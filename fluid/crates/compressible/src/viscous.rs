//! Viscous stresses and heat conduction (compressible Navier–Stokes):
//!   τ = μ(∇u + ∇uᵀ) + (μ_b − 2μ/3)(∇·u) I,   q = −κ∇T,   κ = μ c_p / Pr,
//! added as ∂_d G_d with G_d = −(τ_{·d}, u·τ_{·d} − q_d) moved to the
//! right-hand side. Node gradients and the divergence both use 6th-order
//! central differences (finite-difference form, consistent with the
//! inviscid point-value scheme). Needs 6 ghost cells.

use crate::solver::Mesh;
use crate::system::{System, V};
use rayon::prelude::*;

#[derive(Clone, Copy, Debug)]
pub struct Transport {
    /// dynamic viscosity at the reference temperature
    pub mu: f64,
    /// μ ∝ (T/T_ref)^power (0: constant)
    pub power: f64,
    pub t_ref: f64,
    pub prandtl: f64,
    /// bulk viscosity μ_b
    pub bulk: f64,
}

impl Transport {
    pub fn constant(mu: f64, prandtl: f64) -> Self {
        Transport {
            mu,
            power: 0.0,
            t_ref: 1.0,
            prandtl,
            bulk: 0.0,
        }
    }
    #[inline]
    fn mu_at(&self, t: f64) -> f64 {
        if self.power == 0.0 {
            self.mu
        } else {
            self.mu * (t / self.t_ref).powf(self.power)
        }
    }
}

/// Global frozen-coefficient diffusion estimate. A flux node and the cell
/// receiving its residual can have different densities and heat capacities,
/// so coefficient maxima and capacity minima must be collected separately.
#[derive(Clone, Copy)]
pub(crate) struct DiffusionBounds {
    mu: f64,
    kappa: f64,
    rho: f64,
    capacity: f64,
}

impl Default for DiffusionBounds {
    fn default() -> Self {
        Self {
            mu: 0.0,
            kappa: 0.0,
            rho: f64::INFINITY,
            capacity: f64::INFINITY,
        }
    }
}

impl DiffusionBounds {
    #[inline]
    pub(crate) fn observe<const N: usize, S: System<N>>(
        &mut self,
        sys: &S,
        tr: &Transport,
        w: &V<N>,
    ) {
        self.observe_transport(sys, tr, w);
        let rho = sys.density(w);
        self.rho = self.rho.min(rho);
        self.capacity = self.capacity.min(rho * sys.cv(w));
    }

    /// Ghosts supply flux coefficients but have no evolved heat capacity.
    #[inline]
    pub(crate) fn observe_transport<const N: usize, S: System<N>>(
        &mut self,
        sys: &S,
        tr: &Transport,
        w: &V<N>,
    ) {
        let mu = tr.mu_at(sys.temperature(w));
        self.mu = self.mu.max(mu);
        self.kappa = self.kappa.max(mu * sys.cp(w) / tr.prandtl);
    }

    #[inline]
    pub(crate) fn merge(self, other: Self) -> Self {
        Self {
            mu: self.mu.max(other.mu),
            kappa: self.kappa.max(other.kappa),
            rho: self.rho.min(other.rho),
            capacity: self.capacity.min(other.capacity),
        }
    }

    /// Retain the existing 4/3 thermal safety factor. The momentum rate
    /// includes both the longitudinal shear stress and bulk viscosity.
    pub(crate) fn diffusivity(self, tr: &Transport) -> f64 {
        let momentum = (self.mu * (4.0 / 3.0) + tr.bulk) / self.rho;
        let thermal = self.kappa * (4.0 / 3.0) / self.capacity;
        momentum.max(thermal)
    }
}

const A6: [f64; 3] = [0.75, -0.15, 1.0 / 60.0];

/// Add the viscous terms to `out` (interior cells). `u` must have its
/// ghosts filled (≥ 6 layers on active axes).
pub fn add_viscous<const N: usize, S: System<N>>(
    sys: &S,
    m: &Mesh,
    tr: &Transport,
    u: &[V<N>],
    out: &mut [V<N>],
) {
    assert!(m.ng >= 6, "viscous terms need 6 ghost cells");
    let e_idx = sys
        .energy_index()
        .expect("viscous terms need an energy equation");
    let vi = sys.vel_index();
    let len = m.len();
    let dims: Vec<usize> = (0..3).filter(|&d| m.active(d)).collect();
    let mm = [m.m(0), m.m(1), m.m(2)];
    let stride = [mm[1] * mm[2], mm[2], 1];
    let coord = |q: usize| [q / stride[0], (q / stride[1]) % mm[1], q % mm[2]];
    // primitive fields: u, v, w, T, and μ, c_p
    let prim: Vec<[f64; 6]> = u
        .par_iter()
        .map(|c| {
            let w = sys.to_prim(c);
            let t = sys.temperature(&w);
            [w[vi], w[vi + 1], w[vi + 2], t, tr.mu_at(t), sys.cp(&w)]
        })
        .collect();
    let momentum = std::array::from_fn::<_, 3, _>(|d| sys.momentum_index(d));
    // Cache compact stress/energy fluxes instead of storing gradients. Each
    // node's gradient is consumed here once; neighboring divergence stencils
    // then reuse its flux, with the same arithmetic as the uncached operator.
    let fluxes: Vec<[[f64; 4]; 3]> = (0..len)
        .into_par_iter()
        .map(|q| {
            let c = coord(q);
            let mut f = [[0.0; 4]; 3];
            // Every flux needed by an interior divergence has three valid
            // gradient neighbors in every active direction.
            if dims.iter().any(|&d| c[d] < 3 || c[d] + 3 >= mm[d]) {
                return f;
            }
            let mut g = [[0.0; 3]; 4];
            for &d in &dims {
                for v in 0..4 {
                    let mut s = 0.0;
                    for (l, &a) in A6.iter().enumerate() {
                        let o = (l + 1) * stride[d];
                        s += a * (prim[q + o][v] - prim[q - o][v]);
                    }
                    g[v][d] = s / m.dx[d];
                }
            }
            let p = &prim[q];
            let mu = p[4];
            let div = dims.iter().map(|&d| g[d][d]).sum::<f64>();
            let lam = tr.bulk - 2.0 / 3.0 * mu;
            let kappa = mu * p[5] / tr.prandtl;
            for &d in &dims {
                let mut work = 0.0;
                for i in 0..3 {
                    let mut tau = mu * (g[i][d] + g[d][i]);
                    if i == d {
                        tau += lam * div;
                    }
                    f[d][i] = -tau;
                    work += p[i] * tau;
                }
                f[d][3] = -(work + kappa * g[3][d]);
            }
            f
        })
        .collect();
    // Disjoint output rows avoid an interior-index allocation, a residual
    // staging allocation, and a serial scatter of every cell's result.
    out.par_chunks_mut(stride[0])
        .enumerate()
        .for_each(|(i, row)| {
            if i < m.g(0) || i >= m.g(0) + m.n[0] {
                return;
            }
            for j in m.g(1)..m.g(1) + m.n[1] {
                for k in m.g(2)..m.g(2) + m.n[2] {
                    let r = j * stride[1] + k;
                    let q = i * stride[0] + r;
                    let mut acc = [0.0; N];
                    for &d in &dims {
                        for (l, &a) in A6.iter().enumerate() {
                            let o = (l + 1) * stride[d];
                            let (fp, fm) = (&fluxes[q + o][d], &fluxes[q - o][d]);
                            // Expand only at consumption, preserving the original
                            // per-component accumulation (including zero fields).
                            let mut delta = [0.0; N];
                            for v in 0..3 {
                                delta[momentum[v]] = fp[v] - fm[v];
                            }
                            delta[e_idx] = fp[3] - fm[3];
                            for v in 0..N {
                                acc[v] -= a * delta[v] / m.dx[d];
                            }
                        }
                    }
                    for v in 0..N {
                        row[r][v] += acc[v];
                    }
                }
            }
        });
}
