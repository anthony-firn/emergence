//! Surface tension for the diffuse-interface multiphase models: continuum
//! surface force (Brackbill, Kothe & Zemach 1992) in the form used for the
//! five-equation model (Perigaud & Saurel 2005):
//!   ∂(ρu)/∂t + … = σ κ ∇α,   ∂E/∂t + … = σ κ ∇α · u,
//!   κ = −∇·(∇α̃/|∇α̃|),
//! with α the volume fraction of the reference phase and α̃ a lightly
//! smoothed copy (1-2-1 filter) for the curvature only — the raw diffuse
//! interface of a shock-capturing scheme gives noisy normals. Derivatives
//! are 4th-order central differences on the ghost-filled mesh.

use crate::solver::Mesh;
use crate::system::{System, V};
use rayon::prelude::*;

#[derive(Clone, Copy, Debug)]
pub struct SurfaceTension {
    /// surface tension coefficient σ (N/m)
    pub sigma: f64,
    /// conserved index of the volume fraction α
    pub alpha_index: usize,
    /// smoothing passes for the curvature field
    pub smoothing: usize,
}

const D4: [f64; 2] = [2.0 / 3.0, -1.0 / 12.0];

impl SurfaceTension {
    /// Add σκ∇α (momentum) and σκ∇α·u (energy) to `out` on interior cells.
    /// `u` must have ghosts filled (≥ 4 layers on active axes).
    pub fn add<const N: usize, S: System<N>>(&self, sys: &S, m: &Mesh, u: &[V<N>], out: &mut [V<N>]) {
        let mm = [m.m(0), m.m(1), m.m(2)];
        let stride = [mm[1] * mm[2], mm[2], 1];
        let len = m.len();
        let dims: Vec<usize> = (0..3).filter(|&d| m.active(d)).collect();
        let coord = |q: usize| [q / stride[0], (q / stride[1]) % mm[1], q % mm[2]];
        let fits = |q: usize, d: usize, w: usize| {
            let c = coord(q);
            c[d] >= w && c[d] + w < mm[d]
        };
        let alpha: Vec<f64> = u.par_iter().map(|c| c[self.alpha_index]).collect();
        // smoothed copy for the normals
        let mut sm = alpha.clone();
        for _ in 0..self.smoothing {
            for &d in &dims {
                let prev = sm.clone();
                sm.par_iter_mut().enumerate().for_each(|(q, v)| {
                    if fits(q, d, 1) {
                        *v = 0.25 * prev[q - stride[d]] + 0.5 * prev[q] + 0.25 * prev[q + stride[d]];
                    }
                });
            }
        }
        let grad = |f: &[f64], q: usize, d: usize| -> f64 {
            if !fits(q, d, 2) {
                return 0.0;
            }
            let s = stride[d];
            (D4[0] * (f[q + s] - f[q - s]) + D4[1] * (f[q + 2 * s] - f[q - 2 * s])) / m.dx[d]
        };
        // unit normal field of the smoothed fraction
        let normal: Vec<[f64; 3]> = (0..len)
            .into_par_iter()
            .map(|q| {
                let mut g = [0.0; 3];
                for &d in &dims {
                    g[d] = grad(&sm, q, d);
                }
                let nrm = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
                if nrm > 1e-8 / m.dx.iter().cloned().fold(f64::INFINITY, f64::min) {
                    [g[0] / nrm, g[1] / nrm, g[2] / nrm]
                } else {
                    [0.0; 3]
                }
            })
            .collect();
        let vi = sys.vel_index();
        let e_idx = sys.energy_index();
        let interior: Vec<usize> = (0..m.n[0])
            .flat_map(|i| (0..m.n[1]).flat_map(move |j| (0..m.n[2]).map(move |k| (i, j, k))))
            .map(|(i, j, k)| m.idx(i + m.g(0), j + m.g(1), k + m.g(2)))
            .collect();
        let adds: Vec<(usize, [f64; 3], f64)> = interior
            .par_iter()
            .map(|&q| {
                // κ = −∇·n
                let mut kappa = 0.0;
                for &d in &dims {
                    if fits(q, d, 2) {
                        let s = stride[d];
                        kappa -= (D4[0] * (normal[q + s][d] - normal[q - s][d]) + D4[1] * (normal[q + 2 * s][d] - normal[q - 2 * s][d])) / m.dx[d];
                    }
                }
                let mut f = [0.0; 3];
                for &d in &dims {
                    f[d] = self.sigma * kappa * grad(&alpha, q, d);
                }
                let w = sys.to_prim(&u[q]);
                let work = f[0] * w[vi] + f[1] * w[vi + 1] + f[2] * w[vi + 2];
                (q, f, work)
            })
            .collect();
        for (q, f, work) in adds {
            for d in 0..3 {
                out[q][sys.momentum_index(d)] += f[d];
            }
            if let Some(e) = e_idx {
                out[q][e] += work;
            }
        }
    }
}
