//! Transported scalars (species, dissolved gas, nuclei number densities)
//! carried by the mixture or by one phase.
//!
//! A scalar is stored as its conserved value ρ_c φ, with ρ_c the density of
//! its carrier (the mixture, or one phase's partial density α_kρ_k: dissolved
//! nitrogen lives in the liquid and moves with the liquid's mass flux). Each
//! step it is updated with the carrier's *effective* face mass fluxes of that
//! step (the Runge–Kutta-weighted sum the solver applied to the carrier):
//!   (ρ_c φ)ⁿ⁺¹ = (ρ_c φ)ⁿ − Δt ∇·(F_c φ_upwind),
//! the update the carrier's own mass satisfies with φ = 1, so a uniform φ
//! stays exactly uniform. The face value is a single-step TVD one (second
//! order in space and time: upwind plus the Lax–Wendroff correction
//! ½(1 − ν)ψ(r)Δφ with the monotonized-central limiter ψ and the face
//! Courant number ν), so φ stays within its initial bounds. Then optional diffusion of φ within the
//! carrier (∂(ρ_cφ)/∂t = ∇·(ρ_c D ∇φ), no flux into cells without carrier)
//! and a per-cell kinetics hook (phase change, dissolution) that may change
//! the cell state and its scalars together.

use crate::solver::{Bc, Mesh, Solver};
use crate::system::{System, V};
use rayon::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Carrier {
    /// the mixture density
    Mixture,
    /// one partial density (conserved index)
    Phase(usize),
}

#[derive(Clone, Debug)]
pub struct Scalar {
    pub name: String,
    pub carrier: Carrier,
    /// conserved value ρ_c φ on the ghost-padded mesh (interior is used)
    pub u: Vec<f64>,
    /// diffusivity of φ within the carrier (m²/s)
    pub diffusivity: f64,
}

impl Scalar {
    /// Scalar with specific value φ(x) (per unit carrier mass) on `solver`'s
    /// mesh and current state.
    pub fn new<const N: usize, S: System<N>>(solver: &Solver<N, S>, name: &str, carrier: Carrier, diffusivity: f64, phi: impl Fn([f64; 3]) -> f64) -> Self {
        let m = &solver.mesh;
        let mut u = vec![0.0; m.len()];
        for i in 0..m.n[0] {
            for j in 0..m.n[1] {
                for k in 0..m.n[2] {
                    let q = m.idx(i + m.g(0), j + m.g(1), k + m.g(2));
                    u[q] = carrier_density(&solver.sys, carrier, &solver.u[q]) * phi(m.x(i, j, k));
                }
            }
        }
        Scalar { name: name.into(), carrier, u, diffusivity }
    }

    /// Specific value φ of interior cell (i, j, k).
    pub fn phi<const N: usize, S: System<N>>(&self, solver: &Solver<N, S>, i: usize, j: usize, k: usize) -> f64 {
        let m = &solver.mesh;
        let q = m.idx(i + m.g(0), j + m.g(1), k + m.g(2));
        let rc = carrier_density(&solver.sys, self.carrier, &solver.u[q]);
        if rc > 0.0 {
            self.u[q] / rc
        } else {
            0.0
        }
    }
}

/// Carrier density of a conserved state.
pub fn carrier_density<const N: usize, S: System<N>>(sys: &S, c: Carrier, u: &V<N>) -> f64 {
    match c {
        Carrier::Mixture => sys.density_indices().map(|k| u[k]).sum(),
        Carrier::Phase(k) => u[k],
    }
}

/// Conserved indices summed into a carrier's mass flux.
pub(crate) fn carrier_indices<const N: usize, S: System<N>>(sys: &S, c: Carrier) -> std::ops::Range<usize> {
    match c {
        Carrier::Mixture => sys.density_indices(),
        Carrier::Phase(k) => k..k + 1,
    }
}

/// Face index layout of the flux accumulators: per axis d, pencil-major
/// (pencil (a, b) over the two other axes in increasing order, faces
/// f = 0..n_d between interior cells f−1 and f).
pub(crate) fn face_count(m: &Mesh, d: usize) -> usize {
    (m.n[d] + 1) * m.n[(d + 1) % 3] * m.n[(d + 2) % 3]
}

impl<const N: usize, S: System<N>> Solver<N, S> {
    /// Distinct carriers of the scalars (order of first use).
    pub(crate) fn carriers(&self) -> Vec<Carrier> {
        let mut c = Vec::new();
        for s in &self.scalars {
            if !c.contains(&s.carrier) {
                c.push(s.carrier);
            }
        }
        c
    }

    /// Zero the per-step carrier flux accumulators.
    pub(crate) fn reset_flux_acc(&mut self) {
        let m = self.mesh.clone();
        let nc = self.carriers().len();
        self.flux_acc = (0..nc).map(|_| std::array::from_fn(|d| if m.active(d) { vec![0.0; face_count(&m, d)] } else { Vec::new() })).collect();
    }

    /// Specific values φ of every scalar (ghost-padded; ghosts by the
    /// boundary conditions: periodic copy, else zero gradient) from state u.
    pub(crate) fn scalar_phis(&mut self, u: &[V<N>]) -> Vec<Vec<f64>> {
        let m = self.mesh.clone();
        let carriers = self.carriers();
        self.carrier_old = carriers
            .iter()
            .map(|c| {
                let mut r: Vec<f64> = u.iter().map(|x| carrier_density(&self.sys, *c, x)).collect();
                self.fill_phi_ghosts(&mut r);
                r
            })
            .collect();
        let m = &m;
        self.scalars
            .iter()
            .map(|s| {
                let mut phi: Vec<f64> = (0..m.len())
                    .map(|q| {
                        let rc = carrier_density(&self.sys, s.carrier, &u[q]);
                        if rc > 1e-300 {
                            s.u[q] / rc
                        } else {
                            0.0
                        }
                    })
                    .collect();
                self.fill_phi_ghosts(&mut phi);
                phi
            })
            .collect()
    }

    fn fill_phi_ghosts(&self, f: &mut [f64]) {
        let m = &self.mesh;
        for d in 0..3 {
            if !m.active(d) {
                continue;
            }
            let (n, ng) = (m.n[d], m.ng);
            let others: Vec<usize> = (0..3).filter(|&e| e != d).collect();
            let (e0, e1) = (others[0], others[1]);
            for a in 0..m.m(e0) {
                for b in 0..m.m(e1) {
                    let at = |c: usize| {
                        let mut ijk = [0; 3];
                        ijk[d] = c;
                        ijk[e0] = a;
                        ijk[e1] = b;
                        m.idx(ijk[0], ijk[1], ijk[2])
                    };
                    for gi in 0..ng {
                        let (lo, hi) = m.periodic_ghost_sources(d, gi);
                        f[at(gi)] = if self.bc[d][0] == Bc::Periodic { f[at(lo)] } else { f[at(ng)] };
                        f[at(n + ng + gi)] = if self.bc[d][1] == Bc::Periodic { f[at(hi)] } else { f[at(n + ng - 1)] };
                    }
                }
            }
        }
    }

    /// Transport, diffusion and kinetics of the scalars over a step of Δt
    /// whose effective carrier face fluxes are in `flux_acc`; `phi_old` are
    /// the specific values at the start of the step. `self.u` is the new
    /// state.
    pub(crate) fn finish_scalars(&mut self, phi_old: &[Vec<f64>], dt: f64) {
        self.finish_scalars_with(phi_old, dt, true);
    }

    /// As `finish_scalars`, optionally without the kinetics hook.
    pub(crate) fn finish_scalars_with(&mut self, phi_old: &[Vec<f64>], dt: f64, kinetics: bool) {
        if self.scalars.is_empty() && (self.kinetics.is_none() || !kinetics) {
            return;
        }
        let m = self.mesh.clone();
        let carriers = self.carriers();
        let dims: Vec<usize> = (0..3).filter(|&d| m.active(d)).collect();
        // TVD transport with the carrier's effective face mass fluxes; the
        // face Courant number needs the old carrier density (from φ_old's
        // companion: ρ_c = (ρ_cφ)/φ is ill-defined, so it is passed along)
        let rc_old = &self.carrier_old;
        for (si, s) in self.scalars.iter_mut().enumerate() {
            let ci = carriers.iter().position(|c| *c == s.carrier).unwrap();
            let phi = &phi_old[si];
            let rco = &rc_old[ci];
            for &d in &dims {
                let others: Vec<usize> = (0..3).filter(|&e| e != d).collect();
                let (e0, e1) = (others[0], others[1]);
                let n = m.n[d];
                let inv = dt / m.dx[d];
                let acc = &self.flux_acc[ci][d];
                for a in 0..m.n[e0] {
                    for b in 0..m.n[e1] {
                        let pidx = a * m.n[e1] + b;
                        let at = |cc: usize| {
                            let mut ijk = [0; 3];
                            ijk[d] = cc;
                            ijk[e0] = a + m.g(e0);
                            ijk[e1] = b + m.g(e1);
                            m.idx(ijk[0], ijk[1], ijk[2])
                        };
                        let ng = m.ng;
                        let flux = |f: usize| -> f64 {
                            let fm = acc[pidx * (n + 1) + f];
                            // upwind cell u, downwind dn, second upwind uu
                            let (u, dn, uu) = if fm > 0.0 { (ng + f - 1, ng + f, ng + f - 2) } else { (ng + f, ng + f - 1, ng + f + 1) };
                            let (pu, pd, puu) = (phi[at(u)], phi[at(dn)], phi[at(uu)]);
                            let dphi = pd - pu;
                            let face = if dphi.abs() > 1e-300 {
                                let r = (pu - puu) / dphi;
                                let psi = (2.0 * r).min(0.5 * (1.0 + r)).min(2.0).max(0.0);
                                let nu = (fm.abs() * inv / rco[at(u)].max(1e-300)).min(1.0);
                                pu + 0.5 * (1.0 - nu) * psi * dphi
                            } else {
                                pu
                            };
                            fm * face
                        };
                        let mut lo = flux(0);
                        for cc in 0..n {
                            let hi = flux(cc + 1);
                            s.u[at(ng + cc)] -= inv * (hi - lo);
                            lo = hi;
                        }
                    }
                }
            }
        }
        // diffusion of φ within the carrier (explicit sub-steps)
        for si in 0..self.scalars.len() {
            let dcoef = self.scalars[si].diffusivity;
            if dcoef <= 0.0 {
                continue;
            }
            let carrier = self.scalars[si].carrier;
            let rc: Vec<f64> = self.u.iter().map(|x| carrier_density(&self.sys, carrier, x)).collect();
            let lim: f64 = dims.iter().map(|&d| 2.0 * dcoef / (m.dx[d] * m.dx[d])).sum();
            let nsub = ((dt * lim / 0.8).ceil() as usize).max(1);
            let h = dt / nsub as f64;
            for _ in 0..nsub {
                let mut phi: Vec<f64> = (0..m.len()).map(|q| if rc[q] > 1e-300 { self.scalars[si].u[q] / rc[q] } else { 0.0 }).collect();
                self.fill_phi_ghosts(&mut phi);
                let mut rc_g = rc.clone();
                self.fill_phi_ghosts(&mut rc_g);
                let s = &mut self.scalars[si];
                let upd: Vec<(usize, f64)> = (0..m.n[0])
                    .into_par_iter()
                    .flat_map_iter(|i| {
                        let mut out = Vec::new();
                        for j in 0..m.n[1] {
                            for k in 0..m.n[2] {
                                let c = [i + m.g(0), j + m.g(1), k + m.g(2)];
                                let q = m.idx(c[0], c[1], c[2]);
                                let mut r = 0.0;
                                for &d in &dims {
                                    for sg in [-1i64, 1] {
                                        let mut cn = c;
                                        cn[d] = (c[d] as i64 + sg) as usize;
                                        let qn = m.idx(cn[0], cn[1], cn[2]);
                                        let interior = cn[d] >= m.g(d) && cn[d] < m.g(d) + m.n[d];
                                        if !interior && self.bc[d][(sg > 0) as usize] != Bc::Periodic {
                                            continue; // no flux through walls / open boundaries
                                        }
                                        // no diffusion into cells without carrier
                                        let rf = rc_g[q].min(rc_g[qn]);
                                        r += rf * dcoef * (phi[qn] - phi[q]) / (m.dx[d] * m.dx[d]);
                                    }
                                }
                                out.push((q, r));
                            }
                        }
                        out
                    })
                    .collect();
                for (q, r) in upd {
                    s.u[q] += h * r;
                }
            }
        }
        // kinetics (phase change, dissolution): cell state and scalars
        if !kinetics {
            return;
        }
        if let Some(kin) = &self.kinetics {
            let ns = self.scalars.len();
            let interior: Vec<usize> = (0..m.n[0])
                .flat_map(|i| (0..m.n[1]).flat_map(move |j| (0..m.n[2]).map(move |k| (i, j, k))))
                .map(|(i, j, k)| m.idx(i + m.g(0), j + m.g(1), k + m.g(2)))
                .collect();
            let res: Vec<(usize, V<N>, Vec<f64>)> = interior
                .par_iter()
                .filter_map(|&q| {
                    let mut uc = self.u[q];
                    let mut sc: Vec<f64> = (0..ns).map(|si| self.scalars[si].u[q]).collect();
                    kin(&self.sys, &mut uc, &mut sc, dt).then_some((q, uc, sc))
                })
                .collect();
            self.kinetics_cells = res.len();
            for (q, uc, sc) in res {
                self.u[q] = uc;
                for si in 0..ns {
                    self.scalars[si].u[q] = sc[si];
                }
            }
        }
    }

    /// Totals of a scalar's conserved value over the interior (× volume).
    pub fn scalar_total(&self, si: usize) -> f64 {
        let m = &self.mesh;
        let vol: f64 = (0..3).filter(|&d| m.active(d)).map(|d| m.dx[d]).product();
        let mut t = 0.0;
        for i in 0..m.n[0] {
            for j in 0..m.n[1] {
                for k in 0..m.n[2] {
                    t += self.scalars[si].u[m.idx(i + m.g(0), j + m.g(1), k + m.g(2))] * vol;
                }
            }
        }
        t
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{eos::Nasg, Euler};

    #[test]
    fn scalar_ghosts_wrap_short_periodic_axes_and_corners() {
        let n = [2, 3, 4];
        for ng in [6, 10] {
            let s = Solver::new(
                Euler {
                    eos: Nasg::ideal(1.4),
                },
                Mesh::new(n, [0.0; 3], [1.0; 3]).with_ghosts(ng),
                |_| [1.0, 0.0, 0.0, 0.0, 1.0],
            );
            let m = &s.mesh;
            let mut phi = vec![f64::NAN; m.len()];
            let value = |i, j, k| (100 * i + 10 * j + k) as f64;
            for i in 0..n[0] {
                for j in 0..n[1] {
                    for k in 0..n[2] {
                        phi[m.idx(i + m.g(0), j + m.g(1), k + m.g(2))] = value(i, j, k);
                    }
                }
            }
            s.fill_phi_ghosts(&mut phi);
            for i in 0..m.m(0) {
                for j in 0..m.m(1) {
                    for k in 0..m.m(2) {
                        let c = [i, j, k];
                        let wrapped: [usize; 3] = std::array::from_fn(|d| {
                            (c[d] as isize - m.g(d) as isize).rem_euclid(n[d] as isize) as usize
                        });
                        assert_eq!(
                            phi[m.idx(i, j, k)],
                            value(wrapped[0], wrapped[1], wrapped[2])
                        );
                    }
                }
            }
        }
    }
}
