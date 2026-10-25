//! Heat conduction for single- and multiphase flows, with phase change.
//!
//!   ∂E/∂t = ∇·(k∇T),  k = Σ α_k k_k (phase conductivities), T the
//!   mixture temperature (phase temperatures weighted by mass; equal in
//!   cells relaxed to phase equilibrium),
//! integrated in explicit sub-steps after the flow step (in the all-speed
//! mode: before the pressure solve, as a volume source). After every
//! sub-step the solver's phase relaxation (`Solver::relax`) runs, so heat
//! reaching a cell at its melting or boiling point turns into latent heat at
//! that temperature (the enthalpy method, here with the full mixture EOS):
//! a melting or freezing front moves exactly as fast as conduction brings
//! or removes the latent heat (Stefan problem, `examples/stefan`).
//! Face conductivities are harmonic means; walls hold a temperature
//! (Dirichlet) or are adiabatic.

use crate::solver::{Bc, Solver};
use crate::system::{System, V};
use rayon::prelude::*;

#[derive(Clone, Debug)]
pub struct Conduction {
    /// conductivity of each phase (order of `System::phases`), W/(m K)
    pub k: Vec<f64>,
    /// wall temperature per axis and side (None: adiabatic; periodic axes
    /// wrap)
    pub wall: [[Option<f64>; 2]; 3],
    /// Sub-step safety factor (0 < safety ≤ 0.5):
    /// h ≤ 2*safety*min(ρc_v / sum of face conductances).
    pub safety: f64,
}

impl Conduction {
    pub fn new(k: Vec<f64>) -> Self {
        Conduction {
            k,
            wall: [[None; 2]; 3],
            safety: 0.4,
        }
    }
}

/// Mixture temperature: phase temperatures weighted by phase mass, phases
/// with a negligible volume fraction ignored.
pub fn mixture_temperature(phases: &[(f64, f64, f64)], partial: &[f64]) -> f64 {
    let (mut s, mut m) = (0.0, 0.0);
    for (k, &(al, t, _)) in phases.iter().enumerate() {
        if al > 1e-8 && t.is_finite() {
            s += partial[k] * t;
            m += partial[k];
        }
    }
    if m > 0.0 {
        s / m
    } else {
        f64::NAN
    }
}

/// Common phase sums for model-specific, allocation-free thermal closures.
pub(crate) fn properties(
    t: f64,
    phases: &[(f64, f64, f64)],
    conductivity: &[f64],
) -> (f64, f64, f64) {
    assert_eq!(
        phases.len(),
        conductivity.len(),
        "one conductivity per phase is required"
    );
    let k: f64 = phases
        .iter()
        .zip(conductivity)
        .map(|(p, k)| p.0.max(0.0) * k)
        .sum();
    let capacity: f64 = phases
        .iter()
        .map(|p| if p.0 > 1e-8 { p.0 * p.2 } else { 0.0 })
        .sum();
    (t, k, capacity)
}

impl<const N: usize, S: System<N>> Solver<N, S> {
    /// Conduct heat over Δt in `u` (ghost-padded layout; interior updated),
    /// relaxing phases after every sub-step. Returns the sub-step count.
    pub(crate) fn conduct(&mut self, u: &mut [V<N>], dt: f64) -> usize {
        let Some(cd) = self.conduction.clone() else {
            return 0;
        };
        assert!(dt.is_finite() && dt >= 0.0);
        assert!(cd.safety.is_finite() && cd.safety > 0.0 && cd.safety <= 0.5);
        assert!(cd.k.iter().all(|k| k.is_finite() && *k >= 0.0));
        assert!(cd
            .wall
            .iter()
            .flatten()
            .flatten()
            .all(|t| t.is_finite() && *t >= 0.0));
        let m = self.mesh.clone();
        let dims: Vec<usize> = (0..3).filter(|&d| m.active(d)).collect();
        let e_idx = self
            .sys
            .energy_index()
            .expect("heat conduction needs an energy equation");
        // Compact interior properties and rates, reused across sub-steps.
        let mut cell = vec![(0.0, 0.0, 0.0); m.cells()];
        let mut rates = vec![(0.0, 0.0); m.cells()];
        let mut t_done = 0.0;
        let mut n = 0;
        while t_done < dt * (1.0 - 1e-12) {
            // Temperature, conductivity and volumetric heat capacity ρc_v.
            let sys = &self.sys;
            cell.par_chunks_mut(m.n[2]).enumerate().for_each(|(ij, row)| {
                let i = ij / m.n[1] + m.g(0);
                let j = ij % m.n[1] + m.g(1);
                let first = m.idx(i, j, m.g(2));
                for (k, value) in row.iter_mut().enumerate() {
                    let w = sys.to_prim(&u[first + k]);
                    let (t, conductivity, cv) = sys.conduction_properties(&w, &cd.k);
                    assert!(
                        t.is_finite() && t >= 0.0 && conductivity.is_finite() && cv.is_finite() && cv > 0.0,
                        "conduction requires finite temperature, conductivity and positive heat capacity"
                    );
                    *value = (t, conductivity, cv);
                }
            });
            let lin = |c: [usize; 3]| (c[0] * m.n[1] + c[1]) * m.n[2] + c[2];
            let bc = self.bc;
            let lim = rates
                .par_iter_mut()
                .enumerate()
                .map(|(l, value)| {
                    let c = [l / (m.n[1] * m.n[2]), (l / m.n[2]) % m.n[1], l % m.n[2]];
                    let (t0, k0, _) = cell[l];
                    let (mut r, mut conductance) = (0.0, 0.0);
                    for &d in &dims {
                        for sg in [-1i64, 1] {
                            let v = c[d] as i64 + sg;
                            let nn = m.n[d] as i64;
                            let side = (sg > 0) as usize;
                            let (tn, kn, half) = if v < 0 || v >= nn {
                                if bc[d][side] == Bc::Periodic {
                                    let mut cc = c;
                                    cc[d] = v.rem_euclid(nn) as usize;
                                    let x = cell[lin(cc)];
                                    (x.0, x.1, false)
                                } else if let Some(tw) = cd.wall[d][side] {
                                    (tw, k0, true)
                                } else {
                                    continue;
                                }
                            } else {
                                let mut cc = c;
                                cc[d] = v as usize;
                                let x = cell[lin(cc)];
                                (x.0, x.1, false)
                            };
                            let kf = if k0 + kn > 0.0 {
                                2.0 * k0 * kn / (k0 + kn)
                            } else {
                                0.0
                            };
                            // a wall face is half a cell away
                            let dist2 = if half {
                                0.5 * m.dx[d] * m.dx[d]
                            } else {
                                m.dx[d] * m.dx[d]
                            };
                            r += kf * (tn - t0) / dist2;
                            conductance += kf / dist2;
                        }
                    }
                    assert!(
                        r.is_finite() && conductance.is_finite() && conductance >= 0.0,
                        "non-finite conductive face rate"
                    );
                    *value = (r, conductance);
                    // h*sum(G)/(rho*cv) <= 2*safety <= 1, so frozen positive-
                    // capacity temperatures update as a convex combination.
                    if conductance > 0.0 {
                        cell[l].2 / conductance
                    } else {
                        f64::INFINITY
                    }
                })
                .reduce(|| f64::INFINITY, f64::min);
            let h = (cd.safety * 2.0 * lim).min(dt - t_done);
            assert!(h.is_finite() && h > 0.0, "conduction timestep collapsed");
            let update = |state: &mut V<N>, l: usize| {
                state[e_idx] += h * rates[l].0;
                // Preserve the relaxation hook's commit-on-true semantics.
                if let Some(rx) = &self.relax {
                    let mut relaxed = *state;
                    if rx(sys, &mut relaxed, h) {
                        *state = relaxed;
                    }
                }
            };
            let update_row = |ij: usize, row: &mut [V<N>]| {
                let i = ij / m.m(1);
                let j = ij % m.m(1);
                if i < m.g(0) || i >= m.g(0) + m.n[0] || j < m.g(1) || j >= m.g(1) + m.n[1] {
                    return;
                }
                let l0 = ((i - m.g(0)) * m.n[1] + j - m.g(1)) * m.n[2];
                for k in 0..m.n[2] {
                    update(&mut row[m.g(2) + k], l0 + k);
                }
            };
            if m.cells() >= 8192 {
                u.par_chunks_mut(m.m(2))
                    .enumerate()
                    .for_each(|(ij, row)| update_row(ij, row));
            } else {
                for (ij, row) in u.chunks_mut(m.m(2)).enumerate() {
                    update_row(ij, row);
                }
            }
            t_done += h;
            n += 1;
        }
        self.conduction_substeps += n;
        n
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{eos::Nasg, multiphase::FiveEq, Euler, Mesh};

    #[test]
    fn half_cell_thermal_wall_does_not_drive_a_hot_cell_negative() {
        let mut s = Solver::new(
            Euler {
                eos: Nasg::ideal(1.4),
            },
            Mesh::new([3, 1, 1], [0.0; 3], [3.0, 1.0, 1.0]),
            |x| {
                [
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.4 * if x[0] > 2.0 { 100.0 } else { 1.0 },
                ]
            },
        );
        s.bc[0] = [Bc::Outflow; 2];
        let mut cd = Conduction::new(vec![1.0]);
        cd.wall[0][1] = Some(1.0);
        s.conduction = Some(cd);
        let mut u = std::mem::take(&mut s.u);
        let count = s.conduct(&mut u, 0.4);
        s.u = u;
        assert!(count >= 2); // old cell-only bound took one step to T=-18.8
        for i in 0..3 {
            let w = s.prim(i, 0, 0);
            let t = s.sys.temperature(&w);
            assert!((1.0 - 1e-12..=100.0 + 1e-12).contains(&t));
            assert!(s.sys.admissible(&w));
        }
    }

    #[test]
    fn conductivity_contrast_stays_positive_and_conserves_closed_energy() {
        let sys = FiveEq {
            eos: [Nasg::ideal(1.4); 2],
            kapila: true,
        };
        let mut s = Solver::new(sys, Mesh::new([3, 1, 1], [0.0; 3], [3.0, 1.0, 1.0]), |x| {
            if (1.0..2.0).contains(&x[0]) {
                [1e-6, 0.0, 0.0, 0.0, 0.0, 0.4 * 1e-6 * 100.0, 1.0]
            } else {
                [0.0, 1.0, 0.0, 0.0, 0.0, 0.4, 0.0]
            }
        });
        s.bc[0] = [Bc::Outflow; 2];
        s.conduction = Some(Conduction::new(vec![1.0, 1000.0]));
        let before = s.totals()[5];
        let mut u = std::mem::take(&mut s.u);
        s.conduct(&mut u, 4e-7);
        s.u = u;
        assert!((s.totals()[5] - before).abs() < 1e-13);
        for i in 0..3 {
            let w = s.prim(i, 0, 0);
            let t = s.sys.temperature(&w);
            assert!((1.0 - 1e-12..=100.0 + 1e-12).contains(&t));
            assert!(s.sys.admissible(&w));
        }
    }

    #[test]
    fn periodic_conduction_retains_second_order_cosine_convergence() {
        let mut errors = Vec::new();
        for n in [16, 32, 64] {
            let exact = |x: f64, t: f64| {
                2.0 + 0.2
                    * (-std::f64::consts::TAU.powi(2) * t).exp()
                    * (std::f64::consts::TAU * x).cos()
            };
            let mut s = Solver::new(
                Euler {
                    eos: Nasg::ideal(1.4),
                },
                Mesh::new([n, 1, 1], [0.0; 3], [1.0; 3]),
                |x| [1.0, 0.0, 0.0, 0.0, 0.4 * exact(x[0], 0.0)],
            );
            s.conduction = Some(Conduction::new(vec![1.0]));
            let before = s.totals()[4];
            let mut u = std::mem::take(&mut s.u);
            s.conduct(&mut u, 0.02);
            s.u = u;
            assert!((s.totals()[4] - before).abs() < 1e-13);
            let error = ((0..n)
                .map(|i| {
                    (s.sys.temperature(&s.prim(i, 0, 0)) - exact(s.mesh.x(i, 0, 0)[0], 0.02))
                        .powi(2)
                })
                .sum::<f64>()
                / n as f64)
                .sqrt();
            println!("conduction n={n} rms={error:.9e}");
            errors.push(error);
        }
        for pair in errors.windows(2) {
            assert!((3.7..4.3).contains(&(pair[0] / pair[1])));
        }
    }

    fn cosine_error<const N: usize, S: System<N>>(
        sys: S,
        init: impl Fn(f64, f64) -> V<N> + Sync,
        n: usize,
        conductivities: Vec<f64>,
    ) -> f64 {
        let exact = |x: f64, t: f64| {
            2.0 + 0.2
                * (-std::f64::consts::TAU.powi(2) * t).exp()
                * (std::f64::consts::TAU * x).cos()
        };
        let mut s = Solver::new(sys, Mesh::new([n, 1, 1], [0.0; 3], [1.0; 3]), |x| {
            init(x[0], exact(x[0], 0.0))
        });
        s.conduction = Some(Conduction::new(conductivities));
        let before = s.totals();
        let mut u = std::mem::take(&mut s.u);
        s.conduct(&mut u, 0.02);
        s.u = u;
        let after = s.totals();
        for (a, b) in before.iter().zip(&after) {
            assert!((a - b).abs() < 1e-12);
        }
        ((0..n)
            .map(|i| {
                (s.sys.temperature(&s.prim(i, 0, 0)) - exact(s.mesh.x(i, 0, 0)[0], 0.02)).powi(2)
            })
            .sum::<f64>()
            / n as f64)
            .sqrt()
    }

    #[test]
    fn multiphase_conduction_retains_cosine_accuracy_with_distinct_heat_capacities() {
        use crate::threephase::Kapila3;
        for family in 0..2 {
            let mut errors = Vec::new();
            for n in [16, 32, 64] {
                let error = if family == 0 {
                    let sys = FiveEq {
                        eos: [
                            Nasg::ideal(1.4),
                            Nasg {
                                cv: 2.0,
                                ..Nasg::ideal(1.4)
                            },
                        ],
                        kapila: true,
                    };
                    // rho_k*cv_k = 1 in both phases, so mixture capacity=1.
                    cosine_error(
                        sys,
                        |_, t| [0.4, 0.3, 0.0, 0.0, 0.0, 0.4 * t, 0.4],
                        n,
                        vec![1.0; 2],
                    )
                } else {
                    let sys = Kapila3 {
                        eos: [
                            Nasg {
                                cv: 0.75,
                                ..Nasg::ideal(1.4)
                            },
                            Nasg::ideal(1.4),
                            Nasg {
                                cv: 1.25,
                                ..Nasg::ideal(1.4)
                            },
                        ],
                    };
                    cosine_error(
                        sys,
                        |_, t| {
                            [
                                0.2 / 0.75,
                                0.3,
                                0.5 / 1.25,
                                0.0,
                                0.0,
                                0.0,
                                0.4 * t,
                                0.2,
                                0.3,
                            ]
                        },
                        n,
                        vec![1.0; 3],
                    )
                };
                errors.push(error);
            }
            for pair in errors.windows(2) {
                assert!(
                    (3.7..4.3).contains(&(pair[0] / pair[1])),
                    "family {family}, errors {errors:?}"
                );
            }
        }
    }

    #[test]
    fn thermal_relaxation_preserves_commit_semantics_and_parallel_reproducibility() {
        use std::sync::{
            atomic::{AtomicUsize, Ordering},
            Arc,
        };
        let one = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        let three = rayon::ThreadPoolBuilder::new()
            .num_threads(3)
            .build()
            .unwrap();
        for shape in [[1, 96, 96], [96, 1, 96], [96, 96, 1], [12, 1, 1]] {
            let evolve = |hook: u8| {
                let mut s = Solver::new(
                    Euler {
                        eos: Nasg::ideal(1.4),
                    },
                    Mesh::new(shape, [0.0; 3], [1.0; 3]),
                    |x| {
                        [
                            1.0,
                            0.0,
                            0.0,
                            0.0,
                            0.4 * (2.0
                                + 0.1
                                    * x.into_iter()
                                        .map(|v| (std::f64::consts::TAU * v).cos())
                                        .sum::<f64>()),
                        ]
                    },
                );
                s.conduction = Some(Conduction::new(vec![1.0]));
                let energy = s.totals()[4];
                let calls = Arc::new(AtomicUsize::new(0));
                if hook > 0 {
                    let calls = calls.clone();
                    s.relax = Some(Box::new(move |_, state, h| {
                        calls.fetch_add(1, Ordering::Relaxed);
                        state[4] *= if hook == 1 { 0.0 } else { (-h).exp() };
                        hook == 2
                    }));
                }
                let mut u = std::mem::take(&mut s.u);
                let count = s.conduct(&mut u, 0.001);
                s.u = u;
                assert_eq!(
                    calls.load(Ordering::Relaxed),
                    if hook > 0 { count * s.mesh.cells() } else { 0 }
                );
                let expected = energy * if hook == 2 { (-0.001_f64).exp() } else { 1.0 };
                assert!((s.totals()[4] / expected - 1.0).abs() < 1e-12);
                for i in 0..shape[0] {
                    for j in 0..shape[1] {
                        for k in 0..shape[2] {
                            assert!(s.sys.admissible(&s.prim(i, j, k)));
                        }
                    }
                }
                s.u.into_iter()
                    .map(|c| c.map(f64::to_bits))
                    .collect::<Vec<_>>()
            };
            assert_eq!(one.install(|| evolve(0)), three.install(|| evolve(1)));
            assert_eq!(one.install(|| evolve(2)), three.install(|| evolve(2)));
        }
    }
}
