//! Cross-crate hand-offs through the emergence engine, each checked against
//! the same chain computed with the exact model everywhere.
//!
//! Mie → radiation: the droplet size distribution of every cloud cell is
//! folded into extinction by the photonics Mie ladder, the result is handed to
//! `radiation::Volume`, and the observable is the transmittance of a set of
//! rays. The engine budgets the error of every (cell, size) evaluation by its
//! adjoint weight `∂T_p/∂Q = −T_p ℓ_pc N_c w_n π r_n²`, so the accuracy
//! requested of each droplet is exactly what the rays need.
//!
//! `cargo test -p testlab --test emergence_links -- --nocapture`.

use std::f64::consts::PI;
use std::time::Instant;

use emergence::Request;
use photonics::mie::GammaSizeDistribution;
use photonics::mie_ladder::{self, MieInput, MieQ};
use photonics::C;
use radiation::math::{Bounds, Ray, Vec3};
use radiation::medium::{Medium, Volume, Voxel};
use rayon::prelude::*;

/// Liquid water near 2.13 µm (estimate from Hale & Querry 1973).
const WATER_2_13UM: C = C::new(1.29, 5.0e-4);
const WAVELENGTH: f64 = 2.13e-6;
const N_CELL: usize = 12;
const SIDE: f64 = 200.0; // m
const N_QUAD: usize = 400;

/// Smooth cloud fields: droplet number density (m⁻³; ~10 cm⁻³, thin
/// stratocumulus scale, estimate) and effective radius 5–15 µm.
fn cell_dist(c: [usize; 3]) -> (f64, GammaSizeDistribution) {
    let s = c.map(|i| (i as f64 + 0.5) / N_CELL as f64);
    let n = 1.0e7 * (1.0 + 0.5 * (2.0 * PI * s[0]).sin() * (PI * s[1]).cos());
    let r_eff = 5e-6 + 10e-6 * s[2] * (0.6 + 0.4 * (3.0 * s[0] + s[1]).cos().abs());
    (n, GammaSizeDistribution { r_eff, v_eff: 0.1 })
}

fn bounds() -> Bounds {
    Bounds {
        lo: Vec3([0.0; 3]),
        hi: Vec3([SIDE; 3]),
    }
}

/// 8 × 8 rays along +x and 16 slanted rays.
fn rays() -> Vec<Ray> {
    let mut r = Vec::new();
    for j in 0..8 {
        for k in 0..8 {
            let y = SIDE * (j as f64 + 0.5) / 8.0;
            let z = SIDE * (k as f64 + 0.5) / 8.0;
            r.push(Ray::new(Vec3([-1.0, y, z]), Vec3([1.0, 0.0, 0.0])));
        }
    }
    for k in 0..16 {
        let z = SIDE * (k as f64 + 0.5) / 16.0;
        r.push(Ray::new(Vec3([-1.0, 1.0, z]), Vec3([1.0, 0.7, 0.2 * (k as f64 / 15.0 - 0.5)])));
    }
    r
}

/// Simpson nodes/weights on [0, r_max] (as `mie_ladder::gamma_averaged`).
fn nodes(d: &GammaSizeDistribution) -> Vec<(f64, f64)> {
    let r_max = d.r_eff * (1.0 + 8.0 * d.v_eff.sqrt());
    let h = r_max / N_QUAD as f64;
    (0..=N_QUAD)
        .map(|i| {
            let w = h / 3.0 * if i == 0 || i == N_QUAD { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 };
            (i as f64 * h, w)
        })
        .collect()
}

/// Build the radiation volume from per-cell size-averaged efficiencies
/// `q[c]` (and number densities) and return the transmittance of every ray.
fn transmittances(q: &[(f64, f64, MieQ)], rays: &[Ray]) -> Vec<f64> {
    let media: Vec<Medium> = q
        .iter()
        .map(|&(n, area, q)| {
            let g = if q.q_sca > 0.0 { q.gq_sca / q.q_sca } else { 0.0 };
            Medium::grey(n * area * (q.q_ext - q.q_sca), n * area * q.q_sca, g)
        })
        .collect();
    let cells = (0..media.len())
        .map(|i| Voxel {
            medium: i,
            scale: 1.0,
            temperature: 0.0,
        })
        .collect();
    let v = Volume::new([N_CELL; 3], bounds(), media, cells);
    rays.iter().map(|&r| v.transmittance(r, f64::INFINITY, WAVELENGTH)).collect()
}

#[test]
fn mie_to_radiation_transmittance_is_certified_end_to_end() {
    let rays = rays();
    let n_cells = N_CELL * N_CELL * N_CELL;
    // Path length of every ray in every cell (the transmittance sensitivity).
    let geometry = Volume::new(
        [N_CELL; 3],
        bounds(),
        vec![Medium::grey(0.0, 0.0, 0.0)],
        vec![
            Voxel {
                medium: 0,
                scale: 1.0,
                temperature: 0.0
            };
            n_cells
        ],
    );
    let mut path = vec![vec![0.0; n_cells]; rays.len()];
    for (p, &r) in rays.iter().enumerate() {
        for s in geometry.segments(r, f64::INFINITY) {
            path[p][s.cell] += s.end - s.start;
        }
    }
    let cells: Vec<(f64, GammaSizeDistribution)> = (0..n_cells)
        .map(|i| cell_dist([i / (N_CELL * N_CELL), (i / N_CELL) % N_CELL, i % N_CELL]))
        .collect();
    let k = 2.0 * PI / WAVELENGTH;
    let x_max = k * cells.iter().map(|c| c.1.r_eff).fold(0.0, f64::max) * (1.0 + 8.0 * 0.1f64.sqrt()) * 1.01;

    let t0 = Instant::now();
    let ladder = mie_ladder::ladder(vec![mie_ladder::closure_rung(WATER_2_13UM, 0.3, x_max, 1e-7, 1 << 15)], 0);
    let t_build = t0.elapsed();

    // Items: every (cell, node) of every cell a ray crosses.
    let touched: Vec<usize> = (0..n_cells).filter(|&c| path.iter().any(|p| p[c] > 0.0)).collect();
    let mut inputs = Vec::new();
    let mut owner = Vec::new(); // (cell, area weight w_n π r² pdf)
    for &c in &touched {
        for (r, w) in nodes(&cells[c].1) {
            inputs.push(MieInput { x: k * r, m: WATER_2_13UM });
            owner.push((c, w * PI * r * r * cells[c].1.pdf(r)));
        }
    }

    // Size-averaged efficiencies per cell from item values.
    let assemble = |values: &[MieQ]| -> Vec<(f64, f64, MieQ)> {
        let mut acc = vec![(0.0, [0.0; 3]); n_cells];
        for ((c, aw), q) in owner.iter().zip(values) {
            acc[*c].0 += aw;
            acc[*c].1[0] += aw * q.q_ext;
            acc[*c].1[1] += aw * q.q_sca;
            acc[*c].1[2] += aw * q.gq_sca;
        }
        (0..n_cells)
            .map(|c| {
                let (a, s) = acc[c];
                let q = if a > 0.0 {
                    MieQ {
                        q_ext: s[0] / a,
                        q_sca: s[1] / a,
                        gq_sca: s[2] / a,
                    }
                } else {
                    MieQ::default()
                };
                (cells[c].0, a, q)
            })
            .collect()
    };

    // Adjoint weights from a pilot transmittance (the linearisation point):
    // W_i = max_p T_p ℓ_pc N_c (w_n π r² pdf), for the bound Σ W_i e_i ≤ tol.
    let tol = 1e-6;
    let t0 = Instant::now();
    let pilot = ladder.solve_functional(&inputs, &vec![1.0 / inputs.len() as f64; inputs.len()], Request::new(1e-3));
    let t_pilot = transmittances(&assemble(&pilot.items.iter().map(|s| s.value).collect::<Vec<_>>()), &rays);
    // Per-cell sensitivity max_p T_p ℓ_pc (margin 1.1 for T_true ≤ T_pilot e^{δτ}).
    let sens: Vec<f64> = (0..n_cells)
        .map(|c| (0..rays.len()).map(|p| t_pilot[p] * path[p][c]).fold(0.0, f64::max) * cells[c].0 * 1.1)
        .collect();
    let weights: Vec<f64> = owner.iter().map(|&(c, aw)| sens[c] * aw).collect();
    let dt_pilot = t0.elapsed();
    let pilot_cost = pilot.cost;
    let f = ladder.refine_functional(pilot, &inputs, &weights, Request::new(tol));
    let dt_final = t0.elapsed() - dt_pilot;
    let values: Vec<MieQ> = f.items.iter().map(|s| s.value).collect();
    let t_engine = transmittances(&assemble(&values), &rays);

    // A posteriori bound per ray from the final per-item claims:
    // δτ_p = Σ_c ℓ_pc N_c Σ_{i∈c} (w π r² pdf)_i e_i, |ΔT_p| ≤ T_p (e^{δτ_p} − 1).
    let mut cell_err = vec![0.0; n_cells];
    for (&(c, aw), s) in owner.iter().zip(&f.items) {
        cell_err[c] += aw * s.claim.error;
    }
    let bound: Vec<f64> = (0..rays.len())
        .map(|p| {
            let dtau: f64 = (0..n_cells).map(|c| path[p][c] * cells[c].0 * cell_err[c]).sum();
            t_engine[p] * dtau.exp_m1()
        })
        .collect();
    let dt_engine = t0.elapsed();

    // Reference: exact Mie at every node, same radiation code.
    let t0 = Instant::now();
    let exact: Vec<MieQ> = inputs.par_iter().map(|i| mie_ladder::exact(i.x, i.m)).collect();
    let t_exact = transmittances(&assemble(&exact), &rays);
    let dt_exact = t0.elapsed();

    // The reference is itself exact only to 2e-8 per node (measured Wiscombe
    // truncation): δτ_ref = Σ_c ℓ_pc N_c <πr²>_c · 2e-8.
    let reference = assemble(&exact);
    let mut worst: f64 = 0.0;
    for p in 0..rays.len() {
        let err = (t_engine[p] - t_exact[p]).abs();
        let dtau_ref: f64 = (0..n_cells).map(|c| path[p][c] * reference[c].0 * reference[c].1 * 2e-8).sum();
        let ref_err = t_exact[p] * dtau_ref.exp_m1();
        assert!(bound[p] <= tol * 1.0001, "ray {p}: bound {:.2e} above tolerance", bound[p]);
        assert!(err <= bound[p] + ref_err + 1e-15, "ray {p}: error {err:.2e} > bound {:.2e}", bound[p]);
        worst = worst.max(err);
    }
    let h = f.histogram(ladder.len());
    let t_range = t_exact.iter().fold((1.0f64, 0.0f64), |(a, b), &t| (a.min(t), b.max(t)));
    let cost_exact: f64 = inputs.iter().map(|i| mie_ladder::exact_cost(i.x)).sum();
    println!(
        "{} rays through {} cells ({} items): T ∈ [{:.3}, {:.3}], worst |ΔT| {worst:.2e} (bound ≤ {tol:e});\n\
         rungs {:?} (rayleigh, closure, exact); predicted cost {:.1}× below exact;\n\
         wall: engine {dt_engine:.2?} (pilot {dt_pilot:.2?}, final {dt_final:.2?}) + closure {t_build:.2?} vs exact {dt_exact:.2?}",
        rays.len(),
        touched.len(),
        inputs.len(),
        t_range.0,
        t_range.1,
        h,
        cost_exact / (pilot_cost + f.cost)
    );
}
