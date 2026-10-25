//! Atomistic-rung benchmarks: mean-field and Monte-Carlo Curie points,
//! spin-wave dispersion against LLG, and the discrete domain-wall fold.
//!
//! Every number is checked against an exact identity or a published result;
//! the measured errors are reported in `analysis/magnetism.md`.

use magnetism::geometry::LatticeKind;
use magnetism::micromag::discrete_wall_chain;
use magnetism::rng::Rng;
use magnetism::spin::{
    measure_frequency, mean_field_tc, spin_wave_omega, CuriePoint, Heisenberg, McConfig,
    SpinChain,
};

fn rel(actual: f64, expected: f64) -> f64 {
    (actual - expected).abs() / expected.abs().max(1e-300)
}

#[test]
fn mean_field_curie_temperature() {
    // k_B T_C = J z S^2 / 3; z = 6, 8, 12 for sc, bcc, fcc.
    assert!(rel(mean_field_tc(LatticeKind::Sc, 1.0, 1.0), 2.0) < 1e-14);
    assert!(rel(mean_field_tc(LatticeKind::Bcc, 1.0, 1.0), 8.0 / 3.0) < 1e-14);
    assert!(rel(mean_field_tc(LatticeKind::Fcc, 1.0, 1.0), 4.0) < 1e-14);
    // Spin magnitude enters as S^2.
    assert!(rel(mean_field_tc(LatticeKind::Sc, 1.0, 0.5), 0.5) < 1e-14);
    // The published classical sc Heisenberg transition is 1.44299(14) J.
    assert!((CuriePoint::MC_SC_HEISENBERG - 1.442_99).abs() < 1e-9);
}

fn u4(kind: LatticeKind, l: usize, t: f64, seed: u64) -> f64 {
    let cfg = McConfig {
        equil_sweeps: 1800,
        samples: 400,
        sweeps_per_sample: 4,
    };
    let mut h = Heisenberg::new(kind, l, 1.0, 1.0, 0.0, [0.0, 0.0, 1.0], [0.0; 3], 1234 + l as u64);
    let mut rng = Rng::new(seed + l as u64 * 101);
    h.binder_cumulant(t, &cfg, &mut rng).u4
}

/// Temperature where two Binder cumulants cross, by linear interpolation.
fn crossing(kind: LatticeKind, l1: usize, l2: usize, seed_a: u64, seed_b: u64) -> Option<f64> {
    let temps: Vec<f64> = (0..13).map(|i| 1.28 + 0.02 * i as f64).collect();
    let mut prev: Option<(f64, f64)> = None;
    for &t in &temps {
        let d = u4(kind, l1, t, seed_a) - u4(kind, l2, t, seed_b);
        if let Some((tp, dp)) = prev {
            if dp * d < 0.0 {
                return Some(tp + (t - tp) * dp.abs() / (dp.abs() + d.abs()));
            }
        }
        prev = Some((t, d));
    }
    None
}

fn pair_crossing(l1: usize, l2: usize) -> f64 {
    let mut values = Vec::new();
    for (sa, sb) in [(11u64, 22u64), (31, 42), (51, 62)] {
        if let Some(tc) = crossing(LatticeKind::Sc, l1, l2, sa, sb) {
            values.push(tc);
        }
    }
    values.iter().sum::<f64>() / values.len() as f64
}

#[test]
fn binder_crossing_reproduces_mc_curie_temperature() {
    // Chen, Ferrenberg & Landau (1993): k_B T_C / J = 1.44299(14).
    // Finite lattices and finite statistics give a noisy crossing; use
    // several small sizes and independent seeds and report the spread.
    let pairs = [(4usize, 6usize), (6, 8), (4, 8)];
    let values: Vec<f64> = pairs.iter().map(|&(a, b)| pair_crossing(a, b)).collect();
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    println!(
        "Binder crossings L=4/6 {:.4}, L=6/8 {:.4}, L=4/8 {:.4}, mean {:.4}",
        values[0], values[1], values[2], mean
    );
    let spread = values.iter().cloned().fold(f64::MIN, f64::max)
        - values.iter().cloned().fold(f64::MAX, f64::min);
    assert!(
        (mean - 1.442_99).abs() < 0.08,
        "Binder crossing {mean:.4} (spread {spread:.3}) vs published 1.44299"
    );
    // The mean-field value 2.0 is the wrong answer by a large margin.
    assert!(mean_field_tc(LatticeKind::Sc, 1.0, 1.0) / 1.442_99 > 1.3);
}

#[test]
fn spin_wave_dispersion_matches_llg_chain() {
    let n = 64;
    let a = 1.0;
    let dt = 0.02;
    for m in [8i32, 16, 24] {
        let k = 2.0 * std::f64::consts::PI * m as f64 / (n as f64 * a);
        // Linear spin-wave theory: omega = 2 J S (z - sum cos k.delta).
        // With S = 1/2 and gamma = 2S/hbar = 1 the LLG (unit spins) has the
        // same dispersion.
        let expected = spin_wave_omega(1.0, 0.5, &[[a, 0.0, 0.0], [-a, 0.0, 0.0]], [k, 0.0, 0.0]);
        let mut chain = SpinChain::new(n, a, 1.0, 0.0, 1.0, 0.0);
        chain.seed_spin_wave(k, 0.01);
        let mut sig = Vec::new();
        for step in 0..6000 {
            chain.step(dt);
            if step % 5 == 0 {
                sig.push(chain.spins[0][0]);
            }
        }
        let measured = measure_frequency(&sig, dt * 5.0, expected * 0.9, expected * 1.1, 4000);
        assert!(
            rel(measured, expected) < 2e-3,
            "k={k:.4} spin wave {expected:.6}, LLG {measured:.6}"
        );
    }
}

#[test]
fn discrete_wall_fold_and_its_lattice_correction() {
    // delta/a = sqrt(J/(2K)) for the J=1, a=1 chain.
    let cases = [(1.0_f64, 6e-2), (2.0, 3e-2), (5.0, 2e-2), (10.0, 2e-3)];
    let mut errors = Vec::new();
    for (target, tol) in cases {
        let k = 1.0 / (2.0 * target * target);
        let w = discrete_wall_chain(1.0, k, 1.0, 201, 200_000);
        let e_err = rel(w.sigma_measured, w.sigma_continuum);
        let d_err = rel(w.delta_measured, w.delta_continuum);
        assert!(e_err < tol, "delta/a={target}: energy error {e_err:.3e}");
        assert!(d_err < tol, "delta/a={target}: width error {d_err:.3e}");
        errors.push((target, e_err));
    }
    // The continuum fold becomes exact as delta/a -> infinity: the error
    // must fall monotonically with the wall width.
    for pair in errors.windows(2) {
        assert!(
            pair[1].1 < pair[0].1,
            "delta/a={} error {:.3e} should beat delta/a={} error {:.3e}",
            pair[1].0,
            pair[1].1,
            pair[0].0,
            pair[0].1
        );
    }
}
