//! Micromagnetic-rung benchmarks: exchange-stiffness fold, Bloch wall,
//! Walker breakdown and Kittel FMR.

use magnetism::geometry::{exchange_stiffness, LatticeKind};
use magnetism::hysteresis::Macrospin;
use magnetism::micromag::{
    anisotropy_field, exchange_length, kittel_omega, walker_field, walker_velocity, wall_delta,
    wall_sigma, Llg1d, WalkerReduced, GAMMA,
};
use magnetism::spin::measure_frequency;

fn rel(actual: f64, expected: f64) -> f64 {
    (actual - expected).abs() / expected.abs().max(1e-300)
}

#[test]
fn exchange_stiffness_fold_relations() {
    // A = J S^2 sum|delta|^2/(12 V_atom); for the conventional cubic cells
    // sum|delta|^2 = 6 a^2, giving A = J S^2 n_atoms/(2a).
    let (j, s, a) = (1.0, 0.5, 1.0);
    assert!(rel(exchange_stiffness(LatticeKind::Sc, j, s, a), 0.125) < 1e-12);
    assert!(rel(exchange_stiffness(LatticeKind::Bcc, j, s, a), 0.25) < 1e-12);
    assert!(rel(exchange_stiffness(LatticeKind::Fcc, j, s, a), 0.5) < 1e-12);
    for kind in [LatticeKind::Sc, LatticeKind::Bcc, LatticeKind::Fcc] {
        let n = kind.n_atoms() as f64;
        let closed = j * s * s * n / (2.0 * a);
        assert!(rel(exchange_stiffness(kind, j, s, a), closed) < 1e-12);
    }
}

#[test]
fn bloch_wall_matches_exact_profile_energy_and_width() {
    let (a, k, ms) = (1.0e-11, 1.0e5, 1.0e6);
    let delta = wall_delta(a, k);
    let dx = delta / 8.0;
    let n = 200;
    let mut s = Llg1d::new(n, dx, a, k, ms, 1.0);
    s.seed_wall();
    let dt = 2.0e-14;
    s.relax(dt, 40_000);

    // sigma = 4 sqrt(A K), delta = sqrt(A/K).
    let sigma_rel = rel(s.wall_energy(), wall_sigma(a, k));
    assert!(sigma_rel < 2e-3, "wall energy error {sigma_rel:.3e}");

    // Width from the maximum slope dtheta/dx = 1/delta at the centre.
    let mut best = 0.0_f64;
    for i in 1..n - 1 {
        let th0 = dot(s.m[i - 1], s.easy).clamp(-1.0, 1.0).acos();
        let th1 = dot(s.m[i + 1], s.easy).clamp(-1.0, 1.0).acos();
        best = best.max((th0 - th1).abs() / (2.0 * dx));
    }
    let delta_rel = rel(1.0 / best, delta);
    assert!(delta_rel < 2e-3, "wall width error {delta_rel:.3e}");
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[test]
fn walker_velocity_and_breakdown() {
    // Exact 1D model: H_W = alpha M_s/2, v = gamma delta H/alpha below it.
    let (a, k, ms, alpha) = (1.0e-11, 1.0e5, 1.0e6, 0.1);
    let hw = walker_field(ms, alpha);
    assert!(rel(hw, 0.5 * alpha * ms) < 1e-14);
    let dt = 1.0e-13;

    let averaged = |h: f64| -> f64 {
        let mut w = WalkerReduced::new(a, k, ms, alpha);
        for _ in 0..20_000 {
            w.step(h, dt);
        }
        let mut sum = 0.0;
        let samples = 20_000;
        for _ in 0..samples {
            w.step(h, dt);
            sum += w.velocity(h);
        }
        sum / samples as f64
    };

    for f in [0.3_f64, 0.5, 0.8, 1.0] {
        let h = f * hw;
        let v = averaged(h);
        let v_lin = walker_velocity(a, k, h, alpha);
        assert!(rel(v, v_lin) < 2e-3, "f={f}: v={v:.4e} vs {v_lin:.4e}");
    }
    // Above H_W the average velocity drops below the linear law.
    let v15 = averaged(1.5 * hw);
    let v20 = averaged(2.0 * hw);
    let vlin = gamma_helper(a, k, alpha, 1.5 * hw);
    assert!(v15 < 0.5 * vlin, "Walker breakdown: {v15:.4e} vs {vlin:.4e}");
    assert!(v20 < v15, "velocity must keep falling above breakdown");
}

fn gamma_helper(a: f64, k: f64, alpha: f64, h: f64) -> f64 {
    GAMMA * wall_delta(a, k) * h / alpha
}

#[test]
fn driven_wall_pde_at_low_field_matches_walker_velocity() {
    // Seed the wall at the predicted steady Walker azimuth so the short
    // integration measures the terminal velocity rather than the transient.
    let (a, k, ms, alpha) = (1.0e-11, 1.0e5, 1.0e6, 0.1);
    let delta = wall_delta(a, k);
    let hw = walker_field(ms, alpha);
    let h = 0.3 * hw;
    let phi = 0.5 * (2.0 * h / (alpha * ms)).asin();
    let dx = delta / 10.0;
    let n = 1200;
    let dt = 2.0e-14;
    let mut s = Llg1d::new(n, dx, a, k, ms, alpha);
    s.easy = [1.0, 0.0, 0.0];
    s.demag = vec![([0.0, 0.0, 1.0], 1.0)];
    s.h_ext = [h, 0.0, 0.0];
    let x0 = 0.5 * (n as f64 - 1.0) * dx;
    for i in 0..n {
        let th = 2.0 * (((s.x(i) - x0) / delta).exp()).atan();
        s.m[i] = [th.cos(), th.sin() * phi.cos(), th.sin() * phi.sin()];
    }
    for _ in 0..3000 {
        s.step(dt);
    }
    let q0 = s.wall_position();
    for _ in 0..5000 {
        s.step(dt);
    }
    let q1 = s.wall_position();
    let v = (q1 - q0) / (5000.0 * dt);
    let v_lin = walker_velocity(a, k, h, alpha);
    assert!(
        rel(v, v_lin) < 0.15,
        "PDE wall velocity {v:.4e} vs Walker {v_lin:.4e}"
    );
}

#[test]
fn kittel_fmr_frequency() {
    // Normalised macrospin (H_K = 1, gamma = 1, M_s = 1) with a thin-film
    // demagnetising factor N = 1: omega = sqrt((H+H_K)(H+H_K+M_s)).
    for h in [0.2_f64, 0.5, 1.0, 2.0] {
        let expected = ((h + 1.0) * (h + 2.0)).sqrt();
        let mut ms = Macrospin::new(1.0e-4);
        ms.easy = [1.0, 0.0, 0.0];
        ms.ms = 1.0;
        ms.demag = vec![([0.0, 0.0, 1.0], 1.0)];
        ms.m = [0.99995, 0.0, 0.01];
        let dt = 0.002;
        let mut sig = Vec::new();
        for st in 0..200_000 {
            ms.step([h, 0.0, 0.0], dt);
            if st % 20 == 0 {
                sig.push(ms.m[2]);
            }
        }
        let measured = measure_frequency(&sig, dt * 20.0, expected * 0.9, expected * 1.1, 4000);
        assert!(
            rel(measured, expected) < 1e-3,
            "Kittel h={h}: {measured:.6} vs {expected:.6}"
        );
    }
    // The physical Kittel formula uses gamma = 2.211e5 m/(A s).
    let omega = kittel_omega(1.0e5, 1.0e5, 8.0e5);
    assert!(omega > 0.0);
    let _ = (anisotropy_field(1.0e5, 1.0e6), exchange_length(1.0e-11, 1.0e6));
}
