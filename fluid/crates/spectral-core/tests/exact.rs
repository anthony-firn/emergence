//! Validation against exact results. These must pass before any speed claim.

use spectral_core::{fft::Fft3, init, Grid, Solver};

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)
}

/// Deterministic pseudo-random field (xorshift) for transform tests.
fn noise(len: usize, seed: u64) -> Vec<f64> {
    let mut s = seed;
    (0..len)
        .map(|_| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s as f64 / u64::MAX as f64) - 0.5
        })
        .collect()
}

#[test]
fn fft_roundtrip_is_exact() {
    for n in [8, 16, 24, 32] {
        let g = Grid::new(n);
        let fft = Fft3::new(&g);
        let x = noise(g.real_len(), 0x9E3779B97F4A7C15 ^ n as u64);
        let mut spec = vec![Default::default(); g.spec_len()];
        let mut back = vec![0.0; g.real_len()];
        fft.forward(&x, &mut spec);
        fft.inverse_inplace(&mut spec, &mut back);
        let err = max_abs_diff(&x, &back);
        assert!(err < 1e-14, "n={n} roundtrip err {err:e}");
    }
}

#[test]
fn fft_coefficients_match_analytic() {
    // f = cos(3x) sin(2y) cos(z): coefficients ±1/8 at (±3, ±2, ±1) with sign pattern.
    let g = Grid::new(16);
    let fft = Fft3::new(&g);
    let f = init::from_fn(&g, |x, y, z| [(3.0 * x).cos() * (2.0 * y).sin() * z.cos(), 0.0, 0.0]);
    let mut spec = vec![Default::default(); g.spec_len()];
    fft.forward(&f[0], &mut spec);
    let mut energy = 0.0;
    for idx in 0..g.spec_len() {
        energy += g.weight(idx) * spec[idx].norm_sqr();
        let k = g.k(idx);
        let on = k[0].abs() == 3.0 && k[1].abs() == 2.0 && k[2] == 1.0;
        if !on {
            assert!(spec[idx].norm() < 1e-15, "spurious mode {k:?}");
        } else {
            assert!((spec[idx].norm() - 0.125).abs() < 1e-15);
        }
    }
    // Parseval: mean of f² = 1/8
    assert!((energy - 0.125).abs() < 1e-15);
}

#[test]
fn abc_beltrami_decays_exactly() {
    // curl u = k u  =>  u × ω = 0  =>  u(t) = u(0) e^{-ν k² t}
    let (n, nu, k) = (32, 0.05, 2.0);
    let g = Grid::new(n);
    let u0 = init::abc(&g, 1.0, 0.7, 0.4, k);
    let mut s = Solver::new(g, nu, u0.clone());
    for _ in 0..100 {
        s.step(0.01);
    }
    let decay = (-nu * k * k * s.t).exp();
    for c in 0..3 {
        let exact: Vec<f64> = u0[c].iter().map(|v| v * decay).collect();
        let err = max_abs_diff(&s.velocity_phys(c), &exact);
        assert!(err < 1e-12, "ABC component {c} err {err:e}");
    }
}

#[test]
fn taylor_green_2d_decays_exactly() {
    // Nonlinear term is a pure gradient: only correct if the projection is right.
    let (n, nu) = (32, 0.1);
    let g = Grid::new(n);
    let u0 = init::taylor_green_2d(&g);
    let mut s = Solver::new(g, nu, u0.clone());
    for _ in 0..100 {
        s.step(0.01);
    }
    let decay = (-2.0 * nu * s.t).exp();
    for c in 0..3 {
        let exact: Vec<f64> = u0[c].iter().map(|v| v * decay).collect();
        let err = max_abs_diff(&s.velocity_phys(c), &exact);
        assert!(err < 1e-12, "TG2D component {c} err {err:e}");
    }
}

#[test]
fn inviscid_energy_is_conserved_and_flow_stays_solenoidal() {
    let g = Grid::new(32);
    let mut s = Solver::new(g, 0.0, init::taylor_green(&Grid::new(32)));
    let e0 = s.diagnostics().energy;
    assert!((e0 - 0.125).abs() < 1e-14, "TG initial energy {e0}");
    for _ in 0..100 {
        s.step(0.01);
    }
    let d = s.diagnostics();
    let drift = ((d.energy - e0) / e0).abs();
    assert!(drift < 1e-9, "inviscid energy drift {drift:e}");
    assert!(d.max_div < 1e-14, "divergence {:e}", d.max_div);
    // something nontrivial happened
    assert!(d.enstrophy > 0.375 * 1.05);
}

#[test]
fn time_integration_is_fourth_order() {
    let run = |dt: f64| {
        let g = Grid::new(16);
        let mut s = Solver::new(g.clone(), 0.01, init::taylor_green(&g));
        let steps = (0.8 / dt).round() as usize;
        for _ in 0..steps {
            s.step(dt);
        }
        (0..3).map(|c| s.velocity_phys(c)).collect::<Vec<_>>()
    };
    let reference = run(0.0025);
    let errs: Vec<f64> = [0.04, 0.02, 0.01]
        .iter()
        .map(|&dt| {
            let u = run(dt);
            (0..3).map(|c| max_abs_diff(&u[c], &reference[c])).fold(0.0, f64::max)
        })
        .collect();
    for w in errs.windows(2) {
        let order = (w[0] / w[1]).log2();
        assert!(order > 3.7, "observed order {order:.2}, errors {errs:?}");
    }
}

#[test]
fn fft_roundtrip_odd_sizes() {
    for n in [9, 21, 43, 63] {
        let g = Grid::new_any(n);
        let fft = Fft3::new(&g);
        let x = noise(g.real_len(), 0xA5A5 ^ n as u64);
        let mut spec = vec![Default::default(); g.spec_len()];
        let mut back = vec![0.0; g.real_len()];
        fft.forward(&x, &mut spec);
        let energy: f64 = (0..g.spec_len()).map(|i| g.weight(i) * spec[i].norm_sqr()).sum();
        let mean_sq: f64 = x.iter().map(|v| v * v).sum::<f64>() / x.len() as f64;
        fft.inverse_inplace(&mut spec, &mut back);
        let err = max_abs_diff(&x, &back);
        assert!(err < 1e-14, "n={n} roundtrip err {err:e}");
        assert!((energy - mean_sq).abs() < 1e-14, "n={n} Parseval {energy} vs {mean_sq}");
    }
}
