//! The certificate must never be smaller than the true error.

use spectral_core::certify::{sym3_min_eig, Certifier};
use spectral_core::grid::wavenumber;
use spectral_core::{init, Grid, Solver};

#[test]
fn closed_form_min_eigenvalue() {
    // diagonal, and a rotated matrix with known eigenvalues {-2, 1, 4}
    assert!((sym3_min_eig(3.0, -1.0, 2.0, 0.0, 0.0, 0.0) + 1.0).abs() < 1e-14);
    let (c, s) = (0.6f64, 0.8f64); // rotation in the x–y plane
    let (l1, l2, l3) = (-2.0, 1.0, 4.0);
    let a = c * c * l1 + s * s * l2;
    let b = s * s * l1 + c * c * l2;
    let d = c * s * (l1 - l2);
    assert!((sym3_min_eig(a, b, l3, d, 0.0, 0.0) - l1).abs() < 1e-12);
}

#[test]
fn abc_flow_certificate_is_tiny() {
    // exact NS solution with no truncation: residual and bound at roundoff
    let g = Grid::new(16);
    let mut s = Solver::new(g.clone(), 0.1, init::abc(&g, 1.0, 0.7, 0.4, 1.0));
    let mut cert = Certifier::new(&g, 0.0);
    for _ in 0..50 {
        s.step_certified(0.02, &mut cert);
    }
    let last = cert.history.last().unwrap();
    eprintln!("ABC: bound {:.2e}, resid trunc {:.1e} time {:.1e}, mu {:.3}", last.bound, last.resid_trunc, last.resid_time, last.mu);
    assert!(last.bound < 1e-9);
}

/// Box-mean L² distance between a coarse and a fine spectral state.
fn spectral_distance(gc: &Grid, uc: &[Vec<num_complex::Complex64>; 3], gf: &Grid, uf: &[Vec<num_complex::Complex64>; 3]) -> f64 {
    let mut map = vec![None; gf.spec_len()];
    for idx in 0..gc.spec_len() {
        let (i, j, l) = gc.spec_ijl(idx);
        let w = |k: i64| k.rem_euclid(gf.n as i64) as usize;
        let fi = (w(wavenumber(i, gc.n)) * gf.n + w(wavenumber(j, gc.n))) * gf.nh + l;
        map[fi] = Some(idx);
    }
    let mut e2 = 0.0;
    for fi in 0..gf.spec_len() {
        for c in 0..3 {
            let coarse = map[fi].map(|ci| uc[c][ci]).unwrap_or_default();
            e2 += gf.weight(fi) * (uf[c][fi] - coarse).norm_sqr();
        }
    }
    e2.sqrt()
}

#[test]
fn certificate_bounds_true_error_on_underresolved_taylor_green() {
    let (re, dt, t_end): (f64, f64, f64) = (200.0, 0.01, 1.5);
    let gc = Grid::new(16);
    let mut coarse = Solver::new(gc.clone(), 1.0 / re, init::taylor_green(&gc));
    let gf = Grid::new(64);
    let mut fine = Solver::new(gf.clone(), 1.0 / re, init::taylor_green(&gf));
    let mut cert = Certifier::new(&gc, 0.0);
    let steps = (t_end / dt).round() as usize;
    let mut worst_ratio = f64::INFINITY;
    for k in 1..=steps {
        let st = coarse.step_certified(dt, &mut cert);
        for _ in 0..4 {
            fine.step(dt / 4.0);
        }
        if k % 25 == 0 {
            let err = spectral_distance(&gc, &coarse.u, &gf, &fine.u);
            eprintln!(
                "t={:.2}: true error {err:.3e}  certified bound {:.3e}  (ratio {:.1})  resid trunc {:.1e} time {:.1e} mu {:.2}",
                st.t, st.bound, st.bound / err, st.resid_trunc, st.resid_time, st.mu
            );
            assert!(st.bound >= err, "bound {} < true error {} at t={}", st.bound, err, st.t);
            worst_ratio = worst_ratio.min(st.bound / err);
        }
    }
    eprintln!("tightest bound/error ratio: {worst_ratio:.2}");
}
