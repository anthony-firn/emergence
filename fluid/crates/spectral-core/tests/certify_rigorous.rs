//! Rigorous certificate: never below the true error; regularity certified.

use spectral_core::certify::{Certifier, RigorousCertifier};
use spectral_core::grid::wavenumber;
use spectral_core::{init, Grid, Solver};

fn distance(gc: &Grid, uc: &[Vec<num_complex::Complex64>; 3], gf: &Grid, uf: &[Vec<num_complex::Complex64>; 3], grad: bool) -> f64 {
    let mut map = vec![None; gf.spec_len()];
    for idx in 0..gc.spec_len() {
        let (i, j, l) = gc.spec_ijl(idx);
        let w = |k: i64| k.rem_euclid(gf.n as i64) as usize;
        map[(w(wavenumber(i, gc.n)) * gf.n + w(wavenumber(j, gc.n))) * gf.nh + l] = Some(idx);
    }
    let mut e2 = 0.0;
    for fi in 0..gf.spec_len() {
        let k = gf.k(fi);
        let wk = if grad { k[0] * k[0] + k[1] * k[1] + k[2] * k[2] } else { 1.0 };
        for c in 0..3 {
            let coarse = map[fi].map(|ci| uc[c][ci]).unwrap_or_default();
            e2 += gf.weight(fi) * wk * (uf[c][fi] - coarse).norm_sqr();
        }
    }
    e2.sqrt()
}

#[test]
fn abc_flow_is_certified_regular() {
    let g = Grid::new(16);
    let mut s = Solver::new(g.clone(), 0.1, init::abc(&g, 1.0, 0.7, 0.4, 1.0));
    let mut rc = RigorousCertifier::new(&g, 0.0, 0.0);
    for _ in 0..50 {
        s.step_certified_rigorous(0.02, &mut rc);
    }
    let st = rc.history.last().unwrap();
    eprintln!("ABC t={:.2}: L2 bound {:.2e}, H1 bound {:?}, mu sup {:.3} (sampled {:.3})", st.t, st.l2_bound, st.h1_bound, st.mu_sup, st.mu_sampled);
    assert!(st.l2_bound < 1e-9);
    assert!(st.h1_bound.unwrap() < 1e-8);
}

/// Runs Taylor–Green (Re = 200) at resolution n with the rigorous
/// certificate against a 64³ reference; checks every bound while active
/// and returns the regularity horizon (last time smoothness is certified).
fn certified_run(n: usize, t_end: f64) -> f64 {
    let (re, dt): (f64, f64) = (200.0, 0.01);
    let gc = Grid::new(n);
    let mut coarse = Solver::new(gc.clone(), 1.0 / re, init::taylor_green(&gc));
    let mut coarse2 = Solver::new(gc.clone(), 1.0 / re, init::taylor_green(&gc));
    let gf = Grid::new(64);
    let mut fine = Solver::new(gf.clone(), 1.0 / re, init::taylor_green(&gf));
    let mut rc = RigorousCertifier::new(&gc, 0.0, 0.0);
    let mut sc = Certifier::new(&gc, 0.0);
    let steps = (t_end / dt).round() as usize;
    let mut horizon = 0.0;
    for k in 1..=steps {
        let st = coarse.step_certified_rigorous(dt, &mut rc);
        let ss = coarse2.step_certified(dt, &mut sc);
        for _ in 0..4 {
            fine.step(dt / 4.0);
        }
        if st.h1_bound.is_some() {
            horizon = st.t;
        }
        if k % 25 == 0 {
            let (e0, e1) = (distance(&gc, &coarse.u, &gf, &fine.u, false), distance(&gc, &coarse.u, &gf, &fine.u, true));
            eprintln!(
                "N={n} t={:.2}: L2 true {e0:.3e} rigorous {:.3e} (sampled {:.3e}) | H1 true {e1:.3e} rigorous {:?} | mu sup {:.2} (sampled {:.2})",
                st.t, st.l2_bound, ss.bound, st.h1_bound.map(|b| format!("{b:.3e}")), st.mu_sup, st.mu_sampled
            );
            assert!(st.l2_bound >= e0);
            if let Some(b) = st.h1_bound {
                assert!(b >= e1);
            }
        }
    }
    eprintln!("N={n}: regularity certified up to t = {horizon:.2}");
    horizon
}

#[test]
fn rigorous_bounds_hold_and_resolution_extends_certified_regularity() {
    let h16 = certified_run(16, 1.0);
    let h32 = certified_run(32, 1.0);
    assert!(h16 >= 0.25, "N=16 horizon {h16}");
    assert!(h32 > h16, "better resolution must certify regularity longer: {h32} vs {h16}");
}
