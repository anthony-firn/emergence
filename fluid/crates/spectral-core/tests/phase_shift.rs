//! Phase-shift dealiasing with spherical truncation.

use spectral_core::grid::wavenumber;
use spectral_core::{init, Grid, Solver};

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)
}

#[test]
fn spherical_grid_keeps_exact_solutions_exact() {
    let g = Grid::new_spherical(32);
    let u0 = init::abc(&g, 1.0, 0.7, 0.4, 2.0);
    let mut s = Solver::new(g, 0.05, u0.clone());
    for _ in 0..100 {
        s.step(0.01);
    }
    let decay = (-0.05f64 * 4.0 * s.t).exp();
    for c in 0..3 {
        let exact: Vec<f64> = u0[c].iter().map(|v| v * decay).collect();
        assert!(max_abs_diff(&s.velocity_phys(c), &exact) < 1e-12);
    }
}

#[test]
fn spherical_retains_more_modes() {
    let (c, s) = (Grid::new(48), Grid::new_spherical(48));
    let count = |g: &Grid| (0..g.spec_len()).filter(|&i| g.kept(i)).count();
    let ratio = count(&s) as f64 / count(&c) as f64;
    eprintln!("kmax cube {} sphere {}; retained modes ratio {ratio:.2}", c.kmax(), s.kmax());
    assert!(s.kmax() as f64 > 1.35 * c.kmax() as f64);
}

/// Box-mean L² distance between a coarse and a fine spectral state.
fn distance(gc: &Grid, uc: &[Vec<num_complex::Complex64>; 3], gf: &Grid, uf: &[Vec<num_complex::Complex64>; 3]) -> f64 {
    let mut map = vec![None; gf.spec_len()];
    for idx in 0..gc.spec_len() {
        if !gc.kept(idx) {
            continue;
        }
        let (i, j, l) = gc.spec_ijl(idx);
        let w = |k: i64| k.rem_euclid(gf.n as i64) as usize;
        map[(w(wavenumber(i, gc.n)) * gf.n + w(wavenumber(j, gc.n))) * gf.nh + l] = Some(idx);
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
fn spherical_24_is_as_accurate_as_cube_32_on_taylor_green() {
    let (re, dt, t_end): (f64, f64, f64) = (400.0, 0.01, 2.5);
    let steps = (t_end / dt).round() as usize;
    let run = |g: Grid| {
        let mut s = Solver::new(g.clone(), 1.0 / re, init::taylor_green(&g));
        for _ in 0..steps {
            s.step(dt);
        }
        (g, s.u)
    };
    let (gf, uf) = run(Grid::new(96));
    let (gc, uc) = run(Grid::new(32));
    let (gs, us) = run(Grid::new_spherical(24));
    let (gc24, uc24) = run(Grid::new(24));
    let (ec, es, ec24) = (distance(&gc, &uc, &gf, &uf), distance(&gs, &us, &gf, &uf), distance(&gc24, &uc24, &gf, &uf));
    eprintln!("t={t_end}: error cube 24^3 {ec24:.3e}, spherical 24^3 {es:.3e}, cube 32^3 {ec:.3e}");
    assert!(es < 1.25 * ec, "spherical 24 ({es:e}) should match cube 32 ({ec:e})");
    assert!(es < 0.6 * ec24, "spherical 24 ({es:e}) should beat cube 24 ({ec24:e})");
}
