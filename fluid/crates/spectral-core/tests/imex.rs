//! IMEX scheme with implicit penalization.

use spectral_core::{init, Grid, Penalization, Solver};

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)
}

#[test]
fn imex_is_third_order_without_objects() {
    let run = |dt: f64, imex: bool| {
        let g = Grid::new(16);
        let mut s = Solver::new(g.clone(), 0.01, init::taylor_green(&g));
        let steps = (0.8 / dt).round() as usize;
        for _ in 0..steps {
            if imex {
                s.step_imex(dt)
            } else {
                s.step_ifrk4(dt)
            }
        }
        (0..3).map(|c| s.velocity_phys(c)).collect::<Vec<_>>()
    };
    let reference = run(0.001, false);
    let errs: Vec<f64> = [0.08, 0.04, 0.02]
        .iter()
        .map(|&dt| {
            let u = run(dt, true);
            (0..3).map(|c| max_abs_diff(&u[c], &reference[c])).fold(0.0, f64::max)
        })
        .collect();
    let orders: Vec<f64> = errs.windows(2).map(|w| (w[0] / w[1]).log2()).collect();
    eprintln!("IMEX errors {errs:?}, observed orders {orders:?}");
    assert!(orders.iter().all(|&o| o > 2.7), "orders {orders:?}");
}

/// Steady Stokes flow past a sphere in a periodic array; returns (Cd, CG iterations).
fn stokes_cd(n: usize, eta: f64, dt: f64, steps: usize, sharp: bool) -> (f64, [usize; 3]) {
    let two_pi = 2.0 * std::f64::consts::PI;
    let (d, u, nu) = (two_pi / 8.0, 1e-3, 1.0);
    let g = Grid::new(n);
    let center = [0.5 * two_pi; 3];
    let pen = if sharp {
        Penalization::sphere_sharp(&g, center, d, eta)
    } else {
        Penalization::sphere(&g, center, d, g.dx(), eta)
    }
    .implicit();
    let mut s = Solver::new(g.clone(), nu, init::from_fn(&g, |_, _, _| [0.0; 3]));
    s.penalization = Some(pen);
    s.set_mean_flow([u, 0.0, 0.0]);
    for _ in 0..steps {
        s.step(dt);
    }
    let f = s.solid_force().unwrap();
    let q = 0.5 * u * u * std::f64::consts::PI * d * d / 4.0;
    (f[0] / q, s.last_cg_iters().unwrap())
}

#[test]
fn implicit_penalization_converges_in_eta_and_is_dt_independent() {
    // Re = U D / ν = 7.854e-4; exact periodic-array value Cd = 24/Re · 1.2139 = 37,094
    let exact = 24.0 / 7.853_981_6e-4 * 1.2139;
    let (cd6, it6) = stokes_cd(32, 1e-6, 1.0, 12, true);
    let (cd8, it8) = stokes_cd(32, 1e-8, 1.0, 12, true);
    let (cd8b, _) = stokes_cd(32, 1e-8, 0.25, 48, true);
    let (sm6, _) = stokes_cd(32, 1e-6, 1.0, 12, false);
    let (sm8, _) = stokes_cd(32, 1e-8, 1.0, 12, false);
    eprintln!("smooth mask (for contrast): Cd(1e-6)={sm6:.1} Cd(1e-8)={sm8:.1}");
    eprintln!("N=32: Cd(eta=1e-6)={cd6:.2} iters {it6:?}; Cd(eta=1e-8)={cd8:.2} iters {it8:?}; dt=0.25: {cd8b:.2}; exact {exact:.1}");
    // sharp mask: remaining η dependence is the √(νη) Brinkman layer
    // (√(νη)/D = 0.13% at η = 1e-6), so 1e-6 → 1e-8 moves Cd by < 0.6%
    assert!(((cd6 - cd8) / cd8).abs() < 6e-3, "eta dependence {cd6} vs {cd8}");
    // the smooth mask is far worse (front creeps like 1/ln(1/η))
    assert!(((sm6 - sm8) / sm8).abs() > 5.0 * ((cd6 - cd8) / cd8).abs());
    // steady state independent of the time step (stiffly accurate IMEX)
    assert!(((cd8 - cd8b) / cd8).abs() < 1e-4, "dt dependence {cd8} vs {cd8b}");
    // 4 cells per diameter: only coarse agreement expected here
    assert!(((cd8 - exact) / exact).abs() < 0.25, "Cd {cd8} vs exact {exact}");
}
