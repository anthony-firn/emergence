//! The core discretisation (no immersed body) against exact solutions and
//! published data.

use std::f64::consts::PI;

use immersed::{Bc, Grid, Solver};

fn periodic() -> [Bc; 2] {
    [Bc::Periodic, Bc::Periodic]
}

/// Plane Poiseuille flow through a channel with a velocity inlet, a
/// pressure outlet and no-slip walls.
///
/// The developed discrete solution is known in closed form: interior second
/// differences of a quadratic are exact and the ghost-cell wall closure
/// requires `u(−h/2) = −u(h/2)`, so `u_j = A (y_j(1 − y_j) + h²/4)` with `A`
/// fixed by the inlet flux, `Σ h·4y_j(1 − y_j) = 4(1/6 + h²/12)` (midpoint
/// rule), i.e. `A = 4(1/6 + h²/12)/(1/6 + h²/3)`, and `dp/dx = −2νA`. The
/// solver must reproduce it to round-off; its error against the exact
/// parabola is then `h² − 3h³ + …` (second order, approached from below).
#[test]
fn channel_poiseuille_reproduces_the_exact_discrete_solution() {
    // Long enough that the inlet mismatch (O(h²)) has decayed like
    // e^{−4.2 x/H} (slowest Stokes eigenmode) to round-off where we look.
    let (len, nu) = (8.0, 0.1);
    let exact = |y: f64| 4.0 * y * (1.0 - y);
    let mut errs = Vec::new();
    for ny in [8usize, 16, 32] {
        let h = 1.0 / ny as f64;
        let grid = Grid::new([8 * ny, ny, 1], h, [0.0; 3]);
        let bcs = [[Bc::inlet(move |x, _| [exact(x[1]), 0.0, 0.0]), Bc::Outlet], [Bc::wall(), Bc::wall()], periodic()];
        let mut s = Solver::new(grid, nu, 1.0, bcs, vec![]);
        s.initialize(|x| [exact(x[1]), 0.0, 0.0]);
        let (steps, rate) = s.run_until_steady(0.8, 200.0, 1e-11);
        assert!(rate < 1e-11, "ny {ny}: not steady after {steps} steps (rate {rate:e})");
        let a = 4.0 * (1.0 / 6.0 + h * h / 12.0) / (1.0 / 6.0 + h * h / 3.0);
        let (mut disc, mut err) = (0.0f64, 0.0f64);
        for j in 0..ny {
            let y = (j as f64 + 0.5) * h;
            let u = s.velocity_at([0.75 * len, y, 0.5 * h])[0];
            disc = disc.max((u - a * (y * (1.0 - y) + h * h / 4.0)).abs());
            err = err.max((u - exact(y)).abs());
        }
        let dpdx = s.pressure_at([6.0, 0.5, 0.0]) - s.pressure_at([5.0, 0.5, 0.0]);
        println!(
            "ny {ny:3}: {steps:5} steps; vs discrete closed form {disc:.1e}, dp/dx {:.1e}; vs exact parabola {err:.3e}; max div {:.1e}",
            (dpdx + 2.0 * nu * a).abs(),
            s.max_divergence()
        );
        assert!(disc < 1e-9 && (dpdx + 2.0 * nu * a).abs() < 1e-9);
        assert!(s.max_divergence() < 1e-10);
        errs.push(err);
    }
    // Ratios approach 4 from below (h² − 3h³): 3.13, 3.59, then 3.80.
    assert!((errs[1] / errs[2]).log2() > 1.8);
}

/// Decaying 2D Taylor–Green vortex (exact solution), periodic box.
#[test]
fn taylor_green_vortex_converges_at_second_order() {
    let nu = 0.1;
    let t_end = 0.5;
    let decay = |t: f64| (-2.0 * nu * t).exp();
    let mut errs = Vec::new();
    for n in [16usize, 32, 64] {
        let h = 2.0 * PI / n as f64;
        let grid = Grid::new([n, n, 1], h, [0.0; 3]);
        let mut s = Solver::new(grid, nu, 1.0, [periodic(), periodic(), periodic()], vec![]);
        s.initialize(|x| [-x[0].cos() * x[1].sin(), x[0].sin() * x[1].cos(), 0.0]);
        s.initialize_pressure(|x| -0.25 * ((2.0 * x[0]).cos() + (2.0 * x[1]).cos()));
        while s.t < t_end - 1e-12 {
            let dt = s.stable_dt(0.2).min(t_end - s.t);
            s.step(dt);
        }
        let mut err: f64 = 0.0;
        for j in 0..n {
            for i in 0..n {
                let x = [i as f64 * h, (j as f64 + 0.5) * h, 0.5 * h];
                let u = s.vel[0].get(i as isize, j as isize, 0);
                err = err.max((u + x[0].cos() * x[1].sin() * decay(s.t)).abs());
            }
        }
        println!("n {n:3}: max |u − exact| {err:.3e}");
        errs.push(err);
    }
    for w in errs.windows(2) {
        let order = (w[0] / w[1]).log2();
        assert!(order > 1.9, "order {order}");
    }
}

/// Lid-driven cavity, Re = 100, against Ghia, Ghia & Shin (J. Comput.
/// Phys. 48, 1982, Tables I–II).
#[test]
fn lid_driven_cavity_matches_ghia_at_re_100() {
    let y_u = [0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547];
    let u_g = [0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717];
    let x_v = [0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625];
    let v_g = [-0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233];
    let n = 64;
    let h = 1.0 / n as f64;
    let grid = Grid::new([n, n, 1], h, [0.0; 3]);
    let bcs = [[Bc::wall(), Bc::wall()], [Bc::wall(), Bc::moving_wall([1.0, 0.0, 0.0])], periodic()];
    let mut s = Solver::new(grid, 0.01, 1.0, bcs, vec![]);
    let (steps, rate) = s.run_until_steady(0.8, 60.0, 1e-7);
    assert!(rate < 1e-7, "not steady: rate {rate:e} after {steps} steps");
    let mut worst: f64 = 0.0;
    for (y, ug) in y_u.iter().zip(&u_g) {
        worst = worst.max((s.velocity_at([0.5, *y, 0.5 * h])[0] - ug).abs());
    }
    for (x, vg) in x_v.iter().zip(&v_g) {
        worst = worst.max((s.velocity_at([*x, 0.5, 0.5 * h])[1] - vg).abs());
    }
    println!("cavity Re=100, 64²: {steps} steps to steady, worst |Δ| vs Ghia {worst:.4} (lid speed 1)");
    assert!(worst < 0.015, "worst deviation from Ghia {worst}");
}
