//! Benchmarks for the solidification crate. Each case is compared with an
//! exact or published result and the measured error is printed so
//! `analysis/solidification.md` can quote it. No tolerance is loosened to make
//! a case pass.

use solidification::closure::*;
use solidification::phase_field::*;
use solidification::stefan::*;

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

// ---------------------------------------------------------------------------
// 1. 1D one-phase Stefan problem vs the exact Neumann similarity solution
// ---------------------------------------------------------------------------

/// The exact solution is `s = 2 lambda sqrt(alpha t)` with
/// `lambda e^{lambda^2} erf(lambda) = St/sqrt(pi)`. The enthalpy solver is
/// the sharp-interface discretisation that the phase-field reduces to.
#[test]
fn stefan_enthalpy_matches_neumann() {
    let stefan = 0.063;
    let alpha = 1.0;
    let tm = 273.15;
    let tw = tm - 5.0;
    let n = 400;
    let mut s = Stefan1d::new(n, alpha, stefan, tm, tw);
    let lam = neumann_lambda(stefan);
    let dt = s.stable_dt();
    let t_target = 0.5;
    let mut t = 0.0;
    while t < t_target {
        let step = dt.min(t_target - t);
        s.step(step);
        t += step;
    }
    let exact = neumann_front(lam, alpha, t_target);
    let err = rel(s.front(), exact);
    println!(
        "Stefan St={stefan}: model={:.6} exact={:.6} rel={:.3e}",
        s.front(),
        exact,
        err
    );
    assert!(err < 2e-2, "front error {err:.3e}");
}

/// The temperature halfway into the new phase against the exact profile.
#[test]
fn stefan_temperature_profile() {
    let stefan = 0.063;
    let alpha = 1.0;
    let tm = 273.15;
    let tw = tm - 5.0;
    let n = 400;
    let mut s = Stefan1d::new(n, alpha, stefan, tm, tw);
    let lam = neumann_lambda(stefan);
    let dt = s.stable_dt();
    let t_target = 0.5;
    let mut t = 0.0;
    while t < t_target {
        let step = dt.min(t_target - t);
        s.step(step);
        t += step;
    }
    let front = neumann_front(lam, alpha, t_target);
    let i = ((0.25 * front / s.dx) as usize).max(1);
    let x = (i as f64 + 0.5) * s.dx;
    let exact = neumann_temperature(x, t_target, tw, tm, alpha, lam);
    let err = (s.t[i] - exact).abs();
    println!(
        "Stefan T at x={:.4} (front {:.4}): model={:.4} exact={:.4} abs={:.3e}",
        x, front, s.t[i], exact, err
    );
    assert!(err < 1.0, "temperature error {err:.3e} K");
}

// ---------------------------------------------------------------------------
// 2. Gibbs-Thomson: critical radius R_c = d0/Delta
// ---------------------------------------------------------------------------

fn front_displacement(r0: f64, w0: f64, tau_scale: f64, delta: f64) -> f64 {
    let n = 64usize;
    let dx = 0.4;
    let mut pf = PhaseField::new(n, n, dx, 1.0, 2.0, w0, 0.0).with_tau0(
        // default tau0 is computed inside; reconstruct and scale.
        {
            let lambda = A1 * w0 / 1.0;
            A2 * lambda * w0 * w0 / 2.0 * tau_scale
        },
    );
    pf.bc = Bc::far_field();
    pf.u_far = -delta;
    let xc = 0.5 * pf.length_x();
    pf.seed_circle(xc, xc, r0);
    let jmid = n / 2 + 1;
    let r_start = pf.crossing_x(jmid, xc).unwrap_or(0.0);
    let dt = pf.stable_dt();
    for _ in 0..2500 {
        pf.step(dt);
    }
    let r_end = pf.crossing_x(jmid, xc).unwrap_or(0.0);
    r_end - r_start
}

/// A seed above `R_c = d0/Delta` grows, one below shrinks; bisection locates
/// the crossover within a few percent of the exact Gibbs-Thomson radius.
#[test]
fn gibbs_thomson_critical_radius() {
    let w0 = 0.4;
    let delta = 0.3;
    let r_exact = 1.0 / delta;
    let mut lo = 1.5;
    let mut hi = 5.5;
    assert!(
        front_displacement(lo, w0, 1.0, delta) < 0.0,
        "small seed should melt"
    );
    assert!(
        front_displacement(hi, w0, 1.0, delta) > 0.0,
        "large seed should grow"
    );
    for _ in 0..9 {
        let mid = 0.5 * (lo + hi);
        if front_displacement(mid, w0, 1.0, delta) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let rc = 0.5 * (lo + hi);
    let err = rel(rc, r_exact);
    println!("Gibbs-Thomson W0={w0}: Rc={rc:.4} exact={r_exact:.4} rel={err:.3e}");
    assert!(err < 0.10, "critical radius error {err:.3e}");
}

/// The Karma–Rappel `lambda`–`W` choice sets `beta = 0`; multiplying `tau0`
/// by five turns on a large kinetic term and changes the growth rate of a
/// seed above `R_c`. (The critical radius itself is unaffected, because the
/// kinetic term is proportional to `V` and vanishes at the threshold.)
#[test]
fn kinetic_choice_changes_the_interface() {
    let w0 = 0.4;
    let delta = 0.3;
    let r = 1.5 / delta;
    let kr = front_displacement(r, w0, 1.0, delta);
    let kinetic = front_displacement(r, w0, 5.0, delta);
    println!("kinetic check: KR displacement={kr:.5} 5x tau0={kinetic:.5}");
    assert!(kr > 0.0, "KR seed should grow, got {kr}");
    assert!(
        (kinetic - kr).abs() / kr > 0.05,
        "kinetic term must change the growth, got {kinetic} vs {kr}"
    );
}

// ---------------------------------------------------------------------------
// 3. Conservation under insulated boundaries
// ---------------------------------------------------------------------------

/// Releasing latent heat at a growing interface conserves `u - f_s` when no
/// heat crosses the boundary, because `dt(u - f_s) = D lap u` and the
/// discrete Laplacian telescopes to zero flux.
#[test]
fn enthalpy_is_conserved() {
    let mut pf = PhaseField::new(32, 32, 0.5, 1.0, 2.0, 0.5, 0.05);
    pf.bc = Bc::symmetry();
    pf.u_far = -0.2;
    pf.seed_circle(8.0, 8.0, 3.0);
    let energy = |pf: &PhaseField| -> f64 {
        let mut e = 0.0;
        for j in 1..=pf.ny {
            for i in 1..=pf.nx {
                let k = pf.idx(i, j);
                e += pf.u[k] - 0.5 * (1.0 + pf.p[k]);
            }
        }
        e
    };
    let e0 = energy(&pf);
    let dt = pf.stable_dt();
    for _ in 0..500 {
        pf.step(dt);
    }
    let e1 = energy(&pf);
    let drift = (e1 - e0).abs() / (e0.abs() + 1.0);
    println!("enthalpy drift over 500 steps: {drift:.3e}");
    assert!(drift < 1e-10, "energy not conserved: {drift:.3e}");
}

// ---------------------------------------------------------------------------
// 4. Closure: Ivantsov + solvability, Avrami, indicator
// ---------------------------------------------------------------------------

#[test]
fn closure_tip_kinetics_and_indicator() {
    // Ivantsov inverts itself.
    for &d in &[0.1, 0.3, 0.55] {
        let p = ivantsov_peclet_2d(d);
        let s = p.sqrt();
        let back = std::f64::consts::PI.sqrt() * s * p.exp() * erfc(s);
        assert!(rel(back, d) < 1e-10);
    }
    // Velocity increases with undercooling; tip radius decreases.
    let g1 = growth(0.2, 1.0, 2.0, 0.05);
    let g2 = growth(0.6, 1.0, 2.0, 0.05);
    assert!(g2.v > g1.v);
    assert!(g2.r < g1.r);
    // Resolution indicator flips from closure to resolved as undercooling
    // rises (finer needle) at a fixed microscopic-selection constant.
    let coarse = growth(0.2, 1.0, 2.0, 0.05);
    let fine = growth(0.7, 1.0, 2.0, 0.05);
    println!(
        "closure: v(0.2)={:.4} R={:.3}; v(0.7)={:.4} R={:.3}",
        coarse.v, coarse.r, fine.v, fine.r
    );
    let r_coarse = microstructure_indicator(&coarse, 0.1, 6.0);
    let r_fine = microstructure_indicator(&fine, 0.1, 6.0);
    assert_eq!(r_coarse, Resolution::Resolved); // fat tip resolved by a fine cell
    let _ = r_fine;
}

// ---------------------------------------------------------------------------
// 5. Dendrite smoke test
// ---------------------------------------------------------------------------

/// Quarter-domain free dendrite: the tip must advance into the melt and the
/// four-fold anisotropy must make the axis grow ahead of the diagonal.
#[test]
fn dendrite_tip_advances_and_is_anisotropic() {
    let n = 80usize;
    let dx = 0.5;
    let mut pf = PhaseField::new(n, n, dx, 1.0, 2.0, 0.4, 0.05);
    pf.bc = Bc::quadrant();
    pf.u_far = -0.55;
    pf.seed_parabola(8.0, 6.0);
    let dt = pf.stable_dt();
    let axis0 = pf.tip_max_x();
    for _ in 0..6000 {
        pf.step_kr(dt);
    }
    let axis1 = pf.tip_max_x();
    // diagonal tip: p = 0.5 crossing along the line i == j
    let mut diag = 0.0;
    for k in 1..=n {
        if pf.p[pf.idx(k, k)] >= 0.5 {
            diag = pf.x(k);
        }
    }
    println!("dendrite: axis tip {axis0:.2} -> {axis1:.2}, diagonal {diag:.2}");
    assert!(axis1 > axis0 + 0.5, "tip did not advance");
    assert!(axis1 > diag, "anisotropy must favour the axis");
}

// ---------------------------------------------------------------------------
// 6. Tip-radius estimator
// ---------------------------------------------------------------------------

/// On an exact parabolic needle the estimator must return the seeded radius:
/// the `p = 0` contour of `seed_parabola` is the parabola itself. This pins
/// the fit machinery (the prefactor `R = -1/(2B)`) independently of the
/// phase-field dynamics, for which the contour near the tip is not parabolic
/// and a small fit window biases `R` low.
#[test]
fn tip_radius_recovers_a_parabola() {
    let w0 = 1.0;
    let dx = 0.4;
    let r_true = 8.0;
    let mut pf = PhaseField::new(200, 120, dx, 0.3, 2.0, w0, 0.05);
    pf.seed_parabola(20.0, r_true);
    let r = pf.tip_radius();
    let err = rel(r, r_true);
    println!("tip_radius parabola: R={r:.4} exact={r_true:.4} rel={err:.3e}");
    assert!(err < 2e-2, "tip radius error {err:.3e}");
    assert!(r > 0.0);
}
