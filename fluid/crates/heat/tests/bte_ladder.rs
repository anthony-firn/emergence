//! The grey-slab conductivity ladder on the emergence engine: every rung's
//! claim against the Richardson-extrapolated integral equation, the
//! reference's own claim against a finer extrapolation, and the cost of
//! choosing rungs over twelve decades of Knudsen number.
//!
//! `cargo test -p heat --test bte_ladder -- --nocapture` prints the tables
//! quoted in `analysis/emergence_engine.md`.

use emergence::Request;
use heat::bte1d::solve_gray_bte_integral;
use heat::bte_ladder::*;
use rayon::prelude::*;

#[test]
fn reference_claim_holds_against_a_finer_extrapolation() {
    // Extrapolate from (2m, 4m) and compare with the rung's (m, 2m) value.
    let rows: Vec<(f64, f64, f64)> = [0.3, 1.0, 3.0, 10.0]
        .par_iter()
        .map(|&kn| {
            let r = integral_equation(kn);
            let m = ie_cells(kn);
            let r2 = solve_gray_bte_integral(kn, 2 * m).kappa_ratio;
            let r4 = solve_gray_bte_integral(kn, 4 * m).kappa_ratio;
            let finer = r4 + (r4 - r2) / 2.41;
            (kn, (r.ratio - finer).abs(), r.estimate)
        })
        .collect();
    for (kn, err, claim) in rows {
        println!("Kn {kn:5}: |extrap(m) − extrap(2m)| {err:.2e}, claim {claim:.2e}");
        assert!(err <= claim, "Kn {kn}: {err:e} > {claim:e}");
    }
}

#[test]
fn analytic_claims_hold_on_an_independent_sweep() {
    // 48 log-spaced Kn in [0.2, 300], none of them used to fit the constants.
    let kns: Vec<f64> = (0..48).map(|i| 0.2 * 1500f64.powf((i as f64 + 0.37) / 48.0)).collect();
    let refs: Vec<KappaRatio> = kns.par_iter().map(|&kn| integral_equation(kn)).collect();
    let ladder = ladder(None, 0);
    let mut worst = vec![0.0f64; 3];
    for (kn, r) in kns.iter().zip(&refs) {
        for (rung, w) in worst.iter_mut().enumerate() {
            let emergence::Assessment::Claim(c) = ladder.rung(rung).assess(kn) else { continue };
            let v = ladder.rung(rung).eval(kn);
            let err = metric(&v, r);
            assert!(err <= c.error + r.estimate, "{} at Kn {kn}: error {err:.3e} > claim {:.3e}", ladder.rung_name(rung), c.error);
            *w = w.max((err - r.estimate).max(0.0) / c.error);
        }
    }
    for (rung, w) in worst.iter().enumerate() {
        println!("{:>20}: worst error/claim {w:.3}", ladder.rung_name(rung));
    }
}

#[test]
fn ladder_spans_twelve_decades_at_a_fraction_of_the_exact_cost() {
    let closure = closure_rung(0.3, 30.0, 1e-6);
    println!(
        "closure: {} patches, {} integral-equation samples, coverage {}",
        closure.closure().patches(),
        closure.closure().samples,
        closure.closure().coverage()
    );
    let ladder = ladder(Some(closure), 0);
    let ie = ladder.len() - 3;
    for &tol in &[1e-3, 1e-5] {
        ladder.reset_counters();
        let kns: Vec<f64> = (0..1500).map(|i| 1e-6 * 1e12f64.powf((i as f64 + 0.5) / 1500.0)).collect();
        let solved = ladder.solve_many(&kns, Request::new(tol));
        let mut cost = 0.0;
        let mut ie_everywhere = 0.0;
        for (kn, s) in kns.iter().zip(&solved) {
            assert!(s.met, "tol {tol} unmet at Kn {kn} (rung {})", ladder.rung_name(s.rung));
            cost += s.cost;
            ie_everywhere += ie_cost(*kn);
        }
        let n_ie = solved.iter().filter(|s| s.rung >= ie).count();
        println!(
            "tol {tol:e}: 1500 Kn in [1e-6, 1e6]; cost {cost:.3e} ns vs {ie_everywhere:.3e} ns for the integral equation everywhere \
             ({:.1e}×); {n_ie} integral-equation solves",
            ie_everywhere / cost
        );
        println!("{}", ladder.report());
    }
    // Audit a sample of the cheap answers against the reference.
    let audited = ladder.solve(&2.0, Request::new(1e-5));
    assert!(audited.met);
    let exact = integral_equation(2.0);
    assert!(metric(&audited.value, &exact) <= audited.claim.error + exact.estimate);
}

#[test]
fn certified_film_conductivity_spends_accuracy_where_the_heat_is() {
    use heat::effective::effective_kappa;
    use heat::materials::{silicon_branches, silicon_params};
    let t = 300.0;
    let branches = silicon_branches();
    let modes = mode_spectrum(&branches, &silicon_params(f64::INFINITY), t, 400);
    let bulk = heat::kinetic::kappa_rta(&branches, &silicon_params(f64::INFINITY), t);
    let sum: f64 = modes.iter().map(|m| m.kappa).sum();
    assert!((sum - bulk).abs() <= 1e-9 * bulk, "mode spectrum {sum} vs kappa_rta {bulk}");

    let ladder = ladder(Some(closure_rung(0.3, 30.0, 1e-6)), 0);
    let tol = 1e-4 * bulk;
    let mut previous = 0.0;
    for &l in &[1e-7, 1e-6, 1e-5, 1.0] {
        let goal = film_kappa(&ladder, &modes, l, Request::new(tol));
        assert!(goal.claim.error <= tol);
        let casimir = effective_kappa(&branches, &silicon_params(f64::INFINITY), t, l).kappa_effective;
        print!(
            "L = {l:.0e} m: κ_film = {:.3} ± {:.1e} W/(m K) (bulk {bulk:.1}; Casimir τ_B = L/v gives {casimir:.2}); rungs {:?}; cost {:.2e}",
            goal.kappa, goal.claim.error, goal.histogram, goal.cost
        );
        // Physics checks: monotone in L, never above bulk, and the deficit
        // bounded per mode by 1 − S(Kn) ≤ min(1, 1.5 Kn) (S ≥ the Milne
        // jump law up to its small claim; 2ζ = 1.42).
        assert!(goal.kappa > previous && goal.kappa <= bulk + goal.claim.error);
        previous = goal.kappa;
        let deficit_bound: f64 = modes.iter().map(|m| m.kappa * (1.5 * m.mfp / l).min(1.0)).sum();
        assert!(bulk - goal.kappa <= deficit_bound + goal.claim.error);
        if l == 1e-6 {
            // The same bound with a uniform per-mode tolerance tol / Σκ_i.
            let uniform: Vec<_> = modes.iter().map(|m| ladder.solve(&(m.mfp / l), Request::new(tol / bulk))).collect();
            let k_uniform: f64 = modes.iter().zip(&uniform).map(|(m, s)| m.kappa * s.value.ratio).sum();
            let claim_uniform: f64 = modes.iter().zip(&uniform).map(|(m, s)| m.kappa * s.claim.error).sum();
            let cost_uniform: f64 = uniform.iter().map(|s| s.cost).sum();
            assert!((goal.kappa - k_uniform).abs() <= goal.claim.error + claim_uniform);
            assert!(cost_uniform > 100.0 * goal.cost);
            print!("; uniform per-mode tolerance costs {cost_uniform:.2e} ({:.1e}×)", cost_uniform / goal.cost);
        }
        println!();
    }
}
