//! Engine tests against exact results.
//!
//! The model problems have exactly known errors, so every claim the engine
//! makes can be checked: Taylor polynomials of `exp` with the Lagrange
//! remainder (a rigorous bound), a deliberately over-optimistic asymptotic
//! indicator (calibration), closure tables of known functions, and the
//! error-budget planner against brute force.

use emergence::*;

/// Deterministic uniform numbers in [0, 1) (SplitMix64).
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn taylor_exp(x: f64, n: usize) -> f64 {
    let (mut term, mut sum) = (1.0, 1.0);
    for k in 1..=n {
        term *= x / k as f64;
        sum += term;
    }
    sum
}

fn factorial(n: usize) -> f64 {
    (1..=n).map(|k| k as f64).product()
}

/// Lagrange remainder bound for the degree-`n` Taylor polynomial of `exp` at
/// `x`: `e^{max(0,x)} |x|^{n+1} / (n+1)!`, plus round-off of the sum.
fn lagrange_bound(x: f64, n: usize) -> f64 {
    x.max(0.0).exp() * x.abs().powi(n as i32 + 1) / factorial(n + 1) + 4.0 * (n as f64 + 1.0) * f64::EPSILON * x.abs().exp()
}

/// Taylor rungs of degree 1..=20 plus `f64::exp` as the reference.
fn exp_ladder() -> Ladder<f64, f64> {
    let mut l = Ladder::new("exp", |a: &f64, b: &f64| (a - b).abs());
    for n in 1..=20usize {
        l.push(FnRung::new(
            format!("taylor-{n}"),
            move |_| n as f64,
            move |&x: &f64| {
                if x.abs() <= 1.0 {
                    Assessment::Claim(Claim::rigorous(lagrange_bound(x, n)))
                } else {
                    Assessment::OutOfDomain
                }
            },
            move |&x| taylor_exp(x, n),
        ));
    }
    let r = l.push(FnRung::new(
        "libm-exp",
        |_| 40.0,
        |&x: &f64| Assessment::Claim(Claim::rigorous(2.0 * f64::EPSILON * x.exp())),
        |&x: &f64| x.exp(),
    ));
    l.reference(r).audit(Audit::Every(1))
}

#[test]
fn rigorous_ladder_is_cost_optimal_and_never_exceeds_its_claim() {
    let ladder = exp_ladder();
    let mut rng = Rng(1);
    for &tol in &[1e-2, 1e-6, 1e-12] {
        for _ in 0..2000 {
            let x = 2.0 * rng.next() - 1.0;
            let s = ladder.solve(&x, Request::with_rigor(tol, Rigor::Rigorous));
            assert!(s.met, "tol {tol} x {x}");
            let err = (s.value - x.exp()).abs();
            assert!(err <= s.claim.error, "x {x}: error {err:e} > claim {:e}", s.claim.error);
            assert!(s.claim.error <= tol);
            // Optimality: no cheaper rung's bound meets the tolerance.
            if s.rung > 0 {
                assert!(lagrange_bound(x, s.rung) > tol, "x {x}: degree {} would do", s.rung);
            }
        }
    }
    let st = ladder.stats();
    assert!(st.iter().all(|s| s.rigorous_defects == 0 && s.violations == 0));
    let audits: u64 = st.iter().map(|s| s.audits).sum();
    assert!(audits > 0);
    println!("{}", ladder.report());
}

#[test]
fn out_of_domain_unfolds_to_the_reference_and_unmeetable_requests_say_so() {
    let ladder = exp_ladder();
    // |x| > 1: every Taylor rung is out of domain, only exp applies.
    let s = ladder.solve(&3.0, Request::new(1e-10));
    assert_eq!(ladder.rung_name(s.rung), "libm-exp");
    assert!(s.met);
    // A tolerance below round-off cannot be met: the best answer comes back
    // flagged, never silently.
    let s = ladder.solve(&0.5, Request::new(1e-30));
    assert!(!s.met);
    assert!((s.value - 0.5f64.exp()).abs() < 1e-15);
}

#[test]
fn audits_calibrate_an_over_optimistic_indicator() {
    // A cubic Taylor rung whose indicator claims 1/10 of the leading remainder
    // term: an asymptotic estimate with the wrong constant.
    let optimistic = FnRung::new(
        "taylor-3-optimistic",
        |_| 3.0,
        |&x: &f64| Assessment::Claim(Claim::asymptotic(0.1 * x.powi(4) / 24.0)),
        |&x: &f64| taylor_exp(x, 3),
    );
    let exact = FnRung::new(
        "libm-exp",
        |_| 40.0,
        |&x: &f64| Assessment::Claim(Claim::rigorous(2.0 * f64::EPSILON * x.exp())),
        |&x: &f64| x.exp(),
    );
    let ladder = Ladder::new("calibration", |a: &f64, b: &f64| (a - b).abs())
        .with(optimistic)
        .with(exact)
        .reference(1)
        .audit(Audit::Every(1));
    let mut rng = Rng(7);
    let mut late_violations = 0;
    for i in 0..4000 {
        let x = 2.0 * rng.next() - 1.0;
        let s = ladder.solve(&x, Request::new(1e-3));
        if i >= 1000 && s.rung == 0 && (s.value - x.exp()).abs() > s.claim.error {
            late_violations += 1;
        }
    }
    let st = &ladder.stats()[0];
    println!("{}", ladder.report());
    // The true constant is e^ξ/24 with ξ ∈ (−1, 1): up to 27× the claim.
    assert!(st.worst_ratio > 8.0 && st.worst_ratio < 30.0, "worst ratio {}", st.worst_ratio);
    assert!(st.safety >= st.worst_ratio, "safety {} below worst ratio", st.safety);
    assert!(st.violations <= 6, "{} violations before calibration settled", st.violations);
    assert_eq!(late_violations, 0);
}

#[test]
fn failed_certificate_unfolds_to_the_next_rung() {
    // A rung that only knows its error after running: certified for x ≥ 0.
    let deferred = FnRung::new("cheap-deferred", |_| 1.0, |_: &f64| Assessment::Deferred, |&x: &f64| x)
        .with_certificate(|&x: &f64, _out: &f64| Some(Claim::rigorous(if x >= 0.0 { 1e-12 } else { 1.0 })));
    let exact = FnRung::new(
        "exact",
        |_| 100.0,
        |_: &f64| Assessment::Claim(Claim::rigorous(0.0)),
        |&x: &f64| x,
    );
    let ladder = Ladder::new("deferred", |a: &f64, b: &f64| (a - b).abs()).with(deferred).with(exact);
    let a = ladder.solve(&0.5, Request::new(1e-6));
    assert_eq!((a.rung, a.cost), (0, 1.0));
    let b = ladder.solve(&-0.5, Request::new(1e-6));
    assert_eq!((b.rung, b.cost), (1, 101.0));
    assert!(b.met);
}

#[test]
fn closure_claims_hold_on_held_out_random_points() {
    let f = |x: f64| [(10.0 * x).sin() * x.exp(), 1.0 / (1.0 + 25.0 * x * x)];
    for &tol in &[1e-6, 1e-10] {
        let c = Closure::build(f, -1.0, 1.0, Axis::Linear, tol, ClosureOptions::default());
        assert_eq!(c.coverage(), 1.0);
        let mut rng = Rng(3);
        let mut worst: f64 = 0.0;
        for _ in 0..20000 {
            let x = 2.0 * rng.next() - 1.0;
            let (v, claim) = c.eval(x).unwrap();
            assert!(claim <= tol);
            let e = f(x).iter().zip(&v).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
            assert!(e <= claim, "x {x}: error {e:e} > claim {claim:e}");
            worst = worst.max(e / claim);
        }
        println!(
            "tol {tol:e}: {} patches, {} samples, worst error/claim {worst:.2e}",
            c.patches(),
            c.samples
        );
    }
}

#[test]
fn log_closure_spans_decades_and_unconverged_patches_unfold() {
    // A power law over 12 decades is smooth in ln x. (Values reach 100, so an
    // absolute 1e-12 would be below round-off and correctly never converge.)
    let c = Closure::build(|x: f64| [x.powf(1.0 / 3.0)], 1e-6, 1e6, Axis::Log, 1e-10, ClosureOptions::default());
    assert_eq!(c.coverage(), 1.0);
    let mut rng = Rng(5);
    for _ in 0..5000 {
        let x = 10f64.powf(12.0 * rng.next() - 6.0);
        let (v, claim) = c.eval(x).unwrap();
        assert!((v[0] - x.cbrt()).abs() <= claim, "x {x}");
    }
    // A jump cannot be interpolated: the patch around it stays unconverged and
    // reports out-of-domain, so a ladder unfolds there.
    let opt = ClosureOptions {
        min_width: 1e-4,
        ..ClosureOptions::default()
    };
    let h = Closure::build(|x: f64| [if x < 0.3 { 0.0 } else { 1.0 }], 0.0, 1.0, Axis::Linear, 1e-8, opt);
    let bad = h.unconverged();
    assert!(!bad.is_empty() && bad.iter().all(|&(a, b)| a <= 0.3 + 1e-4 && b >= 0.3 - 1e-4));
    assert!(h.coverage() > 0.999 && h.coverage() < 1.0);
    assert!(h.eval(0.3).is_none());
    assert!(h.eval(0.1).is_some());
}

/// Brute-force minimum cost over all rung choices with Σ e ≤ tol.
fn brute_force(costs: &[Vec<f64>], errs: &[Vec<f64>], tol: f64) -> f64 {
    fn rec(i: usize, costs: &[Vec<f64>], errs: &[Vec<f64>], left: f64, acc: f64, best: &mut f64) {
        if i == costs.len() {
            if left >= -1e-15 {
                *best = best.min(acc);
            }
            return;
        }
        for r in 0..costs[i].len() {
            rec(i + 1, costs, errs, left - errs[i][r], acc + costs[i][r], best);
        }
    }
    let mut best = f64::INFINITY;
    rec(0, costs, errs, tol, 0.0, &mut best);
    best
}

#[test]
fn budget_plan_is_feasible_and_near_optimal() {
    let mut rng = Rng(11);
    for _trial in 0..200 {
        let (n_items, n_rungs) = (6, 4);
        // Item i, rung r: cost grows and error falls with r (random factors).
        let mut costs = vec![vec![0.0; n_rungs]; n_items];
        let mut errs = vec![vec![0.0; n_rungs]; n_items];
        for i in 0..n_items {
            let (mut c, mut e) = (1.0, 1.0);
            for r in 0..n_rungs {
                costs[i][r] = c;
                errs[i][r] = e;
                c *= 1.5 + 4.0 * rng.next();
                e *= 0.02 + 0.5 * rng.next();
            }
        }
        let costs_c = costs.clone();
        let errs_c = errs.clone();
        let mut ladder: Ladder<usize, f64> = Ladder::new("knapsack", |a: &f64, b: &f64| (a - b).abs());
        for r in 0..n_rungs {
            let (cc, ee) = (costs_c.clone(), errs_c.clone());
            ladder.push(FnRung::new(
                format!("r{r}"),
                move |&i: &usize| cc[i][r],
                move |&i: &usize| Assessment::Claim(Claim::rigorous(ee[i][r])),
                |_| 0.0,
            ));
        }
        let items: Vec<usize> = (0..n_items).collect();
        let weights = vec![1.0; n_items];
        let max_err: f64 = errs.iter().map(|e| e[n_rungs - 1]).sum();
        let tol = max_err * (1.0 + 30.0 * rng.next());
        let plan = ladder.plan(&items, &weights, Request::new(tol));
        assert!(plan.bound <= tol * (1.0 + 1e-12), "infeasible plan");
        let opt = brute_force(&costs, &errs, tol);
        // Greedy on convex hulls: within the largest single upgrade cost.
        let max_step = costs.iter().map(|c| c[n_rungs - 1]).fold(0.0, f64::max);
        assert!(plan.cost <= opt + max_step, "plan {} vs optimum {opt}", plan.cost);
        assert!(plan.cost >= opt * (1.0 - 1e-12));
    }
}

#[test]
fn goal_oriented_quadrature_beats_a_uniform_tolerance() {
    // J = ∫_0^1 e^x dx by composite Simpson on 2001 nodes, each node value from
    // the exp ladder. Simpson's own error is ≤ h^4 max|f''''|/180 = 1.5e-14.
    let ladder = exp_ladder().audit(Audit::Off);
    let n = 2000;
    let h = 1.0 / n as f64;
    let nodes: Vec<f64> = (0..=n).map(|i| i as f64 * h).collect();
    let weights: Vec<f64> = (0..=n)
        .map(|i| h / 3.0 * if i == 0 || i == n { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 })
        .collect();
    // Weight the right half of the interval 100× less by splitting J into
    // J = ∫_0^{1/2} + 0.01 ∫_{1/2}^1: the far half needs less accuracy.
    // The node at x = 1/2 closes a left panel (h/3) and opens a right one.
    let weights: Vec<f64> = (0..=n)
        .map(|i| match i.cmp(&(n / 2)) {
            std::cmp::Ordering::Less => weights[i],
            std::cmp::Ordering::Equal => h / 3.0 * 1.01,
            std::cmp::Ordering::Greater => 0.01 * weights[i],
        })
        .collect();
    let exact = (0.5f64.exp() - 1.0) + 0.01 * (1f64.exp() - 0.5f64.exp());
    let tol = 1e-10;
    let simpson_err = 1.5e-14;

    let f = ladder.solve_functional(&nodes, &weights, Request::with_rigor(tol, Rigor::Rigorous));
    let j = f.sum(&weights, |v| *v);
    assert!(f.met);
    assert!((j - exact).abs() <= f.claim.error + simpson_err, "error {:e} > claim {:e}", (j - exact).abs(), f.claim.error);

    // The same accuracy with a uniform per-node tolerance tol / Σ|w|.
    let w_sum: f64 = weights.iter().sum();
    let uniform: Vec<_> = nodes.iter().map(|x| ladder.solve(x, Request::with_rigor(tol / w_sum, Rigor::Rigorous))).collect();
    let uniform_cost: f64 = uniform.iter().map(|s| s.cost).sum();
    println!(
        "goal-oriented cost {:.0} vs uniform {:.0} ({:.2}×); claim {:.2e}, actual {:.2e}",
        f.cost,
        uniform_cost,
        uniform_cost / f.cost,
        f.claim.error,
        (j - exact).abs()
    );
    assert!(f.cost < uniform_cost);
}

#[test]
fn refine_functional_reuses_unchanged_items_and_matches_a_fresh_solve() {
    let ladder = exp_ladder().audit(Audit::Off);
    let xs: Vec<f64> = (0..500).map(|i| -1.0 + 2.0 * i as f64 / 499.0).collect();
    let w1 = vec![1.0 / 500.0; 500];
    let pilot = ladder.solve_functional(&xs, &w1, Request::with_rigor(1e-6, Rigor::Rigorous));
    // New weights: the right half matters 1000× more.
    let w2: Vec<f64> = xs.iter().map(|&x| if x > 0.0 { 1.0 } else { 1e-3 } / 500.0).collect();
    let fresh = ladder.solve_functional(&xs, &w2, Request::with_rigor(1e-9, Rigor::Rigorous));
    let refined = ladder.refine_functional(pilot.clone(), &xs, &w2, Request::with_rigor(1e-9, Rigor::Rigorous));
    assert!(refined.met && fresh.met);
    let reused = refined.items.iter().zip(&pilot.items).filter(|(r, p)| r.rung == p.rung).count();
    for (r, f) in refined.items.iter().zip(&fresh.items) {
        assert_eq!(r.rung, f.rung);
        assert_eq!(r.value, f.value);
    }
    assert!(reused > 0 && refined.cost < fresh.cost);
    println!("reused {reused}/500 items; refine cost {:.0} vs fresh {:.0}", refined.cost, fresh.cost);
}

#[test]
fn closure_lookup_is_exact_at_patch_and_bucket_boundaries() {
    // Uneven patches (a steep feature near 0.7) and probes exactly on every
    // patch boundary and on the bucket grid: the value must match the
    // function to the claim at every probe.
    let f = |x: f64| [(-(x - 0.7f64).powi(2) * 400.0).exp() + x];
    let c = Closure::build(f, 0.0, 1.0, Axis::Linear, 1e-10, ClosureOptions::default());
    let mut probes: Vec<f64> = (0..=4096).map(|i| i as f64 / 4096.0).collect();
    let mut rng = Rng(17);
    probes.extend((0..4096).map(|_| rng.next()));
    for x in probes {
        let (v, claim) = c.eval(x).unwrap();
        assert!((v[0] - f(x)[0]).abs() <= claim, "x {x}: {:e} > {claim:e}", (v[0] - f(x)[0]).abs());
    }
}
