//! The Mie ladder on the emergence engine: every claim checked against the
//! exact series, the resonance gate on closures, and the size-distribution
//! average with a goal-oriented budget.
//!
//! Run with `cargo test -p photonics --test mie_ladder -- --nocapture` for the
//! tables quoted in `analysis/emergence_engine.md`.

use std::time::Instant;

use emergence::Request;
use photonics::mie::GammaSizeDistribution;
use photonics::mie_ladder::*;
use photonics::C;

/// Bond & Bergstrom (2006) recommended black-carbon index at 550 nm.
const SOOT: C = C::new(1.95, 0.79);
/// Mineral dust at 550 nm (estimate; typical of the d'Almeida/OPAC range).
const DUST: C = C::new(1.53, 0.008);
/// Liquid water near 2.13 µm (estimate from Hale & Querry 1973).
const WATER_2_13UM: C = C::new(1.29, 5.0e-4);
/// Liquid water at 550 nm (Hale & Querry 1973: κ ≈ 2e-9).
const WATER_VISIBLE: C = C::new(1.333, 2.0e-9);

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

/// The exact rung's own claimed truncation error (x ≤ 2000).
const EXACT_ERR: f64 = 2e-8;

#[test]
fn rayleigh2_returns_less_than_it_claims() {
    // R2 against the exact series where Rayleigh is in domain.
    for &m in &[SOOT, DUST, WATER_2_13UM, C::new(3.0, 1.0)] {
        for k in 0..60 {
            let x = 10f64.powf(-3.0 + 3.0 * k as f64 / 60.0);
            if x.max(m.abs() * x) > 1.0 {
                continue;
            }
            let claim = metric(&rayleigh2(x, m), &rayleigh1(x, m));
            let err = metric(&rayleigh2(x, m), &exact(x, m));
            assert!(err <= claim + EXACT_ERR, "m {m:?} x {x}: R2 error {err:e} > claim {claim:e}");
        }
    }
}

#[test]
fn per_point_claims_hold_and_the_ladder_picks_cheap_rungs() {
    let materials = [SOOT, DUST, WATER_2_13UM];
    let t = Instant::now();
    let closures: Vec<_> = materials.iter().map(|&m| closure_rung(m, 0.3, 200.0, 1e-7, 1 << 15)).collect();
    let build = t.elapsed();
    for c in &closures {
        let cl = c.closure();
        println!(
            "{:>34}: {:5} patches, {:7} exact samples, coverage {:.4}",
            emergence::Rung::<MieInput, MieQ>::name(c),
            cl.patches(),
            cl.samples,
            cl.coverage()
        );
        assert!((cl.coverage() - 1.0).abs() < 1e-12);
    }
    println!("closure build: {build:.2?}");
    let ladder = ladder(closures, 1);
    let mut rng = Rng(42);
    for &tol in &[1e-3, 1e-5, 1e-7] {
        ladder.reset_counters();
        let mut cost = 0.0;
        let mut cost_exact = 0.0;
        let mut worst: f64 = 0.0;
        for _ in 0..3000 {
            let m = materials[(rng.next() * 3.0) as usize];
            let x = 10f64.powf(-3.0 + (200f64.log10() + 3.0) * rng.next());
            let s = ladder.solve(&MieInput { x, m }, Request::new(tol));
            let reference = ladder.len() - 1;
            let err = if s.rung == reference { 0.0 } else { s.audited.expect("audit every result") };
            assert!(s.met, "tol {tol} unmet at x {x}");
            assert!(err <= s.claim.error + EXACT_ERR, "m {m:?} x {x}: error {err:e} > claim {:e}", s.claim.error);
            worst = worst.max(err / s.claim.error.max(EXACT_ERR));
            cost += s.cost;
            cost_exact += exact_cost(x);
        }
        println!("tol {tol:e}: cost {:.1}× below exact everywhere, worst error/claim {worst:.2}", cost_exact / cost);
        println!("{}", ladder.report());
        let st = ladder.stats();
        assert!(st.iter().all(|s| s.violations == 0));
    }
    // Below the reference's own validated accuracy the request is flagged
    // unmet (with the best answer), never reported as met.
    let s = ladder.solve(&MieInput { x: 50.0, m: SOOT }, Request::new(1e-10));
    assert!(!s.met);
}

#[test]
fn resonance_finder_is_step_independent_and_widths_are_half_maxima() {
    let m = WATER_VISIBLE;
    let r1 = resonances(m, 60.0, 64.0, 0.01, 1e-9);
    let r2 = resonances(m, 60.0, 64.0, 0.005, 1e-9);
    assert_eq!(r1.len(), r2.len(), "the scan step changed the resonance count");
    for (a, b) in r1.iter().zip(&r2) {
        assert!((a.x - b.x).abs() <= 1e-9 * a.x && a.order == b.order);
    }
    let narrow = r1.iter().map(|r| r.width).fold(f64::INFINITY, f64::min);
    let broad = r1.iter().map(|r| r.width).fold(0.0, f64::max);
    println!("{} resonances in 60 ≤ x ≤ 64, widths {narrow:.2e} … {broad:.2e}", r1.len());
    // |c_n|² is a Lorentzian: half its peak value at x* ± Γ/2.
    for r in r1.iter().filter(|r| r.width < 1e-2) {
        let c2 = |x: f64| {
            let s = photonics::mie::MieScattering::with_terms(x, m, r.order);
            let c = if r.electric { s.a[r.order - 1] } else { s.b[r.order - 1] };
            c.abs2()
        };
        let half = 0.5 * (c2(r.x - 0.5 * r.width) + c2(r.x + 0.5 * r.width)) / c2(r.x);
        assert!((half - 0.5).abs() < 0.02, "order {} at x {}: half-max ratio {half}", r.order, r.x);
    }
}

#[test]
fn resonances_narrower_than_a_table_keep_it_out_of_domain_exactly_there() {
    // Visible water near x = 60: radiative Q of the highest orders reaches the
    // absorption limit n/(2κ) ~ 3e8, so some resonances are far narrower than
    // any table spacing.
    let m = WATER_VISIBLE;
    let (lo, hi) = (60.0, 64.0);
    let peaks = resonances(m, lo, hi, 0.01, 1e-9);
    // A table built without the resonance list passes its held-out check
    // between resonances and is wrong on them.
    let blind = emergence::Closure::build(|x| exact(x, m).to_vec(), lo, hi, emergence::Axis::Linear, 1e-7, emergence::ClosureOptions::default());
    let mut worst_blind: f64 = 0.0;
    for r in &peaks {
        let (v, claim) = blind.eval(r.x).expect("blind table claims everything");
        worst_blind = worst_blind.max(metric(&MieQ::from_slice(v), &exact(r.x, m)) / claim);
    }
    // The resonance-gated closure resolves every resonance it can afford and
    // leaves the rest out of domain. With a generous patch budget it resolves
    // all; with a small one it must exclude the narrowest.
    for budget in [1 << 14, 48] {
        let gated = closure_rung(m, lo, hi, 1e-7, budget);
        let (mut resolved, mut excluded) = (0, 0);
        for r in &peaks {
            match gated.closure().eval(r.x) {
                Some((v, claim)) => {
                    let err = metric(&MieQ::from_slice(v), &exact(r.x, m));
                    assert!(err <= claim + EXACT_ERR, "x {}: error {err:e} > claim {claim:e}", r.x);
                    resolved += 1;
                }
                None => excluded += 1,
            }
        }
        let cov = gated.closure().coverage();
        println!(
            "{} resonances; blind table worst error/claim {worst_blind:.2e}; gated table ({budget} patches max) \
             resolves {resolved}, excludes {excluded}, covers {cov:.6} with {} patches",
            peaks.len(),
            gated.closure().patches()
        );
        if budget == 48 {
            assert!(excluded > 0 && cov < 1.0);
        } else {
            assert_eq!(excluded, 0);
        }
        // Through the ladder every resonance comes back within its claim.
        let ladder = ladder(vec![gated], 1);
        for r in &peaks {
            let s = ladder.solve(&MieInput { x: r.x, m }, Request::new(1e-6));
            assert!(s.met && s.audited.unwrap_or(0.0) <= s.claim.error + EXACT_ERR);
        }
    }
    assert!(worst_blind > 10.0, "expected the blind table to miss a resonance");
}

#[test]
fn resonance_aware_closure_holds_its_claim_on_the_resonances() {
    // Water at 2.13 µm (κ = 5e-4): check the claim exactly where it is
    // hardest, on every resonance peak.
    let c = closure_rung(WATER_2_13UM, 0.3, 40.0, 1e-7, 1 << 15);
    let peaks = resonances(WATER_2_13UM, 5.0, 40.0, 0.01, 1e-9);
    assert!(peaks.len() > 50, "found only {} resonances", peaks.len());
    let mut worst: f64 = 0.0;
    for r in &peaks {
        let (v, claim) = c.closure().eval(r.x).expect("covered");
        let err = metric(&MieQ::from_slice(v), &exact(r.x, WATER_2_13UM));
        assert!(err <= claim + EXACT_ERR, "x {}: error {err:e} > claim {claim:e}", r.x);
        worst = worst.max(err / claim);
    }
    println!("{} resonances in 5 ≤ x ≤ 40 at κ = 5e-4: worst error/claim {worst:.2e}", peaks.len());
}

#[test]
fn cloud_field_averages_match_exact_mie_at_every_node() {
    // 2.13 µm (a MODIS droplet-size band), gamma distributions with
    // r_eff = 4 … 20 µm and v_eff = 0.1 (Hansen & Travis 1974), one per cell.
    let wavelength = 2.13e-6;
    let m = WATER_2_13UM;
    let cells = 200;
    let n_quad = 2000;
    let x_max = 2.0 * std::f64::consts::PI / wavelength * 20e-6 * (1.0 + 8.0 * 0.1f64.sqrt()) * 1.01;
    let t = Instant::now();
    let ladder = ladder(vec![closure_rung(m, 0.3, x_max, 1e-7, 1 << 15)], 0);
    let build = t.elapsed();

    let dists: Vec<GammaSizeDistribution> = (0..cells)
        .map(|c| GammaSizeDistribution {
            r_eff: 4e-6 + 16e-6 * c as f64 / (cells - 1) as f64,
            v_eff: 0.1,
        })
        .collect();
    let req = Request::new(1e-5);
    let t = Instant::now();
    // Both paths run in parallel over cells (and nested over nodes).
    use rayon::prelude::*;
    let engine: Vec<AveragedQ> = dists
        .par_iter()
        .map(|d| gamma_averaged(&ladder, d, m, wavelength, 1.0, n_quad, req).0)
        .collect();
    let t_engine = t.elapsed();

    // Direct: the exact series at every node, one series per node (the fair
    // baseline; photonics::mie::gamma_averaged evaluates it four times).
    let t = Instant::now();
    let direct: Vec<MieQ> = dists
        .par_iter()
        .map(|d| {
            let r_max = d.r_eff * (1.0 + 8.0 * d.v_eff.sqrt());
            let h = r_max / n_quad as f64;
            let k = 2.0 * std::f64::consts::PI / wavelength;
            let (s, a) = (0..=n_quad)
                .into_par_iter()
                .map(|i| {
                    let r = i as f64 * h;
                    let w = h / 3.0 * if i == 0 || i == n_quad { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 };
                    let wa = w * std::f64::consts::PI * r * r * d.pdf(r);
                    let q = exact(k * r, m);
                    ([wa * q.q_ext, wa * q.q_sca, wa * q.gq_sca], wa)
                })
                .reduce(|| ([0.0; 3], 0.0), |(p, a), (q, b)| ([p[0] + q[0], p[1] + q[1], p[2] + q[2]], a + b));
            MieQ { q_ext: s[0] / a, q_sca: s[1] / a, gq_sca: s[2] / a }
        })
        .collect();
    let t_direct = t.elapsed();

    let mut worst: f64 = 0.0;
    let mut cost = 0.0;
    // Predicted cost of the exact series at every node of every cell.
    let cost_exact: f64 = dists
        .iter()
        .map(|d| {
            let r_max = d.r_eff * (1.0 + 8.0 * d.v_eff.sqrt());
            let k = 2.0 * std::f64::consts::PI / wavelength;
            (0..=n_quad).map(|i| exact_cost(k * r_max * i as f64 / n_quad as f64)).sum::<f64>()
        })
        .sum();
    for (e, d) in engine.iter().zip(&direct) {
        let err = metric(&e.q, d);
        assert!(e.claim.error <= req.tol);
        assert!(err <= e.claim.error + EXACT_ERR, "error {err:e} > claim {:e}", e.claim.error);
        worst = worst.max(err);
        cost += e.cost;
    }
    // Consistency with the existing (four-series-per-node) implementation.
    let legacy = photonics::mie::gamma_averaged(&dists[100], m, wavelength, 1.0, n_quad);
    assert!((legacy.q_ext - engine[100].q.q_ext).abs() <= 1e-5 + EXACT_ERR);
    println!(
        "{cells} cells × {} nodes: worst |error| {worst:.2e} (claim ≤ {:.0e}); predicted cost {:.1}× below exact; \
         wall time engine {t_engine:.2?} (+ closure build {build:.2?}) vs direct {t_direct:.2?}",
        n_quad + 1,
        req.tol,
        cost_exact / cost
    );
    println!("{}", ladder.report());
}
