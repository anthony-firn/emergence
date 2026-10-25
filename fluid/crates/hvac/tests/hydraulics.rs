//! Friction, drops and operating points.

use hvac::hydraulics::*;

#[test]
fn friction_factor_against_published_approximations() {
    // Colebrook solved exactly (residual checked); Haaland (1983) is an
    // explicit fit published within about 1.5 % of it, which it meets. The
    // Swamee–Jain fit (1976), often quoted at 1 %, is measured here: 2.83 %
    // at its domain's corner (Re 5000, ε/D 0.01), better than 1 % only away
    // from it. It is recorded as a fact about that fit.
    let (mut wh, mut wsj): (f64, f64) = (0.0, 0.0);
    for re in [4e3, 5e3, 1e4, 1e5, 1e6, 1e7, 1e8] {
        for rr in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2] {
            let (f, regime) = friction(re, rr);
            assert_eq!(regime, Regime::Turbulent);
            // The Colebrook residual itself.
            let r = 1.0 / f.sqrt() + 2.0 * (rr / 3.7 + 2.51 / (re * f.sqrt())).log10();
            assert!(r.abs() < 1e-12, "residual {r}");
            wh = wh.max((haaland(re, rr) / f - 1.0).abs());
            // Swamee–Jain's stated domain starts at Re = 5000.
            if re >= 5e3 {
                wsj = wsj.max((swamee_jain(re, rr) / f - 1.0).abs());
            }
        }
    }
    println!("worst against Colebrook: Haaland {:.2} %, Swamee–Jain {:.2} %", 100.0 * wh, 100.0 * wsj);
    assert!(wh < 0.015, "Haaland {wh}");
    assert!(wsj < 0.03, "Swamee–Jain {wsj}");
    // Laminar is exact; the transitional band is flagged.
    assert_eq!(friction(1000.0, 1e-3), (0.064, Regime::Laminar));
    assert_eq!(friction(3000.0, 1e-3).1, Regime::Transitional);
}

#[test]
fn drops_and_hagen_poiseuille() {
    // Laminar oil in a pipe: Δp = 128 μ L Q / (π D⁴) exactly.
    let run = Run::round(0.02, 10.0, 0.0, 0.0);
    let (rho, mu, q) = (900.0, 0.1, 1e-4);
    let d = run.drop(q, rho, mu);
    let hp = 128.0 * mu * 10.0 * q / (std::f64::consts::PI * 0.02f64.powi(4));
    assert_eq!(d.regime, Regime::Laminar);
    assert!((d.dp / hp - 1.0).abs() < 1e-12);
    // A square duct's equivalent diameter is 1.0932 a (ASHRAE Table 3 lists
    // 109 mm for 100 × 100 mm).
    assert!((equivalent_diameter(0.1, 0.1) - 0.10932).abs() < 1e-4);
}

#[test]
fn operating_point_against_the_exact_intersection() {
    // A fully rough system is a pure quadratic Δp = k Q²; with a quadratic
    // fan curve the intersection is a root of a quadratic.
    let fan = Curve::through([(0.0, 600.0, 0.3), (0.5, 520.0, 0.65), (1.0, 300.0, 0.6)]);
    let k = 400.0;
    let op = operating_point(&fan, &|q| k * q * q, 2.0).unwrap();
    let (a, b, c) = (fan.c[2] - k, fan.c[1], fan.c[0]);
    let exact = (-b - (b * b - 4.0 * a * c).sqrt()) / (2.0 * a);
    println!("fan meets system at {:.4} m³/s, {:.1} Pa, η {:.3}, {:.1} W", op.flow, op.pressure, op.efficiency, op.power);
    assert!((op.flow / exact - 1.0).abs() < 1e-12);
    assert!((fan.pressure(0.5) - 520.0).abs() < 1e-9 && (fan.efficiency(0.5) - 0.65).abs() < 1e-12);
    // A system the fan cannot drive is reported, not extrapolated.
    assert!(operating_point(&fan, &|q| 700.0 + q, 2.0).is_none());
}
