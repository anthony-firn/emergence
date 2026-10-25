//! Fresnel integrals against independent high-order quadrature and the
//! tabulated values of Abramowitz & Stegun (Table 7.7).

use num_complex::Complex64 as Cx;
use waveoptics::fresnel::fresnel;

/// Composite Gauss–Legendre (10 points) on [0, x] with panels short enough
/// to resolve the oscillation (an independent reference).
fn quad(x: f64) -> Cx {
    let (xg, wg) = (
        [-0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.1488743389816312, 0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717],
        [0.0666713443086881, 0.1494513491505806, 0.2190863625159820, 0.2692667193099963, 0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2190863625159820, 0.1494513491505806, 0.0666713443086881],
    );
    let panels = ((x * x * 20.0).ceil() as usize).max(20);
    let h = x / panels as f64;
    let mut s = Cx::new(0.0, 0.0);
    for p in 0..panels {
        let (a, b) = (p as f64 * h, (p + 1) as f64 * h);
        for (t, w) in xg.iter().zip(&wg) {
            let u = 0.5 * (a + b) + 0.5 * (b - a) * t;
            let ph = std::f64::consts::PI / 2.0 * u * u;
            s += Cx::new(ph.cos(), ph.sin()) * (0.5 * (b - a) * w);
        }
    }
    s
}

#[test]
fn matches_tables_and_independent_quadrature() {
    // A&S Table 7.7.
    for (x, c, s) in [(0.5, 0.4923442, 0.0647324), (1.0, 0.7798934, 0.4382591), (2.0, 0.4882534, 0.3434157), (3.0, 0.6057208, 0.4963130)] {
        let f = fresnel(x);
        assert!((f.re - c).abs() < 1e-7 && (f.im - s).abs() < 1e-7, "x {x}: {f}");
    }
    let mut worst: f64 = 0.0;
    for k in 0..400 {
        let x = 0.01 + 12.0 * k as f64 / 400.0;
        let e = (fresnel(x) - quad(x)).norm();
        worst = worst.max(e);
    }
    println!("worst |F − quadrature| on [0, 12]: {worst:.2e}");
    assert!(worst < 1e-12);
    // Continuity across the series / continued-fraction switch, oddness,
    // and the limit F(∞) = (1+i)/2.
    assert!((fresnel(2.5 - 1e-12) - fresnel(2.5 + 1e-12)).norm() < 1e-11);
    assert!((fresnel(-1.3) + fresnel(1.3)).norm() < 1e-16);
    assert!((fresnel(1e4) - Cx::new(0.5, 0.5)).norm() < 1e-4);
}
