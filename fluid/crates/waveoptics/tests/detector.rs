//! Photon-counting detector checks: quadrature, Poisson statistics, pixel
//! integration, silicon QE model, and single-photon build-up.

use std::f64::consts::PI;
use waveoptics::detector::*;

#[test]
fn gauss_legendre_is_exact_for_polynomials() {
    for n in 1..12 {
        let (x, w) = gauss_legendre(n);
        for p in 0..(2 * n) {
            let s: f64 = x.iter().zip(&w).map(|(x, w)| w * x.powi(p as i32)).sum();
            let exact = if p % 2 == 1 { 0.0 } else { 2.0 / (p as f64 + 1.0) };
            assert!((s - exact).abs() < 1e-13, "n {n} p {p}: {s} vs {exact}");
        }
    }
}

#[test]
fn poisson_deviates_have_the_right_mean_variance_and_zero_probability() {
    let mut rng = Rng(11);
    for &m in &[0.3, 4.0, 29.0, 31.0, 400.0, 1e5] {
        let n = 200_000;
        let xs: Vec<f64> = (0..n).map(|_| rng.poisson(m)).collect();
        let mean = xs.iter().sum::<f64>() / n as f64;
        let var = xs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
        // standard errors: mean √(m/n), variance ≈ m √(2/n)
        assert!((mean - m).abs() < 5.0 * (m / n as f64).sqrt(), "mean {mean} vs {m}");
        assert!((var - m).abs() < 5.0 * m * (2.0 / n as f64).sqrt() + 1e-9, "var {var} vs {m}");
        if m < 5.0 {
            let p0 = xs.iter().filter(|&&x| x == 0.0).count() as f64 / n as f64;
            assert!((p0 - (-m).exp()).abs() < 5.0 * ((-m).exp() / n as f64).sqrt());
        }
    }
}

#[test]
fn pixel_integration_meets_tolerance_on_fine_fringes() {
    // E(x) = 1 + cos(2πx/P) with P = 2.2 pixels: the rule must integrate it.
    let s = Sensor::cmos([64, 4], 3.45e-6, 1e-3);
    let period = 2.2 * s.pitch;
    for &tol in &[1e-3, 1e-8] {
        let n = s.points_per_axis(period, tol);
        let (pts, _) = s.quadrature(n);
        let e: Vec<f64> = pts.iter().map(|p| 1.0 + (2.0 * PI * p[0] / period).cos()).collect();
        let got = s.expected_electrons(n, &[(550e-9, e)]);
        let k = s.exposure * 550e-9 / HC * s.qe.at(550e-9);
        let mut worst: f64 = 0.0;
        for j in 0..4 {
            for i in 0..64 {
                let c = s.pixel_center(i, j);
                let h = 0.5 * s.pitch;
                let exact = k * s.pitch * (2.0 * h + period / (2.0 * PI) * ((2.0 * PI * (c[0] + h) / period).sin() - (2.0 * PI * (c[0] - h) / period).sin()));
                worst = worst.max((got[i + 64 * j] - exact).abs() / (k * s.pitch * s.pitch));
            }
        }
        println!("tol {tol:e}: {n} points per axis, worst relative error {worst:.1e}");
        assert!(worst < tol);
    }
}

#[test]
fn silicon_qe_follows_the_absorption_table() {
    // Table points are reproduced exactly by the interpolation (per metre).
    assert!((silicon_alpha(600e-9) - 4.14e5).abs() < 1e-6 * 4.14e5);
    assert!((silicon_alpha(1000e-9) - 64.0e2).abs() < 1e-6 * 64e2);
    let q = Qe::bsi_cmos();
    let manual = |l: f64| {
        let a = silicon_alpha(l);
        0.95 * ((-a * 50e-9).exp() - (-a * 3.05e-6).exp())
    };
    for l in [400e-9, 550e-9, 700e-9, 900e-9] {
        assert!((q.at(l) - manual(l)).abs() < 1e-15);
    }
    println!("BSI QE: 400 nm {:.2}, 550 nm {:.2}, 700 nm {:.2}, 900 nm {:.2}", q.at(400e-9), q.at(550e-9), q.at(700e-9), q.at(900e-9));
}

#[test]
fn single_photon_buildup_converges_to_the_probability_density() {
    // χ² of N detections against the expected distribution: E[χ²] = dof.
    let p: Vec<f64> = (0..50).map(|i| (i as f64 * 0.37).cos().powi(2) + 0.05).collect();
    let total: f64 = p.iter().sum();
    let mut rng = Rng(5);
    for &n in &[1_000usize, 100_000] {
        let ev = single_photon_events(&p, n, &mut rng);
        let mut h = vec![0.0; 50];
        for e in ev {
            h[e] += 1.0;
        }
        let chi2: f64 = h.iter().zip(&p).map(|(o, pi)| {
            let e = n as f64 * pi / total;
            (o - e).powi(2) / e
        }).sum();
        println!("N = {n}: χ² = {chi2:.1} (dof 49)");
        assert!(chi2 < 49.0 + 5.0 * (2.0f64 * 49.0).sqrt(), "χ² {chi2}");
    }
}
