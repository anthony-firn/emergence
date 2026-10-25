//! Spreading resistance against its exact limits.

use project::thermal::spreading;

#[test]
fn whole_plate_source_has_no_spreading() {
    let (r, _) = spreading(0.04, 0.05, 0.04, 0.05, 0.003, 200.0, 500.0);
    assert!(r.abs() < 1e-15, "{r}");
}

#[test]
fn small_square_source_approaches_the_half_space() {
    // A uniform-flux square of side c on a half-space has mean rise
    // R k c = I / (2π), with I = (4/3)(1 − √2 + 3 asinh 1) the self-potential
    // of the unit square: 0.473201 (Yovanovich quotes 0.4732).
    let exact = (4.0 / 3.0) * (1.0 - 2f64.sqrt() + 3.0 * 1f64.asinh()) / (2.0 * std::f64::consts::PI);
    assert!((exact - 0.4732).abs() < 5e-5);
    let c = 1.0;
    let k = 1.0;
    // A plate L × L × L/2 with a weak sink below: the finite plate departs
    // from the half-space at O(c/L), so extrapolate in c/L.
    let r = |l: f64| spreading(l, l, c, c, 0.5 * l, k, 1e-3).0 * k * c;
    let (r1, r2) = (r(20.0), r(40.0));
    let r0 = 2.0 * r2 - r1;
    println!("R k c: L=20c {r1:.6}, L=40c {r2:.6}, extrapolated {r0:.6}, exact {exact:.6}");
    assert!((r0 / exact - 1.0).abs() < 2e-3, "{r0} vs {exact}");
    assert!(r1 < r2 && r2 < exact);
}

#[test]
fn thin_plate_with_strong_cooling_barely_spreads() {
    // Bi = h t / k ≫ 1 under a thin plate: heat goes straight down, so the
    // source sees the 1-D resistance of its own footprint. Spreading is then
    // 1/(h c d) + t/(k c d) minus the plate's 1-D value, to O(t/c).
    let (a, b, c, d, t, k, h) = (0.04, 0.05, 0.01, 0.015, 1e-4, 200.0, 1e7);
    let (r, err) = spreading(a, b, c, d, t, k, h);
    let footprint = (t / k + 1.0 / h) / (c * d) - (t / k + 1.0 / h) / (a * b);
    println!("thin plate: {r:.5e} vs footprint {footprint:.5e} (truncation {err:.1e})");
    assert!((r / footprint - 1.0).abs() < 0.05, "{r} vs {footprint}");
}
