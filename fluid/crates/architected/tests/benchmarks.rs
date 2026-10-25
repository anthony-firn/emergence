//! Benchmarks against closed-form and published results.
//!
//! Every number here is compared with an exact or published expression and
//! the relative error is asserted. Tolerances are physical, not tuned: a
//! convergence study tightens its own tolerance by refining `t/l`.

use architected::buckling::*;
use architected::cell::*;
use architected::lattices::*;
use architected::{Cell2d, Cell3d, Kind};

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} vs {b:e}, tol {tol:e}");
}

/// Gibson & Ashby regular-hexagon bending-only modulus.
fn gibson_ashby(t_over_l: f64) -> (f64, f64) {
    let e = 4.0 / 3.0_f64.sqrt() * t_over_l.powi(3);
    (e, 1.0)
}

/// Masters & Evans / beam-theory axial correction for a regular hexagon:
/// `E*/Es = 2.309 (t/l)^3 / (1 + 3 (t/l)^2)`.
fn regular_honeycomb_with_stretching(t_over_l: f64) -> f64 {
    4.0 / 3.0_f64.sqrt() * t_over_l.powi(3) / (1.0 + 3.0 * t_over_l.powi(2))
}

#[test]
fn honeycomb_bending_limit_matches_gibson_ashby() {
    // As t/l -> 0 the frame lattice must approach the bending-only result.
    let l = 1.0;
    let mut prev_err = f64::INFINITY;
    for &t_over_l in &[0.04, 0.02, 0.01, 0.005] {
        let t = t_over_l * l;
        let (raw, a1, a2) = regular_honeycomb(t, l, Kind::Frame);
        let cell = Cell2d::new(a1, a2, &raw, Kind::Frame, None);
        let c = cell.stiffness();
        let m = moduli2d(&c);
        let (e_ga, nu_ga) = gibson_ashby(t_over_l);
        let e_rel = (m.ex - e_ga).abs() / e_ga;
        let nu_err = (m.nu_xy - nu_ga).abs();
        // A single unit cell with rigid joints differs from the ideal
        // beam-theory cell only at O((t/l)^2); the error must fall with t/l.
        assert!(
            e_rel < 0.12,
            "t/l={t_over_l}: E* {:.4e} vs GA {:.4e}, rel {e_rel:.3}",
            m.ex,
            e_ga
        );
        assert!(nu_err < 0.06, "t/l={t_over_l}: nu {:.4} vs 1", m.nu_xy);
        assert!(e_rel < prev_err, "error did not fall with decreasing t/l");
        prev_err = e_rel;
    }
}

#[test]
fn honeycomb_includes_stretching_like_masters_evans() {
    // For the regular hexagon the stretching-inclusive beam-theory result
    // (Masters & Evans 1996) is
    //   E*/Es = (4/sqrt3) (t/l)^3 / (1 + 3 (t/l)^2),
    //   nu*   = (1 - (t/l)^2) / (1 + 3 (t/l)^2),
    // and the frame unit cell reproduces it to machine precision.
    let l = 1.0;
    for &t_over_l in &[0.05, 0.1, 0.2, 0.25] {
        let t = t_over_l * l;
        let (raw, a1, a2) = regular_honeycomb(t, l, Kind::Frame);
        let cell = Cell2d::new(a1, a2, &raw, Kind::Frame, None);
        let m = moduli2d(&cell.stiffness());
        let pred = regular_honeycomb_with_stretching(t_over_l);
        let nu_pred = (1.0 - t_over_l * t_over_l) / (1.0 + 3.0 * t_over_l * t_over_l);
        let rel = (m.ex - pred).abs() / pred;
        assert!(
            rel < 1e-9,
            "t/l={t_over_l}: E* {:.6e} vs corrected {:.6e}, rel {rel:.3e}",
            m.ex,
            pred
        );
        assert!(
            (m.nu_xy - nu_pred).abs() < 1e-9,
            "t/l={t_over_l}: nu {:.6} vs {nu_pred:.6}",
            m.nu_xy
        );
    }
}

#[test]
fn general_hexagonal_honeycomb_matches_gibson_ashby() {
    // The general Gibson & Ashby bending-only in-plane moduli
    //   E1/Es = (t/l)^3 cos(theta) / ((h/l + sin theta) sin^2 theta)
    //   nu12  = cos^2(theta) / ((h/l + sin theta) sin theta)
    // for a non-regular cell, in the slender limit.
    let l = 1.0;
    let t = 0.005 * l;
    for &(h_over_l, theta_deg) in &[(1.5, 20.0), (1.0, 45.0), (2.0, 15.0)] {
        let theta: f64 = theta_deg * std::f64::consts::PI / 180.0;
        let (raw, a1, a2, _h) = honeycomb(t, l, h_over_l, theta, Kind::Frame);
        let cell = Cell2d::new(a1, a2, &raw, Kind::Frame, None);
        let m = moduli2d(&cell.stiffness());
        let e_pred = (t / l).powi(3) * theta.cos()
            / ((h_over_l + theta.sin()) * theta.sin().powi(2));
        let nu_pred = theta.cos().powi(2) / ((h_over_l + theta.sin()) * theta.sin());
        assert!(
            (m.ex - e_pred).abs() / e_pred < 1e-3,
            "h/l={h_over_l}, theta={theta_deg}: E1 {:.4e} vs {:.4e}",
            m.ex,
            e_pred
        );
        assert!(
            (m.nu_xy - nu_pred).abs() < 1e-3,
            "h/l={h_over_l}, theta={theta_deg}: nu {:.4} vs {:.4}",
            m.nu_xy,
            nu_pred
        );
    }
}

#[test]
fn reentrant_honeycomb_has_poisson_minus_one() {
    // Masters & Evans: nu12 = cos^2(theta) / ((h/l + sin theta) sin theta)
    // is -1 for theta = -30 deg and h/l = 2 in the bending limit.
    let theta = -std::f64::consts::PI / 6.0;
    let (h_over_l, l) = (2.0, 1.0);
    let nu_pred = theta.cos().powi(2) / ((h_over_l + theta.sin()) * theta.sin());
    close(nu_pred, -1.0, 1e-12);
    // Small t/l so bending dominates.
    let t = 0.01 * l;
    let (raw, a1, a2, _h) = honeycomb(t, l, h_over_l, theta, Kind::Frame);
    let cell = Cell2d::new(a1, a2, &raw, Kind::Frame, None);
    let m = moduli2d(&cell.stiffness());
    assert!(m.nu_xy < 0.0, "not auxetic: nu = {}", m.nu_xy);
    let rel = (m.nu_xy - nu_pred).abs() / nu_pred.abs();
    assert!(rel < 0.05, "nu {:+.4} vs {nu_pred:+.4}, rel {rel:.3}", m.nu_xy);
}

#[test]
fn triangular_lattice_is_stretch_dominated_and_isotropic() {
    // Pin-jointed triangular lattice: E*/Es = rho/3, nu = 1/3, isotropic.
    let l = 1.0;
    for &t_over_l in &[0.05, 0.02, 0.01] {
        let t = t_over_l * l;
        let (raw, a1, a2) = triangular(t, l);
        let cell = Cell2d::new(a1, a2, &raw, Kind::Truss, None);
        let rho = cell.relative_density();
        let c = cell.stiffness();
        let m = moduli2d(&c);
        let e_pred = rho / 3.0;
        let rel = (m.ex - e_pred).abs() / e_pred;
        assert!(rel < 1e-9, "t/l={t_over_l}: E* {:.6e} vs {e_pred:.6e}", m.ex);
        assert!((m.ey - m.ex).abs() / m.ex < 1e-9, "not isotropic");
        close(m.nu_xy, 1.0 / 3.0, 1e-9);
        close(m.nu_yx, 1.0 / 3.0, 1e-9);
    }
}

/// Deshpande, Fleck & Ashby (2001) octet-truss result:
/// `E*/Es = (1/9) rho` at low relative density.
#[test]
fn octet_truss_matches_deshpande() {
    let a = 1.0;
    for &t_over_l in &[0.05, 0.03, 0.02] {
        // t/l is strut side over strut length; l = a/sqrt(2).
        let l = a / 2.0f64.sqrt();
        let t = t_over_l * l;
        let (raw, a1, a2, a3) = octet_truss(t, a);
        let cell = Cell3d::new(a1, a2, a3, &raw, Kind::Truss);
        let rho = cell.relative_density();
        let m = moduli3d(&cell.stiffness());
        let e_pred = rho / 9.0;
        let rel = (m.e - e_pred).abs() / e_pred;
        assert!(
            rel < 1e-9,
            "t/l={t_over_l}: E*/Es {:.6e} vs rho/9 {:.6e}, rel {rel:.4}",
            m.e,
            e_pred
        );
        // The octet truss (FCC nearest-neighbour network) has cubic elastic
        // symmetry with C11 = 2 C12 = 2 C44, i.e. Zener factor 2, and
        // nu = 1/3 along <100>. (It is not elastically isotropic.)
        let c = cell.stiffness();
        assert!((c[0][0] - 2.0 * c[0][1]).abs() / c[0][0] < 1e-9);
        assert!((c[0][1] - c[3][3]).abs() / c[0][1] < 1e-9);
        assert!((m.nu - 1.0 / 3.0).abs() < 1e-9, "nu = {} at t/l {t_over_l}", m.nu);
    }
}

#[test]
fn cubic_frame_shear_is_bending_dominated() {
    // A simple-cubic frame carries direct axial load but its shear comes from
    // strut bending, so G*/Es ~ rho^2. Check the log-log slope.
    let a = 1.0;
    let g_of = |t_over_a: f64| {
        let t = t_over_a * a;
        let (raw, a1, a2, a3) = cubic_frame(t, a);
        let cell = Cell3d::new(a1, a2, a3, &raw, Kind::Frame);
        (cell.relative_density(), moduli3d(&cell.stiffness()).g)
    };
    let (rho1, g1) = g_of(0.03);
    let (rho2, g2) = g_of(0.06);
    let slope = (g2.ln() - g1.ln()) / (rho2.ln() - rho1.ln());
    assert!(
        (slope - 2.0).abs() < 0.1,
        "shear exponent {slope:.3}, expected 2"
    );
    // The simple-cubic frame has the exact coefficient G*/Es = rho^2 / 18.
    for (rho, g) in [(rho1, g1), (rho2, g2)] {
        assert!(
            (g - rho * rho / 18.0).abs() / (rho * rho / 18.0) < 1e-9,
            "G {g:.4e} vs rho^2/18 {:.4e}",
            rho * rho / 18.0
        );
    }
}

#[test]
fn buckling_strain_scales_as_density_squared() {
    // Honeycomb Euler buckling strain eps_cr ~ (t/l)^2.
    let l = 1.0;
    let cell_of = |t_over_l: f64| {
        let t = t_over_l * l;
        let (raw, a1, a2) = regular_honeycomb(t, l, Kind::Frame);
        Cell2d::new(a1, a2, &raw, Kind::Frame, None)
    };
    let b1 = uniaxial_2d(&cell_of(0.01), 1, 1.0).critical_strain;
    let b2 = uniaxial_2d(&cell_of(0.02), 1, 1.0).critical_strain;
    let slope = (b2.ln() - b1.ln()) / (0.02f64.ln() - 0.01f64.ln());
    assert!(
        (slope - 2.0).abs() < 0.15,
        "buckling exponent {slope:.3}, expected 2"
    );
    // The linear model holds below this strain and fails above it.
    assert!(b1 > 0.0 && b1.is_finite());
}
