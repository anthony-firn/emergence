//! Benchmarks against exact and published results.

use wetting::capillary::*;
use wetting::constants::*;
use wetting::lattice::{self, Branch, SquareWave};
use wetting::texture::*;
use wetting::young_laplace::*;

fn close(a: f64, b: f64, tol: f64, msg: &str) {
    assert!(
        (a - b).abs() <= tol,
        "{msg}: {a:e} vs {b:e}, |Δ|={:e} > {tol:e}",
        (a - b).abs()
    );
}

// ---------------------------------------------------------------------------
// 1. Exact capillary statics
// ---------------------------------------------------------------------------

#[test]
fn zero_gravity_sessile_is_spherical_cap() {
    // With Bo = 0 the Bashforth–Adams integration must reproduce the exact
    // spherical cap h = R(1-cosθ), a = R sinθ.
    for &deg in &[30.0_f64, 60.0, 90.0, 120.0] {
        let th = deg.to_radians();
        let d = AxisymmetricDrop::sessile_to_contact_angle(0.0, th, 1e-4);
        let e = d.last();
        close(e.z, 1.0 - th.cos(), 3e-4, "cap height");
        close(e.r, th.sin(), 3e-4, "cap base radius");
        let c = th.cos();
        let v_exact = std::f64::consts::PI / 3.0 * (2.0 - 3.0 * c + c.powi(3));
        close(e.volume, v_exact, 3e-3, "cap volume");
    }
}

#[test]
fn spherical_cap_geometry_relations() {
    for &deg in &[20.0_f64, 55.0, 110.0, 160.0] {
        let th = deg.to_radians();
        let cap = SphericalCap::new(th, 1.7);
        close(cap.base_radius(), cap.radius * th.sin(), 1e-14, "base radius");
        close(cap.height(), cap.radius * (1.0 - th.cos()), 1e-14, "height");
        let from_v = SphericalCap::from_volume(th, cap.volume());
        close(from_v.radius, cap.radius, 1e-12, "volume round trip");
        close(cap.height() / cap.base_radius(), (th / 2.0).tan(), 1e-12, "h/a");
    }
}

#[test]
fn pendant_matches_opendrop_reference() {
    // Reference profiles from OpenDrop (Berry, Neeson, Dagastine, Chan &
    // Tabor, J. Colloid Interface Sci. 454, 226 (2015)), test_younglaplace.cpp,
    // at Bond number Bo = ρ g b²/γ = 0.21.
    let d = AxisymmetricDrop::pendant(0.21, 4.0, 1e-4);
    let refs = [
        (0.4, 0.389470759, 0.0787805625),
        (0.8, 0.718911827, 0.301185635),
        (1.6, 1.03671685, 1.01682486),
        (3.2, 0.353242294, 2.40649189),
    ];
    for &(s, r, z) in &refs {
        let e = d.at_s(s);
        close(e.r, r, 1e-4, "pendant r");
        close(e.z, z, 1e-4, "pendant z");
    }
    let e = d.at_s(4.0);
    close(e.volume, 5.53648, 2e-3, "pendant volume");
    close(e.area, 15.9890, 3e-3, "pendant area");
}

#[test]
fn jurin_rise_water() {
    let h = jurin_rise(WATER_GAMMA, 0.0, WATER_RHO, G_STANDARD, 0.5e-3);
    let expect = 2.0 * WATER_GAMMA / (WATER_RHO * G_STANDARD * 0.5e-3);
    close(h, expect, 1e-12, "jurin formula");
    let lc = water_capillary_length();
    close(h, 2.0 * lc * lc / 0.5e-3, 1e-12, "jurin via l_c");
    assert!((h - 0.02975).abs() < 2e-4, "jurin ~29.75 mm, got {h}");
}

#[test]
fn planar_puddle_and_wall_meniscus_exact() {
    let lc = water_capillary_length();
    for &deg in &[30.0_f64, 60.0, 90.0, 120.0] {
        let th = deg.to_radians();
        let hp = planar_puddle_height(WATER_GAMMA, th, WATER_RHO, G_STANDARD);
        close(hp, 2.0 * lc * (th / 2.0).sin(), 1e-12, "puddle height");
        let prof = planar_puddle_profile(WATER_GAMMA, th, WATER_RHO, G_STANDARD, 1e-3);
        let znum = prof.last().unwrap().1;
        close(znum, hp, 1e-5, "puddle RK4");
        let zw = planar_meniscus_wall_height(WATER_GAMMA, th, WATER_RHO, G_STANDARD);
        close(zw, lc * (2.0 * (1.0 - th.sin())).sqrt(), 1e-12, "wall height");
    }
}

#[test]
fn planar_meniscus_profile_matches_ode() {
    let (gamma, rho, g) = (WATER_GAMMA, WATER_RHO, G_STANDARD);
    for &deg in &[30.0_f64, 60.0, 90.0] {
        let th = deg.to_radians();
        let prof = planar_meniscus_profile(gamma, th, rho, g, 40.0, 1e-3);
        let target_z = 0.5 * planar_meniscus_wall_height(gamma, th, rho, g);
        let mut best = prof[0];
        for p in &prof {
            if (p.z - target_z).abs() < (best.z - target_z).abs() {
                best = *p;
            }
        }
        let x_closed = planar_meniscus_x(best.z, gamma, th, rho, g);
        close(best.x, x_closed, 2e-3 * x_closed.max(1.0), "meniscus profile");
    }
}

// ---------------------------------------------------------------------------
// 2. Texture homogenisation
// ---------------------------------------------------------------------------

#[test]
fn square_post_geometry_and_fold() {
    let p = SquarePosts::new(10e-6, 3e-6, 8e-6);
    close(p.solid_fraction(), 0.09, 1e-12, "φ_s");
    close(p.roughness(), 1.0 + 4.0 * 3e-6 * 8e-6 / (10e-6f64).powi(2), 1e-12, "r");
    let th = 110.0_f64.to_radians();
    close(p.wenzel_cos(th), p.roughness() * th.cos(), 1e-12, "wenzel");
    close(p.cassie_cos(th), p.solid_fraction() * (th.cos() + 1.0) - 1.0, 1e-12, "cassie");
    let threshold = p.cassie_transition_cos();
    for &deg in &[20.0_f64, 90.0, 150.0] {
        let t = deg.to_radians();
        let lhs = t.cos() < threshold;
        assert_eq!(lhs, p.cassie_is_stable(t), "transition at {deg}");
    }
}

#[test]
fn hex_post_geometry() {
    let p = HexPosts::new(20e-6, 3e-6, 10e-6);
    close(p.solid_fraction(), 3.0 * 9e-12 / 400e-12, 1e-12, "hex φ_s");
    let expect_r = 1.0 + 4.0 * 3.0f64.sqrt() * 3e-6 * 10e-6 / (20e-6f64).powi(2);
    close(p.roughness(), expect_r, 1e-12, "hex r");
}

#[test]
fn impalement_pressure_and_size() {
    let p = SquarePosts::new(10e-6, 3e-6, 8e-6);
    let th = 110.0_f64.to_radians();
    let pimp = p.impalement_pressure(WATER_GAMMA, th);
    close(pimp, -2.0 * WATER_GAMMA * th.cos() / p.gap(), 1e-12, "p_impale");
    close(p.critical_drop_radius(th), p.gap() / th.cos().abs(), 1e-12, "R_c");
    close(
        drop_laplace_pressure(WATER_GAMMA, p.critical_drop_radius(th)),
        pimp.abs(),
        1e-9,
        "critical balance",
    );
}

#[test]
fn cassie_sag_and_impalement_criterion() {
    // Sagitta of the pinned bridging meniscus at the impalement pressure:
    // R_c = gamma/|p_imp| = w/(2|cos theta_Y|), s = R_c - sqrt(R_c^2 - (w/2)^2).
    let p = SquarePosts::new(10e-6, 3e-6, 8e-6);
    let th = 110.0_f64.to_radians();
    let rc = p.gap() / (2.0 * th.cos().abs());
    let expect = rc - (rc * rc - (p.gap() / 2.0).powi(2)).sqrt();
    close(p.sag_before_impalement(WATER_GAMMA, th), expect, 1e-15, "sag");
    // For this geometry the critical sag is far below the post height, so the
    // Cassie state is not lost by sag; impalement needs the Laplace pressure
    // criterion, p_drop = 2 gamma / R_c = |p_imp|.
    assert!(p.sag_before_impalement(WATER_GAMMA, th) < p.height, "sag exceeds h");
}

// ---------------------------------------------------------------------------
// 3. Validity of the fold (direct free-energy minimisation on the texture)
// ---------------------------------------------------------------------------

#[test]
fn flat_wall_recovers_young_angle() {
    // On a flat substrate the minimiser must reproduce θ_Y for every size,
    // including obtuse angles (the arc is solved with the exact cap area).
    let flat = SquareWave::flat();
    for &deg in &[30.0_f64, 60.0, 90.0, 110.0, 120.0, 140.0] {
        let th = deg.to_radians();
        let d = lattice::minimise(Branch::Wenzel, &flat, th.cos(), 1000.0, 60.0, 3000);
        assert!(
            (d.angle - th).abs() < 5e-3,
            "flat θ_Y={deg}: apparent {:.3} deg",
            d.angle.to_degrees()
        );
    }
}

#[test]
fn wenzel_fold_converges_with_size() {
    // Square posts, θ_Y = 110°: the apparent angle must approach the Wenzel
    // angle as the base radius over the period grows, and the error must be
    // smaller at the larger size.
    let wave = SquareWave::new(10.0, 4.0, 3.0);
    let costh = 110.0f64.to_radians().cos();
    let pred = lattice::wenzel_cos(&wave, costh).acos();
    let small = lattice::minimise(Branch::Wenzel, &wave, costh, 1600.0, 120.0, 4000);
    let large = lattice::minimise(Branch::Wenzel, &wave, costh, 25600.0, 240.0, 8000);
    let e_small = (small.angle - pred).abs();
    let e_large = (large.angle - pred).abs();
    assert!(lattice::size_ratio(&large, &wave) > 4.0, "large drop not large");
    assert!(e_large < e_small, "error grew: {e_small:.3} -> {e_large:.3}");
    assert!(e_large < 1.0f64.to_radians(), "large-drop error {:.3} deg", e_large.to_degrees());
    assert!(e_small > 2.0f64.to_radians(), "small drop unexpectedly converged");
}

#[test]
fn cassie_fold_converges_with_size() {
    let wave = SquareWave::new(10.0, 4.0, 3.0);
    let costh = 110.0f64.to_radians().cos();
    let pred = lattice::cassie_cos(&wave, costh).acos();
    let small = lattice::minimise(Branch::Cassie, &wave, costh, 1600.0, 120.0, 4000);
    let large = lattice::minimise(Branch::Cassie, &wave, costh, 102400.0, 400.0, 8000);
    let e_small = (small.angle - pred).abs();
    let e_large = (large.angle - pred).abs();
    assert!(e_large < e_small, "cassie error grew: {e_small:.3} -> {e_large:.3}");
    assert!(e_large < 1.5f64.to_radians(), "cassie error {:.3} deg", e_large.to_degrees());
}
