//! Optical elements against the paraxial formulas, against real rays
//! traced by `waveoptics`, and against the one case ray optics makes
//! exact: a paraboloid.

use std::f64::consts::PI;

use cad::kits::optic::{Lens, Mirror, Prism, Sag, Window};
use cad::meshing::mesh_solid;

#[test]
fn the_sag_equation_and_its_slope_agree() {
    for surface in [
        Sag::sphere(0.05),
        Sag::sphere(-0.08),
        Sag::parabola(0.1),
        Sag::conic(0.06, -0.5),
        Sag::conic(0.06, 2.0),
        Sag { radius: 0.05, conic: 0.0, aspheric: vec![1.0e3, -2.0e6] },
        Sag::plane(),
    ] {
        let a = 0.01;
        let mut worst: f64 = 0.0;
        for i in 1..100 {
            let r = a * i as f64 / 100.0;
            let h = 1e-8;
            let fd = (surface.at(r + h).unwrap() - surface.at(r - h).unwrap()) / (2.0 * h);
            worst = worst.max((surface.slope(r).unwrap() - fd).abs());
        }
        assert!(worst < 1e-6, "{surface:?}: analytic slope vs differences {worst:e}");
        // A sphere's sag is exactly R - sqrt(R^2 - r^2).
        if surface.conic == 0.0 && surface.aspheric.is_empty() && surface.radius.is_finite() {
            let r = 0.008;
            let big = surface.radius;
            let exact = big.signum() * (big.abs() - (big * big - r * r).sqrt());
            assert!((surface.at(r).unwrap() - exact).abs() < 1e-15, "sphere sag");
        }
        // A paraboloid's is exactly r^2 / 2R.
        if surface.conic == -1.0 {
            let r = 0.008;
            assert!((surface.at(r).unwrap() - r * r / (2.0 * surface.radius)).abs() < 1e-18);
        }
    }
    // A sphere only reaches its own radius; a paraboloid reaches forever.
    assert!((Sag::sphere(0.05).max_radius() - 0.05).abs() < 1e-15);
    assert!(Sag::parabola(0.05).max_radius().is_infinite());
    assert!(Sag::sphere(0.05).at(0.06).is_none(), "beyond the sphere");
}

#[test]
fn a_paraboloid_has_no_spherical_aberration() {
    // The exact statement in ray optics: a paraboloid brings every ray
    // parallel to its axis to one point, R/2, whatever the height. This
    // checks the sag and its normal together.
    let m = Mirror::paraboloid(0.2, 0.05, 0.01);
    let f = m.focal_length();
    let mut worst: f64 = 0.0;
    for i in 1..=50 {
        let h = 0.05 * i as f64 / 50.0;
        let focus = m.marginal_focus(h).expect("reflects");
        worst = worst.max((focus - f).abs() / f);
    }
    assert!(worst < 1e-14, "paraboloid focus spread {worst:e} of f");

    // A sphere of the same vertex radius falls short, and the shortfall
    // grows as h^2: that is the longitudinal spherical aberration.
    let s = Mirror::spherical(0.2, 0.05, 0.01);
    let paraxial = s.focal_length();
    let mut errors = Vec::new();
    for h in [0.01, 0.02, 0.04] {
        let focus = s.marginal_focus(h).expect("reflects");
        assert!(focus < paraxial, "a spherical mirror's marginal focus must fall short");
        errors.push(paraxial - focus);
    }
    // Doubling h quadruples the aberration.
    let ratio = errors[2] / errors[1];
    println!("spherical mirror LSA at h = 0.01, 0.02, 0.04 m: {errors:?}");
    assert!((ratio - 4.0).abs() < 0.3, "LSA should scale as h^2: ratio {ratio:.3}");
    // And it vanishes as h -> 0.
    let tiny = s.marginal_focus(1e-5).unwrap();
    assert!((tiny - paraxial).abs() / paraxial < 1e-8, "paraxial limit");
}

#[test]
fn lens_focal_length_matches_a_paraxial_trace_in_waveoptics() {
    // The generated surfaces are handed to `waveoptics::raytrace`, whose
    // paraxial trace is an independent implementation of first-order
    // optics. It must reproduce the lensmaker's equation.
    use waveoptics::materials::Material;
    use waveoptics::raytrace::{Surface, System};

    // N-BK7 by its Sellmeier coefficients (Schott data sheet), as
    // `waveoptics` carries them.
    let glass = Material::Sellmeier {
        b: [1.039_612_12, 0.231_792_344, 1.010_469_45],
        c: [6.000_698_67e-15, 2.001_791_44e-14, 1.035_606_53e-10],
    };
    let lambda = 587.6e-9;
    let n = glass.index(lambda);

    for lens in [
        Lens::new(0.1, -0.1, 0.006, 0.0125),      // biconvex
        Lens::plano_convex(0.0515, 0.004, 0.0125), // plano-convex
        Lens::new(0.08, 0.3, 0.005, 0.0125),       // meniscus
        Lens::new(-0.1, 0.1, 0.003, 0.0125),       // biconcave
    ] {
        let f_formula = lens.focal_length(n, 1.0);
        let data = lens.surface_data();
        let surfaces = vec![
            Surface {
                radius: data[0].0,
                conic: data[0].1,
                semi_diameter: data[0].2,
                thickness: data[0].3,
                after: glass,
            },
            Surface {
                radius: data[1].0,
                conic: data[1].1,
                semi_diameter: data[1].2,
                thickness: 0.05,
                after: Material::Vacuum,
            },
        ];
        let system = System {
            object_medium: Material::Vacuum,
            surfaces,
            stop_z: 0.0,
            stop_radius: 0.001,
            image_z: 0.05,
        };
        let paraxial = system.paraxial(lambda);
        let rel = (paraxial.efl - f_formula).abs() / f_formula.abs();
        println!(
            "R1 = {:+.4} R2 = {:+.4}: lensmaker {f_formula:+.6} m, paraxial trace {:+.6} m ({rel:.3e})",
            lens.front.radius, lens.back.radius, paraxial.efl
        );
        assert!(rel < 1e-12, "focal length differs by {rel:e}");

        // The back focal distance from the same two formulas.
        let (_, bfd) = lens.focal_distances(n, 1.0);
        let rel = (paraxial.bfd - bfd).abs() / bfd.abs();
        assert!(rel < 1e-10, "back focal distance differs by {rel:e}");
    }
}

#[test]
fn a_real_ray_through_the_generated_lens_reaches_the_paraxial_focus() {
    // A ray at a small height must cross the axis at the paraxial focus:
    // this is the geometry of the surfaces as the tracer sees them, not a
    // formula.
    use waveoptics::materials::Material;
    use waveoptics::raytrace::{Ray, Surface, System};

    let glass = Material::Cauchy { a: 1.5168, b: 0.0, c: 0.0 };
    let lambda = 587.6e-9;
    let n = glass.index(lambda);
    let lens = Lens::plano_convex(0.0515, 0.004, 0.0125);
    let (_, bfd) = lens.focal_distances(n, 1.0);

    let data = lens.surface_data();
    let make = |image_z: f64| System {
        object_medium: Material::Vacuum,
        surfaces: vec![
            Surface { radius: data[0].0, conic: data[0].1, semi_diameter: data[0].2, thickness: data[0].3, after: glass },
            Surface { radius: data[1].0, conic: data[1].1, semi_diameter: data[1].2, thickness: image_z, after: Material::Vacuum },
        ],
        stop_z: 0.0,
        stop_radius: 0.0125,
        image_z: lens.thickness + image_z,
    };
    let system = make(bfd);

    // Where a ray at height h crosses the axis, from its traced state.
    let crossing = |h: f64| -> f64 {
        let ray = Ray { p: [h, 0.0, -0.01], d: [0.0, 0.0, 1.0], opl: 0.0 };
        let out = system.trace(lambda, ray).expect("ray must get through");
        // Advance to the axis.
        let t = -out.p[0] / out.d[0];
        out.p[2] + t * out.d[2]
    };
    let paraxial_focus = lens.thickness + bfd;
    let near_axis = crossing(1e-4);
    let rel = (near_axis - paraxial_focus).abs() / paraxial_focus;
    println!("ray at 0.1 mm crosses at {near_axis:.6} m, paraxial focus {paraxial_focus:.6} m ({rel:.3e})");
    assert!(rel < 1e-5, "a near-axis ray must reach the paraxial focus: {rel:e}");

    // A marginal ray focuses short: spherical aberration, growing as h^2.
    let mut lsa = Vec::new();
    for h in [0.003, 0.006, 0.012] {
        let z = crossing(h);
        assert!(z < paraxial_focus, "the marginal focus must be short of the paraxial one");
        lsa.push(paraxial_focus - z);
    }
    println!("longitudinal spherical aberration at h = 3, 6, 12 mm: {lsa:?}");
    let ratio = lsa[2] / lsa[1];
    assert!((ratio - 4.0).abs() < 0.8, "LSA should grow near h^2: ratio {ratio:.3}");
}

#[test]
fn lens_solids_measure_their_own_sag() {
    let lens = Lens::new(0.1, -0.1, 0.006, 0.0125);
    let solid = lens.solid().expect("solid");
    // On the axis the glass runs from 0 to the centre thickness.
    assert!(solid.sdf([0.0, 0.0, 0.5 * lens.thickness]) < 0.0, "the axis must be inside the glass");
    assert!(solid.sdf([0.0, 0.0, -1e-4]) > 0.0 && solid.sdf([0.0, 0.0, lens.thickness + 1e-4]) > 0.0);
    // At the clear aperture the glass is only the edge thickness thick.
    let et = lens.edge_thickness().unwrap();
    let a = lens.semi_diameter;
    let sf = lens.front.at(a).unwrap();
    assert!(solid.sdf([a * 0.999, 0.0, sf + 0.5 * et]) < 0.0, "the edge must be solid");
    assert!(solid.sdf([a * 0.999, 0.0, sf - 1e-4]) > 0.0, "beyond the front surface");

    // Volume against the closed form: the glass is a cylinder of the edge
    // thickness plus a spherical cap at each end, and a cap of height h on
    // a sphere of radius R has volume pi h^2 (3R - h)/3.
    let cap = |radius: f64| {
        let big = radius.abs();
        let h = big - (big * big - a * a).sqrt();
        PI * h * h * (3.0 * big - h) / 3.0
    };
    let analytic = PI * a * a * et + cap(lens.front.radius) + cap(lens.back.radius);
    let v = lens.volume().unwrap();
    let rel = (v - analytic).abs() / analytic;
    println!("lens volume {v:e} m^3, closed form {analytic:e} ({rel:.3e})");
    assert!(rel < 1e-12, "volume differs from the spherical closed form by {rel:e}");

    // And the meshed solid agrees to the meshing error.
    let mesh = mesh_solid(&solid, lens.semi_diameter / 30.0).expect("mesh");
    assert!(mesh.report().watertight, "{:?}", mesh.report());
    let rel = (mesh.report().volume - v).abs() / v;
    assert!(rel < 0.01, "meshed volume off by {rel:.4}");

    // N-BK7 at 2510 kg/m^3: a 25 mm lens weighs a few grams.
    let mass = lens.mass(2510.0).unwrap();
    assert!(mass > 1e-3 && mass < 1e-1, "{:.4} kg", mass);

    // A lens that would close to a knife edge is refused.
    let knife = Lens::new(0.02, -0.02, 0.001, 0.0125);
    let err = knife.solid().unwrap_err();
    assert!(err.contains("knife edge"), "{err}");
}

#[test]
fn prism_minimum_deviation_is_where_the_trace_says_it_is() {
    let n = 1.5168; // N-BK7 at the d line
    let p = Prism::equilateral(0.025, 0.025);
    let formula = p.minimum_deviation(n).expect("passes");
    // Scan the incidence angle for the minimum of the traced deviation.
    let mut best = (f64::INFINITY, 0.0);
    for i in 0..200_001 {
        let i1 = 0.2 + (1.4 - 0.2) * i as f64 / 200_000.0;
        if let Some(d) = p.deviation(n, i1) {
            if d < best.0 {
                best = (d, i1);
            }
        }
    }
    println!(
        "equilateral prism, n = {n}: traced minimum {:.6} rad at i1 = {:.6}, formula {formula:.6}",
        best.0, best.1
    );
    assert!((best.0 - formula).abs() < 1e-6, "minimum deviation {:.8} vs {formula:.8}", best.0);
    // The minimum is at the symmetric passage.
    let symmetric = p.minimum_deviation_incidence(n).unwrap();
    assert!((best.1 - symmetric).abs() < 1e-4, "{:.6} vs {symmetric:.6}", best.1);
    // For an equilateral prism of N-BK7 that is 38.65 degrees:
    // 2 asin(1.5168 * sin 30) - 60.
    assert!((formula.to_degrees() - 38.647).abs() < 0.01, "{:.3} degrees", formula.to_degrees());

    // A right-angle prism of the same glass cannot pass a ray at all: its
    // apex exceeds 2 asin(1/n).
    let ra = Prism::right_angle(0.02, 0.02);
    assert!(ra.apex > Prism::maximum_apex(n));
    assert!(ra.minimum_deviation(n).is_none(), "the second face totally reflects");
    // Total internal reflection shows up in the trace as well.
    assert!(ra.deviation(n, 0.5).is_none());

    // The solid is a triangular extrusion of the right area.
    let solid = p.solid();
    let mesh = mesh_solid(&solid, 0.0008).expect("mesh");
    let half_base = 0.025 * (PI / 6.0).sin();
    let height = 0.025 * (PI / 6.0).cos();
    let expect = half_base * height * 0.025;
    let rel = (mesh.report().volume - expect).abs() / expect;
    assert!(rel < 0.02, "prism volume {:e} vs {expect:e}", mesh.report().volume);
}

#[test]
fn a_window_shifts_the_focus_by_the_published_amount() {
    let w = Window { thickness: 0.003, semi_diameter: 0.0125 };
    let n = 1.5168;
    // t (1 - 1/n), the standard plane-parallel-plate shift.
    assert!((w.focus_shift(n) - 0.003 * (1.0 - 1.0 / n)).abs() < 1e-18);
    assert!(w.focus_shift(n) > 0.0 && w.focus_shift(n) < w.thickness);
    // A window of index 1 shifts nothing.
    assert!(w.focus_shift(1.0).abs() < 1e-18);
    let mesh = mesh_solid(&w.solid(), 0.0004).expect("mesh");
    let expect = PI * 0.0125f64.powi(2) * 0.003;
    assert!((mesh.report().volume - expect).abs() / expect < 0.01);
}
