//! Ray tracing and materials against catalog data, exact optics and
//! measured human ocular chromatic aberration.

use waveoptics::materials::*;
#[allow(unused_imports)]
use waveoptics::materials::water as _water;
use waveoptics::raytrace::*;

const INF: f64 = f64::INFINITY;

#[test]
fn sellmeier_glasses_reproduce_catalog_nd_and_abbe_numbers() {
    // Schott catalog: N-BK7 517642, N-BAF10 670471, SF10 728284.
    for (name, m, nd, vd) in [("N-BK7", Material::n_bk7(), 1.51680, 64.17), ("N-BAF10", Material::n_baf10(), 1.67003, 47.11), ("SF10", Material::sf10(), 1.72825, 28.41)] {
        let n = m.index(587.56e-9);
        let v = abbe(&m);
        println!("{name}: n_d {n:.5} (catalog {nd}), V_d {v:.2} (catalog {vd})");
        assert!((n - nd).abs() < 3e-5 && (v - vd).abs() < 0.1, "{name}");
    }
}

fn singlet(r1: f64, r2: f64, t: f64, image_z: f64) -> System {
    System {
        object_medium: Material::Vacuum,
        surfaces: vec![
            Surface { radius: r1, conic: 0.0, semi_diameter: 12.7e-3, thickness: t, after: Material::n_bk7() },
            Surface { radius: r2, conic: 0.0, semi_diameter: 12.7e-3, thickness: image_z - t, after: Material::Vacuum },
        ],
        stop_z: 0.0,
        stop_radius: 12.0e-3,
        image_z,
    }
}

#[test]
fn paraxial_power_matches_the_thick_lens_formula() {
    let (r1, r2, t) = (51.5e-3, -51.5e-3, 5.3e-3);
    let s = singlet(r1, r2, t, 0.1);
    for lam in [450e-9, 587.56e-9, 700e-9] {
        let n = Material::n_bk7().index(lam);
        let power = (n - 1.0) * (1.0 / r1 - 1.0 / r2 + (n - 1.0) * t / (n * r1 * r2));
        let p = s.paraxial(lam);
        assert!((p.efl - 1.0 / power).abs() < 1e-12, "{} vs {}", p.efl, 1.0 / power);
    }
}

#[test]
fn cartesian_ellipsoid_focuses_a_plane_wave_perfectly() {
    // A surface from vacuum into glass with conic −1/n² images an axial
    // plane wave without aberration at n R/(n − 1): the traced wavefront
    // aberration must vanish (Fermat), while a sphere shows spherical
    // aberration.
    let lam = 587.56e-9;
    let n = Material::n_bk7().index(lam);
    let r = 20e-3;
    let focus = n * r / (n - 1.0);
    let mk = |conic: f64| System {
        object_medium: Material::Vacuum,
        surfaces: vec![Surface { radius: r, conic, semi_diameter: 10e-3, thickness: focus, after: Material::n_bk7() }],
        stop_z: 0.0,
        stop_radius: 6e-3,
        image_z: focus,
    };
    let (_, w) = mk(-1.0 / (n * n)).wavefront(lam, 64);
    let worst = w.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
    println!("Cartesian ellipsoid: max |W| = {worst:.2e} m ({:.1e} waves)", worst / lam);
    assert!(worst < 1e-12);
    let (_, ws) = mk(0.0).wavefront(lam, 64);
    let edge = ws.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
    println!("sphere: max |W| = {:.2} waves", edge / lam);
    assert!(edge / lam > 1.0);
}

/// Thorlabs AC254-050-A (Thorlabs Zemax file): N-BAF10/SF10 achromat.
fn achromat() -> System {
    let image_z = 9.0e-3 + 2.5e-3 + 43.22018758686e-3;
    System {
        object_medium: Material::Vacuum,
        surfaces: vec![
            Surface { radius: 1.0 / 2.999400119976e1, conic: 0.0, semi_diameter: 12.7e-3, thickness: 9.0e-3, after: Material::n_baf10() },
            Surface { radius: 1.0 / -4.48833034111311e1, conic: 0.0, semi_diameter: 12.7e-3, thickness: 2.5e-3, after: Material::sf10() },
            Surface { radius: 1.0 / -3.4355996839248e0, conic: 0.0, semi_diameter: 12.7e-3, thickness: 43.22018758686e-3, after: Material::Vacuum },
        ],
        stop_z: 0.0,
        stop_radius: 12.7e-3,
        image_z,
    }
}

#[test]
fn achromat_reproduces_the_vendor_design_and_cancels_colour() {
    let a = achromat();
    // The Zemax file's curvatures are in 1/mm: radii (m) = 1e-3 / curv.
    let p = a.paraxial(587.56e-9);
    println!("AC254-050-A: EFL {:.3} mm, BFD {:.3} mm (Zemax image distance 43.220 mm)", p.efl * 1e3, p.bfd * 1e3);
    assert!((p.efl - 50.0e-3).abs() < 0.3e-3, "EFL {}", p.efl);
    assert!((p.bfd - 43.22e-3).abs() < 0.3e-3, "BFD {}", p.bfd);
    // Chromatic focal shift F−C: achromat vs a BK7 singlet of similar power.
    let shift = |s: &System| (s.paraxial(486.13e-9).bfd - s.paraxial(656.27e-9).bfd).abs();
    let sg = singlet(51.5e-3, -51.5e-3, 5.3e-3, 0.06);
    println!("F–C focal shift: achromat {:.1} µm, BK7 singlet {:.1} µm", shift(&a) * 1e6, shift(&sg) * 1e6);
    assert!(shift(&a) < 0.25 * shift(&sg));
    // Paraxial focus equals the exact trace of a near-axis ray.
    for nm in [450.0, 587.56, 700.0] {
        let l = nm * 1e-9;
        let r = a.trace(l, Ray { p: [1e-5, 0.0, -1e-3], d: [0.0, 0.0, 1.0], opl: 0.0 }).unwrap();
        let z = r.p[2] - r.p[0] / r.d[0] * r.d[2] - 11.5e-3;
        assert!((z - a.paraxial(l).bfd).abs() < 1e-8, "{nm} nm");
    }
}

use waveoptics::eye::*;

#[test]
fn model_eye_power_and_chromatic_aberration_match_human_data() {
    let nav = eye(3e-3, EyeDispersion::Navarro);
    // The Cauchy interpolation reproduces the published media indices to
    // their stated precision (four decimals, ±5e-5 rounding; a three-term
    // least-squares fit through four rounded values leaves up to ~6e-5).
    for (m, tab) in [(nav.surfaces[0].after, [1.3828, 1.3777, 1.3760, 1.3747]), (nav.surfaces[2].after, [1.4292, 1.4222, 1.4200, 1.4183])] {
        for (l, n) in [458e-9, 543e-9, 589.3e-9, 632.8e-9].iter().zip(tab) {
            assert!((m.index(*l) - n).abs() < 1e-4);
        }
    }
    let p543 = nav.paraxial(543e-9);
    let power = 1.0 / p543.efl;
    let focus = nav.surfaces.iter().take(3).map(|s| s.thickness).sum::<f64>() + p543.bfd;
    println!("eye: power {power:.2} D, paraxial focus − retina at 543 nm {:.3} mm, EP {:.2} mm behind the cornea, radius {:.3} mm", (focus - nav.image_z) * 1e3, p543.ep_z * 1e3, p543.ep_radius * 1e3);
    assert!((power - 60.4).abs() < 0.5, "power {power}");
    assert!((focus - nav.image_z).abs() < 0.05e-3, "not emmetropic");
    let phys = eye(3e-3, EyeDispersion::Physiological);
    println!("physiological dispersion scale {:.3}", physiological_scale());
    let mut worst_as: f64 = 0.0;
    let mut worst_thibos: f64 = 0.0;
    for nm in [420.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0] {
        let l = nm * 1e-9;
        let human = lca_atchison_smith(l) - lca_atchison_smith(589.3e-9);
        let thibos = lca_thibos(l) - lca_thibos(589.3e-9);
        let (mn, mp) = (model_lca(&nav, l), model_lca(&phys, l));
        println!("  {nm} nm: Navarro {mn:+.3} D, physiological {mp:+.3} D | Atchison & Smith {human:+.3} D, Thibos {thibos:+.3} D");
        if nm >= 450.0 {
            worst_as = worst_as.max((mp - human).abs());
            worst_thibos = worst_thibos.max((mp - thibos).abs());
        }
    }
    println!("physiological eye, 450–700 nm: worst |Δ| vs Atchison & Smith {worst_as:.3} D, vs Thibos {worst_thibos:.3} D");
    assert!(worst_as < 0.15 && worst_thibos < 0.2);
}
