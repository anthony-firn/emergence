//! Image formation against closed-form Fourier optics.

use num_complex::Complex64 as Cx;
use std::f64::consts::PI;
use waveoptics::color::*;
use waveoptics::imaging::*;
use waveoptics::materials::Material;
use waveoptics::raytrace::{Surface, System};

/// Radial profile helpers on a square window centred on the axis.
fn window(n: usize, pitch: f64) -> Window {
    Window { center: [0.0; 2], pixels: [n, n], pitch }
}

#[test]
fn ideal_lens_forms_the_airy_pattern_with_exact_normalisation() {
    let (a, f, lam) = (2e-3, 50e-3, 550e-9);
    let p = PupilFunction::ideal(a, f, lam, 128, 0.0);
    let first_zero = 1.22 * lam * f / (2.0 * a);
    let w = window(121, first_zero / 20.0);
    let u = image_field(&p, &vec![Cx::new(1.0, 0.0); 128 * 128], &w);
    let i: Vec<f64> = u.iter().map(|v| v.norm_sqr()).collect();
    // peak: (power through the pupil)·(πa²)/(λ²f²) = (πa²)²/(λf)²; the
    // sampled pupil area is used so the check isolates the transform
    let open = p.values.iter().filter(|v| v.norm() > 0.0).count() as f64 * p.spacing().powi(2);
    let peak = (open / (lam * f)).powi(2);
    let centre = i[60 + 121 * 60];
    println!("Airy peak {centre:.6e} vs {peak:.6e}");
    assert!((centre - peak).abs() / peak < 1e-9);
    // first zero along x at 1.22 λ f / D (pixel 20 from the centre)
    let row: Vec<f64> = (0..121).map(|k| i[k + 121 * 60] / centre).collect();
    let min_k = (61..100).min_by(|&x, &y| row[x].total_cmp(&row[y])).unwrap();
    println!("first minimum at {} px (expected 20), I/I₀ = {:.2e}", min_k - 60, row[min_k]);
    assert_eq!(min_k - 60, 20);
    assert!(row[min_k] < 2e-3);
    // encircled energy within the first zero: 1 − J0² − J1² at 3.8317 = 0.8378
    let pitch = w.pitch;
    let mut inside = 0.0;
    let total_in = open; // incoming power (|U_in| = 1)
    let fine = window(241, first_zero / 40.0);
    let uf = image_field(&p, &vec![Cx::new(1.0, 0.0); 128 * 128], &fine);
    for q in 0..241 * 241 {
        let (x, y) = ((q % 241) as f64 - 120.0, (q / 241) as f64 - 120.0);
        if (x * x + y * y).sqrt() * fine.pitch <= first_zero {
            inside += uf[q].norm_sqr() * fine.pitch * fine.pitch;
        }
    }
    let _ = pitch;
    let ee = inside / total_in;
    println!("encircled energy within the first zero {ee:.4} (Airy 0.8378)");
    assert!((ee - 0.8378).abs() < 0.004);
}

#[test]
fn defocus_lowers_the_strehl_ratio_as_sinc_squared() {
    for w20 in [0.25, 0.5, 1.0] {
        let lam = 550e-9;
        let p = PupilFunction::ideal(1e-3, 20e-3, lam, 256, w20 * lam);
        let phi = 2.0 * PI * w20;
        let exact = ((phi / 2.0).sin() / (phi / 2.0)).powi(2);
        println!("defocus {w20} λ: Strehl {:.5} (exact {exact:.5})", p.strehl());
        assert!((p.strehl() - exact).abs() < 2e-3);
    }
}

#[test]
fn traced_optics_strehl_matches_marechal_and_perfect_systems() {
    let lam = 587.56e-9;
    let n = Material::n_bk7().index(lam);
    let r = 20e-3;
    let focus = n * r / (n - 1.0);
    let mk = |conic: f64, stop: f64| System {
        object_medium: Material::Vacuum,
        surfaces: vec![Surface { radius: r, conic, semi_diameter: 10e-3, thickness: focus, after: Material::n_bk7() }],
        stop_z: 0.0,
        stop_radius: stop,
        image_z: focus,
    };
    let perfect = PupilFunction::traced(&mk(-1.0 / (n * n), 4e-3), lam, 128);
    println!("Cartesian ellipsoid: Strehl {:.9}", perfect.strehl());
    assert!(perfect.strehl() > 1.0 - 1e-9);
    // a sphere stopped down to small aberration: Maréchal exp(−(2πσ/λ)²)
    let sph = PupilFunction::traced(&mk(0.0, 2.4e-3), lam, 128);
    let sigma = sph.rms_wavefront();
    let marechal = (-(2.0 * PI * sigma / lam).powi(2)).exp();
    println!("sphere: rms {:.3} λ, Strehl {:.4}, Maréchal {marechal:.4}", sigma / lam, sph.strehl());
    assert!(sigma / lam > 0.02 && sigma / lam < 0.1);
    assert!((sph.strehl() - marechal).abs() < 0.01);
}

#[test]
fn pinhole_goes_from_airy_to_geometric_shadow() {
    let lam = 550e-9;
    let f = 50e-3;
    // Fresnel number N = (D/2)²/(λ f): far field for N ≪ 1.
    let d = 40e-6;
    let p = PupilFunction::pinhole(d, f, lam, 64);
    let zero = 1.22 * lam * f / d;
    let w = window(81, zero / 20.0);
    let i: Vec<f64> = image_field(&p, &vec![Cx::new(1.0, 0.0); 64 * 64], &w).iter().map(|v| v.norm_sqr()).collect();
    let row: Vec<f64> = (0..81).map(|k| i[k + 81 * 40] / i[40 + 81 * 40]).collect();
    let min_k = (41..75).min_by(|&x, &y| row[x].total_cmp(&row[y])).unwrap();
    println!("pinhole N = {:.3}: first minimum at {} px (Airy 20)", (d / 2.0).powi(2) / (lam * f), min_k - 40);
    assert!((min_k as i32 - 60).abs() <= 1);
    // Large pinhole (N = 18): geometric shadow of diameter D, edge at half
    // intensity, interior near 1 with Fresnel ripples.
    let d = 2e-3;
    let p = PupilFunction::pinhole(d, f, lam, 512);
    let w = Window { center: [0.0; 2], pixels: [161, 1], pitch: 1.6 * d / 160.0 };
    let i: Vec<f64> = image_field(&p, &vec![Cx::new(1.0, 0.0); 512 * 512], &w).iter().map(|v| v.norm_sqr()).collect();
    let xs = w.coords(0);
    let at = |x: f64| i[xs.iter().position(|&v| (v - x).abs() < 0.5 * w.pitch).unwrap()];
    let inner: f64 = xs.iter().zip(&i).filter(|(x, _)| x.abs() < 0.3 * d).map(|(_, v)| *v).sum::<f64>() / xs.iter().filter(|x| x.abs() < 0.3 * d).count() as f64;
    println!("pinhole N = {:.1}: mean inner irradiance {inner:.3}, at the edge {:.3}, outside {:.4}", (d / 2.0).powi(2) / (lam * f), at(0.5 * d), at(0.75 * d));
    assert!((inner - 1.0).abs() < 0.05 && (at(0.5 * d) - 0.25).abs() < 0.12 && at(0.75 * d) < 0.02);
}

#[test]
fn cie_matrix_is_recovered_exactly_from_the_cone_fundamentals() {
    let (m, worst) = lms_to_xyz();
    println!("LMS→XYZ {:?}, worst residual {worst:.1e}", m);
    assert!(worst < 5e-6);
    // Y is a positive combination of L and M only (CIE 2006 construction).
    assert!(m[1][2].abs() < 1e-6 && m[1][0] > 0.0 && m[1][1] > 0.0);
}
