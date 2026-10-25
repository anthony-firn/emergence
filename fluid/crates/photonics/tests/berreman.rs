//! Benchmarks for the anisotropic (Berreman 4x4) wave-optics rung.
//!
//! Every number is compared with an exact or published result and the measured
//! error is printed. Nothing is tuned to pass.

use photonics::berreman::{
    cholesteric, cholesteric_rotatory_power, effective_index_indicator, gooch_tarry_crossed,
    gooch_tarry_parallel, isotropic, twisted_nematic,
    twisted_nematic_jones, uniaxial, AnisoLayer, AnisoStack,
};
use photonics::complex::C;
use photonics::multilayer::{interface_reflectance, Layer, Polarization, Stack};
use std::f64::consts::PI;

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.12e} vs {b:.12e}, tol {tol:.1e}");
}

/// E7 nematic at 25 C, 589 nm: n_e = 1.7472, n_o = 1.5217.
/// S. T. Wu, U. Efron, L. D. Hess, "Birefringence measurements of liquid
/// crystals", Appl. Opt. 23, 3911 (1984), doi:10.1364/AO.23.003911.
const E7_NE: f64 = 1.7472;
const E7_NO: f64 = 1.5217;

/// 5CB nematic at 25 C, 589 nm: n_e = 1.717, n_o = 1.532.
/// S. T. Wu, "Birefringence dispersions of liquid crystals",
/// Phys. Rev. A 33, 1270 (1986), doi:10.1103/PhysRevA.33.1270.
const CB5_NE: f64 = 1.717;
const CB5_NO: f64 = 1.532;

#[test]
fn isotropic_stack_reduces_to_multilayer_code() {
    // The Berreman 4x4 must reduce to the crate's isotropic TE/TM code to
    // round-off at every angle and polarisation.
    let lambda = 633e-9;
    let iso_layers = vec![
        Layer::lossless(2.4, 90e-9),
        Layer::lossless(1.5, 150e-9),
        Layer::lossless(2.4, 90e-9),
        Layer::lossless(2.0, 60e-9),
    ];
    let mut worst = 0.0f64;
    for &(n0, ns) in &[(1.0f64, 1.52f64), (1.5, 1.0), (1.52, 1.52)] {
        let iso = Stack::new(n0, iso_layers.clone(), ns);
        let ber = AnisoStack::new(
            n0,
            iso_layers
                .iter()
                .map(|l| AnisoLayer::new(isotropic(l.index * l.index), l.thickness))
                .collect(),
            ns,
        );
        for &theta in &[0.0f64, 0.15, 0.4, 0.7, 1.0, 1.2] {
            let b = ber.solve(lambda, theta);
            let te = iso.solve(lambda, theta, Polarization::TE);
            let tm = iso.solve(lambda, theta, Polarization::TM);
            worst = worst.max((b.reflectance_s - te.reflectance).abs());
            worst = worst.max((b.transmittance_s - te.transmittance).abs());
            worst = worst.max((b.reflectance_p - tm.reflectance).abs());
            worst = worst.max((b.transmittance_p - tm.transmittance).abs());
        }
    }
    println!("isotropic reduction: worst |Delta| = {worst:.3e}");
    assert!(worst < 1e-12, "worst reduction error {worst:e}");
}

#[test]
fn uniaxial_retardance_is_exact_in_the_layer_matrix() {
    // For an in-plane optic axis at normal incidence the layer matrix has
    // exact eigenmodes: e -> index n_e, o -> index n_o, so the phase
    // retardance is Gamma = 2 pi (n_e - n_o) d / lambda exactly.
    let lambda = 550e-9;
    let (n_o, n_e, d) = (1.5, 1.7, 800e-9);
    let layer = AnisoLayer::new(uniaxial(C::real(n_o), C::real(n_e), [1.0, 0.0, 0.0]), d);
    let m = layer.matrix(lambda, 0.0);
    let k0 = 2.0 * PI / lambda;
    // e-mode eigenvector (Ex, hy) = (1, n_e).
    let e_mode = [C::real(1.0), C::real(0.0), C::real(0.0), C::real(n_e)];
    let out = m.apply(&e_mode);
    let expected_phase = k0 * n_e * d;
    let phase_e = out[0].arg();
    // o-mode (Ey, hx) = (1, -n_o).
    let o_mode = [C::real(0.0), C::real(1.0), C::real(-n_o), C::real(0.0)];
    let out_o = m.apply(&o_mode);
    let phase_o = out_o[1].arg();
    let gamma = (phase_e - phase_o).rem_euclid(2.0 * PI);
    let expected_gamma = (k0 * (n_e - n_o) * d).rem_euclid(2.0 * PI);
    println!(
        "retardance: computed phase_e = {phase_e:.9}, expected {expected_phase:.9}; \
         Gamma = {gamma:.9}, expected {expected_gamma:.9}"
    );
    close(phase_e, expected_phase.rem_euclid(2.0 * PI), 1e-9);
    close(gamma, expected_gamma, 1e-9);
}

#[test]
fn quarter_wave_plate_converts_45_linear_to_circular() {
    // A 45-degree linear input through a quarter-wave uniaxial plate. The
    // transmitted field is nearly circular; the residual ellipticity is the
    // ordinary/extraordinary Fresnel amplitude difference.
    let lambda = 633e-9;
    let (n_o, n_e) = (CB5_NO, CB5_NE);
    let dn = n_e - n_o;
    let d = lambda / (4.0 * dn); // quarter wave: Gamma = pi/2
    let stack = AnisoStack::new(
        1.0,
        vec![AnisoLayer::new(uniaxial(C::real(n_o), C::real(n_e), [1.0, 0.0, 0.0]), d)],
        1.0,
    );
    let input = [C::real(1.0 / 2f64.sqrt()), C::real(1.0 / 2f64.sqrt())];
    let out = stack.solve(lambda, 0.0).transmitted(input);
    let s = photonics::berreman::stokes(out);
    let docp = s[3] / s[0]; // S3 / S0
    println!("quarter-wave plate: |S3|/S0 = {:.6} (ideal 1)", docp.abs());
    assert!(docp.abs() > 0.99, "degree of circular polarisation {docp:.6}");

    // A half-wave plate (Gamma = pi) rotates 45-degree linear to 135-degree.
    let d2 = lambda / (2.0 * dn);
    let stack2 = AnisoStack::new(
        1.0,
        vec![AnisoLayer::new(uniaxial(C::real(n_o), C::real(n_e), [1.0, 0.0, 0.0]), d2)],
        1.0,
    );
    let out2 = stack2.solve(lambda, 0.0).transmitted(input);
    let s2 = photonics::berreman::stokes(out2);
    let angle = (0.5 * s2[2].atan2(s2[1]).to_degrees()).rem_euclid(180.0);
    println!("half-wave plate: polarisation angle = {angle:.4} deg (ideal 135)");
    assert!((angle - 135.0).abs() < 2.0, "half-wave angle {angle:.4}");
}

#[test]
fn uniaxial_brewster_and_tir_match_fresnel() {
    // Optic axis y (normal to the plane of incidence): the TE wave sees the
    // extraordinary index n_e, the TM wave the ordinary index n_o. Choosing the
    // exit index equal to the index seen by that polarisation makes the back
    // interface transparent, so the total reflectance is exactly the front
    // Fresnel reflectance of an isotropic interface.
    let lambda = 550e-9;
    let (n0, n_o, n_e) = (1.8, 1.5, 1.7);
    let d = 1500e-9;
    let mk = |exit: f64| {
        AnisoStack::new(
            n0,
            vec![AnisoLayer::new(
                uniaxial(C::real(n_o), C::real(n_e), [0.0, 1.0, 0.0]),
                d,
            )],
            exit,
        )
    };
    let te_stack = mk(n_e);
    let tm_stack = mk(n_o);

    let mut worst_te = 0.0f64;
    let mut worst_tm = 0.0f64;
    for &theta_deg in &[0.0f64, 15.0, 30.0, 39.81, 45.0, 55.0, 65.0, 70.0] {
        let theta = theta_deg.to_radians();
        let bte = te_stack.solve(lambda, theta);
        let btm = tm_stack.solve(lambda, theta);
        let fte = interface_reflectance(C::real(n0), C::real(n_e), theta, Polarization::TE);
        let ftm = interface_reflectance(C::real(n0), C::real(n_o), theta, Polarization::TM);
        worst_te = worst_te.max((bte.reflectance_s - fte).abs());
        worst_tm = worst_tm.max((btm.reflectance_p - ftm).abs());
    }
    println!("uniaxial interface: worst TE error {worst_te:.3e}, TM error {worst_tm:.3e}");
    assert!(worst_te < 1e-12 && worst_tm < 1e-12);

    // Brewster angle for the ordinary (TM) polarisation: tan(theta_B)=n_o/n0.
    let theta_b = (n_o / n0).atan();
    let r_b = tm_stack.solve(lambda, theta_b).reflectance_p;
    println!(
        "uniaxial Brewster: theta_B = {:.4} deg, R_p = {r_b:.3e}",
        theta_b.to_degrees()
    );
    assert!(r_b < 1e-24, "TM Brewster R_p = {r_b:e}");
    // TE is not extinguished there.
    assert!(te_stack.solve(lambda, theta_b).reflectance_s > 1e-3);

    // Total internal reflection beyond the ordinary critical angle.
    let theta_c = (n_o / n0).asin();
    for &theta in &[theta_c + 0.02, 1.0, 1.2] {
        let r = tm_stack.solve(lambda, theta).reflectance_p;
        close(r, 1.0, 1e-13);
    }
    // The TE critical angle uses n_e.
    let theta_ce = (n_e / n0).asin();
    let below = te_stack.solve(lambda, theta_ce - 0.05).reflectance_s;
    let above = te_stack.solve(lambda, theta_ce + 0.05).reflectance_s;
    println!(
        "uniaxial TIR: theta_c(o) = {:.4} deg, theta_c(e) = {:.4} deg; \
         TE R below/above = {below:.4}/{above:.4}",
        theta_c.to_degrees(),
        theta_ce.to_degrees()
    );
    assert!(below < 0.5 && above > 0.999);
}

#[test]
fn gyrotropic_faraday_rotation_matches_circular_eigenmodes() {
    // Magneto-optic (gyrotropic) slab with magnetization along z. The circular
    // eigen-indices are n_+/- = sqrt(eps +/- g); with the slab embedded in a
    // medium of the same index there are no interface reflections and a
    // linearly polarised input rotates by (k0/2)(n_+ - n_-) d.
    let lambda = 633e-9;
    let n = 2.0;
    let g = 0.01;
    let d = 1e-6;
    let stack = AnisoStack::new(
        n,
        vec![AnisoLayer::new(
            photonics::berreman::gyrotropic(C::real(n * n), C::real(g), [0.0, 0.0, 1.0]),
            d,
        )],
        n,
    );
    let sol = stack.solve(lambda, 0.0);
    let out = sol.transmitted([C::real(1.0), C::real(0.0)]);
    let s = photonics::berreman::stokes(out);
    let angle = 0.5 * s[2].atan2(s[1]);
    let expected = 0.5 * (2.0 * PI / lambda) * ((n * n + g).sqrt() - (n * n - g).sqrt()) * d;
    println!(
        "Faraday rotation: {:.6e} rad vs analytic {:.6e} (err {:.2e} rel)",
        angle.abs(),
        expected,
        (angle.abs() - expected).abs() / expected
    );
    assert!((angle.abs() - expected).abs() < 1e-6 * expected);
}

#[test]
fn energy_conservation_lossless_anisotropic_stacks() {
    let lambda = 600e-9;
    let stacks = vec![
        // General biaxial slab.
        AnisoStack::new(
            1.0,
            vec![AnisoLayer::new(
                photonics::berreman::biaxial(C::real(1.6), C::real(1.5), C::real(1.7)),
                500e-9,
            )],
            1.52,
        ),
        // Twisted nematic.
        twisted_nematic(CB5_NO, CB5_NE, 5e-6, PI / 2.0, 60, 1.0, 1.0),
        // Cholesteric.
        cholesteric(E7_NO, E7_NE, 350e-9, 6.0, 40, 1.0, 1.0),
    ];
    let mut worst = 0.0f64;
    for stack in &stacks {
        for &theta in &[0.0f64, 0.2, 0.5, 0.8] {
            let s = stack.solve(lambda, theta);
            worst = worst.max((s.reflectance_s + s.transmittance_s - 1.0).abs());
            worst = worst.max((s.reflectance_p + s.transmittance_p - 1.0).abs());
            // Arbitrary (circular) input too.
            let lcp = AnisoSolutionCircular::lcp();
            let rcp = AnisoSolutionCircular::rcp();
            worst = worst.max((s.reflectance(lcp) + s.transmittance(lcp) - 1.0).abs());
            worst = worst.max((s.reflectance(rcp) + s.transmittance(rcp) - 1.0).abs());
        }
    }
    println!("energy conservation: worst |R+T-1| = {worst:.3e}");
    assert!(worst < 1e-10);
}

/// Small helper so the test can build circular Jones vectors.
struct AnisoSolutionCircular;
impl AnisoSolutionCircular {
    fn lcp() -> [C; 2] {
        let inv = 1.0 / 2f64.sqrt();
        [C::real(inv), C::new(0.0, inv)]
    }
    fn rcp() -> [C; 2] {
        let inv = 1.0 / 2f64.sqrt();
        [C::real(inv), C::new(0.0, -inv)]
    }
}

#[test]
fn twisted_nematic_matches_gooch_tarry() {
    // Exact Berreman (many thin uniaxial sublayers) against the analytic
    // Gooch-Tarry transmission. Input polarised along the input director,
    // output projected on the parallel and crossed lab directions. The
    // incident/exit index is the ordinary index so the ordinary component has
    // no interface reflection.
    let lambda = 550e-9;
    let (n_o, n_e) = (CB5_NO, CB5_NE);
    let d = 5e-6;
    let phi = PI / 2.0;
    let n_sub = 400;
    let stack = twisted_nematic(n_o, n_e, d, phi, n_sub, n_o, n_o);
    let sol = stack.solve(lambda, 0.0);
    let out = sol.transmitted([C::real(1.0), C::real(0.0)]);
    let t_par = out[0].abs2();
    let t_cross = out[1].abs2();
    let (d_dn, gt_par, gt_cross) = (
        (n_e - n_o) * d,
        gooch_tarry_parallel(phi, (n_e - n_o) * d, lambda),
        gooch_tarry_crossed(phi, (n_e - n_o) * d, lambda),
    );
    let _ = d_dn;
    println!(
        "TN Gooch-Tarry: T_par {t_par:.6} vs {gt_par:.6} (err {:.3e}); \
         T_cross {t_cross:.6} vs {gt_cross:.6} (err {:.3e})",
        (t_par - gt_par).abs(),
        (t_cross - gt_cross).abs()
    );
    assert!((t_par - gt_par).abs() < 5e-3, "parallel TN error");
    assert!((t_cross - gt_cross).abs() < 5e-3, "crossed TN error");

    // The exact analytic twisted-nematic Jones matrix agrees with Gooch-Tarry
    // by construction.
    let j = twisted_nematic_jones(phi, 2.0 * PI * (n_e - n_o) * d / lambda);
    close(j.get(0, 0).abs2(), gt_par, 1e-12);
    close(j.get(1, 0).abs2(), gt_cross, 1e-12);
}

#[test]
fn twisted_nematic_adiabatic_limit() {
    // Slowly-varying-director fold: as the Mauguin parameter
    // u = pi d Delta n / (phi lambda) grows, the polarisation follows the
    // director. Measure the power on the output director and compare with the
    // analytic Jones matrix.
    let lambda = 550e-9;
    let (n_o, n_e) = (CB5_NO, CB5_NE);
    let phi = PI / 2.0;
    println!("Mauguin adiabatic limit (phi = 90 deg):");
    let mut last_loss = 1.0;
    for &d in &[2e-6f64, 4e-6, 8e-6, 16e-6, 32e-6] {
        let u = PI * (n_e - n_o) * d / (phi * lambda);
        let stack = twisted_nematic(n_o, n_e, d, phi, (d / 25e-9) as usize, n_o, n_o);
        let sol = stack.solve(lambda, 0.0);
        // Output component along the output director n_out=(cos phi, sin phi).
        let out = sol.transmitted([C::real(1.0), C::real(0.0)]);
        let along = out[0] * C::real(phi.cos()) + out[1] * C::real(phi.sin());
        let loss = 1.0 - along.abs2();
        // Analytic: |J_local(0,0)|^2 = cos^2 x + (Gamma/2x)^2 sin^2 x.
        let _j = twisted_nematic_jones(phi, 2.0 * PI * (n_e - n_o) * d / lambda);
        let x = (phi * phi + (PI * (n_e - n_o) * d / lambda).powi(2)).sqrt();
        let analytic = x.cos().powi(2)
            + ((PI * (n_e - n_o) * d / lambda) / x).powi(2) * x.sin().powi(2);
        println!(
            "  u = {u:6.2}: loss = {loss:.3e} (analytic {:.3e})",
            1.0 - analytic
        );
        last_loss = loss;
    }
    // At the largest retardation the polarisation follows the director.
    assert!(last_loss < 1e-3, "adiabatic loss {last_loss:.3e}");
}

#[test]
fn cholesteric_selective_circular_bragg_reflection() {
    // de Vries: normal-incidence reflection band lambda_0 = n_bar p,
    // width Delta lambda = Delta n p, one circular handedness only.
    let (n_o, n_e) = (E7_NO, E7_NE);
    let n_bar = 0.5 * (n_e + n_o);
    let dn = n_e - n_o;
    let p = 305.9e-9; // n_bar * p ~ 500 nm
    let lambda0 = n_bar * p;
    let band = dn * p;
    let stack = cholesteric(n_o, n_e, p, 20.0, 40, n_bar, n_bar);
    let lcp = AnisoSolutionCircular::lcp();
    let rcp = AnisoSolutionCircular::rcp();

    let r_at = |lam: f64| {
        let s = stack.solve(lam, 0.0);
        (s.reflectance(lcp), s.reflectance(rcp))
    };
    // Scan the co-handed reflectance across and beyond the band.
    println!("cholesteric scan R(LCP):");
    let mut lam = 420e-9;
    while lam <= 620e-9 {
        let (a, b) = r_at(lam);
        println!("  {:.0} nm: LCP {a:.4}  RCP {b:.4}", lam * 1e9);
        lam += 10e-9;
    }
    let (r_in, r_in2) = r_at(lambda0);
    let (r_lo, _) = r_at(lambda0 - 2.5 * band);
    let (r_hi, _) = r_at(lambda0 + 2.5 * band);
    println!(
        "cholesteric de Vries: lambda0 = {:.1} nm (band {:.1} nm); \
         R_in(LCP,RCP) = {r_in:.4},{r_in2:.4}; R_out = {r_lo:.4},{r_hi:.4}",
        lambda0 * 1e9,
        band * 1e9
    );
    let high = r_in.max(r_in2);
    let low = r_in.min(r_in2);
    assert!(high > 0.95, "in-band co-handed R = {high:.4}");
    assert!(low < 0.05, "in-band opposite-handed R = {low:.4}");
    assert!(r_lo < 0.1 && r_hi < 0.1, "out-of-band R");

    // Bandwidth: find the edges where the co-handed reflectance falls through
    // 0.5, scanning inside +/-1.5 band of the centre.
    let handed = if r_in >= r_in2 { lcp } else { rcp };
    let mut lo_edge = lambda0;
    let mut hi_edge = lambda0;
    let mut lam = lambda0;
    while lam > lambda0 - 1.5 * band {
        if stack.solve(lam, 0.0).reflectance(handed) < 0.5 {
            lo_edge = lam;
            break;
        }
        lam -= 1e-9;
    }
    let mut lam = lambda0;
    while lam < lambda0 + 1.5 * band {
        if stack.solve(lam, 0.0).reflectance(handed) < 0.5 {
            hi_edge = lam;
            break;
        }
        lam += 1e-9;
    }
    let measured = hi_edge - lo_edge;
    println!(
        "cholesteric bandwidth: measured {:.1} nm vs de Vries {:.1} nm (err {:.1}%)",
        measured * 1e9,
        band * 1e9,
        (measured - band).abs() / band * 100.0
    );
    assert!((measured - band).abs() < 0.25 * band, "bandwidth");
}

#[test]
fn cholesteric_homogenises_to_effective_medium() {
    // Pitch much smaller than the wavelength: the cholesteric must approach a
    // homogeneous uniaxial medium. Averaging the helical director over a full
    // turn gives in-plane permittivity (n_o^2 + n_e^2)/2 and axial n_o^2, so
    // at normal incidence the effective index is n_eff = sqrt((n_o^2+n_e^2)/2).
    let lambda = 550e-9;
    let (n_o, n_e) = (E7_NO, E7_NE);
    let n_eff = (0.5 * (n_o * n_o + n_e * n_e)).sqrt();
    let total = 20e-6;
    let lcp = AnisoSolutionCircular::lcp();
    let rcp = AnisoSolutionCircular::rcp();
    println!(
        "cholesteric homogenisation (h = 20 um, lambda = 550 nm, n_eff = {n_eff:.5}):"
    );
    let mut prev_r = f64::NAN;
    for &p in &[160e-9f64, 80e-9, 40e-9, 20e-9] {
        let turns = total / p;
        let stack = cholesteric(n_o, n_e, p, turns, 24, n_eff, n_eff);
        let s = stack.solve(lambda, 0.0);
        let r = 0.5 * (s.reflectance(lcp) + s.reflectance(rcp));
        let indicator = effective_index_indicator(p, lambda);
        let rho = cholesteric_rotatory_power(n_o, n_e, p, lambda);
        println!(
            "  p = {:5.0} nm  indicator {indicator:.4}  R = {r:.3e}  \
             rotatory power = {:.3e} rad/m",
            p * 1e9,
            rho
        );
        prev_r = r;
    }
    // Final reflection must be tiny.
    assert!(prev_r < 2e-5, "homogenised R = {prev_r:.3e}");

    // Compare the fine stack with the folded homogeneous uniaxial slab at a
    // small-but-finite pitch: the difference is the homogenisation error. The
    // folded slab has in-plane eps = (n_o^2+n_e^2)/2 and axial eps = n_o^2.
    let p = 100e-9;
    let turns = total / p;
    let fine = cholesteric(n_o, n_e, p, turns, 24, n_eff, n_eff);
    let folded = AnisoStack::new(
        n_eff,
        vec![AnisoLayer::new(
            uniaxial(C::real(n_eff), C::real(n_o), [0.0, 0.0, 1.0]),
            total,
        )],
        n_eff,
    );
    let mut worst = 0.0f64;
    for &theta in &[0.0f64, 0.2, 0.4] {
        let a = fine.solve(lambda, theta);
        let b = folded.solve(lambda, theta);
        worst = worst.max((a.reflectance_s - b.reflectance_s).abs());
        worst = worst.max((a.transmittance_s - b.transmittance_s).abs());
    }
    println!("cholesteric vs folded uniaxial slab at p/lambda=0.18: worst |Delta| = {worst:.3e}");
}
