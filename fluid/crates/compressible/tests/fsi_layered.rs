//! Fast validation of the compressible solver's fluid–solid interface against
//! exact layered acoustics (`phononic`). The model and experiment code live in
//! `examples/fsi_layered.rs` (included as a module) so the example and the
//! tests share one implementation.
//!
//! Every assertion states a measured bound, not a tuned number. See
//! `analysis/fsi.md` for the numbers and the convergence table.

#[path = "../examples/fsi_layered.rs"]
mod fsi;

use phononic::materials;

/// The shear-bearing `GprMix` extension must reproduce the single-material
/// longitudinal speed `sqrt(c0^2 + 4 c_s^2/3)`. (Its interface behaviour is
/// reported as a negative result; see `analysis/fsi.md`.)
#[test]
fn gprmix_solid_speed_matches_c_longitudinal() {
    let sys = fsi::mix(fsi::water(), fsi::solid(&materials::STEEL));
    let c = sys.material_c_l(1);
    let v = fsi::solid_pulse_speed(&sys, 800, 1.5e-4);
    assert!((v / c - 1.0).abs() < 5e-3, "measured {v:.1} m/s vs c_L {c:.1} m/s");
}

/// Two ideal gases at a stationary contact are an exact solution of the
/// multi-material system: the pressure must not move.
#[test]
fn ideal_gas_contact_is_exact() {
    let dev = fsi::contact_probe(256, 500);
    assert!(dev < 1e-10, "max |dp|/p = {dev:.3e}");
}

/// The crate's validated multi-material model (`FiveEq`, Allaire closure) on
/// the water/steel contrast; the peak reflection must be within 0.15 of the
/// exact Fresnel coefficient `(Z2-Z1)/(Z2+Z1)`.
#[test]
fn fiveeq_water_steel_half_space_reflection() {
    let (r, rf, _) = fsi::fiveeq_fresnel(256);
    assert!((r - rf).abs() < 0.15, "peak |r| {r:.4} vs Fresnel {rf:.4}");
}

/// Fluid–solid plate: the transferred spectra converge with resolution. The
/// error is set by the diffuse interface, which is a fixed number of cells
/// wide, so this is the resolution at which the 20 mm steel plate is resolved.
/// The transmitted amplitude is small (|T| ~ 0.06-0.1), so the absolute
/// incident-weighted error is the meaningful measure.
#[test]
fn plate_transmission_converges() {
    let (sys, media) = fsi::fiveeq_system();
    let zw = materials::WATER.rho * materials::WATER.c;
    let mut err = Vec::new();
    for n in [256usize, 512] {
        let cfg = fsi::plate_config(n);
        let exact = cfg.stack_media(&media, 0, true);
        let m = fsi::meas5(&sys, &media, 0, &cfg, 20_000.0, 120_000.0);
        let (_re, te) = fsi::exact_curve(&m.f, &exact, zw, zw);
        err.push(fsi::mag_abs_error(&m.t, &te, &m.inc));
    }
    assert!(err[1] < err[0], "|T| error did not fall: {err:?}");
    assert!(err[1] < 0.25, "|T| abs error {:.3} at N=512", err[1]);
}

/// The laminate fold. Below the first band gap the homogenised (Reuss) layer
/// reproduces the exact transmission; inside the gap the exact transmission
/// collapses while the homogenised medium predicts almost total transmission,
/// and the solver follows the exact one.
#[test]
fn laminate_fold_and_gap() {
    let (sys, media) = fsi::fiveeq_system_al_steel();
    let zs = materials::STEEL.rho * materials::STEEL.c_p;
    let cfg = fsi::laminate_config(768);
    let layers: Vec<(f64, f64, f64)> = cfg
        .layers
        .iter()
        .map(|&(m, d)| (media.rho[m], media.c[m], d))
        .collect();
    let (rho_e, c_e, d_e, _, _) = fsi::homogenised(&layers);
    let exact = cfg.stack_media(&media, 1, true);
    let mut homog = vec![exact[0]];
    homog.push((rho_e, c_e, d_e));
    homog.push(*exact.last().unwrap());

    let m = fsi::meas5(&sys, &media, 1, &cfg, 100_000.0, 700_000.0);
    let (_re, te) = fsi::exact_curve(&m.f, &exact, zs, zs);
    let (_rh, th) = fsi::exact_curve(&m.f, &homog, zs, zs);

    // below the gap (exact gap 278..522 kHz): the fold is valid
    let kb = fsi::nearest(&m.f, 150_000.0).unwrap();
    assert!(
        (m.t[kb].abs() - te[kb].abs()).abs() < 0.35,
        "below gap: measured {:.4}, exact {:.4}",
        m.t[kb].abs(),
        te[kb].abs()
    );
    assert!(
        (th[kb].abs() - te[kb].abs()).abs() < 0.35,
        "below gap: homogenised {:.4}, exact {:.4}",
        th[kb].abs(),
        te[kb].abs()
    );

    // inside the gap: the folded medium must fail, the solver must not
    let kg = fsi::nearest(&m.f, 400_000.0).unwrap();
    assert!(te[kg].abs() < 0.1, "exact in-gap |T| = {:.4}", te[kg].abs());
    assert!(th[kg].abs() > 0.5, "homogenised in-gap |T| = {:.4}", th[kg].abs());
    assert!(m.t[kg].abs() < 0.2, "solver in-gap |T| = {:.4}", m.t[kg].abs());
}
