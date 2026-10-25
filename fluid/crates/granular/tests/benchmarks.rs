//! Benchmarks for the granular / suspension rung.
//!
//! Every test checks a published or exact result and reports the error.
//! Runs are kept small (two-sphere collisions and a ~100-particle shear
//! cell for a few thousand steps); longer shear-cell and heap runs live in
//! `examples/large_runs.rs`.

use granular::contact::ContactParams;
use granular::dem::{hertz_contact_time, two_sphere_collision, ShearCell};
use granular::rheology::{granular_regime, Boyer2011, GranularRegime, Jop2006};
use granular::stokes::{
    drag_numerical, simple_shear_e, stresslet_exact, stresslet_numerical,
};
use granular::suspension::{
    batchelor_brownian, batchelor_green, einstein, kd_divergence_exponent,
    kd_second_coefficient, krieger_dougherty, BG_NONBROWNIAN_K2, INTRINSIC_VISCOSITY,
    PHI_M_RCP,
};
use std::f64::consts::PI;

fn velocity() -> [f64; 3] {
    [0.7, -0.3, 0.5]
}

// ---------------------------------------------------------------------------
// Suspension rung
// ---------------------------------------------------------------------------

#[test]
fn stokes_drag_from_the_exact_stress_field() {
    let v = velocity();
    let (a, mu) = (0.8, 1.3);
    let f = drag_numerical(v, a, mu, 24, 48);
    let vmag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let exact = 6.0 * PI * mu * a * vmag;
    let fmag = (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt();
    let err = (fmag - exact).abs() / exact;
    eprintln!("Stokes drag: F_num={fmag:.12e} F_exact={exact:.12e} rel.err={err:.3e}");
    assert!(err < 1e-10, "drag rel. error {err}");
    // The force is antiparallel to v.
    let cosang = -(f[0] * v[0] + f[1] * v[1] + f[2] * v[2]) / (fmag * vmag);
    assert!(cosang > 1.0 - 1e-9, "force not antiparallel: cos={cosang}");
}

#[test]
fn stresslet_matches_the_exact_single_sphere_result() {
    let (a, mu) = (0.8, 1.3);
    // A pure (traceless) extensional rate of strain.
    let mut e = [[0.0; 3]; 3];
    e[0][0] = 1.0;
    e[1][1] = -0.5;
    e[2][2] = -0.5;
    let exact = stresslet_exact(e, a, mu);
    let num = stresslet_numerical(e, a, mu, 24, 48);
    let mut max_err = 0.0_f64;
    let mut scale = 0.0_f64;
    for i in 0..3 {
        for j in 0..3 {
            max_err = max_err.max((exact[i][j] - num[i][j]).abs());
            scale = scale.max(exact[i][j].abs());
        }
    }
    eprintln!(
        "stresslet: max abs err {max_err:.3e} (scale {scale:.3e}, rel {:.3e})",
        max_err / scale
    );
    assert!(max_err / scale < 1e-10, "stresslet rel. err {}", max_err / scale);
}

#[test]
fn einstein_viscosity_from_the_stresslet_exponent() {
    // The dilute bulk stress assembled from exact single-sphere stresslets
    // must give exactly 1 + 5/2 phi. This is the numerics of the fold.
    let phi = 0.1;
    let a = 0.3;
    let e = simple_shear_e(1.0);
    let n = phi / ((4.0 / 3.0) * PI * a * a * a);
    let ratio = granular::suspension::dilute_bulk_stress_ratio(n, a, e);
    let eta_rel = ratio[0][1] / e[0][1];
    let err = (eta_rel - einstein(phi)).abs();
    eprintln!("Einstein: eta/eta0={eta_rel:.12} exact={:.12} err={err:.2e}", einstein(phi));
    assert!(err < 1e-12);
    // The coefficient of phi is 5/2.
    let numeric_coeff = (eta_rel - 1.0) / phi;
    assert!((numeric_coeff - INTRINSIC_VISCOSITY).abs() < 1e-12);
}

#[test]
fn batchelor_green_and_krieger_dougherty_series() {
    // Batchelor & Green (1972) non-Brownian k2 = 5.2; Batchelor (1977)
    // Brownian k2 = 6.17 (as compiled by Mendoza & Santamaría-Holek 2009).
    let phi = 1e-4;
    let k2_bg = (batchelor_green(phi) - 1.0 - 2.5 * phi) / (phi * phi);
    let k2_b = (batchelor_brownian(phi) - 1.0 - 2.5 * phi) / (phi * phi);
    eprintln!("BG k2 (numeric)={k2_bg:.6}  Brownian k2={k2_b:.6}");
    assert!((k2_bg - BG_NONBROWNIAN_K2).abs() < 1e-6);
    assert!((k2_b - 6.17).abs() < 1e-6);

    // Krieger–Dougherty with phi_m = 0.63 has an implied O(phi^2)
    // coefficient close to 5.2; the RCP value is the one that matches.
    let k2_kd = kd_second_coefficient(PHI_M_RCP, INTRINSIC_VISCOSITY);
    eprintln!("KD phi_m=0.63 implied k2={k2_kd:.4}");
    assert!((k2_kd - BG_NONBROWNIAN_K2).abs() < 0.15, "k2={k2_kd}");
    // Its divergence exponent is 2.5 phi_m ~ 1.58, weaker than the ~2 seen
    // near jamming; that is an honest miss, recorded in analysis/granular.md.
    let exp = kd_divergence_exponent(PHI_M_RCP, INTRINSIC_VISCOSITY);
    eprintln!("KD divergence exponent={exp:.4}");
    assert!((exp - 1.575).abs() < 1e-12);
}

#[test]
fn hard_sphere_dense_reference_points() {
    // Published low-shear hard-sphere relative viscosities (Phillies 2004
    // reanalysis, arXiv:cond-mat/0403305): a dynamic crossover at
    // phi~0.41, eta_r~11, and the melting point phi=0.494, eta_r~49.
    // Krieger–Dougherty with the RCP value is a non-Brownian/high-frequency
    // form; report the (large) miss honestly.
    let refs = [(0.41, 11.0), (0.494, 49.0)];
    for (phi, eta_ref) in refs {
        let kd = krieger_dougherty(phi, PHI_M_RCP, INTRINSIC_VISCOSITY);
        let err = (kd - eta_ref).abs() / eta_ref;
        eprintln!("KD(phi={phi})={kd:.3} vs published {eta_ref}: rel.err={err:.3}");
        assert!(kd > 1.0);
    }
    // The KD closure must at least reproduce the exact Einstein slope.
    let phi = 1e-5;
    let slope = (krieger_dougherty(phi, PHI_M_RCP, INTRINSIC_VISCOSITY) - 1.0) / phi;
    // Finite-difference slope; O(phi) correction is ~ k2*phi ~ 5e-5.
    assert!((slope - 2.5).abs() < 1e-3, "slope {slope}");
}

#[test]
fn suspension_indicators_switch_regimes() {
    use granular::suspension::{particle_reynolds, peclet, suspension_regime, SuspensionRegime};
    let (rho_f, a, eta_f, kt) = (1000.0, 1e-3, 1e-3, 4.1e-21);
    // A slow, small particle: dilute Newtonian.
    let re = particle_reynolds(rho_f, 0.01, a, eta_f);
    let pe = peclet(0.01, a, eta_f, kt);
    assert_eq!(
        suspension_regime(0.05, 0.63, re, pe),
        SuspensionRegime::DiluteNewtonian
    );
    // Fast flow: inertial and the Newtonian fold fails.
    let re_fast = particle_reynolds(rho_f, 1e5, a, eta_f);
    assert_eq!(
        suspension_regime(0.05, 0.63, re_fast, pe),
        SuspensionRegime::Inertial
    );
}

// ---------------------------------------------------------------------------
// Granular rung: contact
// ---------------------------------------------------------------------------

fn contact_params(restitution: f64) -> ContactParams {
    ContactParams {
        young: 1000.0,
        poisson: 0.3,
        mu: 0.4,
        mu_rolling: 0.0,
        restitution,
    }
}

#[test]
fn hertz_contact_time_matches_exact_formula() {
    let p = contact_params(1.0);
    for &v in &[0.01_f64, 0.05, 0.2] {
        let exact = hertz_contact_time(&p, 0.5, 1.0, v);
        let res = two_sphere_collision(&p, 0.5, 1.0, v, false);
        let err = (res.contact_time - exact).abs() / exact;
        eprintln!(
            "Hertz contact time v={v}: sim={:.6e} exact={exact:.6e} rel.err={err:.3}",
            res.contact_time
        );
        assert!(err < 5e-3, "v={v} rel. err {err}");
    }
}

#[test]
fn restitution_tracks_the_damping_target() {
    for &e_target in &[0.3_f64, 0.5, 0.7, 0.9] {
        let res = two_sphere_collision(&contact_params(e_target), 0.5, 1.0, 0.1, true);
        let err = (res.restitution - e_target).abs();
        eprintln!(
            "restitution target {e_target}: measured {:.4} err {err:.4}",
            res.restitution
        );
        // The Tsuji mapping is approximate for a Hertzian contact; keep a
        // loose but honest bound.
        assert!(err < 0.08, "e target {e_target} measured {}", res.restitution);
    }
}

// ---------------------------------------------------------------------------
// Granular rung: shear cell and mu(I)
// ---------------------------------------------------------------------------

fn shear_cell(n: usize, phi: f64, gdot: f64, dt: f64, seed: u64) -> ShearCell {
    ShearCell::new_cubic(
        n,
        phi,
        gdot,
        contact_params(0.5),
        0.5,
        1.0,
        dt,
        seed,
    )
}

/// Measure `(I, mu)` at a few shear rates under constant pressure. Kept
/// private so the example can do the longer version.
fn measure_mu_i(n: usize, phi: f64, p0: f64, rates: &[f64], steps: usize, dt: f64) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for (k, &g) in rates.iter().enumerate() {
        let mut cell = shear_cell(n, phi, g, dt, 1000 + k as u64);
        cell.set_pressure(p0, 2.0);
        let (p, tau, mu, i, _tg) = cell.run(steps);
        eprintln!(
            "  gdot={g:.4} P={p:.4} tau={tau:.4} mu={mu:.4} I={i:.4} phi={:.4}",
            cell.volume_fraction()
        );
        out.push((i, mu));
    }
    out
}

#[test]
fn shear_cell_mu_increases_with_inertial_number() {
    // A small constant-pressure shear cell: 108 particles, P0 = 0.05. The
    // point is the trend and the fit, not glass-bead-precise numbers; the
    // interparticle friction (0.4) and stiffness set the numbers.
    let rates = [0.005_f64, 0.01, 0.02, 0.04, 0.08, 0.16];
    let data = measure_mu_i(108, 0.58, 0.05, &rates, 10000, 0.005);
    // mu must be positive and rise from the quasi-static plateau to the
    // collisional regime.
    assert!(data.iter().all(|&(_, m)| m > 0.0 && m.is_finite()));
    let first = data.first().unwrap().1;
    let last = data.last().unwrap().1;
    eprintln!("mu(I): first={first:.4} last={last:.4}");
    assert!(last > first, "mu did not increase: {data:?}");
    // The low-I plateau should sit near the interparticle friction.
    assert!((first - 0.4).abs() < 0.12, "low-I mu {first} far from mu_p=0.4");
    let (fit, rmse) = Jop2006::fit(&data);
    eprintln!(
        "DEM mu(I) fit: mu_s={:.4} mu_2={:.4} I0={:.4} rmse={rmse:.4}",
        fit.mu_s, fit.mu_2, fit.i0
    );
    // Plausible fit: increasing friction and a positive cross-over.
    assert!(fit.mu_2 > fit.mu_s);
    assert!(fit.i0 > 0.0 && fit.i0.is_finite());
    // The fitted quasi-static friction should be close to the particle
    // friction; the DEM glass-bead comparison is recorded in the write-up.
    assert!((fit.mu_s - 0.4).abs() < 0.15, "fitted mu_s {}", fit.mu_s);
    // Compare with the glass-bead closure (Jop et al. 2006), reporting the
    // difference honestly; stiffness/friction shift the fitted numbers.
    let glass = Jop2006::glass_beads();
    eprintln!(
        "vs Jop glass beads: mu_s {:.3} vs {:.3}, mu_2 {:.3} vs {:.3}, I0 {:.3} vs {:.3}",
        fit.mu_s, glass.mu_s, fit.mu_2, glass.mu_2, fit.i0, glass.i0
    );
    assert!((fit.mu_2 - glass.mu_2).abs() < 0.2, "mu_2 {} vs {}", fit.mu_2, glass.mu_2);
}

#[test]
fn granular_regimes_by_inertial_number() {
    assert_eq!(granular_regime(1e-4), GranularRegime::QuasiStatic);
    assert_eq!(granular_regime(0.01), GranularRegime::DenseFlow);
    assert_eq!(granular_regime(0.2), GranularRegime::Transition);
    assert_eq!(granular_regime(0.5), GranularRegime::Collisional);
}

#[test]
fn boyer_viscous_number_fold() {
    // The immersed-grain closure: mu(J) and phi(J) (Boyer et al. 2011).
    let b = Boyer2011::frictional();
    let j_low = 1e-6;
    let j_high = 1.0;
    assert!(b.mu(j_high) > b.mu(j_low));
    assert!(b.phi(j_high) < b.phi(j_low));
    // Frictionless particles jam at a higher packing.
    assert!(Boyer2011::frictionless().phi_m > b.phi_m);
    eprintln!(
        "Boyer: mu(1)={:.4} phi(1)={:.4} phi_m={:.3}",
        b.mu(j_high),
        b.phi(j_high),
        b.phi_m
    );
}
