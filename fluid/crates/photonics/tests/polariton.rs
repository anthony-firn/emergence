//! Benchmarks for the phonon-polariton rung. Every number is checked against an
//! exact identity or a sourced handbook value; measured deviations are printed
//! so `analysis/polariton.md` can quote them. No tolerance is loosened to make a
//! case pass.

use photonics::polariton::{
    born_oscillator_strength, cm_to_rad, film_thickness_indicator, lambda_to_wavenumber,
    lorentz_single, polar_crystal_from_chain, rad_to_cm, spatial_dispersion_indicator,
    wavenumber_to_lambda, DiatomicChain, PolarCrystal, AMU,
};
use photonics::dispersion::Dispersion;
use std::f64::consts::PI;

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.12e} vs {b:.12e}, tol {tol:.1e}");
}

fn crystals() -> Vec<PolarCrystal> {
    vec![
        PolarCrystal::silicon_carbide(),
        PolarCrystal::sodium_chloride(),
        PolarCrystal::gallium_arsenide(),
        PolarCrystal::hbn_inplane(),
    ]
}

#[test]
fn lyddane_sachs_teller_and_handbook_static_permittivity() {
    println!("\n=== LST and static permittivity vs handbook ===");
    for c in crystals() {
        let m = &c.modes[0];
        let lst = (m.lo_cm / m.to_cm).powi(2);
        let model = c.static_epsilon();
        // The single-oscillator form satisfies LST exactly, by construction.
        close(model, c.eps_inf * lst, 1e-12);
        close(c.lst_ratio(), lst, 1e-12);
        let err = model - c.eps0_measured;
        println!(
            "{:14} eps_inf={:6.3} eps0(model)={:7.4} eps0(handbook)={:6.3} \
             rel.err={:+.3}% (LST ratio {:.5})",
            c.name,
            c.eps_inf,
            model,
            c.eps0_measured,
            100.0 * err / c.eps0_measured,
            lst
        );
        // The model static permittivity must reproduce the handbook value to a
        // few percent (the residual is the independent-measurement scatter).
        assert!(
            (err / c.eps0_measured).abs() < 0.035,
            "{}: static eps off by more than 3.5%",
            c.name
        );
    }
}

#[test]
fn lorentz_and_factorised_forms_agree_and_limits_hold() {
    // A single oscillator with equal damping in numerator and denominator.
    let m = PolarCrystal::silicon_carbide().modes[0];
    let gamma_cm = 3.5;
    let mut c = PolarCrystal::silicon_carbide();
    c.modes[0].gamma_to_cm = gamma_cm;
    c.modes[0].gamma_lo_cm = gamma_cm;
    for &nu in &[100.0f64, 500.0, 796.2, 850.0, 972.7, 1200.0, 5000.0] {
        let omega = cm_to_rad(nu);
        let factorised = c.epsilon(omega);
        let lorentz = lorentz_single(c.eps_inf, cm_to_rad(m.to_cm), cm_to_rad(m.lo_cm), cm_to_rad(gamma_cm), omega);
        close((factorised - lorentz).abs(), 0.0, 1e-12);
    }
    // High-frequency limit is eps_inf; static limit is eps(0).
    close(c.epsilon(cm_to_rad(1.0e6)).re, c.eps_inf, 1e-4);
    close(c.epsilon(0.0).re, c.static_epsilon(), 1e-10);
    close(c.epsilon(0.0).im, 0.0, 1e-12);
    // The pole sits at omega_TO: for the lossless form Re eps -> -inf just
    // above it and +inf just below, so its sign changes across the pole.
    let lossless = c.lossless();
    let wto = cm_to_rad(m.to_cm);
    let below = lossless.epsilon(wto - 1.0e9).re;
    let above = lossless.epsilon(wto + 1.0e9).re;
    println!(
        "pole check: Re eps(TO-)= {below:.3e}, Re eps(TO+)= {above:.3e}",
    );
    assert!(below > 0.0 && above < 0.0);
    // The numerator zero sits at omega_LO: the lossless permittivity vanishes
    // there.
    close(
        c.lossless().epsilon_wavenumber(m.lo_cm).abs(),
        0.0,
        1e-9,
    );
}

#[test]
fn lossless_bulk_reststrahlen_band_is_a_perfect_reflector() {
    println!("\n=== Bulk Reststrahlen reflectance (lossless) ===");
    for c in crystals() {
        let lossless = c.lossless();
        let m = lossless.modes[0];
        // Exactly R = 1 wherever eps(omega) < 0, i.e. between TO and LO.
        for &frac in &[0.1f64, 0.3, 0.5, 0.7, 0.9] {
            let nu = m.to_cm + frac * (m.lo_cm - m.to_cm);
            let r = lossless.bulk_reflectance(nu);
            close(r, 1.0, 1e-12);
        }
        // Just outside the band eps > 0 and R < 1.
        let r_below = lossless.bulk_reflectance(0.98 * m.to_cm);
        let r_above = lossless.bulk_reflectance(1.02 * m.lo_cm);
        let r_in_lo = lossless.bulk_reflectance(m.to_cm * 1.001);
        let r_in_hi = lossless.bulk_reflectance(m.lo_cm * 0.999);
        println!(
            "{:14} band {:.1}-{:.1} cm^-1: R(0.98 TO)={:.4}, R(TO+)={:.4}, \
             R(LO-)={:.4}, R(1.02 LO)={:.4}",
            c.name, m.to_cm, m.lo_cm, r_below, r_in_lo, r_in_hi, r_above
        );
        assert!(r_below < 1.0 && r_above < 1.0);
    }
}

#[test]
fn sic_reststrahlen_edges_and_peak_match_published_band() {
    // The task quotes the SiC Reststrahlen band as ~790-970 cm^-1. The band
    // edges are the TO and LO frequencies, where Re eps changes sign; use the
    // lossless form so they are exact. The damped half-power edges are reported
    // for comparison.
    let lossless = PolarCrystal::silicon_carbide().lossless();
    let m = lossless.modes[0];
    let (lo, hi) = (m.to_cm, m.lo_cm);
    // R = 1 exactly across the band (eps < 0) and < 1 on either side.
    close(lossless.bulk_reflectance(0.5 * (lo + hi)), 1.0, 1e-12);
    assert!(lossless.epsilon_wavenumber(lo * 1.001).re < 0.0);
    assert!(lossless.epsilon_wavenumber(hi * 0.999).re < 0.0);
    assert!(lossless.epsilon_wavenumber(lo * 0.999).re > 0.0);
    assert!(lossless.epsilon_wavenumber(hi * 1.001).re > 0.0);

    // Damped spectrum: peak reflectance and half-power edges.
    let c = PolarCrystal::silicon_carbide();
    let cross = |a: f64, b: f64| {
        let f = |x: f64| c.bulk_reflectance(x) - 0.5;
        let (mut a, mut b) = (a, b);
        for _ in 0..100 {
            let mid = 0.5 * (a + b);
            if f(a) * f(mid) <= 0.0 {
                b = mid;
            } else {
                a = mid;
            }
        }
        0.5 * (a + b)
    };
    let edge_low = cross(600.0, 0.5 * (lo + hi));
    let edge_high = cross(0.5 * (lo + hi), 1100.0);
    let mut peak = 0.0_f64;
    let mut peak_nu = 0.0;
    for i in 0..=4000 {
        let nu = lo + (hi - lo) * i as f64 / 4000.0;
        let r = c.bulk_reflectance(nu);
        if r > peak {
            peak = r;
            peak_nu = nu;
        }
    }
    println!(
        "\n=== SiC Reststrahlen band ===\n  published 790-970 cm^-1; \
         lossless band {lo:.1}-{hi:.1} cm^-1 (Re eps < 0); \
         errors {:+.1}/{:+.1} cm^-1\n  damped half-power edges \
         {edge_low:.1}-{edge_high:.1} cm^-1, peak R={peak:.4} at {peak_nu:.1} cm^-1",
        lo - 790.0,
        hi - 970.0
    );
    assert!((lo - 790.0).abs() < 12.0, "low edge {lo}");
    assert!((hi - 970.0).abs() < 12.0, "high edge {hi}");
    assert!(peak > 0.95, "peak {peak}");
}

#[test]
fn thin_film_and_bulk_share_the_reststrahlen_band() {
    // Lossless so that eps < 0 gives exactly R = 1 for any thickness: the film
    // is a perfect mirror across the band, transmitting far below it.
    let c = PolarCrystal::silicon_carbide().lossless();
    let m = c.modes[0];
    let mid = 0.5 * (m.to_cm + m.lo_cm);
    let r_2um = c.film_reflectance(mid, 2.0e-6, 1.5);
    let r_thick = c.film_reflectance(mid, 20.0e-6, 1.5);
    let r_out = c.film_reflectance(0.5 * m.to_cm, 2.0e-6, 1.5);
    println!(
        "\n=== SiC film at {mid:.0} cm^-1 (lossless) ===\n  R(2 um)={r_2um:.4}, \
         R(20 um)={r_thick:.4}, R(0.5 TO)={r_out:.4}"
    );
    // A 2 um film is thinner than the evanescent decay length, so a little
    // light tunnels to the substrate; a 20 um film is opaque. Far outside the
    // band the film is a partial mirror at most.
    assert!((0.98..1.0).contains(&r_2um));
    close(r_thick, 1.0, 1e-12);
    assert!(r_out < 0.7, "below the band R should not be near 1: {r_out}");
}

#[test]
fn bulk_polariton_two_branches_and_a_gap() {
    println!("\n=== Bulk polariton dispersion (lossless SiC) ===");
    let c = PolarCrystal::silicon_carbide().lossless();
    let m = c.modes[0];
    let (wto, wlo) = (cm_to_rad(m.to_cm), cm_to_rad(m.lo_cm));

    // Lower branch: real, positive, rises to the TO asymptote.
    let k_lo = c.bulk_polariton_wavevector(0.5 * wto);
    let k_to = c.bulk_polariton_wavevector(wto * (1.0 - 1e-6));
    assert!(k_lo.im.abs() < 1e-6 && k_lo.re > 0.0);
    assert!(k_to.im.abs() < 1e-3 && k_to.re > k_lo.re);
    // Upper branch: k -> 0 as omega -> LO+ and grows with omega.
    let k_lo_plus = c.bulk_polariton_wavevector(wlo * (1.0 + 1e-6));
    let k_high = c.bulk_polariton_wavevector(3.0 * wlo);
    assert!(k_lo_plus.im.abs() < 1e-3 && k_lo_plus.re < k_high.re);
    // Gap: no propagating mode between TO and LO; k is purely imaginary.
    for &frac in &[0.2f64, 0.5, 0.8] {
        let k = c.bulk_polariton_wavevector(wto + frac * (wlo - wto));
        assert!(k.re.abs() < 1e-9, "gap should be evanescent, k={k:?}");
        assert!(k.im > 0.0);
    }
    println!(
        "  omega(TO)/2pi c = {:.1} cm^-1: k(a.u.) lower {:.4}, upper {:.4}; \
         gap evanescent",
        rad_to_cm(wto),
        k_lo.re,
        k_lo_plus.re
    );
}

#[test]
fn surface_phonon_polariton_exists_where_eps_below_minus_one() {
    println!("\n=== Surface phonon polariton (lossless) ===");
    for c in crystals() {
        let lossless = c.lossless();
        let m = lossless.modes[0];
        let ws = lossless.surface_mode_wavenumber();
        // eps = -1 exactly at the surface-mode frequency.
        close(lossless.epsilon_wavenumber(ws).re, -1.0, 1e-9);
        close(ws, lossless.surface_mode_wavenumber_numeric(), 1e-4);
        // eps < -1 on (TO, ws), the range where a surface mode exists.
        assert!(lossless.epsilon_wavenumber(0.5 * (m.to_cm + ws)).re < -1.0);
        let k = lossless.surface_polariton_wavevector(cm_to_rad(ws * (1.0 - 1e-4)));
        assert!(k.im.abs() < 1e-6 && k.re > 0.0);
        // k_sp -> infinity as omega -> ws from below
        let k_lo = lossless.surface_polariton_wavevector(cm_to_rad(ws * (1.0 - 1e-3)));
        let k_hi = lossless.surface_polariton_wavevector(cm_to_rad(ws * (1.0 - 1e-2)));
        assert!(k_lo.re > k_hi.re);
        // No surface mode on (ws, LO): -1 < eps < 0, so eps/(eps+1) < 0.
        let k_bad = lossless.surface_polariton_wavevector(cm_to_rad(0.5 * (ws + m.lo_cm)));
        assert!(k_bad.im.abs() > 1e-3 && k_bad.re.abs() < 1e-9);
        println!(
            "  {:14} omega_s = {ws:.1} cm^-1 (between TO={:.1} and LO={:.1})",
            c.name, m.to_cm, m.lo_cm
        );
    }
}

#[test]
fn otto_atr_dip_angle_matches_analytic_surface_mode() {
    println!("\n=== Otto ATR vs analytic surface-polariton angle (SiC) ===");
    // With a lossless crystal the reflectance is exactly 1 at every angle
    // (total reflection, no absorption channel), so check that first, then use
    // the damped crystal where the surface mode produces a real dip.
    let lossless = PolarCrystal::silicon_carbide().lossless();
    let prism = 2.4;
    let gap = 4.0e-6;
    let nu = 900.0;
    for &theta in &[20.0f64, 27.0, 30.0, 35.0] {
        close(
            lossless.otto_reflectance(prism, gap, nu, theta.to_radians()),
            1.0,
            1e-9,
        );
    }
    let c = PolarCrystal::silicon_carbide();
    let sine = c.surface_polariton_sine(prism, nu);
    let analytic = sine.asin().to_degrees();
    let (theta, r_min) = c.otto_dip_angle(
        prism,
        gap,
        nu,
        (analytic - 6.0).to_radians(),
        (analytic + 6.0).to_radians(),
    );
    let err = theta.to_degrees() - analytic;
    println!(
        "  nu={nu} cm^-1: analytic {analytic:.4} deg, TMM dip {:.4} deg, \
         err {err:+.4} deg, R_min={r_min:.4}",
        theta.to_degrees()
    );
    assert!(err.abs() < 1.0, "angle error {err}");
    assert!(r_min < 0.9, "dip not resolved, R_min={r_min}");
}

#[test]
fn otto_dip_wavenumber_equals_eps_minus_one_solution() {
    // At a fixed prism angle the reflectance dips where k_sp equals the
    // evanescent tangential wavevector, i.e. Re eps = s^2/(1-s^2). Compare the
    // transfer-matrix dip with that solution for the damped crystal.
    println!("\n=== Otto ATR dip wavenumber (damped SiC) ===");
    let prism = 2.4;
    let gap = 6.0e-6;
    let theta = 30.0_f64.to_radians();
    let s = prism * theta.sin();
    assert!(s > 1.0);
    let c = PolarCrystal::silicon_carbide();
    let analytic = c.surface_mode_wavenumber_at_s(s);
    // Scan a window around the analytic surface mode; the reflectance has
    // narrow gap fringes, so sample finely.
    let (lo, hi) = (analytic - 45.0, analytic + 25.0);
    let n = 8000;
    let mut best = f64::INFINITY;
    let mut best_nu = lo;
    for i in 0..=n {
        let nu = lo + (hi - lo) * i as f64 / n as f64;
        let r = c.otto_reflectance(prism, gap, nu, theta);
        if r < best {
            best = r;
            best_nu = nu;
        }
    }
    let err = best_nu - analytic;
    println!(
        "  analytic {analytic:.2} cm^-1, TMM dip {best_nu:.2} cm^-1, \
         err {err:+.2} cm^-1, R_min={best:.4}"
    );
    assert!(err.abs() < 5.0, "dip off by {err}");
}

#[test]
fn diatomic_chain_omega_to_is_the_dielectric_pole() {
    // NaCl-like masses; spring chosen so the zone-centre optical mode is the
    // measured 164 cm^-1 TO frequency of NaCl.
    let m_na = 22.99 * AMU;
    let m_cl = 35.45 * AMU;
    let target_wto = cm_to_rad(164.0);
    let reduced = m_na * m_cl / (m_na + m_cl);
    let spring = 0.5 * target_wto * target_wto * reduced;
    let chain = DiatomicChain::new(5.64e-10, m_na, m_cl, spring);

    // Dynamical matrix vs its closed form.
    let (ac, op) = chain.omega_sq(0.0);
    // Acoustic sum rule: omega_acoustic(0) = 0 to round-off in the squared
    // frequency (the dynamical matrix cancels two large numbers).
    assert!(ac < 1e-10 * op, "acoustic omega^2 at Gamma: {ac:e}");
    close(op, 2.0 * spring * (1.0 / m_na + 1.0 / m_cl), 1e-3);
    // Product of the two squared frequencies is 4 k^2 sin^2(qa/2)/(m1 m2).
    for &q in &[1.0e8f64, 1.0e9, PI / 5.64e-10] {
        let (a2, o2) = chain.omega_sq(q);
        let want = 4.0 * spring * spring * (q * chain.a / 2.0).sin().powi(2) / (m_na * m_cl);
        close(a2 * o2, want, 1e-9 * want);
    }
    // At q = 0 the product vanishes; the cancellation is only round-off.
    let (a2_0, o2_0) = chain.omega_sq(0.0);
    close(a2_0 * o2_0, 0.0, 1e-12 * o2_0 * o2_0);
    close(chain.omega_to(), target_wto, 1e-6 * target_wto);

    // Fold: the Born charge creates the LO mode and the static permittivity.
    let n_pairs = 4.0 / (5.64e-10_f64).powi(3); // 4 NaCl pairs per cubic cell
    let crystal = polar_crystal_from_chain(&chain, 2.34, 1.0, n_pairs, 0.0);
    // The dielectric pole is exactly the chain's zone-centre TO frequency.
    close(cm_to_rad(crystal.to_cm()), chain.omega_to(), 1e-6);
    // LST is an identity of the construction.
    let s = born_oscillator_strength(1.0, n_pairs, chain.reduced_mass());
    let wto = chain.omega_to();
    let wlo = (wto * wto + s / crystal.eps_inf).sqrt();
    close(cm_to_rad(crystal.lo_cm()), wlo, 1e-6);
    close(crystal.static_epsilon(), crystal.eps_inf + s / (wto * wto), 1e-6);
    close(crystal.static_epsilon(), crystal.eps_inf * (crystal.lo_cm() / crystal.to_cm()).powi(2), 1e-9);
    println!(
        "\n=== Lattice-dynamics fold (NaCl masses) ===\n  chain omega_TO = {:.3} cm^-1 \
         (target 164); dielectric pole = {:.3} cm^-1; S = {s:.4e} rad^2/s^2\n  \
         Born-charge omega_LO = {:.2} cm^-1 (measured 264); model eps0 = {:.3} \
         (handbook 5.90)",
        chain.to_wavenumber(),
        crystal.to_cm(),
        rad_to_cm(wlo),
        crystal.static_epsilon()
    );
}

#[test]
fn spatial_dispersion_indicator_is_quadratic_then_fails() {
    // The chain's optical frequency departs from its zone-centre value as
    // (q a)^2 at small q; at the zone boundary the local eps(omega) closure is
    // simply wrong.
    let m_na = 22.99 * AMU;
    let m_cl = 35.45 * AMU;
    let chain = DiatomicChain::new(5.64e-10, m_na, m_cl, 10.0);
    let w0 = chain.omega_to();
    let rel = |q: f64| (chain.omega_sq(q).1.sqrt() - w0).abs() / w0;
    let e1 = rel(1.0e-2 * PI / chain.a);
    let e2 = rel(0.5e-2 * PI / chain.a);
    let ratio = e1 / e2;
    let zone = rel(PI / chain.a);
    println!(
        "\n=== Spatial-dispersion indicator ===\n  rel. shift at q=0.01 pi/a: {e1:.3e}, \
         at half: {e2:.3e}, ratio {ratio:.3} (expect 4 for O(q^2))\n  at q=pi/a: {zone:.3e}; \
         indicator (q a) = {:.3}",
        spatial_dispersion_indicator(PI / chain.a, chain.a)
    );
    assert!((ratio - 4.0).abs() < 0.05, "not quadratic: ratio {ratio}");
    assert!(e1 < 5.0e-4);
    assert!(zone > 0.1, "zone-boundary shift should be large");
    // Film indicator: a 10-lattice-spacing film is homogeneous, a 2-spacing
    // film is not.
    assert!(film_thickness_indicator(chain.a, 10.0 * chain.a) < 0.11);
    assert!(film_thickness_indicator(chain.a, 2.0 * chain.a) > 0.4);
}

#[test]
fn dispersion_trait_and_wavenumber_conversions() {
    let c = PolarCrystal::silicon_carbide();
    let nu = 850.0;
    let lambda = wavenumber_to_lambda(nu);
    close(lambda_to_wavenumber(lambda), nu, 1e-6);
    // The Dispersion trait must agree with epsilon_wavenumber.
    close((c.index(lambda) - c.epsilon_wavenumber(nu).sqrt()).abs(), 0.0, 1e-12);

}
