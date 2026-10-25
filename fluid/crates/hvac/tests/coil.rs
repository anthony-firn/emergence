//! Coils and exchangers: balances and exact limits.

use hvac::coil::*;
use hvac::psychro::{dew_point, enthalpy, w_from_rh, w_s, P_ATM};

#[test]
fn cooling_coil_balances() {
    // 26.7 °C / 50 % room air (ASHRAE's common rating point) through a coil
    // with a 10 °C apparatus dew point and bypass factor 0.1.
    let p = P_ATM;
    let w_in = w_from_rh(26.7, 0.5, p);
    let m = 1.2;
    let c = cooling_coil(26.7, w_in, 10.0, 0.1, m, p);
    println!(
        "coil: leaves {:.2} °C, {:.5} kg/kg; total {:.0} W, sensible {:.0} W, latent {:.0} W (SHR {:.2}), condensate {:.2} g/s",
        c.t_out, c.w_out, c.total, c.sensible, c.latent, c.sensible / c.total, c.condensate * 1e3
    );
    assert!(c.wet && c.latent > 0.0);
    // Energy: the total is the enthalpy drop; the latent part carries the
    // condensed water's latent heat (2501 kJ/kg at the reference, plus the
    // vapour's sensible term in ASHRAE's enthalpy).
    assert!((c.total - m * (enthalpy(26.7, w_in) - enthalpy(c.t_out, c.w_out))).abs() < 1e-9);
    let latent = m * (w_in - c.w_out) * (2501e3 + 1860.0 * c.t_out);
    assert!((c.latent / latent - 1.0).abs() < 1e-9, "{} vs {latent}", c.latent);
    // The leaving state lies on the line to the ADP: no bypass gives
    // saturated air at the ADP.
    let full = cooling_coil(26.7, w_in, 10.0, 0.0, m, p);
    assert!((full.t_out - 10.0).abs() < 1e-12 && (full.w_out - w_s(10.0, p)).abs() < 1e-15);
    // An ADP above the entering dew point leaves a dry coil.
    let dry = cooling_coil(26.7, w_in, dew_point(w_in, p) + 2.0, 0.1, m, p);
    assert!(!dry.wet && dry.latent.abs() < 1e-9 && dry.w_out == w_in);
}

#[test]
fn effectiveness_limits() {
    // C_r → 0 (condensing or evaporating stream): every arrangement gives
    // 1 − e^−NTU.
    for arr in [Arrangement::Counterflow, Arrangement::Parallel, Arrangement::CrossflowUnmixed] {
        assert!((effectiveness(arr, 1.3, 0.0) - (1.0 - (-1.3f64).exp())).abs() < 1e-15);
    }
    // Balanced counterflow: NTU/(1 + NTU); parallel never exceeds 1/(1 + C_r).
    assert!((effectiveness(Arrangement::Counterflow, 3.0, 1.0) - 0.75).abs() < 1e-15);
    assert!((effectiveness(Arrangement::Counterflow, 3.0, 1.0 - 1e-9) - 0.75).abs() < 1e-8);
    assert!(effectiveness(Arrangement::Parallel, 50.0, 0.5) < 1.0 / 1.5 + 1e-12);
    // Crossflow lies between parallel and counterflow.
    for (ntu, cr) in [(0.5, 0.5), (2.0, 0.8), (4.0, 0.3)] {
        let (p, x, c) = (
            effectiveness(Arrangement::Parallel, ntu, cr),
            effectiveness(Arrangement::CrossflowUnmixed, ntu, cr),
            effectiveness(Arrangement::Counterflow, ntu, cr),
        );
        assert!(p < x && x < c, "{ntu} {cr}: {p} {x} {c}");
    }
}
