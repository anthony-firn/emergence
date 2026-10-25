//! The device law's temperature dependence against the drift–diffusion
//! transistors it was folded from.
//!
//! The TCAD currents below were produced by
//! `cargo run -p semiconductor --release --example ekv_temperature`
//! (150 nm n- and pMOS, full mobility model, |V_ds| = 1.2 V); that run takes
//! ~15 min, so its numbers are recorded here and the fitted coefficients in
//! `gates::{N,P}_TEMPERATURE` are checked against them, including the four
//! temperatures held out of the fit.

use semiconductor::circuit::Fet;
use semiconductor::gates::Ekv;

const T: [f64; 8] = [250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0];
const HELD: [bool; 8] = [false, true, false, true, false, true, false, true];
// (I_off, I_on) in A/µm
const N_TCAD: [(f64, f64); 8] = [
    (1.4625e-11, 3.87934e-4),
    (7.6921e-11, 3.67902e-4),
    (3.1214e-10, 3.49479e-4),
    (1.0373e-9, 3.32605e-4),
    (2.9460e-9, 3.17196e-4),
    (7.3777e-9, 3.03175e-4),
    (1.6677e-8, 2.90466e-4),
    (3.4630e-8, 2.78997e-4),
];

fn check(dev: Ekv, tcad: &[(f64, f64)], on_tol_held: f64, off_tol_held: f64) {
    let pol = dev.pol;
    let i = |d: &Ekv, vgs: f64| (d.id(pol * vgs, pol * 1.2)).abs();
    let base = (i(&dev, 0.0), i(&dev, 1.2));
    let (b_off, b_on) = tcad[2];
    let mut worst = (0.0f64, 0.0f64);
    for k in 0..T.len() {
        let d = dev.at(T[k], 1.2).unwrap();
        let on = (i(&d, 1.2) / base.1) / (tcad[k].1 / b_on);
        let off = (i(&d, 0.0) / base.0) / (tcad[k].0 / b_off);
        println!("{} T {:.0} K{}: I_on ratio {on:.4}, I_off ratio {off:.4}", if pol > 0.0 { "n" } else { "p" }, T[k], if HELD[k] { " (held out)" } else { "" });
        if HELD[k] {
            worst = (worst.0.max((on - 1.0).abs()), worst.1.max((off - 1.0).abs()));
        }
    }
    assert!(worst.0 < on_tol_held && worst.1 < off_tol_held, "held-out errors {worst:?}");
}

#[test]
fn nmos_temperature_fold() {
    check(Ekv::n(), &N_TCAD, 0.035, 0.10);
}

#[test]
fn pmos_temperature_fold() {
    const P_TCAD: [(f64, f64); 8] = [
        (5.8399e-12, 1.71688e-4),
        (3.1196e-11, 1.61298e-4),
        (1.2886e-10, 1.52376e-4),
        (4.3638e-10, 1.44666e-4),
        (1.2636e-9, 1.37962e-4),
        (3.2264e-9, 1.32115e-4),
        (7.4328e-9, 1.27009e-4),
        (1.5721e-8, 1.22550e-4),
    ];
    check(Ekv::p(), &P_TCAD, 0.035, 0.12);
}

#[test]
fn outside_the_calibrated_range_is_refused() {
    assert!(Ekv::n().at(240.0, 1.2).is_err());
    assert!(Ekv::n().at(430.0, 1.2).is_err());
}
