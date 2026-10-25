//! The isobaric fold: a mushy-zone closure and its resolution indicator.
//!
//! Run with `cargo run --release -p solidification --example isobaric_closure`.
//!
//! The macroscopic isobaric model of `analysis/compressible.md` §15 carries
//! one enthalpy per cell. When the dendritic microstructure is finer than the
//! cell, its effect must be a closure: the latent-heat release rate and the
//! mushy-zone fraction as functions of the resolved undercooling. This example
//! evaluates that closure for succinonitrile and ice, gives the dimensional
//! tip velocity `v(Delta T)`, the Avrami solid fraction, and the indicator
//! (tip diameter vs cell size) that says when the macroscopic model must
//! unfold to the phase field.

use solidification::closure::*;
use solidification::materials::{ICE, SUCCINONITRILE};

/// `sigma_star` for the 2D solvability condition `R^2 v = 2 D d0 / sigma_star`.
/// The value is calibrated once to the Karma-Rappel benchmark quoted in the
/// stage brief (`V d0/D ~ 0.017` at `Delta = 0.55`, `eps4 = 0.05`): from
/// `p = 0.2554` this gives `sigma_star = 0.017/(2 p^2)`. It is a *closure
/// parameter*, not a fit to the experimental number below.
fn sigma_star_2d() -> f64 {
    // The 2D Ivantsov Peclet number for Delta = 0.55 is 0.25693, so
    // sigma_star = 0.0170/(2*0.25693^2) = 0.1288.
    sigma_star_2d_kr_benchmark()
}

fn main() {
    let sigma = sigma_star_2d();
    println!("2D solvability closure: sigma_star = {sigma:.4}\n");
    for m in [SUCCINONITRILE, ICE] {
        println!(
            "{}: d0 = {:.3e} m, D = {:.3e} m^2/s",
            m.name, m.d0, m.alpha_liquid
        );
        println!("  dT(K)   Delta    p        v (m/s)      R (m)      tip diam (m)");
        for &dt_k in &[0.5, 1.0, 2.0, 5.0] {
            let delta = dt_k / (m.latent / m.c_vol_solid);
            if delta >= 1.0 {
                continue;
            }
            let gp = growth(delta, m.d0, m.alpha_liquid, sigma);
            println!(
                "  {:>5.1}  {:>6.4}  {:>6.4}  {:>10.3e}  {:>10.3e}  {:>10.3e}",
                dt_k,
                delta,
                gp.pe,
                gp.v,
                gp.r,
                gp.tip_diameter()
            );
        }
        // Indicator for a macroscopic cell of 1 mm (a lake-ice or casting cell).
        let dx_macro = 1.0e-3;
        let gp = growth(
            2.0 / (m.latent / m.c_vol_solid),
            m.d0,
            m.alpha_liquid,
            sigma,
        );
        let res = microstructure_indicator(&gp, dx_macro, 6.0);
        println!(
            "  at dT=2 K, tip diameter/cell = {:.3e} -> {:?}\n",
            gp.tip_diameter() / dx_macro,
            res
        );
    }

    // Avrami mushy-zone fraction for a fixed tip velocity, appropriate when the
    // undercooling is nearly uniform. N = 1e12 m^-3 (a modest nucleation
    // density), v = 50 um/s.
    let v = 50e-6;
    let n0 = 1e12;
    println!("Avrami (3D): N={n0:.1e} m^-3, v={v:.1e} m/s");
    for &t in &[0.0, 0.5, 1.0, 2.0, 4.0] {
        let x = avrami_from_growth(t, n0, v, 3);
        let dx = avrami_rate_from_growth(t, n0, v, 3);
        println!(
            "  t={t:>4.1} s: X={x:.4}  dX/dt={dx:.3e} /s  L dX/dt={:.3e} W/m^3",
            SUCCINONITRILE.latent * dx
        );
    }
}
