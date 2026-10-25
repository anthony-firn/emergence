//! Surface plasmon resonance in the Kretschmann configuration (optional rung).
//!
//! A thin metal film on the flat face of a glass prism couples the evanescent
//! prism field to a surface plasmon polariton at the outer metal/air interface.
//! The transfer matrix gives the reflectance dip; the analytic condition for
//! the plasmon wavevector is
//!
//! ```text
//!   k_x = k0 sqrt(eps_m eps_d / (eps_m + eps_d)),
//!   sin(theta_SP) = sqrt(eps_m eps_d / (eps_m + eps_d)) / n_prism.
//! ```
//!
//! With a real metal permittivity the dip angle equals the analytic value.
//! With the causal Drude-Lorentz gold permittivity the dip is shifted and
//! broadened by the imaginary part and by the finite film thickness; the test
//! reports that shift rather than tuning it away.
use crate::complex::C;
use crate::dispersion::Dispersion;
use crate::multilayer::{Layer, Polarization, Stack};

/// Analytic surface-plasmon sine, using the complex metal permittivity.
pub fn surface_plasmon_sine(prism_index: f64, eps_metal: C, eps_dielectric: f64) -> f64 {
    let ratio = eps_metal * eps_dielectric / (eps_metal + C::real(eps_dielectric));
    ratio.sqrt().re / prism_index
}

/// p-polarised reflectance of prism / metal film / ambient at angle `theta`.
pub fn kretschmann_reflectance(
    prism_index: f64,
    metal: C,
    thickness: f64,
    ambient_index: f64,
    lambda: f64,
    theta: f64,
) -> f64 {
    let stack = Stack::with_complex(
        C::real(prism_index),
        vec![Layer {
            index: metal,
            mu: C::real(1.0),
            thickness,
        }],
        C::real(ambient_index),
    );
    stack.solve(lambda, theta, Polarization::TM).reflectance
}

/// Reflectance of a metal film using any dispersion model.
pub fn kretschmann_reflectance_model(
    prism_index: f64,
    metal: &dyn Dispersion,
    thickness: f64,
    ambient_index: f64,
    lambda: f64,
    theta: f64,
) -> f64 {
    kretschmann_reflectance(
        prism_index,
        metal.index(lambda),
        thickness,
        ambient_index,
        lambda,
        theta,
    )
}

/// Search the reflectance minimum over `[theta_lo, theta_hi]` by ternary
/// search (the dip is unimodal on the tested interval). Returns `(angle, R)`.
pub fn reflectance_dip(
    prism_index: f64,
    metal: C,
    thickness: f64,
    ambient_index: f64,
    lambda: f64,
    theta_lo: f64,
    theta_hi: f64,
) -> (f64, f64) {
    let f =
        |t: f64| kretschmann_reflectance(prism_index, metal, thickness, ambient_index, lambda, t);
    let (mut a, mut b) = (theta_lo, theta_hi);
    for _ in 0..200 {
        let m1 = a + (b - a) / 3.0;
        let m2 = b - (b - a) / 3.0;
        if f(m1) < f(m2) {
            b = m2;
        } else {
            a = m1;
        }
    }
    let t = 0.5 * (a + b);
    (t, f(t))
}
