//! Folding a sub-wavelength layered medium into an effective homogeneous one,
//! and the classical bounds for isotropic two-phase composites.
//!
//! This is the *coarse* rung above the exact transfer matrix. At normal
//! incidence a laminate loaded in series is an isostress problem, so its exact
//! long-wavelength modulus is the Reuss (harmonic) average
//!
//! ```text
//!   1 / M_eff = sum_i f_i / M_i ,      rho_eff = sum_i f_i rho_i ,
//! ```
//!
//! with `f_i = d_i / Lambda` the volume fractions. The exact Bloch dispersion
//! of the periodic stack ([`crate::layered::Periodic`]) must converge to
//! `c_eff = sqrt(M_eff / rho_eff)` as `q Lambda -> 0`; that is the benchmark.
//!
//! For an isotropic two-phase composite the effective bulk and shear moduli
//! are bracketed by the Voigt (arithmetic) and Reuss (harmonic) averages and,
//! more tightly, by the Hashin–Shtrikman variational bounds. The ordering
//! `Reuss <= HS_lower <= HS_upper <= Voigt` is exact for a two-phase mixture
//! and is tested.
//!
//! # Validity indicator
//!
//! The folded medium is valid when the phase cannot resolve one period. The
//! relative phase-velocity error is `O((omega Lambda / c_eff)^2)`, in step with
//! the laminate expansion of the Rytov relation; [`homogenisation_error`] is
//! the measured quadratic envelope. Unfold to [`crate::layered`] when it
//! exceeds the required accuracy.
use crate::complex::C;
use crate::layered::Periodic;

/// Reuss (harmonic / series) average of positive moduli.
pub fn reuss(moduli: &[f64], fractions: &[f64]) -> f64 {
    assert_eq!(moduli.len(), fractions.len());
    let inv: f64 = moduli
        .iter()
        .zip(fractions)
        .map(|(&m, &f)| f / m)
        .sum();
    1.0 / inv
}

/// Voigt (arithmetic / parallel) average.
pub fn voigt(moduli: &[f64], fractions: &[f64]) -> f64 {
    assert_eq!(moduli.len(), fractions.len());
    moduli.iter().zip(fractions).map(|(&m, &f)| f * m).sum()
}

/// Exact long-wavelength normal-incidence modulus of a laminate.
pub fn laminate_effective_modulus(layers: &[(f64, f64)]) -> f64 {
    let fractions: Vec<f64> = layers.iter().map(|&(_, d)| d).collect();
    let total: f64 = fractions.iter().sum();
    let f: Vec<f64> = fractions.iter().map(|d| d / total).collect();
    let m: Vec<f64> = layers.iter().map(|&(md, _)| md).collect();
    reuss(&m, &f)
}

/// Effective density of a laminate, `sum f_i rho_i`.
pub fn laminate_effective_density(layers: &[(f64, f64)]) -> f64 {
    let total: f64 = layers.iter().map(|&(_, d)| d).sum();
    layers.iter().map(|&(r, d)| r * d / total).sum()
}

/// Effective normal-incidence phase speed of a laminate, from parallel
/// modulus, density and thickness slices.
pub fn laminate_effective_speed(moduli: &[f64], densities: &[f64], thicknesses: &[f64]) -> f64 {
    assert_eq!(moduli.len(), thicknesses.len());
    assert_eq!(densities.len(), thicknesses.len());
    let m: Vec<(f64, f64)> = moduli.iter().copied().zip(thicknesses.iter().copied()).collect();
    let r: Vec<(f64, f64)> = densities.iter().copied().zip(thicknesses.iter().copied()).collect();
    (laminate_effective_modulus(&m) / laminate_effective_density(&r)).sqrt()
}

/// Hashin–Shtrikman bounds `(lower, upper)` on the effective bulk modulus of
/// an isotropic two-phase composite. Phase `1` and `2` may be given in either
/// order; the stiffer phase is identified internally.
pub fn hashin_shtrikman_bulk(k1: f64, g1: f64, k2: f64, g2: f64, f1: f64) -> (f64, f64) {
    assert!((0.0..=1.0).contains(&f1));
    let (ks, gs, f_s, kw, gw, f_w) = if k1 >= k2 {
        (k1, g1, f1, k2, g2, 1.0 - f1)
    } else {
        (k2, g2, 1.0 - f1, k1, g1, f1)
    };
    let lower = kw + f_s / (1.0 / (ks - kw) + 3.0 * f_w / (3.0 * kw + 4.0 * gw));
    let upper = ks + f_w / (1.0 / (kw - ks) + 3.0 * f_s / (3.0 * ks + 4.0 * gs));
    (lower, upper)
}

/// Hashin–Shtrikman bounds `(lower, upper)` on the effective shear modulus.
pub fn hashin_shtrikman_shear(k1: f64, g1: f64, k2: f64, g2: f64, f1: f64) -> (f64, f64) {
    assert!((0.0..=1.0).contains(&f1));
    let (ks, gs, f_s, kw, gw, f_w) = if g1 >= g2 {
        (k1, g1, f1, k2, g2, 1.0 - f1)
    } else {
        (k2, g2, 1.0 - f1, k1, g1, f1)
    };
    let lower = gw
        + f_s
            / (1.0 / (gs - gw) + 6.0 * f_w * (kw + 2.0 * gw) / (5.0 * gw * (3.0 * kw + 4.0 * gw)));
    let upper = gs
        + f_w
            / (1.0 / (gw - gs) + 6.0 * f_s * (ks + 2.0 * gs) / (5.0 * gs * (3.0 * ks + 4.0 * gs)));
    (lower, upper)
}

/// Exact Bloch phase velocity of the infinite periodic medium at `omega`,
/// `v = omega / Re q` (folded to the first Brillouin zone).
pub fn bloch_phase_velocity(periodic: &Periodic, omega: f64) -> f64 {
    let q = periodic.bloch_wavevector(omega).re;
    omega / q
}

/// Relative phase-velocity error of the folded model against the exact Bloch
/// dispersion of the periodic medium. Measured, not an estimate: this is the
/// quantity the indicator is meant to bound.
pub fn phase_velocity_error(periodic: &Periodic, omega: f64, c_eff: f64) -> f64 {
    let v = bloch_phase_velocity(periodic, omega);
    (v - c_eff).abs() / c_eff
}

/// Empirical quadratic validity indicator for the folded laminate,
/// `C (omega Lambda / c_eff)^2`.
///
/// The coefficient was measured against the exact Rytov dispersion of several
/// two-layer stacks (see `analysis/phononic.md`); the leading relative
/// phase-velocity error is `-(1/24)(omega Lambda/c_eff)^2` times a contrast
/// factor of order one, so this envelope is a switch rather than a proven
/// bound. It can be exceeded near a constituent resonance or at very high
/// contrast.
pub fn homogenisation_error(omega: f64, period: f64, c_eff: f64) -> f64 {
    let x = omega * period / c_eff;
    x * x / 12.0
}

/// Period-to-wavelength indicator `Lambda / lambda`.
pub fn period_over_wavelength(period: f64, wavelength: f64) -> f64 {
    period / wavelength
}

/// Wavelength of a lossless layer at `omega`.
pub fn wavelength(layer: &crate::layered::Layer, omega: f64) -> f64 {
    2.0 * std::f64::consts::PI / layer.wavenumber(omega).re
}

/// Real part helper for a lossless medium.
pub fn real_part(c: C) -> f64 {
    assert!(c.im.abs() < 1e-12, "expected a real quantity, got {c:?}");
    c.re
}
