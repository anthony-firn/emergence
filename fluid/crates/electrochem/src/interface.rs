//! The electrochemical interface rung: Nernst, Butler–Volmer, and the Tafel
//! and linear limits it folds to.
//!
//! For the one-electron-step reaction `O + n e⁻ ⇌ R` with activities
//! `a_O`, `a_R`,
//!
//! ```text
//!   Nernst:          E = E° + (RT/nF) ln(a_O/a_R)
//!   Butler–Volmer:   i = i₀ [ exp(α_a nFη/RT) − exp(−α_c nFη/RT) ]
//!   Tafel (η≫RT/nF): i ≈ i₀ exp(α_a nFη/RT), slope dη/dlog₁₀i = 2.303RT/(α_a nF)
//!   linear (η≪RT/nF): i ≈ i₀ nFη/RT
//! ```
//!
//! References: Bard & Faulkner, *Electrochemical Methods*, 2nd ed., ch. 1 and
//! 3 (Nernst, Butler–Volmer, Tafel, linear polarisation); the diffusion-limited
//! current uses the Nernst diffusion-layer thickness δ.

use crate::{F_FARADAY, R_GAS};

#[cfg(test)]
use crate::T_REF;

/// Nernst potential of `O + n e⁻ ⇌ R` (V) from the standard potential `e0`
/// and the activities (or concentrations, same units) of the oxidised and
/// reduced species.
pub fn nernst_potential(e0: f64, n: f64, t: f64, a_ox: f64, a_red: f64) -> f64 {
    e0 + (R_GAS * t / (n * F_FARADAY)) * (a_ox / a_red).ln()
}

/// The `RT/F` prefactor in volts, exposed for convenience.
pub fn rt_over_f(t: f64) -> f64 {
    R_GAS * t / F_FARADAY
}

/// Butler–Volmer current density (A m⁻²) for overpotential `eta` (V).
///
/// * `i0` — exchange current density (A m⁻²)
/// * `alpha_a`, `alpha_c` — anodic/cathodic transfer coefficients
/// * `n` — electrons per rate-determining step
pub fn butler_volmer(eta: f64, i0: f64, alpha_a: f64, alpha_c: f64, n: f64, t: f64) -> f64 {
    let f = n * F_FARADAY / (R_GAS * t);
    i0 * (clamped_exp(alpha_a * f * eta) - clamped_exp(-alpha_c * f * eta))
}

/// Invert Butler–Volmer: the overpotential (V) that carries current density
/// `i` (A m⁻²). Newton iteration on `butler_volmer`, warm-started from the
/// symmetric `α=½` solution `η = 2 asinh(i/2i₀)/f`.
///
/// Returns `None` for `i0 <= 0`.
pub fn butler_volmer_overpotential(
    i: f64,
    i0: f64,
    alpha_a: f64,
    alpha_c: f64,
    n: f64,
    t: f64,
) -> Option<f64> {
    if !(i0 > 0.0) {
        return None;
    }
    let f = n * F_FARADAY / (R_GAS * t);
    // Warm start: for α_a = α_c = ½ the inverse is exactly 2 asinh(i/2i0)/f.
    let mut eta = 2.0 * (i / (2.0 * i0)).asinh() / f;
    for _ in 0..100 {
        let ea = clamped_exp(alpha_a * f * eta);
        let ec = clamped_exp(-alpha_c * f * eta);
        let res = i0 * (ea - ec) - i;
        let dres = i0 * f * (alpha_a * ea + alpha_c * ec);
        if dres <= 0.0 || !dres.is_finite() {
            break;
        }
        let step = res / dres;
        eta -= step;
        if step.abs() < 1e-14 * (1.0 + eta.abs()) {
            break;
        }
    }
    Some(eta)
}

/// Tafel slope `dη/dlog₁₀i = 2.303 RT/(α n F)` (V per decade).
pub fn tafel_slope(alpha: f64, n: f64, t: f64) -> f64 {
    (10.0f64).ln() * R_GAS * t / (alpha * n * F_FARADAY)
}

/// Asymptotic Tafel current density (A m⁻²) for large anodic overpotential.
pub fn tafel_current(eta: f64, i0: f64, alpha_a: f64, n: f64, t: f64) -> f64 {
    let f = n * F_FARADAY / (R_GAS * t);
    i0 * clamped_exp(alpha_a * f * eta)
}

/// Linear (small-signal) current density (A m⁻²): `i = i₀ nFη/RT`.
pub fn linear_current(eta: f64, i0: f64, n: f64, t: f64) -> f64 {
    i0 * n * F_FARADAY * eta / (R_GAS * t)
}

/// Overpotential for a target Tafel current, `η = (RT/α_a nF) ln(i/i₀)`.
pub fn tafel_overpotential(i: f64, i0: f64, alpha_a: f64, n: f64, t: f64) -> f64 {
    let f = n * F_FARADAY / (R_GAS * t);
    (i / i0).ln() / (alpha_a * f)
}

/// Diffusion-limited (plateau) current density magnitude (A m⁻²) through a
/// Nernst diffusion layer of thickness `delta` (m):
/// `i_L = n F D c_bulk / δ`.
pub fn limiting_current(n: f64, d: f64, c_bulk: f64, delta: f64) -> f64 {
    n * F_FARADAY * d * c_bulk / delta
}

/// Surface concentration of the reacting species (mol m⁻³) when the total
/// current is `i`, given bulk `c_bulk` and plateau `i_lim`:
/// `c_s = c_bulk (1 − i/i_L)`.
pub fn surface_concentration(c_bulk: f64, i: f64, i_lim: f64) -> f64 {
    c_bulk * (1.0 - i / i_lim)
}

/// Current when a kinetic current `i_k` (Butler–Volmer, no transport limit)
/// is in series with a mass-transport limit `i_L`: `i = i_k i_L/(i_k + i_L)`.
///
/// This is the mass-transport-corrected (Koutecký–Levich) current; it recovers
/// the kinetic current for `|i_k| ≪ i_L` and the plateau for `|i_k| ≫ i_L`.
pub fn mixed_current(i_k: f64, i_lim: f64) -> f64 {
    i_k * i_lim / (i_k + i_lim)
}

#[inline]
fn clamped_exp(x: f64) -> f64 {
    x.clamp(-700.0, 700.0).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    const T: f64 = T_REF;

    #[test]
    fn nernst_equal_activities_is_standard_potential() {
        assert_eq!(nernst_potential(-0.5, 1.0, T, 1.0, 1.0), -0.5);
    }

    #[test]
    fn nernst_decade_is_59_mv_per_n1() {
        // One decade of a_O/a_R at 298.15 K: RT/F ln(10) = 59.16 mV.
        let e = nernst_potential(0.0, 1.0, T, 10.0, 1.0);
        let d = 10.0f64.ln() * rt_over_f(T);
        assert!((e - d).abs() < 1e-15);
        assert!((d * 1e3 - 59.16).abs() < 0.01);
    }

    #[test]
    fn tafel_slope_matches_closed_form() {
        // Numeric slope over two decades of Tafel current vs the closed form.
        let (i0, aa, _ac, n) = (1.0, 0.5, 0.5, 1.0);
        let i1 = 1.0;
        let i2 = 100.0;
        let e1 = tafel_overpotential(i1, i0, aa, n, T);
        let e2 = tafel_overpotential(i2, i0, aa, n, T);
        let numeric = (e2 - e1) / (i2.log10() - i1.log10());
        let exact = tafel_slope(aa, n, T);
        assert!((numeric - exact).abs() < 1e-12);
        // 2.303 RT/(αnF) for α=0.5, n=1 at 298.15 K ≈ 118.3 mV/decade.
        assert!((exact * 1e3 - 118.3).abs() < 0.1);
    }

    #[test]
    fn butler_volmer_reduces_to_linear_and_tafel() {
        let (i0, aa, ac, n) = (1e-3, 0.5, 0.5, 1.0);
        // Small overpotential → linear.
        let eta_small = 1e-6;
        let i_bv = butler_volmer(eta_small, i0, aa, ac, n, T);
        let i_lin = linear_current(eta_small, i0, n, T);
        assert!((i_bv - i_lin).abs() / i_lin.abs() < 1e-9);
        // Large anodic overpotential → Tafel.
        let eta_big = 0.3;
        let i_bv = butler_volmer(eta_big, i0, aa, ac, n, T);
        let i_taf = tafel_current(eta_big, i0, aa, n, T);
        assert!((i_bv - i_taf).abs() / i_taf.abs() < 1e-5);
    }

    #[test]
    fn bv_inverse_round_trips() {
        let (i0, aa, ac, n) = (1e-3, 0.5, 0.5, 1.0);
        for &i in &[-50.0, -1.0, -1e-6, 0.0, 1e-6, 1.0, 50.0] {
            let eta = butler_volmer_overpotential(i, i0, aa, ac, n, T).unwrap();
            let back = butler_volmer(eta, i0, aa, ac, n, T);
            assert!((back - i).abs() < 1e-9 * (1.0 + i.abs()), "i={i} back={back}");
        }
    }

    #[test]
    fn mass_transport_limits() {
        let i_lim = limiting_current(1.0, 1e-9, 1000.0, 1e-4);
        assert!(i_lim > 0.0);
        // Kinetic current far below the plateau → unchanged.
        assert!((mixed_current(1e-6, i_lim) - 1e-6).abs() / 1e-6 < 1e-3);
        // Far above → plateau.
        assert!((mixed_current(1e9, i_lim) - i_lim).abs() / i_lim < 1e-6);
        // Surface depletion at half the plateau.
        assert!((surface_concentration(1000.0, 0.5 * i_lim, i_lim) - 500.0).abs() < 1e-9);
    }
}
