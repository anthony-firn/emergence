//! The isobaric (macroscopic) fold: a mushy-zone closure and its indicator.
//!
//! # Why a closure
//!
//! A macroscopic heat/enthalpy solver at cell size `dx_macro` cannot resolve a
//! dendrite whose tip radius is `R ~ 10 d0 = 10 nm`–`10 µm`. The solid fraction
//! and the latent-heat release rate are then sub-grid. The closure here gives
//! them from the resolved undercooling `Delta = (Tm - T)/(L/c_p)`:
//!
//! 1. the **dendrite tip kinetics** `v(Delta)` from Ivantsov's diffusion
//!    solution plus the microscopic-solvability selection condition;
//! 2. an **Avrami** growth law converting that velocity and a nucleation
//!    density into a solid fraction and a latent-heat release rate.
//!
//! # Ivantsov and selection
//!
//! For a 2D parabolic needle the thermal field gives (Ivantsov 1947)
//!
//! ```text
//!   Delta = sqrt(pi p) e^p erfc(sqrt(p)),      p = v R / (2 D),
//! ```
//!
//! which fixes only `v R`. The missing relation is the solvability condition
//! (Langer & Müller-Krumbhaar 1978; Karma & Rappel 1998)
//!
//! ```text
//!   R^2 v = 2 D d0 / sigma_star,
//! ```
//!
//! with `sigma_star` the (anisotropy-dependent) selection constant. Combining
//! the two, `v = 2 D sigma_star p^2 / d0`, so `Delta` alone determines `p`.
//! `sigma_star` is a *closure parameter*; it can be calibrated once against the
//! phase-field tip or taken from the 2D solvability literature. It is never
//! tuned inside a validation test.
//!
//! # The selection constant for the Karma-Rappel benchmark
//!
//! At `Delta = 0.55`, `eps4 = 0.05`, the Karma-Rappel 2D phase-field /
//! Green-function benchmark has the sharp-interface tip speed
//! `V d0/D = 0.0170` (Karma & Rappel, *Phys. Rev. E* **57**, 4323 (1998),
//! Table II; the phase-field column is `0.0168`). The same simulations show
//! that the operating point is Ivantsov to within the measurement error: the
//! correctly measured `p = 0` tip radius is `R = 2 D p/V`, with
//! `p = 0.2569` the 2D Ivantsov value for `Delta = 0.55`. The selection
//! constant follows from `V` alone,
//!
//! ```text
//!   sigma_star = (V d0/D) / (2 p^2) = 0.0170 / (2 * 0.2569^2) = 0.1288.
//! ```
//!
//! This is the value returned by [`sigma_star_2d_kr_benchmark`]. (An earlier
//! note reported `sigma_star ~ 0.5-0.7` from the phase field; that came from a
//! tip-radius fit over only the diffuse-interface cap, which underestimates `R`
//! by ~2.4x. Correcting the fit gives `R/W0 = 9.3` at `D = 2` and `5.5` at
//! `D = 3`, matching `2 D p/V`, and `sigma_star = 0.12-0.13`.)

use crate::stefan::erfc;

/// Solve the 2D Ivantsov relation `Delta = sqrt(pi p) e^p erfc(sqrt(p))` for
/// the Peclet number `p`. Monotone in `p`; bisection on `[0, 4]`, capped so
/// `e^p` cannot overflow.
pub fn ivantsov_peclet_2d(delta: f64) -> f64 {
    assert!((0.0..1.0).contains(&delta), "Delta must be in (0,1)");
    let f = |p: f64| {
        let s = p.max(0.0).sqrt();
        std::f64::consts::PI.sqrt() * s * p.exp() * erfc(s) - delta
    };
    let mut lo = 0.0;
    let mut hi = 4.0;
    // f(4) is ~1 - tiny; ensure bracketing for extreme Delta.
    while f(hi) < 0.0 && hi < 40.0 {
        hi *= 1.5;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// A point on the closure curve: undercooling -> tip velocity and radius.
#[derive(Clone, Copy, Debug)]
pub struct GrowthPoint {
    /// Dimensionless undercooling `Delta`.
    pub delta: f64,
    /// Ivantsov Peclet number `p = v R/(2D)`.
    pub pe: f64,
    /// Tip velocity in the same units as `D/d0`.
    pub v: f64,
    /// Tip radius in the same units as `d0`.
    pub r: f64,
}

impl GrowthPoint {
    /// Tip diameter (physical length, same units as the `d0` passed to
    /// [`growth`]).
    pub fn tip_diameter(&self) -> f64 {
        2.0 * self.r
    }
}

/// Tip kinetics from Ivantsov + solvability. `d0` and `d_diff` are the
/// capillary length and thermal diffusivity. The returned velocity is in
/// `d_diff/d0` units (e.g. m/s when those are SI) and the radius is a physical
/// length (e.g. m). `sigma_star` is the closure selection constant.
///
/// For the Karma-Rappel 2D benchmark use [`sigma_star_2d_kr_benchmark`].
pub fn growth(delta: f64, d0: f64, d_diff: f64, sigma_star: f64) -> GrowthPoint {
    let pe = ivantsov_peclet_2d(delta);
    let v = 2.0 * d_diff * sigma_star * pe * pe / d0;
    let r = 2.0 * d_diff * pe / v;
    GrowthPoint { delta, pe, v, r }
}

/// The 2D microscopic-solvability selection constant `sigma_star` for the
/// Karma-Rappel benchmark `Delta = 0.55`, `eps4 = 0.05`.
///
/// It is fixed by the published sharp-interface tip speed `V d0/D = 0.0170`
/// (Karma & Rappel 1998, Table II) and the 2D Ivantsov Peclet number
/// `p(0.55) = 0.25693`:
/// `sigma_star = 2 D d0/(R^2 V) = (V d0/D)/(2 p^2) = 0.1288`.
/// This is a *calibration* of the closure, not a fit to a validation number.
pub fn sigma_star_2d_kr_benchmark() -> f64 {
    const V_KR_GREEN: f64 = 0.0170;
    let p = ivantsov_peclet_2d(0.55);
    V_KR_GREEN / (2.0 * p * p)
}

/// Selection constant `sigma_star = 2 D d0/(R^2 v)` of a measured tip operating
/// point, for checking a simulation against the closure.
pub fn selection_sigma(v: f64, r: f64, d0: f64, d_diff: f64) -> f64 {
    2.0 * d_diff * d0 / (r * r * v)
}

/// Avrami solid fraction `X(t) = 1 - exp(-(t/tau)^n)`.
pub fn avrami_solid_fraction(t: f64, tau: f64, n: f64) -> f64 {
    if t <= 0.0 {
        return 0.0;
    }
    -(-(t / tau).powf(n)).exp_m1()
}

/// Avrami release rate `dX/dt` (per unit time), for the latent-heat source.
pub fn avrami_rate(t: f64, tau: f64, n: f64) -> f64 {
    if t <= 0.0 {
        return 0.0;
    }
    n / tau * (t / tau).powf(n - 1.0) * (-(t / tau).powf(n)).exp()
}

/// Growth-based Avrami solid fraction: `N` nuclei per unit area (2D) or volume
/// (3D), each growing at the tip velocity `v` in every direction, give an
/// extended measure `g N (v t)^d`; the actual fraction is
/// `X = 1 - exp(-g N (v t)^d)` (Kolmogorov–Johnson–Mehl–Avrami). `dimension` is
/// 2 or 3; `g = pi` in 2D, `4 pi/3` in 3D.
pub fn avrami_from_growth(t: f64, nucleation_density: f64, v: f64, dimension: u32) -> f64 {
    if t <= 0.0 {
        return 0.0;
    }
    let g = if dimension == 2 {
        std::f64::consts::PI
    } else {
        4.0 * std::f64::consts::PI / 3.0
    };
    let ext = g * nucleation_density * (v * t).powi(dimension as i32);
    -(-ext).exp_m1()
}

/// `dX/dt` of [`avrami_from_growth`], the latent-heat release rate per unit
/// time (multiply by `L` for W/m³).
pub fn avrami_rate_from_growth(t: f64, nucleation_density: f64, v: f64, dimension: u32) -> f64 {
    if t <= 0.0 {
        return 0.0;
    }
    let g = if dimension == 2 {
        std::f64::consts::PI
    } else {
        4.0 * std::f64::consts::PI / 3.0
    };
    let ext = g * nucleation_density * (v * t).powi(dimension as i32);
    let dext = g
        * nucleation_density
        * dimension as f64
        * v.powi(dimension as i32)
        * t.powi(dimension as i32 - 1);
    dext * (-ext).exp()
}

/// Is the microstructure resolved by a macroscopic cell of size `dx_macro`?
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Resolution {
    /// Tip diameter `>= cells_per_feature * dx_macro`: resolve it.
    Resolved,
    /// Tip diameter below that: use the closure.
    Closure,
}

/// Indicator for when to unfold: compare the tip diameter against the
/// macroscopic cell. `cells_per_feature` is the number of cells needed to
/// represent a feature (4–8 is typical).
pub fn microstructure_indicator(
    point: &GrowthPoint,
    dx_macro: f64,
    cells_per_feature: f64,
) -> Resolution {
    if point.tip_diameter() >= cells_per_feature * dx_macro {
        Resolution::Resolved
    } else {
        Resolution::Closure
    }
}

/// A one-line closure for a macroscopic solver: given the resolved undercooling
/// and a total solid fraction `X` from a nucleation model, return the
/// volumetric latent-heat release rate `L dX/dt`. This is the rate the isobaric
/// energy equation needs; `dX/dDelta` comes from differentiating the Avrami law
/// along the known cooling path, which the caller supplies as the Avrami
/// parameters `tau`, `n`.
pub fn latent_release_rate(t: f64, tau: f64, n: f64, latent: f64) -> f64 {
    latent * avrami_rate(t, tau, n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ivantsov_inverts_itself() {
        for &d in &[0.05, 0.2, 0.4, 0.55, 0.8] {
            let p = ivantsov_peclet_2d(d);
            let s = p.sqrt();
            let back = std::f64::consts::PI.sqrt() * s * p.exp() * erfc(s);
            assert!((back - d).abs() < 1e-12, "delta={d} back={back}");
        }
    }

    #[test]
    fn avrami_endpoints_and_monotonicity() {
        assert_eq!(avrami_solid_fraction(0.0, 1.0, 3.0), 0.0);
        assert!(avrami_solid_fraction(1.0, 1.0, 3.0) > 0.6);
        assert!(avrami_solid_fraction(5.0, 1.0, 3.0) > 0.99);
        let mut prev = -1.0;
        for i in 0..50 {
            let t = i as f64 * 0.1;
            let x = avrami_solid_fraction(t, 1.0, 3.0);
            assert!(x >= prev - 1e-15);
            prev = x;
        }
    }

    #[test]
    fn avrami_from_growth_is_consistent() {
        let x = avrami_from_growth(2.0, 1e-3, 0.5, 3);
        let dx = avrami_rate_from_growth(2.0, 1e-3, 0.5, 3);
        let x2 = avrami_from_growth(2.0 + 1e-6, 1e-3, 0.5, 3);
        assert!((dx - (x2 - x) / 1e-6).abs() / dx < 1e-3);
        assert!(x > 0.0 && x < 1.0);
    }

    #[test]
    fn growth_is_physically_scaled() {
        // v d0/D = 2 sigma* p^2 and R/d0 = 1/(sigma* p).
        let gp = growth(0.55, 1.0, 2.0, 0.02);
        assert!(gp.v > 0.0 && gp.r > 0.0);
        let lhs = gp.v * 1.0 / 2.0;
        let rhs = 2.0 * 0.02 * gp.pe * gp.pe;
        assert!((lhs - rhs).abs() < 1e-14);
    }

    /// The closure, calibrated once on Karma-Rappel's published tip speed for
    /// `Delta = 0.55`, `eps4 = 0.05`, reproduces both the velocity and the
    /// independently measured phase-field tip radius of that benchmark.
    ///
    /// The measured radii are the `p = 0` whole-needle parabolic fits of this
    /// work (`analysis/solidification.md` §3(ii)): `R/W0 = 9.31` for the
    /// `D = 2` run at `dx/W0 = 0.4` (`V d0/D = 0.0153`) and `R/W0 = 5.46` for
    /// the `D = 3` run at `dx/W0 = 0.2` (`V d0/D = 0.0172`). The residual
    /// against the closure's sharp-interface point is the finite `dx` and
    /// finite `W0/d0` scatter of the phase-field runs (~10 %).
    #[test]
    fn closure_reproduces_kr_2d_benchmark() {
        let sigma = sigma_star_2d_kr_benchmark();
        assert!((sigma - 0.1288).abs() < 5e-4, "sigma* = {sigma}");

        // Calibration point: KR 1998 Table II (Delta = 0.55, eps4 = 0.05),
        // sharp-interface (Green-function) velocity V d0/D = 0.0170.
        // D = 2 uses the benchmark d0/W0 = 0.2770, D = 3 the KR-convergence
        // value 0.1846.
        for (d0w, d) in [(0.27696_f64, 2.0_f64), (0.18464, 3.0)] {
            let gp = growth(0.55, d0w, d, sigma);
            let v = gp.v * d0w / d; // V d0/D
            assert!((v - 0.0170).abs() < 1e-3, "V d0/D = {v}");
            // R = 2 D p/v exactly, by construction.
            let r_iv = 2.0 * d * gp.pe / gp.v;
            assert!((gp.r - r_iv).abs() < 1e-12);
        }

        // Measured phase-field operating points (p = 0 whole-needle fit):
        // (d0/W0, D, V d0/D, R/W0). Each satisfies the 2D Ivantsov relation
        // to within the fit uncertainty and yields the calibrated selection
        // constant within the finite-dx scatter.
        let measured = [(0.27696, 2.0, 0.01530, 9.31), (0.18464, 3.0, 0.01720, 5.46)];
        let p_iv = ivantsov_peclet_2d(0.55);
        for (d0w, d, vd, rw) in measured {
            let p_meas = vd * rw / (2.0 * d0w);
            assert!(
                (p_meas - p_iv).abs() / p_iv < 0.05,
                "p = {p_meas} vs Ivantsov {p_iv}"
            );
            let sig = vd / (2.0 * p_meas * p_meas);
            assert!(
                (sig - sigma).abs() / sigma < 0.12,
                "sigma* = {sig} vs {sigma}"
            );
            // The closure at the measured tip radius and the measured velocity
            // is the same operating point (round trip through selection_sigma).
            let v_phys = gp_v_units(vd, d, d0w);
            let back = selection_sigma(v_phys, rw, d0w, d);
            assert!((back - sig).abs() / sig < 1e-10);
        }
    }

    /// `V d0/D` -> `v` in `d/d0` units with `W0 = 1`.
    fn gp_v_units(vd: f64, d: f64, d0w: f64) -> f64 {
        vd * d / d0w
    }
}
