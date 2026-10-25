//! Reversible cyclic voltammetry and the Randles–Ševčík peak.
//!
//! For a reversible one-electron couple with equal diffusion coefficients the
//! total concentration `a + b` is uniform (`= 1`) and the Nernst relation at
//! the electrode fixes a **time-dependent Dirichlet** surface value
//!
//! ```text
//!   a(0,τ) = 1 / (1 + e^{−σ}),   σ = nF(E−E°)/RT
//! ```
//!
//! (the ratio `a/b = e^σ` with `a+b=1`). The dimensionless scan `σ = σ₀ − τ`
//! then gives the classical CV, whose peak flux is the Randles–Ševčík number
//!
//! ```text
//!   j_p = 0.4463      (dimensionless, semi-infinite)
//!   i_p = 0.4463 nFAc √(nFvD/RT)   (A)
//! ```
//!
//! Sources: Bard & Faulkner, *Electrochemical Methods*, 2nd ed., eq. 6.2.12
//! and Table 6.2.1; the integral-equation solution dates to Nicholson & Shain,
//! *Anal. Chem.* **36**, 706 (1964).

use crate::diffusion::{Boundary, Diffusion1D};
use crate::{F_FARADAY, R_GAS};

/// Randles–Ševčík dimensionless peak `0.4463` for a reversible n-electron
/// couple with equal diffusion coefficients.
pub const RANDLES_SEVCIK_CONST: f64 = 0.4463;

/// Randles–Ševčík peak current (A): `i_p = 0.4463 n F A c √(n F v D/(R T))`.
pub fn randles_sevcik_peak(n_e: f64, area: f64, c: f64, d: f64, v: f64, t: f64) -> f64 {
    RANDLES_SEVCIK_CONST * n_e * F_FARADAY * area * c * (n_e * F_FARADAY * v * d / (R_GAS * t)).sqrt()
}

/// A reversible CV scan trace: `(τ, σ, j)` with the dimensionless surface flux
/// `j = ∂a/∂X|₀`.
pub fn reversible_cv(n: usize, x_max: f64, d_tau: f64, sigma_start: f64, tau_end: f64) -> Vec<(f64, f64, f64)> {
    // Initial condition: a = 1 everywhere.
    let mut s = Diffusion1D::uniform(n, x_max, 1.0, Boundary::Dirichlet(1.0), Boundary::NoFlux, 1.0);
    let steps = (tau_end / d_tau).round() as usize;
    let mut out = Vec::with_capacity(steps + 1);
    for k in 0..=steps {
        let tau = k as f64 * d_tau;
        let sigma = sigma_start - tau;
        let a0 = 1.0 / (1.0 + (-sigma).exp());
        let c = s.concentrations();
        // Second-order one-sided surface gradient (the boundary value a0 is
        // exact), so the peak converges at O(ΔX²) rather than O(ΔX) as with
        // the simple (a1−a0)/ΔX difference.
        let j = (-3.0 * a0 + 4.0 * c[1] - c[2]) / (2.0 * s.dx);
        out.push((tau, sigma, j));
        if k < steps {
            // Enforce the surface value at the new time; the previous step's
            // value enters the explicit CN side, giving the centred boundary.
            let sigma = sigma_start - (tau + d_tau);
            s.set_left(Boundary::Dirichlet(1.0 / (1.0 + (-sigma).exp())));
            s.step(d_tau);
        }
    }
    out
}

/// Peak dimensionless flux of a reversible CV scan.
pub fn reversible_cv_peak(n: usize, x_max: f64, d_tau: f64, sigma_start: f64, tau_end: f64) -> f64 {
    reversible_cv(n, x_max, d_tau, sigma_start, tau_end)
        .iter()
        .map(|&(_, _, j)| j)
        .fold(0.0f64, f64::max)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peak_is_randles_sevcik() {
        // Semi-infinite limit: start well positive of E°, X_max large compared
        // with the diffusion length at the peak (√τ_p ≈ 3.3).
        let peak = reversible_cv_peak(1200, 60.0, 2e-3, 10.0, 22.0);
        let rel = (peak - RANDLES_SEVCIK_CONST).abs() / RANDLES_SEVCIK_CONST;
        assert!(rel < 5e-3, "CV peak {peak} vs 0.4463 (rel {rel})");
    }

    #[test]
    fn peak_is_quadratically_converged_in_dtau() {
        let p1 = reversible_cv_peak(800, 60.0, 4e-3, 10.0, 22.0);
        let p2 = reversible_cv_peak(800, 60.0, 2e-3, 10.0, 22.0);
        assert!((p1 - p2).abs() < 1e-3, "dtau sensitivity {p1} vs {p2}");
    }

    #[test]
    fn peak_scales_with_sqrt_scan_rate_and_n_to_3_2() {
        // i_p ∝ √v at fixed n, and ∝ n^{3/2} at fixed v.
        let p1 = randles_sevcik_peak(1.0, 1.0, 1.0, 1e-5, 0.1, 298.15);
        let p4 = randles_sevcik_peak(1.0, 1.0, 1.0, 1e-5, 0.4, 298.15);
        assert!((p4 / p1 - 2.0).abs() < 1e-12);
        let pn = randles_sevcik_peak(2.0, 1.0, 1.0, 1e-5, 0.1, 298.15);
        assert!((pn / p1 - 2f64.powf(1.5)).abs() < 1e-10);
    }
}
