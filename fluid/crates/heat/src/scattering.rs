//! Phonon scattering rates: umklapp, normal, point-defect (isotope) and
//! boundary (Casimir), in the Callaway/Holland form.

use crate::phonon::{Branch, HBAR};
use std::f64::consts::PI;

/// Parameters of the Callaway/Holland relaxation-time model.
///
/// The umklapp rate is the branch-resolved Morelli form
/// `τ_U^{-1} = ħ γ² ω² T /(M v² Θ) · exp(−Θ/3T)`
/// (Morelli, Heremans & Slack 2002; Kazan et al. 2010). `γ` and `Θ` live on
/// the [`Branch`]; see [`crate::materials`] for the silicon values.
#[derive(Clone, Copy, Debug)]
pub struct ScatteringParams {
    /// Average atomic mass `M` (kg).
    pub mass: f64,
    /// Normal-process strength: `τ_N^{-1} = β τ_U^{-1}`.
    pub beta: f64,
    /// Mass-fluctuation parameter `Γ = Σ f_i (1 − m_i/M̄)²`.
    pub gamma_iso: f64,
    /// Volume per atom `V_0` (m³), for the Klemens point-defect rate.
    pub volume: f64,
    /// Boundary (Casimir) length `L` (m); `f64::INFINITY` for a bulk crystal.
    pub length: f64,
}

impl ScatteringParams {
    /// Umklapp rate `τ_U^{-1}` (1/s).
    pub fn umklapp(&self, b: &Branch, omega: f64, t: f64) -> f64 {
        if omega <= 0.0 {
            return 0.0;
        }
        let pref = HBAR * b.gamma * b.gamma * omega * omega * t
            / (self.mass * b.v * b.v * b.theta());
        pref * (-b.theta() / (3.0 * t)).exp()
    }

    /// Normal-process rate `τ_N^{-1} = β τ_U^{-1}` (1/s).
    pub fn normal(&self, b: &Branch, omega: f64, t: f64) -> f64 {
        self.beta * self.umklapp(b, omega, t)
    }

    /// Point-defect (isotope) rate, Klemens:
    /// `τ_I^{-1} = V_0 Γ ω⁴ /(4π v³)` (1/s).
    pub fn point_defect(&self, b: &Branch, omega: f64) -> f64 {
        self.volume * self.gamma_iso * omega.powi(4) / (4.0 * PI * b.v.powi(3))
    }

    /// Boundary (Casimir) rate `τ_B^{-1} = v/L` (1/s).
    pub fn boundary(&self, b: &Branch, _omega: f64) -> f64 {
        b.v / self.length
    }

    /// Total (combined) rate `τ_c^{-1}`, including normal processes as a
    /// resistance when `beta > 0` (the Callaway "combined" rate).
    pub fn total(&self, b: &Branch, omega: f64, t: f64) -> f64 {
        self.umklapp(b, omega, t)
            + self.normal(b, omega, t)
            + self.point_defect(b, omega)
            + self.boundary(b, omega)
    }

    /// Resistive rate `τ_R^{-1}` (umklapp + isotope + boundary), excluding the
    /// momentum-conserving normal processes.
    pub fn resistive(&self, b: &Branch, omega: f64, t: f64) -> f64 {
        self.umklapp(b, omega, t) + self.point_defect(b, omega) + self.boundary(b, omega)
    }
}
