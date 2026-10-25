//! A virtual DMA / rheometer and time–temperature superposition (WLF).
//!
//! The linear-viscoelastic response is a Prony (multi-mode Maxwell) series
//! `G(t) = Σ_i g_i e^{-t/λ_i}`. A virtual instrument applies a shear
//! history `γ(t)` and integrates the diagonal Maxwell ODEs
//! `dσ_i/dt = g_i γ̇ - σ_i/λ_i` (plus a solvent term `η_s γ̇`), then extracts
//! `G'` and `G''` from the steady periodic part.
//!
//! Time–temperature superposition uses the Williams–Landel–Ferry equation
//!
//! ```text
//! log₁₀ a_T = − C₁ (T − T_ref) / (C₂ + T − T_ref),
//! ```
//!
//! with the universal constants `C₁ = 17.44`, `C₂ = 51.6 K` relative to the
//! glass transition `T_g` (M. L. Williams, R. F. Landel & J. D. Ferry,
//! *J. Am. Chem. Soc.* **77**, 3701 (1955)). The master-curve relation is
//! `G'(ω; T) = G'(a_T ω; T_ref)`, i.e. each mode time becomes `a_T λ_i`.
//!
//! Polystyrene reference data (sources in `analysis/polymer.md`):
//! `T_g = 373 K`, plateau modulus `G_N⁰ ≈ 2.0×10⁵ Pa`, entanglement molar
//! mass `M_e ≈ 13.3 kg/mol`.

use crate::continuum::PronySeries;

/// Williams–Landel–Ferry shift factor.
#[derive(Clone, Copy, Debug)]
pub struct Wlf {
    pub c1: f64,
    pub c2: f64,
    pub t_ref: f64,
}

impl Wlf {
    /// Universal constants relative to the glass transition `T_g`.
    pub fn universal(t_g: f64) -> Self {
        Self { c1: 17.44, c2: 51.6, t_ref: t_g }
    }

    /// `log₁₀ a_T`.
    pub fn log10_a(&self, t: f64) -> f64 {
        let dt = t - self.t_ref;
        -self.c1 * dt / (self.c2 + dt)
    }

    /// Shift factor `a_T`. `a_T < 1` above `T_ref` (faster relaxation).
    pub fn a(&self, t: f64) -> f64 {
        10f64.powf(self.log10_a(t))
    }
}

/// Shift a Prony series in temperature: `λ_i → a_T λ_i` (no vertical
/// shift).
pub fn shift_prony(series: &[(f64, f64)], a: f64) -> Vec<(f64, f64)> {
    series.iter().map(|&(g, l)| (g, l * a)).collect()
}

/// A multi-mode Maxwell fluid with optional solvent, driven by a shear
/// history. This is a virtual rheometer.
#[derive(Clone, Debug)]
pub struct Rheometer {
    pub prony: PronySeries,
    pub eta_s: f64,
}

impl Rheometer {
    pub fn new(prony: PronySeries, eta_s: f64) -> Self {
        Self { prony, eta_s }
    }

    /// Stress relaxation after a step strain `γ₀`: `σ(t) = γ₀ G(t)`.
    pub fn stress_relaxation(&self, gamma0: f64, times: &[f64]) -> Vec<f64> {
        times.iter().map(|&t| gamma0 * self.prony.g_t(t)).collect()
    }

    /// Integrate the Maxwell ODEs under `γ(t) = γ₀ sin(ωt)` and extract
    /// `(G', G'')` from the last period by projecting the steady stress
    /// onto `sin` and `cos`. Returns `None` if `dt` is not fine enough.
    pub fn oscillatory_from_ode(
        &self,
        gamma0: f64,
        omega: f64,
        dt: f64,
        n_periods: usize,
    ) -> (f64, f64) {
        let period = std::f64::consts::TAU / omega;
        let steps_per_period = (period / dt).round() as usize;
        let n_steps = steps_per_period * n_periods;
        // State: one stress per mode.
        let mut s: Vec<f64> = vec![0.0; self.prony.modes.len()];
        let mut acc_s = 0.0;
        let mut acc_c = 0.0;
        let mut count = 0.0;
        let last_period_start = (n_periods - 1) * steps_per_period;
        for n in 0..n_steps {
            let t = n as f64 * dt;
            // Stress at time t (before the step) so the projection phase is
            // exact — important for stiff modes whose in-phase response is
            // small compared with their out-of-phase response.
            let gamma_dot_now = gamma0 * omega * (omega * t).cos();
            let sigma: f64 = s.iter().sum::<f64>() + self.eta_s * gamma_dot_now;
            if n >= last_period_start {
                acc_s += sigma * (omega * t).sin();
                acc_c += sigma * (omega * t).cos();
                count += 1.0;
            }
            let gdot = gamma_dot_now;
            let gd_next = gamma0 * omega * (omega * (t + dt)).cos();
            // Exact integrating factor for each independent mode, with the
            // shear rate linear over the step: it is exact for piecewise
            // linear γ̇ and A-stable for any dt/λ.
            for (i, &(g, l)) in self.prony.modes.iter().enumerate() {
                let e = (-dt / l).exp();
                let a = gdot;
                let b = (gd_next - gdot) / dt;
                let i1 = l * (1.0 - e);
                let i2 = l * l * (1.0 - (1.0 + dt / l) * e);
                let integral = (a + b * dt) * i1 - b * i2;
                s[i] = e * s[i] + g * integral;
            }
        }
        // (2/P)∫ σ sin = γ₀ G';  (2/P)∫ σ cos = γ₀ G''.
        (2.0 * acc_s / count / gamma0, 2.0 * acc_c / count / gamma0)
    }

    /// Analytic `(G', G'')` including the solvent.
    pub fn complex(&self, omega: f64) -> (f64, f64) {
        let (gp, gpp) = self.prony.complex(omega);
        (gp, gpp + self.eta_s * omega)
    }
}

/// A polystyrene-like melt spectrum at a reference temperature. The times
/// are Rouse-like (`λ_p = λ_1/p²`) and their absolute scale is an estimate;
/// the shape, plateau modulus and WLF shift are what the tests use.
///
/// Constants: `T_g = 373 K`, `G_N⁰ = 2.0×10⁵ Pa`, `M_e = 13.3 kg/mol`
/// (J. D. Ferry, *Viscoelastic Properties of Polymers*, 3rd ed. 1980;
/// J. Brandrup & E. H. Immergut, *Polymer Handbook*).
pub fn polystyrene_melt() -> PronySeries {
    let g_n0 = 2.0e5;
    let n_modes = 5;
    let tau1 = 0.1; // estimate at ~453 K
    let modes = (1..=n_modes)
        .map(|p| (g_n0 / n_modes as f64, tau1 / (p * p) as f64))
        .collect();
    PronySeries::new(modes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wlf_universal_constants() {
        let w = Wlf::universal(373.15);
        // a_T = 1 at T_g.
        assert!((w.a(373.15) - 1.0).abs() < 1e-12);
        // At 50 K above T_g the universal WLF gives log10 a_T ≈ −8.58.
        let l = w.log10_a(423.15);
        assert!((l - (-8.583)).abs() < 0.01, "log10 a = {l}");
        // Monotone decreasing with temperature.
        assert!(w.a(450.0) < w.a(400.0));
    }

    #[test]
    fn time_temperature_superposition_collapses() {
        let w = Wlf::universal(373.15);
        let base = polystyrene_melt();
        let omega = 3.0;
        for &t in &[400.0, 430.0, 470.0] {
            let a = w.a(t);
            let shifted = PronySeries::new(shift_prony(&base.modes, a));
            // G*(ω; T) = G*(a_T ω; T_ref).
            let (gp_t, gpp_t) = shifted.complex(omega);
            let (gp_ref, gpp_ref) = base.complex(a * omega);
            assert!((gp_t - gp_ref).abs() / gp_ref < 1e-12);
            assert!((gpp_t - gpp_ref).abs() / gpp_ref < 1e-12);
        }
    }

    #[test]
    fn virtual_dma_matches_analytic_moduli() {
        let base = polystyrene_melt();
        let r = Rheometer::new(base.clone(), 0.0);
        for &omega in &[0.1, 1.0, 10.0] {
            let (gp_m, gpp_m) = r.oscillatory_from_ode(1e-3, omega, 0.005, 60);
            let (gp_a, gpp_a) = r.complex(omega);
            assert!((gp_m - gp_a).abs() / gp_a < 0.01, "ω={omega} G' {gp_m} vs {gp_a}");
            assert!((gpp_m - gpp_a).abs() / gpp_a < 0.01, "ω={omega} G'' {gpp_m} vs {gpp_a}");
        }
    }

    #[test]
    fn stress_relaxation_is_g_of_t() {
        let base = polystyrene_melt();
        let r = Rheometer::new(base.clone(), 0.0);
        let gamma = 0.02;
        let times = [0.0, 0.05, 0.2];
        let s = r.stress_relaxation(gamma, &times);
        for (i, &t) in times.iter().enumerate() {
            assert!((s[i] - gamma * base.g_t(t)).abs() < 1e-15);
        }
    }
}

