//! Viscoelastic test: a virtual DMA frequency sweep on a polystyrene-like
//! melt, with the Williams–Landel–Ferry time–temperature shift.
//!
//! The coarse model is a multi-mode Maxwell (Prony) spectrum from the Rouse
//! modes of the polymer crate; the virtual instrument integrates the Maxwell
//! ODEs under `γ(t) = γ₀ sin ωt` and projects the steady stress onto
//! `sin`/`cos`, which must reproduce the analytic `G*(ω)`. The Deborah
//! number `De = λ₁ω` is the validity indicator: the linear-viscoelastic
//! closure is certified only where the slowest mode is excited
//! (`0.1 < De < 10`); below it the fluid is Newtonian.
//!
//! Sources: Ferry, *Viscoelastic Properties of Polymers*, 3rd ed. (1980) —
//! polystyrene `T_g = 373 K`, plateau modulus `G_N⁰ ≈ 2.0×10⁵ Pa`;
//! Williams, Landel & Ferry (1955) — `C₁ = 17.44`, `C₂ = 51.6 K`.

use crate::{Level, Measurement};
use polymer::continuum::PronySeries;
use polymer::indicator::{RheologyRegime, ViscoelasticIndicator};
use polymer::rheology::{polystyrene_melt, shift_prony, Rheometer, Wlf};

/// The polystyrene spectrum at the reference temperature, with its
/// temperature-dependent shift.
pub struct Polystyrene {
    pub t_g: f64,
    pub t_ref: f64,
}

impl Polystyrene {
    pub fn new() -> Self {
        Self { t_g: 373.0, t_ref: 453.0 }
    }
    pub fn wlf(&self) -> Wlf {
        // The universal WLF constants are referred to T_g; the reference
        // temperature is where polystyrene_melt()'s times are quoted.
        let mut w = Wlf::universal(self.t_g);
        w.t_ref = self.t_ref;
        w
    }
    pub fn spectrum(&self) -> PronySeries {
        polystyrene_melt()
    }
}

impl Default for Polystyrene {
    fn default() -> Self {
        Self::new()
    }
}

/// The DMA sweep and the WLF/TTS check.
pub fn dma() -> Vec<Measurement> {
    let mut out = Vec::new();
    let ps = Polystyrene::new();
    let base = ps.spectrum();
    let rheo = Rheometer::new(base.clone(), 0.0);
    let lambda1 = base.modes.iter().map(|&(_, l)| l).fold(0.0_f64, f64::max);

    // Frequency sweep. The ODE instrument uses a fine step; the analytic
    // complex modulus is the reference.
    for &omega in &[1.0_f64, 10.0, 100.0] {
        let (gp_meas, gpp_meas) = rheo.oscillatory_from_ode(1e-3, omega, 0.002, 40);
        let (gp_ref, gpp_ref) = rheo.complex(omega);
        let de = lambda1 * omega;
        let regime = ViscoelasticIndicator::oscillatory(lambda1, omega).regime;
        // Report the storage modulus; the loss modulus is the second row.
        out.push(Measurement::new(
            "polystyrene",
            "DMA sweep",
            Level::Macroscopic,
            &format!("G'({omega:.0} rad/s)"),
            gp_meas,
            gp_ref,
            "Pa",
            regime == RheologyRegime::Viscoelastic,
            &format!(
                "virtual DMA: Maxwell ODEs under γ = γ₀ sin ωt projected onto sin/cos, vs the \
                 analytic Prony G*(ω). De = λ₁ω = {de:.2}, regime {regime:?}.",
            ),
        ));
        out.push(Measurement::new(
            "polystyrene",
            "DMA sweep",
            Level::Macroscopic,
            &format!("G''({omega:.0} rad/s)"),
            gpp_meas,
            gpp_ref,
            "Pa",
            regime == RheologyRegime::Viscoelastic,
            "loss modulus from the same virtual DMA; the reference is the analytic Prony G''(ω)",
        ));
    }
    // Plateau modulus G_N⁰ = sum of the mode weights.
    let gn0: f64 = base.modes.iter().map(|&(g, _)| g).sum();
    out.push(Measurement::new(
        "polystyrene",
        "DMA sweep",
        Level::Macroscopic,
        "G_N^0",
        gn0,
        2.0e5,
        "Pa",
        (gn0 / 2.0e5 - 1.0).abs() < 0.05,
        "plateau modulus as the sum of the Prony weights vs Ferry (1980), G_N⁰ ≈ 2.0×10⁵ Pa. \
         The absolute time scale of the spectrum is an estimate; the plateau and the WLF shift \
         are the sourced numbers.",
    ));

    // WLF shift at T_g + 50 K, and the exact TTS collapse identity.
    let w = Wlf::universal(ps.t_g);
    let log_a = w.log10_a(ps.t_g + 50.0);
    out.push(Measurement::new(
        "polystyrene",
        "WLF shift",
        Level::Macroscopic,
        "log10 a_T (T_g+50 K)",
        log_a,
        -8.5827,
        "",
        (log_a + 8.5827).abs() < 0.01,
        "Williams–Landel–Ferry with the universal C₁ = 17.44, C₂ = 51.6 K at T_g = 373 K; the \
         reference is the same universal law quoted in analysis/polymer.md",
    ));
    // TTS: G*(ω; T) = G*(a_T ω; T_ref) by construction; report the residual.
    let a = w.a(ps.t_ref + 30.0);
    let shifted = PronySeries::new(shift_prony(&base.modes, a));
    let omega = 3.0;
    let (gp_t, _gpp_t) = shifted.complex(omega);
    let (gp_ref, _gpp_ref) = base.complex(a * omega);
    let collapse = (gp_t - gp_ref).abs() / gp_ref.abs().max(1e-300);
    out.push(Measurement::new(
        "polystyrene",
        "WLF shift",
        Level::Macroscopic,
        "TTS collapse residual",
        collapse,
        0.0,
        "",
        collapse < 1e-12,
        "time–temperature superposition G*(ω; T) = G*(a_T ω; T_ref), exact by construction for a \
         horizontally shifted Prony series (no vertical shift)",
    ));
    out
}
