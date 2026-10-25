//! Electrochemistry and batteries across the emergent levels of
//! [`analysis/emergence_ladder.md`](../../../analysis/emergence_ladder.md).
//!
//! The ladder for a cell, from the cheapest rung to the reference:
//!
//! ```text
//!  equivalent circuit (Thevenin RC)         <- BMS / system level
//!        ^  indicator: voltage error vs pulse length / C-rate
//!  single-particle model (SPM)              <- cell control level
//!        ^  indicator: electrolyte diffusion time / discharge time
//!  electrolyte-enhanced SPM / DFN           <- the reference cell model
//!        ^
//!  interface kinetics + 1D transport        <- the electrochemical rung
//! ```
//!
//! Each module is a rung or a benchmark of one:
//!
//! * [`interface`] — Nernst potential, Butler–Volmer, its Tafel and linear
//!   folds, and the mass-transport-limited current.
//! * [`diffusion`] — a Crank–Nicolson 1D diffusion solver and the exact
//!   **Cottrell** transient `i(t) = nFAc√(D/(πt))` for a potential step.
//! * [`cv`] — a dimensionless reversible cyclic-voltammetry solver and the
//!   **Randles–Ševčík** peak current `i_p = 0.4463 nFAc√(nFvD/(RT))`.
//! * [`particle`] — spherical solid diffusion (finite volume) against the
//!   exact series for a **constant surface flux**.
//! * [`spm`] — the single-particle cell model with a published parameter set
//!   and Butler–Volmer kinetics; conservation and capacity vs C-rate.
//! * [`spme`] — the electrolyte-enhanced SPM that restores the electrolyte
//!   concentration (the finer rung when the SPM indicator is too large).
//! * [`ecm`] — the Thevenin equivalent-circuit fold, its parameters extracted
//!   from an SPM/DFN pulse response, and the fold error.
//! * [`thermal`] — reversible and irreversible heat generation
//!   `q = I(U−V) − I T dU/dT` and the interface to the `heat` crate.
//!
//! # Conventions
//! Current `I > 0` is **discharge** (lithium leaves the negative particle,
//! enters the positive particle). Potentials are vs Li/Li⁺ unless stated.
//!
//! # References
//! * A. J. Bard, L. R. Faulkner, *Electrochemical Methods*, 2nd ed. (2001) —
//!   Nernst, Butler–Volmer, Tafel, Cottrell, Randles–Ševčík.
//! * J. Crank, *The Mathematics of Diffusion*, 2nd ed. (1975) — sphere series.
//! * M. Doyle, T. F. Fuller, J. Newman, *J. Electrochem. Soc.* **140**, 1526
//!   (1993) — the DFN model.
//! * C.-H. Chen et al., *J. Electrochem. Soc.* **167**, 080534 (2020) — LG M50
//!   parameter set (noted `Chen2020`).

pub mod cv;
pub mod dfn;
pub mod diffusion;
pub mod ecm;
pub mod interface;
pub mod particle;
pub mod spm;
pub mod spme;
pub mod thermal;

/// Universal gas constant (J mol⁻¹ K⁻¹), CODATA 2018.
pub const R_GAS: f64 = 8.314_462_618;
/// Faraday constant (C mol⁻¹), CODATA 2018.
pub const F_FARADAY: f64 = 96_485.332_12;
/// Room temperature used unless a call passes another value (K).
pub const T_REF: f64 = 298.15;

/// Thermal voltage `RT/F` (V) at temperature `t`.
pub fn thermal_voltage(t: f64) -> f64 {
    R_GAS * t / F_FARADAY
}
