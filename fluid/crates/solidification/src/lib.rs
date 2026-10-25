//! Crystal growth from the melt: anisotropic phase-field solidification, its
//! sharp-interface and isobaric folds, and the mushy-zone closure.
//!
//! This crate is one rung of the emergence ladder in
//! `analysis/emergence_ladder.md` and fills the item left open in
//! `analysis/compressible.md` §15 ("dendritic ice and solidification
//! microstructure":
//!
//! ```text
//!   isobaric enthalpy model (macroscopic cell, no interface)
//!        ^  indicator: tip radius / dendrite spacing vs cell size
//!   Karma-Rappel anisotropic phase field  <-- this crate
//!        ^  thin-interface limit W0/d0 -> 0
//!   sharp-interface Stefan + Gibbs-Thomson
//!        ^  exact
//!   1D Neumann similarity solution
//! ```
//!
//! * [`phase_field`] — 2D Karma–Rappel (1996/1998) phase-field solidification
//!   of a pure substance with four-fold surface-energy anisotropy, explicit
//!   time stepping and rayon-parallel rows. [`PhaseField::step_kr`] is the
//!   corrected operator (KR Eq. 88 variational anisotropy term);
//!   [`PhaseField::step`] is the plain operator kept for the GPU = CPU check.
//! * [`stefan`] — the exact Neumann similarity solution and a 1D enthalpy
//!   solver used as its reference.
//! * [`closure`] — Ivantsov + microscopic-solvability tip kinetics, an Avrami
//!   mushy-zone law, and the microstructure-resolution indicator that decides
//!   when a macroscopic solver must unfold.
//! * [`materials`] — succinonitrile and ice/water properties with sources.
//!
//! # Quick start
//!
//! ```
//! use solidification::phase_field::PhaseField;
//! use solidification::stefan::neumann_lambda;
//!
//! // A quarter dendrite: W0/d0 = 0.4, four-fold anisotropy 0.05.
//! let mut pf = PhaseField::new(64, 64, 0.25, 1.0, 2.0, 0.4, 0.05);
//! pf.u_far = -0.55;
//! pf.seed_corner(2.0);
//! let lambda = neumann_lambda(0.5);
//! assert!(lambda > 0.0);
//! ```

pub mod closure;
pub mod materials;
pub mod phase_field;
pub mod stefan;

pub use closure::{
    avrami_from_growth, avrami_rate, avrami_rate_from_growth, avrami_solid_fraction, growth,
    ivantsov_peclet_2d, latent_release_rate, microstructure_indicator, GrowthPoint, Resolution,
};
pub use materials::{MeltMaterial, ICE, SUCCINONITRILE};
pub use phase_field::{Bc, Edge, PhaseField, A1, A2};
pub use stefan::{
    erf, erfc, neumann_front, neumann_lambda, neumann_residual, neumann_temperature, Stefan1d,
};

#[cfg(test)]
mod tests {
    /// The Karma–Rappel constants have the values quoted in the source paper.
    #[test]
    fn karma_rappel_constants() {
        assert!((super::A1 - 5.0 * std::f64::consts::SQRT_2 / 8.0).abs() < 1e-12);
        assert!((super::A2 - 0.6267).abs() < 1e-12);
    }
}
