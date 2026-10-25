//! Elastic and acoustic waves in layered and structured materials, and their
//! homogenised (folded) descriptions.
//!
//! This is the mechanical counterpart of [`photonics`](../photonics/index.html):
//! where that crate keeps the phase of multiply reflected light, this one keeps
//! the phase of multiply reflected sound and the sub-wavelength structure of
//! the matter. Both sit below a geometrical/ray rung and provide the cheap
//! indicator and the fold/unfold for the wave structure.
//!
//! ```text
//!   ray acoustics / effective medium   <-- fold up: period << wavelength
//!          ^   indicator (omega Lambda / c)^2
//!   transfer matrix (layered)          <-- exact, coherent, lossy
//!          ^   indicator N |Im q| Lambda
//!   Bloch band structure               <-- infinite periodic medium, one cell
//!   locally resonant metamaterial      <-- inner resonance, sub-wavelength gap
//! ```
//!
//! # What is exact here
//!
//! * [`layered`]: transfer matrices for 1D elastic/acoustic stacks at normal
//!   incidence, with optional loss through a complex modulus; `r, t, R, T, A`;
//!   and the exact Rytov Bloch dispersion of the infinite periodic medium.
//! * [`homogenize`]: the exact long-wavelength laminate modulus (Reuss) and
//!   density, the symmetric Voigt/Reuss/Hashin–Shtrikman bounds for isotropic
//!   two-phase composites, and the phase-velocity validity indicator.
//! * [`resonant`]: the mass-in-mass lattice, its exact and effective-mass
//!   dispersions and its sub-wavelength gap.
//! * [`fluid_coupling`]: oblique fluid–solid–fluid transmission through an
//!   elastic plate with P and SV waves, showing coincidence and Lamb peaks.
//!
//! # The emergent boundaries
//!
//! * Fold a sub-wavelength laminate into an effective medium; the indicator is
//!   the relative phase-velocity error `~ (omega Lambda / c_eff)^2`.
//! * Fold a chain of locally resonant cells below `omega0` into a heavy chain
//!   of mass `m1 + m2`; unfold to the two-degree-of-freedom cell near `omega0`.
//! * Unfold a finite stack at long wavelength back to [`layered`] when the
//!   Bloch decay over the stack is small.
//!
//! The organising idea is in `analysis/emergence_ladder.md`; the acoustics
//! instance is written up in `analysis/phononic.md`.
pub mod complex;
pub mod fluid_coupling;
pub mod helmholtz;
pub mod homogenize;
pub mod lamb;
pub mod layered;
pub mod materials;
pub mod pwe;
pub mod resonant;

pub use complex::C;
pub use layered::{Layer, Mat2, Periodic, Solution, Stack};
