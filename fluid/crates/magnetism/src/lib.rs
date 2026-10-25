//! Magnetism ladder: atomistic spins -> micromagnetics -> hysteresis
//! (`analysis/emergence_ladder.md`, `analysis/magnetism.md`).
//!
//! The fine rung is the classical Heisenberg model on a crystal lattice:
//! exchange `J`, uniaxial anisotropy `K` and a Zeeman field.  The coarse
//! rung is micromagnetics, where the same physics is carried by the
//! exchange stiffness `A`, the anisotropy constant `K` and the saturation
//! magnetisation `M_s`.  The fold between them is the long-wavelength
//! expansion of the spin Hamiltonian (and, for domain walls, the
//! continuum profile `theta = 2 atan(exp(x/delta))`).  A second fold
//! collapses a particle to a single macrospin, valid while the particle is
//! smaller than the single-domain size set by the exchange length and wall
//! width.
//!
//! - [`geometry`]: sc/bcc/fcc conventional cells, nearest-neighbour tables
//!   and the atomistic exchange stiffness `A(J, S, a)`;
//! - [`spin`]: Metropolis Monte Carlo, the Binder cumulant, mean-field and
//!   Monte-Carlo Curie temperatures, linear spin-wave theory and a spin
//!   chain LLG integrator;
//! - [`micromag`]: the 1D finite-difference Landau-Lifshitz-Gilbert solver
//!   (exchange + uniaxial anisotropy + Zeeman + thin-film demagnetising
//!   field), the exact Bloch wall, Walker breakdown and Kittel FMR;
//! - [`hysteresis`]: the Stoner-Wohlfarth astroid, random-orientation
//!   coercivity and the single-domain validity indicator;
//! - [`demag`]: the Newell-Williams-Dunlop cuboid demagnetising tensor and
//!   a zero-padded FFT demagnetising-field operator (the non-local unfold of
//!   the local-shape closure);
//! - [`llg3d`]: a 3-D finite-difference LLG solver built on that operator;
//! - [`sp4`]: the µMAG standard problem #4 reversal as a reusable run;
//! - [`materials`]: sourced constants for Fe, Ni, permalloy and NdFeB.
//!
//! All benchmarks are checked against exact or published results and the
//! measured errors are reported in `analysis/magnetism.md`; see
//! `tests/atomistic.rs`, `tests/micromag.rs` and `tests/hysteresis.rs`.

pub mod demag;
pub mod geometry;
pub mod hysteresis;
pub mod llg3d;
pub mod materials;
pub mod micromag;
pub mod rng;
pub mod sp4;
pub mod spin;

pub use geometry::{exchange_stiffness, LatticeKind};
pub use spin::{CuriePoint, Heisenberg};
