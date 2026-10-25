//! Materials ladder, atomistic rung: from atoms to continuum elasticity
//! (`analysis/emergence_ladder.md`).
//!
//! The coarse model here is continuum elasticity: a homogeneous strain
//! field, the 6x6 Voigt stiffness `C_ij`, and the Christoffel sound speeds
//! it predicts. The fine model is the atomistic crystal with a pair
//! potential. The fold between them is the Cauchy-Born rule (affine
//! deformation plus, for non-Bravais lattices, internal relaxation of the
//! basis), and the validity indicator is the relative dispersion error
//! `(ω_lattice - v k)/ω_lattice` as a function of `k·a`.
//!
//! - [`lattice`]: primitive vectors + basis, neighbour shells and bonds;
//! - [`potential`]: Lennard-Jones/Morse pair potentials and their exact
//!   static lattice sums;
//! - [`elastic`]: Cauchy-Born 6x6 stiffness with internal relaxation, and a
//!   finite-difference check of the energy under strain;
//! - [`phonon`]: the harmonic dynamical matrix `D(k)` and its frequencies;
//! - [`continuum`]: the Christoffel equation and the dispersion-error
//!   indicator that says when the continuum rung may be used;
//! - [`anisotropic`]: real single-crystal elastic constants, symmetry
//!   classes, Bond rotation, the compliance and directional properties;
//! - [`polycrystal`]: the fold to the isotropic continuum (Voigt/Reuss/Hill,
//!   Hashin-Shtrikman bounds) and the grain-size/wavelength validity
//!   indicator.
//!
//! The rungs are validated against exact lattice sums and exact
//! 1D-chain/dispersion identities; see `tests/atomistic.rs`,
//! `tests/anisotropic.rs`, `tests/polycrystal.rs` and
//! `analysis/materials.md`.

pub mod anisotropic;
pub mod continuum;
pub mod defects;
pub mod elastic;
pub mod fracture;
pub mod lattice;
pub mod linalg;
pub mod phonon;
pub mod plasticity;
pub mod polycrystal;
pub mod potential;
pub mod quasicontinuum;
pub mod quasiharmonic;

pub use anisotropic::{CrystalMaterial, Symmetry, M6};
pub use lattice::{Bond, Crystal, Lattice, Shell, V3};
pub use potential::{Harmonic, LennardJones, Morse, PairPotential};
