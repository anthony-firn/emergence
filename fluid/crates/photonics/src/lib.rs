//! Wave optics of layered and structured materials, and the light-matter
//! coupling that follows from it.
//!
//! This crate is the rung below [`radiation`](../radiation/index.html), which
//! is geometrical optics only. Geometrical optics traces rays and multiplies
//! interface reflectances; it is exactly what remains of wave optics when the
//! coherence of the light and the sub-wavelength structure of the matter are
//! both unresolved. This crate makes those two assumptions explicit, and
//! provides the cheap indicator and the unfold for each.
//!
//! ```text
//!   geometrical optics / incoherent    <-- fold up: phase lost, average intensities
//!          ^   indicator n d / L_c
//!   multilayer transfer matrix         <-- exact, coherent, any angle, lossy,
//!          ^                                including complex magnetic mu
//!   photonic crystal band structure    <-- infinite periodic medium, one unit cell
//! ```
//!
//! # What is exact here
//!
//! * `multilayer`: the transfer-matrix method for planar isotropic stacks at
//!   any angle and either polarisation, with complex `N = n + i kappa`.
//!   Outputs `r, t, R, T, A`, per-layer absorptance, and the field used to
//!   deposit absorbed power as a heat source (W m^-3).
//! * `bloch`: the exact band structure of an infinite 1D photonic crystal from
//!   the unit-cell transfer matrix, `cos(K Lambda) = (1/2) tr M`.
//! * `pwe`: the band structure of an infinite 2D crystal of rods or holes by
//!   the plane-wave expansion method, with the inverse-`eps` (Ho) rule, and the
//!   two published 2D benchmarks (square rods, triangular air holes).
//! * `dispersion`: causal material models (Sellmeier, Drude, Lorentz,
//!   split-ring) with sourced constants.
//! * `metamaterial`: double-negative (negative-index) media with a magnetic
//!   response, the Veselago/Pendry perfect lens and its loss-limited resolution.
//!
//! # The emergent boundaries
//!
//! * `effective_medium`: fold a sub-wavelength laminate or mixture into a
//!   homogeneous (possibly uniaxial) medium; indicator `Lambda/lambda`.
//! * `incoherent`: unfold a thick film whose optical thickness exceeds the
//!   source coherence length to the geometrical-optics intensity sum;
//!   indicator `n d / L_c`.
//! * `multilayer::Solution::power_density`: the closure that hands absorbed
//!   optical power to a thermal solver as a volumetric source.
//!
//! The organising idea is in `analysis/emergence_ladder.md`; the wave-optics
//! instance is written up in `analysis/photonics.md`.
pub mod berreman;
pub mod bloch;
pub mod complex;
pub mod dispersion;
pub mod effective_medium;
pub mod incoherent;
pub mod metamaterial;
pub mod mie;
pub mod mie_ladder;
pub mod multilayer;
pub mod plasmon;
pub mod polariton;
pub mod pwe;

pub use complex::C;
pub use dispersion::Dispersion;
pub use multilayer::{Layer, Polarization, Solution, Stack};

/// Vacuum speed of light [m/s] (exact, SI).
pub const C0: f64 = 299_792_458.0;
/// Vacuum permeability [H/m] (CODATA 2018).
pub const MU0: f64 = 1.256_637_062_12e-6;
/// Vacuum permittivity [F/m] (CODATA 2018).
pub const EPS0: f64 = 8.854_187_812_8e-12;
/// Free-space impedance Z0 in ohms.
pub const Z0: f64 = 376.730_313_668;
