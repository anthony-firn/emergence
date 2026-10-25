//! Molecular dynamics: the statistical rung of the emergence ladder
//! (`analysis/emergence_ladder.md`).
//!
//! The long-term goal of this crate is the fold from atoms to the
//! continuum: the equation of state and the transport coefficients the
//! Navier–Stokes solvers use as closures, together with a cheap indicator
//! that says when those closures are valid.
//!
//! * [`potential`] — Lennard-Jones and Morse pairs, three cutoff schemes
//!   (shifted, shifted-force, tail-corrected);
//! * [`neighbor`] — cell-linked, minimum-image neighbour lists;
//! * [`system`] — the periodic box, force/stress/heat-flux kernels and
//!   velocity-Verlet integration;
//! * [`thermostat`] — Bussi stochastic velocity rescaling, Berendsen and
//!   Nosé–Hoover;
//! * [`eos`] — the Johnson–Zollweg–Gubbins (1993) Lennard-Jones equation
//!   of state used as the thermodynamic reference;
//! * [`transport`] — Green–Kubo shear viscosity and thermal conductivity,
//!   self-diffusion from MSD and VACF, and the Yeh–Hummer finite-size
//!   correction;
//! * [`argon`] — reduced-to-SI closures for argon and measured liquid-argon
//!   values near 90 K;
//! * [`validity`] — Knudsen and scale-separation indicators for the
//!   continuum limit.
//!
//! Everything runs in reduced Lennard-Jones units unless the [`argon`]
//! module is used. See `analysis/md.md` for the benchmarks, with error
//! bars.

pub mod argon;
pub mod barostat;
pub mod eos;
pub mod neighbor;
pub mod potential;
pub mod rng;
pub mod system;
pub mod thermostat;
pub mod transport;
pub mod validity;

pub use barostat::BerendsenBarostat;
pub use eos::JzgEos;
pub use potential::{Cutoff, Potential};
pub use rng::Rng;
pub use system::System;
pub use thermostat::Thermostat;

/// Argon-like Lennard-Jones parameters in reduced units (`sigma = 1`,
/// `epsilon = 1`).
pub fn lj() -> Potential {
    Potential::LennardJones(materials::LennardJones::new(1.0, 1.0))
}
