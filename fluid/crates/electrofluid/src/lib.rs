//! Electrical properties of fluids and how fields change flow and
//! "gloopiness": the electromagnetic rungs of `analysis/emergence_ladder.md`
//! (Maxwell → quasi-static → circuits; Poisson–Nernst–Planck/Debye layer →
//! leaky dielectric → electrorheology).
//!
//! - [`pnp`]: 1D Poisson–Nernst–Planck / Poisson–Boltzmann, Debye length,
//!   Gouy–Chapman and Grahame;
//! - [`eo`]: electro-osmosis, Helmholtz–Smoluchowski slip as the closure that
//!   replaces a resolved Debye layer;
//! - [`leaky`]: leaky-dielectric (Taylor–Melcher) drop, small-deformation
//!   discriminating function and circulation;
//! - [`er`]: electrorheology, point-dipole chains, Bingham fold and Mason
//!   number.

pub mod eo;
pub mod er;
pub mod leaky;
pub mod pnp;

pub use eo::{hs_slip, Channel};
pub use er::{dipole_force, ErSim};
pub use leaky::LeakyDrop;
pub use pnp::{debye_length, grahame_sigma, gouy_chapman, pb_solve};
