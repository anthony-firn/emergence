//! Particles meeting fluids: suspensions and granular matter, and the fold
//! to a continuum rheology (`analysis/emergence_ladder.md`).
//!
//! This crate is the statistical/particulate rung of the ladder. Two fine
//! models live here, both cheap enough to run in a unit test:
//!
//! * **suspensions** — the exact Stokes flow around a sphere
//!   ([`stokes`]), Einstein's effective viscosity obtained from the exact
//!   single-sphere stresslet, the Batchelor–Green `O(phi^2)` coefficient
//!   and the Krieger–Dougherty closure ([`suspension`]);
//! * **granular matter** — a soft-sphere discrete-element method with a
//!   Hertz–Mindlin contact, viscous damping and Coulomb friction
//!   ([`contact`], [`dem`]), benchmarked on the Hertz contact time, the
//!   coefficient of restitution and a Lees–Edwards shear cell.
//!
//! The **fold** is the continuum rheology the particle model is replaced by
//! when the particles are unresolved: the `mu(I)` law of Jop, Forterre &
//! Pouliquen (2006) for dry granular flow, and the `mu(J)` law of Boyer,
//! Guazzelli & Pouliquen (2011) for the fluid-immersed (viscous-number)
//! case. Both the closures and the indicators that bound their validity are
//! in [`rheology`].
//!
//! All quantities are in **SI** unless a module says otherwise. The DEM
//! uses a self-consistent set of units fixed by the particle diameter
//! `d = 1`, density `rho_p = 1`, so stresses are `rho_p d^2/T^2` etc.; the
//! mapping is stated in [`dem`].
//!
//! See `analysis/granular.md` for the benchmarks, sources and the honest
//! misses.

pub mod contact;
pub mod dem;
pub mod linalg;
pub mod rheology;
pub mod rng;
pub mod stokes;
pub mod suspension;

pub use contact::{Contact, ContactParams};
pub use dem::{DemParams, Particle, ShearCell};
pub use rheology::{GranularRegime, Jop2006, MuI, Boyer2011};
pub use suspension::{SuspensionRegime, krieger_dougherty};
