//! CPU kinetic closures, not a quantum field / confinement solver.
//! Momenta and energies use GeV (c=1); positions [m], times [s], fields SI.
//! Every unresolved decay retains its four-momentum and additive charges.
#![doc = include_str!("../README.md")]
pub mod coupling;
pub mod decay;
pub mod fields;
pub mod geodesic;
pub mod kinematics;
pub mod qed;
pub mod species;
pub mod thermal;
pub mod transport;
pub mod weak;

pub use fields::UniformField;
pub use kinematics::FourMomentum;
use radiation::math as radiation_math;
pub use radiation::math::{Rng, Vec3};
pub use species::{Domain, Species};
pub use transport::{Particle, Simulation};
pub const C: f64 = radiation::C;
pub const GEV_J: f64 = 1.602_176_634e-10;
pub const HBAR_GEV_S: f64 = radiation::H / (std::f64::consts::TAU * GEV_J);
pub const HBAR_C: f64 = HBAR_GEV_S * C;
pub const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19;
pub const ALPHA: f64 = 1.0 / 137.035_999_084;
