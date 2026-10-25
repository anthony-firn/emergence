//! Compressible multi-material finite-volume solver.
//!
//! The incompressible spectral solver (`spectral-core`) is the coarse model
//! where Ma ≪ 1; this crate is the model it unfolds to when sound, shocks
//! or cavitation matter.

pub mod elements;
pub mod allspeed;
pub mod allotropes;
pub mod amr;
pub mod boundary;
pub mod capillary;
pub mod eos;
pub mod exact;
pub mod gpr;
pub mod gravity;
pub mod helmholtz;
pub mod convection2d;
pub mod isobaric;
pub mod isobaric2d;
pub mod kinetics;
pub mod multiphase;
pub mod phase_change;
pub mod realflow;
pub mod realfluid;
pub mod radiation;
pub mod reference_data;
pub mod refrigerant_data;
pub mod scalars;
pub mod substances;
pub mod standard_model;
pub mod scheme;
pub mod sharp_interface;
pub mod sharp_interface2d;
pub mod solver;
pub mod system;
pub mod thermal;
pub mod thermal_eos;
pub mod threephase;
pub mod viscous;

pub use scheme::FluxMode;
pub use solver::{Bc, Mesh, Solver};
pub use system::{Euler, System};
