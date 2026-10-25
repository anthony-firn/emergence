//! Incompressible Navier–Stokes on a Cartesian MAC grid with immersed CAD
//! geometry and engineering boundary conditions (velocity inlets, pressure
//! outlets, no-slip/moving walls, slip, periodic).

pub mod case;
pub mod field;
pub mod poisson;
pub mod post;
pub mod solver;

pub use field::{Arr, Grid};
pub use solver::{Bc, Body, Solver};
