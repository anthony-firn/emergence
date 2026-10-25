//! f64 CPU reference solver for incompressible Navier–Stokes on the periodic
//! box. It is the accuracy baseline the GPU backend is checked against.

pub mod backend;
pub mod certify;
pub mod constants;
pub mod diagnostics;
pub mod epidemic;
pub mod fft;
pub mod grid;
pub mod imex;
pub mod init;
pub mod krylov;
pub mod ladder;
pub mod mhd;
pub mod nufft;
pub mod penalization;
pub mod regime;
pub mod solver;
pub mod surface;
pub mod twophase;
pub mod virus;

pub use backend::Backend;
pub use diagnostics::Diagnostics;
pub use grid::Grid;
pub use penalization::{Penalization, Shape};
pub use solver::Solver;
pub use surface::Surface;
