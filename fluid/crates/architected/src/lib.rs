//! Architected mechanical metamaterials: the fold from a beam lattice to an
//! effective elastic continuum, and the unfold back to the struts where the
//! fold fails.
//!
//! This is the mechanical instance of the ladder in
//! `analysis/emergence_ladder.md`. The coarse model is an effective stiffness
//! tensor `C*`; the fine model is a linear frame (or pin-jointed bar) finite
//! element lattice. The fold is periodic homogenisation of one unit cell:
//!
//! ```text
//!   finite N x N lattice of beams  --homogenise-->  C* (continuum)
//!          ^   indicator: cell size / deformation wavelength, and strut stress
//!          |             / Euler buckling stress
//!      unfolds to the discrete lattice when either indicator fails
//! ```
//!
//! Modules:
//!
//! * [`beam`] — 2D/3D linear frame and bar element stiffness matrices.
//! * [`cell`] — periodic unit cells and effective stiffness by the energy
//!   method, for frame or truss lattices.
//! * [`lattices`] — honeycomb, re-entrant honeycomb, triangular, octet truss,
//!   simple-cubic and square beam lattices.
//! * [`frame`] — a small direct frame solver for cantilever and fixed-fixed
//!   benchmark beams.
//! * [`fold`] — finite `N x N` lattice blocks compared with the homogenised
//!   continuum, the size-effect indicator and the boundary-layer correction.
//! * [`buckling`] — Euler strut buckling inside the macroscopic state: where
//!   the linear homogenised model stops holding.
//! * [`bloch`] — elastic-wave dispersion of a 2D beam lattice (Bloch-Floquet),
//!   whose long-wavelength slopes must equal the homogenised stiffness.
//! * [`linalg`] — the small dense solver used throughout.
//!
//! The benchmarks and the fold/unfold map are in `analysis/architected.md`.

pub mod beam;
pub mod bloch;
pub mod buckling;
pub mod cell;
pub mod fold;
pub mod frame;
pub mod lattices;
pub mod linalg;

pub use cell::{Cell2d, Cell3d, Kind, Moduli2d, Moduli3d};
