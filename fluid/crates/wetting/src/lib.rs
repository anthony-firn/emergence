//! Wetting of textured surfaces: capillary statics folded into an effective
//! boundary condition for the flow solvers.
//!
//! This is the wetting instance of the ladder in `analysis/emergence_ladder.md`.
//! The fine model is an explicit drop on a periodic microstructure; the coarse
//! model is an effective contact angle, `cos θ* = r cos θ_Y` (Wenzel) or
//! `cos θ* = φ_s (cos θ_Y + 1) − 1` (Cassie–Baxter). The validity indicator is
//! the drop base radius over the texture period, `a / λ`.
//!
//! Modules:
//!
//! * [`constants`] — sourced water properties (surface tension, density) and
//!   the capillary length.
//! * [`young_laplace`] — axisymmetric sessile and pendant drops from the
//!   arc-length Young–Laplace equations (Bashforth–Adams), plus the exact
//!   zero-gravity spherical cap.
//! * [`capillary`] — exact planar results: Jurin rise, the meniscus at a
//!   vertical wall, and the 2D puddle height `2 l_c sin(θ/2)`.
//! * [`texture`] — pillar/post geometry and the Wenzel / Cassie–Baxter fold,
//!   the Cassie→Wenzel transition, impalement pressure and hysteresis.
//! * [`lattice`] — a direct sharp-interface free-energy minimisation of a 2D
//!   drop on an explicit square-wave substrate (the check of the fold).
//!
//! Benchmarks and the fold/unfold map are in `analysis/wetting.md`.

pub mod capillary;
pub mod constants;
pub mod lattice;
pub mod meniscus;
pub mod texture;
pub mod young_laplace;
