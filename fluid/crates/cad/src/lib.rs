//! Parametric CAD for the simulation crates.
//!
//! The solvers in this workspace all need the same thing from geometry and
//! each needs it in a different form: the flow solver wants a signed
//! distance and a wall normal, the drift–diffusion solver wants doped
//! regions on a tensor mesh, the optics crates want a layer stack, the
//! lattice crate wants a strut graph. This crate is the one place that
//! *makes* geometry, and it hands each solver the form it takes.
//!
//! ```text
//!   parameters (expr.rs)
//!        |
//!   sketches (sketch.rs) --extrude/revolve--+
//!   NURBS curves/surfaces (nurbs.rs) -------+
//!   paths (path.rs) ------sweep/loft--------+--> Solid  (solid.rs)
//!   imported CAD (geometry::Sdf) -----------+     | implicit algebra:
//!                                                 | booleans, offsets,
//!   generators (kits/): parts, lattices,          | shells, arrays,
//!   pcb, silicon, textile, pipe, section,         | polar folds,
//!   airfoil, optic, crystal, sheet  ------------->| periodic level sets
//!                                                 v
//!                      meshing.rs: watertight triangles, anisotropic
//!                      cells, block culling
//!                                                 |
//!    sim.rs: signed-distance grids, material label grids, layer stacks,
//!            strut graphs, doping closures, VTK, STL, case files
//!   laplace.rs: 2D electrostatics on a drawn cross-section
//! ```
//!
//! # The kernel invariant
//!
//! A [`Solid`] is a scalar field, not a boundary representation. Every node
//! of the tree keeps `sdf` 1-Lipschitz and zero exactly on the boundary, so
//! the **sign is exact everywhere** and the magnitude is a **safe lower
//! bound** on the clearance. [`Solid::exact`] says when the bound is tight.
//! [`solid`] explains why this is enough to model without a B-rep kernel,
//! and what it costs.
//!
//! # What is validated
//!
//! Every generator is checked against an exact result or a published one;
//! the numbers and the misses are in `analysis/cad.md`. In short: exact
//! primitive distances to 1e-15, Pappus and Green identities to 1e-12,
//! second-order mesh convergence, ISO 724 thread dimensions, the standard
//! involute gear relations, Hammerstad–Jensen microstrip impedance from a
//! field solve on the generated cross-section, a MOSFET built by a process
//! flow reproducing the validated device in `semiconductor`, Peirce's woven
//! geometry and Munden's knitted relations, and octet-truss stiffness
//! agreeing with the homogenisation in `architected`.
//!
//! # Honest limits
//!
//! * There is **no topological B-rep kernel**. [`nurbs`] has rational
//!   curves and single untrimmed surface patches, but no trimmed faces, no
//!   surface-surface intersection, no exact fillet topology and no feature
//!   history referring to faces and edges — so no STEP export. Exact CAD
//!   still comes in through `geometry::cad` (STEP/IGES by tessellation)
//!   and goes out as triangles.
//! * **Sharp edges are rounded by the grid** when a solid is meshed: the
//!   dual isosurface cannot represent a crease sharper than one cell.
//!   Analytic queries (`sdf`) keep the sharp edge exactly.
//! * **No constraint solver.** Sketches are explicit; `expr` gives
//!   dimension-driven parametrics, not geometric constraints.
//! * **No assemblies, mates, tolerances or drawings.**
//! * **Manufacturing output is unverified** against a CAM tool: the Gerber
//!   and Excellon writers in [`kits::pcb`] are written to their
//!   specifications and structurally self-tested, and nothing more.
//! * The generators are **structure generators, not process simulators**.
//!   Nothing relaxes a fabric or a nanoparticle, predicts a K-factor, or
//!   chooses a glass; each kit says which of its inputs are the caller's
//!   measurements.

pub mod expr;
pub mod laplace;
pub mod meshing;
pub mod model;
pub mod nurbs;
pub mod path;
pub mod sketch;
pub mod solid;

pub mod kits;
pub mod sim;

pub use expr::Params;
pub use meshing::{mesh_solid, mesh_solid_in, surface_nets, Grid3};
pub use nurbs::{Nurbs, Surface};
pub use path::{Path3, Profile};
pub use sketch::{Polygon, Region2, V2};
pub use solid::{Placement, Solid, Tpms};

/// Re-exported so callers need not depend on `geometry` directly.
pub use geometry::{P3, TriMesh};
