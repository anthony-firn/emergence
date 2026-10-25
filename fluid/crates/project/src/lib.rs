//! The shared project model.
//!
//! One file per product, read and written by every discipline. Parameters,
//! materials, CAD parts, analyses and requirements form one dependency graph,
//! and every value that crosses a team boundary carries
//!
//! * its **dimension**, so units cannot be misread between teams;
//! * an **interval** that contains the true value, propagated from the stated
//!   tolerances and every model's own error claim;
//! * its **rigor** (the emergence engine's [`Rigor`](emergence::Rigor)): a
//!   rigorous bound, a validated one, or an asymptotic estimate;
//! * an **owner**, so a change can be traced to the team that made it and to
//!   the requirements of other teams it moves.
//!
//! A requirement is then *met*, *failed* or honestly *undecided*, and the
//! report names the inputs, and their owners, that its uncertainty comes from.
//! See `analysis/project_model.md`.

pub mod analysis;
pub mod chip_design;
pub mod cli;
pub mod diff;
pub mod doc;
pub mod expr;
pub mod fdm;
pub mod graph;
pub mod html;
pub mod hvac_analyses;
pub mod materials;
pub mod render;
pub mod report;
pub mod serve;
pub mod thermal;
pub mod units;
pub mod value;

pub use doc::Doc;
pub use expr::Expr;
pub use graph::{Check, Evaluation, Project, Status};
pub use materials::Library;
pub use units::{parse_quantity, parse_unit, Dim, Unit};
pub use value::Val;
