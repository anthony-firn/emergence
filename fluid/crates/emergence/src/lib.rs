//! The emergence engine.
//!
//! Every accurate shortcut in this repository has one shape
//! (`analysis/emergence_ladder.md`): a cheap model valid in a regime, a cheap
//! indicator of its error, and an unfold to a finer model (or a closure
//! computed once from it) where the indicator fails. This crate makes that
//! shape executable and shared:
//!
//! * [`Claim`] / [`Rigor`]: every rung states its error *and* how far that
//!   statement can be trusted (asymptotic estimate, validated, rigorous bound);
//! * [`Rung`]: a model with a cost, an a priori indicator and an optional a
//!   posteriori certificate;
//! * [`Ladder`]: picks the cheapest rung whose claim meets a [`Request`],
//!   unfolds when a certificate fails, never reports a tolerance it did not
//!   reach, and audits its own choices against the reference rung, learning a
//!   safety factor where an indicator was too optimistic;
//! * [`Closure`] / [`ClosureRung`]: fold a fine model once into an adaptive
//!   Chebyshev table with a held-out error check, which becomes a rung;
//! * [`Ladder::solve_functional`]: goal-oriented budgets, spending work only
//!   where an item's weight in the final answer makes its error matter.

pub mod budget;
pub mod claim;
pub mod closure;
pub mod ladder;
pub mod rung;

pub use budget::{Functional, Plan};
pub use claim::{Assessment, Claim, Request, Rigor};
pub use closure::{Axis, Closure, ClosureOptions, ClosureRung};
pub use ladder::{Audit, Ladder, RungStats, Solved};
pub use rung::{FnRung, Rung};
