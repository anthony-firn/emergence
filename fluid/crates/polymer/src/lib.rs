//! Polymers and viscoelasticity: "gloopiness" from chains to continuum.
//!
//! This crate is the viscoelastic rung of
//! [`emergence_ladder.md`](../../../analysis/emergence_ladder.md): the fine
//! model is a bead-spring chain (a Rouse or Zimm chain, or a single
//! dumbbell) and the coarse model is a linear/nonlinear viscoelastic
//! constitutive law (Maxwell, standard linear solid, Oldroyd-B). The fold
//! between them is the set of normal-mode relaxation times and the stress
//! relaxation modulus `G(t)`, and the validity indicator is the Deborah
//! (or Weissenberg) number.
//!
//! Rungs
//! -----
//! * [`rouse`] — Rouse normal modes, `τ_p`, `G(t)`, `η₀`, `G*(ω)`, and a
//!   fold of the Rouse spectrum onto a few-mode Maxwell / Prony model with
//!   a stated error;
//! * [`zimm`] — the same with a pre-averaged hydrodynamic interaction
//!   (Zimm), solved as a generalized eigenproblem;
//! * [`bd`] — Brownian dynamics of a Rouse chain (Euler–Maruyama) whose
//!   Green–Kubo `G(t)` and end-to-end autocorrelation are checked against
//!   the analytic modes;
//! * [`fene`] — finite extensibility: Hookean dumbbells blow up when the
//!   Weissenberg number passes ½, FENE springs do not;
//! * [`continuum`] — Maxwell, standard linear solid, multi-mode Maxwell
//!   (Prony), and the Oldroyd-B fluid (the exact continuum fold of Hookean
//!   dumbbells), with exact start-up / oscillatory benchmarks;
//! * [`rheology`] — a virtual DMA / rheometer and time–temperature
//!   superposition with the WLF equation;
//! * [`indicator`] — Deborah and Weissenberg numbers classify a flow as
//!   Newtonian, viscoelastic or elastic, and switch the chain model;
//! * [`gpr_fold`] — map a Maxwell mode onto the GPR strain-relaxation time
//!   and shear speed (with GPR's `μ = ρc_s²τ/6` convention).
//!
//! All constants are cited in the module that uses them. See
//! `analysis/polymer.md` for the benchmarks, sources and errors.

pub mod bd;
pub mod continuum;
pub mod fene;
pub mod gpr_fold;
pub mod indicator;
pub mod rheology;
pub mod rng;
pub mod rouse;
pub mod zimm;

pub use rng::Rng;
