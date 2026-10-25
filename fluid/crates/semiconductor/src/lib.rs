//! Silicon device simulation: the lower rungs of the electronics ladder
//! (`analysis/emergence_ladder.md`).
//!
//! - `silicon`: material model (band gap, carrier statistics, mobility,
//!   recombination) with sources;
//! - `device`: drift–diffusion TCAD on tensor meshes (Si + SiO₂);
//! - `mosfet`: planar MOSFET structures;
//! - `compact`: table compact models folded from TCAD, with certificates;
//! - `circuit`: CMOS inverter, noise margins, ring oscillator;
//! - `gates`: standard cells (INV, NAND2, NOR2, XOR2, full adder) at
//!   transistor level, and their timing/energy characterisation;
//! - `timing`: NLDM interpolation, event-driven gate simulation, static
//!   timing analysis and the error certificate;
//! - `analytic`: exact references (Poisson–Boltzmann MOS, Pao–Sah,
//!   Shockley);
//! - `band`, `mesh`: banded LU and meshes.

pub mod analytic;
pub mod band;
pub mod circuit;
pub mod compact;
pub mod device;
pub mod gates;
pub mod mesh;
pub mod mosfet;
pub mod silicon;
pub mod quantum;
pub mod seq;
pub mod thermal;
pub mod timing;
pub mod tunnel;
