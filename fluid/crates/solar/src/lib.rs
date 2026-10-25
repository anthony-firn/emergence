//! Reduced solar models and explicit observational forecast validation.
pub mod corona;
pub mod photosphere;
pub mod render;
pub mod structure;
pub mod surface;
pub mod trace;
pub const R_SUN: f64 = 6.957e8;
pub const GM_SUN: f64 = 1.327_124_4e20;
pub const L_SUN: f64 = 3.828e26;
