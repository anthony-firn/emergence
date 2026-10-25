//! SI-unit radiation transport and a spectral camera sharing optical data.
//! Geometrical optics, unpolarized light, static media (v/c negligible).
pub mod geometry;
pub mod math;
pub mod medium;
pub mod optics;
pub mod render;
pub mod sampling;
pub mod thermal;
pub mod transport;

pub const C: f64 = 299_792_458.0;
pub const H: f64 = 6.626_070_15e-34;
pub const KB: f64 = 1.380_649e-23;
pub const SIGMA_SB: f64 = 5.670_374_419e-8;
pub const A_RAD: f64 = 4.0 * SIGMA_SB / C;
