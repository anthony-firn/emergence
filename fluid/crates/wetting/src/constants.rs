//! Sourced material properties and derived capillary scales.
//!
//! Every constant carries its source. Water at 20 °C:
//! * surface tension γ = 72.8 mN/m — Vargaftik, Volkov & Voljak, *J. Phys.
//!   Chem. Ref. Data* **12**, 817 (1983); also the CRC Handbook.
//! * density ρ = 998.21 kg/m³ — IAPWS-95 / CRC Handbook at 20 °C.
//! * standard gravity g = 9.80665 m/s² (CGPM 1901).

/// Surface tension of water against air at 20 °C (N/m).
pub const WATER_GAMMA: f64 = 72.8e-3;
/// Density of water at 20 °C (kg/m³).
pub const WATER_RHO: f64 = 998.21;
/// Standard gravity (m/s²).
pub const G_STANDARD: f64 = 9.80665;

/// Capillary length `l_c = sqrt(γ/(ρg))` (m).
pub fn capillary_length(gamma: f64, rho: f64, g: f64) -> f64 {
    (gamma / (rho * g)).sqrt()
}

/// Capillary length of water at 20 °C: `l_c = 2.727 mm`.
pub fn water_capillary_length() -> f64 {
    capillary_length(WATER_GAMMA, WATER_RHO, G_STANDARD)
}
