//! Equations of state.
//!
//! Noble–Abel stiffened gas (NASG, Le Métayer & Saurel 2016):
//!   p(ρ, e) = (γ − 1)(e − q)/(v − b) − γ p∞,   v = 1/ρ,
//!   c² = γ (p + p∞) v² / (v − b),
//!   T(p, v) = (p + p∞)(v − b)/((γ − 1) c_v).
//! b = q = 0 is the stiffened gas (liquids at moderate pressure), and with
//! p∞ = 0 as well, the ideal gas. NASG reproduces both liquid and vapour
//! saturation curves of a pure substance with one parameter set per phase,
//! which is what phase-transition (cavitation, boiling) models need.

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Nasg {
    pub gamma: f64,
    pub p_inf: f64,
    /// co-volume, m³/kg
    pub b: f64,
    /// reference energy, J/kg
    pub q: f64,
    /// heat capacity at constant volume, J/(kg K)
    pub cv: f64,
    /// entropy constant q′, J/(kg K) (Gibbs free energy / phase equilibrium)
    pub q_prime: f64,
}

impl Nasg {
    pub const fn ideal(gamma: f64) -> Self {
        Nasg { gamma, p_inf: 0.0, b: 0.0, q: 0.0, cv: 1.0, q_prime: 0.0 }
    }

    pub const fn stiffened(gamma: f64, p_inf: f64) -> Self {
        Nasg { gamma, p_inf, b: 0.0, q: 0.0, cv: 1.0, q_prime: 0.0 }
    }

    /// Pressure from density and specific internal energy.
    #[inline]
    pub fn pressure(&self, rho: f64, e: f64) -> f64 {
        let v = 1.0 / rho;
        (self.gamma - 1.0) * (e - self.q) / (v - self.b) - self.gamma * self.p_inf
    }

    /// Specific internal energy from density and pressure.
    #[inline]
    pub fn energy(&self, rho: f64, p: f64) -> f64 {
        let v = 1.0 / rho;
        (p + self.gamma * self.p_inf) * (v - self.b) / (self.gamma - 1.0) + self.q
    }

    /// Volumetric internal energy ρe from ρ and p.
    #[inline]
    pub fn rho_e(&self, rho: f64, p: f64) -> f64 {
        rho * self.energy(rho, p)
    }

    #[inline]
    pub fn sound_speed_sq(&self, rho: f64, p: f64) -> f64 {
        let v = 1.0 / rho;
        self.gamma * (p + self.p_inf) * v * v / (v - self.b)
    }

    #[inline]
    pub fn temperature(&self, rho: f64, p: f64) -> f64 {
        let v = 1.0 / rho;
        (p + self.p_inf) * (v - self.b) / ((self.gamma - 1.0) * self.cv)
    }

    /// Specific Gibbs free energy g(p, T) (NASG closed form); equal Gibbs
    /// energies of two phases define saturation.
    pub fn gibbs(&self, p: f64, t: f64) -> f64 {
        let (g, pi, cv, b, q, qp) = (self.gamma, self.p_inf, self.cv, self.b, self.q, self.q_prime);
        // log form: T^γ overflows for stiff phases
        (g * cv - qp) * t - cv * t * (g * t.ln() - (g - 1.0) * (p + pi).ln()) + b * p + q
    }

    /// Admissible state: positive density, sound speed real.
    #[inline]
    pub fn admissible(&self, rho: f64, p: f64) -> bool {
        rho > 0.0 && 1.0 / rho > self.b && p + self.p_inf > 0.0
    }
}

/// Reference materials (SI units).
pub mod materials {
    use super::Nasg;
    /// Air / nitrogen as an ideal diatomic gas.
    pub const AIR: Nasg = Nasg { gamma: 1.4, p_inf: 0.0, b: 0.0, q: 0.0, cv: 717.5, q_prime: 0.0 };
    /// Liquid water, stiffened gas commonly used for shock–bubble problems.
    pub const WATER_SG: Nasg = Nasg { gamma: 4.4, p_inf: 6.0e8, b: 0.0, q: 0.0, cv: 1.0, q_prime: 0.0 };
    /// Liquid water, NASG (Le Métayer & Saurel 2016, fitted 300–500 K).
    pub const WATER_LIQUID_NASG: Nasg = Nasg { gamma: 1.19, p_inf: 6.217e8, b: 6.7212e-4, q: -1.177788e6, cv: 3610.0, q_prime: 0.0 };
    /// Water vapour, NASG (same reference).
    pub const WATER_VAPOR_NASG: Nasg = Nasg { gamma: 1.47, p_inf: 0.0, b: 0.0, q: 2.077616e6, cv: 955.0, q_prime: 14317.0 };
}
