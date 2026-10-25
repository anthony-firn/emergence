//! Fold a Maxwell mode onto the GPR (Godunov–Peshkov–Romenski) unified
//! continuum model's shear speed and strain-relaxation time.
//!
//! In GPR (see `compressible/src/gpr.rs`) the shear (elastic) modulus is
//! `G = ρ c_s²` and the emergent Newtonian viscosity in the `De ≪ 1` limit
//! is
//!
//! ```text
//! μ = ρ c_s² τ / 6,
//! ```
//!
//! where `τ` is the strain-relaxation time (the module's own convention;
//! the factor 6 is the linear-relaxation result: the shear strain decays
//! as `e^{-6t/τ}`). A Maxwell mode has `G(t) = G e^{-t/λ}` and
//! `η = G λ`. Matching the modulus and viscosity gives
//!
//! ```text
//! ρ c_s² = G,        τ_GPR = 6 λ.
//! ```
//!
//! That factor of six is exactly where a naive `τ = λ` mapping would put
//! GPR's emergent viscosity at `η/6`.

/// A single Maxwell mode to be represented in GPR.
#[derive(Clone, Copy, Debug)]
pub struct MaxwellMode {
    /// shear modulus `G = n k_B T`
    pub g: f64,
    /// relaxation time `λ`
    pub lambda: f64,
}

impl MaxwellMode {
    /// Maxwell zero-shear viscosity `η = G λ`.
    pub fn eta(&self) -> f64 {
        self.g * self.lambda
    }
}

/// GPR parameters representing a Maxwell mode.
#[derive(Clone, Copy, Debug)]
pub struct GprRheology {
    pub rho: f64,
    pub cs: f64,
    pub tau: f64,
}

impl GprRheology {
    /// GPR elastic shear modulus `G = ρ c_s²`.
    pub fn shear_modulus(&self) -> f64 {
        self.rho * self.cs * self.cs
    }
    /// GPR emergent viscosity `μ = ρ c_s² τ / 6`.
    pub fn viscosity(&self) -> f64 {
        self.rho * self.cs * self.cs * self.tau / 6.0
    }
    /// The equivalent Maxwell time `λ = τ/6`.
    pub fn lambda(&self) -> f64 {
        self.tau / 6.0
    }
}

/// Map a Maxwell mode onto GPR at a chosen density `ρ`.
pub fn maxwell_to_gpr(m: MaxwellMode, rho: f64) -> GprRheology {
    GprRheology {
        rho,
        cs: (m.g / rho).sqrt(),
        tau: 6.0 * m.lambda,
    }
}

/// Inverse map.
pub fn gpr_to_maxwell(r: GprRheology) -> MaxwellMode {
    MaxwellMode { g: r.shear_modulus(), lambda: r.lambda() }
}

/// Shear strain of a linear GPR cell under no forcing decays as
/// `e^{-6 t/τ}`; the stress `G γ(t)` is returned. This is the scalar
/// linearisation of `Gpr::relax` and is what fixes the factor 6.
pub fn linear_gpr_stress(r: &GprRheology, t: f64) -> f64 {
    r.shear_modulus() * (-6.0 * t / r.tau).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn maxwell_map_preserves_viscosity_and_shear_speed() {
        let m = MaxwellMode { g: 2.5e5, lambda: 1e-3 };
        let r = maxwell_to_gpr(m, 1000.0);
        // G = ρ c_s²
        assert!((r.shear_modulus() - m.g).abs() / m.g < 1e-14);
        // η = G λ = μ = ρ c_s² τ/6
        assert!((r.viscosity() - m.eta()).abs() / m.eta() < 1e-14);
        // τ = 6λ
        assert!((r.tau - 6.0 * m.lambda).abs() < 1e-18);
        // Round trip.
        let back = gpr_to_maxwell(r);
        assert!((back.g - m.g).abs() / m.g < 1e-14);
        assert!((back.lambda - m.lambda).abs() / m.lambda < 1e-14);
    }

    #[test]
    fn linear_gpr_stress_has_maxwell_time() {
        let m = MaxwellMode { g: 1.0, lambda: 0.5 };
        let r = maxwell_to_gpr(m, 1.0);
        // σ(t) = G e^{-t/λ}; at t = λ it has fallen to G/e.
        let s = linear_gpr_stress(&r, m.lambda);
        assert!((s - m.g / std::f64::consts::E).abs() < 1e-12);
    }
}
