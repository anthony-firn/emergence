//! Heat generation in a cell and the interface to the thermal solvers.
//!
//! Bernardi's lumped energy balance gives the local heat rate per unit
//! electrode area
//!
//! ```text
//!   q = I (U − V) − I T dU/dT
//!       \_______/   \_________/
//!       irreversible  reversible (entropic)
//! ```
//!
//! (`I > 0` discharge, `U` the open-circuit voltage, `V` the terminal voltage,
//! `dU/dT` the entropic coefficient). The first term is the Joule/overpotential
//! loss and is positive for discharge; the second can heat or cool the cell.
//!
//! Source: D. Bernardi, E. Pawlikowski, J. Newman, *J. Electrochem. Soc.*
//! **132**, 5 (1985).
//!
//! # Interface to `heat`
//! The `heat` crate provides lattice conductivity, `κ(T, L)`
//! ([`heat::effective_kappa`]) and its Fourier-validity flag; this crate
//! provides the *source*. A continuum solver assembles
//!
//! ```text
//!   ρ c_p ∂T/∂t = ∇·(κ ∇T) + q_vol,   q_vol = q_areal / L
//! ```
//!
//! where `L` is the electrode thickness ([`volumetric_source`]). The exchange
//! is a heat-rate density and does not require modifying `heat`.

/// Irreversible (overpotential) heat per unit area (W m⁻²): `i (U − V)`.
pub fn irreversible_heat(i_app: f64, u_ocv: f64, v: f64) -> f64 {
    i_app * (u_ocv - v)
}

/// Reversible (entropic) heat per unit area (W m⁻²): `−i T dU/dT`.
pub fn reversible_heat(i_app: f64, t: f64, dudt: f64) -> f64 {
    -i_app * t * dudt
}

/// Total Bernardi heat per unit area (W m⁻²).
pub fn total_heat(i_app: f64, u_ocv: f64, v: f64, t: f64, dudt: f64) -> f64 {
    irreversible_heat(i_app, u_ocv, v) + reversible_heat(i_app, t, dudt)
}

/// Convert an area-specific heat rate to a volumetric source (W m⁻³) for a
/// region of thickness `l` (m).
pub fn volumetric_source(q_areal: f64, l: f64) -> f64 {
    q_areal / l
}

/// Rough entropic coefficient of a graphite/NMC cell (V K⁻¹).
///
/// `ESTIMATE`: sign-changing over state of charge, magnitude ~0.1–0.5 mV/K;
/// the caller should replace it with a measured `dU/dT(z)`.
pub fn dudt_estimate(state_of_charge: f64) -> f64 {
    let z = state_of_charge.clamp(0.0, 1.0);
    -1.0e-4 + 4.0e-4 * (std::f64::consts::PI * z).sin().abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_current_gives_no_heat() {
        assert_eq!(total_heat(0.0, 4.0, 3.7, 298.15, 1e-4), 0.0);
    }

    #[test]
    fn irreversible_heat_is_positive_on_discharge() {
        let q = irreversible_heat(10.0, 4.0, 3.7);
        assert!(q > 0.0);
        assert!((q - 3.0).abs() < 1e-12);
    }

    #[test]
    fn reversible_heat_flips_with_current_direction() {
        let a = reversible_heat(10.0, 298.15, 2e-4);
        let b = reversible_heat(-10.0, 298.15, 2e-4);
        assert!((a + b).abs() < 1e-12);
    }

    #[test]
    fn volumetric_source_conversion() {
        // 100 W/m² over an 80 µm electrode → 1.25e6 W/m³.
        assert!((volumetric_source(100.0, 80e-6) - 1.25e6).abs() < 1.0);
    }
}
