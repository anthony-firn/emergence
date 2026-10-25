//! The fold to the isotropic continuum: effective polycrystal moduli and the
//! grain-size validity indicator.
//!
//! A single crystal is anisotropic. A random-orientation polycrystal is
//! (statistically) isotropic at scales much larger than a grain, so the
//! directional rung of [`crate::anisotropic`] folds onto two moduli, a bulk
//! modulus `K` and a shear modulus `G`. This module provides:
//!
//! - the Voigt (uniform strain), Reuss (uniform stress) and Hill averages;
//! - the Hashin-Shtrikman bounds for a cubic polycrystal;
//! - the isotropic moduli, Young's modulus, Poisson ratio and sound speeds,
//!   and the parameters the GPR solid model needs (see the note below);
//! - the emergent boundary between the isotropic and single-crystal rungs
//!   as a grain-size/wavelength indicator, and the leading-order Rayleigh
//!   longitudinal attenuation.
//!
//! The averages are exact orientation averages of the stiffness and
//! compliance. The Hashin-Shtrikman treatment uses the orientation average
//! of a cubic tensor with the two shear modes (`C' = (C11−C12)/2`, weight
//! 2/5, and `C44`, weight 3/5) as the two phases of the standard two-phase
//! bound (Watt, Davies & O'Connell 1976; Peselnick & Meister 1975). For a
//! cubic crystal the bulk modulus is exactly `(C11+2C12)/3`, so only the
//! shear bound is nontrivial.
//!
//! Using the result in the compressible solver (do not add a dependency on
//! that crate here): build a NASG with `p_inf = K/γ`, `b = 0`, `q = 0` so
//! that at `p = 0` its sound speed is `c0 = sqrt(γ p_inf/ρ0) = sqrt(K/ρ0)`,
//! set `Gpr { eos, cs: sqrt(G/ρ0), tau: f64::INFINITY, tau_model:
//! TauModel::Constant }` for an elastic solid. The longitudinal speed the
//! solver will then carry is `sqrt(c0² + 4cs²/3) = sqrt((K+4G/3)/ρ0)`.

use crate::anisotropic::{compliance, M6};
use std::f64::consts::PI;

/// Isotropic bulk and shear moduli, in Pa.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsotropicModuli {
    pub bulk: f64,
    pub shear: f64,
}

impl IsotropicModuli {
    pub fn young(&self) -> f64 {
        9.0 * self.bulk * self.shear / (3.0 * self.bulk + self.shear)
    }
    pub fn poisson(&self) -> f64 {
        (3.0 * self.bulk - 2.0 * self.shear) / (2.0 * (3.0 * self.bulk + self.shear))
    }
    /// Longitudinal (P) sound speed `sqrt((K + 4G/3)/ρ)`.
    pub fn longitudinal_speed(&self, rho: f64) -> f64 {
        ((self.bulk + 4.0 * self.shear / 3.0) / rho).sqrt()
    }
    /// Shear (S) sound speed `sqrt(G/ρ)`.
    pub fn shear_speed(&self, rho: f64) -> f64 {
        (self.shear / rho).sqrt()
    }
    /// Bulk (hydrostatic) sound speed `sqrt(K/ρ)`.
    pub fn bulk_speed(&self, rho: f64) -> f64 {
        (self.bulk / rho).sqrt()
    }
}

/// Voigt (uniform-strain) average of the 6x6 stiffness.
pub fn voigt(c: &M6) -> IsotropicModuli {
    let k = (c[0][0] + c[1][1] + c[2][2] + 2.0 * (c[0][1] + c[0][2] + c[1][2])) / 9.0;
    let g = (c[0][0] + c[1][1] + c[2][2] - (c[0][1] + c[0][2] + c[1][2])
        + 3.0 * (c[3][3] + c[4][4] + c[5][5]))
        / 15.0;
    IsotropicModuli { bulk: k, shear: g }
}

/// Reuss (uniform-stress) average, from the compliance.
pub fn reuss(c: &M6) -> IsotropicModuli {
    let s = compliance(c);
    let k = 1.0 / (s[0][0] + s[1][1] + s[2][2] + 2.0 * (s[0][1] + s[0][2] + s[1][2]));
    let g = 15.0
        / (4.0 * (s[0][0] + s[1][1] + s[2][2]) - 4.0 * (s[0][1] + s[0][2] + s[1][2])
            + 3.0 * (s[3][3] + s[4][4] + s[5][5]));
    IsotropicModuli { bulk: k, shear: g }
}

/// Hill average: the arithmetic mean of the Voigt and Reuss estimates.
pub fn hill(c: &M6) -> IsotropicModuli {
    let v = voigt(c);
    let r = reuss(c);
    IsotropicModuli {
        bulk: 0.5 * (v.bulk + r.bulk),
        shear: 0.5 * (v.shear + r.shear),
    }
}

/// Hashin-Shtrikman shear bounds for a cubic polycrystal (Pa).
///
/// Returns `(G_lower, G_upper)`. The bulk modulus is exactly
/// `(C11+2C12)/3` for a cubic crystal. The two shear modes of the cubic
/// tensor are the soft and stiff phases of the standard two-phase bound.
pub fn hashin_shtrikman_cubic(c11: f64, c12: f64, c44: f64) -> (f64, f64) {
    let k = (c11 + 2.0 * c12) / 3.0;
    let cprime = 0.5 * (c11 - c12);
    let modes = [(cprime, 2.0 / 5.0), (c44, 3.0 / 5.0)];
    let (soft, f_soft) = if modes[0].0 <= modes[1].0 {
        modes[0]
    } else {
        modes[1]
    };
    let (stiff, f_stiff) = if modes[0].0 <= modes[1].0 {
        modes[1]
    } else {
        modes[0]
    };
    let lower = soft
        + f_stiff
            / (1.0 / (stiff - soft)
                + 6.0 * f_soft * (k + 2.0 * soft) / (5.0 * soft * (3.0 * k + 4.0 * soft)));
    let upper = stiff
        + f_soft
            / (1.0 / (soft - stiff)
                + 6.0 * f_stiff * (k + 2.0 * stiff) / (5.0 * stiff * (3.0 * k + 4.0 * stiff)));
    (lower, upper)
}

/// Everything the isotropic continuum needs for a polycrystal of density
/// `rho`.
#[derive(Clone, Copy, Debug)]
pub struct IsotropicParameters {
    pub density: f64,
    pub moduli: IsotropicModuli,
    /// Shear sound speed `sqrt(G/ρ)`, the GPR `cs`.
    pub shear_sound_speed: f64,
    /// Bulk sound speed `sqrt(K/ρ)`, the GPR equation-of-state `c0`.
    pub bulk_sound_speed: f64,
    /// Longitudinal sound speed `sqrt((K + 4G/3)/ρ)`.
    pub longitudinal_speed: f64,
}

/// Isotropic parameters from the Hill average of a stiffness (Hill is the
/// usual best single estimate; Voigt and Reuss bracket it).
pub fn isotropic_parameters_stiffness(c: &M6, density: f64) -> IsotropicParameters {
    isotropic_parameters(hill(c), density)
}

/// Isotropic parameters from moduli already chosen (e.g. an HS bound).
pub fn isotropic_parameters(moduli: IsotropicModuli, density: f64) -> IsotropicParameters {
    IsotropicParameters {
        density,
        moduli,
        shear_sound_speed: moduli.shear_speed(density),
        bulk_sound_speed: moduli.bulk_speed(density),
        longitudinal_speed: moduli.longitudinal_speed(density),
    }
}

// ---------------------------------------------------------------------------
// Emergent boundary: grain size vs wavelength
// ---------------------------------------------------------------------------

/// Which scattering regime the grain-size/wavelength ratio selects.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ScatteringRegime {
    /// `d ≪ λ`: Rayleigh scattering, attenuation `∝ d³ f⁴`; the isotropic
    /// continuum with a scattering loss is valid.
    Rayleigh,
    /// `d ~ λ`: stochastic scattering; neither simple limit holds, and the
    /// isotropic constant-moduli model is marginal.
    Stochastic,
    /// `d ≫ λ`: geometric (ray) regime; the wave sees individual grains, so
    /// the single-crystal anisotropy must be unfolded grain by grain.
    Geometric,
}

/// Validity indicator from grain size `d` and wavelength `λ`.
///
/// The thresholds `d/λ = 0.1` and `d/λ = 10` are conventional (Stanke &
/// Kino 1984; Weaver 1990). They are order-of-magnitude boundaries, not
/// sharp transitions: the Rayleigh asymptotics apply for `d/λ ≪ 1`, the ray
/// limit for `d/λ ≫ 1`.
pub fn scattering_regime(grain_size: f64, wavelength: f64) -> ScatteringRegime {
    let ratio = grain_size / wavelength;
    if ratio < 0.1 {
        ScatteringRegime::Rayleigh
    } else if ratio > 10.0 {
        ScatteringRegime::Geometric
    } else {
        ScatteringRegime::Stochastic
    }
}

/// Whether the isotropic polycrystal model is valid at `(d, λ)`: true in the
/// Rayleigh regime, where the only missing physics is a scattering
/// attenuation, and false in the stochastic and geometric regimes.
pub fn isotropic_polycrystal_valid(grain_size: f64, wavelength: f64) -> bool {
    scattering_regime(grain_size, wavelength) == ScatteringRegime::Rayleigh
}

/// Leading-order Rayleigh longitudinal attenuation `α_L` (1/m) for a cubic
/// polycrystal at frequency `f`, grain size `d` and density `rho`.
///
/// Uses the Stanke-Kino / Weaver `α ∝ d³ f⁴ / v_L⁴` scaling with the
/// orientation variance of the cubic longitudinal modulus as the
/// dimensionless anisotropy factor:
/// `α_L = (8π³/375) (d³ f⁴ / v_L⁴) · Var[M(n)]/⟨M⟩²`, where
/// `M(n) = C_{ijkl} n_i n_j n_k n_l` and, for cubic symmetry,
/// `Var[M]/⟨M⟩² = (16/525)(C11−C12−2C44)² / ⟨M⟩²`. This is the leading
/// Born-approximation result; the overall prefactor is approximate to an
/// order-unity factor and the formula is valid only for `d ≪ λ` (the
/// Rayleigh regime). It is not a substitute for the full Stanke-Kino
/// angular integrals.
pub fn rayleigh_longitudinal_attenuation(
    c11: f64,
    c12: f64,
    c44: f64,
    rho: f64,
    grain_size: f64,
    frequency: f64,
) -> f64 {
    let moduli = hill(&crate::anisotropic::cubic(c11, c12, c44));
    let v_l = moduli.longitudinal_speed(rho);
    let delta = c11 - c12 - 2.0 * c44;
    let mean_m = (3.0 * c11 + 2.0 * c12 + 4.0 * c44) / 5.0;
    let anisotropy = (16.0 / 525.0) * delta * delta / (mean_m * mean_m);
    (8.0 * PI.powi(3) / 375.0) * (grain_size.powi(3) * frequency.powi(4) / v_l.powi(4)) * anisotropy
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::anisotropic::{cubic, hexagonal};

    fn close(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol * b.abs().max(1.0), "{a} vs {b}");
    }

    #[test]
    fn cubic_voigt_reuss_have_the_closed_forms() {
        let (c11, c12, c44) = (168.4e9, 121.4e9, 75.4e9);
        let c = cubic(c11, c12, c44);
        let v = voigt(&c);
        close(v.bulk, (c11 + 2.0 * c12) / 3.0, 1e-12);
        close(v.shear, (c11 - c12 + 3.0 * c44) / 5.0, 1e-12);
        let r = reuss(&c);
        close(r.bulk, (c11 + 2.0 * c12) / 3.0, 1e-12);
        close(
            r.shear,
            5.0 * (c11 - c12) * c44 / (4.0 * c44 + 3.0 * (c11 - c12)),
            1e-12,
        );
        // Cubic bulk is exactly orientation-independent.
        close(v.bulk, r.bulk, 1e-12);
    }

    #[test]
    fn hs_bounds_bracket_hill_and_are_tighter_than_voigt_reuss() {
        let (c11, c12, c44) = (168.4e9, 121.4e9, 75.4e9);
        let c = cubic(c11, c12, c44);
        let (lo, hi) = hashin_shtrikman_cubic(c11, c12, c44);
        let h = hill(&c);
        let v = voigt(&c);
        let r = reuss(&c);
        assert!(lo <= h.shear && h.shear <= hi, "HS must bracket Hill");
        assert!(v.shear >= hi - 1e-6, "HS inside Voigt");
        assert!(r.shear <= lo + 1e-6, "HS inside Reuss");
    }

    #[test]
    fn isotropic_wave_speeds_are_ordered() {
        let m = IsotropicModuli {
            bulk: 137.1e9,
            shear: 47.3e9,
        };
        let rho = 8960.0;
        assert!(m.longitudinal_speed(rho) > m.shear_speed(rho));
        assert!(m.bulk_speed(rho) > m.shear_speed(rho));
    }

    #[test]
    fn scattering_indicator_boundaries() {
        assert_eq!(scattering_regime(1e-5, 1e-3), ScatteringRegime::Rayleigh);
        assert_eq!(scattering_regime(1e-3, 1e-3), ScatteringRegime::Stochastic);
        assert_eq!(scattering_regime(1e-1, 1e-3), ScatteringRegime::Geometric);
        assert!(isotropic_polycrystal_valid(1e-5, 1e-3));
        assert!(!isotropic_polycrystal_valid(1e-1, 1e-3));
    }

    #[test]
    fn attenuation_scales_as_rayleigh() {
        let (c11, c12, c44) = (168.4e9, 121.4e9, 75.4e9);
        let rho = 8960.0;
        let a1 = rayleigh_longitudinal_attenuation(c11, c12, c44, rho, 1e-4, 5e6);
        let a2 = rayleigh_longitudinal_attenuation(c11, c12, c44, rho, 2e-4, 5e6);
        let a3 = rayleigh_longitudinal_attenuation(c11, c12, c44, rho, 1e-4, 1e7);
        close(a2 / a1, 8.0, 1e-12); // d^3
        close(a3 / a1, 16.0, 1e-12); // f^4
    }

    #[test]
    fn hexagonal_voigt_bulk_uses_the_five_constants() {
        let c = hexagonal(59.7e9, 25.8e9, 21.7e9, 61.7e9, 16.4e9);
        let v = voigt(&c);
        let expected = (2.0 * 59.7 + 61.7 + 2.0 * (25.8 + 2.0 * 21.7)) / 9.0 * 1e9;
        close(v.bulk, expected, 1e-12);
    }
}
