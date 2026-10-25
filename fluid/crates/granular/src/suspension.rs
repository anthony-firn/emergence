//! Dilute and dense suspension rheology: Einstein, Batchelor–Green and
//! Krieger–Dougherty, with the indicators that say when a Newtonian
//! suspension loses validity.
//!
//! # Einstein's viscosity from the stresslet
//!
//! For `N` non-interacting, force- and torque-free spheres in a volume `V`
//! with an imposed rate of strain `E`, the volume-averaged stress is
//! `sigma = 2 mu E + (N/V) S`, and each sphere contributes the exact
//! single-sphere stresslet `S = (20/3) pi mu a^3 E` ([`crate::stokes`]).
//! With the volume fraction `phi = N (4/3) pi a^3 / V`,
//!
//! ```text
//! sigma = 2 mu E (1 + 5/2 phi),
//! ```
//!
//! i.e. `eta = eta0 (1 + 5/2 phi)` exactly in the dilute limit. The
//! identity `(N/V)(20/3) pi a^3 = (5/2) phi` is checked numerically in
//! [`crate::suspension`] and in the integration tests.
//!
//! # Higher order and dense closures
//!
//! * Batchelor & Green (1972) computed the `O(phi^2)` hydrodynamic
//!   coefficient for a random suspension with **no Brownian motion**:
//!   `k2 = 5.2` (high-frequency viscosity). Batchelor (1977) added the
//!   Brownian stresses and obtained `k2 = 6.17`; the pair-Smoluchowski
//!   result of Cichocki & Felderhof (1988) is `k2 = 5.91`. All three are
//!   provided.
//! * Krieger & Dougherty (1959):
//!   `eta/eta0 = (1 - phi/phi_m)^(-[eta] phi_m)`, with intrinsic viscosity
//!   `[eta] = 5/2` and maximum packing `phi_m`. For hard spheres the
//!   low-shear (random-close-packing) value is `phi_m ~ 0.63`, the
//!   high-shear/fcc value `~ 0.74`; Mendoza & Santamaría-Holek (2009)
//!   report `phi_m = 0.63` fits the low-shear data of de Kruif et al.
//!   (1985) and Krieger (1972). The `O(phi^2)` coefficient implied by the
//!   closure is [`kd_second_coefficient`].
//!
//! Sources: Einstein (1906) *Ann. Phys.* **19**, 289; Batchelor & Green
//! (1972) *J. Fluid Mech.* **56**, 401; Batchelor (1977) *J. Fluid Mech.*
//! **83**, 97; Cichocki & Felderhof (1988) *J. Chem. Phys.* **89**, 3705;
//! Krieger & Dougherty (1959) *Trans. Soc. Rheol.* **3**, 137; Mendoza &
//! Santamaría-Holek (2009) arXiv:0810.5111; de Kruif et al. (1985)
//! *J. Chem. Phys.* **83**, 4717.

use std::f64::consts::PI;

/// Einstein intrinsic viscosity, `[eta] = 5/2`.
pub const INTRINSIC_VISCOSITY: f64 = 2.5;

/// Batchelor & Green (1972) `O(phi^2)` coefficient, non-Brownian
/// (high-frequency) random suspension (confirmed by Mendoza &
/// Santamaría-Holek 2009, Eq. 1).
pub const BG_NONBROWNIAN_K2: f64 = 5.2;

/// Batchelor (1977) `O(phi^2)` coefficient including Brownian stresses
/// (as quoted by Mendoza & Santamaría-Holek 2009, Eq. 2).
pub const BATCHELOR_BROWNIAN_K2: f64 = 6.17;

/// Pair-Smoluchowski (exact to `O(phi^2)`) Brownian coefficient of
/// Cichocki & Felderhof (1988), as quoted by Mendoza & Santamaría-Holek
/// (2009).
pub const CICHOCKI_FELDERHOF_K2: f64 = 5.91;

/// Default maximum packing for the low-shear random-close-packing closure.
pub const PHI_M_RCP: f64 = 0.63;

/// Maximum packing for fcc/colloidal-crystal order (high-shear limit).
pub const PHI_M_FCC: f64 = 0.7404;

/// Einstein dilute-limit relative viscosity `1 + 5/2 phi`.
pub fn einstein(phi: f64) -> f64 {
    1.0 + INTRINSIC_VISCOSITY * phi
}

/// Batchelor–Green second-order relative viscosity, non-Brownian
/// (high-frequency): `1 + 5/2 phi + 5.2 phi^2`.
pub fn batchelor_green(phi: f64) -> f64 {
    1.0 + INTRINSIC_VISCOSITY * phi + BG_NONBROWNIAN_K2 * phi * phi
}

/// Batchelor (1977) second-order relative viscosity including Brownian
/// stresses: `1 + 5/2 phi + 6.17 phi^2`.
pub fn batchelor_brownian(phi: f64) -> f64 {
    1.0 + INTRINSIC_VISCOSITY * phi + BATCHELOR_BROWNIAN_K2 * phi * phi
}

/// Krieger–Dougherty relative viscosity
/// `(1 - phi/phi_m)^(-[eta] phi_m)`.
pub fn krieger_dougherty(phi: f64, phi_m: f64, intrinsic: f64) -> f64 {
    assert!(phi_m > 0.0);
    if phi >= phi_m {
        return f64::INFINITY;
    }
    (1.0 - phi / phi_m).powf(-intrinsic * phi_m)
}

/// The `O(phi^2)` coefficient implied by the Krieger–Dougherty closure,
/// obtained by expanding `(1 - phi/phi_m)^(-[eta] phi_m)`:
/// `k2 = [eta]([eta] phi_m + 1) / (2 phi_m) = 25/8 + [eta]/(2 phi_m)`.
pub fn kd_second_coefficient(phi_m: f64, intrinsic: f64) -> f64 {
    intrinsic * (intrinsic * phi_m + 1.0) / (2.0 * phi_m)
}

/// The divergence exponent of the Krieger–Dougherty closure as
/// `phi -> phi_m`: `eta ~ (1 - phi/phi_m)^(-[eta] phi_m)`.
pub fn kd_divergence_exponent(phi_m: f64, intrinsic: f64) -> f64 {
    intrinsic * phi_m
}

/// Dilute-limit bulk stress from summing `n` per-volume exact single-sphere
/// stresslets. `n` is the number density, `a` the radius and `e` the
/// imposed symmetric traceless rate of strain. Returns `(2 mu E + n S)/(2
/// mu)` as a matrix whose `xy` entry divided by `e_xy` is the relative
/// viscosity. This is the numerical assembly that reproduces [`einstein`].
pub fn dilute_bulk_stress_ratio(
    n: f64,
    a: f64,
    e: [[f64; 3]; 3],
) -> [[f64; 3]; 3] {
    // (N/V) S / (2 mu E) = (n (20/3) pi a^3 / 2) * E / E.
    let coeff = n * (20.0 / 3.0) * PI * a * a * a / 2.0;
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = e[i][j] * (1.0 + coeff);
        }
    }
    r
}

/// Particle Reynolds number `Re_p = rho_f gdot a^2 / eta_f`, the ratio of
/// inertial to viscous stress at the particle scale. Stokes flow (and the
/// Newtonian suspension closure) requires `Re_p << 1`; inertia gives
/// shear thickening and, at larger `Re_p`, turbulence-like effects.
pub fn particle_reynolds(rho_f: f64, gdot: f64, a: f64, eta_f: f64) -> f64 {
    rho_f * gdot * a * a / eta_f
}

/// Péclet number `Pe = 6 pi eta_f a^3 gdot / (k_B T)`, the ratio of the
/// shear rate to the Brownian diffusion rate of a sphere. The Brownian
/// `O(phi^2)` term (Batchelor 1977) is only present for `Pe << 1`.
pub fn peclet(gdot: f64, a: f64, eta_f: f64, kt: f64) -> f64 {
    if kt <= 0.0 {
        return f64::INFINITY;
    }
    6.0 * PI * eta_f * a * a * a * gdot / kt
}

/// Sedimentation Péclet (the same number for the stresslet) is the ratio of
/// the shear to the Brownian stress; the Brownian fold switches at `Pe ~ 1`.
pub fn brownian_stress_fraction(pe: f64) -> f64 {
    1.0 / (1.0 + pe)
}

/// Which regime a suspension is in, from the cheap indicators `phi/phi_m`,
/// `Re_p` and `Pe`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SuspensionRegime {
    /// `phi << phi_m`, `Re_p << 1`, `Pe >> 1`: Einstein/Stokes.
    DiluteNewtonian,
    /// `phi <~ phi_m`, `Re_p << 1`: Batchelor–Green / Krieger–Dougherty.
    DenseNewtonian,
    /// `Pe <~ 1`: Brownian stresses set the viscosity.
    Brownian,
    /// `Re_p >~ 1`: inertial shear thickening; the Stokes closure fails.
    Inertial,
    /// `phi -> phi_m`: jamming / yield stress.
    Jamming,
}

impl SuspensionRegime {
    pub fn label(self) -> &'static str {
        match self {
            SuspensionRegime::DiluteNewtonian => "dilute-Newtonian",
            SuspensionRegime::DenseNewtonian => "dense-Newtonian",
            SuspensionRegime::Brownian => "Brownian",
            SuspensionRegime::Inertial => "inertial shear-thickening",
            SuspensionRegime::Jamming => "jamming",
        }
    }
}

/// Classify a suspension. The thresholds are the standard order-unity
/// crossovers, not sharp constants:
/// `phi/phi_m > 0.95` -> jamming; `Re_p > 1` -> inertial;
/// `Pe < 1` and `Re_p < 1` -> Brownian; otherwise dilute/dense by
/// `phi < 0.1`.
pub fn suspension_regime(
    phi: f64,
    phi_m: f64,
    re_p: f64,
    pe: f64,
) -> SuspensionRegime {
    if phi / phi_m > 0.95 {
        SuspensionRegime::Jamming
    } else if re_p > 1.0 {
        SuspensionRegime::Inertial
    } else if pe < 1.0 {
        SuspensionRegime::Brownian
    } else if phi < 0.1 {
        SuspensionRegime::DiluteNewtonian
    } else {
        SuspensionRegime::DenseNewtonian
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn einstein_stresslet_assembly_matches() {
        // Pick a volume fraction, a radius and a volume; the assembled bulk
        // stress must give exactly 1 + 5/2 phi.
        let phi = 0.08;
        let a = 0.5;
        let v = 1.0;
        // phi = n (4/3) pi a^3 -> n.
        let n = phi / ((4.0 / 3.0) * PI * a * a * a);
        let e = crate::stokes::simple_shear_e(1.0);
        let ratio = dilute_bulk_stress_ratio(n, a, e);
        let eta_rel = ratio[0][1] / e[0][1];
        assert!((eta_rel - einstein(phi)).abs() < 1e-12, "{eta_rel}");
        let _ = v;
    }

    #[test]
    fn kd_second_coefficient_matches_expansion() {
        for phi_m in [0.55, 0.58, 0.63, 0.64, 0.7404] {
            let k2 = kd_second_coefficient(phi_m, INTRINSIC_VISCOSITY);
            // Finite-difference the expansion.
            let f = |x: f64| -> f64 {
                (krieger_dougherty(x, phi_m, INTRINSIC_VISCOSITY)
                    - 1.0
                    - INTRINSIC_VISCOSITY * x)
                    / (x * x)
            };
            let fd = f(1e-4);
            assert!((k2 - fd).abs() < 1e-3 * k2.abs() + 1e-6, "{phi_m}: {k2} vs {fd}");
        }
    }

    #[test]
    fn regimes_switch_where_expected() {
        assert_eq!(
            suspension_regime(0.02, 0.63, 1e-3, 1e4),
            SuspensionRegime::DiluteNewtonian
        );
        assert_eq!(
            suspension_regime(0.4, 0.63, 1e-3, 1e4),
            SuspensionRegime::DenseNewtonian
        );
        assert_eq!(
            suspension_regime(0.4, 0.63, 1e-3, 0.1),
            SuspensionRegime::Brownian
        );
        assert_eq!(
            suspension_regime(0.4, 0.63, 5.0, 1e4),
            SuspensionRegime::Inertial
        );
        assert_eq!(
            suspension_regime(0.62, 0.63, 1e-3, 1e4),
            SuspensionRegime::Jamming
        );
    }
}
