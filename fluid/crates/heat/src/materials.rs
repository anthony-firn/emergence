//! Material models for lattice heat transport.
//!
//! The silicon parameters are the branch-resolved Debye–Callaway set of
//! **Morelli, Heremans & Slack, Phys. Rev. B 66, 195304 (2002)** as tabulated
//! by **Kazan, Guisbiers, Pereira, Correia, Masri, Bruyant, Volz & Royer,
//! J. Appl. Phys. 107, 083503 (2010)**, Table II ("bulk isotopic silicon"):
//! `v_T = 5840 m/s`, `v_L = 8430 m/s`, `Θ_T = 216 K`, `Θ_L = 566 K`,
//! `γ_T = 0.45`, `γ_L = 1.09`, `M = 4.66e-26 kg`, normal fraction `β = 1.5`.
//! These are **fitted** parameters; agreement with the measured κ(T) is a
//! consistency check, not a prediction. The isotope term is the Klemens
//! point-defect rate with the mass-fluctuation parameter computed from the
//! isotopic abundances (not fitted).

use crate::phonon::Branch;
use crate::scattering::ScatteringParams;

/// Natural silicon isotope abundances (IUPAC/CIAAW): ²⁸Si, ²⁹Si, ³⁰Si.
pub const SI_ABUNDANCES: [(f64, f64); 3] =
    [(28.0, 0.922_23), (29.0, 0.046_85), (30.0, 0.030_92)];

/// The mass-fluctuation parameter `Γ = Σ f_i (1 − m_i/M̄)²` for a given set of
/// `(mass, abundance)` pairs. Abundances are normalised internally.
pub fn mass_fluctuation(isotopes: &[(f64, f64)]) -> f64 {
    let tot: f64 = isotopes.iter().map(|x| x.1).sum();
    let mbar: f64 = isotopes.iter().map(|x| x.0 * x.1).sum::<f64>() / tot;
    isotopes
        .iter()
        .map(|(m, f)| (f / tot) * (1.0 - m / mbar).powi(2))
        .sum()
}

/// Natural-silicon mass fluctuation (≈ 2.009e-4, matching Inyushkin et al.
/// 2018, who quote `g² = 200.7e-6`).
pub fn silicon_gamma_natural() -> f64 {
    mass_fluctuation(&SI_ABUNDANCES)
}

/// Branch-resolved silicon model: two transverse and one longitudinal Debye
/// branch (Morelli et al. 2002).
pub fn silicon_branches() -> Vec<Branch> {
    vec![
        Branch::debye("TA", 5840.0, 216.0, 0.45, 2.0),
        Branch::debye("LA", 8430.0, 566.0, 1.09, 1.0),
    ]
}

/// Same branches with a Born–von Kármán sine dispersion, for comparison.
pub fn silicon_sine_branches() -> Vec<Branch> {
    vec![
        Branch::sine("TA", 5840.0, 216.0, 0.45, 2.0),
        Branch::sine("LA", 8430.0, 566.0, 1.09, 1.0),
    ]
}

/// Natural-silicon scattering parameters. `length` is the Casimir boundary
/// length in metres; use `f64::INFINITY` for a bulk crystal. The reference
/// bulk sample length quoted by Kazan et al. is 14 mm.
pub fn silicon_params(length: f64) -> ScatteringParams {
    ScatteringParams {
        // average atomic mass of natural silicon (kg): 28.0855 u
        mass: 4.66e-26,
        beta: 1.5,
        gamma_iso: silicon_gamma_natural(),
        // volume per atom of the diamond lattice: a³/8, a = 0.5431 nm
        volume: 2.0e-29,
        length,
    }
}

/// Enriched ²⁸Si (99.995%), see Inyushkin et al., J. Appl. Phys. 123, 095104
/// (2018). The residual ²⁹/³⁰ content makes Γ tiny but nonzero.
pub fn silicon_params_enriched(length: f64) -> ScatteringParams {
    let isotopes = [(28.0, 0.999_95), (29.0, 0.000_03), (30.0, 0.000_02)];
    ScatteringParams {
        gamma_iso: mass_fluctuation(&isotopes),
        ..silicon_params(length)
    }
}

/// Germanium, branch-resolved (**estimated** parameters, not a published fit;
/// sound speeds and a Debye temperature of 374 K are the standard values, the
/// branch split is the same shape as silicon). `v_T = 3.30e3`,
/// `v_L = 5.40e3 m/s`, `Θ_T = 110 K`, `Θ_L = 340 K`, `γ_T = 0.6`,
/// `γ_L = 1.2`, `M = 1.205e-25 kg`, `V_0 = 4.5e-29 m³`.
pub fn germanium_branches() -> Vec<Branch> {
    vec![
        Branch::debye("TA", 3300.0, 110.0, 0.6, 2.0),
        Branch::debye("LA", 5400.0, 340.0, 1.2, 1.0),
    ]
}

/// Germanium scattering parameters (**estimate**).
pub fn germanium_params(length: f64) -> ScatteringParams {
    let ge = [
        (69.924, 0.2057),
        (71.922, 0.2745),
        (72.923, 0.0775),
        (73.921, 0.3650),
        (75.921, 0.0773),
    ];
    ScatteringParams {
        mass: 1.205e-25,
        beta: 1.5,
        gamma_iso: mass_fluctuation(&ge),
        volume: 4.5e-29,
        length,
    }
}
