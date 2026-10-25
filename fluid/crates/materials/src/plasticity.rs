//! The fold from the defect rung to continuum plasticity.
//!
//! The coarse model is the GPR elasto-plastic solid of
//! `fluid/crates/compressible/src/gpr.rs` (read, not modified): a
//! hyperbolic solid whose excess shear stress relaxes on a time `τ`, with
//! `TauModel::Plastic { tau0, sigma_y, n }` dropping `τ` above a von Mises
//! yield stress `σ_y`. The fine model is the atomistic γ-surface computed
//! in [`crate::defects`]. The bridge is the Peierls-Nabarro (PN) model:
//! a dislocation core is a soliton of the γ-surface, the single-crystal
//! critical resolved shear stress (CRSS) is its Peierls stress `σ_P`, and
//! a Taylor factor folds it onto the polycrystal yield stress.
//!
//! This module deliberately *does not* depend on the `compressible` crate.
//! It returns a plain [`GprPlastic`] whose fields map one-to-one onto the
//! GPR closure:
//!
//! ```text
//! GprPlastic { tau0, sigma_y, n }  <->  compressible::gpr::TauModel::Plastic { tau0, sigma_y, n }
//! ```
//!
//! ## Peierls-Nabarro model
//!
//! For an edge dislocation with Burgers vector `b` on a slip plane of
//! spacing `d`, the classical core width is
//!
//! ```text
//! ζ = d / (2 (1 - ν)),
//! ```
//!
//! and the Peierls stress is exponentially small in the core width,
//!
//! ```text
//! σ_P ≈ (2 G / (1 - ν)) exp(-2 π ζ / b).
//! ```
//!
//! With a sinusoidal γ-surface `γ(φ) = (γ_us/2)(1 - cos(2πφ/b))` the
//! continuum PN equation,
//!
//! ```text
//! (G / (2π(1-ν))) ∫ (dφ/dx')/(x - x') dx' = dγ/dφ,
//! ```
//!
//! has the exact soliton solution `φ(x) = b/2 + (b/π) arctan((x-x0)/ζ)`,
//! whose core width is fixed by the γ amplitude,
//!
//! ```text
//! ζ = G b² / (4 π² (1 - ν) γ_us).
//! ```
//!
//! The natural elastic amplitude `γ_us = G b² / (2 π² d)` reproduces the
//! classical `ζ = d/(2(1-ν))`. The continuum soliton is *translation
//! invariant* (any applied stress moves it), so the Peierls stress comes
//! from the discreteness of the lattice, not from the continuum soliton;
//! the exponential above is the leading lattice result. PN is therefore a
//! **crude closure**: the exponential can be wrong by orders of magnitude
//! and the prefactor is not universally agreed (compare the fault-plane
//! ideal strength `π γ_us / b`). Reference: Nabarro, *Proc. Phys. Soc.*
//! **59**, 256 (1947); Peierls, *Proc. Phys. Soc.* **52**, 34 (1940);
//! Hirth & Lothe, *Theory of Dislocations*; Vitek, *Phil. Mag.* **18**,
//! 773 (1968); Schoeck, *Mater. Sci. Eng. A* **400-401**, 7 (2005).
//!
//! ## The fold to the GPR yield stress
//!
//! `σ_y = M_T τ_CRSS` with the fcc polycrystal Taylor factor
//! `M_T ≈ 3.06` (Taylor 1938; Kocks 1970 gives 3.06 for random fcc
//! textures). The single-crystal `τ_CRSS` is either the PN stress or the
//! Taylor hardening law `τ = α G b √ρ` (Taylor 1934; `α` typically 0.2-0.5
//! for fcc). Above `σ_y` the GPR shear stress relaxes in `τ0`; matching the
//! Orowan plastic strain rate `γ̇ = ρ b² τ_res / B` (drag coefficient `B`)
//! to the GPR rate `σ/(G τ0)` gives the mapping
//!
//! ```text
//! τ0 = M_T B / (G ρ b²),
//! ```
//!
//! documented in [`gpr_tau0_from_drag`].

use std::f64::consts::PI;

/// Taylor factor for a random fcc polycrystal.
pub const TAYLOR_FCC: f64 = 3.06;

/// Burgers vector of a perfect fcc dislocation, `|a/2⟨110⟩| = a/√2`.
pub fn fcc_burgers(a_cubic: f64) -> f64 {
    a_cubic / 2.0_f64.sqrt()
}

/// `{111}` interplanar spacing, `d = a/√3`.
pub fn fcc_slip_plane_spacing(a_cubic: f64) -> f64 {
    a_cubic / 3.0_f64.sqrt()
}

/// Peierls-Nabarro estimate for an edge dislocation.
#[derive(Clone, Copy, Debug)]
pub struct PeierlsNabarro {
    pub b: f64,
    pub d: f64,
    pub nu: f64,
    pub g: f64,
    /// Classical core width `d / (2(1-ν))`.
    pub zeta: f64,
    /// `(2G/(1-ν)) exp(-2πζ/b)`.
    pub sigma_p: f64,
    /// Elastic/sinusoidal γ amplitude `G b²/(2π²d)` that reproduces `ζ`.
    pub gamma_us_sinusoidal: f64,
}

/// PN core width and stress from the isotropic shear modulus `g` and
/// Poisson ratio `nu`.
pub fn peierls_nabarro(g: f64, nu: f64, b: f64, d: f64) -> PeierlsNabarro {
    let zeta = d / (2.0 * (1.0 - nu));
    let sigma_p = (2.0 * g / (1.0 - nu)) * (-2.0 * PI * zeta / b).exp();
    PeierlsNabarro {
        b,
        d,
        nu,
        g,
        zeta,
        sigma_p,
        gamma_us_sinusoidal: g * b * b / (2.0 * PI * PI * d),
    }
}

/// Core width implied by a measured sinusoidal γ amplitude,
/// `ζ = G b² / (4π²(1-ν)γ_us)`.
pub fn sinusoidal_core_width(g: f64, nu: f64, b: f64, gamma_us: f64) -> f64 {
    g * b * b / (4.0 * PI * PI * (1.0 - nu) * gamma_us)
}

/// Exact continuum PN soliton profile
/// `φ(x) = b/2 + (b/π) arctan((x-x0)/ζ)`.
pub fn sinusoidal_core(x: f64, b: f64, zeta: f64) -> f64 {
    b / 2.0 + (b / PI) * (x / zeta).atan()
}

/// Ideal shear strength of the fault plane, `max dγ/dφ = π γ_us / b`, for
/// comparison with the exponentially reduced `σ_P`.
pub fn fault_plane_strength(gamma_us: f64, b: f64) -> f64 {
    PI * gamma_us / b
}

/// Taylor hardening: `τ = α G b √ρ` (Taylor 1934).
pub fn taylor_hardening(alpha: f64, g: f64, b: f64, rho_dislocation: f64) -> f64 {
    alpha * g * b * rho_dislocation.max(0.0).sqrt()
}

/// Mean dislocation spacing `ℓ = 1/√ρ`.
pub fn dislocation_spacing(rho_dislocation: f64) -> f64 {
    if rho_dislocation > 0.0 {
        1.0 / rho_dislocation.sqrt()
    } else {
        f64::INFINITY
    }
}

/// Polycrystal yield stress from a single-crystal CRSS.
pub fn taylor_yield(crss: f64, taylor_factor: f64) -> f64 {
    taylor_factor * crss
}

/// GPR `TauModel::Plastic` parameters, mirroring the compressible crate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GprPlastic {
    /// Strain-relaxation time at the yield surface.
    pub tau0: f64,
    /// von Mises yield stress.
    pub sigma_y: f64,
    /// Hardening exponent (large = nearly ideal plasticity).
    pub n: f64,
}

/// Map a yield stress (Pa) and a relaxation time (s) onto the GPR closure.
pub fn gpr_plastic(sigma_y: f64, tau0: f64, n: f64) -> GprPlastic {
    GprPlastic { tau0, sigma_y, n }
}

/// Dislocation-drag mapping of the relaxation time,
/// `τ0 = M_T B / (G ρ b²)`.
///
/// Derivation: above yield the GPR deviatoric stress relaxes as
/// `dσ/dt ≈ -σ/τ0`, so the plastic strain rate is `γ̇_GPR = σ/(G τ0)`.
/// Orowan glide with density `ρ`, Burgers vector `b` and drag `B` gives
/// `γ̇_Orowan = ρ b² τ_res / B` with the resolved stress `τ_res = σ/M_T`.
/// Equating the two eliminates `σ` and yields the mapping.
pub fn gpr_tau0_from_drag(drag: f64, g: f64, b: f64, rho_dislocation: f64, taylor_factor: f64) -> f64 {
    taylor_factor * drag / (g * rho_dislocation.max(1e-300) * b * b)
}

/// Which plasticity rung applies. The continuum GPR closure is valid when
/// many dislocations populate the sample **and** the stress is at or above
/// yield. Otherwise the discrete dislocation/atomistic rung must be
/// unfolded.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PlasticityRung {
    /// `σ < σ_P`: dislocations are pinned, the response is elastic.
    Elastic,
    /// Sample size `L ≳ ℓ = 1/√ρ`: a statistical dislocation ensemble gives
    /// the continuum yield law.
    ContinuumPlastic,
    /// `L ≲ ℓ`: too few dislocations; individual cores and source-limited
    /// strengthening ("smaller is stronger") must be resolved.
    Discrete,
}

/// Validity indicator for the continuum-plasticity fold.
///
/// `sample_size` is the smallest dimension over which the stress is
/// carried, `rho_dislocation` the dislocation density, `sigma` the applied
/// von Mises stress and `sigma_p` the Peierls stress. The two boundaries
/// are: (i) `σ < σ_P` — no glide, the material is elastic; (ii)
/// `L < ℓ = 1/√ρ` — fewer than one dislocation per sample dimension, so
/// the continuum density is meaningless and the discrete rung applies.
/// The boundary is order-of-magnitude, not sharp.
pub fn plasticity_rung(
    sample_size: f64,
    rho_dislocation: f64,
    sigma: f64,
    sigma_p: f64,
) -> PlasticityRung {
    if sigma < sigma_p {
        return PlasticityRung::Elastic;
    }
    if sample_size < dislocation_spacing(rho_dislocation) {
        PlasticityRung::Discrete
    } else {
        PlasticityRung::ContinuumPlastic
    }
}

/// Source-limited "smaller is stronger" strength of a dislocation-starved
/// crystal of size `L`: `σ ≈ G b / L` (Orowan/Stokes-Nabarro
/// source-truncation estimate; see e.g. Uchic et al., *Science* **305**,
/// 986 (2004); Greer & Nix, *Phys. Rev. B* **73**, 245410 (2006)). The
/// coefficient is order unity and this is an estimate, not a calibrated
/// law.
pub fn source_limited_strength(g: f64, b: f64, sample_size: f64) -> f64 {
    g * b / sample_size.max(1e-300)
}

/// Hall-Petch polycrystal strength `σ = σ0 + k / √d` (Hall 1951; Petch
/// 1953). `k` is the Hall-Petch coefficient (Pa·m^1/2).
pub fn hall_petch(sigma0: f64, k: f64, grain_size: f64) -> f64 {
    sigma0 + k / grain_size.max(1e-300).sqrt()
}

/// Isotropic bulk and shear moduli from cubic constants via the
/// Voigt-Reuss-Hill average, plus the Poisson ratio.
///
/// Returns `(g_hill, nu)`; the bulk modulus `K=(C11+2C12)/3` is
/// orientation-independent for cubic symmetry.
pub fn isotropic_from_cubic(c11: f64, c12: f64, c44: f64) -> (f64, f64) {
    let k = (c11 + 2.0 * c12) / 3.0;
    let gv = (c11 - c12 + 3.0 * c44) / 5.0;
    let gr = 5.0 * (c11 - c12) * c44 / (4.0 * c44 + 3.0 * (c11 - c12));
    let g = 0.5 * (gv + gr);
    let nu = (3.0 * k - 2.0 * g) / (2.0 * (3.0 * k + g));
    (g, nu)
}
