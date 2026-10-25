//! Phonon heat transport: kinetic lattice conductivity, size effects, and the
//! diffusive → ballistic fold.
//!
//! This is the phonon rung of the heat ladder in `analysis/emergence_ladder.md`:
//!
//! ```text
//!   continuum heat equation  ρc ∂T/∂t = ∇·(κ(T,L) ∇T)   <-- coarse model
//!            ^   indicator: Kn = Λ/L  (phonon Knudsen number)
//!   kinetic (Callaway/Holland) κ(T), boundary scattering, ballistic Landauer
//!            ^   exact test case: grey BTE between two diffuse walls
//!   discrete ordinates / integral Boltzmann equation
//! ```
//!
//! * [`phonon`] — per-branch Debye and Born–von Kármán sine dispersions, the
//!   mode heat capacity.
//! * [`scattering`] — umklapp, normal, point-defect (isotope) and Casimir
//!   boundary rates.
//! * [`kinetic`] — the relaxation-time conductivity
//!   `κ = (1/3) Σ ∫ C v_g² τ g dω` and the Callaway normal-drift correction.
//! * [`bte1d`] — the grey Boltzmann equation between two walls, solved by
//!   discrete ordinates and by the exact integral equation; the
//!   Fourier↔ballistic fold.
//! * [`materials`] — silicon (natural and ²⁸Si-enriched) and germanium
//!   parameters with sources.
//! * [`effective`] — `κ(T, L)` and the Fourier-validity flag for a continuum
//!   heat solver.
//!
//! # Quick start
//!
//! ```
//! use heat::materials::{silicon_branches, silicon_params};
//! use heat::kinetic::kappa_rta;
//!
//! let branches = silicon_branches();
//! let params = silicon_params(f64::INFINITY);
//! let k = kappa_rta(&branches, &params, 300.0);
//! assert!(k > 50.0 && k < 250.0);
//! ```

pub mod bte1d;
pub mod bte_ladder;
pub mod effective;
pub mod kinetic;
pub mod linalg;
pub mod materials;
pub mod phonon;
pub mod scattering;

use phonon::Branch;
use scattering::ScatteringParams;

/// Reference bulk sample length used by the silicon benchmarks: the 14 mm
/// crystal quoted by Kazan et al. (2010). At low temperature the intrinsic
/// umklapp rate vanishes and only boundary scattering keeps κ finite, so the
/// κ(T) curve requires a finite sample length.
pub const SI_REFERENCE_LENGTH: f64 = 14.0e-3;

/// Natural-silicon lattice conductivity κ(T) (W/(m K)) from the
/// Morelli/Kazan branch model, for a 14 mm reference sample.
pub fn silicon_kappa(t: f64) -> f64 {
    kinetic::kappa_rta(
        &materials::silicon_branches(),
        &materials::silicon_params(SI_REFERENCE_LENGTH),
        t,
    )
}

/// Natural-silicon lattice conductivity including the Callaway normal-drift
/// correction (14 mm reference sample).
pub fn silicon_kappa_callaway(t: f64) -> f64 {
    kinetic::kappa_callaway(
        &materials::silicon_branches(),
        &materials::silicon_params(SI_REFERENCE_LENGTH),
        t,
    )
}

/// Isotopically enriched ²⁸Si lattice conductivity κ(T) (W/(m K))
/// (14 mm reference sample).
pub fn silicon_kappa_enriched(t: f64) -> f64 {
    kinetic::kappa_rta(
        &materials::silicon_branches(),
        &materials::silicon_params_enriched(SI_REFERENCE_LENGTH),
        t,
    )
}

/// Effective conductivity of silicon at temperature `t` in a sample of size
/// `l` (metres), with the Fourier-validity flag.
pub fn silicon_effective_kappa(t: f64, l: f64) -> effective::EffectiveKappa {
    effective::effective_kappa(
        &materials::silicon_branches(),
        &materials::silicon_params(l),
        t,
        l,
    )
}

/// Generic per-branch shim so callers can reuse a branch list with a single
/// parameter set.
pub fn branch_set<'a>(
    branches: &'a [Branch],
    params: &'a ScatteringParams,
) -> (&'a [Branch], &'a ScatteringParams) {
    (branches, params)
}
