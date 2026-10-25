//! Effective conductivity and the fold/unfold validity flag for a continuum
//! heat solver.
//!
//! A finite sample of size `L` cannot be described by the bulk conductivity
//! once the phonon mean free path `Λ` approaches `L`. The cheap indicator is
//! the phonon Knudsen number `Kn = Λ/L`:
//!
//! | `Kn` | regime | model |
//! |---|---|---|
//! | `≪ 1` | Fourier | bulk `κ(T)`, boundary scattering folded into `τ` |
//! | `~ 1` | transitional | effective `κ(T, L)`, flag as uncertain |
//! | `≫ 1` | ballistic | Landauer/Casimir `κ ≈ (1/4) C v L` |
//!
//! [`effective_kappa`] returns both the bulk and size-limited conductivity,
//! `Kn`, and the validity flag. A continuum solver that only needs a
//! coefficient for `∇·(κ∇T)` should call this and use
//! [`EffectiveKappa::fourier_valid`] to decide whether to trust it.

use crate::kinetic::{kappa_rta, mean_free_path};
use crate::phonon::Branch;
use crate::scattering::ScatteringParams;

/// Transport regime from the phonon Knudsen number.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TransportRegime {
    Fourier,
    Transitional,
    Ballistic,
}

impl TransportRegime {
    pub fn name(&self) -> &'static str {
        match self {
            TransportRegime::Fourier => "Fourier",
            TransportRegime::Transitional => "transitional",
            TransportRegime::Ballistic => "ballistic",
        }
    }
}

/// Effective conductivity of a finite sample and its validity.
#[derive(Clone, Copy, Debug)]
pub struct EffectiveKappa {
    /// Bulk (infinite) conductivity (W/(m K)).
    pub kappa_bulk: f64,
    /// Size-limited conductivity (W/(m K)), Casimir boundary scattering.
    pub kappa_effective: f64,
    /// Heat-capacity-weighted mean free path `Λ` (m).
    pub mean_free_path: f64,
    /// Phonon Knudsen number `Λ/L`.
    pub knudsen: f64,
    pub regime: TransportRegime,
    /// True when Fourier's law with the returned `kappa_effective` is
    /// trustworthy (`Kn < 0.1`); false when the ballistic unfold is needed.
    pub fourier_valid: bool,
}

impl EffectiveKappa {
    /// The coefficient to hand to a continuum heat equation `ρ c ∂T/∂t =
    /// ∇·(κ∇T)`. Equal to [`Self::kappa_effective`]; the flag says whether the
    /// diffusion closure is valid.
    pub fn conductivity(&self) -> f64 {
        self.kappa_effective
    }
}

/// Effective conductivity `κ(T, L)` of a sample of size `L`, with the
/// boundary scattering length set to `L`.
///
/// The bulk value is computed first (Casimir length `= ∞`); the size-limited
/// value adds `τ_B^{-1} = v/L` to the combined rate. The ratio
/// `κ_eff/κ_bulk` is the fold error of the Fourier closure.
pub fn effective_kappa(
    branches: &[Branch],
    base: &ScatteringParams,
    t: f64,
    l: f64,
) -> EffectiveKappa {
    let bulk = ScatteringParams {
        length: f64::INFINITY,
        ..*base
    };
    let bounded = ScatteringParams { length: l, ..*base };
    let kappa_bulk = kappa_rta(branches, &bulk, t);
    let kappa_effective = kappa_rta(branches, &bounded, t);
    let (mfp, _) = mean_free_path(branches, &bulk, t);
    let kn = mfp / l;
    let regime = if kn < 0.1 {
        TransportRegime::Fourier
    } else if kn > 10.0 {
        TransportRegime::Ballistic
    } else {
        TransportRegime::Transitional
    };
    EffectiveKappa {
        kappa_bulk,
        kappa_effective,
        mean_free_path: mfp,
        knudsen: kn,
        regime,
        fourier_valid: kn < 0.1,
    }
}

/// Landauer/Casimir ballistic limit of a grey gas,
/// `κ_ballistic = (1/4) C v L` (W/(m K)).
pub fn kappa_ballistic(cv: f64, v: f64, l: f64) -> f64 {
    0.25 * cv * v * l
}
