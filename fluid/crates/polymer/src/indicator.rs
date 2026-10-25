//! The viscoelastic validity indicator: Deborah and Weissenberg numbers.
//!
//! A flow sees a Newtonian fluid, a viscoelastic fluid or an elastic solid
//! according to how the material relaxation time `λ` compares with the
//! flow's observation time `t_obs` and its strain rate `γ̇`:
//!
//! * `De = λ / t_obs` (or `λ ω` in oscillatory flow);
//! * `Wi = λ γ̇`.
//!
//! For steady shear with `t_obs = 1/γ̇` the two coincide. The coarse model
//! (Navier–Stokes, no memory) is valid for `De ≪ 1`; a linear viscoelastic
//! constitutive law is valid for small `Wi`; a solid-like elastic response
//! appears for `De ≫ 1`.
//!
//! The chain-level validity indicator is separate and richer: a *Hookean*
//! dumbbell (Oldroyd-B) has no steady state in strong extension once
//! `Wi = λ ε̇ > ½` (the coil–stretch transition), so beyond that the coarse
//! Hookean closure must unfold to a finite-extensibility (FENE) spring.

use crate::fene::HOOKEAN_COIL_STRETCH_WI;

/// Which constitutive rung a flow sees.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RheologyRegime {
    /// `De ≪ 1`: memory is irrelevant, use Navier–Stokes with `η₀`.
    Newtonian,
    /// `De ~ 1` and small `Wi`: linear viscoelastic (Maxwell / Prony).
    Viscoelastic,
    /// `De ≫ 1`: the material has not relaxed; use an elastic/viscoelastic
    /// model with a finite modulus (GPR in the solid limit).
    ElasticSolid,
}

/// The indicator and the regime it selects.
#[derive(Clone, Copy, Debug)]
pub struct ViscoelasticIndicator {
    /// Deborah number `λ / t_obs`
    pub de: f64,
    /// Weissenberg number `λ γ̇`
    pub wi: f64,
    pub regime: RheologyRegime,
}

impl ViscoelasticIndicator {
    /// `de = λ/t_obs`; below `de_lo` the flow is Newtonian, above `de_hi`
    /// it is solid-like.
    pub fn new(lambda: f64, t_obs: f64, gamma_dot: f64) -> Self {
        let de = lambda / t_obs;
        let wi = lambda * gamma_dot;
        let regime = if de < 0.1 {
            RheologyRegime::Newtonian
        } else if de > 10.0 {
            RheologyRegime::ElasticSolid
        } else {
            RheologyRegime::Viscoelastic
        };
        Self { de, wi, regime }
    }

    /// Oscillatory variant: `De = λ ω`.
    pub fn oscillatory(lambda: f64, omega: f64) -> Self {
        Self::new(lambda, 1.0 / omega, omega)
    }

    /// Does a Hookean-dumbbell (Oldroyd-B) closure remain valid in
    /// extension? It fails beyond the coil–stretch transition.
    pub fn hookean_extension_valid(&self) -> bool {
        self.wi <= HOOKEAN_COIL_STRETCH_WI
    }

    /// The extension rate at which a Hookean dumbbell of time `λ` blows up.
    pub fn coil_stretch_rate(lambda: f64) -> f64 {
        HOOKEAN_COIL_STRETCH_WI / lambda
    }
}

/// Steady-shear regime with observation time `1/γ̇`.
pub fn shear_indicator(lambda: f64, gamma_dot: f64) -> ViscoelasticIndicator {
    ViscoelasticIndicator::new(lambda, 1.0 / gamma_dot, gamma_dot)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regimes_and_coil_stretch() {
        assert_eq!(shear_indicator(1.0, 1e-3).regime, RheologyRegime::Newtonian);
        assert_eq!(shear_indicator(1.0, 1.0).regime, RheologyRegime::Viscoelastic);
        assert_eq!(shear_indicator(1.0, 1e3).regime, RheologyRegime::ElasticSolid);

        let ind = ViscoelasticIndicator::new(1.0, 1.0 / 0.6, 0.6);
        assert!(!ind.hookean_extension_valid());
        assert!((ViscoelasticIndicator::coil_stretch_rate(1.0) - 0.5).abs() < 1e-14);
        // De and Wi coincide when t_obs = 1/γ̇.
        assert!((ind.de - ind.wi).abs() < 1e-14);
    }
}
