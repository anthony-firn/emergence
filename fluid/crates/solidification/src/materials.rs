//! Material properties for the dimensional example.
//!
//! Two classic dendrite materials are provided: succinonitrile (SCN), the
//! transparent plastic crystal used for essentially all quantitative dendritic
//! growth experiments, and ice/water. Every constant carries its source; the
//! ones that are estimates are labelled.
//!
//! The phase-field benchmarks are all *nondimensional*: they use `d0 = 1`,
//! `D = 2` (Karma–Rappel units) and an anisotropy `eps4`, so they do not depend
//! on these numbers. This module only converts a result to SI (microns,
//! kelvin, m/s) for the worked example.

/// Solid–liquid properties of a pure substance, with provenance.
#[derive(Clone, Copy, Debug)]
pub struct MeltMaterial {
    pub name: &'static str,
    /// Melting temperature at one atmosphere (K).
    pub tm: f64,
    /// Latent heat of fusion per unit volume of the solid (J/m³).
    pub latent: f64,
    /// Solid volumetric heat capacity `rho_s c_p` (J/(m³ K)).
    pub c_vol_solid: f64,
    /// Thermal diffusivity of the liquid (m²/s).
    pub alpha_liquid: f64,
    /// Solid–liquid interfacial free energy (J/m²).
    pub gamma: f64,
    /// Capillary length `d0 = gamma Tm c_p / L^2` (m).
    pub d0: f64,
    /// Four-fold surface-energy anisotropy `eps4` (dimensionless).
    pub eps4: f64,
}

impl MeltMaterial {
    /// Gibbs–Thomson coefficient `Gamma = gamma Tm / L` (K m).
    pub fn gibbs_thomson(&self) -> f64 {
        self.gamma * self.tm / self.latent
    }
}

/// Succinonitrile (NC–CH₂–CH₂–CN), the standard dendrite material.
///
/// * `Tm = 331.24 K`, `L = 4.6e4 J/kg`, `c_p = 2.0e3 J/(kg K)`,
///   `rho_l = 9.85e2 kg/m³`, `gamma = 8.9e-3 J/m²`,
///   `alpha_l = 1.1e-7 m²/s` — Glicksman, *Freezing of a pure melt*,
///   and Glicksman, Winsa, Hahn, Lograsso & Tirmizi, "Dendritic growth of
///   succinonitrile", in *Dendritic Growth* / NASA experiments (1976–1993).
///   The values are the ones tabulated in Glicksman's SCN reviews; `d0` and
///   `Gamma` are derived, not fitted.
/// * `eps4 = 0.0055` — measured surface-energy anisotropy of SCN
///   (Glicksman & Singh), *not* the `0.05` used in the benchmark.
pub const SUCCINONITRILE: MeltMaterial = MeltMaterial {
    name: "succinonitrile",
    tm: 331.24,
    latent: 4.55e7,
    c_vol_solid: 1.97e6,
    alpha_liquid: 1.1e-7,
    gamma: 8.9e-3,
    d0: 2.82e-9,
    eps4: 0.0055,
};

/// Ice growing from supercooled water.
///
/// * `Tm = 273.15 K`, `L = 3.34e8 J/m³` (334 kJ/kg of water),
///   `gamma = 2.91e-2 J/m²` (ice/water, 0 °C; Hardy 1977; IAPWS),
///   `c_p(solid) = 2.1e3 J/(kg K)`, `rho_s = 917 kg/m³` (so
///   `c_vol = 1.93e6 J/(m³ K)`), `alpha_solid = 1.14e-6 m²/s`
///   (k = 2.2 W/(m K)).
/// * `d0 = gamma Tm c_p / L^2 ≈ 1.4e-10 m` is derived. It is small because
///   the volumetric latent heat of water is large; ice dendrites are the
///   textbook case where the microstructure is far finer than any
///   macroscopic cell.
/// * `eps4` for ice: the basal/prism anisotropy is six-fold (~0.01–0.05
///   depending on the face). `0.05` is an *estimate* used here because this
///   crate models four-fold anisotropy; see `analysis/solidification.md`.
pub const ICE: MeltMaterial = MeltMaterial {
    name: "ice/water",
    tm: 273.15,
    latent: 3.34e8,
    c_vol_solid: 1.93e6,
    alpha_liquid: 1.3e-7,
    gamma: 2.91e-2,
    d0: 1.38e-10,
    eps4: 0.05,
};

/// Nondimensional temperature `u = (T - Tm) / (L/c_p)` from a physical
/// temperature.
pub fn temperature_to_u(m: &MeltMaterial, t: f64) -> f64 {
    (t - m.tm) / (m.latent / m.c_vol_solid)
}

/// Physical temperature from the nondimensional `u`.
pub fn u_to_temperature(m: &MeltMaterial, u: f64) -> f64 {
    m.tm + u * (m.latent / m.c_vol_solid)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derived_lengths_are_physical() {
        let s = SUCCINONITRILE;
        let gt = s.gibbs_thomson();
        assert!(gt > 1e-8 && gt < 1e-6, "SCN Gamma = {gt}");
        assert!(s.d0 > 1e-9 && s.d0 < 1e-8, "SCN d0 = {}", s.d0);
        // d0 = Gamma / (L/c_p)
        assert!((gt / (s.latent / s.c_vol_solid) - s.d0).abs() / s.d0 < 1e-2);

        assert!(ICE.d0 > 1e-11 && ICE.d0 < 1e-9, "ice d0 = {}", ICE.d0);
    }

    #[test]
    fn round_trip_temperature() {
        let m = SUCCINONITRILE;
        for &t in &[300.0, 331.24, 350.0] {
            let u = temperature_to_u(&m, t);
            assert!((u_to_temperature(&m, u) - t).abs() < 1e-12);
        }
    }
}
