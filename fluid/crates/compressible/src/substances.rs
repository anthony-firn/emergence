//! Molecular substances as three-phase models (solid, liquid, vapour),
//! fitted like the elements (`elements::fit_inputs`): the 1-atm transition
//! temperatures and latent heats exactly, everything else predicted and
//! checked against data that was not used (`examples/substances`).
//!
//! Inputs are explicit per substance because molecular solids break the
//! elements' generic rules (ice has half the heat capacity of liquid water
//! and is 8% less dense than it; CO₂ sublimes at 1 atm). Sources: CRC
//! Handbook values as quoted on the Wikipedia data pages (fetched
//! 2026-10-04), IAPWS for water, Petrenko & Whitworth (Physics of Ice) for
//! ice elasticity; inputs that are estimates are listed in `estimated`.
//! Not represented: water's density maximum at 4 °C (the NASG liquid has
//! one sign of thermal expansion), critical points (separate NASG
//! branches never merge — that is the reference-EOS rung).

use crate::elements::{fit_inputs, FitInputs, Grade, ThreePhase, R_GAS};

#[derive(Clone, Debug)]
pub struct Substance {
    pub name: &'static str,
    pub formula: &'static str,
    pub inputs: FitInputs,
    pub sources: &'static str,
    /// liquid–vapour saturation pressure data (T K, p Pa), not fitted
    pub vapour_pressure: Vec<(f64, f64)>,
    /// solid–vapour (sublimation) pressure data (T K, p Pa), not fitted
    pub sublimation_pressure: Vec<(f64, f64)>,
    /// triple point (T K, p Pa) data
    pub triple: Option<(f64, f64)>,
    /// melting curve T(p) data (p Pa, T K), not fitted
    pub melting_curve: Vec<(f64, f64)>,
}

impl Substance {
    pub fn model(&self) -> ThreePhase {
        fit_inputs(&self.inputs)
    }
}

/// IAPWS R14-08 (2011) melting pressure of ice Ih, 251.165 K ≤ T ≤ 273.16 K.
pub fn iapws_ice_ih_melting_pressure(t: f64) -> f64 {
    let th = t / 273.16;
    let (a, b) = ([0.119539337e7, 0.808183159e5, 0.333826860e4], [0.3e1, 0.2575e2, 0.10375e3]);
    611.657 * (1.0 + (0..3).map(|i| a[i] * (1.0 - th.powf(b[i]))).sum::<f64>())
}

fn per_kg(j_per_mol: f64, m: f64) -> f64 {
    j_per_mol / m
}

pub fn all() -> Vec<Substance> {
    let mut v = Vec::new();
    // ---------------------------------------------------------------- water
    {
        let m = 18.015e-3;
        let cvr: f64 = 27.5 / R_GAS; // c_v of steam at 100 °C (data page)
        v.push(Substance {
            name: "water",
            formula: "H₂O",
            inputs: FitInputs {
                molar_mass: m,
                t_melt: 273.15,
                t_boil: 373.124,
                sublimes: false,
                p_triple: Some(611.657),
                rho_s: 916.7, // ice Ih at 0 °C
                t_ref_s: 273.15,
                alpha_s: 3.0 * 5.3e-5, // ice, 0 °C
                c_bulk_s: (8.8e9f64 / 916.7).sqrt(), // K = 8.8 GPa
                shear_speed: (3.5e9f64 / 916.7).sqrt(), // G = 3.5 GPa
                rho_l: 999.84, // 0 °C
                t_ref_l: 273.15,
                alpha_l: 2.1e-4, // 20 °C value; negative below 4 °C in reality
                c_l: 1403.0, // 0 °C
                cp_s: 2110.0, // ice at 0 °C
                cp_l: per_kg(75.97, m), // 0 °C
                gamma_v: 1.0 + 1.0 / cvr,
                cv_v: cvr * R_GAS / m,
                h_fus: per_kg(6010.0, m),
                h_vap: per_kg(40680.0, m),
                t_crit: Some(647.096),
                grade: Grade::A,
                estimated: vec!["liquid thermal expansion at 20 °C (no density maximum)"],
            },
            sources: "data page: Δh_fus 6.01 kJ/mol, Δh_vap 40.68 kJ/mol (373 K), c_p 75.97 J/(mol K) and c 1403 m/s at 0 °C, c_v steam 27.5 J/(mol K); IAPWS: T_b, triple, critical; ice: Petrenko & Whitworth",
            // IAPWS-95 saturation pressures
            vapour_pressure: vec![(298.15, 3169.9), (323.15, 12352.0), (348.15, 38595.0), (373.15, 101418.0), (423.15, 476160.0)],
            // Murphy & Koop (2005) vapour pressure over ice
            sublimation_pressure: vec![(263.15, 259.9), (253.15, 103.2), (243.15, 38.0), (233.15, 12.84)],
            triple: Some((273.16, 611.657)),
            melting_curve: [272.0, 270.0, 265.0, 260.0, 255.0].iter().map(|&t| (iapws_ice_ih_melting_pressure(t), t)).collect(),
        });
    }
    // -------------------------------------------------------- carbon dioxide
    {
        let m = 44.0095e-3;
        let cvr: f64 = (32.36 - R_GAS) / R_GAS; // c_p 32.36 J/(mol K) at 200 K
        let c_l = 950.0;
        v.push(Substance {
            name: "carbon dioxide",
            formula: "CO₂",
            inputs: FitInputs {
                molar_mass: m,
                t_melt: 216.58, // triple point (melting needs p > p_t)
                t_boil: 194.686, // sublimation at 1 atm
                sublimes: true,
                p_triple: Some(518.5e3),
                rho_s: 1562.0,
                t_ref_s: 194.686,
                alpha_s: 2.0 * 54.6 / m / (c_l / 0.85f64).powi(2),
                c_bulk_s: c_l / 0.85,
                shear_speed: 0.5 * c_l / 0.85,
                rho_l: 1178.0, // at the triple point
                t_ref_l: 216.58,
                alpha_l: 2.0 * per_kg(85.0, m) / (c_l * c_l),
                c_l,
                cp_s: per_kg(54.6, m),
                cp_l: per_kg(85.0, m),
                gamma_v: 1.0 + 1.0 / cvr,
                cv_v: cvr * R_GAS / m,
                h_fus: per_kg(9019.0, m),
                h_vap: per_kg(25230.0, m), // sublimation at 194.7 K
                t_crit: Some(304.18),
                grade: Grade::B,
                estimated: vec!["c_liquid ≈ 950 m/s", "solid sound speeds from c_liquid", "thermal expansion from Grüneisen Γ ≈ 2", "c_p solid 54.6 / liquid 85 J/(mol K)"],
            },
            sources: "data page: triple 216.58 K / 518.5 kPa, Δh_fus 9.019 kJ/mol, critical 304.18 K; CRC: Δh_sub 25.23 kJ/mol, ρ_solid 1562, ρ_liquid(triple) 1178",
            // NIST saturation pressures of the liquid
            vapour_pressure: vec![(230.0, 0.8929e6), (250.0, 1.7850e6), (273.15, 3.4851e6), (293.15, 5.7291e6)],
            sublimation_pressure: vec![],
            triple: Some((216.58, 518.5e3)),
            melting_curve: vec![],
        });
    }
    // ---------------------------------------------------------------- ammonia
    {
        let m = 17.031e-3;
        let cvr: f64 = (35.06 - R_GAS) / R_GAS;
        let c_l = 1700.0;
        let cp_l = per_kg(80.80, m);
        v.push(Substance {
            name: "ammonia",
            formula: "NH₃",
            inputs: FitInputs {
                molar_mass: m,
                t_melt: 195.42,
                t_boil: 239.82,
                sublimes: false,
                p_triple: Some(6.06e3),
                rho_s: 817.0,
                t_ref_s: 195.42,
                alpha_s: 2.0 * 0.6 * cp_l / (c_l / 0.85f64).powi(2),
                c_bulk_s: c_l / 0.85,
                shear_speed: 0.5 * c_l / 0.85,
                rho_l: 682.0, // at the boiling point
                t_ref_l: 239.82,
                alpha_l: 2.0 * cp_l / (c_l * c_l),
                c_l,
                cp_s: 0.6 * cp_l,
                cp_l,
                gamma_v: 1.0 + 1.0 / cvr,
                cv_v: cvr * R_GAS / m,
                h_fus: per_kg(5653.0, m),
                h_vap: per_kg(23350.0, m),
                t_crit: Some(405.5),
                grade: Grade::B,
                estimated: vec!["c_liquid ≈ 1700 m/s", "solid sound speeds from c_liquid", "c_p,solid ≈ 0.6 c_p,liquid", "thermal expansion from Grüneisen Γ ≈ 2"],
            },
            sources: "data page: Δh_fus 5.653, Δh_vap 23.35 kJ/mol, c_p liquid 80.80 / gas 35.06 J/(mol K), triple 195.4 K / 6.06 kPa, critical 405.5 K; CRC: T_m, T_b, ρ",
            vapour_pressure: vec![(273.15, 429.4e3), (293.15, 857.5e3), (313.15, 1555.0e3)],
            sublimation_pressure: vec![],
            triple: Some((195.4, 6.06e3)),
            melting_curve: vec![],
        });
    }
    // ---------------------------------------------------------------- methane
    {
        let m = 16.043e-3;
        let cvr = 3.0; // translation + rotation; vibration frozen near 110 K
        let c_l = 1400.0;
        let cp_l = per_kg(54.6, m);
        v.push(Substance {
            name: "methane",
            formula: "CH₄",
            inputs: FitInputs {
                molar_mass: m,
                t_melt: 90.69,
                t_boil: 111.65,
                sublimes: false,
                p_triple: Some(11.7e3),
                rho_s: 490.0,
                t_ref_s: 90.69,
                alpha_s: 2.0 * 0.8 * cp_l / (c_l / 0.85f64).powi(2),
                c_bulk_s: c_l / 0.85,
                shear_speed: 0.5 * c_l / 0.85,
                rho_l: 451.5, // at the triple point
                t_ref_l: 90.69,
                alpha_l: 2.0 * cp_l / (c_l * c_l),
                c_l,
                cp_s: 0.8 * cp_l,
                cp_l,
                gamma_v: 1.0 + 1.0 / cvr,
                cv_v: cvr * R_GAS / m,
                h_fus: per_kg(1100.0, m),
                h_vap: per_kg(8170.0, m),
                t_crit: Some(190.6),
                grade: Grade::B,
                estimated: vec!["ρ_solid ≈ 490 kg/m³", "c_liquid ≈ 1400 m/s", "c_p liquid ≈ 54.6 J/(mol K), solid ≈ 0.8 liquid", "solid sound speeds from c_liquid", "thermal expansion from Grüneisen Γ ≈ 2"],
            },
            sources: "data page: Δh_fus 1.1, Δh_vap 8.17 kJ/mol, triple 90.67 K / 11.7 kPa, critical 190.6 K; CRC: T_m, T_b, ρ_liquid",
            vapour_pressure: vec![(100.0, 34.4e3), (120.0, 191.4e3), (150.0, 1040.5e3), (180.0, 3286.1e3)],
            sublimation_pressure: vec![],
            triple: Some((90.67, 11.7e3)),
            melting_curve: vec![],
        });
    }
    v
}

pub fn by_name(name: &str) -> Option<Substance> {
    all().into_iter().find(|s| s.name == name || s.formula == name)
}
