//! The periodic table as phase models.
//!
//! Every element with enough data gets three phases — solid, liquid,
//! vapour — each a stiffened-gas/NASG equation of state, fitted at 1 atm
//! so that the melting and boiling points and both latent heats are
//! reproduced exactly. Everything else is a prediction of the model and is
//! checked against independent data (`examples/periodic_table`):
//! the standard state at 298 K, triple points, and the full vapour-pressure
//! curve (1 Pa … 100 kPa). Critical points are outside what separate-phase
//! NASG models can describe (the liquid and vapour branches never merge);
//! that rung needs a reference equation of state.
//!
//! Data: `data/elements.csv` (built by `fluid/data/build_elements.py`,
//! every value with its source). Missing values are filled by classic
//! correlations and the element's grade records how much was estimated:
//!   A  all inputs measured
//!   B  latent heats or condensed-phase sound speeds from correlations
//!      (Richards' rule ΔS_fus ≈ R, Trouton's rule ΔS_vap ≈ 10.5R,
//!      c_liquid ≈ 0.85 c_bulk,solid, densities ±5–10% between phases)
//!   C  several generic estimates (expect tens of % in derived slopes)
//!   P  melting/boiling point or other inputs are published predictions
//!      or estimates, not measurements (late actinides, superheavy
//!      elements; `data/elements_predicted.csv`, built by
//!      `fluid/data/build_predicted.py` from the quoted predictions)
//!   D  not modelled (no melting/boiling data, measured or predicted)

use crate::eos::Nasg;
use std::sync::OnceLock;

pub const R_GAS: f64 = 8.314_462_618;
pub const P_ATM: f64 = 101_325.0;

#[derive(Clone, Debug)]
pub struct Element {
    pub z: u32,
    pub symbol: String,
    pub name: String,
    pub mass_u: f64,
    pub standard_state: String,
    pub block: String,
    pub radioactive: bool,
    pub rho_solid: Option<f64>,
    pub rho_liquid_mp: Option<f64>,
    pub t_melt: Option<f64>,
    pub t_boil: Option<f64>,
    pub sublimes: bool,
    pub allotrope: String,
    pub h_fus: Option<f64>,
    pub h_vap: Option<f64>,
    pub cp_molar: Option<f64>,
    pub t_crit: Option<f64>,
    pub p_crit: Option<f64>,
    pub t_triple: Option<f64>,
    pub p_triple: Option<f64>,
    pub bulk_modulus: Option<f64>,
    pub shear_modulus: Option<f64>,
    pub c_long_solid: Option<f64>,
    pub c_trans_solid: Option<f64>,
    pub c_liquid: Option<f64>,
    pub c_liquid_t: Option<f64>,
    /// linear thermal expansion coefficient of the solid at 25 °C, 1/K
    pub alpha_lin_solid: Option<f64>,
    /// disagreements between sources found while building the table
    pub data_flags: String,
    pub sources: String,
    /// inputs taken from published predictions/estimates (field:qualifier[refs])
    pub predicted: String,
}

const CSV: &str = include_str!("../data/elements.csv");
const VP_CSV: &str = include_str!("../data/vapor_pressure.csv");
const PRED_CSV: &str = include_str!("../data/elements_predicted.csv");

/// Split one CSV line (fields may be double-quoted).
fn split_csv(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut q = false;
    for ch in line.chars() {
        match ch {
            '"' => q = !q,
            ',' if !q => out.push(std::mem::take(&mut cur)),
            _ => cur.push(ch),
        }
    }
    out.push(cur);
    out
}

/// All 118 elements, by atomic number (index z − 1).
pub fn all() -> &'static [Element] {
    static T: OnceLock<Vec<Element>> = OnceLock::new();
    T.get_or_init(|| {
        let mut lines = CSV.lines();
        let head: Vec<String> = split_csv(lines.next().unwrap());
        let col = |name: &str| head.iter().position(|h| h == name).unwrap_or_else(|| panic!("column {name}"));
        let c: std::collections::HashMap<&str, usize> = [
            "z", "symbol", "name", "mass_u", "standard_state", "block", "radioactive", "rho_solid", "rho_liquid_mp", "t_melt", "t_boil", "sublimes", "allotrope",
            "h_fus", "h_vap", "cp_molar", "t_crit", "p_crit", "t_triple", "p_triple", "bulk_modulus", "shear_modulus", "c_long_solid", "c_trans_solid", "c_liquid",
            "c_liquid_t", "alpha_lin_solid", "data_flags", "sources",
        ]
        .iter()
        .map(|&n| (n, col(n)))
        .collect();
        lines
            .filter(|l| !l.trim().is_empty())
            .map(|l| {
                let f = split_csv(l);
                let s = |k: &str| f[c[k]].clone();
                let o = |k: &str| f[c[k]].trim().parse::<f64>().ok();
                Element {
                    z: s("z").parse().unwrap(),
                    symbol: s("symbol"),
                    name: s("name"),
                    mass_u: s("mass_u").trim_matches(|ch: char| ch == '[' || ch == ']').parse().unwrap_or(f64::NAN),
                    standard_state: s("standard_state"),
                    block: s("block"),
                    radioactive: s("radioactive") == "1",
                    rho_solid: o("rho_solid"),
                    rho_liquid_mp: o("rho_liquid_mp"),
                    t_melt: o("t_melt"),
                    t_boil: o("t_boil"),
                    sublimes: s("sublimes") == "1",
                    allotrope: s("allotrope"),
                    h_fus: o("h_fus"),
                    h_vap: o("h_vap"),
                    cp_molar: o("cp_molar"),
                    t_crit: o("t_crit"),
                    p_crit: o("p_crit"),
                    t_triple: o("t_triple"),
                    p_triple: o("p_triple"),
                    bulk_modulus: o("bulk_modulus"),
                    shear_modulus: o("shear_modulus"),
                    c_long_solid: o("c_long_solid"),
                    c_trans_solid: o("c_trans_solid"),
                    c_liquid: o("c_liquid"),
                    c_liquid_t: o("c_liquid_t"),
                    alpha_lin_solid: o("alpha_lin_solid"),
                    data_flags: s("data_flags"),
                    sources: s("sources"),
                    predicted: String::new(),
                }
            })
            .collect::<Vec<Element>>()
            .into_iter()
            .map(merge_predicted)
            .collect()
    })
}

/// Fill missing inputs from the predicted/estimated table (measured values
/// always win).
fn merge_predicted(mut el: Element) -> Element {
    let mut lines = PRED_CSV.lines();
    let head: Vec<String> = split_csv(lines.next().unwrap());
    let col = |name: &str| head.iter().position(|h| h == name).unwrap();
    for l in lines.filter(|l| !l.trim().is_empty()) {
        let f = split_csv(l);
        if f[col("z")].parse::<u32>().ok() != Some(el.z) {
            continue;
        }
        let o = |k: &str| f[col(k)].trim().parse::<f64>().ok();
        let fill = |dst: &mut Option<f64>, k: &str| {
            if dst.is_none() {
                *dst = o(k);
            }
        };
        fill(&mut el.t_melt, "t_melt");
        fill(&mut el.t_boil, "t_boil");
        fill(&mut el.rho_solid, "rho_solid");
        fill(&mut el.rho_liquid_mp, "rho_liquid_mp");
        fill(&mut el.h_fus, "h_fus");
        fill(&mut el.h_vap, "h_vap");
        fill(&mut el.cp_molar, "cp_molar");
        el.predicted = f[col("notes")].clone();
    }
    el
}

pub fn by_symbol(sym: &str) -> Option<&'static Element> {
    all().iter().find(|e| e.symbol == sym)
}

/// Measured vapour-pressure curve (p in Pa, T in K) of element z.
pub fn vapor_pressure(z: u32) -> Vec<(f64, f64)> {
    VP_CSV
        .lines()
        .skip(1)
        .filter_map(|l| {
            let f: Vec<&str> = l.split(',').collect();
            (f[0].parse::<u32>().ok()? == z).then(|| (f[1].parse().unwrap(), f[2].parse().unwrap()))
        })
        .collect()
}

/// Elements whose tabulated latent heats refer to the diatomic molecule X₂
/// (both sources quote H₂, N₂, O₂, F₂, Cl₂, Br₂, I₂).
fn atomicity(z: u32) -> f64 {
    if [1, 7, 8, 9, 17, 35, 53].contains(&z) {
        2.0
    } else {
        1.0
    }
}

/// Dominant vapour molecule near the boiling point (atoms per molecule)
/// and its molar c_v / R: translation + rotation, plus vibration mostly
/// excited (90%) for the heavy polyatomic vapours; H₂ … I₂ keep 5/2
/// (vibration frozen at their low boiling points; Te₂ is hot and heavy).
fn vapour_species(z: u32) -> (f64, f64) {
    match z {
        1 | 7 | 8 | 9 | 17 | 35 | 53 => (2.0, 2.5),
        // P₄ and As₄: the tabulated heats (per atom) give ΔS_vap ≈ 11R per
        // mole of P₄ (Trouton-consistent), so the P₄ vapour is used
        15 | 33 => (4.0, 3.0 + 0.9 * 6.0),
        // S vapour is an S₈ ⇄ S₆ ⇄ S₂ equilibrium (speciation is a chemistry
        // sub-model of its own); S₈ is used and the element is flagged
        16 => (8.0, 3.0 + 0.9 * 18.0),
        // Se, Sb, Te: the per-atom heats are Trouton-consistent only with
        // a monatomic vapour basis (ΔS/R ≈ 11–12.5 per atom)
        _ => (1.0, 1.5),
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Grade {
    A,
    B,
    C,
    P,
    D,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Phase {
    Solid,
    Liquid,
    Vapor,
}

/// Three-phase model of one element (per kg of the vapour's formula unit).
#[derive(Clone, Debug)]
pub struct ThreePhase {
    pub solid: Nasg,
    pub liquid: Nasg,
    pub vapor: Nasg,
    /// shear wave speed of the solid (for the GPR solid model)
    pub shear_speed: f64,
    /// kg per mole of formula units (X₂ for diatomic elements)
    pub molar_mass: f64,
    pub grade: Grade,
    /// inputs that were estimated
    pub estimated: Vec<&'static str>,
    /// critical temperature (data) above which only one fluid phase exists;
    /// separate liquid/vapour NASG branches cannot represent it, so it is
    /// imposed
    pub t_crit: Option<f64>,
}

/// NASG through (ρ, c, c_p, α) at (T, p) (Le Métayer & Saurel 2016):
///   γ = 1 + c²α²T/c_p,  c_v = c_p/γ,  p∞ = (γ−1)c_v/(αv) − p,  b = v(1 − αT).
/// With b = 0 the same three data force α = 1/T (a solid would double its
/// volume from 300 to 600 K), so the covolume carries the thermal
/// expansion.
fn nasg_fit(rho: f64, c: f64, cp: f64, alpha: f64, t: f64, p: f64) -> Nasg {
    let v = 1.0 / rho;
    let gamma = 1.0 + c * c * alpha * alpha * t / cp;
    let cv = cp / gamma;
    let p_inf = (gamma - 1.0) * cv / (alpha * v) - p;
    Nasg { gamma, p_inf, b: v * (1.0 - alpha * t), q: 0.0, cv, q_prime: 0.0 }
}

/// Volumetric thermal expansion from the thermodynamic Grüneisen relation
/// α = Γ c_p / c² with Γ ≈ 2 (typical of solids and liquid metals; ±30%).
fn alpha_gruneisen(cp: f64, c: f64) -> f64 {
    2.0 * cp / (c * c)
}

/// Gibbs energy without the q′ term: g = a(p, T) − q′ T.
fn gibbs_a(e: &Nasg, p: f64, t: f64) -> f64 {
    e.gibbs(p, t) + e.q_prime * t
}

pub fn fit(el: &Element) -> Result<ThreePhase, String> {
    let (Some(tm), Some(tb)) = (el.t_melt, el.t_boil) else {
        return Err("no melting or boiling point".into());
    };
    if !el.mass_u.is_finite() {
        return Err("no atomic mass".into());
    }
    let mut est: Vec<&'static str> = Vec::new();
    let nat = atomicity(el.z);
    let m = nat * el.mass_u * 1e-3; // kg per mole of formula units
    let gas_at_rt = el.standard_state.contains("Gas");
    // vapour: ideal gas of its dominant molecule
    let (n_mol, cvr) = vapour_species(el.z);
    let gamma_v = 1.0 + 1.0 / cvr;
    let cv_v = cvr * R_GAS / (n_mol * el.mass_u * 1e-3);
    let cp_v = gamma_v * cv_v;
    // condensed-phase heat capacities (per kg)
    let cp_rt = el.cp_molar.map(|c| c * nat / m); // CRC molar value: per mole of the standard-state unit
    let (cp_l, cp_s) = if gas_at_rt {
        est.push("cp_liquid ≈ 2 cp_gas, cp_solid ≈ 0.8 cp_liquid");
        (2.0 * cp_v, 1.6 * cp_v)
    } else if el.standard_state.contains("Liquid") {
        let c = cp_rt.ok_or("no heat capacity")?;
        est.push("cp_solid ≈ 0.9 cp_liquid");
        (c, 0.9 * c)
    } else {
        let c = match cp_rt {
            Some(c) => c,
            None => {
                est.push("cp_solid = 3R (Dulong–Petit)");
                3.0 * R_GAS / (el.mass_u * 1e-3)
            }
        };
        est.push("cp_liquid ≈ 1.15 cp_solid");
        (1.15 * c, c)
    };
    // densities
    let rho_l = match (el.rho_liquid_mp, el.rho_solid) {
        (Some(r), _) => r,
        (None, Some(rs)) => {
            est.push("ρ_liquid ≈ 0.95 ρ_solid");
            0.95 * rs
        }
        _ => return Err("no condensed-phase density".into()),
    };
    let rho_s = match el.rho_solid {
        Some(r) => r,
        None => {
            est.push("ρ_solid ≈ 1.1 ρ_liquid");
            1.1 * rho_l
        }
    };
    // solid bulk and shear sound speeds
    let c_bulk = if let Some(k) = el.bulk_modulus {
        (k / rho_s).sqrt()
    } else if let (Some(cl), Some(ct)) = (el.c_long_solid, el.c_trans_solid) {
        (cl * cl - 4.0 / 3.0 * ct * ct).max(0.1 * cl * cl).sqrt()
    } else if let Some(cl) = el.c_liquid {
        est.push("c_bulk,solid ≈ c_liquid / 0.85");
        cl / 0.85
    } else {
        est.push("c_bulk,solid generic (√(15 R T_m / M))");
        (15.0 * R_GAS * tm / (el.mass_u * 1e-3)).sqrt()
    };
    let shear_speed = if let Some(g) = el.shear_modulus {
        (g / rho_s).sqrt()
    } else if let Some(ct) = el.c_trans_solid {
        ct
    } else {
        est.push("c_shear ≈ 0.5 c_bulk");
        0.5 * c_bulk
    };
    let c_l = match el.c_liquid {
        Some(c) => c,
        None => {
            est.push("c_liquid ≈ 0.85 c_bulk,solid");
            0.85 * c_bulk
        }
    };
    // latent heats (per kg)
    let h_fus = match el.h_fus {
        Some(h) => h / m,
        None => {
            est.push("Δh_fus from Richards' rule (ΔS ≈ R per atom)");
            nat * R_GAS * tm / m
        }
    };
    let h_vap = match el.h_vap {
        Some(h) => h / m,
        None => {
            est.push("Δh_vap from Trouton's rule (ΔS ≈ 10.5 R)");
            10.5 * R_GAS * tb / m
        }
    };
    let t_ref_s = if tm < 298.15 { tm } else { 298.15 };
    let alpha_s = match el.alpha_lin_solid {
        Some(a) if t_ref_s > 200.0 => 3.0 * a,
        _ => {
            est.push("thermal expansion of the solid from Grüneisen Γ ≈ 2");
            alpha_gruneisen(cp_s, c_bulk)
        }
    };
    est.push("thermal expansion of the liquid from Grüneisen Γ ≈ 2");
    // the liquid density is at the melting point, or at the boiling point
    // for the elements that are gases at room temperature (infobox data)
    let t_ref_l = if el.sources.contains("rho_liquid:infobox(bp)") { tb } else { tm };
    if el.sublimes {
        est.push("sublimation: Δh_sub = listed vaporization heat; liquid placed at the triple point");
    }
    let grade = if est.iter().any(|e| e.contains("generic") || e.contains("sublimation")) || est.len() > 5 {
        Grade::C
    } else if est.iter().any(|e| !e.starts_with("cp_liquid") && !e.starts_with("cp_solid") && !e.starts_with("thermal expansion of the liquid")) {
        Grade::B
    } else {
        Grade::A
    };
    if el.z == 16 {
        est.push("vapour speciation (S₈/S₆/S₂) not modelled");
    }
    if [1, 2].contains(&el.z) {
        est.push("quantum liquid: classical EOS outside its range");
    }
    let grade = if el.z == 16 || el.z == 1 || el.z == 2 { Grade::C } else { grade };
    let grade = if el.predicted.is_empty() {
        grade
    } else {
        est.push("melting/boiling point or other inputs are published predictions (grade P)");
        Grade::P
    };
    Ok(fit_inputs(&FitInputs {
        molar_mass: m,
        t_melt: tm,
        t_boil: tb,
        sublimes: el.sublimes,
        p_triple: el.p_triple,
        rho_s,
        t_ref_s,
        alpha_s: alpha_s.min(0.5 / t_ref_s),
        c_bulk_s: c_bulk,
        shear_speed,
        rho_l,
        t_ref_l,
        alpha_l: alpha_gruneisen(cp_l, c_l).min(0.5 / t_ref_l),
        c_l,
        cp_s,
        cp_l,
        gamma_v,
        cv_v,
        h_fus,
        h_vap,
        t_crit: el.t_crit,
        grade,
        estimated: est,
    }))
}

/// Explicit inputs of a three-phase fit (all per kg of the formula unit).
#[derive(Clone, Debug)]
pub struct FitInputs {
    /// kg per mole of formula units
    pub molar_mass: f64,
    /// melting (or triple, for sublimers) and boiling (or sublimation)
    /// temperatures at 1 atm
    pub t_melt: f64,
    pub t_boil: f64,
    pub sublimes: bool,
    pub p_triple: Option<f64>,
    /// solid: density at t_ref_s, volumetric expansion, bulk and shear
    /// sound speeds
    pub rho_s: f64,
    pub t_ref_s: f64,
    pub alpha_s: f64,
    pub c_bulk_s: f64,
    pub shear_speed: f64,
    /// liquid: density at t_ref_l, volumetric expansion, sound speed
    pub rho_l: f64,
    pub t_ref_l: f64,
    pub alpha_l: f64,
    pub c_l: f64,
    /// heat capacities (J/kg/K) and the vapour's ideal-gas γ, c_v
    pub cp_s: f64,
    pub cp_l: f64,
    pub gamma_v: f64,
    pub cv_v: f64,
    /// latent heats (J/kg): fusion, and vaporization (sublimation for
    /// sublimers) at 1 atm
    pub h_fus: f64,
    pub h_vap: f64,
    pub t_crit: Option<f64>,
    pub grade: Grade,
    pub estimated: Vec<&'static str>,
}

/// Three NASG phases reproducing the 1-atm transition temperatures and
/// latent heats exactly from explicit inputs.
pub fn fit_inputs(i: &FitInputs) -> ThreePhase {
    let (tm, tb) = (i.t_melt, i.t_boil);
    let mut liquid = nasg_fit(i.rho_l, i.c_l, i.cp_l, i.alpha_l, i.t_ref_l, P_ATM);
    let mut solid = nasg_fit(i.rho_s, i.c_bulk_s, i.cp_s, i.alpha_s, i.t_ref_s, P_ATM);
    let mut vapor = Nasg { gamma: i.gamma_v, p_inf: 0.0, b: 0.0, q: 0.0, cv: i.cv_v, q_prime: 0.0 };
    let cp_v = i.gamma_v * i.cv_v;
    // energy references: q_l = q'_l = 0
    liquid.q = 0.0;
    liquid.q_prime = 0.0;
    if i.sublimes {
        // at 1 atm the solid goes straight to vapour at t_boil (sublimation);
        // the liquid exists only above the triple point (≈ t_melt)
        solid.q = 0.0;
        solid.q_prime = 0.0;
        vapor.q = i.h_vap + (i.cp_s - cp_v) * tb;
        vapor.q_prime = (gibbs_a(&vapor, P_ATM, tb) - gibbs_a(&solid, P_ATM, tb)) / tb;
        let p_t = i.p_triple.unwrap_or(P_ATM);
        liquid.q = solid.q + i.h_fus + (i.cp_s - i.cp_l) * tm;
        liquid.q_prime = (gibbs_a(&liquid, p_t, tm) - gibbs_a(&solid, p_t, tm)) / tm;
    } else {
        vapor.q = i.h_vap + (i.cp_l - cp_v) * tb;
        vapor.q_prime = (gibbs_a(&vapor, P_ATM, tb) - gibbs_a(&liquid, P_ATM, tb)) / tb;
        solid.q = (i.cp_l - i.cp_s) * tm - i.h_fus;
        solid.q_prime = (gibbs_a(&solid, P_ATM, tm) - gibbs_a(&liquid, P_ATM, tm)) / tm;
    }
    ThreePhase { solid, liquid, vapor, shear_speed: i.shear_speed, molar_mass: i.molar_mass, grade: i.grade, estimated: i.estimated.clone(), t_crit: i.t_crit }
}

impl ThreePhase {
    pub fn eos(&self, ph: Phase) -> &Nasg {
        match ph {
            Phase::Solid => &self.solid,
            Phase::Liquid => &self.liquid,
            Phase::Vapor => &self.vapor,
        }
    }

    /// Stable phase at (p, T): the lowest Gibbs energy; above the critical
    /// temperature the fluid (reported as vapour) unless the solid wins.
    pub fn stable(&self, p: f64, t: f64) -> Phase {
        if self.t_crit.is_some_and(|tc| t > tc) {
            return if self.solid.gibbs(p, t) < self.vapor.gibbs(p, t) { Phase::Solid } else { Phase::Vapor };
        }
        let g = [Phase::Solid, Phase::Liquid, Phase::Vapor].map(|ph| (ph, self.eos(ph).gibbs(p, t)));
        g.iter().cloned().fold(g[0], |a, b| if b.1 < a.1 { b } else { a }).0
    }

    /// Coexistence temperature of phases a and b at pressure p (bisection on
    /// g_a − g_b over [t_lo, t_hi]); None if no sign change.
    pub fn coexistence_t(&self, a: Phase, b: Phase, p: f64, t_lo: f64, t_hi: f64) -> Option<f64> {
        let f = |t: f64| self.eos(a).gibbs(p, t) - self.eos(b).gibbs(p, t);
        let (mut lo, mut hi) = (t_lo, t_hi);
        let (flo, fhi) = (f(lo), f(hi));
        if !(flo * fhi < 0.0) {
            return None;
        }
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if (f(mid) < 0.0) == (flo < 0.0) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some(0.5 * (lo + hi))
    }

    /// Temperature at which the vapour becomes the stable phase at pressure
    /// p (boiling or sublimation): g_v = min(g_s, g_l).
    pub fn vaporization_t(&self, p: f64, t_lo: f64, t_hi: f64) -> Option<f64> {
        let f = |t: f64| {
            let gc = self.solid.gibbs(p, t).min(self.liquid.gibbs(p, t));
            self.vapor.gibbs(p, t) - gc
        };
        let (mut lo, mut hi) = (t_lo, t_hi);
        let (flo, fhi) = (f(lo), f(hi));
        if !(flo * fhi < 0.0) {
            return None;
        }
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if (f(mid) < 0.0) == (flo < 0.0) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some(0.5 * (lo + hi))
    }

    /// Triple point (T, p): where the solid–liquid and liquid–vapour lines
    /// meet. Scans ln p for a sign change of T_sl(p) − T_lv(p) (either line
    /// may not exist at extreme pressures), then bisects.
    pub fn triple_point(&self, t_guess: f64) -> Option<(f64, f64)> {
        let (tl, th) = (0.2 * t_guess, 5.0 * t_guess);
        // coexistence near the guess first: with very different phase heat
        // capacities (ice/water) the Gibbs curves can cross a second,
        // unphysical time far from it, and the wide bracket then has no
        // sign change
        let co = |a: Phase, b: Phase, p: f64| self.coexistence_t(a, b, p, 0.7 * t_guess, 1.4 * t_guess).or_else(|| self.coexistence_t(a, b, p, tl, th));
        let diff = |lp: f64| -> Option<f64> {
            let p = lp.exp();
            let tsl = co(Phase::Solid, Phase::Liquid, p)?;
            let tlv = co(Phase::Liquid, Phase::Vapor, p)?;
            Some(tsl - tlv)
        };
        let n = 240;
        let (a, b) = (1e-4f64.ln(), 1e9f64.ln());
        let grid: Vec<(f64, Option<f64>)> = (0..=n).map(|i| {
            let lp = a + (b - a) * i as f64 / n as f64;
            (lp, diff(lp))
        }).collect();
        for w in grid.windows(2) {
            if let ((l0, Some(d0)), (l1, Some(d1))) = (w[0], w[1]) {
                if d0 * d1 <= 0.0 {
                    let (mut lo, mut hi, dlo) = (l0, l1, d0);
                    for _ in 0..100 {
                        let mid = 0.5 * (lo + hi);
                        match diff(mid) {
                            Some(d) if (d < 0.0) == (dlo < 0.0) => lo = mid,
                            Some(_) => hi = mid,
                            None => break,
                        }
                    }
                    let p = (0.5 * (lo + hi)).exp();
                    let t = co(Phase::Solid, Phase::Liquid, p)?;
                    return Some((t, p));
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_loads_and_fits_reproduce_1_atm_transitions() {
        let els = all();
        assert_eq!(els.len(), 118);
        assert_eq!(els[28].symbol, "Cu");
        let cu = fit(&els[28]).unwrap();
        let tm = cu.coexistence_t(Phase::Solid, Phase::Liquid, P_ATM, 300.0, 3000.0).unwrap();
        let tb = cu.coexistence_t(Phase::Liquid, Phase::Vapor, P_ATM, 1500.0, 5000.0).unwrap();
        assert!((tm - els[28].t_melt.unwrap()).abs() < 1e-6, "{tm}");
        assert!((tb - els[28].t_boil.unwrap()).abs() < 1e-6, "{tb}");
        assert_eq!(cu.stable(P_ATM, 298.15), Phase::Solid);
    }
}
