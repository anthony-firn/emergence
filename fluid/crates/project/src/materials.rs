//! One sourced materials library for every discipline.
//!
//! Electrical reads copper's resistivity from the same record that thermal
//! reads its conductivity from, so a team can never use a different copper.
//! Every property states its source and a relative uncertainty, and the
//! uncertainty travels into every result that uses it. Temperature-dependent
//! properties are tables, interpolated linearly. A query outside the table is
//! an error, never an extrapolation.
//!
//! In expressions a property is `mat.<material>.<property>` (the value at the
//! table's reference point for tables) or `mat.<material>.<property>(T)`.

use std::collections::BTreeMap;

use emergence::Rigor;

use crate::units::{parse_quantity, Dim};
use crate::value::Val;

/// Property names and their dimensions. A project-local material may only use
/// these names, so `k` always means thermal conductivity.
pub const PROPERTIES: &[(&str, Dim, &str)] = &[
    ("k", Dim([1, 1, -3, 0, -1, 0, 0]), "thermal conductivity"),
    ("rho", Dim([-3, 1, 0, 0, 0, 0, 0]), "density"),
    ("cp", Dim([2, 0, -2, 0, -1, 0, 0]), "specific heat"),
    ("mu", Dim([-1, 1, -1, 0, 0, 0, 0]), "dynamic viscosity"),
    ("nu", Dim([2, 0, -1, 0, 0, 0, 0]), "kinematic viscosity"),
    ("pr", Dim::NONE, "Prandtl number"),
    ("resistivity", Dim([3, 1, -3, -2, 0, 0, 0]), "electrical resistivity at 20 degC"),
    ("alpha_r", Dim([0, 0, 0, 0, -1, 0, 0]), "temperature coefficient of resistivity at 20 degC"),
];

pub fn property_dim(name: &str) -> Option<Dim> {
    PROPERTIES.iter().find(|p| p.0 == name).map(|p| p.1)
}

#[derive(Clone, Debug, PartialEq)]
pub enum Data {
    Const(f64),
    /// Values against temperature (K), increasing temperatures.
    Table(Vec<(f64, f64)>),
}

#[derive(Clone, Debug, PartialEq)]
pub struct Property {
    pub data: Data,
    pub dim: Dim,
    /// Relative uncertainty (half-width / value).
    pub rel: f64,
    /// Handbook values are `Validated`; estimates are `Asymptotic`.
    pub rigor: Rigor,
    pub source: String,
}

impl Property {
    /// Value at the reference point (the first table temperature at or above
    /// 300 K, or the constant).
    pub fn reference(&self) -> Val {
        let v = match &self.data {
            Data::Const(v) => *v,
            Data::Table(t) => t.iter().find(|p| p.0 >= 300.0).unwrap_or(&t[t.len() - 1]).1,
        };
        self.bounded(v, v, v, Rigor::Rigorous)
    }

    fn bounded(&self, nom: f64, lo: f64, hi: f64, arg_rigor: Rigor) -> Val {
        let e = self.rel * nom.abs().max(lo.abs()).max(hi.abs());
        Val { nom, lo: lo - e, hi: hi + e, dim: self.dim, rigor: self.rigor.min(arg_rigor) }
    }

    /// Value at temperature `t` (an interval): the exact range of the
    /// piecewise-linear table over the interval, widened by the uncertainty.
    pub fn at(&self, t: &Val) -> Result<Val, String> {
        if t.dim != Dim::TEMPERATURE {
            return Err(format!("material property needs a temperature, got {}", t.dim));
        }
        let tab = match &self.data {
            Data::Const(v) => return Ok(self.bounded(*v, *v, *v, t.rigor)),
            Data::Table(tab) => tab,
        };
        let (t0, t1) = (tab[0].0, tab[tab.len() - 1].0);
        if t.lo < t0 || t.hi > t1 {
            return Err(format!("temperature [{:.2}, {:.2}] K is outside the table {t0}–{t1} K", t.lo, t.hi));
        }
        let interp = |x: f64| -> f64 {
            let i = tab.windows(2).position(|w| x <= w[1].0).unwrap_or(tab.len() - 2);
            let (a, b) = (tab[i], tab[i + 1]);
            a.1 + (b.1 - a.1) * (x - a.0) / (b.0 - a.0)
        };
        let mut lo = interp(t.lo).min(interp(t.hi));
        let mut hi = interp(t.lo).max(interp(t.hi));
        for &(x, v) in tab {
            if x > t.lo && x < t.hi {
                lo = lo.min(v);
                hi = hi.max(v);
            }
        }
        Ok(self.bounded(interp(t.nom), lo, hi, t.rigor))
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Material {
    pub description: String,
    pub props: BTreeMap<String, Property>,
}

/// The library: built-in records plus any project-local ones.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Library {
    pub materials: BTreeMap<String, Material>,
}

fn c(v: f64, name: &str, rel: f64, rigor: Rigor, source: &str) -> (String, Property) {
    let dim = property_dim(name).expect("known property");
    (name.into(), Property { data: Data::Const(v), dim, rel, rigor, source: source.into() })
}

fn tab(rows: &[(f64, f64)], name: &str, rel: f64, rigor: Rigor, source: &str) -> (String, Property) {
    let dim = property_dim(name).expect("known property");
    (name.into(), Property { data: Data::Table(rows.to_vec()), dim, rel, rigor, source: source.into() })
}

const INCROPERA_A1: &str = "Incropera et al., Fundamentals of Heat and Mass Transfer, 7th ed., Table A.1 (300 K)";
const INCROPERA_A4: &str = "Incropera et al., Fundamentals of Heat and Mass Transfer, 7th ed., Table A.4";
const INCROPERA_A6: &str = "Incropera et al., Fundamentals of Heat and Mass Transfer, 7th ed., Table A.6 (300 K)";
const ASM: &str = "ASM Handbook Vol. 2 (1990), wrought aluminium alloy properties";

impl Library {
    /// The built-in records. Validated values cite a handbook; values marked
    /// "estimate" are typical datasheet figures with a wide uncertainty.
    pub fn builtin() -> Library {
        use Rigor::{Asymptotic as Est, Validated as Hb};
        let mut m = BTreeMap::new();
        let mut add = |name: &str, description: &str, props: Vec<(String, Property)>| {
            m.insert(name.to_string(), Material { description: description.into(), props: props.into_iter().collect() });
        };
        add(
            "copper",
            "Copper, annealed (conductors, PCB foil, heat spreaders)",
            vec![
                // IACS: 100 % conductivity is 1/58 Ω·mm²/m at 20 °C by definition
                // (IEC 60028). Electrodeposited PCB foil runs a few % higher:
                // the 3 % uncertainty is an estimate covering it.
                c(1.0 / 58e6, "resistivity", 0.03, Est, "IEC 60028 (IACS) at 20 °C; ±3 % estimate for PCB foil"),
                c(0.00393, "alpha_r", 0.01, Hb, "IEC 60028, temperature coefficient at 20 °C"),
                c(401.0, "k", 0.05, Hb, INCROPERA_A1),
                c(8933.0, "rho", 0.005, Hb, INCROPERA_A1),
                c(385.0, "cp", 0.01, Hb, INCROPERA_A1),
            ],
        );
        add(
            "al6063",
            "Aluminium 6063-T5 (extruded heat sinks)",
            vec![
                c(209.0, "k", 0.05, Hb, ASM),
                c(2700.0, "rho", 0.01, Hb, ASM),
                c(900.0, "cp", 0.03, Hb, ASM),
            ],
        );
        add(
            "al6061",
            "Aluminium 6061-T6 (machined parts, enclosures)",
            vec![
                c(167.0, "k", 0.05, Hb, ASM),
                c(2700.0, "rho", 0.01, Hb, ASM),
                c(896.0, "cp", 0.03, Hb, ASM),
            ],
        );
        add(
            "silicon",
            "Silicon, pure (dies)",
            vec![c(148.0, "k", 0.05, Hb, INCROPERA_A1), c(2330.0, "rho", 0.005, Hb, INCROPERA_A1), c(712.0, "cp", 0.01, Hb, INCROPERA_A1)],
        );
        add(
            "fr4",
            "FR-4 laminate, through-plane",
            vec![
                c(0.3, "k", 0.2, Est, "estimate: typical laminate datasheets 0.25–0.35 W/(m·K) through-plane"),
                c(1850.0, "rho", 0.05, Est, "estimate: typical laminate datasheets"),
                c(1100.0, "cp", 0.15, Est, "estimate: typical laminate datasheets"),
            ],
        );
        add(
            "pla",
            "PLA, as printed (FDM)",
            vec![
                c(0.13, "k", 0.25, Est, "estimate: literature range 0.11–0.20 W/(m·K), print dependent"),
                c(1240.0, "rho", 0.05, Est, "NatureWorks Ingeo 4043D datasheet (solid); printed parts are lighter"),
            ],
        );
        add(
            "abs",
            "ABS, as printed (FDM)",
            vec![
                c(0.17, "k", 0.25, Est, "estimate: typical datasheets"),
                c(1040.0, "rho", 0.05, Est, "estimate: typical datasheets (solid)"),
            ],
        );
        add(
            "air",
            "Dry air at 1 atm",
            vec![
                tab(&[(250.0, 1.3947), (300.0, 1.1614), (350.0, 0.9950), (400.0, 0.8711)], "rho", 0.01, Hb, INCROPERA_A4),
                tab(&[(250.0, 1006.0), (300.0, 1007.0), (350.0, 1009.0), (400.0, 1014.0)], "cp", 0.01, Hb, INCROPERA_A4),
                tab(&[(250.0, 159.6e-7), (300.0, 184.6e-7), (350.0, 208.2e-7), (400.0, 230.1e-7)], "mu", 0.01, Hb, INCROPERA_A4),
                tab(&[(250.0, 11.44e-6), (300.0, 15.89e-6), (350.0, 20.92e-6), (400.0, 26.41e-6)], "nu", 0.01, Hb, INCROPERA_A4),
                tab(&[(250.0, 22.3e-3), (300.0, 26.3e-3), (350.0, 30.0e-3), (400.0, 33.8e-3)], "k", 0.01, Hb, INCROPERA_A4),
                tab(&[(250.0, 0.720), (300.0, 0.707), (350.0, 0.700), (400.0, 0.690)], "pr", 0.01, Hb, INCROPERA_A4),
            ],
        );
        add(
            "water",
            "Liquid water at 300 K",
            vec![
                c(997.0, "rho", 0.005, Hb, INCROPERA_A6),
                c(4179.0, "cp", 0.005, Hb, INCROPERA_A6),
                c(855e-6, "mu", 0.01, Hb, INCROPERA_A6),
                c(0.613, "k", 0.01, Hb, INCROPERA_A6),
                c(5.83, "pr", 0.01, Hb, INCROPERA_A6),
            ],
        );
        Library { materials: m }
    }

    pub fn get(&self, material: &str, prop: &str) -> Result<&Property, String> {
        let m = self.materials.get(material).ok_or_else(|| format!("unknown material `{material}`"))?;
        m.props.get(prop).ok_or_else(|| {
            format!(
                "material `{material}` has no `{prop}` (has: {})",
                m.props.keys().cloned().collect::<Vec<_>>().join(", ")
            )
        })
    }

    /// Resolve `mat.<material>.<prop>` into its parts.
    pub fn split(id: &str) -> Option<(&str, &str)> {
        let rest = id.strip_prefix("mat.")?;
        rest.split_once('.')
    }

    /// Add or replace a project-local material property from text:
    /// `value` like `"3 W/(m*K)"`, `rel` a relative uncertainty.
    pub fn define(&mut self, material: &str, prop: &str, value: &str, rel: f64, source: &str) -> Result<(), String> {
        let dim = property_dim(prop).ok_or_else(|| {
            format!(
                "unknown property `{prop}` (known: {})",
                PROPERTIES.iter().map(|p| p.0).collect::<Vec<_>>().join(", ")
            )
        })?;
        let (q, unit) = parse_quantity(value)?;
        if q.dim != dim || unit.is_affine() {
            return Err(format!("{material}.{prop} = {value}: expected {}", dim.symbol()));
        }
        let m = self
            .materials
            .entry(material.to_string())
            .or_insert_with(|| Material { description: format!("project material `{material}`"), props: BTreeMap::new() });
        // Project values are what a team asserts from a datasheet: validated
        // when sourced, an estimate otherwise.
        let rigor = if source.is_empty() || source.contains("estimate") { Rigor::Asymptotic } else { Rigor::Validated };
        m.props.insert(prop.into(), Property { data: Data::Const(q.si), dim, rel, rigor, source: source.into() });
        Ok(())
    }
}
