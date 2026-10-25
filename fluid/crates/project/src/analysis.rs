//! Analyses: the typed models a project can run.
//!
//! An analysis declares its inputs and outputs with dimensions, so a binding
//! in the wrong unit is refused before anything runs. It also declares
//!
//! * which inputs it is **monotone** in (coordinate-wise, for any values of
//!   the others). When every toleranced input is, the extremes over the
//!   tolerance box lie on its corners, and corner enumeration gives a
//!   rigorous enclosure. Otherwise the corners are a sample and the result is
//!   downgraded to asymptotic;
//! * its **assumptions**, printed with every result. The error claim covers
//!   the numerics *within* those assumptions and nothing else.
//!
//! Inputs may default to a material property through a material slot
//! ([`MatRef`]), which is how electrical and thermal end up reading the same
//! copper.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Mutex;

use emergence::Rigor;

use crate::thermal;
use crate::units::{parse_unit, Dim};

/// A material property an input defaults to.
#[derive(Clone, Debug)]
pub struct MatRef {
    /// Material slot name (bound to a material in the project).
    pub slot: String,
    pub prop: String,
    /// Input supplying the temperature for tabulated properties.
    pub at: Option<String>,
}

#[derive(Clone, Debug)]
pub struct Port {
    pub name: String,
    /// `None`: a CAD parameter, accepted as a length or a pure number.
    pub dim: Option<Dim>,
    /// Display unit.
    pub unit: String,
    pub doc: String,
    pub monotone: bool,
    /// Used only to look up material properties; not an argument of `run`'s
    /// arithmetic, so its tolerance needs no corners of its own.
    pub passive: bool,
    pub from: Option<MatRef>,
    /// When unbound, reuse the binding of this other input.
    pub default: Option<String>,
}

#[derive(Clone, Debug)]
pub struct OutPort {
    pub name: String,
    pub dim: Dim,
    pub unit: String,
    pub doc: String,
}

/// One run at one point.
#[derive(Clone, Debug)]
pub struct Outcome {
    /// Per output, SI.
    pub values: Vec<f64>,
    /// Absolute error claim per output, SI.
    pub errors: Vec<f64>,
    pub rigor: Rigor,
    /// Validity indicators and regime, for the report.
    pub notes: Vec<String>,
}

pub trait Analysis: Send + Sync {
    fn kind(&self) -> &str;
    fn summary(&self) -> &str;
    fn inputs(&self) -> Vec<Port>;
    fn outputs(&self) -> Vec<OutPort>;
    fn assumptions(&self) -> Vec<String>;
    /// Run at one point; `x` follows [`Analysis::inputs`], in SI.
    fn run(&self, x: &[f64]) -> Result<Outcome, String>;
}

fn port(name: &str, unit: &str, doc: &str, monotone: bool) -> Port {
    Port {
        name: name.into(),
        dim: Some(parse_unit(unit).expect("valid unit").dim),
        unit: unit.into(),
        doc: doc.into(),
        monotone,
        passive: false,
        from: None,
        default: None,
    }
}

fn mat(mut p: Port, slot: &str, prop: &str, at: Option<&str>) -> Port {
    p.from = Some(MatRef { slot: slot.into(), prop: prop.into(), at: at.map(String::from) });
    p
}

fn out(name: &str, unit: &str, doc: &str) -> OutPort {
    OutPort { name: name.into(), dim: parse_unit(unit).expect("valid unit").dim, unit: unit.into(), doc: doc.into() }
}

/// A new instance of a built-in analysis kind.
pub fn builtin(kind: &str) -> Option<Box<dyn Analysis>> {
    match kind {
        "pcb_trace" => Some(Box::new(PcbTrace)),
        "conduction_slab" => Some(Box::new(Slab)),
        "plate_fin_heatsink" => Some(Box::new(PlateFinSink)),
        "pcb_microstrip" => Some(Box::new(Microstrip)),
        _ => None,
    }
}

pub const KINDS: &[&str] = &[
    "pcb_trace",
    "pcb_microstrip",
    "conduction_slab",
    "plate_fin_heatsink",
    "chip_design",
    "zone_loads",
    "refrigeration_cycle",
    "duct_run",
    "fan",
    "cooling_coil",
    "thermal_comfort",
];

/// Build an analysis from its kind and settings; `dir` is where relative
/// files in `config` resolve. Returns the model and its cache identity.
pub fn build(kind: &str, config: &std::collections::BTreeMap<String, String>, dir: &std::path::Path) -> Result<(Box<dyn Analysis>, String), String> {
    if kind == "chip_design" {
        let rtl = config.get("rtl").ok_or("chip_design needs config.rtl, the Verilog file")?;
        let d = crate::chip_design::ChipDesign::load(dir, rtl)?;
        let key = d.key();
        return Ok((Box::new(d), key));
    }
    use crate::hvac_analyses as hv;
    let need = |k: &str| config.get(k).ok_or_else(|| format!("`{kind}` needs config.{k}"));
    match kind {
        "zone_loads" => {
            let z = hv::ZoneLoads::load(dir, need("zone")?)?;
            let key = z.key();
            return Ok((Box::new(z), key));
        }
        "refrigeration_cycle" => {
            let r = need("refrigerant")?;
            return Ok((Box::new(hv::Refrigeration::new(r)?), format!("refrigeration:{r}")));
        }
        "fan" => {
            let c = need("curve")?;
            return Ok((Box::new(hv::Fan::new(c)?), format!("fan:{c}")));
        }
        "duct_run" => return Ok((Box::new(hv::DuctRun), kind.into())),
        "cooling_coil" => return Ok((Box::new(hv::Coil), kind.into())),
        "thermal_comfort" => return Ok((Box::new(hv::Comfort), kind.into())),
        _ => {}
    }
    if let Some(k) = config.keys().next() {
        return Err(format!("`{kind}` takes no config (found `{k}`)"));
    }
    builtin(kind).map(|m| (m, kind.to_string())).ok_or_else(|| format!("unknown kind `{kind}` (known: {})", KINDS.join(", ")))
}

/// DC resistance and loss of a PCB trace.
pub struct PcbTrace;

/// Reference temperature of resistivity data (IEC 60028).
const T20: f64 = 293.15;

impl Analysis for PcbTrace {
    fn kind(&self) -> &str {
        "pcb_trace"
    }
    fn summary(&self) -> &str {
        "DC resistance, voltage drop and Joule loss of a straight PCB trace"
    }
    fn inputs(&self) -> Vec<Port> {
        vec![
            port("width", "mm", "trace width", true),
            port("thickness", "um", "copper thickness", true),
            port("length", "mm", "trace length", true),
            port("current", "A", "DC (or RMS) current", true),
            port("temperature", "degC", "copper temperature", true),
            mat(port("resistivity", "Ohm*m", "resistivity at 20 degC", true), "conductor", "resistivity", None),
            mat(port("alpha_r", "1/K", "temperature coefficient at 20 degC", true), "conductor", "alpha_r", None),
        ]
    }
    fn outputs(&self) -> Vec<OutPort> {
        vec![
            out("resistance", "mOhm", "DC resistance"),
            out("voltage_drop", "mV", "voltage drop"),
            out("power", "W", "Joule loss"),
        ]
    }
    fn assumptions(&self) -> Vec<String> {
        vec![
            "DC current (no skin or proximity effect)".into(),
            "uniform rectangular cross-section; etch profile inside the width tolerance".into(),
            "copper at one uniform temperature".into(),
        ]
    }
    fn run(&self, x: &[f64]) -> Result<Outcome, String> {
        let [w, t, l, i, temp, rho, a] = x[..] else { unreachable!() };
        if w <= 0.0 || t <= 0.0 || l <= 0.0 {
            return Err("trace dimensions must be positive".into());
        }
        let r = rho * (1.0 + a * (temp - T20)) * l / (w * t);
        Ok(Outcome { values: vec![r, i * r, i * i * r], errors: vec![0.0; 3], rigor: Rigor::Rigorous, notes: vec![] })
    }
}

/// Characteristic impedance of a surface microstrip, from a field solve of
/// the drawn cross-section (`cad::laplace`), copper thickness included.
pub struct Microstrip;

/// Solve one microstrip cross-section: half the line against a mirror, a
/// grounded box `span` trace-plus-substrate widths across and `height`
/// substrate heights tall, at cell `cell`.
fn microstrip_field(w: f64, h: f64, t: f64, eps_r: f64, cell: f64, box_scale: f64) -> (f64, f64) {
    use cad::laplace::{quasi_tem, Edge, Section};
    use cad::Region2;
    let span = 6.0 * box_scale * (w + 2.0 * h);
    let height = h + 8.0 * box_scale * h;
    let section = Section::new(([0.0, 0.0], [0.5 * span, height]), cell)
        .with_dielectric(Region2::rect([0.25 * span, 0.5 * h], [0.5 * span, 0.5 * h]), eps_r)
        .with_conductor(Region2::rect([0.25 * span, 0.0], [0.5 * span, 1e-15]), 0.0)
        .with_conductor(Region2::rect([0.0, h + 0.5 * t], [0.5 * w, 0.5 * t]), 1.0)
        .with_edges([Edge::Mirror, Edge::Grounded, Edge::Grounded, Edge::Grounded])
        .with_multiplicity(2.0);
    // Capacitance is an energy: its error goes as the residual squared.
    let (z0, eps_eff, _) = quasi_tem(&section, 1.0, 1e-8, 40000);
    (z0, eps_eff)
}

impl Analysis for Microstrip {
    fn kind(&self) -> &str {
        "pcb_microstrip"
    }
    fn summary(&self) -> &str {
        "Characteristic impedance of a surface microstrip from a quasi-TEM field solve of its cross-section"
    }
    fn inputs(&self) -> Vec<Port> {
        // Enlarging a conductor or raising a permittivity can only raise a
        // capacitance (Dirichlet's principle), and Z0 = 1/(c √(C C_air)):
        // monotone in width, thickness and permittivity. Height moves the
        // trace and grows the substrate at once, so it is not declared.
        vec![
            port("width", "mm", "trace width", true),
            port("thickness", "um", "copper thickness", true),
            port("height", "mm", "dielectric height to the reference plane", false),
            port("eps_r", "1", "dielectric relative permittivity at the signal frequency", true),
        ]
    }
    fn outputs(&self) -> Vec<OutPort> {
        vec![
            out("z0", "Ohm", "characteristic impedance"),
            out("eps_eff", "1", "effective permittivity"),
            out("delay", "ps/mm", "propagation delay"),
            out("z0_formula", "Ohm", "Hammerstad–Jensen zero-thickness formula (cross-check)"),
        ]
    }
    fn assumptions(&self) -> Vec<String> {
        vec![
            "quasi-TEM, lossless, uncoated (no solder mask), rectangular copper".into(),
            "single trace over an infinite reference plane, no neighbouring copper".into(),
            "permittivity is the laminate's value at the signal frequency".into(),
        ]
    }
    fn run(&self, x: &[f64]) -> Result<Outcome, String> {
        let [w, t, h, er] = x[..] else { unreachable!() };
        if w <= 0.0 || t <= 0.0 || h <= 0.0 || er < 1.0 {
            return Err("microstrip needs positive dimensions and eps_r ≥ 1".into());
        }
        let c0 = 299_792_458.0;
        let (zf, eff_f) = cad::kits::pcb::microstrip_z0(w, h, er);
        let u = w / h;
        if !(0.05..=20.0).contains(&u) {
            return Err(format!("W/h = {u:.3} is outside 0.05–20, where neither rung is validated"));
        }
        if t / h < 1.0 / 12.0 {
            // Rung 0: copper too thin for the uniform-grid solve to draw.
            // Hammerstad–Jensen (quoted 1 %) assumes zero thickness; Wheeler's
            // effective-width increase Δw = (t/π)(1 + ln(2h/t)) measures what
            // that omits, and is counted as error, not applied.
            let dw = t / std::f64::consts::PI * (1.0 + (2.0 * h / t).ln());
            let (zt, _) = cad::kits::pcb::microstrip_z0(w + dw, h, er);
            let err = 0.01 * zf + (zf - zt).abs();
            return Ok(Outcome {
                values: vec![zf, eff_f, eff_f.sqrt() / c0, zf],
                errors: vec![err, 0.01 * eff_f, 0.005 * eff_f.sqrt() / c0, 0.01 * zf],
                rigor: Rigor::Asymptotic,
                notes: vec![format!(
                    "t/h = {:.3} < 1/12: closed form (Hammerstad–Jensen), thickness effect {:.2} Ω counted as error",
                    t / h,
                    (zf - zt).abs()
                )],
            });
        }
        // Rung 1: the field solve, on a grid that puts the copper's thickness
        // and half-width on grid lines at both resolutions, so the only
        // discretisation error left is the smooth one.
        let misfit = |n: f64| -> f64 {
            [t, 0.5 * w]
                .iter()
                .map(|len| {
                    let (k1, k2) = (len * n / h, len * 0.5 * n / h);
                    (k1 - k1.round()).abs() + (k2 - k2.round()).abs()
                })
                .sum()
        };
        let n = (6..=12).map(|k| 2.0 * k as f64).min_by(|a, b| misfit(*a).total_cmp(&misfit(*b))).unwrap();
        let cell = h / n;
        let ((z, eps_eff), ((zc, _), (zb, _))) = rayon::join(
            || microstrip_field(w, h, t, er, cell, 1.0),
            || rayon::join(|| microstrip_field(w, h, t, er, 2.0 * cell, 1.0), || microstrip_field(w, h, t, er, 2.0 * cell, 1.5)),
        );
        if !(z.is_finite() && zc.is_finite() && zb.is_finite()) {
            return Err("the field solve did not resolve the trace".into());
        }
        // Discretisation (cell halved) plus truncation (box enlarged).
        let err = (z - zc).abs() + (zb - zc).abs();
        let delay = eps_eff.sqrt() / c0;
        Ok(Outcome {
            values: vec![z, eps_eff, delay, zf],
            errors: vec![err, eps_eff * err / z, delay * 0.5 * err / z, 0.01 * zf],
            rigor: Rigor::Asymptotic,
            notes: vec![
                format!("W/h = {u:.2}, t/h = {:.2}; field solve at h/{n:.0} with misfit {:.2} cells", t / h, misfit(n)),
                format!(
                    "field solve {z:.2} Ω ± {err:.2}; zero-thickness formula {zf:.2} Ω ({:+.1} %)",
                    100.0 * (zf / z - 1.0)
                ),
            ],
        })
    }
}

/// Conduction through a flat layer (a thermal interface pad, a board).
pub struct Slab;

impl Analysis for Slab {
    fn kind(&self) -> &str {
        "conduction_slab"
    }
    fn summary(&self) -> &str {
        "1-D conduction resistance of a flat layer"
    }
    fn inputs(&self) -> Vec<Port> {
        vec![
            port("thickness", "mm", "layer thickness", true),
            port("area", "mm^2", "area the heat crosses", true),
            mat(port("k", "W/(m*K)", "thermal conductivity", true), "material", "k", None),
        ]
    }
    fn outputs(&self) -> Vec<OutPort> {
        vec![out("resistance", "K/W", "thermal resistance")]
    }
    fn assumptions(&self) -> Vec<String> {
        vec!["1-D conduction: uniform flux over the area, no spreading, no contact resistance".into()]
    }
    fn run(&self, x: &[f64]) -> Result<Outcome, String> {
        let [t, a, k] = x[..] else { unreachable!() };
        if t <= 0.0 || a <= 0.0 || k <= 0.0 {
            return Err("thickness, area and k must be positive".into());
        }
        Ok(Outcome { values: vec![t / (k * a)], errors: vec![0.0], rigor: Rigor::Rigorous, notes: vec![] })
    }
}

/// Plate-fin heat sink in ducted laminar forced air.
pub struct PlateFinSink;

/// Laminar-turbulent transition in a parallel-plate channel (Incropera §8.1;
/// the usual design threshold).
const RE_TRANSITION: f64 = 2300.0;

/// The developing-channel solution depends only on (Pr, X), and corner
/// enumeration revisits the same pairs many times.
fn channel_cached(pr: f64, x: f64) -> thermal::Channel {
    static CACHE: Mutex<Option<HashMap<(u64, u64), thermal::Channel>>> = Mutex::new(None);
    let key = (pr.to_bits(), x.to_bits());
    if let Some(c) = CACHE.lock().unwrap().get_or_insert_with(HashMap::new).get(&key) {
        return *c;
    }
    let c = thermal::channel(pr, x);
    CACHE.lock().unwrap().get_or_insert_with(HashMap::new).insert(key, c);
    c
}

impl Analysis for PlateFinSink {
    fn kind(&self) -> &str {
        "plate_fin_heatsink"
    }
    fn summary(&self) -> &str {
        "Thermal resistance and pressure drop of a plate-fin heat sink in ducted laminar air"
    }
    fn inputs(&self) -> Vec<Port> {
        let mut t = port("air_temperature", "degC", "inlet air temperature", true);
        t.passive = true;
        vec![
            port("width", "mm", "base width, across the flow", false),
            port("length", "mm", "base length, along the flow", false),
            port("base_thickness", "mm", "base plate thickness (thicker: more conduction, less spreading)", false),
            port("fin_height", "mm", "fin height above the base", false),
            port("fin_thickness", "mm", "fin thickness", false),
            port("fin_count", "1", "number of fins", false),
            port("approach_velocity", "m/s", "air velocity ahead of the sink (duct cross-section = fin envelope)", true),
            t,
            Port { default: Some("width".into()), ..port("source_width", "mm", "heat source width, centred on the base (default: whole base)", false) },
            Port { default: Some("length".into()), ..port("source_length", "mm", "heat source length along the flow (default: whole base)", false) },
            mat(port("fin_k", "W/(m*K)", "fin conductivity", true), "fin", "k", None),
            mat(port("air_k", "W/(m*K)", "air conductivity", true), "coolant", "k", Some("air_temperature")),
            mat(port("air_nu", "m^2/s", "air kinematic viscosity", false), "coolant", "nu", Some("air_temperature")),
            mat(port("air_pr", "1", "air Prandtl number", false), "coolant", "pr", Some("air_temperature")),
            mat(port("air_rho", "kg/m^3", "air density", true), "coolant", "rho", Some("air_temperature")),
            mat(port("air_cp", "J/(kg*K)", "air specific heat", true), "coolant", "cp", Some("air_temperature")),
        ]
    }
    fn outputs(&self) -> Vec<OutPort> {
        vec![
            out("thermal_resistance", "K/W", "base-to-inlet-air resistance"),
            out("h", "W/(m^2*K)", "mean heat transfer coefficient"),
            out("fin_efficiency", "1", "fin efficiency"),
            out("channel_velocity", "m/s", "mean velocity between fins"),
            out("reynolds", "1", "channel Reynolds number on D_h = 2·gap"),
            out("pressure_drop", "Pa", "core friction pressure drop"),
            out("spreading_resistance", "K/W", "spreading from the source into the base (included in thermal_resistance)"),
        ]
    }
    fn assumptions(&self) -> Vec<String> {
        vec![
            "ducted: all air passes between the fins (no bypass, fin tips touch an adiabatic shroud)".into(),
            "laminar, uniform inlet velocity; channels treated as parallel plates (fin height ≫ gap)".into(),
            "fin walls at uniform temperature for h; 1-D fins with adiabatic tips".into(),
            "spreading: a centred uniform-flux source; the fins act as a uniform effective h under the base".into(),
            "pressure drop is core friction only (no inlet contraction or exit expansion losses)".into(),
            "air properties at the inlet temperature".into(),
        ]
    }
    fn run(&self, x: &[f64]) -> Result<Outcome, String> {
        let [w, l, tb, hf, tf, n, v, _t, sw, sl, kf, ka, nu, pr, rho, cp] = x[..] else { unreachable!() };
        if sw <= 0.0 || sl <= 0.0 || sw > w * (1.0 + 1e-12) || sl > l * (1.0 + 1e-12) {
            return Err(format!(
                "heat source {:.1} × {:.1} mm does not fit on the {:.1} × {:.1} mm base",
                sw * 1e3, sl * 1e3, w * 1e3, l * 1e3
            ));
        }
        if n.fract() != 0.0 || n < 2.0 {
            return Err(format!("fin_count must be a whole number ≥ 2, got {n}"));
        }
        let gap = (w - n * tf) / (n - 1.0);
        if gap <= 0.0 || hf <= 0.0 || l <= 0.0 || tb <= 0.0 {
            return Err(format!("fins do not fit: gap = {:.3} mm", gap * 1e3));
        }
        let u = v * w / ((n - 1.0) * gap);
        let dh = 2.0 * gap;
        let re = u * dh / nu;
        if re > RE_TRANSITION {
            return Err(format!(
                "channel Reynolds number {re:.0} exceeds {RE_TRANSITION}: turbulent channel flow has no rung yet"
            ));
        }
        let a = 0.5 * gap;
        let x_end = l * nu / (u * a * a);
        let ch = channel_cached(pr, x_end);
        let h = ch.nu_mean * ka / dh;
        let h_err = ch.nu_error * ka / dh;
        let m = (2.0 * h / (kf * tf)).sqrt();
        let eta = (m * hf).tanh() / (m * hf);
        let bi = h * tf / (2.0 * kf);
        let cap = rho * u * (n - 1.0) * gap * hf * cp;
        let area = (n - 1.0) * l * (gap + 2.0 * eta * hf);
        let ntu = h * area / cap;
        let r_conv = 1.0 / (cap * (1.0 - (-ntu).exp()));
        let r_base = tb / (kf * w * l);
        // The fin array, seen from the base, is a uniform effective h.
        let h_eff = 1.0 / (r_conv * w * l);
        let (r_sp, sp_err) = thermal::spreading(w, l, sw.min(w), sl.min(l), tb, kf, h_eff);
        let r = r_conv + r_base + r_sp;
        // Conductance is linear in h at most (x e^-x / (1 - e^-x) ≤ 1), and the
        // 1-D fin errs at O(Bi): the claim adds both, plus series truncation.
        let r_err = r_conv * (h_err / h + bi) + sp_err;
        let dp = ch.dp * rho * u * u;
        let dp_err = ch.dp_error * rho * u * u;
        let xstar = l / (dh * re * pr);
        let regime = if xstar > 0.05 {
            "thermally fully developed over most of the length"
        } else if xstar < 0.005 {
            "boundary layers still thin (flat-plate-like)"
        } else {
            "developing"
        };
        Ok(Outcome {
            values: vec![r, h, eta, u, re, dp, r_sp],
            errors: vec![r_err, h_err, bi * eta, 0.0, 0.0, dp_err, sp_err],
            rigor: Rigor::Asymptotic,
            notes: vec![
                format!("Re_Dh = {re:.0}, Pe = {:.0}, x* = L/(D_h Re Pr) = {xstar:.4}: {regime}", re * pr),
                format!("fin Biot number {bi:.2e} (1-D fin error ~Bi)"),
                format!("fin height / gap = {:.1} (parallel-plate channel assumed)", hf / gap),
                format!("Nu_mean = {:.3} ± {:.1e} (fine vs coarse march)", ch.nu_mean, ch.nu_error),
            ],
        })
    }
}

/// A CAD part as an analysis: geometry in, mass properties out.
pub struct Part {
    pub model: cad::model::Model,
    /// Directory the model's own file references resolve against.
    pub base: PathBuf,
    /// CAD parameter names bound by the project, in input order.
    pub params: Vec<String>,
    /// Mesh cell size (m).
    pub cell: f64,
    /// Printed by fused deposition.
    pub fdm: bool,
}

impl Part {
    /// The solid at the given CAD parameter values (SI).
    pub fn solid(&self, x: &[f64]) -> Result<cad::Solid, String> {
        let mut p = cad::Params::new();
        for (name, v) in self.params.iter().zip(x) {
            p.set(name, *v);
        }
        self.model.build(&p, &self.base)
    }
}

impl Analysis for Part {
    fn kind(&self) -> &str {
        "part"
    }
    fn summary(&self) -> &str {
        "CAD part: volume, area, mass and envelope from the parametric model"
    }
    fn inputs(&self) -> Vec<Port> {
        let mut v: Vec<Port> = self
            .params
            .iter()
            .map(|n| Port {
                name: n.clone(),
                dim: None,
                unit: "m".into(),
                doc: format!("CAD parameter `{n}`"),
                monotone: false,
                passive: false,
                from: None,
                default: None,
            })
            .collect();
        v.push(mat(port("density", "kg/m^3", "material density", true), "material", "rho", None));
        if self.fdm {
            v.push(port("nozzle", "mm", "extrusion line width (≈ nozzle diameter)", true));
            v.push(port("perimeters", "1", "wall lines (perimeters; skins taken as the same thickness)", true));
            v.push(port("infill", "%", "infill density", true));
            v.push(port("max_flow", "mm^3/s", "hotend's maximum volumetric flow", true));
        }
        v
    }
    fn outputs(&self) -> Vec<OutPort> {
        let mut v = vec![
            out("volume", "cm^3", "solid volume"),
            out("area", "cm^2", "surface area"),
            out("mass", "g", "mass"),
            out("size_x", "mm", "envelope along x"),
            out("size_y", "mm", "envelope along y"),
            out("size_z", "mm", "envelope along z"),
        ];
        if self.fdm {
            v.extend([
                out("printed_mass", "g", "mass as printed (walls + infill)"),
                out("print_time_min", "min", "least possible print time (extrusion-limited)"),
                out("filament_length", "m", "1.75 mm filament used"),
                out("min_wall", "mm", "thinnest wall in the part"),
            ]);
        }
        v
    }
    fn assumptions(&self) -> Vec<String> {
        let mut a = vec!["mass properties from the meshed solid; error estimated from a mesh twice as coarse".into()];
        if self.fdm {
            a.push("printed as walls of perimeters × line width all round plus uniform infill; skins equal to the walls".into());
            a.push("print time is the extrusion-limited lower bound: travel, acceleration and cooling only add".into());
        }
        a
    }
    fn run(&self, x: &[f64]) -> Result<Outcome, String> {
        let n = self.params.len();
        let rho = x[n];
        let solid = self.solid(&x[..n])?;
        let fine = cad::mesh_solid(&solid, self.cell)?;
        let coarse = cad::mesh_solid(&solid, 2.0 * self.cell)?;
        let (vf, vc) = (fine.volume(), coarse.volume());
        let (af, ac) = (fine.area(), coarse.area());
        let (lo, hi) = fine.bbox();
        let (clo, chi) = coarse.bbox();
        let size: Vec<f64> = (0..3).map(|k| hi[k] - lo[k]).collect();
        let size_err: Vec<f64> = (0..3).map(|k| ((hi[k] - lo[k]) - (chi[k] - clo[k])).abs()).collect();
        let mut values = vec![vf, af, rho * vf, size[0], size[1], size[2]];
        let mut errors = vec![(vf - vc).abs(), (af - ac).abs(), rho * (vf - vc).abs(), size_err[0], size_err[1], size_err[2]];
        let mut notes = vec![format!("{} triangles at cell {:.3} mm", fine.triangles.len(), self.cell * 1e3)];
        if self.fdm {
            let [nozzle, perimeters, infill, flow] = x[n + 1..n + 5] else { unreachable!() };
            let wall = nozzle * perimeters;
            let pr = crate::fdm::analyse(&solid, vf, self.cell, wall, infill)?;
            // 1.75 mm filament (the common size).
            let filament_area = std::f64::consts::PI * 0.875e-3 * 0.875e-3;
            // Outputs are SI: seconds, displayed in minutes.
            values.extend([rho * pr.printed_volume, pr.printed_volume / flow, pr.printed_volume / filament_area, pr.min_wall]);
            let pe = pr.printed_error + (vf - vc).abs();
            errors.extend([rho * pe, pe / flow, pe / filament_area, pr.min_wall_error]);
            notes.push(format!(
                "walls {:.2} mm ({perimeters} × {:.2} mm); thinnest wall {:.2} ± {:.2} mm",
                wall * 1e3,
                nozzle * 1e3,
                pr.min_wall * 1e3,
                pr.min_wall_error * 1e3
            ));
        }
        Ok(Outcome { values, errors, rigor: Rigor::Asymptotic, notes })
    }
}
