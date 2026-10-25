//! HVAC analyses for the project graph (`hvac` crate): zone loads, the
//! refrigeration cycle, duct runs, fans, cooling coils and comfort.

use std::collections::BTreeMap;
use std::path::Path;

use emergence::Rigor;
use serde_json::Value as Json;

use crate::analysis::{Analysis, OutPort, Outcome, Port};
use crate::units::parse_unit;

fn port(name: &str, unit: &str, doc: &str) -> Port {
    Port {
        name: name.into(),
        dim: Some(parse_unit(unit).expect("unit").dim),
        unit: unit.into(),
        doc: doc.into(),
        monotone: false,
        passive: false,
        from: None,
        default: None,
    }
}

fn mono(mut p: Port) -> Port {
    p.monotone = true;
    p
}

fn out(name: &str, unit: &str, doc: &str) -> OutPort {
    OutPort { name: name.into(), dim: parse_unit(unit).expect("unit").dim, unit: unit.into(), doc: doc.into() }
}

const C2K: f64 = 273.15;

// ---------------------------------------------------------------------------
// Zone loads
// ---------------------------------------------------------------------------

/// The unit a zone-file field takes when it names an input.
fn field_unit(field: &str) -> &'static str {
    match field {
        "thickness" => "mm",
        "k" => "W/(m*K)",
        "rho" => "kg/m^3",
        "cp" => "J/(kg*K)",
        "area" => "m^2",
        "u" => "W/(m^2*K)",
        "shgc" | "absorptance" | "people" => "1",
        "person_sensible" | "person_latent" | "lighting" | "equipment" => "W",
        "outdoor_air" => "L/s",
        _ => "1",
    }
}

/// A zone file: numbers, or strings naming analysis inputs.
pub struct ZoneLoads {
    spec: Json,
    /// Inputs named in the file, with their units, in port order after the
    /// fixed climate inputs.
    named: Vec<(String, &'static str)>,
    key: String,
}

const FIXED: &[(&str, &str, &str)] = &[
    ("design_dry_bulb", "degC", "summer design dry bulb (ASHRAE 0.4 %)"),
    ("daily_range", "K", "summer design daily range"),
    ("coincident_wet_bulb", "degC", "mean coincident wet bulb"),
    ("winter_design", "degC", "winter design dry bulb (ASHRAE 99.6 %)"),
    ("cooling_setpoint", "degC", "room temperature held in cooling"),
    ("cooling_rh", "%", "room relative humidity held in cooling"),
    ("heating_setpoint", "degC", "room temperature held in heating"),
    ("supply_delta_t", "K", "supply air below the room in cooling"),
];

fn collect(v: &Json, field: &str, out: &mut Vec<(String, &'static str)>) -> Result<(), String> {
    match v {
        Json::String(s) => {
            let u = field_unit(field);
            if let Some((_, prev)) = out.iter().find(|(n, _)| n == s) {
                if *prev != u {
                    return Err(format!("input `{s}` is used as both {prev} and {u}"));
                }
            } else {
                out.push((s.clone(), u));
            }
        }
        Json::Array(a) => {
            for x in a {
                collect(x, field, out)?;
            }
        }
        Json::Object(o) => {
            for (k, x) in o {
                if k != "name" && k != "assembly" && k != "flow" {
                    collect(x, k, out)?;
                }
            }
        }
        _ => {}
    }
    Ok(())
}

impl ZoneLoads {
    pub fn load(dir: &Path, file: &str) -> Result<ZoneLoads, String> {
        let path = dir.join(file);
        let text = std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
        let spec: Json = serde_json::from_str(&text).map_err(|e| format!("{file}: {e}"))?;
        let mut named = Vec::new();
        for key in ["assemblies", "opaque", "windows", "gains", "outdoor_air"] {
            if let Some(v) = spec.get(key) {
                collect(v, key, &mut named)?;
            }
        }
        named.retain(|(n, _)| !FIXED.iter().any(|f| f.0 == n));
        let h = text.bytes().fold(0xcbf29ce484222325u64, |h, b| (h ^ b as u64).wrapping_mul(0x100000001b3));
        let me = ZoneLoads { spec, named, key: format!("zone:{file}:{h:016x}") };
        // Build once with placeholder inputs to catch structural mistakes.
        let x = vec![1.0; FIXED.len() + me.named.len()];
        me.build(&x)?;
        Ok(me)
    }
    pub fn key(&self) -> String {
        self.key.clone()
    }

    fn num(&self, v: &Json, x: &BTreeMap<&str, f64>) -> Result<f64, String> {
        match v {
            Json::Number(n) => Ok(n.as_f64().unwrap()),
            Json::String(s) => x.get(s.as_str()).copied().ok_or_else(|| format!("no input `{s}`")),
            other => Err(format!("expected a number or an input name, found {other}")),
        }
    }

    fn build(&self, xs: &[f64]) -> Result<(hvac::zone::Zone, hvac::solar::Site, f64), String> {
        use hvac::envelope::{Assembly, Flow, Layer};
        use hvac::zone::*;
        let x: BTreeMap<&str, f64> = self.named.iter().enumerate().map(|(i, (n, _))| (n.as_str(), xs[FIXED.len() + i])).collect();
        let s = &self.spec;
        fn get<'a>(v: &'a Json, k: &str) -> Result<&'a Json, String> {
            v.get(k).ok_or_else(|| format!("zone file: missing `{k}`"))
        }
        let site_j = get(s, "site")?;
        let site = hvac::solar::Site {
            latitude: self.num(get(site_j, "latitude")?, &x)?,
            longitude: self.num(get(site_j, "longitude")?, &x)?,
            meridian: self.num(get(site_j, "meridian")?, &x)?,
            tau_b: self.num(get(site_j, "tau_b")?, &x)?,
            tau_d: self.num(get(site_j, "tau_d")?, &x)?,
            ground_reflectance: self.num(get(site_j, "ground_reflectance")?, &x)?,
        };
        let day = self.num(get(site_j, "day_of_year")?, &x)?;
        let mut assemblies = BTreeMap::new();
        for (name, a) in get(s, "assemblies")?.as_object().ok_or("assemblies must be an object")? {
            let flow = match a.get("flow").and_then(|f| f.as_str()).unwrap_or("horizontal") {
                "horizontal" => Flow::Horizontal,
                "upward" => Flow::Upward,
                "downward" => Flow::Downward,
                f => return Err(format!("assembly `{name}`: flow `{f}`")),
            };
            let mut layers = Vec::new();
            for l in get(a, "layers")?.as_array().ok_or("layers must be a list")? {
                // Thickness in mm in the file (as its input unit is).
                layers.push(Layer {
                    thickness: self.num(get(l, "thickness")?, &x)? * if matches!(get(l, "thickness")?, Json::String(_)) { 1.0 } else { 1e-3 },
                    k: self.num(get(l, "k")?, &x)?,
                    rho: self.num(get(l, "rho")?, &x)?,
                    cp: self.num(get(l, "cp")?, &x)?,
                });
            }
            assemblies.insert(name.clone(), Assembly::new(layers, flow));
        }
        let mut opaque = Vec::new();
        for o in get(s, "opaque")?.as_array().ok_or("opaque must be a list")? {
            let an = get(o, "assembly")?.as_str().ok_or("assembly must be a name")?;
            let assembly = assemblies.get(an).ok_or_else(|| format!("no assembly `{an}`"))?.clone();
            let internal = o.get("internal").and_then(|v| v.as_bool()).unwrap_or(false);
            let exterior = if internal {
                None
            } else {
                Some((self.num(get(o, "azimuth")?, &x)?, self.num(get(o, "tilt")?, &x)?, self.num(get(o, "absorptance")?, &x)?))
            };
            opaque.push(Opaque { name: o.get("name").and_then(|v| v.as_str()).unwrap_or("").into(), assembly, area: self.num(get(o, "area")?, &x)?, exterior });
        }
        let mut windows = Vec::new();
        for w in s.get("windows").and_then(|v| v.as_array()).cloned().unwrap_or_default() {
            windows.push(Window {
                name: w.get("name").and_then(|v| v.as_str()).unwrap_or("").into(),
                area: self.num(get(&w, "area")?, &x)?,
                u: self.num(get(&w, "u")?, &x)?,
                shgc: self.num(get(&w, "shgc")?, &x)?,
                azimuth: self.num(get(&w, "azimuth")?, &x)?,
                tilt: self.num(get(&w, "tilt")?, &x)?,
            });
        }
        let g = get(s, "gains")?;
        let rad = get(g, "radiant")?.as_array().ok_or("radiant must be [people, lighting, equipment]")?;
        let sch = get(g, "schedule")?.as_array().ok_or("schedule must be [start, end] hours")?;
        let gains = Gains {
            people: self.num(get(g, "people")?, &x)?,
            person_sensible: self.num(get(g, "person_sensible")?, &x)?,
            person_latent: self.num(get(g, "person_latent")?, &x)?,
            lighting: self.num(get(g, "lighting")?, &x)?,
            equipment: self.num(get(g, "equipment")?, &x)?,
            radiant: (self.num(&rad[0], &x)?, self.num(&rad[1], &x)?, self.num(&rad[2], &x)?),
            schedule: (self.num(&sch[0], &x)?, self.num(&sch[1], &x)?),
        };
        // Outdoor air: a named input is in SI (m³/s); a number in the file
        // is in L/s.
        let oa_j = get(s, "outdoor_air")?;
        let oa = self.num(oa_j, &x)? * if matches!(oa_j, Json::String(_)) { 1.0 } else { 1e-3 };
        Ok((
            Zone { floor_area: self.num(get(s, "floor_area")?, &x)?, volume: self.num(get(s, "volume")?, &x)?, opaque, windows, gains, outdoor_air: oa },
            site,
            day,
        ))
    }
}

impl Analysis for ZoneLoads {
    fn kind(&self) -> &str {
        "zone_loads"
    }
    fn summary(&self) -> &str {
        "Design cooling (sensible, latent) and heating loads of a room by heat balance over the ASHRAE design day"
    }
    fn inputs(&self) -> Vec<Port> {
        let mut v: Vec<Port> = FIXED.iter().map(|(n, u, d)| port(n, u, d)).collect();
        v.extend(self.named.iter().map(|(n, u)| port(n, u, "named in the zone file")));
        v
    }
    fn outputs(&self) -> Vec<OutPort> {
        vec![
            out("cooling_sensible", "W", "peak sensible cooling load"),
            out("cooling_latent", "W", "latent load at that hour"),
            out("cooling_total", "W", "sensible + latent at the peak"),
            out("peak_hour", "h", "hour of the sensible peak (solar time of day)"),
            out("heating", "W", "design heating load"),
            out("supply_airflow", "L/s", "supply air to meet the sensible peak"),
            out("ua", "W/K", "envelope conductance"),
            out("solar_peak", "W", "peak solar transmitted through the windows"),
        ]
    }
    fn assumptions(&self) -> Vec<String> {
        vec![
            "ASHRAE clear-sky design day, repeated to a periodic response; outdoor humidity constant over the day".into(),
            "combined convective-radiative surface coefficients (ISO 6946); radiant gains spread over the room's surfaces by area".into(),
            "windows: U and SHGC at normal incidence, no shading; no thermal bridges".into(),
            "internal elements between identical rooms (adiabatic mid-plane)".into(),
            "heating: steady at the winter design temperature, no sun or internal gains".into(),
        ]
    }
    fn run(&self, xs: &[f64]) -> Result<Outcome, String> {
        let (zone, site, day) = self.build(xs)?;
        let [tmax, range, twb, twin, tset, rh, theat, dts] = xs[..FIXED.len()] else { unreachable!() };
        let (tmax, twb, twin, tset, theat) = (tmax - C2K, twb - C2K, twin - C2K, tset - C2K, theat - C2K);
        let p = hvac::psychro::P_ATM;
        let w_out = hvac::psychro::w_from_wet_bulb(tmax, twb, p);
        let w_in = hvac::psychro::w_from_rh(tset, rh, p);
        let dd = hvac::zone::DesignDay { day_of_year: day, t_max: tmax, range, w_out, pressure: p };
        // Fine and coarse time steps: the difference is the error estimate.
        let fine = zone.cooling_load(&site, &dd, tset, w_in, 4);
        let coarse = zone.cooling_load(&site, &dd, tset, w_in, 2);
        let sens = fine.peak_sensible;
        let err = (fine.peak_sensible - coarse.peak_sensible).abs();
        let heating = zone.heating_load(theat, twin);
        let solar = fine.hours.iter().map(|h| h.solar_transmitted).fold(0.0, f64::max);
        let flow = hvac::zone::supply_airflow(sens.max(0.0), dts);
        Ok(Outcome {
            values: vec![sens, fine.peak_latent, sens + fine.peak_latent, fine.peak_hour * 3600.0, heating, flow, zone.ua(), solar],
            errors: vec![err, 0.0, err, 0.0, 0.0, flow * err / sens.abs().max(1e-9), 0.0, 0.0],
            rigor: Rigor::Asymptotic,
            notes: vec![
                format!("outdoor air {:.1} L/s, humidity ratio {:.4} outdoors, {:.4} in the room", zone.outdoor_air * 1e3, w_out, w_in),
                format!("periodic after {} design days; time-step estimate {err:.1} W", fine.days),
            ],
        })
    }
}

// ---------------------------------------------------------------------------
// Refrigeration cycle
// ---------------------------------------------------------------------------

pub struct Refrigeration {
    eos: compressible::helmholtz::Helmholtz,
}

impl Refrigeration {
    pub fn new(name: &str) -> Result<Refrigeration, String> {
        Ok(Refrigeration { eos: hvac::refrigeration::refrigerant(name)? })
    }
}

impl Analysis for Refrigeration {
    fn kind(&self) -> &str {
        "refrigeration_cycle"
    }
    fn summary(&self) -> &str {
        "Single-stage vapour-compression cycle on a reference refrigerant equation of state"
    }
    fn inputs(&self) -> Vec<Port> {
        vec![
            port("evaporating", "degC", "evaporating temperature"),
            port("condensing", "degC", "condensing temperature"),
            port("superheat", "K", "superheat at the evaporator exit"),
            port("subcooling", "K", "subcooling at the condenser exit"),
            port("isentropic_efficiency", "1", "compressor isentropic efficiency"),
            mono(port("capacity", "W", "cooling capacity to deliver")),
        ]
    }
    fn outputs(&self) -> Vec<OutPort> {
        vec![
            out("cop_cooling", "1", "cooling COP"),
            out("cop_heating", "1", "heating COP"),
            out("compressor_power", "W", "compressor shaft power at the capacity"),
            out("heat_rejected", "W", "condenser heat"),
            out("mass_flow", "g/s", "refrigerant mass flow"),
            out("discharge_temperature", "degC", "compressor discharge temperature"),
            out("p_evaporating", "kPa", "evaporating pressure"),
            out("p_condensing", "kPa", "condensing pressure"),
        ]
    }
    fn assumptions(&self) -> Vec<String> {
        vec!["no pressure drops in exchangers or lines; adiabatic compressor; single stage, subcritical".into()]
    }
    fn run(&self, x: &[f64]) -> Result<Outcome, String> {
        let [te, tc, sh, sc, eta, cap] = x[..] else { unreachable!() };
        let c = hvac::refrigeration::cycle(&self.eos, te, tc, sh, sc, eta)?;
        let m = cap / c.q_evap;
        let w = m * c.w_comp;
        Ok(Outcome {
            values: vec![c.cop_cooling, c.cop_heating, w, cap + w, m, c.t_discharge, c.p_evap, c.p_cond],
            errors: vec![0.0; 8],
            rigor: Rigor::Validated,
            notes: vec![format!("{}: expansion leaves quality {:.3}; Carnot COP {:.2}", self.eos.name, c.quality_in, te / (tc - te))],
        })
    }
}

// ---------------------------------------------------------------------------
// Duct run, fan, coil, comfort
// ---------------------------------------------------------------------------

pub struct DuctRun;

impl Analysis for DuctRun {
    fn kind(&self) -> &str {
        "duct_run"
    }
    fn summary(&self) -> &str {
        "Pressure drop of a duct or pipe run (Darcy–Weisbach, Colebrook); height 0 means round of diameter `width`"
    }
    fn inputs(&self) -> Vec<Port> {
        let mut t = port("fluid_temperature", "degC", "fluid temperature (for properties)");
        t.passive = true;
        vec![
            mono(port("flow", "L/s", "volume flow")),
            port("width", "mm", "width, or diameter for a round duct"),
            port("height", "mm", "height (0 for round)"),
            mono(port("length", "m", "straight length")),
            mono(port("roughness", "mm", "absolute roughness (galvanised steel 0.09 mm, ASHRAE)")),
            mono(port("k_sum", "1", "sum of fitting loss coefficients")),
            t,
            Port { from: Some(crate::analysis::MatRef { slot: "fluid".into(), prop: "rho".into(), at: Some("fluid_temperature".into()) }), ..port("rho", "kg/m^3", "fluid density") },
            Port { from: Some(crate::analysis::MatRef { slot: "fluid".into(), prop: "mu".into(), at: Some("fluid_temperature".into()) }), ..port("mu", "Pa*s", "fluid viscosity") },
        ]
    }
    fn outputs(&self) -> Vec<OutPort> {
        vec![
            out("pressure_drop", "Pa", "total pressure drop"),
            out("velocity", "m/s", "mean velocity"),
            out("reynolds", "1", "Reynolds number"),
            out("friction_factor", "1", "Darcy friction factor"),
        ]
    }
    fn assumptions(&self) -> Vec<String> {
        vec!["fully developed friction plus fitting coefficients; transitional flow (Re 2300–4000) given the turbulent factor".into()]
    }
    fn run(&self, x: &[f64]) -> Result<Outcome, String> {
        let [q, w, h, l, eps, k, _t, rho, mu] = x[..] else { unreachable!() };
        let run = if h <= 0.0 { hvac::hydraulics::Run::round(w, l, eps, k) } else { hvac::hydraulics::Run::rectangular(w, h, l, eps, k) };
        let d = run.drop(q, rho, mu);
        let note = match d.regime {
            hvac::hydraulics::Regime::Transitional => "transitional flow: the turbulent friction factor is used (conservative)".to_string(),
            r => format!("{r:?} flow"),
        };
        Ok(Outcome { values: vec![d.dp, d.velocity, d.reynolds, d.friction], errors: vec![0.0; 4], rigor: Rigor::Validated, notes: vec![note] })
    }
}

/// A fan meeting a quadratic system curve through its design point.
pub struct Fan {
    curve: hvac::hydraulics::Curve,
}

impl Fan {
    /// `curve`: three "flow pressure efficiency" triples in L/s, Pa and
    /// fraction, separated by `;`.
    pub fn new(curve: &str) -> Result<Fan, String> {
        let pts: Vec<(f64, f64, f64)> = curve
            .split(';')
            .map(|p| {
                let v: Vec<f64> = p.split_whitespace().map(|x| x.parse::<f64>().map_err(|_| format!("bad number `{x}` in the fan curve"))).collect::<Result<_, _>>()?;
                if v.len() != 3 {
                    return Err("each fan-curve point is `flow_L/s pressure_Pa efficiency`".to_string());
                }
                Ok((v[0] * 1e-3, v[1], v[2]))
            })
            .collect::<Result<_, _>>()?;
        if pts.len() != 3 {
            return Err("the fan curve needs three points".into());
        }
        Ok(Fan { curve: hvac::hydraulics::Curve::through([pts[0], pts[1], pts[2]]) })
    }
}

impl Analysis for Fan {
    fn kind(&self) -> &str {
        "fan"
    }
    fn summary(&self) -> &str {
        "Fan operating point on a quadratic system curve through the system's design point"
    }
    fn inputs(&self) -> Vec<Port> {
        vec![port("design_flow", "L/s", "system design flow"), port("design_pressure", "Pa", "system pressure drop at the design flow")]
    }
    fn outputs(&self) -> Vec<OutPort> {
        vec![
            out("flow", "L/s", "operating flow"),
            out("pressure", "Pa", "fan pressure rise"),
            out("efficiency", "1", "fan efficiency at the operating point"),
            out("power", "W", "shaft power"),
        ]
    }
    fn assumptions(&self) -> Vec<String> {
        vec!["system pressure drop scales as flow² (turbulent, fixed geometry)".into()]
    }
    fn run(&self, x: &[f64]) -> Result<Outcome, String> {
        let [qd, pd] = x[..] else { unreachable!() };
        let k = pd / (qd * qd);
        let qmax = 10.0 * qd;
        let op = hvac::hydraulics::operating_point(&self.curve, &|q| k * q * q, qmax).ok_or("the fan cannot drive this system (its curve does not cross the system curve)")?;
        Ok(Outcome { values: vec![op.flow, op.pressure, op.efficiency, op.power], errors: vec![0.0; 4], rigor: Rigor::Validated, notes: vec![] })
    }
}

pub struct Coil;

impl Analysis for Coil {
    fn kind(&self) -> &str {
        "cooling_coil"
    }
    fn summary(&self) -> &str {
        "Cooling and dehumidifying coil by apparatus dew point and bypass factor"
    }
    fn inputs(&self) -> Vec<Port> {
        vec![
            mono(port("airflow", "L/s", "air volume flow (at 1.204 kg/m³)")),
            port("entering", "degC", "entering dry bulb"),
            port("entering_rh", "%", "entering relative humidity"),
            port("adp", "degC", "apparatus dew point"),
            port("bypass", "1", "bypass factor"),
        ]
    }
    fn outputs(&self) -> Vec<OutPort> {
        vec![
            out("total", "W", "total capacity"),
            out("sensible", "W", "sensible capacity"),
            out("latent", "W", "latent capacity"),
            out("leaving", "degC", "leaving dry bulb"),
            out("leaving_rh", "%", "leaving relative humidity"),
            out("condensate", "g/s", "condensate"),
        ]
    }
    fn assumptions(&self) -> Vec<String> {
        vec!["ADP/bypass-factor coil model; dry when the ADP is above the entering dew point".into()]
    }
    fn run(&self, x: &[f64]) -> Result<Outcome, String> {
        let [q, t, rh, adp, bf] = x[..] else { unreachable!() };
        let p = hvac::psychro::P_ATM;
        let (t, adp) = (t - C2K, adp - C2K);
        let w = hvac::psychro::w_from_rh(t, rh, p);
        let c = hvac::coil::cooling_coil(t, w, adp, bf, 1.204 * q, p);
        let rh_out = hvac::psychro::rh(c.t_out, c.w_out, p);
        Ok(Outcome {
            values: vec![c.total, c.sensible, c.latent, c.t_out + C2K, rh_out, c.condensate],
            errors: vec![0.0; 6],
            rigor: Rigor::Validated,
            notes: vec![format!("{} coil, SHR {:.2}", if c.wet { "wet" } else { "dry" }, c.sensible / c.total)],
        })
    }
}

pub struct Comfort;

impl Analysis for Comfort {
    fn kind(&self) -> &str {
        "thermal_comfort"
    }
    fn summary(&self) -> &str {
        "Predicted mean vote and percentage dissatisfied (ISO 7730 / ASHRAE 55)"
    }
    fn inputs(&self) -> Vec<Port> {
        vec![
            port("air_temperature", "degC", "air temperature"),
            port("radiant_temperature", "degC", "mean radiant temperature"),
            port("air_speed", "m/s", "relative air speed"),
            port("rh", "%", "relative humidity"),
            port("met", "1", "metabolic rate (met)"),
            port("clo", "1", "clothing insulation (clo)"),
        ]
    }
    fn outputs(&self) -> Vec<OutPort> {
        vec![out("pmv", "1", "predicted mean vote"), out("ppd", "%", "predicted percentage dissatisfied"), out("operative", "degC", "operative temperature")]
    }
    fn assumptions(&self) -> Vec<String> {
        vec!["steady state, uniform environment (no draught or asymmetry)".into()]
    }
    fn run(&self, x: &[f64]) -> Result<Outcome, String> {
        let [ta, tr, v, rh, met, clo] = x[..] else { unreachable!() };
        let (ta, tr) = (ta - C2K, tr - C2K);
        let (pmv, ppd) = hvac::comfort::pmv_ppd(ta, tr, v, rh, met, clo, 0.0);
        Ok(Outcome {
            values: vec![pmv, ppd / 100.0, hvac::comfort::operative(ta, tr, v) + C2K],
            errors: vec![0.0; 3],
            rigor: Rigor::Validated,
            notes: vec![],
        })
    }
}
