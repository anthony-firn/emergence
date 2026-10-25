//! A chip as a project analysis: Verilog in, timing and power out, as
//! functions of the supply and the junction temperature, so the chip can
//! sit in the same graph as the package, the heat sink and the board.
//!
//! At load the design is parsed, synthesized onto the transistor-level
//! cells and checked equivalent to its RTL (`chip`). Each run analyses the
//! netlist at the requested corner through the characterised grid.

use std::collections::HashMap;
use std::path::Path;
use std::sync::Mutex;

use emergence::Rigor;

use crate::analysis::{Analysis, OutPort, Outcome, Port};
use crate::units::parse_unit;

/// Relative error claimed for the power estimate: on the transistor-level
/// check (`chip/tests/power.rs`) it read 14.8 % low (no glitches, no
/// internal-node switching); 25 % is that check's tolerance.
const POWER_REL: f64 = 0.25;
/// Relative error claimed for leakage: the device law's temperature fold
/// reproduces the TCAD leakage ratio within 3 % at held-out temperatures
/// inside its fit and 11 % at 425 K (`semiconductor/tests/temperature.rs`);
/// each cell takes its worst input state.
const LEAK_REL: f64 = 0.12;

pub struct ChipDesign {
    pub source: String,
    pub module: chip::rtl::Module,
    pub netlist: chip::synth::Netlist,
    pub equivalence: chip::equiv::Equivalence,
    activity: Mutex<HashMap<u64, Vec<f64>>>,
}

fn port(name: &str, unit: &str, doc: &str, monotone: bool) -> Port {
    Port {
        name: name.into(),
        dim: Some(parse_unit(unit).expect("unit").dim),
        unit: unit.into(),
        doc: doc.into(),
        monotone,
        passive: false,
        from: None,
        default: None,
    }
}

fn out(name: &str, unit: &str, doc: &str) -> OutPort {
    OutPort { name: name.into(), dim: parse_unit(unit).expect("unit").dim, unit: unit.into(), doc: doc.into() }
}

impl ChipDesign {
    /// Load `rtl` (relative to `dir`), synthesize and check equivalence.
    pub fn load(dir: &Path, rtl: &str) -> Result<ChipDesign, String> {
        let path = dir.join(rtl);
        let source = std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
        let module = chip::rtl::parse(&source).map_err(|e| format!("{rtl}: {e}"))?;
        let netlist = chip::synth::synthesize(&module).map_err(|e| format!("{rtl}: {e}"))?;
        let equivalence = chip::equiv::check(&module, &netlist, 20_000)?;
        if let Some(m) = &equivalence.mismatch {
            return Err(format!("{rtl}: synthesis does not match the RTL: {m}"));
        }
        Ok(ChipDesign { source, module, netlist, equivalence, activity: Mutex::new(HashMap::new()) })
    }

    fn toggles(&self, activity: f64) -> Vec<f64> {
        let key = activity.to_bits();
        if let Some(t) = self.activity.lock().unwrap().get(&key) {
            return t.clone();
        }
        let t = chip::power::activity(&self.netlist, activity, 4000, 1);
        self.activity.lock().unwrap().insert(key, t.clone());
        t
    }

    /// Identity of the design for result caching.
    pub fn key(&self) -> String {
        let h = self.source.bytes().fold(0xcbf29ce484222325u64, |h, b| (h ^ b as u64).wrapping_mul(0x100000001b3));
        format!("chip:{}:{h:016x}", self.netlist.name)
    }
}

impl Analysis for ChipDesign {
    fn kind(&self) -> &str {
        "chip_design"
    }
    fn summary(&self) -> &str {
        "Verilog synthesized onto transistor-characterised cells: maximum clock, slack and power at a supply and junction temperature"
    }
    fn inputs(&self) -> Vec<Port> {
        vec![
            port("vdd", "V", "supply voltage", false),
            port("temperature", "degC", "junction temperature", false),
            port("frequency", "MHz", "clock frequency", true),
            port("activity", "%", "probability an input bit toggles per cycle", false),
            port("wire_cap", "fF", "wire capacitance per fan-out connection", false),
            port("input_delay", "ps", "input arrival after the clock edge (from upstream flip-flops)", false),
            port("instances", "1", "copies of the design on the die (power and leakage scale with it)", true),
        ]
    }
    fn outputs(&self) -> Vec<OutPort> {
        vec![
            out("fmax", "MHz", "highest clock frequency every setup and output check allows"),
            out("setup_slack", "ps", "worst setup slack at the clock frequency"),
            out("hold_slack", "ps", "worst hold slack (independent of frequency)"),
            out("power", "W", "total power of all instances"),
            out("power_switching", "W", "switching and clock power"),
            out("power_leakage", "W", "leakage power"),
            out("energy_per_cycle", "pJ", "supply energy per clock cycle, one instance"),
            out("transistors", "1", "transistors per instance"),
            out("flops", "1", "flip-flops per instance"),
        ]
    }
    fn assumptions(&self) -> Vec<String> {
        vec![
            "cells: lumped-capacitance compact transistor model of the 150 nm technology (analysis/gates.md)".into(),
            "an ideal clock (no skew, no clock tree); wires as a fixed capacitance per connection".into(),
            "power from settled-value activity: glitches and internal-node switching are excluded (measured −15 % against transistors)".into(),
            "leakage: every cell in its worst input state".into(),
            "all instances at one junction temperature".into(),
        ]
    }
    fn run(&self, x: &[f64]) -> Result<Outcome, String> {
        let [vdd, temp, freq, activity, wire, in_delay, copies] = x[..] else { unreachable!() };
        if freq <= 0.0 || !(0.0..=1.0).contains(&activity) || copies < 1.0 {
            return Err("frequency must be positive, activity within 0–100 %, at least one instance".into());
        }
        let toggles = self.toggles(activity);
        let c = chip::sta::Constraints { period: 1.0 / freq, wire_cap: wire, input_delay: in_delay, ..Default::default() };
        let (m, ind) = chip::corners::at(&self.netlist, &c, &toggles, vdd, temp)?;
        // The interpolation indicator was optimistic by up to 1.43× at four
        // held-out corners (chip/examples/corners.rs): it is doubled.
        let t_err = m.min_period_error + 2.0 * ind.min_period;
        let fmax = 1.0 / m.min_period;
        let fmax_err = t_err / (m.min_period * m.min_period);
        let switching = m.energy_per_cycle * freq * copies;
        let sw_err = switching * POWER_REL + ind.energy_per_cycle * freq * copies;
        let leak = m.leakage * copies;
        let leak_err = leak * (LEAK_REL + ind.leakage);
        // The critical path itself, exact at the nearest grid corner.
        let near = |grid: &[f64; 3], x: f64| *grid.iter().min_by(|a, b| (*a - x).abs().total_cmp(&(*b - x).abs())).unwrap();
        let corner = chip::library::Corner { vdd: near(&chip::corners::VDDS, vdd), temp: near(&chip::corners::TEMPS, temp) };
        let path = chip::corners::exact(&self.netlist, &c, &toggles, corner).map(|(_, t)| {
            let cells: Vec<&str> = t.critical.iter().map(|s| s.cell).collect();
            format!(
                "critical path at {:.2} V, {:.1} K: {} stages to {}: {}",
                corner.vdd,
                corner.temp,
                cells.len(),
                t.endpoint,
                cells.join(" → ")
            )
        });
        let eq = &self.equivalence;
        Ok(Outcome {
            values: vec![fmax, c.period - m.min_period, m.hold_slack, switching + leak, switching, leak, m.energy_per_cycle, self.netlist.transistors() as f64, self.netlist.flops.len() as f64],
            errors: vec![fmax_err, t_err, ind.hold_slack, sw_err + leak_err, sw_err, leak_err, m.energy_per_cycle * POWER_REL + ind.energy_per_cycle, 0.0, 0.0],
            rigor: Rigor::Asymptotic,
            notes: vec![
                format!("{}: {:?}", self.netlist.name, self.netlist.cells()),
                format!(
                    "netlist equivalent to the RTL on {} (input, state) points{}",
                    eq.points,
                    if eq.exhaustive { ", every one: a proof" } else { ", sampled" }
                ),
                format!("corner {vdd:.3} V, {temp:.1} K: interpolation indicator {:.2} ps on the period", ind.min_period * 1e12),
                path.unwrap_or_else(|e| format!("critical path unavailable: {e}")),
            ],
        })
    }
}
