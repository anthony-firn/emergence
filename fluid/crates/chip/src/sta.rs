//! Clocked static timing analysis.
//!
//! Paths start at the flip-flops (clock-to-Q) and the primary inputs (an
//! input delay after the clock edge) and end at the flip-flops (setup and
//! hold) and the primary outputs (an output delay before the next edge).
//! Rise and fall are carried separately at every net, through each cell's
//! unateness: INV, NAND2 and NOR2 invert, the full adder's carry is
//! positive-unate, XOR2 and the sum are non-unate (both input edges reach
//! both output edges). Delays and slews come from the NLDM tables at each
//! stage's own input slew and load, late for setup and early for hold.
//!
//! The certificate on the critical path adds, stage by stage, the table's
//! held-out interpolation error times its delay (`semiconductor::timing`).
//! It covers interpolation of the characterisation, not the lumped
//! capacitance model or wires, which are assumptions.

use semiconductor::gates::CellKind;
use semiconductor::timing::TimingTable;

use crate::library::CellLib;
use crate::synth::Netlist;

/// Edges: 0 = rise, 1 = fall.
const RISE: usize = 0;
const FALL: usize = 1;

#[derive(Clone, Copy, Debug)]
pub struct Constraints {
    /// Clock period (s).
    pub period: f64,
    /// Clock slew at the flip-flops (s).
    pub clock_slew: f64,
    /// Slew of the primary inputs (s).
    pub input_slew: f64,
    /// Primary inputs arrive this long after the clock edge (s).
    pub input_delay: f64,
    /// Primary outputs are needed this long before the next edge (s).
    pub output_delay: f64,
    /// Capacitance on each primary output (F).
    pub output_load: f64,
    /// Wire capacitance per fan-out connection (F).
    pub wire_cap: f64,
}

impl Default for Constraints {
    fn default() -> Self {
        Constraints {
            period: 2e-9,
            clock_slew: 80e-12,
            input_slew: 80e-12,
            input_delay: 0.0,
            output_delay: 0.0,
            output_load: 5e-15,
            wire_cap: 0.5e-15,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Stage {
    pub cell: &'static str,
    pub net: usize,
    pub edge: &'static str,
    pub arrival: f64,
}

#[derive(Clone, Debug)]
pub struct Timing {
    /// Shortest clock period every setup and output check allows (s).
    pub min_period: f64,
    /// Bound on the interpolation error of `min_period` (s).
    pub min_period_error: f64,
    /// Worst setup / output slack at the given period (s); negative fails.
    pub setup_slack: f64,
    /// Worst hold slack (s); negative fails at any clock frequency.
    pub hold_slack: f64,
    /// The critical path, launch to capture.
    pub critical: Vec<Stage>,
    pub endpoint: String,
    /// Load on every net (F).
    pub loads: Vec<f64>,
}

#[derive(Clone, Copy)]
struct Arr {
    late: [f64; 2],
    early: [f64; 2],
    slew_late: [f64; 2],
    slew_early: [f64; 2],
    /// accumulated certificate along the late path, per edge
    cert: [f64; 2],
    /// (driving instance, input pin, input edge) of the late arrival
    from: [Option<(usize, usize, usize)>; 2],
}

const NONE: Arr = Arr {
    late: [f64::NEG_INFINITY; 2],
    early: [f64::INFINITY; 2],
    slew_late: [0.0; 2],
    slew_early: [0.0; 2],
    cert: [0.0; 2],
    from: [None; 2],
};

/// Input edges that can cause output edge `out` through pin `pin`.
fn sources(kind: CellKind, o: usize, out: usize) -> &'static [usize] {
    match (kind, o) {
        (CellKind::Inv | CellKind::Nand2 | CellKind::Nor2, _) => {
            if out == RISE { &[FALL] } else { &[RISE] }
        }
        (CellKind::FullAdder, 1) => {
            if out == RISE { &[RISE] } else { &[FALL] }
        }
        _ => &[RISE, FALL],
    }
}

/// Net loads: the input pins each net drives, plus wire and output loads.
pub fn loads(n: &Netlist, lib: &CellLib, c: &Constraints) -> Vec<f64> {
    let mut load = vec![0.0; n.comb.n_nets];
    for inst in &n.comb.insts {
        let ct = lib.comb.cell(inst.kind);
        for (p, &net) in inst.ins.iter().enumerate() {
            load[net] += ct.c_in[p] + c.wire_cap;
        }
    }
    for f in &n.flops {
        load[f.d] += lib.dff.c_d + c.wire_cap;
    }
    for (_, bits) in &n.outputs {
        for &b in bits {
            load[b] += c.output_load;
        }
    }
    load
}

pub fn analyse(n: &Netlist, lib: &CellLib, c: &Constraints) -> Timing {
    let load = loads(n, lib, c);
    let mut arr = vec![NONE; n.comb.n_nets];
    let launch = |a: &mut Arr, t: [f64; 2], s: [f64; 2]| {
        a.late = t;
        a.early = t;
        a.slew_late = s;
        a.slew_early = s;
    };
    for (_, bits) in &n.inputs {
        for &b in bits {
            launch(&mut arr[b], [c.input_delay; 2], [c.input_slew; 2]);
        }
    }
    let ff = &lib.dff;
    for f in &n.flops {
        let (dr, sr) = ff.clk_q_rise.interp(c.clock_slew, load[f.q]);
        let (df, sf) = ff.clk_q_fall.interp(c.clock_slew, load[f.q]);
        launch(&mut arr[f.q], [dr, df], [sr, sf]);
    }
    // Tie nets never switch: they stay at NONE.
    for (ii, inst) in n.comb.insts.iter().enumerate() {
        let ct = lib.comb.cell(inst.kind);
        for (o, &out) in inst.outs.iter().enumerate() {
            let mut a = NONE;
            for edge in [RISE, FALL] {
                for (p, &inp) in inst.ins.iter().enumerate() {
                    let pin = &ct.outs[o][p];
                    let tab: &TimingTable = if edge == RISE { &pin.rise } else { &pin.fall };
                    for &ie in sources(inst.kind, o, edge) {
                        let src = &arr[inp];
                        if src.late[ie].is_finite() {
                            let (d, s) = tab.interp(src.slew_late[ie], load[out]);
                            let t = src.late[ie] + d;
                            if t > a.late[edge] {
                                a.late[edge] = t;
                                a.slew_late[edge] = s;
                                a.cert[edge] = src.cert[ie] + tab.cert_error * d;
                                a.from[edge] = Some((ii, p, ie));
                            }
                        }
                        if src.early[ie].is_finite() {
                            let (d, s) = tab.interp(src.slew_early[ie], load[out]);
                            let t = src.early[ie] + d;
                            if t < a.early[edge] {
                                a.early[edge] = t;
                                a.slew_early[edge] = s;
                            }
                        }
                    }
                }
            }
            arr[out] = a;
        }
    }
    // Endpoints.
    let mut need = f64::NEG_INFINITY; // shortest period the worst endpoint allows
    let mut need_cert = 0.0;
    let mut worst: Option<(usize, usize, String)> = None;
    let mut hold_slack = f64::INFINITY;
    for f in &n.flops {
        let a = &arr[f.d];
        for (edge, setup, hold) in [(RISE, ff.setup_rise, ff.hold_rise), (FALL, ff.setup_fall, ff.hold_fall)] {
            if a.late[edge].is_finite() && a.late[edge] + setup > need {
                need = a.late[edge] + setup;
                need_cert = a.cert[edge];
                worst = Some((f.d, edge, format!("{} (setup)", f.name)));
            }
            if a.early[edge].is_finite() {
                hold_slack = hold_slack.min(a.early[edge] - hold);
            }
        }
    }
    for (name, bits) in &n.outputs {
        for (k, &b) in bits.iter().enumerate() {
            for edge in [RISE, FALL] {
                let a = &arr[b];
                if a.late[edge].is_finite() && a.late[edge] + c.output_delay > need {
                    need = a.late[edge] + c.output_delay;
                    need_cert = a.cert[edge];
                    worst = Some((b, edge, format!("{name}[{k}] (output)")));
                }
            }
        }
    }
    // Walk the critical path back to its launch.
    let mut critical = Vec::new();
    let mut endpoint = String::new();
    if let Some((mut net, mut edge, ep)) = worst.clone() {
        endpoint = ep;
        loop {
            let a = &arr[net];
            match a.from[edge] {
                Some((ii, p, ie)) => {
                    let inst = &n.comb.insts[ii];
                    critical.push(Stage { cell: inst.kind.name(), net, edge: if edge == RISE { "rise" } else { "fall" }, arrival: a.late[edge] });
                    net = inst.ins[p];
                    edge = ie;
                }
                None => {
                    let what = if n.flops.iter().any(|f| f.q == net) { "DFF clk→Q" } else { "input" };
                    critical.push(Stage { cell: what, net, edge: if edge == RISE { "rise" } else { "fall" }, arrival: a.late[edge] });
                    break;
                }
            }
        }
        critical.reverse();
    }
    let min_period = if need.is_finite() { need } else { 0.0 };
    Timing {
        min_period,
        min_period_error: need_cert,
        setup_slack: c.period - min_period,
        hold_slack,
        critical,
        endpoint,
        loads: load,
    }
}
