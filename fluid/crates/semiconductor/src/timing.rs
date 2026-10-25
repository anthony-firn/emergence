//! Gate-level timing rung: non-linear delay model (NLDM) interpolation, an
//! event-driven logic simulator, static timing analysis, and the error
//! certificate.
//!
//! The coarse model is a per-cell table: propagation delay and output slew
//! as functions of input slew and load capacitance, for rise and fall. The
//! validity indicator is the interpolation error, measured on held-out
//! characterisation points and accumulated along the critical path. Where
//! input slews or loads leave the characterised box, or a path runs through
//! a cell that was not characterised, the netlist unfolds to the
//! transistor-level simulator in `gates.rs`.

use crate::gates::{CellKind, CellTiming};

/// One non-linear delay model: a table on a rectangular (input slew, load)
/// grid, interpolated bilinearly (clamped at the edges).
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct TimingTable {
    /// input slew axis (s), 10 %–90 % of V_dd
    pub slew: Vec<f64>,
    /// load capacitance axis (F)
    pub load: Vec<f64>,
    /// [is][il]
    pub delay: Vec<Vec<f64>>,
    pub out_slew: Vec<Vec<f64>>,
    /// relative delay error certificate from held-out characterisation
    pub cert_error: f64,
    /// absolute output-slew error certificate (s)
    pub cert_slew: f64,
}

impl TimingTable {
    pub fn new(slew: Vec<f64>, load: Vec<f64>, delay: Vec<Vec<f64>>, out_slew: Vec<Vec<f64>>) -> Self {
        TimingTable { slew, load, delay, out_slew, cert_error: 0.0, cert_slew: 0.0 }
    }

    /// Constant table (a single point, no slew/load dependence). Useful for
    /// unit tests and for idealised comparisons.
    pub fn constant(delay: f64, out_slew: f64) -> Self {
        TimingTable {
            slew: vec![0.0],
            load: vec![0.0],
            delay: vec![vec![delay]],
            out_slew: vec![vec![out_slew]],
            cert_error: 0.0,
            cert_slew: 0.0,
        }
    }

    fn interp_axis(x: &[f64], t: f64) -> (usize, f64) {
        let n = x.len();
        if n == 1 {
            return (0, 0.0);
        }
        let t = t.clamp(x[0], x[n - 1]);
        let mut i = 0;
        while i + 2 < n && x[i + 1] < t {
            i += 1;
        }
        let h = x[i + 1] - x[i];
        (i, if h > 0.0 { (t - x[i]) / h } else { 0.0 })
    }

    /// Bilinear interpolation of (delay, output slew) at (input slew, load).
    pub fn interp(&self, s: f64, l: f64) -> (f64, f64) {
        let (i, fs) = Self::interp_axis(&self.slew, s);
        let (j, fl) = Self::interp_axis(&self.load, l);
        let (i1, j1) = ((i + 1).min(self.slew.len() - 1), (j + 1).min(self.load.len() - 1));
        let b = |tab: &Vec<Vec<f64>>| {
            let v00 = tab[i][j];
            let v01 = tab[i][j1];
            let v10 = tab[i1][j];
            let v11 = tab[i1][j1];
            (v00 * (1.0 - fs) + v10 * fs) * (1.0 - fl) + (v01 * (1.0 - fs) + v11 * fs) * fl
        };
        (b(&self.delay), b(&self.out_slew))
    }

    /// Largest relative delay error against held-out points
    /// `(input slew, load, transistor-level delay, output slew)`; also
    /// returns the largest absolute slew error.
    pub fn certify(&self, pts: &[(f64, f64, f64, f64)]) -> (f64, f64) {
        let mut ere = 0.0f64;
        let mut ese = 0.0f64;
        for &(s, l, d_true, sl_true) in pts {
            let (d, sl) = self.interp(s, l);
            if d_true.abs() > 0.0 {
                ere = ere.max(((d - d_true) / d_true).abs());
            }
            ese = ese.max((sl - sl_true).abs());
        }
        (ere, ese)
    }

    /// Is (input slew, load) inside the characterised box?
    pub fn contains(&self, s: f64, l: f64) -> bool {
        let (s0, s1) = (self.slew[0], *self.slew.last().unwrap());
        let (l0, l1) = (self.load[0], *self.load.last().unwrap());
        s >= s0 - 1e-18 && s <= s1 + 1e-18 && l >= l0 - 1e-18 && l <= l1 + 1e-18
    }
}

/// The gate library: characterised timing for every standard cell, tagged by
/// the characterisation box so the simulator knows when to unfold.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Library {
    pub cells: Vec<CellTiming>,
    /// VLDM voltage (V)
    pub vdd: f64,
}

impl Library {
    pub fn cell(&self, k: CellKind) -> &CellTiming {
        &self.cells[k as usize]
    }

    /// External (fan-out) load of every net used for timing: the driven
    /// input capacitances plus a wire estimate. The driver's own
    /// drain-diffusion capacitance is already part of the characterisation,
    /// so it must not be added here.
    pub fn loads(&self, net: &GateNet, wire_cap: f64) -> Vec<f64> {
        let mut load = vec![wire_cap; net.n_nets];
        for inst in &net.insts {
            let ct = self.cell(inst.kind);
            for (p, &n) in inst.ins.iter().enumerate() {
                load[n] += ct.c_in[p];
            }
        }
        load
    }

    /// Driver's own output diffusion capacitance of every net (F); zero for
    /// primary inputs.
    pub fn driver_caps(&self, net: &GateNet) -> Vec<f64> {
        let mut c = vec![0.0; net.n_nets];
        for inst in &net.insts {
            let ct = self.cell(inst.kind);
            for (o, &n) in inst.outs.iter().enumerate() {
                c[n] = ct.c_out[o];
            }
        }
        c
    }
}

/// A gate-level instance. `outs.len()` is 1 except for the full adder (2).
#[derive(Clone, Debug)]
pub struct Inst {
    pub kind: CellKind,
    pub ins: Vec<usize>,
    pub outs: Vec<usize>,
}

/// A gate-level netlist. Instances must be in topological order (all inputs
/// produced earlier).
#[derive(Clone, Debug)]
pub struct GateNet {
    pub n_nets: usize,
    pub inputs: Vec<usize>,
    pub insts: Vec<Inst>,
}

impl GateNet {
    /// Evaluate the logic of the whole netlist for a primary-input vector.
    pub fn eval(&self, in_vec: &[bool]) -> Vec<bool> {
        let mut v = vec![false; self.n_nets];
        for (i, &n) in self.inputs.iter().enumerate() {
            v[n] = in_vec[i];
        }
        for inst in &self.insts {
            let iv: Vec<bool> = inst.ins.iter().map(|&n| v[n]).collect();
            let ov = inst.kind.eval(&iv);
            for (o, val) in ov.iter().enumerate() {
                v[inst.outs[o]] = *val;
            }
        }
        v
    }
}

/// One event on a net: time of the 50 % crossing, new logic value, slew (s).
#[derive(Clone, Copy, Debug)]
struct Ev {
    t: f64,
    net: usize,
    val: bool,
    slew: f64,
}

/// Result of a gate-level event simulation.
#[derive(Clone, Debug)]
pub struct SimResult {
    pub value: Vec<bool>,
    pub arrival: Vec<f64>,
    pub slew: Vec<f64>,
    /// total dynamic energy (J)
    pub energy: f64,
    /// total leakage energy over the simulated window (J)
    pub leakage_energy: f64,
    /// (time, net, value) of every scheduled output transition
    pub events: Vec<(f64, usize, bool)>,
}

impl SimResult {
    /// Arrival time of the last transition on `net`.
    pub fn last(&self, net: usize) -> f64 {
        self.events.iter().filter(|(_, n, _)| *n == net).map(|(t, _, _)| *t).fold(f64::NAN, f64::max)
    }
}

/// Event-driven gate-level simulation. `stim` is a list of (time, input
/// vector) steps; each changed primary input enters with slew `in_slew`.
/// The simulation runs to `t_end` or until the event queue empties.
pub fn simulate(lib: &Library, net: &GateNet, stim: &[(f64, Vec<bool>)], in_slew: f64, t_end: f64) -> SimResult {
    let load = lib.loads(net, 0.2e-15);
    let v0 = &stim[0].1;
    let mut value = net.eval(v0);
    let mut arrival = vec![0.0; net.n_nets];
    let mut slew = vec![in_slew; net.n_nets];
    for &n in &net.inputs {
        arrival[n] = stim[0].0;
    }
    let mut energy = 0.0;
    let mut events = Vec::new();
    let mut queue: Vec<Ev> = Vec::new();
    // Seed transitions from the second stimulus step onward.
    for w in stim.windows(2) {
        let (t, ref to) = w[1];
        let ref from = w[0].1;
        for (i, &n) in net.inputs.iter().enumerate() {
            if from[i] != to[i] {
                queue.push(Ev { t, net: n, val: to[i], slew: in_slew });
            }
        }
    }
    queue.sort_by(|a, b| b.t.partial_cmp(&a.t).unwrap());
    let mut processed = 0usize;
    while let Some(ev) = queue.pop() {
        if ev.t > t_end {
            break;
        }
        if value[ev.net] == ev.val {
            continue;
        }
        value[ev.net] = ev.val;
        arrival[ev.net] = ev.t;
        slew[ev.net] = ev.slew;
        events.push((ev.t, ev.net, ev.val));
        processed += 1;
        for inst in &net.insts {
            let pin = match inst.ins.iter().position(|&n| n == ev.net) {
                Some(p) => p,
                None => continue,
            };
            let iv: Vec<bool> = inst.ins.iter().map(|&n| value[n]).collect();
            let ov = inst.kind.eval(&iv);
            for (o, &newout) in ov.iter().enumerate() {
                let out = inst.outs[o];
                if newout == value[out] {
                    continue;
                }
                let ct = lib.cell(inst.kind);
                let pt = &ct.outs[o][pin];
                let tab = if newout { &pt.rise } else { &pt.fall };
                let (d, sl) = tab.interp(ev.slew, load[out]);
                let etab = if newout { &pt.e_rise } else { &pt.e_fall };
                energy += etab.interp(ev.slew, load[out]).0;
                let tev = Ev { t: ev.t + d, net: out, val: newout, slew: sl };
                // insert keeping descending-time order
                let pos = queue.partition_point(|x| x.t > tev.t);
                queue.insert(pos, tev);
            }
        }
        if processed > 1_000_000 {
            break;
        }
    }
    let dt = t_end - stim[0].0;
    let leakage_energy: f64 = net
        .insts
        .iter()
        .map(|inst| lib.cell(inst.kind).i_leak * lib.vdd * dt)
        .sum();
    SimResult { value, arrival, slew, energy, leakage_energy, events }
}

/// Static timing analysis. Returns (arrival, slew) per net and the critical
/// path (list of nets from a primary input to the worst output), using the
/// worst (late) delay of each gate. Instances must be topologically ordered.
pub fn sta(lib: &Library, net: &GateNet, net_loads: &[f64], in_slew: f64) -> (Vec<f64>, Vec<f64>, Vec<usize>) {
    let mut arrival = vec![0.0f64; net.n_nets];
    let mut slew = vec![0.0f64; net.n_nets];
    let mut pred = vec![usize::MAX; net.n_nets];
    for &n in &net.inputs {
        arrival[n] = 0.0;
        slew[n] = in_slew;
    }
    for inst in &net.insts {
        let ct = lib.cell(inst.kind);
        for (o, &out) in inst.outs.iter().enumerate() {
            let mut best = f64::NEG_INFINITY;
            let mut best_pin = 0usize;
            let mut best_slew = 0.0;
            for (p, &inn) in inst.ins.iter().enumerate() {
                let (d, _) = ct.outs[o][p].rise.interp(slew[inn], net_loads[out]);
                let a = arrival[inn] + d;
                let sl = ct.outs[o][p].rise.interp(slew[inn], net_loads[out]).1;
                if a > best {
                    best = a;
                    best_pin = p;
                    best_slew = sl;
                }
            }
            arrival[out] = best;
            slew[out] = best_slew;
            pred[out] = if best_pin < inst.ins.len() { inst.ins[best_pin] } else { usize::MAX };
        }
    }
    // critical path: worst output net, walk pred back to a primary input
    let mut worst = 0usize;
    let mut wt = f64::NEG_INFINITY;
    for i in 0..net.n_nets {
        if arrival[i] > wt {
            wt = arrival[i];
            worst = i;
        }
    }
    let mut path = vec![worst];
    let mut cur = worst;
    while pred[cur] != usize::MAX {
        cur = pred[cur];
        path.push(cur);
    }
    path.reverse();
    (arrival, slew, path)
}

/// Certificate: accumulate the per-cell relative delay error along a path of
/// (cell kind, input pin, output, rising) steps. Returns the conservative
/// bound Σ ε_i d_i (s) and the RSS bound √Σ (ε_i d_i)².
pub fn path_error_bound(lib: &Library, steps: &[(CellKind, usize, usize, bool)], slews: &[f64], loads: &[f64]) -> (f64, f64) {
    let mut sum = 0.0;
    let mut rss = 0.0;
    for (k, (kind, pin, o, rise)) in steps.iter().enumerate() {
        let ct = lib.cell(*kind);
        let tab = if *rise { &ct.outs[*o][*pin].rise } else { &ct.outs[*o][*pin].fall };
        let eps = tab.cert_error;
        let (d, _) = tab.interp(slews[k], loads[k]);
        sum += eps * d;
        rss += (eps * d) * (eps * d);
    }
    (sum, rss.sqrt())
}

/// Build the whole standard-cell library: characterise every cell on the
/// given grid and attach held-out certificates.
pub fn build_library(model: &crate::gates::DeviceModel, slews: &[f64], loads: &[f64], held: &[(f64, f64)]) -> Library {
    use rayon::prelude::*;
    let cells: Vec<CellTiming> = CellKind::ALL
        .par_iter()
        .map(|&k| {
            let mut ct = crate::gates::characterise(model, k, slews, loads);
            if !held.is_empty() {
                crate::gates::set_certificates(model, &mut ct, held);
            }
            ct
        })
        .collect();
    Library { cells, vdd: model.vdd }
}

/// Convenience: the timing table of a pin/output direction.
pub fn table<'a>(ct: &'a CellTiming, o: usize, pin: usize, rise: bool) -> &'a TimingTable {
    if rise { &ct.outs[o][pin].rise } else { &ct.outs[o][pin].fall }
}

/// A chain of `n` inverters. `inputs = [0]`, the output is net `n`.
pub fn inverter_chain(n: usize) -> GateNet {
    let mut insts = Vec::new();
    for i in 0..n {
        insts.push(Inst { kind: CellKind::Inv, ins: vec![i], outs: vec![i + 1] });
    }
    GateNet { n_nets: n + 1, inputs: vec![0], insts }
}

/// A `bits`-bit ripple-carry adder of full-adder cells. Returns the netlist,
/// the primary input nets `[a…, b…, cin]`, the sum nets (LSB first) and the
/// final carry net.
pub fn ripple_carry_adder(bits: usize) -> (GateNet, Vec<usize>, Vec<usize>, usize) {
    let a: Vec<usize> = (0..bits).collect();
    let b: Vec<usize> = (bits..2 * bits).collect();
    let cin = 2 * bits;
    let mut n_nets = 2 * bits + 1;
    let mut carry = cin;
    let mut insts = Vec::new();
    let mut sums = Vec::new();
    for i in 0..bits {
        let sum = n_nets;
        n_nets += 1;
        let cout = n_nets;
        n_nets += 1;
        insts.push(Inst { kind: CellKind::FullAdder, ins: vec![a[i], b[i], carry], outs: vec![sum, cout] });
        sums.push(sum);
        carry = cout;
    }
    let mut inputs = a.clone();
    inputs.extend_from_slice(&b);
    inputs.push(cin);
    (GateNet { n_nets, inputs: inputs.clone(), insts }, inputs, sums, carry)
}

/// Evaluate a ripple-carry adder's expected outputs for an input vector
/// (a bits then b bits then cin), for checking the gate simulation.
pub fn expected_sum(bits: usize, in_vec: &[bool]) -> (Vec<bool>, bool) {
    let a: usize = (0..bits).filter(|&i| in_vec[i]).map(|i| 1 << i).sum();
    let b: usize = (0..bits).filter(|&i| in_vec[bits + i]).map(|i| 1 << i).sum();
    let cin = in_vec[2 * bits] as usize;
    let s = a + b + cin;
    let mut out = Vec::with_capacity(bits);
    for i in 0..bits {
        out.push((s >> i) & 1 == 1);
    }
    (out, (s >> bits) & 1 == 1)
}

/// Critical-path arrival and the accumulated interpolation error bound along
/// it. For each gate and output the later of the rise and fall tables chooses
/// the path and contributes its own certificate, so the bound is conservative
/// over the transition direction.
/// Returns (worst arrival, accumulated bound).
pub fn critical_path_bound(lib: &Library, net: &GateNet, net_loads: &[f64], in_slew: f64) -> (f64, f64) {
    let mut arrival = vec![0.0f64; net.n_nets];
    let mut bound = vec![0.0f64; net.n_nets];
    let mut slew = vec![in_slew; net.n_nets];
    for &n in &net.inputs {
        slew[n] = in_slew;
    }
    for inst in &net.insts {
        let ct = lib.cell(inst.kind);
        for (o, &out) in inst.outs.iter().enumerate() {
            let mut best = f64::NEG_INFINITY;
            let mut best_bound = 0.0;
            let mut best_slew = in_slew;
            for (p, &inn) in inst.ins.iter().enumerate() {
                let pin = &ct.outs[o][p];
                let (dr, sr) = pin.rise.interp(slew[inn], net_loads[out]);
                let (df, sf) = pin.fall.interp(slew[inn], net_loads[out]);
                let (d, sl, eps) = if dr >= df {
                    (dr, sr, pin.rise.cert_error)
                } else {
                    (df, sf, pin.fall.cert_error)
                };
                let a = arrival[inn] + d;
                if a > best {
                    best = a;
                    best_slew = sl;
                    best_bound = bound[inn] + eps * d;
                }
            }
            arrival[out] = best;
            bound[out] = best_bound;
            slew[out] = best_slew;
        }
    }
    let wt = arrival.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let wb = bound.iter().cloned().fold(0.0f64, f64::max);
    (wt, wb)
}
