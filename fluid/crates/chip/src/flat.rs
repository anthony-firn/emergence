//! A synthesized design as one transistor network: every cell's and every
//! flip-flop's transistors, with the clock as a driven node. This is the
//! reference the gate-level timing is checked against.

use semiconductor::gates::{DeviceModel, Mos, Net, Template};
use semiconductor::seq::{dff_template, node};

use crate::synth::Netlist;

pub struct Flat {
    pub net: Net,
    /// Transistor node of every gate-level net.
    pub map: Vec<usize>,
    pub clock: usize,
    /// Per flip-flop: transistor nodes of its template (indexed by
    /// `semiconductor::seq::node`).
    pub flops: Vec<Vec<usize>>,
}

fn place(net: &mut Net, model: &DeviceModel, t: &Template, pins: &[(usize, usize)]) -> Vec<usize> {
    let mut lm = vec![usize::MAX; t.n_nodes];
    lm[0] = 0;
    lm[1] = 1;
    for &(local, global) in pins {
        lm[local] = global;
    }
    for slot in lm.iter_mut().skip(2) {
        if *slot == usize::MAX {
            *slot = net.n_nodes;
            net.n_nodes += 1;
            net.driven.push(false);
            net.cap.push(0.0);
        }
    }
    for d in &t.devs {
        let dev = Mos { p: d.p, g: lm[d.g], d: lm[d.d], s: lm[d.s], scale: d.scale };
        net.cap[dev.d] += model.drain_cap(&dev);
        net.cap[dev.s] += model.source_cap(&dev);
        if !net.driven[dev.g] {
            net.cap[dev.g] += model.gate_cap(&dev);
        }
        net.devs.push(dev);
    }
    lm
}

/// Flatten. Primary inputs and the tie nets are driven nodes; flip-flop
/// outputs are not (the flip-flops drive them).
pub fn flatten(model: &DeviceModel, n: &Netlist) -> Flat {
    let mut net = Net { n_nodes: 2, driven: vec![true, true], cap: vec![0.0, 0.0], devs: vec![] };
    let mut map = vec![0usize; n.comb.n_nets];
    for m in map.iter_mut() {
        *m = net.n_nodes;
        net.n_nodes += 1;
        net.driven.push(false);
        net.cap.push(0.0);
    }
    let q_nets: Vec<usize> = n.flops.iter().map(|f| f.q).collect();
    for &i in &n.comb.inputs {
        if !q_nets.contains(&i) {
            net.driven[map[i]] = true;
        }
    }
    let clock = net.n_nodes;
    net.n_nodes += 1;
    net.driven.push(true);
    net.cap.push(0.0);
    for inst in &n.comb.insts {
        let t = inst.kind.template();
        let mut pins: Vec<(usize, usize)> = inst.ins.iter().enumerate().map(|(i, &x)| (t.inputs[i], map[x])).collect();
        pins.extend(inst.outs.iter().enumerate().map(|(o, &x)| (t.outputs[o], map[x])));
        place(&mut net, model, &t, &pins);
    }
    let t = dff_template();
    let flops = n
        .flops
        .iter()
        .map(|f| place(&mut net, model, &t, &[(node::D, map[f.d]), (node::CLK, clock), (node::Q, map[f.q])]))
        .collect();
    Flat { net, map, clock, flops }
}

/// A transistor-level run of a design over clock cycles.
pub struct CycleRun {
    /// Energy drawn from the supply after the first edge (J).
    pub energy: f64,
    /// Simulated time after the first edge (s).
    pub duration: f64,
    /// Output bits (in `Netlist::outputs` order) sampled just before each
    /// clock edge after the first.
    pub outputs: Vec<Vec<bool>>,
}

/// Simulate every transistor: the clock rises every `period` (first edge at
/// one period), input vector k is applied `input_delay` after edge k, all
/// edges with a linear 10–90 % `slew`; time step `dt`. Starts from the
/// all-zero register state with inputs low.
pub fn run_cycles(model: &DeviceModel, n: &Netlist, vectors: &[Vec<bool>], period: f64, input_delay: f64, slew: f64, dt: f64) -> CycleRun {
    use semiconductor::gates::{settle, supply_current, transient};
    let f = flatten(model, n);
    let vdd = model.vdd;
    let full = slew / 0.8;
    let lin = |t: f64, t50: f64| ((t - (t50 - 0.5 * full)) / full).clamp(0.0, 1.0);
    let clk = move |t: f64| -> f64 {
        let k = (t / period).floor();
        let ph = t - k * period;
        // rising at k·period, falling at k·period + period/2
        if k < 1.0 {
            return 0.0;
        }
        if ph < 0.25 * period {
            vdd * lin(ph, 0.0)
        } else if ph < 0.75 * period {
            vdd * (1.0 - lin(ph, 0.5 * period))
        } else {
            vdd * lin(ph, period)
        }
    };
    let in_nets: Vec<usize> = n.inputs.iter().flat_map(|(_, b)| b.iter().copied()).collect();
    let node_input: std::collections::HashMap<usize, usize> = in_nets.iter().enumerate().map(|(k, &net)| (f.map[net], k)).collect();
    let one = f.map[n.one];
    let value = |k: usize, t: f64| -> f64 {
        // bit k: last applied vector's value, ramping from the previous one
        let idx = ((t - input_delay) / period).floor() as isize; // edge index whose vector applies
        let cur = if idx >= 1 { vectors.get((idx - 1) as usize).map(|v| v[k]).unwrap_or(vectors.last().unwrap()[k]) } else { false };
        let prev = if idx >= 2 { vectors[(idx - 2) as usize][k] } else { false };
        let x = if idx >= 1 { lin(t, idx as f64 * period + input_delay) } else { 0.0 };
        vdd * ((prev as u8 as f64) * (1.0 - x) + (cur as u8 as f64) * x)
    };
    let drive = |node_: usize, t: f64| -> f64 {
        if node_ == 1 || node_ == one {
            return vdd;
        }
        if node_ == f.clock {
            return clk(t);
        }
        match node_input.get(&node_) {
            Some(&k) => value(k, t),
            None => 0.0,
        }
    };
    // Settled all-zero start.
    let mut iv = vec![false; n.comb.inputs.len()];
    iv[n.comb.inputs.iter().position(|&x| x == n.one).unwrap()] = true;
    let logic = n.comb.eval(&iv);
    let mut guess = vec![0.0; f.net.n_nodes];
    guess[1] = vdd;
    for net in 0..n.comb.n_nets {
        guess[f.map[net]] = if logic[net] { vdd } else { 0.0 };
    }
    for (fi, nodes) in f.flops.iter().enumerate() {
        let d = logic[n.flops[fi].d];
        let v = |b: bool| if b { vdd } else { 0.0 };
        guess[nodes[node::CLKB]] = vdd;
        guess[nodes[node::M]] = v(d);
        guess[nodes[node::MB]] = v(!d);
        guess[nodes[node::MF]] = v(d);
        guess[nodes[node::S]] = vdd;
        guess[nodes[node::SF]] = vdd;
        guess[nodes[node::Q]] = 0.0;
    }
    let v0 = settle(&f.net, model, &|i| drive(i, 0.0), &guess);
    let t_end = (vectors.len() as f64 + 1.0) * period + 0.9 * period;
    let times: Vec<f64> = (0..=(t_end / dt).ceil() as usize).map(|k| k as f64 * dt).collect();
    let w = transient(&f.net, model, &drive, &times, &v0, &[]);
    let mut energy = 0.0;
    let mut outputs = Vec::new();
    let mut next_sample = 2.0 * period - 0.05 * period;
    for k in 1..w.t.len() {
        if w.t[k] >= period {
            let i0 = supply_current(&f.net, model, &w.v[k - 1]);
            let i1 = supply_current(&f.net, model, &w.v[k]);
            energy += vdd * 0.5 * (i0 + i1) * (w.t[k] - w.t[k - 1]);
        }
        if w.t[k] >= next_sample {
            outputs.push(n.outputs.iter().flat_map(|(_, b)| b.iter().map(|&net| w.v[k][f.map[net]] > 0.5 * vdd)).collect());
            next_sample += period;
        }
    }
    CycleRun { energy, duration: t_end - period, outputs }
}
