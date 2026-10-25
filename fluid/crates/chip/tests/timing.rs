//! Clocked static timing against the same pipeline simulated transistor by
//! transistor.

use chip::flat::flatten;
use chip::library::{self, Corner};
use chip::rtl::parse;
use chip::sta::{analyse, Constraints};
use chip::synth::synthesize;
use semiconductor::gates::{settle, transient, DeviceModel, Ekv};
use semiconductor::seq::node;

const PIPE: &str = "module pipe(input clk, input [1:0] a, input [1:0] b, output reg [2:0] s);
  reg [1:0] ra, rb;
  always @(posedge clk) begin ra <= a; rb <= b; s <= ra + rb; end
endmodule";

/// Transistor-level: capture ra + rb at an edge `period` after the launch
/// edge. Returns the value in `s` afterwards.
fn captured(period: f64, a: u8, b: u8) -> u8 {
    let m = parse(PIPE).unwrap();
    let n = synthesize(&m).unwrap();
    let (en, ep) = (Ekv::n(), Ekv::p());
    let model = DeviceModel::new(&en, &ep);
    let vdd = model.vdd;
    let f = flatten(&model, &n);
    let slew = 80e-12;
    let full = slew / 0.8;
    let ramp = |t: f64, t50: f64, up: bool| -> f64 {
        let x = ((t - (t50 - 0.5 * full)) / full).clamp(0.0, 1.0);
        vdd * if up { x } else { 1.0 - x }
    };
    // Edges: 0.3 ns (load 0), 1.3 ns (launch a, b), 1.3 ns + period (capture).
    let (e1, e2) = (0.3e-9, 1.3e-9);
    let e3 = e2 + period;
    let clk = move |t: f64| -> f64 {
        if t < e1 + 0.25e-9 {
            ramp(t, e1, true)
        } else if t < e2 - 0.25e-9 {
            ramp(t, e1 + 0.5e-9, false)
        } else if t < e2 + 0.5 * period {
            ramp(t, e2, true)
        } else if t < 0.5 * (e2 + 0.5 * period + e3) {
            ramp(t, e2 + 0.5 * period, false)
        } else {
            ramp(t, e3, true)
        }
    };
    let input_of = |net: usize| -> Option<(u8, usize)> {
        for (name, bits) in &n.inputs {
            if let Some(k) = bits.iter().position(|&b| b == net) {
                return Some((if name == "a" { a } else { b }, k));
            }
        }
        None
    };
    let drive = |node_: usize, t: f64| -> f64 {
        if node_ == 1 {
            return vdd;
        }
        if node_ == f.clock {
            return clk(t);
        }
        if node_ == f.map[n.one] {
            return vdd;
        }
        for net in 0..n.comb.n_nets {
            if f.map[net] == node_ {
                if let Some((v, k)) = input_of(net) {
                    return if (v >> k) & 1 == 1 { ramp(t, 0.8e-9, true) } else { 0.0 };
                }
            }
        }
        0.0
    };
    // Initial state: everything zero, clock low.
    let mut zeros = vec![false; n.comb.inputs.len()];
    zeros[n.comb.inputs.iter().position(|&x| x == n.one).unwrap()] = true;
    let logic = n.comb.eval(&zeros);
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
    let v0 = settle(&f.net, &model, &|i| drive(i, 0.0), &guess);
    let t_end = e3 + 0.6e-9;
    let times: Vec<f64> = (0..=(t_end / 1e-12) as usize).map(|k| k as f64 * 1e-12).collect();
    let w = transient(&f.net, &model, &drive, &times, &v0, &[]);
    let last = w.v.last().unwrap();
    let s = &n.outputs.iter().find(|(name, _)| name == "s").unwrap().1;
    s.iter().enumerate().fold(0u8, |acc, (k, &net)| acc | (((last[f.map[net]] > 0.5 * vdd) as u8) << k))
}

#[test]
fn static_timing_is_safe_and_tight_against_transistors() {
    let lib = library::at(Corner::NOMINAL).unwrap();
    let m = parse(PIPE).unwrap();
    let n = synthesize(&m).unwrap();
    // No wires and no output load: the transistor-level run has neither.
    // Inputs come from upstream flip-flops, a clock-to-Q after the edge
    // (with zero input delay STA rightly reports a 2.8 ps hold violation).
    let c = Constraints { wire_cap: 0.0, output_load: 0.0, input_delay: 100e-12, ..Default::default() };
    let t = analyse(&n, &lib, &c);
    let path: Vec<String> = t.critical.iter().map(|s| format!("{} {} {:.0}", s.cell, s.edge, s.arrival * 1e12)).collect();
    println!("STA min period {:.1} ± {:.1} ps, hold slack {:.1} ps, endpoint {}; path {:?}", t.min_period * 1e12, t.min_period_error * 1e12, t.hold_slack * 1e12, t.endpoint, path);
    assert!(t.hold_slack > 0.0);
    // The worst carry ripple: 00 + 00, then 11 + 01 = 100.
    assert_eq!(captured(1e-9, 3, 1), 4, "captures at a relaxed clock");
    let (mut lo, mut hi) = (0.2 * t.min_period, 1e-9);
    while hi - lo > 2e-12 {
        let mid = 0.5 * (lo + hi);
        if captured(mid, 3, 1) == 4 { hi = mid } else { lo = mid }
    }
    println!("transistor-level shortest period for the carry ripple: {:.1} ps", hi * 1e12);
    // Safe: no pattern fails at the STA period plus its certificate.
    assert!(hi <= t.min_period + t.min_period_error + 2e-12, "STA {:.1} ps is optimistic against {:.1} ps", t.min_period * 1e12, hi * 1e12);
    // And not wildly pessimistic.
    assert!(hi >= 0.7 * t.min_period, "STA {:.1} ps vs transistor {:.1} ps", t.min_period * 1e12, hi * 1e12);
}
