//! Sequential cells: the positive-edge D flip-flop that makes logic into a
//! machine.
//!
//! The cell is the textbook static master–slave flip-flop of four
//! transmission gates and five inverters (18 transistors; e.g. Weste &
//! Harris, *CMOS VLSI Design*, 4th ed., Fig. 10.19): the master latch is
//! transparent while the clock is low, the slave while it is high, so the
//! value at D just before the rising edge appears at Q after it. It is
//! simulated with the same compact transistor law and backward-Euler engine
//! as the combinational cells.
//!
//! Characterisation follows the industry definitions:
//! - **clock-to-Q** delay and output slew, as NLDM tables over clock slew and
//!   load, for a captured 1 (Q rises) and a captured 0 (Q falls);
//! - **setup time**: the smallest data-to-clock separation (50 % to 50 %) at
//!   which the flop still captures with clock-to-Q no more than 10 % longer
//!   than with data long settled (the push-out criterion);
//! - **hold time**: the smallest clock-to-data-change separation at which
//!   the flop still captures the old value with at most 10 % push-out
//!   (negative when the data may change before the edge);
//! - pin capacitances, energy per clock cycle with data steady, and leakage.

use crate::gates::{dc_solve, net_from_template, supply_current, transient, wave_crossing, DeviceModel, Mos, Net, Template};
use crate::timing::TimingTable;

/// Node numbers of the flip-flop template.
pub mod node {
    pub const GND: usize = 0;
    pub const VDD: usize = 1;
    pub const D: usize = 2;
    pub const CLK: usize = 3;
    pub const CLKB: usize = 4;
    pub const M: usize = 5;
    pub const MB: usize = 6;
    pub const MF: usize = 7;
    pub const S: usize = 8;
    pub const Q: usize = 9;
    pub const SF: usize = 10;
}

fn inv(devs: &mut Vec<Mos>, a: usize, y: usize) {
    devs.push(Mos { p: false, g: a, d: y, s: node::GND, scale: 1.0 });
    devs.push(Mos { p: true, g: a, d: y, s: node::VDD, scale: 1.0 });
}

/// Transmission gate between `a` and `b`, on when `on` is high (`onb` low).
fn tg(devs: &mut Vec<Mos>, a: usize, b: usize, on: usize, onb: usize) {
    devs.push(Mos { p: false, g: on, d: b, s: a, scale: 1.0 });
    devs.push(Mos { p: true, g: onb, d: b, s: a, scale: 1.0 });
}

/// The 18-transistor master–slave D flip-flop. Inputs `[D, CLK]`,
/// output `[Q]`.
pub fn dff_template() -> Template {
    use node::*;
    let mut devs = Vec::new();
    inv(&mut devs, CLK, CLKB);
    tg(&mut devs, D, M, CLKB, CLK); // master input, open while CLK is low
    inv(&mut devs, M, MB);
    inv(&mut devs, MB, MF);
    tg(&mut devs, MF, M, CLK, CLKB); // master keeper, closed while CLK is high
    tg(&mut devs, MB, S, CLK, CLKB); // slave input, open while CLK is high
    inv(&mut devs, S, Q);
    inv(&mut devs, Q, SF);
    tg(&mut devs, SF, S, CLKB, CLK); // slave keeper, closed while CLK is low
    Template { n_nodes: 11, inputs: vec![D, CLK], outputs: vec![Q], devs }
}

/// A DC state of the flop with `clk` low or high, `d` at D, and the latches
/// holding `q` (slave) and `m` (master), found from that guess.
fn state(net: &Net, model: &DeviceModel, d: bool, clk: bool, m: bool, q: bool) -> Vec<f64> {
    use node::*;
    let v = |b: bool| if b { model.vdd } else { 0.0 };
    let mut g = vec![0.0; net.n_nodes];
    g[VDD] = model.vdd;
    g[D] = v(d);
    g[CLK] = v(clk);
    g[CLKB] = v(!clk);
    g[M] = v(m);
    g[MB] = v(!m);
    g[MF] = v(m);
    g[S] = v(!q);
    g[Q] = v(q);
    g[SF] = v(!q);
    let gv = g.clone();
    dc_solve(net, model, &|i| gv[i], &g)
}

fn ramp(t: f64, t50: f64, slew: f64, from: f64, to: f64) -> f64 {
    let full = slew / 0.8;
    let (a, b) = (t50 - 0.5 * full, t50 + 0.5 * full);
    if t <= a {
        from
    } else if t >= b {
        to
    } else {
        from + (to - from) * (t - a) / full
    }
}

fn grid(t_end: f64, fine: &[(f64, f64)]) -> Vec<f64> {
    // 2 ps steps, 0.25 ps inside each (centre, half-width) window.
    let mut t = vec![0.0];
    while *t.last().unwrap() < t_end {
        let x = *t.last().unwrap();
        let dt = if fine.iter().any(|&(c, w)| (x - c).abs() < w) { 0.25e-12 } else { 2e-12 };
        t.push((x + dt).min(t_end));
    }
    t
}

/// One capture: the flop holds `!value` at Q with D at `!value`; D moves
/// to `value` with its 50 % point at `t_data`, the clock rises with its 50 %
/// point at `t_clk`; `after` optionally moves D back to `!value` at that
/// time (a hold test). Returns (clock-to-Q delay, output slew, captured
/// value, supply energy), the delay NaN when Q never crossed.
#[allow(clippy::too_many_arguments)]
pub fn capture(
    model: &DeviceModel,
    value: bool,
    t_data: f64,
    t_clk: f64,
    after: Option<f64>,
    clk_slew: f64,
    data_slew: f64,
    load: f64,
) -> (f64, f64, bool, f64) {
    use node::*;
    let tmpl = dff_template();
    let net = net_from_template(model, &tmpl);
    let vdd = model.vdd;
    let v = |b: bool| if b { vdd } else { 0.0 };
    let v0 = state(&net, model, !value, false, !value, !value);
    let t_end = t_clk.max(t_data).max(after.unwrap_or(0.0)) + 400e-12 + 4.0 * load * vdd / 2e-4;
    let mut windows = vec![(t_clk, 0.6 * clk_slew + 80e-12), (t_data, 0.6 * data_slew + 20e-12)];
    if let Some(a) = after {
        windows.push((a, 0.6 * data_slew + 20e-12));
    }
    let times = grid(t_end, &windows);
    let drive = |i: usize, t: f64| -> f64 {
        match i {
            VDD => vdd,
            D => {
                let x = ramp(t, t_data, data_slew, v(!value), v(value));
                match after {
                    Some(a) if t > a - data_slew => ramp(t, a, data_slew, v(value), v(!value)),
                    _ => x,
                }
            }
            CLK => ramp(t, t_clk, clk_slew, 0.0, vdd),
            _ => 0.0,
        }
    };
    let w = transient(&net, model, &drive, &times, &v0, &[(Q, load)]);
    let q_end = w.v.last().unwrap()[Q];
    let captured = q_end > 0.5 * vdd;
    let rising = value;
    let t50 = wave_crossing(&w, Q, 0.5 * vdd, rising, t_clk - 0.5 * clk_slew);
    let slew = if rising {
        wave_crossing(&w, Q, 0.9 * vdd, true, t_clk - clk_slew) - wave_crossing(&w, Q, 0.1 * vdd, true, t_clk - clk_slew)
    } else {
        wave_crossing(&w, Q, 0.1 * vdd, false, t_clk - clk_slew) - wave_crossing(&w, Q, 0.9 * vdd, false, t_clk - clk_slew)
    };
    let mut q = 0.0;
    for k in 1..w.t.len() {
        let i0 = supply_current(&net, model, &w.v[k - 1]);
        let i1 = supply_current(&net, model, &w.v[k]);
        q += 0.5 * (i0 + i1) * (w.t[k] - w.t[k - 1]);
    }
    (t50 - t_clk, slew, captured, vdd * q)
}

/// Characterised flip-flop.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct FlopTiming {
    /// D and CLK pin capacitances (F).
    pub c_d: f64,
    pub c_clk: f64,
    /// Clock-to-Q for a captured 1 (Q rises) and a captured 0 (Q falls), over
    /// (clock slew, load); the `e_*` tables hold supply energy in `delay`.
    pub clk_q_rise: TimingTable,
    pub clk_q_fall: TimingTable,
    pub e_rise: TimingTable,
    pub e_fall: TimingTable,
    /// Setup and hold for a rising and a falling D (s), at the nominal clock
    /// and data slews.
    pub setup_rise: f64,
    pub setup_fall: f64,
    pub hold_rise: f64,
    pub hold_fall: f64,
    /// Energy per clock cycle with D and Q steady (J).
    pub e_clock: f64,
    /// Worst leakage current from V_dd over the four static states (A).
    pub i_leak: f64,
}

impl FlopTiming {
    pub fn setup(&self) -> f64 {
        self.setup_rise.max(self.setup_fall)
    }
    pub fn hold(&self) -> f64 {
        self.hold_rise.max(self.hold_fall)
    }
}

/// Settled clock-to-Q with data long before the edge.
fn nominal(model: &DeviceModel, value: bool, clk_slew: f64, data_slew: f64, load: f64) -> (f64, f64, f64) {
    let t_clk = 300e-12;
    let (d, s, ok, e) = capture(model, value, 60e-12, t_clk, None, clk_slew, data_slew, load);
    assert!(ok == value && d.is_finite(), "flip-flop failed to capture {value} with settled data");
    (d, s, e)
}

/// Setup time for capturing `value`: bisection on the data-to-clock
/// separation for 10 % clock-to-Q push-out.
pub fn setup_time(model: &DeviceModel, value: bool, clk_slew: f64, data_slew: f64, load: f64) -> f64 {
    let (d0, _, _) = nominal(model, value, clk_slew, data_slew, load);
    let t_clk = 400e-12;
    let good = |sep: f64| -> bool {
        let (d, _, ok, _) = capture(model, value, t_clk - sep, t_clk, None, clk_slew, data_slew, load);
        ok == value && d.is_finite() && d <= 1.1 * d0
    };
    let (mut lo, mut hi) = (-150e-12, 300e-12);
    assert!(good(hi), "flip-flop does not capture even with 300 ps of setup");
    if good(lo) {
        return lo;
    }
    for _ in 0..30 {
        let mid = 0.5 * (lo + hi);
        if good(mid) {
            hi = mid;
        } else {
            lo = mid;
        }
        if hi - lo < 0.05e-12 {
            break;
        }
    }
    hi
}

/// Hold time for capturing `value` while D returns to `!value` after the
/// edge: bisection on the clock-to-data separation for a correct capture
/// with at most 10 % push-out.
pub fn hold_time(model: &DeviceModel, value: bool, clk_slew: f64, data_slew: f64, load: f64) -> f64 {
    let (d0, _, _) = nominal(model, value, clk_slew, data_slew, load);
    let t_clk = 400e-12;
    let good = |sep: f64| -> bool {
        let (d, _, ok, _) = capture(model, value, 60e-12, t_clk, Some(t_clk + sep), clk_slew, data_slew, load);
        ok == value && d.is_finite() && d <= 1.1 * d0
    };
    let (mut lo, mut hi) = (-250e-12, 250e-12);
    assert!(good(hi), "flip-flop does not hold even with 250 ps of hold");
    if good(lo) {
        return lo;
    }
    for _ in 0..30 {
        let mid = 0.5 * (lo + hi);
        if good(mid) {
            hi = mid;
        } else {
            lo = mid;
        }
        if hi - lo < 0.05e-12 {
            break;
        }
    }
    hi
}

/// Characterise the flip-flop over clock slews and loads. Setup and hold
/// are taken at the middle clock slew, data slew 50 ps, load 4 fF.
pub fn characterise(model: &DeviceModel, slews: &[f64], loads: &[f64]) -> FlopTiming {
    use node::*;
    let tmpl = dff_template();
    let net = net_from_template(model, &tmpl);
    let data_slew = 50e-12;
    let table = |value: bool| -> (TimingTable, TimingTable) {
        let mut d = vec![vec![0.0; loads.len()]; slews.len()];
        let mut s = d.clone();
        let mut e = d.clone();
        for (i, &cs) in slews.iter().enumerate() {
            for (j, &cl) in loads.iter().enumerate() {
                let (dd, ss, ee) = nominal(model, value, cs, data_slew, cl);
                d[i][j] = dd;
                s[i][j] = ss;
                e[i][j] = ee;
            }
        }
        (TimingTable::new(slews.to_vec(), loads.to_vec(), d, s), TimingTable::new(slews.to_vec(), loads.to_vec(), e, vec![vec![0.0; loads.len()]; slews.len()]))
    };
    let (clk_q_rise, e_rise) = table(true);
    let (clk_q_fall, e_fall) = table(false);
    let cs = slews[slews.len() / 2];
    let cl = 4e-15;
    // Pin capacitances: CLK drives the clock inverter and four transmission
    // gate devices; D sees the input transmission gate's diffusion and,
    // while it is open, the master inverter's gate (the worst case).
    let c_clk = net.devs.iter().filter(|m| m.g == CLK).map(|m| model.gate_cap(m)).sum();
    let c_d = net.devs.iter().filter(|m| m.s == D || m.d == D).map(|m| model.source_cap(m)).sum::<f64>()
        + net.devs.iter().filter(|m| m.g == M).map(|m| model.gate_cap(m)).sum::<f64>();
    // A clock cycle with D and Q steady: rise then fall.
    let e_clock = {
        let v0 = state(&net, model, false, false, false, false);
        let times = grid(1.2e-9, &[(300e-12, 100e-12), (800e-12, 100e-12)]);
        let vdd = model.vdd;
        let drive = |i: usize, t: f64| -> f64 {
            match i {
                VDD => vdd,
                CLK => if t < 550e-12 { ramp(t, 300e-12, cs, 0.0, vdd) } else { ramp(t, 800e-12, cs, vdd, 0.0) },
                _ => 0.0,
            }
        };
        let w = transient(&net, model, &drive, &times, &v0, &[(Q, cl)]);
        let mut q = 0.0;
        for k in 1..w.t.len() {
            q += 0.5 * (supply_current(&net, model, &w.v[k - 1]) + supply_current(&net, model, &w.v[k])) * (w.t[k] - w.t[k - 1]);
        }
        vdd * q
    };
    let mut i_leak = 0.0f64;
    for clk in [false, true] {
        for b in [false, true] {
            let v = state(&net, model, b, clk, b, b);
            i_leak = i_leak.max(supply_current(&net, model, &v));
        }
    }
    FlopTiming {
        c_d,
        c_clk,
        clk_q_rise,
        clk_q_fall,
        e_rise,
        e_fall,
        setup_rise: setup_time(model, true, cs, data_slew, cl),
        setup_fall: setup_time(model, false, cs, data_slew, cl),
        hold_rise: hold_time(model, true, cs, data_slew, cl),
        hold_fall: hold_time(model, false, cs, data_slew, cl),
        e_clock,
        i_leak,
    }
}
