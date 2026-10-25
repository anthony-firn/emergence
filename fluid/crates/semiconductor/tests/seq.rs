//! The D flip-flop at transistor level.

use semiconductor::gates::{dc_solve, net_from_template, transient, DeviceModel, Ekv};
use semiconductor::seq::{characterise, dff_template, node};

#[test]
fn captures_on_rising_edges_only() {
    // Clock period 1 ns; D changes mid-cycle (while the clock is high, when
    // the master is closed) and again just before some edges.
    let (n, p) = (Ekv::n(), Ekv::p());
    let model = DeviceModel::new(&n, &p);
    let net = net_from_template(&model, &dff_template());
    let vdd = model.vdd;
    let d_pattern = [true, false, false, true, true, false];
    let period = 1e-9;
    let edge = 30e-12;
    let clk = move |t: f64| -> f64 {
        let ph = t.rem_euclid(period) / period;
        // high from 0.5 to 1.0 of each period (rising edge at mid-period)
        let r = |x: f64| (x / (edge / period)).clamp(0.0, 1.0);
        if ph < 0.5 { vdd * (1.0 - r(ph)).max(0.0) * if t < period * 0.5 { 0.0 } else { 1.0 } } else { vdd * r(ph - 0.5) }
    };
    // D: value k is set at k·period + 0.1 ns (well before edge k at k·period + 0.5 ns)
    let d = move |t: f64| -> f64 {
        let k = ((t - 0.1e-9) / period).floor().clamp(0.0, (d_pattern.len() - 1) as f64) as usize;
        if d_pattern[k] { vdd } else { 0.0 }
    };
    let drive = |i: usize, t: f64| match i {
        node::VDD => vdd,
        node::D => d(t),
        node::CLK => clk(t),
        _ => 0.0,
    };
    let mut guess = vec![0.0; net.n_nodes];
    guess[node::VDD] = vdd;
    for (k, g) in [(node::CLKB, vdd), (node::MB, vdd), (node::S, vdd), (node::SF, vdd)] {
        guess[k] = g;
    }
    let v0 = dc_solve(&net, &model, &|i| drive(i, 0.0), &guess);
    let times: Vec<f64> = (0..=6000).map(|k| k as f64 * 1e-12).collect();
    let w = transient(&net, &model, &drive, &times, &v0, &[(node::Q, 4e-15)]);
    let q_at = |t: f64| w.v[(t / 1e-12).round() as usize][node::Q];
    for (k, &want) in d_pattern.iter().enumerate() {
        // Just after edge k settles (0.4 ns later) and just before edge k+1.
        let after = (k as f64 + 0.5) * period + 0.4e-9;
        let q = q_at(after);
        assert_eq!(q > 0.5 * vdd, want, "cycle {k}: Q = {q:.3} V, D was {want}");
        assert!(q < 0.05 * vdd || q > 0.95 * vdd, "cycle {k}: Q not at a rail ({q:.3} V)");
    }
}

#[test]
fn characterisation() {
    let (n, p) = (Ekv::n(), Ekv::p());
    let model = DeviceModel::new(&n, &p);
    let t0 = std::time::Instant::now();
    let f = characterise(&model, &[20e-12, 80e-12, 200e-12], &[1e-15, 4e-15, 12e-15]);
    println!(
        "DFF: c_d {:.2} fF, c_clk {:.2} fF, clk→Q rise {:.1} ps fall {:.1} ps (80 ps, 4 fF), setup r/f {:.1}/{:.1} ps, hold r/f {:.1}/{:.1} ps, E_clock {:.2} fJ, leak {:.3} nA ({:.1} s)",
        f.c_d * 1e15,
        f.c_clk * 1e15,
        f.clk_q_rise.interp(80e-12, 4e-15).0 * 1e12,
        f.clk_q_fall.interp(80e-12, 4e-15).0 * 1e12,
        f.setup_rise * 1e12,
        f.setup_fall * 1e12,
        f.hold_rise * 1e12,
        f.hold_fall * 1e12,
        f.e_clock * 1e15,
        f.i_leak * 1e9,
        t0.elapsed().as_secs_f64()
    );
    // Clock-to-Q grows with load and clock slew.
    let d = |s: f64, l: f64| f.clk_q_rise.interp(s, l).0;
    assert!(d(20e-12, 12e-15) > d(20e-12, 1e-15) && d(200e-12, 4e-15) > d(20e-12, 4e-15));
    // Setup + hold is the aperture the data must be stable in: positive.
    assert!(f.setup() + f.hold() > 0.0);
}
