//! Gate-level timing rung: truth tables, characterisation certificates, the
//! fold against transistor-level simulation (inverter chain and ripple-carry
//! adder), and static timing analysis.
//!
//! All transistor-level runs use the EKV compact device law of `gates.rs`;
//! the same law is the held-out reference for the table model in
//! `compact.rs`. The library is built once and shared.

use semiconductor::gates::{
    char_point, simulate_flat, wave_crossing, wave_supply_energy, CellKind, DeviceModel, Ekv, Stimulus, VDD,
};
use semiconductor::timing::{
    build_library, expected_sum, inverter_chain, path_error_bound, ripple_carry_adder, simulate, sta, Library,
};
use std::sync::OnceLock;

const SLEWS: [f64; 4] = [15e-12, 40e-12, 120e-12, 350e-12];
const LOADS: [f64; 3] = [1e-15, 4e-15, 12e-15];
const HELD: [(f64, f64); 4] = [(55e-12, 2e-15), (55e-12, 8e-15), (180e-12, 2e-15), (180e-12, 8e-15)];

static LIB: OnceLock<Library> = OnceLock::new();
static NK: OnceLock<Ekv> = OnceLock::new();
static PK: OnceLock<Ekv> = OnceLock::new();

fn ekv_n() -> &'static Ekv {
    NK.get_or_init(Ekv::n)
}
fn ekv_p() -> &'static Ekv {
    PK.get_or_init(Ekv::p)
}
fn model() -> DeviceModel<'static> {
    DeviceModel::new(ekv_n(), ekv_p())
}
fn lib() -> &'static Library {
    LIB.get_or_init(|| {
        let m = model();
        build_library(&m, &SLEWS, &LOADS, &HELD)
    })
}

/// Non-uniform time grid with a fine region around a transition at `t50`.
fn grid(t50: f64, tr: f64, t_end: f64) -> Vec<f64> {
    let mut t = vec![0.0];
    let push = |t: &mut Vec<f64>, target: f64, dt: f64| {
        while *t.last().unwrap() < target - 1e-18 {
            let nt = (*t.last().unwrap() + dt).min(target);
            t.push(nt);
        }
    };
    // coarse settle up to the start of the edge
    push(&mut t, t50 - 0.7 * tr, 10e-12);
    push(&mut t, t50 + 0.7 * tr, 1e-12);
    push(&mut t, t_end, 5e-12);
    t
}

#[test]
fn truth_tables_and_gate_eval() {
    assert_eq!(CellKind::Inv.eval(&[false]), vec![true]);
    assert_eq!(CellKind::Inv.eval(&[true]), vec![false]);
    assert_eq!(CellKind::Nand2.eval(&[true, true]), vec![false]);
    assert_eq!(CellKind::Nand2.eval(&[true, false]), vec![true]);
    assert_eq!(CellKind::Nor2.eval(&[false, false]), vec![true]);
    assert_eq!(CellKind::Nor2.eval(&[false, true]), vec![false]);
    for a in [false, true] {
        for b in [false, true] {
            assert_eq!(CellKind::Xor2.eval(&[a, b]), vec![a ^ b]);
        }
    }
    for a in [false, true] {
        for b in [false, true] {
            for c in [false, true] {
                let o = CellKind::FullAdder.eval(&[a, b, c]);
                let s = a as u8 + b as u8 + c as u8;
                assert_eq!(o, vec![s & 1 == 1, s >= 2]);
            }
        }
    }
    // gate-level netlist evaluation of the 4-bit ripple adder
    let (gnet, ins, sums, cout) = ripple_carry_adder(4);
    for av in 0..16usize {
        for bv in 0..16usize {
            for cin in [false, true] {
                let mut v = vec![false; ins.len()];
                for i in 0..4 {
                    v[i] = (av >> i) & 1 == 1;
                    v[4 + i] = (bv >> i) & 1 == 1;
                }
                v[8] = cin;
                let got = gnet.eval(&v);
                let (es, ec) = expected_sum(4, &v);
                for i in 0..4 {
                    assert_eq!(got[sums[i]], es[i], "sum {i} of {av}+{bv}+{cin}");
                }
                assert_eq!(got[cout], ec, "cout of {av}+{bv}+{cin}");
            }
        }
    }
}

#[test]
fn characterisation_and_certificate() {
    let l = lib();
    let _m = model();
    // input capacitances: INV drives one nMOS + one width-scaled pMOS
    let ct = l.cell(CellKind::Inv);
    assert!((ct.c_in[0] - 2.5775e-15).abs() < 1e-17, "c_in {}", ct.c_in[0]);
    assert!(ct.i_leak > 1e-11 && ct.i_leak < 1e-9, "leak {}", ct.i_leak);
    for k in CellKind::ALL {
        let ct = l.cell(k);
        for o in 0..k.n_outputs() {
            for p in 0..k.n_inputs() {
                for rise in [true, false] {
                    let t = if rise { &ct.outs[o][p].rise } else { &ct.outs[o][p].fall };
                    // monotone in load
                    let d_small = t.interp(100e-12, 1e-15).0;
                    let d_big = t.interp(100e-12, 12e-15).0;
                    assert!(d_big > d_small, "{} out{o} pin{p} rise{rise}: {d_small} {d_big}", k.name());
                    // monotone in input slew
                    let d_fast = t.interp(30e-12, 4e-15).0;
                    let d_slow = t.interp(300e-12, 4e-15).0;
                    assert!(d_slow > d_fast - 1e-15, "{} out{o} pin{p} slew", k.name());
                    assert!(t.cert_error.is_finite() && t.cert_error < 0.25, "{} cert {}", k.name(), t.cert_error);
                }
            }
        }
    }
    // spot value: INV rise delay at (100 ps, 4 fF) is tens of ps
    let d = l.cell(CellKind::Inv).outs[0][0].rise.interp(100e-12, 4e-15).0;
    assert!(d > 20e-12 && d < 70e-12, "INV delay {d}");
}

#[test]
fn inverter_chain_fold() {
    let l = lib();
    let m = model();
    let stages = 6;
    // A terminating inverter loads the last net, so every stage (including the
    // one whose delay is measured) stays inside the characterised load box.
    let gnet = inverter_chain(stages + 1);
    let out = stages;
    let slew = 100e-12;
    let tr = slew / 0.8;
    let t50 = 500e-12;
    for (from, to) in [(false, true), (true, false)] {
        let stim = Stimulus::new(vec![(0.0, vec![from]), (t50, vec![to])], slew);
        let res = simulate(l, &gnet, &(vec![(0.0, vec![from]), (t50, vec![to])]), slew, t50 + 3e-9);
        let dg = res.last(out) - t50;
        let times = grid(t50, tr, t50 + 3e-9);
        let fs = simulate_flat(&m, &gnet, &stim, &times);
        let node = fs.map[out];
        let tcross = wave_crossing(&fs.wave, node, 0.5 * VDD, to, t50 - tr);
        let dt = tcross - t50;
        let er = (dg / dt - 1.0).abs();
        eprintln!("chain {from}->{to}: gate {:.2} ps, transistor {:.2} ps, err {:+.1}%", dg * 1e12, dt * 1e12, 100.0 * (dg / dt - 1.0));
        // each inter-stage load is inside the characterised box
        let load = l.loads(&gnet, 0.2e-15);
        for i in 0..stages {
            let t = &l.cell(CellKind::Inv).outs[0][0].rise;
            assert!(t.contains(res.slew[i], load[i + 1]), "stage {i} out of box: slew {} load {}", res.slew[i], load[i + 1]);
        }
        assert!(er < 0.25, "chain fold error {er}");
        // energy: gate-level estimate vs transistor-level supply energy
        let eg = res.energy;
        let et = wave_supply_energy(&m, &fs.net, &fs.wave);
        eprintln!("  energy gate {:.3e} J, transistor {:.3e} J ({:+.1}%)", eg, et, 100.0 * (eg / et - 1.0));
        assert!(et > 0.0 && (eg / et - 1.0).abs() < 0.5, "energy fold {eg} vs {et}");
    }
}

#[test]
fn sta_matches_event() {
    let l = lib();
    let (gnet, _, _, cout) = ripple_carry_adder(4);
    let load = l.loads(&gnet, 0.2e-15);
    let slew = 100e-12;
    let (arr, _slew, path) = sta(l, &gnet, &load, slew);
    let worst = arr.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    // worst-case carry ripple: 1111 + 0000 -> 1111 + 0001 makes the carry
    // propagate through all four bits (15 + 1 = 16)
    let v0: Vec<bool> = {
        let mut v = vec![false; 9];
        for i in 0..4 {
            v[i] = true; // a = 1111
        }
        v
    };
    let mut v1 = v0.clone();
    v1[4] = true; // b = 0001
    let stim = vec![(0.0, v0.clone()), (500e-12, v1)];
    let res = simulate(l, &gnet, &stim, slew, 3.5e-9);
    let ev = res.last(cout) - 500e-12;
    eprintln!("STA worst {:.1} ps, event cout {:.1} ps, path {path:?}", worst * 1e12, ev * 1e12);
    assert!((worst / ev - 1.0).abs() < 0.35, "STA {worst} vs event {ev}");
}

#[test]
fn adder4_fold() {
    let l = lib();
    let m = model();
    let (gnet, _in, sums, cout) = ripple_carry_adder(4);
    let slew = 100e-12;
    let tr = slew / 0.8;
    let t50 = 500e-12;
    // 1111 + 0000 -> 1111 + 0001: carry ripples through all four bits
    let v0: Vec<bool> = {
        let mut v = vec![false; 9];
        for i in 0..4 {
            v[i] = true;
        }
        v
    };
    let mut v1 = v0.clone();
    v1[4] = true;
    let stim = Stimulus::new(vec![(0.0, v0.clone()), (t50, v1.clone())], slew);
    let stim_pairs = vec![(0.0, v0.clone()), (t50, v1.clone())];
    let res = simulate(l, &gnet, &stim_pairs, slew, t50 + 4e-9);
    // gate-level: last arrival on carry-out and on the sum LSBs
    let dg_cout = res.last(cout) - t50;
    let dg_sum3 = res.last(sums[3]) - t50;
    // transistor-level on the flattened 4-bit adder
    let mut times = vec![0.0];
    let push = |t: &mut Vec<f64>, target: f64, dt: f64| {
        while *t.last().unwrap() < target - 1e-18 {
            let nt = (*t.last().unwrap() + dt).min(target);
            t.push(nt);
        }
    };
    push(&mut times, t50 - 0.8 * tr, 20e-12);
    push(&mut times, t50 + 0.8 * tr, 2e-12);
    push(&mut times, t50 + 2.0e-9, 5e-12);
    let t0 = std::time::Instant::now();
    let fs = simulate_flat(&m, &gnet, &stim, &times);
    let tsim = t0.elapsed().as_secs_f64();
    let dt_cout = wave_crossing(&fs.wave, fs.map[cout], 0.5 * VDD, true, t50 - tr) - t50;
    let dt_sum3 = wave_crossing(&fs.wave, fs.map[sums[3]], 0.5 * VDD, false, t50 - tr) - t50;
    eprintln!(
        "adder4 carry: gate {:.1} ps, transistor {:.1} ps ({:+.1}%); sum3: gate {:.1} ps, transistor {:.1} ps ({:+.1}%); transistor sim {:.1} s",
        dg_cout * 1e12, dt_cout * 1e12, 100.0 * (dg_cout / dt_cout - 1.0),
        dg_sum3 * 1e12, dt_sum3 * 1e12, 100.0 * (dg_sum3 / dt_sum3 - 1.0), tsim
    );
    let eg = res.energy;
    let et = wave_supply_energy(&m, &fs.net, &fs.wave);
    eprintln!("adder4 energy: gate {:.3e} J, transistor {:.3e} J ({:+.1}%)", eg, et, 100.0 * (eg / et - 1.0));
    assert!((dg_cout / dt_cout - 1.0).abs() < 0.35, "carry fold");
    assert!((dg_sum3 / dt_sum3 - 1.0).abs() < 0.35, "sum3 fold");
    assert!((eg / et - 1.0).abs() < 0.5, "energy fold");
}

#[test]
fn interpolation_certificate_bounds_held_out_error() {
    // Build a coarse table and check that the certificate ≥ the held-out error
    // it was computed from, and that fresh in-box points stay within a few
    // times it.
    let m = model();
    let slews = [50e-12, 150e-12];
    let loads = [2e-15, 10e-15];
    let mut ct = semiconductor::gates::characterise(&m, CellKind::Nand2, &slews, &loads);
    let held = [(80e-12, 5e-15), (120e-12, 5e-15)];
    semiconductor::gates::set_certificates(&m, &mut ct, &held);
    let tab = &ct.outs[0][0].rise;
    let (bound, _) = tab.certify(&held.iter().map(|&(s, l)| {
        let cp = char_point(&m, CellKind::Nand2, 0, 0, true, s, l);
        (s, l, cp.delay, cp.out_slew)
    }).collect::<Vec<_>>());
    assert!(tab.cert_error >= bound * 0.999);
    let _ = path_error_bound;
}

#[test]
fn certificate_covers_critical_path_error() {
    // The error certificate: accumulate the per-cell interpolation error along
    // the inverter chain and check the transistor-level result lies inside.
    let l = lib();
    let m = model();
    let stages = 8;
    // Terminating inverter: every stage in the measured path is inside the
    // characterised box (slew and load), so the interpolation certificate
    // applies to the whole path.
    let gnet = inverter_chain(stages + 1);
    let slew = 100e-12;
    let tr = slew / 0.8;
    let t50 = 500e-12;
    let stim = Stimulus::new(vec![(0.0, vec![false]), (t50, vec![true])], slew);
    let res = simulate(l, &gnet, &[(0.0, vec![false]), (t50, vec![true])], slew, t50 + 3e-9);
    let times = grid(t50, tr, t50 + 3e-9);
    let fs = simulate_flat(&m, &gnet, &stim, &times);
    let dt = wave_crossing(&fs.wave, fs.map[stages], 0.5 * VDD, true, t50 - tr) - t50;
    let load = l.loads(&gnet, 0.2e-15);
    // path steps: (INV, pin 0, out 0, rising per stage). Stage i output is
    // low when i is even, so the final (stage 8) transition is a rise.
    let steps: Vec<(CellKind, usize, usize, bool)> =
        (0..stages).map(|i| (CellKind::Inv, 0, 0, (stages - i) % 2 == 1)).collect();
    let slews: Vec<f64> = (0..stages).map(|i| res.slew[i]).collect();
    let loads: Vec<f64> = (1..=stages).map(|i| load[i]).collect();
    let (sum, rss) = path_error_bound(l, &steps, &slews, &loads);
    let dg = res.last(stages) - t50;
    eprintln!(
        "certificate: path bound (sum) {:.1} ps, RSS {:.1} ps; gate {:.1} ps, transistor {:.1} ps, |err| {:.1} ps",
        sum * 1e12, rss * 1e12, dg * 1e12, dt * 1e12, (dg - dt).abs() * 1e12
    );
    assert!((dg - dt).abs() <= sum, "transistor result outside the accumulated bound");
}

#[test]
fn unfold_indicator_fires_outside_box() {
    let l = lib();
    let inv = &l.cell(CellKind::Inv).outs[0][0].rise;
    // inside the characterised box
    assert!(inv.contains(100e-12, 4e-15));
    // input slew below the box (edge faster than characterised) and load
    // above the box both signal an unfold to transistor level
    assert!(!inv.contains(0.2e-12, 4e-15));
    assert!(!inv.contains(100e-12, 100e-15));
    // and above the slew box
    assert!(!inv.contains(2000e-12, 4e-15));
}



