//! Gate-level standard cells and certified timing: the electronics ladder from
//! transistor circuits to logic, towards a CPU (`analysis/emergence_ladder.md`,
//! `analysis/gates.md`).
//!
//! Run: `cargo run -p semiconductor --example adder`.
//!
//! The transistor rung is the EKV compact MOSFET law (the all-region law used
//! as the held-out reference for the PCHIP table model in `compact.rs`),
//! calibrated to the scaled 150 nm technology of
//! `analysis/semiconductor.md`: I_on,n = 330 µA/µm, I_on,p = 230 µA/µm,
//! V_dd = 1.2 V, gate capacitance 1.031 fF/µm. Standard cells (INV, NAND2,
//! NOR2, XOR2, full adder) are transistor netlists; characterisation fits an
//! NLDM timing and energy table; a gate-level event simulator and static
//! timing analysis consume the tables; and the fold is certified against the
//! transistor-level simulation of an inverter chain and a ripple-carry adder.

use semiconductor::gates::{
    simulate_flat, wave_crossing, wave_supply_energy, CellKind, DeviceModel, Ekv, Stimulus, VDD,
};
use semiconductor::timing::{
    build_library, critical_path_bound, inverter_chain, ripple_carry_adder, simulate, sta, Library,
};

const SLEWS: [f64; 4] = [15e-12, 40e-12, 120e-12, 350e-12];
const LOADS: [f64; 3] = [1e-15, 4e-15, 12e-15];
const HELD: [(f64, f64); 4] = [(55e-12, 2e-15), (55e-12, 8e-15), (180e-12, 2e-15), (180e-12, 8e-15)];

fn time_grid(t50: f64, tr: f64, t_end: f64, dt_post: f64) -> Vec<f64> {
    let mut t = vec![0.0];
    let push = |t: &mut Vec<f64>, target: f64, dt: f64| {
        while *t.last().unwrap() < target - 1e-18 {
            let nt = (*t.last().unwrap() + dt).min(target);
            t.push(nt);
        }
    };
    push(&mut t, t50 - 0.8 * tr, 20e-12);
    push(&mut t, t50 + 0.8 * tr, 2e-12);
    push(&mut t, t_end, dt_post);
    t
}

/// Two logic vectors on the primary inputs of a ripple-carry adder.
fn adder_vectors(bits: usize, a: usize, b: usize, cin: bool) -> Vec<bool> {
    let mut v = vec![false; 2 * bits + 1];
    for i in 0..bits {
        v[i] = (a >> i) & 1 == 1;
        v[bits + i] = (b >> i) & 1 == 1;
    }
    v[2 * bits] = cin;
    v
}

fn print_cell_table(lib: &Library) {
    println!("1. Standard-cell library (characterised at input slew 100 ps, load 4 fF)");
    println!(
        "   {:6} {:>6} {:>9} {:>9} {:>10} {:>10} {:>10} {:>8}",
        "cell", "#tr", "c_in(fF)", "c_out(fF)", "d_rise(ps)", "d_fall(ps)", "E_sw(fJ)", "P_leak(nW)"
    );
    for k in CellKind::ALL {
        let ct = lib.cell(k);
        let tr = k.template();
        let s = 100e-12;
        let l = 4e-15;
        let (dr, _) = ct.outs[0][0].rise.interp(s, l);
        let (df, _) = ct.outs[0][0].fall.interp(s, l);
        let e = ct.outs[0][0].e_rise.interp(s, l).0;
        println!(
            "   {:6} {:>6} {:>9.3} {:>9.3} {:>10.2} {:>10.2} {:>10.3} {:>8.3}",
            k.name(),
            tr.devs.len(),
            ct.c_in[0] * 1e15,
            ct.c_out[0] * 1e15,
            dr * 1e12,
            df * 1e12,
            e * 1e15,
            ct.i_leak * VDD * 1e9
        );
    }
    println!("   (FA shown for input pin 0; all inputs have the same c_in.)");
    println!("   Per-cell worst relative delay error certificate from held-out points:");
    for k in CellKind::ALL {
        let ct = lib.cell(k);
        let mut e = 0.0f64;
        for o in 0..k.n_outputs() {
            for p in 0..k.n_inputs() {
                e = e.max(ct.outs[o][p].rise.cert_error).max(ct.outs[o][p].fall.cert_error);
            }
        }
        println!("     {:6}: {:.2}%", k.name(), 100.0 * e);
    }
}

fn chain_fold(lib: &Library, m: &DeviceModel, stages: usize) {
    println!("\n2. Inverter chain ({stages} stages, terminating inverter load), gate-level vs transistor-level");
    // `stages + 1` inverters: the last one loads net `stages`, so the measured
    // path is entirely inside the characterised slew/load box.
    let gnet = inverter_chain(stages + 1);
    let slew = 100e-12;
    let tr = slew / 0.8;
    let t50 = 500e-12;
    for (from, to) in [(false, true), (true, false)] {
        let stim = Stimulus::new(vec![(0.0, vec![from]), (t50, vec![to])], slew);
        let res = simulate(lib, &gnet, &[(0.0, vec![from]), (t50, vec![to])], slew, t50 + 4e-9);
        let fs = simulate_flat(m, &gnet, &stim, &time_grid(t50, tr, t50 + 4e-9, 5e-12));
        let dg = res.last(stages) - t50;
        let dt = wave_crossing(&fs.wave, fs.map[stages], 0.5 * VDD, to, t50 - tr) - t50;
        let eg = res.energy;
        let et = wave_supply_energy(m, &fs.net, &fs.wave);
        let (_ta, bound) = critical_path_bound(lib, &gnet, &lib.loads(&gnet, 0.2e-15), slew);
        let cover = if (dg - dt).abs() <= bound { "within" } else { "OUTSIDE" };
        println!(
            "   input {from}->{to}: delay gate {:6.1} ps, transistor {:6.1} ps ({:+.1}%); interpolation bound {:4.1} ps ({cover})",
            dg * 1e12,
            dt * 1e12,
            100.0 * (dg / dt - 1.0),
            bound * 1e12
        );
        println!(
            "                      energy gate {:.3e} J, transistor {:.3e} J ({:+.1}%)",
            eg,
            et,
            100.0 * (eg / et - 1.0)
        );
    }
}

/// Run one or more input transitions of a ripple-carry adder through both
/// levels of the ladder and print delay/energy errors.
fn adder_fold(lib: &Library, m: &DeviceModel, bits: usize, transitions: &[(usize, usize, bool, usize, usize, bool)], dt_post: f64) {
    println!("\n{bits}-bit ripple-carry adder, gate-level vs transistor-level");
    let (gnet, _in, sums, cout) = ripple_carry_adder(bits);
    let slew = 100e-12;
    let tr = slew / 0.8;
    let t50 = 500e-12;
    let load = lib.loads(&gnet, 0.2e-15);
    let (sta_arr, _slew, _path) = sta(lib, &gnet, &load, slew);
    let (sta_worst, sta_bound) = critical_path_bound(lib, &gnet, &load, slew);
    let sta_simple = sta_arr.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!(
        "   STA worst arrival {:.1} ps (rise edge), {:.1} ps (both-edge, conservative); accumulated interpolation bound {:.1} ps",
        sta_simple * 1e12,
        sta_worst * 1e12,
        sta_bound * 1e12
    );
    for (i, &(a0, b0, c0, a1, b1, c1)) in transitions.iter().enumerate() {
        let v0 = adder_vectors(bits, a0, b0, c0);
        let v1 = adder_vectors(bits, a1, b1, c1);
        let stim = Stimulus::new(vec![(0.0, v0.clone()), (t50, v1.clone())], slew);
        let res = simulate(lib, &gnet, &[(0.0, v0.clone()), (t50, v1.clone())], slew, t50 + 4e-9);
        let t0 = std::time::Instant::now();
        let fs = simulate_flat(m, &gnet, &stim, &time_grid(t50, tr, t50 + 3e-9, dt_post));
        let tsim = t0.elapsed().as_secs_f64();
        let l1 = gnet.eval(&v1);
        let dir_cout = l1[cout];
        let dir_sum = l1[sums[bits - 1]];
        let dc = res.last(cout) - t50;
        let tc = wave_crossing(&fs.wave, fs.map[cout], 0.5 * VDD, dir_cout, t50 - tr) - t50;
        let ds = res.last(sums[bits - 1]) - t50;
        let ts = wave_crossing(&fs.wave, fs.map[sums[bits - 1]], 0.5 * VDD, dir_sum, t50 - tr) - t50;
        let eg = res.energy;
        let et = wave_supply_energy(m, &fs.net, &fs.wave);
        let fmt = |dg: f64, dt: f64| -> String {
            if dg.is_nan() || dt.is_nan() {
                "      (no transition)".to_string()
            } else {
                format!("gate {:6.1} ps, transistor {:6.1} ps ({:+.1}%)", dg * 1e12, dt * 1e12, 100.0 * (dg / dt - 1.0))
            }
        };
        println!(
            "   [{i}] a {a0:0w$b} b {b0:0w$b} -> a {a1:0w$b} b {b1:0w$b}:",
            w = bits
        );
        println!("       carry  {}{}", fmt(dc, tc), if dc.is_nan() { "" } else if (dc - tc).abs() <= sta_bound { " [within bound]" } else { " [OUTSIDE bound]" });
        println!("       sumMSB {}", fmt(ds, ts));
        println!(
            "       energy gate {:.3e} J, transistor {:.3e} J ({:+.1}%); transistor transient {:.1} s",
            eg,
            et,
            100.0 * (eg / et - 1.0),
            tsim
        );
    }
}

fn main() {
    let n = Ekv::n();
    let p = Ekv::p();
    let m = DeviceModel::new(&n, &p);
    println!(
        "EKV compact device law, scaled 150 nm: I_on,n {:.0} µA, I_on,p {:.0} µA (W = 1 µm), V_dd {:.1} V, W_p/W_n {:.2}",
        n.mag(VDD, VDD) * 1e6,
        p.mag(VDD, VDD) * 1e6,
        VDD,
        m.p_scale
    );
    let t0 = std::time::Instant::now();
    let lib = build_library(&m, &SLEWS, &LOADS, &HELD);
    println!("Characterised the standard-cell library in {:.1} s\n", t0.elapsed().as_secs_f64());

    print_cell_table(&lib);

    chain_fold(&lib, &m, 8);

    // 4-bit, several transitions including the worst-case carry ripple
    adder_fold(
        &lib,
        &m,
        4,
        &[
            // (a0,b0,cin0, a1,b1,cin1)
            (0b0000, 0b0000, false, 0b1111, 0b0000, false), // all sums rise
            (0b1111, 0b0000, false, 0b1111, 0b0001, false), // worst-case carry: 15+1=16
            (0b1010, 0b0101, false, 0b1010, 0b0110, false),
            (0b0000, 0b1111, false, 0b0001, 0b1111, false), // carry into bit 3
        ],
        5e-12,
    );

    // 8-bit, the worst-case carry ripple: 255 + 1 = 256
    adder_fold(
        &lib,
        &m,
        8,
        &[(0xff, 0x00, false, 0xff, 0x01, false)],
        6e-12,
    );

    println!("\nFold map");
    println!("   coarse model : NLDM delay/slew/energy tables per cell (INV, NAND2, NOR2, XOR2, FA)");
    println!("   indicator    : interpolation error on held-out points, accumulated along a path");
    println!("   unfold when  : input slew or load leaves the table box (or glitches/coupling matter)");
    println!("   here         : transistor-level EKV netlist, backward-Euler with a scattered Jacobian");
}
