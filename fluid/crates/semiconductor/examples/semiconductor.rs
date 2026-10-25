//! Silicon and transistors, from material to circuit.
//!
//!   semiconductor mms         manufactured solutions: convergence order of
//!                             the drift–diffusion discretization (1D, 2D)
//!   semiconductor silicon     band gap, intrinsic density vs measurement,
//!                             mobility
//!   semiconductor diode       pn junction: equilibrium, depletion
//!                             capacitance, Shockley diode law
//!   semiconductor moscap      MOS capacitor vs the exact Poisson–Boltzmann
//!                             solution: gate charge and low-frequency C–V
//!   semiconductor mosfet [L]  2D nMOSFET (gate length L µm) vs Pao–Sah
//!   semiconductor short       short-channel effects in a scaled technology:
//!                             threshold roll-off, DIBL, punch-through vs
//!                             gate length; the natural length
//!   semiconductor mobility    inversion-layer mobility (Lombardi): effective
//!                             mobility vs effective field, universality
//!   semiconductor transient   transient TCAD: MOS capacitor under gate ramps
//!                             (deep depletion ↔ low-frequency C–V), CMOS
//!                             inverter switching: quasi-static table model
//!                             vs mixed-mode transient TCAD
//!   semiconductor switching   CMOS inverter switching: quasi-static table
//!                             model vs mixed-mode transient TCAD
//!   semiconductor thermal     heat: solver vs exact (nonlinear k), energy
//!                             balance, diode J₀(T), leakage vs T, self-
//!                             heating, coolant (Graetz), thermal runaway
//!   semiconductor quantum     inversion-layer confinement: Schrödinger–
//!                             Poisson vs density gradient vs classical
//!   semiconductor quantum2d   density gradient inside the 2D solver: MOS
//!                             capacitor vs 1D references, 150 nm MOSFET
//!   semiconductor tunnel      gate tunnelling: Fowler–Nordheim limit,
//!                             leakage vs oxide thickness
//!   semiconductor breakdown   avalanche breakdown of a p⁺n junction (impact
//!                             ionization) vs the ionization integral
//!   semiconductor cmos        fold TCAD nMOS/pMOS into certified table
//!                             models; inverter transfer curve and noise
//!                             margins vs mixed-mode TCAD; ring oscillator
use semiconductor::analytic::{depletion_capacitance, shockley_j0, Mos};
use semiconductor::circuit::{crossing, inverter_transient_qs, inverter_vout, ring_oscillator, vtc};
use semiconductor::compact::{fold, sample_with_charges, tcad_id, QTable, Table};
use semiconductor::device::{Contact, ContactKind, Device, Statistics};
use semiconductor::mesh::{graded, uniform, Material, Mesh};
use semiconductor::mosfet::MosfetSpec;
use semiconductor::quantum::{Electrons, MosStack};
use semiconductor::thermal::{graetz, runaway, DeviceHeat};
use semiconductor::tunnel::{fowler_nordheim, tsu_esaki};
use semiconductor::silicon::{alpha_n, alpha_p, band_gap, bgn_slotboom, degeneracy, nc, ni_from_bands, ni_measured, nv, Silicon, EPS_OX, EPS_SI, Q};

/// Both node rows at x index i of a 1D strip.
fn line_nodes(m: &Mesh, i: usize) -> Vec<usize> {
    (0..m.ny()).map(|j| m.node(i, j)).collect()
}

fn diode() {
    let si = Silicon::new(300.0);
    let vt = si.vt();
    let (na, nd) = (1e23, 1e22); // 1e17 / 1e16 cm⁻³
    let mun = si.mob_n.mobility(na, si.t);
    let mup = si.mob_p.mobility(nd, si.t);
    let (tn, tp) = (si.tau_n(na), si.tau_p(nd));
    let ln = (mun * vt * tn).sqrt();
    let lp = (mup * vt * tp).sqrt();
    let (lpside, lnside) = (10.0 * ln, 10.0 * lp);
    let xj = lpside;
    let x = graded(&[(0.0, lpside / 200.0), (xj, 0.5e-9), (xj + lnside, lnside / 200.0)]);
    let mesh = Mesh::line(x, |_| Material::Silicon);
    let last = mesh.nx() - 1;
    let contacts = vec![
        Contact { name: "anode".into(), nodes: line_nodes(&mesh, 0), kind: ContactKind::Ohmic, v: 0.0 },
        Contact { name: "cathode".into(), nodes: line_nodes(&mesh, last), kind: ContactKind::Ohmic, v: 0.0 },
    ];
    let mut dev = Device::new(mesh.clone(), si.clone(), move |x, _| if x < xj { (0.0, na) } else { (nd, 0.0) }, contacts, 1.0);
    let its = dev.solve().expect("equilibrium");
    let an = dev.contact("anode");
    let (jn, jp) = dev.current(an);
    let fmax = (0..mesh.nodes()).map(|k| dev.phin[k].abs().max(dev.phip[k].abs())).fold(0.0, f64::max);
    let vbi = vt * ((0.5 * na / si.ni).asinh() + (0.5 * nd / si.ni).asinh());
    println!("abrupt pn junction, N_a 1e17 / N_d 1e16 cm⁻³, {} nodes, long base (10 diffusion lengths: L_n {:.1} µm, L_p {:.1} µm)", mesh.nx(), ln * 1e6, lp * 1e6);
    println!("  equilibrium ({its} Newton steps): V_bi {vbi:.4} V, terminal current {:.1e} A/m², max |quasi-Fermi level| {:.1e} V (exactly flat: detailed balance)", jn + jp, fmax * vt);
    // depletion capacitance at reverse bias: net positive charge on the n side
    let qn = |d: &Device| -> f64 { (0..mesh.nodes()).filter(|&k| mesh.x[mesh.ij(k).0] > xj).map(|k| Q * d.semiconductor_volume(k) * (d.p(k) - d.n(k) + d.nd[k] - d.na[k])).sum() };
    println!("  depletion capacitance (dQ/dV, ±10 mV) vs abrupt-junction theory with the 2V_t tail correction:");
    for v in [0.0, -0.5, -1.0, -2.0, -5.0] {
        let mut c = [0.0; 2];
        for (s, dv) in [(0, -0.01), (1, 0.01)] {
            dev.ramp(an, v + dv, 0.25).expect("bias");
            c[s] = qn(&dev);
        }
        let cap = (c[0] - c[1]) / 0.02;
        let th = depletion_capacitance(EPS_SI, na, nd, vbi, v, vt);
        println!("    V {v:+5.1}: C {:8.3} nF/cm² (theory {:8.3}, {:+.2}%)", cap * 1e5, th * 1e5, 100.0 * (cap / th - 1.0));
    }
    let j0 = shockley_j0(si.ni, vt, na, nd, mun, tn, mup, tp);
    println!("  forward bias vs Shockley J₀(e^(V/V_t) − 1), J₀ = {:.3e} A/cm²:", j0 * 1e-4);
    let mut prev: Option<(f64, f64)> = None;
    for v in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8] {
        dev.ramp(an, v, 0.05).expect("forward");
        let (jn, jp) = dev.current(an);
        let j = jn + jp;
        let sh = j0 * ((v / vt).exp() - 1.0);
        let ideal = prev.map(|(v0, j0_)| (v - v0) / (vt * (j / j0_).ln()));
        println!("    V {v:.2}: J {:.4e} A/cm² (Shockley {:.4e}, ratio {:.4}){}", j * 1e-4, sh * 1e-4, j / sh, ideal.map_or(String::new(), |n| format!(", ideality {n:.3}")));
        prev = Some((v, j));
    }
    dev.ramp(an, -1.0, 0.25).expect("reverse");
    let (jn, jp) = dev.current(an);
    println!("    V −1.00: J {:.4e} A/cm² (−J₀ = {:.4e}; the rest is SRH generation in the depletion region)", (jn + jp) * 1e-4, -j0 * 1e-4);
    // Shockley's assumptions made exact: no recombination inside the
    // depletion region (edges from the depletion approximation at each
    // bias), so all current is minority-carrier diffusion in the neutral
    // regions
    println!("  the same diode without depletion-region recombination (Shockley's assumption):");
    let (mun_n, mup_p) = (si.mob_n.mobility(nd, si.t), si.mob_p.mobility(na, si.t));
    for v in [-1.0, 0.1, 0.3, 0.5, 0.55, 0.6] {
        let w = (2.0 * EPS_SI * (vbi - v) / Q * (na + nd) / (na * nd)).sqrt();
        let (xp, xn) = (w * nd / (na + nd), w * na / (na + nd));
        dev.recomb_scale = Some((0..mesh.nodes()).map(|k| { let x = mesh.x[mesh.ij(k).0]; if x > xj - xp && x < xj + xn { 0.0 } else { 1.0 } }).collect());
        dev.ramp(an, v, 0.05).expect("bias");
        dev.solve().expect("bias");
        let (jn, jp) = dev.current(an);
        let j = jn + jp;
        let sh = j0 * ((v / vt).exp() - 1.0);
        // ohmic drop in the neutral regions (series resistance)
        let rs = (lnside - xn) / (Q * nd * mun_n) + (lpside - xp) / (Q * na * mup_p);
        let sh_r = j0 * (((v - j * rs) / vt).exp() - 1.0);
        println!("    V {v:+.2}: J {:+.5e} A/cm² (Shockley {:+.5e}, ratio {:.5}; at V − J·R_s: ratio {:.5})", j * 1e-4, sh * 1e-4, j / sh, j / sh_r);
    }
}

/// Long-channel nMOSFET at constant mobility vs the Pao–Sah double
/// integral (the exact gradual-channel limit of drift–diffusion).
fn mosfet(l_um: f64) {
    let si = Silicon::new(300.0);
    let vt = si.vt();
    let spec = MosfetSpec::nmos(l_um * 1e-6);
    let mut dev = spec.build(si.clone());
    let (mun, mup) = (0.05, 0.02);
    dev.mu_const = Some((mun, mup));
    dev.update_mobility();
    let t0 = std::time::Instant::now();
    let it = dev.solve().expect("equilibrium");
    println!(
        "nMOSFET: gate {l_um} µm (metallurgical channel {:.4} µm), t_ox 10 nm, N_a 1e17, n⁺ 1e20 cm⁻³, constant μ_n 500 cm²/(V s); mesh {}×{} ({} unknowns); equilibrium in {it} Newton steps, {:.2} s",
        spec.l_eff() * 1e6,
        dev.mesh.nx(),
        dev.mesh.ny(),
        3 * dev.mesh.nodes(),
        t0.elapsed().as_secs_f64()
    );
    let mos = Mos { na: spec.n_sub, ni: si.ni, vt, eps_s: EPS_SI, cox: EPS_OX / spec.tox, vfb: -vt * (0.5 * spec.n_sub / si.ni).asinh() - spec.gate_offset(si.t) };
    let (d, g) = (dev.contact("drain"), dev.contact("gate"));
    for vd in [0.05, 1.0] {
        dev.ramp(g, -0.2, 0.2).expect("gate");
        dev.ramp(d, vd, 0.1).expect("drain");
        println!("  V_d {vd} V: I_d vs Pao–Sah (W = 1 µm)");
        let mut prev: Option<(f64, f64, f64)> = None;
        let t0 = std::time::Instant::now();
        for i in 0..=27 {
            let vg = -0.2 + 0.1 * i as f64;
            dev.ramp(g, vg, 0.1).expect("gate");
            let id = dev.current(d).0 + dev.current(d).1;
            let ps = mos.pao_sah(vg, vd, spec.l_eff(), mun) * spec.width;
            let ss = prev.map(|(v0, i0, p0)| (1e3 * (vg - v0) / (id / i0).log10(), 1e3 * (vg - v0) / (ps / p0).log10()));
            if i % 2 == 0 || vg < 0.45 {
                println!(
                    "    V_g {vg:+.1}: I_d {id:.4e} A (Pao–Sah {ps:.4e}, ratio {:.4}){}",
                    id / ps,
                    ss.filter(|_| vg < 0.4).map_or(String::new(), |(a, b)| format!(", slope {a:.1} mV/dec (Pao–Sah {b:.1})"))
                );
            }
            prev = Some((vg, id, ps));
        }
        println!("    ({:.1} s for the sweep)", t0.elapsed().as_secs_f64());
    }
}

/// Threshold (constant-current, I_d = 1e-7 A·W/L) and leakage of one
/// device at drain voltage vd; also returns I_on at V_g = V_d = vdd.
fn vt_ioff(dev: &mut Device, spec: &MosfetSpec, vd: f64, vdd: f64) -> (f64, f64, f64, f64) {
    let (d, g) = (dev.contact("drain"), dev.contact("gate"));
    dev.ramp(g, 0.0, 0.2).expect("gate");
    dev.ramp(d, vd, 0.1).expect("drain");
    let icrit = 1e-7 * spec.width / spec.l_gate;
    let id = |dev: &Device| { let (a, b) = dev.current(d); a + b };
    let ioff = id(dev);
    let mut prev = (0.0, ioff);
    let mut vt = f64::NAN;
    let mut ss = f64::NAN;
    let mut v = 0.0;
    while v < vdd - 1e-9 {
        v += 0.025;
        dev.ramp(g, v, 0.05).expect("gate");
        let i = id(dev);
        if v <= 0.1 + 1e-9 {
            ss = 1e3 * (v - prev.0) / (i / prev.1).log10();
        }
        if vt.is_nan() && i >= icrit {
            let t = (icrit / prev.1).ln() / (i / prev.1).ln();
            vt = prev.0 + t * (v - prev.0);
            break;
        }
        prev = (v, i);
    }
    dev.ramp(g, vdd, 0.1).expect("gate");
    let ion = if (vd - vdd).abs() < 1e-12 { id(dev) } else { f64::NAN };
    (vt, ioff, ion, ss)
}

/// Short-channel effects in a scaled technology: drift–diffusion with
/// doping-dependent mobility and velocity saturation.
fn short_channel() {
    let si = Silicon::new(300.0);
    let vt_ = si.vt();
    let vdd = 1.2;
    let base = MosfetSpec::scaled(1e-6, false);
    let na = base.n_sub;
    let phi_f = vt_ * (0.5 * na / si.ni).asinh();
    let wd = (2.0 * EPS_SI * 2.0 * phi_f / (Q * na)).sqrt();
    let l_nat = (EPS_SI / EPS_OX * base.tox * wd).sqrt();
    println!("short-channel effects: nMOSFET, t_ox 4 nm, N_a 5e17 cm⁻³, junctions 60 nm, V_dd {vdd} V; Masetti mobility + velocity saturation");
    println!("  natural length (Liu et al. 1993, η = 1): l = √(ε_si t_ox W_d/ε_ox) = {:.1} nm (W_d {:.1} nm)", l_nat * 1e9, wd * 1e9);
    println!("  {:>7} {:>8} {:>9} {:>9} {:>8} {:>10} {:>10} {:>9} {:>6}", "L_gate", "L_met", "V_t lin", "V_t sat", "DIBL", "I_off", "I_on", "SS sat", "time");
    let mut rows = Vec::new();
    for l_nm in [1000.0, 300.0, 200.0, 150.0, 120.0, 100.0, 80.0, 65.0] {
        let spec = MosfetSpec::scaled(l_nm * 1e-9, false);
        let mut dev = spec.build(si.clone());
        dev.velocity_saturation = true;
        let t0 = std::time::Instant::now();
        dev.solve().expect("equilibrium");
        let (vtl, _, _, _) = vt_ioff(&mut dev, &spec, 0.05, vdd);
        let d = dev.contact("drain");
        dev.ramp(dev.contact("gate"), 0.0, 0.2).expect("gate");
        dev.ramp(d, 0.0, 0.2).expect("drain");
        let (vts, ioff, ion, ss) = vt_ioff(&mut dev, &spec, vdd, vdd);
        let dibl = (vtl - vts) / (vdd - 0.05) * 1e3;
        println!(
            "  {:5.0}nm {:6.1}nm {:8.4}V {:8.4}V {:5.0}mV/V {:9.2e}A {:9.2e}A {:6.1}mV/dec  {:.0}s",
            l_nm,
            spec.l_eff() * 1e9,
            vtl,
            vts,
            dibl,
            ioff,
            ion,
            ss,
            t0.elapsed().as_secs_f64()
        );
        rows.push((spec.l_eff(), vtl, vts));
    }
    // decay length of the roll-off: ΔV_t ∝ exp(−L/2l), locally between
    // successive gate lengths (Liu's form has a prefactor and a second
    // exponential, so the local length approaches l only at short L)
    let vt_long = rows[0].1;
    let lens: Vec<String> = rows
        .windows(2)
        .skip(1)
        .map(|w| {
            let (d0, d1) = (vt_long - w[0].1, vt_long - w[1].1);
            let l = 0.5 * (w[0].0 - w[1].0) / (d1 / d0).ln();
            format!("{:.0}–{:.0} nm: {:.1}", w[0].0 * 1e9, w[1].0 * 1e9, l * 1e9)
        })
        .collect();
    println!("  local decay length of the roll-off ΔV_t(lin) ∝ exp(−L_met/2l) between successive devices (nm): {}", lens.join(", "));
    println!("  natural length with η = 1: {:.1} nm", l_nat * 1e9);
}

/// Two TCAD devices coupled by Kirchhoff's current law at the inverter
/// output (mixed-mode): regula falsi (Illinois) on V_out.
fn tcad_inverter(n: &mut Device, p: &mut Device, p_scale: f64, vdd: f64, vin: f64) -> (f64, usize) {
    let (nd, ng) = (n.contact("drain"), n.contact("gate"));
    let (pd, pg) = (p.contact("drain"), p.contact("gate"));
    n.ramp(ng, vin, 0.1).expect("n gate");
    p.ramp(pg, vin - vdd, 0.1).expect("p gate");
    let mut evals = 0;
    let mut f = |vout: f64, n: &mut Device, p: &mut Device| -> f64 {
        evals += 1;
        n.ramp(nd, vout, 0.1).expect("n drain");
        p.ramp(pd, vout - vdd, 0.1).expect("p drain");
        -p_scale * tcad_id(p) - tcad_id(n)
    };
    let (mut a, mut b) = (0.0, vdd);
    let (mut fa, mut fb) = (f(a, n, p), f(b, n, p));
    let mut side = 0;
    let mut c = 0.5 * vdd;
    for _ in 0..40 {
        c = (a * fb - b * fa) / (fb - fa);
        let fc = f(c, n, p);
        if (b - a).abs() < 1e-6 || fc.abs() < 1e-14 {
            break;
        }
        if fc * fb > 0.0 {
            b = c;
            fb = fc;
            if side == -1 {
                fa *= 0.5;
            }
            side = -1;
        } else {
            a = c;
            fa = fc;
            if side == 1 {
                fb *= 0.5;
            }
            side = 1;
        }
        if (b - a).abs() < 2e-6 {
            break;
        }
    }
    (c, evals)
}

fn cmos() {
    let si = Silicon::new(300.0);
    let vt = si.vt();
    let vdd = 1.2;
    let l = 150e-9;
    let (sn, sp) = (MosfetSpec::scaled(l, false), MosfetSpec::scaled(l, true));
    let mut n = sn.build(si.clone());
    let mut p = sp.build(si.clone());
    for d in [&mut n, &mut p] {
        d.velocity_saturation = true;
        d.solve().expect("equilibrium");
    }
    let vgs: Vec<f64> = (0..=12).map(|i| 0.1 * i as f64).collect();
    let vds: Vec<f64> = vec![0.0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2];
    println!("CMOS from TCAD: n- and pMOS, gate 150 nm, scaled technology (t_ox 4 nm, 5e17 cm⁻³), V_dd {vdd} V, W = 1 µm");
    let t0 = std::time::Instant::now();
    let (tn, en, dn) = fold(&sn, &mut n, &vgs, &vds);
    let (mut tp, ep, dp) = fold(&sp, &mut p, &vgs, &vds);
    println!("  folded into table models ({}×{} TCAD points each, {:.0} s); certificate at the {} cell centres never used for fitting:", vgs.len(), vds.len(), t0.elapsed().as_secs_f64(), (vgs.len() - 1) * (vds.len() - 1));
    println!("    nMOS: worst relative error {:.2}% above 1e-3 I_max, worst {:.3} decades below", 100.0 * en, dn);
    println!("    pMOS: worst relative error {:.2}% above 1e-3 I_max, worst {:.3} decades below", 100.0 * ep, dp);
    let ion_n = tn.id(vdd, vdd);
    let ion_p = -tp.id(-vdd, -vdd);
    println!("    I_on: nMOS {:.1} µA/µm, pMOS {:.1} µA/µm (ratio {:.2}: electrons are faster)", ion_n * 1e6, ion_p * 1e6, ion_n / ion_p);
    // size the pMOS for a symmetric switching threshold
    let (mut lo, mut hi) = (0.5, 6.0);
    for _ in 0..50 {
        let mid = 0.5 * (lo + hi);
        tp.scale = mid;
        if inverter_vout(&tn, &tp, vdd, 0.5 * vdd) < 0.5 * vdd {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let ratio = 0.5 * (lo + hi);
    tp.scale = ratio;
    let (_curve, c) = vtc(&tn, &tp, vdd, 480);
    println!("  inverter (W_p/W_n = {ratio:.2} for V_M = V_dd/2), table models:");
    println!("    V_M {:.4} V, peak gain {:.1}, V_IL {:.4} V, V_IH {:.4} V, V_OL {:.4} V, V_OH {:.4} V", c.vm, c.gain, c.vil, c.vih, c.vol, c.voh);
    println!("    noise margins NM_L {:.3} V, NM_H {:.3} V (> 0: the Boolean abstraction of the gate level is valid)", c.nml, c.nmh);
    println!("  check against mixed-mode TCAD (two drift–diffusion devices coupled by KCL):");
    for vin in [0.2, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 1.0] {
        let vm = inverter_vout(&tn, &tp, vdd, vin);
        let (vtcad, ev) = tcad_inverter(&mut n, &mut p, ratio, vdd, vin);
        println!("    V_in {vin:.2}: V_out table {vm:.4} V, TCAD {vtcad:.4} V ({:+.1} mV, {ev} TCAD evaluations)", 1e3 * (vm - vtcad));
    }
    // node capacitance from TCAD gate charges and junction depletion
    let qg = |d: &mut Device, vg: f64| -> f64 {
        let (dd, g) = (d.contact("drain"), d.contact("gate"));
        d.ramp(dd, 0.0, 0.1).expect("drain");
        d.ramp(g, vg, 0.1).expect("gate");
        d.charge(g)
    };
    let cg_n = (qg(&mut n, vdd) - qg(&mut n, 0.0)) / vdd;
    let cg_p = (qg(&mut p, -vdd) - qg(&mut p, 0.0)) / -vdd;
    let vbi = vt * ((0.5 * sn.n_sub / si.ni).asinh() + (0.5 * sn.n_sd / si.ni).asinh());
    let cj = depletion_capacitance(EPS_SI, sn.n_sub, sn.n_sd, vbi, -0.5 * vdd, vt) * sn.width * sn.l_sd;
    let c_node = (cg_n + ratio * cg_p) + (1.0 + ratio) * cj;
    println!("  ring oscillator, 5 stages: node capacitance {:.3} fF = gates (TCAD charge, n {:.3} + p {:.3}×{ratio:.2} fF) + drain junctions ({:.3} fF per µm², ×{:.2})", c_node * 1e15, cg_n * 1e15, cg_p * 1e15, cj * 1e15, 1.0 + ratio);
    let (period, tpd, _) = ring_oscillator(&tn, &tp, vdd, 5, c_node, 0.05e-12, 2e-9);
    let i_avg = 0.5 * (ion_n + ratio * ion_p);
    println!("    period {:.2} ps ({:.2} GHz), delay per stage {:.2} ps (C·V_dd/(2 I_on,avg) estimate {:.2} ps)", period * 1e12, 1e-9 / period, tpd * 1e12, c_node * vdd / (2.0 * i_avg) * 1e12);
}

/// Inversion-layer mobility: extract μ_eff(E_eff) from long-channel
/// devices at several substrate dopings; with surface scattering the
/// curves collapse onto one universal curve at high field.
fn mobility() {
    let si = Silicon::new(300.0);
    let vd = 0.02;
    println!("inversion-layer mobility: nMOS, gate 2 µm, t_ox 4 nm, V_d {vd} V, Masetti bulk + Lombardi surface (phonon + roughness) mobility");
    println!("  μ_eff = I_d L/(W Q_inv V_d) at mid-channel; E_eff = (|Q_dep| + |Q_inv|/2)/ε_si");
    for screening in [false, true] {
    println!("  impurity scattering {}:", if screening { "re-screened by the local carriers (Brooks–Herring unfold of Masetti)" } else { "Masetti as fitted to neutral bulk (unscreened by the inversion layer)" });
    let mut curves = Vec::new();
    for n_sub in [3e22, 1e23, 3e23] {
        let mut spec = MosfetSpec::scaled(2e-6, false);
        spec.n_sub = n_sub;
        let mut dev = spec.build(si.clone());
        dev.lombardi = true;
        dev.screening = screening;
        dev.solve().expect("equilibrium");
        let (d, g) = (dev.contact("drain"), dev.contact("gate"));
        dev.ramp(d, vd, 0.02).expect("drain");
        let m = dev.mesh.clone();
        let xm = 0.5 * (spec.edges().0 + spec.edges().1);
        let i = (0..m.nx()).min_by(|&a, &b| (m.x[a] - xm).abs().partial_cmp(&(m.x[b] - xm).abs()).unwrap()).unwrap();
        let dxd = 0.5 * (m.x[i + 1] - m.x[i - 1]);
        let mut rows = Vec::new();
        let mut vg = 0.6;
        while vg <= 2.51 {
            dev.ramp(g, vg, 0.1).expect("gate");
            let id = tcad_id(&dev);
            let (mut qinv, mut qdep) = (0.0, 0.0);
            for j in 0..m.ny() {
                let k = m.node(i, j);
                if !dev.is_semiconductor(k) {
                    continue;
                }
                let v = dev.semiconductor_volume(k) / dxd;
                qinv += Q * dev.n(k) * v;
                qdep += Q * (dev.na[k] - dev.nd[k] - dev.p(k)).max(0.0) * v;
            }
            let mu = id * spec.l_eff() / (spec.width * qinv * vd);
            let eeff = (qdep + 0.5 * qinv) / EPS_SI;
            rows.push((eeff, mu * 1e4, qinv));
            vg += 0.1;
        }
        curves.push((n_sub, rows));
    }
    println!("  {:>12} {:>32} {:>32} {:>32}", "", "N_a 3e16 cm⁻³", "1e17", "3e17");
    println!("  {:>12} {:>32} {:>32} {:>32}", "", "E_eff (MV/cm)  μ_eff (cm²/Vs)", "E_eff  μ_eff", "E_eff  μ_eff");
    let nrow = curves[0].1.len();
    for r in (0..nrow).step_by(2) {
        let cells: Vec<String> = curves.iter().map(|c| format!("{:10.3} {:10.1}", c.1[r].0 * 1e-8, c.1[r].1)).collect();
        println!("  {:>12} {:>32} {:>32} {:>32}", format!("V_g {:.1}", 0.6 + 0.1 * r as f64), cells[0], cells[1], cells[2]);
    }
    // universality: compare the curves at common E_eff (log interpolation)
    let interp = |c: &Vec<(f64, f64, f64)>, e: f64| -> Option<f64> {
        c.windows(2).find(|w| (w[0].0 - e) * (w[1].0 - e) <= 0.0).map(|w| {
            let t = (e.ln() - w[0].0.ln()) / (w[1].0.ln() - w[0].0.ln());
            (w[0].1.ln() + t * (w[1].1.ln() - w[0].1.ln())).exp()
        })
    };
    println!("  universality: spread of μ_eff between the three dopings at equal E_eff");
    for e_mv in [0.3, 0.5, 0.7, 1.0] {
        let vals: Vec<f64> = curves.iter().filter_map(|c| interp(&c.1, e_mv * 1e8)).collect();
        if vals.len() >= 2 {
            let (lo, hi) = (vals.iter().cloned().fold(f64::INFINITY, f64::min), vals.iter().cloned().fold(0.0, f64::max));
            println!("    E_eff {e_mv:.1} MV/cm: μ_eff {:.0}–{:.0} cm²/Vs ({} curves, spread {:.1}%)", lo, hi, vals.len(), 100.0 * (hi / lo - 1.0));
        }
    }
    }
}

/// MOS capacitor under linear gate ramps: transient gate current / ramp
/// rate = capacitance; fast ramps give deep depletion (no time to generate
/// inversion electrons), slow ramps the low-frequency curve.
fn transient_moscap() {
    let si = Silicon::new(300.0);
    let vt = si.vt();
    let (na, tox, lsi) = (1e23, 10e-9, 1e-6);
    let offset = 0.5 * band_gap(300.0);
    let mos = Mos { na, ni: si.ni, vt, eps_s: EPS_SI, cox: EPS_OX / tox, vfb: -vt * (0.5 * na / si.ni).asinh() - offset };
    let x = graded(&[(0.0, 0.5e-9), (tox, 0.05e-9), (tox + 20e-9, 1e-9), (tox + lsi, 20e-9)]);
    let mesh = Mesh::line(x, move |x| if x < tox { Material::Oxide } else { Material::Silicon });
    let last = mesh.nx() - 1;
    println!("MOS capacitor under gate ramps −2 → +2 V (n⁺ poly / 10 nm SiO₂ / p 1e17 cm⁻³); C = I_g/(dV/dt), BDF2, 200 steps per ramp");
    // (slower ramps fold to the quasi-static solve: at 1 mV/s the
    // displacement current, ~3.5e-6 A/m², is ~1e-14 of the accumulation-
    // layer drift/diffusion fluxes — below double precision)
    let ramps = [(1e6, "1 V/µs"), (0.1, "0.1 V/s")];
    let mut results = Vec::new();
    for &(rate, label) in &ramps {
        let contacts = vec![
            Contact { name: "gate".into(), nodes: line_nodes(&mesh, 0), kind: ContactKind::Gate { offset }, v: 0.0 },
            Contact { name: "body".into(), nodes: line_nodes(&mesh, last), kind: ContactKind::Ohmic, v: 0.0 },
        ];
        let mut dev = Device::new(mesh.clone(), si.clone(), move |_, _| (0.0, na), contacts, 1.0);
        dev.solve().expect("equilibrium");
        let (g, b) = (dev.contact("gate"), dev.contact("body"));
        dev.ramp(g, -2.0, 0.2).expect("accumulation");
        dev.begin_transient();
        let steps = 200;
        let dt = 4.0 / rate / steps as f64;
        let mut rows = Vec::new();
        let mut cons: f64 = 0.0;
        let mut substeps = 0;
        for s in 1..=steps {
            let t0 = dt * (s - 1) as f64;
            substeps += dev.advance(dt, &|d: &mut Device, f: f64| d.contacts[g].v = -2.0 + rate * (t0 + f * dt)).expect("transient step");
            if s >= 2 {
                let (ig, ib) = (dev.transient_current(g), dev.transient_current(b));
                cons = cons.max((ig + ib).abs() / ig.abs().max(1e-300));
                rows.push((dev.contacts[g].v, ig / rate / mos.cox));
            }
        }
        println!("  ramp {label}: total current conserved to {cons:.1e} (|I_g + I_b|/|I_g|), {substeps} steps (adaptive)");
        results.push(rows);
    }
    let cexact = |f: &dyn Fn(f64) -> f64, v: f64| (f(v + 1e-3) - f(v - 1e-3)) / 2e-3 / mos.cox;
    // the quasi-static fold: steady solves (gate charge difference)
    let contacts = vec![
        Contact { name: "gate".into(), nodes: line_nodes(&mesh, 0), kind: ContactKind::Gate { offset }, v: 0.0 },
        Contact { name: "body".into(), nodes: line_nodes(&mesh, last), kind: ContactKind::Ohmic, v: 0.0 },
    ];
    let mut qs = Device::new(mesh.clone(), si.clone(), move |_, _| (0.0, na), contacts, 1.0);
    qs.solve().expect("equilibrium");
    let gq = qs.contact("gate");
    let mut c_qs = |v: f64| -> f64 {
        qs.ramp(gq, v - 1e-3, 0.1).expect("qs");
        let a = qs.charge(gq);
        qs.ramp(gq, v + 1e-3, 0.1).expect("qs");
        (qs.charge(gq) - a) / 2e-3 / mos.cox
    };
    println!("  {:>7} {:>12} {:>12} {:>14} {:>14} {:>14}", "V_g", "1 V/µs", "0.1 V/s", "QS fold", "deep depl.", "low-freq.");
    for target in [-1.5, -1.0, -0.6, -0.2, 0.2, 0.5, 1.0, 1.5, 1.9] {
        let pick = |r: &Vec<(f64, f64)>| r.iter().min_by(|a, b| (a.0 - target).abs().partial_cmp(&(b.0 - target).abs()).unwrap()).map(|x| (x.0, x.1)).unwrap();
        let (v0, c0) = pick(&results[0]);
        let (_, c1) = pick(&results[1]);
        let c2 = c_qs(v0);
        // exact curves at the node voltage (the BDF2 current is at the step end)
        println!("  {:7.3} {:12.5} {:12.5} {:14.5} {:14.5} {:14.5}", v0, c0, c1, c2, cexact(&|v| mos.qg_deep(v), v0), cexact(&|v| mos.qg(v), v0));
    }
}

/// Mixed-mode transient: an inverter of two drift–diffusion devices
/// driving a load capacitance, on the time grid `times`; KCL at the output
/// solved each variable-step BDF2 step by a secant on V_out started from
/// the previous step's slope (each evaluation advances both devices from
/// their committed histories), bracketing as a fallback.
#[allow(clippy::too_many_arguments)]
fn tcad_switching(n: &mut Device, p: &mut Device, ps: f64, vdd: f64, cl: f64, vin: &dyn Fn(f64) -> f64, times: &[f64]) -> (Vec<(f64, f64, f64)>, usize) {
    let (vo0, _) = tcad_inverter(n, p, ps, vdd, vin(0.0));
    let (nd, ng) = (n.contact("drain"), n.contact("gate"));
    let (pd, pg) = (p.contact("drain"), p.contact("gate"));
    n.begin_transient();
    p.begin_transient();
    let mut vhist = vec![vo0];
    let mut out = vec![(times[0], vin(times[0]), vo0)];
    let mut slope = f64::NAN;
    let mut evals = 0;
    let mut dt_prev = 0.0;
    // one implicit step to time t from the last committed level; None if a
    // device step fails or KCL cannot be bracketed (the caller subdivides)
    #[allow(clippy::too_many_arguments)]
    fn kcl_step(n: &mut Device, p: &mut Device, ids: (usize, usize, usize, usize), ps: f64, vdd: f64, cl: f64, vi: f64, dt: f64, dt_prev: f64, vhist: &[f64], slope: &mut f64, evals: &mut usize) -> Option<f64> {
        let (nd, ng, pd, pg) = ids;
        let a: [f64; 3] = if vhist.len() >= 2 {
            let w = dt / dt_prev;
            [(1.0 + 2.0 * w) / (1.0 + w), -(1.0 + w), w * w / (1.0 + w)]
        } else {
            [1.0, -1.0, 0.0]
        };
        let hist: f64 = vhist.iter().enumerate().map(|(l, v)| a[l + 1] * v).sum();
        let mut eval = |v: f64, n: &mut Device, p: &mut Device| -> f64 {
            *evals += 1;
            n.contacts[ng].v = vi;
            n.contacts[nd].v = v;
            p.contacts[pg].v = vi - vdd;
            p.contacts[pd].v = v - vdd;
            if n.step(dt).is_none() || p.step(dt).is_none() {
                n.reject_step();
                p.reject_step();
                return f64::NAN;
            }
            cl * (a[0] * v + hist) / dt + n.transient_current(nd) + ps * p.transient_current(pd)
        };
        let v_last = vhist[0];
        let mut v0 = if vhist.len() >= 2 { v_last + (v_last - vhist[1]) * dt / dt_prev } else { v_last };
        let mut f0 = eval(v0, n, p);
        let mut converged = false;
        if f0.is_finite() && slope.is_finite() && *slope > 0.0 {
            for _ in 0..8 {
                let v1 = v0 - f0 / *slope;
                let f1 = eval(v1, n, p);
                if !f1.is_finite() {
                    break;
                }
                if (v1 - v0).abs() < 1e-6 {
                    converged = true;
                    v0 = v1;
                    break;
                }
                if (v1 - v0).abs() > 1e-12 && f1 != f0 {
                    let s2 = (f1 - f0) / (v1 - v0);
                    if s2 > 0.0 {
                        *slope = s2;
                    }
                }
                v0 = v1;
                f0 = f1;
            }
        }
        if !converged {
            // bracket (f increases with V_out), then Illinois
            let (mut lo, mut hi) = (v_last - 0.05, v_last + 0.05);
            let (mut flo, mut fhi) = (eval(lo, n, p), eval(hi, n, p));
            let mut grow = 0;
            while flo.is_finite() && fhi.is_finite() && flo * fhi > 0.0 && grow < 6 {
                if flo > 0.0 {
                    lo -= 0.1 * (1 << grow) as f64;
                    flo = eval(lo, n, p);
                } else {
                    hi += 0.1 * (1 << grow) as f64;
                    fhi = eval(hi, n, p);
                }
                grow += 1;
            }
            if !(flo.is_finite() && fhi.is_finite() && flo * fhi <= 0.0) {
                n.reject_step();
                p.reject_step();
                return None;
            }
            *slope = (fhi - flo) / (hi - lo);
            let mut side = 0;
            for _ in 0..60 {
                v0 = (lo * fhi - hi * flo) / (fhi - flo);
                let fc = eval(v0, n, p);
                if !fc.is_finite() {
                    n.reject_step();
                    p.reject_step();
                    return None;
                }
                if (hi - lo).abs() < 1e-6 {
                    break;
                }
                if fc * fhi > 0.0 {
                    hi = v0;
                    fhi = fc;
                    if side == -1 {
                        flo *= 0.5;
                    }
                    side = -1;
                } else {
                    lo = v0;
                    flo = fc;
                    if side == 1 {
                        fhi *= 0.5;
                    }
                    side = 1;
                }
            }
        }
        // commit the converged step
        if !eval(v0, n, p).is_finite() {
            return None;
        }
        n.accept_step();
        p.accept_step();
        Some(v0)
    }
    let ids = (nd, ng, pd, pg);
    let mut t_now = times[0];
    for &t_target in &times[1..] {
        // advance to t_target, halving on failure
        let mut h = t_target - t_now;
        while t_now < t_target - 1e-21 {
            let t = (t_now + h).min(t_target);
            let dt = t - t_now;
            match kcl_step(n, p, ids, ps, vdd, cl, vin(t), dt, dt_prev, &vhist, &mut slope, &mut evals) {
                Some(v) => {
                    vhist.insert(0, v);
                    vhist.truncate(2);
                    dt_prev = dt;
                    t_now = t;
                    out.push((t, vin(t), v));
                    h = (2.0 * h).min(t_target - t_now).max(1e-21);
                }
                None => {
                    h *= 0.5;
                    assert!(h > (t_target - t_now) * 1e-4 && h > 1e-18, "mixed-mode step failed at t = {t_now:e}");
                }
            }
        }
    }
    n.end_transient();
    p.end_transient();
    (out, evals)
}

/// CMOS inverter switching: quasi-static table model (currents + charges
/// folded from TCAD) vs mixed-mode transient TCAD, for input rise times
/// from slow to comparable with the channel transit time.
fn switching() {
    let si = Silicon::new(300.0);
    let vdd = 1.2;
    let l = 150e-9;
    let (sn, sp) = (MosfetSpec::scaled(l, false), MosfetSpec::scaled(l, true));
    let mut n = sn.build(si.clone());
    let mut p = sp.build(si.clone());
    for d in [&mut n, &mut p] {
        d.velocity_saturation = true;
        d.lombardi = true;
        d.screening = true;
        d.solve().expect("equilibrium");
    }
    let vgs: Vec<f64> = (0..=12).map(|i| 0.1 * i as f64).collect();
    let vds: Vec<f64> = vec![0.0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2];
    let t0 = std::time::Instant::now();
    let fold_dev = |d: &mut Device, pol: f64| {
        let gs: Vec<f64> = vgs.iter().map(|v| pol * v).collect();
        let ds: Vec<f64> = vds.iter().map(|v| pol * v).collect();
        let (id, q) = sample_with_charges(d, &gs, &ds);
        (Table::new(vgs.clone(), vds.clone(), id, pol), QTable { vgs: vgs.clone(), vds: vds.clone(), q, polarity: pol, scale: 1.0 })
    };
    let (tn, qn) = fold_dev(&mut n, 1.0);
    let (mut tp, mut qp) = fold_dev(&mut p, -1.0);
    // symmetric sizing
    let (mut lo, mut hi) = (0.5, 6.0);
    for _ in 0..50 {
        let mid = 0.5 * (lo + hi);
        tp.scale = mid;
        if inverter_vout(&tn, &tp, vdd, 0.5 * vdd) < 0.5 * vdd {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let ratio = 0.5 * (lo + hi);
    tp.scale = ratio;
    qp.scale = ratio;
    let ion = tn.id(vdd, vdd);
    println!("CMOS inverter switching, gate 150 nm, V_dd {vdd} V; mobility: Masetti + screening + Lombardi + velocity saturation");
    println!("  folded I_d and quasi-static (Q_g, Q_d) tables in {:.0} s; I_on,n {:.1} µA/µm, W_p/W_n {ratio:.2}", t0.elapsed().as_secs_f64(), ion * 1e6);
    // channel transit time: L²/(μ_eff V), μ_eff from the linear-region table
    let mu_eff = tn.id(vdd, 0.025) * sn.l_eff() / (0.025 * EPS_OX / sn.tox * (vdd - 0.35) * sn.width);
    let tau_tr = sn.l_eff() * sn.l_eff() / (mu_eff * (vdd - 0.35));
    println!("  channel transit time L²/(μ_eff (V_dd − V_t)) ≈ {:.2} ps (μ_eff {:.0} cm²/Vs): the quasi-static fold needs input edges ≫ this", tau_tr * 1e12, mu_eff * 1e4);
    let cl = 2e-15;
    let rises: Vec<f64> = std::env::var("SW_RISE").map(|v| v.split(',').map(|x| x.parse::<f64>().unwrap() * 1e-12).collect()).unwrap_or(vec![50e-12, 5e-12, 1e-12]);
    for tr in rises {
        let delay = 2e-12;
        let vin = move |t: f64| if t < delay { 0.0 } else if t < delay + tr { vdd * (t - delay) / tr } else { vdd };
        // time grid: tr/20 through the edge (+3 ps), at most 0.25 ps (or
        // the edge step if coarser) through the output fall, then 1 ps
        let fine = (tr / 20.0).min(0.5e-12);
        let mid = fine.max(0.25e-12);
        let t_edge_end = delay + tr + 3e-12;
        let t_fine_end = delay + tr + 25e-12;
        let t_end = t_fine_end + 20e-12;
        let mut times = vec![0.0];
        while *times.last().unwrap() < t_end - 1e-18 {
            let t = *times.last().unwrap();
            let h = if t < t_edge_end { fine } else if t < t_fine_end { mid } else { 1e-12 };
            times.push((t + h).min(t_end));
        }
        let tq = std::time::Instant::now();
        let w_qs = inverter_transient_qs(&tn, &qn, &tp, &qp, vdd, cl, &vin, &times);
        let t_qs = tq.elapsed().as_secs_f64();
        let t1 = std::time::Instant::now();
        let (w_tc, evals) = tcad_switching(&mut n, &mut p, ratio, vdd, cl, &vin, &times);
        let tin = delay + 0.5 * tr;
        let (dq, dtc) = (crossing(&w_qs, 0.5 * vdd, true, false) - tin, crossing(&w_tc, 0.5 * vdd, true, false) - tin);
        let fall = |w: &[(f64, f64, f64)]| crossing(w, 0.1 * vdd, true, false) - crossing(w, 0.9 * vdd, true, false);
        let vmin = |w: &[(f64, f64, f64)]| w.iter().map(|x| x.2).fold(f64::INFINITY, f64::min);
        let vmax = |w: &[(f64, f64, f64)]| w.iter().map(|x| x.2).fold(f64::NEG_INFINITY, f64::max);
        println!(
            "  input rise {:5.1} ps (edge/transit {:5.1}): delay t_pHL table {:6.2} ps, TCAD {:6.2} ps ({:+.1}%); fall 90→10% table {:6.2} ps, TCAD {:6.2} ps; V_out range table [{:+.3}, {:.3}], TCAD [{:+.3}, {:.3}] V; {} steps, {} device-pair evaluations; table {:.3} s, TCAD {:.0} s",
            tr * 1e12,
            tr / tau_tr,
            dq * 1e12,
            dtc * 1e12,
            100.0 * (dq / dtc - 1.0),
            fall(&w_qs) * 1e12,
            fall(&w_tc) * 1e12,
            vmin(&w_qs),
            vmax(&w_qs),
            vmin(&w_tc),
            vmax(&w_tc),
            w_tc.len() - 1,
            evals,
            t_qs,
            t1.elapsed().as_secs_f64()
        );
    }
}

fn thermal() {
    // 1. heat solver vs the exact Kirchhoff-transform solution: a silicon
    //    slab of depth L with uniform heating, insulated top, bottom at T₀,
    //    k = 148 (300/T)^1.3 → U = ∫k dT = q(L² − y²)/2 exactly
    println!("1. heat solver: slab with uniform heating 2e14 W/m³, k_Si(T) ∝ T^−1.3 (Picard), bottom 300 K, exact by Kirchhoff transform");
    let (lth, qv, t0): (f64, f64, f64) = (10e-6, 2e14, 300.0);
    let exact = |y: f64| {
        let u = 0.5 * qv * (lth * lth - y * y);
        (t0.powf(-0.3) - 0.3 * u / (148.0 * 300f64.powf(1.3))).powf(-1.0 / 0.3)
    };
    let mut prev: Option<f64> = None;
    for n in [10, 20, 40, 80] {
        let mesh = Mesh::new(uniform(0.0, 1e-6, 2), uniform(0.0, lth, n), |_, _| Material::Silicon);
        let dev = Device::new(mesh.clone(), Silicon::new(300.0), |_, _| (0.0, 1e21), vec![], 1e-6);
        let q: Vec<f64> = (0..mesh.nodes()).map(|k| qv * dev.semiconductor_volume(k) * dev.depth).collect();
        let h = DeviceHeat::solve(&dev, &q, t0);
        let err = (0..mesh.nodes()).map(|k| (h.t[k] - exact(mesh.y[mesh.ij(k).1])).abs()).fold(0.0, f64::max);
        println!("   {n:3} cells: T_max {:.4} K (exact {:.4}), max error {err:.3e} K{}, {} Picard iterations", h.t.iter().cloned().fold(0.0, f64::max), exact(0.0), prev.map_or(String::new(), |p: f64| format!(", order {:.2}", (p / err).log2())), h.picard_iterations);
        prev = Some(err);
    }
    // 2. diode saturation current vs temperature (Shockley, exact test)
    println!("2. diode J₀(T) without depletion recombination vs Shockley with n_i(T), μ(T), V_t(T):");
    for t in [250.0, 300.0, 350.0, 400.0] {
        let si = Silicon::new(t);
        let vt = si.vt();
        let (na, nd) = (1e23, 1e22);
        let mun = si.mob_n.mobility(na, t);
        let mup = si.mob_p.mobility(nd, t);
        let (tn, tp) = (si.tau_n(na), si.tau_p(nd));
        let ln = (mun * vt * tn).sqrt();
        let lp = (mup * vt * tp).sqrt();
        let (lps, lns) = (10.0 * ln, 10.0 * lp);
        let xj = lps;
        let x = graded(&[(0.0, lps / 200.0), (xj, 0.5e-9), (xj + lns, lns / 200.0)]);
        let mesh = Mesh::line(x, |_| Material::Silicon);
        let last = mesh.nx() - 1;
        let contacts = vec![
            Contact { name: "anode".into(), nodes: line_nodes(&mesh, 0), kind: ContactKind::Ohmic, v: 0.0 },
            Contact { name: "cathode".into(), nodes: line_nodes(&mesh, last), kind: ContactKind::Ohmic, v: 0.0 },
        ];
        let mut dev = Device::new(mesh.clone(), si.clone(), move |x, _| if x < xj { (0.0, na) } else { (nd, 0.0) }, contacts, 1.0);
        let vbi = vt * ((0.5 * na / si.ni).asinh() + (0.5 * nd / si.ni).asinh());
        let v = 0.25;
        let w = (2.0 * EPS_SI * (vbi - v) / Q * (na + nd) / (na * nd)).sqrt();
        let (xp, xn) = (w * nd / (na + nd), w * na / (na + nd));
        dev.recomb_scale = Some((0..mesh.nodes()).map(|k| { let x = mesh.x[mesh.ij(k).0]; if x > xj - xp && x < xj + xn { 0.0 } else { 1.0 } }).collect());
        dev.solve().expect("equilibrium");
        let an = dev.contact("anode");
        dev.ramp(an, v, 0.05).expect("forward");
        let (a, b) = dev.current(an);
        let j0 = shockley_j0(si.ni, vt, na, nd, mun, tn, mup, tp);
        let j0_dd = (a + b) / ((v / vt).exp() - 1.0);
        println!("   T {t:.0} K: J₀ {:.4e} A/cm² (Shockley {:.4e}, ratio {:.5}); n_i {:.3e} cm⁻³", j0_dd * 1e-4, j0 * 1e-4, j0_dd / j0, si.ni * 1e-6);
    }
    // 3. the 150 nm nMOS: leakage vs temperature, energy balance, self-heating
    let vdd = 1.2;
    let spec = MosfetSpec::scaled(150e-9, false);
    let mut dev = spec.build(Silicon::new(300.0));
    dev.velocity_saturation = true;
    dev.lombardi = true;
    dev.screening = true;
    dev.solve().expect("equilibrium");
    let (d, g) = (dev.contact("drain"), dev.contact("gate"));
    println!("3. nMOS 150 nm (scaled technology, full mobility model): leakage vs temperature at V_g = 0, V_d = {vdd} V");
    dev.ramp(d, vdd, 0.1).expect("drain");
    let mut ioff_t = Vec::new();
    for t in [300.0, 325.0, 350.0, 375.0, 400.0] {
        dev.set_temperature(t);
        dev.solve().expect("temperature");
        let i = tcad_id(&dev);
        ioff_t.push((t, i));
        println!("   T {t:.0} K: I_off {:.3e} A/µm", i);
    }
    let fit = (ioff_t[4].1 / ioff_t[0].1).ln() / 100.0;
    println!("   leakage doubles every {:.1} K", 2f64.ln() / fit);
    for t in [375.0, 350.0, 325.0, 300.0] {
        dev.set_temperature(t);
        dev.solve().expect("back to 300 K");
    }
    dev.ramp(g, vdd, 0.1).expect("gate on");
    let ion = tcad_id(&dev);
    let heat = dev.heat();
    let ptot: f64 = heat.iter().sum();
    println!("   on state (V_g = V_d = {vdd} V): I_on {:.1} µA/µm; dissipated heat Σ H {:.6e} W vs I_d V_d {:.6e} W (energy balance {:.1e})", ion * 1e6, ptot, ion * vdd, (ptot / (ion * vdd) - 1.0).abs());
    // electro-thermal: device heat → 2D temperature → device at its mean
    // temperature, iterated; the bottom of the 0.5 µm deep domain at T_die
    let t_die = 300.0;
    let mut t_dev = t_die;
    for it in 0..10 {
        let q = dev.heat();
        let h = DeviceHeat::solve(&dev, &q, t_die);
        let (mean, spread, tmax) = h.summary(&q);
        println!("   electro-thermal iteration {it}: I_on {:.2} µA/µm, mean device temperature {mean:.3} K, spread {spread:.3} K, peak {tmax:.3} K", tcad_id(&dev) * 1e6);
        if (mean - t_dev).abs() < 1e-3 {
            println!("   isothermal-device fold indicator: spread {spread:.3} K ≪ T — valid (unfold to non-isothermal drift–diffusion if the spread reaches ~10 K)");
            break;
        }
        t_dev = mean;
        dev.set_temperature(t_dev);
        dev.solve().expect("electro-thermal");
    }
    // 4. coolant: laminar microchannel (Graetz problem)
    println!("4. coolant: laminar flow between parallel plates, both walls at constant flux (Graetz, marched; axial conduction folded for Pe ≫ 1)");
    let nu = graetz(800, 0.5, 8000);
    for xs in [1e-5f64, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.5] {
        let r = nu.iter().min_by(|a, b| (a.0.ln() - xs.ln()).abs().partial_cmp(&(b.0.ln() - xs.ln()).abs()).unwrap()).unwrap();
        let lev = 1.490 * r.0.powf(-1.0 / 3.0);
        println!("   x* {:.2e}: Nu_x {:.4} (entrance asymptote 1.490 x*^−1/3 = {:.4}; fully developed 140/17 = {:.5})", r.0, r.1, lev, 140.0 / 17.0);
    }
    // 5. chip energy balance with TCAD leakage: thermal runaway
    let leak = |t: f64| -> f64 {
        // log-linear interpolation of I_off(T) per µm, extrapolated
        let (t0_, i0) = ioff_t[0];
        let s = (ioff_t[4].1 / ioff_t[0].1).ln() / (ioff_t[4].0 - t0_);
        i0 * (s * (t - t0_)).exp()
    };
    let width_total = 2e9 * 1.0; // µm of transistor width on the chip (off-state leakage)
    let p_dyn = 40.0;
    let k_w = 0.6;
    let dh = 200e-6;
    let h_fd = 8.235 * k_w / dh;
    println!("5. chip: 1 cm², {:.0e} µm of off transistors leaking I_off(T)·V_dd, P_dyn {p_dyn} W, coolant 300 K; microchannels D_h 200 µm, h = 8.235 k/D_h = {:.0} W/(m² K)", width_total, h_fd);
    for (label, r_th) in [("microchannel cooler", 1.0 / (h_fd * 1e-4) + 0.05), ("air heat sink", 0.6), ("poor heat sink", 1.2)] {
        let pts = runaway(300.0, r_th, p_dyn, &|t| width_total * leak(t) * vdd, 250.0, 700.0);
        let desc: Vec<String> = pts.iter().map(|(t, st)| format!("{t:.1} K ({})", if *st { "stable" } else { "unstable" })).collect();
        println!("   {label}: R_th {r_th:.3} K/W → fixed points {}", if desc.is_empty() { "none: thermal runaway".to_string() } else { desc.join(", ") });
    }
}

fn quantum() {
    let stack = |na: f64, tox: f64| MosStack { na, tox, t: 300.0, depth: 150e-9, lq: 15e-9, m_dg: 0.3 };
    // inversion-charge threshold: 1e12 electrons/cm²
    let q_th = Q * 1e16;
    let vt_of = |st: &MosStack, mode: Electrons| -> f64 {
        let (mut lo, mut hi) = (-0.5, 2.0);
        for _ in 0..40 {
            let mid = 0.5 * (lo + hi);
            if st.solve(mid, mode).q_inv < q_th {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    };
    // calibrate the density-gradient mass on Schrödinger–Poisson: stack A,
    // V_g = 1 V
    let a = stack(1e24, 2e-9);
    let q_sp = a.solve(1.0, Electrons::SchrodingerPoisson).q_inv;
    let (mut lo, mut hi): (f64, f64) = (0.05, 3.0);
    for _ in 0..40 {
        let mid = (lo * hi).sqrt();
        let mut st = a.clone();
        st.m_dg = mid;
        // heavier mass → weaker quantum correction → more charge
        if st.solve(1.0, Electrons::DensityGradient).q_inv < q_sp {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let m_dg = (lo * hi).sqrt();
    println!("quantum confinement in the inversion layer (n⁺ poly / SiO₂ / p-Si, 300 K)");
    println!("  density-gradient mass calibrated on Schrödinger–Poisson (N_a 1e18 cm⁻³, t_ox 2 nm, V_g 1 V): m_dg = {m_dg:.4} m₀");
    for (label, na, tox) in [("A: N_a 1e18, t_ox 2 nm (calibration stack)", 1e24, 2e-9), ("B: N_a 3e17, t_ox 3 nm", 3e23, 3e-9), ("C: N_a 3e18, t_ox 1.5 nm", 3e24, 1.5e-9)] {
        let mut st = stack(na, tox);
        st.m_dg = m_dg;
        println!("  {label}");
        println!("    {:>5} {:>12} {:>12} {:>12} {:>9} {:>9} {:>9} {:>8} {:>10}", "V_g", "Q_inv cl.", "Q_inv SP", "Q_inv DG", "DG err", "z̄ cl.", "z̄ SP", "z̄ DG", "E₀ (SP)");
        for vg in [0.6, 1.0, 1.4, 1.8] {
            let (c, sp, dg) = (st.solve(vg, Electrons::Classical), st.solve(vg, Electrons::SchrodingerPoisson), st.solve(vg, Electrons::DensityGradient));
            println!(
                "    {vg:5.1} {:9.3e} C/m² {:9.3e} {:9.3e} {:+8.2}% {:6.2} nm {:6.2} nm {:6.2} nm {:7.0} meV",
                c.q_inv,
                sp.q_inv,
                dg.q_inv,
                100.0 * (dg.q_inv / sp.q_inv - 1.0),
                c.centroid * 1e9,
                sp.centroid * 1e9,
                dg.centroid * 1e9,
                sp.e0 * 1e3
            );
        }
        let (vc, vs, vd) = (vt_of(&st, Electrons::Classical), vt_of(&st, Electrons::SchrodingerPoisson), vt_of(&st, Electrons::DensityGradient));
        // inversion capacitance at high field → effective oxide thickness
        let cinv = |mode: Electrons| (st.solve(1.8, mode).q_inv - st.solve(1.6, mode).q_inv) / 0.2;
        let (cc, cs) = (cinv(Electrons::Classical), cinv(Electrons::SchrodingerPoisson));
        let dt_eff = EPS_OX * (1.0 / cs - 1.0 / cc);
        let indicator = (st.solve(1.2, Electrons::Classical).q_inv / st.solve(1.2, Electrons::SchrodingerPoisson).q_inv - 1.0).abs();
        println!(
            "    threshold (1e12 cm⁻²): classical {vc:.4} V, SP {vs:.4} V (shift {:+.1} mV), DG {vd:.4} V ({:+.1} mV from SP); inversion C: SP/classical {:.3} → +{:.2} nm effective oxide",
            1e3 * (vs - vc),
            1e3 * (vd - vs),
            cs / cc,
            dt_eff * 1e9
        );
        println!("    fold indicator |Q_cl/Q_SP − 1| at V_g 1.2 V = {:.1}%: {}", 100.0 * indicator, if indicator > 0.01 { "classical drift–diffusion invalid here — unfold to density gradient (or Schrödinger–Poisson)" } else { "classical fold valid" });
    }
}

fn quantum2d() {
    let m_dg = 0.2013; // calibrated on Schrödinger–Poisson (semiconductor quantum)
    let si = Silicon::new(300.0);
    let vt = si.vt();
    // 1. a 2D MOS capacitor (uniform in x) vs the 1D density-gradient and
    //    Schrödinger–Poisson solutions of the same stack
    let (na, tox) = (1e24, 2e-9);
    let mut y = graded(&[(-tox, tox / 8.0), (0.0, 0.02e-9)]);
    y.pop();
    y.extend(graded(&[(0.0, 0.02e-9), (2e-9, 0.1e-9), (15e-9, 0.3e-9), (150e-9, 5e-9)]));
    let mesh = Mesh::new(uniform(0.0, 1e-7, 2), y, |_, y| if y < 0.0 { Material::Oxide } else { Material::Silicon });
    let ny = mesh.ny();
    let gate: Vec<usize> = (0..mesh.nx()).map(|i| mesh.node(i, 0)).collect();
    let body: Vec<usize> = (0..mesh.nx()).map(|i| mesh.node(i, ny - 1)).collect();
    let offset = 0.5 * band_gap(300.0);
    let contacts = vec![
        Contact { name: "gate".into(), nodes: gate, kind: ContactKind::Gate { offset }, v: 0.0 },
        Contact { name: "body".into(), nodes: body, kind: ContactKind::Ohmic, v: 0.0 },
    ];
    let mk = |dg: bool| {
        let mut d = Device::new(mesh.clone(), si.clone(), move |_, _| (0.0, na), contacts.clone(), 1e-7);
        if dg {
            d.dg_mass = Some(m_dg);
        }
        d.solve_quantum().expect("equilibrium");
        d
    };
    let (mut dcl, mut ddg) = (mk(false), mk(true));
    let st = MosStack { na, tox, t: 300.0, depth: 150e-9, lq: 15e-9, m_dg };
    println!("quasi-2D density gradient (columns along y) inside the device solver; m_dg {m_dg} m₀ from Schrödinger–Poisson");
    println!("  1. MOS capacitor N_a 1e18, t_ox 2 nm, as a 2D device vs the 1D solutions (inversion charge, C/m²)");
    let qinv = |d: &Device| -> f64 {
        let m = &d.mesh;
        let i = 1;
        let dxd = 0.5 * (m.x[2] - m.x[0]);
        let nb = si.ni * (-(0.5 * na / si.ni).asinh()).exp();
        (0..m.ny()).map(|j| m.node(i, j)).filter(|&k| d.is_semiconductor(k)).map(|k| Q * (d.n(k) - nb).max(0.0) * d.semiconductor_volume(k) / dxd).sum()
    };
    for vg in [0.6, 1.0, 1.4, 1.8] {
        for d in [&mut dcl, &mut ddg] {
            let g = d.contact("gate");
            let v0 = d.contacts[g].v;
            let steps = ((vg - v0).abs() / 0.1).ceil().max(1.0) as usize;
            for s in 1..=steps {
                d.contacts[g].v = v0 + (vg - v0) * s as f64 / steps as f64;
                d.solve_quantum().expect("gate bias");
            }
        }
        let (a, b) = (qinv(&dcl), qinv(&ddg));
        let (c1, d1, s1) = (st.solve(vg, Electrons::Classical).q_inv, st.solve(vg, Electrons::DensityGradient).q_inv, st.solve(vg, Electrons::SchrodingerPoisson).q_inv);
        println!("     V_g {vg:.1}: 2D classical {a:.4e} (1D {c1:.4e}, {:+.2}%), 2D DG {b:.4e} (1D DG {d1:.4e}, {:+.2}%; SP {s1:.4e}, {:+.2}%)", 100.0 * (a / c1 - 1.0), 100.0 * (b / d1 - 1.0), 100.0 * (b / s1 - 1.0));
    }
    let _ = vt;
    // 2. the 150 nm MOSFET with and without quantum confinement
    println!("  2. nMOS 150 nm (t_ox 4 nm, N_a 5e17): threshold and on-current with the quantum potential");
    for dg in [false, true] {
        let spec = MosfetSpec::scaled(150e-9, false);
        let mut dev = spec.build(si.clone());
        dev.velocity_saturation = true;
        dev.lombardi = true;
        dev.screening = true;
        if dg {
            dev.dg_mass = Some(m_dg);
        }
        dev.solve_quantum().expect("equilibrium");
        let (d, g) = (dev.contact("drain"), dev.contact("gate"));
        let ramp = |dev: &mut Device, c: usize, v: f64| {
            let v0 = dev.contacts[c].v;
            let steps = ((v - v0).abs() / 0.1).ceil().max(1.0) as usize;
            for s in 1..=steps {
                dev.contacts[c].v = v0 + (v - v0) * s as f64 / steps as f64;
                dev.solve_quantum().expect("bias");
            }
        };
        ramp(&mut dev, d, 0.05);
        // threshold by constant current 1e-7 A·W/L
        let icrit = 1e-7 * spec.width / spec.l_gate;
        let mut prev = (0.0, tcad_id(&dev));
        let mut vth = f64::NAN;
        let mut v = 0.0;
        while v < 1.0 {
            v += 0.05;
            ramp(&mut dev, g, v);
            let i = tcad_id(&dev);
            if i >= icrit {
                let t = (icrit / prev.1).ln() / (i / prev.1).ln();
                vth = prev.0 + t * (v - prev.0);
                break;
            }
            prev = (v, i);
        }
        ramp(&mut dev, g, 1.2);
        ramp(&mut dev, d, 1.2);
        let ion = tcad_id(&dev);
        let lmin = dev.lambda.iter().cloned().fold(0.0, f64::min);
        println!("     {}: V_t {vth:.4} V, I_on {:.1} µA/µm (most negative quantum potential {:.1} V_t)", if dg { "density gradient" } else { "classical       " }, ion * 1e6, lmin);
    }
    let st = MosStack { na: 5e23, tox: 4e-9, t: 300.0, depth: 150e-9, lq: 15e-9, m_dg };
    let q_th = Q * 1e16;
    let vt_of = |mode: Electrons| -> f64 {
        let (mut lo, mut hi) = (-0.5, 2.0);
        for _ in 0..40 {
            let mid = 0.5 * (lo + hi);
            if st.solve(mid, mode).q_inv < q_th {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    };
    let (v_cl, v_sp, v_dg) = (vt_of(Electrons::Classical), vt_of(Electrons::SchrodingerPoisson), vt_of(Electrons::DensityGradient));
    println!("     1D references for this stack (charge criterion 1e12 cm⁻²): threshold shift Schrödinger–Poisson − classical {:+.1} mV, density gradient − classical {:+.1} mV", 1e3 * (v_sp - v_cl), 1e3 * (v_dg - v_cl));
}

fn tunnel() {
    println!("gate tunnelling, Tsu–Esaki with WKB through the oxide barrier (Φ_B 3.15 eV, m_ox 0.5 m₀)");
    println!("  1. limit check: metal emitter (E_F 5 eV, 10 K), 10 nm oxide, against the analytic Fowler–Nordheim law");
    for f_mv in [6.0, 8.0, 10.0, 12.0, 14.0] {
        let f = f_mv * 1e8;
        let j = tsu_esaki(10e-9, f, 5.0 + 3.15, 5.0, 5.0 - f * 10e-9, 10.0, 1.0, 0.5);
        let jfn = fowler_nordheim(3.15, f, 0.5);
        println!("     F {f_mv:4.1} MV/cm: J {:.4e} A/cm² (Fowler–Nordheim {:.4e}, ratio {:.4})", j * 1e-4, jfn * 1e-4, j / jfn);
    }
    println!("  2. MOS inversion layer (n⁺ poly / SiO₂ / p 1e18 cm⁻³, V_g 1.0 V): direct-tunnelling gate current vs oxide thickness");
    let si = Silicon::new(300.0);
    let vt = si.vt();
    let chi_c = (nc(300.0) / si.ni).ln();
    let mut prev: Option<(f64, f64)> = None;
    for tox_nm in [3.0, 2.5, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0] {
        let st = MosStack { na: 1e24, tox: tox_nm * 1e-9, t: 300.0, depth: 150e-9, lq: 15e-9, m_dg: 0.3 };
        let vg = 1.0;
        let sol = st.solve(vg, Electrons::Classical);
        let psi_s = sol.psi[0];
        let psi_g = (vg + 0.5 * band_gap(300.0)) / vt;
        let field = (psi_g - psi_s) * vt / st.tox;
        let eta = psi_s - chi_c;
        let jg = tsu_esaki(st.tox, field, 3.15, eta * vt, -vg + vt * (psi_s - chi_c), 300.0, 0.5, 0.5);
        let dec = prev.map(|(t0, j0)| (jg / j0).log10() / ((t0 - tox_nm) / 0.2));
        println!("     t_ox {tox_nm:3.1} nm: oxide field {:.2} MV/cm, J_g {:.3e} A/cm²{}", field * 1e-8, jg * 1e-4, dec.map_or(String::new(), |d| format!(" ({d:.2} decades per 0.2 nm)")));
        prev = Some((tox_nm, jg));
    }
}

/// Avalanche breakdown of a one-sided abrupt p⁺n junction.
fn breakdown() {
    let si = Silicon::new(300.0);
    let (na, nd) = (1e25, 1e22); // 1e19 / 1e16 cm⁻³
    let (lp, ln) = (1e-6, 7e-6);
    let x = graded(&[(0.0, 50e-9), (lp - 0.2e-6, 10e-9), (lp, 2e-9), (lp + 0.3e-6, 10e-9), (lp + 4e-6, 20e-9), (lp + ln, 100e-9)]);
    let mesh = Mesh::line(x, |_| Material::Silicon);
    let last = mesh.nx() - 1;
    let mk = |ii: bool| {
        let contacts = vec![
            Contact { name: "anode".into(), nodes: line_nodes(&mesh, 0), kind: ContactKind::Ohmic, v: 0.0 },
            Contact { name: "cathode".into(), nodes: line_nodes(&mesh, last), kind: ContactKind::Ohmic, v: 0.0 },
        ];
        let mut d = Device::new(mesh.clone(), si.clone(), move |x, _| if x < lp { (0.0, na) } else { (nd, 0.0) }, contacts, 1.0);
        d.impact_ionization = ii;
        d.solve().expect("equilibrium");
        d
    };
    let (mut d0, mut d1) = (mk(false), mk(true));
    let an = d0.contact("anode");
    // ionization integral on a field profile: breakdown when
    // ∫ α_p exp(−∫_x^W (α_p − α_n) dx') dx = 1 (holes travel toward x = 0)
    // (`exclude`: drop each carrier's α where it is majority, > 10% of the
    // doping, as the device model does — measures that fold's effect)
    let ion_integral = |d: &Device, exclude: bool| -> f64 {
        let f = d.field_magnitude();
        let xs: Vec<f64> = (0..mesh.nx()).map(|i| mesh.x[i]).collect();
        let k = |i: usize| mesh.node(i, 0);
        let an_ = |i: usize| if exclude && d.n(k(i)) > 0.1 * (d.nd[k(i)] + d.na[k(i)]) { 0.0 } else { alpha_n(f[k(i)]) };
        let ap_ = |i: usize| if exclude && d.p(k(i)) > 0.1 * (d.nd[k(i)] + d.na[k(i)]) { 0.0 } else { alpha_p(f[k(i)]) };
        let n = xs.len();
        // inner integral from the right: S(x) = ∫_x^L (α_p − α_n)
        let mut s = vec![0.0; n];
        for i in (0..n - 1).rev() {
            let g = |q: usize| ap_(q) - an_(q);
            s[i] = s[i + 1] + 0.5 * (g(i) + g(i + 1)) * (xs[i + 1] - xs[i]);
        }
        (0..n - 1).map(|i| 0.5 * (ap_(i) * (-s[i]).exp() + ap_(i + 1) * (-s[i + 1]).exp()) * (xs[i + 1] - xs[i])).sum()
    };
    let vt = si.vt();
    let vbi = vt * ((0.5 * na / si.ni).asinh() + (0.5 * nd / si.ni).asinh());
    let bv_sze = 60.0 * (band_gap(300.0) / 1.1f64).powf(1.5) * (nd * 1e-6 / 1e16f64).powf(-0.75);
    println!("p⁺n junction 1e19/1e16 cm⁻³ ({} nodes): reverse bias, impact ionization (van Overstraeten–de Man) on vs off", mesh.nx());
    let mut prev: Option<(f64, f64, f64, f64)> = None;
    let mut bv_dd = f64::NAN;
    let mut bv_int = f64::NAN;
    let mut bv_int_x = f64::NAN;
    let mut v = 0.0;
    while v > -80.0 {
        v -= if v > -40.0 { 5.0 } else { 1.0 };
        d0.ramp(an, v, 1.0).expect("no-II bias");
        let i0 = (d0.current(an).0 + d0.current(an).1).abs();
        let integ = ion_integral(&d0, false);
        let integ_x = ion_integral(&d0, true);
        let r = d1.ramp(an, v, 0.5);
        let (m, ii_nodes) = if r.is_some() {
            let i1 = (d1.current(an).0 + d1.current(an).1).abs();
            (i1 / i0, d1.ii_nodes.get())
        } else {
            (f64::INFINITY, 0)
        };
        println!("  V {v:6.1}: |J| no-II {:.3e} A/cm², multiplication M {:9.4}, ionization integral {integ:.5} (majority tails excluded {integ_x:.5}; nodes ionizing {ii_nodes})", i0 * 1e-4, m);
        if let Some((v0, m0, i0_, ix0)) = prev {
            if bv_int.is_nan() && i0_ < 1.0 && integ >= 1.0 {
                bv_int = v0 + (1.0 - i0_) / (integ - i0_) * (v - v0);
            }
            if bv_int_x.is_nan() && ix0 < 1.0 && integ_x >= 1.0 {
                bv_int_x = v0 + (1.0 - ix0) / (integ_x - ix0) * (v - v0);
            }
            if m.is_finite() && m > 2.0 && m0 > 1.5 {
                // 1/M → 0 linearly at breakdown: extrapolate from the last
                // two converged points
                let (a, b) = (1.0 / m0, 1.0 / m);
                bv_dd = v0 + a / (a - b) * (v - v0);
            }
        }
        if !m.is_finite() || m > 200.0 {
            break;
        }
        prev = Some((v, m, integ, integ_x));
    }
    println!("  breakdown voltage: drift–diffusion (1/M → 0) {:.2} V; ionization integral = 1 on the same field {:.2} V ({:.2} V with the majority-tail fold); Sze–Gibbons empirical {:.1} V (V_bi {vbi:.2} V)", -bv_dd, -bv_int, -bv_int_x, bv_sze);
}

fn moscap() {
    let si = Silicon::new(300.0);
    let vt = si.vt();
    let (na, tox, lsi) = (1e23, 10e-9, 1e-6);
    let offset = 0.5 * band_gap(300.0); // n⁺ poly gate
    let mos = Mos { na, ni: si.ni, vt, eps_s: EPS_SI, cox: EPS_OX / tox, vfb: -vt * (0.5 * na / si.ni).asinh() - offset };
    println!("MOS capacitor: n⁺ poly / 10 nm SiO₂ / p-Si 1e17 cm⁻³; V_FB {:.4} V (exact Poisson–Boltzmann reference)", mos.vfb);
    let vgs: Vec<f64> = (0..=60).map(|i| -2.0 + 0.075 * i as f64).collect();
    for refine in [1.0, 2.0, 4.0] {
        let h = 0.1e-9 / refine;
        let x = graded(&[(0.0, 1e-9 / refine), (tox, h), (tox + 20e-9, 2e-9 / refine), (tox + lsi, 40e-9 / refine)]);
        let mesh = Mesh::line(x.clone(), move |x| if x < tox { Material::Oxide } else { Material::Silicon });
        let last = mesh.nx() - 1;
        let contacts = vec![
            Contact { name: "gate".into(), nodes: line_nodes(&mesh, 0), kind: ContactKind::Gate { offset }, v: 0.0 },
            Contact { name: "body".into(), nodes: line_nodes(&mesh, last), kind: ContactKind::Ohmic, v: 0.0 },
        ];
        let mut dev = Device::new(mesh.clone(), si.clone(), |_, _| (0.0, na), contacts, 1.0);
        dev.solve().expect("equilibrium");
        let g = dev.contact("gate");
        let mut worst: f64 = 0.0;
        let mut rows = Vec::new();
        for &vg in &vgs {
            dev.ramp(g, vg, 0.1).expect("gate sweep");
            let q = dev.charge(g);
            let qe = mos.qg(vg);
            worst = worst.max((q - qe).abs());
            rows.push((vg, q, qe));
        }
        let qmax = rows.iter().map(|r| r.2.abs()).fold(0.0, f64::max);
        println!("  {} nodes (interface spacing {:.3} nm): max |Q_g − exact| = {:.3e} of max |Q_g|", mesh.nx(), h * 1e9, worst / qmax);
        if refine == 4.0 {
            println!("    low-frequency C/C_ox (dQ_g/dV_g) vs exact:");
            for w in rows.windows(3).step_by(4) {
                let c = (w[2].1 - w[0].1) / (w[2].0 - w[0].0) / mos.cox;
                let ce = (w[2].2 - w[0].2) / (w[2].0 - w[0].0) / mos.cox;
                println!("      V_g {:+.3} V: {:.5} (exact {:.5})", w[1].0, c, ce);
            }
        }
    }
    let phi_f = vt * (0.5 * na / si.ni).asinh();
    let vth = mos.vfb + 2.0 * phi_f + (4.0 * EPS_SI * Q * na * phi_f).sqrt() / mos.cox;
    println!("  classical threshold V_FB + 2φ_F + √(4ε_s q N_a φ_F)/C_ox = {vth:.4} V");
}

/// Manufactured solution on the unit-scaled domain: returns the maximum
/// nodal errors in (ψ, φ_n, φ_p) (units of V_t).
fn mms(mesh: Mesh, lx: f64, ly: f64, heavy: bool) -> [f64; 3] {
    let si = Silicon::new(300.0);
    let vt = si.vt();
    let ni = si.ni;
    let (mun, mup) = (0.1, 0.04);
    let two_d = mesh.ny() > 2;
    let fy = move |y: f64| if two_d { (std::f64::consts::PI * y / ly).cos() } else { 1.0 };
    // heavy: electrons pushed deep into degeneracy (ψ − φ_n up to ~24 V_t,
    // η > 0 against χ_c = ln(N_c/n_i) ≈ 21.8), Fermi–Dirac + band-gap
    // narrowing on 3e20 cm⁻³ doping
    // heavy: smooth deep potential ψ = 23.5 tanh(·) (degenerate on both
    // sides, η ≈ +2), compensated doping with N_d + N_a = 3e26 m⁻³ held
    // constant (band-gap narrowing smooth) and N_d − N_a = n − p of the
    // manufactured state (exactly neutral, a well-posed problem)
    let psi = move |x: f64, y: f64| {
        if heavy {
            23.5 * ((x - 0.5 * lx) / (0.12 * lx)).tanh() * (0.9 + 0.1 * fy(y)) + 0.3 * (std::f64::consts::PI * x / lx).sin()
        } else {
            8.0 * (std::f64::consts::PI * x / lx).sin() * (0.8 + 0.2 * fy(y)) - 2.0
        }
    };
    let phin = move |x: f64, y: f64| 1.5 * x / lx + 0.3 * (2.0 * std::f64::consts::PI * x / lx).sin() * fy(y);
    let phip = move |x: f64, y: f64| -1.0 + 0.5 * (std::f64::consts::PI * x / lx).cos() + 0.2 * fy(y);
    let (use_bgn, use_fd) = (heavy, heavy);
    let ntot_heavy = 3e26;
    let (chi_c, chi_v) = ((nc(300.0) / ni).ln(), (nv(300.0) / ni).ln());
    let dgap = move |_x: f64, _y: f64| if use_bgn { bgn_slotboom(ntot_heavy) / vt } else { 0.0 };
    let th_n = move |x: f64, y: f64, u: f64| 0.5 * dgap(x, y) + if use_fd { degeneracy(u + 0.5 * dgap(x, y) - chi_c).0 } else { 0.0 };
    let th_p = move |x: f64, y: f64, w: f64| 0.5 * dgap(x, y) + if use_fd { degeneracy(w + 0.5 * dgap(x, y) - chi_v).0 } else { 0.0 };
    let dop = move |x: f64, y: f64| -> (f64, f64) {
        if heavy {
            let (ps, fnn, fpp) = (psi(x, y), phin(x, y), phip(x, y));
            let c = ni * ((ps - fnn + th_n(x, y, ps - fnn)).exp() - (fpp - ps + th_p(x, y, fpp - ps)).exp());
            (0.5 * (ntot_heavy + c), 0.5 * (ntot_heavy - c))
        } else {
            let c = 3e20 * ((x - 0.5 * lx) / (0.15 * lx)).tanh();
            (c.max(0.0), (-c).max(0.0))
        }
    };
    // boundary nodes: fixed to the manufactured solution
    let (nx, ny) = (mesh.nx(), mesh.ny());
    let mut nodes = Vec::new();
    let mut vals = Vec::new();
    for i in 0..nx {
        for j in 0..ny {
            let boundary = i == 0 || i == nx - 1 || (two_d && (j == 0 || j == ny - 1));
            if boundary {
                let (x, y) = (mesh.x[i], mesh.y[j]);
                nodes.push(mesh.node(i, j));
                vals.push([psi(x, y), phin(x, y), phip(x, y)]);
            }
        }
    }
    let contacts = vec![Contact { name: "boundary".into(), nodes, kind: ContactKind::Fixed(vals), v: 0.0 }];
    let mut dev = Device::new(mesh.clone(), si.clone(), dop, contacts, 1.0);
    dev.mu_const = Some((mun, mup));
    dev.update_mobility();
    if heavy {
        if use_fd {
            dev.statistics = Statistics::FermiDirac;
        }
        dev.bgn = use_bgn;
        dev.update_statistics();
    }
    // continuous operators by nested central differences
    let kappa = Q * ni / vt;
    let d = 2e-5 * lx.min(if two_d { ly } else { lx });
    let nf = move |x: f64, y: f64| (psi(x, y) - phin(x, y) + th_n(x, y, psi(x, y) - phin(x, y))).exp();
    let pf = move |x: f64, y: f64| (phip(x, y) - psi(x, y) + th_p(x, y, phip(x, y) - psi(x, y))).exp();
    let grad = |f: &dyn Fn(f64, f64) -> f64, x: f64, y: f64| -> (f64, f64) { ((f(x + d, y) - f(x - d, y)) / (2.0 * d), if two_d { (f(x, y + d) - f(x, y - d)) / (2.0 * d) } else { 0.0 }) };
    let div = |fx: &dyn Fn(f64, f64) -> (f64, f64), x: f64, y: f64| -> f64 {
        let mut s = (fx(x + d, y).0 - fx(x - d, y).0) / (2.0 * d);
        if two_d {
            s += (fx(x, y + d).1 - fx(x, y - d).1) / (2.0 * d);
        }
        s
    };
    // J = −μ n ∇φ (exact for any statistics)
    let jn = |x: f64, y: f64| {
        let g = grad(&phin, x, y);
        let n = nf(x, y);
        (-mun * n * g.0, -mun * n * g.1)
    };
    let jp = |x: f64, y: f64| {
        let g = grad(&phip, x, y);
        let p = pf(x, y);
        (-mup * p * g.0, -mup * p * g.1)
    };
    let dpsi = |x: f64, y: f64| {
        let g = grad(&psi, x, y);
        (EPS_SI * g.0, EPS_SI * g.1)
    };
    let mut src = vec![[0.0; 3]; mesh.nodes()];
    for i in 1..nx - 1 {
        for j in 0..ny {
            if two_d && (j == 0 || j == ny - 1) {
                continue;
            }
            let k = mesh.node(i, j);
            let (x, y) = (mesh.x[i], mesh.y[j]);
            let (dd, da) = dop(x, y);
            let (n, p) = (nf(x, y), pf(x, y));
            let ntot = dd + da;
            let (tn, tp) = (si.tau_n(ntot), si.tau_p(ntot));
            let e = (th_n(x, y, psi(x, y) - phin(x, y)) + th_p(x, y, phip(x, y) - psi(x, y))).exp();
            let r = (n * p - e) / (tp * (n + e.sqrt()) + tn * (p + e.sqrt())) + ni * ni * (si.c_n * n + si.c_p * p) * (n * p - e);
            let f_psi = -div(&dpsi, x, y) - kappa * (p - n + (dd - da) / ni);
            let f_n = div(&jn, x, y) - r / vt;
            let f_p = div(&jp, x, y) + r / vt;
            let v = dev.semiconductor_volume(k);
            src[k] = [-v * f_psi, -v * f_n, -v * f_p];
        }
    }
    dev.source = Some(src);
    // start from a perturbed manufactured state
    for k in 0..mesh.nodes() {
        let (i, j) = mesh.ij(k);
        let (x, y) = (mesh.x[i], mesh.y[j]);
        dev.psi[k] = psi(x, y) + 0.3;
        dev.phin[k] = phin(x, y) - 0.2;
        dev.phip[k] = phip(x, y) + 0.2;
    }
    let its = dev.solve().expect("MMS Newton failed");
    if use_fd {
        assert!(dev.fd_nodes.get() > 0);
    }
    let mut err = [0.0f64; 3];
    for k in 0..mesh.nodes() {
        let (i, j) = mesh.ij(k);
        let (x, y) = (mesh.x[i], mesh.y[j]);
        err[0] = err[0].max((dev.psi[k] - psi(x, y)).abs());
        err[1] = err[1].max((dev.phin[k] - phin(x, y)).abs());
        err[2] = err[2].max((dev.phip[k] - phip(x, y)).abs());
    }
    let _ = its;
    err
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let case = a.get(1).cloned().unwrap_or("mms".into());
    match case.as_str() {
        "mms" => {
            let lx = 1e-6;
            println!("manufactured solution, 1D, uniform mesh (L = 1 µm, densities up to 3e20 m⁻³ doping, ψ amplitude 8 V_t):");
            let mut prev: Option<[f64; 3]> = None;
            for n in [20, 40, 80, 160, 320] {
                let e = mms(Mesh::line(uniform(0.0, lx, n), |_| Material::Silicon), lx, 1.0, false);
                let ord = prev.map(|p| [(p[0] / e[0]).log2(), (p[1] / e[1]).log2(), (p[2] / e[2]).log2()]);
                println!("  {n:4} cells: max error ψ {:.3e}, φn {:.3e}, φp {:.3e} V_t{}", e[0], e[1], e[2], ord.map_or(String::new(), |o| format!("   order {:.2} {:.2} {:.2}", o[0], o[1], o[2])));
                prev = Some(e);
            }
            println!("manufactured solution, 1D, degenerate (Fermi–Dirac + band-gap narrowing, compensated 3e20 cm⁻³ doping, ψ ±23.5 V_t, η ≈ +2):");
            let mut prev: Option<[f64; 3]> = None;
            for n in [20, 40, 80, 160, 320] {
                let e = mms(Mesh::line(uniform(0.0, lx, n), |_| Material::Silicon), lx, 1.0, true);
                let ord = prev.map(|p| [(p[0] / e[0]).log2(), (p[1] / e[1]).log2(), (p[2] / e[2]).log2()]);
                println!("  {n:4} cells: max error ψ {:.3e}, φn {:.3e}, φp {:.3e} V_t{}", e[0], e[1], e[2], ord.map_or(String::new(), |o| format!("   order {:.2} {:.2} {:.2}", o[0], o[1], o[2])));
                prev = Some(e);
            }
            println!("manufactured solution, 2D, graded mesh (1 × 0.5 µm, spacing ratio 4 across the domain):");
            let ly = 0.5e-6;
            let mut prev: Option<[f64; 3]> = None;
            for s in [1.0, 2.0, 4.0, 8.0] {
                let h0 = 0.05e-6 / s;
                let x = graded(&[(0.0, h0 / 2.0), (lx, 2.0 * h0)]);
                let y = graded(&[(0.0, h0), (ly, h0 / 2.0)]);
                let (nx, ny) = (x.len(), y.len());
                let e = mms(Mesh::new(x, y, |_, _| Material::Silicon), lx, ly, false);
                let ord = prev.map(|p| [(p[0] / e[0]).log2(), (p[1] / e[1]).log2(), (p[2] / e[2]).log2()]);
                println!("  {nx:3}×{ny:3} nodes: max error ψ {:.3e}, φn {:.3e}, φp {:.3e} V_t{}", e[0], e[1], e[2], ord.map_or(String::new(), |o| format!("   order {:.2} {:.2} {:.2}", o[0], o[1], o[2])));
                prev = Some(e);
            }
        }
        "silicon" => {
            println!("silicon band structure (Thurmond gap, Green 1990 densities of states) vs measured n_i (Misiakos & Tsamakis 1993):");
            for t in [200.0, 250.0, 300.0, 350.0, 400.0] {
                println!("  T {t:5.0} K: E_g {:.4} eV, n_i from bands {:.3e} cm⁻³, measured {:.3e} ({:+.1}%)", band_gap(t), ni_from_bands(t) * 1e-6, ni_measured(t) * 1e-6, 100.0 * (ni_from_bands(t) / ni_measured(t) - 1.0));
            }
            let si = Silicon::new(300.0);
            println!("  low-field mobility at 300 K (Masetti):");
            for n in [1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20] {
                let mn = si.mob_n.mobility(n * 1e6, 300.0) * 1e4;
                let mp = si.mob_p.mobility(n * 1e6, 300.0) * 1e4;
                let rho_n = 1.0 / (Q * n * mn);
                let rho_p = 1.0 / (Q * n * mp);
                println!("    N {n:7.0e} cm⁻³: μ_n {mn:7.1}, μ_p {mp:6.1} cm²/(V s); resistivity n-type {rho_n:9.4}, p-type {rho_p:9.4} Ω·cm");
            }
        }
        "diode" => diode(),
        "moscap" => moscap(),
        "mosfet" => mosfet(a.get(2).map_or(10.0, |v| v.parse().unwrap())),
        "short" => short_channel(),
        "cmos" => cmos(),
        "breakdown" => breakdown(),
        "mobility" => mobility(),
        "transient" => transient_moscap(),
        "switching" => switching(),
        "thermal" => thermal(),
        "quantum" => quantum(),
        "tunnel" => tunnel(),
        "quantum2d" => quantum2d(),
        _ => panic!("unknown case"),
    }
}
