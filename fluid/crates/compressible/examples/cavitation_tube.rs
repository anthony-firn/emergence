//! Expansion (cavitation) tube, after Saurel, Petitpas & Abgrall (2008):
//! liquid water at 355 K, 1 bar, with a trace of vapour (α_v = 1e-2),
//! pulled apart at ±u. Without mass transfer the expansion drives the
//! pressure far below saturation (the liquid is put in tension); with the
//! equilibrium relaxation the central region boils, pressure stays at
//! p_sat(T) and vapour forms, and the waves travel at the much lower
//! equilibrium sound speed.
//!
//!   cavitation_tube N [u] [t_end]
use compressible::eos::materials::{WATER_LIQUID_NASG, WATER_VAPOR_NASG};
use compressible::multiphase::FiveEq;
use compressible::phase_change::{p_sat, v_pt, PhaseChange};
type V = [f64; 7];
const MAXV: usize = 7;
use compressible::{Bc, Mesh, Solver};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let n: usize = a.get(1).map(|v| v.parse().unwrap()).unwrap_or(1000);
    let u0: f64 = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(2.0);
    let t_end: f64 = a.get(3).map(|v| v.parse().unwrap()).unwrap_or(3.2e-3);
    let (l, v) = (WATER_LIQUID_NASG, WATER_VAPOR_NASG);
    let (p0, t0, av0) = (1e5, 355.0, 1e-2);
    let (rl, rv) = (1.0 / v_pt(&l, p0, t0), 1.0 / v_pt(&v, p0, t0));
    let sys = FiveEq { eos: [l, v], kapila: true };
    println!("liquid ρ {rl:.2} kg/m³, vapour ρ {rv:.4} kg/m³ at {p0} Pa, {t0} K; p_sat({t0} K) = {:.0} Pa (NASG)", p_sat(&l, &v, t0));
    let init = |x: [f64; 3]| -> V {
        let mut w = [0.0; MAXV];
        let al = 1.0 - av0;
        w[0] = al * rl;
        w[1] = av0 * rv;
        w[2] = if x[0] < 0.5 { -u0 } else { u0 };
        w[5] = p0;
        w[6] = al;
        w
    };
    for with_pc in [false, true] {
        let mesh = Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]);
        let mut s = Solver::new(sys, mesh, init);
        s.bc[0] = [Bc::Outflow; 2];
        s.cfl = 0.4;
        if with_pc {
            let pc = PhaseChange { liquid: 0, tension: 0.0, present: 1e-6, t_guess: 355.0 };
            s.relax = Some(Box::new(move |sys, u, _dt| pc.relax(sys, u)));
        }
        let tot0 = s.totals();
        let mut max_active = 0;
        while s.t < t_end {
            let dt = s.dt_cfl().min(t_end - s.t);
            s.step(dt);
            max_active = max_active.max(s.relaxed_last);
        }
        let tot1 = s.totals();
        let mid = n / 2;
        let w = s.prim(mid, 0, 0);
        let rho_l = w[0] / w[6];
        let t_l = l.temperature(rho_l, w[5]);
        let pmin = (0..n).map(|i| s.prim(i, 0, 0)[5]).fold(f64::INFINITY, f64::min);
        let vapour_width = (0..n).filter(|&i| 1.0 - s.prim(i, 0, 0)[6] > 0.1).count() as f64 / n as f64;
        println!(
            "{}: centre p = {:.1} Pa, p_sat(T_l = {:.2} K) = {:.1} Pa, α_v = {:.4}; min p {:.3e} Pa; region with α_v > 0.1: {:.3} m; \
             mass change {:.1e}, energy change {:.1e} (outflow BCs); active cells (max) {max_active}",
            if with_pc { "with phase change   " } else { "without phase change" },
            w[5],
            t_l,
            p_sat(&l, &v, t_l),
            1.0 - w[6],
            pmin,
            vapour_width,
            (tot1[0] + tot1[1]) / (tot0[0] + tot0[1]) - 1.0,
            tot1[5] / tot0[5] - 1.0,
        );
    }
}
