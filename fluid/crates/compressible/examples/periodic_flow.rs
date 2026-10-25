//! Every element's phase transitions inside the flow solver.
//!
//!   periodic_flow [N]
//!
//! 1. Equilibrium round trip, solid–liquid at T_m and liquid–vapour at T_b
//!    (1 atm): build a 30%-transformed mixture state (v, e) and ask the
//!    relaxation solver for (p, T, Y); it must return (1 atm, T, 0.3).
//! 2. Expansion (cavitation) tube with the element's own liquid and vapour
//!    in the compressible five-equation solver with phase-change
//!    relaxation: the liquid, between its melting and boiling points with
//!    a 1% vapour trace, is pulled apart; the centre must settle on the
//!    model's saturation pressure p_sat(T) with no loss of admissibility.
use compressible::elements::{all, fit, Phase, P_ATM};
use compressible::multiphase::FiveEq;
use compressible::phase_change::{e_pt, p_sat, t_sat_from, v_pt, PhaseChange};
use compressible::{Bc, Mesh, Solver};

fn main() {
    let n: usize = std::env::args().nth(1).map(|v| v.parse().unwrap()).unwrap_or(200);
    let (mut ok_eq, mut n_eq, mut ok_tube, mut n_tube) = (0, 0, 0, 0);
    let mut worst_tube = (0.0f64, String::new());
    let mut fails = Vec::new();
    println!("  Z      melt (T, Y) error      boil (T, Y) error      tube: centre p/p_sat(T_l)  α_v   steps");
    for el in all() {
        let Ok(m) = fit(el) else { continue };
        if el.sublimes {
            continue;
        }
        let (tm, tb) = (el.t_melt.unwrap(), el.t_boil.unwrap());
        // 1. equilibrium round trips
        let mut eq_txt = String::new();
        for (a, b, t) in [(Phase::Solid, Phase::Liquid, tm), (Phase::Liquid, Phase::Vapor, tb)] {
            let (ea, eb) = (*m.eos(a), *m.eos(b));
            let sys = FiveEq { eos: [ea, eb], kapila: true };
            let pc = PhaseChange { liquid: 0, tension: 0.0, present: 1e-9, t_guess: t };
            let y = 0.3;
            let v = (1.0 - y) * v_pt(&ea, P_ATM, t) + y * v_pt(&eb, P_ATM, t);
            let e = (1.0 - y) * e_pt(&ea, P_ATM, t) + y * e_pt(&eb, P_ATM, t);
            n_eq += 1;
            match pc.equilibrium(&sys, v, e) {
                Some(q) => {
                    let (dt, dy, dp) = ((q.t - t) / t, q.y - y, (q.p - P_ATM) / P_ATM);
                    let good = dt.abs() < 1e-6 && dy.abs() < 1e-6 && dp.abs() < 1e-5;
                    ok_eq += good as usize;
                    if !good {
                        fails.push(format!("{} {:?}-{:?} equilibrium off (dT/T {dt:.1e}, dY {dy:.1e}, dp/p {dp:.1e})", el.symbol, a, b));
                    }
                    eq_txt += &format!("{:9.1e} {:9.1e}   ", dt, dy);
                }
                None => {
                    fails.push(format!("{} {:?}-{:?} equilibrium not found", el.symbol, a, b));
                    eq_txt += "   not found          ";
                }
            }
        }
        // 2. cavitation tube in the flow solver
        let (l, v) = (m.liquid, m.vapor);
        let t0 = tm + 0.6 * (tb - tm);
        let p0 = P_ATM;
        let (rl, rv) = (1.0 / v_pt(&l, p0, t0), 1.0 / v_pt(&v, p0, t0));
        let sys = FiveEq { eos: [l, v], kapila: true };
        let av0 = 1e-2;
        let mut w0 = [0.0; 7];
        w0[0] = (1.0 - av0) * rl;
        w0[1] = av0 * rv;
        w0[5] = p0;
        w0[6] = 1.0 - av0;
        let c_mix = compressible::System::sound_speed(&sys, &w0);
        let u0 = 0.05 * c_mix;
        let mesh = Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]);
        let mut s = Solver::new(sys, mesh, |x| {
            let mut w = w0;
            w[2] = if x[0] < 0.5 { -u0 } else { u0 };
            w
        });
        s.bc[0] = [Bc::Outflow; 2];
        s.cfl = 0.4;
        let pc = PhaseChange { liquid: 0, tension: 0.0, present: 1e-6, t_guess: t0 };
        s.relax = Some(Box::new(move |sys, u, _dt| pc.relax(sys, u)));
        let t_end = 0.3 / c_mix;
        let mut steps = 0;
        let mut healthy = true;
        while s.t < t_end && steps < 20_000 {
            let dt = s.dt_cfl().min(t_end - s.t);
            if !(dt.is_finite() && dt > 0.0) {
                healthy = false;
                break;
            }
            s.step(dt);
            steps += 1;
        }
        n_tube += 1;
        let w = s.prim(n / 2, 0, 0);
        let rho_l = w[0] / w[6].max(1e-12);
        let t_l = l.temperature(rho_l, w[5]);
        let ps = p_sat(&l, &v, t_l);
        let _ = t_sat_from;
        let ratio = w[5] / ps;
        let good = healthy && w.iter().all(|x| x.is_finite()) && (ratio - 1.0).abs() < 1e-3;
        if good {
            ok_tube += 1;
        } else {
            fails.push(format!("{} tube: p/p_sat {ratio:.4}, healthy {healthy}", el.symbol));
        }
        if (ratio - 1.0).abs() > worst_tube.0 {
            worst_tube = ((ratio - 1.0).abs(), el.symbol.clone());
        }
        println!("{:>3} {:<3}  {eq_txt} {:>10.6}            {:.3}  {steps}", el.z, el.symbol, ratio, 1.0 - w[6]);
    }
    println!();
    println!("equilibrium round trips: {ok_eq}/{n_eq} exact (|dT/T|, |dY| < 1e-6)");
    println!("cavitation tubes: {ok_tube}/{n_tube} settle on p_sat(T_l) within 0.1% (worst |p/p_sat − 1| = {:.2e}, {})", worst_tube.0, worst_tube.1);
    for f in &fails {
        println!("  fail: {f}");
    }
}
