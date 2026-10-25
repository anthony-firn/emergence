//! Solid + liquid + vapour of one substance in one flow.
//!
//!   three_phase [symbol] [N]
//!
//! 1. Maximum-entropy equilibrium round trips: states built at known
//!    (p, T, mass fractions) — single phases, melting, boiling, sublimation
//!    and triple-point mixtures — must be recovered from (v, e) alone.
//! 2. Triple-point tube in the three-phase Kapila flow solver: a uniform
//!    solid–liquid–vapour mixture at (T_t, p_t) is pulled apart (left/right
//!    at ∓u) or pushed together. Three coexisting phases pin p and T to the
//!    triple point: expansion must evaporate liquid and freeze liquid (the
//!    evaporative-cooling route to nitrogen slush), compression the
//!    reverse — with mass and energy conserved.
use compressible::elements::{by_symbol, fit, Phase, P_ATM};
use compressible::phase_change::{e_pt, v_pt};
use compressible::threephase::{Equilibrium3, Kapila3, N3};
use compressible::{Bc, Mesh, Solver, System};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let sym = a.get(1).cloned().unwrap_or("N".into());
    let n: usize = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(200);
    let el = by_symbol(&sym).expect("element");
    let model = fit(el).expect("fit");
    let eq = Equilibrium3::new(model.clone(), el.t_melt.unwrap()).expect("triple point");
    let (tt, pt) = eq.triple;
    println!("{} ({}): model triple point T {tt:.3} K, p {pt:.4e} Pa (data T {:?}, p {:?})", el.name, el.symbol, el.t_triple, el.p_triple);
    let eos = [model.solid, model.liquid, model.vapor];
    // 1. round trips
    let mix = |p: f64, t: f64, y: [f64; 3]| -> (f64, f64) {
        let v: f64 = (0..3).map(|k| y[k] * v_pt(&eos[k], p, t)).sum();
        let e: f64 = (0..3).map(|k| y[k] * e_pt(&eos[k], p, t)).sum();
        (v, e)
    };
    let tm = model.coexistence_t(Phase::Solid, Phase::Liquid, P_ATM, 0.3 * tt, 3.0 * tt).unwrap();
    let tb = model.coexistence_t(Phase::Liquid, Phase::Vapor, P_ATM, 0.3 * tt, 5.0 * tt).unwrap();
    let p_sub = 0.3 * pt;
    let t_sub = model.coexistence_t(Phase::Solid, Phase::Vapor, p_sub, 0.3 * tt, 2.0 * tt).unwrap();
    let cases: Vec<(&str, f64, f64, [f64; 3])> = vec![
        ("solid", P_ATM, 0.8 * tm, [1.0, 0.0, 0.0]),
        ("liquid", P_ATM, 0.5 * (tm + tb), [0.0, 1.0, 0.0]),
        ("vapour", P_ATM, 1.3 * tb, [0.0, 0.0, 1.0]),
        ("melting at 1 atm", P_ATM, tm, [0.4, 0.6, 0.0]),
        ("boiling at 1 atm", P_ATM, tb, [0.0, 0.7, 0.3]),
        ("sublimation at 0.3 p_t", p_sub, t_sub, [0.5, 0.0, 0.5]),
        ("triple point", pt, tt, [0.3, 0.5, 0.2]),
        ("triple point, vapour-rich", pt, tt, [0.1, 0.2, 0.7]),
    ];
    println!("  equilibrium round trips (state recovered from v, e):");
    for (name, p, t, y) in cases {
        let (v, e) = mix(p, t, y);
        match eq.equilibrate(v, e) {
            Some(s) => println!(
                "    {name:28} → {:22} dT/T {:+.1e}  dp/p {:+.1e}  dy {:.1e}",
                Equilibrium3::label(&s),
                (s.t - t) / t,
                (s.p - p) / p,
                (0..3).map(|k| (s.y[k] - y[k]).abs()).fold(0.0, f64::max)
            ),
            None => println!("    {name:28} → not found"),
        }
    }
    // 2. triple-point tube in the flow solver
    let sys = Kapila3 { eos };
    let y0 = [0.3, 0.6, 0.1];
    let rk: [f64; 3] = std::array::from_fn(|k| 1.0 / v_pt(&eos[k], pt, tt));
    let mut w0 = [0.0; N3];
    let mut al = [0.0; 3];
    let rho = 1.0 / (0..3).map(|k| y0[k] / rk[k]).sum::<f64>();
    for k in 0..3 {
        w0[k] = y0[k] * rho;
        al[k] = w0[k] / rk[k];
    }
    w0[6] = pt;
    w0[7] = al[0];
    w0[8] = al[1];
    let c_mix = sys.sound_speed(&w0);
    println!("  mixture: ρ {rho:.3} kg/m³, volume fractions {al:.3?}, Wood sound speed {c_mix:.2} m/s");
    for (label, sign) in [("expansion", 1.0), ("compression", -1.0)] {
        let u0 = 0.05 * c_mix;
        let mesh = Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]);
        let mut s = Solver::new(sys, mesh, |x| {
            let mut w = w0;
            w[3] = sign * if x[0] < 0.5 { -u0 } else { u0 };
            w
        });
        s.bc[0] = [Bc::Outflow; 2];
        s.cfl = 0.4;
        let eqc = eq.clone();
        s.relax = Some(Box::new(move |sys, u, _dt| eqc.relax(sys, u)));
        let tot0 = s.totals();
        let t_end = 0.2 / c_mix;
        let mut steps = 0;
        while s.t < t_end && steps < 50_000 {
            let dt = s.dt_cfl().min(t_end - s.t);
            s.step(dt);
            steps += 1;
        }
        let tot1 = s.totals();
        let wc = s.prim(n / 2, 0, 0);
        let rc = Kapila3::rho(&wc);
        let yc: [f64; 3] = std::array::from_fn(|k| wc[k] / rc);
        let ke = 0.5 * rc * wc[3] * wc[3];
        let st = eq.equilibrate(1.0 / rc, (s.sys.to_cons(&wc)[6] - ke) / rc);
        println!(
            "  {label:11}: centre p/p_t {:.6}, T {:.3} K, mass fractions {:.4?} (start {:?}), Δmass {:.1e}, ΔE {:.1e} (outflow BCs), {steps} steps, active cells {}",
            wc[6] / pt,
            st.map(|x| x.t).unwrap_or(f64::NAN),
            yc,
            y0,
            (tot1[0] + tot1[1] + tot1[2]) / (tot0[0] + tot0[1] + tot0[2]) - 1.0,
            tot1[6] / tot0[6] - 1.0,
            s.relaxed_last
        );
    }
}
