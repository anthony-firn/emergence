//! Finite-rate cavitation with nucleation, and gas dissolution (`kinetics.rs`).
//!
//!   kinetics_validate cnt          cavitation thresholds: CNT (homogeneous) and Blake (nuclei)
//!   kinetics_validate tube [N]     expansion tube from pure liquid: none / equilibrium / finite rate
//!   kinetics_validate henry [N]    nitrogen bubbles dissolving in degassed water (minutes, all-speed)
use compressible::allspeed::AllSpeedOpts;
use compressible::eos::materials::{WATER_LIQUID_NASG, WATER_SG, WATER_VAPOR_NASG};
use compressible::eos::Nasg;
use compressible::kinetics::{blake_threshold, cnt_rate, Cavitation, Dissolution, Nuclei};
use compressible::multiphase::FiveEq;
use compressible::phase_change::{p_sat, v_pt, PhaseChange};
use compressible::scalars::{Carrier, Scalar};
use compressible::{Bc, Mesh, Solver, System};

const SIGMA_300K: f64 = 0.0717;
const M_WATER: f64 = 0.018015 / 6.02214076e23;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let case = a.get(1).cloned().unwrap_or("cnt".into());
    let n: usize = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(500);
    let (l, v) = (WATER_LIQUID_NASG, WATER_VAPOR_NASG);
    match case.as_str() {
        "cnt" => {
            let t = 300.0;
            let pv = p_sat(&l, &v, t);
            let rho = 1.0 / v_pt(&l, 1e5, t);
            // pressure where J V τ = 1 for V = 1 mm³, τ = 1 s (bisection in p)
            let target = 1.0 / (1e-9 * 1.0);
            let (mut lo, mut hi) = (-1e9f64, pv - 1.0);
            for _ in 0..200 {
                let mid = 0.5 * (lo + hi);
                if cnt_rate(rho, M_WATER, SIGMA_300K, t, mid, pv) > target {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            println!("water {t} K: p_v = {pv:.0} Pa (NASG), σ = {SIGMA_300K} N/m");
            println!("  CNT cavitation pressure (J·1 mm³·1 s = 1): {:.1} MPa (experiment, water in quartz inclusions: about −140 MPa, Azouzi et al. 2013)", 0.5 * (lo + hi) / 1e6);
            for (r, label) in [(1e-6, "1 µm"), (1e-5, "10 µm"), (1e-4, "100 µm")] {
                println!("  Blake threshold, {label} nucleus at 1 bar: p_v − {:.0} Pa", pv - blake_threshold(r, 1e5, pv, SIGMA_300K));
            }
        }
        "tube" => {
            // Saurel–Petitpas–Abgrall expansion tube, but starting from pure
            // liquid (no vapour at all): vapour must nucleate
            let (p0, t0, av0, u0, t_end) = (1e5, 355.0, 0.0, 2.0, 3.2e-3);
            let (rl, rv) = (1.0 / v_pt(&l, p0, t0), 1.0 / v_pt(&v, p0, t0));
            let sys = FiveEq { eos: [l, v], kapila: true };
            let pv = p_sat(&l, &v, t0);
            println!("expansion tube from pure liquid at {t0} K (p_v {pv:.0} Pa), ±{u0} m/s, {n} cells, t = {t_end} s");
            let init = move |x: [f64; 3]| -> [f64; 7] { [(1.0 - av0) * rl, av0 * rv, if x[0] < 0.5 { -u0 } else { u0 }, 0.0, 0.0, p0, 1.0 - av0] };
            let pc = PhaseChange { liquid: 0, tension: 0.0, present: 1e-6, t_guess: t0 };
            let mut cases: Vec<(String, Option<f64>, bool)> = vec![("no phase change".into(), None, false), ("equilibrium (HEM)".into(), None, true)];
            for nd in [1e8, 1e10, 1e12, 1e14] {
                cases.push((format!("finite rate, n = {nd:.0e} m⁻³"), Some(nd), false));
            }
            let only: Option<usize> = std::env::var("CASE").ok().map(|v| v.parse().unwrap());
            for (ci, (name, nd, hem)) in cases.into_iter().enumerate() {
                if only.is_some_and(|o| o != ci) {
                    continue;
                }
                let mut s = Solver::new(sys, Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]), init);
                s.bc[0] = [Bc::Outflow; 2];
                s.cfl = 0.4;
                if hem {
                    s.relax = Some(Box::new(move |sys, u, _dt| pc.relax(sys, u)));
                }
                if let Some(nd) = nd {
                    let cav = Cavitation { pc, sigma: 0.0636, nuclei: Some(Nuclei { density: nd, radius: 1e-5, p_ref: 1e5 }), cnt_scalar: None, molecule_mass: M_WATER, condensation: 1.0 };
                    s.kinetics = Some(Box::new(move |sys, u, sc, dt| cav.apply(sys, u, sc, dt)));
                }
                let mut pmin = f64::INFINITY;
                let mut steps = 0;
                while s.t < t_end {
                    let dt = s.dt_cfl().min(t_end - s.t);
                    s.step(dt);
                    steps += 1;
                    if std::env::var("TRACE").is_ok() && steps % 500 == 0 {
                        let bad: Vec<usize> = (0..n).filter(|&i| !sys.admissible(&s.prim(i, 0, 0))).collect();
                        eprintln!("    step {steps}: t {:.3e} dt {dt:.3e} centre p {:.3e}, inadmissible {:?}", s.t, s.prim(n / 2, 0, 0)[5], &bad[..bad.len().min(4)]);
                        if let Some(&i) = bad.first() {
                            eprintln!("      cell {i}: {:?}", s.prim(i, 0, 0));
                        }
                    }
                    pmin = pmin.min((0..n).map(|i| s.prim(i, 0, 0)[5]).fold(f64::INFINITY, f64::min));
                }
                let w = s.prim(n / 2, 0, 0);
                let width = (0..n).filter(|&i| 1.0 - s.prim(i, 0, 0)[6] > 0.1).count() as f64 / n as f64;
                println!("  {name:<28}: centre p {:>10.1} Pa, α_v {:.4}; lowest p during run {:>11.4e} Pa; α_v > 0.1 over {width:.3} m", w[5], 1.0 - w[6], pmin);
            }
        }
        "henry" => {
            // N₂ bubbles (α_g = 3 %, 10⁹ m⁻³) in a slab of degassed water, open
            // ends at 1 bar; Henry k_H = 1.79e-7 kg m⁻³ Pa⁻¹ (18 mg/L at 1 atm),
            // D = 2e-9 m²/s, Sh = 2
            let n2 = Nasg { gamma: 1.4, p_inf: 0.0, b: 0.0, q: 0.0, cv: 743.0, q_prime: 0.0 };
            let sys = FiveEq { eos: [WATER_SG, n2], kapila: true };
            let ag0: f64 = std::env::var("AG0").ok().map(|v| v.parse().unwrap()).unwrap_or(0.03);
            let (rho_w, rho_g, p0) = (1000.0, 1.16, 1e5);
            let kh: f64 = std::env::var("KH").ok().map(|v| v.parse().unwrap()).unwrap_or(1.79e-7);
            let dis = Dissolution { liquid: 0, scalar: 0, henry: kh, diffusivity: 2e-9, sherwood: 2.0, bubble_density: 1e9, present: 1e-9 };
            let len = 0.1;
            let init = move |x: [f64; 3]| -> [f64; 7] {
                let ag = 1e-8 + (ag0 - 1e-8) * 0.5 * (((x[0] - 0.04) / 0.002).tanh() - ((x[0] - 0.06) / 0.002).tanh());
                [(1.0 - ag) * rho_w, ag * rho_g, 0.0, 0.0, 0.0, p0, 1.0 - ag]
            };
            let mut s = Solver::new(sys, Mesh::new([n, 1, 1], [0.0; 3], [len, 1.0, 1.0]), init);
            let far = sys.to_cons(&init([0.0; 3]));
            s.bc[0] = [Bc::Characteristic(far); 2];
            s.allspeed = Some(AllSpeedOpts { min_gain: 0.0, ..Default::default() });
            s.scalars.push(Scalar::new(&s, "dissolved N₂", Carrier::Phase(0), dis.diffusivity, |_| 0.0));
            s.kinetics = Some(Box::new(move |sys, u, sc, dt| dis.apply(sys, u, sc, dt)));
            let n2_total = |s: &Solver<7, FiveEq>| -> f64 {
                let vol = len / n as f64;
                (0..n).map(|i| s.u[s.mesh.idx(i + s.mesh.g(0), 0, 0)][1] * vol).sum::<f64>() + s.scalar_total(0)
            };
            let m0 = n2_total(&s);
            // independent ODE of the same transfer law at the slab centre
            // (dissolved concentration c, gas fraction from the N₂ balance)
            let ode = |t_end: f64| -> f64 {
                let (mut c, mut t) = (0.0f64, 0.0);
                let h = 0.01;
                let rate = |c: f64| -> f64 {
                    let ag = (ag0 * rho_g - (1.0 - ag0) * c) / rho_g;
                    if ag <= 0.0 {
                        return 0.0;
                    }
                    let al = 1.0 - ag;
                    let rb = (3.0 * ag / (4.0 * std::f64::consts::PI * 1e9)).cbrt();
                    let kl = 2.0 * 2e-9 / (2.0 * rb);
                    kl * 3.0 * ag / rb * (1.79e-7 * p0 - c) / al
                };
                while t < t_end - 1e-12 {
                    let k1 = rate(c);
                    let k2 = rate(c + 0.5 * h * k1);
                    let k3 = rate(c + 0.5 * h * k2);
                    let k4 = rate(c + h * k3);
                    c += h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
                    t += h;
                }
                c
            };
            let cs = 1.79e-7 * p0;
            println!("N₂ bubbles (α_g {ag0}, 1e9 m⁻³) in degassed water, open ends at 1 bar; c_s = {:.2} mg/L", cs * 1e3);
            println!("    t (s)   centre c (mg/L)   ODE (mg/L)   centre α_g    centre p (Pa)   N₂ inventory change");
            for t_next in [10.0, 30.0, 100.0, 300.0, 1000.0] {
                while s.t < t_next - 1e-9 {
                    let dt = (t_next - s.t).min(std::env::var("DT").ok().map(|v| v.parse().unwrap()).unwrap_or(2.0));
                    s.step(dt);
                    if std::env::var("TRACE").is_ok() {
                        let (mut pmin, mut pmax, mut umax) = (f64::INFINITY, f64::NEG_INFINITY, 0.0f64);
                        let mut first_bad = None;
                        for i in 0..n {
                            let w = s.prim(i, 0, 0);
                            pmin = pmin.min(w[5]);
                            pmax = pmax.max(w[5]);
                            umax = umax.max(w[2].abs());
                            if first_bad.is_none() && !w[5].is_finite() {
                                first_bad = Some(i);
                            }
                        }
                        let wc = s.prim(n / 2, 0, 0);
                        let we = s.prim(n * 4 / 10, 0, 0);
                        eprintln!("  t {:.3}: p in [{pmin:.6e}, {pmax:.6e}], |u| max {umax:.3e}, first NaN cell {first_bad:?}, fallbacks {}; centre α_g {:.6} p {:.2}; x=0.04 α_g {:.3e} p {:.2}", s.t, s.allspeed_stats.explicit_fallbacks, 1.0 - wc[6], wc[5], 1.0 - we[6], we[5]);
                    }
                }
                let t_out = t_next;
                let i = n / 2;
                let w = s.prim(i, 0, 0);
                let c = s.scalars[0].phi(&s, i, 0, 0) * w[0] / w[6];
                println!("  {t_out:7.0}   {:13.3}     {:9.3}     {:.5}      {:12.1}    {:+.2e}", c * 1e3, ode(t_out) * 1e3, 1.0 - w[6], w[5], n2_total(&s) / m0 - 1.0);
            }
            let st = s.allspeed_stats;
            println!("  all-speed: {} solves, {:.0} CG iterations/solve, max acoustic Courant {:.1e}", st.solves, st.cg_iterations as f64 / st.solves.max(1) as f64, st.max_acoustic_cfl);
        }
        _ => panic!("unknown case"),
    }
}
