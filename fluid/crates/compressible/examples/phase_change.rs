//! State changes driven by heat (`isobaric.rs`): melting, freezing,
//! condensation fronts against the exact Neumann solution, and freezing of
//! supercooled water droplets by homogeneous nucleation.
//!
//!   phase_change stefan [N]     model material, iron, lake ice, nitrogen
//!   phase_change droplet        supercooled water droplets cooled at 1 K/min
//!   phase_change allotropes     iron (bcc→fcc→bcc→liquid→vapour) and tin
//!                               heated and cooled: thermal arrests, volume jumps
//!
//! Stefan problems: a slab at its transition temperature touches a wall
//! held at T_w; the new phase grows from the wall as conduction removes or
//! brings the latent heat. The new phase next to the fixed wall is at rest
//! (density changes move only the old phase), so the one-phase Neumann
//! solution s = 2λ√(αt), λe^{λ²}erf λ = St/√π, St = c_p|T_w − T_t|/L holds
//! with the new phase's properties even when the densities differ.
use compressible::allotropes;
use compressible::elements::{by_symbol, fit, fit_inputs, FitInputs, Grade, ThreePhase, P_ATM};
use compressible::isobaric::{enthalpy, Freezing, Isobaric, PhaseSet};
use compressible::phase_change::v_pt;
use compressible::substances::by_name;

fn erf(x: f64) -> f64 {
    if x.abs() >= 3.0 {
        return 1.0f64.copysign(x);
    }
    let (mut sum, mut term, x2) = (x, x, x * x);
    for n in 1..200 {
        term *= -x2 / n as f64;
        let add = term / (2 * n + 1) as f64;
        sum += add;
        if add.abs() < 1e-17 * sum.abs() {
            break;
        }
    }
    2.0 / std::f64::consts::PI.sqrt() * sum
}

fn neumann_lambda(st: f64) -> f64 {
    let f = |l: f64| l * (l * l).exp() * erf(l) - st / std::f64::consts::PI.sqrt();
    let (mut lo, mut hi) = (1e-8, 5.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

/// One Stefan case: `old` phase at its transition with `new`, wall at
/// T_t + dt_wall.
#[allow(clippy::too_many_arguments)]
fn stefan(name: &str, model: ThreePhase, k: [f64; 3], old: usize, new: usize, dt_wall: f64, n: usize, t_range: (f64, f64)) {
    let eos = [model.solid, model.liquid, model.vapor];
    let lad = compressible::isobaric::Ladder::new(&PhaseSet::from(&model), P_ATM, t_range.0, t_range.1);
    let tt = lad.transitions.iter().find(|x| (x.1 == old && x.2 == new) || (x.1 == new && x.2 == old)).expect("transition").0;
    let (en, eo) = (eos[new], eos[old]);
    let latent = (enthalpy(&en, P_ATM, tt) - enthalpy(&eo, P_ATM, tt)).abs();
    let cp = en.gamma * en.cv;
    let rho = 1.0 / v_pt(&en, P_ATM, tt);
    let alpha = k[new] / (rho * cp);
    let st = cp * dt_wall.abs() / latent;
    let lam = neumann_lambda(st);
    // slab long enough that the front reaches about a third of it
    let t_end = 1.0;
    let s_end = 2.0 * lam * (alpha * t_end).sqrt();
    let rho_old = 1.0 / v_pt(&eo, P_ATM, tt);
    let len_old = 3.0 * s_end * rho / rho_old;
    // start a hair inside the old phase's single-phase range
    let t0 = tt * (1.0 - 1e-9 * (new as f64 - old as f64).signum());
    let mut s = Isobaric::new(PhaseSet::from(&model), P_ATM, k.to_vec(), n, len_old / n as f64, old, t0, t_range);
    s.wall[0] = Some(tt + dt_wall);
    println!(
        "{name}: transition {tt:.2} K ({}→{}), wall {:+.0} K, L {:.1} kJ/kg, ρ_old {rho_old:.4} / ρ_new {rho:.4} kg/m³, α_new {alpha:.3e} m²/s, St {st:.4}, λ {lam:.5}",
        ["solid", "liquid", "vapour"][old],
        ["solid", "liquid", "vapour"][new],
        dt_wall,
        latent / 1e3
    );
    let mut steps = 0;
    for i in 1..=5 {
        let t_out = t_end * i as f64 / 5.0;
        steps += s.advance(t_out - s.t);
        // thickness of the new phase
        let temps = s.temperatures();
        let front: f64 = (0..n).map(|c| s.m[c] * s.y[c][new] * v_pt(&eos[new], P_ATM, temps[c])).sum();
        let exact = 2.0 * lam * (alpha * s.t).sqrt();
        let (w, x) = s.geometry();
        let _ = w;
        // temperature halfway into the new phase vs the exact profile
        let c = x.iter().position(|&xc| xc > 0.5 * exact).unwrap_or(0);
        let tex = tt + dt_wall - dt_wall * erf(x[c] / (2.0 * (alpha * s.t).sqrt())) / erf(lam);
        println!(
            "    t {:5.2} s: front {:8.4} mm (exact {:8.4}, {:+6.2}%); T at {:.4} mm {:.3} K (exact {:.3})",
            s.t,
            front * 1e3,
            exact * 1e3,
            100.0 * (front / exact - 1.0),
            x[c] * 1e3,
            temps[c],
            tex
        );
    }
    println!("    {steps} conduction steps, {n} cells");
}

/// Heat one cell of a substance at constant power from `t0` to `t1` and
/// cool it back: the temperature halts at each transition for as long as
/// the latent heat takes (thermal arrests, the classical way the iron
/// transitions were found), and the specific volume jumps (dilatometry).
/// `dv`: measured volume jumps of the solid–solid transitions (%) to
/// compare.
fn arrests(set: &PhaseSet, molar: f64, range: (f64, f64), t0: f64, t1: f64, dv: &[f64]) {
    let p = P_ATM;
    let lowest = compressible::isobaric::Ladder::new(set, p, range.0, range.1).lowest;
    let mut s = Isobaric::new(set.clone(), p, vec![1.0; set.len()], 1, 1.0, lowest, t0, range);
    s.m[0] = 1.0;
    let lad = s.ladder.clone();
    let mut nss = 0;
    for &(tt, a, b) in &lad.transitions {
        let dh = (set.enthalpy(b, p, tt) - set.enthalpy(a, p, tt)) * molar;
        let (va, vb) = (set.specific_volume(a, p, tt), set.specific_volume(b, p, tt));
        let meas = if set.is_solid(a) && set.is_solid(b) { dv.get(nss) } else { None };
        nss += (set.is_solid(a) && set.is_solid(b)) as usize;
        println!(
            "    {:>10} → {:<10} at {:8.2} K: ΔH {:+8.3} kJ/mol{}{}",
            set.names[a],
            set.names[b],
            tt,
            dh / 1e3,
            if b == set.vapor { String::new() } else { format!(", ΔV/V {:+6.2}%", 100.0 * (vb / va - 1.0)) },
            meas.map_or(String::new(), |x| format!(" (measured ≈ {x:+.2}%)"))
        );
    }
    // heating at 1 kW/kg, then cooling at 1 kW/kg; time at each transition
    let power = 1000.0;
    for (label, sign, t_stop) in [("heating", -1.0, t1), ("cooling", 1.0, t0)] {
        s.sink = sign * power;
        let mut stay = vec![0.0; lad.transitions.len()];
        let dt = 0.05;
        loop {
            let tc = s.temperatures()[0];
            if (sign < 0.0 && tc >= t_stop) || (sign > 0.0 && tc <= t_stop) {
                break;
            }
            s.advance(dt);
            let tn = s.temperatures()[0];
            for (j, &(tt, _, _)) in lad.transitions.iter().enumerate() {
                if (tn - tt).abs() < 1e-6 * tt && (tc - tt).abs() < 1e-6 * tt {
                    stay[j] += dt;
                }
            }
        }
        let arr: Vec<String> = lad.transitions.iter().zip(&stay).map(|(x, &d)| format!("{:.0} K {:.1} s", x.0, d)).collect();
        println!("    {label} at 1 kW/kg — arrests: {} (latent heat = power × arrest)", arr.join(", "));
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let case = a.get(1).cloned().unwrap_or("stefan".into());
    let n: usize = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(400);
    match case.as_str() {
        "stefan" => {
            // model material: iron-like, equal densities, negligible expansion
            let m = 55.845e-3;
            let model = fit_inputs(&FitInputs {
                molar_mass: m,
                t_melt: 1811.0,
                t_boil: 3134.0,
                sublimes: false,
                p_triple: None,
                rho_s: 7000.0,
                t_ref_s: 1811.0,
                alpha_s: 1e-7,
                c_bulk_s: 4000.0,
                shear_speed: 2000.0,
                rho_l: 7000.0,
                t_ref_l: 1811.0,
                alpha_l: 1e-7,
                c_l: 4000.0,
                cp_s: 800.0,
                cp_l: 800.0,
                gamma_v: 5.0 / 3.0,
                cv_v: 1.5 * 8.314 / m,
                h_fus: 247e3,
                h_vap: 6.1e6,
                t_crit: None,
                grade: Grade::C,
                estimated: vec![],
            });
            stefan("model material (equal densities)", model, [35.0, 35.0, 0.05], 0, 1, 189.0, n, (500.0, 5000.0));
            stefan("iron melting (periodic-table fit)", fit(by_symbol("Fe").unwrap()).unwrap(), [33.0, 35.0, 0.05], 0, 1, 189.0, n, (500.0, 5000.0));
            stefan("lake ice (water at 0 °C, wall −10 °C)", by_name("water").unwrap().model(), [2.22, 0.561, 0.025], 1, 0, -10.0, n, (150.0, 1000.0));
            stefan("liquid nitrogen condensing on a wall 10 K below boiling", fit(by_symbol("N").unwrap()).unwrap(), [0.25, 0.14, 0.0075], 2, 1, -10.0, n, (30.0, 300.0));
            // silicon: conductivity 22 W/(m K) solid at the melting point,
            // ~56 liquid (liquid silicon is a metal)
            let si = fit(by_symbol("Si").unwrap()).unwrap();
            stefan("silicon melting (periodic-table fit)", si.clone(), [22.0, 56.0, 0.05], 0, 1, 100.0, n, (500.0, 5000.0));
            // liquid silicon is denser than the solid: the melting point
            // falls with pressure (Clausius–Clapeyron from the two fits)
            let set = PhaseSet::from(&si);
            let tm = |p: f64| compressible::isobaric::Ladder::new(&set, p, 500.0, 5000.0).transitions[0].0;
            let (t1, t2) = (tm(P_ATM), tm(P_ATM + 100e6));
            let (vs, vl) = (v_pt(&si.solid, P_ATM, t1), v_pt(&si.liquid, P_ATM, t1));
            println!(
                "silicon: ρ_solid {:.0} / ρ_liquid {:.0} kg/m³ at the melting point; melting point {t1:.2} K at 1 atm, {t2:.2} K at 100 MPa: slope {:.4} K/MPa (Si–I/Si–II/liquid triple point ≈ 10.5 GPa, ≈ 1000 K: average slope ≈ −0.065 K/MPa)",
                1.0 / vs,
                1.0 / vl,
                (t2 - t1) / 100.0
            );
        }
        "droplet" => {
            // homogeneous ice nucleation: Koop et al. (2000), J from the water
            // activity shift Δa_w = 1 − p_ice(T)/p_liq(T), with both vapour
            // pressures from the fitted water model (sublimation and the
            // supercooled liquid's boiling curve) — no fitted freezing data
            let water = by_name("water").unwrap();
            let model = water.model();
            let eos = [model.solid, model.liquid, model.vapor];
            let psat = move |a: usize, t: f64| -> f64 {
                // coexistence pressure of phase a with the vapour at T
                let (mut lo, mut hi) = (1e-6f64.ln(), 1e6f64.ln());
                for _ in 0..200 {
                    let mid = 0.5 * (lo + hi);
                    let p = mid.exp();
                    if eos[a].gibbs(p, t) < eos[2].gibbs(p, t) {
                        hi = mid; // condensed phase stable: pressure above saturation
                    } else {
                        lo = mid;
                    }
                }
                (0.5 * (lo + hi)).exp()
            };
            let rate = move |t: f64| -> f64 {
                if t >= 273.15 {
                    return 0.0;
                }
                let da = 1.0 - psat(0, t) / psat(1, t);
                let lg = -906.7 + 8502.0 * da - 26924.0 * da * da + 29180.0 * da * da * da;
                10f64.powf(lg.min(15.0)) * 1e6 // cm⁻³ s⁻¹ → m⁻³ s⁻¹
            };
            let rho_l = 1.0 / v_pt(&eos[1], P_ATM, 273.15);
            let cp_l = eos[1].gamma * eos[1].cv;
            let latent = enthalpy(&eos[1], P_ATM, 273.15) - enthalpy(&eos[0], P_ATM, 273.15);
            println!("supercooled water droplets cooled at 1 K/min; homogeneous nucleation (Koop et al. 2000 with the model's ice/liquid vapour pressures)");
            println!("    a_w^i = p_ice/p_liq from the model at 238 K: {:.4} (Murphy & Koop 2005: 0.6955)", psat(0, 238.0) / psat(1, 238.0));
            // the same nucleation law with the reference vapour pressures
            // (Murphy & Koop 2005, eqs. 7 and 10): what the EOS model adds
            let mk_ratio = |t: f64| -> f64 {
                let pi = (9.550426 - 5723.265 / t + 3.53068 * t.ln() - 0.00728332 * t).exp();
                let pl = (54.842763 - 6763.22 / t - 4.210 * t.ln() + 0.000367 * t + (0.0415 * (t - 218.8)).tanh() * (53.878 - 1331.22 / t - 9.44523 * t.ln() + 0.014025 * t)).exp();
                pi / pl
            };
            let freeze_t = |d_um: f64, ratio: &dyn Fn(f64) -> f64| -> f64 {
                // deterministic onset: ∫ J V dt = 1 while cooling at 1 K/min
                let vol = std::f64::consts::PI / 6.0 * (d_um * 1e-6f64).powi(3);
                let mut n_exp = 0.0;
                let mut t = 273.15;
                while n_exp < 1.0 && t > 200.0 {
                    let da = 1.0 - ratio(t);
                    let lg = -906.7 + 8502.0 * da - 26924.0 * da * da + 29180.0 * da * da * da;
                    n_exp += 10f64.powf(lg.min(15.0)) * 1e6 * vol * 0.6; // 0.01 K steps at 1 K/min = 0.6 s
                    t -= 0.01;
                }
                t
            };
            println!("    same law with the Murphy–Koop vapour pressures: freezing at {:.2} / {:.2} / {:.2} K for 1 / 10 / 100 µm", freeze_t(1.0, &mk_ratio), freeze_t(10.0, &mk_ratio), freeze_t(100.0, &mk_ratio));
            for d_um in [1.0, 10.0, 100.0] {
                let vol = std::f64::consts::PI / 6.0 * (d_um * 1e-6f64).powi(3);
                // 0D: one cell of mass ρV (per "unit area" = the droplet)
                let mut s = Isobaric::new(PhaseSet::from(&model), P_ATM, vec![2.22, 0.561, 0.025], 1, vol, 1, 273.0, (100.0, 1000.0));
                s.m[0] = rho_l * vol;
                s.sink = cp_l / 60.0;
                s.freezing = Some(Freezing { rate: Box::new(rate), growth: Box::new(|dt| 1e-3 * dt * dt) });
                let mut t_nuc = f64::NAN;
                let mut t_after = f64::NAN;
                let mut y_after = f64::NAN;
                let mut t_low = f64::INFINITY;
                for _ in 0..200_000 {
                    let before = s.y[0][0];
                    s.advance(0.05);
                    let t = s.temperatures()[0];
                    t_low = t_low.min(t);
                    if before <= 0.0 && s.y[0][0] > 0.0 && t_nuc.is_nan() {
                        t_nuc = t_low;
                    }
                    // recalescence finished: the droplet sits at the melting point
                    if !t_nuc.is_nan() && t_after.is_nan() && (t - 273.15).abs() < 0.05 {
                        t_after = t;
                        y_after = s.y[0][0];
                        break;
                    }
                }
                let y_exact = (enthalpy(&eos[1], P_ATM, 273.15) - enthalpy(&eos[1], P_ATM, t_nuc)) / latent;
                println!(
                    "  {d_um:5.0} µm: freezes at {t_nuc:.2} K; recalescence to {t_after:.2} K with ice fraction {y_after:.4} (energy balance {y_exact:.4}, c_pΔT/L {:.4})",
                    cp_l * (273.15 - t_nuc) / latent
                );
            }
            println!("  measured homogeneous freezing of pure water droplets at ~1 K/min: about 235–237 K for 1–100 µm");
        }
        "allotropes" => {
            let (iron, dh_gd) = allotropes::iron();
            let m = 55.845e-3;
            println!("iron: bcc (α/δ), fcc (γ), liquid, vapour; bcc and fcc use the temperature-dependent Dorogokupets et al. (2017) thermal EOS plus the magnetic term; fitted to T(α→γ), T(γ→δ) and ΔH(α→γ) = 0.90 kJ/mol");
            println!("    fitted data: T(α→γ) 1185 K, T(γ→δ) 1667 K, ΔH(α→γ) 0.90 kJ/mol; the volume jumps and ΔH(γ→δ) are predictions");
            println!("    c_p bcc {:.1} J/(mol K) at 300 K, {:.1} at 1185 K (magnetic peak {:.1} at T_C 1043 K) / fcc {:.1}",
                iron.cp(0, P_ATM, 300.0) * m, iron.cp(0, P_ATM, 1185.0) * m, iron.cp(0, P_ATM, 1043.0) * m, iron.cp(1, P_ATM, 1185.0) * m);
            println!("    predicted ΔH(γ→δ) {:.3} kJ/mol (measured 0.84)", dh_gd / 1e3);
            arrests(&iron, m, (250.0, 4000.0), 300.0, 3300.0, &[-0.98, 0.53]);
            let tin = allotropes::tin();
            println!("tin: grey (α, diamond cubic), white (β, tetragonal), liquid, vapour (densities are inputs; real grey tin forms only slowly on cooling — tin pest — which equilibrium relaxation does not delay)");
            arrests(&tin, 118.71e-3, (100.0, 4000.0), 150.0, 3000.0, &[]);
        }
        _ => panic!("unknown case"),
    }
}
