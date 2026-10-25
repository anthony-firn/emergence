//! Two-phase (five-equation) validation.
//!
//!   multiphase_validate advect     water/air interface in uniform p, u (Abgrall test)
//!   multiphase_validate tube       water (1 GPa) | air (1 bar) shock tube vs exact
//!   multiphase_validate wood       sound speed of the mixture vs Wood's formula
use compressible::eos::materials::{AIR, WATER_SG};
use compressible::exact::{sample, star, Side};
use compressible::multiphase::FiveEq;
type V = [f64; 7];
const MAXV: usize = 7;
use compressible::{Bc, FluxMode, Mesh, Solver, System};
use std::time::Instant;

const EPS: f64 = 1e-8;

/// primitive state: water fraction α, phase densities, velocity, pressure
fn prim(alpha: f64, rho_w: f64, rho_a: f64, u: f64, p: f64) -> V {
    let mut w = [0.0; MAXV];
    w[0] = alpha * rho_w;
    w[1] = (1.0 - alpha) * rho_a;
    w[2] = u;
    w[5] = p;
    w[6] = alpha;
    w
}

fn advect() {
    let sys = FiveEq { eos: [WATER_SG, AIR], kapila: true };
    let (u0, p0) = (100.0, 1e5);
    for mode in [FluxMode::Hybrid, FluxMode::Weno] {
        for n in [100usize, 200] {
            let mesh = Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]);
            let mut s = Solver::new(sys, mesh, |x| {
                let a = if (0.25..0.75).contains(&x[0]) { 1.0 - EPS } else { EPS };
                prim(a, 1000.0, 1.2, u0, p0)
            });
            s.mode = mode;
            s.cfl = 0.4;
            s.run_to(1.0 / u0); // one period
            let (mut dp, mut du, mut da) = (0.0f64, 0.0f64, 0.0f64);
            for i in 0..n {
                let w = s.prim(i, 0, 0);
                dp = dp.max((w[5] - p0).abs() / p0);
                du = du.max((w[2] - u0).abs() / u0);
                let x = s.mesh.x(i, 0, 0)[0];
                let a = if (0.25..0.75).contains(&x) { 1.0 - EPS } else { EPS };
                da += (w[6] - a).abs() / n as f64;
            }
            println!(
                "  advect {:?} N={n}: max |Δp|/p {dp:.2e}, max |Δu|/u {du:.2e}, L1(α) {da:.3e}, shock-path faces {:.1}%",
                mode,
                100.0 * s.stats.weno_faces as f64 / s.stats.faces as f64
            );
        }
    }
}

fn tube() {
    let sys = FiveEq { eos: [WATER_SG, AIR], kapila: true };
    let (rw, pw, ra, pa) = (1000.0, 1e9, 50.0, 1e5);
    let sl = Side { rho: rw, u: 0.0, p: pw, eos: WATER_SG };
    let sr = Side { rho: ra, u: 0.0, p: pa, eos: AIR };
    let st = star(&sl, &sr);
    let (x0, t_end) = (0.7, 2.4e-4);
    println!("water|air tube: p* = {:.4e} Pa, u* = {:.2} m/s", st.p, st.u);
    for mode in [FluxMode::Hybrid, FluxMode::Weno] {
        let mut prev: Option<f64> = None;
        for n in [200usize, 400, 800, 1600] {
            let mesh = Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]);
            let mut s = Solver::new(sys, mesh, |x| if x[0] < x0 { prim(1.0 - EPS, rw, ra, 0.0, pw) } else { prim(EPS, rw, ra, 0.0, pa) });
            s.bc[0] = [Bc::Outflow; 2];
            s.mode = mode;
            s.cfl = 0.4;
            let t0 = Instant::now();
            s.run_to(t_end);
            let (mut e1, mut ep, mut scale) = (0.0, 0.0, 0.0);
            for i in 0..n {
                let x = s.mesh.x(i, 0, 0)[0];
                let (re, _, pe, _) = sample(&sl, &sr, &st, (x - x0) / t_end);
                let w = s.prim(i, 0, 0);
                e1 += (FiveEq::rho(&w) - re).abs() / n as f64;
                ep += (w[5] - pe).abs() / n as f64;
                scale += re / n as f64;
            }
            let rate = prev.map(|p| format!("{:.2}", (p / e1).log2())).unwrap_or("-".into());
            prev = Some(e1);
            println!(
                "  {:?} N={n:5}: L1(ρ)/ρ̄ {:.3e} rate {rate:>5}  L1(p)/p_L {:.3e}  shock-path faces {:5.1}%  fallback {}  {:.2}s",
                mode,
                e1 / scale,
                ep / pw,
                100.0 * s.stats.weno_faces as f64 / s.stats.faces as f64,
                s.stats.fallback_faces,
                t0.elapsed().as_secs_f64()
            );
        }
    }
}

fn wood() {
    let sys = FiveEq { eos: [WATER_SG, AIR], kapila: true };
    let p = 1e5;
    println!("  sound speed of water/air mixture at 1 bar:");
    for a_air in [0.0, 1e-4, 1e-2, 0.1, 0.5, 0.9, 1.0] {
        let alpha = (1.0 - a_air as f64).clamp(EPS, 1.0 - EPS);
        let w = prim(alpha, 1000.0, 1.2, 0.0, p);
        let rho = FiveEq::rho(&w);
        let zw = WATER_SG.gamma * (p + WATER_SG.p_inf);
        let za = AIR.gamma * p;
        let wood = (1.0 / (rho * (alpha / zw + (1.0 - alpha) / za))).sqrt();
        println!("    air fraction {a_air:7}: model c = {:9.3} m/s, Wood = {:9.3} m/s", sys.sound_speed(&w), wood);
    }
}

fn main() {
    match std::env::args().nth(1).as_deref() {
        Some("advect") => advect(),
        Some("tube") => tube(),
        Some("wood") => wood(),
        _ => {
            wood();
            advect();
            tube();
        }
    }
}
