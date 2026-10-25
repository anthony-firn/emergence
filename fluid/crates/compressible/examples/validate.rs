//! Validation of the compressible scheme against exact solutions.
//!
//!   validate riemann   — Sod, Lax, Toro-3 and a water shock tube vs the exact Riemann solution
//!   validate vortex    — isentropic vortex (smooth, exact): order of accuracy
//!   validate shuosher  — shock / entropy-wave interaction vs a fine reference
use compressible::eos::{materials, Nasg};
use compressible::exact::{sample, star, Side};
type V = [f64; 5];
const MAXV: usize = 5;
use compressible::scheme::ShockScheme;
use compressible::{Bc, Euler, FluxMode, Mesh, Solver};

fn shock() -> ShockScheme {
    if std::env::var("SHOCK").as_deref() == Ok("fdweno") { ShockScheme::FdWeno } else { ShockScheme::Wcns }
}
use std::time::Instant;

fn prim(rho: f64, u: f64, v: f64, p: f64) -> V {
    let mut w = [0.0; MAXV];
    w[..5].copy_from_slice(&[rho, u, v, 0.0, p]);
    w
}

fn riemann_case(name: &str, eos: Nasg, l: (f64, f64, f64), r: (f64, f64, f64), x0: f64, t_end: f64) {
    let sl = Side { rho: l.0, u: l.1, p: l.2, eos };
    let sr = Side { rho: r.0, u: r.1, p: r.2, eos };
    let st = star(&sl, &sr);
    println!("{name}: p* = {:.6e}, u* = {:.6}", st.p, st.u);
    for mode in [FluxMode::Weno, FluxMode::Hybrid] {
        let mut prev: Option<f64> = None;
        for n in [100usize, 200, 400, 800, 1600] {
            let mesh = Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]);
            let mut s = Solver::new(Euler { eos }, mesh, |x| if x[0] < x0 { prim(l.0, l.1, 0.0, l.2) } else { prim(r.0, r.1, 0.0, r.2) });
            s.bc[0] = [Bc::Outflow; 2];
            s.mode = mode;
        s.shock = shock();
            s.shock = shock();
            s.cfl = 0.4;
            let t0 = Instant::now();
            let steps = s.run_to(t_end);
            let el = t0.elapsed().as_secs_f64();
            let (mut e1, mut scale) = (0.0, 0.0);
            for i in 0..n {
                let x = s.mesh.x(i, 0, 0)[0];
                let (re, _, _, _) = sample(&sl, &sr, &st, (x - x0) / t_end);
                e1 += (s.prim(i, 0, 0)[0] - re).abs() / n as f64;
                scale += re.abs() / n as f64;
            }
            let rate = prev.map(|p| format!("{:.2}", (p / e1).log2())).unwrap_or("-".into());
            prev = Some(e1);
            println!(
                "  {:?} N={n:5}: L1(ρ)/|ρ| {:.3e}  rate {rate:>5}  WENO faces {:5.1}%  fallback {}  {steps} steps {:.2}s",
                mode,
                e1 / scale,
                100.0 * s.stats.weno_faces as f64 / s.stats.faces as f64,
                s.stats.fallback_faces,
                el
            );
        }
    }
}

fn vortex() {
    // Shu's isentropic vortex, γ = 1.4, strength β = 5, mean flow (1, 1),
    // periodic [−5, 5]²; exact solution = initial condition translated.
    let g = 1.4;
    let beta = 5.0;
    let pi = std::f64::consts::PI;
    let exact = move |x: f64, y: f64, t: f64| -> V {
        let wrap = |a: f64| (a + 5.0).rem_euclid(10.0) - 5.0;
        let (xr, yr) = (wrap(x - t), wrap(y - t));
        let r2 = xr * xr + yr * yr;
        let f = beta / (2.0 * pi) * (0.5 * (1.0 - r2)).exp();
        let tt = 1.0 - (g - 1.0) * beta * beta / (8.0 * g * pi * pi) * (1.0 - r2).exp();
        let rho = tt.powf(1.0 / (g - 1.0));
        prim(rho, 1.0 - f * yr, 1.0 + f * xr, rho * tt)
    };
    let t_end = 2.0;
    for mode in [FluxMode::Central, FluxMode::Hybrid, FluxMode::Weno] {
        let mut prev: Option<f64> = None;
        for n in [20usize, 40, 80, 160] {
            let mesh = Mesh::new([n, n, 1], [-5.0, -5.0, 0.0], [5.0, 5.0, 1.0]);
            let mut s = Solver::new(Euler { eos: Nasg::ideal(g) }, mesh, |x| exact(x[0], x[1], 0.0));
            s.mode = mode;
        s.shock = shock();
            s.shock = shock();
            s.cfl = 0.4;
            let t0 = Instant::now();
            // fixed dt ∝ dx so the time error (3rd order) stays below the space error
            let dt = 0.4 * s.mesh.dx[0] / 3.0 * (20.0 / n as f64).powf(2.0 / 3.0).min(1.0);
            let steps = (t_end / dt).ceil() as usize;
            for _ in 0..steps {
                s.step(t_end / steps as f64);
            }
            let el = t0.elapsed().as_secs_f64();
            let mut e1 = 0.0;
            for i in 0..n {
                for j in 0..n {
                    let x = s.mesh.x(i, j, 0);
                    e1 += (s.prim(i, j, 0)[0] - exact(x[0], x[1], t_end)[0]).abs() / (n * n) as f64;
                }
            }
            let rate = prev.map(|p| format!("{:.2}", (p / e1).log2())).unwrap_or("-".into());
            prev = Some(e1);
            println!(
                "  vortex {:?} N={n:4}: L1(ρ) {:.3e}  rate {rate:>5}  WENO faces {:5.1}%  {:.2}s",
                mode,
                e1,
                100.0 * s.stats.weno_faces as f64 / s.stats.faces as f64,
                el
            );
        }
    }
}

fn shu_osher() {
    let eos = Nasg::ideal(1.4);
    let init = |x: [f64; 3]| {
        if x[0] < -4.0 {
            prim(3.857143, 2.629369, 0.0, 10.33333)
        } else {
            prim(1.0 + 0.2 * (5.0 * x[0]).sin(), 0.0, 0.0, 1.0)
        }
    };
    let run = |n: usize, mode: FluxMode| {
        let mesh = Mesh::new([n, 1, 1], [-5.0, 0.0, 0.0], [5.0, 1.0, 1.0]);
        let mut s = Solver::new(Euler { eos }, mesh, init);
        s.bc[0] = [Bc::Outflow; 2];
        s.mode = mode;
        s.shock = shock();
        s.cfl = 0.4;
        let t0 = Instant::now();
        s.run_to(1.8);
        (s, t0.elapsed().as_secs_f64())
    };
    let (rf, _) = run(6400, FluxMode::Weno);
    let refr: Vec<f64> = (0..6400).map(|i| rf.prim(i, 0, 0)[0]).collect();
    for mode in [FluxMode::Weno, FluxMode::Hybrid] {
        for n in [200usize, 400, 800, 1600] {
            let (s, el) = run(n, mode);
            let r = 6400 / n;
            let e1: f64 = (0..n).map(|i| (s.prim(i, 0, 0)[0] - refr[i * r..(i + 1) * r].iter().sum::<f64>() / r as f64).abs()).sum::<f64>() / n as f64;
            println!(
                "  shu-osher {:?} N={n:5}: L1(ρ) vs 6400-cell reference {:.3e}  WENO faces {:5.1}%  {:.2}s",
                mode,
                e1,
                100.0 * s.stats.weno_faces as f64 / s.stats.faces as f64,
                el
            );
        }
    }
}

fn main() {
    let which = std::env::args().nth(1).unwrap_or("all".into());
    if which == "riemann" || which == "all" {
        let gas = Nasg::ideal(1.4);
        riemann_case("Sod", gas, (1.0, 0.0, 1.0), (0.125, 0.0, 0.1), 0.5, 0.2);
        riemann_case("Lax", gas, (0.445, 0.698, 3.528), (0.5, 0.0, 0.571), 0.5, 0.14);
        riemann_case("Toro test 3 (p ratio 1e5)", gas, (1.0, 0.0, 1000.0), (1.0, 0.0, 0.01), 0.5, 0.012);
        riemann_case("water shock tube (stiffened gas, 1 GPa : 1 bar)", materials::WATER_SG, (1000.0, 0.0, 1e9), (1000.0, 0.0, 1e5), 0.5, 1.2e-4);
    }
    if which == "vortex" || which == "all" {
        vortex();
    }
    if which == "shuosher" || which == "all" {
        shu_osher();
    }
}
