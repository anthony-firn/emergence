//! Resolved viscous shock vs Becker's exact solution (Pr = 3/4, constant μ,
//! no bulk viscosity): with total enthalpy conserved through the shock,
//!   (4μ/3m) u′ = ((γ+1)/2γ)(u − u₁)(u − u₂)/u,
//! integrating to x(u) = C (u₁ ln(u₁ − u) − u₂ ln(u − u₂))/(u₁ − u₂) + x₀.
//!
//!   becker [mach] [mu]
use compressible::eos::Nasg;
type V = [f64; 5];
const MAXV: usize = 5;
use compressible::viscous::Transport;
use compressible::{Bc, Euler, FluxMode, Mesh, Solver};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mach: f64 = a.get(1).map(|v| v.parse().unwrap()).unwrap_or(2.0);
    let mu: f64 = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(0.02);
    let g = 1.4;
    let (rho1, p1) = (1.0, 1.0 / g);
    let u1 = mach; // c₁ = 1
    let m = rho1 * u1;
    let u2 = u1 * ((g - 1.0) / (g + 1.0) + 2.0 / ((g + 1.0) * mach * mach));
    let h = g / (g - 1.0) * p1 / rho1 + 0.5 * u1 * u1;
    let c = 4.0 * mu / (3.0 * m) * 2.0 * g / (g + 1.0);
    let xu = |u: f64| c * (u1 * (u1 - u).ln() - u2 * (u - u2).ln()) / (u1 - u2);
    let x0 = xu(0.5 * (u1 + u2));
    let exact = |x: f64| -> V {
        // invert x(u) by bisection (x decreasing in u)
        let (mut lo, mut hi) = (u2 * (1.0 + 1e-15), u1 * (1.0 - 1e-15));
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if xu(mid) - x0 > x { lo = mid } else { hi = mid }
        }
        let u = 0.5 * (lo + hi);
        let rho = m / u;
        let p = rho * (g - 1.0) / g * (h - 0.5 * u * u);
        let mut w = [0.0; MAXV];
        w[..5].copy_from_slice(&[rho, u, 0.0, 0.0, p]);
        w
    };
    let delta = (u1 - u2) / ((g + 1.0) / (2.0 * g) * (u1 - u2).powi(2) / 4.0 / ((u1 + u2) / 2.0) * 3.0 * m / (4.0 * mu)); // max-slope thickness
    println!("Mach {mach}, mu {mu}: shock thickness (max slope) {delta:.4}");
    let modes = match std::env::var("MODE").as_deref() {
        Ok("weno") => vec![FluxMode::Weno],
        Ok("hybrid") => vec![FluxMode::Hybrid],
        Ok("central") => vec![FluxMode::Central],
        _ => vec![FluxMode::Hybrid, FluxMode::Weno],
    };
    for mode in modes {
        let mut prev: Option<f64> = None;
        for n in [50usize, 100, 200, 400] {
            let len = 20.0 * delta;
            let mesh = Mesh::new([n, 1, 1], [-0.5 * len, 0.0, 0.0], [0.5 * len, 1.0, 1.0]).with_ghosts(6);
            let mut s = Solver::new(Euler { eos: Nasg { cv: 1.0, ..Nasg::ideal(g) } }, mesh, |x| exact(x[0]));
            let sys = Euler { eos: Nasg { cv: 1.0, ..Nasg::ideal(g) } };
            use compressible::System;
            s.bc[0] = [Bc::Dirichlet(sys.to_cons(&exact(-0.5 * len - 3.0 * s.mesh.dx[0]))), Bc::Dirichlet(sys.to_cons(&exact(0.5 * len + 3.0 * s.mesh.dx[0])))];
            s.mode = mode;
            s.viscous = Some(Transport::constant(mu, 0.75));
            s.cfl = 0.3;
            let t_end = 2.0 * len / u1;
            if std::env::var("DEBUG").is_ok() {
                for it in 0..200 {
                    let dt = s.dt_cfl();
                    if it % 20 == 0 || !(dt > 1e-8) {
                        let rmax = (0..n).map(|i| s.prim(i, 0, 0)[0]).fold(f64::NAN, f64::max);
                        let e1: f64 = (0..n).map(|i| (s.prim(i, 0, 0)[0] - exact(s.mesh.x(i, 0, 0)[0])[0]).abs()).sum::<f64>() / n as f64;
                        let (imin, rmin) = (0..n).map(|i| (i, s.prim(i, 0, 0)[0])).fold((0, f64::INFINITY), |a, b| if b.1 < a.1 { b } else { a });
                        let w0 = s.prim(0, 0, 0); let wn = s.prim(n - 1, 0, 0);
                        println!("it {it} dt {dt:e} t {:.4} rho_max {rmax:.4} rho_min {rmin:.4}@{imin} L1 {e1:.3e} left (rho {:.4} u {:.4} p {:.4}) right (rho {:.4} u {:.4} p {:.4})", s.t, w0[0], w0[1], w0[4], wn[0], wn[1], wn[4]);
                        if !(dt > 1e-8) { break; }
                    }
                    s.step(dt);
                }
                return;
            }
            s.run_to(t_end);
            let e1: f64 = (0..n).map(|i| (s.prim(i, 0, 0)[0] - exact(s.mesh.x(i, 0, 0)[0])[0]).abs()).sum::<f64>() / n as f64;
            let rate = prev.map(|p| format!("{:.2}", (p / e1).log2())).unwrap_or("-".into());
            prev = Some(e1);
            println!(
                "  {:?} N={n:4} (dx/δ = {:.3}): L1(ρ) {:.3e} rate {rate:>5}  WENO faces {:5.1}%",
                mode,
                len / n as f64 / delta,
                e1,
                100.0 * s.stats.weno_faces as f64 / s.stats.faces as f64
            );
        }
    }
}
