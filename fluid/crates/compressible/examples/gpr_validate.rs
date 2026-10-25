//! GPR unified continuum model: solid and fluid limits.
//!
//!   gpr_validate waves    elastic solid (τ = ∞): shear waves at c_s,
//!                         longitudinal at √(c₀² + 4c_s²/3)
//!   gpr_validate stokes   Stokes' first problem: GPR → Navier–Stokes
//!                         (μ = ρ₀c_s²τ/6) as the Deborah number → 0
use compressible::eos::Nasg;
use compressible::gpr::{Gpr, NG};
use compressible::viscous::Transport;
use compressible::{Bc, Euler, Mesh, Solver};

fn prim(rho: f64, v: [f64; 3], p: f64) -> [f64; NG] {
    let mut w = [0.0; NG];
    w[0] = rho;
    w[1..4].copy_from_slice(&v);
    for i in 0..3 {
        w[4 + 3 * i + i] = 1.0;
    }
    w[13] = p;
    w
}

fn erf(x: f64) -> f64 {
    // Abramowitz–Stegun 7.1.26 is too coarse; use a series / continued fraction
    let t = x.abs();
    let r = if t < 3.0 {
        let mut sum = t;
        let mut term = t;
        let mut n = 0.0;
        while term.abs() > 1e-17 * sum.abs() {
            n += 1.0;
            term *= -t * t / n;
            sum += term / (2.0 * n + 1.0);
        }
        2.0 / std::f64::consts::PI.sqrt() * sum
    } else {
        // asymptotic continued fraction for erfc
        let mut f = 0.0;
        for k in (1..60).rev() {
            f = (k as f64 / 2.0) / (t + f);
        }
        1.0 - (-t * t).exp() / std::f64::consts::PI.sqrt() / (t + f)
    };
    r.copysign(x)
}

fn waves() {
    let eos = Nasg::ideal(1.4);
    let (rho, p, cs) = (1.0, 1.0 / 1.4, 1.0);
    let sys = Gpr { eos, cs, tau: f64::INFINITY, tau_model: compressible::gpr::TauModel::Constant };
    let c_long = (1.0f64 + 4.0 / 3.0 * cs * cs).sqrt();
    for (name, comp, exact) in [("shear", 2usize, cs), ("longitudinal", 1usize, c_long)] {
        for n in [400usize, 800, 1600] {
            let mesh = Mesh::new([n, 1, 1], [-5.0, 0.0, 0.0], [5.0, 1.0, 1.0]);
            let eps: f64 = std::env::var("EPS").ok().map(|v| v.parse().unwrap()).unwrap_or(1e-4);
            let mut s = Solver::new(sys, mesh, |x| {
                let mut v = [0.0; 3];
                v[comp - 1] = eps * (-(x[0] / 0.3).powi(2)).exp();
                prim(rho, v, p)
            });
            s.cfl = 0.4;
            let t_end = 2.0;
            s.run_to(t_end);
            // right-moving pulse: peak of |v_comp| for x > 0, refined by a parabola
            let vals: Vec<f64> = (0..n).map(|i| s.prim(i, 0, 0)[comp]).collect();
            let half = n / 2;
            let (ip, _) = (half..n - 1).map(|i| (i, vals[i].abs())).fold((half, 0.0), |a, b| if b.1 > a.1 { b } else { a });
            let (ym, y0, yp) = (vals[ip - 1].abs(), vals[ip].abs(), vals[ip + 1].abs());
            let off = 0.5 * (ym - yp) / (ym - 2.0 * y0 + yp);
            let xp = s.mesh.x(ip, 0, 0)[0] + off * s.mesh.dx[0];
            let speed = xp / t_end;
            println!("  {name:12} N={n:5}: speed {speed:.6} vs {exact:.6} ({:+.2e}), peak amplitude {:.3e} (ε/2 = {:.3e})", speed / exact - 1.0, y0, 0.5 * eps);
        }
    }
}

fn stokes() {
    let eos = Nasg::ideal(1.4);
    let (rho, p, v0, t_end): (f64, f64, f64, f64) = (1.0, 1.0 / 1.4, 0.1, 1.0);
    let mu = 1e-2;
    let nu = mu / rho;
    let exact = |x: f64| v0 * erf(x / (2.0 * (nu * t_end).sqrt()));
    let n: usize = std::env::var("N").ok().map(|v| v.parse().unwrap()).unwrap_or(400);
    let err = |vals: &[f64], xs: &[f64]| -> f64 { vals.iter().zip(xs).map(|(v, &x)| (v - exact(x)).abs()).sum::<f64>() / vals.len() as f64 };
    let css: Vec<f64> = std::env::var("CS").ok().map(|v| v.split(',').map(|x| x.parse().unwrap()).collect()).unwrap_or(vec![1.0, 2.0, 4.0, 8.0]);
    for cs in css {
        let tau = 6.0 * mu / (rho * cs * cs);
        let sys = Gpr { eos, cs, tau, tau_model: compressible::gpr::TauModel::Constant };
        let mesh = Mesh::new([n, 1, 1], [-1.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let mut s = Solver::new(sys, mesh, |x| prim(rho, [0.0, if x[0] < 0.0 { -v0 } else { v0 }, 0.0], p));
        s.bc[0] = [Bc::Outflow; 2];
        s.cfl = 0.4;
        s.relax = Some(Box::new(|sys: &Gpr, u: &mut [f64; NG], dt: f64| sys.relax(u, dt)));
        s.strang = std::env::var("LIE").is_err();
        s.run_to(t_end);
        let xs: Vec<f64> = (0..n).map(|i| s.mesh.x(i, 0, 0)[0]).collect();
        let vals: Vec<f64> = (0..n).map(|i| s.prim(i, 0, 0)[2]).collect();
        println!(
            "  GPR c_s = {cs:4}: τ = {tau:.2e}, De = τ/t = {:.1e}, N = {n}: L1(v_y) vs Navier–Stokes erf {:.3e}  (shock-path faces {:.1}%)",
            tau / t_end,
            err(&vals, &xs),
            100.0 * s.stats.weno_faces as f64 / s.stats.faces as f64
        );
    }
    // the coarse rung itself: compressible Navier–Stokes
    let mesh = Mesh::new([n, 1, 1], [-1.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let mut s = Solver::new(Euler { eos }, mesh, |x| [rho, 0.0, if x[0] < 0.0 { -v0 } else { v0 }, 0.0, p]);
    s.bc[0] = [Bc::Outflow; 2];
    s.viscous = Some(Transport { mu, power: 0.0, t_ref: 1.0, prandtl: 1e30, bulk: 0.0 });
    s.cfl = 0.4;
    s.run_to(t_end);
    let xs: Vec<f64> = (0..n).map(|i| s.mesh.x(i, 0, 0)[0]).collect();
    let vals: Vec<f64> = (0..n).map(|i| s.prim(i, 0, 0)[2]).collect();
    println!("  Navier–Stokes solver (same grid):      L1(v_y) vs erf {:.3e}", err(&vals, &xs));
}

/// Symmetric elastic–plastic impact. Linear-elastic theory in uniaxial
/// strain: the elastic precursor runs at c_L = √(c₀² + 4c_s²/3) and carries
/// the Hugoniot elastic limit σ_HEL = (K + 4G/3) σ_Y / (2G) (von Mises
/// yield at shear strain σ_Y/2G), i.e. a velocity jump u_HEL = σ_HEL/(ρ c_L);
/// behind it the plastic wave moves at about the bulk speed c₀.
fn impact() {
    let (rho, c0, cs, sy, u) = (1.0f64, 1.0f64, 1.0f64, 0.01f64, 0.05f64);
    let gamma = 4.0;
    let eos = Nasg::stiffened(gamma, rho * c0 * c0 / gamma);
    let sys = Gpr { eos, cs, tau: f64::INFINITY, tau_model: compressible::gpr::TauModel::Plastic { tau0: 1e-4, sigma_y: sy, n: std::env::var("NPL").ok().map(|v| v.parse().unwrap()).unwrap_or(20.0) } };
    let (k_bulk, g_mod) = (rho * c0 * c0, rho * cs * cs);
    let c_l = (c0 * c0 + 4.0 / 3.0 * cs * cs).sqrt();
    let u_hel = (k_bulk + 4.0 / 3.0 * g_mod) * sy / (2.0 * g_mod) / (rho * c_l);
    let t_end = 0.4;
    println!("  theory (lab frame, target moving at −u): precursor speed c_L − u = {:.4}, jump u_HEL = {u_hel:.5} for ideal plasticity", c_l - u);
    for n in [800usize, 1600, 3200] {
        let mesh = Mesh::new([n, 1, 1], [-1.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let mut s = Solver::new(sys, mesh, |x| prim(rho, [if x[0] < 0.0 { u } else { -u }, 0.0, 0.0], 0.0));
        s.bc[0] = [Bc::Outflow; 2];
        s.cfl = 0.4;
        s.strang = true;
        s.relax = Some(Box::new(|sys: &Gpr, w: &mut [f64; NG], dt: f64| sys.relax(w, dt)));
        s.run_to(t_end);
        let v: Vec<f64> = (0..n).map(|i| s.prim(i, 0, 0)[1]).collect();
        let x: Vec<f64> = (0..n).map(|i| s.mesh.x(i, 0, 0)[0]).collect();
        if std::env::var("PROFILE").is_ok() {
            for i in (n / 2..n).step_by(n / 80) {
                let w = s.prim(i, 0, 0);
                let st = s.sys.stress(&w);
                println!("    x {:.3} v {:+.5} p {:+.5} σ11 {:+.5} σ22 {:+.5} seq {:.5}", x[i], v[i], w[13], st[0][0], st[1][1], (st[0][0] - st[1][1]).abs());
            }
            return;
        }
        // right half: material initially at −u; precursor raises v to −u + u_HEL,
        // the plastic wave brings it to 0. Fronts = half-way crossings.
        let level1 = -u + u_hel;
        let cross = |target: f64| -> f64 {
            for i in n / 2..n - 1 {
                if (v[i] - target) * (v[i + 1] - target) <= 0.0 && v[i] != v[i + 1] {
                    return x[i] + (target - v[i]) / (v[i + 1] - v[i]) * (x[i + 1] - x[i]);
                }
            }
            f64::NAN
        };
        let x_plastic = cross(0.5 * level1);
        let x_elastic = cross(-u + 0.5 * u_hel);
        // plateau between the fronts
        let idx: Vec<usize> = (0..n).filter(|&i| x[i] > x_plastic + 0.05 && x[i] < x_elastic - 0.05).collect();
        let pv = idx.iter().map(|&i| v[i]).sum::<f64>() / idx.len().max(1) as f64;
        // longitudinal stress and von Mises stress actually reached in the plateau
        let (mut sl, mut se) = (0.0, 0.0);
        for &i in &idx {
            let w = s.prim(i, 0, 0);
            let st = s.sys.stress(&w);
            sl += (-w[13] + st[0][0]) / idx.len() as f64;
            se += (st[0][0] - st[1][1]).abs() / idx.len() as f64;
        }
        let sl_elastic = (k_bulk + 4.0 / 3.0 * g_mod) * se / (2.0 * g_mod);
        println!(
            "  N={n:5}: precursor speed {:.4}, jump {:.5} ({:+.1}% vs ideal), plateau σ_eq/σ_Y {:.3}, |σ_xx| {:.5} vs elastic relation {:.5} ({:+.2}%), plastic front {:.4}",
            x_elastic / t_end,
            pv + u,
            100.0 * ((pv + u) / u_hel - 1.0),
            se / sy,
            -sl,
            sl_elastic,
            100.0 * (-sl / sl_elastic - 1.0),
            x_plastic / t_end,
        );
    }
}

fn main() {
    match std::env::args().nth(1).as_deref() {
        Some("waves") => waves(),
        Some("stokes") => stokes(),
        Some("impact") => impact(),
        _ => {
            waves();
            stokes();
        }
    }
}
