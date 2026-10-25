//! Rayleigh collapse of a gas bubble in water: 3D five-equation simulation
//! (one octant, symmetry planes) against the Keller–Miksis and
//! Rayleigh–Plesset equations.
//!
//!   bubble_collapse N [p_ratio] [L/R0]
//!
//! Keller–Miksis (liquid compressibility to first order in Ṙ/c):
//!   (1 − Ṙ/c) R R̈ + (3/2)(1 − Ṙ/3c) Ṙ² = (1 + Ṙ/c)(p_B − p∞)/ρ + R ṗ_B/(ρc),
//!   p_B = p_g0 (R0/R)^{3γ} (adiabatic gas, no surface tension/viscosity).
//! Initial liquid pressure p∞ + (R0/r)(p_g0 − p∞) (the potential-flow field
//! at rest), as in Johnsen & Colonius (2009).
use compressible::eos::materials::{AIR, WATER_SG};
use compressible::multiphase::FiveEq;
type V = [f64; 7];
const MAXV: usize = 7;
use compressible::{Bc, Mesh, Solver};
use std::time::Instant;

fn km(r0: f64, pg0: f64, pinf: f64, rho: f64, c: f64, gamma: f64, t_end: f64, compressible: bool) -> Vec<(f64, f64)> {
    let ci = if compressible { 1.0 / c } else { 0.0 };
    let rhs = |r: f64, v: f64| -> f64 {
        let pb = pg0 * (r0 / r).powf(3.0 * gamma);
        let dpb = -3.0 * gamma * pb * v / r;
        ((1.0 + v * ci) * (pb - pinf) / rho + r * dpb * ci / rho - 1.5 * (1.0 - v * ci / 3.0) * v * v) / ((1.0 - v * ci) * r)
    };
    let (mut r, mut v, mut t) = (r0, 0.0, 0.0);
    let dt = t_end / 200_000.0;
    let mut out = vec![(0.0, r0)];
    while t < t_end {
        let (k1r, k1v) = (v, rhs(r, v));
        let (k2r, k2v) = (v + 0.5 * dt * k1v, rhs(r + 0.5 * dt * k1r, v + 0.5 * dt * k1v));
        let (k3r, k3v) = (v + 0.5 * dt * k2v, rhs(r + 0.5 * dt * k2r, v + 0.5 * dt * k2v));
        let (k4r, k4v) = (v + dt * k3v, rhs(r + dt * k3r, v + dt * k3v));
        r += dt / 6.0 * (k1r + 2.0 * k2r + 2.0 * k3r + k4r);
        v += dt / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
        t += dt;
        out.push((t, r));
    }
    out
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let n: usize = a.get(1).map(|v| v.parse().unwrap()).unwrap_or(48);
    let ratio: f64 = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(100.0);
    let lr: f64 = a.get(3).map(|v| v.parse().unwrap()).unwrap_or(5.0);
    let r0 = 1e-3;
    let pg0 = 1e5;
    let pinf = pg0 * ratio;
    let (rho_w, rho_g) = (1000.0, pg0 / (287.0 * 293.0)); // air at 293 K
    let len = lr * r0;
    let eps = 1e-6;
    // KAPILA=0: Allaire five-equation model (no K∇·u term)
    let sys = FiveEq { eos: [WATER_SG, AIR], kapila: std::env::var("KAPILA").map(|v| v != "0").unwrap_or(true) };
    let mesh = Mesh::new([n, n, n], [0.0; 3], [len, len, len]);
    let dx = mesh.dx[0];
    let init = |x: [f64; 3]| -> V {
        let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        // water fraction: smoothed step over ~1.5 cells
        let alpha = eps + (1.0 - 2.0 * eps) * 0.5 * (1.0 + ((r - r0) / (0.75 * dx)).tanh());
        let p = if r > r0 { pinf + r0 / r * (pg0 - pinf) } else { pg0 };
        let mut w = [0.0; MAXV];
        w[0] = alpha * rho_w;
        w[1] = (1.0 - alpha) * rho_g;
        w[5] = p;
        w[6] = alpha;
        w
    };
    let mut s = Solver::new(sys, mesh, init);
    let far = {
        let mut w = [0.0; MAXV];
        w[0] = (1.0 - eps) * rho_w;
        w[1] = eps * rho_g;
        w[5] = pinf;
        w[6] = 1.0 - eps;
        compressible::System::to_cons(&s.sys, &w)
    };
    let outer = if std::env::var("BC").as_deref() == Ok("zero") { Bc::Outflow } else { Bc::Characteristic(far) };
    for d in 0..3 {
        s.bc[d] = [Bc::Reflect, outer];
    }
    s.cfl = 0.4;
    if let Ok(v) = std::env::var("ALLSPEED") {
        // Mach-adaptive all-speed stepping; ALLSPEED=<max acoustic Courant>
        let cap = v.parse::<f64>().unwrap_or(f64::INFINITY);
        s.allspeed = Some(compressible::allspeed::AllSpeedOpts { cfl: 0.4, max_acoustic_cfl: cap, ..Default::default() });
    }
    let tc = 0.915 * r0 * (rho_w / (pinf - pg0)).sqrt();
    let t_end = 1.3 * tc;
    let c_w = (WATER_SG.gamma * (pinf + WATER_SG.p_inf) / rho_w).sqrt();
    let km_c = km(r0, pg0, pinf, rho_w, c_w, AIR.gamma, t_end, true);
    let rp = km(r0, pg0, pinf, rho_w, c_w, AIR.gamma, t_end, false);
    let rmin = |v: &[(f64, f64)]| v.iter().cloned().fold((0.0, f64::INFINITY), |a, b| if b.1 < a.1 { b } else { a });
    let (tk, rk) = rmin(&km_c);
    let (tr, rr) = rmin(&rp);
    println!("p∞/p_g0 = {ratio}, R0 = {} cells, box {lr} R0: Rayleigh time {tc:.4e} s", r0 / dx);
    println!("  Keller–Miksis: R_min/R0 {:.4} at t/t_c {:.4}; Rayleigh–Plesset: {:.4} at {:.4}", rk / r0, tk / tc, rr / r0, tr / tc);
    let vol = |s: &Solver<7, FiveEq>| -> f64 {
        let mut v = 0.0;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    v += 1.0 - s.prim(i, j, k)[6];
                }
            }
        }
        8.0 * v * dx * dx * dx
    };
    let radius = |v: f64| (3.0 * v / (4.0 * std::f64::consts::PI)).cbrt();
    let t0 = Instant::now();
    let mut hist = vec![(0.0, radius(vol(&s)))];
    let mut next = 0.05 * tc;
    let mut steps = 0usize;
    while s.t < t_end {
        steps += 1;
        let dt = s.dt_cfl().min(next - s.t);
        s.step(dt);
        if std::env::var("DEBUG").is_ok() {
            let mut pmin = (f64::INFINITY, 0usize, 0usize, 0usize);
            for i in 0..n { for j in 0..n { for k in 0..n { let p = s.prim(i,j,k)[5]; if p < pmin.0 { pmin = (p, i, j, k); } }}}
            println!("t/tc {:.4} fallback {} pmin {:.3e} at ({},{},{})", s.t / tc, s.stats.fallback_faces, pmin.0, pmin.1, pmin.2, pmin.3);
            let mut worst: Option<(usize, usize, usize, f64)> = None;
            for i in 0..n { for j in 0..n { for k in 0..n {
                let w = s.prim(i, j, k);
                let bad = !w[5].is_finite() || w[5] < 0.0 || w[6] < -0.01 || w[6] > 1.01;
                if bad { worst = Some((i, j, k, w[5])); }
            }}}
            if let Some((i, j, k, p)) = worst {
                println!("t/tc {:.4}: bad cell ({i},{j},{k}) r/R0 {:.3} p {p:e}", s.t / tc, (s.mesh.x(i,j,k).iter().map(|v| v*v).sum::<f64>()).sqrt() / r0);
                for di in 0..5usize { let ii = (i + di).saturating_sub(2).min(n-1); let w = s.prim(ii, j, k); println!("   i={ii}: a1 {:.3e} a2 {:.3e} u {:.3e} v {:.3e} p {:.3e} alpha {:.3e}", w[0], w[1], w[2], w[3], w[5], w[6]); }
                return;
            }
        }
        if s.t >= next - 1e-15 {
            hist.push((s.t, radius(vol(&s))));
            next += 0.05 * tc;
        }
    }
    let interp = |v: &[(f64, f64)], t: f64| {
        let i = v.partition_point(|p| p.0 < t).min(v.len() - 1).max(1);
        let (a, b) = (v[i - 1], v[i]);
        a.1 + (b.1 - a.1) * (t - a.0) / (b.0 - a.0)
    };
    println!("  t/t_c   R/R0 sim   Keller–Miksis   Rayleigh–Plesset");
    for &(t, r) in &hist {
        println!("  {:.3}   {:.4}     {:.4}          {:.4}", t / tc, r / r0, interp(&km_c, t) / r0, interp(&rp, t) / r0);
    }
    let (ts, rs) = rmin(&hist);
    println!(
        "  simulation: R_min/R0 {:.4} at t/t_c {:.3} (sampled every 0.05 t_c); shock-path faces {:.1}%; fallback {}; {:.0}s",
        rs / r0,
        ts / tc,
        100.0 * s.stats.weno_faces as f64 / s.stats.faces as f64,
        s.stats.fallback_faces,
        t0.elapsed().as_secs_f64()
    );
    let a = s.allspeed_stats;
    println!(
        "  {steps} steps{}",
        if s.allspeed.is_some() {
            format!(
                "; all-speed: {:.0} CG iterations/solve, {:.1}% of cell-steps in the incompressible regime (acoustic Courant > 1), max acoustic Courant {:.1}, {} steps redone explicitly ({} sub-steps)",
                a.cg_iterations as f64 / a.solves.max(1) as f64,
                100.0 * a.incompressible_cells as f64 / a.cells.max(1) as f64,
                a.max_acoustic_cfl,
                a.explicit_fallbacks,
                a.explicit_substeps
            )
        } else {
            String::new()
        }
    );
}
