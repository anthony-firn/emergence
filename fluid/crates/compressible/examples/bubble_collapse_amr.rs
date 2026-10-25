//! Rayleigh collapse of a gas bubble on the adaptive multiresolution grid
//! (octant, symmetry planes) vs Keller–Miksis / Rayleigh–Plesset.
//!
//!   bubble_collapse_amr LEVELS [eps] [p_ratio]
//!
//! Base grid 64³ (4³ blocks of 16³); LEVELS = 2 → effective 128³,
//! 3 → 256³. The interface and the collapse region stay at the finest
//! level; the far field coarsens.
use compressible::amr::{Amr, AmrConfig};
use compressible::eos::materials::{AIR, WATER_SG};
use compressible::multiphase::FiveEq;
use compressible::Bc;
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
    let levels: usize = a.get(1).map(|v| v.parse().unwrap()).unwrap_or(3);
    let eps: f64 = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(1e-3);
    let ratio: f64 = a.get(3).map(|v| v.parse().unwrap()).unwrap_or(100.0);
    let (r0, pg0, lr) = (1e-3, 1e5, 5.0);
    let pinf = pg0 * ratio;
    let (rho_w, rho_g) = (1000.0, pg0 / (287.0 * 293.0));
    let len = lr * r0;
    let sys = FiveEq { eos: [WATER_SG, AIR], kapila: true };
    let b = 16;
    let nb0 = 4;
    let n_fine = nb0 * b << (levels - 1);
    let dx_f = len / n_fine as f64;
    let init = move |x: [f64; 3]| -> [f64; 7] {
        let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        let e = 1e-6;
        let alpha = e + (1.0 - 2.0 * e) * 0.5 * (1.0 + ((r - r0) / (0.75 * dx_f)).tanh());
        let p = if r > r0 { pinf + r0 / r * (pg0 - pinf) } else { pg0 };
        let mut w = [0.0; 7];
        w[0] = alpha * rho_w;
        w[1] = (1.0 - alpha) * rho_g;
        w[5] = p;
        w[6] = alpha;
        w
    };
    let cfg = AmrConfig { b, levels, eps, regrid_every: 4, base_blocks: [nb0, nb0, nb0], ghosts: 5 };
    let t0 = Instant::now();
    let mut s = Amr::new(sys, cfg, [0.0; 3], [len; 3], init);
    let far = {
        let mut w = [0.0; 7];
        w[0] = (1.0 - 1e-6) * rho_w;
        w[1] = 1e-6 * rho_g;
        w[5] = pinf;
        w[6] = 1.0 - 1e-6;
        compressible::System::to_cons(&s.sys, &w)
    };
    for d in 0..3 {
        s.bc[d] = [Bc::Reflect, Bc::Characteristic(far)];
    }
    s.subcycle = std::env::var("SUBCYCLE").is_ok();
    let tc = 0.915 * r0 * (rho_w / (pinf - pg0)).sqrt();
    let t_end = std::env::var("TEND").ok().map(|v| v.parse::<f64>().unwrap()).unwrap_or(1.3) * tc;
    let c_w = (WATER_SG.gamma * (pinf + WATER_SG.p_inf) / rho_w).sqrt();
    let km_c = km(r0, pg0, pinf, rho_w, c_w, AIR.gamma, t_end, true);
    let rp = km(r0, pg0, pinf, rho_w, c_w, AIR.gamma, t_end, false);
    let rmin = |v: &[(f64, f64)]| v.iter().cloned().fold((0.0, f64::INFINITY), |a, b| if b.1 < a.1 { b } else { a });
    let (tk, rk) = rmin(&km_c);
    println!(
        "p∞/p_g0 = {ratio}, effective {n_fine}³ (R0 = {:.1} finest cells), ε = {eps:.0e}; Keller–Miksis R_min/R0 {:.4} at t/t_c {:.4}; init {:.0}s, {} cells ({:.1}% of uniform)",
        r0 / dx_f,
        rk / r0,
        tk / tc,
        t0.elapsed().as_secs_f64(),
        s.cells(),
        100.0 * s.cells() as f64 / s.uniform_cells() as f64
    );
    let radius = |s: &Amr<7, FiveEq>| -> f64 {
        let mut v = 0.0;
        s.for_each_leaf_cell(|_, _, u, dv| v += (1.0 - u[6]) * dv);
        (3.0 * 8.0 * v / (4.0 * std::f64::consts::PI)).cbrt()
    };
    let interp = |v: &[(f64, f64)], t: f64| {
        let i = v.partition_point(|p| p.0 < t).min(v.len() - 1).max(1);
        let (a, b) = (v[i - 1], v[i]);
        a.1 + (b.1 - a.1) * (t - a.0) / (b.0 - a.0)
    };
    let t1 = Instant::now();
    let mut hist = vec![(0.0, radius(&s))];
    let mut next = 0.05 * tc;
    let mut cell_steps = 0.0;
    let mut steps = 0.0;
    println!("  t/t_c   R/R0 sim   Keller–Miksis   Rayleigh–Plesset   cells (% uniform)   levels");
    while s.t < t_end {
        let dt = if s.subcycle { s.dt_cfl_subcycled() } else { s.dt_cfl() }.min(next - s.t);
        if s.subcycle {
            s.step_subcycled(dt);
        } else {
            s.step(dt);
        }
        cell_steps += s.cells() as f64;
        steps += 1.0;
        if s.t >= next - 1e-15 {
            let r = radius(&s);
            hist.push((s.t, r));
            println!(
                "  {:.3}   {:.4}     {:.4}          {:.4}            {:>9} ({:4.1}%)   {:?}",
                s.t / tc,
                r / r0,
                interp(&km_c, s.t) / r0,
                interp(&rp, s.t) / r0,
                s.cells(),
                100.0 * s.cells() as f64 / s.uniform_cells() as f64,
                s.level_counts()
            );
            next += 0.05 * tc;
        }
    }
    let (ts, rs) = rmin(&hist);
    println!(
        "  simulation: R_min/R0 {:.4} at t/t_c {:.3}; mean cells {:.1}% of uniform; {:.0} steps, {:.0}s; profile {:?}",
        rs / r0,
        ts / tc,
        100.0 * cell_steps / steps / s.uniform_cells() as f64,
        steps,
        t1.elapsed().as_secs_f64(),
        s.prof.map(|v| v.round())
    );
}
