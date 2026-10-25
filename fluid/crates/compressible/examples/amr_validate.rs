//! Adaptive multiresolution vs the uniform finest grid.
//!
//!   amr_validate shuosher   1D shock / entropy-wave interaction
//!   amr_validate vortex     2D isentropic vortex (smooth, moving)
//!   amr_validate riemann2d  2D four-shock Riemann problem
//!
//! For each threshold ε: the difference between the adaptive solution and
//! the uniform grid at the finest level (the reference scheme), which the
//! multiresolution analysis bounds by O(ε), next to the discretization
//! error of that uniform grid, the cells used, and the run times.
use compressible::amr::{Amr, AmrConfig};
use compressible::eos::Nasg;
use compressible::{Bc, Euler, FluxMode, Mesh, Solver};
use std::time::Instant;

type V = [f64; 5];

fn prim(rho: f64, u: f64, v: f64, p: f64) -> V {
    [rho, u, v, 0.0, p]
}

struct Case {
    name: &'static str,
    lo: [f64; 3],
    hi: [f64; 3],
    base: [usize; 3],
    levels: usize,
    t_end: f64,
    bc: Bc<5>,
    init: fn([f64; 3]) -> V,
}

fn shu_osher(x: [f64; 3]) -> V {
    if x[0] < -4.0 {
        prim(3.857143, 2.629369, 0.0, 10.33333)
    } else {
        prim(1.0 + 0.2 * (5.0 * x[0]).sin(), 0.0, 0.0, 1.0)
    }
}

fn vortex(x: [f64; 3]) -> V {
    let (g, beta) = (1.4f64, 5.0f64);
    let pi = std::f64::consts::PI;
    let r2 = x[0] * x[0] + x[1] * x[1];
    let f = beta / (2.0 * pi) * (0.5 * (1.0 - r2)).exp();
    let tt = 1.0 - (g - 1.0) * beta * beta / (8.0 * g * pi * pi) * (1.0 - r2).exp();
    let rho = tt.powf(1.0 / (g - 1.0));
    prim(rho, 1.0 - f * x[1], 1.0 + f * x[0], rho * tt)
}

/// Lax–Liu configuration 3 (four shocks).
fn riemann2d(x: [f64; 3]) -> V {
    match (x[0] > 0.8, x[1] > 0.8) {
        (true, true) => prim(1.5, 0.0, 0.0, 1.5),
        (false, true) => prim(0.5323, 1.206, 0.0, 0.3),
        (false, false) => prim(0.138, 1.206, 1.206, 0.029),
        (true, false) => prim(0.5323, 0.0, 1.206, 0.3),
    }
}

fn main() {
    let which = std::env::args().nth(1).unwrap_or("shuosher".into());
    let case = match which.as_str() {
        "vortex" => Case { name: "isentropic vortex", lo: [-5.0, -5.0, 0.0], hi: [5.0, 5.0, 1.0], base: [2, 2, 0], levels: 4, t_end: 2.0, bc: Bc::Periodic, init: vortex },
        "riemann2d" => Case { name: "2D Riemann (config 3)", lo: [0.0, 0.0, 0.0], hi: [1.0, 1.0, 1.0], base: [2, 2, 0], levels: 4, t_end: 0.3, bc: Bc::Outflow, init: riemann2d },
        _ => Case { name: "Shu–Osher", lo: [-5.0, 0.0, 0.0], hi: [5.0, 1.0, 1.0], base: [2, 0, 0], levels: 6, t_end: 1.8, bc: Bc::Outflow, init: shu_osher },
    };
    // optional: fewer levels with the same finest grid (more base blocks)
    let mut case = case;
    if let Some(lv) = std::env::var("LEVELS").ok().and_then(|v| v.parse::<usize>().ok()) {
        let f = 1usize << (case.levels - lv);
        case.base = case.base.map(|v| v * f);
        case.levels = lv;
    }
    let eos = Nasg::ideal(1.4);
    let b = 16;
    let nd = case.base.iter().filter(|&&v| v > 0).count();
    let nf: [usize; 3] = std::array::from_fn(|d| if case.base[d] > 0 { case.base[d] * b << (case.levels - 1) } else { 1 });
    let hi: [f64; 3] = std::array::from_fn(|d| if case.base[d] > 0 { case.hi[d] } else { case.lo[d] + 1.0 });
    // uniform finest grid (reference scheme)
    let t0 = Instant::now();
    let mesh = Mesh::new(nf, case.lo, hi);
    let mode = match std::env::var("MODE").as_deref() {
        Ok("weno") => FluxMode::Weno,
        _ => FluxMode::Hybrid,
    };
    let mut uni = Solver::new(Euler { eos }, mesh, case.init);
    uni.mode = mode;
    for d in 0..3 {
        uni.bc[d] = [case.bc; 2];
    }
    uni.cfl = 0.4;
    uni.run_to(case.t_end);
    let t_uni = t0.elapsed().as_secs_f64();
    let n_uni: usize = nf.iter().product();
    // twice finer uniform grid: discretization error of the reference
    let nf2: [usize; 3] = std::array::from_fn(|d| if case.base[d] > 0 { 2 * nf[d] } else { 1 });
    let mut fine = Solver::new(Euler { eos }, Mesh::new(nf2, case.lo, hi), case.init);
    fine.mode = mode;
    for d in 0..3 {
        fine.bc[d] = [case.bc; 2];
    }
    fine.cfl = 0.4;
    fine.run_to(case.t_end);
    // compare uniform vs fine (average 2^D fine cells onto each uniform cell)
    let avg_fine = |i: usize, j: usize, k: usize| -> f64 {
        let mut s = 0.0;
        let mut c = 0.0;
        for a in 0..if nf2[0] > 1 { 2 } else { 1 } {
            for bb in 0..if nf2[1] > 1 { 2 } else { 1 } {
                for cc in 0..if nf2[2] > 1 { 2 } else { 1 } {
                    let ii = if nf2[0] > 1 { 2 * i + a } else { 0 };
                    let jj = if nf2[1] > 1 { 2 * j + bb } else { 0 };
                    let kk = if nf2[2] > 1 { 2 * k + cc } else { 0 };
                    s += fine.prim(ii, jj, kk)[0];
                    c += 1.0;
                }
            }
        }
        s / c
    };
    let mut disc = 0.0;
    for i in 0..nf[0] {
        for j in 0..nf[1] {
            for k in 0..nf[2] {
                disc += (uni.prim(i, j, k)[0] - avg_fine(i, j, k)).abs();
            }
        }
    }
    disc /= n_uni as f64;
    println!("{}: {nd}D, finest {:?} = {n_uni} cells, uniform run {t_uni:.1}s; discretization error of that grid L1(ρ) ≈ {disc:.3e}", case.name, &nf[..nd]);
    let epss: Vec<f64> = std::env::var("EPS").ok().map(|v| v.split(",").map(|x| x.parse().unwrap()).collect()).unwrap_or(vec![1e-2, 1e-3, 1e-4, 1e-5]);
    for eps in epss {
        let cfg = AmrConfig { b, levels: case.levels, eps, regrid_every: 4, base_blocks: case.base, ghosts: 5 };
        let t0 = Instant::now();
        let mut a = Amr::new(Euler { eos }, cfg, case.lo, hi, case.init);
        for d in 0..3 {
            a.bc[d] = [case.bc; 2];
        }
        a.mode = mode;
        a.subcycle = std::env::var("SUBCYCLE").is_ok();
        let mass0 = a.totals()[0];
        let mut cell_steps = 0usize;
        let mut max_cells = 0usize;
        let mut steps = 0usize;
        while a.t < case.t_end * (1.0 - 1e-14) {
            let dt = if a.subcycle { a.dt_cfl_subcycled() } else { a.dt_cfl() }.min(case.t_end - a.t);
            if a.subcycle {
                a.step_subcycled(dt);
            } else {
                a.step(dt);
            }
            steps += 1;
            cell_steps += a.cells();
            max_cells = max_cells.max(a.cells());
        }
        let t_amr = t0.elapsed().as_secs_f64();
        println!("    profile [restrict, ghosts, fluxes, reflux, regrid] s: {:?}; mass drift {:.1e} ({} steps)", a.prof.map(|v| (v * 10.0).round() / 10.0), a.totals()[0] / mass0 - 1.0, steps);
        let l = case.levels - 1;
        let mut diff = 0.0;
        let mut dmax: f64 = 0.0;
        for i in 0..nf[0] {
            for j in 0..nf[1] {
                for k in 0..nf[2] {
                    let g = [i as i64, j as i64, k as i64];
                    let r = compressible::System::to_prim(&a.sys, &a.sample(l, g))[0];
                    let e = (r - uni.prim(i, j, k)[0]).abs();
                    diff += e;
                    dmax = dmax.max(e);
                }
            }
        }
        diff /= n_uni as f64;
        if std::env::var("WHERE").is_ok() {
            let mut worst = (0.0, 0usize);
            for i in 0..nf[0] {
                let g = [i as i64, 0, 0];
                let r = compressible::System::to_prim(&a.sys, &a.sample(l, g))[0];
                let e = (r - uni.prim(i, 0, 0)[0]).abs();
                if e > worst.0 { worst = (e, i); }
            }
            println!("    worst at x = {:.3}", uni.mesh.x(worst.1, 0, 0)[0]);
        }
        println!(
            "  ε = {eps:.0e}: |AMR − uniform| L1(ρ) {diff:.3e} (max {dmax:.2e}); mean cells {:.0} ({:.1}% of uniform), peak {:.1}%, levels {:?}, {:.1}s ({:.1}× faster)",
            cell_steps as f64 / steps as f64,
            100.0 * cell_steps as f64 / steps as f64 / n_uni as f64,
            100.0 * max_cells as f64 / n_uni as f64,
            a.level_counts(),
            t_amr,
            t_uni / t_amr
        );
    }
}
