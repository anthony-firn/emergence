//! Adaptive grids on the GPU vs the CPU `Amr`.
//!
//!   gpu_amr_validate riemann2d LEVELS T     2D four-shock Riemann problem (Euler, outflow)
//!   gpu_amr_validate bubble2d LEVELS T      cylindrical bubble collapse (five-equation,
//!                                           reflecting axes, characteristic far field)
//!   gpu_amr_validate bubble3d LEVELS T      spherical bubble collapse (timing; T in t_c)
//!
//! Same initial hierarchy on both; regridding runs on the CPU for both. In
//! f64 the GPU must reproduce the CPU to round-off (identical hierarchies
//! are then compared leaf by leaf); f32 is compared with the f64 result.
//! Env CPU=0 skips the CPU run, F64=0 the GPU f64 run.
use compressible::amr::{Amr, AmrConfig};
use compressible::eos::materials::{AIR, WATER_SG};
use compressible::eos::Nasg;
use compressible::multiphase::FiveEq;
use compressible::{Bc, Euler, System};
use compressible_cuda::amr::GpuAmr;
use compressible_cuda::{GpuSystem, Real};
use std::time::Instant;

fn structure<const N: usize, S: System<N>>(a: &Amr<N, S>) -> Vec<(usize, [i64; 3], bool)> {
    a.blocks.iter().map(|b| (b.level, b.idx, b.leaf)).collect()
}

/// max |Δ| (absolute, relative to max |·|) of primitive `comp` over leaves
fn diff<const N: usize, S: System<N>>(a: &Amr<N, S>, b: &Amr<N, S>, comp: usize) -> Option<(f64, f64)> {
    if structure(a) != structure(b) {
        return None;
    }
    let (mut e, mut sc) = (0.0f64, 0.0f64);
    for (x, y) in a.blocks.iter().zip(b.blocks.iter()).filter(|(x, _)| x.leaf) {
        for (u, v) in x.u.iter().zip(y.u.iter()) {
            let (p, q) = (a.sys.to_prim(u)[comp], a.sys.to_prim(v)[comp]);
            e = e.max((p - q).abs());
            sc = sc.max(p.abs());
        }
    }
    Some((e, e / sc))
}

fn run_gpu<R: Real, const N: usize, S: System<N>>(a: Amr<N, S>, sys: GpuSystem, t_end: f64, label: &str) -> Amr<N, S> {
    let t0 = Instant::now();
    let mut g = GpuAmr::<R, N, S>::new(a, sys).expect("gpu init");
    g.profile = std::env::var("PROFILE").is_ok();
    let n = g.run_to(t_end).expect("gpu run");
    g.download().unwrap();
    let el = t0.elapsed().as_secs_f64();
    println!(
        "  {label:<12}: {n} steps, {el:.2}s, {} leaves at end, counters {:?}{}",
        g.leaf_count(),
        g.counters().unwrap(),
        if g.profile {
            format!(
                ", profile [restrict, ghosts, fluxes, reflux, regrid] {:?}, regrid [details, download, refine, upload, coarsen] {:?}",
                g.prof.map(|v| (v * 100.0).round() / 100.0),
                g.rprof.map(|v| (v * 100.0).round() / 100.0)
            )
        } else {
            String::new()
        }
    );
    g.cpu
}

fn compare_all<const N: usize, S: System<N>>(make: impl Fn() -> Amr<N, S>, sys: GpuSystem, t_end: f64, comp: usize, name: &str) {
    let base = make();
    println!("{name}: {} leaves initially, {} cells ({:.1}% of uniform)", base.blocks.iter().filter(|b| b.leaf).count(), base.cells(), 100.0 * base.cells() as f64 / base.uniform_cells() as f64);
    let want = |k: &str| std::env::var(k).map(|v| v != "0").unwrap_or(true);
    let cpu = want("CPU").then(|| {
        let mut c = make();
        let t0 = Instant::now();
        let n = c.run_to(t_end);
        println!("  CPU f64     : {n} steps, {:.2}s, {} leaves at end", t0.elapsed().as_secs_f64(), c.blocks.iter().filter(|b| b.leaf).count());
        c
    });
    let g64 = want("F64").then(|| run_gpu::<f64, N, S>(make(), sys, t_end, "GPU f64"));
    let g32 = run_gpu::<f32, N, S>(make(), sys, t_end, "GPU f32");
    let show = |a: &Amr<N, S>, b: &Amr<N, S>, what: &str| match diff(a, b, comp) {
        Some((e, r)) => println!("  {what}: max |Δ| {e:.3e} (relative {r:.3e})"),
        None => println!("  {what}: hierarchies differ (different regrid decisions)"),
    };
    if let (Some(c), Some(g)) = (&cpu, &g64) {
        show(c, g, "GPU f64 vs CPU");
    }
    if let Some(g) = &g64 {
        show(g, &g32, "GPU f32 vs GPU f64");
    } else if let Some(c) = &cpu {
        show(c, &g32, "GPU f32 vs CPU");
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let case = a.get(1).cloned().unwrap_or("riemann2d".into());
    let levels: usize = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(3);
    let tt: Option<f64> = a.get(3).map(|v| v.parse().unwrap());
    match case.as_str() {
        "riemann2d" => {
            let eos = Nasg::ideal(1.4);
            let make = || {
                let cfg = AmrConfig { b: 16, levels, eps: 1e-3, regrid_every: 4, base_blocks: [4, 4, 0], ghosts: 5 };
                let mut s = Amr::new(Euler { eos }, cfg, [0.0; 3], [1.0, 1.0, 0.0], |x| {
                    // configuration 3 (Lax & Liu)
                    let (r, u, v, p) = match (x[0] > 0.5, x[1] > 0.5) {
                        (true, true) => (1.5, 0.0, 0.0, 1.5),
                        (false, true) => (0.5323, 1.206, 0.0, 0.3),
                        (false, false) => (0.138, 1.206, 1.206, 0.029),
                        (true, false) => (0.5323, 0.0, 1.206, 0.3),
                    };
                    [r, u, v, 0.0, p]
                });
                for d in 0..2 {
                    s.bc[d] = [Bc::Outflow; 2];
                }
                s
            };
            compare_all(make, GpuSystem::Euler(eos), tt.unwrap_or(0.3), 0, "2D Riemann config 3 (Euler)");
        }
        "bubble2d" | "bubble3d" => {
            let d3 = case == "bubble3d";
            let (r0, pg0, ratio, lr) = (1e-3, 1e5, 100.0, 5.0);
            let pinf = pg0 * ratio;
            let (rho_w, rho_g) = (1000.0, pg0 / (287.0 * 293.0));
            let len = lr * r0;
            let sys = FiveEq { eos: [WATER_SG, AIR], kapila: true };
            let (b, nb0) = (16, 4);
            let n_fine = nb0 * b << (levels - 1);
            let dx_f = len / n_fine as f64;
            let eps = if d3 { 1e-2 } else { 1e-3 };
            let make = move || {
                let init = move |x: [f64; 3]| -> [f64; 7] {
                    let r = (x[0] * x[0] + x[1] * x[1] + if d3 { x[2] * x[2] } else { 0.0 }).sqrt();
                    let e = 1e-6;
                    let alpha = e + (1.0 - 2.0 * e) * 0.5 * (1.0 + ((r - r0) / (0.75 * dx_f)).tanh());
                    // 3D: p∞ + R0/r (p_g − p∞); 2D (cylinder): logarithmic profile out to the box
                    let p = if r <= r0 {
                        pg0
                    } else if d3 {
                        pinf + r0 / r * (pg0 - pinf)
                    } else {
                        pinf + (pg0 - pinf) * (1.0 - (r / r0).ln() / (len * 1.5 / r0).ln()).max(0.0)
                    };
                    [alpha * rho_w, (1.0 - alpha) * rho_g, 0.0, 0.0, 0.0, p, alpha]
                };
                let cfg = AmrConfig { b, levels, eps, regrid_every: 4, base_blocks: [nb0, nb0, if d3 { nb0 } else { 0 }], ghosts: 5 };
                let mut s = Amr::new(sys, cfg, [0.0; 3], [len, len, if d3 { len } else { 0.0 }], init);
                let far = {
                    let mut w = [0.0; 7];
                    w[0] = (1.0 - 1e-6) * rho_w;
                    w[1] = 1e-6 * rho_g;
                    w[5] = pinf;
                    w[6] = 1.0 - 1e-6;
                    s.sys.to_cons(&w)
                };
                for d in 0..3 {
                    s.bc[d] = [Bc::Reflect, Bc::Characteristic(far)];
                }
                s
            };
            let tc = 0.915 * r0 * (rho_w / (pinf - pg0)).sqrt();
            let name = format!("{} bubble collapse, effective {n_fine}^{}, ε {eps:.0e}", if d3 { "spherical" } else { "cylindrical" }, if d3 { 3 } else { 2 });
            compare_all(make, GpuSystem::FiveEq { eos: [WATER_SG, AIR], kapila: true }, tt.unwrap_or(if d3 { 0.05 } else { 0.5 }) * tc, 6, &name);
        }
        _ => panic!("unknown case"),
    }
}
