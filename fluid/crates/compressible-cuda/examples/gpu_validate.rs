//! GPU vs CPU for the compressible solver.
//!
//!   gpu_validate riemann2d N     2D four-shock Riemann problem (Euler)
//!   gpu_validate tube N          water | air shock tube (five-equation)
//!   gpu_validate taylor N        3D compressible Taylor–Green, Ma 0.1 (timing)
//!
//! f64 on the GPU must reproduce the CPU solver to round-off; f32 is
//! compared against the f64 result; wall times for all three.
use compressible::eos::materials::{AIR, WATER_SG};
use compressible::eos::Nasg;
use compressible::multiphase::FiveEq;
use compressible::{Bc, Euler, Mesh, Solver, System};
use compressible_cuda::{GpuSolver, GpuSystem};
use std::time::Instant;

fn compare<const N: usize, S: System<N>>(cpu: &Solver<N, S>, gpu: &[[f64; N]], comp: usize) -> (f64, f64) {
    let m = &cpu.mesh;
    let (mut e, mut sc) = (0.0f64, 0.0f64);
    for i in 0..m.n[0] {
        for j in 0..m.n[1] {
            for k in 0..m.n[2] {
                let q = m.idx(i + m.g(0), j + m.g(1), k + m.g(2));
                let a = cpu.sys.to_prim(&cpu.u[q])[comp];
                let b = cpu.sys.to_prim(&gpu[q])[comp];
                e = e.max((a - b).abs());
                sc = sc.max(a.abs());
            }
        }
    }
    (e, e / sc)
}

fn run<const N: usize, S: System<N> + Clone>(name: &str, mk: impl Fn() -> Solver<N, S>, gsys: GpuSystem, t_end: f64, comp: usize) {
    let mut cpu = mk();
    let t0 = Instant::now();
    let gpu_only = std::env::var("GPU_ONLY").is_ok();
    let steps = if gpu_only { 0 } else { cpu.run_to(t_end) };
    let t_cpu = t0.elapsed().as_secs_f64();
    let mut g64 = GpuSolver::<f64>::from_cpu(&mk(), gsys).unwrap();
    let t0 = Instant::now();
    let s64 = if gpu_only { 0 } else { g64.run_to(t_end).unwrap() };
    g64.synchronize().unwrap();
    let t_g64 = t0.elapsed().as_secs_f64();
    let u64v = g64.download::<N>().unwrap();
    let mut g32 = GpuSolver::<f32>::from_cpu(&mk(), gsys).unwrap();
    let t0 = Instant::now();
    let s32 = g32.run_to(t_end).unwrap();
    g32.synchronize().unwrap();
    let t_g32 = t0.elapsed().as_secs_f64();
    let u32v = g32.download::<N>().unwrap();
    let (e64, r64) = compare(&cpu, &u64v, comp);
    let (e32, r32) = compare(&cpu, &u32v, comp);
    println!("{name}: {} cells", cpu.mesh.cells());
    if gpu_only {
        println!("  GPU f32: {s32} steps, {t_g32:.2}s, {:.1} M cell-steps/s, counters {:?}", (cpu.mesh.cells() * s32) as f64 / t_g32 / 1e6, g32.counters().unwrap());
        return;
    }
    println!("  CPU f64      : {steps} steps, {t_cpu:.2}s");
    println!("  GPU f64      : {s64} steps, {t_g64:.2}s ({:.1}× CPU), max |Δ| vs CPU {e64:.2e} (relative {r64:.2e}), counters {:?}", t_cpu / t_g64, g64.counters().unwrap());
    println!("  GPU f32      : {s32} steps, {t_g32:.2}s ({:.1}× CPU), max |Δ| vs CPU {e32:.2e} (relative {r32:.2e})", t_cpu / t_g32);
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let which = a.get(1).cloned().unwrap_or("riemann2d".into());
    let n: usize = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(256);
    match which.as_str() {
        "tube" => {
            let sys = FiveEq { eos: [WATER_SG, AIR], kapila: true };
            let mk = || {
                let mesh = Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]);
                let mut s = Solver::new(sys, mesh, |x| {
                    let (a, rho_w, rho_a, p) = if x[0] < 0.7 { (1.0 - 1e-8, 1000.0, 50.0, 1e9) } else { (1e-8, 1000.0, 50.0, 1e5) };
                    [a * rho_w, (1.0 - a) * rho_a, 0.0, 0.0, 0.0, p, a]
                });
                s.bc[0] = [Bc::Outflow; 2];
                s.cfl = 0.4;
                s
            };
            run("water | air tube (five-equation)", mk, GpuSystem::FiveEq { eos: [WATER_SG, AIR], kapila: true }, 2.4e-4, 5);
        }
        "taylor" => {
            let eos = Nasg::ideal(1.4);
            let pi = std::f64::consts::PI;
            let mk = || {
                let mesh = Mesh::new([n, n, n], [0.0; 3], [2.0 * pi; 3]);
                let (ma, g) = (0.1, 1.4);
                let p0 = 1.0 / (g * ma * ma);
                let mut s = Solver::new(Euler { eos }, mesh, |x| {
                    let (u, v) = (x[0].sin() * x[1].cos() * x[2].cos(), -x[0].cos() * x[1].sin() * x[2].cos());
                    let p = p0 + ((2.0 * x[0]).cos() + (2.0 * x[1]).cos()) * ((2.0 * x[2]).cos() + 2.0) / 16.0;
                    [p / p0, u, v, 0.0, p]
                });
                s.cfl = 0.4;
                s
            };
            run("3D compressible Taylor–Green (Ma 0.1, Euler)", mk, GpuSystem::Euler(eos), 0.5, 1);
        }
        _ => {
            let eos = Nasg::ideal(1.4);
            let mk = || {
                let mesh = Mesh::new([n, n, 1], [0.0; 3], [1.0, 1.0, 1.0]);
                let mut s = Solver::new(Euler { eos }, mesh, |x| {
                    let p = |r: f64, u: f64, v: f64, p: f64| [r, u, v, 0.0, p];
                    match (x[0] > 0.8, x[1] > 0.8) {
                        (true, true) => p(1.5, 0.0, 0.0, 1.5),
                        (false, true) => p(0.5323, 1.206, 0.0, 0.3),
                        (false, false) => p(0.138, 1.206, 1.206, 0.029),
                        (true, false) => p(0.5323, 0.0, 1.206, 0.3),
                    }
                });
                s.bc[0] = [Bc::Outflow; 2];
                s.bc[1] = [Bc::Outflow; 2];
                s.cfl = 0.4;
                s
            };
            run("2D Riemann config 3 (Euler)", mk, GpuSystem::Euler(eos), 0.3, 0);
        }
    }
}
