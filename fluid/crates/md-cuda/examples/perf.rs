//! Performance sweep: CPU (md crate, rayon) vs GPU f64/f32.
//!
//!   cargo run -p md-cuda --release --example perf -- [N ...]
//!
//! Default sizes are 4000, 32 000 and 256 000 atoms (fcc at rho = 0.8,
//! rc = 2.5). We take a few velocity-Verlet steps with the GPU f64/f32
//! engines (which rebuild the cell list every force evaluation) and, unless
//! `SKIP_CPU` is set, with the CPU engine. Each step is timed separately and
//! the *minimum* is reported: the machine is shared, so the mean is inflated
//! by contention and periodic neighbour-list rebuilds.
//!
//! The CPU `md` neighbour-list build for 256 000 atoms takes minutes on this
//! machine, so 256k is usually run with `SKIP_CPU=1`; the CPU 256k number is
//! taken from a separate force-evaluation measurement.
//!
//! The GPU is an RTX 3050 Laptop, software power-capped to 15 W (default
//! 35 W) and observed at its 210 MHz idle clock rather than the 2100 MHz
//! maximum, and it is shared with other CUDA work; treat the GPU numbers as
//! lower bounds.

use md::potential::{Cutoff, Potential};
use md::system::System;
use md_cuda::GpuMd;
use std::time::Instant;

/// Reduced-time unit of the argon Lennard-Jones map, in seconds. Constants
/// from `md::argon` (Rowley & Nicholson 1972): sigma = 3.405 Angstrom,
/// eps/k_B = 119.8 K, m = 39.948 u.
fn argon_time_unit_s() -> f64 {
    const SIGMA: f64 = 3.405e-10;
    const EPS: f64 = 119.8 * 1.380_649e-23;
    const M: f64 = 39.948 * 1.660_539_066_60e-27;
    SIGMA * (M / EPS).sqrt()
}

fn ns_per_day(steps_per_s: f64, dt: f64) -> f64 {
    steps_per_s * dt * argon_time_unit_s() * 86400.0 / 1e-9
}

fn lj() -> Potential {
    Potential::LennardJones(materials::LennardJones::new(1.0, 1.0))
}

fn fcc_positions(n_cells: usize, density: f64) -> (Vec<[f64; 3]>, f64) {
    let a = (4.0 / density).cbrt();
    let l = n_cells as f64 * a;
    let basis = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ];
    let mut pos = Vec::with_capacity(4 * n_cells * n_cells * n_cells);
    for i in 0..n_cells {
        for j in 0..n_cells {
            for k in 0..n_cells {
                for b in &basis {
                    pos.push([(i as f64 + b[0]) * a, (j as f64 + b[1]) * a, (k as f64 + b[2]) * a]);
                }
            }
        }
    }
    (pos, l)
}

fn row(n: usize, k: usize, engine: &str, best_s: f64, mean_s: f64, dt: f64, speedup: Option<f64>) {
    let sps = 1.0 / best_s;
    let sp_str = match speedup {
        Some(x) => format!("{x:>9.2}x"),
        None => format!("{:>10}", "-"),
    };
    println!(
        "{n:>9} {k:>6} {engine:>14} {sps:>12.1} {:>10.1} {sp_str} {:>10.3e} {:>10.3e}",
        ns_per_day(sps, dt),
        best_s,
        mean_s
    );
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let sizes: Vec<usize> = if args.is_empty() {
        vec![4000, 32_000, 256_000]
    } else {
        args.iter().map(|s| s.parse().unwrap()).collect()
    };
    let skip_cpu = std::env::var("SKIP_CPU").is_ok();
    let skip_gpu = std::env::var("SKIP_GPU").is_ok();
    let gpu = !skip_gpu && md_cuda::gpu_available();
    if !gpu {
        eprintln!("no CUDA device available; running CPU only");
    }

    let dt = 0.005;
    println!(
        "argon time unit = {:.4e} s; dt = {dt} reduced = {:.4e} ps; nvidia-smi: {}",
        argon_time_unit_s(),
        dt * argon_time_unit_s() * 1e12,
        gpu_clock()
    );
    println!(
        "{:>9} {:>6} {:>14} {:>12} {:>10} {:>10} {:>10} {:>10}",
        "N", "k", "engine", "steps/s", "ns/day", "speedup", "best s", "mean s"
    );

    for &n_target in &sizes {
        let k = ((n_target as f64 / 4.0).cbrt().round() as usize).max(1);
        let n = 4 * k * k * k;
        let steps = if n >= 200_000 { 5 } else if n >= 20_000 { 20 } else { 100 };

        let (pos, box_l) = fcc_positions(k, 0.8);
        let mut rng = md::Rng::new(3);
        let vel: Vec<[f64; 3]> = (0..n)
            .map(|_| [rng.normal(), rng.normal(), rng.normal()])
            .collect();

        // ---- CPU ----
        let mut cpu_best = f64::INFINITY;
        if !skip_cpu {
            let mut cpu = System::fcc(k, 0.8, lj(), 2.5, 0.3, Cutoff::Shifted, 3);
            cpu.set_velocities_mb(1.0);
            let mut times = Vec::with_capacity(steps);
            // First step warms up and is not counted.
            cpu.step_nve(dt);
            for _ in 0..steps {
                let t = Instant::now();
                cpu.step_nve(dt);
                times.push(t.elapsed().as_secs_f64());
            }
            let best = times.iter().cloned().fold(f64::INFINITY, f64::min);
            let mean = times.iter().sum::<f64>() / times.len() as f64;
            cpu_best = best;
            row(n, k, "CPU f64", best, mean, dt, None);
        }

        if !gpu {
            continue;
        }
        // ---- GPU f64 ----
        {
            let mut g = GpuMd::<f64>::new(&pos, &vel, box_l, 2.5, &lj(), Cutoff::Shifted).unwrap();
            g.step(dt).unwrap();
            let mut times = Vec::with_capacity(steps);
            for _ in 0..steps {
                let t = Instant::now();
                g.step(dt).unwrap();
                g.synchronize().unwrap();
                times.push(t.elapsed().as_secs_f64());
            }
            let best = times.iter().cloned().fold(f64::INFINITY, f64::min);
            let mean = times.iter().sum::<f64>() / times.len() as f64;
            let sp = if cpu_best.is_finite() { Some(cpu_best / best) } else { None };
            row(n, k, "GPU f64", best, mean, dt, sp);
        }
        // ---- GPU f32 ----
        {
            let mut g = GpuMd::<f32>::new(&pos, &vel, box_l, 2.5, &lj(), Cutoff::Shifted).unwrap();
            g.step(dt).unwrap();
            let mut times = Vec::with_capacity(steps);
            for _ in 0..steps {
                let t = Instant::now();
                g.step(dt).unwrap();
                g.synchronize().unwrap();
                times.push(t.elapsed().as_secs_f64());
            }
            let best = times.iter().cloned().fold(f64::INFINITY, f64::min);
            let mean = times.iter().sum::<f64>() / times.len() as f64;
            let sp = if cpu_best.is_finite() { Some(cpu_best / best) } else { None };
            row(n, k, "GPU f32", best, mean, dt, sp);
        }
    }
}

/// Best-effort read of the current SM clock, for the log only.
fn gpu_clock() -> String {
    std::process::Command::new("nvidia-smi")
        .args(["--query-gpu=clocks.sm,power.draw", "--format=csv,noheader"])
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .unwrap_or_else(|| "unknown".into())
}
