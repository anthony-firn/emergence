//! GPU vs CPU for the Karma-Rappel phase-field step, plus a timing/speedup
//! measurement.
//!
//!   cargo run --release -p solidification-cuda --example gpu_validate [N] [STEPS]
//!
//! f64 on the GPU must reproduce the CPU to round-off (interior cells); f32 is
//! compared against f64. Requires the CUDA runtime on the loader path
//! (`source scripts/cuda-env.sh`).

use solidification::phase_field::{Bc, PhaseField};
use solidification_cuda::{PhaseFieldGpu, Real};
use std::time::Instant;

fn make(n: usize, dx: f64, w0: f64) -> PhaseField {
    let mut pf = PhaseField::new(n, n, dx, 1.0, 2.0, w0, 0.05);
    pf.bc = Bc::quadrant();
    pf.u_far = -0.55;
    pf.seed_parabola(0.1 * n as f64 * dx, 4.0);
    pf
}

fn max_diff_interior(a: &[f64], b: &[f64], nx: usize, ny: usize) -> f64 {
    let m = nx + 2;
    let mut e = 0.0f64;
    for j in 1..=ny {
        for i in 1..=nx {
            let k = j * m + i;
            e = e.max((a[k] - b[k]).abs());
        }
    }
    e
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let n: usize = a.get(1).map(|v| v.parse().unwrap()).unwrap_or(128);
    let steps: usize = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(300);
    let dx = 0.4;
    let w0 = 1.0;
    let cpu = make(n, dx, w0);
    let dt = cpu.stable_dt();

    let mut c = make(n, dx, w0);
    let t0 = Instant::now();
    for _ in 0..steps {
        c.step(dt);
    }
    let t_cpu = t0.elapsed().as_secs_f64();

    let mut g64 = PhaseFieldGpu::<f64>::from_cpu(&make(n, dx, w0)).unwrap();
    let t0 = Instant::now();
    g64.run(steps, dt).unwrap();
    g64.synchronize().unwrap();
    let t_g64 = t0.elapsed().as_secs_f64();

    let mut g32 = PhaseFieldGpu::<f32>::from_cpu(&make(n, dx, w0)).unwrap();
    let t0 = Instant::now();
    g32.run(steps, dt).unwrap();
    g32.synchronize().unwrap();
    let t_g32 = t0.elapsed().as_secs_f64();

    let (cp, cu) = (c.p.clone(), c.u.clone());
    let (p64, u64) = (g64.download_p().unwrap(), g64.download_u().unwrap());
    let (p32, u32) = (g32.download_p().unwrap(), g32.download_u().unwrap());

    let cells = (n * n) as f64;
    println!("phase-field step: {n}x{n} = {} cells, {steps} steps, dt={dt:.3e}", n * n);
    println!(
        "  CPU f64  : {t_cpu:8.3}s  {:.3} Gcell/s",
        cells * steps as f64 / t_cpu / 1e9
    );
    println!(
        "  GPU f64  : {t_g64:8.3}s  {:.3} Gcell/s  ({:.1}x CPU)  max|dp| vs CPU {:.3e}",
        cells * steps as f64 / t_g64 / 1e9,
        t_cpu / t_g64,
        max_diff_interior(&cp, &p64, n, n)
    );
    println!(
        "  GPU f32  : {t_g32:8.3}s  {:.3} Gcell/s  ({:.1}x CPU)  max|dp| vs f64 {:.3e}",
        cells * steps as f64 / t_g32 / 1e9,
        t_cpu / t_g32,
        max_diff_interior(&p64, &p32, n, n)
    );
    let _ = cu;
    let _ = u64;
    let _ = u32;
    let _ = <f64 as Real>::SUFFIX;
}
