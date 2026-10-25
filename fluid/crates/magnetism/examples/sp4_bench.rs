//! Timing harness for the µMAG SP4 hot loop (`DemagFft::field_ws` and
//! `Llg3d::step_ws`).  Used to produce the performance table in
//! `analysis/magnetism.md`; not needed for production runs, which use
//! `examples/mumag_sp4.rs`.
use magnetism::demag::{DemagFft, FftWorkspace};
use magnetism::llg3d::Llg3d;
use magnetism::micromag::{GAMMA, MU0};
use magnetism::sp4;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let nx: usize = args.get(1).map(|s| s.parse().unwrap()).unwrap_or(128);
    let ny: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(32);
    let nz: usize = args.get(3).map(|s| s.parse().unwrap()).unwrap_or(1);
    let dx = sp4::LX / nx as f64;
    let dy = sp4::LY / ny as f64;
    let dz = sp4::LZ / nz as f64;
    let mut s = Llg3d::new(nx, ny, nz, dx, dy, dz, sp4::A_EX, 0.0, sp4::MS, sp4::ALPHA);
    s.demag = Some(DemagFft::new(nx, ny, nz, [dx, dy, dz]));
    s.easy = [1.0, 0.0, 0.0];
    let mut ws = s.make_scratch();

    let padded = s.demag.as_ref().unwrap().padded_len();
    let mut out = vec![[0.0; 3]; s.n_cells()];
    let mut fws = FftWorkspace::new(padded);
    let reps = 50;
    let t0 = Instant::now();
    for _ in 0..reps {
        s.demag.as_ref().unwrap().field_ws(&s.m, &mut out, &mut fws);
    }
    let dt_field = t0.elapsed().as_secs_f64() / reps as f64;

    s.set_uniform([1.0, 1.0, 1.0]);
    let dt_relax = 0.1 * dx * dx * MU0 * sp4::MS / (2.0 * sp4::A_EX * GAMMA);
    let max_relax: usize = args.get(4).map(|s| s.parse().unwrap()).unwrap_or(200_000);
    let t0 = Instant::now();
    let steps = s.relax_ws(dt_relax, 1.0e2, max_relax, &mut ws);
    let relax_s = t0.elapsed().as_secs_f64();

    s.h_ext = [-24.6e-3 / MU0, 4.3e-3 / MU0, 0.0];
    let c_ex = 2.0 * sp4::A_EX / (MU0 * sp4::MS);
    let dt = 0.05 * dx * dx / (c_ex * GAMMA);
    let reps = 50;
    let t0 = Instant::now();
    for _ in 0..reps {
        s.step_ws(dt, &mut ws);
    }
    let step = t0.elapsed().as_secs_f64() / reps as f64;

    println!(
        "grid {}x{}x{} padded {}: field_ws {:.3} ms | {} relax steps {:.1} s | step {:.3} ms -> 0.6ns {} steps {:.0} s | total {:.0} s",
        nx, ny, nz, padded, dt_field * 1e3, steps, relax_s, step * 1e3,
        (0.6e-9 / dt).ceil(), (0.6e-9 / dt).ceil() as f64 * step, relax_s + (0.6e-9 / dt).ceil() as f64 * step
    );
}
