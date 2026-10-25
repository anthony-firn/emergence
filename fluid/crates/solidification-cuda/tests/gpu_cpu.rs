//! GPU-vs-CPU checks. These tests skip gracefully if no CUDA device or the
//! CUDA runtime libraries are unavailable. Run with
//!   source fluid/scripts/cuda-env.sh && cargo test -p solidification-cuda

use solidification::phase_field::{Bc, PhaseField};
use solidification_cuda::{PhaseFieldGpu, Real};

fn gpu_available() -> bool {
    cudarc::driver::CudaContext::new(0).is_ok()
}

fn seed(n: usize, dx: f64, w0: f64) -> PhaseField {
    let mut pf = PhaseField::new(n, n, dx, 1.0, 2.0, w0, 0.05);
    pf.bc = Bc::quadrant();
    pf.u_far = -0.55;
    pf.seed_parabola(8.0, 4.0);
    pf
}

/// f64 on the GPU must reproduce the CPU update to round-off.
#[test]
fn gpu_f64_matches_cpu() {
    if !gpu_available() {
        eprintln!("no CUDA device / libraries: skipping GPU test");
        return;
    }
    let mut cpu = seed(48, 0.4, 0.5);
    let dt = cpu.stable_dt();
    let mut gpu = PhaseFieldGpu::<f64>::from_cpu(&cpu).unwrap();

    let steps = 40;
    for _ in 0..steps {
        cpu.step(dt);
        gpu.step(dt).unwrap();
    }
    let gp = gpu.download_p().unwrap();
    let gu = gpu.download_u().unwrap();
    let (mut ep, mut eu) = (0.0f64, 0.0f64);
    // Compare interior cells only: the CPU leaves ghosts from the start-of-step
    // `apply_bc`, while the GPU's swapped buffer keeps its (unused) scratch
    // ghosts until the next step. Ghosts of both are refreshed before use.
    let m = cpu.nx + 2;
    for j in 1..=cpu.ny {
        for i in 1..=cpu.nx {
            let k = j * m + i;
            ep = ep.max((cpu.p[k] - gp[k]).abs());
            eu = eu.max((cpu.u[k] - gu[k]).abs());
        }
    }
    println!("GPU f64 vs CPU after {steps} steps: max|dp|={ep:.3e} max|du|={eu:.3e}");
    assert!(ep < 1e-10, "phase field drift {ep:.3e}");
    assert!(eu < 1e-10, "temperature field drift {eu:.3e}");
}

/// f32 is close to f64 over a short run (the production precision).
#[test]
fn gpu_f32_close_to_f64() {
    if !gpu_available() {
        eprintln!("no CUDA device / libraries: skipping GPU test");
        return;
    }
    let cpu = seed(48, 0.4, 0.5);
    let dt = cpu.stable_dt();
    let mut g64 = PhaseFieldGpu::<f64>::from_cpu(&cpu).unwrap();
    let mut g32 = PhaseFieldGpu::<f32>::from_cpu(&cpu).unwrap();
    for _ in 0..40 {
        g64.step(dt).unwrap();
        g32.step(dt).unwrap();
    }
    let a = g64.download_p().unwrap();
    let b = g32.download_p().unwrap();
    let e = a.iter().zip(b.iter()).map(|(x, y)| (x - y).abs()).fold(0.0f64, f64::max);
    println!("GPU f32 vs f64 after 40 steps: max|dp|={e:.3e}");
    assert!(e < 2e-3, "f32 deviates {e:.3e}");
}

/// The `Real` conversion used by the port.
#[test]
fn real_roundtrip() {
    assert_eq!(<f32 as Real>::from_f64(0.5).to_f64(), 0.5);
    assert_eq!(<f64 as Real>::from_f64(0.5).to_f64(), 0.5);
}

/// The variational anisotropy term (KR 1998 Eq. 88) must actually change the
/// update, and the 9-point Laplacian must differ from the 5-point one. This is
/// what makes the corrected operator select a steady dendrite where the
/// `solidification` crate's operator does not.
#[test]
fn kr_operator_differs_from_plain_operator() {
    if !gpu_available() {
        eprintln!("no CUDA device / libraries: skipping GPU test");
        return;
    }
    // `seed` has eps4 = 0.05, so both the variational and grid terms are live.
    let cpu = seed(32, 0.4, 1.0);
    let dt = cpu.stable_dt();
    let run = |mode: i32| {
        let mut g = PhaseFieldGpu::<f64>::from_cpu(&cpu).unwrap();
        g.run_kr_mode(20, dt, mode).unwrap();
        g.download_p().unwrap()
    };
    let m3 = run(3); // 9-point + variational (production KR)
    let m1 = run(1); // 9-point, no variational
    let m0 = run(0); // 5-point, no variational (the CPU crate's operator)
    let diff = |a: &[f64], b: &[f64]| {
        a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0f64, f64::max)
    };
    let d_variational = diff(&m3, &m1);
    let d_stencil = diff(&m1, &m0);
    println!(
        "KR mode3 vs mode1 (variational): {d_variational:.3e}; mode1 vs mode0 (stencil): {d_stencil:.3e}"
    );
    assert!(d_variational > 1e-6, "variational term inactive");
    assert!(d_stencil > 1e-6, "9-point stencil equals 5-point");
}
