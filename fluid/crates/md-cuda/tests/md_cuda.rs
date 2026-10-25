//! GPU vs CPU checks for the Lennard-Jones engine.
//!
//! Every test skips with a printed note when no CUDA device is present.
//! All runs are small (108 atoms, a few hundred steps) so the suite stays
//! well under the 30 s budget.

use md::potential::{Cutoff, Potential};
use md::system::System;
use md_cuda::gk::running_gk_viscosity;
use md_cuda::{gpu_available, GpuMd};

fn lj() -> Potential {
    Potential::LennardJones(materials::LennardJones::new(1.0, 1.0))
}

fn skip(name: &str) {
    eprintln!("[skip] {name}: no CUDA device/driver available");
}

fn force_rel_err(a: &[[f64; 3]], b: &[[f64; 3]]) -> (f64, f64) {
    let mut md = 0.0f64;
    let mut scale = 0.0f64;
    for (x, y) in a.iter().zip(b) {
        for k in 0..3 {
            md = md.max((x[k] - y[k]).abs());
            scale = scale.max(x[k].abs());
        }
    }
    (md, md / scale.max(f64::MIN_POSITIVE))
}

/// A disordered (liquid) configuration: fcc at rho = 0.8 melted by a short
/// NVE run at T = 1.
fn liquid() -> System {
    let mut sys = System::fcc(3, 0.8, lj(), 2.5, 0.3, Cutoff::Shifted, 20240607);
    sys.set_velocities_mb(1.0);
    for _ in 0..400 {
        sys.step_nve(0.005);
    }
    sys
}

fn report_state(name: &str, cpu: &System, gpu: &mut GpuMd<f64>) {
    let gf = gpu.forces().unwrap();
    let (md, rel) = force_rel_err(&gf, &cpu.force);
    let ep_cpu = cpu.potential_energy();
    let ep_gpu = gpu.potential_energy().unwrap();
    let vir_cpu = cpu.internal_pair_virial();
    let vir_gpu = gpu.internal_pair_virial().unwrap();
    let e_rel = (ep_gpu - ep_cpu).abs() / ep_cpu.abs();
    let v_rel = (vir_gpu - vir_cpu).abs() / vir_cpu.abs();
    println!(
        "{name} N={} f64: force max|d|={md:.3e} rel={rel:.3e}; PE rel={e_rel:.3e} ({ep_cpu:.15}); virial rel={v_rel:.3e} ({vir_cpu:.15})",
        cpu.n()
    );
    assert!(rel < 1e-12, "force relative error {rel:e}");
    assert!(e_rel < 1e-12, "energy relative error {e_rel:e}");
    assert!(v_rel < 1e-12, "virial relative error {v_rel:e}");
}

#[test]
fn fcc_force_energy_virial_match_cpu() {
    if !gpu_available() {
        return skip("fcc_force_energy_virial_match_cpu");
    }
    let mut cpu = System::fcc(3, 0.8, lj(), 2.5, 0.3, Cutoff::Shifted, 7);
    // A perfect fcc lattice has exactly cancelling forces, so the relative
    // force error is dominated by round-off of a near-zero quantity. Give
    // each atom a small deterministic displacement (a thermally displaced
    // fcc) so the forces are O(1) and the relative metric is meaningful.
    let mut rng = md::Rng::new(7);
    for p in cpu.pos.iter_mut() {
        for a in 0..3 {
            p[a] += rng.uniform(-0.05, 0.05);
        }
    }
    cpu.rebuild_neighbors();
    cpu.set_velocities_mb(1.2);
    cpu.compute_forces();
    let mut gpu = GpuMd::<f64>::from_md(&cpu).unwrap();
    report_state("fcc+0.05", &cpu, &mut gpu);
}

#[test]
fn liquid_force_energy_virial_match_cpu() {
    if !gpu_available() {
        return skip("liquid_force_energy_virial_match_cpu");
    }
    let mut cpu = liquid();
    cpu.compute_forces();
    let mut gpu = GpuMd::<f64>::from_md(&cpu).unwrap();
    gpu.compute_forces().unwrap();
    report_state("liquid", &cpu, &mut gpu);
}

#[test]
fn nve_trajectory_tracks_cpu_f64() {
    if !gpu_available() {
        return skip("nve_trajectory_tracks_cpu_f64");
    }
    let mut cpu = System::fcc(3, 0.8, lj(), 2.5, 0.3, Cutoff::Shifted, 99);
    cpu.set_velocities_mb(1.0);
    let mut gpu = GpuMd::<f64>::new(&cpu.pos, &cpu.vel, cpu.box_l, cpu.rc, &cpu.potential, cpu.scheme)
        .unwrap();
    let dt = 0.002;
    let e_cpu0 = cpu.total_energy();
    let e_gpu0 = gpu.total_energy().unwrap();
    println!(
        "initial energy: CPU {e_cpu0:.15}, GPU {e_gpu0:.15}, rel diff {:.3e}",
        (e_gpu0 - e_cpu0).abs() / e_cpu0.abs()
    );
    assert!((e_gpu0 - e_cpu0).abs() / e_cpu0.abs() < 1e-12);

    let steps = 2000;
    let mut report: Vec<usize> = vec![1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000];
    let mut err1 = 0.0f64;
    let (mut dev_cpu, mut dev_gpu) = (0.0f64, 0.0f64);
    for s in 1..=steps {
        cpu.step_nve(dt);
        gpu.step(dt).unwrap();
        if report.first() == Some(&s) {
            report.remove(0);
            let gp = gpu.positions().unwrap();
            let mut md = 0.0f64;
            for (a, b) in gp.iter().zip(&cpu.pos) {
                for k in 0..3 {
                    md = md.max((a[k] - b[k]).abs());
                }
            }
            if s == 1 {
                err1 = md;
            }
            println!("  step {s:4}: max |x_gpu - x_cpu| = {md:.3e}");
        }
        if s % 100 == 0 {
            dev_cpu = dev_cpu.max((cpu.total_energy() - e_cpu0).abs() / e_cpu0.abs());
            dev_gpu = dev_gpu.max((gpu.total_energy().unwrap() - e_gpu0).abs() / e_gpu0.abs());
        }
    }
    println!(
        "  first-step max position error {err1:.3e}; max energy deviation over {steps} steps: CPU {dev_cpu:.3e}, GPU f64 {dev_gpu:.3e}"
    );
    assert!(err1 < 1e-9, "first-step position error {err1:e}");
    assert!(dev_gpu < 1e-4, "GPU f64 energy deviation {dev_gpu:e}");
}

/// The scan kernel walks the cells in blocks of 256, and the neighbour
/// logic special-cases small grids (`nc <= 2`). Cover a grid larger than
/// one scan chunk (343 cells) and the minimal 3-cell image to make sure the
/// cell-list offsets, not just the 8-cell `fcc(3)` case, match the CPU.
#[test]
fn large_cell_grid_matches_cpu() {
    if !gpu_available() {
        return skip("large_cell_grid_matches_cpu");
    }
    for k in [11usize, 5] {
        let mut cpu = System::fcc(k, 0.8, lj(), 2.5, 0.3, Cutoff::Shifted, 4242);
        let mut rng = md::Rng::new(k as u64 + 1);
        for p in cpu.pos.iter_mut() {
            for a in 0..3 {
                p[a] += rng.uniform(-0.05, 0.05);
            }
        }
        cpu.rebuild_neighbors();
        cpu.set_velocities_mb(1.0);
        cpu.compute_forces();
        let mut gpu = GpuMd::<f64>::from_md(&cpu).unwrap();
        report_state(&format!("fcc k={k}"), &cpu, &mut gpu);
    }
}

/// `Cutoff::TailCorrected` uses the bare potential on the device and adds
/// the analytic tail correction on the host; check energy and virial (which
/// must *not* include the tail in the raw pair virial).
#[test]
fn tail_corrected_energy_matches_cpu() {
    if !gpu_available() {
        return skip("tail_corrected_energy_matches_cpu");
    }
    let mut cpu = System::fcc(4, 0.8, lj(), 2.5, 0.3, Cutoff::TailCorrected, 5);
    // Displace so the forces are O(1): a perfect lattice has exactly
    // cancelling forces and the relative-force metric would be round-off.
    let mut rng = md::Rng::new(5);
    for p in cpu.pos.iter_mut() {
        for a in 0..3 {
            p[a] += rng.uniform(-0.05, 0.05);
        }
    }
    cpu.rebuild_neighbors();
    cpu.set_velocities_mb(1.0);
    cpu.compute_forces();
    let mut gpu = GpuMd::<f64>::from_md(&cpu).unwrap();
    let e_cpu = cpu.potential_energy();
    let e_gpu = gpu.potential_energy().unwrap();
    let e_rel = (e_gpu - e_cpu).abs() / e_cpu.abs();
    println!("TailCorrected PE: CPU {e_cpu:.15}, GPU {e_gpu:.15}, rel {e_rel:.3e}");
    assert!(e_rel < 1e-12, "tail-corrected energy relative error {e_rel:e}");
    report_state("fcc tail-corrected", &cpu, &mut gpu);
}

/// The scan walks `ncells` in chunks of 256. Fix a dilute random
/// configuration whose grid has 13 824 cells (54 chunks) and compare the
/// device per-cell occupancy to a CPU binning exactly: a wrong carry would
/// put atoms in the wrong cells and show up here as an integer mismatch.
#[test]
fn cell_list_scan_matches_cpu_binning() {
    if !gpu_available() {
        return skip("cell_list_scan_matches_cpu_binning");
    }
    let n = 4000usize;
    let l = 60.0f64;
    let rc = 2.5f64;
    let nc = (l / rc).floor() as i64; // 24 per axis -> 13 824 cells
    assert!(nc > 2);
    let mut rng = md::Rng::new(7);
    let pos: Vec<[f64; 3]> = (0..n)
        .map(|_| [rng.uniform(0.0, l), rng.uniform(0.0, l), rng.uniform(0.0, l)])
        .collect();
    let vel = vec![[0.0f64; 3]; n];
    let pot = lj();
    let mut gpu = GpuMd::<f64>::new(&pos, &vel, l, rc, &pot, Cutoff::Shifted).unwrap();
    let got = gpu.cell_counts().unwrap();
    let ncells = (nc * nc * nc) as usize;
    assert_eq!(got.len(), ncells);
    let mut want = vec![0i32; ncells];
    for p in &pos {
        let ix = ((p[0] / l * nc as f64).floor() as i64).rem_euclid(nc) as usize;
        let iy = ((p[1] / l * nc as f64).floor() as i64).rem_euclid(nc) as usize;
        let iz = ((p[2] / l * nc as f64).floor() as i64).rem_euclid(nc) as usize;
        want[ix + (nc as usize) * (iy + (nc as usize) * iz)] += 1;
    }
    let mismatches: Vec<usize> = (0..ncells).filter(|&c| got[c] != want[c]).collect();
    assert!(
        mismatches.is_empty(),
        "{} of {ncells} cells differ (first at {:?}: gpu {} vs cpu {})",
        mismatches.len(),
        mismatches.first(),
        mismatches.first().map(|&c| got[c]).unwrap_or(0),
        mismatches.first().map(|&c| want[c]).unwrap_or(0),
    );
    assert_eq!(got.iter().sum::<i32>(), n as i32);
    println!("cell-list scan: {ncells} cells ({nc}^3), occupancy matches CPU binning exactly");
}

/// The `f64`-accumulated pressure tensor must match the CPU engine on the
/// same configuration, including the isotropic tail pressure; the `f32`
/// force path keeps f64 accumulation so its off-diagonals are close too.
#[test]
fn pressure_tensor_matches_cpu() {
    if !gpu_available() {
        return skip("pressure_tensor_matches_cpu");
    }
    let mut cpu = liquid();
    for scheme in [Cutoff::Shifted, Cutoff::TailCorrected] {
        cpu.scheme = scheme;
        cpu.rebuild_neighbors();
        cpu.compute_forces();
        let pc = cpu.pressure_tensor();
        let mut g64 = GpuMd::<f64>::from_md(&cpu).unwrap();
        let p64 = g64.pressure_tensor_f64().unwrap();
        let mut g32 = GpuMd::<f32>::from_md(&cpu).unwrap();
        let p32 = g32.pressure_tensor_f64().unwrap();
        for k in 0..6 {
            let d64 = (p64[k] - pc[k]).abs();
            let d32 = (p32[k] - pc[k]).abs();
            assert!(d64 < 1e-11, "{scheme:?} comp {k}: f64 {d64:e} ({} vs {})", p64[k], pc[k]);
            assert!(d32 < 2e-4, "{scheme:?} comp {k}: f32 {d32:e} ({} vs {})", p32[k], pc[k]);
        }
    }
}

/// End-to-end smoke test of the device pressure sampling and the
/// block-averaged Green–Kubo estimator on a short `f32` trajectory.
#[test]
fn pressure_sampling_and_gk_smoke() {
    if !gpu_available() {
        return skip("pressure_sampling_and_gk_smoke");
    }
    // 500 atoms, fcc k=5, tail-corrected, small but enough to make the
    // estimator exercise a real trajectory.
    let mut cpu = System::fcc(5, 0.8442, lj(), 2.5, 0.3, Cutoff::TailCorrected, 17);
    cpu.set_velocities_mb(0.722);
    let mut gpu = GpuMd::<f32>::new(
        &cpu.pos, &cpu.vel, cpu.box_l, cpu.rc, &cpu.potential, cpu.scheme,
    )
    .unwrap();
    // Melt the metastable fcc and equilibrate at the triple-point liquid
    // state, then run NVE. Rescaling every 20 steps is only an
    // equilibration device; the production trajectory is microcanonical.
    let dt = 0.004;
    for s in 0..1500 {
        gpu.step(dt).unwrap();
        if s % 20 == 0 {
            gpu.rescale_to_temperature(1.5).unwrap();
        }
    }
    for s in 0..4000 {
        gpu.step(dt).unwrap();
        if s % 20 == 0 {
            gpu.rescale_to_temperature(0.722).unwrap();
        }
    }
    let t_eq = gpu.temperature().unwrap();
    let steps = 3000usize;
    gpu.reserve_pressure_samples(steps).unwrap();
    let mut tsum = 0.0;
    for _ in 0..steps {
        gpu.step(dt).unwrap();
        gpu.sample_pressure_tensor().unwrap();
        tsum += gpu.temperature().unwrap();
    }
    let press = gpu.take_pressure_samples().unwrap();
    assert_eq!(press.len(), steps);
    let kbt = tsum / steps as f64;
    let vol = cpu.box_l.powi(3);
    let curve = md_cuda::gk::running_gk_viscosity(&press, dt, vol, kbt, 600, 4);
    let (eta, err, lag) = md_cuda::gk::plateau_value(&curve);
    let (eta_eh, _) = md_cuda::gk::einstein_helfand_viscosity(&press, dt, vol, kbt, 600, 4, 100, 500);
    println!(
        "N=500 f32 short run: T*(eq)={t_eq:.4}, <T*>(prod)={kbt:.4}, eta_GK*={eta:.3}+/-{err:.3} at lag {lag}, eta_EH*={eta_eh:.3}"
    );
    assert!(eta.is_finite() && err.is_finite());
    assert!(eta > 0.5 && eta < 8.0, "implausible eta* {eta}");
    assert!(eta_eh.is_finite());
}

#[test]
fn gk_running_integral_matches_cpu_md_crate() {
    if !gpu_available() {
        return skip("gk_running_integral_matches_cpu_md_crate");
    }
    // Same small system, same f64 trajectory: the GPU pressure-tensor samples
    // must feed the same Green-Kubo running integral as the CPU `md` crate.
    let mut cpu = System::fcc(3, 0.8, lj(), 2.5, 0.3, Cutoff::Shifted, 4242);
    cpu.set_velocities_mb(0.9);
    for _ in 0..100 {
        cpu.step_nve(0.004);
    }
    cpu.compute_forces();
    let mut gpu = GpuMd::<f64>::new(
        &cpu.pos, &cpu.vel, cpu.box_l, cpu.rc, &cpu.potential, cpu.scheme,
    )
    .unwrap();
    let dt = 0.004;
    let steps = 400usize;
    gpu.reserve_pressure_samples(steps).unwrap();
    let mut cpu_press = Vec::with_capacity(steps);
    for _ in 0..steps {
        cpu.step_nve(dt);
        gpu.step(dt).unwrap();
        cpu_press.push(cpu.pressure_tensor());
        gpu.sample_pressure_tensor().unwrap();
    }
    let gpu_press = gpu.take_pressure_samples().unwrap();
    assert_eq!(gpu_press.len(), steps);
    let mut maxd = 0.0f64;
    for (a, b) in cpu_press.iter().zip(&gpu_press) {
        for k in 0..6 {
            maxd = maxd.max((a[k] - b[k]).abs());
        }
    }
    println!("GPU f64 vs CPU pressure tensor, max|d| over {steps} steps = {maxd:.3e}");
    assert!(maxd < 1e-10, "pressure-tensor sample mismatch {maxd:e}");

    // Running GK integral, GPU block estimator (one block) vs the scalar
    // `md::transport::running_gk` applied per off-diagonal component.
    let max_lag = 200usize;
    let vol = cpu.box_l.powi(3);
    let kbt = 0.9;
    let curve = running_gk_viscosity(&gpu_press, dt, vol, kbt, max_lag, 1);
    let mut cpu_run = vec![0.0f64; max_lag + 1];
    for &c in &[3usize, 4, 5] {
        let series: Vec<f64> = cpu_press.iter().map(|p| p[c]).collect();
        let rk = md::transport::running_gk(&series, dt, max_lag);
        for (k, r) in rk.iter().enumerate() {
            cpu_run[k] += vol / (3.0 * kbt) * r;
        }
    }
    for &k in &[10usize, 50, 100, 200] {
        let d = (curve.mean[k] - cpu_run[k]).abs();
        let scale = cpu_run[k].abs().max(1.0);
        assert!(d < 1e-9 * scale, "lag {k}: GPU {} vs CPU {}", curve.mean[k], cpu_run[k]);
    }
    println!(
        "GK running integral GPU vs CPU(md): lag200 {:.10} vs {:.10}",
        curve.mean[max_lag], cpu_run[max_lag]
    );
}

#[test]
fn f32_precision_and_energy_drift_are_reported() {
    if !gpu_available() {
        return skip("f32_precision_and_energy_drift_are_reported");
    }
    let mut cpu0 = System::fcc(3, 0.8, lj(), 2.5, 0.3, Cutoff::Shifted, 1234);
    // Displace slightly (as in the fcc test) so the forces are O(1).
    let mut rng = md::Rng::new(11);
    for p in cpu0.pos.iter_mut() {
        for a in 0..3 {
            p[a] += rng.uniform(-0.05, 0.05);
        }
    }
    cpu0.rebuild_neighbors();
    cpu0.set_velocities_mb(1.0);
    cpu0.compute_forces();
    let e_ref = cpu0.total_energy();

    // Single-configuration force precision of f32 against the CPU f64 force.
    let g32_0 = GpuMd::<f32>::new(&cpu0.pos, &cpu0.vel, cpu0.box_l, cpu0.rc, &cpu0.potential, cpu0.scheme)
        .unwrap();
    let (md32, rel32) = force_rel_err(&g32_0.forces().unwrap(), &cpu0.force);
    println!("f32 vs CPU f64 force: max|d| = {md32:.3e}, rel = {rel32:.3e}");

    let dt = 0.002;
    let steps = 2000;
    for (label, tag) in [("f64", 0), ("f32", 1)] {
        if tag == 0 {
            let mut g64 = GpuMd::<f64>::new(
                &cpu0.pos, &cpu0.vel, cpu0.box_l, cpu0.rc, &cpu0.potential, cpu0.scheme,
            )
            .unwrap();
            let e0 = g64.total_energy().unwrap();
            let mut dev = 0.0f64;
            for s in 0..steps {
                g64.step(dt).unwrap();
                if (s + 1) % 100 == 0 {
                    dev = dev.max((g64.total_energy().unwrap() - e0).abs() / e0.abs());
                }
            }
            println!("{label} NVE {steps} steps: max |dE|/|E| = {dev:.3e}");
        } else {
            let mut g32 = GpuMd::<f32>::new(
                &cpu0.pos, &cpu0.vel, cpu0.box_l, cpu0.rc, &cpu0.potential, cpu0.scheme,
            )
            .unwrap();
            let e0 = g32.total_energy().unwrap();
            let mut dev = 0.0f64;
            for s in 0..steps {
                g32.step(dt).unwrap();
                if (s + 1) % 100 == 0 {
                    dev = dev.max((g32.total_energy().unwrap() - e0).abs() / e0.abs());
                }
            }
            let final_rel = (g32.total_energy().unwrap() - e_ref).abs() / e_ref.abs();
            println!("{label} NVE {steps} steps: max |dE|/|E| = {dev:.3e}, |E_end - E_CPU(0)|/|E| = {final_rel:.3e}");
            assert!(dev < 1e-2, "f32 energy deviation {dev:e} unexpectedly large");
        }
    }
}
