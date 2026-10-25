//! GPU backend vs. exact solutions and vs. the f64 CPU reference.
//! Tolerances reflect f32 storage (~1e-7 relative per operation).

use spectral_core::{init, Backend, Grid, Solver};
use spectral_cuda::CudaSolver;

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)
}

#[test]
fn gpu_abc_beltrami_decays_exactly() {
    let (n, nu, k) = (32, 0.05, 2.0);
    let g = Grid::new(n);
    let u0 = init::abc(&g, 1.0, 0.7, 0.4, k);
    let mut s = CudaSolver::new(g, nu, u0.clone()).expect("cuda init");
    for _ in 0..100 {
        s.step(0.01);
    }
    let decay = (-nu * k * k * s.time()).exp();
    for c in 0..3 {
        let exact: Vec<f64> = u0[c].iter().map(|v| v * decay).collect();
        let err = max_abs_diff(&s.velocity_phys(c), &exact);
        eprintln!("GPU ABC component {c}: max err {err:.2e}");
        assert!(err < 1e-6, "component {c} err {err:e}");
    }
}

#[test]
fn gpu_matches_cpu_reference_on_taylor_green() {
    let (n, re, dt, steps) = (32, 1600.0, 0.01, 200);
    let g = Grid::new(n);
    let mut cpu = Solver::new(g.clone(), 1.0 / re, init::taylor_green(&g));
    let mut gpu = CudaSolver::new(g.clone(), 1.0 / re, init::taylor_green(&g)).expect("cuda init");
    for _ in 0..steps {
        Backend::step(&mut cpu, dt);
        gpu.step(dt);
    }
    let (dc, dg) = (Backend::diagnostics(&mut cpu), gpu.diagnostics());
    let rel = |a: f64, b: f64| ((a - b) / a).abs();
    eprintln!(
        "t={:.2}: E cpu {:.10e} gpu {:.10e} (rel {:.1e}); Z cpu {:.10e} gpu {:.10e} (rel {:.1e}); div gpu {:.1e}",
        dc.t, dc.energy, dg.energy, rel(dc.energy, dg.energy),
        dc.enstrophy, dg.enstrophy, rel(dc.enstrophy, dg.enstrophy), dg.max_div
    );
    // Bounds lock in the unbiased expm1/fma integrating factor (a plain
    // f32 exp factor gave 1.4e-5 here).
    assert!(rel(dc.energy, dg.energy) < 3e-6);
    assert!(rel(dc.enstrophy, dg.enstrophy) < 1e-6);
    assert!(dg.max_div < 1e-5);
    for c in 0..3 {
        let err = max_abs_diff(&Backend::velocity_phys(&mut cpu, c), &gpu.velocity_phys(c));
        eprintln!("field component {c}: max |cpu - gpu| = {err:.2e}");
        assert!(err < 3e-6, "component {c} err {err:e}");
    }
}

#[test]
fn gpu_imex_implicit_sphere_matches_cpu() {
    use spectral_core::Penalization;
    let two_pi = 2.0 * std::f64::consts::PI;
    let (n, nu, d, eta, dt, steps) = (32, 1.0, two_pi / 8.0, 1e-4, 2.0, 15);
    let g = Grid::new(n);
    let pen = Penalization::sphere_sharp(&g, [0.5 * two_pi; 3], d, eta).implicit();
    let bf = [4e-5, 0.0, 0.0];
    let zero = || init::from_fn(&g, |_, _, _| [0.0; 3]);
    let mut cpu = Solver::new(g.clone(), nu, zero());
    cpu.penalization = Some(pen.clone());
    cpu.body_force = Some(bf);
    let mut gpu = CudaSolver::new(g.clone(), nu, zero()).expect("cuda init");
    gpu.set_penalization(&pen).unwrap();
    gpu.body_force = Some(bf);
    for _ in 0..steps {
        cpu.step(dt);
        gpu.step(dt);
    }
    let (uc, ug) = (cpu.mean_velocity()[0], gpu.try_mean_velocity().unwrap()[0]);
    let (fc, fg) = (cpu.solid_force().unwrap()[0], gpu.try_solid_force().unwrap().unwrap()[0]);
    let (it, res) = gpu.last_cg().unwrap();
    eprintln!(
        "IMEX t={}: U cpu {uc:.8e} gpu {ug:.8e} (rel {:.1e}); F cpu {fc:.8e} gpu {fg:.8e} (rel {:.1e}); gpu CG {it:?} {res:?}",
        cpu.t, ((uc - ug) / uc).abs(), ((fc - fg) / fc).abs()
    );
    assert!(((uc - ug) / uc).abs() < 1e-4);
    assert!(((fc - fg) / fc).abs() < 1e-3);
}

#[test]
fn gpu_matches_cpu_with_phase_shift_dealiasing() {
    let (n, re, dt, steps) = (32, 1600.0, 0.01, 200);
    let g = Grid::new_spherical(n);
    let mut cpu = Solver::new(g.clone(), 1.0 / re, init::taylor_green(&g));
    let mut gpu = CudaSolver::new(g.clone(), 1.0 / re, init::taylor_green(&g)).expect("cuda init");
    for _ in 0..steps {
        Backend::step(&mut cpu, dt);
        gpu.step(dt);
    }
    let (dc, dg) = (Backend::diagnostics(&mut cpu), gpu.diagnostics());
    let rel = |a: f64, b: f64| ((a - b) / a).abs();
    eprintln!(
        "spherical t={:.2}: E rel {:.1e}, Z rel {:.1e}",
        dc.t, rel(dc.energy, dg.energy), rel(dc.enstrophy, dg.enstrophy)
    );
    assert!(rel(dc.energy, dg.energy) < 3e-6);
    assert!(rel(dc.enstrophy, dg.enstrophy) < 3e-6);
    for c in 0..3 {
        let err = max_abs_diff(&Backend::velocity_phys(&mut cpu, c), &gpu.velocity_phys(c));
        assert!(err < 3e-6, "component {c} err {err:e}");
    }
}

#[test]
fn gpu_surface_walls_match_cpu() {
    use spectral_core::surface::resolution;
    use spectral_core::Surface;
    let pi = std::f64::consts::PI;
    let (n, nu, dt, steps) = (32, 1.0, 5.0, 10);
    let g = Grid::new(n);
    let surface = Surface::sphere([pi; 3], pi / 3.0, resolution(&g));
    let bf = [5e-5, 0.0, 0.0];
    let zero = || init::from_fn(&g, |_, _, _| [0.0; 3]);
    let mut cpu = Solver::new(g.clone(), nu, zero());
    cpu.surface = Some(surface.clone());
    cpu.body_force = Some(bf);
    let mut gpu = CudaSolver::new(g.clone(), nu, zero()).expect("cuda init");
    gpu.set_surface(&surface).unwrap();
    gpu.body_force = Some(bf);
    for _ in 0..steps {
        cpu.step(dt);
        gpu.step(dt);
    }
    let (uc, ug) = (cpu.mean_velocity()[0], gpu.try_mean_velocity().unwrap()[0]);
    let (fc, fg) = (cpu.solid_force().unwrap()[0], gpu.try_solid_force().unwrap().unwrap()[0]);
    eprintln!(
        "surface t={}: U cpu {uc:.8e} gpu {ug:.8e} (rel {:.1e}); F cpu {fc:.8e} gpu {fg:.8e} (rel {:.1e})",
        cpu.t, ((uc - ug) / uc).abs(), ((fc - fg) / fc).abs()
    );
    assert!(((uc - ug) / uc).abs() < 1e-4);
    assert!(((fc - fg) / fc).abs() < 1e-3);
}

#[test]
fn gpu_penalization_target_field_matches_cpu() {
    use spectral_core::Penalization;
    let pi = std::f64::consts::PI;
    let (n, nu, dt, steps) = (32, 0.05, 0.01, 50);
    let g = Grid::new(n);
    let mut pen = Penalization::sphere(&g, [pi; 3], 2.0, g.dx(), 0.05);
    pen.u_target = Some(init::from_fn(&g, |x, y, z| [0.3 * y.sin(), 0.2 * z.cos(), 0.1 * x.sin()]));
    let mut cpu = Solver::new(g.clone(), nu, init::taylor_green(&g));
    cpu.penalization = Some(pen.clone());
    let mut gpu = CudaSolver::new(g.clone(), nu, init::taylor_green(&g)).expect("cuda init");
    gpu.set_penalization(&pen).unwrap();
    for _ in 0..steps {
        cpu.step(dt);
        gpu.step(dt);
    }
    for c in 0..3 {
        let err = max_abs_diff(&Backend::velocity_phys(&mut cpu, c), &gpu.velocity_phys(c));
        eprintln!("target-field penalization component {c}: max |cpu - gpu| = {err:.2e}");
        assert!(err < 3e-6);
    }
}

#[test]
fn compound_surface_part_forces_add_up() {
    use spectral_core::surface::{resolution, Extension};
    use spectral_core::Surface;
    let pi = std::f64::consts::PI;
    let (n, nu, dt, steps) = (48, 1.0, 5.0, 8);
    let g = Grid::new(n);
    let make = |c: [f64; 3]| {
        let mut s = Surface::sphere(c, 0.9, 0.5 * resolution(&g));
        s.extension = Some(Extension { width: 0.5, sigma_rel: 1.0 });
        s.with_basis(resolution(&g))
    };
    let both = Surface::union(vec![make([pi - 1.4, pi, pi]), make([pi + 1.4, pi, pi])]);
    let bf = [5e-5, 0.0, 0.0];
    let zero = || init::from_fn(&g, |_, _, _| [0.0; 3]);
    let mut cpu = Solver::new(g.clone(), nu, zero());
    cpu.surface = Some(both.clone());
    cpu.body_force = Some(bf);
    let mut gpu = CudaSolver::new(g.clone(), nu, zero()).expect("cuda init");
    gpu.set_surface(&both).unwrap();
    gpu.body_force = Some(bf);
    for _ in 0..steps {
        cpu.step(dt);
        gpu.step(dt);
    }
    let (pc, pg) = (cpu.surface_part_forces().unwrap(), gpu.surface_part_forces().unwrap());
    let (tc, tg) = (cpu.solid_force().unwrap()[0], gpu.try_solid_force().unwrap().unwrap()[0]);
    eprintln!("parts cpu {:.6e} {:.6e} (sum {:.6e}, total {tc:.6e}); gpu {:.6e} {:.6e} (total {tg:.6e})", pc[0][0], pc[1][0], pc[0][0] + pc[1][0], pg[0][0], pg[1][0]);
    assert!(((pc[0][0] + pc[1][0]) / tc - 1.0).abs() < 1e-9);
    assert!(((pg[0][0] + pg[1][0]) / tg - 1.0).abs() < 1e-3);
    // mirror-symmetric pair: equal drag on both
    assert!((pc[0][0] / pc[1][0] - 1.0).abs() < 1e-3);
}

#[test]
fn gpu_certificate_matches_cpu() {
    use spectral_core::certify::Certifier;
    use spectral_cuda::CudaCertifier;
    let (n, re, dt, steps) = (24, 200.0, 0.01, 150);
    let g = Grid::new(n);
    let mut cpu = Solver::new(g.clone(), 1.0 / re, init::taylor_green(&g));
    let mut gpu = CudaSolver::new(g.clone(), 1.0 / re, init::taylor_green(&g)).expect("cuda init");
    let mut cc = Certifier::new(&g, 0.0);
    let mut gc = CudaCertifier::new(&gpu, 0.0).unwrap();
    for k in 1..=steps {
        let a = cpu.step_certified(dt, &mut cc);
        let b = gpu.step_certified(dt, &mut gc).unwrap();
        if k % 50 == 0 {
            eprintln!(
                "t={:.2}: bound cpu {:.4e} gpu {:.4e} | trunc {:.3e}/{:.3e} | time {:.1e}/{:.1e} | mu {:.4}/{:.4}",
                a.t, a.bound, b.bound, a.resid_trunc, b.resid_trunc, a.resid_time, b.resid_time, a.mu, b.mu
            );
        }
    }
    let (a, b) = (cc.history.last().unwrap(), gc.history.last().unwrap());
    assert!((a.bound / b.bound - 1.0).abs() < 1e-3, "bounds {} vs {}", a.bound, b.bound);
    assert!((a.mu / b.mu - 1.0).abs() < 1e-4);
}
