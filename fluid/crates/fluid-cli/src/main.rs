//! fluid: command-line driver.
//!
//!   fluid tg [--backend cpu|cuda] [--spherical 0|1] [--n 64] [--re 1600] [--t-end 10] [--cfl 0.5] [--every 0.1] [--out tg.csv]
//!
//!   fluid sphere [--backend cpu|cuda] [--n 128] [--re 100] [--d 0.785] [--eta 2e-3] [--u 1] [--sponge 1] [--delta 1] [--implicit 0] [--eta-c C] [--walls penal|point|ext] [--ext-width W]
//!                [--t-end 40] [--cfl 0.5] [--every 0.25] [--out sphere.csv]
//!
//!   fluid tg --certify  [--backend cpu|cuda] [--n 32] [--re 400] [--t-end 4]  (certified error bound)
//!   fluid tg --adaptive TARGET [--n 16] [--re 400] [--t-end 4] (certified adaptive resolution)
//!
//!   fluid case path/to/case.json   (engineering run: CAD geometry, inlets/outlets/walls, forces, VTK)
//!
//! `tg` runs the 3D Taylor–Green vortex and writes t, E, Z, ε = 2νZ, -dE/dt,
//! the divergence check and the spectral-tail energy fraction to CSV.
//! `sphere` runs uniform flow (U∞ = 1 along x) past a penalized sphere with a
//! downstream sponge and writes the drag/lift coefficients.

use std::{env, fs::File, io::Write, time::Instant};

use spectral_core::certify::{AdaptPolicy, Adaptive, Certifier};
use spectral_core::surface::{resolution, Extension};
use spectral_core::epidemic::{lattice, seed_virions, Epidemic};
use spectral_core::virus::{Cloud, Medium, Scaling, Virion};
use spectral_core::{init, Backend, Grid, Penalization, Solver, Surface};

fn arg<T: std::str::FromStr>(args: &[String], name: &str, default: T) -> T {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .map(|v| v.parse().unwrap_or_else(|_| panic!("bad value for {name}")))
        .unwrap_or(default)
}

fn taylor_green(args: &[String]) {
    let n: usize = arg(args, "--n", 64);
    let re: f64 = arg(args, "--re", 1600.0);
    let t_end: f64 = arg(args, "--t-end", 10.0);
    let cfl: f64 = arg(args, "--cfl", 0.5);
    let every: f64 = arg(args, "--every", 0.1);
    let out: String = arg(args, "--out", format!("tg_n{n}_re{re}.csv"));

    let backend: String = arg(args, "--backend", "cpu".to_string());
    // --spherical 1: phase-shift dealiasing with spherical truncation (41% larger k_max)
    let spherical: u8 = arg(args, "--spherical", 0);
    let g = if spherical != 0 { Grid::new_spherical(n) } else { Grid::new(n) };
    let u0 = init::taylor_green(&g);
    let mut s = make_backend(&backend, g, 1.0 / re, u0);
    let mut f = File::create(&out).expect("create output");
    writeln!(f, "t,energy,enstrophy,dissipation,minus_dEdt,max_div,tail_fraction").unwrap();

    let start = Instant::now();
    let mut prev = s.diagnostics();
    let mut steps = 0usize;
    let mut next_out = 0.0;
    eprintln!("Taylor–Green [{}] N={n}^3 Re={re} -> t={t_end}, output {out}", s.name());
    while s.time() < t_end - 1e-12 {
        let t = s.time();
        let dt = s.cfl_dt(cfl).min(next_out + every - t).min(t_end - t);
        s.step(dt);
        steps += 1;
        let t = s.time();
        if t >= next_out + every - 1e-12 || t >= t_end - 1e-12 {
            let d = s.diagnostics();
            let dedt = -(d.energy - prev.energy) / (d.t - prev.t);
            writeln!(
                f,
                "{:.6},{:.12e},{:.12e},{:.12e},{:.12e},{:.3e},{:.3e}",
                d.t, d.energy, d.enstrophy, d.dissipation, dedt, d.max_div, d.tail_energy_fraction
            )
            .unwrap();
            eprintln!(
                "t={:7.3} E={:.6e} eps={:.6e} tail={:.1e} steps={} {:.1}s",
                d.t,
                d.energy,
                d.dissipation,
                d.tail_energy_fraction,
                steps,
                start.elapsed().as_secs_f64()
            );
            prev = d;
            next_out = t;
        }
    }
    let secs = start.elapsed().as_secs_f64();
    eprintln!(
        "done: {steps} steps in {secs:.1}s ({:.3} s/step, {:.1} ns per point-step)",
        secs / steps as f64,
        secs / steps as f64 / (n * n * n) as f64 * 1e9
    );
}

/// Certified run on the CPU reference solver: fixed N (`target` = None) or
/// certificate-driven adaptive resolution toward a target L² error at t_end.
fn taylor_green_certified(args: &[String], target: Option<f64>) {
    let backend: String = arg(args, "--backend", "cpu".to_string());
    if backend == "cuda" {
        assert!(target.is_none(), "certified adaptive resolution is CPU-only");
        #[cfg(feature = "cuda")]
        return taylor_green_certified_cuda(args);
        #[cfg(not(feature = "cuda"))]
        panic!("CUDA support requires rebuilding with --features cuda");
    }
    let n: usize = arg(args, "--n", 32);
    let re: f64 = arg(args, "--re", 400.0);
    let t_end: f64 = arg(args, "--t-end", 4.0);
    let cfl: f64 = arg(args, "--cfl", 0.5);
    let every: f64 = arg(args, "--every", 0.1);
    let out: String = arg(args, "--out", format!("tg_cert_n{n}_re{re}.csv"));
    let g = Grid::new(n);
    let mut s = Solver::new(g.clone(), 1.0 / re, init::taylor_green(&g));
    let mut cert = Certifier::new(&g, 0.0);
    let mut ad = target.map(|tg| Adaptive::new(AdaptPolicy::new(tg, t_end)));
    let mut f = File::create(&out).expect("create output");
    writeln!(f, "t,n,energy,bound,resid_trunc,resid_time,mu").unwrap();
    let start = Instant::now();
    let (mut steps, mut point_steps, mut next_out) = (0usize, 0.0f64, 0.0);
    while s.t < t_end - 1e-12 {
        let dt = s.cfl_dt(cfl).min(t_end - s.t);
        let st = s.step_certified(dt, &mut cert);
        steps += 1;
        point_steps += (s.grid.n as f64).powi(3);
        if let Some(a) = ad.as_mut() {
            if let Some(m) = a.decide(&s, &cert, &st, dt) {
                let (ns, dropped) = s.regrid(m);
                s = ns;
                cert.rebind(&s.grid, dropped);
                a.changes.push((s.t, m));
                eprintln!("t={:.3}: regrid -> {m}^3 (dropped {dropped:.1e})", s.t);
            }
        }
        if s.t >= next_out + every - 1e-12 || s.t >= t_end - 1e-12 {
            let d = s.diagnostics();
            writeln!(f, "{:.6},{},{:.12e},{:.6e},{:.6e},{:.6e},{:.6e}", s.t, s.grid.n, d.energy, cert.y.sqrt(), st.resid_trunc, st.resid_time, st.mu).unwrap();
            next_out = s.t;
        }
    }
    let secs = start.elapsed().as_secs_f64();
    eprintln!(
        "done: t={t_end} certified ||u - v|| <= {:.3e} (box-mean L2), final N={}, {steps} steps, {:.3e} point-steps, {secs:.1}s",
        cert.y.sqrt(), s.grid.n, point_steps
    );
}

/// Certified fixed-resolution run on the GPU.
#[cfg(feature = "cuda")]
fn taylor_green_certified_cuda(args: &[String]) {
    let n: usize = arg(args, "--n", 64);
    let re: f64 = arg(args, "--re", 400.0);
    let t_end: f64 = arg(args, "--t-end", 2.5);
    let cfl: f64 = arg(args, "--cfl", 0.5);
    let every: f64 = arg(args, "--every", 0.1);
    let out: String = arg(args, "--out", format!("tg_cert_cuda_n{n}_re{re}.csv"));
    let g = Grid::new(n);
    let mut s = spectral_cuda::CudaSolver::new(g.clone(), 1.0 / re, init::taylor_green(&g)).expect("CUDA init");
    let mut cert = spectral_cuda::CudaCertifier::new(&s, 0.0).expect("certifier");
    let mut f = File::create(&out).expect("create output");
    writeln!(f, "t,n,energy,bound,resid_trunc,resid_time,mu").unwrap();
    let start = Instant::now();
    let (mut steps, mut next_out) = (0usize, 0.0);
    while s.t < t_end - 1e-12 {
        let dt = Backend::cfl_dt(&mut s, cfl).min(t_end - s.t);
        let st = s.step_certified(dt, &mut cert).expect("certified step");
        steps += 1;
        if s.t >= next_out + every - 1e-12 || s.t >= t_end - 1e-12 {
            let d = Backend::diagnostics(&mut s);
            writeln!(f, "{:.6},{n},{:.12e},{:.6e},{:.6e},{:.6e},{:.6e}", s.t, d.energy, st.bound, st.resid_trunc, st.resid_time, st.mu).unwrap();
            next_out = s.t;
        }
    }
    eprintln!(
        "done [cuda]: t={t_end} certified ||u - v|| <= {:.3e} (box-mean L2), N={n}, {steps} steps, {:.1}s",
        cert.y.sqrt(),
        start.elapsed().as_secs_f64()
    );
}

fn make_backend(name: &str, g: Grid, nu: f64, u0: [Vec<f64>; 3]) -> Box<dyn Backend> {
    match name {
        "cpu" => Box::new(Solver::new(g, nu, u0)),
        #[cfg(feature = "cuda")]
        "cuda" => Box::new(spectral_cuda::CudaSolver::new(g, nu, u0).expect("CUDA backend init")),
        #[cfg(not(feature = "cuda"))]
        "cuda" => panic!("CUDA support requires rebuilding with --features cuda"),
        other => panic!("unknown backend {other}"),
    }
}

fn sphere(args: &[String]) {
    let n: usize = arg(args, "--n", 128);
    let re: f64 = arg(args, "--re", 100.0);
    let two_pi = 2.0 * std::f64::consts::PI;
    let d: f64 = arg(args, "--d", two_pi / 8.0);
    let eta: f64 = arg(args, "--eta", 2e-3);
    let t_end: f64 = arg(args, "--t-end", 40.0);
    let cfl: f64 = arg(args, "--cfl", 0.5);
    let every: f64 = arg(args, "--every", 0.25);
    let backend: String = arg(args, "--backend", "cpu".to_string());
    let out: String = arg(args, "--out", format!("sphere_n{n}_re{re}.csv"));
    let sponge: u8 = arg(args, "--sponge", 1);
    let delta_cells: f64 = arg(args, "--delta", 1.0);
    // Free-stream speed; small U with moderate ν keeps the penalization wall
    // layer sqrt(ν η) thin at low Reynolds number.
    let speed: f64 = arg(args, "--u", 1.0);
    let u_inf = [speed, 0.0, 0.0];
    let nu = speed * d / re;

    let g = Grid::new(n);
    // Grid-matched penalization: η = c·dx²/ν puts the Brinkman layer √(νη) at
    // √c·dx. With the spectral discretization this converges monotonically,
    // whereas η → 0 converges to an over-stiffened body (see
    // analysis/stokes_array_check.md). Overrides --eta when given.
    let eta_c: f64 = arg(args, "--eta-c", f64::NAN);
    let eta = if eta_c.is_nan() { eta } else { eta_c * g.dx() * g.dx() / nu };
    let center = [0.3 * two_pi, 0.5 * two_pi, 0.5 * two_pi];
    // --walls penal (volume penalization, default) | point | ext
    let walls: String = arg(args, "--walls", "penal".to_string());
    let implicit: u8 = arg(args, "--implicit", 0);
    let mut pen = if walls != "penal" {
        // surface walls carry the body; penalization only provides the sponge
        Penalization::sphere(&g, center, 0.0, g.dx(), eta)
    } else if implicit != 0 {
        Penalization::sphere_sharp(&g, center, d, eta).implicit()
    } else {
        Penalization::sphere(&g, center, d, delta_cells * g.dx(), eta)
    };
    if sponge != 0 {
        pen = pen.with_sponge(&g, 0.8 * two_pi, 0.15 * two_pi, 0.05, u_inf);
    }
    let dt_pen = pen.max_dt();
    let u0 = init::from_fn(&g, |_, _, _| [0.0; 3]);
    let surface = match walls.as_str() {
        "point" => Some(Surface::sphere(center, 0.5 * d, resolution(&g))),
        "ext" => {
            let width: f64 = arg(args, "--ext-width", 0.6 * 0.5 * d);
            if width / g.dx() < 5.0 {
                // the smooth extension needs a resolved shell (≳ 5 cells);
                // thinner shells do not converge (analysis/surface_walls.md)
                eprintln!(
                    "warning: extension shell is {:.1} cells (< 5); falling back to point-force walls",
                    width / g.dx()
                );
                Some(Surface::sphere(center, 0.5 * d, resolution(&g)))
            } else {
                let mut sf = Surface::sphere(center, 0.5 * d, 0.5 * resolution(&g));
                sf.extension = Some(Extension { width, sigma_rel: 1.0 });
                Some(sf.with_basis(resolution(&g)))
            }
        }
        _ => None,
    };
    let mut s = make_backend(&backend, g, nu, u0);
    s.set_penalization(pen);
    if let Some(sf) = surface {
        s.set_surface(sf);
    }
    s.set_mean_flow(u_inf);

    let area = std::f64::consts::PI * d * d / 4.0;
    let q = 0.5 * u_inf[0] * u_inf[0] * area;
    eprintln!(
        "sphere [{}] N={n}^3 Re={re} D={d:.4} ({:.1} cells) eta={eta:e} (wall layer sqrt(nu*eta)={:.2e} = {:.2}% of D) dt_pen={dt_pen:.2e}",
        s.name(),
        d / (two_pi / n as f64),
        (nu * eta).sqrt(),
        100.0 * (nu * eta).sqrt() / d
    );
    let mut f = File::create(&out).expect("create output");
    writeln!(f, "t,cd,cl_y,cl_z,energy,max_div").unwrap();
    let start = Instant::now();
    let (mut steps, mut next_out) = (0usize, 0.0);
    while s.time() < t_end - 1e-12 {
        let t = s.time();
        let dt = s.cfl_dt(cfl).min(dt_pen).min(next_out + every - t).min(t_end - t);
        s.step(dt);
        steps += 1;
        let t = s.time();
        if t >= next_out + every - 1e-12 || t >= t_end - 1e-12 {
            let fo = s.solid_force().unwrap();
            let dg = s.diagnostics();
            let (cd, cy, cz) = (fo[0] / q, fo[1] / q, fo[2] / q);
            writeln!(f, "{t:.6},{cd:.8e},{cy:.8e},{cz:.8e},{:.10e},{:.3e}", dg.energy, dg.max_div).unwrap();
            eprintln!("t={t:7.3} Cd={cd:.5} Cl=({cy:+.2e},{cz:+.2e}) steps={steps} {:.1}s", start.elapsed().as_secs_f64());
            next_out = t;
        }
    }
    eprintln!("done: {steps} steps in {:.1}s", start.elapsed().as_secs_f64());
}

/// Virus-particle (virion) transport: point particles advected by the flow
/// with Brownian diffusion (Stokes–Einstein) and gravitational settling
/// (Stokes + Cunningham slip). A virion is ~100 nm ≪ any grid cell with a
/// ~1e-7 s momentum relaxation time, so the overdamped Langevin model is
/// exact; see spectral_core::virus for the physics and validation tests.
///
///   fluid virus [--medium air|water] [--radius 60e-9] [--np 20000]
///               [--flow uniform|tg|none] [--vel 0.05] [--box-len 6.2832e-3]
///               [--n N] [--t-end 5] [--every 0.25] [--seed S] [--out FILE]
///
/// --vel maps solver velocity 1 to m/s; --box-len is the physical box side
/// (the 2π box) in metres. Writes t, mean displacement, MSD and the exact
/// theory values (drift U t, settling −v_s t, Einstein 6 D t) to CSV.
fn virus(args: &[String]) {
    let medium = match arg::<String>(args, "--medium", "air".to_string()).as_str() {
        "air" => Medium::air(),
        "water" => Medium::water(),
        other => panic!("unknown medium {other} (air|water)"),
    };
    let mut v = Virion::influenza();
    v.radius = arg(args, "--radius", 60e-9);
    let np: usize = arg(args, "--np", 20_000);
    let flow: String = arg(args, "--flow", "uniform".to_string());
    let sc = Scaling {
        box_len: arg(args, "--box-len", 2.0 * std::f64::consts::PI * 1e-3),
        velocity: arg(args, "--vel", 0.05),
    };
    let t_end_s: f64 = arg(args, "--t-end", 5.0); // seconds
    let every_s: f64 = arg(args, "--every", 0.25); // seconds
    let seed: u64 = arg(args, "--seed", 0x5EED);
    let out: String = arg(args, "--out", format!("virus_{}.csv", medium.name));

    let d_si = v.diffusivity(&medium);
    let vs_si = v.settling_speed(&medium);
    let nu = sc.nu_solver(&medium);
    let d = sc.diffusivity_solver(d_si);
    let vs = sc.settling_solver(vs_si);
    let t_end = t_end_s / sc.t0();
    let every = every_s / sc.t0();
    eprintln!(
        "virus [{}]: r={:.0} nm, D={:.3e} m²/s (C_c={:.2}), v_s={:.3e} m/s, τ={:.1e} s; Pe={:.1}, box={:.2} mm, U={:.3} m/s, Re={:.2}",
        medium.name,
        v.radius * 1e9,
        d_si,
        v.cunningham(&medium),
        vs_si,
        v.relaxation_time(&medium),
        sc.velocity * sc.l0() / d_si,
        sc.box_len * 1e3,
        sc.velocity,
        sc.velocity * sc.l0() / medium.nu(),
    );

    // Background flow: `none` and `uniform` are exact on the coarsest grid;
    // `tg` is the decaying Taylor–Green vortex at the Re printed above.
    let n: usize = arg(args, "--n", if flow == "tg" { 32 } else { 8 });
    let g = Grid::new(n);
    let u0 = match flow.as_str() {
        "tg" => init::taylor_green(&g),
        "uniform" | "none" => init::from_fn(&g, |_, _, _| [0.0; 3]),
        other => panic!("unknown flow {other} (uniform|tg|none)"),
    };
    let mut s = Solver::new(g, nu, u0);
    if flow == "uniform" {
        s.set_mean_flow([1.0, 0.0, 0.0]); // held exactly, = U x̂ in SI
    }

    let mut cloud = Cloud::new(np, d, vs, seed);
    let mut f = File::create(&out).expect("create output");
    writeln!(f, "t_s,mean_x_m,mean_y_m,mean_z_m,msd_m2,exact_6Dt_m2,exact_drift_x_m,exact_settle_z_m").unwrap();

    let start = Instant::now();
    let (mut steps, mut next_out) = (0usize, 0.0);
    while s.time() < t_end - 1e-12 {
        let t = s.time();
        let dt = s
            .cfl_dt(0.5)
            .min(every / 4.0)
            .min(next_out + every - t)
            .min(t_end - t);
        s.step(dt);
        let u: [Vec<f64>; 3] = std::array::from_fn(|c| s.velocity_phys(c));
        cloud.step(Some(&u), n, dt);
        steps += 1;
        let t = s.time();
        if t >= next_out + every - 1e-12 || t >= t_end - 1e-12 {
            let st = cloud.stats();
            let (l0, t0) = (sc.l0(), sc.t0());
            let ts = t * t0;
            let drift = if flow == "uniform" { sc.velocity * ts } else { 0.0 };
            writeln!(
                f,
                "{ts:.6},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{drift:.6e},{:.6e}",
                st.mean[0] * l0,
                st.mean[1] * l0,
                st.mean[2] * l0,
                st.msd_about_mean * l0 * l0,
                6.0 * d_si * ts,
                -vs_si * ts
            )
            .unwrap();
            eprintln!(
                "t={ts:7.3}s drift_x={:+.3e} m settle_z={:+.3e} m MSD={:.3e} m² (6Dt={:.3e}) steps={} {:.1}s",
                st.mean[0] * l0,
                st.mean[2] * l0,
                st.msd_about_mean * l0 * l0,
                6.0 * d_si * ts,
                steps,
                start.elapsed().as_secs_f64()
            );
            next_out = t;
        }
    }

    // Final accuracy summary against the exact theory.
    let st = cloud.stats();
    eprintln!("done: {steps} steps in {:.1}s, {np} particles -> {out}", start.elapsed().as_secs_f64());
    if flow == "tg" {
        // In a non-uniform flow the MSD is dominated by advective dispersion
        // (particles trapped in vortex cells), so 6Dt is only a floor.
        eprintln!("note: MSD includes flow dispersion; 6Dt = {:.3e} m² is the diffusion-only floor", 6.0 * d_si * t_end * sc.t0());
    } else {
        let msd_rel = (st.msd_about_mean - 6.0 * d * t_end) / (6.0 * d * t_end);
        eprintln!("accuracy: MSD {:+.2}% off 6Dt (statistical σ ≈ {:.2}%)", msd_rel * 100.0, (2.0 / (3.0 * np as f64)).sqrt() * 100.0);
        if flow == "uniform" {
            let err = (st.mean[0] - t_end).abs() / t_end;
            eprintln!("accuracy: x-drift {err:.2e} relative to U t");
        }
    }
}

/// Spatially-resolved virus spread + reproduction, coupled to the fluid
/// solver for advection. Virions are Brownian particles that infect a lattice
/// of host cells; infected cells produce new virions and die. The infection
/// front is an emergent boundary: only the sparse active region (cells with
/// virions/infected cells) is simulated, so a large environment runs cheaply.
///
///   fluid epidemic [--cells N] [--np V0] [--d D] [--vs VS]
///                 [--beta B] [--p P] [--c C] [--delta DELTA] [--r-bind R]
///                 [--flow none|uniform|tg] [--vel U] [--box-len L]
///                 [--t-end T] [--every DT] [--seed S] [--naive] [--out FILE]
///
/// Rates are in solver units (box [0, 2π)³). --naive runs the O(V·S) full-
/// domain infection search for comparison; the default spatial hash is exact
/// (same outcome, same seed) but only touches the active region.
fn epidemic(args: &[String]) {
    let n: usize = arg(args, "--cells", 16);
    let np: usize = arg(args, "--np", 200);
    let d: f64 = arg(args, "--d", 0.01);
    let vs: f64 = arg(args, "--vs", 0.0);
    let beta: f64 = arg(args, "--beta", 1e-4);
    let p: f64 = arg(args, "--p", 50.0);
    let c: f64 = arg(args, "--c", 0.1);
    let delta: f64 = arg(args, "--delta", 0.2);
    let r_bind: f64 = arg(args, "--r-bind", 0.5);
    let flow: String = arg(args, "--flow", "none".to_string());
    let t_end: f64 = arg(args, "--t-end", 20.0);
    let every: f64 = arg(args, "--every", 1.0);
    let seed: u64 = arg(args, "--seed", 0x5EED);
    let naive: bool = args.iter().any(|a| a == "--naive");
    let out: String = arg(args, "--out", "epidemic.csv".to_string());

    let vol = (2.0 * std::f64::consts::PI).powi(3);
    let cells = lattice(n);
    let virions = seed_virions(np, seed);
    let h = r_bind;
    let mut e = Epidemic::new(cells, virions, beta * vol, p, c, delta, r_bind, d, vs, vol, h, naive, seed ^ 0x1234);

    // Background flow (optional): uniform is held exactly; tg is the decaying
    // Taylor–Green vortex. `none` is quiescent.
    let grid_n: usize = arg(args, "--n", if flow == "tg" { 32 } else { 8 });
    let g = Grid::new(grid_n);
    let u0 = match flow.as_str() {
        "tg" => init::taylor_green(&g),
        "uniform" | "none" => init::from_fn(&g, |_, _, _| [0.0; 3]),
        other => panic!("unknown flow {other} (none|uniform|tg)"),
    };
    let mut s = Solver::new(g, 0.1, u0);
    if flow == "uniform" {
        s.set_mean_flow([1.0, 0.0, 0.0]);
    }

    let mut f = File::create(&out).expect("create output");
    writeln!(f, "t,S,I,D,V,active_cells,total_cells,active_fraction").unwrap();
    let start = Instant::now();
    let (mut steps, mut next_out) = (0usize, 0.0);
    while e.t < t_end - 1e-12 {
        let dt = s.cfl_dt(0.5).min(every / 4.0).min(next_out + every - e.t).min(t_end - e.t);
        s.step(dt);
        let u: [Vec<f64>; 3] = std::array::from_fn(|c| s.velocity_phys(c));
        e.step(Some(&u), grid_n, dt);
        steps += 1;
        if e.t >= next_out + every - 1e-12 || e.t >= t_end - 1e-12 {
            let (s_, i, d_, v) = e.stats();
            let af = e.active_fraction();
            writeln!(f, "{:.4},{s_},{i},{d_},{v},{},{},{af:.6}", e.t, e.active_cells, e.total_cells).unwrap();
            eprintln!("t={:6.2} S={s_:5} I={i:4} D={d_:4} V={v:6} active={:.1}% (speedup {:.0}x) steps={steps} {:.1}s", e.t, af * 100.0, 1.0 / af.max(1e-9), start.elapsed().as_secs_f64());
            next_out = e.t;
        }
    }
    let (s_, i, d_, v) = e.stats();
    let af = e.active_fraction();
    eprintln!("done: {steps} steps in {:.1}s -> {out}", start.elapsed().as_secs_f64());
    eprintln!("final: S={s_} I={i} D={d_} V={v}; active region {:.2}% of domain (emergent-boundary speedup {:.0}x)", af * 100.0, 1.0 / af.max(1e-9));
}

fn main() {
    let args: Vec<String> = env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("tg") if args.iter().any(|a| a == "--adaptive") => {
            let target: f64 = arg(&args, "--adaptive", 1e-3);
            taylor_green_certified(&args[2..], Some(target))
        }
        Some("tg") if args.iter().any(|a| a == "--certify") => taylor_green_certified(&args[2..], None),
        Some("tg") => taylor_green(&args[2..]),
        Some("case") => run_case(&args[2..]),
        Some("cad") => run_cad(&args[2..]),
        Some("project") => std::process::exit(project::cli::run(&args[2..])),
        Some("sphere") => sphere(&args[2..]),
        Some("virus") => virus(&args[2..]),
        Some("epidemic") => epidemic(&args[2..]),
        _ => {
            eprintln!("usage: fluid tg [--backend cpu|cuda] [--n N] [--re RE] [--t-end T] [--cfl C] [--every DT] [--out FILE]");
            eprintln!("       fluid sphere [--backend cpu|cuda] [--n N] [--re RE] [--d D] [--eta ETA] [--t-end T] [--cfl C] [--every DT] [--out FILE]");
            eprintln!("       fluid virus [--medium air|water] [--radius R] [--np NP] [--flow uniform|tg|none] [--vel U] [--box-len L] [--t-end T] [--out FILE]");
            eprintln!("       fluid epidemic [--cells N] [--np V0] [--d D] [--beta B] [--p P] [--c C] [--delta DELTA] [--r-bind R] [--flow none|uniform|tg] [--t-end T] [--out FILE]");
            eprintln!("       fluid case path/to/case.json");
            eprintln!("       fluid cad path/to/model.json [--cell H] [--stl F] [--vtk F] [--set name=value ...]");
            eprintln!("       fluid project report|check|diff|set|fmt|export|kinds ...   (fluid project for details)");
            std::process::exit(2);
        }
    }
}

/// Build a parametric CAD model: report it, and write STL and VTK.
///
/// `--set name=value` overrides a parameter, so a design sweep is a shell
/// loop. `--cell` is the meshing resolution; without it the model's own
/// `cell` is used, and failing that one two-hundredth of the bounding-box
/// diagonal.
fn run_cad(args: &[String]) {
    let Some(path) = args.first() else {
        eprintln!("usage: fluid cad path/to/model.json [--cell H] [--stl F] [--vtk F] [--set name=value ...]");
        std::process::exit(2);
    };
    if path.starts_with("--") {
        eprintln!("usage: fluid cad path/to/model.json [--cell H] [--stl F] [--vtk F] [--set name=value ...]");
        std::process::exit(2);
    }
    let model = cad::model::Model::read(path).unwrap_or_else(|e| {
        eprintln!("{e}");
        std::process::exit(1)
    });
    let base = std::path::Path::new(path).parent().unwrap_or(std::path::Path::new(".")).to_path_buf();

    let mut overrides = cad::Params::new();
    let mut i = 0;
    while i < args.len() {
        if args[i] == "--set" {
            let Some(pair) = args.get(i + 1) else {
                eprintln!("--set needs name=value");
                std::process::exit(2);
            };
            let Some((name, value)) = pair.split_once('=') else {
                eprintln!("--set needs name=value, found {pair:?}");
                std::process::exit(2);
            };
            let v: f64 = value.parse().unwrap_or_else(|_| {
                eprintln!("--set {name}: {value:?} is not a number");
                std::process::exit(2)
            });
            overrides.set(name, v);
            i += 2;
        } else {
            i += 1;
        }
    }

    let params = model.resolve(&overrides).unwrap_or_else(|e| {
        eprintln!("{e}");
        std::process::exit(1)
    });
    let solid = model.build(&overrides, &base).unwrap_or_else(|e| {
        eprintln!("{e}");
        std::process::exit(1)
    });

    let (lo, hi) = solid.bbox();
    if (0..3).any(|k| !lo[k].is_finite() || !hi[k].is_finite()) {
        eprintln!("the model is unbounded; clip it with a \"clip\" feature before meshing");
        std::process::exit(1);
    }
    let diagonal = ((hi[0] - lo[0]).powi(2) + (hi[1] - lo[1]).powi(2) + (hi[2] - lo[2]).powi(2)).sqrt();
    let cell = if args.iter().any(|a| a == "--cell") {
        arg(args, "--cell", diagonal / 200.0)
    } else {
        model.cell.unwrap_or(diagonal / 200.0)
    };

    println!("model {} ({} parameters)", if model.name.is_empty() { path.clone() } else { model.name.clone() }, params.values.len());
    for (k, v) in &params.values {
        println!("  {k} = {v}");
    }
    println!(
        "bounding box [{:.6}, {:.6}, {:.6}] .. [{:.6}, {:.6}, {:.6}], diagonal {diagonal:.6}",
        lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]
    );
    println!("exact signed distance: {}", solid.exact());

    let started = std::time::Instant::now();
    let mesh = cad::mesh_solid(&solid, cell).unwrap_or_else(|e| {
        eprintln!("{e}");
        std::process::exit(1)
    });
    let report = mesh.report();
    println!(
        "meshed at cell {cell:.6} in {:.2} s: {} triangles, {} vertices, volume {:.6e}, area {:.6e}",
        started.elapsed().as_secs_f64(),
        report.triangles,
        report.vertices,
        report.volume,
        report.area
    );
    println!(
        "  watertight {} (boundary {}, non-manifold {}, inconsistent {}, degenerate {})",
        report.watertight,
        report.boundary_edges,
        report.nonmanifold_edges,
        report.inconsistent_edges,
        report.degenerate
    );
    if report.boundary_edges > 0 {
        println!("  the surface is open: refine --cell, or the model has a feature thinner than a cell");
    } else if report.nonmanifold_edges > 0 {
        println!("  closed but pinched: solids in contact, or a gap thinner than a cell");
    }

    let stl = arg::<String>(args, "--stl", String::new());
    if !stl.is_empty() {
        mesh.write_stl(&stl).unwrap_or_else(|e| {
            eprintln!("{e:?}");
            std::process::exit(1)
        });
        println!("wrote {stl}");
    }
    let vtk = arg::<String>(args, "--vtk", String::new());
    if !vtk.is_empty() {
        let grid = cad::Grid3::covering((lo, hi), cell, 2).unwrap_or_else(|e| {
            eprintln!("{e}");
            std::process::exit(1)
        });
        cad::sim::write_vtk_field(&solid, &grid, &vtk).unwrap_or_else(|e| {
            eprintln!("{e}");
            std::process::exit(1)
        });
        println!("wrote {vtk} ({} x {} x {} nodes)", grid.dims[0], grid.dims[1], grid.dims[2]);
    }
}

fn run_case(args: &[String]) {
    let Some(path) = args.first() else {
        eprintln!("usage: fluid case path/to/case.json");
        std::process::exit(2);
    };
    let (case, base) = immersed::case::Case::load(path).unwrap_or_else(|e| {
        eprintln!("{e}");
        std::process::exit(1)
    });
    match case.run(&base, |line| println!("{line}")) {
        Ok(s) => {
            for (b, body) in s.bodies.iter().enumerate() {
                println!(
                    "body {b} ({}, {} triangles, watertight {}): force {:?} (surface {:?}){}",
                    body.file,
                    body.triangles,
                    body.watertight,
                    body.force,
                    body.force_surface,
                    body.coefficients.map(|c| format!(", coefficients {c:?}")).unwrap_or_default()
                );
            }
            println!(
                "{} steps to t = {:.3} ({}converged, rate {:.2e}) in {:.1} s; results in {}",
                s.steps,
                s.time,
                if s.converged { "" } else { "not " },
                s.final_rate,
                s.wall_seconds,
                base.join(&case.output.directory).display()
            );
        }
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
    }
}
