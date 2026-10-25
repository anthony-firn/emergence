//! Surface-constrained walls on two problems with exact answers.
//!
//! usage: wall_probe channel|sphere N alpha [ext_width [sigma_rel [basis_alpha]]]
//!   alpha: marker spacing in units of the resolution length 2π/(2k_max+1)
//!   ext_width: use smooth forcing extension with this shell thickness
//! channel: plane Poiseuille, fluid |y − π| < π/2; exact fluid mean G h²/(3ν)
//! sphere:  periodic sphere array (Hasimoto / Sangani–Acrivos K = 1.2139)
use spectral_core::surface::resolution;
use spectral_core::surface::Extension;
use spectral_core::{init, Grid, Solver, Surface};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let case = a[1].as_str();
    let n: usize = a[2].parse().unwrap();
    let alpha: f64 = a[3].parse().unwrap();
    let pi = std::f64::consts::PI;
    let vol = (2.0 * pi).powi(3);
    let nu = 1.0;
    let g = Grid::new(n);
    let spacing = alpha * resolution(&g);
    let (surface, gx) = match case {
        "channel" => (Surface::channel_y(pi, 0.5 * pi, spacing), 1e-3),
        _ => {
            let r = sphere_radius();
            (Surface::sphere([pi; 3], r, spacing), 6.0 * pi * nu * r * 1e-3 * k_exact(r) / vol)
        }
    };
    let mut surface = surface;
    if let Some(wd) = a.get(4).map(|v| v.parse::<f64>().unwrap()) {
        let sigma_rel = a.get(5).map(|v| v.parse().unwrap()).unwrap_or(1.0);
        surface.extension = Some(Extension { width: wd, sigma_rel });
        if let Some(ab) = a.get(6).map(|v| v.parse::<f64>().unwrap()) {
            surface = surface.with_basis(ab * resolution(&g));
        }
    }
    let nm = surface.points.len();
    let mut s = Solver::new(g.clone(), nu, init::from_fn(&g, |_, _, _| [0.0; 3]));
    s.surface = Some(surface.clone());
    s.body_force = Some([gx, 0.0, 0.0]);
    let t0 = std::time::Instant::now();
    let dt: f64 = std::env::var("DT").ok().map(|v| v.parse().unwrap()).unwrap_or(25.0);
    let steps: usize = std::env::var("STEPS").ok().map(|v| v.parse().unwrap()).unwrap_or(16);
    for k in 0..steps {
        s.step(dt);
        if std::env::var("VERBOSE").is_ok() {
            eprintln!("  step {k}: krylov {:?} resid {:?}", s.last_cg_iters().unwrap(), s.last_cg_resid().unwrap());
        }
    }
    let box_int = s.mean_velocity()[0] * vol;
    {
        let mut nf = spectral_core::nufft::Nufft::new(&g, 12);
        let uw = nf.interpolate(&s.u[0], &surface.points);
        let umax = uw.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let ucl = nf.interpolate(&s.u[0], &[[pi, pi, pi]])[0];
        eprintln!("  wall slip max|u_x(X_m)| = {umax:.2e}  (centre-line u_x = {ucl:.4e}, ratio {:.1e})", umax / ucl.abs());
        if case == "channel" {
            let h = 0.5 * pi;
            let ys: Vec<[f64; 3]> = (0..=20).map(|i| [1.0, pi - h + h * i as f64 / 10.0, 2.0]).collect();
            let uy = nf.interpolate(&s.u[0], &ys);
            let peak = gx * h * h / (2.0 * nu);
            // independent fluid mean: 48-point Gauss–Legendre of the profile
            let (xs, ws) = gauss_legendre(48);
            let qp: Vec<[f64; 3]> = xs.iter().map(|&x| [1.0, pi + h * x, 2.0]).collect();
            let uq = nf.interpolate(&s.u[0], &qp);
            let mean_q: f64 = uq.iter().zip(&ws).map(|(u, w)| u * w).sum::<f64>() / 2.0;
            eprintln!("  fluid mean by quadrature / exact = {:.9} (err {:+.2e})", mean_q / (gx * h * h / (3.0 * nu)), mean_q / (gx * h * h / (3.0 * nu)) - 1.0);
            for (p, u) in ys.iter().zip(&uy).step_by(2) {
                let yy = p[1] - pi;
                let exact = gx / (2.0 * nu) * (h * h - yy * yy);
                eprintln!("    y-pi={yy:+.3}: u={u:.8e} exact={exact:.8e} err/peak={:+.2e}", (u - exact) / peak);
            }
        }
    }
    let solid_int = surface.solid_integral(&g, &s.u[0]);
    let (iters, resid) = (s.last_cg_iters().unwrap(), s.last_cg_resid().unwrap());
    let report = match case {
        "channel" => {
            let h = 0.5 * pi;
            let fluid_mean = (box_int - solid_int) / (vol - surface.solid_volume());
            let exact = gx * h * h / (3.0 * nu);
            format!("fluid mean / exact = {:.6} (err {:+.4}%)", fluid_mean / exact, 100.0 * (fluid_mean / exact - 1.0))
        }
        _ => {
            let r = sphere_radius();
            let ke = k_exact(r);
            let u_sup = (box_int - solid_int) / vol;
            // Sangani–Acrivos / Zick–Homsy convention: F = τ|∇⟨p⟩| (whole cell, i.e.
            // including the mean-pressure-gradient force G·V_s on the sphere),
            // consistent with Darcy permeability k = τ/(6πaK).
            let k = gx * vol / (6.0 * pi * nu * r * u_sup);
            let f = s.solid_force().unwrap()[0];
            format!("K = {k:.5} (exact {ke:.5}, err {:+.3}%), force check {:.6}", 100.0 * (k / ke - 1.0), f / (gx * vol))
        }
    };
    let tag = surface
        .extension
        .map(|e| format!("ext w={} s={}{}", e.width, e.sigma_rel, surface.basis.as_ref().map(|b| format!(" lsq centres={}", b.len())).unwrap_or_default()))
        .unwrap_or("point".into());
    println!(
        "{case} [{tag}] N={n:3} alpha={alpha} markers={nm}: {report}  CG {iters:?} resid {:?} {:.1}s",
        resid.map(|r| format!("{r:.0e}")),
        t0.elapsed().as_secs_f64()
    );
}

/// Gauss–Legendre nodes and weights on [−1, 1] (Newton on P_n).
fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut xs = vec![0.0; n];
    let mut ws = vec![0.0; n];
    for i in 0..n {
        let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            let dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1e-15 {
                let wt = 2.0 / ((1.0 - x * x) * dp * dp);
                xs[i] = x;
                ws[i] = wt;
                break;
            }
        }
    }
    (xs, ws)
}

/// Sphere radius: π/8 by default, or RADIUS from the environment.
fn sphere_radius() -> f64 {
    std::env::var("RADIUS").ok().map(|v| v.parse().unwrap()).unwrap_or(std::f64::consts::PI / 8.0)
}

/// Sangani & Acrivos (1982) simple-cubic drag K(c), through c^{10/3}.
fn k_exact(r: f64) -> f64 {
    let c = 4.0 / 3.0 * std::f64::consts::PI * r.powi(3) / (2.0 * std::f64::consts::PI).powi(3);
    let c13 = c.cbrt();
    1.0 / (1.0 - 1.7601 * c13 + c - 1.5593 * c * c + 3.9799 * c.powf(8.0 / 3.0) - 3.0734 * c.powf(10.0 / 3.0))
}
