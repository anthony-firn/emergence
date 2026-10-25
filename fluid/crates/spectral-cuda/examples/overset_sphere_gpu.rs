//! Overset spectral refinement, two-way coupled by overlapping Schwarz
//! (Dirichlet–Dirichlet), on the GPU: Stokes flow through a periodic array
//! of small spheres (radius a = π/8) against the exact drag K = 1.2139
//! (Sangani & Acrivos, whole-cell convention).
//!
//! Coarse grid (N_c³, the (2π)³ cell), driven by a uniform force G, with a
//! hole r < R_h around the sphere: point-force walls on the hole surface
//! whose prescribed velocity is the patch solution there. Patch (N_p³, a
//! periodic box of side Lp, scale s = Lp/2π; ν/s, τ = t/s, forcing s·G):
//! the shell between the body (no-slip) and an outer sphere R_o > R_h whose
//! velocity is the coarse solution; both walls are smooth forcing
//! extensions. The overlap R_h < r < R_o makes the exchange converge.
//!
//! Drag two ways: the momentum balance K = G·V/(6πμaU), with the
//! superficial velocity U assembled from the coarse field outside the hole
//! and the patch field in a < r < R_h (exact spectral integrals), and the
//! patch body force (×s²).
//!
//! usage: overset_sphere_gpu N_coarse N_patch [t_end]
//! env: RH, RO (hole / outer radii, physical; default 0.55 / 0.75),
//!      HOLE_EXT=1 (smooth-extension hole wall instead of point forces),
//!      DT (default 0.05)
use num_complex::Complex64;
use spectral_core::nufft::Nufft;
use spectral_core::surface::{resolution, Extension};
use spectral_core::{fft::Fft3, init, Grid, Surface};
use spectral_cuda::CudaSolver;

const PI: f64 = std::f64::consts::PI;

fn smooth_step_down(r: f64, r0: f64, r1: f64) -> f64 {
    if r <= r0 {
        return 1.0;
    }
    if r >= r1 {
        return 0.0;
    }
    let t = (r - r0) / (r1 - r0);
    let e = |y: f64| if y > 0.0 { (-1.0 / y).exp() } else { 0.0 };
    e(1.0 - t) / (e(1.0 - t) + e(t))
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let nc: usize = args[1].parse().unwrap();
    let np: usize = args[2].parse().unwrap();
    let t_end: f64 = args.get(3).map(|v| v.parse().unwrap()).unwrap_or(15.0);
    let (nu, a, cen) = (10.0, PI / 8.0, [PI; 3]);
    let vol = (2.0 * PI).powi(3);
    let k_exact = {
        let c = 4.0 / 3.0 * PI * a.powi(3) / vol;
        1.0 / (1.0 - 1.7601 * c.cbrt() + c - 1.5593 * c * c + 3.9799 * c.powf(8.0 / 3.0) - 3.0734 * c.powf(10.0 / 3.0))
    };
    let u0 = 1e-3;
    let gx = 6.0 * PI * nu * a * u0 * k_exact / vol;
    let stokes = |u: f64| 6.0 * PI * nu * a * u;

    // ---------------- coarse grid with a hole r < R_h (geometry, physical units)
    let env = |k: &str, d: f64| std::env::var(k).ok().map(|v| v.parse().unwrap()).unwrap_or(d);
    let (r_hole, r_outer) = (env("RH", 0.55), env("RO", 0.75));
    let gc = Grid::new(nc);
    let hole = if env("HOLE_EXT", 0.0) != 0.0 {
        // smooth extension inside the hole: the overlap must be resolved by
        // the coarse grid, and first-order hole walls distort it
        let mut h = Surface::sphere(cen, r_hole, 0.5 * resolution(&gc));
        h.extension = Some(Extension { width: (6.0 * gc.dx()).min(0.7 * r_hole), sigma_rel: 1.0 });
        h.with_basis(resolution(&gc))
    } else {
        Surface::sphere(cen, r_hole, resolution(&gc))
    };
    let hole_pts = hole.points.clone();
    let mut coarse = CudaSolver::new(gc.clone(), nu, init::from_fn(&gc, |_, _, _| [0.0; 3])).unwrap();
    coarse.set_surface(&hole).unwrap();
    coarse.body_force = Some([gx, 0.0, 0.0]);
    coarse.set_mean_flow([u0, 0.0, 0.0]).unwrap(); // start near the answer

    // ---------------- patch
    let ro = 2.0; // outer sphere in patch units
    let s = r_outer / ro;
    let gp = Grid::new(np);
    let rb = a / s;
    let to_phys = |xi: [f64; 3]| [cen[0] + (xi[0] - PI) * s, cen[1] + (xi[1] - PI) * s, cen[2] + (xi[2] - PI) * s];
    let mut body = Surface::sphere([PI; 3], rb, 0.5 * resolution(&gp));
    body.extension = Some(Extension { width: (6.0 * gp.dx()).min(0.7 * rb), sigma_rel: 1.0 });
    let body = body.with_basis(resolution(&gp));
    let mut outer = Surface::sphere_inverse([PI; 3], ro, 0.5 * resolution(&gp));
    outer.extension = Some(Extension { width: 0.9, sigma_rel: 1.0 });
    let outer = outer.with_basis(resolution(&gp));
    let patch_surface = Surface::union(vec![body, outer]);
    let outer_range = patch_surface.parts[1].0.clone();
    let outer_phys: Vec<[f64; 3]> = patch_surface.points[outer_range.clone()].iter().map(|&x| to_phys(x)).collect();
    let mut patch = CudaSolver::new(gp.clone(), nu / s, init::from_fn(&gp, |_, _, _| [0.0; 3])).unwrap();
    patch.set_surface(&patch_surface).unwrap();
    // windowed forcing s·G
    let fp = init::from_fn(&gp, |x, y, z| {
        let r = ((x - PI).powi(2) + (y - PI).powi(2) + (z - PI).powi(2)).sqrt();
        [s * gx * smooth_step_down(r, ro + 0.3, ro + 0.6), 0.0, 0.0]
    });
    let fft = Fft3::new(&gp);
    let mut fph: [Vec<Complex64>; 3] = std::array::from_fn(|_| vec![Complex64::default(); gp.spec_len()]);
    for c in 0..3 {
        fft.forward(&fp[c], &mut fph[c]);
    }
    patch.set_forcing(&fph).unwrap();

    let mut nf_c = Nufft::new(&gc, 10);
    let mut nf_p = Nufft::new(&gp, 10);
    let hole_xi: Vec<[f64; 3]> = hole_pts.iter().map(|x| [PI + (x[0] - cen[0]) / s, PI + (x[1] - cen[1]) / s, PI + (x[2] - cen[2]) / s]).collect();
    let hole_c = Surface::sphere(cen, r_hole, 1.0); // integrals only
    let ball_ph = Surface::sphere([PI; 3], r_hole / s, 1.0);
    let ball_pa = Surface::sphere([PI; 3], rb, 1.0);
    let dt: f64 = std::env::var("DT").ok().map(|v| v.parse().unwrap()).unwrap_or(0.05);
    let t0 = std::time::Instant::now();
    let mut step = 0usize;
    let report_every = ((1.0 / dt).round() as usize).max(1);
    while coarse.t < t_end - 1e-12 {
        // coarse → patch outer boundary
        let cspec: Vec<Vec<Complex64>> = (0..3).map(|c| coarse.try_spectrum(c).unwrap().remove(0)).collect();
        let vals: Vec<Vec<f64>> = (0..3).map(|c| nf_c.interpolate(&cspec[c], &outer_phys)).collect();
        let mut targets = vec![[0.0; 3]; patch_surface.points.len()];
        for (q, m) in outer_range.clone().enumerate() {
            targets[m] = [vals[0][q], vals[1][q], vals[2][q]];
        }
        patch.set_surface_targets(&targets).unwrap();
        patch.try_step(dt / s).unwrap();
        // patch → coarse hole boundary
        let pspec: Vec<Vec<Complex64>> = (0..3).map(|c| patch.try_spectrum(c).unwrap().remove(0)).collect();
        let pv: Vec<Vec<f64>> = (0..3).map(|c| nf_p.interpolate(&pspec[c], &hole_xi)).collect();
        let ht: Vec<[f64; 3]> = (0..hole_pts.len()).map(|q| [pv[0][q], pv[1][q], pv[2][q]]).collect();
        coarse.set_surface_targets(&ht).unwrap();
        coarse.try_step(dt).unwrap();
        step += 1;
        if step % report_every == 0 || coarse.t >= t_end - 1e-12 {
            let cs0 = coarse.try_spectrum(0).unwrap().remove(0);
            let ps0 = patch.try_spectrum(0).unwrap().remove(0);
            // ∫ true u_x: coarse outside the hole + patch in a < r < R_h (+ 0 in the body)
            let outside = coarse.try_mean_velocity().unwrap()[0] * vol - hole_c.solid_integral(&gc, &cs0);
            let annulus = s.powi(3) * (ball_ph.solid_integral(&gp, &ps0) - ball_pa.solid_integral(&gp, &ps0));
            let u_sup = (outside + annulus) / vol;
            let k_mom = gx * vol / stokes(u_sup);
            let fb = patch.surface_part_forces().unwrap()[0][0] * s * s;
            let k_patch = fb / stokes(u_sup);
            println!(
                "t={:6.2}: K(momentum) = {k_mom:.6} ({:+.3}%)  K(patch body force) = {k_patch:.6} ({:+.3}%)  solves coarse {:?} patch {:?}  {:.0}s",
                coarse.t,
                100.0 * (k_mom / k_exact - 1.0),
                100.0 * (k_patch / k_exact - 1.0),
                coarse.last_cg().map(|c| c.0),
                patch.last_cg().map(|c| c.0),
                t0.elapsed().as_secs_f64()
            );
        }
    }
    println!(
        "exact K = {k_exact:.6}; coarse {nc}^3 (dx {:.3}), hole R_h = {r_hole} ({}), patch {np}^3 (physical dx {:.4}, {:.1}x finer), outer R_o = {r_outer}",
        gc.dx(),
        if hole.extension.is_some() { "extension" } else { "point forces" },
        gp.dx() * s,
        gc.dx() / (gp.dx() * s)
    );
}
