//! GPU convergence study: Stokes flow through a periodic sphere array with
//! point-force surface walls, against Sangani & Acrivos K(c).
//! usage: sphere_array_gpu N [radius] [alpha] [ext_width basis_alpha]
use spectral_core::surface::{resolution, Extension};
use spectral_core::{init, Grid, Surface};
use spectral_cuda::CudaSolver;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let n: usize = a[1].parse().unwrap();
    let pi = std::f64::consts::PI;
    let r: f64 = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(pi / 3.0);
    let alpha: f64 = a.get(3).map(|v| v.parse().unwrap()).unwrap_or(1.0);
    let (nu, vol) = (1.0, (2.0 * pi).powi(3));
    let c = 4.0 / 3.0 * pi * r.powi(3) / vol;
    let ke = 1.0 / (1.0 - 1.7601 * c.cbrt() + c - 1.5593 * c * c + 3.9799 * c.powf(8.0 / 3.0) - 3.0734 * c.powf(10.0 / 3.0));
    let g = Grid::new(n);
    let mut surface = Surface::sphere([pi; 3], r, alpha * resolution(&g));
    if let (Some(wd), Some(ab)) = (a.get(4), a.get(5)) {
        surface.extension = Some(Extension { width: wd.parse().unwrap(), sigma_rel: 1.0 });
        surface = surface.with_basis(ab.parse::<f64>().unwrap() * resolution(&g));
    }
    let gx = 6.0 * pi * nu * r * 1e-3 * ke / vol;
    let mut s = CudaSolver::new(g.clone(), nu, init::from_fn(&g, |_, _, _| [0.0; 3])).unwrap();
    s.set_surface(&surface).unwrap();
    s.body_force = Some([gx, 0.0, 0.0]);
    let t0 = std::time::Instant::now();
    for _ in 0..16 {
        s.try_step(25.0).unwrap();
    }
    let uh = s.try_spectrum(0).unwrap().remove(0);
    let u_sup = (s.try_mean_velocity().unwrap()[0] * vol - surface.solid_integral(&g, &uh)) / vol;
    // Sangani–Acrivos / Zick–Homsy convention: F = τ|∇⟨p⟩| (whole cell, i.e.
    // including the mean-pressure-gradient force G·V_s on the sphere),
    // consistent with Darcy permeability k = τ/(6πaK).
    let k = gx * vol / (6.0 * pi * nu * r * u_sup);
    println!(
        "N={n:4} R={r:.4} {} markers={}: K = {k:.6} (exact {ke:.6}, err {:+.4}%)  {:.1}s",
        if surface.extension.is_some() { "ext" } else { "point" },
        surface.points.len(),
        100.0 * (k / ke - 1.0),
        t0.elapsed().as_secs_f64()
    );
}
