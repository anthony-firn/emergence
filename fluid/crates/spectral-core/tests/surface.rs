//! Surface-constrained walls on plane Poiseuille flow (exact answer known).

use spectral_core::surface::{resolution, Extension};
use spectral_core::{init, Grid, Solver, Surface};

/// Relative error of the fluid-region mean velocity.
fn channel_error(n: usize, ext: Option<Extension>, alpha: f64) -> f64 {
    let pi = std::f64::consts::PI;
    let (nu, gx, h) = (1.0, 1e-3, 0.5 * pi);
    let vol = (2.0 * pi).powi(3);
    let g = Grid::new(n);
    let mut surface = Surface::channel_y(pi, h, alpha * resolution(&g));
    surface.extension = ext;
    let mut s = Solver::new(g.clone(), nu, init::from_fn(&g, |_, _, _| [0.0; 3]));
    s.surface = Some(surface.clone());
    s.body_force = Some([gx, 0.0, 0.0]);
    for _ in 0..16 {
        s.step(25.0);
    }
    let fluid_mean = (s.mean_velocity()[0] * vol - surface.solid_integral(&g, &s.u[0])) / (vol - surface.solid_volume());
    fluid_mean / (gx * h * h / (3.0 * nu)) - 1.0
}

#[test]
fn point_force_walls_are_first_order() {
    let (e16, e32) = (channel_error(16, None, 1.0), channel_error(32, None, 1.0));
    eprintln!("point forces: N=16 {e16:+.4}, N=32 {e32:+.4}");
    let order = (e16 / e32).abs().log2();
    assert!((0.7..1.6).contains(&order), "order {order}");
}

#[test]
fn smooth_forcing_extension_converges_fast() {
    let ext = Some(Extension { width: 1.0, sigma_rel: 1.0 });
    let (e16, e32) = (channel_error(16, ext, 2.0), channel_error(32, ext, 2.0));
    eprintln!("smooth extension: N=16 {e16:+.2e}, N=32 {e32:+.2e}");
    assert!(e32.abs() < 5e-3, "N=32 error {e32}");
    assert!((e16 / e32).abs() > 8.0, "order too low: {e16} -> {e32}");
}
