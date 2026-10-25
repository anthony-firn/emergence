//! Diagnose implicit penalization on the Stokes sphere-array case, against
//! Hasimoto / Sangani–Acrivos: drag = 6πμaU·K, K(c) = 1.2139 here.
//!
//! usage: stokes_probe N dt steps eta [center_off_cells] [wall_off_cells] [mode]
//!   mode "fixed"  (default): hold the mean flow U = 1e-3, drag from ∫g
//!   mode "driven": prescribe the pressure gradient G, measure U; drag = G·V
use spectral_core::{init, Grid, Penalization, Solver};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let n: usize = a[1].parse().unwrap();
    let dt: f64 = a[2].parse().unwrap();
    let steps: usize = a[3].parse().unwrap();
    let eta: f64 = a[4].parse().unwrap();
    let off: f64 = a.get(5).map(|v| v.parse().unwrap()).unwrap_or(0.0);
    let wall: f64 = a.get(6).map(|v| v.parse().unwrap()).unwrap_or(0.0);
    let driven = a.get(7).map(|v| v == "driven").unwrap_or(false);
    let two_pi = 2.0 * std::f64::consts::PI;
    let (d, nu) = (two_pi / 8.0, 1.0);
    let radius = 0.5 * d;
    let vol = two_pi.powi(3);
    let k_exact = 1.2139;
    let g = Grid::new(n);
    let pen = Penalization::sphere_sharp_offset(&g, [0.5 * two_pi + off * g.dx(); 3], d, eta, wall * g.dx()).implicit();
    let mut s = Solver::new(g.clone(), nu, init::from_fn(&g, |_, _, _| [0.0; 3]));
    s.penalization = Some(pen);
    // G chosen so that U ≈ 1e-3 at the exact K
    let gx = 6.0 * std::f64::consts::PI * nu * radius * 1e-3 * k_exact / vol;
    if driven {
        s.body_force = Some([gx, 0.0, 0.0]);
    } else {
        s.set_mean_flow([1e-3, 0.0, 0.0]);
    }
    for k in 0..steps {
        s.step(dt);
        if k + 1 < steps {
            continue;
        }
        let f_g = s.solid_force().unwrap()[0];
        let u = s.mean_velocity()[0];
        let stokes = 6.0 * std::f64::consts::PI * nu * radius * u;
        let (k_mom, k_g) = (gx * vol / stokes, f_g / stokes);
        println!(
            "N={n} eta={eta:e} {} off={off} wall={wall}: U={u:.6e} K(momentum)={} K(int g)={k_g:.4} (exact {k_exact}) iters {:?} resid {:?}",
            if driven { "driven" } else { "fixed " },
            if driven { format!("{k_mom:.4}") } else { "  -   ".into() },
            s.last_cg_iters().unwrap(),
            s.last_cg_resid().unwrap().map(|r| format!("{r:.0e}"))
        );
    }
}
