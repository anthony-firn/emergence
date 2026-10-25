//! Plane Poiseuille flow between two penalized slabs, driven by a uniform
//! pressure gradient G: an unambiguous exact answer for checking
//! convergence of implicit penalization. Fluid occupies |y − π| < h; the
//! exact box-mean velocity is U = G h³ / (3 π ν).
//!
//! usage: channel_probe N c_eta [wall_offset_cells] [fullgrid]
//!   fullgrid: impose the penalization on the N grid (over-determined)
//!   instead of the solution grid M = 2k_max + 1
//!   c_eta > 0: η = c·dx²/ν (grid-matched);  c_eta < 0: η = |c_eta| (absolute)
//!
//! Also reports the error against the exact flow for the staircase geometry
//! the solver actually sees (walls midway between the last fluid and first
//! solid point of the grid the constraint is imposed on), which separates
//! penalization error from geometry sampling.
use spectral_core::{init, Grid, Penalization, Shape, Solver};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let n: usize = a[1].parse().unwrap();
    let c: f64 = a[2].parse().unwrap();
    let wall: f64 = a.get(3).map(|v| v.parse().unwrap()).unwrap_or(0.0);
    let pi = std::f64::consts::PI;
    let (nu, gx, h) = (1.0, 1e-3, 0.5 * pi);
    let g = Grid::new(n);
    let eta = if c > 0.0 { c * g.dx() * g.dx() / nu } else { -c };
    let mut pen = Penalization::sphere_sharp(&g, [pi; 3], 0.1, eta).implicit();
    let hs = h + wall * g.dx();
    let shape = Shape::ChannelY { center: pi, half_width: hs };
    pen.chi = shape.sharp_mask(n);
    pen.shape = Some(shape);
    pen.on_solution_grid = !a.get(4).map(|v| v == "fullgrid").unwrap_or(false);
    let mut s = Solver::new(g.clone(), nu, init::from_fn(&g, |_, _, _| [0.0; 3]));
    s.penalization = Some(pen);
    s.body_force = Some([gx, 0.0, 0.0]);
    for _ in 0..16 {
        s.step(5.0);
    }
    let u = s.mean_velocity()[0];
    let exact = gx * h.powi(3) / (3.0 * pi * nu);
    let full = a.get(4).map(|v| v == "fullgrid").unwrap_or(false);
    let m = if full { n } else { 2 * g.kmax() + 1 };
    let dxm = 2.0 * pi / m as f64;
    let fluid_rows = (0..m).filter(|&j| (j as f64 * dxm - pi).abs() < hs).count();
    let h_eff = 0.5 * fluid_rows as f64 * dxm;
    let exact_eff = gx * h_eff.powi(3) / (3.0 * pi * nu);
    println!(
        "N={n:4} {} (grid {m:3}) eta={eta:.1e}: err vs true walls {:+7.3}%, vs staircase walls {:+8.4}%  iters {:?} resid {:?}",
        if full { "full-grid" } else { "sol-grid " },
        100.0 * (u / exact - 1.0),
        100.0 * (u / exact_eff - 1.0),
        s.last_cg_iters().unwrap(),
        s.last_cg_resid().unwrap().map(|r| format!("{r:.0e}"))
    );
}
