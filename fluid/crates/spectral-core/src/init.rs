//! Initial conditions, including flows with exact solutions used for validation.

use rayon::prelude::*;

use crate::grid::Grid;

/// Fill three physical components from a pointwise function of (x, y, z).
pub fn from_fn(grid: &Grid, f: impl Fn(f64, f64, f64) -> [f64; 3] + Sync) -> [Vec<f64>; 3] {
    let n = grid.n;
    let pts: Vec<[f64; 3]> = (0..grid.real_len())
        .into_par_iter()
        .map(|p| {
            let (i, j, l) = (p / (n * n), (p / n) % n, p % n);
            f(grid.x(i), grid.x(j), grid.x(l))
        })
        .collect();
    std::array::from_fn(|c| pts.iter().map(|v| v[c]).collect())
}

/// 3D Taylor–Green vortex (the Re = 1600 benchmark initial condition, V0 = L = 1).
pub fn taylor_green(grid: &Grid) -> [Vec<f64>; 3] {
    from_fn(grid, |x, y, z| [x.sin() * y.cos() * z.cos(), -x.cos() * y.sin() * z.cos(), 0.0])
}

/// 2D Taylor–Green vortex: exact NS solution u(t) = u(0) e^{-2νt}.
/// Its nonlinear term is a pure gradient, so this tests the projection.
pub fn taylor_green_2d(grid: &Grid) -> [Vec<f64>; 3] {
    from_fn(grid, |x, y, _| [x.sin() * y.cos(), -x.cos() * y.sin(), 0.0])
}

/// Arnold–Beltrami–Childress flow with wavenumber `k`: curl u = k u, so
/// u × ω = 0 and u(t) = u(0) e^{-νk²t} exactly.
pub fn abc(grid: &Grid, a: f64, b: f64, c: f64, k: f64) -> [Vec<f64>; 3] {
    from_fn(grid, move |x, y, z| {
        [
            a * (k * z).sin() + c * (k * y).cos(),
            b * (k * x).sin() + a * (k * z).cos(),
            c * (k * y).sin() + b * (k * x).cos(),
        ]
    })
}
