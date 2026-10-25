//! Common interface for the CPU reference solver and accelerated backends,
//! so drivers and cross-checks are written once.

use crate::diagnostics::Diagnostics;
use crate::penalization::Penalization;
use crate::surface::Surface;

pub trait Backend {
    fn name(&self) -> &'static str;
    fn n(&self) -> usize;
    fn nu(&self) -> f64;
    fn time(&self) -> f64;
    /// Advance one IF-RK4 step of size `h`.
    fn step(&mut self, h: f64);
    fn diagnostics(&mut self) -> Diagnostics;
    /// Advective CFL step dt = cfl·dx/max|u|.
    fn cfl_dt(&mut self, cfl: f64) -> f64;
    /// Physical velocity component `c` as f64 (host copy).
    fn velocity_phys(&mut self, c: usize) -> Vec<f64>;
    /// Install solid objects / sponge (Brinkman penalization).
    fn set_penalization(&mut self, p: Penalization);
    /// Fix the uniform (k = 0) free-stream velocity.
    fn set_mean_flow(&mut self, m: [f64; 3]);
    /// Force of the fluid on the penalized solid (unit density), if any.
    fn solid_force(&mut self) -> Option<[f64; 3]>;
    /// Install surface-constrained walls (markers on the true geometry).
    fn set_surface(&mut self, s: Surface);
}

impl Backend for crate::Solver {
    fn name(&self) -> &'static str {
        "cpu-f64"
    }
    fn n(&self) -> usize {
        self.grid.n
    }
    fn nu(&self) -> f64 {
        self.nu
    }
    fn time(&self) -> f64 {
        self.t
    }
    fn step(&mut self, h: f64) {
        crate::Solver::step(self, h)
    }
    fn diagnostics(&mut self) -> Diagnostics {
        crate::Solver::diagnostics(self)
    }
    fn cfl_dt(&mut self, cfl: f64) -> f64 {
        crate::Solver::cfl_dt(self, cfl)
    }
    fn velocity_phys(&mut self, c: usize) -> Vec<f64> {
        crate::Solver::velocity_phys(self, c)
    }
    fn set_penalization(&mut self, p: Penalization) {
        self.penalization = Some(p);
    }
    fn set_mean_flow(&mut self, m: [f64; 3]) {
        crate::Solver::set_mean_flow(self, m)
    }
    fn solid_force(&mut self) -> Option<[f64; 3]> {
        crate::Solver::solid_force(self)
    }
    fn set_surface(&mut self, s: Surface) {
        self.surface = Some(s);
    }
}
