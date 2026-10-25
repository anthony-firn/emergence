//! Brinkman volume penalization: solid objects and an outflow sponge.
//!
//! The momentum equation gains the term
//!
//!   -(χ/η)(u - u_s) - (σ/η_s)(u - U∞)
//!
//! where χ ∈ [0,1] marks the solid (moving with velocity u_s), σ marks a
//! sponge layer relaxing the flow to the free stream U∞ so the wake does not
//! re-enter through the periodic boundary. As η → 0 the penalized solution
//! converges to Navier–Stokes with no-slip walls (Angot, Bruneau & Fabrie,
//! Numer. Math. 81, 1999); the wall model error scales like √(νη), which is
//! reported as a separate error term, not hidden in the discretization error.
//!
//! The term is added in physical space inside the right-hand side, where the
//! velocity is already available, so objects cost no extra FFTs. Explicit
//! RK4 is stable for this term when dt ≤ 2.78 η (see [`Penalization::max_dt`]).

use rayon::prelude::*;

use crate::grid::Grid;

/// Analytic body geometry, so a sharp mask can be sampled on any grid (the
/// implicit solve uses the odd-sized solution grid, see `imex.rs`).
#[derive(Clone, Copy, Debug)]
pub enum Shape {
    Sphere { center: [f64; 3], radius: f64 },
    /// Solid where |y − center| ≥ half_width (a plane channel along x, z).
    ChannelY { center: f64, half_width: f64 },
    /// Fluid *inside* the sphere, solid outside (the fluid region of an
    /// overset refinement patch; the outside is a fictitious buffer).
    SphereInverse { center: [f64; 3], radius: f64 },
    /// Solid cylinder along z: solid where the (x, y) distance from
    /// `center` is < radius (an electrode rod, a stirrer shaft).
    CylinderZ { center: [f64; 2], radius: f64 },
    /// Fluid inside the z-cylinder, solid outside (a pipe or a dish wall).
    CylinderZInverse { center: [f64; 2], radius: f64 },
}

impl Shape {
    /// Periodic (x, y) distance from a z-axis through `c`.
    pub fn axis_distance(x: [f64; 3], c: [f64; 2]) -> f64 {
        let tp = 2.0 * std::f64::consts::PI;
        let d = |a: f64, b: f64| {
            let mut v = (a - b).rem_euclid(tp);
            if v > 0.5 * tp {
                v -= tp;
            }
            v
        };
        (d(x[0], c[0]).powi(2) + d(x[1], c[1]).powi(2)).sqrt()
    }
}

impl Shape {
    /// Sharp (0/1) mask sampled at the points of an m³ grid on [0, 2π)³.
    pub fn sharp_mask(&self, m: usize) -> Vec<f64> {
        let dx = 2.0 * std::f64::consts::PI / m as f64;
        (0..m * m * m)
            .into_par_iter()
            .map(|p| {
                let (i, j, l) = (p / (m * m), (p / m) % m, p % m);
                let x = [i as f64 * dx, j as f64 * dx, l as f64 * dx];
                let solid = match *self {
                    Shape::Sphere { center: c, radius } => {
                        let d2: f64 = (0..3).map(|q| (x[q] - c[q]).powi(2)).sum();
                        d2.sqrt() < radius
                    }
                    Shape::ChannelY { center, half_width } => (x[1] - center).abs() >= half_width,
                    Shape::SphereInverse { center: c, radius } => {
                        let d2: f64 = (0..3).map(|q| (x[q] - c[q]).powi(2)).sum();
                        d2.sqrt() > radius
                    }
                    Shape::CylinderZ { center, radius } => Shape::axis_distance(x, center) < radius,
                    Shape::CylinderZInverse { center, radius } => Shape::axis_distance(x, center) > radius,
                };
                if solid { 1.0 } else { 0.0 }
            })
            .collect()
    }
}

#[derive(Clone, Debug)]
pub struct Penalization {
    /// Solid indicator χ (smoothed), N³ physical values.
    pub chi: Vec<f64>,
    pub eta: f64,
    pub u_solid: [f64; 3],
    /// Sponge indicator σ, N³ physical values (all zero = no sponge).
    pub sponge: Vec<f64>,
    pub eta_sponge: f64,
    pub u_inf: [f64; 3],
    /// Treat the solid term implicitly (IMEX step with a CG solve for the
    /// penalty force) instead of explicitly in the right-hand side. Lets η be
    /// tiny at the CFL time step.
    pub implicit: bool,
    /// Analytic geometry when known (sharp masks). With it, the implicit
    /// solve imposes the penalization on the solution grid, where the number
    /// of constraints matches the number of retained Fourier modes.
    pub shape: Option<Shape>,
    /// Impose the implicit penalization on the odd solution grid
    /// M = 2k_max + 1 (one constraint per retained mode) instead of the N
    /// grid. Off by default: measured on plane Poiseuille against the exact
    /// staircase-geometry flow it converges from above (+23.5/+10.5/+7.2% at
    /// N = 32/64/96), while the N grid converges from below and is more
    /// accurate at N = 64 (−4.6%). Kept for experiments.
    pub on_solution_grid: bool,
    /// Spatially varying target velocity for the χ term (physical, N³ per
    /// component), replacing the uniform `u_solid` — e.g. a coarse grid
    /// relaxed toward a refinement patch's solution (explicit mode only).
    pub u_target: Option<[Vec<f64>; 3]>,
}

/// Smoothed Heaviside with compact transition: exactly 1 for s ≤ -δ, exactly
/// 0 for s ≥ δ, C^∞ in between. Compact support matters: a tanh profile has
/// exponential tails that penalize the fluid outside the body with strength
/// χ/η, so the effective body grows as η → 0 (observed: Stokes drag rising
/// 1.55 → 2.07 × theory as η went 1e-3 → 1e-4).
#[inline]
fn smooth_inside(s: f64, delta: f64) -> f64 {
    if s <= -delta {
        return 1.0;
    }
    if s >= delta {
        return 0.0;
    }
    // x ∈ (0,1) from inside to outside; standard C^∞ transition function
    let x = (s + delta) / (2.0 * delta);
    let f = |y: f64| if y > 0.0 { (-1.0 / y).exp() } else { 0.0 };
    let (a, b) = (f(1.0 - x), f(x));
    a / (a + b)
}

impl Penalization {
    /// Sphere of diameter `d` centered at `c`, with mask smoothing width
    /// `delta` (≈ 1–1.5 grid spacings limits Gibbs oscillations).
    ///
    /// The smoothing layer lies *inside* the body, χ falling from 1 at
    /// r = R - 2δ to exactly 0 at r = R: as η → 0 every point with χ > 0
    /// becomes rigid, so this makes the η → 0 limit the true geometry rather
    /// than a sphere inflated by δ.
    pub fn sphere(grid: &Grid, c: [f64; 3], d: f64, delta: f64, eta: f64) -> Self {
        let n = grid.n;
        let r = 0.5 * d;
        let chi = (0..grid.real_len())
            .into_par_iter()
            .map(|p| {
                let (i, j, l) = (p / (n * n), (p / n) % n, p % n);
                let dx = [grid.x(i) - c[0], grid.x(j) - c[1], grid.x(l) - c[2]];
                let dist = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
                smooth_inside(dist - r + delta, delta)
            })
            .collect();
        Self {
            chi,
            eta,
            u_solid: [0.0; 3],
            sponge: vec![0.0; grid.real_len()],
            eta_sponge: f64::INFINITY,
            u_inf: [0.0; 3],
            implicit: false,
            shape: None,
            on_solution_grid: false,
            u_target: None,
        }
    }

    /// Sphere with a sharp (binary) mask, χ = 1 at grid points with r < R.
    ///
    /// For implicit penalization with small η this is the right mask: every
    /// point with χ > 0 is rigid once η/(hγχ) ≪ 1, so any compactly
    /// supported mask converges to this staircase body anyway, but a smooth
    /// one only does so like 1/ln(1/η) (the rigid front creeps through the
    /// region where χ is tiny). Measured: smooth mask drag changed 9% between
    /// η = 1e-6 and 1e-8; the sharp mask is η-independent there.
    pub fn sphere_sharp(grid: &Grid, c: [f64; 3], d: f64, eta: f64) -> Self {
        Self::sphere_sharp_offset(grid, c, d, eta, 0.0)
    }

    /// Sharp sphere whose solid points are those with r < R − `offset`.
    /// With a staircase body the effective no-slip wall sits roughly midway
    /// between the last solid and first fluid point, so offset = dx/2 places
    /// that effective wall at r ≈ R.
    pub fn sphere_sharp_offset(grid: &Grid, c: [f64; 3], d: f64, eta: f64, offset: f64) -> Self {
        let n = grid.n;
        let r = 0.5 * d - offset;
        let mut p = Self::sphere(grid, c, d, grid.dx(), eta);
        p.shape = Some(Shape::Sphere { center: c, radius: r });
        p.chi.par_iter_mut().enumerate().for_each(|(q, x)| {
            let (i, j, l) = (q / (n * n), (q / n) % n, q % n);
            let dx = [grid.x(i) - c[0], grid.x(j) - c[1], grid.x(l) - c[2]];
            *x = if (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt() < r { 1.0 } else { 0.0 };
        });
        p
    }

    /// Add a sponge slab across x ∈ [x0, x0 + width) (periodic), relaxing the
    /// flow toward `u_inf` with time scale `eta_sponge`.
    pub fn with_sponge(mut self, grid: &Grid, x0: f64, width: f64, eta_sponge: f64, u_inf: [f64; 3]) -> Self {
        let n = grid.n;
        let two_pi = 2.0 * std::f64::consts::PI;
        let delta = 2.0 * grid.dx();
        self.sponge.par_iter_mut().enumerate().for_each(|(p, s)| {
            let i = p / (n * n);
            // distance from slab center, periodic
            let xc = x0 + 0.5 * width;
            let mut d = (grid.x(i) - xc).rem_euclid(two_pi);
            if d > 0.5 * two_pi {
                d -= two_pi;
            }
            *s = smooth_inside(d.abs() - 0.5 * width, delta);
        });
        self.eta_sponge = eta_sponge;
        self.u_inf = u_inf;
        self
    }

    /// Use the implicit (IMEX) treatment of the solid term.
    pub fn implicit(mut self) -> Self {
        self.implicit = true;
        self
    }

    /// Largest stable explicit RK4 step for the explicitly treated terms.
    pub fn max_dt(&self) -> f64 {
        if self.implicit {
            2.5 * self.eta_sponge
        } else {
            2.5 * self.eta.min(self.eta_sponge)
        }
    }

    /// Subtract the explicitly treated penalization terms from `n` at every
    /// physical point (the solid term is skipped when it is implicit).
    pub fn apply(&self, u: &[Vec<f64>; 3], n: &mut [Vec<f64>; 3]) {
        let ie = if self.implicit { 0.0 } else { 1.0 / self.eta };
        let is = 1.0 / self.eta_sponge;
        let [n0, n1, n2] = n;
        n0.par_iter_mut()
            .zip(n1.par_iter_mut())
            .zip(n2.par_iter_mut())
            .enumerate()
            .for_each(|(p, ((a, b), c))| {
                let cs = self.chi[p] * ie;
                let cf = self.sponge[p] * is;
                let t = match &self.u_target {
                    Some(t) => [t[0][p], t[1][p], t[2][p]],
                    None => self.u_solid,
                };
                *a -= cs * (u[0][p] - t[0]) + cf * (u[0][p] - self.u_inf[0]);
                *b -= cs * (u[1][p] - t[1]) + cf * (u[1][p] - self.u_inf[1]);
                *c -= cs * (u[2][p] - t[2]) + cf * (u[2][p] - self.u_inf[2]);
            });
    }

    /// Force exerted by the fluid on the solid, F = (1/η) ∫ χ (u - u_s) dx,
    /// over the (2π)³ box (unit density).
    pub fn solid_force(&self, grid: &Grid, u: &[Vec<f64>; 3]) -> [f64; 3] {
        let vol_per_point = grid.dx().powi(3);
        let s = (0..grid.real_len())
            .into_par_iter()
            .map(|p| {
                let w = self.chi[p];
                [w * (u[0][p] - self.u_solid[0]), w * (u[1][p] - self.u_solid[1]), w * (u[2][p] - self.u_solid[2])]
            })
            .reduce(|| [0.0; 3], |a, b| [a[0] + b[0], a[1] + b[1], a[2] + b[2]]);
        [s[0], s[1], s[2]].map(|v| v * vol_per_point / self.eta)
    }
}
