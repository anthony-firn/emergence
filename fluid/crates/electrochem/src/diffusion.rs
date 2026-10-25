//! 1D diffusion and the Cottrell transient.
//!
//! The rung below the interface: after a potential step to a transport-limited
//! value the surface concentration is clamped and the current is a pure
//! diffusion problem,
//!
//! ```text
//!   ∂c/∂t = D ∂²c/∂x²,   c(0,t)=0,   c(∞,t)=c_bulk
//!   i(t)  = nFA D ∂c/∂x|₀ = nFAc_bulk √(D/(πt))         (Cottrell)
//! ```
//!
//! [`Diffusion1D`] integrates it with Crank–Nicolson (2nd order in space and
//! time, unconditionally stable) and a Thomas tridiagonal solve. The exact
//! benchmark is [`cottrell_current`] (Bard & Faulkner eq. 5.2.11), and the
//! order of the scheme is checked on a smooth Neumann eigenmode, for which the
//! exact solution is `exp(−Dπ²t/L²) cos(πx/L)`.

use crate::F_FARADAY;

/// Boundary condition on one end of the 1D domain.
#[derive(Clone, Copy, Debug)]
pub enum Boundary {
    /// Fixed concentration `c`.
    Dirichlet(f64),
    /// Zero flux `∂c/∂x = 0`.
    NoFlux,
}

/// Uniform-grid Crank–Nicolson diffusion solver on `[0, L]`.
pub struct Diffusion1D {
    pub n: usize,
    pub dx: f64,
    pub d: f64,
    left: Boundary,
    right: Boundary,
    c: Vec<f64>,
}

impl Diffusion1D {
    /// Build a solver with `n` cells (`n+1` nodes) on `[0, length]`.
    pub fn new(n: usize, length: f64, d: f64, left: Boundary, right: Boundary, init: Vec<f64>) -> Self {
        assert_eq!(init.len(), n + 1);
        Self { n, dx: length / n as f64, d, left, right, c: init }
    }

    /// Uniform initial condition.
    pub fn uniform(n: usize, length: f64, d: f64, left: Boundary, right: Boundary, c0: f64) -> Self {
        Self::new(n, length, d, left, right, vec![c0; n + 1])
    }

    pub fn concentrations(&self) -> &[f64] {
        &self.c
    }

    /// Replace the left boundary condition (e.g. for a potential ramp).
    pub fn set_left(&mut self, b: Boundary) {
        self.left = b;
    }

    /// Replace the right boundary condition.
    pub fn set_right(&mut self, b: Boundary) {
        self.right = b;
    }

    /// Advance one Crank–Nicolson step of size `dt`.
    pub fn step(&mut self, dt: f64) {
        let r = self.d * dt / (self.dx * self.dx);
        let m = self.n + 1;
        // Tridiagonal system a_i x_{i-1} + b_i x_i + c_i x_{i+1} = rhs_i.
        let mut a = vec![0.0; m];
        let mut b = vec![0.0; m];
        let mut cc = vec![0.0; m];
        let mut rhs = vec![0.0; m];

        // Left boundary row.
        match self.left {
            Boundary::Dirichlet(v) => {
                b[0] = 1.0;
                rhs[0] = v;
            }
            Boundary::NoFlux => {
                // Ghost c_{-1} = c_1: (−r) c_{0} ... interior-like row.
                b[0] = 1.0 + r;
                cc[0] = -r;
                rhs[0] = (1.0 - r) * self.c[0] + r * self.c[1];
            }
        }

        for i in 1..self.n {
            a[i] = -0.5 * r;
            b[i] = 1.0 + r;
            cc[i] = -0.5 * r;
            rhs[i] = 0.5 * r * self.c[i - 1] + (1.0 - r) * self.c[i] + 0.5 * r * self.c[i + 1];
        }

        // Right boundary row (index n).
        match self.right {
            Boundary::Dirichlet(v) => {
                b[self.n] = 1.0;
                rhs[self.n] = v;
            }
            Boundary::NoFlux => {
                a[self.n] = -r;
                b[self.n] = 1.0 + r;
                rhs[self.n] = r * self.c[self.n - 1] + (1.0 - r) * self.c[self.n];
            }
        }

        self.c = thomas(&mut a, &mut b, &mut cc, &mut rhs);
    }

    /// Advance `steps` steps of size `dt`.
    pub fn run(&mut self, dt: f64, steps: usize) {
        for _ in 0..steps {
            self.step(dt);
        }
    }

    /// Diffusion current `i = n F A D (c_1 − c_0)/dx` at the left boundary
    /// (positive when material flows from the bulk toward `x = 0`).
    pub fn left_current(&self, n_e: f64, area: f64) -> f64 {
        n_e * F_FARADAY * area * self.d * (self.c[1] - self.c[0]) / self.dx
    }
}

/// Thomas algorithm for a tridiagonal system; `a` and `cc` are overwritten.
fn thomas(a: &mut [f64], b: &mut [f64], cc: &mut [f64], d: &mut [f64]) -> Vec<f64> {
    let n = b.len();
    for i in 1..n {
        let m = a[i] / b[i - 1];
        b[i] -= m * cc[i - 1];
        d[i] -= m * d[i - 1];
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d[n - 1] / b[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (d[i] - cc[i] * x[i + 1]) / b[i];
    }
    x
}

/// Exact Cottrell current (A) after a potential step to the transport-limited
/// plateau: `i(t) = n F A c_bulk √(D/(πt))`.
///
/// Bard & Faulkner, *Electrochemical Methods*, 2nd ed., eq. 5.2.11.
pub fn cottrell_current(t: f64, n_e: f64, area: f64, c_bulk: f64, d: f64) -> f64 {
    n_e * F_FARADAY * area * c_bulk * (d / (std::f64::consts::PI * t)).sqrt()
}

/// Solve the Cottrell problem in the given *dimensionless* units and return
/// `(t, i_numeric, i_exact)` at each requested time. `D`, `L`, `c_bulk` and the
/// times are in whatever consistent units the caller uses; `n_e * area * F`
/// only scales the current, so the comparison is reported normalised.
pub fn cottrell_study(
    n: usize,
    length: f64,
    d: f64,
    dt: f64,
    sample_times: &[f64],
) -> Vec<(f64, f64, f64)> {
    let mut solver = Diffusion1D::uniform(n, length, d, Boundary::Dirichlet(0.0), Boundary::NoFlux, 1.0);
    let mut out = Vec::new();
    let mut t = 0.0;
    let mut k = 0usize;
    let mut next = 0usize;
    // Number of steps to just reach the last sample time.
    let t_end = *sample_times.last().unwrap();
    let steps = (t_end / dt).ceil() as usize;
    while next < sample_times.len() && k <= steps {
        while next < sample_times.len() && t + 1e-12 >= sample_times[next] {
            // Current uses (c_1 - c_0)/dx; scale to the analytic prefactor.
            let num = (solver.concentrations()[1] - solver.concentrations()[0]) / solver.dx;
            // Exact nondimensional flux for c_bulk = 1 is √(D/(πt)).
            let exact = (d / (std::f64::consts::PI * sample_times[next])).sqrt();
            out.push((sample_times[next], num, exact));
            next += 1;
        }
        solver.step(dt);
        t += dt;
        k += 1;
    }
    out
}

/// L2 error of the Crank–Nicolson solve of `∂c/∂t = D ∂²c/∂x²` against the
/// exact Neumann eigenmode `cos(πx/L) exp(−Dπ²t/L²)` at time `t_end`.
///
/// Used for the order study: the scheme is 2nd order in `dt` and `dx`.
pub fn eigenmode_error(n: usize, length: f64, d: f64, dt: f64, t_end: f64) -> f64 {
    let init: Vec<f64> = (0..=n)
        .map(|i| {
            let x = i as f64 * length / n as f64;
            (std::f64::consts::PI * x / length).cos()
        })
        .collect();
    let mut solver = Diffusion1D::new(n, length, d, Boundary::NoFlux, Boundary::NoFlux, init);
    let steps = (t_end / dt).round() as usize;
    solver.run(dt, steps);
    let decay = (-d * std::f64::consts::PI.powi(2) * (steps as f64) * dt / (length * length)).exp();
    let mut err = 0.0;
    for i in 0..=n {
        let x = i as f64 * length / n as f64;
        let exact = decay * (std::f64::consts::PI * x / length).cos();
        let e = solver.concentrations()[i] - exact;
        err += e * e;
    }
    (err / (n as f64 + 1.0)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cottrell_matches_within_half_a_percent() {
        // Dimensionless: D = 1, L = 6 so the finite length is invisible while
        // √(Dt) ≲ 1. Compare the numerical surface flux to √(D/πt).
        let sample = [0.5, 1.0];
        let res = cottrell_study(800, 6.0, 1.0, 2e-4, &sample);
        for (t, num, exact) in res {
            let rel = (num - exact).abs() / exact;
            assert!(rel < 5e-3, "Cottrell rel error at t={t}: {rel}");
        }
    }

    #[test]
    fn eigenmode_second_order_in_time() {
        // Fix dx, refine dt; the spatial error is a constant offset, so use
        // the differences e(dt) − e(dt/2) to isolate the O(dt²) term.
        let (n, l, d, t_end) = (1024, 1.0, 1.0, 0.1);
        let e1 = eigenmode_error(n, l, d, 1e-2, t_end);
        let e2 = eigenmode_error(n, l, d, 5e-3, t_end);
        let e3 = eigenmode_error(n, l, d, 2.5e-3, t_end);
        let ratio = (e1 - e2) / (e2 - e3);
        assert!(ratio > 3.4 && ratio < 4.6, "time order ratio {ratio} (e={e1},{e2},{e3})");
    }

    #[test]
    fn eigenmode_second_order_in_space() {
        // Fix dt, refine dx; the temporal error is a constant offset.
        let (l, d, dt, t_end) = (1.0, 1.0, 2e-5, 0.05);
        let e1 = eigenmode_error(16, l, d, dt, t_end);
        let e2 = eigenmode_error(32, l, d, dt, t_end);
        let e3 = eigenmode_error(64, l, d, dt, t_end);
        let ratio = (e1 - e2) / (e2 - e3);
        assert!(ratio > 3.0 && ratio < 5.0, "space order ratio {ratio} (e={e1},{e2},{e3})");
    }

    #[test]
    fn mass_is_conserved_with_no_flux_walls() {
        let mut s = Diffusion1D::uniform(50, 1.0, 1.0, Boundary::NoFlux, Boundary::NoFlux, 1.0);
        s.run(1e-3, 100);
        let sum: f64 = s.concentrations().iter().sum();
        assert!((sum - 51.0).abs() < 1e-9);
    }
}
