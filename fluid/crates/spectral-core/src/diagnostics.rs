//! Global diagnostics. All averages are over the box, (2π)⁻³ ∫ · dx, which is
//! the normalization used by the Taylor–Green Re = 1600 reference data.

use rayon::prelude::*;

use crate::solver::Solver;

#[derive(Clone, Copy, Debug, Default)]
pub struct Diagnostics {
    pub t: f64,
    /// ⟨|u|²⟩/2
    pub energy: f64,
    /// ⟨|ω|²⟩/2
    pub enstrophy: f64,
    /// 2ν·enstrophy, equal to -dE/dt for unforced flow
    pub dissipation: f64,
    /// max |k·û| / |k| over modes, a divergence check (should be ~1e-16·|u|)
    pub max_div: f64,
    /// Fraction of energy in the outer third of retained shells (resolution check)
    pub tail_energy_fraction: f64,
}

impl Solver {
    pub fn diagnostics(&self) -> Diagnostics {
        let g = &self.grid;
        let k2 = self.k2();
        let kmax = (g.n / 3) as f64;
        let [a, b, c] = &self.u;
        let (e, z, div, tail) = (0..g.spec_len())
            .into_par_iter()
            .map(|idx| {
                let w = g.weight(idx);
                let m2 = a[idx].norm_sqr() + b[idx].norm_sqr() + c[idx].norm_sqr();
                let k = g.k(idx);
                let kn = k2[idx].sqrt();
                let div = if kn > 0.0 {
                    (a[idx] * k[0] + b[idx] * k[1] + c[idx] * k[2]).norm() / kn
                } else {
                    0.0
                };
                let tail = if kn > 2.0 * kmax / 3.0 { 0.5 * w * m2 } else { 0.0 };
                (0.5 * w * m2, 0.5 * w * k2[idx] * m2, div, tail)
            })
            .reduce(
                || (0.0, 0.0, 0.0, 0.0),
                |x, y| (x.0 + y.0, x.1 + y.1, x.2.max(y.2), x.3 + y.3),
            );
        Diagnostics {
            t: self.t,
            energy: e,
            enstrophy: z,
            dissipation: 2.0 * self.nu * z,
            max_div: div,
            tail_energy_fraction: if e > 0.0 { tail / e } else { 0.0 },
        }
    }

    /// max over the grid of |u|, used for the CFL time step.
    pub fn max_speed(&mut self) -> f64 {
        let comps: Vec<Vec<f64>> = (0..3).map(|c| self.velocity_phys(c)).collect();
        (0..comps[0].len())
            .into_par_iter()
            .map(|p| (comps[0][p].powi(2) + comps[1][p].powi(2) + comps[2][p].powi(2)).sqrt())
            .reduce(|| 0.0, f64::max)
    }

    /// Advective CFL time step dt = cfl·dx/max|u|.
    pub fn cfl_dt(&mut self, cfl: f64) -> f64 {
        cfl * self.grid.dx() / self.max_speed().max(1e-12)
    }
}
