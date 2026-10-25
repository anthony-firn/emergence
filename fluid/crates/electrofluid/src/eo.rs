//! Electro-osmosis: the Helmholtz–Smoluchowski slip as the closure that
//! replaces a resolved Debye layer.
//!
//! Resolved rung: steady Stokes flow with the PNP/PB charge density as a body
//! force,
//!     −μ ∇²u = ρ_e E,   u = 0 on the walls,
//! for a channel of width h with both walls at potential ζ.
//!
//! Fold: when the Debye layers are thin (λ_D/h → 0) the body force is a
//! boundary layer and the interior is a **plug flow** with the
//! Helmholtz–Smoluchowski slip velocity
//!     u_s = −ε ε₀ ζ E / μ
//! as the effective boundary condition (the closure). The validity indicator
//! of the fold is λ_D/h.
//!
//! We verify the fold by resolving the channel for a range of λ_D/h and
//! comparing the centreline velocity u(0) with u_s.

use crate::pnp::{charge_density, debye_length, thomas, EPS0};

/// Helmholtz–Smoluchowski slip velocity (m/s).
pub fn hs_slip(eps_r: f64, zeta: f64, e_field: f64, mu: f64) -> f64 {
    -eps_r * EPS0 * zeta * e_field / mu
}

pub struct Channel {
    pub eps_r: f64,
    pub t: f64,
    pub c_molar: f64,
    pub z: f64,
    pub zeta: f64,
    pub e_field: f64,
    pub mu: f64,
    pub h: f64,
}

impl Channel {
    pub fn lambda(&self) -> f64 {
        debye_length(self.eps_r, self.t, self.c_molar, self.z)
    }

    pub fn slip(&self) -> f64 {
        hs_slip(self.eps_r, self.zeta, self.e_field, self.mu)
    }

    /// Resolve PB potential and charge, then solve μ u'' = −ρ_e E. Returns
    /// (y, u_flow) on the full channel −h/2 … h/2. `n` interior cells.
    pub fn resolve(&self, n: usize) -> (Vec<f64>, Vec<f64>) {
        let _lambda = self.lambda();
        // solve PB on the full channel with both walls at zeta
        let (y, psi) = crate::pnp::pb_solve_dirichlet(
            self.eps_r,
            self.t,
            self.c_molar,
            self.z,
            self.zeta,
            self.zeta,
            self.h,
            n,
        );
        let y: Vec<f64> = y.iter().map(|&v| v - 0.5 * self.h).collect();
        let u: Vec<f64> = psi.iter().map(|&p| crate::pnp::u_of_psi(p, self.t, self.z)).collect();
        // assemble tridiagonal for flow: μ u'' = -ρ_e E, u(±h/2)=0
        let m = n - 1;
        let hh = self.h / n as f64;
        let mut a = vec![0.0; m];
        let mut b = vec![0.0; m];
        let mut c = vec![0.0; m];
        let mut d = vec![0.0; m];
        for i in 1..n {
            let k = i - 1;
            b[k] = -2.0 * self.mu / (hh * hh);
            if k > 0 {
                a[k] = self.mu / (hh * hh);
            }
            if k < m - 1 {
                c[k] = self.mu / (hh * hh);
            }
            // μ v'' = -ρ_e E → v'' = -ρ_e E/μ
            d[k] = -charge_density(u[i], self.c_molar, self.z) * self.e_field;
        }
        let v = thomas(&a, &b, &c, &d);
        let mut flow = vec![0.0];
        flow.extend(v);
        flow.push(0.0);
        (y, flow)
    }

    /// Centreline velocity and the HS prediction.
    pub fn centerline(&self, n: usize) -> (f64, f64) {
        let (_, v) = self.resolve(n);
        (v[n / 2], self.slip())
    }
}
