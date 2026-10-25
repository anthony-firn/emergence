//! Electrorheology (ER): dielectric particles in an insulating oil polarise in
//! an applied field, attract along the field and chain.
//!
//! Point-dipole model. A sphere of radius a in a fluid of relative permittivity
//! ε_f has induced dipole moment p = 4π ε₀ ε_f a³ β E along the field, with the
//! Clausius–Mossotti factor β = (ε_p − ε_f)/(ε_p + 2 ε_f). The force between
//! two aligned dipoles separated by r is
//!     F = (3 p² / 4π ε₀ ε_f r⁴)[2 cosθ ŷ + (1 − 5 cos²θ) r̂]
//! whose contact value (r = 2a, θ = 0) is the familiar
//!     F = (3/16) π ε₀ ε_f d² β² E² · (something of order one)
//! with d = 2a; we use the exact tensor form above. Scaling with φ gives the
//! electrostatic yield stress
//!     τ_y = C φ ε₀ ε_f β² E²,
//! with C a structure factor of order 0.1–1 (Klingenberg 1993; Parthasarathy &
//! Klingenberg 1996; Halsey 1992).
//!
//! Fold: the suspension is a Bingham fluid, τ = τ_y + η γ̇ (τ > τ_y), with the
//! Mason number Mn = η γ̇ /(ε₀ ε_f β² E²) as the validity indicator: chains
//! survive for Mn ≲ 1 and break above.
//!
//! This module runs a small 2D overdamped particle simulation (dipole pair
//! forces + short-range repulsion + Stokes drag, Lees–Edwards shear) to show
//! chain formation and measure the structural shear stress. Simulation units:
//! length a, force F₀ = ε₀ ε_f β² E² a², mobility 1 (so 6πμa = 1). The
//! dimensionless shear coefficient is `shear`; the physical Mason number is
//! Mn = shear/(6π).

use std::f64::consts::PI;

/// Full dipole pair force on i from j, in units of F₀ = ε₀ ε_f β² E² a²,
/// with the separation (dx, dy) = r_j − r_i in units of a and dipoles along
/// +y. Returns the force on i (attractive when aligned along the field).
pub fn dipole_force(dx: f64, dy: f64) -> (f64, f64) {
    let r2 = dx * dx + dy * dy;
    if r2 < 1e-12 {
        return (0.0, 0.0);
    }
    let r = r2.sqrt();
    let rhatx = dx / r;
    let rhaty = dy / r;
    let cos = rhaty; // field is +y
    // The standard tensor formula gives the force on j when r = r_j − r_i;
    // negate it to get the force on i. F/F0 = −12π/r⁴[(1−5cos²)r̂ + 2cos ŷ].
    let c = -12.0 * PI / (r2 * r2);
    let radial = 1.0 - 5.0 * cos * cos;
    let fx = c * (radial * rhatx);
    let fy = c * (radial * rhaty + 2.0 * cos);
    (fx, fy)
}

#[derive(Clone, Debug)]
pub struct ErSim {
    pub n: usize,
    pub lx: f64,
    pub ly: f64,
    /// dimensionless shear coefficient γ̇ in simulation units
    pub shear: f64,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub fx: Vec<f64>,
    pub fy: Vec<f64>,
    pub rep_k: f64,
}

impl ErSim {
    pub fn new(n: usize, lx: f64, ly: f64) -> Self {
        let mut s = Self {
            n,
            lx,
            ly,
            shear: 0.0,
            x: vec![0.0; n],
            y: vec![0.0; n],
            fx: vec![0.0; n],
            fy: vec![0.0; n],
            rep_k: 0.5,
        };
        for i in 0..n {
            let nx = (n as f64).sqrt().round() as usize;
            let sp = lx / nx as f64;
            s.x[i] = ((i % nx) as f64 + 0.5) * sp;
            s.y[i] = ((i / nx) as f64 + 0.5) * sp;
            // break the square-lattice symmetry (a stationary saddle of the
            // dipole potential) so chains can form
            let jx = ((i as f64 * 12.9898).sin() * 43758.5453).fract();
            let jy = ((i as f64 * 78.233).sin() * 43758.5453).fract();
            s.x[i] += (jx - 0.5) * 0.8;
            s.y[i] += (jy - 0.5) * 0.8;
        }
        s
    }

    /// Lees–Edwards minimum-image separation (with the box sheared by `shear`).
    fn separation(&self, i: usize, j: usize) -> (f64, f64) {
        let mut dx = self.x[j] - self.x[i];
        let mut dy = self.y[j] - self.y[i];
        dy -= self.ly * (dy / self.ly).round();
        // shear the x-image: images are displaced by shear * dy
        dx -= self.shear * dy;
        dx -= self.lx * (dx / self.lx).round();
        (dx, dy)
    }

    /// Sum pair forces (dipole + short-range repulsion) into fx, fy.
    pub fn forces(&mut self) {
        for v in self.fx.iter_mut() {
            *v = 0.0;
        }
        for v in self.fy.iter_mut() {
            *v = 0.0;
        }
        for i in 0..self.n {
            for j in (i + 1)..self.n {
                let (dx, dy) = self.separation(i, j);
                let r2 = dx * dx + dy * dy;
                let r = r2.sqrt();
                let (mut fx, mut fy) = dipole_force(dx, dy);
                if r < 2.0 {
                    let rep = self.rep_k * (2.0 / r).powi(24);
                    fx -= rep * dx;
                    fy -= rep * dy;
                }
                let cap = 500.0;
                let fmag = (fx * fx + fy * fy).sqrt();
                if fmag > cap {
                    let s = cap / fmag;
                    fx *= s;
                    fy *= s;
                }
                self.fx[i] += fx;
                self.fy[i] += fy;
                self.fx[j] -= fx;
                self.fy[j] -= fy;
            }
        }
    }

    pub fn step(&mut self, dt: f64) {
        self.forces();
        for i in 0..self.n {
            let vshear = self.shear * self.y[i];
            self.x[i] += (self.fx[i] + vshear) * dt;
            self.y[i] += self.fy[i] * dt;
            if self.y[i] >= self.ly {
                self.y[i] -= self.ly;
            } else if self.y[i] < 0.0 {
                self.y[i] += self.ly;
            }
            self.x[i] -= self.lx * (self.x[i] / self.lx).floor();
        }
        self.enforce_exclusion();
    }

    /// Project overlapping pairs apart, guaranteeing the sheared separation
    /// r ≥ 2a. The push is along the plain separation because the sheared
    /// separation has Jacobian [[1,−γ],[0,1]]: moving Δ by (2/r_s − 1)Δ scales
    /// the sheared separation by exactly 2/r_s.
    pub fn enforce_exclusion(&mut self) {
        for _ in 0..40 {
            let mut moved = false;
            for i in 0..self.n {
                for j in (i + 1)..self.n {
                    let (sx, sy) = self.separation(i, j);
                    let rs = (sx * sx + sy * sy).sqrt();
                    if rs < 2.0 && rs > 1e-9 {
                        // push vector D with J D = s, J = [[1,−γ],[0,1]]
                        let dxp = sx + self.shear * sy;
                        let dyp = sy;
                        let fac = 0.5 * (2.0 / rs - 1.0);
                        self.x[i] -= fac * dxp;
                        self.y[i] -= fac * dyp;
                        self.x[j] += fac * dxp;
                        self.y[j] += fac * dyp;
                        moved = true;
                    }
                }
            }
            if !moved {
                break;
            }
        }
    }

    /// Structural shear stress τ_xy = (1/V) Σ F_ij,x (y_j − y_i), in units of
    /// F₀/a² = ε₀ ε_f β² E², from the dipole pair forces (the hard-core
    /// constraint force is an artefact of the projection and is omitted).
    pub fn shear_stress(&self) -> f64 {
        let mut s = 0.0;
        for i in 0..self.n {
            for j in (i + 1)..self.n {
                let (dx, dy) = self.separation(i, j);
                let (fx, _) = dipole_force(dx, dy);
                s += fx * dy;
            }
        }
        s / (self.lx * self.ly)
    }

    /// Minimum pair separation (diagnostic).
    pub fn min_sep(&self) -> f64 {
        let mut m = f64::INFINITY;
        for i in 0..self.n {
            for j in (i + 1)..self.n {
                let (dx, dy) = self.separation(i, j);
                m = m.min((dx * dx + dy * dy).sqrt());
            }
        }
        m
    }

    /// Chain order parameter: fraction of particles with a neighbour within 3a
    /// whose direction is within 25° of the field.
    pub fn chain_order(&self) -> f64 {
        let mut good = 0;
        for i in 0..self.n {
            let mut ok = false;
            for j in 0..self.n {
                if i == j {
                    continue;
                }
                let (dx, dy) = self.separation(i, j);
                let r = (dx * dx + dy * dy).sqrt();
                if r < 3.0 && dy.abs() / r > 0.906 {
                    ok = true;
                    break;
                }
            }
            if ok {
                good += 1;
            }
        }
        good as f64 / self.n as f64
    }

    /// Average electrostatic (excluded-volume) potential energy per particle:
    /// U_dipole = (p²/4π ε₀ ε_f r³)(1 − 3cos²θ) plus a soft repulsion.
    pub fn energy_per_particle(&self) -> f64 {
        let mut u = 0.0;
        for i in 0..self.n {
            for j in (i + 1)..self.n {
                let (dx, dy) = self.separation(i, j);
                let r = (dx * dx + dy * dy).sqrt();
                if r < 1e-9 {
                    continue;
                }
                let cos = dy / r;
                let p2 = 4.0 * PI; // (p²/4π ε₀ ε_f) / F₀ a² = 4π
                u += p2 / (r * r * r) * (1.0 - 3.0 * cos * cos);
                if r < 2.0 {
                    u += self.rep_k / 23.0 * (2.0 / r).powi(23);
                }
            }
        }
        u / self.n as f64
    }
}

/// Point-dipole yield-stress scale (C = 1): τ_y = φ ε₀ ε_f β² E². The structure
/// factor C is measured by the simulation.
pub fn tau_y_scale(phi: f64) -> f64 {
    phi
}
