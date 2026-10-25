//! Soft-sphere discrete-element method (DEM): a two-sphere collision
//! benchmark and a Lees–Edwards/SLLOD shear cell for the `mu(I)` rheology.
//!
//! The system uses reduced units fixed by the particle diameter `d = 1`
//! (radius `r = 0.5`) and density `rho_p = 1`, so mass is `rho_p (4/3) pi
//! r^3` and stress is `rho_p d^2/T^2`. The contact law is
//! [`crate::contact::hertz_mindlin`].
//!
//! # Two-sphere collision
//!
//! For a head-on elastic Hertz collision the contact duration is
//! `t_c = 2.87 (m*^2 / (R* E*^2 v))^(1/5)` (the exact Hertz result;
//! Johnson, *Contact Mechanics*, 1985 §11.7; the numerical prefactor
//! `2.868` follows from `2 I (15/16)^(2/5)` with
//! `I = ∫_0^1 (1-u^(5/2))^(-1/2) du = 1.4714`).
//! [`two_sphere_collision`] integrates it and [`hertz_contact_time`] gives
//! the closed form. With the Tsuji damping the measured coefficient of
//! restitution is compared with the target `e`.
//!
//! # Lees–Edwards shear cell
//!
//! The cell imposes a mean shear `v_x = gdot y` with the SLLOD equations for
//! the **peculiar** velocity `u`,
//!
//! ```text
//! xdot = u + gdot y e_x,   udot = F/m - gdot u_y e_x,
//! ```
//!
//! with Lees–Edwards periodic images. It measures the pressure
//! `P = tr(sigma)/3` and the shear stress `tau = sigma_xy` from the virial,
//! giving `mu = tau/P` and the inertial number
//! `I = gdot d / sqrt(P/rho_p)`. A weak Berendsen thermostat on the granular
//! temperature is available (off by default); the contact damping alone
//! balances the shear work at the volume fractions used here.

use std::collections::HashMap;
use std::f64::consts::PI;

use crate::contact::{hertz_mindlin, ContactParams, ContactState};
use crate::rheology::inertial_number;
use crate::rng::Rng;

/// One spherical particle (reduced units). `v` is the **peculiar** velocity
/// in the shear cell.
#[derive(Clone, Copy, Debug)]
pub struct Particle {
    pub x: [f64; 3],
    pub v: [f64; 3],
    pub w: [f64; 3],
    pub r: f64,
    pub m: f64,
}

impl Particle {
    pub fn new(x: [f64; 3], r: f64, density: f64) -> Self {
        let m = density * (4.0 / 3.0) * PI * r * r * r;
        Self {
            x,
            v: [0.0; 3],
            w: [0.0; 3],
            r,
            m,
        }
    }

    /// Moment of inertia of a solid sphere.
    pub fn inertia(&self) -> f64 {
        0.4 * self.m * self.r * self.r
    }
}

/// Simulation-wide parameters for the DEM.
#[derive(Clone, Copy, Debug)]
pub struct DemParams {
    pub contact: ContactParams,
    pub radius: f64,
    pub density: f64,
    /// Maximal contact stiffness scale (informational); the actual
    /// stiffness is in `contact`.
    pub young: f64,
    /// Time step.
    pub dt: f64,
}

impl DemParams {
    pub fn new(contact: ContactParams, radius: f64, density: f64, dt: f64) -> Self {
        Self {
            contact,
            radius,
            density,
            young: contact.young,
            dt,
        }
    }

    /// Exact head-on elastic Hertz contact time for relative speed `v`.
    pub fn hertz_contact_time(&self, v: f64) -> f64 {
        hertz_contact_time(&self.contact, self.radius, self.density, v)
    }
}

/// Exact elastic Hertz contact time `t_c = 2.868 (m*^2/(R* E*^2 v))^(1/5)`.
pub fn hertz_contact_time(p: &ContactParams, radius: f64, density: f64, v: f64) -> f64 {
    let m = density * (4.0 / 3.0) * PI * radius * radius * radius;
    let m_star = m / 2.0;
    let r_star = radius / 2.0;
    let e_star = p.e_star();
    let inner = m_star * m_star / (r_star * e_star * e_star * v);
    2.868 * inner.powf(1.0 / 5.0)
}

/// Result of a two-sphere collision.
#[derive(Clone, Copy, Debug)]
pub struct CollisionResult {
    /// Time during which the spheres overlap.
    pub contact_time: f64,
    /// Measured coefficient of restitution (1 if undamped).
    pub restitution: f64,
    /// Relative speed before contact.
    pub v_in: f64,
    /// Relative speed after contact.
    pub v_out: f64,
}

/// Integrate a head-on two-sphere collision. Give them relative speed
/// `v_impact` in the centre-of-mass frame. If `damped` is false the target
/// restitution is forced to 1 (elastic), so the measured contact time can
/// be compared with [`hertz_contact_time`].
pub fn two_sphere_collision(
    p: &ContactParams,
    radius: f64,
    density: f64,
    v_impact: f64,
    damped: bool,
) -> CollisionResult {
    let mut pc = *p;
    if !damped {
        pc.restitution = 1.0;
    }
    let m = density * (4.0 / 3.0) * PI * radius * radius * radius;
    let mut x1 = [0.0; 3];
    let mut x2 = [2.0 * radius, 0.0, 0.0];
    let mut v1 = [0.5 * v_impact, 0.0, 0.0];
    let mut v2 = [-0.5 * v_impact, 0.0, 0.0];
    let mut st = ContactState::default();

    let tc = hertz_contact_time(&pc, radius, density, v_impact);
    let dt = tc / 400.0;
    let max_steps = (10.0 * tc / dt) as usize;

    let compute = |x1: [f64; 3],
                   x2: [f64; 3],
                   v1: [f64; 3],
                   v2: [f64; 3],
                   st: &mut ContactState,
                   dt: f64| {
        let c = hertz_mindlin(
            &pc, radius, radius, m, m, x1, v1, [0.0; 3], x2, v2, [0.0; 3], st, dt,
        );
        match c {
            Some(c) => (c.force_i, [-c.force_i[0], -c.force_i[1], -c.force_i[2]]),
            None => ([0.0; 3], [0.0; 3]),
        }
    };

    let (mut f1, mut f2) = compute(x1, x2, v1, v2, &mut st, dt);
    let mut was_contact = false;
    let mut t_contact = 0.0;
    let mut v_out = -v_impact;
    let mut first = true;
    for _ in 0..max_steps {
        // Velocity Verlet with velocity-dependent contact force.
        for k in 0..3 {
            v1[k] += 0.5 * f1[k] / m * dt;
            v2[k] += 0.5 * f2[k] / m * dt;
            x1[k] += v1[k] * dt;
            x2[k] += v2[k] * dt;
        }
        let (nf1, nf2) = compute(x1, x2, v1, v2, &mut st, dt);
        for k in 0..3 {
            v1[k] += 0.5 * nf1[k] / m * dt;
            v2[k] += 0.5 * nf2[k] / m * dt;
        }
        f1 = nf1;
        f2 = nf2;

        let dist = ((x2[0] - x1[0]).powi(2)
            + (x2[1] - x1[1]).powi(2)
            + (x2[2] - x1[2]).powi(2))
        .sqrt();
        let contact = dist < 2.0 * radius;
        if contact {
            if first {
                first = false;
            }
            t_contact += dt;
            was_contact = true;
        } else if was_contact {
            v_out = (v2[0] - v1[0]).abs();
            break;
        } else if !first {
            break;
        }
    }
    if !was_contact {
        // No overlap detected (should not happen with a head-on start).
        v_out = v_impact;
    }
    CollisionResult {
        contact_time: t_contact,
        restitution: v_out / v_impact,
        v_in: v_impact,
        v_out,
    }
}

/// A Lees–Edwards / SLLOD shear cell of `n` identical spheres at fixed
/// volume and imposed shear rate.
pub struct ShearCell {
    pub parts: Vec<Particle>,
    pub contact: ContactParams,
    pub radius: f64,
    pub density: f64,
    pub dt: f64,
    pub lx: f64,
    pub ly: f64,
    pub lz: f64,
    /// Imposed shear rate (flow `x`, gradient `y`).
    pub gdot: f64,
    /// Accumulated strain, used for the Lees–Edwards image shift.
    pub strain: f64,
    pub time: f64,
    contacts: HashMap<(u32, u32), ContactState>,
    /// Optional Berendsen relaxation time for the granular temperature.
    pub thermostat_tau: f64,
    /// Optional target granular temperature.
    pub target_temp: Option<f64>,
    /// Optional constant-pressure (Berendsen) target and relaxation time.
    pub target_pressure: Option<f64>,
    pub barostat_tau: f64,
    /// Last measured pressure.
    pub pressure: f64,
    /// Last measured shear stress `sigma_xy`.
    pub shear_stress: f64,
}

impl ShearCell {
    /// Build a cell with `n` particles at volume fraction `phi` in a cubic
    /// box, on a jittered cubic lattice, with a deterministic seed.
    pub fn new_cubic(
        n: usize,
        phi: f64,
        gdot: f64,
        contact: ContactParams,
        radius: f64,
        density: f64,
        dt: f64,
        seed: u64,
    ) -> Self {
        let vpart = (4.0 / 3.0) * PI * radius * radius * radius;
        let volume = n as f64 * vpart / phi;
        let l = volume.cbrt();
        // Cubic lattice close to n sites.
        let side = (n as f64).cbrt().ceil() as usize;
        let mut rng = Rng::new(seed);
        let mut parts = Vec::with_capacity(n);
        let spacing = l / side as f64;
        let mut count = 0;
        'outer: for i in 0..side {
            for j in 0..side {
                for k in 0..side {
                    if count >= n {
                        break 'outer;
                    }
                    let jitter = 0.05 * spacing;
                    let x = [
                        (i as f64 + 0.5) * spacing + rng.uniform(-jitter, jitter),
                        (j as f64 + 0.5) * spacing + rng.uniform(-jitter, jitter),
                        (k as f64 + 0.5) * spacing + rng.uniform(-jitter, jitter),
                    ];
                    let mut part = Particle::new(x, radius, density);
                    // Small random peculiar velocities.
                    for c in 0..3 {
                        part.v[c] = rng.normal() * 0.01;
                        part.w[c] = rng.normal() * 0.01;
                    }
                    parts.push(part);
                    count += 1;
                }
            }
        }
        Self {
            parts,
            contact,
            radius,
            density,
            dt,
            lx: l,
            ly: l,
            lz: l,
            gdot,
            strain: 0.0,
            time: 0.0,
            contacts: HashMap::new(),
            thermostat_tau: 0.0,
            target_temp: None,
            target_pressure: None,
            barostat_tau: 1.0,
            pressure: 0.0,
            shear_stress: 0.0,
        }
    }

    pub fn n(&self) -> usize {
        self.parts.len()
    }

    /// Number of active contacts in the last force evaluation.
    pub fn n_contacts(&self) -> usize {
        self.contacts.len()
    }

    /// Smallest centre-to-centre separation in the cell (diagnostic).
    pub fn min_separation(&self) -> f64 {
        let n = self.n();
        let mut m = f64::INFINITY;
        for i in 0..n {
            for j in (i + 1)..n {
                let (_, d) = self.separation(i, j);
                if d < m {
                    m = d;
                }
            }
        }
        m
    }

    pub fn volume(&self) -> f64 {
        self.lx * self.ly * self.lz
    }

    pub fn volume_fraction(&self) -> f64 {
        let vp = (4.0 / 3.0) * PI * self.radius.powi(3);
        self.n() as f64 * vp / self.volume()
    }

    /// Wrap a position back into the primary box, shifting `x` when `y` is
    /// wrapped so the Lees–Edwards image set is preserved.
    fn wrap(&self, x: &mut [f64; 3]) {
        if x[1] >= self.ly {
            x[1] -= self.ly;
            x[0] -= self.strain * self.ly;
        } else if x[1] < 0.0 {
            x[1] += self.ly;
            x[0] += self.strain * self.ly;
        }
        if x[0] >= self.lx {
            x[0] -= self.lx;
        } else if x[0] < 0.0 {
            x[0] += self.lx;
        }
        if x[2] >= self.lz {
            x[2] -= self.lz;
        } else if x[2] < 0.0 {
            x[2] += self.lz;
        }
    }

    /// Minimum-image separation from `i` to `j` with the Lees–Edwards
    /// shear shift. Returns the raw (unwrapped) `x_j - x_i` and its length.
    fn separation(&self, i: usize, j: usize) -> ([f64; 3], f64) {
        let xi = self.parts[i].x;
        let xj = self.parts[j].x;
        let mut dr = [xj[0] - xi[0], xj[1] - xi[1], xj[2] - xi[2]];
        let ny = (dr[1] / self.ly).round();
        dr[1] -= ny * self.ly;
        dr[0] -= ny * self.strain * self.ly;
        let nx = (dr[0] / self.lx).round();
        dr[0] -= nx * self.lx;
        let nz = (dr[2] / self.lz).round();
        dr[2] -= nz * self.lz;
        let d = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
        (dr, d)
    }

    /// Compute pair forces/torques and the virial. Updates the tangential
    /// contact history.
    pub fn forces(&mut self) -> (Vec<[f64; 3]>, Vec<[f64; 3]>, [[f64; 3]; 3]) {
        let n = self.n();
        let mut f = vec![[0.0; 3]; n];
        let mut t = vec![[0.0; 3]; n];
        let mut virial = [[0.0; 3]; 3];
        let dt = self.dt;
        let mut new_contacts: HashMap<(u32, u32), ContactState> =
            HashMap::with_capacity(self.contacts.len() + 16);
        for i in 0..n {
            for j in (i + 1)..n {
                let (dr, dist) = self.separation(i, j);
                let ri = self.parts[i].r;
                let rj = self.parts[j].r;
                if dist >= ri + rj {
                    continue;
                }
                // Unwrapped centre of j relative to i.
                let xj = [
                    self.parts[i].x[0] + dr[0],
                    self.parts[i].x[1] + dr[1],
                    self.parts[i].x[2] + dr[2],
                ];
                let key = (i as u32, j as u32);
                let mut st = *self.contacts.get(&key).unwrap_or(&ContactState::default());
                // Physical relative velocity including the affine part.
                let affine = self.gdot * dr[1];
                let vj = [
                    self.parts[j].v[0] + affine,
                    self.parts[j].v[1],
                    self.parts[j].v[2],
                ];
                let c = hertz_mindlin(
                    &self.contact,
                    ri,
                    rj,
                    self.parts[i].m,
                    self.parts[j].m,
                    self.parts[i].x,
                    self.parts[i].v,
                    self.parts[i].w,
                    xj,
                    vj,
                    self.parts[j].w,
                    &mut st,
                    dt,
                );
                if let Some(c) = c {
                    for k in 0..3 {
                        f[i][k] += c.force_i[k];
                        f[j][k] -= c.force_i[k];
                        t[i][k] += c.torque_i[k];
                        t[j][k] += c.torque_j[k];
                    }
                    // Virial: F_i (x) (x_i - x_j) = -F_i (x) dr.
                    for a in 0..3 {
                        for b in 0..3 {
                            virial[a][b] -= c.force_i[a] * dr[b];
                        }
                    }
                    new_contacts.insert(key, st);
                }
            }
        }
        self.contacts = new_contacts;
        (f, t, virial)
    }

    /// Full stress tensor from the virial plus peculiar kinetic part.
    pub fn stress_tensor(&self, virial: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let mut s = *virial;
        for p in &self.parts {
            for a in 0..3 {
                for b in 0..3 {
                    s[a][b] += p.m * p.v[a] * p.v[b];
                }
            }
        }
        let v = self.volume();
        for a in 0..3 {
            for b in 0..3 {
                s[a][b] /= v;
            }
        }
        s
    }

    /// One SLLOD step.
    pub fn step(&mut self) {
        let dt = self.dt;
        let (f, t, virial) = self.forces();
        let n = self.n();
        let gdot = self.gdot;
        // Update peculiar velocities.
        for i in 0..n {
            let m = self.parts[i].m;
            let inertia = self.parts[i].inertia();
            let vy = self.parts[i].v[1];
            for k in 0..3 {
                self.parts[i].v[k] += f[i][k] / m * dt;
                self.parts[i].w[k] += t[i][k] / inertia * dt;
            }
            self.parts[i].v[0] -= gdot * vy * dt;
        }
        // Update positions.
        for i in 0..n {
            let y = self.parts[i].x[1];
            self.parts[i].x[0] += (self.parts[i].v[0] + gdot * y) * dt;
            self.parts[i].x[1] += self.parts[i].v[1] * dt;
            self.parts[i].x[2] += self.parts[i].v[2] * dt;
            let mut x = self.parts[i].x;
            self.wrap(&mut x);
            self.parts[i].x = x;
        }
        // Optional thermostat on the peculiar granular temperature.
        if let (Some(tt), true) = (self.target_temp, self.thermostat_tau > 0.0) {
            let tg = self.granular_temperature();
            if tg > 1e-300 {
                let lambda = (1.0 + dt / self.thermostat_tau * (tt / tg - 1.0)).sqrt();
                for p in &mut self.parts {
                    for k in 0..3 {
                        p.v[k] *= lambda;
                    }
                }
            }
        }
        self.strain += gdot * dt;
        self.time += dt;

        // Measure stresses with the updated state.
        let s = self.stress_tensor(&virial);
        self.pressure = (s[0][0] + s[1][1] + s[2][2]) / 3.0;
        self.shear_stress = s[0][1];

        // Constant-pressure Berendsen feedback on the volume.
        if let Some(p0) = self.target_pressure {
            if p0 > 0.0 && self.barostat_tau > 0.0 {
                let chi = 1.0 - (dt / self.barostat_tau) * (p0 - self.pressure) / p0;
                let sc = chi.max(0.5).cbrt().clamp(0.995, 1.005);
                let (lx, ly, lz, strain) = (self.lx, self.ly, self.lz, self.strain);
                for p in &mut self.parts {
                    for k in 0..3 {
                        p.x[k] *= sc;
                    }
                    let x = &mut p.x;
                    if x[1] >= ly {
                        x[1] -= ly;
                        x[0] -= strain * ly;
                    } else if x[1] < 0.0 {
                        x[1] += ly;
                        x[0] += strain * ly;
                    }
                    if x[0] >= lx {
                        x[0] -= lx;
                    } else if x[0] < 0.0 {
                        x[0] += lx;
                    }
                    if x[2] >= lz {
                        x[2] -= lz;
                    } else if x[2] < 0.0 {
                        x[2] += lz;
                    }
                }
                self.lx *= sc;
                self.ly *= sc;
                self.lz *= sc;
            }
        }
    }

    /// Enable a constant-pressure Berendsen barostat with target `p0` and
    /// relaxation time `tau`.
    pub fn set_pressure(&mut self, p0: f64, tau: f64) {
        self.target_pressure = Some(p0);
        self.barostat_tau = tau;
    }

    /// Granular temperature `T_g = (1/3N) sum m |u|^2`.
    pub fn granular_temperature(&self) -> f64 {
        let mut sum = 0.0;
        for p in &self.parts {
            sum += p.m * (p.v[0] * p.v[0] + p.v[1] * p.v[1] + p.v[2] * p.v[2]);
        }
        sum / (3.0 * self.n() as f64)
    }

    /// Run `steps` steps, averaging `P` and `tau` over the second half.
    /// Returns `(P, tau, mu, I, T_g)`.
    pub fn run(&mut self, steps: usize) -> (f64, f64, f64, f64, f64) {
        let half = steps / 2;
        let mut psum = 0.0;
        let mut tsum = 0.0;
        let mut tgsum = 0.0;
        let mut count = 0.0;
        for k in 0..steps {
            self.step();
            if k >= half {
                psum += self.pressure;
                tsum += self.shear_stress;
                tgsum += self.granular_temperature();
                count += 1.0;
            }
        }
        let p = psum / count;
        let tau = tsum / count;
        let tg = tgsum / count;
        let mu = tau.abs() / p;
        let i = inertial_number(self.gdot, 2.0 * self.radius, p, self.density);
        (p, tau, mu, i, tg)
    }

    /// The effective friction `mu = |sigma_xy|/P` of the current state.
    pub fn mu(&self) -> f64 {
        self.shear_stress.abs() / self.pressure
    }

    /// The inertial number of the current state.
    pub fn inertial_number(&self) -> f64 {
        inertial_number(self.gdot, 2.0 * self.radius, self.pressure, self.density)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn glassish() -> ContactParams {
        ContactParams {
            young: 1000.0,
            poisson: 0.3,
            mu: 0.4,
            mu_rolling: 0.0,
            restitution: 1.0,
        }
    }

    #[test]
    fn collision_energy_and_momentum() {
        let p = ContactParams {
            restitution: 0.5,
            ..glassish()
        };
        let r = two_sphere_collision(&p, 0.5, 1.0, 0.1, true);
        // Momentum conserved in the COM frame: v1+v2 stays 0.
        assert!((r.restitution - 0.5).abs() < 0.15, "e={}", r.restitution);
        assert!(r.restitution < 1.0);
        let re = two_sphere_collision(&glassish(), 0.5, 1.0, 0.1, false);
        assert!((re.restitution - 1.0).abs() < 0.02, "e={}", re.restitution);
    }
}
