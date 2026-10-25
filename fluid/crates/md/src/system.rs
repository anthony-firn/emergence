//! The MD system: positions, velocities, forces, a periodic cubic box, a
//! cell-linked neighbour list, and velocity-Verlet integration.
//!
//! Reduced Lennard-Jones units are used throughout: length `sigma`, energy
//! `epsilon`, mass `m`, so temperature is in `epsilon/k_B`, pressure in
//! `epsilon/sigma^3` and time in `sigma sqrt(m/epsilon)`. Nothing in this
//! module knows about SI units; the [`crate::argon`] module converts.

use crate::neighbor::{build_neighbor_lists, min_image};
use crate::potential::{pair_energy, pair_force, Cutoff, Potential};
use crate::rng::Rng;
use rayon::prelude::*;

#[derive(Clone, Copy, Default)]
struct AtomAccum {
    force: [f64; 3],
    /// Half the sum of pair energies involving this atom.
    pe: f64,
    /// Sum of `-(u'/r) dr (x) dr` over neighbours; halved when totalled.
    stress: [f64; 6],
    /// Sum of `(F_ij . (v_i+v_j)) r_ij` over neighbours; halved when totalled.
    jv: [f64; 3],
}

/// A molecular-dynamics system in a periodic cubic box.
pub struct System {
    pub pos: Vec<[f64; 3]>,
    pub vel: Vec<[f64; 3]>,
    pub force: Vec<[f64; 3]>,
    pub mass: Vec<f64>,
    pub box_l: f64,
    pub potential: Potential,
    pub rc: f64,
    pub skin: f64,
    pub scheme: Cutoff,
    pub neighbors: Vec<Vec<usize>>,
    acc: Vec<AtomAccum>,
    pub rng: Rng,
    pub time: f64,
    pub steps: u64,
    disp: Vec<f64>,
    /// Nosé–Hoover friction coefficient (only used by the NH thermostat).
    pub nh_eta: f64,
}

impl System {
    /// Build a system from explicit positions. Velocities start at rest;
    /// call [`System::set_velocities_mb`].
    pub fn new(
        pos: Vec<[f64; 3]>,
        box_l: f64,
        potential: Potential,
        rc: f64,
        skin: f64,
        scheme: Cutoff,
        seed: u64,
    ) -> Self {
        let n = pos.len();
        assert!(n > 0);
        let neighbors = build_neighbor_lists(&pos, box_l, rc, skin);
        let mut sys = Self {
            pos,
            vel: vec![[0.0; 3]; n],
            force: vec![[0.0; 3]; n],
            mass: vec![1.0; n],
            box_l,
            potential,
            rc,
            skin,
            scheme,
            neighbors,
            acc: vec![AtomAccum::default(); n],
            rng: Rng::new(seed),
            time: 0.0,
            steps: 0,
            disp: vec![0.0; n],
            nh_eta: 0.0,
        };
        sys.compute_forces();
        sys
    }

    /// Perfect face-centred-cubic crystal: `n_cells^3` conventional cells,
    /// `4 n_cells^3` atoms, at the given number density.
    pub fn fcc(
        n_cells: usize,
        density: f64,
        potential: Potential,
        rc: f64,
        skin: f64,
        scheme: Cutoff,
        seed: u64,
    ) -> Self {
        let a = (4.0 / density).cbrt();
        let l = n_cells as f64 * a;
        let basis = [
            [0.0, 0.0, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.0],
        ];
        let mut pos = Vec::with_capacity(4 * n_cells * n_cells * n_cells);
        for i in 0..n_cells {
            for j in 0..n_cells {
                for k in 0..n_cells {
                    for b in &basis {
                        pos.push([
                            (i as f64 + b[0]) * a,
                            (j as f64 + b[1]) * a,
                            (k as f64 + b[2]) * a,
                        ]);
                    }
                }
            }
        }
        Self::new(pos, l, potential, rc, skin, scheme, seed)
    }

    pub fn n(&self) -> usize {
        self.pos.len()
    }

    pub fn volume(&self) -> f64 {
        self.box_l.powi(3)
    }

    /// Affinely scale the box and all positions by `mu` (isotropic
    /// barostat move). Velocities are left untouched; the barostat is meant
    /// to be paired with a thermostat.
    pub fn scale_box(&mut self, mu: f64) {
        self.box_l *= mu;
        for p in &mut self.pos {
            for a in 0..3 {
                p[a] *= mu;
            }
        }
        for d in &mut self.disp {
            *d *= mu;
        }
        // Forces and neighbour lists are invalidated by the affine move.
        self.rebuild_neighbors();
        self.compute_forces();
    }

    pub fn density(&self) -> f64 {
        self.n() as f64 / self.volume()
    }

    /// Instantaneous kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        self.vel
            .iter()
            .zip(&self.mass)
            .map(|(v, m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum()
    }

    /// Instantaneous temperature `2K/(3N)` (all atoms count).
    pub fn temperature(&self) -> f64 {
        2.0 * self.kinetic_energy() / (3.0 * self.n() as f64)
    }

    /// Temperature after removing the centre-of-mass kinetic energy
    /// (`3(N-1)` degrees of freedom), the value the canonical distribution
    /// of a periodic box actually has.
    pub fn temperature_com(&self) -> f64 {
        let n = self.n() as f64;
        if n <= 1.0 {
            return 0.0;
        }
        let ke = self.kinetic_energy() - self.com_kinetic_energy();
        2.0 * ke / (3.0 * (n - 1.0))
    }

    fn com_kinetic_energy(&self) -> f64 {
        let mut p = [0.0; 3];
        let mut m = 0.0;
        for (v, mi) in self.vel.iter().zip(&self.mass) {
            for a in 0..3 {
                p[a] += mi * v[a];
            }
            m += mi;
        }
        0.5 * (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) / m
    }

    /// Total configurational (potential) energy, including the analytic
    /// tail correction for [`Cutoff::TailCorrected`].
    pub fn potential_energy(&self) -> f64 {
        let mut e: f64 = self.acc.iter().map(|a| a.pe).sum();
        if self.scheme == Cutoff::TailCorrected {
            e += self.n() as f64 * self.potential.tail_energy_per_atom(self.density(), self.rc);
        }
        e
    }

    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy() + self.potential_energy()
    }

    /// Draw velocities from the Maxwell–Boltzmann distribution at `t`.
    pub fn set_velocities_mb(&mut self, t: f64) {
        for i in 0..self.n() {
            let m = self.mass[i];
            self.vel[i] = self.rng.maxwell_boltzmann(t, m);
        }
        self.remove_com_momentum();
    }

    /// Subtract the centre-of-mass velocity.
    pub fn remove_com_momentum(&mut self) {
        let mut p = [0.0; 3];
        let mut m = 0.0;
        for (v, mi) in self.vel.iter().zip(&self.mass) {
            for a in 0..3 {
                p[a] += mi * v[a];
            }
            m += mi;
        }
        for (v, _mi) in self.vel.iter_mut().zip(&self.mass) {
            for a in 0..3 {
                v[a] -= p[a] / m;
            }
        }
    }

    /// Rescale all velocities to exactly the target temperature (3N dof).
    pub fn rescale_to_temperature(&mut self, t: f64) {
        let k = self.kinetic_energy();
        if k <= 0.0 {
            return;
        }
        let target = 1.5 * self.n() as f64 * t;
        let lambda = (target / k).sqrt();
        for v in &mut self.vel {
            for a in 0..3 {
                v[a] *= lambda;
            }
        }
    }

    fn needs_rebuild(&self) -> bool {
        let half_skin = 0.5 * self.skin;
        self.disp.iter().any(|&d| d > half_skin)
    }

    /// Rebuild the neighbour lists and reset the displacement counters.
    pub fn rebuild_neighbors(&mut self) {
        self.neighbors = build_neighbor_lists(&self.pos, self.box_l, self.rc, self.skin);
        for d in &mut self.disp {
            *d = 0.0;
        }
    }

    /// Recompute forces, per-atom energies, stress and heat-flux virial.
    pub fn compute_forces(&mut self) {
        if self.needs_rebuild() {
            self.rebuild_neighbors();
        }
        let rc = self.rc;
        let rc2 = rc * rc;
        let l = self.box_l;
        let scheme = self.scheme;
        let potential = &self.potential;
        let n = self.n();
        let System {
            pos,
            vel,
            neighbors,
            acc,
            ..
        } = self;
        acc.par_iter_mut().enumerate().for_each(|(i, a)| {
            let pi = pos[i];
            let vi = vel[i];
            let mut force = [0.0; 3];
            let mut pe = 0.0;
            let mut stress = [0.0; 6];
            let mut jv = [0.0; 3];
            for &j in &neighbors[i] {
                let dr = min_image([pos[j][0] - pi[0], pos[j][1] - pi[1], pos[j][2] - pi[2]], l);
                let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                if r2 > rc2 {
                    continue;
                }
                let r = r2.sqrt();
                let u = pair_energy(potential, r, rc, scheme);
                let fmag = pair_force(potential, r, rc, scheme);
                let inv_r = 1.0 / r;
                let f = [
                    fmag * dr[0] * inv_r,
                    fmag * dr[1] * inv_r,
                    fmag * dr[2] * inv_r,
                ];
                force[0] += f[0];
                force[1] += f[1];
                force[2] += f[2];
                pe += u;
                let c = -fmag * inv_r;
                stress[0] += c * dr[0] * dr[0];
                stress[1] += c * dr[1] * dr[1];
                stress[2] += c * dr[2] * dr[2];
                stress[3] += c * dr[0] * dr[1];
                stress[4] += c * dr[0] * dr[2];
                stress[5] += c * dr[1] * dr[2];
                let vj = vel[j];
                let fdotv =
                    f[0] * (vi[0] + vj[0]) + f[1] * (vi[1] + vj[1]) + f[2] * (vi[2] + vj[2]);
                // r_ij = r_i - r_j = -dr
                jv[0] -= fdotv * dr[0];
                jv[1] -= fdotv * dr[1];
                jv[2] -= fdotv * dr[2];
            }
            *a = AtomAccum {
                force,
                pe: 0.5 * pe,
                stress,
                jv,
            };
        });
        for (f, a) in self.force.iter_mut().zip(&self.acc) {
            *f = a.force;
        }
        let _ = n;
    }

    /// Instantaneous pressure tensor `[xx, yy, zz, xy, xz, yz]` in reduced
    /// units (pressure = force/area).
    pub fn pressure_tensor(&self) -> [f64; 6] {
        let mut p = [0.0; 6];
        for (v, m) in self.vel.iter().zip(&self.mass) {
            p[0] += m * v[0] * v[0];
            p[1] += m * v[1] * v[1];
            p[2] += m * v[2] * v[2];
            p[3] += m * v[0] * v[1];
            p[4] += m * v[0] * v[2];
            p[5] += m * v[1] * v[2];
        }
        let mut s = [0.0; 6];
        for a in &self.acc {
            for (sk, ak) in s.iter_mut().zip(&a.stress) {
                *sk += 0.5 * ak;
            }
        }
        for (pk, sk) in p.iter_mut().zip(&s) {
            *pk += *sk;
        }
        if self.scheme == Cutoff::TailCorrected {
            // `tail` is a pressure; scale by V before dividing the
            // extensive tensor by V.
            let tail = self.potential.tail_pressure(self.density(), self.rc) * self.volume();
            p[0] += tail;
            p[1] += tail;
            p[2] += tail;
        }
        let v = self.volume();
        std::array::from_fn(|i| p[i] / v)
    }

    /// Diagnostic: trace of the internal pair virial, `0.5 * sum_i
    /// tr(stress_i)`. Should equal `-sum_{i<j} u'(r) r`.
    pub fn internal_pair_virial(&self) -> f64 {
        0.5 * self
            .acc
            .iter()
            .map(|a| a.stress[0] + a.stress[1] + a.stress[2])
            .sum::<f64>()
    }

    /// Scalar pressure `tr(P)/3`.
    pub fn pressure(&self) -> f64 {
        let p = self.pressure_tensor();
        (p[0] + p[1] + p[2]) / 3.0
    }

    /// Volume-integrated heat flux `J V` (no `1/V` factor). The thermal
    /// conductivity follows from `kappa = (1/(V k_B T^2)) int <JV(0)JV(t)>`.
    ///
    /// `acc[i].jv` already symmetrises to `(F_ij . (v_i+v_j)) r_ij`; summing
    /// it over all atoms double counts each pair, and the Irving–Kirkwood
    /// pair flux carries `1/2`, so the total weight is `1/4`.
    pub fn heat_flux(&self) -> [f64; 3] {
        let mut j = [0.0; 3];
        for i in 0..self.n() {
            let v = self.vel[i];
            let m = self.mass[i];
            let e = 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) + self.acc[i].pe;
            for a in 0..3 {
                j[a] += e * v[a] + 0.25 * self.acc[i].jv[a];
            }
        }
        j
    }

    /// One velocity-Verlet step of the microcanonical equations of motion.
    pub fn step_nve(&mut self, dt: f64) {
        let n = self.n();
        for i in 0..n {
            let half = 0.5 * dt / self.mass[i];
            for a in 0..3 {
                self.vel[i][a] += half * self.force[i][a];
            }
        }
        for i in 0..n {
            for a in 0..3 {
                let d = dt * self.vel[i][a];
                self.pos[i][a] += d;
            }
            self.disp[i] += dt
                * (self.vel[i][0] * self.vel[i][0]
                    + self.vel[i][1] * self.vel[i][1]
                    + self.vel[i][2] * self.vel[i][2])
                    .sqrt();
        }
        self.compute_forces();
        for i in 0..n {
            let half = 0.5 * dt / self.mass[i];
            for a in 0..3 {
                self.vel[i][a] += half * self.force[i][a];
            }
        }
        self.time += dt;
        self.steps += 1;
    }

    /// One Nosé–Hoover step for a single thermostat chain of length one.
    /// The friction is integrated with a symmetric split around the
    /// velocity-Verlet step (Allen & Tildesley, *Computer Simulation of
    /// Liquids*, ch. 6).
    pub fn step_nose_hoover(&mut self, dt: f64, target: f64, tau: f64) {
        let nf = 3.0 * self.n() as f64;
        let q = nf * target * tau * tau;
        let g0 = 2.0 * self.kinetic_energy() - nf * target;
        let eta_half = self.nh_eta + 0.5 * dt * g0 / q;
        let n = self.n();
        for i in 0..n {
            let half = 0.5 * dt;
            for a in 0..3 {
                self.vel[i][a] +=
                    half * (self.force[i][a] / self.mass[i] - eta_half * self.vel[i][a]);
            }
        }
        for i in 0..n {
            for a in 0..3 {
                self.pos[i][a] += dt * self.vel[i][a];
            }
            self.disp[i] += dt
                * (self.vel[i][0] * self.vel[i][0]
                    + self.vel[i][1] * self.vel[i][1]
                    + self.vel[i][2] * self.vel[i][2])
                    .sqrt();
        }
        self.compute_forces();
        let g1 = 2.0 * self.kinetic_energy() - nf * target;
        let eta_new = eta_half + 0.5 * dt * g1 / q;
        for i in 0..n {
            let half = 0.5 * dt;
            for a in 0..3 {
                self.vel[i][a] +=
                    half * (self.force[i][a] / self.mass[i] - eta_new * self.vel[i][a]);
            }
        }
        self.nh_eta = eta_new;
        self.time += dt;
        self.steps += 1;
    }

    /// Radial distribution function on `[0, rmax]` with `nbins` bins,
    /// normalised so that `g -> 1` in a homogeneous fluid. Uses min-image
    /// distances over all unordered pairs.
    pub fn rdf(&self, rmax: f64, nbins: usize) -> Vec<f64> {
        let n = self.n();
        let dr = rmax / nbins as f64;
        let mut hist = vec![0.0f64; nbins];
        for i in 0..n {
            for j in (i + 1)..n {
                let d = min_image(
                    [
                        self.pos[i][0] - self.pos[j][0],
                        self.pos[i][1] - self.pos[j][1],
                        self.pos[i][2] - self.pos[j][2],
                    ],
                    self.box_l,
                );
                let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                if r < rmax {
                    let b = (r / dr) as usize;
                    hist[b.min(nbins - 1)] += 1.0;
                }
            }
        }
        let rho = self.density();
        let norm = n as f64 * (n as f64 - 1.0) / self.volume();
        let mut g = vec![0.0; nbins];
        for b in 0..nbins {
            let r = (b as f64 + 0.5) * dr;
            // unordered pairs -> factor 2.
            g[b] = 2.0 * hist[b] / (norm * 4.0 * std::f64::consts::PI * r * r * dr);
        }
        let _ = rho;
        g
    }
}
