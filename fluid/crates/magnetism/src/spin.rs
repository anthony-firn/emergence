//! Atomistic classical Heisenberg spins: Metropolis Monte Carlo, the Binder
//! cumulant and Curie temperature, linear spin-wave theory and a spin-chain
//! Landau-Lifshitz-Gilbert integrator.
//!
//! The Hamiltonian used throughout is
//!
//! ```text
//! H = -J sum_<ij> S_i . S_j - K sum_i (S_i . e)^2 - H_ext . sum_i S_i,
//! ```
//!
//! with unit classical spins (`|S_i| = 1`).  The mean-field Curie
//! temperature is `k_B T_C = J z S^2 / 3` with `S = 1` for unit spins; the
//! coordination `z` is 6, 8, 12 for sc, bcc, fcc.
//!
//! Monte-Carlo reference: Chen, Ferrenberg & Landau, *Phys. Rev. B* **48**,
//! 3249 (1993) give `k_B T_C / J = 1.44299(14)` for the classical
//! nearest-neighbour Heisenberg model on the simple-cubic lattice.  The
//! Binder cumulant `U_4 = 1 - <M^4>/(3 <M^2>^2)` crosses between system
//! sizes at `T_C`.

use crate::geometry::{dot, norm, LatticeKind, V3};
use crate::rng::Rng;

pub struct BinderResult {
    pub u4: f64,
    pub m2: f64,
    pub m4: f64,
    pub samples: usize,
}

pub struct McConfig {
    pub equil_sweeps: usize,
    pub samples: usize,
    pub sweeps_per_sample: usize,
}

impl Default for McConfig {
    fn default() -> Self {
        Self {
            equil_sweeps: 2000,
            samples: 400,
            sweeps_per_sample: 4,
        }
    }
}

/// A finite `l x l x l` conventional-cell cluster with periodic boundaries.
pub struct Heisenberg {
    pub kind: LatticeKind,
    pub l: usize,
    pub a: f64,
    pub j: f64,
    pub k: f64,
    pub easy: V3,
    pub h: V3,
    pub spins: Vec<V3>,
    n_sites: usize,
    neighbors: Vec<Vec<usize>>,
}

impl Heisenberg {
    pub fn new(
        kind: LatticeKind,
        l: usize,
        a: f64,
        j: f64,
        k: f64,
        easy: V3,
        h: V3,
        seed: u64,
    ) -> Self {
        let n_basis = kind.n_atoms();
        let n_sites = l * l * l * n_basis;
        let mut init_rng = Rng::new(seed ^ 0xABCD);
        let spins = (0..n_sites).map(|_| init_rng.unit_vector()).collect();
        let neigh = kind.neighbors(a);
        let neighbors = build_neighbor_indices(kind, l, &neigh, n_basis, n_sites);
        Self {
            kind,
            l,
            a,
            j,
            k,
            easy: normalize(easy),
            h,
            spins,
            n_sites,
            neighbors,
        }
    }

    pub fn n_sites(&self) -> usize {
        self.n_sites
    }

    /// Total energy of the current configuration.
    pub fn energy(&self) -> f64 {
        let mut e = 0.0;
        for i in 0..self.n_sites {
            for &j in &self.neighbors[i] {
                e -= 0.5 * self.j * dot(self.spins[i], self.spins[j]);
            }
            e -= self.k * dot(self.spins[i], self.easy).powi(2);
            e -= dot(self.h, self.spins[i]);
        }
        e
    }

    /// One Metropolis sweep (N single-spin proposals with uniform
    /// proposals on the unit sphere, which is symmetric and satisfies
    /// detailed balance).
    pub fn sweep(&mut self, temperature: f64, rng: &mut Rng) {
        for _ in 0..self.n_sites {
            let i = (rng.uniform() * self.n_sites as f64) as usize % self.n_sites;
            let old = self.spins[i];
            let new = rng.unit_vector();
            let d = [new[0] - old[0], new[1] - old[1], new[2] - old[2]];
            let mut de = 0.0;
            for &j in &self.neighbors[i] {
                de -= self.j * dot(d, self.spins[j]);
            }
            de -= self.k
                * (dot(new, self.easy).powi(2) - dot(old, self.easy).powi(2));
            de -= dot(self.h, d);
            if de <= 0.0 || rng.uniform() < (-de / temperature).exp() {
                self.spins[i] = new;
            }
        }
    }

    pub fn run(&mut self, temperature: f64, sweeps: usize, rng: &mut Rng) {
        for _ in 0..sweeps {
            self.sweep(temperature, rng);
        }
    }

    /// Order parameter magnitude `M = |(1/N) sum S_i|`.
    pub fn magnetization(&self) -> f64 {
        let mut m = [0.0; 3];
        for s in &self.spins {
            m[0] += s[0];
            m[1] += s[1];
            m[2] += s[2];
        }
        norm(m) / self.n_sites as f64
    }

    /// Binder cumulant `U_4 = 1 - <M^4>/(3 <M^2>^2)` and its moments.
    pub fn binder_cumulant(
        &mut self,
        temperature: f64,
        cfg: &McConfig,
        rng: &mut Rng,
    ) -> BinderResult {
        self.run(temperature, cfg.equil_sweeps, rng);
        let mut sum2 = 0.0;
        let mut sum4 = 0.0;
        for _ in 0..cfg.samples {
            self.run(temperature, cfg.sweeps_per_sample, rng);
            let m = self.magnetization();
            let m2 = m * m;
            sum2 += m2;
            sum4 += m2 * m2;
        }
        let m2 = sum2 / cfg.samples as f64;
        let m4 = sum4 / cfg.samples as f64;
        BinderResult {
            u4: 1.0 - m4 / (3.0 * m2 * m2),
            m2,
            m4,
            samples: cfg.samples,
        }
    }
}

fn normalize(v: V3) -> V3 {
    let n = norm(v).max(1e-300);
    [v[0] / n, v[1] / n, v[2] / n]
}

fn build_neighbor_indices(
    _kind: LatticeKind,
    l: usize,
    neigh: &[Vec<crate::geometry::Neighbor>],
    n_basis: usize,
    n_sites: usize,
) -> Vec<Vec<usize>> {
    let mut out = vec![Vec::new(); n_sites];
    for cx in 0..l {
        for cy in 0..l {
            for cz in 0..l {
                for bi in 0..n_basis {
                    let i = ((cx * l + cy) * l + cz) * n_basis + bi;
                    for nbr in &neigh[bi] {
                        let nx = (cx as i32 + nbr.cell[0]).rem_euclid(l as i32) as usize;
                        let ny = (cy as i32 + nbr.cell[1]).rem_euclid(l as i32) as usize;
                        let nz = (cz as i32 + nbr.cell[2]).rem_euclid(l as i32) as usize;
                        let j = ((nx * l + ny) * l + nz) * n_basis + nbr.basis;
                        out[i].push(j);
                    }
                }
            }
        }
    }
    out
}

/// Classical mean-field Curie temperature `k_B T_C = J z S^2 / 3` for a
/// cubic lattice with coordination `z` and spin magnitude `S`.
pub fn mean_field_tc(kind: LatticeKind, j: f64, s: f64) -> f64 {
    j * kind.coordination() as f64 * s * s / 3.0
}

/// Linear spin-wave frequency of the ferromagnetic Heisenberg model,
/// `omega(k) = 2 J S (z - sum_delta cos(k . delta)) / hbar`, in units with
/// `hbar = 1`.  `neighbors` are the bond vectors of one atom.
pub fn spin_wave_omega(j: f64, s: f64, neighbors: &[V3], k: V3) -> f64 {
    let z = neighbors.len() as f64;
    let mut sum = 0.0;
    for delta in neighbors {
        sum += dot(k, *delta).cos();
    }
    2.0 * j * s * (z - sum)
}

/// Curie-point utilities.
pub struct CuriePoint;

impl CuriePoint {
    pub fn mean_field(kind: LatticeKind, j: f64, s: f64) -> f64 {
        mean_field_tc(kind, j, s)
    }

    /// Published classical sc Heisenberg Monte-Carlo result
    /// (Chen, Ferrenberg & Landau 1993).
    pub const MC_SC_HEISENBERG: f64 = 1.442_99;
}

/// A classical spin chain with unit spins and Landau-Lifshitz-Gilbert
/// dynamics `dS/dt = -gamma S x H_eff + alpha S x dS/dt`, integrated with
/// RK4.  `H_eff = J sum_nn S_j + 2 K (S . e) e + H_ext` (energy units).
pub struct SpinChain {
    pub n: usize,
    pub a: f64,
    pub j: f64,
    pub k: f64,
    pub easy: V3,
    pub h: V3,
    pub gamma: f64,
    pub alpha: f64,
    pub spins: Vec<V3>,
}

impl SpinChain {
    pub fn new(n: usize, a: f64, j: f64, k: f64, gamma: f64, alpha: f64) -> Self {
        Self {
            n,
            a,
            j,
            k,
            easy: [0.0, 0.0, 1.0],
            h: [0.0, 0.0, 0.0],
            gamma,
            alpha,
            spins: vec![[0.0, 0.0, 1.0]; n],
        }
    }

    /// Right-hand side of the explicit (Landau-Lifshitz) form of the LLG.
    fn rhs(&self, spins: &[V3]) -> Vec<V3> {
        let im = |i: usize| (i + self.n - 1) % self.n;
        let ip = |i: usize| (i + 1) % self.n;
        let g = self.gamma / (1.0 + self.alpha * self.alpha);
        let mut out = vec![[0.0; 3]; self.n];
        for i in 0..self.n {
            let s = spins[i];
            let a = dot(s, self.easy);
            let h = [
                self.j * (spins[im(i)][0] + spins[ip(i)][0])
                    + 2.0 * self.k * a * self.easy[0]
                    + self.h[0],
                self.j * (spins[im(i)][1] + spins[ip(i)][1])
                    + 2.0 * self.k * a * self.easy[1]
                    + self.h[1],
                self.j * (spins[im(i)][2] + spins[ip(i)][2])
                    + 2.0 * self.k * a * self.easy[2]
                    + self.h[2],
            ];
            let s_cross_h = cross(s, h);
            let s_cross_scross_h = cross(s, s_cross_h);
            for c in 0..3 {
                out[i][c] = -g * (s_cross_h[c] + self.alpha * s_cross_scross_h[c]);
            }
        }
        out
    }

    pub fn step(&mut self, dt: f64) {
        let k1 = self.rhs(&self.spins);
        let s2: Vec<V3> = (0..self.n)
            .map(|i| add(self.spins[i], scale(k1[i], 0.5 * dt)))
            .collect();
        let k2 = self.rhs(&s2);
        let s3: Vec<V3> = (0..self.n)
            .map(|i| add(self.spins[i], scale(k2[i], 0.5 * dt)))
            .collect();
        let k3 = self.rhs(&s3);
        let s4: Vec<V3> = (0..self.n)
            .map(|i| add(self.spins[i], scale(k3[i], dt)))
            .collect();
        let k4 = self.rhs(&s4);
        for i in 0..self.n {
            let mut s = self.spins[i];
            for c in 0..3 {
                s[c] += dt / 6.0 * (k1[i][c] + 2.0 * k2[i][c] + 2.0 * k3[i][c] + k4[i][c]);
            }
            self.spins[i] = normalize(s);
        }
    }

    /// Seed a single small-amplitude spin-wave mode
    /// `S = (eps cos(k x), eps sin(k x), sqrt(1-eps^2))`.
    pub fn seed_spin_wave(&mut self, k: f64, eps: f64) {
        for i in 0..self.n {
            let x = i as f64 * self.a;
            let z = (1.0 - eps * eps).max(0.0).sqrt();
            self.spins[i] = [eps * (k * x).cos(), eps * (k * x).sin(), z];
        }
    }
}

/// Estimate the angular frequency of a uniformly sampled scalar signal by a
/// grid least-squares fit to `A cos(omega t) + B sin(omega t)`.
pub fn measure_frequency(signal: &[f64], dt: f64, omega_lo: f64, omega_hi: f64, steps: usize) -> f64 {
    let mean = signal.iter().sum::<f64>() / signal.len() as f64;
    let mut best = (omega_lo, 0.0_f64);
    for n in 0..steps {
        let omega = omega_lo + (omega_hi - omega_lo) * n as f64 / (steps - 1) as f64;
        let (mut a, mut b) = (0.0, 0.0);
        for (i, &v) in signal.iter().enumerate() {
            let t = i as f64 * dt;
            a += (v - mean) * (omega * t).cos();
            b += (v - mean) * (omega * t).sin();
        }
        let amp = a * a + b * b;
        if amp > best.1 {
            best = (omega, amp);
        }
    }
    best.0
}

fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
