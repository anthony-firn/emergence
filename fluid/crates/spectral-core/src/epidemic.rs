//! Spatially-resolved virus spread and reproduction (agent-based), coupled to
//! the fluid solver for advection.
//!
//! Virions are Brownian particles (see [`crate::virus`]) that infect
//! susceptible host cells; an infected cell produces new virions and dies.
//! This is the target-cell-limited model of viral dynamics, but spatially
//! resolved:
//!
//!   ∂C/∂t + u·∇C = D∇²C − c·C + p·I        (virion concentration)
//!   ∂S/∂t = −(β/vol) S·V                   (susceptible cells)
//!   ∂I/∂t =  (β/vol) S·V − δ·I             (infected cells)
//!
//! where C = virions/vol, S, I = cell densities, V = total virion count,
//! β = infection rate (vol/(virion·time)), p = virions produced per infected
//! cell per time, c = virion clearance rate, δ = infected-cell death rate.
//! In the well-mixed limit (binding radius ≥ domain) the agent-based model
//! reduces *exactly* to the target-cell ODE, which is used for validation
//! (see tests). More generally, a virion within `r_bind` of a susceptible
//! cell infects it with probability `β·dt / bind_vol`, where `bind_vol` is
//! the volume of the periodic ball of radius `r_bind`. Averaged over
//! uniformly mixed virions the expected per-cell rate is therefore
//! `β·V/vol` for *any* binding radius, exactly matching the mean field.
//!
//! # Emergent-boundary optimization
//!
//! Virions exist only near infected cells, so in a large environment the
//! *active region* — the set of cells containing virions or infected cells —
//! is a small fraction of the domain. This module tracks only that sparse
//! active region and uses a spatial hash so each virion checks only host cells
//! within its binding radius (O(V·local)) instead of the whole domain
//! (O(V·S)). The infection front is the emergent boundary we optimize across:
//! [`Epidemic::active_fraction`] reports how much of the domain is live, and
//! the implied speedup over a full-domain simulation is 1/active_fraction.

use std::collections::{HashMap, HashSet};

use crate::virus::trilinear;

const TWO_PI: f64 = 2.0 * std::f64::consts::PI;
const PI: f64 = std::f64::consts::PI;

/// Host-cell infection state.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CellState {
    Susceptible,
    Infected,
    Dead,
}

/// A host cell at a fixed position.
#[derive(Clone, Copy, Debug)]
pub struct HostCell {
    pub pos: [f64; 3],
    pub state: CellState,
}

/// A virion: position plus a per-particle RNG stream (deterministic).
#[derive(Clone, Copy, Debug)]
pub struct Virion {
    pub pos: [f64; 3],
    rng: u64,
}

/// Spatially-resolved viral epidemic. All rates are in solver units
/// (box [0, 2π)³); see [`crate::virus::Scaling`] to map from SI.
pub struct Epidemic {
    pub cells: Vec<HostCell>,
    pub virions: Vec<Virion>,
    /// Infection rate, vol/(virion·time).
    pub beta: f64,
    /// Virions produced per infected cell per time.
    pub p: f64,
    /// Virion clearance rate, 1/time.
    pub c: f64,
    /// Infected-cell death rate, 1/time.
    pub delta: f64,
    /// Binding radius (virion must be within this to infect a cell).
    pub r_bind: f64,
    /// Virion diffusivity, solver units.
    pub d: f64,
    /// Settling speed (gravity along −z), solver units.
    pub vs: f64,
    /// Domain volume, solver units.
    pub vol: f64,
    /// Volume of the periodic region within `r_bind` of a point. The
    /// per-pair infection probability is normalized by this so that the
    /// spatially resolved kernel has the same mean-field rate as the ODE
    /// for *any* binding radius, not just `r_bind >=` the box.
    bind_vol: f64,
    /// Spatial-hash cell size.
    pub h: f64,
    pub t: f64,
    /// If true, use the naive O(V·S) infection search (for comparison).
    pub naive: bool,
    cell_hash: HashMap<[i64; 3], Vec<usize>>,
    infected: Vec<usize>,
    /// Number of cells in the active region (virions ∪ infected cells).
    pub active_cells: usize,
    /// Total cells partitioning the domain.
    pub total_cells: usize,
    rng: u64,
}

impl Epidemic {
    /// Build an epidemic. `cells` and `virions` are the initial state; the
    /// spatial hash is built from the (fixed) host-cell positions.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        cells: Vec<HostCell>,
        virions: Vec<Virion>,
        beta: f64,
        p: f64,
        c: f64,
        delta: f64,
        r_bind: f64,
        d: f64,
        vs: f64,
        vol: f64,
        h: f64,
        naive: bool,
        seed: u64,
    ) -> Self {
        let mut cell_hash: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
        for (i, cell) in cells.iter().enumerate() {
            cell_hash.entry(cell_key(cell.pos, h)).or_default().push(i);
        }
        let infected = cells
            .iter()
            .enumerate()
            .filter(|(_, c)| c.state == CellState::Infected)
            .map(|(i, _)| i)
            .collect();
        let bind_vol = binding_volume(r_bind);
        let mut e = Self {
            cells,
            virions,
            beta,
            p,
            c,
            delta,
            r_bind,
            d,
            vs,
            vol,
            bind_vol,
            h,
            t: 0.0,
            naive,
            cell_hash,
            infected,
            active_cells: 0,
            total_cells: 0,
            rng: seed,
        };
        e.update_active();
        e
    }

    /// Advance one step of size `dt`. `u` is the physical velocity on the
    /// `n`³ grid (or `None` for a quiescent fluid).
    pub fn step(&mut self, u: Option<&[Vec<f64>; 3]>, n: usize, dt: f64) {
        let kick = (2.0 * self.d * dt).sqrt();
        // A susceptible cell at distance <= r_bind of a virion is infected
        // with probability beta*dt/bind_vol. Averaged over uniformly mixed
        // virions, the expected number of virions within r_bind is
        // V*bind_vol/vol, so the per-cell rate is beta*V/vol: exactly the
        // target-cell ODE mean-field rate for any r_bind.
        let beta_dt = self.beta * dt / self.bind_vol;
        let c_dt = self.c * dt;
        let delta_dt = self.delta * dt;
        let p_dt = self.p * dt;
        let (h, r_bind, vs) = (self.h, self.r_bind, self.vs);
        // 1. Virion transport, clearance, infection.
        let mut infections: Vec<usize> = Vec::new();
        let mut survivors: Vec<Virion> = Vec::with_capacity(self.virions.len());
        let mut cand: Vec<usize> = Vec::new();
        for mut v in std::mem::take(&mut self.virions) {
            let vel = match u {
                Some(u) => trilinear(u, n, v.pos),
                None => [0.0; 3],
            };
            for c in 0..3 {
                let settle = if c == 2 { -vs * dt } else { 0.0 };
                let delta = (vel[c] * dt + settle) + kick * gauss(&mut v.rng);
                v.pos[c] = (v.pos[c] + delta).rem_euclid(TWO_PI);
            }
            if uniform(&mut v.rng) < c_dt {
                continue; // virion cleared
            }
            // Candidates (susceptible, within r_bind) in index order, so the
            // hashed and naive searches consume the RNG identically.
            cand.clear();
            if self.naive {
                cand.extend((0..self.cells.len()).filter(|&ci| {
                    self.cells[ci].state == CellState::Susceptible && dist(v.pos, self.cells[ci].pos) <= r_bind
                }));
            } else {
                for nk in cells_within(v.pos, r_bind, h) {
                    if let Some(idxs) = self.cell_hash.get(&nk) {
                        cand.extend(idxs.iter().copied().filter(|&ci| {
                            self.cells[ci].state == CellState::Susceptible && dist(v.pos, self.cells[ci].pos) <= r_bind
                        }));
                    }
                }
                cand.sort_unstable();
            }
            for &ci in &cand {
                if uniform(&mut v.rng) < beta_dt {
                    infections.push(ci);
                }
            }
            survivors.push(v);
        }
        self.virions = survivors;

        // Apply infections (once per cell per step).
        for &ci in &infections {
            if self.cells[ci].state == CellState::Susceptible {
                self.cells[ci].state = CellState::Infected;
                self.infected.push(ci);
            }
        }

        // 2. Infected cells: produce virions, die. Only iterate the (sparse)
        //    infected list — the emergent-boundary optimization.
        let mut produced: Vec<[f64; 3]> = Vec::new();
        let mut keep: Vec<usize> = Vec::with_capacity(self.infected.len());
        for &ci in &self.infected {
            let cell = &self.cells[ci];
            if cell.state != CellState::Infected {
                continue;
            }
            let nprod = poisson(p_dt, &mut self.rng);
            for _ in 0..nprod {
                produced.push(cell.pos);
            }
            if uniform(&mut self.rng) < delta_dt {
                self.cells[ci].state = CellState::Dead;
            } else {
                keep.push(ci);
            }
        }
        self.infected = keep;
        for pos in produced {
            self.virions.push(Virion { pos, rng: child_seed(&mut self.rng) });
        }

        self.t += dt;
        self.update_active();
    }

    /// (susceptible, infected, dead, virions) counts.
    pub fn stats(&self) -> (usize, usize, usize, usize) {
        let (mut s, mut i, mut d) = (0, 0, 0);
        for c in &self.cells {
            match c.state {
                CellState::Susceptible => s += 1,
                CellState::Infected => i += 1,
                CellState::Dead => d += 1,
            }
        }
        (s, i, d, self.virions.len())
    }

    /// Fraction of the domain that is live (virions ∪ infected cells).
    /// The implied speedup over a full-domain simulation is 1/this.
    pub fn active_fraction(&self) -> f64 {
        self.active_cells as f64 / self.total_cells.max(1) as f64
    }

    fn update_active(&mut self) {
        let mut set: HashSet<[i64; 3]> = HashSet::new();
        for v in &self.virions {
            set.insert(cell_key(v.pos, self.h));
        }
        for &ci in &self.infected {
            set.insert(cell_key(self.cells[ci].pos, self.h));
        }
        self.active_cells = set.len();
        let per = bins_per_axis(self.h);
        self.total_cells = (per * per * per) as usize;
    }
}

/// Build a cubic lattice of susceptible host cells with `n` per axis.
pub fn lattice(n: usize) -> Vec<HostCell> {
    let dx = TWO_PI / n as f64;
    let mut cells = Vec::with_capacity(n * n * n);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                cells.push(HostCell {
                    pos: [(i as f64 + 0.5) * dx, (j as f64 + 0.5) * dx, (k as f64 + 0.5) * dx],
                    state: CellState::Susceptible,
                });
            }
        }
    }
    cells
}

/// `np` virions clustered near the box centre.
pub fn seed_virions(np: usize, seed: u64) -> Vec<Virion> {
    let mut rng = seed;
    let mut out = Vec::with_capacity(np);
    for _ in 0..np {
        let mut pos = [0.0; 3];
        for c in 0..3 {
            pos[c] = std::f64::consts::PI + 0.1 * (uniform(&mut rng) - 0.5);
        }
        out.push(Virion { pos, rng: child_seed(&mut rng) });
    }
    out
}

#[inline]
fn bins_per_axis(h: f64) -> i64 {
    ((TWO_PI / h).ceil() as i64).max(1)
}

#[inline]
fn cell_key(p: [f64; 3], h: f64) -> [i64; 3] {
    let m = bins_per_axis(h);
    [
        ((p[0].rem_euclid(TWO_PI) / h).floor() as i64).rem_euclid(m),
        ((p[1].rem_euclid(TWO_PI) / h).floor() as i64).rem_euclid(m),
        ((p[2].rem_euclid(TWO_PI) / h).floor() as i64).rem_euclid(m),
    ]
}

/// All hash-cell keys whose bins can contain a point within periodic distance
/// `r` of `p`. The keys are wrapped with `rem_euclid`, so the search is correct
/// for any `h` — including `h` that does not divide `2π` — and translating the
/// whole system leaves the candidate set unchanged.
fn cells_within(p: [f64; 3], r: f64, h: f64) -> Vec<[i64; 3]> {
    let m = bins_per_axis(h);
    // Per axis: start key and number of consecutive keys to visit. `k` gets a
    // one-bin margin for the fractional last bin when `h` does not divide 2π.
    let mut ranges = [[0i64; 2]; 3];
    for c in 0..3 {
        let base = (p[c].rem_euclid(TWO_PI) / h).floor() as i64;
        let k = (r / h).ceil() as i64 + 1;
        let count = 2 * k + 1;
        if count >= m {
            ranges[c] = [0, m];
        } else {
            ranges[c] = [(base - k).rem_euclid(m), count];
        }
    }
    let mut out = Vec::with_capacity((ranges[0][1] * ranges[1][1] * ranges[2][1]) as usize);
    for a in 0..ranges[0][1] {
        let i = (ranges[0][0] + a).rem_euclid(m);
        for b in 0..ranges[1][1] {
            let j = (ranges[1][0] + b).rem_euclid(m);
            for c in 0..ranges[2][1] {
                let k = (ranges[2][0] + c).rem_euclid(m);
                out.push([i, j, k]);
            }
        }
    }
    out
}

/// Volume of the periodic region within distance `r` of a point on the torus
/// `[0, 2π)³`. Wrapping each coordinate into `[-π, π)` is a volume-preserving
/// bijection, so this is the volume of a ball of radius `r` intersected with
/// the cube `[-π, π)³`. Used to normalize the infection kernel so a virion at
/// distance `<= r_bind` infects with probability `β·dt / bind_vol`, giving the
/// ODE mean-field rate `β·V/vol` for every binding radius.
fn binding_volume(r: f64) -> f64 {
    let a = PI;
    let full = 8.0 * a * a * a;
    if r <= 0.0 {
        return 0.0;
    }
    if r >= a * 3.0_f64.sqrt() {
        return full;
    }
    let ball = 4.0 / 3.0 * PI * r * r * r;
    if r <= a {
        return ball;
    }
    // One spherical cap per cube face; no overlaps while r <= a√2.
    let caps = 2.0 * PI * (r - a) * (r - a) * (2.0 * r + a);
    if r <= a * 2.0_f64.sqrt() {
        return ball - caps;
    }
    // a√2 < r < a√3: adjacent caps overlap along the 12 cube edges (no triple
    // overlaps yet). Add those pairwise overlaps back. In polar coordinates
    // around the edge, E = (4/3)∫_{asin(a/r)}^{π/4} (r² − a²csc²θ)^{3/2} dθ.
    let lo = (a / r).asin();
    let hi = std::f64::consts::FRAC_PI_4;
    let n = 8192;
    let step = (hi - lo) / n as f64;
    let f = |t: f64| {
        let s = t.sin();
        (r * r - a * a / (s * s)).max(0.0).powf(1.5)
    };
    let mut sum = f(lo) + f(hi);
    for i in 1..n {
        let w = if i % 2 == 0 { 2.0 } else { 4.0 };
        sum += w * f(lo + i as f64 * step);
    }
    let edges = 12.0 * (4.0 / 3.0) * (sum * step / 3.0);
    ball - caps + edges
}

#[inline]
fn wrap(d: f64) -> f64 {
    d - TWO_PI * (d / TWO_PI).round()
}

#[inline]
fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let (dx, dy, dz) = (wrap(a[0] - b[0]), wrap(a[1] - b[1]), wrap(a[2] - b[2]));
    (dx * dx + dy * dy + dz * dz).sqrt()
}

// --- deterministic per-particle RNG (xorshift64, same family as the solver) ---

#[inline]
fn next_u64(s: &mut u64) -> u64 {
    let mut x = *s;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *s = x;
    x
}

/// Seed for an independent child stream. A raw `next_u64` output would start
/// the child one step along the parent's own sequence, so consecutive
/// children would replay one stream shifted by one draw; splitmix64 scrambles
/// it to an unrelated point of the xorshift cycle.
fn child_seed(s: &mut u64) -> u64 {
    let mut z = next_u64(s).wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    (z ^ (z >> 31)).max(1)
}

#[inline]
fn uniform(s: &mut u64) -> f64 {
    1.0 - (next_u64(s) >> 11) as f64 / (1u64 << 53) as f64
}

#[inline]
fn gauss(s: &mut u64) -> f64 {
    let u1 = uniform(s);
    let u2 = uniform(s);
    (-2.0 * u1.ln()).sqrt() * (TWO_PI * u2).cos()
}

/// Poisson sample with mean `lambda` (unbiased: E = lambda).
fn poisson(lambda: f64, s: &mut u64) -> usize {
    let base = lambda.floor();
    let frac = lambda - base;
    base as usize + usize::from(uniform(s) < frac)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// RK4 integration of the well-mixed target-cell ODE.
    fn ode(s0: f64, i0: f64, v0: f64, beta: f64, p: f64, c: f64, delta: f64, vol: f64, t_end: f64, dt: f64) -> (f64, f64, f64) {
        let f = |s: f64, i: f64, v: f64| {
            let ds = -(beta / vol) * s * v;
            let di = (beta / vol) * s * v - delta * i;
            let dv = p * i - c * v;
            (ds, di, dv)
        };
        let (mut s, mut i, mut v) = (s0, i0, v0);
        let mut t = 0.0;
        while t < t_end - 1e-12 {
            let (k1s, k1i, k1v) = f(s, i, v);
            let (k2s, k2i, k2v) = f(s + 0.5 * dt * k1s, i + 0.5 * dt * k1i, v + 0.5 * dt * k1v);
            let (k3s, k3i, k3v) = f(s + 0.5 * dt * k2s, i + 0.5 * dt * k2i, v + 0.5 * dt * k2v);
            let (k4s, k4i, k4v) = f(s + dt * k3s, i + dt * k3i, v + dt * k3v);
            s += dt / 6.0 * (k1s + 2.0 * k2s + 2.0 * k3s + k4s);
            i += dt / 6.0 * (k1i + 2.0 * k2i + 2.0 * k3i + k4i);
            v += dt / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
            t += dt;
        }
        (s, i, v)
    }

    /// Well-mixed limit (r_bind ≥ domain): the agent-based model must match
    /// the target-cell ODE. Statistical noise ~1/√(infections) sets the
    /// tolerance.
    #[test]
    fn well_mixed_matches_ode() {
        let vol = TWO_PI * TWO_PI * TWO_PI;
        let n = 10; // 10³ = 1000 cells
        let cells = lattice(n);
        let s0 = cells.len() as f64;
        let v0 = 200.0;
        let virions = seed_virions(v0 as usize, 0xABCD);
        let beta = 1e-4 * vol; // β/vol = 1e-4
        let (p, c, delta) = (2.0, 0.1, 0.2); // R0 = p·(β/vol)·S0/δ = 1.0
        let r_bind = 20.0; // ≥ domain diagonal 2π√3 ≈ 10.9 → well mixed
        let h = r_bind; // one hash cell covers the whole domain (fast well-mixed)
        let mut e = Epidemic::new(cells, virions, beta, p, c, delta, r_bind, 0.0, 0.0, vol, h, false, 0x1234);
        let (dt, t_end) = (0.02, 10.0);
        let mut last = (0, 0, 0, 0);
        for _ in 0..(t_end / dt) as usize {
            e.step(None, 0, dt);
            last = e.stats();
        }
        let (s, i, v) = (last.0 as f64, last.1 as f64, last.3 as f64);
        let (os, oi, ov) = ode(s0, 0.0, v0, beta, p, c, delta, vol, t_end, dt);
        // Relative tolerance ~15% (statistical noise).
        assert!((s - os).abs() / os < 0.15, "S: agent {s:.0} vs ODE {os:.0}");
        assert!((i - oi).abs() / oi.max(1.0) < 0.20, "I: agent {i:.0} vs ODE {oi:.0}");
        assert!((v - ov).abs() / ov.max(1.0) < 0.20, "V: agent {v:.0} vs ODE {ov:.0}");
    }

    /// No virions ⇒ no infection.
    #[test]
    fn no_virions_no_infection() {
        let vol = TWO_PI * TWO_PI * TWO_PI;
        let cells = lattice(8);
        let e = Epidemic::new(cells, vec![], 1e-4 * vol, 50.0, 0.1, 0.2, 0.5, 0.0, 0.0, vol, 0.5, false, 1);
        let (s, i, d, v) = e.stats();
        assert_eq!(s, 512);
        assert_eq!(i, 0);
        assert_eq!(d, 0);
        assert_eq!(v, 0);
    }

    /// Spatial spread: with a small binding radius the infection advances as a
    /// front and the active region is a small fraction of the domain that
    /// grows over time (the emergent boundary).
    ///
    /// Note: the original setup used `beta = 1e-2·vol`, which after the kernel
    /// is normalized to the mean-field rate gives a basic reproduction number
    /// of ~10⁴ — the whole box saturates long before `t_end`, so the "sparse
    /// front" claim could not hold. The parameters below keep R0 = O(1) so the
    /// front really is a small, growing shell.
    #[test]
    fn spatial_front_and_active_region() {
        let vol = TWO_PI * TWO_PI * TWO_PI;
        let n = 24; // 24³ = 13824 cells
        let cells = lattice(n);
        let virions = seed_virions(300, 0x55);
        let beta = 3e-5 * vol; // R0 = (p/c)·(β/vol)·S0/δ ≈ 2
        let (p, c, delta) = (2.0, 0.5, 0.5);
        let r_bind = 0.35;
        let h = r_bind;
        let mut e = Epidemic::new(cells, virions, beta, p, c, delta, r_bind, 0.1, 0.0, vol, h, false, 0x99);
        let (dt, t_end) = (0.05, 8.0);
        let mut fracs = Vec::new();
        let mut infected = 0;
        for _ in 0..(t_end / dt) as usize {
            e.step(None, 0, dt);
            infected = e.stats().1;
            fracs.push(e.active_fraction());
        }
        // Infection spread to a meaningful number of cells.
        assert!(infected > 50, "infection should spread, got {infected} infected");
        // Active region is a small fraction of the domain (emergent boundary).
        let frac = e.active_fraction();
        assert!(frac < 0.5, "active fraction {frac:.3} should be ≪ 1");
        // It grew from a small seed (front expansion).
        assert!(fracs[fracs.len() / 2] > fracs[0], "active region should grow");
    }

    /// The spatial hash and the naive O(V·S) search give the same infection
    /// outcome (deterministic, same seed) — the optimization is exact.
    /// `h = 0.5` deliberately does not divide `2π` and differs from
    /// `r_bind`, exercising the periodic key wrapping. Parameters are
    /// subcritical so the run stays small/fast.
    #[test]
    fn hashed_matches_naive() {
        let vol = TWO_PI * TWO_PI * TWO_PI;
        let n = 10;
        let cells = lattice(n);
        let beta = 5e-4 * vol;
        let (p, c, delta) = (2.0, 2.0, 0.5);
        let r_bind = 0.7;
        let h = 0.5; // does not divide 2π, and differs from r_bind
        let (dt, t_end) = (0.05, 4.0);
        let run = |naive: bool| {
            let virions = seed_virions(200, 0x77);
            let mut e = Epidemic::new(cells.clone(), virions, beta, p, c, delta, r_bind, 0.01, 0.0, vol, h, naive, 0x31);
            for _ in 0..(t_end / dt) as usize {
                e.step(None, 0, dt);
            }
            e.stats()
        };
        let (hs, hi, hd, hv) = run(false);
        let (ns, ni, nd, nv) = run(true);
        assert_eq!((hs, hi, hd, hv), (ns, ni, nd, nv), "hash and naive must agree");
    }

    /// Virions uniformly distributed over the periodic box.
    fn uniform_virions(np: usize, seed: u64) -> Vec<Virion> {
        let mut s = seed | 1;
        let mut out = Vec::with_capacity(np);
        for _ in 0..np {
            let mut pos = [0.0; 3];
            for c in 0..3 {
                pos[c] = uniform(&mut s) * TWO_PI;
            }
            out.push(Virion { pos, rng: child_seed(&mut s) });
        }
        out
    }

    /// Minimum-image distance: opposite sides of the periodic box are close.
    #[test]
    fn periodic_distance_uses_minimum_image() {
        let eps = 0.13;
        let a = [eps, 2.0, 4.0];
        let b = [TWO_PI - eps, 2.0, 4.0];
        assert!((dist(a, b) - 2.0 * eps).abs() < 1e-12, "dist = {}", dist(a, b));
        let off = [0.7, 3.9, 1.1];
        let ta = [a[0] + off[0], a[1] + off[1], a[2] + off[2]];
        let tb = [b[0] + off[0], b[1] + off[1], b[2] + off[2]];
        assert!((dist(ta, tb) - dist(a, b)).abs() < 1e-12);
        assert_eq!(dist([0.0; 3], [0.0; 3]), 0.0);
    }

    /// `cell_key` must be invariant under a full period shift (periodic
    /// bookkeeping) and produce keys in range.
    #[test]
    fn cell_key_is_periodic() {
        for &h in &[0.37, 0.5, 0.7, 1.0, 2.5] {
            let p = [0.3, 1.7, 5.9];
            let shifted = [p[0] + TWO_PI, p[1] - 3.0 * TWO_PI, p[2] + 7.0 * TWO_PI];
            assert_eq!(cell_key(p, h), cell_key(shifted, h));
            let m = bins_per_axis(h);
            for c in cell_key(p, h) {
                assert!((0..m).contains(&c), "key out of range for h = {h}");
            }
        }
    }

    /// The hash neighbour enumeration must return exactly the cells inside the
    /// periodic binding radius, for hash sizes that do and do not divide 2π.
    #[test]
    fn periodic_hash_search_matches_brute_force() {
        let n = 8;
        let cells = lattice(n);
        let r = 0.9;
        for &h in &[0.37, 0.5, 0.7, 1.0, 1.3, 2.1, 5.0] {
            let mut hash: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
            for (i, c) in cells.iter().enumerate() {
                hash.entry(cell_key(c.pos, h)).or_default().push(i);
            }
            let mut s = 0x5EEDu64;
            for _ in 0..300 {
                let p = [uniform(&mut s) * TWO_PI, uniform(&mut s) * TWO_PI, uniform(&mut s) * TWO_PI];
                let mut got: Vec<usize> = Vec::new();
                for key in cells_within(p, r, h) {
                    if let Some(idxs) = hash.get(&key) {
                        got.extend(idxs.iter().copied().filter(|&ci| dist(p, cells[ci].pos) <= r));
                    }
                }
                got.sort_unstable();
                got.dedup();
                let mut want: Vec<usize> =
                    (0..cells.len()).filter(|&ci| dist(p, cells[ci].pos) <= r).collect();
                want.sort_unstable();
                assert_eq!(got, want, "h = {h}, p = {p:?}");
            }
        }
    }

    /// `binding_volume` must equal the volume of the periodic ball. Check it
    /// against an independent midpoint quadrature of the wrapped cube.
    #[test]
    fn binding_volume_matches_quadrature() {
        let a = PI;
        let n = 120usize;
        let dv = (2.0 * a / n as f64).powi(3);
        for &r in &[0.4, 1.5, 3.0, 4.5, 5.2] {
            let mut count = 0usize;
            for i in 0..n {
                let x = -a + (i as f64 + 0.5) * 2.0 * a / n as f64;
                for j in 0..n {
                    let y = -a + (j as f64 + 0.5) * 2.0 * a / n as f64;
                    for k in 0..n {
                        let z = -a + (k as f64 + 0.5) * 2.0 * a / n as f64;
                        if x * x + y * y + z * z <= r * r {
                            count += 1;
                        }
                    }
                }
            }
            let quad = count as f64 * dv;
            let exact = binding_volume(r);
            assert!((quad - exact).abs() / exact < 0.03, "r = {r}: quadrature {quad} vs {exact}");
        }
        // Exact limits and the small-r sphere limit.
        let vol = TWO_PI * TWO_PI * TWO_PI;
        assert_eq!(binding_volume(0.0), 0.0);
        assert_eq!(binding_volume(PI * 3.0_f64.sqrt() + 1.0), vol);
        let r = 0.5;
        assert!((binding_volume(r) - 4.0 / 3.0 * PI * r.powi(3)).abs() < 1e-12);
    }

    /// The per-cell infection rate must be the ODE mean-field rate
    /// `β·V/vol` for *any* binding radius, not just when `r_bind` covers the
    /// box. Ensemble of one-step runs with uniformly mixed virions.
    #[test]
    fn mean_field_rate_independent_of_binding_radius() {
        let vol = TWO_PI * TWO_PI * TWO_PI;
        let n = 8;
        let s0 = (n * n * n) as f64;
        let beta = 1e-3 * vol; // β/vol = 1e-3
        let vcount = 200usize;
        let dt = 0.05;
        let reps = 300usize;
        for &r_bind in &[0.4, 1.0, 2.5] {
            let mut infections = 0usize;
            for rep in 0..reps {
                let mut e = Epidemic::new(
                    lattice(n),
                    uniform_virions(vcount, 0x1000 + rep as u64),
                    beta,
                    0.0,
                    0.0,
                    0.0,
                    r_bind,
                    0.0,
                    0.0,
                    vol,
                    r_bind.max(0.5),
                    false,
                    0x2000 + rep as u64,
                );
                let before = e.stats().0;
                e.step(None, 0, dt);
                infections += before - e.stats().0;
            }
            let expected = (beta / vol) * s0 * vcount as f64 * dt * reps as f64;
            let rel = (infections as f64 - expected).abs() / expected;
            assert!(rel < 0.15, "r_bind = {r_bind}: {infections} vs expected {expected:.1}");
        }
    }

    /// Translating the entire system must not change the outcome. The offset
    /// is not a multiple of `h`, and `h` does not divide `2π`.
    #[test]
    fn translation_invariance() {
        let vol = TWO_PI * TWO_PI * TWO_PI;
        let n = 10;
        let beta = 5e-4 * vol;
        let (p, c, delta) = (2.0, 2.0, 0.5);
        let r_bind = 0.7;
        let h = 0.5; // does not divide 2π
        let (dt, steps) = (0.05, 80);
        let run = |off: [f64; 3]| {
            let mut cells = lattice(n);
            for cell in &mut cells {
                for k in 0..3 {
                    cell.pos[k] = (cell.pos[k] + off[k]).rem_euclid(TWO_PI);
                }
            }
            let mut virions = seed_virions(150, 0x77);
            for v in &mut virions {
                for k in 0..3 {
                    v.pos[k] = (v.pos[k] + off[k]).rem_euclid(TWO_PI);
                }
            }
            // d = 0 keeps positions fixed, so the comparison is free of
            // accumulated Brownian floating-point drift.
            let mut e = Epidemic::new(cells, virions, beta, p, c, delta, r_bind, 0.0, 0.0, vol, h, false, 0x31);
            for _ in 0..steps {
                e.step(None, 0, dt);
            }
            e.stats()
        };
        let base = run([0.0, 0.0, 0.0]);
        let shifted = run([0.371, 1.913, 0.07]);
        assert!(base.1 > 0, "translation test should be non-trivial, got {base:?}");
        assert_eq!(base, shifted, "translation must not change the outcome");
    }
}
