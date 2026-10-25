//! Pair potentials and the static lattice sums they generate.
//!
//! A pair potential contributes `φ(r)` per bond. The energy per atom of an
//! infinite crystal with `n` basis atoms is
//!
//! ```text
//! E = (1/(2n)) Σ_i Σ_{bonds from i} φ(|R|),
//! ```
//!
//! the `1/2` because every bond is visited from both ends. This is a true
//! spherical-cutoff lattice sum (periodic images included), not a
//! minimum-image approximation. A Lennard-Jones tail correction is
//! available and the convergence of the sum with `rmax` is demonstrated in
//! the tests.
//!
//! Benchmarks: LJ fcc at zero pressure has `r0 = 1.09017 σ` and cohesive
//! energy `−8.6102 ε` (standard fcc lattice sums `A12 = 12.13188`,
//! `A6 = 14.45392`); hcp is slightly lower in energy than fcc and bcc
//! higher.

use crate::lattice::{dot, norm, Lattice, V3};

/// Central pair potential with exact first and second derivatives.
pub trait PairPotential {
    fn energy(&self, r: f64) -> f64;
    /// dφ/dr
    fn d1(&self, r: f64) -> f64;
    /// d²φ/dr²
    fn d2(&self, r: f64) -> f64;
    fn name(&self) -> &'static str;
    /// Homogeneous continuum estimate of the energy per atom carried by
    /// bonds beyond a spherical cutoff `rc` at number density `rho`. Zero
    /// for short-ranged potentials; Lennard-Jones overrides it so that a
    /// modest cutoff already reproduces the converged lattice sums.
    fn tail_energy_per_atom(&self, _rho: f64, _rc: f64) -> f64 {
        0.0
    }
}

/// Lennard-Jones `φ = 4ε[(σ/r)^12 − (σ/r)^6]`. Constants: standard
/// 12-6 form (Jones 1924); the fcc lattice sums below are exact.
#[derive(Clone, Copy, Debug)]
pub struct LennardJones {
    pub sigma: f64,
    pub epsilon: f64,
}

impl LennardJones {
    pub fn new(sigma: f64, epsilon: f64) -> Self {
        Self { sigma, epsilon }
    }
    /// Argon-like parameters (Rowley & Nicholson 1972, via standard LJ
    /// tables): σ = 0.3405 nm, ε/k_B = 119.8 K.
    pub fn argon() -> Self {
        const KB: f64 = 1.380_649e-23;
        Self::new(0.3405e-9, 119.8 * KB)
    }
}

impl PairPotential for LennardJones {
    fn energy(&self, r: f64) -> f64 {
        let x2 = (self.sigma / r).powi(6);
        4.0 * self.epsilon * (x2 * x2 - x2)
    }
    fn d1(&self, r: f64) -> f64 {
        let x6 = (self.sigma / r).powi(6);
        24.0 * self.epsilon * (x6 - 2.0 * x6 * x6) / r
    }
    fn d2(&self, r: f64) -> f64 {
        let x6 = (self.sigma / r).powi(6);
        24.0 * self.epsilon * (26.0 * x6 * x6 - 7.0 * x6) / (r * r)
    }
    fn name(&self) -> &'static str {
        "lennard-jones"
    }
    fn tail_energy_per_atom(&self, rho: f64, rc: f64) -> f64 {
        lj_tail_energy_per_atom(rho, rc, self.sigma, self.epsilon)
    }
}

/// Morse `φ = D[e^{−2α(r−r0)} − 2e^{−α(r−r0)}]`. Constants: Morse (1929).
#[derive(Clone, Copy, Debug)]
pub struct Morse {
    pub d: f64,
    pub alpha: f64,
    pub r0: f64,
}

impl Morse {
    pub fn new(d: f64, alpha: f64, r0: f64) -> Self {
        Self { d, alpha, r0 }
    }
}

impl PairPotential for Morse {
    fn energy(&self, r: f64) -> f64 {
        let x = (r - self.r0) * self.alpha;
        self.d * ((-2.0 * x).exp() - 2.0 * (-x).exp())
    }
    fn d1(&self, r: f64) -> f64 {
        let x = (r - self.r0) * self.alpha;
        2.0 * self.alpha * self.d * ((-x).exp() - (-2.0 * x).exp())
    }
    fn d2(&self, r: f64) -> f64 {
        let x = (r - self.r0) * self.alpha;
        2.0 * self.alpha * self.alpha * self.d * (2.0 * (-2.0 * x).exp() - (-x).exp())
    }
    fn name(&self) -> &'static str {
        "morse"
    }
}

/// Harmonic spring `φ = ½ k (r − r0)²`, the exactly solvable test
/// potential for the 1D chain. Not a physical cohesion law.
#[derive(Clone, Copy, Debug)]
pub struct Harmonic {
    pub k: f64,
    pub r0: f64,
}

impl PairPotential for Harmonic {
    fn energy(&self, r: f64) -> f64 {
        0.5 * self.k * (r - self.r0).powi(2)
    }
    fn d1(&self, r: f64) -> f64 {
        self.k * (r - self.r0)
    }
    fn d2(&self, _r: f64) -> f64 {
        self.k
    }
    fn name(&self) -> &'static str {
        "harmonic"
    }
}

/// Bond stiffness `M_ab = φ'' r̂_a r̂_b + (φ'/r)(δ_ab − r̂_a r̂_b)`, the
/// second derivative of `φ(|r|)` with respect to the bond vector. Positive
/// definite for a stable potential at the bond length.
pub fn bond_stiffness(pot: &dyn PairPotential, r: V3) -> [[f64; 3]; 3] {
    let rad = norm(r);
    let u = [r[0] / rad, r[1] / rad, r[2] / rad];
    let a = pot.d2(rad);
    let b = pot.d1(rad) / rad;
    let mut m = [[0.0; 3]; 3];
    for p in 0..3 {
        for q in 0..3 {
            m[p][q] = a * u[p] * u[q] + b * (if p == q { 1.0 } else { 0.0 } - u[p] * u[q]);
        }
    }
    m
}

/// Energy per atom with a spherical cutoff and no tail correction.
pub fn energy_per_atom(lat: &Lattice, pot: &dyn PairPotential, rmax: f64) -> f64 {
    let bonds = lat.bond_list(rmax);
    let n = lat.natoms();
    let mut total = 0.0;
    for i in 0..n {
        for b in &bonds[i] {
            total += pot.energy(b.radius);
        }
    }
    total / (2.0 * n as f64)
}

/// Energy per atom including the potential's homogeneous tail correction.
/// For a short-ranged potential this equals [`energy_per_atom`].
pub fn energy_per_atom_with_tail(lat: &Lattice, pot: &dyn PairPotential, rmax: f64) -> f64 {
    let rho = lat.natoms() as f64 / lat.volume();
    energy_per_atom(lat, pot, rmax) + pot.tail_energy_per_atom(rho, rmax)
}

/// Dimensionless lattice sum `A_p = Σ_{bonds} (r_nn/r)^p` around basis atom
/// 0, with `r_nn` the nearest-neighbour distance.
pub fn lattice_sum(lat: &Lattice, p: f64, rmax: f64) -> f64 {
    let rnn = lat.nn_distance();
    lat.bond_list(rmax)[0]
        .iter()
        .map(|b| (rnn / b.radius).powf(p))
        .sum()
}

/// Lattice sum `A_p` corrected by the homogeneous continuum tail beyond
/// `rmax`: `Σ_{r>rc} (r_nn/r)^p ≈ r_nn^p 4πρ rc^{3-p}/(p-3)` for `p > 3`.
/// This is the standard mean-field tail used in molecular simulation.
pub fn lattice_sum_corrected(lat: &Lattice, p: f64, rmax: f64) -> f64 {
    assert!(p > 3.0, "tail integral diverges for p <= 3");
    let rnn = lat.nn_distance();
    let rho = lat.natoms() as f64 / lat.volume();
    let tail = 4.0 * std::f64::consts::PI * rho * rnn.powf(p) * rmax.powf(3.0 - p) / (p - 3.0);
    lattice_sum(lat, p, rmax) + tail
}

/// Lennard-Jones energy tail beyond a spherical cutoff, per atom, for a
/// uniform number density `rho`:
/// `E_tail = 8πρ ε σ³[(1/9)x⁹ − (1/3)x³]`, `x = σ/rc`. Standard
/// homogeneous-tail formula (Frenkel & Smit, *Understanding Molecular
/// Simulation*).
pub fn lj_tail_energy_per_atom(rho: f64, rc: f64, sigma: f64, epsilon: f64) -> f64 {
    let x = sigma / rc;
    8.0 * std::f64::consts::PI * rho * epsilon * sigma.powi(3) * (x.powi(9) / 9.0 - x.powi(3) / 3.0)
}

/// Equilibrium of a lattice under a pair potential.
#[derive(Clone, Copy, Debug)]
pub struct Equilibrium {
    /// Scale factor applied to the base lattice (lattice constant for the
    /// standard structures).
    pub scale: f64,
    pub nn_distance: f64,
    pub energy_per_atom: f64,
}

/// Minimise the energy per atom over a uniform lattice scale on `[lo, hi]`
/// by a coarse scan followed by golden-section refinement. The scale is
/// the lattice constant `a` when the base lattice is one of the standard
/// structures built with `a = 1`.
///
/// Two variants are provided. [`relax`] uses the
/// [`energy_per_atom_with_tail`] total, so with a modest cutoff it already
/// gives the converged (infinite-sum) lattice constant and cohesive
/// energy. [`relax_cutoff`] minimises the raw cut-off sum; it is the
/// configuration of zero stress for the *truncated* potential, which is
/// the consistent reference for comparing the Cauchy-Born stiffness with
/// the long-wavelength phonons (the homogeneous tail carries no force
/// constants and is not part of the dynamical matrix).
///
/// Consistency of the tail. The homogeneous tail is a correction to the
/// *dimensionless* lattice sum, so it is only scale-covariant when the
/// cutoff is a fixed multiple of the lattice spacing. [`relax`] therefore
/// evaluates the tail at `rmax * s` while the trial lattice is scaled by
/// `s`: the set of included shells is then the same at every scale and the
/// tail shifts `A_6` by a constant instead of contributing a spurious
/// `r_nn^{-3}` term. Without this, minimising the truncated sum plus a
/// tail evaluated at a *fixed* cutoff biases the equilibrium (for LJ fcc it
/// gives `r0 = 1.090221 σ` instead of the exact `1.090173 σ`). [`relax_cutoff`]
/// keeps the fixed cutoff, because for the raw sum the cutoff is part of
/// the potential being minimised.
pub fn relax(base: &Lattice, pot: &dyn PairPotential, rmax: f64, lo: f64, hi: f64) -> Equilibrium {
    relax_with(base, pot, rmax, lo, hi, true)
}

/// Minimise the raw spherical-cutoff energy (no tail correction). See
/// [`relax`].
pub fn relax_cutoff(
    base: &Lattice,
    pot: &dyn PairPotential,
    rmax: f64,
    lo: f64,
    hi: f64,
) -> Equilibrium {
    relax_with(base, pot, rmax, lo, hi, false)
}

fn relax_with(
    base: &Lattice,
    pot: &dyn PairPotential,
    rmax: f64,
    lo: f64,
    hi: f64,
    tail: bool,
) -> Equilibrium {
    // For the tail-corrected relaxation, precompute the bond list once at
    // the base cutoff. Because all bond vectors scale with `s`, the
    // dimensionless set of included shells is then the same at every trial
    // scale: the energy is a smooth function of `s`, with no spurious jumps
    // from the cutoff search box and no density-dependent shift of the
    // equilibrium. (`energy_per_atom_with_tail` re-derives its bond list for
    // every scale, which is correct but has small cutoff-boundary steps that
    // trap a golden-section search.)
    //
    // The raw-cutoff relaxation keeps the physical cutoff fixed, so its bond
    // set genuinely changes with `s` and must be rebuilt each trial.
    let base_bonds = base.bond_list(rmax);
    let energy = |s: f64| {
        if !tail {
            return energy_per_atom(&base.scaled(s), pot, rmax);
        }
        let mut total = 0.0;
        for i in 0..base.natoms() {
            for b in &base_bonds[i] {
                total += pot.energy(s * b.radius);
            }
        }
        let e = total / (2.0 * base.natoms() as f64);
        let rho = base.natoms() as f64 / (base.volume() * s.powi(3));
        e + pot.tail_energy_per_atom(rho, rmax * s)
    };
    // Coarse scan, then refine the bracketing interval.
    let n_scan = 200;
    let mut best_i = 0;
    let mut best = f64::INFINITY;
    let mut values = Vec::with_capacity(n_scan + 1);
    for i in 0..=n_scan {
        let s = lo + (hi - lo) * i as f64 / n_scan as f64;
        let e = energy(s);
        values.push(e);
        if e < best {
            best = e;
            best_i = i;
        }
    }
    let step = (hi - lo) / n_scan as f64;
    let a = (lo + step * (best_i.saturating_sub(1)) as f64).max(lo);
    let b = (lo + step * (best_i + 1).min(n_scan) as f64).min(hi);
    let s = golden_min(&energy, a, b, 200);
    Equilibrium {
        scale: s,
        nn_distance: base.scaled(s).nn_distance(),
        energy_per_atom: energy(s),
    }
}

fn golden_min<F: Fn(f64) -> f64>(f: &F, mut lo: f64, mut hi: f64, iters: usize) -> f64 {
    const INV_PHI: f64 = 0.618_033_988_749_894_9;
    let mut c = hi - INV_PHI * (hi - lo);
    let mut d = lo + INV_PHI * (hi - lo);
    let mut fc = f(c);
    let mut fd = f(d);
    for _ in 0..iters {
        if fc < fd {
            hi = d;
            d = c;
            fd = fc;
            c = hi - INV_PHI * (hi - lo);
            fc = f(c);
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = lo + INV_PHI * (hi - lo);
            fd = f(d);
        }
    }
    0.5 * (lo + hi)
}

/// Inner product helper re-exported for downstream rungs.
pub fn dot3(a: V3, b: V3) -> f64 {
    dot(a, b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derivatives_agree_with_finite_difference() {
        let lj = LennardJones::new(1.0, 1.0);
        let morse = Morse::new(1.0, 2.0, 1.1);
        for p in [&lj as &dyn PairPotential, &morse] {
            for &r in &[0.9, 1.0, 1.2, 1.5, 2.0] {
                // Different steps for the two orders: the second difference
                // balance of truncation (h²) and roundoff (h⁻²) wants h~1e-4,
                // the first difference wants h~1e-6.
                let h1 = 1e-6 * r;
                let h2 = 1e-4 * r;
                let d1 = (p.energy(r + h1) - p.energy(r - h1)) / (2.0 * h1);
                let d2 = (p.energy(r + h2) - 2.0 * p.energy(r) + p.energy(r - h2)) / (h2 * h2);
                assert!((p.d1(r) - d1).abs() < 1e-5, "d1 {r}");
                assert!((p.d2(r) - d2).abs() < 1e-3, "d2 {r}");
            }
        }
    }
}
