//! Harmonic lattice dynamics: the dynamical matrix `D(k)` and phonon
//! frequencies for any lattice and pair potential.
//!
//! With force constants `Φ_{ij}(m) = ∂²E/∂u_i(0)∂u_j(m)`, the dynamical
//! matrix is `D_{ij}(k) = (1/√(m_i m_j)) Σ_m Φ_{ij}(m) e^{ik·R}`. For a
//! pair potential each bond contributes the stiffness `M` of
//! [`crate::potential::bond_stiffness`] to the on-site block and
//! `−M e^{ik·R}` to the bond block, which enforces the acoustic sum rule
//! exactly. Eigenvalues of the Hermitian `D` are `ω²`; in the code they are
//! obtained from the equivalent real symmetric embedding.
//!
//! Benchmarks: the 1D monatomic chain reproduces
//! `ω = 2√(K/m)|sin(ka/2)|` to machine precision; `ω(k=0) = 0` (acoustic
//! sum rule); at the equilibrium lattice constant the dispersion is real
//! and non-negative (dynamically stable).

use crate::lattice::{dot, norm, Bond, Lattice, V3};
use crate::linalg::sym_eigenvalues;
use crate::potential::{bond_stiffness, PairPotential};
use std::f64::consts::PI;

/// Precomputed bond list and masses for repeated dispersion evaluations.
pub struct Harmonics {
    pub bonds: Vec<Vec<Bond>>,
    pub masses: Vec<f64>,
    pub natoms: usize,
}

impl Harmonics {
    pub fn new(lat: &Lattice, rmax: f64, masses: Vec<f64>) -> Self {
        assert_eq!(masses.len(), lat.natoms());
        Self {
            bonds: lat.bond_list(rmax),
            masses,
            natoms: lat.natoms(),
        }
    }

    /// Unit masses, one per basis atom.
    pub fn unit_masses(lat: &Lattice, rmax: f64) -> Self {
        Self::new(lat, rmax, vec![1.0; lat.natoms()])
    }

    fn dim(&self) -> usize {
        3 * self.natoms
    }

    /// Hermitian dynamical matrix at wavevector `k`, returned as
    /// `(real, imag)` flat row-major of size `dim × dim`.
    pub fn dynamical(&self, pot: &dyn PairPotential, k: V3) -> Vec<(f64, f64)> {
        let dim = self.dim();
        let mut d = vec![(0.0, 0.0); dim * dim];
        for i in 0..self.natoms {
            for b in &self.bonds[i] {
                let m = bond_stiffness(pot, b.r);
                let phase = dot(k, b.r);
                let (sin, cos) = phase.sin_cos();
                let amp = (self.masses[i] * self.masses[b.j]).sqrt();
                for a in 0..3 {
                    for c in 0..3 {
                        // On-site block.
                        d[(3 * i + a) * dim + (3 * i + c)].0 += m[a][c] / self.masses[i];
                        // Bond block: −M e^{ik·R}/√(m_i m_j).
                        let idx = (3 * i + a) * dim + (3 * b.j + c);
                        d[idx].0 -= m[a][c] * cos / amp;
                        d[idx].1 -= m[a][c] * sin / amp;
                    }
                }
            }
        }
        d
    }

    /// Phonon frequencies `ω = √(eigenvalue)` at `k`, ascending. Negative
    /// numerical eigenvalues are clipped to zero.
    pub fn frequencies(&self, pot: &dyn PairPotential, k: V3) -> Vec<f64> {
        self.eigenvalues(pot, k)
            .into_iter()
            .map(|x| x.max(0.0).sqrt())
            .collect()
    }

    /// Raw eigenvalues `ω²` of the dynamical matrix at `k`, ascending. These
    /// are the physical (non-duplicated) eigenvalues; a dynamically stable
    /// configuration has all of them non-negative.
    pub fn eigenvalues(&self, pot: &dyn PairPotential, k: V3) -> Vec<f64> {
        let dim = self.dim();
        let d = self.dynamical(pot, k);
        // Real symmetric embedding [[B, -C], [C, B]] of D = B + iC. Every
        // Hermitian eigenvalue of D appears exactly twice, so after sorting
        // the physical spectrum is the even-indexed subsequence (the
        // duplicate block of each eigenvalue has even length and starts at
        // an even index). Taking the first `dim` values would instead grab
        // the transverse zeros of a chain.
        let mut r = vec![0.0; (2 * dim) * (2 * dim)];
        for p in 0..dim {
            for q in 0..dim {
                let (re, im) = d[p * dim + q];
                r[p * 2 * dim + q] = re;
                r[p * 2 * dim + (dim + q)] = -im;
                r[(dim + p) * 2 * dim + q] = im;
                r[(dim + p) * 2 * dim + (dim + q)] = re;
            }
        }
        let e = sym_eigenvalues(&r, 2 * dim);
        e.into_iter().step_by(2).take(dim).collect()
    }

    /// Acoustic-branch velocities along direction `n` at small `k = q n`.
    pub fn acoustic_velocities(&self, pot: &dyn PairPotential, n: V3, q: f64) -> Vec<f64> {
        let len = norm(n);
        let khat = [n[0] / len, n[1] / len, n[2] / len];
        let k = [khat[0] * q, khat[1] * q, khat[2] * q];
        let mut w = self.frequencies(pot, k);
        w.truncate(3.min(w.len()));
        w.into_iter().map(|omega| omega / q).collect()
    }
}

/// 1D monatomic chain dispersion (exact): `ω = 2√(K/m)|sin(ka/2)|`,
/// written out here as the closed-form reference for the general code.
pub fn chain_frequency(k: f64, a: f64, spring: f64, mass: f64) -> f64 {
    2.0 * (spring / mass).sqrt() * (k * a / 2.0).sin().abs()
}

/// High-symmetry path for fcc, in units of `2π/a`, as (label, k) points.
/// Γ–X–W–L–Γ. The Brillouin-zone corner coordinates are standard
/// (Ashcroft & Mermin).
pub fn fcc_path() -> Vec<(&'static str, V3)> {
    vec![
        ("Γ", [0.0, 0.0, 0.0]),
        ("X", [0.0, 1.0, 0.0]),
        ("W", [0.5, 1.0, 0.0]),
        ("L", [0.5, 0.5, 0.5]),
        ("Γ", [0.0, 0.0, 0.0]),
    ]
}

/// Sample a piecewise-linear `k`-path with `points_per_segment` points per
/// segment and return `(distance, frequencies)`.
pub fn sample_path(
    harmonics: &Harmonics,
    pot: &dyn PairPotential,
    path: &[(&str, V3)],
    scale: f64,
    points_per_segment: usize,
) -> Vec<(f64, Vec<f64>)> {
    let mut out = Vec::new();
    let mut dist = 0.0;
    for seg in path.windows(2) {
        let a = seg[0].1;
        let b = seg[1].1;
        let seg_len = norm([b[0] - a[0], b[1] - a[1], b[2] - a[2]]) * scale;
        for s in 0..points_per_segment {
            let t = s as f64 / points_per_segment as f64;
            let k = [
                (a[0] + (b[0] - a[0]) * t) * scale,
                (a[1] + (b[1] - a[1]) * t) * scale,
                (a[2] + (b[2] - a[2]) * t) * scale,
            ];
            out.push((dist + t * seg_len, harmonics.frequencies(pot, k)));
        }
        dist += seg_len;
    }
    let last = path.last().unwrap().1;
    out.push((dist, harmonics.frequencies(pot, last)));
    out
}

/// Angular frequency of a 1D chain built as a general lattice, for
/// cross-checking against [`chain_frequency`].
pub fn chain_lattice(a: f64) -> Lattice {
    const HUGE: f64 = 1.0e6;
    Lattice::custom(
        [[a, 0.0, 0.0], [0.0, HUGE, 0.0], [0.0, 0.0, HUGE]],
        vec![[0.0, 0.0, 0.0]],
    )
}

/// The reciprocal-space length `2π/a` used to scale the path.
pub fn two_pi_over(a: f64) -> f64 {
    2.0 * PI / a
}
