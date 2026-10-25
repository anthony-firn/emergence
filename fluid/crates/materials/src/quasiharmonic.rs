//! The quasi-harmonic rung: thermal properties of crystals from lattice
//! dynamics.
//!
//! The harmonic phonons of [`crate::phonon`] give the vibrational free
//! energy as a Brillouin-zone average over the discrete spectrum
//! `ω_{qν}`:
//!
//! ```text
//! F_vib(V,T) = (1/N_q) Σ_q Σ_ν [ ħω/2 + k_B T ln(1 − e^{−ħω/k_B T}) ]
//! ```
//!
//! per primitive cell (the classical limit replaces the bracket by
//! `k_B T ln(ħω/k_B T)`). Summing `F(V,T) = E_static(V) + F_vib(V,T)`
//! over volumes and minimising at zero pressure is the quasi-harmonic
//! approximation (QHA): the only volume dependence of the vibrations is
//! through the harmonic frequencies, so thermal expansion enters through
//! the mode Grüneisen parameters
//!
//! ```text
//! γ_ν(q) = −∂ln ω_{qν}/∂ln V.
//! ```
//!
//! From the minimiser `V(T)` this module returns the thermal expansion
//! `α`, the isothermal bulk modulus `B_T(T)`, the heat capacity `C_V(T)`,
//! the entropy, `C_P − C_V = α² B_T V T`, and the thermodynamic Grüneisen
//! parameter (both from `α B_T V / C_V` and as a mode heat-capacity
//! weighted average of `γ_ν`).
//!
//! The pattern of [`crate::lib`] continues: the harmonic crystal is cheap
//! but only valid while the vibrations stay small compared with the
//! nearest-neighbour distance (Lindemann ratio ≲ 0.1) and while the
//! temperature is not so high that the neglected cubic/anharmonic terms
//! matter. This module exposes the two indicators of that boundary — the
//! quantum/classical crossover `T/θ_D`, and the Lindemann ratio
//! `√⟨u²⟩/d_nn` — so a finer (molecular-dynamics) rung can take over where
//! they fail.
//!
//! The quantum path uses `ħω/2` (zero-point energy); at `T → 0` the QHA
//! volume is therefore *larger* than the static minimum. That zero-point
//! lattice expansion is a genuine quantum effect reproduced below for
//! solid argon.

use crate::continuum::christoffel_velocities;
use crate::lattice::{cross, dot, scale, Lattice, V3};
use crate::phonon::Harmonics;
use crate::potential::{energy_per_atom_with_tail, PairPotential};
use std::f64::consts::PI;

/// Boltzmann constant, SI (2019 redefinition), J/K.
pub const KB: f64 = 1.380_649e-23;
/// Reduced Planck constant, SI, J·s.
pub const HBAR: f64 = 1.054_571_817e-34;
/// Atomic mass unit, SI, kg (CODATA 2018).
pub const AMU: f64 = 1.660_539_066_60e-27;

/// Which statistics to use for the vibrational free energy.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Statistics {
    /// `ħω/2 + k_B T ln(1 − e^{−ħω/k_B T})`; the physical, quantum rung.
    Quantum,
    /// `k_B T ln(ħω/k_B T)`; the classical fold (no zero-point energy,
    /// `C_V = 3 k_B` per atom by construction).
    Classical,
}

/// Reciprocal lattice vectors `b_i` with `a_i · b_j = 2π δ_ij`.
pub fn reciprocal_vectors(lat: &Lattice) -> [V3; 3] {
    let v = dot(lat.prim[0], cross(lat.prim[1], lat.prim[2]));
    let two_pi = 2.0 * PI;
    [
        scale(cross(lat.prim[1], lat.prim[2]), two_pi / v),
        scale(cross(lat.prim[2], lat.prim[0]), two_pi / v),
        scale(cross(lat.prim[0], lat.prim[1]), two_pi / v),
    ]
}

/// Monkhorst–Pack grid on a `n = (n1, n2, n3)` mesh.
///
/// Fractional coordinates `u_i = (2m − n_i − 1)/(2 n_i)` for
/// `m = 0..n_i−1`, mapped to Cartesian `k = Σ u_i b_i`. The points are
/// equally weighted (`1/Πn_i`); symmetry is not used to reduce the grid.
/// For even `n_i` no point falls on Γ, which avoids the three exact
/// zero acoustic modes; for odd `n_i` Γ is present and those modes are
/// skipped by the thermodynamic sums (they carry zero phase-space
/// measure).
pub fn monkhorst_pack(lat: &Lattice, n: [usize; 3]) -> Vec<V3> {
    let b = reciprocal_vectors(lat);
    let mut out = Vec::with_capacity(n[0] * n[1] * n[2]);
    for i in 0..n[0] {
        let u0 = (2.0 * i as f64 - n[0] as f64 - 1.0) / (2.0 * n[0] as f64);
        for j in 0..n[1] {
            let u1 = (2.0 * j as f64 - n[1] as f64 - 1.0) / (2.0 * n[1] as f64);
            for k in 0..n[2] {
                let u2 = (2.0 * k as f64 - n[2] as f64 - 1.0) / (2.0 * n[2] as f64);
                let mut kv = [0.0; 3];
                for c in 0..3 {
                    kv[c] = u0 * b[0][c] + u1 * b[1][c] + u2 * b[2][c];
                }
                out.push(kv);
            }
        }
    }
    out
}

/// A Brillouin-zone Monkhorst–Pack grid and the per-point weight.
#[derive(Clone, Debug)]
pub struct KGrid {
    pub points: Vec<V3>,
    pub n: [usize; 3],
    pub weight: f64,
}

impl KGrid {
    pub fn monkhorst_pack(lat: &Lattice, n: [usize; 3]) -> Self {
        let npts = n[0] * n[1] * n[2];
        Self {
            points: monkhorst_pack(lat, n),
            n,
            weight: 1.0 / npts as f64,
        }
    }
}

/// Eigenvalues (ascending) of a 3×3 Hermitian matrix. The argument is
/// `[H00, H11, H22, Re H01, Im H01, Re H02, Im H02, Re H12, Im H12]`.
/// Used instead of the generic Jacobi solver because the Brillouin-zone
/// loops call it hundreds of thousands of times.
pub(crate) fn hermitian3_eigenvalues(m: [f64; 9]) -> [f64; 3] {
    let (a, d, f) = (m[0], m[1], m[2]);
    let (br, bi) = (m[3], m[4]);
    let (cr, ci) = (m[5], m[6]);
    let (er, ei) = (m[7], m[8]);
    let tr = a + d + f;
    let b2 = br * br + bi * bi;
    let c2 = cr * cr + ci * ci;
    let e2 = er * er + ei * ei;
    let minors = (a * d - b2) + (a * f - c2) + (d * f - e2);
    // det(H) = a(df−|e|²) − b(b* f − e* c) + c(b* e* − d c*), real part.
    let x = br * f - (er * cr + ei * ci);
    let y = -bi * f - (er * ci - ei * cr);
    let term2 = br * x - bi * y;
    let z = br * er - bi * ei - d * cr;
    let w = -br * ei - bi * er + d * ci;
    let term3 = cr * z - ci * w;
    let det = a * (d * f - e2) - term2 + term3;
    let p = minors - tr * tr / 3.0;
    let q = -2.0 * tr * tr * tr / 27.0 + tr * minors / 3.0 - det;
    if p > -1e-12 * (1.0 + tr * tr) {
        // Isotropic (all eigenvalues equal): p = 0 for a Hermitian matrix.
        return [tr / 3.0; 3];
    }
    let r = 2.0 * (-p / 3.0).sqrt();
    let arg = (3.0 * q / (p * r)).clamp(-1.0, 1.0);
    let theta = arg.acos() / 3.0;
    let mut out = [0.0; 3];
    for (k, o) in out.iter_mut().enumerate() {
        *o = r * (theta - 2.0 * PI * k as f64 / 3.0).cos() + tr / 3.0;
    }
    out.sort_by(|x, y| x.total_cmp(y));
    out
}

/// Raw eigenvalues `ω²` of the dynamical matrix at `k`, using the fast
/// 3×3 path for a monatomic (Bravais) cell and the generic solver
/// otherwise. Negative values signal a dynamical instability.
pub fn eigenvalues(harmonics: &Harmonics, pot: &dyn PairPotential, k: V3) -> Vec<f64> {
    if harmonics.natoms == 1 {
        let d = harmonics.dynamical(pot, k);
        let m = [
            d[0].0, d[4].0, d[8].0, d[1].0, d[1].1, d[2].0, d[2].1, d[5].0, d[5].1,
        ];
        hermitian3_eigenvalues(m).to_vec()
    } else {
        harmonics.eigenvalues(pot, k)
    }
}

/// Phonon frequencies at `k`, using the fast 3×3 path for a monatomic
/// (Bravais) cell and the generic solver otherwise.
pub fn spectrum(harmonics: &Harmonics, pot: &dyn PairPotential, k: V3) -> Vec<f64> {
    eigenvalues(harmonics, pot, k)
        .into_iter()
        .map(|x| x.max(0.0).sqrt())
        .collect()
}

/// Phonon density of states `g(ω)` from a Brillouin-zone grid,
/// normalised so that `Σ g_i Δω = 3` per atom (the three translational
/// modes at Γ are dropped). Returns `(bin centres, g)`.
pub fn density_of_states(
    harmonics: &Harmonics,
    pot: &dyn PairPotential,
    kpoints: &[V3],
    n_bins: usize,
) -> (Vec<f64>, Vec<f64>) {
    let mut all: Vec<f64> = Vec::new();
    for k in kpoints {
        for w in spectrum(harmonics, pot, *k) {
            if w > 0.0 {
                all.push(w);
            }
        }
    }
    let wmax = all.iter().cloned().fold(0.0, f64::max);
    let dw = wmax / n_bins as f64;
    let mut g = vec![0.0; n_bins];
    let norm = 1.0 / (kpoints.len() as f64 * harmonics.natoms as f64 * dw);
    for w in all {
        let mut b = (w / dw) as usize;
        if b >= n_bins {
            b = n_bins - 1;
        }
        g[b] += norm;
    }
    let centres = (0..n_bins).map(|i| (i as f64 + 0.5) * dw).collect();
    (centres, g)
}

/// Debye temperature from the elastic tensor and the Christoffel sound
/// speeds:
///
/// ```text
/// θ_D = (ħ v_D / k_B) (6π² n)^{1/3},
/// 3/v_D³ = (1/4π) ∫ dΩ (1/v_l³ + 1/v_t1³ + 1/v_t2³),
/// ```
///
/// with `n` the number density of atoms. The angular integral is a
/// midpoint quadrature over `n_theta × n_phi` directions. `mass_density`
/// is the mass per unit volume (used by the Christoffel velocities).
pub fn debye_temperature_elastic(
    c: &[[f64; 6]; 6],
    mass_density: f64,
    number_density: f64,
    hbar: f64,
    kb: f64,
    n_theta: usize,
    n_phi: usize,
) -> f64 {
    let dtheta = PI / n_theta as f64;
    let dphi = 2.0 * PI / n_phi as f64;
    let mut acc = 0.0;
    for it in 0..n_theta {
        let theta = (it as f64 + 0.5) * dtheta;
        let (st, ct) = theta.sin_cos();
        for ip in 0..n_phi {
            let phi = (ip as f64 + 0.5) * dphi;
            let (sp, cp) = phi.sin_cos();
            let n = [st * cp, st * sp, ct];
            let v = christoffel_velocities(c, mass_density, n);
            acc += st
                * (1.0 / v[0].powi(3) + 1.0 / v[1].powi(3) + 1.0 / v[2].powi(3));
        }
    }
    // (1/4π) ∫ dΩ f = (π/(2 n_theta n_phi)) Σ f sinθ.
    let avg = PI / (2.0 * n_theta as f64 * n_phi as f64) * acc;
    let v_debye = (3.0 / avg).powf(1.0 / 3.0);
    hbar / kb * v_debye * (6.0 * PI * PI * number_density).powf(1.0 / 3.0)
}

/// Debye temperature implied by a low-temperature heat capacity that obeys
/// `C_V = (12π⁴/5) k_B (T/θ_D)³` per atom. `a3 = C_V/T³` is in
/// `energy/(atom·K⁴)`; with `k_B = 1` and dimensionless `C_V`, pass the
/// dimensionless value.
pub fn debye_temperature_from_cv3(a3: f64, kb: f64) -> f64 {
    let c = 12.0 * PI.powi(4) / 5.0 * kb;
    (c / a3).powf(1.0 / 3.0)
}

/// One Brillouin-zone-average thermodynamic evaluation at volume scale
/// `s` (all lengths are `s ×` the base lattice) and temperature `t`.
#[derive(Clone, Copy, Debug)]
pub struct Thermo {
    /// Vibrational free energy per atom.
    pub f: f64,
    /// Vibrational internal energy per atom.
    pub u: f64,
    /// Heat capacity per atom at constant volume.
    pub cv: f64,
    /// Smallest `ω²` over the grid; negative means dynamically unstable.
    pub min_eigenvalue: f64,
    /// Largest `ω²` over the grid.
    pub max_eigenvalue: f64,
    /// `true` when no mode is soft (`min_eigenvalue > 1e-4 max_eigenvalue`),
    /// the domain where the quasi-harmonic free energy is meaningful.
    pub stable: bool,
}

/// Quasi-harmonic model: a base lattice, a potential, a Brillouin-zone
/// grid and the physical constants. `rmax` is the bond cutoff at scale 1
/// (in the same length units as `base`); it is scaled with the lattice so
/// the bond set — and therefore `F(V,T)` — is smooth in the volume.
pub struct Qha<'a> {
    pub base: Lattice,
    pub pot: &'a dyn PairPotential,
    pub rmax: f64,
    pub masses: Vec<f64>,
    pub grid: KGrid,
    pub hbar: f64,
    pub kb: f64,
    /// Uniform lattice scale that minimises the static energy.
    pub static_scale: f64,
    /// Largest scale with no soft phonon mode (the harmonic spinodal).
    /// Temperature independent, since the harmonic frequencies do not
    /// depend on `T`; beyond it the quasi-harmonic free energy diverges.
    pub spinodal_scale: f64,
}

impl<'a> Qha<'a> {
    /// Build the model. `grid_size` is the Monkhorst–Pack mesh. The static
    /// minimum is found once by [`crate::potential::relax`] (tail-corrected
    /// energy), and the harmonic spinodal by a scan over scales.
    pub fn new(
        base: Lattice,
        pot: &'a dyn PairPotential,
        rmax: f64,
        masses: Vec<f64>,
        grid_size: [usize; 3],
        hbar: f64,
        kb: f64,
    ) -> Self {
        assert_eq!(masses.len(), base.natoms());
        let grid = KGrid::monkhorst_pack(&base, grid_size);
        let static_scale = crate::potential::relax(&base, pot, rmax, 1.3, 1.8).scale;
        let mut q = Self {
            base,
            pot,
            rmax,
            masses,
            grid,
            hbar,
            kb,
            static_scale,
            spinodal_scale: static_scale,
        };
        q.spinodal_scale = q.find_spinodal();
        q
    }

    /// Scan upward from the static minimum for the first scale where a
    /// phonon mode goes soft (`min ω² ≤ 1e-4 max ω²`). A dedicated coarse
    /// 4×4×4 grid is used so the scan is cheap and independent of the
    /// thermodynamic mesh.
    fn find_spinodal(&self) -> f64 {
        let s0 = self.static_scale;
        let cap = 1.25 * s0;
        let steps = 80;
        let coarse = monkhorst_pack(&self.base, [4, 4, 4]);
        let mut last = s0;
        for i in 1..=steps {
            let s = s0 + (cap - s0) * i as f64 / steps as f64;
            let h = self.harmonics(s);
            let mut min_eig = f64::INFINITY;
            let mut max_eig: f64 = 0.0;
            for k in &coarse {
                for e in eigenvalues(&h, self.pot, *k) {
                    min_eig = min_eig.min(e);
                    max_eig = max_eig.max(e);
                }
            }
            if min_eig <= 1e-4 * max_eig {
                break;
            }
            last = s;
        }
        last
    }

    /// Bond list and masses at a uniform volume scale `s`.
    pub fn harmonics(&self, s: f64) -> Harmonics {
        let lat = self.base.scaled(s);
        Harmonics::new(&lat, self.rmax * s, self.masses.clone())
    }

    /// Static (0 K, no zero-point) energy per atom at scale `s`.
    pub fn static_energy_per_atom(&self, s: f64) -> f64 {
        let lat = self.base.scaled(s);
        energy_per_atom_with_tail(&lat, self.pot, self.rmax * s)
    }

    /// Brillouin-zone-averaged vibrational free energy, internal energy and
    /// heat capacity per atom. One pass over the spectrum. Unstable
    /// configurations (any `ω² < 0`) are reported through `stable`.
    pub fn thermo(&self, s: f64, t: f64, stats: Statistics) -> Thermo {
        let h = self.harmonics(s);
        let nq = self.grid.points.len() as f64;
        let na = self.base.natoms() as f64;
        let mut f = 0.0;
        let mut u = 0.0;
        let mut cv = 0.0;
        let mut min_eig = f64::INFINITY;
        let mut max_eig: f64 = 0.0;
        for k in &self.grid.points {
            for e in eigenvalues(&h, self.pot, *k) {
                min_eig = min_eig.min(e);
                max_eig = max_eig.max(e);
                if !(e > 0.0) {
                    continue; // three acoustic zeros at Γ, or an unstable mode
                }
                let w = e.sqrt();
                match stats {
                    Statistics::Classical => {
                        let ln = (self.hbar * w / (self.kb * t)).ln();
                        f += self.kb * t * ln;
                        u += self.kb * t;
                        cv += self.kb;
                    }
                    Statistics::Quantum => {
                        let x = self.hbar * w / (self.kb * t);
                        let ex = (-x).exp();
                        let om = -(-x).exp_m1(); // 1 − e^{−x}, stable
                        f += 0.5 * self.hbar * w + self.kb * t * om.ln();
                        u += 0.5 * self.hbar * w + self.hbar * w * ex / om;
                        cv += self.kb * x * x * ex / (om * om);
                    }
                }
            }
        }
        Thermo {
            f: f / (nq * na),
            u: u / (nq * na),
            cv: cv / (nq * na),
            min_eigenvalue: min_eig,
            max_eigenvalue: max_eig,
            stable: min_eig > 1e-4 * max_eig,
        }
    }

    /// Total free energy per atom `E_static(s) + F_vib(s,T)`, `+∞` outside
    /// the harmonic stability domain (where the free energy is not defined).
    pub fn total_free_energy(&self, s: f64, t: f64, stats: Statistics) -> f64 {
        let th = self.thermo(s, t, stats);
        if !th.stable {
            return f64::INFINITY;
        }
        self.static_energy_per_atom(s) + th.f
    }

    /// Entropy per atom, `S = (U − F)/T`.
    pub fn entropy_per_atom(&self, s: f64, t: f64, stats: Statistics) -> f64 {
        let th = self.thermo(s, t, stats);
        (th.u - th.f) / t
    }

    /// Minimise the total free energy over the uniform scale at zero
    /// pressure, restricted to the harmonic stability domain. A coarse scan
    /// plus golden-section refinement. When the minimum sits on the
    /// stability boundary the quasi-harmonic model is breaking down (see
    /// [`ThermalProperties::near_instability`]).
    pub fn equilibrium_scale(&self, t: f64, stats: Statistics) -> f64 {
        let s0 = self.static_scale;
        let lo = 0.90 * s0;
        let hi = (1.15 * s0).min(self.spinodal_scale).max(s0);
        let f = |s: f64| self.total_free_energy(s, t, stats);
        let n = 40;
        let mut best: usize = 0;
        let mut bestf = f64::INFINITY;
        for i in 0..=n {
            let s = lo + (hi - lo) * i as f64 / n as f64;
            let v = f(s);
            if v < bestf {
                bestf = v;
                best = i;
            }
        }
        let step = (hi - lo) / n as f64;
        let a = (lo + step * (best.saturating_sub(1)) as f64).max(lo);
        let b = (lo + step * (best + 1).min(n) as f64).min(hi);
        golden_min(&f, a, b, 40)
    }

    /// Largest stable scale (the harmonic spinodal).
    pub fn stable_scale_limit(&self, _t: f64, _stats: Statistics) -> f64 {
        self.spinodal_scale
    }

    /// Atom volume at scale `s`.
    pub fn atom_volume(&self, s: f64) -> f64 {
        self.base.volume() / self.base.natoms() as f64 * s * s * s
    }

    /// Scale for an atom volume `v` (cubic scaling).
    pub fn scale_for_volume(&self, v: f64) -> f64 {
        let v1 = self.base.volume() / self.base.natoms() as f64;
        (v / v1).cbrt()
    }

    /// Total free energy as a function of atom volume.
    pub fn total_free_energy_volume(&self, v: f64, t: f64, stats: Statistics) -> f64 {
        self.total_free_energy(self.scale_for_volume(v), t, stats)
    }

    /// Linear thermal expansion `α = (1/V)(∂V/∂T)_P` at zero pressure,
    /// from a central difference of the equilibrium volume.
    pub fn thermal_expansion(&self, t: f64, stats: Statistics) -> f64 {
        let dt = (0.01 * t).max(1e-6);
        let sp = self.equilibrium_scale(t + dt, stats);
        let sm = self.equilibrium_scale((t - dt).max(1e-6), stats);
        let s = self.equilibrium_scale(t, stats);
        // V ∝ s³, so (1/V)dV/dT = 3 (1/s) ds/dT.
        3.0 / s * (sp - sm) / (2.0 * dt)
    }

    /// Isothermal bulk modulus `B_T = V ∂²F/∂V²` at the equilibrium volume.
    pub fn bulk_modulus_isothermal(&self, t: f64, stats: Statistics) -> f64 {
        let s = self.equilibrium_scale(t, stats);
        self.bulk_modulus_at(s, t, stats)
    }

    /// `B_T = V ∂²F/∂V²` at an explicit volume scale. Falls back to a
    /// one-sided difference if a stencil point is outside the stability
    /// domain.
    pub fn bulk_modulus_at(&self, s: f64, t: f64, stats: Statistics) -> f64 {
        let v0 = self.atom_volume(s);
        let dv = 1e-3 * v0;
        let f = |v: f64| self.total_free_energy_volume(v, t, stats);
        let f0 = f(v0);
        let fp = f(v0 + dv);
        let fm = f(v0 - dv);
        if fp.is_finite() && fm.is_finite() {
            v0 * (fp - 2.0 * f0 + fm) / (dv * dv)
        } else if fp.is_finite() {
            let f2 = f(v0 + 2.0 * dv);
            v0 * (f2 - 2.0 * fp + f0) / (dv * dv)
        } else if fm.is_finite() {
            let f2 = f(v0 - 2.0 * dv);
            v0 * (f0 - 2.0 * fm + f2) / (dv * dv)
        } else {
            f64::NAN
        }
    }

    /// Heat capacity at constant pressure minus that at constant volume,
    /// per atom: `C_P − C_V = α² B_T V T` (standard thermodynamic identity).
    pub fn cp_minus_cv(&self, t: f64, stats: Statistics) -> f64 {
        let s = self.equilibrium_scale(t, stats);
        let alpha = self.thermal_expansion(t, stats);
        let b = self.bulk_modulus_isothermal(t, stats);
        let v = self.atom_volume(s);
        alpha * alpha * b * v * t
    }

    /// Thermodynamic Grüneisen parameter `γ = α B_T V / C_V` at zero
    /// pressure.
    pub fn gruneisen_alpha(&self, t: f64, stats: Statistics) -> f64 {
        let s = self.equilibrium_scale(t, stats);
        let cv = self.thermo(s, t, stats).cv;
        let alpha = self.thermal_expansion(t, stats);
        let b = self.bulk_modulus_isothermal(t, stats);
        alpha * b * self.atom_volume(s) / cv
    }

    /// Mode Grüneisen parameters `γ_ν(q) = −∂ln ω/∂ln V` by central
    /// difference in volume, averaged with the mode heat capacities.
    pub fn gruneisen_modes(&self, t: f64, stats: Statistics) -> f64 {
        self.gruneisen_modes_at(self.equilibrium_scale(t, stats), t, stats)
    }

    /// [`gruneisen_modes`](Self::gruneisen_modes) at an explicit volume
    /// scale (useful for the indicator at a fixed configuration).
    pub fn gruneisen_modes_at(&self, s: f64, t: f64, stats: Statistics) -> f64 {
        let eps = 1e-4;
        let s2 = s * (1.0 + eps);
        let s1 = s * (1.0 - eps);
        let h2 = self.harmonics(s2);
        let h1 = self.harmonics(s1);
        let dlnv = 3.0 * (s2 / s1).ln();
        let mut num = 0.0;
        let mut den = 0.0;
        for k in &self.grid.points {
            let w2 = spectrum(&h2, self.pot, *k);
            let w1 = spectrum(&h1, self.pot, *k);
            for (a, b) in w2.iter().zip(w1.iter()) {
                if !(*a > 0.0) || !(*b > 0.0) {
                    continue;
                }
                let gamma = -(a.ln() - b.ln()) / dlnv;
                let c = mode_heat_capacity(*b, t, self.hbar, self.kb, stats);
                num += gamma * c;
                den += c;
            }
        }
        num / den
    }

    /// Mean-square displacement per atom from the harmonic spectrum,
    /// `⟨u²⟩ = (1/N_q m) Σ_{qν} (ħ/2ω) coth(ħω/2k_BT)`. Valid for a
    /// monatomic Bravais lattice (fcc here); for a multi-atom basis the
    /// mass is averaged and this is an estimate.
    pub fn mean_square_displacement(&self, s: f64, t: f64) -> f64 {
        let h = self.harmonics(s);
        let m_avg = self.masses.iter().sum::<f64>() / self.masses.len() as f64;
        let nq = self.grid.points.len() as f64;
        let mut acc = 0.0;
        for k in &self.grid.points {
            for w in spectrum(&h, self.pot, *k) {
                if !(w > 0.0) {
                    continue;
                }
                let coth = if self.hbar * w / (2.0 * self.kb * t) > 30.0 {
                    1.0
                } else {
                    1.0 / (self.hbar * w / (2.0 * self.kb * t)).tanh()
                };
                acc += self.hbar / (2.0 * m_avg * w) * coth;
            }
        }
        acc / nq
    }

    /// Lindemann ratio `√⟨u²⟩ / d_nn` at the equilibrium volume; the
    /// classic indicator of the harmonic/QHA boundary (≈0.1 at melting).
    pub fn lindemann_ratio(&self, t: f64) -> f64 {
        let s = self.equilibrium_scale(t, Statistics::Quantum);
        let rnn = self.base.scaled(s).nn_distance();
        self.mean_square_displacement(s, t).sqrt() / rnn
    }

    /// All zero-pressure thermal properties at temperature `t` in one go
    /// (three free-energy minimisations plus local derivatives).
    pub fn properties(&self, t: f64, stats: Statistics) -> ThermalProperties {
        let s = self.equilibrium_scale(t, stats);
        let dt = (0.01 * t).max(1e-6);
        let sp = self.equilibrium_scale(t + dt, stats);
        let sm = self.equilibrium_scale((t - dt).max(1e-6), stats);
        let alpha = 3.0 / s * (sp - sm) / (2.0 * dt);
        let th = self.thermo(s, t, stats);
        let b = self.bulk_modulus_at(s, t, stats);
        let v = self.atom_volume(s);
        let x = alpha * b * v;
        let rnn = self.base.scaled(s).nn_distance();
        let lindemann = self.mean_square_displacement(s, t).sqrt() / rnn;
        let spinodal = self.stable_scale_limit(t, stats);
        ThermalProperties {
            scale: s,
            volume: v,
            free_energy: self.static_energy_per_atom(s) + th.f,
            cv: th.cv,
            entropy: (th.u - th.f) / t,
            alpha,
            bulk_modulus: b,
            cp_minus_cv: x * alpha * t,
            gruneisen: if th.cv > 0.0 { x / th.cv } else { 0.0 },
            gruneisen_modes: self.gruneisen_modes_at(s, t, stats),
            lindemann,
            spinodal_scale: spinodal,
            near_instability: s > 0.995 * spinodal,
        }
    }
}

/// Zero-pressure thermal properties at one temperature.
#[derive(Clone, Copy, Debug)]
pub struct ThermalProperties {
    /// Equilibrium uniform lattice scale (lengths are `s ×` base).
    pub scale: f64,
    /// Atom volume at that scale.
    pub volume: f64,
    /// Total free energy per atom (static + vibrational).
    pub free_energy: f64,
    /// Heat capacity at constant volume, per atom.
    pub cv: f64,
    /// Entropy per atom.
    pub entropy: f64,
    /// Linear thermal expansion `(1/V)(∂V/∂T)_P`.
    pub alpha: f64,
    /// Isothermal bulk modulus.
    pub bulk_modulus: f64,
    /// `C_P − C_V` per atom.
    pub cp_minus_cv: f64,
    /// Thermodynamic Grüneisen parameter `α B_T V / C_V`.
    pub gruneisen: f64,
    /// Mode heat-capacity weighted average of `−∂lnω/∂lnV`.
    pub gruneisen_modes: f64,
    /// Lindemann ratio `√⟨u²⟩/d_nn`.
    pub lindemann: f64,
    /// Largest stable scale at this temperature (the harmonic spinodal).
    pub spinodal_scale: f64,
    /// `true` when the equilibrium sits on the spinodal boundary, i.e. the
    /// quasi-harmonic description is breaking down.
    pub near_instability: bool,
}

/// Heat capacity per mode (per atom for one branch) at frequency `w` and
/// temperature `t`. Classical statistics gives `k_B`.
pub fn mode_heat_capacity(w: f64, t: f64, hbar: f64, kb: f64, stats: Statistics) -> f64 {
    match stats {
        Statistics::Classical => kb,
        Statistics::Quantum => {
            let x = hbar * w / (kb * t);
            let ex = (-x).exp();
            let om = -(-x).exp_m1();
            kb * x * x * ex / (om * om)
        }
    }
}

/// Quantum/classical crossover indicator: `T/θ_D`. The classical fold is
/// accurate when this is large, the Debye `T³` law when it is small.
pub fn crossover_indicator(t: f64, theta_debye: f64) -> f64 {
    t / theta_debye
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

/// Convenience model for fcc Lennard-Jones argon in *reduced* units
/// (σ = ε = m = 1) with the physical constants baked in as conversion
/// factors below. Working in reduced units keeps every length O(1), which
/// the shared bond-list builder assumes; physical quantities are recovered
/// with the constants [`ARGON_SIGMA_M`], [`ARGON_EPSILON_KB`],
/// [`ARGON_MASS_AMU`] and the reduced [`argon_hbar_natural`].
pub struct ArgonModel {
    /// Dimensionless LJ (σ = ε = 1).
    pub pot: crate::potential::LennardJones,
    pub rmax_scale: f64,
    pub grid: [usize; 3],
}

/// Argon LJ σ in metres (0.3405 nm; Rowley & Nicholson 1972).
pub const ARGON_SIGMA_M: f64 = 0.3405e-9;
/// Argon LJ ε/k_B in kelvin (119.8 K; Rowley & Nicholson 1972).
pub const ARGON_EPSILON_KB: f64 = 119.8;
/// Argon atomic mass in u (39.948 u; standard).
pub const ARGON_MASS_AMU: f64 = 39.948;

/// Frequency scale `√(ε/(m σ²))` of the argon LJ model in rad/s, from the
/// constants above.
pub fn argon_omega_scale() -> f64 {
    let eps = ARGON_EPSILON_KB * KB;
    let m = ARGON_MASS_AMU * AMU;
    (eps / (m * ARGON_SIGMA_M * ARGON_SIGMA_M)).sqrt()
}

/// Reduced Planck constant in the LJ reduced units `(σ, ε, m)`, i.e.
/// `ħ_nat = ħ_phys √(ε/(mσ²)) / ε`. With this value and `k_B = 1` the
/// reduced temperature is in units of `ε/k_B = 119.8 K` for argon.
pub fn argon_hbar_natural() -> f64 {
    HBAR * argon_omega_scale() / (ARGON_EPSILON_KB * KB)
}

impl ArgonModel {
    pub fn new(grid: [usize; 3], rmax_scale: f64) -> Self {
        Self {
            pot: crate::potential::LennardJones::new(1.0, 1.0),
            rmax_scale,
            grid,
        }
    }

    /// A [`Qha`] in reduced units. Natural temperature `T*` maps to
    /// `T = T* · 119.8 K` and natural length to `s · σ`. The reduced `ħ` is
    /// [`argon_hbar_natural`], so energies are in units of `ε`.
    pub fn qha(&self) -> Qha<'_> {
        let base = Lattice::fcc(1.0);
        Qha::new(
            base,
            &self.pot,
            self.rmax_scale,
            vec![1.0],
            self.grid,
            argon_hbar_natural(),
            1.0,
        )
    }
}

/// Debye temperature of the reduced argon model in kelvin from a reduced
/// (natural) value in units of `ε/k_B`: `θ_K = θ* · 119.8 K`.
pub fn argon_theta_kelvin(theta_natural: f64) -> f64 {
    theta_natural * ARGON_EPSILON_KB
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::potential::LennardJones;

    #[test]
    fn fast_hermitian3_matches_generic() {
        let pot = LennardJones::new(1.0, 1.0);
        let s = 1.545;
        let h = Harmonics::unit_masses(&Lattice::fcc(s), 2.5 * s);
        for i in 0..8 {
            for j in 0..8 {
                let k = [0.05 * i as f64, 0.05 * j as f64, 0.05];
                let fast = spectrum(&h, &pot, k);
                let slow = h.frequencies(&pot, k);
                for (a, b) in fast.iter().zip(slow.iter()) {
                    assert!(
                        (a - b).abs() < 1e-6,
                        "k={k:?} fast={fast:?} slow={slow:?}"
                    );
                }
            }
        }
    }
}
