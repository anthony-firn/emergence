//! Locally resonant acoustic metamaterial: the mass-in-mass lattice.
//!
//! A chain of outer masses `m1` joined by springs `k1`, each carrying an inner
//! resonator `m2` on a spring `k2`, is the canonical locally resonant
//! metamaterial. Its unit cell has two degrees of freedom, so there are two
//! branches, and the inner resonance opens a **sub-wavelength** band gap: the
//! gap frequency is set by `k2/m2`, not by the Bragg condition of the lattice.
//!
//! Eliminating the inner mass gives the exact dynamic mass
//!
//! ```text
//!   m_eff(omega) = m1 + m2 omega0^2 / (omega0^2 - omega^2),
//!   omega0^2 = k2 / m2,
//! ```
//!
//! which is **negative** for `omega0 < omega < omega0 sqrt(1 + m2/m1)`. The
//! lattice dispersion is then the monatomic form with `m1 -> m_eff`,
//!
//! ```text
//!   cos(q a) = 1 - omega^2 m_eff(omega) / (2 k1).
//! ```
//!
//! The negative-mass interval sits inside the true lattice gap: the acoustic
//! branch tops out at the frequency where `omega^2 m_eff = 4 k1`, slightly
//! below `omega0`. Both are returned, and the exact 2-DOF unit-cell eigenproblem
//! is used as an independent check of the effective-mass elimination.
//!
//! The cheap coarse model is the homogeneous chain with the static mass
//! `m1 + m2`; the indicator is `omega / omega0` (or the gap width). Unfold to
//! the two-degree-of-freedom cell once the drive approaches `omega0`.
use crate::complex::C;

/// One-dimensional mass-in-mass lattice with lattice constant `a`.
#[derive(Clone, Copy, Debug)]
pub struct MassInMass {
    pub a: f64,
    pub m1: f64,
    pub m2: f64,
    pub k1: f64,
    pub k2: f64,
}

impl MassInMass {
    pub fn new(a: f64, m1: f64, m2: f64, k1: f64, k2: f64) -> Self {
        assert!(a > 0.0 && m1 > 0.0 && m2 > 0.0 && k1 > 0.0 && k2 > 0.0);
        Self { a, m1, m2, k1, k2 }
    }
    /// Inner resonance `omega0 = sqrt(k2/m2)`.
    pub fn omega0(&self) -> f64 {
        (self.k2 / self.m2).sqrt()
    }
    /// Long-wavelength (static) acoustic speed of the chain.
    pub fn acoustic_speed(&self) -> f64 {
        self.a * (self.k1 / (self.m1 + self.m2)).sqrt()
    }
    /// Bragg frequency of the bare chain, `2 sqrt(k1/m1)`.
    pub fn bragg_frequency(&self) -> f64 {
        2.0 * (self.k1 / self.m1).sqrt()
    }
    /// Exact effective dynamic mass.
    pub fn effective_mass(&self, omega: f64) -> f64 {
        let w0sq = self.omega0() * self.omega0();
        self.m1 + self.m2 * w0sq / (w0sq - omega * omega)
    }
    /// Frequency where `m_eff = 0`, the top of the gap and of the
    /// negative-mass interval, `omega0 sqrt(1 + m2/m1)`.
    pub fn gap_top(&self) -> f64 {
        (self.k2 * (1.0 / self.m1 + 1.0 / self.m2)).sqrt()
    }
    /// Interval in which the effective mass is negative,
    /// `(omega0, gap_top)`.
    pub fn negative_mass_interval(&self) -> (f64, f64) {
        (self.omega0(), self.gap_top())
    }
    /// Lower edge of the true lattice gap, where the acoustic branch reaches
    /// the zone boundary: the smaller root of
    /// `m1 x^2 - (m1 omega0^2 + k2 + 4 k1) x + 4 k1 omega0^2 = 0`.
    pub fn gap_bottom(&self) -> f64 {
        let w0sq = self.omega0() * self.omega0();
        let b = self.m1 * w0sq + self.k2 + 4.0 * self.k1;
        let disc = (b * b - 4.0 * self.m1 * 4.0 * self.k1 * w0sq).max(0.0).sqrt();
        let x = (b - disc) / (2.0 * self.m1);
        x.max(0.0).sqrt()
    }
    /// The full lattice band gap `(gap_bottom, gap_top)`. The negative-mass
    /// interval is contained in it.
    pub fn band_gap(&self) -> (f64, f64) {
        (self.gap_bottom(), self.gap_top())
    }
    /// Effective-mass dispersion `q(omega)`, `cos(q a) = 1 - omega^2 m_eff /
    /// (2 k1)`. Real in a pass band, complex in the gap.
    pub fn wavevector(&self, omega: f64) -> C {
        let arg = C::real(1.0 - omega * omega * self.effective_mass(omega) / (2.0 * self.k1));
        if arg.im.abs() < 1e-14 && arg.re.abs() <= 1.0 {
            C::real(arg.re.acos() / self.a)
        } else {
            let mag = arg.abs().max(1.0);
            let re_q = if arg.re < 0.0 {
                std::f64::consts::PI / self.a
            } else {
                0.0
            };
            C::new(re_q, mag.acosh() / self.a)
        }
    }
    /// Exact two-degree-of-freedom unit-cell eigenvalues `(omega_low^2,
    /// omega_high^2)` for a given `cos(q a)`. Solves
    /// `m1 m2 x^2 - (P m2 + m1 k2) x + 2 k1 k2 (1 - cos) = 0` with
    /// `P = 2 k1 (1 - cos) + k2`.
    pub fn eigen_omega_squared(&self, cos_qa: f64) -> (f64, f64) {
        let p = 2.0 * self.k1 * (1.0 - cos_qa) + self.k2;
        let b = p * self.m2 + self.m1 * self.k2;
        let c = 2.0 * self.k1 * self.k2 * (1.0 - cos_qa);
        let disc = (b * b - 4.0 * self.m1 * self.m2 * c).max(0.0).sqrt();
        let lo = (b - disc) / (2.0 * self.m1 * self.m2);
        let hi = (b + disc) / (2.0 * self.m1 * self.m2);
        (lo, hi)
    }
    /// Exact two-degree-of-freedom eigenfrequencies for a given Bloch
    /// wavenumber.
    pub fn eigen_frequencies(&self, q: f64) -> (f64, f64) {
        let (lo, hi) = self.eigen_omega_squared((q * self.a).cos());
        (lo.max(0.0).sqrt(), hi.max(0.0).sqrt())
    }
    /// Indicator that the inner resonance is sub-wavelength: the gap
    /// wavelength divided by the lattice constant, computed at `omega0` with
    /// the static speed.
    pub fn gap_wavelength_over_lattice(&self) -> f64 {
        2.0 * std::f64::consts::PI * self.acoustic_speed() / (self.a * self.omega0())
    }
}
