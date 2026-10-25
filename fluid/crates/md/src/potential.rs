//! Pair potentials and the truncation schemes the MD engine integrates.
//!
//! Three cutoff conventions are supported, all with the same force
//! expression where it applies:
//!
//! * [`Cutoff::Shifted`] — truncated and shifted, `u_ts(r) = u(r) - u(rc)`
//!   for `r <= rc`. Continuous in energy, discontinuous in force, so it is
//!   the right scheme for NVE energy-conservation tests.
//! * [`Cutoff::ShiftedForce`] — also subtracts the force jump,
//!   `u(r) - u(rc) - (r-rc) u'(rc)`. Energy and force both continuous.
//! * [`Cutoff::TailCorrected`] — bare truncation of `u`, with the standard
//!   homogeneous long-range corrections to energy and pressure added
//!   analytically. This is the scheme to compare with the infinite-range
//!   Lennard-Jones equation of state.
//!
//! The Lennard-Jones tail formulas are the textbook ones (Frenkel & Smit,
//! *Understanding Molecular Simulation*): energy per atom
//! `8 pi rho eps sigma^3 [(1/9) x^9 - (1/3) x^3]` and pressure correction
//! `-(16/3) pi rho^2 eps sigma^3 [x^3 - (2/3) x^9]`, `x = sigma/rc`.

use materials::{LennardJones, PairPotential};

/// How the potential is cut off at `rc`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cutoff {
    Shifted,
    ShiftedForce,
    TailCorrected,
}

/// A pair potential usable by the MD engine.
#[derive(Clone, Copy, Debug)]
pub enum Potential {
    LennardJones(LennardJones),
    Morse(materials::Morse),
}

impl Potential {
    /// Bare pair energy `u(r)`.
    #[inline]
    pub fn energy(&self, r: f64) -> f64 {
        match self {
            Potential::LennardJones(p) => p.energy(r),
            Potential::Morse(p) => p.energy(r),
        }
    }

    /// `du/dr`.
    #[inline]
    pub fn d1(&self, r: f64) -> f64 {
        match self {
            Potential::LennardJones(p) => p.d1(r),
            Potential::Morse(p) => p.d1(r),
        }
    }

    /// Homogeneous tail energy per atom for a bare truncation at `rc`.
    pub fn tail_energy_per_atom(&self, rho: f64, rc: f64) -> f64 {
        match self {
            Potential::LennardJones(p) => p.tail_energy_per_atom(rho, rc),
            Potential::Morse(_) => 0.0,
        }
    }

    /// Homogeneous tail pressure correction for a bare truncation at `rc`
    /// (negative: the missing attraction lowers the pressure).
    pub fn tail_pressure(&self, rho: f64, rc: f64) -> f64 {
        match self {
            Potential::LennardJones(p) => {
                let x = p.sigma / rc;
                let x3 = x * x * x;
                let x9 = x3 * x3 * x3;
                -(16.0 / 3.0)
                    * std::f64::consts::PI
                    * rho
                    * rho
                    * p.epsilon
                    * p.sigma.powi(3)
                    * (x3 - (2.0 / 3.0) * x9)
            }
            Potential::Morse(_) => 0.0,
        }
    }

    /// Sigma of a Lennard-Jones potential, or `1.0` otherwise.
    pub fn sigma(&self) -> f64 {
        match self {
            Potential::LennardJones(p) => p.sigma,
            Potential::Morse(_) => 1.0,
        }
    }

    /// Epsilon of a Lennard-Jones potential, or `1.0` otherwise.
    pub fn epsilon(&self) -> f64 {
        match self {
            Potential::LennardJones(p) => p.epsilon,
            Potential::Morse(_) => 1.0,
        }
    }
}

impl From<LennardJones> for Potential {
    fn from(p: LennardJones) -> Self {
        Potential::LennardJones(p)
    }
}

/// Scheme-aware pair energy.
#[inline]
pub fn pair_energy(pot: &Potential, r: f64, rc: f64, scheme: Cutoff) -> f64 {
    match scheme {
        Cutoff::Shifted => pot.energy(r) - pot.energy(rc),
        Cutoff::ShiftedForce => pot.energy(r) - pot.energy(rc) - (r - rc) * pot.d1(rc),
        Cutoff::TailCorrected => pot.energy(r),
    }
}

/// Scheme-aware `du_eff/dr`. The force on atom `i` from `j` is
/// `du_eff/dr * (r_j - r_i)/r`, so this is also the radial force magnitude.
#[inline]
pub fn pair_force(pot: &Potential, r: f64, rc: f64, scheme: Cutoff) -> f64 {
    match scheme {
        Cutoff::Shifted | Cutoff::TailCorrected => pot.d1(r),
        Cutoff::ShiftedForce => pot.d1(r) - pot.d1(rc),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shifted_potential_is_continuous_in_energy() {
        let pot = Potential::LennardJones(LennardJones::new(1.0, 1.0));
        let rc = 2.5;
        let below = pair_energy(&pot, rc - 1e-9, rc, Cutoff::Shifted);
        assert!(below.abs() < 1e-8);
        // Shifted-force has a continuous derivative at rc.
        let eps = 1e-6;
        let f_below = pair_force(&pot, rc - eps, rc, Cutoff::ShiftedForce);
        let f_above = 0.0;
        assert!((f_below - f_above).abs() < 1e-5);
    }

    #[test]
    fn tail_pressure_has_the_expected_sign() {
        let pot = Potential::LennardJones(LennardJones::new(1.0, 1.0));
        // Missing attraction must lower the pressure.
        assert!(pot.tail_pressure(0.8, 2.5) < 0.0);
        // And the energy tail is negative too.
        assert!(pot.tail_energy_per_atom(0.8, 2.5) < 0.0);
    }
}
