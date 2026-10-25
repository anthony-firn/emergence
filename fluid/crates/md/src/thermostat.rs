//! Thermostats for the canonical (NVT) ensemble.
//!
//! * [`Thermostat::Bussi`] — the stochastic velocity-rescaling thermostat
//!   of Bussi, Donadio & Parrinello (2007), which draws kinetic energies
//!   from the exact canonical distribution. This is the workhorse used for
//!   the thermodynamic benchmarks.
//! * [`Thermostat::Berendsen`] — weak-coupling rescaling (Berendsen et al.
//!   1984). Cheap, but samples the wrong fluctuations; provided for
//!   equilibration only.
//! * [`Thermostat::NoseHoover`] — a single Nosé–Hoover chain coupled to
//!   the atomic degrees of freedom, integrated symmetrically.

use crate::system::System;

/// How the temperature is controlled.
#[derive(Clone, Debug)]
pub enum Thermostat {
    None,
    Bussi { tau: f64 },
    Berendsen { tau: f64 },
    NoseHoover { tau: f64 },
}

impl Thermostat {
    pub fn bussi(tau: f64) -> Self {
        Thermostat::Bussi { tau }
    }
    pub fn berendsen(tau: f64) -> Self {
        Thermostat::Berendsen { tau }
    }
    pub fn nose_hoover(tau: f64) -> Self {
        Thermostat::NoseHoover { tau }
    }

    /// Advance the system one time step under this thermostat.
    pub fn step(&mut self, sys: &mut System, dt: f64, target: f64) {
        match self {
            Thermostat::None => sys.step_nve(dt),
            Thermostat::Bussi { .. } => {
                sys.step_nve(dt);
                self.bussi_rescale(sys, dt, target);
            }
            Thermostat::Berendsen { .. } => {
                sys.step_nve(dt);
                self.berendsen_rescale(sys, dt, target);
            }
            Thermostat::NoseHoover { tau } => {
                sys.step_nose_hoover(dt, target, *tau);
            }
        }
    }

    /// Exact stochastic velocity rescaling (Bussi–Donadio–Parrinello 2007,
    /// appendix Eq. (A7)–(A10)). `N_f` is the number of degrees of freedom.
    fn bussi_rescale(&self, sys: &mut System, dt: f64, target: f64) {
        let Thermostat::Bussi { tau } = self else {
            return;
        };
        let nf = 3.0 * sys.n() as f64;
        // Bussi et al. quantify the kinetic energy as K = sum_i m v_i^2
        // = 2 * (kinetic energy) and target Kbar = N_f k_B T. The update
        // only uses the ratio Kbar/(N_f K), so writing it with the ordinary
        // kinetic energy means using Kbar = N_f k_B T.
        let k = sys.kinetic_energy();
        if k <= 0.0 {
            return;
        }
        let kbar = nf * target;
        let c = (-dt / tau).exp();
        let r1 = sys.rng.normal();
        // sum of N_f - 1 squared standard normals ~ Gamma((N_f-1)/2, 1)
        let s = sys.rng.gamma(0.5 * (nf - 1.0));
        let alpha2 = c
            + (kbar / (nf * k)) * (1.0 - c) * (r1 * r1 + s)
            + 2.0 * c.sqrt() * ((kbar * (1.0 - c)) / (nf * k)).sqrt() * r1;
        let alpha = alpha2.max(0.0).sqrt();
        for v in &mut sys.vel {
            for a in 0..3 {
                v[a] *= alpha;
            }
        }
    }

    fn berendsen_rescale(&self, sys: &mut System, dt: f64, target: f64) {
        let Thermostat::Berendsen { tau } = self else {
            return;
        };
        let nf = 3.0 * sys.n() as f64;
        let t_inst = 2.0 * sys.kinetic_energy() / nf;
        if t_inst <= 0.0 {
            return;
        }
        let lambda = (1.0 + (dt / tau) * (target / t_inst - 1.0)).max(0.0).sqrt();
        for v in &mut sys.vel {
            for a in 0..3 {
                v[a] *= lambda;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::potential::{Cutoff, Potential};
    use materials::LennardJones;

    #[test]
    fn bussi_reaches_the_target_temperature() {
        let pot = Potential::LennardJones(LennardJones::new(1.0, 1.0));
        let mut sys = System::fcc(3, 0.8, pot, 2.5, 0.3, Cutoff::TailCorrected, 1);
        sys.set_velocities_mb(1.0);
        let mut thermo = Thermostat::bussi(0.5);
        for _ in 0..2000 {
            thermo.step(&mut sys, 0.004, 1.5);
        }
        let mut sum = 0.0;
        let ns = 3000;
        for _ in 0..ns {
            thermo.step(&mut sys, 0.004, 1.5);
            sum += sys.temperature_com();
        }
        let mean = sum / ns as f64;
        assert!(
            (mean - 1.5).abs() < 0.05,
            "Bussi mean temperature {mean}, target 1.5"
        );
    }
}
