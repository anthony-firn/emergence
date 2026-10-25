//! A simple Berendsen barostat for NPT runs.
//!
//! Berendsen, Postma, van Gunsteren, DiNola & Haak, *J. Chem. Phys.* **81**,
//! 3684 (1984). The box is scaled isotropically by
//!
//! ```text
//! mu = [1 - (dt/tau_p) kappa (P0 - P)]^{1/3},
//! ```
//!
//! with `kappa` the isothermal compressibility. The scheme is a weak
//! coupling: it drives the average pressure to `P0` but, like the
//! Berendsen thermostat, does not sample the correct NPT fluctuations. It
//! is provided for equilibration and for reaching a target density; the
//! production NVT benchmarks use the Bussi thermostat at fixed volume.

use crate::system::System;

/// Weak-coupling isotropic barostat.
#[derive(Clone, Copy, Debug)]
pub struct BerendsenBarostat {
    /// Pressure-coupling time.
    pub tau: f64,
    /// Isothermal compressibility (reduced units).
    pub compressibility: f64,
    /// Target pressure.
    pub target_pressure: f64,
}

impl BerendsenBarostat {
    pub fn new(target_pressure: f64, tau: f64, compressibility: f64) -> Self {
        Self {
            tau,
            compressibility,
            target_pressure,
        }
    }

    /// Advance the box by one step at the current computed pressure.
    pub fn step(&self, sys: &mut System, dt: f64) {
        let p = sys.pressure();
        let arg = 1.0 - (dt / self.tau) * self.compressibility * (self.target_pressure - p);
        if arg <= 0.0 {
            return;
        }
        sys.scale_box(arg.powf(1.0 / 3.0));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::potential::{Cutoff, Potential};
    use crate::thermostat::Thermostat;
    use materials::LennardJones;

    #[test]
    fn berendsen_barostat_drives_the_pressure_toward_the_target() {
        // Start a dense liquid at high pressure and couple to P0 = 0.5.
        let pot = Potential::LennardJones(LennardJones::new(1.0, 1.0));
        let mut sys = System::fcc(4, 0.9, pot, 2.5, 0.3, Cutoff::TailCorrected, 1);
        sys.set_velocities_mb(1.0);
        let mut thermo = Thermostat::bussi(0.5);
        for _ in 0..2000 {
            thermo.step(&mut sys, 0.004, 1.0);
        }
        let p0 = 0.5;
        let baro = BerendsenBarostat::new(p0, 1.0, 0.05);
        let p_start = sys.pressure();
        for _ in 0..4000 {
            thermo.step(&mut sys, 0.004, 1.0);
            baro.step(&mut sys, 0.004);
        }
        let p_end = sys.pressure();
        println!("pressure {p_start:.3} -> {p_end:.3}, target {p0}");
        assert!(
            (p_end - p0).abs() < (p_start - p0).abs(),
            "barostat did not reduce the pressure error"
        );
    }
}
