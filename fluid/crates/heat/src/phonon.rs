//! Phonon dispersion models per branch, and the mode heat capacity.

use std::f64::consts::PI;

/// Reduced Planck constant (J s), CODATA 2022.
pub const HBAR: f64 = 1.054_571_817e-34;
/// Boltzmann constant (J/K), SI 2019.
pub const KB: f64 = 1.380_649e-23;

/// Dispersion relation of a branch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Dispersion {
    /// Linear `ω = v q` with a Debye cutoff.
    Debye,
    /// Born–von Kármán sine: `ω = ω_max sin(π q / (2 q_max))`, so that the
    /// group velocity at `q → 0` is the sound speed `v`.
    Sine,
}

/// One acoustic phonon branch (one polarisation).
#[derive(Clone, Copy, Debug)]
pub struct Branch {
    pub name: &'static str,
    /// Long-wavelength sound speed (m/s).
    pub v: f64,
    /// Cutoff angular frequency `ω_max` (rad/s).
    pub omega_max: f64,
    /// Number of equivalent branches (e.g. 2 for the two transverse modes).
    pub degeneracy: f64,
    /// Grüneisen anharmonicity parameter `γ` for the umklapp rate.
    pub gamma: f64,
    pub dispersion: Dispersion,
}

impl Branch {
    /// Debye branch from a cutoff temperature `Θ = ħ ω_max / k_B`.
    pub fn debye(name: &'static str, v: f64, theta: f64, gamma: f64, degeneracy: f64) -> Self {
        Branch {
            name,
            v,
            omega_max: KB * theta / HBAR,
            degeneracy,
            gamma,
            dispersion: Dispersion::Debye,
        }
    }

    /// Sine (Born–von Kármán) branch from a cutoff temperature.
    pub fn sine(name: &'static str, v: f64, theta: f64, gamma: f64, degeneracy: f64) -> Self {
        Branch {
            name,
            v,
            omega_max: KB * theta / HBAR,
            degeneracy,
            gamma,
            dispersion: Dispersion::Sine,
        }
    }

    /// Debye wavevector of the branch. For the sine dispersion it is chosen so
    /// that `dω/dq|_0 = v`.
    pub fn q_max(&self) -> f64 {
        match self.dispersion {
            Dispersion::Debye => self.omega_max / self.v,
            Dispersion::Sine => PI * self.omega_max / (2.0 * self.v),
        }
    }

    /// Branch frequency at wavevector `q` (rad/s).
    pub fn omega(&self, q: f64) -> f64 {
        match self.dispersion {
            Dispersion::Debye => self.v * q,
            Dispersion::Sine => {
                self.omega_max * (PI * q / (2.0 * self.q_max())).sin()
            }
        }
    }

    /// Group speed `|dω/dq|` at wavevector `q` (m/s).
    pub fn group_velocity(&self, q: f64) -> f64 {
        match self.dispersion {
            Dispersion::Debye => self.v,
            Dispersion::Sine => {
                self.omega_max * (PI / (2.0 * self.q_max()))
                    * (PI * q / (2.0 * self.q_max())).cos()
            }
        }
    }

    /// Debye cutoff temperature `Θ = ħ ω_max / k_B` (K).
    pub fn theta(&self) -> f64 {
        HBAR * self.omega_max / KB
    }
}

/// Heat capacity of one mode at `x = ħω/(k_B T)`:
/// `C = k_B x² e^x/(e^x-1)²` (J/K per mode). Classical limit `k_B`.
pub fn mode_heat_capacity(x: f64) -> f64 {
    if x < 1e-8 {
        return KB;
    }
    let ex = x.exp();
    KB * x * x * ex / (ex - 1.0).powi(2)
}

/// Total volumetric heat capacity (J/(m³ K)) from the model branches,
/// `C_V = ∫ C(ω) g(ω) dω`.
pub fn volumetric_heat_capacity(branches: &[Branch], t: f64) -> f64 {
    let mut total = 0.0;
    for b in branches {
        let qmax = b.q_max();
        let n = 200;
        let (xg, wg) = crate::linalg::gauss_legendre_01(n);
        let mut s = 0.0;
        for i in 0..n {
            let q = xg[i] * qmax;
            let w = b.omega(q);
            let x = HBAR * w / (KB * t);
            let g = q * q / (2.0 * PI * PI); // DOS in q-space
            s += wg[i] * mode_heat_capacity(x) * g;
        }
        total += b.degeneracy * qmax * s;
    }
    total
}
