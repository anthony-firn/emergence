//! Material dispersion: how the complex refractive index depends on wavelength.
//!
//! Wave optics needs the material response as a function of frequency, and the
//! response is causal. Causality ties the real and imaginary parts through the
//! Kramers-Kronig relations, which is why the model functions here are chosen
//! from a small set of analytic, causal forms rather than interpolated:
//!
//! * **Sellmeier** for transparent dielectrics far below their absorption
//!   edge: `n^2(lambda) = 1 + sum_i B_i lambda^2/(lambda^2 - C_i)`.
//! * **Drude** for the free-electron response of a metal:
//!   `eps(omega) = 1 - omega_p^2/(omega^2 + i gamma omega)`, so
//!   `eps -> 1 - omega_p^2/omega^2` at high frequency and `n ~ 0`,
//!   `kappa = omega_p/omega` at low frequency.
//! * **Lorentz** oscillators for bound electrons:
//!   `eps(omega) = 1 + sum_j f_j omega_p^2/(omega_j^2 - omega^2 - i gamma_j omega)`.
//!
//! The Drude-Lorentz fit of Rakić et al. (1998) is used for gold and silver;
//! the parameters are exactly those tabulated with the refractiveindex.info
//! model scripts. All of these forms are causal, and a numerical
//! Kramers-Kronig check on the Drude response is one of the tests.
//!
//! Material constants are cited at the constructor. Wavelengths are metres in
//! the public API; the fits are internally in micrometres or electron-volts.
use crate::complex::C;

/// Photon energy hc/e in eV·m. CODATA 2018: h = 6.62607015e-34 J s,
/// c = 299792458 m/s, e = 1.602176634e-19 C.
pub const HC_EV_M: f64 = 1.2398419843320026e-6;

pub fn lambda_to_ev(lambda_m: f64) -> f64 {
    HC_EV_M / lambda_m
}
pub fn ev_to_lambda(ev: f64) -> f64 {
    HC_EV_M / ev
}

pub trait Dispersion {
    /// Complex refractive index `n + i kappa` at vacuum wavelength `lambda_m`.
    fn index(&self, lambda_m: f64) -> C;
}

// ---------------------------------------------------------------------------
// Sellmeier dielectrics
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct Sellmeier {
    pub name: &'static str,
    /// `B_i` and `C_i` (micrometres squared) pairs.
    pub terms: Vec<(f64, f64)>,
    pub valid_um: (f64, f64),
    pub source: &'static str,
}

impl Sellmeier {
    pub fn n2(&self, lambda_um: f64) -> f64 {
        let l2 = lambda_um * lambda_um;
        1.0 + self
            .terms
            .iter()
            .map(|&(b, c)| b * l2 / (l2 - c))
            .sum::<f64>()
    }
    pub fn n_um(&self, lambda_um: f64) -> f64 {
        self.n2(lambda_um).max(0.0).sqrt()
    }
    /// Fused silica, Malitson (1965), valid 0.21–3.71 µm, at 20 C.
    /// Measured n(589 nm) = 1.4585 (sodium D line).
    pub fn fused_silica() -> Self {
        Self {
            name: "fused silica",
            terms: vec![
                (0.6961663, 0.0684043_f64.powi(2)),
                (0.4079426, 0.1162414_f64.powi(2)),
                (0.8974794, 9.896161_f64.powi(2)),
            ],
            valid_um: (0.21, 3.71),
            source: "I. H. Malitson, JOSA 55, 1205 (1965), doi:10.1364/JOSA.55.001205",
        }
    }
    /// Distilled water at 20 C, Daimon & Masumura (2007), valid 0.182–1.129 µm.
    /// The four-term Sellmeier coefficients are those of the
    /// refractiveindex.info `Daimon-20.0C` entry.
    pub fn water() -> Self {
        Self {
            name: "water (20 C)",
            terms: vec![
                (0.5684027565, 0.005101829712),
                (0.1726177391, 0.01821153936),
                (0.02086189578, 0.02620722293),
                (0.1130748688, 10.69792721),
            ],
            valid_um: (0.182, 1.129),
            source: "M. Daimon & A. Masumura, Appl. Opt. 46, 3811 (2007), doi:10.1364/AO.46.003811",
        }
    }
}

impl Dispersion for Sellmeier {
    fn index(&self, lambda_m: f64) -> C {
        C::real(self.n_um(lambda_m * 1e6))
    }
}

// ---------------------------------------------------------------------------
// Drude-Lorentz metals
// ---------------------------------------------------------------------------

/// A bound-electron oscillator in the Lorentz sum. `omega0`, `gamma` and the
/// free-electron `omega_p` are in electron-volts (the Rakić convention).
#[derive(Clone, Copy, Debug)]
pub struct Oscillator {
    pub f: f64,
    pub gamma: f64,
    pub omega0: f64,
}

#[derive(Clone, Debug)]
pub struct DrudeLorentz {
    pub name: &'static str,
    pub omega_p: f64,
    /// Drude oscillator strength and damping (eV).
    pub f0: f64,
    pub gamma0: f64,
    pub oscillators: Vec<Oscillator>,
    pub source: &'static str,
}

impl DrudeLorentz {
    /// `epsilon(omega)` at photon energy in eV.
    pub fn epsilon_ev(&self, omega: f64) -> C {
        // Drude free-electron term: eps = 1 - (f0 wp^2)/(w(w + i gamma0)).
        let wp2 = self.omega_p * self.omega_p;
        let drude =
            C::real(self.f0 * wp2) / (C::real(omega) * (C::real(omega) + C::i() * self.gamma0));
        let mut eps = C::real(1.0) - drude;
        for o in &self.oscillators {
            let denom = C::real(o.omega0 * o.omega0 - omega * omega) - C::i() * (omega * o.gamma);
            eps = eps + C::real(o.f * wp2) / denom;
        }
        eps
    }
    /// `epsilon(omega)` at angular frequency `omega` in rad/s.
    pub fn epsilon_angular(&self, omega_rad: f64) -> C {
        assert!(omega_rad > 0.0);
        // Convert to the eV energy scale used by the fit.
        // hbar = 6.582119569e-16 eV s (CODATA 2018).
        let ev = 6.582_119_569e-16 * omega_rad;
        self.epsilon_ev(ev)
    }
    /// Lorentz-Drude fit for gold. Rakić et al. (1998), Table I.
    pub fn gold() -> Self {
        Self {
            name: "gold",
            omega_p: 9.03,
            f0: 0.760,
            gamma0: 0.053,
            oscillators: vec![
                Oscillator {
                    f: 0.024,
                    gamma: 0.241,
                    omega0: 0.415,
                },
                Oscillator {
                    f: 0.010,
                    gamma: 0.345,
                    omega0: 0.830,
                },
                Oscillator {
                    f: 0.071,
                    gamma: 0.870,
                    omega0: 2.969,
                },
                Oscillator {
                    f: 0.601,
                    gamma: 2.494,
                    omega0: 4.304,
                },
                Oscillator {
                    f: 4.384,
                    gamma: 2.214,
                    omega0: 13.32,
                },
            ],
            source: "A. D. Rakic et al., Appl. Opt. 37, 5271 (1998), doi:10.1364/AO.37.005271",
        }
    }
    /// Lorentz-Drude fit for silver. Rakić et al. (1998), Table I.
    pub fn silver() -> Self {
        Self {
            name: "silver",
            omega_p: 9.01,
            f0: 0.845,
            gamma0: 0.048,
            oscillators: vec![
                Oscillator {
                    f: 0.065,
                    gamma: 3.886,
                    omega0: 0.816,
                },
                Oscillator {
                    f: 0.124,
                    gamma: 0.452,
                    omega0: 4.481,
                },
                Oscillator {
                    f: 0.011,
                    gamma: 0.065,
                    omega0: 8.185,
                },
                Oscillator {
                    f: 0.840,
                    gamma: 0.916,
                    omega0: 9.083,
                },
                Oscillator {
                    f: 5.646,
                    gamma: 2.419,
                    omega0: 20.29,
                },
            ],
            source: "A. D. Rakic et al., Appl. Opt. 37, 5271 (1998), doi:10.1364/AO.37.005271",
        }
    }
}

impl Dispersion for DrudeLorentz {
    fn index(&self, lambda_m: f64) -> C {
        let eps = self.epsilon_ev(lambda_to_ev(lambda_m));
        eps.sqrt()
    }
}

/// Pure Drude permittivity with `hbar = 1` (energy units of `omega_p`, `gamma`).
pub fn drude(omega: f64, omega_p: f64, gamma: f64) -> C {
    C::real(1.0) - C::real(omega_p * omega_p) / (C::real(omega) * (C::real(omega) + C::i() * gamma))
}

/// Split-ring-resonator effective permeability (Pendry 1999),
///
/// ```text
///   mu(omega) = 1 - F omega^2 / (omega^2 - omega0^2 + i gamma omega),
/// ```
///
/// a causal Lorentz oscillator with resonance `omega0`, oscillator strength
/// `F` and damping `gamma` (all angular frequencies in rad/s; `F` is
/// dimensionless). The real part of `mu` is negative just above `omega0`,
/// between `omega0` and `omega0/sqrt(1 - F)`. The Drude permittivity
/// [`drude`] supplies the electric response of a wire medium; together they
/// give a negative-index band where `Re eps < 0` and `Re mu < 0`. With
/// `gamma > 0` the imaginary part is positive (a passive medium).
///
/// Source: J. B. Pendry, A. J. Holden, D. J. Robbins, W. J. Stewart,
/// "Magnetism from conductors and enhanced nonlinear phenomena",
/// IEEE Trans. Microwave Theory Tech. **47**, 2075 (1999),
/// doi:10.1109/22.798002 (effective-medium SRR form, eqs. 5-6 of the paper).
#[derive(Clone, Copy, Debug)]
pub struct SplitRing {
    pub omega0: f64,
    pub strength: f64,
    pub gamma: f64,
}

impl SplitRing {
    pub fn new(omega0: f64, strength: f64, gamma: f64) -> Self {
        assert!(omega0 > 0.0);
        Self {
            omega0,
            strength,
            gamma,
        }
    }
    pub fn mu(&self, omega: f64) -> C {
        let denom = C::real(omega * omega - self.omega0 * self.omega0)
            + C::i() * (self.gamma * omega);
        C::real(1.0) - C::real(self.strength * omega * omega) / denom
    }
    /// Permeability at vacuum wavelength `lambda_m`.
    pub fn mu_at_lambda(&self, lambda_m: f64) -> C {
        let omega = 2.0 * std::f64::consts::PI * crate::C0 / lambda_m;
        self.mu(omega)
    }
}

/// Single Lorentz oscillator permittivity (energy units of `omega0`, `gamma`).
pub fn lorentz(omega: f64, omega0: f64, gamma: f64, strength: f64) -> C {
    C::real(1.0)
        + C::real(strength) / (C::real(omega0 * omega0 - omega * omega) - C::i() * (omega * gamma))
}
