//! Ideal, zero-chemical-potential quantum gas. This is a kinetic closure,
//! not an interacting QCD / electroweak EOS near transitions.
use crate::{
    species::{Domain, ALL},
    Species, GEV_J, HBAR_C,
};
use std::f64::consts::PI;
#[derive(Clone, Copy, Debug, Default)]
pub struct State {
    /// Natural units: number density [GeV³], pressure and energy [GeV⁴].
    pub number: f64,
    pub pressure: f64,
    pub energy: f64,
}
impl std::ops::Add for State {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            number: self.number + b.number,
            pressure: self.pressure + b.pressure,
            energy: self.energy + b.energy,
        }
    }
}
impl State {
    pub fn number_si(self) -> f64 {
        self.number / HBAR_C.powi(3)
    }
    pub fn energy_si(self) -> f64 {
        self.energy * GEV_J / HBAR_C.powi(3)
    }
    pub fn pressure_si(self) -> f64 {
        self.pressure * GEV_J / HBAR_C.powi(3)
    }
    /// Effective relativistic degrees of freedom, energy-based.
    pub fn g_star(self, temperature: f64) -> f64 {
        30.0 * self.energy / (PI * PI * temperature.powi(4))
    }
}
/// T = k_B T in GeV. Massless species use analytic Bose/Fermi integrals;
/// massive species use Simpson quadrature in p/T, with a Boltzmann tail.
pub fn species(s: Species, temperature: f64) -> Result<State, &'static str> {
    if !temperature.is_finite() || temperature <= 0.0 {
        return Err("temperature must be positive finite k_B T [GeV]");
    }
    let g = s.degeneracy();
    let fermion = s.fermion();
    let a = s.mass() / temperature;
    if a == 0.0 {
        let energy =
            g * PI * PI * temperature.powi(4) / 30.0 * if fermion { 7.0 / 8.0 } else { 1.0 };
        let number = g * 1.202_056_903_159_594 * temperature.powi(3) / (PI * PI)
            * if fermion { 0.75 } else { 1.0 };
        return Ok(State {
            number,
            pressure: energy / 3.0,
            energy,
        });
    }
    if a > 745.0 {
        return Ok(State::default());
    }
    // At small a the analytic limit would lose O(a²) mass corrections. Retain
    // the integral, extending the momentum interval for nonrelativistic gas.
    let hi = 50.0 + 12.0 * a.sqrt();
    let n = 4096;
    let h = hi / n as f64;
    let mut sum = State::default();
    for i in 1..n {
        // p=0 and far-tail endpoint both vanish to stated precision.
        let p = i as f64 * h;
        let e = p.hypot(a);
        let f = if fermion {
            1.0 / (e.exp() + 1.0)
        } else {
            1.0 / e.exp_m1()
        };
        let w = if i % 2 == 0 { 2.0 } else { 4.0 };
        sum.number += w * p * p * f;
        sum.energy += w * p * p * e * f;
        sum.pressure += w * p.powi(4) * f / (3.0 * e);
    }
    let prefactor = g * h / (6.0 * PI * PI);
    Ok(State {
        number: sum.number * prefactor * temperature.powi(3),
        energy: sum.energy * prefactor * temperature.powi(4),
        pressure: sum.pressure * prefactor * temperature.powi(4),
    })
}
/// Vacuum includes elementary colorless species; it omits hadrons entirely.
/// Partonic includes all colors/antiparticles once; valid away from confinement
/// only as a NONINTERACTING approximation. It is not the fluid matter EOS.
pub fn equilibrium(temperature: f64, domain: Domain) -> Result<State, &'static str> {
    if !temperature.is_finite() || temperature <= 0.0 {
        return Err("invalid temperature");
    }
    let mut state = State::default();
    for s in ALL {
        if domain == Domain::Partonic || !s.colored() {
            state = state + species(s, temperature)?;
        }
    }
    Ok(state)
}
