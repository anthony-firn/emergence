//! Elastic neutrino-electron scattering in the low-energy Fermi theory.
//! Finite electron recoil; massless neutrinos; no radiative corrections or
//! inelastic channels. External couplings can calibrate the effective theory.
use crate::{radiation_math, FourMomentum, Rng, Species, HBAR_C};
use std::f64::consts::PI;
pub const FERMI_GEV2: f64 = 1.166_378_7e-5;
#[derive(Clone, Copy, Debug)]
pub struct Couplings {
    pub left: f64,
    pub right: f64,
    pub fermi: f64,
}
impl Couplings {
    /// Tree-level chiral coefficients, with electron-flavor charged current.
    /// sin²(theta_W) and its renormalization convention are caller inputs.
    pub fn tree(neutrino: Species, sin2_theta: f64) -> Result<Self, &'static str> {
        if !neutrino.neutrino() || !sin2_theta.is_finite() || !(0.0..1.0).contains(&sin2_theta) {
            return Err("invalid neutrino/weak angle");
        }
        let electron = neutrino.pdg_id().abs() == 12;
        let (mut left, mut right) = (
            sin2_theta - 0.5 + if electron { 1.0 } else { 0.0 },
            sin2_theta,
        );
        if neutrino.pdg_id() < 0 {
            std::mem::swap(&mut left, &mut right);
        }
        Ok(Self {
            left,
            right,
            fermi: FERMI_GEV2,
        })
    }
    fn validate(self) -> Result<(), &'static str> {
        if !self.left.is_finite()
            || !self.right.is_finite()
            || !self.fermi.is_finite()
            || self.fermi <= 0.0
        {
            return Err("invalid weak effective coefficients");
        }
        Ok(())
    }
}
fn domain(energy: f64) -> Result<(), &'static str> {
    let m = Species::Electron.mass();
    let mw = Species::WPlus.mass();
    if !energy.is_finite() || energy <= 0.0 || 2.0 * m * energy / (mw * mw) > 1e-3 {
        return Err("neutrino energy outside low-transfer Fermi closure");
    }
    Ok(())
}
fn endpoint(energy: f64) -> f64 {
    2.0 * energy * energy / (Species::Electron.mass() + 2.0 * energy)
}
fn shape(c: Couplings, energy: f64, t: f64) -> f64 {
    let m = Species::Electron.mass();
    (c.left * c.left + c.right * c.right * (1.0 - t / energy).powi(2)
        - c.left * c.right * m * t / (energy * energy))
        .max(0.0)
}
/// d sigma / d T [m²/GeV], with target-rest energy and electron kinetic recoil.
pub fn differential(c: Couplings, energy: f64, recoil: f64) -> Result<f64, &'static str> {
    c.validate()?;
    domain(energy)?;
    if !recoil.is_finite() || recoil < 0.0 || recoil > endpoint(energy) {
        return Err("electron recoil outside kinematic support");
    }
    Ok(2.0 * c.fermi * c.fermi * Species::Electron.mass() / PI
        * HBAR_C
        * HBAR_C
        * shape(c, energy, recoil))
}
/// Analytic integrated elastic cross section [m²]. Caller schedules encounters.
pub fn cross_section(c: Couplings, energy: f64) -> Result<f64, &'static str> {
    c.validate()?;
    domain(energy)?;
    let m = Species::Electron.mass();
    let t = endpoint(energy);
    let y = t / energy;
    let integral = t
        * (c.left * c.left + c.right * c.right * (1.0 - y + y * y / 3.0)
            - c.left * c.right * m * t / (2.0 * energy * energy));
    Ok(2.0 * c.fermi * c.fermi * m / PI * HBAR_C * HBAR_C * integral.max(0.0))
}
/// Return outgoing neutrino/electron in the original frame. The energy and
/// recoil angular kernel is sampled in the electron rest frame. Electron
/// flavor/antiflavor selection is encoded by the supplied chiral coefficients.
pub fn scatter(
    neutrino: FourMomentum,
    electron: FourMomentum,
    c: Couplings,
    rng: &mut Rng,
) -> Result<[FourMomentum; 2], &'static str> {
    c.validate()?;
    neutrino.validate_on_shell(0.0)?;
    let m = Species::Electron.mass();
    electron.validate_on_shell(m)?;
    let incoming = neutrino.to_rest_on_shell(0.0, electron, m);
    let energy = incoming.energy;
    domain(energy)?;
    let maximum = endpoint(energy);
    // Quadratic, convex in recoil: max on either endpoint. This remains
    // efficient for both neutrinos and antineutrinos including interference.
    let bound = shape(c, energy, 0.0).max(shape(c, energy, maximum));
    if bound == 0.0 {
        return Err("zero elastic weak coupling");
    }
    let recoil = loop {
        let t = maximum * rng.uniform();
        if rng.uniform() * bound <= shape(c, energy, t) {
            break t;
        }
    };
    let outgoing_energy = energy - recoil;
    let cosine = (1.0 - m * recoil / (energy * outgoing_energy)).clamp(-1.0, 1.0);
    let direction =
        radiation_math::around(incoming.momentum.unit(), cosine, 2.0 * PI * rng.uniform());
    let out = FourMomentum {
        energy: outgoing_energy,
        momentum: direction * outgoing_energy,
    };
    let electron_out = FourMomentum::on_shell(m, incoming.momentum - out.momentum)?;
    Ok([
        out.from_rest_on_shell(0.0, electron, m),
        electron_out.from_rest_on_shell(m, electron, m),
    ])
}
