//! Inclusive decay closures: measured PDG W/Z/tau weights; dominant top
//! bW channel; CERN Report 3 Higgs weights at 125.09 GeV (mass here 125.20).
//! Weights are normalized together, never just over the resolved subset.
use crate::{
    kinematics,
    species::{Charges, Domain},
    FourMomentum, Rng, Species,
};
use Species::*;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ClosureReason {
    HadronicTauDecay,
    OffShellElectroweakDecay,
    /// These partons must be processed jointly by a QCD closure, not heated
    /// independently. The retained four-momentum includes rest energy.
    Hadronization,
}
#[derive(Clone, Debug)]
pub struct Unresolved {
    pub origin: Species,
    pub reason: ClosureReason,
    pub four: FourMomentum,
    pub charges: Charges,
    pub position: crate::Vec3,
    pub time: f64,
}
#[derive(Clone, Debug)]
pub enum Products {
    Two([Species; 2]),
    Leptonic([Species; 3]),
    Closure(ClosureReason),
}
#[derive(Clone, Debug)]
pub struct Channel {
    pub weight: f64,
    pub products: Products,
}
fn two(weight: f64, a: Species, b: Species) -> Channel {
    Channel {
        weight,
        products: Products::Two([a, b]),
    }
}
fn lepton(weight: f64, a: Species, b: Species, c: Species) -> Channel {
    Channel {
        weight,
        products: Products::Leptonic([a, b, c]),
    }
}
fn closure(weight: f64, reason: ClosureReason) -> Channel {
    Channel {
        weight,
        products: Products::Closure(reason),
    }
}
/// Public channel list makes approximations and unresolved probability visible.
pub fn channels(species: Species) -> Vec<Channel> {
    let conjugated = matches!(species, AntiMuon | AntiTau | AntiTop | WMinus);
    let s = if conjugated {
        species.conjugate()
    } else {
        species
    };
    let mut out = match s {
        Muon => vec![lepton(1.0, Electron, AntiElectronNeutrino, MuonNeutrino)],
        Tau => vec![
            lepton(0.1785, Electron, AntiElectronNeutrino, TauNeutrino),
            lepton(0.1737, Muon, AntiMuonNeutrino, TauNeutrino),
            closure(1.0 - 0.1785 - 0.1737, ClosureReason::HadronicTauDecay),
        ],
        Top => vec![two(1.0, Bottom, WPlus)],
        WPlus => {
            let mut v = vec![
                two(0.1071, Positron, ElectronNeutrino),
                two(0.1063, AntiMuon, MuonNeutrino),
                two(0.1138, AntiTau, TauNeutrino),
            ];
            // Approximate CKM magnitudes; inclusive hadronic weight is measured.
            let pairs = [
                (Up, AntiDown, 0.97435_f64),
                (Up, AntiStrange, 0.225),
                (Up, AntiBottom, 0.0037),
                (Charm, AntiDown, 0.221),
                (Charm, AntiStrange, 0.975),
                (Charm, AntiBottom, 0.041),
            ];
            let norm = pairs.iter().map(|p| p.2 * p.2).sum::<f64>();
            for (a, b, ckm) in pairs {
                v.push(two(0.6741 * ckm * ckm / norm, a, b));
            }
            v
        }
        Z => vec![
            two(0.033632, Electron, Positron),
            two(0.033662, Muon, AntiMuon),
            two(0.033696, Tau, AntiTau),
            two(0.2 / 3.0, ElectronNeutrino, AntiElectronNeutrino),
            two(0.2 / 3.0, MuonNeutrino, AntiMuonNeutrino),
            two(0.2 / 3.0, TauNeutrino, AntiTauNeutrino),
            two(0.116, Up, AntiUp),
            two(0.1203, Charm, AntiCharm),
            two(0.1512, Bottom, AntiBottom),
            two((0.69911 - 0.116 - 0.1203 - 0.1512) / 2.0, Down, AntiDown),
            two(
                (0.69911 - 0.116 - 0.1203 - 0.1512) / 2.0,
                Strange,
                AntiStrange,
            ),
        ],
        Higgs => {
            let mut v = vec![
                two(0.575, Bottom, AntiBottom),
                two(0.063, Tau, AntiTau),
                two(0.000219, Muon, AntiMuon),
                two(0.029, Charm, AntiCharm),
                two(0.00228, Photon, Photon),
                two(0.00155, Z, Photon),
                closure(0.216 + 0.0266, ClosureReason::OffShellElectroweakDecay),
            ];
            // Eight adjoint components, color averaged; no shower evolution.
            for g in crate::species::GLUONS {
                v.push(two(0.0856 / 8.0, g, g));
            }
            v
        }
        _ => vec![],
    };
    if conjugated {
        for c in &mut out {
            match &mut c.products {
                Products::Two(p) => p.iter_mut().for_each(|s| *s = s.conjugate()),
                Products::Leptonic(p) => p.iter_mut().for_each(|s| *s = s.conjugate()),
                Products::Closure(_) => {}
            }
        }
    }
    let total = out.iter().map(|c| c.weight).sum::<f64>();
    if total > 0.0 {
        for c in &mut out {
            c.weight /= total;
        }
    }
    out
}
#[derive(Default)]
pub struct Decay {
    pub particles: Vec<(Species, FourMomentum)>,
    pub unresolved: Vec<(Species, ClosureReason, FourMomentum, Charges)>,
}
/// Decays an on-shell, unpolarized parent. Domain controls colored outputs.
pub fn sample(
    species: Species,
    four: FourMomentum,
    domain: Domain,
    rng: &mut Rng,
) -> Result<Decay, &'static str> {
    four.validate_on_shell(species.mass())?;
    let table = channels(species);
    if table.is_empty() {
        return Err("no free-particle decay closure for this species");
    }
    let mut x = rng.uniform();
    let mut chosen = &table[table.len() - 1];
    for c in &table {
        x -= c.weight;
        if x <= 0.0 {
            chosen = c;
            break;
        }
    }
    let mut result = Decay::default();
    let products = match &chosen.products {
        Products::Closure(reason) => {
            result
                .unresolved
                .push((species, *reason, four, species.charges()));
            return Ok(result);
        }
        Products::Two(p) => {
            let k = kinematics::two_body(four, species.mass(), p[0].mass(), p[1].mass(), rng)?;
            p.iter().copied().zip(k).collect::<Vec<_>>()
        }
        Products::Leptonic(p) => {
            let k = kinematics::leptonic(four, species.mass(), p[0].mass(), rng)?;
            p.iter().copied().zip(k).collect::<Vec<_>>()
        }
    };
    for (s, p) in products {
        if s.allowed(domain) {
            result.particles.push((s, p));
        } else {
            result
                .unresolved
                .push((s, ClosureReason::Hadronization, p, s.charges()));
        }
    }
    Ok(result)
}
