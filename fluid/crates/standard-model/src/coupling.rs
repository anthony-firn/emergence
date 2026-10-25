//! Explicit conversion to the existing SI radiation/fluid source interface.
use crate::{FourMomentum, Particle, Species, C, GEV_J};
use radiation::{
    math::Ray,
    transport::{Deposit, Packet},
};
/// One physical photon -> energy packet with photon_count=1. Packets containing
/// many photons must be split/sample weighted histories; never reinterpret their
/// entire energy as one high-energy collision photon.
pub fn photon_packet(p: &Particle, seed: u64) -> Result<Packet, &'static str> {
    if p.species() != Species::Photon {
        return Err("only a photon converts to a radiation packet");
    }
    let energy = p.four().energy * GEV_J;
    let wavelength = radiation::H * C / energy;
    if !wavelength.is_finite() || wavelength <= 0.0 || !energy.is_finite() || energy <= 0.0 {
        return Err("unrepresentable photon energy or wavelength");
    }
    Ok(Packet {
        ray: Ray::new(p.position(), p.four().momentum),
        wavelength,
        energy,
        rng_state: seed,
    })
}
/// Signed transferred TOTAL energy and momentum, not an automatic absorption
/// law. Rest energy, escape of neutrinos and hadronization must be accounted for
/// by the caller's physical interaction closure before applying this source.
pub fn deposit(cell: usize, transfer: FourMomentum) -> Result<Deposit, &'static str> {
    let energy = transfer.energy * GEV_J;
    let momentum = transfer.momentum * (GEV_J / C);
    if !energy.is_finite() || momentum.0.iter().any(|x| !x.is_finite()) {
        return Err("nonfinite particle transfer");
    }
    Ok(Deposit {
        cell,
        energy,
        momentum,
    })
}
