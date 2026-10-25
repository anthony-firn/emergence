use crate::{
    math::{Ray, Rng, Vec3},
    medium::{SpectralMedia, Volume},
    optics::sample_hg,
    C,
};
use rayon::prelude::*;

/// Incoherent energy packet (many photons). No c/dx explicit fluid CFL.
#[derive(Clone, Debug)]
pub struct Packet {
    pub ray: Ray,
    pub wavelength: f64,
    pub energy: f64,
    pub rng_state: u64,
}
impl Packet {
    pub fn photon_count(&self) -> f64 {
        self.energy * self.wavelength / (crate::H * C)
    }
}
#[derive(Clone, Debug)]
pub struct Deposit {
    pub cell: usize,
    pub energy: f64,
    pub momentum: Vec3,
}
#[derive(Clone, Debug)]
pub struct Flight {
    pub deposits: Vec<Deposit>,
    pub escaped_energy: f64,
    pub escaped_momentum: Vec3,
    pub live: Option<Packet>,
    /// Nonzero only when the event budget was reached: resume this time
    /// before advancing to the next fluid step. Energy is never discarded.
    pub pending_time: f64,
    pub events: usize,
}
/// Static straight-ray volume transport, implicit capture and HG scattering.
/// Vacuum speed c; refractive interfaces belong to the geometry transport.
/// Momentum transfer = incoming packet momentum - outgoing packet momentum.
pub fn fly(volume: &Volume, mut packet: Packet, dt: f64, event_budget: usize) -> Flight {
    assert!(
        packet.energy >= 0.0
            && packet.energy.is_finite()
            && packet.wavelength > 0.0
            && packet.wavelength.is_finite()
    );
    assert!(dt >= 0.0 && dt.is_finite() && event_budget > 0);
    let mut rng = Rng(packet.rng_state);
    let media = SpectralMedia::new(&volume.media, packet.wavelength);
    let mut remaining = C * dt;
    let mut out = Flight {
        deposits: vec![],
        escaped_energy: 0.0,
        escaped_momentum: Vec3::default(),
        live: None,
        pending_time: 0.0,
        events: 0,
    };
    while remaining > 0.0 && packet.energy > 0.0 {
        if out.events >= event_budget {
            out.pending_time = remaining / C;
            break;
        }
        let tau = -rng.uniform().ln();
        let mut accumulated = 0.0;
        let mut collision = None;
        for seg in volume.segments(packet.ray, remaining) {
            let voxel = volume.cells[seg.cell];
            let c = media.get(voxel.medium);
            let (a, s) = (c.a * voxel.scale, c.s * voxel.scale);
            let st = a + s;
            let depth = st * (seg.end - seg.start);
            if st > 0.0 && accumulated + depth > tau {
                collision = Some((seg.cell, seg.start + (tau - accumulated) / st, a, s));
                break;
            }
            accumulated += depth;
        }
        if let Some((cell, distance, a, s)) = collision {
            let before = packet.ray.dir * (packet.energy / C);
            let absorbed = packet.energy * a / (a + s);
            packet.energy -= absorbed;
            let medium = &volume.media[volume.cells[cell].medium];
            let dir = if s > 0.0 {
                sample_hg(packet.ray.dir, medium.g, &mut rng)
            } else {
                packet.ray.dir
            };
            out.deposits.push(Deposit {
                cell,
                energy: absorbed,
                momentum: before - dir * (packet.energy / C),
            });
            packet.ray = Ray::new(packet.ray.at(distance), dir);
            remaining -= distance;
            out.events += 1;
        } else {
            let exit = volume
                .bounds
                .interval(packet.ray, f64::INFINITY)
                .map(|p| p.1)
                .unwrap_or(0.0);
            if exit <= remaining {
                out.escaped_energy = packet.energy;
                out.escaped_momentum = packet.ray.dir * (packet.energy / C);
                packet.energy = 0.0;
            } else {
                packet.ray.origin = packet.ray.at(remaining);
            }
            remaining = 0.0;
        }
    }
    packet.rng_state = rng.0;
    if packet.energy > 0.0 {
        out.live = Some(packet);
    }
    out
}
/// Deterministic per-packet streams regardless of Rayon scheduling.
/// Each history writes private tallies; callers reduce/apply them once.
pub fn fly_batch(
    volume: &Volume,
    packets: Vec<Packet>,
    dt: f64,
    event_budget: usize,
) -> Vec<Flight> {
    packets
        .into_par_iter()
        .map(|p| fly(volume, p, dt, event_budget))
        .collect()
}

pub struct Emission {
    pub packets: Vec<Packet>,
    pub withdrawals: Vec<Deposit>,
}
/// Kirchhoff thermal emission in one narrow wavelength bin [m]. Its centre
/// quadrature is a spectral discretization, not a bolometric opacity.
/// Available thermal energy is [J/cell]. Reject, rather than clamp, a step
/// that exhausts it. Apply withdrawals to matter BEFORE flying the packets.
/// Isotropic emission recoil and packet momentum balance exactly per history.
pub fn emit_band(
    volume: &Volume,
    dt: f64,
    l: f64,
    width: f64,
    packets_per_cell: usize,
    available: &[f64],
    seed: u64,
) -> Result<Emission, &'static str> {
    if available.len() != volume.cells.len()
        || packets_per_cell == 0
        || !dt.is_finite()
        || dt < 0.0
        || !l.is_finite()
        || l <= 0.0
        || !width.is_finite()
        || width <= 0.0
    {
        return Err("invalid thermal emission inputs");
    }
    let dv = volume.cell_volume();
    let mut energy = vec![];
    for (q, v) in volume.cells.iter().enumerate() {
        let (a, _) = volume.coefficients(q, l);
        // This packet solver is for n=1 diffuse/gas transport.
        if (volume.media[v.medium].n.at(l) - 1.0).abs() > 1e-6 {
            return Err("thermal packets require n=1; refractive media need geometry transport");
        }
        let e = 4.0
            * std::f64::consts::PI
            * a
            * crate::optics::planck(l, v.temperature)
            * dv
            * dt
            * width;
        if !available[q].is_finite() || available[q] < 0.0 || !e.is_finite() || e > available[q] {
            return Err("emission exceeds available thermal energy; reduce step");
        }
        energy.push(e);
    }
    let mut rng = Rng(seed);
    let mut result = Emission {
        packets: vec![],
        withdrawals: vec![],
    };
    let dx = volume.dx();
    for (q, e) in energy.into_iter().enumerate() {
        if e == 0.0 {
            continue;
        }
        let i = [
            q / (volume.n[1] * volume.n[2]),
            q / volume.n[2] % volume.n[1],
            q % volume.n[2],
        ];
        let mut recoil = Vec3::default();
        for _ in 0..packets_per_cell {
            let p = volume.bounds.lo
                + Vec3(std::array::from_fn(|d| {
                    (i[d] as f64 + rng.uniform()) * dx.0[d]
                }));
            let dir = crate::math::around(
                Vec3([0.0, 0.0, 1.0]),
                2.0 * rng.uniform() - 1.0,
                2.0 * std::f64::consts::PI * rng.uniform(),
            );
            let packet = Packet {
                ray: Ray::new(p, dir),
                wavelength: l,
                energy: e / packets_per_cell as f64,
                rng_state: rng.next_u64(),
            };
            recoil = recoil - dir * (packet.energy / C);
            result.packets.push(packet);
        }
        result.withdrawals.push(Deposit {
            cell: q,
            energy: -e,
            momentum: recoil,
        });
    }
    Ok(result)
}
