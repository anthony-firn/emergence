//! Explicit adapters: read-only snapshots for cameras; conservative source
//! updates for radiation dynamics. Call sources as split fluid steps.
use crate::{Solver, System};
use radiation::{
    math::{Bounds, Vec3},
    medium::{Medium, Volume, Voxel},
    thermal,
    transport::Deposit,
};

/// Map each solver primitive cell to calibrated optical coefficients and T.
/// The callback chooses medium and density/volume-fraction scaling; mechanical
/// EOS parameters alone cannot supply measured optical properties.
pub fn snapshot<const N: usize, S: System<N>>(
    s: &Solver<N, S>,
    media: Vec<Medium>,
    map: impl Fn(&S, &[f64; N]) -> Voxel,
) -> Volume {
    let m = &s.mesh;
    let mut cells = Vec::with_capacity(m.cells());
    for i in 0..m.n[0] {
        for j in 0..m.n[1] {
            for k in 0..m.n[2] {
                cells.push(map(&s.sys, &s.prim(i, j, k)));
            }
        }
    }
    Volume::new(
        m.n,
        Bounds {
            lo: Vec3(m.lo),
            hi: Vec3(std::array::from_fn(|d| m.lo[d] + m.n[d] as f64 * m.dx[d])),
        },
        media,
        cells,
    )
}
fn interior<const N: usize, S: System<N>>(s: &Solver<N, S>, q: usize) -> usize {
    let m = &s.mesh;
    let i = q / (m.n[1] * m.n[2]);
    let j = q / m.n[2] % m.n[1];
    let k = q % m.n[2];
    m.idx(i + m.g(0), j + m.g(1), k + m.g(2))
}
/// Photon energy is TOTAL energy: momentum transfer consumes the corresponding
/// kinetic energy from that same increment, rather than double counting it.
/// All updates are checked before committing any cell. Returns energy [J].
pub fn deposit<const N: usize, S: System<N>>(
    s: &mut Solver<N, S>,
    tallies: &[Deposit],
) -> Result<f64, &'static str> {
    let ie = s.sys.energy_index().ok_or("system has no total energy")?;
    let dv = s.mesh.dx.iter().product::<f64>();
    let mut updates = std::collections::BTreeMap::<usize, [f64; 4]>::new();
    let mut total = 0.0;
    for t in tallies {
        if t.cell >= s.mesh.cells()
            || !t.energy.is_finite()
            || t.momentum.0.iter().any(|v| !v.is_finite())
        {
            return Err("invalid photon tally");
        }
        let v = updates.entry(t.cell).or_default();
        v[0] += t.energy / dv;
        for d in 0..3 {
            v[d + 1] += t.momentum.0[d] / dv;
        }
        total += t.energy;
    }
    let mut prepared = vec![];
    for (q, v) in updates {
        let idx = interior(s, q);
        let mut u = s.u[idx];
        u[ie] += v[0];
        for d in 0..3 {
            u[s.sys.momentum_index(d)] += v[d + 1];
        }
        if !s.sys.admissible(&s.sys.to_prim(&u)) || u.iter().any(|v| !v.is_finite()) {
            return Err("photon source produces inadmissible fluid state; reduce split step");
        }
        prepared.push((idx, u));
    }
    for (idx, u) in prepared {
        s.u[idx] = u;
    }
    s.allspeed_p.clear();
    Ok(total)
}
/// Grey n=1 absorption/emission with frozen volumetric Cv. Valid for thermal
/// equilibrium within each cell, fixed phase fractions, and negligible v/c.
/// Does not transport radiation: pair with diffusion or angular packets.
/// Apply phase relaxation after this source if a transition is crossed.
/// System::phases must supply the true fixed-volume heat capacity; a generic
/// approximate/default closure is not a calibrated radiation thermodynamics.
pub fn exchange_grey<const N: usize, S: System<N>>(
    s: &mut Solver<N, S>,
    e: &mut [f64],
    planck_opacity: &[f64],
    dt: f64,
) -> Result<f64, &'static str> {
    if e.len() != s.mesh.cells() || e.len() != planck_opacity.len() || !dt.is_finite() || dt < 0.0 {
        return Err("invalid radiation field/time step");
    }
    let ie = s.sys.energy_index().ok_or("system has no total energy")?;
    let mut prepared = vec![];
    let mut sum = 0.0;
    for q in 0..e.len() {
        if !e[q].is_finite()
            || e[q] < 0.0
            || !planck_opacity[q].is_finite()
            || planck_opacity[q] < 0.0
        {
            return Err("invalid radiation energy/opacity");
        }
        let idx = interior(s, q);
        let w = s.sys.to_prim(&s.u[idx]);
        let phases = s.sys.phases(&w);
        let t = s.sys.temperature(&w);
        let cv: f64 = phases.iter().map(|p| p.0 * p.2).sum();
        if !t.is_finite() || t < 0.0 || !cv.is_finite() || cv <= 0.0 {
            return Err("system has no valid thermal closure");
        }
        if phases
            .iter()
            .any(|p| p.0 > 1e-6 && (!p.1.is_finite() || (p.1 - t).abs() > 1e-4 * t.max(1.0)))
        {
            return Err("grey exchange requires thermal equilibrium within each cell");
        }
        let r = thermal::exchange(t, cv, e[q], planck_opacity[q], dt);
        let mut u = s.u[idx];
        u[ie] += r.matter_delta;
        let updated = s.sys.to_prim(&u);
        if !s.sys.admissible(&updated) || u.iter().any(|v| !v.is_finite()) {
            return Err("radiative exchange produces inadmissible fluid state");
        }
        let actual_t = s.sys.temperature(&updated);
        if !actual_t.is_finite()
            || (actual_t - r.temperature).abs() > 1e-8 * r.temperature.max(1.0)
            || s.sys.phases(&updated).iter().any(|p| {
                p.0 > 1e-6
                    && (!p.1.is_finite()
                        || (p.1 - r.temperature).abs() > 1e-4 * r.temperature.max(1.0))
            })
        {
            return Err("grey source leaves its fixed-Cv thermal-equilibrium closure; use a coupled EOS/phase source");
        }
        sum += r.matter_delta;
        prepared.push((idx, u, r.radiation));
    }
    for (q, (idx, u, en)) in prepared.into_iter().enumerate() {
        s.u[idx] = u;
        e[q] = en;
    }
    s.allspeed_p.clear();
    Ok(sum * s.mesh.dx.iter().product::<f64>())
}
