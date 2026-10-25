//! Dynamic test: a 1D symmetric plate impact with the compressible GPR
//! elastic–plastic solid.
//!
//! The flyer and target are the same material. In the target's rest frame
//! the flyer moves at `V` and the target is at rest, so the particle
//! velocity behind the plastic shock is exactly `V/2` (symmetric impact).
//! Ahead of the plastic front an elastic precursor carries the Hugoniot
//! elastic limit `σ_HEL = Y(1−ν)/(1−2ν)`. The plastic shock speed is
//! compared with both the measured Hugoniot `Us = c0 + s·up` (Marsh 1980)
//! and the Rankine–Hugoniot speed of the NASG closure actually used by the
//! solver: a stiffened-gas EOS is only approximately linear in `up`, and
//! the report keeps the two errors separate.

use crate::Measurement;
use compressible::eos::Nasg;
use compressible::gpr::{Gpr, TauModel, NG};
use compressible::{Bc, Mesh, Solver};

use crate::{Level, Material};

/// NASG stiffened-gas `γ` that reproduces the measured weak-shock slope
/// `s = (γ+1)/4` at `up → 0`.
pub fn gamma_for_slope(s: f64) -> f64 {
    4.0 * s - 1.0
}

/// Rankine–Hugoniot shock speed of a stiffened gas:
/// `Us² − ((γ+1)/2) up Us − c0² = 0`.
pub fn rh_us(gamma: f64, c0: f64, up: f64) -> f64 {
    let a = 0.5 * (gamma + 1.0);
    let b = a * up;
    0.5 * (b + (b * b + 4.0 * c0 * c0).sqrt())
}

fn prim(rho: f64, v: [f64; 3], p: f64) -> [f64; NG] {
    let mut w = [0.0; NG];
    w[0] = rho;
    w[1..4].copy_from_slice(&v);
    for i in 0..3 {
        w[4 + 3 * i + i] = 1.0;
    }
    w[13] = p;
    w
}

/// First crossing of `level` in the particle velocity on the target side.
fn crossing(x: &[f64], v: &[f64], lo: usize, hi: usize, level: f64) -> f64 {
    for i in lo..hi - 1 {
        if (v[i] - level) * (v[i + 1] - level) <= 0.0 && v[i] != v[i + 1] {
            return x[i] + (level - v[i]) / (v[i + 1] - v[i]) * (x[i + 1] - x[i]);
        }
    }
    f64::NAN
}

fn band_mean<F: Fn(usize) -> f64>(x: &[f64], a: f64, b: f64, f: F) -> (f64, usize) {
    let (mut s, mut n) = (0.0, 0usize);
    for i in 0..x.len() {
        if x[i] > a && x[i] < b {
            s += f(i);
            n += 1;
        }
    }
    (s / n.max(1) as f64, n)
}

/// Everything the plate-impact test measures.
#[derive(Clone, Debug)]
pub struct ImpactResult {
    pub material: &'static str,
    pub impact_velocity: f64,
    pub steps: usize,
    /// Particle velocity behind the plastic shock (should be V/2).
    pub particle_velocity: f64,
    pub particle_velocity_ref: f64,
    /// Elastic-precursor axial stress (should be σ_HEL).
    pub hel_stress: f64,
    pub hel_stress_ref: f64,
    /// Elastic-precursor particle-velocity jump (should be σ_HEL/(ρ c_L)).
    pub hel_velocity: f64,
    pub hel_velocity_ref: f64,
    /// Measured elastic-precursor speed (should be c_L).
    pub precursor_speed: f64,
    pub longitudinal_speed: f64,
    /// Measured plastic shock speed.
    pub shock_speed: f64,
    /// Measured Hugoniot `c0 + s·up`.
    pub shock_speed_ref: f64,
    /// Rankine–Hugoniot speed of the NASG EOS that was actually solved.
    pub shock_speed_rh: f64,
    pub gamma: f64,
    pub bulk_speed: f64,
    pub shear_speed: f64,
    pub poisson: f64,
    pub yield_stress: f64,
}

/// Run the symmetric plate impact. `domain` is the half-length of the
/// target (and of the flyer); the grid has `n` cells in 1D.
pub fn plate_impact(mat: &Material, impact_velocity: f64, n: usize, domain: f64) -> ImpactResult {
    use materials::polycrystal::{hill, isotropic_parameters_stiffness};
    let moduli = hill(&mat.c);
    let ip = isotropic_parameters_stiffness(&mat.c, mat.density);
    let nu = moduli.poisson();
    let c_l = ip.longitudinal_speed;
    let cs = ip.shear_sound_speed;
    let c_bulk = ip.bulk_sound_speed;
    let gamma = gamma_for_slope(mat.hugoniot_s);
    let eos = Nasg::stiffened(gamma, mat.density * c_bulk * c_bulk / gamma);
    let sys = Gpr {
        eos,
        cs,
        tau: f64::INFINITY,
        tau_model: TauModel::Plastic {
            tau0: 1e-9,
            sigma_y: mat.yield_stress,
            n: 400.0,
        },
    };
    let sigma_hel = mat.yield_stress * (1.0 - nu) / (1.0 - 2.0 * nu);
    let u_hel = sigma_hel / (mat.density * c_l);
    let up = 0.5 * impact_velocity;
    let t_end = 0.75 * domain / c_l;

    let mesh = Mesh::new([n, 1, 1], [-domain, 0.0, 0.0], [domain, 1.0, 1.0]);
    let mut s = Solver::new(sys, mesh, |x| {
        if x[0] < 0.0 {
            prim(mat.density, [impact_velocity, 0.0, 0.0], 0.0)
        } else {
            prim(mat.density, [0.0, 0.0, 0.0], 0.0)
        }
    });
    s.bc[0] = [Bc::Outflow; 2];
    s.cfl = 0.4;
    s.strang = true;
    s.relax = Some(Box::new(|sys: &Gpr, w: &mut [f64; NG], dt: f64| sys.relax(w, dt)));
    let steps = s.run_to(t_end);

    let xs: Vec<f64> = (0..n).map(|i| s.mesh.x(i, 0, 0)[0]).collect();
    let vs: Vec<f64> = (0..n).map(|i| s.prim(i, 0, 0)[1]).collect();
    let half = n / 2;
    let x_el = crossing(&xs, &vs, half, n, 0.5 * u_hel);
    let x_pl = crossing(&xs, &vs, half, n, 0.5 * (u_hel + up));
    let shock_speed = x_pl / t_end;
    let precursor_speed = x_el / t_end;

    // Behind the plastic shock: between the interface and the plastic front.
    let (particle_velocity, _) = band_mean(&xs, 0.15 * x_pl, 0.85 * x_pl, |i| vs[i]);
    // Precursor plateau: the middle of the band between the two fronts.
    let (a, b) = (x_pl + 0.25 * (x_el - x_pl), x_pl + 0.75 * (x_el - x_pl));
    let (hel_velocity, _) = band_mean(&xs, a, b, |i| vs[i]);
    // Axial stress in the plateau: the compressive stress is
    // `−σ_xx = p − s_xx` with the GPR sign convention `σ_phys = −p δ + s`.
    let (hel_stress, _) = band_mean(&xs, a, b, |i| {
        let w = s.prim(i, 0, 0);
        let st = s.sys.stress(&w);
        w[13] - st[0][0]
    });

    ImpactResult {
        material: mat.name,
        impact_velocity,
        steps,
        particle_velocity,
        particle_velocity_ref: up,
        hel_stress,
        hel_stress_ref: sigma_hel,
        hel_velocity,
        hel_velocity_ref: u_hel,
        precursor_speed,
        longitudinal_speed: c_l,
        shock_speed,
        shock_speed_ref: mat.hugoniot_c0 + mat.hugoniot_s * up,
        shock_speed_rh: rh_us(gamma, c_bulk, up),
        gamma,
        bulk_speed: c_bulk,
        shear_speed: cs,
        poisson: nu,
        yield_stress: mat.yield_stress,
    }
}

impl ImpactResult {
    /// The dynamic-test measurements for the summary report.
    pub fn measurements(&self) -> Vec<Measurement> {
        let m = self.material;
        let mut out = Vec::new();
        let rel = |a: f64, b: f64| (a - b).abs() / b.abs().max(1e-300);
        out.push(Measurement::new(
            m,
            "plate impact (dynamic)",
            Level::Continuum,
            "particle velocity",
            self.particle_velocity,
            self.particle_velocity_ref,
            "m/s",
            rel(self.particle_velocity, self.particle_velocity_ref) < 0.01,
            "symmetric impact: exactly half the impact velocity",
        ));
        out.push(Measurement::new(
            m,
            "plate impact (dynamic)",
            Level::Continuum,
            "σ_HEL",
            self.hel_stress / 1e6,
            self.hel_stress_ref / 1e6,
            "MPa",
            rel(self.hel_stress, self.hel_stress_ref) < 0.10,
            "elastic precursor axial stress vs Y(1−ν)/(1−2ν); GPR power-law \
             hardening (n = 400) yields a few percent below ideal plasticity",
        ));
        out.push(Measurement::new(
            m,
            "plate impact (dynamic)",
            Level::Continuum,
            "u_HEL",
            self.hel_velocity,
            self.hel_velocity_ref,
            "m/s",
            rel(self.hel_velocity, self.hel_velocity_ref) < 0.10,
            "precursor particle-velocity jump vs σ_HEL/(ρ c_L)",
        ));
        out.push(Measurement::new(
            m,
            "plate impact (dynamic)",
            Level::Continuum,
            "precursor speed",
            self.precursor_speed,
            self.longitudinal_speed,
            "m/s",
            rel(self.precursor_speed, self.longitudinal_speed) < 0.03,
            "elastic precursor vs the polycrystal longitudinal speed; the \
             residual is finite-amplitude steepening plus grid diffusion",
        ));
        out.push(Measurement::new(
            m,
            "plate impact (dynamic)",
            Level::Continuum,
            "shock speed",
            self.shock_speed,
            self.shock_speed_ref,
            "m/s",
            rel(self.shock_speed, self.shock_speed_ref) < 0.02,
            "plastic shock vs the measured Hugoniot Us = c0 + s·up (Marsh \
             1980); the NASG/stiffened-gas EOS reproduces a linear Us–up \
             only approximately",
        ));
        out.push(Measurement::new(
            m,
            "plate impact (dynamic)",
            Level::Continuum,
            "EOS shock speed",
            self.shock_speed_rh,
            self.shock_speed_ref,
            "m/s",
            rel(self.shock_speed_rh, self.shock_speed_ref) < 0.02,
            &format!(
                "Rankine–Hugoniot speed of the NASG closure (γ = {:.3}, \
                 c0 = {:.0} m/s) vs the linear Marsh fit: the curvature of \
                 the stiffened-gas Hugoniot is the EOS-form error",
                self.gamma, self.bulk_speed
            ),
        ));
        out
    }
}
