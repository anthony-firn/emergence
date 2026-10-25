//! Regime map: local indicators of where each model in the hierarchy is
//! valid — the "emergent boundaries" at which more (or less) detail is
//! needed. The global certificate says *whether* a run is resolved; these
//! say *where* it is not, and where the continuum description itself fails.
//!
//! * truncation density: |P(v × ω)| restricted to the modes the 2/3 rule
//!   discards, in physical space — the local source of the certificate's
//!   residual (where the grid stops representing the dynamics);
//! * local Kolmogorov resolution k_max·η(x), η = (ν³/ε)^{1/4} with the
//!   local dissipation ε = 2ν|S|² (DNS practice wants ≳ 1.5);
//! * continuum validity, in physical units of a material (`Material`):
//!   the non-equilibrium parameter λ|∇u|/c (molecular collision time ×
//!   strain rate; a gradient-length Knudsen number λ|∇u|/|u| would blow up
//!   at stagnation points), the Mach number |u|/c (incompressibility), and
//!   λ/Δx (resolving below molecular scales is meaningless);
//! * flagged regions: periodic connected components above a threshold, with
//!   centroids — where to put a refinement patch.

use num_complex::Complex64;
use rayon::prelude::*;

use crate::fft::Fft3;
use crate::grid::{wavenumber, Grid};
use crate::solver::{zeros_spec, Solver, Spec};

const I: Complex64 = Complex64::new(0.0, 1.0);

/// Molecular-scale properties for continuum checks (SI units).
#[derive(Clone, Copy, Debug)]
pub struct Material {
    pub name: &'static str,
    /// mean free path (gas) or molecular spacing (liquid), m
    pub mean_free_path: f64,
    /// speed of sound, m/s
    pub sound_speed: f64,
    pub density: f64,
    /// dynamic viscosity, Pa·s
    pub viscosity: f64,
    /// surface tension against its vapour / air, N/m (0 for gases)
    pub surface_tension: f64,
    /// saturated vapour pressure at the reference temperature, Pa
    pub vapor_pressure: f64,
    /// electrical conductivity, S/m
    pub conductivity: f64,
    /// relative permittivity
    pub permittivity: f64,
}

/// Vacuum permeability and permittivity (SI).
pub const MU_0: f64 = 1.256_637_062e-6;
pub const EPS_0: f64 = 8.854_187_813e-12;

impl Material {
    /// 20 °C, 1 atm. Conductivity of pure water; tap water is ~5e-2 S/m.
    pub const WATER: Material = Material {
        name: "water",
        mean_free_path: 0.3e-9,
        sound_speed: 1482.0,
        density: 998.2,
        viscosity: 1.002e-3,
        surface_tension: 0.0728,
        vapor_pressure: 2339.0,
        conductivity: 5.5e-6,
        permittivity: 80.1,
    };
    pub const NITROGEN_GAS: Material = Material {
        name: "nitrogen gas",
        mean_free_path: 66e-9,
        sound_speed: 349.0,
        density: 1.165,
        viscosity: 1.76e-5,
        surface_tension: 0.0,
        vapor_pressure: f64::INFINITY,
        conductivity: 1e-15,
        permittivity: 1.00055,
    };
    /// at its boiling point, 77 K
    pub const LIQUID_NITROGEN: Material = Material {
        name: "liquid nitrogen",
        mean_free_path: 0.35e-9,
        sound_speed: 853.0,
        density: 807.0,
        viscosity: 1.58e-4,
        surface_tension: 0.00885,
        vapor_pressure: 101_325.0,
        conductivity: 1e-15,
        permittivity: 1.43,
    };
    pub const AIR: Material = Material {
        name: "air",
        mean_free_path: 68e-9,
        sound_speed: 343.0,
        density: 1.204,
        viscosity: 1.81e-5,
        surface_tension: 0.0,
        vapor_pressure: f64::INFINITY,
        conductivity: 1e-14,
        permittivity: 1.00059,
    };
    /// 20 °C. A metal: permittivity is not meaningful (charge relaxes in
    /// ~1e-19 s); 1 is used.
    pub const MERCURY: Material = Material {
        name: "mercury",
        mean_free_path: 0.3e-9,
        sound_speed: 1451.0,
        density: 13_534.0,
        viscosity: 1.526e-3,
        surface_tension: 0.485,
        vapor_pressure: 0.17,
        conductivity: 1.04e6,
        permittivity: 1.0,
    };

    pub fn kinematic_viscosity(&self) -> f64 {
        self.viscosity / self.density
    }
}

/// Dimensionless groups of a flow: each one marks an emergent boundary of
/// some model on the ladder.
#[derive(Clone, Copy, Debug)]
pub struct FlowNumbers {
    pub reynolds: f64,
    /// U/c: incompressible below ~0.3; shocks possible above 1
    pub mach: f64,
    /// ρU²L/γ: surface tension dominates below ~1
    pub weber: f64,
    /// (p − p_v)/(½ρU²): cavitation inception when it drops to O(1)
    pub cavitation: f64,
    /// μ₀σUL: the applied field is unperturbed (quasi-static MHD) when ≪ 1
    pub magnetic_reynolds: f64,
    /// BL√(σ/μ): Hartmann layers have thickness L/Ha
    pub hartmann: f64,
    /// ε₀ε_r/σ ÷ L/U: charges relax instantly (ohmic) when ≪ 1; a perfect
    /// dielectric when ≫ 1; leaky-dielectric EHD in between
    pub charge_relaxation: f64,
}

impl FlowNumbers {
    pub fn new(mat: &Material, length: f64, velocity: f64, b_field: f64, pressure: f64) -> Self {
        let q = 0.5 * mat.density * velocity * velocity;
        FlowNumbers {
            reynolds: velocity * length / mat.kinematic_viscosity(),
            mach: velocity / mat.sound_speed,
            weber: if mat.surface_tension > 0.0 { mat.density * velocity * velocity * length / mat.surface_tension } else { f64::INFINITY },
            cavitation: if mat.vapor_pressure.is_finite() { (pressure - mat.vapor_pressure) / q } else { f64::INFINITY },
            magnetic_reynolds: MU_0 * mat.conductivity * velocity * length,
            hartmann: b_field * length * (mat.conductivity / mat.viscosity).sqrt(),
            charge_relaxation: EPS_0 * mat.permittivity / mat.conductivity * velocity / length,
        }
    }

    /// The models valid at these numbers (coarsest first) and the ones that
    /// must be unfolded.
    pub fn regimes(&self) -> Vec<(&'static str, bool)> {
        vec![
            ("incompressible (Ma < 0.3)", self.mach < 0.3),
            ("no shocks (Ma < 1)", self.mach < 1.0),
            ("no cavitation (σ_c > 1)", self.cavitation > 1.0),
            ("quasi-static low-Rm MHD (Rm < 0.1)", self.magnetic_reynolds < 0.1),
            ("ohmic conductor (τ_e U/L < 0.01)", self.charge_relaxation < 0.01),
            ("perfect dielectric (τ_e U/L > 100)", self.charge_relaxation > 100.0),
        ]
    }
}

/// Connected region where an indicator exceeds its threshold.
#[derive(Clone, Debug)]
pub struct Region {
    pub points: usize,
    pub centroid: [f64; 3],
    pub max: f64,
}

pub struct RegimeMap {
    pub n: usize,
    /// local truncation-residual magnitude (N grid)
    pub truncation: Vec<f64>,
    /// k_max · η(x)
    pub kolmogorov: Vec<f64>,
    /// |ω(x)| (the naive indicator, for comparison)
    pub vorticity: Vec<f64>,
}

fn smooth_eps(x: f64) -> f64 {
    x.max(1e-300)
}

impl Solver {
    /// Compute the regime map of the current state.
    pub fn regime_map(&mut self) -> RegimeMap {
        let g = self.grid.clone();
        let n = g.n;
        let km = g.kmax() as i64;
        // padded grid holding v × ω exactly
        let mut m = (3 * n) / 2;
        if m < (4 * km + 2) as usize {
            m = (4 * km + 2) as usize;
        }
        m += m % 2;
        let gp = Grid::new(m);
        let fp = Fft3::new(&gp);
        let mut pairs = Vec::new();
        for idx in 0..g.spec_len() {
            if !g.kept(idx) {
                continue;
            }
            let (i, j, l) = g.spec_ijl(idx);
            let w = |k: i64| k.rem_euclid(m as i64) as usize;
            pairs.push(((w(wavenumber(i, n)) * m + w(wavenumber(j, n))) * gp.nh + l, idx));
        }
        let mut spec = zeros_spec(gp.spec_len());
        let mut tmp = vec![Complex64::default(); gp.spec_len()];
        let nrp = gp.real_len();
        let mut u: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; nrp]);
        let mut w: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; nrp]);
        for c in 0..3 {
            for &(pi, ni) in &pairs {
                spec[c][pi] = self.u[c][ni];
            }
            fp.inverse(&spec[c], &mut tmp, &mut u[c]);
        }
        for c in 0..3 {
            let (a, b) = ((c + 1) % 3, (c + 2) % 3);
            tmp.iter_mut().for_each(|v| *v = Complex64::default());
            for &(pi, _) in &pairs {
                let k = gp.k(pi);
                tmp[pi] = I * (spec[b][pi] * k[a] - spec[a][pi] * k[b]);
            }
            fp.inverse_inplace(&mut tmp, &mut w[c]);
        }
        let mut prod: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; nrp]);
        for q in 0..nrp {
            prod[0][q] = u[1][q] * w[2][q] - u[2][q] * w[1][q];
            prod[1][q] = u[2][q] * w[0][q] - u[0][q] * w[2][q];
            prod[2][q] = u[0][q] * w[1][q] - u[1][q] * w[0][q];
        }
        let mut ps = zeros_spec(gp.spec_len());
        for c in 0..3 {
            fp.forward(&prod[c], &mut ps[c]);
        }
        // keep only the discarded, solenoidal part
        for idx in 0..gp.spec_len() {
            let (i, j, l) = gp.spec_ijl(idx);
            let (kx, ky, kz) = (wavenumber(i, m), wavenumber(j, m), l as i64);
            let k = gp.k(idx);
            let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
            if (kx.abs() <= km && ky.abs() <= km && kz <= km) || k2 == 0.0 {
                for c in 0..3 {
                    ps[c][idx] = Complex64::default();
                }
                continue;
            }
            let d = (ps[0][idx] * k[0] + ps[1][idx] * k[1] + ps[2][idx] * k[2]) / k2;
            for c in 0..3 {
                ps[c][idx] -= d * k[c];
            }
        }
        let mut tr: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; nrp]);
        for c in 0..3 {
            fp.inverse_inplace(&mut ps[c], &mut tr[c]);
        }
        // sample the padded fields at the N-grid points (nearest padded point)
        let to_p = |i: usize| ((i as f64 * m as f64 / n as f64).round() as usize) % m;
        let nr = g.real_len();
        let mut truncation = vec![0.0; nr];
        let mut vorticity = vec![0.0; nr];
        for p in 0..nr {
            let (i, j, l) = (p / (n * n), (p / n) % n, p % n);
            let q = (to_p(i) * m + to_p(j)) * m + to_p(l);
            truncation[p] = (tr[0][q].powi(2) + tr[1][q].powi(2) + tr[2][q].powi(2)).sqrt();
            vorticity[p] = (w[0][q].powi(2) + w[1][q].powi(2) + w[2][q].powi(2)).sqrt();
        }
        // local Kolmogorov resolution from |S|² on the N grid
        let s2 = self.strain_squared();
        let nu = self.nu;
        let kmaxf = g.kmax() as f64;
        let kolmogorov = s2
            .par_iter()
            .map(|&x| {
                let eps = 2.0 * nu * x;
                kmaxf * (nu.powi(3) / smooth_eps(eps)).powf(0.25)
            })
            .collect();
        RegimeMap { n, truncation, kolmogorov, vorticity }
    }

    /// |S|² at the N-grid points.
    fn strain_squared(&mut self) -> Vec<f64> {
        let g = self.grid.clone();
        let comps = [(0, 0, 1.0), (1, 1, 1.0), (2, 2, 1.0), (0, 1, 2.0), (0, 2, 2.0), (1, 2, 2.0)];
        let mut acc = vec![0.0; g.real_len()];
        let mut sp: Spec = vec![Complex64::default(); g.spec_len()];
        let mut out = vec![0.0; g.real_len()];
        for &(i, j, mult) in &comps {
            for idx in 0..g.spec_len() {
                let k = g.k(idx);
                sp[idx] = 0.5 * I * (self.u[j][idx] * k[i] + self.u[i][idx] * k[j]);
            }
            self.fft.inverse_inplace(&mut sp, &mut out);
            acc.par_iter_mut().zip(out.par_iter()).for_each(|(a, &s)| *a += mult * s * s);
        }
        acc
    }

    /// Continuum validity in physical units: the box side is `length` (m),
    /// velocities are scaled by `velocity` (m/s).
    pub fn continuum_check(&mut self, mat: &Material, length: f64, velocity: f64) -> ContinuumCheck {
        let g = self.grid.clone();
        let scale_x = length / (2.0 * std::f64::consts::PI);
        let s2 = self.strain_squared();
        let comps: Vec<Vec<f64>> = (0..3).map(|c| self.velocity_phys(c)).collect();
        let umax = (0..g.real_len())
            .map(|p| (comps[0][p].powi(2) + comps[1][p].powi(2) + comps[2][p].powi(2)).sqrt())
            .fold(0.0, f64::max);
        let grad_max = s2.iter().cloned().fold(0.0, f64::max).sqrt() * 2f64.sqrt() * velocity / scale_x;
        ContinuumCheck {
            non_equilibrium: mat.mean_free_path * grad_max / mat.sound_speed,
            mach: umax * velocity / mat.sound_speed,
            mfp_over_dx: mat.mean_free_path / (g.dx() * scale_x),
        }
    }
}

/// Where the continuum (incompressible Navier–Stokes) model itself holds.
#[derive(Clone, Copy, Debug)]
pub struct ContinuumCheck {
    /// λ |∇u|_max / c: collision time × strain rate (≪ 1 for continuum)
    pub non_equilibrium: f64,
    /// |u|_max / c (≲ 0.3 for incompressibility)
    pub mach: f64,
    /// λ / Δx (≪ 1, or the grid resolves below molecular scales)
    pub mfp_over_dx: f64,
}

impl ContinuumCheck {
    pub fn valid(&self) -> bool {
        self.non_equilibrium < 0.01 && self.mach < 0.3 && self.mfp_over_dx < 0.1
    }
}

impl RegimeMap {
    /// Periodic connected components of `field > threshold` (6-neighbour).
    pub fn regions(field: &[f64], n: usize, threshold: f64) -> Vec<Region> {
        let dx = 2.0 * std::f64::consts::PI / n as f64;
        let mut seen = vec![false; field.len()];
        let mut out = Vec::new();
        for start in 0..field.len() {
            if seen[start] || field[start] <= threshold {
                continue;
            }
            let mut stack = vec![start];
            seen[start] = true;
            let (mut cnt, mut mx) = (0usize, 0.0f64);
            // centroid on the torus: average of unit-circle embeddings
            let mut cs = [[0.0f64; 2]; 3];
            while let Some(p) = stack.pop() {
                cnt += 1;
                mx = mx.max(field[p]);
                let ijl = [p / (n * n), (p / n) % n, p % n];
                for d in 0..3 {
                    let th = ijl[d] as f64 * dx;
                    cs[d][0] += th.cos();
                    cs[d][1] += th.sin();
                }
                for d in 0..3 {
                    for s in [1usize, n - 1] {
                        let mut q = ijl;
                        q[d] = (q[d] + s) % n;
                        let qi = (q[0] * n + q[1]) * n + q[2];
                        if !seen[qi] && field[qi] > threshold {
                            seen[qi] = true;
                            stack.push(qi);
                        }
                    }
                }
            }
            let centroid = [0, 1, 2].map(|d| cs[d][1].atan2(cs[d][0]).rem_euclid(2.0 * std::f64::consts::PI));
            out.push(Region { points: cnt, centroid, max: mx });
        }
        out.sort_by(|a, b| b.points.cmp(&a.points));
        out
    }
}
