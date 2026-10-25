//! Fluid–solid coupling of the compressible solver, validated against exact
//! 1D layered acoustics (`phononic`).
//!
//! # What exists in the crate
//!
//! * Solids are the GPR unified-continuum model (`compressible::gpr::Gpr`): a
//!   single-material `System<14>` with one NASG equation of state and one
//!   shear speed `c_s`. Its longitudinal linear speed is
//!   `c_L = sqrt(c_0^2 + 4 c_s^2/3)`, with `c_0` the bulk (EOS) speed.
//! * Multi-material flow is the Kapila five-equation model
//!   (`compressible::multiphase::FiveEq`): a `System<7>` carrying two partial
//!   densities and a volume fraction, one NASG EOS per phase, with a
//!   *diffuse* interface and Wood's mixture sound speed.
//!
//! The crate has **no shear-bearing multi-material system**. `threephase.rs`
//! says it explicitly: "its elasticity (GPR distortion per phase) is a
//! separate rung". To test the solver's fluid–solid coupling at normal
//! incidence we therefore implement that missing rung here (this file only)
//! as `GprMix`: the GPR distortion field carried together with the two
//! partial densities and the volume fraction of the Kapila model. Its
//! per-material longitudinal speed is exactly `sqrt(c_0^2 + 4 c_s^2/3)`,
//! which is the speed handed to `phononic`.
//!
//! # Measurements
//!
//! 1. `plate`: a Gaussian-modulated acoustic pulse in water hits a steel
//!    plate. The reflected and transmitted pressure spectra at probes just
//!    outside the plate are compared with the exact normal-incidence
//!    transfer matrix of the stack [water spacer, steel, water spacer]
//!    (`phononic::layered`), reported against grid resolution.
//! 2. `laminate`: 8 steel/aluminium quarter-wave periods. The transmitted
//!    spectrum is compared with (a) the exact finite stack and (b) the
//!    Reuss-homogenised effective layer. Below the first band gap the folded
//!    medium matches the exact one; inside the gap it does not. The
//!    `phononic::homogenize::homogenisation_error` indicator is quoted.
//!
//! Material constants come from `phononic::materials` (sourced there). Each
//! NASG is fitted to reproduce `(rho, c)` at `P0` exactly, because the
//! comparison is a linear-acoustics one; `gamma` only controls the (tiny)
//! nonlinearity.
//!
//! The 2D oblique case is **not** done: the only supported shear-bearing
//! system is single-material, so an oblique fluid–solid test would need the
//! same multi-material rung extended to 2D at a cost outside this task.

#![allow(dead_code)]

use compressible::eos::Nasg;
use compressible::multiphase::FiveEq;
use compressible::system::{Mat, System};
use compressible::{Bc, FluxMode, Mesh, Solver};
use phononic::homogenize;
use phononic::layered::{quarter_wave_band_edges, Layer, Stack};
use phononic::materials;

pub const NG: usize = 16;
type V = [f64; NG];
type M = Mat<NG>;

// --- state layout -------------------------------------------------------
// conserved: (a0 rho0, a1 rho1, rho u, rho v, rho w, A_xx..A_zz, E, alpha)
// primitive:  (a0 rho0, a1 rho1, u, v, w, A_xx..A_zz, p, alpha)
const ID1: usize = 0;
const ID2: usize = 1;
const IU: usize = 2;
const IA: usize = 5;
const IP: usize = 14;
const IAL: usize = 15;

#[inline]
fn amat(w: &V) -> [[f64; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|k| w[IA + 3 * i + k]))
}

#[inline]
fn gram(a: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    std::array::from_fn(|k| std::array::from_fn(|l| (0..3).map(|i| a[i][k] * a[i][l]).sum()))
}

#[inline]
fn dev(g: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let tr = (g[0][0] + g[1][1] + g[2][2]) / 3.0;
    std::array::from_fn(|k| std::array::from_fn(|l| g[k][l] - if k == l { tr } else { 0.0 }))
}

#[inline]
fn det3(a: &[[f64; 3]; 3]) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

/// Two-material GPR: the GPR distortion carried with the Kapila partial
/// densities and volume fraction. Material 0 is the fluid (`cs = 0`),
/// material 1 the solid; `alpha` is the volume fraction of material 1.
#[derive(Clone, Copy, Debug)]
pub struct GprMix {
    pub eos: [Nasg; 2],
    /// rest densities of the two materials
    pub rho: [f64; 2],
    /// shear sound speeds c_s of the two materials
    pub cs: [f64; 2],
    pub kapila: bool,
}

impl GprMix {
    #[inline]
    fn alphas(w: &V) -> [f64; 2] {
        [1.0 - w[IAL], w[IAL]]
    }
    /// Mixture density.
    #[inline]
    pub fn rho_mix(w: &V) -> f64 {
        w[ID1] + w[ID2]
    }
    /// Effective shear-speed-squared, Voigt average of the constituents.
    #[inline]
    fn cs2(&self, w: &V) -> f64 {
        let a = w[IAL];
        (1.0 - a) * self.cs[0] * self.cs[0] + a * self.cs[1] * self.cs[1]
    }

    /// D, B of the linear internal-energy form rho e = p D + B.
    #[inline]
    fn db(&self, w: &V) -> (f64, f64) {
        let al = Self::alphas(w);
        let (mut d, mut b) = (0.0, 0.0);
        for k in 0..2 {
            let e = &self.eos[k];
            let vol = al[k] - e.b * w[k];
            d += vol / (e.gamma - 1.0);
            b += e.gamma * e.p_inf * vol / (e.gamma - 1.0) + e.q * w[k];
        }
        (d, b)
    }

    #[inline]
    fn rhoe_thermal(&self, w: &V) -> f64 {
        let (d, b) = self.db(w);
        w[IP] * d + b
    }

    /// GPR elastic energy per unit mass, `(c_s^2/4)|dev G|^2`.
    #[inline]
    fn e2(&self, w: &V) -> f64 {
        let dg = dev(&gram(&amat(w)));
        let s: f64 = dg.iter().flatten().map(|x| x * x).sum();
        0.25 * self.cs2(w) * s
    }

    /// Cauchy stress `sigma = -rho c_s^2 G dev G`.
    pub fn stress(&self, w: &V) -> [[f64; 3]; 3] {
        let a = amat(w);
        let g = gram(&a);
        let dg = dev(&g);
        let f = -Self::rho_mix(w) * self.cs2(w);
        std::array::from_fn(|i| std::array::from_fn(|k| f * (0..3).map(|m| g[i][m] * dg[m][k]).sum::<f64>()))
    }

    #[inline]
    fn phase_rho_c2(&self, w: &V, k: usize) -> f64 {
        let al = Self::alphas(w);
        let e = &self.eos[k];
        let rho_k = if al[k] > 1e-12 && w[k] > 0.0 { w[k] / al[k] } else { 0.0 };
        e.gamma * (w[IP] + e.p_inf) / (1.0 - e.b * rho_k).max(1e-12)
    }

    fn kapila_k(&self, w: &V) -> f64 {
        if !self.kapila {
            return 0.0;
        }
        let al = [w[IAL].clamp(0.0, 1.0), 1.0 - w[IAL].clamp(0.0, 1.0)];
        let (z1, z2) = (self.phase_rho_c2(w, 0), self.phase_rho_c2(w, 1));
        let den = al[1] * z1 + al[0] * z2;
        if den <= 0.0 {
            return 0.0;
        }
        al[0] * al[1] * (z2 - z1) / den
    }

    /// Wood mixture `rho c_0^2` of the two bulk EOS (with the Kapila term).
    fn rho_c2_bulk(&self, w: &V) -> f64 {
        let (d, _) = self.db(w);
        let p = w[IP];
        let mut num = self.rhoe_thermal(w) + p;
        for k in 0..2 {
            let e = &self.eos[k];
            num += w[k] * (e.b * (p + e.gamma * e.p_inf) / (e.gamma - 1.0) - e.q);
        }
        num += self.kapila_k(w) * self.drhoe_dalpha(p);
        num / d
    }

    #[inline]
    fn drhoe_dalpha(&self, p: f64) -> f64 {
        let (e1, e2) = (&self.eos[0], &self.eos[1]);
        (p + e1.gamma * e1.p_inf) / (e1.gamma - 1.0) - (p + e2.gamma * e2.p_inf) / (e2.gamma - 1.0)
    }

    fn p_inf_mix(&self, w: &V) -> f64 {
        let a = w[IAL].clamp(0.0, 1.0);
        a * self.eos[1].p_inf + (1.0 - a) * self.eos[0].p_inf
    }

    /// Bulk (EOS) longitudinal speed of the mixture, without shear.
    pub fn c0(&self, w: &V) -> f64 {
        (self.rho_c2_bulk(w).max(0.0) / Self::rho_mix(w)).sqrt()
    }

    /// Longitudinal linear speed of pure material `m` at rest: the number
    /// handed to `phononic`.
    pub fn material_c_l(&self, m: usize) -> f64 {
        let e = &self.eos[m];
        let c0sq = e.gamma * (P0 + e.p_inf) / self.rho[m];
        (c0sq + 4.0 / 3.0 * self.cs[m] * self.cs[m]).sqrt()
    }
}

impl System<NG> for GprMix {
    fn to_prim(&self, u: &V) -> V {
        let rho = u[ID1] + u[ID2];
        let v = [u[IU] / rho, u[IU + 1] / rho, u[IU + 2] / rho];
        let mut w = [0.0; NG];
        w[ID1] = u[ID1];
        w[ID2] = u[ID2];
        w[IU..IU + 3].copy_from_slice(&v);
        w[IA..IA + 9].copy_from_slice(&u[IA..IA + 9]);
        w[IAL] = u[IAL];
        let ke = 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        let rhoe = u[IP] - rho * ke - rho * self.e2(&w);
        let (d, b) = self.db(&w);
        w[IP] = (rhoe - b) / d;
        w
    }

    #[inline]
    fn velocity_conserved(&self, u: &V, d: usize) -> f64 {
        u[IU + d] / (u[ID1] + u[ID2])
    }

    fn to_cons(&self, w: &V) -> V {
        let rho = w[ID1] + w[ID2];
        let v = [w[IU], w[IU + 1], w[IU + 2]];
        let mut u = [0.0; NG];
        u[ID1] = w[ID1];
        u[ID2] = w[ID2];
        for i in 0..3 {
            u[IU + i] = rho * v[i];
        }
        u[IA..IA + 9].copy_from_slice(&w[IA..IA + 9]);
        u[IAL] = w[IAL];
        let ke = 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        u[IP] = self.rhoe_thermal(w) + rho * self.e2(w) + rho * ke;
        u
    }

    /// Longitudinal speed: bulk Wood speed plus the GPR shear stiffening.
    fn sound_speed(&self, w: &V) -> f64 {
        let bulk = self.rho_c2_bulk(w).max(0.0) / Self::rho_mix(w);
        (bulk + 4.0 / 3.0 * self.cs2(w)).max(0.0).sqrt()
    }

    fn vel_index(&self) -> usize {
        IU
    }
    fn momentum_index(&self, d: usize) -> usize {
        IU + d
    }
    fn energy_index(&self) -> Option<usize> {
        Some(IP)
    }
    fn density(&self, w: &V) -> f64 {
        Self::rho_mix(w)
    }
    fn density_indices(&self) -> std::ops::Range<usize> {
        0..2
    }
    fn temperature(&self, w: &V) -> f64 {
        let al = Self::alphas(w);
        let mut t = 0.0;
        let mut wsum = 0.0;
        for k in 0..2 {
            let rho_k = if al[k] > 1e-12 && w[k] > 0.0 { w[k] / al[k] } else { 0.0 };
            if rho_k > 0.0 {
                t += w[k] * self.eos[k].temperature(rho_k, w[IP]);
                wsum += w[k];
            }
        }
        if wsum > 0.0 { t / wsum } else { f64::NAN }
    }

    fn eigen(&self, _w: &V, _d: usize, l: &mut M, r: &mut M) {
        for i in 0..NG {
            l[i] = [0.0; NG];
            r[i] = [0.0; NG];
            l[i][i] = 1.0;
            r[i][i] = 1.0;
        }
    }

    fn eigenvalues(&self, w: &V, d: usize) -> V {
        [w[IU + d].abs() + self.sound_speed(w); NG]
    }

    /// Rusanov flux (as the single-material GPR system uses).
    fn riemann_flux(&self, wl: &V, wr: &V, d: usize) -> (V, f64) {
        let s = self.max_speed(wl, d).max(self.max_speed(wr, d));
        let (fl, fr) = (self.flux(wl, d), self.flux(wr, d));
        let (ql, qr) = (self.to_cons(wl), self.to_cons(wr));
        let mut f = [0.0; NG];
        for k in 0..NG {
            f[k] = 0.5 * (fl[k] + fr[k]) - 0.5 * s * (qr[k] - ql[k]);
        }
        (f, 0.5 * (wl[IU + d] + wr[IU + d]))
    }

    fn split_flux(&self, a: &V, b: &V, d: usize) -> V {
        let un = IU + d;
        let rho = 0.5 * (Self::rho_mix(a) + Self::rho_mix(b));
        let u = 0.5 * (a[un] + b[un]);
        let c = rho * u;
        let (sa, sb) = (self.stress(a), self.stress(b));
        let mut f = [0.0; NG];
        f[ID1] = 0.5 * (a[ID1] + b[ID1]) * u;
        f[ID2] = 0.5 * (a[ID2] + b[ID2]) * u;
        for i in 0..3 {
            f[IU + i] = c * 0.5 * (a[IU + i] + b[IU + i]) - 0.5 * (sa[i][d] + sb[i][d]);
        }
        f[un] += 0.5 * (a[IP] + b[IP]);
        for i in 0..3 {
            let mut s = 0.0;
            for m in 0..3 {
                s += 0.25 * (a[IA + 3 * i + m] + b[IA + 3 * i + m]) * (a[IU + m] + b[IU + m]);
            }
            f[IA + 3 * i + d] = s;
        }
        f[IAL] = 0.5 * (a[IAL] + b[IAL]) * u;
        let kin = c * 0.5 * (a[IU] * b[IU] + a[IU + 1] * b[IU + 1] + a[IU + 2] * b[IU + 2]);
        let rhoe = 0.5 * (self.rhoe_thermal(a) + self.rhoe_thermal(b));
        let e2 = 0.5 * (Self::rho_mix(a) * self.e2(a) + Self::rho_mix(b) * self.e2(b));
        let pw = 0.5 * (a[un] * b[IP] + b[un] * a[IP]);
        let sw: f64 = (0..3).map(|i| 0.5 * (sa[i][d] * b[IU + i] + sb[i][d] * a[IU + i])).sum();
        f[IP] = kin + u * (rhoe + e2) + pw - sw;
        f
    }

    fn flux(&self, w: &V, d: usize) -> V {
        let un = IU + d;
        let q = self.to_cons(w);
        let s = self.stress(w);
        let v = [w[IU], w[IU + 1], w[IU + 2]];
        let mut f = [0.0; NG];
        f[ID1] = q[ID1] * v[d];
        f[ID2] = q[ID2] * v[d];
        for i in 0..3 {
            f[IU + i] = q[IU + i] * v[d] - s[i][d];
        }
        f[un] += w[IP];
        for i in 0..3 {
            f[IA + 3 * i + d] = (0..3).map(|m| w[IA + 3 * i + m] * v[m]).sum();
        }
        f[IAL] = w[IAL] * v[d];
        f[IP] = (q[IP] + w[IP]) * v[d] - (0..3).map(|i| s[i][d] * v[i]).sum::<f64>();
        f
    }

    fn dudw(&self, w: &V) -> M {
        let mut j = [[0.0; NG]; NG];
        let u0 = self.to_cons(w);
        for c in 0..NG {
            let h = 1e-7 * w[c].abs().max(1e-7);
            let mut wp = *w;
            wp[c] += h;
            let up = self.to_cons(&wp);
            for r in 0..NG {
                j[r][c] = (up[r] - u0[r]) / h;
            }
        }
        j
    }

    fn admissible(&self, w: &V) -> bool {
        let (d, _) = self.db(w);
        w.iter().all(|x| x.is_finite())
            && Self::rho_mix(w) > 0.0
            && d > 0.0
            && w[IP] + self.p_inf_mix(w) > 0.0
            && det3(&amat(w)) > 0.0
            && (-0.05..=1.05).contains(&w[IAL])
    }

    fn nonconservative(&self) -> Option<usize> {
        Some(IAL)
    }

    fn nc_coef(&self, w: &V) -> f64 {
        w[IAL] + self.kapila_k(w)
    }

    fn nc_general(&self) -> bool {
        true
    }

    fn nc_rhs(&self, w: &V, dw: &V, d: usize) -> V {
        let mut r = [0.0; NG];
        for i in 0..3 {
            let s: f64 = (0..3).map(|j| w[IU + j] * dw[IA + 3 * i + j]).sum();
            for k in 0..3 {
                let mut t = w[IU + d] * dw[IA + 3 * i + k];
                if k == d {
                    t -= s;
                }
                r[IA + 3 * i + k] = -t;
            }
        }
        r
    }

    fn sensor(&self, w: &V, out: &mut [f64; 8]) -> usize {
        out[0] = Self::rho_mix(w).max(1e-300).ln();
        out[1] = (w[IP] + self.p_inf_mix(w)).max(1e-300).ln();
        out[2] = w[IAL];
        3
    }
}

// ---------------------------------------------------------------------
// Material construction
// ---------------------------------------------------------------------

pub const P0: f64 = 101_325.0; // 1 atm
pub const GAMMA: f64 = 4.4; // stiffened-gas exponent (nonlinearity only)

/// Stiffened-gas NASG fitted to `(rho, c)` at pressure `p0` exactly:
/// `p_inf = rho c^2 / gamma - p0`.
pub fn nasg_at(rho: f64, c: f64, p0: f64, gamma: f64) -> Nasg {
    Nasg::stiffened(gamma, rho * c * c / gamma - p0)
}

/// Water as the fluid constituent: `c_s = 0`.
pub fn water() -> (f64, Nasg, f64) {
    let f = materials::WATER;
    (f.rho, nasg_at(f.rho, f.c, P0, GAMMA), 0.0)
}

/// A solid constituent. The NASG bulk speed is set so the GPR longitudinal
/// speed equals the handbook `c_p`: `c_0 = sqrt(c_p^2 - 4 c_s^2/3)`.
pub fn solid(s: &phononic::materials::Solid) -> (f64, Nasg, f64) {
    let c0 = (s.c_p * s.c_p - 4.0 / 3.0 * s.c_s * s.c_s).max(1.0).sqrt();
    (s.rho, nasg_at(s.rho, c0, P0, GAMMA), s.c_s)
}

/// Build a two-material system, fluid first, one solid.
pub fn mix(fluid: (f64, Nasg, f64), sol: (f64, Nasg, f64)) -> GprMix {
    GprMix {
        eos: [fluid.1, sol.1],
        rho: [fluid.0, sol.0],
        cs: [fluid.2, sol.2],
        kapila: false,
    }
}

// ---------------------------------------------------------------------
// Run configuration
// ---------------------------------------------------------------------

/// A 1D run: a pulse in the fluid hits a stack of material layers.
#[derive(Clone, Debug)]
pub struct RunConfig {
    pub n: usize,
    pub lo: f64,
    pub hi: f64,
    /// left and right interface of the test region
    pub x1: f64,
    pub x2: f64,
    /// layers inside `[x1, x2]` as (material index, thickness)
    pub layers: Vec<(usize, f64)>,
    pub xc: f64,
    pub width: f64,
    pub f0: f64,
    pub amp: f64,
    pub t_end: f64,
    pub cfl: f64,
    pub tol: f64,
}

impl RunConfig {
    pub fn dx(&self) -> f64 {
        (self.hi - self.lo) / self.n as f64
    }
    /// First/last plate cells (centres inside `[x1, x2]`), and the fluid
    /// spacer lengths in front of and behind.
    pub fn cell_span(&self) -> (usize, usize, f64, f64) {
        let dx = self.dx();
        let i1 = (0..self.n).find(|&i| self.lo + (i as f64 + 0.5) * dx >= self.x1).unwrap();
        let i2 = (0..self.n).rev().find(|&i| self.lo + (i as f64 + 0.5) * dx < self.x2).unwrap();
        let d_front = (self.x1 - (self.lo + (i1 as f64 - 0.5) * dx)).max(0.0);
        let d_back = ((self.lo + (i2 as f64 + 1.5) * dx) - self.x2).max(0.0);
        (i1, i2, d_front, d_back)
    }
    fn material_at(&self, x: f64) -> usize {
        if x < self.x1 || x >= self.x2 {
            return 0;
        }
        let mut s = self.x1;
        for &(m, d) in &self.layers {
            if x < s + d {
                return m;
            }
            s += d;
        }
        0
    }
}

/// Time series at the two probes plus mesh/stack geometry.
pub struct RunResult {
    pub t: Vec<f64>,
    pub pl: Vec<f64>,
    pub pr: Vec<f64>,
    pub dt: f64,
    pub i1: usize,
    pub i2: usize,
    pub d_front: f64,
    pub d_back: f64,
    pub steps: usize,
}

/// Smooth volume fraction of material `mat` in the stack: a sum of tanh
/// transitions over each layer of that material. Initialising the exact
/// diffuse profile avoids the spurious transient a sharp jump produces.
fn smooth_fraction(cfg: &RunConfig, x: f64, mat: usize) -> f64 {
    let delta = 1.5 * cfg.dx();
    let mut alpha = 0.0;
    let mut pos = cfg.x1;
    for &(m, d) in &cfg.layers {
        if m == mat {
            alpha += 0.5 * ((x - pos) / delta).tanh() - 0.5 * ((x - (pos + d)) / delta).tanh();
        }
        pos += d;
    }
    alpha.clamp(0.0, 1.0)
}

/// Primitive state of the fluid at rest plus a right-going pressure pulse.
fn fluid_at(sys: &GprMix, x: f64, cfg: &RunConfig, amp: f64) -> V {
    let rho = sys.rho[0];
    let c = sys.material_c_l(0);
    let mut w = [0.0; NG];
    w[ID1] = rho;
    w[IA] = 1.0;
    w[IA + 4] = 1.0;
    w[IA + 8] = 1.0;
    let r = x - cfg.xc;
    let env = (-(r / cfg.width).powi(2)).exp();
    let dp = amp * env * (2.0 * std::f64::consts::PI * cfg.f0 * r / c).sin();
    w[IP] = P0 + dp;
    w[IU] = dp / (rho * c);
    w
}

/// Initial state for the run: a smooth interface (if `with_layers`) plus the
/// pulse in the fluid region.
fn init_at(sys: &GprMix, cfg: &RunConfig, with_layers: bool, x: f64) -> V {
    let a = if with_layers { smooth_fraction(cfg, x, 1) } else { 0.0 };
    let mut w = [0.0; NG];
    w[ID1] = (1.0 - a) * sys.rho[0];
    w[ID2] = a * sys.rho[1];
    w[IAL] = a;
    w[IA] = 1.0;
    w[IA + 4] = 1.0;
    w[IA + 8] = 1.0;
    let rho = sys.rho[0];
    let c = sys.material_c_l(0);
    let r = x - cfg.xc;
    let env = (-(r / cfg.width).powi(2)).exp();
    let dp = cfg.amp * env * (2.0 * std::f64::consts::PI * cfg.f0 * r / c).sin();
    w[IP] = P0 + dp;
    w[IU] = dp / (rho * c);
    w
}

/// Build a solver for `cfg` and run to `t_end`, recording the two probes.
pub fn run(sys: &GprMix, cfg: &RunConfig, with_layers: bool, dt_override: Option<f64>) -> RunResult {
    let (i1, i2, d_front, d_back) = cfg.cell_span();
    let sys_local = *sys;
    let cfg_local = cfg.clone();
    let mesh = Mesh::new([cfg.n, 1, 1], [cfg.lo, 0.0, 0.0], [cfg.hi, 1.0, 1.0]);
    let mut s = Solver::new(sys_local, mesh, move |x| init_at(&sys_local, &cfg_local, with_layers, x[0]));
    let far = sys_local.to_cons(&fluid_at(&sys_local, -1e9, cfg, 0.0));
    s.bc[0] = [Bc::Characteristic(far), Bc::Characteristic(far)];
    s.cfl = cfg.cfl;
    s.tol = cfg.tol;
    s.mode = match std::env::var("FSI_MODE").as_deref() {
        Ok("central") => FluxMode::Central,
        Ok("weno") => FluxMode::Weno,
        _ => FluxMode::Hybrid,
    };

    let pl_idx = i1 - 1;
    let pr_idx = i2 + 1;
    let dt = dt_override.unwrap_or_else(|| s.dt_cfl());
    let mut t = Vec::new();
    let mut pl = Vec::new();
    let mut pr = Vec::new();
    let mut steps = 0usize;
    while s.t < cfg.t_end - 1e-12 {
        let d = dt.min(cfg.t_end - s.t);
        s.step(d);
        t.push(s.t);
        pl.push(s.prim(pl_idx, 0, 0)[IP]);
        pr.push(s.prim(pr_idx, 0, 0)[IP]);
        steps += 1;
    }
    RunResult { t, pl, pr, dt, i1, i2, d_front, d_back, steps }
}

// ---------------------------------------------------------------------
// Spectra
// ---------------------------------------------------------------------

/// Iterative radix-2 FFT (in place); `n` a power of two.
pub fn fft(re: &mut [f64], im: &mut [f64]) {
    let n = re.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            re.swap(i, j);
            im.swap(i, j);
        }
    }
    let mut len = 2;
    while len <= n {
        let ang = -2.0 * std::f64::consts::PI / len as f64;
        let (wr, wi) = (ang.cos(), ang.sin());
        let mut i = 0;
        while i < n {
            let (mut wlr, mut wli) = (1.0f64, 0.0f64);
            for k in 0..len / 2 {
                let ur = re[i + k];
                let ui = im[i + k];
                let vr = re[i + k + len / 2] * wlr - im[i + k + len / 2] * wli;
                let vi = re[i + k + len / 2] * wli + im[i + k + len / 2] * wlr;
                re[i + k] = ur + vr;
                im[i + k] = ui + vi;
                re[i + k + len / 2] = ur - vr;
                im[i + k + len / 2] = ui - vi;
                let nwr = wlr * wr - wli * wi;
                wli = wlr * wi + wli * wr;
                wlr = nwr;
            }
            i += len;
        }
        len <<= 1;
    }
}

/// Linear resample of `(ts, vals)` onto `nu` uniform points on `[0, t_end]`.
pub fn resample(ts: &[f64], vals: &[f64], nu: usize, t_end: f64) -> Vec<f64> {
    let dtu = t_end / nu as f64;
    let mut out = vec![0.0; nu];
    let mut j = 0usize;
    for m in 0..nu {
        let t = m as f64 * dtu;
        while j + 1 < ts.len() && ts[j + 1] < t {
            j += 1;
        }
        if ts.is_empty() {
            out[m] = 0.0;
        } else if j + 1 >= ts.len() {
            out[m] = *vals.last().unwrap();
        } else {
            let f = (t - ts[j]) / (ts[j + 1] - ts[j]);
            out[m] = vals[j] * (1.0 - f) + vals[j + 1] * f;
        }
    }
    out
}

/// Complex FFT of a real signal.
pub fn spectrum(x: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut re = x.to_vec();
    let mut im = vec![0.0; x.len()];
    fft(&mut re, &mut im);
    (re, im)
}

/// Small complex pair used by the public spectrum API.
#[derive(Clone, Copy, Debug)]
pub struct Cx {
    pub re: f64,
    pub im: f64,
}
impl Cx {
    pub fn new(re: f64, im: f64) -> Self {
        Cx { re, im }
    }
    pub fn abs(self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }
    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }
    pub fn sub(self, o: Cx) -> Cx {
        Cx::new(self.re - o.re, self.im - o.im)
    }
    pub fn div(self, o: Cx) -> Cx {
        let d = o.re * o.re + o.im * o.im;
        Cx::new((self.re * o.re + self.im * o.im) / d, (self.im * o.re - self.re * o.im) / d)
    }
}

/// Exact `(r, t)` of a layer list between half-spaces of impedance `z0`
/// (incident) and `zs` (exit).
pub fn phononic_rt_z(z0: f64, layers: &[(f64, f64, f64)], zs: f64, f: f64) -> (Cx, Cx) {
    let ls: Vec<Layer> = layers.iter().map(|&(r, c, d)| Layer::lossless(r, c, d)).collect();
    let st = Stack::new(z0, ls, zs);
    let sol = st.solve(2.0 * std::f64::consts::PI * f);
    (Cx::new(sol.r.re, sol.r.im), Cx::new(sol.t.re, sol.t.im))
}

/// Exact `(r, t)` of a layer list between water half-spaces.
pub fn phononic_rt(layers: &[(f64, f64, f64)], f: f64) -> (Cx, Cx) {
    let z = materials::WATER.rho * materials::WATER.c;
    phononic_rt_z(z, layers, z, f)
}

/// The physical layer list of a run, including the probe water spacers.
pub fn stack_of(cfg: &RunConfig, with_layers: bool) -> Vec<(usize, f64)> {
    let (_, _, d_front, d_back) = cfg.cell_span();
    let mut v = vec![(0usize, d_front)];
    if with_layers {
        v.extend(cfg.layers.iter().copied());
    }
    v.push((0, d_back));
    v
}

fn layer_rho_c(sys: &GprMix, m: usize, d: f64) -> (f64, f64, f64) {
    (sys.rho[m], sys.material_c_l(m), d)
}

// ---------------------------------------------------------------------
// Experiments
// ---------------------------------------------------------------------

const NU: usize = 16384;

/// Measured vs exact transfer functions over a frequency grid.
pub struct Spectra {
    pub f: Vec<f64>,
    pub inc: Vec<f64>,
    pub r_meas: Vec<Cx>,
    pub t_meas: Vec<Cx>,
    pub r_exact: Vec<Cx>,
    pub t_exact: Vec<Cx>,
    pub steps: usize,
    pub dt: f64,
}

/// Run the pulse experiment and build the transfer functions.
///
/// `exact_layers` is the reference stack used for `r_exact`/`t_exact` (either
/// the true layer list, or its homogenised replacement). The measured
/// spectra are always from `cfg` (which carries the true layers).
pub fn spectra(
    sys: &GprMix,
    cfg: &RunConfig,
    exact_layers: &[(f64, f64, f64)],
    f_lo: f64,
    f_hi: f64,
) -> Spectra {
    let with = run(sys, cfg, true, None);
    let without = run(sys, cfg, false, Some(with.dt));
    let t_end = cfg.t_end;
    let dtu = t_end / NU as f64;
    let inc = spectrum(&resample(&without.t, &without.pl, NU, t_end));
    let tot = spectrum(&resample(&with.t, &with.pl, NU, t_end));
    let tra = spectrum(&resample(&with.t, &with.pr, NU, t_end));
    let df = 1.0 / (NU as f64 * dtu);
    let k0 = (f_lo / df).round() as usize;
    let k1 = (f_hi / df).round() as usize;
    let mut out = Spectra {
        f: Vec::new(),
        inc: Vec::new(),
        r_meas: Vec::new(),
        t_meas: Vec::new(),
        r_exact: Vec::new(),
        t_exact: Vec::new(),
        steps: with.steps,
        dt: with.dt,
    };
    for k in k0..=k1 {
        let fr = k as f64 * df;
        let (ir, ii) = (inc.0[k], inc.1[k]);
        let inc_c = Cx::new(ir, ii);
        let mag = inc_c.abs();
        if mag < 1e-12 {
            continue;
        }
        let r = Cx::new(tot.0[k] - ir, tot.1[k] - ii).div(inc_c);
        let t = Cx::new(tra.0[k], tra.1[k]).div(inc_c);
        let (re, te) = phononic_rt(exact_layers, fr);
        out.f.push(fr);
        out.inc.push(mag);
        out.r_meas.push(r);
        out.t_meas.push(t);
        out.r_exact.push(re);
        out.t_exact.push(te);
    }
    out
}

/// Plate configuration: 20 mm steel plate at x = 0.25 m in water.
pub fn plate_config(n: usize) -> RunConfig {
    RunConfig {
        n,
        lo: 0.0,
        hi: 0.6,
        x1: 0.25,
        x2: 0.27,
        layers: vec![(1, 0.020)],
        xc: 0.10,
        width: 0.012,
        f0: 80_000.0,
        amp: 100.0,
        t_end: 2.6e-4,
        cfl: 0.4,
        tol: 1e-3,
    }
}

/// Laminate configuration: 8 steel/water quarter-wave periods at 400 kHz.
/// The two-phase model limits the constituents to water and one solid, so
/// the second layer is water rather than epoxy/aluminium.
/// Laminate configuration: 8 aluminium/steel quarter-wave periods at
/// 400 kHz, in a steel matrix (phase 1 is steel, phase 0 aluminium).
pub fn laminate_config(n: usize) -> RunConfig {
    let f0 = 400_000.0;
    let d_al = materials::ALUMINIUM.c_p / (4.0 * f0);
    let d_s = materials::STEEL.c_p / (4.0 * f0);
    let mut layers = Vec::new();
    for _ in 0..8 {
        layers.push((0usize, d_al));
        layers.push((1, d_s));
    }
    let d = 8.0 * (d_al + d_s);
    RunConfig {
        n,
        lo: 0.0,
        hi: 0.6,
        x1: 0.25,
        x2: 0.25 + d,
        layers,
        xc: 0.10,
        width: 0.002,
        f0: 350_000.0,
        amp: 100.0,
        t_end: 1.6e-4,
        cfl: 0.4,
        tol: 1e-3,
    }
}

/// Relative L2 error of the measured magnitude against the exact one,
/// weighted by the incident spectrum over the frequency band.
pub fn mag_error(meas: &[Cx], exact: &[Cx], inc: &[f64]) -> f64 {
    let (mut num, mut den) = (0.0, 0.0);
    for k in 0..meas.len() {
        let w = inc[k] * inc[k];
        num += w * (meas[k].abs() - exact[k].abs()).powi(2);
        den += w * exact[k].abs().powi(2);
    }
    (num / den).sqrt()
}

/// RMS absolute magnitude error, incident-weighted. For a plate in water the
/// transmitted amplitude is small (|T| ~ 0.06-0.1), so this is the natural
/// absolute measure of the interface error.
pub fn mag_abs_error(meas: &[Cx], exact: &[Cx], inc: &[f64]) -> f64 {
    let (mut num, mut den) = (0.0, 0.0);
    for k in 0..meas.len() {
        let w = inc[k] * inc[k];
        num += w * (meas[k].abs() - exact[k].abs()).powi(2);
        den += w;
    }
    (num / den.max(1e-300)).sqrt()
}

/// RMS phase error (radians) between measured and exact, magnitude-weighted.
pub fn phase_error(meas: &[Cx], exact: &[Cx], inc: &[f64]) -> f64 {
    let (mut num, mut den) = (0.0, 0.0);
    for k in 0..meas.len() {
        let w = inc[k] * inc[k] * exact[k].abs().powi(2);
        let mut d = meas[k].arg() - exact[k].arg();
        while d > std::f64::consts::PI {
            d -= 2.0 * std::f64::consts::PI;
        }
        while d < -std::f64::consts::PI {
            d += 2.0 * std::f64::consts::PI;
        }
        num += w * d * d;
        den += w;
    }
    (num / den.max(1e-300)).sqrt()
}

/// Reuss-homogenised layer `(rho_eff, c_eff, thickness)` of a periodic
/// stack, plus the period and speed used by the validity indicator.
pub fn homogenised(layers: &[(f64, f64, f64)]) -> (f64, f64, f64, f64, f64) {
    let total: f64 = layers.iter().map(|l| l.2).sum();
    let mut inv_m = 0.0;
    let mut rho = 0.0;
    for &(r, c, d) in layers {
        let f = d / total;
        inv_m += f / (r * c * c);
        rho += f * r;
    }
    let m = 1.0 / inv_m;
    let c_eff = (m / rho).sqrt();
    (rho, c_eff, total, c_eff, total)
}

// =====================================================================
// Normal-incidence fluid–solid validation with the crate's validated
// multi-material model (`FiveEq`).
//
// At normal incidence the shear degrees of freedom decouple: a solid's
// P-wave reflection is fixed by the constrained modulus M = K + 4G/3, i.e.
// by its density `rho` and longitudinal speed `c_L = sqrt(M/rho)`. A
// two-phase stiffened-gas fluid with the same `(rho, c_L)` therefore has
// exactly the same normal-incidence impedance and wave speed. `FiveEq` is
// the crate's validated multi-material model and is stable on this contrast
// (unlike the shear-bearing `GprMix` extension above). The solid's shear
// speed enters only through `c_L`, which is also what `phononic` receives.
// =====================================================================

/// Longitudinal (constrained-modulus) data of the two phases.
#[derive(Clone, Copy)]
pub struct Media {
    pub rho: [f64; 2],
    pub c: [f64; 2],
}

/// Two phases as stiffened-gas fluids, both at their constrained-modulus
/// longitudinal speed.
pub fn fiveeq_system_media(p0: (f64, f64), p1: (f64, f64)) -> (FiveEq, Media) {
    let e = [nasg_at(p0.0, p0.1, P0, GAMMA), nasg_at(p1.0, p1.1, P0, GAMMA)];
    (
        FiveEq { eos: e, kapila: false },
        Media { rho: [p0.0, p1.0], c: [p0.1, p1.1] },
    )
}

/// Water (phase 0) and steel (phase 1), both as stiffened-gas fluids at the
/// constrained-modulus speed.
pub fn fiveeq_system() -> (FiveEq, Media) {
    fiveeq_system_media(
        (materials::WATER.rho, materials::WATER.c),
        (materials::STEEL.rho, materials::STEEL.c_p),
    )
}

/// Aluminium (phase 0) and steel (phase 1) for the laminate.
pub fn fiveeq_system_al_steel() -> (FiveEq, Media) {
    fiveeq_system_media(
        (materials::ALUMINIUM.rho, materials::ALUMINIUM.c_p),
        (materials::STEEL.rho, materials::STEEL.c_p),
    )
}

impl RunConfig {
    /// Phononic layer list for the two-phase media model.
    pub fn stack_media(&self, media: &Media, bg: usize, with_layers: bool) -> Vec<(f64, f64, f64)> {
        let (_, _, df, db) = self.cell_span();
        let mut v = vec![(media.rho[bg], media.c[bg], df)];
        if with_layers {
            for &(m, d) in &self.layers {
                v.push((media.rho[m], media.c[m], d));
            }
        }
        v.push((media.rho[bg], media.c[bg], db));
        v
    }
}

/// Build and run a `FiveEq` experiment. `bg` is the background phase (0 or
/// 1); the pulse is launched in it and the other phase forms the inclusions.
pub fn run5(
    sys: &FiveEq,
    media: &Media,
    bg: usize,
    cfg: &RunConfig,
    with_layers: bool,
    dt_override: Option<f64>,
) -> RunResult {
    let (i1, i2, d_front, d_back) = cfg.cell_span();
    let sy = *sys;
    let cf = cfg.clone();
    let inc = 1 - bg;
    let (rho_bg, c_bg) = (media.rho[bg], media.c[bg]);
    let mesh = Mesh::new([cfg.n, 1, 1], [cfg.lo, 0.0, 0.0], [cfg.hi, 1.0, 1.0]);
    let build = move |x: f64, amp: f64| {
        let a_inc = if with_layers { smooth_fraction(&cf, x, inc) } else { 0.0 };
        let a_bg = 1.0 - a_inc;
        let a0 = if bg == 0 { a_bg } else { a_inc };
        let r = x - cf.xc;
        let env = (-(r / cf.width).powi(2)).exp();
        let dp = amp * env * (2.0 * std::f64::consts::PI * cf.f0 * r / c_bg).sin();
        let mut w = [0.0; 7];
        w[0] = a0 * media.rho[0];
        w[1] = (1.0 - a0) * media.rho[1];
        w[2] = dp / (rho_bg * c_bg);
        w[5] = P0 + dp;
        w[6] = a0;
        w
    };
    let far = sy.to_cons(&build(-1e9, 0.0));
    let amp = cfg.amp;
    let mut s = Solver::new(sy, mesh, move |x| build(x[0], amp));
    s.bc[0] = [Bc::Characteristic(far), Bc::Characteristic(far)];
    s.cfl = cfg.cfl;
    let pl_idx = i1 - 1;
    let pr_idx = i2 + 1;
    let dt = dt_override.unwrap_or_else(|| s.dt_cfl());
    let mut t = Vec::new();
    let mut pl = Vec::new();
    let mut pr = Vec::new();
    let mut steps = 0usize;
    while s.t < cfg.t_end - 1e-12 {
        let d = dt.min(cfg.t_end - s.t);
        s.step(d);
        t.push(s.t);
        pl.push(s.prim(pl_idx, 0, 0)[5]);
        pr.push(s.prim(pr_idx, 0, 0)[5]);
        steps += 1;
    }
    RunResult { t, pl, pr, dt, i1, i2, d_front, d_back, steps }
}

/// Measured transfer functions of a `FiveEq` run.
pub struct Meas {
    pub f: Vec<f64>,
    pub inc: Vec<f64>,
    pub r: Vec<Cx>,
    pub t: Vec<Cx>,
    pub steps: usize,
    pub dt: f64,
}

/// Run the pulse experiment once and build the transfer functions.
pub fn meas5(sys: &FiveEq, media: &Media, bg: usize, cfg: &RunConfig, f_lo: f64, f_hi: f64) -> Meas {
    let with = run5(sys, media, bg, cfg, true, None);
    let without = run5(sys, media, bg, cfg, false, Some(with.dt));
    let t_end = cfg.t_end;
    let dtu = t_end / NU as f64;
    let inc = spectrum(&resample(&without.t, &without.pl, NU, t_end));
    let tot = spectrum(&resample(&with.t, &with.pl, NU, t_end));
    let tra = spectrum(&resample(&with.t, &with.pr, NU, t_end));
    let df = 1.0 / (NU as f64 * dtu);
    let k0 = (f_lo / df).round() as usize;
    let k1 = (f_hi / df).round() as usize;
    let mut out = Meas { f: Vec::new(), inc: Vec::new(), r: Vec::new(), t: Vec::new(), steps: with.steps, dt: with.dt };
    for k in k0..=k1 {
        let fr = k as f64 * df;
        let inc_c = Cx::new(inc.0[k], inc.1[k]);
        let mag = inc_c.abs();
        if mag < 1e-12 {
            continue;
        }
        out.f.push(fr);
        out.inc.push(mag);
        out.r.push(Cx::new(tot.0[k] - inc.0[k], tot.1[k] - inc.1[k]).div(inc_c));
        out.t.push(Cx::new(tra.0[k], tra.1[k]).div(inc_c));
    }
    out
}

/// Exact `(r, t)` of a phononic layer list on the frequency grid of `meas`,
/// between half-spaces of impedance `z0` and `zs`.
pub fn exact_curve(f: &[f64], layers: &[(f64, f64, f64)], z0: f64, zs: f64) -> (Vec<Cx>, Vec<Cx>) {
    let mut r = Vec::with_capacity(f.len());
    let mut t = Vec::with_capacity(f.len());
    for &fr in f {
        let (rr, tt) = phononic_rt_z(z0, layers, zs, fr);
        r.push(rr);
        t.push(tt);
    }
    (r, t)
}

// ---------------------------------------------------------------------
// main / example output
// ---------------------------------------------------------------------

/// A longitudinal pressure pulse in a pure GPR solid; returns the measured
/// front speed.
pub fn solid_pulse_speed(sys: &GprMix, n: usize, t_end: f64) -> f64 {
    let (lo, hi) = (0.0, 3.0);
    let rho = sys.rho[1];
    let (xc, width, _f0) = (0.4, 0.02, 200_000.0);
    let mesh = Mesh::new([n, 1, 1], [lo, 0.0, 0.0], [hi, 1.0, 1.0]);
    let sy = *sys;
    let mut s = Solver::new(sy, mesh, move |x| {
        let mut w = [0.0; NG];
        w[ID2] = rho;
        w[IAL] = 1.0;
        w[IA] = 1.0;
        w[IA + 4] = 1.0;
        w[IA + 8] = 1.0;
        let r = x[0] - xc;
        let env = (-(r / width).powi(2)).exp();
        // pure velocity pulse, as in `gpr_validate waves` (clean eigenmode)
        w[IP] = P0;
        w[IU] = 10.0 * env;
        w
    });
    s.bc[0] = [Bc::Outflow; 2];
    s.cfl = 0.4;
    s.run_to(t_end);
    let dx = (hi - lo) / n as f64;
    let (mut best, mut ibest) = (0.0f64, 0usize);
    for i in n / 5..n {
        let v = s.prim(i, 0, 0)[IU].abs();
        if v > best {
            best = v;
            ibest = i;
        }
    }
    (lo + (ibest as f64 + 0.5) * dx - xc) / t_end
}

/// Water against a stiff fluid (cs = 0): peak reflection ratio against Fresnel.
pub fn fresnel_ratio(n: usize) -> (f64, f64) {
    let stiff = (
        materials::STEEL.rho,
        nasg_at(materials::STEEL.rho, materials::STEEL.c_p, P0, GAMMA),
        0.0,
    );
    fresnel_ratio_sys(&mix(water(), stiff), n)
}

/// Same, for an arbitrary second fluid.
pub fn fresnel_ratio_sys(sys: &GprMix, n: usize) -> (f64, f64) {
    let cfg = RunConfig {
        n,
        lo: 0.0,
        hi: 0.8,
        x1: 0.30,
        x2: 0.80,
        layers: vec![(1, 0.50)],
        xc: 0.12,
        width: 0.012,
        f0: 80_000.0,
        amp: 100.0,
        t_end: 3.0e-4,
        cfl: 0.4,
        tol: 1e-3,
    };
    let with = run(sys, &cfg, true, None);
    let without = run(sys, &cfg, false, Some(with.dt));
    let pk = |v: &[f64]| v.iter().cloned().fold(0.0f64, |m, x| m.max(x.abs()));
    let inc: Vec<f64> = without.pl.iter().map(|a| a - P0).collect();
    let refl: Vec<f64> = with.pl.iter().zip(&without.pl).map(|(a, b)| a - b).collect();
    let z1 = materials::WATER.rho * materials::WATER.c;
    let z2 = sys.rho[1] * sys.material_c_l(1);
    let r_f = (z2 - z1).abs() / (z2 + z1);
    println!("    [n={n}] pk(inc)={:.4}, pk(refl)={:.4}, pk(plate)={:.4}", pk(&inc), pk(&refl), pk(&with.pl.iter().map(|a| a - P0).collect::<Vec<_>>()));
    (pk(&refl) / pk(&inc), r_f)
}

/// Stationary ideal-gas contact: all p_inf = 0, equal pressure. Returns the
/// maximum relative pressure deviation after `steps` steps at rest.
pub fn contact_probe(n: usize, steps: usize) -> f64 {
    let gamma = 1.4;
    let p0 = 1.0;
    let e0 = nasg_at(1.0, (gamma * p0 / 1.0f64).sqrt(), p0, gamma);
    let e1 = nasg_at(4.0, (gamma * p0 / 4.0f64).sqrt(), p0, gamma);
    let sys = GprMix { eos: [e0, e1], rho: [1.0, 4.0], cs: [0.0, 0.0], kapila: true };
    let mesh = Mesh::new([n, 1, 1], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let mut s = Solver::new(sys, mesh, |x| {
        let a = 0.5 * (1.0 + ((x[0] - 0.5) / 0.02).tanh());
        let mut w = [0.0; NG];
        w[ID1] = (1.0 - a) * 1.0;
        w[ID2] = a * 4.0;
        w[IAL] = a;
        w[IA] = 1.0;
        w[IA + 4] = 1.0;
        w[IA + 8] = 1.0;
        w[IP] = p0;
        w
    });
    s.bc[0] = [Bc::Outflow; 2];
    s.cfl = 0.4;
    let dt = s.dt_cfl();
    for _ in 0..steps {
        s.step(dt);
    }
    let mut dev = 0.0f64;
    for i in n / 4..3 * n / 4 {
        let p = s.prim(i, 0, 0)[IP];
        dev = dev.max(((p - p0) / p0).abs());
    }
    dev
}

/// Stationary water / stiff-fluid contact at rest. Returns the max relative
/// pressure deviation and the max |u| after `steps` steps.
pub fn stiff_contact_probe(n: usize, steps: usize) -> (f64, f64) {
    let stiff = (
        materials::STEEL.rho,
        nasg_at(materials::STEEL.rho, materials::STEEL.c_p, P0, GAMMA),
        0.0,
    );
    let sys = mix(water(), stiff);
    let (rho0, rho1) = (sys.rho[0], sys.rho[1]);
    let mesh = Mesh::new([n, 1, 1], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let mut s = Solver::new(sys, mesh, |x| {
        let a = 0.5 * (1.0 + ((x[0] - 0.5) / 0.02).tanh());
        let mut w = [0.0; NG];
        w[ID1] = (1.0 - a) * rho0;
        w[ID2] = a * rho1;
        w[IAL] = a;
        w[IA] = 1.0;
        w[IA + 4] = 1.0;
        w[IA + 8] = 1.0;
        w[IP] = P0;
        w
    });
    s.bc[0] = [Bc::Outflow; 2];
    s.cfl = 0.4;
    let dt = s.dt_cfl();
    for st in 0..steps {
        s.step(dt);
        if st % 50 == 0 || st + 1 == steps {
            let (mut dp, mut u) = (0.0f64, 0.0f64);
            for i in n / 4..3 * n / 4 {
                let w = s.prim(i, 0, 0);
                dp = dp.max(((w[IP] - P0) / P0).abs());
                u = u.max(w[IU].abs());
            }
            println!("    step {st:4}: max|dp|/p {dp:.3e}, max|u| {u:.3e}, dt {dt:.2e}");
        }
    }
    let (mut dp, mut u) = (0.0f64, 0.0f64);
    for i in n / 4..3 * n / 4 {
        let w = s.prim(i, 0, 0);
        dp = dp.max(((w[IP] - P0) / P0).abs());
        u = u.max(w[IU].abs());
    }
    (dp, u)
}

/// The crate's existing multi-material model (`FiveEq`) on the water/steel
/// contrast, treated as two stiffened-gas fluids with the constrained
/// longitudinal modulus. Returns peak reflection ratio vs Fresnel, the
/// maximum pressure deviation after the run, and the peak transmitted ratio.
pub fn fiveeq_fresnel(n: usize) -> (f64, f64, f64) {
    use compressible::multiphase::FiveEq;
    let eos = [
        nasg_at(materials::WATER.rho, materials::WATER.c, P0, GAMMA),
        nasg_at(materials::STEEL.rho, materials::STEEL.c_p, P0, GAMMA),
    ];
    let sys = FiveEq { eos, kapila: false };
    let (rho_w, rho_s, c_w) = (materials::WATER.rho, materials::STEEL.rho, materials::WATER.c);
    let (lo, hi) = (0.0, 0.8);
    let x1 = 0.30;
    let mesh = Mesh::new([n, 1, 1], [lo, 0.0, 0.0], [hi, 1.0, 1.0]);
    let build = |x: f64, pulse: bool, iface: bool| {
        let a1 = if iface {
            1.0 - 0.5 * (1.0 + ((x - x1) / (1.5 * (hi - lo) / n as f64)).tanh())
        } else {
            1.0
        };
        let r = x - 0.12;
        let env = (-(r / 0.012).powi(2)).exp();
        let dp = if pulse { 100.0 * env * (2.0 * std::f64::consts::PI * 80_000.0 * r / c_w).sin() } else { 0.0 };
        let mut w = [0.0; 7];
        w[0] = a1 * rho_w;
        w[1] = (1.0 - a1) * rho_s;
        w[2] = dp / (rho_w * c_w);
        w[5] = P0 + dp;
        w[6] = a1;
        w
    };
    let run = |pulse: bool, iface: bool| -> (Vec<f64>, Vec<f64>, f64) {
        let mut s = Solver::new(sys, mesh.clone(), move |x| build(x[0], pulse, iface));
        let far = sys.to_cons(&build(-1e9, false, false));
        s.bc[0] = [Bc::Characteristic(far), Bc::Characteristic(far)];
        s.cfl = 0.4;
        let i1 = (0..n).find(|&i| lo + (i as f64 + 0.5) * (hi - lo) / n as f64 >= x1).unwrap();
        let dt = s.dt_cfl();
        let t_end = 3.0e-4;
        let mut pl = Vec::new();
        let mut pr = Vec::new();
        while s.t < t_end - 1e-12 {
            let d = dt.min(t_end - s.t);
            s.step(d);
            pl.push(s.prim(i1 - 1, 0, 0)[5]);
            pr.push(s.prim(i1 + 1, 0, 0)[5]);
        }
        (pl, pr, dt)
    };
    let (pl, _, _) = run(true, true);
    let (pl0, _, _) = run(true, false);
    let pk = |v: &[f64]| v.iter().cloned().fold(0.0f64, |m, x| m.max(x.abs()));
    let inc: Vec<f64> = pl0.iter().map(|a| a - P0).collect();
    let refl: Vec<f64> = pl.iter().zip(&pl0).map(|(a, b)| a - b).collect();
    let z1 = rho_w * c_w;
    let z2 = rho_s * materials::STEEL.c_p;
    let r_f = (z2 - z1).abs() / (z2 + z1);
    (pk(&refl) / pk(&inc), r_f, pk(&refl))
}

fn speed_check(sys: &GprMix) {
    println!("single-material longitudinal speeds:");
    for m in 0..2 {
        let mut w = [0.0; NG];
        if m == 0 {
            w[ID1] = sys.rho[0];
        } else {
            w[ID2] = sys.rho[m];
            w[IAL] = 1.0;
        }
        w[IA] = 1.0;
        w[IA + 4] = 1.0;
        w[IA + 8] = 1.0;
        w[IP] = P0;
        let c_impl = sys.sound_speed(&w);
        let c_analytic = sys.material_c_l(m);
        println!(
            "  material {m} (rho {:.1}): solver sqrt(c0^2+4cs^2/3) = {:.4} m/s, analytic {:.4} m/s ({:+.2e})",
            sys.rho[m],
            c_impl,
            c_analytic,
            c_impl / c_analytic - 1.0
        );
    }
}

fn ns_from_env(def: &str) -> Vec<usize> {
    std::env::var("FSI_N")
        .unwrap_or_else(|_| def.into())
        .split(',')
        .filter_map(|s| s.trim().parse().ok())
        .collect()
}

fn plate_report(sys: &GprMix) {
    println!("\nsingle plate (20 mm steel in water), (set FSI_N, default 256,512):");
    for n in ns_from_env("256,512") {
        let cfg = plate_config(n);
        let exact: Vec<(f64, f64, f64)> = stack_of(&cfg, true)
            .into_iter()
            .map(|(m, d)| layer_rho_c(sys, m, d))
            .collect();
        let s = spectra(sys, &cfg, &exact, 20_000.0, 250_000.0);
        let er = mag_error(&s.r_meas, &s.r_exact, &s.inc);
        let et = mag_error(&s.t_meas, &s.t_exact, &s.inc);
        let pr = phase_error(&s.r_meas, &s.r_exact, &s.inc);
        let pt = phase_error(&s.t_meas, &s.t_exact, &s.inc);
        println!(
            "  N={n:5} steps {:5}: |R| rel err {:.3e}, |T| rel err {:.3e}; phase rms {:.3e}/{:.3e} rad",
            s.steps, er, et, pr, pt
        );
        for target in [40_000.0f64, 80_000.0, 160_000.0] {
            if let Some(k) = nearest(&s.f, target) {
                println!(
                    "      f={:7.0} Hz  |R| {:.5} vs {:.5}  |T| {:.5} vs {:.5}",
                    s.f[k], s.r_meas[k].abs(), s.r_exact[k].abs(), s.t_meas[k].abs(), s.t_exact[k].abs()
                );
            }
        }
    }
}

fn laminate_report(sys: &GprMix) {
    println!("\n8 steel/Al periods (quarter-wave at 400 kHz):");
    let f0 = 400_000.0;
    let z_s = materials::STEEL.rho * materials::STEEL.c_p;
    let z_al = materials::ALUMINIUM.rho * materials::ALUMINIUM.c_p;
    let (gap_lo, gap_hi) = quarter_wave_band_edges(z_s, z_al, f0);
    println!("  exact first gap {:.0} .. {:.0} Hz", gap_lo, gap_hi);
    for n in [1024usize] {
        let cfg = laminate_config(n);
        let layers: Vec<(f64, f64, f64)> = cfg
            .layers
            .iter()
            .map(|&(m, d)| layer_rho_c(sys, m, d))
            .collect();
        let period = 8.0f64.recip() * layers.iter().map(|l| l.2).sum::<f64>();
        let (rho_e, c_e, d_e, _, _) = homogenised(&layers);
        let stack: Vec<(f64, f64, f64)> = stack_of(&cfg, true)
            .into_iter()
            .map(|(m, d)| layer_rho_c(sys, m, d))
            .collect();
        let mut homog = vec![stack[0]];
        homog.push((rho_e, c_e, d_e));
        homog.push(*stack.last().unwrap());
        let s = spectra(sys, &cfg, &stack, 60_000.0, 620_000.0);
        let sh = spectra(sys, &cfg, &homog, 60_000.0, 620_000.0);
        println!(
            "  N={n:5} steps {:5}: exact vs measured |T| rel err {:.3e}; homogenised vs measured {:.3e}; rho_eff {:.1}, c_eff {:.1} m/s",
            s.steps,
            mag_error(&s.t_meas, &s.t_exact, &s.inc),
            mag_error(&s.t_meas, &sh.t_exact, &s.inc),
            rho_e,
            c_e
        );
        println!("  {:>8} {:>9} {:>9} {:>9} {:>10}", "f(kHz)", "|T|meas", "|T|exact", "|T|homog", "indicator");
        for target in [100_000.0f64, 200_000.0, 350_000.0, 450_000.0, 550_000.0] {
            if let Some(k) = nearest(&s.f, target) {
                let ind = homogenize::homogenisation_error(2.0 * std::f64::consts::PI * s.f[k], period, c_e);
                println!(
                    "  {:>8.0} {:>9.5} {:>9.5} {:>9.5} {:>10.2e}",
                    s.f[k] / 1e3,
                    s.t_meas[k].abs(),
                    s.t_exact[k].abs(),
                    sh.t_exact[k].abs(),
                    ind
                );
            }
        }
    }
}

fn plate_report5() {
    let (sys, media) = fiveeq_system();
    let zw = materials::WATER.rho * materials::WATER.c;
    println!("\n[FiveEq] single plate (20 mm steel in water), (set FSI_N, try 256,512,1024,2048):");
    for n in ns_from_env("256,512") {
        let cfg = plate_config(n);
        let exact = cfg.stack_media(&media, 0, true);
        let m = meas5(&sys, &media, 0, &cfg, 20_000.0, 250_000.0);
        let (re, te) = exact_curve(&m.f, &exact, zw, zw);
        println!(
            "  N={n:5} steps {:5}: |R| rel err {:.3e}, |T| rel err {:.3e}; |T| abs err {:.3e}; phase rms {:.3e}/{:.3e} rad",
            m.steps,
            mag_error(&m.r, &re, &m.inc),
            mag_error(&m.t, &te, &m.inc),
            mag_abs_error(&m.t, &te, &m.inc),
            phase_error(&m.r, &re, &m.inc),
            phase_error(&m.t, &te, &m.inc),
        );
        for target in [40_000.0f64, 80_000.0, 160_000.0] {
            if let Some(k) = nearest(&m.f, target) {
                println!(
                    "      f={:7.0} Hz  |R| {:.5} vs {:.5}  |T| {:.5} vs {:.5}",
                    m.f[k], m.r[k].abs(), re[k].abs(), m.t[k].abs(), te[k].abs()
                );
            }
        }
    }
}

fn laminate_report5() {
    let (sys, media) = fiveeq_system_al_steel();
    let z_al = materials::ALUMINIUM.rho * materials::ALUMINIUM.c_p;
    let z_s = materials::STEEL.rho * materials::STEEL.c_p;
    println!("\n[FiveEq] 8 aluminium/steel periods in a steel matrix (quarter-wave at 400 kHz):");
    let f0 = 400_000.0;
    let (gap_lo, gap_hi) = quarter_wave_band_edges(z_s, z_al, f0);
    println!("  exact first gap {:.0} .. {:.0} Hz", gap_lo, gap_hi);
    for n in ns_from_env("1024") {
        let cfg = laminate_config(n);
        let layers: Vec<(f64, f64, f64)> = cfg
            .layers
            .iter()
            .map(|&(m, d)| (media.rho[m], media.c[m], d))
            .collect();
        let period = layers.iter().map(|l| l.2).sum::<f64>() / 8.0;
        let (rho_e, c_e, d_e, _, _) = homogenised(&layers);
        let exact = cfg.stack_media(&media, 1, true);
        let mut homog = vec![exact[0]];
        homog.push((rho_e, c_e, d_e));
        homog.push(*exact.last().unwrap());
        let m = meas5(&sys, &media, 1, &cfg, 100_000.0, 700_000.0);
        let (_re, te) = exact_curve(&m.f, &exact, z_s, z_s);
        let (_rh, th) = exact_curve(&m.f, &homog, z_s, z_s);
        println!(
            "  N={n:5} steps {:5}: measured vs exact |T| rel err {:.3e}; homogenised vs exact {:.3e}; rho_eff {:.1}, c_eff {:.1} m/s",
            m.steps,
            mag_error(&m.t, &te, &m.inc),
            mag_error(&th, &te, &m.inc),
            rho_e,
            c_e
        );
        println!("  {:>8} {:>9} {:>9} {:>9} {:>10}", "f(kHz)", "|T|meas", "|T|exact", "|T|homog", "indicator");
        for target in [150_000.0f64, 250_000.0, 350_000.0, 450_000.0, 650_000.0] {
            if let Some(k) = nearest(&m.f, target) {
                let ind = homogenize::homogenisation_error(2.0 * std::f64::consts::PI * m.f[k], period, c_e);
                println!(
                    "  {:>8.0} {:>9.5} {:>9.5} {:>9.5} {:>10.2e}",
                    m.f[k] / 1e3,
                    m.t[k].abs(),
                    te[k].abs(),
                    th[k].abs(),
                    ind
                );
            }
        }
    }
}

pub fn nearest(f: &[f64], target: f64) -> Option<usize> {
    f.iter()
        .enumerate()
        .min_by(|a, b| (a.1 - target).abs().partial_cmp(&(b.1 - target).abs()).unwrap())
        .map(|(i, _)| i)
}

fn main() {
    let arg = std::env::args().nth(1).unwrap_or_else(|| "all".into());
    let sys = mix(water(), solid(&materials::STEEL));
    let _ = &sys;
    match arg.as_str() {
        "stiffcontact" => {
            for n in [256usize] {
                let (dp, u) = stiff_contact_probe(n, 200);
                println!("stiff contact N={n}: max|dp|/p {dp:.3e}, max|u| {u:.3e}");
            }
        }
        "speed" => speed_check(&sys),
        "solid" => {
            let c = sys.material_c_l(1);
            for n in [800usize, 1600] {
                let v = solid_pulse_speed(&sys, n, 1.5e-4);
                println!("pure steel pulse N={n}: speed {:.1} m/s vs c_L {:.1} ({:+.2e})", v, c, v / c - 1.0);
            }
        }
        "five" => {
            for n in [256usize, 512] {
                let (r, rf, pr) = fiveeq_fresnel(n);
                println!("FiveEq water/steel half-space N={n}: peak |r| {r:.4} vs Fresnel {rf:.4}, pk(refl) {pr:.3e}");
            }
        }
        "scatter" => {
            for (rho, c) in [(2000.0f64, 2000.0f64), (4000.0, 3000.0), (7850.0, 5960.0)] {
                let sys = mix(water(), (rho, nasg_at(rho, c, P0, GAMMA), 0.0));
                let (r, rf) = fresnel_ratio_sys(&sys, 256);
                println!("contrast rho={rho:.0} c={c:.0}: peak |r| {r:.4} vs Fresnel {rf:.4}");
            }
        }
        "fresnel" => {
            for n in [256usize, 512] {
                let (r, rf) = fresnel_ratio(n);
                println!("fluid-fluid interface N={n}: peak |r| {:.4} vs Fresnel {:.4} ({:+.2e})", r, rf, r / rf - 1.0);
            }
        }
        "contact" => {
            for n in [256usize, 512] {
                let dev = contact_probe(n, 500);
                println!("ideal-gas contact N={n}: max |dp|/p after 500 steps = {:.3e}", dev);
            }
        }
        "debug" => {
            for n in [128usize, 256] {
                let cfg = plate_config(n);
                let with = run(&sys, &cfg, true, None);
                let without = run(&sys, &cfg, false, Some(with.dt));
                let stat = |r: &RunResult| {
                    let mn = r.pl.iter().cloned().fold(f64::MAX, f64::min);
                    let mx = r.pl.iter().cloned().fold(f64::MIN, f64::max);
                    (mn, mx)
                };
                let (a, b) = stat(&with);
                let (c, d) = stat(&without);
                println!("N={n}: plate steps {} dt {:.2e}, pl in [{a:.4},{b:.4}]; ref steps {} dt {:.2e}, pl in [{c:.4},{d:.4}]", with.steps, with.dt, without.steps, without.dt);
                if with.pl.len() == without.pl.len() {
                    let refl: Vec<f64> = with.pl.iter().zip(&without.pl).map(|(a, b)| a - b).collect();
                    let inc: Vec<f64> = without.pl.iter().map(|a| a - P0).collect();
                    let tra: Vec<f64> = with.pr.iter().map(|a| a - P0).collect();
                    let pk = |v: &[f64]| v.iter().cloned().fold(0.0f64, |m, x| m.max(x.abs()));
                    println!("   peaks: incident {:.3}, reflected {:.3}, transmitted {:.3} Pa", pk(&inc), pk(&refl), pk(&tra));
                }
                for k in (0..with.pl.len()).step_by((with.pl.len() / 12).max(1)) {
                    println!("   plate t={:.3e} pl={:+.5} pr={:+.5}", with.t[k], with.pl[k], with.pr[k]);
                }
                for k in (0..without.pl.len()).step_by((without.pl.len() / 6).max(1)) {
                    println!("   ref   t={:.3e} pl={:+.5} pr={:+.5}", without.t[k], without.pl[k], without.pr[k]);
                }
            }
        }
        "debug5" => {
            let (sys, media) = fiveeq_system();
            let cfg = plate_config(256);
            let w = run5(&sys, &media, 0, &cfg, true, None);
            let wo = run5(&sys, &media, 0, &cfg, false, Some(w.dt));
            let pk = |v: &[f64]| v.iter().cloned().fold(0.0f64, |m, x| m.max(x.abs()));
            let inc: Vec<f64> = wo.pl.iter().map(|a| a - P0).collect();
            let refl: Vec<f64> = w.pl.iter().zip(&wo.pl).map(|(a, b)| a - b).collect();
            let tra: Vec<f64> = w.pr.iter().map(|a| a - P0).collect();
            println!("FiveEq plate N=256 steps {}: pk inc {:.3}, refl {:.3}, tra {:.3}", w.steps, pk(&inc), pk(&refl), pk(&tra));
            let (mut bi, mut bv) = (0, 0.0);
            for (i, &v) in tra.iter().enumerate() {
                if v.abs() > bv {
                    bv = v.abs();
                    bi = i;
                }
            }
            println!("  transmitted peak at t={:.3e}: {:.4}", w.t[bi], tra[bi]);
        }
        "plate" => plate_report5(),
        "laminate" => laminate_report5(),
        "gprmix" => {
            speed_check(&sys);
            plate_report(&sys);
        }
        _ => {
            plate_report5();
            laminate_report5();
        }
    }
}
