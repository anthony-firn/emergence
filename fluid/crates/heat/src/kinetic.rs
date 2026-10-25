//! Kinetic (relaxation-time) lattice thermal conductivity.
//!
//! The Callaway/Holland expression used here is the branch-resolved RTA
//! integral
//!
//! ```text
//! κ = (1/3) Σ_branches deg ∫ C(ω) v_g(ω)² τ(ω) g(ω) dω ,
//! ```
//!
//! with `g(ω) = q²/(2π²) / v_g` the Debye density of states and
//! `τ^{-1} = τ_U^{-1} + τ_N^{-1} + τ_I^{-1} + τ_B^{-1}`.
//!
//! The optional [`kappa_callaway`] adds the Callaway normal-process drift
//! correction, which is positive and depends on the frequency dependence of
//! the rates.

use crate::linalg::gauss_legendre_01;
use crate::phonon::{mode_heat_capacity, Branch, HBAR, KB};
use crate::scattering::ScatteringParams;
use std::f64::consts::PI;

/// Number of Gauss–Legendre points used for the wavevector integrals.
const NQ: usize = 400;

/// Result of the kinetic integral, split into the RTA part `kappa_1` and the
/// Callaway normal-drift part `kappa_2`.
#[derive(Clone, Copy, Debug, Default)]
pub struct Kinetic {
    pub kappa_1: f64,
    pub kappa_2: f64,
}

impl Kinetic {
    pub fn total(&self) -> f64 {
        self.kappa_1 + self.kappa_2
    }
}

/// Relaxation-time (RTA) conductivity `κ_1` (W/(m K)) at temperature `t`.
pub fn kappa_rta(branches: &[Branch], p: &ScatteringParams, t: f64) -> f64 {
    kappa_split(branches, p, t).kappa_1
}

/// Full Callaway conductivity `κ_1 + κ_2` (W/(m K)).
///
/// `κ_2 = (1/3) (∫ w τ_c/τ_N dq)² / (∫ w τ_c/(τ_N τ_R) dq)` with the same
/// mode weight `w = C v_g² q²/(2π²)`. This is the standard Callaway
/// normal-process correction; it reduces to `κ_1` when `β → 0` and is exact
/// for a constant relaxation time (where normal processes do not resist).
pub fn kappa_callaway(branches: &[Branch], p: &ScatteringParams, t: f64) -> f64 {
    kappa_split(branches, p, t).total()
}

/// Both Callaway parts, for diagnostics and tests.
pub fn kappa_split(branches: &[Branch], p: &ScatteringParams, t: f64) -> Kinetic {
    let mut out = Kinetic::default();
    for b in branches {
        let qmax = b.q_max();
        let (xg, wg) = gauss_legendre_01(NQ);
        let mut k1 = 0.0;
        let mut j2 = 0.0;
        let mut k3 = 0.0;
        for i in 0..NQ {
            let q = xg[i] * qmax;
            let wq = wg[i] * qmax;
            let omega = b.omega(q);
            let vg = b.group_velocity(q);
            let x = HBAR * omega / (KB * t);
            let c = mode_heat_capacity(x);
            let g = q * q / (2.0 * PI * PI);
            let weight = c * vg * vg * g;
            let rtot = p.total(b, omega, t);
            let rn = p.normal(b, omega, t);
            let rr = p.resistive(b, omega, t);
            let tc = if rtot > 0.0 { 1.0 / rtot } else { 0.0 };
            k1 += wq * weight * tc;
            if p.beta > 0.0 && rr > 0.0 {
                j2 += wq * weight * tc * rn;
                k3 += wq * weight * tc * rn * rr;
            }
        }
        out.kappa_1 += b.degeneracy * k1 / 3.0;
        if k3 > 0.0 {
            out.kappa_2 += b.degeneracy * j2 * j2 / k3 / 3.0;
        }
    }
    out
}

/// Spectral mean free path at the dominant (thermal) frequency, weighted by
/// the heat-capacity integrand. Returns `(Λ_eff, ω_dom)` in metres and rad/s.
///
/// The indicator for the diffusive→ballistic fold is the phonon Knudsen
/// number `Kn = Λ_eff / L`.
pub fn mean_free_path(branches: &[Branch], p: &ScatteringParams, t: f64) -> (f64, f64) {
    let mut num = 0.0;
    let mut den = 0.0;
    let mut num_omega = 0.0;
    for b in branches {
        let qmax = b.q_max();
        let (xg, wg) = gauss_legendre_01(NQ);
        for i in 0..NQ {
            let q = xg[i] * qmax;
            let wq = wg[i] * qmax;
            let omega = b.omega(q);
            let vg = b.group_velocity(q);
            let x = HBAR * omega / (KB * t);
            let c = mode_heat_capacity(x);
            let g = q * q / (2.0 * PI * PI);
            let weight = c * vg * vg * g * wq;
            let rtot = p.total(b, omega, t);
            let tau = if rtot > 0.0 { 1.0 / rtot } else { f64::INFINITY };
            let mfp = vg * tau;
            num += b.degeneracy * weight * mfp;
            num_omega += b.degeneracy * weight * omega;
            den += b.degeneracy * weight;
        }
    }
    if den > 0.0 {
        (num / den, num_omega / den)
    } else {
        (f64::INFINITY, 0.0)
    }
}

/// Debye temperature of the branch set, defined as the heat-capacity-weighted
/// branch cutoff temperature (a diagnostic; `≈` the usual Debye temperature).
pub fn effective_debye_temperature(branches: &[Branch], t: f64) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for b in branches {
        let qmax = b.q_max();
        let (xg, wg) = gauss_legendre_01(NQ);
        for i in 0..NQ {
            let q = xg[i] * qmax;
            let wq = wg[i] * qmax;
            let omega = b.omega(q);
            let x = HBAR * omega / (KB * t);
            let c = mode_heat_capacity(x);
            let g = q * q / (2.0 * PI * PI);
            num += b.degeneracy * wq * c * g * b.theta();
            den += b.degeneracy * wq * c * g;
        }
    }
    if den > 0.0 {
        num / den
    } else {
        0.0
    }
}
