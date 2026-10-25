//! Finite-rate interphase mass transfer for the five-equation model:
//! cavitation / condensation with nucleation, and gas dissolution /
//! exsolution. Both act per cell after each step (`Solver::kinetics`),
//! conserve ρ, ρu and E exactly, and move mass only as fast as the
//! sub-grid physics allows.
//!
//! **Cavitation** (liquid ↔ its vapour). A pure liquid below its vapour
//! pressure is metastable; vapour appears only once nuclei grow:
//! * heterogeneous nuclei (gas pockets of radius R_n, number density n per
//!   unit liquid volume) become unstable below the Blake threshold
//!     p_B = p_v − 4σ/(3R_c),  R_c = √(3 p_g0 R_n³ / 2σ),
//!   p_g0 = p_ref − p_v + 2σ/R_n (nucleus in equilibrium at p_ref);
//! * homogeneous nucleation (classical nucleation theory) creates nuclei
//!   at the rate J = N_l √(2σ/πm) exp(−16πσ³ / (3k_BT (p_v − p)²)), counted
//!   in a transported scalar (nuclei per unit liquid mass).
//! Vapour then grows / condenses at the Schnerr–Sauer rate (bubble
//! population of n bubbles per liquid volume, Rayleigh–Plesset inertial
//! growth velocity √(2|p_v − p| / 3ρ_l)),
//!   ṁ = (ρ_lρ_v/ρ) (3α_vα_l/R_b) sign(p_v − p) √(2|p_v − p|/3ρ_l),
//!   R_b = (3α_v / 4πnα_l)^{1/3},
//! and never beyond the thermodynamic equilibrium state
//! (`PhaseChange::equilibrium_state`), which it relaxes toward. As n → ∞
//! the model tends to the equilibrium (HEM) relaxation. Liquid tension
//! before nucleation is represented; once vapour exists the transfer is
//! at least what lifts the mixture pressure to 5 % of p_v (a
//! single-pressure closure cannot hold vapour beside liquid in tension).
//!
//! **Dissolution** (a gas phase ↔ gas dissolved in the liquid, carried by
//! the liquid as a scalar: dissolved mass per unit liquid mass). Henry's
//! law gives the saturation concentration c_s = k_H p; transfer through
//! the interfacial area of the sub-grid bubbles (n_b per unit volume,
//! radius R_b = (3α_g/4πn_b)^{1/3}, area A = 3α_g/R_b) with the
//! mass-transfer coefficient k_L = Sh D / 2R_b:
//!   ∂c/∂t = k_L A (c_s − c) / α_l,
//! integrated exactly over the step for frozen coefficients (so any Δt is
//! stable). The gas loses the mass the liquid gains (and the volume that
//! mass occupied; the volume of the dissolved gas itself is neglected).

use crate::eos::Nasg;
use crate::multiphase::FiveEq;
use crate::phase_change::{p_sat, PhaseChange};
use crate::system::System;

type V = [f64; 7];

pub const K_B: f64 = 1.380649e-23;

/// Heterogeneous nuclei: free gas pockets.
#[derive(Clone, Copy, Debug)]
pub struct Nuclei {
    /// number per unit liquid volume (m⁻³)
    pub density: f64,
    /// radius at p_ref (m)
    pub radius: f64,
    /// pressure at which the nuclei are in equilibrium (Pa)
    pub p_ref: f64,
}

/// Blake threshold of a gas nucleus of radius `r_n` in equilibrium at
/// `p_ref`, for vapour pressure `p_v` and surface tension `sigma`.
pub fn blake_threshold(r_n: f64, p_ref: f64, p_v: f64, sigma: f64) -> f64 {
    let pg0 = (p_ref - p_v + 2.0 * sigma / r_n).max(0.0);
    let rc = (3.0 * pg0 * r_n.powi(3) / (2.0 * sigma)).sqrt();
    p_v - 4.0 * sigma / (3.0 * rc)
}

/// Classical homogeneous nucleation rate (m⁻³ s⁻¹) of vapour in a liquid
/// of density `rho_l` (molecular mass `m`) at temperature `t`, pressure
/// `p` below the vapour pressure `p_v`.
pub fn cnt_rate(rho_l: f64, m: f64, sigma: f64, t: f64, p: f64, p_v: f64) -> f64 {
    if p >= p_v {
        return 0.0;
    }
    let j0 = rho_l / m * (2.0 * sigma / (std::f64::consts::PI * m)).sqrt();
    let dg = 16.0 * std::f64::consts::PI * sigma.powi(3) / (3.0 * (p_v - p).powi(2));
    j0 * (-dg / (K_B * t)).exp()
}

#[derive(Clone, Copy, Debug)]
pub struct Cavitation {
    /// equilibrium machinery (liquid index, saturation-curve guess)
    pub pc: PhaseChange,
    /// surface tension (N/m)
    pub sigma: f64,
    /// heterogeneous nuclei (None: only homogeneous nucleation)
    pub nuclei: Option<Nuclei>,
    /// homogeneous nucleation into scalar `cnt_scalar` (nuclei per unit
    /// liquid mass, carrier = the liquid's partial density); the molecular
    /// mass of the liquid (kg)
    pub cnt_scalar: Option<usize>,
    pub molecule_mass: f64,
    /// condensation rate relative to vaporization (Schnerr–Sauer: 1)
    pub condensation: f64,
}

impl Cavitation {
    /// Advance the cell by Δt; returns true if it changed.
    pub fn apply(&self, sys: &FiveEq, u: &mut V, sc: &mut [f64], dt: f64) -> bool {
        let il = self.pc.liquid;
        let iv = 1 - il;
        let (l, v): (&Nasg, &Nasg) = (&sys.eos[il], &sys.eos[iv]);
        let w = sys.to_prim(u);
        let alpha_l = if il == 0 { w[6] } else { 1.0 - w[6] };
        let alpha_v = 1.0 - alpha_l;
        if alpha_l <= self.pc.present {
            return false; // pure vapour: nothing to condense onto here
        }
        let p = w[5];
        let rho_l = w[il] / alpha_l;
        let t_l = l.temperature(rho_l, p);
        if !(t_l.is_finite() && t_l > 0.0) {
            return false;
        }
        let pv = p_sat(l, v, t_l);
        let mut changed = false;
        // homogeneous nucleation: accumulate nuclei
        let mut n_cnt = 0.0;
        if let Some(k) = self.cnt_scalar {
            if p < pv {
                let j = cnt_rate(rho_l, self.molecule_mass, self.sigma, t_l, p, pv);
                if j > 0.0 {
                    sc[k] += alpha_l * j * dt; // (α_lρ_l)·(nuclei per liquid mass)
                    changed = true;
                }
            }
            n_cnt = sc[k] * rho_l / w[il].max(1e-300);
        }
        let n_het = self.nuclei.map_or(0.0, |nu| nu.density);
        // active bubble population: heterogeneous nuclei only below Blake
        let het_active = self.nuclei.is_some_and(|nu| alpha_v > self.pc.present || p < blake_threshold(nu.radius, nu.p_ref, pv, self.sigma));
        let n = if het_active { n_het } else { 0.0 } + n_cnt;
        if n <= 0.0 {
            return changed;
        }
        if alpha_v <= self.pc.present && p >= pv {
            return changed; // no vapour and nothing to grow it
        }
        // Schnerr–Sauer rate
        let rho = w[0] + w[1];
        let rho_v = if alpha_v > 1e-12 && w[iv] > 0.0 { w[iv] / alpha_v } else { 1.0 / crate::phase_change::v_pt(v, pv.max(1.0), t_l) };
        // seed: a just-activated population starts at its critical radius
        let av_eff = if alpha_v > self.pc.present {
            alpha_v
        } else {
            let rc = 2.0 * self.sigma / (pv - p).max(1e-300);
            (n * 4.0 / 3.0 * std::f64::consts::PI * rc.powi(3) * alpha_l).min(1e-3).max(self.pc.present * 1.0001)
        };
        let rb = (3.0 * av_eff / (4.0 * std::f64::consts::PI * n * alpha_l)).cbrt();
        let dpv = pv - p;
        let mut mdot = rho_l * rho_v / rho * 3.0 * av_eff * alpha_l / rb * dpv.signum() * (2.0 * dpv.abs() / (3.0 * rho_l)).sqrt();
        if mdot < 0.0 {
            mdot *= self.condensation;
        }
        let dy = mdot * dt / rho; // change of the vapour mass fraction
        if dy == 0.0 {
            return changed;
        }
        // target: the equilibrium state (or all liquid when condensing past it)
        let target = match self.pc.equilibrium_state(sys, u) {
            Some(t) => t,
            None => {
                if dy > 0.0 {
                    return changed; // the equilibrium is pure liquid: no vaporization
                }
                let mut wt = w;
                wt[il] = rho;
                wt[iv] = 0.0;
                wt[6] = if il == 0 { 1.0 } else { 0.0 };
                let mut ut = sys.to_cons(&wt);
                ut[5] = u[5];
                ut
            }
        };
        let y = w[iv] / rho;
        let yt = target[iv] / rho;
        if (yt - y) * dy <= 0.0 {
            return changed; // already past equilibrium in this direction
        }
        let mut f = (dy / (yt - y)).clamp(0.0, 1.0);
        // with vapour in the cell the single-pressure closure cannot hold
        // liquid in tension (the bubbles' interiors are at p_v): transfer at
        // least enough to bring the mixture pressure to 5 % of p_v, after
        // which the rate law governs the approach to saturation
        let blend = |f: f64| -> V {
            let mut un = *u;
            for k in [0usize, 1, 6] {
                un[k] += f * (target[k] - u[k]);
            }
            un
        };
        let p_floor = 0.05 * pv;
        if dy > 0.0 && sys.to_prim(&blend(f))[5] < p_floor && sys.to_prim(&blend(1.0))[5] >= p_floor {
            let (mut lo, mut hi) = (f, 1.0);
            for _ in 0..50 {
                let mid = 0.5 * (lo + hi);
                if sys.to_prim(&blend(mid))[5] < p_floor {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            f = hi;
        }
        *u = blend(f);
        true
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Dissolution {
    /// FiveEq phase index of the liquid (the other is the gas)
    pub liquid: usize,
    /// scalar index of the dissolved gas (per unit liquid mass)
    pub scalar: usize,
    /// Henry coefficient: saturation concentration per pressure (kg m⁻³ Pa⁻¹)
    pub henry: f64,
    /// diffusivity of the dissolved gas in the liquid (m²/s)
    pub diffusivity: f64,
    /// Sherwood number of the bubbles (2: stagnant)
    pub sherwood: f64,
    /// sub-grid bubble number density (m⁻³)
    pub bubble_density: f64,
    /// gas volume fraction below which there is no interface
    pub present: f64,
}

impl Dissolution {
    /// Advance the cell by Δt; returns true if it changed.
    pub fn apply(&self, sys: &FiveEq, u: &mut V, sc: &mut [f64], dt: f64) -> bool {
        let il = self.liquid;
        let ig = 1 - il;
        let w = sys.to_prim(u);
        let alpha_l = if il == 0 { w[6] } else { 1.0 - w[6] };
        let alpha_g = 1.0 - alpha_l;
        if alpha_g <= self.present || alpha_l <= self.present {
            return false;
        }
        let c = sc[self.scalar] / alpha_l; // dissolved mass per liquid volume
        let cs = self.henry * w[5].max(0.0);
        let rb = (3.0 * alpha_g / (4.0 * std::f64::consts::PI * self.bubble_density)).cbrt();
        let area = 3.0 * alpha_g / rb;
        let kl = self.sherwood * self.diffusivity / (2.0 * rb);
        let lam = kl * area / alpha_l;
        let dc = (cs - c) * (1.0 - (-lam * dt).exp());
        // mass per mixture volume into the liquid, within what exists
        let dm = (alpha_l * dc).min(w[ig]).max(-sc[self.scalar]);
        if dm == 0.0 {
            return false;
        }
        let rho_g = w[ig] / alpha_g;
        let rho_l = w[il] / alpha_l;
        let mut wn = w;
        wn[ig] = w[ig] - dm;
        wn[il] = w[il] + dm;
        // at fixed cell volume the phases stay in mechanical equilibrium:
        // the volume the transferred gas occupied (dm/ρ_g) is shared in
        // proportion to α_k/(ρ_kc_k²), as the Kapila term shares compression
        // — the remaining gas expands into nearly all of it and the cell
        // pressure drops a little (the flow refills it afterwards). The
        // volume of the dissolved gas itself is neglected (N₂ at saturation
        // changes water's volume by ~2e-5). A gas phase that is used up
        // leaves pure liquid.
        let ag_new = if wn[ig] <= 1e-9 * w[ig].max(1e-300) {
            wn[il] += wn[ig];
            wn[ig] = 0.0;
            0.0
        } else {
            let z = |e: &Nasg, rk: f64| e.gamma * (w[5] + e.p_inf) / (1.0 - e.b * rk).max(1e-12);
            let (zl, zg) = (z(&sys.eos[il], rho_l), z(&sys.eos[ig], rho_g));
            let share = (alpha_g / zg) / (alpha_l / zl + alpha_g / zg);
            (alpha_g - dm / rho_g * (1.0 - share)).clamp(0.0, 1.0)
        };
        wn[6] = if il == 0 { 1.0 - ag_new } else { ag_new };
        let un = sys.to_cons(&wn);
        u[0] = un[0];
        u[1] = un[1];
        u[6] = un[6];
        sc[self.scalar] += dm;
        true
    }
}
