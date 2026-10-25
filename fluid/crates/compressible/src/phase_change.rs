//! Liquid–vapour phase change (cavitation, flash boiling, condensation) for
//! the five-equation model with the two phases of one substance, each with
//! a NASG equation of state (`eos::materials::WATER_*_NASG`).
//!
//! Model ladder. Mass transfer is a relaxation far faster than the flow
//! (molecular kinetics at the interface), so where it is active the mixture
//! sits at full thermodynamic equilibrium: equal pressure, temperature and
//! Gibbs free energy (the homogeneous equilibrium model, HEM). Where it is
//! not — pure liquid above saturation, pure vapour — nothing is done. The
//! emergent boundary is the saturation surface p = p_sat(T), shifted by a
//! nucleation threshold: real liquid sustains tension until nuclei grow,
//! which is sub-grid physics (nuclei content, classical nucleation theory)
//! summarized here by `tension`, the pressure deficit below p_sat needed
//! to start vapour formation in a cell that has none.
//!
//! The relaxation conserves ρ, ρu and E exactly (only α, α_kρ_k change), and
//! the new state is consistent with the five-equation mixture closure: with
//! ρ_k = 1/v_k(p, T), ρe = Σ α_kρ_k e_k(p, ρ_k) reproduces p.

use crate::eos::Nasg;
use crate::multiphase::FiveEq;
use crate::system::System;

type V = [f64; 7];

/// NASG specific volume and energy at (p, T).
#[inline]
pub fn v_pt(e: &Nasg, p: f64, t: f64) -> f64 {
    (e.gamma - 1.0) * e.cv * t / (p + e.p_inf) + e.b
}

#[inline]
pub fn e_pt(e: &Nasg, p: f64, t: f64) -> f64 {
    (p + e.gamma * e.p_inf) / (p + e.p_inf) * e.cv * t + e.q
}

/// ∂g/∂T at constant p (= −s).
#[inline]
fn dgdt(e: &Nasg, p: f64, t: f64) -> f64 {
    -e.q_prime - e.cv * (e.gamma * t.ln() - (e.gamma - 1.0) * (p + e.p_inf).ln())
}

/// Saturation temperature at pressure p (g_l = g_v), Newton from T₀.
pub fn t_sat(l: &Nasg, v: &Nasg, p: f64) -> f64 {
    t_sat_from(l, v, p, 400.0)
}

/// Saturation temperature with a starting guess (Newton; the iteration is
/// safeguarded to stay within a factor 3 per step).
pub fn t_sat_from(l: &Nasg, v: &Nasg, p: f64, t0: f64) -> f64 {
    let mut t = t0;
    for _ in 0..200 {
        let f = l.gibbs(p, t) - v.gibbs(p, t);
        let df = dgdt(l, p, t) - dgdt(v, p, t);
        let tn = (t - f / df).clamp(0.3 * t, 3.0 * t);
        if (tn - t).abs() < 1e-12 * t {
            return tn;
        }
        t = tn;
    }
    t
}

/// Saturation pressure at temperature T (bisection in ln p).
pub fn p_sat(l: &Nasg, v: &Nasg, t: f64) -> f64 {
    let (mut lo, mut hi) = (1e-6f64.ln(), 1e10f64.ln());
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        // T_sat increases with p
        if t_sat_from(l, v, mid.exp(), t) < t {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (0.5 * (lo + hi)).exp()
}

/// Coexistence of phases `l` (Y = 0) and `v` (Y = 1) at specific volume
/// `vol` and energy `energy`: equal p, T and Gibbs energy; Y is the mass
/// fraction of `v`. Scans ln p over 1e-6 … 1e10 Pa for a sign change of the
/// energy residual along the coexistence line, then bisects.
pub fn pair_equilibrium(l: &Nasg, v: &Nasg, vol: f64, energy: f64, t_guess: f64) -> Option<Equilibrium> {
    pair_equilibrium_near(l, v, vol, energy, t_guess, None)
}

/// As `pair_equilibrium`, searching first within a factor e⁴ of `p_guess`
/// (a cell already near equilibrium needs 12 saturation solves to bracket
/// instead of 64); falls back to the full scan.
pub fn pair_equilibrium_near(l: &Nasg, v: &Nasg, vol: f64, energy: f64, t_guess: f64, p_guess: Option<f64>) -> Option<Equilibrium> {
    let state = |p: f64| -> (f64, f64, f64) {
        let t = t_sat_from(l, v, p, t_guess);
        let (vl, vv) = (v_pt(l, p, t), v_pt(v, p, t));
        let y = (vol - vl) / (vv - vl);
        let e = y * e_pt(v, p, t) + (1.0 - y) * e_pt(l, p, t);
        (energy - e, t, y)
    };
    let scan = |a: f64, b: f64, n: usize| -> Option<(f64, f64, f64)> {
        let mut prev: Option<(f64, f64)> = None;
        for i in 0..=n {
            let lp = a + (b - a) * i as f64 / n as f64;
            let (r, t, _) = state(lp.exp());
            if !r.is_finite() || !t.is_finite() || t <= 0.0 {
                prev = None;
                continue;
            }
            if let Some((lp0, r0)) = prev {
                if r0 * r <= 0.0 {
                    return Some((lp0, lp, r0));
                }
            }
            prev = Some((lp, r));
        }
        None
    };
    let (a, b) = (1e-6f64.ln(), 1e10f64.ln());
    let local = p_guess.filter(|p| p.is_finite() && *p > 1e-6).and_then(|p| {
        let lp = p.ln();
        scan((lp - 4.0).max(a), (lp + 4.0).min(b), 12)
    });
    let (mut lo, mut hi, r_lo) = local.or_else(|| scan(a, b, 64))?;
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let (r, _, _) = state(mid.exp());
        if (r > 0.0) == (r_lo > 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let p = (0.5 * (lo + hi)).exp();
    let (_, t, y) = state(p);
    (0.0..=1.0).contains(&y).then_some(Equilibrium { p, t, y })
}

#[derive(Clone, Copy, Debug)]
pub struct PhaseChange {
    /// index of the liquid phase in `FiveEq::eos` (the other is its vapour)
    pub liquid: usize,
    /// pressure deficit below p_sat(T_l) that nucleates vapour in a cell
    /// without any (Pa); 0 = equilibrium (no tension)
    pub tension: f64,
    /// volume fraction above which a phase counts as present
    pub present: f64,
    /// starting temperature for the saturation Newton iteration (K)
    pub t_guess: f64,
}

/// Result of an equilibrium solve.
#[derive(Clone, Copy, Debug)]
pub struct Equilibrium {
    pub p: f64,
    pub t: f64,
    /// vapour mass fraction
    pub y: f64,
}

impl PhaseChange {
    fn eos<'a>(&self, sys: &'a FiveEq) -> (&'a Nasg, &'a Nasg) {
        (&sys.eos[self.liquid], &sys.eos[1 - self.liquid])
    }

    /// Two-phase equilibrium (p, T, Y) at given specific volume and energy,
    /// or `None` if the state is single-phase (Y would leave [0, 1]).
    pub fn equilibrium(&self, sys: &FiveEq, vol: f64, energy: f64) -> Option<Equilibrium> {
        self.equilibrium_near(sys, vol, energy, None)
    }

    /// `equilibrium`, searching near `p_guess` first.
    pub fn equilibrium_near(&self, sys: &FiveEq, vol: f64, energy: f64, p_guess: Option<f64>) -> Option<Equilibrium> {
        let (l, v) = self.eos(sys);
        pair_equilibrium_near(l, v, vol, energy, self.t_guess, p_guess)
    }

    /// Relax the conserved state `u` to thermodynamic equilibrium where the
    /// criterion holds. Returns true if the cell was changed.
    pub fn relax(&self, sys: &FiveEq, u: &mut V) -> bool {
        let w = sys.to_prim(u);
        let (il, iv) = (self.liquid, 1 - self.liquid);
        let alpha_l = if il == 0 { w[6] } else { 1.0 - w[6] };
        let alpha_v = 1.0 - alpha_l;
        let (l, _v) = self.eos(sys);
        let p = w[5];
        let both = alpha_l > self.present && alpha_v > self.present;
        if !both {
            if alpha_l <= self.present {
                return false; // pure vapour: condensation needs liquid nuclei too
            }
            // pure liquid: cavitate only beyond the nucleation threshold
            let rho_l = w[il] / alpha_l.max(1e-300);
            let t_l = l.temperature(rho_l, p);
            if p >= p_sat(l, sys.eos.get(iv).unwrap(), t_l) - self.tension {
                return false;
            }
        }
        let Some(un) = self.equilibrium_state(sys, u) else {
            return false;
        };
        // conserve E exactly (to round-off the closure reproduces it)
        u[0] = un[0];
        u[1] = un[1];
        u[6] = un[6];
        true
    }

    /// The two-phase equilibrium state with the same ρ, ρu and E as `u`
    /// (None if the equilibrium is single-phase).
    pub fn equilibrium_state(&self, sys: &FiveEq, u: &V) -> Option<V> {
        let w = sys.to_prim(u);
        let (il, iv) = (self.liquid, 1 - self.liquid);
        let rho = w[0] + w[1];
        let ke = 0.5 * rho * (w[2] * w[2] + w[3] * w[3] + w[4] * w[4]);
        let energy = (u[5] - ke) / rho;
        let eq = self.equilibrium_near(sys, 1.0 / rho, energy, Some(w[5]))?;
        let (_l, v) = self.eos(sys);
        let av = eq.y * rho * v_pt(v, eq.p, eq.t);
        let mut wn = w;
        wn[iv] = eq.y * rho;
        wn[il] = (1.0 - eq.y) * rho;
        wn[6] = if il == 0 { 1.0 - av } else { av };
        wn[5] = eq.p;
        Some(sys.to_cons(&wn))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eos::materials::{WATER_LIQUID_NASG, WATER_VAPOR_NASG};

    /// NASG water saturation curve vs IAPWS steam tables (the NASG fit is
    /// for 300–500 K; Le Métayer & Saurel report a few % in p_sat).
    #[test]
    fn saturation_curve_matches_steam_tables() {
        let (l, v) = (WATER_LIQUID_NASG, WATER_VAPOR_NASG);
        for (t, p_table) in [(298.15, 3_169.9), (323.15, 12_352.0), (373.15, 101_418.0), (423.15, 476_160.0), (473.15, 1_554_900.0)] {
            let p = p_sat(&l, &v, t);
            let err = p / p_table - 1.0;
            println!("T = {t} K: p_sat NASG {p:.0} Pa vs IAPWS {p_table:.0} Pa ({:+.1}%)", 100.0 * err);
            assert!(err.abs() < 0.25, "p_sat at {t}: {p} vs {p_table}");
        }
    }
}
