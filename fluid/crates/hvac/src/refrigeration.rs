//! The vapour-compression cycle: air conditioners, heat pumps, chillers.
//!
//! Single stage on real refrigerant properties (the reference equations of
//! state in `compressible::refrigerant_data`):
//! 1. evaporator exit: saturation pressure at the evaporating temperature,
//!    superheated by ΔT_sh;
//! 2. compressor exit: isentropic to the condensing pressure, then
//!    h₂ = h₁ + (h₂ₛ − h₁)/η_is;
//! 3. condenser exit: subcooled by ΔT_sc at the condensing pressure;
//! 4. isenthalpic expansion back to the evaporating pressure.
//!
//! Pressure drops in the heat exchangers and lines, and the compressor's
//! heat loss, are neglected. Transcritical operation (condensing above the
//! critical temperature: CO₂ heat pumps) needs a gas-cooler pressure as a
//! design variable and is refused.

use compressible::helmholtz::{Helmholtz, Props};
use compressible::refrigerant_data as rd;

/// The refrigerants with reference equations of state.
pub fn refrigerant(name: &str) -> Result<Helmholtz, String> {
    Ok(match name.to_ascii_uppercase().replace('-', "").as_str() {
        "R32" => rd::r32(),
        "R134A" => rd::r134a(),
        "R410A" => rd::r410a(),
        "R290" | "PROPANE" => rd::r290(),
        "R1234YF" => rd::r1234yf(),
        "R744" | "CO2" => rd::r744(),
        _ => return Err(format!("unknown refrigerant `{name}` (R32, R134a, R410A, R290, R1234yf, R744)")),
    })
}

/// Saturation temperature (K) at pressure `p`.
pub fn t_sat(eos: &Helmholtz, p: f64) -> Result<f64, String> {
    if p >= eos.p_c {
        return Err(format!("{:.0} kPa is above the critical pressure of {}", p / 1e3, eos.name));
    }
    let (mut lo, mut hi) = (eos.t_triple.max(0.6 * eos.t_c), eos.t_c * (1.0 - 1e-6));
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        let (ps, _, _) = eos.saturation(mid).ok_or("saturation failed")?;
        if ps < p {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-10 * hi {
            break;
        }
    }
    Ok(0.5 * (lo + hi))
}

/// Single-phase properties at (p, T): vapour starts from the ideal gas,
/// liquid from the saturated liquid density.
fn at(eos: &Helmholtz, p: f64, t: f64, liquid: bool) -> Result<Props, String> {
    let rho0 = if liquid {
        eos.saturation(t.min(eos.t_c * 0.999)).map(|s| s.1).ok_or("no liquid start")?
    } else {
        p / (eos.r() * t)
    };
    let rho = eos.density(p, t, rho0).ok_or_else(|| format!("no density at {:.0} kPa, {:.2} K", p / 1e3, t))?;
    Ok(eos.props(rho, t))
}

/// Vapour state at pressure `p` with entropy `s`.
fn vapour_ps(eos: &Helmholtz, p: f64, s: f64, t_low: f64) -> Result<Props, String> {
    let (mut lo, mut hi) = (t_low, t_low + 300.0);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if at(eos, p, mid, false)?.s < s {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-14 * hi {
            break;
        }
    }
    at(eos, p, 0.5 * (lo + hi), false)
}

#[derive(Clone, Copy, Debug)]
pub struct Cycle {
    pub p_evap: f64,
    pub p_cond: f64,
    /// Enthalpies at states 1–4 (J/kg).
    pub h: [f64; 4],
    pub t_discharge: f64,
    pub cop_cooling: f64,
    pub cop_heating: f64,
    /// Refrigerating effect and compressor work per kg (J/kg).
    pub q_evap: f64,
    pub w_comp: f64,
    /// Quality after expansion.
    pub quality_in: f64,
}

/// The cycle between evaporating `t_e` and condensing `t_c` (K).
pub fn cycle(eos: &Helmholtz, t_e: f64, t_c: f64, superheat: f64, subcool: f64, eta_is: f64) -> Result<Cycle, String> {
    if t_c >= eos.t_c {
        return Err(format!(
            "condensing at {:.1} °C is above {}'s critical temperature ({:.1} °C): transcritical cycles need a gas-cooler pressure and are not supported",
            t_c - 273.15,
            eos.name,
            eos.t_c - 273.15
        ));
    }
    if !(t_e < t_c) || !(0.0..=1.0).contains(&eta_is) || eta_is == 0.0 || superheat < 0.0 || subcool < 0.0 {
        return Err("need t_e < t_c, 0 < η_is ≤ 1 and non-negative superheat and subcooling".into());
    }
    let (p_e, rl_e, rv_e) = eos.saturation(t_e).ok_or("no saturation at the evaporating temperature")?;
    let (p_c, _, _) = eos.saturation(t_c).ok_or("no saturation at the condensing temperature")?;
    let s1 = if superheat > 0.0 { at(eos, p_e, t_e + superheat, false)? } else { eos.props(rv_e, t_e) };
    let s2s = vapour_ps(eos, p_c, s1.s, t_c)?;
    let h2 = s1.h + (s2s.h - s1.h) / eta_is;
    // Discharge temperature: the vapour at p_c with enthalpy h₂.
    let (mut lo, mut hi) = (s2s.t, s2s.t + 300.0);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if at(eos, p_c, mid, false)?.h < h2 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let t2 = 0.5 * (lo + hi);
    let s3 = if subcool > 0.0 {
        at(eos, p_c, t_c - subcool, true)?
    } else {
        let (_, rl, _) = eos.saturation(t_c).unwrap();
        eos.props(rl, t_c)
    };
    let h4 = s3.h;
    let (hl, hv) = (eos.props(rl_e, t_e).h, eos.props(rv_e, t_e).h);
    let x4 = (h4 - hl) / (hv - hl);
    let (q, w) = (s1.h - h4, h2 - s1.h);
    Ok(Cycle {
        p_evap: p_e,
        p_cond: p_c,
        h: [s1.h, h2, s3.h, h4],
        t_discharge: t2,
        cop_cooling: q / w,
        cop_heating: (h2 - s3.h) / w,
        q_evap: q,
        w_comp: w,
        quality_in: x4,
    })
}
