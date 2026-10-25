//! Single-particle model (SPM) of a lithium-ion cell: the coarse model that
//! remains when electrolyte gradients are negligible.
//!
//! Each porous electrode is folded onto one representative spherical particle
//! with a *uniform* reaction current density. Solid diffusion is resolved in
//! [`crate::particle`], the surface reaction is Butler–Volmer
//! ([`crate::interface`]), and the terminal voltage is
//!
//! ```text
//!   V = U_p(x_p^s) − U_n(x_n^s) + η_p − η_n − i_app R_Ω
//!   i_loc,i = i_app /(a_s,i L_i),   a_s,i = 3 ε_s,i / R_i
//! ```
//!
//! with `i_app` the applied current density (A m⁻², discharge positive).
//! Parameter set: C.-H. Chen et al., *J. Electrochem. Soc.* **167**, 080534
//! (2020), "LG M50" 21700 cell. Values below are taken from that paper (and
//! the PyBaMM `Chen2020` implementation); the entries that are estimates or
//! transcriptions rather than directly published are marked `ESTIMATE`.
//!
//! **Honesty note.** The negative-electrode OCP is the published Chen 2020
//! graphite fit. The positive-electrode OCP here is a *surrogate* sigmoid (the
//! published NMC fit was not available offline); the model therefore
//! reproduces the shape and the conservation laws, not the published voltage
//! digits. See `analysis/electrochem.md`.
//!
//! # The fold
//! SPM is valid when electrolyte transport is quasi-static and its
//! concentration is nearly uniform: indicator `t_e / t_d`, with
//! `t_e = L_e² / D_e`. [`crate::spm::electrolyte_indicator`] returns it.

use crate::interface::{butler_volmer_overpotential, tafel_slope};
use crate::particle::SphereParticle;
use crate::F_FARADAY;

/// One porous electrode folded to a single particle.
#[derive(Clone, Debug)]
pub struct Electrode {
    /// Coating thickness (m).
    pub thickness: f64,
    /// Active-material volume fraction.
    pub eps_s: f64,
    /// Particle radius (m).
    pub radius: f64,
    /// Solid diffusivity (m² s⁻¹).
    pub d_s: f64,
    /// Maximum solid concentration (mol m⁻³).
    pub c_max: f64,
    /// Initial solid concentration (mol m⁻³).
    pub c0: f64,
    /// Reaction rate constant (m^2.5 mol^-0.5 s^-1).
    pub k: f64,
    /// Anodic transfer coefficient.
    pub alpha_a: f64,
    /// Cathodic transfer coefficient.
    pub alpha_c: f64,
}

impl Electrode {
    /// Specific active surface area `a_s = 3 ε_s / R` (m² m⁻³).
    pub fn a_s(&self) -> f64 {
        3.0 * self.eps_s / self.radius
    }

    /// Local reaction current density for applied current density `i_app`
    /// (A m⁻²), positive = anodic (lithium leaves the particle).
    pub fn local_current(&self, i_app: f64) -> f64 {
        i_app / (self.a_s() * self.thickness)
    }
}

/// Positive-electrode OCP surrogate (NMC-like), decreasing with lithium
/// content. `ESTIMATE`: sigmoidal fit chosen to span 4.1→3.5 V over x∈[0.25,0.9].
pub fn ocp_nmc(x: f64) -> f64 {
    3.85 + 0.40 * (-4.0 * (x - 0.45)).tanh()
}

/// Negative-electrode OCP: graphite LG M50 fit of Chen et al. (2020).
pub fn ocp_graphite(x: f64) -> f64 {
    1.9793 * (-39.3631 * x).exp() + 0.2482
        - 0.0909 * (29.8538 * (x - 0.1234)).tanh()
        - 0.04478 * (14.9159 * (x - 0.2769)).tanh()
        - 0.02098 * (30.4444 * (x - 0.6100)).tanh()
}

/// Full cell parameters (area-specific; the electrode area cancels).
#[derive(Clone, Debug)]
pub struct CellParams {
    pub neg: Electrode,
    pub pos: Electrode,
    /// Electrolyte concentration for the exchange current (mol m⁻³).
    pub c_e: f64,
    /// Area-specific ohmic resistance (Ω m²). `ESTIMATE`.
    pub r_ohmic: f64,
    /// Temperature (K).
    pub t: f64,
}

impl CellParams {
    /// LG M50 parameter set (Chen et al. 2020). `ESTIMATE` marks values that
    /// were not available verbatim offline.
    pub fn lg_m50() -> Self {
        Self {
            neg: Electrode {
                thickness: 85.2e-6,
                eps_s: 0.75,
                radius: 5.86e-6,
                d_s: 3.3e-14,
                c_max: 33133.0,
                c0: 29866.0,
                k: 2.334e-11, // ESTIMATE of the Chen2020 rate constant
                alpha_a: 0.5,
                alpha_c: 0.5,
            },
            pos: Electrode {
                thickness: 75.6e-6,
                eps_s: 0.665,
                radius: 5.22e-6,
                d_s: 4.0e-15,
                c_max: 63104.0,
                c0: 17038.0,
                k: 6.48e-11, // ESTIMATE
                alpha_a: 0.5,
                alpha_c: 0.5,
            },
            c_e: 1000.0,
            r_ohmic: 2.0e-3, // ESTIMATE area-specific resistance
            t: 298.15,
        }
    }

    /// 1C current density (A m⁻²) from the positive-electrode usable range
    /// `Δx_p = 0.6`: `i_1C = F ε_s L c_max Δx_p / 3600`.
    pub fn i_ref(&self) -> f64 {
        F_FARADAY * self.pos.eps_s * self.pos.thickness * self.pos.c_max * 0.6 / 3600.0
    }

    /// Electrolyte transport time `L_e²/D_e` (s) with an effective electrolyte
    /// length `L_e` = separator + half-electrodes and `D_e` an estimate.
    pub fn electrolyte_time(&self) -> f64 {
        let l_e = 25e-6 + 0.5 * (self.neg.thickness + self.pos.thickness);
        let d_e = 2.0e-10; // ESTIMATE typical LiPF6/carbonate
        l_e * l_e / d_e
    }
}

/// Indicator for the SPM fold: `t_e / t_d` where `t_d = 3600/C` seconds.
/// The electrolyte is quasi-static (SPM valid) when this is small.
pub fn electrolyte_indicator(c_rate: f64, p: &CellParams) -> f64 {
    p.electrolyte_time() / (3600.0 / c_rate)
}

/// A running SPM cell.
pub struct Spm {
    pub params: CellParams,
    pub neg: SphereParticle,
    pub pos: SphereParticle,
}

impl Spm {
    /// Build from parameters with `nodes` shells per particle.
    pub fn new(params: CellParams, nodes: usize) -> Self {
        let neg = SphereParticle::new(nodes, params.neg.radius, params.neg.d_s, params.neg.c0);
        let pos = SphereParticle::new(nodes, params.pos.radius, params.pos.d_s, params.pos.c0);
        Self { params, neg, pos }
    }

    /// Equilibrium (open-circuit) voltage (V).
    pub fn ocv(&self) -> f64 {
        let xn = self.neg.surface() / self.params.neg.c_max;
        let xp = self.pos.surface() / self.params.pos.c_max;
        ocp_nmc(xp) - ocp_graphite(xn)
    }

    /// Exchange current density (A m⁻²) for an electrode at surface
    /// concentration `c_s`.
    fn i0(&self, e: &Electrode, c_s: f64) -> f64 {
        let cs = c_s.clamp(1e-6 * e.c_max, e.c_max - 1e-6 * e.c_max);
        e.k * F_FARADAY * (self.params.c_e * cs * (e.c_max - cs)).sqrt()
    }

    /// Overpotential (V) of an electrode carrying local current density `i`.
    fn eta(&self, e: &Electrode, c_s: f64, i: f64) -> f64 {
        let i0 = self.i0(e, c_s);
        butler_volmer_overpotential(i, i0, e.alpha_a, e.alpha_c, 1.0, self.params.t)
            .unwrap_or_else(|| {
                // Fall back to the Tafel slope times log(i/i0) away from i0.
                let s = tafel_slope(e.alpha_a, 1.0, self.params.t);
                s * (i.abs() / i0).log10() * i.signum()
            })
    }

    /// Terminal voltage (V) at applied current density `i_app` (A m⁻²).
    pub fn voltage(&self, i_app: f64) -> f64 {
        let p = &self.params;
        let cs_n = self.neg.surface();
        let cs_p = self.pos.surface();
        let i_n = p.neg.local_current(i_app);
        // The positive electrode is the cathode on discharge: its local
        // current is cathodic (negative).
        let i_p = -p.pos.local_current(i_app);
        let eta_n = self.eta(&p.neg, cs_n, i_n);
        let eta_p = self.eta(&p.pos, cs_p, i_p);
        let xn = cs_n / p.neg.c_max;
        let xp = cs_p / p.pos.c_max;
        ocp_nmc(xp) - ocp_graphite(xn) + eta_p - eta_n - i_app * p.r_ohmic
    }

    /// Advance the solid diffusion by `dt` at applied current density `i_app`.
    pub fn step(&mut self, dt: f64, i_app: f64) {
        let p = &self.params;
        let i_n = p.neg.local_current(i_app);
        let i_p = -p.pos.local_current(i_app);
        // Surface molar flux: outward positive for the negative particle;
        // inward (negative outward) for the positive particle.
        self.neg.step(dt, i_n / F_FARADAY);
        self.pos.step(dt, i_p / F_FARADAY);
    }

    /// Total solid lithium per unit electrode area (mol m⁻²), an invariant of
    /// the discharge.
    pub fn solid_lithium(&self) -> f64 {
        self.params.neg.eps_s * self.params.neg.thickness * self.neg.average()
            + self.params.pos.eps_s * self.params.pos.thickness * self.pos.average()
    }

    /// Constant-current discharge to a voltage cutoff.
    ///
    /// Returns `(capacity_Ah_per_m2, trace)` where each trace entry is
    /// `(t_s, V_V, x_n, x_p)`. `c_rate` is relative to
    /// [`CellParams::i_ref`].
    pub fn discharge(&self, c_rate: f64, v_cutoff: f64, dt: f64, max_steps: usize)
        -> (f64, Vec<(f64, f64, f64, f64)>)
    {
        let mut cell = Spm::new(self.params.clone(), self.neg.nodes);
        cell.neg.c.copy_from_slice(&self.neg.c);
        cell.pos.c.copy_from_slice(&self.pos.c);
        let i_app = c_rate * self.params.i_ref();
        let mut trace = Vec::new();
        let mut t = 0.0;
        for _ in 0..max_steps {
            let v = cell.voltage(i_app);
            trace.push((
                t,
                v,
                cell.neg.surface() / self.params.neg.c_max,
                cell.pos.surface() / self.params.pos.c_max,
            ));
            if v < v_cutoff {
                break;
            }
            cell.step(dt, i_app);
            t += dt;
        }
        let capacity = i_app * t / 3600.0 * 1000.0; // mAh m⁻²
        (capacity, trace)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ocps_are_monotone_decreasing() {
        for i in 0..100 {
            let x1 = i as f64 / 100.0;
            let x2 = x1 + 0.01;
            assert!(ocp_graphite(x1) > ocp_graphite(x2));
            assert!(ocp_nmc(x1) > ocp_nmc(x2));
        }
    }

    #[test]
    fn open_circuit_voltage_is_in_the_physical_window() {
        let p = CellParams::lg_m50();
        let cell = Spm::new(p, 30);
        let v = cell.ocv();
        assert!(v > 3.5 && v < 4.3, "initial OCV {v}");
    }

    #[test]
    fn solid_lithium_is_conserved() {
        let p = CellParams::lg_m50();
        let mut cell = Spm::new(p.clone(), 30);
        let before = cell.solid_lithium();
        let i_app = p.i_ref();
        for _ in 0..200 {
            cell.step(10.0, i_app);
        }
        let after = cell.solid_lithium();
        let rel = (after - before).abs() / before;
        assert!(rel < 1e-9, "solid Li drift {rel}");
    }

    #[test]
    fn faradaic_charge_matches_lithium_transfer() {
        // The volumetric-average concentration drop of the negative electrode
        // equals the charge passed per unit volume: Δc̄ = i_app t /(F ε_s L).
        let p = CellParams::lg_m50();
        let mut cell = Spm::new(p.clone(), 60);
        let before = cell.neg.average();
        let i_app = p.i_ref();
        let (dt, steps) = (5.0, 100);
        for _ in 0..steps {
            cell.step(dt, i_app);
        }
        let delta = before - cell.neg.average();
        let expected = i_app * dt * steps as f64 / (F_FARADAY * p.neg.eps_s * p.neg.thickness);
        let rel = (delta - expected).abs() / expected;
        assert!(rel < 1e-9, "charge mismatch {rel} (Δc̄={delta}, expected={expected})");
    }

    #[test]
    fn low_rate_voltage_approaches_ocv() {
        let p = CellParams::lg_m50();
        let mut cell = Spm::new(p.clone(), 40);
        let i_app = 0.001 * p.i_ref(); // C/1000
        let v = cell.voltage(i_app);
        let ocv = cell.ocv();
        assert!((v - ocv).abs() < 5e-3, "V {v} vs OCV {ocv}");
        let _ = &mut cell;
    }

    #[test]
    fn capacity_drops_with_c_rate() {
        let p = CellParams::lg_m50();
        let cell = Spm::new(p.clone(), 40);
        let (c1, _) = cell.discharge(1.0, 2.5, 10.0, 2000);
        let (c3, _) = cell.discharge(3.0, 2.5, 5.0, 4000);
        assert!(c1 > c3, "capacity should fall with C-rate: {c1} vs {c3}");
        assert!(c1 > 0.0 && c3 > 0.0);
    }

    #[test]
    fn electrolyte_indicator_grows_with_c_rate() {
        let p = CellParams::lg_m50();
        assert!(electrolyte_indicator(1.0, &p) > electrolyte_indicator(0.1, &p));
    }
}
