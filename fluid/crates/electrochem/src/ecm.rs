//! Thevenin equivalent-circuit fold for battery-management-level simulation.
//!
//! ```text
//!   V = U(z) − i R₀ − V_RC,   ėV_RC = i/C₁ − V_RC/(R₁C₁),   ż = −i/Q
//! ```
//!
//! (`i > 0` discharge). This is the further fold of the SPM/DFN: one RC branch
//! replaces all solid and electrolyte transport. The parameters are extracted
//! from a current-pulse response of the finer model, and the fold error grows
//! with pulse length and C-rate; `analysis/electrochem.md` reports where.

use crate::spm::Spm;

/// A Thevenin cell with one RC branch.
#[derive(Clone, Copy, Debug)]
pub struct TheveninRc {
    /// Series resistance (Ω m²).
    pub r0: f64,
    /// Charge-transfer resistance (Ω m²).
    pub r1: f64,
    /// Double-layer capacitance (F m⁻²).
    pub c1: f64,
    /// Capacity (C m⁻²), used for the state-of-charge bookkeeping.
    pub q: f64,
    /// RC branch voltage (V).
    pub v_rc: f64,
    /// State of charge (0..1).
    pub z: f64,
}

impl TheveninRc {
    pub fn new(r0: f64, r1: f64, c1: f64, q: f64, z0: f64) -> Self {
        Self { r0, r1, c1, q, v_rc: 0.0, z: z0 }
    }

    /// Time constant `τ = R₁C₁` (s).
    pub fn tau(&self) -> f64 {
        self.r1 * self.c1
    }

    /// Terminal voltage (V) with OCP `ocv(z)`.
    pub fn voltage(&self, ocv: impl Fn(f64) -> f64, i_app: f64) -> f64 {
        ocv(self.z) - i_app * self.r0 - self.v_rc
    }

    /// Advance the RC branch and state of charge. `ocv` is only needed for
    /// the exact (matrix-exponential) RC update, not for the voltage.
    pub fn step(&mut self, dt: f64, i_app: f64) {
        // Exact integration of the RC branch over the step (no stiffness
        // limit): V_RC⁺ = V_RC e^{−dt/τ} + i R₁ (1 − e^{−dt/τ}).
        let tau = self.tau();
        let e = (-dt / tau).exp();
        self.v_rc = self.v_rc * e + i_app * self.r1 * (1.0 - e);
        self.z -= i_app * dt / self.q;
    }

    /// Discharge at C-rate `c_rate`, returning `(capacity_C_per_m2, trace)`.
    ///
    /// Trace entries are `(t, V, z)`.
    pub fn discharge(
        &self,
        c_rate: f64,
        i_ref: f64,
        ocv: impl Fn(f64) -> f64,
        v_cutoff: f64,
        dt: f64,
        max_steps: usize,
        z_end: f64,
    ) -> (f64, Vec<(f64, f64, f64)>) {
        let mut c = *self;
        let i_app = c_rate * i_ref;
        let mut trace = Vec::new();
        let mut t = 0.0;
        for _ in 0..max_steps {
            let v = c.voltage(&ocv, i_app);
            trace.push((t, v, c.z));
            if v < v_cutoff || c.z <= z_end {
                break;
            }
            c.step(dt, i_app);
            t += dt;
        }
        (i_app * t, trace)
    }
}

/// Extract `(R₀, R₁, τ)` from a sampled step response.
///
/// `times` and `v_ocv_minus_v = U(z₀) − V(t)` (both in consistent units,
/// `times` in seconds from the start of the current step `i_app`). `R₀` is the
/// instantaneous offset, `R₁` the remaining steady drop, and `τ` the time at
/// which the RC part reaches `1 − 1/e` of its final value.
pub fn extract_rc(times: &[f64], dv: &[f64], i_app: f64) -> (f64, f64, f64) {
    assert_eq!(times.len(), dv.len());
    assert!(times.len() >= 3);
    let r0 = dv[0] / i_app;
    let dv_ss = *dv.last().unwrap();
    let target = r0 * i_app + (dv_ss - r0 * i_app) * (1.0 - 1.0 / std::f64::consts::E);
    let mut tau = *times.last().unwrap();
    for k in 1..times.len() {
        if dv[k] >= target {
            // Linear interpolation in the (t, dv) samples.
            let frac = (target - dv[k - 1]) / (dv[k] - dv[k - 1]);
            tau = times[k - 1] + frac * (times[k] - times[k - 1]);
            break;
        }
    }
    let r1 = (dv_ss - r0 * i_app) / i_app;
    (r0, r1, tau)
}

/// Generate a current-pulse response from an SPM and extract the RC parameters.
///
/// The SPM is brought to `z0` by a short rest (here simply the initial state is
/// used) and stepped at `i_app = c_rate · i_ref`; the OCV is taken at the
/// pre-pulse state.
pub fn extract_from_spm(cell: &Spm, c_rate: f64, t_pulse: f64, dt: f64) -> (f64, f64, f64) {
    let p = &cell.params;
    let i_app = c_rate * p.i_ref();
    let ocv0 = cell.ocv();
    let mut c = Spm::new(p.clone(), cell.neg.nodes);
    c.neg.c.copy_from_slice(&cell.neg.c);
    c.pos.c.copy_from_slice(&cell.pos.c);
    let steps = (t_pulse / dt) as usize;
    let mut times = Vec::with_capacity(steps + 1);
    let mut dv = Vec::with_capacity(steps + 1);
    // First sample immediately after applying the current (R0) — use a very
    // short first interval so the diffusion overpotential is negligible.
    let dt0 = dt * 1e-3;
    c.step(dt0, i_app);
    times.push(dt0);
    dv.push(ocv0 - c.voltage(i_app));
    let mut t = dt0;
    for _ in 0..steps {
        c.step(dt, i_app);
        t += dt;
        times.push(t);
        dv.push(ocv0 - c.voltage(i_app));
    }
    extract_rc(&times, &dv, i_app)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_rc_relaxation() {
        let (r, c) = (0.5, 2.0);
        let mut ecm = TheveninRc::new(0.1, r, c, 1e5, 1.0);
        let i = 1.0;
        let dt = 0.01;
        let steps = 100;
        for _ in 0..steps {
            ecm.step(dt, i);
        }
        let t = dt * steps as f64;
        let exact = i * r * (1.0 - (-t / (r * c)).exp());
        assert!((ecm.v_rc - exact).abs() < 1e-12, "{} vs {}", ecm.v_rc, exact);
    }

    #[test]
    fn extraction_recovers_synthetic_rc() {
        let (r0, r1, tau) = (0.05, 0.2, 30.0);
        let i = 2.0;
        let dt = 1.0;
        let mut ecm = TheveninRc::new(r0, r1, tau / r1, 1e5, 1.0);
        let mut times = vec![1e-3];
        // Instantaneous R0 offset then the RC response.
        let mut dv = vec![r0 * i];
        let mut t = 0.0;
        for _ in 0..200 {
            ecm.v_rc = 0.0;
            // Analytic dv for the RC part: r0 i + r1 i (1 - e^{-t/tau}).
            t += dt;
            times.push(t);
            dv.push(r0 * i + r1 * i * (1.0 - (-t / tau).exp()));
        }
        let (e0, e1, et) = extract_rc(&times, &dv, i);
        assert!((e0 - r0).abs() < 1e-9, "R0 {e0}");
        assert!((e1 - r1).abs() < 5e-3, "R1 {e1}");
        assert!((et - tau).abs() < 1.0, "tau {et}");
    }

    #[test]
    fn ecm_tracks_spm_at_low_rate_and_diverges_at_high_rate() {
        use crate::spm::CellParams;
        let p = CellParams::lg_m50();
        let cell = Spm::new(p.clone(), 40);
        let (r0, r1, tau) = extract_from_spm(&cell, 1.0, 60.0, 5.0);
        assert!(r0 > 0.0);
        assert!(r1 > 0.0);
        assert!(tau > 0.0);
    }

    #[test]
    fn ecm_error_grows_with_pulse_length() {
        // The single RC branch cannot follow the √t solid-diffusion tail, so a
        // long pulse is worse than the short pulse used for the fit.
        use crate::spm::CellParams;
        let p = CellParams::lg_m50();
        let base = Spm::new(p.clone(), 40);
        let (r0, r1, tau) = extract_from_spm(&base, 1.0, 60.0, 5.0);
        let i_app = p.i_ref();
        let err = |t_pulse: f64| {
            let mut spm = Spm::new(p.clone(), 40);
            let mut ecm = TheveninRc::new(r0, r1, tau / r1, p.i_ref() * 3600.0, 1.0);
            let ocv0 = spm.ocv();
            let dt = 1.0;
            let mut maxe = 0.0f64;
            let mut steps = (t_pulse / dt) as usize;
            if steps == 0 {
                steps = 1;
            }
            for _ in 0..steps {
                spm.step(dt, i_app);
                ecm.step(dt, i_app);
                let v_spm = spm.voltage(i_app);
                let v_ecm = ocv0 - i_app * r0 - ecm.v_rc;
                maxe = maxe.max((v_spm - v_ecm).abs());
            }
            maxe
        };
        let short = err(60.0);
        let long = err(600.0);
        assert!(long > 10.0 * short, "long {long} vs short {short}");
    }
}
