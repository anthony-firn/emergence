//! Electrolyte-enhanced single-particle model (SPMe): the fold of the SPM that
//! restores the *electrolyte* concentration.
//!
//! The SPM assumes the electrolyte is uniform. SPMe keeps a 1D electrolyte
//! concentration across negative electrode | separator | positive electrode
//! with the reaction current smeared uniformly (so the solid and electrolyte
//! problems stay one-dimensional in the radial and x directions respectively):
//!
//! ```text
//!   ε_e ∂c_e/∂t = ∂_x(D_e^eff ∂_x c_e) + (1−t₊) a_s i_loc / F
//!   V_SPMe = V_SPM + (2RT/F)(1−t₊) ln(c_e⁺/c_e⁻)
//! ```
//!
//! References for the DFN/SPMe reduction: M. Doyle, T. Fuller, J. Newman,
//! *J. Electrochem. Soc.* **140**, 1526 (1993); Marquis et al.,
//! *J. Electrochem. Soc.* **166**, A3693 (2019) (SPMe error analysis).
//!
//! The true reference is the pseudo-2D DFN, which resolves the **position
//! dependence** of the reaction; SPMe is the intermediate rung. See
//! `analysis/electrochem.md`.

use crate::spm::{CellParams, Spm};
use crate::{F_FARADAY, R_GAS};

/// 1D electrolyte concentration on `[0, length]` with zero-flux boundaries and
/// a per-node volumetric source (Crank–Nicolson + Thomas).
pub struct Electrolyte1D {
    pub n: usize,
    pub length: f64,
    /// x coordinates of the face-centred nodes.
    pub x: Vec<f64>,
    /// Electrolyte concentration (mol m⁻³).
    pub c: Vec<f64>,
    dx: f64,
    /// Effective diffusivity `D_e ε_e^{b−1}` (m² s⁻¹).
    d_eff: f64,
}

impl Electrolyte1D {
    pub fn new(n: usize, length: f64, c0: f64, d_eff: f64) -> Self {
        let dx = length / n as f64;
        let x = (0..n).map(|i| (i as f64 + 0.5) * dx).collect();
        Self { n, length, x, c: vec![c0; n], dx, d_eff }
    }

    /// Concentration at the node nearest `x0`.
    pub fn at(&self, x0: f64) -> f64 {
        let i = ((x0 / self.dx) - 0.5).round().clamp(0.0, (self.n - 1) as f64) as usize;
        self.c[i]
    }

    /// Advance one step with source `s` (mol m⁻³ s⁻¹) per node.
    pub fn step(&mut self, dt: f64, s: &[f64]) {
        assert_eq!(s.len(), self.n);
        let r = self.d_eff * dt / (self.dx * self.dx);
        let n = self.n;
        let mut a = vec![0.0; n];
        let mut b = vec![0.0; n];
        let mut cc = vec![0.0; n];
        let mut rhs = vec![0.0; n];
        for i in 0..n {
            // No-flux ends: the missing neighbour mirrors the interior node,
            // so the end rows get the full off-diagonal coefficient −r.
            if i == 0 {
                a[i] = 0.0;
                b[i] = 1.0 + r;
                cc[i] = -r;
            } else if i == n - 1 {
                a[i] = -r;
                b[i] = 1.0 + r;
                cc[i] = 0.0;
            } else {
                a[i] = -0.5 * r;
                b[i] = 1.0 + r;
                cc[i] = -0.5 * r;
            }
            let left = if i == 0 { self.c[1] } else { self.c[i - 1] };
            let right = if i == n - 1 { self.c[i - 1] } else { self.c[i + 1] };
            rhs[i] = 0.5 * r * left + (1.0 - r) * self.c[i] + 0.5 * r * right + dt * s[i];
        }
        self.c = thomas(&mut a, &mut b, &mut cc, &mut rhs);
    }

    /// Total moles per unit area (mol m⁻²).
    pub fn total(&self) -> f64 {
        self.c.iter().sum::<f64>() * self.dx
    }
}

fn thomas(a: &mut [f64], b: &mut [f64], cc: &mut [f64], d: &mut [f64]) -> Vec<f64> {
    let n = b.len();
    for i in 1..n {
        let m = a[i] / b[i - 1];
        b[i] -= m * cc[i - 1];
        d[i] -= m * d[i - 1];
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d[n - 1] / b[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (d[i] - cc[i] * x[i + 1]) / b[i];
    }
    x
}

/// SPMe cell: SPM + electrolyte.
pub struct Spme {
    pub spm: Spm,
    pub elec: Electrolyte1D,
    /// Transference number `t₊`.
    pub t_plus: f64,
    x_sep_n: f64,
    x_sep_p: f64,
    l_n: f64,
    l_p: f64,
}

impl Spme {
    pub fn new(params: CellParams, solid_nodes: usize, elec_nodes: usize, t_plus: f64) -> Self {
        let l_n = params.neg.thickness;
        let l_sep = 12.0e-6; // ESTIMATE separator thickness
        let l_p = params.pos.thickness;
        let length = l_n + l_sep + l_p;
        let eps_e: f64 = 0.4; // ESTIMATE electrode porosity
        let bruggeman: f64 = 1.5; // ESTIMATE
        let d_e = 2.0e-10; // ESTIMATE
        let d_eff = d_e * eps_e.powf(bruggeman - 1.0);
        let elec = Electrolyte1D::new(elec_nodes, length, params.c_e, d_eff);
        Self {
            spm: Spm::new(params, solid_nodes),
            elec,
            t_plus,
            x_sep_n: l_n,
            x_sep_p: l_n + l_sep,
            l_n,
            l_p,
        }
    }

    /// Volumetric electrolyte source (mol m⁻³ s⁻¹) for applied current density
    /// `i_app` (A m⁻², discharge positive).
    pub fn source(&self, i_app: f64) -> Vec<f64> {
        let p = &self.spm.params;
        let amp = (1.0 - self.t_plus) * i_app / F_FARADAY;
        let mut raw: Vec<f64> = self
            .elec
            .x
            .iter()
            .map(|&x| {
                if x < self.l_n {
                    amp / p.neg.thickness
                } else if x >= self.x_sep_p {
                    -amp / p.pos.thickness
                } else {
                    0.0
                }
            })
            .collect();
        // The exact source integrates to zero (produced Li+ = consumed Li+).
        // Node/region misalignment leaves a tiny discrete residual; remove it
        // so the electrolyte inventory is conserved exactly.
        let mean = raw.iter().sum::<f64>() * self.elec.dx / self.elec.length;
        for s in raw.iter_mut() {
            *s -= mean;
        }
        raw
    }

    /// Terminal voltage (V): SPM voltage plus the electrolyte concentration
    /// overpotential `(2RT/F)(1−t₊) ln(c_e⁺/c_e⁻)`.
    pub fn voltage(&self, i_app: f64) -> f64 {
        let base = self.spm.voltage(i_app);
        let ce_n = self.elec.at(self.x_sep_n);
        let ce_p = self.elec.at(self.x_sep_p);
        let eta_c = 2.0 * R_GAS * self.spm.params.t / F_FARADAY * (1.0 - self.t_plus) * (ce_p / ce_n).ln();
        base + eta_c
    }

    pub fn step(&mut self, dt: f64, i_app: f64) {
        self.spm.step(dt, i_app);
        let s = self.source(i_app);
        self.elec.step(dt, &s);
    }

    /// Constant-current discharge to cutoff. Returns `(capacity_mAh_m2, trace)`
    /// with trace `(t, V_spme, V_spm, ce_n, ce_p)`.
    pub fn discharge(
        &self,
        c_rate: f64,
        v_cutoff: f64,
        dt: f64,
        max_steps: usize,
    ) -> (f64, Vec<(f64, f64, f64, f64, f64)>) {
        let p = self.spm.params.clone();
        let mut spm = Spm::new(p.clone(), self.spm.neg.nodes);
        spm.neg.c.copy_from_slice(&self.spm.neg.c);
        spm.pos.c.copy_from_slice(&self.spm.pos.c);
        let mut elec = Electrolyte1D::new(self.elec.n, self.elec.length, p.c_e, self.elec.d_eff);
        let t_plus = self.t_plus;
        let (x_sep_n, x_sep_p, l_n, l_p) = (self.x_sep_n, self.x_sep_p, self.l_n, self.l_p);
        let i_app = c_rate * p.i_ref();
        let mut trace = Vec::new();
        let mut t = 0.0;
        for _ in 0..max_steps {
            let v_spm = spm.voltage(i_app);
            let ce_n = elec.at(x_sep_n);
            let ce_p = elec.at(x_sep_p);
            let eta_c = 2.0 * R_GAS * p.t / F_FARADAY * (1.0 - t_plus) * (ce_p / ce_n).ln();
            let v = v_spm + eta_c;
            trace.push((t, v, v_spm, ce_n, ce_p));
            if v < v_cutoff {
                break;
            }
            spm.step(dt, i_app);
            let amp = (1.0 - t_plus) * i_app / F_FARADAY;
            let mut s: Vec<f64> = elec
                .x
                .iter()
                .map(|&x| {
                    if x < l_n {
                        amp / l_n
                    } else if x >= x_sep_p {
                        -amp / l_p
                    } else {
                        0.0
                    }
                })
                .collect();
            let mean = s.iter().sum::<f64>() * elec.dx / elec.length;
            for si in s.iter_mut() {
                *si -= mean;
            }
            elec.step(dt, &s);
            t += dt;
        }
        (i_app * t / 3600.0 * 1000.0, trace)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn electrolyte_source_is_charge_neutral() {
        // Produced Li+ in the negative equals consumed Li+ in the positive, so
        // the total electrolyte inventory is conserved.
        let p = CellParams::lg_m50();
        let s = Spme::new(p, 30, 60, 0.2594);
        let total: f64 = s.source(50.0).iter().sum::<f64>() * s.elec.dx;
        assert!(total.abs() < 1e-9, "net source {total}");
    }

    #[test]
    fn zero_current_leaves_electrolyte_uniform() {
        let p = CellParams::lg_m50();
        let mut s = Spme::new(p, 20, 40, 0.2594);
        let src = s.source(0.0);
        for _ in 0..50 {
            s.elec.step(10.0, &src);
        }
        let c0 = s.elec.c[0];
        for &c in &s.elec.c {
            assert!((c - c0).abs() < 1e-9);
        }
    }

    #[test]
    fn discharge_depletes_positive_electrolyte() {
        // On discharge Li+ is produced at the negative and consumed at the
        // positive, so c_e is higher at the negative/separator interface.
        let p = CellParams::lg_m50();
        let cell = Spme::new(p, 30, 60, 0.2594);
        let (_, trace) = cell.discharge(1.0, 2.5, 10.0, 2000);
        let (_, _, _, ce_n, ce_p) = *trace.last().unwrap();
        assert!(ce_n > ce_p, "ce_n {ce_n} ce_p {ce_p}");
    }

    #[test]
    fn spm_error_grows_with_c_rate() {
        let p = CellParams::lg_m50();
        let cell = Spme::new(p, 30, 80, 0.2594);
        let max_diff = |c_rate: f64| {
            let (_, trace) = cell.discharge(c_rate, 2.5, 5.0, 4000);
            trace
                .iter()
                .map(|&(_, v, vspm, _, _)| (v - vspm).abs())
                .fold(0.0f64, f64::max)
        };
        let d1 = max_diff(1.0);
        let d5 = max_diff(5.0);
        assert!(d5 > d1, "electrolyte correction should grow: {d1} vs {d5}");
    }

    #[test]
    fn spme_voltage_is_close_to_spm_at_low_rate() {
        let p = CellParams::lg_m50();
        let cell = Spme::new(p, 30, 80, 0.2594);
        let (_, trace) = cell.discharge(0.2, 2.5, 20.0, 4000);
        let max_diff = trace
            .iter()
            .map(|&(_, v, vspm, _, _)| (v - vspm).abs())
            .fold(0.0f64, f64::max);
        assert!(max_diff < 5e-3, "low-rate correction too large: {max_diff}");
    }
}
