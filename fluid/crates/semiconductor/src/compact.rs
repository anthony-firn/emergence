//! Compact transistor models folded from TCAD: the circuit rung of the
//! ladder.
//!
//! A table model I_d(V_gs, V_ds) sampled from the drift–diffusion device,
//! factored as I_d = I_sat(V_gs) · s(V_gs, V_ds), with I_sat the current at
//! the largest V_ds and s ∈ [0, 1] the output-characteristic shape.
//! Both are interpolated with monotone piecewise-cubic Hermite polynomials
//! (PCHIP, Fritsch–Carlson): I_sat in z = asinh(I/I₀), where exponential
//! subthreshold current becomes nearly linear in V_gs, and s directly (it
//! is bounded and smooth). Monotonicity in both voltages is preserved (no
//! spurious negative conductance for the circuit solver). Interpolating
//! asinh(I_d) along V_ds instead fails: it jumps from 0 at V_ds = 0 to ~25
//! within one cell.
//!
//! Certificate: the table is checked against fresh TCAD solutions at the
//! centres of its cells (points never used to build it); the worst error
//! there bounds the model's error inside its domain, and outside the
//! domain the model refuses (the circuit unfolds to TCAD). This is the
//! "learned surrogate with a validity certificate" of the emergence ladder,
//! with a table as the simplest learner.

use crate::device::Device;
use crate::mosfet::MosfetSpec;

/// PCHIP slopes (Fritsch–Carlson) for data y on nodes x.
fn pchip_slopes(x: &[f64], y: &[f64]) -> Vec<f64> {
    let n = x.len();
    let h: Vec<f64> = (0..n - 1).map(|i| x[i + 1] - x[i]).collect();
    let del: Vec<f64> = (0..n - 1).map(|i| (y[i + 1] - y[i]) / h[i]).collect();
    let mut d = vec![0.0; n];
    if n == 2 {
        return vec![del[0]; 2];
    }
    for i in 1..n - 1 {
        if del[i - 1] * del[i] <= 0.0 {
            d[i] = 0.0;
        } else {
            let w1 = 2.0 * h[i] + h[i - 1];
            let w2 = h[i] + 2.0 * h[i - 1];
            d[i] = (w1 + w2) / (w1 / del[i - 1] + w2 / del[i]);
        }
    }
    let end = |h0: f64, h1: f64, d0: f64, d1: f64| -> f64 {
        let s = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
        if s * d0 <= 0.0 {
            0.0
        } else if d0 * d1 <= 0.0 && s.abs() > 3.0 * d0.abs() {
            3.0 * d0
        } else {
            s
        }
    };
    d[0] = end(h[0], h[1], del[0], del[1]);
    d[n - 1] = end(h[n - 2], h[n - 3], del[n - 2], del[n - 3]);
    d
}

/// Evaluate the PCHIP interpolant at t.
fn pchip_eval(x: &[f64], y: &[f64], d: &[f64], t: f64) -> f64 {
    let n = x.len();
    let i = match x.iter().position(|&xi| xi > t) {
        Some(0) => 0,
        Some(k) => k - 1,
        None => n - 2,
    };
    let h = x[i + 1] - x[i];
    let s = (t - x[i]) / h;
    let (h00, h10, h01, h11) = (2.0 * s * s * s - 3.0 * s * s + 1.0, s * s * s - 2.0 * s * s + s, -2.0 * s * s * s + 3.0 * s * s, s * s * s - s * s);
    h00 * y[i] + h10 * h * d[i] + h01 * y[i + 1] + h11 * h * d[i + 1]
}

/// Table model of one transistor. Voltages are for an n-type device; a
/// p-type device is stored with both voltages and the current negated
/// (`polarity` −1), so the same monotone table serves both.
#[derive(Clone, Debug)]
pub struct Table {
    pub vgs: Vec<f64>,
    pub vds: Vec<f64>,
    /// z = asinh(I_sat/I₀) per V_gs
    pub zsat: Vec<f64>,
    /// shape I_d/I_sat at [vgs][vds]
    pub shape: Vec<Vec<f64>>,
    pub i0: f64,
    pub polarity: f64,
    /// width scale applied on evaluation (tables are built for the TCAD
    /// width)
    pub scale: f64,
    dzsat: Vec<f64>,
    dshape: Vec<Vec<f64>>,
}

impl Table {
    pub fn new(vgs: Vec<f64>, vds: Vec<f64>, id: Vec<Vec<f64>>, polarity: f64) -> Self {
        let i0 = 1e-15;
        let last = vds.len() - 1;
        let zsat: Vec<f64> = id.iter().map(|row| (polarity * row[last] / i0).asinh()).collect();
        let shape: Vec<Vec<f64>> = id.iter().map(|row| row.iter().map(|&i| i / row[last]).collect()).collect();
        let dzsat = pchip_slopes(&vgs, &zsat);
        let dshape = shape.iter().map(|row| pchip_slopes(&vds, row)).collect();
        Table { vgs, vds, zsat, shape, i0, polarity, scale: 1.0, dzsat, dshape }
    }

    /// Is (vgs, vds) (n-type sense) inside the sampled domain?
    pub fn contains(&self, vgs: f64, vds: f64) -> bool {
        let (g0, g1) = (self.vgs[0], *self.vgs.last().unwrap());
        let (d0, d1) = (self.vds[0], *self.vds.last().unwrap());
        vgs >= g0 - 1e-9 && vgs <= g1 + 1e-9 && vds.abs() >= d0 - 1e-9 && vds.abs() <= d1 + 1e-9
    }

    /// Drain current (A) at terminal voltages in the device's own polarity
    /// (p-type: V_gs, V_ds ≤ 0). V_ds of the wrong sign swaps source and
    /// drain (symmetric device, body effect of the swap neglected).
    pub fn id(&self, vgs: f64, vds: f64) -> f64 {
        let (g, d) = (self.polarity * vgs, self.polarity * vds);
        let (g, d, sg) = if d >= 0.0 { (g, d, 1.0) } else { (g - d, -d, -1.0) };
        let g = g.clamp(self.vgs[0], *self.vgs.last().unwrap());
        let d = d.min(*self.vds.last().unwrap());
        // shape along V_ds in every V_gs row, then along V_gs
        let col: Vec<f64> = (0..self.vgs.len()).map(|k| pchip_eval(&self.vds, &self.shape[k], &self.dshape[k], d)).collect();
        let dcol = pchip_slopes(&self.vgs, &col);
        let sh = pchip_eval(&self.vgs, &col, &dcol, g);
        let isat = self.i0 * pchip_eval(&self.vgs, &self.zsat, &self.dzsat, g).sinh();
        self.polarity * sg * self.scale * isat * sh
    }
}

/// Drain current of a TCAD device (source and body at 0).
pub fn tcad_id(dev: &Device) -> f64 {
    let (a, b) = dev.current(dev.contact("drain"));
    a + b
}

/// Sample I_d on the grid with bias continuation (snake order over V_ds).
/// Voltages in the device's polarity. Returns [vgs][vds].
pub fn sample(dev: &mut Device, vgs: &[f64], vds: &[f64]) -> Vec<Vec<f64>> {
    let (d, g) = (dev.contact("drain"), dev.contact("gate"));
    let mut out = vec![vec![0.0; vds.len()]; vgs.len()];
    for (k, &vg) in vgs.iter().enumerate() {
        dev.ramp(g, vg, 0.1).expect("gate ramp");
        let order: Vec<usize> = if k % 2 == 0 { (0..vds.len()).collect() } else { (0..vds.len()).rev().collect() };
        for j in order {
            dev.ramp(d, vds[j], 0.1).expect("drain ramp");
            out[k][j] = tcad_id(dev);
        }
    }
    out
}

/// Build the table from TCAD and certify it at the cell centres. Returns
/// the table and (worst relative error where |I_d| > 1e-3 of the maximum,
/// worst error in decades below that).
pub fn fold(spec: &MosfetSpec, dev: &mut Device, vgs: &[f64], vds: &[f64]) -> (Table, f64, f64) {
    let pol = if spec.pmos { -1.0 } else { 1.0 };
    let gs: Vec<f64> = vgs.iter().map(|v| pol * v).collect();
    let ds: Vec<f64> = vds.iter().map(|v| pol * v).collect();
    let id = sample(dev, &gs, &ds);
    let table = Table::new(vgs.to_vec(), vds.to_vec(), id.clone(), pol);
    let imax = id.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
    let mid = |v: &[f64]| -> Vec<f64> { v.windows(2).map(|w| 0.5 * (w[0] + w[1])).collect() };
    let (mg, md) = (mid(vgs), mid(vds));
    let gm: Vec<f64> = mg.iter().map(|v| pol * v).collect();
    let dm: Vec<f64> = md.iter().map(|v| pol * v).collect();
    let check = sample(dev, &gm, &dm);
    let (mut rel, mut dec) = (0.0f64, 0.0f64);
    for (k, row) in check.iter().enumerate() {
        for (j, &it) in row.iter().enumerate() {
            let im = table.id(gm[k], dm[j]);
            if it.abs() > 1e-3 * imax {
                rel = rel.max((im / it - 1.0).abs());
            } else if it.abs() > 0.0 && im.abs() > 0.0 {
                dec = dec.max((im / it).abs().log10().abs());
            }
        }
    }
    (table, rel, dec)
}

#[cfg(test)]
mod tests {
    use super::*;
    /// The interpolant reproduces a smooth synthetic transistor law between
    /// its nodes.
    #[test]
    fn table_interpolates_synthetic_law() {
        // EKV all-region law: I = I_s[ln²(1 + e^{v_p/2}) − ln²(1 + e^{(v_p − v_ds)/2})]
        let law = |g: f64, d: f64| {
            let (vt, n) = (0.0259, 1.3);
            let vp = (g - 0.35) / (n * vt);
            let f = |x: f64| (1.0 + (0.5 * x).exp()).ln().powi(2);
            2e-7 * (f(vp) - f(vp - d / vt))
        };
        let vgs: Vec<f64> = (0..=12).map(|i| 0.1 * i as f64).collect();
        let vds = vec![0.0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2];
        let id: Vec<Vec<f64>> = vgs.iter().map(|&g| vds.iter().map(|&d| law(g, d)).collect()).collect();
        let t = Table::new(vgs.clone(), vds.clone(), id, 1.0);
        for (k, &g) in vgs.iter().enumerate() {
            for (j, &d) in vds.iter().enumerate() {
                let m = t.id(g, d);
                let e = law(g, d);
                assert!((m - e).abs() <= 1e-9 * e.abs() + 1e-30, "node ({k},{j}) g {g} d {d}: model {m:e} law {e:e}");
            }
        }
        let mut worst: f64 = 0.0;
        for w in vgs.windows(2) {
            for v in vds.windows(2) {
                let (g, d) = (0.5 * (w[0] + w[1]), 0.5 * (v[0] + v[1]));
                worst = worst.max((t.id(g, d) / law(g, d) - 1.0).abs());
            }
        }
        eprintln!("worst midpoint relative error {worst:.3e}");
        assert!(worst < 0.05);
    }
}

/// Quasi-static terminal charges (gate, drain) on the same grid as the
/// current table, from `Device::qs_charges` (charge-conserving: the four
/// terminal charges sum to zero). Interpolated with PCHIP in both voltages.
#[derive(Clone, Debug)]
pub struct QTable {
    pub vgs: Vec<f64>,
    pub vds: Vec<f64>,
    /// [vgs][vds] → (Q_g, Q_d) in coulombs (TCAD width)
    pub q: Vec<Vec<(f64, f64)>>,
    pub polarity: f64,
    pub scale: f64,
}

impl QTable {
    /// (Q_g, Q_d) at terminal voltages in the device's own polarity
    /// (domain clamped; V_ds < 0 is clamped to 0, a small overshoot range).
    pub fn q(&self, vgs: f64, vds: f64) -> (f64, f64) {
        let (g, d) = (self.polarity * vgs, (self.polarity * vds).max(0.0));
        let g = g.clamp(self.vgs[0], *self.vgs.last().unwrap());
        let d = d.min(*self.vds.last().unwrap());
        let mut out = [0.0; 2];
        for (c, o) in out.iter_mut().enumerate() {
            let col: Vec<f64> = (0..self.vgs.len())
                .map(|k| {
                    let row: Vec<f64> = self.q[k].iter().map(|x| if c == 0 { x.0 } else { x.1 }).collect();
                    let dr = pchip_slopes(&self.vds, &row);
                    pchip_eval(&self.vds, &row, &dr, d)
                })
                .collect();
            let dc = pchip_slopes(&self.vgs, &col);
            *o = self.scale * pchip_eval(&self.vgs, &col, &dc, g);
        }
        (out[0], out[1])
    }
}

/// Sample drain current and quasi-static (Q_g, Q_d) on the grid (device
/// polarity voltages), with bias continuation.
#[allow(clippy::type_complexity)]
pub fn sample_with_charges(dev: &mut Device, vgs: &[f64], vds: &[f64]) -> (Vec<Vec<f64>>, Vec<Vec<(f64, f64)>>) {
    let (d, g) = (dev.contact("drain"), dev.contact("gate"));
    let mut id = vec![vec![0.0; vds.len()]; vgs.len()];
    let mut q = vec![vec![(0.0, 0.0); vds.len()]; vgs.len()];
    for (k, &vg) in vgs.iter().enumerate() {
        dev.ramp(g, vg, 0.1).expect("gate ramp");
        let order: Vec<usize> = if k % 2 == 0 { (0..vds.len()).collect() } else { (0..vds.len()).rev().collect() };
        for j in order {
            dev.ramp(d, vds[j], 0.1).expect("drain ramp");
            id[k][j] = tcad_id(dev);
            let qs = dev.qs_charges();
            q[k][j] = (qs[g], qs[d]);
        }
    }
    (id, q)
}
