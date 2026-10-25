//! CMOS circuits on compact models: inverter transfer curve, noise margins
//! (the indicator that the logic rung above is valid), ring oscillator.

/// A transistor model for circuits: drain current (into the drain) from
/// V_gs, V_ds in the device's own polarity.
pub trait Fet: Sync {
    fn id(&self, vgs: f64, vds: f64) -> f64;
}

impl Fet for crate::compact::Table {
    fn id(&self, vgs: f64, vds: f64) -> f64 {
        crate::compact::Table::id(self, vgs, vds)
    }
}

/// Net current into the output node of an inverter (pMOS from V_dd, nMOS
/// to ground): −I_d,p − I_d,n.
pub fn inverter_current(n: &dyn Fet, p: &dyn Fet, vdd: f64, vin: f64, vout: f64) -> f64 {
    -p.id(vin - vdd, vout - vdd) - n.id(vin, vout)
}

/// Inverter output at input vin (bisection: the node current falls
/// monotonically with V_out).
pub fn inverter_vout(n: &dyn Fet, p: &dyn Fet, vdd: f64, vin: f64) -> f64 {
    let (mut lo, mut hi) = (0.0, vdd);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if inverter_current(n, p, vdd, vin, mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Voltage transfer curve characteristics: switching threshold V_M and the
/// noise margins NM_L = V_IL − V_OL, NM_H = V_OH − V_IH (unity-gain points).
#[derive(Clone, Copy, Debug)]
pub struct Vtc {
    pub vm: f64,
    pub vil: f64,
    pub vih: f64,
    pub vol: f64,
    pub voh: f64,
    pub nml: f64,
    pub nmh: f64,
    pub gain: f64,
}

pub fn vtc(n: &dyn Fet, p: &dyn Fet, vdd: f64, points: usize) -> (Vec<(f64, f64)>, Vtc) {
    let curve: Vec<(f64, f64)> = (0..=points).map(|i| {
        let vin = vdd * i as f64 / points as f64;
        (vin, inverter_vout(n, p, vdd, vin))
    }).collect();
    let slope: Vec<(f64, f64, f64)> = curve.windows(2).map(|w| (0.5 * (w[0].0 + w[1].0), 0.5 * (w[0].1 + w[1].1), (w[1].1 - w[0].1) / (w[1].0 - w[0].0))).collect();
    let cross = |from_left: bool| -> (f64, f64) {
        let it: Box<dyn Iterator<Item = &(f64, f64, f64)>> = if from_left { Box::new(slope.iter()) } else { Box::new(slope.iter().rev()) };
        let mut prev: Option<&(f64, f64, f64)> = None;
        for s in it {
            if s.2 <= -1.0 {
                if let Some(p) = prev {
                    let t = (-1.0 - p.2) / (s.2 - p.2);
                    return (p.0 + t * (s.0 - p.0), p.1 + t * (s.1 - p.1));
                }
                return (s.0, s.1);
            }
            prev = Some(s);
        }
        (f64::NAN, f64::NAN)
    };
    let (vil, voh) = cross(true);
    let (vih, vol) = cross(false);
    let mut vm = f64::NAN;
    for w in curve.windows(2) {
        let (a, b) = (w[0].1 - w[0].0, w[1].1 - w[1].0);
        if a >= 0.0 && b < 0.0 {
            vm = w[0].0 + a / (a - b) * (w[1].0 - w[0].0);
        }
    }
    let gain = slope.iter().map(|s| -s.2).fold(0.0, f64::max);
    (curve, Vtc { vm, vil, vih, vol, voh, nml: vil - vol, nmh: voh - vih, gain })
}

/// Ring oscillator of `stages` (odd) inverters with node capacitance c (F):
/// RK4 from a perturbed state; returns the period (s) measured on node 0
/// after `settle` periods, and the propagation delay per stage.
pub fn ring_oscillator(n: &dyn Fet, p: &dyn Fet, vdd: f64, stages: usize, c: f64, dt: f64, t_end: f64) -> (f64, f64, Vec<(f64, f64)>) {
    let mut v: Vec<f64> = (0..stages).map(|i| if i % 2 == 0 { 0.0 } else { vdd }).collect();
    v[0] = 0.5 * vdd;
    let f = |v: &[f64]| -> Vec<f64> { (0..stages).map(|i| inverter_current(n, p, vdd, v[(i + stages - 1) % stages], v[i]) / c).collect() };
    let mut t = 0.0;
    let mut crossings = Vec::new();
    let mut trace = Vec::new();
    let mut prev = v[0];
    while t < t_end {
        let k1 = f(&v);
        let v2: Vec<f64> = v.iter().zip(&k1).map(|(a, k)| a + 0.5 * dt * k).collect();
        let k2 = f(&v2);
        let v3: Vec<f64> = v.iter().zip(&k2).map(|(a, k)| a + 0.5 * dt * k).collect();
        let k3 = f(&v3);
        let v4: Vec<f64> = v.iter().zip(&k3).map(|(a, k)| a + dt * k).collect();
        let k4 = f(&v4);
        for i in 0..stages {
            v[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        t += dt;
        // rising crossings of V_dd/2 on node 0
        if prev < 0.5 * vdd && v[0] >= 0.5 * vdd {
            let tc = t - dt * (v[0] - 0.5 * vdd) / (v[0] - prev);
            crossings.push(tc);
        }
        prev = v[0];
        trace.push((t, v[0]));
    }
    let k = crossings.len();
    let period = if k >= 4 { (crossings[k - 1] - crossings[k - 3]) / 2.0 } else { f64::NAN };
    (period, period / (2.0 * stages as f64), trace)
}

/// Quasi-static transient of an inverter driving a load capacitance:
/// d/dt[Q_d,n + Q_d,p + C_L V_out] = −I_d,n − I_d,p, with the drain charges
/// depending on (V_in, V_out). Variable-step BDF2 in the node charge on the
/// time grid, Newton on V_out (numerical derivative). Returns (t, V_in,
/// V_out).
#[allow(clippy::too_many_arguments)]
pub fn inverter_transient_qs(
    n: &dyn Fet,
    nq: &crate::compact::QTable,
    p: &dyn Fet,
    pq: &crate::compact::QTable,
    vdd: f64,
    cl: f64,
    vin: &dyn Fn(f64) -> f64,
    times: &[f64],
) -> Vec<(f64, f64, f64)> {
    let charge = |vi: f64, vo: f64| nq.q(vi, vo).1 + pq.q(vi - vdd, vo - vdd).1 + cl * vo;
    let current = |vi: f64, vo: f64| -p.id(vi - vdd, vo - vdd) - n.id(vi, vo);
    // start from the DC solution at the first time
    let mut vo = inverter_vout(n, p, vdd, vin(times[0]));
    let mut q_hist = vec![charge(vin(times[0]), vo)];
    let mut out = vec![(times[0], vin(times[0]), vo)];
    let mut dt_prev = 0.0;
    for k in 1..times.len() {
        let (t, dt) = (times[k], times[k] - times[k - 1]);
        let vi = vin(t);
        let a: [f64; 3] = if q_hist.len() >= 2 {
            let w = dt / dt_prev;
            [(1.0 + 2.0 * w) / (1.0 + w), -(1.0 + w), w * w / (1.0 + w)]
        } else {
            [1.0, -1.0, 0.0]
        };
        dt_prev = dt;
        let hist: f64 = q_hist.iter().enumerate().map(|(l, q)| a[l + 1] * q).sum();
        let resid = |v: f64| (a[0] * charge(vi, v) + hist) / dt - current(vi, v);
        for _ in 0..50 {
            let r = resid(vo);
            let h = 1e-6;
            let dr = (resid(vo + h) - resid(vo - h)) / (2.0 * h);
            let step = r / dr;
            vo -= step.clamp(-0.2, 0.2);
            if step.abs() < 1e-12 {
                break;
            }
        }
        q_hist.insert(0, charge(vi, vo));
        q_hist.truncate(2);
        out.push((t, vi, vo));
    }
    out
}

/// 50 % crossing time of a waveform (falling if `falling`).
pub fn crossing(w: &[(f64, f64, f64)], level: f64, falling: bool, input: bool) -> f64 {
    for p in w.windows(2) {
        let (a, b) = if input { (p[0].1, p[1].1) } else { (p[0].2, p[1].2) };
        let hit = if falling { a > level && b <= level } else { a < level && b >= level };
        if hit {
            return p[0].0 + (a - level) / (a - b) * (p[1].0 - p[0].0);
        }
    }
    f64::NAN
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compact::{QTable, Table};
    /// Cost of the quasi-static table transient at production table size
    /// (13 × 16, EKV-law currents and smooth charges).
    #[test]
    fn qs_transient_timing() {
        let vgs: Vec<f64> = (0..=12).map(|i| 0.1 * i as f64).collect();
        let vds = vec![0.0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2];
        let law = |g: f64, d: f64| {
            let vp = (g - 0.35) / (1.3 * 0.0259);
            let f = |x: f64| (1.0 + (0.5 * x).exp()).ln().powi(2);
            2e-6 * (f(vp) - f(vp - d / 0.0259))
        };
        let id: Vec<Vec<f64>> = vgs.iter().map(|&g| vds.iter().map(|&d| law(g, d)).collect()).collect();
        let idp: Vec<Vec<f64>> = id.iter().map(|r| r.iter().map(|x| -x).collect()).collect();
        let q: Vec<Vec<(f64, f64)>> = vgs.iter().map(|&g| vds.iter().map(|&d| (1e-15 * g, -0.3e-15 * g - 0.2e-15 * d)).collect()).collect();
        let tn = Table::new(vgs.clone(), vds.clone(), id, 1.0);
        let mut tp = Table::new(vgs.clone(), vds.clone(), idp, -1.0);
        tp.scale = 1.0;
        let qn = QTable { vgs: vgs.clone(), vds: vds.clone(), q: q.clone(), polarity: 1.0, scale: 1.0 };
        let qp = QTable { vgs, vds, q, polarity: -1.0, scale: 1.0 };
        let times: Vec<f64> = (0..=175).map(|k| k as f64 * 0.5e-12).collect();
        let vin = |t: f64| if t < 2e-12 { 0.0 } else if t < 52e-12 { 1.2 * (t - 2e-12) / 50e-12 } else { 1.2 };
        let t0 = std::time::Instant::now();
        let w = inverter_transient_qs(&tn, &qn, &tp, &qp, 1.2, 2e-15, &vin, &times);
        eprintln!("QS table transient: {} steps in {:.4} s", w.len() - 1, t0.elapsed().as_secs_f64());
        assert!(w.last().unwrap().2 < 0.1);
    }
}
