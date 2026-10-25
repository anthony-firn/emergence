//! Heat: the thermal rungs below and above the device.
//!
//! - `DeviceHeat`: steady 2D conduction on the device mesh (box method,
//!   k_Si(T) by Picard iteration, oxide k constant), sources from
//!   `Device::heat`, the bottom of the silicon held at the die temperature,
//!   other boundaries adiabatic. Fold: the device is simulated isothermal at
//!   its mean temperature while the temperature spread inside it stays
//!   small (the indicator); otherwise a non-isothermal drift–diffusion model
//!   would be needed.
//! - `graetz`: the coolant. Laminar flow between parallel plates heated at
//!   constant flux, u ∂T/∂x = α ∂²T/∂y² (Pe ≫ 1: axial conduction folded
//!   away), marched with variable-step BDF2. Gives the local Nusselt number from the
//!   entrance to the fully developed limit Nu = 8.235 — the downstream fold.
//! - `runaway`: chip energy balance T = T_c + R_th (P_dyn + P_leak(T)):
//!   fixed points, their stability, and thermal runaway where none exist.

use crate::band::Band;
use crate::device::Device;
use crate::mesh::Material;
use crate::silicon::{k_si, K_OX};

/// Steady heat conduction on the device mesh.
pub struct DeviceHeat {
    /// node temperatures (K)
    pub t: Vec<f64>,
    pub picard_iterations: usize,
}

impl DeviceHeat {
    /// Solve for the temperature with sources `q` (W per node, per device
    /// depth) and the bottom silicon row held at `t_die`.
    pub fn solve(dev: &Device, q: &[f64], t_die: f64) -> Self {
        let m = &dev.mesh;
        let (nx, ny) = (m.nx(), m.ny());
        let nn = m.nodes();
        let mut t = vec![t_die; nn];
        let mut its = 0;
        for it in 0..50 {
            its = it + 1;
            let mut a = Band::new(nn, ny, ny);
            let mut b = vec![0.0; nn];
            // element conductances (k at the element's mean temperature)
            for i in 0..nx - 1 {
                for j in 0..ny - 1 {
                    let (hx, hy) = (m.x[i + 1] - m.x[i], m.y[j + 1] - m.y[j]);
                    let ks: [usize; 4] = [m.node(i, j), m.node(i + 1, j), m.node(i, j + 1), m.node(i + 1, j + 1)];
                    let tm = 0.25 * ks.iter().map(|&k| t[k]).sum::<f64>();
                    let k = if m.material(i, j) == Material::Silicon { k_si(tm) } else { K_OX };
                    // each element contributes half its height to its two x edges
                    // and half its width to its two y edges
                    let gx = k * 0.5 * hy / hx;
                    let gy = k * 0.5 * hx / hy;
                    for (p, q2, g) in [(ks[0], ks[1], gx), (ks[2], ks[3], gx), (ks[0], ks[2], gy), (ks[1], ks[3], gy)] {
                        a.add(p, p, g);
                        a.add(p, q2, -g);
                        a.add(q2, q2, g);
                        a.add(q2, p, -g);
                    }
                }
            }
            for k in 0..nn {
                b[k] = q[k] / dev.depth;
            }
            // Dirichlet: bottom row
            for i in 0..nx {
                let k = m.node(i, ny - 1);
                for c in k.saturating_sub(ny)..=(k + ny).min(nn - 1) {
                    if a.get(k, c) != 0.0 {
                        a.add(k, c, -a.get(k, c));
                    }
                }
                a.add(k, k, 1.0);
                b[k] = t_die;
            }
            a.solve(&mut b).expect("heat solve");
            let change = b.iter().zip(&t).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max);
            t = b;
            if change < 1e-9 {
                break;
            }
        }
        DeviceHeat { t, picard_iterations: its }
    }

    /// Heat-weighted mean temperature and the spread (max − min) over the
    /// nodes carrying 99% of the heat.
    pub fn summary(&self, q: &[f64]) -> (f64, f64, f64) {
        let total: f64 = q.iter().map(|x| x.max(0.0)).sum();
        let mean = q.iter().zip(&self.t).map(|(qi, ti)| qi.max(0.0) * ti).sum::<f64>() / total.max(1e-300);
        let mut idx: Vec<usize> = (0..q.len()).collect();
        idx.sort_by(|&a, &b| q[b].partial_cmp(&q[a]).unwrap());
        let (mut acc, mut lo, mut hi) = (0.0, f64::INFINITY, f64::NEG_INFINITY);
        for k in idx {
            if acc > 0.99 * total {
                break;
            }
            acc += q[k].max(0.0);
            lo = lo.min(self.t[k]);
            hi = hi.max(self.t[k]);
        }
        let tmax = self.t.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (mean, hi - lo, tmax)
    }
}

/// Graetz problem, parallel plates (gap 2H), both walls at heat flux q'',
/// Poiseuille flow of mean velocity U: march u ∂T/∂x = α ∂²T/∂y² from a
/// uniform inlet temperature. Returns (x* = x/(D_h Re Pr), Nu_x) with
/// D_h = 4H, Nu = q'' D_h / (k (T_w − T_b)); dimensionless (only x* and Nu
/// matter), solved on y ∈ [0, H] by symmetry with n cells.
pub fn graetz(n: usize, x_end: f64, steps: usize) -> Vec<(f64, f64)> {
    // nondimensional: y ∈ [0, 1] (wall at 1), u = 1.5 (1 − y²) (mean 1),
    // θ = (T − T_in) k/(q'' H); ∂θ/∂ξ = (1/u) ∂²θ/∂y², ξ = α x/(U H²);
    // wall flux ∂θ/∂y = 1 at y = 1, symmetry at y = 0.
    // x* = x/(D_h Re Pr) = α x/(U D_h²) = ξ/16.
    let h = 1.0 / n as f64;
    let y: Vec<f64> = (0..=n).map(|i| i as f64 * h).collect();
    let u: Vec<f64> = y.iter().map(|&yy| 1.5 * (1.0 - yy * yy)).collect();
    let mut th = vec![0.0; n + 1];
    let xi_end = 16.0 * x_end;
    let mut out = Vec::new();
    // geometric steps in ξ (resolve the entrance)
    let r = (xi_end / 1e-7f64).powf(1.0 / steps as f64);
    let mut xi = 0.0;
    let mut dxi = 1e-7 * (r - 1.0).max(1e-12);
    // finite volumes of width h (half at the ends); wall node volume has
    // u = 0: a capacity-free node there, handled by the flux condition
    let vol: Vec<f64> = (0..=n).map(|i| if i == 0 || i == n { 0.5 * h } else { h }).collect();
    // variable-step BDF2 (L-stable: Crank–Nicolson rings on the stiff modes
    // once the steps grow), BDF1 for the first step
    let mut th_old: Option<Vec<f64>> = None;
    let mut dxi_old = 0.0;
    for _ in 0..steps {
        let (a0, a1, a2) = match &th_old {
            Some(_) => {
                let w = dxi / dxi_old;
                ((1.0 + 2.0 * w) / (1.0 + w), -(1.0 + w), w * w / (1.0 + w))
            }
            None => (1.0, -1.0, 0.0),
        };
        let mut a = Band::new(n + 1, 1, 1);
        let mut rhs = vec![0.0; n + 1];
        for i in 0..=n {
            let cap = u[i] * vol[i] / dxi;
            let flux = if i == n { 1.0 } else { 0.0 };
            let hist = a1 * th[i] + th_old.as_ref().map_or(0.0, |o| a2 * o[i]);
            rhs[i] = -cap * hist + flux;
            a.add(i, i, cap * a0);
            if i > 0 {
                a.add(i, i, 1.0 / h);
                a.add(i, i - 1, -1.0 / h);
            }
            if i < n {
                a.add(i, i, 1.0 / h);
                a.add(i, i + 1, -1.0 / h);
            }
        }
        a.solve(&mut rhs).expect("graetz");
        th_old = Some(std::mem::replace(&mut th, rhs));
        dxi_old = dxi;
        xi += dxi;
        dxi *= r;
        // bulk temperature: ∫u θ dy / ∫u dy (trapezoid on the same volumes)
        let num: f64 = (0..=n).map(|i| u[i] * th[i] * vol[i]).sum();
        let den: f64 = (0..=n).map(|i| u[i] * vol[i]).sum();
        let tb = num / den;
        // Nu = q'' D_h/(k (T_w − T_b)) = 4/(θ_w − θ_b)
        out.push((xi / 16.0, 4.0 / (th[n] - tb)));
    }
    out
}

/// Fixed points of T = T_c + R_th (P_dyn + P_leak(T)) on [t_lo, t_hi]:
/// (T, stable) where stable means d/dT [R_th P_leak] < 1.
pub fn runaway(t_c: f64, r_th: f64, p_dyn: f64, p_leak: &dyn Fn(f64) -> f64, t_lo: f64, t_hi: f64) -> Vec<(f64, bool)> {
    let f = |t: f64| t_c + r_th * (p_dyn + p_leak(t)) - t;
    let n = 4000;
    let mut out = Vec::new();
    let mut prev = (t_lo, f(t_lo));
    for i in 1..=n {
        let t = t_lo + (t_hi - t_lo) * i as f64 / n as f64;
        let ft = f(t);
        if prev.1 * ft <= 0.0 && prev.1 != ft {
            let (mut a, mut b) = (prev.0, t);
            for _ in 0..100 {
                let mid = 0.5 * (a + b);
                if f(a) * f(mid) <= 0.0 {
                    b = mid;
                } else {
                    a = mid;
                }
            }
            let tr = 0.5 * (a + b);
            let h = 1e-3;
            let slope = r_th * (p_leak(tr + h) - p_leak(tr - h)) / (2.0 * h);
            out.push((tr, slope < 1.0));
        }
        prev = (t, ft);
    }
    out
}

#[cfg(test)]
mod tests {
    #[test]
    fn graetz_converges() {
        for (n, steps) in [(200usize, 2000usize), (400, 4000), (800, 8000)] {
            let nu = super::graetz(n, 0.5, steps);
            let at = |xs: f64| nu.iter().min_by(|a, b| (a.0.ln() - xs.ln()).abs().partial_cmp(&(b.0.ln() - xs.ln()).abs()).unwrap()).unwrap().1;
            eprintln!("n {n}: Nu(1e-5) {:.4} Nu(1e-4) {:.4} Nu(1e-3) {:.4} Nu(0.1) {:.5} Nu(0.5) {:.5}", at(1e-5), at(1e-4), at(1e-3), at(0.1), at(0.5));
            // fully developed: Nu = 140/17 exactly (parallel plates, H1)
            assert!((at(0.5) / (140.0 / 17.0) - 1.0).abs() < 2e-5);
        }
    }
}
