//! Quantum confinement in MOS inversion layers (1D, depth z from the
//! Si/SiO₂ interface): the rung below drift–diffusion.
//!
//! - Schrödinger–Poisson (the unfold): effective-mass subbands of the two
//!   (100) valley families — Δ2 (m_z 0.916, m_d 0.19, g 2) and Δ4 (m_z
//!   0.19, m_d 0.417, g 4) — in the self-consistent potential, occupied
//!   with 2D Fermi–Dirac statistics, n(z) = Σ g m_d kT/(πħ²)
//!   ln(1 + e^{(E_F−E_i)/kT}) |ψ_i(z)|²; wavefunctions vanish at the oxide
//!   (barrier 3.15 eV) and at the end of the quantum box. Poisson with the
//!   oxide as a Robin condition; Trellakis predictor–corrector iteration.
//! - Density gradient (the fold): n = n_i e^{ψ + Λ}, Λ = 2b (√n)''/√n with
//!   b = ħ²/(12 q m_dg V_t), √n = 0 at the oxide; Newton on (ψ, √n). One
//!   parameter, m_dg, calibrated on Schrödinger–Poisson and certified on
//!   other biases, dopings and oxide thicknesses.
//! - Classical (Boltzmann) electrons: the drift–diffusion limit.
//!
//! Equilibrium (φ = 0), p-type substrate, n⁺ poly gate (offset +E_g/2).

use crate::band::Band;
use crate::mesh::graded;
use crate::silicon::{band_gap, nc, Silicon, EPS_OX, EPS_SI, KB, Q};
use crate::tunnel::{HBAR, M0};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Electrons {
    Classical,
    SchrodingerPoisson,
    DensityGradient,
}

#[derive(Clone, Debug)]
pub struct MosStack {
    pub na: f64,
    pub tox: f64,
    pub t: f64,
    /// depth of the silicon domain (m)
    pub depth: f64,
    /// quantum box (m): electrons inside are quantized
    pub lq: f64,
    /// density-gradient mass (units of m₀)
    pub m_dg: f64,
}

#[derive(Clone, Debug)]
pub struct MosSolution {
    pub z: Vec<f64>,
    /// potential (units of V_t, intrinsic-referenced)
    pub psi: Vec<f64>,
    /// electron density (m⁻³)
    pub n: Vec<f64>,
    /// inversion charge per area (C/m², positive number)
    pub q_inv: f64,
    /// electron centroid (m)
    pub centroid: f64,
    /// lowest subband energy above the surface band edge (eV), SP only
    pub e0: f64,
    pub iterations: usize,
}

impl MosStack {
    fn mesh(&self) -> Vec<f64> {
        graded(&[(0.0, 0.02e-9), (2e-9, 0.1e-9), (self.lq, 0.3e-9), (self.depth, 5e-9)])
    }

    /// Solve at gate voltage vg.
    pub fn solve(&self, vg: f64, mode: Electrons) -> MosSolution {
        let si = Silicon::new(self.t);
        let vt = si.vt();
        let ni = si.ni;
        let z = self.mesh();
        let nz = z.len();
        let psi_b = -(0.5 * self.na / ni).asinh();
        let psi_g = (vg + 0.5 * band_gap(self.t)) / vt;
        let cox = EPS_OX / self.tox;
        // dual widths
        let w: Vec<f64> = (0..nz)
            .map(|i| {
                let l = if i > 0 { z[i] - z[i - 1] } else { 0.0 };
                let r = if i + 1 < nz { z[i + 1] - z[i] } else { 0.0 };
                0.5 * (l + r)
            })
            .collect();
        let chi_c = (nc(self.t) / ni).ln();
        let kt = KB * self.t;
        // initial potential: classical surface potential guess, linear decay
        let mut psi: Vec<f64> = z.iter().map(|&zz| psi_b + (psi_g - psi_b) * 0.3 * (-zz / 20e-9).exp()).collect();
        let iq = z.iter().position(|&zz| zz >= self.lq).unwrap_or(nz - 1);
        // electron density (in units of n_i) for the present potential
        let mut n_hat = vec![0.0; nz];
        let mut e0 = f64::NAN;
        let mut u: Vec<f64> = psi.iter().map(|&p| (0.5 * p).exp()).collect();
        u[0] = 0.0;
        let mut its = 0;
        for outer in 0..200 {
            its = outer + 1;
            match mode {
                Electrons::Classical => {
                    for i in 0..nz {
                        n_hat[i] = psi[i].exp();
                    }
                }
                Electrons::SchrodingerPoisson => {
                    let (nq, e_low) = self.subband_density(&z, &w, &psi, iq, chi_c, kt, ni);
                    e0 = e_low;
                    for i in 0..nz {
                        n_hat[i] = if i < iq { nq[i] } else { psi[i].exp() };
                    }
                }
                Electrons::DensityGradient => {}
            }
            // Poisson (and DG) Newton with the predictor n(ψ) = n_old e^{ψ − ψ_old}
            let psi_old = psi.clone();
            let n_old = n_hat.clone();
            let b = HBAR * HBAR / (12.0 * Q * self.m_dg * M0 * vt);
            let nv = if mode == Electrons::DensityGradient { 2 } else { 1 };
            let mut converged_inner = false;
            for _ in 0..100 {
                let mut a = Band::new(nv * nz, 2 * nv, 2 * nv);
                let mut f = vec![0.0; nv * nz];
                for i in 0..nz {
                    let r = nv * i;
                    // Poisson: Σ ε (ψ_j − ψ_i)/h + q n_i w (p − n − N_a/n_i)/V_t·(1/ε…) = 0
                    let mut lap = 0.0;
                    for (j, sgn) in [(i.wrapping_sub(1), -1i32), (i + 1, 1)] {
                        if j < nz {
                            let h = (z[j] - z[i]).abs();
                            let g = EPS_SI / h;
                            lap += g * (psi[j] - psi[i]);
                            a.add(r, nv * j, g);
                            a.add(r, r, -g);
                            let _ = sgn;
                        }
                    }
                    if i == 0 {
                        // oxide: ε_ox (ψ_g − ψ_s)/t_ox flux into the first box
                        lap += cox * (psi_g - psi[0]);
                        a.add(r, r, -cox);
                    }
                    let p_hat = (-psi[i]).exp();
                    let (nn, dn_dpsi, dn_du) = match mode {
                        Electrons::DensityGradient => (u[i] * u[i], 0.0, 2.0 * u[i]),
                        _ => {
                            let v = n_old[i] * (psi[i] - psi_old[i]).exp();
                            (v, v, 0.0)
                        }
                    };
                    let kappa = Q * ni / vt;
                    let rho = p_hat - nn - self.na / ni;
                    f[r] = lap + kappa * w[i] * rho;
                    a.add(r, r, kappa * w[i] * (-p_hat - dn_dpsi));
                    if mode == Electrons::DensityGradient {
                        a.add(r, r + 1, -kappa * w[i] * dn_du);
                    }
                    if i == nz - 1 {
                        // bulk: neutral potential
                        for c in [r.saturating_sub(nv), r] {
                            let v = a.get(r, c);
                            if v != 0.0 {
                                a.add(r, c, -v);
                            }
                        }
                        a.add(r, r, 1.0);
                        f[r] = psi[i] - psi_b;
                    }
                    if mode == Electrons::DensityGradient {
                        // 2b u'' − u (2 ln u − ψ) = 0 (box form), u(0) = 0,
                        // u(end) = bulk value
                        let r2 = r + 1;
                        if i == 0 || i == nz - 1 {
                            let target = if i == 0 { 0.0 } else { (0.5 * psi_b).exp() };
                            a.add(r2, r2, 1.0);
                            f[r2] = u[i] - target;
                        } else {
                            let mut lapu = 0.0;
                            for j in [i - 1, i + 1] {
                                let h = (z[j] - z[i]).abs();
                                lapu += 2.0 * b * (u[j] - u[i]) / h;
                                a.add(r2, nv * j + 1, 2.0 * b / h);
                                a.add(r2, r2, -2.0 * b / h);
                            }
                            let ui = u[i].max(1e-300);
                            let src = ui * (2.0 * ui.ln() - psi[i]);
                            f[r2] = lapu - w[i] * src;
                            a.add(r2, r2, -w[i] * (2.0 * ui.ln() - psi[i] + 2.0));
                            a.add(r2, r, w[i] * ui);
                        }
                    }
                }
                f.iter_mut().for_each(|v| *v = -*v);
                a.solve(&mut f).expect("Poisson solve");
                let dmax = (0..nz).map(|i| f[nv * i].abs()).fold(0.0, f64::max);
                let s = (2.0 / dmax).min(1.0);
                for i in 0..nz {
                    psi[i] += s * f[nv * i];
                    if mode == Electrons::DensityGradient {
                        // keep √n positive: limit relative decrease
                        let du = s * f[nv * i + 1];
                        u[i] = if i == 0 { 0.0 } else { (u[i] + du).max(0.1 * u[i]) };
                    }
                }
                let du_max = if mode == Electrons::DensityGradient { (0..nz).map(|i| (f[nv * i + 1] / u[i].max(1e-300)).abs()).fold(0.0, f64::max) } else { 0.0 };
                if dmax * s < 1e-11 && du_max < 1e-9 {
                    converged_inner = true;
                    break;
                }
            }
            let _ = converged_inner;
            let change = psi.iter().zip(&psi_old).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
            if mode != Electrons::SchrodingerPoisson || change < 1e-8 {
                if mode == Electrons::DensityGradient {
                    for i in 0..nz {
                        n_hat[i] = u[i] * u[i];
                    }
                } else if mode == Electrons::Classical {
                    for i in 0..nz {
                        n_hat[i] = psi[i].exp();
                    }
                }
                break;
            }
        }
        let n: Vec<f64> = n_hat.iter().map(|v| v * ni).collect();
        // inversion charge (electrons above the bulk value) and centroid
        let nb = ni * (psi_b).exp();
        let (mut qi, mut zc) = (0.0, 0.0);
        for i in 0..nz {
            let ex = (n[i] - nb).max(0.0) * w[i];
            qi += ex;
            zc += ex * z[i];
        }
        MosSolution { centroid: zc / qi.max(1e-300), q_inv: Q * qi, z, psi, n, e0, iterations: its }
    }

    /// Electron density (units of n_i) from the subbands of both valley
    /// families on [0, z_iq), and the lowest subband energy (eV) above the
    /// surface band edge.
    #[allow(clippy::too_many_arguments)]
    fn subband_density(&self, z: &[f64], w: &[f64], psi: &[f64], iq: usize, chi_c: f64, kt: f64, ni: f64) -> (Vec<f64>, f64) {
        let nz = z.len();
        let mut n = vec![0.0; nz];
        let ec_s = kt * (chi_c - psi[0]);
        let mut e_low = f64::INFINITY;
        // interior nodes 1..iq−1 (ψ_wf = 0 at 0 and iq)
        let m = iq - 1;
        for (mz, md, g) in [(0.916, 0.19, 2.0), (0.19, 0.417, 4.0)] {
            // −ħ²/(2 m_z) d²/dz² + V, V = kT(χ_c − ψ): symmetric form with the
            // dual widths, H̃ = W^{−1/2} K W^{−1/2} + V
            let c = HBAR * HBAR / (2.0 * mz * M0);
            let mut diag = vec![0.0; m];
            let mut off = vec![0.0; m.saturating_sub(1)];
            for k in 0..m {
                let i = k + 1;
                let (hl, hr) = (z[i] - z[i - 1], z[i + 1] - z[i]);
                diag[k] = c * (1.0 / hl + 1.0 / hr) / w[i] + kt * (chi_c - psi[i]);
                if k + 1 < m {
                    off[k] = -c / hr / (w[i] * w[i + 1]).sqrt();
                }
            }
            let nev = 12;
            let evals = tridiag_eigenvalues(&diag, &off, nev);
            let dos = g * md * M0 * kt / (std::f64::consts::PI * HBAR * HBAR);
            for &e in &evals {
                let occ_arg = -e / kt; // E_F = 0
                if occ_arg < -40.0 {
                    continue;
                }
                e_low = e_low.min(e);
                let v = tridiag_eigenvector(&diag, &off, e);
                // normalize ∫|φ|² dz = 1 with φ = W^{−1/2} ζ → Σ ζ² = 1
                let norm: f64 = v.iter().map(|x| x * x).sum();
                let occ = if occ_arg > 40.0 { occ_arg } else { occ_arg.exp().ln_1p() };
                for k in 0..m {
                    let i = k + 1;
                    n[i] += dos * occ * v[k] * v[k] / norm / w[i] / ni;
                }
            }
        }
        (n, (e_low - ec_s) / Q)
    }
}

/// Lowest `count` eigenvalues of a symmetric tridiagonal matrix (Sturm
/// bisection).
pub fn tridiag_eigenvalues(d: &[f64], e: &[f64], count: usize) -> Vec<f64> {
    let n = d.len();
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for i in 0..n {
        let r = if i > 0 { e[i - 1].abs() } else { 0.0 } + if i + 1 < n { e[i].abs() } else { 0.0 };
        lo = lo.min(d[i] - r);
        hi = hi.max(d[i] + r);
    }
    // number of eigenvalues below x
    let below = |x: f64| -> usize {
        let mut cnt = 0;
        let mut q = d[0] - x;
        if q < 0.0 {
            cnt += 1;
        }
        for i in 1..n {
            let qq = if q.abs() < 1e-300 { 1e-300 } else { q };
            q = d[i] - x - e[i - 1] * e[i - 1] / qq;
            if q < 0.0 {
                cnt += 1;
            }
        }
        cnt
    };
    (0..count.min(n))
        .map(|k| {
            let (mut a, mut b) = (lo, hi);
            for _ in 0..200 {
                let mid = 0.5 * (a + b);
                if below(mid) > k {
                    b = mid;
                } else {
                    a = mid;
                }
            }
            0.5 * (a + b)
        })
        .collect()
}

/// Eigenvector for eigenvalue λ by inverse iteration (tridiagonal solves).
pub fn tridiag_eigenvector(d: &[f64], e: &[f64], lambda: f64) -> Vec<f64> {
    let n = d.len();
    let shift = lambda + 1e-12 * lambda.abs().max(1e-30);
    let mut v = vec![1.0; n];
    for _ in 0..3 {
        // solve (T − shift) x = v (Thomas with partial safety)
        let mut cp = vec![0.0; n];
        let mut dp = vec![0.0; n];
        let mut b0 = d[0] - shift;
        if b0.abs() < 1e-300 {
            b0 = 1e-300;
        }
        cp[0] = if n > 1 { e[0] / b0 } else { 0.0 };
        dp[0] = v[0] / b0;
        for i in 1..n {
            let mut den = d[i] - shift - e[i - 1] * cp[i - 1];
            if den.abs() < 1e-300 {
                den = 1e-300;
            }
            cp[i] = if i + 1 < n { e[i] / den } else { 0.0 };
            dp[i] = (v[i] - e[i - 1] * dp[i - 1]) / den;
        }
        let mut x = vec![0.0; n];
        x[n - 1] = dp[n - 1];
        for i in (0..n - 1).rev() {
            x[i] = dp[i] - cp[i] * x[i + 1];
        }
        let norm = x.iter().map(|a| a * a).sum::<f64>().sqrt();
        v = x.iter().map(|a| a / norm).collect();
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    /// Infinite square well: E_k = (kπħ)²/(2 m L²) on the same discrete
    /// operator (graded mesh), second order.
    #[test]
    fn square_well() {
        let l = 5e-9;
        let m = 0.19 * M0;
        let z = graded(&[(0.0, 0.02e-9), (l, 0.1e-9)]);
        let nz = z.len();
        let w: Vec<f64> = (0..nz).map(|i| 0.5 * (if i > 0 { z[i] - z[i - 1] } else { 0.0 } + if i + 1 < nz { z[i + 1] - z[i] } else { 0.0 })).collect();
        let c = HBAR * HBAR / (2.0 * m);
        let mm = nz - 2;
        let mut d = vec![0.0; mm];
        let mut e = vec![0.0; mm - 1];
        for k in 0..mm {
            let i = k + 1;
            let (hl, hr) = (z[i] - z[i - 1], z[i + 1] - z[i]);
            d[k] = c * (1.0 / hl + 1.0 / hr) / w[i];
            if k + 1 < mm {
                e[k] = -c / hr / (w[i] * w[i + 1]).sqrt();
            }
        }
        let ev = tridiag_eigenvalues(&d, &e, 3);
        for (k, val) in ev.iter().enumerate() {
            let exact = ((k + 1) as f64 * std::f64::consts::PI * HBAR).powi(2) / (2.0 * m * l * l);
            assert!((val / exact - 1.0).abs() < 2e-3, "level {k}: {val} vs {exact}");
        }
    }
}
