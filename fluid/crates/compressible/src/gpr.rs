//! Unified continuum mechanics: the Godunov–Peshkov–Romenski (GPR) model
//! (Peshkov & Romenski 2016; Dumbser, Peshkov, Romenski & Zanotti 2016).
//! One hyperbolic system describes elastic solids, viscoelastic media and
//! viscous fluids. The state carries the distortion field A (a local
//! "basis" of the material's elastic memory):
//!   ∂ρ/∂t + ∂_k(ρv_k) = 0
//!   ∂(ρv_i)/∂t + ∂_k(ρv_iv_k + pδ_ik − σ_ik) = 0
//!   ∂A_ik/∂t + ∂_k(A_im v_m) + v_j(∂_j A_ik − ∂_k A_ij) = −ψ_ik/θ(τ)
//!   ∂(ρE)/∂t + ∂_k(ρE v_k + (pδ_ik − σ_ik)v_i) = 0
//! with E = e(ρ, p) + (c_s²/4)|dev G|² + |v|²/2, G = AᵀA,
//! σ = −ρc_s² G dev G, ψ = c_s² A dev G, θ = τc_s²|A|^{−5/3}/3.
//!
//! Emergence. τ is the time the material takes to forget a shear strain.
//! Against the flow time t, the Deborah number De = τ/t picks the rung:
//! De ≫ 1 is an elastic solid (τ = ∞ gives finite-strain elasticity, shear
//! waves at c_s), De ≪ 1 is a Newtonian fluid with μ = ρ₀c_s²τ/6 (the
//! relaxation collapses onto the slow manifold dev G ≈ 0 and the
//! Navier–Stokes stress emerges), De ~ 1 is viscoelastic. Melting is τ
//! collapsing; plasticity is τ dropping where the stress passes the yield
//! stress.
//!
//! Relaxation is integrated per cell after each step (`relax`). It acts
//! only on the singular values of A: with A = UΣVᵀ the source keeps U and V
//! and det A (so ρ/ρ₀ = det A stays consistent), and the squared singular
//! values y_i (eigenvalues of G) obey dy_i/dt = −2k y_i (y_i − ȳ),
//! k = 3 det(A)^{5/3}/τ, integrated exponentially in ln y (exact for small
//! strain at any Δt, L-stable). Use Strang splitting (`Solver::strang`):
//! first-order splitting leaves an O(Δt/τ) error in the emergent viscosity.

use crate::eos::Nasg;
use crate::system::{Mat, System};

pub const NG: usize = 14;
type V = [f64; NG];
type M = Mat<NG>;
const IA: usize = 4;
const IP: usize = 13;

#[derive(Clone, Copy, Debug)]
pub struct Gpr {
    pub eos: Nasg,
    /// shear sound speed c_s
    pub cs: f64,
    /// strain relaxation time τ (f64::INFINITY: elastic solid); used by
    /// `TauModel::Constant`
    pub tau: f64,
    pub tau_model: TauModel,
}

/// Closures for the strain-relaxation time: the material's emergent
/// rheology.
#[derive(Clone, Copy, Debug)]
pub enum TauModel {
    /// τ = `Gpr::tau`
    Constant,
    /// Elasto-plasticity: τ = τ₀ (σ_Y/σ_eq)ⁿ with the von Mises stress
    /// σ_eq. Below yield τ is huge (elastic), above it the shear stress
    /// relaxes within τ₀; n → ∞ is ideal plasticity.
    Plastic { tau0: f64, sigma_y: f64, n: f64 },
    /// Melting: τ goes from `solid` to `liquid` (log-linear) across
    /// [t_melt − width, t_melt + width] in temperature.
    Melting { solid: f64, liquid: f64, t_melt: f64, width: f64 },
}

#[inline]
fn amat(w: &V) -> [[f64; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|k| w[IA + 3 * i + k]))
}

#[inline]
fn gram(a: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    std::array::from_fn(|k| std::array::from_fn(|l| (0..3).map(|i| a[i][k] * a[i][l]).sum()))
}

#[inline]
fn dev(g: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let tr = (g[0][0] + g[1][1] + g[2][2]) / 3.0;
    std::array::from_fn(|k| std::array::from_fn(|l| g[k][l] - if k == l { tr } else { 0.0 }))
}

#[inline]
fn det3(a: &[[f64; 3]; 3]) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

/// Eigen-decomposition of a symmetric 3×3 matrix (cyclic Jacobi):
/// returns eigenvalues and eigenvectors as columns of V.
pub fn sym_eig3(m: &[[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut a = *m;
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..50 {
        let off = a[0][1].abs() + a[0][2].abs() + a[1][2].abs();
        if off < 1e-15 * (a[0][0].abs() + a[1][1].abs() + a[2][2].abs() + 1e-300) {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if a[p][q].abs() < 1e-300 {
                continue;
            }
            let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let t = if theta == 0.0 { 1.0 } else { t };
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            for k in 0..3 {
                let (akp, akq) = (a[k][p], a[k][q]);
                a[k][p] = c * akp - s * akq;
                a[k][q] = s * akp + c * akq;
            }
            for k in 0..3 {
                let (apk, aqk) = (a[p][k], a[q][k]);
                a[p][k] = c * apk - s * aqk;
                a[q][k] = s * apk + c * aqk;
            }
            for k in 0..3 {
                let (vkp, vkq) = (v[k][p], v[k][q]);
                v[k][p] = c * vkp - s * vkq;
                v[k][q] = s * vkp + c * vkq;
            }
        }
    }
    ([a[0][0], a[1][1], a[2][2]], v)
}

impl Gpr {
    /// Elastic energy per unit mass E₂ = (c_s²/4)|dev G|².
    #[inline]
    fn e2(&self, a: &[[f64; 3]; 3]) -> f64 {
        let dg = dev(&gram(a));
        let s: f64 = dg.iter().flatten().map(|x| x * x).sum();
        0.25 * self.cs * self.cs * s
    }

    /// Stress σ = −ρc_s² G dev G.
    pub fn stress(&self, w: &V) -> [[f64; 3]; 3] {
        let a = amat(w);
        let g = gram(&a);
        let dg = dev(&g);
        let f = -w[0] * self.cs * self.cs;
        std::array::from_fn(|i| std::array::from_fn(|k| f * (0..3).map(|m| g[i][m] * dg[m][k]).sum::<f64>()))
    }

    /// Relax the distortion over dt (see module docs). Returns true.
    pub fn relax(&self, u: &mut V, dt: f64) -> bool {
        if matches!(self.tau_model, TauModel::Constant) && !self.tau.is_finite() {
            return false;
        }
        let a: [[f64; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|k| u[IA + 3 * i + k]));
        let det = det3(&a);
        if det <= 0.0 {
            return false;
        }
        // keep the internal energy fixed: elastic energy released by
        // relaxation becomes heat (total energy is conserved)
        let g = gram(&a);
        let (y, vv) = sym_eig3(&g);
        // Exponential integrator in x_i = ln y_i: dx_i/dt = −2k(y_i − ȳ) is
        // linearized as −2kȳ(x_i − x̄), whose exact solution
        //   x_i ← x̄ + (x_i − x̄) e^{−2kȳh}
        // is exact for small strain at any step size and L-stable (h → ∞
        // lands on the slow manifold dev G = 0); Σx_i = ln det G is kept.
        // Sub-steps bound the nonlinear (large-strain) error.
        let mut x: [f64; 3] = std::array::from_fn(|i| y[i].max(1e-300).ln());
        let rho = u[0];
        let temp = {
            let w = self.to_prim(u);
            self.eos.temperature(w[0], w[IP])
        };
        // relaxation rate for the current strain (k = 3 det^{5/3} / τ)
        let rate = |x: &[f64; 3]| -> f64 {
            let tau = match self.tau_model {
                TauModel::Constant => self.tau,
                TauModel::Plastic { tau0, sigma_y, n } => {
                    // principal stresses in G's eigenbasis: −ρc_s² y_i (y_i − ȳ)
                    let yv: [f64; 3] = std::array::from_fn(|i| x[i].exp());
                    let yb = (yv[0] + yv[1] + yv[2]) / 3.0;
                    let sg: [f64; 3] = std::array::from_fn(|i| -rho * self.cs * self.cs * yv[i] * (yv[i] - yb));
                    let sm = (sg[0] + sg[1] + sg[2]) / 3.0;
                    let seq = (1.5 * sg.iter().map(|v| (v - sm).powi(2)).sum::<f64>()).sqrt();
                    if seq <= 0.0 { f64::INFINITY } else { tau0 * (sigma_y / seq).powf(n) }
                }
                TauModel::Melting { solid, liquid, t_melt, width } => {
                    let f = ((temp - (t_melt - width)) / (2.0 * width)).clamp(0.0, 1.0);
                    (solid.ln() * (1.0 - f) + liquid.ln() * f).exp()
                }
            };
            if tau.is_finite() { 3.0 * det.powf(5.0 / 3.0) / tau } else { 0.0 }
        };
        let k = rate(&x);
        let xm = (x[0] + x[1] + x[2]) / 3.0;
        let spread = x.iter().map(|v| (v - xm).abs()).fold(0.0, f64::max);
        let nsub = ((4.0 * dt * k * spread).ceil() as usize).clamp(1, 256);
        let h = dt / nsub as f64;
        // von Mises stress of a log-strain state (principal stresses in G's
        // eigenbasis are −ρc_s² y_i (y_i − ȳ))
        let seq_of = |x: &[f64; 3]| -> f64 {
            let yv: [f64; 3] = std::array::from_fn(|i| x[i].exp());
            let yb = (yv[0] + yv[1] + yv[2]) / 3.0;
            let sg: [f64; 3] = std::array::from_fn(|i| -rho * self.cs * self.cs * yv[i] * (yv[i] - yb));
            let sm = (sg[0] + sg[1] + sg[2]) / 3.0;
            (1.5 * sg.iter().map(|v| (v - sm).powi(2)).sum::<f64>()).sqrt()
        };
        for _ in 0..nsub {
            let k = rate(&x);
            let yb = x.iter().map(|v| v.exp()).sum::<f64>() / 3.0;
            let xm = (x[0] + x[1] + x[2]) / 3.0;
            let f = (-2.0 * k * yb * h).exp();
            let x0 = x;
            for v in x.iter_mut() {
                *v = xm + (*v - xm) * f;
            }
            // Radial return: relaxation cannot carry the stress below σ*,
            // the level whose relaxation time equals the sub-step (below it
            // the material is effectively elastic over the step). As n → ∞,
            // σ* → σ_Y: projection onto the yield surface.
            if let TauModel::Plastic { tau0, sigma_y, n } = self.tau_model {
                let s_star = sigma_y * (tau0 / h).powf(1.0 / n);
                if seq_of(&x0) > s_star && seq_of(&x) < s_star {
                    let (mut lo, mut hi) = (0.0f64, 1.0f64); // scale of the deviation
                    for _ in 0..60 {
                        let mid = 0.5 * (lo + hi);
                        let xt: [f64; 3] = std::array::from_fn(|i| xm + (x0[i] - xm) * mid);
                        if seq_of(&xt) < s_star { lo = mid } else { hi = mid }
                    }
                    let sc = 0.5 * (lo + hi);
                    x = std::array::from_fn(|i| xm + (x0[i] - xm) * sc);
                }
            }
        }
        let z: [f64; 3] = std::array::from_fn(|i| x[i].exp());
        // restore det G = det(A)² exactly
        let scale = (det * det / (z[0] * z[1] * z[2])).cbrt();
        let fac: [f64; 3] = std::array::from_fn(|i| (z[i] * scale / y[i]).sqrt());
        // A_new = A V diag(fac) Vᵀ
        let av: [[f64; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|c| (0..3).map(|m| a[i][m] * vv[m][c]).sum()));
        for i in 0..3 {
            for kk in 0..3 {
                u[IA + 3 * i + kk] = (0..3).map(|c| av[i][c] * fac[c] * vv[kk][c]).sum();
            }
        }
        true
    }

    fn c0(&self, w: &V) -> f64 {
        self.eos.sound_speed_sq(w[0], w[IP]).max(0.0).sqrt()
    }
}

impl System<NG> for Gpr {
    fn to_prim(&self, u: &V) -> V {
        let rho = u[0];
        let mut w = *u;
        let v = [u[1] / rho, u[2] / rho, u[3] / rho];
        w[1] = v[0];
        w[2] = v[1];
        w[3] = v[2];
        let a = amat(u);
        let e1 = u[IP] / rho - self.e2(&a) - 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        w[IP] = self.eos.pressure(rho, e1);
        w
    }

    #[inline]
    fn velocity_conserved(&self, u: &V, d: usize) -> f64 {
        u[1 + d] / u[0]
    }
    fn to_cons(&self, w: &V) -> V {
        let rho = w[0];
        let mut u = *w;
        u[1] = rho * w[1];
        u[2] = rho * w[2];
        u[3] = rho * w[3];
        let a = amat(w);
        let ke = 0.5 * (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
        u[IP] = rho * (self.eos.energy(rho, w[IP]) + self.e2(&a) + ke);
        u
    }

    /// Upper bound of the characteristic speeds (longitudinal waves,
    /// stiffened by finite strain).
    fn sound_speed(&self, w: &V) -> f64 {
        let g = gram(&amat(w));
        let s = ((g[0][0] + g[1][1] + g[2][2]) / 3.0).max(1.0);
        (self.c0(w).powi(2) + 4.0 / 3.0 * self.cs * self.cs * s * s * 1.2).sqrt()
    }

    fn vel_index(&self) -> usize {
        1
    }

    fn momentum_index(&self, d: usize) -> usize {
        1 + d
    }

    fn energy_index(&self) -> Option<usize> {
        Some(IP)
    }

    /// Component-wise reconstruction (no characteristic decomposition yet).
    fn eigen(&self, _w: &V, _d: usize, l: &mut M, r: &mut M) {
        for i in 0..NG {
            l[i] = [0.0; NG];
            r[i] = [0.0; NG];
            l[i][i] = 1.0;
            r[i][i] = 1.0;
        }
    }

    fn eigenvalues(&self, w: &V, d: usize) -> V {
        [w[1 + d].abs() + self.sound_speed(w); NG]
    }

    /// Rusanov flux with the speed bound; face velocity = mean.
    fn riemann_flux(&self, wl: &V, wr: &V, d: usize) -> (V, f64) {
        let s = self.max_speed(wl, d).max(self.max_speed(wr, d));
        let (fl, fr) = (self.flux(wl, d), self.flux(wr, d));
        let (ql, qr) = (self.to_cons(wl), self.to_cons(wr));
        let mut f = [0.0; NG];
        for k in 0..NG {
            f[k] = 0.5 * (fl[k] + fr[k]) - 0.5 * s * (qr[k] - ql[k]);
        }
        (f, 0.5 * (wl[1 + d] + wr[1 + d]))
    }

    /// Symmetric two-point flux: arithmetic means of the transported
    /// quantities, products split symmetrically.
    fn split_flux(&self, a: &V, b: &V, d: usize) -> V {
        let un = 1 + d;
        let rho = 0.5 * (a[0] + b[0]);
        let u = 0.5 * (a[un] + b[un]);
        let c = rho * u;
        let (sa, sb) = (self.stress(a), self.stress(b));
        let mut f = [0.0; NG];
        f[0] = c;
        for i in 0..3 {
            f[1 + i] = c * 0.5 * (a[1 + i] + b[1 + i]) - 0.5 * (sa[i][d] + sb[i][d]);
        }
        f[un] += 0.5 * (a[IP] + b[IP]);
        // A column d: A_im v_m
        for i in 0..3 {
            let mut s = 0.0;
            for m in 0..3 {
                s += 0.25 * (a[IA + 3 * i + m] + b[IA + 3 * i + m]) * (a[1 + m] + b[1 + m]);
            }
            f[IA + 3 * i + d] = s;
        }
        let kin = c * 0.5 * (a[1] * b[1] + a[2] * b[2] + a[3] * b[3]);
        let rhoe = 0.5 * (self.eos.rho_e(a[0], a[IP]) + self.eos.rho_e(b[0], b[IP]));
        let e2 = 0.5 * (a[0] * self.e2(&amat(a)) + b[0] * self.e2(&amat(b)));
        let pw = 0.5 * (a[un] * b[IP] + b[un] * a[IP]);
        let sw: f64 = (0..3).map(|i| 0.5 * (sa[i][d] * b[1 + i] + sb[i][d] * a[1 + i])).sum();
        f[IP] = kin + u * (rhoe + e2) + pw - sw;
        f
    }

    fn flux(&self, w: &V, d: usize) -> V {
        let un = 1 + d;
        let q = self.to_cons(w);
        let s = self.stress(w);
        let v = [w[1], w[2], w[3]];
        let mut f = [0.0; NG];
        f[0] = q[0] * v[d];
        for i in 0..3 {
            f[1 + i] = q[1 + i] * v[d] - s[i][d];
        }
        f[un] += w[IP];
        for i in 0..3 {
            f[IA + 3 * i + d] = (0..3).map(|m| w[IA + 3 * i + m] * v[m]).sum();
        }
        f[IP] = (q[IP] + w[IP]) * v[d] - (0..3).map(|i| s[i][d] * v[i]).sum::<f64>();
        f
    }

    /// Numerical ∂u/∂w (only used by the FD-WENO path).
    fn dudw(&self, w: &V) -> M {
        let mut j = [[0.0; NG]; NG];
        let u0 = self.to_cons(w);
        for c in 0..NG {
            let h = 1e-7 * w[c].abs().max(1e-7);
            let mut wp = *w;
            wp[c] += h;
            let up = self.to_cons(&wp);
            for r in 0..NG {
                j[r][c] = (up[r] - u0[r]) / h;
            }
        }
        j
    }

    fn admissible(&self, w: &V) -> bool {
        w.iter().all(|x| x.is_finite()) && w[0] > 0.0 && w[IP] + self.eos.p_inf > 0.0 && det3(&amat(w)) > 0.0
    }

    fn nc_general(&self) -> bool {
        true
    }

    /// −[v_d ∂_d A_ik − δ_kd Σ_j v_j ∂_d A_ij]
    fn nc_rhs(&self, w: &V, dw: &V, d: usize) -> V {
        let mut r = [0.0; NG];
        for i in 0..3 {
            let s: f64 = (0..3).map(|j| w[1 + j] * dw[IA + 3 * i + j]).sum();
            for k in 0..3 {
                let mut t = w[1 + d] * dw[IA + 3 * i + k];
                if k == d {
                    t -= s;
                }
                r[IA + 3 * i + k] = -t;
            }
        }
        r
    }

    fn sensor(&self, w: &V, out: &mut [f64; 8]) -> usize {
        out[0] = w[0].max(1e-300).ln();
        out[1] = (w[IP] + self.eos.p_inf).max(1e-300).ln();
        let pairs = [(0, 1), (0, 2), (1, 2), (1, 0), (2, 0), (2, 1)];
        for (n, (i, k)) in pairs.iter().enumerate() {
            out[2 + n] = w[IA + 3 * i + k];
        }
        8
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn relaxation_keeps_det_and_kills_shear_strain() {
        let g = Gpr { eos: Nasg::ideal(1.4), cs: 1.0, tau: 1e-3, tau_model: TauModel::Constant };
        let mut u = [0.0; NG];
        u[0] = 1.0;
        let a = [[1.0, 0.2, 0.0], [0.05, 0.9, 0.1], [0.0, 0.0, 1.1]];
        for i in 0..3 {
            for k in 0..3 {
                u[IA + 3 * i + k] = a[i][k];
            }
        }
        let d0 = det3(&a);
        g.relax(&mut u, 1.0);
        let an: [[f64; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|k| u[IA + 3 * i + k]));
        let dg = dev(&gram(&an));
        let s: f64 = dg.iter().flatten().map(|x| x.abs()).sum();
        assert!((det3(&an) - d0).abs() < 1e-12, "det changed");
        assert!(s < 1e-8, "dev G not relaxed: {s}");
    }

    #[test]
    fn eigen3_reconstructs() {
        let m = [[2.0, 0.3, -0.1], [0.3, 1.0, 0.2], [-0.1, 0.2, 0.5]];
        let (l, v) = sym_eig3(&m);
        for i in 0..3 {
            for j in 0..3 {
                let r: f64 = (0..3).map(|c| v[i][c] * l[c] * v[j][c]).sum();
                assert!((r - m[i][j]).abs() < 1e-12);
            }
        }
    }
}
