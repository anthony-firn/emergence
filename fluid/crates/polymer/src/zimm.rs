//! Zimm model: Rouse chain with a pre-averaged hydrodynamic interaction.
//!
//! The Oseen tensor couples the beads. Pre-averaging over the Gaussian
//! equilibrium distribution replaces it by the scalar
//! `Ω_ij = (6π η_s)⁻¹ ⟨1/R_ij⟩ = (6π η_s)⁻¹ √(6/(π b²|i-j|))`, so with
//! `H_ij = δ_ij + h* / √|i-j|` (`h*` the hydrodynamic-interaction
//! parameter, `h* = 0` recovers Rouse) the equations of motion become
//! `ζ H dR/dt = -k A R + f`. Normal modes solve the generalized
//! eigenproblem `A φ = λ H φ`, and
//!
//! ```text
//! τ_p = ζ / (k λ_p).
//! ```
//!
//! For `h* > 0` the longest time crosses over from the free-draining
//! `τ_1 ∝ N²` to the Zimm `τ_1 ∝ N^{3/2}` (the `N^{3ν}` law with `ν = ½`);
//! in a good solvent `ν ≈ 0.588` gives `N^{1.76}`.
//!
//! References: B. H. Zimm, *J. Chem. Phys.* **24**, 269 (1956);
//! Doi & Edwards (1986), ch. 4; Rubinstein & Colby (2003), ch. 8.

/// A Zimm chain with pre-averaged hydrodynamic interaction.
#[derive(Clone, Copy, Debug)]
pub struct Zimm {
    /// number of beads `N`
    pub n_beads: usize,
    /// Kuhn length `b`
    pub b: f64,
    /// bead friction coefficient `ζ`
    pub zeta: f64,
    /// thermal energy `k_B T`
    pub kbt: f64,
    /// pre-averaged hydrodynamic-interaction parameter `h*` (`0` = Rouse)
    pub h_star: f64,
}

impl Zimm {
    pub fn spring_k(&self) -> f64 {
        3.0 * self.kbt / (self.b * self.b)
    }

    /// Rouse matrix `A` (N×N): diagonal 2 (1 at the free ends), -1 for
    /// neighbours.
    fn rouse_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.n_beads;
        let mut a = vec![vec![0.0; n]; n];
        for i in 0..n {
            a[i][i] = if i == 0 || i == n - 1 { 1.0 } else { 2.0 };
            if i + 1 < n {
                a[i][i + 1] = -1.0;
                a[i + 1][i] = -1.0;
            }
        }
        a
    }

    /// Pre-averaged friction matrix `H`.
    fn friction_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.n_beads;
        let mut h = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                h[i][j] = if i == j {
                    1.0
                } else {
                    self.h_star / ((i as f64 - j as f64).abs()).sqrt()
                };
            }
        }
        h
    }

    /// Generalized eigenproblem of the *mobility* form: the relaxation
    /// rates are the eigenvalues of `H A` (mobility × connectivity),
    /// where `H = I + h* K` is the dimensionless friction matrix. This is
    /// the Zimm problem; its small-`p` eigenvalues scale as `(p/N)^{3/2}`,
    /// giving `τ₁ ∝ N^{3/2}`. Computed symmetrically as `Lᵀ A L` with
    /// `H = L Lᵀ` (Cholesky), which is similar to `H A`.
    pub fn eigenvalues(&self) -> Vec<f64> {
        let n = self.n_beads;
        let a = self.rouse_matrix();
        let h = self.friction_matrix();
        let l = cholesky(&h);
        // M = Lᵀ A L (symmetric, eigenvalues = eigenvalues of H A).
        // First T = A L.
        let mut t = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for k in 0..n {
                    s += a[i][k] * l[k][j];
                }
                t[i][j] = s;
            }
        }
        let mut m = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for k in 0..n {
                    s += l[k][i] * t[k][j];
                }
                m[i][j] = s;
            }
        }
        let (mut eigs, _) = jacobi_eigen(&mut m);
        eigs.sort_by(|x, y| x.partial_cmp(y).unwrap());
        eigs
    }

    /// Relaxation times `τ_p = ζ/(k λ_p)` for the `N-1` nonzero modes,
    /// slowest first.
    pub fn taus(&self) -> Vec<f64> {
        let k = self.spring_k();
        let mut t: Vec<f64> = self
            .eigenvalues()
            .into_iter()
            .filter(|&l| l > 1e-9)
            .map(|l| self.zeta / (k * l))
            .collect();
        t.sort_by(|x, y| y.partial_cmp(x).unwrap());
        t
    }

    pub fn tau1(&self) -> f64 {
        self.taus()[0]
    }
}

/// Cholesky factorization `H = L Lᵀ` (lower triangular `L`).
fn cholesky(h: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = h.len();
    let mut l = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = h[i][j];
            for k in 0..j {
                s -= l[i][k] * l[j][k];
            }
            if i == j {
                l[i][j] = s.max(1e-300).sqrt();
            } else {
                l[i][j] = s / l[j][j];
            }
        }
    }
    l
}

/// Cyclic Jacobi eigensolver for a symmetric matrix. Returns eigenvalues
/// and the eigenvector matrix (columns).
pub fn jacobi_eigen(a: &mut [Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut v = vec![vec![0.0; n]; n];
    for i in 0..n {
        v[i][i] = 1.0;
    }
    for _ in 0..100 {
        let mut off = 0.0;
        for i in 0..n {
            for j in i + 1..n {
                off += a[i][j] * a[i][j];
            }
        }
        if off < 1e-24 {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let (vkp, vkq) = (v[k][p], v[k][q]);
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    ((0..n).map(|i| a[i][i]).collect(), v)
}

/// The Rouse limit of [`Zimm`] must reproduce the Rouse times.
#[cfg(test)]
mod tests {
    use super::*;
    use crate::rouse::Rouse;

    #[test]
    fn zero_hydrodynamic_interaction_is_rouse() {
        let z = Zimm { n_beads: 8, b: 1.0, zeta: 1.0, kbt: 1.0, h_star: 0.0 };
        let r = Rouse::new(8, 1.0, 1.0, 1.0, 1.0);
        let tz = z.taus();
        assert_eq!(tz.len(), 7);
        for p in 0..7 {
            let tr = r.tau(p + 1);
            assert!((tz[p] - tr).abs() / tr < 1e-9, "p={} zimm {} rouse {}", p + 1, tz[p], tr);
        }
    }

    #[test]
    fn zimm_longest_time_scales_as_three_halves() {
        let mk = |n| Zimm { n_beads: n, b: 1.0, zeta: 1.0, kbt: 1.0, h_star: 0.25 };
        let t20 = mk(20).tau1();
        let t40 = mk(40).tau1();
        let exponent = (t40 / t20).ln() / 2.0_f64.ln();
        assert!(exponent > 1.35 && exponent < 1.75, "exponent {exponent}");

        // Free draining gives ≈ N².
        let r = |n| Zimm { n_beads: n, b: 1.0, zeta: 1.0, kbt: 1.0, h_star: 0.0 };
        let e = (r(40).tau1() / r(20).tau1()).ln() / 2.0_f64.ln();
        assert!(e > 1.8, "Rouse exponent {e}");
    }
}
