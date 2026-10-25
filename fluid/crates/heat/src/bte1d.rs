//! The diffusive → ballistic fold for a one-dimensional grey phonon gas.
//!
//! Set-up: a slab `0 ≤ x ≤ L` with two diffuse walls at temperatures `T₁` and
//! `T₂`, a grey (frequency-independent) phonon gas with volumetric heat
//! capacity `C`, group speed `v` and relaxation time `τ`. In the steady,
//! linearised relaxation-time approximation the intensity `i(x, μ)` obeys
//!
//! ```text
//! μ ∂i/∂x = (i₀ − i)/Λ ,   i₀(x) = ½ ∫_{-1}^{1} i(x, μ) dμ ,
//! ```
//!
//! with `Λ = vτ` the mean free path and `μ = cos θ`. Energy conservation fixes
//! `i₀`; the walls emit isotropically (`i(0,μ>0) = C T₁/2`, etc.).
//!
//! The **indicator** is the phonon Knudsen number `Kn = Λ/L`:
//! * `Kn ≪ 1` → Fourier: `q = (1/3) C v Λ (T₁−T₂)/L`, i.e.
//!   `κ_eff/κ₀ → 1`;
//! * `Kn ≫ 1` → ballistic (Landauer/Casimir): `q = (1/4) C v (T₁−T₂)`, i.e.
//!   `κ_eff/κ₀ → (3/4)/Kn`.
//!
//! Both an exact integral-equation solution (the benchmark) and a
//! discrete-ordinates (DOM) solver are provided.

use crate::linalg::{expn, gauss_legendre_01, Lu};

/// Result of a grey-BTE solve.
#[derive(Clone, Debug)]
pub struct GraySolution {
    /// `Kn = Λ/L`.
    pub kn: f64,
    /// Dimensionless flux `q/(C v ΔT)`.
    pub q: f64,
    /// Effective conductivity ratio `κ_eff/κ₀ = 3 q /(2 Kn)`.
    pub kappa_ratio: f64,
    /// The local equilibrium profile `i₀(x)/Δi` (integral-equation solver).
    pub profile: Vec<f64>,
}

/// Exact asymptotic effective-conductivity ratio:
/// Fourier `1` for `Kn→0`, ballistic `(3/4)/Kn` for `Kn→∞`.
/// A smooth interpolation is `(1 + (4/3)Kn)⁻¹`, which matches both limits.
pub fn ratio_interpolation(kn: f64) -> f64 {
    1.0 / (1.0 + 4.0 / 3.0 * kn)
}

/// Ballistic (Casimir/Landauer) limit of the effective-conductivity ratio.
pub fn ratio_ballistic(kn: f64) -> f64 {
    0.75 / kn
}

/// Fourier limit of the effective-conductivity ratio.
pub fn ratio_fourier() -> f64 {
    1.0
}

// ---------------------------------------------------------------------------
// Discrete ordinates (DOM)
// ---------------------------------------------------------------------------

/// Directional quadrature on `μ ∈ [-1, 1]`.
fn angular_quadrature(n: usize) -> (Vec<f64>, Vec<f64>) {
    let (x, w) = gauss_legendre_01(n);
    let mu: Vec<f64> = x.iter().map(|&xi| 2.0 * xi - 1.0).collect();
    let ww: Vec<f64> = w.iter().map(|&wi| 2.0 * wi).collect();
    (mu, ww)
}

/// Block-tridiagonal matrix for the upwind DOM.
///
/// Unknown ordering per cell: the `n_ord` angular intensities. The streaming
/// term couples the same angle in neighbouring cells, so the lower and upper
/// blocks are diagonal; the closure `i₀ = ½ Σ w i` makes the diagonal block
/// dense. Solved by block-Thomas elimination.
pub fn solve_gray_bte_dom(kn: f64, n_ord: usize, n_cell: usize) -> GraySolution {
    let m = n_cell;
    let n = n_ord;
    let h = 1.0 / m as f64;
    let lambda = kn;
    let (mu, w) = angular_quadrature(n);
    // diagonal blocks A_j, lower L_j and upper U_j (diagonal), rhs b_j
    let mut a = vec![vec![vec![0.0; n]; n]; m];
    let mut l = vec![vec![0.0; n]; m];
    let mut u = vec![vec![0.0; n]; m];
    let mut b = vec![vec![0.0; n]; m];
    for j in 0..m {
        for i in 0..n {
            let absmu = mu[i].abs();
            a[j][i][i] = absmu / h + 1.0 / lambda;
            for k in 0..n {
                a[j][i][k] -= w[k] / (2.0 * lambda);
            }
            if mu[i] > 0.0 {
                l[j][i] = -mu[i] / h;
                if j == 0 {
                    b[j][i] = mu[i] / h; // left wall i₁ = 1
                }
            } else {
                u[j][i] = -absmu / h;
                if j == m - 1 {
                    b[j][i] = absmu / h * 0.0; // right wall i₂ = 0
                }
            }
        }
    }
    // forward block sweep
    let mut lu: Vec<Option<Lu>> = Vec::with_capacity(m);
    let mut bprime: Vec<Vec<f64>> = vec![vec![0.0; n]; m];
    let a0 = &a[0];
    lu.push(Lu::new(a0));
    bprime[0] = b[0].clone();
    for j in 1..m {
        let (x, y) = {
            let prev = lu[j - 1].as_ref().unwrap();
            // X = A_{j-1}^{-1} U_{j-1}, U diagonal: column k scaled by U[k][k]
            let mut x = vec![vec![0.0; n]; n];
            for k in 0..n {
                let mut e = vec![0.0; n];
                e[k] = u[j - 1][k];
                let col = prev.solve(&e);
                for i in 0..n {
                    x[i][k] = col[i];
                }
            }
            let y = prev.solve(&bprime[j - 1]);
            (x, y)
        };
        let mut aj = a[j].clone();
        for i in 0..n {
            for k in 0..n {
                aj[i][k] -= l[j][i] * x[i][k];
            }
        }
        lu.push(Lu::new(&aj));
        for i in 0..n {
            bprime[j][i] = b[j][i] - l[j][i] * y[i];
        }
    }
    // back substitution
    let mut ijet = vec![vec![0.0; n]; m];
    for j in (0..m).rev() {
        let lu_j = lu[j].as_ref().unwrap();
        let mut rhs = bprime[j].clone();
        if j + 1 < m {
            for i in 0..n {
                rhs[i] -= u[j][i] * ijet[j + 1][i];
            }
        }
        ijet[j] = lu_j.solve(&rhs);
    }
    // flux, averaged over cells
    let mut q = 0.0;
    for j in 0..m {
        for i in 0..n {
            q += w[i] * mu[i] * ijet[j][i];
        }
    }
    q /= m as f64;
    // physical flux q_phys = C v ΔT · q, with the (i₁=1, i₂=0) normalisation
    // used here T₁ = 2/C, T₂ = 0 so ΔT = 2/C and q_phys = C v (2/C) q /? -- we
    // keep the dimensionless q = q_phys/(C v ΔT). Then κ_eff/κ₀ = 3q/(2Kn).
    let kappa_ratio = 3.0 * q / (2.0 * kn);
    GraySolution {
        kn,
        q,
        kappa_ratio,
        profile: vec![],
    }
}

// ---------------------------------------------------------------------------
// Exact integral equation (benchmark)
// ---------------------------------------------------------------------------

/// Antiderivative pieces for the E₁ cell integral.
fn e1_left(x: f64, c: f64, kn: f64) -> (f64, f64) {
    let z = (c - x) / kn;
    (kn * expn(2, z), x * kn * expn(2, z) - kn * kn * expn(3, z))
}
fn e1_right(x: f64, c: f64, kn: f64) -> (f64, f64) {
    let z = (x - c) / kn;
    (
        -kn * expn(2, z),
        -x * kn * expn(2, z) - kn * kn * expn(3, z),
    )
}

/// `∫_a^b φ(s) E₁(|c−s|/Kn) ds` for the linear hat `φ` with values
/// `(hat_a, hat_b)` at the cell ends.
fn cell_e1(c: f64, a: f64, b: f64, hat_a: f64, hat_b: f64, kn: f64) -> f64 {
    let h = b - a;
    let alpha = hat_a - (hat_b - hat_a) * a / h;
    let beta = (hat_b - hat_a) / h;
    let mut total = 0.0;
    if a < c {
        let x0 = a;
        let x1 = b.min(c);
        let (f0a, f1a) = e1_left(x0, c, kn);
        let (f0b, f1b) = e1_left(x1, c, kn);
        total += alpha * (f0b - f0a) + beta * (f1b - f1a);
    }
    if b > c {
        let x0 = a.max(c);
        let x1 = b;
        let (f0a, f1a) = e1_right(x0, c, kn);
        let (f0b, f1b) = e1_right(x1, c, kn);
        total += alpha * (f0b - f0a) + beta * (f1b - f1a);
    }
    total
}

/// `∫_0^1 i(x) E₂(x/Kn) dx` for piecewise-linear `i` on a uniform grid.
fn flux_integral(i: &[f64], kn: f64) -> f64 {
    let m = i.len() - 1;
    let h = 1.0 / m as f64;
    let mut total = 0.0;
    for j in 0..m {
        let a = j as f64 * h;
        let b = a + h;
        let ia = i[j];
        let ib = i[j + 1];
        let alpha = ia - (ib - ia) * a / h;
        let beta = (ib - ia) / h;
        let f0 = |x: f64| -kn * expn(3, x / kn);
        let f1 =
            |x: f64| -x * kn * expn(3, x / kn) - kn * kn * expn(4, x / kn);
        total += alpha * (f0(b) - f0(a)) + beta * (f1(b) - f1(a));
    }
    total
}

/// Exact (numerical) grey-BTE solution from the integral equation
///
/// ```text
/// i₀(x) = ½[i₁ E₂(x/Λ) + i₂ E₂((L−x)/Λ)] + (1/2Λ) ∫_0^L i₀ E₁(|x−x'|/Λ) dx' ,
/// ```
///
/// solved by a Nyström method with piecewise-linear `i₀` and exact cell
/// integrals. This is the published/derived benchmark for the DOM; it is
/// accurate while the grid resolves the `O(Kn)` boundary layer
/// (`m ≳ 20/Kn`).
pub fn solve_gray_bte_integral(kn: f64, m: usize) -> GraySolution {
    let lambda = kn;
    let n = m + 1;
    let h = 1.0 / m as f64;
    let mut a = vec![vec![0.0; n]; n];
    for i in 0..n {
        a[i][i] = 1.0;
    }
    for i in 0..n {
        let c = i as f64 * h;
        for j in 0..m {
            let lo = j as f64 * h;
            let hi = lo + h;
            a[i][j] -= 0.5 * cell_e1(c, lo, hi, 1.0, 0.0, lambda) / lambda;
            a[i][j + 1] -= 0.5 * cell_e1(c, lo, hi, 0.0, 1.0, lambda) / lambda;
        }
    }
    let mut b = vec![0.0; n];
    for i in 0..n {
        let c = i as f64 * h;
        b[i] = 0.5 * expn(2, c / lambda); // i₁ = 1, i₂ = 0
    }
    let mut lu = Lu::new(&a).expect("singular integral operator");
    let _ = &mut lu;
    let i0 = lu.solve(&b);
    let q = 0.5 - flux_integral(&i0, lambda) / lambda;
    let kappa_ratio = 3.0 * q / (2.0 * kn);
    GraySolution {
        kn,
        q,
        kappa_ratio,
        profile: i0,
    }
}
