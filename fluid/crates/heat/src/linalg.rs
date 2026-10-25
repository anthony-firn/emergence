//! Small numerical helpers: Gauss–Legendre quadrature, dense LU solve, and
//! exponential integrals. No external dependencies.


/// Gauss–Legendre nodes and weights on `[0, 1]`.
///
/// Newton iteration on the Legendre polynomial `P_n`, starting from the
/// Chebyshev approximation of the roots (Abramowitz & Stegun 25.4.29).
pub fn gauss_legendre_01(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    for i in 0..n {
        // initial guess for the i-th root of P_n on [-1, 1]
        let mut z = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        for _ in 0..100 {
            // evaluate P_n(z) and P_n'(z) by the three-term recurrence
            let (mut p0, mut p1) = (1.0, 0.0);
            for j in 0..n {
                let p2 = p1;
                p1 = p0;
                p0 = ((2.0 * j as f64 + 1.0) * z * p1 - j as f64 * p2) / (j as f64 + 1.0);
            }
            let dp = n as f64 * (z * p0 - p1) / (z * z - 1.0);
            let dz = p0 / dp;
            z -= dz;
            if dz.abs() < 1e-15 {
                break;
            }
        }
        // weights w_i = 2 / ((1 - z^2) P_n'(z)^2)
        let (mut p0, mut p1) = (1.0, 0.0);
        for j in 0..n {
            let p2 = p1;
            p1 = p0;
            p0 = ((2.0 * j as f64 + 1.0) * z * p1 - j as f64 * p2) / (j as f64 + 1.0);
        }
        let dp = n as f64 * (z * p0 - p1) / (z * z - 1.0);
        // nodes are found in descending order; store ascending on [0, 1]
        x[n - 1 - i] = 0.5 * (z + 1.0);
        w[n - 1 - i] = 1.0 / ((1.0 - z * z) * dp * dp);
    }
    (x, w)
}

/// Solve `A x = b` by LU decomposition with partial pivoting. `a` is an
/// `n x n` matrix in row-major `Vec<Vec<f64>>`, `b` is overwritten by the
/// solution. Returns `false` if the matrix is singular.
pub fn solve_dense(a: &mut [Vec<f64>], b: &mut [f64]) -> bool {
    match Lu::new(a) {
        Some(lu) => {
            let x = lu.solve(b);
            b.copy_from_slice(&x);
            true
        }
        None => false,
    }
}

/// A reusable dense LU factorisation with partial pivoting.
pub struct Lu {
    lu: Vec<Vec<f64>>,
    piv: Vec<usize>,
}

impl Lu {
    /// Factorise `a` in place (cloned into the struct). `None` if singular.
    pub fn new(a: &[Vec<f64>]) -> Option<Lu> {
        let n = a.len();
        let mut lu = a.to_vec();
        let mut piv: Vec<usize> = (0..n).collect();
        for k in 0..n {
            let mut p = k;
            let mut best = lu[k][k].abs();
            for i in k + 1..n {
                if lu[i][k].abs() > best {
                    best = lu[i][k].abs();
                    p = i;
                }
            }
            if best < 1e-300 {
                return None;
            }
            if p != k {
                lu.swap(k, p);
                piv.swap(k, p);
            }
            let pivv = lu[k][k];
            for i in k + 1..n {
                let f = lu[i][k] / pivv;
                lu[i][k] = f;
                if f != 0.0 {
                    for j in k + 1..n {
                        lu[i][j] -= f * lu[k][j];
                    }
                }
            }
        }
        Some(Lu { lu, piv })
    }

    /// Solve `A x = b`.
    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = b.len();
        let mut y = vec![0.0; n];
        for i in 0..n {
            let mut s = b[self.piv[i]];
            for j in 0..i {
                s -= self.lu[i][j] * y[j];
            }
            y[i] = s;
        }
        let mut x = y;
        for i in (0..n).rev() {
            let mut s = x[i];
            for j in i + 1..n {
                s -= self.lu[i][j] * x[j];
            }
            x[i] = s / self.lu[i][i];
        }
        x
    }
}

/// Exponential integral `E_n(z) = ∫_0^1 μ^{n-2} e^{-z/μ} dμ`, `n >= 1`.
///
/// This is the standard `E_n` of radiative transfer. For `0 < z ≤ 1` the
/// power series for `E_1` is summed and the recurrence
/// `E_{n+1}(z) = (e^{-z} − z E_n(z))/n` used for the higher orders; for
/// `z > 1` the same recurrence is seeded from a modified-Lentz evaluation of
/// the continued fraction (Abramowitz & Stegun 5.1.22).
pub fn expn(n: usize, z: f64) -> f64 {
    if n == 0 {
        return (-z).exp() / z.max(1e-300);
    }
    if z <= 0.0 {
        return 1.0 / (n as f64 - 1.0).max(1.0);
    }
    let e1 = if z <= 1.0 {
        // E_1(z) = -γ - ln z + Σ_{k≥1} (-1)^{k+1} z^k /(k k!)
        const GAMMA: f64 = 0.577_215_664_901_532_9;
        let mut fact = 1.0; // z^k/k!
        let mut sum = 0.0;
        for k in 1..=200 {
            fact *= z / k as f64;
            let t = fact / k as f64;
            if k % 2 == 1 {
                sum += t;
            } else {
                sum -= t;
            }
            if t < 1e-18 {
                break;
            }
        }
        -GAMMA - z.ln() + sum
    } else if z > 700.0 {
        0.0
    } else {
        // modified Lentz for E_1(z) = e^{-z} / (z + 1/(1 + 1/(z + 2/(1 + ...))))
        let tiny = 1e-300;
        let mut f = z;
        let mut c = f;
        let mut d = 0.0;
        for k in 1..=300 {
            let a = ((k + 1) / 2) as f64;
            let b = if k % 2 == 0 { z } else { 1.0 };
            d = b + a * d;
            if d.abs() < tiny {
                d = tiny;
            }
            c = b + a / c;
            if c.abs() < tiny {
                c = tiny;
            }
            d = 1.0 / d;
            let delta = c * d;
            f *= delta;
            if (delta - 1.0).abs() < 1e-15 {
                break;
            }
        }
        (-z).exp() / f
    };
    // recurrence upward: E_{n+1} = (e^{-z} - z E_n)/n
    let mut en = e1;
    for k in 1..n {
        en = ((-z).exp() - z * en) / k as f64;
    }
    en
}
