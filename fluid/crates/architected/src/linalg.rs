//! Small dense linear algebra for unit-cell and frame solves.
//!
//! Unit cells and benchmark frames are small (tens to a few thousand degrees
//! of freedom), so a plain row-major dense LU with partial pivoting is used.
//! The only reason this exists rather than a dependency is the house rule that
//! the physics crates stay self-contained; the routines are textbook and are
//! checked by the tests in this module.

/// Row-major dense square matrix.
#[derive(Clone, Debug)]
pub struct Dense {
    pub n: usize,
    pub a: Vec<f64>,
}

impl Dense {
    pub fn zeros(n: usize) -> Self {
        Self {
            n,
            a: vec![0.0; n * n],
        }
    }
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.a[i * self.n + j]
    }
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.a[i * self.n + j] = v;
    }
    #[inline]
    pub fn add(&mut self, i: usize, j: usize, v: f64) {
        self.a[i * self.n + j] += v;
    }
    pub fn add_outer_block(&mut self, rows: &[usize], k: &[f64]) {
        // k is a dense m*m matrix in the same order as rows.
        let m = rows.len();
        for (ii, &ri) in rows.iter().enumerate() {
            for (jj, &rj) in rows.iter().enumerate() {
                self.a[ri * self.n + rj] += k[ii * m + jj];
            }
        }
    }
    pub fn mul_vec(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.n];
        for i in 0..self.n {
            let mut s = 0.0;
            for j in 0..self.n {
                s += self.a[i * self.n + j] * x[j];
            }
            y[i] = s;
        }
        y
    }
}

/// Solve `A x = b` in place. `b` is overwritten with the solution.
/// Returns `None` if `A` is singular to working precision.
pub fn solve(a: &mut Dense, b: &mut [f64]) -> Option<()> {
    let n = a.n;
    assert_eq!(b.len(), n);
    // Partial-pivoting LU.
    for k in 0..n {
        let mut piv = k;
        let mut best = a.get(k, k).abs();
        for i in (k + 1)..n {
            let v = a.get(i, k).abs();
            if v > best {
                best = v;
                piv = i;
            }
        }
        if best <= 0.0 {
            return None;
        }
        if piv != k {
            for j in 0..n {
                let t = a.get(k, j);
                a.set(k, j, a.get(piv, j));
                a.set(piv, j, t);
            }
            b.swap(k, piv);
        }
        let akk = a.get(k, k);
        for i in (k + 1)..n {
            let f = a.get(i, k) / akk;
            if f == 0.0 {
                continue;
            }
            a.set(i, k, 0.0);
            for j in (k + 1)..n {
                a.add(i, j, -f * a.get(k, j));
            }
            b[i] -= f * b[k];
        }
    }
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in (i + 1)..n {
            s -= a.get(i, j) * b[j];
        }
        b[i] = s / a.get(i, i);
    }
    Some(())
}

/// Inverse of a small dense matrix by Gauss-Jordan.
pub fn invert(a: &Dense) -> Option<Dense> {
    let n = a.n;
    let mut m = Dense::zeros(n);
    let mut inv = Dense::zeros(n);
    for i in 0..n {
        for j in 0..n {
            m.set(i, j, a.get(i, j));
        }
        inv.set(i, i, 1.0);
    }
    for k in 0..n {
        let mut piv = k;
        let mut best = m.get(k, k).abs();
        for i in (k + 1)..n {
            if m.get(i, k).abs() > best {
                best = m.get(i, k).abs();
                piv = i;
            }
        }
        if best <= 0.0 {
            return None;
        }
        if piv != k {
            for j in 0..n {
                let t = m.get(k, j);
                m.set(k, j, m.get(piv, j));
                m.set(piv, j, t);
                let t = inv.get(k, j);
                inv.set(k, j, inv.get(piv, j));
                inv.set(piv, j, t);
            }
        }
        let d = m.get(k, k);
        for j in 0..n {
            m.set(k, j, m.get(k, j) / d);
            inv.set(k, j, inv.get(k, j) / d);
        }
        for i in 0..n {
            if i == k {
                continue;
            }
            let f = m.get(i, k);
            if f == 0.0 {
                continue;
            }
            for j in 0..n {
                m.add(i, j, -f * m.get(k, j));
                inv.add(i, j, -f * inv.get(k, j));
            }
        }
    }
    Some(inv)
}

/// Cholesky factor `L` with `A = L Lᵀ` (lower triangular). `None` if not SPD.
pub fn cholesky(a: &Dense) -> Option<Dense> {
    let n = a.n;
    let mut l = Dense::zeros(n);
    for i in 0..n {
        for j in 0..=i {
            let mut s = a.get(i, j);
            for k in 0..j {
                s -= l.get(i, k) * l.get(j, k);
            }
            if i == j {
                if s <= 0.0 {
                    return None;
                }
                l.set(i, j, s.sqrt());
            } else {
                l.set(i, j, s / l.get(j, j));
            }
        }
    }
    Some(l)
}

/// Eigenvalues and eigenvectors of a real symmetric matrix by the cyclic
/// Jacobi method. Returns `(eigenvalues, eigenvectors)` with eigenvectors as
/// columns of the returned `Dense`.
pub fn jacobi_eigen(a: &Dense) -> (Vec<f64>, Dense) {
    let n = a.n;
    let mut m = a.clone();
    let mut v = Dense::zeros(n);
    for i in 0..n {
        v.set(i, i, 1.0);
    }
    for _sweep in 0..100 {
        let mut off = 0.0;
        for p in 0..n {
            for q in (p + 1)..n {
                off += m.get(p, q) * m.get(p, q);
            }
        }
        if off < 1e-30 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = m.get(p, q);
                if apq.abs() < 1e-300 {
                    continue;
                }
                let app = m.get(p, p);
                let aqq = m.get(q, q);
                let theta = 0.5 * (aqq - app) / apq;
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let mkp = m.get(k, p);
                    let mkq = m.get(k, q);
                    m.set(k, p, c * mkp - s * mkq);
                    m.set(k, q, s * mkp + c * mkq);
                }
                for k in 0..n {
                    let mpk = m.get(p, k);
                    let mqk = m.get(q, k);
                    m.set(p, k, c * mpk - s * mqk);
                    m.set(q, k, s * mpk + c * mqk);
                }
                for k in 0..n {
                    let vkp = v.get(k, p);
                    let vkq = v.get(k, q);
                    v.set(k, p, c * vkp - s * vkq);
                    v.set(k, q, s * vkp + c * vkq);
                }
            }
        }
    }
    let vals = (0..n).map(|i| m.get(i, i)).collect();
    (vals, v)
}

/// Solve the symmetric generalised eigenproblem `K x = λ M x` for `M` SPD.
/// Returns `(λ, X)` sorted ascending with eigenvectors as columns of `X`,
/// normalised so that `Xᵀ M X = I`.
pub fn generalized_eigen(k: &Dense, m: &Dense) -> Option<(Vec<f64>, Dense)> {
    let l = cholesky(m)?;
    // B = L^-1 K L^-T.
    let n = k.n;
    // Solve L Y = K  -> Y = L^-1 K.
    let mut y = Dense::zeros(n);
    for col in 0..n {
        let mut rhs: Vec<f64> = (0..n).map(|i| k.get(i, col)).collect();
        // forward substitution
        for i in 0..n {
            let mut s = rhs[i];
            for j in 0..i {
                s -= l.get(i, j) * rhs[j];
            }
            rhs[i] = s / l.get(i, i);
        }
        for i in 0..n {
            y.set(i, col, rhs[i]);
        }
    }
    // B = Y L^-T  (i.e. B = L^-1 K L^-T).
    let mut b = Dense::zeros(n);
    for row in 0..n {
        let mut rhs: Vec<f64> = (0..n).map(|j| y.get(row, j)).collect();
        // solve L z = y_row  => z = L^-1 y_row, then B[row][j] = z[j]
        for i in 0..n {
            let mut s = rhs[i];
            for j in 0..i {
                s -= l.get(i, j) * rhs[j];
            }
            rhs[i] = s / l.get(i, i);
        }
        for j in 0..n {
            b.set(row, j, rhs[j]);
        }
    }
    // symmetrise to kill round-off asymmetry
    for i in 0..n {
        for j in (i + 1)..n {
            let avg = 0.5 * (b.get(i, j) + b.get(j, i));
            b.set(i, j, avg);
            b.set(j, i, avg);
        }
    }
    let (vals, z) = jacobi_eigen(&b);
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&i, &j| vals[i].partial_cmp(&vals[j]).unwrap());
    let mut vals_sorted = vec![0.0; n];
    let mut x = Dense::zeros(n);
    for (newc, &oldc) in idx.iter().enumerate() {
        vals_sorted[newc] = vals[oldc];
        // x = L^-T z_column
        let mut rhs: Vec<f64> = (0..n).map(|i| z.get(i, oldc)).collect();
        // solve L^T x = z: back substitution
        for i in (0..n).rev() {
            let mut s = rhs[i];
            for j in (i + 1)..n {
                s -= l.get(j, i) * rhs[j];
            }
            rhs[i] = s / l.get(i, i);
        }
        for i in 0..n {
            x.set(i, newc, rhs[i]);
        }
    }
    Some((vals_sorted, x))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn close(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "{a:e} vs {b:e}");
    }
    #[test]
    fn solve_and_invert() {
        let mut a = Dense {
            n: 2,
            a: vec![3.0, 1.0, 1.0, 2.0],
        };
        let mut b = vec![5.0, 5.0];
        solve(&mut a, &mut b).unwrap();
        close(b[0], 1.0, 1e-14);
        close(b[1], 2.0, 1e-14);
        let a = Dense {
            n: 2,
            a: vec![3.0, 1.0, 1.0, 2.0],
        };
        let inv = invert(&a).unwrap();
        close(inv.get(0, 0), 2.0 / 5.0, 1e-14);
        close(inv.get(0, 1), -1.0 / 5.0, 1e-14);
    }
    #[test]
    fn generalised_eigen_matches_hand_case() {
        // diag(1,2) x = λ I x  -> λ = 1, 2.
        let k = Dense {
            n: 2,
            a: vec![1.0, 0.0, 0.0, 2.0],
        };
        let m = Dense {
            n: 2,
            a: vec![1.0, 0.0, 0.0, 1.0],
        };
        let (vals, _x) = generalized_eigen(&k, &m).unwrap();
        close(vals[0], 1.0, 1e-12);
        close(vals[1], 2.0, 1e-12);
    }
}
