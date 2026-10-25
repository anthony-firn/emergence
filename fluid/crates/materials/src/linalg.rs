//! Small dense linear algebra for the lattice rungs. Everything here is
//! tiny (at most a few dozen rows), so a cyclic Jacobi eigensolver and a
//! pivoted Gaussian solve are enough; there are no external dependencies.

/// Eigenvalues of a real symmetric matrix given row-major flat.
/// Returns them in the order the Jacobi sweeps leave them (unsorted).
pub fn jacobi_eigenvalues(a: &[f64], n: usize) -> Vec<f64> {
    assert_eq!(a.len(), n * n);
    let mut a = a.to_vec();
    for _ in 0..100 {
        let mut off = 0.0;
        for p in 0..n {
            for q in p + 1..n {
                off += a[p * n + q] * a[p * n + q];
            }
        }
        if off.sqrt() < 1e-14 {
            break;
        }
        for p in 0..n.saturating_sub(1) {
            for q in p + 1..n {
                let apq = a[p * n + q];
                if apq == 0.0 {
                    continue;
                }
                let app = a[p * n + p];
                let aqq = a[q * n + q];
                let theta = (aqq - app) / (2.0 * apq);
                let t = if theta >= 0.0 {
                    1.0 / (theta + (theta * theta + 1.0).sqrt())
                } else {
                    -1.0 / (-theta + (theta * theta + 1.0).sqrt())
                };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let akp = a[k * n + p];
                    let akq = a[k * n + q];
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let apk = a[p * n + k];
                    let aqk = a[q * n + k];
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
            }
        }
    }
    (0..n).map(|i| a[i * n + i]).collect()
}

/// Eigenvalues of a real symmetric matrix, ascending.
pub fn sym_eigenvalues(a: &[f64], n: usize) -> Vec<f64> {
    let mut e = jacobi_eigenvalues(a, n);
    e.sort_by(|x, y| x.total_cmp(y));
    e
}

/// Eigen-decomposition of a real symmetric matrix: eigenvalues ascending
/// and, for each, a unit eigenvector. `a` is row-major flat; the returned
/// eigenvector `vecs[i]` has `n` components and belongs to `eval[i]`.
pub fn sym_eigen(a: &[f64], n: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    assert_eq!(a.len(), n * n);
    let mut m = a.to_vec();
    // Columns of `v` accumulate the eigenvectors (V <- V J).
    let mut v = vec![vec![0.0; n]; n];
    for i in 0..n {
        v[i][i] = 1.0;
    }
    for _ in 0..100 {
        let mut off = 0.0;
        for p in 0..n {
            for q in p + 1..n {
                off += m[p * n + q] * m[p * n + q];
            }
        }
        if off.sqrt() < 1e-14 {
            break;
        }
        for p in 0..n.saturating_sub(1) {
            for q in p + 1..n {
                let mpq = m[p * n + q];
                if mpq == 0.0 {
                    continue;
                }
                let app = m[p * n + p];
                let aqq = m[q * n + q];
                let theta = (aqq - app) / (2.0 * mpq);
                let t = if theta >= 0.0 {
                    1.0 / (theta + (theta * theta + 1.0).sqrt())
                } else {
                    -1.0 / (-theta + (theta * theta + 1.0).sqrt())
                };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let mkp = m[k * n + p];
                    let mkq = m[k * n + q];
                    m[k * n + p] = c * mkp - s * mkq;
                    m[k * n + q] = s * mkp + c * mkq;
                }
                for k in 0..n {
                    let mpk = m[p * n + k];
                    let mqk = m[q * n + k];
                    m[p * n + k] = c * mpk - s * mqk;
                    m[q * n + k] = s * mpk + c * mqk;
                }
                for k in 0..n {
                    let vkp = v[k][p];
                    let vkq = v[k][q];
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut pairs: Vec<(f64, Vec<f64>)> = (0..n)
        .map(|i| (m[i * n + i], (0..n).map(|k| v[k][i]).collect::<Vec<f64>>()))
        .collect();
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
    pairs.into_iter().unzip()
}

/// Solve `a x = b` for a small dense system by partial-pivoted Gauss
/// elimination. `a` is row-major and is consumed.
pub fn solve(mut a: Vec<f64>, mut b: Vec<f64>, n: usize) -> Vec<f64> {
    assert_eq!(a.len(), n * n);
    assert_eq!(b.len(), n);
    for col in 0..n {
        let mut piv = col;
        let mut best = a[col * n + col].abs();
        for r in col + 1..n {
            let v = a[r * n + col].abs();
            if v > best {
                best = v;
                piv = r;
            }
        }
        assert!(best > 1e-300, "singular internal-coordinate Hessian");
        if piv != col {
            for k in 0..n {
                a.swap(col * n + k, piv * n + k);
            }
            b.swap(col, piv);
        }
        let d = a[col * n + col];
        for r in col + 1..n {
            let factor = a[r * n + col] / d;
            if factor == 0.0 {
                continue;
            }
            for k in col..n {
                a[r * n + k] -= factor * a[col * n + k];
            }
            b[r] -= factor * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let mut s = b[r];
        for k in r + 1..n {
            s -= a[r * n + k] * x[k];
        }
        x[r] = s / a[r * n + r];
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jacobi_recovers_known_spectrum() {
        let a = [2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0];
        let e = sym_eigenvalues(&a, 3);
        assert!((e[0] - 1.0).abs() < 1e-12, "{e:?}");
        assert!((e[1] - 1.0).abs() < 1e-12, "{e:?}");
        assert!((e[2] - 4.0).abs() < 1e-12, "{e:?}");
    }

    #[test]
    fn solve_matches_residual() {
        let a = vec![3.0, 1.0, 1.0, 2.0];
        let x = solve(a.clone(), vec![5.0, 5.0], 2);
        // 3x+y=5, x+2y=5 -> x=1, y=2
        assert!((x[0] - 1.0).abs() < 1e-14 && (x[1] - 2.0).abs() < 1e-14);
    }
}
