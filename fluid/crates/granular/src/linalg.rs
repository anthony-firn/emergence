//! Tiny dense linear algebra (partial-pivoted Gauss solve) for the small
//! least-squares systems used to fit the `mu(I)` closure and the small
//! boundary-integral checks. No external dependency.

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
        assert!(best > 1e-300, "singular system");
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

/// Least-squares solve of the (small) over-determined system `A x = b`
/// built from normal equations. `cols` is the number of columns of `A`.
pub fn least_squares(a: &[Vec<f64>], b: &[f64], cols: usize) -> Vec<f64> {
    let mut ata = vec![0.0; cols * cols];
    let mut atb = vec![0.0; cols];
    for (row, bi) in a.iter().zip(b) {
        for p in 0..cols {
            atb[p] += row[p] * bi;
            for q in 0..cols {
                ata[p * cols + q] += row[p] * row[q];
            }
        }
    }
    solve(ata, atb, cols)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solve_residual() {
        let x = solve(vec![3.0, 1.0, 1.0, 2.0], vec![5.0, 5.0], 2);
        assert!((x[0] - 1.0).abs() < 1e-13 && (x[1] - 2.0).abs() < 1e-13);
    }
}
