//! Restarted GMRES for small dense-free linear systems (the smooth forcing
//! extension's boundary system is not symmetric).

/// Solve A x = b with GMRES(restart). `x` holds the initial guess and the
/// result. Returns (matrix-vector products, final relative residual).
pub fn gmres(
    mut apply: impl FnMut(&[f64], &mut [f64]),
    b: &[f64],
    x: &mut [f64],
    restart: usize,
    tol: f64,
    max_matvecs: usize,
) -> (usize, f64) {
    let n = b.len();
    let dot = |u: &[f64], v: &[f64]| u.iter().zip(v).map(|(p, q)| p * q).sum::<f64>();
    let bnorm = dot(b, b).sqrt().max(1e-300);
    let mut ax = vec![0.0; n];
    let mut matvecs = 0;
    loop {
        apply(x, &mut ax);
        matvecs += 1;
        let r: Vec<f64> = b.iter().zip(&ax).map(|(p, q)| p - q).collect();
        let beta = dot(&r, &r).sqrt();
        if beta / bnorm < tol || matvecs >= max_matvecs {
            return (matvecs, beta / bnorm);
        }
        // Arnoldi with modified Gram–Schmidt and Givens rotations
        let mut v: Vec<Vec<f64>> = vec![r.iter().map(|q| q / beta).collect()];
        let mut h = vec![vec![0.0; restart]; restart + 1];
        let (mut cs, mut sn) = (vec![0.0; restart], vec![0.0; restart]);
        let mut g = vec![0.0; restart + 1];
        g[0] = beta;
        let mut k_used = 0;
        for k in 0..restart {
            let mut wv = vec![0.0; n];
            apply(&v[k], &mut wv);
            matvecs += 1;
            for (j, vj) in v.iter().enumerate() {
                h[j][k] = dot(&wv, vj);
                wv.iter_mut().zip(vj).for_each(|(p, q)| *p -= h[j][k] * q);
            }
            h[k + 1][k] = dot(&wv, &wv).sqrt();
            for j in 0..k {
                let t = cs[j] * h[j][k] + sn[j] * h[j + 1][k];
                h[j + 1][k] = -sn[j] * h[j][k] + cs[j] * h[j + 1][k];
                h[j][k] = t;
            }
            let den = (h[k][k] * h[k][k] + h[k + 1][k] * h[k + 1][k]).sqrt().max(1e-300);
            cs[k] = h[k][k] / den;
            sn[k] = h[k + 1][k] / den;
            h[k][k] = den;
            h[k + 1][k] = 0.0;
            g[k + 1] = -sn[k] * g[k];
            g[k] *= cs[k];
            k_used = k + 1;
            let hk1 = dot(&wv, &wv).sqrt();
            if (g[k + 1] / bnorm).abs() < tol || matvecs >= max_matvecs || hk1 < 1e-300 {
                break;
            }
            v.push(wv.iter().map(|q| q / hk1).collect());
        }
        // back substitution and update
        let mut y = vec![0.0; k_used];
        for i in (0..k_used).rev() {
            let s: f64 = (i + 1..k_used).map(|j| h[i][j] * y[j]).sum();
            y[i] = (g[i] - s) / h[i][i];
        }
        for (j, yj) in y.iter().enumerate() {
            x.iter_mut().zip(&v[j]).for_each(|(p, q)| *p += yj * q);
        }
        if matvecs >= max_matvecs {
            apply(x, &mut ax);
            let r2: f64 = b.iter().zip(&ax).map(|(p, q)| (p - q) * (p - q)).sum();
            return (matvecs + 1, r2.sqrt() / bnorm);
        }
    }
}
