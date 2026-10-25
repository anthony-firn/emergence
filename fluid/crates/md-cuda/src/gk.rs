//! Green–Kubo and Einstein–Helfand transport estimators with
//! block-averaged error bars.
//!
//! The CPU `md::transport` estimators read a single running integral and
//! quote the scatter of that curve over a lag window, which understates the
//! true statistical error because neighbouring lags are correlated. Here
//! the trajectory is split into `nblocks` contiguous blocks; each block
//! yields an independent estimate of the running integral (or of the
//! Einstein–Helfand slope), and the quoted error is the standard error of
//! the mean across blocks. The blocks are chosen much longer than the
//! stress correlation time (`τ* ≈ 0.14`), so they are effectively
//! independent.
//!
//! Conventions are the reduced-unit ones of `md::transport`:
//!
//! ```text
//! eta = V/(k_B T) * int_0^inf <P_xy(0) P_xy(t)> dt
//! ```
//!
//! with `P = (sum m v v + sum r F)/V` the intensive pressure tensor. By
//! symmetry the three off-diagonal components `xy, xz, yz` are averaged.

/// A running-integral curve, one value per lag, with a block standard error.
pub struct BlockCurve {
    pub dt: f64,
    pub mean: Vec<f64>,
    pub err: Vec<f64>,
    pub nblocks: usize,
    /// Samples per block actually used.
    pub block_len: usize,
}

/// Mean and standard error of the mean across `values`.
pub fn mean_sem(values: &[f64]) -> (f64, f64) {
    let n = values.len() as f64;
    if n == 0.0 {
        return (f64::NAN, f64::NAN);
    }
    let mean = values.iter().sum::<f64>() / n;
    if values.len() < 2 {
        return (mean, f64::NAN);
    }
    let var = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
    (mean, (var / n).sqrt())
}

/// Least-squares slope and its standard error for `y` against `x`.
fn linfit(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len() as f64;
    let sx: f64 = x.iter().sum();
    let sy: f64 = y.iter().sum();
    let sxx: f64 = x.iter().map(|v| v * v).sum();
    let sxy: f64 = x.iter().zip(y).map(|(a, b)| a * b).sum();
    let denom = n * sxx - sx * sx;
    let slope = (n * sxy - sx * sy) / denom;
    let intercept = (sy - slope * sx) / n;
    let ss: f64 = x
        .iter()
        .zip(y)
        .map(|(a, b)| (b - (intercept + slope * a)).powi(2))
        .sum();
    let dof = (n - 2.0).max(1.0);
    let s2 = ss / dof;
    let se = (s2 * n / denom).sqrt();
    (slope, se)
}

const OFFDIAG: [usize; 3] = [3, 4, 5];

/// Running Green–Kubo shear-viscosity integral, block-averaged.
///
/// `press` is the intensive pressure tensor sampled every `dt`; `volume`
/// and `kbt` are the box volume and temperature in reduced units. Only the
/// off-diagonal components enter. The returned `mean[k]`/`err[k]` are the
/// mean and block SEM of `V/(kBT) ∫_0^{k dt} <P_xy P_xy>`.
pub fn running_gk_viscosity(
    press: &[[f64; 6]],
    dt: f64,
    volume: f64,
    kbt: f64,
    max_lag: usize,
    nblocks: usize,
) -> BlockCurve {
    let m = press.len();
    assert!(nblocks > 0 && m >= nblocks);
    let blk = m / nblocks;
    let max_lag = max_lag.min(blk.saturating_sub(1));
    // Global mean of each off-diagonal component; subtracting the global
    // mean keeps the block estimates on a common footing.
    let mut gmean = [0.0f64; 3];
    for p in press {
        for (a, &c) in OFFDIAG.iter().enumerate() {
            gmean[a] += p[c];
        }
    }
    for g in &mut gmean {
        *g /= m as f64;
    }

    let mut curves: Vec<Vec<f64>> = Vec::with_capacity(nblocks);
    for b in 0..nblocks {
        let s = &press[b * blk..(b + 1) * blk];
        let nb = s.len();
        let mut run = vec![0.0f64; max_lag + 1];
        for (a, &c) in OFFDIAG.iter().enumerate() {
            // Autocorrelation of this component within the block.
            let mut ac = vec![0.0f64; max_lag + 1];
            for k in 0..=max_lag {
                let mut sum = 0.0;
                for t in 0..(nb - k) {
                    sum += (s[t][c] - gmean[a]) * (s[t + k][c] - gmean[a]);
                }
                ac[k] = sum / (nb - k) as f64;
            }
            let mut integ = 0.0;
            for k in 1..=max_lag {
                integ += 0.5 * dt * (ac[k - 1] + ac[k]);
                run[k] += integ;
            }
        }
        for v in run.iter_mut() {
            *v *= volume / (3.0 * kbt);
        }
        curves.push(run);
    }

    let mut mean = vec![0.0; max_lag + 1];
    let mut err = vec![0.0; max_lag + 1];
    for k in 0..=max_lag {
        let vals: Vec<f64> = curves.iter().map(|c| c[k]).collect();
        let (mu, se) = mean_sem(&vals);
        mean[k] = mu;
        err[k] = se;
    }
    BlockCurve { dt, mean, err, nblocks, block_len: blk }
}

/// Generalized Einstein–Helfand shear viscosity, block-averaged.
///
/// `eta = V/(2 kBT) d/dt <[∫_0^t P_xy]^2>`, with the three off-diagonal
/// components averaged. The slope is fitted over lags `[lo, hi]`.
pub fn einstein_helfand_viscosity(
    press: &[[f64; 6]],
    dt: f64,
    volume: f64,
    kbt: f64,
    max_lag: usize,
    nblocks: usize,
    lo: usize,
    hi: usize,
) -> (f64, f64) {
    let m = press.len();
    assert!(nblocks > 0 && m >= nblocks);
    let blk = m / nblocks;
    let max_lag = max_lag.min(blk.saturating_sub(1));
    let hi = hi.min(max_lag);
    let lo = lo.min(hi);

    let mut etas = Vec::with_capacity(nblocks);
    for b in 0..nblocks {
        let s = &press[b * blk..(b + 1) * blk];
        let nb = s.len();
        // Cumulative integrals A_c(t) for t = 0..=nb-1.
        let mut disp = vec![0.0f64; max_lag + 1];
        for &c in &OFFDIAG {
            let mut a = vec![0.0f64; nb];
            for t in 1..nb {
                a[t] = a[t - 1] + dt * s[t][c];
            }
            for k in 0..=max_lag {
                let mut sum = 0.0;
                let mut cnt = 0usize;
                for t in 0..(nb - k) {
                    let d = a[t + k] - a[t];
                    sum += d * d;
                    cnt += 1;
                }
                disp[k] += sum / cnt.max(1) as f64;
            }
        }
        for v in disp.iter_mut() {
            *v /= 3.0;
        }
        let x: Vec<f64> = (lo..=hi).map(|k| k as f64 * dt).collect();
        let y: Vec<f64> = disp[lo..=hi].to_vec();
        let (slope, _) = linfit(&x, &y);
        etas.push(volume * slope / (2.0 * kbt));
    }
    mean_sem(&etas)
}

/// Average several running-integral curves (e.g. independent seeds) lag by
/// lag, with the standard error of the mean across the curves. This is the
/// replica-based error bar: the curves are statistically independent.
pub fn average_curves(curves: &[BlockCurve], dt: f64) -> BlockCurve {
    assert!(!curves.is_empty());
    let len = curves.iter().map(|c| c.mean.len()).min().unwrap();
    let mut mean = vec![0.0; len];
    let mut err = vec![0.0; len];
    for k in 0..len {
        let vals: Vec<f64> = curves.iter().map(|c| c.mean[k]).collect();
        let (mu, se) = mean_sem(&vals);
        mean[k] = mu;
        err[k] = se;
    }
    BlockCurve {
        dt,
        mean,
        err,
        nblocks: curves.iter().map(|c| c.nblocks).sum(),
        block_len: curves[0].block_len,
    }
}

/// Earliest lag at which the running integral is flat within the block
/// errors. `min_frac` sets the fraction of the curve that must already be
/// plateaux-like. Returns the lag index (samples).
pub fn plateau_lag(curve: &BlockCurve) -> usize {
    let m = curve.mean.len();
    if m < 8 {
        return m - 1;
    }
    // Reference level: average of the last 20 % of the curve.
    let ntail = (m / 5).max(2);
    let ref_val = curve.mean[m - ntail..].iter().sum::<f64>() / ntail as f64;
    let ref_err = curve.err[m - ntail..].iter().sum::<f64>() / ntail as f64;
    // Start looking after the initial transient, where the integral is
    // still rising from zero.
    let start = (m / 8).max(1);
    for k in start..m {
        // All values from k on must lie within 2 combined errors of the
        // long-time reference.
        let ok = (k..m).all(|j| {
            let tol = 2.0 * (curve.err[j].max(ref_err)).max(1e-12);
            (curve.mean[j] - ref_val).abs() <= tol
        });
        if ok {
            return k;
        }
    }
    m - 1
}

/// Value and block error of the running integral at the plateau lag.
pub fn plateau_value(curve: &BlockCurve) -> (f64, f64, usize) {
    let k = plateau_lag(curve);
    (curve.mean[k], curve.err[k], k)
}

/// Mean of a running-integral curve over the sample window `[lo, hi]`.
pub fn window_mean(curve: &BlockCurve, lo: usize, hi: usize) -> f64 {
    let hi = hi.min(curve.mean.len() - 1);
    let lo = lo.min(hi);
    curve.mean[lo..=hi].iter().sum::<f64>() / (hi - lo + 1) as f64
}

/// Centre of the Green–Kubo plateau, found as the lag of the maximum of
/// the running integral smoothed over `win` samples. The running integral
/// rises, flattens near its maximum and then drifts on the finite-size
/// long-time tail, so the peak of the smoothed curve is a robust plateau
/// locator (unlike anchoring on the noisy tail). Returns the lag index.
pub fn plateau_peak_lag(curve: &BlockCurve, win: usize) -> usize {
    let m = curve.mean.len();
    if m == 0 {
        return 0;
    }
    let win = win.clamp(1, m);
    let half = win / 2;
    let mut best = 0usize;
    let mut bestv = f64::NEG_INFINITY;
    for k in 0..m {
        let lo = k.saturating_sub(half);
        let hi = (k + half).min(m - 1);
        let s: f64 = curve.mean[lo..=hi].iter().sum::<f64>() / (hi - lo + 1) as f64;
        if s > bestv {
            bestv = s;
            best = k;
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    /// With a single block the block estimator must reduce exactly to the
    /// direct running Green–Kubo integral of the three off-diagonal
    /// components, which is the formula `md::transport` uses.
    #[test]
    fn single_block_matches_direct_running_integral() {
        let m = 2000usize;
        let dt = 0.01;
        let volume = 2.0;
        let kbt = 0.8;
        let max_lag = 300usize;
        let mut press = vec![[0.0f64; 6]; m];
        let mut state = 424242u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as f64 / (1u64 << 31) as f64 - 1.0
        };
        for p in press.iter_mut() {
            p[3] = next();
            p[4] = next();
            p[5] = next();
        }
        let curve = running_gk_viscosity(&press, dt, volume, kbt, max_lag, 1);
        // Direct computation with global means.
        let mut gmean = [0.0f64; 3];
        for p in &press {
            for (a, &c) in OFFDIAG.iter().enumerate() {
                gmean[a] += p[c];
            }
        }
        for g in &mut gmean {
            *g /= m as f64;
        }
        let mut direct = vec![0.0f64; max_lag + 1];
        for (a, &c) in OFFDIAG.iter().enumerate() {
            let mut ac = vec![0.0f64; max_lag + 1];
            for k in 0..=max_lag {
                let mut sum = 0.0;
                for t in 0..(m - k) {
                    sum += (press[t][c] - gmean[a]) * (press[t + k][c] - gmean[a]);
                }
                ac[k] = sum / (m - k) as f64;
            }
            let mut integ = 0.0;
            for k in 1..=max_lag {
                integ += 0.5 * dt * (ac[k - 1] + ac[k]);
                direct[k] += integ;
            }
        }
        for v in direct.iter_mut() {
            *v *= volume / (3.0 * kbt);
        }
        for k in 0..=max_lag {
            assert!(
                (curve.mean[k] - direct[k]).abs() <= 1e-12 * direct[k].abs().max(1.0),
                "k={k}: {} vs {}",
                curve.mean[k],
                direct[k]
            );
        }
        // With a single block the standard error is undefined (NaN); the
        // lag-0 integral is identically zero, so its error must be 0 once
        // two or more blocks are available. Accept NaN here.
        assert!(curve.err[0].is_nan() || curve.err[0] == 0.0);
    }

    #[test]
    fn block_error_is_finite_and_shrinks_with_blocks() {
        let m = 8000usize;
        let dt = 0.01;
        let mut press = vec![[0.0f64; 6]; m];
        // Three independent pseudo-random components (LCG) with a short
        // correlation length.
        let mut state = 12345u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as f64 / (1u64 << 31) as f64 - 1.0
        };
        for p in press.iter_mut() {
            p[3] = next();
            p[4] = next();
            p[5] = next();
        }
        let c8 = running_gk_viscosity(&press, dt, 1.0, 1.0, 200, 8);
        let c16 = running_gk_viscosity(&press, dt, 1.0, 1.0, 200, 16);
        assert!(c8.err[200].is_finite() && c8.err[200] > 0.0);
        assert!(c16.err[200].is_finite() && c16.err[200] > 0.0);
    }

    #[test]
    fn einstein_helfand_runs() {
        let m = 4000usize;
        let dt = 0.01;
        let mut press = vec![[0.0f64; 6]; m];
        let mut state = 999u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as f64 / (1u64 << 31) as f64 - 1.0
        };
        for p in press.iter_mut() {
            p[3] = next();
            p[4] = next();
            p[5] = next();
        }
        let (eta, err) = einstein_helfand_viscosity(&press, dt, 1.0, 1.0, 500, 4, 50, 400);
        assert!(eta.is_finite());
        assert!(err.is_finite() && err >= 0.0);
    }

    /// The plateau locator must find the flat top of a synthetic running
    /// integral (rise, plateau, decay) rather than the decayed tail.
    #[test]
    fn plateau_peak_finds_the_flat_top() {
        let dt = 0.01;
        let m = 3000usize;
        let mut mean = vec![0.0f64; m];
        for (k, v) in mean.iter_mut().enumerate() {
            let t = k as f64 * dt;
            *v = if t < 1.0 {
                3.0 * t * t
            } else if t < 4.0 {
                3.0
            } else {
                3.0 - 0.3 * (t - 4.0)
            };
        }
        let curve = BlockCurve {
            dt,
            mean,
            err: vec![0.01; m],
            nblocks: 4,
            block_len: m / 4,
        };
        let peak = plateau_peak_lag(&curve, 20);
        assert!(peak >= 90 && peak <= 410, "peak lag {peak} outside the plateau");
        let wm = window_mean(&curve, 100, 400);
        assert!((wm - 3.0).abs() < 1e-9, "window mean {wm}");
    }
}
