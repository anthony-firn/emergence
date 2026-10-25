//! Time-correlation analysis: Green–Kubo transport coefficients and
//! self-diffusion from the mean-squared displacement and the velocity
//! autocorrelation.
//!
//! Conventions (reduced units, `k_B = 1`):
//!
//! * shear viscosity `eta = (V/(kT)) int_0^inf <P_xy(0) P_xy(t)> dt`,
//!   with `P` the *intensive* pressure/stress tensor `(sum m v v + sum r F)/V`;
//! * thermal conductivity `kappa = (1/(V T^2)) int_0^inf <J_x(0) J_x(t)> dt`,
//!   with `J` the volume-integrated heat flux of
//!   [`crate::system::System::heat_flux`];
//! * self-diffusion `D = lim MSD/(6t) = (1/3) int <v(0).v(t)> dt`.
//!
//! Green–Kubo integrals are read from the plateau of the running integral.
//! For short trajectories the running integral has not converged: the
//! reported scatter is the standard deviation of the running integral over
//! the chosen lag window, which is a lower bound on the true statistical
//! error (independent time origins are not available in a short run). This
//! is stated honestly wherever it is used.
//!
//! The finite-size correction to `D` is that of Yeh & Hummer (2004),
//! `D_0 = D_PBC - T xi/(6 pi eta L)` with `xi = 2.837297` (reduced units).

/// A sampled time series and its normalised autocorrelation.
pub struct Correlation {
    /// Sampling interval.
    pub dt: f64,
    /// `c[k]` for lags `0..=max_lag`.
    pub c: Vec<f64>,
    /// Number of origin pairs contributing to each lag.
    pub count: Vec<usize>,
}

/// Mean-subtracted autocorrelation of a scalar series.
pub fn autocorrelate(series: &[f64], max_lag: usize) -> Correlation {
    let m = series.len();
    let mean = if m == 0 {
        0.0
    } else {
        series.iter().sum::<f64>() / m as f64
    };
    let max_lag = max_lag.min(m.saturating_sub(1));
    let mut c = vec![0.0; max_lag + 1];
    let mut count = vec![0usize; max_lag + 1];
    for (k, (ck, nk)) in c.iter_mut().zip(count.iter_mut()).enumerate() {
        let mut acc = 0.0;
        for t in 0..(m - k) {
            acc += (series[t] - mean) * (series[t + k] - mean);
        }
        *ck = acc;
        *nk = m - k;
    }
    for (ck, nk) in c.iter_mut().zip(&count) {
        if *nk > 0 {
            *ck /= *nk as f64;
        }
    }
    Correlation { dt: 1.0, c, count }
}

/// Trapezoidal running integral `int_0^{t_k} C dt` of an autocorrelation.
pub fn running_integral(corr: &Correlation) -> Vec<f64> {
    let c = &corr.c;
    let mut out = vec![0.0; c.len()];
    for k in 1..c.len() {
        out[k] = out[k - 1] + 0.5 * corr.dt * (c[k - 1] + c[k]);
    }
    out
}

/// Running Green–Kubo integral `I(t_k) = int_0^{t_k} <x(0)x(t)> dt`.
pub fn running_gk(series: &[f64], dt: f64, max_lag: usize) -> Vec<f64> {
    let mut corr = autocorrelate(series, max_lag);
    corr.dt = dt;
    running_integral(&corr)
}

/// Green–Kubo integral of a scalar series at a fixed lag cutoff.
pub fn gk_integral(series: &[f64], dt: f64, max_lag: usize) -> f64 {
    *running_gk(series, dt, max_lag).last().unwrap()
}

/// Mean and scatter of a running integral over the lag window `[lo, hi]`.
/// The scatter is the standard deviation of the running curve over the
/// window; it is a rough lower bound on the true statistical error.
pub fn plateau(running: &[f64], lo: usize, hi: usize) -> (f64, f64) {
    let hi = hi.min(running.len() - 1);
    let lo = lo.min(hi);
    let s = &running[lo..=hi];
    let n = s.len() as f64;
    let mean = s.iter().sum::<f64>() / n;
    if s.len() < 2 {
        return (mean, f64::NAN);
    }
    let var = s.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
    (mean, var.sqrt())
}

/// Average the running integrals of several series.
fn average_running(series: &[Vec<f64>], dt: f64, max_lag: usize) -> Vec<f64> {
    let mut acc: Vec<f64> = Vec::new();
    for (i, s) in series.iter().enumerate() {
        let r = running_gk(s, dt, max_lag);
        if i == 0 {
            acc = r;
        } else {
            for (a, b) in acc.iter_mut().zip(&r) {
                *a += b;
            }
        }
    }
    let inv = 1.0 / series.len() as f64;
    for a in &mut acc {
        *a *= inv;
    }
    acc
}

/// Shear viscosity from the three off-diagonal stress components.
///
/// `stress[i]` holds the *intensive* pressure tensor
/// `[xx, yy, zz, xy, xz, yz] = Pi/V` at sample `i`. The plateau window
/// `[lo, hi]` is in sample indices; `lo*dt` to `hi*dt` should span the
/// region where the running integral is flat.
pub fn green_kubo_viscosity(
    stress: &[[f64; 6]],
    dt: f64,
    volume: f64,
    kbt: f64,
    max_lag: usize,
    lo: usize,
    hi: usize,
) -> (f64, f64) {
    let comps = [3usize, 4, 5];
    let series: Vec<Vec<f64>> = comps
        .iter()
        .map(|&c| stress.iter().map(|s| s[c]).collect())
        .collect();
    let run = average_running(&series, dt, max_lag);
    let scaled: Vec<f64> = run.iter().map(|i| volume * i / kbt).collect();
    plateau(&scaled, lo, hi)
}

/// Generalized Einstein (Helfand) estimate of the shear viscosity:
/// `eta = V/(2 kT) d<[int_0^t P_xy]^2>/dt` in the linear regime. This is
/// the same integral as [`green_kubo_viscosity`] but read from the slope of
/// a smoother displacement function; it is less sensitive to the noisy
/// long-time tail of the stress autocorrelation.
pub fn green_kubo_viscosity_einstein(
    stress: &[[f64; 6]],
    dt: f64,
    volume: f64,
    kbt: f64,
    max_lag: usize,
    lo: usize,
    hi: usize,
) -> (f64, f64) {
    let m = stress.len();
    let max_lag = max_lag.min(m.saturating_sub(1));
    let mut sum = vec![0.0; max_lag + 1];
    let mut count = vec![0usize; max_lag + 1];
    for &c in &[3usize, 4, 5] {
        let mut a = vec![0.0; m];
        for t in 1..m {
            a[t] = a[t - 1] + dt * stress[t][c];
        }
        for t in 0..m {
            for k in 0..=max_lag {
                if t + k >= m {
                    break;
                }
                let d = a[t + k] - a[t];
                sum[k] += d * d;
                count[k] += 1;
            }
        }
    }
    for (s, c) in sum.iter_mut().zip(&count) {
        if *c > 0 {
            *s /= *c as f64;
        }
    }
    let hi = hi.min(max_lag);
    let x: Vec<f64> = (lo..=hi).map(|i| i as f64 * dt).collect();
    let y: Vec<f64> = sum[lo..=hi].to_vec();
    let (slope, se) = linfit(&x, &y);
    (volume * slope / (2.0 * kbt), volume * se / (2.0 * kbt))
}

/// Thermal conductivity from the three heat-flux components.
///
/// `heat[i]` holds the volume-integrated heat flux `J V` at sample `i`, so
/// the coefficient is `int/(V T^2)`.
pub fn green_kubo_conductivity(
    heat: &[[f64; 3]],
    dt: f64,
    volume: f64,
    kbt: f64,
    max_lag: usize,
    lo: usize,
    hi: usize,
) -> (f64, f64) {
    let series: Vec<Vec<f64>> = (0..3)
        .map(|c| heat.iter().map(|h| h[c]).collect())
        .collect();
    let run = average_running(&series, dt, max_lag);
    let scaled: Vec<f64> = run.iter().map(|i| i / (volume * kbt * kbt)).collect();
    plateau(&scaled, lo, hi)
}

/// Mean-squared displacement from unwrapped trajectory frames.
///
/// `frames` is flattened as `frame * n + atom`.
pub fn msd(frames: &[[f64; 3]], n: usize, max_lag: usize) -> Vec<f64> {
    let m = frames.len() / n;
    let max_lag = max_lag.min(m.saturating_sub(1));
    let mut out = vec![0.0; max_lag + 1];
    let mut count = vec![0usize; max_lag + 1];
    for f in 0..m {
        for k in 0..=max_lag {
            let g = f + k;
            if g >= m {
                break;
            }
            let mut acc = 0.0;
            for i in 0..n {
                let a = frames[f * n + i];
                let b = frames[g * n + i];
                acc += (b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2) + (b[2] - a[2]).powi(2);
            }
            out[k] += acc / n as f64;
            count[k] += 1;
        }
    }
    for (o, c) in out.iter_mut().zip(&count) {
        if *c > 0 {
            *o /= *c as f64;
        }
    }
    out
}

/// Velocity autocorrelation from trajectory frames.
pub fn vacf(frames: &[[f64; 3]], n: usize, max_lag: usize) -> Vec<f64> {
    let m = frames.len() / n;
    let max_lag = max_lag.min(m.saturating_sub(1));
    let mut out = vec![0.0; max_lag + 1];
    let mut count = vec![0usize; max_lag + 1];
    for f in 0..m {
        for k in 0..=max_lag {
            let g = f + k;
            if g >= m {
                break;
            }
            let mut acc = 0.0;
            for i in 0..n {
                let a = frames[f * n + i];
                let b = frames[g * n + i];
                acc += a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
            }
            out[k] += acc / n as f64;
            count[k] += 1;
        }
    }
    for (o, c) in out.iter_mut().zip(&count) {
        if *c > 0 {
            *o /= *c as f64;
        }
    }
    out
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

/// Self-diffusion from the MSD over the window `[lo, hi]` (indices into the
/// MSD array). Returns `(D, sigma_D)`.
pub fn diffusion_from_msd(msd: &[f64], dt: f64, lo: usize, hi: usize) -> (f64, f64) {
    let hi = hi.min(msd.len() - 1);
    let x: Vec<f64> = (lo..=hi).map(|i| i as f64 * dt).collect();
    let y: Vec<f64> = msd[lo..=hi].to_vec();
    let (slope, se) = linfit(&x, &y);
    (slope / 6.0, se / 6.0)
}

/// Self-diffusion from MSD block averaging: split the trajectory into
/// `nblocks` contiguous pieces, fit each, and return the mean and standard
/// error. This is the honest statistical error of the MSD estimator.
pub fn diffusion_from_msd_blocks(
    frames: &[[f64; 3]],
    n: usize,
    dt: f64,
    lo: usize,
    hi: usize,
    nblocks: usize,
) -> (f64, f64) {
    let m = frames.len() / n;
    if nblocks < 2 || m < 4 * nblocks {
        let full = msd(frames, n, hi);
        return diffusion_from_msd(&full, dt, lo, hi);
    }
    let blk = m / nblocks;
    let mut vals = Vec::with_capacity(nblocks);
    for b in 0..nblocks {
        let s = &frames[b * blk * n..(b + 1) * blk * n];
        let ms = msd(s, n, hi.min(blk - 1));
        let (d, _) = diffusion_from_msd(&ms, dt, lo, hi.min(blk - 1));
        vals.push(d);
    }
    let mean = vals.iter().sum::<f64>() / nblocks as f64;
    let var = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (nblocks - 1) as f64;
    (mean, (var / nblocks as f64).sqrt())
}

/// Self-diffusion from the VACF, `D = (1/3) int C_v dt`, integrated to
/// `hi`.
pub fn diffusion_from_vacf(vacf: &[f64], dt: f64, hi: usize) -> (f64, f64) {
    let hi = hi.min(vacf.len() - 1);
    let mut integ = 0.0;
    for k in 1..=hi {
        integ += 0.5 * dt * (vacf[k - 1] + vacf[k]);
    }
    let tail: f64 = (hi + 1..vacf.len())
        .map(|k| dt * vacf[k].abs())
        .sum::<f64>()
        .min(dt * vacf[hi].abs() * 10.0);
    ((1.0 / 3.0) * integ, (1.0 / 3.0) * tail)
}

/// Yeh & Hummer (2004) constant for the cubic periodic-box finite-size
/// correction, `xi = 2.837297`.
pub const YEH_HUMMER_XI: f64 = 2.837_297_477_735_76;

/// Subtract the finite-size enhancement from a periodic-box `D`.
/// `D_0 = D_PBC - T xi / (6 pi eta L)`, reduced units.
pub fn finite_size_corrected_d(d_pbc: f64, temperature: f64, eta: f64, box_l: f64) -> f64 {
    d_pbc - temperature * YEH_HUMMER_XI / (6.0 * std::f64::consts::PI * eta * box_l)
}

/// Integrated autocorrelation time `tau = int_0^inf C(t)/C(0) dt`, the
/// time scale a continuum closure must resolve.
pub fn correlation_time(corr: &Correlation) -> f64 {
    if corr.c.is_empty() || corr.c[0] == 0.0 {
        return 0.0;
    }
    let mut tau = 0.0;
    for k in 1..corr.c.len() {
        tau += 0.5 * corr.dt * (corr.c[k - 1] + corr.c[k]) / corr.c[0];
    }
    tau
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn autocorrelation_of_constant_is_zero() {
        let s = vec![3.0; 100];
        let c = autocorrelate(&s, 10);
        assert!(c.c.iter().all(|v| v.abs() < 1e-12));
    }

    #[test]
    fn msd_and_vacf_of_linear_motion_are_exact() {
        let n = 2;
        let v = 0.5;
        let dtf = 0.1;
        let m = 50;
        let mut pos = Vec::new();
        let mut vel = Vec::new();
        for f in 0..m {
            for i in 0..n {
                let sign = if i == 0 { 1.0 } else { -1.0 };
                pos.push([sign * v * f as f64 * dtf, 0.0, 0.0]);
                vel.push([sign * v, 0.0, 0.0]);
            }
        }
        let ms = msd(&pos, n, 20);
        for k in 0..=20 {
            let t = k as f64 * dtf;
            assert!((ms[k] - (v * t).powi(2)).abs() < 1e-10, "k {k}");
        }
        let vc = vacf(&vel, n, 20);
        for k in 0..=20 {
            assert!((vc[k] - v * v).abs() < 1e-10);
        }
    }

    #[test]
    fn finite_size_correction_lowers_d() {
        let d0 = finite_size_corrected_d(0.03, 0.722, 3.2, 8.5);
        assert!(d0 < 0.03);
        assert!(
            (0.03 - d0 - 0.722 * YEH_HUMMER_XI / (6.0 * std::f64::consts::PI * 3.2 * 8.5)).abs()
                < 1e-12
        );
    }
}
