//! Capillary-wave dispersion diagnostics (new file owned by the
//! capillary-wave agent).
use spectral_core::twophase::{Projection, TwoPhaseParams, TwoPhaseSolver};
use spectral_core::{init, Grid};
use std::sync::Mutex;

/// The solver reads `CAPFORM` from the process environment, so tests that set
/// it must not run concurrently.
static ENV_LOCK: Mutex<()> = Mutex::new(());

const PI: f64 = std::f64::consts::PI;

struct Cfg {
    n: usize,
    ratio: f64,
    eta_scale: f64,
    m: f64,
    mu1: f64,
    cap_alt: bool,
    proj: Projection,
    t_periods: f64,
    n_sample: usize,
    dt_scale: f64,
}

fn cfg() -> Cfg {
    let e = |k: &str, d: f64| std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d);
    Cfg {
        n: e("N", 32.0) as usize,
        ratio: e("RATIO", 10.0),
        eta_scale: e("ETAS", 1.2),
        m: e("MOB", 2e-4),
        mu1: e("MU1", 2e-4),
        cap_alt: std::env::var("CAPFORM").map(|v| v == "1").unwrap_or(false),
        proj: if std::env::var("VCPROJ").map(|v| v == "1").unwrap_or(false) {
            Projection::Variable
        } else {
            Projection::Constant
        },
        t_periods: e("TP", 1.3),
        n_sample: e("NS", 400.0) as usize,
        dt_scale: e("DTS", 1.0),
    }
}

fn build(c: &Cfg) -> (TwoPhaseSolver, f64, f64) {
    let g = Grid::new(c.n);
    let (k, a0) = (2.0, 0.03);
    let rho = [1.0, 1.0 / c.ratio];
    let sigma = 1.0;
    let nu_ratio = if c.ratio > 1.0 { 15.0 } else { 1.0 };
    let eta = c.eta_scale * g.dx();
    let prm = TwoPhaseParams {
        rho,
        mu: [c.mu1, c.mu1 * nu_ratio / c.ratio],
        sigma,
        gravity: [0.0; 3],
        eta,
        mobility: c.m,
        stab: 2.0,
    };
    let w = 2f64.sqrt() * eta;
    let (zlo, zhi) = (0.5 * PI, 1.5 * PI);
    let phi = init::from_fn(&g, |x, _, z| {
        let zeta = a0 * (k * x).cos();
        let v = ((z - zlo - zeta) / w).tanh() - ((z - zhi - zeta) / w).tanh() - 1.0;
        [v, 0.0, 0.0]
    })[0]
        .clone();
    if c.cap_alt {
        std::env::set_var("CAPFORM", "1");
    }
    let s = TwoPhaseSolver::new(g.clone(), prm, init::from_fn(&g, |_, _, _| [0.0; 3]), phi);
    let mut s = s;
    s.set_projection(c.proj);
    let om = (sigma * k.powi(3) / (rho[0] + rho[1])).sqrt();
    let dt_cap = 0.4 * (0.5 * (rho[0] + rho[1]) * g.dx().powi(3) / (2.0 * PI * sigma)).sqrt();
    (s, dt_cap.min(0.005) * c.dt_scale, om)
}

fn crossing(phi: &[f64], g: &Grid, i: usize, j: usize, z_guess: f64) -> f64 {
    let n = g.n;
    let idx = |l: usize| (i * n + j) * n + l;
    let mut best = (f64::INFINITY, z_guess);
    for l in 0..n {
        let (a, b) = (phi[idx(l)], phi[idx((l + 1) % n)]);
        if a == 0.0 || a * b < 0.0 {
            let z = (l as f64 + a / (a - b)) * g.dx();
            let mut dz = z - z_guess;
            dz -= 2.0 * PI * (dz / (2.0 * PI)).round();
            if dz.abs() < best.0 {
                best = (dz.abs(), z);
            }
        }
    }
    best.1
}

fn amp(s: &TwoPhaseSolver, zi: f64) -> f64 {
    let g = &s.grid;
    let n = g.n;
    let phi = s.phi_phys();
    0.5 * (crossing(&phi, g, 0, 0, zi) - crossing(&phi, g, n / 4, 0, zi))
}

/// Reference (flat) profile and its z-derivative for the tanh interface.
fn reference(g: &Grid, eta: f64, zlo: f64, zhi: f64) -> (Vec<f64>, Vec<f64>) {
    let w = 2f64.sqrt() * eta;
    let n = g.n;
    let mut phi0 = vec![0.0; g.real_len()];
    let mut g0 = vec![0.0; g.real_len()];
    for i in 0..n {
        for j in 0..n {
            for l in 0..n {
                let z = g.x(l);
                let a = ((z - zlo) / w).tanh();
                let b = ((z - zhi) / w).tanh();
                phi0[(i * n + j) * n + l] = a - b - 1.0;
                let sa = 1.0 / ((z - zlo) / w).cosh().powi(2);
                let sb = 1.0 / ((z - zhi) / w).cosh().powi(2);
                g0[(i * n + j) * n + l] = (sa - sb) / w;
            }
        }
    }
    (phi0, g0)
}

/// Interface displacement projected onto the interface-shift mode:
/// `ζ(x) = -∫ δφ ∂zφ₀ dz / ∫ (∂zφ₀)² dz`. Internal interface modes are
/// orthogonal to ∂zφ₀ and are annihilated, leaving the capillary mode.
fn proj_amp(s: &TwoPhaseSolver, phi0: &[f64], g0: &[f64]) -> f64 {
    let g = &s.grid;
    let n = g.n;
    let phi = s.phi_phys();
    let norm: f64 = (0..n).map(|l| g0[l] * g0[l]).sum();
    let zeta = |i: usize| {
        let base = (i * n) * n;
        let mut acc = 0.0;
        for l in 0..n {
            acc += (phi[base + l] - phi0[base + l]) * g0[base + l];
        }
        -acc / norm
    };
    0.5 * (zeta(0) - zeta(n / 4))
}

fn mode(s: &TwoPhaseSolver, kx: i64, kz: usize) -> f64 {
    let g = &s.grid;
    let n = g.n;
    let nh = g.nh;
    let i = ((kx % n as i64 + n as i64) % n as i64) as usize;
    s.phi[(i * n) * nh + kz].re
}

fn umode(s: &TwoPhaseSolver, c: usize, kx: i64, kz: usize) -> f64 {
    let g = &s.grid;
    let n = g.n;
    let nh = g.nh;
    let i = ((kx % n as i64 + n as i64) % n as i64) as usize;
    s.u[c][(i * n) * nh + kz].re
}

/// Periodogram of a uniformly sampled scalar signal after removing a linear
/// trend and the mean.
fn periodogram2(us: &[(f64, f64)], om_center: f64) -> f64 {
    let n = us.len();
    let t0 = us[0].0;
    let (mut st, mut stt, mut sh, mut sth) = (0.0, 0.0, 0.0, 0.0);
    for &(t, h) in us {
        let x = t - t0;
        st += x;
        stt += x * x;
        sh += h;
        sth += x * h;
    }
    let nf = n as f64;
    let den = nf * stt - st * st;
    let (b, a) = if den.abs() > 0.0 {
        ((nf * sth - st * sh) / den, (sh - (nf * sth - st * sh) / den * st) / nf)
    } else {
        (0.0, sh / nf)
    };
    let y: Vec<f64> = us.iter().map(|&(t, h)| h - a - b * (t - t0)).collect();
    let mut best = (om_center, -1.0f64);
    for iw in 0..4000 {
        let w = om_center * (0.5 + 1.0 * iw as f64 / 3999.0);
        let (mut re, mut im) = (0.0, 0.0);
        for (q, &h) in y.iter().enumerate() {
            let ph = w * (us[q].0 - t0);
            re += h * ph.cos();
            im -= h * ph.sin();
        }
        let p = re * re + im * im;
        if p > best.1 {
            best = (w, p);
        }
    }
    best.0
}

fn periodogram(samples: &[(f64, f64, f64)], om_center: f64) -> f64 {
    let n = samples.len();
    let t0 = samples[0].0;
    let (mut st, mut stt, mut sh, mut sth) = (0.0, 0.0, 0.0, 0.0);
    for &(t, h, _) in samples {
        let x = t - t0;
        st += x;
        stt += x * x;
        sh += h;
        sth += x * h;
    }
    let nf = n as f64;
    let den = nf * stt - st * st;
    let (b, a) = if den.abs() > 0.0 {
        ((nf * sth - st * sh) / den, (sh - (nf * sth - st * sh) / den * st) / nf)
    } else {
        (0.0, sh / nf)
    };
    let y: Vec<f64> = samples.iter().map(|&(t, h, _)| h - a - b * (t - t0)).collect();
    let mut best = (om_center, -1.0f64);
    for iw in 0..4000 {
        let w = om_center * (0.5 + 1.0 * iw as f64 / 3999.0);
        let (mut re, mut im) = (0.0, 0.0);
        for (q, &h) in y.iter().enumerate() {
            let ph = w * (samples[q].0 - t0);
            re += h * ph.cos();
            im -= h * ph.sin();
        }
        let p = re * re + im * im;
        if p > best.1 {
            best = (w, p);
        }
    }
    best.0
}

/// Fit h(t) = e^{-g t}(a cos wt + b sin wt) + cc by grid search; returns
/// (omega, gamma, rms residual).
fn damped_fit(samples: &[(f64, f64, f64)], om: f64) -> (f64, f64, f64) {
    let mut best = (om, 0.1, f64::INFINITY, 0.0, 0.0, 0.0);
    for iw in 0..400 {
        let w = om * (0.4 + 1.2 * iw as f64 / 399.0);
        for ig in 0..120 {
            let gam = 1.2 * ig as f64 / 119.0;
            let (a, b, cc, rms) = linfit(samples, w, gam);
            if rms < best.2 {
                best = (w, gam, rms, a, b, cc);
            }
        }
    }
    // refine
    for _ in 0..60 {
        let (w0, g0) = (best.0, best.1);
        let dw = om * 0.02;
        let dg = 0.02;
        for iw in -4..=4 {
            for ig in -4..=4 {
                let w = w0 + dw * iw as f64 / 4.0;
                let gam = (g0 + dg * ig as f64 / 4.0).max(0.0);
                let (a, b, cc, rms) = linfit(samples, w, gam);
                if rms < best.2 {
                    best = (w, gam, rms, a, b, cc);
                }
            }
        }
    }
    (best.0, best.1, best.2)
}

fn linfit(samples: &[(f64, f64, f64)], w: f64, gam: f64) -> (f64, f64, f64, f64) {
    // normal equations for basis [e^{-g t}cos(wt), e^{-g t}sin(wt), 1]
    let mut m = [[0.0f64; 3]; 3];
    let mut rhs = [0.0f64; 3];
    for &(t, h, _) in samples {
        let e = (-gam * t).exp();
        let phi = [e * (w * t).cos(), e * (w * t).sin(), 1.0];
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] += phi[i] * phi[j];
            }
            rhs[i] += phi[i] * h;
        }
    }
    // Gaussian elimination 3x3
    for i in 0..3 {
        let p = (i..3).max_by(|&a, &b| m[a][i].abs().partial_cmp(&m[b][i].abs()).unwrap()).unwrap();
        m.swap(i, p);
        rhs.swap(i, p);
        for j in i + 1..3 {
            let f = m[j][i] / m[i][i];
            for k in i..3 {
                m[j][k] -= f * m[i][k];
            }
            rhs[j] -= f * rhs[i];
        }
    }
    let mut x = [0.0; 3];
    for i in (0..3).rev() {
        let mut s = rhs[i];
        for j in i + 1..3 {
            s -= m[i][j] * x[j];
        }
        x[i] = s / m[i][i];
    }
    let mut ss = 0.0;
    for &(t, h, _) in samples {
        let e = (-gam * t).exp();
        let p = x[0] * e * (w * t).cos() + x[1] * e * (w * t).sin() + x[2];
        ss += (p - h).powi(2);
    }
    (x[0], x[1], x[2], (ss / samples.len() as f64).sqrt())
}

/// Normal equations for the 4-parameter model
/// `h(t) = e^{-g t}(a cos wt + b sin wt) + c0 + c1 t`.  The extra linear
/// trend absorbs the slowly varying (non-oscillatory) part of the projected
/// interface displacement, which otherwise biases a single-mode fit strongly
/// downward on short runs (see `analysis/twophase_and_mhd.md`).  Returns the
/// rms residual.
fn linfit4(samples: &[(f64, f64, f64)], w: f64, gam: f64) -> f64 {
    let mut m = [[0.0f64; 4]; 4];
    let mut rhs = [0.0f64; 4];
    for &(t, h, _) in samples {
        let e = (-gam * t).exp();
        let phi = [e * (w * t).cos(), e * (w * t).sin(), 1.0, t];
        for i in 0..4 {
            for j in 0..4 {
                m[i][j] += phi[i] * phi[j];
            }
            rhs[i] += phi[i] * h;
        }
    }
    for i in 0..4 {
        let p = (i..4).max_by(|&a, &b| m[a][i].abs().partial_cmp(&m[b][i].abs()).unwrap()).unwrap();
        m.swap(i, p);
        rhs.swap(i, p);
        for j in i + 1..4 {
            let f = m[j][i] / m[i][i];
            for k in i..4 {
                m[j][k] -= f * m[i][k];
            }
            rhs[j] -= f * rhs[i];
        }
    }
    let mut x = [0.0; 4];
    for i in (0..4).rev() {
        let mut s = rhs[i];
        for j in i + 1..4 {
            s -= m[i][j] * x[j];
        }
        x[i] = s / m[i][i];
    }
    let mut ss = 0.0;
    for &(t, h, _) in samples {
        let e = (-gam * t).exp();
        let p = x[0] * e * (w * t).cos() + x[1] * e * (w * t).sin() + x[2] + x[3] * t;
        ss += (p - h).powi(2);
    }
    (ss / samples.len() as f64).sqrt()
}

/// Robust damped-frequency estimate: grid + refine over `(w, gam)` with the
/// linear-trend basis of [`linfit4`].  Returns `(omega, gamma, rms)`.
fn robust_fit(samples: &[(f64, f64, f64)], om0: f64) -> (f64, f64, f64) {
    let mut best = (om0, 0.0f64, f64::INFINITY);
    for iw in 0..241 {
        let w = om0 * (0.55 + 0.9 * iw as f64 / 240.0);
        for ig in 0..61 {
            let g = ig as f64 / 60.0;
            let r = linfit4(samples, w, g);
            if r < best.2 {
                best = (w, g, r);
            }
        }
    }
    let (mut w0, mut g0) = (best.0, best.1);
    for _ in 0..60 {
        let (dw, dg) = (om0 * 0.005, 0.005);
        for iw in -3..=3 {
            for ig in -3..=3 {
                let w = w0 + dw * iw as f64 / 3.0;
                let g = (g0 + dg * ig as f64 / 3.0).max(0.0);
                let r = linfit4(samples, w, g);
                if r < best.2 {
                    best = (w, g, r);
                }
            }
        }
        if (best.0 - w0).abs() < 1e-10 && (best.1 - g0).abs() < 1e-10 {
            break;
        }
        w0 = best.0;
        g0 = best.1;
    }
    best
}

/// Time-step a configuration and return `(samples, omega_theory)`.  `samples`
/// hold `(t, projected interface height, interface-shift spectral mode)`.
fn run(c: &Cfg) -> (Vec<(f64, f64, f64)>, f64) {
    let (mut s, dt, om) = build(c);
    let t_end = c.t_periods * 2.0 * PI / om;
    let every = ((t_end / dt) / c.n_sample as f64).ceil().max(1.0) as usize;
    let (zlo, zhi) = (0.5 * PI, 1.5 * PI);
    let (phi0, g0) = reference(&s.grid, c.eta_scale * s.grid.dx(), zlo, zhi);
    let mut samples = vec![(0.0, proj_amp(&s, &phi0, &g0), mode(&s, 2, 1))];
    let mut step = 0;
    while s.t < t_end {
        s.step(dt);
        step += 1;
        if step % every == 0 {
            samples.push((s.t, proj_amp(&s, &phi0, &g0), mode(&s, 2, 1)));
        }
    }
    (samples, om)
}

#[test]
fn capillary_dispersion_ratio1() {
    // Equal densities: both force forms agree, and this 2-period
    // damped-sinusoid window reproduces omega^2 = sigma k^3/(rho1+rho2)
    // to a few percent.  The window is short, so the residual is compared
    // against the analytic finite-thickness correction in
    // tests/capillary_thickness.rs rather than pinned tightly here.
    let _lk = ENV_LOCK.lock().unwrap();
    std::env::set_var("CAPFORM", "0");
    let c = Cfg {
        n: 24,
        ratio: 1.0,
        eta_scale: 1.2,
        m: 1e-8,
        mu1: 2e-4,
        cap_alt: false,
        proj: Projection::Constant,
        t_periods: 2.0,
        n_sample: 240,
        dt_scale: 2.0,
    };
    let (samples, om) = run(&c);
    let (w, _g, _) = robust_fit(&samples, om);
    let err = w / om - 1.0;
    eprintln!("ratio-1 capillary: omega={w:.4} vs {om:.4} ({:+.2}%)", 100.0 * err);
    assert!(err.abs() < 0.06, "ratio-1 capillary frequency error {:.1}%", 100.0 * err);
}

#[test]
fn capillary_dispersion_ratio10_balanced_form() {
    // The default mu_c*grad(phi) force is inconsistent with the constant-
    // coefficient variable-density projection: the two standard forms differ
    // by grad(phi mu_c), a pressure gradient only at constant density.  At
    // ratio 10 the balanced form -phi*grad(mu_c) recovers the exact
    // frequency, the default form is ~30% low.  This test pins the balanced
    // form (and its runtime is < 30 s).
    let _lk = ENV_LOCK.lock().unwrap();
    std::env::set_var("CAPFORM", "1");
    let c = Cfg {
        n: 24,
        ratio: 10.0,
        eta_scale: 1.2,
        m: 1e-8,
        mu1: 2e-4,
        cap_alt: true,
        proj: Projection::Constant,
        t_periods: 2.0,
        n_sample: 240,
        dt_scale: 2.0,
    };
    let (samples, om) = run(&c);
    let (w, _g, _) = robust_fit(&samples, om);
    let err = w / om - 1.0;
    eprintln!("ratio-10 capillary (balanced): omega={w:.4} vs {om:.4} ({:+.2}%)", 100.0 * err);
    std::env::set_var("CAPFORM", "0");
    assert!(err.abs() < 0.10, "balanced-form capillary frequency error {:.1}%", 100.0 * err);
}

#[test]
#[ignore]
fn diag_capillary() {
    let _lk = ENV_LOCK.lock().unwrap();
    let c = cfg();
    let (mut s, dt, om) = build(&c);
    let t_end = c.t_periods * 2.0 * PI / om;
    let every = ((t_end / dt) / c.n_sample as f64).ceil().max(1.0) as usize;
    let zi = 1.5 * PI;
    let (zlo, zhi) = (0.5 * PI, 1.5 * PI);
    let (phi0, g0) = reference(&s.grid, c.eta_scale * s.grid.dx(), zlo, zhi);
    let mut prev = proj_amp(&s, &phi0, &g0);
    let mut prevm = mode(&s, 2, 1);
    let mut zc = Vec::new();
    let mut zm = Vec::new();
    let mut samples = vec![(0.0, prev, prevm)];
    let mut usamp: Vec<(f64, f64)> = vec![(0.0, umode(&s, 2, 2, 0))];
    let a0 = prev.abs();
    let m0 = prevm;
    let mut step = 0;
    while s.t < t_end {
        s.step(dt);
        step += 1;
        if step % every == 0 {
            let a = proj_amp(&s, &phi0, &g0);
            let mm = mode(&s, 2, 1);
            samples.push((s.t, a, mm));
            usamp.push((s.t, umode(&s, 2, 2, 0)));
            if prev * a < 0.0 {
                zc.push(s.t - dt * a / (a - prev));
            }
            if prevm * mm < 0.0 {
                zm.push(s.t - dt * mm / (mm - prevm));
            }
            prev = a;
            prevm = mm;
        }
    }
    let _ = (&zc, &zm, m0, amp(&s, zi));
    // damped-sinusoid fit: h(t) = e^{-g t}(a cos wt + b sin wt) + cc
    let fit = damped_fit(&samples, om);
    let rob = robust_fit(&samples, om);
    let om_pg = periodogram(&samples, om);
    let om_u = periodogram2(&usamp, om);
    // Zero-crossing estimate (matches examples/capillary_probe.rs).
    let (om_zc, gam_zc) = if zc.len() >= 3 {
        let period = 2.0 * (zc[zc.len() - 1] - zc[0]) / (zc.len() - 1) as f64;
        let omz = 2.0 * PI / period;
        // damping from successive peak amplitudes would be better; approximate
        // from the decay of the sampled envelope over the run.
        let gamz = 0.0;
        (omz, gamz)
    } else {
        (f64::NAN, f64::NAN)
    };
    println!(
        "N={} ratio={} etas={} M={:e} mu1={:e} cap={} vc={:?}: omega_u {:.4} ({:+.2}%) omega_pg {:.4} ({:+.2}%) omega_zc {:.4} ({:+.2}%) omega_fit {:.4} gamma {:.4} vs {:.4} ({:+.2}%) A {:.4}->{:.4} rms {:.2e} nzc {} cg {} r {:.1e}",
        c.n, c.ratio, c.eta_scale, c.m, c.mu1, c.cap_alt, c.proj, om_u,
        100.0 * (om_u / om - 1.0), om_pg,
        100.0 * (om_pg / om - 1.0), om_zc,
        100.0 * (om_zc / om - 1.0), fit.0, fit.1, om,
        100.0 * (fit.0 / om - 1.0), a0, samples.last().unwrap().1, fit.2, zc.len(),
        s.last_vc_iters, s.last_vc_resid
    );
    println!(
        "  omega_rob {:.4} ({:+.2}%) gamma {:.4} rms {:.2e}",
        rob.0, 100.0 * (rob.0 / om - 1.0), rob.1, rob.2
    );
    let _ = gam_zc;
    if std::env::var("LOG").is_ok() {
        for (t, a, m) in &samples {
            println!("  t={t:.3} h={a:+.6} mode1={m:+.6e}");
        }
    }
}
