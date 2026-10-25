//! Convection in the channels of a plate-fin heat sink.
//!
//! Two exact limits bracket the flow between two fins:
//!
//! * near the leading edge, each fin is an isolated flat plate, and the
//!   Blasius–Pohlhausen similarity solution ([`similarity`]) is exact for the
//!   boundary-layer equations;
//! * far downstream, the flow is fully developed, with `Nu_Dh = 7.5407` for
//!   uniform wall temperature (Shah & London 1978, *Laminar Flow Forced
//!   Convection in Ducts*, Table 104).
//!
//! A real sink usually sits between them. [`march`] solves the parabolized
//! (boundary-layer) equations for simultaneously developing flow between
//! parallel plates (the method of Bodoia & Osterle 1961, *Appl. Sci. Res. A*
//! 10:265): implicit in the streamwise direction, with the pressure gradient
//! fixed at each station by mass conservation. It recovers both limits (see
//! the tests), so it is the rung used for every sink, and its own
//! discretisation error is estimated by repeating it on a grid twice as fine.
//!
//! Scaling: half-gap `a`, mean velocity `U`, `X = x ν / (U a²)`, `Y = y / a`
//! (wall at 0, mid-plane at 1), `u/U`, `V = v a / ν`, `P = p / (ρ U²)`,
//! `θ = (T − T_wall) / (T_inlet − T_wall)`. Then `x / (D_h Re_Dh) = X / 16`.

/// Mean velocity times the half-gap, in the scaled variables, is 1.
#[derive(Clone, Copy, Debug)]
pub struct Station {
    pub x: f64,
    /// Local Nusselt number on the hydraulic diameter `D_h = 4a`.
    pub nu_local: f64,
    /// Bulk (mixed-mean) temperature.
    pub theta_b: f64,
    /// Pressure relative to the inlet.
    pub p: f64,
    /// Mid-plane velocity.
    pub u_center: f64,
}

impl Station {
    /// Mean Nusselt number over `[0, x]`, defined by the log-mean temperature
    /// difference so that `θ_b = exp(−NTU)` exactly.
    pub fn nu_mean(&self, pr: f64) -> f64 {
        -4.0 * pr * self.theta_b.ln() / self.x
    }
}

fn thomas(a: &[f64], b: &[f64], c: &[f64], r: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut cp = vec![0.0; n];
    let mut rp = vec![0.0; n];
    cp[0] = c[0] / b[0];
    rp[0] = r[0] / b[0];
    for i in 1..n {
        let m = b[i] - a[i] * cp[i - 1];
        cp[i] = c[i] / m;
        rp[i] = (r[i] - a[i] * rp[i - 1]) / m;
    }
    let mut x = vec![0.0; n];
    x[n - 1] = rp[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = rp[i] - cp[i] * x[i + 1];
    }
    x
}

/// March from the inlet to `x_end` (scaled) with `ny` cells across the half
/// gap (clustered at the wall) and `nx` streamwise steps (clustered at the
/// inlet). Returns every station.
pub fn march(pr: f64, x_end: f64, ny: usize, nx: usize) -> Vec<Station> {
    assert!(ny >= 4 && nx >= 2 && x_end > 0.0 && pr > 0.0);
    let y: Vec<f64> = (0..=ny).map(|j| 1.0 - (0.5 * std::f64::consts::PI * j as f64 / ny as f64).cos()).collect();
    let n = ny; // unknowns at j = 1..=ny
    // Finite-difference weights at interior nodes.
    let mut w = vec![[0.0f64; 3]; ny + 1]; // first derivative
    let mut d = vec![[0.0f64; 3]; ny + 1]; // second derivative
    for j in 1..ny {
        let (hm, hp) = (y[j] - y[j - 1], y[j + 1] - y[j]);
        w[j] = [-hp / (hm * (hm + hp)), (hp - hm) / (hm * hp), hm / (hp * (hm + hp))];
        let (dm, dp) = (2.0 / ((hm + hp) * hm), 2.0 / ((hm + hp) * hp));
        d[j] = [dm, -(dm + dp), dp];
    }
    let hn = y[ny] - y[ny - 1];
    d[ny] = [2.0 / (hn * hn), -2.0 / (hn * hn), 0.0]; // mirror ghost at the mid-plane
    let trapz = |f: &dyn Fn(usize) -> f64| -> f64 { (1..=ny).map(|j| 0.5 * (f(j) + f(j - 1)) * (y[j] - y[j - 1])).sum() };

    let mut u = vec![1.0; ny + 1];
    u[0] = 0.0;
    let mut th = u.clone();
    let flow = trapz(&|j| u[j]);
    let mut v = vec![0.0; ny + 1];
    let mut p = 0.0;
    let wall_grad = |f: &[f64]| -> f64 {
        let (h1, h2) = (y[1], y[2]);
        f[1] * h2 / (h1 * (h2 - h1)) - f[2] * h1 / (h2 * (h2 - h1))
    };
    let mut out = Vec::with_capacity(nx);
    let mut x_prev = 0.0;
    for k in 1..=nx {
        let x = x_end * (k as f64 / nx as f64).powi(2);
        let dx = x - x_prev;
        let (u_old, th_old) = (u.clone(), th.clone());
        let (mut ub, mut vb) = (u_old.clone(), v.clone());
        let mut g = 0.0;
        for _ in 0..3 {
            // Momentum: A u = u_old ū/Δ − G, solved as u = u_a − G u_b.
            let (mut a, mut b, mut c) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
            let mut r = vec![0.0; n];
            for j in 1..=ny {
                let i = j - 1;
                a[i] = vb[j] * w[j][0] - d[j][0];
                b[i] = ub[j] / dx + vb[j] * w[j][1] - d[j][1];
                c[i] = vb[j] * w[j][2] - d[j][2];
                r[i] = ub[j] * u_old[j] / dx;
            }
            a[0] = 0.0;
            let ua = thomas(&a, &b, &c, &r);
            let ubb = thomas(&a, &b, &c, &vec![1.0; n]);
            let q = |s: &[f64]| -> f64 { (1..=ny).map(|j| 0.5 * (s[j - 1] + if j >= 2 { s[j - 2] } else { 0.0 }) * (y[j] - y[j - 1])).sum() };
            g = (q(&ua) - flow) / q(&ubb);
            for j in 1..=ny {
                u[j] = ua[j - 1] - g * ubb[j - 1];
            }
            // Continuity: V_Y = −(u − u_old)/Δ, V(0) = 0.
            v[0] = 0.0;
            for j in 1..=ny {
                let s = 0.5 * ((u[j] - u_old[j]) + (u[j - 1] - u_old[j - 1])) / dx;
                v[j] = v[j - 1] - s * (y[j] - y[j - 1]);
            }
            ub = u.clone();
            vb = v.clone();
        }
        p += g * dx;
        // Energy with the converged velocity.
        let (mut a, mut b, mut c) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        let mut r = vec![0.0; n];
        for j in 1..=ny {
            let i = j - 1;
            a[i] = v[j] * w[j][0] - d[j][0] / pr;
            b[i] = u[j] / dx + v[j] * w[j][1] - d[j][1] / pr;
            c[i] = v[j] * w[j][2] - d[j][2] / pr;
            r[i] = u[j] * th_old[j] / dx;
        }
        a[0] = 0.0;
        let t = thomas(&a, &b, &c, &r);
        for j in 1..=ny {
            th[j] = t[j - 1];
        }
        let theta_b = trapz(&|j| u[j] * th[j]) / trapz(&|j| u[j]);
        out.push(Station { x, nu_local: 4.0 * wall_grad(&th) / theta_b, theta_b, p, u_center: u[ny] });
        x_prev = x;
    }
    out
}

/// The developing-channel result used by a heat sink, with its estimated
/// discretisation error.
#[derive(Clone, Copy, Debug)]
pub struct Channel {
    pub nu_mean: f64,
    /// Estimated absolute error of `nu_mean`: the change from a grid half as
    /// fine in both directions (conservative while the error is asymptotic).
    pub nu_error: f64,
    /// `(p_inlet − p_outlet) / (ρ U²)`: core friction only.
    pub dp: f64,
    pub dp_error: f64,
}

pub fn channel(pr: f64, x_end: f64) -> Channel {
    let coarse = *march(pr, x_end, 40, 200).last().unwrap();
    let fine = *march(pr, x_end, 80, 400).last().unwrap();
    Channel {
        nu_mean: fine.nu_mean(pr),
        nu_error: (fine.nu_mean(pr) - coarse.nu_mean(pr)).abs(),
        dp: -fine.p,
        dp_error: (fine.p - coarse.p).abs(),
    }
}

/// The Blasius–Pohlhausen flat-plate boundary layer, the exact leading-edge
/// limit. Writes the similarity equations `f''' + f f''/2 = 0` and
/// `θ'' + Pr f θ'/2 = 0` with `η = y √(U/(ν x))`.
pub mod similarity {
    fn rk4(s: f64, eta_max: f64, steps: usize) -> Vec<[f64; 4]> {
        // state: f, f', f'', F = ∫f
        let rhs = |z: [f64; 4]| [z[1], z[2], -0.5 * z[0] * z[2], z[0]];
        let h = eta_max / steps as f64;
        let mut z = [0.0, 0.0, s, 0.0];
        let mut out = vec![z];
        for _ in 0..steps {
            let k1 = rhs(z);
            let k2 = rhs(std::array::from_fn(|i| z[i] + 0.5 * h * k1[i]));
            let k3 = rhs(std::array::from_fn(|i| z[i] + 0.5 * h * k2[i]));
            let k4 = rhs(std::array::from_fn(|i| z[i] + h * k3[i]));
            z = std::array::from_fn(|i| z[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]));
            out.push(z);
        }
        out
    }

    const ETA: f64 = 12.0;
    const STEPS: usize = 4800;

    /// Wall shear `f''(0)` of the Blasius solution, by shooting.
    pub fn blasius() -> f64 {
        let (mut lo, mut hi) = (0.1, 1.0);
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if rk4(mid, ETA, STEPS).last().unwrap()[1] < 1.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// `θ'(0)`, so that the local Nusselt number is `θ'(0) Re_x^{1/2}`.
    /// From `θ' = θ'(0) exp(−Pr F/2)`, `F = ∫ f`, and `θ(∞) = 1`.
    pub fn theta_wall(pr: f64) -> f64 {
        let z = rk4(blasius(), ETA, STEPS);
        let h = ETA / STEPS as f64;
        // Simpson over η; beyond ETA the integrand is below exp(−Pr·60/2·…).
        let g: Vec<f64> = z.iter().map(|s| (-0.5 * pr * s[3]).exp()).collect();
        let mut sum = g[0] + g[STEPS];
        for (i, v) in g.iter().enumerate().take(STEPS).skip(1) {
            sum += if i % 2 == 1 { 4.0 * v } else { 2.0 * v };
        }
        1.0 / (sum * h / 3.0)
    }
}

/// Spreading resistance of a centred `c × d` isoflux source on an `a × b`
/// plate of thickness `t` and conductivity `k`, with edges adiabatic and a
/// uniform heat transfer coefficient `h` on the opposite face: the mean
/// source temperature minus the 1-D value `t/(kab) + 1/(hab)`, per watt.
///
/// Exact cosine-series solution of Laplace's equation. Mode `(m, n)` with
/// `β = π √((m/a)² + (n/b)²)` responds at the heated face with
/// `G = (1 + ε tanh βt) / (kβ (tanh βt + ε))`, `ε = h/(kβ)`, and the
/// source's coefficients `s_m = 4 (−1)^{m/2} sin(mπc/2a) / (mπ)` (even `m`
/// only, by symmetry). Returns the resistance and the change from summing
/// half as many modes (the truncation estimate).
pub fn spreading(a: f64, b: f64, c: f64, d: f64, t: f64, k: f64, h: f64) -> (f64, f64) {
    assert!(c > 0.0 && d > 0.0 && c <= a && d <= b && t > 0.0 && k > 0.0 && h > 0.0);
    // Enough modes to resolve the source edges many times over.
    let modes = |len: f64, src: f64| -> usize { (40.0 * len / src).ceil() as usize + 40 };
    let (nx, ny) = (modes(a, c), modes(b, d));
    let weight = |m: usize, len: f64, src: f64| -> f64 {
        if m == 0 {
            src / len
        } else {
            let mf = m as f64;
            let s = 4.0 * (0.5 * mf * std::f64::consts::PI * src / len).sin() / (mf * std::f64::consts::PI);
            len / (2.0 * src) * s * s
        }
    };
    let wx: Vec<f64> = (0..=nx).map(|i| weight(2 * i, a, c)).collect();
    let wy: Vec<f64> = (0..=ny).map(|j| weight(2 * j, b, d)).collect();
    let sum = |kx: usize, ky: usize| -> f64 {
        let mut s = 0.0;
        for i in 0..=kx {
            if wx[i] == 0.0 {
                continue;
            }
            let lx = 2.0 * i as f64 * std::f64::consts::PI / a;
            for j in 0..=ky {
                if (i == 0 && j == 0) || wy[j] == 0.0 {
                    continue;
                }
                let ly = 2.0 * j as f64 * std::f64::consts::PI / b;
                let beta = (lx * lx + ly * ly).sqrt();
                let eps = h / (k * beta);
                let tau = (beta * t).tanh();
                s += wx[i] * wy[j] * (1.0 + eps * tau) / (k * beta * (tau + eps));
            }
        }
        s / (c * d)
    };
    let full = sum(nx, ny);
    let half = sum(nx / 2, ny / 2);
    (full, (full - half).abs())
}
