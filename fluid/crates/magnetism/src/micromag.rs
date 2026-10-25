//! Micromagnetics: the coarse rung below the atomistic spins.
//!
//! The continuum energy density is
//!
//! ```text
//! e = A (grad m)^2 - K (m . e)^2 - mu0 M_s H . m
//!     + (mu0/2) M_s^2 sum_k N_k (m . n_k)^2,
//! ```
//!
//! with the exchange stiffness `A` folded up from `J` in [`crate::geometry`],
//! uniaxial anisotropy `K`, saturation magnetisation `M_s` and local-shape
//! demagnetising factors `N_k`.  The Landau-Lifshitz-Gilbert equation
//!
//! ```text
//! dm/dt = -gamma/(1+alpha^2) [ m x H_eff + alpha m x (m x H_eff) ]
//! ```
//!
//! is integrated with RK4 on a 1D finite-difference grid.  The exact Bloch
//! wall `theta(x) = 2 atan(exp(x/delta))` has width `delta = sqrt(A/K)` and
//! energy per area `sigma = 4 sqrt(A K)`; the field-driven wall has the
//! Walker velocity `v = gamma delta H / alpha` below the Walker field
//! `H_W = alpha M_s / 2` (Schryer & Walker 1974; Malozemoff & Slonczewski
//! 1979; the 1D model result quoted in Bertotti, *Hysteresis in Magnetism*).
//! The thin-film Kittel resonance is
//! `omega = gamma sqrt((H+H_K)(H+H_K+M_s))`.

use crate::geometry::{dot, norm, V3};

/// `mu0 gamma_e`, the SI LLG gyromagnetic ratio, in m/(A s):
/// `mu0 * 1.760859644e11 rad/(s T) = 2.211e5 m/(A s)` (CODATA 2018).
pub const GAMMA: f64 = 2.211e5;
pub const MU0: f64 = 1.256_637_062_12e-6;

pub fn normalize(v: V3) -> V3 {
    let n = norm(v).max(1e-300);
    [v[0] / n, v[1] / n, v[2] / n]
}

/// Domain-wall width parameter `delta = sqrt(A/K)` (m).
pub fn wall_delta(a: f64, k: f64) -> f64 {
    (a / k).sqrt()
}

/// Domain-wall energy per unit area `sigma = 4 sqrt(A K)` (J/m^2).
pub fn wall_sigma(a: f64, k: f64) -> f64 {
    4.0 * (a * k).sqrt()
}

/// Anisotropy field `H_K = 2K/(mu0 M_s)` (A/m).
pub fn anisotropy_field(k: f64, ms: f64) -> f64 {
    2.0 * k / (MU0 * ms)
}

/// Exchange length `l_ex = sqrt(2A/(mu0 M_s^2))` (m).
pub fn exchange_length(a: f64, ms: f64) -> f64 {
    (2.0 * a / (MU0 * ms * ms)).sqrt()
}

/// Exact Bloch-wall profile: `m = (cos theta, sin theta, 0)` with
/// `theta = 2 atan(exp(x/delta))`, easy axis x, wall normal x.
pub fn wall_profile(x: f64, delta: f64) -> V3 {
    let theta = 2.0 * ((x / delta).exp()).atan();
    [theta.cos(), theta.sin(), 0.0]
}

/// Discrete atomistic 1D chain wall: exchange `J` between neighbours,
/// easy-axis anisotropy `K`, lattice constant `a`, unit spins.  The fold
/// maps to `A = J a/2`, `K_1d = K/a`.  Relaxes the interior angles by
/// nonlinear over-relaxation with fixed end spins `theta_0=0`,
/// `theta_{n-1}=pi`, and returns the measured width parameter and the
/// energy cost relative to the uniform state.
pub struct DiscreteWall {
    pub delta_measured: f64,
    pub sigma_measured: f64,
    pub delta_continuum: f64,
    pub sigma_continuum: f64,
    pub ratio_delta_over_a: f64,
}

pub fn discrete_wall_chain(j: f64, k: f64, a: f64, n: usize, iterations: usize) -> DiscreteWall {
    // Continuum coefficients for this 1D chain.
    let a_ex = 0.5 * j * a;
    let k_1d = k / a;
    let delta_c = wall_delta(a_ex, k_1d);
    let sigma_c = wall_sigma(a_ex, k_1d);

    // Relax theta by damped gradient descent (accept/reject line search) on
    // E = sum_a [-J cos(dtheta) - K cos^2 theta], starting from the
    // continuum profile.
    let mut theta = vec![0.0; n];
    let x0 = 0.5 * (n - 1) as f64;
    for (i, t) in theta.iter_mut().enumerate() {
        *t = 2.0 * (((i as f64 - x0) / delta_c).exp()).atan();
    }
    let energy = |t: &[f64]| -> f64 {
        let mut e = 0.0;
        for i in 0..n - 1 {
            e -= j * (t[i] - t[i + 1]).cos();
        }
        for &ti in t {
            e -= k * ti.cos().powi(2);
        }
        e
    };
    let mut e = energy(&theta);
    let mut step = 0.001;
    for _ in 0..iterations {
        let mut grad = vec![0.0; n];
        for i in 1..n - 1 {
            grad[i] = j * ((theta[i] - theta[i - 1]).sin() + (theta[i] - theta[i + 1]).sin())
                + k * (2.0 * theta[i]).sin();
        }
        let trial: Vec<f64> = theta
            .iter()
            .enumerate()
            .map(|(i, &t)| if i == 0 || i == n - 1 { t } else { t - step * grad[i] })
            .collect();
        let e_trial = energy(&trial);
        if e_trial < e {
            theta = trial;
            e = e_trial;
            step *= 1.1;
        } else {
            step *= 0.5;
            if step < 1e-14 {
                break;
            }
        }
    }

    let e_wall = (energy(&theta) - energy(&vec![0.0; n])) / a;
    let sigma = e_wall.abs();

    // Width from the maximum slope dtheta/dx at theta = pi/2.
    let mut best = 0.0_f64;
    let mut idx = 0;
    for i in 1..n - 1 {
        let slope = (theta[i + 1] - theta[i - 1]).abs() / (2.0 * a);
        if slope > best {
            best = slope;
            idx = i;
        }
    }
    let _ = idx;
    let delta_measured = 1.0 / best;

    DiscreteWall {
        delta_measured,
        sigma_measured: sigma,
        delta_continuum: delta_c,
        sigma_continuum: sigma_c,
        ratio_delta_over_a: delta_c / a,
    }
}

/// A 1D finite-difference micromagnetic LLG solver.
pub struct Llg1d {
    pub n: usize,
    pub dx: f64,
    pub a: f64,
    pub k: f64,
    pub ms: f64,
    pub easy: V3,
    pub h_ext: V3,
    pub alpha: f64,
    pub demag: Vec<(V3, f64)>,
    pub periodic: bool,
    pub m: Vec<V3>,
}

impl Llg1d {
    pub fn new(n: usize, dx: f64, a: f64, k: f64, ms: f64, alpha: f64) -> Self {
        Self {
            n,
            dx,
            a,
            k,
            ms,
            easy: [1.0, 0.0, 0.0],
            h_ext: [0.0, 0.0, 0.0],
            alpha,
            demag: Vec::new(),
            periodic: false,
            m: vec![[1.0, 0.0, 0.0]; n],
        }
    }

    pub fn x(&self, i: usize) -> f64 {
        i as f64 * self.dx
    }

    pub fn effective_field(&self, m: &[V3], i: usize) -> V3 {
        let (im, ip) = if self.periodic {
            ((i + self.n - 1) % self.n, (i + 1) % self.n)
        } else {
            (i.saturating_sub(1), (i + 1).min(self.n - 1))
        };
        let ex = 2.0 * self.a / (MU0 * self.ms);
        let ani = 2.0 * self.k / (MU0 * self.ms);
        let am = dot(m[i], self.easy);
        let mut h = [0.0; 3];
        for c in 0..3 {
            h[c] = ex * (m[im][c] - 2.0 * m[i][c] + m[ip][c]) / (self.dx * self.dx)
                + ani * am * self.easy[c]
                + self.h_ext[c];
        }
        for (axis, nfac) in &self.demag {
            let an = dot(m[i], *axis);
            for c in 0..3 {
                h[c] -= nfac * self.ms * an * axis[c];
            }
        }
        h
    }

    fn rhs(&self, m: &[V3]) -> Vec<V3> {
        let g = GAMMA / (1.0 + self.alpha * self.alpha);
        let mut out = vec![[0.0; 3]; self.n];
        for i in 0..self.n {
            let s = m[i];
            let h = self.effective_field(m, i);
            let s_cross_h = cross(s, h);
            let s_cross_scross_h = cross(s, s_cross_h);
            for c in 0..3 {
                out[i][c] = -g * (s_cross_h[c] + self.alpha * s_cross_scross_h[c]);
            }
        }
        out
    }

    /// One RK4 step with re-normalisation.
    pub fn step(&mut self, dt: f64) {
        let k1 = self.rhs(&self.m);
        let s2: Vec<V3> = (0..self.n)
            .map(|i| add(self.m[i], scale(k1[i], 0.5 * dt)))
            .collect();
        let k2 = self.rhs(&s2);
        let s3: Vec<V3> = (0..self.n)
            .map(|i| add(self.m[i], scale(k2[i], 0.5 * dt)))
            .collect();
        let k3 = self.rhs(&s3);
        let s4: Vec<V3> = (0..self.n)
            .map(|i| add(self.m[i], scale(k3[i], dt)))
            .collect();
        let k4 = self.rhs(&s4);
        for i in 0..self.n {
            let mut s = self.m[i];
            for c in 0..3 {
                s[c] += dt / 6.0 * (k1[i][c] + 2.0 * k2[i][c] + 2.0 * k3[i][c] + k4[i][c]);
            }
            self.m[i] = normalize(s);
        }
    }

    /// Initialise a single Bloch wall centred in the grid.
    pub fn seed_wall(&mut self) {
        let delta = wall_delta(self.a, self.k);
        let x0 = 0.5 * (self.n as f64 - 1.0) * self.dx;
        for i in 0..self.n {
            self.m[i] = wall_profile(self.x(i) - x0, delta);
        }
    }

    /// Total energy (J/m^2) of the current configuration.
    pub fn energy(&self) -> f64 {
        let mut e = 0.0;
        for i in 0..self.n {
            let mi = self.m[i];
            // Exchange with forward differences (last cell reuses the
            // backward difference so every cell contributes).
            let grad = if i + 1 < self.n {
                scale(sub(self.m[i + 1], mi), 1.0 / self.dx)
            } else if i > 0 {
                scale(sub(mi, self.m[i - 1]), 1.0 / self.dx)
            } else {
                [0.0; 3]
            };
            e += self.dx * (self.a * dot(grad, grad) - self.k * dot(mi, self.easy).powi(2));
            e -= self.dx * MU0 * self.ms * dot(self.h_ext, mi);
            for (axis, nfac) in &self.demag {
                e += self.dx * 0.5 * MU0 * self.ms * self.ms * nfac * dot(mi, *axis).powi(2);
            }
        }
        e
    }

    /// Wall energy relative to the uniform state (J/m^2).
    pub fn wall_energy(&self) -> f64 {
        let uniform = -self.k * (self.n as f64) * self.dx;
        self.energy() - uniform
    }

    /// Position (m) of the wall centre, where `m . easy` crosses zero.
    pub fn wall_position(&self) -> f64 {
        for i in 0..self.n - 1 {
            let a0 = dot(self.m[i], self.easy);
            let a1 = dot(self.m[i + 1], self.easy);
            if a0 == 0.0 {
                return self.x(i);
            }
            if a0 * a1 < 0.0 {
                let f = a0 / (a0 - a1);
                return self.x(i) + f * self.dx;
            }
        }
        f64::NAN
    }

    /// Relax to the static wall with heavy damping.
    pub fn relax(&mut self, dt: f64, steps: usize) {
        let saved = self.alpha;
        self.alpha = 1.0;
        for _ in 0..steps {
            self.step(dt);
        }
        self.alpha = saved;
    }
}

fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Field-driven wall velocity for the 1D model, `v = gamma delta H / alpha`.
pub fn walker_velocity(a: f64, k: f64, h: f64, alpha: f64) -> f64 {
    GAMMA * wall_delta(a, k) * h / alpha
}

/// Walker breakdown field of the 1D model, `H_W = alpha M_s / 2`.
pub fn walker_field(ms: f64, alpha: f64) -> f64 {
    0.5 * alpha * ms
}

/// Kittel thin-film resonance angular frequency
/// `omega = gamma sqrt((H + H_K)(H + H_K + M_s))`.
pub fn kittel_omega(h: f64, h_k: f64, ms: f64) -> f64 {
    GAMMA * ((h + h_k) * (h + h_k + ms)).sqrt()
}

/// The standard 1D Walker collective-coordinate model for the driven
/// domain wall, in wall position `q` and azimuth `phi`:
///
/// ```text
/// dq/dt  = gamma delta/(2(1+alpha^2)) [ 2 alpha H - M_s sin 2 phi ]
/// dphi/dt = -gamma/(1+alpha^2) [ H + (alpha M_s/2) sin 2 phi ]
/// ```
///
/// Below the Walker field `H_W = alpha M_s/2` the azimuth settles at
/// `sin 2 phi = -2H/(alpha M_s)` and the wall moves at
/// `v = gamma delta H / alpha`; above it the azimuth precesses and the
/// average velocity drops.
pub struct WalkerReduced {
    pub a: f64,
    pub k: f64,
    pub ms: f64,
    pub alpha: f64,
    pub q: f64,
    pub phi: f64,
}

impl WalkerReduced {
    pub fn new(a: f64, k: f64, ms: f64, alpha: f64) -> Self {
        Self {
            a,
            k,
            ms,
            alpha,
            q: 0.0,
            phi: 0.0,
        }
    }

    pub fn delta(&self) -> f64 {
        wall_delta(self.a, self.k)
    }

    pub fn step(&mut self, h: f64, dt: f64) {
        let g = GAMMA / (1.0 + self.alpha * self.alpha);
        let d = self.delta();
        let dq = g * d / 2.0 * (2.0 * self.alpha * h - self.ms * (2.0 * self.phi).sin());
        let dphi = -g * (h + 0.5 * self.alpha * self.ms * (2.0 * self.phi).sin());
        self.q += dq * dt;
        self.phi += dphi * dt;
    }

    /// Instantaneous wall velocity.
    pub fn velocity(&self, h: f64) -> f64 {
        let g = GAMMA / (1.0 + self.alpha * self.alpha);
        g * self.delta() / 2.0
            * (2.0 * self.alpha * h - self.ms * (2.0 * self.phi).sin())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wall_profile_is_the_exact_stationary_solution() {
        let (a, k) = (1.0e-11, 1.0e5);
        let delta = wall_delta(a, k);
        // dtheta/dx = sech(x/delta)/delta, so the exchange+anisotropy
        // integrand is flat: A theta'^2 = K sin^2 theta.
        for x in [-2.0, -0.5, 0.0, 0.5, 2.0] {
            let m = wall_profile(x * delta, delta);
            let theta = 2.0 * (x.exp()).atan();
            let dtheta = 1.0 / (x.cosh()) / delta;
            let lhs = a * dtheta * dtheta;
            let rhs = k * theta.sin().powi(2);
            assert!((lhs - rhs).abs() / rhs < 1e-10, "x={x}");
            assert!((m[0] - theta.cos()).abs() < 1e-12);
        }
    }
}
