//! Mie scattering by a homogeneous sphere: the exact wave-optics closure for a
//! particle, droplet or bubble larger than a molecule and smaller than the
//! coherence/geometrical-optics limit.
//!
//! `multilayer` and `pwe` treat planar and periodic matter. A single sphere of
//! arbitrary size is the remaining canonical scatterer, and it is the closure the
//! geometrical-optics `radiation` crate needs for a participating medium:
//! extinction, scattering, absorption and the asymmetry parameter `g` become the
//! per-metre coefficients `sigma_ext`, `sigma_sca`, `sigma_a` and `g` of a
//! [`radiation::Medium`](../../radiation/index.html).
//!
//! # Exact model
//!
//! For a sphere of radius `a` and relative refractive index `m = m_particle /
//! n_medium` (`N = n + i kappa`, positive imaginary part is absorption, matching
//! `complex::C`), the size parameter is `x = 2 pi n_medium a / lambda`. The Mie
//! coefficients are
//!
//! ```text
//!   a_n = [ (D_n(mx)/m + n/x) psi_n(x) - psi_{n-1}(x) ]
//!       / [ (D_n(mx)/m + n/x) xi_n(x) - xi_{n-1}(x) ]
//!   b_n = [ (m D_n(mx) + n/x) psi_n(x) - psi_{n-1}(x) ]
//!       / [ (m D_n(mx) + n/x) xi_n(x) - xi_{n-1}(x) ]
//! ```
//!
//! with the Riccati–Bessel functions `psi_n = x j_n`, `chi_n = -x y_n`,
//! `xi_n = psi_n + i chi_n`, and the logarithmic derivative
//! `D_n(z) = psi_n'(z)/psi_n(z)` (Bohren & Huffman 1983, eqs. 4.52, 4.88).
//!
//! Stability follows Wiscombe (1980):
//! * `D_n(mx)` from the **downward** recurrence (started from the Lentz continued
//!   fraction), which is the minimal solution and never overflows;
//! * `psi_n(x)` from Miller's **downward** recurrence with rescaling;
//! * `chi_n(x)` from the **upward** recurrence, the growing solution;
//! * the series is truncated at `n_stop = x + 4.05 x^{1/3} + 2` terms.
//!
//! # Efficiencies
//!
//! `Q_ext`, `Q_sca`, `Q_abs = Q_ext - Q_sca`, `Q_back` and the asymmetry
//! parameter `g` (Bohren & Huffman, eqs. 4.62–4.64, 4.71):
//!
//! ```text
//!   Q_ext = (2/x^2) sum (2n+1) Re(a_n + b_n)
//!   Q_sca = (2/x^2) sum (2n+1) (|a_n|^2 + |b_n|^2)
//!   Q_back = (1/x^2) |sum (2n+1) (-1)^n (a_n - b_n)|^2
//!   g = (4/(x^2 Q_sca)) sum [ n(n+2)/(n+1) Re(a_n a*_{n+1} + b_n b*_{n+1})
//!                            + (2n+1)/(n(n+1)) Re(a_n b*_n) ]
//! ```
//!
//! The angular amplitudes `S_1`, `S_2` use the `pi_n`, `tau_n` recurrence
//! (Bohren & Huffman, eq. 4.47). The unpolarised phase function
//! `p(theta) = 2 (|S_1|^2 + |S_2|^2) / (x^2 Q_sca)` is normalised to
//! `(1/4 pi) integral p dOmega = 1`.
//!
//! # Folds and their indicators
//!
//! The exact series is cheap for a single particle, but a size distribution or a
//! coupled transport solve wants a closed form. Two folds are provided, both with
//! the indicator that bounds them:
//!
//! | fold | valid when | indicator | error |
//! |---|---|---|---|
//! | Rayleigh | `x << 1` and `|m|x << 1` | `x` (and `|m|x`) | relative O(x^2) |
//! | anomalous diffraction | `x >> 1`, `|m-1| << 1` | `rho = 2x|m-1|` | O(1/x) |
//! | extinction paradox | `x >> 1` | `1/x` | `Q_ext -> 2` |
//!
//! # Material constants
//!
//! Water and silica indices come from [`crate::dispersion`] (Daimon–Masumura
//! 2007 and Malitson 1965). Particle indices not from that table are cited at the
//! call site or in the tests and labelled "estimate".

use crate::complex::C;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Riccati-Bessel machinery
// ---------------------------------------------------------------------------

/// Wiscombe (1980) series truncation `n_stop = x + 4.05 x^{1/3} + 2`.
///
/// The truncation error of the summed efficiencies is then about `1e-6`.
pub fn wiscombe_terms(x: f64) -> usize {
    if x <= 0.0 {
        return 0;
    }
    (x + 4.05 * x.cbrt() + 2.0) as usize
}

/// Logarithmic derivative `D_n(z) = psi_n'(z)/psi_n(z)` by the Lentz continued
/// fraction (Lentz, Appl. Opt. **15**, 668 (1976)).
fn lentz_dn(z: C, n: usize) -> C {
    let mut zinv = C::real(2.0) / z;
    let alpha = zinv * (n as f64 + 0.5);
    let mut aj = zinv * -(n as f64 + 1.5);
    let mut alpha_j1 = alpha.recip() + aj;
    let mut alpha_j2 = aj;
    let mut ratio = alpha_j1 / alpha_j2;
    let mut runratio = alpha * ratio;
    while (ratio.abs() - 1.0).abs() > 1e-14 {
        aj = zinv - aj;
        alpha_j1 = alpha_j1.recip() + aj;
        alpha_j2 = alpha_j2.recip() + aj;
        ratio = alpha_j1 / alpha_j2;
        zinv = -zinv;
        runratio = ratio * runratio;
    }
    -C::real(n as f64) / z + runratio
}

/// `D_0 .. D_nmax` at `z`, by the downward recurrence seeded with [`lentz_dn`].
///
/// The downward recurrence is stable for any complex `z` (Wiscombe 1980); the
/// result is shifted so that `d[k] = D_k`.
fn log_derivatives(z: C, n_max: usize) -> Vec<C> {
    let mut d = vec![C::default(); n_max + 1];
    d[n_max] = lentz_dn(z, n_max);
    for n in (1..=n_max).rev() {
        let nz = C::real(n as f64) / z;
        d[n - 1] = nz - (d[n] + nz).recip();
    }
    d
}

/// `psi_0 .. psi_n` at the real argument `x` by Miller's downward recurrence.
///
/// The downward (minimal) solution is stable; the arbitrary scale is fixed
/// against the analytic `psi_0 = sin x`, `psi_1 = sin x / x - cos x`.
fn riccati_psi(x: f64, n_terms: usize) -> Vec<f64> {
    let n_start = n_terms.max(x as usize) + 25 + (1.5 * x.sqrt()) as usize;
    let mut p = vec![0.0f64; n_start + 2];
    p[n_start] = 1e-50;
    for n in (1..=n_start).rev() {
        p[n - 1] = (2 * n + 1) as f64 / x * p[n] - p[n + 1];
        if p[n - 1].abs() > 1e200 {
            for v in p[n - 1..=n_start + 1].iter_mut() {
                *v /= 1e200;
            }
        }
    }
    let psi0 = x.sin();
    let psi1 = x.sin() / x - x.cos();
    let scale = if psi1.abs() > psi0.abs() {
        psi1 / p[1]
    } else {
        psi0 / p[0]
    };
    for v in p.iter_mut() {
        *v *= scale;
    }
    p.truncate(n_terms + 1);
    p
}

// ---------------------------------------------------------------------------
// Mie coefficients and efficiencies
// ---------------------------------------------------------------------------

/// Mie coefficients `a_n`, `b_n` for a homogeneous sphere, plus the `x`, `m`
/// used. Indices are zero-based: `a[0]` is `a_1`.
#[derive(Clone, Debug)]
pub struct MieScattering {
    pub x: f64,
    /// Relative complex index, `N = n + i kappa` (absorption is `kappa > 0`).
    pub m: C,
    pub a: Vec<C>,
    pub b: Vec<C>,
}

/// Integrated Mie efficiencies for a homogeneous sphere.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Efficiencies {
    pub q_ext: f64,
    pub q_sca: f64,
    pub q_abs: f64,
    /// Radar/backscattering efficiency `Q_back`.
    pub q_back: f64,
    /// Asymmetry parameter `g = <cos theta>`.
    pub g: f64,
}

impl Efficiencies {
    /// Radiation-pressure efficiency `Q_pr = Q_ext - g Q_sca`.
    pub fn q_pr(&self) -> f64 {
        self.q_ext - self.g * self.q_sca
    }
}

impl MieScattering {
    /// Build the coefficients for size parameter `x` and relative index `m`.
    ///
    /// `m` follows the crate convention `N = n + i kappa`. The kernel internally
    /// uses the `N = n - i kappa` convention of Bohren & Huffman, so the
    /// imaginary part is negated here; the physical efficiencies are unchanged.
    pub fn new(x: f64, m: C) -> Self {
        Self::with_terms(x, m, wiscombe_terms(x))
    }

    /// As [`MieScattering::new`] with an explicit series length `n_terms`
    /// (for truncation studies; [`wiscombe_terms`] is the default).
    pub fn with_terms(x: f64, m: C, n_terms: usize) -> Self {
        let kernel_m = C::new(m.re, -m.im);
        if x <= 0.0 || n_terms == 0 {
            return Self {
                x,
                m,
                a: Vec::new(),
                b: Vec::new(),
            };
        }
        let inv_x = 1.0 / x;
        let psi = riccati_psi(x, n_terms);
        let d = log_derivatives(kernel_m * x, n_terms);

        let mut a = vec![C::default(); n_terms];
        let mut b = vec![C::default(); n_terms];
        let mut chi_nm1 = x.cos();
        let mut chi_n = chi_nm1 * inv_x + x.sin();
        for n in 1..=n_terms {
            let n_over_x = n as f64 * inv_x;
            let xi_nm1 = C::new(psi[n - 1], chi_nm1);
            let xi_n = C::new(psi[n], chi_n);
            let dn = d[n];
            let temp = dn / kernel_m + n_over_x;
            a[n - 1] = (temp * psi[n] - psi[n - 1]) / (temp * xi_n - xi_nm1);
            let temp = dn * kernel_m + n_over_x;
            b[n - 1] = (temp * psi[n] - psi[n - 1]) / (temp * xi_n - xi_nm1);
            let chi = (2 * n + 1) as f64 * inv_x * chi_n - chi_nm1;
            chi_nm1 = chi_n;
            chi_n = chi;
        }
        Self { x, m, a, b }
    }

    /// Build the coefficients for a **perfectly conducting** sphere (`|m| -> inf`),
    /// the limit in which the internal field vanishes. The coefficients are
    /// `a_n = (n/x psi_n - psi_{n-1})/(n/x xi_n - xi_{n-1})`, `b_n = psi_n/xi_n`.
    ///
    /// This is distinct from `new(x, C::new(0, 100))`, which is a *finite* index:
    /// the general series then differs from this limit by `O(1/|m|)` (about 0.3%
    /// for `|m| = 100` at `x = 100`).
    pub fn perfect_conductor(x: f64) -> Self {
        let n_terms = wiscombe_terms(x);
        if x <= 0.0 || n_terms == 0 {
            return Self {
                x,
                m: C::new(f64::INFINITY, 0.0),
                a: Vec::new(),
                b: Vec::new(),
            };
        }
        let inv_x = 1.0 / x;
        let psi = riccati_psi(x, n_terms);
        let mut a = vec![C::default(); n_terms];
        let mut b = vec![C::default(); n_terms];
        let mut chi_nm1 = x.cos();
        let mut chi_n = chi_nm1 * inv_x + x.sin();
        for n in 1..=n_terms {
            let n_over_x = n as f64 * inv_x;
            let xi_nm1 = C::new(psi[n - 1], chi_nm1);
            let xi_n = C::new(psi[n], chi_n);
            a[n - 1] = (C::real(psi[n]) * n_over_x - psi[n - 1]) / (xi_n * n_over_x - xi_nm1);
            b[n - 1] = C::real(psi[n]) / xi_n;
            let chi = (2 * n + 1) as f64 * inv_x * chi_n - chi_nm1;
            chi_nm1 = chi_n;
            chi_n = chi;
        }
        Self {
            x,
            m: C::new(f64::INFINITY, 0.0),
            a,
            b,
        }
    }

    /// Build from a physical radius and wavelength (metres).
    ///
    /// `m_particle` is the particle index `N = n + i kappa`; `n_medium` is the
    /// (real) ambient index. The size parameter is `x = 2 pi n_medium a / lambda`.
    pub fn from_radius(radius: f64, wavelength: f64, m_particle: C, n_medium: f64) -> Self {
        let x = 2.0 * PI * n_medium * radius / wavelength;
        let m = m_particle / n_medium;
        Self::new(x, m)
    }

    /// Number of Mie terms in the truncated series.
    pub fn n_terms(&self) -> usize {
        self.a.len()
    }

    /// `Q_ext`, `Q_sca`, `Q_abs`, `Q_back` and `g`.
    pub fn efficiencies(&self) -> Efficiencies {
        let x2 = self.x * self.x;
        if self.a.is_empty() || x2 == 0.0 {
            return Efficiencies::default();
        }
        let n = self.a.len();
        let mut q_ext = 0.0;
        let mut q_sca = 0.0;
        let mut g = 0.0;
        let mut back = C::default();
        for i in 0..n {
            let nn = (i + 1) as f64;
            let cn = 2.0 * nn + 1.0;
            q_ext += cn * (self.a[i].re + self.b[i].re);
            q_sca += cn * (self.a[i].abs2() + self.b[i].abs2());
            if i + 1 < n {
                let cross = self.a[i] * self.a[i + 1].conj() + self.b[i] * self.b[i + 1].conj();
                g += nn * (nn + 2.0) / (nn + 1.0) * cross.re;
            }
            g += cn / (nn * (nn + 1.0)) * (self.a[i] * self.b[i].conj()).re;
            let sign = if (i + 1) % 2 == 0 { 1.0 } else { -1.0 };
            back = back + (self.a[i] - self.b[i]) * (sign * cn);
        }
        q_ext *= 2.0 / x2;
        q_sca *= 2.0 / x2;
        let q_back = back.abs2() / x2;
        let g = if q_sca > 0.0 {
            4.0 / (x2 * q_sca) * g
        } else {
            0.0
        };
        Efficiencies {
            q_ext,
            q_sca,
            q_abs: q_ext - q_sca,
            q_back,
            g,
        }
    }

    /// Angular amplitudes `(S_1, S_2)` at `mu = cos theta`.
    ///
    /// These are the Bohren & Huffman amplitudes (the `S1_S2` "wiscombe"
    /// normalisation of `miepython`), so the unpolarised differential cross
    /// section is `(|S_1|^2 + |S_2|^2)/(2 k^2)`.
    pub fn s1_s2(&self, mu: f64) -> (C, C) {
        let mut s1 = C::default();
        let mut s2 = C::default();
        let mut pi_nm1 = 0.0f64;
        let mut pi_n = 1.0f64;
        for i in 0..self.a.len() {
            let n = (i + 1) as f64;
            let tau_n = n * mu * pi_n - (n + 1.0) * pi_nm1;
            let scale = (2.0 * n + 1.0) / (n * (n + 1.0));
            s1 = s1 + (self.a[i] * pi_n + self.b[i] * tau_n) * scale;
            s2 = s2 + (self.a[i] * tau_n + self.b[i] * pi_n) * scale;
            let pi_np1 = ((2.0 * n + 1.0) * mu * pi_n - (n + 1.0) * pi_nm1) / n;
            pi_nm1 = pi_n;
            pi_n = pi_np1;
        }
        (s1, s2)
    }

    /// Unpolarised phase function at scattering angle `theta` (radians),
    /// normalised to `(1/4 pi) integral p dOmega = 1`.
    pub fn phase_function(&self, theta: f64) -> f64 {
        self.phase_function_with(theta, self.efficiencies().q_sca)
    }

    /// Phase function with `Q_sca` already known (avoids recomputing the series).
    pub fn phase_function_with(&self, theta: f64, q_sca: f64) -> f64 {
        if q_sca <= 0.0 {
            return 0.0;
        }
        let (s1, s2) = self.s1_s2(theta.cos());
        2.0 * (s1.abs2() + s2.abs2()) / (self.x * self.x * q_sca)
    }
}

/// Convenience: efficiencies for a homogeneous sphere.
pub fn efficiencies(x: f64, m: C) -> Efficiencies {
    MieScattering::new(x, m).efficiencies()
}

// ---------------------------------------------------------------------------
// Folds: Rayleigh and anomalous diffraction
// ---------------------------------------------------------------------------

/// Rayleigh (small-particle) scattering efficiency
/// `Q_sca = (8/3) x^4 |(m^2-1)/(m^2+2)|^2`.
///
/// Validity indicator: `x << 1` and `|m| x << 1`; the relative error of the fold
/// against exact Mie is `O(x^2)`.
pub fn rayleigh_qsca(x: f64, m: C) -> f64 {
    let m2 = m * m;
    let k = (m2 - 1.0) / (m2 + 2.0);
    (8.0 / 3.0) * x.powi(4) * k.abs2()
}

/// Rayleigh extinction efficiency, including absorption (Bohren & Huffman
/// eq. 5.9): `Q_ext = 4 x Im[(m^2-1)/(m^2+2)] + (8/3) x^4 |K|^2`.
pub fn rayleigh_qext(x: f64, m: C) -> f64 {
    let m2 = m * m;
    let k = (m2 - 1.0) / (m2 + 2.0);
    4.0 * x * k.im + rayleigh_qsca(x, m)
}

/// Rayleigh unpolarised phase function `(3/4)(1 + cos^2 theta)`.
pub fn rayleigh_phase(theta: f64) -> f64 {
    let mu = theta.cos();
    0.75 * (1.0 + mu * mu)
}

/// Rayleigh validity indicator: the largest of `x` and `|m| x`. The fold error
/// is roughly this squared.
pub fn rayleigh_indicator(x: f64, m: C) -> f64 {
    x.max(m.abs() * x)
}

/// van de Hulst anomalous-diffraction extinction efficiency for a non-absorbing
/// sphere, `Q_ext = 2 - (4/rho) sin rho + (4/rho^2)(1 - cos rho)` with the phase
/// delay `rho = 2 x (m - 1)`.
///
/// Intermediate fold: `x >> 1`, `|m - 1| << 1`. Indicator `rho`; for a strongly
/// absorbing particle pass `rho = 2 x (m - 1)` as a complex number and take the
/// real part of the analytic continuation, or use the exact series.
pub fn anomalous_diffraction_qext_real(x: f64, m: f64) -> f64 {
    let rho = 2.0 * x * (m - 1.0);
    anomalous_diffraction_qext_from_rho(rho)
}

/// Anomalous-diffraction `Q_ext` from the phase delay `rho` directly.
pub fn anomalous_diffraction_qext_from_rho(rho: f64) -> f64 {
    if rho == 0.0 {
        return 0.0;
    }
    2.0 - 4.0 / rho * rho.sin() + 4.0 / (rho * rho) * (1.0 - rho.cos())
}

// ---------------------------------------------------------------------------
// Size distributions and participating-medium properties
// ---------------------------------------------------------------------------

/// Gamma size distribution with effective radius `r_eff` and effective variance
/// `v_eff` (Hansen & Travis 1974):
/// `n(r) ~ r^alpha exp(-r/(v_eff r_eff))`, `alpha = 1/v_eff - 3`, normalised so
/// that `integral n dr = 1`. The area-weighted radius `r_eff = <r^3>/<r^2>`.
///
/// Cloud droplets use `v_eff ~ 0.1` (dispersion ~ 0.35).
#[derive(Clone, Copy, Debug)]
pub struct GammaSizeDistribution {
    pub r_eff: f64,
    pub v_eff: f64,
}

impl GammaSizeDistribution {
    /// Stratocumulus-like distribution: `r_eff = 10 µm`, `v_eff = 0.1`.
    pub fn stratocumulus() -> Self {
        Self {
            r_eff: 10.0e-6,
            v_eff: 0.1,
        }
    }
    pub fn alpha(&self) -> f64 {
        1.0 / self.v_eff - 3.0
    }
    pub fn beta(&self) -> f64 {
        1.0 / (self.v_eff * self.r_eff)
    }
    /// Normalised number density per metre of radius, at radius `r` (m).
    pub fn pdf(&self, r: f64) -> f64 {
        let a = self.alpha();
        let b = self.beta();
        if r <= 0.0 {
            return if a == 0.0 { b } else { 0.0 };
        }
        (b.powf(a + 1.0) / gamma_fn(a + 1.0)) * r.powf(a) * (-b * r).exp()
    }
    /// `p`-th raw moment `~ beta^{-p} Gamma(alpha+1+p)/Gamma(alpha+1)`.
    pub fn moment(&self, p: f64) -> f64 {
        (ln_gamma(self.alpha() + 1.0 + p) - ln_gamma(self.alpha() + 1.0)).exp()
            / self.beta().powf(p)
    }
}

/// Lanczos log-gamma, `g = 7`, `N = 9` (Numerical Recipes). Accurate to ~1e-14
/// for the positive half-line used here.
pub fn ln_gamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_93,
        676.520_368_121_885_1,
        -1259.139_216_722_402_8,
        771.323_428_777_653_13,
        -176.615_029_162_140_59,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // Reflection: Gamma(x) Gamma(1-x) = pi / sin(pi x).
        return (PI).ln() - (PI * x).sin().abs().ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = C[0];
    let t = x + G + 0.5;
    for (i, &c) in C.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// Gamma function on the positive half-line.
pub fn gamma_fn(x: f64) -> f64 {
    ln_gamma(x).exp()
}

/// Size-averaged optical properties. `c_*` are per-particle mean cross sections
/// (m^2); `q_*` are their area-weighted dimensionless efficiencies.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SizeAveraged {
    pub c_ext: f64,
    pub c_sca: f64,
    pub c_abs: f64,
    pub q_ext: f64,
    pub q_sca: f64,
    pub q_abs: f64,
    pub g: f64,
}

/// Composite Simpson quadrature on `[a, b]`, `n` even sub-intervals.
fn simpson<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, n: usize) -> f64 {
    let n = if n % 2 == 0 { n } else { n + 1 };
    let h = (b - a) / n as f64;
    let mut s = f(a) + f(b);
    for i in 1..n {
        let x = a + i as f64 * h;
        s += if i % 2 == 0 { 2.0 } else { 4.0 } * f(x);
    }
    s * h / 3.0
}

/// Mie properties averaged over a gamma size distribution.
///
/// `m` is the relative index at `wavelength` (metres), `n_medium` the ambient
/// index. `n_quad` is the number of Simpson sub-intervals (2000 is ample). The
/// returned `c_*` are per-particle mean cross sections.
pub fn gamma_averaged(
    dist: &GammaSizeDistribution,
    m: C,
    wavelength: f64,
    n_medium: f64,
    n_quad: usize,
) -> SizeAveraged {
    let r_max = dist.r_eff * (1.0 + 8.0 * dist.v_eff.sqrt());
    let k = 2.0 * PI * n_medium / wavelength;
    let area = |r: f64| PI * r * r;
    let pdf = |r: f64| dist.pdf(r);

    let c_ext = simpson(
        |r| efficiencies(k * r, m).q_ext * area(r) * pdf(r),
        0.0,
        r_max,
        n_quad,
    );
    let c_sca = simpson(
        |r| efficiencies(k * r, m).q_sca * area(r) * pdf(r),
        0.0,
        r_max,
        n_quad,
    );
    let c_abs = simpson(
        |r| efficiencies(k * r, m).q_abs * area(r) * pdf(r),
        0.0,
        r_max,
        n_quad,
    );
    let mean_area = simpson(|r| area(r) * pdf(r), 0.0, r_max, n_quad);
    let g_num = simpson(
        |r| {
            let e = efficiencies(k * r, m);
            e.g * e.q_sca * area(r) * pdf(r)
        },
        0.0,
        r_max,
        n_quad,
    );

    SizeAveraged {
        c_ext,
        c_sca,
        c_abs,
        q_ext: c_ext / mean_area,
        q_sca: c_sca / mean_area,
        q_abs: c_abs / mean_area,
        g: if c_sca > 0.0 { g_num / c_sca } else { 0.0 },
    }
}

/// Per-metre radiative properties of a cloud of particles, exactly the inputs
/// `radiation::Medium` consumes: `sigma_ext`, `sigma_sca`, `sigma_abs` (m^-1)
/// and the asymmetry parameter `g`. The single-scattering albedo is
/// `sigma_sca / sigma_ext`; the ambient index is `n_medium`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ParticipatingMedium {
    pub sigma_ext: f64,
    pub sigma_sca: f64,
    pub sigma_abs: f64,
    pub g: f64,
    pub n_medium: f64,
}

impl ParticipatingMedium {
    pub fn albedo(&self) -> f64 {
        if self.sigma_ext > 0.0 {
            self.sigma_sca / self.sigma_ext
        } else {
            0.0
        }
    }
}

/// Radiative properties of a monodisperse population of number density `n_number`
/// (m^-3) at `wavelength`.
pub fn monodisperse_medium(
    x: f64,
    m: C,
    wavelength: f64,
    n_medium: f64,
    n_number: f64,
) -> ParticipatingMedium {
    let e = efficiencies(x, m);
    let radius = x * wavelength / (2.0 * PI * n_medium);
    let c = PI * radius * radius;
    ParticipatingMedium {
        sigma_ext: n_number * c * e.q_ext,
        sigma_sca: n_number * c * e.q_sca,
        sigma_abs: n_number * c * e.q_abs,
        g: e.g,
        n_medium,
    }
}

/// Radiative properties of a gamma-distributed cloud of number density `n_number`
/// (m^-3). This is the cloud/droplet/bubble closure for `radiation`.
pub fn cloud_medium(
    dist: &GammaSizeDistribution,
    m: C,
    wavelength: f64,
    n_medium: f64,
    n_number: f64,
    n_quad: usize,
) -> ParticipatingMedium {
    let s = gamma_averaged(dist, m, wavelength, n_medium, n_quad);
    ParticipatingMedium {
        sigma_ext: n_number * s.c_ext,
        sigma_sca: n_number * s.c_sca,
        sigma_abs: n_number * s.c_abs,
        g: s.g,
        n_medium,
    }
}

// ---------------------------------------------------------------------------
// Optical forces
// ---------------------------------------------------------------------------

/// Optical-force estimates for a Rayleigh-regime sphere in a focused Gaussian
/// beam.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct OpticalForce {
    /// Exact Mie radiation-pressure efficiency `Q_pr = Q_ext - g Q_sca`.
    pub q_pr: f64,
    /// Radiation-pressure (scattering) force on the beam axis, newtons.
    pub scattering_force: f64,
    /// Transverse gradient force at the steepest part of the Gaussian, newtons.
    pub gradient_force: f64,
    /// Trap stiffness near focus `dF/dr`, N/m (Rayleigh estimate).
    pub stiffness: f64,
    /// Potential depth `(1/2) k w0^2` in units of `k_B T`.
    pub depth_over_kbt: f64,
}

/// `k_B` [J/K] (CODATA 2018 exact).
pub const K_B: f64 = 1.380_649e-23;

/// Estimate optical forces on a sphere in a Gaussian trap.
///
/// * `radius`, `wavelength`, `n_medium`: geometry and ambient index.
/// * `m`: relative index (Rayleigh Clausius–Mossotti factor `K = (m^2-1)/(m^2+2)`).
/// * `power`: beam power (W).
/// * `w0`: beam waist (m). The Gaussian has `I(r) = I0 exp(-2 r^2/w0^2)`,
///   `I0 = 2 P/(pi w0^2)`.
/// * `temperature`: for the `k_B T` comparison.
///
/// The scattering force uses the exact Mie `Q_pr`. The gradient force and
/// stiffness use the standard Rayleigh result
/// `F_grad = (2 pi n_medium a^3/c) K grad I` (Ashkin 1986; Neuman & Block 2004).
/// For a 1 µm bead this is a marginal approximation (`x ~ 4`) and is reported as
/// an order-of-magnitude estimate; `Q_pr` is exact.
pub fn optical_force(
    radius: f64,
    wavelength: f64,
    n_medium: f64,
    m: C,
    power: f64,
    w0: f64,
    temperature: f64,
) -> OpticalForce {
    let x = 2.0 * PI * n_medium * radius / wavelength;
    let e = efficiencies(x, m);
    let q_pr = e.q_pr();
    let i0 = 2.0 * power / (PI * w0 * w0);
    // Radiation pressure on-axis: (n_medium/c) Q_pr times the power intercepted
    // by the sphere. Integrating the Gaussian over the sphere gives
    // P (1 - exp(-2 a^2/w0^2)); for a << w0 this is the usual I0 * pi a^2.
    let overlap = 1.0 - (-2.0 * radius * radius / (w0 * w0)).exp();
    let scattering_force = n_medium / crate::C0 * q_pr * power * overlap;
    // Rayleigh gradient force at the steepest Gaussian slope (~1.21 I0/w0).
    let m2 = m * m;
    let k_cm = ((m2 - 1.0) / (m2 + 2.0)).re;
    let grad_i = 2.0 * i0 / w0; // conservative |grad I| scale
    let pref = 2.0 * PI * n_medium * radius.powi(3) * k_cm / crate::C0;
    let gradient_force = pref * grad_i;
    // d/dr of -(4 r/w0^2) I0 at the centre gives the centre stiffness.
    let stiffness = pref * 4.0 * i0 / (w0 * w0);
    let depth = 0.5 * stiffness * w0 * w0;
    OpticalForce {
        q_pr,
        scattering_force,
        gradient_force,
        stiffness,
        depth_over_kbt: depth / (K_B * temperature),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dispersion::Sellmeier;
    use crate::Dispersion;

    fn close(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "{a} vs {b} (tol {tol})");
    }

    /// Wiscombe MIEV0 test cases, as reproduced by the `miepython` project's
    /// reference suite. The second column is the published value.
    struct Case {
        m: C,
        x: f64,
        q_ext: f64,
        q_sca: f64,
        q_back: f64,
        g: f64,
    }

    const CASES: &[Case] = &[
        // m = 1.55 (Bohren & Huffman dielectric sphere), x = 2 pi 0.525/0.6328
        Case {
            m: C::new(1.55, 0.0),
            x: 5.212819668567135,
            q_ext: 3.10543,
            q_sca: 3.10543,
            q_back: 2.92534,
            g: 0.63314,
        },
        // m = 1.55 - 0.1i
        Case {
            m: C::new(1.55, 0.1),
            x: 5.212819668567135,
            q_ext: 2.86165188243,
            q_sca: 1.66424911991,
            q_back: 0.20599534080,
            g: 0.80128972639,
        },
        // MIEV0 case 5, m = 0.75, x = 0.099
        Case {
            m: C::new(0.75, 0.0),
            x: 0.099,
            q_ext: 0.000007,
            q_sca: 0.000007,
            q_back: 0.0,
            g: 0.001448,
        },
        // MIEV0 case 7, m = 0.75, x = 10
        Case {
            m: C::new(0.75, 0.0),
            x: 10.0,
            q_ext: 2.232265,
            q_sca: 2.232265,
            q_back: 0.0,
            g: 0.896473,
        },
        // MIEV0 case 8, m = 0.75, x = 1000
        Case {
            m: C::new(0.75, 0.0),
            x: 1000.0,
            q_ext: 1.997908,
            q_sca: 1.997908,
            q_back: 0.0,
            g: 0.844944,
        },
        // MIEV0 water absorbing case 10, m = 1.33 - 1e-5 i, x = 100
        // (the published case fixes Q_sca and g; Q_ext is not tabulated)
        Case {
            m: C::new(1.33, 0.00001),
            x: 100.0,
            q_ext: 0.0,
            q_sca: 2.096594,
            q_back: 0.0,
            g: 0.868959,
        },
        // MIEV0 water absorbing case 11, x = 10000 (stability)
        Case {
            m: C::new(1.33, 0.00001),
            x: 10000.0,
            q_ext: 0.0,
            q_sca: 1.723857,
            q_back: 0.0,
            g: 0.907840,
        },
        // MIEV0 case 14, m = 1.5 - 1i, x = 1
        // (the published case fixes Q_ext and Q_sca)
        Case {
            m: C::new(1.5, 1.0),
            x: 1.0,
            q_ext: 2.336321,
            q_sca: 0.6634538,
            q_back: 0.0,
            g: 0.192136,
        },
        // MIEV0 case 15, m = 1.5 - 1i, x = 100
        Case {
            m: C::new(1.5, 1.0),
            x: 100.0,
            q_ext: 2.097502,
            q_sca: 1.283697,
            q_back: 0.0,
            g: 0.850252,
        },
        // MIEV0 case 16, m = 1.5 - 1i, x = 10000
        Case {
            m: C::new(1.5, 1.0),
            x: 10000.0,
            q_ext: 2.004368,
            q_sca: 1.236575,
            q_back: 0.0,
            g: 0.846309,
        },
        // MIEV0 case 17, m = 10 - 10i, x = 1
        Case {
            m: C::new(10.0, 10.0),
            x: 1.0,
            q_ext: 0.0,
            q_sca: 2.049405,
            q_back: 0.0,
            g: -0.110664,
        },
        // MIEV0 case 18, m = 10 - 10i, x = 100
        Case {
            m: C::new(10.0, 10.0),
            x: 100.0,
            q_ext: 0.0,
            q_sca: 1.836785,
            q_back: 0.0,
            g: 0.556215,
        },
        // MIEV0 case 19, m = 10 - 10i, x = 10000
        Case {
            m: C::new(10.0, 10.0),
            x: 10000.0,
            q_ext: 0.0,
            q_sca: 1.795393,
            q_back: 0.0,
            g: 0.548194,
        },
    ];

    #[test]
    fn wiscombe_miev0_benchmarks() {
        let mut worst = 0.0f64;
        for c in CASES {
            let e = efficiencies(c.x, c.m);
            // Compare the quantities that the published case fixes.
            if c.q_ext != 0.0 {
                let err = (e.q_ext - c.q_ext).abs();
                worst = worst.max(err);
                assert!(
                    err < 1e-5,
                    "Qext m={:?} x={} got {} want {}",
                    c.m,
                    c.x,
                    e.q_ext,
                    c.q_ext
                );
            }
            let err = (e.q_sca - c.q_sca).abs();
            worst = worst.max(err);
            assert!(
                err < 1e-5,
                "Qsca m={:?} x={} got {} want {}",
                c.m,
                c.x,
                e.q_sca,
                c.q_sca
            );
            let err = (e.g - c.g).abs();
            worst = worst.max(err);
            assert!(
                err < 2e-5,
                "g m={:?} x={} got {} want {}",
                c.m,
                c.x,
                e.g,
                c.g
            );
            if c.q_back != 0.0 {
                let err = (e.q_back - c.q_back).abs();
                worst = worst.max(err);
                assert!(
                    err < 1e-5,
                    "Qback m={:?} x={} got {} want {}",
                    c.m,
                    c.x,
                    e.q_back,
                    c.q_back
                );
            }
        }
        // The benchmark suite spans x = 0.099 to 10000 and |m| up to 14.
        assert!(worst < 1e-5, "worst error {worst}");
    }

    #[test]
    fn perfectly_conducting_sphere() {
        // Wiscombe MIEV0 perfectly-conducting cases (the `|m| -> inf` limit).
        for &(x, qsca, g) in &[
            (0.101, 0.000348, -0.397262),
            (100.0, 2.008102, 0.500926),
            (10000.0, 2.000289, 0.500070),
        ] {
            let e = MieScattering::perfect_conductor(x).efficiencies();
            close(e.q_sca, qsca, 1e-5);
            close(e.g, g, 1e-5);
            // A conductor is lossless in the scattering sense: Q_ext = Q_sca.
            close(e.q_ext, e.q_sca, 1e-10 * e.q_sca.max(1.0));
        }
    }

    #[test]
    fn bohren_huffman_backscatter() {
        // m = 1.55, x = 2 pi 0.525/0.6328: Q_back = 2.92534.
        let x = 2.0 * PI * 0.525 / 0.6328;
        let e = efficiencies(x, C::new(1.55, 0.0));
        close(e.q_back, 2.92534, 1e-5);
        // m = 1.55 - 0.1i: Q_back = 0.20599534080.
        let e = efficiencies(x, C::new(1.55, 0.1));
        close(e.q_back, 0.20599534080, 1e-7);
    }

    #[test]
    fn optical_theorem_equals_forward_amplitude() {
        // Q_ext = (4/x^2) Re S(0), with S1(0) = S2(0).
        for &(m, x) in &[
            (C::new(1.5, 0.0), 2.0),
            (C::new(1.5, 1.0), 5.0),
            (C::new(1.33, 0.00001), 20.0),
            (C::new(0.75, 0.0), 50.0),
        ] {
            let s = MieScattering::new(x, m);
            let e = s.efficiencies();
            let (s1, _) = s.s1_s2(1.0);
            let q_opt = 4.0 / (x * x) * s1.re;
            close(e.q_ext, q_opt, 1e-12 * e.q_ext.abs().max(1.0));
        }
    }

    #[test]
    fn energy_conservation_lossless() {
        for &m in &[1.05, 1.2, 1.33, 1.5, 2.0, 0.75] {
            for &x in &[0.1, 1.0, 5.0, 20.0, 100.0, 500.0] {
                let e = efficiencies(x, C::new(m, 0.0));
                close(e.q_ext, e.q_sca, 1e-11 * e.q_ext.max(1.0));
                close(e.q_abs, 0.0, 1e-11 * e.q_ext.max(1.0));
            }
        }
    }

    #[test]
    fn phase_function_normalisation() {
        // (1/4 pi) integral p dOmega = 1, checked with Simpson in mu.
        let s = MieScattering::new(5.0, C::new(1.5, 0.1));
        let q = s.efficiencies().q_sca;
        let n = 20000;
        let integral = simpson(
            |mu| {
                let theta = mu.acos();
                2.0 * PI * s.phase_function_with(theta, q)
            },
            -1.0,
            1.0,
            n,
        );
        close(integral / (4.0 * PI), 1.0, 1e-6);
    }

    #[test]
    fn angular_amplitudes_match_published_values() {
        // miepython / Bohren & Huffman, m = 1.5 - 1i, x = 1, "wiscombe" norm.
        let s = MieScattering::new(1.0, C::new(1.5, 1.0));
        let theta = [0.0f64, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0];
        for (i, &t) in theta.iter().enumerate() {
            let (s1, s2) = s.s1_s2(t.to_radians().cos());
            let (re1, im1, re2, im2) = match i {
                0 => (0.584080, 0.190515, 0.584080, 0.190515),
                1 => (0.565702, 0.187200, 0.500161, 0.145611),
                2 => (0.517525, 0.178443, 0.287964, 0.041054),
                3 => (0.456340, 0.167167, 0.0362285, -0.0618265),
                4 => (0.400212, 0.156643, -0.174875, -0.122959),
                5 => (0.362157, 0.149391, -0.305682, -0.143846),
                _ => (0.348844, 0.146829, -0.348844, -0.146829),
            };
            close(s1.re, re1, 1e-5);
            close(s1.im, im1, 1e-5);
            close(s2.re, re2, 1e-5);
            close(s2.im, im2, 1e-5);
        }
    }

    #[test]
    fn rayleigh_fold_converges_quadratically() {
        let m = C::new(1.05, 0.0);
        let mut err_prev = f64::INFINITY;
        // Error falls as O(x^2) as x -> 0.
        let mut ratios = Vec::new();
        for &x in &[0.1, 0.05, 0.025, 0.0125] {
            let exact = efficiencies(x, m).q_sca;
            let approx = rayleigh_qsca(x, m);
            let rel = (exact - approx).abs() / exact;
            ratios.push(rel);
            if err_prev.is_finite() {
                let _ = err_prev;
            }
            err_prev = rel;
        }
        for w in ratios.windows(2) {
            // halving x should cut the relative error by about 4.
            assert!(w[0] / w[1] > 3.0 && w[0] / w[1] < 5.0, "{:?}", ratios);
        }
    }

    #[test]
    fn anomalous_diffraction_fold_behaviour() {
        // For x >> 1 and |m-1| << 1 the ADT fold approaches the exact Q_ext.
        let m = 1.05;
        let x = 200.0;
        let exact = efficiencies(x, C::new(m, 0.0)).q_ext;
        let adt = anomalous_diffraction_qext_real(x, m);
        assert!((exact - adt).abs() < 0.08, "exact {exact} adt {adt}");
        // Extinction paradox: Q_ext -> 2 as x -> infinity.
        let big = efficiencies(100_000.0, C::new(1.5, 0.0)).q_ext;
        assert!((big - 2.0).abs() < 0.01, "Q_ext(1e5) = {big}");
    }

    #[test]
    fn small_x_limit_and_phase_function() {
        // Rayleigh phase function: p(0) = p(pi) = 3/2, p(pi/2) = 3/4.
        close(rayleigh_phase(0.0), 1.5, 1e-14);
        close(rayleigh_phase(PI), 1.5, 1e-14);
        close(rayleigh_phase(PI / 2.0), 0.75, 1e-14);
        // Mie phase function tends to the Rayleigh form for x << 1.
        let s = MieScattering::new(1e-3, C::new(1.5, 0.0));
        let q = s.efficiencies().q_sca;
        for &t in &[0.0, 0.7, PI / 2.0, PI] {
            close(s.phase_function_with(t, q), rayleigh_phase(t), 1e-4);
        }
    }

    #[test]
    fn gamma_distribution_moments() {
        let d = GammaSizeDistribution::stratocumulus();
        // r_eff = <r^3>/<r^2>.
        close(d.moment(3.0) / d.moment(2.0), d.r_eff, 1e-12 * d.r_eff);
        // v_eff = <r^4>/(r_eff^2 <r^2>) - 1.
        let v = d.moment(4.0) / (d.r_eff * d.r_eff * d.moment(2.0)) - 1.0;
        close(v, d.v_eff, 1e-12);
        // Simpson at 4000 points and the truncated tail agree with 1 to ~1e-8.
        let norm = simpson(
            |r| d.pdf(r),
            0.0,
            d.r_eff * (1.0 + 8.0 * d.v_eff.sqrt()),
            4000,
        );
        close(norm, 1.0, 1e-7);
    }

    #[test]
    fn cloud_extinction_matches_textbook_scaling() {
        // Stratocumulus: r_eff = 10 µm, v_eff = 0.1, N = 100 cm^-3.
        // The large-particle law sigma_ext = 3 LWC / (2 rho_w r_eff) should hold
        // with Q_ext -> 2; the textbook value at LWC = 0.3 g/m^3 is ~0.045 m^-1.
        let d = GammaSizeDistribution::stratocumulus();
        let n_number = 100e6; // 100 cm^-3
                              // Water index from the crate's Daimon-Masumura Sellmeier fit.
        let n_water = Sellmeier::water().index(550e-9).re;
        let m = C::new(n_water, 0.0);
        let mu = cloud_medium(&d, m, 550e-9, 1.0, n_number, 4000);
        // Mass from the distribution and the sphere density.
        let rho = 1000.0; // water, kg/m^3
        let mass_per_particle = rho * 4.0 / 3.0 * PI * d.moment(3.0);
        let lwc = n_number * mass_per_particle;
        let sigma_large = 3.0 * lwc / (2.0 * rho * d.r_eff);
        assert!(
            mu.sigma_ext > 0.02 && mu.sigma_ext < 0.09,
            "sigma_ext={}",
            mu.sigma_ext
        );
        // Extinction is close to the Q_ext -> 2 law (within 10%).
        assert!(
            (mu.sigma_ext - sigma_large).abs() / sigma_large < 0.1,
            "mie {} large {}",
            mu.sigma_ext,
            sigma_large
        );
        // Visible cloud droplets are almost conservative scatterers, g ~ 0.85.
        assert!(mu.albedo() > 0.99, "albedo={}", mu.albedo());
        assert!((mu.g - 0.85).abs() < 0.05, "g={}", mu.g);
    }

    #[test]
    fn optical_force_scales_and_traps() {
        // 1 µm silica bead in water at 1064 nm, 100 mW, w0 = 0.5 µm, 300 K.
        let radius = 0.5e-6;
        let n_silica = Sellmeier::fused_silica().index(1064e-9).re;
        let n_water = Sellmeier::water().index(1064e-9).re;
        let m = C::new(n_silica / n_water, 0.0);
        let f = optical_force(radius, 1064e-9, n_water, m, 0.1, 0.5e-6, 300.0);
        assert!(f.q_pr > 0.01 && f.q_pr < 0.5, "q_pr={}", f.q_pr);
        // Scattering force is a fraction of a pN at 100 mW.
        assert!(
            f.scattering_force > 1e-12 && f.scattering_force < 1e-9,
            "Fscat={}",
            f.scattering_force
        );
        // The trap is deep compared with k_B T.
        assert!(f.depth_over_kbt > 10.0, "U/kT = {}", f.depth_over_kbt);
        // Both forces scale linearly with power.
        let f2 = optical_force(radius, 1064e-9, n_water, m, 0.2, 0.5e-6, 300.0);
        close(f2.scattering_force / f.scattering_force, 2.0, 1e-12);
        close(f2.gradient_force / f.gradient_force, 2.0, 1e-12);
    }

    #[test]
    fn virion_cross_section_is_submicron_squared() {
        // 100 nm virion in water, index estimate m = 1.1 (protein/lipid envelope).
        let radius = 50e-9;
        let s = MieScattering::from_radius(radius, 550e-9, C::new(1.1 * 1.3334, 0.0), 1.3334);
        let e = s.efficiencies();
        let c_sca = e.q_sca * PI * radius * radius;
        // Weak scatterer: Q_sca small and the cross-section below the geometric one.
        assert!(e.q_sca < 0.1, "Q_sca={}", e.q_sca);
        assert!(c_sca > 1e-18 && c_sca < 1e-15, "C_sca={c_sca}");
        // Positive, and forward-scattering dominated.
        assert!(e.g > 0.0 && e.g < 1.0, "g={}", e.g);
        // Compare with the Rayleigh fold: within a factor few at x ~ 0.76.
        let rayleigh = rayleigh_qsca(s.x, s.m) * PI * radius * radius;
        assert!(c_sca / rayleigh > 0.5 && c_sca / rayleigh < 2.0);
    }
}
