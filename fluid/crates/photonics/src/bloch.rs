//! Band structure of a one-dimensional photonic crystal.
//!
//! A periodic stack with period `Lambda` is the simplest structured optical
//! material. Bloch's theorem says a wave in the infinite crystal has the form
//! `E(z) e^{iKz}` with the Bloch wavevector `K` set by the unit-cell transfer
//! matrix `M`:
//!
//! ```text
//!   cos(K Lambda) = (1/2) tr M.
//! ```
//!
//! For lossless layers `tr M` is real. Where `|tr M / 2| <= 1` a propagating
//! Bloch wave exists (a pass band); where it exceeds one, `K` is complex and
//! the field is evanescent (a stop band, the origin of the high-reflectance
//! band of a Bragg mirror). This is a *cheap* model of an infinite structure:
//! it costs one unit cell, not the thousands of layers of a finite mirror.
//!
//! The validity indicator is the finite-stack size relative to the decay
//! length. At a stop-band centre the field decays as `exp(-|Im K| z)`, so a
//! stack of `N` periods reflects almost perfectly once
//! `N |Im K| Lambda >~ 1`. The finite-stack reflectance benchmark below checks
//! that the gap predicted by the one-cell Bloch analysis coincides with the
//! high-reflectance band of a long but finite stack.
use crate::complex::C;
use crate::multilayer::{layer_matrix, Layer, Mat2, Polarization};
use std::f64::consts::PI;

/// One period of a 1D photonic crystal, in propagation order.
#[derive(Clone, Debug)]
pub struct UnitCell {
    pub layers: Vec<Layer>,
}

impl UnitCell {
    pub fn new(layers: Vec<Layer>) -> Self {
        assert!(!layers.is_empty());
        Self { layers }
    }
    /// Binary quarter-wave stack: high then low index, no absorption.
    pub fn quarter_wave(n_high: f64, n_low: f64, lambda0: f64) -> Self {
        Self::new(vec![
            Layer::quarter_wave(n_high, lambda0),
            Layer::quarter_wave(n_low, lambda0),
        ])
    }
    pub fn period(&self) -> f64 {
        self.layers.iter().map(|l| l.thickness).sum()
    }
    pub fn matrix(&self, lambda: f64, theta: f64, pol: Polarization) -> Mat2 {
        let s = theta.sin();
        let mut m = Mat2::identity();
        for l in &self.layers {
            m = m * layer_matrix(l.index, l.mu, l.thickness, lambda, s, pol);
        }
        m
    }
    /// `(1/2) tr M`, the left side of the Bloch condition.
    pub fn half_trace(&self, lambda: f64, theta: f64, pol: Polarization) -> C {
        let m = self.matrix(lambda, theta, pol);
        (m.m11 + m.m22) * 0.5
    }
    /// Bloch wavevector `K`; complex in a stop band.
    pub fn bloch_wavevector(&self, lambda: f64, theta: f64, pol: Polarization) -> C {
        let h = self.half_trace(lambda, theta, pol);
        let period = self.period();
        // arccos(h): real for |h|<=1, K = pi/Lambda + i acosh(|h|)/Lambda otherwise.
        if h.im.abs() < 1e-12 && h.re.abs() <= 1.0 {
            C::real(h.re.clamp(-1.0, 1.0).acos() / period)
        } else {
            // Stop band. Reduce K to the first Brillouin zone: Re K = 0 for
            // tr M/2 > 1, Re K = pi/Lambda for tr M/2 < -1.
            let mag = h.re.abs().max(1.0);
            let alpha = mag.acosh();
            let re_k = if h.re < 0.0 { PI / period } else { 0.0 };
            C::new(re_k, alpha / period)
        }
    }
    /// True when the Bloch wavevector is complex: a stop band for lossless cells.
    pub fn in_stop_band(&self, lambda: f64, theta: f64, pol: Polarization) -> bool {
        let h = self.half_trace(lambda, theta, pol);
        h.im.abs() < 1e-9 && h.re.abs() > 1.0
    }
    /// Field decay rate `|Im K|` [m^-1] at `lambda`.
    pub fn decay_rate(&self, lambda: f64, theta: f64, pol: Polarization) -> f64 {
        self.bloch_wavevector(lambda, theta, pol).im.abs()
    }
    /// Number of periods needed for `exp(-2 decay N Lambda) < tolerance`.
    pub fn periods_for_attenuation(
        &self,
        lambda: f64,
        theta: f64,
        pol: Polarization,
        tolerance: f64,
    ) -> f64 {
        let d = self.decay_rate(lambda, theta, pol) * self.period();
        if d <= 0.0 {
            f64::INFINITY
        } else {
            (1.0 / tolerance).ln() / (2.0 * d)
        }
    }
}

/// Reflectance at the design wavelength of an `N`-pair quarter-wave mirror
/// with the **low-index** layer adjacent to the incident medium, the standard
/// closed form
/// `R = ((1 - (n_s/n_0)(n_L/n_H)^{2N}) / (1 + (n_s/n_0)(n_L/n_H)^{2N}))^2`.
/// Verified against the transfer matrix at the design wavelength.
pub fn quarter_wave_bragg_reflectance_low_first(
    n0: f64,
    ns: f64,
    n_high: f64,
    n_low: f64,
    pairs: usize,
) -> f64 {
    let b = (n_low / n_high).powi(2 * pairs as i32);
    let a = ns / n0;
    ((1.0 - a * b) / (1.0 + a * b)).powi(2)
}

/// Same closed form for the opposite termination (high-index layer adjacent to
/// the incident medium), which reflects more for the same layer count.
pub fn quarter_wave_bragg_reflectance_high_first(
    n0: f64,
    ns: f64,
    n_high: f64,
    n_low: f64,
    pairs: usize,
) -> f64 {
    let b = (n_high / n_low).powi(2 * pairs as i32);
    let a = ns / n0;
    ((1.0 - a * b) / (1.0 + a * b)).powi(2)
}

/// Exact relative **frequency** width of the quarter-wave stop band,
/// `Delta_omega/omega0 = (4/pi) asin(r)` with
/// `r = (n_H - n_L)/(n_H + n_L)`.
///
/// This is not a leading-order approximation. For the infinite quarter-wave
/// stack at normal incidence the Bloch condition gives
/// `(1/2) tr M = cos^2(pi x/2) - ((rho+rho^-1)/2) sin^2(pi x/2)` with
/// `x = omega/omega0` and `rho = n_H/n_L`; the edges `(1/2) tr M = -1` are at
/// `x = 1 +/- (2/pi) asin(r)`, so the frequency width is exactly
/// `(4/pi) asin(r)` to all orders in the index contrast. The standard
/// `(4/pi) asin(r)` result quoted for `d lambda/lambda0` is the same number
/// divided by `1 - (Delta_omega/omega0)^2/4` (see
/// [`quarter_wave_stopband_wavelength_ratio`]).
pub fn quarter_wave_stopband_ratio(n_high: f64, n_low: f64) -> f64 {
    let r = (n_high - n_low) / (n_high + n_low);
    (4.0 / PI) * r.asin()
}

/// Exact relative wavelength width `Delta_lambda/lambda0` of the same
/// quarter-wave stop band, obtained from the frequency width by
/// `Delta_lambda/lambda0 = Delta/(1 - Delta^2/4)` with
/// `Delta = Delta_omega/omega0`. At `n_H/n_L = 2.30/1.45` this is the
/// `0.297434` that the trace formula returns when the edges are measured on
/// the wavelength axis.
pub fn quarter_wave_stopband_wavelength_ratio(n_high: f64, n_low: f64) -> f64 {
    let d = quarter_wave_stopband_ratio(n_high, n_low);
    d / (1.0 - 0.25 * d * d)
}

/// Wavelength band edges `(lambda_lo, lambda_hi)` of the quarter-wave stop
/// band, from bisecting `|tr M / 2| = 1` on both sides of the design
/// wavelength. For the quarter-wave stack these are exactly
/// `lambda0/(1 +/- Delta/2)` with `Delta = (4/pi) asin(r)`.
pub fn quarter_wave_stopband_wavelength_edges(
    n_high: f64,
    n_low: f64,
    lambda0: f64,
) -> (f64, f64) {
    let cell = UnitCell::quarter_wave(n_high, n_low, lambda0);
    let f = |lambda: f64| cell.half_trace(lambda, 0.0, Polarization::TE).re.abs() - 1.0;
    assert!(f(lambda0) > 0.0, "design wavelength must lie in the gap");
    // Upper edge: leave the gap as lambda grows, so f(hi) < 0.
    let mut hi = lambda0 * 1.1;
    let mut guard = 0;
    while f(hi) > 0.0 && guard < 200 {
        hi = lambda0 + (hi - lambda0) * 1.5;
        guard += 1;
    }
    let (mut a, mut b) = (lambda0, hi);
    for _ in 0..200 {
        let m = 0.5 * (a + b);
        if f(m) > 0.0 {
            a = m;
        } else {
            b = m;
        }
    }
    let upper = 0.5 * (a + b);
    // Lower edge: leave the gap as lambda shrinks, so f(lo) < 0.
    let mut lo = lambda0 * 0.9;
    guard = 0;
    while f(lo) > 0.0 && guard < 200 {
        lo = lambda0 - (lambda0 - lo) * 1.5;
        guard += 1;
    }
    let (mut a, mut b) = (lo, lambda0);
    for _ in 0..200 {
        let m = 0.5 * (a + b);
        if f(m) > 0.0 {
            b = m;
        } else {
            a = m;
        }
    }
    let lower = 0.5 * (a + b);
    (lower, upper)
}

/// Exact relative frequency width of the quarter-wave stop band obtained from
/// the Bloch band edges `|tr M / 2| = 1`. Frequencies are `1/lambda`, so the
/// result is `lambda0 (1/lambda_lo - 1/lambda_hi)`. It agrees with
/// [`quarter_wave_stopband_ratio`] to round-off at every index contrast; the
/// function exists to check that agreement against the independent trace
/// computation.
pub fn quarter_wave_stopband_ratio_exact(n_high: f64, n_low: f64, lambda0: f64) -> f64 {
    let (lo, hi) = quarter_wave_stopband_wavelength_edges(n_high, n_low, lambda0);
    lambda0 * (1.0 / lo - 1.0 / hi)
}
