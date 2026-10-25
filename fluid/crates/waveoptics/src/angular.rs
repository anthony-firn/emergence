//! Band-limited angular-spectrum propagation of the openings' fields to a
//! plane: the fast rung for whole planes of points (a device's entrance
//! pupil).
//!
//! The first Rayleigh–Sommerfeld integral with given boundary values is
//! *exactly* the angular-spectrum propagation of those values (Goodman,
//! *Introduction to Fourier Optics*, §3.10):
//!
//! ```text
//! U(u, z) = ∫∫ Ũ(f) e^{i k z √(1 − λ²|f|²)} e^{i 2π f·u} d²f,
//! Ũ(f) = ∫∫_opening t(ξ) e^{−i 2π f·ξ} d²ξ,
//! ```
//!
//! with `t` the Kirchhoff boundary value (the incident wave on the
//! opening). For a rectangle and a plane incident wave `Ũ` is a product of
//! closed forms (exact); for a point source the incident phase is expanded
//! about the opening centre to second order in each axis (closed-form
//! Fresnel moments), and the dropped cross and cubic terms and the `1/R`
//! amplitude variation are bounded ([`OpeningSpectrum::phase_bound`]).
//!
//! A device only records the directions that reach its detector window, so
//! the integral is restricted to a band of spatial frequencies (with a
//! raised-cosine roll-off) and evaluated as a separable sum on a frequency
//! grid whose spacing makes the periodic replicas of the band-limited field
//! fall clear of the evaluation region (Poisson summation). The power of
//! each opening outside the band is reported, so the truncation is visible.

use std::f64::consts::PI;

use num_complex::Complex64 as Cx;
use rayon::prelude::*;

use crate::diffraction::{moments, Incident, Rect};

/// The angular spectrum of one opening's boundary values.
#[derive(Clone, Copy, Debug)]
pub struct OpeningSpectrum {
    center: [f64; 2],
    half: [f64; 2],
    /// Incident field at the opening centre.
    tc: Cx,
    /// Phase gradient and second derivatives of the incident wave at the
    /// centre.
    g: [f64; 2],
    h: [f64; 2],
    /// Bound on the dropped incident phase (rad) and relative amplitude
    /// variation over the opening.
    pub phase_bound: f64,
}

impl OpeningSpectrum {
    pub fn new(rect: &Rect, incident: &Incident, wavelength: f64) -> Self {
        let k = 2.0 * PI / wavelength;
        let c = rect.center;
        let half = [0.5 * rect.size[0], 0.5 * rect.size[1]];
        match *incident {
            Incident::Plane { amplitude, sines } => Self {
                center: c,
                half,
                tc: amplitude * Cx::from_polar(1.0, k * (sines[0] * c[0] + sines[1] * c[1])),
                g: [k * sines[0], k * sines[1]],
                h: [0.0; 2],
                phase_bound: 0.0,
            },
            Incident::Point { amplitude, position, distance } => {
                let w = [c[0] - position[0], c[1] - position[1]];
                let r = (w[0] * w[0] + w[1] * w[1] + distance * distance).sqrt();
                let hxy = k * w[0] * w[1] / r.powi(3);
                // third derivative of R along any transverse direction
                // ≤ min(3ρ/R³, 1.155/R²) (as in the diffraction module)
                let diag = half[0].hypot(half[1]);
                let rho = (w[0].abs() + half[0]).hypot(w[1].abs() + half[1]);
                let rmin = ((w[0].abs() - half[0]).max(0.0).powi(2) + (w[1].abs() - half[1]).max(0.0).powi(2) + distance * distance).sqrt();
                let d3 = (3.0 * rho / rmin.powi(3)).min(1.155 / (rmin * rmin));
                let amp = rho / (rmin * rmin) * diag; // |∇ ln(1/R)| |Δ|
                Self {
                    center: c,
                    half,
                    tc: amplitude * Cx::from_polar(1.0 / r, k * r),
                    g: [k * w[0] / r, k * w[1] / r],
                    h: [k * (1.0 / r - w[0] * w[0] / r.powi(3)), k * (1.0 / r - w[1] * w[1] / r.powi(3))],
                    phase_bound: hxy.abs() * half[0] * half[1] + k / 6.0 * d3 * diag.powi(3) + amp,
                }
            }
        }
    }

    /// `Ũ(f)` (spatial frequency `f` in cycles/m).
    pub fn at(&self, f: [f64; 2]) -> Cx {
        self.tc * Cx::from_polar(1.0, -2.0 * PI * (f[0] * self.center[0] + f[1] * self.center[1])) * self.axis(0, f[0]) * self.axis(1, f[1])
    }

    /// One axis factor `∫ e^{i((g − 2πf)s + h s²/2)} ds` over the half-width.
    pub fn axis(&self, a: usize, f: f64) -> Cx {
        moments(self.g[a] - 2.0 * PI * f, self.h[a], self.half[a], true).0
    }

    /// Power through the opening, `∫|t|²` (to the accuracy of the local
    /// model).
    pub fn power(&self) -> f64 {
        self.tc.norm_sqr() * 4.0 * self.half[0] * self.half[1]
    }

    /// Transverse support of the opening along an axis.
    pub fn extent(&self, a: usize) -> (f64, f64) {
        (self.center[a] - self.half[a], self.center[a] + self.half[a])
    }
}

/// A rectangular band of spatial frequencies with a raised-cosine edge:
/// weight 1 within `inner` of `center`, falling to 0 at `outer`.
#[derive(Clone, Copy, Debug)]
pub struct Band {
    pub center: [f64; 2],
    pub inner: [f64; 2],
    pub outer: [f64; 2],
}

impl Band {
    pub fn weight(&self, a: usize, f: f64) -> f64 {
        let d = (f - self.center[a]).abs();
        if d <= self.inner[a] {
            1.0
        } else if d >= self.outer[a] {
            0.0
        } else {
            0.5 * (1.0 + (PI * (d - self.inner[a]) / (self.outer[a] - self.inner[a])).cos())
        }
    }
}

/// Fields of each opening on the grid `xs × ys` of the plane `z`:
/// `out[opening][i + xs.len() j]`, plus the fraction of each opening's
/// power outside the band. `period[a]` is the replica spacing (≥ the span
/// of the band-limited field together with the grid).
pub fn plane_fields(spectra: &[OpeningSpectrum], wavelength: f64, z: f64, band: &Band, period: [f64; 2], xs: &[f64], ys: &[f64]) -> (Vec<Vec<Cx>>, Vec<f64>) {
    let k = 2.0 * PI / wavelength;
    let freqs = |a: usize| -> Vec<f64> {
        let df = 1.0 / period[a];
        let m = (band.outer[a] / df).ceil() as i64;
        (-m..=m).map(|i| band.center[a] + i as f64 * df).collect()
    };
    let (fx, fy) = (freqs(0), freqs(1));
    let (mx, my) = (fx.len(), fy.len());
    let dfa = 1.0 / (period[0] * period[1]);
    let wx: Vec<f64> = fx.iter().map(|&f| band.weight(0, f)).collect();
    let wy: Vec<f64> = fy.iter().map(|&f| band.weight(1, f)).collect();
    // transfer function (evanescent components dropped: z ≫ λ)
    let hz: Vec<Cx> = (0..mx * my)
        .into_par_iter()
        .map(|q| {
            let (i, j) = (q % mx, q / mx);
            let s2 = wavelength * wavelength * (fx[i] * fx[i] + fy[j] * fy[j]);
            if s2 < 1.0 {
                Cx::from_polar(wx[i] * wy[j], k * z * (1.0 - s2).sqrt())
            } else {
                Cx::new(0.0, 0.0)
            }
        })
        .collect();
    let ex: Vec<Vec<Cx>> = xs.iter().map(|&x| fx.iter().map(|&f| Cx::from_polar(1.0, 2.0 * PI * f * x)).collect()).collect();
    let ey: Vec<Vec<Cx>> = ys.iter().map(|&y| fy.iter().map(|&f| Cx::from_polar(1.0, 2.0 * PI * f * y)).collect()).collect();
    let mut fields = Vec::with_capacity(spectra.len());
    let mut losses = Vec::with_capacity(spectra.len());
    for s in spectra {
        // separable spectrum: Ũ = tc e^{−i2πf·c} X(fx) Y(fy)
        let ax: Vec<Cx> = fx.iter().map(|&f| s.axis(0, f) * Cx::from_polar(1.0, -2.0 * PI * f * s.center[0])).collect();
        let ay: Vec<Cx> = fy.iter().map(|&f| s.axis(1, f) * Cx::from_polar(1.0, -2.0 * PI * f * s.center[1])).collect();
        let inband: f64 = (0..my).map(|j| (0..mx).map(|i| (ax[i] * ay[j] * s.tc).norm_sqr() * (wx[i] * wy[j]).powi(2)).sum::<f64>()).sum::<f64>() * dfa;
        losses.push((1.0 - inband / s.power()).max(0.0));
        // stage 1: A[j][x] = Σ_i Ũ H e^{i2π fx x}
        let a: Vec<Vec<Cx>> = (0..my)
            .into_par_iter()
            .map(|j| {
                let row: Vec<Cx> = (0..mx).map(|i| ax[i] * ay[j] * hz[i + mx * j]).collect();
                ex.iter().map(|e| row.iter().zip(e).map(|(r, e)| r * e).sum()).collect()
            })
            .collect();
        let pref = s.tc * dfa;
        let u: Vec<Cx> = (0..xs.len() * ys.len())
            .into_par_iter()
            .map(|q| {
                let (xi, yi) = (q % xs.len(), q / xs.len());
                pref * (0..my).map(|j| a[j][xi] * ey[yi][j]).sum::<Cx>()
            })
            .collect();
        fields.push(u);
    }
    (fields, losses)
}
