//! Transfer matrices for 1D elastic/acoustic layered media, and the Bloch
//! dispersion of the infinite periodic stack.
//!
//! Longitudinal waves at normal incidence in a layer are described by a
//! pressure `p` and a particle velocity `v` that are continuous across every
//! interface. For a homogeneous layer of thickness `d`, density `rho` and
//! complex modulus `M` the impedance is `Z = rho c = sqrt(rho M)`, the
//! wavenumber is `k = omega / c = omega sqrt(rho / M)`, and the pair at the far
//! face follows from the pair at the near face through the exact characteristic
//! matrix
//!
//! ```text
//!   [ p(d) ]     [  cos(k d)   -i Z sin(k d) ] [ p(0) ]
//!   [      ]  =  [                            ] [      ].
//!   [ v(d) ]     [ -i sin(k d)/Z   cos(k d)   ] [ v(0) ]
//! ```
//!
//! The product over the layers gives the response of a finite stack between two
//! half-spaces, and `r`, `t`, `R`, `T`, `A` follow. The same unit-cell matrix
//! gives the exact band structure of the infinite periodic medium through the
//! Rytov relation,
//!
//! ```text
//!   cos(q Lambda) = cos(k1 d1) cos(k2 d2)
//!                 - (1/2)(Z1/Z2 + Z2/Z1) sin(k1 d1) sin(k2 d2).
//! ```
//!
//! # Exactness and the emergent boundaries
//!
//! The transfer matrix is exact for planar, isotropic, perfectly bonded
//! layers. It resolves the interference of multiply reflected partial waves,
//! which is a **spatial** sub-grid structure: when the period is much smaller
//! than the wavelength the stack must be folded *up* to the effective
//! medium ([`crate::homogenize`]), and when a layer becomes much thicker than
//! the pulse it must be unfolded *down* to a ray/intensity description. The
//! cheap indicators are the period-to-wavelength ratio and the number of
//! periods times the per-period Bloch decay.
use crate::complex::C;
use crate::materials::{Fluid, Solid};
use std::f64::consts::PI;

/// One homogeneous layer. `modulus` is the relevant elastic modulus at normal
/// incidence: the bulk modulus for a fluid, or the constrained (P-wave)
/// modulus `lambda + 2 mu` for a solid.
#[derive(Clone, Copy, Debug)]
pub struct Layer {
    pub rho: f64,
    /// Complex modulus `M = M0 (1 + i eta)`; `eta >= 0` is absorption.
    pub modulus: C,
    pub thickness: f64,
}

impl Layer {
    /// Lossless layer from density and phase speed `c`.
    pub fn lossless(rho: f64, c: f64, thickness: f64) -> Self {
        assert!(rho > 0.0 && c > 0.0 && thickness > 0.0);
        Self {
            rho,
            modulus: C::real(rho * c * c),
            thickness,
        }
    }
    /// Layer from an explicit (possibly complex) modulus.
    pub fn from_modulus(rho: f64, modulus: C, thickness: f64) -> Self {
        assert!(rho > 0.0 && thickness > 0.0 && modulus.re > 0.0);
        Self {
            rho,
            modulus,
            thickness,
        }
    }
    /// Lossy layer, `M = rho c^2 (1 + i loss)`.
    pub fn with_loss(rho: f64, c: f64, thickness: f64, loss: f64) -> Self {
        assert!(loss >= 0.0);
        let m = rho * c * c;
        Self::from_modulus(rho, C::new(m, m * loss), thickness)
    }
    /// Layer from a fluid or a solid (P-wave, normal incidence).
    pub fn from_fluid(f: &Fluid, thickness: f64) -> Self {
        Self::lossless(f.rho, f.c, thickness)
    }
    /// Normal-incidence layer of a solid: uses `rho c_p^2 = lambda + 2 mu`.
    pub fn from_solid(s: &Solid, thickness: f64) -> Self {
        Self::from_modulus(s.rho, C::real(s.rho * s.c_p * s.c_p), thickness)
    }
    /// Acoustic impedance `Z = sqrt(rho M)`.
    pub fn impedance(&self) -> C {
        (self.modulus * self.rho).sqrt()
    }
    /// Wavenumber `k = omega sqrt(rho / M)`, principal branch so a passive
    /// modulus gives `Im(k) >= 0`.
    pub fn wavenumber(&self, omega: f64) -> C {
        (C::real(self.rho) / self.modulus).sqrt() * omega
    }
    /// Phase speed `c = sqrt(M / rho)`.
    pub fn speed(&self) -> C {
        (self.modulus / self.rho).sqrt()
    }
    /// Quarter-wave thickness at design frequency `f0` (lossless convention).
    pub fn quarter_wave(rho: f64, c: f64, f0: f64) -> Self {
        Self::lossless(rho, c, c / (4.0 * f0))
    }
}

/// 2x2 characteristic matrix, row-major.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Mat2 {
    pub m11: C,
    pub m12: C,
    pub m21: C,
    pub m22: C,
}

impl Mat2 {
    pub fn identity() -> Self {
        Self {
            m11: C::real(1.0),
            m12: C::real(0.0),
            m21: C::real(0.0),
            m22: C::real(1.0),
        }
    }
    /// Characteristic matrix from impedance, wavenumber and thickness.
    pub fn layer(z: C, k: C, thickness: f64) -> Self {
        let theta = k * thickness;
        let (s, c) = (theta.sin_complex(), theta.cos_complex());
        let i = C::i();
        Self {
            m11: c,
            m12: -(i * z * s),
            m21: -(i * s / z),
            m22: c,
        }
    }
    /// Characteristic matrix of a layer at angular frequency `omega`.
    pub fn from_layer(l: &Layer, omega: f64) -> Self {
        Self::layer(l.impedance(), l.wavenumber(omega), l.thickness)
    }
    pub fn apply(self, v: (C, C)) -> (C, C) {
        (
            self.m11 * v.0 + self.m12 * v.1,
            self.m21 * v.0 + self.m22 * v.1,
        )
    }
}

impl std::ops::Mul for Mat2 {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self {
            m11: self.m11 * b.m11 + self.m12 * b.m21,
            m12: self.m11 * b.m12 + self.m12 * b.m22,
            m21: self.m21 * b.m11 + self.m22 * b.m21,
            m22: self.m21 * b.m12 + self.m22 * b.m22,
        }
    }
}

/// Helpers on `C` for trigonometric functions of a complex argument. Kept here
/// so `complex.rs` stays a pure arithmetic type.
trait Trig {
    fn cos_complex(self) -> C;
    fn sin_complex(self) -> C;
}
impl Trig for C {
    fn cos_complex(self) -> C {
        let e = (C::i() * self).exp();
        let em = (C::i() * (-self)).exp();
        (e + em) * 0.5
    }
    fn sin_complex(self) -> C {
        let e = (C::i() * self).exp();
        let em = (C::i() * (-self)).exp();
        (e - em) / (C::i() * 2.0)
    }
}

/// A finite stack between two half-spaces. Half-spaces are specified by their
/// acoustic impedance `Z` (real for a lossless fluid or at normal incidence).
#[derive(Clone, Debug)]
pub struct Stack {
    pub incident_z: C,
    pub layers: Vec<Layer>,
    pub exit_z: C,
}

/// Solution of a finite stack at one frequency.
#[derive(Clone, Copy, Debug)]
pub struct Solution {
    /// Pressure reflection coefficient.
    pub r: C,
    /// Pressure transmission coefficient.
    pub t: C,
    pub reflectance: f64,
    /// Power transmittance `(Z_in/Z_out) |t|^2`.
    pub transmittance: f64,
    /// `1 - R - T`; positive for a lossy stack.
    pub absorptance: f64,
    pub incident_z: C,
    pub exit_z: C,
    pub omega: f64,
}

impl Solution {
    /// Particle-velocity transmission amplitude `t_v = t Z_in / Z_out`. With
    /// this convention the lossless energy identity takes the symmetric form
    /// `|r|^2 + (Z_out/Z_in) |t_v|^2 = 1`.
    pub fn velocity_transmission(&self) -> C {
        self.t * self.incident_z / self.exit_z
    }
    /// `R + T`, equal to one for a lossless stack.
    pub fn energy(&self) -> f64 {
        self.reflectance + self.transmittance
    }
}

impl Stack {
    pub fn new(incident_z: f64, layers: Vec<Layer>, exit_z: f64) -> Self {
        assert!(incident_z > 0.0 && exit_z > 0.0);
        Self {
            incident_z: C::real(incident_z),
            layers,
            exit_z: C::real(exit_z),
        }
    }
    pub fn from_media(incident: &Fluid, layers: Vec<Layer>, exit: &Fluid) -> Self {
        Self::new(incident.impedance(), layers, exit.impedance())
    }
    pub fn single(incident_z: f64, layer: Layer, exit_z: f64) -> Self {
        Self::new(incident_z, vec![layer], exit_z)
    }
    /// Total characteristic matrix of the finite layers.
    pub fn matrix(&self, omega: f64) -> Mat2 {
        let mut m = Mat2::identity();
        for l in &self.layers {
            m = m * Mat2::from_layer(l, omega);
        }
        m
    }
    /// Solve the stack at angular frequency `omega`.
    pub fn solve(&self, omega: f64) -> Solution {
        assert!(omega > 0.0 && omega.is_finite());
        let m = self.matrix(omega);
        let z0 = self.incident_z;
        let zs = self.exit_z;
        // Input impedance of the terminated stack.
        let z_in = (m.m22 - m.m12 / zs) / (C::real(-1.0) * m.m21 + m.m11 / zs);
        let r = (z_in - z0) / (z_in + z0);
        let t = (C::real(1.0) + r) * (m.m11 + m.m12 / z_in);
        let reflectance = r.abs2();
        let transmittance = (z0.re / zs.re) * t.abs2();
        Solution {
            r,
            t,
            reflectance,
            transmittance,
            absorptance: 1.0 - reflectance - transmittance,
            incident_z: z0,
            exit_z: zs,
            omega,
        }
    }
    pub fn reflectance(&self, omega: f64) -> f64 {
        self.solve(omega).reflectance
    }
}

/// One period of an infinite 1D periodic medium: two lossless (or lossy)
/// layers, in propagation order.
#[derive(Clone, Copy, Debug)]
pub struct Periodic {
    pub a: Layer,
    pub b: Layer,
}

impl Periodic {
    pub fn new(a: Layer, b: Layer) -> Self {
        Self { a, b }
    }
    /// Binary quarter-wave stack at design frequency `f0`.
    pub fn quarter_wave(rho1: f64, c1: f64, rho2: f64, c2: f64, f0: f64) -> Self {
        Self::new(
            Layer::quarter_wave(rho1, c1, f0),
            Layer::quarter_wave(rho2, c2, f0),
        )
    }
    pub fn period(&self) -> f64 {
        self.a.thickness + self.b.thickness
    }
    /// Left side of the Rytov relation,
    /// `cos(k1 d1) cos(k2 d2) - (1/2)(Z1/Z2 + Z2/Z1) sin(k1 d1) sin(k2 d2)`.
    pub fn rytov_cos(&self, omega: f64) -> C {
        let z1 = self.a.impedance();
        let z2 = self.b.impedance();
        let t1 = self.a.wavenumber(omega) * self.a.thickness;
        let t2 = self.b.wavenumber(omega) * self.b.thickness;
        t1.cos_complex() * t2.cos_complex()
            - (z1 / z2 + z2 / z1) * 0.5 * t1.sin_complex() * t2.sin_complex()
    }
    /// Unit-cell transfer matrix (same as the Rytov left side is `tr M / 2`).
    pub fn cell_matrix(&self, omega: f64) -> Mat2 {
        Mat2::from_layer(&self.a, omega) * Mat2::from_layer(&self.b, omega)
    }
    /// Exact Bloch wavevector `q` of the infinite periodic medium. Real in a
    /// pass band, complex in a band gap (reduced to `Re q in [0, pi/Lambda]`).
    pub fn bloch_wavevector(&self, omega: f64) -> C {
        let h = self.rytov_cos(omega);
        let period = self.period();
        if h.im.abs() < 1e-12 && h.re.abs() <= 1.0 {
            C::real(h.re.clamp(-1.0, 1.0).acos() / period)
        } else {
            let mag = h.abs().max(1.0);
            let alpha = mag.acosh();
            let re_q = if h.re < 0.0 { PI / period } else { 0.0 };
            C::new(re_q, alpha / period)
        }
    }
    pub fn in_band_gap(&self, omega: f64) -> bool {
        let h = self.rytov_cos(omega);
        h.im.abs() < 1e-9 && h.re.abs() > 1.0
    }
    /// Field decay rate `|Im q|` [m^-1].
    pub fn decay_rate(&self, omega: f64) -> f64 {
        self.bloch_wavevector(omega).im.abs()
    }
    /// Number of periods for `exp(-2 decay N Lambda) < tolerance`.
    pub fn periods_for_attenuation(&self, omega: f64, tolerance: f64) -> f64 {
        let d = self.decay_rate(omega) * self.period();
        if d <= 0.0 {
            f64::INFINITY
        } else {
            (1.0 / tolerance).ln() / (2.0 * d)
        }
    }
}

/// Exact normal-incidence band edges of a lossless quarter-wave stack, in
/// frequency. The gap is centred on `f0` with relative width
/// `(4/pi) asin((Z1 - Z2)/(Z1 + Z2))`; the exact edges solve
/// `tr M / 2 = -1`.
pub fn quarter_wave_band_edges(z1: f64, z2: f64, f0: f64) -> (f64, f64) {
    assert!(z1 > 0.0 && z2 > 0.0 && f0 > 0.0);
    let asin = ((z1 - z2) / (z1 + z2)).asin();
    let half = (2.0 / PI) * asin;
    (f0 * (1.0 - half), f0 * (1.0 + half))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "{a:e} vs {b:e}");
    }

    #[test]
    fn single_interface_fresnel_energy() {
        // Two half-spaces, no layer: r = (Z2 - Z1)/(Z2 + Z1).
        let stack = Stack::new(1.0, vec![], 2.0);
        let s = stack.solve(1.0);
        close(s.r.re, 1.0 / 3.0, 1e-14);
        close(s.energy(), 1.0, 1e-14);
        let tv = s.velocity_transmission();
        close(s.reflectance + (s.exit_z.re / s.incident_z.re) * tv.abs2(), 1.0, 1e-14);
    }
}
