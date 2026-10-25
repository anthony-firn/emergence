//! Minimal complex arithmetic for wave optics.
//!
//! Wave optics is complex-valued for two independent reasons, and the crate
//! needs both. A passive material has a complex refractive index
//! `N = n + i kappa`; the imaginary part is absorption. And even in a lossless
//! stack a wave above the critical angle, or inside a photonic band gap, has an
//! imaginary normal wavevector: the field is evanescent, not zero.
//!
//! This is a deliberately small type (about a hundred lines) so the physics
//! modules have no external dependency. The branch of `sqrt` is the physically
//! relevant one: for `N^2 - s^2` with a positive imaginary part (an absorbing
//! medium) or a negative real part (an evanescent wave), the principal square
//! root has a non-negative imaginary part, giving decay along the propagation
//! direction. See `multilayer::normal_wavevector`.
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct C {
    pub re: f64,
    pub im: f64,
}

impl C {
    pub const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    pub const fn real(re: f64) -> Self {
        Self { re, im: 0.0 }
    }
    pub const fn i() -> Self {
        Self { re: 0.0, im: 1.0 }
    }
    pub fn is_finite(self) -> bool {
        self.re.is_finite() && self.im.is_finite()
    }
    pub fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }
    pub fn abs2(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }
    pub fn recip(self) -> Self {
        let d = self.abs2();
        Self::new(self.re / d, -self.im / d)
    }
    /// Principal square root (branch cut on the negative real axis, upper half).
    pub fn sqrt(self) -> Self {
        let r = self.abs();
        if r == 0.0 {
            return Self::default();
        }
        let re = (0.5 * (r + self.re)).sqrt();
        let im = (0.5 * (r - self.re)).sqrt();
        Self::new(re, if self.im < 0.0 { -im } else { im })
    }
    pub fn exp(self) -> Self {
        let e = self.re.exp();
        Self::new(e * self.im.cos(), e * self.im.sin())
    }
    pub fn powi(self, n: i32) -> Self {
        let mut out = Self::real(1.0);
        let mut b = self;
        let mut k = n.unsigned_abs();
        while k > 0 {
            if k & 1 == 1 {
                out = out * b;
            }
            b = b * b;
            k >>= 1;
        }
        if n < 0 {
            out.recip()
        } else {
            out
        }
    }
    pub fn scale(self, s: f64) -> Self {
        Self::new(self.re * s, self.im * s)
    }
}

impl Add for C {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self::new(self.re + b.re, self.im + b.im)
    }
}
impl Sub for C {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self::new(self.re - b.re, self.im - b.im)
    }
}
impl Mul for C {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self::new(
            self.re * b.re - self.im * b.im,
            self.re * b.im + self.im * b.re,
        )
    }
}
impl Div for C {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        let d = b.abs2();
        Self::new(
            (self.re * b.re + self.im * b.im) / d,
            (self.im * b.re - self.re * b.im) / d,
        )
    }
}
impl Neg for C {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.re, -self.im)
    }
}
impl Mul<f64> for C {
    type Output = Self;
    fn mul(self, b: f64) -> Self {
        self.scale(b)
    }
}
impl Div<f64> for C {
    type Output = Self;
    fn div(self, b: f64) -> Self {
        self.scale(1.0 / b)
    }
}
impl Add<f64> for C {
    type Output = Self;
    fn add(self, b: f64) -> Self {
        Self::new(self.re + b, self.im)
    }
}
impl Sub<f64> for C {
    type Output = Self;
    fn sub(self, b: f64) -> Self {
        Self::new(self.re - b, self.im)
    }
}
impl From<f64> for C {
    fn from(v: f64) -> Self {
        Self::real(v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: C, b: C, tol: f64) {
        assert!((a - b).abs() <= tol, "{a:?} vs {b:?}");
    }

    #[test]
    fn arithmetic_and_principal_sqrt() {
        let z = C::new(3.0, 4.0);
        close(z * z.recip(), C::real(1.0), 1e-15);
        close(z.sqrt() * z.sqrt(), z, 1e-14);
        // Negative real axis: principal branch is +i sqrt(|x|).
        close(C::real(-4.0).sqrt(), C::new(0.0, 2.0), 1e-15);
        // Absorbing medium: upper half plane.
        close(
            C::new(1.0, 2.0).sqrt() * C::new(1.0, 2.0).sqrt(),
            C::new(1.0, 2.0),
            1e-14,
        );
        close(C::real(2.0).exp(), C::real(2.0_f64.exp()), 1e-13);
        let quarter = (C::i() * std::f64::consts::FRAC_PI_2).exp();
        assert!(quarter.re.abs() < 1e-15 && (quarter.im - 1.0).abs() < 1e-15);
    }
}
