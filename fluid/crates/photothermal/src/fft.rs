//! Minimal in-place radix-2 FFT, used to launch a photoacoustic pulse into a
//! layered stack in the frequency domain.
//!
//! The transform convention is the standard signal-processing one,
//! `X_k = sum_n x_n exp(-2 pi i k n / N)` and its inverse with `+i` and a
//! `1/N`, so a real, even input has a real, even spectrum and a pure time
//! shift is a phase ramp. The implementation is iterative Cooley-Tukey; it is
//! deliberately tiny and dependency-free, since the photoacoustic benchmark
//! only needs a few thousand points.

use std::f64::consts::PI;

/// A small complex number for the transform.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Z {
    pub re: f64,
    pub im: f64,
}

impl Z {
    pub const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    pub const fn real(re: f64) -> Self {
        Self { re, im: 0.0 }
    }
    pub fn abs2(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    pub fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }
    pub fn add(self, b: Self) -> Self {
        Self::new(self.re + b.re, self.im + b.im)
    }
    pub fn sub(self, b: Self) -> Self {
        Self::new(self.re - b.re, self.im - b.im)
    }
    pub fn mul(self, b: Self) -> Self {
        Self::new(
            self.re * b.re - self.im * b.im,
            self.re * b.im + self.im * b.re,
        )
    }
    pub fn scale(self, s: f64) -> Self {
        Self::new(self.re * s, self.im * s)
    }
}

/// In-place FFT. `a.len()` must be a power of two. `inverse` selects the
/// opposite sign and applies the `1/N`.
pub fn fft(a: &mut [Z], inverse: bool) {
    let n = a.len();
    assert!(n.is_power_of_two(), "FFT length must be a power of two");
    if n <= 1 {
        return;
    }
    // Bit-reversal permutation.
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            a.swap(i, j);
        }
    }
    // Butterflies.
    let sign = if inverse { 1.0 } else { -1.0 };
    let mut len = 2usize;
    while len <= n {
        let ang = sign * 2.0 * PI / len as f64;
        let wlen = Z::new(ang.cos(), ang.sin());
        let mut i = 0;
        while i < n {
            let mut w = Z::real(1.0);
            for k in 0..len / 2 {
                let u = a[i + k];
                let v = a[i + k + len / 2].mul(w);
                a[i + k] = u.add(v);
                a[i + k + len / 2] = u.sub(v);
                w = w.mul(wlen);
            }
            i += len;
        }
        len <<= 1;
    }
    if inverse {
        let inv = 1.0 / n as f64;
        for x in a.iter_mut() {
            *x = x.scale(inv);
        }
    }
}

/// Forward transform of a real signal.
pub fn fft_real(x: &[f64]) -> Vec<Z> {
    let mut a: Vec<Z> = x.iter().map(|&v| Z::real(v)).collect();
    fft(&mut a, false);
    a
}

/// Inverse transform taking the real part.
pub fn ifft_real(x: &[Z]) -> Vec<f64> {
    let mut a = x.to_vec();
    fft(&mut a, true);
    a.iter().map(|z| z.re).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_and_known_spectrum() {
        let n = 16;
        let x: Vec<f64> = (0..n).map(|i| (i as f64 * 0.4).sin() + 0.5).collect();
        let back = ifft_real(&fft_real(&x));
        for (a, b) in x.iter().zip(&back) {
            assert!((a - b).abs() < 1e-12);
        }
        // A constant has all its energy at k = 0.
        let flat = vec![1.0f64; n];
        let spec = fft_real(&flat);
        assert!((spec[0].re - n as f64).abs() < 1e-12);
        for z in &spec[1..] {
            assert!(z.abs() < 1e-12);
        }
    }
}
