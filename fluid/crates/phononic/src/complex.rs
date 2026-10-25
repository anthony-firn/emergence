//! Minimal complex arithmetic for wave mechanics.
//!
//! Two things make this crate complex-valued. A lossy layer has a complex
//! modulus `M = M0 (1 + i eta)`; the imaginary part is mechanical absorption.
//! And even in a lossless stack a wave below its critical angle, or inside a
//! band gap, has an imaginary normal wavenumber: the field is evanescent, not
//! zero.
//!
//! The type is deliberately small (no dependency) and matches the convention of
//! the wave-optics crate: principal `sqrt` with a branch cut on the negative
//! real axis, upper half plane. For a passive modulus `M = M0(1 + i eta)` and
//! time dependence `exp(+i omega t)`, the principal root gives
//! `Im(k) >= 0`, so the forward wave `exp(-i k z)` decays along `+z`.
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
    /// Principal square root with a non-negative real part (propagating-wave
    /// branch). Opposite half-plane choice used for a wave vector that must
    /// point in a chosen direction.
    pub fn sqrt_nonneg_re(self) -> Self {
        let s = self.sqrt();
        if s.re >= 0.0 {
            s
        } else {
            -s
        }
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

impl From<f64> for C {
    fn from(v: f64) -> Self {
        Self::real(v)
    }
}

macro_rules! binop {
    ($trait:ident, $method:ident, $op:tt) => {
        impl $trait for C {
            type Output = C;
            fn $method(self, b: C) -> C {
                C::new(self.re $op b.re, self.im $op b.im)
            }
        }
    };
}
binop!(Add, add, +);
binop!(Sub, sub, -);

impl Mul for C {
    type Output = C;
    fn mul(self, b: C) -> C {
        C::new(
            self.re * b.re - self.im * b.im,
            self.re * b.im + self.im * b.re,
        )
    }
}
impl Div for C {
    type Output = C;
    fn div(self, b: C) -> C {
        let d = b.abs2();
        C::new(
            (self.re * b.re + self.im * b.im) / d,
            (self.im * b.re - self.re * b.im) / d,
        )
    }
}
impl Neg for C {
    type Output = C;
    fn neg(self) -> C {
        C::new(-self.re, -self.im)
    }
}
impl Mul<f64> for C {
    type Output = C;
    fn mul(self, b: f64) -> C {
        self.scale(b)
    }
}
impl Div<f64> for C {
    type Output = C;
    fn div(self, b: f64) -> C {
        self.scale(1.0 / b)
    }
}
impl Add<f64> for C {
    type Output = C;
    fn add(self, b: f64) -> C {
        C::new(self.re + b, self.im)
    }
}
impl Sub<f64> for C {
    type Output = C;
    fn sub(self, b: f64) -> C {
        C::new(self.re - b, self.im)
    }
}

/// Solve a small dense complex linear system `A x = b` by Gaussian elimination
/// with partial pivoting. `a` is row-major `n x n`. Returns `None` if singular.
pub fn solve(a: &[C], b: &[C], n: usize) -> Option<Vec<C>> {
    let mut m = a.to_vec();
    let mut x = b.to_vec();
    for col in 0..n {
        // Pivot on the largest magnitude in this column.
        let mut piv = col;
        for row in col + 1..n {
            if m[row * n + col].abs() > m[piv * n + col].abs() {
                piv = row;
            }
        }
        if m[piv * n + col].abs() < 1e-300 {
            return None;
        }
        if piv != col {
            for k in 0..n {
                m.swap(col * n + k, piv * n + k);
            }
            x.swap(col, piv);
        }
        let d = m[col * n + col];
        for k in col..n {
            m[col * n + k] = m[col * n + k] / d;
        }
        x[col] = x[col] / d;
        for row in 0..n {
            if row == col {
                continue;
            }
            let f = m[row * n + col];
            if f.abs() == 0.0 {
                continue;
            }
            for k in col..n {
                m[row * n + k] = m[row * n + k] - f * m[col * n + k];
            }
            x[row] = x[row] - f * x[col];
        }
    }
    Some(x)
}

/// Particular solution of a possibly rank-deficient complex system, with free
/// variables set to zero. Used for boundary systems that contain an unforced
/// free vibration (for example the shear mode at exactly normal incidence),
/// where the physical reflection/transmission is independent of that free
/// amplitude. The system is assumed consistent.
pub fn solve_particular(a: &[C], b: &[C], n: usize) -> Vec<C> {
    let mut m = a.to_vec();
    let mut rhs = b.to_vec();
    let mut pivot_col = vec![usize::MAX; n];
    let mut row = 0usize;
    for col in 0..n {
        if row >= n {
            break;
        }
        let mut best = row;
        for r in row + 1..n {
            if m[r * n + col].abs() > m[best * n + col].abs() {
                best = r;
            }
        }
        // Skip the column only if it is exactly zero below the current row.
        // A magnitude-relative threshold is wrong here because the boundary
        // matrix mixes stresses (~1e13) with displacements (~1), so a small
        // but valid pivot can be many orders below the global scale.
        let mut colmax = 0.0f64;
        for r in row..n {
            colmax = colmax.max(m[r * n + col].abs());
        }
        if colmax <= 0.0 {
            continue; // free variable
        }
        if best != row {
            for k in 0..n {
                m.swap(row * n + k, best * n + k);
            }
            rhs.swap(row, best);
        }
        let d = m[row * n + col];
        for k in col..n {
            m[row * n + k] = m[row * n + k] / d;
        }
        rhs[row] = rhs[row] / d;
        for r in 0..n {
            if r == row {
                continue;
            }
            let f = m[r * n + col];
            if f.abs() == 0.0 {
                continue;
            }
            for k in col..n {
                m[r * n + k] = m[r * n + k] - f * m[row * n + k];
            }
            rhs[r] = rhs[r] - f * rhs[row];
        }
        pivot_col[row] = col;
        row += 1;
    }
    let mut x = vec![C::real(0.0); n];
    for r in 0..row {
        x[pivot_col[r]] = rhs[r];
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: C, b: C, tol: f64) {
        assert!((a - b).abs() <= tol, "{a:?} vs {b:?}");
    }

    #[test]
    fn particular_solution_residual() {
        let a = [
            C::real(0.0),
            C::new(0.0, -2.0),
            C::real(1.0),
            C::real(5.0),
            C::real(0.0),
            C::real(0.0),
            C::real(0.0),
            C::new(-3.0, 1.0),
            C::real(0.8),
        ];
        let b = [C::real(-2.0), C::real(0.0), C::real(0.0)];
        let x = solve_particular(&a, &b, 3);
        for r in 0..3 {
            let mut acc = C::real(0.0);
            for c in 0..3 {
                acc = acc + a[r * 3 + c] * x[c];
            }
            close(acc, b[r], 1e-12);
        }
    }

    #[test]
    fn arithmetic_and_principal_sqrt() {
        let z = C::new(3.0, 4.0);
        close(z * z.recip(), C::real(1.0), 1e-15);
        close(z.sqrt() * z.sqrt(), z, 1e-14);
        // Passive loss: M = M0(1 + i eta), forward exp(-i k z) must decay.
        let m = C::new(1.0, 0.1);
        let k = m.sqrt(); // proportional to wavenumber
        assert!(k.im > 0.0, "passive loss must give Im k > 0, got {k:?}");
    }

    #[test]
    fn solves_a_small_system() {
        // [[1, i],[2, 1]] x = [1+2i, 3].
        let a = [C::real(1.0), C::i(), C::real(2.0), C::real(1.0)];
        let b = [C::new(1.0, 2.0), C::real(3.0)];
        let x = solve(&a, &b, 2).unwrap();
        let r0 = a[0] * x[0] + a[1] * x[1];
        let r1 = a[2] * x[0] + a[3] * x[1];
        close(r0, b[0], 1e-14);
        close(r1, b[1], 1e-14);
    }
}
