//! Fresnel integrals `F(x) = C(x) + i S(x) = ∫₀ˣ exp(iπt²/2) dt`.
//!
//! `F(x) = ((1+i)/2) erf(√π (1−i) x/2)`. For `|x| ≤ 2.5` the power series
//! `Σ (iπ/2)ⁿ x^{2n+1}/(n!(2n+1))` is summed directly (largest term
//! `~e^{πx²/2} ≈ 2·10³`, so cancellation costs < 1e-12); beyond, the
//! complementary error function is evaluated by its continued fraction
//! (modified Lentz), which converges for `Re z > 0` and here has
//! `arg z = −π/4`. `C` and `S` are odd.

use num_complex::Complex64 as Cx;
use std::f64::consts::PI;

/// `C(x) + i S(x)`.
pub fn fresnel(x: f64) -> Cx {
    if x < 0.0 {
        return -fresnel(-x);
    }
    if x <= 2.5 {
        let mut term = Cx::new(x, 0.0);
        let mut sum = term;
        let z = Cx::new(0.0, PI / 2.0 * x * x);
        for n in 1..200 {
            // term_n = (iπx²/2)ⁿ x / n!, contribution term_n/(2n+1)
            term = term * z / n as f64;
            let c = term / (2 * n + 1) as f64;
            sum += c;
            if c.norm() < 1e-17 * sum.norm() {
                break;
            }
        }
        return sum;
    }
    let z = Cx::new(1.0, -1.0) * (PI.sqrt() * x / 2.0);
    let one_plus_i_half = Cx::new(0.5, 0.5);
    one_plus_i_half * (Cx::new(1.0, 0.0) - erfc_cf(z))
}

/// `erfc(z)` for `Re z > 0` by the continued fraction
/// `erfc z = e^{−z²}/√π · 1/(z + (1/2)/(z + 1/(z + (3/2)/(z + …))))`.
fn erfc_cf(z: Cx) -> Cx {
    let tiny = 1e-300;
    // f = z + a1/(z + a2/(z + …)), a_n = n/2 (modified Lentz).
    let mut f = z;
    if f.norm() < tiny {
        f = Cx::new(tiny, 0.0);
    }
    let mut c = f;
    let mut d = Cx::new(0.0, 0.0);
    for n in 1..5000 {
        let a = n as f64 / 2.0;
        d = z + a * d;
        if d.norm() < tiny {
            d = Cx::new(tiny, 0.0);
        }
        c = z + a / c;
        if c.norm() < tiny {
            c = Cx::new(tiny, 0.0);
        }
        d = d.inv();
        let delta = c * d;
        f *= delta;
        if (delta - 1.0).norm() < 1e-16 {
            break;
        }
    }
    (-z * z).exp() / (PI.sqrt() * f)
}

/// `∫_{x1}^{x2} exp(iπt²/2) dt = F(x2) − F(x1)`, evaluated without
/// cancellation for long intervals far from the origin.
pub fn fresnel_span(x1: f64, x2: f64) -> Cx {
    fresnel(x2) - fresnel(x1)
}
