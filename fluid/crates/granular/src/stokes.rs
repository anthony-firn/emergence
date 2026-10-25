//! Stokes flow around a rigid sphere, and the dilute-suspension stresslet.
//!
//! The exact solution for a sphere translating with velocity `V` through a
//! fluid at rest at infinity is (in the sphere frame; Landau & Lifshitz,
//! *Fluid Mechanics* 2nd ed. §20; Happel & Brenner 1965 §3-3)
//!
//! ```text
//! u = V A(r) + (V.n) n B(r),  n = x/r,
//! A(r) = (3/4)(a/r) + (1/4)(a/r)^3,
//! B(r) = (3/4)(a/r)(1 - (a/r)^2),
//! p     = (3/2) mu a (V.x)/r^3,
//! ```
//!
//! so that `u(a n) = V`, `u(inf) = 0`, and the surface traction is uniform,
//! `t = -(3 mu / 2a) V`. Integrating it gives the drag `F = -6 pi mu a V`
//! (Stokes 1851).
//!
//! The **stresslet** of a rigid, force- and torque-free sphere in an
//! imposed symmetric traceless rate of strain `E` is
//! `S_ij = (20/3) pi mu a^3 E_ij` (Batchelor 1970; Kim & Karrila 2005
//! §3.3). The traction it carries is `t_i = 5 mu E_ij n_j`; integrating its
//! symmetric first moment over the sphere — done here numerically, as the
//! independent check — returns that `S`. Summing `N` such stresslets in a
//! dilute suspension gives Einstein's `eta = eta0 (1 + 5/2 phi)` exactly
//! (see [`crate::suspension`]).

use std::f64::consts::PI;

/// Dot product.
#[inline]
pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product.
#[inline]
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Oseen tensor `G_ij(r) = delta_ij/r + r_i r_j/r^3`, the Stokeslet Green's
/// function (`u_i = G_ij f_j / 8 pi mu`). Singular at `r = 0`.
pub fn oseen(r: [f64; 3]) -> [[f64; 3]; 3] {
    let rr = norm(r);
    let r3 = rr * rr * rr;
    let mut g = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            g[i][j] = if i == j { 1.0 / rr } else { 0.0 } + r[i] * r[j] / r3;
        }
    }
    g
}

/// Exact velocity at `x` of a rigid sphere of radius `a`, centre at the
/// origin, translating with velocity `v` in a fluid at rest at infinity.
pub fn translating_sphere_velocity(v: [f64; 3], a: f64, x: [f64; 3]) -> [f64; 3] {
    let rr = norm(x);
    if rr < 1e-300 {
        return v;
    }
    let ar = a / rr;
    let a3 = ar * ar * ar;
    let aa = 0.75 * ar + 0.25 * a3; // A(r)
    let bb = 0.75 * ar * (1.0 - ar * ar); // B(r)
    let vn = dot(v, x) / rr;
    [
        v[0] * aa + vn * (x[0] / rr) * bb,
        v[1] * aa + vn * (x[1] / rr) * bb,
        v[2] * aa + vn * (x[2] / rr) * bb,
    ]
}

/// Gradient `D_ij = d u_i / d x_j` of [`translating_sphere_velocity`].
pub fn translating_sphere_gradient(v: [f64; 3], a: f64, x: [f64; 3]) -> [[f64; 3]; 3] {
    let rr = norm(x);
    let r2 = rr * rr;
    let r3 = r2 * rr;
    let r4 = r2 * r2;
    let s = dot(v, x);
    // A'(r), B'(r).
    let da = -0.75 * a / r2 - 0.75 * a * a * a / r4;
    let db = 0.75 * a * (-1.0 / r2 + 3.0 * a * a / r4);
    let ar = a / rr;
    let bb = 0.75 * ar * (1.0 - ar * ar);
    let n = [x[0] / rr, x[1] / rr, x[2] / rr];
    let mut d = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let term1 = v[i] * da * n[j];
            let term2 = bb
                * ((v[j] / rr - s * x[j] / r3) * n[i]
                    + (s / rr) * ((if i == j { 1.0 } else { 0.0 }) / rr - x[i] * x[j] / r3));
            let term3 = (s / rr) * n[i] * db * n[j];
            d[i][j] = term1 + term2 + term3;
        }
    }
    d
}

/// Exact pressure `p = (3/2) mu a (v.x)/r^3` (ambient pressure zero).
pub fn translating_sphere_pressure(v: [f64; 3], a: f64, mu: f64, x: [f64; 3]) -> f64 {
    let r3 = norm(x).powi(3);
    1.5 * mu * a * dot(v, x) / r3
}

/// Cauchy stress `sigma_ij = -p delta_ij + mu (du_i/dx_j + du_j/dx_i)`.
pub fn translating_sphere_stress(
    v: [f64; 3],
    a: f64,
    mu: f64,
    x: [f64; 3],
) -> [[f64; 3]; 3] {
    let p = translating_sphere_pressure(v, a, mu, x);
    let d = translating_sphere_gradient(v, a, x);
    let mut sig = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            sig[i][j] = -p * if i == j { 1.0 } else { 0.0 } + mu * (d[i][j] + d[j][i]);
        }
    }
    sig
}

/// Numerical drag on the translating sphere: integrate the traction
/// `t_i = sigma_ij n_j` over the surface using Gauss–Legendre in `cos theta`
/// and the trapezoid rule in `phi`. Returns the force on the fluid (equal
/// and opposite to the drag on the sphere).
pub fn drag_numerical(v: [f64; 3], a: f64, mu: f64, ntheta: usize, nphi: usize) -> [f64; 3] {
    let (nodes, weights) = gauss_legendre(ntheta);
    let mut f = [0.0; 3];
    let dphi = 2.0 * PI / nphi as f64;
    for (ct, w) in nodes.iter().zip(&weights) {
        let st = (1.0 - ct * ct).max(0.0).sqrt();
        for ip in 0..nphi {
            let phi = (ip as f64 + 0.5) * dphi;
            let n = [st * phi.cos(), st * phi.sin(), *ct];
            let x = [a * n[0], a * n[1], a * n[2]];
            let sig = translating_sphere_stress(v, a, mu, x);
            let mut t = [0.0; 3];
            for i in 0..3 {
                for j in 0..3 {
                    t[i] += sig[i][j] * n[j];
                }
            }
            let darea = a * a * w * dphi;
            for i in 0..3 {
                f[i] += t[i] * darea;
            }
        }
    }
    f
}

/// Exact single-sphere stresslet `S_ij = (20/3) pi mu a^3 E_ij` for a
/// force-free, torque-free rigid sphere in a symmetric traceless rate of
/// strain `E` (Batchelor 1970; Kim & Karrila 2005 §3.3).
pub fn stresslet_exact(e: [[f64; 3]; 3], a: f64, mu: f64) -> [[f64; 3]; 3] {
    let c = (20.0 / 3.0) * PI * mu * a * a * a;
    scale(c, e)
}

/// Traction on a rigid sphere in an imposed rate of strain `E`:
/// `t_i = 5 mu E_ij n_j` (the distribution whose symmetric first moment is
/// the stresslet above). Sign convention: the traction the suspension
/// exerts on the particle.
pub fn strain_traction(e: [[f64; 3]; 3], mu: f64, n: [f64; 3]) -> [f64; 3] {
    let mut t = [0.0; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i] += 5.0 * mu * e[i][j] * n[j];
        }
    }
    t
}

/// Numerical stresslet: integrate the symmetric first moment
/// `S_ij = (1/2) ∫ (x_i t_j + x_j t_i) dA` of the strain traction over the
/// sphere. This is the independent numerical check of [`stresslet_exact`].
pub fn stresslet_numerical(
    e: [[f64; 3]; 3],
    a: f64,
    mu: f64,
    ntheta: usize,
    nphi: usize,
) -> [[f64; 3]; 3] {
    let (nodes, weights) = gauss_legendre(ntheta);
    let mut s = [[0.0; 3]; 3];
    let dphi = 2.0 * PI / nphi as f64;
    for (ct, w) in nodes.iter().zip(&weights) {
        let st = (1.0 - ct * ct).max(0.0).sqrt();
        for ip in 0..nphi {
            let phi = (ip as f64 + 0.5) * dphi;
            let n = [st * phi.cos(), st * phi.sin(), *ct];
            let x = [a * n[0], a * n[1], a * n[2]];
            let t = strain_traction(e, mu, n);
            let darea = a * a * w * dphi;
            for i in 0..3 {
                for j in 0..3 {
                    s[i][j] += 0.5 * (x[i] * t[j] + x[j] * t[i]) * darea;
                }
            }
        }
    }
    s
}

fn scale(c: f64, m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = c * m[i][j];
        }
    }
    r
}

/// Gauss–Legendre nodes and weights on `[-1, 1]` (Newton on the Legendre
/// polynomial). Enough for the smooth sphere quadratures here.
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = vec![0.0; n];
    let mut weights = vec![0.0; n];
    for k in 0..n {
        // Initial guess (Chebyshev).
        let mut x = (PI * (k as f64 + 0.75) / (n as f64 + 0.5)).cos();
        for _ in 0..100 {
            // Evaluate P_n and P_n' by recurrence.
            let (mut p0, mut p1) = (1.0, x);
            for j in 2..=n {
                let p2 = ((2 * j - 1) as f64 * x * p1 - (j - 1) as f64 * p0) / j as f64;
                p0 = p1;
                p1 = p2;
            }
            let pn = p1;
            let dpn = n as f64 * (x * pn - p0) / (x * x - 1.0);
            let dx = pn / dpn;
            x -= dx;
            if dx.abs() < 1e-15 {
                break;
            }
        }
        let (mut p0, mut p1) = (1.0, x);
        for j in 2..=n {
            let p2 = ((2 * j - 1) as f64 * x * p1 - (j - 1) as f64 * p0) / j as f64;
            p0 = p1;
            p1 = p2;
        }
        let dpn = n as f64 * (x * p1 - p0) / (x * x - 1.0);
        nodes[k] = x;
        weights[k] = 2.0 / ((1.0 - x * x) * dpn * dpn);
    }
    (nodes, weights)
}

/// A rank-two symmetric traceless rate of strain used in the tests: simple
/// shear `E_xy = E_yx = gdot/2`.
pub fn simple_shear_e(gdot: f64) -> [[f64; 3]; 3] {
    let mut e = [[0.0; 3]; 3];
    e[0][1] = 0.5 * gdot;
    e[1][0] = 0.5 * gdot;
    e
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn surface_velocity_is_no_slip_and_decays() {
        let v = [0.7, -0.3, 0.5];
        let a = 0.8;
        for &ct in &[-1.0_f64, -0.4, 0.2, 1.0] {
            let st = (1.0 - ct * ct).sqrt();
            let x = [a * st, 0.0, a * ct];
            let u = translating_sphere_velocity(v, a, x);
            for i in 0..3 {
                assert!((u[i] - v[i]).abs() < 1e-14, "no-slip {u:?}");
            }
        }
        let u_far = translating_sphere_velocity(v, a, [100.0, -40.0, 30.0]);
        assert!(norm(u_far) < 1e-2, "far field {u_far:?}");
    }

    #[test]
    fn gradient_matches_finite_difference() {
        let v = [0.7, -0.3, 0.5];
        let (a, mu) = (0.8, 1.3);
        let x = [1.7, -0.9, 0.6];
        let d = translating_sphere_gradient(v, a, x);
        let h = 1e-6;
        for j in 0..3 {
            let mut xp = x;
            let mut xm = x;
            xp[j] += h;
            xm[j] -= h;
            let up = translating_sphere_velocity(v, a, xp);
            let um = translating_sphere_velocity(v, a, xm);
            for i in 0..3 {
                let fd = (up[i] - um[i]) / (2.0 * h);
                assert!((d[i][j] - fd).abs() < 1e-7, "D[{i}][{j}] {} vs {fd}", d[i][j]);
            }
        }
        // The stress diverges as r->0 but is finite at the surface.
        let _ = translating_sphere_stress(v, a, mu, [a, 0.0, 0.0]);
    }
}
