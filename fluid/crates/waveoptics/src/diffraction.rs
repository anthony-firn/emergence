//! Scalar diffraction by rectangular apertures: a ladder of models with
//! rigorous error claims.
//!
//! The aperture lies in `z = 0`; the observation point `P = (x, y, L)` has
//! `L > 0`. With Kirchhoff boundary values (the incident field in the
//! opening, zero on the screen), the exact scalar field is the first
//! Rayleigh–Sommerfeld integral
//!
//! ```text
//! U(P) = 1/(iλ) ∫∫ U_inc(Q) (L/r) (1 − 1/(ikr)) e^{ikr}/r dA,   r = |P − Q|.
//! ```
//!
//! | rung | model | claim |
//! |---|---|---|
//! | Fraunhofer | linear phase across the aperture (sinc × sinc) | rigorous: Fresnel bound + `k q (a_x² + a_y²)/8` |
//! | Fresnel | quadratic phase, closed form in Fresnel integrals, exact phase and amplitude at the aperture centre | rigorous: variation of the neglected terms across the aperture |
//! | Rayleigh–Sommerfeld | the integral above by Gauss–Legendre panels resolving the phase | a posteriori (refined rule) |
//!
//! **Why centre corrections.** The textbook Fresnel condition bounds the
//! whole neglected quartic phase `k ρ⁴/(8L³)`. For an off-axis point most of
//! it is a constant across a small aperture: it moves the phase, not the
//! pattern. It is computed exactly at the aperture centre; only its
//! variation across the aperture is an error. Without it the two slits of a
//! double slit get different, uncorrected phase offsets off axis and the
//! fringes move.
//!
//! **Bounds used.** `g(ρ) = √(z²+ρ²) − z − ρ²/(2z)` satisfies
//! `−ρ⁴/(8z³) ≤ g ≤ −ρ⁴/(8z³) + ρ⁶/(16z⁵)`, so over the aperture
//! `|g − g(centre)| ≤ (ρ⁴_max − ρ⁴_min)/(8z³) + ρ⁶_max/(16z⁵)`, with
//! `ρ²` ranges exact for a rectangle. The amplitude factor's variation is
//! bounded from the same ranges, and `|e^{iδ} − 1| ≤ |δ|`.
//!
//! **Scalar validity.** Kirchhoff boundary values are wrong within a few
//! wavelengths of the edges; for an opening of width `a` the relative error
//! of the far field is `O(λ/a)`. [`Rect::scalar_validity`] reports `a/λ`.
//! A slit only a few wavelengths wide needs a Maxwell solver, which is not
//! built; such apertures are flagged, not silently computed.

use std::f64::consts::PI;

use emergence::{Assessment, Audit, Claim, Ladder, Rung};
use num_complex::Complex64 as Cx;

use crate::fresnel::fresnel_span;

/// A rectangular opening in the plane `z = 0`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rect {
    pub center: [f64; 2],
    pub size: [f64; 2],
}

impl Rect {
    pub fn area(&self) -> f64 {
        self.size[0] * self.size[1]
    }
    fn lo(&self, a: usize) -> f64 {
        self.center[a] - 0.5 * self.size[a]
    }
    fn hi(&self, a: usize) -> f64 {
        self.center[a] + 0.5 * self.size[a]
    }
    /// Smallest opening width in wavelengths (`a/λ`).
    pub fn scalar_validity(&self, wavelength: f64) -> f64 {
        self.size[0].min(self.size[1]) / wavelength
    }
}

/// The wave arriving at the aperture plane.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Incident {
    /// `A e^{ik(s_x ξ + s_y η)}`: a plane wave with direction sines `s`.
    Plane { amplitude: Cx, sines: [f64; 2] },
    /// `A e^{ikR}/R` from a point `(x_s, y_s)` at distance `z_s` before the
    /// aperture plane.
    Point { amplitude: Cx, position: [f64; 2], distance: f64 },
}

/// One evaluation: the field of one aperture at one point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Query {
    pub aperture: Rect,
    pub incident: Incident,
    pub point: [f64; 3],
    pub wavelength: f64,
}

impl Query {
    fn k(&self) -> f64 {
        2.0 * PI / self.wavelength
    }
}

/// `[min, max]` of `(p − u)²` for `u ∈ [lo, hi]`.
fn sq_range(p: f64, lo: f64, hi: f64) -> (f64, f64) {
    let (a, b) = ((p - lo).powi(2), (p - hi).powi(2));
    let max = a.max(b);
    let min = if p >= lo && p <= hi { 0.0 } else { a.min(b) };
    (min, max)
}

/// Range of `ρ²` between a point and the rectangle (transverse).
fn rho2_range(p: [f64; 2], ap: &Rect) -> (f64, f64) {
    let (x0, x1) = sq_range(p[0], ap.lo(0), ap.hi(0));
    let (y0, y1) = sq_range(p[1], ap.lo(1), ap.hi(1));
    (x0 + y0, x1 + y1)
}

/// Bound on the variation of `g(ρ) = √(z²+ρ²) − z − ρ²/(2z)` over the
/// aperture.
fn g_variation(z: f64, (r2min, r2max): (f64, f64)) -> f64 {
    (r2max * r2max - r2min * r2min) / (8.0 * z.powi(3)) + r2max.powi(3) / (16.0 * z.powi(5))
}

/// Exact scalar kernel factors `(1/R)·(L/r)(1 − 1/(ikr))/r` magnitude
/// pieces, used for centre corrections and bounds.
struct Geometry {
    k: f64,
    l: f64,
    /// Quadratic coefficient of the Fresnel phase: `k q/2`, `q = 1/z_s + 1/L`.
    q: f64,
    /// Fresnel phase `Φ_F(ξ, η)`; also gives the linear coefficients.
    xs: [f64; 2],
    zs: Option<f64>,
    sines: [f64; 2],
    amplitude: Cx,
}

impl Geometry {
    fn new(qr: &Query) -> Self {
        let k = qr.k();
        let l = qr.point[2];
        let (amplitude, xs, zs, sines) = match qr.incident {
            Incident::Plane { amplitude, sines } => (amplitude, [0.0; 2], None, sines),
            Incident::Point { amplitude, position, distance } => (amplitude, position, Some(distance), [0.0; 2]),
        };
        let q = 1.0 / l + zs.map_or(0.0, |z| 1.0 / z);
        Self {
            k,
            l,
            q,
            xs,
            zs,
            sines,
            amplitude,
        }
    }

    /// Fresnel phase at aperture point `(ξ, η)` (without the constant
    /// `k(z_s + L)`).
    fn phi_f(&self, p: [f64; 3], u: [f64; 2]) -> f64 {
        let mut s = 0.0;
        for a in 0..2 {
            s += (p[a] - u[a]).powi(2) / (2.0 * self.l);
            s += match self.zs {
                Some(z) => (u[a] - self.xs[a]).powi(2) / (2.0 * z),
                None => self.sines[a] * u[a],
            };
        }
        self.k * s
    }

    /// Exact phase `k(R + r)` (minus the same constant).
    fn phi_exact(&self, p: [f64; 3], u: [f64; 2]) -> f64 {
        let r = ((p[0] - u[0]).powi(2) + (p[1] - u[1]).powi(2) + self.l * self.l).sqrt() - self.l;
        let inc = match self.zs {
            Some(z) => ((u[0] - self.xs[0]).powi(2) + (u[1] - self.xs[1]).powi(2) + z * z).sqrt() - z,
            None => self.sines[0] * u[0] + self.sines[1] * u[1],
        };
        self.k * (r + inc)
    }

    /// Exact amplitude factor `A/R · (L/r)(1 − 1/(ikr))/r` (complex).
    fn amp_exact(&self, p: [f64; 3], u: [f64; 2]) -> Cx {
        let r = ((p[0] - u[0]).powi(2) + (p[1] - u[1]).powi(2) + self.l * self.l).sqrt();
        let inc = match self.zs {
            Some(z) => 1.0 / ((u[0] - self.xs[0]).powi(2) + (u[1] - self.xs[1]).powi(2) + z * z).sqrt(),
            None => 1.0,
        };
        self.amplitude * inc * (self.l / r) / r * (Cx::new(1.0, 0.0) - Cx::new(0.0, -1.0 / (self.k * r)))
    }

    /// Constant phase `e^{ik(z_s + L)}` common to every rung.
    fn carrier(&self) -> Cx {
        Cx::from_polar(1.0, self.k * (self.l + self.zs.unwrap_or(0.0)))
    }

    fn prefactor(&self, wavelength: f64) -> Cx {
        Cx::new(0.0, -1.0 / wavelength) // 1/(iλ)
    }
}

/// Exact phase gradient `∇_ξ k(r + R)` at aperture point `u`.
fn exact_gradient(g: &Geometry, p: [f64; 3], u: [f64; 2]) -> [f64; 2] {
    let r = ((p[0] - u[0]).powi(2) + (p[1] - u[1]).powi(2) + g.l * g.l).sqrt();
    let rs = g.zs.map(|z| ((u[0] - g.xs[0]).powi(2) + (u[1] - g.xs[1]).powi(2) + z * z).sqrt());
    std::array::from_fn(|a| {
        g.k * (-(p[a] - u[a]) / r
            + match rs {
                Some(rs) => (u[a] - g.xs[a]) / rs,
                None => g.sines[a],
            })
    })
}

/// Fraunhofer rung (tilted): exact phase and phase gradient at the aperture
/// centre, linear phase across the aperture: `A_c e^{iΦ(c)} a_x a_y
/// sinc(α_x a_x/2) sinc(α_y a_y/2)`.
pub fn fraunhofer(qr: &Query) -> Cx {
    let g = Geometry::new(qr);
    let c = qr.aperture.center;
    let p = qr.point;
    let alpha = exact_gradient(&g, p, c);
    let sinc = |x: f64| if x.abs() < 1e-8 { 1.0 - x * x / 6.0 } else { x.sin() / x };
    let integral = qr.aperture.area() * sinc(alpha[0] * qr.aperture.size[0] / 2.0) * sinc(alpha[1] * qr.aperture.size[1] / 2.0);
    g.prefactor(qr.wavelength) * g.carrier() * g.amp_exact(p, c) * Cx::from_polar(integral, g.phi_exact(p, c))
}

/// Fresnel rung: exact integral of the quadratic (Fresnel) phase over the
/// rectangle, corrected to the exact phase and amplitude at the centre.
pub fn fresnel(qr: &Query) -> Cx {
    let g = Geometry::new(qr);
    let c = qr.aperture.center;
    let p = qr.point;
    // Per axis Φ_F = (kq/2)(ξ − ξ₀)² + const, ξ₀ = (x/L + x_s/z_s)/q
    // (plane wave: x/L − s_x).
    let kq2 = g.k * g.q / 2.0;
    let scale = (kq2 / (PI / 2.0)).sqrt();
    let mut integral = Cx::new(1.0, 0.0);
    let mut phi_vertex = g.phi_f(p, c);
    for a in 0..2 {
        let lin = p[a] / g.l + g.zs.map_or(-g.sines[a], |z| g.xs[a] / z);
        let xi0 = lin / g.q;
        let (u1, u2) = ((qr.aperture.lo(a) - xi0) * scale, (qr.aperture.hi(a) - xi0) * scale);
        integral *= fresnel_span(u1, u2) / scale;
        phi_vertex -= kq2 * (c[a] - xi0).powi(2);
    }
    let delta = g.phi_exact(p, c) - g.phi_f(p, c);
    g.prefactor(qr.wavelength) * g.carrier() * g.amp_exact(p, c) * Cx::from_polar(1.0, delta + phi_vertex) * integral
}

/// Rigorous error bounds `(fresnel, fraunhofer)` on the field.
pub fn bounds(qr: &Query) -> (f64, f64) {
    let g = Geometry::new(qr);
    let p = qr.point;
    let ap = &qr.aperture;
    let ro = rho2_range([p[0], p[1]], ap);
    let mut dphi = g.k * g_variation(g.l, ro);
    let (rmin, rmax) = ((g.l * g.l + ro.0).sqrt(), (g.l * g.l + ro.1).sqrt());
    let (mut smin, mut smax, mut inv_rs_min) = (1.0, 1.0, 0.0);
    if let Some(z) = g.zs {
        let rs = rho2_range(g.xs, ap);
        dphi += g.k * g_variation(z, rs);
        smin = 1.0 / (z * z + rs.1).sqrt();
        smax = 1.0 / (z * z + rs.0).sqrt();
        inv_rs_min = smax;
    }
    // amplitude factor (1/R)(L/r²)(1 − 1/(ikr)): range over the aperture
    let m_max = smax * g.l / (rmin * rmin);
    let m_min = smin * g.l / (rmax * rmax);
    let amp_var = (m_max - m_min) + 2.0 * m_max / (g.k * rmin);
    let scale = g.amplitude.norm() * ap.area() / qr.wavelength;
    let fres = scale * (amp_var + m_max * dphi);
    // Tilted Fraunhofer: Hessian of k(r + R) has eigenvalues in
    // [0, k(1/r + 1/R)], so the neglected quadratic is ≤ k(1/r_min +
    // 1/R_min)(a_x² + a_y²)/8.
    let quad = g.k * (1.0 / rmin + inv_rs_min) * (ap.size[0].powi(2) + ap.size[1].powi(2)) / 8.0;
    (fres, scale * (amp_var + m_max * quad))
}

/// Order of the local expansion of the phase about the aperture centre.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Order {
    /// Linear phase (tilted Fraunhofer).
    Linear,
    /// Linear + diagonal quadratic phase (tilted Fresnel).
    Quadratic,
}

/// A closed-form approximation and its rigorous error bound.
#[derive(Clone, Copy, Debug)]
pub struct Approx {
    pub u: Cx,
    pub bound: f64,
}

/// `(∫ e^{iφ}, ∫ s e^{iφ})` over `s ∈ [−c, c]` with `φ = g s + h s²/2`,
/// and the phase error committed when `h` is dropped (`0` otherwise).
pub(crate) fn moments(g: f64, h: f64, c: f64, keep_quadratic: bool) -> (Cx, Cx, f64) {
    let i = Cx::new(0.0, 1.0);
    if keep_quadratic && h * c * c > 1e-4 {
        // Complete the square: φ = (h/2)(s + g/h)² − g²/(2h).
        let sq = (h / PI).sqrt();
        let v = |s: f64| (s + g / h) * sq;
        let i0 = Cx::from_polar(1.0, -g * g / (2.0 * h)) * fresnel_span(v(-c), v(c)) / sq;
        let phi = |s: f64| g * s + 0.5 * h * s * s;
        let i1 = (Cx::from_polar(1.0, phi(c)) - Cx::from_polar(1.0, phi(-c))) / (i * h) - i0 * (g / h);
        return (i0, i1, 0.0);
    }
    // Linear phase: exact moments; dropping h costs ≤ h c²/2.
    let x = g * c;
    let (i0, i1) = if x.abs() < 1e-3 {
        let x2 = x * x;
        (Cx::new(2.0 * c * (1.0 - x2 / 6.0 + x2 * x2 / 120.0), 0.0), i * (2.0 * g * c * c * c / 3.0) * (1.0 - x2 / 10.0 + x2 * x2 / 280.0))
    } else {
        (Cx::new(2.0 * x.sin() / g, 0.0), i * (2.0 * (x.sin() - x * x.cos()) / (g * g)))
    };
    (i0, i1, 0.5 * h * c * c)
}

/// Local expansion about the aperture centre: exact phase, gradient and
/// diagonal Hessian of `k(r + R)`, exact linear amplitude variation; the
/// integral over the rectangle is then a product of closed forms. The
/// bound covers everything dropped:
///
/// * the Hessian cross term `|H_ξη(c)| (a_x/2)(a_y/2)`;
/// * the third-order phase remainder `(k/6)(D₃r + D₃R)|Δ|³`, using
///   `|d³r/ds³| = 3|c|(1 − c²)/r² ≤ min(3ρ/r³, 1.155/r²)` along any
///   transverse direction (`c` the direction cosine);
/// * the quadratic itself for [`Order::Linear`] (`Σ H_aa a_a²/8`);
/// * the second-order amplitude remainder `max|m''| |Δ|²/2` and the
///   variation of the `1 − 1/(ikr)` factor.
pub fn local_expansion(qr: &Query, order: Order) -> Approx {
    local_impl(qr, order, true)
}

/// The bound of [`local_expansion`] without computing the field (cheap).
pub fn local_bound(qr: &Query, order: Order) -> f64 {
    local_impl(qr, order, false).bound
}

fn local_impl(qr: &Query, order: Order, want_field: bool) -> Approx {
    let g = Geometry::new(qr);
    let c = qr.aperture.center;
    let p = qr.point;
    let ap = &qr.aperture;
    let w = [c[0] - p[0], c[1] - p[1]];
    let r = (w[0] * w[0] + w[1] * w[1] + g.l * g.l).sqrt();
    let (ws, rs) = match g.zs {
        Some(z) => {
            let ws = [c[0] - g.xs[0], c[1] - g.xs[1]];
            (ws, Some((ws[0] * ws[0] + ws[1] * ws[1] + z * z).sqrt()))
        }
        None => ([0.0; 2], None),
    };
    let grad = exact_gradient(&g, p, c);
    // Hessian of k(r + R): k[(δ − ŵŵᵀ)/r] on the transverse block.
    let hess = |a: usize, b: usize| -> f64 {
        let d = if a == b { 1.0 } else { 0.0 };
        let mut h = (d - w[a] * w[b] / (r * r)) / r;
        if let Some(rs) = rs {
            h += (d - ws[a] * ws[b] / (rs * rs)) / rs;
        }
        g.k * h
    };
    // Amplitude m = A (1/R)(L/r²), linear part β = ∇ ln m.
    let beta: [f64; 2] = std::array::from_fn(|a| -2.0 * w[a] / (r * r) - rs.map_or(0.0, |rs| ws[a] / (rs * rs)));
    let keep = order == Order::Quadratic;
    let half = [0.5 * ap.size[0], 0.5 * ap.size[1]];
    let mut drop_q = 0.0;
    let mut u = Cx::new(0.0, 0.0);
    if want_field {
        let (mut i0, mut i1) = ([Cx::new(0.0, 0.0); 2], [Cx::new(0.0, 0.0); 2]);
        for a in 0..2 {
            let (m0, m1, dq) = moments(grad[a], hess(a, a), half[a], keep);
            i0[a] = m0;
            i1[a] = m1;
            drop_q += dq;
        }
        let integral = i0[0] * i0[1] + i1[0] * i0[1] * beta[0] + i0[0] * i1[1] * beta[1];
        u = g.prefactor(qr.wavelength) * g.carrier() * g.amp_exact(p, c) * Cx::from_polar(1.0, g.phi_exact(p, c)) * integral;
    } else if keep {
        for a in 0..2 {
            let h = hess(a, a);
            if h * half[a] * half[a] <= 1e-4 {
                drop_q += 0.5 * h * half[a] * half[a];
            }
        }
    }

    // ---- bound
    let ro = rho2_range([p[0], p[1]], ap);
    let (rmin, rmax) = ((g.l * g.l + ro.0).sqrt(), (g.l * g.l + ro.1).sqrt());
    let rho_max = ro.1.sqrt();
    let d3 = |rho: f64, rmin: f64| (3.0 * rho / rmin.powi(3)).min(1.155 / (rmin * rmin));
    let mut d3sum = d3(rho_max, rmin);
    let (mut rsmin, mut rhos_max) = (f64::INFINITY, 0.0);
    if let Some(z) = g.zs {
        let rr = rho2_range(g.xs, ap);
        rsmin = (z * z + rr.0).sqrt();
        rhos_max = rr.1.sqrt();
        d3sum += d3(rhos_max, rsmin);
    }
    let diag = (half[0] * half[0] + half[1] * half[1]).sqrt();
    let mut dphi = hess(0, 1).abs() * half[0] * half[1] + g.k / 6.0 * d3sum * diag.powi(3);
    if keep {
        dphi += drop_q; // only nonzero where h c² ≤ 1e-4
    } else {
        dphi += 0.5 * (hess(0, 0) * half[0] * half[0] + hess(1, 1) * half[1] * half[1]);
    }
    // amplitude: m = A L f g, f = r⁻², g = R⁻¹
    let rf = rho_max / rmin;
    let (f0, f1, f2) = (1.0 / (rmin * rmin), 2.0 * rf / rmin.powi(3), (2.0 + 6.0 * rf * rf) / rmin.powi(4));
    let (g0, g1, g2) = if g.zs.is_some() {
        let rq = rhos_max / rsmin;
        (1.0 / rsmin, rq / (rsmin * rsmin), (1.0 + 2.0 * rq * rq) / rsmin.powi(3))
    } else {
        (1.0, 0.0, 0.0)
    };
    let a_mag = g.amplitude.norm() * g.l;
    let m_max = a_mag * f0 * g0;
    let m2 = a_mag * (f2 * g0 + 2.0 * f1 * g1 + f0 * g2);
    // |m̃| over the aperture: |m_c (1 + β·Δ)| ≤ m_c(1 + |β||Δ|)
    let m_c = g.amp_exact(p, c).norm();
    let mt_max = m_c * (1.0 + (beta[0].abs() * half[0] + beta[1].abs() * half[1]));
    let damp = 0.5 * m2 * diag * diag + m_max * (1.0 / rmin - 1.0 / rmax) / g.k * 2.0;
    let bound = ap.area() / qr.wavelength * (damp + mt_max * dphi);
    Approx { u, bound }
}

/// Semi-analytic line integral for elongated openings: along the long
/// axis, Gauss–Legendre panels on the *exact* phase (≤ 1 rad per 10-point
/// panel); across the narrow axis, at each node, the local expansion of
/// [`local_expansion`] restricted to that axis (exact phase, gradient,
/// Hessian and linear amplitude in the narrow coordinate, so no cross term
/// is dropped). Returns the field and the rigorous bound on the
/// narrow-axis expansion (third-order phase and second-order amplitude
/// remainders in the narrow coordinate only); the long-axis quadrature
/// error of a 10-point rule on ≤ 1 rad panels is below 1e-13 relative and
/// is checked against the full integral in the tests.
pub fn line_integral(qr: &Query) -> Approx {
    line_impl(qr, true)
}

/// The bound of [`line_integral`] without the quadrature (cheap).
pub fn line_bound(qr: &Query) -> f64 {
    line_impl(qr, false).bound
}

fn line_impl(qr: &Query, want_field: bool) -> Approx {
    let g = Geometry::new(qr);
    let p = qr.point;
    let ap = &qr.aperture;
    let (na, lb) = if ap.size[0] <= ap.size[1] { (0, 1) } else { (1, 0) };
    let hn = 0.5 * ap.size[na];
    // panels along the long axis
    let ro = sq_range(p[lb], ap.lo(lb), ap.hi(lb)).1.sqrt();
    let so = ro / (ro * ro + g.l * g.l).sqrt();
    let si = match g.zs {
        Some(z) => {
            let rs = sq_range(g.xs[lb], ap.lo(lb), ap.hi(lb)).1.sqrt();
            rs / (rs * rs + z * z).sqrt()
        }
        None => g.sines[lb].abs(),
    };
    let panels = if want_field { ((g.k * (so + si) * ap.size[lb]).ceil() as usize).max(1) } else { 0 };
    let (xg, wg) = gl10();
    let mut sum = Cx::new(0.0, 0.0);
    for j in 0..panels {
        let (b0, b1) = (ap.lo(lb) + ap.size[lb] * j as f64 / panels as f64, ap.lo(lb) + ap.size[lb] * (j + 1) as f64 / panels as f64);
        for (t, wt) in xg.iter().zip(&wg) {
            let mut u = ap.center;
            u[lb] = 0.5 * (b0 + b1) + 0.5 * (b1 - b0) * t;
            // narrow-axis expansion at this node
            let w = [u[0] - p[0], u[1] - p[1]];
            let r = (w[0] * w[0] + w[1] * w[1] + g.l * g.l).sqrt();
            let mut grad = -(p[na] - u[na]) / r;
            let mut hess = (1.0 - w[na] * w[na] / (r * r)) / r;
            let mut beta = -2.0 * w[na] / (r * r);
            if let Some(z) = g.zs {
                let ws = [u[0] - g.xs[0], u[1] - g.xs[1]];
                let rs = (ws[0] * ws[0] + ws[1] * ws[1] + z * z).sqrt();
                grad += ws[na] / rs;
                hess += (1.0 - ws[na] * ws[na] / (rs * rs)) / rs;
                beta -= ws[na] / (rs * rs);
            } else {
                grad += g.sines[na];
            }
            let (m0, m1, _) = moments(g.k * grad, g.k * hess, hn, true);
            let val = g.amp_exact(p, u) * Cx::from_polar(1.0, g.phi_exact(p, u)) * (m0 + m1 * beta);
            sum += val * (0.5 * (b1 - b0) * wt);
        }
    }
    let u = g.prefactor(qr.wavelength) * g.carrier() * sum;
    // ---- bound: third-order phase and second-order amplitude in the
    // narrow coordinate (all other dependence is exact).
    let ro2 = rho2_range([p[0], p[1]], ap);
    let rmin = (g.l * g.l + ro2.0).sqrt();
    let wn_max = sq_range(p[na], ap.lo(na), ap.hi(na)).1.sqrt();
    // ∂³r/∂ξ³ = −3w/r³ + 3w³/r⁵ (narrow axis): |·| ≤ 3w/r³ (1 + w²/r²)
    let mut d3 = 3.0 * wn_max / rmin.powi(3) * (1.0 + wn_max * wn_max / (rmin * rmin));
    let (mut g0, mut g1, mut g2) = (1.0, 0.0, 0.0);
    if let Some(z) = g.zs {
        let rr = rho2_range(g.xs, ap);
        let rsmin = (z * z + rr.0).sqrt();
        let wsn = sq_range(g.xs[na], ap.lo(na), ap.hi(na)).1.sqrt();
        d3 += 3.0 * wsn / rsmin.powi(3) * (1.0 + wsn * wsn / (rsmin * rsmin));
        let q = wsn / rsmin;
        g0 = 1.0 / rsmin;
        g1 = q / (rsmin * rsmin);
        g2 = (1.0 + 2.0 * q * q) / rsmin.powi(3);
    }
    let rf = wn_max / rmin;
    let (f0, f1, f2) = (1.0 / (rmin * rmin), 2.0 * rf / rmin.powi(3), (2.0 + 6.0 * rf * rf) / rmin.powi(4));
    let a_mag = g.amplitude.norm() * g.l;
    let m_max = a_mag * f0 * g0;
    let m2 = a_mag * (f2 * g0 + 2.0 * f1 * g1 + f0 * g2);
    let dphi = g.k / 6.0 * d3 * hn.powi(3);
    let beta_max = 2.0 * wn_max / (rmin * rmin) + g1 / g0.max(1e-300);
    let bound = ap.area() / qr.wavelength * (0.5 * m2 * hn * hn + m_max * (1.0 + beta_max * hn) * dphi + 2.0 * m_max * (hn / (rmin * rmin)) / g.k);
    Approx { u, bound }
}

/// Integrand evaluations of [`line_integral`].
fn line_points(qr: &Query) -> f64 {
    let g = Geometry::new(qr);
    let ap = &qr.aperture;
    let lb = if ap.size[0] <= ap.size[1] { 1 } else { 0 };
    let ro = sq_range(qr.point[lb], ap.lo(lb), ap.hi(lb)).1.sqrt();
    let so = ro / (ro * ro + g.l * g.l).sqrt();
    10.0 * ((g.k * (so + 1.0) * ap.size[lb]).ceil() + 1.0)
}

struct LineRung;

impl Rung<Query, Field> for LineRung {
    fn name(&self) -> &str {
        "line integral"
    }
    fn cost(&self, q: &Query) -> f64 {
        150.0 * line_points(q)
    }
    fn assess(&self, q: &Query) -> Assessment {
        if !scalar_ok(q) {
            return Assessment::OutOfDomain;
        }
        Assessment::Claim(Claim::validated(line_bound(q)))
    }
    fn eval(&self, q: &Query) -> Field {
        Field { u: line_integral(q).u, estimate: 0.0 }
    }
}

/// Rayleigh–Sommerfeld integral by Gauss–Legendre panels; `refine`
/// multiplies the panel count. Returns the field.
pub fn rayleigh_sommerfeld(qr: &Query, refine: usize) -> Cx {
    let g = Geometry::new(qr);
    let p = qr.point;
    let ap = &qr.aperture;
    // Panels so that the phase varies by ≲ 2 rad per panel: |∂Φ/∂ξ| ≤ k(|sinθ_o| + |sinθ_s|) ≤ 2k.
    let panels: [usize; 2] = std::array::from_fn(|a| {
        let ro = sq_range(p[a], ap.lo(a), ap.hi(a)).1.sqrt();
        let so = ro / (ro * ro + g.l * g.l).sqrt();
        let si = match g.zs {
            Some(z) => {
                let rs = sq_range(g.xs[a], ap.lo(a), ap.hi(a)).1.sqrt();
                rs / (rs * rs + z * z).sqrt()
            }
            None => g.sines[a].abs(),
        };
        refine * (((g.k * (so + si) * ap.size[a] / 2.0).ceil() as usize).max(1) + 1)
    });
    let (xg, wg) = gl10();
    let mut sum = Cx::new(0.0, 0.0);
    for i in 0..panels[0] {
        let (a0, a1) = (ap.lo(0) + ap.size[0] * i as f64 / panels[0] as f64, ap.lo(0) + ap.size[0] * (i + 1) as f64 / panels[0] as f64);
        for (tx, wx) in xg.iter().zip(&wg) {
            let xi = 0.5 * (a0 + a1) + 0.5 * (a1 - a0) * tx;
            let wxx = 0.5 * (a1 - a0) * wx;
            for j in 0..panels[1] {
                let (b0, b1) = (ap.lo(1) + ap.size[1] * j as f64 / panels[1] as f64, ap.lo(1) + ap.size[1] * (j + 1) as f64 / panels[1] as f64);
                for (ty, wy) in xg.iter().zip(&wg) {
                    let eta = 0.5 * (b0 + b1) + 0.5 * (b1 - b0) * ty;
                    let w = wxx * 0.5 * (b1 - b0) * wy;
                    let u = [xi, eta];
                    sum += g.amp_exact(p, u) * Cx::from_polar(w, g.phi_exact(p, u));
                }
            }
        }
    }
    g.prefactor(qr.wavelength) * g.carrier() * sum
}

fn gl10() -> ([f64; 10], [f64; 10]) {
    (
        [-0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.1488743389816312, 0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717],
        [0.0666713443086881, 0.1494513491505806, 0.2190863625159820, 0.2692667193099963, 0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2190863625159820, 0.1494513491505806, 0.0666713443086881],
    )
}

/// Output of the ladder: the field and, for the quadrature rung, its a
/// posteriori estimate.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Field {
    pub u: Cx,
    pub estimate: f64,
}

pub fn metric(a: &Field, b: &Field) -> f64 {
    (a.u - b.u).norm()
}

/// Estimated cost of the quadrature rung (integrand evaluations).
fn rs_points(qr: &Query, refine: usize) -> f64 {
    let g = Geometry::new(qr);
    let p = qr.point;
    let ap = &qr.aperture;
    let n: f64 = (0..2)
        .map(|a| {
            let ro = sq_range(p[a], ap.lo(a), ap.hi(a)).1.sqrt();
            let so = ro / (ro * ro + g.l * g.l).sqrt();
            10.0 * refine as f64 * ((g.k * (so + 1.0) * ap.size[a] / 2.0).ceil() + 1.0)
        })
        .product();
    n
}

struct FraunhoferRung;
struct FresnelRung;
struct RsRung;

/// Scalar (Kirchhoff) model domain: openings at least this many
/// wavelengths wide.
pub const MIN_WIDTH_IN_WAVELENGTHS: f64 = 10.0;

fn scalar_ok(q: &Query) -> bool {
    q.aperture.scalar_validity(q.wavelength) >= MIN_WIDTH_IN_WAVELENGTHS && q.point[2] > 0.0
}

impl Rung<Query, Field> for FraunhoferRung {
    fn name(&self) -> &str {
        "tilted fraunhofer"
    }
    fn cost(&self, _: &Query) -> f64 {
        60.0
    }
    fn assess(&self, q: &Query) -> Assessment {
        if !scalar_ok(q) {
            return Assessment::OutOfDomain;
        }
        Assessment::Claim(Claim::rigorous(local_bound(q, Order::Linear)))
    }
    fn eval(&self, q: &Query) -> Field {
        Field { u: local_expansion(q, Order::Linear).u, estimate: 0.0 }
    }
}

impl Rung<Query, Field> for FresnelRung {
    fn name(&self) -> &str {
        "tilted fresnel"
    }
    fn cost(&self, _: &Query) -> f64 {
        250.0
    }
    fn assess(&self, q: &Query) -> Assessment {
        if !scalar_ok(q) {
            return Assessment::OutOfDomain;
        }
        Assessment::Claim(Claim::rigorous(local_bound(q, Order::Quadratic)))
    }
    fn eval(&self, q: &Query) -> Field {
        Field { u: local_expansion(q, Order::Quadratic).u, estimate: 0.0 }
    }
}

impl Rung<Query, Field> for RsRung {
    fn name(&self) -> &str {
        "rayleigh-sommerfeld"
    }
    fn cost(&self, q: &Query) -> f64 {
        40.0 * (rs_points(q, 1) + rs_points(q, 2))
    }
    fn assess(&self, q: &Query) -> Assessment {
        if !scalar_ok(q) {
            return Assessment::OutOfDomain;
        }
        Assessment::Deferred
    }
    fn eval(&self, q: &Query) -> Field {
        let a = rayleigh_sommerfeld(q, 1);
        let b = rayleigh_sommerfeld(q, 2);
        Field { u: b, estimate: (a - b).norm() + 1e-15 * b.norm() }
    }
    fn certify(&self, _: &Query, out: &Field) -> Option<Claim> {
        Some(Claim::asymptotic(out.estimate))
    }
}

/// The diffraction ladder (reference: Rayleigh–Sommerfeld).
pub fn ladder(audit_every: u64) -> Ladder<Query, Field> {
    let mut l = Ladder::new("rectangle diffraction", metric);
    l.push(FraunhoferRung);
    l.push(FresnelRung);
    l.push(LineRung);
    let r = l.push(RsRung);
    l.reference(r).audit(if audit_every == 0 { Audit::Off } else { Audit::Every(audit_every) })
}
