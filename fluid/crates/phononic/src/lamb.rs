//! Free-plate Rayleigh–Lamb dispersion, and the coincidence angles at which
//! its leaky modes couple to a fluid load.
//!
//! For an isotropic, traction-free plate of thickness `d` the displacement is
//! the same P + SV potential sum used by [`crate::fluid_coupling`]. The two
//! traction-free faces give a 4x4 homogeneous system whose determinant
//! vanishes on the guided modes. Splitting by mirror symmetry yields the
//! classical Rayleigh–Lamb equations (Graff, *Wave Motion in Elastic Solids*,
//! ch. 6; Rose, *Ultrasonic Waves in Solid Media*, ch. 8). With
//! `alpha^2 = omega^2/c_p^2 - k^2`, `beta^2 = omega^2/c_s^2 - k^2`,
//! `delta = k^2 - beta^2`, `h = d/2`, the pole-free forms used here are
//!
//! ```text
//!   symmetric (S):     delta^2 sin(beta h) cos(alpha h)
//!                    + 4 alpha beta k^2 sin(alpha h) cos(beta h) = 0,
//!   antisymmetric (A): 4 alpha beta k^2 sin(beta h) cos(alpha h)
//!                    + delta^2 sin(alpha h) cos(beta h) = 0.
//! ```
//!
//! These are equivalent to the usual `tan` ratio forms but have no poles.
//! Because `alpha` and `beta` turn imaginary below the P and S speeds, the
//! residual is evaluated with the real functions appropriate to each of the
//! three speed intervals (`c > c_p`, `c_s < c < c_p`, `c < c_s`); the explicit
//! hyperbolic forms avoid returnung the wrong (identically zero) component of
//! a complex expression. The cross-multiplication introduces one non-physical
//! root wherever `sin(beta h) = 0` (the `beta -> 0` degeneracy at `c_s`) or
//! `sin(alpha h) = 0` on the antisymmetric branch at `c_p`; callers should
//! skip speeds within a small tolerance of the two body-wave speeds.
//!
//! The two branches recover the thin-plate limits: `S0 -> sqrt(E/(rho(1-nu^2)))`
//! and `A0 -> 0` (flexural) as `omega d -> 0`.
//!
//! A mode with phase speed `c > c_f` has a coincidence angle
//! `theta = asin(c_f/c)`; the immersed-plate transmission peaks sit at these
//! angles, shifted slightly by the fluid loading.
use crate::materials::Solid;
use std::f64::consts::PI;

/// Symmetric (extensional) or antisymmetric (flexural) Rayleigh–Lamb branch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LambKind {
    /// Symmetric modes: `S0, S1, ...` (extensional).
    Symmetric,
    /// Antisymmetric modes: `A0, A1, ...` (flexural).
    Antisymmetric,
}

/// One Rayleigh–Lamb mode at a given frequency: phase speed and branch.
#[derive(Clone, Copy, Debug)]
pub struct LambMode {
    pub speed: f64,
    pub kind: LambKind,
}

/// Pole-free real Rayleigh–Lamb residual at phase speed `c`; zero on a mode.
///
/// `alpha` and `beta` are continued below their body-wave speeds by replacing
/// `sin(i x) -> i sinh(x)` (and likewise the cosines), which is what makes the
/// product `4 alpha beta k^2 sin(alpha h) cos(beta h)` real in every interval.
pub fn residual(solid: &Solid, d: f64, omega: f64, c: f64, kind: LambKind) -> f64 {
    let k = omega / c;
    let k2 = k * k;
    let h = 0.5 * d;
    let kp2 = (omega / solid.c_p).powi(2);
    let ks2 = (omega / solid.c_s).powi(2);
    let a2 = kp2 - k2;
    let b2 = ks2 - k2;
    // `delta = k^2 - beta^2` is `2 k^2 - ks^2` in every interval.
    let delta = 2.0 * k2 - ks2;
    let d2 = delta * delta;
    if b2 >= 0.0 {
        // Both body waves propagating (c >= c_s). `alpha` real when c >= c_p.
        let (sa, ca) = if a2 >= 0.0 {
            let al = a2.sqrt();
            ((al * h).sin(), (al * h).cos())
        } else {
            let al = (-a2).sqrt();
            ((al * h).sinh(), (al * h).cosh())
        };
        let be = b2.sqrt();
        let (sb, cb) = ((be * h).sin(), (be * h).cos());
        // `4 alpha beta k^2` becomes `-4|alpha| beta k^2` when alpha is
        // imaginary, but the sign differs between the two branches: the `i`
        // from `sin(alpha h)` sits in the second term of S and the first
        // term of A, so only S picks up the minus sign in this interval.
        let (four_s, four_a) = if a2 >= 0.0 {
            let m = 4.0 * a2.sqrt() * be * k2;
            (m, m)
        } else {
            let m = 4.0 * (-a2).sqrt() * be * k2;
            (-m, m)
        };
        match kind {
            LambKind::Symmetric => d2 * sb * ca + four_s * sa * cb,
            LambKind::Antisymmetric => four_a * sb * ca + d2 * sa * cb,
        }
    } else {
        // Both body waves evanescent (c < c_s); all functions hyperbolic and
        // `4 alpha beta k^2 = -4 |alpha| |beta| k^2` is real.
        let al = (-a2).sqrt();
        let be = (-b2).sqrt();
        let (sa, ca) = ((al * h).sinh(), (al * h).cosh());
        let (sb, cb) = ((be * h).sinh(), (be * h).cosh());
        let four = -4.0 * al * be * k2;
        match kind {
            LambKind::Symmetric => d2 * sb * ca + four * sa * cb,
            LambKind::Antisymmetric => four * sb * ca + d2 * sa * cb,
        }
    }
}

/// Find all Rayleigh–Lamb modes with phase speed in `[c_lo, c_hi]` at angular
/// frequency `omega`, for one branch. The residual is smooth there, so a
/// sign-change scan plus bisection is robust. Speeds within `1e-6` relative of
/// `c_s` or `c_p` are skipped: those are the removable roots produced by
/// cross-multiplying the `tan` ratio equations.
pub fn modes_branch(
    solid: &Solid,
    d: f64,
    omega: f64,
    c_lo: f64,
    c_hi: f64,
    kind: LambKind,
    steps: usize,
) -> Vec<f64> {
    assert!(c_lo > 0.0 && c_hi > c_lo && steps > 0);
    let f = |c: f64| residual(solid, d, omega, c, kind);
    let skip = |c: f64| {
        let tol = 1e-6;
        (c - solid.c_s).abs() < tol * solid.c_s || (c - solid.c_p).abs() < tol * solid.c_p
    };
    let mut out = Vec::new();
    let mut c_prev = c_lo;
    let mut v_prev = f(c_prev);
    for i in 1..=steps {
        let c = c_lo + (c_hi - c_lo) * i as f64 / steps as f64;
        let v = f(c);
        if v_prev == 0.0 {
            if !skip(c_prev) {
                out.push(c_prev);
            }
        } else if v_prev * v < 0.0 {
            let (mut a, mut b) = (c_prev, c);
            let (mut va, mut _vb) = (v_prev, v);
            for _ in 0..200 {
                let m = 0.5 * (a + b);
                let vm = f(m);
                if va * vm <= 0.0 {
                    b = m;
                } else {
                    a = m;
                    va = vm;
                }
                if b - a < 1e-12 * (1.0 + 0.5 * (a + b)) {
                    break;
                }
            }
            let root = 0.5 * (a + b);
            if !skip(root) {
                out.push(root);
            }
        }
        c_prev = c;
        v_prev = v;
    }
    out
}

/// All Rayleigh–Lamb modes with phase speed in `[c_lo, c_hi]`, both branches,
/// sorted by speed.
pub fn modes(
    solid: &Solid,
    d: f64,
    omega: f64,
    c_lo: f64,
    c_hi: f64,
    steps: usize,
) -> Vec<LambMode> {
    let mut out: Vec<LambMode> = Vec::new();
    for kind in [LambKind::Symmetric, LambKind::Antisymmetric] {
        for speed in modes_branch(solid, d, omega, c_lo, c_hi, kind, steps) {
            out.push(LambMode { speed, kind });
        }
    }
    out.sort_by(|a, b| a.speed.partial_cmp(&b.speed).unwrap());
    out
}

/// Coincidence angle `asin(c_f / c)` for a guided mode with phase speed `c`;
/// `None` when the mode is slower than the fluid (no coupling angle).
pub fn coincidence_angle(c_fluid: f64, mode_speed: f64) -> Option<f64> {
    let s = c_fluid / mode_speed;
    if s <= 1.0 {
        Some(s.asin())
    } else {
        None
    }
}

/// Frequency-thickness product [Hz m].
pub fn fd(f: f64, d: f64) -> f64 {
    f * d
}

/// Two pi, for callers building angular frequencies.
pub const TWO_PI: f64 = 2.0 * PI;
