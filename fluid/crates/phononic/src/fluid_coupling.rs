//! Oblique sound transmission through an elastic plate immersed in a fluid
//! (fluid–solid–fluid), with both P and SV waves in the plate.
//!
//! At oblique incidence a fluid cannot support shear, so the plate is loaded
//! only by a normal traction and the tangential traction must vanish on both
//! faces. Inside the plate the displacement is the sum of a compressional (P)
//! and a shear-vertical (SV) potential,
//!
//! ```text
//!   u = grad(phi) + curl(0, psi, 0),   phi, psi ~ exp(i(k_x x +- k_z z)),
//!   k_zP = sqrt(k_P^2 - k_x^2),  k_zS = sqrt(k_S^2 - k_x^2),
//! ```
//!
//! and the four stress/displacement components `[u_x, u_z, sigma_zz, sigma_xz]`
//! are propagated across the plate by a 4x4 transfer matrix. The fluid enters
//! through the impedance of a plane wave for the normal velocity,
//! `Z_f = rho_f c_f / cos(theta)`, and the conditions
//! `sigma_zz = -p`, `sigma_xz = 0`, with `p + Z_f v_z` fixed by the incident
//! wave at the front face and `p - Z_f v_z = 0` at the back face.
//!
//! The transmission shows the plate's **coincidence** and **Lamb-wave**
//! resonances: sharp peaks where the trace speed `c_f / sin(theta)` matches a
//! guided mode of the plate. At normal incidence the shear motion decouples and
//! the result must reduce exactly to the longitudinal stack of
//! [`crate::layered`]; that limit and energy conservation are the benchmarks.
//!
//! # Two solvers
//!
//! [`solve`] uses an unconditionally stable **global (stiffness) matrix**
//! formulation: the four partial-wave amplitudes are normalised at the face
//! they originate from, so every exponential factor is bounded by one
//! (`|exp(i k_z d)| <= 1`). This is the single-layer instance of the
//! Rokhlin-Wang global/stiffness-matrix method and conserves energy to
//! round-off even where `|Im k_zS| d` is large (past the shear critical angle
//! at high `f d`). [`solve_transfer_matrix`] is the historical 4x4 propagator,
//! kept for comparison; it is exact in the well-conditioned regime but mixes
//! growing and decaying evanescent partial waves and loses precision once
//! `|Im k_zS| d` is large.
use crate::complex::{solve as solve_linear, solve_particular, C};
use crate::layered::Stack;
use crate::materials::{Fluid, Solid};
use std::f64::consts::PI;

/// Solution of the fluid–plate–fluid problem at one frequency and angle.
#[derive(Clone, Copy, Debug)]
pub struct PlateSolution {
    /// Pressure reflection coefficient.
    pub r: C,
    /// Pressure transmission coefficient.
    pub t: C,
    pub reflectance: f64,
    /// Power transmittance (same fluid on both sides).
    pub transmittance: f64,
    /// `1 - R - T`; zero for a lossless plate.
    pub absorptance: f64,
    pub theta: f64,
    pub omega: f64,
}

impl PlateSolution {
    /// `R + T`; equal to one for a lossless plate.
    pub fn energy(&self) -> f64 {
        self.reflectance + self.transmittance
    }
}

/// Complex vertical wavenumber `sqrt(k^2 - k_x^2)`, principal branch, so an
/// evanescent wave decays along `+z`.
pub fn vertical_wavenumber(k: f64, kx: f64) -> C {
    (C::real(k * k - kx * kx)).sqrt()
}

/// The solid's 4x4 fundamental solution matrix at depth `z`, columns
/// `[A, B, C, D]` (down/up P, down/up SV) and rows `[u_x, u_z, sigma_zz,
/// sigma_xz]`. Row-major.
#[allow(clippy::too_many_arguments)]
fn state_matrix(z: f64, kx: f64, alpha: C, beta: C, lam: f64, mu: f64, kp2: f64) -> [C; 16] {
    let ea = (C::i() * alpha * z).exp();
    let eb = (C::i() * (-alpha) * z).exp();
    let ec = (C::i() * beta * z).exp();
    let ed = (C::i() * (-beta) * z).exp();
    let ikx = C::i() * kx;
    let ia = C::i() * alpha;
    let ib = C::i() * beta;
    let szz_p = -(C::real(lam * kp2) + alpha * alpha * (2.0 * mu));
    let sxp = C::real(-2.0 * mu * kx) * alpha;
    let sxv = (beta * beta - C::real(kx * kx)) * mu;
    let mut m = [C::real(0.0); 16];
    // Column A: down P.
    m[0] = ikx * ea;
    m[4] = ia * ea;
    m[8] = szz_p * ea;
    m[12] = sxp * ea;
    // Column B: up P.
    m[1] = ikx * eb;
    m[5] = -ia * eb;
    m[9] = szz_p * eb;
    m[13] = -sxp * eb;
    // Column C: down SV.
    m[2] = -ib * ec;
    m[6] = ikx * ec;
    m[10] = C::real(-2.0 * mu * kx) * beta * ec;
    m[14] = sxv * ec;
    // Column D: up SV.
    m[3] = ib * ed;
    m[7] = ikx * ed;
    m[11] = C::real(2.0 * mu * kx) * beta * ed;
    m[15] = sxv * ed;
    m
}

/// The solid's four partial-wave base vectors at `z = 0`, rows
/// `[u_x, u_z, sigma_zz, sigma_xz]`, columns `[A, B, C, D]` (down/up P,
/// down/up SV). Row-major. Multiplying column `B` by `qa = exp(i alpha d)`
/// gives the up-P state at `z = d`, and likewise for the others, which is the
/// scaling that keeps the global matrix unconditionally stable.
fn base_vectors(kx: f64, alpha: C, beta: C, lam: f64, mu: f64, kp2: f64) -> [C; 16] {
    let ikx = C::i() * kx;
    let ia = C::i() * alpha;
    let ib = C::i() * beta;
    let szz_p = -(C::real(lam * kp2) + alpha * alpha * (2.0 * mu));
    let sxp = C::real(-2.0 * mu * kx) * alpha;
    let sxv = (beta * beta - C::real(kx * kx)) * mu;
    let two_mu_kx_beta = C::real(2.0 * mu * kx) * beta;
    let mut m = [C::real(0.0); 16];
    // Column A: down P at z = 0.
    m[0] = ikx;
    m[4] = ia;
    m[8] = szz_p;
    m[12] = sxp;
    // Column B: up P measured at z = 0 (base, before the qa factor).
    m[1] = ikx;
    m[5] = -ia;
    m[9] = szz_p;
    m[13] = -sxp;
    // Column C: down SV at z = 0.
    m[2] = -ib;
    m[6] = ikx;
    m[10] = -two_mu_kx_beta;
    m[14] = sxv;
    // Column D: up SV measured at z = d (base, before the qb factor).
    m[3] = ib;
    m[7] = ikx;
    m[11] = two_mu_kx_beta;
    m[15] = sxv;
    m
}

/// Scaled state matrices at the two faces. Columns are ordered
/// `[pf, pb, sf, sb]`, where `pf`, `sf` are the down-P / down-SV amplitudes at
/// `z = 0` and `pb`, `sb` are the up-P / up-SV amplitudes at `z = d`. Every
/// exponential factor is `qa` or `qb`, both bounded by one.
fn scaled_states(
    kx: f64,
    alpha: C,
    beta: C,
    lam: f64,
    mu: f64,
    kp2: f64,
    qa: C,
    qb: C,
) -> ([C; 16], [C; 16]) {
    let base = base_vectors(kx, alpha, beta, lam, mu, kp2);
    let col = |m: &[C; 16], c: usize| [m[c], m[4 + c], m[8 + c], m[12 + c]];
    let (a0, b0, c0, d0) = (col(&base, 0), col(&base, 1), col(&base, 2), col(&base, 3));
    let mut m0 = [C::real(0.0); 16];
    let mut md = [C::real(0.0); 16];
    let set = |m: &mut [C; 16], c: usize, v: [C; 4]| {
        for r in 0..4 {
            m[r * 4 + c] = v[r];
        }
    };
    // At z = 0: down waves base, up waves carry their origin factor.
    set(&mut m0, 0, a0);
    set(&mut m0, 1, [qa * b0[0], qa * b0[1], qa * b0[2], qa * b0[3]]);
    set(&mut m0, 2, c0);
    set(&mut m0, 3, [qb * d0[0], qb * d0[1], qb * d0[2], qb * d0[3]]);
    // At z = d: up waves base, down waves carry their origin factor.
    set(&mut md, 0, [qa * a0[0], qa * a0[1], qa * a0[2], qa * a0[3]]);
    set(&mut md, 1, b0);
    set(&mut md, 2, [qb * c0[0], qb * c0[1], qb * c0[2], qb * c0[3]]);
    set(&mut md, 3, d0);
    (m0, md)
}

/// Unconditionally stable global (stiffness) matrix solver. The four
/// partial-wave amplitudes are the unknowns; the boundary conditions are the
/// vanishing tangential traction on both faces and the fluid impedance
/// conditions at each face. Because each unknown is normalised at its face of
/// origin, no factor `exp(|Im k_z| d)` ever appears and the system stays well
/// conditioned past the shear critical angle.
pub fn solve_stable(fluid: &Fluid, solid: &Solid, d: f64, omega: f64, theta: f64) -> PlateSolution {
    assert!(d > 0.0 && omega > 0.0 && theta.is_finite());
    let kx = (omega / fluid.c) * theta.sin();
    let cos_t = theta.cos();
    let zf = fluid.rho * fluid.c / cos_t;
    let (lam, mu) = (solid.lambda(), solid.mu());
    let (kp, ks) = (omega / solid.c_p, omega / solid.c_s);
    let alpha = vertical_wavenumber(kp, kx);
    let beta = vertical_wavenumber(ks, kx);
    let qa = (C::i() * alpha * d).exp();
    let qb = (C::i() * beta * d).exp();
    let (m0, md) = scaled_states(kx, alpha, beta, lam, mu, kp * kp, qa, qb);
    let iwzf = C::i() * (omega * zf);
    // Equations, in order:
    //   0: sigma_xz(0) = 0
    //   1: sigma_zz(0) - i omega Zf u_z(0) = -2
    //   2: sigma_xz(d) = 0
    //   3: sigma_zz(d) + i omega Zf u_z(d) = 0
    let mut a = [C::real(0.0); 16];
    let mut rhs = [C::real(0.0); 4];
    for c in 0..4 {
        a[c] = m0[12 + c];
        a[4 + c] = m0[8 + c] - iwzf * m0[4 + c];
        a[8 + c] = md[12 + c];
        a[12 + c] = md[8 + c] + iwzf * md[4 + c];
    }
    rhs[1] = C::real(-2.0);
    // Row equilibration: the equations mix displacements (~1) and stresses
    // (~1e13), so scale each row before elimination.
    for r in 0..4 {
        let mut mx = 0.0f64;
        for c in 0..4 {
            mx = mx.max(a[r * 4 + c].abs());
        }
        if mx > 0.0 {
            for c in 0..4 {
                a[r * 4 + c] = a[r * 4 + c] / mx;
            }
            rhs[r] = rhs[r] / mx;
        }
    }
    let x = solve_particular(&a, &rhs, 4);
    let mut szz0 = C::real(0.0);
    let mut szzd = C::real(0.0);
    for c in 0..4 {
        szz0 = szz0 + m0[8 + c] * x[c];
        szzd = szzd + md[8 + c] * x[c];
    }
    let r = -szz0 - C::real(1.0);
    let t = -szzd;
    PlateSolution {
        r,
        t,
        reflectance: r.abs2(),
        transmittance: t.abs2(),
        absorptance: 1.0 - r.abs2() - t.abs2(),
        theta,
        omega,
    }
}

/// Solve the fluid–plate–fluid problem for a plate of thickness `d` at
/// angular frequency `omega` and incidence angle `theta`. Uses the
/// unconditionally stable global matrix formulation.
pub fn solve(fluid: &Fluid, solid: &Solid, d: f64, omega: f64, theta: f64) -> PlateSolution {
    solve_stable(fluid, solid, d, omega, theta)
}

/// Determinant of the free-plate boundary system at plate phase speed `c`:
/// rows `[sigma_zz(0), sigma_xz(0), sigma_zz(d), sigma_xz(d)]`, columns the
/// four scaled partial-wave amplitudes. It vanishes exactly on a
/// Rayleigh–Lamb mode and is the independent check for [`crate::lamb`].
/// The same scaled (bounded-exponential) matrices are used, so it is stable at
/// any `f d`.
pub fn free_plate_determinant(solid: &Solid, d: f64, omega: f64, c: f64) -> C {
    let kx = omega / c;
    let (lam, mu) = (solid.lambda(), solid.mu());
    let (kp, ks) = (omega / solid.c_p, omega / solid.c_s);
    let alpha = vertical_wavenumber(kp, kx);
    let beta = vertical_wavenumber(ks, kx);
    let qa = (C::i() * alpha * d).exp();
    let qb = (C::i() * beta * d).exp();
    let (m0, md) = scaled_states(kx, alpha, beta, lam, mu, kp * kp, qa, qb);
    let mut m = [C::real(0.0); 16];
    for c in 0..4 {
        m[c] = m0[8 + c];
        m[4 + c] = m0[12 + c];
        m[8 + c] = md[8 + c];
        m[12 + c] = md[12 + c];
    }
    det4(&m)
}

/// 4x4 complex determinant (row-major).
fn det4(a: &[C; 16]) -> C {
    let s = |i: usize, j: usize| a[i * 4 + j];
    let m0 = s(0, 0)
        * (s(1, 1) * (s(2, 2) * s(3, 3) - s(2, 3) * s(3, 2))
            - s(1, 2) * (s(2, 1) * s(3, 3) - s(2, 3) * s(3, 1))
            + s(1, 3) * (s(2, 1) * s(3, 2) - s(2, 2) * s(3, 1)));
    let m1 = s(0, 1)
        * (s(1, 0) * (s(2, 2) * s(3, 3) - s(2, 3) * s(3, 2))
            - s(1, 2) * (s(2, 0) * s(3, 3) - s(2, 3) * s(3, 0))
            + s(1, 3) * (s(2, 0) * s(3, 2) - s(2, 2) * s(3, 0)));
    let m2 = s(0, 2)
        * (s(1, 0) * (s(2, 1) * s(3, 3) - s(2, 3) * s(3, 1))
            - s(1, 1) * (s(2, 0) * s(3, 3) - s(2, 3) * s(3, 0))
            + s(1, 3) * (s(2, 0) * s(3, 1) - s(2, 1) * s(3, 0)));
    let m3 = s(0, 3)
        * (s(1, 0) * (s(2, 1) * s(3, 2) - s(2, 2) * s(3, 1))
            - s(1, 1) * (s(2, 0) * s(3, 2) - s(2, 2) * s(3, 0))
            + s(1, 2) * (s(2, 0) * s(3, 1) - s(2, 1) * s(3, 0)));
    m0 - m1 + m2 - m3
}

/// Historical 4x4 transfer-matrix solver, kept for comparison. It is exact in
/// the well-conditioned regime (small `|Im k_zS| d`) but its propagator mixes
/// growing and decaying evanescent partial waves, so it loses precision past
/// the shear critical angle once `|Im k_zS| d` is large.
pub fn solve_transfer_matrix(
    fluid: &Fluid,
    solid: &Solid,
    d: f64,
    omega: f64,
    theta: f64,
) -> PlateSolution {
    assert!(d > 0.0 && omega > 0.0 && theta.is_finite());
    let kx = (omega / fluid.c) * theta.sin();
    // Normal-velocity impedance of the fluid half-spaces.
    let cos_t = theta.cos();
    let zf = fluid.rho * fluid.c / cos_t;
    let (lam, mu) = (solid.lambda(), solid.mu());
    let (kp, ks) = (omega / solid.c_p, omega / solid.c_s);
    let alpha = vertical_wavenumber(kp, kx);
    let beta = vertical_wavenumber(ks, kx);

    let e0 = state_matrix(0.0, kx, alpha, beta, lam, mu, kp * kp);
    let ed = state_matrix(d, kx, alpha, beta, lam, mu, kp * kp);
    // Propagator P = E(d) E(0)^-1: invert E(0) column by column, then
    // multiply by E(d). (Solving E(0) x = E(d)[:,j] would give
    // E(0)^-1 E(d), the wrong direction.)
    let mut inv = [C::real(0.0); 16];
    for col in 0..4 {
        let mut rhs = [C::real(0.0); 4];
        rhs[col] = C::real(1.0);
        let x = solve_linear(&e0, &rhs, 4).expect("ideal plate propagator is singular");
        for row in 0..4 {
            inv[row * 4 + col] = x[row];
        }
    }
    let mut b = [C::real(0.0); 16];
    for row in 0..4 {
        for col in 0..4 {
            let mut acc = C::real(0.0);
            for k in 0..4 {
                acc = acc + ed[row * 4 + k] * inv[k * 4 + col];
            }
            b[row * 4 + col] = acc;
        }
    }

    // Unknowns x = [u_x(0), u_z(0), sigma_zz(0)], with sigma_xz(0) = 0.
    let iwzf = C::i() * (omega * zf);
    let a33 = [
        C::real(0.0),
        -iwzf,
        C::real(1.0),
        b[12],
        b[13],
        b[14],
        b[8] + iwzf * b[4],
        b[9] + iwzf * b[5],
        b[10] + iwzf * b[6],
    ];
    let rhs = [C::real(-2.0), C::real(0.0), C::real(0.0)];
    let x = solve_particular(&a33, &rhs, 3);
    let szz0 = x[2];
    let r = -szz0 - C::real(1.0);
    // sigma_zz(d) = row 2 of B times s0 (s0[3] = 0).
    let szzd = b[8] * x[0] + b[9] * x[1] + b[10] * x[2];
    let t = -szzd;
    let reflectance = r.abs2();
    let transmittance = t.abs2();
    PlateSolution {
        r,
        t,
        reflectance,
        transmittance,
        absorptance: 1.0 - reflectance - transmittance,
        theta,
        omega,
    }
}

/// Transmittance only.
pub fn transmittance(fluid: &Fluid, solid: &Solid, d: f64, omega: f64, theta: f64) -> f64 {
    solve(fluid, solid, d, omega, theta).transmittance
}

/// Trace (phase) speed along the plate, `c_f / sin(theta)`.
pub fn trace_speed(fluid: &Fluid, theta: f64) -> f64 {
    fluid.c / theta.sin()
}

/// Coincidence angle where the trace speed equals the plate's compressional
/// speed, `asin(c_f / c_p)`; `None` when the fluid is faster than the plate.
pub fn coincidence_angle(fluid: &Fluid, solid: &Solid) -> Option<f64> {
    coincidence_angle_for_speed(fluid, solid.c_p)
}

/// Coincidence angle for an arbitrary guided-plate speed,
/// `asin(c_f / c_plate)`; `None` when the plate is slower than the fluid.
pub fn coincidence_angle_for_speed(fluid: &Fluid, speed: f64) -> Option<f64> {
    let s = fluid.c / speed;
    if s <= 1.0 {
        Some(s.asin())
    } else {
        None
    }
}

/// Normal-incidence reference through the longitudinal stack.
pub fn normal_incidence_reference(fluid: &Fluid, solid: &Solid, d: f64, omega: f64) -> f64 {
    let layer = crate::layered::Layer::from_solid(solid, d);
    Stack::from_media(fluid, vec![layer], fluid)
        .solve(omega)
        .transmittance
}

/// Extensional (S0) plate speed in the thin-plate limit,
/// `sqrt(E / (rho (1 - nu^2)))`. This is the guided speed whose coincidence
/// with the fluid trace produces the strong S0 transmission peak.
pub fn plate_speed_s0(solid: &Solid) -> f64 {
    let nu = (solid.lambda() / (2.0 * (solid.lambda() + solid.mu()))).clamp(-0.99, 0.49);
    (solid.young() / (solid.rho * (1.0 - nu * nu))).sqrt()
}

/// Number of half-wavelengths of a compressional wave in the plate; used to
/// label thickness resonances.
pub fn plate_halfwaves(solid: &Solid, f: f64, d: f64) -> f64 {
    2.0 * f * d / solid.c_p
}

/// Design frequency of the `n`-th thickness resonance `f = n c_p / (2 d)`.
pub fn thickness_resonance(solid: &Solid, d: f64, n: u32) -> f64 {
    n as f64 * solid.c_p / (2.0 * d)
}

/// Two times pi, for callers building angular frequencies.
pub const TWO_PI: f64 = 2.0 * PI;
