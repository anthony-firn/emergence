//! The Berreman 4x4 transfer matrix for anisotropic (crystalline and
//! liquid-crystalline) layered media.
//!
//! `multilayer` handles isotropic stacks with a 2x2 (Abelès) matrix per
//! polarisation. When the layers are anisotropic the two polarisations couple,
//! and the correct local model is Berreman's 4x4 matrix (D. W. Berreman,
//! "Optics in stratified and anisotropic media: 4x4-matrix formulation",
//! J. Opt. Soc. Am. **62**, 502 (1972), doi:10.1364/JOSA.62.000502).
//!
//! Maxwell's equations in a uniform layer are written for the tangential
//! field state
//!
//! ```text
//!   psi = (E_x, E_y, H_x, H_y)^T          (H scaled by the vacuum impedance)
//!   d psi / d z = i k0 Delta psi,          k0 = 2 pi / lambda,
//! ```
//!
//! with the conserved tangential wavevector `s = n0 sin(theta0)` and the
//! relative permittivity tensor `eps = eps_ij` (complex, so the form is
//! lossy; `-i g`-type off-diagonal elements make it gyrotropic). The layer
//! matrix is `P = exp(i k0 d Delta)`; the stack is the ordered product. This
//! is exactly Maxwell's equations, so it reduces to `multilayer` for
//! isotropic `eps` and is the reference for every anisotropic case.
//!
//! Reflection and transmission are returned as full complex 2x2 Jones
//! matrices in the tangential-electric-field basis `(E_x, E_y)` (the `p` and
//! `s` amplitudes). Reflectance and transmittance follow from the normal
//! Poynting flux weights `Re(q)` (`s`) and `Re(eps/q)` (`p`), so
//! `R + T = 1` for every lossless stack.
//!
//! # The fold
//!
//! The exact layered model is the fine rung. Two coarse closures are the
//! subject of this module's benchmarks:
//!
//! * **Sub-wavelength periodic medium** (cholesteric pitch `p << lambda`):
//!   the stack homogenises to a uniaxial effective medium (plus a small
//!   optical activity). Indicator `p / lambda`; the reflected power falls as
//!   `O((p/lambda)^2)`. Recover it with [`effective_index_convergence`].
//! * **Slowly-varying director** (Mauguin/adiabatic limit): when the
//!   retardation `pi d (n_e - n_o) / lambda` much exceeds the twist angle
//!   `phi`, the polarisation follows the local director (a waveguide). The
//!   indicator is the inverse Mauguin number `phi lambda / (pi d Delta n)`;
//!   [`twisted_nematic_jones`] gives the exact Jones matrix for comparison.
use crate::complex::C;
use std::f64::consts::PI;

/// A 3x3 complex relative permittivity tensor `eps_ij`.
pub type Tensor3 = [[C; 3]; 3];

/// Identity 3x3 scaled by `eps`.
pub fn isotropic(eps: C) -> Tensor3 {
    let z = C::real(0.0);
    let mut t = [[z; 3]; 3];
    for i in 0..3 {
        t[i][i] = eps;
    }
    t
}

/// Isotropic permittivity from a refractive index `eps = n^2`.
pub fn isotropic_index(n: C) -> Tensor3 {
    isotropic(n * n)
}

/// Uniaxial permittivity with ordinary index `n_o`, extraordinary index
/// `n_e`, and optic axis `axis` (need not be unit; it is normalised):
/// `eps = n_o^2 I + (n_e^2 - n_o^2) a a^T`.
pub fn uniaxial(n_o: C, n_e: C, axis: [f64; 3]) -> Tensor3 {
    let len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    assert!(len > 0.0 && len.is_finite());
    let a = [axis[0] / len, axis[1] / len, axis[2] / len];
    let no2 = n_o * n_o;
    let dn = n_e * n_e - no2;
    let mut t = isotropic(no2);
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = t[i][j] + C::real(a[i] * a[j]) * dn;
        }
    }
    t
}

/// Biaxial permittivity from principal indices, with the principal axes along
/// `x`, `y`, `z` (apply [`rotate`] for another orientation).
pub fn biaxial(n_x: C, n_y: C, n_z: C) -> Tensor3 {
    let z = C::real(0.0);
    let mut t = [[z; 3]; 3];
    t[0][0] = n_x * n_x;
    t[1][1] = n_y * n_y;
    t[2][2] = n_z * n_z;
    t
}

/// Rotate a permittivity tensor, `eps' = R eps R^T`, for a real 3x3 `R`.
pub fn rotate(eps: &Tensor3, r: &[[f64; 3]; 3]) -> Tensor3 {
    // tmp = eps R^T, out = R tmp.
    let mut tmp = [[C::real(0.0); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut acc = C::real(0.0);
            for k in 0..3 {
                acc = acc + eps[i][k] * r[j][k];
            }
            tmp[i][j] = acc;
        }
    }
    let mut out = [[C::real(0.0); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut acc = C::real(0.0);
            for k in 0..3 {
                acc = acc + C::real(r[i][k]) * tmp[k][j];
            }
            out[i][j] = acc;
        }
    }
    out
}

/// Rotation matrix about a unit axis by `angle` (Rodrigues).
pub fn rotation(axis: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
    let len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    let (x, y, z) = (axis[0] / len, axis[1] / len, axis[2] / len);
    let c = angle.cos();
    let s = angle.sin();
    let t = 1.0 - c;
    [
        [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
        [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
        [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
    ]
}

/// Gyrotropic (magneto-optic) permittivity, magnetization along a unit axis:
/// `eps = eps0 I - i g [a]_x` in the plane normal to `a` (the standard
/// Faraday/Voigt form `eps_xy = -i g`, `eps_yx = +i g` for `a = z`).
pub fn gyrotropic(eps0: C, g: C, axis: [f64; 3]) -> Tensor3 {
    let len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    let a = [axis[0] / len, axis[1] / len, axis[2] / len];
    // [a]_x is the cross-product matrix; eps = eps0 I - i g [a]_x.
    let cross = [
        [C::real(0.0), C::real(-a[2]), C::real(a[1])],
        [C::real(a[2]), C::real(0.0), C::real(-a[0])],
        [C::real(-a[1]), C::real(a[0]), C::real(0.0)],
    ];
    let mut t = isotropic(eps0);
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = t[i][j] - C::i() * g * cross[i][j];
        }
    }
    t
}

// ---------------------------------------------------------------------------
// 4x4 complex matrix
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
pub struct Mat4 {
    pub m: [[C; 4]; 4],
}

impl Mat4 {
    pub fn zero() -> Self {
        Self {
            m: [[C::real(0.0); 4]; 4],
        }
    }
    pub fn identity() -> Self {
        let mut a = Self::zero();
        for i in 0..4 {
            a.m[i][i] = C::real(1.0);
        }
        a
    }
    pub fn mul(&self, b: &Self) -> Self {
        let mut out = Self::zero();
        for i in 0..4 {
            for j in 0..4 {
                let mut acc = C::real(0.0);
                for k in 0..4 {
                    acc = acc + self.m[i][k] * b.m[k][j];
                }
                out.m[i][j] = acc;
            }
        }
        out
    }
    pub fn add(&self, b: &Self) -> Self {
        let mut out = Self::zero();
        for i in 0..4 {
            for j in 0..4 {
                out.m[i][j] = self.m[i][j] + b.m[i][j];
            }
        }
        out
    }
    pub fn scale(&self, s: C) -> Self {
        let mut out = Self::zero();
        for i in 0..4 {
            for j in 0..4 {
                out.m[i][j] = self.m[i][j] * s;
            }
        }
        out
    }
    /// Apply the matrix to a 4-vector.
    pub fn apply(&self, v: &[C; 4]) -> [C; 4] {
        let mut out = [C::real(0.0); 4];
        for i in 0..4 {
            let mut acc = C::real(0.0);
            for k in 0..4 {
                acc = acc + self.m[i][k] * v[k];
            }
            out[i] = acc;
        }
        out
    }

    /// Matrix exponential by scaling and squaring with a truncated Taylor
    /// series. The argument here is `A = i k0 d Delta`, whose norm may be
    /// large (thick or strongly birefringent layers), hence the scaling.
    pub fn expm(a: &Mat4) -> Self {
        // Infinity norm (max absolute row sum).
        let mut norm = 0.0f64;
        for i in 0..4 {
            let mut row = 0.0;
            for j in 0..4 {
                row += a.m[i][j].abs();
            }
            norm = norm.max(row);
        }
        let mut s = 0u32;
        let mut scaled = norm;
        while scaled > 0.5 {
            scaled *= 0.5;
            s += 1;
        }
        let factor = 2f64.powi(-(s as i32));
        let b = a.scale(C::real(factor));
        // Taylor series for exp(B).
        let mut term = Self::identity();
        let mut sum = Self::identity();
        for k in 1..=14u32 {
            term = term.mul(&b).scale(C::real(1.0 / k as f64));
            sum = sum.add(&term);
        }
        for _ in 0..s {
            sum = sum.mul(&sum);
        }
        sum
    }
}

/// Solve the 4x4 complex linear system `A x = b` by Gaussian elimination with
/// partial pivoting. Returns `None` for a singular system.
pub fn solve4(a: &[[C; 4]; 4], b: &[C; 4]) -> Option<[C; 4]> {
    let mut m = *a;
    let mut x = *b;
    for col in 0..4 {
        // Pivot.
        let mut piv = col;
        for r in col + 1..4 {
            if m[r][col].abs() > m[piv][col].abs() {
                piv = r;
            }
        }
        if m[piv][col].abs() < 1e-300 {
            return None;
        }
        if piv != col {
            m.swap(piv, col);
            x.swap(piv, col);
        }
        let d = m[col][col];
        for r in col + 1..4 {
            let f = m[r][col] / d;
            for c in col..4 {
                m[r][c] = m[r][c] - f * m[col][c];
            }
            x[r] = x[r] - f * x[col];
        }
    }
    let mut out = [C::real(0.0); 4];
    for i in (0..4).rev() {
        let mut acc = x[i];
        for j in i + 1..4 {
            acc = acc - m[i][j] * out[j];
        }
        out[i] = acc / m[i][i];
    }
    Some(out)
}

// ---------------------------------------------------------------------------
// Berreman layer and stack
// ---------------------------------------------------------------------------

/// The Berreman differential matrix `Delta` for a uniform medium with
/// relative permittivity `eps`, scalar relative permeability `mu`, and
/// conserved tangential wavevector `s`.
pub fn berreman_delta(eps: &Tensor3, mu: C, s: f64) -> Mat4 {
    let e = eps;
    let ezz = e[2][2];
    let sc = C::real(s);
    let s2 = C::real(s * s);
    let mut d = Mat4::zero();
    d.m[0][0] = -(sc * e[2][0]) / ezz;
    d.m[0][1] = -(sc * e[2][1]) / ezz;
    d.m[0][3] = mu - s2 / ezz;
    d.m[1][2] = -mu;
    d.m[2][0] = -(e[1][0] - e[1][2] * e[2][0] / ezz);
    d.m[2][1] = -(e[1][1] - e[1][2] * e[2][1] / ezz - s2 / mu);
    d.m[2][3] = sc * e[1][2] / ezz;
    d.m[3][0] = e[0][0] - e[0][2] * e[2][0] / ezz;
    d.m[3][1] = e[0][1] - e[0][2] * e[2][1] / ezz;
    d.m[3][3] = -(sc * e[0][2]) / ezz;
    d
}

#[derive(Clone, Copy, Debug)]
pub struct AnisoLayer {
    pub eps: Tensor3,
    /// Scalar relative permeability (non-magnetic optical media: 1).
    pub mu: C,
    pub thickness: f64,
}

impl AnisoLayer {
    pub fn new(eps: Tensor3, thickness: f64) -> Self {
        assert!(thickness > 0.0 && thickness.is_finite());
        Self {
            eps,
            mu: C::real(1.0),
            thickness,
        }
    }
    pub fn with_mu(eps: Tensor3, mu: C, thickness: f64) -> Self {
        assert!(thickness > 0.0 && thickness.is_finite());
        Self { eps, mu, thickness }
    }
    /// 4x4 layer transfer matrix at conserved tangential wavevector `s`.
    pub fn matrix(&self, lambda: f64, s: f64) -> Mat4 {
        let k0 = 2.0 * PI / lambda;
        let delta = berreman_delta(&self.eps, self.mu, s);
        Mat4::expm(&delta.scale(C::i() * (k0 * self.thickness)))
    }
}

/// Full 2x2 complex Jones matrix, rows = output `(E_x, E_y)`, columns =
/// input `(E_x, E_y)` in the tangential-electric-field basis.
#[derive(Clone, Copy, Debug)]
pub struct Jones {
    pub m: [[C; 2]; 2],
}

impl Jones {
    pub fn get(&self, i: usize, j: usize) -> C {
        self.m[i][j]
    }
    pub fn apply(&self, v: [C; 2]) -> [C; 2] {
        [
            self.m[0][0] * v[0] + self.m[0][1] * v[1],
            self.m[1][0] * v[0] + self.m[1][1] * v[1],
        ]
    }
}

impl std::ops::Mul for Jones {
    type Output = Jones;
    fn mul(self, b: Jones) -> Jones {
        let mut m = [[C::real(0.0); 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                m[i][j] = self.m[i][0] * b.m[0][j] + self.m[i][1] * b.m[1][j];
            }
        }
        Jones { m }
    }
}

#[derive(Clone, Debug)]
pub struct AnisoSolution {
    pub r: Jones,
    pub t: Jones,
    /// Reflectance/transmittance for a pure `s` input (`E_y`).
    pub reflectance_s: f64,
    pub transmittance_s: f64,
    /// Reflectance/transmittance for a pure `p` input (`E_x`).
    pub reflectance_p: f64,
    pub transmittance_p: f64,
    /// Transmitted tangential field for an input Jones vector.
    pub exit_forward_phase: C,
    /// Incident tangential `s` flux weight `Re(q0)`.
    pub incident: f64,
    /// Exit-medium refractive index (isotropic).
    pub exit_i: f64,
    pub q0: C,
    pub q_out: C,
    pub lambda: f64,
    pub s: f64,
}

impl AnisoSolution {
    pub fn absorptance_s(&self) -> f64 {
        1.0 - self.reflectance_s - self.transmittance_s
    }
    pub fn absorptance_p(&self) -> f64 {
        1.0 - self.reflectance_p - self.transmittance_p
    }
    /// Transmitted tangential field for an input `(E_x, E_y)`.
    pub fn transmitted(&self, input: [C; 2]) -> [C; 2] {
        self.t.apply(input)
    }
    /// Reflected tangential `s`/`p` flux weights in the incident medium.
    fn incident_weights(&self) -> (f64, f64) {
        let w_s = self.q0.re;
        let w_p = (C::real(self.incident * self.incident) / self.q0).re;
        (w_s, w_p)
    }
    /// Reflectance of the stack for an arbitrary incident tangential Jones
    /// vector `(E_x, E_y)`, using the normal Poynting flux weights of the
    /// incident medium.
    pub fn reflectance(&self, input: [C; 2]) -> f64 {
        let r = self.r.apply(input);
        let (w_s, w_p) = self.incident_weights();
        let num = w_s * r[0].abs2() + w_p * r[1].abs2();
        let den = w_s * input[0].abs2() + w_p * input[1].abs2();
        num / den
    }
    /// Transmittance of the stack for an arbitrary incident tangential Jones
    /// vector `(E_x, E_y)`, using the normal Poynting flux weights of the
    /// incident and exit (isotropic, non-magnetic) media.
    pub fn transmittance(&self, input: [C; 2]) -> f64 {
        let t = self.t.apply(input);
        let (w_s, w_p) = self.incident_weights();
        let we_s = self.q_out.re;
        let we_p = (C::real(self.exit_i * self.exit_i) / self.q_out).re;
        let num = we_s * t[0].abs2() + we_p * t[1].abs2();
        let den = w_s * input[0].abs2() + w_p * input[1].abs2();
        num / den
    }
    /// Transmitted field for a unit Stokes-`S3` circular input. `handedness`
    /// `+1` is `(x + i y)/sqrt(2)`, `-1` is `(x - i y)/sqrt(2)`.
    pub fn circular_input(handedness: f64) -> [C; 2] {
        let inv = 1.0 / 2f64.sqrt();
        [C::real(inv), C::new(0.0, handedness * inv)]
    }
}

#[derive(Clone, Debug)]
pub struct AnisoStack {
    pub incident: f64,
    pub layers: Vec<AnisoLayer>,
    pub exit: f64,
}

fn iso_q(n: f64, s: f64) -> C {
    // Physical branch: evanescent waves decay along +z (Im q >= 0).
    let arg = C::real(n * n - s * s);
    arg.sqrt()
}

/// Mode vectors in an isotropic medium, in the order `(s+, s-, p+, p-)`, at
/// the reference plane `z = 0` (phases applied separately). The `s` mode has
/// unit tangential `E_y`, the `p` mode unit tangential `E_x`.
fn iso_modes(n: f64, _s: f64, q: C) -> [[C; 4]; 4] {
    let n2q = C::real(n * n) / q;
    [
        // s+ : (Ex, Ey, hx, hy)
        [C::real(0.0), C::real(1.0), -q, C::real(0.0)],
        // s-
        [C::real(0.0), C::real(1.0), q, C::real(0.0)],
        // p+
        [C::real(1.0), C::real(0.0), C::real(0.0), n2q],
        // p-
        [C::real(1.0), C::real(0.0), C::real(0.0), -n2q],
    ]
}

impl AnisoStack {
    pub fn new(incident: f64, layers: Vec<AnisoLayer>, exit: f64) -> Self {
        assert!(incident > 0.0 && exit > 0.0);
        Self {
            incident,
            layers,
            exit,
        }
    }

    /// Total 4x4 transfer matrix of the finite layers at conserved `s`. This
    /// convention maps the front fields to the back: `psi(d) = M psi(0)`, so
    /// the first layer acts first and the product is `P_N ... P_1`.
    pub fn transfer_matrix(&self, lambda: f64, s: f64) -> Mat4 {
        let mut m = Mat4::identity();
        for l in &self.layers {
            m = l.matrix(lambda, s).mul(&m);
        }
        m
    }

    pub fn solve(&self, lambda: f64, theta: f64) -> AnisoSolution {
        self.solve_s(lambda, self.incident * theta.sin())
    }

    /// Solve at a fixed conserved tangential wavevector `s = k_x / k0`.
    pub fn solve_s(&self, lambda: f64, s: f64) -> AnisoSolution {
        assert!(lambda > 0.0 && lambda.is_finite() && s.is_finite());
        let k0 = 2.0 * PI / lambda;
        let q0 = iso_q(self.incident, s);
        let q_out = iso_q(self.exit, s);
        let m = self.transfer_matrix(lambda, s);

        // Columns of V_in: s+, s-, p+, p-.
        let v_in = iso_modes(self.incident, s, q0);
        // Apply M to the incident-side modes.
        let mut u = [[C::real(0.0); 4]; 4];
        for c in 0..4 {
            u[c] = m.apply(&v_in[c]);
        }
        // Exit forward modes at z = d_total (common propagation phase).
        let d_total: f64 = self.layers.iter().map(|l| l.thickness).sum();
        let phase = (C::i() * (k0 * q_out.re * d_total)).exp();
        let v_out = iso_modes(self.exit, s, q_out);
        let v_out_f0 = v_out[0]; // s+
        let v_out_f1 = v_out[2]; // p+

        // Unknowns x = (r_s, r_p, t_s, t_p); A x = b.
        // A = [ U(:,s-), U(:,p-), -Vout_s, -Vout_p ].
        let mut a = [[C::real(0.0); 4]; 4];
        for i in 0..4 {
            a[i][0] = u[1][i];
            a[i][1] = u[3][i];
            a[i][2] = -(v_out_f0[i] * phase);
            a[i][3] = -(v_out_f1[i] * phase);
        }

        // For incident basis s and p.
        let solve_basis = |a_s: C, a_p: C| -> ([C; 2], [C; 2]) {
            let mut b = [C::real(0.0); 4];
            for i in 0..4 {
                b[i] = -(u[0][i] * a_s + u[2][i] * a_p);
            }
            let x = solve4(&a, &b).expect("Berreman system singular");
            ([x[0], x[1]], [x[2], x[3]])
        };
        let (r0, t0) = solve_basis(C::real(1.0), C::real(0.0)); // s input
        let (r1, t1) = solve_basis(C::real(0.0), C::real(1.0)); // p input

        let r = Jones {
            m: [[r0[0], r1[0]], [r0[1], r1[1]]],
        };
        let t = Jones {
            m: [[t0[0], t1[0]], [t0[1], t1[1]]],
        };

        // Flux weights: s -> Re(q); p -> Re(n^2 / q).
        let w_s0 = q0.re;
        let w_p0 = (C::real(self.incident * self.incident) / q0).re;
        let w_sout = q_out.re;
        let w_pout = (C::real(self.exit * self.exit) / q_out).re;

        // s input: reflected (r0[0], r0[1]), transmitted (t0[0], t0[1]).
        let refl_s = w_s0 * r0[0].abs2() + w_p0 * r0[1].abs2();
        let inc_s = w_s0;
        let trans_s = w_sout * t0[0].abs2() + w_pout * t0[1].abs2();
        // p input.
        let refl_p = w_s0 * r1[0].abs2() + w_p0 * r1[1].abs2();
        let inc_p = w_p0;
        let trans_p = w_sout * t1[0].abs2() + w_pout * t1[1].abs2();

        AnisoSolution {
            r,
            t,
            reflectance_s: refl_s / inc_s,
            transmittance_s: trans_s / inc_s,
            reflectance_p: refl_p / inc_p,
            transmittance_p: trans_p / inc_p,
            exit_forward_phase: phase,
            incident: self.incident,
            exit_i: self.exit,
            q0,
            q_out,
            lambda,
            s,
        }
    }
}

/// Stokes parameters `(S0, S1, S2, S3)` of a tangential Jones vector, with
/// `S3 = 2 Im(E_x E_y*)` (positive for the `+z`-propagating left-handed
/// convention of Born & Wolf). Normalise by `S0` for the Poincaré vector.
pub fn stokes(e: [C; 2]) -> [f64; 4] {
    let (ex, ey) = (e[0], e[1]);
    let s0 = ex.abs2() + ey.abs2();
    let s1 = ex.abs2() - ey.abs2();
    let exy = ex * ey.conj();
    [s0, s1, 2.0 * exy.re, 2.0 * exy.im]
}

// ---------------------------------------------------------------------------
// Liquid-crystal closures
// ---------------------------------------------------------------------------

/// The exact Jones matrix of a twisted nematic cell at normal incidence in
/// the input-director frame, from the coupled-mode (Mauguin) solution. With
/// twist `phi` (radians), retardation `gamma_r = 2 pi d (n_e - n_o)/lambda`
/// and `x = sqrt(phi^2 + (gamma_r/2)^2)`:
///
/// ```text
///   J = R(phi) [ cos x - i (gamma_r/(2x)) sin x ,  (phi/x) sin x
///                -(phi/x) sin x,  cos x + i (gamma_r/(2x)) sin x ]
/// ```
///
/// The Gooch–Tarry transmission between polarisers parallel/crossed to the
/// input director is `T_par = |J_11|^2`, `T_cross = |J_21|^2`. These are the
/// (approximate, slowly-varying) analytic benchmarks for the exact Berreman
/// solution in `tests/berreman.rs`.
pub fn twisted_nematic_jones(phi: f64, gamma_r: f64) -> Jones {
    let x = (phi * phi + 0.25 * gamma_r * gamma_r).sqrt();
    let (sx, cx) = (x.sin(), x.cos());
    let p = C::new(cx, -(gamma_r / (2.0 * x)) * sx);
    let q = C::real((phi / x) * sx);
    let local = Jones {
        m: [[p, q], [-q, p.conj()]],
    };
    let (cp, sp) = (phi.cos(), phi.sin());
    let rot = Jones {
        m: [
            [C::real(cp), C::real(-sp)],
            [C::real(sp), C::real(cp)],
        ],
    };
    rot * local
}

/// Gooch–Tarry transmission of a twisted nematic between polarisers crossed
/// to the input director: `T = 1 - (phi/x)^2 sin^2 x`,
/// `x = sqrt(phi^2 + (pi d Delta n / lambda)^2)` (Gooch & Tarry, J. Phys. D
/// **8**, 1575 (1975), doi:10.1088/0022-3727/8/13/007).
pub fn gooch_tarry_crossed(phi: f64, d_dn: f64, lambda: f64) -> f64 {
    let gamma_r = 2.0 * PI * d_dn / lambda;
    let x = (phi * phi + 0.25 * gamma_r * gamma_r).sqrt();
    if x == 0.0 {
        return 0.0;
    }
    let t = (phi / x) * x.sin();
    1.0 - t * t
}

/// Parallel-polariser Gooch–Tarry transmission, `(phi/x)^2 sin^2 x`.
pub fn gooch_tarry_parallel(phi: f64, d_dn: f64, lambda: f64) -> f64 {
    let gamma_r = 2.0 * PI * d_dn / lambda;
    let x = (phi * phi + 0.25 * gamma_r * gamma_r).sqrt();
    if x == 0.0 {
        return 0.0;
    }
    let t = (phi / x) * x.sin();
    t * t
}

/// Build a twisted nematic cell as a stack of `n_sub` thin uniaxial layers
/// whose optic axis rotates linearly from `0` to `twist` over `thickness`.
pub fn twisted_nematic(
    n_o: f64,
    n_e: f64,
    thickness: f64,
    twist: f64,
    n_sub: usize,
    incident: f64,
    exit: f64,
) -> AnisoStack {
    assert!(n_sub > 0);
    let dz = thickness / n_sub as f64;
    let mut layers = Vec::with_capacity(n_sub);
    for i in 0..n_sub {
        let theta = twist * (i as f64 + 0.5) / n_sub as f64;
        layers.push(AnisoLayer::new(
            uniaxial(C::real(n_o), C::real(n_e), [theta.cos(), theta.sin(), 0.0]),
            dz,
        ));
    }
    AnisoStack::new(incident, layers, exit)
}

/// Build a cholesteric (chiral nematic) layer: optic axis in the `x-y` plane
/// rotating by `2 pi` over the pitch `p`. `turns` full revolutions, `n_sub`
/// sublayers per revolution.
pub fn cholesteric(
    n_o: f64,
    n_e: f64,
    pitch: f64,
    turns: f64,
    n_sub: usize,
    incident: f64,
    exit: f64,
) -> AnisoStack {
    assert!(n_sub > 0 && turns > 0.0);
    let total = pitch * turns;
    let dz = total / (n_sub as f64 * turns);
    let n_layers = (n_sub as f64 * turns).round() as usize;
    let mut layers = Vec::with_capacity(n_layers);
    for i in 0..n_layers {
        let z = (i as f64 + 0.5) * dz;
        let theta = 2.0 * PI * z / pitch;
        layers.push(AnisoLayer::new(
            uniaxial(C::real(n_o), C::real(n_e), [theta.cos(), theta.sin(), 0.0]),
            dz,
        ));
    }
    AnisoStack::new(incident, layers, exit)
}

/// Leading-order optical-activity (rotatory power) of a short-pitch
/// cholesteric, `rho ~ pi (Delta n)^2 p^3 / (4 n_bar^2 lambda^2)` (de Vries,
/// Acta Cryst. **4**, 219 (1951)), used as the homogenised-medium diagnostic.
pub fn cholesteric_rotatory_power(n_o: f64, n_e: f64, pitch: f64, lambda: f64) -> f64 {
    let dn = n_e - n_o;
    let nbar = 0.5 * (n_e + n_o);
    PI * dn * dn * pitch.powi(3) / (4.0 * nbar * nbar * lambda * lambda)
}

/// Validity indicator of the cholesteric homogenisation (pitch-to-wavelength
/// ratio). Unfold to the exact Berreman stack when this exceeds the required
/// accuracy; the reflected power of a short-pitch cholesteric scales as
/// `O((p/lambda)^2)`.
pub fn effective_index_indicator(pitch: f64, lambda: f64) -> f64 {
    pitch / lambda
}

/// Simple convergence helper: the difference of two successive reflectance
/// values, so a caller can stop refining when it falls below a tolerance.
pub fn effective_index_convergence(r_fine: f64, r_coarse: f64) -> f64 {
    (r_fine - r_coarse).abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "{a:.12e} vs {b:.12e}");
    }

    #[test]
    fn isotropic_reduces_to_multilayer_normal_incidence() {
        use crate::multilayer::{Layer, Polarization, Stack};
        let lambda = 550e-9;
        // Glass / TiO2 / glass at normal incidence.
        let iso = Stack::new(
            1.0,
            vec![Layer::lossless(2.4, 80e-9), Layer::lossless(1.5, 120e-9)],
            1.52,
        );
        let ber = AnisoStack::new(
            1.0,
            vec![
                AnisoLayer::new(isotropic_index(C::real(2.4)), 80e-9),
                AnisoLayer::new(isotropic_index(C::real(1.5)), 120e-9),
            ],
            1.52,
        );
        let b = ber.solve(lambda, 0.0);
        let te = iso.solve(lambda, 0.0, Polarization::TE);
        let tm = iso.solve(lambda, 0.0, Polarization::TM);
        close(b.reflectance_s, te.reflectance, 1e-13);
        close(b.transmittance_s, te.transmittance, 1e-13);
        close(b.reflectance_p, tm.reflectance, 1e-13);
        close(b.transmittance_p, tm.transmittance, 1e-13);
    }

    #[test]
    fn isotropic_reduces_to_multilayer_oblique() {
        use crate::multilayer::{Layer, Polarization, Stack};
        let lambda = 633e-9;
        let layers = vec![
            Layer::lossless(2.4, 90e-9),
            Layer::lossless(1.5, 150e-9),
            Layer::lossless(2.4, 90e-9),
        ];
        let iso = Stack::new(1.5, layers.clone(), 1.0);
        let ber = AnisoStack::new(
            1.5,
            layers
                .iter()
                .map(|l| AnisoLayer::new(isotropic(l.index * l.index), l.thickness))
                .collect(),
            1.0,
        );
        for &theta in &[0.0f64, 0.2, 0.5, 0.8, 1.0] {
            let b = ber.solve(lambda, theta);
            let te = iso.solve(lambda, theta, Polarization::TE);
            let tm = iso.solve(lambda, theta, Polarization::TM);
            close(b.reflectance_s, te.reflectance, 1e-12);
            close(b.transmittance_s, te.transmittance, 1e-12);
            close(b.reflectance_p, tm.reflectance, 1e-12);
            close(b.transmittance_p, tm.transmittance, 1e-12);
        }
    }

    #[test]
    fn lossless_energy_conservation_anisotropic() {
        let lambda = 600e-9;
        let stack = cholesteric(1.5, 1.7, 400e-9, 2.0, 40, 1.0, 1.0);
        for &theta in &[0.0f64, 0.3, 0.6] {
            let s = stack.solve(lambda, theta);
            close(s.reflectance_s + s.transmittance_s, 1.0, 1e-10);
            close(s.reflectance_p + s.transmittance_p, 1.0, 1e-10);
        }
    }

    #[test]
    fn rotation_aligns_biaxial_tensor() {
        // biaxial(n_o, n_e, n_o) has its extraordinary axis along y; rotating
        // it onto +x must equal the uniaxial tensor with optic axis +x.
        let b = biaxial(C::real(1.5), C::real(1.7), C::real(1.5));
        let r = rotation([0.0, 0.0, 1.0], -PI / 2.0);
        let rb = rotate(&b, &r);
        let u = uniaxial(C::real(1.5), C::real(1.7), [1.0, 0.0, 0.0]);
        for i in 0..3 {
            for j in 0..3 {
                assert!((rb[i][j] - u[i][j]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn gyrotropic_faraday_rotation() {
        // eps_xy = -i g gives circular eigen-indices n_pm ~ n +/- g/(2n);
        // the transmitted linear polarisation rotates by
        // (k0/2)(n_+ - n_-) d ~ pi g d / (lambda n).
        let lambda = 633e-9;
        let n = 2.0;
        let g = 0.01;
        let d = 1e-6;
        let stack = AnisoStack::new(
            2.0,
            vec![AnisoLayer::new(
                gyrotropic(C::real(n * n), C::real(g), [0.0, 0.0, 1.0]),
                d,
            )],
            2.0,
        );
        let sol = stack.solve(lambda, 0.0);
        let out = sol.transmitted([C::real(1.0), C::real(0.0)]);
        let st = stokes(out);
        let angle = 0.5 * st[2].atan2(st[1]);
        let expected = PI * g * d / (lambda * n);
        assert!(
            (angle.abs() - expected).abs() < 1e-3 * expected,
            "Faraday rotation {:.6e} vs {:.6e}",
            angle.abs(),
            expected
        );
    }
}
