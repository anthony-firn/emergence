//! Plane-wave expansion (PWE) band structure of two-dimensional phononic
//! crystals.
//!
//! A square lattice of parallel circular cylinders (period `a`, radius `r`)
//! with isotropic elastic constituents. Two polarisations are implemented:
//!
//! * **Anti-plane / SH** (`u = w z_hat`), the scalar equation
//!   `div(mu grad w) + rho omega^2 w = 0`. Expanding `w` and the material
//!   functions in the reciprocal lattice gives the generalised eigenproblem
//!   `M_mu (K.K') w = omega^2 M_rho w`, which is **mathematically identical to
//!   TM photonics** `div((1/eps) grad H) + (omega^2/c^2) H = 0`: the same
//!   solver with `A = mu`, `B = rho` reduces to `A = 1/eps`, `B = 1/c^2`.
//! * **In-plane / mixed P-SV** (`u` in the `x-z` plane), the full anisotropic
//!   elastic equation `rho omega^2 u_i = K_j (C_ijkl)_{G-G'} K'_k u_l` with the
//!   isotropic stiffness.
//!
//! The Fourier coefficients of a material function `P(r) = P_m + (P_i-P_m)
//! chi(r)` are `P_0 = P_m + f (P_i-P_m)` and, for `G != 0`, `P_G =
//! (P_i-P_m) F(G)` with the disk structure factor
//! `F(G) = 2 f J_1(G r)/(G r)` and filling fraction `f = pi r^2/a^2`
//! (Kushwaha et al., Phys. Rev. Lett. **71**, 2022 (1993)).
//!
//! # Conditioning
//!
//! Internally the problem is made dimensionless: lengths are measured in `a`,
//! moduli in `mu_matrix` and densities in `rho_matrix`. The eigenvalues are
//! then `Omega^2 = omega^2 rho_m a^2 / mu_m` of order one, which keeps the
//! dense eigensolver well scaled.
//!
//! # Convergence
//!
//! PWE converges slowly for high elastic contrast (solid/fluid or stiff
//! inclusion in a soft matrix), because the displacement is discontinuous
//! across the interface and a Fourier series of a step converges only as
//! `1/G`. The number of plane waves is reported with every benchmark;
//! convergence is demonstrated by refining `N`. The low-frequency end needs
//! few waves and is used to check the folded effective medium.
use crate::materials::Solid;
use std::f64::consts::PI;

/// A square-lattice phononic crystal of circular cylinders.
#[derive(Clone, Copy, Debug)]
pub struct SquareLattice {
    /// Lattice constant [m].
    pub a: f64,
    /// Cylinder radius [m].
    pub radius: f64,
    pub inclusion: Solid,
    pub matrix: Solid,
}

impl SquareLattice {
    pub fn new(a: f64, radius: f64, inclusion: Solid, matrix: Solid) -> Self {
        assert!(a > 0.0 && radius > 0.0 && radius < 0.5 * a);
        Self {
            a,
            radius,
            inclusion,
            matrix,
        }
    }
    /// Filling fraction `f = pi r^2 / a^2`.
    pub fn filling_fraction(&self) -> f64 {
        PI * self.radius * self.radius / (self.a * self.a)
    }
    /// Disk structure factor `F(G)`: `f` at `G = 0`, `2 f J_1(Gr)/(Gr)`
    /// otherwise. `g` is the physical wavenumber.
    pub fn structure_factor(&self, g: f64) -> f64 {
        disk_structure_factor(self.radius, self.filling_fraction(), g)
    }
}

/// Bessel function of the first kind, order 1. Power series for `|x| < 8`
/// (accurate to round-off there) and the Abramowitz & Stegun 9.4.7
/// large-argument expansion for `|x| >= 8` (relative accuracy ~1e-8).
pub fn bessel_j1(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 8.0 {
        let h = 0.5 * x;
        let mut term = h;
        let mut sum = term;
        let mut m = 1.0;
        loop {
            term *= -(h * h) / (m * (m + 1.0));
            sum += term;
            if term.abs() <= 1e-17 * sum.abs().max(1e-300) {
                break;
            }
            m += 1.0;
            if m > 200.0 {
                break;
            }
        }
        sum
    } else {
        let z = 8.0 / ax;
        let y = z * z;
        let xx = ax - 2.356194491;
        let ans1 = 1.0
            + y * (0.183105e-2
                + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
        let ans2 = 0.04687499995
            + y * (-0.2002690873e-3
                + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
        let ans = (0.636619772 / ax).sqrt() * (xx.cos() * ans1 - z * xx.sin() * ans2);
        if x < 0.0 {
            -ans
        } else {
            ans
        }
    }
}

/// Disk structure factor `F(G)`: `f` at `G = 0`, `2 f J_1(Gr)/(Gr)`
/// otherwise.
fn disk_structure_factor(radius: f64, f: f64, g: f64) -> f64 {
    if g == 0.0 {
        f
    } else {
        let x = g * radius;
        2.0 * f * bessel_j1(x) / x
    }
}

/// Dimensionless plane-wave vectors `G a = 2 pi (n_x, n_y)`,
/// `n_x, n_y in [-m, m]`.
fn plane_waves(m: i32) -> Vec<(f64, f64)> {
    let mut v = Vec::with_capacity(((2 * m + 1) * (2 * m + 1)) as usize);
    for nx in -m..=m {
        for ny in -m..=m {
            v.push((2.0 * PI * nx as f64, 2.0 * PI * ny as f64));
        }
    }
    v
}

/// Dimensionless convolution matrix `c[p][q] = P_{G_p - G_q}`, with `P`
/// scaled by its matrix value.
fn conv_matrix(pw: &[(f64, f64)], radius_over_a: f64, f: f64, p_inc: f64, p_mat: f64) -> Vec<f64> {
    let n = pw.len();
    let mut c = vec![0.0; n * n];
    for p in 0..n {
        for q in 0..n {
            let dnx = ((pw[p].0 - pw[q].0) / (2.0 * PI)).round();
            let dny = ((pw[p].1 - pw[q].1) / (2.0 * PI)).round();
            let val = if dnx == 0.0 && dny == 0.0 {
                p_mat + (p_inc - p_mat) * f
            } else {
                let g = 2.0 * PI * (dnx * dnx + dny * dny).sqrt();
                (p_inc - p_mat) * disk_structure_factor(radius_over_a, f, g)
            };
            c[p * n + q] = val;
        }
    }
    c
}

/// Cholesky factor `L` (lower, row-major) of a real symmetric positive
/// definite matrix. Returns `None` if not positive definite.
pub fn cholesky(a: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut l = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = a[i * n + j];
            for k in 0..j {
                s -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if s <= 0.0 {
                    return None;
                }
                l[i * n + i] = s.sqrt();
            } else {
                l[i * n + j] = s / l[j * n + j];
            }
        }
    }
    Some(l)
}

/// Forward substitution `L y = b` for a lower-triangular `L` (row-major).
fn forward_solve(l: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut y = vec![0.0; n];
    for i in 0..n {
        let mut s = b[i];
        for k in 0..i {
            s -= l[i * n + k] * y[k];
        }
        y[i] = s / l[i * n + i];
    }
    y
}

/// Eigenvalues (ascending) of the real symmetric generalised problem
/// `A x = lambda B x`, `B` symmetric positive definite. `A` and `B` are
/// row-major `n x n`.
pub fn generalized_eigenvalues(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let l = cholesky(b, n).expect("mass matrix must be positive definite");
    // C = L^{-1} A L^{-T}. First Y = L^{-1} A (column by column).
    let mut y = vec![0.0; n * n];
    let mut col = vec![0.0; n];
    for j in 0..n {
        for i in 0..n {
            col[i] = a[i * n + j];
        }
        let sol = forward_solve(&l, &col, n);
        for i in 0..n {
            y[i * n + j] = sol[i];
        }
    }
    // Then C = (L^{-1} Y^T)^T.
    let mut z = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            col[i] = y[j * n + i]; // column j of Y^T is row j of Y
        }
        let sol = forward_solve(&l, &col, n);
        for i in 0..n {
            z[i * n + j] = sol[i];
        }
    }
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            c[i * n + j] = z[j * n + i];
        }
    }
    symmetric_eigenvalues(&mut c, n)
}

/// Eigenvalues (ascending) of a real symmetric matrix by Householder
/// tridiagonalisation followed by the implicit QL algorithm. `a` (row-major,
/// `n x n`) is overwritten. This is the fast path used by
/// [`generalized_eigenvalues`].
pub fn symmetric_eigenvalues(a: &mut [f64], n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];
    tridiagonalize(a, n, &mut d, &mut e);
    tql1(&mut d, &mut e, n);
    d.sort_by(|x, y| x.partial_cmp(y).unwrap());
    d
}

/// Householder reduction of a real symmetric matrix (row-major `a`) to
/// tridiagonal form: on return `d[i]` is the diagonal and `e[i]` the
/// subdiagonal (`e[0] = 0`).
fn tridiagonalize(a: &mut [f64], n: usize, d: &mut [f64], e: &mut [f64]) {
    for k in 0..n.saturating_sub(2) {
        let mut norm = 0.0;
        for i in k + 1..n {
            norm += a[i * n + k] * a[i * n + k];
        }
        let norm = norm.sqrt();
        if norm == 0.0 {
            continue;
        }
        let alpha = if a[(k + 1) * n + k] >= 0.0 {
            -norm
        } else {
            norm
        };
        let mut v = vec![0.0; n];
        v[k + 1] = a[(k + 1) * n + k] - alpha;
        for i in k + 2..n {
            v[i] = a[i * n + k];
        }
        let mut vtv = 0.0;
        for i in k + 1..n {
            vtv += v[i] * v[i];
        }
        if vtv == 0.0 {
            continue;
        }
        let v2 = 2.0 / vtv;
        // p = v2 A v.
        let mut p = vec![0.0; n];
        for i in 0..n {
            let mut s = 0.0;
            for j in k + 1..n {
                s += a[i * n + j] * v[j];
            }
            p[i] = v2 * s;
        }
        let mut pv = 0.0;
        for i in 0..n {
            pv += p[i] * v[i];
        }
        // A <- A - v w^T - w v^T.
        for i in 0..n {
            let w_i = p[i] - 0.5 * v2 * pv * v[i];
            for j in 0..n {
                a[i * n + j] -= v[i] * (p[j] - 0.5 * v2 * pv * v[j]) + w_i * v[j];
            }
        }
    }
    for i in 0..n {
        d[i] = a[i * n + i];
        e[i] = if i == 0 { 0.0 } else { a[i * n + (i - 1)] };
    }
}

/// Eigenvalues of a real symmetric tridiagonal matrix by the implicit QL
/// algorithm with Wilkinson shift (EISPACK/Numerical Recipes `tql1`).
/// `d` is the diagonal, `e` the subdiagonal with `e[0]` ignored.
fn tql1(d: &mut [f64], e: &mut [f64], n: usize) {
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;
    for l in 0..n {
        let mut iter = 0;
        loop {
            let mut m = l;
            while m < n - 1 {
                let dd = d[m].abs() + d[m + 1].abs();
                if e[m].abs() + dd == dd {
                    break;
                }
                m += 1;
            }
            if m == l {
                break;
            }
            iter += 1;
            assert!(iter <= 50, "tql1 did not converge");
            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let mut r = g.hypot(1.0);
            g = d[m] - d[l] + e[l] / (g + sign(r, g));
            let mut s = 1.0;
            let mut c = 1.0;
            let mut p = 0.0;
            let mut i = m;
            let mut broke = false;
            while i > l {
                i -= 1;
                let f = s * e[i];
                let b = c * e[i];
                r = f.hypot(g);
                e[i + 1] = r;
                if r == 0.0 {
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    broke = true;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;
            }
            if broke {
                continue;
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }
    }
}

/// `sign(magnitude, sign_of)` helper.
fn sign(a: f64, b: f64) -> f64 {
    if b >= 0.0 {
        a.abs()
    } else {
        -a.abs()
    }
}

/// Eigenvalues (ascending) of a real symmetric matrix by the cyclic Jacobi
/// method. `a` (row-major, `n x n`) is overwritten. Convergence is to
/// machine precision relative to the matrix norm.
pub fn jacobi_eigenvalues(a: &mut [f64], n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    let mut frob = 0.0f64;
    for v in a.iter() {
        frob += v * v;
    }
    let tol = 1e-30 * frob;
    for _sweep in 0..100 {
        let mut off = 0.0f64;
        for p in 0..n {
            for q in p + 1..n {
                off += a[p * n + q] * a[p * n + q];
            }
        }
        if off <= tol {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                let apq = a[p * n + q];
                if apq == 0.0 {
                    continue;
                }
                let app = a[p * n + p];
                let aqq = a[q * n + q];
                let theta = (aqq - app) / (2.0 * apq);
                let t = if theta >= 0.0 {
                    1.0 / (theta + (theta * theta + 1.0).sqrt())
                } else {
                    -1.0 / (-theta + (theta * theta + 1.0).sqrt())
                };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let akp = a[k * n + p];
                    let akq = a[k * n + q];
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let apk = a[p * n + k];
                    let aqk = a[q * n + k];
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                a[p * n + p] = app - t * apq;
                a[q * n + q] = aqq + t * apq;
                a[p * n + q] = 0.0;
                a[q * n + p] = 0.0;
            }
        }
    }
    let mut e: Vec<f64> = (0..n).map(|i| a[i * n + i]).collect();
    e.sort_by(|x, y| x.partial_cmp(y).unwrap());
    e
}

/// Dimensionless scalar eigenvalue squares `Omega^2` at dimensionless Bloch
/// wavevector `(kxa, kya)`.
fn scalar_eigenvalues_sq(
    radius_over_a: f64,
    f: f64,
    m: i32,
    kxa: f64,
    kya: f64,
    a_inc: f64,
    a_mat: f64,
    b_inc: f64,
    b_mat: f64,
) -> Vec<f64> {
    let pw = plane_waves(m);
    let n = pw.len();
    let ca = conv_matrix(&pw, radius_over_a, f, a_inc, a_mat);
    let cb = conv_matrix(&pw, radius_over_a, f, b_inc, b_mat);
    let mut amat = vec![0.0; n * n];
    for p in 0..n {
        let kxp = kxa + pw[p].0;
        let kyp = kya + pw[p].1;
        for q in 0..n {
            let kxq = kxa + pw[q].0;
            let kyq = kya + pw[q].1;
            amat[p * n + q] = ca[p * n + q] * (kxp * kxq + kyp * kyq);
        }
    }
    generalized_eigenvalues(&amat, &cb, n)
}

/// Anti-plane (SH) band frequencies [rad/s] at physical `(kx, ky)`, ascending.
/// `A` is the shear modulus `mu` and `B` the density `rho`.
pub fn sh_frequencies(lat: &SquareLattice, m: i32, kx: f64, ky: f64) -> Vec<f64> {
    let (mu_m, rho_m, a) = (lat.matrix.mu(), lat.matrix.rho, lat.a);
    let e = scalar_eigenvalues_sq(
        lat.radius / a,
        lat.filling_fraction(),
        m,
        kx * a,
        ky * a,
        lat.inclusion.mu() / mu_m,
        lat.matrix.mu() / mu_m,
        lat.inclusion.rho / rho_m,
        lat.matrix.rho / rho_m,
    );
    let scale = (mu_m / (rho_m * a * a)).sqrt();
    e.into_iter().map(|x| scale * x.max(0.0).sqrt()).collect()
}

/// TM-photonics reduction check: the scalar solver with `A = 1/eps` and a
/// constant `B = 1/c^2` returns the photonic frequencies `omega`. With
/// `mu <-> 1/eps` and `rho` constant the SH and TM band structures coincide.
pub fn tm_frequencies(
    eps_inc: f64,
    eps_mat: f64,
    a: f64,
    radius: f64,
    m: i32,
    kx: f64,
    ky: f64,
) -> Vec<f64> {
    let f = PI * radius * radius / (a * a);
    // Dimensionless: A = (1/eps)/(1/eps_mat) = eps_mat/eps, B = 1.
    let e = scalar_eigenvalues_sq(
        radius / a,
        f,
        m,
        kx * a,
        ky * a,
        eps_mat / eps_inc,
        1.0,
        1.0,
        1.0,
    );
    // Physical photonic frequency for c = 1: omega = Omega / a.
    let scale = 1.0 / a;
    e.into_iter().map(|x| scale * x.max(0.0).sqrt()).collect()
}

/// Dimensionless in-plane mixed P-SV eigenvalue squares `Omega^2` at
/// `(kxa, kya)`. The elastic operator is `M_il = lambda_{G-G'}(K.K') delta_il
/// + mu_{G-G'}(K_i K'_l + K_l K'_i)`; the mass is `rho_{G-G'} delta_il`.
fn inplane_eigenvalues_sq(lat: &SquareLattice, m: i32, kxa: f64, kya: f64) -> Vec<f64> {
    let pw = plane_waves(m);
    let n = pw.len();
    let a = lat.a;
    let f = lat.filling_fraction();
    let mu_m = lat.matrix.mu();
    let rho_m = lat.matrix.rho;
    let cl = conv_matrix(
        &pw,
        lat.radius / a,
        f,
        lat.inclusion.lambda() / mu_m,
        lat.matrix.lambda() / mu_m,
    );
    let cm = conv_matrix(
        &pw,
        lat.radius / a,
        f,
        lat.inclusion.mu() / mu_m,
        lat.matrix.mu() / mu_m,
    );
    let cr = conv_matrix(
        &pw,
        lat.radius / a,
        f,
        lat.inclusion.rho / rho_m,
        lat.matrix.rho / rho_m,
    );
    let dim = 2 * n;
    let mut mat = vec![0.0; dim * dim];
    let idx = |block: usize, p: usize| block * n + p;
    for p in 0..n {
        let kxp = kxa + pw[p].0;
        let kyp = kya + pw[p].1;
        for q in 0..n {
            let kxq = kxa + pw[q].0;
            let kyq = kya + pw[q].1;
            let kk = kxp * kxq + kyp * kyq;
            let l = cl[p * n + q];
            let mu = cm[p * n + q];
            // M_il = mu (K.K') delta_il + lambda K_i K'_l + mu K_l K'_i.
            mat[idx(0, p) * dim + idx(0, q)] = mu * kk + (l + mu) * kxp * kxq;
            mat[idx(1, p) * dim + idx(1, q)] = mu * kk + (l + mu) * kyp * kyq;
            mat[idx(0, p) * dim + idx(1, q)] = l * kxp * kyq + mu * kyp * kxq;
            mat[idx(1, p) * dim + idx(0, q)] = l * kyp * kxq + mu * kxp * kyq;
        }
    }
    let mut b = vec![0.0; dim * dim];
    for p in 0..n {
        for q in 0..n {
            b[idx(0, p) * dim + idx(0, q)] = cr[p * n + q];
            b[idx(1, p) * dim + idx(1, q)] = cr[p * n + q];
        }
    }
    generalized_eigenvalues(&mat, &b, dim)
}

/// In-plane band frequencies [rad/s] at physical `(kx, ky)`, ascending.
pub fn inplane_frequencies(lat: &SquareLattice, m: i32, kx: f64, ky: f64) -> Vec<f64> {
    let (mu_m, rho_m, a) = (lat.matrix.mu(), lat.matrix.rho, lat.a);
    let e = inplane_eigenvalues_sq(lat, m, kx * a, ky * a);
    let scale = (mu_m / (rho_m * a * a)).sqrt();
    e.into_iter().map(|x| scale * x.max(0.0).sqrt()).collect()
}

/// Wave polarisation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Polarization {
    /// Anti-plane (SH): `u` parallel to the cylinders.
    SH,
    /// In-plane mixed P-SV.
    InPlane,
}

/// Lowest-band phase speed along a direction, from the PWE eigenfrequency at
/// `(kx, ky)`. Used for the effective-medium check.
pub fn lowest_band_speed(lat: &SquareLattice, m: i32, kx: f64, ky: f64, kind: Polarization) -> f64 {
    let w = match kind {
        Polarization::SH => sh_frequencies(lat, m, kx, ky),
        Polarization::InPlane => inplane_frequencies(lat, m, kx, ky),
    };
    let k = (kx * kx + ky * ky).sqrt();
    w[0] / k
}

/// High-symmetry points of the square-lattice irreducible Brillouin zone.
/// Returns `(kx, ky, segment label)` along `Gamma-X-M-Gamma`.
pub fn band_path(a: f64, points_per_segment: usize) -> Vec<(f64, f64, &'static str)> {
    let x = PI / a;
    let m = PI / a;
    let mut path = Vec::new();
    for i in 0..=points_per_segment {
        let t = i as f64 / points_per_segment as f64;
        path.push((x * t, 0.0, "Gamma-X"));
    }
    for i in 1..=points_per_segment {
        let t = i as f64 / points_per_segment as f64;
        path.push((x, x * t, "X-M"));
    }
    for i in 1..=points_per_segment {
        let t = i as f64 / points_per_segment as f64;
        path.push((m * (1.0 - t), m * (1.0 - t), "M-Gamma"));
    }
    path
}

/// Normalised frequency `Omega = f a / c0` with `c0` a reference speed.
pub fn normalized_frequency(omega: f64, a: f64, c0: f64) -> f64 {
    omega * a / (2.0 * PI * c0)
}

/// Voigt (arithmetic) effective transverse speed, the reference used by
/// Kushwaha et al. to normalise the frequency axis.
pub fn voigt_transverse_speed(lat: &SquareLattice) -> f64 {
    let f = lat.filling_fraction();
    let mu = f * lat.inclusion.mu() + (1.0 - f) * lat.matrix.mu();
    let rho = f * lat.inclusion.rho + (1.0 - f) * lat.matrix.rho;
    (mu / rho).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "{a:e} vs {b:e}");
    }

    #[test]
    fn bessel_matches_reference() {
        // Reference values from tables (J1).
        close(bessel_j1(1.0), 0.4400505857449335, 1e-8);
        close(bessel_j1(5.0), -0.3275791375914653, 1e-8);
        close(bessel_j1(10.0), 0.0434727461688614, 1e-8);
        close(bessel_j1(20.0), 0.0668331241758502, 1e-8);
        close(bessel_j1(40.0), 0.1260383180375850, 1e-8);
        close(bessel_j1(60.0), 0.0465983837581662, 1e-8);
    }

    #[test]
    fn householder_matches_jacobi() {
        // An arbitrary symmetric 4x4 matrix.
        #[rustfmt::skip]
        let base = [
            4.0, 1.0, -2.0, 0.5,
            1.0, 3.0, 0.0, -1.0,
            -2.0, 0.0, 5.0, 0.25,
            0.5, -1.0, 0.25, 2.0,
        ];
        let mut a1 = base;
        let mut a2 = base;
        let fast = symmetric_eigenvalues(&mut a1, 4);
        let slow = jacobi_eigenvalues(&mut a2, 4);
        for i in 0..4 {
            close(fast[i], slow[i], 1e-12);
        }
    }

    #[test]
    fn jacobi_diagonal_matrix() {
        // Already diagonal: eigenvalues are the diagonal.
        let mut a = vec![3.0, 0.0, 0.0, 2.0];
        let e = jacobi_eigenvalues(&mut a, 2);
        close(e[0], 2.0, 1e-14);
        close(e[1], 3.0, 1e-14);
    }
}
