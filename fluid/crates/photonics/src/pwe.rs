//! Two-dimensional photonic crystals by the plane-wave expansion (PWE) method.
//!
//! The 2D crystal is a periodic array of circular rods (or holes) in a
//! background, translationally invariant along the rod axis. For a wave in the
//! plane, the two polarisations decouple and each gives a Hermitian eigenvalue
//! problem in the reciprocal lattice `G` (Ho, Chan & Soukoulis, PRL **65**,
//! 3152 (1990)):
//!
//! * **TM**, `E` parallel to the rods: the scalar field obeys
//!   `div(grad E) + (omega/c)^2 eps E = 0`, a *generalised* problem
//!   `diag(|k+G|^2) e = (omega/c)^2 eps e`.
//! * **TE**, `H` parallel to the rods: the master operator is
//!   `-div(eps^-1 grad H) = (omega/c)^2 H`.
//!
//! Both are discretised with the **inverse-`eps` (Ho–Chan–Soukoulis) rule**:
//! the truncated dielectric matrix `E_{GG'} = eps_{G-G'}` is inverted once per
//! `k` and enters symmetrically,
//!
//! ```text
//!   TM:  A_{GG'} = |k+G| |k+G'| [E^-1]_{GG'}
//!   TE:  A_{GG'} = (k+G).(k+G') [E^-1]_{GG'}
//! ```
//!
//! Both matrices are real symmetric and share the same `(omega/c)^2`
//! eigenvalues. The inverse rule is used for *both* polarisations because the
//! dielectric function is discontinuous at the rod surface: the direct Fourier
//! coefficients of `1/eps` decay only as `1/|G|`, so a direct plane-wave
//! truncation converges slowly (the ``staircase'' problem), whereas the matrix
//! inverse reproduces the exact answer with a handful of plane waves, while
//! the direct `1/eps` rule is still drifting at tens of plane waves. The
//! convergence is measured in `tests/optics.rs` against both published
//! benchmarks.
//!
//! Frequencies are reported as `omega a / 2 pi c`. The long-wavelength limit is
//! the effective-medium rung (`crate::effective_medium`); the validity
//! indicator is the lattice-constant-to-wavelength ratio, exactly as for the
//! 1D crystal.
//!
//! Benchmark (Joannopoulos et al., *Photonic Crystals: Molding the Flow of
//! Light*, 2nd ed., Princeton 2008, Fig. 5.2, 5.10 and Appendix C):
//!
//! * a square lattice of alumina rods, `eps = 8.9`, `r = 0.2 a`, has a TM gap
//!   `0.3224-0.4426`, a gap-midgap ratio of `31.4%` (the book quotes 31.4%),
//!   and no TE gap;
//! * a triangular lattice of air holes in `eps = 13`, `r = 0.48 a`, has a
//!   complete (TE and TM) gap `0.430-0.520`, a gap-midgap ratio of `19.0%`
//!   (the book quotes 18.6%).
//!
use std::f64::consts::PI;

/// Two-dimensional Bravais lattice of circular inclusions.
#[derive(Clone, Copy, Debug)]
pub enum Lattice {
    /// Square lattice with lattice constant `a` [m].
    Square { a: f64 },
    /// Triangular lattice with nearest-neighbour spacing `a` [m].
    Triangular { a: f64 },
}

impl Lattice {
    pub fn a(&self) -> f64 {
        match *self {
            Lattice::Square { a } | Lattice::Triangular { a } => a,
        }
    }
    /// Reciprocal basis `(b1, b2)` in units of `2 pi / a`.
    pub fn reciprocal_basis(&self) -> ([f64; 2], [f64; 2]) {
        match *self {
            Lattice::Square { .. } => ([1.0, 0.0], [0.0, 1.0]),
            Lattice::Triangular { .. } => (
                [1.0, 1.0 / 3.0_f64.sqrt()],
                [0.0, 2.0 / 3.0_f64.sqrt()],
            ),
        }
    }
    /// Area fraction of a circular inclusion of radius `r` (same units as `a`).
    pub fn area_fraction(&self, r: f64) -> f64 {
        let f = r / self.a();
        match *self {
            Lattice::Square { .. } => PI * f * f,
            Lattice::Triangular { .. } => 2.0 * PI * f * f / 3.0_f64.sqrt(),
        }
    }
    /// Cartesian `G = m b1 + n b2` in units of `2 pi / a`.
    pub fn g_vector(&self, m: i32, n: i32) -> (f64, f64) {
        let (b1, b2) = self.reciprocal_basis();
        (
            m as f64 * b1[0] + n as f64 * b2[0],
            m as f64 * b1[1] + n as f64 * b2[1],
        )
    }
    /// High-symmetry path for the band structure, in units of `2 pi / a`:
    /// `Gamma-X-M-Gamma` (square) or `Gamma-M-K-Gamma` (triangular).
    pub fn symmetry_path(&self, points_per_segment: usize) -> Vec<(f64, f64)> {
        let g = (0.0, 0.0);
        let corners: [(f64, f64); 3] = match *self {
            Lattice::Square { .. } => [g, (0.5, 0.0), (0.5, 0.5)],
            Lattice::Triangular { .. } => {
                let (b1, b2) = self.reciprocal_basis();
                let m = (0.5 * b1[0], 0.5 * b1[1]);
                let k = ((b1[0] + b2[0]) / 3.0, (b1[1] + b2[1]) / 3.0);
                [g, m, k]
            }
        };
        let mut path = Vec::with_capacity(3 * points_per_segment + 1);
        for s in 0..3 {
            let a = corners[s];
            let b = corners[(s + 1) % 3];
            for i in 0..points_per_segment {
                let t = i as f64 / points_per_segment as f64;
                path.push((a.0 + (b.0 - a.0) * t, a.1 + (b.1 - a.1) * t));
            }
        }
        path.push(g);
        path
    }
}

/// Polarisation of a 2D crystal mode.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Polarization2D {
    /// `E` parallel to the rods (out of plane); `A = |k+G| |k+G'| [E^-1]`.
    TM,
    /// `H` parallel to the rods; `A = (k+G).(k+G') [E^-1]`.
    TE,
}

/// A 2D photonic crystal and its plane-wave truncation.
#[derive(Clone, Debug)]
pub struct PweCrystal {
    pub lattice: Lattice,
    /// Permittivity of the circular inclusions (rods or holes).
    pub eps_inclusion: f64,
    /// Permittivity of the background.
    pub eps_background: f64,
    /// Inclusion radius (same length units as `a`).
    pub radius: f64,
    /// Plane-wave cutoff `N`: the grid runs over `m, n in [-N, N]`.
    pub n_grid: i32,
}

/// Bessel function `J1(x)` to about `1e-8` (better for most `x`): the power
/// series for `|x| < 8` and the asymptotic expansion
/// `J1 ~ sqrt(2/(pi x)) (cos(chi) P - sin(chi) Q)`, `chi = x - 3 pi/4`, for
/// larger arguments. The series is summed until the term is below `1e-18`
/// relative, so it is accurate to round-off in its range.
pub fn bessel_j1(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 1e-12 {
        return 0.0;
    }
    if ax < 8.0 {
        // J1(x) = sum_k (-1)^k (x/2)^(2k+1) / (k! (k+1)!).
        let half = x / 2.0;
        let mut term = half;
        let mut sum = term;
        let mut k = 0.0;
        while term.abs() > 1e-18 * sum.abs().max(1e-300) && k < 300.0 {
            k += 1.0;
            term *= -(half * half) / (k * (k + 1.0));
            sum += term;
        }
        sum
    } else {
        // a_k(1) = prod_{j=1}^k (4 - (2j-1)^2) / (8 j), a_0 = 1.
        let mut a = [0.0f64; 10];
        a[0] = 1.0;
        for k in 1..10 {
            let kf = k as f64;
            a[k] = a[k - 1] * (4.0 - (2.0 * kf - 1.0).powi(2)) / (8.0 * kf);
        }
        let (mut p, mut q) = (0.0, 0.0);
        let mut k = 0usize;
        while k < 10 {
            p += if (k / 2) % 2 == 0 { a[k] } else { -a[k] } / ax.powi(k as i32);
            k += 2;
        }
        k = 1;
        while k < 10 {
            q += if ((k - 1) / 2) % 2 == 0 { a[k] } else { -a[k] } / ax.powi(k as i32);
            k += 2;
        }
        let chi = ax - 0.75 * PI;
        let j = (2.0 / (PI * ax)).sqrt() * (chi.cos() * p - chi.sin() * q);
        if x < 0.0 {
            -j
        } else {
            j
        }
    }
}

/// Fourier coefficient of the piecewise-constant circular inclusion.
///
/// `dg` is the magnitude of `G - G'` in units of `2 pi / a`, `r` the radius in
/// the same units as `a`, and `fill` the area fraction. For `dg = 0` the
/// coefficient is the area average; otherwise it is
/// `(val_inc - val_bg) * fill * 2 J1(2 pi dg r) / (2 pi dg r)`.
fn inclusion_coefficient(dg: f64, val_inc: f64, val_bg: f64, fill: f64, r: f64) -> f64 {
    if dg < 1e-12 {
        val_bg + (val_inc - val_bg) * fill
    } else {
        let x = 2.0 * PI * dg * r;
        (val_inc - val_bg) * fill * 2.0 * bessel_j1(x) / x
    }
}

impl PweCrystal {
    pub fn new(
        lattice: Lattice,
        eps_inclusion: f64,
        eps_background: f64,
        radius: f64,
        n_grid: i32,
    ) -> Self {
        assert!(eps_inclusion > 0.0 && eps_background > 0.0);
        assert!(radius > 0.0);
        assert!(n_grid >= 1);
        Self {
            lattice,
            eps_inclusion,
            eps_background,
            radius,
            n_grid,
        }
    }

    pub fn area_fraction(&self) -> f64 {
        self.lattice.area_fraction(self.radius)
    }

    /// Plane-wave `G` vectors in units of `2 pi / a`, as Cartesian pairs.
    pub fn plane_waves(&self) -> Vec<(f64, f64)> {
        let n = self.n_grid;
        let mut out = Vec::with_capacity(((2 * n + 1) * (2 * n + 1)) as usize);
        for m in -n..=n {
            for l in -n..=n {
                out.push(self.lattice.g_vector(m, l));
            }
        }
        out
    }

    /// Number of plane waves retained.
    pub fn number_of_plane_waves(&self) -> usize {
        let n = (2 * self.n_grid + 1) as usize;
        n * n
    }

    fn k_vectors(&self, k: (f64, f64)) -> (Vec<f64>, Vec<f64>) {
        let pw = self.plane_waves();
        let mut kx = Vec::with_capacity(pw.len());
        let mut ky = Vec::with_capacity(pw.len());
        for (gx, gy) in pw {
            kx.push(k.0 + gx);
            ky.push(k.1 + gy);
        }
        (kx, ky)
    }

    /// The real symmetric plane-wave matrix whose eigenvalues are
    /// `(omega a / 2 pi c)^2`.
    pub fn matrix(&self, k: (f64, f64), pol: Polarization2D) -> Vec<f64> {
        let pw = self.plane_waves();
        let n = pw.len();
        let (kx, ky) = self.k_vectors(k);
        let fill = self.area_fraction();
        let r = self.radius / self.lattice.a();
        // The dielectric matrix is shared by both polarisations; invert it
        // once with the Ho rule. The diagonal term `eps_{G-G'}` is the area
        // average, the off-diagonal the 2 J1(x)/x Bessel coefficient.
        let mut b = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                let dgx = pw[i].0 - pw[j].0;
                let dgy = pw[i].1 - pw[j].1;
                let dg = (dgx * dgx + dgy * dgy).sqrt();
                b[i * n + j] = inclusion_coefficient(
                    dg,
                    self.eps_inclusion,
                    self.eps_background,
                    fill,
                    r,
                );
            }
        }
        let binv = cholesky_inverse(&mut b, n);
        let mut a = vec![0.0; n * n];
        match pol {
            Polarization2D::TE => {
                // A_{ij} = (k+G_i).(k+G_j) [eps^-1]_{ij}.
                for i in 0..n {
                    for j in 0..n {
                        a[i * n + j] =
                            binv[i * n + j] * (kx[i] * kx[j] + ky[i] * ky[j]);
                    }
                }
            }
            Polarization2D::TM => {
                // A_{ij} = |k+G_i| |k+G_j| [eps^-1]_{ij}.
                for i in 0..n {
                    let ki = (kx[i] * kx[i] + ky[i] * ky[i]).sqrt();
                    for j in 0..n {
                        let kj = (kx[j] * kx[j] + ky[j] * ky[j]).sqrt();
                        a[i * n + j] = ki * kj * binv[i * n + j];
                    }
                }
            }
        }
        a
    }

    /// Lowest `n_bands` normalised frequencies `omega a / 2 pi c` at `k`
    /// (Cartesian, in units of `2 pi / a`).
    pub fn bands(&self, k: (f64, f64), pol: Polarization2D, n_bands: usize) -> Vec<f64> {
        let a = self.matrix(k, pol);
        let n = self.number_of_plane_waves();
        lanczos_lowest(&a, n, n_bands)
            .into_iter()
            .map(|lam| lam.max(0.0).sqrt())
            .collect()
    }

    /// Long-wavelength slope `d(omega a/2 pi c)/d(k a/2 pi)` at `Gamma`,
    /// evaluated by finite difference at `k = dk (1, 0)`.
    pub fn lowest_band_slope(&self, pol: Polarization2D, dk: f64) -> f64 {
        let f = self.bands((dk, 0.0), pol, 1)[0];
        f / dk
    }

    /// Band edges of the lowest `n_bands` bands along the high-symmetry path.
    /// Returns, for each band index, `(min, max)`.
    pub fn band_extrema(
        &self,
        pol: Polarization2D,
        n_bands: usize,
        points_per_segment: usize,
    ) -> Vec<(f64, f64)> {
        let path = self.lattice.symmetry_path(points_per_segment);
        let mut ext = vec![(f64::INFINITY, f64::NEG_INFINITY); n_bands];
        for k in path {
            let f = self.bands(k, pol, n_bands);
            for (i, &fi) in f.iter().enumerate() {
                if fi < ext[i].0 {
                    ext[i].0 = fi;
                }
                if fi > ext[i].1 {
                    ext[i].1 = fi;
                }
            }
        }
        ext
    }

    /// All polarisation gaps among the lowest `n_bands` bands along the
    /// high-symmetry path, as `(lower_band, gap_lo, gap_hi)`. A gap exists
    /// when the maximum of band `lower` is below the minimum of band
    /// `lower + 1`.
    pub fn band_gaps(
        &self,
        pol: Polarization2D,
        n_bands: usize,
        points_per_segment: usize,
    ) -> Vec<(usize, f64, f64)> {
        let ext = self.band_extrema(pol, n_bands, points_per_segment);
        let mut out = Vec::new();
        for lower in 0..n_bands.saturating_sub(1) {
            let lo = ext[lower].1;
            let hi = ext[lower + 1].0;
            if lo < hi {
                out.push((lower, lo, hi));
            }
        }
        out
    }

    /// A polarisation gap: `(gap_lo, gap_hi)` from the maximum of band `lower`
    /// and the minimum of band `lower + 1`; `None` if they overlap. Midgap
    /// width is `(hi - lo) / ((hi + lo) / 2)`.
    pub fn band_gap(
        &self,
        pol: Polarization2D,
        lower: usize,
        points_per_segment: usize,
    ) -> Option<(f64, f64)> {
        let ext = self.band_extrema(pol, lower + 2, points_per_segment);
        let lo = ext[lower].1;
        let hi = ext[lower + 1].0;
        if lo < hi {
            Some((lo, hi))
        } else {
            None
        }
    }

    /// Widest complete (TE and TM) gap among adjacent bands: the intersection
    /// of every TE polarisation gap with every TM polarisation gap, keeping the
    /// one with the largest gap-midgap ratio. For an air-hole crystal the two
    /// polarisations can have their gaps between different band indices, so all
    /// pairings are searched rather than assuming a fixed one.
    pub fn complete_gap(
        &self,
        n_bands: usize,
        points_per_segment: usize,
    ) -> Option<(f64, f64)> {
        let te = self.band_gaps(Polarization2D::TE, n_bands, points_per_segment);
        let tm = self.band_gaps(Polarization2D::TM, n_bands, points_per_segment);
        let mut best: Option<(f64, f64)> = None;
        let mut best_ratio = -1.0;
        for &(_, tlo, thi) in &te {
            for &(_, mlo, mhi) in &tm {
                let lo = tlo.max(mlo);
                let hi = thi.min(mhi);
                if lo < hi {
                    let ratio = (hi - lo) / (0.5 * (hi + lo));
                    if ratio > best_ratio {
                        best_ratio = ratio;
                        best = Some((lo, hi));
                    }
                }
            }
        }
        best
    }
}

/// Gap-midgap ratio `(hi - lo) / ((hi + lo)/2)` of a gap.
pub fn gap_midgap_ratio(gap: (f64, f64)) -> f64 {
    (gap.1 - gap.0) / (0.5 * (gap.0 + gap.1))
}

/// In-place Cholesky factorisation then inverse of a symmetric positive
/// definite matrix given row-major.
pub fn cholesky_inverse(a: &mut [f64], n: usize) -> Vec<f64> {
    // Lower triangular L with A = L L^T.
    for i in 0..n {
        for j in 0..=i {
            let mut s = a[i * n + j];
            for k in 0..j {
                s -= a[i * n + k] * a[j * n + k];
            }
            if i == j {
                assert!(s > 0.0, "dielectric matrix is not positive definite");
                a[i * n + i] = s.sqrt();
            } else {
                a[i * n + j] = s / a[j * n + j];
            }
        }
    }
    // Inverse of L by forward substitution, stored in `linv` (lower).
    let mut linv = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            if i < j {
                continue;
            }
            let mut s = if i == j { 1.0 } else { 0.0 };
            for k in j..i {
                s -= a[i * n + k] * linv[k * n + j];
            }
            linv[i * n + j] = s / a[i * n + i];
        }
    }
    // A^-1 = L^-T L^-1.
    let mut inv = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for k in i.max(j)..n {
                s += linv[k * n + i] * linv[k * n + j];
            }
            inv[i * n + j] = s;
        }
    }
    inv
}

fn matvec(a: &[f64], n: usize, x: &[f64], y: &mut [f64]) {
    for i in 0..n {
        let row = &a[i * n..i * n + n];
        let mut s = 0.0;
        for j in 0..n {
            s += row[j] * x[j];
        }
        y[i] = s;
    }
}

fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).map(|(a, b)| a * b).sum()
}

/// Eigenvalues of a small real symmetric matrix, ascending (cyclic Jacobi).
pub fn jacobi_eigenvalues_sym(a: &mut [f64], n: usize) -> Vec<f64> {
    for _sweep in 0..100 {
        let mut off = 0.0;
        for p in 0..n {
            for q in p + 1..n {
                off += a[p * n + q] * a[p * n + q];
            }
        }
        if off.sqrt() < 1e-14 {
            break;
        }
        for p in 0..n.saturating_sub(1) {
            for q in p + 1..n {
                let apq = a[p * n + q];
                if apq == 0.0 {
                    continue;
                }
                let theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
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
            }
        }
    }
    (0..n).map(|i| a[i * n + i]).collect()
}

/// Lowest `n_bands` eigenvalues of a dense real symmetric matrix by Lanczos
/// iteration with full reorthogonalisation. The tridiagonal Ritz problem is
/// diagonalised with [`jacobi_eigenvalues_sym`].
///
/// The start vector is a deterministic linear-congruential pseudo-random
/// vector. A structured start (e.g. a low-frequency sinusoid) can be
/// orthogonal to the ground state at particular plane-wave counts, which makes
/// Lanczos return the next eigenvalue instead of the lowest; a pseudo-random
/// start has non-negligible overlap with every eigenvector.
pub fn lanczos_lowest(a: &[f64], n: usize, n_bands: usize) -> Vec<f64> {
    let n_bands = n_bands.max(1).min(n);
    let m = (12 * n_bands + 120).min(n).max(n_bands + 1);
    let mut state: u64 =
        0x9E37_79B9_7F4A_7C15 ^ (n as u64) ^ ((n_bands as u64) << 32);
    let mut q: Vec<f64> = (0..n)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 33) as f64) / ((1u64 << 31) as f64) - 1.0
        })
        .collect();
    normalize(&mut q);
    let mut qs: Vec<Vec<f64>> = Vec::with_capacity(m);
    let mut alpha = vec![0.0; m];
    let mut beta = vec![0.0; m];
    qs.push(q.clone());
    let mut v = vec![0.0; n];
    let mut steps = m;
    for j in 0..m {
        matvec(a, n, &qs[j], &mut v);
        alpha[j] = dot(&qs[j], &v);
        for i in 0..n {
            v[i] -= alpha[j] * qs[j][i];
            if j > 0 {
                v[i] -= beta[j - 1] * qs[j - 1][i];
            }
        }
        // Full reorthogonalisation (twice for numerical stability).
        for _pass in 0..2 {
            for qi in qs.iter() {
                let c = dot(qi, &v);
                for i in 0..n {
                    v[i] -= c * qi[i];
                }
            }
        }
        let b = dot(&v, &v).sqrt();
        if j + 1 >= m {
            break;
        }
        beta[j] = b;
        if b < 1e-14 {
            steps = j + 1;
            break;
        }
        let next: Vec<f64> = v.iter().map(|x| x / b).collect();
        qs.push(next);
    }
    // Tridiagonal Ritz matrix.
    let d = steps.min(qs.len());
    let mut t = vec![0.0; d * d];
    for i in 0..d {
        t[i * d + i] = alpha[i];
        if i + 1 < d {
            t[i * d + i + 1] = beta[i];
            t[(i + 1) * d + i] = beta[i];
        }
    }
    let mut e = jacobi_eigenvalues_sym(&mut t, d);
    e.sort_by(|x, y| x.total_cmp(y));
    e.truncate(n_bands);
    e
}

fn normalize(v: &mut [f64]) {
    let norm = dot(v, v).sqrt();
    assert!(norm > 0.0, "zero Lanczos start vector");
    for x in v.iter_mut() {
        *x /= norm;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bessel_j1_matches_series_and_recursion() {
        // Values from the standard tables: J1(1) = 0.4400505857,
        // J1(2) = 0.5767248078, J1(5) = -0.3275791376,
        // J1(10) = 0.0434727462.
        let table = [(1.0, 0.4400505857), (2.0, 0.5767248078), (5.0, -0.3275791376), (10.0, 0.0434727462)];
        for (x, j) in table {
            assert!((bessel_j1(x) - j).abs() < 5e-8, "J1({x}) = {}", bessel_j1(x));
        }
    }

    #[test]
    fn empty_lattice_is_free_photon_dispersion() {
        // eps = 1 everywhere: f = |k + G| in units of 2 pi / a.
        let c = PweCrystal::new(Lattice::Square { a: 1.0 }, 1.0, 1.0, 0.2, 4);
        let f = c.bands((0.25, 0.0), Polarization2D::TM, 3);
        // The lowest folded free-photon frequency at k=(0.25,0) is 0.25.
        assert!((f[0] - 0.25).abs() < 1e-10, "{f:?}");
        let g = c.bands((0.25, 0.0), Polarization2D::TE, 3);
        assert!((g[0] - 0.25).abs() < 1e-10, "{g:?}");
    }

    #[test]
    fn cholesky_inverse_is_inverse() {
        let mut a = vec![4.0, 1.0, 1.0, 1.0, 3.0, 0.5, 1.0, 0.5, 2.0];
        let orig = a.clone();
        let inv = cholesky_inverse(&mut a, 3);
        for i in 0..3 {
            for j in 0..3 {
                let mut s = 0.0;
                for k in 0..3 {
                    s += orig[i * 3 + k] * inv[k * 3 + j];
                }
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((s - want).abs() < 1e-12, "({i},{j}) = {s}");
            }
        }
    }

    #[test]
    fn lanczos_finds_lowest_of_small_matrix() {
        let a = [2.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 2.0];
        let all = jacobi_eigenvalues_sym(&mut a.clone(), 3);
        let low = lanczos_lowest(&a, 3, 2);
        let mut all = all;
        all.sort_by(|x, y| x.total_cmp(y));
        assert!((low[0] - all[0]).abs() < 1e-10, "{low:?} {all:?}");
        assert!((low[1] - all[1]).abs() < 1e-10, "{low:?} {all:?}");
    }
}
