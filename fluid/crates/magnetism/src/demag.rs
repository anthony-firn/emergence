//! FFT demagnetising field on a finite-difference grid.
//!
//! The local-shape demagnetising field used by [`crate::micromag`] is a
//! one-cell closure: it ignores the non-local magnetostatic interaction.
//! The exact field is the convolution of the magnetisation with the
//! demagnetising tensor `N` of a rectangular cell,
//!
//! ```text
//! H_demag(r_i) = - sum_j N(r_i - r_j) M(r_j),
//! ```
//!
//! where `N` is the Newell-Williams-Dunlop tensor.  For cuboid cells it has
//! a closed form (Newell, Williams & Dunlop, *J. Geophys. Res.* **98**,
//! 9551 (1993)): for an observer at relative position `(x,y,z)` from the
//! centre of a cell of size `(dx,dy,dz)`, with the relative source
//! coordinates `xi = x - x'`, etc. evaluated at the eight corners
//! `x' = ±dx/2`, `y' = ±dy/2`, `z' = ±dz/2`,
//!
//! ```text
//! N_xx = -(1/4pi)[ CS_eta,zeta(atan(eta zeta/(xi r))) ]_{xi=+dx/2}^{xi=-dx/2} ...
//! ```
//!
//! and the off-diagonal components use `-ln(zeta + r)` type primitives
//! (see the private helpers).  The self term (`r = 0`) is the analytic
//! shape factor of the cell, e.g. `N_xx = N_yy = N_zz = 1/3` for a cube.
//! The implementation below was checked against direct numerical
//! quadrature of `∂'_i ∂'_j (1/r)` to machine precision (see
//! `tests/demag.rs::newell_tensor_matches_numerical_quadrature`).
//!
//! The non-local field is evaluated with a zero-padded FFT convolution:
//! the magnetisation is padded to twice the grid size so the circular
//! convolution of the padded arrays equals the linear convolution and the
//! open boundary is exact.  Grid sizes are rounded up to powers of two so a
//! plain radix-2 Cooley-Tukey transform suffices.
//!
//! The transform plan (twiddles and bit-reversal permutation) is built once
//! per grid and reused, the six complex working arrays are owned by a
//! caller-supplied [`FftWorkspace`] so the hot loop is allocation-free, and
//! the three transform axes and the tensor contraction run in parallel with
//! rayon.  The physics is unchanged; `field` returns bit-identical values
//! before and after the optimisation (the per-element arithmetic is the
//! same).
//!
//! The tensor is dimensionless and `H_demag` comes out in the same units as
//! `M` (A/m); the vacuum permeability enters only the energy.

use crate::geometry::{dot, V3};
use rayon::prelude::*;
use std::f64::consts::PI;

/// Complex pair `[re, im]`.
pub type C = [f64; 2];

const ZERO: C = [0.0, 0.0];

#[inline]
fn cadd(a: C, b: C) -> C {
    [a[0] + b[0], a[1] + b[1]]
}
#[inline]
fn csub(a: C, b: C) -> C {
    [a[0] - b[0], a[1] - b[1]]
}
#[inline]
fn cmul(a: C, b: C) -> C {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

/// A reusable radix-2 Cooley-Tukey FFT plan: the bit-reversal permutation
/// is precomputed so it is not rebuilt on every transform.
pub struct FftPlan {
    n: usize,
    rev: Vec<usize>,
}

impl FftPlan {
    pub fn new(n: usize) -> Self {
        assert!(n.is_power_of_two(), "FFT length {n} is not a power of two");
        let mut rev = vec![0usize; n];
        if n > 1 {
            let bits = n.trailing_zeros();
            for (i, r) in rev.iter_mut().enumerate() {
                *r = i.reverse_bits() >> (usize::BITS - bits);
            }
        }
        Self { n, rev }
    }

    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// In-place transform of `data[..n]`; `inverse` selects the sign of the
    /// exponent and applies the `1/n` normalisation.
    pub fn execute(&self, data: &mut [C], inverse: bool) {
        let n = self.n;
        if n <= 1 {
            return;
        }
        for (i, &j) in self.rev.iter().enumerate() {
            if i < j {
                data.swap(i, j);
            }
        }
        let sign = if inverse { 1.0 } else { -1.0 };
        let mut len = 2usize;
        while len <= n {
            let ang = sign * 2.0 * PI / len as f64;
            let wlen = [ang.cos(), ang.sin()];
            let half = len / 2;
            for start in (0..n).step_by(len) {
                let mut w = [1.0, 0.0];
                for k in 0..half {
                    let u = data[start + k];
                    let v = cmul(data[start + k + half], w);
                    data[start + k] = cadd(u, v);
                    data[start + k + half] = csub(u, v);
                    w = cmul(w, wlen);
                }
            }
            len <<= 1;
        }
        if inverse {
            let s = 1.0 / n as f64;
            for x in data.iter_mut() {
                x[0] *= s;
                x[1] *= s;
            }
        }
    }
}

/// In-place radix-2 Cooley-Tukey FFT.  `data.len()` must be a power of two.
/// Kept as a convenience wrapper; the hot path uses [`FftPlan::execute`].
pub fn fft(data: &mut [C], inverse: bool) {
    FftPlan::new(data.len()).execute(data, inverse);
}

/// Parallel 3-D FFT on an x-major flat array of length `mx*my*mz`.
///
/// `plans` holds the transform plans for the `x`, `y` and `z` axes and
/// `tmp` must have the same length as `data`; it is scratch for the x-axis
/// transpose.  All three axes are transformed with rayon.
pub fn fft3d_ws(
    data: &mut [C],
    mx: usize,
    my: usize,
    mz: usize,
    inverse: bool,
    plans: &[FftPlan; 3],
    tmp: &mut [C],
) {
    let m = mx * my * mz;
    debug_assert_eq!(data.len(), m);
    debug_assert_eq!(tmp.len(), m);

    // z (contiguous)
    if mz > 1 {
        data.par_chunks_mut(mz)
            .for_each(|row| plans[2].execute(row, inverse));
    }
    // y (stride mz)
    if my > 1 {
        data.par_chunks_mut(my * mz).for_each_init(
            || vec![ZERO; my],
            |buf, block| {
                for z in 0..mz {
                    for (y, b) in buf.iter_mut().enumerate() {
                        *b = block[y * mz + z];
                    }
                    plans[1].execute(buf, inverse);
                    for (y, b) in buf.iter().enumerate() {
                        block[y * mz + z] = *b;
                    }
                }
            },
        );
    }
    // x (stride my*mz): transpose to (y,z,x) so x is contiguous, transform,
    // transpose back.
    if mx > 1 {
        let nm = my * mz;
        tmp.par_chunks_mut(mx).enumerate().for_each(|(col, dst)| {
            let y = col / mz;
            let z = col % mz;
            for (x, d) in dst.iter_mut().enumerate() {
                *d = data[(x * my + y) * mz + z];
            }
        });
        tmp.par_chunks_mut(mx)
            .for_each(|row| plans[0].execute(row, inverse));
        data.par_chunks_mut(nm).enumerate().for_each(|(x, block)| {
            for y in 0..my {
                for z in 0..mz {
                    block[y * mz + z] = tmp[(y * mz + z) * mx + x];
                }
            }
        });
    }
}

/// Convenience 3-D FFT that builds its own plan and scratch.
pub fn fft3d(data: &mut [C], mx: usize, my: usize, mz: usize, inverse: bool) {
    let plans = [FftPlan::new(mx), FftPlan::new(my), FftPlan::new(mz)];
    let mut tmp = vec![ZERO; mx * my * mz];
    fft3d_ws(data, mx, my, mz, inverse, &plans, &mut tmp);
}

#[inline]
fn p_xx(xi: f64, eta: f64, zeta: f64) -> f64 {
    let r = (xi * xi + eta * eta + zeta * zeta).sqrt();
    (eta * zeta / (xi * r)).atan()
}
#[inline]
fn p_yy(xi: f64, eta: f64, zeta: f64) -> f64 {
    let r = (xi * xi + eta * eta + zeta * zeta).sqrt();
    (xi * zeta / (eta * r)).atan()
}
#[inline]
fn p_zz(xi: f64, eta: f64, zeta: f64) -> f64 {
    let r = (xi * xi + eta * eta + zeta * zeta).sqrt();
    (xi * eta / (zeta * r)).atan()
}
#[inline]
fn q_xy(xi: f64, eta: f64, zeta: f64) -> f64 {
    let r = (xi * xi + eta * eta + zeta * zeta).sqrt();
    -(zeta + r).ln()
}
#[inline]
fn q_xz(xi: f64, eta: f64, zeta: f64) -> f64 {
    let r = (xi * xi + eta * eta + zeta * zeta).sqrt();
    -(eta + r).ln()
}
#[inline]
fn q_yz(xi: f64, eta: f64, zeta: f64) -> f64 {
    let r = (xi * xi + eta * eta + zeta * zeta).sqrt();
    -(xi + r).ln()
}

/// Corner sum over the two remaining source coordinates:
/// `f(b2,c2) - f(b1,c2) - f(b2,c1) + f(b1,c1)`.
#[inline]
fn corner_sum<F: Fn(f64, f64, f64) -> f64>(f: F, fixed: f64, r1: [f64; 2], r2: [f64; 2]) -> f64 {
    let s = [-1.0, 1.0];
    let mut acc = 0.0;
    for j in 0..2 {
        for k in 0..2 {
            acc += s[j] * s[k] * f(fixed, r1[j], r2[k]);
        }
    }
    acc
}

/// Newell-Williams-Dunlop demagnetising tensor of a cuboid cell.
///
/// `offset` is the integer cell displacement from the observer to the
/// source; `cell` the cell dimensions.  Returns the symmetric tensor `N`
/// (dimensionless) such that the field at the observer is
/// `H_demag = - N . M_source`.
pub fn newell_tensor(offset: [i64; 3], cell: [f64; 3]) -> [[f64; 3]; 3] {
    let hx = 0.5 * cell[0];
    let hy = 0.5 * cell[1];
    let hz = 0.5 * cell[2];
    let x = offset[0] as f64 * cell[0];
    let y = offset[1] as f64 * cell[1];
    let z = offset[2] as f64 * cell[2];
    // Source corners at -h (index 0) and +h (index 1); relative coord = obs - source.
    let xi = [x + hx, x - hx];
    let eta = [y + hy, y - hy];
    let zeta = [z + hz, z - hz];

    let ixx = corner_sum(p_xx, xi[1], eta, zeta) - corner_sum(p_xx, xi[0], eta, zeta);
    let iyy = corner_sum(|f, a, b| p_yy(a, f, b), eta[1], xi, zeta)
        - corner_sum(|f, a, b| p_yy(a, f, b), eta[0], xi, zeta);
    let izz = corner_sum(|f, a, b| p_zz(a, b, f), zeta[1], xi, eta)
        - corner_sum(|f, a, b| p_zz(a, b, f), zeta[0], xi, eta);
    let ixy = corner_sum(q_xy, xi[1], eta, zeta) - corner_sum(q_xy, xi[0], eta, zeta);
    let ixz = corner_sum(q_xz, xi[1], eta, zeta) - corner_sum(q_xz, xi[0], eta, zeta);
    let iyz = corner_sum(|f, a, b| q_yz(a, f, b), eta[1], xi, zeta)
        - corner_sum(|f, a, b| q_yz(a, f, b), eta[0], xi, zeta);

    let c = -1.0 / (4.0 * PI);
    let (nxx, nyy, nzz, nxy, nxz, nyz) = (c * ixx, c * iyy, c * izz, c * ixy, c * ixz, c * iyz);
    [[nxx, nxy, nxz], [nxy, nyy, nyz], [nxz, nyz, nzz]]
}

#[inline]
fn next_pow2(n: usize) -> usize {
    let mut p = 1usize;
    while p < n {
        p <<= 1;
    }
    p
}

/// Reusable complex working arrays for [`DemagFft::field_ws`], sized to the
/// zero-padded grid.  Owning them here keeps the hot loop allocation-free.
pub struct FftWorkspace {
    fx: Vec<C>,
    fy: Vec<C>,
    fz: Vec<C>,
    hx: Vec<C>,
    hy: Vec<C>,
    hz: Vec<C>,
    tmp: Vec<C>,
    tmp2: Vec<C>,
    tmp3: Vec<C>,
}

impl FftWorkspace {
    pub fn new(padded_len: usize) -> Self {
        Self {
            fx: vec![ZERO; padded_len],
            fy: vec![ZERO; padded_len],
            fz: vec![ZERO; padded_len],
            hx: vec![ZERO; padded_len],
            hy: vec![ZERO; padded_len],
            hz: vec![ZERO; padded_len],
            tmp: vec![ZERO; padded_len],
            tmp2: vec![ZERO; padded_len],
            tmp3: vec![ZERO; padded_len],
        }
    }

    pub fn len(&self) -> usize {
        self.fx.len()
    }

    pub fn is_empty(&self) -> bool {
        self.fx.is_empty()
    }
}

/// Zero-padded FFT demagnetising-field operator for an `nx x ny x nz` grid.
///
/// The cell ordering is x-major: `index = (i*ny + j)*nz + k`.
pub struct DemagFft {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub mx: usize,
    pub my: usize,
    pub mz: usize,
    pub cell: [f64; 3],
    /// Fourier transforms of the six independent tensor components
    /// (xx, xy, xz, yy, yz, zz) on the padded grid.
    kernels: [Vec<C>; 6],
    /// Transform plans for the x, y and z axes.
    plans: [FftPlan; 3],
}

impl DemagFft {
    pub fn new(nx: usize, ny: usize, nz: usize, cell: [f64; 3]) -> Self {
        assert!(nx > 0 && ny > 0 && nz > 0);
        let mx = next_pow2(2 * nx - 1);
        let my = next_pow2(2 * ny - 1);
        let mz = next_pow2(2 * nz - 1);
        let m = mx * my * mz;

        let mut k = [
            vec![ZERO; m],
            vec![ZERO; m],
            vec![ZERO; m],
            vec![ZERO; m],
            vec![ZERO; m],
            vec![ZERO; m],
        ];
        // Place the kernel at wrapped offsets so that the circular
        // convolution on the doubled grid is the open-boundary linear
        // convolution.
        for ix in -(nx as i64 - 1)..=(nx as i64 - 1) {
            for iy in -(ny as i64 - 1)..=(ny as i64 - 1) {
                for iz in -(nz as i64 - 1)..=(nz as i64 - 1) {
                    let t = newell_tensor([ix, iy, iz], cell);
                    let ii = (ix.rem_euclid(mx as i64)) as usize;
                    let jj = (iy.rem_euclid(my as i64)) as usize;
                    let kk = (iz.rem_euclid(mz as i64)) as usize;
                    let idx = (ii * my + jj) * mz + kk;
                    k[0][idx] = [t[0][0], 0.0];
                    k[1][idx] = [t[0][1], 0.0];
                    k[2][idx] = [t[0][2], 0.0];
                    k[3][idx] = [t[1][1], 0.0];
                    k[4][idx] = [t[1][2], 0.0];
                    k[5][idx] = [t[2][2], 0.0];
                }
            }
        }
        let plans = [FftPlan::new(mx), FftPlan::new(my), FftPlan::new(mz)];
        for kk in k.iter_mut() {
            fft3d_ws(kk, mx, my, mz, false, &plans, &mut vec![ZERO; m]);
        }
        Self {
            nx,
            ny,
            nz,
            mx,
            my,
            mz,
            cell,
            kernels: k,
            plans,
        }
    }

    pub fn n_cells(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// Number of complex cells in the zero-padded transform grid.
    pub fn padded_len(&self) -> usize {
        self.mx * self.my * self.mz
    }

    /// Demagnetising field `H_demag = -N * M` (A/m) for unit-cell
    /// magnetisation `m` (the source `M = M_s m`).  `out` has the grid
    /// ordering of [`DemagFft`].
    ///
    /// This allocating entry point is kept for tests and one-off calls; the
    /// LLG hot loop uses [`DemagFft::field_ws`].
    pub fn field(&self, m: &[V3], out: &mut [V3]) {
        let mut ws = FftWorkspace::new(self.padded_len());
        self.field_ws(m, out, &mut ws);
    }

    /// Allocation-free demagnetising field: `ws` must be sized by
    /// [`DemagFft::padded_len`].
    pub fn field_ws(&self, m: &[V3], out: &mut [V3], ws: &mut FftWorkspace) {
        assert_eq!(m.len(), self.n_cells());
        assert_eq!(out.len(), self.n_cells());
        let mx = self.mx;
        let my = self.my;
        let mz = self.mz;
        let mp = mx * my * mz;
        debug_assert_eq!(ws.len(), mp);

        // Pack the magnetisation into the padded grid (zero elsewhere).
        ws.fx[..mp].fill(ZERO);
        ws.fy[..mp].fill(ZERO);
        ws.fz[..mp].fill(ZERO);
        for i in 0..self.nx {
            for j in 0..self.ny {
                for k in 0..self.nz {
                    let src = (i * self.ny + j) * self.nz + k;
                    let dst = (i * my + j) * mz + k;
                    ws.fx[dst] = [m[src][0], 0.0];
                    ws.fy[dst] = [m[src][1], 0.0];
                    ws.fz[dst] = [m[src][2], 0.0];
                }
            }
        }
        let plans = &self.plans;
        {
            let (fx, fy, fz) = (&mut ws.fx, &mut ws.fy, &mut ws.fz);
            let (t1, t2, t3) = (&mut ws.tmp, &mut ws.tmp2, &mut ws.tmp3);
            rayon::join(
                || fft3d_ws(fx, mx, my, mz, false, plans, t1),
                || {
                    rayon::join(
                        || fft3d_ws(fy, mx, my, mz, false, plans, t2),
                        || fft3d_ws(fz, mx, my, mz, false, plans, t3),
                    )
                },
            );
        }

        // H = -N M: accumulate N M into (hx,hy,hz) then negate.
        {
            let (fx, fy, fz) = (&ws.fx, &ws.fy, &ws.fz);
            let (hx, hy, hz) = (&mut ws.hx, &mut ws.hy, &mut ws.hz);
            let kernels = &self.kernels;
            hx.par_iter_mut()
                .zip(hy.par_iter_mut())
                .zip(hz.par_iter_mut())
                .enumerate()
                .for_each(|(idx, ((x, y), z))| {
                    let (a, b, c) = (fx[idx], fy[idx], fz[idx]);
                    let (nxx, nxy, nxz) = (kernels[0][idx], kernels[1][idx], kernels[2][idx]);
                    let (nyy, nyz, nzz) = (kernels[3][idx], kernels[4][idx], kernels[5][idx]);
                    *x = cadd(cadd(cmul(nxx, a), cmul(nxy, b)), cmul(nxz, c));
                    *y = cadd(cadd(cmul(nxy, a), cmul(nyy, b)), cmul(nyz, c));
                    *z = cadd(cadd(cmul(nxz, a), cmul(nyz, b)), cmul(nzz, c));
                });
        }
        let plans = &self.plans;
        {
            let (hx, hy, hz) = (&mut ws.hx, &mut ws.hy, &mut ws.hz);
            let (t1, t2, t3) = (&mut ws.tmp, &mut ws.tmp2, &mut ws.tmp3);
            rayon::join(
                || fft3d_ws(hx, mx, my, mz, true, plans, t1),
                || {
                    rayon::join(
                        || fft3d_ws(hy, mx, my, mz, true, plans, t2),
                        || fft3d_ws(hz, mx, my, mz, true, plans, t3),
                    )
                },
            );
        }
        for i in 0..self.nx {
            for j in 0..self.ny {
                for k in 0..self.nz {
                    let src = (i * my + j) * mz + k;
                    let dst = (i * self.ny + j) * self.nz + k;
                    out[dst] = [-ws.hx[src][0], -ws.hy[src][0], -ws.hz[src][0]];
                }
            }
        }
    }

    /// Magnetostatic energy (J) of `m` for saturation `ms`, using
    /// `E = -(mu0 V_cell ms^2 / 2) sum_i m_i . n_i` with `n = -N m`.
    pub fn energy(&self, m: &[V3], ms: f64, mu0: f64) -> f64 {
        let mut h = vec![[0.0; 3]; self.n_cells()];
        self.field(m, &mut h);
        let v = self.cell[0] * self.cell[1] * self.cell[2];
        let mut acc = 0.0;
        for i in 0..m.len() {
            acc += dot(m[i], h[i]);
        }
        // h = -N m is in units of M_s, so the physical field is ms*h.
        -0.5 * mu0 * v * ms * ms * acc
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn self_tensor_of_a_cube_is_one_third() {
        let t = newell_tensor([0, 0, 0], [1.0, 1.0, 1.0]);
        for a in 0..3 {
            assert!((t[a][a] - 1.0 / 3.0).abs() < 1e-14, "{a}");
            for b in 0..3 {
                if a != b {
                    assert!(t[a][b].abs() < 1e-14);
                }
            }
        }
    }

    #[test]
    fn tensor_trace_is_zero_off_self() {
        for off in [[1, 0, 0], [1, 1, 0], [2, 3, 4], [-1, 2, -3]] {
            let t = newell_tensor(off, [1.0, 1.0, 1.0]);
            let tr = t[0][0] + t[1][1] + t[2][2];
            assert!(tr.abs() < 1e-12, "off {off:?} trace {tr}");
        }
    }

    #[test]
    fn fft_uniform_cube_gives_minus_one_third() {
        let n = 5;
        let d = DemagFft::new(n, n, n, [1.0, 1.0, 1.0]);
        let m = vec![[1.0, 0.0, 0.0]; n * n * n];
        let mut h = vec![[0.0; 3]; n * n * n];
        d.field(&m, &mut h);
        let c = (n / 2 * n + n / 2) * n + n / 2;
        assert!((h[c][0] + 1.0 / 3.0).abs() < 1e-12, "Hx {}", h[c][0]);
        // Odd grid: the centre cell sits at the symmetry point, so the
        // transverse field is exactly zero for a cubic body.
        assert!(
            h[c][1].abs() < 1e-12 && h[c][2].abs() < 1e-12,
            "h={:?}",
            h[c]
        );
    }

    #[test]
    fn workspace_and_allocating_paths_agree() {
        // The reusable-workspace path must reproduce the convenience path.
        let d = DemagFft::new(7, 5, 3, [1.0, 1.2, 0.7]);
        let mut m = vec![[0.0; 3]; d.n_cells()];
        for (i, mi) in m.iter_mut().enumerate() {
            *mi = [
                (i as f64 * 0.37).sin(),
                (i as f64 * 0.11).cos(),
                (i as f64 * 0.53).sin(),
            ];
        }
        let mut h1 = vec![[0.0; 3]; d.n_cells()];
        d.field(&m, &mut h1);
        let mut h2 = vec![[0.0; 3]; d.n_cells()];
        let mut ws = FftWorkspace::new(d.padded_len());
        d.field_ws(&m, &mut h2, &mut ws);
        for (a, b) in h1.iter().zip(&h2) {
            assert!((a[0] - b[0]).abs() < 1e-12);
            assert!((a[1] - b[1]).abs() < 1e-12);
            assert!((a[2] - b[2]).abs() < 1e-12);
        }
    }

    #[test]
    fn tensor_values_match_quadrature_reference() {
        // Values cross-checked against direct numerical quadrature of
        // d'_i d'_j(1/r) in /tmp (see analysis/magnetism.md).
        let t = newell_tensor([1, 0, 0], [1.0, 1.0, 1.0]);
        assert!((t[0][0] + 0.1347824).abs() < 1e-6);
        let t = newell_tensor([1, 1, 0], [1.0, 1.0, 1.0]);
        assert!((t[0][0] + 0.0137462).abs() < 1e-6);
        assert!((t[0][1] + 0.0429224).abs() < 1e-6);
        assert!((t[2][2] - 0.0274923).abs() < 1e-6);
        let t = newell_tensor([0, 1, 0], [1.0, 1.0, 1.0]);
        assert!((t[0][0] - 0.0673912).abs() < 1e-6);
        assert!((t[1][1] + 0.1347824).abs() < 1e-6);
    }
}
