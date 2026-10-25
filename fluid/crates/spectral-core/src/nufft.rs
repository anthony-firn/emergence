//! Non-uniform FFT between the solver's retained Fourier modes and arbitrary
//! points, by Gaussian gridding (Greengard & Lee, SIAM Rev. 46, 2004).
//!
//! * interpolation (type 2): u(X_m) = Σ_k û_k e^{ik·X_m}
//! * spreading (type 1):     f̂_k  = (2π)⁻³ Σ_m F_m e^{−ik·X_m}
//!
//! i.e. f = Σ_m F_m δ(x − X_m) projected onto the retained modes. The two
//! are adjoint up to the factor (2π)³, so J·H⁻¹P·S is symmetric positive
//! semi-definite (used for the surface-constraint solve in `imex.rs`).
//!
//! Method (1D notation, products in 3D): with g_τ(x) = e^{−x²/(4τ)} and
//! φ(k) = 2π ĝ_τ(k) = 2√(πτ) e^{−k²τ},
//!   interpolate: W = û/Φ on an oversampled grid, w = F⁻¹W, u(X) = h³ Σ_j w_j g(X − x_j)
//!   spread:      F_τ(x_j) = Σ_m F_m g(x_j − X_m), f̂ = F[F_τ]/Φ
//! where Φ = φ(kx)φ(ky)φ(kz) and the Gaussian is truncated to ±`half_width`
//! oversampled points. τ balances the truncation error exp(−π²w²/(M²τ))
//! against the aliasing error exp(−M²τ(1 − 1/R)) for oversampling
//! R = M/(2k_max+1): M²τ = πw/√(1 − 1/R), error ≈ exp(−πw√(1 − 1/R)) times
//! the deconvolution gain exp(k_max²τ).

use num_complex::Complex64;
use rayon::prelude::*;

use crate::fft::Fft3;
use crate::grid::{wavenumber, Grid};
use crate::solver::Spec;

pub struct Nufft {
    /// oversampled grid
    pub m: usize,
    pub tau: f64,
    pub half_width: usize,
    fft: Fft3,
    /// oversampled spectral index ↔ N-grid spectral index for retained modes
    pairs: Vec<(usize, usize)>,
    /// 1/Φ(k) for each pair
    inv_phi: Vec<f64>,
    work_r: Vec<f64>,
    work_c: Spec,
}

impl Nufft {
    /// For an N grid (any truncation): oversampling ≥ 2, Gaussian half width
    /// `half_width` oversampled points (8 gives ~1e-7, 12 ~1e-10 relative accuracy).
    pub fn new(gn: &Grid, half_width: usize) -> Self {
        let kmax = gn.kmax();
        let modes = 2 * kmax + 1;
        let mut m = 2 * modes;
        m += m % 2;
        let r = m as f64 / modes as f64;
        let tau = std::f64::consts::PI * half_width as f64 / ((1.0 - 1.0 / r).sqrt() * (m * m) as f64);
        let gm = Grid::new(m.max(8));
        let phi = |k: f64| 2.0 * (std::f64::consts::PI * tau).sqrt() * (-k * k * tau).exp();
        let mut pairs = Vec::new();
        let mut inv_phi = Vec::new();
        for ni in 0..gn.spec_len() {
            if !gn.kept(ni) {
                continue;
            }
            let (i, j, l) = gn.spec_ijl(ni);
            let (kx, ky, kz) = (wavenumber(i, gn.n), wavenumber(j, gn.n), l as i64);
            let w = |k: i64| k.rem_euclid(gm.n as i64) as usize;
            let mi = (w(kx) * gm.n + w(ky)) * gm.nh + l;
            pairs.push((mi, ni));
            inv_phi.push(1.0 / (phi(kx as f64) * phi(ky as f64) * phi(kz as f64)));
        }
        let (ns, nr) = (gm.spec_len(), gm.real_len());
        Self {
            m: gm.n,
            tau,
            half_width,
            fft: Fft3::new(&gm),
            pairs,
            inv_phi,
            work_r: vec![0.0; nr],
            work_c: vec![Complex64::default(); ns],
        }
    }

    fn h(&self) -> f64 {
        2.0 * std::f64::consts::PI / self.m as f64
    }

    /// For point coordinate x: the first stencil index and the 1D Gaussian
    /// weights at the 2w consecutive grid points (periodic).
    fn stencil(&self, x: f64) -> (i64, Vec<f64>) {
        let h = self.h();
        let w = self.half_width as i64;
        let base = (x / h).floor() as i64 - w + 1;
        let weights = (0..2 * w)
            .map(|q| {
                let d = x - (base + q) as f64 * h;
                (-d * d / (4.0 * self.tau)).exp()
            })
            .collect();
        (base, weights)
    }

    /// Type 2: values of the field with coefficients `uh` (N-grid spectrum,
    /// retained modes) at `points`.
    pub fn interpolate(&mut self, uh: &Spec, points: &[[f64; 3]]) -> Vec<f64> {
        self.work_c.iter_mut().for_each(|v| *v = Complex64::default());
        for (q, &(mi, ni)) in self.pairs.iter().enumerate() {
            self.work_c[mi] = uh[ni] * self.inv_phi[q];
        }
        self.fft.inverse_inplace(&mut self.work_c, &mut self.work_r);
        let (m, w2) = (self.m as i64, 2 * self.half_width);
        let h3 = self.h().powi(3);
        let wr = &self.work_r;
        points
            .par_iter()
            .map(|x| {
                let (bx, wx) = self.stencil(x[0]);
                let (by, wy) = self.stencil(x[1]);
                let (bz, wz) = self.stencil(x[2]);
                let mut acc = 0.0;
                for a in 0..w2 {
                    let i = (bx + a as i64).rem_euclid(m) as usize;
                    for b in 0..w2 {
                        let j = (by + b as i64).rem_euclid(m) as usize;
                        let row = (i * m as usize + j) * m as usize;
                        let wab = wx[a] * wy[b];
                        for c in 0..w2 {
                            let l = (bz + c as i64).rem_euclid(m) as usize;
                            acc += wab * wz[c] * wr[row + l];
                        }
                    }
                }
                acc * h3
            })
            .collect()
    }

    /// Type 1: retained-mode coefficients of Σ_m F_m δ(x − X_m), written into
    /// `out` (N-grid spectrum; other entries zeroed).
    pub fn spread(&mut self, forces: &[f64], points: &[[f64; 3]], out: &mut Spec) {
        let (m, w2) = (self.m as i64, 2 * self.half_width);
        let mu = self.m;
        // serial accumulation (markers are few relative to grid points)
        self.work_r.iter_mut().for_each(|v| *v = 0.0);
        for (x, &f) in points.iter().zip(forces) {
            let (bx, wx) = self.stencil(x[0]);
            let (by, wy) = self.stencil(x[1]);
            let (bz, wz) = self.stencil(x[2]);
            for a in 0..w2 {
                let i = (bx + a as i64).rem_euclid(m) as usize;
                for b in 0..w2 {
                    let j = (by + b as i64).rem_euclid(m) as usize;
                    let row = (i * mu + j) * mu;
                    let wab = f * wx[a] * wy[b];
                    for c in 0..w2 {
                        let l = (bz + c as i64).rem_euclid(m) as usize;
                        self.work_r[row + l] += wab * wz[c];
                    }
                }
            }
        }
        self.fft.forward(&self.work_r, &mut self.work_c);
        out.iter_mut().for_each(|v| *v = Complex64::default());
        for (q, &(mi, ni)) in self.pairs.iter().enumerate() {
            out[ni] = self.work_c[mi] * self.inv_phi[q];
        }
    }
}
