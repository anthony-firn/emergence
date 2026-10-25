//! Periodic grid on [0, 2π)³ and its Fourier dual.
//!
//! Layouts (shared with the CUDA backend):
//! - physical: `N*N*N` reals, index `(i*N + j)*N + l`, `l` (z) fastest.
//! - spectral: `N*N*NH` complex, `NH = N/2 + 1`, index `(i*N + j)*NH + l`,
//!   the real-to-complex half-spectrum stored along z.

/// Integer wavenumber of FFT index `i` on an `n`-point axis.
#[inline]
pub fn wavenumber(i: usize, n: usize) -> i64 {
    if i <= n / 2 {
        i as i64
    } else {
        i as i64 - n as i64
    }
}

#[derive(Clone, Debug)]
pub struct Grid {
    pub n: usize,
    pub nh: usize,
    /// Wavenumbers along x and y (full axes) and z (half axis).
    pub kx: Vec<f64>,
    pub kz: Vec<f64>,
    /// 2/3-rule dealiasing: modes with 3|k_i| >= N (and all Nyquist modes) are zeroed.
    pub keep_x: Vec<bool>,
    pub keep_z: Vec<bool>,
    /// Spherical truncation |k|² < kc² with kc = √2·N/3 (used with
    /// phase-shift dealiasing) instead of the 2/3 cube.
    pub spherical: bool,
    pub kc2: f64,
}

impl Grid {
    pub fn new(n: usize) -> Self {
        assert!(n >= 8 && n % 2 == 0, "grid size must be even and >= 8");
        Self::new_any(n)
    }

    /// Grid for phase-shift dealiasing (Rogallo 1981): spherical truncation
    /// at kc = √2·N/3. Single and triple aliases cancel between a random grid
    /// shift Δ and Δ + half a cell; double aliases need |p|,|q| < kc with
    /// p + q aliased in two directions, which this truncation excludes.
    /// Retains ~41% larger k_max than the 2/3 cube at the same N.
    pub fn new_spherical(n: usize) -> Self {
        let mut g = Self::new(n);
        let kc = std::f64::consts::SQRT_2 * n as f64 / 3.0;
        g.spherical = true;
        g.kc2 = kc * kc;
        let nyq = (n / 2) as f64;
        g.keep_x = g.kx.iter().map(|&k| k.abs() < nyq && k.abs() < kc).collect();
        g.keep_z = g.kz.iter().map(|&k| k.abs() < nyq && k.abs() < kc).collect();
        g
    }

    /// Grid of any size ≥ 3, including odd sizes (used for the penalization
    /// solution grid). Odd sizes have no Nyquist mode.
    pub fn new_any(n: usize) -> Self {
        assert!(n >= 3);
        let nh = n / 2 + 1;
        let kx: Vec<f64> = (0..n).map(|i| wavenumber(i, n) as f64).collect();
        let kz: Vec<f64> = (0..nh).map(|l| l as f64).collect();
        let keep = |k: f64| 3.0 * k.abs() < n as f64 && k.abs() < (n / 2) as f64;
        Self {
            n,
            nh,
            keep_x: kx.iter().map(|&k| keep(k)).collect(),
            keep_z: kz.iter().map(|&k| keep(k)).collect(),
            spherical: false,
            kc2: 0.0,
            kx,
            kz,
        }
    }

    pub fn real_len(&self) -> usize {
        self.n * self.n * self.n
    }

    pub fn spec_len(&self) -> usize {
        self.n * self.n * self.nh
    }

    pub fn dx(&self) -> f64 {
        2.0 * std::f64::consts::PI / self.n as f64
    }

    /// Coordinates of physical point `(i, j, l)`.
    #[inline]
    pub fn x(&self, i: usize) -> f64 {
        i as f64 * self.dx()
    }

    /// Split a spectral index into `(i, j, l)`.
    #[inline]
    pub fn spec_ijl(&self, idx: usize) -> (usize, usize, usize) {
        let l = idx % self.nh;
        let ij = idx / self.nh;
        (ij / self.n, ij % self.n, l)
    }

    /// Wavevector of spectral index.
    #[inline]
    pub fn k(&self, idx: usize) -> [f64; 3] {
        let (i, j, l) = self.spec_ijl(idx);
        [self.kx[i], self.kx[j], self.kz[l]]
    }

    #[inline]
    pub fn kept(&self, idx: usize) -> bool {
        let (i, j, l) = self.spec_ijl(idx);
        let cube = self.keep_x[i] && self.keep_x[j] && self.keep_z[l];
        if !self.spherical {
            return cube;
        }
        let (a, b, c) = (self.kx[i], self.kx[j], self.kz[l]);
        cube && a * a + b * b + c * c < self.kc2
    }

    /// Largest retained |k| per axis (2/3 rule: 3|k| < N; spherical: < kc).
    pub fn kmax(&self) -> usize {
        if self.spherical {
            let kc = self.kc2.sqrt();
            ((kc - 1e-9).floor() as usize).min((self.n - 1) / 2)
        } else {
            (self.n - 1) / 3
        }
    }

    /// Multiplicity of a half-spectrum entry in full-spectrum sums
    /// (the conjugate partner at -k is not stored for 0 < kz < N/2).
    #[inline]
    pub fn weight(&self, idx: usize) -> f64 {
        let l = idx % self.nh;
        if l == 0 || (self.n % 2 == 0 && l == self.n / 2) {
            1.0
        } else {
            2.0
        }
    }
}
