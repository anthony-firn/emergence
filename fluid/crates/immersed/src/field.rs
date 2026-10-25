//! Grid geometry and padded arrays for the staggered (MAC) layout.

use geometry::P3;

/// A uniform Cartesian grid of cubic cells, `n[0] × n[1] × n[2]`, low
/// corner at `origin`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Grid {
    pub n: [usize; 3],
    pub h: f64,
    pub origin: P3,
}

impl Grid {
    /// Grid covering `[lo, lo + n h]`.
    pub fn new(n: [usize; 3], h: f64, origin: P3) -> Self {
        Self { n, h, origin }
    }

    /// Node counts of velocity component `c`: one extra along its own axis.
    pub fn comp_dims(&self, c: usize) -> [usize; 3] {
        std::array::from_fn(|a| self.n[a] + usize::from(a == c))
    }

    /// Position of node `(i, j, k)` of component `c` (`c = 3`: cell centre).
    #[inline]
    pub fn node(&self, c: usize, i: isize, j: isize, k: isize) -> P3 {
        let idx = [i, j, k];
        std::array::from_fn(|a| self.origin[a] + self.h * (idx[a] as f64 + if a == c { 0.0 } else { 0.5 }))
    }

    pub fn cells(&self) -> usize {
        self.n.iter().product()
    }

    pub fn hi(&self) -> P3 {
        std::array::from_fn(|a| self.origin[a] + self.h * self.n[a] as f64)
    }

    /// Spatial dimension (axes with more than one cell).
    pub fn dim(&self) -> usize {
        self.n.iter().filter(|&&n| n > 1).count()
    }
}

/// A 3D array with one ghost layer on every side, `x` fastest.
#[derive(Clone, Debug, PartialEq)]
pub struct Arr {
    /// Interior extents.
    pub d: [usize; 3],
    pub data: Vec<f64>,
    sx: usize,
    sxy: usize,
}

impl Arr {
    pub fn new(d: [usize; 3]) -> Self {
        let sx = d[0] + 2;
        let sxy = sx * (d[1] + 2);
        Self {
            d,
            data: vec![0.0; sxy * (d[2] + 2)],
            sx,
            sxy,
        }
    }

    /// Flat index of `(i, j, k)`, `−1 ≤ i ≤ d[0]` etc.
    #[inline(always)]
    pub fn at(&self, i: isize, j: isize, k: isize) -> usize {
        (i + 1) as usize + self.sx * (j + 1) as usize + self.sxy * (k + 1) as usize
    }

    #[inline(always)]
    pub fn get(&self, i: isize, j: isize, k: isize) -> f64 {
        self.data[self.at(i, j, k)]
    }

    #[inline(always)]
    pub fn set(&mut self, i: isize, j: isize, k: isize, v: f64) {
        let q = self.at(i, j, k);
        self.data[q] = v;
    }

    /// Strides `(1, sx, sxy)` of the flat layout.
    #[inline(always)]
    pub fn strides(&self) -> [usize; 3] {
        [1, self.sx, self.sxy]
    }

    /// Inverse of [`Arr::at`].
    pub fn ijk(&self, q: usize) -> [isize; 3] {
        [(q % self.sx) as isize - 1, ((q / self.sx) % (self.d[1] + 2)) as isize - 1, (q / self.sxy) as isize - 1]
    }

    /// Trilinear interpolation at fractional index `f` (node units),
    /// clamped to the padded array. Returns the 8 (index, weight) pairs.
    pub fn stencil(&self, f: P3) -> [(usize, f64); 8] {
        let mut base = [0isize; 3];
        let mut t = [0.0; 3];
        for a in 0..3 {
            if self.d[a] == 1 {
                // A single node along this axis (2D slab): the coordinate is
                // irrelevant; never let round-off pick the ghost layer.
                base[a] = 0;
                t[a] = 0.0;
                continue;
            }
            let lo = -1.0;
            let hi = self.d[a] as f64; // last ghost index
            let x = f[a].clamp(lo, hi);
            let mut b = x.floor() as isize;
            if b as f64 >= hi {
                b = hi as isize - 1;
            }
            base[a] = b;
            t[a] = x - b as f64;
        }
        let mut out = [(0usize, 0.0); 8];
        for (m, o) in out.iter_mut().enumerate() {
            let (di, dj, dk) = ((m & 1) as isize, ((m >> 1) & 1) as isize, ((m >> 2) & 1) as isize);
            let w = (if di == 1 { t[0] } else { 1.0 - t[0] }) * (if dj == 1 { t[1] } else { 1.0 - t[1] }) * (if dk == 1 { t[2] } else { 1.0 - t[2] });
            *o = (self.at(base[0] + di, base[1] + dj, base[2] + dk), w);
        }
        out
    }

    #[inline]
    pub fn interp(&self, stencil: &[(usize, f64); 8]) -> f64 {
        stencil.iter().map(|&(q, w)| w * self.data[q]).sum()
    }
}
