//! Exact fast Poisson solver for the cell-centred box Laplacian.
//!
//! With the body immersed, the pressure equation is the plain 7-point
//! Laplacian on the box, with one condition per box face:
//!
//! * Neumann (`∂φ/∂n = 0`: walls, velocity inlets, slip), ghost `φ₋₁ = φ₀`;
//! * Dirichlet (`φ = 0` on the face: pressure outlets), ghost `φ₋₁ = −φ₀`;
//! * periodic.
//!
//! The 1D second-difference matrix with each pair of end conditions is
//! diagonalised by a real trigonometric transform, so the 3D operator is
//! inverted exactly (to round-off) by one transform per axis, a diagonal
//! division and the inverse transforms, in `O(N log N)` for any grid size:
//!
//! | ends | basis `v_k(i)` | eigenvalue `× h²` |
//! |---|---|---|
//! | N–N | DCT-II `cos(πk(i+½)/n)` | `−4 sin²(πk/2n)` |
//! | D–D | DST-II `sin(π(k+1)(i+½)/n)` | `−4 sin²(π(k+1)/2n)` |
//! | N–D | DCT-IV `cos(π(k+½)(i+½)/n)` | `−4 sin²(π(2k+1)/4n)` |
//! | D–N | DCT-IV on the reversed line | same |
//! | periodic | Hartley `cas(2πik/n)` | `−4 sin²(πk/n)` |
//!
//! Every transform is computed with a complex FFT of length `n` or `2n`.

use std::f64::consts::PI;
use std::sync::Arc;

use rayon::prelude::*;
use rustfft::num_complex::Complex64 as Cx;
use rustfft::{Fft, FftPlanner};

/// Condition at one end of an axis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum End {
    Neumann,
    Dirichlet,
    Periodic,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    NN,
    DD,
    ND,
    DN,
    PP,
}

/// Per-thread FFT work buffers.
#[derive(Default)]
struct Work {
    buf: Vec<Cx>,
    scratch: Vec<Cx>,
}

struct Axis1 {
    n: usize,
    kind: Kind,
    fwd: Arc<dyn Fft<f64>>,
    inv: Arc<dyn Fft<f64>>,
    /// Eigenvalues of the 1D operator (already divided by h²).
    lambda: Vec<f64>,
}

impl Axis1 {
    fn new(n: usize, lo: End, hi: End, h: f64, planner: &mut FftPlanner<f64>) -> Self {
        let kind = match (lo, hi) {
            (End::Neumann, End::Neumann) => Kind::NN,
            (End::Dirichlet, End::Dirichlet) => Kind::DD,
            (End::Neumann, End::Dirichlet) => Kind::ND,
            (End::Dirichlet, End::Neumann) => Kind::DN,
            (End::Periodic, End::Periodic) => Kind::PP,
            _ => panic!("periodic must be paired with periodic"),
        };
        let len = if kind == Kind::PP { n } else { 2 * n };
        let s = |x: f64| x.sin().powi(2);
        let lambda = (0..n)
            .map(|k| {
                let k = k as f64;
                let n = n as f64;
                -4.0 / (h * h)
                    * match kind {
                        Kind::NN => s(PI * k / (2.0 * n)),
                        Kind::DD => s(PI * (k + 1.0) / (2.0 * n)),
                        Kind::ND | Kind::DN => s(PI * (2.0 * k + 1.0) / (4.0 * n)),
                        Kind::PP => s(PI * k / n),
                    }
            })
            .collect();
        Self {
            n,
            kind,
            fwd: planner.plan_fft_forward(len),
            inv: planner.plan_fft_inverse(len),
            lambda,
        }
    }

    fn scratch_len(&self) -> usize {
        self.fwd.get_inplace_scratch_len().max(self.inv.get_inplace_scratch_len())
    }

    fn forward(&self, x: &mut [f64], buf: &mut Work) {
        match self.kind {
            Kind::NN => self.dct2(x, buf),
            Kind::DD => {
                for (i, v) in x.iter_mut().enumerate() {
                    if i % 2 == 1 {
                        *v = -*v;
                    }
                }
                self.dct2(x, buf);
                x.reverse();
            }
            Kind::ND => self.dct4(x, buf),
            Kind::DN => {
                x.reverse();
                self.dct4(x, buf);
            }
            Kind::PP => self.dht(x, buf),
        }
    }

    fn inverse(&self, x: &mut [f64], buf: &mut Work) {
        let n = self.n as f64;
        match self.kind {
            Kind::NN => self.idct2(x, buf),
            Kind::DD => {
                x.reverse();
                self.idct2(x, buf);
                for (i, v) in x.iter_mut().enumerate() {
                    if i % 2 == 1 {
                        *v = -*v;
                    }
                }
            }
            Kind::ND => {
                self.dct4(x, buf);
                x.iter_mut().for_each(|v| *v *= 2.0 / n);
            }
            Kind::DN => {
                self.dct4(x, buf);
                x.iter_mut().for_each(|v| *v *= 2.0 / n);
                x.reverse();
            }
            Kind::PP => {
                self.dht(x, buf);
                x.iter_mut().for_each(|v| *v /= n);
            }
        }
    }

    /// `y_k = Σ x_i cos(πk(2i+1)/2n)` = `Re[e^{−iπk/2n} FFT_2n(x, 0…)_k]`.
    fn dct2(&self, x: &mut [f64], w: &mut Work) {
        let Work { buf, scratch } = w;
        let n = self.n;
        buf.clear();
        buf.extend(x.iter().map(|&v| Cx::new(v, 0.0)));
        buf.resize(2 * n, Cx::new(0.0, 0.0));
        self.fwd.process_with_scratch(buf, scratch);
        for (k, y) in x.iter_mut().enumerate() {
            let a = -PI * k as f64 / (2.0 * n as f64);
            *y = (Cx::new(a.cos(), a.sin()) * buf[k]).re;
        }
    }

    /// Inverse of [`Self::dct2`]: `x_i = (1/n)[y₀ + 2 Σ_{k≥1} y_k cos(πk(2i+1)/2n)]`.
    fn idct2(&self, x: &mut [f64], w: &mut Work) {
        let Work { buf, scratch } = w;
        let n = self.n;
        buf.clear();
        buf.extend(x.iter().enumerate().map(|(k, &y)| {
            let a = PI * k as f64 / (2.0 * n as f64);
            Cx::new(a.cos(), a.sin()) * (if k == 0 { y } else { 2.0 * y })
        }));
        buf.resize(2 * n, Cx::new(0.0, 0.0));
        self.inv.process_with_scratch(buf, scratch);
        for (i, v) in x.iter_mut().enumerate() {
            *v = buf[i].re / n as f64;
        }
    }

    /// `y_k = Σ x_i cos(π(2k+1)(2i+1)/4n)` (self-inverse up to `2/n`).
    fn dct4(&self, x: &mut [f64], w: &mut Work) {
        let Work { buf, scratch } = w;
        let n = self.n;
        buf.clear();
        buf.extend(x.iter().enumerate().map(|(i, &v)| {
            let a = -PI * i as f64 / (2.0 * n as f64);
            Cx::new(a.cos(), a.sin()) * v
        }));
        buf.resize(2 * n, Cx::new(0.0, 0.0));
        self.fwd.process_with_scratch(buf, scratch);
        for (k, y) in x.iter_mut().enumerate() {
            let a = -PI * (2 * k + 1) as f64 / (4.0 * n as f64);
            *y = (Cx::new(a.cos(), a.sin()) * buf[k]).re;
        }
    }

    /// Discrete Hartley transform `H_k = Σ x_i cas(2πik/n)` = `Re X_k − Im X_k`.
    fn dht(&self, x: &mut [f64], w: &mut Work) {
        let Work { buf, scratch } = w;
        buf.clear();
        buf.extend(x.iter().map(|&v| Cx::new(v, 0.0)));
        self.fwd.process_with_scratch(buf, scratch);
        for (k, y) in x.iter_mut().enumerate() {
            *y = buf[k].re - buf[k].im;
        }
    }
}

/// The box Laplacian inverse on an `n[0] × n[1] × n[2]` cell grid (x fastest).
pub struct Poisson {
    pub n: [usize; 3],
    pub h: f64,
    axes: [Axis1; 3],
    ends: [[End; 2]; 3],
}

impl Poisson {
    pub fn new(n: [usize; 3], h: f64, ends: [[End; 2]; 3]) -> Self {
        let mut planner = FftPlanner::new();
        let axes = std::array::from_fn(|a| Axis1::new(n[a], ends[a][0], ends[a][1], h, &mut planner));
        Self { n, h, axes, ends }
    }

    /// True when the operator has a null space (constants).
    pub fn singular(&self) -> bool {
        self.ends.iter().all(|e| e[0] != End::Dirichlet && e[1] != End::Dirichlet)
    }

    /// Apply one transform along axis `a` to every line of `f`, in parallel
    /// over lines (x lines in place; y and z lines gathered, transformed
    /// and scattered back).
    fn sweep(&self, f: &mut [f64], a: usize, forward: bool) {
        let [nx, ny, nz] = self.n;
        let ax = &self.axes[a];
        let init = || Work {
            buf: Vec::new(),
            scratch: vec![Cx::new(0.0, 0.0); ax.scratch_len()],
        };
        let run = |line: &mut [f64], w: &mut Work| {
            if forward {
                ax.forward(line, w)
            } else {
                ax.inverse(line, w)
            }
        };
        if a == 0 {
            f.par_chunks_mut(nx).for_each_init(init, |w, line| run(line, w));
            return;
        }
        // Lines along y (stride nx) or z (stride nx·ny), indexed by their
        // starting offset.
        let (len, stride) = if a == 1 { (ny, nx) } else { (nz, nx * ny) };
        let starts: Vec<usize> = if a == 1 {
            (0..nz).flat_map(|k| (0..nx).map(move |i| i + nx * ny * k)).collect()
        } else {
            (0..nx * ny).collect()
        };
        let src: &[f64] = f;
        let lines: Vec<Vec<f64>> = starts
            .par_iter()
            .map_init(init, |w, &s0| {
                let mut line: Vec<f64> = (0..len).map(|m| src[s0 + stride * m]).collect();
                run(&mut line, w);
                line
            })
            .collect();
        for (s0, line) in starts.iter().zip(&lines) {
            for (m, v) in line.iter().enumerate() {
                f[s0 + stride * m] = *v;
            }
        }
    }

    /// Solve `L φ = rhs` in place. For a singular operator the constant
    /// mode of the right-hand side is dropped (its mean is the
    /// incompatibility; it is returned) and `φ` has zero mean.
    pub fn solve(&self, f: &mut [f64]) -> f64 {
        let [nx, ny, nz] = self.n;
        assert_eq!(f.len(), nx * ny * nz);
        let mean = f.iter().sum::<f64>() / f.len() as f64;
        for a in 0..3 {
            if self.n[a] > 1 || self.axes[a].kind != Kind::PP {
                self.sweep(f, a, true);
            }
        }
        let (lx, ly, lz) = (&self.axes[0].lambda, &self.axes[1].lambda, &self.axes[2].lambda);
        f.par_chunks_mut(nx * ny).enumerate().for_each(|(k, pl)| {
            for j in 0..ny {
                for i in 0..nx {
                    let l = lx[i] + ly[j] + lz[k];
                    let v = &mut pl[i + nx * j];
                    *v = if l != 0.0 { *v / l } else { 0.0 };
                }
            }
        });
        for a in (0..3).rev() {
            if self.n[a] > 1 || self.axes[a].kind != Kind::PP {
                self.sweep(f, a, false);
            }
        }
        if self.singular() {
            mean
        } else {
            0.0
        }
    }

    /// Apply the discrete operator (for tests and residual checks).
    pub fn apply(&self, phi: &[f64]) -> Vec<f64> {
        let [nx, ny, nz] = self.n;
        let h2 = self.h * self.h;
        let get = |i: isize, j: isize, k: isize, c: f64| -> f64 {
            // Reflect into the box according to the end conditions; `c` is
            // the value at the cell the ghost mirrors.
            let idx = [i, j, k];
            let n = [nx as isize, ny as isize, nz as isize];
            let mut q = idx;
            let mut sign = 1.0;
            for a in 0..3 {
                if q[a] < 0 || q[a] >= n[a] {
                    let end = self.ends[a][if q[a] < 0 { 0 } else { 1 }];
                    match end {
                        End::Periodic => q[a] = q[a].rem_euclid(n[a]),
                        End::Neumann => return c,
                        End::Dirichlet => {
                            sign = -1.0;
                            return sign * c;
                        }
                    }
                }
            }
            sign * phi[(q[0] + n[0] * (q[1] + n[1] * q[2])) as usize]
        };
        let mut out = vec![0.0; phi.len()];
        out.par_iter_mut().enumerate().for_each(|(q, o)| {
            let (i, j, k) = ((q % nx) as isize, ((q / nx) % ny) as isize, (q / (nx * ny)) as isize);
            let c = phi[q];
            let mut s = -6.0 * c;
            for (di, dj, dk) in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)] {
                s += get(i + di, j + dj, k + dk, c);
            }
            *o = s / h2;
        });
        out
    }
}
