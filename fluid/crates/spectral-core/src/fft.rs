//! Parallel 3D real-to-complex FFT on the layouts described in `grid.rs`.
//!
//! The forward transform is normalized by 1/N³, so spectral arrays hold
//! Fourier coefficients: u(x) = Σ_k û_k e^{ik·x}.

use std::sync::Arc;

use num_complex::Complex64;
use rayon::prelude::*;
use realfft::{ComplexToReal, RealFftPlanner, RealToComplex};
use rustfft::{Fft, FftPlanner};

use crate::grid::Grid;

/// Raw pointer that may be shared across rayon tasks touching disjoint indices.
#[derive(Clone, Copy)]
struct SyncPtr(*mut Complex64);
unsafe impl Send for SyncPtr {}
unsafe impl Sync for SyncPtr {}

pub struct Fft3 {
    n: usize,
    nh: usize,
    r2c: Arc<dyn RealToComplex<f64>>,
    c2r: Arc<dyn ComplexToReal<f64>>,
    fwd: Arc<dyn Fft<f64>>,
    inv: Arc<dyn Fft<f64>>,
}

impl Fft3 {
    pub fn new(grid: &Grid) -> Self {
        let n = grid.n;
        let mut rp = RealFftPlanner::<f64>::new();
        let mut cp = FftPlanner::<f64>::new();
        Self {
            n,
            nh: grid.nh,
            r2c: rp.plan_fft_forward(n),
            c2r: rp.plan_fft_inverse(n),
            fwd: cp.plan_fft_forward(n),
            inv: cp.plan_fft_inverse(n),
        }
    }

    /// Complex transform along y (stride NH) inside each x-plane, then along x
    /// (stride N*NH). `forward` selects the direction.
    fn xy_pass(&self, spec: &mut [Complex64], forward: bool) {
        let (n, nh) = (self.n, self.nh);
        let plan = if forward { &self.fwd } else { &self.inv };
        let scratch_len = plan.get_inplace_scratch_len();

        // y: each x-plane is a contiguous block of N*NH.
        spec.par_chunks_mut(n * nh).for_each_init(
            || (vec![Complex64::default(); n], vec![Complex64::default(); scratch_len]),
            |(line, scratch), plane| {
                for l in 0..nh {
                    for j in 0..n {
                        line[j] = plane[j * nh + l];
                    }
                    plan.process_with_scratch(line, scratch);
                    for j in 0..n {
                        plane[j * nh + l] = line[j];
                    }
                }
            },
        );

        // x: columns (j, l) are disjoint across tasks.
        let ptr = SyncPtr(spec.as_mut_ptr());
        (0..n).into_par_iter().for_each_init(
            || (vec![Complex64::default(); n], vec![Complex64::default(); scratch_len]),
            |(line, scratch), j| {
                let p = ptr;
                for l in 0..nh {
                    for i in 0..n {
                        // SAFETY: index (i, j, l) is touched only by the task owning j.
                        line[i] = unsafe { *p.0.add((i * n + j) * nh + l) };
                    }
                    plan.process_with_scratch(line, scratch);
                    for i in 0..n {
                        unsafe { *p.0.add((i * n + j) * nh + l) = line[i] };
                    }
                }
            },
        );
    }

    /// physical -> spectral (normalized by 1/N³).
    pub fn forward(&self, real: &[f64], spec: &mut [Complex64]) {
        let (n, nh) = (self.n, self.nh);
        debug_assert_eq!(real.len(), n * n * n);
        debug_assert_eq!(spec.len(), n * n * nh);
        let norm = 1.0 / (n * n * n) as f64;
        let scratch_len = self.r2c.get_scratch_len();
        real.par_chunks(n).zip(spec.par_chunks_mut(nh)).for_each_init(
            || (vec![0.0; n], vec![Complex64::default(); scratch_len]),
            |(buf, scratch), (src, dst)| {
                buf.copy_from_slice(src);
                self.r2c.process_with_scratch(buf, dst, scratch).expect("r2c");
                for v in dst.iter_mut() {
                    *v *= norm;
                }
            },
        );
        self.xy_pass(spec, true);
    }

    /// spectral -> physical. `spec` is used as workspace and left modified.
    pub fn inverse_inplace(&self, spec: &mut [Complex64], real: &mut [f64]) {
        let (n, nh) = (self.n, self.nh);
        self.xy_pass(spec, false);
        let scratch_len = self.c2r.get_scratch_len();
        spec.par_chunks_mut(nh).zip(real.par_chunks_mut(n)).for_each_init(
            || vec![Complex64::default(); scratch_len],
            |scratch, (src, dst)| {
                // A real field has real DC (and, for even n, Nyquist) entries
                // along z; drop roundoff. For odd n the last entry is an
                // ordinary complex mode and must be left alone.
                src[0].im = 0.0;
                if n % 2 == 0 {
                    src[nh - 1].im = 0.0;
                }
                self.c2r.process_with_scratch(src, dst, scratch).expect("c2r");
            },
        );
    }

    /// spectral -> physical, preserving the input (copies into `work`).
    pub fn inverse(&self, spec: &[Complex64], work: &mut [Complex64], real: &mut [f64]) {
        work.par_iter_mut().zip(spec.par_iter()).for_each(|(w, s)| *w = *s);
        self.inverse_inplace(work, real);
    }
}
