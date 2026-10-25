//! Three-dimensional finite-difference Landau-Lifshitz-Gilbert solver with
//! the FFT demagnetising field.
//!
//! This is the unfold of the local-shape demagnetising approximation in
//! [`crate::micromag`]: the exchange, anisotropy and Zeeman terms are the
//! same, but the magnetostatic field is the exact non-local convolution
//! computed by [`crate::demag::DemagFft`].  The LLG
//!
//! ```text
//! dm/dt = -gamma/(1+alpha^2) [ m x H_eff + alpha m x (m x H_eff) ]
//! ```
//!
//! is integrated with RK4 and renormalisation.  Exchange uses mirror
//! (zero-derivative) boundaries, the natural condition for a free surface.
//!
//! All per-step buffers (the RK4 stages, the effective-field scratch and the
//! FFT workspace) live in a caller-owned [`LlgScratch`], so a long reversal
//! does not allocate in its hot loop; the exchange and Gilbert terms are
//! evaluated with rayon.  The `*_ws` methods are the allocation-free path;
//! the plain methods are thin convenience wrappers used by tests.

use crate::demag::{DemagFft, FftWorkspace};
use crate::geometry::{dot, V3};
use crate::micromag::{normalize, GAMMA, MU0};
use rayon::prelude::*;

fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Reusable per-step buffers for [`Llg3d`].  Sized for one grid; the FFT
/// workspace is present only when the solver has a [`DemagFft`] attached.
pub struct LlgScratch {
    h: Vec<V3>,
    k1: Vec<V3>,
    k2: Vec<V3>,
    k3: Vec<V3>,
    k4: Vec<V3>,
    s2: Vec<V3>,
    s3: Vec<V3>,
    s4: Vec<V3>,
    fft: Option<FftWorkspace>,
}

impl LlgScratch {
    pub fn new(n_cells: usize, fft: Option<FftWorkspace>) -> Self {
        Self {
            h: vec![[0.0; 3]; n_cells],
            k1: vec![[0.0; 3]; n_cells],
            k2: vec![[0.0; 3]; n_cells],
            k3: vec![[0.0; 3]; n_cells],
            k4: vec![[0.0; 3]; n_cells],
            s2: vec![[0.0; 3]; n_cells],
            s3: vec![[0.0; 3]; n_cells],
            s4: vec![[0.0; 3]; n_cells],
            fft,
        }
    }
}

/// 3-D micromagnetic state on an `nx x ny x nz` grid of rectangular cells.
pub struct Llg3d {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub dx: f64,
    pub dy: f64,
    pub dz: f64,
    pub a: f64,
    pub k: f64,
    pub ms: f64,
    pub alpha: f64,
    pub easy: V3,
    pub h_ext: V3,
    pub demag: Option<DemagFft>,
    pub m: Vec<V3>,
}

impl Llg3d {
    pub fn new(
        nx: usize,
        ny: usize,
        nz: usize,
        dx: f64,
        dy: f64,
        dz: f64,
        a: f64,
        k: f64,
        ms: f64,
        alpha: f64,
    ) -> Self {
        Self {
            nx,
            ny,
            nz,
            dx,
            dy,
            dz,
            a,
            k,
            ms,
            alpha,
            easy: [1.0, 0.0, 0.0],
            h_ext: [0.0, 0.0, 0.0],
            demag: None,
            m: vec![[1.0, 0.0, 0.0]; nx * ny * nz],
        }
    }

    /// Allocate the reusable step buffers for this grid.
    pub fn make_scratch(&self) -> LlgScratch {
        let fft = self
            .demag
            .as_ref()
            .map(|d| FftWorkspace::new(d.padded_len()));
        LlgScratch::new(self.n_cells(), fft)
    }

    pub fn n_cells(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    #[inline]
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        (i * self.ny + j) * self.nz + k
    }

    pub fn set_uniform(&mut self, dir: V3) {
        let d = normalize(dir);
        for x in self.m.iter_mut() {
            *x = d;
        }
    }

    /// Exchange coefficient `2A/(mu0 M_s)` in m^2.
    pub fn exchange_coeff(&self) -> f64 {
        2.0 * self.a / (MU0 * self.ms)
    }

    /// Effective field (A/m) at every cell.  Exchange uses mirror
    /// boundaries; `out` gets the grid ordering of [`Llg3d`].
    pub fn effective_field(&self, m: &[V3], out: &mut [V3]) {
        let mut fft = self
            .demag
            .as_ref()
            .map(|d| FftWorkspace::new(d.padded_len()));
        self.effective_field_ws(m, out, &mut fft);
    }

    /// Allocation-free effective field using a caller-supplied FFT
    /// workspace.  The demagnetising field is written into `out` first and
    /// then the exchange, anisotropy and Zeeman terms are added in a single
    /// parallel pass, so no separate demag buffer is needed.
    pub fn effective_field_ws(&self, m: &[V3], out: &mut [V3], fft: &mut Option<FftWorkspace>) {
        let demag_scale = if let Some(d) = &self.demag {
            let ws = fft
                .as_mut()
                .expect("demag workspace missing while demag field is attached");
            d.field_ws(m, out, ws);
            self.ms
        } else {
            0.0
        };
        let ex = self.exchange_coeff();
        let ani = 2.0 * self.k / (MU0 * self.ms);
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);
        let cx = ex / (self.dx * self.dx);
        let cy = ex / (self.dy * self.dy);
        let cz = ex / (self.dz * self.dz);

        out.par_chunks_mut(ny * nz)
            .enumerate()
            .for_each(|(i, orow)| {
                let im = i.saturating_sub(1);
                let ip = (i + 1).min(nx - 1);
                for j in 0..ny {
                    let jm = j.saturating_sub(1);
                    let jp = (j + 1).min(ny - 1);
                    for k in 0..nz {
                        let km = k.saturating_sub(1);
                        let kp = (k + 1).min(nz - 1);
                        let c = self.index(i, j, k);
                        let cl = j * nz + k;
                        let mc = m[c];
                        let lapx = add(
                            add(m[self.index(im, j, k)], m[self.index(ip, j, k)]),
                            scale(mc, -2.0),
                        );
                        let lapy = add(
                            add(m[self.index(i, jm, k)], m[self.index(i, jp, k)]),
                            scale(mc, -2.0),
                        );
                        let lapz = add(
                            add(m[self.index(i, j, km)], m[self.index(i, j, kp)]),
                            scale(mc, -2.0),
                        );
                        let am = dot(mc, self.easy);
                        let demag = scale(orow[cl], demag_scale);
                        orow[cl] = [
                            demag[0]
                                + cx * lapx[0]
                                + cy * lapy[0]
                                + cz * lapz[0]
                                + ani * am * self.easy[0]
                                + self.h_ext[0],
                            demag[1]
                                + cx * lapx[1]
                                + cy * lapy[1]
                                + cz * lapz[1]
                                + ani * am * self.easy[1]
                                + self.h_ext[1],
                            demag[2]
                                + cx * lapx[2]
                                + cy * lapy[2]
                                + cz * lapz[2]
                                + ani * am * self.easy[2]
                                + self.h_ext[2],
                        ];
                    }
                }
            });
    }

    fn rhs_ws(&self, m: &[V3], out: &mut [V3], tmp_h: &mut [V3], fft: &mut Option<FftWorkspace>) {
        self.effective_field_ws(m, tmp_h, fft);
        let g = GAMMA / (1.0 + self.alpha * self.alpha);
        let alpha = self.alpha;
        out.par_iter_mut().enumerate().for_each(|(c, o)| {
            let s = m[c];
            let hh = tmp_h[c];
            let s_cross_h = cross(s, hh);
            let s_cross_scross_h = cross(s, s_cross_h);
            *o = [
                -g * (s_cross_h[0] + alpha * s_cross_scross_h[0]),
                -g * (s_cross_h[1] + alpha * s_cross_scross_h[1]),
                -g * (s_cross_h[2] + alpha * s_cross_scross_h[2]),
            ];
        });
    }

    /// One RK4 step with renormalisation of the final state, allocating a
    /// scratch buffer (convenience; use [`Llg3d::step_ws`] in a loop).
    pub fn step(&mut self, dt: f64) {
        let mut ws = self.make_scratch();
        self.step_ws(dt, &mut ws);
    }

    /// One allocation-free RK4 step.
    pub fn step_ws(&mut self, dt: f64, ws: &mut LlgScratch) {
        let n = self.n_cells();
        self.rhs_ws(&self.m, &mut ws.k1, &mut ws.h, &mut ws.fft);
        for i in 0..n {
            ws.s2[i] = normalize(add(self.m[i], scale(ws.k1[i], 0.5 * dt)));
        }
        self.rhs_ws(&ws.s2, &mut ws.k2, &mut ws.h, &mut ws.fft);
        for i in 0..n {
            ws.s3[i] = normalize(add(self.m[i], scale(ws.k2[i], 0.5 * dt)));
        }
        self.rhs_ws(&ws.s3, &mut ws.k3, &mut ws.h, &mut ws.fft);
        for i in 0..n {
            ws.s4[i] = normalize(add(self.m[i], scale(ws.k3[i], dt)));
        }
        self.rhs_ws(&ws.s4, &mut ws.k4, &mut ws.h, &mut ws.fft);
        for i in 0..n {
            let mut s = self.m[i];
            for c in 0..3 {
                s[c] +=
                    dt / 6.0 * (ws.k1[i][c] + 2.0 * ws.k2[i][c] + 2.0 * ws.k3[i][c] + ws.k4[i][c]);
            }
            self.m[i] = normalize(s);
        }
    }

    /// Spatially averaged magnetisation.
    pub fn average(&self) -> V3 {
        let n = self.n_cells() as f64;
        let mut a = [0.0; 3];
        for m in &self.m {
            for c in 0..3 {
                a[c] += m[c];
            }
        }
        [a[0] / n, a[1] / n, a[2] / n]
    }

    /// Maximum `|m x H_eff|` over the grid (rad/s once multiplied by
    /// `gamma`); a relaxation/convergence measure.
    pub fn max_torque(&self) -> f64 {
        let mut ws = self.make_scratch();
        self.max_torque_ws(&mut ws)
    }

    /// Allocation-free maximum torque.
    pub fn max_torque_ws(&self, ws: &mut LlgScratch) -> f64 {
        self.effective_field_ws(&self.m, &mut ws.h, &mut ws.fft);
        self.m
            .par_iter()
            .zip(ws.h.par_iter())
            .map(|(m, h)| {
                let t = cross(*m, *h);
                (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt()
            })
            .reduce(|| 0.0, f64::max)
    }

    /// Demagnetising energy (J) of the current state.
    pub fn demag_energy(&self) -> f64 {
        match &self.demag {
            Some(d) => d.energy(&self.m, self.ms, MU0),
            None => 0.0,
        }
    }

    /// Relax with heavy damping until the maximum torque falls below
    /// `tol` or `max_steps` is reached.  Returns the number of steps used.
    pub fn relax(&mut self, dt: f64, tol: f64, max_steps: usize) -> usize {
        let mut ws = self.make_scratch();
        self.relax_ws(dt, tol, max_steps, &mut ws)
    }

    /// Allocation-free relaxation.
    pub fn relax_ws(&mut self, dt: f64, tol: f64, max_steps: usize, ws: &mut LlgScratch) -> usize {
        let saved = self.alpha;
        self.alpha = 1.0;
        let mut used = 0;
        for _ in 0..max_steps {
            used += 1;
            self.step_ws(dt, ws);
            if self.max_torque_ws(ws) < tol {
                break;
            }
        }
        self.alpha = saved;
        used
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_state_is_stationary_without_demag() {
        let mut s = Llg3d::new(4, 4, 1, 1e-8, 1e-8, 3e-9, 1.3e-11, 0.0, 8e5, 0.02);
        s.set_uniform([1.0, 0.0, 0.0]);
        s.step(1e-13);
        let a = s.average();
        assert!((a[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn workspace_step_matches_allocating_step() {
        let (nx, ny, nz) = (8, 6, 2);
        let (dx, dy, dz) = (2e-8, 2.2e-8, 3e-9);
        let mut a = Llg3d::new(nx, ny, nz, dx, dy, dz, 1.3e-11, 0.0, 8e5, 0.02);
        a.demag = Some(DemagFft::new(nx, ny, nz, [dx, dy, dz]));
        a.set_uniform([1.0, 1.0, 1.0]);
        let mut b = Llg3d::new(nx, ny, nz, dx, dy, dz, 1.3e-11, 0.0, 8e5, 0.02);
        b.demag = Some(DemagFft::new(nx, ny, nz, [dx, dy, dz]));
        b.m = a.m.clone();
        a.step(1e-14);
        let mut ws = b.make_scratch();
        b.step_ws(1e-14, &mut ws);
        for (x, y) in a.m.iter().zip(&b.m) {
            for c in 0..3 {
                assert!((x[c] - y[c]).abs() < 1e-13, "{x:?} vs {y:?}");
            }
        }
    }
}
