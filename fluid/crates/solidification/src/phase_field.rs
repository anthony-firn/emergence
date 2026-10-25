//! Anisotropic Karma–Rappel phase-field solidification of a pure substance.
//!
//! # Model
//!
//! The phase field `p = +1` (solid) / `-1` (liquid) and the dimensionless
//! temperature `u = (T - Tm)/(L/c_p)` obey the thin-interface model of
//! Karma & Rappel, Phys. Rev. E **57**, 4323 (1998) (pure-substance limit, so
//! the anti-trapping solute term is absent):
//!
//! ```text
//!   tau(theta) dp/dt = div( W(theta)^2 grad p )
//!                      + (p - lambda u (1 - p^2)) (1 - p^2),
//!   du/dt            = D lap u + (1/2) dp/dt.
//! ```
//!
//! The four-fold anisotropy is carried by
//! `W(theta) = W0 (1 + eps4 cos 4 theta)`, `theta` being the angle of the
//! interface normal `grad p`. With the Karma–Rappel parameter choice
//!
//! ```text
//!   lambda = a1 W0 / d0,        tau0 = a2 lambda W0^2 / D,
//!   a1 = 5 sqrt(2)/8 = 0.883 883,   a2 = 0.6267,
//! ```
//!
//! the interface kinetic coefficient vanishes to leading order in `W0/d0`, so
//! the interface sits on the local Gibbs–Thomson temperature
//! `u_i = -d0(theta) kappa` (+ higher-order terms in `W0`). The `tau ∝ W^2`
//! choice (rather than the mean `tau0`) is what cancels the *anisotropy* of the
//! kinetic coefficient as well; the isotropic statement is the classical
//! Karma–Rappel result.
//!
//! # Discretisation
//!
//! Uniform grid, cell-centred, explicit forward Euler, second-order central
//! differences. The anisotropic Laplacian is expanded as
//!
//! ```text
//!   div(W^2 grad p) = W^2 lap p + 2 W (dW/dtheta) (grad theta . grad p),
//! ```
//!
//! with `theta_x = (p_x p_xy - p_y p_xx)/|grad p|^2` (and the symmetric
//! expression for `theta_y`). Ghost cells impose either symmetry (zero normal
//! derivative) or the far field. Rows are updated in parallel with rayon.
//!
//! Two operators are provided. [`PhaseField::step`] is the plain
//! `div(W^2 grad p)` expansion above; it passes through the right dendrite tip
//! speed and then decays because it omits the variational anisotropy term of
//! KR Eq. (88). [`PhaseField::step_kr`] adds that term,
//! `p_x d_y Q - p_y d_x Q` with `Q = W dW/dtheta`, and uses KR's 9-point
//! isotropic Laplacian; it is the operator whose thin-interface limit selects
//! the KR dendrite, and it matches `solidification-cuda`'s `step_kr_mode` to
//! round-off. `step_kr` falls back to the 5-point, variational-free
//! discretisation at `eps4 = 0`, so the isotropic limit is bit-for-bit the old
//! `step`.
//!
//! Scratch buffers are swapped, so the solver allocates once.

use rayon::prelude::*;

/// First Karma–Rappel constant, `5 sqrt(2)/8`.
pub const A1: f64 = 0.883_883_476_483_184_4;
/// Second Karma–Rappel constant, from the `O(W)` interface-stretching term.
pub const A2: f64 = 0.626_7;

/// Which physical boundary is applied on a pair of opposite edges.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Edge {
    /// Zero normal derivative for `p` and `u` (symmetry plane).
    Symmetry,
    /// `p = -1`, `u = u_far`.
    FarField,
}

/// Boundary condition for the four edges of the rectangle.
#[derive(Clone, Copy, Debug)]
pub struct Bc {
    pub left: Edge,
    pub right: Edge,
    pub bottom: Edge,
    pub top: Edge,
}

impl Bc {
    /// Symmetry in `x = 0, y = 0`, far field on `x = Lx, y = Ly` (dendrite
    /// octant/quadrant).
    pub fn quadrant() -> Self {
        Self {
            left: Edge::Symmetry,
            right: Edge::FarField,
            bottom: Edge::Symmetry,
            top: Edge::FarField,
        }
    }
    /// Far field on every edge (isolated seed in an infinite melt).
    pub fn far_field() -> Self {
        Self {
            left: Edge::FarField,
            right: Edge::FarField,
            bottom: Edge::FarField,
            top: Edge::FarField,
        }
    }
    /// Symmetry on every edge (periodic-like reflecting box).
    pub fn symmetry() -> Self {
        Self {
            left: Edge::Symmetry,
            right: Edge::Symmetry,
            bottom: Edge::Symmetry,
            top: Edge::Symmetry,
        }
    }
}

/// A two-dimensional phase-field state with ghost cells.
///
/// Arrays are laid out row-major with a one-cell ghost border:
/// `idx(i, j) = j * (nx + 2) + i`, `i, j = 0 ..= nx+1 / ny+1`. Interior cells
/// are `i = 1..=nx`, `j = 1..=ny`.
#[derive(Clone)]
pub struct PhaseField {
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub d0: f64,
    pub d: f64,
    pub w0: f64,
    pub eps4: f64,
    pub lambda: f64,
    pub tau0: f64,
    pub bc: Bc,
    pub u_far: f64,
    /// Phase field, `+1` solid, `-1` liquid.
    pub p: Vec<f64>,
    pub u: Vec<f64>,
    p_new: Vec<f64>,
    u_new: Vec<f64>,
    /// `Q = W dW/dtheta` at the interior cells; the ghost cells stay zero and
    /// reproduce `solidification-cuda`'s `pf_q` buffer exactly.
    q: Vec<f64>,
}

impl PhaseField {
    /// Build a solver with `W0 = w0`, capillary length `d0`, diffusivity `d`
    /// and four-fold anisotropy `eps4`. `lambda` and `tau0` follow the
    /// Karma–Rappel kinetic-free choice.
    pub fn new(nx: usize, ny: usize, dx: f64, d0: f64, d: f64, w0: f64, eps4: f64) -> Self {
        let lambda = A1 * w0 / d0;
        let tau0 = A2 * lambda * w0 * w0 / d;
        let n = (nx + 2) * (ny + 2);
        Self {
            nx,
            ny,
            dx,
            d0,
            d,
            w0,
            eps4,
            lambda,
            tau0,
            bc: Bc::quadrant(),
            u_far: 0.0,
            p: vec![-1.0; n],
            u: vec![0.0; n],
            p_new: vec![-1.0; n],
            u_new: vec![0.0; n],
            q: vec![0.0; n],
        }
    }

    /// Override `tau0` (used only to probe kinetic effects; the benchmarks use
    /// the Karma–Rappel value).
    pub fn with_tau0(mut self, tau0: f64) -> Self {
        self.tau0 = tau0;
        self
    }

    #[inline]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        j * (self.nx + 2) + i
    }

    /// Physical `x` of interior column `i` (`i = 1..=nx`).
    #[inline]
    pub fn x(&self, i: usize) -> f64 {
        (i as f64 - 0.5) * self.dx
    }

    #[inline]
    pub fn y(&self, j: usize) -> f64 {
        (j as f64 - 0.5) * self.dx
    }

    /// The simulation box side to the centre of the last cell.
    pub fn length_x(&self) -> f64 {
        self.nx as f64 * self.dx
    }
    pub fn length_y(&self) -> f64 {
        self.ny as f64 * self.dx
    }

    /// A smoothed circular seed of radius `r0` at `(xc, yc)`; temperature set
    /// to `u_far` everywhere.
    pub fn seed_circle(&mut self, xc: f64, yc: f64, r0: f64) {
        for j in 1..=self.ny {
            for i in 1..=self.nx {
                let r = ((self.x(i) - xc).powi(2) + (self.y(j) - yc).powi(2)).sqrt();
                let v = -((r - r0) / (std::f64::consts::SQRT_2 * self.w0)).tanh();
                let k = self.idx(i, j);
                self.p[k] = v;
                self.u[k] = self.u_far;
            }
        }
    }

    /// A quarter-disk seed at the symmetry corner `(0, 0)`.
    pub fn seed_corner(&mut self, r0: f64) {
        self.seed_circle(0.0, 0.0, r0);
    }

    /// A parabolic needle seed for the dendrite: solid where
    /// `x < x_tip - y^2/(2 r_tip)`, smoothed over the interface width.
    /// `x_tip` and `r_tip` are physical lengths; the tip sits on `y = 0`.
    pub fn seed_parabola(&mut self, x_tip: f64, r_tip: f64) {
        for j in 1..=self.ny {
            for i in 1..=self.nx {
                let x = self.x(i);
                let y = self.y(j);
                let xi = x_tip - y * y / (2.0 * r_tip);
                let v = -((x - xi) / (std::f64::consts::SQRT_2 * self.w0)).tanh();
                let k = self.idx(i, j);
                self.p[k] = v;
                self.u[k] = self.u_far;
            }
        }
    }

    /// Largest `x` interior cell with `p >= 0.5` (the actual tip), over all
    /// rows.
    pub fn tip_max_x(&self) -> f64 {
        let mut best = 0.0;
        for j in 1..=self.ny {
            for i in 1..=self.nx {
                if self.p[self.idx(i, j)] >= 0.5 {
                    let x = self.x(i);
                    if x > best {
                        best = x;
                    }
                }
            }
        }
        best
    }

    /// Interpolated `x` where `p` crosses `level` going outward along row `j`
    /// (from `x = 0` toward `+x`), or `None`.
    pub fn crossing_level_x(&self, j: usize, level: f64, xc: f64) -> Option<f64> {
        let mut prev = self.p[self.idx(1, j)];
        for i in 2..=self.nx {
            let cur = self.p[self.idx(i, j)];
            if prev >= level && cur < level {
                let x0 = self.x(i - 1);
                let x1 = self.x(i);
                let t = (prev - level) / (prev - cur);
                return Some(x0 + t * (x1 - x0) - xc);
            }
            prev = cur;
        }
        None
    }

    /// Interpolated `x` where `p` crosses `0.5` going outward along row `j`
    /// (from `x = 0` toward `+x`), or `None`.
    pub fn crossing_x(&self, j: usize, xc: f64) -> Option<f64> {
        self.crossing_level_x(j, 0.5, xc)
    }

    /// Tip radius from a least-squares fit `x(y) = x0 - y^2/(2R)` to the
    /// `p = 0` interface contour.
    ///
    /// The diffuse interface rounds the very tip on the scale `W0`: the local
    /// radius of curvature of the `p = 0` contour there is only a few `W0`,
    /// independent of the macroscopic tip radius. The first ~`1.5 W0` of the
    /// contour is therefore excluded, and the fit runs over the single-valued
    /// monotone needle out to `40 W0` (the paraboloid is monotone in `y`; a
    /// sidebranch would break monotonicity and stop the fit). The `40 W0` cap
    /// is the benchmark scale: the selected dendrite is only approximately
    /// parabolic, and fitting the near-parallel flanks of a very long needle
    /// overshoots `R` by ~25 %. The old estimator fitted the `p = 0.5` contour
    /// over `ny.min(8)` rows, only `~2.8 W0` at `dx = 0.4 W0`, i.e. entirely
    /// inside the rounded cap; that biased `R` low by ~2.4x and made the
    /// Karma-Rappel benchmark look like `sigma* ~ 0.6` instead of the true
    /// ~0.13. See `analysis/solidification.md` §3(ii).
    pub fn tip_radius(&self) -> f64 {
        let y0 = 1.5 * self.w0;
        let y1 = 40.0 * self.w0;
        let mut ys = Vec::new();
        let mut xs = Vec::new();
        let mut prev = f64::INFINITY;
        for j in 1..=self.ny {
            let Some(x) = self.crossing_level_x(j, 0.0, 0.0) else {
                break;
            };
            let y = self.y(j);
            if y < y0 {
                prev = x;
                continue;
            }
            if y > y1 {
                break;
            }
            // The needle contour is monotone decreasing; a sidebranch or the
            // flank turning back stops the single-valued window.
            if x > prev + 0.25 * self.dx {
                break;
            }
            ys.push(y * y);
            xs.push(x);
            prev = x;
        }
        if ys.len() < 4 {
            return f64::NAN;
        }
        // Ordinary least squares x = A + B * (y^2); then R = -1/(2B).
        let n = ys.len() as f64;
        let sy: f64 = ys.iter().sum();
        let sx: f64 = xs.iter().sum();
        let syy: f64 = ys.iter().map(|v| v * v).sum();
        let sxy: f64 = ys.iter().zip(xs.iter()).map(|(a, b)| a * b).sum();
        let den = n * syy - sy * sy;
        if den.abs() < 1e-18 {
            return f64::NAN;
        }
        let b = (n * sxy - sy * sx) / den;
        if b >= -1e-12 {
            return f64::NAN;
        }
        -1.0 / (2.0 * b)
    }

    /// Apply ghost-cell boundary conditions.
    fn apply_bc(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let m = nx + 2;
        let u_far = self.u_far;
        // left/right
        for j in 0..(ny + 2) {
            let base = j * m;
            match self.bc.left {
                Edge::Symmetry => {
                    let pv = self.p[base + 1];
                    let uv = self.u[base + 1];
                    self.p[base] = pv;
                    self.u[base] = uv;
                }
                Edge::FarField => {
                    self.p[base] = -1.0;
                    self.u[base] = u_far;
                }
            }
            match self.bc.right {
                Edge::Symmetry => {
                    let pv = self.p[base + nx];
                    let uv = self.u[base + nx];
                    self.p[base + nx + 1] = pv;
                    self.u[base + nx + 1] = uv;
                }
                Edge::FarField => {
                    self.p[base + nx + 1] = -1.0;
                    self.u[base + nx + 1] = u_far;
                }
            }
        }
        // bottom/top (row 0 is the ghost of row 1; row ny+1 of row ny)
        for i in 0..m {
            match self.bc.bottom {
                Edge::Symmetry => {
                    let pv = self.p[m + i];
                    let uv = self.u[m + i];
                    self.p[i] = pv;
                    self.u[i] = uv;
                }
                Edge::FarField => {
                    self.p[i] = -1.0;
                    self.u[i] = u_far;
                }
            }
            match self.bc.top {
                Edge::Symmetry => {
                    let pv = self.p[ny * m + i];
                    let uv = self.u[ny * m + i];
                    self.p[(ny + 1) * m + i] = pv;
                    self.u[(ny + 1) * m + i] = uv;
                }
                Edge::FarField => {
                    self.p[(ny + 1) * m + i] = -1.0;
                    self.u[(ny + 1) * m + i] = u_far;
                }
            }
        }
    }

    /// Stability limit `dt <= 0.4 min(dx^2 tau0/(4 W0^2), dx^2/(4D))`.
    pub fn stable_dt(&self) -> f64 {
        let dt_p = self.dx * self.dx * self.tau0 / (4.0 * self.w0 * self.w0);
        let dt_u = self.dx * self.dx / (4.0 * self.d);
        // tau is at most tau0 (W/W0)^2; use the smallest.
        0.4 * dt_p.min(dt_u)
    }

    /// Advance one explicit step.
    pub fn step(&mut self, dt: f64) {
        self.apply_bc();
        let nx = self.nx;
        let ny = self.ny;
        let m = nx + 2;
        let dx = self.dx;
        let dx2 = dx * dx;
        let inv_twodx = 1.0 / (2.0 * dx);
        let inv_fourd2 = 1.0 / (4.0 * dx2);
        let w0 = self.w0;
        let eps4 = self.eps4;
        let lambda = self.lambda;
        let tau0 = self.tau0;
        let d = self.d;

        // Split borrows: read `p`/`u`, write the scratch buffers.
        let p = &self.p;
        let u = &self.u;
        let p_new = &mut self.p_new;
        let u_new = &mut self.u_new;

        // Parallel over rows. Each row writes a disjoint stretch of the
        // scratch buffers, so no atomics are needed.
        p_new
            .par_chunks_exact_mut(m)
            .zip(u_new.par_chunks_exact_mut(m))
            .enumerate()
            .for_each(|(j, (pn, un))| {
                if j == 0 || j > ny {
                    return;
                }
                for i in 1..=nx {
                    let k = j * m + i;
                    let p_x = (p[k + 1] - p[k - 1]) * inv_twodx;
                    let p_y = (p[k + m] - p[k - m]) * inv_twodx;
                    let p_xx = (p[k + 1] - 2.0 * p[k] + p[k - 1]) / dx2;
                    let p_yy = (p[k + m] - 2.0 * p[k] + p[k - m]) / dx2;
                    let p_xy =
                        (p[k + m + 1] - p[k + m - 1] - p[k - m + 1] + p[k - m - 1]) * inv_fourd2;
                    let g2 = p_x * p_x + p_y * p_y;
                    let one_minus = 1.0 - p[k] * p[k];
                    let drive = (p[k] - lambda * u[k] * one_minus) * one_minus;
                    let (lap_anis, tau) = if g2 > 1e-18 {
                        let theta = p_y.atan2(p_x);
                        let c4 = (4.0 * theta).cos();
                        let s4 = (4.0 * theta).sin();
                        let w = w0 * (1.0 + eps4 * c4);
                        let dw = -4.0 * w0 * eps4 * s4;
                        let theta_x = (p_x * p_xy - p_y * p_xx) / g2;
                        let theta_y = (p_x * p_yy - p_y * p_xy) / g2;
                        let lap =
                            w * w * (p_xx + p_yy) + 2.0 * w * dw * (theta_x * p_x + theta_y * p_y);
                        let ratio = w / w0;
                        (lap, tau0 * ratio * ratio)
                    } else {
                        (w0 * w0 * (p_xx + p_yy), tau0)
                    };
                    let ds = (lap_anis + drive) / tau;
                    let mut dp = dt * ds;
                    if !dp.is_finite() {
                        dp = 0.0;
                    }
                    pn[i] = p[k] + dp;
                    let lap_u = (u[k + 1] - 2.0 * u[k] + u[k - 1]) / dx2
                        + (u[k + m] - 2.0 * u[k] + u[k - m]) / dx2;
                    un[i] = u[k] + d * lap_u * dt + 0.5 * dp;
                }
            });
        // Copy the interior scratch into the state; ghost cells of `p_new`
        // are never read, and `apply_bc` refreshes the state's ghosts next
        // step.
        for j in 1..=ny {
            let s = j * m + 1;
            let e = j * m + nx + 1;
            self.p[s..e].copy_from_slice(&self.p_new[s..e]);
            self.u[s..e].copy_from_slice(&self.u_new[s..e]);
        }
    }

    /// Advance one explicit step with the **corrected** Karma–Rappel
    /// operator (KR 1998 Eq. 88): [`Self::step`] plus the variational
    /// anisotropy term, optionally with KR's 9-point isotropic Laplacian.
    /// With `Q = W dW/dtheta`, the variational term is `p_x d_y Q - p_y d_x Q`;
    /// it is what makes the interface normal relax toward the anisotropic
    /// equilibrium shape and gives a selected dendrite (the plain
    /// [`Self::step`] passes through the right tip speed and then decays).
    ///
    /// `mode` mirrors `solidification-cuda`'s `step_kr_mode` bit for bit so
    /// the two agree to round-off (`Q` is built from the algebraic `c4`/`s4`
    /// rather than `atan2`, exactly as in the device kernel):
    /// * bit 0 — use KR's 9-point isotropic Laplacian (otherwise 5-point);
    /// * bit 1 — include the variational anisotropy term.
    ///
    /// [`Self::step_kr`] is `mode = 3` (the GPU production operator); at
    /// `eps4 = 0` it falls back to `mode = 2` so the isotropic results stay
    /// bit-for-bit [`Self::step`]. The GPU equality check lives in
    /// `solidification-cuda`.
    pub fn step_kr_mode(&mut self, dt: f64, mode: i32) {
        self.apply_bc();
        let nx = self.nx;
        let ny = self.ny;
        let m = nx + 2;
        let dx = self.dx;
        let dx2 = dx * dx;
        let inv_twodx = 1.0 / (2.0 * dx);
        let inv_fourd2 = 1.0 / (4.0 * dx2);
        let w0 = self.w0;
        let eps4 = self.eps4;
        let lambda = self.lambda;
        let tau0 = self.tau0;
        let d = self.d;

        // Q = W dW/dtheta on the interior (identical expression to the GPU
        // `pf_q` kernel: c4 and s4 are evaluated algebraically rather than via
        // atan2, so the two agree to round-off). Ghost q are never written, so
        // they stay zero, as in the GPU buffer.
        {
            let p = &self.p;
            let q = &mut self.q;
            q.par_chunks_exact_mut(m).enumerate().for_each(|(j, qrow)| {
                if j == 0 || j > ny {
                    return;
                }
                for (i, qv) in qrow.iter_mut().enumerate().take(nx + 1).skip(1) {
                    let k = j * m + i;
                    let px = (p[k + 1] - p[k - 1]) * inv_twodx;
                    let py = (p[k + m] - p[k - m]) * inv_twodx;
                    let g2 = px * px + py * py;
                    *qv = if g2 > 1e-18 {
                        let g2i = 1.0 / g2;
                        let c4 = (px * px * px * px - 6.0 * px * px * py * py
                            + py * py * py * py)
                            * g2i
                            * g2i;
                        let s4 = 4.0 * px * py * (px * px - py * py) * g2i * g2i;
                        let w = w0 * (1.0 + eps4 * c4);
                        let dw = -4.0 * w0 * eps4 * s4;
                        w * dw
                    } else {
                        0.0
                    };
                }
            });
        }

        let p = &self.p;
        let u = &self.u;
        let q = &self.q;
        let p_new = &mut self.p_new;
        let u_new = &mut self.u_new;

        p_new
            .par_chunks_exact_mut(m)
            .zip(u_new.par_chunks_exact_mut(m))
            .enumerate()
            .for_each(|(j, (pn, un))| {
                if j == 0 || j > ny {
                    return;
                }
                for i in 1..=nx {
                    let k = j * m + i;
                    let p_x = (p[k + 1] - p[k - 1]) * inv_twodx;
                    let p_y = (p[k + m] - p[k - m]) * inv_twodx;
                    let p_xx = (p[k + 1] - 2.0 * p[k] + p[k - 1]) / dx2;
                    let p_yy = (p[k + m] - 2.0 * p[k] + p[k - m]) / dx2;
                    let p_xy =
                        (p[k + m + 1] - p[k + m - 1] - p[k - m + 1] + p[k - m - 1]) * inv_fourd2;
                    let g2 = p_x * p_x + p_y * p_y;
                    // Isotropic 9-point Laplacian (KR 1998 Sec. III), which
                    // suppresses the grid anisotropy of the 5-point stencil.
                    // The reciprocal multiply matches the device kernel's
                    // `inv_sixdx2` exactly.
                    let inv_sixdx2 = 1.0 / (6.0 * dx2);
                    let lap9 = (4.0
                        * (p[k + 1] + p[k - 1] + p[k + m] + p[k - m])
                        + (p[k + m + 1] + p[k + m - 1] + p[k - m + 1] + p[k - m - 1])
                        - 20.0 * p[k])
                        * inv_sixdx2;
                    let lap = if mode & 1 != 0 { lap9 } else { p_xx + p_yy };
                    let one_minus = 1.0 - p[k] * p[k];
                    let drive = (p[k] - lambda * u[k] * one_minus) * one_minus;
                    let (lap_anis, tau) = if g2 > 1e-18 {
                        let g2i = 1.0 / g2;
                        let c4 = (p_x * p_x * p_x * p_x - 6.0 * p_x * p_x * p_y * p_y
                            + p_y * p_y * p_y * p_y)
                            * g2i
                            * g2i;
                        let s4 = 4.0 * p_x * p_y * (p_x * p_x - p_y * p_y) * g2i * g2i;
                        let w = w0 * (1.0 + eps4 * c4);
                        let dw = -4.0 * w0 * eps4 * s4;
                        let theta_x = (p_x * p_xy - p_y * p_xx) / g2;
                        let theta_y = (p_x * p_yy - p_y * p_xy) / g2;
                        let base =
                            w * w * lap + 2.0 * w * dw * (theta_x * p_x + theta_y * p_y);
                        // Variational anisotropy: p_x d_y Q - p_y d_x Q.
                        let mut anis = base;
                        if mode & 2 != 0 {
                            let q_y = (q[k + m] - q[k - m]) * inv_twodx;
                            let q_x = (q[k + 1] - q[k - 1]) * inv_twodx;
                            anis += p_x * q_y - p_y * q_x;
                        }
                        let ratio = w / w0;
                        (anis, tau0 * ratio * ratio)
                    } else {
                        (w0 * w0 * lap, tau0)
                    };
                    let ds = (lap_anis + drive) / tau;
                    let mut dp = dt * ds;
                    if !dp.is_finite() {
                        dp = 0.0;
                    }
                    pn[i] = p[k] + dp;
                    let lap_u = (u[k + 1] - 2.0 * u[k] + u[k - 1]) / dx2
                        + (u[k + m] - 2.0 * u[k] + u[k - m]) / dx2;
                    un[i] = u[k] + d * lap_u * dt + 0.5 * dp;
                }
            });
        for j in 1..=ny {
            let s = j * m + 1;
            let e = j * m + nx + 1;
            self.p[s..e].copy_from_slice(&self.p_new[s..e]);
            self.u[s..e].copy_from_slice(&self.u_new[s..e]);
        }
    }

    /// The corrected Karma–Rappel operator used by the benchmarks: mode 3
    /// (9-point Laplacian plus the variational term), falling back to mode 2
    /// at `eps4 = 0` so the isotropic results stay exactly [`Self::step`].
    pub fn step_kr(&mut self, dt: f64) {
        let mode = if self.eps4 == 0.0 { 2 } else { 3 };
        self.step_kr_mode(dt, mode);
    }

    /// Total solid fraction, used for conservation checks.
    pub fn solid_fraction(&self) -> f64 {
        let n = (self.nx * self.ny) as f64;
        let s: f64 = (1..=self.ny)
            .flat_map(|j| (1..=self.nx).map(move |i| self.idx(i, j)))
            .map(|k| 0.5 * (1.0 + self.p[k]))
            .sum();
        s / n
    }

    /// Mean dimensionless temperature over the interior.
    pub fn mean_u(&self) -> f64 {
        let n = (self.nx * self.ny) as f64;
        let s: f64 = (1..=self.ny)
            .flat_map(|j| (1..=self.nx).map(move |i| self.idx(i, j)))
            .map(|k| self.u[k])
            .sum();
        s / n
    }

    /// The rightmost interface crossing `p = 0` along the symmetry row `j0`
    /// (for a dendrite growing in `+x`). Returns the physical `x`, or `None`.
    pub fn tip_x(&self, j0: usize) -> Option<f64> {
        let mut last = None;
        for i in 1..=self.nx {
            let k = self.idx(i, j0);
            if self.p[k] >= 0.0 {
                last = Some(self.x(i));
            }
        }
        last
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seed_is_solid_and_far_field_liquid() {
        let mut pf = PhaseField::new(32, 32, 0.25, 1.0, 2.0, 0.5, 0.05);
        pf.u_far = -0.55;
        pf.seed_circle(4.0, 4.0, 1.5);
        assert!(pf.p[pf.idx(17, 17)] > 0.9);
        assert!(pf.p[pf.idx(1, 1)] < 0.0);
    }

    /// With `eps4 = 0` the model is isotropic: a circular seed must stay
    /// circular (the field is a function of radius only).
    #[test]
    fn isotropic_seed_stays_radial() {
        let mut pf = PhaseField::new(48, 48, 0.25, 1.0, 2.0, 0.5, 0.0);
        pf.bc = Bc::far_field();
        pf.u_far = -0.3;
        pf.seed_circle(6.0, 6.0, 2.0);
        let dt = pf.stable_dt();
        for _ in 0..20 {
            pf.step(dt);
        }
        // Cells at the same radius on the x and y axes differ only by
        // discretisation; check two symmetry-equivalent cells.
        let a = pf.p[pf.idx(24, 24)];
        let b = pf.p[pf.idx(24, 24)];
        assert_eq!(a, b);
        // Compare a ring sample: (18,24) and (24,18) are mirror images.
        let c = pf.p[pf.idx(18, 24)];
        let d = pf.p[pf.idx(24, 18)];
        assert!((c - d).abs() < 1e-10, "isotropy broken: {c} vs {d}");
    }

    #[test]
    fn a_flat_interface_stays_flat() {
        // A planar interface in 1D (y-independent) must not develop structure.
        let mut pf = PhaseField::new(8, 32, 0.5, 1.0, 2.0, 0.5, 0.05);
        pf.bc = Bc {
            left: Edge::FarField,
            right: Edge::FarField,
            bottom: Edge::Symmetry,
            top: Edge::Symmetry,
        };
        pf.u_far = -0.55;
        for j in 1..=pf.ny {
            for i in 1..=pf.nx {
                let k = pf.idx(i, j);
                if i <= 2 {
                    pf.p[k] = 1.0;
                }
                pf.u[k] = -0.55;
            }
        }
        let dt = pf.stable_dt();
        for _ in 0..10 {
            pf.step(dt);
        }
        for i in 1..=pf.nx {
            let a = pf.p[pf.idx(i, 3)];
            let b = pf.p[pf.idx(i, 25)];
            assert!((a - b).abs() < 1e-12, "row divergence at i={i}");
        }
    }

    /// Independent reference for one cell of the corrected operator: the
    /// interface-normal angle is taken with `atan2` and `Q = W dW/dtheta` is
    /// rebuilt from trig, a different numerical route from the production
    /// algebraic `c4`/`s4`. A sign or factor error in the variational term is
    /// therefore caught.
    fn hand_update(pf: &PhaseField, dt: f64, i: usize, j: usize) -> (f64, f64) {
        let m = pf.nx + 2;
        let dx = pf.dx;
        let inv2 = 1.0 / (2.0 * dx);
        let inv4 = 1.0 / (4.0 * dx * dx);
        let p = &pf.p;
        let u = &pf.u;
        let q_at = |ci: usize, cj: usize| -> f64 {
            let k = cj * m + ci;
            let px = (p[k + 1] - p[k - 1]) * inv2;
            let py = (p[k + m] - p[k - m]) * inv2;
            let g2 = px * px + py * py;
            if g2 <= 1e-18 {
                return 0.0;
            }
            let theta = py.atan2(px);
            let w = pf.w0 * (1.0 + pf.eps4 * (4.0 * theta).cos());
            let dw = -4.0 * pf.w0 * pf.eps4 * (4.0 * theta).sin();
            w * dw
        };
        let k = j * m + i;
        let px = (p[k + 1] - p[k - 1]) * inv2;
        let py = (p[k + m] - p[k - m]) * inv2;
        let pxx = (p[k + 1] - 2.0 * p[k] + p[k - 1]) / (dx * dx);
        let pyy = (p[k + m] - 2.0 * p[k] + p[k - m]) / (dx * dx);
        let pxy = (p[k + m + 1] - p[k + m - 1] - p[k - m + 1] + p[k - m - 1]) * inv4;
        let g2 = px * px + py * py;
        let qx = (q_at(i + 1, j) - q_at(i - 1, j)) * inv2;
        let qy = (q_at(i, j + 1) - q_at(i, j - 1)) * inv2;
        let (lap_anis, tau) = if g2 > 1e-18 {
            let theta = py.atan2(px);
            let w = pf.w0 * (1.0 + pf.eps4 * (4.0 * theta).cos());
            let dw = -4.0 * pf.w0 * pf.eps4 * (4.0 * theta).sin();
            let tx = (px * pxy - py * pxx) / g2;
            let ty = (px * pyy - py * pxy) / g2;
            // Production `step_kr` is mode 3: 9-point Laplacian plus the
            // variational term.
            let lap9 = (4.0 * (p[k + 1] + p[k - 1] + p[k + m] + p[k - m])
                + (p[k + m + 1] + p[k + m - 1] + p[k - m + 1] + p[k - m - 1])
                - 20.0 * p[k])
                / (6.0 * dx * dx);
            let base = w * w * lap9 + 2.0 * w * dw * (tx * px + ty * py);
            (base + px * qy - py * qx, pf.tau0 * (w / pf.w0).powi(2))
        } else {
            (pf.w0 * pf.w0 * (pxx + pyy), pf.tau0)
        };
        let om = 1.0 - p[k] * p[k];
        let drive = (p[k] - pf.lambda * u[k] * om) * om;
        let dp = dt * ((lap_anis + drive) / tau);
        let lap_u = (u[k + 1] - 2.0 * u[k] + u[k - 1]) / (dx * dx)
            + (u[k + m] - 2.0 * u[k] + u[k - m]) / (dx * dx);
        (p[k] + dp, u[k] + pf.d * lap_u * dt + 0.5 * dp)
    }

    /// The corrected step must reproduce a hand-computed update at a cell on a
    /// curved interface (where the variational term is live).
    #[test]
    fn kr_step_matches_hand_computed_cell() {
        let mut pf = PhaseField::new(5, 5, 0.5, 1.0, 2.0, 0.5, 0.05);
        pf.bc = Bc::far_field();
        pf.u_far = -0.3;
        // A curved interface passing through the centre cell (3,3).
        pf.seed_circle(0.0, 0.0, 1.6);
        // Vary u so the drive and the latent-heat update are exercised too.
        for j in 1..=pf.ny {
            for i in 1..=pf.nx {
                let k = pf.idx(i, j);
                pf.u[k] = -0.3 + 0.01 * i as f64 - 0.02 * j as f64;
            }
        }
        pf.apply_bc();
        let dt = 1e-3;
        let (pe, ue) = hand_update(&pf, dt, 3, 3);
        pf.step_kr(dt);
        let k = pf.idx(3, 3);
        assert!(
            (pf.p[k] - pe).abs() < 1e-12,
            "p mismatch: {} vs hand {} ",
            pf.p[k],
            pe
        );
        assert!(
            (pf.u[k] - ue).abs() < 1e-12,
            "u mismatch: {} vs hand {} ",
            pf.u[k],
            ue
        );
    }

    /// `eps4 = 0` is the isotropic limit: the variational term is exactly zero
    /// (`Q ∝ eps4`), so the corrected step must be bit-for-bit the old step.
    #[test]
    fn kr_step_is_exactly_isotropic_at_eps4_zero() {
        let mut a = PhaseField::new(32, 32, 0.25, 1.0, 2.0, 0.5, 0.0);
        a.bc = Bc::far_field();
        a.u_far = -0.3;
        a.seed_circle(4.0, 4.0, 1.5);
        let mut b = a.clone();
        let dt = a.stable_dt();
        for _ in 0..20 {
            a.step(dt);
            b.step_kr(dt);
        }
        for k in 0..a.p.len() {
            assert_eq!(a.p[k], b.p[k], "p diverged at {k}");
            assert_eq!(a.u[k], b.u[k], "u diverged at {k}");
        }
    }

    /// The variational term is active for `eps4 != 0` and changes the update.
    #[test]
    fn kr_step_differs_from_plain_when_anisotropic() {
        let mut a = PhaseField::new(48, 48, 0.4, 1.0, 2.0, 0.5, 0.05);
        a.bc = Bc::quadrant();
        a.u_far = -0.55;
        a.seed_parabola(8.0, 4.0);
        let mut b = a.clone();
        let dt = a.stable_dt();
        for _ in 0..20 {
            a.step(dt);
            b.step_kr(dt);
        }
        let dmax = (1..=a.ny)
            .flat_map(|j| (1..=a.nx).map(move |i| (i, j)))
            .map(|(i, j)| (a.p[a.idx(i, j)] - b.p[b.idx(i, j)]).abs())
            .fold(0.0f64, f64::max);
        assert!(dmax > 1e-8, "variational term inactive, max|dp|={dmax:.3e}");
    }
}
