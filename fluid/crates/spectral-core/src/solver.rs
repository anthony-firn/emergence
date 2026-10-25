//! Incompressible Navier–Stokes on the periodic box, pseudo-spectral.
//!
//!   ∂t û = P[ (u × ω)^ ] - ν|k|² û + P f̂,        P = I - k kᵀ/|k|²
//!
//! Rotational form (the gradient part of (u·∇)u is absorbed into pressure,
//! which P removes), 2/3-rule dealiasing, and the viscous term integrated
//! exactly by an integrating factor inside classical RK4 (Lawson's IF-RK4).
//! The semi-discrete Galerkin system conserves energy exactly when ν = 0.

use num_complex::Complex64;
use rayon::prelude::*;

use crate::fft::Fft3;
use crate::grid::Grid;
use crate::penalization::Penalization;

pub type Spec = Vec<Complex64>;

const I: Complex64 = Complex64::new(0.0, 1.0);

/// Smooth body force, given in spectral space (already divergence-free or
/// projected on use). `None` means unforced.
pub type Forcing = Option<[Spec; 3]>;

pub struct Solver {
    pub grid: Grid,
    pub nu: f64,
    pub t: f64,
    pub u: [Spec; 3],
    pub forcing: Forcing,
    /// Solid objects / sponge (Brinkman penalization), applied in physical space.
    pub penalization: Option<Penalization>,
    /// Surface-constrained walls (markers on the true geometry), solved in
    /// the IMEX implicit stages; takes the place of solid penalization.
    pub surface: Option<crate::surface::Surface>,
    pub(crate) imex: Option<Box<crate::imex::ImexWork>>,
    /// Uniform body force (a mean pressure gradient −∇P̄, per unit mass).
    /// `None`: the mean flow is held fixed. `Some(G)`: the mean is free and
    /// driven by G; at steady state the force on any solid is exactly G·(2π)³.
    pub body_force: Option<[f64; 3]>,
    /// Conducting fluid in a magnetic field (quasi-static, low Rm).
    pub mhd: Option<crate::mhd::Mhd>,
    rng: u64,
    pub(crate) fft: Fft3,
    pub(crate) k2: Vec<f64>,
    // workspaces
    pub(crate) phys_u: [Vec<f64>; 3],
    pub(crate) phys_w: [Vec<f64>; 3],
    pub(crate) phys_n: [Vec<f64>; 3],
    pub(crate) tmp: Spec,
    stage: [Spec; 3],
    k_a: [Spec; 3],
    k_b: [Spec; 3],
    k_c: [Spec; 3],
    k_d: [Spec; 3],
}

/// Zero dealiased modes and remove the gradient part of a spectral vector
/// field: v ← mask·(I - k kᵀ/|k|²) v. The k = 0 (mean flow) mode is zeroed.
pub fn project(g: &Grid, v: &mut [Spec; 3]) {
    project_with(g, v, false)
}

/// As [`project`], but with `keep_mean` the k = 0 mode is left untouched
/// (used when the mean flow is free and driven by a uniform body force).
pub fn project_with(g: &Grid, v: &mut [Spec; 3], keep_mean: bool) {
    let saved = [v[0][0], v[1][0], v[2][0]];
    project_inner(g, v);
    if keep_mean {
        for c in 0..3 {
            v[c][0] = saved[c];
        }
    }
}

fn project_inner(g: &Grid, v: &mut [Spec; 3]) {
    let [a, b, c] = v;
    a.par_iter_mut()
        .zip(b.par_iter_mut())
        .zip(c.par_iter_mut())
        .enumerate()
        .for_each(|(idx, ((x, y), z))| {
            let k = g.k(idx);
            let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
            if !g.kept(idx) || k2 == 0.0 {
                *x = Complex64::default();
                *y = Complex64::default();
                *z = Complex64::default();
                return;
            }
            let kdot = (*x * k[0] + *y * k[1] + *z * k[2]) / k2;
            *x -= kdot * k[0];
            *y -= kdot * k[1];
            *z -= kdot * k[2];
        });
}

/// Spectrum of a physical velocity field, dealiased and made divergence-free.
pub fn initial_spectrum(g: &Grid, u_phys: &[Vec<f64>; 3]) -> [Spec; 3] {
    let fft = Fft3::new(g);
    let mut u = zeros_spec(g.spec_len());
    for c in 0..3 {
        fft.forward(&u_phys[c], &mut u[c]);
    }
    project(g, &mut u);
    u
}

pub(crate) fn zeros_spec(n: usize) -> [Spec; 3] {
    std::array::from_fn(|_| vec![Complex64::default(); n])
}

impl Solver {
    /// Build from a physical velocity field (each component N³ reals).
    /// The field is dealiased and projected to be exactly divergence-free.
    pub fn new(grid: Grid, nu: f64, u_phys: [Vec<f64>; 3]) -> Self {
        let fft = Fft3::new(&grid);
        let ns = grid.spec_len();
        let nr = grid.real_len();
        let k2 = (0..ns)
            .map(|idx| {
                let k = grid.k(idx);
                k[0] * k[0] + k[1] * k[1] + k[2] * k[2]
            })
            .collect();
        let mut s = Self {
            nu,
            t: 0.0,
            u: zeros_spec(ns),
            forcing: None,
            penalization: None,
            surface: None,
            imex: None,
            body_force: None,
            mhd: None,
            rng: 0x9E37_79B9_7F4A_7C15,
            fft,
            k2,
            phys_u: std::array::from_fn(|_| vec![0.0; nr]),
            phys_w: std::array::from_fn(|_| vec![0.0; nr]),
            phys_n: std::array::from_fn(|_| vec![0.0; nr]),
            tmp: vec![Complex64::default(); ns],
            stage: zeros_spec(ns),
            k_a: zeros_spec(ns),
            k_b: zeros_spec(ns),
            k_c: zeros_spec(ns),
            k_d: zeros_spec(ns),
            grid,
        };
        for c in 0..3 {
            s.fft.forward(&u_phys[c], &mut s.u[c]);
        }
        let mut u = std::mem::take(&mut s.u);
        s.dealias_and_project(&mut u);
        s.u = u;
        s
    }

    /// Zero dealiased modes and remove the gradient (non-solenoidal) part.
    pub fn dealias_and_project(&self, v: &mut [Spec; 3]) {
        project(&self.grid, v);
    }

    /// Nonlinear (+ forcing) tendency N(û) = P[(u × ω)^] + P f̂, written into `out`.
    pub(crate) fn rhs(&mut self, uh: &[Spec; 3], out: &mut [Spec; 3]) {
        self.rhs_shifted(uh, out, None)
    }

    /// Per-axis phase factors e^{i k Δ} for a grid shift Δ.
    fn shift_phases(&self, d: [f64; 3]) -> [Vec<Complex64>; 3] {
        let g = &self.grid;
        [
            g.kx.iter().map(|&k| Complex64::from_polar(1.0, k * d[0])).collect(),
            g.kx.iter().map(|&k| Complex64::from_polar(1.0, k * d[1])).collect(),
            g.kz.iter().map(|&k| Complex64::from_polar(1.0, k * d[2])).collect(),
        ]
    }

    /// Tendency with the nonlinear term evaluated on a grid shifted by Δ
    /// (phase-shift dealiasing); `None` means no shift.
    pub(crate) fn rhs_shifted(&mut self, uh: &[Spec; 3], out: &mut [Spec; 3], shift: Option<[f64; 3]>) {
        let ph = shift.map(|d| self.shift_phases(d));
        let g = &self.grid;
        let phase = |idx: usize| -> Complex64 {
            match &ph {
                None => Complex64::new(1.0, 0.0),
                Some([ex, ey, ez]) => {
                    let (i, j, l) = g.spec_ijl(idx);
                    ex[i] * ey[j] * ez[l]
                }
            }
        };
        // velocity to physical space (on the shifted grid: û e^{ik·Δ})
        for c in 0..3 {
            let uc = &uh[c];
            self.tmp.par_iter_mut().enumerate().for_each(|(idx, t)| *t = uc[idx] * phase(idx));
            self.fft.inverse_inplace(&mut self.tmp, &mut self.phys_u[c]);
        }
        // vorticity ω̂ = i k × û, to physical space
        for c in 0..3 {
            let (a, b) = ((c + 1) % 3, (c + 2) % 3);
            let (ua, ub) = (&uh[a], &uh[b]);
            self.tmp.par_iter_mut().enumerate().for_each(|(idx, w)| {
                let k = g.k(idx);
                *w = I * (ub[idx] * k[a] - ua[idx] * k[b]) * phase(idx);
            });
            self.fft.inverse_inplace(&mut self.tmp, &mut self.phys_w[c]);
        }
        // n = u × ω in physical space
        let [u0, u1, u2] = &self.phys_u;
        let [w0, w1, w2] = &self.phys_w;
        let [n0, n1, n2] = &mut self.phys_n;
        n0.par_iter_mut()
            .zip(n1.par_iter_mut())
            .zip(n2.par_iter_mut())
            .enumerate()
            .for_each(|(p, ((a, b), c))| {
                *a = u1[p] * w2[p] - u2[p] * w1[p];
                *b = u2[p] * w0[p] - u0[p] * w2[p];
                *c = u0[p] * w1[p] - u1[p] * w0[p];
            });
        if let Some(pen) = &self.penalization {
            pen.apply(&self.phys_u, &mut self.phys_n);
        }
        for c in 0..3 {
            self.fft.forward(&self.phys_n[c], &mut out[c]);
            if ph.is_some() {
                out[c].par_iter_mut().enumerate().for_each(|(idx, v)| *v *= phase(idx).conj());
            }
        }
        if let Some(f) = &self.forcing {
            for c in 0..3 {
                out[c].par_iter_mut().zip(f[c].par_iter()).for_each(|(o, fv)| *o += *fv);
            }
        }
        if let Some(m) = &self.mhd {
            m.add_force(&self.grid, &self.fft, uh, out);
        }
        match self.body_force {
            None => self.dealias_and_project(out),
            Some(bf) => {
                // free mean: keep the k = 0 tendency (e.g. explicit drag) and add G
                project_with(&self.grid, out, true);
                for c in 0..3 {
                    out[c][0] += Complex64::new(bf[c], 0.0);
                }
            }
        }
    }

    /// Integrating factor e^{-ν|k|² h} applied in place: v ← E(h) v.
    fn apply_if(&self, h: f64, v: &mut [Spec; 3]) {
        let nu = self.nu;
        if let Some(m) = self.mhd.as_ref().filter(|m| m.joule_in_if && m.b_var.is_none() && m.fluid_mask.is_none()) {
            let g = &self.grid;
            for comp in v.iter_mut() {
                comp.par_iter_mut()
                    .zip(self.k2.par_iter())
                    .enumerate()
                    .for_each(|(i, (x, &k2))| *x *= (-(nu * k2 + m.joule_rate(g.k(i))) * h).exp());
            }
            return;
        }
        for comp in v.iter_mut() {
            comp.par_iter_mut()
                .zip(self.k2.par_iter())
                .for_each(|(x, &k2)| *x *= (-nu * k2 * h).exp());
        }
    }

    /// out ← a + s·b (componentwise).
    pub(crate) fn axpy(out: &mut [Spec; 3], a: &[Spec; 3], s: f64, b: &[Spec; 3]) {
        for c in 0..3 {
            out[c]
                .par_iter_mut()
                .zip(a[c].par_iter().zip(b[c].par_iter()))
                .for_each(|(o, (x, y))| *o = *x + *y * s);
        }
    }

    /// Advance one step of size h: IF-RK4, or the IMEX scheme when a solid is
    /// treated implicitly.
    pub fn step(&mut self, h: f64) {
        if self.surface.is_some() || self.penalization.as_ref().is_some_and(|p| p.implicit) {
            assert!(!self.mhd.as_ref().is_some_and(|m| m.joule_in_if), "joule_in_if needs the IF-RK4 path");
            self.step_imex(h)
        } else {
            self.step_ifrk4(h)
        }
    }

    /// Final relative CG residuals of the last IMEX step's implicit stages.
    pub fn last_cg_resid(&self) -> Option<[f64; 3]> {
        self.imex.as_ref().map(|w| [w.last_resid[1], w.last_resid[2], w.last_resid[3]])
    }

    /// CG iterations of the last IMEX step's implicit stages (2..4).
    pub fn last_cg_iters(&self) -> Option<[usize; 3]> {
        self.imex.as_ref().map(|w| [w.last_iters[1], w.last_iters[2], w.last_iters[3]])
    }

    /// One Lawson IF-RK4 step of size h.
    pub fn step_ifrk4(&mut self, h: f64) {
        let u0 = std::mem::take(&mut self.u);
        let mut st = std::mem::take(&mut self.stage);
        let (mut ka, mut kb, mut kc, mut kd) = (
            std::mem::take(&mut self.k_a),
            std::mem::take(&mut self.k_b),
            std::mem::take(&mut self.k_c),
            std::mem::take(&mut self.k_d),
        );

        // Phase-shift dealiasing on spherical grids: stages a, c use a random
        // shift Δ and b, d use Δ + half a cell, so the shifted groups carry
        // equal weight (1/6 + 1/3 each) and single/triple aliases cancel.
        let shifts: [Option<[f64; 3]>; 4] = if self.grid.spherical {
            let d = self.next_shift();
            let hc = 0.5 * self.grid.dx();
            let d2 = [d[0] + hc, d[1] + hc, d[2] + hc];
            [Some(d), Some(d2), Some(d), Some(d2)]
        } else {
            [None; 4]
        };
        // a = N(u0);  stage1 = E(h/2)(u0 + h/2 a)
        self.rhs_shifted(&u0, &mut ka, shifts[0]);
        Self::axpy(&mut st, &u0, 0.5 * h, &ka);
        self.apply_if(0.5 * h, &mut st);
        // b = N(stage1);  stage2 = E(h/2) u0 + h/2 b
        self.rhs_shifted(&st, &mut kb, shifts[1]);
        let mut e_half_u0 = u0.clone();
        self.apply_if(0.5 * h, &mut e_half_u0);
        Self::axpy(&mut st, &e_half_u0, 0.5 * h, &kb);
        // c = N(stage2);  stage3 = E(h) u0 + h E(h/2) c  = E(h/2)(E(h/2) u0 + h c)
        self.rhs_shifted(&st, &mut kc, shifts[2]);
        Self::axpy(&mut st, &e_half_u0, h, &kc);
        self.apply_if(0.5 * h, &mut st);
        // d = N(stage3)
        self.rhs_shifted(&st, &mut kd, shifts[3]);
        // u1 = E(h)(u0 + h/6 a) + E(h/2) h/3 (b + c) + h/6 d
        //    = E(h/2)[ E(h/2)(u0 + h/6 a) + h/3 (b + c) ] + h/6 d
        let mut u1 = u0;
        Self::axpy(&mut st, &u1, h / 6.0, &ka); // st = u0 + h/6 a
        self.apply_if(0.5 * h, &mut st);
        Self::axpy(&mut u1, &st, h / 3.0, &kb);
        Self::axpy(&mut st, &u1, h / 3.0, &kc);
        self.apply_if(0.5 * h, &mut st);
        Self::axpy(&mut u1, &st, h / 6.0, &kd);

        self.u = u1;
        self.stage = st;
        self.k_a = ka;
        self.k_b = kb;
        self.k_c = kc;
        self.k_d = kd;
        self.t += h;
    }

    /// Physical velocity component c (fresh array).
    pub fn velocity_phys(&mut self, c: usize) -> Vec<f64> {
        let mut out = vec![0.0; self.grid.real_len()];
        self.fft.inverse(&self.u[c], &mut self.tmp, &mut out);
        out
    }

    /// Next random grid shift in [0, dx)³ (deterministic xorshift sequence).
    fn next_shift(&mut self) -> [f64; 3] {
        let dx = self.grid.dx();
        let mut r = || {
            self.rng ^= self.rng << 13;
            self.rng ^= self.rng >> 7;
            self.rng ^= self.rng << 17;
            (self.rng >> 11) as f64 / (1u64 << 53) as f64 * dx
        };
        [r(), r(), r()]
    }

    /// Set the uniform (k = 0) velocity. The right-hand side never changes
    /// this mode, so it acts as a fixed free stream (the implied uniform body
    /// force balances the drag of any penalized objects).
    pub fn set_mean_flow(&mut self, m: [f64; 3]) {
        for c in 0..3 {
            self.u[c][0] = Complex64::new(m[c], 0.0);
        }
    }

    /// The uniform (k = 0) velocity: the mean flow over the box.
    pub fn mean_velocity(&self) -> [f64; 3] {
        [self.u[0][0].re, self.u[1][0].re, self.u[2][0].re]
    }

    /// Force of the fluid on the penalized solid, if any (unit density).
    pub fn solid_force(&mut self) -> Option<[f64; 3]> {
        if self.surface.is_some() {
            return self.imex.as_ref().map(|w| w.last_force);
        }
        if self.penalization.as_ref()?.implicit {
            return self.imex.as_ref().map(|w| w.last_force);
        }
        let u: [Vec<f64>; 3] = std::array::from_fn(|c| self.velocity_phys(c));
        self.penalization.as_ref().map(|p| p.solid_force(&self.grid, &u))
    }

    pub fn k2(&self) -> &[f64] {
        &self.k2
    }
}
