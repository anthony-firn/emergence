//! GPU (CUDA, NVRTC) port of the Karma-Rappel anisotropic phase-field
//! solidification step from `fluid/crates/solidification/src/phase_field.rs`.
//!
//! The device kernel is a line-by-line transcription of the CPU update, so
//! `f64` on the GPU reproduces the CPU solver to round-off (NVRTC is compiled
//! with `--fmad=false` and IEEE division/sqrt to avoid FMA contraction, which
//! Rust does not perform at the default `x86-64` target). `f32` is the fast
//! consumer-GPU path used for the long convergence runs and is compared
//! against the f64 result.
//!
//! The state is cell-centred with a one-cell ghost border, exactly as the CPU
//! solver; boundary conditions are the four edges `Symmetry` / `FarField`.

use cudarc::driver::{
    CudaContext, CudaFunction, CudaSlice, CudaStream, DeviceRepr, LaunchConfig, PushKernelArg,
    ValidAsZeroBits,
};
use cudarc::nvrtc::{compile_ptx_with_opts, CompileOptions};
use solidification::phase_field::{Bc, Edge};
use std::sync::Arc;

pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

const BLOCK: u32 = 256;

/// Real type on the device.
pub trait Real: DeviceRepr + ValidAsZeroBits + Copy + Default + Send + Sync + 'static {
    const SUFFIX: &'static str;
    fn from_f64(v: f64) -> Self;
    fn to_f64(self) -> f64;
}
impl Real for f64 {
    const SUFFIX: &'static str = "f64";
    fn from_f64(v: f64) -> Self {
        v
    }
    fn to_f64(self) -> f64 {
        self
    }
}
impl Real for f32 {
    const SUFFIX: &'static str = "f32";
    fn from_f64(v: f64) -> Self {
        v as f32
    }
    fn to_f64(self) -> f64 {
        self as f64
    }
}

struct Kern {
    bc_lr: CudaFunction,
    bc_bt: CudaFunction,
    step: CudaFunction,
    q: CudaFunction,
    step_kr: CudaFunction,
}

fn load_kernels<R: Real>() -> Result<(Arc<CudaStream>, Kern)> {
    let ctx = CudaContext::new(0)?;
    let stream = ctx.default_stream();
    let ptx = compile_ptx_with_opts(
        include_str!("kernels.cu"),
        CompileOptions {
            arch: Some("sm_86"),
            fmad: Some(false),
            ..Default::default()
        },
    )?;
    let module = ctx.load_module(ptx)?;
    let nm = |base: &str| format!("{base}_{}", R::SUFFIX);
    let k = Kern {
        bc_lr: module.load_function(&nm("pf_bc_lr"))?,
        bc_bt: module.load_function(&nm("pf_bc_bt"))?,
        step: module.load_function(&nm("pf_step"))?,
        q: module.load_function(&nm("pf_q"))?,
        step_kr: module.load_function(&nm("pf_step_kr"))?,
    };
    Ok((stream, k))
}

fn cfg(len: usize) -> LaunchConfig {
    LaunchConfig {
        grid_dim: ((len as u64).div_ceil(BLOCK as u64).max(1) as u32, 1, 1),
        block_dim: (BLOCK, 1, 1),
        shared_mem_bytes: 0,
    }
}

fn edge_code(e: Edge) -> i32 {
    match e {
        Edge::Symmetry => 0,
        Edge::FarField => 1,
    }
}

/// GPU phase-field solver.
pub struct PhaseFieldGpu<R: Real> {
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
    stream: Arc<CudaStream>,
    k: Kern,
    p: CudaSlice<R>,
    u: CudaSlice<R>,
    pn: CudaSlice<R>,
    un: CudaSlice<R>,
    q: CudaSlice<R>,
}

impl<R: Real> PhaseFieldGpu<R> {
    /// Allocate a GPU solver with the same parameters and initial fields as a
    /// CPU [`solidification::phase_field::PhaseField`]. The CPU solver's `p`,
    /// `u`, boundary conditions and `u_far` are copied.
    pub fn from_cpu(cpu: &solidification::phase_field::PhaseField) -> Result<Self> {
        let (stream, k) = load_kernels::<R>()?;
        let n = (cpu.nx + 2) * (cpu.ny + 2);
        let host_p: Vec<R> = cpu.p.iter().map(|v| R::from_f64(*v)).collect();
        let host_u: Vec<R> = cpu.u.iter().map(|v| R::from_f64(*v)).collect();
        let gw = Self {
            nx: cpu.nx,
            ny: cpu.ny,
            dx: cpu.dx,
            d0: cpu.d0,
            d: cpu.d,
            w0: cpu.w0,
            eps4: cpu.eps4,
            lambda: cpu.lambda,
            tau0: cpu.tau0,
            bc: cpu.bc,
            u_far: cpu.u_far,
            p: stream.clone_htod(&host_p)?,
            u: stream.clone_htod(&host_u)?,
            pn: stream.alloc_zeros::<R>(n)?,
            un: stream.alloc_zeros::<R>(n)?,
            q: stream.alloc_zeros::<R>(n)?,
            stream,
            k,
        };
        Ok(gw)
    }

    fn fill_bc(&mut self) -> Result<()> {
        let (nx, ny) = (self.nx, self.ny);
        let ufar = R::from_f64(self.u_far);
        let (bl, br) = (edge_code(self.bc.left), edge_code(self.bc.right));
        let (bb, bt) = (edge_code(self.bc.bottom), edge_code(self.bc.top));
        let (nxr, nyr) = (nx as i32, ny as i32);
        {
            let mut lb = self.stream.launch_builder(&self.k.bc_lr);
            lb.arg(&nxr).arg(&nyr).arg(&bl).arg(&br).arg(&ufar).arg(&mut self.p).arg(&mut self.u);
            unsafe { lb.launch(cfg(ny + 2))? };
        }
        {
            let mut lb = self.stream.launch_builder(&self.k.bc_bt);
            lb.arg(&nxr).arg(&nyr).arg(&bb).arg(&bt).arg(&ufar).arg(&mut self.p).arg(&mut self.u);
            unsafe { lb.launch(cfg(nx + 2))? };
        }
        Ok(())
    }

    /// One explicit forward-Euler step with time step `dt` (the operator of
    /// `solidification/src/phase_field.rs`, which omits the variational
    /// anisotropy term of KR Eq. 88). Kept for the GPU = CPU check.
    pub fn step(&mut self, dt: f64) -> Result<()> {
        self.fill_bc()?;
        let nxr = self.nx as i32;
        let nyr = self.ny as i32;
        let dx = R::from_f64(self.dx);
        let inv_twodx = R::from_f64(1.0 / (2.0 * self.dx));
        let inv_fourd2 = R::from_f64(1.0 / (4.0 * self.dx * self.dx));
        let w0 = R::from_f64(self.w0);
        let eps4 = R::from_f64(self.eps4);
        let lambda = R::from_f64(self.lambda);
        let tau0 = R::from_f64(self.tau0);
        let d = R::from_f64(self.d);
        let dtr = R::from_f64(dt);
        let mut lb = self.stream.launch_builder(&self.k.step);
        lb.arg(&nxr).arg(&nyr).arg(&dx).arg(&inv_twodx).arg(&inv_fourd2).arg(&w0).arg(&eps4)
            .arg(&lambda).arg(&tau0).arg(&d).arg(&dtr)
            .arg(&self.p).arg(&self.u).arg(&mut self.pn).arg(&mut self.un);
        unsafe { lb.launch(cfg(self.nx * self.ny))? };
        std::mem::swap(&mut self.p, &mut self.pn);
        std::mem::swap(&mut self.u, &mut self.un);
        Ok(())
    }

    /// One explicit step with the **full** Karma-Rappel anisotropic operator
    /// (KR 1998 Eq. 88), which adds the variational term
    /// `p_x d_y Q - p_y d_x Q`, `Q = W dW/dtheta`, to the 9-point isotropic
    /// Laplacian. `solidification`'s CPU solver omits both, so this step does
    /// *not* equal the CPU step; it is the operator whose thin-interface limit
    /// is the free-boundary problem that KR benchmarked. The benchmark example
    /// uses this.
    pub fn step_kr(&mut self, dt: f64) -> Result<()> {
        self.step_kr_mode(dt, 3)
    }

    /// As [`step_kr`], but `mode` selects the discretisation for comparisons:
    /// bit 0 = use the 9-point isotropic Laplacian (else 5-point),
    /// bit 1 = include the variational anisotropy term. The CPU solver's
    /// operator is mode = 0 (and differs by using rayons it is not on device).
    pub fn step_kr_mode(&mut self, dt: f64, mode: i32) -> Result<()> {
        self.fill_bc()?;
        let nxr = self.nx as i32;
        let nyr = self.ny as i32;
        let inv_twodx = R::from_f64(1.0 / (2.0 * self.dx));
        let w0 = R::from_f64(self.w0);
        let eps4 = R::from_f64(self.eps4);
        {
            let mut lb = self.stream.launch_builder(&self.k.q);
            lb.arg(&nxr).arg(&nyr).arg(&inv_twodx).arg(&w0).arg(&eps4).arg(&self.p).arg(&mut self.q);
            unsafe { lb.launch(cfg(self.nx * self.ny))? };
        }
        let dx = R::from_f64(self.dx);
        let inv_fourd2 = R::from_f64(1.0 / (4.0 * self.dx * self.dx));
        let lambda = R::from_f64(self.lambda);
        let tau0 = R::from_f64(self.tau0);
        let d = R::from_f64(self.d);
        let dtr = R::from_f64(dt);
        let mut lb = self.stream.launch_builder(&self.k.step_kr);
        lb.arg(&nxr).arg(&nyr).arg(&dx).arg(&inv_twodx).arg(&inv_fourd2).arg(&w0).arg(&eps4)
            .arg(&lambda).arg(&tau0).arg(&d).arg(&dtr)
            .arg(&mode)
            .arg(&self.p).arg(&self.u).arg(&self.q).arg(&mut self.pn).arg(&mut self.un);
        unsafe { lb.launch(cfg(self.nx * self.ny))? };
        std::mem::swap(&mut self.p, &mut self.pn);
        std::mem::swap(&mut self.u, &mut self.un);
        Ok(())
    }

    /// Advance `n` corrected-KR steps of size `dt`.
    pub fn run_kr(&mut self, n: usize, dt: f64) -> Result<()> {
        self.run_kr_mode(n, dt, 3)
    }

    /// Advance `n` steps with an explicit stencil `mode` (see [`step_kr_mode`]).
    pub fn run_kr_mode(&mut self, n: usize, dt: f64, mode: i32) -> Result<()> {
        for _ in 0..n {
            self.step_kr_mode(dt, mode)?;
        }
        Ok(())
    }

    /// Advance `n` steps of size `dt`.
    pub fn run(&mut self, n: usize, dt: f64) -> Result<()> {
        for _ in 0..n {
            self.step(dt)?;
        }
        Ok(())
    }

    /// Stability limit, identical to the CPU `PhaseField::stable_dt`.
    pub fn stable_dt(&self) -> f64 {
        let dt_p = self.dx * self.dx * self.tau0 / (4.0 * self.w0 * self.w0);
        let dt_u = self.dx * self.dx / (4.0 * self.d);
        0.4 * dt_p.min(dt_u)
    }

    /// Download `p` (with ghosts), in the CPU layout (`idx = j*(nx+2)+i`).
    pub fn download_p(&self) -> Result<Vec<f64>> {
        Ok(self.stream.clone_dtoh(&self.p)?.into_iter().map(R::to_f64).collect())
    }
    /// Download `u` (with ghosts).
    pub fn download_u(&self) -> Result<Vec<f64>> {
        Ok(self.stream.clone_dtoh(&self.u)?.into_iter().map(R::to_f64).collect())
    }

    /// Upload `p` (with ghosts); used for restart.
    pub fn upload_p(&mut self, p: &[f64]) -> Result<()> {
        let host: Vec<R> = p.iter().map(|v| R::from_f64(*v)).collect();
        self.p = self.stream.clone_htod(&host)?;
        Ok(())
    }
    /// Upload `u` (with ghosts).
    pub fn upload_u(&mut self, u: &[f64]) -> Result<()> {
        let host: Vec<R> = u.iter().map(|v| R::from_f64(*v)).collect();
        self.u = self.stream.clone_htod(&host)?;
        Ok(())
    }

    pub fn synchronize(&self) -> Result<()> {
        self.stream.synchronize()?;
        Ok(())
    }
}

/// Host-side interface analysis on a downloaded phase field (CPU layout).
pub struct HostView<'a> {
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub p: &'a [f64],
}

impl<'a> HostView<'a> {
    #[inline]
    fn idx(&self, i: usize, j: usize) -> usize {
        j * (self.nx + 2) + i
    }
    #[inline]
    fn x(&self, i: usize) -> f64 {
        (i as f64 - 0.5) * self.dx
    }
    #[inline]
    fn y(&self, j: usize) -> f64 {
        (j as f64 - 0.5) * self.dx
    }

    /// Largest `x` interior cell with `p >= 0.5`, over all rows.
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

    /// Interpolated `x` where `p` crosses 0.5 going outward along row `j`.
    pub fn crossing_x(&self, j: usize) -> Option<f64> {
        let mut prev = self.p[self.idx(1, j)];
        for i in 2..=self.nx {
            let cur = self.p[self.idx(i, j)];
            if prev >= 0.5 && cur < 0.5 {
                let x0 = self.x(i - 1);
                let x1 = self.x(i);
                let t = (prev - 0.5) / (prev - cur);
                return Some(x0 + t * (x1 - x0));
            }
            prev = cur;
        }
        None
    }

    /// Tip radius from a least-squares fit `x(y) = x0 - y^2/(2R)` to the
    /// `p = 0.5` contour near the axis (same estimator as the CPU solver).
    pub fn tip_radius(&self) -> f64 {
        let mut ys = Vec::new();
        let mut xs = Vec::new();
        let kmax = self.ny.min(8);
        for j in 1..=kmax {
            if let Some(x) = self.crossing_x(j) {
                let y = self.y(j) - self.y(1);
                ys.push(y * y);
                xs.push(x);
            } else {
                break;
            }
        }
        if ys.len() < 4 {
            return f64::NAN;
        }
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

    /// Total solid fraction over the interior.
    pub fn solid_fraction(&self) -> f64 {
        let n = (self.nx * self.ny) as f64;
        let mut s = 0.0;
        for j in 1..=self.ny {
            for i in 1..=self.nx {
                s += 0.5 * (1.0 + self.p[self.idx(i, j)]);
            }
        }
        s / n
    }
}
