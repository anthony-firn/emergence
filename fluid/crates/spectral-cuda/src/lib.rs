//! CUDA + cuFFT backend: the same algorithm as `spectral_core::Solver`
//! (rotational form, 2/3 dealiasing, Lawson IF-RK4), with f32 storage and
//! f64 accumulation in reductions. Everything stays on the device; the host
//! only sees diagnostics and explicit downloads.
//!
//! Runtime libraries are loaded dynamically: put the directories holding
//! `libcufft.so.11` and `libnvrtc.so.12` on `LD_LIBRARY_PATH`.

use std::sync::Arc;

use cudarc::cufft::{sys::cufftType, sys::float2, CudaFft};
use cudarc::driver::{CudaContext, CudaFunction, CudaSlice, CudaStream, LaunchConfig, PushKernelArg};
use cudarc::nvrtc::{compile_ptx_with_opts, CompileOptions};
use spectral_core::{solver::initial_spectrum, Backend, Diagnostics, Grid, Penalization, Surface};

mod certify_gpu;
mod imex;
mod surface_gpu;

pub use certify_gpu::CudaCertifier;

pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

const BLOCK: u32 = 256;
const REDUCE_BLOCKS: u32 = 1024;

struct Kernels {
    curl: CudaFunction,
    cross: CudaFunction,
    project: CudaFunction,
    combine5: CudaFunction,
    diag: CudaFunction,
    max_speed2: CudaFunction,
    phase_copy: CudaFunction,
    accum_if: CudaFunction,
    penalize: CudaFunction,
    penalize_field: CudaFunction,
    solid_force: CudaFunction,
}

/// Device copy of a [`Penalization`].
struct DevPenal {
    chi: CudaSlice<f32>,
    sponge: CudaSlice<f32>,
    eta: f64,
    eta_sponge: f64,
    u_solid: [f32; 3],
    u_inf: [f32; 3],
    implicit: bool,
    /// spatially varying χ-term target (see Penalization::u_target)
    target: Option<[CudaSlice<f32>; 3]>,
}

/// Device workspaces used inside one right-hand-side evaluation.
struct Work {
    stream: Arc<CudaStream>,
    k: Kernels,
    r2c: CudaFft,
    c2r: CudaFft,
    n: i32,
    nh: i32,
    tmp: [CudaSlice<float2>; 3],
    up: [CudaSlice<f32>; 3],
    wp: [CudaSlice<f32>; 3],
    reduce_f64: CudaSlice<f64>,
    reduce_f32: CudaSlice<f32>,
    penal: Option<DevPenal>,
    /// spherical truncation radius² (0 = 2/3 cube)
    kc2: f32,
}

pub struct CudaSolver {
    pub n: usize,
    pub nu: f64,
    pub t: f64,
    w: Work,
    u: [CudaSlice<float2>; 3],
    st: [CudaSlice<float2>; 3],
    /// current stage tendency (low-storage IF-RK4 reuses one buffer)
    k: [CudaSlice<float2>; 3],
    /// running uⁿ⁺¹ accumulator
    acc: [CudaSlice<float2>; 3],
    spherical: bool,
    dx: f64,
    rng: u64,
    forcing: Option<[CudaSlice<float2>; 3]>,
    /// Uniform body force driving a free mean flow; `None` holds the mean fixed
    /// (see `spectral_core::Solver::body_force`).
    pub body_force: Option<[f64; 3]>,
    /// max |u| measured at the start of the latest step (lags by ≤ 1 step)
    last_umax: Option<f64>,
    ik: imex::ImexKernels,
    imex: Option<Box<imex::ImexDev>>,
    nk: surface_gpu::NufftKernels,
    surf: Option<Box<surface_gpu::DevSurface>>,
    grid: Grid,
    module: Arc<cudarc::driver::CudaModule>,
}

fn cfg(len: usize) -> LaunchConfig {
    LaunchConfig {
        grid_dim: ((len as u64).div_ceil(BLOCK as u64) as u32, 1, 1),
        block_dim: (BLOCK, 1, 1),
        shared_mem_bytes: 0,
    }
}

fn alloc3<T: cudarc::driver::DeviceRepr + cudarc::driver::ValidAsZeroBits>(
    s: &Arc<CudaStream>,
    len: usize,
) -> Result<[CudaSlice<T>; 3]> {
    Ok([s.alloc_zeros(len)?, s.alloc_zeros(len)?, s.alloc_zeros(len)?])
}

impl Work {
    fn ns(&self) -> usize {
        (self.n * self.n * self.nh) as usize
    }
    fn nr(&self) -> usize {
        (self.n * self.n * self.n) as usize
    }

    /// spectral component `src` -> physical `dst` (via a copy; c2r clobbers its input).
    fn to_phys(
        stream: &Arc<CudaStream>,
        c2r: &CudaFft,
        src: &CudaSlice<float2>,
        tmp: &mut CudaSlice<float2>,
        dst: &mut CudaSlice<f32>,
    ) -> Result<()> {
        stream.memcpy_dtod(src, tmp)?;
        c2r.exec_c2r(tmp, dst)?;
        Ok(())
    }

    /// out = P[(u × ω)^] + P f̂ for state `uh`.
    fn rhs(
        &mut self,
        uh: &[CudaSlice<float2>; 3],
        out: &mut [CudaSlice<float2>; 3],
        forcing: Option<&[CudaSlice<float2>; 3]>,
        body_force: Option<[f64; 3]>,
        shift: Option<[f64; 3]>,
    ) -> Result<()> {
        let (n, nh) = (self.n, self.nh);
        let (ns, nr) = (self.ns(), self.nr());
        let s = self.stream.clone();
        let sh = shift.unwrap_or([0.0; 3]).map(|v| v as f32);
        for c in 0..3 {
            if shift.is_some() {
                unsafe {
                    s.launch_builder(&self.k.phase_copy)
                        .arg(&mut self.tmp[c]).arg(&uh[c])
                        .arg(&sh[0]).arg(&sh[1]).arg(&sh[2]).arg(&n).arg(&nh)
                        .launch(cfg(ns))?;
                }
                self.c2r.exec_c2r(&mut self.tmp[c], &mut self.up[c])?;
            } else {
                Self::to_phys(&s, &self.c2r, &uh[c], &mut self.tmp[c], &mut self.up[c])?;
            }
        }
        {
            let [t0, t1, t2] = &mut self.tmp;
            unsafe {
                s.launch_builder(&self.k.curl)
                    .arg(&uh[0]).arg(&uh[1]).arg(&uh[2])
                    .arg(t0).arg(t1).arg(t2)
                    .arg(&n).arg(&nh)
                    .arg(&sh[0]).arg(&sh[1]).arg(&sh[2])
                    .launch(cfg(ns))?;
            }
        }
        for c in 0..3 {
            self.c2r.exec_c2r(&mut self.tmp[c], &mut self.wp[c])?;
        }
        {
            let [w0, w1, w2] = &mut self.wp;
            let nr64 = nr as i64;
            unsafe {
                s.launch_builder(&self.k.cross)
                    .arg(&self.up[0]).arg(&self.up[1]).arg(&self.up[2])
                    .arg(&mut *w0).arg(&mut *w1).arg(&mut *w2)
                    .arg(&nr64)
                    .launch(cfg(nr))?;
            }
            if let Some(pen) = self.penal.as_ref().filter(|p| p.target.is_some()) {
                let tg = pen.target.as_ref().unwrap();
                let ie = (1.0 / pen.eta) as f32;
                let is = (1.0 / pen.eta_sponge) as f32;
                let ui = pen.u_inf;
                unsafe {
                    s.launch_builder(&self.k.penalize_field)
                        .arg(&self.up[0]).arg(&self.up[1]).arg(&self.up[2])
                        .arg(&mut *w0).arg(&mut *w1).arg(&mut *w2)
                        .arg(&pen.chi).arg(&pen.sponge)
                        .arg(&tg[0]).arg(&tg[1]).arg(&tg[2])
                        .arg(&ie).arg(&is)
                        .arg(&ui[0]).arg(&ui[1]).arg(&ui[2])
                        .arg(&nr64)
                        .launch(cfg(nr))?;
                }
            } else if let Some(pen) = &self.penal {
                // the solid term is handled by the implicit solve in IMEX mode
                let ie = if pen.implicit { 0.0f32 } else { (1.0 / pen.eta) as f32 };
                let is = (1.0 / pen.eta_sponge) as f32;
                let (us, ui) = (pen.u_solid, pen.u_inf);
                unsafe {
                    s.launch_builder(&self.k.penalize)
                        .arg(&self.up[0]).arg(&self.up[1]).arg(&self.up[2])
                        .arg(&mut *w0).arg(&mut *w1).arg(&mut *w2)
                        .arg(&pen.chi).arg(&pen.sponge)
                        .arg(&ie).arg(&is)
                        .arg(&us[0]).arg(&us[1]).arg(&us[2])
                        .arg(&ui[0]).arg(&ui[1]).arg(&ui[2])
                        .arg(&nr64)
                        .launch(cfg(nr))?;
                }
            }
        }
        for c in 0..3 {
            self.r2c.exec_r2c(&self.wp[c], &mut out[c])?;
        }
        let scale = 1.0f32 / nr as f32;
        let has_f: i32 = forcing.is_some() as i32;
        let [o0, o1, o2] = out;
        // When unforced, pass the up/tmp buffers as dummies (never read).
        let f = forcing.unwrap_or(&self.tmp);
        let keep: i32 = body_force.is_some() as i32;
        let bf = body_force.unwrap_or([0.0; 3]).map(|v| v as f32);
        unsafe {
            s.launch_builder(&self.k.project)
                .arg(o0).arg(o1).arg(o2)
                .arg(&f[0]).arg(&f[1]).arg(&f[2])
                .arg(&has_f).arg(&scale).arg(&n).arg(&nh)
                .arg(&keep).arg(&bf[0]).arg(&bf[1]).arg(&bf[2])
                .arg(&self.kc2).arg(&sh[0]).arg(&sh[1]).arg(&sh[2])
                .launch(cfg(ns))?;
        }
        Ok(())
    }

    /// max |u|² over the physical velocity currently held in `up`.
    fn max_speed2_of_up(&mut self) -> Result<f64> {
        let nr64 = self.nr() as i64;
        let s = self.stream.clone();
        let lc = LaunchConfig {
            grid_dim: (REDUCE_BLOCKS, 1, 1),
            block_dim: (BLOCK, 1, 1),
            shared_mem_bytes: BLOCK * 4,
        };
        unsafe {
            s.launch_builder(&self.k.max_speed2)
                .arg(&self.up[0]).arg(&self.up[1]).arg(&self.up[2])
                .arg(&mut self.reduce_f32)
                .arg(&nr64)
                .launch(lc)?;
        }
        let parts = s.clone_dtoh(&self.reduce_f32)?;
        Ok(parts.into_iter().fold(0.0f32, f32::max) as f64)
    }
}

impl CudaSolver {
    pub fn new(grid: Grid, nu: f64, u_phys: [Vec<f64>; 3]) -> Result<Self> {
        let ctx = CudaContext::new(0)?;
        let stream = ctx.default_stream();
        let ptx = compile_ptx_with_opts(
            include_str!("kernels.cu"),
            CompileOptions { arch: Some("sm_86"), ..Default::default() },
        )?;
        let module = ctx.load_module(ptx)?;
        let ik = imex::ImexKernels::load(&module)?;
        let nk = surface_gpu::NufftKernels::load(&module)?;
        let k = Kernels {
            curl: module.load_function("curl")?,
            cross: module.load_function("cross_inplace")?,
            project: module.load_function("project_scale")?,
            combine5: module.load_function("combine5")?,
            diag: module.load_function("diag_partial")?,
            max_speed2: module.load_function("max_speed2_partial")?,
            phase_copy: module.load_function("phase_copy")?,
            accum_if: module.load_function("accum_if")?,
            penalize: module.load_function("penalize")?,
            penalize_field: module.load_function("penalize_field")?,
            solid_force: module.load_function("solid_force_partial")?,
        };
        let n = grid.n;
        let (ns, nr) = (grid.spec_len(), grid.real_len());
        let ni = n as i32;
        let r2c = CudaFft::plan_3d(ni, ni, ni, cufftType::CUFFT_R2C, stream.clone())?;
        let c2r = CudaFft::plan_3d(ni, ni, ni, cufftType::CUFFT_C2R, stream.clone())?;

        // Initial spectrum on the host in f64 (dealiased, projected), then upload.
        let grid_keep = grid.clone();
        let u0 = initial_spectrum(&grid, &u_phys);
        drop(u_phys);
        let mut u = alloc3::<float2>(&stream, ns)?;
        for c in 0..3 {
            let host: Vec<float2> =
                u0[c].iter().map(|z| float2 { x: z.re as f32, y: z.im as f32 }).collect();
            stream.memcpy_htod(&host, &mut u[c])?;
        }

        let w = Work {
            k,
            r2c,
            c2r,
            n: ni,
            nh: grid.nh as i32,
            tmp: alloc3(&stream, ns)?,
            up: alloc3(&stream, nr)?,
            wp: alloc3(&stream, nr)?,
            reduce_f64: stream.alloc_zeros(3 * REDUCE_BLOCKS as usize)?,
            reduce_f32: stream.alloc_zeros(REDUCE_BLOCKS as usize)?,
            penal: None,
            kc2: if grid.spherical { grid.kc2 as f32 } else { 0.0 },
            stream: stream.clone(),
        };
        Ok(Self {
            n,
            nu,
            t: 0.0,
            u,
            st: alloc3(&stream, ns)?,
            k: alloc3(&stream, ns)?,
            acc: alloc3(&stream, ns)?,
            spherical: grid.spherical,
            dx: grid.dx(),
            rng: 0x9E37_79B9_7F4A_7C15,
            forcing: None,
            body_force: None,
            last_umax: None,
            ik,
            imex: None,
            nk,
            surf: None,
            grid: grid_keep,
            module: module.clone(),
            w,
        })
    }

    /// out[c] = Σ_m coef[m]·E(h[m])·x_m[c] for c = 0..3, E(h) = e^{-ν|k|²h}.
    /// Terms with coef 0 are not read.
    #[allow(clippy::too_many_arguments)]
    fn combine(
        w: &Work,
        nu: f64,
        out: &mut [CudaSlice<float2>; 3],
        xs: [&[CudaSlice<float2>; 3]; 5],
        coef: [f64; 5],
        h: [f64; 5],
    ) -> Result<()> {
        let (n, nh) = (w.n, w.nh);
        let nu = nu as f32;
        let c: Vec<f32> = coef.iter().map(|&v| v as f32).collect();
        let hh: Vec<f32> = h.iter().map(|&v| v as f32).collect();
        for comp in 0..3 {
            unsafe {
                w.stream
                    .launch_builder(&w.k.combine5)
                    .arg(&mut out[comp])
                    .arg(&xs[0][comp]).arg(&xs[1][comp]).arg(&xs[2][comp])
                    .arg(&xs[3][comp]).arg(&xs[4][comp])
                    .arg(&c[0]).arg(&c[1]).arg(&c[2]).arg(&c[3]).arg(&c[4])
                    .arg(&hh[0]).arg(&hh[1]).arg(&hh[2]).arg(&hh[3]).arg(&hh[4])
                    .arg(&nu).arg(&n).arg(&nh)
                    .launch(cfg(w.ns()))?;
            }
        }
        Ok(())
    }

    /// out[c] += coef·E(h)·x[c], E(h) = e^{-ν|k|²h} (unbiased f32 form).
    fn accum_if(w: &Work, nu: f64, out: &mut [CudaSlice<float2>; 3], x: &[CudaSlice<float2>; 3], coef: f64, h: f64) -> Result<()> {
        let (n, nh, c, hh, nu) = (w.n, w.nh, coef as f32, h as f32, nu as f32);
        for comp in 0..3 {
            unsafe {
                w.stream
                    .launch_builder(&w.k.accum_if)
                    .arg(&mut out[comp]).arg(&x[comp])
                    .arg(&c).arg(&hh).arg(&nu).arg(&n).arg(&nh)
                    .launch(cfg(w.ns()))?;
            }
        }
        Ok(())
    }

    /// Next random grid shift in [0, dx)³ — the same xorshift sequence as
    /// `spectral_core::Solver`, so both backends shift identically.
    fn next_shift(&mut self) -> [f64; 3] {
        let dx = self.dx;
        let mut r = || {
            self.rng ^= self.rng << 13;
            self.rng ^= self.rng >> 7;
            self.rng ^= self.rng << 17;
            (self.rng >> 11) as f64 / (1u64 << 53) as f64 * dx
        };
        [r(), r(), r()]
    }

    /// One Lawson IF-RK4 step, low storage: the stage tendency reuses one
    /// buffer and uⁿ⁺¹ is accumulated as stages finish (4 spectral sets
    /// instead of 6, which is what lets 384³ fit in 6 GB).
    pub fn try_step(&mut self, h: f64) -> Result<()> {
        if self.surf.is_some() || self.w.penal.as_ref().is_some_and(|p| p.implicit) {
            return self.try_step_imex(h);
        }
        // phase-shift dealiasing on spherical grids (see spectral_core::Grid::new_spherical)
        let shifts: [Option<[f64; 3]>; 4] = if self.spherical {
            let d = self.next_shift();
            let hc = 0.5 * self.dx;
            let d2 = [d[0] + hc, d[1] + hc, d[2] + hc];
            [Some(d), Some(d2), Some(d), Some(d2)]
        } else {
            [None; 4]
        };
        let nu = self.nu;
        let f = self.forcing.as_ref();
        let bf = self.body_force;
        // a = N(u); record max|u| from the physical velocity just built
        self.w.rhs(&self.u, &mut self.k, f, bf, shifts[0])?;
        self.last_umax = Some(self.w.max_speed2_of_up()?.sqrt());
        let u = &self.u;
        // acc = E(h) u + h/6 E(h) a ;  st1 = E(h/2) u + h/2 E(h/2) a
        Self::combine(&self.w, nu, &mut self.acc, [u, &self.k, u, u, u], [1.0, h / 6.0, 0.0, 0.0, 0.0], [h, h, 0.0, 0.0, 0.0])?;
        Self::combine(&self.w, nu, &mut self.st, [u, &self.k, u, u, u], [1.0, 0.5 * h, 0.0, 0.0, 0.0], [0.5 * h, 0.5 * h, 0.0, 0.0, 0.0])?;
        // b: acc += h/3 E(h/2) b ;  st2 = E(h/2) u + h/2 b
        self.w.rhs(&self.st, &mut self.k, f, bf, shifts[1])?;
        Self::accum_if(&self.w, nu, &mut self.acc, &self.k, h / 3.0, 0.5 * h)?;
        Self::combine(&self.w, nu, &mut self.st, [u, &self.k, u, u, u], [1.0, 0.5 * h, 0.0, 0.0, 0.0], [0.5 * h, 0.0, 0.0, 0.0, 0.0])?;
        // c: acc += h/3 E(h/2) c ;  st3 = E(h) u + h E(h/2) c
        self.w.rhs(&self.st, &mut self.k, f, bf, shifts[2])?;
        Self::accum_if(&self.w, nu, &mut self.acc, &self.k, h / 3.0, 0.5 * h)?;
        Self::combine(&self.w, nu, &mut self.st, [u, &self.k, u, u, u], [1.0, h, 0.0, 0.0, 0.0], [h, 0.5 * h, 0.0, 0.0, 0.0])?;
        // d: acc += h/6 d
        self.w.rhs(&self.st, &mut self.k, f, bf, shifts[3])?;
        Self::accum_if(&self.w, nu, &mut self.acc, &self.k, h / 6.0, 0.0)?;
        std::mem::swap(&mut self.u, &mut self.acc);
        self.t += h;
        Ok(())
    }

    pub fn try_diagnostics(&mut self) -> Result<Diagnostics> {
        let w = &mut self.w;
        let (n, nh) = (w.n, w.nh);
        let lc = LaunchConfig {
            grid_dim: (REDUCE_BLOCKS, 1, 1),
            block_dim: (BLOCK, 1, 1),
            shared_mem_bytes: 3 * BLOCK * 8,
        };
        unsafe {
            w.stream
                .launch_builder(&w.k.diag)
                .arg(&self.u[0]).arg(&self.u[1]).arg(&self.u[2])
                .arg(&mut w.reduce_f64)
                .arg(&n).arg(&nh)
                .launch(lc)?;
        }
        let parts = w.stream.clone_dtoh(&w.reduce_f64)?;
        let (mut e, mut z, mut d) = (0.0, 0.0, 0.0f64);
        for p in parts.chunks(3) {
            e += p[0];
            z += p[1];
            d = d.max(p[2]);
        }
        Ok(Diagnostics {
            t: self.t,
            energy: e,
            enstrophy: z,
            dissipation: 2.0 * self.nu * z,
            max_div: d,
            tail_energy_fraction: f64::NAN, // not computed on device yet
        })
    }

    pub fn try_velocity_phys(&mut self, c: usize) -> Result<Vec<f64>> {
        let s = self.w.stream.clone();
        Work::to_phys(&s, &self.w.c2r, &self.u[c], &mut self.w.tmp[c], &mut self.w.up[c])?;
        Ok(s.clone_dtoh(&self.w.up[c])?.into_iter().map(f64::from).collect())
    }

    pub fn try_max_speed(&mut self) -> Result<f64> {
        let s = self.w.stream.clone();
        for c in 0..3 {
            Work::to_phys(&s, &self.w.c2r, &self.u[c], &mut self.w.tmp[c], &mut self.w.up[c])?;
        }
        Ok(self.w.max_speed2_of_up()?.sqrt())
    }

    /// Upload solid/sponge masks; replaces any previous penalization.
    pub fn set_penalization(&mut self, p: &Penalization) -> Result<()> {
        let s = &self.w.stream;
        let to32 = |v: &[f64]| v.iter().map(|&x| x as f32).collect::<Vec<f32>>();
        let f3 = |v: [f64; 3]| v.map(|x| x as f32);
        self.w.penal = Some(DevPenal {
            chi: s.clone_htod(&to32(&p.chi))?,
            sponge: s.clone_htod(&to32(&p.sponge))?,
            eta: p.eta,
            eta_sponge: p.eta_sponge,
            u_solid: f3(p.u_solid),
            u_inf: f3(p.u_inf),
            implicit: p.implicit,
            target: match &p.u_target {
                None => None,
                Some(t) => Some([s.clone_htod(&to32(&t[0]))?, s.clone_htod(&to32(&t[1]))?, s.clone_htod(&to32(&t[2]))?]),
            },
        });
        Ok(())
    }

    /// Replace the penalization's target field (see `Penalization::u_target`).
    pub fn set_penalization_target(&mut self, t: &[Vec<f64>; 3]) -> Result<()> {
        let s = self.w.stream.clone();
        let pen = self.w.penal.as_mut().ok_or("no penalization installed")?;
        let tg = match pen.target.as_mut() {
            Some(tg) => tg,
            None => {
                pen.target = Some(alloc3(&s, t[0].len())?);
                pen.target.as_mut().unwrap()
            }
        };
        for c in 0..3 {
            let v: Vec<f32> = t[c].iter().map(|&x| x as f32).collect();
            s.memcpy_htod(&v, &mut tg[c])?;
        }
        Ok(())
    }

    /// Force on each part of the installed surface (last implicit stage).
    pub fn surface_part_forces(&self) -> Option<Vec<[f64; 3]>> {
        self.surf.as_ref().map(|sd| sd.part_forces())
    }

    /// Replace the per-marker prescribed velocities of the installed surface.
    pub fn set_surface_targets(&mut self, t: &[[f64; 3]]) -> Result<()> {
        let sd = self.surf.as_mut().ok_or("no surface installed")?;
        sd.set_targets(t.to_vec());
        Ok(())
    }

    /// Set a steady body force from its spectral coefficients (f64 host arrays,
    /// N-grid layout); it is projected with the nonlinear term.
    pub fn set_forcing(&mut self, f: &[spectral_core::solver::Spec; 3]) -> Result<()> {
        let s = &self.w.stream;
        let mut dev = alloc3::<float2>(s, f[0].len())?;
        for c in 0..3 {
            let host: Vec<float2> = f[c].iter().map(|z| float2 { x: z.re as f32, y: z.im as f32 }).collect();
            s.memcpy_htod(&host, &mut dev[c])?;
        }
        self.forcing = Some(dev);
        Ok(())
    }

    /// Install surface-constrained walls (point forces at the markers).
    pub fn set_surface(&mut self, surface: &Surface) -> Result<()> {
        self.surf = Some(Box::new(surface_gpu::DevSurface::new(&self.w.stream, &self.grid, surface)?));
        Ok(())
    }

    /// Set the fixed uniform (k = 0) velocity; see `Solver::set_mean_flow`.
    pub fn set_mean_flow(&mut self, m: [f64; 3]) -> Result<()> {
        for c in 0..3 {
            let v = [float2 { x: m[c] as f32, y: 0.0 }];
            self.w.stream.memcpy_htod(&v, &mut self.u[c].slice_mut(0..1))?;
        }
        Ok(())
    }

    /// Spectral velocity component `c` (f64 host copy of the coefficients).
    pub fn try_spectrum(&self, c: usize) -> Result<Vec<spectral_core::solver::Spec>> {
        let v = self.w.stream.clone_dtoh(&self.u[c])?;
        Ok(vec![v.into_iter().map(|z| num_complex::Complex64::new(z.x as f64, z.y as f64)).collect()])
    }

    /// The uniform (k = 0) velocity: the mean flow over the box.
    pub fn try_mean_velocity(&self) -> Result<[f64; 3]> {
        let mut m = [0.0; 3];
        for c in 0..3 {
            let v = self.w.stream.clone_dtoh(&self.u[c].slice(0..1))?;
            m[c] = v[0].x as f64;
        }
        Ok(m)
    }

    /// Force of the fluid on the penalized solid, F = (1/η) ∫ χ (u - u_s) dx.
    pub fn try_solid_force(&mut self) -> Result<Option<[f64; 3]>> {
        if let Some(sd) = &self.surf {
            return Ok(Some(sd.last_force));
        }
        match &self.w.penal {
            None => return Ok(None),
            Some(p) if p.implicit => return Ok(self.imex.as_ref().map(|d| d.last_force)),
            _ => {}
        }
        let s = self.w.stream.clone();
        for c in 0..3 {
            Work::to_phys(&s, &self.w.c2r, &self.u[c], &mut self.w.tmp[c], &mut self.w.up[c])?;
        }
        let w = &mut self.w;
        let pen = w.penal.as_ref().unwrap();
        let nr64 = (w.n as i64).pow(3);
        let lc = LaunchConfig {
            grid_dim: (REDUCE_BLOCKS, 1, 1),
            block_dim: (BLOCK, 1, 1),
            shared_mem_bytes: 3 * BLOCK * 8,
        };
        let us = pen.u_solid;
        unsafe {
            s.launch_builder(&w.k.solid_force)
                .arg(&w.up[0]).arg(&w.up[1]).arg(&w.up[2])
                .arg(&pen.chi)
                .arg(&us[0]).arg(&us[1]).arg(&us[2])
                .arg(&mut w.reduce_f64)
                .arg(&nr64)
                .launch(lc)?;
        }
        let parts = s.clone_dtoh(&w.reduce_f64)?;
        let mut f = [0.0f64; 3];
        for p in parts.chunks(3) {
            for c in 0..3 {
                f[c] += p[c];
            }
        }
        let dv = (2.0 * std::f64::consts::PI / w.n as f64).powi(3);
        Ok(Some(f.map(|v| v * dv / pen.eta)))
    }

    pub fn synchronize(&self) -> Result<()> {
        self.w.stream.synchronize()?;
        Ok(())
    }
}

impl Backend for CudaSolver {
    fn name(&self) -> &'static str {
        "cuda-f32"
    }
    fn n(&self) -> usize {
        self.n
    }
    fn nu(&self) -> f64 {
        self.nu
    }
    fn time(&self) -> f64 {
        self.t
    }
    fn step(&mut self, h: f64) {
        self.try_step(h).expect("CUDA step failed")
    }
    fn diagnostics(&mut self) -> Diagnostics {
        self.try_diagnostics().expect("CUDA diagnostics failed")
    }
    fn cfl_dt(&mut self, cfl: f64) -> f64 {
        let umax = match self.last_umax {
            Some(v) => v,
            None => self.try_max_speed().expect("CUDA max speed failed"),
        };
        cfl * 2.0 * std::f64::consts::PI / self.n as f64 / umax.max(1e-12)
    }
    fn velocity_phys(&mut self, c: usize) -> Vec<f64> {
        self.try_velocity_phys(c).expect("CUDA download failed")
    }
    fn set_penalization(&mut self, p: Penalization) {
        CudaSolver::set_penalization(self, &p).expect("CUDA mask upload failed")
    }
    fn set_mean_flow(&mut self, m: [f64; 3]) {
        CudaSolver::set_mean_flow(self, m).expect("CUDA mean flow failed")
    }
    fn solid_force(&mut self) -> Option<[f64; 3]> {
        self.try_solid_force().expect("CUDA solid force failed")
    }
    fn set_surface(&mut self, s: Surface) {
        CudaSolver::set_surface(self, &s).expect("CUDA surface setup failed")
    }
}
