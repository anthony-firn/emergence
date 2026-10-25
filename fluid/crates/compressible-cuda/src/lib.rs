//! GPU (CUDA, NVRTC) port of the compressible hybrid central/WCNS solver
//! on uniform grids: Euler and the five-equation multiphase model, in f64
//! (bit-for-bit comparable with the CPU solver up to summation order) or
//! f32 (fast on consumer GPUs).
//!
//! The device state is built from a CPU `compressible::Solver` (same cell
//! layout with ghosts, converted to structure-of-arrays) and can be
//! downloaded back for comparison. Kernels: `src/kernels.cu`.

use compressible::eos::Nasg;
use compressible::scheme::FluxMode;
use compressible::solver::{Bc, Mesh};
use compressible::{Solver, System};
use cudarc::driver::{CudaContext, CudaFunction, CudaSlice, CudaStream, DeviceRepr, LaunchConfig, PushKernelArg, ValidAsZeroBits};
use cudarc::nvrtc::{compile_ptx_with_opts, CompileOptions};
use std::sync::Arc;

pub mod amr;

pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

const BLOCK: u32 = 128;
const RATE_BLOCKS: u32 = 256;

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

/// Which system the kernels solve.
#[derive(Clone, Copy, Debug)]
pub enum GpuSystem {
    Euler(Nasg),
    FiveEq { eos: [Nasg; 2], kapila: bool },
}

impl GpuSystem {
    pub(crate) fn name(&self) -> &'static str {
        match self {
            GpuSystem::Euler(_) => "euler",
            GpuSystem::FiveEq { .. } => "five",
        }
    }
    pub(crate) fn nv(&self) -> usize {
        match self {
            GpuSystem::Euler(_) => 5,
            GpuSystem::FiveEq { .. } => 7,
        }
    }
    pub(crate) fn mom(&self, d: usize) -> i32 {
        match self {
            GpuSystem::Euler(_) => 1 + d as i32,
            GpuSystem::FiveEq { .. } => 2 + d as i32,
        }
    }
    pub(crate) fn eos_args(&self) -> ([f64; 8], i32) {
        match self {
            GpuSystem::Euler(e) => ([e.gamma, e.p_inf, e.b, e.q, 1.4, 0.0, 0.0, 0.0], 0),
            GpuSystem::FiveEq { eos: [a, b], kapila } => ([a.gamma, a.p_inf, a.b, a.q, b.gamma, b.p_inf, b.b, b.q], *kapila as i32),
        }
    }
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct Geo {
    pub m: [i32; 3],
    pub n: [i32; 3],
    pub g: [i32; 3],
    _pad: i32,
    pub stride: [i64; 3],
    /// cells per block (with ghosts)
    pub ncell: i64,
    /// variable stride of the field arrays (all blocks)
    pub ntot: i64,
    /// variable stride of the right-hand side
    pub rtot: i64,
    /// blocks in a launch, first right-hand-side slot
    pub nb: i32,
    pub s0: i32,
}
unsafe impl DeviceRepr for Geo {}

pub(crate) struct Kern {
    pub(crate) theta: CudaFunction,
    pub(crate) midflux: CudaFunction,
    pub(crate) flux: CudaFunction,
    pub(crate) accum: CudaFunction,
    pub(crate) prim: CudaFunction,
    pub(crate) rate: CudaFunction,
    pub(crate) rk: CudaFunction,
    pub(crate) ghost: CudaFunction,
    pub(crate) amr_ghost: CudaFunction,
    pub(crate) amr_rgather: CudaFunction,
    pub(crate) amr_pgather: CudaFunction,
    pub(crate) amr_ppass: CudaFunction,
    pub(crate) amr_rpass: CudaFunction,
    pub(crate) amr_faces: CudaFunction,
    pub(crate) amr_detail: CudaFunction,
    pub(crate) amr_reflux: CudaFunction,
}

/// Compile the kernels and load the entry points of one system / precision.
pub(crate) fn load_kernels<R: Real>(sys: GpuSystem) -> Result<(Arc<CudaStream>, Kern)> {
    let ctx = CudaContext::new(0)?;
    let stream = ctx.default_stream();
    let ptx = compile_ptx_with_opts(include_str!("kernels.cu"), CompileOptions { arch: Some("sm_86"), ..Default::default() })?;
    let module = ctx.load_module(ptx)?;
    let nm = |base: &str| format!("{base}_{}_{}", sys.name(), R::SUFFIX);
    let k = Kern {
        theta: module.load_function(&nm("theta"))?,
        midflux: module.load_function(&nm("midflux"))?,
        flux: module.load_function(&nm("flux"))?,
        accum: module.load_function(&nm("accum"))?,
        prim: module.load_function(&nm("prim"))?,
        rate: module.load_function(&nm("rate"))?,
        rk: module.load_function(&nm("rk"))?,
        ghost: module.load_function(&nm("ghost"))?,
        amr_ghost: module.load_function(&nm("amr_ghost"))?,
        amr_rgather: module.load_function(&nm("amr_rgather"))?,
        amr_pgather: module.load_function(&nm("amr_pgather"))?,
        amr_ppass: module.load_function(&nm("amr_ppass"))?,
        amr_rpass: module.load_function(&nm("amr_rpass"))?,
        amr_faces: module.load_function(&nm("amr_faces"))?,
        amr_detail: module.load_function(&nm("amr_detail"))?,
        amr_reflux: module.load_function(&nm("amr_reflux"))?,
    };
    Ok((stream, k))
}

pub(crate) fn cfg(len: usize) -> LaunchConfig {
    LaunchConfig { grid_dim: ((len as u64).div_ceil(BLOCK as u64).max(1) as u32, 1, 1), block_dim: (BLOCK, 1, 1), shared_mem_bytes: 0 }
}

pub struct GpuSolver<R: Real> {
    pub sys: GpuSystem,
    pub mesh: Mesh,
    pub geo: Geo,
    pub bc: [[Bc<7>; 2]; 3],
    pub cfl: f64,
    pub tol: f64,
    pub mode: FluxMode,
    pub positivity: bool,
    pub t: f64,
    stream: Arc<CudaStream>,
    k: Kern,
    nv: usize,
    u: CudaSlice<R>,
    u0: CudaSlice<R>,
    w: CudaSlice<R>,
    rhs: CudaSlice<R>,
    f: CudaSlice<R>,
    uf: CudaSlice<R>,
    th: CudaSlice<R>,
    mf: CudaSlice<R>,
    mu: CudaSlice<R>,
    counters: CudaSlice<u64>,
    rates: CudaSlice<R>,
    dir: [[CudaSlice<R>; 2]; 3],
    bl0: CudaSlice<i32>,
}

impl<R: Real> GpuSolver<R> {
    /// Upload the state of a CPU solver (its mesh, boundary conditions and
    /// scheme settings are copied too).
    pub fn from_cpu<const N: usize, S: System<N>>(cpu: &Solver<N, S>, sys: GpuSystem) -> Result<Self> {
        assert_eq!(N, sys.nv(), "system size mismatch");
        let (stream, k) = load_kernels::<R>(sys)?;
        let m = cpu.mesh.clone();
        let mm = [m.m(0), m.m(1), m.m(2)];
        let geo = Geo {
            m: mm.map(|v| v as i32),
            n: m.n.map(|v| v as i32),
            g: [m.g(0), m.g(1), m.g(2)].map(|v| v as i32),
            _pad: 0,
            stride: [(mm[1] * mm[2]) as i64, mm[2] as i64, 1],
            ncell: m.len() as i64,
            ntot: m.len() as i64,
            rtot: m.cells() as i64,
            nb: 1,
            s0: 0,
        };
        let ncell = m.len();
        let nint = m.cells();
        let nv = N;
        // SoA upload
        let mut host = vec![R::default(); nv * ncell];
        for (c, v) in cpu.u.iter().enumerate() {
            for q in 0..nv {
                host[q * ncell + c] = R::from_f64(v[q]);
            }
        }
        let u = stream.clone_htod(&host)?;
        let u0 = stream.alloc_zeros::<R>(nv * ncell)?;
        let w = stream.alloc_zeros::<R>(nv * ncell)?;
        let rhs = stream.alloc_zeros::<R>(nv * nint)?;
        let maxf = (0..3).filter(|&d| m.active(d)).map(|d| (m.n[d] + 1) * m.n[(d + 1) % 3] * m.n[(d + 2) % 3]).max().unwrap_or(1);
        let f = stream.alloc_zeros::<R>(nv * maxf)?;
        let uf = stream.alloc_zeros::<R>(maxf)?;
        let th = stream.alloc_zeros::<R>(maxf)?;
        let maxm = (0..3).filter(|&d| m.active(d)).map(|d| (m.n[d] + 5) * m.n[(d + 1) % 3] * m.n[(d + 2) % 3]).max().unwrap_or(1);
        let mf = stream.alloc_zeros::<R>(nv * maxm)?;
        let mu = stream.alloc_zeros::<R>(maxm)?;
        let counters = stream.alloc_zeros::<u64>(4)?;
        let rates = stream.alloc_zeros::<R>(2 * RATE_BLOCKS as usize)?;
        let mut bc = [[Bc::Periodic; 2]; 3];
        let mut dir: Vec<[CudaSlice<R>; 2]> = Vec::new();
        for d in 0..3 {
            let mut pair = Vec::new();
            for s in 0..2 {
                let (b, vals) = match cpu.bc[d][s] {
                    Bc::Periodic => (Bc::Periodic, vec![R::default(); nv]),
                    Bc::Outflow => (Bc::Outflow, vec![R::default(); nv]),
                    Bc::Reflect => (Bc::Reflect, vec![R::default(); nv]),
                    Bc::Characteristic(v) => {
                        let mut arr = [0.0; 7];
                        arr[..N].copy_from_slice(&v);
                        (Bc::Characteristic(arr), v.iter().map(|x| R::from_f64(*x)).collect())
                    }
                    Bc::Dirichlet(v) => {
                        let mut arr = [0.0; 7];
                        arr[..N].copy_from_slice(&v);
                        (Bc::Dirichlet(arr), v.iter().map(|x| R::from_f64(*x)).collect())
                    }
                };
                bc[d][s] = b;
                pair.push(stream.clone_htod(&vals)?);
            }
            let b1 = pair.pop().unwrap();
            let b0 = pair.pop().unwrap();
            dir.push([b0, b1]);
        }
        let dir: [[CudaSlice<R>; 2]; 3] = [dir.remove(0), dir.remove(0), dir.remove(0)];
        let bl0 = stream.clone_htod(&[0i32])?;
        Ok(GpuSolver {
            sys,
            mesh: m,
            geo,
            bc,
            cfl: cpu.cfl,
            tol: cpu.tol,
            mode: cpu.mode,
            positivity: cpu.positivity,
            t: cpu.t,
            stream,
            k,
            nv,
            u,
            u0,
            w,
            rhs,
            f,
            uf,
            th,
            mf,
            mu,
            counters,
            rates,
            dir,
            bl0,
        })
    }

    fn eos(&self) -> ([R; 8], i32) {
        let (a, kap) = self.sys.eos_args();
        (a.map(R::from_f64), kap)
    }

    fn fill_ghosts(&mut self) -> Result<()> {
        for d in 0..3 {
            if !self.mesh.active(d) {
                continue;
            }
            let code = |b: &Bc<7>| -> i32 {
                match b {
                    Bc::Periodic => 0,
                    Bc::Outflow => 1,
                    Bc::Reflect => 2,
                    Bc::Dirichlet(_) => 3,
                    Bc::Characteristic(_) => 4,
                }
            };
            let (lo, hi) = (code(&self.bc[d][0]), code(&self.bc[d][1]));
            let (e0, e1) = { let mut o: Vec<usize> = (0..3).filter(|&e| e != d).collect(); o.sort(); (o[0], o[1]) };
            let np = self.mesh.m(e0) * self.mesh.m(e1) * self.mesh.g(d);
            let mom = self.sys.mom(d);
            let di = d as i32;
            let geo = self.geo;
            let (e, kap) = self.eos();
            let mut lb = self.stream.launch_builder(&self.k.ghost);
            for v in e.iter() {
                lb.arg(v);
            }
            lb.arg(&kap).arg(&mut self.u).arg(&geo).arg(&di).arg(&lo).arg(&hi).arg(&mom).arg(&self.dir[d][0]).arg(&self.dir[d][1]);
            unsafe { lb.launch(cfg(np))? };
        }
        Ok(())
    }

    fn prim(&mut self) -> Result<()> {
        let (e, kap) = self.eos();
        let geo = self.geo;
        let mut lb = self.stream.launch_builder(&self.k.prim);
        for v in e.iter() {
            lb.arg(v);
        }
        lb.arg(&kap).arg(&self.u).arg(&mut self.w).arg(&geo).arg(&self.bl0);
        unsafe { lb.launch(cfg(self.mesh.len()))? };
        Ok(())
    }

    fn rhs(&mut self, dt: f64) -> Result<()> {
        self.fill_ghosts()?;
        self.prim()?;
        let (e, kap) = self.eos();
        let ndim = (0..3).filter(|&d| self.mesh.active(d)).count() as f64;
        let mode: i32 = match self.mode {
            FluxMode::Hybrid => 0,
            FluxMode::Weno => 1,
            FluxMode::Central => 2,
        };
        let tol = R::from_f64(self.tol);
        let geo = self.geo;
        let mut first: i32 = 1;
        for d in 0..3 {
            if !self.mesh.active(d) {
                continue;
            }
            let dx = self.mesh.dx[d];
            let lambda = R::from_f64(if self.positivity { ndim * dt / dx } else { 0.0 });
            let di = d as i32;
            let nf = (self.mesh.n[d] + 1) * self.mesh.n[(d + 1) % 3] * self.mesh.n[(d + 2) % 3];
            let nm = (self.mesh.n[d] + 5) * self.mesh.n[(d + 1) % 3] * self.mesh.n[(d + 2) % 3];
            {
                let mut lb = self.stream.launch_builder(&self.k.theta);
                for v in e.iter() {
                    lb.arg(v);
                }
                lb.arg(&kap).arg(&self.w).arg(&geo).arg(&self.bl0).arg(&di).arg(&mode).arg(&tol).arg(&mut self.th);
                unsafe { lb.launch(cfg(nf))? };
            }
            if mode != 2 {
                let mut lb = self.stream.launch_builder(&self.k.midflux);
                for v in e.iter() {
                    lb.arg(v);
                }
                lb.arg(&kap).arg(&self.w).arg(&geo).arg(&self.bl0).arg(&di).arg(&self.th).arg(&mut self.mf).arg(&mut self.mu).arg(&mut self.counters);
                unsafe { lb.launch(cfg(nm))? };
            }
            {
                let mut lb = self.stream.launch_builder(&self.k.flux);
                for v in e.iter() {
                    lb.arg(v);
                }
                lb.arg(&kap).arg(&self.w).arg(&geo).arg(&self.bl0).arg(&di).arg(&lambda).arg(&self.th).arg(&self.mf).arg(&self.mu).arg(&mut self.f).arg(&mut self.uf).arg(&mut self.counters);
                unsafe { lb.launch(cfg(nf))? };
            }
            {
                let dxr = R::from_f64(dx);
                let mut lb = self.stream.launch_builder(&self.k.accum);
                for v in e.iter() {
                    lb.arg(v);
                }
                lb.arg(&kap).arg(&self.w).arg(&geo).arg(&self.bl0).arg(&di).arg(&dxr).arg(&self.f).arg(&self.uf).arg(&mut self.rhs).arg(&first);
                unsafe { lb.launch(cfg(self.mesh.cells()))? };
            }
            first = 0;
        }
        Ok(())
    }

    fn rk(&mut self, al: f64, be: f64, dt: f64) -> Result<()> {
        let geo = self.geo;
        let (a, b, t) = (R::from_f64(al), R::from_f64(be), R::from_f64(dt));
        let mut lb = self.stream.launch_builder(&self.k.rk);
        lb.arg(&mut self.u).arg(&self.u0).arg(&self.rhs).arg(&geo).arg(&self.bl0).arg(&a).arg(&b).arg(&t);
        unsafe { lb.launch(cfg(self.mesh.cells()))? };
        Ok(())
    }

    /// One SSP-RK3 step.
    pub fn step(&mut self, dt: f64) -> Result<()> {
        self.stream.memcpy_dtod(&self.u, &mut self.u0)?;
        self.rhs(dt)?;
        self.rk(0.0, 1.0, dt)?;
        self.rhs(dt)?;
        self.rk(0.75, 0.25, dt)?;
        self.rhs(dt)?;
        self.rk(1.0 / 3.0, 2.0 / 3.0, dt)?;
        self.t += dt;
        Ok(())
    }

    /// CFL (and multiphase dilatation) time step.
    pub fn dt_cfl(&mut self) -> Result<f64> {
        self.fill_ghosts()?;
        self.prim()?;
        let (e, kap) = self.eos();
        let geo = self.geo;
        let dx = self.mesh.dx.map(R::from_f64);
        {
            let lc = LaunchConfig { grid_dim: (RATE_BLOCKS, 1, 1), block_dim: (256, 1, 1), shared_mem_bytes: 0 };
            let mut lb = self.stream.launch_builder(&self.k.rate);
            for v in e.iter() {
                lb.arg(v);
            }
            lb.arg(&kap).arg(&self.w).arg(&geo).arg(&self.bl0).arg(&dx[0]).arg(&dx[1]).arg(&dx[2]).arg(&mut self.rates);
            unsafe { lb.launch(lc)? };
        }
        let r = self.stream.clone_dtoh(&self.rates)?;
        let (mut r1, mut r2) = (0.0f64, 0.0f64);
        for b in 0..RATE_BLOCKS as usize {
            r1 = r1.max(r[2 * b].to_f64());
            r2 = r2.max(r[2 * b + 1].to_f64());
        }
        let mut dt = self.cfl / r1;
        if r2 > 0.0 {
            dt = dt.min(0.5 / r2);
        }
        Ok(dt)
    }

    pub fn run_to(&mut self, t_end: f64) -> Result<usize> {
        let mut n = 0;
        while self.t < t_end * (1.0 - 1e-14) {
            let dt = self.dt_cfl()?.min(t_end - self.t);
            if !(dt.is_finite() && dt > 0.0) {
                return Err("time step collapsed".into());
            }
            self.step(dt)?;
            n += 1;
        }
        Ok(n)
    }

    /// Conserved state with ghosts, AoS in f64 (CPU layout).
    pub fn download<const N: usize>(&self) -> Result<Vec<[f64; N]>> {
        let h = self.stream.clone_dtoh(&self.u)?;
        let ncell = self.mesh.len();
        Ok((0..ncell).map(|c| std::array::from_fn(|q| h[q * ncell + c].to_f64())).collect())
    }

    /// [shock-path faces, first-order fallbacks, positivity-limited faces]
    pub fn counters(&self) -> Result<[u64; 3]> {
        let c = self.stream.clone_dtoh(&self.counters)?;
        Ok([c[0], c[1], c[2]])
    }

    pub fn synchronize(&self) -> Result<()> {
        self.stream.synchronize()?;
        Ok(())
    }

    pub fn nv(&self) -> usize {
        self.nv
    }
}
