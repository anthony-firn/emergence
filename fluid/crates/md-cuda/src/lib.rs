//! GPU (CUDA, NVRTC) molecular dynamics for Lennard-Jones pair potentials.
//!
//! The device state is built from a CPU [`md::System`] and stays resident on
//! the GPU across steps: only diagnostics are copied back. Kernels live in
//! `src/kernels.cu` and cover the cell-list construction (bin, scan, fill),
//! Lennard-Jones forces / energies / virial with a truncated-shifted (or
//! shifted-force / bare) cutoff, velocity-Verlet integration, and the
//! diagnostic reductions.
//!
//! Both `f64` and `f32` are instantiated. On a consumer GPU (RTX 3050) the
//! `f64` path runs at 1/64 of the `f32` rate, so `f32` is the useful one
//! for production; the `f64` path exists to check the CPU engine to
//! round-off.

use cudarc::driver::{
    CudaContext, CudaFunction, CudaSlice, CudaStream, DeviceRepr, LaunchConfig,
    PushKernelArg, ValidAsZeroBits,
};
use cudarc::nvrtc::{compile_ptx_with_opts, CompileOptions};
use md::potential::{Cutoff, Potential};
use std::sync::Arc;

pub mod gk;

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
    counts: CudaFunction,
    scan: CudaFunction,
    fill: CudaFunction,
    force: CudaFunction,
    kick: CudaFunction,
    drift: CudaFunction,
    diag: CudaFunction,
    press: CudaFunction,
    scale: CudaFunction,
}

fn load_kernels<R: Real>() -> Result<(Arc<CudaStream>, Kern)> {
    let ctx = CudaContext::new(0)?;
    let stream = ctx.default_stream();
    let ptx = compile_ptx_with_opts(
        include_str!("kernels.cu"),
        CompileOptions { arch: Some("sm_86"), ..Default::default() },
    )?;
    let module = ctx.load_module(ptx)?;
    let nm = |base: &str| format!("md_{base}_{}", R::SUFFIX);
    let k = Kern {
        counts: module.load_function(&nm("cell_counts"))?,
        scan: module.load_function(&nm("cell_scan"))?,
        fill: module.load_function(&nm("cell_fill"))?,
        force: module.load_function(&nm("force"))?,
        kick: module.load_function(&nm("kick"))?,
        drift: module.load_function(&nm("drift"))?,
        diag: module.load_function(&nm("diag"))?,
        press: module.load_function(&nm("press"))?,
        scale: module.load_function(&nm("scale"))?,
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

/// True if a CUDA context can be created (a GPU and driver are present).
///
/// With cudarc's `dynamic-loading` feature a missing driver library makes
/// `CudaContext::new` *panic* rather than return `Err`, so probe for the
/// library first and treat any panic as "no GPU".
pub fn gpu_available() -> bool {
    if !unsafe { cudarc::driver::sys::is_culib_present() } {
        return false;
    }
    std::panic::catch_unwind(|| CudaContext::new(0).is_ok()).unwrap_or(false)
}

fn tail_energy(pot: &Potential, rho: f64, rc: f64) -> f64 {
    pot.tail_energy_per_atom(rho, rc)
}

fn tail_pressure(pot: &Potential, rho: f64, rc: f64) -> f64 {
    pot.tail_pressure(rho, rc)
}

/// A GPU-resident molecular-dynamics system.
pub struct GpuMd<R: Real> {
    stream: Arc<CudaStream>,
    k: Kern,
    n: usize,
    box_l: f64,
    rc: f64,
    sigma: f64,
    epsilon: f64,
    /// 0 = shifted, 1 = shifted-force, 2 = bare (tail-corrected).
    scheme: i32,
    tail: bool,
    nc: [i32; 3],
    ncells: usize,
    pos: CudaSlice<R>,
    vel: CudaSlice<R>,
    force: CudaSlice<R>,
    pe: CudaSlice<R>,
    stress: CudaSlice<R>,
    counts: CudaSlice<i32>,
    starts: CudaSlice<i32>,
    order: CudaSlice<i32>,
    cell_of: CudaSlice<i32>,
    diag: CudaSlice<R>,
    /// Double-precision pressure-tensor samples for Green–Kubo analysis:
    /// `6 * press_cap` values, `press_len` slots filled.
    press_samples: CudaSlice<f64>,
    press_cap: usize,
    press_len: usize,
}

fn scheme_code(c: Cutoff) -> i32 {
    match c {
        Cutoff::Shifted => 0,
        Cutoff::ShiftedForce => 1,
        Cutoff::TailCorrected => 2,
    }
}

impl<R: Real> GpuMd<R> {
    /// Build a system from an explicit CPU configuration.
    pub fn new(
        pos: &[[f64; 3]],
        vel: &[[f64; 3]],
        box_l: f64,
        rc: f64,
        potential: &Potential,
        scheme: Cutoff,
    ) -> Result<Self> {
        let n = pos.len();
        assert!(n > 0, "empty system");
        assert!(box_l > 2.0 * rc, "box {box_l} too small for min-image rc={rc}");
        // Cell width must be at least rc so that all pairs within rc lie in
        // adjacent cells. `skin = 0` because the list is rebuilt each step.
        let nc: [i32; 3] = std::array::from_fn(|_| ((box_l / rc).floor() as i32).max(1));
        let ncells = (nc[0] * nc[1] * nc[2]) as usize;

        let (stream, k) = load_kernels::<R>()?;
        let mut pflat = vec![R::default(); 3 * n];
        let mut vflat = vec![R::default(); 3 * n];
        for (i, (p, v)) in pos.iter().zip(vel).enumerate() {
            for a in 0..3 {
                pflat[3 * i + a] = R::from_f64(p[a]);
                vflat[3 * i + a] = R::from_f64(v[a]);
            }
        }
        let gpos = stream.clone_htod(&pflat)?;
        let gvel = stream.clone_htod(&vflat)?;
        let force = stream.alloc_zeros::<R>(3 * n)?;
        let pe = stream.alloc_zeros::<R>(n)?;
        let stress = stream.alloc_zeros::<R>(6 * n)?;
        let counts = stream.alloc_zeros::<i32>(ncells)?;
        let starts = stream.alloc_zeros::<i32>(ncells)?;
        let order = stream.alloc_zeros::<i32>(n)?;
        let cell_of = stream.alloc_zeros::<i32>(n)?;
        let diag = stream.alloc_zeros::<R>(13)?;
        let press_samples = stream.alloc_zeros::<f64>(0)?;
        let mut sys = Self {
            stream,
            k,
            n,
            box_l,
            rc,
            sigma: potential.sigma(),
            epsilon: potential.epsilon(),
            scheme: scheme_code(scheme),
            tail: scheme == Cutoff::TailCorrected,
            nc,
            ncells,
            pos: gpos,
            vel: gvel,
            force,
            pe,
            stress,
            counts,
            starts,
            order,
            cell_of,
            diag,
            press_samples,
            press_cap: 0,
            press_len: 0,
        };
        sys.compute_forces()?;
        Ok(sys)
    }

    /// Upload a CPU [`md::System`] (positions, velocities and potential
    /// parameters are copied; masses are assumed unity as in `md`).
    pub fn from_md(sys: &md::System) -> Result<Self> {
        // A fresh copy of the neighbour lists is irrelevant on the GPU.
        Self::new(
            &sys.pos,
            &sys.vel,
            sys.box_l,
            sys.rc,
            &sys.potential,
            sys.scheme,
        )
    }

    pub fn n(&self) -> usize {
        self.n
    }

    pub fn box_l(&self) -> f64 {
        self.box_l
    }

    pub fn scheme(&self) -> i32 {
        self.scheme
    }

    /// Build the cell list (bin, scan, fill).
    fn build_cells(&mut self) -> Result<()> {
        let n = self.n;
        let (ncx, ncy, ncz) = (self.nc[0], self.nc[1], self.nc[2]);
        let ncells = self.ncells;
        let invbox = R::from_f64(1.0 / self.box_l);
        self.stream.memset_zeros(&mut self.counts)?;
        {
            let mut lb = self.stream.launch_builder(&self.k.counts);
            lb.arg(&self.pos).arg(&mut self.counts).arg(&mut self.cell_of).arg(&n);
            lb.arg(&ncx).arg(&ncy).arg(&ncz).arg(&invbox);
            unsafe { lb.launch(cfg(n))? };
        }
        {
            let lc = LaunchConfig { grid_dim: (1, 1, 1), block_dim: (256, 1, 1), shared_mem_bytes: 0 };
            let mut lb = self.stream.launch_builder(&self.k.scan);
            lb.arg(&mut self.counts).arg(&mut self.starts).arg(&ncells);
            unsafe { lb.launch(lc)? };
        }
        {
            let mut lb = self.stream.launch_builder(&self.k.fill);
            lb.arg(&self.cell_of).arg(&mut self.starts).arg(&mut self.order).arg(&n);
            unsafe { lb.launch(cfg(n))? };
        }
        Ok(())
    }

    /// Recompute forces, per-atom energies and virials on the device.
    pub fn compute_forces(&mut self) -> Result<()> {
        self.build_cells()?;
        let n = self.n;
        let (ncx, ncy, ncz) = (self.nc[0], self.nc[1], self.nc[2]);
        let boxr = R::from_f64(self.box_l);
        let rc = R::from_f64(self.rc);
        let rc2 = R::from_f64(self.rc * self.rc);
        let sigma = R::from_f64(self.sigma);
        let eps = R::from_f64(self.epsilon);
        let scheme = self.scheme;
        let mut lb = self.stream.launch_builder(&self.k.force);
        lb.arg(&self.pos)
            .arg(&self.starts)
            .arg(&self.counts)
            .arg(&self.order)
            .arg(&mut self.force)
            .arg(&mut self.pe)
            .arg(&mut self.stress)
            .arg(&n)
            .arg(&ncx)
            .arg(&ncy)
            .arg(&ncz)
            .arg(&boxr)
            .arg(&rc)
            .arg(&rc2)
            .arg(&sigma)
            .arg(&eps)
            .arg(&scheme);
        unsafe { lb.launch(cfg(n))? };
        Ok(())
    }

    /// One velocity-Verlet step of the microcanonical dynamics.
    pub fn step(&mut self, dt: f64) -> Result<()> {
        let n = self.n;
        let half = R::from_f64(0.5 * dt);
        {
            let mut lb = self.stream.launch_builder(&self.k.kick);
            lb.arg(&mut self.vel).arg(&self.force).arg(&n).arg(&half);
            unsafe { lb.launch(cfg(n))? };
        }
        let dtr = R::from_f64(dt);
        {
            let mut lb = self.stream.launch_builder(&self.k.drift);
            lb.arg(&mut self.pos).arg(&self.vel).arg(&n).arg(&dtr);
            unsafe { lb.launch(cfg(n))? };
        }
        self.compute_forces()?;
        {
            let mut lb = self.stream.launch_builder(&self.k.kick);
            lb.arg(&mut self.vel).arg(&self.force).arg(&n).arg(&half);
            unsafe { lb.launch(cfg(n))? };
        }
        Ok(())
    }

    /// Run `steps` velocity-Verlet steps.
    pub fn run(&mut self, steps: usize, dt: f64) -> Result<()> {
        for _ in 0..steps {
            self.step(dt)?;
        }
        Ok(())
    }

    /// Rescale all velocities to exactly `target` temperature (`3N` degrees
    /// of freedom). Used for equilibration only; production runs are NVE so
    /// the dynamics are not perturbed. Costs a device→host round trip.
    pub fn rescale_to_temperature(&mut self, target: f64) -> Result<()> {
        let ke = self.kinetic_energy()?; // 0.5 sum v^2
        if ke <= 0.0 {
            return Ok(());
        }
        let want = 1.5 * self.n as f64 * target;
        let lambda = R::from_f64((want / ke).sqrt());
        let n = self.n;
        let mut lb = self.stream.launch_builder(&self.k.scale);
        lb.arg(&mut self.vel).arg(&n).arg(&lambda);
        unsafe { lb.launch(cfg(n))? };
        Ok(())
    }

    /// Allocate room for `max` pressure-tensor samples (each 6 f64). This
    /// keeps the per-step sampling on the device so the hot loop never
    /// synchronises; call [`Self::take_pressure_samples`] at the end.
    pub fn reserve_pressure_samples(&mut self, max: usize) -> Result<()> {
        self.press_samples = self.stream.alloc_zeros::<f64>(6 * max)?;
        self.press_cap = max;
        self.press_len = 0;
        Ok(())
    }

    /// Recompute the instantaneous [xx,yy,zz,xy,xz,yz] pressure-tensor
    /// numerator and append it to the sample buffer. The reduction is done
    /// in `f64` (`atomicAdd` on `double`), so the `f32` force path still
    /// accumulates the virial sum in double precision.
    pub fn sample_pressure_tensor(&mut self) -> Result<()> {
        if self.press_len >= self.press_cap {
            return Err("pressure sample buffer full".into());
        }
        let off = (6 * self.press_len) as i32;
        let n = self.n;
        let mut lb = self.stream.launch_builder(&self.k.press);
        lb.arg(&self.vel)
            .arg(&self.stress)
            .arg(&mut self.press_samples)
            .arg(&n)
            .arg(&off);
        unsafe { lb.launch(cfg(n))? };
        self.press_len += 1;
        Ok(())
    }

    /// Synchronise and return the collected intensive pressure-tensor
    /// samples `[xx,yy,zz,xy,xz,yz]/V` (with the isotropic tail pressure on
    /// the diagonal for the `TailCorrected` scheme). The buffer is zeroed
    /// again so it can be reused.
    pub fn take_pressure_samples(&mut self) -> Result<Vec<[f64; 6]>> {
        self.stream.synchronize()?;
        let all = self.stream.clone_dtoh(&self.press_samples)?;
        let v = self.box_l.powi(3);
        let mut tail = 0.0;
        if self.tail {
            let rho = self.n as f64 / v;
            let pot = Potential::LennardJones(materials::LennardJones::new(self.sigma, self.epsilon));
            tail = tail_pressure(&pot, rho, self.rc);
        }
        let m = self.press_len;
        let out = (0..m)
            .map(|s| {
                let mut p: [f64; 6] = std::array::from_fn(|k| all[6 * s + k] / v);
                p[0] += tail;
                p[1] += tail;
                p[2] += tail;
                p
            })
            .collect();
        self.stream.memset_zeros(&mut self.press_samples)?;
        self.press_len = 0;
        Ok(out)
    }

    /// Number of pressure samples currently buffered.
    pub fn pressure_sample_count(&self) -> usize {
        self.press_len
    }

    /// Recompute the 13 raw diagnostic reductions.
    fn diag_values(&mut self) -> Result<[f64; 13]> {
        let n = self.n;
        self.stream.memset_zeros(&mut self.diag)?;
        {
            let mut lb = self.stream.launch_builder(&self.k.diag);
            lb.arg(&self.vel).arg(&self.pe).arg(&self.stress).arg(&mut self.diag).arg(&n);
            unsafe { lb.launch(cfg(n))? };
        }
        let h = self.stream.clone_dtoh(&self.diag)?;
        Ok(std::array::from_fn(|i| h[i].to_f64()))
    }

    pub fn positions(&self) -> Result<Vec<[f64; 3]>> {
        let h = self.stream.clone_dtoh(&self.pos)?;
        Ok(self.flat_to_vec(&h))
    }

    pub fn velocities(&self) -> Result<Vec<[f64; 3]>> {
        let h = self.stream.clone_dtoh(&self.vel)?;
        Ok(self.flat_to_vec(&h))
    }

    pub fn forces(&self) -> Result<Vec<[f64; 3]>> {
        let h = self.stream.clone_dtoh(&self.force)?;
        Ok(self.flat_to_vec(&h))
    }

    fn flat_to_vec(&self, h: &[R]) -> Vec<[f64; 3]> {
        (0..self.n).map(|i| std::array::from_fn(|a| h[3 * i + a].to_f64())).collect()
    }

    /// Total configurational energy (including the analytic tail correction
    /// for the bare-truncated scheme).
    pub fn potential_energy(&mut self) -> Result<f64> {
        let d = self.diag_values()?;
        let rho = self.n as f64 / self.box_l.powi(3);
        let mut e = d[0];
        if self.tail {
            let pot = Potential::LennardJones(materials::LennardJones::new(self.sigma, self.epsilon));
            e += tail_energy(&pot, rho, self.rc) * self.n as f64;
        }
        Ok(e)
    }

    pub fn kinetic_energy(&mut self) -> Result<f64> {
        let d = self.diag_values()?;
        Ok(0.5 * (d[7] + d[8] + d[9]))
    }

    /// Instantaneous temperature `2K/(3N)` in reduced units.
    pub fn temperature(&mut self) -> Result<f64> {
        Ok(2.0 * self.kinetic_energy()? / (3.0 * self.n as f64))
    }

    pub fn total_energy(&mut self) -> Result<f64> {
        let d = self.diag_values()?;
        let rho = self.n as f64 / self.box_l.powi(3);
        let mut pe = d[0];
        if self.tail {
            let pot = Potential::LennardJones(materials::LennardJones::new(self.sigma, self.epsilon));
            pe += tail_energy(&pot, rho, self.rc) * self.n as f64;
        }
        Ok(pe + 0.5 * (d[7] + d[8] + d[9]))
    }

    /// Instantaneous pressure tensor `[xx, yy, zz, xy, xz, yz]`, reduced
    /// units, matching `md::System::pressure_tensor`.
    pub fn pressure_tensor(&mut self) -> Result<[f64; 6]> {
        let d = self.diag_values()?;
        let v = self.box_l.powi(3);
        // kinetic tensor (indices 7..12) + half the pair-virial sum.
        let mut p = [0.0f64; 6];
        for k in 0..6 {
            p[k] = d[7 + k] + 0.5 * d[1 + k];
        }
        if self.tail {
            let rho = self.n as f64 / v;
            let pot = Potential::LennardJones(materials::LennardJones::new(self.sigma, self.epsilon));
            let t = tail_pressure(&pot, rho, self.rc) * v;
            p[0] += t;
            p[1] += t;
            p[2] += t;
        }
        for x in p.iter_mut() {
            *x /= v;
        }
        Ok(p)
    }

    /// Instantaneous pressure tensor `[xx,yy,zz,xy,xz,yz]` in reduced units,
    /// with the reduction accumulated in `f64` on the device (so the `f32`
    /// path keeps double-precision pressure statistics). Matches
    /// `md::System::pressure_tensor` including the isotropic tail pressure.
    pub fn pressure_tensor_f64(&mut self) -> Result<[f64; 6]> {
        let n = self.n;
        let mut one = self.stream.alloc_zeros::<f64>(6)?;
        {
            let mut lb = self.stream.launch_builder(&self.k.press);
            lb.arg(&self.vel).arg(&self.stress).arg(&mut one).arg(&n).arg(&0i32);
            unsafe { lb.launch(cfg(n))? };
        }
        let h = self.stream.clone_dtoh(&one)?;
        let v = self.box_l.powi(3);
        let mut p: [f64; 6] = std::array::from_fn(|k| h[k] / v);
        if self.tail {
            let rho = self.n as f64 / v;
            let pot = Potential::LennardJones(materials::LennardJones::new(self.sigma, self.epsilon));
            let t = tail_pressure(&pot, rho, self.rc);
            p[0] += t;
            p[1] += t;
            p[2] += t;
        }
        Ok(p)
    }

    /// Scalar pressure `tr(P)/3`.
    pub fn pressure(&mut self) -> Result<f64> {
        let p = self.pressure_tensor()?;
        Ok((p[0] + p[1] + p[2]) / 3.0)
    }

    /// Internal pair virial `0.5 sum_i tr(stress_i)`, matching
    /// `md::System::internal_pair_virial`.
    pub fn internal_pair_virial(&mut self) -> Result<f64> {
        let d = self.diag_values()?;
        Ok(0.5 * (d[1] + d[2] + d[3]))
    }

    pub fn synchronize(&self) -> Result<()> {
        self.stream.synchronize()?;
        Ok(())
    }

    /// Diagnostic: per-cell occupancy after a rebuild, length `ncells`.
    /// Integer and exact, so it can be compared to a CPU binning.
    pub fn cell_counts(&mut self) -> Result<Vec<i32>> {
        self.build_cells()?;
        self.stream.synchronize()?;
        Ok(self.stream.clone_dtoh(&self.counts)?)
    }

    /// Diagnostic: cell-grid size, maximum occupancy and mean occupancy
    /// after a rebuild. Used by tests/examples, not on the hot path.
    pub fn cell_stats(&mut self) -> Result<(usize, i32, f64)> {
        self.build_cells()?;
        self.stream.synchronize()?;
        let h = self.stream.clone_dtoh(&self.counts)?;
        let max = h.iter().copied().max().unwrap_or(0);
        let sum: i64 = h.iter().map(|&x| x as i64).sum();
        Ok((self.ncells, max, sum as f64 / self.ncells as f64))
    }
}
