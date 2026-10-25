//! CUDA backend for the implicit renderer.
//!
//! The CAD field is **compiled**, not interpreted: [`codegen`] emits a
//! `__device__ double scene_sdf(double3)` whose arithmetic is the same
//! expressions, in the same order, as `cad::solid`'s, and nvrtc builds it
//! together with a fixed sphere-tracing kernel. So the device and host
//! fields are the same function, which the tests check over tens of
//! thousands of points, and the GPU image is not an approximation of the
//! CPU one — it is the same image.
//!
//! ```text
//!   cad::Solid --codegen--> CUDA source --nvrtc--> module
//!                                                    |
//!   camera ------------------------------------------>| trace kernel
//!                                                    v
//!                              depth, normal and step buffers
//! ```
//!
//! # Precision, and what it does and does not buy
//!
//! `Double` is the default because it is the mode in which the device and
//! host fields agree *to the last bit*: any disagreement is then a bug in
//! the generated code, where in `float` it could be either a bug or
//! rounding and the test would have to tolerate it and so would stop
//! catching bugs.
//!
//! `Single` is available and its deviation is measured, not assumed:
//! 4.4e-8 absolute on a field of order one, a few units in the last place.
//! It is only **1.3x faster**, though, on a card whose `float` throughput
//! is thirty-two times its `double` throughput. That is worth stating
//! plainly: a sphere-tracing kernel is not floating-point-throughput
//! bound. Neighbouring pixels take wildly different numbers of steps, so
//! a warp waits for its slowest lane, and the cost is branch divergence
//! and latency rather than arithmetic. Precision is the wrong knob;
//! per-warp work balancing would be the right one, and is not done here.
//!
//! # Runtime requirements
//!
//! `libnvrtc` and the CUDA driver, loaded dynamically:
//! `source fluid/scripts/cuda-env.sh`. Every entry point returns a
//! `Result`, and [`available`] says whether a device can be opened at all,
//! so a machine with no GPU gets a clear answer rather than a panic.

use std::sync::Arc;

use cad::solid::Solid;
use cudarc::driver::{CudaContext, CudaFunction, CudaStream, LaunchConfig, PushKernelArg};
use cudarc::nvrtc::{compile_ptx_with_opts, CompileOptions};
use geometry::P3;
use viz::camera::Camera;

pub mod codegen;

pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

const BLOCK: u32 = 128;

/// Whether a CUDA device can be opened.
pub fn available() -> bool {
    CudaContext::new(0).is_ok()
}

/// A compiled scene: one solid, on the device.
pub struct GpuTracer {
    _ctx: Arc<CudaContext>,
    stream: Arc<CudaStream>,
    trace: CudaFunction,
    sample: CudaFunction,
    bbox: (P3, P3),
    bounded: bool,
    /// Device functions the generator emitted, for reporting.
    pub functions: usize,
    /// The precision the field is evaluated in.
    pub precision: codegen::Precision,
    /// The generated source, kept so a caller can inspect or dump it.
    pub source: String,
}

/// What a GPU trace produced.
#[derive(Clone, Debug)]
pub struct GpuFrame {
    pub width: usize,
    pub height: usize,
    /// Distance to the first hit; `+inf` for a miss and **negative** for a
    /// march that ran out of steps, carrying how far it had gone.
    pub depth: Vec<f64>,
    /// Unit normals, three per pixel.
    pub normals: Vec<f64>,
    pub steps: Vec<i32>,
}

impl GpuFrame {
    pub fn hits(&self) -> usize {
        self.depth.iter().filter(|d| d.is_finite() && **d > 0.0).count()
    }

    pub fn exhausted(&self) -> usize {
        self.depth.iter().filter(|d| **d < 0.0).count()
    }

    pub fn total_steps(&self) -> u64 {
        self.steps.iter().map(|&s| s as u64).sum()
    }
}

impl GpuTracer {
    /// Compile a solid for the device in double precision: the mode whose
    /// field is checked bit for bit against the host.
    pub fn new(solid: &Solid) -> Result<Self> {
        Self::new_with(solid, codegen::Precision::Double)
    }

    /// Compile a solid at a chosen precision. `Single` is several times
    /// faster on a consumer card; its deviation from the exact path is
    /// measured in `tests/gpu.rs`, not assumed.
    pub fn new_with(solid: &Solid, precision: codegen::Precision) -> Result<Self> {
        let generated = codegen::emit_with(solid, precision)?;
        let source = format!(
            "{}\n{}",
            generated.source,
            include_str!("kernels.cu")
        );
        let ctx = CudaContext::new(0)?;
        let stream = ctx.default_stream();
        let ptx = compile_ptx_with_opts(
            &source,
            CompileOptions {
                arch: Some("sm_86"),
                ..Default::default()
            },
        )?;
        let module = ctx.load_module(ptx)?;
        let bbox = solid.bbox();
        let bounded = (0..3).all(|k| bbox.0[k].is_finite() && bbox.1[k].is_finite());
        Ok(Self {
            trace: module.load_function("trace")?,
            sample: module.load_function("sample_field")?,
            _ctx: ctx,
            stream,
            bbox,
            bounded,
            functions: generated.functions,
            precision,
            source,
        })
    }

    /// Evaluate the device field at a list of points.
    ///
    /// This is the test of the code generator: the values must equal
    /// `Solid::sdf` at the same points.
    pub fn sample(&self, points: &[P3]) -> Result<Vec<f64>> {
        let n = points.len();
        if n == 0 {
            return Ok(Vec::new());
        }
        let flat: Vec<f64> = points.iter().flat_map(|p| p.iter().copied()).collect();
        let d_points = self.stream.clone_htod(&flat)?;
        let mut d_out = self.stream.alloc_zeros::<f64>(n)?;
        let count = n as i32;
        unsafe {
            self.stream
                .launch_builder(&self.sample)
                .arg(&d_points)
                .arg(&mut d_out)
                .arg(&count)
                .launch(cfg(n))?;
        }
        Ok(self.stream.clone_dtoh(&d_out)?)
    }

    /// Sphere-trace a frame.
    pub fn trace(
        &self,
        camera: &Camera,
        width: usize,
        height: usize,
        options: &viz::MarchOptions,
    ) -> Result<GpuFrame> {
        let (w, h) = (width.max(1), height.max(1));
        let n = w * h;
        let frame = camera.frame();
        let aspect = w as f64 / h as f64;
        let half_h = (0.5 * camera.fov_y).tan();
        let t_max = if options.t_max.is_finite() { options.t_max } else { 1e30 };

        let mut d_depth = self.stream.alloc_zeros::<f64>(n)?;
        let mut d_normals = self.stream.alloc_zeros::<f64>(3 * n)?;
        let mut d_steps = self.stream.alloc_zeros::<i32>(n)?;

        let eye = d3(camera.eye);
        let right = d3(frame.right);
        let up = d3(frame.up);
        let forward = d3(frame.forward);
        let lo = d3(self.bbox.0);
        let hi = d3(self.bbox.1);
        let bounded = i32::from(self.bounded);
        let (wi, hi_i) = (w as i32, h as i32);
        let max_steps = options.max_steps as i32;
        unsafe {
            self.stream
                .launch_builder(&self.trace)
                .arg(&eye)
                .arg(&right)
                .arg(&up)
                .arg(&forward)
                .arg(&half_h)
                .arg(&aspect)
                .arg(&wi)
                .arg(&hi_i)
                .arg(&lo)
                .arg(&hi)
                .arg(&bounded)
                .arg(&options.tolerance)
                .arg(&max_steps)
                .arg(&t_max)
                .arg(&options.normal_epsilon)
                .arg(&mut d_depth)
                .arg(&mut d_normals)
                .arg(&mut d_steps)
                .launch(cfg(n))?;
        }
        Ok(GpuFrame {
            width: w,
            height: h,
            depth: self.stream.clone_dtoh(&d_depth)?,
            normals: self.stream.clone_dtoh(&d_normals)?,
            steps: self.stream.clone_dtoh(&d_steps)?,
        })
    }
}

fn cfg(len: usize) -> LaunchConfig {
    LaunchConfig {
        grid_dim: (((len as u32) + BLOCK - 1) / BLOCK, 1, 1),
        block_dim: (BLOCK, 1, 1),
        shared_mem_bytes: 0,
    }
}

/// `double3` in the layout the driver expects.
#[repr(C)]
#[derive(Clone, Copy)]
struct Double3 {
    x: f64,
    y: f64,
    z: f64,
}

unsafe impl cudarc::driver::DeviceRepr for Double3 {}

fn d3(p: P3) -> Double3 {
    Double3 { x: p[0], y: p[1], z: p[2] }
}

/// As [`cpu_frame`], across all cores: the fair comparison for a
/// throughput claim.
///
/// The single-threaded version is the correctness reference, because its
/// result does not depend on how the work was divided; this one is the
/// performance reference.
pub fn cpu_frame_parallel(
    solid: &Solid,
    camera: &Camera,
    width: usize,
    height: usize,
    options: &viz::MarchOptions,
) -> GpuFrame {
    use radiation::math::Rng;
    use rayon::prelude::*;
    let (w, h) = (width.max(1), height.max(1));
    let aspect = w as f64 / h as f64;
    let rows: Vec<(Vec<f64>, Vec<f64>, Vec<i32>)> = (0..h)
        .into_par_iter()
        .map(|y| {
            let mut depth = vec![0.0; w];
            let mut normals = vec![0.0; 3 * w];
            let mut steps = vec![0i32; w];
            for x in 0..w {
                let mut rng = Rng(1);
                let [u, v] = camera.pixel(x, y, w, h, [0.5, 0.5]);
                let ray = camera.ray(u, v, aspect, &mut rng);
                match viz::trace::march(solid, ray, options) {
                    viz::Trace::Hit(hit) => {
                        depth[x] = hit.t;
                        steps[x] = hit.steps as i32;
                        for k in 0..3 {
                            normals[3 * x + k] = hit.normal[k];
                        }
                    }
                    viz::Trace::Miss => depth[x] = f64::INFINITY,
                    viz::Trace::Exhausted { t, steps: s, .. } => {
                        depth[x] = -t;
                        steps[x] = s as i32;
                    }
                }
            }
            (depth, normals, steps)
        })
        .collect();
    let mut frame = GpuFrame {
        width: w,
        height: h,
        depth: Vec::with_capacity(w * h),
        normals: Vec::with_capacity(3 * w * h),
        steps: Vec::with_capacity(w * h),
    };
    for (d, n, s) in rows {
        frame.depth.extend(d);
        frame.normals.extend(n);
        frame.steps.extend(s);
    }
    frame
}

/// The same trace on the CPU, for comparison: the reference the GPU is
/// checked against.
pub fn cpu_frame(
    solid: &Solid,
    camera: &Camera,
    width: usize,
    height: usize,
    options: &viz::MarchOptions,
) -> GpuFrame {
    use radiation::math::{Rng, Vec3};
    let (w, h) = (width.max(1), height.max(1));
    let aspect = w as f64 / h as f64;
    let mut depth = vec![0.0; w * h];
    let mut normals = vec![0.0; 3 * w * h];
    let mut steps = vec![0i32; w * h];
    for y in 0..h {
        for x in 0..w {
            let i = y * w + x;
            let mut rng = Rng(1);
            let [u, v] = camera.pixel(x, y, w, h, [0.5, 0.5]);
            let ray = camera.ray(u, v, aspect, &mut rng);
            let _ = Vec3(ray.dir.0);
            match viz::trace::march(solid, ray, options) {
                viz::Trace::Hit(hit) => {
                    depth[i] = hit.t;
                    steps[i] = hit.steps as i32;
                    for k in 0..3 {
                        normals[3 * i + k] = hit.normal[k];
                    }
                }
                viz::Trace::Miss => {
                    depth[i] = f64::INFINITY;
                }
                viz::Trace::Exhausted { t, steps: s, .. } => {
                    depth[i] = -t;
                    steps[i] = s as i32;
                }
            }
        }
    }
    GpuFrame { width: w, height: h, depth, normals, steps }
}
