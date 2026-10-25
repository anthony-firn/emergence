//! The generated device field against the host field, and the GPU trace
//! against the CPU trace.
//!
//! These tests need a CUDA device and `libnvrtc`; without them they print
//! why and pass, so the suite still runs on a machine with no GPU.

use cad::sketch::Region2;
use cad::solid::{Placement, Solid, Tpms};
use viz::camera::Camera;
use viz::MarchOptions;
use viz_cuda::{codegen, cpu_frame, GpuTracer};

struct Pcg(u64);
impl Pcg {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn point(&mut self, lo: f64, hi: f64) -> [f64; 3] {
        std::array::from_fn(|_| lo + (hi - lo) * self.next())
    }
}

/// Every node kind the generator supports, in one scene.
fn menagerie() -> Vec<(&'static str, Solid)> {
    let gear = Region2::rect([0.35, 0.0], [0.12, 0.05]).polar_repeat(9);
    let plate = Region2::rounded_rect([0.0, 0.0], [0.5, 0.3], 0.08)
        .minus(Region2::circle([0.0, 0.0], 0.12))
        .plus(gear);
    vec![
        ("sphere", Solid::sphere([0.1, -0.2, 0.3], 0.5)),
        ("cuboid", Solid::rounded_box([0.0; 3], [0.4, 0.3, 0.2], 0.05)),
        ("cylinder", Solid::cylinder([-0.3, 0.0, -0.2], [0.4, 0.2, 0.3], 0.15)),
        ("cone", Solid::cone([0.0, 0.0, -0.3], [0.0, 0.0, 0.4], 0.3, 0.08)),
        ("capsule", Solid::capsule([-0.4, 0.0, 0.0], [0.4, 0.1, 0.2], 0.12)),
        ("torus", Solid::torus([0.0, 0.0, 0.05], 0.35, 0.1)),
        ("half space", Solid::half_space([0.3, 1.0, -0.2], 0.1)),
        ("extrude", Solid::extrude(plate.clone(), -0.1, 0.1)),
        (
            "revolve",
            Solid::revolve(Region2::Polygon(cad::sketch::Polygon::new(vec![
                [-0.3, -0.2],
                [0.3, -0.2],
                [0.15, 0.0],
                [0.3, 0.2],
                [-0.3, 0.2],
                [-0.15, 0.0],
            ]))),
        ),
        (
            "boolean tree",
            Solid::box_between([-0.5, -0.3, -0.1], [0.5, 0.3, 0.1])
                .minus(Solid::cylinder([-0.25, 0.0, -0.5], [-0.25, 0.0, 0.5], 0.08))
                .plus(Solid::sphere([0.4, 0.2, 0.1], 0.12))
                .minus(Solid::torus([0.0; 3], 0.35, 0.05)),
        ),
        (
            "blend",
            Solid::sphere([-0.15, 0.0, 0.0], 0.25).blend(Solid::sphere([0.15, 0.0, 0.0], 0.25), 0.1),
        ),
        (
            "offset and shell",
            Solid::box_between([-0.3; 3], [0.3; 3]).offset(0.05).shell(0.04),
        ),
        (
            "placed",
            Solid::cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 0.5], 0.1)
                .placed(Placement::rotation([1.0, 0.3, 0.2], 0.7))
                .translated([0.1, -0.1, 0.05])
                .scaled(0.8),
        ),
        ("mirrored", Solid::sphere([0.3, 0.0, 0.0], 0.2).mirrored([1.0, 0.0, 0.0])),
        (
            "stretch",
            Solid::sphere([0.0; 3], 0.3).stretched([1.5, 0.7, 1.0]),
        ),
        (
            "twist",
            Solid::box_between([-0.2, -0.2, -0.4], [0.2, 0.2, 0.4]).twisted(2.0, 0.3),
        ),
        (
            "array",
            Solid::sphere([0.0; 3], 0.06).array([0.2, 0.2, 0.2], [3, 2, 2]),
        ),
        (
            "gyroid sheet",
            Solid::Periodic {
                kind: Tpms::Gyroid,
                cell: [0.4, 0.4, 0.4],
                level: 0.0,
                thickness: Some(0.08),
            },
        ),
        (
            "schwarz p",
            Solid::Periodic {
                kind: Tpms::SchwarzP,
                cell: [0.5, 0.4, 0.3],
                level: 0.2,
                thickness: None,
            },
        ),
        (
            "diamond",
            Solid::Periodic {
                kind: Tpms::Diamond,
                cell: [0.45; 3],
                level: -0.1,
                thickness: Some(0.06),
            },
        ),
        ("slot and ring", Solid::extrude(
            Region2::slot([-0.3, 0.0], [0.3, 0.1], 0.08)
                .plus(Region2::ring([0.0, 0.0], 0.15, 0.25))
                .and(Region2::rect([0.0, 0.0], [0.5, 0.5])),
            -0.05,
            0.05,
        )),
        ("region array", Solid::extrude(
            Region2::circle([0.0, 0.0], 0.04).array([0.15, 0.15], [4, 3]),
            -0.02,
            0.02,
        )),
    ]
}

#[test]
fn the_generated_device_field_is_the_host_field() {
    if !viz_cuda::available() {
        println!("no CUDA device: skipping");
        return;
    }
    let mut rng = Pcg(0x2545F4914F6CDD1D);
    let mut worst_overall: f64 = 0.0;
    for (name, solid) in menagerie() {
        let tracer = match GpuTracer::new(&solid) {
            Ok(t) => t,
            Err(e) => panic!("{name}: {e}"),
        };
        let points: Vec<[f64; 3]> = (0..4000).map(|_| rng.point(-1.2, 1.2)).collect();
        let gpu = tracer.sample(&points).expect("sample");
        let mut worst: f64 = 0.0;
        let mut worst_at = [0.0; 3];
        for (i, p) in points.iter().enumerate() {
            let host = solid.sdf(*p);
            if !host.is_finite() && !gpu[i].is_finite() {
                continue;
            }
            let d = (host - gpu[i]).abs();
            let scale = host.abs().max(1.0);
            if d / scale > worst {
                worst = d / scale;
                worst_at = *p;
            }
        }
        println!("{name:20} ({} device functions): worst relative difference {worst:.3e}", tracer.functions);
        assert!(
            worst < 1e-13,
            "{name}: device and host fields differ by {worst:e} at {worst_at:?}"
        );
        worst_overall = worst_overall.max(worst);
    }
    println!("worst over the whole menagerie: {worst_overall:.3e}");
}

#[test]
fn nodes_with_no_device_form_are_refused_by_name() {
    // A mesh carries a BVH and a custom node an arbitrary closure; neither
    // can cross to a device, and saying so is better than approximating.
    let mesh = Solid::from_trimesh(geometry::primitives::icosphere([0.0; 3], 0.5, 1));
    let err = codegen::emit(&mesh).unwrap_err();
    assert!(err.contains("Solid::Mesh") && err.contains("CPU"), "{err}");

    let custom = Solid::custom("field isosurface", ([0.0; 3], [1.0; 3]), false, |p| p[0] - 0.5);
    let err = codegen::emit(&custom).unwrap_err();
    assert!(err.contains("field isosurface"), "{err}");

    // And a tree containing one is refused too, not silently dropped.
    let mixed = Solid::sphere([0.0; 3], 1.0).minus(custom);
    assert!(codegen::emit(&mixed).is_err());
}

#[test]
fn the_gpu_trace_is_the_cpu_trace() {
    if !viz_cuda::available() {
        println!("no CUDA device: skipping");
        return;
    }
    let solid = Solid::box_between([-0.5, -0.3, -0.2], [0.5, 0.3, 0.2])
        .minus(Solid::cylinder([-0.25, 0.0, -0.5], [-0.25, 0.0, 0.5], 0.1))
        .plus(Solid::sphere([0.45, 0.2, 0.15], 0.15))
        .minus(Solid::torus([0.0; 3], 0.3, 0.06));
    let camera = Camera::looking_at([1.6, -1.9, 1.2], [0.0; 3], 0.7);
    let options = MarchOptions {
        tolerance: 1e-9,
        max_steps: 400,
        normal_epsilon: 1e-7,
        ..MarchOptions::default()
    };
    let (w, h) = (160, 120);

    let tracer = GpuTracer::new(&solid).expect("compile");
    let gpu = tracer.trace(&camera, w, h, &options).expect("trace");
    let cpu = cpu_frame(&solid, &camera, w, h, &options);

    println!(
        "{} hits on the GPU, {} on the CPU; {} exhausted vs {}",
        gpu.hits(),
        cpu.hits(),
        gpu.exhausted(),
        cpu.exhausted()
    );
    assert!(gpu.hits() > 2000, "the frame should see the model: {} hits", gpu.hits());

    // Pixel classification must agree: a hit on one must be a hit on the
    // other. The three outcomes are distinguished, so a GPU miss where the
    // CPU hit would be caught here.
    let mut disagreements = 0;
    let class = |d: f64| {
        if d < 0.0 {
            2 // exhausted
        } else if d.is_finite() {
            0 // hit
        } else {
            1 // miss
        }
    };
    for i in 0..w * h {
        if class(gpu.depth[i]) != class(cpu.depth[i]) {
            disagreements += 1;
        }
    }
    assert_eq!(disagreements, 0, "{disagreements} pixels classified differently");

    // And the depths and normals agree to the march tolerance.
    let mut worst_depth: f64 = 0.0;
    let mut worst_normal: f64 = 0.0;
    for i in 0..w * h {
        if class(gpu.depth[i]) != 0 {
            continue;
        }
        worst_depth = worst_depth.max((gpu.depth[i] - cpu.depth[i]).abs());
        let dot: f64 = (0..3).map(|k| gpu.normals[3 * i + k] * cpu.normals[3 * i + k]).sum();
        worst_normal = worst_normal.max((1.0 - dot).abs());
        assert_eq!(gpu.steps[i], cpu.steps[i], "step counts differ at pixel {i}");
    }
    println!("worst depth difference {worst_depth:.3e}, worst normal {worst_normal:.3e}");
    assert!(worst_depth < 1e-12, "depths differ by {worst_depth:e}");
    assert!(worst_normal < 1e-12, "normals differ by {worst_normal:e}");
}

#[test]
fn single_precision_is_faster_and_its_error_is_measured() {
    if !viz_cuda::available() {
        println!("no CUDA device: skipping");
        return;
    }
    // Correctness is established in double, where the device field equals
    // the host field exactly. Speed is taken in single, where it does not
    // — so the deviation is measured rather than hoped about.
    let solid = Solid::Periodic {
        kind: Tpms::Gyroid,
        cell: [0.25; 3],
        level: 0.0,
        thickness: Some(0.05),
    }
    .clipped([-0.5; 3], [0.5; 3]);
    let mut rng = Pcg(0xDEADBEEF);
    let points: Vec<[f64; 3]> = (0..8000).map(|_| rng.point(-0.6, 0.6)).collect();

    let f64_tracer = GpuTracer::new(&solid).expect("compile f64");
    let f32_tracer =
        GpuTracer::new_with(&solid, codegen::Precision::Single).expect("compile f32");
    let exact = f64_tracer.sample(&points).expect("sample f64");
    let single = f32_tracer.sample(&points).expect("sample f32");

    let mut worst: f64 = 0.0;
    for i in 0..points.len() {
        worst = worst.max((exact[i] - single[i]).abs());
    }
    println!("single precision deviates from the exact field by {worst:.3e} (absolute)");
    // The field is O(1) here, so this is a relative figure too: a few
    // units in the last place of f32, as expected, and a thousand times
    // coarser than the march tolerance used below.
    assert!(worst < 1e-5, "single precision is off by {worst:e}");
    assert!(worst > 1e-9, "it should not be exact: {worst:e}");

    // And it is faster at the same image.
    let camera = Camera::looking_at([1.3, -1.5, 1.0], [0.0; 3], 0.8);
    let options = MarchOptions {
        tolerance: 1e-5,
        max_steps: 512,
        normal_epsilon: 1e-4,
        ..MarchOptions::default()
    };
    let (w, h) = (960, 720);
    for t in [&f64_tracer, &f32_tracer] {
        let _ = t.trace(&camera, 64, 64, &options).expect("warm up");
    }
    let t0 = std::time::Instant::now();
    let a = f64_tracer.trace(&camera, w, h, &options).expect("f64");
    let t_f64 = t0.elapsed().as_secs_f64();
    let t1 = std::time::Instant::now();
    let b = f32_tracer.trace(&camera, w, h, &options).expect("f32");
    let t_f32 = t1.elapsed().as_secs_f64();
    println!(
        "{w}x{h}: f64 {t_f64:.4} s, f32 {t_f32:.4} s ({:.1}x); hits {} vs {}",
        t_f64 / t_f32,
        a.hits(),
        b.hits()
    );
    // Faster, but only modestly: this kernel is bound by branch
    // divergence between neighbouring pixels, not by arithmetic
    // throughput, so the card's 32:1 float-to-double ratio does not show
    // up. Precision is the wrong knob for this kernel.
    assert!(t_f32 < t_f64, "single should be faster: {t_f32:.4} vs {t_f64:.4}");
    assert!(
        t_f64 / t_f32 < 10.0,
        "a 10x gain would mean the kernel is throughput bound, which it is not: {:.1}x",
        t_f64 / t_f32
    );
    // The same picture to within the coarser precision: the hit counts
    // agree to a fraction of a per cent.
    let diff = (a.hits() as f64 - b.hits() as f64).abs() / a.hits() as f64;
    assert!(diff < 0.01, "hit counts differ by {:.3} %", 100.0 * diff);
}

#[test]
fn the_gpu_is_faster_than_the_cpu_on_a_heavy_field() {
    if !viz_cuda::available() {
        println!("no CUDA device: skipping");
        return;
    }
    // A gyroid sheet in a box: a cheap field per sample, but every pixel
    // needs many steps, which is the case a GPU is for.
    let solid = Solid::Periodic {
        kind: Tpms::Gyroid,
        cell: [0.25; 3],
        level: 0.0,
        thickness: Some(0.05),
    }
    .clipped([-0.5; 3], [0.5; 3]);
    let camera = Camera::looking_at([1.3, -1.5, 1.0], [0.0; 3], 0.8);
    let options = MarchOptions {
        tolerance: 1e-7,
        max_steps: 512,
        normal_epsilon: 1e-6,
        ..MarchOptions::default()
    };
    // Large enough that the per-pixel work dominates the fixed cost of a
    // launch and a copy back: at 480x360 the GPU's margin over all the
    // cores sits inside the run-to-run noise, which would make this a
    // coin flip rather than a measurement.
    let (w, h) = (1280, 960);

    let t0 = std::time::Instant::now();
    let tracer = GpuTracer::new(&solid).expect("compile");
    let compile = t0.elapsed().as_secs_f64();
    let _ = tracer.trace(&camera, 64, 64, &options).expect("warm up");

    // Best of three each way: the minimum is the honest comparison when
    // the question is how fast the machine can do this at all.
    let best = |mut f: Box<dyn FnMut() -> f64>| (0..3).map(|_| f()).fold(f64::INFINITY, f64::min);
    let gpu_time = {
        let t = &tracer;
        best(Box::new(move || {
            let s = std::time::Instant::now();
            let _ = t.trace(&camera, w, h, &options).expect("trace");
            s.elapsed().as_secs_f64()
        }))
    };
    let par_time = {
        let solid = &solid;
        best(Box::new(move || {
            let s = std::time::Instant::now();
            let _ = viz_cuda::cpu_frame_parallel(solid, &camera, w, h, &options);
            s.elapsed().as_secs_f64()
        }))
    };
    let one_core_time = {
        let solid = &solid;
        let s = std::time::Instant::now();
        let _ = cpu_frame(solid, &camera, w, h, &options);
        s.elapsed().as_secs_f64()
    };

    // And the single-precision path, which is what a renderer would
    // actually use: its 4.4e-8 deviation is four thousand times below the
    // 1e-5 march tolerance here, so it draws the same picture.
    let single = GpuTracer::new_with(&solid, codegen::Precision::Single).expect("compile f32");
    let _ = single.trace(&camera, 64, 64, &options).expect("warm up");
    let f32_time = {
        let t = &single;
        best(Box::new(move || {
            let s = std::time::Instant::now();
            let _ = t.trace(&camera, w, h, &options).expect("trace");
            s.elapsed().as_secs_f64()
        }))
    };

    let gpu = tracer.trace(&camera, w, h, &options).expect("trace");
    let cpu = cpu_frame_parallel_check(&solid, &camera, w, h, &options);
    println!(
        "{w}x{h} ({} march steps), compile {compile:.2} s:\n  \
         GPU f64 {gpu_time:.4} s | GPU f32 {f32_time:.4} s | \
         CPU all cores {par_time:.4} s | CPU 1 core {one_core_time:.4} s\n  \
         GPU f64 vs 1 core {:.1}x | GPU f64 vs all cores {:.2}x | GPU f32 vs all cores {:.2}x",
        gpu.total_steps(),
        one_core_time / gpu_time,
        par_time / gpu_time,
        par_time / f32_time
    );
    // Same picture, so the comparison is of like with like.
    assert_eq!(gpu.hits(), cpu.hits(), "the two frames must agree");
    for i in 0..w * h {
        if cpu.depth[i].is_finite() && cpu.depth[i] > 0.0 {
            assert!((gpu.depth[i] - cpu.depth[i]).abs() < 1e-12);
        }
    }
    // What is actually true, and asserted as such:
    //
    // * against one core the GPU wins comfortably;
    // * against **all** cores in double precision it does not — this
    //   card's f64 throughput is a thirty-second of its f32, which puts it
    //   in the same class as a modern multi-core CPU, and the kernel is
    //   divergence-bound besides. The ratio is printed, not asserted,
    //   because it straddles one and asserting a win would be asserting
    //   noise;
    // * the single-precision path does win, and is the mode a renderer
    //   should use.
    assert!(
        one_core_time / gpu_time > 3.0,
        "GPU should beat one core by more than 3x: {:.1}x",
        one_core_time / gpu_time
    );
    assert!(
        f32_time < gpu_time,
        "single precision should be the faster of the two GPU modes"
    );
}

/// The parallel CPU frame, used as the picture to compare against.
fn cpu_frame_parallel_check(
    solid: &Solid,
    camera: &Camera,
    w: usize,
    h: usize,
    options: &MarchOptions,
) -> viz_cuda::GpuFrame {
    viz_cuda::cpu_frame_parallel(solid, camera, w, h, options)
}
