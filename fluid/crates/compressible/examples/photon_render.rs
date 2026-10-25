//! CPU spectral render of an evolved water/air droplet with a trapped bubble.
//! photon_render [output_directory] [grid=40] [spp=256] [width=640] [steps=12] [fixed|exact|adaptive]
use compressible::{
    capillary::SurfaceTension, elements::fit_inputs, eos::Nasg, multiphase::FiveEq,
    phase_change::v_pt, substances::by_name, Bc, Mesh, Solver, System,
};
use radiation::{
    geometry::{isosurface, quad, Bvh, Primitive, Shape},
    math::Vec3,
    medium::{Medium, Voxel},
    optics::{OpticalMaterial, Spectrum},
    render::{AdaptiveOptions, AreaLight, Camera, RenderOptions, Scene, Surface},
};
use std::{fs, path::PathBuf, time::Instant};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let output = PathBuf::from(args.get(1).map(String::as_str).unwrap_or("renders"));
    let parse = |i: usize, default: usize| {
        args.get(i)
            .map(|s| s.parse::<usize>().expect("integer argument"))
            .unwrap_or(default)
    };
    let n = parse(2, 40);
    let samples = parse(3, 256);
    let width = parse(4, 640);
    let steps = parse(5, 12);
    assert!(n >= 12 && width >= 32 && samples >= 16 && samples % 16 == 0);
    fs::create_dir_all(&output)?;
    let start = Instant::now();
    let water_eos = fit_inputs(&by_name("water").unwrap().inputs).liquid;
    let air = Nasg {
        gamma: 1.4,
        p_inf: 0.0,
        b: 0.0,
        q: 0.0,
        cv: 717.5,
        q_prime: 0.0,
    };
    let (p, t) = (101325.0, 298.0);
    let rhow = 1.0 / v_pt(&water_eos, p, t);
    let rhog = 1.0 / v_pt(&air, p, t);
    let sys = FiveEq {
        eos: [water_eos, air],
        kapila: true,
    };
    let len = 0.064;
    let mesh = Mesh::new([n; 3], [-len / 2.0; 3], [len / 2.0; 3]);
    let dx = mesh.dx[0];
    let center = Vec3([-0.010, 0.001, 0.0]);
    let bubble = Vec3([-0.010, 0.004, 0.003]);
    let mut solver = Solver::new(sys, mesh, |x| {
        let q = Vec3(x) - center;
        let r = (q.0[0].powi(2) / 1.08_f64.powi(2)
            + q.0[1].powi(2) / 0.93_f64.powi(2)
            + q.0[2].powi(2))
        .sqrt();
        let drop = 0.5 * (1.0 + ((0.013 - r) / (0.6 * dx)).tanh());
        let b = 0.5 * (1.0 + (((Vec3(x) - bubble).length() - 0.004) / (0.5 * dx)).tanh());
        let alpha = (drop * b).clamp(1e-7, 1.0 - 1e-7);
        // Smooth divergence-free strain supplies a genuine evolved flow.
        [
            alpha * rhow,
            (1.0 - alpha) * rhog,
            4.0 * x[0],
            -2.0 * x[1],
            -2.0 * x[2],
            p,
            alpha,
        ]
    });
    solver.bc = [[Bc::Outflow; 2]; 3];
    solver.cfl = 0.25;
    solver.surface_tension = Some(SurfaceTension {
        sigma: 0.072,
        alpha_index: 6,
        smoothing: 2,
    });
    for _ in 0..steps {
        solver.step(solver.dt_cfl());
    }
    let mut alpha = Vec::with_capacity(n * n * n);
    let mut max_speed: f64 = 0.0;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let w = solver.prim(i, j, k);
                assert!(solver.sys.admissible(&w) && w.iter().all(|v| v.is_finite()));
                alpha.push(w[6]);
                max_speed = max_speed.max((w[2] * w[2] + w[3] * w[3] + w[4] * w[4]).sqrt());
            }
        }
    }
    let water = OpticalMaterial::water();
    let iron = OpticalMaterial::iron();
    let absorption = Spectrum::new(
        water
            .k
            .0
            .iter()
            .map(|&(l, _)| (l, water.absorption(l)))
            .collect(),
    );
    let media = vec![Medium {
        absorption,
        scattering: Spectrum::constant(0.0),
        g: 0.0,
        n: water.n.clone(),
    }];
    let volume = compressible::radiation::snapshot(&solver, media, |sys, w| Voxel {
        medium: 0,
        scale: w[6].clamp(0.0, 1.0),
        temperature: sys.temperature(w),
    });
    let mut primitives = isosurface(
        solver.mesh.n,
        Vec3(solver.mesh.lo),
        Vec3(solver.mesh.dx),
        &alpha,
        0.5,
        0,
    );
    let interface_triangles = primitives.len();
    let materials = vec![
        Surface::Dielectric {
            optics: water,
            outside_n: 1.0,
        },
        Surface::Conductor {
            optics: iron,
            temperature: 1800.0,
        },
        Surface::Checker {
            a: Spectrum::constant(0.64),
            b: Spectrum::constant(0.12),
            period: 0.008,
        },
        Surface::Emitter {
            temperature: 6500.0,
            scale: 8e-5,
        },
        Surface::Emitter {
            temperature: 4500.0,
            scale: 1e-4,
        },
    ];
    primitives.push(Primitive {
        shape: Shape::Sphere {
            center: Vec3([0.017, -0.006, 0.0]),
            radius: 0.009,
        },
        material: 1,
    });
    primitives.extend(quad(
        Vec3([-0.20, -0.016, -0.20]),
        Vec3([0.40, 0.0, 0.0]),
        Vec3([0.0, 0.0, 0.40]),
        2,
    ));
    let lights = vec![
        AreaLight {
            origin: Vec3([-0.05, 0.065, -0.025]),
            u: Vec3([0.07, 0.0, 0.0]),
            v: Vec3([0.0, 0.0, 0.04]),
            temperature: 6500.0,
            scale: 8e-5,
            material: 3,
        },
        AreaLight {
            origin: Vec3([0.05, 0.0, -0.02]),
            u: Vec3([0.0, 0.05, 0.0]),
            v: Vec3([0.0, 0.0, 0.04]),
            temperature: 4500.0,
            scale: 1e-4,
            material: 4,
        },
    ];
    for (i, l) in lights.iter().enumerate() {
        primitives.extend(quad(l.origin, l.u, l.v, 3 + i));
    }
    let scene = Scene {
        geometry: Bvh::new(primitives),
        materials,
        volume: Some(volume),
        lights,
        environment_temperature: 6500.0,
        environment_scale: 5e-6,
        epsilon: len * 1e-7,
    };
    let camera = Camera::new(
        Vec3([0.051, 0.031, 0.085]),
        Vec3([-0.003, -0.003, 0.0]),
        Vec3([0.0, 1.0, 0.0]),
        38.0,
    );
    let sim_seconds = start.elapsed().as_secs_f64();
    println!("Evolved {n}³ water/air cells, {steps} fluid steps, t={:.6e} s, max |u|={max_speed:.4} m/s; {interface_triangles} interface triangles; {:.2}s preparation",solver.t,sim_seconds);
    let before = solver.u.clone();
    let render_start = Instant::now();
    let opts = RenderOptions {
        width,
        height: width * 5 / 8,
        samples,
        max_depth: 128,
        seed: 42,
    };
    let mode = args.get(6).map(String::as_str).unwrap_or("fixed");
    let image = match mode {
        "fixed" => scene.render(&camera, opts),
        "exact" => scene.render_with_budget(&camera, opts, 0.0),
        "adaptive" => scene.render_adaptive(
            &camera,
            opts,
            AdaptiveOptions {
                min_samples: 64.min(samples),
                ..Default::default()
            },
        ),
        _ => panic!("render mode must be fixed, exact or adaptive"),
    };
    assert_eq!(
        solver.u, before,
        "camera must leave simulation state unchanged"
    );
    assert!(image.rgb.iter().flatten().all(|v| v.is_finite()));
    image.write_png(output.join("photons.png"), 0.018)?;
    image.write_pfm(output.join("photons.pfm"))?;
    let mean_error =
        image.standard_error_y.iter().sum::<f64>() / image.standard_error_y.len() as f64;
    let elapsed = render_start.elapsed().as_secs_f64();
    let summary=format!("fluid_grid={n}\nfluid_steps={steps}\nfluid_time_s={}\ninterface_triangles={interface_triangles}\nwidth={}\nheight={}\nsamples_per_pixel={samples}\nseed=42\nmax_depth=128\ntruncated_paths={}\nmean_standard_error_y={mean_error}\npreparation_wall_s={sim_seconds}\nrender_wall_s={elapsed}\noptical_water_temperature_K=298\nillustrative_iron_temperature_K=1800\niron_optics=room_temperature_Johnson_Christy_not_melt_calibration\n",solver.t,opts.width,opts.height,image.truncated_paths);
    let total_paths = image.samples_per_pixel.iter().sum::<usize>();
    let summary=format!("{summary}render_mode={mode}\nmain_wavelength_paths={total_paths}\npilot_wavelength_paths={}\nmean_samples_per_pixel={}\nthermal_spectral_radiance_budget={}\n",image.pilot_paths,total_paths as f64/image.rgb.len() as f64,if mode=="exact"{0.0}else{1e-6});
    fs::write(output.join("photons.txt"), summary)?;
    println!(
        "{}×{}, {samples} spp, {:.2}s rendering, {} truncated paths, mean luminance SE {:.3e}; {}",
        opts.width,
        opts.height,
        elapsed,
        image.truncated_paths,
        mean_error,
        output.display()
    );
    Ok(())
}
