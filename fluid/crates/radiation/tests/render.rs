use radiation::{geometry::*, math::*, medium::*, optics::*, render::*};
fn empty_scene() -> Scene {
    Scene {
        geometry: Bvh::new(vec![]),
        materials: vec![],
        volume: None,
        lights: vec![],
        environment_temperature: 500.0,
        environment_scale: 1.0,
        epsilon: 1e-8,
    }
}
#[test]
fn slab_camera_matches_exact_emitting_transfer_solution() {
    let mut scene = empty_scene();
    scene.volume = Some(Volume::new(
        [8, 1, 1],
        Bounds {
            lo: Vec3([0.0, -1.0, -1.0]),
            hi: Vec3([1.0, 1.0, 1.0]),
        },
        vec![Medium::grey(2.0, 0.0, 0.0)],
        vec![
            Voxel {
                medium: 0,
                scale: 1.0,
                temperature: 1200.0
            };
            8
        ],
    ));
    let l = 700e-9;
    let expected =
        planck(l, 1200.0) * (1.0 - (-2.0_f64).exp()) + planck(l, 500.0) * (-2.0_f64).exp();
    let (actual, cap) = scene.radiance(
        Ray::new(Vec3([-1.0, 0.0, 0.0]), Vec3([1.0, 0.0, 0.0])),
        l,
        &mut Rng(42),
        64,
    );
    assert!(!cap);
    assert!((actual / expected - 1.0).abs() < 1e-14);
}
#[test]
fn isosurface_plane_position_orientation_and_bvh() {
    let n = [8; 3];
    let dx = Vec3([0.125; 3]);
    let alpha = (0..8)
        .flat_map(|i| (0..64).map(move |_| 1.0 - (i as f64 + 0.5) * 0.125))
        .collect::<Vec<_>>();
    let mesh = isosurface(n, Vec3([0.0; 3]), dx, &alpha, 0.63, 0);
    let tree = Bvh::new(mesh);
    let r = Ray::new(Vec3([1.0, 0.371, 0.419]), Vec3([-1.0, 0.0, 0.0]));
    let h = tree.hit(r, 100.0).unwrap();
    assert!((h.point.0[0] - 0.37).abs() < 1e-14);
    assert!(h.normal.0[0] > 0.999999);
}
#[test]
fn renderer_is_reproducible_across_thread_counts_and_writes_hdr() {
    let scene = empty_scene();
    let camera = Camera::new(
        Vec3([0.0, 0.0, 2.0]),
        Vec3([0.0; 3]),
        Vec3([0.0, 1.0, 0.0]),
        40.0,
    );
    let opts = RenderOptions {
        width: 4,
        height: 3,
        samples: 32,
        max_depth: 64,
        seed: 42,
    };
    let a = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(|| scene.render(&camera, opts));
    let b = rayon::ThreadPoolBuilder::new()
        .num_threads(3)
        .build()
        .unwrap()
        .install(|| scene.render(&camera, opts));
    assert_eq!(a.rgb, b.rgb);
    assert!(a.rgb.iter().flatten().all(|v| v.is_finite()));
    assert_eq!(a.truncated_paths, 0);
    let path = std::env::temp_dir().join(format!("radiation-test-{}.pfm", std::process::id()));
    a.write_pfm(&path).unwrap();
    let data = std::fs::read(&path).unwrap();
    assert!(data.starts_with(b"PF\n4 3\n-1.0\n"));
    std::fs::remove_file(path).unwrap();
}

#[test]
fn diffuse_camera_matches_lambert_under_uniform_environment() {
    let mut scene = empty_scene();
    scene.geometry = Bvh::new(quad(
        Vec3([-10.0, -10.0, 0.0]),
        Vec3([20.0, 0.0, 0.0]),
        Vec3([0.0, 20.0, 0.0]),
        0,
    ));
    scene.materials = vec![Surface::Diffuse {
        albedo: Spectrum::constant(0.5),
    }];
    let l = 600e-9;
    for seed in 0..100 {
        let (v, cap) = scene.radiance(
            Ray::new(Vec3([0.0, 0.0, 1.0]), Vec3([0.0, 0.0, -1.0])),
            l,
            &mut Rng(seed),
            64,
        );
        assert!(!cap);
        assert!((v / (0.5 * planck(l, scene.environment_temperature)) - 1.0).abs() < 1e-14);
    }
}

#[test]
fn opaque_thermal_surface_obeys_kirchhoff() {
    let mut scene = empty_scene();
    scene.environment_scale = 0.0;
    scene.geometry = Bvh::new(quad(
        Vec3([-10.0, -10.0, 0.0]),
        Vec3([20.0, 0.0, 0.0]),
        Vec3([0.0, 20.0, 0.0]),
        0,
    ));
    let iron = OpticalMaterial::iron();
    let l = 650e-9;
    let expected = (1.0 - fresnel_conductor(1.0, iron.n.at(l), iron.k.at(l))) * planck(l, 1800.0);
    scene.materials = vec![Surface::Conductor {
        optics: iron,
        temperature: 1800.0,
    }];
    let (v, cap) = scene.radiance(
        Ray::new(Vec3([0.0, 0.0, 1.0]), Vec3([0.0, 0.0, -1.0])),
        l,
        &mut Rng(2),
        64,
    );
    assert!(!cap);
    assert!((v / expected - 1.0).abs() < 1e-14);
}

#[test]
fn closed_dielectric_returns_vacuum_radiance_without_eta_drift() {
    let mut scene = empty_scene();
    scene.geometry = Bvh::new(vec![Primitive {
        shape: Shape::Sphere {
            center: Vec3([0.0; 3]),
            radius: 1.0,
        },
        material: 0,
    }]);
    scene.materials = vec![Surface::Dielectric {
        optics: OpticalMaterial::water(),
        outside_n: 1.0,
    }];
    let l = 550e-9;
    let expected = planck(l, scene.environment_temperature);
    for seed in 0..1000 {
        let (v, cap) = scene.radiance(
            Ray::new(Vec3([0.0, 0.0, 2.0]), Vec3([0.0, 0.0, -1.0])),
            l,
            &mut Rng(seed),
            64,
        );
        assert!(!cap);
        assert!((v / expected - 1.0).abs() < 1e-12);
    }
}
