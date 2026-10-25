use radiation::{
    geometry::*, math::*, medium::*, optics::*, render::*, sampling::WavelengthSampler,
};
#[test]
fn tiny_and_large_triangles_accept_shared_edges_and_vertices() {
    for scale in [1e-12, 1e-9, 1.0, 1e9] {
        let b = Bvh::new(quad(
            Vec3([0.0; 3]),
            Vec3([scale, 0.0, 0.0]),
            Vec3([0.0, scale, 0.0]),
            0,
        ));
        for x in [0.0, 0.125, 0.5, 0.875, 1.0] {
            // Exactly shared diagonal, including both shared vertices.
            let r = Ray::new(Vec3([x * scale, x * scale, scale]), Vec3([0.0, 0.0, -1.0]));
            let h = b.hit(r, 2.0 * scale).unwrap();
            assert!((h.distance / scale - 1.0).abs() < 1e-14);
            assert!(b.occluded(r, 2.0 * scale));
            assert!(!b.occluded(r, 0.5 * scale));
        }
    }
}
#[test]
fn accelerated_hits_and_occlusion_agree_with_exhaustive_search() {
    let primitives: Vec<_> = (0..80)
        .map(|i| Primitive {
            shape: Shape::Sphere {
                center: Vec3([(i % 8) as f64 * 0.3, (i / 8) as f64 * 0.3, 0.0]),
                radius: 0.12,
            },
            material: i,
        })
        .collect();
    let brute: Vec<_> = primitives
        .iter()
        .map(|p| Bvh::new(vec![p.clone()]))
        .collect();
    let b = Bvh::new(primitives);
    let mut rng = Rng(4);
    for _ in 0..500 {
        let ray = Ray::new(
            Vec3([rng.uniform() * 3.0 - 0.3, rng.uniform() * 3.0 - 0.3, 2.0]),
            Vec3([
                0.15 * (rng.uniform() - 0.5),
                0.15 * (rng.uniform() - 0.5),
                -1.0,
            ]),
        );
        let expected = brute
            .iter()
            .filter_map(|b| b.hit(ray, 10.0))
            .min_by(|a, b| a.distance.total_cmp(&b.distance));
        let actual = b.hit(ray, 10.0);
        assert_eq!(actual.map(|h| h.material), expected.map(|h| h.material));
        for limit in [1.0, 1.95, 10.0] {
            assert_eq!(
                b.occluded(ray, limit),
                expected.is_some_and(|h| h.distance < limit)
            );
        }
    }
}
fn scene() -> Scene {
    Scene {
        geometry: Bvh::new(vec![]),
        materials: vec![],
        lights: vec![],
        volume: None,
        environment_temperature: 6500.0,
        environment_scale: 1e-5,
        epsilon: 1e-9,
    }
}
#[test]
fn omitted_thermal_emission_respects_the_budget_per_weighted_path() {
    let mut scene = scene();
    scene.environment_scale = 0.0;
    let mut medium = Medium::grey(3.0, 10.0, 0.7);
    medium.n = Spectrum::constant(1.3);
    scene.volume = Some(Volume::new(
        [16, 1, 1],
        Bounds {
            lo: Vec3([0.0, -10.0, -10.0]),
            hi: Vec3([1.0, 10.0, 10.0]),
        },
        vec![medium],
        (0..16)
            .map(|i| Voxel {
                medium: 0,
                scale: 1.0,
                temperature: 500.0 + i as f64 * 50.0,
            })
            .collect(),
    ));
    let ray = Ray::new(Vec3([0.0; 3]), Vec3([1.0, 0.0, 0.0]));
    for seed in 0..100 {
        for budget in [0.0, 1e-6, 1e6, 1e12] {
            let (exact, cap) = scene.radiance(ray, 2e-6, &mut Rng(seed), 128);
            let (fast, cap_fast) =
                scene.radiance_with_budget(ray, 2e-6, &mut Rng(seed), 128, budget);
            assert_eq!(cap, cap_fast);
            let roundoff = exact.abs() * 1e-13;
            assert!(
                exact - fast >= -roundoff && exact - fast <= budget + roundoff,
                "error {} budget {budget}",
                exact - fast
            );
            if budget == 0.0 {
                assert_eq!(exact, fast);
            }
        }
    }
}
#[test]
fn spectral_importance_sampling_retains_support_and_correct_pdf() {
    let sampler = WavelengthSampler::default();
    let n = 200000;
    let mut width = 0.0;
    let mut y = 0.0;
    for i in 0..n {
        let (l, p) = sampler.sample((i as f64 + 0.5) / n as f64);
        assert!((380e-9..780e-9).contains(&l) && p > 0.0);
        width += 1.0 / p;
        y += cie_xyz(l)[1] / p;
    }
    assert!((width / n as f64 / 400e-9 - 1.0).abs() < 1e-4);
    let exact = (0..10000)
        .map(|i| cie_xyz((380.0 + 400.0 * (i as f64 + 0.5) / 10000.0) * 1e-9)[1])
        .sum::<f64>()
        * 400e-9
        / 10000.0;
    assert!((y / n as f64 / exact - 1.0).abs() < 1e-4);
    assert!(sampler.sample(1e-12).0 < 380.01e-9);
    assert!(sampler.sample(1.0 - 1e-12).0 > 779.99e-9);
}
#[test]
fn adaptive_independent_main_samples_match_blackbody_colour_integral() {
    let scene = scene();
    let camera = Camera::new(
        Vec3([0.0, 0.0, 1.0]),
        Vec3([0.0; 3]),
        Vec3([0.0, 1.0, 0.0]),
        40.0,
    );
    let options = RenderOptions {
        width: 32,
        height: 16,
        samples: 512,
        max_depth: 64,
        seed: 15,
    };
    let adaptive = AdaptiveOptions {
        pilot_samples: 32,
        min_samples: 64,
        ..Default::default()
    };
    let image = scene.render_adaptive(&camera, options, adaptive);
    let mut expected = [0.0; 3];
    let n = 20000;
    for i in 0..n {
        let l = (380.0 + 400.0 * (i as f64 + 0.5) / n as f64) * 1e-9;
        let xyz = cie_xyz(l);
        for d in 0..3 {
            expected[d] += xyz[d] * planck(l, 6500.0) * 1e-5 * 400e-9 / n as f64;
        }
    }
    let expected = xyz_to_rgb(expected);
    for (d, e) in expected.into_iter().enumerate() {
        let mean = image.rgb.iter().map(|p| p[d]).sum::<f64>() / image.rgb.len() as f64;
        assert!(
            (mean / e - 1.0).abs() < 0.005,
            "channel {d}: {mean} expected {e}"
        );
    }
    assert!(image
        .samples_per_pixel
        .iter()
        .all(|&n| (64..=512).contains(&n) && n.is_multiple_of(16)));
    assert_eq!(image.pilot_paths, 32 * 16 * 32);
    assert_eq!(image.truncated_paths, 0);
    let other = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(|| scene.render_adaptive(&camera, options, adaptive));
    assert_eq!(image.rgb, other.rgb);
    assert_eq!(image.samples_per_pixel, other.samples_per_pixel);
}
