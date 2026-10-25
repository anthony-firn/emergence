//! Independent exact/published checks for the native photographic renderer.
use geometry::{bvh::Bvh, primitives};
use radiation::math::{Ray, Rng, Vec3};
use std::f64::consts::PI;
use viz::{
    pbr::{ggx_d, AreaLight, Bsdf, Medium, Mesh, Object, Shape, Texture, World},
    Image,
};

#[test]
fn mesh_rays_match_exact_box_and_do_not_crack_at_shared_edges() {
    let mesh = primitives::cuboid([-1.; 3], [1.; 3]);
    let bvh = Bvh::build(&mesh);
    let mut rng = Rng(97);
    let mut worst = 0_f64;
    for _ in 0..2000 {
        let x = -1. + 2. * rng.uniform();
        let y = -1. + 2. * rng.uniform();
        let hit = bvh.ray(&mesh, [x, y, 3.], [0., 0., -1.], 0., 10.).unwrap();
        worst = worst.max((hit.distance - 2.).abs());
        assert!((hit.barycentric.iter().sum::<f64>() - 1.).abs() < 1e-14);
    }
    for p in [[0., 0., 3.], [-1., -1., 3.], [1., 1., 3.], [0.5, 0.5, 3.]] {
        assert!((bvh.ray(&mesh, p, [0., 0., -1.], 0., 10.).unwrap().distance - 2.).abs() < 1e-14);
    }
    assert!(bvh
        .ray(&mesh, [1.001, 0., 3.], [0., 0., -1.], 0., 10.)
        .is_none());
    eprintln!("mesh vs exact box: max distance error {worst:e} m");
    assert!(worst < 1e-14);
}

#[test]
fn ggx_projected_density_normalizes_and_brdf_is_reciprocal() {
    let mut worst = 0_f64;
    for alpha in [[0.1, 0.3], [0.3, 0.7], [0.8, 0.8]] {
        let mut integral = 0.;
        let (nz, np) = (4096, 256);
        // Integrate in theta; concentration near the pole is resolved by the
        // angle grid, rather than assuming uniform solid angle is sufficient.
        for i in 0..nz {
            let t = (i as f64 + 0.5) * PI / (2. * nz as f64);
            for j in 0..np {
                let p = (j as f64 + 0.5) * 2. * PI / np as f64;
                integral += ggx_d(Vec3([t.sin() * p.cos(), t.sin() * p.sin(), t.cos()]), alpha)
                    * t.cos()
                    * t.sin()
                    * PI
                    / (2. * nz as f64)
                    * 2.
                    * PI
                    / np as f64;
            }
        }
        let error = (integral - 1.).abs();
        eprintln!("GGX {alpha:?} projected density {integral:.12}, error {error:e}");
        worst = worst.max(error);
        let bsdf = Bsdf::Coated {
            albedo: [0.7, 0.6, 0.5],
            index: 1.5,
            alpha,
        };
        let a = Vec3([0.4, -0.2, 1.]).unit();
        let b = Vec3([-0.1, 0.7, 1.]).unit();
        let x = bsdf.evaluate(a, b).0;
        let y = bsdf.evaluate(b, a).0;
        for k in 0..3 {
            assert!((x[k] - y[k]).abs() < 1e-14);
        }
    }
    assert!(worst < 5e-6);
}

#[test]
fn ggx_sampling_matches_independent_quadrature_and_conserves_energy() {
    let bsdf = Bsdf::Coated {
        albedo: [0.7; 3],
        index: 1.5,
        alpha: [0.23, 0.41],
    };
    let wo = Vec3([0.5, 0., 1.]).unit();
    let mut exact = 0.;
    let (nz, np) = (800, 512);
    for i in 0..nz {
        let z = (i as f64 + 0.5) / nz as f64;
        for j in 0..np {
            let p = (j as f64 + 0.5) * 2. * PI / np as f64;
            let wi = Vec3([
                (1. - z * z).sqrt() * p.cos(),
                (1. - z * z).sqrt() * p.sin(),
                z,
            ]);
            exact += bsdf.evaluate(wo, wi).0[0] * z * 2. * PI / (nz * np) as f64;
        }
    }
    let mut measured = 0.;
    let mut rng = Rng(578102);
    let n = 250000;
    for _ in 0..n {
        if let Some((_, weight, _)) = bsdf.sample(wo, &mut rng) {
            measured += weight[0] / n as f64;
        }
    }
    let error = (measured - exact).abs();
    eprintln!("coated GGX integral {exact:.9}, sampled {measured:.9}, error {error:e}");
    assert!(exact > 0. && exact <= 1.);
    assert!(error < 0.004);
}

#[test]
fn furnace_and_absorption_recover_exact_results() {
    let world = World::new(
        vec![Object::new(
            Shape::Quad {
                corner: Vec3([-1., -1., 0.]),
                u: Vec3([2., 0., 0.]),
                v: Vec3([0., 2., 0.]),
            },
            Bsdf::Diffuse {
                albedo: [0.7, 0.5, 0.3],
            },
        )],
        vec![],
        [1.; 3],
    );
    let mut rng = Rng(456);
    let ray = Ray::new(Vec3([0., 0., 0.5]), Vec3([0., 0., -1.]));
    for _ in 0..1000 {
        let value = world.radiance(ray, &mut rng, 8);
        for (v, e) in value.into_iter().zip([0.7, 0.5, 0.3]) {
            assert!((v - e).abs() < 1e-14);
        }
    }
    let medium = Medium {
        index: 1.333,
        absorption: [2., 5., 10.],
    };
    for (v, exact) in
        medium
            .transmission(0.2)
            .into_iter()
            .zip([(-0.4_f64).exp(), (-1_f64).exp(), (-2_f64).exp()])
    {
        assert!((v - exact).abs() < 1e-15);
    }
    // Clear plane-parallel slab in a uniform furnace: entry/exit eta²
    // factors cancel, and reflection plus transmission returns unity.
    let slab = World::new(
        vec![Object::new(
            Shape::Mesh(Mesh::new(
                primitives::cuboid([-1., -1., -0.02], [1., 1., 0.]),
                false,
            )),
            Bsdf::Dielectric {
                outside: Medium::AIR,
                inside: Medium {
                    index: 1.333,
                    absorption: [0.; 3],
                },
            },
        )],
        vec![],
        [1.; 3],
    );
    for _ in 0..1000 {
        for v in slab.radiance(ray, &mut rng, 24) {
            assert!((v - 1.).abs() < 1e-14);
        }
    }
}

#[test]
fn rectangle_transport_matches_exact_irradiance_with_mis() {
    let light = AreaLight {
        corner: Vec3([-0.475, 0.5, 1.]),
        u: Vec3([0.95, 0., 0.]),
        v: Vec3([0., -1., 0.]),
        radiance: [1.; 3],
    };
    let world = World::new(
        vec![Object::new(
            Shape::Quad {
                corner: Vec3([-2., -2., 0.]),
                u: Vec3([4., 0., 0.]),
                v: Vec3([0., 4., 0.]),
            },
            Bsdf::Diffuse { albedo: [0.5; 3] },
        )],
        vec![light],
        [0.; 3],
    );
    let (a, b, d) = (0.475_f64, 0.5_f64, 1.0_f64);
    let exact = 0.5 * 2. / PI
        * (a / (a * a + d * d).sqrt() * (b / (a * a + d * d).sqrt()).atan()
            + b / (b * b + d * d).sqrt() * (a / (b * b + d * d).sqrt()).atan());
    let ray = Ray::new(Vec3([0., 0., 0.5]), Vec3([0., 0., -1.]));
    let mut rng = Rng(78279);
    let n = 150000;
    let mut measured = 0.;
    for _ in 0..n {
        measured += world.radiance(ray, &mut rng, 4)[0] / n as f64;
    }
    let error = (measured - exact).abs();
    eprintln!("MIS rectangle: {measured:.9}, exact {exact:.9}, absolute error {error:e}");
    assert!(error < 0.001);
}

#[test]
fn linear_texture_filtering_and_pfm_rows_round_trip() {
    let t = Texture::from_linear(2, 2, vec![[0.; 3], [1.; 3], [1.; 3], [0.; 3]]);
    for v in t.sample(0.5, 0.5, 1.) {
        assert!((v - 0.5).abs() < 1e-14);
    }
    assert_eq!(t.sample(0.25, 0.25, 0.), [0.; 3]);
    assert_eq!(t.sample(1.25, -0.75, 0.), [0.; 3]);
    let mut image = Image::new(2, 2);
    image.pixels = vec![[0., 0.5, 1.], [2., 3., 4.], [5., 6., 7.], [8., 9., 10.]];
    let path = std::env::temp_dir().join(format!("native-pbr-pfm-{}.pfm", std::process::id()));
    image.write_pfm(&path).unwrap();
    let read = viz::pbr::read_pfm(&path).unwrap();
    std::fs::remove_file(path).unwrap();
    assert_eq!(read.pixels, image.pixels);
}

#[test]
fn hollow_bubble_interfaces_preserve_radiance_in_a_furnace() {
    let water = Medium {
        index: 1.333,
        absorption: [0.; 3],
    };
    let world = World::new(
        vec![
            Object::new(
                Shape::Sphere {
                    centre: Vec3([0.; 3]),
                    radius: 0.05,
                },
                Bsdf::Dielectric {
                    outside: Medium::AIR,
                    inside: water,
                },
            ),
            Object::new(
                Shape::Sphere {
                    centre: Vec3([0.; 3]),
                    radius: 0.05 - 8e-6,
                },
                Bsdf::Dielectric {
                    outside: water,
                    inside: Medium::AIR,
                },
            ),
        ],
        vec![],
        [1.; 3],
    );
    let ray = Ray::new(Vec3([0.01, 0., 0.3]), Vec3([0., 0., -1.]));
    let mut rng = Rng(778122);
    for _ in 0..10000 {
        for v in world.radiance(ray, &mut rng, 24) {
            assert!((v - 1.).abs() < 1e-13);
        }
    }
}

#[test]
fn restarting_rows_is_identical_to_a_single_render() {
    let world = World::new(
        vec![Object::new(
            Shape::Sphere {
                centre: Vec3([0.; 3]),
                radius: 0.24,
            },
            Bsdf::Coated {
                albedo: [0.4, 0.2, 0.1],
                index: 1.5,
                alpha: [0.2, 0.4],
            },
        )],
        vec![AreaLight {
            corner: Vec3([-0.5, 0.3, 0.8]),
            u: Vec3([1., 0., 0.]),
            v: Vec3([0., -0.5, 0.]),
            radiance: [4.; 3],
        }],
        [0.03; 3],
    );
    let camera =
        viz::Camera::looking_at([0., -1., 0.5], [0.; 3], 0.5).with_depth_of_field(0.01, 1.);
    let full = viz::pbr::render_rows(&world, &camera, 16, 12, 8, 0, 12);
    let a = viz::pbr::render_rows(&world, &camera, 16, 12, 8, 0, 5);
    let b = viz::pbr::render_rows(&world, &camera, 16, 12, 8, 5, 12);
    assert_eq!(&full.pixels[..80], &a.pixels);
    assert_eq!(&full.pixels[80..], &b.pixels);
}
