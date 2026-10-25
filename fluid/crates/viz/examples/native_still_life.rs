//! Native still life. Imagegen supplies wood and linen albedo maps; every
//! final pixel is computed by the Rust geometry, material and light models.
//! CLI: OUT WIDTH SPP [ROW_START ROW_END | assemble]
use geometry::TriMesh;
use radiation::math::Vec3;
use std::{
    f64::consts::PI,
    hash::{Hash, Hasher},
    path::PathBuf,
    sync::Arc,
};
use viz::{
    pbr::{self, AreaLight, Bsdf, Medium, Mesh, Object, Shape, Texture, World},
    Camera, Image,
};
use wetting::meniscus::Meniscus;

fn noise(x: f64, y: f64) -> f64 {
    let hash = |i: i64, j: i64| {
        let mut h = (i as u64).wrapping_mul(0x9e3779b97f4a7c15)
            ^ (j as u64).wrapping_mul(0xbf58476d1ce4e5b9);
        h = (h ^ (h >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        h = (h ^ (h >> 27)).wrapping_mul(0x94d049bb133111eb);
        ((h ^ (h >> 31)) >> 11) as f64 / (1_u64 << 53) as f64
    };
    let (i, j) = (x.floor() as i64, y.floor() as i64);
    let smooth = |t: f64| t * t * (3. - 2. * t);
    let (u, v) = (smooth(x - x.floor()), smooth(y - y.floor()));
    (hash(i, j) * (1. - u) + hash(i + 1, j) * u) * (1. - v)
        + (hash(i, j + 1) * (1. - u) + hash(i + 1, j + 1) * u) * v
}
fn perturb(n: Vec3, dx: f64, dy: f64) -> Vec3 {
    (n - Vec3([dx, dy, 0.])).unit()
}
fn mesh(mesh: TriMesh, bsdf: Bsdf) -> Object {
    Object::new(Shape::Mesh(Mesh::new(mesh, true)), bsdf)
}
fn box_object(lo: [f64; 3], hi: [f64; 3], bsdf: Bsdf) -> Object {
    Object::new(
        Shape::Mesh(Mesh::new(geometry::primitives::cuboid(lo, hi), false)),
        bsdf,
    )
}

/// Smooth meridian revolved into shared triangle geometry, not a normal-only
/// shape illusion. The full profile includes the interior and the foot.
fn lathe(points: &[(f64, f64)], angles: usize, subdiv: usize) -> TriMesh {
    let mut profile = Vec::new();
    for i in 0..points.len() - 1 {
        let p0 = points[i.saturating_sub(1)];
        let p1 = points[i];
        let p2 = points[i + 1];
        let p3 = points[(i + 2).min(points.len() - 1)];
        for j in 0..subdiv {
            let t = j as f64 / subdiv as f64;
            let cat = |a: f64, b: f64, c: f64, d: f64| {
                0.5 * (2. * b
                    + (-a + c) * t
                    + (2. * a - 5. * b + 4. * c - d) * t * t
                    + (-a + 3. * b - 3. * c + d) * t * t * t)
            };
            profile.push((
                cat(p0.0, p1.0, p2.0, p3.0).max(0.),
                cat(p0.1, p1.1, p2.1, p3.1),
            ));
        }
    }
    profile.push(*points.last().unwrap());
    let mut vertices = Vec::new();
    let mut triangles = Vec::new();
    for &(r, z) in &profile {
        for i in 0..angles {
            let a = 2. * PI * i as f64 / angles as f64;
            vertices.push([r * a.cos(), r * a.sin(), z]);
        }
    }
    for j in 0..profile.len() - 1 {
        for i in 0..angles {
            let a = (j * angles + i) as u32;
            let b = (j * angles + (i + 1) % angles) as u32;
            let c = a + angles as u32;
            let d = b + angles as u32;
            triangles.push([a, b, c]);
            triangles.push([b, d, c]);
        }
    }
    TriMesh::new(vertices, triangles)
}
fn handle() -> TriMesh {
    let (nu, nv) = (160, 32);
    let mut vertices = Vec::new();
    let mut triangles = Vec::new();
    for i in 0..nu {
        let u = 2. * PI * i as f64 / nu as f64;
        let c = Vec3([0.059 + 0.025 * u.cos(), 0., 0.066 + 0.034 * u.sin()]);
        let tangent = Vec3([-0.025 * u.sin(), 0., 0.034 * u.cos()]).unit();
        let radial = Vec3([0., 1., 0.]).cross(tangent).unit();
        for j in 0..nv {
            let v = 2. * PI * j as f64 / nv as f64;
            let r = 0.0065 * (1. + 0.045 * (3. * u).sin());
            let p = c + radial * (r * v.cos()) + Vec3([0., r * v.sin(), 0.]);
            vertices.push(p.0);
        }
    }
    for i in 0..nu {
        for j in 0..nv {
            let a = (i * nv + j) as u32;
            let b = (((i + 1) % nu) * nv + j) as u32;
            let c = (i * nv + (j + 1) % nv) as u32;
            let d = (((i + 1) % nu) * nv + (j + 1) % nv) as u32;
            triangles.push([a, b, c]);
            triangles.push([b, d, c]);
        }
    }
    let mut mesh = TriMesh::new(vertices, triangles);
    if mesh.volume() < 0. {
        mesh.flip();
    }
    mesh
}
fn cloth_height(u: f64, v: f64, layer: usize) -> f64 {
    // Estimate: prescribed folded linen, thickness 0.55 mm, not a mechanical
    // prediction. Two separate closed sheets, with lifted folds and hems.
    let ridge = 0.004 * (-((u + 0.053 + 0.003 * (v * 24.).sin()) / 0.008).powi(2)).exp();
    let curl = 0.0018 * (-((u - 0.064) / 0.005).powi(2)).exp() * (0.7 + 0.3 * (v * 35.).sin());
    let wrinkle =
        0.0012 * (u * 52. + 0.7 * (v * 18.).sin()).sin() + 0.0007 * (v * 61. + u * 13.).sin();
    0.0027 + layer as f64 * 0.0011 + ridge + curl + wrinkle
}
fn cloth_mesh(layer: usize) -> TriMesh {
    let (nu, nv) = (120, 184);
    let mut vertices = Vec::new();
    let mut triangles = Vec::new();
    let transform = |u: f64, v: f64, z: f64| {
        let (s, c) = (-0.16_f64).sin_cos();
        [-0.145 + c * u - s * v, 0.012 + s * u + c * v, z]
    };
    for side in 0..2 {
        for j in 0..=nv {
            for i in 0..=nu {
                let u = -0.069 + 0.138 * i as f64 / nu as f64;
                let v = -0.109 + 0.218 * j as f64 / nv as f64;
                let shift = if layer == 1 {
                    0.002 + 0.001 * (v * 30.).sin()
                } else {
                    0.
                };
                let hem = 0.00018 * (-((u.abs() - 0.066) / 0.0015).powi(2)).exp();
                let z = cloth_height(u, v, layer) + hem - side as f64 * 0.00055;
                vertices.push(transform(u + shift, v + layer as f64 * 0.001, z));
            }
        }
    }
    let count = (nu + 1) * (nv + 1);
    for side in 0..2 {
        for j in 0..nv {
            for i in 0..nu {
                let a = (side * count + j * (nu + 1) + i) as u32;
                let b = a + 1;
                let c = a + (nu + 1) as u32;
                let d = c + 1;
                if side == 0 {
                    triangles.push([a, b, c]);
                    triangles.push([b, d, c]);
                } else {
                    triangles.push([a, c, b]);
                    triangles.push([b, c, d]);
                }
            }
        }
    }
    let mut edge = Vec::new();
    for i in 0..=nu {
        edge.push(i);
    }
    for j in 1..=nv {
        edge.push(j * (nu + 1) + nu);
    }
    for i in (0..nu).rev() {
        edge.push(nv * (nu + 1) + i);
    }
    for j in (1..nv).rev() {
        edge.push(j * (nu + 1));
    }
    for i in 0..edge.len() {
        let a = edge[i] as u32;
        let b = edge[(i + 1) % edge.len()] as u32;
        let c = a + count as u32;
        let d = b + count as u32;
        triangles.push([a, c, b]);
        triangles.push([b, c, d]);
    }
    TriMesh::new(vertices, triangles)
}

fn build(texture: Arc<Texture>, fabric_texture: Arc<Texture>) -> World {
    // All appearance parameters below are estimates, not measured spectra.
    let ceramic = Bsdf::Coated {
        albedo: [0.63, 0.59, 0.49],
        index: 1.5,
        alpha: [0.028, 0.028],
    };
    let ceramic_shader = move |p: Vec3, n: Vec3, _footprint: f64| {
        let speck = noise(p.0[0] * 3800., p.0[1] * 3100. + p.0[2] * 4800.).powi(14);
        let variation = 0.95 + 0.05 * noise(p.0[0] * 130., p.0[2] * 220.);
        let albedo = [0.63, 0.59, 0.49].map(|v| v * variation * (1. - 0.35 * speck));
        (
            Bsdf::Coated {
                albedo,
                index: 1.5,
                alpha: [0.025 + 0.020 * speck; 2],
            },
            n,
        )
    };
    let mut objects = Vec::new();
    objects.push(
        box_object(
            [-1., -1., -0.055],
            [1., 1.2, 0.],
            Bsdf::Diffuse { albedo: [0.3; 3] },
        )
        .shader(move |p, n, footprint| {
            // Albedo only is generated. Rust computes mapping, linear filtering,
            // micro-relief, roughness, all illumination, reflections and shadows.
            let (u, v) = ((p.0[0] + 0.3) / 0.55, (p.0[1] + 0.4) / 0.9);
            let rgb = texture.sample(u, v, footprint / 0.55);
            let h = |u: f64, v: f64| pbr::luminance(texture.sample(u, v, 0.0003));
            let eps = 0.00015;
            let dx = 0.000016 * (h(u + eps, v) - h(u - eps, v)) / (2. * eps * 0.55);
            let dy = 0.000016 * (h(u, v + eps) - h(u, v - eps)) / (2. * eps * 0.9);
            (
                Bsdf::Coated {
                    albedo: rgb.map(|v| v * 0.88),
                    index: 1.47,
                    alpha: [0.29, 0.19],
                },
                if n.0[2] > 0.9 { perturb(n, dx, dy) } else { n },
            )
        }),
    );
    objects.push(box_object(
        [-1., 0.90, 0.],
        [1., 0.94, 1.3],
        Bsdf::Diffuse {
            albedo: [0.56, 0.54, 0.49],
        },
    ));
    let cup = lathe(
        &[
            (0., 0.009),
            (0.029, 0.009),
            (0.032, 0.0052),
            (0.038, 0.0052),
            (0.041, 0.010),
            (0.042, 0.026),
            (0.046, 0.075),
            (0.0495, 0.114),
            (0.050, 0.118),
            (0.049, 0.121),
            (0.046, 0.121),
            (0.045, 0.117),
            (0.0425, 0.087),
            (0.040, 0.042),
            (0.036, 0.027),
            (0., 0.027),
        ],
        192,
        8,
    );
    let cup_geometry = cup.clone();
    let cup_bvh = geometry::bvh::Bvh::build(&cup_geometry);
    objects.push(mesh(cup, ceramic).shader(ceramic_shader));
    objects.push(mesh(handle(), ceramic).shader(ceramic_shader));
    let saucer = lathe(
        &[
            (0., 0.0013),
            (0.045, 0.0013),
            (0.074, 0.003),
            (0.087, 0.008),
            (0.091, 0.012),
            (0.089, 0.015),
            (0.085, 0.015),
            (0.076, 0.010),
            (0.050, 0.0052),
            (0.036, 0.0052),
            (0., 0.0052),
        ],
        192,
        8,
    );
    objects.push(mesh(saucer, ceramic).shader(ceramic_shader));
    let fabric = Bsdf::Coated {
        albedo: [0.63, 0.59, 0.50],
        index: 1.53,
        alpha: [0.6, 0.85],
    };
    for layer in 0..2 {
        let woven = fabric_texture.clone();
        objects.push(
            mesh(cloth_mesh(layer), fabric).shader(move |p, n, footprint| {
                let dx = p.0[0] + 0.145;
                let dy = p.0[1] - 0.012;
                let (s, c) = 0.16_f64.sin_cos();
                let u = c * dx - s * dy;
                let v = s * dx + c * dy;
                // Estimate: 10 cm texture extent, luminance-derived yarn
                // relief up to 120 um. Generated albedo is not measured height.
                let (uu, vv) = ((u + 0.069) / 0.10, (v + 0.109) / 0.10);
                let albedo = woven.sample(uu, vv, footprint / 0.10).map(|v| v * 0.88);
                let height = |u: f64, v: f64| pbr::luminance(woven.sample(u, v, footprint / 0.10));
                let eps = 0.0004;
                let gu =
                    0.00012 * (height(uu + eps, vv) - height(uu - eps, vv)) / (2. * eps * 0.10);
                let gv =
                    0.00012 * (height(uu, vv + eps) - height(uu, vv - eps)) / (2. * eps * 0.10);
                let nn = if n.0[2] > 0.3 {
                    perturb(n, c * gu + s * gv, -s * gu + c * gv)
                } else {
                    n
                };
                (
                    Bsdf::Coated {
                        albedo,
                        index: 1.53,
                        alpha: [0.6, 0.85],
                    },
                    nn,
                )
            }),
        );
    }
    let steel = Bsdf::Conductor {
        index: [2.5; 3],
        extinction: [3.3; 3],
        alpha: [0.035, 0.075],
    }; // Estimate: visible steel-like optical constants.
    let mut spoon = lathe(
        &[
            (0., 0.0035),
            (0.3, 0.004),
            (0.6, 0.006),
            (0.9, 0.010),
            (1., 0.012),
            (1.025, 0.0116),
            (0.94, 0.009),
            (0.6, 0.0053),
            (0.3, 0.0033),
            (0., 0.0029),
        ],
        128,
        6,
    );
    spoon.transform(|p| [0.121 + p[0] * 0.013, -0.075 + p[1] * 0.023, p[2]]);
    objects.push(mesh(spoon, steel));
    let mut stem = geometry::primitives::cuboid([-0.002, -0.125, -0.0008], [0.002, 0., 0.0008]);
    stem.transform(|p| {
        [
            0.122 + 0.08 * p[1] + p[0],
            -0.094 + p[1],
            0.006 - 0.01 * p[1] + p[2],
        ]
    });
    objects.push(mesh(stem, steel));
    for (x, y, a) in [
        (-0.096, -0.105, 0.2_f64),
        (-0.077, -0.121, -0.5),
        (-0.097, -0.130, 0.8),
    ] {
        let mut bean = geometry::primitives::icosphere([0., 0., 0.], 1., 4);
        bean.transform(|p| {
            let (s, c) = a.sin_cos();
            let groove = 0.00065 * (-(p[0] / 0.11_f64).powi(2)).exp() * p[2].max(0.);
            let x0 = p[0] * 0.0058;
            let y0 = p[1] * 0.0037;
            [
                x + c * x0 - s * y0,
                y + s * x0 + c * y0,
                0.003 + p[2] * 0.0028 - groove,
            ]
        });
        objects.push(
            mesh(
                bean,
                Bsdf::Coated {
                    albedo: [0.055, 0.020, 0.008],
                    index: 1.45,
                    alpha: [0.28; 2],
                },
            )
            .shader(move |p, n, _| {
                let t = 0.75 + 0.25 * noise(p.0[0] * 2900., p.0[1] * 4700.);
                (
                    Bsdf::Coated {
                        albedo: [0.065 * t, 0.022 * t, 0.009 * t],
                        index: 1.45,
                        alpha: [0.28; 2],
                    },
                    n,
                )
            }),
        );
    }
    // Water constants are sourced in wetting::constants (IAPWS/CRC,
    // Vargaftik et al. 1983, CGPM). Coffee surface tension/density are folded
    // here into a water-like estimate, not measured coffee composition.
    let base = 0.109;
    let wall_height = 0.00065;
    let radius = cup_bvh
        .ray(
            &cup_geometry,
            [0., 0., base + wall_height],
            [1., 0., 0.],
            0.,
            1.,
        )
        .unwrap()
        .distance;
    let meniscus = Arc::new(
        Meniscus::solve(
            radius,
            wetting::constants::water_capillary_length(),
            wall_height,
            2048,
        )
        .unwrap(),
    );
    // Estimate: dark coffee absorption in m^-1 and visible refractive index.
    let coffee = Medium {
        index: 1.333,
        absorption: [80., 190., 350.],
    };
    let mut bubbles = Vec::new();
    for i in 0..44 {
        let a = 0.15 + i as f64 * 0.078;
        let r = 0.00025 + 0.00048 * noise(i as f64 * 2.13, 7.8);
        let radial = 0.0434 - 1.1 * r;
        let (x, y) = (radial * a.cos(), radial * a.sin());
        let centre = Vec3([x, y, base + meniscus.height(radial) - 0.16 * r]);
        bubbles.push((centre, r));
    }
    // Estimated sparse cluster on the front side, with non-overlapping shells.
    for (x, y, r) in [
        (-0.023_f64, -0.018_f64, 0.00098),
        (-0.0252, -0.0188, 0.00075),
        (-0.027, -0.0198, 0.00060),
        (-0.0245, -0.0165, 0.00055),
    ] {
        let centre = Vec3([x, y, base + meniscus.height(x.hypot(y)) - 0.16 * r]);
        bubbles.push((centre, r));
    }
    let bubbles = Arc::new(bubbles);
    let profile: Vec<_> = (0..=160)
        .map(|i| {
            let r = meniscus.radius * i as f64 / 160.;
            (r, base + meniscus.height(r))
        })
        .collect();
    // Reverse radial traversal for an upward outward normal, then close the
    // liquid at the cup bottom. These meshes remain ordinary project TriMesh.
    // Conform the liquid volume to the actual hollow cup, rather than letting
    // a straight cylinder pierce the tapered ceramic wall. The sub-micron
    // inset resolves coincident mesh interfaces. Below the free surface, this
    // boundary uses a wet opaque ceramic BRDF and retains the coffee medium.
    let mut volume = vec![(0., 0.028)];
    for i in 0..=96 {
        let z = 0.028 + (base + wall_height - 0.028) * i as f64 / 96.;
        let r = cup_bvh
            .ray(&cup_geometry, [0., 0., z], [1., 0., 0.], 0., 1.)
            .unwrap()
            .distance;
        volume.push((r - 8e-7, z));
    }
    volume.extend(profile.iter().rev().copied());
    let liquid = lathe(&volume, 256, 1);
    let cuts = bubbles.clone();
    let surface = meniscus.clone();
    objects.push(
        mesh(
            liquid,
            Bsdf::Dielectric {
                outside: Medium::AIR,
                inside: coffee,
            },
        )
        .accept(move |p| cuts.iter().all(|(c, r)| (p - *c).length() >= *r))
        .shader(move |p, n, _| {
            let r = p.0[0].hypot(p.0[1]);
            if p.0[2] < base + surface.height(r) - 2e-6 {
                (
                    Bsdf::Coated {
                        albedo: [0.63, 0.59, 0.49],
                        index: 1.5 / coffee.index,
                        alpha: [0.055; 2],
                    },
                    n,
                )
            } else {
                (
                    Bsdf::Dielectric {
                        outside: Medium::AIR,
                        inside: coffee,
                    },
                    n,
                )
            }
        }),
    );
    for &(centre, radius) in bubbles.iter() {
        let profile = meniscus.clone();
        // Explicit shell, not painted circles. Radius/film thickness are
        // estimates; 8 um film has real front/back interfaces. Bo < 0.14,
        // spherical curvature is a capillary-dominated geometric closure.
        objects.push(
            Object::new(
                Shape::Sphere { centre, radius },
                Bsdf::Dielectric {
                    outside: Medium::AIR,
                    inside: coffee,
                },
            )
            .accept(move |p| p.0[2] > base + profile.height(p.0[0].hypot(p.0[1]))),
        );
        objects.push(Object::new(
            Shape::Sphere {
                centre,
                radius: radius - 8e-6,
            },
            Bsdf::Dielectric {
                outside: coffee,
                inside: Medium::AIR,
            },
        ));
    }
    // An actual eight-pane window: finite luminous sky patches, opaque wooden
    // mullions and reveal. Radiance is an estimate in relative RGB units.
    let mut lights = Vec::new();
    for j in 0..2 {
        for i in 0..4 {
            lights.push(AreaLight {
                corner: Vec3([-0.56, -0.55 + i as f64 * 0.36, 0.20 + j as f64 * 0.43]),
                u: Vec3([0., 0.33, 0.]),
                v: Vec3([0., 0., 0.40]),
                radiance: if j == 1 {
                    [8.0, 8.5, 9.2]
                } else {
                    [5.7, 5.8, 5.9]
                },
            });
        }
    }
    let frame = Bsdf::Diffuse {
        albedo: [0.38, 0.35, 0.29],
    };
    for i in 1..4 {
        let y = -0.55 + i as f64 * 0.36;
        objects.push(box_object(
            [-0.575, y - 0.03, 0.17],
            [-0.54, y, 1.07],
            frame,
        ));
    }
    objects.push(box_object(
        [-0.575, -0.58, 0.60],
        [-0.54, 0.89, 0.63],
        frame,
    ));
    objects.push(box_object(
        [-0.575, -0.58, 0.17],
        [-0.54, 0.89, 0.20],
        frame,
    ));
    objects.push(box_object(
        [-0.575, -0.58, 1.03],
        [-0.54, 0.89, 1.07],
        frame,
    ));
    // Estimate: very weak room fill; all visible key light comes from the panes.
    World::new(objects, lights, [0.025, 0.028, 0.032])
}

fn display(image: &Image) -> Image {
    // Exposure and highlight shoulder only; no generative editing.
    Image {
        width: image.width,
        height: image.height,
        pixels: image
            .pixels
            .iter()
            .map(|p| p.map(|v| 1. - (-2.1 * v.max(0.)).exp()))
            .collect(),
    }
}
fn finish(
    out: &std::path::Path,
    image: &Image,
    texture: Arc<Texture>,
    fabric: Arc<Texture>,
    detail: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    image.write_pfm(out.join("native-linear.pfm"))?;
    display(image).write_png(out.join("native-unfiltered.png"), 1.)?;
    let world = build(texture, fabric);
    let camera = camera(detail);
    let developed = pbr::denoise(&world, &camera, image);
    display(&developed).write_png(out.join("native.png"), 1.)?;
    Ok(())
}
fn camera(detail: bool) -> Camera {
    if detail {
        let camera = Camera::looking_at([0.15, -0.24, 0.26], [0., 0., 0.110], 20_f64.to_radians());
        return camera.with_depth_of_field(0.0015, camera.subject_distance());
    }
    let camera = Camera::looking_at(
        [0.28, -0.43, 0.31],
        [-0.020, -0.016, 0.05],
        35_f64.to_radians(),
    );
    // Estimate: 3.5 mm pupil; an ideal thin-lens camera, not a measured lens.
    camera.with_depth_of_field(0.0035, camera.subject_distance())
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args: Vec<_> = std::env::args().skip(1).collect();
    let detail = args.last().is_some_and(|s| s == "--detail");
    if detail {
        args.pop();
    }
    let out = PathBuf::from(
        args.first()
            .map(String::as_str)
            .unwrap_or("../analysis/photographic_scene/native-v2"),
    );
    let width = args
        .get(1)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(1000usize);
    let height = width * 2 / 3;
    let samples = args
        .get(2)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(512usize);
    if width < 3 || samples == 0 {
        return Err("width >= 3 and positive samples required".into());
    }
    let wood = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../analysis/photographic_scene/assets/oak-albedo.png");
    let linen = wood.with_file_name("linen-albedo.png");
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    include_str!("native_still_life.rs").hash(&mut hasher);
    include_str!("../src/pbr.rs").hash(&mut hasher);
    include_str!("../../wetting/src/meniscus.rs").hash(&mut hasher);
    include_str!("../../geometry/src/bvh.rs").hash(&mut hasher);
    std::fs::read(&wood)?.hash(&mut hasher);
    std::fs::read(&linen)?.hash(&mut hasher);
    let key = format!("{:016x}", hasher.finish());
    let manifest = serde_json::json!({"width":width,"height":height,"samples":samples,"render_key":key,
        "view":if detail {"liquid detail"}else{"still life"},"final_pixels":"Rust only","generated_inputs":["oak-albedo.png","linen-albedo.png"],"postprocess":"albedo/geometry-guided spatial denoising, exposure, highlight shoulder, sRGB; original radiance retained"});
    std::fs::create_dir_all(&out)?;
    let manifest_path = out.join("render.json");
    if manifest_path.exists() {
        let old: serde_json::Value = serde_json::from_slice(&std::fs::read(&manifest_path)?)?;
        if old != manifest {
            return Err("output has different settings/source; choose a new directory".into());
        }
    } else {
        std::fs::write(&manifest_path, serde_json::to_string_pretty(&manifest)?)?;
    }
    if args.get(3).is_some_and(|s| s == "assemble") {
        let mut tiles = Vec::new();
        for entry in std::fs::read_dir(&out)? {
            let p = entry?.path();
            let stem = p.file_stem().and_then(|s| s.to_str()).unwrap_or("");
            if let Some(start) = stem
                .strip_prefix("rows-")
                .and_then(|s| s.split('-').next())
                .and_then(|s| s.parse::<usize>().ok())
            {
                if p.extension().is_some_and(|s| s == "pfm") {
                    tiles.push((start, p));
                }
            }
        }
        tiles.sort_by_key(|p| p.0);
        let mut image = Image::new(width, height);
        let mut cursor = 0;
        for (start, path) in tiles {
            let tile = pbr::read_pfm(path)?;
            if start != cursor || tile.width != width || start + tile.height > height {
                return Err("tiles overlap or leave a gap".into());
            }
            image.pixels[start * width..(start + tile.height) * width]
                .copy_from_slice(&tile.pixels);
            cursor += tile.height;
        }
        if cursor != height {
            return Err("missing render rows".into());
        }
        finish(
            &out,
            &image,
            Arc::new(Texture::from_png(&wood)?),
            Arc::new(Texture::from_png(&linen)?),
            detail,
        )?;
        return Ok(());
    }
    let start = args
        .get(3)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(0usize);
    let end = args
        .get(4)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(height);
    if start >= end || end > height {
        return Err("invalid row range".into());
    }
    let texture = Arc::new(Texture::from_png(&wood)?);
    let fabric = Arc::new(Texture::from_png(&linen)?);
    let world = build(texture.clone(), fabric.clone());
    let camera = camera(detail);
    let clock = std::time::Instant::now();
    let image = pbr::render_rows(&world, &camera, width, height, samples, start, end);
    image.write_pfm(out.join(format!("rows-{start}-{end}.pfm")))?;
    if start == 0 && end == height {
        finish(&out, &image, texture, fabric, detail)?;
    }
    eprintln!("rows {start}..{end}, {:.3}s", clock.elapsed().as_secs_f64());
    Ok(())
}
