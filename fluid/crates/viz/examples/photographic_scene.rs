//! A reproducible still life using this project's CAD, camera, ray marcher,
//! Fresnel optics, cosine sampler, and PNG writer. No generated image assets.
//! RGB textures are artistic reflectance estimates, not measured spectra.
//! Run: cargo run -p viz --release --example photographic_scene -- OUT WIDTH SPP

use cad::solid::Solid;
use radiation::{
    math::{around, Ray, Rng, Vec3},
    optics::{fresnel_conductor, fresnel_dielectric},
};
use rayon::prelude::*;
use std::{
    f64::consts::PI,
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
};
use viz::{Camera, Image, MarchOptions, Material, Object, Scene};

struct Studio {
    scene: Scene,
    march: MarchOptions,
}

fn studio() -> Studio {
    let mut scene = Scene::new();
    // All reflectances and styling dimensions below are estimates for an
    // illustrative tabletop, not calibrated samples of specific materials.
    scene.materials = vec![Material::Diffuse { albedo: 0.7 }];
    let mut add = |name: &str, solid: Solid| scene.objects.push(Object::new(name, solid, 0));
    add(
        "oak",
        Solid::box_between([-1.5, -1.5, -0.08], [1.5, 1.8, 0.0]),
    );
    add(
        "wall",
        Solid::box_between([-1.5, 0.95, 0.0], [1.5, 1.0, 1.8]),
    );
    // Hollow body, rounded lip, and loop handle; metres throughout.
    let body = Solid::cone([0., 0., 0.012], [0., 0., 0.117], 0.041, 0.050).minus(Solid::cone(
        [0., 0., 0.022],
        [0., 0., 0.14],
        0.036,
        0.047,
    ));
    let lip = Solid::torus([0., 0., 0.117], 0.048, 0.0028);
    let handle = Solid::torus([0., 0., 0.], 0.029, 0.007)
        .rotated([1., 0., 0.], PI / 2.)
        .translated([0.060, 0., 0.065]);
    add("ceramic", Solid::union(vec![body, lip, handle]));
    add("coffee", Solid::cylinder_z([0., 0., 0.102], 0.0455, 0.001));
    add(
        "ceramic",
        Solid::union(vec![
            Solid::cylinder_z([0., 0., 0.006], 0.086, 0.008),
            Solid::torus([0., 0., 0.010], 0.081, 0.004),
        ]),
    );
    // Flattened ellipsoid as a spoon bowl, with a slightly raised handle.
    let spoon = Solid::sphere([0., 0., 0.], 1.)
        .stretched([0.015, 0.022, 0.004])
        .translated([0.115, -0.045, 0.008])
        .minus(
            Solid::sphere([0., 0., 0.], 1.)
                .stretched([0.014, 0.021, 0.004])
                .translated([0.115, -0.045, 0.010]),
        )
        .plus(Solid::capsule(
            [0.115, -0.06, 0.005],
            [0.128, -0.165, 0.008],
            0.003,
        ));
    add("steel", spoon);
    add(
        "linen",
        // The cloth height has |grad h| < 0.35: sinusoidal terms contribute
        // <0.08 and the Gaussian ridge <0.25. Dividing the max-field by
        // sqrt(1+0.35²) gives a conservative 1-Lipschitz marching field.
        // This is an estimated cloth shape, not a simulated cloth model.
        Solid::custom(
            "folded cloth",
            ([-0.25, -0.13, 0.], [-0.05, 0.16, 0.012]),
            false,
            |p| {
                let dx = p[0] + 0.15;
                let dy = p[1] - 0.015;
                let (s, c) = 0.13_f64.sin_cos();
                let u = c * dx - s * dy;
                let v = s * dx + c * dy;
                let h = 0.003
                    + 0.0009 * (u * 42. + 0.7 * (v * 12.).sin()).sin()
                    + 0.0007 * (v * 37.).sin()
                    + 0.0025
                        * (-((u + 0.055) / 0.009).powi(2)).exp()
                        * (0.65 + 0.35 * (v * 15.).sin());
                ((p[2] - h).abs() - 0.0006)
                    .max(u.abs() - 0.068)
                    .max(v.abs() - 0.115)
                    / (1.0 + 0.35_f64.powi(2)).sqrt()
            },
        ),
    );
    // Small irregular offsets keep the composition from being perfectly tiled.
    for (x, y, a) in [
        (-0.096, -0.098, 0.2),
        (-0.076, -0.113, -0.5),
        (-0.098, -0.122, 0.8),
    ] {
        let bean = Solid::sphere([0., 0., 0.], 1.)
            .stretched([0.006, 0.004, 0.003])
            .minus(Solid::box_between(
                [-0.00045, -0.006, 0.001],
                [0.00045, 0.006, 0.004],
            ))
            .rotated([0., 0., 1.], a)
            .translated([x, y, 0.003]);
        add("bean", bean);
    }
    // Large luminous window, represented as geometry so reflections see it.
    add(
        "window",
        Solid::box_between([-0.65, -0.3, 0.15], [-0.645, 0.65, 1.15]),
    );
    Studio {
        scene,
        march: MarchOptions {
            tolerance: 2e-6,
            normal_epsilon: 2e-7,
            max_steps: 4096,
            t_max: 8.,
            ..MarchOptions::default()
        },
    }
}

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
    let a = hash(i, j) * (1. - u) + hash(i + 1, j) * u;
    let b = hash(i, j + 1) * (1. - u) + hash(i + 1, j + 1) * u;
    a * (1. - v) + b * v
}

fn texture(name: &str, p: Vec3) -> [f64; 3] {
    let [x, y, z] = p.0;
    match name {
        "oak" => {
            // Estimate: procedural grain is appearance, not a wood model.
            let bend = 0.009 * noise(x * 12., y * 3.) + 0.003 * noise(x * 27., y * 11.);
            let grain = noise((x + bend) * 700., y * 20.);
            let pore = noise((x + bend) * 2000., y * 45.).powi(9);
            let broad = noise(x * 40., y * 4.);
            let t =
                0.65 + 0.18 * broad + 0.12 * grain - 0.16 * pore + 0.03 * (160. * (x + bend)).sin();
            [0.38 * t, 0.205 * t, 0.095 * t]
        }
        "linen" => {
            let weave = (x * 12500.).sin() * (y * 12500.).sin();
            let t = 0.85 + 0.09 * weave;
            [0.68 * t, 0.64 * t, 0.53 * t]
        }
        "ceramic" => {
            let t = 0.96 + 0.02 * (x * 6300. + y * 5100. + z * 9700.).sin();
            [0.76 * t, 0.72 * t, 0.63 * t]
        }
        "wall" => [0.64, 0.62, 0.56],
        "coffee" => [0.035, 0.014, 0.006],
        "bean" => [0.075, 0.029, 0.012],
        _ => [0.5; 3],
    }
}

fn multiply(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|k| a[k] * b[k])
}
fn accumulate(sum: &mut [f64; 3], a: [f64; 3]) {
    for k in 0..3 {
        sum[k] += a[k];
    }
}

// Rectangular area light estimator: Le * cos(surface) * cos(light) * A/r².
// Dividing by pi is the Lambert BRDF (exact diffuse closure).
fn window_weight(cs: f64, cl: f64, r2: f64) -> f64 {
    cs.max(0.) * cl.max(0.) * (0.95 * 1.0) / (PI * r2)
}

fn trace(studio: &Studio, mut ray: Ray, rng: &mut Rng, exhausted: &AtomicUsize) -> [f64; 3] {
    let mut sum = [0.; 3];
    let mut throughput = [1.; 3];
    let mut specular = true;
    // Estimate: relative RGB radiance, arbitrary exposure units, daylight-like.
    let window = [5.8, 5.6, 5.2];
    for depth in 0..8 {
        let (hit, ex, _) = studio.scene.intersect(ray, &studio.march);
        if ex {
            exhausted.fetch_add(1, Ordering::Relaxed);
        }
        let Some(h) = hit else {
            // Estimate: weak uniform room/environment illumination.
            accumulate(&mut sum, multiply(throughput, [0.012, 0.015, 0.019]));
            break;
        };
        let name = studio.scene.objects[h.object].name.as_str();
        if name == "window" {
            if specular {
                accumulate(&mut sum, multiply(throughput, window));
            }
            break;
        }
        let p = Vec3(h.hit.point);
        let mut n = Vec3(h.hit.normal);
        if ray.dir.dot(n) > 0. {
            n = -n;
        }
        let origin = p + n * (4. * studio.march.tolerance);
        let cosine = (-ray.dir).dot(n).clamp(0., 1.);
        // Estimates: clear ceramic glaze n=1.5; coffee surface n=1.333.
        // These are ideal smooth coats, not a measured roughness model.
        let coat = match name {
            "ceramic" => fresnel_dielectric(cosine, 1., 1.5),
            "coffee" => fresnel_dielectric(cosine, 1., 1.333),
            "steel" => 1.,
            _ => 0.,
        };
        if rng.uniform() < coat {
            if name == "steel" {
                // Estimate: representative visible conductor n=2.5, k=3.3;
                // constant across RGB, not a steel spectral data set.
                let r = fresnel_conductor(cosine, 2.5, 3.3);
                throughput = throughput.map(|t| t * r);
            }
            ray = Ray::new(origin, ray.dir.reflect(n));
            specular = true;
        } else {
            let albedo = texture(name, p);
            throughput = multiply(throughput, albedo);
            let light = Vec3([-0.6449, -0.3 + 0.95 * rng.uniform(), 0.15 + rng.uniform()]);
            let delta = light - origin;
            let r2 = delta.dot(delta);
            let dir = delta.unit();
            let weight = window_weight(n.dot(dir), -dir.0[0], r2);
            if weight > 0. {
                let shadow = Ray::new(origin, dir);
                let mut opts = studio.march;
                opts.t_max = r2.sqrt() - 4. * studio.march.tolerance;
                let (block, ex, _) = studio.scene.intersect(shadow, &opts);
                if ex {
                    exhausted.fetch_add(1, Ordering::Relaxed);
                }
                if block.is_none() && !ex {
                    accumulate(&mut sum, multiply(throughput, window.map(|v| v * weight)));
                }
            }
            let dir = around(n, rng.uniform().sqrt(), 2. * PI * rng.uniform());
            ray = Ray::new(origin, dir);
            specular = false;
        }
        if depth >= 3 {
            let survival = throughput
                .iter()
                .copied()
                .fold(0., f64::max)
                .clamp(0.05, 0.95);
            if rng.uniform() > survival {
                break;
            }
            throughput = throughput.map(|t| t / survival);
        }
    }
    sum
}

// Edge-aware display filtering only. Retain the original linear PFM for
// measurements: this filter reduces noise but deliberately introduces bias.
fn develop(image: &Image) -> Image {
    let (w, h) = (image.width, image.height);
    let mut pixels = image.pixels.clone();
    for step in [1_isize, 2] {
        pixels = (0..pixels.len())
            .into_par_iter()
            .map(|i| {
                let centre = pixels[i];
                let luminance = |p: [f64; 3]| 0.2126 * p[0] + 0.7152 * p[1] + 0.0722 * p[2];
                let lum = luminance(centre);
                let sigma = 0.012 + 0.08 * lum;
                let mut sum = [0.; 3];
                let mut total = 0.;
                for dy in -2_isize..=2 {
                    for dx in -2_isize..=2 {
                        let x = i % w;
                        let y = i / w;
                        let xx = (x as isize + dx * step).clamp(0, w as isize - 1) as usize;
                        let yy = (y as isize + dy * step).clamp(0, h as isize - 1) as usize;
                        let sample = pixels[yy * w + xx];
                        let difference = luminance(sample) - lum;
                        let weight = (-0.5 * (dx * dx + dy * dy) as f64
                            - difference * difference / (2. * sigma * sigma))
                            .exp();
                        accumulate(&mut sum, sample.map(|v| v * weight));
                        total += weight;
                    }
                }
                sum.map(|v| v / total)
            })
            .collect();
    }
    Image {
        width: w,
        height: h,
        pixels,
    }
}

fn read_linear(path: &std::path::Path, width: usize, height: usize) -> std::io::Result<Image> {
    let bytes = std::fs::read(path)?;
    let header = format!("PF\n{width} {height}\n-1.0\n");
    if !bytes.starts_with(header.as_bytes()) || bytes.len() != header.len() + width * height * 12 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "unexpected PFM dimensions or format",
        ));
    }
    let mut image = Image::new(width, height);
    for (i, chunk) in bytes[header.len()..].chunks_exact(12).enumerate() {
        let y = height - 1 - i / width;
        image.pixels[y * width + i % width] = std::array::from_fn(|k| {
            f32::from_le_bytes(chunk[k * 4..k * 4 + 4].try_into().unwrap()) as f64
        });
    }
    Ok(image)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let out = PathBuf::from(
        args.first()
            .map(String::as_str)
            .unwrap_or("../analysis/photographic_scene"),
    );
    let width = args
        .get(1)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(900usize);
    let samples = args
        .get(2)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(256usize);
    if width == 0 || samples == 0 {
        return Err("width and samples must be positive".into());
    }
    let height = width * 2 / 3;
    if height == 0 {
        return Err("width must be at least 2".into());
    }
    std::fs::create_dir_all(&out)?;
    let studio = studio();
    let camera = Camera::looking_at(
        [0.28, -0.43, 0.31],
        [-0.017, -0.01, 0.048],
        35_f64.to_radians(),
    );
    // Estimate: 5 mm entrance pupil; focus at the cup, in metres.
    let camera = camera.with_depth_of_field(0.005, camera.subject_distance());
    let exhausted = AtomicUsize::new(0);
    let rows = AtomicUsize::new(0);
    let start = std::time::Instant::now();
    let mut image = Image::new(width, height);
    let develop_only = args.get(3).is_some_and(|s| s == "--develop");
    if develop_only {
        image = read_linear(&out.join("tabletop-linear.pfm"), width, height)?;
    } else {
        image
            .pixels
            .par_chunks_mut(width)
            .enumerate()
            .for_each(|(y, row)| {
                for (x, pixel) in row.iter_mut().enumerate() {
                    let mut sum = [0.; 3];
                    for s in 0..samples {
                        let mut rng = Rng(((y * width + x) as u64)
                            .wrapping_mul(0x9e3779b97f4a7c15)
                            ^ (s as u64).wrapping_mul(0xbf58476d1ce4e5b9)
                            ^ 0x123456789abcdef);
                        let jitter = [rng.uniform(), rng.uniform()];
                        let [u, v] = camera.pixel(x, y, width, height, jitter);
                        let ray = camera.ray(u, v, width as f64 / height as f64, &mut rng);
                        accumulate(&mut sum, trace(&studio, ray, &mut rng, &exhausted));
                    }
                    *pixel = sum.map(|v| v / samples as f64);
                }
                let done = rows.fetch_add(1, Ordering::Relaxed) + 1;
                if done.is_multiple_of(50) {
                    eprintln!(
                        "{done}/{height} rows, {:.1}s",
                        start.elapsed().as_secs_f64()
                    );
                }
            });
        image.write_pfm(out.join("tabletop-linear.pfm"))?;
    }
    // Keep both presentations so filtering can be inspected visually.
    Image {
        width,
        height,
        pixels: image
            .pixels
            .iter()
            .map(|p| p.map(|v| 1. - (-2.2 * v.max(0.)).exp()))
            .collect(),
    }
    .write_png(out.join("tabletop-raw.png"), 1.)?;
    image = develop(&image);
    // Exposure/shoulder are photographic presentation, not transport physics.
    let exposed = image
        .pixels
        .iter()
        .map(|p| p.map(|v| 1. - (-2.2 * v.max(0.)).exp()))
        .collect();
    let display = Image {
        width,
        height,
        pixels: exposed,
    };
    display.write_png(out.join("tabletop.png"), 1.)?;
    let mut metadata = if develop_only {
        serde_json::from_slice(&std::fs::read(out.join("metadata.json"))?)?
    } else {
        serde_json::json!({"width":width,"height":height,"samples_per_pixel":samples,
        "seconds":start.elapsed().as_secs_f64(),"exhausted_queries":exhausted.load(Ordering::Relaxed),
        "renderer":"project CAD sphere tracing with example-local RGB Lambert/Fresnel transport",
        "limits":"estimated RGB reflectances; ideal smooth coats; eight-bounce truncation; no measured scene validation"})
    };
    metadata["display_filter"] =
        serde_json::json!("two-pass bilateral; original linear PFM retained unfiltered");
    std::fs::write(
        out.join("metadata.json"),
        serde_json::to_string_pretty(&metadata)?,
    )?;
    eprintln!("{}", metadata);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rectangular_light_matches_on_axis_solid_angle_integral() {
        // Independent exact irradiance from a parallel rectangle, integrating
        // cos²/r² (area form factor); rectangle centred at distance 1.
        let (a, b, d) = (0.475_f64, 0.5_f64, 1.0_f64);
        let exact = 2. / PI
            * (a / (a * a + d * d).sqrt() * (b / (a * a + d * d).sqrt()).atan()
                + b / (b * b + d * d).sqrt() * (a / (b * b + d * d).sqrt()).atan());
        let n = 200;
        let mut measured = 0.;
        for i in 0..n {
            for j in 0..n {
                let x = -a + 2. * a * (i as f64 + 0.5) / n as f64;
                let y = -b + 2. * b * (j as f64 + 0.5) / n as f64;
                let r2 = d * d + x * x + y * y;
                measured += window_weight(d / r2.sqrt(), d / r2.sqrt(), r2) / (n * n) as f64;
            }
        }
        let error = (measured - exact).abs();
        eprintln!(
            "rectangle irradiance: {measured:.12}, exact {exact:.12}, absolute error {error:e}"
        );
        assert!(error < 2e-6);
    }
    #[test]
    fn textures_conserve_diffuse_energy() {
        for name in ["oak", "linen", "ceramic", "coffee", "bean", "wall"] {
            for i in 0..1000 {
                let p = Vec3([i as f64 * 0.0017, i as f64 * 0.0031, i as f64 * 0.00013]);
                assert!(texture(name, p).iter().all(|v| (0.0..=1.0).contains(v)));
            }
        }
    }
}
