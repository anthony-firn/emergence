//! The integrator: what a ray is worth.
//!
//! Two modes, because they answer different questions.
//!
//! * **Direct lighting** evaluates the lights analytically and shadows
//!   them with one ray each. It is deterministic, converges immediately,
//!   and is exactly right for the laws a render should be checked
//!   against — the cosine law, the inverse-square law, a disc's form
//!   factor. It has no interreflection.
//! * **Path tracing** follows a ray through scattering events, sampling
//!   one wavelength per path so dispersion and thin-film colour come out
//!   right. It converges as `1/sqrt(N)` and conserves energy, which the
//!   furnace test checks exactly.
//!
//! Both use `radiation::optics` for the physics: `planck` for emission,
//! `fresnel_dielectric` and `fresnel_conductor` for reflectance, `refract`
//! for transmission, and `cie_xyz` with `xyz_to_rgb` for colour. Nothing
//! about the light transport is reimplemented here; what is new is that
//! the geometry is an implicit field rather than a triangle list.
//!
//! # Sampling and determinism
//!
//! Each sample's generator is seeded from the pixel index and the sample
//! number, so the image does not depend on how work was split across
//! threads. Stratified jitter within the pixel gives the usual
//! variance reduction without breaking that.

use std::collections::BTreeMap;

use geometry::{add, dot, norm, scale, sub, P3};
use radiation::math::{around, Ray, Rng, Vec3};
use radiation::optics::{cie_xyz, fresnel_conductor, fresnel_dielectric, refract, xyz_to_rgb};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::camera::Camera;
use crate::field::Field;
use crate::image::{Channel, Image, Pass, Render};
use crate::scene::{Material, Scene, SceneHit};
use crate::trace::MarchOptions;

/// Which integrator to use.
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Integrator {
    /// Analytic lights, one shadow ray each, no interreflection.
    Direct,
    /// Spectral path tracing.
    Path,
    /// Surface normal shaded against the view direction: no lights needed,
    /// useful for looking at geometry.
    Normals,
}

/// Render settings.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct RenderOptions {
    pub width: usize,
    pub height: usize,
    /// Samples per pixel. One is enough for `Direct` on a scene with no
    /// aperture, since nothing in it is stochastic.
    pub samples: usize,
    pub integrator: Integrator,
    pub max_depth: usize,
    /// Wavelength range sampled by the path integrator, metres.
    pub spectrum: (f64, f64),
    pub passes: Vec<Pass>,
    /// Base seed, so a sequence of frames can be decorrelated.
    pub seed: u64,
}

impl Default for RenderOptions {
    fn default() -> Self {
        Self {
            width: 320,
            height: 240,
            samples: 16,
            integrator: Integrator::Path,
            max_depth: 8,
            spectrum: (380e-9, 780e-9),
            passes: Vec::new(),
            seed: 0x9E3779B97F4A7C15,
        }
    }
}

impl RenderOptions {
    pub fn preview(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            samples: 1,
            integrator: Integrator::Direct,
            ..Self::default()
        }
    }
}

/// Deterministic per-sample seed: mixes the pixel and sample indices so
/// the result is independent of thread scheduling.
fn seed_for(base: u64, pixel: usize, sample: usize) -> Rng {
    let mut x = base
        ^ (pixel as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
        ^ (sample as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    // Two rounds of splitmix to decorrelate neighbouring pixels.
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    Rng(x ^ (x >> 31) | 1)
}

/// Radiance along a ray, in the grey (wavelength-independent) sense that
/// `Direct` works in.
fn direct_radiance(scene: &Scene, ray: Ray, opts: &MarchOptions) -> (f64, Option<SceneHit>, bool, usize) {
    let (hit, exhausted, steps) = scene.intersect(ray, opts);
    let Some(h) = hit else {
        return (scene.background, None, exhausted, steps);
    };
    let Some(material) = scene.material_of(&h) else {
        return (0.0, Some(h), exhausted, steps);
    };
    // Face the normal towards the ray.
    let n = if dot(h.hit.normal, ray.dir.0) > 0.0 {
        scale(h.hit.normal, -1.0)
    } else {
        h.hit.normal
    };
    let mut radiance = match material {
        Material::Emitter { radiance } => *radiance,
        Material::Blackbody { .. } => material.emission(550e-9),
        _ => 0.0,
    };
    let albedo = material.albedo_at(h.hit.point);
    if albedo > 0.0 {
        let offset = add(h.hit.point, scale(n, 4.0 * opts.tolerance));
        for light in &scene.lights {
            let (irradiance, to_light, distance) = light.sample(offset, n);
            if irradiance <= 0.0 {
                continue;
            }
            if scene.occluded(offset, to_light, distance * (1.0 - 1e-6), opts) {
                continue;
            }
            // Lambertian BRDF: rho/pi.
            radiance += albedo / std::f64::consts::PI * irradiance;
        }
    }
    (radiance, Some(h), exhausted, steps)
}

/// Spectral radiance along a ray at one wavelength.
fn path_radiance(
    scene: &Scene,
    mut ray: Ray,
    lambda: f64,
    rng: &mut Rng,
    opts: &MarchOptions,
    max_depth: usize,
) -> (f64, Option<SceneHit>, bool, usize) {
    let mut throughput = 1.0;
    let mut radiance = 0.0;
    let mut first: Option<SceneHit> = None;
    let mut exhausted = false;
    let mut steps = 0usize;
    for depth in 0..max_depth {
        let (hit, ex, st) = scene.intersect(ray, opts);
        exhausted |= ex;
        steps += st;
        let Some(h) = hit else {
            radiance += throughput * scene.background;
            break;
        };
        if depth == 0 {
            first = Some(h);
        }
        let Some(material) = scene.material_of(&h) else {
            break;
        };
        let front = dot(h.hit.normal, ray.dir.0) < 0.0;
        let n = if front { h.hit.normal } else { scale(h.hit.normal, -1.0) };
        radiance += throughput * material.emission(lambda);

        match *material {
            Material::Diffuse { .. } | Material::Checker { .. } => {
                let albedo = material.albedo_at(h.hit.point);
                if albedo <= 0.0 {
                    break;
                }
                throughput *= albedo;
                // Cosine-weighted hemisphere: the pdf cancels the cosine
                // and the 1/pi exactly, so the estimator is unbiased with
                // no leftover factor.
                let cos = rng.uniform().sqrt();
                let dir = around(Vec3(n), cos, std::f64::consts::TAU * rng.uniform());
                ray = Ray::new(
                    Vec3(add(h.hit.point, scale(n, 4.0 * opts.tolerance))),
                    dir,
                );
            }
            Material::Conductor { index, extinction } => {
                let cos = dot(scale(ray.dir.0, -1.0), n).abs();
                throughput *= fresnel_conductor(cos, index, extinction);
                let dir = ray.dir.reflect(Vec3(n));
                ray = Ray::new(
                    Vec3(add(h.hit.point, scale(n, 4.0 * opts.tolerance))),
                    dir,
                );
            }
            Material::Mirror => {
                let dir = ray.dir.reflect(Vec3(n));
                ray = Ray::new(
                    Vec3(add(h.hit.point, scale(n, 4.0 * opts.tolerance))),
                    dir,
                );
            }
            Material::Dielectric { .. } | Material::Dispersive { .. } => {
                let index = material.refractive_index(lambda).unwrap_or(1.5);
                let (ni, nt) = if front { (1.0, index) } else { (index, 1.0) };
                let cos = dot(scale(ray.dir.0, -1.0), n).abs();
                let r = fresnel_dielectric(cos, ni, nt);
                let reflect = rng.uniform() < r;
                let dir = if reflect {
                    ray.dir.reflect(Vec3(n))
                } else {
                    match refract(ray.dir, Vec3(n), ni, nt) {
                        Some(d) => d,
                        None => ray.dir.reflect(Vec3(n)),
                    }
                };
                // Russian-roulette choice between the two lobes is already
                // weighted by R, so the throughput is unchanged.
                let side = if dot(dir.0, n) > 0.0 { 1.0 } else { -1.0 };
                ray = Ray::new(
                    Vec3(add(h.hit.point, scale(n, side * 4.0 * opts.tolerance))),
                    dir,
                );
            }
            Material::Emitter { .. } | Material::Blackbody { .. } => break,
        }
        if throughput < 1e-6 {
            break;
        }
    }
    (radiance, first, exhausted, steps)
}

/// Render a scene.
pub fn render(scene: &Scene, camera: &Camera, opts: &RenderOptions) -> Render {
    render_with_field(scene, camera, opts, None)
}

/// Render a scene, optionally sampling a field at each hit for the
/// [`Pass::Field`] channel.
pub fn render_with_field(
    scene: &Scene,
    camera: &Camera,
    opts: &RenderOptions,
    field: Option<&Field>,
) -> Render {
    let (w, h) = (opts.width.max(1), opts.height.max(1));
    let aspect = w as f64 / h as f64;
    let march = MarchOptions::for_scale(scene.scale());
    let samples = opts.samples.max(1);
    let wants = |p: Pass| opts.passes.contains(&p);

    struct Pixel {
        colour: [f64; 3],
        depth: f64,
        normal: [f64; 3],
        position: [f64; 3],
        object: f64,
        steps: f64,
        field: f64,
        exhausted: usize,
    }

    let rows: Vec<Vec<Pixel>> = (0..h)
        .into_par_iter()
        .map(|y| {
            (0..w)
                .map(|x| {
                    let index = y * w + x;
                    let mut colour = [0.0; 3];
                    let mut depth = f64::INFINITY;
                    let mut normal = [0.0; 3];
                    let mut position = [0.0; 3];
                    let mut object = 0.0;
                    let mut steps = 0.0;
                    let mut field_value = f64::NAN;
                    let mut exhausted = 0usize;
                    for s in 0..samples {
                        let mut rng = seed_for(opts.seed, index, s);
                        // Stratified jitter: one sample sits at the pixel
                        // centre so a single-sample render is a clean
                        // point sample rather than a random one.
                        let jitter = if samples == 1 {
                            [0.5, 0.5]
                        } else {
                            let side = (samples as f64).sqrt().ceil() as usize;
                            let (sx, sy) = (s % side, (s / side) % side);
                            [
                                (sx as f64 + rng.uniform()) / side as f64,
                                (sy as f64 + rng.uniform()) / side as f64,
                            ]
                        };
                        let [u, v] = camera.pixel(x, y, w, h, jitter);
                        let ray = camera.ray(u, v, aspect, &mut rng);
                        let (value, hit, ex, st) = match opts.integrator {
                            Integrator::Direct => direct_radiance(scene, ray, &march),
                            Integrator::Normals => {
                                let (hit, ex, st) = scene.intersect(ray, &march);
                                let shade = hit.map_or(0.0, |hh| {
                                    dot(hh.hit.normal, scale(ray.dir.0, -1.0)).abs()
                                });
                                (shade, hit, ex, st)
                            }
                            Integrator::Path => {
                                let lambda = opts.spectrum.0
                                    + (opts.spectrum.1 - opts.spectrum.0) * rng.uniform();
                                let (l, hit, ex, st) =
                                    path_radiance(scene, ray, lambda, &mut rng, &march, opts.max_depth);
                                // One wavelength per path, weighted by the
                                // CIE observer and the sampling width.
                                let xyz = cie_xyz(lambda);
                                let rgb = xyz_to_rgb([
                                    xyz[0] * l * (opts.spectrum.1 - opts.spectrum.0),
                                    xyz[1] * l * (opts.spectrum.1 - opts.spectrum.0),
                                    xyz[2] * l * (opts.spectrum.1 - opts.spectrum.0),
                                ]);
                                for k in 0..3 {
                                    colour[k] += rgb[k] / samples as f64;
                                }
                                (f64::NAN, hit, ex, st)
                            }
                        };
                        if opts.integrator != Integrator::Path {
                            for k in 0..3 {
                                colour[k] += value / samples as f64;
                            }
                        }
                        if ex {
                            exhausted += 1;
                        }
                        steps += st as f64 / samples as f64;
                        if let Some(hh) = hit {
                            if hh.hit.t < depth {
                                depth = hh.hit.t;
                                normal = hh.hit.normal;
                                position = hh.hit.point;
                                object = hh.object as f64 + 1.0;
                                if let Some(f) = field {
                                    field_value = f.at(hh.hit.point);
                                }
                            }
                        }
                    }
                    Pixel {
                        colour,
                        depth,
                        normal,
                        position,
                        object,
                        steps,
                        field: field_value,
                        exhausted,
                    }
                })
                .collect()
        })
        .collect();

    let mut out = Render::new(w, h);
    out.rays = w * h * samples;
    let mut depth_ch = Channel::new(w, h, f64::INFINITY);
    let mut object_ch = Channel::new(w, h, 0.0);
    let mut steps_ch = Channel::new(w, h, 0.0);
    let mut field_ch = Channel::new(w, h, f64::NAN);
    let mut normal_img = Image::new(w, h);
    let mut position_img = Image::new(w, h);
    for (y, row) in rows.iter().enumerate() {
        for (x, p) in row.iter().enumerate() {
            let i = y * w + x;
            out.beauty.pixels[i] = p.colour;
            out.exhausted += p.exhausted;
            depth_ch.values[i] = p.depth;
            object_ch.values[i] = p.object;
            steps_ch.values[i] = p.steps;
            field_ch.values[i] = p.field;
            normal_img.pixels[i] = p.normal;
            position_img.pixels[i] = p.position;
        }
    }
    let mut channels = BTreeMap::new();
    let mut vectors = BTreeMap::new();
    if wants(Pass::Depth) {
        channels.insert(Pass::Depth, depth_ch);
    }
    if wants(Pass::ObjectId) {
        channels.insert(Pass::ObjectId, object_ch);
    }
    if wants(Pass::Steps) {
        channels.insert(Pass::Steps, steps_ch);
    }
    if wants(Pass::Field) {
        channels.insert(Pass::Field, field_ch);
    }
    if wants(Pass::Normal) {
        vectors.insert(Pass::Normal, normal_img);
    }
    if wants(Pass::Position) {
        vectors.insert(Pass::Position, position_img);
    }
    out.channels = channels;
    out.vectors = vectors;
    out
}

/// Irradiance on a surface patch from the scene's lights alone, with
/// shadowing: the quantity a thermal or photovoltaic model wants, in
/// W/m^2.
pub fn irradiance_at(scene: &Scene, point: P3, normal: P3) -> f64 {
    let opts = MarchOptions::for_scale(scene.scale());
    let n = geometry::normalize(normal);
    let offset = add(point, scale(n, 4.0 * opts.tolerance));
    let mut total = 0.0;
    for light in &scene.lights {
        let (e, dir, distance) = light.sample(offset, n);
        if e <= 0.0 {
            continue;
        }
        if scene.occluded(offset, dir, distance * (1.0 - 1e-6), &opts) {
            continue;
        }
        total += e;
    }
    total
}

/// The view direction a camera has at a pixel, without any ray marching:
/// used by the tests to check the projection round trip.
pub fn view_direction(camera: &Camera, u: f64, v: f64, aspect: f64) -> P3 {
    let mut rng = Rng(1);
    let ray = camera.ray(u, v, aspect, &mut rng);
    let d = ray.dir.0;
    let n = norm(d);
    [d[0] / n, d[1] / n, d[2] / n]
}

/// Distance from `a` to `b`, for tests and tools.
pub fn distance(a: P3, b: P3) -> f64 {
    norm(sub(a, b))
}
