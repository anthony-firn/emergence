use crate::geometry::{Bvh, Hit};
use crate::math::{around, Ray, Rng, Vec3};
use crate::medium::{SpectralMedia, Volume};
use crate::optics::{
    cie_xyz, fresnel_conductor, fresnel_dielectric, hg, planck, refract, sample_hg, xyz_to_rgb,
    OpticalMaterial, Spectrum,
};
use crate::sampling::WavelengthSampler;
use rayon::prelude::*;
use std::{
    f64::consts::PI,
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

#[derive(Clone, Debug)]
pub enum Surface {
    Diffuse {
        albedo: Spectrum,
    },
    Checker {
        a: Spectrum,
        b: Spectrum,
        period: f64,
    },
    Dielectric {
        optics: OpticalMaterial,
        outside_n: f64,
    },
    /// Opaque polished conductor; directional emissivity = 1 - Fresnel R.
    Conductor {
        optics: OpticalMaterial,
        temperature: f64,
    },
    Emitter {
        temperature: f64,
        scale: f64,
    },
}
#[derive(Clone, Debug)]
pub struct AreaLight {
    pub origin: Vec3,
    pub u: Vec3,
    pub v: Vec3,
    pub temperature: f64,
    pub scale: f64,
    pub material: usize,
}
impl AreaLight {
    fn area(&self) -> f64 {
        self.u.cross(self.v).length()
    }
    fn radiance(&self, l: f64) -> f64 {
        planck(l, self.temperature) * self.scale
    }
}
pub struct Scene {
    pub geometry: Bvh,
    pub materials: Vec<Surface>,
    pub volume: Option<Volume>,
    pub lights: Vec<AreaLight>,
    pub environment_temperature: f64,
    pub environment_scale: f64,
    /// SI distance, chosen relative to scene scale; introduces this tolerance.
    pub epsilon: f64,
}
#[derive(Clone, Copy, Debug)]
pub struct Camera {
    origin: Vec3,
    forward: Vec3,
    right: Vec3,
    up: Vec3,
    tan_half_fov: f64,
}
impl Camera {
    pub fn new(origin: Vec3, target: Vec3, up: Vec3, vertical_fov_degrees: f64) -> Self {
        let forward = (target - origin).unit();
        let right = forward.cross(up).unit();
        let up = right.cross(forward);
        assert!(vertical_fov_degrees > 0.0 && vertical_fov_degrees < 180.0);
        Self {
            origin,
            forward,
            right,
            up,
            tan_half_fov: (vertical_fov_degrees * PI / 360.0).tan(),
        }
    }
    pub fn ray(&self, x: f64, y: f64, aspect: f64) -> Ray {
        Ray::new(
            self.origin,
            self.forward
                + self.right * ((2.0 * x - 1.0) * aspect * self.tan_half_fov)
                + self.up * ((1.0 - 2.0 * y) * self.tan_half_fov),
        )
    }
}
#[derive(Clone, Copy, Debug)]
pub struct RenderOptions {
    pub width: usize,
    pub height: usize,
    pub samples: usize,
    pub max_depth: usize,
    pub seed: u64,
}
impl Default for RenderOptions {
    fn default() -> Self {
        Self {
            width: 640,
            height: 400,
            samples: 256,
            max_depth: 128,
            seed: 42,
        }
    }
}
/// Pilot samples are discarded; independent main samples avoid optional-stopping bias.
#[derive(Clone, Copy, Debug)]
pub struct AdaptiveOptions {
    pub pilot_samples: usize,
    pub min_samples: usize,
    pub relative_error: f64,
    pub absolute_error: f64,
    /// Absolute upper bound on omitted thermal spectral radiance PER history.
    /// Units W m^-2 sr^-1 m^-1. Set zero for exact volume emission.
    pub thermal_radiance_budget: f64,
}
impl Default for AdaptiveOptions {
    fn default() -> Self {
        Self {
            pilot_samples: 64,
            min_samples: 64,
            relative_error: 0.03,
            absolute_error: 0.02,
            thermal_radiance_budget: 1e-6,
        }
    }
}
struct Pixel {
    rgb: [f64; 3],
    mean_y: f64,
    error: f64,
    truncated: usize,
}
pub struct Image {
    pub width: usize,
    pub height: usize,
    pub rgb: Vec<[f64; 3]>,
    pub standard_error_y: Vec<f64>,
    pub truncated_paths: usize,
    pub samples_per_pixel: Vec<usize>,
    pub pilot_paths: usize,
}
impl Image {
    /// Linear RGB, unexposed, lossless float output for scientific reuse.
    pub fn write_pfm(&self, path: impl AsRef<Path>) -> io::Result<()> {
        let mut f = BufWriter::new(File::create(path)?);
        write!(f, "PF\n{} {}\n-1.0\n", self.width, self.height)?;
        for y in (0..self.height).rev() {
            for x in 0..self.width {
                for v in self.rgb[y * self.width + x] {
                    f.write_all(&(v as f32).to_le_bytes())?;
                }
            }
        }
        Ok(())
    }
    /// Display transform only: exposure, Reinhard, sRGB. The PFM retains HDR.
    pub fn write_png(
        &self,
        path: impl AsRef<Path>,
        exposure: f64,
    ) -> Result<(), Box<dyn std::error::Error>> {
        assert!(exposure > 0.0 && exposure.is_finite());
        let mut encoder = png::Encoder::new(
            BufWriter::new(File::create(path)?),
            self.width as u32,
            self.height as u32,
        );
        encoder.set_color(png::ColorType::Rgb);
        encoder.set_depth(png::BitDepth::Eight);
        encoder.set_source_srgb(png::SrgbRenderingIntent::Perceptual);
        let data: Vec<u8> = self
            .rgb
            .iter()
            .flat_map(|rgb| {
                rgb.map(|v| {
                    let a = (v * exposure).max(0.0);
                    let linear = a / (1.0 + a);
                    let srgb = if linear <= 0.0031308 {
                        12.92 * linear
                    } else {
                        1.055 * linear.powf(1.0 / 2.4) - 0.055
                    };
                    (255.0 * srgb.clamp(0.0, 1.0) + 0.5) as u8
                })
            })
            .collect();
        encoder.write_header()?.write_image_data(&data)?;
        Ok(())
    }
}
impl Scene {
    fn shadow_transmittance(
        &self,
        p: Vec3,
        dir: Vec3,
        distance: f64,
        media: &Option<SpectralMedia>,
    ) -> f64 {
        let r = Ray::new(p + dir * self.epsilon, dir);
        if self
            .geometry
            .occluded(r, (distance - 2.0 * self.epsilon).max(0.0))
        {
            return 0.0;
        }
        self.volume.as_ref().map_or(1.0, |v| {
            v.transmittance_cached(r, distance, media.as_ref().unwrap())
        })
    }
    /// Next-event estimator for area sources. Their matching geometry must
    /// have an Emitter surface with the same spectrum. Emission after a
    /// diffuse/volume continuation is excluded, avoiding double counting.
    fn direct(
        &self,
        p: Vec3,
        l: f64,
        rng: &mut Rng,
        media: &Option<SpectralMedia>,
        weight: impl Fn(Vec3) -> f64,
    ) -> f64 {
        if self.lights.is_empty() {
            return 0.0;
        }
        let light = &self.lights[(rng.uniform() * self.lights.len() as f64) as usize];
        let q = light.origin + light.u * rng.uniform() + light.v * rng.uniform();
        let delta = q - p;
        let distance = delta.length();
        let dir = delta / distance;
        let cos = light.u.cross(light.v).unit().dot(-dir).abs();
        light.radiance(l) * weight(dir) * cos * light.area() * self.lights.len() as f64
            / (distance * distance)
            * self.shadow_transmittance(p, dir, distance, media)
    }
    /// Spectral radiance estimator at one vacuum wavelength. Absorption is
    /// integrated exactly along each voxel segment, scattering is sampled.
    pub fn radiance(&self, ray: Ray, l: f64, rng: &mut Rng, max_depth: usize) -> (f64, bool) {
        self.radiance_prepared(ray, l, rng, max_depth, 0.0, 0.0)
    }
    pub fn radiance_with_budget(
        &self,
        ray: Ray,
        l: f64,
        rng: &mut Rng,
        max_depth: usize,
        budget: f64,
    ) -> (f64, bool) {
        self.radiance_prepared(ray, l, rng, max_depth, self.max_temperature(), budget)
    }
    fn max_temperature(&self) -> f64 {
        self.volume.as_ref().map_or(0.0, |v| {
            v.cells.iter().map(|v| v.temperature).fold(0.0, f64::max)
        })
    }
    fn radiance_prepared(
        &self,
        mut ray: Ray,
        l: f64,
        rng: &mut Rng,
        max_depth: usize,
        max_temperature: f64,
        budget: f64,
    ) -> (f64, bool) {
        assert!(budget.is_finite() && budget >= 0.0);
        let media = self
            .volume
            .as_ref()
            .map(|v| SpectralMedia::new(&v.media, l));
        let source_bound = if budget > 0.0 {
            planck(l, max_temperature)
        } else {
            0.0
        };
        let mut remaining_budget = budget;
        let mut beta = 1.0;
        let mut result = 0.0;
        let mut previous_delta = true;
        for depth in 0..max_depth {
            let hit = self.geometry.hit(ray, f64::INFINITY);
            let limit = hit.map_or(f64::INFINITY, |h| h.distance);
            let mut scatter = None;
            if let Some(volume) = &self.volume {
                let mut tau = -rng.uniform().ln();
                for seg in volume.segments(ray, limit) {
                    let voxel = volume.cells[seg.cell];
                    let c = media.as_ref().unwrap().get(voxel.medium);
                    let (a, s) = (c.a * voxel.scale, c.s * voxel.scale);
                    let ds = seg.end - seg.start;
                    let collision = s > 0.0 && s * ds > tau;
                    let length = if collision { tau / s } else { ds };
                    let optical_depth = a * length;
                    if optical_depth > 0.0 {
                        // 1-exp(-x) <= min(x,1), and B(T) <= B(T_max).
                        // Track the bound including the current Monte Carlo weight:
                        // all skipped contributions together remain below the budget.
                        let bound = beta * optical_depth.min(1.0) * source_bound * c.n2;
                        if budget > 0.0 && bound <= remaining_budget {
                            remaining_budget -= bound;
                        } else {
                            result += beta
                                * (-(-optical_depth).exp_m1())
                                * planck(l, voxel.temperature)
                                * c.n2;
                        }
                    }
                    beta *= (-a * length).exp();
                    if collision {
                        scatter = Some((ray.at(seg.start + length), c.g));
                        break;
                    }
                    tau -= s * ds;
                }
            }
            if beta == 0.0 {
                return (result, false);
            }
            if let Some((p, g)) = scatter {
                result += beta * self.direct(p, l, rng, &media, |wi| hg(ray.dir.dot(wi), g));
                ray = Ray::new(p, sample_hg(ray.dir, g, rng));
                previous_delta = false;
            } else if let Some(h) = hit {
                let front = ray.dir.dot(h.normal) < 0.0;
                let normal = if front { h.normal } else { -h.normal };
                match &self.materials[h.material] {
                    Surface::Emitter { temperature, scale } => {
                        if previous_delta
                            || !self.lights.iter().any(|light| light.material == h.material)
                        {
                            result += beta * planck(l, *temperature) * scale;
                        }
                        return (result, false);
                    }
                    Surface::Dielectric { optics, outside_n } => {
                        let inside = optics.n.at(l);
                        let (ni, nt) = if front {
                            (*outside_n, inside)
                        } else {
                            (inside, *outside_n)
                        };
                        let cos = -ray.dir.dot(normal);
                        let f = fresnel_dielectric(cos, ni, nt);
                        let dir = if rng.uniform() < f {
                            ray.dir.reflect(normal)
                        } else {
                            let eta = ni / nt;
                            beta *= eta * eta;
                            refract(ray.dir, normal, ni, nt)
                                .expect("Fresnel excludes total internal transmission")
                        };
                        ray = Ray::new(h.point + dir * self.epsilon, dir);
                        previous_delta = true;
                    }
                    Surface::Conductor {
                        optics,
                        temperature,
                    } => {
                        let f =
                            fresnel_conductor(-ray.dir.dot(normal), optics.n.at(l), optics.k.at(l));
                        result += beta * (1.0 - f) * planck(l, *temperature);
                        beta *= f;
                        let dir = ray.dir.reflect(normal);
                        ray = Ray::new(h.point + dir * self.epsilon, dir);
                        previous_delta = true;
                    }
                    Surface::Diffuse { albedo } => {
                        let a = albedo.at(l);
                        self.diffuse(
                            &mut ray,
                            h,
                            normal,
                            l,
                            a,
                            rng,
                            &media,
                            &mut beta,
                            &mut result,
                        );
                        previous_delta = false;
                    }
                    Surface::Checker { a, b, period } => {
                        let even = ((h.point.0[0] / period).floor() as i64
                            + (h.point.0[2] / period).floor() as i64)
                            .rem_euclid(2)
                            == 0;
                        let a = if even { a.at(l) } else { b.at(l) };
                        self.diffuse(
                            &mut ray,
                            h,
                            normal,
                            l,
                            a,
                            rng,
                            &media,
                            &mut beta,
                            &mut result,
                        );
                        previous_delta = false;
                    }
                }
            } else {
                result += beta * planck(l, self.environment_temperature) * self.environment_scale;
                return (result, false);
            }
            if depth >= 6 {
                let survival = beta.clamp(0.05, 0.95);
                if rng.uniform() > survival {
                    return (result, false);
                }
                beta /= survival;
            }
        }
        (result, true)
    }
    #[allow(clippy::too_many_arguments)]
    fn diffuse(
        &self,
        ray: &mut Ray,
        h: Hit,
        n: Vec3,
        l: f64,
        a: f64,
        rng: &mut Rng,
        media: &Option<SpectralMedia>,
        beta: &mut f64,
        result: &mut f64,
    ) {
        assert!((0.0..=1.0).contains(&a));
        *result += *beta * self.direct(h.point, l, rng, media, |wi| a / PI * n.dot(wi).max(0.0));
        *beta *= a;
        let dir = around(n, rng.uniform().sqrt(), 2.0 * PI * rng.uniform());
        *ray = Ray::new(h.point + dir * self.epsilon, dir);
    }
    fn validate_options(&self, opts: RenderOptions) {
        assert!(
            opts.width > 0
                && opts.height > 0
                && opts.samples >= 16
                && opts.samples.is_multiple_of(16)
                && opts.max_depth > 0
        );
        assert!(self.epsilon > 0.0 && self.epsilon.is_finite());
    }
    fn pixel(
        &self,
        camera: &Camera,
        opts: RenderOptions,
        q: usize,
        sampler: &WavelengthSampler,
        max_t: f64,
        budget: f64,
    ) -> Pixel {
        let mut rng = Rng(opts
            .seed
            .wrapping_add((q as u64).wrapping_mul(0xd1342543de82ef95)));
        let mut xyz = [0.0; 3];
        let mut truncated = 0;
        let mut mean = 0.0;
        let mut m2 = 0.0;
        let blocks = opts.samples / 4;
        for block in 0..blocks {
            let mut block_y = 0.0;
            let path_seed = rng.next_u64();
            let x = (q % opts.width) as f64 + rng.uniform();
            let y = (q / opts.width) as f64 + rng.uniform();
            let r = camera.ray(
                x / opts.width as f64,
                y / opts.height as f64,
                opts.width as f64 / opts.height as f64,
            );
            for stratum in 0..4 {
                let (l, pdf) = sampler.sample((stratum as f64 + rng.uniform()) / 4.0);
                let (value, cap) = self.radiance_prepared(
                    r,
                    l,
                    &mut Rng(path_seed),
                    opts.max_depth,
                    max_t,
                    budget,
                );
                truncated += cap as usize;
                let cmf = cie_xyz(l);
                for d in 0..3 {
                    xyz[d] += value * cmf[d] / (pdf * opts.samples as f64);
                }
                block_y += value * cmf[1] / (pdf * 4.0);
            }
            let delta = block_y - mean;
            mean += delta / (block + 1) as f64;
            m2 += delta * (block_y - mean);
        }
        Pixel {
            rgb: xyz_to_rgb(xyz),
            mean_y: mean,
            error: (m2 / ((blocks - 1) * blocks) as f64).sqrt(),
            truncated,
        }
    }
    fn collect_image(
        &self,
        opts: RenderOptions,
        pixels: Vec<(Pixel, usize)>,
        pilot_paths: usize,
    ) -> Image {
        Image {
            width: opts.width,
            height: opts.height,
            rgb: pixels.iter().map(|p| p.0.rgb).collect(),
            standard_error_y: pixels.iter().map(|p| p.0.error).collect(),
            truncated_paths: pixels.iter().map(|p| p.0.truncated).sum(),
            samples_per_pixel: pixels.iter().map(|p| p.1).collect(),
            pilot_paths,
        }
    }
    /// Full-support CIE importance sampling, correlated dispersive wavelength
    /// bundles, independent pixel streams. Thermal omission is bounded by
    /// 1e-6 spectral radiance per history; use render_with_budget(...,0) for exact.
    pub fn render(&self, camera: &Camera, opts: RenderOptions) -> Image {
        self.render_with_budget(camera, opts, 1e-6)
    }
    pub fn render_with_budget(&self, camera: &Camera, opts: RenderOptions, budget: f64) -> Image {
        self.validate_options(opts);
        assert!(budget >= 0.0 && budget.is_finite());
        let sampler = WavelengthSampler::default();
        let max_t = self.max_temperature();
        let pixels = (0..opts.width * opts.height)
            .into_par_iter()
            .map(|q| {
                (
                    self.pixel(camera, opts, q, &sampler, max_t, budget),
                    opts.samples,
                )
            })
            .collect();
        self.collect_image(opts, pixels, 0)
    }
    /// Discard an independent pilot, select a fixed main budget per pixel,
    /// then estimate using only the main samples. Noise targets are estimates,
    /// not guarantees for rare caustics that the pilot may never encounter.
    pub fn render_adaptive(
        &self,
        camera: &Camera,
        opts: RenderOptions,
        a: AdaptiveOptions,
    ) -> Image {
        self.validate_options(opts);
        assert!(
            a.pilot_samples >= 16
                && a.pilot_samples.is_multiple_of(16)
                && a.min_samples >= 16
                && a.min_samples.is_multiple_of(16)
                && a.min_samples <= opts.samples
        );
        assert!(
            a.relative_error >= 0.0
                && a.relative_error.is_finite()
                && a.absolute_error > 0.0
                && a.absolute_error.is_finite()
        );
        assert!(a.thermal_radiance_budget >= 0.0 && a.thermal_radiance_budget.is_finite());
        let sampler = WavelengthSampler::default();
        let max_t = self.max_temperature();
        let pixels = (0..opts.width * opts.height)
            .into_par_iter()
            .map(|q| {
                let pilot = self.pixel(
                    camera,
                    RenderOptions {
                        samples: a.pilot_samples,
                        seed: opts.seed ^ 0xa0761d6478bd642f,
                        ..opts
                    },
                    q,
                    &sampler,
                    max_t,
                    a.thermal_radiance_budget,
                );
                let target = a.absolute_error + a.relative_error * pilot.mean_y.abs();
                let need = (pilot.error / target).powi(2) * a.pilot_samples as f64;
                let mut samples =
                    ((need.ceil() as usize).clamp(a.min_samples, opts.samples).div_ceil(16) * 16)
                        .min(opts.samples);
                if pilot.truncated > 0 {
                    samples = opts.samples;
                }
                (
                    self.pixel(
                        camera,
                        RenderOptions { samples, ..opts },
                        q,
                        &sampler,
                        max_t,
                        a.thermal_radiance_budget,
                    ),
                    samples,
                )
            })
            .collect();
        self.collect_image(opts, pixels, opts.width * opts.height * a.pilot_samples)
    }
}
