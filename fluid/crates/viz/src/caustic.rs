//! Caustics by light tracing, and why they are the renderer's best
//! accuracy test.
//!
//! A caustic is light that specular geometry has concentrated: the bright
//! cusp under a wine glass, the bands on a pool floor, the rainbow. An
//! eye-path tracer renders them badly for a structural reason — the path
//! eye → diffuse → specular → light has measure zero under light sampling,
//! so the estimator never finds it — and no amount of sampling fixes that.
//! The cure is to go the other way: trace photons **from** the light,
//! follow them through the specular chain, and deposit their flux where
//! they land (Arvo 1986; Jensen, *Realistic Image Synthesis Using Photon
//! Mapping*, 2001).
//!
//! That makes caustics the most demanding thing a renderer draws, and
//! therefore the most useful thing to check it with. A caustic's position
//! is set by refraction and reflection *exactly*, so it has closed-form
//! references that a diffuse render never does:
//!
//! | caustic | exact reference |
//! |---|---|
//! | a sphere, one internal reflection | the rainbow angle, `cos θ₁ = √((n²−1)/3)` |
//! | a sphere, two internal reflections | the secondary rainbow, `cos θ₁ = √((n²−1)/8)` |
//! | a spherical mirror | the catacaustic cusp at `R/2`, and its nephroid envelope |
//! | a paraboloid | a point, with no caustic at all |
//! | any of them | the flux that arrives equals the flux that left |
//!
//! The tests use every one of those, which is how the specular transport —
//! Snell's law, Fresnel's coefficients, total internal reflection, the
//! march that finds the interfaces — is checked as a whole rather than
//! one formula at a time.
//!
//! # Where geometrical optics gives up
//!
//! At a caustic, ray optics predicts **infinite** irradiance: the
//! Jacobian of the ray map vanishes, and the flux of a tube of rays is
//! squeezed into zero area. That is not a numerical difficulty but the
//! model breaking, and it is the sharpest possible demonstration of the
//! ladder in [`crate::wave`]. The finite answer is the Airy function, and
//! [`fold`] gives it: the caustic's peak is proportional to `λ^(-1/3)` in
//! irradiance, and the fringes beside it to `λ^(2/3)` in width. A photon
//! map's caustic peak rises without limit as the kernel narrows; the wave
//! peak does not.

use geometry::{add, cross, dot, norm, normalize, scale, sub, P3};
use radiation::math::{Ray, Rng, Vec3};
use radiation::optics::{fresnel_conductor, fresnel_dielectric, refract};
use serde::{Deserialize, Serialize};

use crate::scene::{Material, Scene};
use crate::trace::MarchOptions;

/// Where photons are launched from.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
#[serde(tag = "emitter", rename_all = "snake_case")]
pub enum Emitter {
    /// A collimated beam: a disc of radius `radius` at `origin`, all
    /// photons travelling along `direction`, carrying `power` watts in
    /// total.
    ///
    /// This is the exact-reference case. The rainbow angle and the
    /// catacaustic are both defined for parallel light, and a beam has no
    /// solid-angle sampling to get wrong.
    Beam {
        origin: P3,
        direction: P3,
        radius: f64,
        power: f64,
    },
    /// An isotropic point source of `power` watts, with photons aimed into
    /// the cone that subtends `target_radius` at `target`.
    ///
    /// Aiming is an importance decision, not an approximation: the pdf of
    /// the cone is divided out, so the estimate is unbiased, and photons
    /// that would have missed the optic are never traced.
    Point {
        position: P3,
        power: f64,
        target: P3,
        target_radius: f64,
    },
}

impl Emitter {
    pub fn power(&self) -> f64 {
        match *self {
            Emitter::Beam { power, .. } | Emitter::Point { power, .. } => power,
        }
    }

    /// One photon: a ray and the flux it carries.
    pub fn emit(&self, rng: &mut Rng, count: usize) -> (Ray, f64) {
        match *self {
            Emitter::Beam {
                origin,
                direction,
                radius,
                power,
            } => {
                let d = normalize(direction);
                let (u, v) = basis(d);
                let (r, phi) = (radius * rng.uniform().sqrt(), std::f64::consts::TAU * rng.uniform());
                let p = add(origin, add(scale(u, r * phi.cos()), scale(v, r * phi.sin())));
                (Ray::new(Vec3(p), Vec3(d)), power / count as f64)
            }
            Emitter::Point {
                position,
                power,
                target,
                target_radius,
            } => {
                let axis = sub(target, position);
                let distance = norm(axis);
                let axis = scale(axis, 1.0 / distance);
                // Cone that just contains the target sphere.
                let sin_max = (target_radius / distance).min(1.0);
                let cos_max = (1.0 - sin_max * sin_max).max(0.0).sqrt();
                let cos = 1.0 - rng.uniform() * (1.0 - cos_max);
                let dir = radiation::math::around(
                    Vec3(axis),
                    cos,
                    std::f64::consts::TAU * rng.uniform(),
                );
                // Fraction of the full sphere the cone covers: the pdf
                // correction that keeps the estimate unbiased.
                let fraction = 0.5 * (1.0 - cos_max);
                (
                    Ray::new(Vec3(position), dir),
                    power * fraction / count as f64,
                )
            }
        }
    }

    /// Solid-angle fraction of the full sphere that photons are aimed
    /// into; one for a beam.
    pub fn aimed_fraction(&self) -> f64 {
        match *self {
            Emitter::Beam { .. } => 1.0,
            Emitter::Point {
                position,
                target,
                target_radius,
                ..
            } => {
                let distance = norm(sub(target, position));
                let sin_max = (target_radius / distance).min(1.0);
                let cos_max = (1.0 - sin_max * sin_max).max(0.0).sqrt();
                0.5 * (1.0 - cos_max)
            }
        }
    }
}

fn basis(d: P3) -> (P3, P3) {
    let k = (0..3).min_by(|&i, &j| d[i].abs().partial_cmp(&d[j].abs()).unwrap()).unwrap();
    let mut e = [0.0; 3];
    e[k] = 1.0;
    let u = normalize(cross(d, e));
    (u, cross(d, u))
}

/// A plane that photons land on, with a flux estimate per bin.
///
/// Deposition is a bilinear (tent) splat, whose weights sum to one, so the
/// total flux on the receiver is **exactly** the flux of the photons that
/// hit it, whatever the bin size. That is what makes the energy test a
/// test of the transport rather than of the estimator.
#[derive(Clone, Debug)]
pub struct Receiver {
    pub origin: P3,
    /// In-plane axes; their lengths set the extent.
    pub u: P3,
    pub v: P3,
    pub bins: [usize; 2],
    /// Flux in each bin, watts.
    pub flux: Vec<f64>,
    /// CIE XYZ accumulated per bin, for the colour of a dispersive
    /// caustic.
    pub xyz: Vec<[f64; 3]>,
    pub photons: u64,
    /// Flux that landed on the plane but outside the binned extent.
    pub spilled: f64,
}

impl Receiver {
    pub fn new(origin: P3, u: P3, v: P3, bins: [usize; 2]) -> Self {
        let n = bins[0].max(1) * bins[1].max(1);
        Self {
            origin,
            u,
            v,
            bins: [bins[0].max(1), bins[1].max(1)],
            flux: vec![0.0; n],
            xyz: vec![[0.0; 3]; n],
            photons: 0,
            spilled: 0.0,
        }
    }

    /// Unit normal of the plane.
    pub fn normal(&self) -> P3 {
        normalize(cross(self.u, self.v))
    }

    /// Area of one bin.
    pub fn bin_area(&self) -> f64 {
        norm(cross(self.u, self.v)) / (self.bins[0] * self.bins[1]) as f64
    }

    /// Normalised in-plane coordinates of a point, each in `[0, 1]`.
    pub fn coordinates(&self, p: P3) -> [f64; 2] {
        let d = sub(p, self.origin);
        [
            dot(d, self.u) / dot(self.u, self.u),
            dot(d, self.v) / dot(self.v, self.v),
        ]
    }

    /// Where a ray meets the plane, if it does ahead of itself.
    pub fn hit(&self, ray: Ray) -> Option<(f64, P3)> {
        let n = self.normal();
        let denom = dot(ray.dir.0, n);
        if denom.abs() < 1e-300 {
            return None;
        }
        let t = dot(sub(self.origin, ray.origin.0), n) / denom;
        if t <= 1e-12 {
            return None;
        }
        Some((t, ray.at(t).0))
    }

    /// Splat `flux` at `point`, bilinearly. Returns false if it fell
    /// outside the binned extent.
    pub fn deposit(&mut self, point: P3, flux: f64, xyz: [f64; 3]) -> bool {
        let [s, t] = self.coordinates(point);
        if !(0.0..1.0).contains(&s) || !(0.0..1.0).contains(&t) {
            self.spilled += flux;
            return false;
        }
        let (nx, ny) = (self.bins[0], self.bins[1]);
        // Bin centres at (i + 1/2)/n, so a tent between neighbours.
        let (fx, fy) = (s * nx as f64 - 0.5, t * ny as f64 - 0.5);
        let (i0, j0) = (fx.floor(), fy.floor());
        let (ax, ay) = (fx - i0, fy - j0);
        let mut placed = 0.0;
        for (di, wx) in [(0i64, 1.0 - ax), (1, ax)] {
            for (dj, wy) in [(0i64, 1.0 - ay), (1, ay)] {
                let (i, j) = (i0 as i64 + di, j0 as i64 + dj);
                if i < 0 || j < 0 || i >= nx as i64 || j >= ny as i64 {
                    continue;
                }
                let w = wx * wy;
                let k = j as usize * nx + i as usize;
                self.flux[k] += flux * w;
                for c in 0..3 {
                    self.xyz[k][c] += xyz[c] * flux * w;
                }
                placed += w;
            }
        }
        // Weight that fell off the edge of the grid is spill, not loss.
        self.spilled += flux * (1.0 - placed);
        true
    }

    /// Total flux held by the bins.
    pub fn total_flux(&self) -> f64 {
        self.flux.iter().sum()
    }

    /// Irradiance per bin, W/m².
    pub fn irradiance(&self) -> Vec<f64> {
        let a = self.bin_area();
        self.flux.iter().map(|f| f / a).collect()
    }

    /// Brightest bin's irradiance and its indices.
    pub fn peak(&self) -> (f64, [usize; 2]) {
        let a = self.bin_area();
        let mut best = (0.0, [0, 0]);
        for (k, f) in self.flux.iter().enumerate() {
            if *f / a > best.0 {
                best = (*f / a, [k % self.bins[0], k / self.bins[0]]);
            }
        }
        best
    }

    /// Irradiance as a channel for colouring or writing out.
    pub fn channel(&self) -> crate::image::Channel {
        crate::image::Channel {
            width: self.bins[0],
            height: self.bins[1],
            values: self.irradiance(),
        }
    }

    /// The caustic as a linear-RGB image, from the accumulated XYZ.
    pub fn image(&self) -> crate::image::Image {
        let a = self.bin_area();
        let mut img = crate::image::Image::new(self.bins[0], self.bins[1]);
        for (k, xyz) in self.xyz.iter().enumerate() {
            let scaled = [xyz[0] / a, xyz[1] / a, xyz[2] / a];
            img.pixels[k] = radiation::optics::xyz_to_rgb(scaled);
        }
        img
    }
}

/// A histogram of exit directions, binned by deviation from the incident
/// one. This is the receiver a rainbow is measured with.
#[derive(Clone, Debug)]
pub struct Angular {
    pub incident: P3,
    /// Deviation range in radians, `[lo, hi]`.
    pub range: [f64; 2],
    pub flux: Vec<f64>,
    pub photons: u64,
    pub outside: f64,
}

impl Angular {
    pub fn new(incident: P3, range: [f64; 2], bins: usize) -> Self {
        Self {
            incident: normalize(incident),
            range,
            flux: vec![0.0; bins.max(1)],
            photons: 0,
            outside: 0.0,
        }
    }

    /// Deviation of a direction from the incident one, radians.
    pub fn deviation(&self, dir: P3) -> f64 {
        dot(normalize(dir), self.incident).clamp(-1.0, 1.0).acos()
    }

    pub fn deposit(&mut self, dir: P3, flux: f64) {
        let d = self.deviation(dir);
        let n = self.flux.len();
        let f = (d - self.range[0]) / (self.range[1] - self.range[0]);
        if !(0.0..1.0).contains(&f) {
            self.outside += flux;
            return;
        }
        // Flux per unit solid angle needs the sin of the deviation; the
        // raw histogram is kept and divided in `radiance`.
        self.flux[(f * n as f64) as usize] += flux;
    }

    /// Bin centre deviations.
    pub fn centres(&self) -> Vec<f64> {
        let n = self.flux.len();
        (0..n)
            .map(|i| self.range[0] + (self.range[1] - self.range[0]) * (i as f64 + 0.5) / n as f64)
            .collect()
    }

    /// Flux per unit solid angle: the histogram divided by the annulus
    /// `2π sin(θ) Δθ` each bin covers. Without this an isotropic
    /// scatterer would appear to peak at 90°.
    pub fn radiance(&self) -> Vec<f64> {
        let n = self.flux.len();
        let dtheta = (self.range[1] - self.range[0]) / n as f64;
        self.centres()
            .iter()
            .zip(&self.flux)
            .map(|(theta, f)| {
                let solid = 2.0 * std::f64::consts::PI * theta.sin().abs() * dtheta;
                if solid > 0.0 {
                    f / solid
                } else {
                    0.0
                }
            })
            .collect()
    }

    /// Deviation of the brightest bin, by a parabolic fit through its two
    /// neighbours — so the answer is not limited to the bin width.
    pub fn peak_deviation(&self) -> Option<f64> {
        let r = self.radiance();
        let n = r.len();
        let (mut best, mut at) = (0.0, 0usize);
        for (i, v) in r.iter().enumerate() {
            if *v > best {
                best = *v;
                at = i;
            }
        }
        if best <= 0.0 || at == 0 || at + 1 >= n {
            return None;
        }
        let (a, b, c) = (r[at - 1], r[at], r[at + 1]);
        let denom = a - 2.0 * b + c;
        let shift = if denom.abs() > 0.0 { 0.5 * (a - c) / denom } else { 0.0 };
        let dtheta = (self.range[1] - self.range[0]) / n as f64;
        Some(self.range[0] + dtheta * (at as f64 + 0.5 + shift))
    }

    pub fn total_flux(&self) -> f64 {
        self.flux.iter().sum()
    }
}

/// How photons are traced.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct PhotonOptions {
    /// Longest specular chain followed.
    pub max_bounces: usize,
    /// Wavelength range sampled, metres. A single-wavelength range makes
    /// the run monochromatic, which is what the closed-form references
    /// are stated for.
    pub spectrum: (f64, f64),
    /// Stop following a photon whose flux has fallen below this fraction
    /// of its original.
    pub min_weight: f64,
    /// Whether to deposit photons that reach the receiver without ever
    /// meeting the specular geometry. Off by default: the caustic is the
    /// interesting part, and the direct beam swamps it.
    pub include_direct: bool,
    /// Keep only photons with exactly this many specular events.
    ///
    /// This is how one caustic is separated from another, and it is a
    /// selection rather than an approximation. A sphere in parallel light
    /// makes a different caustic for each path: two events (straight
    /// through) is the bright forward beam, three (in, one internal
    /// reflection, out) is the primary rainbow, four is the secondary.
    /// Without this the forward beam buries the bows, which is exactly
    /// what it does in the sky.
    pub select_events: Option<usize>,
}

impl Default for PhotonOptions {
    fn default() -> Self {
        Self {
            max_bounces: 12,
            spectrum: (550e-9, 550e-9),
            min_weight: 1e-4,
            include_direct: false,
            select_events: None,
        }
    }
}

/// What a photon pass did.
#[derive(Clone, Copy, Debug, Default)]
pub struct PhotonStats {
    pub emitted: u64,
    /// Photons that reached the receiver after at least one specular
    /// event.
    pub deposited: u64,
    /// Photons that never met the geometry.
    pub missed: u64,
    /// Photons absorbed, or lost to the weight cutoff or the bounce
    /// limit.
    pub lost: u64,
    pub flux_emitted: f64,
    pub flux_deposited: f64,
    /// Flux that reached the receiver plane outside its binned extent.
    pub flux_spilled: f64,
    /// Flux lost to Fresnel absorption, cutoffs and the bounce limit.
    pub flux_lost: f64,
    pub specular_events: u64,
}

impl PhotonStats {
    /// Flux accounted for, as a fraction of what was emitted: the
    /// bookkeeping check that no energy appeared or vanished.
    pub fn accounted(&self) -> f64 {
        if self.flux_emitted <= 0.0 {
            return 1.0;
        }
        (self.flux_deposited + self.flux_spilled + self.flux_lost) / self.flux_emitted
    }
}

/// Index of refraction of a material at a wavelength, if it has one.
fn index_of(material: &Material, lambda: f64) -> Option<f64> {
    material.refractive_index(lambda)
}

/// Trace `count` photons from `emitter` through `scene` onto `receiver`.
///
/// Photons follow specular events only: a dielectric refracts or reflects
/// by Fresnel's coefficients, a conductor reflects. A photon that meets a
/// diffuse surface is absorbed there, because a caustic is by definition
/// what arrives through the specular chain. Flux is accounted for at every
/// step, and [`PhotonStats::accounted`] must come out at one.
pub fn trace_to_plane(
    scene: &Scene,
    emitter: &Emitter,
    receiver: &mut Receiver,
    count: usize,
    seed: u64,
    options: &PhotonOptions,
) -> PhotonStats {
    let march = MarchOptions::for_scale(scene.scale());
    let mut stats = PhotonStats::default();
    for i in 0..count {
        let mut rng = Rng(mix(seed, i as u64));
        let (ray, flux) = emitter.emit(&mut rng, count);
        let lambda = sample_lambda(&mut rng, options.spectrum);
        let xyz = radiation::optics::cie_xyz(lambda);
        stats.emitted += 1;
        stats.flux_emitted += flux;
        match walk(scene, ray, flux, lambda, &march, options, &mut rng, &mut stats) {
            Walk::Escaped { ray, flux, bounces } => {
                if !selected(bounces, options) {
                    stats.missed += 1;
                    stats.flux_lost += flux;
                    continue;
                }
                match receiver.hit(ray) {
                    Some((_, point)) => {
                        let before = receiver.spilled;
                        if receiver.deposit(point, flux, xyz) {
                            stats.deposited += 1;
                        }
                        stats.flux_deposited += flux - (receiver.spilled - before);
                        stats.flux_spilled += receiver.spilled - before;
                        receiver.photons += 1;
                    }
                    None => {
                        stats.missed += 1;
                        stats.flux_lost += flux;
                    }
                }
            }
            Walk::Absorbed { flux } => {
                stats.lost += 1;
                stats.flux_lost += flux;
            }
        }
    }
    stats
}

/// As [`trace_to_plane`], binning the exit *direction* instead: the
/// measurement a rainbow angle is read from.
pub fn trace_to_angles(
    scene: &Scene,
    emitter: &Emitter,
    receiver: &mut Angular,
    count: usize,
    seed: u64,
    options: &PhotonOptions,
) -> PhotonStats {
    let march = MarchOptions::for_scale(scene.scale());
    let mut stats = PhotonStats::default();
    for i in 0..count {
        let mut rng = Rng(mix(seed, i as u64));
        let (ray, flux) = emitter.emit(&mut rng, count);
        let lambda = sample_lambda(&mut rng, options.spectrum);
        stats.emitted += 1;
        stats.flux_emitted += flux;
        match walk(scene, ray, flux, lambda, &march, options, &mut rng, &mut stats) {
            Walk::Escaped { ray, flux, bounces } => {
                if !selected(bounces, options) {
                    stats.missed += 1;
                    stats.flux_lost += flux;
                    continue;
                }
                let before = receiver.outside;
                receiver.deposit(ray.dir.0, flux);
                receiver.photons += 1;
                let spilled = receiver.outside - before;
                if spilled <= 0.0 {
                    stats.deposited += 1;
                }
                stats.flux_deposited += flux - spilled;
                stats.flux_spilled += spilled;
            }
            Walk::Absorbed { flux } => {
                stats.lost += 1;
                stats.flux_lost += flux;
            }
        }
    }
    stats
}

/// Whether a photon with this many specular events is kept.
fn selected(bounces: usize, options: &PhotonOptions) -> bool {
    match options.select_events {
        Some(k) => bounces == k,
        None => bounces > 0 || options.include_direct,
    }
}

enum Walk {
    Escaped { ray: Ray, flux: f64, bounces: usize },
    Absorbed { flux: f64 },
}

/// Follow one photon through the specular chain until it leaves the scene
/// or is absorbed.
#[allow(clippy::too_many_arguments)]
fn walk(
    scene: &Scene,
    mut ray: Ray,
    mut flux: f64,
    lambda: f64,
    march: &MarchOptions,
    options: &PhotonOptions,
    rng: &mut Rng,
    stats: &mut PhotonStats,
) -> Walk {
    let start = flux;
    let mut bounces = 0usize;
    for _ in 0..options.max_bounces {
        let (hit, _, _) = scene.intersect(ray, march);
        let Some(h) = hit else {
            return Walk::Escaped { ray, flux, bounces };
        };
        let Some(material) = scene.material_of(&h) else {
            return Walk::Absorbed { flux };
        };
        let front = dot(h.hit.normal, ray.dir.0) < 0.0;
        let n = if front { h.hit.normal } else { scale(h.hit.normal, -1.0) };
        let offset = 4.0 * march.tolerance;

        match material {
            Material::Dielectric { .. } | Material::Dispersive { .. } => {
                let index = index_of(material, lambda).unwrap_or(1.5);
                let (ni, nt) = if front { (1.0, index) } else { (index, 1.0) };
                let cos = dot(scale(ray.dir.0, -1.0), n).abs();
                let r = fresnel_dielectric(cos, ni, nt);
                // Russian roulette between the two lobes: unbiased, and
                // it keeps every photon at full weight so the caustic is
                // not noisier than it has to be.
                let reflect = rng.uniform() < r;
                let dir = if reflect {
                    ray.dir.reflect(Vec3(n))
                } else {
                    match refract(ray.dir, Vec3(n), ni, nt) {
                        Some(d) => d,
                        // Total internal reflection: not a special case to
                        // be approximated, just what Snell gives.
                        None => ray.dir.reflect(Vec3(n)),
                    }
                };
                let side = if dot(dir.0, n) > 0.0 { 1.0 } else { -1.0 };
                ray = Ray::new(Vec3(add(h.hit.point, scale(n, side * offset))), dir);
            }
            Material::Conductor { index, extinction } => {
                let cos = dot(scale(ray.dir.0, -1.0), n).abs();
                flux *= fresnel_conductor(cos, *index, *extinction);
                let dir = ray.dir.reflect(Vec3(n));
                ray = Ray::new(Vec3(add(h.hit.point, scale(n, offset))), dir);
            }
            // A mirror that reflects everything: the catacaustic reference
            // wants no Fresnel falloff confusing the geometry.
            Material::Mirror => {
                let dir = ray.dir.reflect(Vec3(n));
                ray = Ray::new(Vec3(add(h.hit.point, scale(n, offset))), dir);
            }
            _ => return Walk::Absorbed { flux },
        }
        bounces += 1;
        stats.specular_events += 1;
        if flux < options.min_weight * start {
            return Walk::Absorbed { flux };
        }
    }
    Walk::Absorbed { flux }
}

fn sample_lambda(rng: &mut Rng, spectrum: (f64, f64)) -> f64 {
    if spectrum.1 <= spectrum.0 {
        spectrum.0
    } else {
        spectrum.0 + (spectrum.1 - spectrum.0) * rng.uniform()
    }
}

fn mix(seed: u64, i: u64) -> u64 {
    let mut x = seed ^ i.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    (x ^ (x >> 31)) | 1
}

// ---------------------------------------------------------------------
// Closed-form references
// ---------------------------------------------------------------------

/// Deviation of a ray through a sphere with `k` internal reflections, as a
/// function of the incidence angle (Descartes' construction).
///
/// ```text
///   D(θ₁) = 2(θ₁ − θ₂) + k(π − 2θ₂)
/// ```
///
/// with `sin θ₁ = n sin θ₂`. `None` when the ray cannot enter.
pub fn sphere_deviation(incidence: f64, index: f64, internal: usize) -> Option<f64> {
    let s = incidence.sin() / index;
    if s.abs() > 1.0 {
        return None;
    }
    let refracted = s.asin();
    Some(2.0 * (incidence - refracted) + internal as f64 * (std::f64::consts::PI - 2.0 * refracted))
}

/// The incidence angle of least deviation, exactly:
/// `cos θ₁ = √((n² − 1) / (k² + 2k))`.
///
/// For one internal reflection this is `√((n²−1)/3)`, the primary
/// rainbow; for two, `√((n²−1)/8)`, the secondary (Descartes 1637; Born &
/// Wolf §14.1).
pub fn rainbow_incidence(index: f64, internal: usize) -> Option<f64> {
    let k = internal as f64;
    let denom = k * k + 2.0 * k;
    if denom <= 0.0 {
        return None;
    }
    let c2 = (index * index - 1.0) / denom;
    if !(0.0..=1.0).contains(&c2) {
        return None;
    }
    Some(c2.sqrt().acos())
}

/// Deviation at the rainbow: the stationary value of
/// [`sphere_deviation`].
pub fn rainbow_deviation(index: f64, internal: usize) -> Option<f64> {
    sphere_deviation(rainbow_incidence(index, internal)?, index, internal)
}

/// The angle a rainbow is *seen* at, measured from the antisolar point:
/// `π − D` for the primary, `D − π` for the secondary.
pub fn rainbow_angle(index: f64, internal: usize) -> Option<f64> {
    let d = rainbow_deviation(index, internal)?;
    Some((std::f64::consts::PI - d).abs())
}

/// Where a ray parallel to the axis, at height `h`, crosses the axis after
/// reflecting in a spherical mirror of radius `R`.
///
/// `R/2` in the limit `h → 0`, and shorter for every finite `h`: the
/// caustic of a spherical mirror is the locus of these crossings, whose
/// cusp is the paraxial focus. Exact.
pub fn mirror_axis_crossing(radius: f64, h: f64) -> Option<f64> {
    if h.abs() >= radius.abs() {
        return None;
    }
    // The surface z = R - sqrt(R^2 - h^2) measured from the vertex, with
    // the normal along the radius.
    let cos = (1.0 - (h / radius).powi(2)).sqrt();
    // Reflected ray makes angle 2*asin(h/R) with the axis.
    let theta = (h / radius).asin();
    let z = radius * (1.0 - cos);
    Some(z + h / (2.0 * theta).tan())
}

/// The catacaustic of a spherical mirror in parallel light: the envelope
/// of the reflected rays, at the parameter `t` where the ray strikes at
/// height `R sin t`.
///
/// Derived here rather than quoted, because the derivation is short and a
/// misremembered nephroid is worse than none. With the vertex at the
/// origin and the concave side facing the light, the ray family is
///
/// ```text
///   P(t) = (R sin t, R(1 − cos t)),   r(t) = (−sin 2t, cos 2t)
/// ```
///
/// and the envelope condition `(P' + u r') ∥ r` gives `u = R cos t / 2`,
/// so
///
/// ```text
///   x_c = R sin³t,   z_c = R(1 − (3/2) cos t + cos³t).
/// ```
///
/// At `t = 0` that is `(0, R/2)`: the cusp of the caustic is the paraxial
/// focus. The tests check the formula against the intersection of
/// neighbouring rays, which is what an envelope is.
pub fn mirror_catacaustic(radius: f64, t: f64) -> [f64; 2] {
    let (s, c) = t.sin_cos();
    [
        radius * s * s * s,
        radius * (1.0 - 1.5 * c + c * c * c),
    ]
}

/// The distance along the reflected ray at which it touches the caustic,
/// `R cos t / 2`.
pub fn mirror_caustic_parameter(radius: f64, t: f64) -> f64 {
    0.5 * radius * t.cos()
}

/// The reflected ray at parameter `t`: its origin on the mirror and its
/// unit direction, for the concave orientation used by
/// [`mirror_catacaustic`].
pub fn mirror_reflected_ray(radius: f64, t: f64) -> ([f64; 2], [f64; 2]) {
    let (s, c) = t.sin_cos();
    (
        [radius * s, radius * (1.0 - c)],
        [-(2.0 * t).sin(), (2.0 * t).cos()],
    )
}

// ---------------------------------------------------------------------
// The wave caustic
// ---------------------------------------------------------------------

/// The Airy function `Ai(x)`, by its Taylor series about zero with the
/// oscillatory asymptotic beyond.
///
/// `Ai'' = x Ai`, so the series coefficients follow from
/// `a_{n+2} = a_{n−1} / ((n+2)(n+1))` with
/// `Ai(0) = 3^{−2/3}/Γ(2/3) = 0.3550280539` and
/// `Ai'(0) = −3^{−1/3}/Γ(1/3) = −0.2588194038`. The tests check the
/// published values and the first three zeros.
///
/// The series is entire, but in double precision it cancels: the terms
/// grow like `x^m` while `Ai(x)` decays, so the asymptotic expansions take
/// over where they are the more accurate of the two. Those crossovers were
/// placed by comparing both branches against the same series evaluated in
/// 60-digit arithmetic, not guessed. The worst relative error anywhere is
/// about `1e-7`, just above `x = 6`, where `Ai` is already `1e-6`.
pub fn airy(x: f64) -> f64 {
    /// `Ai(0) = 3^{-2/3} / Gamma(2/3)`.
    const AI0: f64 = 0.355_028_053_887_817_24;
    /// `Ai'(0) = -3^{-1/3} / Gamma(1/3)`.
    const AIP0: f64 = -0.258_819_403_792_806_80;
    // Asymmetric on purpose. For x > 0 the series cancels badly, because
    // Ai decays exponentially; for x < 0 it holds much further, because
    // the answer stays O(1).
    if (-7.0..=6.0).contains(&x) {
        // a_0 = Ai(0), a_1 = Ai'(0), a_2 = 0, a_m = a_{m-3} / (m (m-1)).
        let mut coeff = [AI0, AIP0, 0.0];
        let mut power = [1.0, x, x * x];
        let mut sum = AI0 + AIP0 * x;
        let x3 = x * x * x;
        // `a_2 = 0`, so every third coefficient vanishes exactly. A break
        // that looks at one term therefore stops on a zero and truncates
        // the series: the last three terms have to be small together.
        let mut recent = [f64::INFINITY; 3];
        for m in 3..120usize {
            let a = coeff[0] / ((m * (m - 1)) as f64);
            let p = power[0] * x3;
            sum += a * p;
            coeff = [coeff[1], coeff[2], a];
            power = [power[1], power[2], p];
            recent = [recent[1], recent[2], (a * p).abs()];
            let worst = recent[0].max(recent[1]).max(recent[2]);
            if m > 12 && worst < 1e-19 * sum.abs().max(1e-300) {
                break;
            }
        }
        sum
    } else if x < 0.0 {
        asymptotic_negative(-x)
    } else {
        asymptotic_positive(x)
    }
}

/// The coefficients `u_k` of the Airy asymptotic expansions (DLMF 9.7.1),
/// from `u_k = (6k-5)(6k-3)(6k-1) / (216 k (2k-1)) * u_{k-1}`.
const AIRY_U: [f64; 11] = [
    1.0, // u_0
    6.94444444444444475e-02, // u_1 = 5/72,
    3.71334876543209860e-02, // u_2 = 385/10368,
    3.79930591278006408e-02, // u_3 = 85085/2239488,
    5.76491904126697208e-02, // u_4 = 37182145/644972544,
    1.16099064025515411e-01, // u_5 = 5391411025/46438023168,
    2.91591399230750503e-01, // u_6 = 5849680962125/20061226008576,
    8.77666969510016881e-01, // u_7 = 1267709431363375/1444408272617472,
    3.07945303017316707e+00, // u_8 = 2562040760785380875/831979165027663872,
    1.23415733323452379e+01, // u_9 = 6653619855759634132375/539122498937926189056,
    5.56227853659170819e+01, // u_10 = 4318199286388002551911375/77633639847061371224064
];

/// `Ai(-z)` for large `z > 0`, by the oscillatory asymptotic expansion
/// (DLMF 9.7.9):
///
/// ```text
///   Ai(-z) ~ z^{-1/4} / sqrt(pi) * [ cos(zeta - pi/4) * sum (-1)^k u_2k / zeta^2k
///                                  + sin(zeta - pi/4) * sum (-1)^k u_2k+1 / zeta^2k+1 ]
/// ```
///
/// with `zeta = (2/3) z^{3/2}`. The `sin` series is not optional: near a
/// zero of the cosine the leading term alone is relatively wrong by
/// `u_1 / (zeta tan(zeta - pi/4))`, which is 3% at `z = 8` — the error
/// that `the_airy_function_matches_its_published_values` caught by
/// comparing the two branches where they meet.
pub fn asymptotic_negative(z: f64) -> f64 {
    let zeta = 2.0 / 3.0 * z.powf(1.5);
    let (mut even, mut odd) = (0.0, 0.0);
    for k in 0..5 {
        let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
        even += sign * AIRY_U[2 * k] / zeta.powi(2 * k as i32);
        odd += sign * AIRY_U[2 * k + 1] / zeta.powi(2 * k as i32 + 1);
    }
    let phase = zeta - std::f64::consts::FRAC_PI_4;
    z.powf(-0.25) / std::f64::consts::PI.sqrt() * (phase.cos() * even + phase.sin() * odd)
}

/// `Ai(z)` for large `z > 0`, by the exponentially small asymptotic
/// expansion (DLMF 9.7.5).
///
/// Public so that the handover in [`airy`] can be tested for what it is —
/// two branches evaluated at the *same* point — rather than by nudging the
/// argument across the crossover, which measures `Ai'` instead.
pub fn asymptotic_positive(z: f64) -> f64 {
    let zeta = 2.0 / 3.0 * z.powf(1.5);
    let mut sum = 0.0;
    for k in 0..7 {
        let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
        sum += sign * AIRY_U[k] / zeta.powi(k as i32);
    }
    (-zeta).exp() / (2.0 * std::f64::consts::PI.sqrt() * z.powf(0.25)) * sum
}

/// Irradiance across a **fold** caustic — the generic caustic, the one a
/// lens or a water surface makes — normalised to one far on the lit side.
///
/// Geometrical optics gives `1/√s` on the lit side and nothing on the
/// dark: infinite at the caustic and discontinuous across it. The exact
/// wave answer is `Ai²`, which is finite everywhere and decays into the
/// shadow (Airy 1838; Berry & Upstill, *Progress in Optics* 18, 1980):
///
/// ```text
///   I(s) = (2π / √w) Ai²(−s / w)
/// ```
///
/// where `s` is the distance from the caustic along its normal and `w` is
/// the caustic's wave scale, [`fold_scale`]. The constant is fixed by
/// requiring the mean of the oscillations to match `1/√s` far on the lit
/// side, which the tests check: `Ai²(−z)` averages to `1/(2π√z)`, so
/// `(2π/√w) · 1/(2π√(s/w)) = 1/√s`.
///
/// Two consequences are what make a real caustic bright but finite. The
/// peak irradiance is `≈ 1.80 / √w`, so it grows only as `λ^{-1/3}`; and
/// the fringes beside it are spaced as `w`, so as `λ^{2/3}`.
pub fn fold(distance: f64, scale: f64) -> f64 {
    if scale <= 0.0 {
        return fold_geometric(distance);
    }
    let a = airy(-distance / scale);
    2.0 * std::f64::consts::PI / scale.sqrt() * a * a
}

/// The peak irradiance of a fold caustic: `2π Ai²(−1.0188) / √w`, about
/// `1.8028 / √w`.
pub fn fold_peak(scale: f64) -> f64 {
    fold(1.018_792_971_647_471 * scale, scale)
}

/// The wave scale of a fold caustic: `(λ² ρ / 2)^{1/3}` for a caustic of
/// radius of curvature `ρ` at wavelength `λ`.
///
/// This is the width over which the Airy pattern replaces the geometric
/// singularity, and the renderer's indicator for caustics: a caustic
/// thinner than a pixel is a ray-optics line, and one wider is a fringed
/// band that ray optics cannot draw at all.
pub fn fold_scale(wavelength: f64, curvature_radius: f64) -> f64 {
    (wavelength * wavelength * curvature_radius.abs() / 2.0).cbrt()
}

/// Geometrical optics' answer for the same fold: `1/√s` on the lit side,
/// zero on the dark, infinite at the caustic.
pub fn fold_geometric(distance: f64) -> f64 {
    if distance > 0.0 {
        distance.sqrt().recip()
    } else if distance == 0.0 {
        // Not a convenient value to return but the honest one: this is
        // where the model fails, and returning a finite number would hide
        // that.
        f64::INFINITY
    } else {
        0.0
    }
}
