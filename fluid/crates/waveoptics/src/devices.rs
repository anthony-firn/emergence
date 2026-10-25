//! Rendering devices: cameras with real lenses, pinhole cameras and the
//! human eye, looking at an interference experiment.
//!
//! A device is placed with its front (first lens vertex, cornea, or
//! pinhole plane) at transverse position `(x, y)` and distance `z` behind
//! the openings, facing along +z. For every incoherent source mode the
//! band-limited angular spectrum of the openings ([`crate::angular`]) gives
//! the complex field on the device's entrance pupil (in object space),
//! restricted to the directions that can reach the detector window; the
//! device's pupil function at that
//! wavelength — the traced aberrations of the real optics, or the Fresnel
//! phase of a pinhole — turns it into the detector field
//! ([`crate::imaging::image_field`]). Mode intensities add. The result is
//! spectral irradiance on the detector, which a [`crate::detector::Sensor`]
//! or a [`Retina`] then turns into counts.
//!
//! Image coordinates are scene-oriented (as displayed after the device's
//! usual inversion).

use std::f64::consts::PI;

use num_complex::Complex64 as Cx;

use crate::angular::{plane_fields, Band, OpeningSpectrum};
use crate::detector::{Rng, Sensor, HC};
use crate::diffraction::Incident;
use crate::polarization::Jones;
use crate::experiment::Experiment;
use crate::imaging::{image_field, PupilFunction, Window};
use crate::raytrace::System;

/// Largest pupil grid per axis (memory: n² complex samples per opening).
pub const MAX_PUPIL: usize = 1024;

/// The optics of a device.
#[derive(Clone, Debug)]
pub enum Optics {
    /// A real refracting system (camera lens or eye) with its detector at
    /// `System::image_z`.
    Refractive(System),
    /// A lensless pinhole of `diameter` with the detector `distance` behind.
    Pinhole { diameter: f64, distance: f64 },
    /// An aberration-free thin lens (pupil radius, focal length) with the
    /// detector in its focal plane: the reference optic.
    Ideal { radius: f64, focal: f64 },
}

/// A device placed in an experiment.
#[derive(Clone, Debug)]
pub struct Device {
    pub optics: Optics,
    /// Front of the device: `(x, y, z)` with `z` measured from the openings.
    pub front: [f64; 3],
}

/// Spectral irradiance on a detector window: `(λ, E_λ)` per wavelength
/// (each already weighted by its share of the source power), plus how the
/// pupil was sampled.
#[derive(Clone, Debug)]
pub struct SpectralImage {
    pub window: Window,
    pub bands: Vec<(f64, Vec<f64>)>,
    pub pupil_samples: usize,
    pub modes: usize,
    /// Largest fraction of an opening's power outside the recorded band
    /// (light that cannot reach the window, up to pupil diffraction).
    pub band_loss: f64,
    /// Largest bound on the incident-wave model's phase error (rad).
    pub phase_bound: f64,
}

impl SpectralImage {
    /// Total irradiance (W/m²) per window sample.
    pub fn irradiance(&self) -> Vec<f64> {
        let mut out = vec![0.0; self.window.pixels[0] * self.window.pixels[1]];
        for (_, e) in &self.bands {
            for (o, v) in out.iter_mut().zip(e) {
                *o += v;
            }
        }
        out
    }

    /// CIE 2006 XYZ (W/m² weighted by the colour matching functions) per
    /// sample: the colour an ideal colorimetric camera records.
    pub fn xyz(&self) -> Vec<[f64; 3]> {
        let mut out = vec![[0.0; 3]; self.window.pixels[0] * self.window.pixels[1]];
        for (l, e) in &self.bands {
            let c = crate::color::xyz(*l);
            for (o, v) in out.iter_mut().zip(e) {
                for k in 0..3 {
                    o[k] += c[k] * v;
                }
            }
        }
        out
    }
}

impl SpectralImage {
    /// Expected photo-electrons per pixel of `sensor`, the window being the
    /// sensor's pixel grid subdivided `sub × sub` ([`sensor_window`]):
    /// midpoint rule over each pixel, silicon QE per band, fill factor as a
    /// uniform area fraction.
    pub fn electrons(&self, sensor: &Sensor, sub: usize) -> Vec<f64> {
        let [px, py] = sensor.pixels;
        assert_eq!(self.window.pixels, [px * sub, py * sub]);
        let area = self.window.pitch * self.window.pitch * sensor.fill_factor;
        let mut out = vec![0.0; px * py];
        for (l, e) in &self.bands {
            let k = sensor.exposure * l / HC * sensor.qe.at(*l) * area;
            for (q, v) in e.iter().enumerate() {
                let (i, j) = ((q % (px * sub)) / sub, (q / (px * sub)) / sub);
                out[i + px * j] += k * v;
            }
        }
        out
    }
}

/// The detector window of a sensor's pixels subdivided `sub × sub`.
pub fn sensor_window(sensor: &Sensor, sub: usize) -> Window {
    Window { center: sensor.center, pixels: [sensor.pixels[0] * sub, sensor.pixels[1] * sub], pitch: sensor.pitch / sub as f64 }
}

/// Subdivisions per pixel so that the midpoint rule over a pixel of
/// `pitch` resolves a pattern of period `period` with relative error
/// `(π pitch/(sub · period))²/24 < tol`.
pub fn subdivisions(pitch: f64, period: f64, tol: f64) -> usize {
    ((PI * pitch / period) / (24.0 * tol).sqrt()).ceil().max(1.0) as usize
}

impl Device {
    fn pupil(&self, wavelength: f64, n: usize) -> PupilFunction {
        match &self.optics {
            Optics::Refractive(s) => PupilFunction::traced(s, wavelength, n),
            Optics::Pinhole { diameter, distance } => PupilFunction::pinhole(*diameter, *distance, wavelength, n),
            Optics::Ideal { radius, focal } => PupilFunction::ideal(*radius, *focal, wavelength, n, 0.0),
        }
    }

    /// Focal length (detector distance for a pinhole) at a wavelength.
    pub fn focal(&self, wavelength: f64) -> f64 {
        match &self.optics {
            Optics::Refractive(s) => s.paraxial(wavelength).efl,
            Optics::Pinhole { distance, .. } => *distance,
            Optics::Ideal { focal, .. } => *focal,
        }
    }

    /// Entrance-pupil radius and plane (from the front) at a wavelength.
    fn pupil_geometry(&self, wavelength: f64) -> (f64, f64) {
        match &self.optics {
            Optics::Refractive(s) => {
                let p = s.paraxial(wavelength);
                (p.ep_radius, p.ep_z)
            }
            Optics::Pinhole { diameter, .. } => (0.5 * diameter, 0.0),
            Optics::Ideal { radius, .. } => (*radius, 0.0),
        }
    }

    /// Largest phase step between neighbouring open pupil samples (rad).
    fn pupil_phase_step(p: &PupilFunction) -> f64 {
        let n = p.n;
        let mut worst: f64 = 0.0;
        for j in 0..n {
            for i in 0..n {
                let a = p.values[i + n * j];
                if a.norm() == 0.0 {
                    continue;
                }
                for (ii, jj) in [(i + 1, j), (i, j + 1)] {
                    if ii < n && jj < n {
                        let b = p.values[ii + n * jj];
                        if b.norm() > 0.0 {
                            worst = worst.max((b * a.conj()).arg().abs());
                        }
                    }
                }
            }
        }
        worst
    }

    /// The device's pupil at one wavelength with the largest ray-angle
    /// error of its own wavefront, `max |∇W|` (rad), resolved to < π/2 per
    /// sample.
    fn resolved_pupil(&self, wavelength: f64, n_min: usize) -> (PupilFunction, f64) {
        let mut n = n_min.clamp(64, MAX_PUPIL);
        loop {
            let p = self.pupil(wavelength, n);
            let step = Self::pupil_phase_step(&p);
            if step <= 0.5 * PI || n >= MAX_PUPIL {
                let slope = step / (2.0 * PI / wavelength * p.spacing());
                return (p, slope);
            }
            n = (2 * n).min(MAX_PUPIL);
        }
    }

    /// Spectral irradiance on `window` (detector coordinates) of `exp`.
    ///
    /// Per wavelength, the band of spatial frequencies that can reach the
    /// window is: the window's direction sines, widened by the device's own
    /// ray-angle errors (`max|∇W|`) and by six diffraction widths
    /// `λ/(2r)` of the pupil, with a raised-cosine roll-off over six more.
    /// The pupil step then keeps every replica of the pupil-weighted field
    /// out of the window (`1/h ≥ f_outer + max|∇W|/λ + w/λ`), and the
    /// frequency grid keeps the replicas of the band-limited field clear of
    /// the pupil ([`crate::angular`]).
    pub fn spectral_image(&self, exp: &Experiment, window: &Window) -> SpectralImage {
        // Spectral and source sampling: two openings reach a window point at
        // direction sine θ with path difference ≤ b(|θ| + max|∇W|) (b the
        // baseline; the device's own wavefront adds at most its slope times
        // the opening separation in the pupil), plus the source side.
        let shortest = exp.source.spectrum.span().0;
        let f0 = self.focal(exp.source.spectrum.center());
        let (_, slope0) = self.resolved_pupil(shortest, 64);
        let w0 = (window.center[0].abs() + 0.5 * window.pitch * window.pixels[0] as f64).hypot(window.center[1].abs() + 0.5 * window.pitch * window.pixels[1] as f64) / f0;
        let max_path = exp.baseline() * (w0 + slope0) + exp.source_path();
        let modes = exp.source.modes(exp.baseline(), max_path, exp.tol);
        let mut lambdas: Vec<f64> = Vec::new();
        for m in &modes {
            if !lambdas.contains(&m.wavelength) {
                lambdas.push(m.wavelength);
            }
        }
        let npx = window.pixels[0] * window.pixels[1];
        let scalar = exp.polarization_free();
        let half = |a: usize| 0.5 * window.pitch * window.pixels[a] as f64;
        let mut bands = Vec::with_capacity(lambdas.len());
        let (mut n_max, mut band_loss, mut phase_bound): (usize, f64, f64) = (0, 0.0, 0.0);
        for &lambda in &lambdas {
            let f = self.focal(lambda);
            let (r, _) = self.pupil_geometry(lambda);
            let (_, slope) = self.resolved_pupil(lambda, 64);
            let diff = 0.5 / r; // pupil diffraction width (cycles/m)
            let band = Band {
                center: [window.center[0] / (lambda * f), window.center[1] / (lambda * f)],
                inner: std::array::from_fn(|a| (half(a) / f + slope) / lambda + 6.0 * diff),
                outer: std::array::from_fn(|a| (half(a) / f + slope) / lambda + 12.0 * diff),
            };
            let need = (0..2).map(|a| band.outer[a] + slope / lambda + half(a) / (lambda * f)).fold(0.0, f64::max);
            let (pupil, _) = self.resolved_pupil(lambda, (((2.0 * r * need).ceil() as usize) + 1) & !1);
            n_max = n_max.max(pupil.n);
            let z = self.front[2] + pupil.plane;
            let h = pupil.spacing();
            let coords = |a: usize| -> Vec<f64> { (0..pupil.n).map(|i| self.front[a] + pupil.center[a] - pupil.radius + h * (i as f64 + 0.5)).collect() };
            let (xs, ys) = (coords(0), coords(1));
            // replica spacing: span of the band-limited field (openings
            // displaced by λ z f over the band, plus the roll-off's ringing)
            // together with the pupil
            let period: [f64; 2] = std::array::from_fn(|a| {
                let ring = 4.0 / (band.outer[a] - band.inner[a]);
                let lo = exp.openings.iter().map(|o| o.rect.center[a] - 0.5 * o.rect.size[a]).fold(f64::INFINITY, f64::min) + lambda * z * (band.center[a] - band.outer[a]) - ring;
                let hi = exp.openings.iter().map(|o| o.rect.center[a] + 0.5 * o.rect.size[a]).fold(f64::NEG_INFINITY, f64::max) + lambda * z * (band.center[a] + band.outer[a]) + ring;
                let (plo, phi) = (self.front[a] - pupil.radius, self.front[a] + pupil.radius);
                (hi - plo).max(phi - lo)
            });
            let mut acc = vec![0.0; npx];
            // group the modes of this wavelength by illumination
            let mut groups: Vec<(Incident, Vec<(Jones, f64)>)> = Vec::new();
            for m in modes.iter().filter(|m| m.wavelength == lambda) {
                match groups.iter_mut().find(|g| g.0 == m.incident) {
                    Some(g) => g.1.push((m.jones, m.weight)),
                    None => groups.push((m.incident, vec![(m.jones, m.weight)])),
                }
            }
            for (incident, pols) in &groups {
                let spectra: Vec<OpeningSpectrum> = exp.openings.iter().map(|o| OpeningSpectrum::new(&o.rect, incident, lambda)).collect();
                phase_bound = spectra.iter().map(|s| s.phase_bound).fold(phase_bound, f64::max);
                let (fields, loss) = plane_fields(&spectra, lambda, z, &band, period, &xs, &ys);
                band_loss = loss.into_iter().fold(band_loss, f64::max);
                if scalar {
                    let total_w: f64 = pols.iter().map(|p| p.1).sum();
                    let u: Vec<Cx> = (0..xs.len() * ys.len()).map(|q| fields.iter().map(|f| f[q]).sum()).collect();
                    let img = image_field(&pupil, &u, window);
                    for (a, v) in acc.iter_mut().zip(&img) {
                        *a += total_w * v.norm_sqr();
                    }
                } else {
                    for (jones, w) in pols {
                        for comp in 0..2 {
                            let u: Vec<Cx> = (0..xs.len() * ys.len())
                                .map(|q| {
                                    let mut e = Jones([Cx::new(0.0, 0.0); 2]);
                                    for (o, f) in exp.openings.iter().zip(&fields) {
                                        e = e.add(&o.element.apply(jones).scale(f[q]));
                                    }
                                    if let Some(a) = &exp.analyzer {
                                        e = a.apply(&e);
                                    }
                                    e.0[comp]
                                })
                                .collect();
                            let img = image_field(&pupil, &u, window);
                            for (a, v) in acc.iter_mut().zip(&img) {
                                *a += w * v.norm_sqr();
                            }
                        }
                    }
                }
            }
            bands.push((lambda, acc));
        }
        SpectralImage { window: *window, bands, pupil_samples: n_max, modes: modes.len(), band_loss, phase_bound }
    }
}

/// The cone mosaic of the central retina.
#[derive(Clone, Debug)]
pub struct Retina {
    /// Centre-to-centre cone spacing (m); ~2.5 µm at the foveal centre
    /// (Curcio et al., J. Comp. Neurol. 292, 497 (1990): peak density
    /// ~1.6–2 × 10⁵ cones/mm²).
    pub spacing: f64,
    /// L : M : S proportions (estimate: 0.62 : 0.31 : 0.07; S cones absent
    /// within `s_free_radius`).
    pub proportions: [f64; 3],
    /// Radius of the S-cone-free foveola (m), ~0.1 mm (≈ 0.35°).
    pub s_free_radius: f64,
    /// Effective collecting area at peak sensitivity (m²): aperture area ×
    /// peak absorptance × isomerisation efficiency (estimate:
    /// 3.1 µm² × 0.68 × 0.67).
    pub peak_area: f64,
    /// Temporal integration time (s), ~0.1 s (estimate).
    pub integration: f64,
}

impl Default for Retina {
    fn default() -> Self {
        Self {
            spacing: 2.5e-6,
            proportions: [0.62, 0.31, 0.07],
            s_free_radius: 0.1e-3,
            peak_area: 3.1e-12 * 0.68 * 0.67,
            integration: 0.1,
        }
    }
}

/// One cone and its response.
#[derive(Clone, Copy, Debug)]
pub struct Cone {
    pub position: [f64; 2],
    /// 0 = L, 1 = M, 2 = S.
    pub kind: usize,
    /// Expected and sampled photoisomerisations in the integration time.
    pub expected: f64,
    pub isomerisations: f64,
    /// Energy-weighted excitation (CIE 2006 fundamental units, W/m²
    /// weighted), scaled by the sampled/expected ratio (photon noise).
    pub excitation: f64,
}

impl Retina {
    /// Sample the retinal image with a triangular cone lattice over the
    /// window (random type assignment by proportion; no S cones within
    /// `s_free_radius` of the fovea at the origin), photon noise included.
    pub fn sample(&self, image: &SpectralImage, rng: &mut Rng) -> Vec<Cone> {
        let w = &image.window;
        let (xs, ys) = (w.coords(0), w.coords(1));
        let (x0, x1) = (xs[0], xs[xs.len() - 1]);
        let (y0, y1) = (ys[0], ys[ys.len() - 1]);
        let row = self.spacing * 3f64.sqrt() / 2.0;
        let mut cones = Vec::new();
        let mut j = 0;
        let mut y = y0;
        while y <= y1 {
            let mut x = x0 + if j % 2 == 1 { 0.5 * self.spacing } else { 0.0 };
            while x <= x1 {
                // retinal coordinates: the fovea is at the origin (device axis)
                let r = x.hypot(y);
                let mut p = self.proportions;
                if r < self.s_free_radius {
                    p[2] = 0.0;
                }
                let u = rng.uniform() * (p[0] + p[1] + p[2]);
                let kind = if u < p[0] { 0 } else if u < p[0] + p[1] { 1 } else { 2 };
                cones.push(Cone { position: [x, y], kind, expected: 0.0, isomerisations: 0.0, excitation: 0.0 });
                x += self.spacing;
            }
            y += row;
            j += 1;
        }
        // bilinear interpolation of each band at the cone centres
        let at = |e: &[f64], p: [f64; 2]| -> f64 {
            let fx = ((p[0] - xs[0]) / w.pitch).clamp(0.0, (xs.len() - 1) as f64);
            let fy = ((p[1] - ys[0]) / w.pitch).clamp(0.0, (ys.len() - 1) as f64);
            let (i, jj) = ((fx as usize).min(xs.len().saturating_sub(2)), (fy as usize).min(ys.len().saturating_sub(2)));
            let (tx, ty) = (fx - i as f64, fy - jj as f64);
            let g = |a: usize, b: usize| e[a + xs.len() * b];
            g(i, jj) * (1.0 - tx) * (1.0 - ty) + g(i + 1, jj) * tx * (1.0 - ty) + g(i, jj + 1) * (1.0 - tx) * ty + g(i + 1, jj + 1) * tx * ty
        };
        // quantal sensitivity ∝ energy sensitivity / λ, normalised to peak 1
        let peak_q: [f64; 3] = std::array::from_fn(|k| (380..=780).map(|nm| crate::color::lms(nm as f64 * 1e-9)[k] / nm as f64).fold(0.0, f64::max));
        for c in cones.iter_mut() {
            let (mut q, mut e) = (0.0, 0.0);
            for (l, band) in &image.bands {
                let irr = at(band, c.position);
                let s = crate::color::lms(*l)[c.kind];
                let sq = s / (l * 1e9) / peak_q[c.kind];
                q += irr * l / HC * self.peak_area * sq * self.integration;
                e += irr * s;
            }
            c.expected = q;
            c.isomerisations = rng.poisson(q);
            c.excitation = if q > 0.0 { e * c.isomerisations / q } else { 0.0 };
        }
        cones
    }

    /// The colour image the mosaic encodes: at each window sample, each
    /// cone class's excitation is estimated from that class's cones nearby
    /// (Gaussian weights, σ = 1.5 spacings of the class), and LMS → XYZ.
    pub fn perceived_xyz(&self, cones: &[Cone], window: &Window) -> Vec<[f64; 3]> {
        let (xs, ys) = (window.coords(0), window.coords(1));
        let m = crate::color::lms_to_xyz().0;
        let sigma: [f64; 3] = std::array::from_fn(|k| 1.5 * self.spacing / self.proportions[k].max(1e-3).sqrt());
        // bucket cones by coarse cell for neighbour search
        let cell = 3.0 * sigma.iter().cloned().fold(0.0, f64::max);
        let (gx0, gy0) = (xs[0] - cell, ys[0] - cell);
        let gw = (((xs[xs.len() - 1] - gx0) / cell).ceil() as usize) + 2;
        let gh = (((ys[ys.len() - 1] - gy0) / cell).ceil() as usize) + 2;
        let mut grid = vec![Vec::new(); gw * gh];
        for (i, c) in cones.iter().enumerate() {
            let gx = ((c.position[0] - gx0) / cell) as usize;
            let gy = ((c.position[1] - gy0) / cell) as usize;
            if gx < gw && gy < gh {
                grid[gx + gw * gy].push(i);
            }
        }
        use rayon::prelude::*;
        (0..xs.len() * ys.len())
            .into_par_iter()
            .map(|q| {
                let p = [xs[q % xs.len()], ys[q / xs.len()]];
                let (gx, gy) = (((p[0] - gx0) / cell) as isize, ((p[1] - gy0) / cell) as isize);
                let mut num = [0.0; 3];
                let mut den = [0.0; 3];
                for dy in -1..=1 {
                    for dx in -1..=1 {
                        let (x, y) = (gx + dx, gy + dy);
                        if x < 0 || y < 0 || x as usize >= gw || y as usize >= gh {
                            continue;
                        }
                        for &i in &grid[x as usize + gw * y as usize] {
                            let c = &cones[i];
                            let d2 = (c.position[0] - p[0]).powi(2) + (c.position[1] - p[1]).powi(2);
                            let wgt = (-0.5 * d2 / (sigma[c.kind] * sigma[c.kind])).exp();
                            num[c.kind] += wgt * c.excitation;
                            den[c.kind] += wgt;
                        }
                    }
                }
                let lms: [f64; 3] = std::array::from_fn(|k| if den[k] > 0.0 { num[k] / den[k] } else { 0.0 });
                std::array::from_fn(|r| (0..3).map(|k| m[r][k] * lms[k]).sum())
            })
            .collect()
    }
}

/// Convert a window of XYZ to linear sRGB.
pub fn xyz_to_rgb(xyz: &[[f64; 3]]) -> Vec<[f64; 3]> {
    xyz.iter().map(|v| crate::color::xyz_to_linear_srgb(*v)).collect()
}
