//! Deterministic CPU render of grey-LTE continuum and a separate magnetic
//! diagnostic. Magnetic colors/curves have no radiometric interpretation.
use crate::{photosphere, surface::Surface, trace::Line};
use gravity::Vec3;
use radiation::render::Image;
use rayon::prelude::*;
use std::f64::consts::TAU;

#[derive(Clone, Copy)]
pub struct View {
    pub longitude: f64,
    pub latitude: f64,
    /// Half the vertical field of view, in photospheric radii.
    pub extent: f64,
}
impl View {
    fn basis(self) -> Result<[Vec3; 3], &'static str> {
        if !self.longitude.is_finite()
            || !self.latitude.is_finite()
            || self.latitude.abs() > std::f64::consts::FRAC_PI_2
            || !self.extent.is_finite()
            || self.extent <= 0.0
        {
            return Err("invalid solar view");
        }
        let (s, c) = self.longitude.sin_cos();
        let (sb, cb) = self.latitude.sin_cos();
        Ok([
            Vec3([-s, c, 0.0]),
            Vec3([-sb * c, -sb * s, cb]),
            Vec3([cb * c, cb * s, sb]),
        ])
    }
    /// Orthographic projected coordinates and depth, in photospheric radii.
    pub fn project(self, point: Vec3) -> Result<[f64; 3], &'static str> {
        if point.0.iter().any(|v| !v.is_finite()) {
            return Err("invalid projected point");
        }
        Ok(self.basis()?.map(|v| v.dot(point)))
    }
    /// The opaque photosphere hides both rear footpoints and far-side loops.
    pub fn visible(self, point: Vec3) -> Result<bool, &'static str> {
        let [u, v, depth] = self.project(point)?;
        Ok(u * u + v * v >= 1.0 || depth >= (1.0 - u * u - v * v).max(0.0).sqrt() - 1e-12)
    }
}
#[derive(Clone, Copy)]
pub struct Options {
    pub width: usize,
    pub height: usize,
    pub samples_per_axis: usize,
    pub view: View,
}
impl Options {
    fn validate(self) -> Result<(), &'static str> {
        self.view.basis()?;
        let pixels = self
            .width
            .checked_mul(self.height)
            .ok_or("solar image overflow")?;
        if self.width == 0
            || self.height == 0
            || pixels > 20_000_000
            || self.samples_per_axis == 0
            || self.samples_per_axis > 8
        {
            return Err("invalid solar image dimensions/sampling");
        }
        Ok(())
    }
    fn scale(self) -> f64 {
        self.height as f64 / (2.0 * self.view.extent)
    }
}
fn image(options: Options, rgb: Vec<[f64; 3]>) -> Image {
    Image {
        width: options.width,
        height: options.height,
        rgb,
        // Deterministic quadrature has no Monte Carlo standard error estimate.
        standard_error_y: vec![f64::NAN; options.width * options.height],
        truncated_paths: 0,
        samples_per_pixel: vec![options.samples_per_axis.pow(2); options.width * options.height],
        pilot_paths: 0,
    }
}
fn disk(
    options: Options,
    value: impl Fn(Vec3, f64) -> [f64; 3] + Sync,
) -> Result<Image, &'static str> {
    options.validate()?;
    let [right, up, normal] = options.view.basis()?;
    let scale = options.scale();
    let n = options.samples_per_axis;
    let rgb = (0..options.width * options.height)
        .into_par_iter()
        .map(|index| {
            let mut out = [0.0; 3];
            for sy in 0..n {
                for sx in 0..n {
                    let x = (index % options.width) as f64 + (sx as f64 + 0.5) / n as f64;
                    let y = (index / options.width) as f64 + (sy as f64 + 0.5) / n as f64;
                    let u = (x - options.width as f64 * 0.5) / scale;
                    let v = (options.height as f64 * 0.5 - y) / scale;
                    if u * u + v * v < 1.0 {
                        let mu = (1.0 - u * u - v * v).sqrt();
                        let rgb = value(right * u + up * v + normal * mu, mu);
                        for (a, b) in out.iter_mut().zip(rgb) {
                            *a += b / (n * n) as f64;
                        }
                    }
                }
            }
            out
        })
        .collect();
    Ok(image(options, rgb))
}
/// Spectrally integrated continuum radiance. Lookup interpolation is in angle,
/// not an invented magnetic-field-to-temperature or sunspot conversion.
pub fn continuum(effective_temperature: f64, options: Options) -> Result<Image, &'static str> {
    const INTERVALS: usize = 512;
    let table = (0..=INTERVALS)
        .map(|i| photosphere::rgb(i as f64 / INTERVALS as f64, effective_temperature))
        .collect::<Result<Vec<_>, _>>()?;
    disk(options, |_, mu| {
        let q = mu * INTERVALS as f64;
        let i = (q as usize).min(INTERVALS - 1);
        let t = q - i as f64;
        std::array::from_fn(|k| table[i][k] * (1.0 - t) + table[i + 1][k] * t)
    })
}
/// Display the actual equal-area cell value. The asinh scale is symmetric and
/// explicitly supplied in tesla. Shading assists shape perception only.
pub fn magnetic(
    surface: &Surface,
    lines: &[Line],
    field_scale: f64,
    options: Options,
) -> Result<Image, &'static str> {
    if !field_scale.is_finite() || field_scale <= 0.0 {
        return Err("invalid magnetic display scale");
    }
    let (n_mu, n_phi) = surface.dimensions();
    let mut out = disk(options, |p, mu| {
        let i = (((p.0[2] + 1.0) * 0.5 * n_mu as f64) as usize).min(n_mu - 1);
        let j =
            ((p.0[1].atan2(p.0[0]).rem_euclid(TAU) / TAU * n_phi as f64) as usize).min(n_phi - 1);
        let b = surface.field()[i * n_phi + j];
        let a = (5.0 * b.abs() / field_scale).asinh() / 5.0_f64.asinh();
        let a = a.min(1.0);
        let color = if b >= 0.0 {
            [1.0, 0.06, 0.02]
        } else {
            [0.02, 0.12, 1.0]
        };
        std::array::from_fn(|k| (0.25 * (1.0 - a) + color[k] * a) * (0.4 + 0.6 * mu))
    })?;
    let radius = surface.parameters().radius;
    for line in lines {
        let color = if line.closed() {
            [1.8, 0.65, 0.10]
        } else if line.open() {
            [0.06, 1.4, 0.85]
        } else {
            [0.4; 3]
        };
        for pair in line.points.windows(2) {
            let a = pair[0] / radius;
            let b = pair[1] / radius;
            let pa = options.view.project(a)?;
            let pb = options.view.project(b)?;
            let length = (pa[0] - pb[0]).hypot(pa[1] - pb[1]) * options.scale();
            let steps = (length * 2.5).ceil().max(1.0) as usize;
            for sample in 0..=steps {
                let t = sample as f64 / steps as f64;
                let p = a * (1.0 - t) + b * t;
                if !options.view.visible(p)? {
                    continue;
                }
                let [u, v, _] = options.view.project(p)?;
                let x = options.width as f64 * 0.5 + u * options.scale() - 0.5;
                let y = options.height as f64 * 0.5 - v * options.scale() - 0.5;
                splat(&mut out, x, y, color);
            }
        }
    }
    Ok(out)
}
fn splat(image: &mut Image, x: f64, y: f64, color: [f64; 3]) {
    for j in y.floor() as isize - 1..=y.floor() as isize + 2 {
        for i in x.floor() as isize - 1..=x.floor() as isize + 2 {
            if i < 0 || j < 0 || i >= image.width as isize || j >= image.height as isize {
                continue;
            }
            let distance = (i as f64 - x).hypot(j as f64 - y);
            let a = ((1.3 - distance) / 1.3).clamp(0.0, 1.0) * 0.45;
            let pixel = &mut image.rgb[j as usize * image.width + i as usize];
            for (v, c) in pixel.iter_mut().zip(color) {
                *v = *v * (1.0 - a) + c * a;
            }
        }
    }
}
