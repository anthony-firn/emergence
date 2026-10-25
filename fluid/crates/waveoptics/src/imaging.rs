//! Image formation by imaging devices: from the complex field arriving at
//! the device's entrance pupil to the field on its detector.
//!
//! For a device focused at infinity (object-space focal length `f`), a
//! plane-wave component at angle `θ` comes to focus at `x' = f θ`, and the
//! detector field is the Fraunhofer transform of the pupil-weighted
//! incoming field (Goodman, *Introduction to Fourier Optics*, §6.1):
//!
//! ```text
//! U(x') = 1/(iλf) ∫∫ P(u) U_in(u) e^{−i 2π u·x'/(λ f)} d²u,
//! P(u) = A(u) e^{i k W(u)},
//! ```
//!
//! with `W` the wavefront aberration traced through the real optics at the
//! same wavelength (so chromatic focus errors, spherical aberration and
//! defocus of the detector all appear), and the incoming field's own
//! curvature carrying the object distance. A lensless pinhole is the same
//! integral with `W(u) = |u|²/(2f)` (the Fresnel kernel to a detector at
//! distance `f`), which reproduces both its geometric shadow and its
//! diffraction.
//!
//! The transform is a matrix (non-FFT) DFT onto an arbitrary detector
//! window (Soummer et al., Opt. Express 15, 15935 (2007)), separable in
//! x and y. Sampling: the pupil grid must resolve the phase of `P·U_in`,
//! and the window must lie within the alias-free width `λ f/h`
//! ([`PupilFunction::alias_free_width`]).

use std::f64::consts::PI;

use num_complex::Complex64 as Cx;
use rayon::prelude::*;

use crate::raytrace::System;

/// A sampled pupil: grid of `n × n` cell centres over `[−r, r]²` of the
/// entrance pupil, centred at `center` (transverse, object space).
#[derive(Clone, Debug)]
pub struct PupilFunction {
    pub n: usize,
    pub radius: f64,
    pub center: [f64; 2],
    /// Distance of the entrance-pupil plane from the device's front (first
    /// vertex, or the pinhole plane), positive inward (m).
    pub plane: f64,
    pub values: Vec<Cx>,
    /// Object-space focal length (detector distance for a pinhole).
    pub focal: f64,
    pub wavelength: f64,
}

impl PupilFunction {
    pub fn spacing(&self) -> f64 {
        2.0 * self.radius / self.n as f64
    }

    /// Transverse sample positions.
    pub fn points(&self) -> Vec<[f64; 2]> {
        let h = self.spacing();
        (0..self.n * self.n)
            .map(|q| {
                let (i, j) = (q % self.n, q / self.n);
                [self.center[0] - self.radius + h * (i as f64 + 0.5), self.center[1] - self.radius + h * (j as f64 + 0.5)]
            })
            .collect()
    }

    /// Width of detector window free of DFT aliasing, `λ f/h`.
    pub fn alias_free_width(&self) -> f64 {
        self.wavelength * self.focal / self.spacing()
    }

    /// Ideal aberration-free circular pupil (radius `a`) with an optional
    /// defocus `W(ρ) = w20 (r/a)²`.
    pub fn ideal(a: f64, focal: f64, wavelength: f64, n: usize, w20: f64) -> Self {
        let mut p = Self { n, radius: a, center: [0.0; 2], plane: 0.0, values: vec![], focal, wavelength };
        let k = 2.0 * PI / wavelength;
        p.values = p
            .points()
            .iter()
            .map(|u| {
                let r2 = u[0] * u[0] + u[1] * u[1];
                if r2 <= a * a {
                    Cx::from_polar(1.0, k * w20 * r2 / (a * a))
                } else {
                    Cx::new(0.0, 0.0)
                }
            })
            .collect();
        p
    }

    /// Lensless pinhole of diameter `d` with the detector at distance `f`:
    /// the Fresnel phase `k|u|²/(2f)` over the opening.
    pub fn pinhole(d: f64, f: f64, wavelength: f64, n: usize) -> Self {
        let a = 0.5 * d;
        let mut p = Self::ideal(a, f, wavelength, n, 0.0);
        let k = 2.0 * PI / wavelength;
        let pts = p.points();
        for (v, u) in p.values.iter_mut().zip(pts) {
            if v.norm() > 0.0 {
                *v = Cx::from_polar(1.0, k * (u[0] * u[0] + u[1] * u[1]) / (2.0 * f));
            }
        }
        p
    }

    /// Pupil of a real system at one wavelength, from its traced wavefront
    /// (vignetted samples are opaque).
    pub fn traced(system: &System, wavelength: f64, n: usize) -> Self {
        let (par, w) = system.wavefront(wavelength, n);
        let k = 2.0 * PI / wavelength;
        let values = w.iter().map(|w| w.map_or(Cx::new(0.0, 0.0), |w| Cx::from_polar(1.0, k * w))).collect();
        Self { n, radius: par.ep_radius, center: [0.0; 2], plane: par.ep_z, values, focal: par.efl, wavelength }
    }

    /// RMS wavefront error (m) over the open pupil, piston removed.
    pub fn rms_wavefront(&self) -> f64 {
        let k = 2.0 * PI / self.wavelength;
        let open: Vec<f64> = self.values.iter().filter(|v| v.norm() > 0.0).map(|v| v.arg() / k).collect();
        let m = open.iter().sum::<f64>() / open.len() as f64;
        (open.iter().map(|w| (w - m).powi(2)).sum::<f64>() / open.len() as f64).sqrt()
    }

    /// Strehl ratio `|Σ P|²/(Σ |P|)²` (on-axis intensity relative to the
    /// aberration-free pupil).
    pub fn strehl(&self) -> f64 {
        let s: Cx = self.values.iter().sum();
        let a: f64 = self.values.iter().map(|v| v.norm()).sum();
        s.norm_sqr() / (a * a)
    }
}

/// A rectangular window of detector samples.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Window {
    pub center: [f64; 2],
    pub pixels: [usize; 2],
    pub pitch: f64,
}

impl Window {
    pub fn coords(&self, axis: usize) -> Vec<f64> {
        (0..self.pixels[axis]).map(|i| self.center[axis] + self.pitch * (i as f64 + 0.5 - 0.5 * self.pixels[axis] as f64)).collect()
    }
}

/// Detector field `U(x')` from the pupil function and the incoming field
/// sampled at [`PupilFunction::points`]. Power is conserved:
/// `Σ|U|² pitch² → Σ|P U_in|² h²` when the window holds the image.
pub fn image_field(pupil: &PupilFunction, incoming: &[Cx], window: &Window) -> Vec<Cx> {
    let n = pupil.n;
    let h = pupil.spacing();
    let scale = 1.0 / (pupil.wavelength * pupil.focal);
    let g: Vec<Cx> = pupil.values.iter().zip(incoming).map(|(p, u)| p * u).collect();
    let u: Vec<f64> = (0..n).map(|i| -pupil.radius + h * (i as f64 + 0.5)).collect();
    let xs = window.coords(0);
    let ys = window.coords(1);
    let ex: Vec<Vec<Cx>> = xs.iter().map(|&x| u.iter().map(|&ui| Cx::from_polar(1.0, -2.0 * PI * ui * x * scale)).collect()).collect();
    let ey: Vec<Vec<Cx>> = ys.iter().map(|&y| u.iter().map(|&vj| Cx::from_polar(1.0, -2.0 * PI * vj * y * scale)).collect()).collect();
    // A[j][xi] = Σ_i g[j][i] e^{...x}
    let a: Vec<Vec<Cx>> = (0..n)
        .into_par_iter()
        .map(|j| xs.iter().enumerate().map(|(xi, _)| (0..n).map(|i| g[i + n * j] * ex[xi][i]).sum()).collect())
        .collect();
    let pref = Cx::new(0.0, -1.0) * (h * h * scale);
    (0..ys.len() * xs.len())
        .into_par_iter()
        .map(|q| {
            let (xi, yi) = (q % xs.len(), q / xs.len());
            pref * (0..n).map(|j| a[j][xi] * ey[yi][j]).sum::<Cx>()
        })
        .collect()
}

/// Zero-padded 2-D FFT workspace for the isoplanatic imagers.
struct Fft2 {
    px: usize,
    py: usize,
    fwd: [std::sync::Arc<dyn rustfft::Fft<f64>>; 2],
    inv: [std::sync::Arc<dyn rustfft::Fft<f64>>; 2],
}

impl Fft2 {
    /// Padding: a quarter of the grid each way holds the spread functions.
    fn new(nx: usize, ny: usize) -> Self {
        let (px, py) = ((nx + nx / 4).next_power_of_two(), (ny + ny / 4).next_power_of_two());
        let mut planner = rustfft::FftPlanner::new();
        Self { px, py, fwd: [planner.plan_fft_forward(px), planner.plan_fft_forward(py)], inv: [planner.plan_fft_inverse(px), planner.plan_fft_inverse(py)] }
    }

    fn run(&self, g: &mut [Cx], inverse: bool) {
        let (ax, ay) = if inverse { (&self.inv[0], &self.inv[1]) } else { (&self.fwd[0], &self.fwd[1]) };
        let (px, py) = (self.px, self.py);
        g.par_chunks_mut(px).for_each(|row| ax.process(row));
        // columns, transposed in blocks
        let mut t = vec![Cx::new(0.0, 0.0); px * py];
        t.par_chunks_mut(py).enumerate().for_each(|(i, col)| {
            for j in 0..py {
                col[j] = g[i + px * j];
            }
            ay.process(col);
        });
        g.par_chunks_mut(px).enumerate().for_each(|(j, row)| {
            for i in 0..px {
                row[i] = t[j + py * i];
            }
        });
    }

    /// Pupil function sampled at `u = λ D f` of every grid frequency.
    fn transfer(&self, delta: f64, wavelength: f64, distance: f64, pupil: &(impl Fn([f64; 2]) -> Cx + Sync)) -> Vec<Cx> {
        let (px, py) = (self.px, self.py);
        let freq = |i: usize, n: usize| (if i < n / 2 { i as f64 } else { i as f64 - n as f64 }) / (n as f64 * delta);
        (0..px * py).into_par_iter().map(|q| pupil([wavelength * distance * freq(q % px, px), wavelength * distance * freq(q / px, py)])).collect()
    }
}

/// Coherent image of a sampled object field by an isoplanatic imager at
/// distance `distance` (paraxial; Goodman §5.3): the object field, with its
/// quadratic phase `k|ξ|²/(2D)`, is filtered by the pupil function
/// `pupil(u)` evaluated at the pupil point `u = λ D f_ξ` of each spatial
/// frequency `f_ξ`. Returns the image intensity in object coordinates
/// (the caller applies the magnification). The grid is `nx × ny` with
/// spacing `delta`, zero-padded to a power of two in each axis.
///
/// A camera focused on the object has `pupil(u) = A(u) e^{ikW(u)}` with its
/// aberrations `W`; a lensless pinhole of radius `a` with the detector at
/// `f` is the same with `W = |u|²(1/D + 1/f)/2` (the Fresnel phases of both
/// propagations over the opening).
pub fn isoplanatic_image(object: &[Cx], nx: usize, ny: usize, delta: f64, wavelength: f64, distance: f64, pupil: impl Fn([f64; 2]) -> Cx + Sync) -> Vec<f64> {
    let fft = Fft2::new(nx, ny);
    let px = fft.px;
    let k = 2.0 * PI / wavelength;
    let mut g = vec![Cx::new(0.0, 0.0); fft.px * fft.py];
    for j in 0..ny {
        for i in 0..nx {
            let x = (i as f64 - 0.5 * nx as f64) * delta;
            let y = (j as f64 - 0.5 * ny as f64) * delta;
            g[i + px * j] = object[i + nx * j] * Cx::from_polar(1.0, k * (x * x + y * y) / (2.0 * distance));
        }
    }
    fft.run(&mut g, false);
    let h = fft.transfer(delta, wavelength, distance, &pupil);
    g.par_iter_mut().zip(&h).for_each(|(v, h)| *v *= h);
    fft.run(&mut g, true);
    let norm = 1.0 / (fft.px * fft.py) as f64;
    (0..nx * ny).map(|q| (g[q % nx + px * (q / nx)] * norm).norm_sqr()).collect()
}

/// The mean of [`isoplanatic_image`] over random-phase objects of
/// intensity `intensity` (fully developed speckle averaged out): the
/// intensity convolved with the discrete `|h|²` of the same imager — the
/// incoherent image.
pub fn incoherent_image(intensity: &[f64], nx: usize, ny: usize, delta: f64, wavelength: f64, distance: f64, pupil: impl Fn([f64; 2]) -> Cx + Sync) -> Vec<f64> {
    let fft = Fft2::new(nx, ny);
    let (px, py) = (fft.px, fft.py);
    let n = (px * py) as f64;
    // amplitude spread function h = IFFT(P)/N, then |h|² and its transform
    let mut psf = fft.transfer(delta, wavelength, distance, &pupil);
    fft.run(&mut psf, true);
    psf.par_iter_mut().for_each(|v| *v = Cx::new((*v / n).norm_sqr(), 0.0));
    fft.run(&mut psf, false);
    let mut g = vec![Cx::new(0.0, 0.0); px * py];
    for j in 0..ny {
        for i in 0..nx {
            g[i + px * j] = Cx::new(intensity[i + nx * j], 0.0);
        }
    }
    fft.run(&mut g, false);
    g.par_iter_mut().zip(&psf).for_each(|(v, o)| *v *= o);
    fft.run(&mut g, true);
    (0..nx * ny).map(|q| g[q % nx + px * (q / nx)].re / n).collect()
}
