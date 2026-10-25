//! Partially coherent sources: spectrum, extent and polarization.
//!
//! A real source is an incoherent superposition of mutually incoherent
//! modes: wavelengths (temporal coherence), points on the emitter (spatial
//! coherence) and orthogonal polarization states. Each mode is propagated
//! coherently; intensities add. The number of modes is set by coherence,
//! not by a fixed resolution:
//!
//! * **Spectrum.** Sampled uniformly in wavenumber `σ = 1/λ`. A uniform
//!   sampling `δσ` makes the computed interference pattern periodic in path
//!   difference with period `1/δσ` (spurious "revivals"), so `δσ` is chosen
//!   so that this period exceeds twice the largest path difference `Δ_max`
//!   in the scene. When the coherence length exceeds `Δ_max` by enough that
//!   the visibility loss is below tolerance, one wavelength suffices (the
//!   monochromatic fold).
//! * **Extent.** A source of width `s` at distance `z_s` from two openings
//!   `d` apart shifts the fringes by up to `π s d/(λ z_s)`; the source is
//!   sampled finely enough to resolve that phase (`≤ π/8` per sample), or
//!   treated as a point when the visibility loss `(π s d/(λ z_s))²/6` is
//!   below tolerance (the coherent fold).
//! * **Polarization.** Polarized sources have one Jones vector; unpolarized
//!   light is the incoherent sum of two orthogonal states.

use std::f64::consts::PI;

use num_complex::Complex64 as Cx;

use crate::polarization::Jones;

/// Spectral power distribution.
#[derive(Clone, Debug, PartialEq)]
pub enum Spectrum {
    /// A single line (laser), wavelength in metres.
    Line(f64),
    /// Gaussian in wavelength: centre and full width at half maximum.
    Gaussian { center: f64, fwhm: f64 },
    /// Planck blackbody at temperature `kelvin` through a Gaussian
    /// band-pass filter (centre, FWHM).
    FilteredBlackbody { kelvin: f64, center: f64, fwhm: f64 },
    /// Sum of Gaussian bands `(relative power, centre, FWHM)`, e.g. a
    /// phosphor-converted white LED.
    Mixture(Vec<(f64, f64, f64)>),
}

/// One spectral sample: wavelength and power fraction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralSample {
    pub wavelength: f64,
    pub weight: f64,
}

/// Planck spectral radiance shape `λ⁻⁵/(e^{hc/λkT} − 1)` (arbitrary units).
pub fn planck(lambda: f64, kelvin: f64) -> f64 {
    // CODATA 2018: h, c, k_B
    const HC_OVER_K: f64 = 6.626_070_15e-34 * 2.997_924_58e8 / 1.380_649e-23;
    lambda.powi(-5) / ((HC_OVER_K / (lambda * kelvin)).exp() - 1.0)
}

impl Spectrum {
    /// Mean wavelength (power weighted, approximately the centre).
    pub fn center(&self) -> f64 {
        match self {
            Spectrum::Line(l) => *l,
            Spectrum::Gaussian { center, .. } | Spectrum::FilteredBlackbody { center, .. } => *center,
            Spectrum::Mixture(b) => b.iter().map(|x| x.0 * x.1).sum::<f64>() / b.iter().map(|x| x.0).sum::<f64>(),
        }
    }

    /// Shortest and longest wavelengths carrying power (±5 standard
    /// deviations of the bands).
    pub fn span(&self) -> (f64, f64) {
        let sd = |w: f64| w / (2.0 * (2.0 * 2f64.ln()).sqrt());
        match self {
            Spectrum::Line(l) => (*l, *l),
            Spectrum::Gaussian { center, fwhm } | Spectrum::FilteredBlackbody { center, fwhm, .. } => (center - 5.0 * sd(*fwhm), center + 5.0 * sd(*fwhm)),
            Spectrum::Mixture(b) => (
                b.iter().map(|x| x.1 - 5.0 * sd(x.2)).fold(f64::INFINITY, f64::min).max(380e-9),
                b.iter().map(|x| x.1 + 5.0 * sd(x.2)).fold(0.0, f64::max).min(780e-9),
            ),
        }
    }

    /// Narrowest band (σ standard deviation in wavenumber) for sampling.
    fn narrowest_sigma_std(&self) -> f64 {
        let sd = |w: f64, c: f64| w / (2.0 * (2.0 * 2f64.ln()).sqrt()) / (c * c);
        match self {
            Spectrum::Line(_) => 0.0,
            Spectrum::Gaussian { center, fwhm } | Spectrum::FilteredBlackbody { center, fwhm, .. } => sd(*fwhm, *center),
            Spectrum::Mixture(b) => b.iter().map(|x| sd(x.2, x.1)).fold(f64::INFINITY, f64::min),
        }
    }

    /// Spectral power density per unit wavelength (relative).
    pub fn density(&self, lambda: f64) -> f64 {
        let gauss = |c: f64, w: f64| {
            let s = w / (2.0 * (2.0 * 2f64.ln()).sqrt());
            (-0.5 * ((lambda - c) / s).powi(2)).exp() / s
        };
        match self {
            Spectrum::Line(_) => 0.0,
            Spectrum::Gaussian { center, fwhm } => gauss(*center, *fwhm),
            Spectrum::FilteredBlackbody { kelvin, center, fwhm } => planck(lambda, *kelvin) * gauss(*center, *fwhm),
            Spectrum::Mixture(b) => b.iter().map(|x| x.0 * gauss(x.1, x.2)).sum(),
        }
    }

    /// Wavenumber FWHM `Δσ ≈ Δλ/λ²` (0 for a line).
    pub fn sigma_width(&self) -> f64 {
        match self {
            Spectrum::Line(_) => 0.0,
            Spectrum::Gaussian { center, fwhm } | Spectrum::FilteredBlackbody { center, fwhm, .. } => fwhm / (center * center),
            Spectrum::Mixture(_) => {
                let (lo, hi) = self.span();
                1.0 / lo - 1.0 / hi
            }
        }
    }

    /// Visibility of fringes at path difference `Δ` for a Gaussian line
    /// (the coherence function): `exp(−(π Δ Δσ)²/(4 ln 2))`.
    pub fn gaussian_visibility(&self, path_difference: f64) -> f64 {
        let x = PI * path_difference * self.sigma_width();
        (-(x * x) / (4.0 * 2f64.ln())).exp()
    }

    /// Spectral samples resolving path differences up to `max_path` with
    /// visibility error below `tol` (see module docs).
    pub fn samples(&self, max_path: f64, tol: f64) -> Vec<SpectralSample> {
        if let Spectrum::Line(l) = *self {
            return vec![SpectralSample { wavelength: l, weight: 1.0 }];
        }
        // Monochromatic fold (single Gaussian lines only).
        if matches!(self, Spectrum::Gaussian { .. } | Spectrum::FilteredBlackbody { .. }) && 1.0 - self.gaussian_visibility(max_path) < tol {
            return vec![SpectralSample { wavelength: self.center(), weight: 1.0 }];
        }
        let (lo, hi) = self.span();
        let (s_lo, s_hi) = (1.0 / hi, 1.0 / lo);
        // Two conditions: revival period 1/δσ ≥ 4 Δ_max, and every band
        // resolved (δσ ≤ σ_std/2, where the trapezoid rule on a smooth line
        // is exponentially accurate).
        let step = (1.0 / (4.0 * max_path.max(1e-12))).min(0.5 * self.narrowest_sigma_std());
        let n = (((s_hi - s_lo) / step).ceil() as usize).max(3) | 1;
        let mut out: Vec<SpectralSample> = (0..n)
            .map(|i| {
                let sig = s_lo + (s_hi - s_lo) * i as f64 / (n - 1) as f64;
                let lambda = 1.0 / sig;
                // density per unit σ: ρ(λ) |dλ/dσ| = ρ(λ) λ²
                SpectralSample { wavelength: lambda, weight: self.density(lambda) * lambda * lambda }
            })
            .filter(|s| s.weight > 0.0)
            .collect();
        let total: f64 = out.iter().map(|s| s.weight).sum();
        out.iter_mut().for_each(|s| s.weight /= total);
        out
    }
}

/// Emitter geometry, transverse to the axis at distance `distance` before
/// the aperture plane, or a collimated beam.
#[derive(Clone, Debug, PartialEq)]
pub enum Extent {
    /// Collimated beam of uniform irradiance (W/m²) at normal incidence.
    Collimated { irradiance: f64 },
    /// Point source of total power (W), isotropic.
    Point { power: f64, distance: f64 },
    /// Incoherent uniform slit source of width `width` along x (and point-
    /// like along y), total power.
    Slit { power: f64, distance: f64, width: f64 },
    /// Incoherent uniform disk of diameter `diameter`, total power.
    Disk { power: f64, distance: f64, diameter: f64 },
}

/// A fully described source.
#[derive(Clone, Debug, PartialEq)]
pub struct Source {
    pub spectrum: Spectrum,
    pub extent: Extent,
    /// `None`: unpolarized.
    pub polarization: Option<Jones>,
}

/// One mutually incoherent mode: an illumination with its power fraction.
#[derive(Clone, Copy, Debug)]
pub struct Mode {
    pub wavelength: f64,
    pub incident: crate::diffraction::Incident,
    pub jones: Jones,
    pub weight: f64,
}

impl Source {
    /// The incoherent modes needed to resolve fringes from openings up to
    /// `baseline` apart over path differences up to `max_path`, with
    /// visibility errors below `tol`.
    pub fn modes(&self, baseline: f64, max_path: f64, tol: f64) -> Vec<Mode> {
        use crate::diffraction::Incident;
        let spectral = self.spectrum.samples(max_path, tol);
        let pols: Vec<(Jones, f64)> = match self.polarization {
            Some(j) => vec![(j.normalized(), 1.0)],
            None => vec![(Jones::horizontal(), 0.5), (Jones::vertical(), 0.5)],
        };
        let mut out = Vec::new();
        for sp in &spectral {
            let lambda = sp.wavelength;
            let points: Vec<([f64; 2], f64, f64, f64)> = match self.extent {
                Extent::Collimated { irradiance } => vec![([0.0; 2], 1.0, irradiance, 0.0)],
                Extent::Point { power, distance } => vec![([0.0; 2], 1.0, power / (4.0 * PI), distance)],
                Extent::Slit { power, distance, width } | Extent::Disk { power, distance, diameter: width } => {
                    // Midpoint sampling of the source: the visibility error of
                    // a step δφ in fringe phase is ≈ δφ²/24, so δφ = √(24 tol)
                    // (at most π/8).
                    let phase = PI * width * baseline / (lambda * distance);
                    let step = (24.0 * tol).sqrt().min(PI / 8.0);
                    let n = if phase * phase / 6.0 < tol { 1 } else { ((phase / step).ceil() as usize).max(2) };
                    let disk = matches!(self.extent, Extent::Disk { .. });
                    let mut pts = Vec::new();
                    for i in 0..n {
                        for j in 0..(if disk { n } else { 1 }) {
                            let x = if n == 1 { 0.0 } else { -0.5 * width + width * (i as f64 + 0.5) / n as f64 };
                            let y = if disk && n > 1 { -0.5 * width + width * (j as f64 + 0.5) / n as f64 } else { 0.0 };
                            if disk && x * x + y * y > 0.25 * width * width && n > 1 {
                                continue;
                            }
                            pts.push(([x, y], 1.0, power / (4.0 * PI), distance));
                        }
                    }
                    let m = pts.len() as f64;
                    pts.into_iter().map(|(p, _, a2, d)| (p, 1.0 / m, a2, d)).collect()
                }
            };
            for (pos, frac, a2, distance) in points {
                for &(jones, pf) in &pols {
                    let weight = sp.weight * frac * pf;
                    let amplitude = Cx::new(a2.sqrt(), 0.0);
                    let incident = if distance == 0.0 {
                        Incident::Plane { amplitude, sines: [0.0; 2] }
                    } else {
                        Incident::Point { amplitude, position: pos, distance }
                    };
                    out.push(Mode { wavelength: lambda, incident, jones, weight });
                }
            }
        }
        out
    }
}
