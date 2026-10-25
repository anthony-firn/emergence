//! Interference experiments: a source, openings with polarization elements,
//! and the irradiance (W/m²) they produce on a plane.
//!
//! Each incoherent source mode (wavelength × emitter point × polarization)
//! is propagated coherently through every opening by the diffraction
//! ladder (certified to the requested field accuracy); the Jones fields of
//! the openings add, an optional analyzer acts, and the intensities of the
//! modes add with their power weights.

use std::f64::consts::PI;

use emergence::{Ladder, Request};
use num_complex::Complex64 as Cx;
use rayon::prelude::*;

use crate::diffraction::{self, Field, Incident, Query, Rect};
use crate::polarization::{Jones, JonesMatrix};
use crate::source::{Mode, Source};

/// One opening and the polarization element covering it.
#[derive(Clone, Copy, Debug)]
pub struct Opening {
    pub rect: Rect,
    pub element: JonesMatrix,
}

/// An interference experiment observed on the plane `z = distance`.
pub struct Experiment {
    pub source: Source,
    pub openings: Vec<Opening>,
    pub distance: f64,
    /// Polarizer in front of the detector, if any.
    pub analyzer: Option<JonesMatrix>,
    /// Relative accuracy of each opening's field (and of the coherence
    /// sampling).
    pub tol: f64,
    ladder: Ladder<Query, Field>,
}

/// Irradiance on a set of points, with how it was computed.
#[derive(Clone, Debug)]
pub struct Irradiance {
    pub values: Vec<f64>,
    pub modes: usize,
    /// Field evaluations per ladder rung (Fraunhofer, Fresnel, line, RS).
    pub rungs: Vec<usize>,
    /// Openings narrower than the scalar domain (≥ 10 λ).
    pub scalar_ok: bool,
}

impl Experiment {
    pub fn new(source: Source, openings: Vec<Opening>, distance: f64, tol: f64) -> Self {
        Self {
            source,
            openings,
            distance,
            analyzer: None,
            tol,
            ladder: diffraction::ladder(0),
        }
    }

    /// Young's double slit: two `width × height` slits whose centres are
    /// `separation` apart along x.
    pub fn double_slit(source: Source, width: f64, separation: f64, height: f64, distance: f64, tol: f64) -> Self {
        let open = |x: f64| Opening {
            rect: Rect { center: [x, 0.0], size: [width, height] },
            element: JonesMatrix::identity(),
        };
        Self::new(source, vec![open(-0.5 * separation), open(0.5 * separation)], distance, tol)
    }

    /// The same experiment with only the openings in `keep` (blocking the
    /// others), as an experimenter does to measure each beam alone.
    pub fn with_only(&self, keep: &[usize]) -> Self {
        let openings = keep.iter().map(|&i| self.openings[i]).collect();
        let mut e = Self::new(self.source.clone(), openings, self.distance, self.tol);
        e.analyzer = self.analyzer;
        e
    }

    /// Largest distance between opening centres.
    pub fn baseline(&self) -> f64 {
        let mut b: f64 = 0.0;
        for a in &self.openings {
            for c in &self.openings {
                let d = ((a.rect.center[0] - c.rect.center[0]).powi(2) + (a.rect.center[1] - c.rect.center[1]).powi(2)).sqrt();
                b = b.max(d);
            }
        }
        b
    }

    /// Largest path difference between openings over the points.
    fn max_path(&self, points: &[[f64; 2]]) -> f64 {
        let l = self.distance;
        let mut m: f64 = 0.0;
        for p in points {
            let r: Vec<f64> = self
                .openings
                .iter()
                .map(|o| ((p[0] - o.rect.center[0]).powi(2) + (p[1] - o.rect.center[1]).powi(2) + l * l).sqrt())
                .collect();
            let (lo, hi) = r.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| (a.min(v), b.max(v)));
            m = m.max(hi - lo);
        }
        m + self.source_path()
    }

    /// Source side of the path difference: an emitter point at transverse
    /// offset ≤ s/2 adds ≤ baseline·(s/2)/z_s.
    pub fn source_path(&self) -> f64 {
        use crate::source::Extent;
        match self.source.extent {
            Extent::Slit { width, distance, .. } | Extent::Disk { diameter: width, distance, .. } => self.baseline() * 0.5 * width / distance,
            _ => 0.0,
        }
    }

    /// The incoherent source modes this experiment needs at these points.
    pub fn modes_for(&self, points: &[[f64; 2]]) -> Vec<Mode> {
        self.source.modes(self.baseline(), self.max_path(points), self.tol)
    }

    /// Irradiance (W/m²) at transverse points `(x, y)` of the observation
    /// plane.
    pub fn irradiance(&self, points: &[[f64; 2]]) -> Irradiance {
        self.irradiance_modes(points, &self.modes_for(points))
    }

    /// Irradiance from given source modes (so that runs with openings
    /// blocked use exactly the same source sampling).
    pub fn irradiance_modes(&self, points: &[[f64; 2]], modes: &[Mode]) -> Irradiance {
        // Group modes that share the scalar field (same wavelength and
        // illumination): the scalar diffraction is computed once.
        let mut groups: Vec<(f64, Incident, Vec<(Jones, f64)>)> = Vec::new();
        for m in modes {
            match groups.iter_mut().find(|g| g.0 == m.wavelength && g.1 == m.incident) {
                Some(g) => g.2.push((m.jones, m.weight)),
                None => groups.push((m.wavelength, m.incident, vec![(m.jones, m.weight)])),
            }
        }
        let n_rungs = self.ladder.len();
        let results: Vec<(f64, Vec<usize>)> = points
            .par_iter()
            .map(|&p| {
                let mut total = 0.0;
                let mut hist = vec![0usize; n_rungs];
                for (lambda, incident, pols) in &groups {
                    let fields: Vec<Cx> = self
                        .openings
                        .iter()
                        .map(|o| {
                            let q = Query {
                                aperture: o.rect,
                                incident: *incident,
                                point: [p[0], p[1], self.distance],
                                wavelength: *lambda,
                            };
                            let s = self.ladder.solve(&q, Request::new(self.tol * field_scale(&q)));
                            hist[s.rung] += 1;
                            s.value.u
                        })
                        .collect();
                    for (jones, w) in pols {
                        let mut e = Jones([Cx::new(0.0, 0.0); 2]);
                        for (o, u) in self.openings.iter().zip(&fields) {
                            e = e.add(&o.element.apply(jones).scale(*u));
                        }
                        if let Some(a) = &self.analyzer {
                            e = a.apply(&e);
                        }
                        total += w * e.intensity();
                    }
                }
                (total, hist)
            })
            .collect();
        let mut rungs = vec![0usize; n_rungs];
        let mut values = Vec::with_capacity(points.len());
        for (v, h) in results {
            values.push(v);
            for (r, c) in rungs.iter_mut().zip(h) {
                *r += c;
            }
        }
        let lambda_min = modes.iter().map(|m| m.wavelength).fold(f64::INFINITY, f64::min);
        Irradiance {
            values,
            modes: modes.len(),
            rungs,
            scalar_ok: self.openings.iter().all(|o| o.rect.scalar_validity(lambda_min) >= diffraction::MIN_WIDTH_IN_WAVELENGTHS),
        }
    }

    /// Scalar field of every opening at points of the plane `z` for one
    /// illumination (wavelength + incident wave): `fields[point][opening]`.
    pub fn opening_fields(&self, z: f64, points: &[[f64; 2]], wavelength: f64, incident: &Incident) -> Vec<Vec<Cx>> {
        points
            .par_iter()
            .map(|p| {
                self.openings
                    .iter()
                    .map(|o| {
                        let q = Query { aperture: o.rect, incident: *incident, point: [p[0], p[1], z], wavelength };
                        self.ladder.solve(&q, Request::new(self.tol * field_scale(&q))).value.u
                    })
                    .collect()
            })
            .collect()
    }

    /// True when every opening transmits both polarizations equally and no
    /// analyzer is present: the two Jones components then carry the same
    /// scalar field and polarization can be ignored.
    pub fn polarization_free(&self) -> bool {
        self.analyzer.is_none() && self.openings.iter().all(|o| o.element == JonesMatrix::identity())
    }

    /// Fields of every opening on the grid `xs × ys` of the plane `z`
    /// (`out[opening][i + xs.len() j]`) by the band-limited angular
    /// spectrum ([`crate::angular`]): the fast rung for whole planes.
    ///
    /// The band holds every direction from an opening to the grid,
    /// widened by four Fresnel zones `√(λz)` (the stationary region of the
    /// frequency integral), with a raised-cosine roll-off over four more.
    /// `audits` grid points (evenly strided) are recomputed with the
    /// certified ladder; the worst difference relative to the field scale
    /// is returned with the fields.
    pub fn grid_fields(&self, z: f64, xs: &[f64], ys: &[f64], wavelength: f64, incident: &Incident, audits: usize) -> (Vec<Vec<Cx>>, f64) {
        use crate::angular::{plane_fields, Band, OpeningSpectrum};
        let zone = (wavelength * z).sqrt();
        let sine = |d: f64| d / (d * d + z * z).sqrt();
        let grid = [xs, ys];
        let mut band = Band { center: [0.0; 2], inner: [0.0; 2], outer: [0.0; 2] };
        let mut period = [0.0; 2];
        for a in 0..2 {
            let (gmin, gmax) = grid[a].iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| (lo.min(v), hi.max(v)));
            let lo_o = self.openings.iter().map(|o| o.rect.center[a] - 0.5 * o.rect.size[a]).fold(f64::INFINITY, f64::min);
            let hi_o = self.openings.iter().map(|o| o.rect.center[a] + 0.5 * o.rect.size[a]).fold(f64::NEG_INFINITY, f64::max);
            let f_lo = sine(gmin - hi_o - 4.0 * zone) / wavelength;
            let f_hi = sine(gmax - lo_o + 4.0 * zone) / wavelength;
            let taper = sine(4.0 * zone) / wavelength;
            band.center[a] = 0.5 * (f_lo + f_hi);
            band.inner[a] = 0.5 * (f_hi - f_lo);
            band.outer[a] = band.inner[a] + taper;
            // replicas of the band-limited field clear of the grid
            let ring = 4.0 / taper;
            let reach = |f: f64| wavelength * z * f / (1.0 - (wavelength * f).powi(2)).max(1e-6).sqrt();
            let (flo, fhi) = (band.center[a] - band.outer[a], band.center[a] + band.outer[a]);
            let (slo, shi) = (lo_o + reach(flo) - ring, hi_o + reach(fhi) + ring);
            period[a] = (shi - gmin).max(gmax - slo);
        }
        let spectra: Vec<OpeningSpectrum> = self.openings.iter().map(|o| OpeningSpectrum::new(&o.rect, incident, wavelength)).collect();
        let (fields, _) = plane_fields(&spectra, wavelength, z, &band, period, xs, ys);
        let n = xs.len() * ys.len();
        let mut worst: f64 = 0.0;
        for t in 0..audits.min(n) {
            let q = (t * n) / audits.min(n) + (n / audits.min(n).max(1)) / 2;
            let q = q.min(n - 1);
            let p = [xs[q % xs.len()], ys[q / xs.len()], z];
            for (o, f) in self.openings.iter().zip(&fields) {
                let qr = Query { aperture: o.rect, incident: *incident, point: p, wavelength };
                let exact = self.ladder.solve(&qr, Request::new(0.1 * self.tol * field_scale(&qr))).value.u;
                worst = worst.max((f[q] - exact).norm() / field_scale(&qr));
            }
        }
        (fields, worst)
    }

    /// Spectral irradiance on the grid `xs × ys` of the observation plane
    /// from [`Experiment::grid_fields`], per wavelength, with the worst
    /// audit error over all illuminations.
    pub fn grid_spectral_irradiance(&self, xs: &[f64], ys: &[f64], modes: &[Mode], audits: usize) -> (Vec<(f64, Vec<f64>)>, f64) {
        let n = xs.len() * ys.len();
        let mut out: Vec<(f64, Vec<f64>)> = Vec::new();
        let mut worst: f64 = 0.0;
        let mut groups: Vec<(f64, Incident, Vec<(Jones, f64)>)> = Vec::new();
        for m in modes {
            match groups.iter_mut().find(|g| g.0 == m.wavelength && g.1 == m.incident) {
                Some(g) => g.2.push((m.jones, m.weight)),
                None => groups.push((m.wavelength, m.incident, vec![(m.jones, m.weight)])),
            }
        }
        for (lambda, incident, pols) in &groups {
            let (fields, err) = self.grid_fields(self.distance, xs, ys, *lambda, incident, audits);
            worst = worst.max(err);
            let acc: Vec<f64> = (0..n)
                .into_par_iter()
                .map(|q| {
                    pols.iter()
                        .map(|(jones, w)| {
                            let mut e = Jones([Cx::new(0.0, 0.0); 2]);
                            for (o, f) in self.openings.iter().zip(&fields) {
                                e = e.add(&o.element.apply(jones).scale(f[q]));
                            }
                            if let Some(a) = &self.analyzer {
                                e = a.apply(&e);
                            }
                            w * e.intensity()
                        })
                        .sum()
                })
                .collect();
            match out.iter_mut().find(|b| b.0 == *lambda) {
                Some(b) => b.1.iter_mut().zip(&acc).for_each(|(a, v)| *a += v),
                None => out.push((*lambda, acc)),
            }
        }
        (out, worst)
    }

    /// Irradiance per wavelength: `(λ, E_λ(points))`, each already weighted
    /// by its share of the source power (for detectors, whose response
    /// depends on photon energy and quantum efficiency).
    pub fn spectral_irradiance(&self, points: &[[f64; 2]], modes: &[Mode]) -> Vec<(f64, Vec<f64>)> {
        let mut lambdas: Vec<f64> = Vec::new();
        for m in modes {
            if !lambdas.contains(&m.wavelength) {
                lambdas.push(m.wavelength);
            }
        }
        lambdas
            .into_iter()
            .map(|l| {
                let sub: Vec<Mode> = modes.iter().filter(|m| m.wavelength == l).copied().collect();
                (l, self.irradiance_modes(points, &sub).values)
            })
            .collect()
    }

    /// Shortest fringe period on the observation plane (`λ_min L/baseline`).
    pub fn fringe_scale(&self, modes: &[Mode]) -> f64 {
        let lmin = modes.iter().map(|m| m.wavelength).fold(f64::INFINITY, f64::min);
        let b = self.baseline();
        if b > 0.0 {
            lmin * self.distance / b
        } else {
            // a single opening: its own diffraction scale
            let a = self.openings.iter().map(|o| o.rect.size[0].min(o.rect.size[1])).fold(f64::INFINITY, f64::min);
            lmin * self.distance / a
        }
    }

    /// Normalised mutual coherence at points, measured the experimenter's
    /// way: `γ = (I − ΣI_k)/(2√(I₁I₂))` from the full pattern and each
    /// opening alone (two openings).
    pub fn coherence(&self, points: &[[f64; 2]]) -> Vec<f64> {
        assert_eq!(self.openings.len(), 2);
        let modes = self.modes_for(points);
        let i = self.irradiance_modes(points, &modes).values;
        let i1 = self.with_only(&[0]).irradiance_modes(points, &modes).values;
        let i2 = self.with_only(&[1]).irradiance_modes(points, &modes).values;
        i.iter().zip(i1.iter().zip(&i2)).map(|(t, (a, b))| (t - a - b) / (2.0 * (a * b).sqrt())).collect()
    }
}

/// Field scale of one opening: `|A| area/(λ L)` (`/z_s` for a point
/// source), the on-axis far-field magnitude.
pub fn field_scale(q: &Query) -> f64 {
    let a = match q.incident {
        Incident::Plane { amplitude, .. } => amplitude.norm(),
        Incident::Point { amplitude, distance, .. } => amplitude.norm() / distance,
    };
    a * q.aperture.area() / (q.wavelength * q.point[2])
}

/// Wavelength-to-wavenumber helper.
pub fn k_of(lambda: f64) -> f64 {
    2.0 * PI / lambda
}

impl Mode {
    pub fn k(&self) -> f64 {
        k_of(self.wavelength)
    }
}
