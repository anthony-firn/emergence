//! Wave optics: the rungs above rays.
//!
//! A ray tracer is a model of light with a stated domain of validity, not a
//! description of it. Light diffracts round an aperture, interferes with
//! itself after a path difference, picks up colour from a film a quarter of
//! a wavelength thick, and arrives in countable quanta. This module makes
//! all four available to the renderer, and — more importantly — makes the
//! *choice* between rays and waves an explicit, calibrated one rather than
//! an assumption.
//!
//! ```text
//!   rung 0  geometrical optics (the shadow)        cost 1
//!           indicator: Fresnel number F = a^2/(lambda z)
//!   rung 1  Fraunhofer                             cost ~1
//!   rung 2  Fresnel (exact quadratic phase)        cost ~1
//!   rung 3  Rayleigh-Sommerfeld (reference)        cost ~N^2
//! ```
//!
//! This is the structure of `analysis/emergence_ladder.md`, and the ladder
//! is an actual [`emergence::Ladder`]: ask for a tolerance and it evaluates
//! the cheapest rung whose claim meets it. Rungs 1 and 2 carry the
//! **rigorous** error bounds `waveoptics::diffraction::bounds` computes, so
//! where they apply the ladder answers with a proof rather than an
//! estimate — and where they do not, those same bounds blow up and the
//! ladder refuses them, which is the machinery earning its keep. The
//! Fraunhofer and Fresnel rungs are *far-field* models: their bounds are
//! tiny at a Fresnel number of 0.04 and astronomical at 360, so a near-field
//! query falls through to Rayleigh-Sommerfeld or back to rays.
//!
//! # The indicator that matters to a renderer
//!
//! Not the Fresnel number but the **penumbra width in pixels**,
//! `sqrt(lambda z) / pixel`. Geometrical optics replaces a band of width
//! `sqrt(lambda z)` with a step; whether that matters depends on how wide
//! the band is compared with a pixel.
//!
//! Calibrating it against the *exact* half-plane solution
//! ([`knife_edge`], two Fresnel integrals, affordable at any distance
//! where Rayleigh-Sommerfeld is not) gives a result worth stating
//! plainly: **a sharp shadow is a far better approximation than it
//! looks.** Pointwise at the edge it is wrong by 3/4, always. But averaged
//! over a pixel the step error cancels against the symmetric part of the
//! fringe pattern, and what survives is only the fringe *asymmetry* — so
//! the pixel error is 0.03 at a quarter-pixel blur and 1e-4 at a
//! thousandth. The bound is `min(0.25 blur, 3/4)`: the ceiling is derived
//! (the pointwise error at the edge is exactly `1 - 1/4`), the slope is
//! fitted, and [`shadow_misplaced_irradiance`] records why the obvious
//! derivation of the slope fails.
//!
//! None of the wave physics is reimplemented. `waveoptics::diffraction`
//! supplies the propagators and their bounds, `waveoptics::source` the
//! coherence of a real spectrum, `waveoptics::detector` the photon
//! statistics of a real sensor, and `photonics::multilayer` the complex
//! amplitude reflectance of a layered surface. What is here is the
//! geometric rung, the indicator that decides between them, and the
//! bridges into the renderer.

use emergence::{Assessment, Claim, FnRung, Ladder, Request, Rigor};
use num_complex::Complex64 as Cx;
use waveoptics::diffraction::{bounds, fraunhofer, fresnel, rayleigh_sommerfeld, Incident, Query, Rect};

/// Fresnel number of an aperture seen from a point: `a^2 / (lambda z)`,
/// with `a` the half-width of the smaller opening dimension.
///
/// Large means the aperture is many Fresnel zones across and the shadow is
/// sharp — geometrical optics applies. Of order one means a single zone
/// fills the opening and the pattern is all diffraction.
pub fn fresnel_number(aperture: &Rect, point: [f64; 3], wavelength: f64) -> f64 {
    let a = 0.5 * aperture.size[0].min(aperture.size[1]);
    let z = point[2].abs().max(1e-30);
    a * a / (wavelength * z)
}

/// The aperture's angular half-width of the first diffraction lobe,
/// `lambda / a`: the blur geometrical optics omits.
pub fn diffraction_angle(aperture: &Rect, wavelength: f64) -> f64 {
    let a = aperture.size[0].min(aperture.size[1]);
    wavelength / a.max(1e-30)
}

/// The Fresnel scale `sqrt(lambda z)`: the width of the penumbra that
/// geometrical optics replaces with a step.
///
/// This, not the Fresnel number, is the indicator that matters to a
/// renderer. A sharp shadow is wrong over a band of this width; whether
/// that matters depends on how wide the band is compared with a pixel.
/// At 550 nm and 10 mm it is 74 um; at 550 nm and 10 m it is 2.3 mm.
pub fn fresnel_scale(distance: f64, wavelength: f64) -> f64 {
    (wavelength * distance.abs()).sqrt()
}

/// The renderer's indicator: the penumbra width in pixels.
///
/// Below one, geometrical optics is right at the image's own resolution
/// and no wave computation can change a pixel. Above one, the shadow edge
/// is visibly soft and fringed, and rays are wrong in a way the image
/// shows.
pub fn shadow_blur_pixels(distance: f64, wavelength: f64, pixel_footprint: f64) -> f64 {
    fresnel_scale(distance, wavelength) / pixel_footprint.max(1e-300)
}

/// Exact scalar irradiance behind a straight edge, normalised to one in
/// the open, as a function of the Fresnel coordinate
/// `v = x sqrt(2 / (lambda z))` (positive into the lit region).
///
/// ```text
///   U(v) = integral_{-v}^{inf} exp(i pi t^2 / 2) dt = (1+i)/2 + F(v)
///   I(v) = |U(v)|^2 / |1 + i|^2 = |(1+i)/2 + F(v)|^2 / 2
/// ```
///
/// with `F` the Fresnel integral, which is odd. Deep in the light
/// `F -> (1+i)/2` and `I -> 1`; deep in the shadow `F -> -(1+i)/2` and
/// `I -> 0`.
///
/// This is the exact solution of the parabolic wave equation for a
/// half-plane (Born & Wolf, *Principles of Optics*, 7th ed., §8.7.1) and
/// the canonical answer to "how wrong is a sharp shadow". `I(0) = 1/4` at
/// the geometric edge — a quarter of the open irradiance, not a half — and
/// the first maximum is 1.37 at `v = 1.217`. Both are checked in the
/// tests.
///
/// It costs two Fresnel integrals, so unlike Rayleigh-Sommerfeld it is
/// affordable at any distance, which is what makes the calibration below
/// possible.
pub fn knife_edge(v: f64) -> f64 {
    (Cx::new(0.5, 0.5) + waveoptics::fresnel::fresnel(v)).norm_sqr() / 2.0
}

/// Geometrical optics' answer for the same edge: one in the light, zero in
/// the shadow, nothing in between.
pub fn knife_edge_geometric(v: f64) -> f64 {
    if v > 0.0 {
        1.0
    } else {
        0.0
    }
}

/// Slope of the renderer's error claim: the pixel-averaged irradiance
/// error of a sharp shadow is bounded by `SHADOW_ERROR_SLOPE * blur` for
/// blur up to a pixel.
///
/// **Fitted**, against the exact [`knife_edge`] over four decades of blur
/// (`tests/wave.rs`): the measured worst ratio below one pixel is 0.16, so
/// this carries about 50 % headroom. See
/// [`shadow_misplaced_irradiance`] for why the obvious derivation of this
/// slope does not work.
pub const SHADOW_ERROR_SLOPE: f64 = 0.25;

/// Where the claim saturates, and why it is exactly three quarters.
///
/// However wide the penumbra, the *worst* a pixel can be wrong is set by
/// the pointwise error at the geometric edge: geometrical optics says the
/// irradiance there is 1, and the exact half-plane solution says
/// [`knife_edge(0.0)`](knife_edge) `= 1/4`. So the error can never exceed
/// `1 - 1/4`. That is a derived bound, not a fit.
pub const SHADOW_ERROR_CEILING: f64 = 0.75;

/// The renderer's error claim for geometrical optics at a shadow edge:
/// `min(SHADOW_ERROR_SLOPE * blur, SHADOW_ERROR_CEILING)`.
///
/// This is the number that decides whether a wave computation can change
/// the image. It is an asymptotic estimate below the ceiling and a bound
/// at it, and [`shadow_claim`] labels it accordingly.
pub fn shadow_error_claim(blur_pixels: f64) -> Claim {
    let linear = SHADOW_ERROR_SLOPE * blur_pixels.abs();
    if linear >= SHADOW_ERROR_CEILING {
        Claim::rigorous(SHADOW_ERROR_CEILING)
    } else {
        Claim::asymptotic(linear)
    }
}

/// [`shadow_error_claim`] from the physical quantities.
pub fn shadow_claim(distance: f64, wavelength: f64, pixel_footprint: f64) -> Claim {
    shadow_error_claim(shadow_blur_pixels(distance, wavelength, pixel_footprint))
}

/// `|integral_{-reach}^{reach} (I_exact - I_geometric) dv|` in Fresnel
/// units: the net irradiance a sharp shadow misplaces at an edge.
///
/// **It tends to zero**, halving as the reach doubles. Diffraction moves
/// irradiance across the edge — a 37 % overshoot on the lit side, a tail
/// into the shadow — but creates none, so the signed error integrates
/// away. The tests check the decay.
///
/// That is also a negative result worth keeping. The obvious way to derive
/// the pixel error — "a pixel of width `dv` containing the edge misplaces
/// a fixed amount `K`, so its mean error is `K/dv`" — fails exactly
/// because `K = 0`. The pixel error therefore falls *faster* than `1/dv`
/// for a pixel centred on the edge, and the worst case is an offset pixel
/// that catches the cumulative error near its peak. The slope in
/// [`SHADOW_ERROR_SLOPE`] is fitted for that reason, not derived.
pub fn shadow_misplaced_irradiance(reach: f64, step: f64) -> f64 {
    let n = ((2.0 * reach / step).ceil() as usize).max(16);
    let h = 2.0 * reach / n as f64;
    let mut acc = 0.0;
    for i in 0..n {
        let v = -reach + h * (i as f64 + 0.5);
        acc += (knife_edge(v) - knife_edge_geometric(v)) * h;
    }
    acc.abs()
}

/// Mean absolute irradiance error of a sharp shadow over a pixel of width
/// `pixel` whose centre sits `offset` from the geometric edge.
///
/// The average is taken in Fresnel coordinates, where the integrand
/// `I_exact - I_geometric` is smooth and decays away from the edge. That
/// makes the cost independent of how small the blur is — integrating in
/// `x` instead needs a sample per fringe across the whole pixel, which is
/// unbounded as the blur shrinks and is how this was first written.
pub fn shadow_error_over_pixel(offset: f64, pixel: f64, distance: f64, wavelength: f64) -> f64 {
    let scale = (2.0 / (wavelength * distance.abs())).sqrt();
    let dv = pixel * scale;
    if !(dv > 0.0) {
        return 0.0;
    }
    let v0 = offset * scale;
    // The profiles agree to 1/v, so the cumulative difference keeps
    // oscillating with a slowly shrinking envelope; 200 Fresnel units is
    // where truncating stops moving the answer at the per-cent level.
    const REACH: f64 = 200.0;
    let (lo, hi) = ((v0 - 0.5 * dv).max(-REACH), (v0 + 0.5 * dv).min(REACH));
    if hi <= lo {
        return 0.0;
    }
    let n = (((hi - lo) / 0.05).ceil() as usize).clamp(16, 20_000);
    let h = (hi - lo) / n as f64;
    let mut acc = 0.0;
    for i in 0..n {
        let v = lo + h * (i as f64 + 0.5);
        acc += (knife_edge(v) - knife_edge_geometric(v)) * h;
    }
    acc.abs() / dv
}

/// The worst pixel-averaged shadow error anywhere near the edge: the
/// quantity the renderer's claim has to bound.
pub fn worst_shadow_error(pixel: f64, distance: f64, wavelength: f64) -> f64 {
    let blur = fresnel_scale(distance, wavelength);
    let reach = 4.0 * blur.max(pixel);
    let mut worst: f64 = 0.0;
    for i in 0..=120 {
        let offset = -reach + 2.0 * reach * i as f64 / 120.0;
        worst = worst.max(shadow_error_over_pixel(offset, pixel, distance, wavelength));
    }
    worst
}

/// Geometrical optics: the irradiance of the aperture's shadow.
///
/// Inside the projected opening the incident irradiance passes; outside it
/// is zero. There is no edge, no fringe and no spreading — that is exactly
/// the content of the approximation, and what the penumbra width measures
/// the cost of.
pub fn geometric_irradiance(qr: &Query) -> f64 {
    let p = qr.point;
    let (shift, amplitude) = match qr.incident {
        Incident::Plane { amplitude, sines } => {
            ([sines[0] * p[2], sines[1] * p[2]], amplitude.norm_sqr())
        }
        Incident::Point {
            amplitude,
            position,
            distance,
        } => {
            // A point source magnifies the opening by (z_s + z)/z_s.
            let m = (distance + p[2]) / distance;
            let scale = amplitude.norm_sqr() / (distance + p[2]).powi(2);
            (
                [
                    position[0] * (1.0 - m) + qr.aperture.center[0] * (m - 1.0),
                    position[1] * (1.0 - m) + qr.aperture.center[1] * (m - 1.0),
                ],
                scale,
            )
        }
    };
    let inside = (0..2).all(|a| {
        let centre = qr.aperture.center[a] + shift[a];
        (p[a] - centre).abs() <= 0.5 * qr.aperture.size[a]
    });
    if inside {
        amplitude
    } else {
        0.0
    }
}

/// Irradiance from a complex field amplitude.
/// Irradiance from a complex field amplitude.
pub fn irradiance(field: Cx) -> f64 {
    field.norm_sqr()
}

/// The optical ladder for the irradiance behind a rectangular aperture.
///
/// Rung 0 is geometrical optics; rungs 1 and 2 are Fraunhofer and Fresnel
/// with the rigorous field bounds `waveoptics` proves for them; rung 3 is
/// Rayleigh–Sommerfeld, the reference. Costs are in units of one field
/// evaluation.
///
/// `refine` sets the reference's quadrature refinement.
pub fn irradiance_ladder(refine: usize) -> Ladder<Query, f64> {
    // Relative irradiance error from a bound on the field amplitude:
    // |I - I0|/I0 <= 2 db/|A| + (db/|A|)^2, with |A| the reference
    // amplitude. The bound is on the field, so it is converted here.
    fn field_bound_to_irradiance(db: f64, amplitude: f64) -> f64 {
        if amplitude <= 0.0 {
            return f64::INFINITY;
        }
        let r = db / amplitude;
        2.0 * r + r * r
    }

    let geometric = FnRung::new(
        "geometric",
        |_: &Query| 1.0,
        |qr: &Query| {
            // The geometric rung's error is set by the penumbra, which it
            // omits entirely: over a region the size of the aperture the
            // mean error is the Fresnel scale divided by the aperture.
            // This is the same statement as the pixel indicator, with the
            // aperture playing the part of the pixel.
            let a = qr.aperture.size[0].min(qr.aperture.size[1]);
            if !(a > 0.0) || !(qr.point[2].abs() > 0.0) {
                return Assessment::OutOfDomain;
            }
            let blur = fresnel_scale(qr.point[2], qr.wavelength) / a;
            Assessment::Claim(shadow_error_claim(blur))
        },
        geometric_irradiance,
    );

    let fraunhofer_rung = FnRung::new(
        "fraunhofer",
        |_: &Query| 2.0,
        |qr: &Query| {
            let (_, fh) = bounds(qr);
            let a = fraunhofer(qr).norm();
            Assessment::Claim(Claim::new(
                field_bound_to_irradiance(fh, a),
                Rigor::Rigorous,
            ))
        },
        |qr: &Query| irradiance(fraunhofer(qr)),
    );

    let fresnel_rung = FnRung::new(
        "fresnel",
        |_: &Query| 6.0,
        |qr: &Query| {
            let (fr, _) = bounds(qr);
            let a = fresnel(qr).norm();
            Assessment::Claim(Claim::new(
                field_bound_to_irradiance(fr, a),
                Rigor::Rigorous,
            ))
        },
        |qr: &Query| irradiance(fresnel(qr)),
    );

    let reference = FnRung::new(
        "rayleigh-sommerfeld",
        move |_: &Query| 400.0 * (refine.max(1) as f64).powi(2),
        |_: &Query| Assessment::Claim(Claim::validated(1e-10)),
        move |qr: &Query| irradiance(rayleigh_sommerfeld(qr, refine)),
    );

    let mut ladder = Ladder::new("aperture irradiance", |a: &f64, b: &f64| {
        let scale = a.abs().max(b.abs()).max(1e-300);
        (a - b).abs() / scale
    });
    ladder.push(geometric);
    ladder.push(fraunhofer_rung);
    ladder.push(fresnel_rung);
    let r = ladder.push(reference);
    ladder.reference(r)
}

/// Solve for the irradiance at a tolerance, letting the ladder choose.
pub fn irradiance_at(qr: &Query, tolerance: f64, refine: usize) -> emergence::Solved<f64> {
    irradiance_ladder(refine).solve(qr, Request::new(tolerance))
}

/// Analytic single-slit irradiance in the far field, normalised to one on
/// axis: `sinc^2(pi a sin(theta) / lambda)`.
///
/// Its first zero is at `sin(theta) = lambda / a` exactly, which is what
/// the tests check the renderer's diffraction against.
pub fn single_slit(sin_theta: f64, width: f64, wavelength: f64) -> f64 {
    let x = std::f64::consts::PI * width * sin_theta / wavelength;
    if x.abs() < 1e-12 {
        1.0
    } else {
        (x.sin() / x).powi(2)
    }
}

/// Analytic two-slit irradiance: the single-slit envelope times
/// `cos^2(pi d sin(theta)/lambda)`, so the fringes are spaced
/// `lambda / d` in `sin(theta)`.
pub fn double_slit(sin_theta: f64, width: f64, separation: f64, wavelength: f64) -> f64 {
    let c = (std::f64::consts::PI * separation * sin_theta / wavelength).cos();
    single_slit(sin_theta, width, wavelength) * c * c
}

/// First zero of a circular aperture's Airy pattern:
/// `sin(theta) = 1.2196699 lambda / D`, the first zero of `2 J1(x)/x`
/// divided by pi (Born & Wolf, *Principles of Optics*, 7th ed., §8.5.2).
pub fn airy_first_zero(diameter: f64, wavelength: f64) -> f64 {
    1.219_669_891 * wavelength / diameter
}

/// A layered surface whose reflectance is wave-optical: a coating, an
/// oxide, a soap film, an anti-reflection stack.
///
/// The reflectance comes from `photonics::multilayer`, so it is the full
/// transfer-matrix result — every internal reflection summed with its
/// phase, both polarisations, complex indices allowed. A ray tracer's
/// Fresnel term is the special case of a single interface.
#[derive(Clone, Debug)]
pub struct Film {
    pub stack: photonics::multilayer::Stack,
}

impl Film {
    /// A single non-absorbing layer of index `n` and thickness `d` on a
    /// substrate of index `n_sub`, in air.
    pub fn single(n: f64, thickness: f64, n_sub: f64) -> Self {
        Self {
            stack: photonics::multilayer::Stack::new(
                1.0,
                vec![photonics::multilayer::Layer::new(n, 0.0, thickness)],
                n_sub,
            ),
        }
    }

    pub fn from_stack(stack: photonics::multilayer::Stack) -> Self {
        Self { stack }
    }

    /// Power reflectance at a wavelength, incidence and polarisation.
    pub fn reflectance(
        &self,
        wavelength: f64,
        theta: f64,
        pol: photonics::multilayer::Polarization,
    ) -> f64 {
        self.stack.solve(wavelength, theta, pol).reflectance
    }

    /// Complex amplitude reflection coefficient: the phase is what makes
    /// interference between paths possible.
    pub fn amplitude(
        &self,
        wavelength: f64,
        theta: f64,
        pol: photonics::multilayer::Polarization,
    ) -> Cx {
        let r = self.stack.solve(wavelength, theta, pol).r;
        Cx::new(r.re, r.im)
    }

    /// Unpolarised reflectance: the mean of the two polarisations, which is
    /// what an unpolarised ray tracer should use.
    pub fn unpolarised(&self, wavelength: f64, theta: f64) -> f64 {
        use photonics::multilayer::Polarization::{TE, TM};
        0.5 * (self.reflectance(wavelength, theta, TE) + self.reflectance(wavelength, theta, TM))
    }

    /// The colour of the film in reflection under an equal-energy
    /// illuminant, as linear RGB.
    ///
    /// This is the whole point of a wave-optical surface: a 300 nm soap
    /// film is not grey, and no reflectance scalar can make it the right
    /// colour. The spectrum is integrated against the CIE observer.
    pub fn reflected_rgb(&self, theta: f64, samples: usize) -> [f64; 3] {
        let (lo, hi) = (380e-9, 780e-9);
        let n = samples.max(8);
        let mut xyz = [0.0; 3];
        let mut weight = 0.0;
        for i in 0..n {
            let lambda = lo + (hi - lo) * (i as f64 + 0.5) / n as f64;
            let r = self.unpolarised(lambda, theta);
            let bar = radiation::optics::cie_xyz(lambda);
            for k in 0..3 {
                xyz[k] += bar[k] * r;
            }
            weight += bar[1];
        }
        if weight > 0.0 {
            for v in &mut xyz {
                *v /= weight;
            }
        }
        radiation::optics::xyz_to_rgb(xyz)
    }
}

/// Reflectance of a single non-absorbing film from the two-beam (Airy)
/// formula, the closed form the transfer matrix must reproduce:
///
/// ```text
///   R = (r1^2 + r2^2 + 2 r1 r2 cos d) / (1 + r1^2 r2^2 + 2 r1 r2 cos d)
///   d = 4 pi n1 d cos(theta1) / lambda
/// ```
///
/// with `r1`, `r2` the amplitude Fresnel coefficients of the two
/// interfaces at normal incidence (Born & Wolf §1.6.4).
pub fn airy_reflectance(n0: f64, n1: f64, n2: f64, thickness: f64, wavelength: f64) -> f64 {
    let r1 = (n0 - n1) / (n0 + n1);
    let r2 = (n1 - n2) / (n1 + n2);
    let delta = 4.0 * std::f64::consts::PI * n1 * thickness / wavelength;
    let num = r1 * r1 + r2 * r2 + 2.0 * r1 * r2 * delta.cos();
    let den = 1.0 + r1 * r1 * r2 * r2 + 2.0 * r1 * r2 * delta.cos();
    num / den
}

/// Brewster angle: the incidence at which p-polarised light is not
/// reflected at all, `atan(n2/n1)`. An exact zero, and one of the sharpest
/// checks available on a Fresnel implementation.
pub fn brewster_angle(n1: f64, n2: f64) -> f64 {
    (n2 / n1).atan()
}

/// Fringe visibility of two interfering paths differing by
/// `path_difference`, for a source of the given spectral width.
///
/// Partial coherence is why a renderer may add two paths' intensities
/// rather than their amplitudes: beyond the coherence length the cross
/// term averages away and the wave answer becomes the ray answer. This is
/// `waveoptics::source`'s Gaussian visibility, which is the exact Fourier
/// transform of a Gaussian spectrum.
pub fn visibility(spectrum: &waveoptics::source::Spectrum, path_difference: f64) -> f64 {
    spectrum.gaussian_visibility(path_difference)
}

/// Coherence length of a Gaussian spectrum of centre `lambda` and full
/// width at half maximum `fwhm`: `lambda^2 / fwhm` to within the usual
/// factor, defined here as the path difference at which the visibility
/// falls to `1/e`.
pub fn coherence_length(spectrum: &waveoptics::source::Spectrum) -> f64 {
    // Bisect the visibility for 1/e.
    let target = std::f64::consts::E.recip();
    let (mut lo, mut hi) = (0.0, 1.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if spectrum.gaussian_visibility(mid) > target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Two interfering paths, combined with the coherence their source allows.
///
/// With full coherence this is `|a1 + a2|^2`; with none it is
/// `|a1|^2 + |a2|^2`; in between the cross term carries the visibility.
/// Taking the ray-optics answer is the `visibility = 0` limit, which is why
/// a renderer that only ever adds intensities is a *partially coherent*
/// renderer with the visibility silently set to zero.
pub fn combine_paths(a1: Cx, a2: Cx, visibility: f64) -> f64 {
    let i1 = a1.norm_sqr();
    let i2 = a2.norm_sqr();
    let cross = 2.0 * (a1 * a2.conj()).re;
    i1 + i2 + visibility.clamp(0.0, 1.0) * cross
}

/// Photon statistics: turn an expected electron count into a measured one.
///
/// A render is a measurement, and a measurement of light is a count. At
/// high exposure the Poisson noise is invisible and the image is the
/// radiance; at low exposure the image *is* the noise, and no amount of
/// path tracing removes it, because it is not sampling error. This uses
/// `waveoptics::detector`'s generator.
pub fn photon_count(expected: f64, rng: &mut waveoptics::detector::Rng) -> f64 {
    rng.poisson(expected.max(0.0))
}

/// Expected electrons per pixel from a radiance, aperture solid angle,
/// pixel area, exposure and quantum efficiency.
///
/// `radiance` is spectral, W/(m^2 sr m); the result is a count.
pub fn expected_electrons(
    radiance: f64,
    wavelength: f64,
    bandwidth: f64,
    solid_angle: f64,
    pixel_area: f64,
    exposure: f64,
    qe: f64,
) -> f64 {
    let photon_energy = radiation::H * radiation::C / wavelength;
    radiance * bandwidth * solid_angle * pixel_area * exposure * qe / photon_energy
}
