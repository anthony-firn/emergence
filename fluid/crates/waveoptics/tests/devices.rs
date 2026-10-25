//! Rendering devices against exact optics: a camera focused at infinity
//! looking through a double slit at a distant source records the slits'
//! angular spectrum (Fraunhofer pattern) whatever the slit distance, and
//! its fringe visibility is the slits' mutual coherence.

use std::f64::consts::PI;
use std::time::Instant;
use waveoptics::devices::{Device, Optics};
use waveoptics::experiment::Experiment;
use waveoptics::imaging::Window;
use waveoptics::polarization::Jones;
use waveoptics::source::{Extent, Source, Spectrum};

const HENE: f64 = 632.8e-9;
const A: f64 = 100e-6; // slit width
const D: f64 = 500e-6; // separation
const H: f64 = 300e-6; // slit height
const F: f64 = 50e-3; // focal length
const Z: f64 = 20e-3; // slits to lens

fn sinc(u: f64) -> f64 {
    if u == 0.0 {
        1.0
    } else {
        u.sin() / u
    }
}

fn camera() -> Device {
    Device { optics: Optics::Ideal { radius: 1.5e-3, focal: F }, front: [0.0, 0.0, Z] }
}

#[test]
fn ideal_camera_records_the_exact_fraunhofer_pattern() {
    let e0 = 1.0;
    let src = Source { spectrum: Spectrum::Line(HENE), extent: Extent::Collimated { irradiance: e0 }, polarization: Some(Jones::horizontal()) };
    let exp = Experiment::double_slit(src, A, D, H, 1.0, 1e-4);
    let window = Window { center: [0.0, 0.0], pixels: [96, 32], pitch: 4e-6 };
    let t = Instant::now();
    let img = camera().spectral_image(&exp, &window);
    let irr = img.irradiance();
    // Angular spectrum of the Kirchhoff boundary values, focused at
    // sin θ = x'/f: I = E₀ |a h sinc sinc 2cos|²/(λ f)².
    let (xs, ys) = (window.coords(0), window.coords(1));
    let model = |x: f64, y: f64| {
        let (fx, fy) = (x / (HENE * F), y / (HENE * F));
        let u = A * H * sinc(PI * A * fx) * sinc(PI * H * fy) * 2.0 * (PI * D * fx).cos();
        e0 * (u / (HENE * F)).powi(2)
    };
    let i0 = model(0.0, 0.0);
    let mut worst: f64 = 0.0;
    for (q, v) in irr.iter().enumerate() {
        worst = worst.max((v - model(xs[q % xs.len()], ys[q / xs.len()])).abs() / i0);
    }
    println!("band loss {:.2e}, pupil {}², {} mode(s), {:.2} s; peak {:.4e} W/m² (exact {:.4e}); max |I − exact|/I₀ = {worst:.2e}", img.band_loss, img.pupil_samples, img.modes, t.elapsed().as_secs_f64(), irr.iter().cloned().fold(0.0, f64::max), i0);
    // Only directions near the window are propagated (band loss: light
    // that lands elsewhere on the detector); what remains differs from the
    // exact pattern by the pupil's coherent spread of the band edge.
    assert!(worst < 1e-3);
    // Fringe period on the detector, from zero crossings of the
    // interference term I − I_envelope (exact: f λ/d).
    let row = &irr[xs.len() * (ys.len() / 2)..xs.len() * (ys.len() / 2 + 1)];
    let y = ys[ys.len() / 2];
    let g: Vec<f64> = xs.iter().zip(row).map(|(&x, v)| v - 0.5 * model(x, y) / (PI * D * x / (HENE * F)).cos().powi(2).max(1e-300) * 1.0).collect();
    let zeros: Vec<f64> = (0..g.len() - 1).filter(|&i| g[i] * g[i + 1] < 0.0).map(|i| xs[i] - g[i] * (xs[i + 1] - xs[i]) / (g[i + 1] - g[i])).collect();
    let period = 2.0 * (zeros.last().unwrap() - zeros[0]) / (zeros.len() - 1) as f64;
    println!("fringe period {:.3} µm from {} zeros (f λ/d = {:.3} µm)", period * 1e6, zeros.len(), F * HENE / D * 1e6);
    assert!((period - F * HENE / D).abs() / (F * HENE / D) < 1e-3);
}

/// Fringe visibility near the image centre: fit `a + b cos + c sin` at the
/// known period to the image divided by the single-slit envelope.
fn visibility(xs: &[f64], row: &[f64]) -> f64 {
    let k = 2.0 * PI * D / (HENE * F);
    let mut m = [[0.0; 3]; 3];
    let mut r = [0.0; 3];
    for (&x, &v) in xs.iter().zip(row) {
        let env = sinc(PI * A * x / (HENE * F)).powi(2);
        let phi = [1.0, (k * x).cos(), (k * x).sin()];
        for i in 0..3 {
            r[i] += phi[i] * v / env;
            for j in 0..3 {
                m[i][j] += phi[i] * phi[j];
            }
        }
    }
    // 3×3 solve (Cramer)
    let det = |m: [[f64; 3]; 3]| m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    let d0 = det(m);
    let c: Vec<f64> = (0..3)
        .map(|col| {
            let mut mm = m;
            for row in 0..3 {
                mm[row][col] = r[row];
            }
            det(mm) / d0
        })
        .collect();
    c[1].hypot(c[2]) / c[0]
}

#[test]
fn camera_fringe_visibility_is_the_slits_mutual_coherence() {
    let zs = 1.0;
    let window = Window { center: [0.0, 0.0], pixels: [64, 1], pitch: 4e-6 };
    for &s in &[0.4e-3, 0.9e-3] {
        let src = Source { spectrum: Spectrum::Line(HENE), extent: Extent::Slit { power: 1e-3, distance: zs, width: s }, polarization: Some(Jones::horizontal()) };
        let exp = Experiment::double_slit(src, A, D, H, 1.0, 1e-4);
        let gamma = exp.coherence(&[[0.0, 0.0]])[0].abs();
        let t = Instant::now();
        let img = camera().spectral_image(&exp, &window);
        let v = visibility(&window.coords(0), &img.irradiance());
        let vcz = sinc(PI * s * D / (HENE * zs)).abs();
        println!("source {:.1} mm: camera visibility {v:.4}, slit-plane |γ| {gamma:.4} (VCZ {vcz:.4}); {} modes, pupil {}², {:.2} s", s * 1e3, img.modes, img.pupil_samples, t.elapsed().as_secs_f64());
        assert!((v - gamma).abs() < 5e-3);
    }
}

use waveoptics::angular::{plane_fields, Band, OpeningSpectrum};
use waveoptics::diffraction::{ladder, Incident, Query, Rect};
use waveoptics::experiment::field_scale;
use num_complex::Complex64 as Cx;

#[test]
fn angular_spectrum_fields_match_the_certified_ladder() {
    // One 100 × 300 µm opening, 20 mm in front of the evaluation plane
    // (Fresnel number ≈ 0.8 across, 7 along): a tilted plane wave and a
    // point source 0.5 m away, off axis. Full band to sines of 0.1 (the
    // points lie within 0.06 of the opening's direction).
    let rect = Rect { center: [0.1e-3, 0.0], size: [A, H] };
    let z = Z;
    let lad = ladder(0);
    let pts = [[0.0, 0.0], [0.15e-3, 0.05e-3], [0.4e-3, -0.1e-3], [-0.6e-3, 0.2e-3], [1.0e-3, 0.0]];
    let (xs, ys): (Vec<f64>, Vec<f64>) = (pts.iter().map(|p| p[0]).collect(), pts.iter().map(|p| p[1]).collect());
    for incident in [
        Incident::Plane { amplitude: Cx::new(1.0, 0.0), sines: [0.01, -0.005] },
        Incident::Point { amplitude: Cx::new(0.5, 0.0), position: [2e-3, -1e-3], distance: 0.5 },
    ] {
        let sines = match incident {
            Incident::Plane { sines, .. } => sines,
            Incident::Point { position, distance, .. } => [(rect.center[0] - position[0]) / distance, (rect.center[1] - position[1]) / distance],
        };
        let spec = OpeningSpectrum::new(&rect, &incident, HENE);
        let mut worst = [0.0f64; 2];
        for (b, inner) in [0.10, 0.14].iter().enumerate() {
            let band = Band { center: [sines[0] / HENE, sines[1] / HENE], inner: [inner / HENE; 2], outer: [(inner + 0.04) / HENE; 2] };
            let period = [2.0 * z * (inner + 0.04) + 4.0 * HENE / 0.04 + 4e-3; 2];
            let t = Instant::now();
            let (fields, loss) = plane_fields(&[spec], HENE, z, &band, period, &xs, &ys);
            for (i, p) in pts.iter().enumerate() {
                let q = Query { aperture: rect, incident, point: [p[0], p[1], z], wavelength: HENE };
                let exact = lad.solve(&q, emergence::Request::new(1e-5 * field_scale(&q))).value.u;
                let got = fields[0][i + xs.len() * i];
                worst[b] = worst[b].max((got - exact).norm() / field_scale(&q));
            }
            println!("{incident:?}: band ±{inner}: max |U − ladder|/scale {:.2e}, power outside band {:.2e}, incident model bound {:.1e} rad, {:.2} s", worst[b], loss[0], spec.phase_bound, t.elapsed().as_secs_f64());
        }
        // the ladder is certified to 1e-5 of the field scale
        assert!(worst[0] < 2e-5 && worst[1] < 2e-5);
    }
}

use waveoptics::detector::Rng;
use waveoptics::imaging::isoplanatic_image;

#[test]
fn rough_screen_gives_fully_developed_speckle_and_edges_ring_as_coherent_theory() {
    // Random-phase screen (unit amplitude, phase uniform on [0, 2π) per
    // sample) seen by an f/8, 50 mm lens 0.5 m away: intensity statistics
    // are negative-exponential (contrast 1) and the mean is the pupil's
    // share of the white angular spectrum, π f_c² δ² with f_c = a/(λD).
    let (nx, ny, delta, lam, dist, a) = (512, 256, 20e-6, HENE, 0.5, 3.125e-3);
    let mut rng = Rng(7);
    let obj: Vec<Cx> = (0..nx * ny).map(|_| Cx::from_polar(1.0, 2.0 * PI * rng.uniform())).collect();
    let pupil = |u: [f64; 2]| if u[0] * u[0] + u[1] * u[1] <= a * a { Cx::new(1.0, 0.0) } else { Cx::new(0.0, 0.0) };
    let img = isoplanatic_image(&obj, nx, ny, delta, lam, dist, pupil);
    // interior only (the padded FFT sees zero field outside)
    let inner: Vec<f64> = (0..nx * ny).filter(|q| (q % nx).abs_diff(nx / 2) < nx / 3 && (q / nx).abs_diff(ny / 2) < ny / 3).map(|q| img[q]).collect();
    let mean = inner.iter().sum::<f64>() / inner.len() as f64;
    let sd = (inner.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / inner.len() as f64).sqrt();
    let fc = a / (lam * dist);
    let expect = PI * fc * fc * delta * delta;
    // speckle cells in the sample: area/(λD/2a)² → relative standard error
    let cells = inner.len() as f64 * delta * delta / (lam * dist / (2.0 * a)).powi(2);
    println!("speckle: contrast {:.4}, mean {:.4e} (π f_c² δ² = {:.4e}), ~{:.0} speckles", sd / mean, mean, expect, cells);
    let se = 3.0 / cells.sqrt();
    assert!((sd / mean - 1.0).abs() < se && (mean / expect - 1.0).abs() < se);
    // A coherent straight edge through a circular pupil: with a smooth
    // profile along the edge (Gaussian, σ = 0.15 mm, spectrum ≪ f_c) only
    // f_y ≈ 0 matters, so the amplitude is 1/2 + Si(2π f_c x)/π and the
    // intensity overshoots to (1/2 + Si(π)/π)² = 1.1873 (Gibbs).
    let (nx, ny, delta) = (1024, 1024, 2e-6);
    let edge: Vec<Cx> = (0..nx * ny)
        .map(|q| {
            let y = (q / nx) as f64 * delta - 0.5 * ny as f64 * delta;
            let x = (q % nx) as f64 * delta - 0.5 * nx as f64 * delta;
            // far end rolled off smoothly (its own ringing would reach the
            // first maximum at ~0.5% in amplitude)
            let tail = if x < 0.4e-3 { 1.0 } else if x < 0.9e-3 { (0.5 * PI * (x - 0.4e-3) / 0.5e-3).cos().powi(2) } else { 0.0 };
            if x >= 0.0 { Cx::new(tail * (-0.5 * (y / 0.15e-3).powi(2)).exp(), 0.0) } else { Cx::new(0.0, 0.0) }
        })
        .collect();
    let img = isoplanatic_image(&edge, nx, ny, delta, lam, dist, pupil);
    let row = &img[nx * (ny / 2)..nx * (ny / 2 + 1)];
    let i = (nx / 2..nx / 2 + 200).max_by(|&a, &b| row[a].total_cmp(&row[b])).unwrap();
    // parabolic refinement of the sampled maximum
    let (l, c, r) = (row[i - 1], row[i], row[i + 1]);
    let peak = c + 0.125 * (r - l).powi(2) / (2.0 * c - l - r);
    println!("coherent edge overshoot {peak:.4} (Gibbs 1.1873)");
    assert!((peak - 1.1873).abs() < 5e-3);
}

use waveoptics::devices::{Retina, SpectralImage};

#[test]
fn retina_catches_photons_as_cone_fundamentals_predict() {
    let window = Window { center: [0.3e-3, 0.0], pixels: [40, 40], pitch: 2e-6 };
    let retina = Retina::default();
    // Uniform 550 nm, 1 W/m² on the retina: each cone's mean catch is
    // E λ/(hc) A s_Q(λ) t with s_Q the quantal sensitivity (peak 1).
    let lam = 550e-9;
    let img = SpectralImage { window, bands: vec![(lam, vec![1.0; 1600])], pupil_samples: 0, modes: 1, band_loss: 0.0, phase_bound: 0.0 };
    let cones = retina.sample(&img, &mut Rng(3));
    let hc = 6.62607015e-34 * 299792458.0;
    let peak_q: [f64; 3] = std::array::from_fn(|k| (380..=780).map(|nm| waveoptics::color::lms(nm as f64 * 1e-9)[k] / nm as f64).fold(0.0, f64::max));
    for kind in 0..3 {
        let c: Vec<&_> = cones.iter().filter(|c| c.kind == kind).collect();
        let sq = waveoptics::color::lms(lam)[kind] / 550.0 / peak_q[kind];
        let expect = lam / hc * retina.peak_area * sq * retina.integration;
        let mean = c.iter().map(|c| c.isomerisations).sum::<f64>() / c.len() as f64;
        let var = c.iter().map(|c| (c.isomerisations - mean).powi(2)).sum::<f64>() / (c.len() - 1) as f64;
        println!("cone {kind}: {} cones, mean catch {mean:.4e} (expected {expect:.4e}), variance/mean {:.3}", c.len(), var / mean);
        assert!((c[0].expected - expect).abs() < 1e-9 * expect);
        let se = 4.0 * (expect / c.len() as f64).sqrt();
        assert!((mean - expect).abs() < se && (var / mean - 1.0).abs() < 4.0 * (2.0 / c.len() as f64).sqrt());
    }
    // A uniform broadband field: the mosaic's reconstruction returns the
    // field's CIE 2006 XYZ (photon noise negligible at this level).
    let bands: Vec<(f64, Vec<f64>)> = (400..=700).step_by(10).map(|nm| (nm as f64 * 1e-9, vec![0.05; 1600])).collect();
    let img = SpectralImage { window, bands: bands.clone(), pupil_samples: 0, modes: 1, band_loss: 0.0, phase_bound: 0.0 };
    let cones = retina.sample(&img, &mut Rng(5));
    let xyz = retina.perceived_xyz(&cones, &window);
    let truth: [f64; 3] = std::array::from_fn(|k| bands.iter().map(|(l, e)| waveoptics::color::xyz(*l)[k] * e[0]).sum());
    let centre = xyz[20 + 40 * 20];
    println!("perceived XYZ {centre:.4?} vs CIE {truth:.4?}");
    for k in 0..3 {
        assert!((centre[k] - truth[k]).abs() < 1e-2 * truth[1]);
    }
}
