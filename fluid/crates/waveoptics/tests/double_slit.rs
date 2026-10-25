//! Young's double slit against closed-form optics: absolute irradiance,
//! fringe spacing and envelope, spatial coherence (van Cittert–Zernike),
//! temporal coherence, and which-path marking / erasure with polarizers.

use num_complex::Complex64 as Cx;
use std::f64::consts::PI;
use waveoptics::experiment::Experiment;
use waveoptics::fresnel::fresnel_span;
use waveoptics::polarization::{Jones, JonesMatrix};
use waveoptics::source::{Extent, Source, Spectrum};

const HENE: f64 = 632.8e-9;
const A: f64 = 50e-6; // slit width
const D: f64 = 250e-6; // slit separation (centre to centre)
const H: f64 = 3e-3; // slit height
const L: f64 = 1.0; // slits to screen

fn laser(irradiance: f64) -> Source {
    Source {
        spectrum: Spectrum::Line(HENE),
        extent: Extent::Collimated { irradiance },
        polarization: Some(Jones::horizontal()),
    }
}

#[test]
fn laser_double_slit_absolute_irradiance_spacing_and_envelope() {
    let e0 = 100.0; // W/m² after a beam expander (estimate of a 1 mW He–Ne spread over ~1 cm²)
    let exp = Experiment::double_slit(laser(e0), A, D, H, L, 1e-4);
    // Absolute centre irradiance: each slit contributes, in the Fresnel
    // regime along its height and the Fraunhofer regime across it,
    // U = (1/(iλL)) e^{ikL} a √(λL/2) [F(h/√(2λL)) − F(−h/√(2λL))].
    let w = H / (2.0 * HENE * L).sqrt();
    let u1 = A * (HENE * L / 2.0).sqrt() * fresnel_span(-w, w) / (HENE * L);
    let i_centre = e0 * (2.0 * u1.norm()).powi(2);
    let got = exp.irradiance(&[[0.0, 0.0]]);
    println!("centre irradiance {:.4} W/m² (closed form {:.4})", got.values[0], i_centre);
    assert!((got.values[0] - i_centre).abs() / i_centre < 2e-3);
    // Fringes along x: cos²(π d x/(λL)) sinc²(π a x/(λL)) relative to the centre.
    let xs: Vec<[f64; 2]> = (0..801).map(|i| [-0.02 + 0.04 * i as f64 / 800.0, 0.0]).collect();
    let irr = exp.irradiance(&xs);
    let mut worst: f64 = 0.0;
    for (p, v) in xs.iter().zip(&irr.values) {
        let x = p[0];
        let r = (x * x + L * L).sqrt();
        let sinc = |u: f64| if u == 0.0 { 1.0 } else { u.sin() / u };
        // far-field model incl. obliquity (L/r)² and 1/r² (cancel to ~1e-4 here)
        let model = i_centre * (PI * D * x / (HENE * r)).cos().powi(2) * sinc(PI * A * x / (HENE * r)).powi(2) * (L / r).powi(4);
        worst = worst.max((v - model).abs() / i_centre);
    }
    println!("max |I − cos²·sinc²|/I₀ over ±20 mm: {worst:.2e}; rungs {:?}", irr.rungs);
    assert!(worst < 3e-3);
    // Fringe spacing, measured the experimenter's way: block one slit, then
    // the other; the interference term (I − I₁ − I₂) has equally spaced
    // zeros (the maxima of I are pulled inward by the sinc² envelope).
    let pts: Vec<[f64; 2]> = (0..801).map(|i| [-0.008 + 0.016 * i as f64 / 800.0, 0.0]).collect();
    let g = exp.coherence(&pts);
    let zeros: Vec<f64> = (0..800)
        .filter(|&i| g[i] * g[i + 1] < 0.0)
        .map(|i| pts[i][0] - g[i] * (pts[i + 1][0] - pts[i][0]) / (g[i + 1] - g[i]))
        .collect();
    // two zeros per fringe period
    let spacing = 2.0 * (zeros.last().unwrap() - zeros[0]) / (zeros.len() - 1) as f64;
    println!("fringe spacing {:.5} mm from {} zeros (λL/d = {:.5} mm)", spacing * 1e3, zeros.len(), HENE * L / D * 1e3);
    assert!((spacing - HENE * L / D).abs() / (HENE * L / D) < 1e-4);
}

#[test]
fn source_width_sets_visibility_as_van_cittert_zernike_predicts() {
    // Monochromatic incoherent slit source of width s at z_s = 0.5 m:
    // |γ| = |sinc(π s d/(λ z_s))|, with a sign flip past the first zero.
    let zs = 0.5;
    for &s in &[0.2e-3, 0.6e-3, 1.0e-3, 1.6e-3] {
        let src = Source {
            spectrum: Spectrum::Line(HENE),
            extent: Extent::Slit { power: 1e-3, distance: zs, width: s },
            polarization: Some(Jones::horizontal()),
        };
        let exp = Experiment::double_slit(src, A, D, H, L, 1e-4);
        let g = exp.coherence(&[[0.0, 0.0]])[0];
        // Pinhole formula (textbook): sinc(π s d/(λ z_s)).
        let u = PI * s * D / (HENE * zs);
        let pinholes = u.sin() / u;
        // Finite slits: a source point at x_s tilts the light on each slit,
        // shifting its own sinc envelope: with q = 1/z_s + 1/L and slit
        // centres c = ∓d/2, the slit factors are sinc(k a/2 (q c − x_s/z_s))
        // and γ = ∫ cos(k d x_s/z_s) S₁S₂ / √(∫S₁² ∫S₂²).
        let k = 2.0 * PI / HENE;
        let q = 1.0 / zs + 1.0 / L;
        let sinc = |v: f64| if v == 0.0 { 1.0 } else { v.sin() / v };
        let n = 20000;
        let (mut num, mut s1, mut s2) = (0.0, 0.0, 0.0);
        for i in 0..n {
            let xs = -0.5 * s + s * (i as f64 + 0.5) / n as f64;
            let f1 = sinc(k * A / 2.0 * (q * (-0.5 * D) - xs / zs));
            let f2 = sinc(k * A / 2.0 * (q * (0.5 * D) - xs / zs));
            num += (k * D * xs / zs).cos() * f1 * f2;
            s1 += f1 * f1;
            s2 += f2 * f2;
        }
        let finite = num / (s1 * s2).sqrt();
        println!("s = {:.1} mm: γ(0) = {g:+.5}; finite-slit closed form {finite:+.5}; pinhole sinc {pinholes:+.5}", s * 1e3);
        assert!((g - finite).abs() < 1e-3, "s {s}: {g} vs {finite}");
    }
}

#[test]
fn bandwidth_washes_out_high_order_fringes() {
    // Gaussian spectrum (600 nm, 20 nm FWHM), collimated. Reference: the
    // closed-form two-beam model (each slit Fraunhofer across, Fresnel along
    // its height, at every wavelength) integrated over the spectrum with
    // 4001 samples, independent of the source's coherence-based sampling.
    let (center, fwhm) = (600e-9, 20e-9);
    let spec = Spectrum::Gaussian { center, fwhm };
    let src = Source { spectrum: spec.clone(), extent: Extent::Collimated { irradiance: 1.0 }, polarization: Some(Jones::horizontal()) };
    let exp = Experiment::double_slit(src, A, D, H, L, 1e-4);
    let two_beam = |x: f64| -> f64 {
        let std = fwhm / (2.0 * (2.0 * 2f64.ln()).sqrt());
        let (mut cross, mut p1, mut p2) = (0.0, 0.0, 0.0);
        for i in 0..4001 {
            let lam = center - 6.0 * std + 12.0 * std * i as f64 / 4000.0;
            let w = (-0.5 * ((lam - center) / std).powi(2)).exp(); // density in λ, uniform λ steps
            let k = 2.0 * PI / lam;
            let wy = H / (2.0 * lam * L).sqrt();
            let fy = (lam * L / 2.0).sqrt() * fresnel_span(-wy, wy);
            let u = |c: f64| {
                let r = ((x - c).powi(2) + L * L).sqrt();
                let arg = k * A * (x - c) / (2.0 * r);
                let s = if arg == 0.0 { 1.0 } else { arg.sin() / arg };
                Cx::from_polar(1.0, k * r) * fy * (A * s / (lam * r) * (L / r))
            };
            let (u1, u2) = (u(-0.5 * D), u(0.5 * D));
            cross += w * (u1 * u2.conj()).re;
            p1 += w * u1.norm_sqr();
            p2 += w * u2.norm_sqr();
        }
        cross / (p1 * p2).sqrt()
    };
    // (d/a = 5: orders 5, 10, 15 fall on single-slit envelope zeros)
    for &order in &[3.0, 7.0, 12.0] {
        let xm = order * center * L / D;
        let period = center * L / D;
        let pts: Vec<[f64; 2]> = (0..40).map(|i| [xm - 0.5 * period + period * i as f64 / 40.0, 0.0]).collect();
        let g = exp.coherence(&pts);
        let worst = pts.iter().zip(&g).map(|(p, gi)| (gi - two_beam(p[0])).abs()).fold(0.0f64, f64::max);
        let env = g.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        println!(
            "order {order}: |γ| envelope {env:.4} (ideal coherence function {:.4}); max |γ − two-beam reference| {worst:.1e}; {} spectral modes",
            spec.gaussian_visibility(D * xm / L),
            exp.modes_for(&pts).len()
        );
        assert!(worst < 2e-3, "order {order}: {worst}");
    }
}

#[test]
fn which_path_marking_kills_fringes_and_the_eraser_restores_them_only_for_coherent_light() {
    let mut marked = Experiment::double_slit(
        Source { spectrum: Spectrum::Line(HENE), extent: Extent::Collimated { irradiance: 1.0 }, polarization: Some(Jones::linear(PI / 4.0)) },
        A,
        D,
        H,
        L,
        1e-4,
    );
    marked.openings[0].element = JonesMatrix::polarizer(0.0);
    marked.openings[1].element = JonesMatrix::polarizer(PI / 2.0);
    let pts: Vec<[f64; 2]> = (0..30).map(|i| [-0.002 + 0.004 * i as f64 / 30.0, 0.0]).collect();
    let g = marked.coherence(&pts).iter().fold(0.0f64, |m, v| m.max(v.abs()));
    println!("H/V marked slits: max |γ| {g:.2e}");
    assert!(g < 1e-9);
    // Eraser: analyzer at 45° makes the paths indistinguishable again.
    marked.analyzer = Some(JonesMatrix::polarizer(PI / 4.0));
    let g = marked.coherence(&pts).iter().fold(0.0f64, |m, v| m.max(v.abs()));
    println!("with 45° analyzer (diagonal source): max |γ| {g:.4}");
    assert!(g > 0.999);
    // Unpolarized light: the H and V parts were never coherent; erasure fails.
    marked.source.polarization = None;
    let g = marked.coherence(&pts).iter().fold(0.0f64, |m, v| m.max(v.abs()));
    println!("with 45° analyzer (unpolarized source): max |γ| {g:.2e}");
    assert!(g < 1e-9);
    let _ = Cx::new(0.0, 0.0);
}

#[test]
fn angular_spectrum_grid_rung_agrees_with_the_ladder_over_a_sensor() {
    // A 4 × 3 mm sensor-sized grid 0.3 m behind the slits (Fresnel number
    // ≈ 47 along the slit height, so the pattern has near-field edges),
    // narrowband partially coherent light: every grid point compared with
    // the certified ladder.
    let src = Source { spectrum: Spectrum::Gaussian { center: HENE, fwhm: 2e-9 }, extent: Extent::Slit { power: 1e-3, distance: 0.5, width: 0.3e-3 }, polarization: None };
    let exp = Experiment::double_slit(src, A, D, H, 0.3, 1e-3);
    let xs: Vec<f64> = (0..41).map(|i| -2e-3 + 4e-3 * i as f64 / 40.0).collect();
    let ys: Vec<f64> = (0..16).map(|j| -1.5e-3 + 3e-3 * j as f64 / 15.0).collect();
    let pts: Vec<[f64; 2]> = ys.iter().flat_map(|&y| xs.iter().map(move |&x| [x, y])).collect();
    let modes = exp.modes_for(&pts);
    let t = std::time::Instant::now();
    let (grid, audit) = exp.grid_spectral_irradiance(&xs, &ys, &modes, 8);
    let t_grid = t.elapsed().as_secs_f64();
    let t = std::time::Instant::now();
    let ladder = exp.irradiance_modes(&pts, &modes);
    let t_ladder = t.elapsed().as_secs_f64();
    let total: Vec<f64> = (0..pts.len()).map(|q| grid.iter().map(|b| b.1[q]).sum()).collect();
    let peak = ladder.values.iter().cloned().fold(0.0, f64::max);
    let worst = total.iter().zip(&ladder.values).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max) / peak;
    println!("{} modes on {} points: grid rung {t_grid:.2} s (audit {audit:.1e}), ladder {t_ladder:.2} s (rungs {:?}); max |ΔI|/I_peak {worst:.2e}", modes.len(), pts.len(), ladder.rungs);
    assert!(worst < 1e-3 && audit < 1e-3);
}
