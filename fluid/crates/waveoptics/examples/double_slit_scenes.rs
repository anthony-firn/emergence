//! Young's double slit, recorded by physical devices.
//!
//! * `a_bare_sensor`: a CMOS sensor in the screen plane, laser light;
//!   raw frame plus a single-photon build-up sequence.
//! * `b_eye`: a human eye (Escudero-Sanz & Navarro, physiological and
//!   published dispersion) looking through the slits held 20 mm in front
//!   of the cornea at a white LED 2 m away; the retinal image, the cone
//!   mosaic's photoisomerisations and the colour they encode.
//! * `c_camera`: a camera with a Thorlabs AC254-050-A achromat at f/8 and a
//!   CMOS sensor, in the same position as the eye.
//! * `d_screen`: laser fringes on a matte white screen, photographed by
//!   an f/8 lens camera and by a pinhole camera (speckle included).
//!
//! Run: `cargo run --release -p waveoptics --example double_slit_scenes`.
//! Output: `fluid/renders/double_slit/` (PNG + report.txt).

use std::f64::consts::PI;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

use num_complex::Complex64 as Cx;
use waveoptics::color::{write_gray16_png, write_srgb_png, xyz, xyz_to_linear_srgb};
use waveoptics::detector::{single_photon_events, Rng, Sensor, HC};
use waveoptics::devices::{sensor_window, subdivisions, Device, Optics, Retina};
use waveoptics::experiment::Experiment;
use waveoptics::eye::{eye, EyeDispersion};
use waveoptics::imaging::{incoherent_image, isoplanatic_image, Window};
use waveoptics::materials::Material;
use waveoptics::polarization::Jones;
use waveoptics::raytrace::{Surface, System};
use waveoptics::source::{Extent, Source, Spectrum};

const A: f64 = 50e-6; // slit width
const D: f64 = 250e-6; // slit separation
const H: f64 = 3e-3; // slit height
const HENE: f64 = 632.8e-9;

/// He–Ne laser, 1 mW expanded to a 5 mm beam: ≈ 51 W/m² (estimate of a
/// typical teaching setup), horizontally polarized.
fn laser() -> Source {
    Source { spectrum: Spectrum::Line(HENE), extent: Extent::Collimated { irradiance: 1e-3 / (PI * 2.5e-3 * 2.5e-3) }, polarization: Some(Jones::horizontal()) }
}

/// Phosphor-converted white LED (estimate of a typical cool-white
/// spectrum: InGaN pump at 450 nm, 20 nm FWHM; YAG:Ce phosphor at 560 nm,
/// 110 nm FWHM; 1 : 3 by power), 1 mm emitter 2 m away. A ~1 W LED emits
/// ~0.3 W of light; Lambertian on-axis intensity Φ/π equals that of an
/// isotropic 4Φ = 1.2 W source.
fn led() -> Source {
    Source {
        spectrum: Spectrum::Mixture(vec![(0.25, 450e-9, 20e-9), (0.75, 560e-9, 110e-9)]),
        extent: Extent::Disk { power: 1.2, distance: 2.0, diameter: 1e-3 },
        polarization: None,
    }
}

/// Thorlabs AC254-050-A (vendor Zemax prescription), stopped to f/8.
fn achromat() -> System {
    let image_z = 9.0e-3 + 2.5e-3 + 43.22018758686e-3;
    System {
        object_medium: Material::Air,
        surfaces: vec![
            Surface { radius: 1.0 / 2.999400119976e1, conic: 0.0, semi_diameter: 12.7e-3, thickness: 9.0e-3, after: Material::n_baf10() },
            Surface { radius: 1.0 / -4.48833034111311e1, conic: 0.0, semi_diameter: 12.7e-3, thickness: 2.5e-3, after: Material::sf10() },
            Surface { radius: 1.0 / -3.4355996839248e0, conic: 0.0, semi_diameter: 12.7e-3, thickness: 43.22018758686e-3, after: Material::Air },
        ],
        stop_z: 0.0,
        stop_radius: 50e-3 / 16.0,
        image_z,
    }
}

fn out_dir() -> PathBuf {
    let p = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../renders/double_slit");
    std::fs::create_dir_all(&p).unwrap();
    p
}

/// Nearest-neighbour upscale for viewing.
fn upscale<T: Copy>(v: &[T], w: usize, h: usize, k: usize) -> Vec<T> {
    (0..w * k * h * k).map(|q| v[(q % (w * k)) / k + w * ((q / (w * k)) / k)]).collect()
}

/// Colour image of monochromatic light of irradiance `e` (display white =
/// brightest sample).
fn mono_rgb(e: &[f64], lambda: f64) -> (Vec<[f64; 3]>, f64) {
    let c = xyz_to_linear_srgb(xyz(lambda));
    let rgb: Vec<[f64; 3]> = e.iter().map(|v| c.map(|k| k * v)).collect();
    let white = e.iter().cloned().fold(0.0, f64::max) * xyz(lambda)[1];
    (rgb, white)
}

/// Fringe period measured the experimenter's way: locate the dark
/// fringes (minima of the lightly smoothed profile relative to its local
/// mean, parabolic refinement)
/// walking out from the centre one period at a time while the fringes stay
/// bright, and fit a straight line to their positions.
fn period(profile: &[f64], dx: f64, guess: f64) -> f64 {
    let n = profile.len();
    let g = guess / dx; // samples per period
    let k = ((g / 10.0).round() as usize).max(0);
    let sm: Vec<f64> = (0..n).map(|i| {
        let (a, b) = (i.saturating_sub(k), (i + k + 1).min(n));
        profile[a..b].iter().sum::<f64>() / (b - a) as f64
    }).collect();
    let peak = sm[n / 4..3 * n / 4].iter().cloned().fold(0.0, f64::max);
    // divide out the envelope (running mean over one period) so a sloping
    // envelope does not pull the minima of low-visibility fringes outward
    let half = (g / 2.0).round() as usize;
    let norm: Vec<f64> = (0..n)
        .map(|i| {
            let (a, b) = (i.saturating_sub(half), (i + half + 1).min(n));
            let mean = sm[a..b].iter().sum::<f64>() / (b - a) as f64;
            if mean > 0.0 { sm[i] / mean } else { 1.0 }
        })
        .collect();
    // the darkest sample within ±half of c, if it is interior to that
    // window and the fringes around it are bright enough to measure
    let local_min = |c: f64, half: f64| -> Option<f64> {
        let (a, b) = ((c - half).max(1.0) as usize, ((c + half) as usize).min(n - 2));
        if a + 2 >= b {
            return None;
        }
        let i = (a..=b).min_by(|&x, &y| norm[x].total_cmp(&norm[y]))?;
        if i == a || i == b {
            return None;
        }
        let near = (i.saturating_sub(g as usize)..(i + g as usize).min(n)).map(|q| sm[q]).fold(0.0, f64::max);
        // stop where the envelope fades or the fringes wash out (white
        // light beyond its coherence length): a dark fringe must be ≥ 25%
        // below the local mean
        if near < 0.2 * peak || norm[i] > 0.75 {
            return None;
        }
        let (l, c0, r) = (norm[i - 1], norm[i], norm[i + 1]);
        let den = l - 2.0 * c0 + r;
        Some(i as f64 + if den > 0.0 { 0.5 * (l - r) / den } else { 0.0 })
    };
    let Some(start) = local_min(n as f64 / 2.0, 0.75 * g) else { return f64::NAN };
    let mut found = vec![(0.0, start)];
    for dir in [-1.0, 1.0] {
        let (mut idx, mut pos) = (0.0, start);
        while let Some(p) = local_min(pos + dir * g, g / 3.0) {
            idx += dir;
            pos = p;
            found.push((idx, p));
        }
    }
    if found.len() < 3 {
        return f64::NAN;
    }
    let m = found.len() as f64;
    let (sx, sy) = found.iter().fold((0.0, 0.0), |(a, b), (x, y)| (a + x, b + y));
    let (mx, my) = (sx / m, sy / m);
    let (num, den) = found.iter().fold((0.0, 0.0), |(a, b), (x, y)| (a + (x - mx) * (y - my), b + (x - mx) * (x - mx)));
    num / den * dx
}

/// Fringe visibility at the centre of a profile: (max − min)/(max + min)
/// over one period either side.
fn visibility(profile: &[f64], samples_per_period: usize) -> f64 {
    let c = profile.len() / 2;
    let s = &profile[c.saturating_sub(samples_per_period)..(c + samples_per_period).min(profile.len())];
    let (lo, hi) = s.iter().fold((f64::INFINITY, 0.0f64), |(a, b), &v| (a.min(v), b.max(v)));
    (hi - lo) / (hi + lo)
}

/// Mean of rows `j0..j1` of a `w`-wide image.
fn rows(v: &[f64], w: usize, j0: usize, j1: usize) -> Vec<f64> {
    (0..w).map(|i| (j0..j1).map(|j| v[i + w * j]).sum::<f64>() / (j1 - j0) as f64).collect()
}

fn to_f64(v: &[u32]) -> Vec<f64> {
    v.iter().map(|&x| x as f64).collect()
}

/// Auto-exposure: scale the exposure so the brightest pixel's expected
/// charge is `fill` of full well.
fn auto_expose(sensor: &mut Sensor, expected_at_1s: &[f64], fill: f64) -> Vec<f64> {
    let peak = expected_at_1s.iter().cloned().fold(0.0, f64::max);
    sensor.exposure = fill * sensor.full_well / peak;
    expected_at_1s.iter().map(|v| v * sensor.exposure).collect()
}

fn scene_bare_sensor(out: &Path, rng: &mut Rng, rep: &mut String) {
    let t = Instant::now();
    let l = 0.3;
    let exp = Experiment::double_slit(laser(), A, D, H, l, 1e-3);
    // 1/2.3"-class sensor: 640 × 480, 5.86 µm pixels
    let mut sensor = Sensor::cmos([640, 480], 5.86e-6, 1.0);
    sensor.full_well = 30_000.0; // estimate for 5.86 µm pixels
    sensor.gain = sensor.full_well / 4095.0;
    let modes = exp.modes_for(&[[0.0, 0.0], [2e-3, 1.5e-3]]);
    let n = sensor.points_per_axis(exp.fringe_scale(&modes), 1e-3);
    let ([xs, ys], _) = sensor.axis_quadrature(n);
    let (spectral, audit) = exp.grid_spectral_irradiance(&xs, &ys, &modes, 16);
    let at_1s = sensor.expected_electrons_grid(n, &spectral);
    let expected = auto_expose(&mut sensor, &at_1s, 0.8);
    let dn = sensor.expose(&expected, rng);
    let [w, h] = sensor.pixels;
    write_gray16_png(out.join("a_bare_sensor_raw16.png"), w, h, &dn, 4095).unwrap();
    let (rgb, white) = mono_rgb(&to_f64(&dn), HENE);
    write_srgb_png(out.join("a_bare_sensor_srgb.png"), w, h, &rgb, white / 0.8).unwrap();
    // single-photon build-up: detections drawn from the expected-charge
    // map (photon counting, QE included in the map)
    let mut counts = vec![0u32; w * h];
    let mut drawn = 0;
    for (k, total) in [20usize, 200, 2_000, 20_000, 200_000, 2_000_000].iter().enumerate() {
        for p in single_photon_events(&expected, total - drawn, rng) {
            counts[p] += 1;
        }
        drawn = *total;
        let max = *counts.iter().max().unwrap();
        // dilate single events so they stay visible when downscaled
        let vis: Vec<u32> = (0..w * h)
            .map(|q| {
                let (i, j) = ((q % w) as isize, (q / w) as isize);
                let mut m = 0;
                for dj in -1..=1isize {
                    for di in -1..=1isize {
                        let (x, y) = (i + di, j + dj);
                        if x >= 0 && y >= 0 && (x as usize) < w && (y as usize) < h {
                            m = m.max(counts[x as usize + w * y as usize]);
                        }
                    }
                }
                if *total <= 2000 { m.min(1) } else { m }
            })
            .collect();
        write_gray16_png(out.join(format!("a_photons_{k}_{total}.png")), w, h, &vis, if *total <= 2000 { 1 } else { max }).unwrap();
    }
    let prof = rows(&to_f64(&dn), w, h / 2 - 40, h / 2 + 40);
    let p = period(&prof, sensor.pitch, HENE * l / D);
    let _ = writeln!(
        rep,
        "a) bare CMOS at L = {l} m: {}×{} px of {:.2} µm, {n}² Gauss points/pixel (angular-spectrum grid rung, ladder audit {audit:.1e} of the field scale), exposure {:.3} ms; fringe period {:.4} mm (λL/d = {:.4} mm, {:+.2}%), central visibility {:.3}; {:.1} s",
        w,
        h,
        sensor.pitch * 1e6,
        sensor.exposure * 1e3,
        p * 1e3,
        HENE * l / D * 1e3,
        (p / (HENE * l / D) - 1.0) * 100.0,
        visibility(&prof, (HENE * l / D / sensor.pitch) as usize),
        t.elapsed().as_secs_f64()
    );
}

fn scene_eye(out: &Path, rng: &mut Rng, rep: &mut String) {
    let exp = Experiment::double_slit(led(), A, D, H, 1.0, 3e-3);
    for (name, disp) in [("physiological", EyeDispersion::Physiological), ("navarro", EyeDispersion::Navarro)] {
        let t = Instant::now();
        let e = eye(3e-3, disp);
        let dev = Device { optics: Optics::Refractive(e.clone()), front: [0.0, 0.0, 0.02] };
        let window = Window { center: [0.0, 0.0], pixels: [400, 60], pitch: 1.5e-6 };
        let img = dev.spectral_image(&exp, &window);
        let [w, h] = window.pixels;
        // the optical retinal image (colorimetric, noise-free)
        let xyz_img = img.xyz();
        let white = xyz_img.iter().map(|v| v[1]).fold(0.0, f64::max);
        let rgb: Vec<[f64; 3]> = xyz_img.iter().map(|v| xyz_to_linear_srgb(*v)).collect();
        write_srgb_png(out.join(format!("b_eye_{name}_retinal_image.png")), w * 3, h * 3, &upscale(&rgb, w, h, 3), white).unwrap();
        // the cone mosaic: photoisomerisations in 0.1 s and the colour they
        // encode (no perceptual filling-in: the central fovea has no S cones,
        // so the encoded image there carries no blue signal)
        let retina = Retina::default();
        let cones = retina.sample(&img, rng);
        let per = xyz_img.iter().map(|v| v[1]).sum::<f64>() / xyz_img.len() as f64;
        let perceived = retina.perceived_xyz(&cones, &window);
        let rgb: Vec<[f64; 3]> = perceived.iter().map(|v| xyz_to_linear_srgb(*v)).collect();
        write_srgb_png(out.join(format!("b_eye_{name}_cone_encoded.png")), w * 3, h * 3, &upscale(&rgb, w, h, 3), white).unwrap();
        // mosaic view: each cone drawn in its class colour, brightness ∝ catch
        let max_q = cones.iter().map(|c| c.isomerisations).fold(0.0, f64::max);
        let mut mosaic = vec![[0.0; 3]; w * h];
        for c in &cones {
            let i = ((c.position[0] - window.coords(0)[0]) / window.pitch).round() as isize;
            let j = ((c.position[1] - window.coords(1)[0]) / window.pitch).round() as isize;
            let col = [[1.0, 0.15, 0.1], [0.15, 1.0, 0.1], [0.2, 0.3, 1.0]][c.kind];
            for dj in 0..1isize {
                for di in 0..1isize {
                    let (x, y) = (i + di, j + dj);
                    if x >= 0 && y >= 0 && (x as usize) < w && (y as usize) < h {
                        mosaic[x as usize + w * y as usize] = col.map(|k| k * c.isomerisations / max_q);
                    }
                }
            }
        }
        write_srgb_png(out.join(format!("b_eye_{name}_cone_catches.png")), w * 3, h * 3, &upscale(&mosaic, w, h, 3), 1.0).unwrap();
        let f = e.paraxial(555e-9).efl;
        // integrate over all rows: the disk's coherence (van Cittert–
        // Zernike) is a property of the whole source image
        let prof = rows(&img.irradiance(), w, 0, h);
        let expect = f * 550e-9 / D;
        let reference = two_beam(&window.coords(0), f, &exp.source.spectrum, [1e-3, 2.0], |_| 1.0);
        if name == "physiological" {
            // the same retinal image through an aberration-free lens with
            // the eye's entrance pupil and focal length (isolates the eye's
            // own aberrations)
            let par = e.paraxial(555e-9);
            let ideal = Device { optics: Optics::Ideal { radius: par.ep_radius, focal: f }, front: [0.0, 0.0, 0.02 + par.ep_z] };
            let ii = ideal.spectral_image(&exp, &window);
            let pi = rows(&ii.irradiance(), w, 0, h);
            let xyz_i = ii.xyz();
            let wi = xyz_i.iter().map(|v| v[1]).fold(0.0, f64::max);
            let rgb: Vec<[f64; 3]> = xyz_i.iter().map(|v| xyz_to_linear_srgb(*v)).collect();
            write_srgb_png(out.join("b_ideal_lens_retinal_image.png"), w * 3, h * 3, &upscale(&rgb, w, h, 3), wi).unwrap();
            let _ = writeln!(
                rep,
                "b) aberration-free lens with the eye's pupil ({:.3} mm radius) and focal length ({:.3} mm): fringe profile within {:.2e} of the two-beam closed form (relative to peak); centre visibility {:.3} (closed form {:.3})",
                par.ep_radius * 1e3,
                f * 1e3,
                shape_error(&pi, &reference),
                visibility(&pi, (expect / window.pitch) as usize),
                visibility(&reference, (expect / window.pitch) as usize)
            );
        }
        let lm: Vec<f64> = cones.iter().filter(|c| c.kind < 2).map(|c| c.isomerisations).collect();
        let _ = writeln!(
            rep,
            "b) eye ({name} dispersion), slits 20 mm before the cornea, LED 2 m away: {} modes, pupil {}², band loss {:.2}; fringe spacing f·550 nm/d = {:.2} µm = {:.2} arcmin; profile differs from the aberration-free closed form by {:.3} of peak (the eye's own aberrations); centre visibility {:.3}; peak retinal irradiance {:.3e} W/m², mean {per:.3e} (Y-weighted); {} cones, L+M catches up to {:.0} per 0.1 s; {:.1} s",
            img.modes,
            img.pupil_samples,
            img.band_loss,
            expect * 1e6,
            550e-9 / D * 180.0 / PI * 60.0,
            shape_error(&prof, &reference),
            visibility(&prof, (expect / window.pitch) as usize),
            img.irradiance().iter().cloned().fold(0.0, f64::max),
            cones.len(),
            lm.iter().cloned().fold(0.0, f64::max),
            t.elapsed().as_secs_f64()
        );
    }
}

fn scene_camera(out: &Path, rng: &mut Rng, rep: &mut String) {
    let t = Instant::now();
    let exp = Experiment::double_slit(led(), A, D, H, 1.0, 3e-3);
    let lens = achromat();
    let f = lens.paraxial(550e-9).efl;
    let dev = Device { optics: Optics::Refractive(lens), front: [0.0, 0.0, 0.02] };
    let mut sensor = Sensor::cmos([480, 64], 3.45e-6, 1.0);
    let sub = subdivisions(sensor.pitch, f * 450e-9 / D, 1e-3);
    let window = sensor_window(&sensor, sub);
    let img = dev.spectral_image(&exp, &window);
    let at_1s = img.electrons(&sensor, sub);
    let expected = auto_expose(&mut sensor, &at_1s, 0.8);
    let dn = sensor.expose(&expected, rng);
    let [w, h] = sensor.pixels;
    write_gray16_png(out.join("c_camera_raw16.png"), w * 2, h * 2, &upscale(&dn, w, h, 2), 4095).unwrap();
    // ideal colorimetric camera (noise-free XYZ), pixel-averaged
    let xyz_sub = img.xyz();
    let ws = w * sub;
    let xyz_px: Vec<[f64; 3]> = (0..w * h)
        .map(|q| {
            let (i, j) = (q % w, q / w);
            let mut s = [0.0; 3];
            for b in 0..sub {
                for a in 0..sub {
                    let v = xyz_sub[i * sub + a + ws * (j * sub + b)];
                    (0..3).for_each(|k| s[k] += v[k]);
                }
            }
            s
        })
        .collect();
    let white = xyz_px.iter().map(|v| v[1]).fold(0.0, f64::max);
    let rgb: Vec<[f64; 3]> = xyz_px.iter().map(|v| xyz_to_linear_srgb(*v)).collect();
    write_srgb_png(out.join("c_camera_colorimetric.png"), w * 2, h * 2, &upscale(&rgb, w, h, 2), white).unwrap();
    // electron-weighted mean wavelength sets the white-light fringe period
    let (mut sl, mut s) = (0.0, 0.0);
    for (l, e) in &img.bands {
        let q = e.iter().sum::<f64>() * l / HC * sensor.qe.at(*l);
        sl += q / l;
        s += q;
    }
    let lam_eff = s / sl;
    let prof = rows(&to_f64(&dn), w, 0, h);
    let clean = rows(&expected, w, 0, h);
    let expect = f * lam_eff / D;
    let xs: Vec<f64> = (0..w).map(|i| sensor.pixel_center(i, 0)[0]).collect();
    let qe = sensor.qe;
    let reference = two_beam(&xs, f, &exp.source.spectrum, [1e-3, 2.0], |l| l * qe.at(l));
    let _ = writeln!(
        rep,
        "c) AC254-050-A at f/8 + CMOS ({}×{} px, {:.2} µm, {sub}² samples/px), slits 20 mm in front: {} modes, pupil {}², exposure {:.1} ms; expected-charge profile within {:.2e} of the two-beam closed form (relative to peak; f·λ_eff/d = {:.2} µm, λ_eff = {:.1} nm), centre visibility in the raw frame {:.3} (closed form {:.3}); {:.1} s",
        w,
        h,
        sensor.pitch * 1e6,
        img.modes,
        img.pupil_samples,
        sensor.exposure * 1e3,
        shape_error(&clean, &reference),
        expect * 1e6,
        lam_eff * 1e9,
        visibility(&prof, (expect / sensor.pitch) as usize),
        visibility(&reference, (expect / sensor.pitch) as usize),
        t.elapsed().as_secs_f64()
    );
}

/// Largest deviation of `profile` from `reference` after least-squares
/// scaling, relative to the scaled reference's peak.
fn shape_error(profile: &[f64], reference: &[f64]) -> f64 {
    let s = profile.iter().zip(reference).map(|(a, b)| a * b).sum::<f64>() / reference.iter().map(|b| b * b).sum::<f64>();
    let peak = reference.iter().cloned().fold(0.0, f64::max) * s;
    profile.iter().zip(reference).map(|(a, b)| (a - s * b).abs()).fold(0.0, f64::max) / peak
}

/// Independent closed-form reference for a device focused at infinity:
/// the two-beam Fraunhofer pattern summed over the source spectrum and
/// integrated over the slit-height direction (`|U|² ∝ (a h/λf)²` times the
/// height pattern's width `λ f/h`),
/// `Σ w/λ · det(λ) · sinc²(π a θ/λ) [1 + γ(λ) cos(2π d θ/λ)]` at
/// `θ = x'/f`, with the uniform-disk (van Cittert–Zernike) coherence
/// `γ = 2J₁(v)/v`, `v = π D_s d/(λ z_s)`, and `det(λ)` the detector's
/// weighting (1 for irradiance, `λ·QE` for photo-electrons).
fn two_beam(xs: &[f64], f: f64, spectrum: &Spectrum, disk: [f64; 2], det: impl Fn(f64) -> f64) -> Vec<f64> {
    let samples = spectrum.samples(D * 0.05, 1e-4);
    let sinc = |u: f64| if u.abs() < 1e-12 { 1.0 } else { u.sin() / u };
    xs.iter()
        .map(|&x| {
            let th = x / f;
            samples
                .iter()
                .map(|sp| {
                    let l = sp.wavelength;
                    let v = PI * disk[0] * D / (l * disk[1]);
                    let gamma = 2.0 * j1(v) / v;
                    sp.weight / l * det(l) * sinc(PI * A * th / l).powi(2) * (1.0 + gamma * (2.0 * PI * D * th / l).cos())
                })
                .sum()
        })
        .collect()
}

/// Bessel J₁ by its integral representation (64-point midpoint rule;
/// exact to ~1e-12 for the small arguments used here).
fn j1(x: f64) -> f64 {
    (0..64).map(|i| {
        let t = PI * (i as f64 + 0.5) / 64.0;
        (t - x * t.sin()).cos()
    }).sum::<f64>() / 64.0
}

/// Contrast of fully developed speckle integrated over a square pixel of
/// side `p` (Goodman, Speckle Phenomena in Optics, §4.6):
/// C² = ∫∫ (1 − |Δx|/p)(1 − |Δy|/p) |μ(Δ)|² dΔ / p², μ = 2J₁(v)/v,
/// v = 2π a |Δ|/(λ z) for a pupil of radius `a` at distance `z` from the
/// detector.
fn speckle_contrast(p: f64, a: f64, z: f64, lambda: f64) -> f64 {
    let n = 200;
    let h = 2.0 * p / n as f64;
    let mut s = 0.0;
    for j in 0..n {
        for i in 0..n {
            let (dx, dy) = (-p + h * (i as f64 + 0.5), -p + h * (j as f64 + 0.5));
            let v = 2.0 * PI * a * dx.hypot(dy) / (lambda * z);
            let mu = if v < 1e-9 { 1.0 } else { 2.0 * j1(v) / v };
            s += (1.0 - dx.abs() / p) * (1.0 - dy.abs() / p) * mu * mu * h * h;
        }
    }
    (s / (p * p)).sqrt()
}

fn scene_screen(out: &Path, rng: &mut Rng, rep: &mut String) {
    let l = 1.0;
    let exp = Experiment::double_slit(laser(), A, D, H, l, 1e-3);
    let mode = exp.modes_for(&[[0.0, 0.0]])[0];
    // matte white paint: albedo 0.9 (estimate); surface relief ≫ λ, so the
    // reflected phase is uniformly random on scales finer than either
    // camera resolves (fully developed speckle; Goodman, Speckle Phenomena
    // in Optics, ch. 4)
    let rho: f64 = 0.9;
    let _ = writeln!(rep, "d) laser fringes on a matte screen (albedo {rho}) at L = {l} m");
    // (name, camera distance, aperture radius, detector distance, lens?,
    //  pixel pitch, screen field of view)
    for (name, dist, radius, f_det, lens, pitch, fov) in [
        ("lens_f8", 0.5, 50e-3 / 16.0, 1.0 / (1.0 / 50e-3 - 1.0 / 0.5), true, 3.45e-6, [24e-3, 6e-3]),
        ("pinhole", 0.15, 0.125e-3, 50e-3, false, 5.86e-6, [40e-3, 8e-3]),
    ] {
        let t = Instant::now();
        let k = 2.0 * PI / HENE;
        let m = f_det / dist;
        let fc = radius / (HENE * dist);
        // screen sampling: k_px × k_px samples per pixel (no interpolation),
        // fine enough for the image intensity (band 2 f_c) to be alias-free
        let pixel_obj = pitch / m;
        // ≥ 6 samples per finest intensity period 1/(2 f_c), at most 8
        let kpx = ((pixel_obj * 12.0 * fc).ceil() as usize).clamp(1, 8);
        let delta = pixel_obj / kpx as f64;
        let (px, py) = ((fov[0] / pixel_obj) as usize & !1, (fov[1] / pixel_obj) as usize & !1);
        let (nx, ny) = (px * kpx, py * kpx);
        let xs: Vec<f64> = (0..nx).map(|i| (i as f64 + 0.5 - 0.5 * nx as f64) * delta).collect();
        let ys: Vec<f64> = (0..ny).map(|j| (j as f64 + 0.5 - 0.5 * ny as f64) * delta).collect();
        let (fields, audit) = exp.grid_fields(l, &xs, &ys, HENE, &mode.incident, 16);
        let screen: Vec<Cx> = (0..nx * ny).map(|q| fields.iter().map(|f| f[q]).sum::<Cx>() * mode.weight.sqrt()).collect();
        drop(fields);
        let e_screen: Vec<f64> = screen.iter().map(|u| rho * u.norm_sqr()).collect();
        let rough: Vec<Cx> = screen.iter().map(|u| u * rho.sqrt() * Cx::from_polar(1.0, 2.0 * PI * rng.uniform())).collect();
        drop(screen);
        let pupil = |u: [f64; 2]| {
            let r2 = u[0] * u[0] + u[1] * u[1];
            if r2 > radius * radius {
                Cx::new(0.0, 0.0)
            } else if lens {
                Cx::new(1.0, 0.0) // focused: object and detector curvatures cancel the lens
            } else {
                Cx::from_polar(1.0, k * r2 * (1.0 / dist + 1.0 / f_det) / 2.0)
            }
        };
        // Lambertian camera equation: E_img = ρ E_s (a/D)²/m²; the random
        // phasor field puts the share π f_c² δ² of its power in the pupil.
        let norm = (radius / dist).powi(2) / (m * m) / (PI * fc * fc * delta * delta);
        let speckled = isoplanatic_image(&rough, nx, ny, delta, HENE, dist, pupil);
        let mean_img = incoherent_image(&e_screen, nx, ny, delta, HENE, dist, pupil);
        let block = |v: &[f64]| -> Vec<f64> {
            (0..px * py)
                .map(|q| {
                    let (i, j) = (q % px, q / px);
                    let mut s = 0.0;
                    for b in 0..kpx {
                        for a in 0..kpx {
                            s += v[i * kpx + a + nx * (j * kpx + b)];
                        }
                    }
                    s * norm / (kpx * kpx) as f64
                })
                .collect()
        };
        let (e_px, smooth) = (block(&speckled), block(&mean_img));
        let mut sensor = Sensor::cmos([px, py], pitch, 1.0);
        let per_s: Vec<f64> = e_px.iter().map(|e| e * pitch * pitch * HENE / HC * sensor.qe.at(HENE)).collect();
        let expected = auto_expose(&mut sensor, &per_s, 0.8);
        let dn = sensor.expose(&expected, rng);
        write_gray16_png(out.join(format!("d_screen_{name}_raw16.png")), px, py, &dn, 4095).unwrap();
        let (rgb, white) = mono_rgb(&to_f64(&dn), HENE);
        write_srgb_png(out.join(format!("d_screen_{name}_srgb.png")), px, py, &rgb, white / 0.8).unwrap();
        let (rgb, white) = mono_rgb(&smooth, HENE);
        write_srgb_png(out.join(format!("d_screen_{name}_speckle_free.png")), px, py, &rgb, white).unwrap();
        let prof = rows(&to_f64(&dn), px, py / 2 - py / 8, py / 2 + py / 8);
        let expect = m * HENE * l / D;
        let p = period(&prof, pitch, expect);
        // speckle contrast: speckled / speckle-free (its ensemble mean) over
        // the bright fringes, against theory for fully developed speckle
        // integrated over a square pixel
        let top = smooth.iter().cloned().fold(0.0, f64::max);
        let ratio: Vec<f64> = e_px.iter().zip(&smooth).filter(|(_, s)| **s > 0.5 * top).map(|(e, s)| e / s).collect();
        let mean = ratio.iter().sum::<f64>() / ratio.len() as f64;
        let c = (ratio.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / ratio.len() as f64).sqrt() / mean;
        let c_theory = speckle_contrast(pitch, radius, f_det, HENE);
        let speckle = HENE * f_det / (2.0 * radius);
        let n_speckles = ratio.len() as f64 * (pitch / speckle).powi(2).min(1.0);
        let p_mean = period(&rows(&smooth, px, py / 2 - py / 8, py / 2 + py / 8), pitch, expect);
        let _ = writeln!(
            rep,
            "   {name}: D = {dist} m, aperture radius {:.3} mm, magnification {m:.4}, {px}×{py} px of {:.2} µm ({kpx}² screen samples/px, ladder audit {audit:.1e}), exposure {:.1} ms; fringe period {:.2} µm in the frame, {:.2} µm in the speckle-free mean (m λL/d = {:.2} µm, {:+.2}% / {:+.2}%); speckle contrast {:.3} ± {:.3} over ~{:.0} speckles, mean/ensemble mean {:.3} (theory {:.3} for speckle {:.1} µm on {:.2} µm pixels); {:.1} s",
            radius * 1e3,
            pitch * 1e6,
            sensor.exposure * 1e3,
            p * 1e6,
            p_mean * 1e6,
            expect * 1e6,
            (p / expect - 1.0) * 100.0,
            (p_mean / expect - 1.0) * 100.0,
            c,
            c / n_speckles.sqrt(),
            n_speckles,
            mean,
            c_theory,
            speckle * 1e6,
            pitch * 1e6,
            t.elapsed().as_secs_f64()
        );
    }
}

fn main() {
    let out = out_dir();
    let mut rng = Rng(2026);
    let mut rep = String::new();
    let only: Vec<String> = std::env::args().skip(1).collect();
    let run = |s: &str| only.is_empty() || only.iter().any(|o| o == s);
    if run("a") {
        scene_bare_sensor(&out, &mut rng, &mut rep);
    }
    if run("b") {
        scene_eye(&out, &mut rng, &mut rep);
    }
    if run("c") {
        scene_camera(&out, &mut rng, &mut rep);
    }
    if run("d") {
        scene_screen(&out, &mut rng, &mut rep);
    }
    print!("{rep}");
    std::fs::write(out.join("report.txt"), rep).unwrap();
}
