//! The optical ladder: is the indicator calibrated, do the wave rungs
//! reproduce the analytic diffraction, and do the wave-optical surfaces
//! reproduce the closed forms?

use std::f64::consts::PI;

use emergence::{Request, Rigor};
use num_complex::Complex64 as Cx;
use photonics::multilayer::Polarization::{TE, TM};
use viz::wave::{
    airy_first_zero, airy_reflectance, brewster_angle, coherence_length, combine_paths,
    double_slit, expected_electrons, fresnel_number, geometric_irradiance, irradiance,
    irradiance_at, irradiance_ladder, knife_edge, knife_edge_geometric, photon_count,
    shadow_blur_pixels, shadow_claim, shadow_misplaced_irradiance, single_slit, visibility,
    worst_shadow_error, Film,
};
use waveoptics::diffraction::{Incident, Query, Rect};

fn plane_wave_query(half_width: f64, z: f64, x: f64, lambda: f64) -> Query {
    Query {
        aperture: Rect {
            center: [0.0, 0.0],
            size: [2.0 * half_width, 2.0 * half_width],
        },
        incident: Incident::Plane {
            amplitude: Cx::new(1.0, 0.0),
            sines: [0.0, 0.0],
        },
        point: [x, 0.0, z],
        wavelength: lambda,
    }
}

#[test]
fn the_exact_knife_edge_matches_its_published_values() {
    // The canonical scalar diffraction solution, and the reference the
    // renderer's indicator is calibrated against. Two published numbers
    // pin it: the irradiance at the geometric edge is exactly a quarter of
    // the open value, and the first maximum is 1.370 at v = 1.217
    // (Born & Wolf, Principles of Optics, 7th ed., table 8.1).
    assert!((knife_edge(0.0) - 0.25).abs() < 1e-12, "I(0) = {:.9}", knife_edge(0.0));
    let mut best = (0.0f64, 0.0f64);
    for i in 0..20000 {
        let v = 4.0 * i as f64 / 20000.0;
        let x = knife_edge(v);
        if x > best.1 {
            best = (v, x);
        }
    }
    println!("first maximum {:.6} at v = {:.4}", best.1, best.0);
    assert!((best.1 - 1.3704).abs() < 1e-3, "first maximum {:.6}", best.1);
    assert!((best.0 - 1.217).abs() < 2e-3, "at v = {:.4}", best.0);
    // Deep in the light it tends to one, deep in the shadow to zero, and
    // it is monotone in neither — that is the fringe pattern.
    assert!((knife_edge(40.0) - 1.0).abs() < 0.02);
    assert!(knife_edge(-40.0) < 1e-3);
    assert!(knife_edge(-1.0) < 0.1 && knife_edge(1.0) > 1.0);
}

#[test]
fn diffraction_misplaces_no_net_irradiance() {
    // The signed error of a sharp shadow integrates away: diffraction
    // moves irradiance across the edge — a 37 % overshoot on the lit side,
    // a tail into the shadow — but creates none. The integral decays with
    // the range, halving as the range doubles, because the profiles agree
    // only as 1/v.
    let mut last = f64::INFINITY;
    for reach in [40.0, 80.0, 160.0, 320.0] {
        let k = shadow_misplaced_irradiance(reach, 0.004);
        println!("net misplaced irradiance over +-{reach}: {k:.6}");
        assert!(k < last, "it should keep falling: {k} after {last}");
        last = k;
    }
    assert!(last < 5e-4, "the net should be heading to zero: {last:.2e}");
}

#[test]
fn the_shadow_indicator_bounds_the_pixel_error() {
    // The claim the renderer rests on: the pixel-averaged error of a sharp
    // shadow is bounded by 0.25 * blur below one pixel, with a ceiling of
    // 3/4. Measured against the exact knife edge over four decades.
    let lambda = 550e-9;
    let pixel = 1e-3;
    let mut rows = Vec::new();
    for &z in &[1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2] {
        let blur = shadow_blur_pixels(z, lambda, pixel);
        let err = worst_shadow_error(pixel, z, lambda);
        let claim = shadow_claim(z, lambda, pixel);
        println!(
            "blur {blur:9.5} px: worst pixel error {err:9.6}, claimed {:9.6} ({})",
            claim.error,
            claim.rigor.name()
        );
        assert!(
            err <= claim.error * 1.001,
            "blur {blur}: measured {err} exceeds the claim {}",
            claim.error
        );
        rows.push((blur, err));
    }
    // Monotone in the blur, as the indicator asserts.
    for w in rows.windows(2) {
        assert!(w[1].1 > w[0].1, "error should grow with blur: {rows:?}");
    }
    // A sharp shadow is far better than it looks once averaged: a quarter
    // of a pixel of blur costs three per cent, a thousandth of a pixel
    // costs under 1e-4.
    let fine = rows.iter().find(|(b, _)| *b < 0.01).expect("a fine case");
    assert!(fine.1 < 1e-4, "blur {:.5} px gives {:.2e}", fine.0, fine.1);
    let (b, e) = rows
        .iter()
        .find(|(b, _)| *b > 0.2 && *b < 0.3)
        .expect("a quarter-pixel case");
    assert!(*e < 0.05, "blur {b:.3} px gives {e:.4}");

    // Pointwise, by contrast, the error at the edge never shrinks at all:
    // it is 3/4 however small the blur, which is exactly the ceiling.
    let pointwise = (knife_edge_geometric(1e-12) - knife_edge(0.0)).abs();
    assert!((pointwise - 0.75).abs() < 1e-12, "pointwise edge error {pointwise}");
    assert!((viz::wave::SHADOW_ERROR_CEILING - 0.75).abs() < 1e-18);
    // The claim saturates there, as a rigorous bound rather than an
    // estimate, and the measured error never exceeds it.
    let huge = shadow_claim(1e4, lambda, pixel);
    assert_eq!(huge.rigor, Rigor::Rigorous);
    assert!((huge.error - 0.75).abs() < 1e-18);
    assert!(rows.last().unwrap().1 < 0.75);
    // The Fresnel scale itself: 74 um at 10 mm, 2.3 mm at 10 m.
    assert!((viz::wave::fresnel_scale(1e-2, lambda) - 74.16e-6).abs() < 1e-7);
    assert!((viz::wave::fresnel_scale(1e1, lambda) - 2.345e-3).abs() < 1e-6);
}

#[test]
fn the_ladder_picks_the_cheapest_rung_that_can_meet_the_tolerance() {
    let lambda = 550e-9;
    // Deep in the geometric regime: a 10 mm opening at 10 mm, where the
    // penumbra is a thousandth of the aperture.
    let near = plane_wave_query(5e-3, 1e-2, 0.0, lambda);
    let f = fresnel_number(&near.aperture, near.point, lambda);
    let ladder = irradiance_ladder(1);
    let loose = ladder.solve(&near, Request::new(0.02));
    println!(
        "F = {f:.0}: tolerance 2 % -> rung {} ({}), claim {:.2e} {}",
        loose.rung,
        ladder.rung_name(loose.rung),
        loose.claim.error,
        loose.claim.rigor.name()
    );
    assert_eq!(loose.rung, 0, "a 2 % tolerance at F = {f:.0} needs only rays");
    assert!(loose.met);

    // The far-field rungs are *refused* here, and by their own rigorous
    // bounds rather than by a hand-written domain test: Fraunhofer and
    // Fresnel assume a paraxial aperture, and `waveoptics` proves how
    // badly that fails in the near field.
    for rung in 1..=2 {
        match ladder.assess(rung, &near) {
            emergence::Assessment::Claim(c) => assert!(
                c.error > 1.0,
                "rung {} should be hopeless at F = {f:.0}, claims {:.2e}",
                ladder.rung_name(rung),
                c.error
            ),
            other => panic!("rung {rung}: {other:?}"),
        }
    }
    // So a tolerance rays cannot meet has nowhere to go but the
    // reference. That is asserted from the claims rather than by running
    // it: Rayleigh-Sommerfeld in the near field is exactly the expensive
    // case the ladder exists to avoid — it is unaffordable at F = 4500,
    // which is why the geometric rung is worth having.
    let request = Request::new(1e-6);
    for rung in 0..3 {
        match ladder.assess(rung, &near) {
            emergence::Assessment::Claim(c) => assert!(
                !c.meets(&request),
                "rung {} should not meet 1e-6 at F = {f:.0}",
                ladder.rung_name(rung)
            ),
            other => panic!("rung {rung}: {other:?}"),
        }
    }
    match ladder.assess(3, &near) {
        emergence::Assessment::Claim(c) => assert!(c.meets(&request)),
        other => panic!("the reference must serve: {other:?}"),
    }
    // The reference is the expensive rung, which is the whole point.
    assert!(loose.cost < 10.0, "the geometric rung is cheap: {}", loose.cost);

    // Deep in the far field the picture inverts: now Fraunhofer and
    // Fresnel are cheap *and* proved, and rays are the hopeless one.
    let far = plane_wave_query(100e-6, 2.0, 0.0, lambda);
    let f = fresnel_number(&far.aperture, far.point, lambda);
    assert!(f < 0.1, "F = {f:e} should be small");
    let solved = ladder.solve(&far, Request::new(1e-3));
    println!(
        "F = {f:.4}: tolerance 1e-3 -> rung {} ({}), claim {:.2e} {}",
        solved.rung,
        ladder.rung_name(solved.rung),
        solved.claim.error,
        solved.claim.rigor.name()
    );
    assert!(solved.rung > 0 && solved.rung < 3, "a far-field wave rung should serve");
    assert_eq!(solved.claim.rigor, Rigor::Rigorous, "and prove its own error");
    assert!(solved.met);
    // The geometric rung is at its ceiling there.
    match ladder.assess(0, &far) {
        emergence::Assessment::Claim(c) => {
            assert!(c.error >= 0.74, "rays should be at the 3/4 ceiling: {:.3}", c.error)
        }
        other => panic!("{other:?}"),
    }

    // And the convenience wrapper agrees with the ladder it builds.
    let direct = irradiance_at(&far, 1e-3, 1);
    assert_eq!(direct.rung, solved.rung);
}

#[test]
fn the_wave_rungs_reproduce_analytic_diffraction() {
    // A slit's first minimum is at sin(theta) = lambda/a, exactly. Check
    // the ladder's far-field answer against the closed form.
    let lambda = 633e-9;
    let a = 40e-6;
    let z = 2.0; // far field: a^2/(lambda z) = 1.3e-3
    let first_min = lambda / a;

    // Scan across the pattern and find the minimum nearest the axis.
    let mut best = (f64::INFINITY, 0.0);
    for i in 1..400 {
        let sin_theta = first_min * 0.4 + first_min * 1.2 * i as f64 / 400.0;
        let x = z * sin_theta;
        // A slit: wide in y, narrow in x.
        let qr = Query {
            aperture: Rect { center: [0.0, 0.0], size: [a, 1.0] },
            incident: Incident::Plane { amplitude: Cx::new(1.0, 0.0), sines: [0.0, 0.0] },
            point: [x, 0.0, z],
            wavelength: lambda,
        };
        let value = irradiance(waveoptics::diffraction::fraunhofer(&qr));
        if value < best.0 {
            best = (value, sin_theta);
        }
    }
    let rel = (best.1 - first_min).abs() / first_min;
    println!("slit first minimum at sin(theta) = {:.6e}, analytic {first_min:.6e} ({rel:.2e})", best.1);
    assert!(rel < 0.01, "first minimum off by {rel:.3e}");

    // The analytic helpers themselves: a sinc^2 with its zeros in the
    // right places, and two-slit fringes spaced lambda/d.
    assert!((single_slit(0.0, a, lambda) - 1.0).abs() < 1e-15);
    assert!(single_slit(lambda / a, a, lambda) < 1e-25, "the first zero must vanish");
    assert!(single_slit(2.0 * lambda / a, a, lambda) < 1e-25);
    let d = 200e-6;
    // Fringe maxima at sin(theta) = m lambda / d.
    for m in 0..4 {
        let s = m as f64 * lambda / d;
        let here = double_slit(s, a, d, lambda);
        let between = double_slit(s + 0.5 * lambda / d, a, d, lambda);
        assert!(here > between, "fringe {m} should be a maximum");
        assert!(between < 1e-20, "the midpoint between fringes must vanish");
    }
    // Airy's first zero is 1.2197 lambda/D, not lambda/D.
    let airy = airy_first_zero(1e-3, lambda);
    assert!((airy / (lambda / 1e-3) - 1.219_669_891).abs() < 1e-9);
}

#[test]
fn a_thin_film_reproduces_the_two_beam_closed_form() {
    // The transfer matrix must agree with the Airy formula for a single
    // non-absorbing layer, exactly — every internal reflection summed.
    let (n0, n1, n2) = (1.0, 1.33, 1.0); // a soap film in air
    for &thickness in &[100e-9, 300e-9, 550e-9, 1200e-9] {
        let film = Film::single(n1, thickness, n2);
        let mut worst: f64 = 0.0;
        for i in 0..80 {
            let lambda = 380e-9 + (780e-9 - 380e-9) * i as f64 / 80.0;
            let got = film.reflectance(lambda, 0.0, TE);
            let want = airy_reflectance(n0, n1, n2, thickness, lambda);
            worst = worst.max((got - want).abs());
        }
        assert!(worst < 1e-12, "d = {thickness:e}: transfer matrix vs Airy {worst:e}");
    }

    // The reflectance maxima of a film in air sit where
    // 2 n d = (m + 1/2) lambda, because the two interfaces' amplitude
    // coefficients have opposite signs.
    let thickness = 300e-9;
    let film = Film::single(n1, thickness, n2);
    let m = 0.0;
    let lambda_peak = 4.0 * n1 * thickness / (2.0 * m + 1.0);
    assert!(
        lambda_peak > 380e-9 && lambda_peak < 1800e-9,
        "peak at {lambda_peak:e}"
    );
    let at_peak = film.reflectance(lambda_peak, 0.0, TE);
    for delta in [-0.15, 0.15] {
        let off = film.reflectance(lambda_peak * (1.0 + delta), 0.0, TE);
        assert!(at_peak > off, "the quarter-wave peak must be a maximum");
    }

    // Energy: a non-absorbing film reflects and transmits everything.
    for &lambda in &[450e-9, 550e-9, 650e-9] {
        let s = film.stack.solve(lambda, 0.3, TM);
        assert!(
            (s.reflectance + s.transmittance - 1.0).abs() < 1e-12,
            "R + T = {} at {lambda:e}",
            s.reflectance + s.transmittance
        );
        assert!(s.absorptance.abs() < 1e-12);
    }

    // And the film has a colour that no scalar reflectance can produce:
    // a 300 nm soap film is not grey.
    let rgb = film.reflected_rgb(0.0, 64);
    let (lo, hi) = (
        rgb.iter().copied().fold(f64::INFINITY, f64::min),
        rgb.iter().copied().fold(0.0, f64::max),
    );
    println!("300 nm soap film in reflection: RGB {rgb:?}");
    assert!(hi > 0.0 && (hi - lo) / hi > 0.2, "the film must be coloured, not grey: {rgb:?}");
    // A film a whole wave thick at mid-visible is much less saturated.
    let thick = Film::single(n1, 10e-6, n2).reflected_rgb(0.0, 64);
    let span = |c: [f64; 3]| {
        let hi = c.iter().copied().fold(0.0, f64::max);
        let lo = c.iter().copied().fold(f64::INFINITY, f64::min);
        (hi - lo) / hi.max(1e-30)
    };
    assert!(span(thick) < span(rgb), "a thick film washes out: {thick:?} vs {rgb:?}");
}

#[test]
fn brewsters_angle_is_an_exact_zero() {
    // p-polarised light is not reflected at all at tan(theta) = n2/n1.
    // One of the sharpest checks there is on a Fresnel implementation,
    // and the renderer gets it from `photonics` rather than its own.
    let (n1, n2) = (1.0, 1.5);
    let theta_b = brewster_angle(n1, n2);
    assert!((theta_b.to_degrees() - 56.3099).abs() < 1e-3, "{:.4} deg", theta_b.to_degrees());
    // A bare interface as a film of negligible thickness.
    let bare = Film::from_stack(photonics::multilayer::Stack::new(
        n1,
        vec![photonics::multilayer::Layer::new(n2, 0.0, 1e-9)],
        n2,
    ));
    let r_p = bare.reflectance(550e-9, theta_b, TM);
    println!("Rp at Brewster: {r_p:e}");
    assert!(r_p < 1e-6, "p-polarised reflectance should vanish: {r_p:e}");
    // s-polarised light is strongly reflected at the same angle.
    let r_s = bare.reflectance(550e-9, theta_b, TE);
    assert!(r_s > 0.1, "Rs at Brewster {r_s}");
    // Away from Brewster, Rp recovers.
    assert!(bare.reflectance(550e-9, theta_b * 0.5, TM) > r_p);
    assert!(bare.reflectance(550e-9, 1.4, TM) > 0.1);
    // Unpolarised is the mean of the two.
    let un = bare.unpolarised(550e-9, theta_b);
    assert!((un - 0.5 * (r_s + r_p)).abs() < 1e-15);
}

#[test]
fn partial_coherence_interpolates_between_waves_and_rays() {
    // Two paths of equal amplitude, in phase. Fully coherent they give
    // four times one path; fully incoherent, twice. A renderer that only
    // adds intensities has silently chosen the second.
    let a = Cx::new(1.0, 0.0);
    assert!((combine_paths(a, a, 1.0) - 4.0).abs() < 1e-15);
    assert!((combine_paths(a, a, 0.0) - 2.0).abs() < 1e-15);
    assert!((combine_paths(a, a, 0.5) - 3.0).abs() < 1e-15);
    // In antiphase, coherence cancels them completely.
    let b = Cx::new(-1.0, 0.0);
    assert!(combine_paths(a, b, 1.0).abs() < 1e-15);
    assert!((combine_paths(a, b, 0.0) - 2.0).abs() < 1e-15);

    // The visibility of a real source falls with the path difference, and
    // its coherence length is the scale on which that happens. A narrow
    // line stays coherent much further than a broad one.
    let narrow = waveoptics::source::Spectrum::Gaussian {
        center: 633e-9,
        fwhm: 0.01e-9,
    };
    let broad = waveoptics::source::Spectrum::Gaussian {
        center: 550e-9,
        fwhm: 100e-9,
    };
    let (lc_n, lc_b) = (coherence_length(&narrow), coherence_length(&broad));
    println!("coherence length: narrow {:.3e} m, broad {:.3e} m", lc_n, lc_b);
    assert!(lc_n > 1e-3, "a 0.01 nm line should stay coherent over millimetres");
    assert!(lc_b < 1e-5, "a 100 nm band should decohere within microns");
    assert!(lc_n > 100.0 * lc_b);
    // Visibility is one at zero path difference and falls monotonically.
    assert!((visibility(&broad, 0.0) - 1.0).abs() < 1e-12);
    let mut last = 1.0;
    for i in 1..20 {
        let v = visibility(&broad, lc_b * i as f64 / 5.0);
        assert!(v <= last + 1e-12, "visibility must not rise");
        last = v;
    }
    assert!(last < 0.05, "far beyond the coherence length nothing interferes: {last}");
}

#[test]
fn photon_counting_has_poisson_statistics() {
    // A render is a measurement, and a measurement of light is a count.
    // At low exposure the image *is* the noise, and that is physics, not
    // sampling error.
    let mut rng = waveoptics::detector::Rng(0x1234_5678_9ABC_DEF0);
    for &mean in &[0.5, 5.0, 50.0, 500.0] {
        let n = 20000;
        let mut sum = 0.0;
        let mut sq = 0.0;
        for _ in 0..n {
            let c = photon_count(mean, &mut rng);
            assert!(c >= 0.0 && c.fract() == 0.0, "a count must be a non-negative integer");
            sum += c;
            sq += c * c;
        }
        let m = sum / n as f64;
        let var = sq / n as f64 - m * m;
        // Poisson: variance equals the mean. The sampling error of the
        // variance itself is about mean*sqrt(2/n).
        let tol = 6.0 * mean * (2.0 / n as f64).sqrt() + 0.05 * mean;
        println!("mean {mean}: measured {m:.3}, variance {var:.3}");
        assert!((m - mean).abs() < 5.0 * (mean / n as f64).sqrt() + 1e-3);
        assert!((var - mean).abs() < tol, "variance {var} vs mean {mean}");
    }

    // The exposure that makes a photon count meaningful. A 550 nm photon
    // carries 3.6e-19 J. Daylight is about 1e7 W/(m^2 sr m) of spectral
    // radiance, so a 1 um pixel over a 100 nm band and a steradian,
    // exposed for a millisecond at half quantum efficiency, collects a few
    // thousand electrons — which is what a real sensor does.
    let n = expected_electrons(1e7, 550e-9, 100e-9, 1.0, 1e-12, 1e-3, 0.5);
    println!("expected electrons: {n:.1}");
    assert!(n > 100.0 && n < 1e6, "{n}");
    // Twice the exposure, twice the count, exactly.
    let n2 = expected_electrons(1e7, 550e-9, 100e-9, 1.0, 1e-12, 2e-3, 0.5);
    assert!((n2 / n - 2.0).abs() < 1e-12);
    // And the shot-noise signal-to-noise goes as the square root.
    assert!(((n2 / n).sqrt() - 2.0f64.sqrt()).abs() < 1e-12);
}

#[test]
fn the_geometric_rung_is_the_shadow_and_nothing_more() {
    let lambda = 550e-9;
    let qr = plane_wave_query(1e-3, 0.1, 0.0, lambda);
    // On axis, inside the opening: the full incident irradiance.
    assert!((geometric_irradiance(&qr) - 1.0).abs() < 1e-15);
    // Outside the projected opening: nothing, with no penumbra at all.
    let outside = plane_wave_query(1e-3, 0.1, 1.5e-3, lambda);
    assert!(geometric_irradiance(&outside).abs() < 1e-18);
    // A tilted plane wave shifts the shadow by z tan(theta), and the
    // geometric rung follows it exactly.
    let sines = [0.1, 0.0];
    let tilted = Query {
        incident: Incident::Plane { amplitude: Cx::new(1.0, 0.0), sines },
        point: [sines[0] * 0.1, 0.0, 0.1],
        ..qr
    };
    assert!((geometric_irradiance(&tilted) - 1.0).abs() < 1e-15);
    let missed = Query {
        point: [sines[0] * 0.1 + 1.5e-3, 0.0, 0.1],
        ..tilted
    };
    assert!(geometric_irradiance(&missed).abs() < 1e-18);
    let _ = PI;
}
