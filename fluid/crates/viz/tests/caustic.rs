//! Caustics as the renderer's accuracy test.
//!
//! A caustic's position is fixed by refraction and reflection exactly, so
//! it has closed-form references a diffuse render never does. Each test
//! here checks the whole specular transport — Snell, Fresnel, total
//! internal reflection, and the march that finds the interfaces — against
//! one of them.

use std::f64::consts::PI;

use cad::solid::Solid;
use viz::caustic::{
    airy, asymptotic_negative, asymptotic_positive, mirror_caustic_parameter,
    mirror_reflected_ray, fold, fold_geometric, fold_peak, fold_scale, mirror_axis_crossing,
    mirror_catacaustic, rainbow_angle, rainbow_deviation, rainbow_incidence, sphere_deviation,
    trace_to_angles, trace_to_plane, Angular, Emitter, PhotonOptions, Receiver,
};
use viz::scene::{Material, Object, Scene};

#[test]
fn the_rainbow_angle_comes_out_of_the_photon_trace() {
    // Descartes' construction, and the sharpest test of refraction there
    // is: parallel light through a sphere with one internal reflection
    // piles up at a least-deviation angle that has a closed form. The
    // photons must find it without being told.
    let n = 1.333;
    let exact_incidence = rainbow_incidence(n, 1).unwrap();
    let exact_deviation = rainbow_deviation(n, 1).unwrap();
    let exact_angle = rainbow_angle(n, 1).unwrap();
    println!(
        "primary rainbow for n = {n}: incidence {:.4}°, deviation {:.4}°, seen at {:.4}°",
        exact_incidence.to_degrees(),
        exact_deviation.to_degrees(),
        exact_angle.to_degrees()
    );
    // The published figure for water is about 42°.
    assert!(
        (exact_angle.to_degrees() - 42.0).abs() < 0.5,
        "primary rainbow at {:.3}°",
        exact_angle.to_degrees()
    );
    // It really is the stationary point of the deviation.
    for delta in [-0.05, 0.05] {
        let d = sphere_deviation(exact_incidence + delta, n, 1).unwrap();
        assert!(d > exact_deviation, "the rainbow is the *least* deviation");
    }

    // Now trace it. A water sphere, a collimated beam wider than it, and
    // a histogram of exit directions.
    let radius = 1.0;
    let scene = Scene::new()
        .with_material(Material::Dielectric { index: n })
        .with_object(Object::new("drop", Solid::sphere([0.0; 3], radius), 0));
    let direction = [0.0, 0.0, 1.0];
    let emitter = Emitter::Beam {
        origin: [0.0, 0.0, -4.0],
        direction,
        radius: 0.999 * radius,
        power: 1.0,
    };
    let mut angles = Angular::new(direction, [0.0, PI], 900);
    // Three specular events: in, one internal reflection, out. Without
    // this selection the two-event forward beam is far brighter and the
    // peak lands on it — which is also why a rainbow needs rain and not a
    // window pane.
    let options = PhotonOptions {
        select_events: Some(3),
        ..PhotonOptions::default()
    };
    let stats = trace_to_angles(&scene, &emitter, &mut angles, 400_000, 0x51ED, &options);
    println!(
        "  {} photons, {} deposited, flux accounted {:.6}",
        stats.emitted,
        stats.deposited,
        stats.accounted()
    );
    // A sum of 400 000 terms, so the bookkeeping is checked to rounding
    // and not beyond it.
    assert!(
        (stats.accounted() - 1.0).abs() < 1e-9,
        "flux must be accounted for: {}",
        stats.accounted()
    );

    // The brightest deviation, by a parabolic fit, against the closed
    // form. The peak is the caustic: a fold in the ray map.
    let measured = angles.peak_deviation().expect("a peak");
    let error = (measured - exact_deviation).to_degrees().abs();
    println!(
        "  measured deviation {:.4}°, exact {:.4}°, error {error:.4}°",
        measured.to_degrees(),
        exact_deviation.to_degrees()
    );
    assert!(error < 0.25, "rainbow deviation off by {error:.4}°");
}

#[test]
fn the_secondary_rainbow_is_there_too_and_in_the_right_place() {
    // Two internal reflections, dimmer and reversed: the same closed form
    // with k = 2, and a test that the photon walk follows a longer
    // specular chain correctly.
    let n = 1.333;
    let primary = rainbow_angle(n, 1).unwrap().to_degrees();
    let secondary = rainbow_angle(n, 2).unwrap().to_degrees();
    println!("primary {primary:.3}°, secondary {secondary:.3}°");
    assert!((secondary - 51.0).abs() < 1.0, "secondary at {secondary:.3}°");
    assert!(secondary > primary, "the secondary is the wider bow");

    let radius = 1.0;
    let scene = Scene::new()
        .with_material(Material::Dielectric { index: n })
        .with_object(Object::new("drop", Solid::sphere([0.0; 3], radius), 0));
    let direction = [0.0, 0.0, 1.0];
    let emitter = Emitter::Beam {
        origin: [0.0, 0.0, -4.0],
        direction,
        radius: 0.999 * radius,
        power: 1.0,
    };
    // The secondary sits between 120° and 135° of deviation; binning only
    // that window keeps the much brighter primary out of the fit.
    let window = [125.0_f64.to_radians(), 135.0_f64.to_radians()];
    let mut angles = Angular::new(direction, window, 400);
    // Four events: in, two internal reflections, out.
    let options = PhotonOptions {
        select_events: Some(4),
        ..PhotonOptions::default()
    };
    let stats = trace_to_angles(&scene, &emitter, &mut angles, 600_000, 0x5EC0, &options);
    assert!((stats.accounted() - 1.0).abs() < 1e-9, "{}", stats.accounted());
    let measured = angles.peak_deviation().expect("a peak");
    let exact = rainbow_deviation(n, 2).unwrap();
    // The secondary's deviation exceeds pi, so the angle *seen* is
    // D - pi; the histogram bins the unsigned deviation, which folds it.
    let folded = 2.0 * PI - exact;
    let error = (measured - folded).to_degrees().abs();
    println!(
        "  secondary: measured {:.4}°, folded exact {:.4}°, error {error:.4}°",
        measured.to_degrees(),
        folded.to_degrees()
    );
    assert!(error < 0.3, "secondary deviation off by {error:.4}°");
}

#[test]
fn dispersion_spreads_the_bow_by_the_published_amount() {
    // The width of a real rainbow is the dispersion of water across the
    // visible: the red bow sits outside the violet one. Using
    // `waveoptics`' water dispersion (Daimon & Masumura 2007), the
    // closed form gives the published spread of about 1.8°.
    let red = waveoptics::materials::water(656.3e-9);
    let violet = waveoptics::materials::water(404.7e-9);
    println!("water index: {red:.5} at 656 nm, {violet:.5} at 405 nm");
    assert!(violet > red, "water is normally dispersive");
    let a_red = rainbow_angle(red, 1).unwrap().to_degrees();
    let a_violet = rainbow_angle(violet, 1).unwrap().to_degrees();
    let spread = a_red - a_violet;
    println!("primary bow: red at {a_red:.3}°, violet at {a_violet:.3}°, width {spread:.3}°");
    // Published: the primary runs from about 40.6° (violet) to 42.4°
    // (red), a width near 1.8°.
    assert!((a_red - 42.4).abs() < 0.3, "red bow at {a_red:.3}°");
    assert!((a_violet - 40.6).abs() < 0.4, "violet bow at {a_violet:.3}°");
    assert!((spread - 1.8).abs() < 0.4, "bow width {spread:.3}°");

    // And a dispersive glass sphere puts the colours in different places,
    // which is what makes a traced caustic coloured rather than white.
    let glass = Material::bk7();
    let n_red = glass.refractive_index(656.3e-9).unwrap();
    let n_blue = glass.refractive_index(486.1e-9).unwrap();
    assert!(n_blue > n_red);
    // The Abbe number of N-BK7 is about 64.
    let n_d = glass.refractive_index(587.6e-9).unwrap();
    let abbe = (n_d - 1.0) / (n_blue - n_red);
    println!("N-BK7: nd = {n_d:.5}, Abbe {abbe:.1}");
    assert!((abbe - 64.17).abs() < 0.5, "Abbe number {abbe:.2}");
}

#[test]
fn a_spherical_mirrors_caustic_has_its_cusp_at_the_paraxial_focus() {
    // The catacaustic of a circle in parallel light is a curve with a cusp
    // at the paraxial focus, R/2 from the vertex. Rays near the axis cross
    // there; marginal rays cross short of it, and the locus of crossings
    // is the caustic.
    let r = 1.0;
    let paraxial = 0.5 * r;
    let near_axis = mirror_axis_crossing(r, 1e-6).unwrap();
    assert!(
        (near_axis - paraxial).abs() < 1e-10,
        "paraxial crossing {near_axis:.9}, expected {paraxial}"
    );
    for h in [0.1, 0.3, 0.6] {
        let z = mirror_axis_crossing(r, h).unwrap();
        assert!(z < paraxial, "a marginal ray at h = {h} must cross short of R/2");
    }
    let cusp = mirror_catacaustic(r, 0.0);
    println!("catacaustic cusp at {cusp:?}, paraxial focus at z = {paraxial}");
    assert!(cusp[0].abs() < 1e-12 && (cusp[1] - paraxial).abs() < 1e-12);

    // The envelope, checked against what an envelope is: the limit of the
    // intersection of neighbouring rays. No tracing, no tolerance to
    // choose but the step.
    for t in [0.05, 0.2, 0.5, 0.9, 1.2] {
        let (p0, d0) = mirror_reflected_ray(r, t);
        let eps = 1e-6;
        let (p1, d1) = mirror_reflected_ray(r, t + eps);
        // Solve p0 + a d0 = p1 + b d1, i.e. [d0, -d1] (a, b)^T = p1 - p0.
        let det = d1[0] * d0[1] - d1[1] * d0[0];
        let a = ((p1[1] - p0[1]) * d1[0] - (p1[0] - p0[0]) * d1[1]) / det;
        let meet = [p0[0] + a * d0[0], p0[1] + a * d0[1]];
        let env = mirror_catacaustic(r, t);
        let err = ((meet[0] - env[0]).powi(2) + (meet[1] - env[1]).powi(2)).sqrt();
        println!(
            "t = {t}: neighbours meet at [{:.9}, {:.9}], envelope [{:.9}, {:.9}], {err:.2e} apart",
            meet[0], meet[1], env[0], env[1]
        );
        // First order in the step, so the agreement is first order too.
        assert!(err < 20.0 * eps, "the envelope formula is wrong at t = {t}: {err:.3e}");
        // And the distance along the ray matches the closed form.
        assert!((a - mirror_caustic_parameter(r, t)).abs() < 20.0 * eps);
    }

    // Traced: a perfect mirror and a collimated beam. A ball of glass
    // always presents a *convex* face, so a concave mirror is the
    // complement — a slab with the ball cut out of it. The vertex is at
    // the origin, the centre of curvature at z = R, and the light arrives
    // along -z, which is the frame both closed forms above are written in.
    let bowl = Solid::cuboid([0.0, 0.0, 0.2 * r], [1.3 * r, 1.3 * r, 0.4 * r])
        .minus(Solid::sphere([0.0, 0.0, r], r));
    let scene = Scene::new()
        .with_material(Material::Mirror)
        .with_object(Object::new("mirror", bowl, 0));
    let options = PhotonOptions {
        max_bounces: 4,
        select_events: Some(1),
        ..PhotonOptions::default()
    };

    // The cusp is where the caustic *ends*, so it is the narrow-beam
    // limit, not the brightest plane of a wide beam. A beam out to h = a
    // spreads its crossings over R/2 down to R - R/(2 cos asin(a/R)),
    // which is about a^2/(4R); the peak of a wide beam therefore sits
    // short of the cusp by that much, and the test is that it converges
    // like a^2.
    let mut errors = Vec::new();
    for &a in &[0.4 * r, 0.2 * r, 0.1 * r] {
        // An aperture small against the beam, or every plane catches
        // everything and there is no peak to find.
        let half = 0.05 * a;
        let mut best = (0.0, f64::NAN);
        let steps = 60;
        for i in 0..=steps {
            let z = paraxial - 0.12 * r + 0.14 * r * i as f64 / steps as f64;
            let mut axial = Receiver::new(
                [-half, -half, z],
                [2.0 * half, 0.0, 0.0],
                [0.0, 2.0 * half, 0.0],
                [8, 8],
            );
            let emitter = Emitter::Beam {
                origin: [0.0, 0.0, 3.0],
                direction: [0.0, 0.0, -1.0],
                radius: a,
                power: 1.0,
            };
            let stats =
                trace_to_plane(&scene, &emitter, &mut axial, 20_000, 0xCA05 + i as u64, &options);
            assert!((stats.accounted() - 1.0).abs() < 1e-9, "{}", stats.accounted());
            let flux = axial.total_flux();
            if flux > best.0 {
                best = (flux, z);
            }
        }
        let spread = a * a / (4.0 * r);
        let error = (best.1 - paraxial).abs();
        println!(
            "beam h <= {a:.3}: peak flux at z = {:.4}, cusp {paraxial}, off by {error:.4}              (caustic spreads {spread:.4})",
            best.1
        );
        // The peak must sit *short* of the cusp, inside the caustic.
        assert!(best.1 <= paraxial + 1e-9, "the peak overshot the cusp: {:.4}", best.1);
        // And within the spread of the caustic, or the scan step, whichever
        // is coarser: a 0.14 m scan in 60 steps cannot resolve better than
        // 2.4 mm.
        let floor = 0.14 * r / steps as f64;
        assert!(
            error < 2.0 * spread.max(floor),
            "the peak is {error:.4} from the cusp, caustic spread only {spread:.4}"
        );
        errors.push(error);
    }
    // Narrowing the beam walks the peak onto the cusp.
    println!("peak-to-cusp distance as the beam narrows: {errors:?}");
    assert!(
        errors[2] < errors[0],
        "a narrower beam must land closer to the cusp: {errors:?}"
    );
    assert!(errors[2] < 0.01 * r, "the narrowest beam missed the cusp: {:.4}", errors[2]);
}

#[test]
fn a_paraboloid_has_no_caustic_at_all() {
    // The one specular surface with no caustic: every ray parallel to the
    // axis passes through one point, so `marginal_focus` is R/2 at every
    // height. That part is exact.
    let r = 1.0;
    let focus = 0.5 * r;
    let exact = cad::kits::optic::Mirror::paraboloid(r, 0.5, 0.05);
    for h in [1e-6, 0.05, 0.2, 0.35, 0.49] {
        let f = exact.marginal_focus(h).expect("a real focus");
        assert!(
            (f - focus).abs() < 1e-12,
            "a paraboloid has no spherical aberration: {f:.15} at h = {h}"
        );
    }

    // Traced, it does make a spot — but the spot is the *mesh*, not the
    // optics. `Mirror::meridian` revolves a polygon of `points` segments,
    // so each facet is a flat mirror whose normal is wrong by about
    // `dz/dr` across half a facet, and the focus it makes is displaced by
    // twice that times the focal length. The displacement is therefore
    // proportional to the facet width, 1/points.
    //
    // So the test is the scaling, which is the statement that the error
    // is discretisation: shrink the aperture in step with 1/points and
    // the caught fraction should hold still.
    let mut caught = Vec::new();
    for points in [40usize, 80, 160] {
        let mut dish = cad::kits::optic::Mirror::paraboloid(r, 0.5, 0.05);
        dish.points = points;
        let scene = Scene::new()
            .with_material(Material::Mirror)
            .with_object(Object::new("dish", dish.solid().expect("paraboloid"), 0));
        let emitter = Emitter::Beam {
            origin: [0.0, 0.0, 3.0],
            direction: [0.0, 0.0, -1.0],
            radius: 0.45,
            power: 1.0,
        };
        let options = PhotonOptions {
            max_bounces: 3,
            select_events: Some(1),
            ..PhotonOptions::default()
        };
        // Half-width proportional to the facet width.
        let half = 32.0 / points as f64 * 0.025;
        let mut window = Receiver::new(
            [-half, -half, focus],
            [2.0 * half, 0.0, 0.0],
            [0.0, 2.0 * half, 0.0],
            [4, 4],
        );
        let stats = trace_to_plane(&scene, &emitter, &mut window, 40_000, 0xABCD, &options);
        let fraction = window.total_flux() / stats.flux_emitted;
        println!(
            "points = {points:4}: {:.4} of the flux through a {:.2} mm aperture",
            fraction,
            2000.0 * half
        );
        caught.push(fraction);
    }
    // The fraction holds still as both the facets and the aperture shrink
    // together: the spot is the facet width.
    let (lo, hi) = (
        caught.iter().cloned().fold(f64::INFINITY, f64::min),
        caught.iter().cloned().fold(0.0, f64::max),
    );
    assert!(
        hi - lo < 0.01,
        "the spot does not scale with the facet width: {caught:?}"
    );
    assert!(lo > 0.95, "too little flux caught at all: {caught:?}");

    // And with the aperture left alone, more facets means more flux
    // through it — the limit is a point.
    let mut better = Vec::new();
    for points in [40usize, 160, 640] {
        let mut dish = cad::kits::optic::Mirror::paraboloid(r, 0.5, 0.05);
        dish.points = points;
        let scene = Scene::new()
            .with_material(Material::Mirror)
            .with_object(Object::new("dish", dish.solid().expect("paraboloid"), 0));
        let emitter = Emitter::Beam {
            origin: [0.0, 0.0, 3.0],
            direction: [0.0, 0.0, -1.0],
            radius: 0.45,
            power: 1.0,
        };
        let options = PhotonOptions {
            max_bounces: 3,
            select_events: Some(1),
            ..PhotonOptions::default()
        };
        let half = 0.0025;
        let mut window = Receiver::new(
            [-half, -half, focus],
            [2.0 * half, 0.0, 0.0],
            [0.0, 2.0 * half, 0.0],
            [4, 4],
        );
        let stats = trace_to_plane(&scene, &emitter, &mut window, 20_000, 0xBEE5, &options);
        let fraction = window.total_flux() / stats.flux_emitted;
        println!("points = {points:4}: {fraction:.4} through a fixed 5 mm aperture");
        better.push(fraction);
    }
    assert!(
        better[1] > better[0] + 0.1 && better[2] > better[1] + 0.1,
        "refining the mesh must tighten the focus: {better:?}"
    );
    assert!(better[2] > 0.9, "640 facets should focus to 5 mm: {:.4}", better[2]);
}

#[test]
fn flux_is_conserved_through_a_glass_sphere() {
    // The bookkeeping test: every photon's flux is deposited, spilled or
    // lost, and the three must sum to what was emitted. Without this a
    // caustic can look right and be the wrong brightness.
    let scene = Scene::new()
        .with_material(Material::Dielectric { index: 1.5 })
        .with_object(Object::new("ball", Solid::sphere([0.0; 3], 0.5), 0));
    let emitter = Emitter::Beam {
        origin: [0.0, 0.0, -2.0],
        direction: [0.0, 0.0, 1.0],
        radius: 0.45,
        power: 3.7,
    };
    // At the paraxial focus of a ball lens, f = nR/(2(n-1)) measured from
    // the centre (Born & Wolf 4.4): 0.75 m for n = 1.5, R = 0.5 m. Farther
    // out the light has already diverged and there is no caustic to see.
    let mut receiver = Receiver::new(
        [-0.1, -0.1, 0.75],
        [0.2, 0.0, 0.0],
        [0.0, 0.2, 0.0],
        [64, 64],
    );
    let stats = trace_to_plane(
        &scene,
        &emitter,
        &mut receiver,
        100_000,
        0xF10C,
        &PhotonOptions::default(),
    );
    println!(
        "emitted {:.4} W, deposited {:.4}, spilled {:.4}, lost {:.4}; accounted {:.12}",
        stats.flux_emitted, stats.flux_deposited, stats.flux_spilled, stats.flux_lost,
        stats.accounted()
    );
    // 100 000 photons of 3.7e-5 W each: the sum carries the beam's power
    // to within the rounding of accumulating that many terms.
    assert!(
        (stats.flux_emitted - 3.7).abs() < 1e-11,
        "the beam carries its power: {}",
        stats.flux_emitted
    );
    assert!(
        (stats.accounted() - 1.0).abs() < 1e-9,
        "flux is not accounted for: {}",
        stats.accounted()
    );
    // The bins hold exactly the flux that landed on them: the tent splat
    // conserves, so this is a statement about the estimator too.
    assert!(
        (receiver.total_flux() - stats.flux_deposited).abs() < 1e-11,
        "the bins hold {:.9} but {:.9} was deposited",
        receiver.total_flux(),
        stats.flux_deposited
    );
    // A caustic is there: the brightest bin is far above the mean.
    let irradiance = receiver.irradiance();
    let mean = irradiance.iter().sum::<f64>() / irradiance.len() as f64;
    let (peak, at) = receiver.peak();
    println!("peak irradiance {peak:.4} W/m^2 at bin {at:?}, mean over the plate {mean:.4}");
    assert!(peak > 10.0 * mean, "a sphere should make a bright caustic: {:.1}x", peak / mean);
    // And it is on the axis, where the paraxial focus is.
    assert!(
        at[0].abs_diff(32) <= 1 && at[1].abs_diff(32) <= 1,
        "the caustic should be on the axis, found bin {at:?}"
    );
}

#[test]
fn the_caustic_estimate_converges_as_one_over_root_n() {
    // Monte Carlo, so the noise must fall as 1/sqrt(N) and the answer must
    // not move. Checked on the flux through a fixed window, which has a
    // definite value.
    //
    // The statement is that `sd * sqrt(N)` is *constant*, which is what
    // 1/sqrt(N) means, and it needs enough runs to say: a sample standard
    // deviation from n runs is itself uncertain by about 1/sqrt(2(n-1)),
    // so eight runs give 27% on each sd and 38% on a ratio of two. An
    // earlier version of this test used eight and reported the noise
    // falling 7.9x for 16x the photons where 4x is predicted — which was
    // the error bar on the estimate, not a violation of 1/sqrt(N).
    let scene = Scene::new()
        .with_material(Material::Dielectric { index: 1.5 })
        .with_object(Object::new("ball", Solid::sphere([0.0; 3], 0.5), 0));
    let emitter = Emitter::Beam {
        origin: [0.0, 0.0, -2.0],
        direction: [0.0, 0.0, 1.0],
        radius: 0.45,
        power: 1.0,
    };
    let window = |count: usize, seed: u64| -> f64 {
        let mut r = Receiver::new(
            [-0.1, -0.1, 1.2],
            [0.2, 0.0, 0.0],
            [0.0, 0.2, 0.0],
            [16, 16],
        );
        trace_to_plane(&scene, &emitter, &mut r, count, seed, &PhotonOptions::default());
        r.total_flux()
    };
    const RUNS: u64 = 24;
    let mut scaled = Vec::new();
    let mut means = Vec::new();
    for &count in &[2_000usize, 8_000, 32_000] {
        let runs: Vec<f64> = (0..RUNS)
            // Well-separated seeds, so the runs are independent samples
            // and not neighbouring states of one stream.
            .map(|k| window(count, 0x9E37_79B9_7F4A_7C15u64.wrapping_mul(k + 1)))
            .collect();
        let mean = runs.iter().sum::<f64>() / runs.len() as f64;
        let var = runs.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (runs.len() - 1) as f64;
        let sd = var.sqrt();
        println!(
            "N = {count:6}: mean {mean:.6}, sd {sd:.6}, sd*sqrt(N) {:.4}",
            sd * (count as f64).sqrt()
        );
        scaled.push(sd * (count as f64).sqrt());
        means.push(mean);
    }
    // The mean does not drift.
    let drift = (means[2] - means[0]).abs() / means[0];
    println!("mean drifted {:.2}% from N = 2000 to N = 32000", 100.0 * drift);
    assert!(drift < 0.03, "the estimate drifted: {:?}", means);
    // And sd * sqrt(N) holds still. With 24 runs each sd carries about 15%,
    // so a ratio of two carries about 21%; a factor of 1.6 is a real
    // constraint and not a tautology.
    let lo = scaled.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = scaled.iter().cloned().fold(0.0, f64::max);
    println!("sd*sqrt(N) spans {lo:.4} to {hi:.4}, a factor of {:.2}", hi / lo);
    assert!(hi / lo < 1.6, "the noise does not fall as 1/sqrt(N): {scaled:?}");
}

#[test]
fn the_airy_function_matches_its_published_values() {
    // The wave caustic rests on this, so it is pinned to the tables
    // (Abramowitz & Stegun, table 10.11; DLMF 9.9).
    assert!((airy(0.0) - 0.355_028_053_887_817).abs() < 1e-13, "Ai(0) = {}", airy(0.0));
    assert!((airy(1.0) - 0.135_292_416_312_881).abs() < 1e-13, "Ai(1) = {}", airy(1.0));
    assert!((airy(-1.0) - 0.535_560_883_292_352).abs() < 1e-13, "Ai(-1) = {}", airy(-1.0));
    // Ai(5) = 1.0834442808e-4. My first reference figure here was wrong
    // in the fifth digit and the series caught it, not the other way
    // round — so it was re-derived in 60-digit arithmetic before being
    // believed.
    assert!((airy(5.0) - 1.083_444_280_849e-4).abs() < 1e-15, "Ai(5) = {:e}", airy(5.0));

    // The first three zeros.
    let zero = |near: f64| {
        let (mut lo, mut hi) = (near - 0.4, near + 0.4);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if airy(mid).signum() == airy(lo).signum() {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    };
    for (near, published) in [
        (-2.34, -2.338_107_410_459_767),
        (-4.09, -4.087_949_444_130_971),
        (-5.52, -5.520_559_828_095_551),
    ] {
        let z = zero(near);
        assert!((z - published).abs() < 1e-9, "zero near {near}: {z:.12} vs {published:.12}");
    }

    // The first maximum of Ai(-x), which sets the caustic's peak.
    let mut best = (0.0f64, 0.0f64);
    for i in 0..200_000 {
        let x = 2.0 * i as f64 / 200_000.0;
        let v = airy(-x);
        if v > best.1 {
            best = (x, v);
        }
    }
    println!("first maximum of Ai(-x): {:.9} at x = {:.9}", best.1, best.0);
    assert!((best.0 - 1.018_792_971_647).abs() < 1e-5);
    assert!((best.1 - 0.535_656_656_9).abs() < 1e-9);

    // Far from the origin, where the Taylor series is no longer the
    // branch that runs. These are the same series evaluated in 60-digit
    // arithmetic, not a table — an independent computation of the same
    // definition, which is what is needed to test the asymptotic.
    for (x, exact) in [
        (-6.0, -0.329_145_173_629_823_105),
        (-8.0, -0.052_705_050_356_386_203),
        (-10.0, 0.040_241_238_486_443_191),
        (-16.0, -0.143_057_931_669_099_698),
        (8.0, 4.692_207_616_099_232e-8),
        (10.0, 1.104_753_255_289_869e-10),
    ] {
        let got = airy(x);
        let rel = (got - exact).abs() / exact.abs();
        println!("Ai({x:5}) = {got:.15e}, exact {exact:.15e}, relative {rel:.1e}");
        assert!(rel < 2e-7, "Ai({x}) = {got:e} vs {exact:e}, relative {rel:.2e}");
    }

    // The series and the asymptotic agree where they meet. This is the
    // check that found the missing `sin` series in the oscillatory
    // asymptotic: the two branches differed by 3% at x = -8.
    // Both branches at the *same* argument. Nudging x across the
    // crossover instead measures Ai', which is O(1) here and swamps the
    // mismatch being looked for.
    for (x, asymptotic) in [(-7.0, asymptotic_negative(7.0)), (6.0, asymptotic_positive(6.0))] {
        let series = airy(x);
        let rel = (series - asymptotic).abs() / series.abs();
        println!("branches at {x}: series {series:e}, asymptotic {asymptotic:e}, relative {rel:.1e}");
        assert!(rel < 2e-7, "the branches disagree at {x}: {series:e} vs {asymptotic:e}");
    }
}

#[test]
fn the_wave_caustic_is_finite_where_rays_say_infinity() {
    // The sharpest statement of the ladder in this crate: at a caustic,
    // geometrical optics predicts infinite irradiance, which is the model
    // failing rather than a numerical difficulty. The exact answer is
    // Airy, and it is finite.
    let lambda = 550e-9;
    let rho = 0.1; // a caustic with a 100 mm radius of curvature
    let w = fold_scale(lambda, rho);
    println!("fold scale at 550 nm and rho = 0.1 m: {:.3} um", w * 1e6);
    assert!(w > 1e-6 && w < 1e-4, "{w:e}");

    // Geometrical optics diverges; the wave answer does not.
    assert!(fold_geometric(0.0).is_infinite());
    let peak = fold_peak(w);
    assert!(peak.is_finite() && peak > 0.0, "peak {peak:e}");
    println!("wave peak {peak:.4e}, geometric at one scale out {:.4e}", fold_geometric(w));
    // The peak is 2 pi Ai^2(-1.0188) / sqrt(w) = 1.8028 / sqrt(w).
    assert!(
        (peak * w.sqrt() - 1.802_8).abs() < 1e-3,
        "peak x sqrt(w) = {:.6}",
        peak * w.sqrt()
    );
    // On the dark side the wave answer is small but not zero: light
    // reaches where rays say it cannot.
    assert!(fold(-3.0 * w, w) > 0.0 && fold(-3.0 * w, w) < 0.05 * peak);
    assert!(fold_geometric(-3.0 * w) == 0.0);

    // Far on the lit side the two agree in the mean, which is what fixes
    // the constant: average over a few fringes.
    for &s in &[40.0, 90.0, 160.0] {
        let centre = s * w;
        let span = 8.0 * w;
        let n = 4000;
        let mut wave = 0.0;
        let mut geom = 0.0;
        for i in 0..n {
            let x = centre - 0.5 * span + span * (i as f64 + 0.5) / n as f64;
            wave += fold(x, w);
            geom += fold_geometric(x);
        }
        let ratio = wave / geom;
        println!("at {s} scales out, wave/geometric = {ratio:.4}");
        assert!((ratio - 1.0).abs() < 0.05, "they should agree in the mean: {ratio:.4}");
    }
}

#[test]
fn the_caustic_brightens_as_one_over_cube_root_of_wavelength() {
    // The two scalings that make a caustic bright but finite: the peak
    // irradiance goes as lambda^(-1/3) and the fringe spacing as
    // lambda^(2/3). Both follow from the Airy form and neither is
    // available to ray optics at all.
    let rho = 0.05;
    let mut peaks = Vec::new();
    for &lambda in &[400e-9, 800e-9] {
        let w = fold_scale(lambda, rho);
        peaks.push((lambda, w, fold_peak(w)));
    }
    let (l0, w0, p0) = peaks[0];
    let (l1, w1, p1) = peaks[1];
    // Doubling the wavelength widens the scale by 2^(2/3) = 1.587 ...
    let widened = w1 / w0;
    println!("doubling lambda: scale x {widened:.4} (2^(2/3) = {:.4})", 2f64.powf(2.0 / 3.0));
    assert!((widened - 2f64.powf(2.0 / 3.0)).abs() < 1e-9);
    // ... and dims the peak by 2^(1/3) = 1.26.
    let dimmed = p0 / p1;
    println!("  peak x {dimmed:.4} (2^(1/3) = {:.4})", 2f64.powf(1.0 / 3.0));
    assert!((dimmed - 2f64.powf(1.0 / 3.0)).abs() < 1e-6);
    assert!(l1 > l0);

    // And the indicator: a caustic narrower than a pixel is a line that
    // ray optics draws correctly; one wider is a fringed band it cannot.
    let pixel = 20e-6;
    let narrow = fold_scale(550e-9, 1e-3) / pixel;
    let wide = fold_scale(550e-9, 10.0) / pixel;
    println!("fold scale in pixels: {narrow:.3} at rho = 1 mm, {wide:.3} at rho = 10 m");
    assert!(narrow < 1.0 && wide > 1.0);
}
