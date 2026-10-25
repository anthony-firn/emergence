//! Fabrics against the published geometric models: Peirce's woven
//! relations, Munden's knitted relations, and exact yarn arithmetic.

use std::f64::consts::PI;

use cad::kits::textile::{
    fibre_density, peirce_plain, Knit, Pattern, Relaxation, Weave, Yarn,
};
use cad::meshing::mesh_solid;

fn cotton(tex: f64) -> Yarn {
    Yarn::new(tex, fibre_density::COTTON)
}

/// A shirting: 30 ends and picks per centimetre of 30 tex cotton, with the
/// yarns flattened as they are in a real woven fabric.
fn shirting() -> Weave {
    let p = 1.0 / 3000.0;
    let yarn = cotton(30.0).flattened(1.3);
    Weave::new(Pattern::Plain, yarn.clone(), yarn, p, p)
}

#[test]
fn yarn_diameter_from_count_is_the_standard_relation() {
    // d = sqrt(4 T / (pi rho phi)): a 30 tex cotton yarn at 0.6 packing.
    let y = cotton(30.0);
    let d = y.diameter();
    let expect = (4.0 * 30e-6 / (PI * 1540.0 * 0.6)).sqrt();
    assert!((d - expect).abs() / expect < 1e-15);
    // About 0.2 mm, as a 30 tex cotton yarn measures.
    assert!(d > 1.8e-4 && d < 2.1e-4, "{:.4} mm", d * 1e3);

    // Flattening preserves the yarn's area exactly.
    let flat = cotton(30.0).flattened(2.0);
    assert!((flat.area() - y.area()).abs() / y.area() < 1e-15);
    let (w, t) = flat.section();
    assert!((w / t - 4.0).abs() < 1e-12, "aspect ratio {}", w / t);

    // Fibre area is the count over the fibre density, exactly.
    assert!((y.fibre_area() - 30e-6 / 1540.0).abs() / y.fibre_area() < 1e-15);
    // And the packing factor is the ratio of fibre area to yarn area.
    assert!((y.fibre_area() / y.area() - 0.6).abs() < 1e-12);
}

#[test]
fn weave_patterns_interlace_correctly() {
    // Plain weave alternates every end and every pick.
    let plain = Pattern::Plain;
    assert_eq!(plain.repeat_size(), (2, 2));
    for i in 0..4 {
        for j in 0..4 {
            assert_eq!(plain.over(i, j), (i + j) % 2 == 0);
        }
    }
    assert!((plain.mean_float() - 1.0).abs() < 1e-12, "plain float {}", plain.mean_float());

    // A 2/2 twill: each end floats over two picks, under two, advancing
    // by one each end.
    let twill = Pattern::Twill { over: 2, under: 2, step: 1 };
    assert_eq!(twill.repeat_size(), (4, 4));
    assert!((twill.mean_float() - 2.0).abs() < 1e-12, "2/2 twill float {}", twill.mean_float());
    // The diagonal: end i+1's interlacing is end i's shifted by one pick.
    for i in 0..8 {
        for j in 0..4 {
            assert_eq!(
                twill.over(i + 1, j),
                twill.over(i, (j + 4 - 1) % 4),
                "twill diagonal at end {i}, pick {j}"
            );
        }
    }

    // A 5-end satin with step 2: one interlacing per end per repeat.
    let satin = Pattern::Satin { repeat: 5, step: 2 };
    for i in 0..5 {
        assert_eq!((0..5).filter(|&j| satin.over(i, j)).count(), 1);
    }
    assert!((satin.mean_float() - 2.5).abs() < 1e-12, "satin float {}", satin.mean_float());
}

#[test]
fn plain_weave_thickness_and_cover_match_peirce() {
    let w = shirting();
    let p = w.warp_spacing;
    let (yw, yt) = w.warp.section();
    // Peirce's constraint h1 + h2 = d1 + d2 sets the thickness, in terms
    // of the yarns' thickness in the fabric.
    assert!(
        (w.thickness() - 2.0 * yt).abs() / yt < 1e-15,
        "thickness {:e} vs t1 + t2 {:e}",
        w.thickness(),
        2.0 * yt
    );

    // Cover factor K = width/p, and the fabric cover K1 + K2 - K1 K2.
    let (k1, k2, cover) = w.cover();
    assert!((k1 - yw / p).abs() < 1e-18 && (k2 - yw / p).abs() < 1e-18);
    assert!((cover - (k1 + k2 - k1 * k2)).abs() < 1e-15);
    assert!(cover > 0.0 && cover < 1.0, "cover {cover:.4}");
    assert!(w.interference() < 0.0, "this fabric should be open");

    // The flattening is inside the range the geometry allows.
    let (fmin, fmax) = w.consistent_flattening();
    assert!(
        w.warp.flattening > fmin && w.warp.flattening < fmax,
        "flattening {} outside the consistent range ({fmin:.3}, {fmax:.3})",
        w.warp.flattening
    );

    // The generated solid's thickness must match: the swept yarns span
    // exactly t1 + t2 in z.
    let solid = w.clone().with_repeats([1, 1]).solid();
    let (lo, hi) = solid.bbox();
    let measured = hi[2] - lo[2];
    assert!(
        (measured - w.thickness()).abs() / w.thickness() < 1e-9,
        "measured thickness {measured:e} vs {:e}",
        w.thickness()
    );

    // And the yarns sweep cleanly.
    w.check().expect("open plain weave should be geometrically consistent");
}

#[test]
fn a_jammed_weave_is_refused_not_drawn() {
    // 80 ends/cm of a yarn 0.2 mm across cannot fit: K1 = 1.6.
    let p = 1.0 / 8000.0;
    let w = Weave::new(Pattern::Plain, cotton(30.0), cotton(30.0), p, p);
    assert!(w.interference() > 0.0);
    let err = w.check().unwrap_err();
    assert!(err.contains("interfere"), "{err}");

    // A round yarn set at 40 ends/cm fits side by side but its centreline
    // then curves more tightly than the yarn is thick, so sweeping it
    // would self-intersect. No flattening fixes it: the allowed range is
    // empty, because p is below d sqrt(pi/2).
    let tight = Weave::new(Pattern::Plain, cotton(30.0), cotton(30.0), 1.0 / 4000.0, 1.0 / 4000.0);
    assert!(tight.interference() < 0.0, "this setting does fit side by side");
    let (fmin, fmax) = tight.consistent_flattening();
    assert!(fmin > fmax, "range should be empty: ({fmin:.3}, {fmax:.3})");
    let err = tight.check().unwrap_err();
    assert!(err.contains("self-intersects"), "{err}");
}

#[test]
fn crimp_and_swept_volume_are_consistent() {
    let w = shirting().with_repeats([1, 1]);
    let p = w.warp_spacing;
    let (c1, c2) = w.crimp();
    // A symmetric plain weave crimps both sets equally.
    assert!((c1 - c2).abs() < 1e-9, "warp crimp {c1:e}, weft {c2:e}");
    // Positive, and in the range a woven fabric's crimp occupies
    // (a few per cent to about 15 %).
    assert!(c1 > 0.0 && c1 < 0.2, "crimp {c1:.5}");

    // Against the exact arc length of the sinusoidal centreline: the path
    // is a raised cosine of amplitude a over a half period p, so its
    // length per crossing is the elliptic integral evaluated here by
    // high-order quadrature, independent of the generated polyline.
    let amp = 0.5 * w.crimp_height();
    let n = 200_000;
    let mut exact = 0.0;
    for i in 0..n {
        let t = (i as f64 + 0.5) / n as f64;
        // z(x) = amp * cos(pi x / p) over one crossing, so
        // dz/dx = -amp pi/p sin(pi x/p).
        let dz = -amp * PI / p * (PI * t).sin();
        exact += (1.0 + dz * dz).sqrt() * p / n as f64;
    }
    let exact_crimp = exact / p - 1.0;
    let rel = (c1 - exact_crimp).abs() / exact_crimp;
    println!("crimp: generated {c1:.6}, exact arc length {exact_crimp:.6} ({rel:.3e})");
    assert!(rel < 0.02, "generated crimp differs from the exact arc length by {rel:.3e}");

    // The swept yarn volume must equal section area times path length.
    let (mesh, check) = w.yarns().into_iter().next().unwrap();
    assert!(check.watertight);
    let path_length = w.warp_path(0, 2).length();
    // The swept volume is the *polygonal* section's area times the path
    // length; the section is an inscribed n-gon, whose area is below the
    // ellipse's by a known factor.
    let (yt, yw) = (w.warp.section().1, w.warp.section().0);
    let polygon_area = cad::path::Profile::ellipse(0.5 * yt, 0.5 * yw, w.sides).area().abs();
    let expect = polygon_area * path_length;
    let rel = (check.volume - expect).abs() / expect;
    let inscribed = w.sides as f64 / 2.0 * (2.0 * PI / w.sides as f64).sin() / PI;
    println!(
        "yarn volume: swept {:e}, polygon area x length {expect:e} ({rel:.3e}); \
         the {}-gon holds {:.4} of the ellipse",
        check.volume, w.sides, inscribed
    );
    assert!(rel < 2e-3, "swept volume differs from area x length by {rel:.3e}");
    // And the polygon really is short of the ellipse by that factor.
    let ellipse_rel = (polygon_area - w.warp.area()).abs() / w.warp.area();
    assert!((ellipse_rel - (1.0 - inscribed)).abs() < 1e-12);
    assert!(mesh.report().watertight);
}

#[test]
fn areal_density_matches_the_textile_calculation() {
    // 24 ends/cm and 20 picks/cm of 30 tex cotton: the standard sum
    // (ends/m * tex * (1+crimp) + picks/m * tex * (1+crimp)) / 1e6 kg/m^2.
    let yarn = cotton(30.0).flattened(1.3);
    let w = Weave::new(Pattern::Plain, yarn.clone(), yarn, 1.0 / 2400.0, 1.0 / 2000.0);
    let (c1, c2) = w.crimp();
    let expect = (2400.0 * 30e-6 * (1.0 + c1)) + (2000.0 * 30e-6 * (1.0 + c2));
    assert!((w.areal_density() - expect).abs() / expect < 1e-15);
    // About 135 g/m^2, a shirting weight.
    let gsm = w.areal_density() * 1e3;
    assert!(gsm > 120.0 && gsm < 150.0, "{gsm:.1} g/m^2");

    // Fibre volume fraction is well under the packing limit for an open
    // fabric, and rises when the fabric is set more closely.
    let vf = w.fibre_volume_fraction();
    assert!(vf > 0.05 && vf < 0.6, "fibre volume fraction {vf:.4}");
    let yarn = cotton(30.0).flattened(1.3);
    let tight = Weave::new(Pattern::Plain, yarn.clone(), yarn, 1.0 / 3000.0, 1.0 / 3000.0);
    assert!(tight.fibre_volume_fraction() > vf);
}

#[test]
fn peirce_plain_weave_solves_its_own_relations() {
    let (d1, d2) = (0.2e-3, 0.2e-3);
    let (p1, p2) = (0.42e-3, 0.42e-3);
    let (th1, th2, h1, h2, l1, l2) = peirce_plain(d1, d2, p1, p2).expect("solution");
    let dd = d1 + d2;
    // The constraint that fixes the fabric thickness.
    assert!((h1 + h2 - dd).abs() / dd < 1e-12, "h1 + h2 = {:e}, D = {dd:e}", h1 + h2);
    // Each direction must satisfy both of Peirce's relations.
    for (p, h, th, l) in [(p2, h1, th1, l1), (p1, h2, th2, l2)] {
        let (s, c) = th.sin_cos();
        let r1 = ((l - dd * th) * c + dd * s - p).abs() / p;
        let r2 = ((l - dd * th) * s + dd * (1.0 - c) - h).abs() / h;
        assert!(r1 < 1e-9, "spacing relation residual {r1:e}");
        assert!(r2 < 1e-9, "crimp-height relation residual {r2:e}");
    }
    // A symmetric fabric splits the thickness equally.
    assert!((h1 - h2).abs() / dd < 1e-6, "h1 = {h1:e}, h2 = {h2:e}");
    assert!(th1 > 0.0 && th1 < 0.5 * PI);
    // Crimp from Peirce: l/p - 1.
    let crimp = l1 / p2 - 1.0;
    assert!(crimp > 0.0 && crimp < 0.5, "Peirce crimp {crimp:.4}");

    // A jammed fabric has no solution, and says so.
    assert!(peirce_plain(d1, d2, 0.15e-3, 0.15e-3).is_err());
}

#[test]
fn knit_geometry_matches_mundens_relations() {
    // Munden's constants, as published.
    assert_eq!(Relaxation::Dry.constants(), (5.0, 3.8));
    assert_eq!(Relaxation::Wet.constants(), (5.3, 4.1));
    assert!((Relaxation::Dry.stitch_constant() - 19.0).abs() < 1e-12);
    assert!((Relaxation::Wet.stitch_constant() - 21.73).abs() < 0.01);

    let yarn = Yarn::new(25.0, fibre_density::WOOL);
    let l = 4e-3;
    for state in [Relaxation::Dry, Relaxation::Wet] {
        let k = Knit::relaxed(yarn.clone(), l, state);
        let (kc, kw) = state.constants();
        // Courses and wales per metre follow Munden's law.
        assert!((1.0 / k.course_spacing - kc / l).abs() / (kc / l) < 1e-15);
        assert!((1.0 / k.wale_spacing - kw / l).abs() / (kw / l) < 1e-15);
        // Stitch density is Ks / l^2.
        let expect = state.stitch_constant() / (l * l);
        assert!((k.stitch_density() - expect).abs() / expect < 1e-12);

        // The amplitude is solved against the continuous arc length, so
        // the generated polyline is short by its chord error, which must
        // fall as the stations are refined.
        let mut errors = Vec::new();
        for stations in [48usize, 192] {
            let fine = Knit { stations, ..k.clone() };
            let measured = fine.measured_loop_length();
            errors.push((measured - l).abs() / l);
        }
        println!(
            "{state:?}: loop length {l:e}; polyline error {:.2e} at 48 stations, {:.2e} at 192; \
             loop reaches {:.2} course spacings",
            errors[0],
            errors[1],
            k.reach_in_courses()
        );
        let order = (errors[0] / errors[1]).log2() / (192f64 / 48.0).log2();
        // Second order, as the chords on the head and sinker arcs
        // require: with m chords on an arc the shortfall is theta^2/24,
        // and the arcs are 40 % of the stitch.
        assert!(order > 1.8, "chord error order {order:.2}");
        assert!(errors[1] < 2e-4, "refined loop length off by {:e}", errors[1]);
        // A plain-knit loop interlooks with the course below, so it must
        // reach more than one course spacing but not many.
        assert!(
            k.reach_in_courses() > 1.0 && k.reach_in_courses() < 2.5,
            "loop reaches {:.3} course spacings",
            k.reach_in_courses()
        );
        // The legs are solved against the three-dimensional length, so
        // the planar length falls short of the loop length by the
        // out-of-plane excursion — half a per cent for this yarn.
        let legs = k.leg_length().unwrap();
        let planar = k.planar_length(legs);
        assert!((k.spatial_length(legs) - l).abs() / l < 1e-9, "3D stitch length must be l");
        let out_of_plane = (l - planar) / l;
        assert!(
            out_of_plane > 0.0 && out_of_plane < 0.03,
            "out-of-plane share {out_of_plane:.4}"
        );

        // Areal density from loop length and stitch density.
        let gsm = k.areal_density() * 1e3;
        assert!(gsm > 50.0 && gsm < 400.0, "{state:?}: {gsm:.1} g/m^2");
        k.check().expect("a relaxed loop should sweep cleanly");
    }

    // A loop too tight for its yarn is refused.
    let fat = Yarn::new(400.0, fibre_density::WOOL);
    let tight = Knit::relaxed(fat, 1.5e-3, Relaxation::Dry);
    let err = tight.check().unwrap_err();
    assert!(err.contains("self-intersects") || err.contains("shorter than"), "{err}");
}

#[test]
fn fabric_solids_mesh_and_feed_a_solver() {
    let p = 1.0 / 3000.0;
    let yarn = cotton(30.0).flattened(1.3);
    let w = Weave::new(Pattern::Twill { over: 2, under: 2, step: 1 }, yarn.clone(), yarn, p, p)
        .with_repeats([1, 1]);
    w.check().expect("consistent");
    let solid = w.solid();
    // A window of the fabric meshes watertight, which is what the flow and
    // thermal solvers need.
    let (lo, hi) = solid.bbox();
    let window = solid.clipped(
        [lo[0] + 0.3 * (hi[0] - lo[0]), lo[1] + 0.3 * (hi[1] - lo[1]), lo[2] - 1e-6],
        [lo[0] + 0.7 * (hi[0] - lo[0]), lo[1] + 0.7 * (hi[1] - lo[1]), hi[2] + 1e-6],
    );
    let d = w.warp.diameter();
    let mesh = mesh_solid(&window, d / 6.0).expect("mesh");
    let report = mesh.report();
    // Peirce's constraint makes warp and weft exactly tangent at every
    // crossing, and the boundary of two tangent solids is non-manifold
    // there, so the mesh is closed but need not be manifold. Closed is
    // what matters: the inside/outside test is still exact.
    assert_eq!(report.boundary_edges, 0, "the fabric surface must be closed: {report:?}");
    assert!(
        report.nonmanifold_edges < report.triangles / 100,
        "too many pinch points for tangent contacts alone: {report:?}"
    );
    assert!(report.volume > 0.0);

    // The knitted fabric likewise.
    let k = Knit::relaxed(Yarn::new(25.0, fibre_density::WOOL), 4e-3, Relaxation::Wet).with_repeats([2, 2]);
    let solid = k.solid();
    assert!(solid.sdf([0.0, 0.0, 0.0]).is_finite());
    let (swept, check) = k.loops().into_iter().next().unwrap();
    assert!(check.watertight, "a swept course must be watertight: {check:?}");
    let one_course = cad::Solid::from_trimesh(swept);
    let mesh = mesh_solid(&one_course, k.yarn.diameter() / 5.0).expect("course mesh");
    assert_eq!(mesh.report().boundary_edges, 0, "remeshed course must be closed");
}
