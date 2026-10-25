//! Mechanical generators against published standards: ISO 724 thread
//! dimensions, the standard involute gear relations, and the octet-truss
//! relative density and stiffness.

use std::f64::consts::PI;

use cad::kits::lattice::{octet_relative_density, LatticeSpec, Topology, TpmsSpec};
use cad::kits::part::{SpurGear, Thread};
use cad::meshing::mesh_solid;
use cad::sketch::Region2;
use cad::solid::Solid;

/// Angular half-width of a region's boundary at radius `r`, measured by
/// bisection on the field itself (not from the formula that built it).
fn measured_half_angle(region: &Region2, r: f64) -> f64 {
    let at = |a: f64| region.sdf([r * a.cos(), r * a.sin()]);
    assert!(at(0.0) < 0.0, "radius {r} is not inside the tooth at angle 0");
    let mut hi = 0.0;
    let step = 1e-4;
    while at(hi) < 0.0 && hi < 1.0 {
        hi += step;
    }
    let mut lo = hi - step;
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if at(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn gear_tooth_thickness_at_the_pitch_circle() {
    // Measured on the generated polygon, which samples the pitch circle
    // exactly, so the only error is the bisection tolerance.
    // The standard circular tooth thickness is s = m (pi/2 + 2 x tan a),
    // which for an unshifted gear is exactly pi m / 2 (ISO 21771).
    for (m, z, shift) in [(2.0, 20usize, 0.0), (1.0, 37, 0.0), (3.0, 12, 0.3), (2.5, 48, -0.2)] {
        let mut gear = SpurGear::new(m, z, 5.0);
        gear.shift = shift;
        let region = gear.region();
        let rp = gear.pitch_radius();
        let psi = measured_half_angle(&region, rp);
        let measured = 2.0 * rp * psi;
        let expected = gear.report().pitch_thickness;
        let rel = (measured - expected).abs() / expected;
        assert!(
            rel < 1e-6,
            "m={m} z={z} x={shift}: measured tooth thickness {measured:.6} vs standard {expected:.6} ({rel:.2e})"
        );
    }
}

#[test]
fn gear_flank_normals_are_tangent_to_the_base_circle() {
    // The defining property of an involute flank: its normal is tangent to
    // the base circle, so the distance from the centre to the flank's
    // tangent line is exactly r_b. This is independent of the formula used
    // to place the points.
    let gear = SpurGear { flank_points: 200, ..SpurGear::new(2.0, 24, 5.0) };
    let region = gear.region();
    let rb = gear.base_radius();
    let (ra, rf) = (gear.tip_radius(), gear.root_radius());
    let mut worst: f64 = 0.0;
    let mut samples = 0;
    for i in 1..20 {
        // Sample radii strictly between the base circle and the tip.
        let r = rb.max(rf) * 1.02 + (ra * 0.98 - rb.max(rf) * 1.02) * i as f64 / 20.0;
        let psi = measured_half_angle(&region, r);
        let p = [r * psi.cos(), r * psi.sin()];
        // Outward normal by central differences of the exact region field.
        let e = 1e-6;
        let g = [
            (region.sdf([p[0] + e, p[1]]) - region.sdf([p[0] - e, p[1]])) / (2.0 * e),
            (region.sdf([p[0], p[1] + e]) - region.sdf([p[0], p[1] - e])) / (2.0 * e),
        ];
        let gn = (g[0] * g[0] + g[1] * g[1]).sqrt();
        let n = [g[0] / gn, g[1] / gn];
        // The involute's normal line passes tangent to the base circle, so
        // the perpendicular distance from the centre to that line, |p x n|,
        // must equal r_b.
        let lever = (p[0] * n[1] - p[1] * n[0]).abs();
        worst = worst.max((lever - rb).abs() / rb);
        samples += 1;
    }
    assert!(samples > 15);
    // The flank is a chordal polyline: with 200 points over the flank its
    // normal differs from the involute's by the chord's own slope error.
    assert!(
        worst < 2e-3,
        "flank normal lever arm differs from r_b by {worst:.2e} relative"
    );

    // The chordal error of the flank falls with the sampling, which is
    // what sets the tooth-thickness accuracy away from the pitch circle.
    let coarse = SpurGear { flank_points: 6, ..SpurGear::new(2.0, 24, 5.0) };
    // Not a round fraction, so the radius does not land on a sampled
    // vertex of either gear.
    let lo = gear.base_radius().max(gear.root_radius());
    let r_mid = lo + 0.3717 * (gear.tip_radius() - lo);
    let err = |g: &SpurGear| {
        let psi_measured = measured_half_angle(&g.region(), r_mid);
        (psi_measured - g.half_angle(r_mid).unwrap()).abs()
    };
    let (e6, e200) = (err(&coarse), err(&gear));
    println!("flank chord error: {e6:.3e} at 6 points, {e200:.3e} at 200");
    assert!(e6 > 1e-8, "the coarse flank should show a chord error, got {e6:.2e}");
    assert!(e200 < e6 / 100.0, "chord error {e6:.2e} -> {e200:.2e} with 6 -> 200 points");
}

#[test]
fn gear_diameters_and_mesh_geometry() {
    let gear = SpurGear { bore: 10.0, ..SpurGear::new(2.0, 20, 6.0) };
    let r = gear.report();
    assert!((r.pitch_diameter - 40.0).abs() < 1e-12);
    assert!((r.base_diameter - 40.0 * 20f64.to_radians().cos()).abs() < 1e-12);
    assert!((r.tip_diameter - 44.0).abs() < 1e-12);
    assert!((r.root_diameter - 35.0).abs() < 1e-12);
    assert!((r.base_pitch - PI * 2.0 * 20f64.to_radians().cos()).abs() < 1e-12);
    assert!(!r.tip_truncated);
    // 20 teeth at 20 degrees is just above the undercut limit 2/sin^2 a.
    assert!(r.undercut_limit > 17.0 && r.undercut_limit < 17.2, "{}", r.undercut_limit);

    // The solid's own extent must match the tip and bore.
    let region = gear.region();
    assert!(region.sdf([22.1, 0.0]) > 0.0 && region.sdf([21.8, 0.0]) < 0.0);
    assert!(region.sdf([4.9, 0.0]) > 0.0, "the bore must be open");
    assert!(region.sdf([5.2, 0.0]) < 0.0);

    // A standard pair: the contact ratio from the ISO 21771 formula.
    let pinion = SpurGear::new(2.0, 20, 6.0);
    let wheel = SpurGear::new(2.0, 20, 6.0);
    assert!((pinion.centre_distance(&wheel) - 40.0).abs() < 1e-12);
    let cr = pinion.contact_ratio(&wheel);
    assert!((cr - 1.557).abs() < 2e-3, "contact ratio {cr:.4}");
}

#[test]
fn thread_dimensions_match_iso_724() {
    // ISO 724 basic dimensions, millimetres: (nominal, pitch, d2, d1).
    let table = [
        (6.0, 1.0, 5.350, 4.917),
        (8.0, 1.25, 7.188, 6.647),
        (10.0, 1.5, 9.026, 8.376),
        (16.0, 2.0, 14.701, 13.835),
        (24.0, 3.0, 22.051, 20.752),
    ];
    for (d, p, d2, d1) in table {
        let t = Thread::new(d, p, 10.0);
        let r = t.report();
        assert!(
            (r.pitch_diameter - d2).abs() < 5e-4,
            "M{d}x{p}: d2 = {:.4}, ISO 724 gives {d2}",
            r.pitch_diameter
        );
        assert!(
            (r.minor - d1).abs() < 5e-4,
            "M{d}x{p}: d1 = {:.4}, ISO 724 gives {d1}",
            r.minor
        );
        // The coarse-pitch table must agree with the entry.
        let coarse = Thread::coarse_mm(d, 10.0).expect("coarse pitch");
        assert!((coarse.pitch - p).abs() < 1e-12, "ISO 261 pitch for M{d}");
    }
}

#[test]
fn thread_solid_has_the_right_profile_and_lead() {
    let t = Thread::new(8.0, 1.25, 6.0);
    let r = t.report();
    let solid = t.solid();

    // Radial extent: crest at d/2, root at d1/2, measured along a ray at
    // the height of a crest.
    let z = 3.0;
    let inside = |rad: f64, theta: f64, z: f64| solid.sdf([rad * theta.cos(), rad * theta.sin(), z]) < 0.0;
    assert!(!inside(r.major / 2.0 + 0.02, 0.0, z), "material beyond the major diameter");
    assert!(inside(r.minor / 2.0 - 0.02, 0.0, z), "the core must be solid");

    // Axial thread width at a radius, as a fraction of the pitch.
    let width_fraction = |rad: f64| {
        let n = 4000;
        let mut c = 0;
        for i in 0..n {
            let zz = 2.0 + 1.25 * i as f64 / n as f64;
            if inside(rad, 0.0, zz) {
                c += 1;
            }
        }
        c as f64 / n as f64
    };
    // At the pitch diameter the basic profile is half metal, half groove.
    let f2 = width_fraction(r.pitch_diameter / 2.0);
    assert!((f2 - 0.5).abs() < 0.01, "width fraction at d2 is {f2:.4}, expected 0.5");
    // At the crest it is the P/8 flat.
    let fc = width_fraction(r.major / 2.0 - 1e-4);
    assert!((fc - 0.125).abs() < 0.01, "width fraction at the crest is {fc:.4}, expected 0.125");

    // Lead: the crest advances exactly one pitch per turn, so the solid is
    // invariant under (theta, z) -> (theta + 2 pi, z + P).
    let rad = r.pitch_diameter / 2.0;
    let mut mismatches = 0;
    for i in 0..200 {
        let theta = 2.0 * PI * i as f64 / 200.0;
        let zz = 2.5 + 0.3 * (i % 7) as f64 / 7.0;
        if inside(rad, theta, zz) != inside(rad, theta + 2.0 * PI, zz + 1.25) {
            mismatches += 1;
        }
    }
    assert_eq!(mismatches, 0, "the helix does not close on itself");

    // The helical field keeps the kernel bound: its Lipschitz factor is
    // small but must be divided out.
    assert!(r.helical_factor > 1.0 && r.helical_factor < 1.01, "{}", r.helical_factor);
}

#[test]
fn octet_truss_density_matches_the_published_relation() {
    let (a, rad) = (0.01, 0.0005);
    let spec = LatticeSpec::new(Topology::Octet, a, rad);
    let graph = spec.graph();
    assert_eq!(graph.struts.len(), Topology::Octet.struts_per_cell());
    assert_eq!(graph.nodes.len(), 4);
    // 24 struts of length a/sqrt(2).
    let expected_length = 24.0 * a / 2.0f64.sqrt();
    assert!(
        (graph.length() - expected_length).abs() < 1e-15,
        "strut length {} vs {expected_length}",
        graph.length()
    );
    // And the Deshpande-Fleck-Ashby relative density.
    let published = octet_relative_density(rad, a);
    assert!(
        (graph.relative_density_upper() - published).abs() / published < 1e-14,
        "{} vs {published}",
        graph.relative_density_upper()
    );

    // The meshed solid's volume fraction is lower, by the material shared
    // at the nodes.
    let solid = graph.solid([1, 1, 1]);
    let mesh = mesh_solid(&solid, a / 80.0).expect("mesh");
    let fraction = mesh.report().volume / a.powi(3);
    println!(
        "octet r/l = {:.4}: measured volume fraction {fraction:.5}, no-overlap bound {published:.5}, \
         overlap {:.1} %",
        rad / (a / 2.0f64.sqrt()),
        100.0 * (1.0 - fraction / published)
    );
    assert!(
        fraction < published,
        "measured fraction {fraction:.5} must stay under the no-overlap bound {published:.5}"
    );
    assert!(
        fraction > 0.8 * published,
        "measured fraction {fraction:.5} is far below the bound {published:.5}"
    );

    for topology in [Topology::SimpleCubic, Topology::BodyCentred, Topology::CubicDiagonal] {
        let g = LatticeSpec::new(topology, a, rad).graph();
        assert_eq!(g.struts.len(), topology.struts_per_cell(), "{topology:?}");
    }
}

#[test]
fn octet_truss_stiffness_agrees_with_architected() {
    // The same lattice through two independent paths: the CAD strut graph
    // and the generator inside `architected`. Their homogenised stiffness
    // must agree, which checks the graph's topology, node reduction and
    // periodic images.
    use architected::cell::{Cell3d, Kind, RawMember3d};

    let (a, side) = (1.0, 0.05);
    let graph = LatticeSpec::new(Topology::Octet, a, side / 2.0).graph();
    let section = |t: f64| (t * t, t.powi(4) / 12.0, 0.1406 * t.powi(4));
    let (area, i, j) = section(side);
    let raw: Vec<RawMember3d> = graph
        .members()
        .into_iter()
        .map(|(p, q)| RawMember3d { a: p, b: q, e: 1.0, area, iy: i, iz: i, j })
        .collect();
    let ours = Cell3d::new([a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a], &raw, Kind::Truss);

    let (reference_raw, a1, a2, a3) = architected::lattices::octet_truss(side, a);
    let theirs = Cell3d::new(a1, a2, a3, &reference_raw, Kind::Truss);

    assert_eq!(ours.members.len(), theirs.members.len(), "member count");
    let (c1, c2) = (ours.stiffness(), theirs.stiffness());
    let scale = c2[0][0].abs().max(1e-30);
    let mut worst: f64 = 0.0;
    for r in 0..6 {
        for c in 0..6 {
            worst = worst.max((c1[r][c] - c2[r][c]).abs() / scale);
        }
    }
    assert!(worst < 1e-10, "stiffness differs by {worst:e} of C11");
    assert!(c1[0][0] > 0.0);
}

#[test]
fn periodic_level_sets_fill_half_the_cell() {
    // Each nodal function is odd under a symmetry of its own lattice, so
    // the level-zero solid fills exactly half of a whole number of cells.
    for kind in ["gyroid", "schwarz_p", "diamond"] {
        let a = 0.01;
        let solid = TpmsSpec::new(kind, a).block([2, 2, 2]).solid();
        let bounds = ([0.0; 3], [2.0 * a, 2.0 * a, 2.0 * a]);
        let v = solid.volume_quadrature(bounds, a / 40.0);
        let fraction = v / (2.0 * a).powi(3);
        assert!(
            (fraction - 0.5).abs() < 5e-3,
            "{kind}: volume fraction {fraction:.5}, expected 0.5 by symmetry"
        );
    }
    // A gyroid sheet is watertight, and its wall thickness is
    // proportional to the nominal one but larger, because the nodal field
    // is normalised by a conservative gradient bound.
    let a = 0.01;
    let mut factors = Vec::new();
    for t in [0.0002, 0.0004] {
        let spec = TpmsSpec::new("gyroid", a).sheet(t).block([1, 1, 1]);
        let mesh = mesh_solid(&spec.solid(), a / 80.0).expect("mesh");
        assert!(mesh.report().watertight, "gyroid sheet at t = {t} is not watertight");
        let achieved = spec.wall_thickness(a / 80.0).expect("thickness");
        println!("gyroid sheet: nominal t = {t}, achieved {achieved:.6}, factor {:.3}", achieved / t);
        assert!(achieved > t, "the achieved wall must exceed the nominal field thickness");
        factors.push(achieved / t);
    }
    // The factor is 1/|grad g| averaged over the surface, so it tends to a
    // limit as the sheet thins; it is not exactly constant.
    for f in &factors {
        assert!((1.5..2.2).contains(f), "wall factor {f:.3} outside the expected range");
    }
    assert!(
        factors[0] > factors[1],
        "the factor should rise toward its thin-sheet limit: {factors:?}"
    );
}

#[test]
fn plate_and_boss_build_as_expected() {
    use cad::kits::part::{boss, rib, Plate};
    let outline = Region2::rounded_rect([0.0, 0.0], [0.05, 0.03], 0.008);
    let plate = Plate::new(outline, 0.004)
        .with_hole(Region2::circle([-0.04, 0.0], 0.0025))
        .with_hole(Region2::circle([0.04, 0.0], 0.0025));
    let solid = plate
        .solid()
        .plus(boss([0.0, 0.0], 0.008, 0.004, 0.01, 0.005))
        .plus(rib([-0.03, 0.0], [0.03, 0.0], 0.003, 0.004, 0.006));
    let mesh = mesh_solid(&solid, 0.0005).expect("mesh");
    let report = mesh.report();
    assert!(report.watertight, "{report:?}");
    // The holes are open, the boss bore is open, the rib is attached.
    assert!(solid.sdf([-0.04, 0.0, 0.002]) > 0.0);
    // Above the rib, the boss bore is open and its wall is solid.
    assert!(solid.sdf([0.0, 0.0, 0.012]) > 0.0, "the boss bore must be open");
    assert!(solid.sdf([0.0, 0.007, 0.012]) < 0.0, "the boss wall");
    assert!(solid.sdf([0.02, 0.0, 0.007]) < 0.0, "the rib");
    assert!(Solid::Empty.sdf([0.0; 3]).is_infinite());
}
