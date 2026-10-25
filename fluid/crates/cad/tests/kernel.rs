//! Kernel validation: exact distances, the Lipschitz bound, exact area and
//! volume identities, and mesh convergence.

use std::f64::consts::PI;

use cad::meshing::{distance_to_mesh, mesh_solid, Grid3};
use cad::path::{loft, Path3, Profile};
use cad::sketch::{Polygon, Region2};
use cad::solid::Solid;
use cad::Params;

/// Deterministic pseudo-random points (xorshift), so a failure reproduces.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn point(&mut self, lo: f64, hi: f64) -> [f64; 3] {
        std::array::from_fn(|_| lo + (hi - lo) * self.next())
    }
}

#[test]
fn cuboid_distance_matches_validated_mesh_sdf() {
    // geometry::Sdf on a polyhedron is exact and is itself validated
    // against the analytic box in the geometry crate, so it is an
    // independent reference for our analytic box field.
    let mesh = geometry::primitives::cuboid([-0.3, -0.2, -0.1], [0.3, 0.2, 0.1]);
    let reference = geometry::Sdf::new(mesh);
    let ours = Solid::box_between([-0.3, -0.2, -0.1], [0.3, 0.2, 0.1]);
    let mut rng = Rng(0x2545F4914F6CDD1D);
    let mut worst: f64 = 0.0;
    for _ in 0..4000 {
        let p = rng.point(-0.6, 0.6);
        worst = worst.max((ours.sdf(p) - reference.distance(p)).abs());
    }
    assert!(worst < 1e-14, "box distance differs from the mesh SDF by {worst:e}");
}

#[test]
fn cylinder_and_cone_distances_are_exact() {
    // Independent closed form for a z-aligned capped cylinder.
    let (r, h) = (0.25, 0.4);
    let ours = Solid::cylinder_z([0.0; 3], r, h);
    let reference = |p: [f64; 3]| {
        let dr = (p[0] * p[0] + p[1] * p[1]).sqrt() - r;
        let dz = p[2].abs() - 0.5 * h;
        if dr <= 0.0 && dz <= 0.0 {
            dr.max(dz)
        } else {
            (dr.max(0.0).powi(2) + dz.max(0.0).powi(2)).sqrt()
        }
    };
    let mut rng = Rng(12345);
    let mut worst: f64 = 0.0;
    for _ in 0..4000 {
        let p = rng.point(-0.6, 0.6);
        worst = worst.max((ours.sdf(p) - reference(p)).abs());
    }
    assert!(worst < 1e-14, "cylinder distance error {worst:e}");

    // A lofted frustum surface is an independent reference for the cone.
    let (ra, rb, len) = (0.2, 0.08, 0.3);
    let path = Path3::line([0.0, 0.0, 0.0], [0.0, 0.0, len], 1);
    let n = 2048;
    let frustum = loft(&path, &[Profile::circle(ra, n), Profile::circle(rb, n)]);
    let cone = Solid::cone([0.0, 0.0, 0.0], [0.0, 0.0, len], ra, rb);
    let mut worst_cone: f64 = 0.0;
    for _ in 0..200 {
        let p = rng.point(-0.4, 0.5);
        let d = distance_to_mesh(&frustum, p);
        worst_cone = worst_cone.max((cone.sdf(p).abs() - d).abs());
    }
    // The polygonal loft inscribes the circle; its sagitta at n = 2048 is
    // r (1 - cos(pi/n)) = 2.4e-7 r.
    assert!(worst_cone < 1e-6, "cone distance error {worst_cone:e} against a {n}-gon loft");
}

#[test]
fn boolean_tree_keeps_the_sign_and_bounds_the_distance() {
    // A plate with two cross-drilled holes, a sphere welded on, and a
    // rounded pocket: unions, differences and offsets together.
    let plate = Solid::box_between([-0.2, -0.12, -0.02], [0.2, 0.12, 0.02]);
    let hole = |x: f64| Solid::cylinder([x, 0.0, -0.1], [x, 0.0, 0.1], 0.03);
    let pocket = Solid::rounded_box([0.0, 0.0, 0.025], [0.05, 0.05, 0.025], 0.01);
    let solid = plate
        .minus(hole(-0.1))
        .minus(hole(0.1))
        .minus(pocket)
        .plus(Solid::sphere([0.15, 0.08, 0.02], 0.03));
    assert!(!solid.exact(), "a boolean tree must not claim exactness");

    let h = 0.004;
    let mesh = mesh_solid(&solid, h).expect("mesh");
    let report = mesh.report();
    assert!(report.watertight, "meshed CSG not watertight: {report:?}");
    assert_eq!(report.nonmanifold_edges, 0);
    assert_eq!(report.inconsistent_edges, 0);

    let reference = geometry::Sdf::new(mesh.clone());
    let mut rng = Rng(777);
    let (mut checked, mut worst_over) = (0usize, 0.0f64);
    for _ in 0..1500 {
        let p = rng.point(-0.25, 0.25);
        let d_true = distance_to_mesh(&mesh, p);
        let f = solid.sdf(p);
        if d_true > 2.0 * h {
            worst_over = worst_over.max(f.abs() - d_true);
            // The sign must agree with the independent mesh sign.
            assert_eq!(
                f < 0.0,
                reference.distance(p) < 0.0,
                "sign disagreement at {p:?}: kernel {f:e}, mesh {:e}",
                reference.distance(p)
            );
            checked += 1;
        }
    }
    assert!(checked > 800, "only {checked} points were far enough from the surface");
    // The reference here is the *meshed* boundary, whose own position error
    // is O(h^2) = 1.6e-5 at this resolution, so the comparison can only
    // bound the excess at that level. The invariant itself is checked
    // directly below, where it is exact.
    assert!(
        worst_over < 5.0 * h * h,
        "kernel distance exceeded the meshed distance by {worst_over:e} (h^2 = {:e})",
        h * h
    );

    // The invariant that |sdf| <= dist is a consequence of the field being
    // 1-Lipschitz and vanishing on the boundary. The Lipschitz property is
    // directly checkable, and is where the proof actually rests.
    let mut worst_ratio: f64 = 0.0;
    for _ in 0..4000 {
        let a = rng.point(-0.3, 0.3);
        let b = rng.point(-0.3, 0.3);
        let d = geometry::norm(geometry::sub(b, a));
        if d > 1e-6 {
            worst_ratio = worst_ratio.max((solid.sdf(b) - solid.sdf(a)).abs() / d);
        }
    }
    assert!(
        worst_ratio <= 1.0 + 1e-12,
        "field is not 1-Lipschitz: worst |df|/|dx| = {worst_ratio}"
    );
}

#[test]
fn rounding_needs_a_remesh_and_offsets_alone_cancel() {
    let cube = Solid::box_between([0.0; 3], [0.1, 0.1, 0.1]);
    // Offsetting down and back up is the identity in a distance algebra:
    // this is why round_convex has to re-distance.
    let cancelled = cube.clone().offset(-0.01).offset(0.01);
    let p = [0.1, 0.1, 0.1];
    assert!(
        (cancelled.sdf(p) - cube.sdf(p)).abs() < 1e-15,
        "erode-then-dilate should cancel exactly"
    );
    // The real operation rounds the corner: the cube's corner is at
    // distance 0 from its vertex; after rounding by r the vertex is
    // outside the solid by r (sqrt 3 - 1).
    let rounded = cube.clone().round_convex(0.01, 0.002).expect("round");
    let d = rounded.sdf(p);
    let expect = 0.01 * (3.0f64.sqrt() - 1.0);
    assert!(
        (d - expect).abs() < 0.0015,
        "rounded corner distance {d:e}, expected {expect:e}"
    );
    // A face centre is untouched by the rounding.
    let face = [0.05, 0.05, 0.1];
    assert!(rounded.sdf(face).abs() < 0.0015, "face moved by {:e}", rounded.sdf(face));
}

#[test]
fn pappus_and_green_identities() {
    // Revolving a circle of radius a whose centre is at r = R gives a
    // torus of volume 2 pi R * pi a^2 (Pappus).
    let (big, small) = (0.3, 0.08);
    let region = Region2::circle([big, 0.0], small);
    let solid = Solid::revolve(region);
    assert!(solid.exact());
    let exact = 2.0 * PI * big * PI * small * small;
    for h in [0.01, 0.005] {
        let v = solid.volume_quadrature(solid.bbox(), h);
        let rel = (v - exact).abs() / exact;
        assert!(rel < 0.02, "Pappus volume error {rel:.4} at h = {h}");
    }

    // Extruding a polygon multiplies its exact Green's-theorem area.
    let poly = Polygon::new(vec![[0.0, 0.0], [0.4, 0.0], [0.4, 0.1], [0.15, 0.1], [0.15, 0.25], [0.0, 0.25]])
        .with_hole(vec![[0.05, 0.05], [0.1, 0.05], [0.1, 0.1], [0.05, 0.1]]);
    let area = poly.area();
    let expected = 0.4 * 0.1 + 0.15 * 0.15 - 0.05 * 0.05;
    assert!((area - expected).abs() < 1e-15, "polygon area {area} vs {expected}");
    let thickness = 0.03;
    let solid = Solid::extrude(Region2::Polygon(poly), 0.0, thickness);
    let v = solid.volume_quadrature(solid.bbox(), 0.0015);
    let rel = (v - area * thickness).abs() / (area * thickness);
    assert!(rel < 0.01, "extruded volume error {rel:.4}");

    // Second moments of a rectangle are exact: I_xx = b h^3 / 12.
    let (b, hh) = (0.06, 0.02);
    let rect = Polygon::new(vec![[0.0, 0.0], [b, 0.0], [b, hh], [0.0, hh]]);
    let (ixx, iyy, ixy) = rect.second_moments();
    assert!((ixx - b * hh.powi(3) / 12.0).abs() < 1e-20, "Ixx {ixx:e}");
    assert!((iyy - hh * b.powi(3) / 12.0).abs() < 1e-20, "Iyy {iyy:e}");
    assert!(ixy.abs() < 1e-20, "Ixy {ixy:e} should vanish by symmetry");
}

#[test]
fn surface_nets_is_second_order() {
    let r = 0.5;
    let sphere = Solid::sphere([0.0; 3], r);
    let (v_exact, a_exact) = (4.0 / 3.0 * PI * r.powi(3), 4.0 * PI * r * r);
    let mut errors = Vec::new();
    for n in [8usize, 16, 32] {
        let h = r / n as f64;
        let mesh = mesh_solid(&sphere, h).expect("mesh");
        let report = mesh.report();
        assert!(report.watertight, "sphere mesh at h = {h} is not watertight: {report:?}");
        let (ev, ea) = (
            (report.volume - v_exact).abs() / v_exact,
            (report.area - a_exact).abs() / a_exact,
        );
        errors.push((h, ev, ea));
    }
    for w in errors.windows(2) {
        let order_v = (w[0].1 / w[1].1).log2();
        let order_a = (w[0].2 / w[1].2).log2();
        assert!(
            order_v > 1.6,
            "volume convergence order {order_v:.2} (errors {:e} -> {:e})",
            w[0].1,
            w[1].1
        );
        assert!(order_a > 1.3, "area convergence order {order_a:.2}");
    }
    // Measured at h = r/32: 3.9e-4 relative in volume, 1.8e-4 in area.
    let (_, ev, ea) = errors[2];
    assert!(ev < 5e-4 && ea < 3e-4, "finest errors: volume {ev:e}, area {ea:e}");
}

#[test]
fn offsets_and_shells_are_exact_on_exact_solids() {
    let s = Solid::sphere([0.0; 3], 0.2);
    let grown = s.clone().offset(0.05);
    assert!(grown.exact());
    assert!((grown.sdf([0.3, 0.0, 0.0]) - (0.3 - 0.25)).abs() < 1e-15);
    let wall = s.offset(0.05).shell(0.02);
    // The shell straddles r = 0.25 with 0.01 either side.
    assert!(wall.sdf([0.25, 0.0, 0.0]) < 0.0);
    assert!(wall.sdf([0.235, 0.0, 0.0]) > 0.0);
    assert!(wall.sdf([0.265, 0.0, 0.0]) > 0.0);
}

#[test]
fn placement_round_trips_and_mirrors() {
    let s = Solid::cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 0.1], 0.02);
    let moved = s
        .clone()
        .rotated([1.0, 0.0, 0.0], PI / 2.0)
        .translated([0.3, 0.0, 0.0]);
    // The cylinder now runs along -y from (0.3, 0, 0).
    assert!(moved.sdf([0.3, -0.05, 0.0]) < 0.0);
    assert!(moved.sdf([0.3, 0.0, 0.05]) > 0.0);
    let mirrored = s.clone().translated([0.3, 0.0, 0.0]).mirrored([1.0, 0.0, 0.0]);
    assert!(mirrored.sdf([-0.3, 0.0, 0.05]) < 0.0);
    // A scaled solid's field scales with it.
    let scaled = Solid::sphere([0.0; 3], 1.0).scaled(0.25);
    assert!((scaled.sdf([0.5, 0.0, 0.0]) - 0.25).abs() < 1e-15);
}

#[test]
fn arrays_are_exact_when_the_shape_fits_its_cell() {
    let hole = Solid::sphere([0.0; 3], 0.01);
    let grid = hole.array([0.05, 0.05, 0.0], [4, 3, 1]);
    // Centres sit at +-0.075, +-0.025 in x and 0, +-0.05 in y.
    assert!(grid.sdf([0.075, 0.05, 0.0]) < 0.0);
    assert!((grid.sdf([0.075, 0.05, 0.02]) - 0.01).abs() < 1e-14);
    assert!(grid.sdf([0.125, 0.0, 0.0]) > 0.0, "the array must be finite");
    let mesh = mesh_solid(&grid, 0.001).expect("mesh");
    let report = mesh.report();
    assert!(report.watertight);
    let v_exact = 12.0 * 4.0 / 3.0 * PI * 0.01f64.powi(3);
    let rel = (report.volume - v_exact).abs() / v_exact;
    // The sphere study gives 3.9e-4 at h = r/32, so 4e-3 at r/10.
    assert!(rel < 5e-3, "12 spheres: volume error {rel:.4}");
}

#[test]
fn expressions_and_parameters() {
    let p = Params::from_defs(&[
        ("wall", "0.002"),
        ("bore", "0.01"),
        ("od", "bore + 2*wall"),
        ("area", "pi * (od/2)^2"),
    ])
    .expect("defs");
    assert!((p.get("od").unwrap() - 0.014).abs() < 1e-18);
    assert!((p.get("area").unwrap() - PI * 0.007 * 0.007).abs() < 1e-18);
    // Precedence and associativity.
    assert!((p.eval("2 + 3 * 4").unwrap() - 14.0).abs() < 1e-15);
    assert!((p.eval("2 ^ 3 ^ 2").unwrap() - 512.0).abs() < 1e-12);
    assert!((p.eval("-2^2").unwrap() + 4.0).abs() < 1e-15);
    assert!((p.eval("min(3, 1, 2) + max(1, 5)").unwrap() - 6.0).abs() < 1e-15);
    assert!((p.eval("rad(180)").unwrap() - PI).abs() < 1e-15);
    assert!((p.eval("1.5e-3 * 2").unwrap() - 3e-3).abs() < 1e-18);
    assert!(p.eval("nope + 1").is_err());
    assert!(p.eval("2 +").is_err());
}

#[test]
fn json_model_is_parametric() {
    let text = r#"{
      "name": "flange",
      "params": { "bolts": 6, "bolt_circle": 0.08, "hole": 0.0033,
                  "thickness": 0.006, "od": "bolt_circle + 6*hole" },
      "shape": {
        "op": "difference",
        "a": { "op": "cylinder_z", "center": [0,0,"thickness/2"],
               "radius": "od/2", "height": "thickness" },
        "b": { "op": "union", "parts": [
                 { "op": "cylinder_z", "center": [0,0,"thickness/2"],
                   "radius": "od/4", "height": "2*thickness" },
                 { "op": "polar_array", "axis": [0,0,1], "count": "bolts",
                   "shape": { "op": "cylinder_z",
                              "center": ["bolt_circle/2",0,"thickness/2"],
                              "radius": "hole", "height": "2*thickness" } } ] }
      }
    }"#;
    let model = cad::model::Model::from_json(text).expect("parse");
    let params = model.resolve(&Params::new()).expect("resolve");
    assert!((params.get("od").unwrap() - 0.0998).abs() < 1e-15);

    let analytic = |bolts: usize, od: f64| {
        let t = 0.006;
        PI * (od / 2.0).powi(2) * t - PI * (od / 4.0).powi(2) * t - bolts as f64 * PI * 0.0033f64.powi(2) * t
    };
    for (bolts, expect_od) in [(6usize, 0.0998), (12, 0.0998)] {
        let mut over = Params::new();
        over.set("bolts", bolts as f64);
        let solid = model.build(&over, ".").expect("build");
        let v = solid.volume_quadrature(solid.bbox(), 0.0006);
        let want = analytic(bolts, expect_od);
        let rel = (v - want).abs() / want;
        assert!(rel < 0.02, "{bolts} bolts: volume error {rel:.4} ({v:e} vs {want:e})");
    }
    // A sweep of the bolt circle moves the holes, so the volume is
    // unchanged but the geometry is not.
    let mut over = Params::new();
    over.set("bolt_circle", 0.07);
    let a = model.build(&Params::new(), ".").unwrap();
    let b = model.build(&over, ".").unwrap();
    let z = 0.003;
    // A bolt hole sits on the bolt circle and solid metal sits off it, and
    // the two models swap which radius is which.
    assert!(a.sdf([0.040, 0.0, z]) > 0.0, "hole expected at r = 0.040 in model a");
    assert!(a.sdf([0.035, 0.0, z]) < 0.0, "metal expected at r = 0.035 in model a");
    assert!(b.sdf([0.035, 0.0, z]) > 0.0, "hole expected at r = 0.035 in model b");
    assert!(b.sdf([0.040, 0.0, z]) < 0.0, "metal expected at r = 0.040 in model b");
}

#[test]
fn unbounded_solids_are_refused_not_guessed() {
    let gyroid = cad::kits::lattice::TpmsSpec::new("gyroid", 0.01).solid();
    let err = mesh_solid(&gyroid, 0.001).unwrap_err();
    assert!(err.contains("unbounded"), "unexpected error: {err}");
    let block = cad::kits::lattice::TpmsSpec::new("gyroid", 0.01).block([1, 1, 1]).solid();
    assert!(mesh_solid(&block, 0.0008).is_ok());
}

#[test]
fn grid_covering_pads_the_boundary() {
    let s = Solid::sphere([0.0; 3], 1.0);
    let g = Grid3::covering(s.bbox(), 0.25, 2).unwrap();
    let f = cad::meshing::field_on(&s, &g);
    let [nx, ny, nz] = g.dims;
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let on_boundary = i == 0 || j == 0 || k == 0 || i == nx - 1 || j == ny - 1 || k == nz - 1;
                if on_boundary {
                    assert!(f[g.index(i, j, k)] > 0.0, "node ({i},{j},{k}) is inside the padded boundary");
                }
            }
        }
    }
}

#[test]
fn polar_repeat_matches_the_explicit_union() {
    // Folding the angle must give the same region as unioning every copy,
    // which is what makes a many-toothed gear affordable.
    let tooth = Region2::rect([0.04, 0.0], [0.008, 0.003]);
    let count = 17;
    let folded = tooth.clone().polar_repeat(count);
    let explicit = Region2::Union {
        parts: (0..count)
            .map(|i| tooth.clone().rotated(2.0 * PI * i as f64 / count as f64))
            .collect(),
    };
    let mut rng = Rng(99);
    let mut worst: f64 = 0.0;
    let mut sign_mismatches = 0;
    for _ in 0..20000 {
        let p = [rng.next() * 0.12 - 0.06, rng.next() * 0.12 - 0.06];
        let (a, b) = (folded.sdf(p), explicit.sdf(p));
        worst = worst.max((a - b).abs());
        if (a < 0.0) != (b < 0.0) {
            sign_mismatches += 1;
        }
    }
    assert_eq!(sign_mismatches, 0, "folding changed the region");
    assert!(worst < 1e-12, "folded distance differs by {worst:e}");

    // And the exact area of the repeat is the copy's times the count.
    let a = folded.area_exact().unwrap();
    assert!((a - count as f64 * 4.0 * 0.008 * 0.003).abs() < 1e-18);
}

#[test]
fn capped_cylinder_boxes_are_tight() {
    // A ball at each end would inflate the box along the axis by the
    // radius; a short wide disc is the case where that matters.
    let disc = Solid::cylinder_z([0.0; 3], 0.05, 0.004);
    let (lo, hi) = disc.bbox();
    assert!((lo[2] + 0.002).abs() < 1e-15 && (hi[2] - 0.002).abs() < 1e-15, "z extent {lo:?}..{hi:?}");
    assert!((lo[0] + 0.05).abs() < 1e-15 && (hi[0] - 0.05).abs() < 1e-15);
    // A tilted cylinder's box still contains it.
    let tilted = Solid::cylinder([0.0; 3], [0.1, 0.1, 0.1], 0.02);
    let (lo, hi) = tilted.bbox();
    let mut rng = Rng(4242);
    for _ in 0..20000 {
        let p = rng.point(-0.2, 0.3);
        if tilted.sdf(p) < 0.0 {
            assert!(
                (0..3).all(|k| p[k] >= lo[k] - 1e-12 && p[k] <= hi[k] + 1e-12),
                "{p:?} is inside the solid but outside the box {lo:?}..{hi:?}"
            );
        }
    }
    // The cone's box too.
    let cone = Solid::cone([0.0, 0.0, 0.0], [0.0, 0.0, 0.01], 0.06, 0.01);
    let (lo, hi) = cone.bbox();
    assert!((hi[0] - 0.06).abs() < 1e-15 && (hi[2] - 0.01).abs() < 1e-15);
    assert!(lo[2].abs() < 1e-15);
}
