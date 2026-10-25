//! NURBS against the exact properties that define them: partition of
//! unity, the rational quadratic circle, exact surfaces of revolution,
//! knot insertion, and interpolation.

use std::f64::consts::PI;

use cad::nurbs::{clamped_knots, Nurbs, Surface};
use cad::solid::Solid;

#[test]
fn the_basis_is_a_partition_of_unity() {
    // Sum_i N_i,p(u) = 1 for every u in the domain, for every degree.
    for degree in 1..=5usize {
        let n = degree + 4;
        let control: Vec<[f64; 3]> = (0..n).map(|i| [i as f64, 0.0, 0.0]).collect();
        let curve = Nurbs::clamped(degree, control).expect("curve");
        let (lo, hi) = curve.domain();
        let mut worst: f64 = 0.0;
        for k in 0..=500 {
            let u = lo + (hi - lo) * k as f64 / 500.0;
            let sum: f64 = (0..n).map(|i| curve.basis(i, u)).sum();
            worst = worst.max((sum - 1.0).abs());
            // And every basis function is non-negative.
            assert!((0..n).all(|i| curve.basis(i, u) >= -1e-15));
        }
        assert!(worst < 1e-14, "degree {degree}: partition of unity off by {worst:e}");
    }
    // The clamped knot vector has the right multiplicities.
    let k = clamped_knots(3, 7);
    assert_eq!(k.len(), 11);
    assert_eq!(&k[..4], &[0.0, 0.0, 0.0, 0.0]);
    assert_eq!(&k[7..], &[1.0, 1.0, 1.0, 1.0]);
}

#[test]
fn the_rational_quadratic_circle_is_exact() {
    // This is the property that makes NURBS worth having: a rational
    // quadratic with weights 1, sqrt(2)/2, 1 on each quadrant is a
    // circular arc exactly, not approximately.
    let (c, r) = ([0.3, -0.2, 0.1], 0.47);
    let circle = Nurbs::circle(c, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], r);
    let mut worst: f64 = 0.0;
    for k in 0..=2000 {
        let u = k as f64 / 2000.0;
        let p = circle.eval(u);
        let d = ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2)).sqrt();
        worst = worst.max((d - r).abs() / r);
        assert!((p[2] - c[2]).abs() < 1e-14, "the circle must stay in its plane");
    }
    assert!(worst < 1e-14, "radius error {worst:e}");

    // Its circumference, by Gauss quadrature of the exact derivative.
    let length = circle.length(16);
    let rel = (length - 2.0 * PI * r).abs() / (2.0 * PI * r);
    assert!(rel < 1e-12, "circumference {length:e} vs 2 pi r ({rel:e})");

    // In a tilted plane too. The second direction is orthogonalised, so
    // a non-perpendicular pair still gives a circle, not an ellipse.
    let tilted = Nurbs::circle(c, [1.0, 1.0, 0.0], [0.0, 1.0, 1.0], r);
    for k in 0..=500 {
        let p = tilted.eval(k as f64 / 500.0);
        let d = ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2) + (p[2] - c[2]).powi(2)).sqrt();
        assert!((d - r).abs() / r < 1e-14, "tilted circle radius error");
    }

    // Partial arcs: a quarter has a quarter of the circumference.
    for quarters in 1..=4usize {
        let arc = Nurbs::arc(c, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], r, quarters).unwrap();
        let expect = 2.0 * PI * r * quarters as f64 / 4.0;
        let rel = (arc.length(16) - expect).abs() / expect;
        assert!(rel < 1e-12, "{quarters} quarters: length error {rel:e}");
    }
}

#[test]
fn derivatives_and_endpoints_are_right() {
    // A clamped curve starts and ends at its outer control points, with
    // the end tangent along the first and last control edges.
    let control = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [3.0, 2.0, 1.0], [4.0, 0.0, 1.0]];
    let curve = Nurbs::bezier(control.clone()).unwrap();
    let a = curve.eval(0.0);
    let b = curve.eval(1.0);
    for k in 0..3 {
        assert!((a[k] - control[0][k]).abs() < 1e-15);
        assert!((b[k] - control[3][k]).abs() < 1e-15);
    }
    // A cubic Bezier's end derivative is 3 (P1 - P0).
    let d0 = curve.derivative(0.0);
    for k in 0..3 {
        let want = 3.0 * (control[1][k] - control[0][k]);
        assert!((d0[k] - want).abs() < 1e-12, "end derivative {d0:?} vs {want}");
    }

    // The analytic derivative against central differences.
    let mut worst: f64 = 0.0;
    for i in 1..100 {
        let u = i as f64 / 100.0;
        let h = 1e-6;
        let fd: Vec<f64> = (0..3)
            .map(|k| (curve.eval(u + h)[k] - curve.eval(u - h)[k]) / (2.0 * h))
            .collect();
        let an = curve.derivative(u);
        for k in 0..3 {
            worst = worst.max((an[k] - fd[k]).abs());
        }
    }
    assert!(worst < 1e-7, "derivative vs finite differences: {worst:e}");

    // Same for a rational curve, where the quotient rule matters.
    let circle = Nurbs::circle([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1.0);
    let mut worst: f64 = 0.0;
    for i in 1..200 {
        let u = i as f64 / 200.0;
        if (u * 4.0).fract().abs() < 1e-9 {
            continue; // skip the knots, where the derivative jumps direction
        }
        let h = 1e-7;
        let an = circle.derivative(u);
        for k in 0..3 {
            let fd = (circle.eval(u + h)[k] - circle.eval(u - h)[k]) / (2.0 * h);
            worst = worst.max((an[k] - fd).abs());
        }
        // The tangent of a circle is perpendicular to its radius.
        let p = circle.eval(u);
        let t = circle.tangent(u);
        assert!((p[0] * t[0] + p[1] * t[1]).abs() < 1e-10, "tangent not perpendicular at u = {u}");
    }
    assert!(worst < 1e-6, "rational derivative vs finite differences: {worst:e}");
}

#[test]
fn knot_insertion_does_not_move_the_curve() {
    // The strongest single check of the evaluator: refining the knot
    // vector adds control points but must leave the curve identical.
    let control = vec![
        [0.0, 0.0, 0.0],
        [1.0, 2.0, 0.5],
        [3.0, 2.5, -0.5],
        [4.0, 0.0, 1.0],
        [6.0, 1.0, 0.0],
        [7.0, -1.0, 0.5],
    ];
    let curve = Nurbs::clamped(3, control).unwrap();
    let refined = curve.insert_knot(0.37, 2).unwrap().insert_knot(0.8, 1).unwrap();
    assert_eq!(refined.control.len(), curve.control.len() + 3);
    let mut worst: f64 = 0.0;
    for i in 0..=1000 {
        let u = i as f64 / 1000.0;
        let (a, b) = (curve.eval(u), refined.eval(u));
        for k in 0..3 {
            worst = worst.max((a[k] - b[k]).abs());
        }
    }
    assert!(worst < 1e-13, "knot insertion moved the curve by {worst:e}");

    // And on a rational curve, where the weights must be carried too.
    let circle = Nurbs::circle([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 2.0);
    let refined = circle.insert_knot(0.1, 1).unwrap();
    let mut worst: f64 = 0.0;
    for i in 0..=1000 {
        let u = i as f64 / 1000.0;
        let p = refined.eval(u);
        worst = worst.max(((p[0] * p[0] + p[1] * p[1]).sqrt() - 2.0).abs());
    }
    assert!(worst < 1e-13, "refined circle radius error {worst:e}");
}

#[test]
fn interpolation_passes_through_its_points() {
    let points: Vec<[f64; 3]> = (0..9)
        .map(|i| {
            let t = i as f64 / 8.0;
            [t, (3.0 * t).sin() * 0.4, t * t * 0.2]
        })
        .collect();
    for degree in [2usize, 3, 4] {
        let curve = Nurbs::interpolate(&points, degree).expect("interpolate");
        // Chord-length parameters, recomputed the same way.
        let mut chords = vec![0.0; points.len()];
        for i in 1..points.len() {
            chords[i] = chords[i - 1]
                + ((points[i][0] - points[i - 1][0]).powi(2)
                    + (points[i][1] - points[i - 1][1]).powi(2)
                    + (points[i][2] - points[i - 1][2]).powi(2))
                .sqrt();
        }
        let total = chords[points.len() - 1];
        let mut worst: f64 = 0.0;
        for (i, p) in points.iter().enumerate() {
            let q = curve.eval(chords[i] / total);
            for k in 0..3 {
                worst = worst.max((q[k] - p[k]).abs());
            }
        }
        assert!(worst < 1e-10, "degree {degree}: interpolation misses by {worst:e}");
    }
}

#[test]
fn surfaces_of_revolution_are_exact_spheres_and_tori() {
    // A semicircular profile revolved is a sphere, exactly.
    let r = 0.73;
    let profile = Nurbs::arc([0.0; 3], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0], r, 2).unwrap();
    let sphere = Surface::revolve(&profile, 4).expect("revolve");
    let mut worst: f64 = 0.0;
    for i in 0..=120 {
        for j in 0..=60 {
            let p = sphere.eval(i as f64 / 120.0, j as f64 / 60.0);
            let d = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            worst = worst.max((d - r).abs() / r);
        }
    }
    assert!(worst < 1e-13, "sphere radius error {worst:e}");

    // Its area, from the tessellation, converges to 4 pi r^2.
    let exact = 4.0 * PI * r * r;
    let mut errors = Vec::new();
    for n in [16usize, 32, 64] {
        let a = sphere.area(2 * n, n, true);
        errors.push((a - exact).abs() / exact);
    }
    for w in errors.windows(2) {
        let order = (w[0] / w[1]).log2();
        assert!(order > 1.8, "sphere area order {order:.2} ({errors:?})");
    }
    // Measured 5.1e-4 at 128 x 64; the inscribed tessellation's O(h^2).
    assert!(errors[2] < 1e-3, "finest sphere area error {:e}", errors[2]);

    // A circle offset from the axis gives an exact torus.
    let (big, small) = (0.5, 0.12);
    let ring = Nurbs::circle([big, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], small);
    let torus = Surface::revolve(&ring, 4).expect("revolve");
    let mut worst: f64 = 0.0;
    for i in 0..=120 {
        for j in 0..=120 {
            let p = torus.eval(i as f64 / 120.0, j as f64 / 120.0);
            let rho = (p[0] * p[0] + p[1] * p[1]).sqrt();
            let d = ((rho - big).powi(2) + p[2] * p[2]).sqrt();
            worst = worst.max((d - small).abs() / small);
        }
    }
    assert!(worst < 1e-13, "torus error {worst:e}");

    // The torus area converges to 4 pi^2 R r (Pappus).
    let exact = 4.0 * PI * PI * big * small;
    let a = torus.area(160, 120, true);
    assert!((a - exact).abs() / exact < 1e-3, "torus area {a:e} vs {exact:e}");

    // The profile must lie in the xz half-plane, and say so if it does not.
    let bad = Nurbs::circle([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.1);
    assert!(Surface::revolve(&bad, 4).is_err());
}

#[test]
fn a_tessellated_surface_becomes_a_solid() {
    // The point of the module: a smooth surface closes the loop back into
    // the implicit kernel, so it can be cut, shelled and meshed.
    let r = 0.2;
    let profile = Nurbs::arc([0.0; 3], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0], r, 2).unwrap();
    let sphere = Surface::revolve(&profile, 4).unwrap();
    let mesh = sphere.tessellate(64, 32, true);
    let report = mesh.report();
    assert!(report.watertight, "a closed surface of revolution must tessellate watertight: {report:?}");
    let rel = (report.volume - 4.0 / 3.0 * PI * r.powi(3)).abs() / (4.0 / 3.0 * PI * r.powi(3));
    // Measured 4.1e-3 at 64 x 32: the inscribed polyhedron's deficit.
    assert!(rel < 6e-3, "tessellated volume error {rel:e}");

    let solid = Solid::from_trimesh(mesh);
    assert!(solid.exact(), "a mesh solid carries an exact signed distance");
    // Cut it in half and the result is still a solid the kernel can mesh.
    let half = solid.minus(Solid::half_space([0.0, 0.0, -1.0], 0.0));
    let cut = cad::mesh_solid(&half, r / 20.0).expect("mesh");
    let v = cut.report().volume;
    let expect = 2.0 / 3.0 * PI * r.powi(3);
    assert!((v - expect).abs() / expect < 0.01, "half-sphere volume {v:e} vs {expect:e}");
}

#[test]
fn curves_feed_sweeps_and_sketches() {
    // A NURBS path swept with a circular profile: a bent tube.
    let spine = Nurbs::interpolate(
        &[[0.0, 0.0, 0.0], [0.05, 0.02, 0.0], [0.1, 0.0, 0.01], [0.15, -0.02, 0.0]],
        3,
    )
    .unwrap();
    let path = spine.to_path_by_length(60, false);
    // Equal arc-length spacing: every segment the same length.
    let segs: Vec<f64> = path
        .points
        .windows(2)
        .map(|w| {
            ((w[1][0] - w[0][0]).powi(2) + (w[1][1] - w[0][1]).powi(2) + (w[1][2] - w[0][2]).powi(2)).sqrt()
        })
        .collect();
    let mean = segs.iter().sum::<f64>() / segs.len() as f64;
    let spread = segs.iter().map(|s| (s - mean).abs()).fold(0.0, f64::max) / mean;
    assert!(spread < 0.02, "arc-length sampling is uneven by {spread:.4}");
    // And the path length matches the curve's own quadrature.
    let rel = (path.length() - spine.length(32)).abs() / spine.length(32);
    assert!(rel < 1e-4, "sampled length off by {rel:e}");

    let radius = 0.004;
    let (mesh, check) = cad::path::sweep_checked(&path, &cad::path::Profile::circle(radius, 24));
    assert!(check.watertight, "swept tube: {check:?}");
    assert!(check.curvature_ratio < 1.0, "tube self-intersects: {check:?}");
    // The tube's area is the lateral area plus two caps. A 24-sided
    // profile is inscribed, so the mesh is slightly under the smooth
    // tube's area, by the usual 1 - sin(pi/n)/(pi/n).
    let smooth = 2.0 * PI * radius * path.length() + 2.0 * PI * radius * radius;
    let inscribed = (PI / 24.0).sin() / (PI / 24.0);
    let rel = (mesh.area() - smooth * inscribed).abs() / (smooth * inscribed);
    assert!(rel < 0.02, "swept tube area off by {rel:e}");

    // A planar curve becomes a sketch region.
    let outline = Nurbs::circle([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.03);
    let exact = PI * 0.03 * 0.03;
    // Uniform in parameter: inscribed, so below the circle, but not the
    // equal-angle polygon, because a rational circle's parameter is not
    // its angle.
    let by_param = outline.to_polygon(128).area().abs();
    assert!(by_param < exact && by_param / exact > 0.999, "{:.8}", by_param / exact);
    // Equally spaced in arc length, it is the equal-angle inscribed
    // n-gon, whose area is (n/2) sin(2 pi/n) / pi of the circle's.
    let by_length = outline.to_polygon_by_length(128).area().abs();
    let inscribed = 128.0 / 2.0 * (2.0 * PI / 128.0).sin() / PI;
    assert!(
        (by_length / exact - inscribed).abs() < 1e-6,
        "arc-length polygon holds {:.9} of the circle, expected {inscribed:.9}",
        by_length / exact
    );
    assert!(by_length > by_param, "even spacing must enclose more area");
    let poly = outline.to_polygon_by_length(128);
    let plate = Solid::extrude(cad::Region2::Polygon(poly), 0.0, 0.002);
    assert!(plate.sdf([0.0, 0.0, 0.001]) < 0.0);
}

#[test]
fn malformed_curves_are_refused() {
    assert!(Nurbs::new(2, vec![[0.0; 3]], vec![1.0], vec![0.0, 0.0, 0.0]).is_err());
    assert!(Nurbs::new(0, vec![[0.0; 3]; 3], vec![1.0; 3], vec![0.0; 4]).is_err());
    assert!(Nurbs::new(2, vec![[0.0; 3]; 4], vec![1.0; 3], clamped_knots(2, 4)).is_err());
    assert!(Nurbs::new(2, vec![[0.0; 3]; 4], vec![1.0, 1.0, -1.0, 1.0], clamped_knots(2, 4)).is_err());
    let mut bad = clamped_knots(2, 4);
    bad[3] = -1.0;
    assert!(Nurbs::new(2, vec![[0.0; 3]; 4], vec![1.0; 4], bad).is_err());
    assert!(Nurbs::interpolate(&[[0.0; 3], [0.0; 3]], 3).is_err());
    assert!(Nurbs::arc([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1.0, 5).is_err());
}
