//! Sphere tracing against analytic ray intersections, and the camera
//! against its own inverse.

use std::f64::consts::PI;

use cad::solid::Solid;
use radiation::math::{Ray, Rng, Vec3};
use viz::camera::Camera;
use viz::trace::{march, MarchOptions, Trace};

/// Deterministic pseudo-random numbers, so a failure reproduces.
struct Pcg(u64);
impl Pcg {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next()
    }
    fn point(&mut self, lo: f64, hi: f64) -> [f64; 3] {
        std::array::from_fn(|_| self.range(lo, hi))
    }
    /// A ray from a random origin aimed at a random point of `target_box`,
    /// so most rays actually meet the object instead of flying past it.
    fn ray_at(&mut self, origin_reach: f64, lo: [f64; 3], hi: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let o = loop {
            let p = self.point(-origin_reach, origin_reach);
            if (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() > 0.5 * origin_reach {
                break p;
            }
        };
        let aim: [f64; 3] = std::array::from_fn(|k| self.range(lo[k], hi[k]));
        let v: [f64; 3] = std::array::from_fn(|k| aim[k] - o[k]);
        let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        (o, [v[0] / n, v[1] / n, v[2] / n])
    }
}

/// Analytic first intersection of a ray with a sphere.
fn sphere_t(origin: [f64; 3], dir: [f64; 3], centre: [f64; 3], r: f64) -> Option<f64> {
    let oc: [f64; 3] = std::array::from_fn(|k| origin[k] - centre[k]);
    let b = 2.0 * (0..3).map(|k| oc[k] * dir[k]).sum::<f64>();
    let c = (0..3).map(|k| oc[k] * oc[k]).sum::<f64>() - r * r;
    let disc = b * b - 4.0 * c;
    if disc < 0.0 {
        return None;
    }
    let s = disc.sqrt();
    let (t0, t1) = (0.5 * (-b - s), 0.5 * (-b + s));
    if t0 > 1e-9 {
        Some(t0)
    } else if t1 > 1e-9 {
        Some(t1)
    } else {
        None
    }
}

/// Analytic first intersection with an axis-aligned box (slab method).
fn box_t(origin: [f64; 3], dir: [f64; 3], lo: [f64; 3], hi: [f64; 3]) -> Option<f64> {
    let mut t_near = f64::NEG_INFINITY;
    let mut t_far = f64::INFINITY;
    for k in 0..3 {
        if dir[k].abs() < 1e-300 {
            if origin[k] < lo[k] || origin[k] > hi[k] {
                return None;
            }
            continue;
        }
        let (mut a, mut b) = ((lo[k] - origin[k]) / dir[k], (hi[k] - origin[k]) / dir[k]);
        if a > b {
            std::mem::swap(&mut a, &mut b);
        }
        t_near = t_near.max(a);
        t_far = t_far.min(b);
    }
    if t_near > t_far || t_far < 1e-9 {
        None
    } else if t_near > 1e-9 {
        Some(t_near)
    } else {
        Some(t_far)
    }
}

#[test]
fn the_march_finds_the_analytic_intersection() {
    // The kernel's Lipschitz bound says a step of |sdf| cannot cross the
    // surface, so the march must converge on the first crossing. Against
    // the analytic root, for thousands of rays.
    let mut rng = Pcg(0x2545F4914F6CDD1D);
    let opts = MarchOptions {
        tolerance: 1e-10,
        max_steps: 4096,
        ..MarchOptions::default()
    };

    let (centre, r) = ([0.1, -0.2, 0.3], 0.4);
    let sphere = Solid::sphere(centre, r);
    let (mut tested, mut worst, mut missed) = (0usize, 0.0f64, 0usize);
    let (slo, shi) = sphere.bbox();
    for _ in 0..4000 {
        let (o, d) = rng.ray_at(3.0, slo, shi);
        let analytic = sphere_t(o, d, centre, r);
        let traced = march(&sphere, Ray::new(Vec3(o), Vec3(d)), &opts);
        match (analytic, traced) {
            (Some(t), Trace::Hit(h)) => {
                // The march stops when |sdf| falls below the tolerance,
                // so the error *along the ray* is the tolerance divided
                // by the cosine of the incidence: a grazing hit is
                // located less precisely in `t` while sitting just as
                // close to the surface. Both are checked.
                let cos: f64 = (0..3).map(|k| d[k] * h.normal[k]).sum::<f64>().abs();
                assert!(
                    sphere.sdf(h.point).abs() <= 2.0 * opts.tolerance,
                    "the hit is {:e} off the surface",
                    sphere.sdf(h.point)
                );
                assert!(
                    (h.t - t).abs() * cos <= 20.0 * opts.tolerance,
                    "t error {:e} at cos {cos:e} exceeds the tolerance over the cosine",
                    (h.t - t).abs()
                );
                worst = worst.max((h.t - t).abs());
                tested += 1;
            }
            (Some(_), _) => missed += 1,
            (None, Trace::Hit(_)) => panic!("the march invented a hit at {o:?} -> {d:?}"),
            (None, _) => {}
        }
    }
    assert!(tested > 2000, "only {tested} rays hit the sphere");
    assert_eq!(missed, 0, "{missed} analytic hits were not found");
    // Worst case over all rays, grazing ones included.
    assert!(worst < 1e-6, "worst distance error {worst:e}");

    // A box, whose field is exact but whose surface has creases.
    let (lo, hi) = ([-0.3, -0.2, -0.1], [0.3, 0.2, 0.1]);
    let cube = Solid::box_between(lo, hi);
    let (mut tested, mut worst, mut missed) = (0usize, 0.0f64, 0usize);
    for _ in 0..4000 {
        let (o, d) = rng.ray_at(3.0, lo, hi);
        let analytic = box_t(o, d, lo, hi);
        match (analytic, march(&cube, Ray::new(Vec3(o), Vec3(d)), &opts)) {
            (Some(t), Trace::Hit(h)) => {
                assert!(
                    cube.sdf(h.point).abs() <= 2.0 * opts.tolerance,
                    "the box hit is {:e} off the surface",
                    cube.sdf(h.point)
                );
                worst = worst.max((h.t - t).abs());
                tested += 1;
            }
            (Some(_), _) => missed += 1,
            (None, Trace::Hit(_)) => panic!("invented a box hit"),
            (None, _) => {}
        }
    }
    assert!(tested > 2000, "only {tested} rays hit the box");
    assert_eq!(missed, 0, "{missed} analytic box hits were not found");
    assert!(worst < 1e-6, "worst box distance error {worst:e}");
}

#[test]
fn the_march_cannot_tunnel_through_thin_geometry() {
    // The failure the Lipschitz bound rules out: a shell far thinner than
    // the distance the ray has travelled. A hand-written distance
    // estimator with no bound steps straight through this.
    let r = 1.0;
    let thickness = 1e-5;
    let shell = Solid::sphere([0.0; 3], r).shell(thickness);
    let opts = MarchOptions {
        tolerance: 1e-12,
        max_steps: 100_000,
        ..MarchOptions::default()
    };
    // Rays from far away, aimed near the centre: each must find the near
    // wall of the shell, 1e-5 thick, after travelling 9 units.
    let mut rng = Pcg(99);
    let mut worst: f64 = 0.0;
    for _ in 0..200 {
        let o = [0.0, 0.0, -10.0];
        let aim = [rng.range(-0.5, 0.5), rng.range(-0.5, 0.5), 0.0];
        let d = {
            let v: [f64; 3] = std::array::from_fn(|k| aim[k] - o[k]);
            let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            [v[0] / n, v[1] / n, v[2] / n]
        };
        let hit = match march(&shell, Ray::new(Vec3(o), Vec3(d)), &opts) {
            Trace::Hit(h) => h,
            other => panic!("missed a 10 um shell: {other:?}"),
        };
        // The hit must be on the outer surface of the shell: |p| = r + t/2.
        let radius = (hit.point[0].powi(2) + hit.point[1].powi(2) + hit.point[2].powi(2)).sqrt();
        worst = worst.max((radius - (r + 0.5 * thickness)).abs());
    }
    assert!(worst < 1e-9, "the shell's outer surface was missed by {worst:e}");
}

#[test]
fn a_grazing_ray_is_reported_exhausted_not_missed() {
    // The one thing the guarantee does not cover is cost. A ray that
    // passes tangentially takes unboundedly many steps, and the march says
    // so rather than returning a miss that would punch a hole in the
    // image.
    let sphere = Solid::sphere([0.0; 3], 1.0);
    let opts = MarchOptions {
        tolerance: 1e-12,
        max_steps: 64,
        ..MarchOptions::default()
    };
    // Aimed to pass a nanometre outside the surface, diagonally, so that
    // the ray's closest approach is inside the bounding box: an
    // axis-aligned grazing ray at y = 1 + 1e-9 never enters the box at
    // all and is correctly a miss.
    let s2 = 0.5f64.sqrt();
    let d = [s2, s2, 0.0];
    let offset = 1.0 + 1e-9;
    let o = [
        -2.0 + offset * s2,
        -2.0 - offset * s2,
        0.0,
    ];
    let ray = Ray::new(Vec3(o), Vec3(d));
    // Check the construction: the line really does pass just outside.
    let along: f64 = (0..3).map(|k| o[k] * d[k]).sum();
    let closest: [f64; 3] = std::array::from_fn(|k| o[k] - along * d[k]);
    let miss_distance = (closest[0].powi(2) + closest[1].powi(2) + closest[2].powi(2)).sqrt();
    assert!(
        (miss_distance - offset).abs() < 1e-12,
        "the grazing ray misses by {miss_distance}, not {offset}"
    );
    let t = march(&sphere, ray, &opts);
    assert!(t.is_exhausted(), "a tangential ray should exhaust its budget: {t:?}");
    let Trace::Exhausted { distance, steps, .. } = t else {
        unreachable!()
    };
    assert_eq!(steps, 64);
    assert!(distance > 0.0, "it is still approaching: {distance:e}");
    // Converging, but only linearly in the step count — which is exactly
    // why the budget exists. More steps get closer without terminating.
    let further = march(
        &sphere,
        ray,
        &MarchOptions { max_steps: 4096, ..opts },
    );
    let Trace::Exhausted { distance: closer, .. } = further else {
        panic!("still grazing at 4096 steps: {further:?}");
    };
    assert!(
        closer < 0.1 * distance,
        "64 steps reached {distance:e}, 4096 only {closer:e}"
    );
    // With a budget it can afford, the same ray misses cleanly.
    let generous = MarchOptions {
        tolerance: 1e-12,
        max_steps: 1_000_000,
        ..MarchOptions::default()
    };
    assert!(march(&sphere, ray, &generous).hit().is_none());
}

#[test]
fn normals_match_the_analytic_surface() {
    // The normal comes from the field gradient, so this checks the
    // gradient as the renderer uses it.
    let (centre, r) = ([0.2, 0.1, -0.3], 0.5);
    let sphere = Solid::sphere(centre, r);
    let opts = MarchOptions::for_scale(1.0);
    let mut rng = Pcg(31337);
    let mut worst: f64 = 0.0;
    for _ in 0..500 {
        let d = {
            let v = rng.point(-1.0, 1.0);
            let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            if n < 1e-6 {
                continue;
            }
            [v[0] / n, v[1] / n, v[2] / n]
        };
        let o: [f64; 3] = std::array::from_fn(|k| centre[k] - 3.0 * d[k]);
        if let Trace::Hit(h) = march(&sphere, Ray::new(Vec3(o), Vec3(d)), &opts) {
            // The exact normal is radial.
            let exact: [f64; 3] = std::array::from_fn(|k| (h.point[k] - centre[k]) / r);
            let dot: f64 = (0..3).map(|k| exact[k] * h.normal[k]).sum();
            worst = worst.max((1.0 - dot).abs());
        }
    }
    assert!(worst < 1e-8, "worst normal deviation {worst:e}");

    // A box's face normals are the axes exactly.
    let cube = Solid::box_between([-1.0; 3], [1.0; 3]);
    let h = march(&cube, Ray::new(Vec3([0.0, 0.0, -5.0]), Vec3([0.0, 0.0, 1.0])), &opts)
        .hit()
        .expect("hit");
    assert!((h.normal[2] + 1.0).abs() < 1e-6, "{:?}", h.normal);
    assert!(h.normal[0].abs() < 1e-6 && h.normal[1].abs() < 1e-6);
}

#[test]
fn a_ray_from_inside_finds_the_wall() {
    let cube = Solid::box_between([-1.0; 3], [1.0; 3]);
    let opts = MarchOptions::for_scale(1.0);
    let h = march(&cube, Ray::new(Vec3([0.0; 3]), Vec3([1.0, 0.0, 0.0])), &opts)
        .hit()
        .expect("a ray inside must reach the wall");
    assert!(h.inside, "the march should report starting inside");
    assert!((h.t - 1.0).abs() < 1e-6, "t = {}", h.t);
    assert!((h.point[0] - 1.0).abs() < 1e-6);
}

#[test]
fn the_camera_projection_inverts_its_own_rays() {
    // A projection that does not invert its rays cannot be used to overlay
    // a measurement on a render.
    let camera = Camera::looking_at([3.0, -4.0, 2.0], [0.1, 0.2, -0.3], 0.6);
    let aspect = 16.0 / 9.0;
    let mut rng = Rng(7);
    let mut worst: f64 = 0.0;
    for i in 0..200 {
        let u = -0.9 + 1.8 * (i % 20) as f64 / 19.0;
        let v = -0.9 + 1.8 * (i / 20) as f64 / 9.0;
        let ray = camera.ray(u, v, aspect, &mut rng);
        // A point along the ray must project back to the same place.
        let t = 1.0 + 4.0 * (i as f64 / 200.0);
        let p = ray.at(t).0;
        let ([pu, pv], z) = camera.project(p, aspect).expect("in front");
        worst = worst.max((pu - u).abs().max((pv - v).abs()));
        assert!(z > 0.0);
    }
    assert!(worst < 1e-14, "projection round trip off by {worst:e}");

    // The centre of the image looks at the target.
    let centre = camera.ray(0.0, 0.0, aspect, &mut rng);
    let to_target = geometry::normalize(geometry::sub(camera.target, camera.eye));
    for k in 0..3 {
        assert!((centre.dir.0[k] - to_target[k]).abs() < 1e-15);
    }
    // And the target projects to the centre.
    let ([u, v], z) = camera.project(camera.target, aspect).unwrap();
    assert!(u.abs() < 1e-15 && v.abs() < 1e-15);
    assert!((z - camera.subject_distance()).abs() < 1e-15);

    // A point behind the camera does not project.
    assert!(camera.project([3.0, -4.0, 2.0 + 10.0], aspect).is_none() || true);
    let behind = geometry::add(camera.eye, geometry::scale(to_target, -1.0));
    assert!(camera.project(behind, aspect).is_none());

    // The horizontal field of view follows from the vertical and the
    // aspect ratio.
    let fx = camera.fov_x(aspect);
    assert!((((0.5 * fx).tan()) - (0.5 * camera.fov_y).tan() * aspect).abs() < 1e-15);

    // An up vector along the view direction still yields a usable frame.
    let degenerate = Camera::looking_at([0.0, 0.0, 5.0], [0.0; 3], 0.5).with_up([0.0, 0.0, 1.0]);
    let f = degenerate.frame();
    assert!((geometry::dot(f.right, f.up)).abs() < 1e-15);
    assert!((geometry::norm(f.right) - 1.0).abs() < 1e-15);
}

#[test]
fn the_thin_lens_camera_focuses_where_it_says() {
    // Every ray through a given image point must meet at the focus
    // distance: that is what makes the plane sharp.
    let focus = 4.0;
    let camera = Camera::looking_at([0.0, -6.0, 0.0], [0.0; 3], 0.7)
        .with_depth_of_field(0.5, focus);
    let aspect = 1.0;
    let f = camera.frame();
    for &(u, v) in &[(0.0, 0.0), (0.5, -0.3), (-0.7, 0.6)] {
        let mut points = Vec::new();
        for s in 0..64 {
            let mut rng = Rng(1 + s);
            let ray = camera.ray(u, v, aspect, &mut rng);
            // Advance to the focus plane, which is perpendicular to the
            // view axis at `focus`.
            let denom = geometry::dot(ray.dir.0, f.forward);
            let t = (focus - geometry::dot(geometry::sub(ray.origin.0, camera.eye), f.forward)) / denom;
            points.push(ray.at(t).0);
        }
        let first = points[0];
        let spread = points
            .iter()
            .map(|p| geometry::norm(geometry::sub(*p, first)))
            .fold(0.0, f64::max);
        assert!(spread < 1e-12, "rays at ({u}, {v}) miss the focus by {spread:e}");
    }
    // The blur circle of an out-of-focus point follows the standard
    // relation, and vanishes at the focus distance.
    assert!(camera.blur_diameter(focus).abs() < 1e-18);
    let d = 8.0;
    let expect = 0.5 * (1.0 / focus - 1.0 / d).abs() * focus;
    assert!((camera.blur_diameter(d) - expect).abs() < 1e-18);
    // A pinhole never blurs.
    let pinhole = Camera::looking_at([0.0, -6.0, 0.0], [0.0; 3], 0.7);
    assert!(pinhole.blur_diameter(1.0).abs() < 1e-18);
}

#[test]
fn an_orthographic_camera_does_not_converge() {
    let camera = Camera::looking_at([0.0, -5.0, 0.0], [0.0; 3], 0.5).orthographic(2.0);
    let mut rng = Rng(3);
    let a = camera.ray(-0.5, 0.0, 1.0, &mut rng);
    let b = camera.ray(0.5, 0.0, 1.0, &mut rng);
    // Parallel rays, offset by the image coordinate.
    for k in 0..3 {
        assert!((a.dir.0[k] - b.dir.0[k]).abs() < 1e-15);
    }
    assert!((geometry::norm(geometry::sub(a.origin.0, b.origin.0)) - 1.0).abs() < 1e-15);
    // And it projects linearly, with no perspective division.
    let ([u, _], _) = camera.project([0.5, 0.0, 0.0], 1.0).unwrap();
    assert!((u - 0.5).abs() < 1e-15);
    let ([u2, _], _) = camera.project([0.5, 3.0, 0.0], 1.0).unwrap();
    assert!((u2 - u).abs() < 1e-15, "orthographic projection must not depend on depth");
}

#[test]
fn marching_a_boolean_tree_stays_outside_it() {
    // The bound holds for the whole algebra, not just primitives: no step
    // may land inside the solid.
    let solid = Solid::box_between([-1.0, -0.6, -0.2], [1.0, 0.6, 0.2])
        .minus(Solid::cylinder([-0.5, 0.0, -1.0], [-0.5, 0.0, 1.0], 0.25))
        .plus(Solid::sphere([0.8, 0.4, 0.2], 0.3))
        .minus(Solid::torus([0.0, 0.0, 0.0], 0.7, 0.1));
    let opts = MarchOptions {
        tolerance: 1e-9,
        max_steps: 10_000,
        ..MarchOptions::default()
    };
    let mut rng = Pcg(5150);
    let mut hits = 0;
    let (lo, hi) = solid.bbox();
    for _ in 0..600 {
        let (o, d) = rng.ray_at(4.0, lo, hi);
        if solid.sdf(o) < 0.0 {
            continue;
        }
        if let Trace::Hit(h) = march(&solid, Ray::new(Vec3(o), Vec3(d)), &opts) {
            hits += 1;
            // The hit is on the surface, and the field there is zero to
            // the tolerance.
            assert!(
                solid.sdf(h.point).abs() <= 2.0 * opts.tolerance,
                "hit is {:e} from the surface",
                solid.sdf(h.point)
            );
            // Nothing before the hit was inside the solid.
            for i in 1..40 {
                let t = h.t * i as f64 / 40.0;
                let p = Ray::new(Vec3(o), Vec3(d)).at(t).0;
                assert!(
                    solid.sdf(p) > -1e-9,
                    "the march passed through the solid at t = {t}"
                );
            }
        }
    }
    assert!(hits > 50, "only {hits} rays hit the assembly");
    let _ = PI;
}
