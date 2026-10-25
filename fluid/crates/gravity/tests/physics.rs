use gravity::{
    direct,
    metric::{Metric, Schwarzschild, WeakField},
    norm, Body, Direct, Field, NBody, PointMass, Spherical, Tree, Vec3, C, G,
};
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t, "{a:e} != {b:e}, tol {t:e}");
}
#[test]
fn uniform_sphere_has_exact_interior_exterior_and_shell_theorem() {
    let rho = 3.0 / (4.0 * std::f64::consts::PI * G);
    let field = Spherical::new(vec![0.0, 0.2, 0.7, 1.0], vec![rho; 3], Vec3::default()).unwrap();
    for r in [0.0, 0.01, 0.2, 0.5, 1.0, 2.0, 10.0] {
        let s = field.sample(Vec3([r, 0.0, 0.0])).unwrap();
        close(
            s.acceleration.0[0],
            if r <= 1.0 { -r } else { -1.0 / (r * r) },
            1e-14,
        );
        close(
            s.potential,
            if r <= 1.0 {
                -0.5 * (3.0 - r * r)
            } else {
                -1.0 / r
            },
            1e-14,
        );
    }
    let hollow = Spherical::new(vec![0.0, 1.0, 2.0], vec![0.0, rho], Vec3::default()).unwrap();
    let a = hollow.sample(Vec3::default()).unwrap();
    let b = hollow.sample(Vec3([0.5, 0.0, 0.0])).unwrap();
    assert_eq!(a.acceleration, Vec3::default());
    assert_eq!(b.acceleration, Vec3::default());
    assert_eq!(a.potential, b.potential);
}
fn cloud(n: usize) -> Vec<Body> {
    let mut rng = radiation::math::Rng(942);
    (0..n)
        .map(|_| Body {
            position: Vec3(std::array::from_fn(|_| 2.0 * rng.uniform() - 1.0)),
            mass: (0.5 + rng.uniform()) / G,
        })
        .collect()
}
#[test]
fn quadrupole_tree_matches_direct_and_excludes_sources() {
    let bodies = cloud(2048);
    let direct_tree = Tree::new(bodies.clone(), 0.02, 0.0).unwrap();
    let tree = Tree::new(bodies.clone(), 0.02, 0.3).unwrap();
    let mut square = 0.0;
    let mut total = 0.0;
    let mut work = 0;
    for i in (0..bodies.len()).step_by(16) {
        let x = bodies[i].position;
        let reference = direct(&bodies, x, 0.02, Some(i)).unwrap();
        let exact = direct_tree.evaluate(x, Some(i)).unwrap().0;
        close(
            norm(exact.acceleration - reference.acceleration),
            0.0,
            2e-11,
        );
        let (approx, w) = tree.evaluate(x, Some(i)).unwrap();
        square += norm(approx.acceleration - reference.acceleration).powi(2);
        total += norm(reference.acceleration).powi(2);
        work += w.direct_pairs;
    }
    let relative = (square / total).sqrt();
    assert!(relative < 0.001, "tree relative RMS {relative}");
    assert!(
        work < 2048 * 128 * 4 / 5,
        "tree did not reduce direct interactions"
    );
    let x = Vec3([30.0, 0.0, 0.0]);
    let a = tree.evaluate(x, Some(0)).unwrap().0;
    let b = direct(&bodies, x, 0.02, Some(0)).unwrap();
    close(a.potential, b.potential, 1e-5 * b.potential.abs());
}
#[test]
fn tree_potential_gradient_matches_its_quadrupole_force() {
    let tree = Tree::new(cloud(200), 0.1, 0.5).unwrap();
    let x = Vec3([20.0, 3.0, 5.0]);
    let s = tree.sample(x).unwrap();
    for d in 0..3 {
        let mut lo = x;
        let mut hi = x;
        lo.0[d] -= 1e-4;
        hi.0[d] += 1e-4;
        let gradient =
            (tree.sample(hi).unwrap().potential - tree.sample(lo).unwrap().potential) / 2e-4;
        close(gradient, -s.acceleration.0[d], 1e-9);
    }
}
#[test]
fn leapfrog_binary_preserves_momentum_and_bounded_orbital_energy() {
    let mut sim = NBody {
        bodies: vec![
            Body {
                position: Vec3([-0.5, 0.0, 0.0]),
                mass: 0.5 / G,
            },
            Body {
                position: Vec3([0.5, 0.0, 0.0]),
                mass: 0.5 / G,
            },
        ],
        velocities: vec![Vec3([0.0, -0.5, 0.0]), Vec3([0.0, 0.5, 0.0])],
        time: 0.0,
        softening: 0.0,
        theta: 0.0,
    };
    let initial = sim.invariants().unwrap();
    let dt = std::f64::consts::TAU / 1000.0;
    for _ in 0..20000 {
        sim.step(dt).unwrap();
        let current = sim.invariants().unwrap();
        close(current.0, initial.0, 1e-8 * initial.0.abs());
        assert_eq!(current.1, initial.1);
    }
    let before = sim.bodies[0].position;
    assert!(sim.step(f64::NAN).is_err());
    assert_eq!(sim.bodies[0].position, before);
}
#[test]
fn metric_and_source_domains_are_checked() {
    let metric = Schwarzschild {
        gm: 1e20,
        center: Vec3::default(),
    };
    let r = 1e9;
    let g = metric.geometry(Vec3([r, 0.0, 0.0])).unwrap();
    let weak = WeakField(PointMass {
        gm: 1e20,
        center: Vec3::default(),
        softening: 0.0,
    })
    .geometry(Vec3([r, 0.0, 0.0]))
    .unwrap();
    close(g.lapse, weak.lapse, 5e-12);
    close(g.scale, weak.scale, 5e-12);
    assert!(metric.geometry(Vec3::default()).is_err());
    assert!(metric
        .geometry(Vec3([1e20 / (4.0 * C * C), 0.0, 0.0]))
        .is_err());
    assert!(Spherical::new(vec![0.0, 1.0, 0.5], vec![1.0; 2], Vec3::default()).is_err());
    assert!(Tree::new(cloud(10), -1.0, 0.5).is_err());
}

#[test]
fn direct_field_preserves_reference_and_scaled_norm_extremes() {
    let sources = cloud(64);
    let field = Direct::new(sources.clone(), 0.02).unwrap();
    for (i, body) in sources.iter().enumerate() {
        let actual = field.evaluate(body.position, Some(i)).unwrap();
        let expected = direct(&sources, body.position, 0.02, Some(i)).unwrap();
        assert_eq!(actual.acceleration, expected.acceleration);
        assert_eq!(actual.potential, expected.potential);
    }
    assert!(field.evaluate(Vec3::default(), Some(64)).is_err());
    assert!(direct(&sources, Vec3::default(), 0.02, Some(64)).is_err());
    for magnitude in [1e-200, 1.0, 1e200] {
        close(
            norm(Vec3([3.0 * magnitude, 4.0 * magnitude, 0.0])) / magnitude,
            5.0,
            2e-15,
        );
    }
}
