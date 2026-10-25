use gravity::{
    metric::{Geometry, Metric, Minkowski, Schwarzschild, WeakField},
    norm, PointMass, Vec3, C,
};
use standard_model::{
    geodesic::{self, Options},
    Domain, FourMomentum, Particle, Simulation, Species,
};
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t, "{a:e} != {b:e}, tol {t:e}");
}
#[test]
fn flat_metric_reproduces_ballistic_motion_and_proper_clock() {
    let m = Species::Muon.mass();
    let p = FourMomentum::on_shell(m, Vec3([3.0 * m, 2.0 * m, 0.0])).unwrap();
    let dt = 1e-7;
    let o = geodesic::propagate(
        &Minkowski,
        Vec3::default(),
        p,
        m,
        dt,
        f64::INFINITY,
        Options::default(),
    )
    .unwrap();
    close(
        norm(o.position - p.momentum * (C * dt / p.energy)),
        0.0,
        1e-12,
    );
    close(o.proper_time, dt * m / p.energy, 1e-22);
    assert_eq!(o.four, p);
}
#[test]
fn equivalence_principle_holds_for_slow_massive_species() {
    let r = 1e9;
    let gm = 1e20;
    let metric = WeakField(PointMass {
        gm,
        center: Vec3::default(),
        softening: 0.0,
    });
    for species in [
        Species::Electron,
        Species::Muon,
        Species::Tau,
        Species::WPlus,
        Species::Higgs,
    ] {
        let m = species.mass();
        let p = FourMomentum::on_shell(m, Vec3::default()).unwrap();
        let o = geodesic::propagate(
            &metric,
            Vec3([r, 0.0, 0.0]),
            p,
            m,
            1.0,
            f64::INFINITY,
            Options::default(),
        )
        .unwrap();
        close(
            C * o.four.momentum.0[0] / o.four.energy,
            -gm / (r * r),
            0.001,
        );
        close(o.position.0[0] - r, -0.5 * gm / (r * r), 0.001);
    }
}
#[test]
fn solar_limb_light_bending_has_general_relativistic_factor_two() {
    let gm = 1.3271244e20;
    let b = 6.957e8;
    let distance = 200.0 * b;
    let metric = Schwarzschild {
        gm,
        center: Vec3::default(),
    };
    let p = FourMomentum::on_shell(0.0, Vec3([1.0, 0.0, 0.0])).unwrap();
    let opts = Options {
        relative: 2e-12,
        position_absolute: 0.001,
        max_step: b / (10.0 * C),
        ..Options::default()
    };
    let o = geodesic::propagate(
        &metric,
        Vec3([-distance, b, 0.0]),
        p,
        0.0,
        2.0 * distance / C,
        f64::INFINITY,
        opts,
    )
    .unwrap();
    let deflection = (-o.four.momentum.0[1] / o.four.momentum.0[0]).atan();
    let expected = 4.0 * gm / (b * C * C);
    close(deflection, expected, 5e-5 * expected);
    assert_eq!(o.proper_time, 0.0);
}
#[test]
fn radial_photon_redshifts_and_conserves_stationary_killing_energy() {
    let metric = Schwarzschild {
        gm: 1.3271244e20,
        center: Vec3::default(),
    };
    let start = Vec3([6.957e8, 0.0, 0.0]);
    let p = FourMomentum::on_shell(0.0, Vec3([1.0, 0.0, 0.0])).unwrap();
    let initial = metric.geometry(start).unwrap().lapse * p.energy;
    let o = geodesic::propagate(
        &metric,
        start,
        p,
        0.0,
        100.0,
        f64::INFINITY,
        Options::default(),
    )
    .unwrap();
    let final_energy = metric.geometry(o.position).unwrap().lapse * o.four.energy;
    assert!(o.four.energy < p.energy);
    close(final_energy, initial, 5e-10 * initial);
}
struct ConstantMetric(f64);
impl Metric for ConstantMetric {
    fn geometry(&self, _: Vec3) -> Result<Geometry, &'static str> {
        Ok(Geometry {
            lapse: self.0,
            scale: 1.0,
            grad_lapse: Vec3::default(),
            grad_scale: Vec3::default(),
        })
    }
}
#[test]
fn gravitational_clock_dilation_and_decay_cascades_keep_energy() {
    let metric = ConstantMetric(0.9);
    let mu = Particle::new(
        Species::Muon,
        Vec3::default(),
        Vec3::default(),
        47,
        Domain::Vacuum,
    )
    .unwrap();
    let event = mu.proper_time_to_decay() / 0.9;
    let mut sim = Simulation::new(Domain::Vacuum);
    sim.insert(mu).unwrap();
    let initial = sim.killing_energy(&metric).unwrap();
    sim.advance_gravity(event * 0.99, &metric, Options::default())
        .unwrap();
    assert_eq!(sim.decays(), 0);
    sim.advance_gravity(event * 0.02, &metric, Options::default())
        .unwrap();
    assert_eq!(sim.decays(), 1);
    close(
        sim.killing_energy(&metric).unwrap(),
        initial,
        1e-12 * initial,
    );
    for p in sim.particles() {
        close(p.time(), event * 1.01, 1e-18);
    }
    let before = sim.totals();
    let time = sim.time();
    assert!(sim
        .advance_gravity(1.0, &ConstantMetric(f64::NAN), Options::default())
        .is_err());
    assert_eq!(sim.totals(), before);
    assert_eq!(sim.time(), time);
}
#[test]
fn geodesic_tolerances_and_step_partition_converge() {
    let metric = Schwarzschild {
        gm: 1e20,
        center: Vec3::default(),
    };
    let m = Species::Electron.mass();
    let start = Vec3([1e9, 0.0, 0.0]);
    let p = FourMomentum::on_shell(m, Vec3([0.0, m * 1e-3, 0.0])).unwrap();
    let opts = Options {
        relative: 1e-12,
        position_absolute: 1e-5,
        ..Options::default()
    };
    let direct = geodesic::propagate(&metric, start, p, m, 1000.0, f64::INFINITY, opts).unwrap();
    let half = geodesic::propagate(&metric, start, p, m, 500.0, f64::INFINITY, opts).unwrap();
    let end = geodesic::propagate(
        &metric,
        half.position,
        half.four,
        m,
        500.0,
        f64::INFINITY,
        opts,
    )
    .unwrap();
    close(norm(direct.position - end.position), 0.0, 0.1);
    close(
        norm(direct.four.momentum - end.four.momentum),
        0.0,
        1e-10 * m,
    );
    close(direct.proper_time, half.proper_time + end.proper_time, 1e-7);
}
