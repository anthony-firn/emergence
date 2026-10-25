use standard_model::{
    fields::UniformField,
    kinematics::{self, norm},
    transport, Domain, FourMomentum, Particle, Rng, Simulation, Species, Vec3, C,
    ELEMENTARY_CHARGE, GEV_J,
};
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}, tol {tol:e}");
}
fn vector(a: Vec3, b: Vec3, tol: f64) {
    for d in 0..3 {
        close(a.0[d], b.0[d], tol);
    }
}
fn four(a: FourMomentum, b: FourMomentum, tol: f64) {
    close(a.energy, b.energy, tol);
    vector(a.momentum, b.momentum, tol);
}
#[test]
fn light_cone_boosts_and_pair_invariants_survive_large_gamma() {
    let m = Species::Electron.mass();
    for gamma in [1.0_f64, 10.0, 1e4, 1e8, 1e12] {
        let frame =
            FourMomentum::on_shell(m, Vec3([m * (gamma * gamma - 1.0).sqrt(), 0.0, 0.0])).unwrap();
        let rest = FourMomentum::on_shell(m, Vec3::default()).unwrap();
        let p = rest.from_rest_on_shell(m, frame, m);
        four(p.to_rest_on_shell(m, frame, m), rest, 2e-15 * m);
        close(p.pair_mass(m, p, m).unwrap(), 2.0 * m, 2e-15 * m);
        let photon = FourMomentum::on_shell(0.0, Vec3([m, 0.0, 0.0])).unwrap();
        let boosted = photon.from_rest_on_shell(0.0, frame, m);
        four(boosted.to_rest_on_shell(0.0, frame, m), photon, 2e-15 * m);
        close(
            p.dot_on_shell(m, boosted, 0.0).unwrap(),
            m * m,
            3e-15 * m * m,
        );
    }
    let a = FourMomentum::on_shell(0.0, Vec3([1e8, 0.0, 0.0])).unwrap();
    let b = FourMomentum::on_shell(0.0, Vec3([1e8, 1.0, 0.0])).unwrap();
    close(a.pair_mass(0.0, b, 0.0).unwrap(), 1.0, 1e-15);
}
#[test]
fn constant_electric_field_worldline_matches_lab_solution_and_proper_time() {
    let m = Species::Positron.mass();
    let dt = 1e-8;
    let field = UniformField {
        electric: Vec3([1e6, 0.0, 0.0]),
        magnetic: Vec3::default(),
    };
    let p = FourMomentum::on_shell(m, Vec3::default()).unwrap();
    let (orbit, event) = field.until(p, m, 1.0, dt, f64::INFINITY).unwrap();
    assert!(!event);
    let rate = ELEMENTARY_CHARGE * C * 1e6 / GEV_J;
    let momentum = rate * dt;
    let energy = m.hypot(momentum);
    close(orbit.four.momentum.0[0], momentum, 2e-14 * momentum);
    close(orbit.four.energy, energy, 2e-14 * energy);
    close(
        orbit.displacement.0[0],
        C / rate * (momentum * momentum / (energy + m)),
        2e-13,
    );
    close(
        orbit.proper_time,
        m / rate * (momentum / m).asinh(),
        2e-14 * dt,
    );
    close(orbit.lab_time, dt, 2e-14 * dt);
    let chi = field.quantum_chi(p, m, 1.0).unwrap();
    close(chi, 1e6 / 1.323_285_474e18, 1e-20);
}
#[test]
fn electric_deceleration_can_return_an_ultrarelativistic_particle_to_rest() {
    let m = Species::Electron.mass();
    let gamma = 1e8_f64;
    let p = FourMomentum::on_shell(m, Vec3([m * (gamma * gamma - 1.0).sqrt(), 0.0, 0.0])).unwrap();
    let field = UniformField {
        electric: Vec3([1e8, 0.0, 0.0]),
        magnetic: Vec3::default(),
    };
    let tau =
        (gamma + (gamma * gamma - 1.0).sqrt()).ln() * m * GEV_J / (ELEMENTARY_CHARGE * C * 1e8);
    let orbit = field.propagate(p, m, -1.0, tau).unwrap();
    close(orbit.four.energy, m, 1e-12 * m);
    close(orbit.four.momentum.0[0], 0.0, 1e-12 * m);
    assert!(orbit.lab_time > tau);
}
#[test]
fn gyro_orbit_is_exact_without_resolving_a_period() {
    let m = Species::Electron.mass();
    let p = FourMomentum::on_shell(m, Vec3([2.0 * m, 0.0, m])).unwrap();
    let field = UniformField {
        electric: Vec3::default(),
        magnetic: Vec3([0.0, 0.0, 1.0]),
    };
    let omega = -ELEMENTARY_CHARGE * C * C / (m * GEV_J);
    for turns in [0.25, 1.0, 10.0, 1000.0] {
        let theta = -std::f64::consts::TAU * turns;
        let tau = theta / omega;
        let orbit = field.propagate(p, m, -1.0, tau).unwrap();
        close(orbit.four.energy, p.energy, 1e-14 * p.energy);
        vector(
            orbit.four.momentum,
            Vec3([2.0 * m * theta.cos(), -2.0 * m * theta.sin(), m]),
            // Constructing tau from theta incurs phase error proportional to turns.
            8.0 * f64::EPSILON * theta.abs().max(1.0) * m,
        );
        vector(
            orbit.displacement,
            Vec3([
                C * 2.0 * theta.sin() / omega,
                C * 2.0 * (theta.cos() - 1.0) / omega,
                C * tau,
            ]),
            1e-12,
        );
        close(orbit.lab_time, tau * p.energy / m, 1e-12 * tau);
    }
}
#[test]
fn parallel_fields_remain_invariant_at_ultrarelativistic_gamma() {
    let field = UniformField {
        electric: Vec3([1e6, 0.0, 0.0]),
        magnetic: Vec3([0.5, 0.0, 0.0]),
    };
    for gamma in [1e4, 1e8, 1e12] {
        let frame = FourMomentum::on_shell(1.0, Vec3([gamma, 0.0, 0.0])).unwrap();
        let rest = field.to_rest(frame, 1.0).unwrap();
        assert_eq!(rest.electric, field.electric);
        assert_eq!(rest.magnetic, field.magnetic);
    }
}
#[test]
fn crossed_fields_matrix_solution_is_covariant_and_composes() {
    let m = Species::Muon.mass();
    let p = FourMomentum::on_shell(m, Vec3([m, 0.2 * m, 0.3 * m])).unwrap();
    let field = UniformField {
        electric: Vec3([3e7, -2e7, 1e7]),
        magnetic: Vec3([0.1, -0.2, 0.3]),
    };
    let tau = 2e-9;
    let direct = field.propagate(p, m, -1.0, tau).unwrap();
    close(direct.four.mass2(), m * m, 1e-11 * m * m);
    let half = field.propagate(p, m, -1.0, tau / 2.0).unwrap();
    let end = field.propagate(half.four, m, -1.0, tau / 2.0).unwrap();
    four(direct.four, end.four, 1e-12);
    vector(
        direct.displacement,
        half.displacement + end.displacement,
        1e-12,
    );
    close(direct.lab_time, half.lab_time + end.lab_time, 1e-20);
    let frame = FourMomentum::on_shell(m, Vec3([0.4 * m, 0.2 * m, -0.1 * m])).unwrap();
    let transformed = field.to_rest(frame, m).unwrap();
    let transformed_p = p.to_rest_on_shell(m, frame, m);
    let transformed_orbit = transformed.propagate(transformed_p, m, -1.0, tau).unwrap();
    four(
        transformed_orbit.four,
        direct.four.to_rest_on_shell(m, frame, m),
        1e-12,
    );
    let displacement = FourMomentum {
        energy: C * direct.lab_time,
        momentum: direct.displacement,
    }
    .to_rest(frame, m);
    close(C * transformed_orbit.lab_time, displacement.energy, 1e-12);
    vector(transformed_orbit.displacement, displacement.momentum, 1e-12);
    close(
        norm(transformed.electric).powi(2) - C * C * norm(transformed.magnetic).powi(2),
        norm(field.electric).powi(2) - C * C * norm(field.magnetic).powi(2),
        100.0,
    );
    close(
        transformed.electric.dot(transformed.magnetic),
        field.electric.dot(field.magnetic),
        1e-8,
    );
}
#[test]
fn accelerating_muon_decay_occurs_at_proper_clock_and_preserves_field_ledger() {
    let mu = Particle::new(
        Species::Muon,
        Vec3::default(),
        Vec3::default(),
        11,
        Domain::Vacuum,
    )
    .unwrap();
    let clock = mu.proper_time_to_decay();
    let m = Species::Muon.mass();
    // Set rapidity at the sampled decay to 2, safely within classical field limits.
    let field = UniformField {
        electric: Vec3([-2.0 * m * GEV_J / (ELEMENTARY_CHARGE * C * clock), 0.0, 0.0]),
        magnetic: Vec3::default(),
    };
    let expected = field.propagate(mu.four(), m, -1.0, clock).unwrap();
    assert!(expected.lab_time > 1.5 * clock);
    let before =
        transport::fly_uniform_field(mu.clone(), expected.lab_time * 0.9, Domain::Vacuum, field)
            .unwrap();
    assert_eq!(before.decays, 0);
    assert_eq!(before.live[0].species(), Species::Muon);
    let at =
        transport::fly_uniform_field(mu.clone(), expected.lab_time, Domain::Vacuum, field).unwrap();
    assert_eq!(at.decays, 1);
    assert_eq!(at.live.len(), 3);
    for p in &at.live {
        vector(p.position(), expected.displacement, 1e-7);
        close(p.time(), expected.lab_time, 1e-14 * expected.lab_time);
    }
    let mut sim = Simulation::new(Domain::Vacuum);
    sim.insert(mu).unwrap();
    let initial = sim.totals();
    let exchange = sim
        .advance_uniform_field(expected.lab_time * 1.2, field)
        .unwrap();
    let final_state = sim.totals();
    four(
        initial.0 + exchange,
        final_state.0,
        1e-10 * final_state.0.energy,
    );
    assert_eq!(initial.1, final_state.1);
    assert_eq!(sim.decays(), 1);
}
#[test]
fn field_steps_partition_correctly_and_failed_steps_are_atomic() {
    let mut a = Simulation::new(Domain::Vacuum);
    let mut b = Simulation::new(Domain::Vacuum);
    let p = Particle::new(
        Species::Positron,
        Vec3::default(),
        Vec3([0.001, 0.0001, 0.0]),
        3,
        Domain::Vacuum,
    )
    .unwrap();
    a.insert(p.clone()).unwrap();
    b.insert(p).unwrap();
    let field = UniformField {
        electric: Vec3([1e6, 2e6, 0.0]),
        magnetic: Vec3([0.0, 0.0, 0.005]),
    };
    let dt = 1e-9;
    let work = a.advance_uniform_field(dt, field).unwrap();
    let mut parts = FourMomentum::default();
    for _ in 0..20 {
        parts = parts + b.advance_uniform_field(dt / 20.0, field).unwrap();
    }
    four(a.particles()[0].four(), b.particles()[0].four(), 1e-13);
    four(work, parts, 1e-13);
    vector(
        a.particles()[0].position(),
        b.particles()[0].position(),
        1e-11,
    );
    close(
        a.particles()[0].proper_age(),
        b.particles()[0].proper_age(),
        1e-20,
    );
    let before = a.totals();
    let bad = UniformField {
        electric: Vec3([f64::NAN, 0.0, 0.0]),
        magnetic: Vec3::default(),
    };
    assert!(a.advance_uniform_field(dt, bad).is_err());
    assert_eq!(before, a.totals());
    close(a.time(), dt, 1e-25);
}
#[test]
fn high_gamma_decay_kinematics_conserve_without_losing_invariant_pair_mass() {
    let mut rng = Rng(5);
    let m = Species::Z.mass();
    let p = FourMomentum::on_shell(m, Vec3([m * 1e8, 0.0, 0.0])).unwrap();
    let pair = kinematics::two_body(p, m, 0.0, 0.0, &mut rng).unwrap();
    four(pair[0] + pair[1], p, 3e-6 * p.energy);
    close(pair[0].pair_mass(0.0, pair[1], 0.0).unwrap(), m, 2e-6 * m);
}
