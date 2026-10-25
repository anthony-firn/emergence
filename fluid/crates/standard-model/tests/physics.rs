use standard_model::{
    coupling,
    decay::{self, Products},
    kinematics, qed,
    species::{Charges, Domain, ALL},
    thermal,
    transport::{self, Particle, Simulation},
    FourMomentum, Rng, Species, Vec3, C,
};
use Species::*;
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}, tolerance {tol:e}");
}
fn conserve(a: FourMomentum, b: FourMomentum) {
    close(a.energy, b.energy, 2e-10 * a.energy.max(1e-3));
    for d in 0..3 {
        close(a.momentum.0[d], b.momentum.0[d], 2e-10 * a.energy.max(1e-3));
    }
}
fn sum(v: impl Iterator<Item = FourMomentum>) -> FourMomentum {
    v.fold(FourMomentum::default(), |a, b| a + b)
}
#[test]
fn complete_catalog_quantum_numbers_and_degrees_of_freedom() {
    assert_eq!(ALL.len(), 37);
    let mut fermion = 0.0;
    let mut boson = 0.0;
    for s in ALL {
        assert_eq!(s.conjugate().conjugate(), s);
        let a = s.charges();
        let b = s.conjugate().charges();
        assert_eq!(a + b, Charges::default());
        assert_eq!(s.mass(), s.conjugate().mass());
        if s.fermion() {
            fermion += s.degeneracy();
        } else {
            boson += s.degeneracy();
        }
        assert_eq!(Species::from_pdg(s.pdg_id()).unwrap().pdg_id(), s.pdg_id());
    }
    assert_eq!(fermion, 90.0);
    assert_eq!(boson, 28.0);
    assert_eq!(boson + 0.875 * fermion, 106.75);
    assert!(Species::from_pdg(39).is_none()); // no hypothetical graviton
    assert_eq!(Down.charges().electric3, -1);
    assert_eq!(Up.charges().electric3, 2);
    assert_eq!(WMinus.charges().electric3, -3);
}
#[test]
fn vacuum_rejects_unconfined_partons() {
    for s in ALL {
        let p = Particle::new(s, Vec3::default(), Vec3([1.0, 0.0, 0.0]), 0, Domain::Vacuum);
        assert_eq!(p.is_ok(), !s.needs_qcd_closure());
        assert!(Particle::new(
            s,
            Vec3::default(),
            Vec3([1.0, 0.0, 0.0]),
            0,
            Domain::Partonic
        )
        .is_ok());
    }
    assert!(decay::channels(Charm).is_empty());
}
#[test]
fn massless_flight_and_massive_speed_are_relativistic() {
    for s in [
        Photon,
        ElectronNeutrino,
        AntiMuonNeutrino,
        standard_model::species::GLUONS[7],
    ] {
        let p = Particle::new(
            s,
            Vec3::default(),
            Vec3([0.0, 0.0, 1.0]),
            1,
            Domain::Partonic,
        )
        .unwrap();
        let f = transport::fly(p, 1e-9, Domain::Partonic).unwrap();
        close(f.live[0].position().0[2], C * 1e-9, 1e-15);
    }
    let p = Particle::new(
        Electron,
        Vec3::default(),
        Vec3([Electron.mass(), 0.0, 0.0]),
        1,
        Domain::Vacuum,
    )
    .unwrap();
    let f = transport::fly(p, 1e-9, Domain::Vacuum).unwrap();
    close(f.live[0].position().0[0], C * 1e-9 / 2f64.sqrt(), 1e-15);
}
#[test]
fn two_body_threshold_shells_and_boosts() {
    let mut rng = Rng(3);
    for (m, a, b) in [
        (10.0, 3.0, 4.0),
        (80.369, 1.77693, 0.0),
        (1.0, 0.0, 0.0),
        (7.0, 3.0, 4.0),
    ] {
        let p = FourMomentum::on_shell(m, Vec3([m * 4.0, m * 2.0, m])).unwrap();
        let out = kinematics::two_body(p, m, a, b, &mut rng).unwrap();
        conserve(p, out[0] + out[1]);
        close(out[0].mass2(), a * a, 1e-9 * m * m);
        close(out[1].mass2(), b * b, 1e-9 * m * m);
        conserve(out[0].to_rest(p, m).from_rest(p, m), out[0]);
    }
    assert!(kinematics::two_body(FourMomentum::default(), 1.0, 2.0, 0.0, &mut rng).is_err());
}
#[test]
fn every_decay_channel_is_open_and_conserves_additive_charges() {
    for s in ALL {
        let table = decay::channels(s);
        assert_eq!(!table.is_empty(), s.width() > 0.0);
        if table.is_empty() {
            continue;
        }
        close(table.iter().map(|c| c.weight).sum(), 1.0, 3e-15);
        for c in table {
            let daughters = match c.products {
                Products::Two(p) => p.to_vec(),
                Products::Leptonic(p) => p.to_vec(),
                Products::Closure(_) => continue,
            };
            assert!(daughters.iter().map(|s| s.mass()).sum::<f64>() < s.mass());
            assert_eq!(
                daughters
                    .iter()
                    .fold(Charges::default(), |q, s| q + s.charges()),
                s.charges()
            );
        }
    }
}
#[test]
fn all_decay_modes_keep_energy_momentum_and_unresolved_probability() {
    let mut rng = Rng(4);
    for s in ALL.into_iter().filter(|s| s.width() > 0.0) {
        for domain in [Domain::Vacuum, Domain::Partonic] {
            let p = FourMomentum::on_shell(s.mass(), Vec3([2.0 * s.mass(), 0.3 * s.mass(), 0.0]))
                .unwrap();
            let mut unresolved = 0;
            for _ in 0..1500 {
                let d = decay::sample(s, p, domain, &mut rng).unwrap();
                conserve(
                    p,
                    sum(d
                        .particles
                        .iter()
                        .map(|v| v.1)
                        .chain(d.unresolved.iter().map(|v| v.2))),
                );
                let q = d
                    .particles
                    .iter()
                    .fold(Charges::default(), |q, p| q + p.0.charges())
                    + d.unresolved.iter().fold(Charges::default(), |q, p| q + p.3);
                assert_eq!(q, s.charges());
                unresolved += usize::from(!d.unresolved.is_empty());
            }
            if matches!(s, Tau | AntiTau | Higgs)
                || domain == Domain::Vacuum && matches!(s, WPlus | WMinus | Z | Top | AntiTop)
            {
                assert!(unresolved > 100);
            }
        }
    }
}
#[test]
fn muon_spectrum_matches_michel_and_neutrino_marginals() {
    let mut rng = Rng(6);
    let m = Muon.mass();
    let parent = FourMomentum::on_shell(m, Vec3::default()).unwrap();
    let mut means = [0.0; 3];
    let mut second = 0.0;
    let n = 30000;
    for _ in 0..n {
        let p = kinematics::leptonic(parent, m, Electron.mass(), &mut rng).unwrap();
        for i in 0..3 {
            means[i] += 2.0 * p[i].energy / m;
        }
        second += (2.0 * p[0].energy / m).powi(2);
    }
    // Integral moments: e and nu_mu 2x²(3-2x), anti_nu_e 12x²(1-x).
    close(means[0] / n as f64, 0.7, 0.005);
    close(means[1] / n as f64, 0.6, 0.005);
    close(means[2] / n as f64, 0.7, 0.005);
    close(second / n as f64, 8.0 / 15.0, 0.006);
}
#[test]
fn lifetime_distribution_and_time_dilation() {
    let tau = Muon.lifetime().unwrap();
    let n = 20000;
    let mut rest = 0;
    let mut boosted = 0;
    for i in 0..n {
        let a = Particle::new(Muon, Vec3::default(), Vec3::default(), i, Domain::Vacuum).unwrap();
        let b = Particle::new(
            Muon,
            Vec3::default(),
            Vec3([Muon.mass() * 99f64.sqrt(), 0.0, 0.0]),
            i,
            Domain::Vacuum,
        )
        .unwrap();
        rest += transport::fly(a, tau, Domain::Vacuum).unwrap().decays;
        boosted += transport::fly(b, tau, Domain::Vacuum).unwrap().decays;
    }
    close(rest as f64 / n as f64, 1.0 - (-1.0_f64).exp(), 0.009);
    close(boosted as f64 / n as f64, 1.0 - (-0.1_f64).exp(), 0.009);
}
fn fixture() -> Simulation {
    let mut sim = Simulation::new(Domain::Partonic);
    for i in 0..33 {
        for (j, s) in ALL.into_iter().enumerate() {
            let p = Particle::new(
                s,
                Vec3::default(),
                Vec3([s.mass() + 0.1, 0.05, 0.01]),
                (i * 37 + j) as u64,
                Domain::Partonic,
            )
            .unwrap();
            sim.insert(p).unwrap();
        }
    }
    sim
}
#[test]
fn cascades_conserve_totals_and_thread_count_does_not_change_histories() {
    let mut a = fixture();
    let mut b = fixture();
    let before = a.totals();
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let three = rayon::ThreadPoolBuilder::new()
        .num_threads(3)
        .build()
        .unwrap();
    one.install(|| a.advance(1e-4).unwrap());
    three.install(|| b.advance(1e-4).unwrap());
    conserve(before.0, a.totals().0);
    assert_eq!(before.1, a.totals().1);
    assert_eq!(a.decays(), b.decays());
    assert_eq!(a.particles().len(), b.particles().len());
    for (x, y) in a.particles().iter().zip(b.particles()) {
        assert_eq!(x.species(), y.species());
        assert_eq!(x.four(), y.four());
        assert_eq!(x.position(), y.position());
    }
    assert_eq!(a.unresolved().len(), b.unresolved().len());
    for (x, y) in a.unresolved().iter().zip(b.unresolved()) {
        assert_eq!(x.four, y.four);
        assert_eq!(x.time, y.time);
    }
    let state = a.totals();
    assert!(a.advance(f64::NAN).is_err());
    assert_eq!(state, a.totals());
}
#[test]
fn splitting_steps_preserves_decay_clock_and_daughter_birth_positions() {
    let p = Particle::new(
        Muon,
        Vec3::default(),
        Vec3([Muon.mass() * 3.0, 0.0, 0.0]),
        8,
        Domain::Vacuum,
    )
    .unwrap();
    let end = p.proper_time_to_decay() * 10.0;
    let mut a = Simulation::new(Domain::Vacuum);
    let mut b = Simulation::new(Domain::Vacuum);
    a.insert(p.clone()).unwrap();
    b.insert(p).unwrap();
    a.advance(end).unwrap();
    for _ in 0..10 {
        b.advance(end / 10.0).unwrap();
    }
    assert_eq!(a.decays(), 1);
    assert_eq!(b.decays(), 1);
    let mut x = a.particles().iter().collect::<Vec<_>>();
    let mut y = b.particles().iter().collect::<Vec<_>>();
    x.sort_by_key(|p| p.species().pdg_id());
    y.sort_by_key(|p| p.species().pdg_id());
    for (x, y) in x.iter().zip(y) {
        assert_eq!(x.four(), y.four());
        for d in 0..3 {
            close(x.position().0[d], y.position().0[d], 1e-10);
        }
    }
}
#[test]
fn boris_magnetic_push_preserves_energy_and_electric_push_reports_work() {
    for s in ALL.into_iter().filter(|s| s.charges().electric3 != 0) {
        let mut p = Particle::new(
            s,
            Vec3::default(),
            Vec3([s.mass(), s.mass(), 0.0]),
            12,
            Domain::Partonic,
        )
        .unwrap();
        let e = p.four().energy;
        for _ in 0..500 {
            p.push_uniform(Vec3::default(), Vec3([0.0, 0.0, 2.0]), 1e-13)
                .unwrap();
        }
        close(p.four().energy, e, 1e-12 * e);
        let before = p.four();
        let dt = 1e-10;
        let field = Vec3([1e8, 0.0, 0.0]);
        let work = p.push_uniform(field, Vec3::default(), dt).unwrap();
        close(
            p.four().momentum.0[0] - before.momentum.0[0],
            s.charges().electric3 as f64 / 3.0 * standard_model::ELEMENTARY_CHARGE * C * dt * 1e8
                / standard_model::GEV_J,
            1e-11,
        );
        close(work, p.four().energy - before.energy, 1e-15);
    }
}
#[test]
fn compton_scattering_matches_energy_shift_and_conserves_recoil() {
    let mut rng = Rng(9);
    for s in [Electron, Positron, Muon, AntiTau] {
        for ratio in [1e-8, 0.01, 1.0, 1000.0] {
            for boost in [0.0, 3.0] {
                let m = s.mass();
                let target = FourMomentum::on_shell(m, Vec3([0.0, m * boost, 0.0])).unwrap();
                let incoming = FourMomentum::on_shell(0.0, Vec3([m * ratio, 0.0, 0.0]))
                    .unwrap()
                    .from_rest(target, m);
                for _ in 0..100 {
                    let [out, recoil] = qed::compton(incoming, target, s, &mut rng).unwrap();
                    conserve(incoming + target, out + recoil);
                    let rest = out.to_rest(target, m);
                    let cosine = rest.momentum.unit().0[0];
                    close(
                        rest.energy,
                        m * ratio / (1.0 + ratio * (1.0 - cosine)),
                        2e-8 * m * ratio,
                    );
                    close(
                        recoil.energy,
                        m.hypot(recoil.momentum.length()),
                        1e-9 * recoil.energy,
                    );
                }
            }
        }
    }
    let thomson = qed::compton_cross_section(Electron, 0.0).unwrap();
    close(thomson, 6.652_458_73e-29, 1e-36);
    assert!(qed::compton_cross_section(Electron, Electron.mass()).unwrap() < thomson);
    assert!(qed::compton_cross_section(Up, 1.0).is_err());
}
#[test]
fn qed_total_cross_section_matches_integrated_differential_kernel() {
    let m = Electron.mass();
    let r = standard_model::ALPHA * standard_model::HBAR_C / m;
    for x in [1e-5, 1e-3, 0.01, 0.1, 1.0, 10.0] {
        let n = 10000;
        let h = 2.0 / n as f64;
        let mut integral = 0.0;
        for i in 0..=n {
            let c = -1.0 + i as f64 * h;
            let y = 1.0 / (1.0 + x * (1.0 - c));
            let w = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            integral += w * y * y * (y + 1.0 / y - 1.0 + c * c);
        }
        let expected = std::f64::consts::PI * r * r * integral * h / 3.0;
        let actual = qed::compton_cross_section(Electron, x * m).unwrap();
        close(actual / expected, 1.0, 2e-6);
    }
}
#[test]
fn threshold_annihilation_and_external_elastic_kernel_conserve() {
    let mut rng = Rng(13);
    let m = Electron.mass();
    let p = FourMomentum::on_shell(m, Vec3([2.0 * m, 0.0, 0.0])).unwrap();
    let out = qed::annihilate_at_threshold(p, p, &mut rng).unwrap();
    conserve(p + p, out[0] + out[1]);
    assert!(qed::annihilate_at_threshold(
        p,
        FourMomentum::on_shell(m, Vec3::default()).unwrap(),
        &mut rng
    )
    .is_err());
    let a = FourMomentum::on_shell(Up.mass(), Vec3([3.0, 0.0, 0.0])).unwrap();
    let b = FourMomentum::on_shell(Charm.mass(), Vec3([-1.0, 0.5, 0.0])).unwrap();
    let out = qed::elastic(a, b, Up.mass(), Charm.mass(), 0.2, &mut rng).unwrap();
    conserve(a + b, out[0] + out[1]);
}
#[test]
fn quantum_gas_analytic_limits_thermodynamic_identity_and_units() {
    let photon = thermal::species(Photon, 0.1).unwrap();
    let nu = thermal::species(ElectronNeutrino, 0.1).unwrap();
    close(nu.energy / photon.energy, 7.0 / 16.0, 1e-15);
    close(photon.pressure / photon.energy, 1.0 / 3.0, 1e-15);
    let high = thermal::equilibrium(1e5, Domain::Partonic).unwrap();
    close(high.g_star(1e5), 106.75, 2e-4);
    let t = Electron.mass() / 100.0;
    let gas = thermal::species(Electron, t).unwrap();
    close(gas.pressure / (gas.number * t), 1.0, 1e-9);
    close(gas.energy / (gas.number * Electron.mass()), 1.015, 3e-4);
    let t = Muon.mass() / 2.0;
    let eps = 1e-4;
    let a = thermal::species(Muon, t * (1.0 - eps)).unwrap();
    let b = thermal::species(Muon, t * (1.0 + eps)).unwrap();
    let mid = thermal::species(Muon, t).unwrap();
    close(
        ((b.pressure - a.pressure) / (2.0 * eps) - mid.pressure) / mid.energy,
        1.0,
        2e-7,
    );
    let kelvin = 300.0;
    let t = radiation::KB * kelvin / standard_model::GEV_J;
    let gas = thermal::species(Photon, t).unwrap();
    close(
        gas.energy_si() / (radiation::A_RAD * kelvin.powi(4)),
        1.0,
        3e-10,
    );
    assert!(thermal::equilibrium(0.0, Domain::Partonic).is_err());
}
#[test]
fn photon_and_source_adapters_preserve_si_energy_and_momentum() {
    let p = Particle::new(
        Photon,
        Vec3([1.0, 2.0, 3.0]),
        Vec3([0.0, 0.0, 0.002]),
        3,
        Domain::Vacuum,
    )
    .unwrap();
    let packet = coupling::photon_packet(&p, 18).unwrap();
    close(packet.photon_count(), 1.0, 3e-16);
    let tally = coupling::deposit(2, p.four()).unwrap();
    assert_eq!(tally.cell, 2);
    close(tally.energy, packet.energy, 1e-25);
    close(tally.momentum.0[2], packet.energy / C, 1e-32);
    assert!(coupling::photon_packet(
        &Particle::new(
            Electron,
            Vec3::default(),
            Vec3::default(),
            3,
            Domain::Vacuum
        )
        .unwrap(),
        0
    )
    .is_err());
}

#[test]
fn scheduled_qed_events_update_live_particles_without_resetting_decay_clocks() {
    let mut sim = Simulation::new(Domain::Vacuum);
    let photon = Particle::new(
        Photon,
        Vec3::default(),
        Vec3([0.01, 0.0, 0.0]),
        1,
        Domain::Vacuum,
    )
    .unwrap();
    let muon = Particle::new(Muon, Vec3::default(), Vec3::default(), 2, Domain::Vacuum).unwrap();
    let clock = muon.proper_time_to_decay();
    sim.insert(photon).unwrap();
    sim.insert(muon).unwrap();
    let before = sim.totals();
    sim.compton(0, 1).unwrap();
    conserve(before.0, sim.totals().0);
    assert_eq!(before.1, sim.totals().1);
    assert_eq!(clock, sim.particles()[1].proper_time_to_decay());
    let before = sim.totals();
    assert!(sim.compton(1, 0).is_err());
    assert_eq!(before, sim.totals());
    assert!(sim.compton(0, 0).is_err());
    assert!(sim.compton(0, 2).is_err());
    let mut sim = Simulation::new(Domain::Vacuum);
    for s in [Electron, Positron] {
        sim.insert(
            Particle::new(
                s,
                Vec3::default(),
                Vec3([s.mass(), 0.0, 0.0]),
                4,
                Domain::Vacuum,
            )
            .unwrap(),
        )
        .unwrap();
    }
    let before = sim.totals();
    sim.annihilate_at_threshold(0, 1).unwrap();
    conserve(before.0, sim.totals().0);
    assert_eq!(before.1, sim.totals().1);
    assert_eq!(sim.particles().len(), 2);
    assert!(sim.particles().iter().all(|p| p.species() == Photon));
}
#[test]
fn invalid_decay_inputs_are_rejected_and_impulses_preserve_clock() {
    let bad = FourMomentum {
        energy: f64::NAN,
        momentum: Vec3::default(),
    };
    let mut rng = Rng(1);
    assert!(kinematics::leptonic(bad, Muon.mass(), Electron.mass(), &mut rng).is_err());
    assert!(kinematics::two_body(bad, 1.0, 0.0, 0.0, &mut rng).is_err());
    assert!(decay::sample(Higgs, bad, Domain::Partonic, &mut rng).is_err());
    let mut p = Particle::new(Muon, Vec3::default(), Vec3::default(), 4, Domain::Vacuum).unwrap();
    let clock = p.proper_time_to_decay();
    p.set_momentum(Vec3([0.01, 0.0, 0.0])).unwrap();
    assert_eq!(clock, p.proper_time_to_decay());
    let before = p.four();
    assert!(p.set_momentum(Vec3([f64::INFINITY, 0.0, 0.0])).is_err());
    assert_eq!(before, p.four());
    assert_eq!(Higgs.twice_spin(), 0);
    assert_eq!(Electron.twice_spin(), 1);
    assert_eq!(Photon.twice_spin(), 2);
}

#[test]
fn unrepresentable_flights_reject_atomically_and_closure_handoff_retains_ledger() {
    let mut sim = Simulation::new(Domain::Vacuum);
    sim.insert(
        Particle::new(
            Photon,
            Vec3::default(),
            Vec3([1.0, 0.0, 0.0]),
            7,
            Domain::Vacuum,
        )
        .unwrap(),
    )
    .unwrap();
    let before = sim.totals();
    assert!(sim.advance(f64::MAX).is_err());
    assert_eq!(sim.time(), 0.0);
    assert_eq!(sim.totals(), before);
    let mut sim = Simulation::new(Domain::Vacuum);
    for i in 0..100 {
        sim.insert(
            Particle::new(Higgs, Vec3::default(), Vec3::default(), i, Domain::Vacuum).unwrap(),
        )
        .unwrap();
    }
    sim.advance(1e-8).unwrap();
    let before = sim.totals();
    let unresolved = sim.take_unresolved();
    assert!(!unresolved.is_empty());
    assert!(sim.unresolved().is_empty());
    let external = sum(unresolved.iter().map(|p| p.four));
    conserve(before.0, sim.totals().0 + external);
    let charges = unresolved
        .iter()
        .fold(Charges::default(), |q, p| q + p.charges);
    assert_eq!(before.1, sim.totals().1 + charges);
}
