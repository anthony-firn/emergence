use standard_model::{
    qed, weak, Domain, FourMomentum, Particle, Rng, Simulation, Species, Vec3, ALPHA, C, HBAR_C,
};
use std::f64::consts::PI;
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t, "{a:e} != {b:e}, tol {t:e}");
}
fn four(a: FourMomentum, b: FourMomentum, t: f64) {
    close(a.energy, b.energy, t);
    for d in 0..3 {
        close(a.momentum.0[d], b.momentum.0[d], t);
    }
}
fn integrate(f: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    let n = 16000;
    let h = (b - a) / n as f64;
    (0..n).map(|i| f(a + (i as f64 + 0.5) * h)).sum::<f64>() * h
}
#[test]
fn annihilation_matches_heitler_rest_target_formula_and_threshold_flux() {
    let m = Species::Electron.mass();
    let r = ALPHA * HBAR_C / m;
    // Independent laboratory rest-target formula from the Geant4/Heitler
    // reference, rather than an integral of our implementation's kernel.
    for gamma in [1.0001_f64, 1.01, 1.5, 2.0, 10.0, 1000.0, 1e6] {
        let p = (gamma * gamma - 1.0).sqrt();
        let expected = PI * r * r / (gamma + 1.0)
            * ((gamma * gamma + 4.0 * gamma + 1.0) / (gamma * gamma - 1.0) * (gamma + p).ln()
                - (gamma + 3.0) / p);
        let s = 2.0 * m * m * (gamma + 1.0);
        close(
            qed::annihilation_cross_section(s).unwrap(),
            expected,
            5e-9 * expected,
        );
        let beta2 = 1.0 - 4.0 * m * m / s;
        close(
            qed::pair_cross_section(Species::Electron, s).unwrap(),
            2.0 * beta2 * expected,
            5e-9 * 2.0 * beta2 * expected,
        );
    }
    let limit = PI * r * r * C;
    close(
        qed::annihilation_sigma_v(4.0 * m * m).unwrap(),
        limit,
        1e-15 * limit,
    );
    for epsilon in [1e-12, 1e-9, 1e-6] {
        close(
            qed::annihilation_sigma_v(4.0 * m * m * (1.0 + epsilon)).unwrap(),
            limit,
            3.0 * epsilon * limit,
        );
    }
}
#[test]
fn finite_mass_pair_differential_kernel_integrates_to_breit_wheeler_total() {
    for species in [Species::Electron, Species::Muon, Species::Tau] {
        let m = species.mass();
        for beta in [0.0001_f64, 0.1, 0.5, 0.9, 0.99] {
            let s = 4.0 * m * m / (1.0 - beta * beta);
            let integrated = PI * ALPHA * ALPHA * HBAR_C * HBAR_C * beta / s
                * integrate(|mu| qed::pair_angular_kernel(s, m, mu).unwrap(), -1.0, 1.0);
            let total = qed::pair_cross_section(species, s).unwrap();
            close(integrated, total, 2e-6 * total);
        }
        assert_eq!(qed::pair_cross_section(species, 3.0 * m * m).unwrap(), 0.0);
    }
}
#[test]
fn pair_angular_samples_match_finite_mass_moments_and_conserve() {
    let m = Species::Electron.mass();
    let mut rng = Rng(292);
    for beta in [0.01_f64, 0.6, 0.99] {
        let energy = m / (1.0 - beta * beta).sqrt();
        let s = 4.0 * energy * energy;
        let photon = FourMomentum::on_shell(0.0, Vec3([0.0, 0.0, energy])).unwrap();
        let opposite = FourMomentum::on_shell(0.0, -photon.momentum).unwrap();
        let electron = FourMomentum::on_shell(m, Vec3([0.0, 0.0, energy * beta])).unwrap();
        let positron = FourMomentum::on_shell(m, -electron.momentum).unwrap();
        let normalization = integrate(|mu| qed::pair_angular_kernel(s, m, mu).unwrap(), -1.0, 1.0);
        let expected = integrate(
            |mu| mu * mu * qed::pair_angular_kernel(s, m, mu).unwrap(),
            -1.0,
            1.0,
        ) / normalization;
        let n = 25000;
        let mut pair_moment = 0.0;
        let mut ann_moment = 0.0;
        for _ in 0..n {
            let pair = qed::breit_wheeler(photon, opposite, Species::Electron, &mut rng).unwrap();
            four(pair[0] + pair[1], photon + opposite, 2e-12 * energy);
            for p in pair {
                p.validate_on_shell(m).unwrap();
            }
            pair_moment += (pair[0].momentum.0[2] / (energy * beta)).powi(2);
            let ann = qed::annihilate(electron, positron, &mut rng).unwrap();
            four(ann[0] + ann[1], electron + positron, 2e-12 * energy);
            ann_moment += (ann[0].momentum.0[2] / energy).powi(2);
        }
        close(pair_moment / n as f64, expected, 0.008);
        close(ann_moment / n as f64, expected, 0.008);
    }
}
#[test]
fn collision_kernels_follow_a_common_boost_and_invariant_threshold() {
    let m = Species::Electron.mass();
    let frame = FourMomentum::on_shell(1.0, Vec3([1e6, 0.0, 0.0])).unwrap();
    let e = FourMomentum::on_shell(m, Vec3([0.0, 0.0, 2.0 * m])).unwrap();
    let p = FourMomentum::on_shell(m, -e.momentum).unwrap();
    let mut rng = Rng(71);
    let mut boosted_rng = Rng(rng.0);
    let rest = qed::annihilate(e, p, &mut rng).unwrap();
    let boosted = qed::annihilate(
        e.from_rest_on_shell(m, frame, 1.0),
        p.from_rest_on_shell(m, frame, 1.0),
        &mut boosted_rng,
    )
    .unwrap();
    for i in 0..2 {
        four(
            boosted[i],
            rest[i].from_rest_on_shell(0.0, frame, 1.0),
            1e-9 * boosted[i].energy,
        );
    }
    let a = FourMomentum::on_shell(0.0, Vec3([0.0, 0.0, 3.0 * m])).unwrap();
    let b = FourMomentum::on_shell(0.0, -a.momentum).unwrap();
    let mut rng = Rng(29);
    let mut boosted_rng = Rng(rng.0);
    let rest = qed::breit_wheeler(a, b, Species::Electron, &mut rng).unwrap();
    let boosted = qed::breit_wheeler(
        a.from_rest_on_shell(0.0, frame, 1.0),
        b.from_rest_on_shell(0.0, frame, 1.0),
        Species::Electron,
        &mut boosted_rng,
    )
    .unwrap();
    for i in 0..2 {
        four(
            boosted[i],
            rest[i].from_rest_on_shell(m, frame, 1.0),
            1e-9 * boosted[i].energy,
        );
    }
    assert!(qed::breit_wheeler(a, a, Species::Electron, &mut rng).is_err());
    assert!(qed::breit_wheeler(a, b, Species::Muon, &mut rng).is_err());
    let invalid = FourMomentum {
        energy: 1.0,
        momentum: Vec3::default(),
    };
    assert!(qed::annihilate(invalid, p, &mut rng).is_err());
}
#[test]
fn neutrino_flavor_couplings_and_recoil_integrals_match_fermi_theory() {
    use Species::*;
    let m = Electron.mass();
    let angle = 0.23112;
    for species in [
        ElectronNeutrino,
        AntiElectronNeutrino,
        MuonNeutrino,
        AntiMuonNeutrino,
        TauNeutrino,
        AntiTauNeutrino,
    ] {
        let c = weak::Couplings::tree(species, angle).unwrap();
        let conjugate = weak::Couplings::tree(species.conjugate(), angle).unwrap();
        assert_eq!(c.left, conjugate.right);
        assert_eq!(c.right, conjugate.left);
        for energy in [1e-6, 1e-3, 1.0] {
            let maximum = 2.0 * energy * energy / (m + 2.0 * energy);
            let integrated = integrate(|t| weak::differential(c, energy, t).unwrap(), 0.0, maximum);
            let total = weak::cross_section(c, energy).unwrap();
            close(integrated, total, 1e-9 * total);
        }
        // E >> m limit, calculated directly from chiral coefficients.
        let energy = 10.0;
        let expected = 2.0 * weak::FERMI_GEV2.powi(2) * m * energy / PI
            * HBAR_C.powi(2)
            * (c.left * c.left + c.right * c.right / 3.0);
        close(
            weak::cross_section(c, energy).unwrap(),
            expected,
            1e-4 * expected,
        );
    }
    assert!(weak::Couplings::tree(Photon, angle).is_err());
    assert!(
        weak::cross_section(weak::Couplings::tree(ElectronNeutrino, angle).unwrap(), 1e9).is_err()
    );
}
#[test]
fn weak_recoil_samples_match_spectrum_kinematics_and_boosts() {
    let m = Species::Electron.mass();
    let electron = FourMomentum::on_shell(m, Vec3::default()).unwrap();
    let mut rng = Rng(529);
    for species in [
        Species::ElectronNeutrino,
        Species::AntiElectronNeutrino,
        Species::MuonNeutrino,
        Species::AntiMuonNeutrino,
    ] {
        let c = weak::Couplings::tree(species, 0.23112).unwrap();
        for energy in [1e-5, 1e-3, 1.0] {
            let incoming = FourMomentum::on_shell(0.0, Vec3([0.0, 0.0, energy])).unwrap();
            let maximum = 2.0 * energy * energy / (m + 2.0 * energy);
            let mean = integrate(
                |t| t * weak::differential(c, energy, t).unwrap(),
                0.0,
                maximum,
            ) / weak::cross_section(c, energy).unwrap();
            let n = 12000;
            let mut sampled = 0.0;
            for _ in 0..n {
                let out = weak::scatter(incoming, electron, c, &mut rng).unwrap();
                four(out[0] + out[1], incoming + electron, 2e-12 * (energy + m));
                let recoil = out[1].energy - m;
                assert!(recoil >= 0.0 && recoil <= maximum * (1.0 + 1e-10));
                sampled += recoil;
                let cosine = out[0].momentum.0[2] / out[0].energy;
                close(cosine, 1.0 - m * recoil / (energy * out[0].energy), 1e-9);
            }
            close(sampled / n as f64, mean, 0.012 * maximum);
            let frame = FourMomentum::on_shell(1.0, Vec3([30.0, 0.0, 0.0])).unwrap();
            let mut a = Rng(123);
            let mut b = Rng(a.0);
            let rest = weak::scatter(incoming, electron, c, &mut a).unwrap();
            let boosted = weak::scatter(
                incoming.from_rest_on_shell(0.0, frame, 1.0),
                electron.from_rest_on_shell(m, frame, 1.0),
                c,
                &mut b,
            )
            .unwrap();
            four(
                boosted[0],
                rest[0].from_rest_on_shell(0.0, frame, 1.0),
                1e-10 * (energy + m),
            );
            four(
                boosted[1],
                rest[1].from_rest_on_shell(m, frame, 1.0),
                1e-10 * (energy + m),
            );
        }
    }
}
fn particle(species: Species, p: Vec3) -> Particle {
    Particle::new(species, Vec3::default(), p, 14, Domain::Vacuum).unwrap()
}
#[test]
fn ultrarelativistic_compton_preserves_forward_recoil_and_cross_section() {
    let m = Species::Electron.mass();
    let target = FourMomentum::on_shell(m, Vec3::default()).unwrap();
    let r = ALPHA * HBAR_C / m;
    let mut rng = Rng(425);
    for x in [1e8_f64, 1e12, 1e20, 1e100] {
        let incoming = FourMomentum::on_shell(0.0, Vec3([0.0, 0.0, m * x])).unwrap();
        let sigma = qed::compton_cross_section(Species::Electron, incoming.energy).unwrap();
        let asymptotic = PI * r * r / x * ((2.0 * x).ln() + 0.5);
        close(sigma, asymptotic, 5e-7 * asymptotic);
        for _ in 0..2000 {
            let out = qed::compton(incoming, target, Species::Electron, &mut rng).unwrap();
            let transfer = incoming.dot_on_shell(0.0, out[0], 0.0).unwrap();
            let expected = m * (incoming.energy - out[0].energy);
            close(transfer, expected, 2e-10 * expected);
            four(incoming + target, out[0] + out[1], 2e-15 * incoming.energy);
        }
    }
}
#[test]
fn scheduled_interactions_preserve_ledger_birth_clocks_and_atomic_rejection() {
    let mut sim = Simulation::new(Domain::Vacuum);
    sim.insert(particle(Species::Electron, Vec3([0.1, 0.0, 0.0])))
        .unwrap();
    sim.insert(particle(Species::Positron, Vec3([-0.1, 0.0, 0.0])))
        .unwrap();
    sim.advance(1e-9).unwrap();
    let initial = sim.totals();
    sim.annihilate(0, 1).unwrap();
    four(sim.totals().0, initial.0, 1e-12);
    assert_eq!(sim.totals().1, initial.1);
    for p in sim.particles() {
        assert_eq!(p.species(), Species::Photon);
        assert_eq!(p.time(), sim.time());
    }
    let before_rejection = sim.totals();
    assert!(sim.pair_produce(0, 1, Species::Tau).is_err());
    assert_eq!(sim.totals(), before_rejection);
    sim.pair_produce(0, 1, Species::Electron).unwrap();
    four(sim.totals().0, initial.0, 1e-12);
    assert_eq!(sim.totals().1, initial.1);
    assert_eq!(sim.particles()[0].species(), Species::Electron);
    assert_eq!(sim.particles()[1].species(), Species::Positron);
    let mut weak_sim = Simulation::new(Domain::Vacuum);
    weak_sim
        .insert(particle(Species::ElectronNeutrino, Vec3([0.001, 0.0, 0.0])))
        .unwrap();
    weak_sim
        .insert(particle(Species::Electron, Vec3::default()))
        .unwrap();
    let before = weak_sim.totals();
    let c = weak::Couplings::tree(Species::ElectronNeutrino, 0.23112).unwrap();
    weak_sim.scatter_neutrino(0, 1, c).unwrap();
    four(before.0, weak_sim.totals().0, 1e-14);
    assert_eq!(before.1, weak_sim.totals().1);
    let before = weak_sim.totals();
    assert!(weak_sim.scatter_neutrino(1, 0, c).is_err());
    assert_eq!(before, weak_sim.totals());
    let mut unstable = Simulation::new(Domain::Vacuum);
    unstable
        .insert(particle(Species::Photon, Vec3([2.0, 0.0, 0.0])))
        .unwrap();
    unstable
        .insert(particle(Species::Photon, Vec3([-2.0, 0.0, 0.0])))
        .unwrap();
    unstable.pair_produce(0, 1, Species::Tau).unwrap();
    for p in unstable.particles() {
        assert!(p.proper_time_to_decay().is_finite());
        assert_eq!(p.proper_age(), 0.0);
    }
}
