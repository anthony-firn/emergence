//! cargo run -p standard-model --release --example relativistic
use standard_model::{
    kinematics::norm, qed, weak, Domain, FourMomentum, Particle, Rng, Simulation, Species,
    UniformField, Vec3, ALPHA, C, ELEMENTARY_CHARGE, GEV_J, HBAR_C,
};
fn main() -> Result<(), &'static str> {
    let m = Species::Electron.mass();
    let electron = FourMomentum::on_shell(m, Vec3([m * 1e12, 0.0, 0.0]))?;
    let recovered = electron.to_rest_on_shell(m, electron, m);
    let mass = electron.pair_mass(m, electron, m)?;
    assert!((mass / (2.0 * m) - 1.0).abs() < 1e-12);
    println!(
        "boost: gamma=1e12 recovered_rest_GeV={:.12e} pair_mass_GeV={:.12e}",
        recovered.energy, mass
    );

    let muon = Particle::new(
        Species::Muon,
        Vec3::default(),
        Vec3::default(),
        11,
        Domain::Vacuum,
    )?;
    let clock = muon.proper_time_to_decay();
    let mu_mass = Species::Muon.mass();
    let field = UniformField {
        electric: Vec3([
            -2.0 * mu_mass * GEV_J / (ELEMENTARY_CHARGE * C * clock),
            0.0,
            0.0,
        ]),
        magnetic: Vec3::default(),
    };
    let event = field.propagate(muon.four(), mu_mass, -1.0, clock)?;
    let mut sim = Simulation::new(Domain::Vacuum);
    sim.insert(muon)?;
    let before = sim.totals();
    let exchange = sim.advance_uniform_field(1.2 * event.lab_time, field)?;
    let after = sim.totals();
    let energy_error = (after.0.energy - before.0.energy - exchange.energy).abs() / after.0.energy;
    let momentum_error =
        norm(after.0.momentum - before.0.momentum - exchange.momentum) / after.0.energy;
    assert_eq!(sim.decays(), 1);
    assert_eq!(before.1, after.1);
    assert!(energy_error < 1e-10 && momentum_error < 1e-10);
    println!("accelerating_muon: proper_decay_s={clock:.9e} lab_decay_s={:.9e} birth_x_m={:.9e} field_energy_GeV={:.9e} energy_ledger_rel={energy_error:.3e} momentum_ledger_rel={momentum_error:.3e}",event.lab_time,event.displacement.0[0],exchange.energy);

    let beta = 0.8_f64;
    let energy = m / (1.0 - beta * beta).sqrt();
    let negative = FourMomentum::on_shell(m, Vec3([0.0, 0.0, beta * energy]))?;
    let positive = FourMomentum::on_shell(m, -negative.momentum)?;
    let s = 4.0 * energy * energy;
    let mut rng = Rng(47);
    let photons = qed::annihilate(negative, positive, &mut rng)?;
    let pair = qed::breit_wheeler(photons[0], photons[1], Species::Electron, &mut rng)?;
    let error = norm(pair[0].momentum + pair[1].momentum) / (2.0 * energy);
    assert!(error < 1e-12);
    println!("qed: beta_CM={beta} sigma_ann_m2={:.9e} sigma_pair_m2={:.9e} threshold_sigma_v_m3_s={:.9e} roundtrip_energy_rel={:.3e}",qed::annihilation_cross_section(s)?,qed::pair_cross_section(Species::Electron,s)?,qed::annihilation_sigma_v(4.0*m*m)?,((pair[0].energy+pair[1].energy)/(2.0*energy)-1.0).abs());

    let energy = 0.001;
    let neutrino = FourMomentum::on_shell(0.0, Vec3([0.0, 0.0, energy]))?;
    let target = FourMomentum::on_shell(m, Vec3::default())?;
    let couplings = weak::Couplings::tree(Species::ElectronNeutrino, 0.23112)?;
    let maximum = 2.0 * energy * energy / (m + 2.0 * energy);
    let mut sum = 0.0;
    let samples = 20000;
    for _ in 0..samples {
        let out = weak::scatter(neutrino, target, couplings, &mut rng)?;
        sum += out[1].energy - m;
    }
    println!("weak: flavor=nu_e sin2_theta_input=0.23112 E_GeV={energy} sigma_m2={:.9e} recoil_max_GeV={maximum:.9e} sampled_mean_recoil_GeV={:.9e} samples={samples}",weak::cross_section(couplings,energy)?,sum/samples as f64);
    let threshold = std::f64::consts::PI * (ALPHA * HBAR_C / m).powi(2) * C;
    assert!((qed::annihilation_sigma_v(4.0 * m * m)? / threshold - 1.0).abs() < 1e-12);
    Ok(())
}
