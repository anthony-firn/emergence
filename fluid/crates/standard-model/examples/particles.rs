//! cargo run -p standard-model --release --example particles -- [copies=1000] [vacuum|partonic]
use standard_model::{
    coupling, qed,
    species::{Domain, ALL},
    thermal, FourMomentum, Particle, Rng, Simulation, Species, Vec3,
};
use std::time::Instant;
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let copies: usize = args
        .get(1)
        .map(|s| s.parse().expect("integer copies"))
        .unwrap_or(1000);
    let domain = match args.get(2).map(String::as_str).unwrap_or("partonic") {
        "vacuum" => Domain::Vacuum,
        "partonic" => Domain::Partonic,
        _ => panic!("expected vacuum|partonic"),
    };
    assert!(copies > 0);
    println!("PDG 2025 species (GeV, lifetimes in s), domain={domain:?}");
    for s in ALL {
        println!(
            "{s:?}: id={} mass={:.10e} charge={:+.3} g={} lifetime={:?} QCD_closure={}",
            s.pdg_id(),
            s.mass(),
            s.charges().electric3 as f64 / 3.0,
            s.degeneracy(),
            s.lifetime(),
            s.needs_qcd_closure()
        );
    }
    let mut sim = Simulation::new(domain);
    for i in 0..copies {
        for (j, s) in ALL.into_iter().enumerate() {
            if s.allowed(domain) {
                sim.insert(
                    Particle::new(
                        s,
                        Vec3::default(),
                        Vec3([s.mass() + 0.1, 0.05, 0.01]),
                        (i * 37 + j) as u64,
                        domain,
                    )
                    .unwrap(),
                )
                .unwrap();
            }
        }
    }
    let initial = sim.particles().len();
    let before = sim.totals();
    let start = Instant::now();
    sim.advance(1e-4).unwrap();
    let seconds = start.elapsed().as_secs_f64();
    let after = sim.totals();
    let error = (after.0.energy - before.0.energy).abs() / before.0.energy;
    let dp = (after.0.momentum - before.0.momentum).length() / before.0.energy;
    assert!(error < 1e-10 && dp < 1e-10);
    assert_eq!(before.1, after.1);
    println!("cascade: initial={initial} live={} unresolved={} decays={} wall_s={seconds:.6} energy_rel={error:.3e} momentum_rel={dp:.3e}",sim.particles().len(),sim.unresolved().len(),sim.decays());
    println!(
        "Unresolved records require hadronic or off-shell closure; they retain energy and charges."
    );
    let photon = sim
        .particles()
        .iter()
        .find(|p| p.species() == Species::Photon)
        .unwrap();
    let packet = coupling::photon_packet(photon, 3).unwrap();
    println!(
        "photon adapter: wavelength_m={:.6e} packet_J={:.6e} photon_count={:.6}",
        packet.wavelength,
        packet.energy,
        packet.photon_count()
    );
    let mut rng = Rng(14);
    let m = Species::Electron.mass();
    let incoming = FourMomentum::on_shell(0.0, Vec3([m, 0.0, 0.0])).unwrap();
    let target = FourMomentum::on_shell(m, Vec3::default()).unwrap();
    let out = qed::compton(incoming, target, Species::Electron, &mut rng).unwrap();
    println!(
        "Compton: in_GeV={:.6e} out_GeV={:.6e} sigma_m2={:.6e}",
        incoming.energy,
        out[0].energy,
        qed::compton_cross_section(Species::Electron, incoming.energy).unwrap()
    );
    for t in [1e-4, 0.1, 1.0, 100.0, 1e5] {
        let gas = thermal::equilibrium(t, domain).unwrap();
        println!(
            "ideal quantum gas: kBT_GeV={t:.3e} g_star={:.6} pressure_Pa={:.6e}",
            gas.g_star(t),
            gas.pressure_si()
        );
    }
}
