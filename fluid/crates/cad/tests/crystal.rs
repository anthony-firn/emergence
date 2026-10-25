//! Atomistic structures against published crystallography: coordination
//! numbers, packing fractions, nearest-neighbour distances, nanotube
//! geometry, and a cross-check against `materials::lattice`.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use cad::kits::crystal::{Crystal, Sheet, Structure};
use cad::solid::Solid;

#[test]
fn coordination_numbers_are_counted_not_assumed() {
    // Count the first-shell neighbours in a periodic supercell and compare
    // with the published coordination number of each structure.
    for structure in [Structure::Sc, Structure::Bcc, Structure::Fcc, Structure::Hcp, Structure::Diamond] {
        let a = 4e-10;
        let c = Crystal::new(structure, a, &["A"]);
        let cell = c.supercell([3, 3, 3]);
        let nn = c.nn_distance();
        // A cutoff between the first and second shell.
        let counts = cell.coordination(nn * 1.1);
        let first = counts[0];
        assert!(
            counts.iter().all(|&n| n == first),
            "{structure:?}: coordination varies across the cell: {counts:?}"
        );
        assert_eq!(
            first,
            structure.coordination(),
            "{structure:?}: counted {first} neighbours, published {}",
            structure.coordination()
        );
        // The measured minimum separation is the nearest-neighbour
        // distance the structure claims.
        let measured = cell.min_separation();
        assert!(
            (measured - nn).abs() / nn < 1e-12,
            "{structure:?}: minimum separation {measured:e} vs nn {nn:e}"
        );
    }
}

#[test]
fn packing_fractions_match_the_published_values() {
    // Ashcroft & Mermin ch. 4.
    let cases = [
        (Structure::Fcc, PI / (3.0 * 2.0f64.sqrt())),
        (Structure::Hcp, PI / (3.0 * 2.0f64.sqrt())),
        (Structure::Bcc, PI * 3.0f64.sqrt() / 8.0),
        (Structure::Sc, PI / 6.0),
        (Structure::Diamond, PI * 3.0f64.sqrt() / 16.0),
    ];
    for (structure, expect) in cases {
        let f = structure.packing_fraction().expect("monatomic");
        assert!(
            (f - expect).abs() / expect < 1e-14,
            "{structure:?}: packing {f:.6} vs published {expect:.6}"
        );
    }
    // Close packing is 0.740480, body-centred 0.680175, diamond 0.340087.
    assert!((Structure::Fcc.packing_fraction().unwrap() - 0.740480).abs() < 1e-6);
    assert!((Structure::Bcc.packing_fraction().unwrap() - 0.680175).abs() < 1e-6);
    assert!((Structure::Diamond.packing_fraction().unwrap() - 0.340087).abs() < 1e-6);
    // A compound structure has no single-sphere packing fraction.
    assert!(Structure::RockSalt.packing_fraction().is_none());
}

#[test]
fn supercells_have_the_right_count_and_density() {
    // Copper: FCC, a = 3.615 A, so 4 atoms per cell and a mass density of
    // 8.93 g/cm^3 at an atomic mass of 63.546 u.
    let cu = Crystal::copper();
    let cell = cu.supercell([4, 4, 4]);
    assert_eq!(cell.len(), 4 * 4 * 4 * 4);
    assert_eq!(cell.composition().get("Cu"), Some(&256));
    let v = cell.volume().unwrap();
    assert!((v - (4.0 * 3.615e-10f64).powi(3)).abs() / v < 1e-15);

    let u = 1.660_539_066_60e-27; // CODATA 2018
    let mut masses = BTreeMap::new();
    masses.insert("Cu".to_string(), 63.546 * u);
    let rho = cell.mass_density(&masses).unwrap();
    println!("copper: {:.1} kg/m^3", rho);
    assert!((rho - 8930.0).abs() < 60.0, "copper density {rho:.1} kg/m^3");

    // Silicon: diamond cubic, 8 atoms per cell, 2329 kg/m^3.
    let si = Crystal::silicon();
    let cell = si.supercell([3, 3, 3]);
    assert_eq!(cell.len(), 8 * 27);
    let mut masses = BTreeMap::new();
    masses.insert("Si".to_string(), 28.0855 * u);
    let rho = cell.mass_density(&masses).unwrap();
    println!("silicon: {:.1} kg/m^3", rho);
    assert!((rho - 2329.0).abs() < 10.0, "silicon density {rho:.1} kg/m^3");
    // And its nearest-neighbour distance is the 2.35 A Si-Si bond.
    assert!((si.nn_distance() * 1e10 - 2.352).abs() < 0.002, "{:.4} A", si.nn_distance() * 1e10);

    // Iron: BCC, 2 atoms per cell, 7874 kg/m^3.
    let fe = Crystal::iron();
    let mut masses = BTreeMap::new();
    masses.insert("Fe".to_string(), 55.845 * u);
    let rho = fe.supercell([3, 3, 3]).mass_density(&masses).unwrap();
    println!("iron: {:.1} kg/m^3", rho);
    assert!((rho - 7874.0).abs() < 40.0, "iron density {rho:.1} kg/m^3");
}

#[test]
fn compound_structures_have_the_right_stoichiometry() {
    let nacl = Crystal::new(Structure::RockSalt, 5.64e-10, &["Na", "Cl"]);
    let c = nacl.supercell([2, 2, 2]).composition();
    assert_eq!(c.get("Na"), c.get("Cl"), "rock salt is 1:1");
    assert_eq!(c.values().sum::<usize>(), 8 * 8);

    let gaas = Crystal::new(Structure::Zincblende, 5.653e-10, &["Ga", "As"]);
    let c = gaas.supercell([2, 2, 2]).composition();
    assert_eq!(c.get("Ga"), c.get("As"), "zincblende is 1:1");

    // Cubic perovskite ABX3.
    let p = Crystal::new(Structure::Perovskite, 3.905e-10, &["Sr", "Ti", "O"]);
    let c = p.supercell([2, 2, 2]).composition();
    assert_eq!(c["Sr"], 8);
    assert_eq!(c["Ti"], 8);
    assert_eq!(c["O"], 24, "three oxygens per formula unit");
}

#[test]
fn the_kernel_cuts_nanoparticles_out_of_a_crystal() {
    // Any solid the kernel can describe becomes a cluster: this is what
    // makes the atomistic kit part of the same CAD system.
    let cu = Crystal::copper();
    let d = 3e-9;
    let particle = cu.nanoparticle(d).expect("particle");
    // Every atom is inside the sphere and none beyond it.
    let (lo, hi) = particle.bbox();
    for k in 0..3 {
        assert!(lo[k] > -0.5 * d - 1e-12 && hi[k] < 0.5 * d + 1e-12);
    }
    // The count follows the number density to the surface-to-volume ratio.
    let expected = 4.0 / (cu.a.powi(3)) * (PI / 6.0 * d.powi(3));
    let ratio = particle.len() as f64 / expected;
    println!("3 nm Cu particle: {} atoms, bulk estimate {:.0} (ratio {ratio:.3})", particle.len(), expected);
    assert!(
        (0.9..1.1).contains(&ratio),
        "{} atoms against a bulk estimate of {expected:.0}",
        particle.len()
    );
    // The lattice is intact: the minimum separation is still the bulk
    // nearest-neighbour distance.
    assert!((particle.min_separation() - cu.nn_distance()).abs() / cu.nn_distance() < 1e-12);
    // Surface atoms have fewer neighbours than the 12 of the interior.
    let coord = particle.coordination(cu.nn_distance() * 1.1);
    assert!(coord.iter().any(|&c| c == 12), "the core must be bulk-coordinated");
    assert!(coord.iter().any(|&c| c < 12), "the surface must be under-coordinated");
    let surface = coord.iter().filter(|&&c| c < 12).count();
    println!("  {surface} of {} atoms are under-coordinated", particle.len());

    // A nanowire, and an arbitrary shape: a sphere with a hole.
    let wire = cu.nanowire(2e-9, 6e-9).expect("wire");
    let (lo, hi) = wire.bbox();
    assert!(hi[2] - lo[2] <= 6e-9 + 1e-12);
    assert!((hi[0] - lo[0]) <= 2e-9 + 1e-12);
    let shell = cu
        .carve(&Solid::sphere([0.0; 3], 1.5e-9).minus(Solid::sphere([0.0; 3], 1.0e-9)))
        .expect("shell");
    assert!(shell.len() < particle.len(), "a hollow shell holds fewer atoms");
    assert!(shell.atoms.iter().all(|a| {
        let r = (a.position[0].powi(2) + a.position[1].powi(2) + a.position[2].powi(2)).sqrt();
        r > 1.0e-9 - 1e-12
    }));

    // An unbounded solid is refused.
    assert!(cu.carve(&Solid::half_space([0.0, 0.0, 1.0], 0.0)).is_err());
}

#[test]
fn graphene_and_nanotubes_match_the_published_geometry() {
    let g = Sheet::graphene();
    // a = 2.46 A gives a 1.42 A bond.
    assert!((g.bond_length() * 1e10 - 1.4203).abs() < 1e-3, "{:.4} A", g.bond_length() * 1e10);
    let flat = g.flat([4, 4]);
    assert_eq!(flat.len(), 32, "two atoms per hexagonal cell");
    assert!(
        (flat.min_separation() - g.bond_length()).abs() / g.bond_length() < 1e-12,
        "the sheet's shortest bond must be the C-C bond"
    );
    // Every atom has three neighbours inside the sheet; count on the
    // interior where the periodic image applies cleanly.
    let coord = flat.coordination(g.bond_length() * 1.1);
    assert!(coord.iter().any(|&c| c == 3), "graphene is three-fold coordinated");

    // Nanotube geometry, against Saito et al. ch. 3.
    for (n, m) in [(10usize, 0usize), (5, 5), (10, 5), (9, 3)] {
        let d = g.tube_diameter(n, m);
        let expect = g.a * ((n * n + n * m + m * m) as f64).sqrt() / PI;
        assert!((d - expect).abs() / d < 1e-15);
        let tube = g.nanotube(n, m, 1).expect("tube");

        // Every atom sits on the cylinder of that diameter.
        let mut worst: f64 = 0.0;
        for a in &tube.atoms {
            let r = (a.position[0].powi(2) + a.position[1].powi(2)).sqrt();
            worst = worst.max((r - 0.5 * d).abs() / (0.5 * d));
        }
        assert!(worst < 1e-12, "({n},{m}): atoms off the cylinder by {worst:e}");

        // The atom count of one translational unit cell.
        let expect_count = g.tube_atoms_per_cell(n, m);
        assert_eq!(
            tube.len(),
            expect_count,
            "({n},{m}): {} atoms, published 4(n^2+nm+m^2)/d_R = {expect_count}",
            tube.len()
        );

        // And the cell length.
        let (lo, hi) = tube.bbox();
        let measured_length = g.tube_cell_length(n, m);
        assert!(
            hi[2] - lo[2] < measured_length + 1e-15,
            "({n},{m}): atoms outside the translational cell"
        );

        // Chiral angle: 0 for zigzag (n, 0), 30 degrees for armchair
        // (n, n).
        let angle = g.chiral_angle(n, m).to_degrees();
        if m == 0 {
            assert!(angle.abs() < 1e-12, "({n},0) is zigzag");
        } else if n == m {
            assert!((angle - 30.0).abs() < 1e-12, "({n},{n}) is armchair");
        } else {
            assert!(angle > 0.0 && angle < 30.0, "({n},{m}) chiral angle {angle:.3}");
        }
        println!(
            "({n},{m}): diameter {:.4} nm, {} atoms per cell, cell {:.4} nm, chiral angle {angle:.2} deg",
            d * 1e9,
            tube.len(),
            measured_length * 1e9
        );
    }

    // A (10,10) armchair tube is 1.356 nm across.
    assert!((g.tube_diameter(10, 10) * 1e9 - 1.356).abs() < 2e-3, "{:.4} nm", g.tube_diameter(10, 10) * 1e9);
    // The bond length survives the rolling to within the curvature's
    // effect on a chord: the shortest separation is the chord of the
    // 1.42 A arc, slightly less.
    let tube = g.nanotube(10, 10, 2).unwrap();
    let min = tube.min_separation();
    assert!(min < g.bond_length() && min > 0.99 * g.bond_length(), "{:e}", min);

    assert!(g.nanotube(0, 0, 1).is_err());
}

#[test]
fn a_cad_lattice_agrees_with_the_materials_crate() {
    // The same crystal through two independent paths: the CAD kit's
    // supercell, and `materials::lattice`, which computes neighbour
    // shells from the primitive vectors. Their nearest-neighbour distances
    // and first-shell counts must agree.
    use materials::lattice::Lattice;

    let a = 3.615e-10;
    let cases: [(Structure, Lattice); 4] = [
        (Structure::Fcc, Lattice::fcc(a)),
        (Structure::Bcc, Lattice::bcc(a)),
        (Structure::Sc, Lattice::sc(a)),
        (Structure::Diamond, Lattice::diamond(a)),
    ];
    for (structure, reference) in cases {
        let ours = Crystal::new(structure, a, &["A"]);
        let theirs = reference.nn_distance();
        let rel = (ours.nn_distance() - theirs).abs() / theirs;
        assert!(
            rel < 1e-12,
            "{structure:?}: CAD nn {:e} vs materials {theirs:e} ({rel:e})",
            ours.nn_distance()
        );
        // First-shell count from the materials crate's own shell analysis.
        let shells = reference.shells(0, theirs * 1.2, 1e-6);
        assert_eq!(
            shells[0].count,
            structure.coordination(),
            "{structure:?}: materials counts {} in the first shell",
            shells[0].count
        );
    }
}

#[test]
fn structures_export_as_extended_xyz() {
    let si = Crystal::silicon().supercell([2, 2, 2]);
    let xyz = si.to_xyz();
    let lines: Vec<&str> = xyz.lines().collect();
    assert_eq!(lines[0].trim(), "64", "the header counts the atoms");
    assert!(lines[1].contains("Lattice=") && lines[1].contains("Properties=species:S:1:pos:R:3"));
    assert_eq!(lines.len(), 66);
    // Positions are in angstroms, so a 2x2x2 silicon cell spans ~10.9 A.
    let coords: Vec<f64> = lines[2]
        .split_whitespace()
        .skip(1)
        .map(|v| v.parse().unwrap())
        .collect();
    assert_eq!(coords.len(), 3);
    let max: f64 = lines[2..]
        .iter()
        .flat_map(|l| l.split_whitespace().skip(1).map(|v| v.parse::<f64>().unwrap()))
        .fold(0.0, f64::max);
    assert!((max - 9.504).abs() < 0.01, "largest coordinate {max:.3} A");

    // A cluster has no lattice line.
    let cluster = Crystal::copper().nanoparticle(1.5e-9).unwrap();
    assert!(!cluster.to_xyz().lines().nth(1).unwrap().contains("Lattice="));
}
