//! Benchmarks for the atomistic rung of the materials ladder (`AGENTS.md`,
//! `analysis/emergence_ladder.md`). Every test checks an exact identity,
//! a published lattice sum, or an independent finite difference; no
//! tolerance is widened to make a model agree with itself.

use materials::continuum::{christoffel_velocities, dispersion_error, k_norm};
use materials::elastic::{
    elastic_tensor, elastic_tensor_extrapolated, elastic_tensor_fd, elastic_tensor_frozen,
};
use materials::lattice::{Crystal, Lattice, IDEAL_C_OVER_A};
use materials::phonon::{
    chain_frequency, chain_lattice, fcc_path, sample_path, two_pi_over, Harmonics,
};
use materials::potential::{
    lattice_sum, lattice_sum_corrected, relax, relax_cutoff, Harmonic, LennardJones, Morse,
    PairPotential,
};

/// Exact fcc Lennard-Jones lattice sums (standard values, e.g. Kittel /
/// Ashcroft–Mermin tables of the fcc lattice sums).
const A6_FCC: f64 = 14.45392;
const A12_FCC: f64 = 12.13188;

fn close(actual: f64, expected: f64, tol: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tol,
        "{what}: {actual:.12e} vs {expected:.12e}, tolerance {tol:.3e}"
    );
}

fn rel_close(actual: f64, expected: f64, tol: f64, what: &str) {
    let scale = expected.abs().max(1e-300);
    assert!(
        (actual - expected).abs() / scale <= tol,
        "{what}: {actual:.12e} vs {expected:.12e}, rel tolerance {tol:.3e}"
    );
}

#[test]
fn shell_counts_are_exact() {
    // fcc: 12, 6, 24, 12 (Ashcroft–Mermin neighbour shells).
    let counts = |lat: &Lattice, rmax: f64| -> Vec<usize> {
        lat.shells(0, rmax, 1e-6)
            .iter()
            .map(|s| s.count)
            .take(4)
            .collect()
    };
    assert_eq!(counts(&Lattice::fcc(1.0), 1.45), vec![12, 6, 24, 12]);
    // bcc: 8, 6, 12.
    assert_eq!(counts(&Lattice::bcc(1.0), 1.45), vec![8, 6, 12]);
    // diamond: 4, 12.
    let dia = counts(&Lattice::diamond(1.0), 0.75);
    assert_eq!(dia, vec![4, 12]);
    // Ideal hcp is close packed: 12 first neighbours at the lattice constant.
    let hcp = Lattice::hcp(1.0);
    let first = hcp.shells(0, 1.001, 1e-6);
    assert_eq!(first[0].count, 12);
    close(first[0].radius, 1.0, 1e-12, "ideal hcp nn");

    // Sc / Bravais nearest-neighbour distances.
    close(Lattice::sc(1.0).nn_distance(), 1.0, 1e-12, "sc nn");
    close(
        Lattice::fcc(1.0).nn_distance(),
        0.5_f64.sqrt(),
        1e-12,
        "fcc nn",
    );
    close(
        Lattice::bcc(1.0).nn_distance(),
        3.0_f64.sqrt() / 2.0,
        1e-12,
        "bcc nn",
    );
    close(
        Lattice::diamond(1.0).nn_distance(),
        3.0_f64.sqrt() / 4.0,
        1e-12,
        "diamond nn",
    );

    // An explicit, non-ideal c/a is honoured.
    let squashed = Lattice::hcp_ca(1.0, 1.5);
    close(
        squashed.nn_distance(),
        (1.0_f64 / 3.0 + 0.5625).sqrt(),
        1e-12,
        "hcp c/a=1.5 nn",
    );
    assert!(Lattice::new(&Crystal::Hcp { c_over_a: 1.5 }, 1.0)
        .label
        .contains("1.5"));
    close(IDEAL_C_OVER_A, (8.0 / 3.0_f64).sqrt(), 1e-15, "ideal c/a");
}

#[test]
fn lj_lattice_sums_converge_to_the_exact_values() {
    let fcc = Lattice::fcc(1.0);
    // The r^-12 sum converges at the first few shells.
    close(lattice_sum(&fcc, 12.0, 3.0), A12_FCC, 5e-5, "A12 rmax=3");
    // The raw r^-6 sum converges slowly; the homogeneous tail accelerates it.
    let mut prev = f64::INFINITY;
    for rmax in [4.0, 6.0, 8.0, 10.0] {
        let corrected = lattice_sum_corrected(&fcc, 6.0, rmax);
        let err = (corrected - A6_FCC).abs();
        assert!(err < prev, "corrected A6 not converging at rmax={rmax}");
        prev = err;
        println!("A6 corrected rmax={rmax}: {corrected:.8} err={err:.3e}");
    }
    close(
        lattice_sum_corrected(&fcc, 6.0, 10.0),
        A6_FCC,
        5e-5,
        "A6 rmax=10",
    );
}

#[test]
fn lj_fcc_relax_reproduces_the_exact_cohesive_energy() {
    let lj = LennardJones::new(1.0, 1.0);
    let fcc = Lattice::fcc(1.0);
    let eq = relax(&fcc, &lj, 10.0, 0.9, 1.8);
    // Exact minimum of the infinite-sum energy,
    // r0 = (2 A12 / A6)^(1/6) sigma, E0 = 2(A12/r0^12 - A6/r0^6) eps.
    let r0 = (2.0 * A12_FCC / A6_FCC).powf(1.0 / 6.0);
    let e_exact = 2.0 * (A12_FCC / r0.powi(12) - A6_FCC / r0.powi(6));
    println!(
        "fcc LJ r0={:.9} (exact {:.9}), E={:.9} (exact {:.9})",
        eq.nn_distance, r0, eq.energy_per_atom, e_exact
    );
    // The tail-consistent relaxation reproduces the exact lattice-sum
    // minimum, not just the truncated-sum one.
    rel_close(eq.nn_distance, r0, 1e-6, "fcc r0");
    rel_close(eq.energy_per_atom, e_exact, 1e-6, "fcc cohesive energy");
}

#[test]
fn lj_fcc_converged_elastic_constants_match_the_infinite_sum_reference() {
    // The raw 8-sigma-cutoff constants (C11 = 104.927, C12 = C44 =
    // 60.019 eps/sigma^3) carry an O(rc^-3) truncation tail. Evaluate the
    // Cauchy-Born stiffness at two cutoffs and eliminate the tail to get the
    // converged (infinite-sum) values at the exact r0. Independent
    // reference (infinite fcc lattice sums at zero pressure, +/-0.02):
    // C11 ~ 105.19, C12 = C44 ~ 60.19 eps/sigma^3.
    let lj = LennardJones::new(1.0, 1.0);
    let r0 = (2.0 * A12_FCC / A6_FCC).powf(1.0 / 6.0);
    let lat = Lattice::fcc(r0 * 2.0_f64.sqrt());
    let c = elastic_tensor_extrapolated(&lat, &lj, 12.0, 20.0);
    println!(
        "converged fcc LJ C11={:.5} C12={:.5} C44={:.5} (ref 105.19/60.19/60.19)",
        c[0][0], c[0][1], c[3][3]
    );
    close(c[0][0], 105.19, 0.02, "converged C11");
    close(c[0][1], 60.19, 0.02, "converged C12");
    close(c[3][3], 60.19, 0.02, "converged C44");
    // Cauchy relation of a central-force crystal at zero pressure is
    // preserved by the extrapolation.
    rel_close(c[0][1], c[3][3], 1e-4, "converged Cauchy C12=C44");
}

#[test]
fn lj_phase_order_is_hcp_then_fcc_then_bcc() {
    let lj = LennardJones::new(1.0, 1.0);
    let fcc = relax(&Lattice::fcc(1.0), &lj, 10.0, 0.9, 1.8);
    let hcp = relax(&Lattice::hcp(1.0), &lj, 10.0, 0.9, 1.8);
    let bcc = relax(&Lattice::bcc(1.0), &lj, 10.0, 0.9, 2.0);
    println!(
        "fcc E={:.6}, hcp E={:.6}, bcc E={:.6}",
        fcc.energy_per_atom, hcp.energy_per_atom, bcc.energy_per_atom
    );
    assert!(
        hcp.energy_per_atom < fcc.energy_per_atom,
        "hcp should lie below fcc"
    );
    assert!(
        bcc.energy_per_atom > fcc.energy_per_atom,
        "bcc should lie above fcc"
    );
    // The known fcc/hcp difference is ~0.01%.
    let rel = (fcc.energy_per_atom - hcp.energy_per_atom).abs() / fcc.energy_per_atom.abs();
    assert!(
        (1e-5..2e-4).contains(&rel),
        "fcc-hcp relative difference {rel:e} outside the expected ~0.01%"
    );
}

#[test]
fn morse_and_harmonic_derivatives_are_exact() {
    // Independent central differences at well-separated orders.
    for (pot, r0) in [
        (&Morse::new(1.0, 2.0, 1.1) as &dyn PairPotential, 1.1),
        (&Harmonic { k: 3.0, r0: 1.2 }, 1.2),
    ] {
        for &r in &[0.9, 1.15, 1.4, 2.0] {
            let h1 = 1e-6 * r;
            let h2 = 1e-4 * r;
            let d1 = (pot.energy(r + h1) - pot.energy(r - h1)) / (2.0 * h1);
            let d2 = (pot.energy(r + h2) - 2.0 * pot.energy(r) + pot.energy(r - h2)) / (h2 * h2);
            close(pot.d1(r), d1, 1e-5 * pot.d2(r0).abs().max(1.0), "d1");
            close(pot.d2(r), d2, 1e-3 * pot.d2(r0).abs().max(1.0), "d2");
        }
    }
}

#[test]
fn cauchy_relation_holds_for_pair_potentials_at_zero_pressure() {
    let lj = LennardJones::new(1.0, 1.0);
    let raw = relax_cutoff(&Lattice::fcc(1.0), &lj, 8.0, 0.9, 1.8);
    let lat = Lattice::fcc(raw.scale);
    let c = elastic_tensor(&lat, &lj, 8.0);
    println!(
        "fcc LJ C11={:.5} C12={:.5} C44={:.5} (eps/sigma^3)",
        c[0][0], c[0][1], c[3][3]
    );

    // Cubic symmetry: C11=C22=C33, C12=C13=C23, C44=C55=C66, and the
    // 1-4/1-5/... couplings vanish.
    close(c[0][0], c[1][1], 1e-9, "C11=C22");
    close(c[1][1], c[2][2], 1e-9, "C22=C33");
    close(c[0][1], c[0][2], 1e-9, "C12=C13");
    close(c[0][1], c[1][2], 1e-9, "C12=C23");
    close(c[3][3], c[4][4], 1e-9, "C44=C55");
    close(c[4][4], c[5][5], 1e-9, "C55=C66");
    for i in 0..6 {
        for j in 0..6 {
            let is_diag_block = matches!(
                (i, j),
                (0, 0)
                    | (1, 1)
                    | (2, 2)
                    | (0, 1)
                    | (1, 0)
                    | (0, 2)
                    | (2, 0)
                    | (1, 2)
                    | (2, 1)
                    | (3, 3)
                    | (4, 4)
                    | (5, 5)
            );
            if !is_diag_block {
                close(c[i][j], 0.0, 1e-9, "cubic zero entry");
            }
        }
    }
    // Exact Cauchy relation for central forces on a Bravais lattice.
    rel_close(c[0][1], c[3][3], 1e-6, "Cauchy C12=C44");

    // Independent finite difference of the relaxed energy under strain.
    let cfd = elastic_tensor_fd(&lat, &lj, 8.0, 1e-4);
    for (i, j) in [(0, 0), (0, 1), (3, 3)] {
        rel_close(c[i][j], cfd[i][j], 1e-4, "analytic vs FD");
    }
}

#[test]
fn diamond_floppy_shear_is_cured_by_internal_relaxation() {
    // Diamond with only nearest-neighbour central springs is isostatic
    // (z = 4 in d = 3): the affine (frozen) C44 is finite but sublattice
    // relaxation removes it entirely, and the finite difference confirms
    // the relaxed value.
    let a = 4.0 / 3.0_f64.sqrt();
    let lat = Lattice::diamond(a);
    let spring = Harmonic { k: 2.0, r0: 1.0 };
    let frozen = elastic_tensor_frozen(&lat, &spring, 1.2);
    let relaxed = elastic_tensor(&lat, &spring, 1.2);
    let fd = elastic_tensor_fd(&lat, &spring, 1.2, 1e-4);
    println!(
        "diamond spring C44 frozen={:.6} relaxed={:.6} fd={:.6}",
        frozen[3][3], relaxed[3][3], fd[3][3]
    );
    assert!(frozen[3][3] > 0.2, "frozen C44 should be finite");
    close(relaxed[3][3], 0.0, 1e-9, "relaxed C44 floppy mode");
    close(fd[3][3], 0.0, 1e-8, "fd C44 floppy mode");
}

#[test]
fn hcp_internal_relaxation_matches_the_finite_difference() {
    // Nearest-neighbour central springs at the ideal lattice constant: every
    // included bond sits at its relaxed length, so the configuration is at
    // zero stress and the relaxed stiffness must equal the finite difference.
    let a = 1.0;
    let lat = Lattice::hcp(a);
    let spring = Harmonic { k: 2.0, r0: a };
    let rmax = 1.2 * a; // first shell (12) only; second shell is at sqrt(2) a
    let frozen = elastic_tensor_frozen(&lat, &spring, rmax);
    let relaxed = elastic_tensor(&lat, &spring, rmax);
    let fd = elastic_tensor_fd(&lat, &spring, rmax, 1e-4);
    println!(
        "hcp spring C11 frozen={:.6} relaxed={:.6} fd={:.6}",
        frozen[0][0], relaxed[0][0], fd[0][0]
    );
    println!(
        "hcp spring C12 frozen={:.6} relaxed={:.6} fd={:.6}",
        frozen[0][1], relaxed[0][1], fd[0][1]
    );
    println!(
        "hcp spring C44 frozen={:.6} relaxed={:.6} fd={:.6}",
        frozen[3][3], relaxed[3][3], fd[3][3]
    );
    assert!(
        (relaxed[0][0] - frozen[0][0]).abs() > 1e-3,
        "internal relaxation should change C11"
    );
    for (i, j) in [(0, 0), (0, 1), (3, 3)] {
        rel_close(relaxed[i][j], fd[i][j], 1e-4, "hcp spring relaxed vs FD");
    }

    // Lennard-Jones hcp at its isotropic equilibrium: relaxation changes the
    // constants substantially; C11 and C12 match the finite difference. (C44
    // is sensitive to the small deviatoric stress of the ideal c/a, so the
    // zero-stress spring above is the sharp test for it.)
    let lj = LennardJones::new(1.0, 1.0);
    let raw = relax_cutoff(&Lattice::hcp(1.0), &lj, 8.0, 0.9, 1.8);
    let ljat = Lattice::hcp(raw.scale);
    let lj_frozen = elastic_tensor_frozen(&ljat, &lj, 8.0);
    let lj_relaxed = elastic_tensor(&ljat, &lj, 8.0);
    let lj_fd = elastic_tensor_fd(&ljat, &lj, 8.0, 1e-4);
    println!(
        "hcp LJ C11 frozen={:.5} relaxed={:.5} fd={:.5}",
        lj_frozen[0][0], lj_relaxed[0][0], lj_fd[0][0]
    );
    assert!(
        (lj_relaxed[0][0] - lj_frozen[0][0]).abs() > 5.0,
        "internal relaxation should matter for LJ hcp"
    );
    for (i, j) in [(0, 0), (0, 1)] {
        rel_close(lj_relaxed[i][j], lj_fd[i][j], 1e-4, "hcp LJ relaxed vs FD");
    }
}

#[test]
fn monatomic_chain_dispersion_is_exact() {
    let a = 1.0;
    let spring = Harmonic { k: 2.0, r0: a };
    let harms = Harmonics::unit_masses(&chain_lattice(a), 1.5 * a);
    for frac in [0.02, 0.05, 0.1, 0.25, 0.5] {
        let k = frac * std::f64::consts::PI / a;
        let w = harms.frequencies(&spring, [k, 0.0, 0.0]);
        // Only the longitudinal mode has non-zero stiffness; the two
        // transverse modes are exactly zero.
        assert!(w[0].abs() < 1e-12 && w[1].abs() < 1e-12);
        close(
            w[2],
            chain_frequency(k, a, spring.k, 1.0),
            1e-12,
            "chain omega",
        );
    }
}

#[test]
fn acoustic_sum_rule_and_dynamical_stability_of_fcc_lj() {
    let lj = LennardJones::new(1.0, 1.0);
    let raw = relax_cutoff(&Lattice::fcc(1.0), &lj, 8.0, 0.9, 1.8);
    let harms = Harmonics::unit_masses(&Lattice::fcc(raw.scale), 8.0);

    // Acoustic sum rule: all three Gamma eigenvalues vanish.
    let gamma = harms.eigenvalues(&lj, [0.0, 0.0, 0.0]);
    assert_eq!(gamma.len(), 3);
    for e in &gamma {
        assert!(e.abs() < 1e-10, "Gamma eigenvalue {e:e} not zero");
    }

    // Stability: at the equilibrium lattice constant every eigenvalue is
    // non-negative, and the acoustic branches rise from zero.
    let a = raw.scale;
    let mut seen_positive = false;
    for frac in [0.05, 0.1, 0.25, 0.5, 0.75] {
        for dir in [[1.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.5, 0.5]] {
            let k = [
                dir[0] * frac * std::f64::consts::PI / a,
                dir[1] * frac * std::f64::consts::PI / a,
                dir[2] * frac * std::f64::consts::PI / a,
            ];
            for e in harms.eigenvalues(&lj, k) {
                assert!(e > -1e-9, "unstable eigenvalue {e:e}");
                if e > 1e-6 {
                    seen_positive = true;
                }
            }
        }
    }
    assert!(seen_positive);
}

#[test]
fn fcc_phonon_path_is_sampled_monotonically() {
    let lj = LennardJones::new(1.0, 1.0);
    let raw = relax_cutoff(&Lattice::fcc(1.0), &lj, 8.0, 0.9, 1.8);
    let harms = Harmonics::unit_masses(&Lattice::fcc(raw.scale), 8.0);
    let points = sample_path(&harms, &lj, &fcc_path(), two_pi_over(raw.scale), 6);
    assert_eq!(points.len(), 4 * 6 + 1);
    for w in points.windows(2) {
        assert!(w[1].0 > w[0].0, "path distance must increase");
    }
    for (_, freqs) in &points {
        assert_eq!(freqs.len(), 3);
        assert!(freqs.iter().all(|w| w.is_finite() && *w >= 0.0));
    }
    // The path ends back at Gamma, so the final frequencies vanish.
    let last = points.last().unwrap();
    assert!(last.1.iter().all(|w| *w < 1e-10));
}

#[test]
fn dispersion_indicator_is_quadratic_at_long_wavelength() {
    // The relative dispersion error (ω_lat - v k)/ω_lat of the longitudinal
    // [100] branch of fcc LJ must scale as (k a)^2 as k -> 0.
    let lj = LennardJones::new(1.0, 1.0);
    let raw = relax_cutoff(&Lattice::fcc(1.0), &lj, 8.0, 0.9, 1.8);
    let lat = Lattice::fcc(raw.scale);
    let c = elastic_tensor(&lat, &lj, 8.0);
    let rho = lat.natoms() as f64 / lat.volume();
    let harms = Harmonics::unit_masses(&lat, 8.0);
    let n = [1.0, 0.0, 0.0];
    let v = christoffel_velocities(&c, rho, n);
    let long_error = |q: f64| {
        let w = harms.acoustic_velocities(&lj, n, q);
        dispersion_error(w[2] * q, v[2], q).abs()
    };
    let e_half = long_error(5e-3);
    let e_full = long_error(1e-2);
    let ratio = e_full / e_half;
    println!("dispersion error q=5e-3: {e_half:.3e}, q=1e-2: {e_full:.3e}, ratio={ratio:.3}");
    assert!((3.5..4.5).contains(&ratio), "not quadratic: ratio {ratio}");
}

#[test]
fn christoffel_sound_speeds_are_the_long_wavelength_phonon_limit() {
    let lj = LennardJones::new(1.0, 1.0);
    let raw = relax_cutoff(&Lattice::fcc(1.0), &lj, 8.0, 0.9, 1.8);
    let lat = Lattice::fcc(raw.scale);
    let c = elastic_tensor(&lat, &lj, 8.0);
    let rho = lat.natoms() as f64 / lat.volume();
    let harms = Harmonics::unit_masses(&lat, 8.0);

    for n in [
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [1.0, 1.0, 1.0],
        [0.5, 0.5, 1.0],
    ] {
        let v_cont = christoffel_velocities(&c, rho, n);
        let q_small = 1e-3;
        let v_lat = harms.acoustic_velocities(&lj, n, q_small);
        assert_eq!(v_lat.len(), 3);
        for mode in 0..3 {
            let err = dispersion_error(v_lat[mode] * q_small, v_cont[mode], q_small).abs();
            println!(
                "n={n:?} mode={mode} v_lat={:.8} v_cont={:.8} rel={err:.3e}",
                v_lat[mode], v_cont[mode]
            );
            assert!(
                err < 1e-5,
                "Christoffel/phonon mismatch n={n:?} mode={mode}: {err:e}"
            );
        }

        // The indicator must vanish as k -> 0 like (k a)^2.
        let mut prev = f64::INFINITY;
        for q in [1e-2, 1e-3] {
            let v = harms.acoustic_velocities(&lj, n, q);
            let e: f64 = (0..3)
                .map(|m| dispersion_error(v[m] * q, v_cont[m], q).abs())
                .fold(0.0, f64::max);
            assert!(e < prev, "dispersion error did not shrink with q");
            prev = e;
        }
    }
    // Check the helper on a known case: unit direction and length.
    close(k_norm([3.0, 4.0, 0.0]), 5.0, 1e-15, "k_norm");
}
