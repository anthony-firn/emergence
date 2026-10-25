//! The atomistic rung and its folds, for the fully self-consistent material
//! (Lennard-Jones argon).
//!
//! Three folds are checked against exact or published references:
//! * the Cauchy–Born stiffness of the LJ fcc crystal against the converged
//!   infinite-lattice-sum values of `analysis/materials.md`;
//! * the Christoffel sound speeds against the long-wavelength harmonic
//!   phonon slopes (the Cauchy–Born fold);
//! * an NVT molecular-dynamics state point against the Johnson–Zollweg–
//!   Gubbins (1993) reference equation of state.

use crate::Measurement;
use materials::elastic::elastic_tensor_extrapolated;
use materials::lattice::Lattice;
use materials::phonon::Harmonics;
use materials::potential::{relax_cutoff, LennardJones};
use materials::V3;

use crate::Level;

/// Converged (infinite-sum) LJ fcc elastic constants, in `ε/σ³`
/// (`analysis/materials.md`: C11 = 105.182 ± 0.02, C12 = C44 = 60.187 ± 0.02).
pub const LJ_C11_REF: f64 = 105.182;
pub const LJ_C12_REF: f64 = 60.187;

/// The argon Cauchy–Born stiffness in SI, the mass density and the
/// zero-pressure lattice constant.
pub fn argon_cauchy_born() -> (materials::M6, f64, f64) {
    let lj = LennardJones::new(1.0, 1.0);
    // Zero-stress (raw-cutoff) equilibrium, consistent with the phonons.
    let raw = relax_cutoff(&Lattice::fcc(1.0), &lj, 12.0, 0.9, 1.8);
    let a_red = raw.scale;
    let lat = Lattice::fcc(a_red);
    let c_red = elastic_tensor_extrapolated(&lat, &lj, 12.0, 20.0);
    let p_unit = md::argon::pressure_unit();
    let c = c_red.map(|row| row.map(|x| x * p_unit));
    let number_density = 4.0 / a_red.powi(3); // atoms per σ³
    let density = number_density * md::argon::mass_kg() / md::argon::SIGMA_M.powi(3);
    (c, density, a_red * md::argon::SIGMA_M)
}

/// The LJ fcc sound-speed fold: Christoffel (from the Cauchy–Born tensor)
/// versus the harmonic phonon slope at `q = 1e-3`. Returns the largest
/// relative error over the three directions and modes.
pub fn argon_christoffel_vs_phonon() -> f64 {
    let lj = LennardJones::new(1.0, 1.0);
    let raw = relax_cutoff(&Lattice::fcc(1.0), &lj, 8.0, 0.9, 1.8);
    let lat = Lattice::fcc(raw.scale);
    let c = materials::elastic::elastic_tensor(&lat, &lj, 8.0);
    let rho = 4.0 / raw.scale.powi(3);
    let harms = Harmonics::unit_masses(&lat, 8.0);
    let dirs: [V3; 3] = [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0]];
    let mut worst: f64 = 0.0;
    for n in dirs {
        let vc = materials::continuum::christoffel_velocities(&c, rho, n);
        let vl = harms.acoustic_velocities(&lj, n, 1e-3);
        for mode in 0..3 {
            let e = (vl[mode] - vc[mode]).abs() / vc[mode];
            worst = worst.max(e);
        }
    }
    worst
}

/// A short NVT state point and its JZG reference. Returns
/// `(T, P, U_res, T_ref, P_ref, U_ref)`.
pub fn argon_md_state_point(eq: usize, prod: usize, seed: u64) -> (f64, f64, f64, f64, f64, f64) {
    use md::eos::JzgEos;
    use md::potential::Cutoff;
    use md::system::System;
    use md::thermostat::Thermostat;

    let lj = || md::potential::Potential::LennardJones(LennardJones::new(1.0, 1.0));
    let (t, rho) = (2.0_f64, 0.5_f64);
    let mut sys = System::fcc(4, rho, lj(), 2.5, 0.3, Cutoff::TailCorrected, seed);
    sys.set_velocities_mb(t);
    let mut th = Thermostat::bussi(0.5);
    let dt = 0.005;
    for _ in 0..eq {
        th.step(&mut sys, dt, t);
    }
    let (mut p, mut u, mut tt) = (0.0, 0.0, 0.0);
    for _ in 0..prod {
        th.step(&mut sys, dt, t);
        p += sys.pressure();
        u += sys.potential_energy() / sys.n() as f64;
        tt += sys.temperature_com();
    }
    let n = prod as f64;
    let eos = JzgEos::new();
    (
        tt / n,
        p / n,
        u / n,
        t,
        eos.pressure(rho, t),
        eos.residual_energy(rho, t),
    )
}

/// All argon atomistic-fold measurements.
pub fn argon_fold(eq: usize, prod: usize, seed: u64) -> Vec<Measurement> {
    let mut out = Vec::new();
    let lj = LennardJones::new(1.0, 1.0);
    // Converged Cauchy–Born constants.
    let raw = relax_cutoff(&Lattice::fcc(1.0), &lj, 12.0, 0.9, 1.8);
    let lat = Lattice::fcc(raw.scale);
    let c_red = elastic_tensor_extrapolated(&lat, &lj, 12.0, 20.0);
    out.push(Measurement::new(
        "Ar",
        "atomistic fold",
        Level::Atomistic,
        "C11",
        c_red[0][0],
        LJ_C11_REF,
        "ε/σ³",
        (c_red[0][0] / LJ_C11_REF - 1.0).abs() < 0.01,
        "Cauchy–Born from the LJ pair potential (rc = 12σ, 20σ, extrapolated) \
         vs the converged infinite-lattice-sum reference in \
         analysis/materials.md",
    ));
    out.push(Measurement::new(
        "Ar",
        "atomistic fold",
        Level::Atomistic,
        "C12",
        c_red[0][1],
        LJ_C12_REF,
        "ε/σ³",
        (c_red[0][1] / LJ_C12_REF - 1.0).abs() < 0.01,
        "Cauchy relation C12 = C44 holds for central forces",
    ));
    out.push(Measurement::new(
        "Ar",
        "atomistic fold",
        Level::Atomistic,
        "C44",
        c_red[3][3],
        LJ_C12_REF,
        "ε/σ³",
        (c_red[3][3] / LJ_C12_REF - 1.0).abs() < 0.01,
        "Cauchy relation C12 = C44 holds for central forces",
    ));

    let err = argon_christoffel_vs_phonon();
    out.push(Measurement::new(
        "Ar",
        "atomistic fold",
        Level::Continuum,
        "max rel. phonon/Christoffel error",
        err,
        0.0,
        "",
        err < 1e-5,
        "Cauchy–Born C_ij → Christoffel speeds are the q→0 harmonic phonon \
         slopes (same raw-cutoff configuration); error is the (ka)² \
         dispersion residual at q = 1e-3",
    ));

    let (t, p, u, t_ref, p_ref, u_ref) = argon_md_state_point(eq, prod, seed);
    out.push(Measurement::new(
        "Ar",
        "atomistic fold",
        Level::Atomistic,
        "T*",
        t,
        t_ref,
        "",
        (t - t_ref).abs() < 0.05,
        "NVT Bussi state point (256 atoms, tail-corrected LJ)",
    ));
    out.push(Measurement::new(
        "Ar",
        "atomistic fold",
        Level::Atomistic,
        "P*",
        p,
        p_ref,
        "",
        (p - p_ref).abs() / p_ref.abs() < 0.05,
        "MD pressure vs JZG (1993) MBWR reference EOS, ρ* = 0.5, T* = 2.0",
    ));
    out.push(Measurement::new(
        "Ar",
        "atomistic fold",
        Level::Atomistic,
        "U_res*",
        u,
        u_ref,
        "",
        (u - u_ref).abs() < 0.05,
        "MD residual energy vs JZG (1993)",
    ));
    out
}
