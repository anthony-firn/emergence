//! Discrete energy budget for the Cahn–Hilliard + variable-density NS solver.
//!
//! Inviscid, zero mobility: the capillary force must be the exact variational
//! derivative of the discrete free energy so that kinetic energy and free
//! energy exchange without drift (up to the time-stepping dissipation).
use num_complex::Complex64;
use spectral_core::fft::Fft3;
use spectral_core::twophase::{Projection, TwoPhaseParams, TwoPhaseSolver};
use spectral_core::{init, Grid};

const PI: f64 = std::f64::consts::PI;
const I: Complex64 = Complex64::new(0.0, 1.0);

fn build(n: usize, ratio: f64, etas: f64, mobility: f64, mu1: f64) -> TwoPhaseSolver {
    let g = Grid::new(n);
    let eta = etas * g.dx();
    let prm = TwoPhaseParams {
        rho: [1.0, 1.0 / ratio],
        mu: [mu1, mu1],
        sigma: 1.0,
        gravity: [0.0; 3],
        eta,
        mobility,
        stab: 2.0,
    };
    let a0 = 0.03;
    let k = 2.0;
    let w = 2f64.sqrt() * eta;
    let (zlo, zhi) = (0.5 * PI, 1.5 * PI);
    let phi = init::from_fn(&g, |x, _, z| {
        let zeta = a0 * (k * x).cos();
        [((z - zlo - zeta) / w).tanh() - ((z - zhi - zeta) / w).tanh() - 1.0, 0.0, 0.0]
    })[0]
        .clone();
    let u0 = init::from_fn(&g, |_, _, _| [0.0; 3]);
    TwoPhaseSolver::new(g.clone(), prm, u0, phi)
}

fn free_energy(s: &TwoPhaseSolver) -> f64 {
    let g = &s.grid;
    let fft = Fft3::new(g);
    let phi = s.phi_phys();
    // spectral gradient magnitude squared: sum_d |D_d phi|^2
    let nr = g.real_len();
    let mut grad2 = vec![0.0f64; nr];
    for d in 0..3 {
        let mut ds: Vec<Complex64> = s
            .phi
            .iter()
            .enumerate()
            .map(|(i, v)| I * g.k(i)[d] * *v)
            .collect();
        let mut w = vec![0.0; nr];
        fft.inverse_inplace(&mut ds, &mut w);
        for q in 0..nr {
            grad2[q] += w[q] * w[q];
        }
    }
    let eta2 = s.prm.eta * s.prm.eta;
    let dens: f64 = (0..nr)
        .map(|q| 0.5 * grad2[q] + (phi[q] * phi[q] - 1.0).powi(2) / (4.0 * eta2))
        .sum::<f64>()
        / nr as f64;
    let lambda = 3.0 * s.prm.sigma * s.prm.eta / (2.0 * 2f64.sqrt());
    lambda * dens
}

fn kinetic(s: &TwoPhaseSolver) -> f64 {
    s.kinetic_energy()
}

#[test]
fn inviscid_energy_conserved() {
    for &n in &[24usize, 32] {
        let mut s = build(n, 1.0, 1.2, 0.0, 0.0);
        s.set_projection(Projection::Constant);
        let dt = 0.002;
        let e0 = free_energy(&s) + kinetic(&s);
        let mut emin = f64::INFINITY;
        let mut emax = f64::NEG_INFINITY;
        let fe0 = free_energy(&s);
        let ke0 = kinetic(&s);
        for step in 0..300 {
            s.step(dt);
            let e = free_energy(&s) + kinetic(&s);
            emin = emin.min(e);
            emax = emax.max(e);
            if step % 100 == 99 {
                eprintln!(
                    "n={n} step={step} E={e:.8e} dE={:+.3e} F={:.6e} K={:.6e}",
                    e - e0,
                    free_energy(&s),
                    kinetic(&s)
                );
            }
        }
        let e1 = free_energy(&s) + kinetic(&s);
        let rel = (e1 - e0) / e0.abs();
        eprintln!(
            "n={n}: E0={e0:.8e} E1={e1:.8e} rel drift={rel:+.3e} min={:.8e} max={:.8e} F: {fe0:.4e}->{:.4e} K: {ke0:.4e}->{:.4e}",
            emin,
            emax,
            free_energy(&s),
            kinetic(&s)
        );
        // The capillary force must be the variational derivative of the
        // discrete free energy: the inviscid, M = 0 budget can only decrease
        // (time-stepping dissipation), and it must not drift upward.
        assert!(e1 <= e0 + 1e-12 * e0.abs(), "total energy increased by {rel:+.3e}");
        assert!(rel.abs() < 1e-4, "discrete energy drift {rel:+.3e}");
    }
}
