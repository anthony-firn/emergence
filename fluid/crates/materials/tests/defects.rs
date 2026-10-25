//! Benchmarks for the defect rung of the materials ladder.
//!
//! References are in `analysis/defects.md`. Every test checks an exact
//! identity (the unrelaxed vacancy energy, the broken-bond lattice sum), a
//! published value, or an ordering; no tolerance is widened to make a
//! model agree with itself.

use materials::defects::{
    broken_bond_surface_energy, fcc_100_slab, fcc_111_slab, gamma_surface_111,
    stacking_fault_hcp_fcc, surface_energy, unstable_stacking_fault, vacancy_formation,
};
use materials::elastic::elastic_tensor_extrapolated;
use materials::lattice::Lattice;
use materials::plasticity::*;
use materials::potential::{energy_per_atom, energy_per_atom_with_tail, relax, LennardJones};
use std::f64::consts::PI;

/// fcc LJ at zero pressure: `a = 1.5417373 σ`, `r_nn = 1.0901729 σ`,
/// `e_coh = -8.6102077 ε` (analysis/materials.md).
fn lj_fcc() -> (LennardJones, f64, f64) {
    let pot = LennardJones::new(1.0, 1.0);
    let a = relax(&Lattice::fcc(1.0), &pot, 8.0, 1.4, 1.7).scale;
    (pot, a, a / 2.0_f64.sqrt())
}

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
fn unrelaxed_vacancy_energy_is_minus_cohesive_energy() {
    // Exact pair-potential identity: removing one atom from the perfect
    // crystal removes a per-atom bond sum 2 e_coh; referring the remaining
    // N-1 atoms to the bulk gives E_v^unrel = -e_coh.
    let (pot, a, _) = lj_fcc();
    let rc = 6.0;
    let e_coh = energy_per_atom_with_tail(&Lattice::fcc(a), &pot, rc);
    // The rc = 6 tail correction is converged to ~2e-3 ε; the exact
    // converged LJ value is -8.6102077 ε (analysis/materials.md).
    close(e_coh, -8.6102077, 3e-3, "LJ fcc cohesive energy");
    let v = vacancy_formation(&Lattice::fcc(a), &pot, rc, [4, 4, 4]);
    // The unrelaxed formation energy must equal -e_coh evaluated at the
    // *same* cutoff and tail convention; the finite supercell shifts it by
    // a few 1e-3 ε.
    close(v.e_unrelaxed, -e_coh, 2e-2, "unrelaxed E_v = -e_coh");
    assert!(
        v.e_relaxed < v.e_unrelaxed,
        "relaxation must lower the vacancy energy"
    );
    assert!(
        v.relaxation_energy > 0.08 && v.relaxation_energy < 0.16,
        "LJ vacancy relaxation energy {} ε outside expected range",
        v.relaxation_energy
    );
    close(v.e_relaxed, 8.494, 1e-2, "relaxed LJ fcc vacancy energy");
}

#[test]
fn surface_energy_equals_broken_bond_lattice_sum() {
    // The rigid (unrelaxed) slab surface energy must equal the exact
    // broken-bond lattice sum gamma = 1/2 rho sum_{r.n>0} (r.n)|phi|.
    let (pot, a, r_nn) = lj_fcc();
    let base = Lattice::fcc(a);
    for rc in [4.0, 5.0, 6.0] {
        let e_bulk = energy_per_atom(&base, &pot, rc);
        let s111 = fcc_111_slab(r_nn, 10, rc + 3.0);
        let (u111, r111) = surface_energy(&s111, &pot, rc, e_bulk, 1e-4);
        let bb111 = broken_bond_surface_energy(&base, &pot, rc, [1.0, 1.0, 1.0]);
        rel_close(u111, bb111, 1e-6, "111 broken-bond vs slab");

        let s100 = fcc_100_slab(a, 12, rc + 3.0);
        let (u100, r100) = surface_energy(&s100, &pot, rc, e_bulk, 1e-4);
        let bb100 = broken_bond_surface_energy(&base, &pot, rc, [0.0, 0.0, 1.0]);
        rel_close(u100, bb100, 1e-6, "100 broken-bond vs slab");

        // Relaxation lowers the surface energy modestly (~0.4 % for LJ).
        assert!(r111 < u111 && r100 < u100, "relaxation must lower gamma");
        assert!(
            (1.0 - r111 / u111).abs() < 0.02,
            "LJ surface relaxation {:.3} too large",
            1.0 - r111 / u111
        );
        // Close-packed (111) is the low-energy face.
        assert!(r111 < r100, "gamma(111) must be below gamma(100)");
    }
}

#[test]
fn surface_energies_are_in_the_literature_range() {
    // Converged broken-bond lattice sums (rc = 12 σ): the classical
    // zero-temperature LJ values in reduced units. These are exact for the
    // rigid lattice of a central pair potential.
    let (pot, a, r_nn) = lj_fcc();
    let base = Lattice::fcc(a);
    let bb111 = broken_bond_surface_energy(&base, &pot, 12.0, [1.0, 1.0, 1.0]);
    let bb100 = broken_bond_surface_energy(&base, &pot, 12.0, [0.0, 0.0, 1.0]);
    close(bb111, 3.100, 5e-3, "gamma(111) broken-bond rc=12");
    close(bb100, 3.220, 5e-3, "gamma(100) broken-bond rc=12");

    // Relaxed slab values at rc = 6 (the same-cutoff reference).
    let e_bulk = energy_per_atom(&base, &pot, 6.0);
    let (_, r111) = surface_energy(&fcc_111_slab(r_nn, 16, 9.0), &pot, 6.0, e_bulk, 1e-4);
    let (_, r100) = surface_energy(&fcc_100_slab(a, 20, 9.0), &pot, 6.0, e_bulk, 1e-4);
    close(r111, 3.010, 1e-2, "relaxed gamma(111) rc=6");
    close(r100, 3.111, 1e-2, "relaxed gamma(100) rc=6");
    assert!(r111 > 2.5 && r111 < 3.3 && r100 > 2.6 && r100 < 3.4);
}

#[test]
fn intrinsic_stacking_fault_is_near_zero_and_negative() {
    // For a pair potential gamma_ISF = 2 (e_hcp - e_fcc) / A_111. LJ puts
    // hcp just below fcc, so the predicted ISF energy is a small negative
    // number -- the well-known pair-potential failure.
    let (pot, a, r_nn) = lj_fcc();
    let sf = stacking_fault_hcp_fcc(a, &pot, 12.0);
    assert!(
        sf.e_hcp < sf.e_fcc,
        "LJ hcp must lie below fcc (analysis/materials.md)"
    );
    close(sf.area_per_atom, 1.0293, 1e-3, "A_111 per atom");
    close(sf.gamma, -1.68e-3, 5e-4, "gamma_ISF hcp-fcc epoch");
    assert!(sf.gamma.abs() < 5e-3, "LJ ISF energy must be tiny");

    // Direct constrained (normal-only) fault at t = 1 is also near zero,
    // confirming the two independent routes agree in magnitude.
    let u = [r_nn, 0.0, 0.0];
    let v = [0.5 * r_nn, 0.5 * 3.0_f64.sqrt() * r_nn, 0.0];
    let partial = [(u[0] + v[0]) / 3.0, (u[1] + v[1]) / 3.0, 0.0];
    let curve = gamma_surface_111(r_nn, 9, 9.0, 4, 6.0, &pot, partial, 12, 1e-4);
    let isf = curve.last().unwrap().gamma_relaxed;
    assert!(isf.abs() < 5e-3, "direct ISF gamma {isf} not near zero");
}

#[test]
fn gamma_surface_and_unstable_stacking_fault() {
    let (pot, _, r_nn) = lj_fcc();
    let u = [r_nn, 0.0, 0.0];
    let v = [0.5 * r_nn, 0.5 * 3.0_f64.sqrt() * r_nn, 0.0];
    let partial = [(u[0] + v[0]) / 3.0, (u[1] + v[1]) / 3.0, 0.0];

    // <112> path to the Shockley partial.
    let c = gamma_surface_111(r_nn, 9, 9.0, 4, 6.0, &pot, partial, 12, 1e-4);
    let u112 = unstable_stacking_fault(&c, false);
    let r112 = unstable_stacking_fault(&c, true);
    close(u112, 0.5514, 2e-3, "unrelaxed gamma_us<112>");
    close(r112, 0.4226, 3e-3, "relaxed gamma_us<112>");
    assert!(r112 > 0.0 && r112 < u112, "relaxation must lower gamma_us");
    // The path starts at the perfect crystal and ends at the ISF.
    close(c[0].gamma_relaxed, 0.0, 1e-12, "gamma(0)");
    assert!(c.last().unwrap().gamma_relaxed.abs() < 5e-3);

    // <110> path (period r_nn); much larger barrier.
    let c2 = gamma_surface_111(r_nn, 9, 9.0, 4, 6.0, &pot, u, 12, 1e-4);
    let u110 = unstable_stacking_fault(&c2, false);
    let r110 = unstable_stacking_fault(&c2, true);
    close(u110, 3.226, 5e-3, "unrelaxed gamma_us<110>");
    close(r110, 1.251, 5e-3, "relaxed gamma_us<110>");
    assert!(r110 > r112, "the <110> barrier should exceed <112>");
}

#[test]
fn peierls_nabarro_and_the_fold_to_elasticity() {
    // LJ fcc elastic constants (infinite-cutoff extrapolation,
    // analysis/materials.md): C11 = 105.182, C12 = C44 = 60.187.
    let (pot, a, r_nn) = lj_fcc();
    let c = elastic_tensor_extrapolated(&Lattice::fcc(a), &pot, 12.0, 20.0);
    let (g, nu) = isotropic_from_cubic(c[0][0], c[0][1], c[3][3]);
    close(g, 40.5746, 1e-3, "Hill shear modulus");
    close(nu, 0.27131, 1e-4, "Poisson ratio");
    close(c[3][3], c[0][1], 1e-3, "central-force Cauchy relation C44=C12");

    let b = fcc_burgers(a);
    let d = fcc_slip_plane_spacing(a);
    close(b, r_nn, 1e-12, "Burgers vector = nn distance");
    let pn = peierls_nabarro(g, nu, b, d);
    close(pn.zeta, d / (2.0 * (1.0 - nu)), 1e-12, "core width");
    close(pn.zeta, 0.61077, 1e-4, "LJ core width");
    close(
        pn.sigma_p,
        (2.0 * g / (1.0 - nu)) * (-2.0 * PI * pn.zeta / b).exp(),
        1e-12,
        "Peierls stress definition",
    );
    close(pn.sigma_p / g, 0.0812, 1e-3, "sigma_P / G");
    // Elastic/sinusoidal amplitude that reproduces the core width.
    close(pn.gamma_us_sinusoidal, 2.7445, 1e-3, "sinusoidal gamma_us");
    close(
        sinusoidal_core_width(g, nu, b, pn.gamma_us_sinusoidal),
        pn.zeta,
        1e-12,
        "sinusoidal width consistency",
    );
    // Exact continuum soliton profile.
    close(sinusoidal_core(0.0, b, pn.zeta), b / 2.0, 1e-12, "soliton centre");
    assert!(sinusoidal_core(-50.0 * pn.zeta, b, pn.zeta) < 0.01 * b);
    close(sinusoidal_core(50.0 * pn.zeta, b, pn.zeta), b, 0.01 * b, "soliton tail");
    // The measured <112> gamma_us is far below the elastic amplitude, so the
    // naive exponential PN closure is inconsistent -- documented as crude.
    let measured = 0.4226;
    assert!(
        measured < pn.gamma_us_sinusoidal / 4.0,
        "LJ gamma_us is far from the elastic sinusoidal estimate"
    );
    assert!(fault_plane_strength(measured, b) < pn.sigma_p);

    // Fold: Taylor factor and the GPR plastic closure.
    let sigma_y = taylor_yield(pn.sigma_p, TAYLOR_FCC);
    close(sigma_y, TAYLOR_FCC * pn.sigma_p, 1e-12, "Taylor yield");
    close(TAYLOR_FCC, 3.06, 1e-12, "fcc Taylor factor");
    let hard = taylor_hardening(0.3, g, b, 1e-3);
    close(hard, 0.3 * g * b * (1e-3f64).sqrt(), 1e-12, "Taylor hardening");
    let tau0 = gpr_tau0_from_drag(1.0, g, b, 1e-3, TAYLOR_FCC);
    close(tau0, TAYLOR_FCC / (g * 1e-3 * b * b), 1e-12, "gpr tau0 mapping");
    let gp = gpr_plastic(sigma_y, tau0, 20.0);
    close(gp.sigma_y, sigma_y, 1e-12, "GPR sigma_y");
    close(gp.tau0, tau0, 1e-12, "GPR tau0");
    assert_eq!(gp.n, 20.0);
}

#[test]
fn plasticity_validity_indicator() {
    let sigma_p = 3.2958;
    close(dislocation_spacing(1e-4), 100.0, 1e-12, "1/sqrt(rho)");
    // Below the Peierls stress the response is elastic.
    assert_eq!(plasticity_rung(1.0, 1e-2, 1.0, sigma_p), PlasticityRung::Elastic);
    // One dislocation per 10 length units: spacing 10 > sample 1 -> discrete.
    assert_eq!(
        plasticity_rung(1.0, 1e-2, 10.0, sigma_p),
        PlasticityRung::Discrete
    );
    // Dense ensemble and large sample: the continuum yield law applies.
    assert_eq!(
        plasticity_rung(100.0, 1e-2, 10.0, sigma_p),
        PlasticityRung::ContinuumPlastic
    );
    // Smaller is stronger: G b / L grows as L shrinks.
    assert!(source_limited_strength(1.0, 1.0, 0.5) > source_limited_strength(1.0, 1.0, 2.0));
    // Hall-Petch: sqrt(d) scaling.
    let s1 = hall_petch(1.0, 1.0, 1.0);
    let s4 = hall_petch(1.0, 1.0, 4.0);
    close(s1, 2.0, 1e-12, "Hall-Petch d=1");
    close(s4, 1.5, 1e-12, "Hall-Petch d=4");
}
