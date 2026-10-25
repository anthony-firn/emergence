//! Hysteresis-fold benchmarks: the Stoner-Wohlfarth astroid, the
//! random-orientation coercivity, and the single-domain indicator.

use magnetism::hysteresis::{
    astroid_hsw, random_assembly_coercivity, single_domain_diameter, Macrospin,
};
use magnetism::micromag::{exchange_length, wall_delta, wall_sigma};
use magnetism::materials::{IRON, NDFEB, NICKEL, PERMALLOY};

fn rel(actual: f64, expected: f64) -> f64 {
    (actual - expected).abs() / expected.abs().max(1e-300)
}

#[test]
fn astroid_matches_llg_macrospin() {
    // h_sw(psi) = 1/(cos^{2/3} psi + sin^{2/3} psi)^{3/2}, in units of H_K.
    for deg in [0.0_f64, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0] {
        let psi = deg.to_radians();
        let exact = astroid_hsw(psi);
        let measured = Macrospin::switching_field_llg(psi, 1.0, 0.02, 40_000);
        assert!(
            rel(measured, exact) < 5e-2,
            "psi={deg}: LLG {measured:.5} vs astroid {exact:.5}"
        );
    }
    // The astroid touches the hard axis at 1 and has its minimum 1/2
    // at 45 degrees: h = 1/(2 (1/sqrt2)^{2/3})^{3/2} = 1/2.
    assert!(rel(astroid_hsw(0.0), 1.0) < 1e-12);
    assert!(rel(astroid_hsw(std::f64::consts::FRAC_PI_2), 1.0) < 1e-9);
    assert!(rel(astroid_hsw(std::f64::consts::FRAC_PI_4), 0.5) < 1e-9);
}

#[test]
fn random_orientation_coercivity_is_0479_hk() {
    // Stoner-Wohlfarth (1948): H_c = 0.479 H_K for an isotropic assembly.
    let hc = random_assembly_coercivity(300, 400, 1.0, 0.02, 300, 7).abs();
    assert!(
        rel(hc, 0.479) < 5e-2,
        "assembly coercivity {hc:.5} vs 0.479"
    );
}

#[test]
fn single_domain_indicator_for_real_materials() {
    // Iron: wall width ~ 21 nm, exchange length ~ 3.4 nm, single-domain
    // diameter ~ 10 nm (the numbers quoted in analysis/magnetism.md).
    let delta_fe = wall_delta(IRON.exchange_a, IRON.anisotropy_k);
    let lex_fe = exchange_length(IRON.exchange_a, IRON.ms);
    let dc_fe = single_domain_diameter(IRON.exchange_a, IRON.anisotropy_k, IRON.ms);
    assert!((delta_fe - 2.09e-8).abs() / 2.09e-8 < 5e-2, "Fe wall {delta_fe:.3e}");
    assert!((lex_fe - 3.38e-9).abs() / 3.38e-9 < 5e-2, "Fe l_ex {lex_fe:.3e}");
    assert!((dc_fe - 9.8e-9).abs() / 9.8e-9 < 1e-1, "Fe D_c {dc_fe:.3e}");

    // The indicator ordering for the material table.  Soft magnets have a
    // wall wider than the exchange length; the hard magnet NdFeB has a very
    // narrow wall (strong anisotropy), so only positivity is asserted there.
    for m in [IRON, NICKEL, PERMALLOY, NDFEB] {
        let k = m.anisotropy_k.max(1e3);
        let d = wall_delta(m.exchange_a, k);
        let lex = exchange_length(m.exchange_a, m.ms);
        let dc = single_domain_diameter(m.exchange_a, k, m.ms);
        assert!(d > 0.0 && lex > 0.0 && dc > 0.0, "{}", m.name);
        if !m.non_cubic {
            assert!(d > lex, "{} wall should exceed exchange length", m.name);
        }
        let _ = wall_sigma(m.exchange_a, k);
    }
}
