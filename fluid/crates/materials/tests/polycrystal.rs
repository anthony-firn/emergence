//! Benchmarks for the polycrystal fold of the materials ladder
//! (`analysis/materials.md`): the isotropic averages against measured
//! polycrystalline values, the Hashin-Shtrikman bounds, the parameters the
//! isotropic GPR solid model consumes, and the grain-size validity indicator.

use materials::anisotropic::{cubic, material};
use materials::polycrystal::{
    hashin_shtrikman_cubic, hill, isotropic_parameters_stiffness, isotropic_polycrystal_valid,
    rayleigh_longitudinal_attenuation, reuss, scattering_regime, voigt, ScatteringRegime,
};

fn close_rel(actual: f64, expected: f64, tol: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tol * expected.abs().max(1e-300),
        "{what}: {actual:.6e} vs {expected:.6e}, rel tolerance {tol:.2e}"
    );
}

/// Measured polycrystalline moduli (Ashby & Jones, *Engineering Materials*;
/// CRC Handbook; typical room-temperature quasi-static values): copper
/// shear ~48 GPa, Young's ~128 GPa; aluminium shear ~26 GPa, Young's
/// ~70 GPa. The Hill average of the single-crystal constants should land
/// within a few percent.
#[test]
fn hill_average_matches_measured_polycrystal_moduli() {
    for (name, g_meas, e_meas) in [("Cu", 48.0, 128.0), ("Al", 26.0, 70.0)] {
        let m = material(name).expect("material in table");
        let h = hill(&m.c);
        println!(
            "{name}: Hill G={:.3} GPa (meas {g_meas}), E={:.3} GPa (meas {e_meas})",
            h.shear / 1e9,
            h.young() / 1e9
        );
        close_rel(h.shear / 1e9, g_meas, 5e-2, "polycrystal shear modulus");
        close_rel(h.young() / 1e9, e_meas, 5e-2, "polycrystal Young's modulus");
    }
}

/// Voigt and Reuss bracket the Hill average, and for a cubic crystal the
/// bulk modulus is orientation-independent (Voigt = Reuss).
#[test]
fn voigt_reuss_bracket_hill() {
    let m = material("Cu").expect("Cu in table");
    let v = voigt(&m.c);
    let r = reuss(&m.c);
    let h = hill(&m.c);
    assert!(v.shear >= h.shear && h.shear >= r.shear);
    close_rel(v.bulk, r.bulk, 1e-9, "cubic bulk Voigt = Reuss");
    close_rel(h.bulk, 0.5 * (v.bulk + r.bulk), 1e-12, "Hill bulk");
}

/// Hashin-Shtrikman shear bounds bracket the Hill average and lie inside the
/// Voigt-Reuss interval; the bounds collapse for a nearly isotropic crystal
/// (tungsten).
#[test]
fn hashin_shtrikman_bounds_are_consistent() {
    for name in ["Cu", "Al", "Ni", "W"] {
        let m = material(name).expect("material in table");
        let (lo, hi) = hashin_shtrikman_cubic(m.c[0][0], m.c[0][1], m.c[3][3]);
        let v = voigt(&m.c).shear;
        let r = reuss(&m.c).shear;
        let h = hill(&m.c).shear;
        println!(
            "{name}: HS=[{:.3},{:.3}] GPa, Hill={:.3}",
            lo / 1e9,
            hi / 1e9,
            h / 1e9
        );
        assert!(lo <= hi, "{name}: HS lower above upper");
        assert!(lo + 1e-6 >= r, "{name}: HS lower inside Reuss");
        assert!(hi - 1e-6 <= v, "{name}: HS upper inside Voigt");
        // Hill is not guaranteed to lie inside the (tighter) HS interval,
        // so allow a small relative slack.
        assert!(
            h >= lo - 1e-3 * lo && h <= hi + 1e-3 * hi,
            "{name}: Hill outside the HS bracket by more than 0.1%"
        );
    }
    // Tungsten's anisotropy is tiny, so the bounds are essentially collapsed.
    let w = material("W").expect("W in table");
    let (lo, hi) = hashin_shtrikman_cubic(w.c[0][0], w.c[0][1], w.c[3][3]);
    let h = hill(&w.c).shear;
    assert!(
        (hi - lo) / h < 1e-3,
        "W HS width {:.3e} is not small",
        hi - lo
    );
}

/// The GPR solid model consumes density, shear sound speed and bulk sound
/// speed. Check those are the isotropic moduli's `sqrt(G/ρ)`, `sqrt(K/ρ)`
/// and that the implied longitudinal speed is `sqrt((K+4G/3)/ρ)`.
#[test]
fn gpr_isotropic_parameters_are_consistent() {
    let cu = material("Cu").expect("Cu in table");
    let iso = isotropic_parameters_stiffness(&cu.c, cu.density);
    close_rel(iso.density, 8960.0, 1e-12, "density");
    close_rel(
        iso.shear_sound_speed,
        (iso.moduli.shear / iso.density).sqrt(),
        1e-12,
        "shear sound speed",
    );
    close_rel(
        iso.bulk_sound_speed,
        (iso.moduli.bulk / iso.density).sqrt(),
        1e-12,
        "bulk sound speed",
    );
    close_rel(
        iso.longitudinal_speed.powi(2),
        iso.bulk_sound_speed.powi(2) + 4.0 / 3.0 * iso.shear_sound_speed.powi(2),
        1e-12,
        "longitudinal speed identity",
    );
    println!(
        "Cu polycrystal: rho={:.0} kg/m3, cs={:.1} m/s, cb={:.1} m/s, vp={:.1} m/s",
        iso.density, iso.shear_sound_speed, iso.bulk_sound_speed, iso.longitudinal_speed
    );
}

/// Grain size vs wavelength selects the regime; the isotropic polycrystal is
/// valid only in the Rayleigh regime, `d ≪ λ`.
#[test]
fn grain_size_validity_indicator() {
    // 1 MHz shear wave in Cu (v ~ 2300 m/s): lambda ~ 2.3 mm.
    let lambda = 2.3e-3;
    assert_eq!(scattering_regime(1e-4, lambda), ScatteringRegime::Rayleigh);
    assert_eq!(
        scattering_regime(1e-3, lambda),
        ScatteringRegime::Stochastic
    );
    assert_eq!(scattering_regime(1e-1, lambda), ScatteringRegime::Geometric);
    assert!(isotropic_polycrystal_valid(1e-4, lambda));
    assert!(!isotropic_polycrystal_valid(1e-1, lambda));
}

/// Rayleigh attenuation scales as `d³ f⁴` and grows with the anisotropy:
/// copper (A ~ 3.2) attenuates more than aluminium (A ~ 1.2) at the same
/// grain size and frequency.
#[test]
fn rayleigh_attenuation_scaling_and_anisotropy() {
    let cu = material("Cu").expect("Cu in table");
    let al = material("Al").expect("Al in table");
    let d = 1e-4;
    let f = 5e6;
    let a_cu =
        rayleigh_longitudinal_attenuation(cu.c[0][0], cu.c[0][1], cu.c[3][3], cu.density, d, f);
    let a_al =
        rayleigh_longitudinal_attenuation(al.c[0][0], al.c[0][1], al.c[3][3], al.density, d, f);
    println!("Rayleigh alpha_L at d=0.1mm, f=5MHz: Cu={a_cu:.3e}, Al={a_al:.3e} 1/m");
    assert!(a_cu > a_al, "copper should attenuate more than aluminium");
    // A genuinely isotropic crystal (C11 = C12 + 2 C44) has no cubic
    // anisotropy and so no Rayleigh scattering within this leading-order
    // model.
    let iso = cubic(200e9, 100e9, 50e9);
    let a_iso = rayleigh_longitudinal_attenuation(iso[0][0], iso[0][1], iso[3][3], 8960.0, d, f);
    close_rel(
        a_iso,
        0.0,
        1e-12,
        "isotropic crystal has no Rayleigh scattering",
    );
}
