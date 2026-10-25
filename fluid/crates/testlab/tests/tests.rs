//! Standard-test validations for the virtual materials lab.
//!
//! Each test checks a measured value against a published or exact reference
//! and reports the error. The dynamic (plate-impact) tests are kept to small
//! 1D grids and a few hundred steps; each runs well under 30 s.

use testlab::dynamic::{gamma_for_slope, plate_impact, rh_us};
use testlab::elastic::{scattering_fold, static_elastic, ultrasonic};
use testlab::{aluminium, argon, copper, material, Measurement};

fn close(rel: f64, tol: f64, what: &str) {
    assert!(rel.abs() < tol, "{what}: relative error {rel:+.3e} > {tol:.3e}");
}

fn by_q<'a>(rows: &'a [Measurement], q: &str) -> &'a Measurement {
    rows.iter().find(|r| r.quantity == q).unwrap_or_else(|| panic!("no row {q}"))
}

#[test]
fn static_tensile_single_crystal_matches_the_closed_form() {
    for m in [copper(), aluminium(), argon()] {
        let rows = static_elastic(&m);
        for q in ["E [100]", "E [110]", "E [111]", "ν [100]", "ν [110]", "ν [111]"] {
            let r = by_q(&rows, q);
            close(r.rel_error(), 1e-9, &format!("{} {q}", m.name));
        }
    }
}

#[test]
fn polycrystal_hill_matches_handbook_for_cu_and_al() {
    for m in [copper(), aluminium()] {
        let rows = static_elastic(&m);
        let e = by_q(&rows, "E (Hill)");
        let nu = by_q(&rows, "ν (Hill)");
        close(e.rel_error(), 0.02, &format!("{} E Hill", m.name));
        close(nu.error(), 0.02, &format!("{} ν Hill", m.name));
    }
}

#[test]
fn christoffel_speeds_match_the_independent_solver_and_bracket_the_polycrystal() {
    for m in [copper(), aluminium(), argon()] {
        let rows = ultrasonic(&m);
        for dir in ["[100]", "[110]", "[111]"] {
            let l = by_q(&rows, &format!("v_L {dir}"));
            let t = by_q(&rows, &format!("v_T {dir}"));
            close(l.rel_error(), 1e-9, &format!("{} v_L {dir}", m.name));
            close(t.rel_error(), 1e-9, &format!("{} v_T {dir}", m.name));
        }
        // The polycrystal fold must sit inside the single-crystal bracket.
        assert!(by_q(&rows, "v_L (Hill)").fold_valid, "{} v_L bracket", m.name);
        assert!(by_q(&rows, "v_T (Hill)").fold_valid, "{} v_T bracket", m.name);
    }
}

#[test]
fn rayleigh_fold_is_valid_and_geometric_fold_is_not() {
    for m in [copper(), aluminium(), argon()] {
        let fine = scattering_fold(&m, 1.0e-5, 5.0e6); // d/λ ≪ 0.1
        let coarse = scattering_fold(&m, 1.0e-2, 5.0e6); // d/λ ≫ 0.1
        assert!(fine.fold_valid, "{} fine grain should be Rayleigh: {}", m.name, fine.note);
        assert!(!coarse.fold_valid, "{} coarse grain fold should be invalid", m.name);
    }
}

#[test]
fn argon_cauchy_born_reproduces_the_converged_lattice_sums() {
    let fold = testlab::atomistic::argon_fold(200, 400, 11);
    for q in ["C11", "C12", "C44"] {
        let r = by_q(&fold, q);
        close(r.rel_error(), 0.01, &format!("LJ {q}"));
    }
}

#[test]
fn argon_christoffel_speeds_are_the_phonon_slopes() {
    let err = testlab::atomistic::argon_christoffel_vs_phonon();
    println!("argon max phonon/Christoffel relative error {err:.3e}");
    assert!(err < 1e-5, "Cauchy–Born fold error {err:.3e} too large");
}

#[test]
fn argon_md_state_point_matches_the_jzg_eos() {
    // 256 atoms, a short NVT trajectory: the error is the statistical scatter
    // plus the tail correction, not a tuned tolerance.
    let fold = testlab::atomistic::argon_fold(1000, 2000, 11);
    let p = by_q(&fold, "P*");
    let u = by_q(&fold, "U_res*");
    println!("MD {} = {:.4} vs JZG {:.4}", p.quantity, p.value, p.reference);
    close(p.rel_error(), 0.08, "MD P*");
    close(u.rel_error(), 0.08, "MD U_res*");
}

/// Run a plate impact and return its measurements plus the result.
fn run_impact(name: &str) -> (testlab::dynamic::ImpactResult, Vec<Measurement>) {
    let m = material(name).unwrap();
    let res = plate_impact(&m, 200.0, 600, 0.4);
    let rows = res.measurements();
    (res, rows)
}

#[test]
fn plate_impact_cu() {
    let (res, rows) = run_impact("Cu");
    // Symmetric impact: particle velocity is exactly half the impact speed.
    let up = by_q(&rows, "particle velocity");
    close(up.rel_error(), 0.01, "Cu particle velocity");
    // Elastic precursor at the Hugoniot elastic limit.
    let hel = by_q(&rows, "σ_HEL");
    close(hel.rel_error(), 0.10, "Cu σ_HEL");
    // Plastic shock against the measured Hugoniot.
    let us = by_q(&rows, "shock speed");
    close(us.rel_error(), 0.02, "Cu shock speed");
    // The NASG/RH closure itself.
    let us_rh = by_q(&rows, "EOS shock speed");
    close(us_rh.rel_error(), 0.02, "Cu EOS shock speed");
    println!(
        "Cu: u_p {:.2}/{:.2} m/s, σ_HEL {:.1}/{:.1} MPa, Us {:.0}/{:.0} m/s, RH {:.0}",
        up.value, up.reference, hel.value, hel.reference, us.value, us.reference, res.shock_speed_rh
    );
}

#[test]
fn plate_impact_al() {
    let (res, rows) = run_impact("Al");
    let up = by_q(&rows, "particle velocity");
    close(up.rel_error(), 0.01, "Al particle velocity");
    let hel = by_q(&rows, "σ_HEL");
    close(hel.rel_error(), 0.10, "Al σ_HEL");
    let us = by_q(&rows, "shock speed");
    close(us.rel_error(), 0.02, "Al shock speed");
    let us_rh = by_q(&rows, "EOS shock speed");
    close(us_rh.rel_error(), 0.02, "Al EOS shock speed");
    println!(
        "Al: u_p {:.2}/{:.2} m/s, σ_HEL {:.1}/{:.1} MPa, Us {:.0}/{:.0} m/s, RH {:.0}",
        up.value, up.reference, hel.value, hel.reference, us.value, us.reference, res.shock_speed_rh
    );
}

#[test]
fn nasg_closure_reproduces_the_hugoniot_slope_only_at_weak_shocks() {
    // The stiffened-gas closure with γ = 4s − 1 matches c0 and the weak-shock
    // slope exactly by construction; the curvature is the EOS-form error.
    for m in [copper(), aluminium()] {
        let gamma = gamma_for_slope(m.hugoniot_s);
        let c0 = m.hugoniot_c0;
        let weak = rh_us(gamma, c0, 100.0);
        let strong = rh_us(gamma, c0, 500.0);
        let weak_ref = c0 + m.hugoniot_s * 100.0;
        let strong_ref = c0 + m.hugoniot_s * 500.0;
        println!(
            "{}: γ = {:.3}, Us(100) = {:.0} vs {:.0} ({:+.2}%), Us(500) = {:.0} vs {:.0} ({:+.2}%)",
            m.name,
            gamma,
            weak,
            weak_ref,
            100.0 * (weak / weak_ref - 1.0),
            strong,
            strong_ref,
            100.0 * (strong / strong_ref - 1.0)
        );
        close(weak / weak_ref - 1.0, 0.01, "weak-shock slope");
        close(strong / strong_ref - 1.0, 0.03, "strong-shock curvature");
    }
}

// ---------------------------------------------------------------------------
// Extended battery: mechanical, thermal, viscoelastic, magnetic, optical,
// interfaces and granular. Each test checks a published or exact reference.
// Heavy paths (quasi-harmonic argon, the small DEM cell) are cached inside
// the library, so repeated calls are cheap.
// ---------------------------------------------------------------------------

use testlab::{interfaces, magnetic, mechanical, optical, thermal, viscoelastic};

#[test]
fn griffith_fold_is_valid_for_brittle_argon_and_invalid_for_ductile_metals() {
    // Argon: atomistic surface energy + Cauchy-Born moduli vs the
    // independent anisotropic Stroh toughness.
    let ar = argon();
    let rows = mechanical::fracture(&ar);
    let kic = rows.iter().find(|r| r.quantity == "K_IC (atomistic γ_s)").unwrap();
    close(kic.rel_error(), 0.10, "Ar Griffith vs anisotropic");
    assert!(kic.fold_valid, "Ar Griffith fold should be certified: {}", kic.note);
    // 2D atomistic Griffith threshold: exact.
    let k2d = rows.iter().find(|r| r.quantity == "K_G (2D sheet)").unwrap();
    close(k2d.rel_error(), 1e-5, "2D atomistic Griffith");
    // Peierls stress reproduces the tabulated value.
    let sp = rows.iter().find(|r| r.quantity == "σ_P").unwrap();
    close(sp.rel_error(), 0.02, "Ar Peierls stress");

    // Metals: the Griffith value is two orders of magnitude below the
    // measured ductile toughness, and the fold is flagged invalid.
    for name in ["Cu", "Al"] {
        let m = material(name).unwrap();
        let rows = mechanical::fracture(&m);
        let kic = rows.iter().find(|r| r.quantity == "K_IC (Griffith)").unwrap();
        assert!(!kic.fold_valid, "{name} Griffith must be unfolded");
        assert!(kic.value < 0.05 * kic.reference, "{name} Griffith {:.3} vs {:}", kic.value, kic.reference);
        println!("{name}: Griffith K_IC {:.3} vs measured {:.1} MPa m^0.5", kic.value, kic.reference);
        // Yield: the handbook value is the bottom rung.
        let y = mechanical::metal_yield(&m);
        let sy = y.iter().find(|r| r.quantity == "σ_y").unwrap();
        assert!(sy.value.is_finite() && sy.value > 0.0);
    }
}

#[test]
fn argon_quasiharmonic_thermal_matches_experiment() {
    let rows = thermal::argon_thermal();
    // Thermal expansion at 4.2, 10 and 20 K within 2%.
    let mut checked = 0;
    for r in &rows {
        if r.test.contains("thermal expansion") {
            close(r.rel_error(), 0.02, &format!("Ar {}", r.quantity));
            checked += 1;
        }
    }
    assert!(checked >= 3, "expected at least three a(T) points");
    // Heat capacity: the well-resolved 40 K point is certified against the
    // Debye model with the calorimetric θ_D; the low-T points are an honest
    // miss (the fcc spectrum is not Debye-like) and are flagged invalid.
    let cv40 = rows.iter().find(|r| r.quantity == "C_V(40 K)/k_B").unwrap();
    close(cv40.rel_error(), 0.15, "Ar C_V(40 K) vs Debye");
    assert!(cv40.fold_valid, "40 K heat-capacity fold should be certified");
    let cv20 = rows.iter().find(|r| r.quantity == "C_V(20 K)/k_B").unwrap();
    assert!(!cv20.fold_valid, "20 K heat-capacity fold must be unfolded");
    println!("Ar C_V(40 K)/k_B = {:.3} vs Debye {:.3}", cv40.value, cv40.reference);
}

#[test]
fn silicon_conductivity_bulk_is_fourier_and_thin_film_is_not() {
    let rows = thermal::silicon_thermal();
    let bulk = rows.iter().find(|r| r.quantity == "κ(300 K, bulk)").unwrap();
    close(bulk.rel_error(), 0.20, "Si bulk kappa");
    assert!(bulk.fold_valid, "bulk must be Fourier");
    let film = rows.iter().find(|r| r.quantity == "κ(300 K, 100 nm film)").unwrap();
    assert!(!film.fold_valid, "100 nm film must be unfolded");
    assert!(film.value < bulk.value, "size effect must lower kappa");
    let kn = rows.iter().find(|r| r.quantity.starts_with("photon Knudsen")).unwrap();
    assert!(kn.value > 1.0, "Kn at 100 nm = {}", kn.value);
    println!("Si: bulk {:.1} W/mK (Fourier), 100 nm {:.1} (Kn={:.2})", bulk.value, film.value, kn.value);
}

#[test]
fn polystyrene_dma_matches_the_analytic_moduli_and_wlf() {
    let rows = viscoelastic::dma();
    let mut n_g = 0;
    for r in &rows {
        if r.test == "DMA sweep" && r.quantity.starts_with("G'") {
            close(r.rel_error(), 0.02, &format!("PS {}", r.quantity));
            n_g += 1;
        }
    }
    assert!(n_g >= 3, "expected three G' points");
    let gn0 = rows.iter().find(|r| r.quantity == "G_N^0").unwrap();
    close(gn0.rel_error(), 0.05, "G_N0");
    let wlf = rows.iter().find(|r| r.quantity == "log10 a_T (T_g+50 K)").unwrap();
    close(wlf.rel_error(), 0.01, "WLF");
    let tts = rows.iter().find(|r| r.quantity == "TTS collapse residual").unwrap();
    assert!(tts.value < 1e-12, "TTS residual {}", tts.value);
}

#[test]
fn magnetism_wall_widths_and_astroid() {
    let rows = magnetic::hysteresis();
    // Iron wall width inside the published 20-40 nm range.
    let fe = rows.iter().find(|r| r.material == "Fe (bcc)" && r.quantity == "δ").unwrap();
    assert!(fe.value >= 20.0 && fe.value <= 40.0, "Fe wall {:.1} nm", fe.value);
    assert!(fe.fold_valid, "Fe wall fold");
    // Permalloy exchange length ≈5 nm.
    let py = rows
        .iter()
        .find(|r| r.material == "Ni80Fe20 (permalloy)" && r.quantity == "l_ex")
        .unwrap();
    close(py.rel_error(), 0.20, "permalloy exchange length");
    // Astroid minimum at 45°.
    let hs = rows.iter().find(|r| r.quantity == "h_sw(45°)/H_K").unwrap();
    close(hs.rel_error(), 5e-3, "SW astroid 45 deg");
    // The macrospin coercivity must be flagged invalid for the bulk magnets.
    for name in ["Fe (bcc)", "Ni80Fe20 (permalloy)"] {
        let hc = rows
            .iter()
            .find(|r| r.material == name && r.quantity == "H_c (macrospin)")
            .unwrap();
        assert!(!hc.fold_valid, "{name} macrospin coercivity must be unfolded");
    }
}

#[test]
fn gold_reflectance_and_sic_reststrahlen() {
    let rows = optical::optical();
    let au = rows.iter().find(|r| r.material == "Au").unwrap();
    close(au.rel_error(), 0.05, "gold R(633 nm)");
    let mid = rows.iter().find(|r| r.quantity == "R(mid-band)").unwrap();
    assert!((mid.value - 1.0).abs() < 0.05, "SiC in-band R = {}", mid.value);
    let to = rows.iter().find(|r| r.quantity == "ω_TO").unwrap();
    let lo = rows.iter().find(|r| r.quantity == "ω_LO").unwrap();
    close(to.rel_error(), 0.02, "SiC TO");
    close(lo.rel_error(), 0.02, "SiC LO");
    assert!(to.value < lo.value, "TO must sit below LO");
    println!("Au R = {:.4} at 633 nm; SiC band {:.1}-{:.1} cm^-1", au.value, to.value, lo.value);
}

#[test]
fn wetting_fold_and_granular_mu_i() {
    let w = interfaces::wetting_contact_angle();
    let wen = w
        .iter()
        .find(|r| r.test == "contact angle (Wenzel)")
        .unwrap();
    // The fine model is within a few degrees of the Wenzel fold at large a/λ.
    close(wen.rel_error(), 0.03, "Wenzel fold");
    assert!(wen.fold_valid, "Wenzel fold should be certified: {}", wen.note);

    let g = interfaces::glass_beads();
    for r in g.iter().filter(|r| r.test == "μ(I) shear cell" && r.quantity.starts_with("μ(")) {
        close(r.rel_error(), 0.30, &format!("DEM {}", r.quantity));
        println!("glass beads {}: mu={:.3} vs Jop {:.3}", r.quantity, r.value, r.reference);
    }
}

#[test]
fn fold_summary_lists_valid_and_invalid_folds() {
    let rows = testlab::report::battery_report();
    let summary = testlab::report::format_fold_summary(&rows);
    assert!(summary.contains("cross-material fold summary"));
    assert!(summary.contains("Cu"));
    assert!(summary.contains("TOTAL"));
    // The report must include one table per material.
    let report = testlab::report::format_report(&rows);
    for name in ["Cu", "Al", "Ar", "Si", "polystyrene", "Fe (bcc)", "Au", "SiC", "glass beads"] {
        assert!(report.contains(name), "report missing {name}");
    }
    println!("{summary}");
}
