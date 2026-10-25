//! Benchmarks for the quasi-harmonic rung of the materials ladder.
//!
//! Every test checks an exact identity, a limiting law (Dulong–Petit,
//! Debye `T³`), a published rare-gas measurement, or a thermodynamic
//! identity; no tolerance is widened to make a model agree with itself.
//! Temperatures, lengths and energies are in the natural units of the
//! model unless a quantity is explicitly converted to SI. The reduced LJ
//! runs use σ = ε = m = 1 and ħ = k_B = 1, so `T` is in units of ε/k_B;
//! for argon the conversions in `quasiharmonic` map these to kelvin,
//! metres and joules.

use materials::elastic::elastic_tensor;
use materials::lattice::Lattice;
use materials::potential::LennardJones;
use materials::quasiharmonic::*;

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

/// Monkhorst–Pack mesh, reciprocal vectors and the DOS normalisation.
#[test]
fn bz_grid_reciprocal_vectors_and_dos_normalisation() {
    let lat = Lattice::fcc(1.0);
    let b = reciprocal_vectors(&lat);
    // a_i · b_j = 2π δ_ij.
    for i in 0..3 {
        for j in 0..3 {
            let d = materials::lattice::dot(lat.prim[i], b[j]);
            let want = if i == j { 2.0 * std::f64::consts::PI } else { 0.0 };
            close(d, want, 1e-10, "reciprocal orthogonality");
        }
    }
    // Even mesh: no Γ point. Odd mesh: Γ present.
    let even = KGrid::monkhorst_pack(&lat, [4, 4, 4]);
    assert_eq!(even.points.len(), 64);
    assert!(even
        .points
        .iter()
        .all(|k| materials::continuum::k_norm(*k) > 1e-9));
    let odd = KGrid::monkhorst_pack(&lat, [3, 3, 3]);
    assert_eq!(odd.points.len(), 27);
    assert!(odd
        .points
        .iter()
        .any(|k| materials::continuum::k_norm(*k) < 1e-12));

    // DOS normalisation: ∫ g dω = 3 per atom (the three Γ zeros carry no
    // measure and are dropped).
    let pot = LennardJones::new(1.0, 1.0);
    let q = Qha::new(lat.clone(), &pot, 2.5, vec![1.0], [8, 8, 8], 1.0, 1.0);
    let h = q.harmonics(q.static_scale);
    let (centres, g) = density_of_states(&h, &pot, &q.grid.points, 80);
    let dw = centres[1] - centres[0];
    let total: f64 = g.iter().sum::<f64>() * dw;
    close(total, 3.0, 0.02, "DOS normalisation");
    assert!(g.iter().all(|x| *x >= 0.0));
}

/// Classical high-temperature limit: `C_V → 3 k_B` per atom (Dulong–Petit).
#[test]
fn dulong_petit_limit() {
    let pot = LennardJones::new(1.0, 1.0);
    let base = Lattice::fcc(1.0);
    let q = Qha::new(base, &pot, 2.5, vec![1.0], [12, 12, 12], 1.0, 1.0);
    let s0 = q.static_scale;
    for t in [80.0, 160.0, 320.0] {
        let cv = q.thermo(s0, t, Statistics::Quantum).cv;
        println!("Dulong-Petit T={t}: Cv={cv:.8} (3-Cv={:.3e})", 3.0 - cv);
        close(cv, 3.0, 0.02, "C_V → 3 k_B");
        // The classical branch is exactly 3 k_B.
        let cvc = q.thermo(s0, t, Statistics::Classical).cv;
        close(cvc, 3.0, 1e-12, "classical C_V = 3 k_B");
    }
}

/// Low-temperature Debye law `C_V = (12π⁴/5) k_B (T/θ_D)³`; extract θ_D
/// and compare with the Debye temperature from the elastic tensor and the
/// directional-average sound speed at the same configuration.
#[test]
fn debye_t3_extracts_the_elastic_debye_temperature() {
    let pot = LennardJones::new(1.0, 1.0);
    let base = Lattice::fcc(1.0);
    let q = Qha::new(base.clone(), &pot, 2.5, vec![1.0], [32, 32, 32], 1.0, 1.0);
    let s0 = q.static_scale;

    // Fit C_V/T³ over T ≈ 0.1 θ_D (where the Brillouin-zone sum resolves
    // the acoustic modes and the Debye law holds).
    let mut sum = 0.0;
    let n = 3;
    for t in [2.5, 3.0, 3.5] {
        let cv = q.thermo(s0, t, Statistics::Quantum).cv;
        let a3 = cv / t.powi(3);
        println!("Debye T={t}: Cv/T3={a3:.8e} theta={:.5}", debye_temperature_from_cv3(a3, 1.0));
        sum += a3;
    }
    let theta_cv = debye_temperature_from_cv3(sum / n as f64, 1.0);

    // Elastic (low-frequency) Debye temperature at the same lattice.
    let c = elastic_tensor(&Lattice::fcc(s0), &pot, 2.5 * s0);
    let nd = 4.0 / (s0 * s0 * s0);
    let theta_el = debye_temperature_elastic(&c, nd, nd, 1.0, 1.0, 40, 80);
    println!("theta_Cv={theta_cv:.5} theta_elastic={theta_el:.5}");

    // The two independent routes must agree to a few percent.
    rel_close(theta_cv, theta_el, 0.05, "theta_D from Cv vs elastic");
    // The Debye law must hold: C_V/T³ nearly constant over the range.
    let a3_a = q.thermo(s0, 2.5, Statistics::Quantum).cv / 2.5_f64.powi(3);
    let a3_b = q.thermo(s0, 3.5, Statistics::Quantum).cv / 3.5_f64.powi(3);
    rel_close(a3_a, a3_b, 0.05, "Debye T^3 plateau");
}

/// Classical-limit QHA: at low `T` the equilibrium volume is linear in `T`
/// (constant α), and the QHA volume exceeds the static minimum.
#[test]
fn classical_qha_expansion_is_linear_at_low_temperature() {
    let pot = LennardJones::new(1.0, 1.0);
    let base = Lattice::fcc(1.0);
    let q = Qha::new(base, &pot, 2.5, vec![1.0], [6, 6, 6], 1.0, 1.0);
    let s0 = q.static_scale;
    let _ = s0;
    let mut alphas = Vec::new();
    for t in [0.02, 0.05, 0.1, 0.15] {
        let p = q.properties(t, Statistics::Classical);
        assert!(p.alpha > 0.0, "positive thermal expansion at T={t}");
        alphas.push((t, p.alpha));
        println!("classical T={t}: s={:.6} alpha={:.6e}", p.scale, p.alpha);
    }
    // α is the low-T plateau of the classical model; it varies slowly.
    let (_, a_lo) = alphas[0];
    let (_, a_hi) = alphas[3];
    rel_close(a_lo, a_hi, 0.5, "classical alpha plateau");

    // Two independent estimates of the thermodynamic Grüneisen parameter
    // agree: α B V / C_V and the mode heat-capacity average of −∂lnω/∂lnV.
    let p = q.properties(0.2, Statistics::Classical);
    rel_close(p.gruneisen, p.gruneisen_modes, 0.05, "Gruneisen alpha vs modes");
}

/// Grüneisen parameter of LJ fcc at zero pressure versus the commonly
/// quoted value for rare-gas solids, and its cutoff insensitivity.
#[test]
fn lj_fcc_gruneisen_parameter() {
    let pot = LennardJones::new(1.0, 1.0);
    let base = Lattice::fcc(1.0);
    for rmax in [2.0, 2.5, 4.0] {
        let q = Qha::new(base.clone(), &pot, rmax, vec![1.0], [6, 6, 6], 1.0, 1.0);
        let p = q.properties(0.3, Statistics::Classical);
        // At high T all modes are equally excited: the plain average.
        let g_modes = q.gruneisen_modes_at(p.scale, 1e3, Statistics::Classical);
        println!("rmax={rmax}: gamma_alpha={:.4} gamma_modes={:.4}", p.gruneisen, g_modes);
        // The literature value for rare-gas/LJ fcc at zero pressure is
        // ~2.6–3.0 (see analysis/quasiharmonic.md); we allow 2.5–3.5.
        assert!(
            (2.5..=3.5).contains(&g_modes),
            "Grüneisen parameter {g_modes} outside 2.5–3.5"
        );
        rel_close(p.gruneisen, g_modes, 0.05, "alpha-B-V/Cv vs mode average");
    }
}

/// Solid argon: zero-point lattice expansion, QHA lattice constant vs
/// temperature against Peterson, Batchelder & Simmons (1966), and the
/// Lindemann/spinodal breakdown approaching melting.
#[test]
fn argon_qha_lattice_constant_and_breakdown() {
    let ar = ArgonModel::new([6, 6, 6], 2.5);
    let aq = ar.qha();
    let sigma = ARGON_SIGMA_M;
    let to_ang = |s: f64| s * sigma * 1e10;

    // Static (classical, no zero-point) equilibrium of the same energy the
    // QHA minimises: golden-section search on static_energy_per_atom.
    let (mut lo, mut hi) = (0.97 * aq.static_scale, 1.03 * aq.static_scale);
    let g = 0.5 * (5f64.sqrt() - 1.0);
    for _ in 0..100 {
        let (x1, x2) = (hi - g * (hi - lo), lo + g * (hi - lo));
        if aq.static_energy_per_atom(x1) < aq.static_energy_per_atom(x2) { hi = x2 } else { lo = x1 }
    }
    let s_static = 0.5 * (lo + hi);
    println!("static minimum of the QHA energy s={s_static:.8}, relax() s={:.8}", aq.static_scale);
    let a_static = to_ang(s_static);
    // Quantum QHA at T → 0: zero-point energy expands the lattice.
    let a0 = to_ang(aq.equilibrium_scale(1e-6, Statistics::Quantum));
    println!(
        "Ar static a={a_static:.4} A, a(T->0)={a0:.4} A, ZPE expansion={:.4} A",
        a0 - a_static
    );
    assert!(
        a0 > a_static,
        "zero-point expansion must be positive: {a0} vs {a_static}"
    );
    // Zero-point expansion of solid argon is about 1% (quantum parameter
    // Λ* ≈ 0.19); the quantitative check is a(T) against experiment below.
    let zpe = (a0 - a_static) / a_static;
    assert!((0.005..0.015).contains(&zpe), "zero-point expansion {zpe:.4} outside 0.5–1.5%");

    // Debye temperature of argon from the elastic tensor and the
    // directional-average sound speed, converted to kelvin. Experimental
    // low-temperature heat-capacity value ≈ 93.5 K (standard rare-gas
    // Debye table; see analysis/quasiharmonic.md).
    let c_ar = elastic_tensor(&Lattice::fcc(aq.static_scale), aq.pot, 2.5 * aq.static_scale);
    let nd_ar = 4.0 / aq.static_scale.powi(3);
    let theta_ar = argon_theta_kelvin(debye_temperature_elastic(
        &c_ar,
        nd_ar,
        nd_ar,
        argon_hbar_natural(),
        1.0,
        40,
        80,
    ));
    println!("Ar Debye temperature from elasticity = {theta_ar:.3} K");
    rel_close(theta_ar, 93.5, 0.02, "Ar Debye temperature");

    // Experimental lattice constants (rounded from Peterson, Batchelder &
    // Simmons 1966; see analysis/quasiharmonic.md). The low-temperature
    // points are the anchors; near melting the QHA is expected to fail.
    let exp = [
        (4.2, 5.311),
        (10.0, 5.312),
        (20.0, 5.317),
        (40.0, 5.342),
        (50.0, 5.362),
        (60.0, 5.389),
        (80.0, 5.462),
    ];
    let mut at50 = None;
    let mut at60 = None;
    for (tk, a_exp) in exp {
        let t = tk / ARGON_EPSILON_KB;
        let p = aq.properties(t, Statistics::Quantum);
        let a = to_ang(p.scale);
        let err = a - a_exp;
        let alpha_si = p.alpha / ARGON_EPSILON_KB;
        let b_si = p.bulk_modulus * ARGON_EPSILON_KB * KB / sigma.powi(3);
        println!(
            "Ar T={tk:4}K a={a:.4} exp={a_exp:.4} err={err:+.4} ({:+.2}%) alpha={alpha_si:.3e}/K B={b_si:.3e}Pa L={:.4} near={}",
            100.0 * err / a_exp,
            p.lindemann,
            p.near_instability
        );
        if tk == 50.0 {
            at50 = Some(p);
        }
        if tk == 60.0 {
            at60 = Some(p);
        }
        // QHA is accurate at low T and degrades systematically on heating.
        if tk <= 20.0 {
            assert!(err.abs() < 0.02, "Ar a({tk} K) error {err}");
        } else if tk <= 40.0 {
            assert!(err.abs() < 0.05, "Ar a({tk} K) error {err}");
        }
        // C_P − C_V = α² B V T holds exactly by construction.
        rel_close(
            p.cp_minus_cv,
            p.alpha * p.alpha * p.bulk_modulus * p.volume * t,
            1e-9,
            "C_P−C_V identity",
        );
        if !p.near_instability {
            // Thermal expansion is positive away from the breakdown.
            assert!(alpha_si > 0.0);
            assert!(p.cp_minus_cv > 0.0 && p.cp_minus_cv < p.cv);
        }
    }

    // Breakdown indicators approaching melting: the Lindemann ratio
    // reaches ≈0.1 and the equilibrium hits the harmonic spinodal.
    let p50 = at50.unwrap();
    let p60 = at60.unwrap();
    println!("Ar 50 K: L={:.4} near={}", p50.lindemann, p50.near_instability);
    println!("Ar 60 K: L={:.4} near={}", p60.lindemann, p60.near_instability);
    assert!(
        (0.08..0.12).contains(&p50.lindemann),
        "Lindemann at 50 K: {}",
        p50.lindemann
    );
    assert!(
        p60.near_instability,
        "60 K should be flagged as quasi-harmonic breakdown"
    );
    // The harmonic spinodal in argon sits below the experimental triple
    // point, so the QHA must not be trusted through melting.
    assert!(to_ang(aq.spinodal_scale) < 5.462, "spinodal above exp a(80 K)");
}
