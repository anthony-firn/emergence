//! Allotrope tests: the temperature-dependent thermal EOS of bcc and fcc
//! iron and the Inden / Hillert–Jarl magnetic contribution, and the
//! unchanged tin model.

use compressible::allotropes;
use compressible::elements::P_ATM;
use compressible::isobaric::{Ladder, Magnetic, PhaseSet};

const R: f64 = 8.314_462_618;
const M: f64 = 55.845e-3;

/// The bcc thermal EOS reproduces the published ambient-pressure data it is
/// taken from: at 1185 K, V = 7.377 cm³/mol (ρ = 7.570 g/cm³), volumetric
/// expansion 5.34e-5 /K and bulk sound speed 4314 m/s (Dorogokupets et al.
/// 2017, Table 1). fcc at 1185 K is 7.320 cm³/mol (ρ = 7.629 g/cm³).
#[test]
fn thermal_eos_reproduces_published_ambient_data() {
    let (iron, _) = allotropes::iron();
    let p = P_ATM;
    let t = 1185.0;
    let vb = iron.specific_volume(0, p, t) * M;
    let vf = iron.specific_volume(1, p, t) * M;
    assert!((vb * 1e6 - 7.377).abs() < 0.002, "bcc V = {} cm³/mol", vb * 1e6);
    assert!((vf * 1e6 - 7.320).abs() < 0.002, "fcc V = {} cm³/mol", vf * 1e6);
    let b = iron.thermal[0].unwrap();
    let f = iron.thermal[1].unwrap();
    let alpha_b = b.alpha(p, t);
    let alpha_f = f.alpha(p, t);
    assert!((alpha_b * 1e6 - 53.36).abs() < 0.1, "bcc α = {} /K", alpha_b);
    assert!((alpha_f * 1e6 - 70.67).abs() < 0.1, "fcc α = {} /K", alpha_f);
    let cb = {
        // K_S uses the total c_v, including the magnetic part (which is a
        // function of T alone, so c_p = c_v for it).
        let v = b.volume(p, t);
        let kt = b.k_t(p, t);
        let alpha = b.alpha(p, t);
        let cv = b.cv(v, t) - b.cv_off + iron.magnetic[0].unwrap().cp(t) * M;
        // the published 4314 m/s is for the raw EOS; the fitted c_v offset
        // (+2.345 J/mol K, for ΔH(α→γ)) lowers it by 0.2%, to 4305 m/s
        (kt + v * t * (alpha * kt).powi(2) / cv).sqrt() / (M / v).sqrt()
    };
    assert!((cb - 4314.0).abs() < 5.0, "bcc c = {cb} m/s");
    // the bcc heat capacity has the magnetic peak at the Curie point
    let cp_b = |t: f64| iron.cp(0, p, t) * M;
    println!("bcc cp: 300 {:.2}, 1000 {:.2}, 1043 {:.2}, 1185 {:.2}, 1667 {:.2}", cp_b(300.0), cp_b(1000.0), cp_b(1043.0), cp_b(1185.0), cp_b(1667.0));
    println!("fcc cp: 1185 {:.2}, 1667 {:.2}", iron.cp(1, p, 1185.0) * M, iron.cp(1, p, 1667.0) * M);
    assert!(cp_b(1043.0) > cp_b(300.0) + 20.0, "no magnetic peak");
    assert!(cp_b(1043.0) > cp_b(1667.0) + 10.0, "peak decays");
}

/// The fit still reproduces the three calibration data: both transition
/// temperatures and ΔH(α→γ) = 0.90 kJ/mol. ΔH(γ→δ) is a prediction.
#[test]
fn iron_fit_reproduces_the_three_data() {
    let (iron, dh_gd) = allotropes::iron();
    let lad = Ladder::new(&iron, P_ATM, 250.0, 4000.0);
    let tt: Vec<f64> = lad.transitions.iter().map(|x| x.0).collect();
    assert_eq!(tt.len(), 4, "transitions {tt:?}");
    for (got, want) in tt.iter().zip(&[1185.0, 1667.0, 1811.0, 3134.0]) {
        assert!((got - want).abs() < 0.5, "transition {got} vs {want}");
    }
    let dh_ag = (iron.enthalpy(1, P_ATM, 1185.0) - iron.enthalpy(0, P_ATM, 1185.0)) * M;
    assert!((dh_ag - 900.0).abs() < 1.0, "ΔH(α→γ) = {dh_ag} J/mol");
    let dh_gd_direct = (iron.enthalpy(0, P_ATM, 1667.0) - iron.enthalpy(1, P_ATM, 1667.0)) * M;
    assert!((dh_gd - dh_gd_direct).abs() < 1e-9, "reported prediction inconsistent");
    // genuinely not the measured 0.84 kJ/mol and not fitted
    assert!(dh_gd.is_finite() && dh_gd > 0.0 && dh_gd < 900.0);
}

/// The volume jumps and ΔH(γ→δ) are predictions of the published EOS; the
/// measured values are −0.98%, +0.53% and 0.84 kJ/mol. The γ→δ volume is
/// essentially zero, a genuine miss (see analysis §15).
#[test]
fn iron_predictions() {
    let (iron, dh_gd) = allotropes::iron();
    let p = P_ATM;
    let dv_ag = iron.specific_volume(1, p, 1185.0) / iron.specific_volume(0, p, 1185.0) - 1.0;
    let dv_gd = iron.specific_volume(0, p, 1667.0) / iron.specific_volume(1, p, 1667.0) - 1.0;
    println!("ΔV/V(α→γ) = {:.3}%, ΔV/V(γ→δ) = {:.3}%, ΔH(γ→δ) = {:.1} J/mol", dv_ag * 100.0, dv_gd * 100.0, dh_gd);
    assert!((dv_ag * 100.0 + 0.77).abs() < 0.02, "ΔV/V(α→γ) = {}", dv_ag * 100.0);
    assert!((dv_gd * 100.0 + 0.031).abs() < 0.02, "ΔV/V(γ→δ) = {}", dv_gd * 100.0);
    assert!((dh_gd - 554.1).abs() < 2.0, "ΔH(γ→δ) = {dh_gd} J/mol");
}

/// The thermal EOS is thermodynamically consistent: h = g − T ∂g/∂T and
/// c_p = ∂h/∂T.
#[test]
fn thermal_eos_consistency() {
    let (iron, _) = allotropes::iron();
    for k in [0, 1] {
        let e = iron.thermal[k].unwrap();
        for &t in &[300.0, 800.0, 1185.0, 1667.0, 1811.0] {
            let dt = 1e-4 * t;
            let dgdt = (iron.gibbs(k, P_ATM, t + dt) - iron.gibbs(k, P_ATM, t - dt)) / (2.0 * dt);
            let h = iron.gibbs(k, P_ATM, t) - t * dgdt;
            assert!((h - iron.enthalpy(k, P_ATM, t)).abs() < 1e-4 * (h.abs() + 1.0), "h mismatch phase {k} at {t}");
            let dhdt = (iron.enthalpy(k, P_ATM, t + dt) - iron.enthalpy(k, P_ATM, t - dt)) / (2.0 * dt);
            assert!((dhdt - iron.cp(k, P_ATM, t)).abs() < 1e-4 * (dhdt.abs() + 1.0), "c_p mismatch phase {k} at {t}");
            let _ = e;
        }
    }
}

/// The bcc magnetic term is carried into the Gibbs energy, enthalpy and
/// heat capacity, and h = g − T ∂g/∂T holds.
#[test]
fn magnetic_gibbs_enthalpy_consistency() {
    let m = Magnetic { tc: 1043.0, p: 0.4, beta: 2.22, mm: 55.845e-3 };
    for &t in &[300.0, 900.0, 1000.0, 1100.0, 1185.0, 1667.0, 2500.0] {
        let dt = 1e-4 * t;
        let dgdt = (m.gibbs(t + dt) - m.gibbs(t - dt)) / (2.0 * dt);
        let h = m.gibbs(t) - t * dgdt;
        assert!((h - m.enthalpy(t)).abs() < 1e-6 * (m.enthalpy(t).abs() + 1.0), "h mismatch at {t}: {} vs {}", h, m.enthalpy(t));
        let dhdt = (m.enthalpy(t + dt) - m.enthalpy(t - dt)) / (2.0 * dt);
        assert!((dhdt - m.cp(t)).abs() < 1e-5 * (m.cp(t).abs() + 1.0), "c_p mismatch at {t}: {} vs {}", dhdt, m.cp(t));
    }
}

/// The magnetic heat capacity peaks near the Curie temperature and decays on
/// both sides; the paramagnetic tail is small.
#[test]
fn magnetic_heat_capacity_peak_at_curie_point() {
    let m = Magnetic { tc: 1043.0, p: 0.4, beta: 2.22, mm: 55.845e-3 };
    let cp = |t: f64| m.cp(t) * m.mm; // J/(mol K)
    // the excess magnetic entropy is R ln(β+1)
    let expected_entropy = R * (m.beta + 1.0).ln();
    // numerical ∫ c_p/T dT from the low-T plateau to the high-T tail
    let n = 200_000;
    let (lo, hi): (f64, f64) = (1.0, 20_000.0);
    let mut s = 0.0;
    for i in 0..n {
        let t = lo * (hi / lo).powf((i as f64 + 0.5) / n as f64);
        s += cp(t) / t * (t * (hi / lo).ln() / n as f64);
    }
    assert!((s - expected_entropy).abs() < 0.05, "magnetic entropy {s} vs R ln(β+1) {expected_entropy}");
    // peak
    assert!(cp(1043.0) > 20.0, "peak too small: {}", cp(1043.0));
    assert!(cp(1043.0) > 3.0 * cp(300.0), "not peaked vs 300 K");
    assert!(cp(1043.0) > 10.0 * cp(2500.0), "not peaked vs 2500 K");
}

/// The bcc thermal phase carries the magnetic term; the other phases do not.
#[test]
fn only_bcc_iron_is_magnetic() {
    let (iron, _) = allotropes::iron();
    assert!(iron.magnetic[0].is_some());
    assert!(iron.magnetic[1..].iter().all(|m| m.is_none()));
    assert!(iron.thermal[0].is_some() && iron.thermal[1].is_some());
    assert!(iron.thermal[2..].iter().all(|e| e.is_none()));
}

/// The two physics anchors are consistent: the bcc and fcc Gibbs energies
/// are equal at the two transition temperatures.
#[test]
fn iron_gibbs_equal_at_transitions() {
    let (iron, _) = allotropes::iron();
    for &t in &[1185.0, 1667.0] {
        assert!((iron.gibbs(0, P_ATM, t) - iron.gibbs(1, P_ATM, t)).abs() < 1e-6 * iron.gibbs(0, P_ATM, t).abs().max(1.0), "g mismatch at {t}");
    }
}

/// Tin keeps its old behaviour: no magnetic term, same transitions and
/// latent heats.
#[test]
fn tin_is_unchanged() {
    let tin: PhaseSet = allotropes::tin();
    assert!(tin.magnetic.iter().all(|m| m.is_none()));
    assert!(tin.thermal.iter().all(|e| e.is_none()));
    let lad = Ladder::new(&tin, P_ATM, 100.0, 4000.0);
    let tt: Vec<f64> = lad.transitions.iter().map(|x| x.0).collect();
    assert!((tt[0] - 286.35).abs() < 0.5, "{tt:?}");
    assert!((tt[1] - 505.08).abs() < 0.5, "{tt:?}");
    let m = 118.71e-3;
    let dh_ab = (tin.enthalpy(1, P_ATM, tt[0]) - tin.enthalpy(0, P_ATM, tt[0])) * m;
    assert!((dh_ab - 2090.0).abs() < 1.0, "ΔH(grey→white) = {dh_ab} J/mol");
    let dh_m = (tin.enthalpy(2, P_ATM, tt[1]) - tin.enthalpy(1, P_ATM, tt[1])) * m;
    assert!((dh_m - 7030.0).abs() < 1.0, "ΔH(white→liquid) = {dh_m} J/mol");
}
