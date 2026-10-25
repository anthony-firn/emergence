//! Solid–solid transitions (allotropes) as extra solid phases, for the
//! isobaric phase-change model (`isobaric::PhaseSet`).
//!
//! **Iron.** At 1 atm iron is bcc (α) up to 1185 K, fcc (γ) to 1667 K, bcc
//! again (δ) to its melting point 1811 K. α and δ are one crystal structure,
//! so one bcc phase must be stable below *and* above the fcc phase:
//! g_bcc − g_fcc changes sign twice.
//!
//! The bcc and fcc phases now carry the full temperature-dependent thermal
//! EOS of Dorogokupets, Dymshits, Litasov & Sokolova, Sci. Rep. 7, 41863
//! (2017), Table 1: a Vinet cold curve, an Einstein thermal term with an
//! Altshuler Grüneisen parameter, and a free-electron term (`thermal_eos.rs`).
//! Volume, bulk modulus, thermal expansion, the temperature dependence of
//! c_p and the electronic term are therefore taken from the published EOS.
//! The bcc phase also carries the Inden / Hillert–Jarl magnetic Gibbs energy
//! (CALPHAD; Dinsdale 1991) with T_C = 1043 K, bcc p = 0.4 and β = 2.22 μ_B.
//!
//! Only the reference of the two phases is fitted, exactly as the NASG model
//! fitted bcc c_v and q, q′: a gauge q − q′T plus a constant c_v offset is
//! chosen so that (i) g_bcc = g_fcc at 1185 K and at 1667 K and (ii)
//! ΔH(α→γ) = 0.90 kJ/mol. The volume jumps and ΔH(γ→δ) are then predictions.
//!
//! The published bcc/fcc `U₀` and the isobaric `Magnetic` form differ from
//! the paper's magnetic term only by a quantity linear in T, which the gauge
//! absorbs; the heat capacity is the paper's λ-shaped one either way.
//!
//! **Tin.** Grey tin (α, diamond cubic, 5.77 g/cm³) turns into white tin
//! (β, 7.27 g/cm³) at 286.35 K on heating (2.09 kJ/mol); the reverse
//! transformation on cooling ("tin pest") is famously slow — it needs
//! nucleation, which the equilibrium model does not delay.

use crate::elements::P_ATM;
use crate::eos::Nasg;
use crate::isobaric::{enthalpy, Magnetic, PhaseSet};
use crate::thermal_eos::{ThermalEos, T_REF};

const R: f64 = 8.314_462_618;

/// NASG phase from density, sound speed, c_p and volumetric expansion at
/// (T, p) (same construction as the element fits).
fn nasg_fit(rho: f64, c: f64, cp: f64, alpha: f64, t: f64, p: f64) -> Nasg {
    let v = 1.0 / rho;
    let gamma = 1.0 + c * c * alpha * alpha * t / cp;
    let cv = cp / gamma;
    let p_inf = (gamma - 1.0) * cv / (alpha * v) - p;
    Nasg { gamma, p_inf, b: v * (1.0 - alpha * t), q: 0.0, cv, q_prime: 0.0 }
}

/// NASG phase with a prescribed specific volume, volumetric expansion and
/// bulk sound speed at (t, p) and a prescribed c_v. The last NASG parameter
/// is free: from c² = γ(γ−1)c_v/(α²t) and γ(p+p∞) = c²αt/v,
///   b = v(1 − αt),  γ = (1 + √(1 + 4c²α²t/c_v))/2,  p∞ = c²αt/(γv) − p,
/// so v, α and c are reproduced exactly whatever c_v is. Choosing c_v below
/// therefore fixes the mechanical EOS (and the sound speed) from measured
/// data, leaving c_v to carry the latent-heat fit.
fn nasg_fit_cv(v: f64, c: f64, alpha: f64, cv: f64, t: f64, p: f64) -> Nasg {
    let b = v * (1.0 - alpha * t);
    let gg = c * c * alpha * alpha * t / cv;
    let gamma = 0.5 * (1.0 + (1.0 + 4.0 * gg).sqrt());
    let p_inf = c * c * alpha * t / (v * gamma) - p;
    Nasg { gamma, p_inf, b, q: 0.0, cv, q_prime: 0.0 }
}

/// Gibbs energy of a NASG phase plus an optional magnetic term.
fn g_of(e: &Nasg, mag: Option<Magnetic>, p: f64, t: f64) -> f64 {
    e.gibbs(p, t) + mag.map_or(0.0, |m| m.gibbs(t))
}

/// Specific enthalpy of a NASG phase plus an optional magnetic term.
fn h_of(e: &Nasg, mag: Option<Magnetic>, p: f64, t: f64) -> f64 {
    enthalpy(e, p, t) + mag.map_or(0.0, |m| m.enthalpy(t))
}

/// q, q′ of phase `a` so that its Gibbs energy equals that of `b` at T1 and
/// T2 (g is linear in q and q′: g = ĝ + q − q′T) including any magnetic
/// contributions.
fn match_two(a: &mut Nasg, mag_a: Option<Magnetic>, b: &Nasg, mag_b: Option<Magnetic>, p: f64, t1: f64, t2: f64) {
    a.q = 0.0;
    a.q_prime = 0.0;
    let d1 = g_of(b, mag_b, p, t1) - g_of(a, mag_a, p, t1);
    let d2 = g_of(b, mag_b, p, t2) - g_of(a, mag_a, p, t2);
    // q − q′ t = d
    a.q_prime = (d1 - d2) / (t2 - t1);
    a.q = d1 + a.q_prime * t1;
}

/// q, q′ of phase `a` so that it coexists with `b` at T with latent heat
/// h_a − h_b = dh, including any magnetic contributions.
fn match_one(a: &mut Nasg, mag_a: Option<Magnetic>, b: &Nasg, mag_b: Option<Magnetic>, p: f64, t: f64, dh: f64) {
    a.q = 0.0;
    a.q_prime = 0.0;
    a.q = h_of(b, mag_b, p, t) + dh - h_of(a, mag_a, p, t);
    a.q_prime = (g_of(a, mag_a, p, t) - g_of(b, mag_b, p, t)) / t;
}

/// Anchor a NASG phase `a` (whose ungauged Gibbs/enthalpy at T are `g_a0`
/// and `h_a0`) to a target Gibbs/enthalpy `(g_b, h_b)` at T with latent heat
/// h_a − h_b = dh. Works for any target phase, thermal or NASG.
fn anchor(a: &mut Nasg, t: f64, g_a0: f64, h_a0: f64, g_b: f64, h_b: f64, dh: f64) {
    a.q = 0.0;
    a.q_prime = 0.0;
    a.q = h_b + dh - h_a0;
    // g_a0 is the ungauged Gibbs; after the shift a.gibbs = g_a0 + q − q′t,
    // so q′ = (g_a0 + q − g_b)/t (a.q is already set).
    a.q_prime = (g_a0 + a.q - g_b) / t;
}

/// bcc-Fe thermal EOS, Dorogokupets et al. 2017, Table 1.
fn bcc_thermal(m: f64) -> ThermalEos {
    ThermalEos {
        v0: 7.092e-6,
        k0: 164.0e9,
        k0p: 5.50,
        theta0: 303.0,
        gamma0: 1.736,
        beta: 1.125,
        ginf: 0.0,
        e0: 198e-6,
        ge: 1.0,
        u0: 0.0,
        q: 0.0,
        q_prime: 0.0,
        cv_off: 0.0,
        t_max: 1811.0,
        mm: m,
    }
}

/// fcc-Fe thermal EOS, Dorogokupets et al. 2017, Table 1.
fn fcc_thermal(m: f64) -> ThermalEos {
    ThermalEos {
        v0: 6.9285e-6,
        k0: 146.2e9,
        k0p: 4.67,
        theta0: 222.5,
        gamma0: 2.203,
        beta: 0.01,
        ginf: 0.0,
        e0: 198e-6,
        ge: 1.0,
        u0: 4470.0,
        q: 0.0,
        q_prime: 0.0,
        cv_off: 0.0,
        t_max: 1811.0,
        mm: m,
    }
}

/// Fit the bcc reference gauge (q, q′, c_v offset) so that the α→γ and γ→δ
/// transitions sit at `t1` and `t2` and ΔH(α→γ) = `dh` J/mol. The published
/// EOS shape (volume, c_p, electronic and magnetic terms) is untouched.
fn fit_gauge(bcc: &mut ThermalEos, fcc: &ThermalEos, mag: &Magnetic, p: f64, t1: f64, t2: f64, dh: f64) {
    let gb1 = bcc.gibbs(p, t1) + mag.gibbs(t1) * bcc.mm;
    let gb2 = bcc.gibbs(p, t2) + mag.gibbs(t2) * bcc.mm;
    let hb1 = bcc.enthalpy(p, t1) + mag.enthalpy(t1) * bcc.mm;
    let gf1 = fcc.gibbs(p, t1);
    let gf2 = fcc.gibbs(p, t2);
    let hf1 = fcc.enthalpy(p, t1);
    // Gauge: G += q − q′T + c[(T−T₀) − T ln(T/T₀)], H += q + c(T−T₀).
    let f = |t: f64| (t - T_REF) - t * (t / T_REF).ln();
    let m = [[1.0, -t1, f(t1)], [1.0, -t2, f(t2)], [1.0, 0.0, t1 - T_REF]];
    let rhs = [gf1 - gb1, gf2 - gb2, hf1 - hb1 - dh];
    let det = |a: [[f64; 3]; 3]| {
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    };
    let d = det(m);
    let dq = det([[rhs[0], m[0][1], m[0][2]], [rhs[1], m[1][1], m[1][2]], [rhs[2], m[2][1], m[2][2]]]);
    let dqp = det([[m[0][0], rhs[0], m[0][2]], [m[1][0], rhs[1], m[1][2]], [m[2][0], rhs[2], m[2][2]]]);
    let dc = det([[m[0][0], m[0][1], rhs[0]], [m[1][0], m[1][1], rhs[1]], [m[2][0], m[2][1], rhs[2]]]);
    bcc.q = dq / d;
    bcc.q_prime = dqp / d;
    bcc.cv_off = dc / d;
}

/// The iron model and the predicted γ→δ latent heat (J/mol).
pub fn iron() -> (PhaseSet, f64) {
    let m = 55.845e-3;
    let p = P_ATM;
    let (t_ag, t_gd, t_m, t_b) = (1185.0, 1667.0, 1811.0, 3134.0);
    let dh_ag = 900.0; // J/mol, α → γ on heating
    let mag = Magnetic { tc: 1043.0, p: 0.4, beta: 2.22, mm: m };

    // Full temperature-dependent solid EOSs (published), gauge fitted.
    let mut bcc_th = bcc_thermal(m);
    let fcc_th = fcc_thermal(m);
    fit_gauge(&mut bcc_th, &fcc_th, &mag, p, t_ag, t_gd, dh_ag);

    // Legacy NASG solids, kept only for code that reads `PhaseSet::eos`
    // directly (the phase-change example and the 2D solver); the phase-set
    // methods above dispatch to the thermal EOS for these two phases.
    let fcc = nasg_fit(7630.0, 3800.0, 34.0 / m, 3.0 * 2.3e-5, t_ag, p);
    let (v_bcc, c_bcc, alpha_bcc) = (1.0 / 7570.0, 4314.0, 5.34e-5);
    let make_bcc = |cv: f64| -> Nasg {
        let mut b = nasg_fit_cv(v_bcc, c_bcc, alpha_bcc, cv, t_ag, p);
        match_two(&mut b, Some(mag), &fcc, None, p, t_ag, t_gd);
        b
    };
    let dh_at = |b: &Nasg, t: f64| h_of(&fcc, None, p, t) - h_of(b, Some(mag), p, t);
    let (mut lo, mut hi) = (fcc.cv, 3.0 * fcc.cv);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if dh_at(&make_bcc(mid), t_ag) < dh_ag / m {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let bcc = make_bcc(0.5 * (lo + hi));

    // liquid: ρ 6.98 g/cm³ at the melting point, c 3820 m/s, c_p 46 J/(mol K);
    // melts from δ (bcc) at 1811 K with 13.81 kJ/mol
    let liquid = nasg_fit(6980.0, 3820.0, 46.0 / m, 1.2e-4, t_m, p);
    // vapour: monatomic ideal gas, boils at 3134 K with 340 kJ/mol
    let cv_v = 1.5 * R / m;
    let vapor = Nasg { gamma: 5.0 / 3.0, p_inf: 0.0, b: 0.0, q: 0.0, cv: cv_v, q_prime: 0.0 };

    let mut set = PhaseSet {
        names: vec!["bcc (α/δ)".into(), "fcc (γ)".into(), "liquid".into(), "vapour".into()],
        eos: vec![bcc, fcc, liquid, vapor],
        liquid: 2,
        vapor: 3,
        t_crit: Some(9250.0),
        magnetic: vec![Some(mag), None, None, None],
        thermal: vec![Some(bcc_th), Some(fcc_th), None, None],
    };

    // Anchor the liquid to the actual (thermal) δ-bcc at the melting point,
    // then the vapour to the liquid at the boiling point.
    let (g_l0, h_l0) = (set.eos[2].gibbs(p, t_m), enthalpy(&set.eos[2], p, t_m));
    let (g_b, h_b) = (set.gibbs(0, p, t_m), set.enthalpy(0, p, t_m));
    anchor(&mut set.eos[2], t_m, g_l0, h_l0, g_b, h_b, 13810.0 / m);
    let g_l = set.gibbs(2, p, t_b);
    let h_l = set.enthalpy(2, p, t_b);
    let (g_v0, h_v0) = (set.eos[3].gibbs(p, t_b), enthalpy(&set.eos[3], p, t_b));
    anchor(&mut set.eos[3], t_b, g_v0, h_v0, g_l, h_l, 340e3 / m);
    let dh_gd_pred = (set.enthalpy(0, p, t_gd) - set.enthalpy(1, p, t_gd)) * m; // h_δ − h_γ at 1667 K
    (set, dh_gd_pred)
}

/// The tin model: grey (α) and white (β) solids, liquid, vapour.
pub fn tin() -> PhaseSet {
    let m = 118.71e-3;
    let p = P_ATM;
    let (t_ab, t_m, t_b) = (286.35, 505.08, 2875.0);
    // white tin: ρ 7.265 g/cm³, c_p 27.1 J/(mol K), linear expansion 2.2e-5,
    // bulk sound speed √(K/ρ) with K ≈ 58 GPa
    let white = nasg_fit(7265.0, 2830.0, 27.1 / m, 3.0 * 2.2e-5, 298.15, p);
    // grey tin: ρ 5.769 g/cm³, c_p 25.8 J/(mol K), linear expansion 4.7e-6,
    // K ≈ 42.6 GPa
    let mut grey = nasg_fit(5769.0, 2720.0, 25.8 / m, 3.0 * 4.7e-6, t_ab, p);
    // grey is the low-temperature phase: h_grey − h_white = −2.09 kJ/mol
    match_one(&mut grey, None, &white, None, p, t_ab, -2090.0 / m);
    let mut liquid = nasg_fit(6990.0, 2464.0, 30.0 / m, 1.0e-4, t_m, p);
    match_one(&mut liquid, None, &white, None, p, t_m, 7030.0 / m);
    let cv_v = 1.5 * R / m;
    let mut vapor = Nasg { gamma: 5.0 / 3.0, p_inf: 0.0, b: 0.0, q: 0.0, cv: cv_v, q_prime: 0.0 };
    match_one(&mut vapor, None, &liquid, None, p, t_b, 296.1e3 / m);
    PhaseSet {
        names: vec!["grey (α)".into(), "white (β)".into(), "liquid".into(), "vapour".into()],
        eos: vec![grey, white, liquid, vapor],
        liquid: 2,
        vapor: 3,
        t_crit: None,
        magnetic: vec![None; 4],
        thermal: vec![None; 4],
    }
}
