//! Magnetic tests: domain-wall width and coercivity for iron and permalloy,
//! at the micromagnetics and macrospin (Stoner–Wohlfarth) rungs.
//!
//! The coarse models are the continuum Bloch wall
//! `δ = sqrt(A/K), σ = 4 sqrt(AK)` and the single-domain Stoner–Wohlfarth
//! astroid. The validity indicator is the particle size against the
//! single-domain diameter `D_c` (Kittel 1946): only a particle smaller than
//! `D_c` reverses by coherent rotation. Bulk iron and permalloy are far
//! larger than `D_c`, so the macrospin coercivity is reported with the fold
//! flagged invalid; their real coercivity is set by domain-wall pinning.
//!
//! Sources: Coey, *Magnetism and Magnetic Materials* (2010); Cullity &
//! Graham, *Introduction to Magnetic Materials*, 2nd ed.; Kittel (1946);
//! Stoner & Wohlfarth (1948).

use crate::{Level, Measurement};
use magnetism::hysteresis::{astroid_hsw, single_domain_diameter, Macrospin};
use magnetism::materials::{MagneticMaterial, IRON, PERMALLOY};
use magnetism::micromag::{exchange_length, wall_delta, wall_sigma};

/// Measured coercivities and wall widths used as references.
struct MagneticRef {
    /// Published domain-wall width range (m); NaN when undefined.
    delta_lo: f64,
    delta_hi: f64,
    /// Measured bulk coercivity (A/m); NaN when not available.
    hc_measured: f64,
    hc_source: &'static str,
}

fn magnetic_ref(name: &str) -> MagneticRef {
    match name {
        "Fe (bcc)" => MagneticRef {
            // Commonly quoted 20–40 nm wall and a soft-iron coercivity of
            // order 80 A/m (Cullity & Graham; Coey).
            delta_lo: 20.0e-9,
            delta_hi: 40.0e-9,
            hc_measured: 80.0,
            hc_source: "soft iron H_c ≈ 80 A/m (Cullity & Graham, Introduction to Magnetic Materials)",
        },
        _ => MagneticRef {
            // Permalloy's wall width is set by shape/thickness; the
            // magnetocrystalline formula is undefined at K → 0.
            delta_lo: f64::NAN,
            delta_hi: f64::NAN,
            hc_measured: 1.0,
            hc_source: "permalloy H_c ≈ 1 A/m (soft magnetic alloy; Coey 2010)",
        },
    }
}

fn wall_width_measurement(m: &MagneticMaterial) -> Measurement {
    let r = magnetic_ref(m.name);
    if m.anisotropy_k > 0.0 {
        let delta = wall_delta(m.exchange_a, m.anisotropy_k);
        let inside = delta >= r.delta_lo && delta <= r.delta_hi;
        Measurement::new(
            m.name,
            "domain wall",
            Level::Polycrystal,
            "δ",
            delta * 1e9,
            r.delta_lo * 1e9,
            "nm",
            inside,
            &format!(
                "continuum Bloch wall δ = sqrt(A/K) with A = {:.1e} J/m, K = {:.1e} J/m³; the \
                 published range is 20–40 nm, so the fold is valid. The wall energy is σ = 4 sqrt(AK) \
                 = {:.3e} J/m².",
                m.exchange_a,
                m.anisotropy_k,
                wall_sigma(m.exchange_a, m.anisotropy_k)
            ),
        )
    } else {
        // Permalloy: K ≈ 0, so the magnetocrystalline δ diverges. Report the
        // exchange length, which is well defined, against the quoted ≈5 nm.
        let lex = exchange_length(m.exchange_a, m.ms);
        Measurement::new(
            m.name,
            "domain wall",
            Level::Polycrystal,
            "l_ex",
            lex * 1e9,
            5.0,
            "nm",
            (lex * 1e9 / 5.0 - 1.0).abs() < 0.20,
            &format!(
                "K ≈ 0, so δ = sqrt(A/K) is undefined: the magnetocrystalline wall-width fold is \
                 invalid for permalloy and its wall is set by shape/thickness. The exchange length \
                 l_ex = sqrt(2A/(μ₀M_s²)) = {:.2} nm is well defined and matches the quoted ≈5 nm.",
                lex * 1e9
            ),
        )
    }
}

fn coercivity_measurement(m: &MagneticMaterial) -> Measurement {
    let r = magnetic_ref(m.name);
    let hk = m.anisotropy_field();
    // Stoner–Wohlfarth single-particle switching field at 45°: the astroid
    // minimum h_sw = 1/2 in units of H_K. Measured by the LLG integrator.
    let h_sw = Macrospin::switching_field_llg(std::f64::consts::FRAC_PI_4, 1.0, 0.02, 40_000);
    let h_sw_exact = astroid_hsw(std::f64::consts::FRAC_PI_4);
    let hc_sw = h_sw * hk;
    let dc = single_domain_diameter(m.exchange_a, m.anisotropy_k.max(1.0), m.ms);
    Measurement::new(
        m.name,
        "hysteresis (Stoner–Wohlfarth)",
        Level::Polycrystal,
        "H_c (macrospin)",
        hc_sw,
        r.hc_measured,
        "A/m",
        false,
        &format!(
            "coherent-rotation switching field at 45° (astroid minimum {h_sw_exact:.4} H_K, LLG \
             measured {h_sw:.5} H_K), H_K = {hk:.3e} A/m. The reference is the measured bulk \
             coercivity: {} The macrospin fold is invalid for a bulk sample (D ≫ D_c = {:.1} nm, \
             reversal is by domain-wall motion and pinning), so the two differ by orders of \
             magnitude.",
            r.hc_source,
            dc * 1e9
        ),
    )
}

/// All magnetic measurements for iron and permalloy.
pub fn hysteresis() -> Vec<Measurement> {
    let mut out = Vec::new();
    for m in [IRON, PERMALLOY] {
        out.push(wall_width_measurement(&m));
        out.push(coercivity_measurement(&m));
    }
    // Astroid identity as a self-check of the LLG macrospin rung.
    let h_sw = Macrospin::switching_field_llg(std::f64::consts::FRAC_PI_4, 1.0, 0.02, 40_000);
    out.push(Measurement::new(
        "Fe (bcc)",
        "hysteresis (Stoner–Wohlfarth)",
        Level::Polycrystal,
        "h_sw(45°)/H_K",
        h_sw,
        0.5,
        "",
        (h_sw - 0.5).abs() < 5e-3,
        "LLG macrospin switching field at the astroid minimum vs the exact Stoner–Wohlfarth value \
         1/2 (Stoner & Wohlfarth 1948)",
    ));
    out
}
