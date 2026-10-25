//! Summary report: every test at every level, grouped by material, plus a
//! cross-material summary of which folds were certified.

use crate::dynamic::plate_impact;
use crate::{
    elastic, interfaces, magnetic, materials, mechanical, optical, thermal, viscoelastic,
    Measurement,
};

/// Run the cheap static and ultrasonic tests for the three elastic materials.
pub fn static_report() -> Vec<Measurement> {
    let mut out = Vec::new();
    for m in materials() {
        out.extend(elastic::static_elastic(&m));
        out.extend(elastic::ultrasonic(&m));
        // Rayleigh case: 10 µm grains at 5 MHz; geometric case: 10 mm grains.
        out.push(elastic::scattering_fold(&m, 1.0e-5, 5.0e6));
        out.push(elastic::scattering_fold(&m, 1.0e-2, 5.0e6));
    }
    out
}

/// The full virtual-test battery: the mechanical, thermal, viscoelastic,
/// magnetic, optical, wetting and granular standard tests on top of the
/// static/ultrasonic suite. Heavy paths (quasi-harmonic argon, the small DEM
/// shear cell) are cached internally, so calling this is cheap after the
/// first call.
pub fn battery_report() -> Vec<Measurement> {
    let mut out = static_report();
    // Mechanical: fracture for the elastic materials, yield where available.
    for m in materials() {
        out.extend(mechanical::fracture(&m));
        out.extend(mechanical::metal_yield(&m));
    }
    // Thermal.
    out.extend(thermal::argon_thermal());
    out.extend(thermal::silicon_thermal());
    // Viscoelastic.
    out.extend(viscoelastic::dma());
    // Magnetic.
    out.extend(magnetic::hysteresis());
    // Optical.
    out.extend(optical::optical());
    // Interfaces / granular.
    out.extend(interfaces::wetting_contact_angle());
    out.extend(interfaces::glass_beads());
    out
}

/// The complete report including the MD state point and the plate impacts.
pub fn full_report() -> Vec<Measurement> {
    let mut out = battery_report();
    out.extend(crate::atomistic::argon_fold(1500, 2500, 11));
    for name in ["Cu", "Al"] {
        let m = materials().into_iter().find(|m| m.name == name).unwrap();
        out.extend(plate_impact(&m, 200.0, 600, 0.4).measurements());
    }
    out
}

/// The three elastic materials plus the new substances of the battery, in the
/// order their tables should appear.
const MATERIAL_ORDER: [&str; 11] = [
    "Cu",
    "Al",
    "Ar",
    "Si",
    "polystyrene",
    "Fe (bcc)",
    "Ni80Fe20 (permalloy)",
    "Au",
    "SiC",
    "water (textured)",
    "glass beads",
];

fn material_note(name: &str) -> String {
    match name {
        "Si" => "silicon: diamond lattice, a = 0.5431 nm; the phonon heat-transport tests use the \
                Morelli–Kazan parameter set".to_string(),
        "polystyrene" => "polystyrene-like melt: T_g = 373 K, G_N⁰ ≈ 2.0×10⁵ Pa (Ferry 1980)".to_string(),
        "Fe (bcc)" => "bcc α-iron: A = 2.1e-11 J/m, M_s = 1.71e6 A/m, K₁ = 4.8e4 J/m³ (Coey 2010)"
            .to_string(),
        "Ni80Fe20 (permalloy)" => {
            "permalloy: A = 1.3e-11 J/m, M_s = 8.0e5 A/m, K ≈ 0 (Coey 2010)".to_string()
        }
        "Au" => "gold: Rakić et al. (1998) Lorentz–Drude fit for the optical constants".to_string(),
        "SiC" => "3C silicon carbide: ω_TO = 796.2, ω_LO = 972.7 cm⁻¹, ε_∞ = 6.52 (Palik 1985)"
            .to_string(),
        "water (textured)" => {
            "water at 20 °C: γ = 72.8 mN/m, ρ = 998.21 kg/m³ (Vargaftik et al. 1983; CRC)"
                .to_string()
        }
        "glass beads" => "glass beads: μ(I) constants μ_s = tan 20.9°, μ_2 = tan 32.8°, I₀ = 0.279 \
                          (Jop et al. 2006)".to_string(),
        _ => {
            // The elastic materials carry their bottom-rung note.
            match materials().into_iter().find(|m| m.name == name) {
                Some(m) => match &m.bottom {
                    crate::BottomRung::PairPotential { potential, .. } => {
                        format!("atomistic model: {potential} (self-consistent)")
                    }
                    crate::BottomRung::MeasuredSingleCrystal { note } => {
                        format!("atomistic model: {note}")
                    }
                },
                None => String::new(),
            }
        }
    }
}

/// A fixed-width table for one group of measurements.
pub fn format_table(rows: &[Measurement]) -> String {
    let mut s = String::new();
    s.push_str(&format!(
        "{:<28} {:<16} {:<24} {:>12} {:>12} {:>10} {:<6}\n",
        "test", "level", "quantity", "value", "reference", "error", "fold"
    ));
    s.push_str(&"-".repeat(116));
    s.push('\n');
    for r in rows {
        let err = if r.reference.abs() > 1e-12 {
            format!("{:+.2e}", r.rel_error())
        } else {
            format!("{:+.2e}", r.error())
        };
        s.push_str(&format!(
            "{:<28} {:<16} {:<24} {:>12.5} {:>12.5} {:>10} {:<6}\n",
            truncate(&r.test, 28),
            r.level.label(),
            truncate(&r.quantity, 24),
            r.value,
            r.reference,
            err,
            if r.fold_valid { "ok" } else { "UNFOLD" }
        ));
        if !r.note.is_empty() {
            s.push_str(&format!("    note: {}\n", r.note));
        }
    }
    s
}

fn truncate(s: &str, n: usize) -> String {
    if s.chars().count() <= n {
        s.to_string()
    } else {
        s.chars().take(n.saturating_sub(1)).collect::<String>() + "…"
    }
}

/// Material names present in `rows`, ordered by the canonical list, then any
/// unknown names in first-appearance order.
fn ordered_materials(rows: &[Measurement]) -> Vec<String> {
    let mut present: Vec<String> = Vec::new();
    for r in rows {
        if !present.contains(&r.material) {
            present.push(r.material.clone());
        }
    }
    let mut out: Vec<String> = Vec::new();
    for name in MATERIAL_ORDER {
        if !name.is_empty() && present.iter().any(|p| p == name) && !out.contains(&name.to_string()) {
            out.push(name.to_string());
        }
    }
    for p in present {
        if !out.contains(&p) {
            out.push(p);
        }
    }
    out
}

/// The whole report as a printable string: one table per material, then the
/// cross-material fold summary.
pub fn format_report(rows: &[Measurement]) -> String {
    let mut s = String::new();
    for name in ordered_materials(rows) {
        let group: Vec<Measurement> = rows
            .iter()
            .filter(|r| r.material == name)
            .cloned()
            .collect();
        if group.is_empty() {
            continue;
        }
        s.push_str(&format!("\n=== {} ===\n", name));
        let note = material_note(&name);
        if !note.is_empty() {
            s.push_str(&format!("{note}\n"));
        }
        s.push_str(&format_table(&group));
    }
    s.push_str(&format_fold_summary(rows));
    s
}

/// Cross-material summary: how many folds were certified, and which tests
/// must be unfolded.
pub fn format_fold_summary(rows: &[Measurement]) -> String {
    let mut s = String::new();
    s.push_str("\n=== cross-material fold summary ===\n");
    s.push_str(&format!(
        "{:<28} {:>6} {:>7} {:>9}  {}\n",
        "material", "tests", "valid", "invalid", "folds to unfold"
    ));
    s.push_str(&"-".repeat(116));
    s.push('\n');
    let mut total = (0usize, 0usize);
    for name in ordered_materials(rows) {
        let group: Vec<&Measurement> = rows.iter().filter(|r| r.material == name).collect();
        if group.is_empty() {
            continue;
        }
        let valid = group.iter().filter(|r| r.fold_valid).count();
        let invalid = group.len() - valid;
        total.0 += valid;
        total.1 += invalid;
        let unfolded: Vec<String> = group
            .iter()
            .filter(|r| !r.fold_valid)
            .map(|r| format!("{}:{}", r.test, r.quantity))
            .collect();
        let list = if unfolded.is_empty() {
            "—".to_string()
        } else {
            unfolded.join(", ")
        };
        s.push_str(&format!(
            "{:<28} {:>6} {:>7} {:>9}  {}\n",
            name,
            group.len(),
            valid,
            invalid,
            truncate(&list, 78)
        ));
    }
    s.push_str(&"-".repeat(116));
    s.push('\n');
    s.push_str(&format!(
        "{:<28} {:>6} {:>7} {:>9}\n",
        "TOTAL",
        total.0 + total.1,
        total.0,
        total.1
    ));
    s
}
