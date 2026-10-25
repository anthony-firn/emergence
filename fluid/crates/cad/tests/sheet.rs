//! Sheet metal against the standard bend formulas and against the
//! generated geometry, and springs against the published stiffness.

use std::f64::consts::PI;

use cad::kits::part::Spring;
use cad::kits::sheet::{k_factor, Bend, SheetPart};
use cad::meshing::mesh_solid;

#[test]
fn the_bend_formulas_reproduce_a_worked_example() {
    // 1 mm sheet, 1 mm inside radius, 90 degree bend, K = 0.44: the
    // textbook numbers are BA = 2.262 mm, OSSB = 2.000 mm,
    // BD = 1.738 mm.
    let t = 1e-3;
    let part = SheetPart::new(t, 0.05, vec![0.03, 0.02], vec![Bend::at_right_angle(t)])
        .unwrap()
        .with_k_factor(0.44);
    let n = part.numbers(0).unwrap();
    assert!((n.allowance * 1e3 - 2.2619).abs() < 1e-3, "BA {:.4} mm", n.allowance * 1e3);
    assert!((n.outside_setback * 1e3 - 2.0).abs() < 1e-9, "OSSB {:.4} mm", n.outside_setback * 1e3);
    assert!((n.deduction * 1e3 - 1.7381).abs() < 1e-3, "BD {:.4} mm", n.deduction * 1e3);
    // The neutral axis sits K T above the inside surface.
    assert!((n.neutral_radius - (t + 0.44 * t)).abs() < 1e-18);
    // BD = 2 OSSB - BA identically.
    assert!((n.deduction - (2.0 * n.outside_setback - n.allowance)).abs() < 1e-18);

    // The flat length, both ways: from the apex dimensions less the
    // deductions, and from the straights plus the allowances.
    let a = part.flat_length().unwrap();
    let b = part.flat_length_from_geometry().unwrap();
    assert!((a - b).abs() / a < 1e-15, "{a:e} vs {b:e}");
    assert!((a * 1e3 - (30.0 + 20.0 - 1.7381)).abs() < 1e-3, "flat {:.4} mm", a * 1e3);

    // At K = 0.5 the bend allowance is the mid-surface arc length exactly.
    let half = part.clone().with_k_factor(0.5);
    let m = half.numbers(0).unwrap();
    assert!((m.allowance - 0.5 * PI * (t + 0.5 * t)).abs() < 1e-18);

    // The rule-of-thumb table is monotone in R/T and reaches 0.5.
    assert!(k_factor(0.5) < k_factor(2.0) && k_factor(2.0) < k_factor(5.0));
    assert!((k_factor(5.0) - 0.5).abs() < 1e-18);
}

#[test]
fn the_generated_section_has_the_developed_length() {
    // A U-channel: two 90 degree bends.
    let t = 1.5e-3;
    let part = SheetPart::new(
        t,
        0.06,
        vec![0.025, 0.04, 0.025],
        vec![Bend::at_right_angle(2.0 * t), Bend::at_right_angle(2.0 * t)],
    )
    .unwrap()
    .with_k_factor(0.42);

    // The neutral-axis path, measured off the generated polyline: the
    // chords of the arcs make it short, so refine and watch it converge to
    // the developed length.
    let developed = part.flat_length().unwrap();
    let mut errors = Vec::new();
    for per_bend in [4usize, 16, 64] {
        let axis = part.neutral_axis(per_bend).unwrap();
        let measured: f64 = axis
            .windows(2)
            .map(|w| ((w[1][0] - w[0][0]).powi(2) + (w[1][1] - w[0][1]).powi(2)).sqrt())
            .sum();
        errors.push((developed - measured).abs() / developed);
    }
    println!("neutral-axis chord error: {errors:?} (developed {:.4} mm)", developed * 1e3);
    let order = (errors[0] / errors[2]).log2() / (64f64 / 4.0).log2();
    assert!(order > 1.8, "chord error order {order:.2}");
    assert!(errors[2] < 1e-5, "refined path length off by {:e}", errors[2]);

    // The formed part is a U: its bounding box is the outside dimensions.
    let (lo, hi) = part.formed_bbox(32).unwrap();
    let width = hi[0] - lo[0];
    let height = hi[1] - lo[1];
    println!("formed U: {:.3} x {:.3} mm", width * 1e3, height * 1e3);
    // The run starts along +x for the first 25 mm flange, turns up for the
    // 40 mm web, and turns back for the second flange: so x spans a
    // flange and y the web.
    assert!((width * 1e3 - 25.0).abs() < 0.2, "the flanges should be 25 mm: {:.3}", width * 1e3);
    assert!((height * 1e3 - 40.0).abs() < 0.2, "the web should be 40 mm: {:.3}", height * 1e3);

    // And it meshes.
    let solid = part.solid(32).unwrap();
    let mesh = mesh_solid(&solid, t / 4.0).expect("mesh");
    assert!(mesh.report().watertight, "{:?}", mesh.report());
    let expect = part.formed_volume(32).unwrap();
    let rel = (mesh.report().volume - expect).abs() / expect;
    assert!(rel < 0.05, "meshed volume {:e} vs section area x width {expect:e}", mesh.report().volume);
}

#[test]
fn volume_is_conserved_only_at_k_one_half() {
    // The neutral-axis model stretches the outside of a bend and
    // compresses the inside; the two cancel exactly only when the neutral
    // axis is at mid-thickness. This is a property of the model, and the
    // test pins it rather than hiding it.
    let t = 1e-3;
    let make = |k: f64| {
        SheetPart::new(t, 0.05, vec![0.03, 0.03], vec![Bend::at_right_angle(t)])
            .unwrap()
            .with_k_factor(k)
    };
    let at_half = make(0.5).volume_change(200).unwrap();
    assert!(at_half.abs() < 2e-4, "at K = 0.5 the volume change should vanish: {at_half:e}");
    for k in [0.33, 0.40, 0.44] {
        let change = make(k).volume_change(200).unwrap();
        println!("K = {k}: volume change {:+.4} %", 100.0 * change);
        assert!(change.abs() > at_half.abs(), "K = {k} should not conserve volume");
        assert!(change.abs() < 0.01, "but the change should be small: {change:e}");
    }
    // The direction is fixed: a neutral axis below mid-thickness means the
    // blank is shorter than the mid-surface, so the formed part has more
    // material than the blank.
    assert!(make(0.33).volume_change(200).unwrap() > 0.0);
}

#[test]
fn impossible_parts_are_refused() {
    let t = 2e-3;
    // A flange shorter than its setbacks cannot be formed.
    let part = SheetPart::new(t, 0.05, vec![0.02, 0.003], vec![Bend::at_right_angle(3.0 * t)]).unwrap();
    let err = part.straight(1).unwrap_err();
    assert!(err.contains("shorter than its setbacks"), "{err}");
    // The shop minimum for that flange is R + T.
    assert!((part.minimum_flange(1).unwrap() - (3.0 * t + t)).abs() < 1e-18);

    // Flange and bend counts must match.
    assert!(SheetPart::new(t, 0.05, vec![0.02, 0.02], vec![]).is_err());
    // A bend the other way is fine: only the magnitude costs material.
    let right = SheetPart::new(t, 0.05, vec![0.02, 0.02], vec![Bend { angle: -0.5 * PI, inside_radius: t }]).unwrap();
    let left = SheetPart::new(t, 0.05, vec![0.02, 0.02], vec![Bend::at_right_angle(t)]).unwrap();
    assert!((right.flat_length().unwrap() - left.flat_length().unwrap()).abs() < 1e-18);

    // A bend angle of zero or pi is not a bend.
    let flat = SheetPart::new(t, 0.05, vec![0.02, 0.02], vec![Bend { angle: 0.0, inside_radius: t }]).unwrap();
    assert!(flat.numbers(0).is_err());
    let folded = SheetPart::new(t, 0.05, vec![0.02, 0.02], vec![Bend { angle: PI, inside_radius: t }]).unwrap();
    assert!(folded.numbers(0).is_err());
}

#[test]
fn blank_and_formed_part_are_both_solids() {
    let t = 1.2e-3;
    let part = SheetPart::new(t, 0.04, vec![0.02, 0.03, 0.02], vec![
        Bend { angle: 0.5 * PI, inside_radius: 2.0 * t },
        Bend { angle: -0.5 * PI, inside_radius: 2.0 * t },
    ])
    .unwrap();
    // A Z-section: one bend each way, so the two flanges run the same way
    // and the overall length is both of them.
    let (lo, hi) = part.formed_bbox(32).unwrap();
    println!("formed Z: {:.3} x {:.3} mm", (hi[0] - lo[0]) * 1e3, (hi[1] - lo[1]) * 1e3);
    assert!(((hi[1] - lo[1]) * 1e3 - 30.0).abs() < 1.0, "the web sets the height");
    assert!((hi[0] - lo[0]) * 1e3 > 35.0, "both flanges run the same way");

    let blank = part.blank().unwrap();
    let mesh = mesh_solid(&blank, t / 3.0).expect("mesh");
    let expect = part.blank_volume().unwrap();
    assert!((mesh.report().volume - expect).abs() / expect < 0.02);

    let formed = mesh_solid(&part.solid(32).unwrap(), t / 3.0).expect("mesh");
    assert!(formed.report().watertight);
    // The two differ only by the model's volume change.
    let rel = (formed.report().volume - expect).abs() / expect;
    assert!(rel < 0.1, "formed {:e} vs blank {expect:e}", formed.report().volume);
}

#[test]
fn spring_stiffness_and_wire_length_are_the_published_forms() {
    // A music-wire compression spring: d = 1.6 mm, D = 12 mm, 8 coils,
    // free length 40 mm. G = 79.3 GPa for music wire (Shigley & Mischke,
    // table 10-5).
    let s = Spring::new(1.6e-3, 12e-3, 40e-3, 8.0);
    let g = 79.3e9;
    let k = s.stiffness(g);
    let expect = g * 1.6e-3f64.powi(4) / (8.0 * 12e-3f64.powi(3) * 8.0);
    assert!((k - expect).abs() / k < 1e-15);
    println!("spring: k = {:.1} N/m, index {:.2}, Wahl {:.4}", k, s.index(), s.wahl_factor());
    // 4.70 N/mm for this geometry.
    assert!((k - 4699.3).abs() < 1.0, "{k:.1} N/m");

    // The spring index, and the Wahl factor it implies.
    assert!((s.index() - 7.5).abs() < 1e-15);
    let c = s.index();
    assert!((s.wahl_factor() - ((4.0 * c - 1.0) / (4.0 * c - 4.0) + 0.615 / c)).abs() < 1e-18);
    assert!(s.wahl_factor() > 1.1 && s.wahl_factor() < 1.3, "{:.4}", s.wahl_factor());
    // Stress rises with the Wahl factor above the plain torsion value.
    let plain = 8.0 * 50.0 * 12e-3 / (PI * 1.6e-3f64.powi(3));
    assert!((s.shear_stress(50.0) / plain - s.wahl_factor()).abs() < 1e-12);

    // The wire length: an exact helix, n sqrt((pi D)^2 + p^2).
    let p = s.pitch();
    let exact = 8.0 * ((PI * 12e-3f64).powi(2) + p * p).sqrt();
    assert!((s.wire_length() - exact).abs() / exact < 1e-15);
    // And the swept solid's volume matches length times section.
    let (solid, check) = s.solid().expect("solid");
    assert!(check.watertight, "{check:?}");
    let rel = (check.volume - s.volume()).abs() / s.volume();
    println!("spring volume: swept {:e}, length x section {:e} ({rel:.3e})", check.volume, s.volume());
    // The swept section is a 16-gon, which holds 0.9936 of the circle.
    let inscribed = 16.0 / 2.0 * (2.0 * PI / 16.0).sin() / PI;
    assert!(rel < 1.0 - inscribed + 0.01, "swept volume off by {rel:.4}");
    assert!(solid.sdf([6e-3, 0.0, 0.0]).abs() < 1e-3);

    // Coils that clash are refused.
    let tight = Spring::new(4e-3, 12e-3, 20e-3, 8.0);
    assert!(tight.coil_clearance() < 0.0);
    let err = tight.solid().unwrap_err();
    assert!(err.contains("clash"), "{err}");

    // Drawable is not the same as makeable. A spring of index 1.6 sweeps
    // perfectly well — the wire clears the coil axis — and could not be
    // wound, so the geometry succeeds and the manufacturability check
    // fails.
    let fat = Spring::new(5e-3, 8e-3, 60e-3, 6.0);
    assert!((fat.index() - 1.6).abs() < 1e-12, "spring index {:.2}", fat.index());
    assert!(fat.solid().is_ok(), "the geometry is self-consistent");
    let err = fat.manufacturability().unwrap_err();
    assert!(err.contains("below 4"), "{err}");
    // And a very slender one tangles.
    let slender = Spring::new(0.5e-3, 10e-3, 40e-3, 10.0);
    assert!(slender.manufacturability().unwrap_err().contains("above 15"));
    s.manufacturability().expect("the music-wire spring is windable");

    // Solid length: the coils stacked.
    assert!((s.solid_length() - 9.0 * 1.6e-3).abs() < 1e-18);
    assert!(s.solid_length() < s.free_length);
}
