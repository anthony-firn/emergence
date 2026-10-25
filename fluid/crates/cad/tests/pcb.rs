//! The PCB kit: stackup, copper geometry, electrical answers about what
//! was drawn, and manufacturing output.

use std::f64::consts::PI;

use cad::kits::pcb::{
    excellon, flatten, gerber_layer, microstrip_width_for, microstrip_z0, Board, Feature, PcbLayer, Via, OZ_COPPER,
    RHO_CU,
};
use cad::sketch::{Polygon, Region2};

fn two_layer(board_mm: f64) -> Board {
    let half = 0.5 * board_mm * 1e-3;
    Board::new(
        Region2::rounded_rect([0.0, 0.0], [half, half], 1e-3),
        vec![
            PcbLayer::Copper { ounces: 1.0 },
            PcbLayer::Dielectric {
                material: "FR-4".into(),
                thickness: 1.5e-3,
                eps_r: 4.4,
                loss_tangent: 0.02,
                kappa: 0.29,
            },
            PcbLayer::Copper { ounces: 1.0 },
        ],
    )
}

#[test]
fn stackup_thickness_and_sheet_resistance() {
    let board = two_layer(10.0);
    let t = board.thickness();
    assert!((t - (1.5e-3 + 2.0 * OZ_COPPER)).abs() < 1e-15, "thickness {t:e}");
    assert_eq!(board.copper_layers(), vec![0, 2]);

    // 1 oz copper is 34.8 um, so the sheet resistance is
    // 1.678e-8 / 34.8e-6 = 0.482 mohm per square, the familiar
    // "half a milliohm per square" of a 1 oz layer.
    let rs = board.sheet_resistance(0).unwrap();
    assert!((rs - RHO_CU / OZ_COPPER).abs() < 1e-18);
    assert!((rs * 1e3 - 0.4822).abs() < 1e-3, "{:.4} mohm/sq", rs * 1e3);

    // The dielectric between the copper layers is what a microstrip
    // formula needs.
    let (h, eps) = board.dielectric_below(1).expect("dielectric");
    assert!((h - 1.5e-3).abs() < 1e-15 && (eps - 4.4).abs() < 1e-15);

    // Series thermal resistance of the stack, m^2 K/W.
    let r = board.thermal_resistance().unwrap();
    let expect = 1.5e-3 / 0.29 + 2.0 * OZ_COPPER / 401.0;
    assert!((r - expect).abs() / expect < 1e-12, "{r:e} vs {expect:e}");

    // And the stack hands straight to the optics/thermal crates.
    let stack = board.layer_stack();
    assert_eq!(stack.layers.len(), 3);
    assert!((stack.total_thickness() - t).abs() < 1e-15);
}

#[test]
fn trace_resistance_from_the_drawn_geometry() {
    // A 1 mm wide, 100 mm long trace in 1 oz copper: R = rho L / (w t).
    let board = two_layer(120.0).with(Feature::Trace {
        layer: 0,
        points: vec![[-50e-3, 0.0], [50e-3, 0.0]],
        width: 1e-3,
    });
    let r20 = board.trace_resistance(0, 20.0).unwrap();
    let expect = RHO_CU * 100e-3 / (1e-3 * OZ_COPPER);
    assert!((r20 - expect).abs() / expect < 1e-12, "{r20:e} vs {expect:e}");
    assert!((r20 * 1e3 - 48.2).abs() < 0.2, "{:.2} mohm", r20 * 1e3);

    // Copper gains 0.393 % per kelvin.
    let r85 = board.trace_resistance(0, 85.0).unwrap();
    let ratio = r85 / r20;
    assert!((ratio - (1.0 + 3.93e-3 * 65.0)).abs() < 1e-12, "{ratio}");

    // A bent trace's length is the sum of its segments.
    let bent = two_layer(120.0).with(Feature::Trace {
        layer: 0,
        points: vec![[0.0, 0.0], [30e-3, 0.0], [30e-3, 40e-3]],
        width: 0.5e-3,
    });
    let r = bent.trace_resistance(0, 20.0).unwrap();
    let expect = RHO_CU * 70e-3 / (0.5e-3 * OZ_COPPER);
    assert!((r - expect).abs() / expect < 1e-12);
}

#[test]
fn via_resistance_and_drill_geometry() {
    let board = two_layer(10.0).with_via(Via {
        center: [0.0, 0.0],
        drill: 0.3e-3,
        annular_ring: 0.15e-3,
        plating: 25e-6,
    });
    let r = board.via_resistance(0).unwrap();
    let r_out = 0.15e-3 + 25e-6;
    let area = PI * (r_out * r_out - 0.15e-3 * 0.15e-3);
    let expect = RHO_CU * board.thickness() / area;
    assert!((r - expect).abs() / expect < 1e-12);
    // A 0.3 mm via through 1.57 mm of board is about a milliohm.
    assert!(r > 0.4e-3 && r < 2e-3, "{:e} ohm", r);

    // The drill is open through the board and the ring is copper.
    let solid = board.solid();
    assert!(solid.sdf([0.0, 0.0, 0.5e-3]) > 0.0, "the drill must go through");
    assert!(solid.sdf([0.0, 0.25e-3, 0.5 * OZ_COPPER]) < 0.0, "the annular ring");
}

#[test]
fn copper_coverage_matches_the_drawn_area() {
    let board = two_layer(10.0).with(Feature::Trace {
        layer: 0,
        points: vec![[-4e-3, 0.0], [4e-3, 0.0]],
        width: 1e-3,
    });
    // A capsule: a rectangle plus two half-discs.
    let copper_area = 8e-3 * 1e-3 + PI * 0.5e-3 * 0.5e-3;
    let board_area = board.outline.area_quadrature(10e-6);
    let expect = copper_area / board_area;
    let coverage = board.coverage(0, 10e-6);
    assert!(
        (coverage - expect).abs() / expect < 0.01,
        "coverage {coverage:.5} vs {expect:.5}"
    );
    // Nothing was drawn on the bottom layer.
    assert!(board.coverage(1, 20e-6) < 1e-12);
}

#[test]
fn board_solid_has_copper_where_it_was_drawn() {
    let board = two_layer(10.0)
        .with(Feature::Trace {
            layer: 1,
            points: vec![[-3e-3, 0.0], [3e-3, 0.0]],
            width: 0.6e-3,
        })
        .with(Feature::Pour {
            layer: 0,
            region: Region2::rect([0.0, 0.0], [4e-3, 4e-3]),
        });
    let solid = board.solid();
    let t = board.thickness();
    let z_top = t - 0.5 * OZ_COPPER;
    // Top-layer trace: copper on the centreline, none 1 mm away.
    assert!(solid.sdf([0.0, 0.0, z_top]) < 0.0, "trace copper missing");
    assert!(solid.sdf([0.0, 1e-3, z_top]) > 0.0, "copper where none was drawn");
    // Bottom pour.
    assert!(solid.sdf([0.0, 1e-3, 0.5 * OZ_COPPER]) < 0.0, "pour missing");
    assert!(solid.sdf([4.5e-3, 0.0, 0.5 * OZ_COPPER]) > 0.0, "pour too big");
    // The dielectric fills the middle everywhere inside the outline.
    assert!(solid.sdf([4.5e-3, 0.0, 0.5 * t]) < 0.0, "dielectric missing");
    assert!(solid.sdf([0.0, 0.0, t + 1e-4]) > 0.0, "material above the board");

    // A small window of the board meshes watertight.
    let window = solid.clipped([-2e-3, -2e-3, -1e-4], [2e-3, 2e-3, t + 1e-4]);
    let mesh = cad::meshing::mesh_solid(&window, 15e-6).expect("mesh");
    assert!(mesh.report().watertight, "{:?}", mesh.report());
}

#[test]
fn impedance_synthesis_round_trips() {
    let (h, eps_r) = (1.5e-3, 4.4);
    for target in [50.0, 75.0, 90.0] {
        let w = microstrip_width_for(target, h, eps_r);
        let (z, _) = microstrip_z0(w, h, eps_r);
        assert!((z - target).abs() < 1e-6, "{target} ohm -> {w:e} m -> {z} ohm");
    }
    // 50 ohms on 1.5 mm FR-4 wants roughly 2.7 mm of copper.
    let w50 = microstrip_width_for(50.0, h, eps_r);
    assert!(w50 > 2.4e-3 && w50 < 3.0e-3, "{:.3} mm", w50 * 1e3);
}

#[test]
fn boundary_tracing_recovers_area_and_perimeter() {
    // flatten() is what the Gerber writer emits, so its accuracy is the
    // accuracy of the manufacturing output.
    let r = 2e-3;
    let circle = Region2::circle([0.0, 0.0], r);
    let mut last_err = f64::INFINITY;
    for cell in [r / 10.0, r / 20.0, r / 40.0] {
        let loops = flatten(&circle, cell);
        assert_eq!(loops.len(), 1, "a disc has one boundary loop");
        let poly = Polygon::new(loops[0].clone());
        let area_err = (poly.area().abs() - PI * r * r).abs() / (PI * r * r);
        let perim_err = (poly.perimeter() - 2.0 * PI * r).abs() / (2.0 * PI * r);
        assert!(perim_err < 0.02, "perimeter error {perim_err:.4} at cell {cell:e}");
        assert!(area_err < last_err, "area error not falling: {area_err:e}");
        last_err = area_err;
    }
    assert!(last_err < 1e-3, "finest area error {last_err:e}");

    // An annulus has two loops.
    let ring = Region2::ring([0.0, 0.0], 1e-3, 2e-3);
    assert_eq!(flatten(&ring, 5e-5).len(), 2);
}

#[test]
fn manufacturing_files_are_structurally_sound() {
    let board = two_layer(10.0)
        .with(Feature::Trace {
            layer: 0,
            points: vec![[-3e-3, 0.0], [3e-3, 0.0]],
            width: 0.6e-3,
        })
        .with_via(Via { center: [2e-3, 1e-3], drill: 0.3e-3, annular_ring: 0.15e-3, plating: 25e-6 })
        .with_via(Via { center: [-2e-3, 1e-3], drill: 0.3e-3, annular_ring: 0.15e-3, plating: 25e-6 })
        .with_via(Via { center: [0.0, -2e-3], drill: 0.8e-3, annular_ring: 0.2e-3, plating: 25e-6 });

    let gerber = gerber_layer(&board.copper(0), 25e-6, "top copper");
    assert!(gerber.starts_with("G04 top copper"));
    assert!(gerber.contains("%FSLAX46Y46*%") && gerber.contains("%MOMM*%"));
    assert!(gerber.trim_end().ends_with("M02*"));
    // Every region opens and closes.
    assert_eq!(gerber.matches("G36*").count(), gerber.matches("G37*").count());
    assert!(gerber.matches("G36*").count() >= 1);
    // Coordinates are in 1e-6 mm: a 3 mm feature reads as 3000000.
    let max_coord: i64 = gerber
        .lines()
        .filter(|l| l.starts_with('X'))
        .filter_map(|l| l[1..].split('Y').next().and_then(|v| v.parse::<i64>().ok()))
        .map(i64::abs)
        .max()
        .unwrap();
    assert!(
        (3_300_000..6_000_000).contains(&max_coord),
        "largest X coordinate {max_coord} should be a few million units (mm * 1e6)"
    );

    let drills = excellon(&board);
    assert!(drills.starts_with("M48") && drills.trim_end().ends_with("M30"));
    assert!(drills.contains("METRIC,TZ"));
    // Two tool sizes for the two drill diameters, three hits in total.
    assert_eq!(drills.matches("T01C").count() + drills.matches("T02C").count(), 2);
    assert!(drills.contains("C0.300") && drills.contains("C0.800"));
    assert_eq!(drills.lines().filter(|l| l.starts_with('X')).count(), 3);
}
