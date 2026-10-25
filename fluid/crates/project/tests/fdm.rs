//! FDM printing against exact geometry.

use cad::Solid;
use project::fdm::analyse;

#[test]
fn printed_volume_of_a_box() {
    // 40 × 30 × 20 mm, 2 perimeters of 0.4 mm, 20 % infill: the walls are
    // V − (40 − 1.6)(30 − 1.6)(20 − 1.6) mm³ exactly.
    let mm = 1e-3;
    let b = Solid::box_between([0.0; 3], [40.0 * mm, 30.0 * mm, 20.0 * mm]);
    let cell = 0.25 * mm;
    let p = analyse(&b, 40.0 * 30.0 * 20.0 * mm * mm * mm, cell, 0.8 * mm, 0.2).unwrap();
    let inner = 38.4 * 28.4 * 18.4 * mm * mm * mm;
    let exact = (40.0 * 30.0 * 20.0 * mm * mm * mm - inner) + 0.2 * inner;
    println!("printed {:.4} cm³ ± {:.4} vs exact {:.4} cm³", p.printed_volume * 1e6, p.printed_error * 1e6, exact * 1e6);
    assert!((p.printed_volume - exact).abs() <= p.printed_error.max(1e-9 * exact), "{p:?} vs {exact}");
    // A plate is 20 mm at its thinnest here (the box's smallest side).
    assert!((p.min_wall - 20.0 * mm).abs() <= p.min_wall_error, "{p:?}");
}

#[test]
fn thinnest_wall_is_found() {
    let mm = 1e-3;
    // A 1.2 mm plate standing on a 5 mm block: the plate is the thinnest wall.
    let s = Solid::box_between([0.0; 3], [20.0 * mm, 20.0 * mm, 5.0 * mm])
        .plus(Solid::box_between([9.4 * mm, 0.0, 5.0 * mm], [10.6 * mm, 20.0 * mm, 15.0 * mm]));
    let cell = 0.2 * mm;
    let p = analyse(&s, 1.0, cell, 0.8 * mm, 0.2).unwrap();
    println!("plate: thinnest {:.3} ± {:.3} mm (exact 1.2)", p.min_wall / mm, p.min_wall_error / mm);
    assert!((p.min_wall - 1.2 * mm).abs() <= p.min_wall_error, "{p:?}");
    // A tube with a 1 mm wall, curved: same claim.
    let tube = Solid::cylinder_z([0.0; 3], 6.0 * mm, 10.0 * mm).minus(Solid::cylinder_z([0.0; 3], 5.0 * mm, 12.0 * mm));
    let p = analyse(&tube, 1.0, cell, 0.8 * mm, 0.2).unwrap();
    println!("tube: thinnest {:.3} ± {:.3} mm (exact 1.0)", p.min_wall / mm, p.min_wall_error / mm);
    assert!((p.min_wall - 1.0 * mm).abs() <= p.min_wall_error, "{p:?}");
}

#[test]
fn the_printed_duct_in_the_example() {
    // Walls of 2 mm with 3 × 0.4 mm perimeters print solid, so the printed
    // mass is the full PLA mass; the least time is volume / flow.
    use project::Project;
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../projects/motor-driver/project.json");
    let p = Project::load(&path).unwrap();
    let e = p.evaluate();
    let v = |id: &str| e.value(id).unwrap();
    assert!((v("duct_body.printed_mass").nom / v("duct_body.mass").nom - 1.0).abs() < 1e-3);
    let t = v("duct_body.print_time_min").nom;
    // PLA 1240 kg/m³ (the library's value), 12 mm³/s.
    let printed_volume = v("duct_body.printed_mass").nom / 1240.0;
    assert!((t - printed_volume / 12e-9).abs() < 1e-9 * t, "seconds, from SI volume and flow: {t}");
    assert!(t > 25.0 * 60.0 && t < 35.0 * 60.0);
    // Exact outer section 44 × 27 mm minus 40 × 23 mm, 80 mm long.
    let exact = (44.0 * 27.0 - 40.0 * 23.0) * 80.0 * 1e-9;
    assert!((v("duct_body.volume").nom / exact - 1.0).abs() < 2e-3);
    assert!((v("duct_body.min_wall").nom - 2e-3).abs() <= 0.75e-3);
}
