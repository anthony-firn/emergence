//! The microstrip analysis: its error claim against a finer, larger solve,
//! its rung choice, and its declared monotonicity.

use cad::laplace::{quasi_tem, Edge, Section};
use cad::Region2;
use project::analysis::{Analysis, Microstrip};

/// An independent reference: finer than either of the analysis's grids and
/// in a box 1.5 times as large as its standard one.
fn reference(w: f64, h: f64, t: f64, er: f64, cell: f64) -> f64 {
    let span = 9.0 * (w + 2.0 * h);
    let section = Section::new(([0.0, 0.0], [0.5 * span, 13.0 * h]), cell)
        .with_dielectric(Region2::rect([0.25 * span, 0.5 * h], [0.5 * span, 0.5 * h]), er)
        .with_conductor(Region2::rect([0.25 * span, 0.0], [0.5 * span, 1e-15]), 0.0)
        .with_conductor(Region2::rect([0.0, h + 0.5 * t], [0.5 * w, 0.5 * t]), 1.0)
        .with_edges([Edge::Mirror, Edge::Grounded, Edge::Grounded, Edge::Grounded])
        .with_multiplicity(2.0);
    quasi_tem(&section, 1.0, 1e-9, 60000).0
}

#[test]
fn the_error_claim_covers_a_finer_larger_solve() {
    // Geometry on grid lines at h/12, h/24 and h/36: 50 um copper, 0.3 mm
    // trace on 0.2 mm laminate.
    let (w, h, t, er) = (0.3e-3, 0.2e-3, 50e-6, 4.2);
    let out = Microstrip.run(&[w, t, h, er]).unwrap();
    let (z, err) = (out.values[0], out.errors[0]);
    let z_ref = reference(w, h, t, er, h / 36.0);
    println!("analysis {z:.3} ± {err:.3} Ω, reference {z_ref:.3} Ω (miss {:.3}); {:?}", (z - z_ref).abs(), out.notes);
    assert!((z - z_ref).abs() <= err, "{z} ± {err} misses {z_ref}");
    // The zero-thickness formula is the leading-order check: it must sit
    // above the thick-copper answer (copper adds capacitance).
    assert!(out.values[3] > z);
}

#[test]
fn thin_copper_takes_the_closed_form_rung() {
    // 35 um on 1.6 mm FR-4: t/h = 0.022, too thin for the uniform grid.
    let out = Microstrip.run(&[2.9e-3, 35e-6, 1.6e-3, 4.4]).unwrap();
    assert!(out.notes[0].contains("closed form"), "{:?}", out.notes);
    // The 50 ohm rule of thumb for this stack, within the stated error.
    assert!((out.values[0] - 50.0).abs() < out.errors[0] + 2.0);
    assert!(Microstrip.run(&[0.01e-3, 35e-6, 1.6e-3, 4.4]).is_err(), "W/h outside the validated range must be refused");
}

#[test]
fn impedance_is_monotone_where_declared() {
    // Dirichlet's principle: wider, thicker copper and higher permittivity
    // can only lower Z0.
    let base = [0.3e-3, 50e-6, 0.2e-3, 4.2];
    let z = |x: [f64; 4]| Microstrip.run(&x).unwrap().values[0];
    let z0 = z(base);
    for (i, step) in [(0, 0.05e-3), (1, 10e-6), (3, 0.3)] {
        let mut x = base;
        x[i] += step;
        assert!(z(x) < z0, "input {i}: {} vs {z0}", z(x));
    }
}
