//! The cross-section solver against exact capacitances, and PCB impedance
//! against the published microstrip formulas.

use cad::laplace::{quasi_tem, Edge, Section, EPS0};
use cad::sketch::Region2;

#[test]
fn parallel_plate_capacitance_is_exact() {
    // With mirror side walls the field is exactly linear, so the discrete
    // capacitance must equal eps0 eps_r W / d to round-off.
    let (w, d, eps_r) = (3e-3, 1e-3, 4.4);
    for cell in [d / 8.0, d / 20.0] {
        let section = Section::new(([0.0, 0.0], [w, d]), cell)
            .with_dielectric(Region2::rect([0.5 * w, 0.5 * d], [w, 0.5 * d]), eps_r)
            .with_conductor(Region2::rect([0.5 * w, 0.0], [w, 1e-15]), 0.0)
            .with_conductor(Region2::rect([0.5 * w, d], [w, 1e-15]), 1.0)
            .with_edges([Edge::Mirror, Edge::Mirror, Edge::Grounded, Edge::Grounded]);
        let sol = section.solve(1e-12, 20000);
        let c = sol.capacitance(1.0);
        let exact = EPS0 * eps_r * w / d;
        let rel = (c - exact).abs() / exact;
        assert!(
            rel < 1e-9,
            "cell {cell:e}: C = {c:e} F/m vs exact {exact:e} ({rel:e}), residual {:e}",
            sol.residual
        );
    }
}

#[test]
fn two_dielectric_layers_in_series_are_exact() {
    // Stacked dielectrics between plates: 1/C = d1/(eps1 W) + d2/(eps2 W).
    let (w, d1, d2, e1, e2) = (2e-3, 0.4e-3, 0.6e-3, 4.4, 1.0);
    let d = d1 + d2;
    let section = Section::new(([0.0, 0.0], [w, d]), d / 25.0)
        .with_dielectric(Region2::rect([0.5 * w, 0.5 * d1], [w, 0.5 * d1]), e1)
        .with_dielectric(Region2::rect([0.5 * w, d1 + 0.5 * d2], [w, 0.5 * d2]), e2)
        .with_conductor(Region2::rect([0.5 * w, 0.0], [w, 1e-15]), 0.0)
        .with_conductor(Region2::rect([0.5 * w, d], [w, 1e-15]), 1.0)
        .with_edges([Edge::Mirror, Edge::Mirror, Edge::Grounded, Edge::Grounded]);
    let c = section.solve(1e-12, 20000).capacitance(1.0);
    let exact = 1.0 / (d1 / (EPS0 * e1 * w) + d2 / (EPS0 * e2 * w));
    let rel = (c - exact).abs() / exact;
    assert!(rel < 1e-9, "series C = {c:e} vs exact {exact:e} ({rel:e})");
}

#[test]
fn microstrip_impedance_matches_hammerstad_jensen() {
    use cad::kits::pcb::microstrip_z0;
    // A 1.6 mm FR-4 board, 35 um copper, for a range of trace widths.
    let (h, t, eps_r) = (1.6e-3, 35e-6, 4.4);
    let mut worst: f64 = 0.0;
    for w in [0.8e-3, 1.6e-3, 2.9e-3] {
        let (z_formula, _) = microstrip_z0(w, h, eps_r);
        // The same geometry, solved. Half the problem by symmetry.
        let span = 6.0 * (w + 2.0 * h);
        let height = h + 8.0 * h;
        let cell = h / 16.0;
        let section = Section::new(([0.0, 0.0], [0.5 * span, height]), cell)
            .with_dielectric(Region2::rect([0.25 * span, 0.5 * h], [0.5 * span, 0.5 * h]), eps_r)
            // Ground plane along the bottom.
            .with_conductor(Region2::rect([0.25 * span, 0.0], [0.5 * span, 1e-15]), 0.0)
            // Half the trace, against the mirror at x = 0.
            .with_conductor(Region2::rect([0.0, h + 0.5 * t], [0.5 * w, 0.5 * t]), 1.0)
            .with_edges([Edge::Mirror, Edge::Grounded, Edge::Grounded, Edge::Grounded])
            // The mirror cuts through the trace, so this is half the line.
            .with_multiplicity(2.0);
        let (z_field, eps_eff, sol) = quasi_tem(&section, 1.0, 1e-10, 40000);
        let rel = (z_field - z_formula).abs() / z_formula;
        println!(
            "w/h = {:.3}: Hammerstad-Jensen {z_formula:.2} ohm, field solve {z_field:.2} ohm \
             ({:+.2} %), eps_eff {eps_eff:.3}, {} iterations",
            w / h,
            100.0 * (z_field - z_formula) / z_formula,
            sol.iterations
        );
        worst = worst.max(rel);
    }
    // Hammerstad & Jensen quote about 1 % for the formula itself and
    // assume zero conductor thickness; the field solve carries the real
    // 35 um copper, which raises the capacitance and so lowers Z0, and
    // adds its own discretisation error at the trace corners. The sign of
    // the disagreement is the same for every width, as the thickness
    // correction requires.
    assert!(worst < 0.05, "worst disagreement {:.2} %", 100.0 * worst);

    // The 50 ohm rule of thumb for this stack: about 2.9 mm wide.
    let (z50, _) = microstrip_z0(2.9e-3, h, eps_r);
    assert!((z50 - 50.0).abs() < 2.0, "2.9 mm on 1.6 mm FR-4 gives {z50:.2} ohm");
}

#[test]
fn field_solve_converges_with_resolution() {
    use cad::kits::pcb::microstrip_z0;
    let (h, t, eps_r, w) = (1.6e-3, 35e-6, 4.4, 1.6e-3);
    let (z_formula, _) = microstrip_z0(w, h, eps_r);
    let mut zs = Vec::new();
    for div in [8.0, 16.0, 32.0] {
        let span = 6.0 * (w + 2.0 * h);
        let section = Section::new(([0.0, 0.0], [0.5 * span, 9.0 * h]), h / div)
            .with_dielectric(Region2::rect([0.25 * span, 0.5 * h], [0.5 * span, 0.5 * h]), eps_r)
            .with_conductor(Region2::rect([0.25 * span, 0.0], [0.5 * span, 1e-15]), 0.0)
            .with_conductor(Region2::rect([0.0, h + 0.5 * t], [0.5 * w, 0.5 * t]), 1.0)
            .with_edges([Edge::Mirror, Edge::Grounded, Edge::Grounded, Edge::Grounded])
            .with_multiplicity(2.0);
        let (z, _, _) = quasi_tem(&section, 1.0, 1e-10, 40000);
        zs.push(z);
    }
    println!("microstrip Z0 vs h/cell: {zs:?} (formula {z_formula:.2})");
    // Successive refinements must move less than the previous step.
    let (d1, d2) = ((zs[1] - zs[0]).abs(), (zs[2] - zs[1]).abs());
    assert!(d2 < d1, "not converging: steps {d1:.3} then {d2:.3}");
}

#[test]
fn conductor_edges_on_grid_lines_are_not_lost_to_round_off() {
    // A trace sitting exactly on the dielectric: its bottom row of nodes
    // lies on the region's edge, where the field can evaluate to +1e-20.
    // Classification must not depend on that, so growing the region by
    // 1 pm (far below a cell) must change nothing. Before the fix the
    // bottom row went to air and Z0 read 64.1 ohm instead of 49.2 ohm.
    let (w, h, t, er) = (0.38e-3, 0.2e-3, 70e-6, 4.2);
    let solve = |grow: f64| {
        let span = 6.0 * (w + 2.0 * h);
        let section = Section::new(([0.0, 0.0], [0.5 * span, 9.0 * h]), h / 8.0)
            .with_dielectric(Region2::rect([0.25 * span, 0.5 * h], [0.5 * span, 0.5 * h]), er)
            .with_conductor(Region2::rect([0.25 * span, 0.0], [0.5 * span, 1e-15]), 0.0)
            .with_conductor(Region2::rect([0.0, h + 0.5 * t], [0.5 * w + grow, 0.5 * t + grow]), 1.0)
            .with_edges([Edge::Mirror, Edge::Grounded, Edge::Grounded, Edge::Grounded])
            .with_multiplicity(2.0);
        quasi_tem(&section, 1.0, 1e-10, 40000).0
    };
    let (exact, grown) = (solve(0.0), solve(1e-12));
    assert!((exact - grown).abs() < 1e-9 * exact, "{exact} vs {grown}");
}
