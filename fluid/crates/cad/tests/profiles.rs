//! Structural sections against the textbook composite formulas and a beam
//! solve, and NACA sections against their published closed forms.

use std::f64::consts::PI;

use cad::kits::airfoil::{Blade, Naca4, A4_CLOSED, A4_OPEN};
use cad::kits::section::{cantilever_tip_deflection, euler_buckling, Section};
use cad::meshing::mesh_solid;

#[test]
fn section_moments_match_the_composite_formulas() {
    // Green's theorem on the drawn polygon against the composite-rectangle
    // formula: two independent computation paths over the same shape.
    let (d, bf, tw, tf): (f64, f64, f64, f64) = (0.3, 0.15, 0.008, 0.012);
    let p = Section::IBeam { d, bf, tw, tf }.properties().unwrap();

    let hw = d - 2.0 * tf;
    let area = 2.0 * bf * tf + hw * tw;
    // 1e-14, not 1e-15: the shoelace sum over twelve vertices of an
    // I-beam accumulates a couple of units in the last place.
    assert!((p.area - area).abs() / area < 1e-14, "area {:e} vs {area:e}", p.area);

    // Ixx: two flanges about their own centroids plus their offsets, and
    // the web.
    let ixx = 2.0 * (bf * tf.powi(3) / 12.0 + bf * tf * (0.5 * (d - tf)).powi(2)) + tw * hw.powi(3) / 12.0;
    assert!((p.ixx - ixx).abs() / ixx < 1e-14, "Ixx {:e} vs {ixx:e}", p.ixx);
    // Iyy: the flanges bending the other way, plus the web.
    let iyy = 2.0 * (tf * bf.powi(3) / 12.0) + hw * tw.powi(3) / 12.0;
    assert!((p.iyy - iyy).abs() / iyy < 1e-14, "Iyy {:e} vs {iyy:e}", p.iyy);
    assert!(p.ixy.abs() < 1e-20, "a doubly symmetric shape has no product of inertia");
    assert!(p.centroid[0].abs() < 1e-15 && p.centroid[1].abs() < 1e-15);
    // Section modulus and radius of gyration follow.
    assert!((p.sx - ixx / (0.5 * d)).abs() / p.sx < 1e-14);
    assert!((p.rx - (ixx / area).sqrt()).abs() / p.rx < 1e-15);
    assert!(!p.exact_j, "an open section's torsion constant is an estimate");

    // A rectangle: bh^3/12 and hb^3/12, exactly.
    let (w, h): (f64, f64) = (0.05, 0.02);
    let r = Section::Rectangle { width: w, height: h }.properties().unwrap();
    assert!((r.ixx - w * h.powi(3) / 12.0).abs() / r.ixx < 1e-15);
    assert!((r.iyy - h * w.powi(3) / 12.0).abs() / r.iyy < 1e-15);
    assert!((r.area - w * h).abs() / r.area < 1e-15);
    assert!((r.j > 0.0) && !r.exact_j, "a rectangle's torsion constant is the series estimate");

    // A box section: the outside less the inside.
    let t: f64 = 0.004;
    let b = Section::Box { width: w, height: h, thickness: t }.properties().unwrap();
    let (iw, ih) = (w - 2.0 * t, h - 2.0 * t);
    let ixx = (w * h.powi(3) - iw * ih.powi(3)) / 12.0;
    assert!((b.ixx - ixx).abs() / ixx < 1e-14, "box Ixx {:e} vs {ixx:e}", b.ixx);
    assert!((b.area - (w * h - iw * ih)).abs() / b.area < 1e-14);
    assert!(b.exact_j, "Bredt's formula is a closed form");

    // Round and tubular shapes keep their closed forms rather than a
    // polygon's approximation.
    let dia: f64 = 0.04;
    let c = Section::Round { diameter: dia }.properties().unwrap();
    assert!((c.ixx - PI * (0.5 * dia).powi(4) / 4.0).abs() / c.ixx < 1e-15);
    assert!((c.area - PI * (0.5 * dia).powi(2)).abs() / c.area < 1e-15);
    assert!((c.j - 2.0 * c.ixx).abs() / c.j < 1e-15, "J = 2 I for a circle");
    let tube = Section::Tube { diameter: dia, thickness: 0.003 }.properties().unwrap();
    assert!(tube.ixx < c.ixx && tube.area < c.area);
}

#[test]
fn asymmetric_sections_find_their_own_centroid() {
    // A channel's centroid is off the web, and a tee's is up in the
    // flange: both come out of Green's theorem, not an assumption.
    let (d, bf, tw, tf): (f64, f64, f64, f64) = (0.2, 0.08, 0.006, 0.01);
    let ch = Section::Channel { d, bf, tw, tf }.properties().unwrap();
    let area = 2.0 * bf * tf + (d - 2.0 * tf) * tw;
    assert!((ch.area - area).abs() / area < 1e-14);
    // First moment about x = 0: two flanges at bf/2 plus the web at tw/2.
    let xbar = (2.0 * bf * tf * (0.5 * bf) + (d - 2.0 * tf) * tw * (0.5 * tw)) / area;
    assert!((ch.centroid[0] - xbar).abs() / xbar < 1e-14, "{:e} vs {xbar:e}", ch.centroid[0]);
    assert!(ch.centroid[1].abs() < 1e-15, "a channel is symmetric about x");

    let tee = Section::Tee { d: 0.15, bf: 0.1, tw: 0.008, tf: 0.012 }.properties().unwrap();
    let (web_a, fl_a) = ((0.15 - 0.012) * 0.008, 0.1 * 0.012);
    let ybar = (web_a * 0.5 * (0.15 - 0.012) + fl_a * (0.15 - 0.006)) / (web_a + fl_a);
    assert!((tee.centroid[1] - ybar).abs() / ybar < 1e-14, "{:e} vs {ybar:e}", tee.centroid[1]);
    // An angle has a product of inertia: it is not symmetric about either
    // centroidal axis.
    let ang = Section::Angle { leg_x: 0.06, leg_y: 0.06, thickness: 0.006 }.properties().unwrap();
    assert!(ang.ixy.abs() > 0.0, "an angle must have a product of inertia");
    assert!((ang.ixx - ang.iyy).abs() / ang.ixx < 1e-14, "an equal angle is symmetric about the diagonal");

    // Malformed shapes are refused.
    assert!(Section::IBeam { d: 0.1, bf: 0.05, tw: 0.004, tf: 0.06 }.properties().is_err());
    assert!(Section::Box { width: 0.02, height: 0.02, thickness: 0.011 }.properties().is_err());
    assert!(Section::Tube { diameter: 0.02, thickness: 0.02 }.properties().is_err());
}

#[test]
fn a_cad_section_in_the_beam_solver_matches_euler_bernoulli() {
    // The section properties drive `architected`'s frame elements, and the
    // tip deflection of a cantilever must be P L^3 / 3 E I exactly.
    use architected::frame::Frame;

    let section = Section::IBeam { d: 0.2, bf: 0.1, tw: 0.006, tf: 0.01 };
    let p = section.properties().unwrap();
    let (e, length, load) = (200e9, 3.0, 5000.0);

    let frame = Frame::straight_beam(length, 8, e, p.area, p.ixx);
    let mut frame = frame;
    frame.supports[0] = [true, true, true]; // fixed at the root
    let mut loads = vec![[0.0; 3]; frame.nodes.len()];
    let tip = frame.nodes.len() - 1;
    loads[tip] = [0.0, -load, 0.0];
    let (u, _) = frame.solve(&loads);

    let exact = cantilever_tip_deflection(load, length, e, p.ixx);
    let rel = (u[tip][1].abs() - exact).abs() / exact;
    println!("cantilever tip: {:e} m, exact {exact:e} ({rel:.3e})", u[tip][1].abs());
    assert!(rel < 1e-10, "tip deflection off by {rel:e}");

    // Euler buckling of the same section as a pinned column.
    let pcr = euler_buckling(e, p.iyy, length);
    assert!((pcr - PI * PI * e * p.iyy / (length * length)).abs() / pcr < 1e-15);
    // Weak-axis buckling governs.
    assert!(pcr < euler_buckling(e, p.ixx, length));

    // Plastic modulus exceeds the elastic one for the same shape.
    let zx = section.plastic_modulus_x().unwrap();
    assert!(zx > p.sx, "Zx {zx:e} should exceed Sx {:e}", p.sx);
    assert!(zx / p.sx < 1.3, "the shape factor of an I-beam is near 1.1: {:.3}", zx / p.sx);
}

#[test]
fn sections_extrude_into_meshable_beams() {
    for section in [
        Section::IBeam { d: 0.2, bf: 0.1, tw: 0.006, tf: 0.01 },
        Section::Channel { d: 0.15, bf: 0.06, tw: 0.005, tf: 0.008 },
        Section::Box { width: 0.08, height: 0.04, thickness: 0.004 },
        Section::Tube { diameter: 0.05, thickness: 0.004 },
        Section::Angle { leg_x: 0.05, leg_y: 0.05, thickness: 0.005 },
    ] {
        let beam = section.beam(0.4).expect("beam");
        let mesh = mesh_solid(&beam, 0.0015).expect("mesh");
        let report = mesh.report();
        assert!(report.watertight, "{section:?}: {report:?}");
        // Volume is area times length, to the meshing error.
        let expect = section.properties().unwrap().area * 0.4;
        let rel = (report.volume - expect).abs() / expect;
        assert!(rel < 0.03, "{section:?}: volume {:e} vs {expect:e} ({rel:.4})", report.volume);
    }
}

#[test]
fn naca_sections_match_their_published_closed_forms() {
    // NACA 2412: 2 % camber at 40 % chord, 12 % thick.
    let a = Naca4::designation(2412, 1.5).unwrap();
    assert!((a.camber - 0.02).abs() < 1e-15);
    assert!((a.camber_position - 0.4).abs() < 1e-15);
    assert!((a.thickness - 0.12).abs() < 1e-15);
    assert!(Naca4::designation(2400, 1.0).is_err(), "zero thickness is not a section");
    assert!(Naca4::designation(2012, 1.0).is_err(), "camber needs a position");

    for code in [12u32, 2412, 4415, 6309] {
        let f = Naca4::designation(code, 2.0).unwrap();
        // The maximum half-thickness is t/2 chords, at x = 0.2998.
        let xm = f.max_thickness_position();
        assert!((xm - 0.2998).abs() < 1e-3, "NACA {code}: thickest at {xm:.4}");
        // The polynomial's own maximum is 1.00027 t, not exactly t: the
        // series is nominally "t at 30 % chord", and the published
        // ordinates are, but the closed form overshoots by 0.027 %. That
        // is a property of the formula, so it is asserted rather than
        // tolerated.
        let tmax = 2.0 * f.half_thickness(xm);
        let overshoot = tmax / f.thickness;
        assert!(
            (1.0002..1.0004).contains(&overshoot),
            "NACA {code}: the polynomial's maximum is {overshoot:.6} t"
        );
        // The maximum camber is m, at p.
        if f.camber > 0.0 {
            let yc = f.camber_at(f.camber_position);
            assert!((yc - f.camber).abs() / f.camber < 1e-15, "NACA {code}: camber {yc:e}");
            // And the camber line is flat there.
            assert!(f.camber_slope(f.camber_position).abs() < 1e-15);
        }
        // Leading-edge radius: 1.1019 t^2 chords.
        let r = f.leading_edge_radius() / f.chord;
        assert!(
            (r - 1.1019 * f.thickness * f.thickness).abs() / r < 1e-3,
            "NACA {code}: LE radius {r:e} chords"
        );
        // The drawn polygon's area against the integral of the section as
        // constructed: thickness laid off perpendicular to the camber
        // line, so the area exceeds the closed-form 2 integral y_t dx.
        let drawn = f.polygon().area().abs();
        let constructed = f.area();
        let vertical = f.area_vertical();
        let rel = (drawn - constructed).abs() / constructed;
        assert!(
            rel < 2e-3,
            "NACA {code}: drawn area {drawn:e} vs constructed {constructed:e} ({rel:.3e})"
        );
        assert!(drawn < constructed, "a sampled outline is inscribed");
        if f.camber > 0.0 {
            assert!(
                constructed > vertical,
                "NACA {code}: the perpendicular construction must enclose more than the vertical one"
            );
        } else {
            // With no camber the two coincide exactly.
            assert!(
                (constructed - vertical).abs() / vertical < 1e-12,
                "NACA {code}: symmetric section areas differ by {:e}",
                (constructed - vertical).abs() / vertical
            );
        }
    }

    // The two trailing-edge conventions, and what each does.
    let open = Naca4 { closed: false, ..Naca4::symmetric(0.12, 1.0) };
    let closed = Naca4 { closed: true, ..Naca4::symmetric(0.12, 1.0) };
    // y_t(1) = 5t(a0 + a1 + a2 + a3 - a4), and the closed form's a4 is
    // chosen to make that zero, so the open form's gap is
    // 2 * 5t (a4_closed - a4_open).
    assert!((open.trailing_edge_gap() - 2.0 * 5.0 * 0.12 * (A4_CLOSED - A4_OPEN)).abs() < 1e-9);
    // The vertical-offset closed form for a 12 % symmetric section.
    assert!((open.area_vertical() - 0.685083 * 0.12).abs() < 1e-6, "{:e}", open.area_vertical());
    assert!((closed.area_vertical() - 0.680883 * 0.12).abs() < 1e-6, "{:e}", closed.area_vertical());
    assert!(open.trailing_edge_gap() > 1e-3, "the original series leaves a gap");
    assert!(closed.trailing_edge_gap() < 1e-15, "the closed form meets exactly");
    // Closing the trailing edge removes a sliver of area, exactly the a4
    // coefficient's own integral.
    assert!(closed.area_vertical() < open.area_vertical());
    let coefficient_difference = 2.0 * 5.0 * 0.12 * (A4_CLOSED - A4_OPEN) / 5.0;
    assert!(
        (open.area_vertical() - closed.area_vertical() - coefficient_difference).abs() < 1e-12,
        "the area difference is the a4 coefficient's integral"
    );
}

#[test]
fn naca_surfaces_are_perpendicular_to_the_camber_line() {
    // The standard construction offsets the thickness along the camber
    // normal, so the two surfaces are equidistant from the camber line.
    let f = Naca4::designation(4412, 1.0).unwrap();
    let mut worst: f64 = 0.0;
    for i in 1..200 {
        let x = i as f64 / 200.0;
        let (up, lo) = f.surfaces(x);
        let yt = f.half_thickness(x);
        let (cx, cy) = (x, f.camber_at(x));
        let du = ((up[0] - cx).powi(2) + (up[1] - cy).powi(2)).sqrt();
        let dl = ((lo[0] - cx).powi(2) + (lo[1] - cy).powi(2)).sqrt();
        worst = worst.max((du - yt).abs().max((dl - yt).abs()) / f.thickness);
    }
    assert!(worst < 1e-14, "surfaces are not at the half-thickness: {worst:e}");

    // A symmetric section's surfaces mirror exactly.
    let s = Naca4::symmetric(0.1, 1.0);
    for i in 0..100 {
        let x = i as f64 / 100.0;
        let (up, lo) = s.surfaces(x);
        assert!((up[0] - lo[0]).abs() < 1e-15 && (up[1] + lo[1]).abs() < 1e-15);
    }
}

#[test]
fn blades_loft_and_their_planform_is_exact() {
    let section = Naca4 { closed: true, points: 60, ..Naca4::designation(2412, 1.0).unwrap() };
    let blade = Blade::new(section, 2.0, 0.4, 0.2)
        .twisted(0.0, -4.0_f64.to_radians())
        .swept(10.0_f64.to_radians());

    // Trapezoidal planform, exactly.
    assert!((blade.planform_area() - 0.5 * 2.0 * 0.6).abs() < 1e-15);
    assert!((blade.taper() - 0.5).abs() < 1e-15);
    assert!((blade.aspect_ratio() - 4.0 / 0.6).abs() < 1e-14);
    // Mean aerodynamic chord of a taper-0.5 wing: (2/3) c_r (1+l+l^2)/(1+l).
    let mac = blade.mean_aerodynamic_chord();
    assert!((mac - 2.0 / 3.0 * 0.4 * 1.75 / 1.5).abs() < 1e-15, "MAC {mac:e}");
    assert!(mac < blade.root_chord && mac > blade.tip_chord);

    // The loft is a closed surface the kernel can take.
    let mesh = blade.mesh();
    let report = mesh.report();
    assert_eq!(report.boundary_edges, 0, "the blade surface must be closed: {report:?}");
    assert!(report.volume > 0.0);
    // Its volume is near the sum of the section areas times the span: a
    // straight loft of a 12 % section.
    let mean_area = 0.5 * section.area() * ((0.4 / 1.0f64).powi(2) + (0.2 / 1.0f64).powi(2));
    let rel = (report.volume - mean_area * 2.0).abs() / (mean_area * 2.0);
    assert!(rel < 0.1, "blade volume {:e} vs the trapezoid estimate ({rel:.3})", report.volume);

    let solid = blade.solid();
    assert!(solid.sdf([0.0, 0.0, 1.0]).is_finite());
}
