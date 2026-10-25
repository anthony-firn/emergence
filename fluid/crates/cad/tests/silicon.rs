//! The silicon process kit: exact profile identities, the Deal-Grove
//! limits, and a device built through the CAD path against the validated
//! structure in `semiconductor`.

use std::f64::consts::PI;

use cad::kits::silicon::{
    arrhenius, deal_grove, deal_grove_offset, Depth, Dopant, Film, Implant, MosProcess, Process,
};

#[test]
fn implant_dose_is_conserved() {
    // The analytic integral of each profile, against quadrature of the
    // profile that was actually built.
    let cases: Vec<Depth> = vec![
        Depth::Implanted { dose: 1e19, range: 50e-9, straggle: 20e-9 },
        Depth::Surface { peak: 1e26, length: 80e-9 },
        Depth::Erfc { surface: 1e26, length: 60e-9 },
    ];
    for depth in cases {
        let h = 1e-11;
        let mut sum = 0.0;
        for i in 0..100_000 {
            let y = (i as f64 + 0.5) * h;
            sum += depth.at(y) * h;
        }
        let analytic = depth.dose();
        let rel = (sum - analytic).abs() / analytic;
        assert!(rel < 1e-6, "{depth:?}: quadrature {sum:e} vs analytic {analytic:e} ({rel:e})");
    }

    // A surface Gaussian's dose is N0 L sqrt(pi) / 2, exactly.
    let d = Depth::Surface { peak: 1e26, length: 80e-9 };
    assert!((d.dose() - 0.5 * 1e26 * 80e-9 * PI.sqrt()).abs() / d.dose() < 1e-15);
}

#[test]
fn anneal_broadens_a_gaussian_and_keeps_its_dose() {
    // Fick's equation: sigma^2 -> sigma^2 + 2 D t, with the dose fixed.
    let implanted = Depth::Implanted { dose: 1e19, range: 100e-9, straggle: 20e-9 };
    let dt = 1e-16; // m^2
    let after = implanted.annealed(dt);
    match after {
        Depth::Implanted { straggle, dose, .. } => {
            let expect = (20e-9f64.powi(2) + 2.0 * dt).sqrt();
            assert!((straggle - expect).abs() / expect < 1e-15, "{straggle:e} vs {expect:e}");
            assert!((dose - 1e19).abs() / 1e19 < 1e-15, "dose must not change");
        }
        _ => panic!("type changed"),
    }
    // Dose is conserved for the surface profile too, where the peak must
    // fall as the profile spreads.
    let surface = Depth::Surface { peak: 1e26, length: 80e-9 };
    let a = surface.annealed(2e-16);
    let rel = (a.dose() - surface.dose()).abs() / surface.dose();
    assert!(rel < 1e-15, "surface dose changed by {rel:e}");
    // And it really did spread.
    match (surface, a) {
        (Depth::Surface { length: l0, .. }, Depth::Surface { length: l1, peak: p1 }) => {
            assert!(l1 > l0 * 1.05 && p1 < 1e26);
        }
        _ => panic!("type changed"),
    }

    // Arrhenius is just the exponential; check one decade by hand.
    let d = arrhenius(1e-4, 3.5, 1273.15);
    let expect: f64 = 1e-4 * f64::exp(-3.5 * 1.602176634e-19 / (1.380649e-23 * 1273.15));
    assert!((d - expect).abs() / expect < 1e-15);
}

#[test]
fn deal_grove_satisfies_its_own_law_and_limits() {
    // Rate constants chosen to exercise both regimes, not taken from a
    // table: the law is what is being checked.
    let (a, b) = (0.165e-6, 0.0117e-12 / 3600.0 * 1e0);
    for time in [1.0, 60.0, 3600.0, 36000.0] {
        let x = deal_grove(a, b, 0.0, time);
        // x^2 + A x = B t, exactly.
        let residual = (x * x + a * x - b * time).abs() / (b * time);
        assert!(residual < 1e-12, "t = {time}: residual {residual:e}");
    }
    // Short times: linear, x -> (B/A) t.
    let t_short = 1e-3;
    let x = deal_grove(a, b, 0.0, t_short);
    let linear = b / a * t_short;
    assert!((x - linear).abs() / linear < 1e-6, "{x:e} vs {linear:e}");
    // Long times: parabolic. The two-term asymptotic is
    // x -> sqrt(B t) - A/2, and that is what the root must approach; the
    // leading term alone is only good to A/(2 sqrt(B t)).
    let t_long = 1e9;
    let x = deal_grove(a, b, 0.0, t_long);
    let parabolic = (b * t_long).sqrt();
    // The expansion is x = sqrt(Bt) - A/2 + A^2/(8 sqrt(Bt)) + ..., so the
    // residual after two terms must be the third term and no more.
    let third = a * a / (8.0 * parabolic);
    let residual = (x - (parabolic - 0.5 * a)).abs();
    assert!(
        residual < 1.01 * third && residual > 0.9 * third,
        "residual {residual:e} against the next asymptotic term {third:e}"
    );
    // An existing oxide is accounted for by tau: growing from x0 for zero
    // time leaves x0.
    let x0 = 20e-9;
    let tau = deal_grove_offset(a, b, x0);
    assert!((deal_grove(a, b, tau, 0.0) - x0).abs() / x0 < 1e-12);
}

#[test]
fn masked_films_and_junctions_land_where_drawn() {
    let proc = MosProcess::nmos(0.25e-6);
    let p = proc.process();
    let (xs, xd) = proc.edges();
    let mid = 0.5 * (xs + xd);

    // Oxide everywhere, poly only over the gate.
    assert_eq!(p.material_at(mid, -5e-9), Some("oxide"));
    assert_eq!(p.material_at(mid, -(proc.oxide + 50e-9)), Some("poly"));
    assert_eq!(p.material_at(0.1e-6, -(proc.oxide + 50e-9)), None, "poly outside the gate mask");
    assert_eq!(p.material_at(0.1e-6, -5e-9), Some("oxide"));
    assert_eq!(p.material_at(mid, 10e-9), Some("silicon"));
    assert!((p.surface_height(mid) - (proc.oxide + proc.gate_thickness)).abs() < 1e-18);
    assert!((p.surface_height(0.1e-6) - proc.oxide).abs() < 1e-18);

    // Under the source, the junction sits at the specified depth.
    let depth = p.junction_depth(0.2e-6).expect("junction under the source");
    assert!(
        (depth - proc.junction).abs() / proc.junction < 0.02,
        "junction at {depth:e}, specified {:e}",
        proc.junction
    );
    // Mid-channel there is no junction: the substrate is undisturbed.
    assert!(p.junction_depth(mid).is_none(), "the channel must not be inverted by doping");

    // The source/drain doping dominates at the surface outside the gate
    // and the substrate dominates in the channel.
    let (d_src, a_src) = p.doping_at(0.2e-6, 0.0);
    assert!(d_src > 100.0 * a_src, "source should be strongly n-type");
    let (d_ch, a_ch) = p.doping_at(mid, 0.0);
    assert!(a_ch > 10.0 * d_ch, "channel should stay p-type");

    // The blanket part of the stack is a layer stack; the gate is masked
    // and is reported as left out rather than silently dropped.
    let (stack, skipped) = p.layer_stack();
    assert_eq!(skipped, vec!["poly".to_string()]);
    assert_eq!(stack.layers.len(), 2);
    assert_eq!(stack.layers[0].material, "oxide");

    // The 3D solid has the gate stack standing above the field oxide.
    let solid = p.solid(1e-6);
    let z = 0.5e-6;
    assert!(solid.sdf([mid, -(proc.oxide + 50e-9), z]) < 0.0, "gate missing");
    assert!(solid.sdf([0.1e-6, -(proc.oxide + 50e-9), z]) > 0.0, "gate too wide");
    assert!(solid.sdf([0.1e-6, -5e-9, z]) < 0.0, "oxide missing");
    assert!(solid.sdf([mid, 0.5e-6, z]) < 0.0, "silicon body missing");
}

#[test]
fn process_reproduces_the_validated_mosfet_structure() {
    // An interface test, not a physics one: the point is that the process
    // flow can express the structure `semiconductor::mosfet` validates
    // against published drift-diffusion results, so a CAD-authored device
    // inherits that validation.
    use semiconductor::mosfet::MosfetSpec;

    let l_gate = 0.25e-6;
    let spec = MosfetSpec::nmos(l_gate);
    let mos = MosProcess {
        gate_length: l_gate,
        oxide: spec.tox,
        substrate: spec.n_sub,
        source_drain: spec.n_sd,
        lateral: spec.sigma_x,
        junction: spec.y_j,
        access: spec.l_sd,
        depth: spec.depth,
        pmos: spec.pmos,
        gate_thickness: 100e-9,
    };
    let p = mos.process();

    // The lateral junction and effective channel length agree with the
    // reference implementation.
    assert!(
        (mos.lateral_junction() - spec.lateral_junction()).abs() < 1e-15,
        "{:e} vs {:e}",
        mos.lateral_junction(),
        spec.lateral_junction()
    );
    assert!((mos.effective_length() - spec.l_eff()).abs() < 1e-15);

    // And the doping closures agree everywhere on the cross-section.
    let sigma_y = spec.y_j / (spec.n_sd / spec.n_sub).ln().sqrt();
    let (xs, xd) = spec.edges();
    let reference = |x: f64, y: f64| -> (f64, f64) {
        let vert = (-(y / sigma_y).powi(2)).exp();
        let lat = 0.5 * semiconductor::mosfet::erfc((x - xs) / spec.sigma_x)
            + 0.5 * semiconductor::mosfet::erfc((xd - x) / spec.sigma_x);
        let sd = spec.n_sd * vert * lat.min(1.0);
        (sd, spec.n_sub)
    };
    let mut worst: f64 = 0.0;
    for i in 0..=120 {
        let x = mos.total_length() * i as f64 / 120.0;
        for j in 0..=60 {
            let y = spec.depth * j as f64 / 60.0;
            let (d1, a1) = p.doping_at(x, y);
            let (d0, a0) = reference(x, y);
            let scale = spec.n_sd;
            worst = worst.max(((d1 - d0).abs() + (a1 - a0).abs()) / scale);
        }
    }
    assert!(worst < 1e-14, "doping closures differ by {worst:e} of the peak");
}

#[test]
fn cad_built_device_solves_like_the_reference() {
    // End to end: a drift-diffusion solve on a device whose doping comes
    // from the process flow, against the same device built by
    // `semiconductor::mosfet`. The currents must match to the solver's own
    // convergence tolerance.
    use semiconductor::device::{Contact, Device};
    use semiconductor::mesh::Material;
    use semiconductor::mosfet::MosfetSpec;
    use semiconductor::silicon::Silicon;

    let l_gate = 0.25e-6;
    // A coarser mesh than the default: the comparison between the two
    // construction paths is exact, so resolution only costs time.
    let spec = MosfetSpec { refine: 0.5, ..MosfetSpec::nmos(l_gate) };
    let si = Silicon::new(300.0);
    let mut reference = spec.build(si.clone());

    // Rebuild the same device, taking the structure from the CAD process.
    let mos = MosProcess {
        gate_length: l_gate,
        oxide: spec.tox,
        substrate: spec.n_sub,
        source_drain: spec.n_sd,
        lateral: spec.sigma_x,
        junction: spec.y_j,
        access: spec.l_sd,
        depth: spec.depth,
        pmos: false,
        gate_thickness: 100e-9,
    };
    let process = mos.process();
    let mesh = reference.mesh.clone();
    // The mesh carries its own materials; check the process agrees with
    // them before reusing it.
    for i in 0..mesh.nx() - 1 {
        for j in 0..mesh.ny() - 1 {
            let x = 0.5 * (mesh.x[i] + mesh.x[i + 1]);
            let y = 0.5 * (mesh.y[j] + mesh.y[j + 1]);
            let want = match mesh.elem[i * (mesh.ny() - 1) + j] {
                Material::Oxide => "oxide",
                Material::Silicon => "silicon",
            };
            let got = process.material_at(x, y).unwrap_or("none");
            // Above the gate the process has polysilicon where the
            // reference mesh models only oxide; everywhere else they must
            // agree.
            if got != want && !(want == "oxide" && got == "poly") {
                panic!("material at ({x:e}, {y:e}): process says {got}, mesh says {want}");
            }
        }
    }

    let contacts: Vec<Contact> = reference
        .contacts
        .iter()
        .map(|c| Contact {
            name: c.name.clone(),
            nodes: c.nodes.clone(),
            kind: c.kind.clone(),
            v: 0.0,
        })
        .collect();
    let mut ours = Device::new(mesh, si, process.doping(), contacts, spec.width);

    // Bias both to the same point and compare the drain current.
    for dev in [&mut reference, &mut ours] {
        dev.initial_guess();
        assert!(dev.solve().is_some(), "equilibrium solve failed");
        let gate = dev.contact("gate");
        let drain = dev.contact("drain");
        assert!(dev.ramp(gate, 1.0, 0.1).is_some(), "gate ramp failed");
        assert!(dev.ramp(drain, 0.1, 0.05).is_some(), "drain ramp failed");
    }
    let drain = reference.contact("drain");
    let (ne, pe) = reference.current(drain);
    let (no, po) = ours.current(drain);
    let i_ref = (ne + pe).abs();
    let i_ours = (no + po).abs();
    assert!(i_ref > 0.0, "no reference current to compare against");
    let rel = (i_ours - i_ref).abs() / i_ref;
    println!("drain current: reference {i_ref:e} A, CAD-built {i_ours:e} A ({rel:e})");
    assert!(rel < 1e-6, "drain currents differ by {rel:e}");
}

#[test]
fn a_bipolar_style_flow_stacks_implants() {
    // Three implants in sequence: a deep n well, a p base, an n emitter.
    // Each junction must land where the profiles cross.
    let p = Process::new(Dopant::acceptor("boron"), 1e21, [0.0, 10e-6], 5e-6)
        .film(Film::new("oxide", 50e-9).with_eps(3.9))
        .implant(Implant::blanket(
            Dopant::phosphorus(),
            Depth::Erfc { surface: 1e24, length: 2e-6 },
        ))
        .implant(Implant::through(
            Dopant::boron(),
            Depth::Erfc { surface: 5e24, length: 0.5e-6 },
            vec![[2e-6, 8e-6]],
            0.2e-6,
        ))
        .implant(Implant::through(
            Dopant::arsenic(),
            Depth::Erfc { surface: 1e26, length: 0.15e-6 },
            vec![[4e-6, 6e-6]],
            0.1e-6,
        ));

    // Under the emitter: n at the surface, then p, then n.
    let x = 5e-6;
    let net = |y: f64| {
        let (d, a) = p.doping_at(x, y);
        d - a
    };
    assert!(net(0.0) > 0.0, "emitter should be n-type");
    let junctions = p.junctions(x, 4000);
    assert_eq!(junctions.len(), 3, "emitter-base, base-well and well-substrate: {junctions:?}");
    let je = junctions[0];
    assert!(je > 0.1e-6 && je < 0.6e-6, "emitter junction at {je:e}");
    assert!(net(je + 0.1e-6) < 0.0, "base should be p-type below the emitter");
    // Deeper still, the well takes over again, and then the substrate.
    assert!(net(3e-6) > 0.0, "the n well should return");
    assert!(junctions[1] > je && junctions[2] > junctions[1]);
    assert!((p.junction_depth(x).unwrap() - je).abs() < 1e-12);

    // Outside the base window there is no base implant, so the well runs
    // to the surface.
    assert!(net_at(&p, 0.5e-6, 0.0) > 0.0, "the well should reach the surface outside the base");
    assert!(p.junction_depth(0.5e-6).is_some(), "well-substrate junction");
}

fn net_at(p: &Process, x: f64, y: f64) -> f64 {
    let (d, a) = p.doping_at(x, y);
    d - a
}
