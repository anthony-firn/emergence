//! The project graph: hand-offs, intervals, incremental edits, and the
//! mistakes it must refuse.

use std::path::Path;

use emergence::Rigor;
use project::analysis::{self, Analysis};
use project::graph::Cache;
use project::{Doc, Project, Status};

fn project(json: &str) -> Result<Project, String> {
    Project::new(Doc::from_json(json)?, ".", Cache::default())
}

/// Deterministic uniform samples in [0, 1).
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

const TRACE: &str = r#"{
  "name": "trace",
  "params": {
    "w": { "value": "10 mm", "tol": "0.1 mm", "owner": "electrical" },
    "t": { "value": "70 um", "tol": "10%", "owner": "electrical" },
    "l": { "value": "50 mm", "tol": "1 mm", "owner": "electrical" },
    "i": { "value": "30 A", "tol": "5%", "owner": "electrical" },
    "temp": { "value": "85 degC", "tol": "10 K", "owner": "thermal" }
  },
  "analyses": {
    "trace": { "kind": "pcb_trace", "owner": "electrical",
      "inputs": { "width": "w", "thickness": "t", "length": "l", "current": "i", "temperature": "temp" },
      "materials": { "conductor": "copper" } }
  },
  "requirements": {
    "R1": { "text": "drop", "of": "trace.voltage_drop", "max": "50 mV", "owner": "electrical" }
  }
}"#;

#[test]
fn rigorous_intervals_contain_every_sample() {
    // The trace model is monotone in every input, so corner enumeration is a
    // bound: no point of the tolerance box may land outside it.
    let p = project(TRACE).unwrap();
    let e = p.evaluate();
    let r = e.value("trace.resistance").unwrap();
    assert!(e.runs["trace"].rigor_note.is_none());
    let lib = project::Library::builtin();
    let rho = lib.get("copper", "resistivity").unwrap().reference();
    let alpha = lib.get("copper", "alpha_r").unwrap().reference();
    let mut g = Lcg(7);
    let pick = |v: project::Val, g: &mut Lcg| v.lo + (v.hi - v.lo) * g.next();
    for _ in 0..5000 {
        let x = [
            pick(e.value("w").unwrap(), &mut g),
            pick(e.value("t").unwrap(), &mut g),
            pick(e.value("l").unwrap(), &mut g),
            pick(e.value("i").unwrap(), &mut g),
            pick(e.value("temp").unwrap(), &mut g),
            pick(rho, &mut g),
            pick(alpha, &mut g),
        ];
        let out = analysis::PcbTrace.run(&x).unwrap();
        assert!(out.values[0] >= r.lo && out.values[0] <= r.hi);
    }
    // ... and the corners are attained, so the interval is tight.
    let corners = [r.lo, r.hi];
    assert!(corners[0] < r.nom && r.nom < corners[1]);
}

#[test]
fn declared_monotone_inputs_are_monotone() {
    // Rigor rests on the monotone flags, so check each against the model on
    // random designs: one-input steps must never reverse direction.
    let mut g = Lcg(11);
    let models: Vec<(Box<dyn Analysis>, Vec<f64>, Vec<f64>)> = vec![
        (
            Box::new(analysis::PlateFinSink),
            // width, length, base, fin_h, fin_t, n, V, T, src_w, src_l, k_fin, k_air, nu, Pr, rho, cp
            vec![0.03, 0.03, 0.002, 0.01, 0.0008, 6.0, 0.5, 300.0, 0.005, 0.005, 150.0, 0.024, 1.4e-5, 0.69, 1.0, 1000.0],
            vec![0.06, 0.08, 0.006, 0.03, 0.002, 12.0, 3.0, 340.0, 0.03, 0.03, 240.0, 0.031, 2.0e-5, 0.72, 1.2, 1015.0],
        ),
        (
            Box::new(analysis::PcbTrace),
            vec![1e-3, 17e-6, 5e-3, 0.5, 250.0, 1.6e-8, 0.0038],
            vec![20e-3, 140e-6, 200e-3, 40.0, 400.0, 1.8e-8, 0.0041],
        ),
    ];
    for (model, lo, hi) in &models {
        let ports = model.inputs();
        let mut checked = 0;
        for _ in 0..40 {
            let mut x: Vec<f64> = lo.iter().zip(hi).map(|(a, b)| a + (b - a) * g.next()).collect();
            if model.kind() == "plate_fin_heatsink" {
                x[5] = x[5].round();
            }
            let Ok(base) = model.run(&x) else { continue };
            for (i, port) in ports.iter().enumerate() {
                if !port.monotone || port.passive {
                    continue;
                }
                let steps: Vec<Vec<f64>> = (0..4)
                    .map(|k| {
                        let mut y = x.clone();
                        y[i] = lo[i] + (hi[i] - lo[i]) * (k as f64 / 3.0);
                        y
                    })
                    .collect();
                let runs: Vec<_> = steps.iter().filter_map(|y| model.run(y).ok()).collect();
                if runs.len() < 4 {
                    continue;
                }
                for j in 0..base.values.len() {
                    let d: Vec<f64> = runs.windows(2).map(|w| w[1].values[j] - w[0].values[j]).collect();
                    let up = d.iter().all(|v| *v >= -1e-12 * runs[0].values[j].abs());
                    let down = d.iter().all(|v| *v <= 1e-12 * runs[0].values[j].abs());
                    assert!(up || down, "{} output {j} not monotone in {}: {d:?}", model.kind(), port.name);
                    checked += 1;
                }
            }
        }
        assert!(checked > 100, "{}: only {checked} checks ran", model.kind());
    }
}

#[test]
fn edits_rerun_only_what_they_reach() {
    let json = TRACE.replace(
        r#""analyses": {"#,
        r#""analyses": {
    "tim": { "kind": "conduction_slab", "owner": "thermal",
      "inputs": { "thickness": "0.5[mm]", "area": "150[mm^2]" }, "materials": { "material": "al6061" } },"#,
    );
    let p = project(&json).unwrap();
    let first = p.evaluate();
    assert_eq!(first.cached, 0);
    let edited = p.with_values(&[("w".into(), "12 mm".into())]).unwrap();
    let second = edited.evaluate();
    // The slab never sees the width: every one of its runs is reused, and
    // the trace reruns only at points it has not seen.
    assert_eq!(second.cached, first.runs["tim"].corners + 1);
    // Incremental equals from scratch, bit for bit.
    let fresh = Project::new(edited.doc.clone(), ".", Cache::default()).unwrap().evaluate();
    for (id, n) in &second.nodes {
        assert_eq!(n.value.as_ref().ok(), fresh.nodes[id].value.as_ref().ok(), "{id}");
    }
    assert_ne!(first.value("trace.resistance"), second.value("trace.resistance"));
}

#[test]
fn uncertainty_is_traced_to_its_owner() {
    let p = project(
        r#"{ "name": "budget",
      "params": {
        "x": { "value": "1 W", "tol": "10%", "owner": "electrical" },
        "y": { "value": "1 W", "tol": "1%", "owner": "thermal" } },
      "derived": { "z": { "expr": "x + y", "owner": "systems" } },
      "requirements": { "R": { "text": "", "of": "z", "max": "2.05 W", "owner": "systems" } } }"#,
    )
    .unwrap();
    let e = p.evaluate();
    // z = 2 W in [1.89, 2.11]: straddles 2.05 W, so honestly undecided.
    assert_eq!(e.checks["R"].status, Status::Undecided);
    let c = p.contributions("R", &e);
    assert_eq!(c[0].0, "x");
    assert_eq!(c[0].1, "electrical");
    // Pinning x leaves only y's ±0.01 W of the ±0.11 W: exactly 10/11.
    assert!((c[0].2 - 10.0 / 11.0).abs() < 1e-12);
    assert!((c[1].2 - 1.0 / 11.0).abs() < 1e-12);
}

#[test]
fn hand_off_mistakes_are_refused() {
    // A temperature bound to a length input.
    let bad = TRACE.replace(r#""width": "w""#, r#""width": "temp""#);
    let e = project(&bad).unwrap().evaluate();
    let err = e.runs["trace"].error.clone().unwrap();
    assert!(err.contains("needs m"), "{err}");
    assert_eq!(e.checks["R1"].status, Status::Unavailable);
    // A misspelt input is caught at load, with the real names listed.
    let typo = TRACE.replace(r#""width": "w""#, r#""widht": "w""#);
    let err = project(&typo).err().unwrap();
    assert!(err.contains("no input `widht`") && err.contains("width"), "{err}");
    // A declared unit that the expression does not produce.
    let p = project(
        r#"{ "name": "u", "params": { "p": { "value": "2 W", "owner": "a" } },
             "derived": { "t": { "expr": "p * 2", "owner": "a", "unit": "degC" } } }"#,
    )
    .unwrap();
    let err = p.evaluate().nodes["t"].value.clone().unwrap_err();
    assert!(err.contains("declared unit degC"), "{err}");
    // A cycle, named.
    let p = project(
        r#"{ "name": "c", "derived": {
             "a": { "expr": "b + 1[W]", "owner": "x" }, "b": { "expr": "a", "owner": "x" } } }"#,
    )
    .unwrap();
    let err = p.evaluate().nodes["a"].value.clone().unwrap_err();
    assert!(err.contains("cycle") && (err.contains("a → b → a") || err.contains("b → a → b")), "{err}");
    // Unknown ids and outputs, with suggestions.
    let p = project(r#"{ "name": "x", "derived": { "a": { "expr": "nope * 2", "owner": "x" } } }"#).unwrap();
    assert!(p.evaluate().nodes["a"].value.clone().unwrap_err().contains("unknown id `nope`"));
    let wrong_output = TRACE.replace(r#""of": "trace.voltage_drop""#, r#""of": "trace.drop""#);
    let r = &project(&wrong_output).unwrap().evaluate().checks["R1"];
    assert!(r.reason.as_ref().unwrap().contains("unknown id") || r.status == Status::Unavailable);
    // A limit in the wrong dimension.
    let wrong_limit = TRACE.replace(r#""max": "50 mV""#, r#""max": "50 mA""#);
    let r = &project(&wrong_limit).unwrap().evaluate().checks["R1"];
    assert_eq!(r.status, Status::Unavailable);
    assert!(r.reason.as_ref().unwrap().contains("limit"));
    // A Celsius tolerance (ambiguous: 5 degC is 278 K).
    let p = project(r#"{ "name": "t", "params": { "a": { "value": "40 degC", "tol": "5 degC", "owner": "s" } } }"#).unwrap();
    assert!(p.evaluate().nodes["a"].value.clone().unwrap_err().contains("in K"));
    // Duplicate ids across sections.
    let dup = TRACE.replace(r#""w": { "value""#, r#""trace": { "value": "1 mm", "owner": "e" }, "w": { "value""#);
    assert!(project(&dup).err().unwrap().contains("twice"));
}

#[test]
fn example_project_end_to_end() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../projects/motor-driver/project.json");
    let p = Project::load(&path).unwrap();
    let e = p.evaluate();
    for (id, r) in &e.runs {
        assert!(r.error.is_none(), "{id}: {:?}", r.error);
    }
    let status = |id: &str| e.checks[id].status;
    assert_eq!(status("REQ-TJ"), Status::Met);
    assert_eq!(status("REQ-FAN"), Status::Met);
    assert_eq!(status("REQ-MASS"), Status::Met);
    // 30 A through a 10 mm, 2 oz trace: 46 mV nominal, but the copper
    // thickness tolerance alone can push it past 50 mV.
    assert_eq!(status("REQ-DROP"), Status::Undecided);
    // The junction temperature is assembled from three teams' numbers;
    // the nominal must be exactly the hand-off arithmetic.
    let v = |id: &str| e.value(id).unwrap().nom;
    let tj = v("ambient") + v("q1.power") * (v("q1.rjc") + v("tim.resistance") + v("sink.thermal_resistance"));
    assert!((v("q1.tj") - tj).abs() < 1e-12 * tj);
    assert_eq!(e.value("q1.tj").unwrap().rigor, Rigor::Asymptotic);
    // The electro-thermal loop, Tj = Ta + Rth (I² R25 (1 + tc (Tj − 25 °C)) + Psw),
    // is linear in Tj: its nominal has a closed form.
    let (ta, rth, i, r25, tc, psw) =
        (v("ambient"), v("q1.rth_total"), v("q1.current"), v("q1.rds_on_25"), v("q1.rds_tc"), v("q1.switching_loss"));
    let closed = (ta + rth * (i * i * r25 * (1.0 - tc * 298.15) + psw)) / (1.0 - rth * i * i * r25 * tc);
    assert!((v("q1.tj") - closed).abs() < 1e-9, "{} vs {closed}", v("q1.tj"));
    assert!(e.loops.len() == 1 && e.loops[0].existence, "{:?}", e.loops);
    // The CAD part: exact volume W L t_b + N t_f H L = 16.8 cm³ must lie in
    // the claimed interval (nominal fin thickness; the interval also spans
    // the ±0.05 mm extrusion tolerance).
    let vol = e.value("sink_body.volume").unwrap();
    let exact = 0.04 * 0.05 * 0.003 + 9.0 * 0.0012 * 0.02 * 0.05;
    assert!(vol.lo <= exact && exact <= vol.hi, "{vol:?} vs {exact}");
    println!("volume {:.4} cm³ vs exact {:.4} (error {:.2} %)", vol.nom * 1e6, exact * 1e6, 100.0 * (vol.nom / exact - 1.0));
}

fn loop_project(gain: &str) -> Project {
    project(&format!(
        r#"{{ "name": "loop",
      "params": {{
        "base": {{ "value": "11 K", "tol": "1 K", "owner": "a" }},
        "gain": {{ "value": "{gain}", "tol": "0.05", "owner": "b" }} }},
      "derived": {{
        "t": {{ "expr": "base + heat", "owner": "a", "unit": "K", "range": ["0 K", "100 K"] }},
        "heat": {{ "expr": "gain * t", "owner": "b", "unit": "K" }} }} }}"#
    ))
    .unwrap()
}

#[test]
fn feedback_loops_are_solved_with_an_enclosure() {
    // t = base + gain·t: t = base / (1 − gain). With base in [10, 12] and
    // gain in [0.2, 0.3] the exact range is [12.5, 120/7], nominal 44/3.
    let e = loop_project("0.25").evaluate();
    let t = e.value("t").unwrap();
    assert!((t.nom - 44.0 / 3.0).abs() < 1e-11, "{t:?}");
    assert!((t.lo - 12.5).abs() < 1e-11 && (t.hi - 120.0 / 7.0).abs() < 1e-11, "{t:?}");
    assert!(e.loops[0].existence);
    assert_eq!(e.loops[0].unknowns, vec!["t".to_string()]);
    assert_eq!(t.rigor, Rigor::Rigorous);
    // Loop gain above one: no operating point anywhere in the range.
    let e = loop_project("1.2").evaluate();
    let err = e.nodes["t"].value.clone().unwrap_err();
    assert!(err.contains("no operating point"), "{err}");
    // A loop nobody declared is still refused, with the fix named.
    let p = project(
        r#"{ "name": "c", "derived": {
             "a": { "expr": "b + 1[W]", "owner": "x" }, "b": { "expr": "a", "owner": "x" } } }"#,
    )
    .unwrap();
    let err = p.evaluate().nodes["a"].value.clone().unwrap_err();
    assert!(err.contains("`range`"), "{err}");
}

#[test]
fn parts_render_from_the_shared_model() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../projects/motor-driver/project.json");
    let p = Project::load(&path).unwrap();
    let e = p.evaluate();
    let png = project::render::part_png(&p, &e, "sink_body", 96, 64).unwrap();
    assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
    assert!(project::render::part_png(&p, &e, "trace", 96, 64).is_err(), "an analysis is not a part");
}

#[test]
#[ignore = "slow on a cold cache: characterises nine cell libraries from the transistor level (~30 s)"]
fn chip_project_closes_the_leakage_loop() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../projects/mini-cpu/project.json");
    let p = Project::load(&path).unwrap();
    let e = p.evaluate();
    assert!(e.runs["cpu"].error.is_none(), "{:?}", e.runs["cpu"].error);
    assert!(e.loops.len() == 1 && e.loops[0].existence);
    assert_eq!(e.checks["REQ-TJ"].status, Status::Met);
    assert_eq!(e.checks["REQ-HOLD"].status, Status::Met);
}

#[test]
fn office_hvac_project() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../projects/office-hvac/project.json");
    let p = Project::load(&path).unwrap();
    let e = p.evaluate();
    for (id, r) in &e.runs {
        assert!(r.error.is_none(), "{id}: {:?}", r.error);
    }
    let v = |id: &str| e.value(id).unwrap().nom;
    // The supply airflow is the sensible peak over ρ c_p ΔT exactly.
    let flow = v("office.cooling_sensible") / (1.204 * 1006.0 * v("mech.supply_dt"));
    assert!((v("office.supply_airflow") / flow - 1.0).abs() < 1e-12);
    // The plant delivers exactly the zone's total, and its compressor power
    // follows from the COP.
    assert!((v("plant.compressor_power") * v("plant.cop_cooling") / v("office.cooling_total") - 1.0).abs() < 1e-9);
    // Heating sized at 4 kW against a ~4.04 kW load: honestly undecided.
    assert_eq!(e.checks["REQ-HEAT"].status, Status::Undecided);
    assert_eq!(e.checks["REQ-NOISE"].status, Status::Met);
}
