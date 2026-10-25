//! Reviewing a change: semantic diff, the CI gate, and edits from the CLI.

use std::path::{Path, PathBuf};
use std::process::Command;

use project::graph::Cache;
use project::{diff, Doc, Project, Status};

const BOARD: &str = r#"{
  "format": "fluid-project/1",
  "name": "board",
  "params": {
    "w": { "value": "10 mm", "tol": "0.1 mm", "owner": "electrical" },
    "i": { "value": "15 A", "owner": "electrical" },
    "temp": { "value": "85 degC", "owner": "thermal" },
    "rth": { "value": "30 K/W", "owner": "thermal" }
  },
  "derived": {
    "rise": { "expr": "trace.power * rth", "owner": "thermal", "unit": "K" }
  },
  "analyses": {
    "trace": { "kind": "pcb_trace", "owner": "electrical",
      "inputs": { "width": "w", "thickness": "35[um]", "length": "40[mm]", "current": "i", "temperature": "temp" },
      "materials": { "conductor": "copper" } }
  },
  "requirements": {
    "DROP": { "text": "drop", "of": "trace.voltage_drop", "max": "50 mV", "owner": "electrical" },
    "RISE": { "text": "local rise", "of": "rise", "max": "30 K", "owner": "thermal" }
  }
}"#;

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("project-test-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn load(json: &str) -> Project {
    Project::new(Doc::from_json(json).unwrap(), ".", Cache::default()).unwrap()
}

#[test]
fn a_change_is_traced_to_the_requirements_it_moves() {
    let old = load(BOARD);
    // 37 mV and 17 K rise at 10 mm; 93 mV and 42 K at 4 mm.
    let new = load(&BOARD.replace(r#""value": "10 mm""#, r#""value": "4 mm""#));
    assert_eq!(old.evaluate().checks["DROP"].status, Status::Met);
    let d = diff::diff(&old, &new);
    assert_eq!(d.changes.len(), 1);
    assert_eq!((d.changes[0].id.as_str(), d.changes[0].owner.as_str()), ("w", "electrical"));
    // Electrical narrowed its trace: its own drop limit fails, and so does
    // thermal's local temperature rise, which electrical never looked at.
    let flip = |id: &str| d.verdicts.iter().find(|v| v.id == id).unwrap();
    assert_eq!(flip("DROP").after, Status::Failed);
    assert_eq!(flip("RISE").before, Status::Met);
    assert_eq!(flip("RISE").after, Status::Failed);
    assert_eq!(flip("RISE").causes, vec![("w".to_string(), "electrical".to_string())]);
    let text = diff::text(&d);
    assert!(text.contains("RISE") && text.contains("met → FAILED") && text.contains("caused by: w (electrical)"), "{text}");
    // An unrelated change moves nothing.
    let same = load(&BOARD.replace(r#""name": "board""#, r#""name": "board v2""#));
    assert!(diff::diff(&old, &same).verdicts.is_empty());
}

fn cli(args: &[&str]) -> i32 {
    project::cli::run(&args.iter().map(|s| s.to_string()).collect::<Vec<_>>())
}

#[test]
fn check_is_a_ci_gate_and_set_refuses_bad_edits() {
    let dir = scratch("cli");
    let file = dir.join("project.json");
    std::fs::write(&file, BOARD).unwrap();
    let f = file.to_str().unwrap();
    assert_eq!(cli(&["check", f]), 0);
    // A what-if from the command line does not touch the file.
    assert_eq!(cli(&["check", f, "--set", "i=40 A"]), 1);
    // A limit at the nominal value lies inside the interval: undecided, exit 3.
    let nominal = load(BOARD).evaluate().value("trace.voltage_drop").unwrap().nom;
    let tight = BOARD.replace(r#""max": "50 mV""#, &format!(r#""max": "{} mV""#, nominal * 1e3));
    std::fs::write(&file, &tight).unwrap();
    assert_eq!(cli(&["check", f]), 3);
    // `set` writes canonical JSON, and refuses values that do not parse.
    assert_eq!(cli(&["set", f, "i", "10 A"]), 0);
    let doc = Doc::read(&file).unwrap();
    assert_eq!(doc.params["i"].value, "10 A");
    assert_eq!(std::fs::read_to_string(&file).unwrap(), doc.to_json());
    assert_eq!(cli(&["set", f, "i", "ten amps"]), 2);
    assert_eq!(cli(&["set", f, "nope", "1 A"]), 2);
    assert_eq!(Doc::read(&file).unwrap().params["i"].value, "10 A");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn diff_against_a_git_revision() {
    if Command::new("git").arg("--version").output().is_err() {
        return;
    }
    let dir = scratch("git");
    let git = |args: &[&str]| {
        let ok = Command::new("git").arg("-C").arg(&dir).args(args).output().unwrap().status.success();
        assert!(ok, "git {args:?}");
    };
    git(&["init", "-q"]);
    std::fs::write(dir.join("project.json"), BOARD).unwrap();
    git(&["add", "project.json"]);
    git(&["-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "v1"]);
    let f = dir.join("project.json");
    assert_eq!(cli(&["set", f.to_str().unwrap(), "rth", "60 K/W"]), 0);
    // The library call behind `fluid project diff FILE --rev HEAD`.
    let new = Project::load(&f).unwrap();
    let out = Command::new("git").arg("-C").arg(&dir).args(["show", "HEAD:./project.json"]).output().unwrap();
    let old = Project::new(Doc::from_json(&String::from_utf8_lossy(&out.stdout)).unwrap(), Path::new(&dir), new.cache()).unwrap();
    let d = diff::diff(&old, &new);
    assert_eq!(d.changes[0].id, "rth");
    assert!(d.verdicts.iter().any(|v| v.id == "RISE" && v.after == Status::Failed));
    assert_eq!(cli(&["diff", f.to_str().unwrap(), "--rev", "HEAD"]), 0);
    let _ = std::fs::remove_dir_all(&dir);
}
