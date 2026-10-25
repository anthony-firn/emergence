//! `fluid project …`: the command line for teams and CI.

use std::path::{Path, PathBuf};
use std::process::Command;

use crate::diff;
use crate::doc::Doc;
use crate::graph::{Project, Status};
use crate::{analysis, html, report};

pub const USAGE: &str = "\
usage: fluid project report FILE [--why] [--html OUT.html] [--set ID=VALUE ...]
       fluid project check  FILE [--set ID=VALUE ...]      exit 0 all met, 1 any failed, 3 undecided/unavailable
       fluid project diff   OLD NEW                         what changed, what it moved, who caused it
       fluid project diff   FILE --rev REV                  against FILE as committed at git revision REV
       fluid project set    FILE ID VALUE                   edit a parameter (rewrites FILE canonically)
       fluid project fmt    FILE                            rewrite FILE in canonical (sorted) form
       fluid project export FILE PART OUT.stl [--set ID=VALUE ...]
       fluid project render FILE PART OUT.png [--width 900] [--height 600] [--set ID=VALUE ...]
       fluid project serve  FILE [--port 8080] [--host 127.0.0.1]   live shared editor in the browser
       fluid project kinds                                  analysis kinds, their inputs and outputs";

/// Run a subcommand; returns the process exit code.
pub fn run(args: &[String]) -> i32 {
    match dispatch(args) {
        Ok(code) => code,
        Err(e) => {
            eprintln!("error: {e}");
            2
        }
    }
}

fn sets(args: &[String]) -> Result<Vec<(String, String)>, String> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < args.len() {
        if args[i] == "--set" {
            let pair = args.get(i + 1).ok_or("--set needs ID=VALUE")?;
            let (id, v) = pair.split_once('=').ok_or_else(|| format!("--set needs ID=VALUE, found {pair:?}"))?;
            out.push((id.to_string(), v.to_string()));
            i += 1;
        }
        i += 1;
    }
    Ok(out)
}

fn flag<'a>(args: &'a [String], name: &str) -> Option<&'a str> {
    args.iter().position(|a| a == name).and_then(|i| args.get(i + 1)).map(String::as_str)
}

fn load(path: &str, args: &[String]) -> Result<Project, String> {
    let p = Project::load(path)?;
    let s = sets(args)?;
    if s.is_empty() {
        Ok(p)
    } else {
        p.with_values(&s)
    }
}

fn dispatch(args: &[String]) -> Result<i32, String> {
    let cmd = args.first().map(String::as_str).unwrap_or("");
    let file = || args.get(1).filter(|a| !a.starts_with("--")).ok_or_else(|| USAGE.to_string());
    match cmd {
        "report" => {
            let p = load(file()?, args)?;
            let e = p.evaluate();
            print!("{}", report::text(&p, &e, args.iter().any(|a| a == "--why")));
            if let Some(out) = flag(args, "--html") {
                std::fs::write(out, html::page(&p, &e)).map_err(|e| format!("{out}: {e}"))?;
                eprintln!("wrote {out}");
            }
            Ok(0)
        }
        "check" => {
            let p = load(file()?, args)?;
            let e = p.evaluate();
            let mut code = 0;
            for c in e.checks.values() {
                println!("{} {:<12} {}", report::status_mark(c.status), c.id, c.status.name());
                code = match (code, c.status) {
                    (_, Status::Failed) | (1, _) => 1,
                    (_, Status::Undecided | Status::Unavailable) => 3,
                    (c, Status::Met) => c,
                };
            }
            Ok(code)
        }
        "diff" => {
            let a = file()?;
            let (old, new) = if let Some(rev) = flag(args, "--rev") {
                let new = Project::load(a)?;
                let dir = checkout(a, rev, &new.doc)?;
                let name = Path::new(a).file_name().unwrap();
                let old = Project::new(Doc::read(dir.join(name))?, &dir, new.cache());
                let _ = std::fs::remove_dir_all(&dir);
                (old?, new)
            } else {
                let b = args.get(2).ok_or_else(|| USAGE.to_string())?;
                let old = Project::load(a)?;
                let new = Project::new(Doc::read(b)?, Path::new(b).parent().unwrap_or(Path::new(".")), old.cache())?;
                (old, new)
            };
            print!("{}", diff::text(&diff::diff(&old, &new)));
            Ok(0)
        }
        "set" => {
            let (path, id, value) = match &args[1..] {
                [p, i, v, ..] => (p, i, v),
                _ => return Err(USAGE.into()),
            };
            let mut doc = Doc::read(path)?;
            let param = doc.params.get_mut(id).ok_or_else(|| format!("`{id}` is not a parameter"))?;
            let old = std::mem::replace(&mut param.value, value.clone());
            // Refuse an edit that does not load.
            let dir = Path::new(path).parent().unwrap_or(Path::new(".")).to_path_buf();
            Project::new(doc.clone(), &dir, Default::default())?;
            crate::graph::param_value(value, doc.params[id].tol.as_ref())?;
            std::fs::write(path, doc.to_json()).map_err(|e| format!("{path}: {e}"))?;
            eprintln!("{id}: {old} → {value}");
            Ok(0)
        }
        "fmt" => {
            let path = file()?;
            let doc = Doc::read(path)?;
            std::fs::write(path, doc.to_json()).map_err(|e| format!("{path}: {e}"))?;
            Ok(0)
        }
        "export" => {
            let (path, part, out) = match &args[1..] {
                [p, part, o, ..] => (p, part, o),
                _ => return Err(USAGE.into()),
            };
            let p = load(path, args)?;
            export(&p, part, out)?;
            Ok(0)
        }
        "render" => {
            let (path, part, out) = match &args[1..] {
                [p, part, o, ..] => (p, part, o),
                _ => return Err(USAGE.into()),
            };
            let p = load(path, args)?;
            let e = p.evaluate();
            let w: usize = flag(args, "--width").unwrap_or("900").parse().map_err(|_| "--width needs a number")?;
            let h: usize = flag(args, "--height").unwrap_or("600").parse().map_err(|_| "--height needs a number")?;
            let png = crate::render::part_png(&p, &e, part, w, h)?;
            std::fs::write(out, png).map_err(|e| format!("{out}: {e}"))?;
            eprintln!("wrote {out}");
            Ok(0)
        }
        "serve" => {
            let path = file()?;
            let host = flag(args, "--host").unwrap_or("127.0.0.1");
            let port = flag(args, "--port").unwrap_or("8080");
            crate::serve::serve(Path::new(path), &format!("{host}:{port}"))?;
            Ok(0)
        }
        "kinds" => {
            for k in analysis::KINDS {
                let a = analysis::builtin(k).unwrap();
                println!("{k}: {}", a.summary());
                for port in a.inputs() {
                    let from = port.from.map(|f| format!(" (default: material slot `{}`.{})", f.slot, f.prop)).unwrap_or_default();
                    println!("  in  {:<18} [{}] {}{from}", port.name, port.unit, port.doc);
                }
                for o in a.outputs() {
                    println!("  out {:<18} [{}] {}", o.name, o.unit, o.doc);
                }
                for s in a.assumptions() {
                    println!("  assumes: {s}");
                }
                println!();
            }
            Ok(0)
        }
        _ => Err(USAGE.into()),
    }
}

/// Write a part's solid, at the project's nominal values, as STL.
pub fn export(p: &Project, part: &str, out: &str) -> Result<(), String> {
    let model = p.part(part).ok_or_else(|| {
        format!("`{part}` is not a part (parts: {})", p.doc.parts.keys().cloned().collect::<Vec<_>>().join(", "))
    })?;
    let e = p.evaluate();
    let run = &e.runs[part];
    if let Some(err) = &run.error {
        return Err(err.clone());
    }
    let x: Vec<f64> = run.inputs.iter().take(model.params.len()).map(|(_, _, v)| v.map(|v| v.nom).unwrap_or(f64::NAN)).collect();
    let mesh = cad::sim::write_stl(&model.solid(&x)?, model.cell, out)?;
    eprintln!("wrote {out}: {} triangles", mesh.triangles.len());
    Ok(())
}

/// Recreate the project file and every part's CAD model as committed at
/// `rev`, in a fresh temporary directory.
fn checkout(path: &str, rev: &str, doc: &Doc) -> Result<PathBuf, String> {
    let file = Path::new(path);
    let dir = file.parent().filter(|d| !d.as_os_str().is_empty()).unwrap_or(Path::new("."));
    let tmp = std::env::temp_dir().join(format!("fluid-project-{}-{}", std::process::id(), rev.replace('/', "_")));
    std::fs::create_dir_all(&tmp).map_err(|e| e.to_string())?;
    let show = |rel: &str| -> Result<Vec<u8>, String> {
        let out = Command::new("git")
            .arg("-C")
            .arg(dir)
            .arg("show")
            .arg(format!("{rev}:./{rel}"))
            .output()
            .map_err(|e| format!("git: {e}"))?;
        if !out.status.success() {
            return Err(format!("git show {rev}:{rel}: {}", String::from_utf8_lossy(&out.stderr).trim()));
        }
        Ok(out.stdout)
    };
    let name = file.file_name().unwrap().to_string_lossy().to_string();
    let text = show(&name)?;
    std::fs::write(tmp.join(&name), &text).map_err(|e| e.to_string())?;
    let old = Doc::from_json(&String::from_utf8_lossy(&text))?;
    for part in old.parts.values().chain(doc.parts.values()) {
        if let Ok(bytes) = show(&part.model) {
            let target = tmp.join(&part.model);
            if let Some(parent) = target.parent() {
                std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
            }
            std::fs::write(target, bytes).map_err(|e| e.to_string())?;
        }
    }
    Ok(tmp)
}
