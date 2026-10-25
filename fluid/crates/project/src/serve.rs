//! `fluid project serve`: the project as a live, shared web page.
//!
//! Everyone on the team opens the same address. A parameter edited in one
//! browser is checked (units, tolerance, a project that still loads),
//! evaluated, written to the project file, and pushed to every other viewer
//! with what it did to the requirements. The file stays the source of truth:
//! edits made with an editor or a `git pull` are picked up the same way, and
//! the CLI, git and CI keep working alongside.
//!
//! Standard library only: one thread, one request at a time, so edits are
//! serialised and two people can never interleave half-applied changes.
//! There is no login. Bind to loopback unless the network is trusted.

use std::io::{BufRead, BufReader, Read, Write};
use std::net::{SocketAddr, TcpListener, TcpStream};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::doc::Doc;
use crate::graph::{Evaluation, Project};
use crate::html::{self, Activity, Live};
use crate::report;

struct State {
    path: PathBuf,
    project: Project,
    eval: Evaluation,
    version: u64,
    activity: Vec<Activity>,
    mtime: Option<SystemTime>,
    /// Part renders for the current version.
    renders: std::collections::HashMap<(String, u64), Vec<u8>>,
}

fn mtime(path: &Path) -> Option<SystemTime> {
    std::fs::metadata(path).and_then(|m| m.modified()).ok()
}

fn clock() -> String {
    let s = SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0) % 86_400;
    format!("{:02}:{:02}:{:02} UTC", s / 3600, s / 60 % 60, s % 60)
}

/// What a change did to the requirements, one line each: status changes
/// first, then moved values.
pub fn impact(before: &Evaluation, after: &Evaluation) -> Vec<String> {
    let mut flips = Vec::new();
    let mut moved = Vec::new();
    for (id, a) in &after.checks {
        let Some(b) = before.checks.get(id) else { continue };
        let f = |c: &crate::graph::Check| c.value.map(|v| report::val(&v, c.unit.as_deref())).unwrap_or_else(|| "unavailable".into());
        if a.status != b.status {
            let why = match (&a.reason, a.value) {
                (Some(r), None) => format!(": {r}"),
                _ => String::new(),
            };
            flips.push(format!("{id} [{}] {} → {}: {} → {}{why}", a.owner, b.status.name(), a.status.name(), f(b), f(a)));
        } else if a.value.map(|v| v.nom) != b.value.map(|v| v.nom) {
            moved.push(format!("{id} [{}] still {}: {} → {}", a.owner, a.status.name(), f(b), f(a)));
        }
    }
    flips.extend(moved);
    flips
}

impl State {
    fn load(path: &Path, cache: crate::graph::Cache) -> Result<(Project, Evaluation), String> {
        let doc = Doc::read(path)?;
        let p = Project::new(doc, path.parent().unwrap_or(Path::new(".")), cache)?;
        let e = p.evaluate();
        Ok((p, e))
    }

    /// Pick up edits made outside the server.
    fn sync(&mut self) {
        let now = mtime(&self.path);
        if now == self.mtime {
            return;
        }
        self.mtime = now;
        match State::load(&self.path, self.project.cache()) {
            Ok((p, e)) => {
                let what = if p.doc == self.project.doc { None } else { Some("changed the project file outside the editor") };
                if let Some(what) = what {
                    let imp = impact(&self.eval, &e);
                    self.activity.push(Activity { time: clock(), who: "someone".into(), what: what.into(), impact: imp });
                    self.project = p;
                    self.eval = e;
                    self.version += 1;
                }
            }
            Err(err) => {
                self.activity.push(Activity {
                    time: clock(),
                    who: "someone".into(),
                    what: format!("saved a project file that does not load (still showing the last good one): {err}"),
                    impact: vec![],
                });
                self.version += 1;
            }
        }
    }

    fn set(&mut self, id: &str, value: &str, who: &str) -> Result<String, String> {
        let def = self.project.doc.params.get(id).ok_or_else(|| format!("`{id}` is not a parameter"))?;
        let old = def.value.clone();
        crate::graph::param_value(value, def.tol.as_ref())?;
        let next = self.project.with_values(&[(id.to_string(), value.to_string())])?;
        let e = next.evaluate();
        // Refuse an edit that breaks the unit of what it feeds.
        if let Some(Err(err)) = e.nodes.get(id).map(|n| &n.value) {
            return Err(err.clone());
        }
        std::fs::write(&self.path, next.doc.to_json()).map_err(|e| format!("{}: {e}", self.path.display()))?;
        self.mtime = mtime(&self.path);
        let imp = impact(&self.eval, &e);
        let owner = &next.doc.params[id].owner;
        let summary = if imp.is_empty() { format!("{id} = {value}: no requirement moved") } else { format!("{id} = {value}: {}", imp.join("; ")) };
        self.activity.push(Activity { time: clock(), who: who.into(), what: format!("set {id} ({owner}) {old} → {value}"), impact: imp });
        self.project = next;
        self.eval = e;
        self.version += 1;
        Ok(summary)
    }
}

struct Request {
    method: String,
    path: String,
    query: String,
    body: Vec<u8>,
}

fn read_request(stream: &mut TcpStream) -> Option<Request> {
    let mut r = BufReader::new(stream.try_clone().ok()?);
    let mut line = String::new();
    r.read_line(&mut line).ok()?;
    let mut parts = line.split_whitespace();
    let method = parts.next()?.to_string();
    let target = parts.next()?.to_string();
    let mut length = 0usize;
    loop {
        let mut h = String::new();
        if r.read_line(&mut h).ok()? == 0 || h == "\r\n" || h == "\n" {
            break;
        }
        if let Some((k, v)) = h.split_once(':') {
            if k.trim().eq_ignore_ascii_case("content-length") {
                length = v.trim().parse().ok()?;
            }
        }
    }
    if length > 1 << 20 {
        return None;
    }
    let mut body = vec![0; length];
    r.read_exact(&mut body).ok()?;
    let (path, query) = target.split_once('?').map(|(a, b)| (a.to_string(), b.to_string())).unwrap_or((target, String::new()));
    Some(Request { method, path, query, body })
}

fn respond(stream: &mut TcpStream, status: &str, kind: &str, body: &[u8]) {
    let head = format!(
        "HTTP/1.1 {status}\r\nContent-Type: {kind}\r\nContent-Length: {}\r\nCache-Control: no-store\r\nConnection: close\r\n\r\n",
        body.len()
    );
    let _ = stream.write_all(head.as_bytes());
    let _ = stream.write_all(body);
}

fn decode(s: &str) -> String {
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        match b[i] {
            b'+' => out.push(b' '),
            b'%' if i + 2 < b.len() => {
                if let Ok(v) = u8::from_str_radix(&s[i + 1..i + 3], 16) {
                    out.push(v);
                    i += 2;
                }
            }
            c => out.push(c),
        }
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

fn handle(state: &mut State, mut stream: TcpStream) {
    let Some(req) = read_request(&mut stream) else { return };
    state.sync();
    let json = |v: serde_json::Value| serde_json::to_vec(&v).unwrap();
    match (req.method.as_str(), req.path.as_str()) {
        ("GET", "/") | ("GET", "/index.html") => {
            let live = Live { version: state.version, activity: &state.activity };
            let page = html::render(&state.project, &state.eval, Some(&live));
            respond(&mut stream, "200 OK", "text/html; charset=utf-8", page.as_bytes());
        }
        ("GET", "/api/state") => {
            respond(&mut stream, "200 OK", "application/json", &json(serde_json::json!({ "version": state.version })));
        }
        ("GET", "/api/why") => {
            let req_id = req.query.split('&').find_map(|kv| kv.strip_prefix("req=")).map(decode).unwrap_or_default();
            let body = if state.project.doc.requirements.contains_key(&req_id) {
                html::contributions_html(&state.project, &state.eval, &req_id)
            } else {
                String::new()
            };
            respond(&mut stream, "200 OK", "text/html; charset=utf-8", body.as_bytes());
        }
        ("GET", "/api/render") => {
            let part = req.query.split('&').find_map(|kv| kv.strip_prefix("part=")).map(decode).unwrap_or_default();
            let key = (part.clone(), state.version);
            if !state.renders.contains_key(&key) {
                state.renders.retain(|k, _| k.1 == state.version);
                match crate::render::part_png(&state.project, &state.eval, &part, 480, 320) {
                    Ok(png) => {
                        state.renders.insert(key.clone(), png);
                    }
                    Err(e) => return respond(&mut stream, "404 Not Found", "text/plain", e.as_bytes()),
                }
            }
            respond(&mut stream, "200 OK", "image/png", &state.renders[&key]);
        }
        ("POST", "/api/set") => {
            let body: serde_json::Value = serde_json::from_slice(&req.body).unwrap_or_default();
            let field = |k: &str| body.get(k).and_then(|v| v.as_str()).unwrap_or("").to_string();
            let (id, value, who) = (field("id"), field("value"), field("who"));
            let who = if who.trim().is_empty() { "someone".to_string() } else { who };
            let reply = match state.set(&id, value.trim(), &who) {
                Ok(summary) => serde_json::json!({ "ok": true, "version": state.version, "summary": summary }),
                Err(e) => serde_json::json!({ "ok": false, "error": e }),
            };
            respond(&mut stream, "200 OK", "application/json", &json(reply));
        }
        _ => respond(&mut stream, "404 Not Found", "text/plain", b"not found"),
    }
}

/// Serve `path` on an already-bound listener, in a background thread.
/// Returns the address it listens on.
pub fn start(path: &Path, listener: TcpListener) -> Result<SocketAddr, String> {
    let (project, eval) = State::load(path, Default::default())?;
    let addr = listener.local_addr().map_err(|e| e.to_string())?;
    let mut state = State {
        path: path.to_path_buf(),
        project,
        eval,
        version: 0,
        activity: vec![],
        mtime: mtime(path),
        renders: Default::default(),
    };
    std::thread::spawn(move || {
        for stream in listener.incoming().flatten() {
            handle(&mut state, stream);
        }
    });
    Ok(addr)
}

/// Serve until killed.
pub fn serve(path: &Path, addr: &str) -> Result<(), String> {
    let listener = TcpListener::bind(addr).map_err(|e| format!("{addr}: {e}"))?;
    let bound = start(path, listener)?;
    eprintln!("serving {} at http://{bound}/", path.display());
    if !bound.ip().is_loopback() {
        eprintln!("warning: there is no login; anyone who can reach {bound} can edit {}", path.display());
    }
    loop {
        std::thread::park();
    }
}
