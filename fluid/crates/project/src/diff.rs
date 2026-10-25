//! Semantic diff between two versions of a project.
//!
//! A text diff of the JSON says that `sink.fin_count` went from 9 to 6. This
//! says what that *did*: which results moved, which requirements changed
//! status, and which teams' changes caused it. It is the review a design
//! change needs before it is merged.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write;

use serde_json::Value as Json;

use crate::expr::Expr;
use crate::graph::{Evaluation, Kind, Project, Status};
use crate::report;
use crate::value::Val;

#[derive(Clone, Debug)]
pub struct Change {
    pub id: String,
    pub owner: String,
    pub what: String,
}

#[derive(Clone, Debug)]
pub struct Effect {
    pub id: String,
    pub owner: String,
    pub before: Option<Val>,
    pub after: Option<Val>,
    pub unit: Option<String>,
}

#[derive(Clone, Debug)]
pub struct Verdict {
    pub id: String,
    pub owner: String,
    pub text: String,
    pub before: Status,
    pub after: Status,
    pub before_value: Option<Val>,
    pub after_value: Option<Val>,
    pub unit: Option<String>,
    /// Why the new version cannot be judged, if it cannot.
    pub reason: Option<String>,
    /// Changed ids upstream of this requirement, with their owners.
    pub causes: Vec<(String, String)>,
}

#[derive(Clone, Debug, Default)]
pub struct Diff {
    pub changes: Vec<Change>,
    pub effects: Vec<Effect>,
    /// Every requirement whose status or value moved.
    pub verdicts: Vec<Verdict>,
}

fn describe(section: &str, a: Option<&Json>, b: Option<&Json>) -> Option<String> {
    let short = |v: &Json| -> String {
        match v {
            Json::Object(o) => {
                if let Some(Json::String(s)) = o.get("value").or_else(|| o.get("expr")) {
                    let tol = o.get("tol").map(|t| format!(" ± {}", t.to_string().trim_matches('"'))).unwrap_or_default();
                    return format!("{s}{tol}");
                }
                v.to_string()
            }
            _ => v.to_string(),
        }
    };
    match (a, b) {
        (None, Some(b)) => Some(format!("added to {section}: {}", short(b))),
        (Some(a), None) => Some(format!("removed from {section} (was {})", short(a))),
        (Some(a), Some(b)) if a != b => {
            let (oa, ob) = (a.as_object(), b.as_object());
            if let (Some(oa), Some(ob)) = (oa, ob) {
                let keys: BTreeSet<&String> = oa.keys().chain(ob.keys()).collect();
                let parts: Vec<String> = keys
                    .into_iter()
                    .filter(|k| oa.get(*k) != ob.get(*k))
                    .map(|k| {
                        let f = |v: Option<&Json>| v.map(|v| v.to_string()).unwrap_or_else(|| "—".into());
                        format!("{k} {} → {}", f(oa.get(k)), f(ob.get(k)))
                    })
                    .collect();
                Some(parts.join("; "))
            } else {
                Some(format!("{} → {}", short(a), short(b)))
            }
        }
        _ => None,
    }
}

/// Definition changes, section by section, with owners.
pub fn changes(old: &Project, new: &Project) -> Vec<Change> {
    let a = serde_json::to_value(&old.doc).unwrap();
    let b = serde_json::to_value(&new.doc).unwrap();
    let mut out = Vec::new();
    for section in ["materials", "params", "derived", "analyses", "parts", "requirements"] {
        let empty = serde_json::Map::new();
        let sa = a.get(section).and_then(|v| v.as_object()).unwrap_or(&empty);
        let sb = b.get(section).and_then(|v| v.as_object()).unwrap_or(&empty);
        let ids: BTreeSet<&String> = sa.keys().chain(sb.keys()).collect();
        for id in ids {
            if let Some(what) = describe(section, sa.get(id), sb.get(id)) {
                let key = if section == "materials" { format!("mat.{id}") } else { id.clone() };
                let owner = new.doc.owner(&key).or_else(|| old.doc.owner(&key)).unwrap_or("?").to_string();
                out.push(Change { id: key, owner, what });
            }
        }
    }
    // A part's CAD model can change without the project file changing.
    for (id, p) in &new.doc.parts {
        if let Some(q) = old.doc.parts.get(id) {
            let read = |proj: &Project, f: &str| std::fs::read_to_string(proj.dir.join(f)).ok();
            if q.model == p.model && read(old, &q.model) != read(new, &p.model) {
                out.push(Change { id: id.clone(), owner: p.owner.clone(), what: format!("CAD model {} edited", p.model) });
            }
        }
    }
    out
}

/// Every id upstream of a requirement (outputs bring in their analysis).
fn upstream(p: &Project, e: &Evaluation, req: &str) -> BTreeSet<String> {
    let mut seen = BTreeSet::new();
    let Some(r) = p.doc.requirements.get(req) else { return seen };
    let mut stack: Vec<String> = Expr::parse(&r.of).map(|e| e.deps()).unwrap_or_default();
    seen.insert(req.to_string());
    while let Some(id) = stack.pop() {
        if !seen.insert(id.clone()) {
            continue;
        }
        if let Some((head, _)) = id.rsplit_once('.') {
            if p.model(head).is_some() {
                seen.insert(head.to_string());
            }
            if let Some(m) = id.strip_prefix("mat.") {
                seen.insert(format!("mat.{}", m.split('.').next().unwrap_or("")));
            }
        }
        if let Some(n) = e.nodes.get(&id) {
            stack.extend(n.deps.iter().cloned());
        }
    }
    seen
}

fn moved(a: Option<Val>, b: Option<Val>) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => {
            let scale = a.nom.abs().max(b.nom.abs()).max(1e-300);
            ((a.nom - b.nom).abs() + (a.lo - b.lo).abs() + (a.hi - b.hi).abs()) > 1e-12 * scale
        }
        (None, None) => false,
        _ => true,
    }
}

pub fn diff(old: &Project, new: &Project) -> Diff {
    let (ea, eb) = (old.evaluate(), new.evaluate());
    let changes = changes(old, new);
    let changed: BTreeMap<&str, &str> = changes.iter().map(|c| (c.id.as_str(), c.owner.as_str())).collect();
    let mut effects = Vec::new();
    let ids: BTreeSet<&String> = ea.nodes.keys().chain(eb.nodes.keys()).collect();
    for id in ids {
        let (na, nb) = (ea.nodes.get(id), eb.nodes.get(id));
        let n = nb.or(na).unwrap();
        if !matches!(n.kind, Kind::Derived | Kind::Output) {
            continue;
        }
        let (va, vb) = (na.and_then(|n| n.value.clone().ok()), nb.and_then(|n| n.value.clone().ok()));
        if moved(va, vb) {
            effects.push(Effect { id: id.clone(), owner: n.owner.clone(), before: va, after: vb, unit: n.unit.clone() });
        }
    }
    let mut verdicts = Vec::new();
    let reqs: BTreeSet<&String> = ea.checks.keys().chain(eb.checks.keys()).collect();
    for id in reqs {
        let (ca, cb) = (ea.checks.get(id), eb.checks.get(id));
        let c = cb.or(ca).unwrap();
        let before = ca.map(|c| c.status).unwrap_or(Status::Unavailable);
        let after = cb.map(|c| c.status).unwrap_or(Status::Unavailable);
        let (va, vb) = (ca.and_then(|c| c.value), cb.and_then(|c| c.value));
        if before == after && !moved(va, vb) {
            continue;
        }
        let up = upstream(new, &eb, id);
        let causes = changed
            .iter()
            .filter(|(cid, _)| up.contains(**cid) || **cid == id.as_str())
            .map(|(cid, o)| (cid.to_string(), o.to_string()))
            .collect();
        verdicts.push(Verdict {
            id: id.clone(),
            owner: c.owner.clone(),
            text: c.text.clone(),
            before,
            after,
            before_value: va,
            after_value: vb,
            unit: c.unit.clone(),
            reason: cb.and_then(|c| c.reason.clone()),
            causes,
        });
    }
    Diff { changes, effects, verdicts }
}

pub fn text(d: &Diff) -> String {
    let mut s = String::new();
    if d.changes.is_empty() {
        let _ = writeln!(s, "No changes.");
        return s;
    }
    let _ = writeln!(s, "Changes");
    for c in &d.changes {
        let _ = writeln!(s, "  {:<24} [{}] {}", c.id, c.owner, c.what);
    }
    let flips: Vec<&Verdict> = d.verdicts.iter().filter(|v| v.before != v.after).collect();
    if !flips.is_empty() {
        let _ = writeln!(s, "\nRequirements that changed status");
        for v in &flips {
            let _ = writeln!(
                s,
                "  {} {:<10} {} → {}   [{}] {}",
                report::status_mark(v.after),
                v.id,
                v.before.name(),
                v.after.name(),
                v.owner,
                v.text
            );
            write_values(&mut s, v);
        }
    }
    let shifted: Vec<&Verdict> = d.verdicts.iter().filter(|v| v.before == v.after).collect();
    if !shifted.is_empty() {
        let _ = writeln!(s, "\nRequirements whose value moved (status unchanged)");
        for v in &shifted {
            let _ = writeln!(s, "  {} {:<10} {}   [{}]", report::status_mark(v.after), v.id, v.after.name(), v.owner);
            write_values(&mut s, v);
        }
    }
    if !d.effects.is_empty() {
        let _ = writeln!(s, "\nResults that moved");
        for e in &d.effects {
            let f = |v: &Option<Val>| v.map(|v| report::val(&v, e.unit.as_deref())).unwrap_or_else(|| "error".into());
            let _ = writeln!(s, "  {:<28} [{}] {} → {}", e.id, e.owner, f(&e.before), f(&e.after));
        }
    }
    s
}

fn write_values(s: &mut String, v: &Verdict) {
    let f = |x: &Option<Val>| x.map(|x| report::val(&x, v.unit.as_deref())).unwrap_or_else(|| "unavailable".into());
    let _ = writeln!(s, "      {} → {}", f(&v.before_value), f(&v.after_value));
    if let (None, Some(r)) = (v.after_value, &v.reason) {
        let _ = writeln!(s, "      why: {r}");
    }
    if !v.causes.is_empty() {
        let c: Vec<String> = v.causes.iter().map(|(id, o)| format!("{id} ({o})")).collect();
        let _ = writeln!(s, "      caused by: {}", c.join(", "));
    }
}
