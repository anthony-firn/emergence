//! Plain-text rendering of an evaluation: what each team needs to read.

use std::fmt::Write;

use emergence::Rigor;

use crate::graph::{Evaluation, Kind, Project, Status};
use crate::units::{parse_unit, sig, Dim};
use crate::value::Val;

/// Convert an SI value for display in `unit` (or plain SI).
fn shown(si: f64, dim: Dim, unit: Option<&str>) -> (f64, String) {
    match unit.map(|u| (u, parse_unit(u))) {
        Some((u, Ok(p))) if p.dim == dim => (p.from_si(si), if u == "1" { String::new() } else { u.to_string() }),
        _ => (si, dim.symbol()),
    }
}

/// `nominal unit [lo, hi]`.
pub fn val(v: &Val, unit: Option<&str>) -> String {
    let (n, u) = shown(v.nom, v.dim, unit);
    let space = if u.is_empty() { "" } else { " " };
    if v.is_exact() {
        return format!("{}{space}{u}", sig(n, 4));
    }
    let (lo, _) = shown(v.lo, v.dim, unit);
    let (hi, _) = shown(v.hi, v.dim, unit);
    format!("{}{space}{u} [{}, {}]", sig(n, 4), sig(lo, 4), sig(hi, 4))
}

/// A difference (no affine offset) in `unit`'s scale.
pub fn delta(d: f64, dim: Dim, unit: Option<&str>) -> String {
    match unit.map(|u| (u, parse_unit(u))) {
        Some((u, Ok(p))) if p.dim == dim => {
            let u = if p.is_affine() { "K" } else { u };
            format!("{} {u}", sig(d / p.scale, 3))
        }
        _ => format!("{} {}", sig(d, 3), dim.symbol()),
    }
}

pub fn rigor(r: Rigor) -> &'static str {
    r.name()
}

pub fn status_mark(s: Status) -> &'static str {
    match s {
        Status::Met => "✔",
        Status::Failed => "✘",
        Status::Undecided => "?",
        Status::Unavailable => "–",
    }
}

/// The full text report. `contributions` adds, per requirement with an
/// interval, which inputs (and owners) its uncertainty comes from.
pub fn text(p: &Project, e: &Evaluation, contributions: bool) -> String {
    let mut s = String::new();
    let _ = writeln!(s, "{}", p.doc.name);
    let _ = writeln!(s, "{}\n", "=".repeat(p.doc.name.chars().count()));
    let _ = writeln!(s, "Requirements");
    for c in e.checks.values() {
        let _ = write!(s, "  {} {:<10} {:<11} ", status_mark(c.status), c.id, c.status.name());
        match &c.value {
            Some(v) => {
                let limit = match (c.min, c.max) {
                    (Some(a), Some(b)) => format!("{} … {}", sig(shown(a, v.dim, c.unit.as_deref()).0, 4), sig(shown(b, v.dim, c.unit.as_deref()).0, 4)),
                    (None, Some(b)) => format!("≤ {}", sig(shown(b, v.dim, c.unit.as_deref()).0, 4)),
                    (Some(a), None) => format!("≥ {}", sig(shown(a, v.dim, c.unit.as_deref()).0, 4)),
                    _ => String::new(),
                };
                let margin = c.margin.map(|m| format!(", margin {}", delta(m, v.dim, c.unit.as_deref()))).unwrap_or_default();
                let _ = writeln!(s, "{} (limit {limit}{margin}; {})", val(v, c.unit.as_deref()), rigor(v.rigor));
            }
            None => {
                let _ = writeln!(s, "{}", c.reason.as_deref().unwrap_or("not evaluated"));
            }
        }
        let _ = writeln!(s, "      {} [{}]", c.text, c.owner);
        if contributions && c.value.is_some_and(|v| v.width() > 0.0) {
            let parts = p.contributions(&c.id, e);
            let shown: Vec<String> = parts
                .iter()
                .filter(|x| x.2 > 0.005)
                .take(4)
                .map(|(id, owner, f)| format!("{id} ({owner}) {:.0}%", 100.0 * f))
                .collect();
            if !shown.is_empty() {
                let _ = writeln!(s, "      uncertainty from: {}", shown.join(", "));
            }
        }
    }
    let mut owners: Vec<&str> = e.nodes.values().filter(|n| n.kind != Kind::Material).map(|n| n.owner.as_str()).collect();
    owners.sort();
    owners.dedup();
    for owner in owners {
        let _ = writeln!(s, "\n{owner}");
        for n in e.nodes.values().filter(|n| n.owner == owner && matches!(n.kind, Kind::Param | Kind::Derived | Kind::Output)) {
            let v = match &n.value {
                Ok(v) => format!("{}  {}", val(v, n.unit.as_deref()), if v.is_exact() { "" } else { rigor(v.rigor) }),
                Err(err) => format!("error: {err}"),
            };
            let tag = match n.kind {
                Kind::Param => "",
                _ => "=",
            };
            let _ = writeln!(s, "  {tag:1} {:<30} {v}", n.id);
        }
    }
    for l in &e.loops {
        let _ = writeln!(s, "\nFeedback loop {} (unknown {}):", l.nodes.join(", "), l.unknowns.join(", "));
        match &l.error {
            Some(err) => {
                let _ = writeln!(s, "  error: {err}");
            }
            None => {
                let _ = writeln!(
                    s,
                    "  settled in {} iterations; {}",
                    l.iterations,
                    if l.existence { "the interval maps into itself, so an operating point exists in it" } else { "existence not proven" }
                );
            }
        }
    }
    let _ = writeln!(s, "\nModels");
    for (id, r) in &e.runs {
        let _ = writeln!(s, "  {id} ({}, {}): {}", r.kind, r.owner, r.summary);
        if let Some(err) = &r.error {
            let _ = writeln!(s, "    error: {err}");
        }
        for n in &r.notes {
            let _ = writeln!(s, "    {n}");
        }
        if let Some(n) = &r.rigor_note {
            let _ = writeln!(s, "    {} tolerance corners; {n}", r.corners);
        } else if r.corners > 0 {
            let _ = writeln!(s, "    {} tolerance corners", r.corners);
        }
        for a in &r.assumptions {
            let _ = writeln!(s, "    assumes: {a}");
        }
    }
    let _ = writeln!(s, "\n{} model runs computed, {} reused", e.computed, e.cached);
    s
}
