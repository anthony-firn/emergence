//! Evaluating a project: the dependency graph.
//!
//! Every parameter, derived value, analysis, part and requirement is a node;
//! edges come from the expressions that reference them. The graph is sorted
//! once (a cycle is reported with its path), then evaluated in order with
//! intervals. Analyses and parts are evaluated at the corners of their
//! toleranced inputs ([`crate::analysis`] explains when that is rigorous).
//! Results are memoised by exact input bits, so after an edit only what the
//! edit reaches runs again, and the answer is identical to a full rerun.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use emergence::Rigor;

use crate::analysis::{self, Analysis, Outcome, Port};
use crate::doc::{Doc, Tol};
use crate::expr::{Expr, Scope};
use crate::materials::{Library, Property};
use crate::units::{parse_quantity, parse_unit, split_quantity, Dim};
use crate::value::Val;

/// Memoised analysis runs, shareable between versions of a project.
pub type Cache = Arc<Mutex<HashMap<(String, Vec<u64>), Result<Outcome, String>>>>;

/// Corner enumeration is exhaustive up to this many toleranced inputs.
const MAX_CORNER_INPUTS: usize = 10;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Kind {
    Param,
    Material,
    Derived,
    Output,
    Requirement,
}

#[derive(Clone, Debug)]
pub struct Node {
    pub id: String,
    pub kind: Kind,
    pub owner: String,
    pub label: String,
    /// Display unit.
    pub unit: Option<String>,
    /// Node ids read directly (for outputs: the analysis's inputs).
    pub deps: Vec<String>,
    pub value: Result<Val, String>,
}

/// How an analysis or part ran.
#[derive(Clone, Debug)]
pub struct Run {
    pub kind: String,
    pub summary: String,
    pub owner: String,
    pub label: String,
    /// (input, binding expression, value)
    pub inputs: Vec<(String, String, Option<Val>)>,
    pub corners: usize,
    pub rigor_note: Option<String>,
    pub notes: Vec<String>,
    pub assumptions: Vec<String>,
    pub error: Option<String>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    Met,
    Failed,
    Undecided,
    Unavailable,
}

impl Status {
    pub fn name(self) -> &'static str {
        match self {
            Status::Met => "met",
            Status::Failed => "FAILED",
            Status::Undecided => "undecided",
            Status::Unavailable => "unavailable",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Check {
    pub id: String,
    pub text: String,
    pub owner: String,
    pub of: String,
    pub status: Status,
    pub value: Option<Val>,
    pub max: Option<f64>,
    pub min: Option<f64>,
    /// Worst-case distance inside the limits (negative: outside), SI.
    pub margin: Option<f64>,
    pub unit: Option<String>,
    pub reason: Option<String>,
}

#[derive(Clone, Debug)]
pub struct Evaluation {
    pub nodes: BTreeMap<String, Node>,
    pub runs: BTreeMap<String, Run>,
    pub checks: BTreeMap<String, Check>,
    /// Model runs actually computed / served from the cache.
    pub computed: usize,
    pub cached: usize,
    /// Feedback loops solved as fixed points.
    pub loops: Vec<Loop>,
}

/// A feedback loop and how it was solved.
#[derive(Clone, Debug)]
pub struct Loop {
    /// Every node on the loop.
    pub nodes: Vec<String>,
    /// The unknowns (nodes with a `range`).
    pub unknowns: Vec<String>,
    pub iterations: usize,
    /// `F(X) ⊆ X` held at the solution, so a fixed point provably exists
    /// in the reported interval (Brouwer).
    pub existence: bool,
    pub error: Option<String>,
}

impl Evaluation {
    pub fn value(&self, id: &str) -> Option<Val> {
        self.nodes.get(id).and_then(|n| n.value.as_ref().ok().copied())
    }
}

/// A loaded project, ready to evaluate.
pub struct Project {
    pub doc: Doc,
    pub dir: PathBuf,
    pub lib: Library,
    models: BTreeMap<String, Box<dyn Analysis>>,
    /// Cache key prefix per model (identity of its definition).
    model_keys: BTreeMap<String, String>,
    cache: Cache,
}

fn valid_id(id: &str) -> bool {
    let mut c = id.chars();
    c.next().is_some_and(|f| f.is_ascii_alphabetic())
        && id.chars().all(|ch| ch.is_ascii_alphanumeric() || ch == '_' || ch == '.' || ch == '-')
        && !id.ends_with('.')
        && !id.starts_with("mat.")
}

/// Parse `"20%"` as a relative tolerance.
fn percent(text: &str) -> Option<f64> {
    text.trim().strip_suffix('%').and_then(|n| n.trim().parse::<f64>().ok()).map(|p| p / 100.0)
}

impl Project {
    pub fn load(path: impl AsRef<Path>) -> Result<Project, String> {
        let path = path.as_ref();
        let doc = Doc::read(path)?;
        Project::new(doc, path.parent().unwrap_or(Path::new(".")), Cache::default())
    }

    pub fn new(doc: Doc, dir: impl AsRef<Path>, cache: Cache) -> Result<Project, String> {
        let dir = dir.as_ref().to_path_buf();
        let mut lib = Library::builtin();
        for (m, def) in &doc.materials {
            for (p, prop) in &def.props {
                let rel = match &prop.tol {
                    None => 0.0,
                    Some(t) => percent(t).ok_or_else(|| format!("material {m}.{p}: tol must be relative, e.g. \"10%\""))?,
                };
                lib.define(m, p, &prop.value, rel, &prop.source)?;
            }
            if !def.description.is_empty() {
                lib.materials.get_mut(m).unwrap().description = def.description.clone();
            }
        }
        // Ids: valid, and unique across every section.
        let mut seen = BTreeSet::new();
        let sections: [(&str, Vec<&String>); 5] = [
            ("params", doc.params.keys().collect()),
            ("derived", doc.derived.keys().collect()),
            ("analyses", doc.analyses.keys().collect()),
            ("parts", doc.parts.keys().collect()),
            ("requirements", doc.requirements.keys().collect()),
        ];
        for (section, ids) in &sections {
            for id in ids {
                if !valid_id(id) {
                    return Err(format!("{section}: `{id}` is not a valid id (letters, digits, _ . -; not starting with mat.)"));
                }
                if !seen.insert(id.to_string()) {
                    return Err(format!("{section}: id `{id}` is used twice"));
                }
            }
        }
        let mut models: BTreeMap<String, Box<dyn Analysis>> = BTreeMap::new();
        let mut model_keys = BTreeMap::new();
        for (id, a) in &doc.analyses {
            let (model, key) = analysis::build(&a.kind, &a.config, &dir).map_err(|e| format!("analysis `{id}`: {e}"))?;
            model_keys.insert(id.clone(), key);
            models.insert(id.clone(), model);
        }
        for (id, p) in &doc.parts {
            for key in p.fdm.iter().flat_map(|f| f.keys()) {
                if !["nozzle", "perimeters", "infill", "max_flow"].contains(&key.as_str()) {
                    return Err(format!("part `{id}`: unknown fdm setting `{key}` (nozzle, perimeters, infill, max_flow)"));
                }
                if p.params.contains_key(key) {
                    return Err(format!("part `{id}`: `{key}` is both a CAD parameter and an fdm setting"));
                }
            }
            let part = build_part(&dir, id, p)?;
            model_keys.insert(
                id.clone(),
                format!("part:{}:{}:{}", dir.join(&p.model).display(), part.params.join(","), part.cell),
            );
            models.insert(id.clone(), Box::new(part));
        }
        // Every node's own fields are checked here; references are checked
        // when the graph is built, so one bad node does not hide the others.
        for (id, a) in &doc.analyses {
            let ports = models[id].inputs();
            for name in a.inputs.keys() {
                if !ports.iter().any(|p| &p.name == name) {
                    return Err(format!(
                        "analysis `{id}`: `{}` has no input `{name}` (inputs: {})",
                        a.kind,
                        ports.iter().map(|p| p.name.as_str()).collect::<Vec<_>>().join(", ")
                    ));
                }
            }
            for slot in a.materials.keys() {
                if !ports.iter().any(|p| p.from.as_ref().is_some_and(|f| &f.slot == slot)) {
                    return Err(format!("analysis `{id}`: `{}` has no material slot `{slot}`", a.kind));
                }
            }
        }
        Ok(Project { doc, dir, lib, models, model_keys, cache })
    }

    /// The CAD part `id`, rebuilt from its definition.
    pub fn part(&self, id: &str) -> Option<analysis::Part> {
        self.doc.parts.get(id).and_then(|p| build_part(&self.dir, id, p).ok())
    }

    pub fn cache(&self) -> Cache {
        self.cache.clone()
    }

    pub fn model(&self, id: &str) -> Option<&dyn Analysis> {
        self.models.get(id).map(|m| m.as_ref())
    }

    /// The input bindings of analysis or part `id`: (port, expression text).
    pub fn bindings(&self, id: &str) -> Result<Vec<(Port, String)>, String> {
        let model = &self.models[id];
        let (given, slots): (BTreeMap<String, String>, BTreeMap<String, String>) =
            if let Some(a) = self.doc.analyses.get(id) {
                (a.inputs.clone(), a.materials.clone())
            } else {
                let p = &self.doc.parts[id];
                let mut given = p.params.clone();
                if let Some(f) = &p.fdm {
                    for (k, v) in f {
                        given.insert(k.clone(), v.clone());
                    }
                }
                (given, BTreeMap::from([("material".to_string(), p.material.clone())]))
            };
        let ports = model.inputs();
        let mut out = Vec::new();
        for port in &ports {
            let text = if let Some(t) = given.get(&port.name) {
                t.clone()
            } else if let Some(t) = port.default.as_ref().and_then(|d| given.get(d)) {
                t.clone()
            } else if let Some(f) = &port.from {
                let m = slots.get(&f.slot).ok_or_else(|| {
                    format!("`{id}`: input `{}` is not bound and material slot `{}` is not set", port.name, f.slot)
                })?;
                match &f.at {
                    None => format!("mat.{m}.{}", f.prop),
                    Some(at) => {
                        let t = given.get(at).ok_or_else(|| format!("`{id}`: input `{at}` is not bound"))?;
                        format!("mat.{m}.{}({t})", f.prop)
                    }
                }
            } else {
                return Err(format!("`{id}`: input `{}` ({}) is not bound", port.name, port.doc));
            };
            out.push((port.clone(), text));
        }
        Ok(out)
    }

    /// Map a referenced id to the graph node that produces it.
    fn producer(&self, id: &str) -> Result<Option<String>, String> {
        if Library::split(id).is_some() {
            return Ok(None);
        }
        if self.doc.params.contains_key(id) || self.doc.derived.contains_key(id) {
            return Ok(Some(id.to_string()));
        }
        if self.models.contains_key(id) {
            let outs: Vec<String> = self.models[id].outputs().iter().map(|o| format!("{id}.{}", o.name)).collect();
            return Err(format!("`{id}` is an analysis; refer to one of its outputs: {}", outs.join(", ")));
        }
        if let Some((head, tail)) = id.rsplit_once('.') {
            if let Some(m) = self.models.get(head) {
                if m.outputs().iter().any(|o| o.name == tail) {
                    return Ok(Some(head.to_string()));
                }
                return Err(format!(
                    "`{head}` has no output `{tail}` (outputs: {})",
                    m.outputs().iter().map(|o| o.name.as_str()).collect::<Vec<_>>().join(", ")
                ));
            }
        }
        Err(format!("unknown id `{id}`"))
    }

    /// Direct dependencies (as referenced ids) of every graph node.
    fn references(&self) -> BTreeMap<String, Result<Vec<String>, String>> {
        let mut out = BTreeMap::new();
        for id in self.doc.params.keys() {
            out.insert(id.clone(), Ok(vec![]));
        }
        for (id, d) in &self.doc.derived {
            out.insert(id.clone(), Expr::parse(&d.expr).map(|e| e.deps()));
        }
        for id in self.models.keys() {
            let deps = self.bindings(id).and_then(|b| {
                let mut all = Vec::new();
                for (_, text) in b {
                    for d in Expr::parse(&text)?.deps() {
                        if !all.contains(&d) {
                            all.push(d);
                        }
                    }
                }
                Ok(all)
            });
            out.insert(id.clone(), deps);
        }
        out
    }

    /// Evaluate everything. `pins` names leaves (params, `mat.<m>.<p>`)
    /// whose uncertainty is removed, which is how contributions are measured.
    pub fn evaluate_pinned(&self, pins: &BTreeSet<String>) -> Evaluation {
        let refs = self.references();
        // Graph edges between producer nodes.
        let mut edges: BTreeMap<String, Result<Vec<String>, String>> = BTreeMap::new();
        for (id, r) in &refs {
            let e = r.clone().and_then(|ids| {
                let mut p = Vec::new();
                for d in ids {
                    if let Some(n) = self.producer(&d)? {
                        if !p.contains(&n) {
                            p.push(n);
                        }
                    }
                }
                Ok(p)
            });
            edges.insert(id.clone(), e);
        }
        let mut values: BTreeMap<String, Result<Val, String>> = BTreeMap::new();
        let mut nodes = BTreeMap::new();
        let mut runs = BTreeMap::new();
        let mut stats = (0usize, 0usize);
        let mut loops = Vec::new();
        for scc in components(&edges) {
            let looped = scc.len() > 1 || edges[&scc[0]].as_ref().is_ok_and(|d| d.contains(&scc[0]));
            if !looped {
                self.eval_unit(&scc[0], &edges, &mut values, &mut runs, pins, &mut stats);
            } else {
                loops.push(self.solve_loop(&scc, &edges, &mut values, &mut runs, pins, &mut stats));
            }
        }
        let (computed, cached) = stats;
        // Nodes for the report.
        for (id, p) in &self.doc.params {
            let unit = p.unit.clone().or_else(|| split_quantity(&p.value).ok().map(|(_, u)| u.to_string()).filter(|u| !u.is_empty()));
            nodes.insert(id.clone(), node(id, Kind::Param, &p.owner, &p.label, unit, vec![], &values));
        }
        for (id, d) in &self.doc.derived {
            let deps = refs[id].clone().unwrap_or_default();
            nodes.insert(id.clone(), node(id, Kind::Derived, &d.owner, &d.label, d.unit.clone(), deps, &values));
        }
        for (id, m) in &self.models {
            let (owner, label) = match self.doc.analyses.get(id) {
                Some(a) => (a.owner.clone(), a.label.clone()),
                None => (self.doc.parts[id].owner.clone(), self.doc.parts[id].label.clone()),
            };
            for o in m.outputs() {
                let oid = format!("{id}.{}", o.name);
                let deps = refs[id].clone().unwrap_or_default();
                let label = if label.is_empty() { o.doc.clone() } else { format!("{label}: {}", o.doc) };
                let mut n = node(&oid, Kind::Output, &owner, &label, Some(o.unit.clone()), deps, &values);
                if !values.contains_key(&oid) {
                    n.value = Err(values.get(id).cloned().and_then(|v| v.err()).unwrap_or_else(|| "not evaluated".into()));
                }
                nodes.insert(oid, n);
            }
        }
        // Materials that were referenced.
        for r in refs.values().flatten() {
            for d in r {
                if let Some((m, p)) = Library::split(d) {
                    let v = self.lib.get(m, p).map(|prop| pinned_prop(prop, pins.contains(d)).reference());
                    let owner = self.doc.owner(d).unwrap_or("materials library").to_string();
                    let label = self.lib.get(m, p).map(|prop| prop.source.clone()).unwrap_or_default();
                    nodes.insert(
                        d.clone(),
                        Node { id: d.clone(), kind: Kind::Material, owner, label, unit: None, deps: vec![], value: v },
                    );
                }
            }
        }
        let checks = self.doc.requirements.keys().map(|id| (id.clone(), self.check(id, &values, pins))).collect();
        for (id, r) in &self.doc.requirements {
            let deps = Expr::parse(&r.of).map(|e| e.deps()).unwrap_or_default();
            let mut n = node(id, Kind::Requirement, &r.owner, &r.text, None, deps, &values);
            n.value = Err("requirement".into());
            nodes.insert(id.clone(), n);
        }
        Evaluation { nodes, runs, checks, computed, cached, loops }
    }

    /// Evaluate one acyclic node (a param, a derived value, or a model and
    /// its outputs) into `values`.
    fn eval_unit(
        &self,
        id: &str,
        edges: &BTreeMap<String, Result<Vec<String>, String>>,
        values: &mut BTreeMap<String, Result<Val, String>>,
        runs: &mut BTreeMap<String, Run>,
        pins: &BTreeSet<String>,
        stats: &mut (usize, usize),
    ) {
        if let Some(p) = self.doc.params.get(id) {
            let v = param_value(&p.value, p.tol.as_ref()).map(|v| if pins.contains(id) { v.pinned() } else { v });
            values.insert(id.to_string(), v);
        } else if let Some(d) = self.doc.derived.get(id) {
            let scope = EvalScope { values, lib: &self.lib, pins };
            let v = match (&edges[id], Expr::parse(&d.expr)) {
                (Err(e), _) => Err(e.clone()),
                (_, Err(e)) => Err(e),
                (Ok(_), Ok(e)) => e.eval(&scope).and_then(|v| check_unit(v, d.unit.as_deref())),
            };
            values.insert(id.to_string(), v);
        } else if self.models.contains_key(id) {
            let (run, outs, c, h) = self.run_model(id, &edges[id], values, pins);
            stats.0 += c;
            stats.1 += h;
            for (name, v) in outs {
                values.insert(format!("{id}.{name}"), v);
            }
            runs.insert(id.to_string(), run);
        }
    }

    /// Solve a strongly connected component as a fixed point `X = F(X)` in
    /// the unknowns' ranges, by `X ← F(X) ∩ X`. Every fixed point inside the
    /// range survives each step (it equals its own image), so an empty
    /// intersection proves there is none there (thermal runaway, for
    /// example). When the final box maps into itself, one exists (Brouwer).
    /// The nominal is iterated from the bottom of the range, which for an
    /// increasing loop converges to the lowest operating point.
    fn solve_loop(
        &self,
        scc: &[String],
        edges: &BTreeMap<String, Result<Vec<String>, String>>,
        values: &mut BTreeMap<String, Result<Val, String>>,
        runs: &mut BTreeMap<String, Run>,
        pins: &BTreeSet<String>,
        stats: &mut (usize, usize),
    ) -> Loop {
        let mut info = Loop { nodes: scc.to_vec(), unknowns: vec![], iterations: 0, existence: false, error: None };
        let cycle = format!("{} → {}", scc.join(" → "), scc[0]);
        let fail = |values: &mut BTreeMap<String, Result<Val, String>>, info: &mut Loop, e: String| {
            for id in scc {
                values.insert(id.clone(), Err(e.clone()));
                if let Some(m) = self.models.get(id) {
                    for o in m.outputs() {
                        values.insert(format!("{id}.{}", o.name), Err(e.clone()));
                    }
                }
            }
            info.error = Some(e);
        };
        let unknowns: Vec<String> = scc.iter().filter(|id| self.doc.derived.get(*id).is_some_and(|d| d.range.is_some())).cloned().collect();
        info.unknowns = unknowns.clone();
        if unknowns.is_empty() {
            fail(values, &mut info, format!("dependency cycle: {cycle}; give one of its derived values a `range` to solve it as a feedback loop"));
            return info;
        }
        // The rest of the loop, in order, with the edges into the unknowns cut.
        let inner: BTreeMap<String, Result<Vec<String>, String>> = scc
            .iter()
            .filter(|id| !unknowns.contains(id))
            .map(|id| {
                let d = edges[id].clone().map(|d| d.into_iter().filter(|x| scc.contains(x) && !unknowns.contains(x)).collect());
                (id.clone(), d)
            })
            .collect();
        let (order, cyclic) = topo(&inner);
        if let Some(c) = cyclic.values().next() {
            fail(values, &mut info, format!("feedback loop {c} has no node with a `range`"));
            return info;
        }
        let mut x: BTreeMap<String, Val> = BTreeMap::new();
        for u in &unknowns {
            let d = &self.doc.derived[u];
            let [a, b] = d.range.as_ref().unwrap();
            let r = (|| -> Result<Val, String> {
                let (qa, _) = parse_quantity(a)?;
                let (qb, _) = parse_quantity(b)?;
                if qa.dim != qb.dim {
                    return Err(format!("`{u}`: range ends differ in dimension"));
                }
                let (lo, hi) = (qa.si.min(qb.si), qa.si.max(qb.si));
                check_unit(Val::new(lo, lo, hi, qa.dim, Rigor::Rigorous)?, d.unit.as_deref())
            })();
            match r {
                Ok(v) => {
                    x.insert(u.clone(), v);
                }
                Err(e) => {
                    fail(values, &mut info, e);
                    return info;
                }
            }
        }
        const MAX_ITERATIONS: usize = 400;
        let image = |x: &BTreeMap<String, Val>,
                     values: &mut BTreeMap<String, Result<Val, String>>,
                     runs: &mut BTreeMap<String, Run>,
                     stats: &mut (usize, usize)|
         -> Result<BTreeMap<String, Val>, String> {
            for (u, v) in x {
                values.insert(u.clone(), Ok(*v));
            }
            for id in &order {
                self.eval_unit(id, edges, values, runs, pins, stats);
            }
            let mut fx = BTreeMap::new();
            for u in x.keys() {
                let d = &self.doc.derived[u];
                let scope = EvalScope { values, lib: &self.lib, pins };
                let v = Expr::parse(&d.expr).and_then(|e| e.eval(&scope)).and_then(|v| check_unit(v, d.unit.as_deref()))?;
                fx.insert(u.clone(), v);
            }
            Ok(fx)
        };
        for it in 1..=MAX_ITERATIONS {
            info.iterations = it;
            let fx = match image(&x, values, runs, stats) {
                Ok(f) => f,
                Err(e) => {
                    fail(values, &mut info, format!("{e} (in feedback loop {cycle})"));
                    return info;
                }
            };
            let mut next = BTreeMap::new();
            let mut moved = false;
            for (u, xv) in &x {
                let f = fx[u];
                let (lo, hi) = (xv.lo.max(f.lo), xv.hi.min(f.hi));
                if lo > hi {
                    let unit = self.doc.derived[u].unit.as_deref();
                    let e = format!(
                        "no operating point for `{u}` in {}: the loop's image {} lies outside it (runaway, or the range is too small)",
                        crate::report::val(xv, unit),
                        crate::report::val(&f, unit)
                    );
                    fail(values, &mut info, e);
                    return info;
                }
                let nom = f.nom.clamp(lo, hi);
                let v = Val { nom, lo, hi, dim: xv.dim, rigor: f.rigor };
                let scale = xv.hi.abs().max(xv.lo.abs()).max(1e-300);
                moved |= (v.lo - xv.lo).abs() > 1e-13 * scale || (v.hi - xv.hi).abs() > 1e-13 * scale || (v.nom - xv.nom).abs() > 1e-13 * scale;
                next.insert(u.clone(), v);
            }
            x = next;
            if !moved {
                // Final pass on the converged box, and the existence check.
                match image(&x, values, runs, stats) {
                    Ok(fx) => {
                        info.existence = x.iter().all(|(u, v)| fx[u].lo >= v.lo && fx[u].hi <= v.hi);
                        for (u, v) in &x {
                            let rigor = if info.existence { v.rigor } else { v.rigor.min(Rigor::Asymptotic) };
                            values.insert(u.clone(), Ok(Val { rigor, ..*v }));
                        }
                    }
                    Err(e) => fail(values, &mut info, e),
                }
                return info;
            }
        }
        fail(values, &mut info, format!("feedback loop {cycle} did not settle in {MAX_ITERATIONS} iterations (loop gain near or above 1?)"));
        info
    }

    pub fn evaluate(&self) -> Evaluation {
        self.evaluate_pinned(&BTreeSet::new())
    }

    fn run_model(
        &self,
        id: &str,
        edges: &Result<Vec<String>, String>,
        values: &BTreeMap<String, Result<Val, String>>,
        pins: &BTreeSet<String>,
    ) -> (Run, Vec<(String, Result<Val, String>)>, usize, usize) {
        let model = &self.models[id];
        let (owner, label) = match self.doc.analyses.get(id) {
            Some(a) => (a.owner.clone(), a.label.clone()),
            None => (self.doc.parts[id].owner.clone(), self.doc.parts[id].label.clone()),
        };
        let mut run = Run {
            kind: model.kind().into(),
            summary: model.summary().into(),
            owner,
            label,
            inputs: vec![],
            corners: 0,
            rigor_note: None,
            notes: vec![],
            assumptions: model.assumptions(),
            error: None,
        };
        let outputs = model.outputs();
        let fail = |run: &mut Run, e: String| {
            run.error = Some(e.clone());
            outputs.iter().map(|o| (o.name.clone(), Err(e.clone()))).collect::<Vec<_>>()
        };
        if let Err(e) = edges {
            let o = fail(&mut run, e.clone());
            return (run, o, 0, 0);
        }
        let bindings = match self.bindings(id) {
            Ok(b) => b,
            Err(e) => {
                let o = fail(&mut run, e);
                return (run, o, 0, 0);
            }
        };
        let scope = EvalScope { values, lib: &self.lib, pins };
        let mut xs = Vec::new();
        let mut bad = None;
        for (port, text) in &bindings {
            let v = Expr::parse(text).and_then(|e| e.eval(&scope)).and_then(|v| {
                match port.dim {
                    Some(d) if d != v.dim => {
                        Err(format!("input `{}` needs {} ({}), the binding gives {}", port.name, d.symbol(), port.unit, v.dim.symbol()))
                    }
                    None if v.dim != Dim::LENGTH && !v.dim.is_none() => {
                        Err(format!("CAD parameter `{}` must be a length or a number, got {}", port.name, v.dim.symbol()))
                    }
                    _ => Ok(v),
                }
            });
            run.inputs.push((port.name.clone(), text.clone(), v.as_ref().ok().copied()));
            match v {
                Ok(v) => xs.push(v),
                Err(e) => {
                    bad.get_or_insert(format!("input `{}`: {e}", port.name));
                }
            }
        }
        if let Some(e) = bad {
            let o = fail(&mut run, e);
            return (run, o, 0, 0);
        }
        let ports: Vec<&Port> = bindings.iter().map(|b| &b.0).collect();
        let nominal: Vec<f64> = xs.iter().map(|v| v.nom).collect();
        let toleranced: Vec<usize> = (0..xs.len()).filter(|&i| xs[i].width() > 0.0 && !ports[i].passive).collect();
        let mut points = vec![nominal.clone()];
        let exhaustive = toleranced.len() <= MAX_CORNER_INPUTS;
        if exhaustive {
            for mask in 0..(1usize << toleranced.len()) {
                let mut x = nominal.clone();
                for (b, &i) in toleranced.iter().enumerate() {
                    x[i] = if mask >> b & 1 == 1 { xs[i].hi } else { xs[i].lo };
                }
                points.push(x);
            }
        } else {
            for &i in &toleranced {
                for end in [xs[i].lo, xs[i].hi] {
                    let mut x = nominal.clone();
                    x[i] = end;
                    points.push(x);
                }
            }
        }
        run.corners = points.len() - 1;
        let key = &self.model_keys[id];
        // Corners are independent runs: evaluate them in parallel, each
        // through the cache.
        let outcomes: Vec<(Result<Outcome, String>, bool)> = {
            use rayon::prelude::*;
            points
                .par_iter()
                .map(|x| {
                    let k = (key.clone(), x.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
                    let hit = self.cache.lock().unwrap().get(&k).cloned();
                    match hit {
                        Some(r) => (r, true),
                        None => {
                            let r = model.run(x);
                            self.cache.lock().unwrap().insert(k, r.clone());
                            (r, false)
                        }
                    }
                })
                .collect()
        };
        let hits = outcomes.iter().filter(|o| o.1).count();
        let computed = outcomes.len() - hits;
        let mut results = Vec::with_capacity(points.len());
        for (x, (r, _)) in points.iter().zip(outcomes) {
            match r {
                Ok(o) => results.push(o),
                Err(e) => {
                    let at = if x == &nominal {
                        "at nominal".to_string()
                    } else {
                        let moved: Vec<String> = toleranced
                            .iter()
                            .filter(|&&i| x[i] != nominal[i])
                            .map(|&i| format!("{} {}", ports[i].name, if x[i] > nominal[i] { "high" } else { "low" }))
                            .collect();
                        format!("at the tolerance corner ({})", moved.join(", "))
                    };
                    let o = fail(&mut run, format!("{e} [{at}]"));
                    return (run, o, computed, hits);
                }
            }
        }
        let mut rigor = xs.iter().map(|v| v.rigor).min().unwrap_or(Rigor::Rigorous);
        rigor = rigor.min(results.iter().map(|o| o.rigor).min().unwrap());
        let loose: Vec<&str> = toleranced.iter().filter(|&&i| !ports[i].monotone).map(|&i| ports[i].name.as_str()).collect();
        if !exhaustive {
            rigor = rigor.min(Rigor::Asymptotic);
            run.rigor_note = Some(format!("{} toleranced inputs: tolerance ends sampled one at a time, not a bound", toleranced.len()));
        } else if !loose.is_empty() {
            rigor = rigor.min(Rigor::Asymptotic);
            run.rigor_note = Some(format!("not proven monotone in {}: corners are a sample, not a bound", loose.join(", ")));
        }
        run.notes = results[0].notes.clone();
        let outs = outputs
            .iter()
            .enumerate()
            .map(|(j, o)| {
                let lo = results.iter().map(|r| r.values[j] - r.errors[j]).fold(f64::INFINITY, f64::min);
                let hi = results.iter().map(|r| r.values[j] + r.errors[j]).fold(f64::NEG_INFINITY, f64::max);
                (o.name.clone(), Val::new(results[0].values[j], lo, hi, o.dim, rigor))
            })
            .collect();
        (run, outs, computed, hits)
    }

    fn check(&self, id: &str, values: &BTreeMap<String, Result<Val, String>>, pins: &BTreeSet<String>) -> Check {
        let r = &self.doc.requirements[id];
        let mut c = Check {
            id: id.into(),
            text: r.text.clone(),
            owner: r.owner.clone(),
            of: r.of.clone(),
            status: Status::Unavailable,
            value: None,
            max: None,
            min: None,
            margin: None,
            unit: None,
            reason: None,
        };
        c.unit = [&r.max, &r.min]
            .into_iter()
            .flatten()
            .find_map(|t| split_quantity(t).ok().map(|(_, u)| u.to_string()).filter(|u| !u.is_empty()));
        let scope = EvalScope { values, lib: &self.lib, pins };
        let v = match Expr::parse(&r.of).and_then(|e| e.eval(&scope)) {
            Ok(v) => v,
            Err(e) => {
                c.reason = Some(e);
                return c;
            }
        };
        let limit = |text: &Option<String>| -> Result<Option<f64>, String> {
            let Some(t) = text else { return Ok(None) };
            let (q, _) = parse_quantity(t)?;
            if q.dim != v.dim {
                return Err(format!("limit `{t}` is {}, the value is {}", q.dim.symbol(), v.dim.symbol()));
            }
            Ok(Some(q.si))
        };
        match (limit(&r.max), limit(&r.min)) {
            (Ok(max), Ok(min)) => {
                c.max = max;
                c.min = min;
            }
            (Err(e), _) | (_, Err(e)) => {
                c.reason = Some(e);
                return c;
            }
        }
        let mut margin = f64::INFINITY;
        let mut failed = false;
        if let Some(m) = c.max {
            margin = margin.min(m - v.hi);
            failed |= v.lo > m;
        }
        if let Some(m) = c.min {
            margin = margin.min(v.lo - m);
            failed |= v.hi < m;
        }
        c.status = if failed {
            Status::Failed
        } else if margin >= 0.0 {
            Status::Met
        } else {
            Status::Undecided
        };
        c.margin = margin.is_finite().then_some(margin);
        c.value = Some(v);
        c
    }

    /// Leaves (toleranced params and material properties) upstream of a
    /// requirement.
    pub fn leaves(&self, requirement: &str, eval: &Evaluation) -> Vec<String> {
        let mut out = BTreeSet::new();
        let mut stack: Vec<String> = Expr::parse(&self.doc.requirements[requirement].of).map(|e| e.deps()).unwrap_or_default();
        let mut seen = BTreeSet::new();
        while let Some(id) = stack.pop() {
            if !seen.insert(id.clone()) {
                continue;
            }
            if Library::split(&id).is_some() || self.doc.params.contains_key(&id) {
                if eval.value(&id).is_some_and(|v| v.width() > 0.0) {
                    out.insert(id);
                }
                continue;
            }
            if let Some(n) = eval.nodes.get(&id) {
                stack.extend(n.deps.iter().cloned());
            }
        }
        out.into_iter().collect()
    }

    /// How much of a requirement's uncertainty each upstream leaf causes:
    /// the fraction of the interval width that disappears when that leaf is
    /// pinned at its nominal value. Sorted, largest first.
    pub fn contributions(&self, requirement: &str, eval: &Evaluation) -> Vec<(String, String, f64)> {
        let Some(base) = eval.checks.get(requirement).and_then(|c| c.value) else { return vec![] };
        if base.width() == 0.0 {
            return vec![];
        }
        use rayon::prelude::*;
        let mut out: Vec<(String, String, f64)> = self
            .leaves(requirement, eval)
            .into_par_iter()
            .filter_map(|leaf| {
                let pinned = self.evaluate_pinned(&BTreeSet::from([leaf.clone()]));
                let w = pinned.checks[requirement].value?.width();
                let owner = self.doc.owner(&leaf).unwrap_or("?").to_string();
                Some((leaf, owner, 1.0 - w / base.width()))
            })
            .collect();
        out.sort_by(|a, b| b.2.total_cmp(&a.2));
        out
    }
}

fn build_part(dir: &Path, id: &str, p: &crate::doc::PartDef) -> Result<analysis::Part, String> {
    let path = dir.join(&p.model);
    let model = cad::model::Model::read(&path).map_err(|e| format!("part `{id}`: {e}"))?;
    for name in p.params.keys() {
        if !model.params.contains_key(name) {
            return Err(format!(
                "part `{id}`: model has no parameter `{name}` (has: {})",
                model.params.keys().cloned().collect::<Vec<_>>().join(", ")
            ));
        }
    }
    let cell = match &p.cell {
        Some(c) => {
            let (q, _) = parse_quantity(c)?;
            if q.dim != Dim::LENGTH {
                return Err(format!("part `{id}`: cell must be a length"));
            }
            q.si
        }
        None => model.cell.ok_or_else(|| format!("part `{id}`: give a mesh `cell`"))?,
    };
    let params: Vec<String> = p.params.keys().cloned().collect();
    let base = path.parent().unwrap_or(Path::new(".")).to_path_buf();
    Ok(analysis::Part { model, base, params, cell, fdm: p.fdm.is_some() })
}

fn node(
    id: &str,
    kind: Kind,
    owner: &str,
    label: &str,
    unit: Option<String>,
    deps: Vec<String>,
    values: &BTreeMap<String, Result<Val, String>>,
) -> Node {
    Node {
        id: id.into(),
        kind,
        owner: owner.into(),
        label: label.into(),
        unit,
        deps,
        value: values.get(id).cloned().unwrap_or_else(|| Err("not evaluated".into())),
    }
}

fn check_unit(v: Val, unit: Option<&str>) -> Result<Val, String> {
    if let Some(u) = unit {
        let du = parse_unit(u)?.dim;
        if du != v.dim {
            return Err(format!("declared unit {u} is {}, but the expression gives {}", du.symbol(), v.dim.symbol()));
        }
    }
    Ok(v)
}

/// A parameter's value and tolerance as an interval. A tolerance is a
/// stated bound, so the interval is rigorous by construction.
pub fn param_value(value: &str, tol: Option<&Tol>) -> Result<Val, String> {
    let (q, unit) = parse_quantity(value)?;
    let v = q.si;
    let (lo, hi) = match tol {
        None => (v, v),
        Some(Tol::Sym(t)) => {
            if let Some(p) = percent(t) {
                (v - p * v.abs(), v + p * v.abs())
            } else {
                let (dq, du) = parse_quantity(t)?;
                if du.is_affine() {
                    return Err(format!("tolerance `{t}`: write temperature differences in K"));
                }
                if dq.dim != q.dim {
                    return Err(format!("tolerance `{t}` is {}, the value is {}", dq.dim.symbol(), q.dim.symbol()));
                }
                (v - dq.si.abs(), v + dq.si.abs())
            }
        }
        Some(Tol::Range([a, b])) => {
            let (qa, _) = parse_quantity(a)?;
            let (qb, _) = parse_quantity(b)?;
            if qa.dim != q.dim || qb.dim != q.dim {
                return Err(format!("range [{a}, {b}] does not match the value's unit"));
            }
            (qa.si.min(qb.si), qa.si.max(qb.si))
        }
    };
    let _ = unit;
    Val::new(v, lo, hi, q.dim, Rigor::Rigorous).map_err(|_| format!("value {value} lies outside its tolerance range"))
}

fn pinned_prop(p: &Property, pinned: bool) -> Property {
    let mut p = p.clone();
    if pinned {
        p.rel = 0.0;
    }
    p
}

struct EvalScope<'a> {
    values: &'a BTreeMap<String, Result<Val, String>>,
    lib: &'a Library,
    pins: &'a BTreeSet<String>,
}

impl Scope for EvalScope<'_> {
    fn var(&self, id: &str) -> Result<Val, String> {
        if let Some((m, p)) = Library::split(id) {
            let prop = self.lib.get(m, p)?;
            return Ok(pinned_prop(prop, self.pins.contains(id)).reference());
        }
        match self.values.get(id) {
            Some(Ok(v)) => Ok(*v),
            // Keep the root cause first, then the path it travelled.
            Some(Err(e)) => Err(format!("{e} (via {id})")),
            None => Err(format!("unknown id `{id}`")),
        }
    }
    fn table(&self, id: &str, arg: &Val) -> Result<Val, String> {
        let (m, p) = Library::split(id).ok_or_else(|| format!("`{id}(…)`: only material properties can be called"))?;
        pinned_prop(self.lib.get(m, p)?, self.pins.contains(id)).at(arg)
    }
}

/// Strongly connected components of `edges` (node → nodes it reads), each
/// listed after every component it reads from: evaluation order.
fn components(edges: &BTreeMap<String, Result<Vec<String>, String>>) -> Vec<Vec<String>> {
    struct T<'a> {
        edges: &'a BTreeMap<String, Result<Vec<String>, String>>,
        index: BTreeMap<&'a str, usize>,
        low: BTreeMap<&'a str, usize>,
        stack: Vec<&'a str>,
        on: BTreeSet<&'a str>,
        next: usize,
        out: Vec<Vec<String>>,
    }
    fn visit<'a>(t: &mut T<'a>, v: &'a str) {
        t.index.insert(v, t.next);
        t.low.insert(v, t.next);
        t.next += 1;
        t.stack.push(v);
        t.on.insert(v);
        if let Ok(deps) = &t.edges[v] {
            for w in deps {
                let Some((w, _)) = t.edges.get_key_value(w) else { continue };
                let w = w.as_str();
                if !t.index.contains_key(w) {
                    visit(t, w);
                    let lw = t.low[w];
                    let lv = t.low.get_mut(v).unwrap();
                    *lv = (*lv).min(lw);
                } else if t.on.contains(w) {
                    let iw = t.index[w];
                    let lv = t.low.get_mut(v).unwrap();
                    *lv = (*lv).min(iw);
                }
            }
        }
        if t.low[v] == t.index[v] {
            let mut c = Vec::new();
            loop {
                let w = t.stack.pop().unwrap();
                t.on.remove(w);
                c.push(w.to_string());
                if w == v {
                    break;
                }
            }
            c.sort();
            t.out.push(c);
        }
    }
    let mut t = T { edges, index: BTreeMap::new(), low: BTreeMap::new(), stack: vec![], on: BTreeSet::new(), next: 0, out: vec![] };
    for k in edges.keys() {
        if !t.index.contains_key(k.as_str()) {
            visit(&mut t, k);
        }
    }
    t.out
}

/// Topological order of `edges` (node → nodes it reads). Nodes on a cycle
/// are returned separately with the cycle spelled out.
fn topo(edges: &BTreeMap<String, Result<Vec<String>, String>>) -> (Vec<String>, BTreeMap<String, String>) {
    #[derive(Clone, Copy, PartialEq)]
    enum Mark {
        New,
        Open,
        Done,
    }
    let mut mark: BTreeMap<&str, Mark> = edges.keys().map(|k| (k.as_str(), Mark::New)).collect();
    let mut order = Vec::new();
    let mut cyclic = BTreeMap::new();
    fn visit<'a>(
        n: &'a str,
        edges: &'a BTreeMap<String, Result<Vec<String>, String>>,
        mark: &mut BTreeMap<&'a str, Mark>,
        path: &mut Vec<&'a str>,
        order: &mut Vec<String>,
        cyclic: &mut BTreeMap<String, String>,
    ) {
        match mark.get(n) {
            Some(Mark::Done) | None => return,
            Some(Mark::Open) => {
                let start = path.iter().position(|p| *p == n).unwrap();
                let cycle: Vec<&str> = path[start..].iter().copied().chain([n]).collect();
                for c in &path[start..] {
                    cyclic.insert(c.to_string(), cycle.join(" → "));
                }
                return;
            }
            Some(Mark::New) => {}
        }
        mark.insert(n, Mark::Open);
        path.push(n);
        if let Ok(deps) = &edges[n] {
            for d in deps {
                visit(d, edges, mark, path, order, cyclic);
            }
        }
        path.pop();
        mark.insert(n, Mark::Done);
        order.push(n.to_string());
    }
    for k in edges.keys() {
        visit(k, edges, &mut mark, &mut Vec::new(), &mut order, &mut cyclic);
    }
    (order, cyclic)
}

impl Project {
    /// A copy with some parameter values replaced (`id`, `"value unit"`),
    /// sharing this project's cache: the what-if behind `--set`.
    pub fn with_values(&self, changes: &[(String, String)]) -> Result<Project, String> {
        let mut doc = self.doc.clone();
        for (id, v) in changes {
            let p = doc.params.get_mut(id).ok_or_else(|| format!("`{id}` is not a parameter"))?;
            p.value = v.clone();
        }
        Project::new(doc, &self.dir, self.cache.clone())
    }
}
