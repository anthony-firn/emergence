//! A self-contained HTML report: one page a cross-functional review can
//! open, share and attach. No scripts, no external assets.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write;

use crate::expr::Expr;
use crate::graph::{Evaluation, Kind, Project, Status};
use crate::report;
use crate::units::{parse_unit, sig, Dim};

fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
}

/// Team colours (hue, chroma chosen to stay distinct in light and dark).
const HUES: [&str; 8] = ["#2f6fdf", "#d9622b", "#2e9e6b", "#9a52c7", "#c43d6a", "#1f9bb3", "#a8871c", "#6b7280"];

fn shown(si: f64, dim: Dim, unit: Option<&str>) -> f64 {
    match unit.map(parse_unit) {
        Some(Ok(u)) if u.dim == dim => u.from_si(si),
        _ => si,
    }
}

const STYLE: &str = r#"
:root{--bg:#f7f7f5;--card:#fff;--ink:#1d1f23;--muted:#62666d;--line:#e3e3df;--met:#1f8a4c;--fail:#c2332b;--und:#b7791f;--na:#7a7f87;--band:#2f6fdf}
@media (prefers-color-scheme:dark){:root{--bg:#141518;--card:#1c1e22;--ink:#e9e9e6;--muted:#a0a4ab;--line:#2c2f35;--met:#3fb873;--fail:#ef5d55;--und:#e0a33a;--na:#8b9099;--band:#6c9cf0}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif}
main{max-width:1080px;margin:0 auto;padding:28px 16px 60px}
h1{font-size:26px;margin:0 0 4px}h2{font-size:18px;margin:36px 0 12px}p.lead{color:var(--muted);margin:0 0 16px;max-width:760px}
.chips{display:flex;gap:8px;flex-wrap:wrap;margin:12px 0}.chip{border:1px solid var(--line);border-radius:999px;padding:2px 10px;font-size:13px;background:var(--card)}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(320px,1fr));gap:14px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px}
.req .head{display:flex;justify-content:space-between;align-items:baseline;gap:8px}.req .id{font-weight:650}
.badge{font-size:12px;font-weight:650;letter-spacing:.02em;border-radius:6px;padding:1px 8px;color:#fff;white-space:nowrap}
.met{background:var(--met)}.failed{background:var(--fail)}.undecided{background:var(--und)}.unavailable{background:var(--na)}
.text{color:var(--muted);font-size:13.5px;margin:4px 0 8px}.value{font-variant-numeric:tabular-nums;font-size:15px}
.small{font-size:12.5px;color:var(--muted)}.owner{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:5px;vertical-align:baseline}
svg text{fill:var(--ink)}.contrib{margin-top:8px}.contrib div{display:flex;align-items:center;gap:6px;font-size:12.5px;margin:2px 0}
.contrib .bar{height:7px;border-radius:4px;flex:none}
table{width:100%;border-collapse:collapse;font-size:13.5px}td,th{padding:5px 8px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{color:var(--muted);font-weight:600}td.num{font-variant-numeric:tabular-nums;white-space:nowrap}
.scroll{overflow-x:auto}ul{margin:4px 0 0 18px;padding:0}li{margin:1px 0}.err{color:var(--fail)}
input,select,button{font:inherit;color:var(--ink);background:var(--bg);border:1px solid var(--line);border-radius:6px;padding:2px 6px}
input.edit{font-variant-numeric:tabular-nums;width:auto}input.edit:focus{outline:2px solid var(--band)}input.bad{border-color:var(--fail)}input.foreign{border-style:dashed}
button.why{margin-top:8px;font-size:12.5px;cursor:pointer}.live{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin-top:8px}
#toast{flex-basis:100%}#toast:empty{display:none}img.part{display:block;width:100%;height:auto;border-radius:8px;margin:8px 0;background:var(--bg)}#toast.err{color:var(--fail)}ul.feed{list-style:none;margin:0}ul.feed li{padding:6px 0;border-bottom:1px solid var(--line)}.impact{margin-left:12px}
"#;

const LIVE_SCRIPT: &str = r#"
const $ = (q, el) => (el || document).querySelector(q);
const $$ = (q) => [...document.querySelectorAll(q)];
const remember = (id) => { const el = $('#' + id); try { el.value = localStorage.getItem('fp-' + id) || ''; } catch (e) {} el.onchange = () => { try { localStorage.setItem('fp-' + id, el.value); } catch (e) {} mark(); }; };
function toast(text, bad) { const t = $('#toast'); t.textContent = text; t.className = 'small' + (bad ? ' err' : ''); }
function mark() { const team = $('#team').value; $$('input.edit').forEach(i => i.classList.toggle('foreign', !!team && i.dataset.owner !== team)); }
async function edit(el) {
  const team = $('#team').value;
  if (team && el.dataset.owner !== team && !confirm(el.dataset.id + ' is owned by ' + el.dataset.owner + '. Change it anyway?')) { el.value = el.defaultValue; return; }
  el.disabled = true; toast('evaluating ' + el.dataset.id + ' = ' + el.value + ' …');
  const who = ($('#who').value || 'someone') + (team ? ' (' + team + ')' : '');
  const r = await fetch('api/set', { method: 'POST', body: JSON.stringify({ id: el.dataset.id, value: el.value, who }) });
  const j = await r.json();
  el.disabled = false;
  if (!j.ok) { el.classList.add('bad'); toast(j.error, true); return; }
  version = j.version; await refresh(); toast(j.summary);
}
async function why(b) {
  b.disabled = true; b.textContent = 'working out the contributions …';
  const r = await fetch('api/why?req=' + encodeURIComponent(b.dataset.req));
  b.nextElementSibling.innerHTML = await r.text(); b.remove();
}
function bind() {
  remember('who'); remember('team'); mark();
  $$('input.edit').forEach(i => { i.onchange = () => edit(i); i.onkeydown = (e) => { if (e.key === 'Escape') { i.value = i.defaultValue; i.blur(); } }; });
  $$('button.why').forEach(b => b.onclick = () => why(b));
}
async function refresh() {
  const t = await (await fetch('./')).text();
  const doc = new DOMParser().parseFromString(t, 'text/html');
  const keep = $('#toast').textContent;
  $('main').replaceWith(doc.querySelector('main'));
  $('#toast').textContent = keep; bind();
}
async function poll() {
  try {
    const j = await (await fetch('api/state')).json();
    const busy = document.activeElement && document.activeElement.classList.contains('edit');
    if (j.version !== version && !busy) { version = j.version; await refresh(); }
  } catch (e) {}
  setTimeout(poll, 1500);
}
bind(); poll();
"#;

/// One entry in the live editor's activity feed.
#[derive(Clone, Debug)]
pub struct Activity {
    pub time: String,
    pub who: String,
    pub what: String,
    /// Requirements whose status the change moved, already formatted.
    pub impact: Vec<String>,
}

/// What a live page adds to the static report.
pub struct Live<'a> {
    pub version: u64,
    pub activity: &'a [Activity],
}

/// The static report.
pub fn page(p: &Project, e: &Evaluation) -> String {
    render(p, e, None)
}

/// Team colours, stable for a project.
pub fn colours(p: &Project, e: &Evaluation) -> BTreeMap<String, &'static str> {
    let mut teams: Vec<String> = e.nodes.values().filter(|n| n.kind != Kind::Material).map(|n| n.owner.clone()).collect();
    teams.extend(p.doc.teams.keys().cloned());
    teams.sort();
    teams.dedup();
    teams.iter().enumerate().map(|(i, t)| (t.clone(), HUES[i % HUES.len()])).collect()
}

/// Where a requirement's uncertainty comes from, as bars.
pub fn contributions_html(p: &Project, e: &Evaluation, req: &str) -> String {
    let colour = colours(p, e);
    let parts: Vec<_> = p.contributions(req, e).into_iter().filter(|x| x.2 > 0.005).take(5).collect();
    let mut s = String::new();
    if parts.is_empty() {
        return s;
    }
    let _ = write!(s, r#"<div class="contrib small">Uncertainty comes from"#);
    for (id, owner, f) in parts {
        let _ = write!(
            s,
            r#"<div><span class="bar" style="width:{:.0}px;background:{}"></span>{} <span>({}) {:.0}%</span></div>"#,
            (90.0 * f).max(3.0),
            colour.get(&owner).copied().unwrap_or("#888"),
            esc(&id),
            esc(&owner),
            100.0 * f
        );
    }
    let _ = write!(s, "</div>");
    s
}

pub fn render(p: &Project, e: &Evaluation, live: Option<&Live>) -> String {
    let mut teams: Vec<String> = e.nodes.values().filter(|n| n.kind != Kind::Material).map(|n| n.owner.clone()).collect();
    teams.extend(p.doc.teams.keys().cloned());
    teams.sort();
    teams.dedup();
    let colour: BTreeMap<String, &str> = teams.iter().enumerate().map(|(i, t)| (t.clone(), HUES[i % HUES.len()])).collect();
    let dot = |owner: &str| format!(r#"<span class="owner" style="background:{}"></span>"#, colour.get(owner).copied().unwrap_or("#888"));

    let mut s = String::new();
    let _ = write!(
        s,
        r#"<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{}</title><style>{STYLE}</style></head><body><main>"#,
        esc(&p.doc.name)
    );
    let _ = write!(s, "<h1>{}</h1>", esc(&p.doc.name));
    if !p.doc.description.is_empty() {
        let _ = write!(s, r#"<p class="lead">{}</p>"#, esc(&p.doc.description));
    }
    let count = |st: Status| e.checks.values().filter(|c| c.status == st).count();
    let _ = write!(s, r#"<div class="chips">"#);
    for st in [Status::Met, Status::Failed, Status::Undecided, Status::Unavailable] {
        let n = count(st);
        if n > 0 {
            let _ = write!(s, r#"<span class="chip"><span class="badge {}">{}</span> {n}</span>"#, st.name().to_lowercase(), st.name());
        }
    }
    for t in &teams {
        let _ = write!(s, r#"<span class="chip">{}{}</span>"#, dot(t), esc(t));
    }
    let _ = write!(s, "</div>");
    if let Some(l) = live {
        let _ = write!(
            s,
            r#"<div class="card live"><label>You are <input id="who" placeholder="name" size="12"></label> <label>on <select id="team"><option value="">any team</option>{}</select></label> <span class="small">Edit any parameter below; everyone viewing this page sees the result. Saved to the project file.</span><div id="toast" class="small"></div></div>"#,
            teams.iter().map(|t| format!(r#"<option>{}</option>"#, esc(t))).collect::<String>()
        );
        if !l.activity.is_empty() {
            let _ = write!(s, r#"<h2>Activity</h2><div class="card"><ul class="feed">"#);
            for a in l.activity.iter().rev().take(12) {
                let _ = write!(s, r#"<li><span class="small">{}</span> <b>{}</b> {}"#, esc(&a.time), esc(&a.who), esc(&a.what));
                for i in &a.impact {
                    let _ = write!(s, r#"<div class="small impact">{}</div>"#, esc(i));
                }
                let _ = write!(s, "</li>");
            }
            let _ = write!(s, "</ul></div>");
        }
    }

    // Requirements.
    let _ = write!(s, r#"<h2>Requirements</h2><div class="grid">"#);
    for c in e.checks.values() {
        let _ = write!(
            s,
            r#"<div class="card req"><div class="head"><span class="id">{}</span><span class="badge {}">{}</span></div><div class="text">{}{}</div>"#,
            esc(&c.id),
            c.status.name().to_lowercase(),
            c.status.name(),
            dot(&c.owner),
            esc(&c.text)
        );
        match c.value {
            Some(v) => {
                let u = c.unit.as_deref();
                let _ = write!(s, r#"<div class="value">{}</div>"#, esc(&report::val(&v, u)));
                let margin = c.margin.map(|m| format!(" · worst-case margin {}", report::delta(m, v.dim, u))).unwrap_or_default();
                let _ = write!(s, r#"<div class="small">{} · {}{}</div>"#, esc(&c.of), v.rigor.name(), esc(&margin));
                // Interval against the limit(s).
                let (lo, hi, nom) = (shown(v.lo, v.dim, u), shown(v.hi, v.dim, u), shown(v.nom, v.dim, u));
                let lims: Vec<f64> = [c.min, c.max].iter().flatten().map(|l| shown(*l, v.dim, u)).collect();
                let mut a = lims.iter().copied().fold(lo, f64::min);
                let mut b = lims.iter().copied().fold(hi, f64::max);
                let pad = 0.12 * (b - a).max(1e-12 * nom.abs().max(1.0));
                a -= pad;
                b += pad;
                let x = |t: f64| 8.0 + 284.0 * (t - a) / (b - a);
                let tint = match c.status {
                    Status::Met => "var(--met)",
                    Status::Failed => "var(--fail)",
                    Status::Undecided => "var(--und)",
                    Status::Unavailable => "var(--na)",
                };
                let _ = write!(
                    s,
                    r#"<svg viewBox="0 0 300 46" width="100%" role="img" aria-label="interval against limit"><line x1="8" y1="22" x2="292" y2="22" stroke="var(--line)" stroke-width="2"/><rect x="{:.1}" y="15" width="{:.1}" height="14" rx="4" fill="{tint}" opacity="0.35"/><line x1="{:.1}" y1="11" x2="{:.1}" y2="33" stroke="{tint}" stroke-width="2.5"/>"#,
                    x(lo),
                    (x(hi) - x(lo)).max(2.0),
                    x(nom),
                    x(nom)
                );
                for l in &lims {
                    let _ = write!(
                        s,
                        r#"<line x1="{0:.1}" y1="6" x2="{0:.1}" y2="38" stroke="var(--ink)" stroke-dasharray="3 2"/><text x="{0:.1}" y="46" font-size="10" text-anchor="middle">limit {1}</text>"#,
                        x(*l),
                        esc(&sig(*l, 4))
                    );
                }
                let _ = write!(s, "</svg>");
                if v.width() > 0.0 {
                    if live.is_some() {
                        let _ = write!(s, r#"<button class="why" data-req="{}">Where does the uncertainty come from?</button><div class="why-out"></div>"#, esc(&c.id));
                    } else {
                        s.push_str(&contributions_html(p, e, &c.id));
                    }
                }
            }
            None => {
                let _ = write!(s, r#"<div class="err small">{}</div>"#, esc(c.reason.as_deref().unwrap_or("not evaluated")));
            }
        }
        let _ = write!(s, "</div>");
    }
    let _ = write!(s, "</div>");

    // Dependency graph.
    let _ = write!(s, r#"<h2>How the teams' numbers connect</h2><div class="card scroll">{}</div>"#, graph_svg(p, e, &colour));

    // Team tables.
    let _ = write!(s, "<h2>Values by team</h2>");
    for t in &teams {
        let rows: Vec<_> = e.nodes.values().filter(|n| &n.owner == t && matches!(n.kind, Kind::Param | Kind::Derived | Kind::Output)).collect();
        if rows.is_empty() {
            continue;
        }
        let _ = write!(
            s,
            r#"<div class="card scroll" style="margin-bottom:12px"><div style="font-weight:650;margin-bottom:4px">{}{}</div><table><tr><th>id</th><th>what</th><th>value</th><th>rigor</th></tr>"#,
            dot(t),
            esc(t)
        );
        for n in rows {
            let (val, rig) = match &n.value {
                Ok(v) => (report::val(v, n.unit.as_deref()), if v.is_exact() { String::new() } else { v.rigor.name().to_string() }),
                Err(err) => (format!("error: {err}"), String::new()),
            };
            let kind = if n.kind == Kind::Param { "" } else { "= " };
            let cell = match (live, p.doc.params.get(&n.id)) {
                (Some(_), Some(def)) => {
                    let tol = match &def.tol {
                        Some(crate::doc::Tol::Sym(t)) => format!(" ± {t}"),
                        Some(crate::doc::Tol::Range([a, b])) => format!(" in [{a}, {b}]"),
                        None => String::new(),
                    };
                    format!(
                        r#"<input class="edit" data-id="{}" data-owner="{}" value="{}" size="{}"><span class="small">{}</span><div class="small">{}</div>"#,
                        esc(&n.id),
                        esc(&n.owner),
                        esc(&def.value),
                        def.value.chars().count().max(6) + 2,
                        esc(&tol),
                        esc(&val)
                    )
                }
                _ => esc(&val),
            };
            let _ = write!(
                s,
                r#"<tr><td>{kind}{}</td><td class="small">{}</td><td class="num">{cell}</td><td class="small">{}</td></tr>"#,
                esc(&n.id),
                esc(&n.label),
                rig
            );
        }
        let _ = write!(s, "</table></div>");
    }

    if !e.loops.is_empty() {
        let _ = write!(s, r#"<h2>Feedback loops</h2><div class="grid">"#);
        for l in &e.loops {
            let _ = write!(s, r#"<div class="card"><div style="font-weight:650">{}</div><div class="small">unknown: {}</div>"#, esc(&l.nodes.join(" → ")), esc(&l.unknowns.join(", ")));
            match &l.error {
                Some(err) => {
                    let _ = write!(s, r#"<div class="err small">{}</div>"#, esc(err));
                }
                None => {
                    let _ = write!(
                        s,
                        r#"<div class="small">settled in {} iterations · {}</div>"#,
                        l.iterations,
                        if l.existence { "the solution interval maps into itself: an operating point provably exists in it" } else { "existence not proven" }
                    );
                }
            }
            let _ = write!(s, "</div>");
        }
        let _ = write!(s, "</div>");
    }
    // Models.
    let _ = write!(s, r#"<h2>Models, regimes and assumptions</h2><div class="grid">"#);
    for (id, r) in &e.runs {
        let _ = write!(
            s,
            r#"<div class="card"><div style="font-weight:650">{}{} <span class="small">{}</span></div><div class="small">{}</div>"#,
            dot(&r.owner),
            esc(id),
            esc(&r.kind),
            esc(&r.summary)
        );
        if let Some(err) = &r.error {
            let _ = write!(s, r#"<div class="err small">{}</div>"#, esc(err));
        } else if r.kind == "part" {
            let src = match live {
                Some(l) => Some(format!("api/render?part={}&v={}", esc(id), l.version)),
                None => crate::render::part_png(p, e, id, 480, 320).ok().map(|png| format!("data:image/png;base64,{}", crate::render::base64(&png))),
            };
            if let Some(src) = src {
                let _ = write!(s, r#"<img class="part" src="{src}" alt="{} at the project's values" width="480" height="320">"#, esc(id));
            }
        }
        if !r.notes.is_empty() || r.rigor_note.is_some() {
            let _ = write!(s, r#"<ul class="small">"#);
            for n in &r.notes {
                let _ = write!(s, "<li>{}</li>", esc(n));
            }
            if let Some(n) = &r.rigor_note {
                let _ = write!(s, "<li>{} tolerance corners; {}</li>", r.corners, esc(n));
            }
            let _ = write!(s, "</ul>");
        }
        let _ = write!(s, r#"<div class="small" style="margin-top:6px">Assumes</div><ul class="small">"#);
        for a in &r.assumptions {
            let _ = write!(s, "<li>{}</li>", esc(a));
        }
        let _ = write!(s, "</ul></div>");
    }
    let _ = write!(s, "</div>");
    let _ = write!(
        s,
        r#"<p class="small" style="margin-top:28px">Intervals contain the true value given the stated tolerances and model assumptions. <b>rigorous</b>: a proven bound; <b>validated</b>: handbook or checked data; <b>asymptotic</b>: an estimate that can be exceeded outside its regime. Generated by <code>fluid project</code>.</p></main>"#
    );
    if let Some(l) = live {
        let _ = write!(s, "<script>let version = {};{LIVE_SCRIPT}</script>", l.version);
    }
    s.push_str("</body></html>");
    s
}

/// Layered drawing of the graph: parameters on the left, requirements on
/// the right, each analysis or part as one box.
fn graph_svg(p: &Project, e: &Evaluation, colour: &BTreeMap<String, &str>) -> String {
    // Graph units: params, derived, models, requirements.
    let unit_of = |id: &str| -> Option<String> {
        if p.doc.params.contains_key(id) || p.doc.derived.contains_key(id) {
            return Some(id.to_string());
        }
        let (head, _) = id.rsplit_once('.')?;
        p.model(head).map(|_| head.to_string())
    };
    let mut deps: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    let mut owner: BTreeMap<String, String> = BTreeMap::new();
    for (id, d) in &p.doc.params {
        deps.insert(id.clone(), BTreeSet::new());
        owner.insert(id.clone(), d.owner.clone());
    }
    for (id, d) in &p.doc.derived {
        let ds = Expr::parse(&d.expr).map(|x| x.deps()).unwrap_or_default();
        deps.insert(id.clone(), ds.iter().filter_map(|d| unit_of(d)).collect());
        owner.insert(id.clone(), d.owner.clone());
    }
    for (id, r) in &e.runs {
        let ds: BTreeSet<String> = r
            .inputs
            .iter()
            .flat_map(|(_, text, _)| Expr::parse(text).map(|x| x.deps()).unwrap_or_default())
            .filter_map(|d| unit_of(&d))
            .collect();
        deps.insert(id.clone(), ds);
        owner.insert(id.clone(), r.owner.clone());
    }
    for (id, r) in &p.doc.requirements {
        let ds = Expr::parse(&r.of).map(|x| x.deps()).unwrap_or_default();
        deps.insert(id.clone(), ds.iter().filter_map(|d| unit_of(d)).collect());
        owner.insert(id.clone(), r.owner.clone());
    }
    // Column = longest path from a leaf; requirements in the last column.
    let mut depth: BTreeMap<String, usize> = BTreeMap::new();
    fn d_of(id: &str, deps: &BTreeMap<String, BTreeSet<String>>, memo: &mut BTreeMap<String, usize>, guard: usize) -> usize {
        if let Some(d) = memo.get(id) {
            return *d;
        }
        if guard > 64 {
            return 0;
        }
        // Provisional depth: an edge back into a feedback loop reads 0, so a
        // loop is drawn once instead of being unrolled column after column.
        memo.insert(id.to_string(), 0);
        let d = deps.get(id).map(|s| s.iter().map(|x| 1 + d_of(x, deps, memo, guard + 1)).max().unwrap_or(0)).unwrap_or(0);
        memo.insert(id.to_string(), d);
        d
    }
    for id in deps.keys() {
        d_of(id, &deps, &mut depth, 0);
    }
    let last = depth.values().copied().max().unwrap_or(0) + 1;
    for id in p.doc.requirements.keys() {
        depth.insert(id.clone(), last);
    }
    let mut cols: BTreeMap<usize, Vec<String>> = BTreeMap::new();
    for (id, d) in &depth {
        cols.entry(*d).or_default().push(id.clone());
    }
    for v in cols.values_mut() {
        v.sort_by(|a, b| (&owner[a], a).cmp(&(&owner[b], b)));
    }
    let (cw, rh, bw, bh) = (190.0, 30.0, 160.0, 22.0);
    let rows = cols.values().map(|v| v.len()).max().unwrap_or(1);
    let (w, h) = (cols.len() as f64 * cw + 20.0, rows as f64 * rh + 20.0);
    let mut pos: BTreeMap<String, (f64, f64)> = BTreeMap::new();
    for (ci, (_, ids)) in cols.iter().enumerate() {
        let off = (rows - ids.len()) as f64 * rh / 2.0;
        for (ri, id) in ids.iter().enumerate() {
            pos.insert(id.clone(), (10.0 + ci as f64 * cw, 10.0 + off + ri as f64 * rh));
        }
    }
    let mut s = format!(r#"<svg viewBox="0 0 {w:.0} {h:.0}" width="{w:.0}" height="{h:.0}" role="img" aria-label="dependency graph">"#);
    for (id, ds) in &deps {
        let (x1, y1) = pos[id];
        for d in ds {
            if let Some(&(x0, y0)) = pos.get(d) {
                let (sx, sy, ex, ey) = (x0 + bw, y0 + bh / 2.0, x1, y1 + bh / 2.0);
                let c = colour.get(&owner[d]).copied().unwrap_or("#888");
                let _ = write!(
                    s,
                    r#"<path d="M{sx:.1},{sy:.1} C{:.1},{sy:.1} {:.1},{ey:.1} {ex:.1},{ey:.1}" fill="none" stroke="{c}" stroke-opacity="0.45" stroke-width="1.4"/>"#,
                    sx + 30.0,
                    ex - 30.0
                );
            }
        }
    }
    for (id, &(x, y)) in &pos {
        let c = colour.get(&owner[id]).copied().unwrap_or("#888");
        let is_model = e.runs.contains_key(id);
        let is_req = p.doc.requirements.contains_key(id);
        let fill = if is_req {
            match e.checks[id].status {
                Status::Met => "var(--met)",
                Status::Failed => "var(--fail)",
                Status::Undecided => "var(--und)",
                Status::Unavailable => "var(--na)",
            }
        } else {
            "var(--card)"
        };
        let label: String = if id.chars().count() > 24 { id.chars().take(23).collect::<String>() + "…" } else { id.clone() };
        let _ = write!(
            s,
            r#"<g><title>{}</title><rect x="{x:.1}" y="{y:.1}" width="{bw}" height="{bh}" rx="{}" fill="{fill}" stroke="{c}" stroke-width="{}"/><text x="{:.1}" y="{:.1}" font-size="11.5" {}>{}</text></g>"#,
            esc(&format!("{id} [{}]", owner[id])),
            if is_model { 3 } else { 11 },
            if is_model { 2.4 } else { 1.4 },
            x + 9.0,
            y + 15.0,
            if is_req { r#"style="fill:#fff;font-weight:600""# } else { "" },
            esc(&label)
        );
    }
    s.push_str("</svg>");
    s
}
