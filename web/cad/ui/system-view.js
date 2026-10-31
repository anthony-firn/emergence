// The System studio: the open project's linked solve, its headline outputs,
// every link and check, the stated assumptions and live certificates,
// coupled sensitivities, the project's optimisation, its reference designs and
// the command log any surface can replay. Everything specific to a project
// comes from its manifest.
import { project, onSolved, mechanical, scheduleSolve } from "../app.js";
import { DOMAIN_INFO } from "../domains/index.js";
import { compatible, label as unitLabel } from "../core/units.js";
import { getPath, setPath, numericPaths, slug, STUDIOS } from "../core/project.js";
import { buildScript } from "../domains/rebuild.js";
import { validParams } from "../repository.js";
import { node, svg, button, format, checkCards, certificates, assumptionList, section, run, download, toast, live, lineChart, miniForm } from "./kit.js";
const root = document.getElementById("system-studio");
let solution = null, study = null, optimum = null;
const manifest = () => project.store.docs.project;
const docsNow = () => ({ ...structuredClone(project.store.docs), mech: structuredClone(mechanical.state()) });
const portLabel = (ports, id) => { const p = ports?.outputs.get(id) || ports?.inputs.get(id); return p ? `${DOMAIN_INFO[p.domain]?.title || p.domain} · ${p.label}` : id; };
const outputLabel = (s, port) => manifest().headline.find(h => h.port === port)?.label || s?.ports.outputs.get(port)?.label || port;
// A what-if copy must still be a valid design: each domain's own validator decides.
const validate = (docs, domain) => { if (domain === "mech") { if (!validParams(docs.mech.params)) throw new Error("outside the mechanical model's ranges"); } else project.store.domains[domain].validate(docs[domain]); };
const quiet = p => p.catch(() => {});

function render() {
  const s = solution, m = manifest(); root.replaceChildren();
  const head = node("div", undefined, "studio-heading"), title = node("div");
  title.append(node("p", `SYSTEM / ${m.name.toUpperCase()}`, "panel-label"), node("h2", "One model across every domain."),
    node("p", "Each domain publishes typed outputs and accepts typed inputs. Links convert units and refuse mismatched dimensions. The solver runs the project's domains in dependency order, iterates the coupled loops, and reports whether they settle."));
  const tools = node("div", undefined, "studio-file-tools");
  tools.append(button("＋ Domain", () => document.getElementById("studio-add").click(), "", "system-add-domain"), button("↻ Solve now", () => scheduleSolve([]), "", "system-solve"));
  head.append(title, tools); root.append(head);
  if (!m.domains.length) { root.append(emptyProject()); log(); return; }
  if (!s) { root.append(node("p", "Solving the linked system on the Rust kernels…", "small muted")); return; }
  // Headline outputs and convergence.
  const hero = node("div", undefined, "system-hero");
  for (const h of m.headline) {
    const o = s.outputs[h.port], card = node("div", undefined, "hero-stat");
    card.append(node("span", outputLabel(s, h.port)), node("b", format(o?.value, o?.unit, 4)), node("small", s.ports.outputs.has(h.port) ? h.port : `${h.port} · not in this project`));
    const x = button("×", () => quiet(run("project.unpin", { port: h.port })), "unpin"); x.setAttribute("aria-label", `Remove ${outputLabel(s, h.port)} from the headline`); card.append(x);
    hero.append(card);
  }
  const conv = node("div", undefined, "hero-stat convergence");
  conv.append(node("span", "Coupled solve"), node("b", `${s.iterations} pass${s.iterations > 1 ? "es" : ""} · ${s.converged ? "settled" : "not settled"}`), node("small", `residual ${Number.isFinite(s.residual) ? s.residual.toExponential(1) : "—"} · ${s.ms.toFixed(0)} ms · order ${s.order.map(d => DOMAIN_INFO[d]?.title || d).join(" → ")}`));
  const hist = s.history.filter(Number.isFinite);
  if (hist.length > 1) conv.append(lineChart({ series: [{ points: hist.map((r, i) => [i + 2, Math.log10(Math.max(r, 1e-16))]) }], x: { label: "pass", min: 2 }, y: { label: "log₁₀ residual" }, width: 260, height: 110 }));
  hero.append(conv); root.append(hero);
  const pin = node("div", undefined, "rtl-actions pin-row"), pick = node("select"); pick.id = "system-pin"; pick.setAttribute("aria-label", "Output to headline");
  for (const [id, p] of s.ports.outputs) if (!m.headline.some(h => h.port === id)) { const o = node("option", `${DOMAIN_INFO[p.domain]?.title} · ${p.label}`); o.value = id; pick.append(o); }
  const named = node("input"); named.id = "system-pin-label"; named.placeholder = "Label (optional)"; named.maxLength = 60; named.setAttribute("aria-label", "Headline label");
  pin.append(pick, named, button("＋ Headline output", () => pick.value && quiet(run("project.pin", { port: pick.value, ...(named.value && { label: named.value }) })), "", "system-pin-add")); root.append(pin);
  root.append(requirementsPanel(s));
  // The link graph.
  const g = section(root, "Links", `${project.store.docs.system.links.filter(l => l.enabled && !s.dormant.includes(l.id)).length} enabled · ${Object.keys(s.applied).length} inputs driven${s.dormant.length ? ` · ${s.dormant.length} dormant` : ""}`);
  g.append(graph(s), linkTable(s));
  // Checks by domain.
  const c = section(root, "Every check", `${s.checks.filter(x => x.status === "met").length} met · ${s.checks.filter(x => x.status === "undecided").length} undecided · ${s.checks.filter(x => x.status === "unmet" || x.status === "unavailable").length} failing or unavailable`);
  const byDomain = {}; for (const x of s.checks) (byDomain[x.domain] ||= []).push(x);
  for (const [d, list] of Object.entries(byDomain)) { c.append(node("h3", `${DOMAIN_INFO[d]?.icon || ""} ${DOMAIN_INFO[d]?.title || d}`, "domain-heading")); checkCards(c, list); certificates(c, s.results[d]?.certificates); }
  if (s.issues.length) { const ul = node("ul", undefined, "link-issues"); for (const i of s.issues) ul.append(node("li", i.message)); c.append(node("h3", "Link and convergence notes", "domain-heading"), ul); }
  root.append(studyPanel(s), optimizer(m.optimize, s), tourPanel(s));
  const r = section(root, "Public reference designs", "What this project's documents are built from"), ul = node("ul", undefined, "credits");
  m.credits.forEach((x, i) => { const li = node("li"), a = node("a", x.title); a.href = x.url; a.rel = "noopener"; li.append(node("b", DOMAIN_INFO[x.domain]?.title || "Project"), a); if (x.licence) li.append(node("span", x.licence, "pill")); li.append(button("×", () => quiet(run("project.uncredit", { index: i + 1 })), "unpin")); if (x.note) li.append(node("small", x.note)); ul.append(li); });
  r.append(ul);
  authoring(r, "Credit a design", d => miniForm(d, { id: "credit", action: "＋ Credit", fields: [{ name: "title", label: "Design", wide: true }, { name: "url", label: "Link (https)", wide: true }, { name: "licence", label: "Licence" }, { name: "domain", label: "Domain", options: [["", "Whole project"], ...m.domains.map(x => [x, DOMAIN_INFO[x].title])] }, { name: "note", label: "Note", wide: true }], submit: a => run("project.credit", a) }));
  // Assumptions.
  const all = project.assumptions(), confirmed = all.filter(a => a.confirmed).length;
  const a = section(root, "Fully defined design", `${confirmed} of ${all.length} assumptions confirmed by you`);
  const meter = node("div", undefined, "meter"); const bar = node("i"); bar.style.width = `${100 * confirmed / Math.max(1, all.length)}%`; meter.append(bar); a.append(meter);
  a.append(node("p", "Nominal, illustrative and prescribed inputs are listed per domain. Checking one records that you have verified it against your own data; it does not change a number. Unchecked items mean the answers above inherit that assumption.", "small muted"));
  assumptionList(a);
  log();
}
function emptyProject() {
  const wrap = node("div", undefined, "empty-project");
  wrap.append(node("h3", "This project has no domains yet."), node("p", "Add the domains your design spans, empty to start from nothing or with a starter design to edit. Then link their outputs to each other's inputs here, state your requirements, and study what drives them.", "small muted"));
  const cards = node("div", undefined, "project-cards");
  for (const [d, info] of Object.entries(DOMAIN_INFO)) if (info.blurb) {
    const c = node("article", undefined, "project-card"); c.append(node("h4", `${info.icon} ${info.title}`), node("p", info.blurb));
    const row = node("div", undefined, "project-card-actions"); row.append(button("Add empty", () => quiet(run("project.add_domain", { domain: d, start: "empty" })), "accent", `system-add-${d}`), button("Add starter", () => quiet(run("project.add_domain", { domain: d, start: "starter" })), "", `system-starter-${d}`)); c.append(row); cards.append(c);
  }
  wrap.append(cards); return wrap;
}
function log() {
  const log = section(root, "Everything is a command", `${project.bus.log.length} commands this session`);
  const actions = node("div", undefined, "rtl-actions");
  actions.append(button("Copy edits as script", () => { navigator.clipboard?.writeText(JSON.stringify(project.bus.script(), null, 2)).then(() => toast("Script copied."), () => toast("Clipboard unavailable; use Download.")); }),
    button("↓ Download script", () => download(JSON.stringify(project.bus.script(), null, 2) + "\n", `${slug(manifest().name)}-script.json`, "application/json")),
    button("↓ Build script", () => { const script = buildScript(project.file().documents); download(JSON.stringify(script, null, 1) + "\n", `${slug(manifest().name)}-build.json`, "application/json"); toast(`${script.length} commands that build this project from a blank one.`); }, "", "system-build-script"));
  const replay = node("label", "↑ Replay script", "file-button"), file = node("input"); file.type = "file"; file.accept = ".json,application/json"; file.id = "system-replay"; replay.htmlFor = file.id; replay.append(file);
  file.addEventListener("change", async () => { try { const list = JSON.parse(await file.files[0].text()); const r = await project.bus.batch(list, { source: "script", label: "Replay script" }); toast(`Replayed ${r.length} commands as one undoable step.`); } catch (e) { toast(`Script rejected, nothing changed: ${e.message}`); } finally { file.value = ""; } });
  actions.append(replay); log.append(actions);
  const table = node("table", undefined, "data-table log-table"), thead = node("thead"), tr = node("tr"); for (const h of ["Time", "Source", "Command", "Result"]) tr.append(node("th", h)); thead.append(tr); table.append(thead);
  const body = node("tbody"); for (const e of project.bus.log.slice(-40).reverse()) { const row = node("tr"); row.dataset.ok = e.ok; for (const v of [new Date(e.time).toLocaleTimeString(), e.source, `${e.op} ${Object.entries(e.args || {}).map(([k, v]) => `${k}=${JSON.stringify(v)}`).join(" ")}`.slice(0, 160), e.ok ? e.summary || "✓" : `✗ ${e.error}`]) row.append(node("td", v)); body.append(row); }
  if (!project.bus.log.length) { const row = node("tr"); const td = node("td", "No commands yet. Edit any studio, or type one below."); td.colSpan = 4; row.append(td); body.append(row); }
  table.append(body); log.append(table);
  const agent = node("details", undefined, "agent-details"); agent.append(node("summary", "Drive it from code or an AI agent"));
  agent.append(node("p", "Browser console: window.emergence.run(\"project.add_domain\", { domain: \"optics\" }), .batch([...]), .commands(), .solution, .documents(), .file()."),
    node("p", `Headless, with the same commands and the same WebAssembly physics, on a downloaded project file: node web/cad/agent/cli.mjs solve --project ${slug(manifest().name)}.emergence.json · cli.mjs new my.emergence.json --example blank · node web/cad/agent/mcp.mjs --project my.emergence.json (an MCP server: claude mcp add emergence -- node web/cad/agent/mcp.mjs --project my.emergence.json).`));
  log.append(agent);
}
function graph(s) {
  const ids = project.active(), n = ids.length, W = 760, H = 440, cx = W / 2, cy = H / 2 + 6, R = 165;
  const pos = Object.fromEntries(ids.map((d, i) => { const a = -Math.PI / 2 + 2 * Math.PI * i / n; return [d, n === 1 ? [cx, cy] : [cx + 1.45 * R * Math.cos(a), cy + R * Math.sin(a)]]; }));
  const root = svg("svg", { viewBox: `0 0 ${W} ${H}`, class: "link-graph", role: "img", "aria-label": "Cross-domain link graph" });
  const issues = new Map(s.issues.filter(i => i.link).map(i => [i.link, i]));
  const pairs = {};
  for (const l of project.store.docs.system.links) { const a = l.from.split(".")[0], b = l.to.split(".")[0]; if (!pos[a] || !pos[b]) continue; const k = a + ">" + b; (pairs[k] ||= []).push(l); }
  for (const [k, links] of Object.entries(pairs)) {
    const [a, b] = k.split(">"), [x1, y1] = pos[a], [x2, y2] = pos[b], mx = (x1 + x2) / 2, my = (y1 + y2) / 2, dx = x2 - x1, dy = y2 - y1, bend = .18;
    const qx = mx - dy * bend, qy = my + dx * bend, state = links.some(l => issues.has(l.id)) ? "issue" : links.every(l => !l.enabled) ? "off" : "on";
    const path = svg("path", { d: `M${x1},${y1} Q${qx},${qy} ${x2},${y2}`, class: `edge edge-${state}`, "marker-end": "url(#arrow)" });
    path.append(svg("title", {}, links.map(l => `${l.from} → ${l.to}${s.applied[l.to] ? ` = ${format(s.applied[l.to].value, s.applied[l.to].unit, 4)}` : ""}${issues.get(l.id) ? ` · ${issues.get(l.id).message}` : ""}${l.enabled ? "" : " · disabled"}`).join("\n")));
    root.append(path, svg("text", { x: qx, y: qy, class: "edge-label", "text-anchor": "middle" }, links.length > 1 ? `${links.length} links` : links[0].to.split(".")[1]));
  }
  const defs = svg("defs"), marker = svg("marker", { id: "arrow", viewBox: "0 0 10 10", refX: 34, refY: 5, markerWidth: 6, markerHeight: 6, orient: "auto-start-reverse" });
  marker.append(svg("path", { d: "M0,0 L10,5 L0,10 z", class: "arrow" })); defs.append(marker); root.prepend(defs);
  for (const d of ids) {
    const [x, y] = pos[d], checks = s.checks.filter(c => c.domain === d), bad = checks.some(c => c.status === "unmet" || c.status === "unavailable"), open = checks.some(c => c.status === "undecided");
    const g = svg("g", { class: `domain-node ${bad ? "bad" : open ? "open" : "ok"}`, role: "button", tabindex: 0, "aria-label": `${DOMAIN_INFO[d].title}: ${checks.length} checks` });
    g.append(svg("circle", { cx: x, cy: y, r: 30 }), svg("text", { x, y: y - 3, "text-anchor": "middle", class: "node-icon" }, DOMAIN_INFO[d].icon), svg("text", { x, y: y + 13, "text-anchor": "middle", class: "node-label" }, DOMAIN_INFO[d].title));
    g.addEventListener("click", () => document.getElementById(`studio-tab-${DOMAIN_INFO[d].studio}`)?.click());
    g.addEventListener("keydown", e => { if (e.key === "Enter") g.dispatchEvent(new Event("click")); });
    root.append(g);
  }
  return root;
}
function linkTable(s) {
  const wrap = node("div", undefined, "link-table-wrap"), table = node("table", undefined, "data-table link-table"), thead = node("thead"), tr = node("tr");
  for (const h of ["On", "From", "To", "Value", "Note", ""]) tr.append(node("th", h)); thead.append(tr); table.append(thead);
  const body = node("tbody"), issues = new Map(s.issues.filter(i => i.link).map(i => [i.link, i]));
  for (const l of project.store.docs.system.links) {
    const dormant = s.dormant.includes(l.id), row = node("tr"); row.dataset.state = dormant ? "off" : issues.has(l.id) ? "issue" : l.enabled ? "on" : "off";
    const on = node("input"); on.type = "checkbox"; on.checked = l.enabled; on.disabled = dormant; on.setAttribute("aria-label", `Enable ${l.from} to ${l.to}`); on.addEventListener("change", () => quiet(run("system.toggle", { id: l.id, enabled: on.checked })));
    const td = node("td"); td.append(on); row.append(td);
    const v = s.applied[l.to], o = s.outputs[l.from], missing = [l.from, l.to].map(p => p.split(".")[0]).find(d => !project.active().includes(d));
    for (const text of [portLabel(s.ports, l.from), portLabel(s.ports, l.to), dormant ? `Dormant: ${DOMAIN_INFO[missing]?.title || missing} is not in this project` : issues.get(l.id)?.message || (v ? format(v.value, v.unit, 4) : o ? `${format(o.value, o.unit, 4)} (not applied)` : "—"), l.note || ""]) row.append(node("td", text));
    const del = node("td"); del.append(button("Remove", () => quiet(run("system.unlink", { id: l.id })))); row.append(del); body.append(row);
  }
  if (!project.store.docs.system.links.length) { const row = node("tr"), td = node("td", "No links yet. Choose an output and an input with matching units below."); td.colSpan = 6; row.append(td); body.append(row); }
  table.append(body); wrap.append(table);
  // New link: outputs, then the inputs whose dimensions agree.
  const form = node("div", undefined, "link-form"), from = node("select"), to = node("select"), note = node("input"); note.placeholder = "Why these are the same quantity"; note.maxLength = 160;
  from.setAttribute("aria-label", "Output port"); to.setAttribute("aria-label", "Input port"); note.setAttribute("aria-label", "Link note");
  for (const [id, p] of s.ports.outputs) { const o = node("option", `${DOMAIN_INFO[p.domain]?.title} · ${p.label} (${unitLabel(p.unit) || "—"})`); o.value = id; from.append(o); }
  const fill = () => { to.replaceChildren(); const out = s.ports.outputs.get(from.value); for (const [id, p] of s.ports.inputs) if (out && compatible(out.unit, p.unit) && p.domain !== out.domain) { const o = node("option", `${DOMAIN_INFO[p.domain]?.title} · ${p.label} (${unitLabel(p.unit) || "—"})`); o.value = id; to.append(o); } if (!to.options.length) to.append(Object.assign(node("option", "No input with matching units"), { value: "" })); };
  from.addEventListener("change", fill); fill();
  form.append(node("span", "New link", "panel-label"), from, node("span", "→"), to, note, button("＋ Link", () => to.value && quiet(run("system.link", { from: from.value, to: to.value, ...(note.value && { note: note.value }) })), "accent", "system-add-link"));
  wrap.append(form);
  return wrap;
}
function studyPanel(s) {
  const m = manifest(), sec = node("section", undefined, "domain-section study"), top = node("div", undefined, "section-top"); top.append(node("h2", "What drives a requirement?"), node("span", "Finite differences through the full coupled solve")); sec.append(top);
  const pick = node("select"); pick.id = "study-output"; pick.setAttribute("aria-label", "Output to study");
  const ports = [...m.headline.map(h => h.port).filter(p => s.ports.outputs.has(p)), ...[...s.ports.outputs.keys()].filter(p => !m.headline.some(h => h.port === p))];
  for (const port of ports) { const o = node("option", `${outputLabel(s, port)} (${port})`); o.value = port; pick.append(o); }
  pick.value = study?.port && ports.includes(study.port) ? study.port : ports[0] ?? "";
  const go = button("Run sensitivity study", async () => {
    if (!m.studies.length || !pick.value) return;
    go.disabled = true; const port = pick.value, base = docsNow(), b = await project.system.fork(base).solve(base.system), y0 = b.outputs[port]?.value, rows = [];
    for (const [i, x] of m.studies.entries()) {
      go.textContent = `Solving ${i + 1} of ${m.studies.length}…`;
      const docs = structuredClone(base), x0 = getPath(docs[x.domain], x.path);
      try {
        if (!Number.isFinite(x0)) throw new Error(`${x.path} is not a number`);
        setPath(docs[x.domain], x.path, Math.min(x.max ?? Infinity, Math.max(x.min ?? -Infinity, x0 * 1.05))); validate(docs, x.domain);
        const r = await project.system.fork(docs).solve(docs.system), y = r.outputs[port]?.value, dx = getPath(docs[x.domain], x.path) / x0 - 1;
        rows.push({ label: x.label, x0, s: dx ? ((y - y0) / y0) / dx : NaN, y });
      } catch (e) { rows.push({ label: x.label, x0, s: NaN, error: e.message }); }
    }
    study = { port, y0, rows: rows.sort((a, b) => (Math.abs(b.s) || 0) - (Math.abs(a.s) || 0)) }; go.disabled = false; go.textContent = "Run sensitivity study"; render();
  }, "accent", "study-run");
  go.disabled = !m.studies.length;
  const controls = node("div", undefined, "rtl-actions"); controls.append(pick, go); sec.append(controls);
  sec.append(node("p", "Each input is raised 5% in a copy of the project and the whole linked system is re-solved. Elasticity is (Δy/y)/(Δx/x): 1 means proportional. Paths that cross domains appear here without being programmed anywhere.", "small muted"));
  // The inputs to study are part of the project.
  const inputs = node("ul", undefined, "study-inputs");
  for (const x of m.studies) { const li = node("li"); li.append(node("span", x.label), node("code", `${x.domain} ${x.path}`), button("×", () => quiet(run("project.unstudy", { id: x.id })), "unpin")); inputs.append(li); }
  if (!m.studies.length) inputs.append(node("li", "No inputs to study yet. Add any numeric field of a domain's document.", "small muted"));
  sec.append(inputs, fieldPaths("study-paths"));
  const form = miniForm(sec, { id: "study", action: "＋ Study input", cls: "study-add", fields: [{ name: "domain", label: "Domain", options: domains() }, { name: "path", label: "Field", list: "study-paths", placeholder: "battery.capacity", wide: true }, { name: "label", label: "Label" }, { name: "min", label: "Lowest", number: true }, { name: "max", label: "Highest", number: true }], submit: a => run("project.study", a) });
  linkPaths(form.domain, "study-paths");
  if (study) {
    const max = Math.max(...study.rows.map(r => Math.abs(r.s) || 0), 1e-9), list = node("div", undefined, "tornado");
    for (const r of study.rows) { const row = node("div"), bar = node("i"); bar.style.width = `${50 * Math.abs(r.s || 0) / max}%`; bar.dataset.sign = r.s >= 0 ? "pos" : "neg";
      row.append(node("span", r.label), node("b", Number.isFinite(r.s) ? r.s.toFixed(3) : r.error || "—"), node("div", undefined, "tornado-track")); row.lastChild.append(bar); list.append(row); }
    sec.append(node("p", `Baseline ${study.port} = ${format(study.y0, "", 5)}`, "small"), list);
  }
  return sec;
}
const domains = () => project.active().map(d => [d, DOMAIN_INFO[d].title]);
const outputs = s => [...s.ports.outputs].map(([id, p]) => [id, `${DOMAIN_INFO[p.domain]?.title} · ${p.label}`]);
// Studio re-renders keep an editor open once a person has opened it.
const opened = new Set();
const authoring = (root, title, fill, open) => { const d = node("details", undefined, "authoring"); d.open = !!open || opened.has(title); d.append(node("summary", title)); fill(d); d.addEventListener("toggle", () => { if (d.open) opened.add(title); else opened.delete(title); }); root.append(d); return d; };
const fieldPaths = id => Object.assign(node("datalist"), { id });
// Suggest every numeric field of the chosen domain's document.
const linkPaths = (select, id) => { const fill = () => { const list = document.getElementById(id); if (!list) return; list.replaceChildren(); const doc = select.value === "mech" ? mechanical.state() : project.store.docs[select.value]; if (doc) for (const p of numericPaths(doc)) list.append(Object.assign(node("option"), { value: p })); };
  select.addEventListener("change", fill); queueMicrotask(fill); };
function requirementsPanel(s) {
  const m = manifest(), sec = section(node("div"), "Requirements", `${m.requirements.length} of your own, checked on every solve`), ul = node("ul", undefined, "item-list");
  for (const q of m.requirements) { const li = node("li"), c = s.checks.find(x => x.domain === "project" && x.id === q.id), o = s.outputs[q.port];
    li.append(node("span", c?.text ?? q.port), node("code", q.port), node("b", format(o?.value, o?.unit, 4)), node("span", c?.status ?? "—", "pill"), button("Remove", () => quiet(run("project.unrequire", { id: q.id })))); ul.append(li); }
  if (!m.requirements.length) ul.append(node("li", "None yet. State what the design must do, on any output of any domain.", "small muted"));
  sec.append(ul); miniForm(sec, { id: "requirement", action: "＋ Requirement", fields: [{ name: "port", label: "Output", options: outputs(s), wide: true }, { name: "max", label: "At most", number: true }, { name: "min", label: "At least", number: true }, { name: "label", label: "Label", wide: true }], submit: a => run("project.require", a) });
  return sec;
}
function tourPanel(s) {
  const m = manifest(), sec = section(node("div"), "Guided tour", `${m.tour.length} step${m.tour.length === 1 ? "" : "s"}`), ul = node("ul", undefined, "item-list");
  m.tour.forEach((t, i) => { const li = node("li"); li.append(node("span", `${i + 1}. ${t.title}`), node("code", t.command || t.studio), button("Remove", () => quiet(run("project.tour_remove", { index: i + 1 })))); ul.append(li); });
  if (!m.tour.length) ul.append(node("li", "No tour yet. A tour shows others around this project with live commands.", "small muted"));
  sec.append(ul);
  authoring(sec, "Add a tour step", d => miniForm(d, { id: "tour-step-form", action: "＋ Step", fields: [{ name: "studio", label: "Studio", options: [["system", "System"], ...m.domains.map(x => [DOMAIN_INFO[x].studio, DOMAIN_INFO[x].title]), ["render", "Render"], ["history", "History"]].filter(([k]) => STUDIOS.includes(k)) }, { name: "title", label: "Title", wide: true }, { name: "text", label: "What to say", wide: true }, { name: "command", label: "Command to run", placeholder: "flight.mass gimbal 1.6", wide: true }, { name: "watch", label: "Outputs to report", placeholder: "flight.total_mass,mech.q1_tj", wide: true }, { name: "checks", label: "Checks to report", placeholder: "mech/REQ-TJ" }], submit: a => run("project.tour_step", a) }), !m.tour.length);
  return sec;
}
// Compass search over the project's variables. Feasible beats infeasible;
// between feasible trials the objective decides, between infeasible ones the
// smaller violation. Every trial is a full coupled solve on a copy.
function optimizer(spec, s) {
  const sec = node("section", undefined, "domain-section optimizer"), top = node("div", undefined, "section-top");
  const define = root => miniForm(root, { id: "optimize", action: spec ? "Update" : "＋ Define optimisation", fields: [{ name: "title", label: "Title", value: spec?.title, wide: true }, { name: "objective", label: "Objective", options: outputs(s), value: spec?.objective.port, wide: true }, { name: "sense", label: "Goal", options: [["min", "Minimise"], ["max", "Maximise"]], value: spec?.objective.sense }, { name: "trials", label: "Most trials", number: true, value: spec?.trials }, { name: "summary", label: "What it does", value: spec?.summary, wide: true }], submit: a => run("project.optimize", a) });
  if (!spec) { top.append(node("h2", "Optimisation"), node("span", "None defined")); sec.append(top, node("p", "Choose an output to minimise or maximise, then the fields it may change and the limits it must respect. Each trial is a full coupled solve.", "small muted")); define(sec); return sec; }
  top.append(node("h2", spec.title), node("span", `${spec.objective.sense === "min" ? "Minimise" : "Maximise"} ${spec.objective.port}`)); sec.append(top);
  sec.append(node("p", `${spec.summary ? spec.summary + " " : ""}Every linked model stays in its validated range.`, "small muted"));
  const vars = node("ul", undefined, "item-list");
  spec.variables.forEach((v, i) => { const li = node("li"); li.append(node("span", `${v.label} ${v.min}–${v.max}, step ${v.step}`), node("code", `${v.domain} ${v.path}`), button("Remove", () => quiet(run("project.optimize_remove", { kind: "variable", index: i + 1 })))); vars.append(li); });
  spec.constraints.forEach((c, i) => { const li = node("li"); li.append(node("span", c.label || `${c.port} ${c.max !== undefined ? "≤ " + c.max : ""} ${c.min !== undefined ? "≥ " + c.min : ""}`), node("code", c.port), button("Remove", () => quiet(run("project.optimize_remove", { kind: "constraint", index: i + 1 })))); vars.append(li); });
  if (!spec.variables.length) vars.append(node("li", "Add at least one variable to run it.", "small muted"));
  sec.append(vars);
  authoring(sec, "Edit this optimisation", d => {
    define(d); d.append(fieldPaths("optimize-paths"));
    const form = miniForm(d, { id: "optimize-variable", action: "＋ Variable", fields: [{ name: "domain", label: "Domain", options: domains() }, { name: "path", label: "Field", list: "optimize-paths", wide: true }, { name: "min", label: "From", number: true }, { name: "max", label: "To", number: true }, { name: "step", label: "Smallest step", number: true }, { name: "label", label: "Label" }, { name: "integer", label: "Whole numbers", check: true }], submit: a => run("project.optimize_variable", a) });
    linkPaths(form.domain, "optimize-paths");
    miniForm(d, { id: "optimize-constraint", action: "＋ Constraint", fields: [{ name: "port", label: "Output", options: outputs(s), wide: true }, { name: "max", label: "At most", number: true }, { name: "min", label: "At least", number: true }, { name: "label", label: "Label", wide: true }], submit: a => run("project.optimize_constraint", a) });
    d.append(button("Remove this optimisation", () => quiet(run("project.optimize_clear")), "", "optimize-clear"));
  }, !spec.variables.length);
  const go = button(`Run: ${spec.title}`, async () => {
    go.disabled = true; const base = docsNow(), vars = spec.variables, eps = 1e-9;
    const snap = (v, x) => { const k = Math.round((Math.min(v.max, Math.max(v.min, x)) - v.min) / v.step); return Math.min(v.max, v.min + k * v.step); };
    const evaluate = async values => {
      const docs = structuredClone(base);
      try { vars.forEach((v, i) => { setPath(docs[v.domain], v.path, values[i]); validate(docs, v.domain); }); } catch (e) { return { values, feasible: false, violation: Infinity, error: e.message }; }
      const r = await project.system.fork(docs).solve(docs.system), y = r.outputs[spec.objective.port]?.value;
      let violation = r.issues.some(i => i.kind === "range") ? 1 : 0;
      for (const c of spec.constraints) { const v = r.outputs[c.port]?.value; if (!Number.isFinite(v)) violation += 10; else { if (c.max !== undefined && v > c.max) violation += (v - c.max) / Math.max(1, Math.abs(c.max)); if (c.min !== undefined && v < c.min) violation += (c.min - v) / Math.max(1, Math.abs(c.min)); } }
      return { values, y, feasible: Number.isFinite(y) && violation === 0, violation, r };
    };
    const better = (a, b) => a.feasible !== b.feasible ? a.feasible : a.feasible ? (spec.objective.sense === "min" ? a.y < b.y - eps : a.y > b.y + eps) : a.violation < b.violation - eps;
    let best = await evaluate(vars.map(v => snap(v, getPath(base[v.domain], v.path)))), trials = 1;
    const steps = vars.map(v => Math.max(v.step, snap({ ...v, min: 0 }, (v.max - v.min) / 8)));
    while (trials < spec.trials && steps.some((st, i) => st >= vars[i].step)) {
      let improved = false;
      for (const [i, v] of vars.entries()) for (const dir of [-1, 1]) {
        if (steps[i] < v.step || trials >= spec.trials) continue;
        const x = snap(v, best.values[i] + dir * steps[i]); if (x === best.values[i]) continue;
        go.textContent = `Trial ${++trials}…`; const t = await evaluate(best.values.map((b, k) => k === i ? x : b));
        if (better(t, best)) { best = t; improved = true; }
      }
      if (!improved) steps.forEach((st, i) => { const half = Math.round(st / 2 / vars[i].step) * vars[i].step; steps[i] = half >= vars[i].step && half < st ? half : st > vars[i].step ? vars[i].step : 0; });
    }
    optimum = { ...best, trials, spec }; go.disabled = false; go.textContent = `Run: ${spec.title}`; render();
  }, "accent", "optimize-run");
  go.disabled = !spec.variables.length;
  sec.append(go);
  if (optimum && optimum.spec.title === spec.title) {
    const o = optimum, vals = spec.variables.map((v, i) => `${v.label} ${o.values[i]}`).join(", "), cons = spec.constraints.map(c => `${outputLabel(solution, c.port)} ${format(o.r?.outputs[c.port]?.value, o.r?.outputs[c.port]?.unit, 4)}`).join(", ");
    sec.append(node("p", o.feasible ? `Best after ${o.trials} coupled solves: ${vals} → ${outputLabel(solution, spec.objective.port)} ${format(o.y, o.r.outputs[spec.objective.port]?.unit, 5)}${cons ? `; ${cons}` : ""}.`
      : `No feasible design found in ${o.trials} coupled solves; the closest was ${vals}${o.error ? ` (${o.error})` : cons ? ` (${cons})` : ""}.`, "optimum"));
    if (o.feasible) sec.append(button("Apply as an undoable edit", () => project.bus.batch(spec.variables.map((v, i) => v.domain === "mech" ? { op: "mech.set", args: { name: v.path.slice(7), value: o.values[i] } } : { op: "project.set_value", args: { domain: v.domain, path: v.path, value: o.values[i] } }), { source: "optimizer", label: `Apply: ${spec.title}` })
      .then(() => toast("Optimum applied. Undo restores the previous design.")).catch(e => toast(e.message)), "", "optimize-apply"));
  }
  return sec;
}
onSolved(s => { solution = s; });
const refresh = live(["system", "project"], render);
document.addEventListener("emergence:studio", e => { if (e.detail === "system") render(); });
project.bus.on(() => { if (!root.hidden) refresh(); });
render();
