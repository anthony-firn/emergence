// Testing studio: what has been predicted, what has been checked against
// measurements, and what has actually been tested, kept apart.
import { project, onSolved } from "../app.js";
import { node, button, format, section, live, lineChart, run, miniForm } from "./kit.js";
const root = document.getElementById("testing-studio");
let solution = null;
const TITLE = { simulation: "Simulated", correlation: "Correlated with measurements", software: "Real software", physical: "Physical" };
const NOTE = { simulation: "Model predictions. They show the design is consistent with itself.", correlation: "Models against measured data they were not fitted to. They show the model can be trusted this far.", software: "Real flight code and real CAD tools run on this design. No hardware is involved, so these are not physical tests.", physical: "Measured on real hardware. Only these show that the thing works." };
const rows3 = text => text.trim().split(/\r?\n/).map(l => l.trim().split(/[\s,;]+/).map(Number)).filter(r => r.length === 3 && r.every(Number.isFinite));
const csv = text => text.trim().split(/\r?\n/).map(l => l.trim().split(/[\s,;]+/).map(Number)).filter(r => r.every(Number.isFinite) && r.length >= 2).map(([rpm, grams, torque]) => ({ rpm, grams, ...(torque !== undefined && { torque }) }));
function render() {
  const d = project.store.get("testing"), r = solution?.results.testing?.details, status = node("p", "", "studio-status"); status.id = "testing-status";
  root.replaceChildren();
  const head = node("div", undefined, "studio-heading"), title = node("div");
  title.append(node("p", `TESTING / ${project.store.docs.project.name.toUpperCase()}`, "panel-label"), node("h2", "What do we actually know?"),
    node("p", "Each result is kept as the kind of evidence it is. A simulation can never count as a test, and a physical test stays open until someone records a measurement on real hardware."));
  head.append(title); root.append(head);
  if (r) {
    const hero = node("div", undefined, "system-hero");
    for (const [k, s] of Object.entries(r.summary)) { const c = node("div", undefined, "hero-stat"); c.append(node("span", TITLE[k]), node("b", `${s.passed} of ${s.total} passed`), node("small", [s.failed && `${s.failed} failed`, s.stale && `${s.stale} stale`, `${s.open} open`].filter(Boolean).join(" · "))); c.dataset.kind = k; hero.append(c); }
    root.append(hero);
    if (r.summary.physical.total && r.summary.physical.passed < r.summary.physical.total) root.append(node("p", `${r.summary.physical.open} of ${r.summary.physical.total} physical tests have not been run. Nothing here has been proven on hardware until they are.`, "test-banner"));
    for (const kind of ["simulation", "correlation", "software", "physical"]) {
      const items = r.items.filter(i => i.kind === kind); if (!items.length) continue;
      const s = section(root, TITLE[kind], NOTE[kind]), table = node("table", undefined, "data-table test-table"), body = node("tbody"), th = node("tr"), head = node("thead");
      for (const h of ["Test", "Limit", "Result", "Status", ""]) th.append(node("th", h)); head.append(th); table.append(head);
      for (const i of items) { const row = node("tr"), name = node("td"); row.dataset.status = i.status; name.append(node("strong", i.title)); if (i.method) name.append(node("small", i.method)); if (i.detail) name.append(node("small", i.detail));
        const limit = [i.min !== undefined && `≥ ${i.min}`, i.max !== undefined && `≤ ${i.max}`].filter(Boolean).join(" ") || "—", del = node("td"); del.append(button("Remove", () => run("testing.remove", { id: i.id }, status).catch(() => {})));
        if (i.result) del.prepend(button("Clear", () => run("testing.clear_result", { id: i.id }, status).catch(() => {})));
        row.append(name, node("td", limit), node("td", i.value !== undefined ? format(i.value, i.unit, 4) : "—"), node("td", i.status), del); body.append(row); }
      table.append(body); s.append(table);
      if (kind === "physical" || kind === "software") miniForm(s, { id: `testing-record-${kind}`, action: "Record a measurement", status, fields: [{ name: "id", label: "Test", options: items.map(i => [i.id, i.title]) }, { name: "value", label: "Measured value", number: true }, { name: "date", label: "Date" }, { name: "by", label: "By" }, { name: "note", label: "Note", wide: true }], submit: a => run("testing.record", a, status) });
      for (const c of items.filter(i => i.correlation)) { const o = section(root, c.dataset.title, `model error ${c.correlation.rms.toFixed(2)}% RMS, worst ${c.correlation.worst.toFixed(2)}%`), pts = c.correlation.rows;
        o.append(lineChart({ series: [{ points: pts.map(p => [p.rpm, p.predicted]) }, { points: pts.map(p => [p.rpm, p.grams]), cls: "system" }], x: { label: "rotor speed, rpm" }, y: { label: "thrust, g", min: 0 }, marks: pts.filter((_, i) => i % 3 === 0).map(p => ({ x: p.rpm, y: p.grams })) }), node("p", "Line: the model fitted to the design's own table. Second line and marks: the independent measurements.", "small muted")); if (c.dataset.source) { const l = node("a", "Data source"); l.href = c.dataset.source; l.rel = "noopener"; o.append(l); } }
      if (kind === "software") for (const t of d.traces ?? []) { const o = section(root, t.title, t.note ?? ""), unit = t.unit ? `, ${t.unit}` : ""; o.append(lineChart({ series: [{ points: t.rows.map(r => [r[0], r[2]]) }, { points: t.rows.map(r => [r[0], r[1]]), cls: "system" }], x: { label: "time, s" }, y: { label: `response${unit}` } }), node("p", "Line: the model's prediction. Second line: the real run.", "small muted")); if (t.source) { const a = node("a", "Data source"); a.href = t.source; a.rel = "noopener"; o.append(a); } }
    }
  } else root.append(node("p", "Solving…", "small muted"));
  const add = section(root, "Add to the plan");
  miniForm(add, { id: "testing-add", action: "＋ Test", status, fields: [{ name: "id", label: "Id", placeholder: "hover-bench" }, { name: "title", label: "Title", wide: true }, { name: "kind", label: "Kind", options: [["physical", "Physical"], ["software", "Real software"], ["simulation", "Simulation"], ["correlation", "Correlation"]] }, { name: "covers", label: "Documents it depends on", placeholder: "controls, flight", when: v => v.kind === "software" || v.kind === "physical" }, { name: "port", label: "Output (simulation)", placeholder: "flight.hover_throttle", when: v => v.kind === "simulation" }, { name: "dataset", label: "Dataset", options: d.datasets.map(s => [s.id, s.title]), when: v => v.kind === "correlation" }, { name: "min", label: "At least", number: true }, { name: "max", label: "At most", number: true }, { name: "unit", label: "Unit" }, { name: "method", label: "Method", wide: true }], submit: a => run("testing.add", a, status) });
  miniForm(add, { id: "testing-dataset", action: "＋ Measured data", status, fields: [{ name: "id", label: "Id", placeholder: "my-bench" }, { name: "title", label: "Title", wide: true }, { name: "source", label: "Source (https)", wide: true }, { name: "rows", label: "Rows: rpm grams [torque], one per line", long: true, wide: true }], submit: a => { const rows = csv(a.rows ?? ""); if (rows.length < 2) throw new Error("Paste at least two rows of rpm and grams."); return run("testing.dataset", { ...a, rows }, status); } });
  if (d.datasets.length) { const ul = node("ul", undefined, "item-list"); for (const s of d.datasets) { const li = node("li"); li.append(node("span", `${s.title}: ${s.rows.length} rows`), button("Remove", () => run("testing.remove_dataset", { id: s.id }, status).catch(() => {}))); ul.append(li); } add.append(ul); }
  miniForm(add, { id: "testing-trace", action: "＋ Measured trace", status, fields: [{ name: "id", label: "Id", placeholder: "roll-step" }, { name: "title", label: "Title", wide: true }, { name: "unit", label: "Unit" }, { name: "note", label: "Note", wide: true }, { name: "source", label: "Source (https)", wide: true }, { name: "rows", label: "Rows: time measured model, one per line", long: true, wide: true }], submit: a => { const rows = rows3(a.rows ?? ""); if (rows.length < 2) throw new Error("Paste at least two rows of time, measured and model."); return run("testing.trace", { ...a, rows }, status); } });
  if (d.traces?.length) { const ul = node("ul", undefined, "item-list"); for (const t of d.traces) { const li = node("li"); li.append(node("span", `${t.title}: ${t.rows.length} rows`), button("Remove", () => run("testing.remove_trace", { id: t.id }, status).catch(() => {}))); ul.append(li); } add.append(ul); }
  add.append(status);
}
onSolved(s => { solution = s; });
live(["testing"], render);
render();
