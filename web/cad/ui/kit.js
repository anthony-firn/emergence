// Small DOM helpers shared by the domain studios. Every editing control runs
// a registry command, so the UI is one client of the same API agents use.
import { project, onSolved } from "../app.js";
import { format, label as unitLabel } from "../core/units.js";
export { format };
export const node = (tag, text, cls) => { const e = document.createElement(tag); if (text !== undefined && text !== null) e.textContent = text; if (cls) e.className = cls; return e; };
export const svg = (tag, attrs = {}, text) => { const e = document.createElementNS("http://www.w3.org/2000/svg", tag); for (const [k, v] of Object.entries(attrs)) if (v !== undefined) e.setAttribute(k, v); if (text !== undefined) e.textContent = text; return e; };
export function button(text, fn, cls, id) { const b = node("button", text, cls); b.type = "button"; if (id) b.id = id; b.addEventListener("click", fn); return b; }
export function download(data, name, type) { const url = URL.createObjectURL(new Blob([data], { type })), a = node("a"); a.href = url; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(url), 1500); }
let toastTimer;
export function toast(message) { const t = document.getElementById("toast"); if (!t) return; t.textContent = message; t.hidden = false; clearTimeout(toastTimer); toastTimer = setTimeout(() => t.hidden = true, 6500); }
// Run a command; report its message or error in the given status element.
export function run(op, args, status) {
  return project.bus.run(op, args).then(r => { if (status && r?.message) status.textContent = r.message; return r; })
    .catch(e => { if (status) status.textContent = e.message; else toast(e.message); throw e; });
}
const quiet = p => p.catch(() => {});
export function numberField(root, { label, id, value, min, max, step = "any", unit, op, args, status, hint }) {
  const l = node("label", label), i = node("input"); i.type = "number"; if (id) { i.id = id; l.htmlFor = id; }
  if (min !== undefined) i.min = min; if (max !== undefined) i.max = max; i.step = step; i.value = Number.isFinite(value) ? value : "";
  if (unit) l.append(node("small", unitLabel(unit) || unit, "field-unit"));
  l.append(i); if (hint) l.append(node("small", hint, "field-hint")); root.append(l);
  i.addEventListener("change", () => {
    if (i.value === "" || !i.checkValidity()) { i.value = value; if (status) status.textContent = `Enter ${min ?? "−∞"}–${max ?? "∞"}${unit ? " " + (unitLabel(unit) || unit) : ""} for ${label}.`; return; }
    quiet(run(op, args(Number(i.value)), status).catch(() => { i.value = value; }));
  });
  return i;
}
export function textField(root, { label, id, value, max = 60, op, args, status }) {
  const l = node("label", label), i = node("input"); if (id) { i.id = id; l.htmlFor = id; } i.value = value; i.maxLength = max; l.append(i); root.append(l);
  i.addEventListener("change", () => quiet(run(op, args(i.value), status).catch(() => { i.value = value; })));
  return i;
}
export function selectField(root, { label, id, value, options, op, args, status }) {
  const l = node("label", label), s = node("select"); if (id) { s.id = id; l.htmlFor = id; }
  for (const [v, t] of options) { const o = node("option", t); o.value = v; s.append(o); }
  s.value = value; l.append(s); root.append(l);
  s.addEventListener("change", () => quiet(run(op, args(s.value), status).catch(() => { s.value = value; })));
  return s;
}
// A compact form that runs one command. Empty fields are left out, so the
// command's own defaults apply; a failed command keeps what was typed.
// What has been typed outlives any re-render, until the form is submitted.
const drafts = new Map();
export function miniForm(root, { id, fields, action, submit, status, cls = "" }) {
  const form = node("div", undefined, `mini-form ${cls}`), inputs = {}, labels = {}, key = f => `${id}-${f.name}`;
  for (const f of fields) {
    const l = node("label", f.label), i = node(f.options ? "select" : f.long ? "textarea" : "input"); i.id = `${id}-${f.name}`; if (f.long) i.rows = 5; l.htmlFor = i.id; if (f.wide) l.classList.add("wide");
    if (f.options) for (const [v, t] of f.options) { const o = node("option", t); o.value = v; i.append(o); }
    else { if (f.number) { i.type = "number"; i.step = f.step ?? "any"; } if (f.check) i.type = "checkbox"; if (f.placeholder) i.placeholder = f.placeholder; if (f.list) i.setAttribute("list", f.list); i.spellcheck = false; }
    if (f.value !== undefined) i[f.check ? "checked" : "value"] = f.value;
    if (drafts.has(key(f))) i[f.check ? "checked" : "value"] = drafts.get(key(f));
    i.addEventListener(f.options || f.check ? "change" : "input", () => drafts.set(key(f), f.check ? i.checked : i.value));
    if (f.check) l.classList.add("check-label");
    l.append(i); form.append(l); inputs[f.name] = i; labels[f.name] = l;
  }
  // A field with `when` shows only while it applies to the form's other values.
  const current = () => Object.fromEntries(fields.map(f => [f.name, inputs[f.name].value]));
  const refresh = () => { for (const f of fields) if (f.when) labels[f.name].hidden = !f.when(current()); };
  for (const f of fields) if (f.options) inputs[f.name].addEventListener("change", refresh);
  refresh();
  form.append(button(action, () => {
    const values = {};
    for (const f of fields) { const i = inputs[f.name]; if (labels[f.name].hidden) continue; if (f.check) { if (i.checked) values[f.name] = true; } else if (i.value !== "") values[f.name] = f.number ? Number(i.value) : i.value; }
    Promise.resolve().then(() => submit(values)).then(() => { for (const f of fields) drafts.delete(key(f)); for (const f of fields) if (f.value === undefined && !f.options) inputs[f.name][f.check ? "checked" : "value"] = f.check ? false : ""; })
      .catch(e => { if (status) status.textContent = e.message; else toast(e.message); });
  }, "accent", `${id}-submit`));
  root.append(form); return inputs;
}
export function checkCards(root, checks) {
  const grid = node("div", undefined, "requirements domain-checks");
  for (const c of checks) {
    const card = node("article", undefined, "requirement"); card.dataset.status = c.status;
    const top = node("div", undefined, "requirement-top"); top.append(node("h3", c.text), node("span", c.status, "status-badge")); card.append(top);
    if (c.value?.nominal !== undefined && Number.isFinite(c.value.nominal)) card.append(node("p", format(c.value.nominal, c.unit, 4), "value"));
    if (c.reason) card.append(node("p", c.reason, "reason"));
    grid.append(card);
  }
  if (!checks.length) grid.append(node("p", "No checks for this domain yet.", "small muted"));
  root.append(grid); return grid;
}
export function certificates(root, list) {
  if (!list?.length) return;
  const ul = node("ul", undefined, "certificates");
  for (const c of list) { const li = node("li", undefined, c.ok ? "cert-ok" : "cert-fail"); li.append(node("b", c.ok ? "✓" : "✗"), node("span", c.label)); ul.append(li); }
  root.append(ul);
}
export function assumptionList(root, domain) {
  const items = project.assumptions().filter(a => !domain || a.domain === domain);
  const ul = node("ul", undefined, "assumption-list");
  for (const a of items) {
    const li = node("li"); li.dataset.kind = a.kind; li.dataset.confirmed = a.confirmed;
    const c = node("input"); c.type = "checkbox"; c.checked = a.confirmed; c.setAttribute("aria-label", `I have checked: ${a.label}`);
    c.addEventListener("change", () => quiet(run("system.confirm", { id: a.id, confirmed: c.checked })));
    const text = node("div"); text.append(node("strong", a.label), node("small", `${a.kind} · ${a.source}`));
    li.append(c, text); ul.append(li);
  }
  if (!items.length) ul.append(node("li", "No stated assumptions.", "small muted"));
  root.append(ul); return ul;
}
export function section(root, title, aside) {
  const s = node("section", undefined, "domain-section"), top = node("div", undefined, "section-top");
  top.append(node("h2", title)); if (aside) top.append(typeof aside === "string" ? node("span", aside) : aside); s.append(top); root.append(s); return s;
}
// Linked inputs: a small "driven by" chip with an unlink action.
export function linkedChip(solution, port) {
  const a = solution?.applied[port]; if (!a) return null;
  const chip = node("p", undefined, "link-note");
  chip.append(node("span", `⛓ ${format(a.value, a.unit, 4)} from ${a.from.join(" + ")}`), button("Unlink", () => quiet(Promise.all(a.links.map(id => run("system.toggle", { id, enabled: false }))))));
  return chip;
}
// Minimal line chart: series [{ points: [[x, y]], cls }], with axes labels.
export function lineChart({ series, x, y, width = 420, height = 220, marks = [] }) {
  const pad = { l: 46, r: 12, t: 12, b: 34 }, all = series.flatMap(s => s.points).concat(marks.map(m => [m.x, m.y])).filter(p => p.every(Number.isFinite));
  const x0 = x.min ?? Math.min(...all.map(p => p[0])), x1 = x.max ?? Math.max(...all.map(p => p[0])), y0 = y.min ?? Math.min(0, ...all.map(p => p[1])), y1 = y.max ?? (Math.max(...all.map(p => p[1])) * 1.08 || 1);
  const sx = v => pad.l + (v - x0) / (x1 - x0 || 1) * (width - pad.l - pad.r), sy = v => height - pad.b - (v - y0) / (y1 - y0 || 1) * (height - pad.t - pad.b);
  const root = svg("svg", { viewBox: `0 0 ${width} ${height}`, class: "chart", role: "img", "aria-label": `${y.label} against ${x.label}` });
  for (let i = 0; i <= 4; i++) { const v = y0 + (y1 - y0) * i / 4; root.append(svg("line", { x1: pad.l, x2: width - pad.r, y1: sy(v), y2: sy(v), class: "chart-grid" }), svg("text", { x: pad.l - 6, y: sy(v) + 3, "text-anchor": "end", class: "chart-tick" }, Number(v.toPrecision(3)))); }
  for (let i = 0; i <= 4; i++) { const v = x0 + (x1 - x0) * i / 4; root.append(svg("text", { x: sx(v), y: height - pad.b + 14, "text-anchor": "middle", class: "chart-tick" }, Number(v.toPrecision(3)))); }
  root.append(svg("text", { x: (pad.l + width - pad.r) / 2, y: height - 4, "text-anchor": "middle", class: "chart-label" }, x.label), svg("text", { x: 10, y: pad.t + 2, class: "chart-label" }, y.label));
  for (const s of series) { const pts = s.points.filter(p => p.every(Number.isFinite)); if (pts.length > 1) root.append(svg("polyline", { points: pts.map(p => `${sx(p[0])},${sy(p[1])}`).join(" "), class: `chart-line ${s.cls || ""}` })); }
  for (const m of marks) if ([m.x, m.y].every(Number.isFinite)) root.append(svg("circle", { cx: sx(m.x), cy: sy(m.y), r: 4, class: "chart-mark" }), svg("text", { x: sx(m.x) + 7, y: sy(m.y) - 7, class: "chart-tick" }, m.label || ""));
  return root;
}
// Re-render a studio when its documents change or a system solve lands, but
// never under someone's fingers: a solve can land while they are typing.
const typing = () => document.activeElement?.matches?.("input:not([type=checkbox]):not([type=range]),textarea,select");
export function live(domains, render) {
  let pending = false;
  const go = () => { if (pending) return; pending = true; requestAnimationFrame(() => { pending = false;
    if (typing()) { document.activeElement.addEventListener("focusout", go, { once: true }); return; }
    render(); }); };
  project.store.on(({ ids }) => { if (ids.some(id => domains.includes(id))) go(); });
  onSolved(go);
  return go;
}
