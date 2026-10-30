// The page shell: the open project's heading and studio tabs, the project
// and domain menus, the command line, the global save bar, live status badges
// and the project's guided tour. Studios load as independent modules.
import { project, onSolved, onSolving, solution } from "./app.js";
import { DOMAIN_INFO, DOMAINS } from "./domains/index.js";
import { current, EXAMPLES, projects, open, openExample, newBlank, createFromFile, createFromBundle, remove } from "./projects.js";
import { readProjectFile, slug } from "./core/project.js";
import { unseal } from "./core/integrity.js";
import "./workspace.js";
import "./studios.js";
import "./ui/system-view.js";
import "./ui/flight-view.js";
import "./ui/parts-view.js";
import "./ui/optics-view.js";
import "./ui/fluids-view.js";
import "./ui/architecture-view.js";
import "./ui/render-view.js";
import { observe } from "./ui/tour.js";

const $ = id => document.getElementById(id);
const el = (tag, text, cls) => { const e = document.createElement(tag); if (text !== undefined) e.textContent = text; if (cls) e.className = cls; return e; };
const allTabs = [...document.querySelectorAll("[data-studio]")];
const domainOf = studio => studio === "mechanical" ? "mech" : studio;
const ALWAYS = new Set(["system", "render", "history"]);
// Only the project's domains have tabs.
const tabs = () => allTabs.filter(b => !b.hidden);
function syncTabs() {
  const on = new Set(project.active());
  for (const b of allTabs) b.hidden = !ALWAYS.has(b.dataset.studio) && !on.has(domainOf(b.dataset.studio));
  $("studio-add").textContent = on.size ? "＋ Domain" : "＋ Add a domain";
  if (!tabs().some(b => b.dataset.studio === activeStudio)) selectStudio("system");
}
export let activeStudio = "system";
export function selectStudio(name, focus = false) {
  if (!tabs().some(b => b.dataset.studio === name)) return;
  activeStudio = name;
  for (const b of allTabs) { const on = b.dataset.studio === name; b.setAttribute("aria-selected", String(on)); b.tabIndex = on ? 0 : -1; $(b.dataset.studio + "-studio").hidden = !on; if (on && focus) b.focus(); }
  document.dispatchEvent(new CustomEvent("emergence:studio", { detail: name }));
}
allTabs.forEach(b => {
  b.addEventListener("click", () => { selectStudio(b.dataset.studio); history.replaceState(null, "", "#" + b.dataset.studio); });
  b.addEventListener("keydown", e => {
    const list = tabs(), i = list.indexOf(b), n = e.key === "ArrowRight" ? (i + 1) % list.length : e.key === "ArrowLeft" ? (i + list.length - 1) % list.length : e.key === "Home" ? 0 : e.key === "End" ? list.length - 1 : -1;
    if (n >= 0) { e.preventDefault(); selectStudio(list[n].dataset.studio, true); history.replaceState(null, "", "#" + list[n].dataset.studio); }
  });
});
// The project heading, from its manifest.
function heading() {
  const m = project.store.docs.project, example = current.meta.source?.example;
  $("project-name").textContent = m.name; $("project-summary").textContent = m.summary || `${m.domains.length ? m.domains.map(d => DOMAIN_INFO[d].title).join(", ") : "No domains yet"}. Add domains, link them and solve them together.`;
  $("project-kind").textContent = example ? "EXAMPLE PROJECT" : "PROJECT"; $("project-slug").textContent = slug(m.name);
  $("start-tour").hidden = !m.tour.length; document.title = `${m.name} · Emergence CAD`;
}
project.store.on(({ ids }) => { if (ids.includes("project")) { syncTabs(); heading(); } });
syncTabs(); heading();
selectStudio(location.hash.slice(1) || "system");
addEventListener("hashchange", () => selectStudio(location.hash.slice(1)));

// Live status per domain on the tabs.
onSolved(s => {
  const by = {};
  for (const c of s.checks) (by[c.domain] ||= []).push(c.status);
  for (const [d, list] of Object.entries(by)) {
    const badge = $(`badge-${d === "mech" ? "mechanical" : d}`); if (!badge) continue;
    const bad = list.filter(x => x === "unmet" || x === "unavailable").length, open = list.filter(x => x === "undecided").length;
    badge.textContent = bad ? `✗${bad}` : open ? `◌${open}` : `✓${list.length}`; badge.dataset.state = bad ? "bad" : open ? "open" : "ok";
    badge.title = `${list.length} checks: ${bad} failing or unavailable, ${open} undecided`;
  }
  const sys = $("badge-system"); sys.textContent = s.converged ? `${s.iterations}×` : "≁"; sys.dataset.state = s.converged ? "ok" : "bad"; sys.title = `${s.iterations} coupled passes, ${s.converged ? "settled" : "not settled"}`;
});
onSolving.add(busy => document.body.classList.toggle("solving", !!busy));
onSolved(() => document.body.classList.remove("solving"));

// Global revision bar: one save covers every domain.
document.addEventListener("emergence:draft", ({ detail }) => {
  $("global-save").disabled = !detail.dirty;
  $("global-draft").textContent = detail.dirty ? `Unsaved draft${detail.studios.length ? " · " + detail.studios.join(", ") : ""}` : "No uncommitted changes";
});
$("global-save").addEventListener("click", () => { $("save-summary").textContent = $("draft-status").textContent; $("save-message").value = $("message").value; $("save-dialog").showModal(); $("save-message").focus(); });
$("save-form").addEventListener("submit", e => { e.preventDefault(); $("message").value = $("save-message").value; $("save-dialog").close(); $("commit").click(); });
$("global-reviews").addEventListener("click", () => { selectStudio("history"); $("tab-review").click(); });

// Projects: examples, this device's projects, a blank one, or a file.
const card = (title, text, actions, cls = "project-card") => { const c = el("article", undefined, cls); c.append(el("h4", title), el("p", text)); const row = el("div", undefined, "project-card-actions"); row.append(...actions); c.append(row); return c; };
const act = (text, fn, cls) => { const b = el("button", text, cls); b.type = "button"; b.addEventListener("click", fn); return b; };
function renderProjects() {
  const ex = $("project-examples"); ex.replaceChildren();
  for (const e of EXAMPLES) {
    const mine = projects().find(p => p.source?.example === e.id), here = current.meta.source?.example === e.id;
    ex.append(card(e.name, e.summary, [here ? el("span", "Open now", "pill") : act(mine ? "Open your copy" : "Open example", () => openExample(e.id), "accent")]));
  }
  const list = $("project-list"); list.replaceChildren();
  for (const p of projects()) {
    const row = el("div", undefined, "project-row"), when = new Date(p.opened || p.created);
    row.append(el("strong", p.name), el("span", `${p.source?.example ? "from an example · " : p.source?.legacy ? "from before projects · " : ""}opened ${when.toLocaleDateString()}`, "small muted"));
    if (p.id === current.id) row.append(el("span", "Open now", "pill"));
    else row.append(act("Open", () => open(p.id)), act("Delete", () => { if (confirm(`Delete “${p.name}” and its revision history from this device? Export its history first to keep a copy.`)) { remove(p.id); renderProjects(); } }));
    list.append(row);
  }
}
$("projects-open").addEventListener("click", () => { renderProjects(); $("projects-dialog").showModal(); });
$("project-new").addEventListener("click", () => newBlank());
$("project-download").addEventListener("click", () => {
  const file = window.emergence.file(), url = URL.createObjectURL(new Blob([JSON.stringify(file, null, 1) + "\n"], { type: "application/json" })), a = el("a");
  a.href = url; a.download = `${slug(file.documents.project.name)}.emergence.json`; a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
});
$("project-file").addEventListener("change", async () => {
  const f = $("project-file").files[0]; $("project-file").value = ""; if (!f) return;
  try {
    if (f.size > 3_000_000) throw new Error("Project files are limited to 3 MB here.");
    const check = await unseal(JSON.parse(await f.text())); if (!check.ok) throw new Error("Its SHA-256 does not match its content.");
    const body = check.body;
    if (body.format === "emergence-cad/1") { const r = createFromBundle(body); open(r.id); return; }
    readProjectFile(body);
    // Every document must pass its own domain's validation before the project exists.
    for (const id of project.store.ids()) if (body.documents[id] !== undefined) project.store.domains[id].validate(structuredClone(body.documents[id]));
    open(createFromFile(body));
  } catch (e) { out.textContent = `✗ ${f.name} was not opened: ${e.message}`; out.dataset.state = "error"; $("projects-dialog").close(); }
});
// Domains: add or remove studios from this project.
function renderDomains() {
  const root = $("domain-choices"), on = project.active(); root.replaceChildren();
  for (const d of DOMAINS) {
    const info = DOMAIN_INFO[d], has = on.includes(d);
    const c = card(`${info.icon} ${info.title}`, info.blurb, [act(has ? "Remove" : "Add", () => project.bus.run(has ? "project.remove_domain" : "project.add_domain", { domain: d }).then(() => { renderDomains(); if (!has) selectStudio(info.studio); }).catch(e => { out.textContent = `✗ ${e.message}`; }), has ? "" : "accent")]);
    c.dataset.on = has; root.append(c);
  }
}
$("studio-add").addEventListener("click", () => { renderDomains(); $("domain-dialog").showModal(); });

// The command line: the registry's parser, completions and the shared undo.
const input = $("command-input"), out = $("command-output"), list = $("command-suggestions"), past = []; let cursor = 0, picks = [], pick = -1;
const context = () => domainOf(activeStudio);
function suggest() {
  picks = input.value.trim() ? project.bus.complete(input.value, context()) : []; pick = -1; list.replaceChildren();
  for (const [i, s] of picks.entries()) { const li = el("li"); li.setAttribute("role", "option"); li.id = `cmd-opt-${i}`; li.append(el("b", s.text), el("span", s.hint || "")); li.addEventListener("mousedown", e => { e.preventDefault(); accept(i); }); list.append(li); }
  list.hidden = !picks.length;
}
function accept(i) {
  const s = picks[i]; if (!s) return;
  const tokens = input.value.split(/\s+/);
  input.value = s.insert ? input.value.replace(/\S*$/, "") + s.text : s.text + " ";
  if (!s.insert && tokens.length > 1) input.value = s.text + " " + tokens.slice(1).join(" ");
  input.focus(); suggest();
}
input.addEventListener("input", suggest);
input.addEventListener("blur", () => setTimeout(() => list.hidden = true, 120));
input.addEventListener("keydown", e => {
  if (e.key === "Tab" && picks.length) { e.preventDefault(); accept(pick >= 0 ? pick : 0); }
  else if (e.key === "ArrowDown" && !list.hidden) { e.preventDefault(); pick = (pick + 1) % picks.length; highlight(); }
  else if (e.key === "ArrowUp" && !list.hidden && pick >= 0) { e.preventDefault(); pick = pick ? pick - 1 : picks.length - 1; highlight(); }
  else if (e.key === "ArrowUp" && past.length) { e.preventDefault(); cursor = Math.max(0, cursor - 1); input.value = past[cursor] || ""; }
  else if (e.key === "ArrowDown" && past.length && list.hidden) { e.preventDefault(); cursor = Math.min(past.length, cursor + 1); input.value = past[cursor] || ""; }
  else if (e.key === "Escape") { list.hidden = true; input.blur(); }
});
function highlight() { [...list.children].forEach((li, i) => li.setAttribute("aria-selected", String(i === pick))); input.setAttribute("aria-activedescendant", pick >= 0 ? `cmd-opt-${pick}` : ""); }
$("command-form").addEventListener("submit", async e => {
  e.preventDefault(); if (pick >= 0 && !list.hidden) { accept(pick); return; }
  const text = input.value.trim(); if (!text) return;
  past.push(text); cursor = past.length; list.hidden = true;
  const [head, ...rest] = text.split(/\s+/);
  if (head === "help" || head === "?") {
    if (!rest.length) { openReference(); out.textContent = `${project.bus.commands.size} commands. Type a domain prefix (project., pcb., optics.…) for completions.`; return; }
    try { const spec = project.bus.resolve(rest[0], context()); out.textContent = `${project.bus.usage(spec.op)} — ${spec.summary || spec.title}${spec.example ? ` · e.g. ${spec.example}` : ""}`; } catch (err) { out.textContent = err.message; }
    input.value = ""; return;
  }
  if (head === "undo" || head === "redo") { out.textContent = project.bus[head]() ? `${head === "undo" ? "Undid" : "Redid"} the last edit.` : `Nothing to ${head}.`; input.value = ""; return; }
  if (head === "solve") { out.textContent = "Solving…"; const s = await project.solve(); out.textContent = `${s.iterations} passes, residual ${s.residual.toExponential(1)}, ${s.ms.toFixed(0)} ms.`; input.value = ""; return; }
  if (head === "script") { navigator.clipboard?.writeText(JSON.stringify(project.bus.script(), null, 2)).catch(() => {}); out.textContent = `Copied ${project.bus.script().length} edits as a JSON script (also in the System studio).`; input.value = ""; return; }
  try {
    const { op, args } = project.bus.parse(text, context()), r = await project.bus.run(op, args, { source: "command line" });
    out.textContent = `✓ ${op}${r?.message ? " · " + r.message : ""}`; out.dataset.state = "ok"; input.value = "";
    const studio = op.split(".")[0] === "mech" ? "mechanical" : op.split(".")[0] === "litho" ? "silicon" : op.split(".")[0];
    if (studio !== activeStudio && tabs().some(t => t.dataset.studio === studio) && studio !== "system") out.textContent += ` · see the ${DOMAIN_INFO[domainOf(studio)]?.title || studio} studio`;
  } catch (err) { out.textContent = `✗ ${err.message}`; out.dataset.state = "error"; }
});
addEventListener("keydown", e => {
  if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === "k") { e.preventDefault(); input.focus(); input.select(); }
  const typing = /INPUT|TEXTAREA|SELECT/.test(document.activeElement?.tagName) || document.activeElement?.isContentEditable;
  if (!typing && (e.ctrlKey || e.metaKey) && e.key.toLowerCase() === "z") { e.preventDefault(); e.shiftKey ? project.bus.redo() : project.bus.undo(); }
});
$("command-undo").addEventListener("click", () => project.bus.undo());
$("command-redo").addEventListener("click", () => project.bus.redo());
const syncUndo = () => { $("command-undo").disabled = !project.bus.undoStack.length; $("command-redo").disabled = !project.bus.redoStack.length; $("command-undo").title = project.bus.undoStack.at(-1) ? `Undo: ${project.bus.undoStack.at(-1).label}` : "Nothing to undo"; };
project.bus.on(syncUndo); syncUndo();
function openReference() {
  const render = () => {
    const f = $("commands-filter").value.toLowerCase(), root = $("commands-list"); root.replaceChildren();
    const groups = {}; for (const c of project.bus.list()) if (!f || (c.op + c.title + (c.summary || "")).toLowerCase().includes(f)) (groups[c.domain] ||= []).push(c);
    for (const [d, cs] of Object.entries(groups)) {
      const sec = el("section"); sec.append(el("h3", `${DOMAIN_INFO[d]?.icon || "◇"} ${DOMAIN_INFO[d]?.title || d}`));
      for (const c of cs) { const row = el("div", undefined, "command-row"); row.append(el("code", project.bus.usage(c.op)), el("span", c.summary || c.title)); if (c.example) row.append(el("small", "e.g. " + c.example));
        row.addEventListener("click", () => { input.value = (c.example || c.op) + (c.example ? "" : " "); $("commands-dialog").close(); input.focus(); suggest(); }); sec.append(row); }
      root.append(sec);
    }
  };
  $("commands-filter").oninput = render; render(); $("commands-dialog").showModal();
}
$("command-reference").addEventListener("click", openReference);
for (const b of document.querySelectorAll("[data-close]")) if (!b.dataset.bound) { b.dataset.bound = 1; b.addEventListener("click", () => $(b.dataset.close).close()); }

// The project's guided tour: each step is a real command against the live
// model, with the design restored afterwards.
let tourIndex = 0, tourSnapshot = null, tourUndo = 0, tourRan = false, TOUR = [];
const ids = () => ["mech", ...project.store.ids()];
function showStep() {
  const step = TOUR[tourIndex];
  $("tour-step").textContent = `STEP ${tourIndex + 1} OF ${TOUR.length} · ${DOMAIN_INFO[domainOf(step.studio)]?.title || step.studio[0].toUpperCase() + step.studio.slice(1)}`;
  $("tour-title").textContent = step.title; $("tour-text").textContent = step.text; $("tour-command").textContent = step.command || "";
  $("tour-command").hidden = !step.command; $("tour-result").textContent = "";
  $("tour-next").textContent = step.command ? "Run step" : tourIndex === TOUR.length - 1 ? "Finish" : "Next";
  $("tour-back").disabled = tourIndex === 0; selectStudio(step.studio);
}
$("start-tour").addEventListener("click", () => {
  TOUR = structuredClone(project.store.docs.project.tour); if (!TOUR.length) return;
  tourSnapshot = Object.fromEntries(ids().map(d => [d, project.bus.domains.get(d).snapshot()])); tourUndo = project.bus.undoStack.length;
  tourIndex = 0; tourRan = false; $("tour").hidden = false; showStep();
});
$("tour-next").addEventListener("click", async () => {
  const step = TOUR[tourIndex];
  if (step.command && !tourRan) {
    $("tour-next").disabled = true; $("tour-result").textContent = "Running and re-solving the linked system…";
    try {
      const before = solution;
      const { op, args } = project.bus.parse(step.command); await project.bus.run(op, args, { source: "tour" });
      const after = await new Promise(resolve => { const off = onSolved(s => { if (s !== before) { off(); resolve(s); } }); });
      $("tour-result").textContent = observe(step, after, before);
      tourRan = true; $("tour-next").textContent = tourIndex === TOUR.length - 1 ? "Finish" : "Next";
    } catch (e) { $("tour-result").textContent = e.message; }
    finally { $("tour-next").disabled = false; }
    return;
  }
  if (tourIndex === TOUR.length - 1) { endTour(); return; }
  tourIndex++; tourRan = false; showStep();
});
$("tour-back").addEventListener("click", () => { if (tourIndex) { tourIndex--; tourRan = false; showStep(); } });
function endTour() {
  if (tourSnapshot) { for (const [d, s] of Object.entries(tourSnapshot)) project.bus.domains.get(d).restore(s); project.bus.undoStack.length = Math.min(tourUndo, project.bus.undoStack.length); project.bus.redoStack = []; project.bus.emit({ type: "reset" }); }
  tourSnapshot = null; $("tour").hidden = true; selectStudio("system");
}
$("tour-close").addEventListener("click", endTour);
