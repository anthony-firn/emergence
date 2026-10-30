// Parts studio: solid components as feature trees built by the Rust kernel.
import { project, onSolved } from "../app.js";
import { TEMPLATES, MATERIALS, instances } from "../domains/parts.js";
import { Renderer, M } from "./gl.js";
import { modelMesh, stl } from "./meshes.js";
import { node, button, format, numberField, selectField, textField, assumptionList, section, run, download, toast } from "./kit.js";
const root = document.getElementById("parts-studio");
let solution = null, selected = null, viewer = null;
const gpu = new Map(), canvas = node("canvas", undefined, "part-canvas");
canvas.id = "parts-canvas"; canvas.tabIndex = 0; canvas.setAttribute("aria-label", "Part preview. Drag to orbit, scroll to zoom, shift-drag to pan.");
async function preview(p, info) {
  try {
    viewer ||= new Renderer(canvas);
    info.textContent = "Meshing in the Rust kernel…";
    const t0 = performance.now(), mesh = await modelMesh(p.model, p.overrides);
    if (selected !== p.id) return;
    const g = gpu.get(mesh) || viewer.upload(mesh); gpu.set(mesh, g);
    const [lo, hi] = mesh.bbox, colour = { cfrp: [.22, .23, .25, 1], "al-6061": [.72, .76, .8, 1], abs: [.28, .36, .44, 1], brass: [.78, .62, .3, 1], "steel-304": [.62, .64, .66, 1] }[p.material] || [.6, .65, .7, 1];
    viewer.set([{ gpu: g, matrix: M.identity(), color: colour, part: p }], { size: Math.max(...hi.map((v, i) => v - lo[i])), step: Math.max(...hi.map((v, i) => v - lo[i])) / 10, z: lo[2] });
    viewer.fit(lo, hi);
    const mass = mesh.volume * MATERIALS[p.material].density / 1000;
    info.textContent = `${mesh.triangle_count.toLocaleString()} triangles · ${mesh.watertight ? "watertight" : "NOT watertight"} · cell ${mesh.cell.toFixed(3)} mm${mesh.coarsened ? " (coarsened for the browser budget)" : ""} · volume ${format(mesh.volume, "mm3", 5)} · ${format(mass, "g", 4)} each × ${p.count} = ${format(mass * p.count, "g", 4)} · ${(performance.now() - t0).toFixed(0)} ms`;
    info.dataset.mesh = "ready";
  } catch (e) { info.textContent = `Mesh unavailable: ${e.message}`; info.dataset.mesh = "error"; }
}
function render() {
  const d = project.store.get("parts"), status = node("p", "", "studio-status"); status.id = "parts-status";
  if (!d.parts.some(p => p.id === selected)) selected = d.parts[0]?.id ?? null;
  root.replaceChildren();
  const head = node("div", undefined, "studio-heading"), title = node("div");
  title.append(node("p", `PARTS / ${project.store.docs.project.name.toUpperCase()}`, "panel-label"), node("h2", "Feature trees on a field kernel."),
    node("p", "Each part is a JSON feature tree in the fluid cad format: primitives, booleans, shells, arrays and kits (gears, threads, sheet metal, pipe runs, lattices, wings, lenses), where every number may be an expression of the part's parameters. The Rust kernel builds an exact signed-distance field and meshes it watertight in the browser."));
  const tools = node("div", undefined, "studio-file-tools"), tpl = node("select"); tpl.id = "parts-template"; tpl.setAttribute("aria-label", "Template");
  for (const [k, t] of Object.entries(TEMPLATES)) { const o = node("option", t.title); o.value = k; tpl.append(o); }
  tools.append(tpl, button("＋ Add part", () => run("parts.add", { template: tpl.value }, status).then(r => { selected = r.id; render(); }).catch(() => {}), "accent", "parts-add"));
  head.append(title, tools); root.append(head);
  const vols = solution?.results.parts?.details?.volumes || {}, total = solution?.outputs["parts.mass"]?.value;
  const grid = node("div", undefined, "parts-grid"), list = node("aside", undefined, "parts-list");
  list.append(node("p", `PARTS · ${total !== undefined ? format(total, "kg", 4) + " total" : "…"}`, "panel-label"));
  for (const p of d.parts) {
    const b = button(`${p.name} ×${p.count}`, () => { selected = p.id; render(); }, "studio-tree-item"); b.setAttribute("aria-pressed", String(selected === p.id));
    const m = vols[p.id] !== undefined ? vols[p.id] * MATERIALS[p.material].density / 1000 * p.count : null;
    b.append(node("small", ` · ${MATERIALS[p.material].title}${m !== null ? ` · ${format(m, "g", 3)}` : ""}`)); list.append(b);
  }
  const centre = node("div", undefined, "parts-centre");
  const info = node("p", "", "mesh-info"); info.id = "parts-mesh-info"; centre.append(canvas, info);
  const insp = node("aside", undefined, "parts-inspector"), p = d.parts.find(x => x.id === selected);
  if (p) {
    insp.append(node("p", "PART", "panel-label"));
    textField(insp, { label: "Name", id: "parts-name", value: p.name, op: "parts.update", args: v => ({ id: p.id, name: v }), status });
    selectField(insp, { label: "Material", id: "parts-material", value: p.material, options: Object.entries(MATERIALS).map(([k, m]) => [k, `${m.title} · ${m.density} g/cm³`]), op: "parts.update", args: v => ({ id: p.id, material: v }), status });
    insp.append(node("p", "PARAMETERS", "panel-label studio-section-label"));
    for (const [k, v] of Object.entries(p.model.params)) {
      const value = p.overrides[k] ?? v;
      if (typeof value === "number") numberField(insp, { label: k, id: `parts-param-${k}`, value, step: "any", op: "parts.param", args: x => ({ id: p.id, name: k, value: x }), status });
      else insp.append(node("p", `${k} = ${value}`, "expr-param"));
    }
    insp.append(node("p", `PLACEMENT · ${p.pattern === "polar" ? `polar ×${p.count} about z` : p.pattern === "stack" ? `stack ×${p.count}` : "single"}`, "panel-label studio-section-label"));
    for (const k of ["x", "y", "z"]) numberField(insp, { label: k.toUpperCase(), id: `parts-${k}`, unit: "mm", value: p[k], min: -2000, max: 2000, step: .5, op: "parts.update", args: v => ({ id: p.id, [k]: v }), status });
    const tree = node("details", undefined, "tree-editor"); tree.append(node("summary", "Feature tree (JSON)"));
    const area = node("textarea"); area.id = "parts-tree"; area.rows = 12; area.spellcheck = false; area.value = JSON.stringify(p.model, null, 1); tree.append(area);
    tree.append(button("Apply feature tree", () => { try { run("parts.model", { id: p.id, model: JSON.parse(area.value) }, status).catch(() => {}); } catch (e) { status.textContent = `JSON: ${e.message}`; } }, "accent", "parts-apply-tree"));
    insp.append(tree, button("↓ STL of this part", async () => { try { const m = await modelMesh(p.model, p.overrides); download(stl(m, p.name), `${p.id}-mm.stl`, "model/stl"); } catch (e) { toast(e.message); } }, "", "parts-stl"),
      button("Remove part", () => run("parts.remove", { id: p.id }, status).catch(() => {})), status);
  } else insp.append(node("p", "Add a part from a template.", "small muted"));
  grid.append(list, centre, insp); root.append(grid);
  const a = section(root, "Assumptions"); assumptionList(a, "parts");
  if (p && !root.hidden) requestAnimationFrame(() => preview(p, info));
}
onSolved(s => { solution = s; const list = root.querySelector(".parts-list .panel-label"); if (list && s.outputs["parts.mass"]) list.textContent = `PARTS · ${format(s.outputs["parts.mass"].value, "kg", 4)} total`; });
project.store.on(({ ids }) => { if (ids.includes("parts")) render(); });
document.addEventListener("emergence:studio", e => { if (e.detail === "parts") render(); });
render();
