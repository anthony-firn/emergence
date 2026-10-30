// Render studio: the open project from every domain's own geometry. Each
// domain builds its geometry in its own frame from the live documents and
// the kernels' meshes; the project's scene places it (at, rx, rz, and count
// copies spaced by phase about z). Domains the scene does not place are laid
// out side by side. With an Architecture domain, the site scope puts the
// design inside the building at true scale.
import { project, onSolved, mechanical } from "../app.js";
import { instances, MATERIALS } from "../domains/parts.js";
import { opticsMesh } from "../domains/optics.js";
import { reference } from "../core/project.js";
import { DOMAIN_INFO } from "../domains/index.js";
import { Renderer, M, boxMesh, discMesh } from "./gl.js";
import { modelMesh, sinkMesh } from "./meshes.js";
import { node, button, format, download, toast } from "./kit.js";
const root = document.getElementById("render-studio");
const COLOURS = { Mechanical: [.66, .72, .78, 1], cfrp: [.2, .21, .23, 1], "al-6061": [.75, .78, .82, 1], abs: [.24, .32, .4, 1], brass: [.78, .62, .3, 1], "steel-304": [.62, .64, .66, 1], copper: [.85, .55, .3, 1], "ti-64": [.6, .6, .63, 1], peek: [.7, .62, .45, 1], pa12: [.85, .85, .82, 1],
  sink: [.78, .8, .84, 1], board: [.12, .42, .26, 1], trace: [.85, .55, .3, 1], glass: [.55, .8, .95, .45], battery: [.85, .65, .2, 1], rotor: [.95, .55, .25, .2], wall: [.82, .8, .76, 1], window: [.55, .8, .95, .35], floor: [.42, .44, .46, 1],
  die: [.2, .2, .26, 1], logic: [.48, .4, .85, 1], memory: [.3, .6, .85, 1], io: [.85, .6, .3, 1], analog: [.4, .75, .5, 1] };
const boxes = new Map();
const box = (w, d, h) => { const k = `${w}|${d}|${h}`; if (!boxes.has(k)) boxes.set(k, boxMesh(w, d, h)); return boxes.get(k); };
const deg = Math.PI / 180;
let viewer = null, scope = "design", section = false, solution = null, generation = 0;
const canvas = node("canvas", undefined, "render-canvas"); canvas.id = "render-canvas"; canvas.tabIndex = 0; canvas.setAttribute("aria-label", "Whole-system 3D view. Drag to orbit, scroll to zoom, shift-drag to pan, click to identify a part.");
const status = node("p", "", "mesh-info"); status.id = "render-status";
const legend = node("div", undefined, "render-legend");
// Each domain's geometry in its own frame, millimetres.
const GEOMETRY = {
  async parts(docs) {
    const out = [];
    for (const p of docs.parts.parts) { const mesh = await modelMesh(p.model, p.overrides);
      for (const inst of instances(p)) out.push({ mesh, matrix: M.mul(M.rz(inst.rz * deg), M.translate(...inst.by)), color: COLOURS[p.material] || COLOURS.Mechanical, label: p.name, detail: `${MATERIALS[p.material].title} · ${format(mesh.volume * MATERIALS[p.material].density / 1000, "g", 3)} each` }); }
    return out;
  },
  async mech() {
    // The live mechanical design with its linked inputs applied.
    const params = { ...mechanical.state().params, ...Object.fromEntries(Object.entries(solution?.applied || {}).filter(([k]) => k.startsWith("mech.")).map(([k, v]) => [k.slice(5), v.value])) };
    if (params.fins !== undefined) params.fins = Math.round(params.fins);
    const sink = await sinkMesh(params).catch(() => null);
    return sink ? [{ mesh: sink, matrix: M.identity(), color: COLOURS.sink, label: "Heat sink", detail: `${params.fins} fins · Q1 ${format(solution?.outputs["mech.q1_tj"]?.value, "degC", 3)}` }] : [];
  },
  pcb(docs) {
    const b = docs.pcb, t = b.stackup.thickness, out = [{ mesh: box(b.width, b.height, t), matrix: M.identity(), color: COLOURS.board, label: "Board", detail: `${b.width} × ${b.height} mm · ${b.routes.length} traces` }];
    for (const c of b.components) out.push({ mesh: box(5, 4, 3), matrix: M.translate(c.x - b.width / 2, b.height / 2 - c.y, t), color: c.footprint === "dip8" ? [.15, .15, .17, 1] : c.footprint === "transistor" ? [.3, .3, .32, 1] : [.6, .55, .45, 1], label: c.id, detail: c.name });
    for (const r of b.routes) for (let k = 1; k < r.points.length; k++) { const [x0, y0] = r.points[k - 1], [x1, y1] = r.points[k], len = Math.hypot(x1 - x0, y1 - y0); if (!len) continue;
      out.push({ mesh: box(len, r.width, .1), matrix: M.mul(M.translate((x0 + x1) / 2 - b.width / 2, b.height / 2 - (y0 + y1) / 2, r.layer === "top" ? t : -.1), M.rz(-Math.atan2(y1 - y0, x1 - x0))), color: COLOURS.trace, label: `${r.from} → ${r.to}`, detail: `${r.layer} · ${r.width} mm` }); }
    return out;
  },
  silicon(docs) {
    // The die at true size: the floorplan is in micrometres.
    const d = docs.silicon, k = 1e-3, out = [{ mesh: box(d.width * k, d.height * k, .3), matrix: M.identity(), color: COLOURS.die, label: "Die", detail: `${(d.width * k).toFixed(1)} × ${(d.height * k).toFixed(1)} mm` }];
    for (const b of d.blocks) out.push({ mesh: box(b.width * k, b.height * k, .05), matrix: M.translate((b.x + b.width / 2 - d.width / 2) * k, (d.height / 2 - b.y - b.height / 2) * k, .3), color: COLOURS[b.kind], label: b.name, detail: `${b.kind} · ${b.power} mW planned` });
    return out;
  },
  async optics(docs) {
    const out = [];
    for (const e of docs.optics.elements) { const m = opticsMesh(e); if (!m) continue; out.push({ mesh: await modelMesh(m, {}, 1.5e5), matrix: M.identity(), color: COLOURS.glass, label: `${e.id} ${e.name}`, detail: `${e.glass} · R ${e.r1 ?? "∞"} / ${e.r2 ?? "∞"}` }); }
    return out;
  },
  flight(docs) {
    // The battery under the body and the momentum-theory actuator discs, X layout.
    const f = docs.flight, s = Math.cbrt(f.battery.mass / 3.9), disc = discMesh(f.prop / 2), load = solution?.outputs["flight.disk_loading"]?.value, arm = Math.max(f.prop * .97, 120);
    const out = [{ mesh: box(180 * s, 75 * s, 55 * s), matrix: M.translate(0, -10, -2.5 - 55 * s), color: COLOURS.battery, label: "Battery", detail: `${f.battery.cells}S ${f.battery.capacity} Ah · ${format(f.battery.mass, "kg", 3)}` }];
    for (let i = 0; i < f.rotors; i++) out.push({ mesh: disc, matrix: M.mul(M.rz((45 + 360 * i / f.rotors) * deg), M.translate(arm, 0, 48)), color: COLOURS.rotor, label: "Rotor actuator disc", detail: `${(f.prop / 25.4).toFixed(1)}″ · disk loading ${format(load, "Pa", 3)}` });
    return out;
  },
};
function building(docs, out) {
  const a = docs.architecture, mm = 1000, wall = 200;
  for (const r of a.rooms) {
    const x0 = r.x * mm, y0 = -r.y * mm, w = r.width * mm, d = r.depth * mm, h = r.height * mm, z = solution?.results.architecture?.details?.rooms.find(q => q.id === r.id);
    out.push({ mesh: box(w, d, 60), matrix: M.translate(x0 + w / 2, y0 - d / 2, -60), color: COLOURS.floor, label: r.name, detail: z ? `${z.temperature.toFixed(1)} °C · ${r.mode}` : r.mode });
    for (const [side, cx, cy, len, rot] of [["north", x0 + w / 2, y0, w, 0], ["south", x0 + w / 2, y0 - d, w, 0], ["west", x0, y0 - d / 2, d, Math.PI / 2], ["east", x0 + w, y0 - d / 2, d, Math.PI / 2]]) {
      const wl = r.walls[side]; if (!wl.exterior && side !== "south" && side !== "east") continue;
      const glass = wl.window > 0 ? Math.min(.9, wl.window / (len / mm * r.height)) : 0;
      out.push({ mesh: box(len * (1 - glass), wall, h), matrix: M.mul(M.translate(cx, cy, 0), M.rz(rot), M.translate(len * glass / 2, 0, 0)), color: COLOURS.wall, label: `${r.name} · ${side} wall`, detail: wl.exterior ? `U ${wl.U} W/m²K` : "interior" });
      if (glass) out.push({ mesh: box(len * glass, wall / 3, h * .6), matrix: M.mul(M.translate(cx, cy, h * .25), M.rz(rot), M.translate(-len * (1 - glass) / 2, 0, 0)), color: COLOURS.window, label: `${r.name} · ${side} window`, detail: `${wl.window} m² · SHGC ${a.glazing.shgc}` });
    }
  }
}
// Axis-aligned bounds of meshes under their matrices, from each mesh's own box.
function bounds(items) {
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (const it of items) {
    const P = it.mesh.positions; if (!P?.length) continue;
    if (!it.mesh.box) { const a = [Infinity, Infinity, Infinity], b = [-Infinity, -Infinity, -Infinity]; for (let i = 0; i < P.length; i += 3) for (let k = 0; k < 3; k++) { a[k] = Math.min(a[k], P[i + k]); b[k] = Math.max(b[k], P[i + k]); } it.mesh.box = [a, b]; }
    const [a, b] = it.mesh.box, m = it.matrix;
    for (const x of [a[0], b[0]]) for (const y of [a[1], b[1]]) for (const z of [a[2], b[2]]) { const p = [0, 1, 2].map(k => m[k] * x + m[4 + k] * y + m[8 + k] * z + m[12 + k]); for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], p[k]); hi[k] = Math.max(hi[k], p[k]); } }
  }
  return [lo, hi];
}
const geometryDomains = () => project.active().filter(d => GEOMETRY[d]);
async function build() {
  const gen = ++generation; viewer ||= new Renderer(canvas);
  viewer.onpick = it => { legend.textContent = it ? `${it.label}${it.detail ? " · " + it.detail : ""}` : "Click a part to identify it."; };
  status.textContent = "Meshing every domain in the Rust kernel…";
  const docs = { ...project.store.docs, mech: mechanical.state() }, scene = docs.project.scene, items = [], t0 = performance.now();
  try {
    // Placed domains first; the rest side by side beyond them.
    let shelf = null;
    const groups = await Promise.all(geometryDomains().map(async d => ({ d, items: await GEOMETRY[d](docs) })));
    for (const { d, items: local } of groups.sort((a, b) => !!scene.placements[b.d] - !!scene.placements[a.d])) {
      if (!local.length) continue;
      const p = scene.placements[d];
      if (!p) {
        const [lo, hi] = bounds(local), placed = bounds(items);
        shelf ??= Number.isFinite(placed[1][0]) ? placed[1][0] + 40 : 0;
        const by = M.translate(shelf - lo[0], -(lo[1] + hi[1]) / 2, -lo[2]); shelf += hi[0] - lo[0] + 40;
        for (const it of local) items.push({ ...it, matrix: M.mul(by, it.matrix), domain: d });
        continue;
      }
      const count = Math.max(1, Math.round(reference(docs, p.count ?? 1) ?? 1));
      for (let i = 0; i < count; i++) {
        const T = M.mul(M.rz(((p.phase ?? 0) + 360 * i / count) * deg), M.translate(...p.at), M.rz((p.rz ?? 0) * deg), M.rx((p.rx ?? 0) * deg));
        for (const it of local) items.push({ ...it, matrix: M.mul(T, it.matrix), domain: d });
      }
    }
    const site = scope === "site" && project.active().includes("architecture");
    if (site) {
      // The design in its room at true scale, then the building around it.
      const a = docs.architecture, room = a.rooms.find(r => r.id === scene.site?.room) || a.rooms.find(r => r.equipment) || a.rooms[0], mm = 1000;
      const by = M.translate((room.x + room.width / 2) * mm, -(room.y + room.depth / 2) * mm, scene.site?.z ?? 900);
      for (const it of items) it.matrix = M.mul(by, it.matrix);
      building(docs, items);
    }
    if (gen !== generation) return;
    if (!items.length) { viewer.set([], { size: 100, step: 10, z: 0 }); status.textContent = project.active().length ? "None of this project's domains has 3D geometry yet: parts, mechanical, PCB, silicon, optics, flight and architecture do." : "Add a domain to see it here."; return; }
    for (const it of items) it.gpu = it.mesh.gpu ||= viewer.upload(it.mesh);
    const [lo, hi] = bounds(items), span = Math.max(hi[0] - lo[0], hi[1] - lo[1], 1);
    viewer.plane = section ? [0, 1, 0, -(lo[1] + hi[1]) / 2] : null;
    viewer.set(items, { size: span * .75, step: 10 ** Math.floor(Math.log10(span / 4)), z: lo[2] - span * .01 });
    const key = `${scope}|${project.active().join()}`;
    if (viewer.fitted !== key) { viewer.fit(lo, hi); viewer.fitted = key; }
    const triangles = items.reduce((s, it) => s + (it.mesh.indices?.length ?? 0) / 3, 0);
    status.textContent = `${items.length} objects · ${Math.round(triangles).toLocaleString()} triangles · from ${[...new Set(items.map(i => i.domain).filter(Boolean))].map(d => DOMAIN_INFO[d].title).join(", ")}${site ? " and the building" : ""} · ${(performance.now() - t0).toFixed(0)} ms`;
  } catch (e) { status.textContent = `Scene unavailable: ${e.message}`; }
}
function render() {
  root.replaceChildren();
  const has = project.active().includes("architecture"), only = has && !geometryDomains().some(d => d !== "architecture");
  if (only) scope = "site"; else if (!has) scope = "design";
  const head = node("div", undefined, "studio-heading"), title = node("div");
  title.append(node("p", `RENDER / ${project.store.docs.project.name.toUpperCase()}`, "panel-label"), node("h2", "The whole design, from each domain's own geometry."),
    node("p", "Parts, heat sinks and lenses are meshed by the Rust kernel from the live documents; boards, dies, batteries, rotor discs and rooms come straight from their parameters. Where each domain sits comes from the project's scene."));
  const tools = node("div", undefined, "studio-file-tools");
  if (has && !only) tools.append(button(scope === "design" ? "⌂ Site scope" : "◭ Design scope", () => { scope = scope === "design" ? "site" : "design"; render(); }, "", "render-scope"));
  const xray = button("◌ X-ray", () => { viewer.xray = !viewer.xray; xray.setAttribute("aria-pressed", String(viewer.xray)); viewer.draw(); }, "", "render-xray"); xray.setAttribute("aria-pressed", String(!!viewer?.xray));
  const cut = button("◧ Section", () => { section = !section; cut.setAttribute("aria-pressed", String(section)); build(); }, "", "render-section"); cut.setAttribute("aria-pressed", String(section));
  tools.append(xray, cut, button("⌖ Fit", () => { viewer.fitted = null; build(); }, "", "render-fit"), button("↓ PNG", async () => { const blob = await viewer.png(); download(blob, `emergence-${scope}.png`, "image/png"); toast("Rendered view saved."); }, "accent", "render-png"));
  head.append(title, tools); root.append(head);
  const stage = node("div", undefined, "render-stage"); stage.append(canvas, legend); root.append(stage, status);
  legend.textContent = "Click a part to identify it.";
  build();
}
onSolved(s => { solution = s; if (!root.hidden) build(); });
project.store.on(({ ids }) => { if (root.hidden) return; if (ids.includes("project")) render(); else build(); });
document.addEventListener("emergence:studio", e => { if (e.detail === "render") try { render(); } catch (err) { root.replaceChildren(node("p", `3D view unavailable: ${err.message}`, "conflict")); } });
// Built on first view: hidden studios do not queue kernel work at start-up.
if (!root.hidden) try { render(); } catch (e) { root.append(node("p", `3D view unavailable: ${e.message}`, "conflict")); }
