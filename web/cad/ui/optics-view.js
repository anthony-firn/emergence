// Optics studio: element prescription, meridional ray fan, spot diagrams.
import { project, onSolved } from "../app.js";
import { surfaces, sag, GLASSES, refractiveIndex, edgeThickness, elementMass, COOKE_SOURCE } from "../domains/optics.js";
import { node, svg, button, format, numberField, checkCards, certificates, assumptionList, section, linkedChip, live, run } from "./kit.js";
const root = document.getElementById("optics-studio");
let solution = null;
function layout(d, r) {
  const ss = surfaces(d), zmin = Math.min(...ss.map(s => s.z)) - 12, zmax = d.target.z + 6, ymax = Math.max(d.target.diameter / 2, ...d.elements.map(e => e.diameter / 2)) + 3;
  const s = svg("svg", { viewBox: `${zmin} ${-ymax} ${zmax - zmin} ${2 * ymax}`, class: "optics-layout", role: "img", "aria-label": "Meridional optical layout with traced rays", preserveAspectRatio: "xMidYMid meet" });
  s.append(svg("line", { x1: zmin, x2: zmax, y1: 0, y2: 0, class: "optic-axis" }));
  for (const e of d.elements) {
    const a = e.diameter / 2;
    if (e.kind === "stop") { s.append(svg("line", { x1: e.z, x2: e.z, y1: -ymax, y2: -a, class: "optic-stop" }), svg("line", { x1: e.z, x2: e.z, y1: a, y2: ymax, class: "optic-stop" })); continue; }
    const front = [], back = [];
    for (let i = 0; i <= 40; i++) { const y = -a + 2 * a * i / 40; front.push([e.z + (e.kind === "lens" ? sag(e.r1, Math.abs(y)) : 0), y]); back.push([e.z + e.thickness + (e.kind === "lens" ? sag(e.r2, Math.abs(y)) : 0), y]); }
    const path = [...front, ...back.reverse()].map((p, i) => `${i ? "L" : "M"}${p[0]},${-p[1]}`).join(" ") + "Z";
    const g = svg("path", { d: path, class: `optic-glass glass-${e.glass}` }); g.append(svg("title", {}, `${e.id} ${e.name} · ${e.glass}`)); s.append(g);
    s.append(svg("text", { x: e.z + (e.thickness || 0) / 2, y: -a - 1.2, "text-anchor": "middle", class: "optic-label" }, e.id));
  }
  s.append(svg("line", { x1: d.target.z, x2: d.target.z, y1: -d.target.diameter / 2, y2: d.target.diameter / 2, class: "optic-target" }));
  const colours = ["ray-0", "ray-1", "ray-2", "ray-3", "ray-4"];
  for (const [fi, f] of (r?.fields || []).entries()) {
    // Meridional fan: rays with x ≈ 0 from the traced beam.
    const fan = f.rays.filter((ray, i) => Math.abs(ray.path[0][0]) < d.source.diameter * .04).slice(0, 13);
    for (const ray of fan) s.append(svg("polyline", { points: ray.path.map(p => `${p[2]},${-p[1]}`).join(" "), class: `optic-ray ${colours[fi]}`, opacity: ray.lost ? .25 : .8 }));
  }
  return s;
}
function spots(d, r) {
  const wrap = node("div", undefined, "spot-row");
  for (const [fi, f] of (r?.fields || []).entries()) {
    const fig = node("figure", undefined, "spot"), hits = f.hits, [cx, cy] = f.centroid, R = Math.max(.002, 3 * (f.rms || .001));
    const s = svg("svg", { viewBox: `${-R} ${-R} ${2 * R} ${2 * R}`, role: "img", "aria-label": `Spot diagram at ${f.field} degrees` });
    s.append(svg("circle", { cx: 0, cy: 0, r: R * .98, class: "spot-frame" }), svg("circle", { cx: 0, cy: 0, r: f.rms, class: "spot-rms" }));
    for (const [x, y] of hits) s.append(svg("circle", { cx: x - cx, cy: -(y - cy), r: R / 90, class: `spot-dot ray-${fi}` }));
    fig.append(s, node("figcaption", `${f.field}° · ${(f.rms * 1000).toFixed(1)} µm RMS · y ${f.centroid[1].toFixed(3)} mm`)); wrap.append(fig);
  }
  return wrap;
}
function render() {
  const d = project.store.get("optics"), r = solution?.results.optics, x = r?.details, status = node("p", "", "studio-status"); status.id = "optics-status";
  root.replaceChildren();
  const head = node("div", undefined, "studio-heading"), title = node("div");
  title.append(node("p", `OPTICS / ${project.store.docs.project.name.toUpperCase()}`, "panel-label"), node("h2", "Every ray, traced for real."),
    node("p", `3D vector Snell refraction through spherical surfaces with per-ray Fresnel loss, chief rays aimed through the stop, paraxial focal lengths checked against the thick-lens equation. Reference design: ${COOKE_SOURCE}.`));
  const tools = node("div", undefined, "studio-file-tools");
  tools.append(button("⌖ Best focus", () => run("optics.focus", {}, status).catch(() => {}), "accent", "optics-focus"));
  head.append(title, tools); root.append(head);
  if (x) {
    const hero = node("div", undefined, "system-hero");
    for (const [label, v] of [["Focal length", format(x.paraxial.efl, "mm", 6)], ["Working f-number", `f/${x.fnumber.toFixed(2)}`], ["On-axis RMS spot", `${(x.trace.rms * 1000).toFixed(2)} µm`], ["Ground sample distance", `${x.gsd.toFixed(3)} cm/px at ${x.altitude} m`],
      ["Transmission (uncoated)", `${(100 * x.trace.transmitted).toFixed(1)}%`], ["Glass mass", format(x.mass, "g", 3)], ["Paraxial focus", format(x.paraxial.focus, "mm", 6)], ["Image plane", format(d.target.z, "mm", 6)]]) { const c = node("div", undefined, "hero-stat"); c.append(node("span", label), node("b", v)); hero.append(c); }
    root.append(hero);
  }
  const lay = section(root, "Meridional layout", `λ ${d.wavelength} nm · fields ${d.fields.join(", ")}° · ${d.rays} rays per field`);
  lay.append(layout(d, x), spots(d, x));
  // Prescription table.
  const pres = section(root, "Prescription");
  const table = node("table", undefined, "data-table prescription"), head2 = node("thead"), tr = node("tr");
  for (const h of ["Element", "Kind", "z vertex", "R1", "R2", "Thickness", "Clear Ø", "Glass", "n", "Edge", "Mass", ""]) tr.append(node("th", h)); head2.append(tr); table.append(head2);
  const body = node("tbody");
  for (const e of [...d.elements].sort((a, b) => a.z - b.z)) {
    const row = node("tr"); row.append(node("td", `${e.id} · ${e.name}`), node("td", e.kind));
    for (const [key, step] of [["z", .01], ["r1", .001], ["r2", .001], ["thickness", .001], ["diameter", .1]]) {
      const td = node("td"); if (e[key] === undefined) { td.textContent = "—"; row.append(td); continue; }
      const i = node("input"); i.type = "number"; i.step = step; i.value = e[key]; i.id = `optics-${e.id}-${key}`; i.setAttribute("aria-label", `${e.id} ${key}`);
      i.addEventListener("change", () => run("optics.update", { id: e.id, [key]: Number(i.value) }, status).catch(() => { i.value = e[key]; })); td.append(i); row.append(td);
    }
    const g = node("td"); if (e.glass) { const sel = node("select"); sel.setAttribute("aria-label", `${e.id} glass`); for (const [k, v] of Object.entries(GLASSES)) { const o = node("option", k); o.value = k; o.title = v.title; sel.append(o); } sel.value = e.glass; sel.addEventListener("change", () => run("optics.update", { id: e.id, glass: sel.value }, status).catch(() => {})); g.append(sel); } else g.textContent = "—";
    row.append(g, node("td", e.glass ? refractiveIndex(e.glass, d.wavelength).toFixed(5) : "—"), node("td", e.kind === "lens" ? `${edgeThickness(e).toFixed(2)} mm` : "—"), node("td", e.kind === "stop" ? "—" : `${elementMass(e).toFixed(2)} g`));
    const del = node("td"); del.append(button("Remove", () => run("optics.remove", { id: e.id }, status).catch(() => {}))); row.append(del); body.append(row);
  }
  table.append(body); pres.append(table);
  const add = node("div", undefined, "rtl-actions");
  add.append(button("＋ Lens", () => run("optics.add", { kind: "lens" }, status).catch(() => {})), button("＋ Window", () => run("optics.add", { kind: "window" }, status).catch(() => {})), button("＋ Stop", () => run("optics.add", { kind: "stop" }, status).catch(() => {})));
  pres.append(add);
  const sys = section(root, "Source, sensor and mission");
  const f = node("div", undefined, "field-grid");
  for (const [key, label, unit, min, max, step, value] of [["wavelength", "Wavelength", "nm", 300, 2500, 1, d.wavelength], ["diameter", "Entrance beam Ø", "mm", .1, 100, .1, d.source.diameter],
    ["target_z", "Image plane z", "mm", -1000, 5000, .001, d.target.z], ["target_diameter", "Sensor diagonal", "mm", .01, 500, .1, d.target.diameter], ["rays", "Rays per field", "", 16, 4000, 1, d.rays]])
    numberField(f, { label, id: `optics-${key}`, unit, value, min, max, step, op: "optics.set", args: v => ({ [key]: v }), status });
  sys.append(f);
  const alt = linkedChip(solution, "optics.altitude"); if (alt) sys.append(alt);
  sys.append(node("p", `Pixel pitch ${d.pixel} µm; ${d.requiredGsd > 0 ? `required ground sample distance ≤ ${d.requiredGsd} cm` : "no ground-sample requirement set"}. Glass indices from Sellmeier fits; N-SK16 n_d = ${refractiveIndex("N-SK16", 587.56).toFixed(5)}, F2 n_d = ${refractiveIndex("F2", 587.56).toFixed(5)} (catalogue 1.62041, 1.62004).`, "small muted"), status);
  if (r) { const c = section(root, "Optical checks"); checkCards(c, r.checks || []); certificates(c, r.certificates); }
  const a = section(root, "Assumptions"); assumptionList(a, "optics");
}
onSolved(s => { solution = s; });
live(["optics"], render);
render();
