// Board-level additions to the PCB editor: the kernel board (SI, y up) for
// cad::kits::pcb, current density, an obstacle-aware 45° router, KiCad export
// and a Gerber reader that checks the fabrication files against the design.
import { FOOTPRINTS, pads, bodyBounds, pcbChecks, segmentDistance } from "../studio-model.js";
export const OZ = 34.8e-3; // mm of copper per oz/ft² (IPC-4562A), as in the kit
export const DRILL = .9;   // illustrative finished drill for every through-hole pad, mm

export function kernelBoard(d) {
  const m = v => v * 1e-3, flip = y => m(d.height - y), copper = { kind: "copper", ounces: d.stackup.copper };
  const features = [], currents = [];
  for (const r of d.routes) {
    features.push({ kind: "trace", layer: r.layer === "top" ? 1 : 0, points: r.points.map(([x, y]) => [m(x), flip(y)]), width: m(r.width) });
    currents.push(d.currents[netOf(d, r)] ?? 0);
  }
  const vias = pads(d).map(p => ({ center: [m(p.x), flip(p.y)], drill: m(DRILL), annular_ring: m(p.radius - DRILL / 2) }));
  return {
    board: { outline: { op: "rect", center: [m(d.width / 2), m(d.height / 2)], half: [m(d.width / 2), m(d.height / 2)] },
      stackup: [copper, { kind: "dielectric", thickness: m(Math.max(.05, d.stackup.thickness - 2 * d.stackup.copper * OZ)), eps_r: d.stackup.epsr }, copper],
      features, vias },
    currents,
  };
}
export const netOf = (d, r) => pads(d).find(p => p.id === r.from)?.net;
export function densityIssues(d) {
  const t = d.stackup.copper * OZ;
  return d.routes.flatMap(r => {
    const a = d.currents[netOf(d, r)] ?? 0, j = a / (r.width * t);
    return j > d.rules.density ? [{ code: "density", object: r.id, message: `${r.id} carries ${a} A at ${j.toFixed(0)} A/mm², above ${d.rules.density} A/mm². Widen it or use a pour.` }] : [];
  });
}
// A* on the snap grid with obstacles inflated by clearance and half the new
// track width (plus half a cell, so diagonal steps cannot clip a corner).
export function autoroute(d, from, to, layer, width) {
  const all = pads(d), a = all.find(p => p.id === from), b = all.find(p => p.id === to);
  if (!a || !b || a.id === b.id || !a.net || a.net !== b.net) throw new Error("Choose two different pads on the same assigned net.");
  const g = Math.max(.1, d.grid), nx = Math.floor(d.width / g) + 1, ny = Math.floor(d.height / g) + 1;
  if (nx * ny > 400000) throw new Error("This board and grid exceed the router's 400,000-cell budget. Use a coarser snap grid.");
  const blocked = new Uint8Array(nx * ny), margin = d.rules.clearance + width / 2 + g / 2;
  const mark = (x0, y0, x1, y1, test) => {
    for (let j = Math.max(0, Math.floor(y0 / g)); j <= Math.min(ny - 1, Math.ceil(y1 / g)); j++)
      for (let i = Math.max(0, Math.floor(x0 / g)); i <= Math.min(nx - 1, Math.ceil(x1 / g)); i++) if (test(i * g, j * g)) blocked[j * nx + i] = 1;
  };
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const x = i * g, y = j * g; if (Math.min(x, y, d.width - x, d.height - y) - width / 2 < d.rules.edge) blocked[j * nx + i] = 1; }
  for (const p of all) if (p.net !== a.net) { const r = p.radius + margin; mark(p.x - r, p.y - r, p.x + r, p.y + r, (x, y) => Math.hypot(x - p.x, y - p.y) < r); }
  for (const r of d.routes) {
    if (r.layer !== layer || netOf(d, r) === a.net) continue;
    for (let k = 1; k < r.points.length; k++) {
      const [u, v] = [r.points[k - 1], r.points[k]], rr = r.width / 2 + margin;
      mark(Math.min(u[0], v[0]) - rr, Math.min(u[1], v[1]) - rr, Math.max(u[0], v[0]) + rr, Math.max(u[1], v[1]) + rr, (x, y) => segmentDistance([x, y], [x, y], u, v) < rr);
    }
  }
  const node = p => [Math.round(p.x / g), Math.round(p.y / g)], [si, sj] = node(a), [ti, tj] = node(b), start = sj * nx + si, goal = tj * nx + ti;
  blocked[start] = 0; blocked[goal] = 0;
  const dirs = [[1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1], [1, -1]], S = nx * ny * 8;
  const cost = new Float64Array(S).fill(Infinity), prev = new Int32Array(S).fill(-1), heap = [];
  const h = k => { const i = k % nx, j = Math.floor(k / nx), dx = Math.abs(i - ti), dy = Math.abs(j - tj); return Math.max(dx, dy) + (Math.SQRT2 - 1) * Math.min(dx, dy); };
  const push = (s, f) => { heap.push([f, s]); let i = heap.length - 1; while (i) { const p = (i - 1) >> 1; if (heap[p][0] <= heap[i][0]) break; [heap[p], heap[i]] = [heap[i], heap[p]]; i = p; } };
  const pop = () => { const top = heap[0], last = heap.pop(); if (heap.length) { heap[0] = last; let i = 0; for (;;) { const l = 2 * i + 1, r = l + 1; let m = i; if (l < heap.length && heap[l][0] < heap[m][0]) m = l; if (r < heap.length && heap[r][0] < heap[m][0]) m = r; if (m === i) break; [heap[m], heap[i]] = [heap[i], heap[m]]; i = m; } } return top; };
  for (let k = 0; k < 8; k++) { cost[start * 8 + k] = 0; push(start * 8 + k, h(start)); }
  let found = -1, expanded = 0;
  while (heap.length) {
    const [, s] = pop(), n = Math.floor(s / 8), dir = s % 8;
    if (n === goal) { found = s; break; }
    if (++expanded > 3e6) break;
    const i = n % nx, j = Math.floor(n / nx);
    for (let k = 0; k < 8; k++) {
      const ii = i + dirs[k][0], jj = j + dirs[k][1];
      if (ii < 0 || jj < 0 || ii >= nx || jj >= ny || blocked[jj * nx + ii]) continue;
      const turn = k === dir ? 0 : Math.min((k - dir + 8) % 8, (dir - k + 8) % 8) > 1 ? .9 : .3, c = cost[s] + (k % 2 ? Math.SQRT2 : 1) + turn, t = (jj * nx + ii) * 8 + k;
      if (c < cost[t]) { cost[t] = c; prev[t] = s; push(t, c + h(jj * nx + ii)); }
    }
  }
  if (found < 0) throw new Error(`No clear ${layer} path from ${from} to ${to} at ${d.rules.clearance} mm clearance. Try the other layer, a narrower track or move a component.`);
  const cells = [];
  for (let s = found; s >= 0; s = prev[s]) { const n = Math.floor(s / 8); if (cells.at(-1) !== n) cells.push(n); }
  cells.reverse();
  const pts = [[a.x, a.y], ...cells.map(n => [(n % nx) * g, Math.floor(n / nx) * g]), [b.x, b.y]].map(p => p.map(v => Math.round(v * 1e6) / 1e6));
  const simple = pts.filter((p, i) => {
    if (i === 0 || i === pts.length - 1) return true;
    const [u, v] = [pts[i - 1], pts[i + 1]];
    return Math.abs((p[0] - u[0]) * (v[1] - u[1]) - (p[1] - u[1]) * (v[0] - u[0])) > 1e-9;
  }).filter((p, i, all) => i === 0 || p[0] !== all[i - 1][0] || p[1] !== all[i - 1][1]);
  if (simple.length > 32) throw new Error("That path needs more than 32 corners. Move components to open a straighter channel.");
  return { id: "route-" + crypto.randomUUID().slice(0, 12), from, to, layer, width, points: simple };
}
// KiCad 8 board file. Pads are written at their rotated offsets on a
// zero-angle footprint so rotation conventions cannot disagree. Written to
// the documented S-expression format; not yet checked by kicad-cli.
export function kicadPCB(d) {
  const q = s => '"' + String(s).replace(/\\/g, "\\\\").replace(/"/g, '\\"') + '"', n = v => Number(v.toFixed(6));
  let seq = 0; const uuid = () => { const h = (++seq).toString(16).padStart(12, "0"); return `00000000-0000-4000-8000-${h}`; };
  const nets = ["", ...[...new Set(pads(d).map(p => p.net).filter(Boolean))].sort()], netId = name => nets.indexOf(name ?? "");
  const out = ["(kicad_pcb (version 20240108) (generator \"emergence\") (generator_version \"1.0\")",
    `  (general (thickness ${n(d.stackup.thickness)}) (legacy_teardrops no))`, "  (paper \"A4\")",
    "  (layers (0 \"F.Cu\" signal) (31 \"B.Cu\" signal) (38 \"B.Mask\" user) (39 \"F.Mask\" user) (36 \"B.SilkS\" user \"B.Silkscreen\") (37 \"F.SilkS\" user \"F.Silkscreen\") (44 \"Edge.Cuts\" user) (47 \"F.CrtYd\" user \"F.Courtyard\") (49 \"F.Fab\" user))",
    "  (setup (pad_to_mask_clearance 0))", ...nets.map((name, i) => `  (net ${i} ${q(name)})`)];
  const byComponent = Map.groupBy ? Map.groupBy(pads(d), p => p.component) : pads(d).reduce((m, p) => m.set(p.component, [...(m.get(p.component) || []), p]), new Map());
  for (const c of d.components) {
    const b = bodyBounds(c), font = "(effects (font (size 1 1) (thickness 0.15)))";
    out.push(`  (footprint ${q("Emergence:" + c.footprint)} (layer "F.Cu") (uuid ${q(uuid())}) (at ${n(c.x)} ${n(c.y)})`,
      `    (property "Reference" ${q(c.id)} (at 0 ${n(-b.height / 2 - 1.2)} 0) (layer "F.SilkS") (uuid ${q(uuid())}) ${font})`,
      `    (property "Value" ${q(c.name)} (at 0 ${n(b.height / 2 + 1.2)} 0) (layer "F.Fab") (uuid ${q(uuid())}) ${font})`,
      `    (fp_rect (start ${n(b.x - c.x)} ${n(b.y - c.y)}) (end ${n(b.x + b.width - c.x)} ${n(b.y + b.height - c.y)}) (stroke (width 0.05) (type solid)) (fill none) (layer "F.CrtYd") (uuid ${q(uuid())}))`);
    for (const p of byComponent.get(c.id) || [])
      out.push(`    (pad ${q(p.number)} thru_hole circle (at ${n(p.x - c.x)} ${n(p.y - c.y)}) (size ${n(2 * p.radius)} ${n(2 * p.radius)}) (drill ${DRILL}) (layers "*.Cu" "*.Mask") (net ${netId(p.net)} ${q(p.net || "")}) (uuid ${q(uuid())}))`);
    out.push("  )");
  }
  for (const r of d.routes) for (let k = 1; k < r.points.length; k++)
    out.push(`  (segment (start ${n(r.points[k - 1][0])} ${n(r.points[k - 1][1])}) (end ${n(r.points[k][0])} ${n(r.points[k][1])}) (width ${n(r.width)}) (layer ${r.layer === "top" ? '"F.Cu"' : '"B.Cu"'}) (net ${netId(netOf(d, r))}) (uuid ${q(uuid())}))`);
  out.push(`  (gr_rect (start 0 0) (end ${n(d.width)} ${n(d.height)}) (stroke (width 0.1) (type default)) (fill none) (layer "Edge.Cuts") (uuid ${q(uuid())}))`, ")");
  return out.join("\n") + "\n";
}
// Reads back the region-only RS-274X the kernel writes: format 4.6, mm,
// G36/G37 regions of D02/D01 moves. Returns loops in millimetres, y up.
export function parseGerber(text) {
  const loops = []; let x = 0, y = 0, current = null, scale = 1e-6;
  for (const raw of text.split(/\*|\n/)) {
    const s = raw.trim(); if (!s) continue;
    if (/%FSLAX(\d)(\d)Y/.test(s)) { const m = s.match(/%FSLAX(\d)(\d)Y/); scale = 10 ** -Number(m[2]); continue; }
    if (s.startsWith("G36")) { current = []; continue; }
    if (s.startsWith("G37")) { if (current?.length > 2) loops.push(current); current = null; continue; }
    const m = s.match(/^(?:X(-?\d+))?(?:Y(-?\d+))?D0([12])$/);
    if (!m) continue;
    if (m[1] !== undefined) x = Number(m[1]) * scale; if (m[2] !== undefined) y = Number(m[2]) * scale;
    if (current) { if (m[3] === "2" && current.length > 2) { loops.push(current); current = []; } current.push([x, y]); }
  }
  return loops;
}
// Even-odd fill of the read-back loops against a raster of the drawn copper.
export function gerberCheck(d, gerber, layer, cell = .1) {
  const loops = parseGerber(gerber), inside = (px, py) => {
    let k = false; for (const l of loops) for (let i = 0, j = l.length - 1; i < l.length; j = i++) {
      const [xi, yi] = l[i], [xj, yj] = l[j]; if ((yi > py) !== (yj > py) && px < (xj - xi) * (py - yi) / (yj - yi) + xi) k = !k;
    } return k;
  };
  const ps = pads(d), rs = d.routes.filter(r => r.layer === layer);
  const drawn = (x, y) => ps.some(p => { const r = Math.hypot(x - p.x, y - p.y); return r <= p.radius && r >= DRILL / 2; })
    || rs.some(r => r.points.slice(1).some((v, k) => segmentDistance([x, y], [x, y], r.points[k], v) <= r.width / 2) && !ps.some(p => Math.hypot(x - p.x, y - p.y) < DRILL / 2));
  let both = 0, onlyG = 0, onlyD = 0;
  for (let y = cell / 2; y < d.height; y += cell) for (let x = cell / 2; x < d.width; x += cell) {
    const g = inside(x, d.height - y), w = drawn(x, y);
    if (g && w) both++; else if (g) onlyG++; else if (w) onlyD++;
  }
  const area = c => c * cell * cell;
  return { loops: loops.length, drawnArea: area(both + onlyD), gerberArea: area(both + onlyG), mismatch: area(onlyG + onlyD), agreement: both / Math.max(1, both + onlyG + onlyD) };
}
export const pcbAssumptions = d => [
  { id: "pcb.epsr", label: `Dielectric εr ${d.stackup.epsr} (FR-4 nominal)`, kind: "nominal", source: "Real laminate runs 4.2–4.8 with resin content and weave; controlled impedance needs the fabricator's value." },
  { id: "pcb.footprints", label: "Footprints and pad sizes", kind: "illustrative", source: "Estimates, not IPC-7351 land patterns. Replace with the manufacturer's package drawing." },
  { id: "pcb.currents", label: "Every trace on a net carries that net's full RMS current", kind: "nominal", source: "Conservative for copper loss; branch currents are not solved." },
  { id: "pcb.drill", label: `${DRILL} mm finished drill on every pad`, kind: "illustrative", source: "Sets via barrels in the drill file and kernel board." },
];
export { FOOTPRINTS, pcbChecks };
