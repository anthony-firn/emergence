// Board-level additions to the PCB editor: the kernel board (SI, y up) for
// cad::kits::pcb, current density, an obstacle-aware 45° router, KiCad export
// and a Gerber reader that checks the fabrication files against the design.
import { FOOTPRINTS, NO_CONNECT, pads, bodyBounds, pcbChecks, segmentDistance, inPolygon, padDistance } from "../studio-model.js";
export const OZ = 34.8e-3; // mm of copper per oz/ft² (IPC-4562A), as in the kit
export const DRILL = .9;   // illustrative finished drill for every through-hole pad, mm

// A net with a pour carries its current there; its traces are stubs and branches
// with small local currents, taken as zero.
const poured = (d, net) => (d.pours ?? []).some(p => p.net === net);
export const traceCurrent = (d, r) => poured(d, netOf(d, r)) ? 0 : d.currents[netOf(d, r)] ?? 0;
export function kernelBoard(d) {
  const m = v => v * 1e-3, flip = y => m(d.height - y), copper = { kind: "copper", ounces: d.stackup.copper };
  const features = [], currents = [];
  for (const r of d.routes) {
    features.push({ kind: "trace", layer: r.layer === "top" ? 1 : 0, points: r.points.map(([x, y]) => [m(x), flip(y)]), width: m(r.width) });
    currents.push(traceCurrent(d, r));
  }
  // Surface pads are copper on the top layer; only through pads are drilled.
  const all = pads(d), vias = all.filter(p => !p.smd).map(p => { const drill = p.drill ?? DRILL; return { center: [m(p.x), flip(p.y)], drill: m(drill), annular_ring: m(p.radius - drill / 2) }; });
  for (const p of all) if (p.smd) features.push({ kind: "pad", layer: 1, center: [m(p.x), flip(p.y)], size: [m(p.w), m(p.h)] });
  // A pour is its outline less a clearance around every other net's copper.
  for (const pour of d.pours ?? []) {
    const c = d.rules.clearance, keep = [];
    for (const p of all) if (p.net !== pour.net && (!p.smd || pour.layer === "top")) keep.push(p.w ? { op: "rect", center: [m(p.x), flip(p.y)], half: [m(p.w / 2 + c), m(p.h / 2 + c)] } : { op: "circle", center: [m(p.x), flip(p.y)], radius: m(p.radius + c) });
    for (const r of d.routes) if (r.layer === pour.layer && netOf(d, r) !== pour.net) for (let k = 1; k < r.points.length; k++) keep.push({ op: "slot", a: [m(r.points[k - 1][0]), flip(r.points[k - 1][1])], b: [m(r.points[k][0]), flip(r.points[k][1])], radius: m(r.width / 2 + c) });
    const outline = { op: "polygon", loops: [pour.points.map(([x, y]) => [m(x), flip(y)])] };
    features.push({ kind: "pour", layer: pour.layer === "top" ? 1 : 0, region: keep.length ? { op: "difference", a: outline, b: { op: "union", parts: keep } } : outline });
  }
  return {
    board: { outline: { op: "rect", center: [m(d.width / 2), m(d.height / 2)], half: [m(d.width / 2), m(d.height / 2)] },
      stackup: [copper, { kind: "dielectric", thickness: m(Math.max(.05, d.stackup.thickness - 2 * d.stackup.copper * OZ)), eps_r: d.stackup.epsr }, copper],
      features, vias },
    currents,
  };
}
export const netOf = (d, r) => { const [id, n] = r.from.split(":"), c = d.components.find(c => c.id === id); return c?.nets[Number(n) - 1]; };
// Each component's dissipation is a port of its own.
export const powerPort = id => `${id.toLowerCase().replace(/[^a-z0-9_]/g, "_")}_power`;
export function densityIssues(d) {
  const t = d.stackup.copper * OZ;
  const traces = d.routes.flatMap(r => {
    const a = traceCurrent(d, r), j = a / (r.width * t);
    return j > d.rules.density ? [{ code: "density", object: r.id, message: `${r.id} carries ${a} A at ${j.toFixed(0)} A/mm², above ${d.rules.density} A/mm². Widen it or use a pour.` }] : [];
  });
  const fills = (d.pours ?? []).flatMap(p => { const a = d.currents[p.net] ?? 0, { width } = pourSize(p), j = a / (width * t); return j > d.rules.density ? [{ code: "density", object: p.id, message: `Pour ${p.id} carries ${a} A at ${j.toFixed(0)} A/mm², above ${d.rules.density} A/mm².` }] : []; });
  return [...traces, ...fills];
}
// A pour's current path, as a strip: its long side, and its area over that length.
export function pourSize(p) {
  const xs = p.points.map(q => q[0]), ys = p.points.map(q => q[1]), length = Math.max(Math.max(...xs) - Math.min(...xs), Math.max(...ys) - Math.min(...ys));
  const area = Math.abs(p.points.reduce((s, [x, y], i) => { const [u, v] = p.points[(i + 1) % p.points.length]; return s + x * v - u * y; }, 0)) / 2;
  return { length, area, width: area / Math.max(1e-6, length) };
}
// Copper loss in the pours, I²R of that strip at 40 °C. An estimate: real current
// spreads and crowds; only a field solve of the drawn pour would say exactly.
export function pourLoss(d) { const t = d.stackup.copper * OZ, rho = 1.678e-8 * (1 + 3.93e-3 * 20); return (d.pours ?? []).reduce((s, p) => { const a = d.currents[p.net] ?? 0, { length, width } = pourSize(p); return s + a * a * rho * (length * 1e-3) / (width * t * 1e-6); }, 0); }
// A* on the snap grid with obstacles inflated by clearance and half the new
// track width (plus half a cell, so diagonal steps cannot clip a corner).
export function autoroute(d, from, to, layer, width) {
  const all = pads(d), a = all.find(p => p.id === from), b = all.find(p => p.id === to);
  if (!a || !b || a.id === b.id || !a.net || a.net === NO_CONNECT || a.net !== b.net) throw new Error("Choose two different pads on the same assigned net.");
  if ((a.smd || b.smd) && layer !== "top") throw new Error("A surface pad is reached on the top layer; use a via for the bottom.");
  const g = Math.max(.1, d.grid), nx = Math.floor(d.width / g) + 1, ny = Math.floor(d.height / g) + 1;
  if (nx * ny > 400000) throw new Error("This board and grid exceed the router's 400,000-cell budget. Use a coarser snap grid.");
  const blocked = new Uint8Array(nx * ny), margin = d.rules.clearance + width / 2 + g / 2;
  const mark = (x0, y0, x1, y1, test) => {
    for (let j = Math.max(0, Math.floor(y0 / g)); j <= Math.min(ny - 1, Math.ceil(y1 / g)); j++)
      for (let i = Math.max(0, Math.floor(x0 / g)); i <= Math.min(nx - 1, Math.ceil(x1 / g)); i++) if (test(i * g, j * g)) blocked[j * nx + i] = 1;
  };
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const x = i * g, y = j * g; if (Math.min(x, y, d.width - x, d.height - y) - width / 2 < d.rules.edge) blocked[j * nx + i] = 1; }
  for (const p of all) if ((p.net !== a.net || p.net === NO_CONNECT) && (!p.smd || layer === "top")) { const r = p.radius + margin; mark(p.x - r, p.y - r, p.x + r, p.y + r, (x, y) => padDistance(p, x, y) < margin); }
  for (const r of d.routes) {
    if (r.layer !== layer || netOf(d, r) === a.net) continue;
    for (let k = 1; k < r.points.length; k++) {
      const [u, v] = [r.points[k - 1], r.points[k]], rr = r.width / 2 + margin;
      mark(Math.min(u[0], v[0]) - rr, Math.min(u[1], v[1]) - rr, Math.max(u[0], v[0]) + rr, Math.max(u[1], v[1]) + rr, (x, y) => segmentDistance([x, y], [x, y], u, v) < rr);
    }
  }
  // Copper yields to a trace, but a trace through another net's pour can cut it in two, so
  // the router goes round one when it can and pays to cross when it must.
  const foreign = new Uint8Array(nx * ny);
  for (const q of d.pours ?? []) if (q.layer === layer && q.net !== a.net) { const xs = q.points.map(v => v[0]), ys = q.points.map(v => v[1]); mark(Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys), (x, y) => { if (inPolygon(q.points, x, y)) foreign[Math.round(y / g) * nx + Math.round(x / g)] = 1; return false; }); }
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
      const turn = k === dir ? 0 : Math.min((k - dir + 8) % 8, (dir - k + 8) % 8) > 1 ? .9 : .3, c = cost[s] + (k % 2 ? Math.SQRT2 : 1) * (foreign[jj * nx + ii] ? 6 : 1) + turn, t = (jj * nx + ii) * 8 + k;
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
  const nets = ["", ...[...new Set(pads(d).map(p => p.net).filter(n => n && n !== NO_CONNECT))].sort()], netId = name => Math.max(0, nets.indexOf(name ?? ""));
  const out = ["(kicad_pcb (version 20240108) (generator \"emergence\") (generator_version \"1.0\")",
    `  (general (thickness ${n(d.stackup.thickness)}) (legacy_teardrops no))`, "  (paper \"A4\")",
    "  (layers (0 \"F.Cu\" signal) (31 \"B.Cu\" signal) (38 \"B.Mask\" user) (39 \"F.Mask\" user) (36 \"B.SilkS\" user \"B.Silkscreen\") (37 \"F.SilkS\" user \"F.Silkscreen\") (44 \"Edge.Cuts\" user) (47 \"F.CrtYd\" user \"F.Courtyard\") (49 \"F.Fab\" user))",
    "  (setup (pad_to_mask_clearance 0))", ...nets.map((name, i) => `  (net ${i} ${q(name)})`)];
  const byComponent = Map.groupBy ? Map.groupBy(pads(d), p => p.component) : pads(d).reduce((m, p) => m.set(p.component, [...(m.get(p.component) || []), p]), new Map());
  for (const c of d.components) {
    if (c.footprint === "via") continue;
    const b = bodyBounds(c, d), font = "(effects (font (size 1 1) (thickness 0.15)))";
    out.push(`  (footprint ${q("Emergence:" + c.footprint)} (layer "F.Cu") (uuid ${q(uuid())}) (at ${n(c.x)} ${n(c.y)})`,
      `    (property "Reference" ${q(c.id)} (at 0 ${n(-b.height / 2 - 1.2)} 0) (layer "F.Fab") (uuid ${q(uuid())}) ${font})`,
      `    (property "Value" ${q(c.name)} (at 0 ${n(b.height / 2 + 1.2)} 0) (layer "F.Fab") (uuid ${q(uuid())}) ${font})`,
      `    (fp_rect (start ${n(b.x - c.x)} ${n(b.y - c.y)}) (end ${n(b.x + b.width - c.x)} ${n(b.y + b.height - c.y)}) (stroke (width 0.05) (type solid)) (fill none) (layer "F.CrtYd") (uuid ${q(uuid())}))`);
    for (const p of byComponent.get(c.id) || [])
      if (p.smd) out.push(`    (pad ${q(p.number)} smd rect (at ${n(p.x - c.x)} ${n(p.y - c.y)}) (size ${n(p.w)} ${n(p.h)}) (layers "F.Cu" "F.Paste" "F.Mask") (net ${netId(p.net)} ${q(p.net || "")}) (uuid ${q(uuid())}))`);
      else out.push(`    (pad ${q(p.number)} thru_hole circle (at ${n(p.x - c.x)} ${n(p.y - c.y)}) (size ${n(2 * p.radius)} ${n(2 * p.radius)}) (drill ${p.drill ?? DRILL}) (layers "*.Cu" "*.Mask") (net ${netId(p.net)} ${q(p.net || "")}) (uuid ${q(uuid())}))`);
    out.push("  )");
  }
  for (const r of d.routes) for (let k = 1; k < r.points.length; k++)
    out.push(`  (segment (start ${n(r.points[k - 1][0])} ${n(r.points[k - 1][1])}) (end ${n(r.points[k][0])} ${n(r.points[k][1])}) (width ${n(r.width)}) (layer ${r.layer === "top" ? '"F.Cu"' : '"B.Cu"'}) (net ${netId(netOf(d, r))}) (uuid ${q(uuid())}))`);
  // A via is a via to KiCad, not a part to place.
  for (const c of d.components) if (c.footprint === "via") for (const p of byComponent.get(c.id) || []) out.push(`  (via (at ${n(p.x)} ${n(p.y)}) (size ${n(2 * p.radius)}) (drill ${p.drill ?? DRILL}) (layers "F.Cu" "B.Cu") (net ${netId(p.net)}) (uuid ${q(uuid())}))`);
  // Pours as zones: KiCad refills them against its own design rules.
  for (const p of d.pours ?? []) out.push(`  (zone (net ${netId(p.net)}) (net_name ${q(p.net)}) (layer ${p.layer === "top" ? '"F.Cu"' : '"B.Cu"'}) (uuid ${q(uuid())}) (hatch edge 0.5) (connect_pads yes (clearance ${n(d.rules.clearance)})) (min_thickness ${n(d.rules.minWidth)}) (filled_areas_thickness no) (fill yes (thermal_gap 0.5) (thermal_bridge_width 0.5)) (polygon (pts ${p.points.map(([x, y]) => `(xy ${n(x)} ${n(y)})`).join(" ")})))`);
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
// Even-odd fill of the read-back loops against a raster of the drawn copper. Both
// are scan-converted, so a board's cost follows its copper rather than cells × edges.
export function gerberCheck(d, gerber, layer, cell = .1) {
  const loops = parseGerber(gerber), nx = Math.max(1, Math.ceil(d.width / cell - .5)), ny = Math.max(1, Math.ceil(d.height / cell - .5)), X = i => (i + .5) * cell, Y = j => (j + .5) * cell;
  const read = new Uint8Array(nx * ny), rows = Array.from({ length: ny }, () => []);
  for (const l of loops) for (let i = 0, j = l.length - 1; i < l.length; j = i++) {
    const [xi, yi] = l[i], [xj, yj] = l[j], from = Math.max(0, Math.floor((d.height - Math.max(yi, yj)) / cell - .5)), to = Math.min(ny - 1, Math.ceil((d.height - Math.min(yi, yj)) / cell - .5));
    for (let r = from; r <= to; r++) { const py = d.height - Y(r); if ((yi > py) !== (yj > py)) rows[r].push((xj - xi) * (py - yi) / (yj - yi) + xi); }
  }
  rows.forEach((c, r) => { c.sort((a, b) => a - b); let k = 0; for (let i = 0; i < nx; i++) { while (k < c.length && c[k] <= X(i)) k++; read[r * nx + i] = (c.length - k) & 1; } });
  // Stamp a box of the raster wherever a test holds.
  const stamp = (grid, x0, y0, x1, y1, test, value = 1) => {
    for (let j = Math.max(0, Math.floor(y0 / cell - .5)); j <= Math.min(ny - 1, Math.ceil(y1 / cell - .5)); j++) for (let i = Math.max(0, Math.floor(x0 / cell - .5)); i <= Math.min(nx - 1, Math.ceil(x1 / cell - .5)); i++)
      if (test(X(i), Y(j), j * nx + i)) grid[j * nx + i] = value;
  };
  const ps = pads(d).filter(p => !p.smd || layer === "top"), nets = new Map(pads(d).map(p => [p.id, p.net])), rs = d.routes.filter(r => r.layer === layer), c = d.rules.clearance;
  const pad = (grid, p, grow) => { const r = (p.w ? Math.hypot(p.w, p.h) / 2 : p.radius) + grow; stamp(grid, p.x - r, p.y - r, p.x + r, p.y + r, (x, y) => padDistance(p, x, y) <= grow); };
  const track = (grid, r, grow) => r.points.slice(1).forEach((v, k) => { const u = r.points[k], w = r.width / 2 + grow; stamp(grid, Math.min(u[0], v[0]) - w, Math.min(u[1], v[1]) - w, Math.max(u[0], v[0]) + w, Math.max(u[1], v[1]) + w, (x, y) => segmentDistance([x, y], [x, y], u, v) <= w); });
  const drawn = new Uint8Array(nx * ny);
  for (const p of (d.pours ?? []).filter(p => p.layer === layer)) {
    const clear = new Uint8Array(nx * ny), xs = p.points.map(q => q[0]), ys = p.points.map(q => q[1]);
    for (const q of ps) if (q.net !== p.net) pad(clear, q, c);
    for (const r of rs) if (nets.get(r.from) !== p.net) track(clear, r, c);
    stamp(drawn, Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys), (x, y, at) => !clear[at] && inPolygon(p.points, x, y));
  }
  for (const p of ps) pad(drawn, p, 0);
  for (const r of rs) track(drawn, r, 0);
  for (const p of pads(d).filter(p => !p.smd)) { const r = (p.drill ?? DRILL) / 2; stamp(drawn, p.x - r, p.y - r, p.x + r, p.y + r, (x, y) => Math.hypot(x - p.x, y - p.y) < r, 0); }
  let both = 0, onlyG = 0, onlyD = 0;
  for (let k = 0; k < read.length; k++) if (read[k] && drawn[k]) both++; else if (read[k]) onlyG++; else if (drawn[k]) onlyD++;
  const area = n => n * cell * cell;
  return { loops: loops.length, drawnArea: area(both + onlyD), gerberArea: area(both + onlyG), mismatch: area(onlyG + onlyD), agreement: both / Math.max(1, both + onlyG + onlyD) };
}
export const pcbAssumptions = d => [
  { id: "pcb.epsr", label: `Dielectric εr ${d.stackup.epsr} (FR-4 nominal)`, kind: "nominal", source: "Real laminate runs 4.2–4.8 with resin content and weave; controlled impedance needs the fabricator's value." },
  { id: "pcb.footprints", label: "Footprints and pad sizes", kind: "illustrative", source: "Estimates, not IPC-7351 land patterns. Replace with the manufacturer's package drawing." },
  { id: "pcb.currents", label: "Every trace carries its net's full RMS current, unless the net has a pour, which then carries it", kind: "nominal", source: "Conservative on bare traces. A trace on a poured net counts as a stub at zero amps, so a heavy path beside a pour is not flagged. Branch currents are not solved." },
  { id: "pcb.pour", label: "Pour loss estimated as a uniform strip; current in a pour is not solved", kind: "nominal", source: "A field solve of the drawn copper, or a thermal image, shows where it crowds." },
  { id: "pcb.drill", label: `${DRILL} mm finished drill on through pads unless a footprint says otherwise`, kind: "illustrative", source: "Sets via barrels in the drill file and kernel board." },
];
export { FOOTPRINTS, pcbChecks };
