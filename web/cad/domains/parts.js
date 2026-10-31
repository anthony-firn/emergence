// Solid components as feature trees: the JSON model format of `fluid cad`,
// built by the Rust kernel (cad::model) in the browser. Every number may be an
// expression over the part's parameters; one model is a family of parts.
// Templates are the repository's worked models, rescaled to millimetres.
export const MATERIALS = {
  "al-6061": { title: "Aluminium 6061-T6", density: 2.70 }, "steel-304": { title: "Stainless steel 304", density: 8.00 },
  "brass": { title: "Brass C360", density: 8.50 }, "copper": { title: "Copper C110", density: 8.96 }, "ti-64": { title: "Titanium Ti-6Al-4V", density: 4.43 },
  "cfrp": { title: "Carbon fibre / epoxy laminate", density: 1.55 }, "abs": { title: "ABS", density: 1.04 }, "peek": { title: "PEEK", density: 1.32 }, "pa12": { title: "Nylon PA12 (printed)", density: 1.01 },
};
// Two plain starts for a part built feature by feature; the rest are finished models.
export const BASE_TEMPLATES = ["block", "disc"];
export const TEMPLATES = {
  block: { title: "Block (build from here)", material: "al-6061", model: { name: "block", cell: .5, params: { w: 60, d: 40, h: 20 }, shape: { op: "box", center: [0, 0, "h/2"], half: ["w/2", "d/2", "h/2"] } } },
  disc: { title: "Disc (build from here)", material: "al-6061", model: { name: "disc", cell: .5, params: { diameter: 60, h: 12 }, shape: { op: "cylinder_z", center: [0, 0, "h/2"], radius: "diameter/2", height: "h" } } },
  flange: { title: "Bolted flange", material: "steel-304", model: { name: "bolted flange", cell: .6,
    params: { bolts: 6, bolt_circle: 80, bolt: 3.3, thickness: 6, bore: 25, od: "bolt_circle + 5*bolt", hub: 12 },
    shape: { op: "difference",
      a: { op: "union", parts: [{ op: "cylinder_z", center: [0, 0, "thickness/2"], radius: "od/2", height: "thickness" }, { op: "cylinder_z", center: [0, 0, "thickness + hub/2"], radius: "bore/2 + 5", height: "hub" }] },
      b: { op: "union", parts: [{ op: "cylinder_z", center: [0, 0, 0], radius: "bore/2", height: "4*(thickness + hub)" },
        { op: "polar_array", axis: [0, 0, 1], count: "bolts", shape: { op: "cylinder_z", center: ["bolt_circle/2", 0, "thickness/2"], radius: "bolt", height: "4*thickness" } }] } } } },
  gear: { title: "Involute spur gear (ISO 21771)", material: "brass", model: { name: "spur gear", cell: .2,
    params: { module: 2, teeth: 24, face: 8, bore: 10 },
    shape: { op: "kit", kit: "spur_gear", args: { module: "module", teeth: "teeth", thickness: "face", bore: "bore", pressure_angle: .3490658503988659, flank_points: 24 } } } },
  bracket: { title: "Formed sheet-metal bracket", material: "al-6061", model: { name: "sheet bracket", cell: .45,
    params: { t: 1.5, radius: 3, width: 40, web: 50, flange: 25 },
    shape: { op: "kit", kit: "sheet_part", args: { thickness: "t", k_factor: .42, width: "width", flanges: ["flange", "web", "flange"],
      bends: [{ angle: 1.5707963267948966, inside_radius: "radius" }, { angle: 1.5707963267948966, inside_radius: "radius" }], per_bend: 16 } } } },
  elbow: { title: "NPS 2 Sch 40 elbow run (ASME B36.10M)", material: "steel-304", model: { name: "elbow run", cell: 1.4,
    params: { od: 60.3, wall: 3.91, leg_length: 120, bend: "1.5 * od" },
    shape: { op: "kit", kit: "pipe_run", args: { pipe: { outer: "od", wall: "wall" }, part: "wall", sides: 40,
      run: { start: [0, 0, 0], direction: [1, 0, 0], legs: [{ leg: "straight", length: "leg_length" }, { leg: "elbow", radius: "bend", angle: "pi/2", axis: [0, 0, 1] }, { leg: "straight", length: "leg_length" }] } } } } },
  pinfin: { title: "Pin-fin heat sink", material: "al-6061", model: { name: "pin-fin heat sink", cell: .4,
    params: { base: 40, plate: 4, pin: 2.5, height: 20, pitch: 6, count: 6 },
    shape: { op: "union", parts: [{ op: "box", center: [0, 0, "plate/2"], half: ["base/2", "base/2", "plate/2"] },
      { op: "array", shape: { op: "cylinder_z", center: [0, 0, "plate + height/2"], radius: "pin/2", height: "height" }, period: ["pitch", "pitch", 0], counts: ["count", "count", 1] }] } } },
  thread: { title: "ISO metric thread", material: "steel-304", model: { name: "metric thread", cell: .1,
    params: { diameter: 8, pitch: 1.25, length: 20 }, shape: { op: "kit", kit: "thread", args: { diameter: "diameter", pitch: "pitch", length: "length" } } } },
  spring: { title: "Compression spring", material: "steel-304", model: { name: "compression spring", cell: .2,
    params: { wire: 1.2, mean_diameter: 10, free_length: 30, coils: 6 }, shape: { op: "kit", kit: "spring", args: { wire: "wire", mean_diameter: "mean_diameter", free_length: "free_length", coils: "coils", steps: 240 } } } },
  gyroid: { title: "Gyroid sheet lattice", material: "pa12", model: { name: "gyroid lattice", cell: .3,
    params: { cell: 8, wall: .6, nx: 3, ny: 3, nz: 2 },
    shape: { op: "clip", shape: { op: "periodic", kind: "gyroid", cell: ["cell", "cell", "cell"], level: 0, thickness: "wall" }, min: [0, 0, 0], max: ["nx*cell", "ny*cell", "nz*cell"] } } },
  wing: { title: "NACA 2412 wing section", material: "abs", model: { name: "NACA 2412 wing", cell: 1.8,
    params: { chord: 200, span: 600 },
    shape: { op: "kit", kit: "wing", args: { camber: .02, camber_position: .4, thickness: .12, chord: "chord", closed: true, points: 100, span: "span" } } } },
  arm: { title: "Carbon arm tube", material: "cfrp", model: { name: "arm tube", cell: .9,
    params: { od: 25, wall: 1.5, start: 70, length: 310 },
    shape: { op: "kit", kit: "pipe_run", args: { pipe: { outer: "od", wall: "wall" }, part: "wall", sides: 32, run: { start: ["start", 0, 0], direction: [1, 0, 0], legs: [{ leg: "straight", length: "length" }] } } } } },
  plate: { title: "Centre plate with lightening holes", material: "cfrp", model: { name: "centre plate", cell: .9,
    params: { size: 190, t: 2.5, hole: 22, pitch: 48, holes: 3 },
    shape: { op: "difference", a: { op: "box", center: [0, 0, "t/2"], half: ["size/2", "size/2", "t/2"], round: 12 },
      b: { op: "array", shape: { op: "cylinder_z", center: [0, 0, "t/2"], radius: "hole/2", height: "4*t" }, period: ["pitch", "pitch", 0], counts: ["holes", "holes", 1] } } } },
  mount: { title: "Motor mount", material: "al-6061", model: { name: "motor mount", cell: .35,
    params: { bolts: 4, bolt_circle: 25, bolt: 1.6, thickness: 3, bore: 8, od: "bolt_circle + 6*bolt", hub: 5 },
    shape: { op: "difference",
      a: { op: "union", parts: [{ op: "cylinder_z", center: [0, 0, "thickness/2"], radius: "od/2", height: "thickness" }, { op: "cylinder_z", center: [0, 0, "thickness + hub/2"], radius: "bore/2 + 3", height: "hub" }] },
      b: { op: "union", parts: [{ op: "cylinder_z", center: [0, 0, 0], radius: "bore/2", height: "4*(thickness + hub)" },
        { op: "polar_array", axis: [0, 0, 1], count: "bolts", shape: { op: "cylinder_z", center: ["bolt_circle/2", 0, "thickness/2"], radius: "bolt", height: "4*thickness" } }] } } } },
  enclosure: { title: "Shelled enclosure with vents", material: "abs", model: { name: "enclosure", cell: .7,
    params: { w: 90, d: 70, h: 40, wall: 2, round: 4, vent: 3, vents: 6 },
    shape: { op: "difference", a: { op: "shell", thickness: "wall", shape: { op: "box", center: [0, 0, "h/2"], half: ["w/2", "d/2", "h/2"], round: "round" } },
      b: { op: "union", parts: [{ op: "box", center: [0, 0, "h"], half: ["w/2 - wall", "d/2 - wall", "wall"] },
        { op: "array", shape: { op: "box", center: [0, "-d/2", "h/2"], half: ["vent/2", "wall*2", "h/4"] }, period: ["vent*2.5", 0, 0], counts: ["vents", 1, 1] }] } } } },
};
// The starter: one bolted flange. Projects bring their own parts lists.
export function partsExample() {
  return { parts: [{ id: "flange1", name: TEMPLATES.flange.title, template: "flange", model: structuredClone(TEMPLATES.flange.model), overrides: {}, x: 0, y: 0, z: 0, rz: 0, count: 1, pattern: "single", material: TEMPLATES.flange.material, domain: "Mechanical" }] };
}
const finite = (n, lo, hi) => Number.isFinite(n) && n >= lo && n <= hi;
function refuseFiles(v) {
  if (Array.isArray(v)) return v.forEach(refuseFiles);
  if (v && typeof v === "object") { if (v.op === "mesh") throw new Error("File-backed mesh features are not available in the browser."); Object.values(v).forEach(refuseFiles); }
}
export function validateParts(d) {
  if (!d || !Array.isArray(d.parts) || d.parts.length > 24 || new Set(d.parts.map(p => p?.id)).size !== d.parts.length) throw new Error("Invalid parts list.");
  return { parts: d.parts.map(p => {
    if (!p || typeof p.id !== "string" || !/^[a-zA-Z][a-zA-Z0-9_-]{0,23}$/.test(p.id) || typeof p.name !== "string" || !p.name.trim() || p.name.length > 60
      || !(p.template === "custom" || TEMPLATES[p.template]) || !MATERIALS[p.material] || !["x", "y", "z"].every(k => finite(p[k], -2000, 2000))
      || typeof p.domain !== "string" || p.domain.length > 24 || !finite(p.rz ?? 0, -360, 360) || !Number.isInteger(p.count ?? 1) || (p.count ?? 1) < 1 || (p.count ?? 1) > 64
      || !["single", "polar", "stack"].includes(p.pattern ?? "single") || !finite(p.stack ?? 0, -2000, 2000)) throw new Error(`Invalid part ${p?.id ?? ""}.`);
    const m = p.model;
    if (!m || typeof m !== "object" || !m.shape || typeof m.shape !== "object" || !m.params || typeof m.params !== "object" || Array.isArray(m.params) || JSON.stringify(m).length > 24000)
      throw new Error(`${p.name}: a model needs params and a shape, under 24 KB.`);
    refuseFiles(m.shape);
    for (const [k, v] of Object.entries(m.params)) if (!/^[a-zA-Z_][a-zA-Z0-9_]{0,31}$/.test(k) || !(typeof v === "number" && Number.isFinite(v) || typeof v === "string" && v.length <= 200)) throw new Error(`${p.name}: parameter ${k}.`);
    if (m.cell !== undefined && !finite(m.cell, .01, 50)) throw new Error(`${p.name}: mesh cell.`);
    const overrides = p.overrides && typeof p.overrides === "object" ? p.overrides : {};
    for (const [k, v] of Object.entries(overrides)) if (!Object.hasOwn(m.params, k) || !Number.isFinite(v)) throw new Error(`${p.name}: override ${k} is not a numeric model parameter.`);
    return { id: p.id, name: p.name, template: p.template, model: { name: typeof m.name === "string" ? m.name.slice(0, 60) : "", params: { ...m.params }, shape: structuredClone(m.shape), ...(m.cell !== undefined && { cell: m.cell }) },
      overrides: { ...overrides }, x: p.x, y: p.y, z: p.z, rz: p.rz ?? 0, count: p.count ?? 1, pattern: p.pattern ?? "single", stack: p.stack ?? 0, material: p.material, domain: p.domain };
  }) };
}
// Placement of every instance: polar patterns rotate about z, stacks repeat along z.
export function instances(p) {
  return Array.from({ length: p.count }, (_, i) => p.pattern === "polar" ? { rz: p.rz + 360 * i / p.count, by: [p.x, p.y, p.z] } : { rz: p.rz, by: [p.x, p.y, p.z + (p.pattern === "stack" ? i * p.stack : 0)] });
}
export const placed = (p, inst = instances(p)[0]) => ({ ...p.model, shape: { op: "rotate", axis: [0, 0, 1], angle: inst.rz * Math.PI / 180, shape: { op: "translate", by: inst.by, shape: p.model.shape } } });
// Building a solid feature by feature. Each step wraps the tree so far, so the
// feature tree reads as the order a person worked in. Every number may be an
// expression of the part's parameters. Positions are the feature's centre, mm.
const half = v => typeof v === "number" ? v / 2 : `(${v})/2`;
export const FEATURE_SHAPES = ["box", "cylinder", "sphere", "cone", "torus", "extrude", "revolve"];
export function primitive({ shape, x = 0, y = 0, z = 0, width = 20, depth = 20, height = 20, radius = 10, radius2 = 5, round = 0, points, angle }) {
  const c = [x, y, z], at = v => typeof z === "number" && typeof v === "number" ? z + v : `(${z})+(${v})`;
  if (shape === "box") return { op: "box", center: c, half: [half(width), half(depth), half(height)], ...(round && { round }) };
  if (shape === "cylinder") return { op: "cylinder_z", center: c, radius, height };
  if (shape === "sphere") return { op: "sphere", center: c, radius };
  if (shape === "cone") return { op: "cone", a: [x, y, at(typeof height === "number" ? -height / 2 : `-(${height})/2`)], b: [x, y, at(half(height))], radius_a: radius, radius_b: radius2 };
  if (shape === "torus") return { op: "torus", center: c, major: radius, minor: radius2 };
  if (shape === "extrude" || shape === "revolve") {
    if (!Array.isArray(points) || points.length < 3 || !points.every(p => Array.isArray(p) && p.length === 2 && p.every(v => Number.isFinite(v) || typeof v === "string"))) throw new Error(`A ${shape} needs a sketch of at least three [x, y] points.`);
    const sketch = { op: "polygon", loops: [points] };
    return shape === "extrude" ? { op: "translate", by: [x, y, 0], shape: { op: "extrude", sketch, z0: at(typeof height === "number" ? -height / 2 : `-(${height})/2`), z1: at(half(height)) } }
      : { op: "translate", by: c, shape: { op: "revolve", sketch, ...(angle !== undefined && { angle: typeof angle === "number" ? angle * Math.PI / 180 : `(${angle})*pi/180` }) } };
  }
  throw new Error(`Shapes: ${FEATURE_SHAPES.join(", ")}.`);
}
export const COMBINE = ["add", "cut", "keep", "blend"];
export function combine(tree, how, other, k = 2) {
  if (how === "add") return tree.op === "union" ? { ...tree, parts: [...tree.parts, other] } : { op: "union", parts: [tree, other] };
  if (how === "cut") return { op: "difference", a: tree, b: other };
  if (how === "keep") return { op: "intersection", parts: [tree, other] };
  if (how === "blend") return { op: "blend", a: tree, b: other, k };
  throw new Error(`Combine by ${COMBINE.join(", ")}.`);
}
export const MODIFY = ["shell", "offset", "round", "fillet", "pattern", "polar", "mirror", "move", "rotate", "scale", "clip"];
const AXES = { x: [1, 0, 0], y: [0, 1, 0], z: [0, 0, 1] };
export function modify(tree, m, cell) {
  const { op, value = 1, axis = "z", x = 0, y = 0, z = 0, nx = 1, ny = 1, nz = 1, count = 6, min, max } = m, a = AXES[axis];
  if (!a) throw new Error("Axes: x, y, z.");
  if (op === "shell") return { op: "shell", shape: tree, thickness: value };
  if (op === "offset") return { op: "offset", shape: tree, distance: value };
  if (op === "round") return { op: "round_convex", shape: tree, radius: value, cell: cell ?? .5 };
  if (op === "fillet") return { op: "round_concave", shape: tree, radius: value, cell: cell ?? .5 };
  if (op === "pattern") return { op: "array", shape: tree, period: [x, y, z], counts: [nx, ny, nz] };
  if (op === "polar") return { op: "polar_array", axis: a, count, shape: tree };
  if (op === "mirror") return { op: "union", parts: [tree, { op: "mirror", normal: a, shape: tree }] };
  if (op === "move") return { op: "translate", by: [x, y, z], shape: tree };
  if (op === "rotate") return { op: "rotate", axis: a, angle: typeof value === "number" ? value * Math.PI / 180 : `(${value})*pi/180`, shape: tree };
  if (op === "scale") return { op: "scale", factor: value, shape: tree };
  if (op === "clip") { if (!min || !max) throw new Error("A clip needs min and max corners."); return { op: "clip", shape: tree, min, max }; }
  throw new Error(`Operations: ${MODIFY.join(", ")}.`);
}
// The tree as lines a person can read, newest feature first.
export function outline(node, depth = 0, out = []) {
  if (!node || typeof node !== "object" || depth > 12) return out;
  const op = node.op; if (!op) return out;
  out.push({ depth, op, detail: op === "kit" ? node.kit : ["thickness", "distance", "radius", "count", "factor"].filter(k => node[k] !== undefined).map(k => `${k} ${node[k]}`).join(" · ") });
  for (const k of ["a", "b", "shape"]) if (node[k] && typeof node[k] === "object") outline(node[k], depth + 1, out);
  for (const n of node.parts ?? []) outline(n, depth + 1, out);
  return out;
}
export const partsAssumptions = d => d.parts.length ? [{ id: "parts.materials", label: "Handbook densities for mass", kind: "nominal", source: "Alloy temper and print infill change mass; printed parts are taken as fully dense." }] : [];
