// Mechanical commands over the revisioned parameter/feature/drawing state.
// The browser workspace and the headless agent pass a host that owns the
// state; both get the same commands, ports and Rust evaluation.
import { FIELDS, DEFAULTS, values, validParams, validFeatures, validDrawings } from "../repository.js";
import { drawingDefaults } from "../drawings.js";
import { assemblyParts } from "../scene.js";
export const PRESETS = {
  reference: {}, compact: { width: 32, length: 35, finHeight: 16, fins: 10, velocity: 2, current: 20 },
  optical: { width: 55, length: 65, fins: 12, opticalPower: 4, absorption: 65, current: 20 },
};
const UNIT = { "": "", "m/s": "m/s", A: "A", "°C": "degC", W: "W", mm: "mm", "%": "%" };
export const FIELD_UNIT = key => key === "fins" ? "count" : UNIT[FIELDS[key].unit] ?? "";
export class MechanicalState {
  constructor() { this.params = values({}); this.features = []; this.drawings = {}; }
  state() { return { params: this.params, features: this.features, drawings: this.drawings }; }
  snapshot() { return structuredClone(this.state()); }
  restore(s) { ({ params: this.params, features: this.features, drawings: this.drawings } = structuredClone(s)); this.changed?.(); }
  set(next) {
    if (!validParams(next.params) || !validFeatures(next.features) || !validDrawings(next.drawings)) throw new Error("That change leaves the mechanical design outside its supported ranges.");
    this.params = next.params; this.features = next.features; this.drawings = next.drawings; this.changed?.();
  }
}
// Released drawings return to review whenever the geometry they show changes.
const review = (drawings, id) => { for (const [k, d] of Object.entries(drawings)) if ((!id || k === id) && d.status === "Released") drawings[k] = { ...d, status: "For review" }; };
const DRAWING_FIELDS = ["number", "title", "revision", "material", "finish", "tolerance", "projection", "prepared", "checked", "status", "notes"];
export function mechanicalCommands(host) {
  const edit = fn => { const s = structuredClone(host.state()); const r = fn(s); host.set(s); return r; };
  const part = (s, id) => { const p = s.features.find(f => f.id === id); if (!p) throw new Error(`No custom part ${id}. Parts: ${s.features.map(f => f.id).join(", ") || "none"}.`); return p; };
  const dims = [["x", -100, 100], ["y", -100, 100], ["z", -100, 100], ["width", 1, 50], ["depth", 1, 50], ["height", 1, 50]];
  return [
    { op: "mech.set", title: "Set a mechanical parameter", summary: "Set one parameter of the heat sink, airflow, electrical or optical inputs.", example: "mech.set fins 11", coalesce: a => "mech.set:" + a.name,
      args: [{ name: "name", type: "enum", values: Object.keys(FIELDS), required: true }, { name: "value", type: "number", required: true }],
      run: ({ name, value }) => edit(s => { const f = FIELDS[name]; if (value < f.min || value > f.max || name === "fins" && !Number.isInteger(value)) throw new Error(`${f.label} must be ${f.min}–${f.max} ${f.unit}.`); s.params[name] = value; review(s.drawings); return { message: `${f.label} = ${value} ${f.unit}` }; }) },
    { op: "mech.preset", title: "Apply a mechanical preset", args: [{ name: "name", type: "enum", values: Object.keys(PRESETS), required: true }],
      run: ({ name }) => edit(s => { s.params = { ...values({}), ...PRESETS[name] }; s.features = []; s.drawings = {}; return { message: `Applied the ${name} preset as a draft.` }; }) },
    { op: "mech.add_part", title: "Sketch and extrude a part", summary: "Origin-centred rectangle or circle on XY, extruded along Z.", example: 'mech.add_part "Mounting block" box width=12 depth=8 height=6',
      args: [{ name: "name", type: "string", maxLength: 60, required: true }, { name: "shape", type: "enum", values: ["box", "cylinder"], default: "box" },
        { name: "width", type: "number", min: 1, max: 50, unit: "mm", default: 12 }, { name: "depth", type: "number", min: 1, max: 50, unit: "mm", default: 8, named: true },
        { name: "height", type: "number", min: 1, max: 50, unit: "mm", default: 6, named: true }, { name: "x", type: "number", min: -100, max: 100, unit: "mm", named: true },
        { name: "y", type: "number", min: -100, max: 100, unit: "mm", named: true, default: 0 }, { name: "z", type: "number", min: -100, max: 100, unit: "mm", named: true, default: 0 },
        { name: "domain", type: "enum", values: ["Mechanical", "Electrical", "Thermal", "Photonics", "Radiation"], default: "Mechanical", named: true }],
      run: a => edit(s => { if (s.features.length >= 16) throw new Error("The prototype supports sixteen custom solids.");
        const p = { id: `part-${crypto.randomUUID().slice(0, 12)}`, name: a.name.trim() || "Part", shape: a.shape, domain: a.domain, x: a.x ?? Math.min(100, s.params.width / 2 + 20), y: a.y, z: a.z,
          width: a.width, depth: a.shape === "cylinder" ? a.width : a.depth, height: a.height };
        s.features.push(p); return { id: p.id, message: `Added ${p.name} (${p.id}).` }; }) },
    { op: "mech.update_part", title: "Edit a part's placement or dimensions", example: "mech.update_part part-1 x=40 height=8", coalesce: a => "mech.update_part:" + a.id + Object.keys(a).sort().join(),
      args: [{ name: "id", type: "string", required: true }, { name: "name", type: "string", maxLength: 60, named: true }, ...dims.map(([name, min, max]) => ({ name, type: "number", min, max, unit: "mm", named: true }))],
      run: ({ id, ...changes }) => edit(s => { const p = part(s, id);
        for (const [k, v] of Object.entries(changes)) if (v !== undefined) p[k] = k === "name" ? (v.trim() || p.name) : v;
        if (p.shape === "cylinder") p.depth = p.width; review(s.drawings, id); return { message: `Updated ${p.name}.` }; }) },
    { op: "mech.duplicate_part", title: "Duplicate a part", args: [{ name: "id", type: "string", required: true }],
      run: ({ id }) => edit(s => { if (s.features.length >= 16) throw new Error("The prototype supports sixteen custom solids."); const p = part(s, id);
        const copy = { ...p, id: "part-" + crypto.randomUUID().slice(0, 12), name: (p.name + " copy").slice(0, 60), x: Math.min(100, p.x + 5) }; s.features.push(copy); return { id: copy.id, message: `Duplicated ${p.name}.` }; }) },
    { op: "mech.remove_part", title: "Delete a part", args: [{ name: "id", type: "string", required: true }],
      run: ({ id }) => edit(s => { const p = part(s, id); s.features = s.features.filter(f => f.id !== id); delete s.drawings[id]; return { message: `Deleted ${p.name}.` }; }) },
    { op: "mech.drawing", title: "Set a shop-drawing field", example: 'mech.drawing sink material "6061-T6 aluminium"',
      args: [{ name: "part", type: "string", required: true }, { name: "field", type: "enum", values: DRAWING_FIELDS, required: true }, { name: "value", type: "string", maxLength: 600, required: true }],
      run: ({ part: id, field, value }) => edit(s => {
        const target = assemblyParts(s.params, s.features).find(p => p.id === id && p.shape !== "context");
        if (!target) throw new Error(`No drawable part ${id}.`);
        const next = { ...drawingDefaults(target), ...(s.drawings[id] || {}), [field]: value };
        if (field !== "status" && next.status === "Released") next.status = "For review";
        s.drawings[id] = next; return { message: `${id} ${field} updated.` }; }) },
  ];
}
// Typed ports onto the original electrothermal project graph in Rust.
export function mechanicalAdapter(host, request) {
  const inputs = Object.entries(FIELDS).filter(([k]) => k !== "opticalSize").map(([id, f]) => ({ id, label: f.label, unit: FIELD_UNIT(id), min: f.min, max: f.max, base: () => host.state().params[id] }));
  const outputs = [
    ["q1_power", "Q1 dissipation", "W"], ["q1_tj", "Q1 junction temperature", "degC"], ["sink_resistance", "Sink thermal resistance", "K/W"],
    ["trace_resistance", "Power trace resistance", "ohm"], ["sink_dp", "Sink pressure drop (Rust graph)", "Pa"], ["gap", "Fin channel gap", "mm"],
    ["sink_volume", "Sink metal volume", "mm3"], ["current_out", "MOSFET RMS current", "A"], ["width_out", "Sink width", "mm"], ["length_out", "Sink length", "mm"],
    ["finHeight_out", "Fin height", "mm"], ["fins_out", "Fin count", "count"], ["finThickness_out", "Fin thickness", "mm"], ["sink_mass", "Heat-sink mass (6063 Al)", "kg"],
  ].map(([id, label, unit]) => ({ id, label, unit }));
  return { title: "Mechanical", ports: () => ({ inputs, outputs }),
    async evaluate(doc, overrides) {
      const p = { ...(doc?.params ?? host.state().params), ...overrides };
      if (overrides.fins !== undefined) p.fins = Math.round(overrides.fins);
      const { result: r } = await request("workspace", p);
      const node = id => r.nodes[id] ? { value: r.nodes[id].nominal, lo: r.nodes[id].lo, hi: r.nodes[id].hi } : { value: NaN };
      const tj = node("q1.tj");
      const names = { "REQ-TJ": "Junction temperature ≤ 100 °C", "REQ-DROP": "PCB trace drop ≤ 50 mV", "REQ-FAN": "Heat-sink pressure drop ≤ 20 Pa" };
      return { outputs: {
        q1_power: { ...node("q1.power"), unit: "W" }, q1_tj: { value: tj.value - 273.15, lo: tj.lo - 273.15, hi: tj.hi - 273.15, unit: "degC" },
        sink_resistance: { ...node("sink.thermal_resistance"), unit: "K/W" }, trace_resistance: { ...node("trace.resistance"), unit: "ohm" },
        sink_dp: { ...node("sink.pressure_drop"), unit: "Pa" }, gap: { value: r.handoffs.gap_mm, unit: "mm" }, sink_volume: { value: r.handoffs.sink_volume_mm3, unit: "mm3" },
        current_out: { value: p.current, unit: "A" }, width_out: { value: p.width, unit: "mm" }, length_out: { value: p.length, unit: "mm" },
        finHeight_out: { value: p.finHeight, unit: "mm" }, fins_out: { value: p.fins, unit: "count" }, finThickness_out: { value: p.finThickness, unit: "mm" },
        sink_mass: { value: r.handoffs.sink_volume_mm3 * 2.70e-6, unit: "kg" } },
        // The graph reports SI: kelvin and volts. Show °C and mV as the requirements state them.
        checks: r.checks.map(c => { const [offset, scale, unit] = { "REQ-TJ": [273.15, 1, "degC"], "REQ-DROP": [0, 1000, "mV"], "REQ-FAN": [0, 1, "Pa"] }[c.id] || [0, 1, c.unit];
          return { id: c.id, text: names[c.id] || c.text, status: c.status, reason: c.reason, unit, ...(c.value && { value: { nominal: (c.value.nominal - offset) * scale, lo: (c.value.lo - offset) * scale, hi: (c.value.hi - offset) * scale } }) }; }), details: r, params: p,
        certificates: [{ label: "Interval contains nominal for every shared node", ok: Object.values(r.nodes).every(n => !n || n.lo <= n.nominal + 1e-12 && n.nominal <= n.hi + 1e-12) }] };
    } };
}
export const mechanicalAssumptions = () => [
  { id: "mech.density", label: "Heat-sink density 2.70 g/cm³ (6063 aluminium)", kind: "nominal", source: "Handbook value; anodising and fasteners not included." },
  { id: "mech.channel", label: "Laminar developing channel correlations for the fin array", kind: "nominal", source: "Rust project graph; out-of-regime inputs report unavailable." },
  { id: "mech.switching", label: "Q1 switching loss 3 W ± 30%", kind: "nominal", source: "Prescribed loss interval from the motor-driver project." },
];
export { DEFAULTS };
