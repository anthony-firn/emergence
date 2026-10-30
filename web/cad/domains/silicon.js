// Floorplan power, the chip-flow stage summary, and GDSII. The stream file is
// written per the GDSII Stream Format (Calma, release 6): 1 nm database unit,
// 1 µm user unit, one BOUNDARY per block on a layer by block kind, the core
// keepout and die outline, and TEXT labels. It is a concept floorplan, with no
// PDK layer map, and is read back to check every boundary before download.
export const GDS_LAYERS = { die: 0, logic: 1, memory: 2, io: 3, analog: 4, core: 10, text: 63 };
export const STAGES = [
  { id: "rtl", title: "RTL", detail: "Verilog subset parse" },
  { id: "synth", title: "Synthesis", detail: "Gate netlist on characterised cells" },
  { id: "equiv", title: "Equivalence", detail: "Netlist vs RTL on (input, state) points" },
  { id: "timing", title: "Timing", detail: "STA interpolated over 9 corners" },
  { id: "power", title: "Power & Tj", detail: "Switching + leakage at junction temperature" },
  { id: "floorplan", title: "Floorplan", detail: "Blocks, keepout, halos, density" },
  { id: "litho", title: "Lithography", detail: "Abbe aerial image of a layout clip" },
  { id: "route", title: "Place & route", detail: "Not in this browser flow", planned: true },
  { id: "signoff", title: "DRC / LVS sign-off", detail: "Requires a foundry PDK", planned: true },
];
// flows: { blockId: chip-flow result }; an RTL block's power replaces its estimate.
export function blockPowers(d, flows = {}) {
  return d.blocks.map(b => { const f = flows?.[b.id]; return { id: b.id, name: b.name, area: b.width * b.height * 1e-6, power: f?.power ? f.power.total : b.power / 1e3, fromRtl: !!f?.power }; });
}
export function powerSummary(d, flows) {
  const blocks = blockPowers(d, flows), total = blocks.reduce((s, b) => s + b.power, 0);
  const peak = blocks.reduce((m, b) => b.area > 0 && b.power / b.area > m.density ? { density: b.power / b.area, id: b.id } : m, { density: 0, id: null });
  return { blocks, total, peak, junction: d.chip.ambient + d.chip.theta * total };
}
// GDSII: big-endian records [length u16][type u8][datatype u8][payload].
const REC = { HEADER: [0x00, 2], BGNLIB: [0x01, 2], LIBNAME: [0x02, 6], UNITS: [0x03, 5], ENDLIB: [0x04, 0], BGNSTR: [0x05, 2], STRNAME: [0x06, 6], ENDSTR: [0x07, 0],
  BOUNDARY: [0x08, 0], TEXT: [0x0c, 0], LAYER: [0x0d, 2], DATATYPE: [0x0e, 2], XY: [0x10, 3], ENDEL: [0x11, 0], TEXTTYPE: [0x16, 2], STRING: [0x19, 6] };
function real8(v) {
  // Excess-64 base-16 floating point.
  if (v === 0) return new Uint8Array(8);
  const sign = v < 0 ? 0x80 : 0; v = Math.abs(v);
  let exp = 64; while (v >= 1) { v /= 16; exp++; } while (v < 1 / 16) { v *= 16; exp--; }
  const out = new Uint8Array(8); out[0] = sign | exp;
  let mant = v; for (let i = 1; i < 8; i++) { mant *= 256; const b = Math.floor(mant); out[i] = b; mant -= b; }
  return out;
}
export function gdsii(d, name = "FLOORPLAN") {
  const parts = [];
  const rec = (type, payload = new Uint8Array(0)) => { const [t, dt] = REC[type], b = new Uint8Array(4 + payload.length); new DataView(b.buffer).setUint16(0, b.length); b[2] = t; b[3] = dt; b.set(payload, 4); parts.push(b); };
  const i16 = (...v) => { const b = new Uint8Array(2 * v.length), dv = new DataView(b.buffer); v.forEach((x, i) => dv.setInt16(2 * i, x)); return b; };
  const i32 = v => { const b = new Uint8Array(4 * v.length), dv = new DataView(b.buffer); v.forEach((x, i) => dv.setInt32(4 * i, Math.round(x))); return b; };
  const str = s => { let b = new TextEncoder().encode(s.replace(/[^\x20-\x7e]/g, "_")); if (b.length % 2) b = Uint8Array.from([...b, 0]); return b; };
  const date = i16(2026, 1, 1, 0, 0, 0, 2026, 1, 1, 0, 0, 0);
  const box = (layer, x, y, w, h) => { rec("BOUNDARY"); rec("LAYER", i16(layer)); rec("DATATYPE", i16(0)); const nm = v => v * 1000;
    rec("XY", i32([nm(x), nm(y), nm(x + w), nm(y), nm(x + w), nm(y + h), nm(x), nm(y + h), nm(x), nm(y)])); rec("ENDEL"); };
  rec("HEADER", i16(600)); rec("BGNLIB", date); rec("LIBNAME", str("EMERGENCE.DB"));
  rec("UNITS", Uint8Array.from([...real8(1e-3), ...real8(1e-9)]));
  rec("BGNSTR", date); rec("STRNAME", str(name));
  // GDS y is up; the editor's y is down. Flip about the die height.
  const fy = (y, h) => d.height - y - h;
  box(GDS_LAYERS.die, 0, 0, d.width, d.height);
  box(GDS_LAYERS.core, d.margin, d.margin, d.width - 2 * d.margin, d.height - 2 * d.margin);
  for (const b of d.blocks) {
    box(GDS_LAYERS[b.kind], b.x, fy(b.y, b.height), b.width, b.height);
    rec("TEXT"); rec("LAYER", i16(GDS_LAYERS.text)); rec("TEXTTYPE", i16(0)); rec("XY", i32([(b.x + b.width / 2) * 1000, (fy(b.y, b.height) + b.height / 2) * 1000])); rec("STRING", str(b.name.slice(0, 60))); rec("ENDEL");
  }
  rec("ENDSTR"); rec("ENDLIB");
  const out = new Uint8Array(parts.reduce((s, p) => s + p.length, 0)); let o = 0; for (const p of parts) { out.set(p, o); o += p.length; }
  return out;
}
export function readGDS(bytes) {
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength), shapes = [], texts = [];
  let o = 0, cur = null, units = null;
  const real = at => { const s = dv.getUint8(at) & 0x80 ? -1 : 1, e = (dv.getUint8(at) & 0x7f) - 64; let m = 0; for (let i = 7; i >= 1; i--) m = (m + dv.getUint8(at + i)) / 256; return s * m * 16 ** e; };
  while (o + 4 <= bytes.length) {
    const len = dv.getUint16(o), type = dv.getUint8(o + 2), body = o + 4;
    if (len < 4) throw new Error("Corrupt GDSII record.");
    if (type === 0x03) units = [real(body), real(body + 8)];
    else if (type === 0x08 || type === 0x0c) cur = { kind: type === 0x08 ? "boundary" : "text" };
    else if (type === 0x0d && cur) cur.layer = dv.getInt16(body);
    else if (type === 0x10 && cur) cur.xy = Array.from({ length: (len - 4) / 4 }, (_, i) => dv.getInt32(body + 4 * i));
    else if (type === 0x19 && cur) cur.text = new TextDecoder().decode(bytes.subarray(body, o + len)).replace(/\0+$/, "");
    else if (type === 0x11 && cur) { (cur.kind === "boundary" ? shapes : texts).push(cur); cur = null; }
    else if (type === 0x04) break;
    o += len;
  }
  return { units, shapes, texts };
}
// Read the stream back and compare every block rectangle exactly (in nm).
export function verifyGDS(d, bytes) {
  const r = readGDS(bytes), problems = [];
  if (!r.units || Math.abs(r.units[0] - 1e-3) > 1e-12 || Math.abs(r.units[1] - 1e-9) > 1e-21) problems.push("Units did not round-trip.");
  for (const b of d.blocks) {
    const y0 = (d.height - b.y - b.height) * 1000, want = [b.x * 1000, y0, (b.x + b.width) * 1000, y0 + b.height * 1000].map(Math.round);
    const hit = r.shapes.find(s => s.layer === GDS_LAYERS[b.kind] && s.xy[0] === want[0] && s.xy[1] === want[1] && s.xy[4] === want[2] && s.xy[5] === want[3]);
    if (!hit) problems.push(`${b.name} boundary did not round-trip.`);
  }
  return { ok: !problems.length, problems, boundaries: r.shapes.length, labels: r.texts.length, bytes: bytes.length };
}
export const siliconAssumptions = d => [
  { id: "silicon.blocks", label: "Block sizes and non-RTL block power", kind: "illustrative", source: "Planning estimates; no PDK, macro abstracts or measured power." },
  { id: "silicon.cells", label: "150 nm EKV cells characterised from the repository's transistor law", kind: "nominal", source: "analysis/chip.md: ideal clock, wires as fixed capacitance, glitches excluded (−15% power)." },
  { id: "silicon.package", label: `Package θJA ${d.chip.theta} K/W`, kind: "nominal", source: "Depends on package, board copper and airflow; use the datasheet value." },
  { id: "silicon.litho", label: "Scalar thin-mask lithography, no defocus or resist blur", kind: "nominal", source: "Kirchhoff mask; real OPC needs rigorous 3D mask and calibrated resist models." },
];
