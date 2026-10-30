// Sequential optics on the z axis, millimetres. Real rays are traced in 3D
// through spherical and flat surfaces with vector Snell refraction and the
// unpolarized Fresnel transmittance at each surface's actual incidence angle.
// Paraxial (y, u) tracing gives focal lengths. Assumptions are listed below and
// in assumptions(): glass indices from published Sellmeier fits, air n = 1,
// no bulk absorption, no ghost reflections, a uniform collimated beam.
import { clone } from "../studio-model.js";
export const GLASSES = {
  // Sellmeier B1..B3, C1..C3 (µm²). Schott N-BK7 and N-SF11 catalogue fits;
  // fused silica and CaF2 from Malitson (1965, 1963).
  "N-BK7": { title: "N-BK7 crown", B: [1.03961212, .231792344, 1.01046945], C: [.00600069867, .0200179144, 103.560653], density: 2.51 },
  "N-SF11": { title: "N-SF11 dense flint", B: [1.73759695, .313747346, 1.89878101], C: [.013188707, .0623068142, 155.23629], density: 4.74 },
  "fused-silica": { title: "Fused silica", B: [.6961663, .4079426, .8974794], C: [.0684043 ** 2, .1162414 ** 2, 9.896161 ** 2], density: 2.20 },
  "CaF2": { title: "Calcium fluoride", B: [.5675888, .4710914, 3.8484723], C: [.050263605 ** 2, .1003909 ** 2, 34.649040 ** 2], density: 3.18 },
  "N-SK16": { title: "N-SK16 dense crown", B: [1.34317774, .241144399, .994317969], C: [.00704687339, .0229005, 92.7508526], density: 3.58 },
  "F2": { title: "F2 flint", B: [1.34533359, .209073176, .937357162], C: [.00997743871, .0470450767, 111.886764], density: 3.61 },
};
export function refractiveIndex(glass, nm) {
  const g = GLASSES[glass]; if (!g) throw new Error(`Unknown glass ${glass}.`);
  const l2 = (nm / 1000) ** 2;
  return Math.sqrt(1 + g.B.reduce((s, b, i) => s + b * l2 / (l2 - g.C[i]), 0));
}
export const KINDS = ["lens", "window", "stop"];
// The starter lens: the Cooke triplet of the Optiland gallery (MIT),
// after Kingslake: EPD 10 mm, fields 0/14/20°, f ≈ 50 mm. Surface 4 is the
// stop; it is modelled as a plane at that vertex. Clear diameters are ours.
export const COOKE_SOURCE = "Optiland gallery, Cooke triplet (MIT licence), after R. Kingslake";
export function opticsExample() {
  const z1 = 10, z2 = z1 + 3.25896 + 6.00755, zs = z2 + .99997, z3 = zs + 4.75041;
  return {
    wavelength: 550, source: { power: 0, diameter: 10 }, fields: [0, 14, 20], pixel: 3.45, requiredGsd: 0,
    elements: [
      { id: "L1", name: "Front crown", kind: "lens", z: z1, r1: 22.01359, r2: -435.76044, thickness: 3.25896, diameter: 17, glass: "N-SK16" },
      { id: "L2", name: "Flint negative", kind: "lens", z: z2, r1: -22.21328, r2: 20.29192, thickness: .99997, diameter: 12, glass: "F2" },
      // The stop sits on surface 4; as a plane it goes just past that surface's edge sag.
      { id: "S1", name: "Aperture stop", kind: "stop", z: Math.round((zs + sag(20.29192, 6) + .01) * 1e5) / 1e5, diameter: 7.8 },
      { id: "L3", name: "Rear crown", kind: "lens", z: z3, r1: 79.68360, r2: -18.39533, thickness: 2.95208, diameter: 16, glass: "N-SK16" },
    ],
    target: { name: "Image sensor", z: Math.round((z3 + 2.95208 + 42.20778) * 1e5) / 1e5, diameter: 36.4 }, rays: 400,
  };
}
const finite = (n, lo, hi) => Number.isFinite(n) && n >= lo && n <= hi;
export function validateOptics(d) {
  if (!d || !finite(d.wavelength, 300, 2500) || !d.source || !finite(d.source.power, 0, 100) || !finite(d.source.diameter, .1, 100)
    || !d.target || !finite(d.target.z, -1000, 5000) || !finite(d.target.diameter, .01, 500) || typeof d.target.name !== "string" || d.target.name.length > 60
    || !Number.isInteger(d.rays) || d.rays < 16 || d.rays > 4000 || !Array.isArray(d.fields ?? [0]) || (d.fields ?? [0]).length > 5 || !(d.fields ?? [0]).every(f => finite(f, 0, 45))
    || !finite(d.pixel ?? 3.45, .5, 50) || !(d.requiredGsd === 0 || finite(d.requiredGsd ?? 1, .001, 1000)) || !Array.isArray(d.elements) || d.elements.length > 24 || new Set(d.elements.map(e => e?.id)).size !== d.elements.length)
    throw new Error("Invalid optical system settings.");
  for (const e of d.elements) {
    if (!e || typeof e.id !== "string" || !/^[a-zA-Z][a-zA-Z0-9_-]{0,23}$/.test(e.id) || typeof e.name !== "string" || !e.name.trim() || e.name.length > 60 || !KINDS.includes(e.kind)
      || !finite(e.z, -1000, 5000) || !finite(e.diameter, .1, 500)) throw new Error(`Invalid element ${e?.id ?? ""}.`);
    if (e.kind !== "stop" && (!finite(e.thickness, .05, 200) || !GLASSES[e.glass])) throw new Error(`${e.id} needs a thickness and a catalogue glass.`);
    if (e.kind === "lens" && (![e.r1, e.r2].every(r => finite(r, -1e5, 1e5)) || [e.r1, e.r2].some(r => r !== 0 && Math.abs(r) < e.diameter / 2)))
      throw new Error(`${e.id}: each radius must be 0 (flat) or at least the semi-diameter.`);
  }
  return {
    wavelength: d.wavelength, source: { power: d.source.power, diameter: d.source.diameter }, rays: d.rays,
    fields: [...new Set(d.fields ?? [0])].sort((a, b) => a - b), pixel: d.pixel ?? 3.45, requiredGsd: d.requiredGsd ?? 1,
    target: { name: d.target.name, z: d.target.z, diameter: d.target.diameter },
    elements: d.elements.map(e => e.kind === "stop" ? { id: e.id, name: e.name, kind: e.kind, z: e.z, diameter: e.diameter }
      : { id: e.id, name: e.name, kind: e.kind, z: e.z, thickness: e.thickness, diameter: e.diameter, glass: e.glass, ...(e.kind === "lens" ? { r1: e.r1, r2: e.r2 } : {}) }),
  };
}
// Ordered refracting/limiting surfaces: { z vertex, R (0 flat), aperture, n after }.
export function surfaces(d) {
  const out = [];
  for (const e of [...d.elements].sort((a, b) => a.z - b.z)) {
    if (e.kind === "stop") { out.push({ id: e.id, z: e.z, R: 0, a: e.diameter / 2, n: 1, stop: true }); continue; }
    const n = refractiveIndex(e.glass, d.wavelength);
    out.push({ id: e.id, z: e.z, R: e.kind === "lens" ? e.r1 : 0, a: e.diameter / 2, n });
    out.push({ id: e.id, z: e.z + e.thickness, R: e.kind === "lens" ? e.r2 : 0, a: e.diameter / 2, n: 1 });
  }
  return out;
}
export const sag = (R, r) => R === 0 ? 0 : R - Math.sign(R) * Math.sqrt(Math.max(0, R * R - r * r));
export function edgeThickness(e) { return e.kind !== "lens" ? e.thickness : e.thickness + sag(e.r2, e.diameter / 2) - sag(e.r1, e.diameter / 2); }
// Paraxial marginal ray from infinity: system focal length and back focus.
export function paraxial(d) {
  const ss = surfaces(d).filter(s => !s.stop);
  if (!ss.length) return { efl: Infinity, bfl: Infinity, focus: Infinity };
  let y = 1, u = 0, n = 1;
  for (let i = 0; i < ss.length; i++) {
    const s = ss[i], c = s.R ? 1 / s.R : 0;
    u = (n * u - y * c * (s.n - n)) / s.n; n = s.n;
    if (i + 1 < ss.length) y += u * (ss[i + 1].z - s.z);
  }
  const efl = -1 / u, bfl = -y / u;
  return { efl, bfl, focus: ss.at(-1).z + bfl };
}
// Thick-lens lensmaker equation in air, the same expression as cad::kits::optic::Lens::focal_length.
export function lensmaker(e, nm) {
  const n = refractiveIndex(e.glass, nm), c1 = e.r1 ? 1 / e.r1 : 0, c2 = e.r2 ? 1 / e.r2 : 0;
  return 1 / ((n - 1) * (c1 - c2 + (n - 1) * e.thickness * c1 * c2 / n));
}
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
function intersect(p, v, s) {
  if (s.R === 0) { if (Math.abs(v[2]) < 1e-15) return null; const t = (s.z - p[2]) / v[2]; return t > 1e-12 ? t : null; }
  const c = [0, 0, s.z + s.R], o = [p[0] - c[0], p[1] - c[1], p[2] - c[2]], b = dot(o, v), q = dot(o, o) - s.R * s.R, disc = b * b - q;
  if (disc < 0) return null;
  // The vertex-side cap: the hit lies on the hemisphere that contains the vertex.
  for (const t of [-b - Math.sqrt(disc), -b + Math.sqrt(disc)].sort((a, b) => a - b)) {
    if (t <= 1e-12) continue;
    if (Math.sign(p[2] + t * v[2] - c[2]) === Math.sign(-s.R)) return t;
  }
  return null;
}
export function fresnel(n1, n2, cosi, cost) {
  const rs = (n1 * cosi - n2 * cost) / (n1 * cosi + n2 * cost), rp = (n2 * cosi - n1 * cost) / (n2 * cosi + n1 * cost);
  return 1 - (rs * rs + rp * rp) / 2;
}
// One ray through the system. Returns its path, final state and why it stopped.
export function traceRay(ss, p0, v0, weight = 1) {
  let p = [...p0], v = [...v0], n = 1, w = weight;
  const path = [[...p]];
  for (const s of ss) {
    const t = intersect(p, v, s);
    if (t === null) return { path, p, v, w, lost: `misses ${s.id}` };
    p = [p[0] + t * v[0], p[1] + t * v[1], p[2] + t * v[2]]; path.push([...p]);
    if (Math.hypot(p[0], p[1]) > s.a + 1e-9) return { path, p, v, w: 0, lost: `clipped by ${s.id}`, blocked: s.id };
    if (s.stop || s.n === n) continue;
    let N = s.R === 0 ? [0, 0, -1] : [p[0] / s.R, p[1] / s.R, (p[2] - (s.z + s.R)) / s.R];
    if (dot(N, v) > 0) N = N.map(x => -x);
    const eta = n / s.n, cosi = -dot(N, v), k = 1 - eta * eta * (1 - cosi * cosi);
    if (k < 0) return { path, p, v, w: 0, lost: `total internal reflection at ${s.id}`, blocked: s.id };
    const cost = Math.sqrt(k);
    w *= fresnel(n, s.n, cosi, cost);
    v = v.map((x, i) => eta * x + (eta * cosi - cost) * N[i]); const len = Math.hypot(...v); v = v.map(x => x / len);
    n = s.n;
  }
  return { path, p, v, w };
}
// Equal-area (Vogel spiral) samples of a uniform circular beam.
export function pupil(diameter, count) {
  const golden = Math.PI * (3 - Math.sqrt(5)), r = diameter / 2;
  return Array.from({ length: count }, (_, i) => { const rr = r * Math.sqrt((i + .5) / count), t = i * golden; return [rr * Math.cos(t), rr * Math.sin(t)]; });
}
export function chiefOffset(d, deg) {
  const ss = surfaces(d), stop = ss.find(s => s.stop), start = Math.min(0, ...ss.map(s => s.z - 1)) - 5, a = deg * Math.PI / 180, v = [0, Math.sin(a), Math.cos(a)];
  if (!deg || !stop) return 0;
  const upto = ss.slice(0, ss.indexOf(stop)), at = y => { const r = traceRay(upto, [0, y, start], v); if (r.lost) return NaN; const t = (stop.z - r.p[2]) / r.v[2]; return r.p[1] + t * r.v[1]; };
  // Scan for a bracket where the chief height at the stop changes sign
  // (vignetted heights are skipped), then bisect inside it.
  let prev = null;
  for (let y = -40; y <= 40; y += .2) {
    const h = at(y); if (!Number.isFinite(h)) { prev = null; continue; }
    if (prev && Math.sign(h) !== Math.sign(prev[1])) {
      let lo = prev[0], hi = y, hlo = prev[1];
      for (let i = 0; i < 60; i++) { const m = (lo + hi) / 2, hm = at(m); if (!Number.isFinite(hm)) break; if (Math.sign(hm) === Math.sign(hlo)) { lo = m; hlo = hm; } else hi = m; }
      return (lo + hi) / 2;
    }
    prev = [y, h];
  }
  return 0;
}
export function trace(d, targetZ = d.target.z, field = 0) {
  const ss = surfaces(d), start = Math.min(0, ...ss.map(s => s.z - 1)) - 5, rays = [], a = field * Math.PI / 180, dir = [0, Math.sin(a), Math.cos(a)], y0 = chiefOffset(d, field);
  let onTarget = 0, transmitted = 0, sx = 0, sy = 0, sw = 0;
  const hits = [];
  for (const [x, y] of pupil(d.source.diameter, d.rays)) {
    const r = traceRay(ss, [x, y + y0, start], dir, 1 / d.rays);
    let landed = null;
    if (!r.lost && r.v[2] > 0) {
      const t = (targetZ - r.p[2]) / r.v[2];
      landed = [r.p[0] + t * r.v[0], r.p[1] + t * r.v[1], targetZ];
      transmitted += r.w;
      if (Math.hypot(landed[0], landed[1]) <= d.target.diameter / 2) onTarget += r.w;
      hits.push([landed[0], landed[1], r.w]); sx += landed[0] * r.w; sy += landed[1] * r.w; sw += r.w;
    }
    rays.push({ path: landed ? [...r.path, landed] : r.path, lost: r.lost, w: r.w });
  }
  const cx = sw ? sx / sw : 0, cy = sw ? sy / sw : 0;
  const rms = sw ? Math.sqrt(hits.reduce((s, [x, y, w]) => s + w * ((x - cx) ** 2 + (y - cy) ** 2), 0) / sw) : NaN;
  return { rays, hits, onTarget, transmitted, rms, centroid: [cx, cy] };
}
// Golden-section search for the plane of least RMS spot radius past the last surface.
export function bestFocus(d) {
  const ss = surfaces(d), last = ss.length ? ss.at(-1).z : 0, p = paraxial(d);
  let lo = last + .01, hi = Number.isFinite(p.focus) && p.focus > last ? p.focus + Math.max(5, (p.focus - last)) : last + 500;
  const f = z => trace({ ...d, rays: Math.min(d.rays, 300) }, z).rms, g = (Math.sqrt(5) - 1) / 2;
  let a = hi - g * (hi - lo), b = lo + g * (hi - lo), fa = f(a), fb = f(b);
  for (let i = 0; i < 60 && hi - lo > 1e-4; i++) {
    if (fa < fb) { hi = b; b = a; fb = fa; a = hi - g * (hi - lo); fa = f(a); } else { lo = a; a = b; fa = fb; b = lo + g * (hi - lo); fb = f(b); }
  }
  return (lo + hi) / 2;
}
export function geometryIssues(d) {
  const issues = [], sorted = [...d.elements].sort((a, b) => a.z - b.z);
  for (const e of sorted) if (e.kind === "lens" && edgeThickness(e) < 0.5) issues.push({ code: "edge", object: e.id, message: `${e.name} edge thickness ${edgeThickness(e).toFixed(2)} mm is below 0.5 mm.` });
  for (let i = 1; i < sorted.length; i++) {
    const a = sorted[i - 1], b = sorted[i], aEnd = a.z + (a.thickness || 0), r = Math.min(a.diameter, b.diameter) / 2;
    const gap = b.z + (b.kind === "lens" ? Math.min(0, sag(b.r1, r)) : 0) - (aEnd + (a.kind === "lens" ? Math.max(0, sag(a.r2, r)) : 0));
    if (gap < 0) issues.push({ code: "overlap", object: b.id, message: `${b.name} intersects ${a.name}.` });
  }
  const last = sorted.at(-1);
  if (last && d.target.z <= last.z + (last.thickness || 0)) issues.push({ code: "target", object: "target", message: "The target plane must lie after the last element." });
  return issues;
}
// Glass volume by quadrature of the axial thickness over the clear aperture.
export function elementVolume(e) {
  if (e.kind === "stop") return 0;
  const a = e.diameter / 2, n = 400; let v = 0;
  for (let i = 0; i < n; i++) { const r = (i + .5) * a / n, t = e.kind === "lens" ? e.thickness + sag(e.r2, r) - sag(e.r1, r) : e.thickness; v += 2 * Math.PI * r * Math.max(0, t) * a / n; }
  return v;
}
export const elementMass = e => e.kind === "stop" ? 0 : elementVolume(e) * GLASSES[e.glass].density / 1000; // g
export function opticsResults(d, altitude = 60) {
  const t = trace(d), p = paraxial(d), issues = geometryIssues(d), area = Math.PI * (d.target.diameter / 2) ** 2;
  const fields = d.fields.map(f => { const r = f ? trace({ ...d, rays: Math.min(d.rays, 300) }, d.target.z, f) : t; return { field: f, rms: r.rms, centroid: r.centroid, transmitted: r.transmitted, rays: r.rays, hits: r.hits }; });
  const blocked = t.rays.filter(r => r.lost).length, mass = d.elements.reduce((s, e) => s + elementMass(e), 0);
  const gsd = d.pixel * 1e-6 * altitude / (p.efl * 1e-3) * 100, fnumber = p.efl / d.source.diameter;
  return { trace: t, paraxial: p, issues, fields, mass, gsd, fnumber, powerOnTarget: d.source.power * t.onTarget, transmitted: t.transmitted, irradiance: d.source.power * t.onTarget / area, blocked };
}
export function opticsAssumptions(d) {
  return [
    ...[...new Set(d.elements.filter(e => e.glass).map(e => e.glass))].map(g => ({ id: `optics.glass.${g}`, label: `${GLASSES[g].title} index ${refractiveIndex(g, d.wavelength).toFixed(5)} at ${d.wavelength} nm`, kind: "nominal", source: "Published Sellmeier fit; melt data varies by lot (±0.0005 typical)." })),
    { id: "optics.air", label: "Air index taken as 1.00000", kind: "nominal", source: "Real air is ≈1.00027 at STP; shifts focus by parts in 10⁴." },
    { id: "optics.beam", label: `Object at infinity, ${d.source.diameter} mm entrance beam, fields ${d.fields.join("/")}°`, kind: "prescribed", source: "Each field's chief ray is aimed through the stop centre; one wavelength at a time; the stop is a plane just past its surface." },
    { id: "optics.losses", label: "Fresnel losses only: no bulk absorption, scattering or ghosts", kind: "nominal", source: "Uncoated surfaces; coatings would raise transmission." },
  ];
}
export const opticsMesh = (e) => e.kind === "stop" ? null : ({
  name: e.name, cell: Math.max(.15, e.diameter / 60),
  shape: { op: "translate", by: [0, 0, e.z], shape: e.kind === "lens"
    ? { op: "kit", kit: "lens", args: { front: { radius: e.r1 || 1e12 }, back: { radius: e.r2 || 1e12 }, thickness: e.thickness, semi_diameter: e.diameter / 2, points: 64 } }
    : { op: "cylinder_z", center: [0, 0, e.thickness / 2], radius: e.diameter / 2, height: e.thickness } },
});
export const copy = clone;
