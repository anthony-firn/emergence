// Series air path: fan, grille, filter, ducts, bends and the heat-sink channels.
// Laminar friction in rectangular passages uses Shah & London (1978)'s exact
// fit fRe(α); turbulent friction uses Colebrook–White on the hydraulic
// diameter. Between Re 2300 and 4000 neither law applies, so the operating
// point is bounded between them and reported as an interval, as the Rust pipe
// kit refuses to invent a transition value. Fan curves are illustrative fits.
export const KINDS = ["duct", "fitting", "filter", "channels"];
export const FANS = {
  axial60: { title: "60 mm axial fan (illustrative)", shutoff: 30, free: .0085 },
  axial80: { title: "80 mm axial fan (illustrative)", shutoff: 45, free: .0125 },
  axial80hp: { title: "80 mm high-pressure axial (illustrative)", shutoff: 120, free: .0165 },
  axial92: { title: "92 mm axial fan (illustrative)", shutoff: 60, free: .021 },
  blower60: { title: "60 mm radial blower (illustrative)", shutoff: 300, free: .006 },
  // Rotor downwash captured by the duct inlet: a pressure source with no fan
  // curve. The pressure comes from the Flight domain's wake.
  wake: { title: "Rotor wake (linked from Flight)", shutoff: 30, free: 10 },
};
export function fluidsExample() {
  return {
    temperature: 25, pressure: 101325, fan: { preset: "axial60", shutoff: 30, free: .0085 }, heatLimit: 15,
    elements: [
      { id: "inlet", name: "Inlet grille", kind: "fitting", K: 1.5, width: 80, height: 80 },
      { id: "filter", name: "Dust filter", kind: "filter", coefficient: 25, width: 80, height: 80 },
      { id: "duct", name: "Supply duct", kind: "duct", width: 60, height: 30, length: 180, roughness: .0015 },
      { id: "bend", name: "90° bend", kind: "fitting", K: .9, width: 60, height: 30 },
      { id: "sink", name: "Heat-sink channels", kind: "channels", fins: 9, gap: 3.65, height: 20, length: 50, width: 40, Kin: .5, Kout: 1 },
    ],
  };
}
const finite = (n, lo, hi) => Number.isFinite(n) && n >= lo && n <= hi;
export function validateFluids(d) {
  if (!d || !finite(d.temperature, -40, 150) || !finite(d.pressure, 5e4, 2e5) || !d.fan || !finite(d.fan.shutoff, 1, 5000) || !finite(d.fan.free, 1e-5, 10)
    || !(d.fan.preset === "custom" || FANS[d.fan.preset]) || !finite(d.heatLimit, .1, 200) || !Array.isArray(d.elements) || d.elements.length > 24
    || new Set(d.elements.map(e => e?.id)).size !== d.elements.length) throw new Error("Invalid air-path settings.");
  const out = d.elements.map(e => {
    if (!e || typeof e.id !== "string" || !/^[a-zA-Z][a-zA-Z0-9_-]{0,23}$/.test(e.id) || typeof e.name !== "string" || !e.name.trim() || e.name.length > 60 || !KINDS.includes(e.kind)
      || !finite(e.width, .5, 2000) || !finite(e.height, .5, 2000)) throw new Error(`Invalid element ${e?.id ?? ""}.`);
    const base = { id: e.id, name: e.name, kind: e.kind, width: e.width, height: e.height };
    if (e.kind === "duct") { if (!finite(e.length, .1, 1e5) || !finite(e.roughness, 0, 10)) throw new Error(`${e.id}: length and roughness.`); return { ...base, length: e.length, roughness: e.roughness }; }
    if (e.kind === "fitting") { if (!finite(e.K, 0, 100)) throw new Error(`${e.id}: loss coefficient K.`); return { ...base, K: e.K }; }
    if (e.kind === "filter") { if (!finite(e.coefficient, 0, 1e4)) throw new Error(`${e.id}: filter coefficient.`); return { ...base, coefficient: e.coefficient }; }
    if (!Number.isInteger(e.fins) || e.fins < 2 || e.fins > 200 || !finite(e.gap, .05, 200) || !finite(e.length, .1, 5000) || !finite(e.Kin, 0, 10) || !finite(e.Kout, 0, 10))
      throw new Error(`${e.id}: channel count, gap, length and K values.`);
    return { ...base, fins: e.fins, gap: e.gap, length: e.length, Kin: e.Kin, Kout: e.Kout };
  });
  if (!finite(d.heatLoad ?? 0, 0, 1e6)) throw new Error("Invalid heat carried by the air.");
  return { temperature: d.temperature, pressure: d.pressure, heatLimit: d.heatLimit, heatLoad: d.heatLoad ?? 0, fan: { preset: d.fan.preset, shutoff: d.fan.shutoff, free: d.fan.free }, elements: out };
}
// Dry air: ideal gas density, Sutherland viscosity, constant cp.
export function air(celsius, pressure = 101325) {
  const T = celsius + 273.15;
  return { rho: pressure / (287.05 * T), mu: 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4), cp: 1006 };
}
// Darcy fRe for laminar flow in a rectangle of aspect ratio α = short/long;
// 96 for parallel plates, 56.91 for a square (Shah & London 1978).
export const shahLondon = a => 96 * (1 - 1.3553 * a + 1.9467 * a ** 2 - 1.7012 * a ** 3 + .9564 * a ** 4 - .2537 * a ** 5);
export function colebrook(re, relRough) {
  let x = -2 * Math.log10(relRough / 3.7 + 5.74 / re ** .9); // Swamee–Jain start
  for (let i = 0; i < 50; i++) { const next = -2 * Math.log10(relRough / 3.7 + 2.51 * x / re); if (Math.abs(next - x) < 1e-13) { x = next; break; } x = next; }
  return 1 / (x * x);
}
// Friction factor bounds: equal outside the transition band.
export function friction(re, alpha, relRough) {
  if (re <= 0) return { lo: 0, hi: 0, regime: "still" };
  const lam = shahLondon(alpha) / re;
  if (re < 2300) return { lo: lam, hi: lam, regime: "laminar" };
  const turb = colebrook(re, relRough);
  if (re > 4000) return { lo: turb, hi: turb, regime: "turbulent" };
  return { lo: Math.min(lam, turb), hi: Math.max(lam, turb), regime: "transition" };
}
const rect = (w, h) => ({ area: w * h, dh: 2 * w * h / (w + h), alpha: Math.min(w, h) / Math.max(w, h) });
export function sinkGeometry(e) {
  const channels = e.fins - 1, approach = e.width * e.height * 1e-6;
  return { channels, approach, channel: rect(e.gap * 1e-3, e.height * 1e-3) };
}
// Pressure drop of one element at flow Q (m³/s); returns { lo, hi, ... }.
export function elementDrop(e, Q, props) {
  const { rho, mu } = props, dyn = v => .5 * rho * v * v;
  if (e.kind === "fitting") { const v = Q / (e.width * e.height * 1e-6); return { lo: e.K * dyn(v), hi: e.K * dyn(v), velocity: v, regime: "fitting" }; }
  if (e.kind === "filter") { const v = Q / (e.width * e.height * 1e-6); return { lo: e.coefficient * v, hi: e.coefficient * v, velocity: v, regime: "filter" }; }
  if (e.kind === "duct") {
    const g = rect(e.width * 1e-3, e.height * 1e-3), v = Q / g.area, re = rho * v * g.dh / mu, f = friction(re, g.alpha, e.roughness * 1e-3 / g.dh), k = e.length * 1e-3 / g.dh * dyn(v);
    return { lo: f.lo * k, hi: f.hi * k, velocity: v, re, regime: f.regime };
  }
  const s = sinkGeometry(e), v = Q / s.channels / s.channel.area, re = rho * v * s.channel.dh / mu, f = friction(re, s.channel.alpha, 1e-6 / s.channel.dh);
  const k = e.length * 1e-3 / s.channel.dh * dyn(v), minor = (e.Kin + e.Kout) * dyn(v);
  return { lo: f.lo * k + minor, hi: f.hi * k + minor, velocity: v, re, regime: f.regime, approach: Q / s.approach };
}
export const fanPressure = (fan, Q) => fan.shutoff * (1 - (Q / fan.free) ** 2);
const system = (d, Q, props, side) => d.elements.reduce((s, e) => s + elementDrop(e, Q, props)[side], 0);
function root(fn, lo, hi) {
  if (fn(lo) < 0 || fn(hi) > 0) return null;
  for (let i = 0; i < 200 && hi - lo > 1e-14; i++) { const m = (lo + hi) / 2; fn(m) > 0 ? lo = m : hi = m; }
  return (lo + hi) / 2;
}
// Operating point: the fan curve meets the system curve. The low-friction
// bound gives the high flow; the high-friction bound gives the low flow.
export function operatingPoint(d, temperature = d.temperature) {
  const props = air(temperature, d.pressure), fan = d.fan;
  const high = root(Q => fanPressure(fan, Q) - system(d, Q, props, "lo"), 1e-12, fan.free);
  const low = root(Q => fanPressure(fan, Q) - system(d, Q, props, "hi"), 1e-12, fan.free);
  if (high === null || low === null) return { props, error: "The fan cannot move air through this path." };
  const Q = (low + high) / 2, elements = d.elements.map(e => ({ id: e.id, name: e.name, kind: e.kind, ...elementDrop(e, Q, props) }));
  return { props, Q, low, high, pressure: fanPressure(fan, Q), elements, total: { lo: system(d, low, props, "hi"), hi: system(d, high, props, "lo") } };
}
export function fluidsResults(d, overrides = {}) {
  const doc = structuredClone(d);
  const sink = doc.elements.find(e => e.kind === "channels");
  if (sink) for (const [k, key] of [["fin_count", "fins"], ["sink_gap", "gap"], ["sink_height", "height"], ["sink_length", "length"], ["sink_width", "width"]]) if (overrides[k] !== undefined) sink[key] = key === "fins" ? Math.round(overrides[k]) : overrides[k];
  const t = overrides.air_temperature ?? doc.temperature, op = operatingPoint(doc, t), heat = overrides.heat_load ?? doc.heatLoad;
  return { doc, op, sink, temperature: t, heat };
}
export function fluidsAssumptions(d) {
  return [
    d.fan.preset === "wake" ? { id: "fluids.wake", label: "Rotor wake as a constant pressure source at the duct inlet", kind: "nominal", source: "The captured share of the hover wake's dynamic pressure, linked from Flight; gusts, forward flight and rotor-to-rotor interaction are not modelled." }
      : { id: "fluids.fan", label: `Fan curve ${d.fan.shutoff} Pa shut-off, ${(d.fan.free * 1000).toFixed(1)} L/s free delivery`, kind: d.fan.preset === "custom" ? "user" : "illustrative", source: "Quadratic fit of a typical axial-fan curve. Use the manufacturer's measured curve." },
    { id: "fluids.losses", label: "Fitting K and filter coefficients", kind: "nominal", source: "Tabulated values (Crane TP-410 style); real fittings vary ±30%." },
    { id: "fluids.air", label: "Dry air: ideal gas, Sutherland viscosity, cp 1006 J/kg·K", kind: "nominal", source: "Humidity and compressibility neglected (Mach < 0.05)." },
    { id: "fluids.split", label: "Equal flow in every heat-sink channel", kind: "nominal", source: "Manifold maldistribution and bypass around the fins are not modelled." },
  ];
}
