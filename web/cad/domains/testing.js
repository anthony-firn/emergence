// A test plan, kept honest about what kind of evidence each result is:
//   simulation   a model's prediction against a limit
//   correlation  a model against a measured dataset it was not fitted to
//   software     real software or a real tool run on the design and recorded: flight
//                code against a simulated aircraft, a CAD tool's design-rule check
//   physical     a measurement on real hardware, recorded by a person
// A recorded test stays "not run" until someone records it, so the plan always
// shows what has been predicted, what has been checked against measurements, and
// what has actually been tested. A result that names the documents it covers goes
// stale when one of them changes: it describes a design that no longer exists.
export const KINDS = ["simulation", "correlation", "software", "physical"];
export const testingExample = () => ({ datasets: [], traces: [], plan: [] });
const canonical = v => Array.isArray(v) ? v.map(canonical) : v && typeof v === "object" ? Object.fromEntries(Object.keys(v).sort().map(k => [k, canonical(v[k])])) : v;
export function digest(docs, covers) { let h = 0x811c9dc5; for (const c of JSON.stringify(covers.map(k => canonical(docs[k] ?? null)))) { h ^= c.charCodeAt(0); h = Math.imul(h, 16777619); } return (h >>> 0).toString(16).padStart(8, "0"); }
const ID = /^[a-zA-Z][a-zA-Z0-9_-]{0,31}$/, PORT = /^[a-z]+\.[a-zA-Z0-9_]+$/, finite = v => Number.isFinite(v);
const text = (v, n, what) => { if (v !== undefined && (typeof v !== "string" || v.length > n)) throw new Error(`${what} is up to ${n} characters.`); return v; };
export function validateTesting(d) {
  if (!d || !Array.isArray(d.plan) || d.plan.length > 48 || !Array.isArray(d.datasets) || d.datasets.length > 16) throw new Error("Invalid test plan.");
  const datasets = d.datasets.map(s => {
    if (!s || !ID.test(s.id ?? "") || typeof s.title !== "string" || !s.title.trim() || s.title.length > 120 || !Array.isArray(s.rows) || s.rows.length < 2 || s.rows.length > 200 || !s.rows.every(r => finite(r?.rpm) && r.rpm > 0 && finite(r.grams) && r.grams > 0)) throw new Error("A dataset has an id, a title and rows of rpm and grams.");
    if (s.source !== undefined && !/^https:\/\/\S{4,300}$/.test(s.source)) throw new Error("A dataset's source is an https link.");
    return { id: s.id, title: s.title, ...(s.source && { source: s.source }), rows: s.rows.map(r => ({ rpm: r.rpm, grams: r.grams, ...(finite(r.torque) && { torque: r.torque }) })) };
  });
  const traces = (d.traces ?? []).map(s => {
    if (!s || !ID.test(s.id ?? "") || typeof s.title !== "string" || !s.title.trim() || s.title.length > 120 || !Array.isArray(s.rows) || s.rows.length < 2 || s.rows.length > 400 || !s.rows.every(r => Array.isArray(r) && r.length === 3 && r.every(finite))) throw new Error("A trace has an id, a title, and rows of time, measured and model.");
    if (s.source !== undefined && !/^https:\/\/\S{4,300}$/.test(s.source)) throw new Error("A trace's source is an https link.");
    return { id: s.id, title: s.title, ...(text(s.unit, 12, "A unit") && { unit: s.unit }), ...(text(s.note, 300, "A note") && { note: s.note }), ...(s.source && { source: s.source }), rows: s.rows.map(r => [...r]) };
  });
  const plan = d.plan.map(t => {
    if (!t || !ID.test(t.id ?? "") || typeof t.title !== "string" || !t.title.trim() || t.title.length > 120 || !KINDS.includes(t.kind)) throw new Error("A test has an id, a title and a kind.");
    if (t.port !== undefined && !PORT.test(t.port)) throw new Error(`${t.id}: a port is domain.output.`);
    for (const k of ["min", "max"]) if (t[k] !== undefined && !finite(t[k])) throw new Error(`${t.id}: ${k} is a number.`);
    if (t.kind === "correlation" && !datasets.some(s => s.id === t.dataset)) throw new Error(`${t.id}: choose a dataset.`);
    const r = t.result; if (r !== undefined && (!r || !finite(r.value))) throw new Error(`${t.id}: a result has a value.`);
    if (r?.digest !== undefined && !/^[0-9a-f]{8}$/.test(r.digest)) throw new Error(`${t.id}: invalid result digest.`);
    if (t.covers !== undefined && (!Array.isArray(t.covers) || t.covers.length > 12 || !t.covers.every(k => typeof k === "string" && /^[a-z]{1,16}$/.test(k)))) throw new Error(`${t.id}: covers is a list of documents.`);
    return { id: t.id, title: t.title, kind: t.kind, ...(text(t.method, 400, "A method") && { method: t.method }), ...(t.port && { port: t.port }), ...(t.dataset && { dataset: t.dataset }), ...(t.covers?.length && { covers: [...t.covers] }), ...(t.min !== undefined && { min: t.min }), ...(t.max !== undefined && { max: t.max }), ...(text(t.unit, 12, "A unit") && { unit: t.unit }),
      ...(r && { result: { value: r.value, ...(text(r.date, 24, "A date") && { date: r.date }), ...(text(r.by, 60, "A name") && { by: r.by }), ...(text(r.note, 300, "A note") && { note: r.note }), ...(r.digest && { digest: r.digest }) } }) };
  });
  if (new Set(plan.map(t => t.id)).size !== plan.length || new Set(datasets.map(s => s.id)).size !== datasets.length || new Set(traces.map(s => s.id)).size !== traces.length) throw new Error("Test, dataset and trace ids are unique.");
  return { datasets, traces, plan };
}
// Rotor similarity fitted to the design's bench table, then asked to predict
// measurements it never saw: T = c·n², and Q = k·n² where torque was measured.
export function correlate(set, table) {
  const fit = (key) => { const rows = table.filter(r => r[key] > 0), num = rows.reduce((s, r) => s + r[key] * r.rpm ** 2, 0), den = rows.reduce((s, r) => s + r.rpm ** 4, 0); return num / den; };
  const c = fit("grams"), k = fit("torque");
  const rows = set.rows.map(r => ({ ...r, predicted: c * r.rpm ** 2, error: 100 * (c * r.rpm ** 2 / r.grams - 1), ...(r.torque !== undefined && { torqueError: 100 * (k * r.rpm ** 2 / r.torque - 1) }) }));
  const rms = key => { const v = rows.filter(r => r[key] !== undefined).map(r => r[key]); return v.length ? Math.sqrt(v.reduce((s, x) => s + x * x, 0) / v.length) : NaN; };
  return { c, rows, rms: rms("error"), worst: Math.max(...rows.map(r => Math.abs(r.error))), torqueRms: rms("torqueError") };
}
const within = (t, v) => (t.min === undefined || v >= t.min) && (t.max === undefined || v <= t.max);
export function evaluatePlan(d, { outputs, table, docs }) {
  const items = d.plan.map(t => {
    const base = { id: t.id, title: t.title, kind: t.kind, method: t.method, min: t.min, max: t.max };
    if (t.kind === "physical" || t.kind === "software") {
      if (!t.result) return { ...base, status: "not run", unit: t.unit };
      const stale = t.result.digest && t.covers?.length && docs && digest(docs, t.covers) !== t.result.digest;
      return { ...base, status: stale ? "stale" : within(t, t.result.value) ? "passed" : "failed", value: t.result.value, unit: t.unit, result: t.result, ...(stale && { detail: `Recorded before ${t.covers.join(", ")} changed. Run it again.` }) };
    }
    if (t.kind === "correlation") {
      const set = d.datasets.find(s => s.id === t.dataset);
      if (!table) return { ...base, status: "not run", detail: "The project has no measured propulsion table to correlate." };
      const c = correlate(set, table); return { ...base, status: within(t, c.rms) ? "passed" : "failed", value: c.rms, unit: "%", correlation: c, dataset: set };
    }
    const o = outputs[t.port]; if (!o || !Number.isFinite(o.value)) return { ...base, status: "not run", detail: `${t.port} has no value in this project.` };
    return { ...base, status: within(t, o.value) ? "passed" : "failed", value: o.value, unit: o.unit };
  });
  const count = (k, s) => items.filter(i => i.kind === k && i.status === s).length, total = k => items.filter(i => i.kind === k).length;
  return { items, summary: Object.fromEntries(KINDS.map(k => [k, { total: total(k), passed: count(k, "passed"), failed: count(k, "failed"), stale: count(k, "stale"), open: count(k, "not run") }])) };
}
