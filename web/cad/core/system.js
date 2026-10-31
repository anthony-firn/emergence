// Cross-domain links. Each domain declares typed input and output ports; a link
// carries one output into one input with unit conversion. The solver runs the
// domains in dependency order and iterates coupled loops (heat raises the room
// temperature, which raises the inlet air temperature, which raises the heat)
// until the linked values stop changing, and reports when they do not.
import { canonical } from "../repository.js";
import { compatible, convert, label } from "./units.js";

export const SYSTEM_DEFAULT = () => ({ links: [], confirmed: {} });
const portDomain = id => id.split(".")[0];
const portName = id => id.slice(id.indexOf(".") + 1);

export function validateSystem(doc, known) {
  if (!doc || typeof doc !== "object" || !Array.isArray(doc.links) || doc.links.length > 64 || !doc.confirmed || typeof doc.confirmed !== "object" || Array.isArray(doc.confirmed))
    throw new Error("Invalid system links document.");
  const ids = new Set();
  const links = doc.links.map(l => {
    if (!l || typeof l.id !== "string" || !/^[a-zA-Z0-9_-]{1,48}$/.test(l.id) || ids.has(l.id)) throw new Error("Each link needs a unique id.");
    ids.add(l.id);
    if (typeof l.from !== "string" || typeof l.to !== "string" || !/^[a-z]+\.[a-zA-Z0-9_.]+$/.test(l.from) || !/^[a-z]+\.[a-zA-Z0-9_.]+$/.test(l.to)) throw new Error(`Link ${l.id} has an invalid port.`);
    if (known && (!known.outputs.has(l.from) || !known.inputs.has(l.to))) throw new Error(`Link ${l.id} refers to an unknown port.`);
    return { id: l.id, from: l.from, to: l.to, enabled: l.enabled !== false, ...(typeof l.note === "string" && l.note && { note: l.note.slice(0, 160) }) };
  });
  const confirmed = Object.fromEntries(Object.entries(doc.confirmed).filter(([k, v]) => /^[a-z]+\.[a-zA-Z0-9_.-]+$/.test(k) && v === true).slice(0, 256));
  return { links, confirmed };
}

export class System {
  // `enabled` names the domains in the project; links to the others are dormant.
  constructor(adapters, docs, enabled = () => true) { this.adapters = adapters; this.docs = docs; this.enabled = enabled; this.memo = new Map(); this.warm = {}; }
  // A what-if solve on modified documents that shares the evaluation cache
  // but not the warm start, so studies never disturb the live solution.
  fork(docs) { const s = new System(this.adapters, () => docs, this.enabled); s.memo = this.memo; s.warm = { ...this.warm }; s.applied = this.applied; return s; }
  ports() {
    const docs = this.docs(), inputs = new Map(), outputs = new Map();
    for (const [id, a] of this.adapters) {
      if (!this.enabled(id)) continue;
      const p = a.ports(docs[id]);
      for (const i of p.inputs) inputs.set(`${id}.${i.id}`, { ...i, domain: id, port: `${id}.${i.id}` });
      for (const o of p.outputs) outputs.set(`${id}.${o.id}`, { ...o, domain: id, port: `${id}.${o.id}` });
    }
    return { inputs, outputs };
  }
  // A link is usable when both ports exist and their dimensions agree.
  check(link, ports = this.ports()) {
    const from = ports.outputs.get(link.from), to = ports.inputs.get(link.to);
    if (!from) return `Output ${link.from} does not exist.`;
    if (!to) return `Input ${link.to} does not exist.`;
    if (!compatible(from.unit, to.unit)) return `${from.label} (${label(from.unit) || "dimensionless"}) cannot drive ${to.label} (${label(to.unit) || "dimensionless"}).`;
    if (portDomain(link.from) === portDomain(link.to) && portName(link.from) === portName(link.to)) return "A port cannot drive itself.";
    return null;
  }
  async evaluate(id, doc, overrides) {
    const key = id + canonical({ doc, overrides });
    if (this.memo.has(key)) return this.memo.get(key);
    const r = await this.adapters.get(id).evaluate(doc, overrides);
    if (this.memo.size > 400) this.memo.clear();
    this.memo.set(key, r); return r;
  }
  order(links) {
    const ids = [...this.adapters.keys()].filter(this.enabled), edges = new Map(ids.map(id => [id, new Set()])), seen = new Set(), done = [], stack = new Set();
    for (const l of links) if (portDomain(l.from) !== portDomain(l.to) && edges.has(portDomain(l.from))) edges.get(portDomain(l.from)).add(portDomain(l.to));
    let cyclic = false;
    const visit = id => { if (stack.has(id)) { cyclic = true; return; } if (seen.has(id)) return; seen.add(id); stack.add(id); for (const n of edges.get(id) || []) visit(n); stack.delete(id); done.push(id); };
    ids.forEach(visit);
    return { order: done.reverse(), cyclic };
  }
  async solve(systemDoc, { maxIterations = 16, tolerance = 1e-4, changed = [] } = {}) {
    const started = performance.now(), docs = this.docs(), ports = this.ports(), issues = [];
    const on = l => this.enabled(portDomain(l.from)) && this.enabled(portDomain(l.to));
    const links = systemDoc.links.filter(l => l.enabled && on(l)), dormant = systemDoc.links.filter(l => !on(l)).map(l => l.id);
    const usable = links.filter(l => { const problem = this.check(l, ports); if (problem) issues.push({ link: l.id, kind: "invalid", message: problem }); return !problem; });
    let { order, cyclic } = this.order(usable);
    // In a loop, start from the edited domain so its change travels the
    // cycle once before the pass that confirms it.
    const first = cyclic && order.findIndex(id => changed.includes(id));
    if (first > 0) order = [...order.slice(first), ...order.slice(0, first)];
    const outputs = {}, results = {}, applied = {};
    // Warm start: links whose source runs later in the order begin from the
    // previous solution instead of the unlinked value, so an edit near the
    // last fixed point usually settles in two passes.
    let previous = Object.fromEntries(Object.entries(this.applied || {}).map(([k, v]) => [k, v.value])), residual = Infinity, iterations = 0, omega = 1, last = Infinity;
    const history = [];
    for (iterations = 1; iterations <= maxIterations; iterations++) {
      const values = {}, passIssues = [];
      for (const id of order) {
        const overrides = {};
        for (const [port, input] of ports.inputs) {
          if (input.domain !== id) continue;
          const source = l => outputs[l.from] ?? (iterations === 1 ? this.warm[l.from] : undefined);
          const incoming = usable.filter(l => l.to === port && Number.isFinite(source(l)?.value));
          if (!incoming.length) continue;
          if (incoming.length > 1 && input.aggregate !== "sum") passIssues.push({ link: incoming[1].id, kind: "multiple", message: `${input.label} accepts one link; using ${incoming[0].from}.` });
          const used = input.aggregate === "sum" ? incoming : incoming.slice(0, 1);
          let v = used.reduce((s, l) => s + convert(source(l).value, source(l).unit, input.unit), 0);
          if (previous[port] !== undefined && omega < 1) v = previous[port] + omega * (v - previous[port]);
          if (input.min !== undefined && v < input.min || input.max !== undefined && v > input.max) {
            passIssues.push({ link: used[0].id, kind: "range", port, message: `${input.label} would be ${round(v)} ${label(input.unit)}, outside its validated ${input.min}–${input.max} ${label(input.unit)}. The ${id} design keeps its own value.` });
            continue;
          }
          overrides[portName(port)] = v; values[port] = { value: v, unit: input.unit, from: used.map(l => l.from), links: used.map(l => l.id) };
        }
        let r;
        try { r = await this.evaluate(id, docs[id], overrides); }
        catch (e) { r = { outputs: {}, checks: [{ id: `${id}-model`, text: `${this.adapters.get(id).title} model`, status: "unavailable", reason: e.message }], error: e.message }; }
        results[id] = r;
        for (const [k, o] of Object.entries(r.outputs || {})) outputs[`${id}.${k}`] = { ...o, domain: id };
      }
      residual = 0;
      for (const [port, v] of Object.entries(values)) {
        const before = previous[port]?.value ?? previous[port];
        residual = Math.max(residual, before === undefined ? Infinity : Math.abs(v.value - before) / Math.max(1e-9, Math.abs(v.value)));
      }
      if (!Object.keys(values).length) residual = 0;
      history.push(residual);
      Object.assign(applied, values);
      for (const k of Object.keys(applied)) if (!values[k]) delete applied[k];
      previous = Object.fromEntries(Object.entries(values).map(([k, v]) => [k, v.value]));
      if (iterations === maxIterations || residual <= tolerance || !cyclic && iterations >= 2) { issues.push(...passIssues); break; }
      if (iterations > 3 && residual > last) omega = .5;
      last = residual;
    }
    const converged = residual <= tolerance || !cyclic;
    this.warm = outputs; this.applied = applied;
    // Work that reads the settled system and feeds nothing back into it: flight
    // tests, the test plan. It runs once, after the coupled loop.
    for (const [id, a] of this.adapters) {
      if (!a.after || !this.enabled(id)) continue;
      try { const r = await a.after(docs[id], { outputs, results, docs }); results[id] = { ...results[id], ...r, outputs: { ...results[id]?.outputs, ...r.outputs }, checks: [...(results[id]?.checks ?? []), ...(r.checks ?? [])] }; for (const [k, o] of Object.entries(r.outputs ?? {})) outputs[`${id}.${k}`] = { ...o, domain: id }; }
      catch (e) { results[id] = { ...results[id], checks: [...(results[id]?.checks ?? []), { id: `${id}-after`, text: `${a.title} evaluation`, status: "unavailable", reason: e.message }] }; }
    }
    if (!converged) issues.push({ kind: "convergence", message: `Coupled links did not settle in ${maxIterations} passes (relative change ${residual.toExponential(1)}). Values shown are the last pass.` });
    const checks = Object.entries(results).flatMap(([domain, r]) => (r.checks || []).map(c => ({ ...c, domain })));
    // The project's own requirements, on any output. An interval that straddles
    // a limit is undecided, not met.
    for (const q of docs.project?.requirements ?? []) {
      const o = outputs[q.port], p = ports.outputs.get(q.port), v = o?.value, lo = o?.lo ?? v, hi = o?.hi ?? v, bound = [q.min !== undefined && `≥ ${q.min}`, q.max !== undefined && `≤ ${q.max}`].filter(Boolean).join(" and ");
      const inside = x => (q.min === undefined || x >= q.min) && (q.max === undefined || x <= q.max);
      checks.push({ id: q.id, domain: "project", text: q.label || `${p?.label || q.port} ${bound}${o ? ` ${label(o.unit)}` : ""}`.trim(), ...(Number.isFinite(v) ? { status: inside(lo) && inside(hi) ? "met" : inside(v) || inside(lo) || inside(hi) ? "undecided" : "unmet", value: { nominal: v }, unit: o.unit } : { status: "unavailable", reason: o ? `${q.port} has no value.` : `${q.port} is not in this project.` }) });
    }
    return { outputs, applied, results, checks, issues, iterations, residual, converged, cyclic, order, history, dormant, ms: performance.now() - started, ports };
  }
}
const round = v => Number(v.toPrecision(4));
