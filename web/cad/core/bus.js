// Every edit is a command. The UI, the command line, the guided demo, scripts,
// the CLI and the MCP server all call the same registry, so anything a person
// can do an agent can do, and every change has one undo path and one log.
const finite = Number.isFinite;
const describeArg = a => ({ name: a.name, type: a.type, ...(a.unit !== undefined && { unit: a.unit }), ...(a.min !== undefined && { min: a.min }),
  ...(a.max !== undefined && { max: a.max }), ...(a.values && { values: a.values }), ...(a.integer && { integer: true }),
  ...(a.required && { required: true }), ...(a.default !== undefined && { default: a.default }), ...(a.help && { help: a.help }) });

export function coerce(arg, raw) {
  const where = `${arg.name}`;
  if (raw === undefined || raw === null || raw === "") {
    if (arg.required) throw new Error(`Missing ${where}.`);
    return arg.default;
  }
  switch (arg.type) {
    case "number": {
      const n = typeof raw === "number" ? raw : Number(String(raw).trim());
      if (!finite(n)) throw new Error(`${where} must be a number.`);
      if (arg.integer && !Number.isInteger(n)) throw new Error(`${where} must be a whole number.`);
      if (arg.min !== undefined && n < arg.min || arg.max !== undefined && n > arg.max) throw new Error(`${where} must be ${arg.min}–${arg.max}${arg.unit ? " " + arg.unit : ""}.`);
      return n;
    }
    case "string": {
      const s = String(raw);
      if (s.length > (arg.maxLength ?? 120)) throw new Error(`${where} is longer than ${arg.maxLength ?? 120} characters.`);
      return s;
    }
    case "enum": {
      const s = String(raw);
      if (!arg.values.includes(s)) throw new Error(`${where} must be one of ${arg.values.join(", ")}.`);
      return s;
    }
    case "boolean": {
      if (typeof raw === "boolean") return raw;
      const s = String(raw).toLowerCase();
      if (["true", "yes", "on", "1"].includes(s)) return true;
      if (["false", "no", "off", "0"].includes(s)) return false;
      throw new Error(`${where} must be true or false.`);
    }
    case "point": {
      const p = Array.isArray(raw) ? raw.map(Number) : String(raw).replace(/^@/, "").split(",").map(Number);
      if (p.length !== 2 || !p.every(finite)) throw new Error(`${where} must be a point such as @40,20.`);
      return p;
    }
    case "json": {
      const v = typeof raw === "string" ? JSON.parse(raw) : raw;
      if (v === null || typeof v !== "object") throw new Error(`${where} must be a JSON object or array.`);
      return structuredClone(v);
    }
    default: throw new Error(`Unsupported argument type ${arg.type}.`);
  }
}

// `pcb.place dip8 @40,20 name="Gate driver"`: named key=value pairs, an @x,y
// point, then positional values in declaration order.
export function tokenize(text) {
  const out = [], re = /([^\s"=]+=)?"((?:[^"\\]|\\.)*)"|\S+/g;
  for (const m of text.matchAll(re)) out.push(m[2] !== undefined ? (m[1] || "") + m[2].replace(/\\(.)/g, "$1") : m[0]);
  return out;
}

export class CommandBus {
  commands = new Map(); domains = new Map(); undoStack = []; redoStack = []; log = []; listeners = new Set(); queue = Promise.resolve();
  domain(id, handlers) { this.domains.set(id, { id, ...handlers }); return this; }
  register(spec) {
    if (!/^[a-z]+\.[a-z_]+$/.test(spec.op)) throw new Error(`Invalid command name ${spec.op}.`);
    this.commands.set(spec.op, { args: [], undoable: true, ...spec }); return this;
  }
  on(fn) { this.listeners.add(fn); return () => this.listeners.delete(fn); }
  emit(event) { for (const fn of this.listeners) try { fn(event); } catch (e) { console.error(e); } }
  list(domain) {
    return [...this.commands.values()].filter(c => !domain || c.domain === domain).map(c => ({ op: c.op, domain: c.domain, title: c.title,
      summary: c.summary, undoable: c.undoable && (!!c.touches || !!this.domains.get(c.domain)?.snapshot),
      args: c.args.map(describeArg), ...(c.example && { example: c.example }) }));
  }
  resolve(op, context) {
    if (this.commands.has(op)) return this.commands.get(op);
    if (context && this.commands.has(`${context}.${op}`)) return this.commands.get(`${context}.${op}`);
    const near = [...this.commands.keys()].filter(k => k.endsWith("." + op) || k.startsWith(op)).slice(0, 4);
    throw new Error(`Unknown command "${op}".${near.length ? " Did you mean " + near.join(", ") + "?" : " Type help to list commands."}`);
  }
  args(spec, raw = {}) {
    const known = new Set(spec.args.map(a => a.name));
    for (const k of Object.keys(raw)) if (!known.has(k)) throw new Error(`${spec.op} has no argument "${k}". Arguments: ${[...known].join(", ") || "none"}.`);
    return Object.fromEntries(spec.args.map(a => [a.name, coerce(a, raw[a.name])]).filter(([, v]) => v !== undefined));
  }
  parse(text, context) {
    const [head, ...tokens] = tokenize(text.trim());
    if (!head) throw new Error("Type a command, or help to list them.");
    const spec = this.resolve(head, context), raw = {}, positional = spec.args.filter(a => !a.named);
    for (const t of tokens) {
      const eq = t.indexOf("=");
      if (eq > 0 && /^[a-zA-Z_][a-zA-Z0-9_]*$/.test(t.slice(0, eq))) raw[t.slice(0, eq)] = t.slice(eq + 1);
      else if (t.startsWith("@")) { const a = spec.args.find(a => a.type === "point" && raw[a.name] === undefined); if (!a) throw new Error(`${spec.op} takes no point.`); raw[a.name] = t; }
      else { const a = positional.find(a => raw[a.name] === undefined); if (!a) throw new Error(`Unexpected value "${t}" for ${spec.op}.`); raw[a.name] = t; }
    }
    return { op: spec.op, args: raw };
  }
  complete(text, context) {
    const [head = "", ...rest] = tokenize(text);
    if (!rest.length && !text.endsWith(" ")) {
      const all = [...this.commands.values()];
      return all.filter(c => c.op.startsWith(head) || context && c.domain === context && c.op.slice(context.length + 1).startsWith(head))
        .sort((a, b) => (b.domain === context) - (a.domain === context) || a.op.localeCompare(b.op)).slice(0, 12)
        .map(c => ({ text: c.op, hint: c.title, usage: this.usage(c.op) }));
    }
    try { const spec = this.resolve(head, context), used = new Set(rest.map(t => t.split("=")[0]));
      return spec.args.filter(a => !used.has(a.name)).map(a => ({ text: `${a.name}=`, hint: `${a.type}${a.unit ? " · " + a.unit : ""}${a.values ? " · " + a.values.join("|") : ""}`, insert: true }));
    } catch { return []; }
  }
  usage(op) {
    const c = this.resolve(op);
    return `${c.op} ${c.args.map(a => (a.required ? "" : "[") + (a.type === "point" ? "@x,y" : a.named ? `${a.name}=…` : a.name) + (a.required ? "" : "]")).join(" ")}`.trim();
  }
  // Serialized: snapshots and the undo stack always describe one edit at a time.
  run(op, raw = {}, options = {}) { const p = this.queue.then(() => this.#run(op, raw, options)); this.queue = p.catch(() => {}); return p; }
  async #run(op, raw, { source = "ui", record = true, context } = {}) {
    const spec = this.resolve(op, context), args = this.args(spec, raw), touched = this.touched(spec, args);
    const undoable = record && spec.undoable && touched.length > 0, before = undoable ? touched.map(id => [id, this.domains.get(id).snapshot()]) : undefined;
    const entry = { time: new Date().toISOString(), op: spec.op, args: structuredClone(args), source, domain: spec.domain };
    try {
      const result = await spec.run(args, { bus: this, source });
      if (undoable) {
        // Slider drags arrive as many commands; one undo step covers a burst.
        const key = spec.coalesce?.(args), top = this.undoStack.at(-1), now = Date.now();
        if (key && top?.key === key && now - top.time < 1500) top.time = now;
        else { this.undoStack.push({ label: spec.title, snaps: before, key, time: now }); if (this.undoStack.length > 100) this.undoStack.shift(); }
        this.redoStack = [];
      }
      this.#record({ ...entry, ok: true, summary: summarize(result), mutates: !!undoable });
      this.emit({ type: "run", op: spec.op, domain: spec.domain, args, result, source });
      return result;
    } catch (e) {
      if (undoable) for (const [id, s] of before) try { this.domains.get(id).restore(s); } catch { /* the domain reported its own failure */ }
      this.#record({ ...entry, ok: false, error: e.message });
      this.emit({ type: "error", op: spec.op, domain: spec.domain, error: e, source });
      throw e;
    }
  }
  // All-or-nothing across domains: any failure restores every touched domain.
  batch(list, options = {}) { const p = this.queue.then(() => this.#batch(list, options)); this.queue = p.catch(() => {}); return p; }
  async #batch(list, { source = "script", label = "Script" } = {}) {
    if (!Array.isArray(list) || list.length > 500) throw new Error("A script is a list of up to 500 commands.");
    const snaps = new Map(), results = [];
    const take = domain => { const d = this.domains.get(domain); if (d?.snapshot && !snaps.has(domain)) snaps.set(domain, d.snapshot()); };
    try {
      for (const [i, step] of list.entries()) {
        const spec = this.resolve(step.op ?? "", step.context), args = this.args(spec, step.args || {});
        this.touched(spec, args).forEach(take);
        try { results.push({ op: spec.op, result: await spec.run(args, { bus: this, source }) }); }
        catch (e) { throw new Error(`Step ${i + 1} (${spec.op}): ${e.message}`); }
        this.#record({ time: new Date().toISOString(), op: spec.op, args, source, domain: spec.domain, ok: true, summary: summarize(results.at(-1).result), mutates: spec.undoable });
      }
    } catch (e) {
      for (const [domain, s] of snaps) this.domains.get(domain).restore(s);
      this.emit({ type: "error", op: "batch", error: e, source }); throw e;
    }
    if (snaps.size) { this.undoStack.push({ label, snaps: [...snaps] }); this.redoStack = []; }
    this.emit({ type: "batch", domains: [...snaps.keys()], results, source });
    return results;
  }
  // A multi-domain edit outside the registry (discard, checkout-free reverts)
  // with one undo step.
  async transaction(label, domains, fn, source = "ui") {
    const snaps = domains.filter(d => this.domains.get(d)?.snapshot).map(d => [d, this.domains.get(d).snapshot()]);
    try { const r = await fn(); this.undoStack.push({ label, snaps }); this.redoStack = []; this.#record({ time: new Date().toISOString(), op: label, args: {}, source, ok: true, mutates: false }); this.emit({ type: "transaction", label, domains }); return r; }
    catch (e) { for (const [d, s] of snaps) this.domains.get(d).restore(s); throw e; }
  }
  // The documents a command edits: its own domain, or those it names.
  touched(spec, args) { return (spec.touches?.(args) ?? [spec.domain]).filter(id => this.domains.get(id)?.snapshot); }
  #record(entry) { this.log.push(entry); if (this.log.length > 300) this.log.shift(); }
  #swap(from, to, type) {
    const e = from.pop(); if (!e) return false;
    const back = e.snaps.map(([domain]) => [domain, this.domains.get(domain).snapshot()]);
    for (const [domain, s] of [...e.snaps].reverse()) this.domains.get(domain).restore(s);
    to.push({ label: e.label, snaps: back });
    this.emit({ type, label: e.label, domains: e.snaps.map(([d]) => d) }); return true;
  }
  undo() { return this.#swap(this.undoStack, this.redoStack, "undo"); }
  redo() { return this.#swap(this.redoStack, this.undoStack, "redo"); }
  clearHistory() { this.undoStack = []; this.redoStack = []; this.emit({ type: "reset" }); }
  script() { return this.log.filter(e => e.ok && e.mutates).map(e => ({ op: e.op, args: e.args })); }
}
function summarize(result) {
  if (result === undefined || result === null) return "";
  if (typeof result === "string") return result.slice(0, 200);
  if (typeof result === "object" && typeof result.message === "string") return result.message.slice(0, 200);
  return "";
}
