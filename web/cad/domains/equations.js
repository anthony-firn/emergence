// Your own relations. Inputs are named values a link can drive; outputs are
// expressions over the inputs and earlier outputs, published as ordinary typed
// ports. This is how a project says what no domain model knows: a rail current
// from a die power and a voltage, a cost, a margin. Units are as you declare
// them; the arithmetic does not check them.
import { parse, names, evaluate, CONSTANTS } from "../core/expr.js";
import { UNITS } from "../core/units.js";
const NAME = /^[a-zA-Z_][a-zA-Z0-9_]{0,31}$/, finite = (v, lo = -1e15, hi = 1e15) => Number.isFinite(v) && v >= lo && v <= hi;
export const equationsExample = () => ({ inputs: [], outputs: [] });
const label = (v, id) => { if (v !== undefined && (typeof v !== "string" || v.length > 60)) throw new Error(`${id}: a label is up to 60 characters.`); return v?.trim() || id; };
const unit = (v, id) => { if (!Object.hasOwn(UNITS, v ?? "")) throw new Error(`${id}: unknown unit ${v}. Units: ${Object.keys(UNITS).filter(Boolean).join(", ")}.`); return v ?? ""; };
export function validateEquations(d) {
  if (!d || !Array.isArray(d.inputs) || d.inputs.length > 64 || !Array.isArray(d.outputs) || d.outputs.length > 64) throw new Error("Invalid equations document.");
  const seen = new Set(), name = id => { if (typeof id !== "string" || !NAME.test(id) || Object.hasOwn(CONSTANTS, id)) throw new Error(`${id} cannot name a value: use letters, digits and _, starting with a letter, and not pi or e.`); if (seen.has(id)) throw new Error(`${id} is defined twice.`); seen.add(id); return id; };
  const inputs = d.inputs.map(x => { name(x?.id); if (!finite(x.value)) throw new Error(`${x.id}: value must be a number.`);
    for (const k of ["min", "max"]) if (x[k] !== undefined && !finite(x[k])) throw new Error(`${x.id}: ${k} must be a number.`);
    return { id: x.id, label: label(x.label, x.id), unit: unit(x.unit, x.id), value: x.value, ...(x.min !== undefined && { min: x.min }), ...(x.max !== undefined && { max: x.max }) }; });
  const outputs = d.outputs.map(x => { const known = new Set(seen), tree = parse(x?.expr); name(x.id);
    for (const r of names(tree)) if (!known.has(r) && !Object.hasOwn(CONSTANTS, r)) throw new Error(`${x.id}: ${r} is not an input or an earlier output.`);
    return { id: x.id, label: label(x.label, x.id), unit: unit(x.unit, x.id), expr: x.expr }; });
  return { inputs, outputs };
}
export function equationsResults(d, overrides = {}) {
  const scope = Object.fromEntries(d.inputs.map(x => [x.id, overrides[x.id] ?? x.value])), values = {}, failed = [];
  for (const x of d.outputs) { const v = evaluate(parse(x.expr), scope); scope[x.id] = v; values[x.id] = v; if (!Number.isFinite(v)) failed.push(x); }
  return { values, scope, failed };
}
export const equationsAssumptions = d => d.outputs.length ? [{ id: "equations.units", label: "Equation units are as declared", kind: "prescribed", source: "Values convert at links by the declared unit; the arithmetic inside an equation is not dimension-checked." }] : [];
