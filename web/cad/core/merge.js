// Three-way merge of design documents. Objects merge field by field; arrays of
// objects with unique string ids merge by id, so two people can move different
// components or add different rooms on separate branches without a conflict.
const sorted = v => Array.isArray(v) ? v.map(sorted) : v && typeof v === "object" ? Object.fromEntries(Object.keys(v).sort().map(k => [k, sorted(v[k])])) : v;
const same = (a, b) => JSON.stringify(sorted(a)) === JSON.stringify(sorted(b));
const plain = v => v !== null && typeof v === "object" && !Array.isArray(v);
const keyed = v => Array.isArray(v) && v.every(x => plain(x) && typeof x.id === "string") && new Set(v.map(x => x.id)).size === v.length;

export function merge3(base, target, source, choices = {}, path = "") {
  const conflicts = [];
  const value = walk(base, target, source, path, choices, conflicts);
  return { value, conflicts };
}
function walk(base, target, source, path, choices, conflicts) {
  if (same(source, base) || same(target, source)) return clone(target);
  if (same(target, base)) return clone(source);
  if (plain(base ?? {}) && plain(target) && plain(source)) {
    const b = base ?? {}, out = {};
    for (const k of new Set([...Object.keys(target), ...Object.keys(source), ...Object.keys(b)])) {
      const v = walk(b[k], target[k], source[k], path ? `${path}.${k}` : k, choices, conflicts);
      if (v !== undefined) out[k] = v;
    }
    return out;
  }
  if (keyed(base ?? []) && keyed(target) && keyed(source)) {
    const index = a => new Map(a.map(x => [x.id, x])), b = index(base ?? []), t = index(target), s = index(source), out = [];
    for (const id of [...t.keys(), ...[...s.keys()].filter(id => !t.has(id))]) {
      const v = walk(b.get(id), t.get(id), s.get(id), `${path}[${id}]`, choices, conflicts);
      if (v !== undefined) out.push(v);
    }
    return out;
  }
  const choice = choices[path];
  conflicts.push({ path, base: clone(base), target: clone(target), source: clone(source), resolved: !!choice });
  return clone(choice === "source" ? source : target);
}
const clone = v => v === undefined ? undefined : structuredClone(v);
