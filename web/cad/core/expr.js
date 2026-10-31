// Arithmetic over named values, parsed by hand: a design file must never run
// code, and the page policy forbids eval. + - * / ^, unary minus, parentheses,
// the constants pi and e, and the functions below.
const FUNCTIONS = { sqrt: Math.sqrt, abs: Math.abs, exp: Math.exp, log: Math.log, log10: Math.log10, sin: Math.sin, cos: Math.cos, tan: Math.tan, asin: Math.asin, acos: Math.acos, atan: Math.atan,
  floor: Math.floor, ceil: Math.ceil, round: Math.round, min: Math.min, max: Math.max, pow: Math.pow, hypot: Math.hypot };
export const CONSTANTS = { pi: Math.PI, e: Math.E };
export const FUNCTION_NAMES = Object.keys(FUNCTIONS);
const TOKEN = /\s*(?:(\d+\.?\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?)|([a-zA-Z_][a-zA-Z0-9_]*)|(.))/y;
export function parse(text) {
  if (typeof text !== "string" || !text.trim() || text.length > 400) throw new Error("An expression is up to 400 characters.");
  const tokens = []; TOKEN.lastIndex = 0;
  while (TOKEN.lastIndex < text.length) { const m = TOKEN.exec(text); if (!m || m[0].trim() === "" ) break; tokens.push(m[1] !== undefined ? { num: Number(m[1]) } : m[2] !== undefined ? { name: m[2] } : { op: m[3] }); }
  let i = 0, depth = 0;
  const peek = op => tokens[i]?.op === op, take = op => peek(op) && ++i;
  const need = op => { if (!take(op)) throw new Error(`Expected "${op}" in ${text}`); };
  const guard = fn => (...a) => { if (++depth > 64) throw new Error("That expression is nested too deeply."); const r = fn(...a); depth--; return r; };
  const atom = guard(() => {
    const t = tokens[i++]; if (!t) throw new Error(`Incomplete expression: ${text}`);
    if (t.num !== undefined) return t.num;
    if (t.op === "(") { const v = sum(); need(")"); return v; }
    if (t.op === "-") return { neg: power() };
    if (t.op === "+") return power();
    if (t.op) throw new Error(`Unexpected "${t.op}" in ${text}`);
    if (!take("(")) return { ref: t.name };
    if (!FUNCTIONS[t.name]) throw new Error(`Unknown function ${t.name}. Functions: ${FUNCTION_NAMES.join(", ")}.`);
    const args = [sum()]; while (take(",")) args.push(sum()); need(")"); return { call: t.name, args };
  });
  const power = guard(() => { const base = atom(); return take("^") ? { op: "^", a: base, b: power() } : base; });
  const product = guard(() => { let v = power(); for (;;) { const op = peek("*") ? "*" : peek("/") ? "/" : null; if (!op) return v; i++; v = { op, a: v, b: power() }; } });
  const sum = guard(() => { let v = product(); for (;;) { const op = peek("+") ? "+" : peek("-") ? "-" : null; if (!op) return v; i++; v = { op, a: v, b: product() }; } });
  const tree = sum(); if (i < tokens.length) throw new Error(`Unexpected text after the expression: ${text}`);
  return tree;
}
export function names(tree, out = new Set()) {
  if (typeof tree !== "object") return out;
  if (tree.ref) out.add(tree.ref); if (tree.neg !== undefined) names(tree.neg, out); if (tree.op) { names(tree.a, out); names(tree.b, out); } if (tree.args) tree.args.forEach(a => names(a, out));
  return out;
}
export function evaluate(tree, scope) {
  if (typeof tree === "number") return tree;
  if (tree.ref) { const v = Object.hasOwn(scope, tree.ref) ? scope[tree.ref] : CONSTANTS[tree.ref]; if (v === undefined) throw new Error(`Unknown name ${tree.ref}.`); return v; }
  if (tree.neg !== undefined) return -evaluate(tree.neg, scope);
  if (tree.call) return FUNCTIONS[tree.call](...tree.args.map(a => evaluate(a, scope)));
  const a = evaluate(tree.a, scope), b = evaluate(tree.b, scope);
  return tree.op === "+" ? a + b : tree.op === "-" ? a - b : tree.op === "*" ? a * b : tree.op === "/" ? a / b : a ** b;
}
