// Structures studio: an arm as a cantilevered tube, against the rotor's excitation.
import { project, onSolved } from "../app.js";
import { armResults } from "../domains/structures.js";
import { node, svg, format, numberField, textField, checkCards, assumptionList, section, linkedChip, live } from "./kit.js";
const root = document.getElementById("structures-studio");
let solution = null;
const FIELDS = [["od", "Tube outer diameter", "mm", .5, "od", "od"], ["wall", "Wall", "mm", .1, "wall", "wall"], ["length", "Free length", "mm", 1, "length", "length"], ["load_thrust", "Peak thrust", "N", 1, "loadThrust", "thrust"], ["load_torque", "Peak torque", "Nm", .05, "loadTorque", "torque"], ["tip_mass", "Mass at the tip", "kg", .01, "tipMass", "tip_mass"]];
const MATERIAL = [["modulus", "Axial modulus", "GPa", 1, "modulus"], ["density", "Density (g/cm³)", "", .01, "density"], ["strength", "Allowable stress", "MPa", 10, "strength"], ["safety_factor", "Safety factor", "", .1, "safetyFactor"], ["deflection_limit", "Deflection limit", "mm", .1, "deflectionLimit"], ["separation", "Mode separation", "%", 1, "separation"], ["blades", "Blades per rotor", "", 1, "blades"]];
// Campbell view: the mode against the bands the rotors sweep.
function campbell(x) {
  const top = Math.max(x.frequency * 1.3, ...x.bands.map(b => b.hi * 1.2), 10), W = 640, H = 92, px = v => 20 + (W - 40) * v / top, s = svg("svg", { viewBox: `0 0 ${W} ${H}`, class: "campbell", role: "img", "aria-label": `First mode ${x.frequency.toFixed(0)} Hz against rotor excitation` });
  x.bands.forEach((b, i) => { s.append(svg("rect", { x: px(b.lo), y: 16 + i * 26, width: px(b.hi) - px(b.lo), height: 20, class: x.conflicts.includes(b) ? "band-hit" : "band" }), svg("text", { x: px(b.hi) + 6, y: 30 + i * 26, class: "chart-tick" }, `${b.name} ${b.lo.toFixed(0)}–${b.hi.toFixed(0)} Hz`)); });
  s.append(svg("line", { x1: px(x.frequency), x2: px(x.frequency), y1: 8, y2: H - 8, class: "mode-line" }), svg("text", { x: px(x.frequency) + 5, y: H - 8, class: "chart-tick" }, `mode ${x.frequency.toFixed(0)} Hz`));
  return s;
}
function render() {
  const d = project.store.get("structures"), r = solution?.results.structures, x = r?.details ?? armResults(d), status = node("p", "", "studio-status"); status.id = "structures-status";
  root.replaceChildren();
  const head = node("div", undefined, "studio-heading"), title = node("div");
  title.append(node("p", `STRUCTURES / ${project.store.docs.project.name.toUpperCase()}`, "panel-label"), node("h2", "Strength is rarely what breaks an arm."),
    node("p", "A cantilevered tube carrying thrust and torque at its tip: bending and torsion combined by von Mises, and the first bending mode by Rayleigh's method. The mode must clear the rotor's once-per-rev and blade-pass across the whole flying speed range. The root is taken as fully clamped, an upper bound: a tap test on the real joint gives the true mode."));
  head.append(title); root.append(head);
  const hero = node("div", undefined, "system-hero");
  for (const [label, v] of [["Peak stress", format(x.vonMises, "MPa", 3)], ["Strength margin", x.margin.toFixed(1)], ["Tip deflection", format(x.deflection, "mm", 3)], ["First mode", format(x.frequency, "Hz", 4)], ["Arm mass", format(x.armMass, "kg", 3)]]) { const c = node("div", undefined, "hero-stat"); c.append(node("span", label), node("b", v)); hero.append(c); }
  root.append(hero);
  if (x.bands.length) section(root, "Mode against rotor excitation", x.conflicts.length ? "inside a band" : "clear of every band").append(campbell(x));
  const g = section(root, d.name, "linked values override these"), f = node("div", undefined, "field-grid");
  textField(f, { label: "Name", id: "structures-name", value: d.name, op: "structures.set", args: v => ({ name: v }), status });
  for (const [arg, label, unit, step, key, port] of FIELDS) { const i = numberField(f, { label, id: `structures-${arg}`, unit, value: d[key], step, op: "structures.set", args: v => ({ [arg]: v }), status }), chip = linkedChip(solution, `structures.${port}`); if (chip) { i.disabled = true; i.parentElement.append(chip); } }
  for (const [arg, label, unit, step, key] of MATERIAL) numberField(f, { label, id: `structures-${arg}`, unit, value: d[key], step, op: "structures.set", args: v => ({ [arg]: v }), status });
  textField(f, { label: "Material source", id: "structures-source", value: d.source, max: 300, op: "structures.set", args: v => ({ source: v }), status });
  g.append(f, status);
  if (r) { const c = section(root, "Structural checks"); checkCards(c, r.checks || []); }
  const a = section(root, "Assumptions"); assumptionList(a, "structures");
}
onSolved(s => { solution = s; });
live(["structures"], render);
render();
