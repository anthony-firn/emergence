// Power studio: a three-phase inverter from its MOSFET's datasheet.
import { project, onSolved } from "../app.js";
import { MOSFET_FIELDS, CAPACITOR_FIELDS } from "../domains/power.js";
import { node, format, numberField, textField, checkCards, certificates, assumptionList, section, linkedChip, live } from "./kit.js";
const root = document.getElementById("power-studio");
let solution = null;
const STAGE = [["pwm", "PWM frequency", "Hz", 1e3, 1e6, 1000, "pwm"], ["deadtime", "Dead time (ns)", "", 0, 1e4, 10, "deadtime"], ["gate_drive", "Gate drive", "V", 4, 20, .5, "gateDrive"], ["gate_resistor", "External gate resistor (Ω)", "", 0, 100, .1, "gateResistor"], ["driver_peak", "Driver peak", "A", .1, 100, .1, "driverPeak"], ["ringing", "Ringing allowance (× bus)", "", 1, 3, .05, "ringing"], ["efficiency_limit", "Efficiency limit", "%", 0, 100, .5, "efficiencyLimit"]];
const OPERATING = [["phase", "Phase current", "A", 0, 1e4, .1, "phase", "phase_current"], ["bus", "Bus voltage", "V", 1, 2000, .1, "bus", "bus_voltage"], ["duty", "Average duty", "%", 0, 100, 1, "duty", "duty"], ["junction", "Junction", "degC", -40, 200, 1, "junction", "junction"], ["peak_bus", "Peak bus", "V", 1, 2000, .1, "peakBus", "peak_voltage"]];
const FET = { vds: ["Drain-source rating (V)", 1], rds25: ["RDS(on) at 25 °C (mΩ)", .1], rdsTc: ["RDS(on) coefficient (1/K)", .0005], qg: ["Total gate charge (nC)", 1], qsw: ["Switching charge (nC)", 1], qoss: ["Output charge (nC)", 1], qrr: ["Reverse recovery charge (nC)", 1], vplateau: ["Gate plateau (V)", .1], rg: ["Internal gate resistance (Ω)", .1], vsd: ["Body diode drop (V)", .01], rthJc: ["Junction to case (K/W)", .01], tjMax: ["Maximum junction (°C)", 1] };
const CAP = { count: ["Capacitors in the bank", 1], capacitance: ["Capacitance each (µF)", 1], esr: ["ESR each (mΩ)", 1], rating: ["Ripple rating each (A)", .1], voltage: ["Rated voltage (V)", 1], grams: ["Mass each (g)", .1], limit: ["Ripple limit (% of bus)", .1] };
const ARG = { rdsTc: "rds_tc", rthJc: "rth_jc", tjMax: "tj_max" };
function render() {
  const d = project.store.get("power"), r = solution?.results.power, x = r?.details, status = node("p", "", "studio-status"); status.id = "power-status";
  root.replaceChildren();
  const head = node("div", undefined, "studio-heading"), title = node("div");
  title.append(node("p", `POWER / ${project.store.docs.project.name.toUpperCase()}`, "panel-label"), node("h2", "Where the watts go in the inverter."),
    node("p", "Block commutation with synchronous rectification: two MOSFETs carry the phase current at any moment and one leg hard-switches each PWM cycle. Conduction uses RDS(on) at the junction temperature; switching uses the gate charges and the gate current they imply. It is an analytic estimate from the datasheet: a double-pulse test on the real board measures the truth."));
  head.append(title); root.append(head);
  if (x) { const hero = node("div", undefined, "system-hero");
    for (const [label, v] of [["Inverter loss", format(x.total, "W", 3)], ["Efficiency", format(x.efficiency, "%", 4)], ["Conduction", format(x.conduction, "W", 3)], ["Switching", format(x.switching, "W", 3)], ["Gate drive", format(x.gate, "W", 3)], ["Dead time", format(x.dead, "W", 3)], [`Hottest FET (${x.hottest.device})`, format(x.hottest.rms, "A", 3)], ["Bus ripple", format(x.ripple, "A", 3)]]) { const c = node("div", undefined, "hero-stat"); c.append(node("span", label), node("b", v)); hero.append(c); }
    root.append(hero); }
  const op = section(root, "Operating point", "linked values override these"), f = node("div", undefined, "field-grid");
  for (const [arg, label, unit, min, max, step, key, port] of OPERATING) { const i = numberField(f, { label, id: `power-${arg}`, unit, value: d[key], min, max, step, op: "power.set", args: v => ({ [arg]: v }), status }), chip = linkedChip(solution, `power.${port}`); if (chip) { i.disabled = true; i.parentElement.append(chip); } }
  op.append(f);
  const s = section(root, d.name, "gate drive and switching"), g = node("div", undefined, "field-grid");
  textField(g, { label: "Name", id: "power-name", value: d.name, op: "power.set", args: v => ({ name: v }), status });
  for (const [arg, label, unit, min, max, step, key] of STAGE) numberField(g, { label, id: `power-${arg}`, unit, value: d[key], min, max, step, op: "power.set", args: v => ({ [arg]: v }), status });
  s.append(g);
  const m = section(root, "MOSFET", d.mosfet.part), h = node("div", undefined, "field-grid");
  textField(h, { label: "Part", id: "power-part", value: d.mosfet.part, max: 80, op: "power.mosfet", args: v => ({ part: v }), status });
  textField(h, { label: "Datasheet", id: "power-source", value: d.mosfet.source ?? "", max: 300, op: "power.mosfet", args: v => ({ source: v }), status });
  for (const [key, [label, step]] of Object.entries(FET)) numberField(h, { label, id: `power-fet-${key}`, value: d.mosfet[key], min: MOSFET_FIELDS[key][0], max: MOSFET_FIELDS[key][1], step, op: "power.mosfet", args: v => ({ [ARG[key] ?? key]: v }), status });
  m.append(h);
  const cap = section(root, "Bus capacitors", d.capacitor.part), cf = node("div", undefined, "field-grid");
  textField(cf, { label: "Part", id: "power-cap-part", value: d.capacitor.part, max: 80, op: "power.capacitor", args: v => ({ part: v }), status });
  textField(cf, { label: "Datasheet", id: "power-cap-source", value: d.capacitor.source ?? "", max: 300, op: "power.capacitor", args: v => ({ source: v }), status });
  for (const [key, [label, step]] of Object.entries(CAP)) numberField(cf, { label, id: `power-cap-${key}`, value: d.capacitor[key], min: CAPACITOR_FIELDS[key][0], max: CAPACITOR_FIELDS[key][1], step, op: "power.capacitor", args: v => ({ [key]: v }), status });
  cap.append(cf, status);
  if (r) { const c = section(root, "Power checks"); checkCards(c, r.checks || []); certificates(c, r.certificates); }
  const a = section(root, "Assumptions"); assumptionList(a, "power");
}
onSolved(s => { solution = s; });
live(["power"], render);
render();
