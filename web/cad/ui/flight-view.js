// Flight studio: the vehicle's mass budget, rotor power and endurance.
import { project, onSolved } from "../app.js";
import { rotor, G } from "../domains/flight.js";
import { node, svg, button, format, numberField, selectField, textField, checkCards, certificates, assumptionList, section, linkedChip, live, lineChart, run, miniForm } from "./kit.js";
const root = document.getElementById("flight-studio");
let solution = null;
function render() {
  const d = project.store.get("flight"), r = solution?.results.flight, x = r?.details, status = node("p", "", "studio-status"); status.id = "flight-status";
  root.replaceChildren();
  const head = node("div", undefined, "studio-heading"), title = node("div");
  title.append(node("p", "FLIGHT / " + d.name.toUpperCase(), "panel-label"), node("h2", "Mass makes thrust, thrust makes current."),
    node("p", "Momentum theory sizes each rotor at hover and at the design thrust-to-weight point. The design-point ESC current drives the MOSFET model, and the ESC's sink and board come back as mass."));
  head.append(title); root.append(head);
  if (x) {
    const hero = node("div", undefined, "system-hero");
    for (const [label, v, u] of [["Take-off mass", x.mass, "kg"], ["Hover power", x.hover.P, "W"], ["Endurance", x.time, "min"], ["Design ESC current", x.design.esc, "A"], ["Hover ESC current", x.hover.esc, "A"],
      ...(d.propulsion ? [["Hover throttle", 100 * x.hover.throttle, "%"], ["Hover speed", x.hover.rpm, "rpm"], ["Thrust/weight, full", x.twFull, ""], ["Thrust/weight, empty", x.twEnd, ""], ["Implied figure of merit", x.hover.fm, ""]] : [["Induced velocity", x.hover.vi, "m/s"]]), ["Disk loading", x.diskLoading, "Pa"], ["Air density", x.rho, ""]]) { const c = node("div", undefined, "hero-stat"); c.append(node("span", label), node("b", u ? format(v, u, 4) : `${v.toFixed(4)} kg/m³`)); hero.append(c); }
    root.append(hero);
  }
  const grid = node("div", undefined, "domain-grid");
  // Plan view of the vehicle, to scale.
  const plan = section(grid, "Plan view", `${d.rotors} × ${(d.prop / 25.4).toFixed(1)}″ rotors`);
  const parts = project.store.get("parts").parts, arm = parts.find(p => p.id === "mounts")?.x ?? 368, span = arm + d.prop / 2 + 20;
  const s = svg("svg", { viewBox: `${-span} ${-span} ${2 * span} ${2 * span}`, class: "vehicle-plan", role: "img", "aria-label": "Vehicle plan view to scale" });
  for (let i = 0; i < d.rotors; i++) { const a = (45 + 360 * i / d.rotors) * Math.PI / 180, xx = arm * Math.cos(a), yy = arm * Math.sin(a);
    s.append(svg("line", { x1: 0, y1: 0, x2: xx, y2: yy, class: "plan-arm" }), svg("circle", { cx: xx, cy: yy, r: d.prop / 2, class: "plan-disc" }), svg("circle", { cx: xx, cy: yy, r: 18, class: "plan-motor" }),
      svg("rect", { x: .52 * xx - 17, y: .52 * yy - 9, width: 34, height: 18, class: "plan-esc", transform: `rotate(${a * 180 / Math.PI} ${.52 * xx} ${.52 * yy})` })); }
  s.append(svg("rect", { x: -95, y: -95, width: 190, height: 190, rx: 12, class: "plan-body" }), svg("circle", { cx: 0, cy: 70, r: 22, class: "plan-camera" }), svg("text", { x: 0, y: 5, "text-anchor": "middle", class: "plan-label" }, d.name));
  plan.append(s, node("p", `Rotor discs are to scale, with ${format(x?.clearance, "mm", 3)} between tips. ESCs sit on the arms in the downwash that cools them.`, "small muted"));
  // Power against thrust-to-weight.
  if (x) {
    const pts = [], amps = [];
    for (let tw = .6; tw <= Math.max(2.6, d.designTW + .4); tw += .05) { const T = tw * x.weight / d.rotors, rr = rotor(T, x.rho, x.area, d.fm), P = d.rotors * rr.shaft / d.motorEfficiency / d.escEfficiency + d.avionicsPower; pts.push([tw, P / 1000]); amps.push([tw, P / x.V]); }
    const ch = section(grid, "Power and current", d.propulsion ? "momentum theory for comparison: measured data sets the numbers" : "momentum theory · constant FM and efficiencies");
    ch.append(lineChart({ series: [{ points: pts }], x: { label: "thrust / weight" }, y: { label: "electrical power, kW" }, marks: [{ x: 1, y: x.hover.P / 1000, label: "hover" }, { x: d.designTW, y: x.design.P / 1000, label: "design" }] }));
    ch.append(node("p", `Pack: ${d.battery.cells}S, ${format(x.V, "V", 3)} nominal, ${format(x.design.Vt, "V", 3)} under the design load (${format(x.design.I, "A", 3)}). Ideal hover power is T·vᵢ; the figure of merit ${d.fm} and efficiencies turn it into battery draw.`, "small muted"));
  }
  root.append(grid);
  // Mass budget.
  const mb = section(root, "Mass budget", x ? `${format(x.mass, "kg", 4)} take-off` : "");
  if (x) {
    const linked = [["ESCs (sink + board) × " + d.rotors, (solution.applied["flight.esc_mass"]?.value ?? 0) * d.rotors, "link"], ["Camera optics", solution.applied["flight.camera_mass"]?.value ?? 0, "link"], ["Airframe parts", solution.applied["flight.parts_mass"]?.value ?? 0, "link"]];
    const items = [...d.masses.map(m => [m.name, m.mass, "fixed"]), ["Battery", d.battery.mass, "battery"], ...linked], bar = node("div", undefined, "mass-bar");
    for (const [name, m, kind] of items) { const seg = node("i"); seg.style.width = `${100 * m / x.mass}%`; seg.dataset.kind = kind; seg.title = `${name}: ${format(m, "kg", 4)}`; bar.append(seg); }
    mb.append(bar);
    const table = node("table", undefined, "data-table"), body = node("tbody");
    for (const m of d.masses) { const tr = node("tr"), name = node("td", m.name), mass = node("td"), del = node("td");
      const i = node("input"); i.type = "number"; i.min = 0; i.max = 500; i.step = .01; i.value = m.mass; i.setAttribute("aria-label", `${m.name} mass in kg`); i.id = `flight-mass-${m.id}`;
      i.addEventListener("change", () => i.checkValidity() && run("flight.mass", { id: m.id, mass: Number(i.value) }, status).catch(() => {}));
      mass.append(i, node("small", " kg")); del.append(button("Remove", () => run("flight.remove_mass", { id: m.id }, status).catch(() => {}))); if (m.part) { const small = node("small", m.part); if (m.source) { const a = node("a", " source"); a.href = m.source; a.rel = "noopener"; small.append(a); } name.append(small); } tr.append(name, mass, del); body.append(tr); }
    for (const [name, m] of linked) { const tr = node("tr", undefined, "linked-row"); tr.append(node("td", "⛓ " + name), node("td", format(m, "kg", 4)), node("td", "linked")); body.append(tr); }
    table.append(body); mb.append(table);
    const add = node("div", undefined, "rtl-actions"), nm = node("input"); nm.placeholder = "New item name"; nm.maxLength = 60; nm.setAttribute("aria-label", "New mass item name");
    const kg = node("input"); kg.type = "number"; kg.min = 0; kg.step = .01; kg.placeholder = "kg"; kg.setAttribute("aria-label", "New mass item kg");
    const pn = node("input"); pn.placeholder = "Maker and part number"; pn.maxLength = 80; pn.setAttribute("aria-label", "New mass item part number");
    const sr = node("input"); sr.placeholder = "Source link (https)"; sr.maxLength = 300; sr.setAttribute("aria-label", "New mass item source");
    add.append(nm, kg, pn, sr, button("＋ Add mass", () => { const id = "m" + Date.now().toString(36); run("flight.mass", { id, mass: Number(kg.value || 0), name: nm.value || "Item", ...(pn.value && { part: pn.value }), ...(sr.value && { source: sr.value }) }, status).catch(() => {}); }));
    mb.append(add);
  }
  propulsion(d, x, status);
  // Parameters.
  const p = section(root, "Vehicle and mission");
  const fields = node("div", undefined, "field-grid");
  selectField(fields, { label: "Rotors", id: "flight-rotors", value: String(d.rotors), options: [["3", "Tricopter"], ["4", "Quad"], ["6", "Hexa"], ["8", "Octo"]], op: "flight.set", args: v => ({ rotors: v }), status });
  for (const [key, label, unit, min, max, step, value] of [["prop", "Propeller diameter", "mm", 50, 2000, 1, d.prop], ["fm", "Figure of merit", "", .3, .85, .01, d.fm], ["motor_efficiency", "Motor efficiency", "", .3, .98, .01, d.motorEfficiency],
    ["esc_efficiency", "ESC efficiency", "", .5, .995, .005, d.escEfficiency], ["design_tw", "Design thrust/weight", "", 1.1, 5, .05, d.designTW], ["min_tw", "Minimum thrust/weight", "", 1, 5, .05, d.minTW],
    ["avionics", "Avionics power", "W", 0, 5000, 1, d.avionicsPower], ["altitude", "Mission altitude", "m", 1, 5000, 1, d.altitude], ["required_time", "Required endurance", "min", 0, 600, 1, d.requiredTime], ["capture", "Wake capture (ESC duct)", "", 0, 1, .01, d.capture]])
    numberField(fields, { label, id: `flight-${key}`, unit, value, min, max, step, op: "flight.set", args: v => ({ [key]: v }), status });
  numberField(fields, { label: "Rotor arm from centre", id: "flight-arm", unit: "mm", value: d.arm, min: 20, max: 5000, step: 1, op: "flight.set", args: v => ({ arm: v }), status });
  textField(fields, { label: "Vehicle name", id: "flight-name", value: d.name, op: "flight.set", args: v => ({ name: v }), status });
  numberField(fields, { label: "Air pressure", id: "flight-pressure", unit: "Pa", value: d.pressure, min: 5e4, max: 1.1e5, step: 100, op: "flight.set", args: v => ({ pressure: v }), status });
  const t = numberField(fields, { label: "Air temperature", id: "flight-temperature", unit: "degC", value: d.temperature, min: -40, max: 55, step: 1, op: "flight.set", args: v => ({ temperature: v }), status });
  const chip = linkedChip(solution, "flight.air_temperature"); if (chip) { t.disabled = true; t.parentElement.append(chip); }
  p.append(fields);
  const b = section(root, "Battery");
  const bf = node("div", undefined, "field-grid");
  for (const [key, label, unit, min, max, step, value] of [["cells", "Cells in series", "", 1, 24, 1, d.battery.cells], ["capacity", "Capacity", "", .1, 200, .5, d.battery.capacity], ["mass", "Pack mass", "kg", .01, 100, .01, d.battery.mass], ["resistance", "Pack resistance", "ohm", 0, 1, .001, d.battery.resistance], ["usable", "Usable fraction", "", .3, 1, .01, d.battery.usable]])
    numberField(bf, { label: key === "capacity" ? "Capacity (Ah)" : label, id: `flight-battery-${key}`, unit, value, min, max, step, op: "flight.battery", args: v => ({ [key]: v }), status });
  b.append(bf, status);
  if (r) { const c = section(root, "Flight checks"); checkCards(c, r.checks || []); certificates(c, r.certificates); }
  const a = section(root, "Assumptions"); assumptionList(a, "flight");
}
// The measured table, shown as it was published, and how to use your own.
function propulsion(d, x, status) {
  const s = section(root, "Propulsion", d.propulsion ? `${d.propulsion.motor} · ${d.propulsion.propeller}` : "momentum theory with a stated figure of merit"), p = d.propulsion;
  if (p) { const table = node("table", undefined, "data-table"), head = node("thead"), tr = node("tr"), body = node("tbody"); for (const h of ["Throttle", "Volts", "Amps", "Watts", "rpm", "Torque, N·m", "Thrust, g", "g/W"]) tr.append(node("th", h)); head.append(tr); table.append(head);
    for (const r of p.table) { const row = node("tr"); for (const v of [`${Math.round(r.throttle * 100)}%`, r.volts, r.amps, r.watts, r.rpm, r.torque, r.grams, (r.grams / r.watts).toFixed(2)]) row.append(node("td", String(v))); body.append(row); } table.append(body); s.append(table);
    s.append(node("p", `KV ${p.kv}, ${p.resistance} Ω, ${p.peakCurrent} A and ${p.maxPower} W limits. Bench air taken as ${p.density} kg/m³; ${Math.round(p.blockage * 100)}% of thrust lost to the arm. Speed, torque and power at any thrust come from this table, corrected for air density by similarity.`, "small muted"));
    if (p.source) { const a = node("a", "Measured data source"); a.href = p.source; a.rel = "noopener"; s.append(a); } s.append(button("Use momentum theory instead", () => run("flight.clear_propulsion", {}, status).catch(() => {}), "", "flight-clear-propulsion")); }
  miniForm(s, { id: "flight-propulsion", action: p ? "Replace the table" : "Use measured data", status, fields: [{ name: "motor", label: "Motor", value: p?.motor }, { name: "propeller", label: "Propeller", value: p?.propeller }, { name: "kv", label: "KV (rpm/V)", number: true, value: p?.kv }, { name: "resistance", label: "Resistance (Ω)", number: true, value: p?.resistance }, { name: "peak_current", label: "Peak current (A)", number: true, value: p?.peakCurrent }, { name: "max_power", label: "Maximum power (W)", number: true, value: p?.maxPower }, { name: "rotor_inertia", label: "Rotor inertia (kg·m²)", number: true, value: p?.rotorInertia }, { name: "density", label: "Test air density", number: true, value: p?.density }, { name: "blockage", label: "Thrust lost to the arm (0–0.5)", number: true, value: p?.blockage }, { name: "source", label: "Source (https)", wide: true, value: p?.source },
    { name: "rows", label: "Rows: throttle% volts amps watts rpm torque grams, one per line", long: true, wide: true, value: p?.table.map(r => [Math.round(r.throttle * 100), r.volts, r.amps, r.watts, r.rpm, r.torque, r.grams].join(" ")).join("\n") }],
    submit: ({ rows, ...a }) => { const table = (rows ?? "").trim().split(/\r?\n/).map(l => l.trim().split(/[\s,;]+/).map(Number)).filter(r => r.length >= 7 && r.every(Number.isFinite)).map(([t, volts, amps, watts, rpm, torque, grams]) => ({ throttle: t / 100, volts, amps, watts, rpm, torque, grams })); if (table.length < 2) throw new Error("Paste at least two rows of seven numbers."); return run("flight.propulsion", { ...a, table }, status); } });
}
onSolved(s => { solution = s; });
live(["flight", "parts"], render);
render();
