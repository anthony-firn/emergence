// Controls studio: PX4's cascade on this airframe, in frequency and in flight.
import { project, onSolved } from "../app.js";
import { paramFile, PX4_DEFAULTS, PX4_SOURCE } from "../domains/controls.js";
import { slug } from "../core/project.js";
import { node, button, format, numberField, checkCards, assumptionList, section, live, lineChart, download, run, toast } from "./kit.js";
const root = document.getElementById("controls-studio");
let solution = null;
const GROUPS = [["Attitude", /^MC_(ROLL|PITCH|YAW)_P$/], ["Rate", /RATE_(P|I|D|K)$/], ["Rate limits", /RATE_MAX$/], ["Position and velocity", /^MPC_/], ["Filters", /^IMU_/]];
const PLANT = [["mass", "Mass", "kg", "mass", "mass"], ["ixx", "Roll inertia", "kgm2", "ixx", "ixx"], ["iyy", "Pitch inertia", "kgm2", "iyy", "iyy"], ["izz", "Yaw inertia", "kgm2", "izz", "izz"], ["arm", "Arm", "mm", "arm", "arm"], ["rotors", "Rotors", "", "rotors", "rotors"], ["hover_throttle", "Hover throttle", "%", "hoverThrottle", "hover_throttle"], ["max_thrust", "Thrust per rotor", "N", "maxThrust", "max_thrust"], ["tau", "Motor lag", "s", "tau", "motor_tau"], ["kappa", "Torque per thrust", "m", "kappa", "torque_ratio"], ["drag_area", "Drag area", "m2", "dragArea", null]];
const series = (rows, scale = 1) => rows.map(r => [r.t, r.y * scale]);
function render() {
  const d = project.store.get("controls"), r = solution?.results.controls, x = r?.details, T = x?.tests, status = node("p", "", "studio-status"); status.id = "controls-status";
  root.replaceChildren();
  const head = node("div", undefined, "studio-heading"), title = node("div");
  title.append(node("p", `CONTROLS / ${project.store.docs.project.name.toUpperCase()}`, "panel-label"), node("h2", "Will it actually fly?"),
    node("p", "PX4's multicopter cascade, with its parameter names and defaults, flown in a 6-DoF simulation at its 400 Hz rate loop with its gyro and D-term filters. Motors lag by a time constant from real motor constants. This is a model of PX4's controller, not PX4 itself; export the parameters and confirm in software-in-the-loop and then in flight."));
  const tools = node("div", undefined, "studio-file-tools");
  tools.append(button("Tune for this airframe", () => run("controls.tune", {}, status).catch(() => {}), "accent", "controls-tune"), button("↓ PX4 parameters", () => { const name = slug(project.store.docs.project.name); download(paramFile(d, project.store.docs.project.name), `${name}.params`, "text/plain"); toast("A QGroundControl parameter file. Review it before loading it onto an aircraft."); }, "", "controls-export"));
  head.append(title, tools); root.append(head);
  if (x) { const a = x.roll, hero = node("div", undefined, "system-hero");
    for (const [label, v] of [["Rate bandwidth", format(a.bandwidth, "Hz", 3)], ["Phase margin", `${a.phaseMargin?.toFixed(1)}°`], ["Gain margin", `${a.gainMargin?.toFixed(1)} dB`], ["Roll step rise", format(T?.rollStep.rise, "s", 2)], ["Roll overshoot", format(T?.rollStep.overshoot, "%", 2)], ["Gust drift", format(T?.gust.drift, "m", 2)]]) { const c = node("div", undefined, "hero-stat"); c.append(node("span", label), node("b", v)); hero.append(c); }
    root.append(hero); }
  if (T) {
    const grid = node("div", undefined, "domain-grid"), step = section(grid, `Roll step, ${d.stepAngle}°`, "simulated"), alt = section(grid, "Altitude step, 2 m", "simulated"), gust = section(grid, `${d.wind} m/s wind, ${d.gust} m/s gust`, "simulated"), fail = section(grid, "One motor stops at 2 s", T.failure.survives ? "holds attitude" : "not survivable");
    step.append(lineChart({ series: [{ points: series(T.rollStep.trace, 180 / Math.PI) }, { points: T.rollStep.trace.map(s => [s.t, s.r * 180 / Math.PI]), cls: "system" }], x: { label: "time, s", min: 0 }, y: { label: "roll, °", min: 0 } }));
    alt.append(lineChart({ series: [{ points: series(T.altitude.trace) }, { points: T.altitude.trace.map(s => [s.t, s.r]), cls: "system" }], x: { label: "time, s", min: 0 }, y: { label: "height, m", min: 0 } }));
    gust.append(lineChart({ series: [{ points: series(T.gust.trace) }], x: { label: "time, s", min: 0 }, y: { label: "drift, m", min: 0 } }));
    fail.append(lineChart({ series: [{ points: series(T.failure.trace) }], x: { label: "time, s", min: 0 }, y: { label: "tilt, °", min: 0 } }), node("p", T.failure.survives ? "The remaining rotors hold attitude." : `Tilt reaches ${T.failure.tilt.toFixed(0)}° and the aircraft falls ${T.failure.altitudeLoss.toFixed(0)} m in three seconds. A quadcopter cannot survive losing a rotor: use a parachute or six rotors if that matters.`, "small muted"));
    root.append(grid);
  }
  if (x) { const c = section(root, "Rate loop, roll", `crossover ${format(x.roll.bandwidth, "Hz", 3)}`); c.append(lineChart({ series: [{ points: x.roll.curve.map(p => [Math.log10(p.w), p.gain]) }], x: { label: "log₁₀ frequency, rad/s" }, y: { label: "loop gain, dB", min: Math.min(-40, ...x.roll.curve.map(p => p.gain)) } })); }
  for (const [name, pattern] of GROUPS) { const s = section(root, name, "PX4 parameters"), f = node("div", undefined, "field-grid");
    for (const k of Object.keys(PX4_DEFAULTS).filter(k => pattern.test(k))) numberField(f, { label: `${k}${d.params[k] === PX4_DEFAULTS[k] ? "" : " (changed)"}`, id: `controls-${k}`, value: d.params[k], min: 0, max: 1e4, step: "any", op: "controls.gain", args: v => ({ name: k, value: v }), status });
    s.append(f); }
  const src = node("a", "PX4 parameter reference"); src.href = PX4_SOURCE; src.rel = "noopener"; root.append(node("p", "", "small muted")); root.lastChild.append(document.createTextNode("Defaults are from the "), src, document.createTextNode(", which changes between releases."));
  const p = section(root, "Plant", "values linked from other domains override these"), f = node("div", undefined, "field-grid");
  for (const [arg, label, unit, key, port] of PLANT) { const on = port && solution?.applied[`controls.${port}`], i = numberField(f, { label, id: `controls-plant-${arg}`, unit, value: on ? on.value : d.plant[key], step: "any", op: "controls.plant", args: v => ({ [arg]: v }), status }); if (on) i.disabled = true; }
  p.append(f);
  const c = section(root, "Flight-test conditions and limits"), g = node("div", undefined, "field-grid");
  for (const [arg, label, unit, key] of [["loop_rate", "Rate-loop frequency", "Hz", "loopRate"], ["gyro_noise", "Gyro noise (°/s/√Hz)", "", "gyroNoise"], ["gust", "Gust", "m/s", "gust"], ["wind", "Steady wind", "m/s", "wind"], ["step_angle", "Step angle (°)", "", "stepAngle"]]) numberField(g, { label, id: `controls-${arg}`, unit, value: d[key], step: "any", op: "controls.set", args: v => ({ [arg]: v }), status });
  for (const [key, label] of [["phaseMargin", "Phase margin at least (°)"], ["gainMargin", "Gain margin at least (dB)"], ["overshoot", "Overshoot at most (%)"], ["settle", "Settling at most (s)"], ["drift", "Gust drift at most (m)"]]) numberField(g, { label, id: `controls-limit-${key}`, value: d.limits[key], step: "any", op: "controls.limit", args: v => ({ name: key, value: v }), status });
  c.append(g, status);
  if (r) { const k = section(root, "Control checks"); checkCards(k, r.checks || []); }
  const a = section(root, "Assumptions"); assumptionList(a, "controls");
}
onSolved(s => { solution = s; });
live(["controls"], render);
render();
