// Any project as the commands that build it from a blank project, in order.
// Replaying this script and getting the same documents back is the proof that
// a person can make the design from nothing with the tools the system gives:
// a field no command can set cannot be reproduced.
import { DEFAULTS } from "../repository.js";
import { TEMPLATES } from "./parts.js";
import { PX4_DEFAULTS } from "./controls.js";
const has = (v, k) => v[k] !== undefined;
const same = (a, b) => JSON.stringify(a) === JSON.stringify(b);
export function buildScript(documents) {
  const out = [], run = (op, args = {}) => out.push({ op, args }), D = documents, M = D.project, on = new Set(M.domains);
  run("project.rename", { name: M.name }); if (M.summary) run("project.describe", { summary: M.summary });
  for (const domain of M.domains) run("project.add_domain", { domain, start: "empty" });
  if (on.has("mech") && D.mech) {
    for (const [name, value] of Object.entries(D.mech.params ?? {})) if (value !== DEFAULTS[name]) run("mech.set", { name, value });
    for (const f of D.mech.features ?? []) run("mech.add_part", { id: f.id, name: f.name, shape: f.shape, domain: f.domain, x: f.x, y: f.y, z: f.z, width: f.width, depth: f.depth, height: f.height });
    for (const [part, sheet] of Object.entries(D.mech.drawings ?? {})) for (const [field, value] of Object.entries(sheet)) run("mech.drawing", { part, field, value });
  }
  if (on.has("flight")) { const f = D.flight;
    run("flight.set", { name: f.name, rotors: String(f.rotors), prop: f.prop, arm: f.arm, fm: f.fm, motor_efficiency: f.motorEfficiency, esc_efficiency: f.escEfficiency, design_tw: f.designTW, avionics: f.avionicsPower, temperature: f.temperature, pressure: f.pressure, altitude: f.altitude, required_time: f.requiredTime, min_tw: f.minTW, capture: f.capture });
    run("flight.battery", { ...f.battery }); for (const m of f.masses) run("flight.mass", { id: m.id, mass: m.mass, name: m.name, ...(m.part && { part: m.part }), ...(m.source && { source: m.source }) });
    if (f.propulsion) { const p = f.propulsion; run("flight.propulsion", { motor: p.motor, propeller: p.propeller, kv: p.kv, resistance: p.resistance, peak_current: p.peakCurrent, max_power: p.maxPower, rotor_inertia: p.rotorInertia, density: p.density, blockage: p.blockage, ...(p.source && { source: p.source }), table: p.table }); }
  }
  if (on.has("parts")) for (const p of D.parts.parts) {
    const template = TEMPLATES[p.template] ? p.template : "block";
    run("parts.add", { template, id: p.id, name: p.name, x: p.x, y: p.y, z: p.z, material: p.material, rz: p.rz, count: p.count, pattern: p.pattern, stack: p.stack });
    if (p.template === "custom" || !same(p.model, TEMPLATES[template].model)) run("parts.model", { id: p.id, model: p.model, keep_template: p.template !== "custom" });
    for (const [name, value] of Object.entries(p.overrides)) run("parts.param", { id: p.id, name, value });
    if (p.domain !== "Mechanical") run("parts.update", { id: p.id, domain: p.domain });
  }
  if (on.has("pcb")) { const b = D.pcb;
    run("pcb.board", { width: b.width, height: b.height, grid: b.grid, clearance: b.rules.clearance, min_width: b.rules.minWidth, edge: b.rules.edge, density: b.rules.density, ...b.stackup });
    for (const [name, f] of Object.entries(b.footprints ?? {})) run("pcb.footprint", { name, kind: "custom", pads: f.pads, width: f.width, height: f.height, title: f.title, mass: f.mass, pad: f.pad, drill: f.drill, ...(f.smd && { smd: true }), ...(f.padW !== undefined && { pad_w: f.padW, pad_h: f.padH }) });
    for (const p of b.pours ?? []) run("pcb.pour", { net: p.net, layer: p.layer, points: p.points, id: p.id });
    for (const c of b.components) { if (c.footprint === "via") { run("pcb.via", { at: [c.x, c.y], net: c.nets[0], ref: c.id }); continue; } run("pcb.place", { footprint: c.footprint, at: [c.x, c.y], name: c.name, ref: c.id }); if (c.rotation) run("pcb.move", { id: c.id, rotation: String(c.rotation) });
      c.nets.forEach((net, i) => { if (net) run("pcb.net", { pad: `${c.id}:${i + 1}`, net }); }); }
    for (const [net, amps] of Object.entries(b.currents)) run("pcb.current", { net, amps });
    for (const [id, watts] of Object.entries(b.powers)) run("pcb.power", { id, watts });
    for (const net of b.loadNets) run("pcb.load_net", { net });
    for (const r of b.routes) run("pcb.trace", { from: r.from, to: r.to, layer: r.layer, width: r.width, points: r.points, id: r.id });
  }
  if (on.has("silicon")) { const s = D.silicon, c = s.chip, l = D.litho;
    run("silicon.core", { width: s.width, height: s.height, grid: s.grid, margin: s.margin, spacing: s.spacing, target: s.target });
    for (const b of s.blocks) { run("silicon.add_block", { kind: b.kind, id: b.id, name: b.name, at: [b.x, b.y], width: b.width, height: b.height, power: b.power }); if (b.fixed) run("silicon.update_block", { id: b.id, fixed: true }); }
    run("silicon.chip", { frequency: c.frequency, vdd: c.vdd, activity: c.activity, instances: c.instances, theta: c.theta, ambient: c.ambient });
    for (const r of c.rtl) run("silicon.rtl", { block: r.block, source: r.source, ...(r.verilog && { verilog: r.verilog }) });
    if (l) { run("litho.set", { wavelength: l.wavelength, na: l.na, sigma_outer: l.sigmaOuter, sigma_inner: l.sigmaInner, threshold: l.threshold, pixel: l.pixel, source_points: l.sourcePoints, steepness: l.steepness, grid: String(l.n) }); for (const s of l.shapes) run("litho.shape", { ...s }); }
  }
  if (on.has("optics")) { const o = D.optics;
    run("optics.set", { wavelength: o.wavelength, power: o.source.power, diameter: o.source.diameter, target_z: o.target.z, target_diameter: o.target.diameter, target_name: o.target.name, rays: o.rays, pixel: o.pixel, required_gsd: o.requiredGsd, object_distance: o.objectDistance, fields: o.fields.join(",") });
    for (const e of o.elements) run("optics.add", { kind: e.kind, id: e.id, name: e.name, z: e.z, diameter: e.diameter, ...(has(e, "thickness") && { thickness: e.thickness }), ...(has(e, "glass") && { glass: e.glass }), ...(has(e, "r1") && { r1: e.r1, r2: e.r2 }) });
  }
  if (on.has("fluids")) { const f = D.fluids;
    run("fluids.fan", { preset: f.fan.preset, ...(f.fan.preset === "custom" && { shutoff: f.fan.shutoff, free: f.fan.free }) }); run("fluids.set", { temperature: f.temperature, pressure: f.pressure, heat_limit: f.heatLimit, heat: f.heatLoad });
    f.elements.forEach((e, i) => { const { kind, id, ...rest } = e; run("fluids.add", { kind, id, ...rest, ...(i && { after: f.elements[i - 1].id }) }); });
  }
  if (on.has("architecture")) { const a = D.architecture;
    run("architecture.site", { name: a.site.name, outdoor: a.site.outdoor, solar: a.site.solar, glazing_u: a.glazing.U, shgc: a.glazing.shgc, units: a.units, limit: a.limit, equipment: a.equipmentPower });
    for (const r of a.rooms) { run("architecture.add_room", { label: r.name, id: r.id, at: [r.x, r.y], width: r.width, depth: r.depth, height: r.height, people: r.people, lighting: r.lighting, supply: r.supply, supply_temp: r.supplyTemp, setpoint: r.setpoint, roof_u: r.roofU, mode: r.mode, ...(r.equipment && { equipment: true }) });
      for (const [side, w] of Object.entries(r.walls)) run("architecture.wall", { id: r.id, side, exterior: w.exterior, u: w.U, window: w.window }); }
  }
  if (on.has("equations")) { const e = D.equations;
    for (const x of e.inputs) run("equations.input", { ...x }); for (const x of e.outputs) run("equations.output", { ...x });
  }
  if (on.has("power")) { const s = D.power, m = s.mosfet;
    run("power.set", { name: s.name, pwm: s.pwm, deadtime: s.deadtime, gate_drive: s.gateDrive, gate_resistor: s.gateResistor, driver_peak: s.driverPeak, ringing: s.ringing, efficiency_limit: s.efficiencyLimit, phase: s.phase, bus: s.bus, duty: s.duty, junction: s.junction, peak_bus: s.peakBus });
    run("power.mosfet", { part: m.part, vds: m.vds, rds25: m.rds25, rds_tc: m.rdsTc, qg: m.qg, qsw: m.qsw, qoss: m.qoss, qrr: m.qrr, vplateau: m.vplateau, rg: m.rg, vsd: m.vsd, rth_jc: m.rthJc, tj_max: m.tjMax, ...(m.source && { source: m.source }) });
    run("power.capacitor", { ...s.capacitor });
  }
  if (on.has("structures")) { const s = D.structures; run("structures.set", { name: s.name, od: s.od, wall: s.wall, length: s.length, modulus: s.modulus, density: s.density, strength: s.strength, safety_factor: s.safetyFactor, deflection_limit: s.deflectionLimit, separation: s.separation, blades: s.blades, tip_mass: s.tipMass, load_thrust: s.loadThrust, load_torque: s.loadTorque, ...(s.source && { source: s.source }) }); }
  if (on.has("controls")) { const c = D.controls, p = c.plant;
    run("controls.set", { loop_rate: c.loopRate, gyro_noise: c.gyroNoise, gust: c.gust, wind: c.wind, step_angle: c.stepAngle });
    run("controls.plant", { mass: p.mass, ixx: p.ixx, iyy: p.iyy, izz: p.izz, arm: p.arm, rotors: p.rotors, hover_throttle: p.hoverThrottle, max_thrust: p.maxThrust, tau: p.tau, kappa: p.kappa, drag_area: p.dragArea });
    for (const [name, value] of Object.entries(c.params)) if (value !== PX4_DEFAULTS[name]) run("controls.gain", { name, value });
    for (const [name, value] of Object.entries(c.limits)) run("controls.limit", { name, value });
  }
  if (on.has("testing")) { const t = D.testing;
    for (const s of t.datasets) run("testing.dataset", { ...s });
    for (const s of t.traces ?? []) run("testing.trace", { ...s });
    for (const { result, covers, ...x } of t.plan) { run("testing.add", { ...x, ...(covers && { covers: covers.join(",") }) }); if (result) { const { digest, ...r } = result; run("testing.record", { id: x.id, ...r }); } }
  }
  for (const l of D.system.links) { run("system.link", { from: l.from, to: l.to, id: l.id, ...(l.note && { note: l.note }) }); if (!l.enabled) run("system.toggle", { id: l.id, enabled: false }); }
  for (const id of Object.keys(D.system.confirmed)) run("system.confirm", { id });
  for (const h of M.headline) run("project.pin", { port: h.port, ...(h.label && { label: h.label }) });
  for (const q of M.requirements) run("project.require", { ...q });
  for (const s of M.studies) run("project.study", { ...s });
  if (M.optimize) { const o = M.optimize;
    run("project.optimize", { title: o.title, objective: o.objective.port, sense: o.objective.sense, ...(o.summary && { summary: o.summary }), trials: o.trials });
    for (const v of o.variables) run("project.optimize_variable", { ...v }); for (const c of o.constraints) run("project.optimize_constraint", { ...c });
  }
  for (const t of M.tour) run("project.tour_step", { studio: t.studio, title: t.title, text: t.text, ...(t.command && { command: t.command }), ...(t.watch && { watch: t.watch.join(",") }), ...(t.checks && { checks: t.checks.join(",") }) });
  for (const [domain, p] of Object.entries(M.scene.placements)) { const [x, y, z] = p.at; run("project.place", { domain, x, y, z, ...(has(p, "rx") && { rx: p.rx }), ...(has(p, "rz") && { rz: p.rz }), ...(has(p, "phase") && { phase: p.phase }), ...(has(p, "count") && { count: String(p.count) }) }); }
  if (M.scene.site) run("project.site", { ...M.scene.site });
  for (const c of M.credits) run("project.credit", { ...c });
  return out;
}
