// The design environment: every domain's documents, typed ports, evaluators,
// live certificates, assumptions and commands. What is being designed comes
// from a project file (core/project.js); the browser and the headless agent
// both build the environment here, with their own model transport and
// mechanical host, so every surface runs the same commands and physics.
import { CommandBus } from "../core/bus.js";
import { Store } from "../core/store.js";
import { System, SYSTEM_DEFAULT, validateSystem } from "../core/system.js";
import { DOMAINS, STUDIOS, projectDefault, validateProject, readProjectFile, projectFile, getPath, setPath } from "../core/project.js";
import { FOOTPRINTS, footprintOf, footprintPads, pcbExample, siliconExample, documentFor, validateDocument, pcbChecks, siliconChecks, connectPads, moveComponent, clone, BLOCK_POWER } from "../studio-model.js";
import { kernelBoard, densityIssues, autoroute, pcbAssumptions, netOf, powerPort, pourLoss, pourSize } from "./pcb.js";
import { powerSummary, siliconAssumptions } from "./silicon.js";
import { lithoExample, validateLitho } from "./litho.js";
import { opticsExample, validateOptics, opticsResults, opticsAssumptions, paraxial, lensmaker, bestFocus, GLASSES, refractiveIndex, KINDS as OPTIC_KINDS } from "./optics.js";
import { fluidsExample, validateFluids, fluidsResults, fluidsAssumptions, FANS, fanPressure, KINDS as FLUID_KINDS } from "./fluids.js";
import { architectureExample, validateArchitecture, architectureResults, architectureAssumptions, overlaps, SIDES, AIR } from "./architecture.js";
import { partsExample, validateParts, partsAssumptions, TEMPLATES, MATERIALS, primitive, combine, modify, FEATURE_SHAPES, COMBINE, MODIFY } from "./parts.js";
import { equationsExample, validateEquations, equationsResults, equationsAssumptions } from "./equations.js";
import { powerExample, validatePower, stageLoss, busRipple, powerAssumptions, MOSFET_FIELDS, CAPACITOR_FIELDS, CAPACITOR_DERATE } from "./power.js";
import { structuresExample, validateStructures, armResults, structuresAssumptions, STRUCTURE_FIELDS } from "./structures.js";
import { controlsExample, validateControls, plantOf, rateLoop, flightTests, tune, controlsAssumptions, PX4_DEFAULTS } from "./controls.js";
import { testingExample, validateTesting, evaluatePlan, digest, KINDS as TEST_KINDS } from "./testing.js";
import { mechanicalCommands, mechanicalAdapter, mechanicalAssumptions } from "./mechanical.js";
import { flightExample, validateFlight, flightResults, flightAssumptions, validatePropulsion } from "./flight.js";
import { instances } from "./parts.js";
import { elementMass } from "./optics.js";
import { FOOTPRINTS as FP } from "../studio-model.js";
import { values } from "../repository.js";
// Illustrative component masses for the board budget, grams.
const FOOTPRINT_MASS = { connector: 1.5, resistor: .3, capacitor: 1, transistor: 2, dip8: .55 };

// Each domain's studio tab and what it adds to a project.
export const DOMAIN_INFO = {
  system: { title: "System", icon: "◎", color: "var(--text)" },
  flight: { title: "Flight", icon: "✈", color: "var(--ember)", studio: "flight", blurb: "Multirotor mass budget, rotor power by momentum theory, battery and endurance." },
  mech: { title: "Mechanical", icon: "⬡", color: "var(--cyan)", studio: "mechanical", blurb: "Power-stage heat sink on the Rust electrothermal graph, assembly, STL and shop drawings." },
  parts: { title: "Parts", icon: "⚙", color: "var(--cyan)", studio: "parts", blurb: "Solid feature trees on the Rust field kernel: gears, sheet metal, pipes, lattices, wings." },
  pcb: { title: "PCB", icon: "▦", color: "var(--ember)", studio: "pcb", blurb: "Placement, obstacle-aware routing, kernel trace resistance and impedance, Gerber and KiCad." },
  silicon: { title: "Silicon", icon: "◈", color: "var(--violet)", studio: "silicon", blurb: "Verilog synthesis, equivalence, nine-corner timing and power, floorplan, GDSII and lithography." },
  optics: { title: "Optics", icon: "◐", color: "var(--violet)", studio: "optics", blurb: "Real 3D ray traces through catalogue glass, spots, focus and image scale." },
  fluids: { title: "Fluids", icon: "≋", color: "var(--cyan)", studio: "fluids", blurb: "Air paths of ducts, fittings, filters and fin channels against a fan or pressure source." },
  architecture: { title: "Architecture", icon: "⌂", color: "var(--green)", studio: "architecture", blurb: "Rooms on a plan with envelopes, glazing, people, equipment and supply air in heat balance." },
  structures: { title: "Structures", icon: "⌖", color: "var(--cyan)", studio: "structures", blurb: "Arm stress, deflection and bending mode against the rotor's once-per-rev and blade-pass." },
  power: { title: "Power", icon: "⚡", color: "var(--ember)", studio: "power", blurb: "A three-phase inverter from its MOSFET's datasheet: conduction, switching, gate and dead-time loss." },
  controls: { title: "Controls", icon: "◉", color: "var(--violet)", studio: "controls", blurb: "PX4's cascade flown in a 6-DoF simulation, with loop margins, step responses, gusts and motor failure." },
  testing: { title: "Testing", icon: "✓", color: "var(--green)", studio: "testing", blurb: "A test plan that separates simulation, correlation against measurements, and physical tests." },
  equations: { title: "Equations", icon: "ƒ", color: "var(--green)", studio: "equations", blurb: "Your own relations: named inputs and formula outputs that link to any other domain." },
  project: { title: "Requirements", icon: "☑", color: "var(--text)" },
};
// Each domain with nothing in it. A starter is an example to edit; empty is a
// design that has not been started.
const EMPTY = {
  pcb: () => ({ ...pcbExample(), components: [], routes: [], currents: {}, powers: {}, loadNets: [] }),
  silicon: () => { const s = siliconExample(); return { ...s, blocks: [], chip: { ...s.chip, rtl: [] } }; },
  litho: () => ({ ...lithoExample(), shapes: [] }),
  optics: () => ({ ...opticsExample(), elements: [], requiredGsd: 0, target: { name: "Image plane", z: 100, diameter: 25 } }),
  fluids: () => ({ ...fluidsExample(), elements: [] }),
  architecture: () => ({ ...architectureExample(), rooms: [], units: 1 }),
  parts: () => ({ parts: [] }),
  flight: () => ({ ...flightExample(), name: "Vehicle", masses: [] }),
  equations: equationsExample, power: powerExample, structures: structuresExample, controls: controlsExample, testing: testingExample,
};
const ID = /^[a-zA-Z][a-zA-Z0-9_-]{0,23}$/;
const partPort = (p, k) => `${p.id}_${k}`.replace(/[^a-zA-Z0-9_]/g, "_");
export { DOMAINS };
const studio = kind => ({ example: () => kind === "pcb" ? pcbExample() : siliconExample(), validate: d => validateDocument({ format: "emergence-studio/1", kind, design: d }).design });
const finiteIn = (v, lo, hi) => Number.isFinite(v) && v >= lo && v <= hi;
const check = (id, text, ok, reason, value, unit) => ({ id, text, status: ok === null ? "undecided" : ok ? "met" : "unmet", ...(reason && { reason }), ...(value !== undefined && { value: { nominal: value }, unit }) });

export const PART_BUDGET = 2.5e5;
export function createProject({ request, mech, loadCorners, verilog = "", mesh }) {
  // The browser passes its cached binary mesher; headless callers mesh directly.
  mesh ||= (model, overrides) => request("model", { input: { model, overrides, budget: PART_BUDGET } }).then(r => r.result);
  const store = new Store({
    pcb: studio("pcb"), silicon: studio("silicon"), litho: { example: lithoExample, validate: validateLitho },
    optics: { example: opticsExample, validate: validateOptics }, fluids: { example: fluidsExample, validate: validateFluids },
    architecture: { example: architectureExample, validate: validateArchitecture }, parts: { example: partsExample, validate: validateParts },
    flight: { example: flightExample, validate: validateFlight },
    system: { example: SYSTEM_DEFAULT, validate: d => validateSystem(d) },
    project: { example: projectDefault, validate: validateProject },
    equations: { example: equationsExample, validate: validateEquations },
    power: { example: powerExample, validate: validatePower },
    structures: { example: structuresExample, validate: validateStructures },
    controls: { example: controlsExample, validate: validateControls },
    testing: { example: testingExample, validate: validateTesting },
  });
  const bus = new CommandBus();
  for (const id of store.ids()) bus.domain(id, { snapshot: () => store.snapshot(id), restore: s => store.restore(id, s) });
  bus.domain("mech", { snapshot: () => mech.snapshot(), restore: s => mech.restore(s) });
  const m = (id, fn) => store.mutate(id, fn);

  // Chip corner libraries are sent to the model worker once per Worker: a
  // restarted Worker has lost them, and a failed load is retried.
  let corners = null;
  const chip = async (d, rtl, temperature) => {
    if (corners && (await corners.catch(() => null))?.generation !== request.generation) corners = null;
    corners ||= (async () => { const data = await loadCorners(); await request("chip", { input: { libs: data.corners, verilog: "module t(input a, output y); assign y = a; endmodule" } }); return { data, generation: request.generation }; })();
    await corners;
    return (await request("chip", { input: { verilog: rtl.source === "custom" ? rtl.verilog : RTL[rtl.source], frequency: d.chip.frequency * 1e6, vdd: d.chip.vdd, activity: d.chip.activity / 100, instances: d.chip.instances, temperature, samples: 3000 } })).result;
  };
  const RTL = typeof verilog === "string" ? { mini8: verilog } : { ...verilog };
  const adapters = new Map([
    ["mech", mechanicalAdapter(mech, request)],
    ["silicon", { title: "Silicon", ports: () => ({
        inputs: [{ id: "ambient", label: "Package ambient", unit: "degC", min: -40, max: 125 }, { id: "frequency", label: "Clock frequency", unit: "", min: .001, max: 5000 }],
        outputs: [{ id: "power", label: "Die power", unit: "W" }, { id: "junction", label: "Die junction temperature", unit: "degC" }, { id: "fmax", label: "RTL maximum clock (MHz)", unit: "" },
          { id: "utilization", label: "Requested block area", unit: "%" }, { id: "peak_density", label: "Peak block power density", unit: "W/mm2" }] }),
      async evaluate(d, ov) {
        d = { ...d, chip: { ...d.chip, ...(ov.ambient !== undefined && { ambient: ov.ambient }), ...(ov.frequency !== undefined && { frequency: ov.frequency }) } };
        const geometry = siliconChecks(d), checks = [check("floorplan", "Floorplan placement checks", !geometry.issues.length, geometry.issues.map(i => i.message).join(" "))];
        let flows = {}, iterations = 0, error = null;
        if (d.chip.rtl.length) {
          // Die electrothermal loop: Tj = Ta + θ·P(Tj), leakage rising with Tj.
          let tj = d.chip.ambient + 273.15;
          try {
            for (iterations = 1; iterations <= 8; iterations++) {
              if (tj < 250 || tj > 425) throw new Error(`Junction ${(tj - 273.15).toFixed(1)} °C is outside the 250–425 K characterised corners.`);
              flows = {}; for (const r of d.chip.rtl) flows[r.block] = { ...(await chip(d, r, tj)), source: r.source };
              const next = d.chip.ambient + 273.15 + d.chip.theta * powerSummary(d, flows).total;
              if (Math.abs(next - tj) < 1e-3) { tj = next; break; } tj = next;
            }
          } catch (e) { error = e.message; }
          for (const [block, f] of Object.entries(flows)) {
            const name = d.blocks.find(b => b.id === block)?.name || block;
            if (f.equivalence) checks.push(check(`equiv-${block}`, `${name}: netlist equals its RTL`, !f.equivalence.mismatch, f.equivalence.mismatch ?? `${f.equivalence.points.toLocaleString()} (input, state) points${f.equivalence.exhaustive ? ", exhaustive: a proof" : ", sampled"}`));
            if (f.timing) checks.push(check(`setup-${block}`, `${name}: setup slack ≥ 0 at ${d.chip.frequency} MHz`, f.timing.setup_slack >= 0, null, f.timing.setup_slack * 1e12, "ps"), check(`hold-${block}`, `${name}: hold slack ≥ 0`, f.timing.hold_slack >= 0, null, f.timing.hold_slack * 1e12, "ps"));
          }
          if (error) checks.push({ id: "chip", text: "Chip flow", status: "unavailable", reason: error });
        }
        const flow = Object.values(flows)[0] ?? null, fmax = Math.min(...Object.values(flows).filter(f => f.timing).map(f => f.timing.fmax));
        const p = powerSummary(d, flows), junction = d.chip.ambient + d.chip.theta * p.total;
        checks.push(check("tj", "Die junction ≤ 125 °C", junction <= 125, null, junction, "degC"));
        return { outputs: { power: { value: p.total, unit: "W" }, junction: { value: junction, unit: "degC" }, fmax: { value: Number.isFinite(fmax) ? fmax / 1e6 : NaN, unit: "" },
          utilization: { value: geometry.utilization, unit: "%" }, peak_density: { value: p.peak.density, unit: "W/mm2" } },
          checks, details: { flow, flows, power: p, geometry, iterations }, certificates: Object.entries(flows).filter(([, f]) => f.equivalence).map(([b, f]) => ({ label: `${b}: netlist ≡ RTL on ${f.equivalence.points.toLocaleString()} points`, ok: !f.equivalence.mismatch })) };
      } }],
    ["pcb", { title: "PCB", ports: d => ({
        inputs: [...d.components.map(c => ({ id: powerPort(c.id), label: `${c.id} dissipation`, unit: "W", min: 0, max: 500 })), { id: "load_current", label: "Load net RMS current", unit: "A", min: 0, max: 500 }],
        outputs: [{ id: "copper_loss", label: "Copper I²R loss", unit: "W" }, { id: "dissipation", label: "Board + component dissipation", unit: "W" }, { id: "unrouted", label: "Unrouted connections", unit: "count" },
          { id: "board_mass", label: "Populated board mass", unit: "kg" }] }),
      async evaluate(d, ov) {
        d = clone(d);
        for (const c of d.components) if (ov[powerPort(c.id)] !== undefined) d.powers[c.id] = ov[powerPort(c.id)];
        if (ov.load_current !== undefined) for (const n of d.loadNets) d.currents[n] = ov.load_current;
        const geo = pcbChecks(d), density = densityIssues(d), k = kernelBoard(d);
        let kernel = null, error = null;
        if (d.routes.length || d.pours.length) try { kernel = (await request("pcb", { input: { board: k.board, currents: k.currents, temperature: 40 } })).result; } catch (e) { error = e.message; }
        const copper = (kernel?.copper_loss ?? 0) + pourLoss(d), components = Object.values(d.powers).reduce((s, w) => s + w, 0);
        const checks = [check("geometry", "Copper geometry checks", !geo.issues.length, geo.issues.slice(0, 3).map(i => i.message).join(" ") + (geo.issues.length > 3 ? ` (+${geo.issues.length - 3} more)` : "")),
          check("density", `Track current density ≤ ${d.rules.density} A/mm²`, !density.length, density.map(i => i.message).join(" ")),
          check("routed", "Every net connected", geo.unrouted === 0 ? true : null, geo.unrouted ? `${geo.unrouted} connection${geo.unrouted > 1 ? "s" : ""} remain` : "")];
        if (error) checks.push({ id: "kernel", text: "Kernel board", status: "unavailable", reason: error });
        // FR-4 at 1.85 g/cm³ over the outline, drawn copper at 8.96 g/cm³, components by footprint.
        const t = d.stackup.copper * .0348, copperArea = d.routes.reduce((s, r) => s + r.points.slice(1).reduce((a, p, i) => a + Math.hypot(p[0] - r.points[i][0], p[1] - r.points[i][1]) * r.width, 0), 0) + d.pours.reduce((s, p) => s + pourSize(p).area, 0);
        const mass = (d.width * d.height * d.stackup.thickness * 1.85e-3 + copperArea * t * 8.96e-3 + d.components.reduce((s, c) => s + (footprintOf(d, c.footprint).mass ?? FOOTPRINT_MASS[c.footprint]), 0)) / 1000;
        return { outputs: { copper_loss: { value: copper, unit: "W" }, dissipation: { value: copper + components, unit: "W" }, unrouted: { value: geo.unrouted, unit: "count" }, board_mass: { value: mass, unit: "kg" } },
          checks, details: { kernel, geometry: geo, density, board: d },
          certificates: kernel ? [{ label: "Kernel trace R = ρL/(wt) at 40 °C", ok: kernel.traces.every(t => Math.abs(t.resistance / (1.678e-8 * (1 + 3.93e-3 * 20) * t.length / (t.width * d.stackup.copper * 34.8e-6)) - 1) < 1e-9) }] : [] };
      } }],
    ["optics", { title: "Optics", ports: () => ({
        inputs: [{ id: "source_power", label: "Source power", unit: "W", min: 0, max: 100 }, { id: "altitude", label: "Object distance (altitude)", unit: "m", min: 1, max: 5000 }],
        outputs: [{ id: "power_on_target", label: "Power on target", unit: "W" }, { id: "efl", label: "System focal length", unit: "mm" }, { id: "spot_rms", label: "RMS spot radius", unit: "um" },
          { id: "transmission", label: "Transmission to target", unit: "%" }, { id: "irradiance", label: "Mean irradiance on target", unit: "W/mm2" },
          { id: "mass", label: "Lens glass mass", unit: "kg" }, { id: "gsd", label: "Ground sample distance", unit: "cm" }, { id: "fnumber", label: "Working f-number", unit: "" }] }),
      async evaluate(d, ov) {
        if (ov.source_power !== undefined) d = { ...d, source: { ...d.source, power: ov.source_power } };
        const altitude = ov.altitude ?? d.objectDistance, r = opticsResults(d, altitude), lens = d.elements.find(e => e.kind === "lens");
        const checks = !d.elements.length ? [] : [check("geometry", "Elements clear each other and the target", !r.issues.length, r.issues.map(i => i.message).join(" ")),
          check("capture", `Every field images onto the ${d.target.name.toLowerCase()}`, r.fields.every(f => f.transmitted > 0 && Math.hypot(...f.centroid) <= d.target.diameter / 2), r.fields.map(f => `${f.field}°: ${(f.rms * 1000).toFixed(1)} µm RMS`).join(" · ")),
          ...(d.requiredGsd > 0 ? [check("gsd", `Ground sample distance ≤ ${d.requiredGsd} cm at ${altitude} m`, r.gsd <= d.requiredGsd, null, r.gsd, "cm")] : [])];
        const certs = [{ label: "Paraxial EFL of a lone lens equals the thick-lens lensmaker equation", ok: !lens || d.elements.filter(e => e.kind === "lens").length > 1 || Math.abs(paraxial({ ...d, elements: [lens] }).efl / lensmaker(lens, d.wavelength) - 1) < 1e-10 }];
        return { outputs: { power_on_target: { value: r.powerOnTarget, unit: "W" }, efl: { value: r.paraxial.efl, unit: "mm" }, spot_rms: { value: r.trace.rms * 1000, unit: "um" },
          transmission: { value: 100 * r.trace.onTarget, unit: "%" }, irradiance: { value: r.irradiance, unit: "W/mm2" }, mass: { value: r.mass / 1000, unit: "kg" },
          gsd: { value: r.gsd, unit: "cm" }, fnumber: { value: r.fnumber, unit: "" } }, checks, details: { ...r, altitude }, certificates: certs };
      } }],
    ["fluids", { title: "Fluids", ports: () => ({
        inputs: [{ id: "air_temperature", label: "Inlet air temperature", unit: "degC", min: -40, max: 150 }, { id: "sink_gap", label: "Sink channel gap", unit: "mm", min: .05, max: 200 },
          { id: "sink_height", label: "Sink channel height", unit: "mm", min: .5, max: 2000 }, { id: "sink_length", label: "Sink channel length", unit: "mm", min: .1, max: 5000 },
          { id: "sink_width", label: "Sink approach width", unit: "mm", min: .5, max: 2000 }, { id: "fin_count", label: "Fin count", unit: "count", min: 2, max: 200 }, { id: "heat_load", label: "Heat carried by the air", unit: "W", min: 0, max: 1e5 },
          { id: "drive_pressure", label: "Wake driving pressure", unit: "Pa", min: 0, max: 5000 }],
        outputs: [{ id: "flow", label: "Volume flow", unit: "L/s" }, { id: "duct_velocity", label: "Sink approach velocity", unit: "m/s" }, { id: "sink_dp", label: "Sink pressure drop", unit: "Pa" },
          { id: "total_dp", label: "Path pressure drop", unit: "Pa" }, { id: "exhaust_rise", label: "Air temperature rise", unit: "K" }, { id: "channel_re", label: "Channel Reynolds number", unit: "" }] }),
      async evaluate(d, ov) {
        if (d.fan.preset === "wake" && ov.drive_pressure !== undefined) d = { ...d, fan: { ...d.fan, shutoff: ov.drive_pressure } };
        const r = fluidsResults(d, ov), op = r.op;
        if (op.error) return { outputs: {}, checks: [check("fan", "Fan overcomes the path", false, op.error)], details: r };
        const sink = op.elements.find(e => e.kind === "channels"), rho = op.props.rho, rise = r.heat !== undefined ? r.heat / (rho * op.Q * op.props.cp) : NaN;
        const transition = op.elements.filter(e => e.regime === "transition");
        const checks = [check("fan", "Fan overcomes the path", true, `${(op.Q * 1000).toFixed(2)} L/s at ${op.pressure.toFixed(1)} Pa`),
          check("regime", "Every passage in a defined flow regime", transition.length ? null : true, transition.length ? `${transition.map(e => e.name).join(", ")} in transition (Re 2300–4000): flow bounded ${(op.low * 1000).toFixed(2)}–${(op.high * 1000).toFixed(2)} L/s` : ""),
          check("rise", `Air temperature rise ≤ ${d.heatLimit} K`, Number.isFinite(rise) ? rise <= d.heatLimit : null, null, rise, "K")];
        const residual = Math.abs(fanPressure(d.fan, op.Q) - op.elements.reduce((s, e) => s + (e.lo + e.hi) / 2, 0)) / Math.max(1, op.pressure);
        return { outputs: { flow: { value: op.Q * 1000, lo: op.low * 1000, hi: op.high * 1000, unit: "L/s" }, duct_velocity: { value: sink ? sink.approach : NaN, unit: "m/s" },
          sink_dp: { value: sink ? (sink.lo + sink.hi) / 2 : NaN, unit: "Pa" }, total_dp: { value: op.pressure, unit: "Pa" }, exhaust_rise: { value: rise, unit: "K" }, channel_re: { value: sink?.re ?? NaN, unit: "" } },
          checks, details: r, certificates: [{ label: "Fan pressure balances the path at the operating point", ok: transition.length > 0 || residual < 1e-6 }] };
      } }],
    ["architecture", { title: "Architecture", ports: () => ({
        inputs: [{ id: "equipment_power", label: "Dissipation per equipment unit", unit: "W", min: 0, max: 1e5, aggregate: "sum" }],
        outputs: [{ id: "room_temperature", label: "Equipment room air temperature", unit: "degC" }, { id: "cooling_load", label: "Building cooling load", unit: "W" },
          { id: "equipment_load", label: "Equipment heat in the room", unit: "W" }, { id: "floor_area", label: "Floor area", unit: "m2" }, { id: "outdoor_temperature", label: "Site design air temperature", unit: "degC" }] }),
      async evaluate(d, ov) {
        const r = architectureResults(d, ov), room = r.equipmentRoom, clash = overlaps(d);
        const checks = [...(room ? [check("room", `Equipment room ≤ ${d.limit} °C`, room.temperature <= d.limit, null, room.temperature, "degC")] : []),
          ...r.rooms.map((x, i) => check(`vent-${x.id}`, `${x.name}: outdoor air ≥ ASHRAE 62.1 office rate`, d.rooms[i].supply / 3.6 >= x.ventilation, `${(d.rooms[i].supply / 3.6).toFixed(0)} L/s supplied, ${x.ventilation.toFixed(0)} L/s required (if supply is outdoor air)`)),
          ...(!d.rooms.length ? [] : [check("plan", "Rooms do not overlap", !clash.length, clash.map(([a, b]) => `${a.name} overlaps ${b.name}`).join("; "))])];
        const balance = d.rooms.map((x, i) => { const z = r.rooms[i]; return x.mode === "supplied" ? Math.abs(AIR.rho * AIR.cp * x.supply / 3600 * (z.temperature - x.supplyTemp) - (z.internal + z.solar + z.UA * (d.site.outdoor - z.temperature))) : 0; });
        return { outputs: { room_temperature: { value: room?.temperature ?? NaN, unit: "degC" }, cooling_load: { value: r.cooling, unit: "W" }, equipment_load: { value: r.equipmentW, unit: "W" }, floor_area: { value: r.floor, unit: "m2" }, outdoor_temperature: { value: d.site.outdoor, unit: "degC" } },
          checks, details: r, certificates: [{ label: "Supplied-room heat balances close", ok: balance.every(b => b < 1e-9) }] };
      } }],
    ["flight", { title: "Flight", ports: () => ({
        inputs: [{ id: "air_temperature", label: "Air temperature", unit: "degC", min: -40, max: 55 }, { id: "esc_mass", label: "Mass per ESC", unit: "kg", min: 0, max: 50, aggregate: "sum" },
          { id: "camera_mass", label: "Camera optics mass", unit: "kg", min: 0, max: 50, aggregate: "sum" }, { id: "parts_mass", label: "Airframe parts mass", unit: "kg", min: 0, max: 500, aggregate: "sum" }],
        outputs: [{ id: "phase_current", label: "Design-point motor phase current", unit: "A" }, { id: "bus_voltage", label: "Bus voltage at the design point", unit: "V" }, { id: "peak_voltage", label: "Full-charge pack voltage", unit: "V" }, { id: "design_duty", label: "Design-point duty", unit: "%" },
          { id: "hover_throttle", label: "Hover throttle", unit: "%" }, { id: "hover_rpm", label: "Hover rotor speed", unit: "rpm" }, { id: "max_rpm", label: "Full-throttle rotor speed", unit: "rpm" }, { id: "max_thrust", label: "Full-throttle thrust per rotor", unit: "N" }, { id: "max_torque", label: "Full-throttle rotor torque", unit: "Nm" },
          { id: "tw_full", label: "Thrust/weight at full charge", unit: "" }, { id: "tw_end", label: "Thrust/weight at end of discharge", unit: "" }, { id: "motor_tau", label: "Motor time constant", unit: "s" }, { id: "torque_ratio", label: "Rotor torque per thrust", unit: "m" }, { id: "implied_fm", label: "Figure of merit the data implies", unit: "" },
          { id: "ixx", label: "Roll inertia", unit: "kgm2" }, { id: "iyy", label: "Pitch inertia", unit: "kgm2" }, { id: "izz", label: "Yaw inertia", unit: "kgm2" }, { id: "arm_length", label: "Rotor arm from centre", unit: "mm" }, { id: "rotor_count", label: "Rotors", unit: "count" }, { id: "tip_mass", label: "Mass at the arm tip", unit: "kg" },
          { id: "total_mass", label: "Take-off mass", unit: "kg" }, { id: "hover_power", label: "Hover electrical power", unit: "W" }, { id: "esc_current", label: "Design-point ESC current", unit: "A" },
          { id: "hover_current", label: "Hover ESC current", unit: "A" }, { id: "flight_time", label: "Hover endurance", unit: "min" }, { id: "disk_loading", label: "Disk loading", unit: "Pa" },
          { id: "downwash", label: "Induced velocity at hover", unit: "m/s" }, { id: "wake_pressure", label: "Captured wake pressure", unit: "Pa" }, { id: "altitude", label: "Mission altitude", unit: "m" },
          { id: "charger_heat", label: "Charger heat per pack (1 h charge, 92%)", unit: "W" }] }),
      async evaluate(d, ov) {
        const r = flightResults(d, ov), h = r.hover, x = r.design, p = d.propulsion;
        const energy = d.battery.capacity * r.V, charger = energy * (1 / .92 - 1), tip = d.masses.filter(m => /motor|prop/i.test(`${m.id} ${m.name}`)).reduce((s, m) => s + m.mass, 0) / d.rotors;
        const checks = [check("endurance", `Hover endurance ≥ ${d.requiredTime} min`, r.time >= d.requiredTime, null, r.time, "min"),
          ...(p ? [check("hover-margin", "Hover throttle ≤ 70% (headroom to manoeuvre)", h.feasible && h.throttle <= .7, null, 100 * h.throttle, "%"),
            check("tw", `Thrust/weight at end of discharge ≥ ${d.minTW}`, r.twEnd >= d.minTW, `${r.twFull.toFixed(2)} at full charge, ${r.twEnd.toFixed(2)} at ${r.Vend.toFixed(1)} V`, r.twEnd, ""),
            check("design-point", `Rotors can reach the design thrust-to-weight ${d.designTW}`, x.feasible, x.feasible ? `${(100 * x.throttle).toFixed(0)}% throttle, ${x.rpm.toFixed(0)} rpm` : "Beyond the top of the bench table"),
            check("motor-limit", `Motor within ${p.peakCurrent} A and ${p.maxPower} W at full throttle`, p.table.at(-1).amps <= p.peakCurrent && p.table.at(-1).watts <= p.maxPower, `${p.table.at(-1).amps} A, ${p.table.at(-1).watts} W on the bench`),
            check("clearance", "Rotors clear each other by at least 25 mm", r.clearance >= 25, null, r.clearance, "mm")] : [check("tw", `Design thrust-to-weight ${d.designTW} ≥ ${d.minTW} for control authority`, d.designTW >= d.minTW)]),
          check("pack", "Pack can deliver the design-point power", Number.isFinite(x.I), Number.isFinite(x.I) ? `${x.I.toFixed(1)} A at ${x.Vt.toFixed(1)} V` : "Power exceeds V²/4R: the pack voltage collapses"),
          check("crate", "Design current ≤ 10C continuous", x.I <= 10 * d.battery.capacity, null, x.I / d.battery.capacity, "")];
        const kt = p ? 60 / (2 * Math.PI * p.kv) : NaN, omega = h.rpm * Math.PI / 30, tau = p ? p.rotorInertia / (kt * kt / Math.max(1e-6, p.resistance) + 2 * h.torque / omega) : NaN;
        return { outputs: { phase_current: { value: x.phase, unit: "A" }, bus_voltage: { value: x.Vt, unit: "V" }, peak_voltage: { value: d.battery.cells * 4.2, unit: "V" }, design_duty: { value: 100 * x.duty, unit: "%" }, hover_throttle: { value: 100 * h.throttle, unit: "%" }, hover_rpm: { value: h.rpm, unit: "rpm" }, max_rpm: { value: p ? p.table.at(-1).rpm : NaN, unit: "rpm" },
          max_thrust: { value: r.thrustFull, unit: "N" }, max_torque: { value: p ? p.table.at(-1).torque : NaN, unit: "Nm" }, tw_full: { value: r.twFull, unit: "" }, tw_end: { value: r.twEnd, unit: "" }, motor_tau: { value: tau, unit: "s" }, torque_ratio: { value: h.torque / h.T, unit: "m" }, implied_fm: { value: h.fm, unit: "" },
          ixx: { value: r.Ixx, unit: "kgm2" }, iyy: { value: r.Ixx, unit: "kgm2" }, izz: { value: r.Izz, unit: "kgm2" }, arm_length: { value: d.arm, unit: "mm" }, rotor_count: { value: d.rotors, unit: "count" }, tip_mass: { value: tip, unit: "kg" },
          total_mass: { value: r.mass, unit: "kg" }, hover_power: { value: h.P, unit: "W" }, esc_current: { value: x.esc, unit: "A" }, hover_current: { value: h.esc, unit: "A" },
          flight_time: { value: r.time, unit: "min" }, disk_loading: { value: r.diskLoading, unit: "Pa" }, downwash: { value: h.vi, unit: "m/s" }, wake_pressure: { value: r.wakePressure, unit: "Pa" },
          altitude: { value: d.altitude, unit: "m" }, charger_heat: { value: charger, unit: "W" } }, checks, details: r,
          certificates: [{ label: "Hover power equals T·vᵢ/(FM·η) for every rotor", ok: Math.abs(h.input - h.T * h.vi / d.fm / d.motorEfficiency / d.escEfficiency) < 1e-9 * h.input }] };
      } }],
    ["equations", { title: "Equations", ports: d => ({
        inputs: d.inputs.map(x => ({ id: x.id, label: x.label, unit: x.unit, ...(x.min !== undefined && { min: x.min }), ...(x.max !== undefined && { max: x.max }) })),
        outputs: d.outputs.map(x => ({ id: x.id, label: x.label, unit: x.unit })) }),
      async evaluate(d, ov) {
        const r = equationsResults(d, ov);
        return { outputs: Object.fromEntries(d.outputs.map(x => [x.id, { value: r.values[x.id], unit: x.unit }])), checks: r.failed.map(x => ({ id: `value-${x.id}`, text: `${x.label} has a value`, status: "unavailable", reason: `${x.expr} is not a finite number with these inputs.` })), details: r };
      } }],
    ["power", { title: "Power", ports: () => ({
        inputs: [{ id: "phase_current", label: "Motor phase current", unit: "A", min: 0, max: 1e4 }, { id: "bus_voltage", label: "Bus voltage", unit: "V", min: 1, max: 2000 }, { id: "duty", label: "Average duty", unit: "%", min: 0, max: 100 }, { id: "junction", label: "Junction temperature", unit: "degC", min: -40, max: 200 }, { id: "peak_voltage", label: "Peak bus voltage", unit: "V", min: 1, max: 2000 }],
        outputs: [{ id: "fet_rms", label: "RMS current in the hottest FET", unit: "A" }, { id: "fet_switching", label: "Hottest FET's switching loss", unit: "W" }, { id: "other_heat", label: "Heat from the other FETs", unit: "W" }, { id: "rds_on", label: "MOSFET RDS(on) at 25 °C", unit: "mohm" }, { id: "esc_loss", label: "Inverter loss", unit: "W" }, { id: "efficiency", label: "Inverter efficiency", unit: "%" }, { id: "ripple_current", label: "Bus capacitor ripple current", unit: "A" }, { id: "bus_ripple", label: "Bus ripple voltage", unit: "V" }, { id: "capacitor_mass", label: "Bus capacitor mass", unit: "kg" }, { id: "gate_power", label: "Gate drive power", unit: "W" }] }),
      async evaluate(d, ov) {
        const phase = ov.phase_current ?? d.phase, bus = ov.bus_voltage ?? d.bus, duty = (ov.duty ?? d.duty) / 100, tj = ov.junction ?? d.junction, peak = ov.peak_voltage ?? d.peakBus, r = stageLoss(d, { phase, bus, duty, tj }), c = d.capacitor, b = busRipple(d, { phase, duty }), input = phase * duty * bus, efficiency = input > 0 ? 100 * (1 - (r.total + b.loss) / input) : NaN, m = d.mosfet;
        const checks = [check("voltage", `MOSFET ${m.vds} V ≥ ${d.ringing} × ${peak} V peak bus`, m.vds >= d.ringing * peak, null, m.vds / peak, ""), check("gate", `Gate peak ≤ driver's ${d.driverPeak} A`, r.gatePeak <= d.driverPeak, null, r.gatePeak, "A"), check("efficiency", `Inverter efficiency ≥ ${d.efficiencyLimit}%`, efficiency >= d.efficiencyLimit, null, efficiency, "%"), check("junction", `Junction ≤ ${m.tjMax} °C`, tj <= m.tjMax, null, tj, "degC"),
          check("capacitor-voltage", `Bus capacitors ${c.voltage} V ≥ ${peak} V peak bus`, c.voltage >= peak, null, c.voltage / peak, ""), check("capacitor-ripple", `Ripple per capacitor ≤ ${CAPACITOR_DERATE * 100}% of its ${c.rating} A rating`, b.share <= CAPACITOR_DERATE * c.rating, null, b.share, "A"), check("bus-ripple", `Bus ripple ≤ ${c.limit}% of the bus`, b.volts <= c.limit / 100 * bus, null, b.volts, "V")];
        return { outputs: { fet_rms: { value: r.hottest.rms, unit: "A" }, fet_switching: { value: r.hottest.extra, unit: "W" }, other_heat: { value: r.other, unit: "W" }, rds_on: { value: m.rds25, unit: "mohm" }, esc_loss: { value: r.total, unit: "W" }, efficiency: { value: efficiency, unit: "%" }, ripple_current: { value: r.ripple, unit: "A" }, bus_ripple: { value: b.volts, unit: "V" }, capacitor_mass: { value: b.mass, unit: "kg" }, gate_power: { value: r.gate, unit: "W" } }, checks, details: { ...r, phase, bus, duty, tj, input, efficiency, capacitors: b } };
      } }],
    ["structures", { title: "Structures", ports: () => ({
        inputs: [{ id: "od", label: "Tube outer diameter", unit: "mm", min: 1, max: 1000 }, { id: "wall", label: "Tube wall", unit: "mm", min: .05, max: 100 }, { id: "length", label: "Free arm length", unit: "mm", min: 5, max: 5000 }, { id: "thrust", label: "Peak thrust at the tip", unit: "N", min: 0, max: 1e5 }, { id: "torque", label: "Peak torque at the tip", unit: "Nm", min: 0, max: 1e5 }, { id: "tip_mass", label: "Mass at the tip", unit: "kg", min: 0, max: 100 }, { id: "rpm_low", label: "Lowest flying rotor speed", unit: "rpm", min: 0, max: 1e6 }, { id: "rpm_high", label: "Highest rotor speed", unit: "rpm", min: 0, max: 1e6 }],
        outputs: [{ id: "stress", label: "Peak von Mises stress", unit: "MPa" }, { id: "margin", label: "Strength margin", unit: "" }, { id: "deflection", label: "Tip deflection at peak thrust", unit: "mm" }, { id: "frequency", label: "First bending mode", unit: "Hz" }, { id: "mass", label: "Arm mass", unit: "kg" }] }),
      async evaluate(d, ov) {
        const r = armResults(d, ov), checks = [check("strength", `Stress within ${d.strength} MPa at safety factor ${d.safetyFactor}`, r.margin >= 0, null, r.margin, ""), check("stiffness", `Tip deflection ≤ ${d.deflectionLimit} mm`, r.deflection <= d.deflectionLimit, null, r.deflection, "mm"),
          check("resonance", `First mode clears rotor excitation by ${d.separation}%`, r.bands.length ? !r.conflicts.length : null, r.bands.length ? (r.conflicts.length ? `${r.frequency.toFixed(0)} Hz sits in the ${r.conflicts.map(b => `${b.name} band (${b.lo.toFixed(0)}–${b.hi.toFixed(0)} Hz)`).join(" and ")}` : r.bands.map(b => `${b.name} ${b.lo.toFixed(0)}–${b.hi.toFixed(0)} Hz`).join(", ")) : "Link the rotor speeds to check", r.frequency, "Hz")];
        return { outputs: { stress: { value: r.vonMises, unit: "MPa" }, margin: { value: r.margin, unit: "" }, deflection: { value: r.deflection, unit: "mm" }, frequency: { value: r.frequency, unit: "Hz" }, mass: { value: r.armMass, unit: "kg" } }, checks, details: r };
      } }],
    ["controls", { title: "Controls", ports: () => ({
        inputs: [{ id: "mass", label: "Take-off mass", unit: "kg", min: .01, max: 1e4 }, { id: "ixx", label: "Roll inertia", unit: "kgm2", min: 1e-6, max: 1e6 }, { id: "iyy", label: "Pitch inertia", unit: "kgm2", min: 1e-6, max: 1e6 }, { id: "izz", label: "Yaw inertia", unit: "kgm2", min: 1e-6, max: 1e6 }, { id: "arm", label: "Rotor arm", unit: "mm", min: 10, max: 1e5 }, { id: "rotors", label: "Rotors", unit: "count", min: 4, max: 8 }, { id: "hover_throttle", label: "Hover throttle", unit: "%", min: 1, max: 99 }, { id: "max_thrust", label: "Full-throttle thrust per rotor", unit: "N", min: .01, max: 1e5 }, { id: "motor_tau", label: "Motor time constant", unit: "s", min: 1e-4, max: 5 }, { id: "torque_ratio", label: "Rotor torque per thrust", unit: "m", min: 0, max: 5 }],
        outputs: [{ id: "bandwidth", label: "Rate-loop bandwidth", unit: "Hz" }, { id: "phase_margin", label: "Rate-loop phase margin", unit: "" }, { id: "gain_margin", label: "Rate-loop gain margin", unit: "" }, { id: "overshoot", label: "Roll step overshoot", unit: "%" }, { id: "settle", label: "Roll step settling", unit: "s" }, { id: "rise", label: "Roll step rise", unit: "s" }, { id: "gust_drift", label: "Largest drift in the gust", unit: "m" }, { id: "altitude_settle", label: "Altitude step settling", unit: "s" }] }),
      async evaluate(d, ov) {
        const P = lastPlant = plantOf(d, ov), roll = rateLoop(d, P, 0), pitch = rateLoop(d, P, 1), pm = Math.min(roll.phaseMargin, pitch.phaseMargin), gm = Math.min(roll.gainMargin, pitch.gainMargin), L = d.limits;
        return { outputs: { bandwidth: { value: Math.min(roll.bandwidth, pitch.bandwidth), unit: "Hz" }, phase_margin: { value: pm, unit: "" }, gain_margin: { value: gm, unit: "" } },
          checks: [check("phase-margin", `Rate-loop phase margin ≥ ${L.phaseMargin}°`, pm >= L.phaseMargin, null, pm, ""), check("gain-margin", `Rate-loop gain margin ≥ ${L.gainMargin} dB`, gm >= L.gainMargin, null, gm, "")], details: { plant: P, roll, pitch, overrides: ov } };
      },
      // Flown once the system has settled; they feed nothing back.
      async after(d, { results }) {
        const r = results.controls?.details; if (!r) return {};
        const key = `flight${JSON.stringify([d, r.overrides])}`; if (!flown.has(key)) { if (flown.size > 12) flown.delete(flown.keys().next().value); flown.set(key, flightTests(d, r.plant)); }
        const f = flown.get(key), L = d.limits;
        return { outputs: { overshoot: { value: f.rollStep.overshoot, unit: "%" }, settle: { value: f.rollStep.settle, unit: "s" }, rise: { value: f.rollStep.rise, unit: "s" }, gust_drift: { value: f.gust.drift, unit: "m" }, altitude_settle: { value: f.altitude.settle, unit: "s" } },
          checks: [check("overshoot", `Roll step overshoot ≤ ${L.overshoot}%`, f.rollStep.overshoot <= L.overshoot, null, f.rollStep.overshoot, "%"), check("settle", `Roll step settles in ${L.settle} s`, f.rollStep.settle <= L.settle, null, f.rollStep.settle, "s"), check("gust", `Gust drift ≤ ${L.drift} m`, f.gust.drift <= L.drift, null, f.gust.drift, "m")], details: { ...r, tests: f } };
      } }],
    ["testing", { title: "Testing", ports: () => ({ inputs: [], outputs: [{ id: "passed", label: "Tests passed", unit: "count" }, { id: "failed", label: "Tests failed", unit: "count" }, { id: "not_run", label: "Tests not run", unit: "count" }, { id: "stale", label: "Results the design has since changed under", unit: "count" }] }),
      async evaluate() { return { outputs: {}, checks: [] }; },
      async after(d, { outputs, docs }) {
        const r = evaluatePlan(d, { outputs, table: docs.flight?.propulsion?.table, docs }), count = s => r.items.filter(i => i.status === s).length;
        return { outputs: { passed: { value: count("passed"), unit: "count" }, failed: { value: count("failed"), unit: "count" }, not_run: { value: count("not run"), unit: "count" }, stale: { value: count("stale"), unit: "count" } }, details: r,
          checks: r.items.filter(i => i.status === "failed").map(i => ({ id: `test-${i.id}`, text: `${i.title} (${i.kind})`, status: "unmet", ...(i.value !== undefined && { value: { nominal: i.value }, unit: i.unit }) })) };
      } }],
    ["parts", { title: "Parts", ports: d => ({ inputs: [], outputs: [{ id: "mass", label: "Parts mass", unit: "kg" },
        // Every numeric dimension of a part is an output another domain can use.
        ...d.parts.flatMap(p => Object.entries({ ...p.model.params, ...p.overrides }).filter(([, v]) => Number.isFinite(v)).map(([k]) => ({ id: partPort(p, k), label: `${p.name}: ${k}`, unit: "mm" })))] }),
      // Volumes come from kernel meshes, cached by model; mass multiplies by count.
      async evaluate(d) {
        let mass = 0; const volumes = {};
        for (const p of d.parts) {
          const key = JSON.stringify([p.model, p.overrides]);
          if (!partVolumes.has(key)) {
            const volume = mesh(p.model, p.overrides).then(r => r.volume); partVolumes.set(key, volume); volume.catch(() => partVolumes.delete(key));
            if (partVolumes.size > 64) partVolumes.delete(partVolumes.keys().next().value);
          }
          volumes[p.id] = await partVolumes.get(key);
          mass += volumes[p.id] * MATERIALS[p.material].density * 1e-6 * p.count;
        }
        const dims = d.parts.flatMap(p => Object.entries({ ...p.model.params, ...p.overrides }).filter(([, v]) => Number.isFinite(v)).map(([k, v]) => [partPort(p, k), { value: v, unit: "mm" }]));
        return { outputs: { mass: { value: mass, unit: "kg" }, ...Object.fromEntries(dims) }, checks: [], details: { volumes } };
      } }],
  ]);
  const partVolumes = new Map(), flown = new Map();
  let lastPlant = null;
  const active = () => store.docs.project.domains;
  const system = new System(adapters, () => ({ ...store.docs, mech: mech.state() }), id => active().includes(id));
  // A domain's commands run only while it is in the project.
  const owner = op => ({ litho: "silicon" })[op.split(".")[0]] ?? op.split(".")[0];
  const register = specs => specs.forEach(s => { const d = owner(s.op);
    bus.register({ domain: s.op.split(".")[0], ...s, ...(DOMAINS.includes(d) && { run: (a, c) => { if (!active().includes(d)) throw new Error(`${DOMAIN_INFO[d].title} is not in this project. Add it with project.add_domain ${d}.`); return s.run(a, c); } }) }); });
  register(mechanicalCommands(mech));

  // PCB
  const comp = (d, id) => { const c = d.components.find(c => c.id === id); if (!c) throw new Error(`No component ${id}. Components: ${d.components.map(c => c.id).join(", ")}.`); return c; };
  const nextId = (items, prefix) => { let i = 1; while (items.some(x => x.id === prefix + i)) i++; return prefix + i; };
  // A person may name anything they add; studies, links and tours refer to it by id.
  const idArg = { name: "id", type: "string", maxLength: 24, named: true, help: "Your own id; one is assigned otherwise" };
  const chooseId = (id, items, prefix) => { if (!id) return nextId(items, prefix); if (!ID.test(id)) throw new Error("Ids start with a letter and use letters, digits, _ and -."); if (items.some(x => x.id === id)) throw new Error(`${id} already exists.`); return id; };
  register([
    { op: "pcb.place", title: "Place a footprint", example: "pcb.place dip8 @40,20 name=\"Gate driver\"", args: [{ name: "footprint", type: "string", maxLength: 48, required: true, help: `${Object.keys(FOOTPRINTS).join(", ")}, or one made with pcb.footprint` },
      { name: "at", type: "point", help: "Board position in mm" }, { name: "name", type: "string", maxLength: 60, named: true }, { name: "ref", type: "string", maxLength: 12, named: true }],
      run: a => m("pcb", d => { if (d.components.length >= 256) throw new Error("This editor supports 256 components.");
        const fp = footprintOf(d, a.footprint); if (!fp) throw new Error(`No footprint ${a.footprint}. Footprints: ${[...Object.keys(FOOTPRINTS), ...Object.keys(d.footprints)].join(", ")}.`);
        const ref = a.ref || nextId(d.components, { connector: "J", resistor: "R", capacitor: "C", transistor: "Q", dip8: "U", via: "V" }[a.footprint] ?? "U");
        if (d.components.some(c => c.id === ref)) throw new Error(`${ref} already exists.`);
        const [x, y] = a.at || [d.width / 2, d.height / 2], snap = v => Math.round(v / d.grid) * d.grid;
        d.components.push({ id: ref, name: a.name || fp.title, footprint: a.footprint, x: snap(x), y: snap(y), rotation: 0, nets: fp.pads.map(() => "") });
        return { id: ref, message: `Placed ${ref}.` }; }) },
    { op: "pcb.footprint", title: "Make a footprint", summary: "A through-hole pattern from a package's numbers (single, dual or quad), or your own pad list. Generated footprints follow the numbers you give; check them against the part's datasheet.", example: "pcb.footprint soic14 dual pins=14 pitch=1.27 span=5.4",
      args: [{ name: "name", type: "string", maxLength: 48, required: true }, { name: "kind", type: "enum", values: ["single", "dual", "quad", "custom"], default: "dual" }, { name: "pins", type: "number", min: 1, max: 256, integer: true, named: true },
        { name: "pitch", type: "number", min: .1, max: 50, unit: "mm", named: true, default: 2.54 }, { name: "span", type: "number", min: 0, max: 200, unit: "mm", named: true, default: 7.62 }, { name: "pads", type: "json", named: true, help: "[[x, y], …] in mm, for kind custom" },
        { name: "width", type: "number", min: .5, max: 300, unit: "mm", named: true }, { name: "height", type: "number", min: .5, max: 300, unit: "mm", named: true }, { name: "title", type: "string", maxLength: 60, named: true }, { name: "mass", type: "number", min: 0, max: 1000, unit: "g", named: true }, { name: "pad", type: "number", min: .2, max: 20, unit: "mm", named: true, help: "Pad diameter" }, { name: "drill", type: "number", min: .1, max: 20, unit: "mm", named: true }, { name: "smd", type: "boolean", named: true, help: "Surface pads on the top layer" }, { name: "pad_w", type: "number", min: .05, max: 100, unit: "mm", named: true }, { name: "pad_h", type: "number", min: .05, max: 100, unit: "mm", named: true }],
      run: a => m("pcb", d => { if (Object.hasOwn(FOOTPRINTS, a.name)) throw new Error(`${a.name} is a built-in footprint.`); if (a.kind !== "custom" && !a.pins) throw new Error("Give pins, or kind custom with a pad list.");
        const list = a.kind === "custom" ? a.pads : footprintPads(a.kind, a.pins, a.pitch, a.span); if (!list?.length) throw new Error("A footprint needs pads.");
        const used = d.components.find(c => c.footprint === a.name && c.nets.length !== list.length); if (used) throw new Error(`${used.id} uses ${a.name} with ${used.nets.length} pads; remove it before changing the pad count.`);
        const span = k => Math.max(...list.map(p => Math.abs(p[k]))) * 2 + 3;
        d.footprints[a.name] = { title: a.title || a.name, width: a.width ?? Number(span(0).toFixed(3)), height: a.height ?? Number(span(1).toFixed(3)), mass: a.mass ?? 1, pad: a.pad ?? 1.6, drill: a.drill ?? .9, ...(a.smd && { smd: true }), ...(a.pad_w !== undefined && { padW: a.pad_w, padH: a.pad_h }), pads: list }; return { message: `Footprint ${a.name}: ${list.length} pads.` }; }) },
    { op: "pcb.via", title: "Add a via", summary: "A plated hole joining the layers on a net.", example: "pcb.via @30,20 GND", args: [{ name: "at", type: "point", required: true }, { name: "net", type: "string", maxLength: 32, required: true }, { name: "ref", type: "string", maxLength: 12, named: true }],
      run: a => m("pcb", d => { if (d.components.length >= 256) throw new Error("This editor supports 256 components."); const ref = a.ref || nextId(d.components, "V"); if (d.components.some(c => c.id === ref)) throw new Error(`${ref} already exists.`); d.components.push({ id: ref, name: "Via", footprint: "via", x: a.at[0], y: a.at[1], rotation: 0, nets: [a.net] }); return { id: ref, message: `Via ${ref} on ${a.net}.` }; }) },
    { op: "pcb.pour", title: "Pour copper for a net", summary: "A polygon on one layer. It joins that net's pads inside it and clears every other net's copper by the clearance rule.", example: "pcb.pour GND bottom points=[[1,1],[69,1],[69,49],[1,49]]", args: [{ name: "net", type: "string", maxLength: 32, required: true }, { name: "layer", type: "enum", values: ["top", "bottom"], default: "top" }, { name: "points", type: "json", named: true, required: true }, { name: "id", type: "string", maxLength: 24, named: true }],
      run: a => m("pcb", d => { const id = a.id || nextId(d.pours, "pour"); if (d.pours.some(p => p.id === id)) throw new Error(`${id} already exists.`); d.pours.push({ id, net: a.net, layer: a.layer, points: a.points }); return { id, message: `Poured ${a.net} on ${a.layer}.` }; }) },
    { op: "pcb.remove_pour", title: "Remove a pour", args: [{ name: "id", type: "string", required: true }], run: ({ id }) => m("pcb", d => { if (!d.pours.some(p => p.id === id)) throw new Error(`No pour ${id}.`); d.pours = d.pours.filter(p => p.id !== id); return { message: `Removed ${id}.` }; }) },
    { op: "pcb.remove_footprint", title: "Delete a footprint of your own", args: [{ name: "name", type: "string", required: true }],
      run: ({ name }) => m("pcb", d => { if (!d.footprints[name]) throw new Error(`No footprint ${name} of your own.`); const c = d.components.find(c => c.footprint === name); if (c) throw new Error(`${c.id} uses ${name}.`); delete d.footprints[name]; return { message: `Deleted ${name}.` }; }) },
    { op: "pcb.trace", title: "Draw a trace through your own points", summary: "Two pads on one net, and every point of the copper between them in mm.", example: "pcb.trace J1:1 U1:1 top 0.5 points=[[10,16],[20,16],[24,21.5]]",
      args: [{ name: "from", type: "string", required: true }, { name: "to", type: "string", required: true }, { name: "layer", type: "enum", values: ["top", "bottom"], default: "top" }, { name: "width", type: "number", min: .1, max: 5, unit: "mm", default: .5 }, { name: "points", type: "json", named: true, required: true }, { name: "id", type: "string", maxLength: 48, named: true }],
      run: a => m("pcb", d => { const r = connectPads(d, a.from, a.to, a.layer, a.width); r.points = a.points; if (a.id) r.id = a.id; if (d.routes.some(x => x.id === r.id)) throw new Error(`${r.id} already exists.`); d.routes.push(r); return { id: r.id, message: `Drew ${a.from} → ${a.to}.` }; }) },
    { op: "pcb.load_net", title: "Mark a net as carrying the linked load current", args: [{ name: "net", type: "string", required: true }, { name: "enabled", type: "boolean", default: true }],
      run: ({ net, enabled }) => m("pcb", d => { d.loadNets = d.loadNets.filter(n => n !== net); if (enabled) d.loadNets.push(net); return { message: `${net} ${enabled ? "follows" : "no longer follows"} the linked load current.` }; }) },
    { op: "pcb.move", title: "Move or rotate a component", example: "pcb.move U1 @36,26 rotation=90", args: [{ name: "id", type: "string", required: true }, { name: "at", type: "point" },
      { name: "rotation", type: "enum", values: ["0", "90", "180", "270"], named: true }],
      run: a => m("pcb", d => { const c = comp(d, a.id), changes = {}; if (a.at) [changes.x, changes.y] = a.at.map(v => Math.round(v / d.grid) * d.grid); if (a.rotation) changes.rotation = Number(a.rotation);
        const removed = moveComponent(d, c, changes); return { message: `Moved ${a.id}${removed ? `; removed ${removed} attached connection${removed > 1 ? "s" : ""}` : ""}.` }; }) },
    { op: "pcb.rename", title: "Rename a component", args: [{ name: "id", type: "string", required: true }, { name: "name", type: "string", maxLength: 60, required: true }],
      run: ({ id, name }) => m("pcb", d => { if (!name.trim()) throw new Error("Enter a name."); comp(d, id).name = name.trim(); return { message: `Renamed ${id}.` }; }) },
    { op: "pcb.load", title: "Replace the board with a design document", summary: "An emergence-studio/1 PCB document or its design object. Undo restores the previous board.", args: [{ name: "design", type: "json", required: true }],
      run: ({ design }) => { const d = validateDocument(design.format ? design : { format: "emergence-studio/1", kind: "pcb", design }); if (d.kind !== "pcb") throw new Error("Open this file in the " + d.kind + " workspace."); store.replace("pcb", d.design, "edit"); return { message: "Board loaded." }; } },
    { op: "pcb.remove", title: "Delete a component", args: [{ name: "id", type: "string", required: true }],
      run: ({ id }) => m("pcb", d => { comp(d, id); d.components = d.components.filter(c => c.id !== id); d.routes = d.routes.filter(r => !r.from.startsWith(id + ":") && !r.to.startsWith(id + ":")); delete d.powers[id]; return { message: `Deleted ${id}.` }; }) },
    { op: "pcb.net", title: "Assign a pad's net", summary: "NC marks a pin deliberately left unconnected; - clears it.", example: "pcb.net U1:4 GATE", args: [{ name: "pad", type: "string", required: true }, { name: "net", type: "string", maxLength: 32, required: true }],
      run: ({ pad, net }) => m("pcb", d => { const [id, n] = pad.split(":"), c = comp(d, id), i = Number(n) - 1;
        if (!Number.isInteger(i) || i < 0 || i >= c.nets.length) throw new Error(`${id} has pads 1–${c.nets.length}.`);
        if (!/^[a-zA-Z0-9_+./-]*$/.test(net)) throw new Error("Net names use letters, digits and _+./-.");
        c.nets[i] = net === "-" ? "" : net; d.routes = d.routes.filter(r => r.from !== pad && r.to !== pad); return { message: `${pad} → ${c.nets[i] || "no net"}.` }; }) },
    { op: "pcb.connect", title: "Connect two pads", summary: "Straight-and-45° connection between two pads on the same net.", example: "pcb.connect J1:1 U1:1 top 0.5",
      args: [{ name: "from", type: "string", required: true }, { name: "to", type: "string", required: true }, { name: "layer", type: "enum", values: ["top", "bottom"], default: "top" }, { name: "width", type: "number", min: .1, max: 5, unit: "mm", default: .5 }],
      run: a => m("pcb", d => { const r = connectPads(d, a.from, a.to, a.layer, a.width); d.routes.push(r); return { id: r.id, message: `Connected ${a.from} → ${a.to}.` }; }) },
    { op: "pcb.route", title: "Route around obstacles", summary: "A* on the snap grid, clear of other nets' copper by the clearance rule.", example: "pcb.route U1:4 R1:1 top 0.5",
      args: [{ name: "from", type: "string", required: true }, { name: "to", type: "string", required: true }, { name: "layer", type: "enum", values: ["top", "bottom"], default: "top" }, { name: "width", type: "number", min: .1, max: 5, unit: "mm", default: .5 }],
      run: a => m("pcb", d => { if (d.routes.length >= 512) throw new Error("This editor supports 512 connections."); const r = autoroute(d, a.from, a.to, a.layer, a.width); d.routes.push(r); return { id: r.id, message: `Routed ${a.from} → ${a.to} with ${r.points.length - 1} segments.` }; }) },
    { op: "pcb.route_all", title: "Route every open connection", summary: "Routes each remaining airwire, trying top then bottom.", args: [{ name: "width", type: "number", min: .1, max: 5, unit: "mm", default: .4 }],
      run: ({ width }) => m("pcb", d => { let done = 0; const failed = [];
        for (const net of pcbChecks(d).nets) for (const [a, b] of net.airwires) {
          let ok = false; for (const layer of ["top", "bottom"]) try { d.routes.push(autoroute(d, a.id, b.id, layer, width)); ok = true; done++; break; } catch { /* try the next layer */ }
          if (!ok) failed.push(`${a.id}–${b.id}`);
        }
        return { message: `Routed ${done} connection${done === 1 ? "" : "s"}${failed.length ? `; no clear path for ${failed.join(", ")}` : ""}.` }; }) },
    { op: "pcb.net_width", title: "Set the width of every trace on a net", example: "pcb.net_width OUT 4", args: [{ name: "net", type: "string", required: true }, { name: "width", type: "number", min: .1, max: 5, unit: "mm", required: true }],
      run: ({ net, width }) => m("pcb", d => { const rs = d.routes.filter(r => netOf(d, r) === net); if (!rs.length) throw new Error(`No traces on ${net}.`); rs.forEach(r => r.width = width); return { message: `${rs.length} ${net} trace${rs.length > 1 ? "s" : ""} now ${width} mm.` }; }) },
    { op: "pcb.unroute", title: "Remove a connection", args: [{ name: "id", type: "string", required: true }],
      run: ({ id }) => m("pcb", d => { if (!d.routes.some(r => r.id === id)) throw new Error(`No connection ${id}.`); d.routes = d.routes.filter(r => r.id !== id); return { message: `Removed ${id}.` }; }) },
    { op: "pcb.current", title: "Set a net's RMS current", example: "pcb.current GATE 0.4", args: [{ name: "net", type: "string", required: true }, { name: "amps", type: "number", min: 0, max: 500, unit: "A", required: true }],
      run: ({ net, amps }) => m("pcb", d => { d.currents[net] = amps; return { message: `${net} carries ${amps} A RMS.` }; }) },
    { op: "pcb.power", title: "Set a component's dissipation", example: "pcb.power R1 0.05", args: [{ name: "id", type: "string", required: true }, { name: "watts", type: "number", min: 0, max: 500, unit: "W", required: true }],
      run: ({ id, watts }) => m("pcb", d => { comp(d, id); d.powers[id] = watts; return { message: `${id} dissipates ${watts} W.` }; }) },
    { op: "pcb.board", title: "Set board size, rules and stackup", example: "pcb.board width=90 clearance=0.2 epsr=4.3", args: [
      ...[["width", 10, 300, "mm"], ["height", 10, 300, "mm"], ["grid", .1, 10, "mm"], ["clearance", .05, 5, "mm"], ["min_width", .1, 5, "mm"], ["edge", 0, 10, "mm"], ["density", 1, 1000, "A/mm²"], ["thickness", .2, 6, "mm"], ["copper", .25, 6, "oz"], ["epsr", 1, 20, ""]]
        .map(([name, min, max, unit]) => ({ name, type: "number", min, max, unit, named: true }))],
      run: a => m("pcb", d => { for (const k of ["width", "height", "grid"]) if (a[k] !== undefined) d[k] = a[k];
        for (const [k, key] of [["clearance", "clearance"], ["min_width", "minWidth"], ["edge", "edge"], ["density", "density"]]) if (a[k] !== undefined) d.rules[key] = a[k];
        for (const k of ["thickness", "copper", "epsr"]) if (a[k] !== undefined) d.stackup[k] = a[k]; return { message: "Board settings updated." }; }) },
  ]);
  // Silicon
  const block = (d, id) => { const b = d.blocks.find(b => b.id === id); if (!b) throw new Error(`No block ${id}. Blocks: ${d.blocks.map(b => b.id).join(", ")}.`); return b; };
  register([
    { op: "silicon.add_block", title: "Add a floorplan block", example: "silicon.add_block logic \"Crypto engine\" @600,500 width=200 height=150", args: [{ name: "kind", type: "enum", values: ["logic", "memory", "io", "analog"], required: true }, idArg,
      { name: "name", type: "string", maxLength: 60 }, { name: "at", type: "point" }, { name: "width", type: "number", min: 1, max: 10000, unit: "µm", named: true, default: 200 },
      { name: "height", type: "number", min: 1, max: 10000, unit: "µm", named: true, default: 160 }, { name: "power", type: "number", min: 0, max: 1e6, unit: "mW", named: true }],
      run: a => m("silicon", d => { if (d.blocks.length >= 64) throw new Error("This editor supports 64 blocks."); const id = chooseId(a.id, d.blocks, a.kind), [x, y] = a.at || [d.margin, d.margin];
        d.blocks.push({ id, name: a.name || { logic: "Logic block", memory: "Memory macro", io: "I/O block", analog: "Analog block" }[a.kind], kind: a.kind, x, y, width: a.width, height: a.height, fixed: false, power: a.power ?? BLOCK_POWER[a.kind] });
        return { id, message: `Added ${id}.` }; }) },
    { op: "silicon.update_block", title: "Edit a block", example: "silicon.update_block dsp power=450 x=480", args: [{ name: "id", type: "string", required: true },
      ...["x", "y", "width", "height"].map(name => ({ name, type: "number", min: name === "x" || name === "y" ? -10000 : 1, max: name === "x" || name === "y" ? 20000 : 10000, unit: "µm", named: true })),
      { name: "power", type: "number", min: 0, max: 1e6, unit: "mW", named: true }, { name: "name", type: "string", maxLength: 60, named: true }, { name: "fixed", type: "boolean", named: true }],
      run: ({ id, ...c }) => m("silicon", d => { const b = block(d, id);
        if (b.fixed && c.fixed !== false && ["x", "y", "width", "height"].some(k => c[k] !== undefined)) throw new Error("Unlock this block before changing its placement or size.");
        for (const [k, v] of Object.entries(c)) if (v !== undefined) b[k] = v; return { message: `Updated ${b.name}.` }; }) },
    { op: "silicon.remove_block", title: "Delete a block", args: [{ name: "id", type: "string", required: true }],
      run: ({ id }) => m("silicon", d => { if (block(d, id).fixed) throw new Error("Unlock this block before deleting it."); d.blocks = d.blocks.filter(b => b.id !== id); d.chip.rtl = d.chip.rtl.filter(r => r.block !== id); return { message: `Deleted ${id}.` }; }) },
    { op: "silicon.load", title: "Replace the floorplan with a design document", args: [{ name: "design", type: "json", required: true }],
      run: ({ design }) => { const d = validateDocument(design.format ? design : { format: "emergence-studio/1", kind: "silicon", design }); if (d.kind !== "silicon") throw new Error("Open this file in the " + d.kind + " workspace."); store.replace("silicon", d.design, "edit"); return { message: "Floorplan loaded." }; } },
    { op: "silicon.core", title: "Set die and core constraints", args: [...[["width", 10, 10000], ["height", 10, 10000], ["grid", .1, 100], ["margin", 0, 1000], ["spacing", 0, 1000], ["target", 1, 100]].map(([name, min, max]) => ({ name, type: "number", min, max, named: true }))],
      run: a => m("silicon", d => { for (const [k, v] of Object.entries(a)) if (v !== undefined) d[k] = v; return { message: "Core constraints updated." }; }) },
    { op: "silicon.chip", title: "Set the chip flow conditions", summary: "Clock, supply, input activity, instances and package θJA.", example: "silicon.chip frequency=150 vdd=1.32",
      args: [{ name: "frequency", type: "number", min: .001, max: 5000, unit: "MHz", named: true }, { name: "vdd", type: "number", min: 1.08, max: 1.32, unit: "V", named: true },
        { name: "activity", type: "number", min: 0, max: 100, unit: "%", named: true }, { name: "instances", type: "number", min: 1, max: 100000, integer: true, named: true },
        { name: "theta", type: "number", min: .1, max: 500, unit: "K/W", named: true }, { name: "ambient", type: "number", min: -40, max: 125, unit: "°C", named: true }],
      run: a => m("silicon", d => { for (const [k, v] of Object.entries(a)) if (v !== undefined) d.chip[k] = v; return { message: "Chip flow settings updated." }; }) },
    { op: "silicon.rtl", title: "Implement a block with RTL", summary: "Source mini8, tt-motor (public Tiny Tapeout motor controller) or custom Verilog in the native subset; none removes it.", example: "silicon.rtl pwm tt-motor",
      args: [{ name: "block", type: "string", required: true }, { name: "source", type: "enum", values: ["mini8", "tt-motor", "custom", "none"], required: true }, { name: "verilog", type: "string", maxLength: 65536, named: true }],
      run: ({ block: id, source, verilog: text }) => m("silicon", d => { block(d, id); d.chip.rtl = d.chip.rtl.filter(r => r.block !== id);
        if (source !== "none") { if (source === "custom" && !text?.trim()) throw new Error("Custom RTL needs Verilog source."); d.chip.rtl.push({ block: id, source, verilog: source === "custom" ? text : null }); }
        return { message: source === "none" ? `${id} uses its estimated power.` : `${id} is implemented by ${source} RTL; synthesis reruns on the next solve.` }; }) },
    { op: "litho.set", title: "Set lithography optics", example: "litho.set na=0.85 sigma_outer=0.8", args: [...[["wavelength", 13, 500, "nm"], ["na", .1, 1.5, ""], ["sigma_outer", 0, 1, ""], ["sigma_inner", 0, 1, ""], ["threshold", .05, .95, ""], ["pixel", 1, 200, "nm"]]
      .map(([name, min, max, unit]) => ({ name, type: "number", min, max, unit, named: true })), { name: "source_points", type: "number", min: 1, max: 9, integer: true, named: true },
      { name: "steepness", type: "number", min: 1, max: 200, named: true }, { name: "grid", type: "enum", values: ["32", "64", "128"], named: true }],
      run: a => m("litho", d => { if (a.grid) d.n = Number(a.grid); for (const [k, key] of [["wavelength", "wavelength"], ["na", "na"], ["sigma_outer", "sigmaOuter"], ["sigma_inner", "sigmaInner"], ["threshold", "threshold"], ["pixel", "pixel"], ["source_points", "sourcePoints"], ["steepness", "steepness"]]) if (a[k] !== undefined) d[key] = a[k]; return { message: "Lithography settings updated." }; }) },
    { op: "litho.shape", title: "Add a rectangle to the mask clip", example: "litho.shape 200 200 160 600", args: ["x", "y", "w", "h"].map(name => ({ name, type: "number", min: 0, max: 1e5, unit: "nm", required: true })),
      run: a => m("litho", d => { d.shapes.push(a); return { message: "Added a mask rectangle." }; }) },
    { op: "litho.remove_shape", title: "Remove a mask rectangle", args: [{ name: "index", type: "number", min: 1, integer: true, required: true }],
      run: ({ index }) => m("litho", d => { if (!d.shapes[index - 1]) throw new Error(`The clip has ${d.shapes.length} rectangles.`); d.shapes.splice(index - 1, 1); return { message: `Removed rectangle ${index}.` }; }) },
    { op: "litho.clear", title: "Clear the mask clip", run: () => m("litho", d => { d.shapes = []; return { message: "Cleared the clip." }; }) },
  ]);
  // Optics
  const element = (d, id) => { const e = d.elements.find(e => e.id === id); if (!e) throw new Error(`No element ${id}. Elements: ${d.elements.map(e => e.id).join(", ")}.`); return e; };
  register([
    { op: "optics.set", title: "Set the beam, wavelength and target", example: "optics.set power=3 wavelength=980 target_z=120", args: [
      ...[["wavelength", 300, 2500, "nm"], ["power", 0, 100, "W"], ["diameter", .1, 100, "mm"], ["target_z", -1000, 5000, "mm"], ["target_diameter", .01, 500, "mm"]].map(([name, min, max, unit]) => ({ name, type: "number", min, max, unit, named: true })),
      { name: "rays", type: "number", min: 16, max: 4000, integer: true, named: true }, { name: "pixel", type: "number", min: .5, max: 50, unit: "µm", named: true }, { name: "required_gsd", type: "number", min: 0, max: 1000, unit: "cm", named: true, help: "0 sets no requirement" },
      { name: "object_distance", type: "number", min: .001, max: 1e6, unit: "m", named: true }, { name: "target_name", type: "string", maxLength: 60, named: true }, { name: "fields", type: "string", maxLength: 60, named: true, help: "Field angles in degrees, e.g. 0,14,20" }],
      run: a => m("optics", d => { if (a.pixel !== undefined) d.pixel = a.pixel; if (a.required_gsd !== undefined) d.requiredGsd = a.required_gsd; if (a.object_distance !== undefined) d.objectDistance = a.object_distance; if (a.target_name) d.target.name = a.target_name;
        if (a.fields !== undefined) { d.fields = a.fields.split(",").map(Number); if (!d.fields.length || !d.fields.every(Number.isFinite)) throw new Error("Fields are angles such as 0,14,20."); }
        if (a.wavelength !== undefined) d.wavelength = a.wavelength; if (a.power !== undefined) d.source.power = a.power; if (a.diameter !== undefined) d.source.diameter = a.diameter;
        if (a.target_z !== undefined) d.target.z = a.target_z; if (a.target_diameter !== undefined) d.target.diameter = a.target_diameter; if (a.rays !== undefined) d.rays = a.rays; return { message: "Optical settings updated." }; }) },
    { op: "optics.add", title: "Add a lens, window or stop", example: "optics.add lens z=70 r1=80 r2=-80 thickness=5 diameter=25 glass=N-BK7", args: [{ name: "kind", type: "enum", values: OPTIC_KINDS, required: true }, idArg, { name: "name", type: "string", maxLength: 60, named: true },
      ...[["z", -1000, 5000], ["r1", -1e5, 1e5], ["r2", -1e5, 1e5], ["thickness", .05, 200], ["diameter", .1, 500]].map(([name, min, max]) => ({ name, type: "number", min, max, unit: "mm", named: true })),
      { name: "glass", type: "enum", values: Object.keys(GLASSES), named: true, default: "N-BK7" }],
      run: a => m("optics", d => { const id = chooseId(a.id, d.elements, { lens: "L", window: "W", stop: "S" }[a.kind]), last = Math.max(0, ...d.elements.map(e => e.z + (e.thickness || 0)));
        const e = { id, name: a.name || { lens: "Lens", window: "Window", stop: "Aperture stop" }[a.kind], kind: a.kind, z: a.z ?? last + 10, diameter: a.diameter ?? 25 };
        if (a.kind !== "stop") Object.assign(e, { thickness: a.thickness ?? 4, glass: a.glass });
        if (a.kind === "lens") Object.assign(e, { r1: a.r1 ?? 60, r2: a.r2 ?? -60 });
        d.elements.push(e); return { id, message: `Added ${e.name} (${id}).` }; }) },
    { op: "optics.update", title: "Edit an optical element", example: "optics.update L1 r1=48 z=22", args: [{ name: "id", type: "string", required: true }, { name: "name", type: "string", maxLength: 60, named: true },
      ...[["z", -1000, 5000], ["r1", -1e5, 1e5], ["r2", -1e5, 1e5], ["thickness", .05, 200], ["diameter", .1, 500]].map(([name, min, max]) => ({ name, type: "number", min, max, unit: "mm", named: true })),
      { name: "glass", type: "enum", values: Object.keys(GLASSES), named: true }],
      run: ({ id, ...c }) => m("optics", d => { const e = element(d, id); for (const [k, v] of Object.entries(c)) if (v !== undefined) { if (e.kind !== "lens" && ["r1", "r2"].includes(k)) throw new Error(`${id} has no curvature.`); if (e.kind === "stop" && ["thickness", "glass"].includes(k)) throw new Error("A stop has no thickness or glass."); e[k] = v; } return { message: `Updated ${e.name}.` }; }) },
    { op: "optics.remove", title: "Remove an optical element", args: [{ name: "id", type: "string", required: true }],
      run: ({ id }) => m("optics", d => { element(d, id); d.elements = d.elements.filter(e => e.id !== id); return { message: `Removed ${id}.` }; }) },
    { op: "optics.focus", title: "Move the target to best focus", summary: "Golden-section search for the least RMS spot radius on real rays.",
      run: () => m("optics", d => { d.target.z = Math.round(bestFocus(d) * 1000) / 1000; return { message: `Target moved to best focus at z = ${d.target.z} mm.` }; }) },
    { op: "optics.design_lens", title: "Solve a plano-convex radius for a focal length", example: "optics.design_lens L1 100", args: [{ name: "id", type: "string", required: true }, { name: "efl", type: "number", min: 1, max: 1e5, unit: "mm", required: true }],
      run: ({ id, efl }) => m("optics", d => { const e = element(d, id); if (e.kind !== "lens") throw new Error(`${id} is not a lens.`);
        // Curved side first; thickness has no effect on f for a plano-convex lens.
        const n = refractiveIndex(e.glass, d.wavelength); e.r1 = Math.round(efl * (n - 1) * 1e4) / 1e4; e.r2 = 0;
        if (e.r1 < e.diameter / 2) throw new Error("That focal length needs a radius smaller than the semi-diameter. Reduce the diameter or increase the focal length.");
        return { message: `${id}: R1 = ${e.r1} mm, flat back, in ${e.glass} at ${d.wavelength} nm (n = ${n.toFixed(5)}).` }; }) },
  ]);
  // Fluids
  const fel = (d, id) => { const e = d.elements.find(e => e.id === id); if (!e) throw new Error(`No element ${id}. Elements: ${d.elements.map(e => e.id).join(", ")}.`); return e; };
  register([
    { op: "fluids.fan", title: "Choose the fan", example: "fluids.fan axial80hp", args: [{ name: "preset", type: "enum", values: [...Object.keys(FANS), "custom"], required: true },
      { name: "shutoff", type: "number", min: 1, max: 5000, unit: "Pa", named: true }, { name: "free", type: "number", min: .00001, max: 10, unit: "m³/s", named: true }],
      run: a => m("fluids", d => { if (a.preset === "custom") d.fan = { preset: "custom", shutoff: a.shutoff ?? d.fan.shutoff, free: a.free ?? d.fan.free }; else d.fan = { preset: a.preset, shutoff: FANS[a.preset].shutoff, free: FANS[a.preset].free };
        return { message: `Fan: ${a.preset === "custom" ? "custom curve" : FANS[a.preset].title}.` }; }) },
    { op: "fluids.set", title: "Set air conditions and the rise limit", args: [{ name: "temperature", type: "number", min: -40, max: 150, unit: "°C", named: true }, { name: "pressure", type: "number", min: 5e4, max: 2e5, unit: "Pa", named: true }, { name: "heat_limit", type: "number", min: .1, max: 200, unit: "K", named: true }, { name: "heat", type: "number", min: 0, max: 1e6, unit: "W", named: true, help: "Heat the air carries, unless a link drives it" }],
      run: a => m("fluids", d => { if (a.heat !== undefined) d.heatLoad = a.heat; if (a.temperature !== undefined) d.temperature = a.temperature; if (a.pressure !== undefined) d.pressure = a.pressure; if (a.heat_limit !== undefined) d.heatLimit = a.heat_limit; return { message: "Air conditions updated." }; }) },
    { op: "fluids.add", title: "Add a duct, fitting, filter or fin channels", example: "fluids.add duct name=Riser width=60 height=30 length=120", args: [{ name: "kind", type: "enum", values: FLUID_KINDS, required: true }, idArg, { name: "name", type: "string", maxLength: 60, named: true },
      ...[["width", .5, 2000, "mm"], ["height", .5, 2000, "mm"], ["length", .1, 1e5, "mm"], ["roughness", 0, 10, "mm"], ["K", 0, 100, ""], ["coefficient", 0, 1e4, "Pa·s/m"], ["gap", .05, 200, "mm"], ["Kin", 0, 10, ""], ["Kout", 0, 10, ""]].map(([name, min, max, unit]) => ({ name, type: "number", min, max, unit, named: true })),
      { name: "fins", type: "number", min: 2, max: 200, integer: true, named: true }, { name: "after", type: "string", named: true }],
      run: a => m("fluids", d => { if (a.kind === "channels" && d.elements.some(e => e.kind === "channels")) throw new Error("An air path holds one set of fin channels."); const id = chooseId(a.id, d.elements, a.kind), e = { id, name: a.name || { duct: "Duct", fitting: "Fitting", filter: "Filter", channels: "Fin channels" }[a.kind], kind: a.kind, width: a.width ?? (a.kind === "channels" ? 40 : 60), height: a.height ?? (a.kind === "channels" ? 20 : 30) };
        if (a.kind === "duct") Object.assign(e, { length: a.length ?? 100, roughness: a.roughness ?? .0015 }); if (a.kind === "fitting") e.K = a.K ?? .5; if (a.kind === "filter") e.coefficient = a.coefficient ?? 25;
        if (a.kind === "channels") Object.assign(e, { fins: a.fins ?? 9, gap: a.gap ?? 3.65, length: a.length ?? 50, Kin: a.Kin ?? .5, Kout: a.Kout ?? 1 });
        const at = a.after ? d.elements.findIndex(x => x.id === a.after) : d.elements.findIndex(x => x.kind === "channels");
        d.elements.splice(at < 0 ? d.elements.length : a.after ? at + 1 : at, 0, e); return { id, message: `Added ${e.name} (${id}).` }; }) },
    { op: "fluids.update", title: "Edit an air-path element", example: "fluids.update duct length=240", args: [{ name: "id", type: "string", required: true }, { name: "name", type: "string", maxLength: 60, named: true },
      ...[["width", .5, 2000], ["height", .5, 2000], ["length", .1, 1e5], ["roughness", 0, 10], ["K", 0, 100], ["coefficient", 0, 1e4], ["gap", .05, 200], ["Kin", 0, 10], ["Kout", 0, 10]].map(([name, min, max]) => ({ name, type: "number", min, max, named: true })),
      { name: "fins", type: "number", min: 2, max: 200, integer: true, named: true }],
      run: ({ id, ...c }) => m("fluids", d => { const e = fel(d, id); for (const [k, v] of Object.entries(c)) if (v !== undefined) { if (!(k in e)) throw new Error(`${e.name} has no ${k}.`); e[k] = v; } return { message: `Updated ${e.name}.` }; }) },
    { op: "fluids.remove", title: "Remove an air-path element", args: [{ name: "id", type: "string", required: true }],
      run: ({ id }) => m("fluids", d => { const e = fel(d, id); d.elements = d.elements.filter(x => x.id !== id); return { message: `Removed ${e.name}.` }; }) },
  ]);
  // Architecture
  const room = (d, id) => { const r = d.rooms.find(r => r.id === id); if (!r) throw new Error(`No room ${id}. Rooms: ${d.rooms.map(r => r.id).join(", ")}.`); return r; };
  const roomArgs = [...[["width", .5, 200, "m"], ["depth", .5, 200, "m"], ["height", 2, 20, "m"], ["people", 0, 2000, ""], ["lighting", 0, 50, "W/m²"], ["supply", 0, 1e6, "m³/h"], ["supply_temp", 5, 40, "°C"], ["setpoint", 10, 35, "°C"], ["roof_u", 0, 5, "W/m²K"]]
    .map(([name, min, max, unit]) => ({ name, type: "number", min, max, unit, named: true })), { name: "mode", type: "enum", values: ["supplied", "conditioned"], named: true }, { name: "equipment", type: "boolean", named: true }, { name: "name", type: "string", maxLength: 60, named: true }];
  const applyRoom = (d, r, a) => { for (const [k, key] of [["width", "width"], ["depth", "depth"], ["height", "height"], ["people", "people"], ["lighting", "lighting"], ["supply", "supply"], ["supply_temp", "supplyTemp"], ["setpoint", "setpoint"], ["roof_u", "roofU"], ["mode", "mode"], ["name", "name"]]) if (a[k] !== undefined) r[key] = a[k];
    // One room holds the equipment; r is already in d.rooms.
    if (a.at) [r.x, r.y] = a.at; if (a.equipment !== undefined) { if (a.equipment) for (const x of d.rooms) x.equipment = x === r; else r.equipment = false; if (!d.rooms.some(x => x.equipment)) d.rooms[0].equipment = true; } };
  register([
    { op: "architecture.add_room", title: "Add a room", example: "architecture.add_room Workshop @13,0 width=6 depth=9", args: [{ name: "label", type: "string", maxLength: 60, required: true }, { name: "at", type: "point", help: "Plan position in m" }, idArg, ...roomArgs.filter(a => a.name !== "name")],
      run: a => m("architecture", d => { if (d.rooms.length >= 32) throw new Error("This plan supports 32 rooms."); const id = chooseId(a.id, d.rooms, "room"), right = Math.max(0, ...d.rooms.map(r => r.x + r.width));
        const w = (ext = true) => ({ exterior: ext, U: .35, window: 0 }), r = { id, name: a.label, x: right, y: 0, width: 5, depth: 5, height: 3, people: 2, lighting: 9, mode: "conditioned", supply: 200, supplyTemp: 16, setpoint: 24, equipment: false, roofU: .25, walls: { north: w(), east: w(), south: w(), west: w() } };
        d.rooms.push(r); applyRoom(d, r, a); return { id, message: `Added ${r.name} (${id}).` }; }) },
    { op: "architecture.update_room", title: "Edit a room", example: "architecture.update_room lab supply=1200 people=5", args: [{ name: "id", type: "string", required: true }, { name: "at", type: "point" }, ...roomArgs],
      run: ({ id, ...a }) => m("architecture", d => { const r = room(d, id); applyRoom(d, r, a); return { message: `Updated ${r.name}.` }; }) },
    { op: "architecture.wall", title: "Set a wall: exterior, U-value and window", example: "architecture.wall lab west window=4", args: [{ name: "id", type: "string", required: true }, { name: "side", type: "enum", values: SIDES, required: true },
      { name: "exterior", type: "boolean", named: true }, { name: "u", type: "number", min: .05, max: 6, unit: "W/m²K", named: true }, { name: "window", type: "number", min: 0, max: 5000, unit: "m²", named: true }],
      run: ({ id, side, exterior, u, window }) => m("architecture", d => { const w = room(d, id).walls[side]; if (exterior !== undefined) w.exterior = exterior; if (u !== undefined) w.U = u; if (window !== undefined) w.window = window; return { message: `${id} ${side} wall updated.` }; }) },
    { op: "architecture.remove_room", title: "Remove a room", args: [{ name: "id", type: "string", required: true }],
      run: ({ id }) => m("architecture", d => { room(d, id); d.rooms = d.rooms.filter(r => r.id !== id); if (d.rooms.length && !d.rooms.some(r => r.equipment)) d.rooms[0].equipment = true; return { message: `Removed ${id}.` }; }) },
    { op: "architecture.site", title: "Set climate, glazing, unit count and limit", example: "architecture.site outdoor=35 units=40", args: [...[["outdoor", -40, 55, "°C"], ["solar", 0, 1200, "W/m²"], ["glazing_u", .3, 7, "W/m²K"], ["shgc", .05, .9, ""], ["limit", 5, 80, "°C"], ["equipment", 0, 1e6, "W"]].map(([name, min, max, unit]) => ({ name, type: "number", min, max, unit, named: true })),
      { name: "units", type: "number", min: 0, max: 10000, integer: true, named: true }, { name: "name", type: "string", maxLength: 60, named: true }],
      run: a => m("architecture", d => { if (a.outdoor !== undefined) d.site.outdoor = a.outdoor; if (a.solar !== undefined) d.site.solar = a.solar; if (a.name !== undefined) d.site.name = a.name;
        if (a.glazing_u !== undefined) d.glazing.U = a.glazing_u; if (a.shgc !== undefined) d.glazing.shgc = a.shgc; if (a.units !== undefined) d.units = a.units; if (a.limit !== undefined) d.limit = a.limit; if (a.equipment !== undefined) d.equipmentPower = a.equipment; return { message: "Site updated." }; }) },
  ]);
  // Parts
  const pt = (d, id) => { const p = d.parts.find(p => p.id === id); if (!p) throw new Error(`No part ${id}. Parts: ${d.parts.map(p => p.id).join(", ") || "none"}.`); return p; };
  const placementArgs = [{ name: "rz", type: "number", min: -360, max: 360, unit: "°", named: true }, { name: "count", type: "number", min: 1, max: 64, integer: true, named: true }, { name: "pattern", type: "enum", values: ["single", "polar", "stack"], named: true }, { name: "stack", type: "number", min: -2000, max: 2000, unit: "mm", named: true, help: "Spacing along z" }];
  register([
    { op: "parts.add", title: "Add a solid from a template", summary: "Block and disc are plain starts to build on with parts.solid and parts.modify.", example: "parts.add flange x=0 y=-90", args: [{ name: "template", type: "enum", values: Object.keys(TEMPLATES), required: true }, idArg, { name: "name", type: "string", maxLength: 60, named: true },
      { name: "x", type: "number", min: -2000, max: 2000, unit: "mm", named: true, default: 0 }, { name: "y", type: "number", min: -2000, max: 2000, unit: "mm", named: true, default: 0 }, { name: "z", type: "number", min: -2000, max: 2000, unit: "mm", named: true, default: 0 },
      { name: "material", type: "enum", values: Object.keys(MATERIALS), named: true }, ...placementArgs],
      run: a => m("parts", d => { if (d.parts.length >= 24) throw new Error("This project supports 24 parts."); const id = chooseId(a.id, d.parts, a.template), t = TEMPLATES[a.template];
        d.parts.push({ id, name: a.name || t.title, template: a.template, model: structuredClone(t.model), overrides: {}, x: a.x, y: a.y, z: a.z, rz: a.rz ?? 0, count: a.count ?? 1, pattern: a.pattern ?? "single", stack: a.stack ?? 0, material: a.material || t.material, domain: "Mechanical" }); return { id, message: `Added ${t.title} (${id}).` }; }) },
    { op: "parts.solid", title: "Combine a shape with a part", summary: "Add, cut or keep the overlap with a box, cylinder, sphere, cone, torus, or a sketched extrude or revolve. Positions are the shape's centre; any dimension may be an expression of the part's parameters.", example: "parts.solid block1 cut cylinder z=10 radius=4 height=30",
      args: [{ name: "id", type: "string", required: true }, { name: "how", type: "enum", values: COMBINE, required: true }, { name: "shape", type: "enum", values: FEATURE_SHAPES, required: true },
        ...["x", "y", "z", "width", "depth", "height", "radius", "radius2", "round", "angle", "k"].map(name => ({ name, type: "expr", named: true })), { name: "points", type: "json", named: true, help: "Sketch points [[x, y], …] in mm" }],
      run: ({ id, how, k, ...a }) => m("parts", d => { const p = pt(d, id); p.model.shape = combine(p.model.shape, how, primitive(a), k); p.template = "custom"; return { message: `${p.name}: ${how} ${a.shape}.` }; }) },
    { op: "parts.modify", title: "Apply a feature to a whole part", summary: "Shell, offset, round, fillet, pattern, polar pattern, mirror, move, rotate, scale or clip what has been built so far.", example: "parts.modify block1 pattern nx=3 x=25",
      args: [{ name: "id", type: "string", required: true }, { name: "op", type: "enum", values: MODIFY, required: true }, { name: "value", type: "expr" }, { name: "axis", type: "enum", values: ["x", "y", "z"], named: true },
        ...["x", "y", "z", "nx", "ny", "nz", "count"].map(name => ({ name, type: "expr", named: true })), { name: "min", type: "json", named: true }, { name: "max", type: "json", named: true }],
      run: ({ id, ...a }) => m("parts", d => { const p = pt(d, id); p.model.shape = modify(p.model.shape, a, p.model.cell); p.template = "custom"; return { message: `${p.name}: ${a.op}.` }; }) },
    { op: "parts.define", title: "Define a parameter of a part", summary: "A number or an expression of other parameters, usable in any dimension.", example: "parts.define block1 wall 2.5",
      args: [{ name: "id", type: "string", required: true }, { name: "name", type: "string", maxLength: 32, required: true }, { name: "value", type: "expr", required: true }],
      run: ({ id, name, value }) => m("parts", d => { const p = pt(d, id); if (!/^[a-zA-Z_][a-zA-Z0-9_]{0,31}$/.test(name)) throw new Error("Parameter names use letters, digits and _."); p.model.params[name] = value; delete p.overrides[name]; return { message: `${p.name}: ${name} = ${value}.` }; }) },
    { op: "parts.param", title: "Set a part's model parameter", example: "parts.param gear1 teeth 30", args: [{ name: "id", type: "string", required: true }, { name: "name", type: "string", required: true }, { name: "value", type: "number", required: true }],
      run: ({ id, name, value }) => m("parts", d => { const p = pt(d, id); if (!Object.hasOwn(p.model.params, name)) throw new Error(`${p.name} has parameters ${Object.keys(p.model.params).join(", ")}.`); p.overrides[name] = value; return { message: `${p.name}: ${name} = ${value}.` }; }) },
    { op: "parts.update", title: "Place, pattern or rename a part", summary: "Position, material, and copies: polar about z, or a stack along z.", example: "parts.update arm1 pattern=polar count=4 rz=45", args: [{ name: "id", type: "string", required: true }, ...["x", "y", "z"].map(name => ({ name, type: "number", min: -2000, max: 2000, unit: "mm", named: true })),
      { name: "material", type: "enum", values: Object.keys(MATERIALS), named: true }, { name: "name", type: "string", maxLength: 60, named: true }, { name: "domain", type: "string", maxLength: 24, named: true }, ...placementArgs, { name: "cell", type: "number", min: .01, max: 50, unit: "mm", named: true, help: "Mesh resolution" }],
      run: ({ id, cell, ...c }) => m("parts", d => { const p = pt(d, id); for (const [k, v] of Object.entries(c)) if (v !== undefined) p[k] = v; if (cell !== undefined) p.model.cell = cell; return { message: `Updated ${p.name}.` }; }) },
    { op: "parts.model", title: "Replace a part's feature tree", summary: "JSON in the fluid cad model format: params, shape (ops and kits), optional cell.", args: [{ name: "id", type: "string", required: true }, { name: "model", type: "json", required: true }, { name: "keep_template", type: "boolean", named: true }],
      run: ({ id, model, keep_template }) => m("parts", d => { const p = pt(d, id); p.model = model; if (!keep_template) p.template = "custom"; p.overrides = {}; return { message: `${p.name} now uses a custom feature tree.` }; }) },
    { op: "parts.remove", title: "Remove a part", args: [{ name: "id", type: "string", required: true }], run: ({ id }) => m("parts", d => { pt(d, id); d.parts = d.parts.filter(p => p.id !== id); return { message: `Removed ${id}.` }; }) },
  ]);
  // Flight
  register([
    { op: "flight.set", title: "Set vehicle and mission parameters", example: "flight.set altitude=120 design_tw=2.2", args: [
      ...[["prop", 50, 2000, "mm"], ["fm", .3, .85, ""], ["motor_efficiency", .3, .98, ""], ["esc_efficiency", .5, .995, ""], ["design_tw", 1.1, 5, ""], ["avionics", 0, 5000, "W"],
        ["temperature", -40, 55, "°C"], ["pressure", 5e4, 1.1e5, "Pa"], ["altitude", 1, 5000, "m"], ["required_time", 0, 600, "min"], ["min_tw", 1, 5, ""], ["capture", 0, 1, ""]].map(([name, min, max, unit]) => ({ name, type: "number", min, max, unit, named: true })),
      { name: "arm", type: "number", min: 20, max: 5000, unit: "mm", named: true }, { name: "rotors", type: "enum", values: ["3", "4", "6", "8"], named: true }, { name: "name", type: "string", maxLength: 60, named: true }],
      run: a => m("flight", d => { const keys = { arm: "arm", prop: "prop", fm: "fm", motor_efficiency: "motorEfficiency", esc_efficiency: "escEfficiency", design_tw: "designTW", avionics: "avionicsPower", temperature: "temperature", pressure: "pressure", altitude: "altitude", required_time: "requiredTime", min_tw: "minTW", capture: "capture", name: "name" };
        for (const [k, key] of Object.entries(keys)) if (a[k] !== undefined) d[key] = a[k]; if (a.rotors) d.rotors = Number(a.rotors); return { message: "Vehicle updated." }; }) },
    { op: "flight.propulsion", title: "Use a measured motor and propeller table", summary: "Rows of throttle (0–1), volts, amps, watts, rpm, torque (N·m) and grams, from a thrust stand or a manufacturer. The vehicle's power, speed and limits then come from the data.", example: 'flight.propulsion "MN5008 KV170" "P18x6.1" kv=170 resistance=0.27 peak_current=15 max_power=720 rotor_inertia=0.0007 table=[…]',
      args: [{ name: "motor", type: "string", maxLength: 80, required: true }, { name: "propeller", type: "string", maxLength: 80, required: true }, ...[["kv", 1, 1e5], ["resistance", 0, 100], ["peak_current", 0, 1e4], ["max_power", 0, 1e6], ["rotor_inertia", 0, 10]].map(([name, min, max]) => ({ name, type: "number", min, max, named: true, required: true })), { name: "density", type: "number", min: .5, max: 1.5, unit: "kg/m³", named: true, default: 1.2 }, { name: "blockage", type: "number", min: 0, max: .5, named: true, default: 0, help: "Fraction of thrust lost to the arm" }, { name: "source", type: "string", maxLength: 300, named: true }, { name: "table", type: "json", named: true, required: true }],
      run: a => m("flight", d => { d.propulsion = validatePropulsion({ motor: a.motor, propeller: a.propeller, kv: a.kv, resistance: a.resistance, peakCurrent: a.peak_current, maxPower: a.max_power, rotorInertia: a.rotor_inertia, density: a.density, blockage: a.blockage, ...(a.source && { source: a.source }), table: a.table }); return { message: `Using ${a.table.length} measured points for ${a.motor} with ${a.propeller}.` }; }) },
    { op: "flight.clear_propulsion", title: "Return to momentum theory", run: () => m("flight", d => { d.propulsion = null; return { message: "Momentum theory with a stated figure of merit." }; }) },
    { op: "flight.mass", title: "Set or add a mass budget item", example: "flight.mass gimbal 1.4", args: [{ name: "id", type: "string", required: true }, { name: "mass", type: "number", min: 0, max: 500, unit: "kg", required: true }, { name: "name", type: "string", maxLength: 60, named: true }, { name: "part", type: "string", maxLength: 80, named: true, help: "Maker and part number" }, { name: "source", type: "string", maxLength: 300, named: true, help: "Where the mass comes from" }],
      run: ({ id, mass, name, part, source }) => m("flight", d => { if (!/^[a-zA-Z][a-zA-Z0-9_-]{0,23}$/.test(id)) throw new Error("Item ids start with a letter."); const item = d.masses.find(x => x.id === id), extra = { ...(part && { part }), ...(source && { source }) };
        if (item) { item.mass = mass; if (name) item.name = name; Object.assign(item, extra); } else { if (d.masses.length >= 24) throw new Error("The budget holds 24 items."); d.masses.push({ id, name: name || id, mass, ...extra }); }
        return { message: `${item?.name || name || id}: ${mass} kg.` }; }) },
    { op: "flight.remove_mass", title: "Remove a mass budget item", args: [{ name: "id", type: "string", required: true }],
      run: ({ id }) => m("flight", d => { if (!d.masses.some(x => x.id === id)) throw new Error(`No item ${id}.`); d.masses = d.masses.filter(x => x.id !== id); return { message: `Removed ${id}.` }; }) },
    { op: "flight.battery", title: "Set the battery pack", example: "flight.battery cells=12 capacity=22 mass=5.2", args: [{ name: "cells", type: "number", min: 1, max: 24, integer: true, named: true },
      ...[["capacity", .1, 200, "Ah"], ["mass", .01, 100, "kg"], ["resistance", 0, 1, "Ω"], ["usable", .3, 1, ""]].map(([name, min, max, unit]) => ({ name, type: "number", min, max, unit, named: true }))],
      run: a => m("flight", d => { for (const k of ["cells", "capacity", "mass", "resistance", "usable"]) if (a[k] !== undefined) d.battery[k] = a[k]; return { message: "Battery updated." }; }) },
  ]);
  // Power, structures, controls and testing
  const numbers = (fields, rename = {}) => Object.entries(fields).map(([key, [min, max]]) => ({ name: rename[key] ?? key, type: "number", min, max, named: true }));
  const snake = { rds25: "rds25", rdsTc: "rds_tc", tjMax: "tj_max", rthJc: "rth_jc" };
  const STAGE = { pwm: [1e3, 1e6], deadtime: [0, 1e4], gateDrive: [4, 20], gateResistor: [0, 100], driverPeak: [.1, 100], ringing: [1, 3], efficiencyLimit: [0, 100], phase: [0, 1e4], bus: [1, 2000], duty: [0, 100], junction: [-40, 200], peakBus: [1, 2000] }, stageNames = { gateDrive: "gate_drive", gateResistor: "gate_resistor", driverPeak: "driver_peak", efficiencyLimit: "efficiency_limit", peakBus: "peak_bus" };
  const apply = (target, fields, names, a) => { for (const k of Object.keys(fields)) { const v = a[names[k] ?? k]; if (v !== undefined) target[k] = v; } };
  const plantArgs = { mass: [.01, 1e4], ixx: [1e-6, 1e6], iyy: [1e-6, 1e6], izz: [1e-6, 1e6], arm: [10, 1e5], rotors: [4, 8], hoverThrottle: [1, 99], maxThrust: [.01, 1e5], tau: [1e-4, 5], kappa: [0, 5], dragArea: [0, 100] }, plantNames = { hoverThrottle: "hover_throttle", maxThrust: "max_thrust", dragArea: "drag_area" };
  const controlArgs = { loopRate: [50, 8000], gyroNoise: [0, 1], gust: [0, 50], wind: [0, 50], stepAngle: [1, 45] }, controlNames = { loopRate: "loop_rate", gyroNoise: "gyro_noise", stepAngle: "step_angle" };
  const test = (d, id) => { const t = d.plan.find(t => t.id === id); if (!t) throw new Error(`No test ${id}. Tests: ${d.plan.map(t => t.id).join(", ") || "none"}.`); return t; };
  register([
    { op: "power.set", title: "Set the power stage's operating point and gate drive", example: "power.set pwm=24000 deadtime=100 gate_resistor=3", args: [...numbers(STAGE, stageNames), { name: "name", type: "string", maxLength: 60, named: true }],
      run: a => m("power", d => { apply(d, STAGE, stageNames, a); if (a.name) d.name = a.name; return { message: "Power stage updated." }; }) },
    { op: "power.mosfet", title: "Set the MOSFET from its datasheet", summary: "RDS(on) in mΩ at 25 °C, charges in nC, plateau and diode in V, gate resistance in Ω.", example: 'power.mosfet "Infineon BSC040N10NS5" vds=100 rds25=3.4 qg=58 qsw=19', args: [{ name: "part", type: "string", maxLength: 80 }, ...numbers(MOSFET_FIELDS, snake), { name: "source", type: "string", maxLength: 300, named: true }],
      run: a => m("power", d => { if (a.part) d.mosfet.part = a.part; if (a.source) d.mosfet.source = a.source; apply(d.mosfet, MOSFET_FIELDS, snake, a); return { message: `MOSFET: ${d.mosfet.part}.` }; }) },
    { op: "power.capacitor", title: "Set the bus capacitors from their datasheet", summary: "A bank of count parts: capacitance in µF, ESR in mΩ, ripple rating in A, rated voltage in V, grams each, and the ripple limit as a percentage of the bus.", example: 'power.capacitor "Rubycon 63ZLJ470M12.5X20" count=3 capacitance=470 esr=33 rating=2.4 voltage=63', args: [{ name: "part", type: "string", maxLength: 80 }, ...numbers(CAPACITOR_FIELDS, {}), { name: "source", type: "string", maxLength: 300, named: true }],
      run: a => m("power", d => { if (a.part) d.capacitor.part = a.part; if (a.source) d.capacitor.source = a.source; apply(d.capacitor, CAPACITOR_FIELDS, {}, a); return { message: `Bus capacitors: ${d.capacitor.count} × ${d.capacitor.part}.` }; }) },
    { op: "structures.set", title: "Set the arm, its material and its limits", example: "structures.set od=30 wall=2 modulus=90", args: [...numbers(STRUCTURE_FIELDS, { safetyFactor: "safety_factor", deflectionLimit: "deflection_limit", tipMass: "tip_mass", loadThrust: "load_thrust", loadTorque: "load_torque" }), { name: "name", type: "string", maxLength: 60, named: true }, { name: "source", type: "string", maxLength: 300, named: true }],
      run: a => m("structures", d => { apply(d, STRUCTURE_FIELDS, { safetyFactor: "safety_factor", deflectionLimit: "deflection_limit", tipMass: "tip_mass", loadThrust: "load_thrust", loadTorque: "load_torque" }, a); if (a.name) d.name = a.name; if (a.source) d.source = a.source; return { message: "Structure updated." }; }) },
    { op: "controls.set", title: "Set the flight-test conditions", example: "controls.set gust=10 step_angle=20", args: numbers(controlArgs, controlNames), run: a => m("controls", d => { apply(d, controlArgs, controlNames, a); return { message: "Flight-test conditions updated." }; }) },
    { op: "controls.plant", title: "Set the plant, where nothing links it", args: numbers(plantArgs, plantNames), run: a => m("controls", d => { apply(d.plant, plantArgs, plantNames, a); return { message: "Plant updated." }; }) },
    { op: "controls.gain", title: "Set a PX4 parameter", example: "controls.gain MC_ROLLRATE_P 0.18", args: [{ name: "name", type: "enum", values: Object.keys(PX4_DEFAULTS), required: true }, { name: "value", type: "number", min: 0, max: 1e4, required: true }], run: ({ name, value }) => m("controls", d => { d.params[name] = value; return { message: `${name} = ${value}.` }; }) },
    { op: "controls.limit", title: "Set a control requirement", args: [{ name: "name", type: "enum", values: ["phaseMargin", "gainMargin", "overshoot", "settle", "drift"], required: true }, { name: "value", type: "number", min: 0, max: 1e4, required: true }], run: ({ name, value }) => m("controls", d => { d.limits[name] = value; return { message: `${name} limit ${value}.` }; }) },
    { op: "controls.tune", title: "Tune the rate loops for this airframe", summary: "The fastest rate loop with at least 55° of phase margin and 8 dB of gain margin.", run: () => m("controls", d => { const t = tune(d, lastPlant ?? plantOf(d)); if (!t) throw new Error("No gains in PX4's range meet the margins for this airframe."); const { bandwidth, ...gains } = t; Object.assign(d.params, gains); return { message: `Tuned: rate P ${gains.MC_ROLLRATE_P}, I ${gains.MC_ROLLRATE_I}, D ${gains.MC_ROLLRATE_D}, about ${bandwidth.toFixed(1)} Hz.` }; }) },
    { op: "testing.add", title: "Add a test to the plan", summary: "simulation reads an output against limits; correlation compares a dataset with the model; software and physical wait for a recorded result. covers lists the documents a recorded result depends on, comma separated: if one changes the result goes stale.", example: 'testing.add hover "Hover throttle" simulation port=flight.hover_throttle max=65', args: [{ name: "id", type: "string", maxLength: 32, required: true }, { name: "title", type: "string", maxLength: 120, required: true }, { name: "kind", type: "enum", values: TEST_KINDS, required: true }, { name: "method", type: "string", maxLength: 400, named: true }, { name: "port", type: "string", maxLength: 80, named: true }, { name: "dataset", type: "string", maxLength: 32, named: true }, { name: "covers", type: "string", maxLength: 120, named: true }, { name: "min", type: "number", named: true }, { name: "max", type: "number", named: true }, { name: "unit", type: "string", maxLength: 12, named: true }],
      run: a => m("testing", d => { if (d.plan.some(t => t.id === a.id)) throw new Error(`${a.id} already exists.`); if (a.kind === "simulation" && !a.port) throw new Error("A simulation test reads an output port."); d.plan.push(Object.fromEntries(Object.entries({ ...a, covers: a.covers ? a.covers.split(",").map(k => k.trim()).filter(Boolean) : undefined }).filter(([, v]) => v !== undefined))); return { message: `Added ${a.id}.` }; }) },
    { op: "testing.record", title: "Record a test result", summary: "A measurement on real hardware, or a real tool's result. Only this makes a physical or software test pass or fail.", example: 'testing.record thrust-stand 2.1 date=2026-11-02 by="A. Engineer"', args: [{ name: "id", type: "string", required: true }, { name: "value", type: "number", required: true }, { name: "date", type: "string", maxLength: 24, named: true }, { name: "by", type: "string", maxLength: 60, named: true }, { name: "note", type: "string", maxLength: 300, named: true }],
      run: ({ id, value, date, by, note }) => m("testing", d => { const t = test(d, id); if (t.kind !== "physical" && t.kind !== "software") throw new Error("Only physical and software tests take a recorded result."); t.result = { value, ...(date && { date }), ...(by && { by }), ...(note && { note }), ...(t.covers?.length && { digest: digest(store.docs, t.covers) }) }; return { message: `Recorded ${id}.` }; }) },
    { op: "testing.clear_result", title: "Clear a recorded result", args: [{ name: "id", type: "string", required: true }], run: ({ id }) => m("testing", d => { delete test(d, id).result; return { message: `${id} is not run.` }; }) },
    { op: "testing.remove", title: "Remove a test", args: [{ name: "id", type: "string", required: true }], run: ({ id }) => m("testing", d => { test(d, id); d.plan = d.plan.filter(t => t.id !== id); return { message: `Removed ${id}.` }; }) },
    { op: "testing.dataset", title: "Add measured data", summary: "Rows of rpm and grams, with torque if measured: yours, or a published set the model was not fitted to.", args: [{ name: "id", type: "string", maxLength: 32, required: true }, { name: "title", type: "string", maxLength: 120, required: true }, { name: "rows", type: "json", required: true }, { name: "source", type: "string", maxLength: 300, named: true }],
      run: ({ id, title, rows, source }) => m("testing", d => { if (d.datasets.some(s => s.id === id)) throw new Error(`${id} already exists.`); d.datasets.push({ id, title, rows, ...(source && { source }) }); return { message: `${rows.length} measured rows in ${id}.` }; }) },
    { op: "testing.trace", title: "Add a measured trace against the model", summary: "Rows of [time, measured, model]: a real run beside what the model predicted for it.", args: [{ name: "id", type: "string", maxLength: 32, required: true }, { name: "title", type: "string", maxLength: 120, required: true }, { name: "rows", type: "json", required: true }, { name: "unit", type: "string", maxLength: 12, named: true }, { name: "note", type: "string", maxLength: 300, named: true }, { name: "source", type: "string", maxLength: 300, named: true }],
      run: ({ id, title, rows, unit, note, source }) => m("testing", d => { if ((d.traces ?? []).some(s => s.id === id)) throw new Error(`A trace ${id} exists.`); (d.traces ??= []).push({ id, title, ...(unit && { unit }), ...(note && { note }), ...(source && { source }), rows }); return { message: `Added ${title}.` }; }) },
    { op: "testing.remove_trace", title: "Remove a trace", args: [{ name: "id", type: "string", required: true }], run: ({ id }) => m("testing", d => { if (!(d.traces ?? []).some(s => s.id === id)) throw new Error(`No trace ${id}.`); d.traces = d.traces.filter(s => s.id !== id); return { message: `Removed ${id}.` }; }) },
    { op: "testing.remove_dataset", title: "Remove measured data", args: [{ name: "id", type: "string", required: true }], run: ({ id }) => m("testing", d => { if (d.plan.some(t => t.dataset === id)) throw new Error(`A test uses ${id}.`); if (!d.datasets.some(s => s.id === id)) throw new Error(`No dataset ${id}.`); d.datasets = d.datasets.filter(s => s.id !== id); return { message: `Removed ${id}.` }; }) },
  ]);
  // Equations
  const unitArg = { name: "unit", type: "string", maxLength: 12, named: true, default: "" }, labelArg = { name: "label", type: "string", maxLength: 60, named: true };
  const put = (list, item) => { const i = list.findIndex(x => x.id === item.id); if (i >= 0) list[i] = item; else list.push(item); };
  register([
    { op: "equations.input", title: "Add or set a named input", summary: "A value of your own that a link can also drive.", example: "equations.input rail_voltage 0.9 unit=V", args: [{ name: "id", type: "string", maxLength: 32, required: true }, { name: "value", type: "number", required: true }, unitArg, labelArg, { name: "min", type: "number", named: true }, { name: "max", type: "number", named: true }],
      run: ({ id, value, unit, label, min, max }) => m("equations", d => { put(d.inputs, { id, value, unit, ...(label && { label }), ...(min !== undefined && { min }), ...(max !== undefined && { max }) }); return { message: `${id} = ${value}.` }; }) },
    { op: "equations.output", title: "Add or set a formula", summary: "An expression over inputs and earlier outputs, with + - * / ^ and functions such as sqrt, min, max, exp, log, sin.", example: 'equations.output rail_current "power / rail_voltage" unit=A', args: [{ name: "id", type: "string", maxLength: 32, required: true }, { name: "expr", type: "string", maxLength: 400, required: true }, unitArg, labelArg],
      run: ({ id, expr, unit, label }) => m("equations", d => { put(d.outputs, { id, expr, unit, ...(label && { label }) }); return { message: `${id} = ${expr}.` }; }) },
    { op: "equations.remove", title: "Remove an input or formula", args: [{ name: "id", type: "string", required: true }],
      run: ({ id }) => m("equations", d => { if (![...d.inputs, ...d.outputs].some(x => x.id === id)) throw new Error(`No ${id}.`); d.inputs = d.inputs.filter(x => x.id !== id); d.outputs = d.outputs.filter(x => x.id !== id); return { message: `Removed ${id}.` }; }) },
  ]);
  // System links and confirmations
  register([
    { op: "system.link", title: "Link an output to an input", example: "system.link optics.power_on_target mech.heat", args: [{ name: "from", type: "string", required: true }, { name: "to", type: "string", required: true }, { name: "note", type: "string", maxLength: 160, named: true }, { name: "id", type: "string", maxLength: 48, named: true }],
      run: a => m("system", d => { const link = { id: a.id || `link-${crypto.randomUUID().slice(0, 8)}`, from: a.from, to: a.to, enabled: true, ...(a.note && { note: a.note }) }, problem = system.check(link);
        if (problem) throw new Error(problem); if (d.links.some(l => l.id === link.id || l.from === a.from && l.to === a.to)) throw new Error("That link already exists."); d.links.push(link); return { id: link.id, message: `Linked ${a.from} → ${a.to}.` }; }) },
    { op: "system.toggle", title: "Enable or disable a link", args: [{ name: "id", type: "string", required: true }, { name: "enabled", type: "boolean" }],
      run: ({ id, enabled }) => m("system", d => { const l = d.links.find(l => l.id === id); if (!l) throw new Error(`No link ${id}.`); l.enabled = enabled ?? !l.enabled; return { message: `${id} ${l.enabled ? "enabled" : "disabled"}.` }; }) },
    { op: "system.unlink", title: "Remove a link", args: [{ name: "id", type: "string", required: true }], run: ({ id }) => m("system", d => { if (!d.links.some(l => l.id === id)) throw new Error(`No link ${id}.`); d.links = d.links.filter(l => l.id !== id); return { message: `Removed ${id}.` }; }) },
    { op: "system.confirm", title: "Confirm or reopen an assumption", summary: "Records that you have checked a nominal or illustrative input.", args: [{ name: "id", type: "string", required: true }, { name: "confirmed", type: "boolean", default: true }],
      run: ({ id, confirmed }) => m("system", d => { if (confirmed) d.confirmed[id] = true; else delete d.confirmed[id]; return { message: `${id} ${confirmed ? "confirmed" : "reopened"}.` }; }) },
  ]);
  // The project manifest: which domains are in use and what to watch.
  const domainArg = { name: "domain", type: "enum", values: DOMAINS, required: true }, startArg = { name: "start", type: "enum", values: ["empty", "starter"], named: true };
  const touchesStart = a => a.start ? ["project", a.domain, ...(a.domain === "silicon" ? ["litho"] : [])] : ["project"];
  const startDomain = (domain, start) => {
    if (domain === "mech") return mech.set({ params: values({}), features: [], drawings: {} });
    const doc = id => start === "empty" ? EMPTY[id]() : store.domains[id].example();
    store.replace(domain, doc(domain), "edit"); if (domain === "silicon") store.replace("litho", doc("litho"), "edit");
  };
  const outputPort = port => { const all = system.ports().outputs; if (!all.has(port)) throw new Error(`No output ${port} in this project. Outputs include ${[...all.keys()].slice(0, 6).join(", ") || "none yet"}.`); };
  const numericField = (domain, path) => { if (!Number.isFinite(getPath(domain === "mech" ? mech.state() : store.docs[domain], path))) throw new Error(`${path} is not a numeric field of the ${DOMAIN_INFO[domain].title} document.`); };
  const portArg = name => ({ name, type: "string", required: true, help: "An output port such as flight.total_mass" });
  register([
    { op: "project.rename", title: "Rename the project", example: 'project.rename "Survey drone"', args: [{ name: "name", type: "string", maxLength: 80, required: true }],
      run: ({ name }) => m("project", d => { if (!name.trim()) throw new Error("Enter a name."); d.name = name.trim(); return { message: `Project renamed ${d.name}.` }; }) },
    { op: "project.describe", title: "Describe the project", args: [{ name: "summary", type: "string", maxLength: 600, required: true }],
      run: ({ summary }) => m("project", d => { d.summary = summary.trim(); return { message: "Summary updated." }; }) },
    { op: "project.add_domain", title: "Add a domain to the project", summary: "start=empty begins with nothing in it, start=starter with an editable example; without start, a domain added back keeps its design.", example: "project.add_domain optics start=empty", args: [domainArg, startArg], touches: touchesStart,
      run: ({ domain, start }) => { if (active().includes(domain)) throw new Error(`${DOMAIN_INFO[domain].title} is already in this project.`); const r = m("project", d => { d.domains = DOMAINS.filter(x => x === domain || d.domains.includes(x)); return { message: `Added ${DOMAIN_INFO[domain].title}${start ? ` (${start})` : ""}.` }; }); if (start) startDomain(domain, start); return r; } },
    { op: "project.reset_domain", title: "Empty a domain or restore its starter", example: "project.reset_domain pcb start=empty", args: [domainArg, { ...startArg, required: true }], touches: touchesStart,
      run: ({ domain, start }) => { if (!active().includes(domain)) throw new Error(`${DOMAIN_INFO[domain].title} is not in this project.`); startDomain(domain, start); return { message: `${DOMAIN_INFO[domain].title}: ${start === "empty" ? "emptied" : "starter restored"}.` }; } },
    { op: "project.remove_domain", title: "Remove a domain from the project", summary: "Its document is kept, and links to it go dormant until it is added back.", args: [domainArg],
      run: ({ domain }) => m("project", d => { if (!d.domains.includes(domain)) throw new Error(`${DOMAIN_INFO[domain].title} is not in this project.`); d.domains = d.domains.filter(x => x !== domain); return { message: `Removed ${DOMAIN_INFO[domain].title}; its design and links are kept.` }; }) },
    { op: "project.pin", title: "Headline an output", example: 'project.pin optics.spot_rms label="Spot size"', args: [portArg("port"), { name: "label", type: "string", maxLength: 60, named: true }],
      run: ({ port, label }) => m("project", d => { outputPort(port); if (d.headline.some(h => h.port === port)) throw new Error(`${port} is already headlined.`); if (d.headline.length >= 12) throw new Error("Headline up to 12 outputs.");
        d.headline.push({ port, ...(label && { label }) }); return { message: `Headlined ${port}.` }; }) },
    { op: "project.unpin", title: "Remove a headline output", args: [portArg("port")],
      run: ({ port }) => m("project", d => { if (!d.headline.some(h => h.port === port)) throw new Error(`${port} is not headlined.`); d.headline = d.headline.filter(h => h.port !== port); return { message: `Removed ${port} from the headline.` }; }) },
    { op: "project.study", title: "Add a sensitivity input", summary: "Any numeric field of a domain's document, by path: object keys or list items by id.", example: 'project.study flight battery.capacity label="Battery capacity"',
      args: [domainArg, { name: "path", type: "string", maxLength: 120, required: true }, { name: "label", type: "string", maxLength: 60, named: true }, { name: "id", type: "string", maxLength: 32, named: true }, { name: "min", type: "number", named: true }, { name: "max", type: "number", named: true }],
      run: ({ domain, path, label, id, min, max }) => m("project", d => { numericField(domain, path);
        if (d.studies.length >= 16) throw new Error("Study up to 16 inputs."); const key = id || nextId(d.studies, "input"); if (d.studies.some(s => s.id === key)) throw new Error(`${key} already exists.`);
        d.studies.push({ id: key, label: label || `${DOMAIN_INFO[domain].title} ${path}`, domain, path, ...(min !== undefined && { min }), ...(max !== undefined && { max }) }); return { id: key, message: `Studying ${domain} ${path}.` }; }) },
    { op: "project.unstudy", title: "Remove a sensitivity input", args: [{ name: "id", type: "string", required: true }],
      run: ({ id }) => m("project", d => { if (!d.studies.some(s => s.id === id)) throw new Error(`No study input ${id}.`); d.studies = d.studies.filter(s => s.id !== id); return { message: `Removed ${id}.` }; }) },
    { op: "project.require", title: "Add a requirement on any output", summary: "Checked on every solve, next to the domains' own checks. An interval that straddles the limit reads undecided.", example: 'project.require mech.q1_tj max=90 label="Q1 stays below 90 °C"',
      args: [portArg("port"), { name: "max", type: "number", named: true }, { name: "min", type: "number", named: true }, { name: "label", type: "string", maxLength: 80, named: true }, { name: "id", type: "string", maxLength: 32, named: true }],
      run: ({ port, max, min, label, id }) => m("project", d => { outputPort(port); if (max === undefined && min === undefined) throw new Error("Give a max, a min or both."); if (d.requirements.length >= 32) throw new Error("A project holds 32 requirements."); const key = id || nextId(d.requirements, "req"); if (d.requirements.some(q => q.id === key)) throw new Error(`${key} already exists.`);
        d.requirements.push({ id: key, port, ...(min !== undefined && { min }), ...(max !== undefined && { max }), ...(label && { label }) }); return { id: key, message: `Requirement on ${port}.` }; }) },
    { op: "project.unrequire", title: "Remove a requirement", args: [{ name: "id", type: "string", required: true }], run: ({ id }) => m("project", d => { if (!d.requirements.some(q => q.id === id)) throw new Error(`No requirement ${id}.`); d.requirements = d.requirements.filter(q => q.id !== id); return { message: `Removed ${id}.` }; }) },
    { op: "project.optimize", title: "Define the project's optimisation", summary: "An objective output to minimise or maximise. Add variables and constraints next.", example: 'project.optimize "Lightest frame" flight.total_mass min',
      args: [{ name: "title", type: "string", maxLength: 80, required: true }, portArg("objective"), { name: "sense", type: "enum", values: ["min", "max"], default: "min" }, { name: "summary", type: "string", maxLength: 400, named: true }, { name: "trials", type: "number", min: 4, max: 200, integer: true, named: true }],
      run: ({ title, objective, sense, summary, trials }) => m("project", d => { outputPort(objective); d.optimize = { title, ...(summary && { summary }), objective: { port: objective, sense }, constraints: d.optimize?.constraints ?? [], variables: d.optimize?.variables ?? [], trials: trials ?? d.optimize?.trials ?? 40 }; return { message: `Optimisation: ${title}.` }; }) },
    { op: "project.optimize_variable", title: "Add an optimisation variable", summary: "Any numeric field, with its range and smallest step.", example: "project.optimize_variable mech params.fins 5 14 1 integer=true",
      args: [domainArg, { name: "path", type: "string", maxLength: 120, required: true }, { name: "min", type: "number", required: true }, { name: "max", type: "number", required: true }, { name: "step", type: "number", required: true }, { name: "label", type: "string", maxLength: 60, named: true }, { name: "integer", type: "boolean", named: true }],
      run: ({ domain, path, min, max, step, label, integer }) => m("project", d => { if (!d.optimize) throw new Error("Define the optimisation first with project.optimize."); numericField(domain, path); if (d.optimize.variables.length >= 6) throw new Error("Up to six variables."); d.optimize.variables.push({ domain, path, label: label || path, min, max, step, ...(integer && { integer: true }) }); return { message: `Variable ${path}.` }; }) },
    { op: "project.optimize_constraint", title: "Add an optimisation constraint", example: "project.optimize_constraint mech.q1_tj max=100",
      args: [portArg("port"), { name: "max", type: "number", named: true }, { name: "min", type: "number", named: true }, { name: "label", type: "string", maxLength: 80, named: true }],
      run: ({ port, max, min, label }) => m("project", d => { outputPort(port); if (max === undefined && min === undefined) throw new Error("Give a max, a min or both."); if (!d.optimize) throw new Error("Define the optimisation first with project.optimize."); if (d.optimize.constraints.length >= 8) throw new Error("Up to eight constraints."); d.optimize.constraints.push({ port, ...(min !== undefined && { min }), ...(max !== undefined && { max }), ...(label && { label }) }); return { message: `Constraint on ${port}.` }; }) },
    { op: "project.optimize_remove", title: "Remove an optimisation variable or constraint", args: [{ name: "kind", type: "enum", values: ["variable", "constraint"], required: true }, { name: "index", type: "number", min: 1, integer: true, required: true }],
      run: ({ kind, index }) => m("project", d => { const list = kind === "variable" ? d.optimize?.variables : d.optimize?.constraints; if (!list?.[index - 1]) throw new Error(`No ${kind} ${index}.`); list.splice(index - 1, 1); return { message: `Removed ${kind} ${index}.` }; }) },
    { op: "project.optimize_clear", title: "Remove the project's optimisation", run: () => m("project", d => { d.optimize = null; return { message: "Optimisation removed." }; }) },
    { op: "project.tour_step", title: "Add a step to the project's tour", summary: "A studio to show and a text; optionally a command to run, then the outputs (watch) and checks (domain/id) to report.", example: 'project.tour_step flight "Heavier camera" "Mass moves everything." command="flight.mass gimbal 1.6" watch=flight.total_mass,mech.q1_tj',
      args: [{ name: "studio", type: "enum", values: STUDIOS, required: true }, { name: "title", type: "string", maxLength: 80, required: true }, { name: "text", type: "string", maxLength: 600, required: true }, { name: "command", type: "string", maxLength: 400, named: true }, { name: "watch", type: "string", maxLength: 300, named: true }, { name: "checks", type: "string", maxLength: 200, named: true }],
      run: ({ studio, title, text, command, watch, checks }) => m("project", d => { if (d.tour.length >= 16) throw new Error("A tour holds 16 steps."); if (command) bus.parse(command); const list = s => s.split(",").map(x => x.trim()).filter(Boolean);
        d.tour.push({ studio, title, text, ...(command && { command }), ...(watch && { watch: list(watch) }), ...(checks && { checks: list(checks) }) }); return { message: `Tour step ${d.tour.length}.` }; }) },
    { op: "project.tour_remove", title: "Remove a tour step", args: [{ name: "index", type: "number", min: 1, integer: true, required: true }],
      run: ({ index }) => m("project", d => { if (!d.tour[index - 1]) throw new Error(`The tour has ${d.tour.length} steps.`); d.tour.splice(index - 1, 1); return { message: `Removed step ${index}.` }; }) },
    { op: "project.place", title: "Place a domain's geometry in the 3D scene", summary: "Position in mm, tilt and turn in degrees, and copies spaced about z. count may read a document, such as flight:rotors.", example: "project.place mech x=265 z=25 rz=90 phase=45 count=flight:rotors",
      args: [domainArg, ...["x", "y", "z"].map(name => ({ name, type: "number", min: -1e6, max: 1e6, unit: "mm", named: true, default: 0 })), ...["rx", "rz", "phase"].map(name => ({ name, type: "number", min: -360, max: 360, unit: "°", named: true })), { name: "count", type: "string", maxLength: 80, named: true }],
      run: ({ domain, x, y, z, rx, rz, phase, count }) => m("project", d => { d.scene.placements[domain] = { at: [x, y, z], ...(rx !== undefined && { rx }), ...(rz !== undefined && { rz }), ...(phase !== undefined && { phase }), ...(count !== undefined && { count: /^\d+$/.test(count) ? Number(count) : count }) }; return { message: `Placed ${DOMAIN_INFO[domain].title}.` }; }) },
    { op: "project.unplace", title: "Return a domain to automatic layout", args: [domainArg], run: ({ domain }) => m("project", d => { delete d.scene.placements[domain]; return { message: `${DOMAIN_INFO[domain].title} is laid out automatically.` }; }) },
    { op: "project.site", title: "Set where the design sits in the building", args: [{ name: "room", type: "string", maxLength: 32 }, { name: "z", type: "number", min: -1e5, max: 1e5, unit: "mm", named: true }],
      run: ({ room, z }) => m("project", d => { d.scene.site = { ...d.scene.site, ...(room && { room }), ...(z !== undefined && { z }) }; return { message: "Site placement updated." }; }) },
    { op: "project.credit", title: "Credit a public reference design", example: 'project.credit "Cooke triplet" https://example.org/lens domain=optics licence=MIT',
      args: [{ name: "title", type: "string", maxLength: 120, required: true }, { name: "url", type: "string", maxLength: 300, required: true }, { name: "domain", type: "enum", values: DOMAINS, named: true }, { name: "licence", type: "string", maxLength: 40, named: true }, { name: "note", type: "string", maxLength: 300, named: true }],
      run: ({ title, url, domain, licence, note }) => m("project", d => { d.credits.push({ title, url, ...(domain && { domain }), ...(licence && { licence }), ...(note && { note }) }); return { message: `Credited ${title}.` }; }) },
    { op: "project.uncredit", title: "Remove a credit", args: [{ name: "index", type: "number", min: 1, integer: true, required: true }],
      run: ({ index }) => m("project", d => { if (!d.credits[index - 1]) throw new Error(`The project has ${d.credits.length} credits.`); d.credits.splice(index - 1, 1); return { message: `Removed credit ${index}.` }; }) },
    { op: "project.set_value", title: "Set any numeric field by path", summary: "The domain's own validator accepts or refuses the result; mechanical parameters go through mech.set.", example: "project.set_value flight battery.capacity 18",
      args: [domainArg, { name: "path", type: "string", maxLength: 120, required: true }, { name: "value", type: "number", required: true }], touches: a => [a.domain],
      run: ({ domain, path, value }) => {
        if (domain === "mech") { if (!/^params\.[a-zA-Z]+$/.test(path)) throw new Error("Mechanical fields are params.<name>."); return bus.resolve("mech.set").run({ name: path.slice(7), value }); }
        return m(domain, d => { setPath(d, path, value); return { message: `${DOMAIN_INFO[domain].title} ${path} = ${value}.` }; });
      } },
  ]);
  const assumptions = () => {
    const docs = store.docs, confirmed = docs.system.confirmed, on = new Set(active());
    return [...(on.has("mech") ? mechanicalAssumptions() : []), ...(on.has("pcb") ? pcbAssumptions(docs.pcb) : []), ...(on.has("silicon") ? siliconAssumptions(docs.silicon) : []),
      ...(on.has("optics") ? opticsAssumptions(docs.optics) : []), ...(on.has("fluids") ? fluidsAssumptions(docs.fluids) : []), ...(on.has("architecture") ? architectureAssumptions(docs.architecture) : []),
      ...(on.has("parts") ? partsAssumptions(docs.parts) : []), ...(on.has("flight") ? flightAssumptions(docs.flight) : []), ...(on.has("equations") ? equationsAssumptions(docs.equations) : []), ...(on.has("power") ? powerAssumptions(docs.power) : []), ...(on.has("structures") ? structuresAssumptions(docs.structures) : []), ...(on.has("controls") ? controlsAssumptions(docs.controls) : [])]
      .map(a => ({ ...a, domain: a.id.split(".")[0], confirmed: !!confirmed[a.id] }));
  };
  // Open a project file: every document is validated by its domain before any
  // replaces the current one; domains the file omits start from their starter.
  const load = (file, reason = "load") => {
    const { documents } = readProjectFile(file), docs = {};
    for (const id of store.ids()) docs[id] = store.domains[id].validate(structuredClone(documents[id] ?? store.domains[id].example()));
    mech.set({ params: values(documents.mech?.params ?? {}), features: documents.mech?.features ?? [], drawings: documents.mech?.drawings ?? {} });
    store.load(docs, reason); bus.clearHistory();
  };
  const file = () => projectFile({ ...store.all(), mech: mech.snapshot() });
  return { store, bus, system, adapters, assumptions, active, load, file, rtlSources: RTL, solve: options => system.solve(store.docs.system, options) };
}
