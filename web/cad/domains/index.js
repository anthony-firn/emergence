// The design environment: every domain's documents, typed ports, evaluators,
// live certificates, assumptions and commands. What is being designed comes
// from a project file (core/project.js); the browser and the headless agent
// both build the environment here, with their own model transport and
// mechanical host, so every surface runs the same commands and physics.
import { CommandBus } from "../core/bus.js";
import { Store } from "../core/store.js";
import { System, SYSTEM_DEFAULT, validateSystem } from "../core/system.js";
import { DOMAINS, projectDefault, validateProject, readProjectFile, projectFile, getPath, setPath } from "../core/project.js";
import { FOOTPRINTS, pcbExample, siliconExample, documentFor, validateDocument, pcbChecks, siliconChecks, connectPads, moveComponent, clone, BLOCK_POWER } from "../studio-model.js";
import { kernelBoard, densityIssues, autoroute, pcbAssumptions, netOf } from "./pcb.js";
import { powerSummary, siliconAssumptions } from "./silicon.js";
import { lithoExample, validateLitho } from "./litho.js";
import { opticsExample, validateOptics, opticsResults, opticsAssumptions, paraxial, lensmaker, bestFocus, GLASSES, refractiveIndex, KINDS as OPTIC_KINDS } from "./optics.js";
import { fluidsExample, validateFluids, fluidsResults, fluidsAssumptions, FANS, fanPressure, KINDS as FLUID_KINDS } from "./fluids.js";
import { architectureExample, validateArchitecture, architectureResults, architectureAssumptions, overlaps, SIDES, AIR } from "./architecture.js";
import { partsExample, validateParts, partsAssumptions, TEMPLATES, MATERIALS } from "./parts.js";
import { mechanicalCommands, mechanicalAdapter, mechanicalAssumptions } from "./mechanical.js";
import { flightExample, validateFlight, flightResults, flightAssumptions } from "./flight.js";
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
};
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
  });
  const bus = new CommandBus();
  for (const id of store.ids()) bus.domain(id, { snapshot: () => store.snapshot(id), restore: s => store.restore(id, s) });
  bus.domain("mech", { snapshot: () => mech.snapshot(), restore: s => mech.restore(s) });
  const m = (id, fn) => store.mutate(id, fn);

  // Chip corner libraries are sent to the model worker once.
  let corners = null;
  const chip = async (d, rtl, temperature) => {
    if (!corners) corners = loadCorners().then(async data => { await request("chip", { input: { libs: data.corners, verilog: "module t(input a, output y); assign y = a; endmodule" } }); return data; });
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
    ["pcb", { title: "PCB", ports: () => ({
        inputs: [{ id: "u1_power", label: "U1 package dissipation", unit: "W", min: 0, max: 500 }, { id: "load_current", label: "Load net RMS current", unit: "A", min: 0, max: 500 }],
        outputs: [{ id: "copper_loss", label: "Copper I²R loss", unit: "W" }, { id: "dissipation", label: "Board + component dissipation", unit: "W" }, { id: "unrouted", label: "Unrouted connections", unit: "count" },
          { id: "board_mass", label: "Populated board mass", unit: "kg" }] }),
      async evaluate(d, ov) {
        d = clone(d);
        if (ov.u1_power !== undefined && d.components.some(c => c.id === "U1")) d.powers.U1 = ov.u1_power;
        if (ov.load_current !== undefined) for (const n of d.loadNets) d.currents[n] = ov.load_current;
        const geo = pcbChecks(d), density = densityIssues(d), k = kernelBoard(d);
        let kernel = null, error = null;
        if (d.routes.length) try { kernel = (await request("pcb", { input: { board: k.board, currents: k.currents, temperature: 40 } })).result; } catch (e) { error = e.message; }
        const copper = kernel?.copper_loss ?? 0, components = Object.values(d.powers).reduce((s, w) => s + w, 0);
        const checks = [check("geometry", "Copper geometry checks", !geo.issues.length, geo.issues.slice(0, 3).map(i => i.message).join(" ") + (geo.issues.length > 3 ? ` (+${geo.issues.length - 3} more)` : "")),
          check("density", `Track current density ≤ ${d.rules.density} A/mm²`, !density.length, density.map(i => i.message).join(" ")),
          check("routed", "Every net connected", geo.unrouted === 0 ? true : null, geo.unrouted ? `${geo.unrouted} connection${geo.unrouted > 1 ? "s" : ""} remain` : "")];
        if (error) checks.push({ id: "kernel", text: "Kernel board", status: "unavailable", reason: error });
        // FR-4 at 1.85 g/cm³ over the outline, drawn copper at 8.96 g/cm³, components by footprint.
        const t = d.stackup.copper * .0348, copperArea = d.routes.reduce((s, r) => s + r.points.slice(1).reduce((a, p, i) => a + Math.hypot(p[0] - r.points[i][0], p[1] - r.points[i][1]) * r.width, 0), 0);
        const mass = (d.width * d.height * d.stackup.thickness * 1.85e-3 + copperArea * t * 8.96e-3 + d.components.reduce((s, c) => s + FOOTPRINT_MASS[c.footprint], 0)) / 1000;
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
        const altitude = ov.altitude ?? 60, r = opticsResults(d, altitude), lens = d.elements.find(e => e.kind === "lens");
        const checks = [check("geometry", "Elements clear each other and the target", !r.issues.length, r.issues.map(i => i.message).join(" ")),
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
        const checks = [check("room", `Equipment room ≤ ${d.limit} °C`, room.temperature <= d.limit, null, room.temperature, "degC"),
          ...r.rooms.map((x, i) => check(`vent-${x.id}`, `${x.name}: outdoor air ≥ ASHRAE 62.1 office rate`, d.rooms[i].supply / 3.6 >= x.ventilation, `${(d.rooms[i].supply / 3.6).toFixed(0)} L/s supplied, ${x.ventilation.toFixed(0)} L/s required (if supply is outdoor air)`)),
          check("plan", "Rooms do not overlap", !clash.length, clash.map(([a, b]) => `${a.name} overlaps ${b.name}`).join("; "))];
        const balance = d.rooms.map((x, i) => { const z = r.rooms[i]; return x.mode === "supplied" ? Math.abs(AIR.rho * AIR.cp * x.supply / 3600 * (z.temperature - x.supplyTemp) - (z.internal + z.solar + z.UA * (d.site.outdoor - z.temperature))) : 0; });
        return { outputs: { room_temperature: { value: room.temperature, unit: "degC" }, cooling_load: { value: r.cooling, unit: "W" }, equipment_load: { value: r.equipmentW, unit: "W" }, floor_area: { value: r.floor, unit: "m2" }, outdoor_temperature: { value: d.site.outdoor, unit: "degC" } },
          checks, details: r, certificates: [{ label: "Supplied-room heat balances close", ok: balance.every(b => b < 1e-9) }] };
      } }],
    ["flight", { title: "Flight", ports: () => ({
        inputs: [{ id: "air_temperature", label: "Air temperature", unit: "degC", min: -40, max: 55 }, { id: "esc_mass", label: "Mass per ESC", unit: "kg", min: 0, max: 50, aggregate: "sum" },
          { id: "camera_mass", label: "Camera optics mass", unit: "kg", min: 0, max: 50, aggregate: "sum" }, { id: "parts_mass", label: "Airframe parts mass", unit: "kg", min: 0, max: 500, aggregate: "sum" }],
        outputs: [{ id: "total_mass", label: "Take-off mass", unit: "kg" }, { id: "hover_power", label: "Hover electrical power", unit: "W" }, { id: "esc_current", label: "Design-point ESC current", unit: "A" },
          { id: "hover_current", label: "Hover ESC current", unit: "A" }, { id: "flight_time", label: "Hover endurance", unit: "min" }, { id: "disk_loading", label: "Disk loading", unit: "Pa" },
          { id: "downwash", label: "Induced velocity at hover", unit: "m/s" }, { id: "wake_pressure", label: "Captured wake pressure", unit: "Pa" }, { id: "altitude", label: "Mission altitude", unit: "m" },
          { id: "charger_heat", label: "Charger heat per pack (1 h charge, 92%)", unit: "W" }] }),
      async evaluate(d, ov) {
        const r = flightResults(d, ov), h = r.hover, x = r.design;
        const energy = d.battery.capacity * r.V, charger = energy * (1 / .92 - 1);
        const checks = [check("endurance", `Hover endurance ≥ ${d.requiredTime} min`, r.time >= d.requiredTime, null, r.time, "min"),
          check("tw", `Design thrust-to-weight ${d.designTW} ≥ ${d.minTW} for control authority`, d.designTW >= d.minTW),
          check("pack", "Pack can deliver the design-point power", Number.isFinite(x.I), Number.isFinite(x.I) ? `${x.I.toFixed(1)} A at ${x.Vt.toFixed(1)} V` : "Power exceeds V²/4R: the pack voltage collapses"),
          check("crate", "Design current ≤ 10C continuous", x.I <= 10 * d.battery.capacity, null, x.I / d.battery.capacity, "")];
        return { outputs: { total_mass: { value: r.mass, unit: "kg" }, hover_power: { value: h.P, unit: "W" }, esc_current: { value: x.esc, unit: "A" }, hover_current: { value: h.esc, unit: "A" },
          flight_time: { value: r.time, unit: "min" }, disk_loading: { value: r.diskLoading, unit: "Pa" }, downwash: { value: h.vi, unit: "m/s" }, wake_pressure: { value: r.wakePressure, unit: "Pa" },
          altitude: { value: d.altitude, unit: "m" }, charger_heat: { value: charger, unit: "W" } }, checks, details: r,
          certificates: [{ label: "Hover power equals T·vᵢ/(FM·η) for every rotor", ok: Math.abs(h.input - h.T * h.vi / d.fm / d.motorEfficiency / d.escEfficiency) < 1e-9 * h.input }] };
      } }],
    ["parts", { title: "Parts", ports: () => ({ inputs: [], outputs: [{ id: "mass", label: "Parts mass", unit: "kg" }] }),
      // Volumes come from kernel meshes, cached by model; mass multiplies by count.
      async evaluate(d) {
        let mass = 0; const volumes = {};
        for (const p of d.parts) {
          const key = JSON.stringify([p.model, p.overrides]);
          if (!partVolumes.has(key)) partVolumes.set(key, mesh(p.model, p.overrides).then(r => r.volume));
          volumes[p.id] = await partVolumes.get(key);
          mass += volumes[p.id] * MATERIALS[p.material].density * 1e-6 * p.count;
        }
        return { outputs: { mass: { value: mass, unit: "kg" } }, checks: [], details: { volumes } };
      } }],
  ]);
  const partVolumes = new Map();
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
  register([
    { op: "pcb.place", title: "Place a footprint", example: "pcb.place dip8 @40,20 name=\"Gate driver\"", args: [{ name: "footprint", type: "enum", values: Object.keys(FOOTPRINTS), required: true },
      { name: "at", type: "point", help: "Board position in mm" }, { name: "name", type: "string", maxLength: 60, named: true }, { name: "ref", type: "string", maxLength: 12, named: true }],
      run: a => m("pcb", d => { if (d.components.length >= 64) throw new Error("This editor supports 64 components.");
        const ref = a.ref || nextId(d.components, { connector: "J", resistor: "R", capacitor: "C", transistor: "Q", dip8: "U" }[a.footprint]);
        if (d.components.some(c => c.id === ref)) throw new Error(`${ref} already exists.`);
        const [x, y] = a.at || [d.width / 2, d.height / 2], snap = v => Math.round(v / d.grid) * d.grid;
        d.components.push({ id: ref, name: a.name || FOOTPRINTS[a.footprint].title, footprint: a.footprint, x: snap(x), y: snap(y), rotation: 0, nets: FOOTPRINTS[a.footprint].pads.map(() => "") });
        return { id: ref, message: `Placed ${ref}.` }; }) },
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
    { op: "pcb.net", title: "Assign a pad's net", example: "pcb.net U1:4 GATE", args: [{ name: "pad", type: "string", required: true }, { name: "net", type: "string", maxLength: 32, required: true }],
      run: ({ pad, net }) => m("pcb", d => { const [id, n] = pad.split(":"), c = comp(d, id), i = Number(n) - 1;
        if (!Number.isInteger(i) || i < 0 || i >= c.nets.length) throw new Error(`${id} has pads 1–${c.nets.length}.`);
        if (!/^[a-zA-Z0-9_+./-]*$/.test(net)) throw new Error("Net names use letters, digits and _+./-.");
        c.nets[i] = net === "-" ? "" : net; d.routes = d.routes.filter(r => r.from !== pad && r.to !== pad); return { message: `${pad} → ${c.nets[i] || "no net"}.` }; }) },
    { op: "pcb.connect", title: "Connect two pads", summary: "Straight-and-45° connection between two pads on the same net.", example: "pcb.connect J1:1 U1:1 top 0.5",
      args: [{ name: "from", type: "string", required: true }, { name: "to", type: "string", required: true }, { name: "layer", type: "enum", values: ["top", "bottom"], default: "top" }, { name: "width", type: "number", min: .1, max: 5, unit: "mm", default: .5 }],
      run: a => m("pcb", d => { const r = connectPads(d, a.from, a.to, a.layer, a.width); d.routes.push(r); return { id: r.id, message: `Connected ${a.from} → ${a.to}.` }; }) },
    { op: "pcb.route", title: "Route around obstacles", summary: "A* on the snap grid, clear of other nets' copper by the clearance rule.", example: "pcb.route U1:4 R1:1 top 0.5",
      args: [{ name: "from", type: "string", required: true }, { name: "to", type: "string", required: true }, { name: "layer", type: "enum", values: ["top", "bottom"], default: "top" }, { name: "width", type: "number", min: .1, max: 5, unit: "mm", default: .5 }],
      run: a => m("pcb", d => { if (d.routes.length >= 128) throw new Error("This editor supports 128 connections."); const r = autoroute(d, a.from, a.to, a.layer, a.width); d.routes.push(r); return { id: r.id, message: `Routed ${a.from} → ${a.to} with ${r.points.length - 1} segments.` }; }) },
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
    { op: "silicon.add_block", title: "Add a floorplan block", example: "silicon.add_block logic \"Crypto engine\" @600,500 width=200 height=150", args: [{ name: "kind", type: "enum", values: ["logic", "memory", "io", "analog"], required: true },
      { name: "name", type: "string", maxLength: 60 }, { name: "at", type: "point" }, { name: "width", type: "number", min: 1, max: 10000, unit: "µm", named: true, default: 200 },
      { name: "height", type: "number", min: 1, max: 10000, unit: "µm", named: true, default: 160 }, { name: "power", type: "number", min: 0, max: 1e6, unit: "mW", named: true }],
      run: a => m("silicon", d => { if (d.blocks.length >= 64) throw new Error("This editor supports 64 blocks."); const id = nextId(d.blocks, a.kind), [x, y] = a.at || [d.margin, d.margin];
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
      .map(([name, min, max, unit]) => ({ name, type: "number", min, max, unit, named: true })), { name: "source_points", type: "number", min: 1, max: 9, integer: true, named: true }],
      run: a => m("litho", d => { for (const [k, key] of [["wavelength", "wavelength"], ["na", "na"], ["sigma_outer", "sigmaOuter"], ["sigma_inner", "sigmaInner"], ["threshold", "threshold"], ["pixel", "pixel"], ["source_points", "sourcePoints"]]) if (a[k] !== undefined) d[key] = a[k]; return { message: "Lithography settings updated." }; }) },
    { op: "litho.shape", title: "Add a rectangle to the mask clip", example: "litho.shape 200 200 160 600", args: ["x", "y", "w", "h"].map(name => ({ name, type: "number", min: 0, max: 1e5, unit: "nm", required: true })),
      run: a => m("litho", d => { d.shapes.push(a); return { message: "Added a mask rectangle." }; }) },
    { op: "litho.clear", title: "Clear the mask clip", run: () => m("litho", d => { d.shapes = []; return { message: "Cleared the clip." }; }) },
  ]);
  // Optics
  const element = (d, id) => { const e = d.elements.find(e => e.id === id); if (!e) throw new Error(`No element ${id}. Elements: ${d.elements.map(e => e.id).join(", ")}.`); return e; };
  register([
    { op: "optics.set", title: "Set the beam, wavelength and target", example: "optics.set power=3 wavelength=980 target_z=120", args: [
      ...[["wavelength", 300, 2500, "nm"], ["power", 0, 100, "W"], ["diameter", .1, 100, "mm"], ["target_z", -1000, 5000, "mm"], ["target_diameter", .01, 500, "mm"]].map(([name, min, max, unit]) => ({ name, type: "number", min, max, unit, named: true })),
      { name: "rays", type: "number", min: 16, max: 4000, integer: true, named: true }],
      run: a => m("optics", d => { if (a.wavelength !== undefined) d.wavelength = a.wavelength; if (a.power !== undefined) d.source.power = a.power; if (a.diameter !== undefined) d.source.diameter = a.diameter;
        if (a.target_z !== undefined) d.target.z = a.target_z; if (a.target_diameter !== undefined) d.target.diameter = a.target_diameter; if (a.rays !== undefined) d.rays = a.rays; return { message: "Optical settings updated." }; }) },
    { op: "optics.add", title: "Add a lens, window or stop", example: "optics.add lens z=70 r1=80 r2=-80 thickness=5 diameter=25 glass=N-BK7", args: [{ name: "kind", type: "enum", values: OPTIC_KINDS, required: true }, { name: "name", type: "string", maxLength: 60, named: true },
      ...[["z", -1000, 5000], ["r1", -1e5, 1e5], ["r2", -1e5, 1e5], ["thickness", .05, 200], ["diameter", .1, 500]].map(([name, min, max]) => ({ name, type: "number", min, max, unit: "mm", named: true })),
      { name: "glass", type: "enum", values: Object.keys(GLASSES), named: true, default: "N-BK7" }],
      run: a => m("optics", d => { const id = nextId(d.elements, { lens: "L", window: "W", stop: "S" }[a.kind]), last = Math.max(0, ...d.elements.map(e => e.z + (e.thickness || 0)));
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
    { op: "fluids.set", title: "Set air conditions and the rise limit", args: [{ name: "temperature", type: "number", min: -40, max: 150, unit: "°C", named: true }, { name: "pressure", type: "number", min: 5e4, max: 2e5, unit: "Pa", named: true }, { name: "heat_limit", type: "number", min: .1, max: 200, unit: "K", named: true }],
      run: a => m("fluids", d => { if (a.temperature !== undefined) d.temperature = a.temperature; if (a.pressure !== undefined) d.pressure = a.pressure; if (a.heat_limit !== undefined) d.heatLimit = a.heat_limit; return { message: "Air conditions updated." }; }) },
    { op: "fluids.add", title: "Add a duct, fitting or filter", example: "fluids.add duct name=Riser width=60 height=30 length=120", args: [{ name: "kind", type: "enum", values: FLUID_KINDS.filter(k => k !== "channels"), required: true }, { name: "name", type: "string", maxLength: 60, named: true },
      ...[["width", .5, 2000, "mm"], ["height", .5, 2000, "mm"], ["length", .1, 1e5, "mm"], ["roughness", 0, 10, "mm"], ["K", 0, 100, ""], ["coefficient", 0, 1e4, "Pa·s/m"]].map(([name, min, max, unit]) => ({ name, type: "number", min, max, unit, named: true })),
      { name: "after", type: "string", named: true }],
      run: a => m("fluids", d => { const id = nextId(d.elements, a.kind), e = { id, name: a.name || { duct: "Duct", fitting: "Fitting", filter: "Filter" }[a.kind], kind: a.kind, width: a.width ?? 60, height: a.height ?? 30 };
        if (a.kind === "duct") Object.assign(e, { length: a.length ?? 100, roughness: a.roughness ?? .0015 }); if (a.kind === "fitting") e.K = a.K ?? .5; if (a.kind === "filter") e.coefficient = a.coefficient ?? 25;
        const at = a.after ? d.elements.findIndex(x => x.id === a.after) : d.elements.findIndex(x => x.kind === "channels");
        d.elements.splice(at < 0 ? d.elements.length : a.after ? at + 1 : at, 0, e); return { id, message: `Added ${e.name} (${id}).` }; }) },
    { op: "fluids.update", title: "Edit an air-path element", example: "fluids.update duct length=240", args: [{ name: "id", type: "string", required: true }, { name: "name", type: "string", maxLength: 60, named: true },
      ...[["width", .5, 2000], ["height", .5, 2000], ["length", .1, 1e5], ["roughness", 0, 10], ["K", 0, 100], ["coefficient", 0, 1e4], ["gap", .05, 200], ["Kin", 0, 10], ["Kout", 0, 10]].map(([name, min, max]) => ({ name, type: "number", min, max, named: true })),
      { name: "fins", type: "number", min: 2, max: 200, integer: true, named: true }],
      run: ({ id, ...c }) => m("fluids", d => { const e = fel(d, id); for (const [k, v] of Object.entries(c)) if (v !== undefined) { if (!(k in e)) throw new Error(`${e.name} has no ${k}.`); e[k] = v; } return { message: `Updated ${e.name}.` }; }) },
    { op: "fluids.remove", title: "Remove an air-path element", args: [{ name: "id", type: "string", required: true }],
      run: ({ id }) => m("fluids", d => { const e = fel(d, id); if (e.kind === "channels") throw new Error("The heat-sink channels are linked to the mechanical design and stay in the path."); d.elements = d.elements.filter(x => x.id !== id); return { message: `Removed ${e.name}.` }; }) },
  ]);
  // Architecture
  const room = (d, id) => { const r = d.rooms.find(r => r.id === id); if (!r) throw new Error(`No room ${id}. Rooms: ${d.rooms.map(r => r.id).join(", ")}.`); return r; };
  const roomArgs = [...[["width", .5, 200, "m"], ["depth", .5, 200, "m"], ["height", 2, 20, "m"], ["people", 0, 2000, ""], ["lighting", 0, 50, "W/m²"], ["supply", 0, 1e6, "m³/h"], ["supply_temp", 5, 40, "°C"], ["setpoint", 10, 35, "°C"], ["roof_u", 0, 5, "W/m²K"]]
    .map(([name, min, max, unit]) => ({ name, type: "number", min, max, unit, named: true })), { name: "mode", type: "enum", values: ["supplied", "conditioned"], named: true }, { name: "equipment", type: "boolean", named: true }, { name: "name", type: "string", maxLength: 60, named: true }];
  const applyRoom = (d, r, a) => { for (const [k, key] of [["width", "width"], ["depth", "depth"], ["height", "height"], ["people", "people"], ["lighting", "lighting"], ["supply", "supply"], ["supply_temp", "supplyTemp"], ["setpoint", "setpoint"], ["roof_u", "roofU"], ["mode", "mode"], ["name", "name"]]) if (a[k] !== undefined) r[key] = a[k];
    if (a.at) [r.x, r.y] = a.at; if (a.equipment !== undefined) { for (const x of d.rooms) x.equipment = false; r.equipment = a.equipment; if (!d.rooms.some(x => x.equipment)) d.rooms[0].equipment = true; } };
  register([
    { op: "architecture.add_room", title: "Add a room", example: "architecture.add_room Workshop @13,0 width=6 depth=9", args: [{ name: "label", type: "string", maxLength: 60, required: true }, { name: "at", type: "point", help: "Plan position in m" }, ...roomArgs.filter(a => a.name !== "name")],
      run: a => m("architecture", d => { if (d.rooms.length >= 32) throw new Error("This plan supports 32 rooms."); const id = nextId(d.rooms, "room"), right = Math.max(...d.rooms.map(r => r.x + r.width));
        const w = (ext = true) => ({ exterior: ext, U: .35, window: 0 }), r = { id, name: a.label, x: right, y: 0, width: 5, depth: 5, height: 3, people: 2, lighting: 9, mode: "conditioned", supply: 200, supplyTemp: 16, setpoint: 24, equipment: false, roofU: .25, walls: { north: w(), east: w(), south: w(), west: w() } };
        applyRoom(d, r, a); d.rooms.push(r); return { id, message: `Added ${r.name} (${id}).` }; }) },
    { op: "architecture.update_room", title: "Edit a room", example: "architecture.update_room lab supply=1200 people=5", args: [{ name: "id", type: "string", required: true }, { name: "at", type: "point" }, ...roomArgs],
      run: ({ id, ...a }) => m("architecture", d => { const r = room(d, id); applyRoom(d, r, a); return { message: `Updated ${r.name}.` }; }) },
    { op: "architecture.wall", title: "Set a wall: exterior, U-value and window", example: "architecture.wall lab west window=4", args: [{ name: "id", type: "string", required: true }, { name: "side", type: "enum", values: SIDES, required: true },
      { name: "exterior", type: "boolean", named: true }, { name: "u", type: "number", min: .05, max: 6, unit: "W/m²K", named: true }, { name: "window", type: "number", min: 0, max: 5000, unit: "m²", named: true }],
      run: ({ id, side, exterior, u, window }) => m("architecture", d => { const w = room(d, id).walls[side]; if (exterior !== undefined) w.exterior = exterior; if (u !== undefined) w.U = u; if (window !== undefined) w.window = window; return { message: `${id} ${side} wall updated.` }; }) },
    { op: "architecture.remove_room", title: "Remove a room", args: [{ name: "id", type: "string", required: true }],
      run: ({ id }) => m("architecture", d => { room(d, id); if (d.rooms.length === 1) throw new Error("Keep at least one room."); d.rooms = d.rooms.filter(r => r.id !== id); if (!d.rooms.some(r => r.equipment)) d.rooms[0].equipment = true; return { message: `Removed ${id}.` }; }) },
    { op: "architecture.site", title: "Set climate, glazing, unit count and limit", example: "architecture.site outdoor=35 units=40", args: [...[["outdoor", -40, 55, "°C"], ["solar", 0, 1200, "W/m²"], ["glazing_u", .3, 7, "W/m²K"], ["shgc", .05, .9, ""], ["limit", 5, 80, "°C"]].map(([name, min, max, unit]) => ({ name, type: "number", min, max, unit, named: true })),
      { name: "units", type: "number", min: 0, max: 10000, integer: true, named: true }, { name: "name", type: "string", maxLength: 60, named: true }],
      run: a => m("architecture", d => { if (a.outdoor !== undefined) d.site.outdoor = a.outdoor; if (a.solar !== undefined) d.site.solar = a.solar; if (a.name !== undefined) d.site.name = a.name;
        if (a.glazing_u !== undefined) d.glazing.U = a.glazing_u; if (a.shgc !== undefined) d.glazing.shgc = a.shgc; if (a.units !== undefined) d.units = a.units; if (a.limit !== undefined) d.limit = a.limit; return { message: "Site updated." }; }) },
  ]);
  // Parts
  const pt = (d, id) => { const p = d.parts.find(p => p.id === id); if (!p) throw new Error(`No part ${id}. Parts: ${d.parts.map(p => p.id).join(", ") || "none"}.`); return p; };
  register([
    { op: "parts.add", title: "Add a solid from a template", example: "parts.add flange @0,-90,0", args: [{ name: "template", type: "enum", values: Object.keys(TEMPLATES), required: true }, { name: "name", type: "string", maxLength: 60, named: true },
      { name: "x", type: "number", min: -2000, max: 2000, unit: "mm", named: true, default: 0 }, { name: "y", type: "number", min: -2000, max: 2000, unit: "mm", named: true, default: 0 }, { name: "z", type: "number", min: -2000, max: 2000, unit: "mm", named: true, default: 0 },
      { name: "material", type: "enum", values: Object.keys(MATERIALS), named: true }],
      run: a => m("parts", d => { if (d.parts.length >= 24) throw new Error("This project supports 24 parts."); const id = nextId(d.parts, a.template), t = TEMPLATES[a.template];
        d.parts.push({ id, name: a.name || t.title, template: a.template, model: structuredClone(t.model), overrides: {}, x: a.x, y: a.y, z: a.z, material: a.material || t.material, domain: "Mechanical" }); return { id, message: `Added ${t.title} (${id}).` }; }) },
    { op: "parts.param", title: "Set a part's model parameter", example: "parts.param gear1 teeth 30", args: [{ name: "id", type: "string", required: true }, { name: "name", type: "string", required: true }, { name: "value", type: "number", required: true }],
      run: ({ id, name, value }) => m("parts", d => { const p = pt(d, id); if (!Object.hasOwn(p.model.params, name)) throw new Error(`${p.name} has parameters ${Object.keys(p.model.params).join(", ")}.`); p.overrides[name] = value; return { message: `${p.name}: ${name} = ${value}.` }; }) },
    { op: "parts.update", title: "Move a part or change its material", args: [{ name: "id", type: "string", required: true }, ...["x", "y", "z"].map(name => ({ name, type: "number", min: -2000, max: 2000, unit: "mm", named: true })),
      { name: "material", type: "enum", values: Object.keys(MATERIALS), named: true }, { name: "name", type: "string", maxLength: 60, named: true }],
      run: ({ id, ...c }) => m("parts", d => { const p = pt(d, id); for (const [k, v] of Object.entries(c)) if (v !== undefined) p[k] = v; return { message: `Updated ${p.name}.` }; }) },
    { op: "parts.model", title: "Replace a part's feature tree", summary: "JSON in the fluid cad model format: params, shape (ops and kits), optional cell.", args: [{ name: "id", type: "string", required: true }, { name: "model", type: "json", required: true }],
      run: ({ id, model }) => m("parts", d => { const p = pt(d, id); p.model = model; p.template = "custom"; p.overrides = {}; return { message: `${p.name} now uses a custom feature tree.` }; }) },
    { op: "parts.remove", title: "Remove a part", args: [{ name: "id", type: "string", required: true }], run: ({ id }) => m("parts", d => { pt(d, id); d.parts = d.parts.filter(p => p.id !== id); return { message: `Removed ${id}.` }; }) },
  ]);
  // Flight
  register([
    { op: "flight.set", title: "Set vehicle and mission parameters", example: "flight.set altitude=120 design_tw=2.2", args: [
      ...[["prop", 50, 2000, "mm"], ["fm", .3, .85, ""], ["motor_efficiency", .3, .98, ""], ["esc_efficiency", .5, .995, ""], ["design_tw", 1.1, 5, ""], ["avionics", 0, 5000, "W"],
        ["temperature", -40, 55, "°C"], ["altitude", 1, 5000, "m"], ["required_time", 0, 600, "min"], ["min_tw", 1, 5, ""], ["capture", 0, 1, ""]].map(([name, min, max, unit]) => ({ name, type: "number", min, max, unit, named: true })),
      { name: "rotors", type: "enum", values: ["3", "4", "6", "8"], named: true }, { name: "name", type: "string", maxLength: 60, named: true }],
      run: a => m("flight", d => { const keys = { prop: "prop", fm: "fm", motor_efficiency: "motorEfficiency", esc_efficiency: "escEfficiency", design_tw: "designTW", avionics: "avionicsPower", temperature: "temperature", altitude: "altitude", required_time: "requiredTime", min_tw: "minTW", capture: "capture", name: "name" };
        for (const [k, key] of Object.entries(keys)) if (a[k] !== undefined) d[key] = a[k]; if (a.rotors) d.rotors = Number(a.rotors); return { message: "Vehicle updated." }; }) },
    { op: "flight.mass", title: "Set or add a mass budget item", example: "flight.mass gimbal 1.4", args: [{ name: "id", type: "string", required: true }, { name: "mass", type: "number", min: 0, max: 500, unit: "kg", required: true }, { name: "name", type: "string", maxLength: 60, named: true }],
      run: ({ id, mass, name }) => m("flight", d => { if (!/^[a-zA-Z][a-zA-Z0-9_-]{0,23}$/.test(id)) throw new Error("Item ids start with a letter."); const item = d.masses.find(x => x.id === id);
        if (item) { item.mass = mass; if (name) item.name = name; } else { if (d.masses.length >= 24) throw new Error("The budget holds 24 items."); d.masses.push({ id, name: name || id, mass }); }
        return { message: `${item?.name || name || id}: ${mass} kg.` }; }) },
    { op: "flight.remove_mass", title: "Remove a mass budget item", args: [{ name: "id", type: "string", required: true }],
      run: ({ id }) => m("flight", d => { if (!d.masses.some(x => x.id === id)) throw new Error(`No item ${id}.`); d.masses = d.masses.filter(x => x.id !== id); return { message: `Removed ${id}.` }; }) },
    { op: "flight.battery", title: "Set the battery pack", example: "flight.battery cells=12 capacity=22 mass=5.2", args: [{ name: "cells", type: "number", min: 1, max: 24, integer: true, named: true },
      ...[["capacity", .1, 200, "Ah"], ["mass", .01, 100, "kg"], ["resistance", 0, 1, "Ω"], ["usable", .3, 1, ""]].map(([name, min, max, unit]) => ({ name, type: "number", min, max, unit, named: true }))],
      run: a => m("flight", d => { for (const k of ["cells", "capacity", "mass", "resistance", "usable"]) if (a[k] !== undefined) d.battery[k] = a[k]; return { message: "Battery updated." }; }) },
  ]);
  // System links and confirmations
  register([
    { op: "system.link", title: "Link an output to an input", example: "system.link optics.power_on_target mech.heat", args: [{ name: "from", type: "string", required: true }, { name: "to", type: "string", required: true }, { name: "note", type: "string", maxLength: 160, named: true }],
      run: a => m("system", d => { const link = { id: `link-${crypto.randomUUID().slice(0, 8)}`, from: a.from, to: a.to, enabled: true, ...(a.note && { note: a.note }) }, problem = system.check(link);
        if (problem) throw new Error(problem); if (d.links.some(l => l.from === a.from && l.to === a.to)) throw new Error("That link already exists."); d.links.push(link); return { id: link.id, message: `Linked ${a.from} → ${a.to}.` }; }) },
    { op: "system.toggle", title: "Enable or disable a link", args: [{ name: "id", type: "string", required: true }, { name: "enabled", type: "boolean" }],
      run: ({ id, enabled }) => m("system", d => { const l = d.links.find(l => l.id === id); if (!l) throw new Error(`No link ${id}.`); l.enabled = enabled ?? !l.enabled; return { message: `${id} ${l.enabled ? "enabled" : "disabled"}.` }; }) },
    { op: "system.unlink", title: "Remove a link", args: [{ name: "id", type: "string", required: true }], run: ({ id }) => m("system", d => { if (!d.links.some(l => l.id === id)) throw new Error(`No link ${id}.`); d.links = d.links.filter(l => l.id !== id); return { message: `Removed ${id}.` }; }) },
    { op: "system.confirm", title: "Confirm or reopen an assumption", summary: "Records that you have checked a nominal or illustrative input.", args: [{ name: "id", type: "string", required: true }, { name: "confirmed", type: "boolean", default: true }],
      run: ({ id, confirmed }) => m("system", d => { if (confirmed) d.confirmed[id] = true; else delete d.confirmed[id]; return { message: `${id} ${confirmed ? "confirmed" : "reopened"}.` }; }) },
  ]);
  // The project manifest: which domains are in use and what to watch.
  const domainArg = { name: "domain", type: "enum", values: DOMAINS, required: true };
  const portArg = name => ({ name, type: "string", required: true, help: "An output port such as flight.total_mass" });
  register([
    { op: "project.rename", title: "Rename the project", example: 'project.rename "Survey drone"', args: [{ name: "name", type: "string", maxLength: 80, required: true }],
      run: ({ name }) => m("project", d => { if (!name.trim()) throw new Error("Enter a name."); d.name = name.trim(); return { message: `Project renamed ${d.name}.` }; }) },
    { op: "project.describe", title: "Describe the project", args: [{ name: "summary", type: "string", maxLength: 600, required: true }],
      run: ({ summary }) => m("project", d => { d.summary = summary.trim(); return { message: "Summary updated." }; }) },
    { op: "project.add_domain", title: "Add a domain to the project", summary: "Opens its studio with a starter design; links to it can then be drawn in the System studio.", example: "project.add_domain optics", args: [domainArg],
      run: ({ domain }) => m("project", d => { if (d.domains.includes(domain)) throw new Error(`${DOMAIN_INFO[domain].title} is already in this project.`); d.domains = DOMAINS.filter(x => x === domain || d.domains.includes(x)); return { message: `Added ${DOMAIN_INFO[domain].title}.` }; }) },
    { op: "project.remove_domain", title: "Remove a domain from the project", summary: "Its document is kept, and links to it go dormant until it is added back.", args: [domainArg],
      run: ({ domain }) => m("project", d => { if (!d.domains.includes(domain)) throw new Error(`${DOMAIN_INFO[domain].title} is not in this project.`); d.domains = d.domains.filter(x => x !== domain); return { message: `Removed ${DOMAIN_INFO[domain].title}; its design and links are kept.` }; }) },
    { op: "project.pin", title: "Headline an output", example: 'project.pin optics.spot_rms label="Spot size"', args: [portArg("port"), { name: "label", type: "string", maxLength: 60, named: true }],
      run: ({ port, label }) => m("project", d => { if (!system.ports().outputs.has(port)) throw new Error(`No output ${port} in this project.`); if (d.headline.some(h => h.port === port)) throw new Error(`${port} is already headlined.`); if (d.headline.length >= 12) throw new Error("Headline up to 12 outputs.");
        d.headline.push({ port, ...(label && { label }) }); return { message: `Headlined ${port}.` }; }) },
    { op: "project.unpin", title: "Remove a headline output", args: [portArg("port")],
      run: ({ port }) => m("project", d => { if (!d.headline.some(h => h.port === port)) throw new Error(`${port} is not headlined.`); d.headline = d.headline.filter(h => h.port !== port); return { message: `Removed ${port} from the headline.` }; }) },
    { op: "project.study", title: "Add a sensitivity input", summary: "Any numeric field of a domain's document, by path: object keys or list items by id.", example: 'project.study flight battery.capacity label="Battery capacity"',
      args: [domainArg, { name: "path", type: "string", maxLength: 120, required: true }, { name: "label", type: "string", maxLength: 60, named: true }],
      run: ({ domain, path, label }) => m("project", d => { const doc = domain === "mech" ? mech.state() : store.docs[domain];
        if (!Number.isFinite(getPath(doc, path))) throw new Error(`${path} is not a numeric field of the ${DOMAIN_INFO[domain].title} document.`);
        if (d.studies.length >= 16) throw new Error("Study up to 16 inputs."); let i = 1; while (d.studies.some(s => s.id === `input${i}`)) i++;
        d.studies.push({ id: `input${i}`, label: label || `${DOMAIN_INFO[domain].title} ${path}`, domain, path }); return { id: `input${i}`, message: `Studying ${domain} ${path}.` }; }) },
    { op: "project.unstudy", title: "Remove a sensitivity input", args: [{ name: "id", type: "string", required: true }],
      run: ({ id }) => m("project", d => { if (!d.studies.some(s => s.id === id)) throw new Error(`No study input ${id}.`); d.studies = d.studies.filter(s => s.id !== id); return { message: `Removed ${id}.` }; }) },
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
      ...(on.has("parts") ? partsAssumptions(docs.parts) : []), ...(on.has("flight") ? flightAssumptions(docs.flight) : [])]
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
