// Multirotor sizing: mass budget, hover and design-point thrust, rotor power,
// battery current and endurance. With a measured propulsion table (a
// manufacturer's or your own thrust stand) power, speed and torque are looked
// up, corrected for air density by similarity; momentum theory then stays as a
// cross-check, as the figure of merit the data implies. Without one, ideal
// induced power T^1.5 / √(2ρA) is divided by a stated figure of merit and
// efficiencies. Endurance integrates a LiPo discharge curve under pack sag, and
// the maximum thrust follows the sagging voltage.
export const G = 9.80665;
// Open-circuit volts per cell against state of charge: the widely used LiPo rest
// voltage table. Typical of the chemistry, not of a particular cell.
const OCV = [[0, 3.27], [.05, 3.61], [.1, 3.69], [.15, 3.71], [.2, 3.73], [.25, 3.75], [.3, 3.77], [.35, 3.79], [.4, 3.8], [.45, 3.82], [.5, 3.84], [.55, 3.85], [.6, 3.87], [.65, 3.91], [.7, 3.95], [.75, 3.98], [.8, 4.02], [.85, 4.08], [.9, 4.11], [.95, 4.15], [1, 4.2]];
export function cellVoltage(soc) { const s = Math.min(1, Math.max(0, soc)), i = Math.min(OCV.length - 2, Math.floor(s / .05)), [x0, y0] = OCV[i], [x1, y1] = OCV[i + 1]; return y0 + (y1 - y0) * (s - x0) / (x1 - x0); }
// Current drawn for power P from a source of open-circuit V and resistance R.
const sag = (V, R, P) => { if (R <= 0) return P / V; const disc = V * V - 4 * R * P; return disc >= 0 ? (V - Math.sqrt(disc)) / (2 * R) : NaN; };
export function flightExample() {
  return {
    name: "Quadcopter", rotors: 4, prop: 381, fm: .65, motorEfficiency: .82, escEfficiency: .97, designTW: 2,
    battery: { cells: 12, capacity: 16, usable: .8, mass: 3.9, resistance: .02 },
    masses: [
      { id: "frame", name: "Carbon frame and arms", mass: 1.45 }, { id: "motors", name: "Motors", mass: .92 }, { id: "props", name: "Propellers", mass: .2 },
      { id: "gimbal", name: "Gimbal and camera body", mass: .55 }, { id: "avionics", name: "Flight controller, GNSS, radio", mass: .18 }, { id: "wiring", name: "Wiring and connectors", mass: .3 },
    ],
    avionicsPower: 14, temperature: 15, pressure: 101325, altitude: 60, requiredTime: 18, minTW: 1.8, capture: .15, arm: 368, propulsion: null,
  };
}
const finite = (n, lo, hi) => Number.isFinite(n) && n >= lo && n <= hi;
export function validateFlight(d) {
  if (!d || typeof d.name !== "string" || d.name.length > 60 || !Number.isInteger(d.rotors) || ![3, 4, 6, 8].includes(d.rotors) || !finite(d.prop, 50, 2000)
    || !finite(d.fm, .3, .85) || !finite(d.motorEfficiency, .3, .98) || !finite(d.escEfficiency, .5, .995) || !finite(d.designTW, 1.1, 5) || !d.battery
    || !Number.isInteger(d.battery.cells) || d.battery.cells < 1 || d.battery.cells > 24 || !finite(d.battery.capacity, .1, 200) || !finite(d.battery.usable, .3, 1)
    || !finite(d.battery.mass, .01, 100) || !finite(d.battery.resistance, 0, 1) || !Array.isArray(d.masses) || d.masses.length > 24 || new Set(d.masses.map(m => m?.id)).size !== d.masses.length
    || !d.masses.every(m => m && typeof m.id === "string" && /^[a-zA-Z][a-zA-Z0-9_-]{0,23}$/.test(m.id) && typeof m.name === "string" && m.name.trim() && m.name.length <= 60 && finite(m.mass, 0, 500))
    || !finite(d.avionicsPower, 0, 5000) || !finite(d.temperature, -40, 55) || !finite(d.pressure, 5e4, 1.1e5) || !finite(d.altitude, 1, 5000)
    || !finite(d.requiredTime, 0, 600) || !finite(d.minTW, 1, 5) || !finite(d.capture, 0, 1) || !finite(d.arm ?? 368, 20, 5000)) throw new Error("Invalid vehicle settings.");
  const text = (v, n) => v === undefined || typeof v === "string" && v.length <= n;
  if (!d.masses.every(m => text(m.part, 80) && (m.source === undefined || typeof m.source === "string" && /^https:\/\/\S{4,300}$/.test(m.source)))) throw new Error("A mass item's part is text and its source an https link.");
  const P = d.propulsion ? validatePropulsion(d.propulsion) : null;
  return { name: d.name, rotors: d.rotors, prop: d.prop, fm: d.fm, motorEfficiency: d.motorEfficiency, escEfficiency: d.escEfficiency, designTW: d.designTW,
    battery: { cells: d.battery.cells, capacity: d.battery.capacity, usable: d.battery.usable, mass: d.battery.mass, resistance: d.battery.resistance },
    masses: d.masses.map(m => ({ id: m.id, name: m.name, mass: m.mass, ...(m.part && { part: m.part }), ...(m.source && { source: m.source }) })), avionicsPower: d.avionicsPower, temperature: d.temperature, pressure: d.pressure,
    altitude: d.altitude, requiredTime: d.requiredTime, minTW: d.minTW, capture: d.capture, arm: d.arm ?? 368, propulsion: P };
}
// A measured table: rows of throttle, volts, amps, watts, rpm, torque (N·m) and
// grams of thrust, rising in thrust, with the motor's own constants.
const ROW = ["throttle", "volts", "amps", "watts", "rpm", "torque", "grams"];
export function validatePropulsion(p) {
  const f = (v, lo, hi, what) => { if (!finite(v, lo, hi)) throw new Error(`Propulsion: ${what}.`); return v; };
  if (!p || typeof p.motor !== "string" || !p.motor.trim() || p.motor.length > 80 || typeof p.propeller !== "string" || !p.propeller.trim() || p.propeller.length > 80) throw new Error("Propulsion: name the motor and the propeller.");
  if (p.source !== undefined && (typeof p.source !== "string" || !/^https:\/\/\S{4,300}$/.test(p.source))) throw new Error("Propulsion: the source is an https link.");
  if (!Array.isArray(p.table) || p.table.length < 2 || p.table.length > 64) throw new Error("Propulsion: a table has 2 to 64 rows.");
  const table = p.table.map(r => Object.fromEntries(ROW.map(k => [k, f(r?.[k], 0, 1e6, `every row needs ${ROW.join(", ")}`)])));
  if (!table.every((r, i) => !i || r.grams > table[i - 1].grams)) throw new Error("Propulsion: table rows rise in thrust.");
  return { motor: p.motor, propeller: p.propeller, ...(p.source && { source: p.source }), kv: f(p.kv, 1, 1e5, "KV"), resistance: f(p.resistance, 0, 100, "motor resistance in ohms"), peakCurrent: f(p.peakCurrent, 0, 1e4, "peak current"), maxPower: f(p.maxPower, 0, 1e6, "maximum power"),
    rotorInertia: f(p.rotorInertia, 0, 10, "rotor inertia in kg·m²"), density: f(p.density ?? 1.2, .5, 1.5, "test air density"), blockage: f(p.blockage ?? 0, 0, .5, "thrust lost to the arm under the disc"), table };
}
// Bench data at thrust T (N) in air of density ρ. At one rpm, thrust, torque and
// power all scale with ρ (similarity), so look up the equivalent bench thrust.
// Below the table, T and Q go as n² and P as n³; above it the rotor cannot do it.
export function propulsionAt(p, T, rho) {
  const k = rho / p.density, g = T / (1 - p.blockage) / G * 1000 / k, t = p.table, top = t.at(-1);
  if (g > top.grams) return null;
  if (g < t[0].grams) { const r = Math.sqrt(g / t[0].grams), a = t[0]; return { throttle: a.throttle * r, rpm: a.rpm * r, torque: a.torque * r * r * k, watts: a.watts * r ** 3 * k }; }
  const i = Math.min(t.length - 2, t.findIndex((r, n) => g <= t[n + 1]?.grams)), a = t[i], b = t[i + 1], f = (g - a.grams) / (b.grams - a.grams), mid = key => a[key] + (b[key] - a[key]) * f;
  return { throttle: mid("throttle"), rpm: mid("rpm"), torque: mid("torque") * k, watts: mid("watts") * k };
}
// The most one rotor gives at pack voltage V: speed follows voltage, thrust n².
export const maxThrust = (p, V, rho) => { const top = p.table.at(-1); return top.grams / 1000 * G * Math.min(1, V / top.volts) ** 2 * (rho / p.density) * (1 - p.blockage); };
// Actuator-disc rotor at thrust T (N): induced velocity, ideal and shaft power.
export function rotor(T, rho, area, fm) {
  const vi = Math.sqrt(T / (2 * rho * area)), ideal = T * vi;
  return { vi, ideal, shaft: ideal / fm, wake: 2 * vi };
}
export function flightResults(d, ov = {}) {
  const t = ov.air_temperature ?? d.temperature, rho = d.pressure / (287.05 * (t + 273.15)), p = d.propulsion, Kt = p ? 60 / (2 * Math.PI * p.kv) : NaN;
  const linked = (ov.esc_mass ?? 0) * d.rotors + (ov.camera_mass ?? 0) + (ov.parts_mass ?? 0);
  const fixed = d.masses.reduce((s, m) => s + m.mass, 0), mass = fixed + d.battery.mass + linked, weight = mass * G;
  const area = Math.PI * (d.prop / 2000) ** 2, V = d.battery.cells * cellVoltage(.5), R = d.battery.resistance, full = d.battery.cells * cellVoltage(1);
  const point = (tw, Vbus = V) => {
    const T = tw * weight / d.rotors, r = rotor(T, rho, area, d.fm), bench = p ? propulsionAt(p, T, rho) : null;
    const shaft = bench ? bench.torque * bench.rpm * Math.PI / 30 : r.shaft, input = bench ? bench.watts : r.shaft / d.motorEfficiency / d.escEfficiency, motor = bench ? shaft : r.shaft / d.motorEfficiency;
    const P = input * d.rotors + d.avionicsPower, I = sag(Vbus, R, P), Vt = Vbus - I * R, phase = bench ? bench.torque / Kt : NaN;
    // Motor back-EMF and IR drop against the bus: the ESC's average duty.
    const duty = bench ? Math.min(1, (bench.rpm / p.kv + phase * p.resistance) / Vt) : NaN;
    return { T, ...r, shaft, motor, input, P, I, Vt, esc: input / Vt, escLoss: input - motor, feasible: !p || !!bench, throttle: bench?.throttle ?? NaN, rpm: bench?.rpm ?? NaN, torque: bench?.torque ?? NaN, phase, duty, fm: r.ideal / shaft };
  };
  const hover = point(1), design = point(d.designTW);
  // Discharge at hover power: the pack's open-circuit voltage falls with charge
  // and the current rises to hold the power.
  let seconds = 0, soc = 1, Vend = NaN; const steps = 80, dQ = d.battery.capacity * 3600 * d.battery.usable / steps;
  for (let i = 0; i < steps; i++) { const Voc = d.battery.cells * cellVoltage(soc - dQ / (2 * d.battery.capacity * 3600)), I = sag(Voc, R, hover.P); seconds += dQ / I; Vend = Voc - I * R; soc -= dQ / (d.battery.capacity * 3600); }
  const time = seconds / 60, thrustFull = p ? maxThrust(p, full - hover.I * R, rho) : NaN, thrustEnd = p ? maxThrust(p, Vend, rho) : NaN;
  const twFull = thrustFull * d.rotors / weight, twEnd = thrustEnd * d.rotors / weight, clearance = 2 * d.arm * Math.sin(Math.PI / d.rotors) - d.prop;
  // A rigid-body estimate of the inertia. What sits at the rotors is the motors,
  // propellers and ESCs; parts are spread along the arms; the rest is a body box.
  const R2 = (d.arm / 1000) ** 2, atRotors = (ov.esc_mass ?? 0) * d.rotors + d.masses.filter(m => /motor|prop/i.test(`${m.id} ${m.name}`)).reduce((s, m) => s + m.mass, 0), arms = ov.parts_mass ?? 0, body = Math.max(0, mass - atRotors - arms);
  const Izz = atRotors * R2 + arms * R2 / 3 + body * (.12 ** 2 + .09 ** 2) / 12 * 2, Ixx = atRotors * R2 / 2 + arms * R2 / 6 + body * (.12 ** 2 + .08 ** 2) / 12;
  return { rho, mass, weight, fixed, linked, area, V, Vend, hover, design, time, diskLoading: weight / (d.rotors * area), wakePressure: .5 * rho * hover.wake ** 2 * d.capture, thrustFull, thrustEnd, twFull, twEnd, clearance, Ixx, Izz, atRotors };
}
export const flightAssumptions = d => [
  ...(d.propulsion ? [
    { id: "flight.bench", label: `${d.propulsion.motor} with ${d.propulsion.propeller}: measured bench table`, kind: "nominal", source: "Manufacturer's static bench data, marked for reference only. Confirm on your own thrust stand; the Testing studio correlates a measured sweep." },
    { id: "flight.density", label: `Bench air density taken as ${d.propulsion.density} kg/m³`, kind: "nominal", source: "The table does not state its conditions; a 5% error here moves hover power by about 5%." },
    { id: "flight.blockage", label: `${Math.round(d.propulsion.blockage * 100)}% of thrust lost to the arm under the disc`, kind: "nominal", source: "Static bench tables have no airframe below the rotor; measure on the assembled arm." },
    { id: "flight.inertia", label: `Rotor inertia ${d.propulsion.rotorInertia} kg·m²`, kind: "nominal", source: "Estimated from blade mass and diameter; a spin-down test measures it." },
  ] : [
    { id: "flight.fm", label: `Rotor figure of merit ${d.fm}`, kind: "nominal", source: "Small rotors measure 0.5–0.75; use a thrust-stand test of the chosen propeller." },
    { id: "flight.efficiency", label: `Motor ${d.motorEfficiency} and ESC ${d.escEfficiency} efficiency`, kind: "nominal", source: "Constant over the operating range; real maps vary with load and speed." },
  ]),
  { id: "flight.battery", label: `${d.battery.cells}S ${d.battery.capacity} Ah pack, ${Math.round(d.battery.usable * 100)}% usable, ${d.battery.resistance * 1000} mΩ`, kind: "nominal", source: "Typical LiPo rest-voltage curve; pack resistance is an estimate unless measured under load." },
  { id: "flight.capture", label: `ESC duct captures ${Math.round(d.capture * 100)}% of the wake dynamic pressure`, kind: "nominal", source: "Placement under the rotor; measure with a pitot in the wake." },
  { id: "flight.hover", label: "Hover in still air, no ground effect, no forward flight", kind: "nominal", source: "Momentum theory of an isolated actuator disc." },
];
