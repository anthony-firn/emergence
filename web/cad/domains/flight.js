// Multirotor sizing: mass budget, hover and design-point thrust, rotor power
// by momentum theory with a figure of merit, battery current and endurance.
// Ideal induced power P = T^1.5 / √(2ρA) is exact for an actuator disc; the
// figure of merit, motor and ESC efficiencies are stated assumptions. The
// design thrust-to-weight point sets the continuous ESC current, which is the
// MOSFET RMS current of the mechanical/electrothermal power stage.
export const G = 9.80665;
export function flightExample() {
  return {
    name: "Quadcopter", rotors: 4, prop: 381, fm: .65, motorEfficiency: .82, escEfficiency: .97, designTW: 2,
    battery: { cells: 12, capacity: 16, usable: .8, mass: 3.9, resistance: .02 },
    masses: [
      { id: "frame", name: "Carbon frame and arms", mass: 1.45 }, { id: "motors", name: "Motors", mass: .92 }, { id: "props", name: "Propellers", mass: .2 },
      { id: "gimbal", name: "Gimbal and camera body", mass: .55 }, { id: "avionics", name: "Flight controller, GNSS, radio", mass: .18 }, { id: "wiring", name: "Wiring and connectors", mass: .3 },
    ],
    avionicsPower: 14, temperature: 15, pressure: 101325, altitude: 60, requiredTime: 18, minTW: 1.8, capture: .15,
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
    || !finite(d.requiredTime, 0, 600) || !finite(d.minTW, 1, 5) || !finite(d.capture, 0, 1)) throw new Error("Invalid vehicle settings.");
  return { name: d.name, rotors: d.rotors, prop: d.prop, fm: d.fm, motorEfficiency: d.motorEfficiency, escEfficiency: d.escEfficiency, designTW: d.designTW,
    battery: { cells: d.battery.cells, capacity: d.battery.capacity, usable: d.battery.usable, mass: d.battery.mass, resistance: d.battery.resistance },
    masses: d.masses.map(m => ({ id: m.id, name: m.name, mass: m.mass })), avionicsPower: d.avionicsPower, temperature: d.temperature, pressure: d.pressure,
    altitude: d.altitude, requiredTime: d.requiredTime, minTW: d.minTW, capture: d.capture };
}
// Actuator-disc rotor at thrust T (N): induced velocity, ideal and shaft power.
export function rotor(T, rho, area, fm) {
  const vi = Math.sqrt(T / (2 * rho * area)), ideal = T * vi;
  return { vi, ideal, shaft: ideal / fm, wake: 2 * vi };
}
export function flightResults(d, ov = {}) {
  const t = ov.air_temperature ?? d.temperature, rho = d.pressure / (287.05 * (t + 273.15));
  const linked = (ov.esc_mass ?? 0) * d.rotors + (ov.camera_mass ?? 0) + (ov.parts_mass ?? 0);
  const fixed = d.masses.reduce((s, m) => s + m.mass, 0), mass = fixed + d.battery.mass + linked, weight = mass * G;
  const area = Math.PI * (d.prop / 2000) ** 2, V = d.battery.cells * 3.7;
  const point = tw => {
    const T = tw * weight / d.rotors, r = rotor(T, rho, area, d.fm), motor = r.shaft / d.motorEfficiency, input = motor / d.escEfficiency;
    // Pack sag at this draw: V_t = V − I·R, solved I = (V − √(V² − 4RP)) / 2R.
    const P = input * d.rotors + d.avionicsPower, disc = V * V - 4 * d.battery.resistance * P;
    const I = d.battery.resistance > 0 ? (disc >= 0 ? (V - Math.sqrt(disc)) / (2 * d.battery.resistance) : NaN) : P / V;
    const Vt = V - I * d.battery.resistance;
    return { T, ...r, motor, input, P, I, Vt, esc: input / Vt, escLoss: input - motor };
  };
  const hover = point(1), design = point(d.designTW);
  const energy = d.battery.capacity * V * d.battery.usable, time = 60 * energy / hover.P;
  return { rho, mass, weight, fixed, linked, area, V, hover, design, time, diskLoading: weight / (d.rotors * area),
    wakePressure: .5 * rho * hover.wake ** 2 * d.capture, maxTW: NaN };
}
export const flightAssumptions = d => [
  { id: "flight.fm", label: `Rotor figure of merit ${d.fm}`, kind: "nominal", source: "Small rotors measure 0.5–0.75; use a thrust-stand test of the chosen propeller." },
  { id: "flight.efficiency", label: `Motor ${d.motorEfficiency} and ESC ${d.escEfficiency} efficiency`, kind: "nominal", source: "Constant over the operating range; real maps vary with load and speed." },
  { id: "flight.battery", label: `${d.battery.cells}S ${d.battery.capacity} Ah pack at 3.7 V/cell, ${Math.round(d.battery.usable * 100)}% usable`, kind: "nominal", source: "Nominal cell voltage; sag uses the stated pack resistance." },
  { id: "flight.capture", label: `ESC duct captures ${Math.round(d.capture * 100)}% of the wake dynamic pressure`, kind: "nominal", source: "Placement under the rotor; measure with a pitot in the wake." },
  { id: "flight.hover", label: "Hover in still air, no ground effect, no forward flight", kind: "nominal", source: "Momentum theory of an isolated actuator disc." },
];
