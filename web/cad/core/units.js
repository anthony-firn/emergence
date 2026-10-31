// Dimensioned quantities for cross-domain links. A unit is a scale and offset
// onto SI base dimensions [m, kg, s, A, K]; links refuse incompatible pairs.
const D = (m = 0, kg = 0, s = 0, A = 0, K = 0) => [m, kg, s, A, K];
const U = (dim, scale = 1, offset = 0, label) => ({ dim, scale, offset, label });
const NONE = D(), LENGTH = D(1), AREA = D(2), VOLUME = D(3), TIME = D(0, 0, 1), POWER = D(2, 1, -3), TEMP = D(0, 0, 0, 0, 1);
const PRESSURE = D(-1, 1, -2), VELOCITY = D(1, 0, -1), FLOW = D(3, 0, -1), CURRENT = D(0, 0, 0, 1), VOLT = D(2, 1, -3, -1), OHM = D(2, 1, -3, -2);
export const UNITS = {
  "": U(NONE, 1, 0, ""), "%": U(NONE, .01, 0, "%"), count: U(NONE, 1, 0, ""),
  kg: U(D(0, 1)), g: U(D(0, 1), 1e-3), N: U(D(1, 1, -2)), min: U(TIME, 60),
  m: U(LENGTH), cm: U(LENGTH, 1e-2), mm: U(LENGTH, 1e-3), "µm": U(LENGTH, 1e-6), um: U(LENGTH, 1e-6, 0, "µm"), nm: U(LENGTH, 1e-9),
  m2: U(AREA, 1, 0, "m²"), mm2: U(AREA, 1e-6, 0, "mm²"), m3: U(VOLUME, 1, 0, "m³"), mm3: U(VOLUME, 1e-9, 0, "mm³"), L: U(VOLUME, 1e-3),
  s: U(TIME), ms: U(TIME, 1e-3), h: U(TIME, 3600),
  W: U(POWER), mW: U(POWER, 1e-3), kW: U(POWER, 1e3),
  K: U(TEMP), degC: U(TEMP, 1, 273.15, "°C"), "°C": U(TEMP, 1, 273.15, "°C"),
  Pa: U(PRESSURE), kPa: U(PRESSURE, 1e3), bar: U(PRESSURE, 1e5),
  "m/s": U(VELOCITY), "m3/s": U(FLOW, 1, 0, "m³/s"), "L/s": U(FLOW, 1e-3), "L/min": U(FLOW, 1e-3 / 60), "m3/h": U(FLOW, 1 / 3600, 0, "m³/h"),
  A: U(CURRENT), mA: U(CURRENT, 1e-3), V: U(VOLT), mV: U(VOLT, 1e-3), ohm: U(OHM, 1, 0, "Ω"), mohm: U(OHM, 1e-3, 0, "mΩ"),
  Hz: U(D(0, 0, -1)), rpm: U(D(0, 0, -1), 1 / 60), MPa: U(PRESSURE, 1e6), GPa: U(PRESSURE, 1e9), kgm2: U(D(2, 1), 1, 0, "kg·m²"), Nm: U(D(2, 1, -2), 1, 0, "N·m"),
  "K/W": U(D(-2, -1, 3, 0, 1)), "W/m2": U(D(0, 1, -3), 1, 0, "W/m²"), "W/mm2": U(D(0, 1, -3), 1e6, 0, "W/mm²"), "W/K": U(D(2, 1, -3, 0, -1)),
  "W/m2K": U(D(0, 1, -3, 0, -1), 1, 0, "W/m²K"), "1/h": U(D(0, 0, -1), 1 / 3600, 0, "1/h"),
};
export function unit(name) {
  const u = UNITS[name ?? ""];
  if (!u) throw new Error(`Unknown unit "${name}".`);
  return u;
}
export const compatible = (a, b) => unit(a).dim.every((n, i) => n === unit(b).dim[i]);
export const toSI = (value, name) => { const u = unit(name); return value * u.scale + u.offset; };
export const fromSI = (value, name) => { const u = unit(name); return (value - u.offset) / u.scale; };
export function convert(value, from, to) {
  if (!compatible(from, to)) throw new Error(`Cannot convert ${label(from) || "a dimensionless value"} to ${label(to) || "a dimensionless value"}.`);
  return fromSI(toSI(value, from), to);
}
// Temperature differences convert by scale only.
export const convertDelta = (value, from, to) => { if (!compatible(from, to)) throw new Error("Incompatible units."); return value * unit(from).scale / unit(to).scale; };
// Display only: unknown units (ps, pJ…) print as written; links still refuse them.
export const label = name => UNITS[name ?? ""] ? UNITS[name ?? ""].label ?? name : name;
export function format(value, name, digits = 3) {
  if (value === undefined || value === null || !Number.isFinite(value)) return "—";
  const a = Math.abs(value), text = a !== 0 && (a >= 1e5 || a < 1e-3) ? value.toExponential(2) : Number(value.toPrecision(digits)).toLocaleString("en", { maximumFractionDigits: 6 });
  const l = label(name);
  return l ? `${text} ${l}` : text;
}
