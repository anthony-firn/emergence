// One three-phase inverter from its MOSFET's datasheet, in block (120°)
// commutation with synchronous rectification: two FETs carry the phase current
// at any moment, and one leg hard-switches each PWM cycle. Conduction uses
// RDS(on) at the junction temperature; switching uses the gate charges and the
// gate current they imply, the standard estimate. It is an analytic estimate,
// not a device simulation: ringing, layout inductance and reverse recovery
// beyond the datasheet figures need a double-pulse test on the real board.
export const POWER_SOURCE = "https://www.infineon.com/assets/row/public/documents/24/49/infineon-bsc040n10ns5-datasheet-en.pdf";
const finite = (v, lo, hi) => Number.isFinite(v) && v >= lo && v <= hi;
export const CAPACITOR_SOURCE = "https://www.digikey.com.mx/en/products/detail/rubycon/63ZLJ470M12-5X20/15667867";
export const CAPACITOR_DERATE = .8; // of the ripple rating, for life at temperature
export const CAPACITOR_FIELDS = { count: [1, 12], capacitance: [.1, 1e5], esr: [.1, 1e4], rating: [.01, 1000], voltage: [1, 5000], grams: [0, 1000], limit: [.01, 100] };
// Listed at 470 µF, 63 V, 33 mΩ and 2.4 A at 100 kHz; each weighs about 4 g (an estimate).
const capacitorExample = () => ({ part: "Rubycon 63ZLJ470M12.5X20", count: 4, capacitance: 470, esr: 33, rating: 2.4, voltage: 63, grams: 4, limit: 2, source: CAPACITOR_SOURCE });
export const MOSFET_FIELDS = { vds: [10, 1000], rds25: [.1, 1000], rdsTc: [0, .02], qg: [1, 2000], qsw: [1, 1000], qoss: [0, 2000], qrr: [0, 5000], vplateau: [1, 15], rg: [0, 50], vsd: [.1, 3], rthJc: [.01, 50], tjMax: [100, 200] };
export const powerExample = () => ({
  name: "ESC power stage", pwm: 24000, deadtime: 100, gateDrive: 12, gateResistor: 3, driverPeak: 3.7, ringing: 1.5, efficiencyLimit: 95, phase: 15, bus: 44, duty: 60, junction: 60, peakBus: 50.4,
  // Infineon BSC040N10NS5, final data sheet Rev. 2.4, 2020-05-15. rdsTc is read
  // from diagram 9 (about +0.6 %/K of the 25 °C value at the typical curve).
  mosfet: { part: "Infineon BSC040N10NS5", vds: 100, rds25: 3.4, rdsTc: .006, qg: 58, qsw: 19, qoss: 75, qrr: 90, vplateau: 4.6, rg: 1.3, vsd: .86, rthJc: .5, tjMax: 150, source: POWER_SOURCE },
  capacitor: capacitorExample(),
});
export function validatePower(d) {
  if (!d || typeof d.name !== "string" || !d.name.trim() || d.name.length > 60 || !finite(d.pwm, 1e3, 1e6) || !finite(d.deadtime, 0, 1e4) || !finite(d.gateDrive, 4, 20) || !finite(d.gateResistor, 0, 100) || !finite(d.driverPeak, .1, 100) || !finite(d.ringing, 1, 3) || !finite(d.efficiencyLimit, 0, 100) || !finite(d.phase, 0, 1e4) || !finite(d.bus, 1, 2000) || !finite(d.duty, 0, 100) || !finite(d.junction, -40, 200) || !finite(d.peakBus, 1, 2000) || !d.mosfet) throw new Error("Invalid power stage settings.");
  const m = d.mosfet; if (typeof m.part !== "string" || !m.part.trim() || m.part.length > 80 || m.source !== undefined && (typeof m.source !== "string" || !/^https:\/\/\S{4,300}$/.test(m.source))) throw new Error("Name the MOSFET; its source is an https link.");
  for (const [k, [lo, hi]] of Object.entries(MOSFET_FIELDS)) if (!finite(m[k], lo, hi)) throw new Error(`MOSFET ${k} must be ${lo}–${hi}.`);
  const c = d.capacitor ?? capacitorExample(); if (typeof c.part !== "string" || !c.part.trim() || c.part.length > 80 || c.source !== undefined && (typeof c.source !== "string" || !/^https:\/\/\S{4,300}$/.test(c.source))) throw new Error("Name the capacitor; its source is an https link.");
  for (const [k, [lo, hi]] of Object.entries(CAPACITOR_FIELDS)) if (!finite(c[k], lo, hi)) throw new Error(`Capacitor ${k} must be ${lo}–${hi}.`);
  if (!Number.isInteger(c.count)) throw new Error("A capacitor bank has a whole number of parts.");
  return { name: d.name, pwm: d.pwm, deadtime: d.deadtime, gateDrive: d.gateDrive, gateResistor: d.gateResistor, driverPeak: d.driverPeak, ringing: d.ringing, efficiencyLimit: d.efficiencyLimit, phase: d.phase, bus: d.bus, duty: d.duty, junction: d.junction, peakBus: d.peakBus,
    mosfet: { part: m.part, ...Object.fromEntries(Object.keys(MOSFET_FIELDS).map(k => [k, m[k]])), ...(m.source && { source: m.source }) },
    capacitor: { part: c.part, ...Object.fromEntries(Object.keys(CAPACITOR_FIELDS).map(k => [k, c[k]])), ...(c.source && { source: c.source }) } };
}
// phase: motor current (A); bus: DC volts; duty: 0–1; tj: junction °C.
export function stageLoss(d, { phase, bus, duty, tj }) {
  const m = d.mosfet, f = d.pwm, I = phase, D = Math.min(1, Math.max(0, duty)), R = m.rds25 * 1e-3 * (1 + m.rdsTc * (tj - 25));
  // The high side is on for duty D of its third of a cycle; the low side for the
  // rest of that third (rectifying) and all of its own third.
  const hiRms = I * Math.sqrt(D / 3), loRms = I * Math.sqrt((2 - D) / 3), conduction = 2 * I * I * R;
  const Rtotal = d.gateResistor + m.rg, Ion = Math.max(1e-6, (d.gateDrive - m.vplateau) / Rtotal), Ioff = m.vplateau / Rtotal, tOn = m.qsw * 1e-9 / Ion, tOff = m.qsw * 1e-9 / Ioff;
  const Esw = .5 * bus * I * (tOn + tOff), Eoss = m.qoss * 1e-9 * bus, Err = m.qrr * 1e-9 * bus * Math.min(1, I / 50), switching = (Esw + Eoss + Err) * f;
  const gate = 2 * m.qg * 1e-9 * d.gateDrive * f, dead = 2 * m.vsd * I * d.deadtime * 1e-9 * f, total = conduction + switching + gate + dead;
  // Per device: commutation shares the switching among the three high sides and
  // the dead-time diode loss among the three low sides.
  const high = hiRms * hiRms * R + switching / 3, low = loRms * loRms * R + dead / 3;
  const hottest = high >= low ? { rms: hiRms, extra: switching / 3, device: "high side" } : { rms: loRms, extra: dead / 3, device: "low side" };
  return { R, hiRms, loRms, conduction, switching, gate, dead, total, high, low, hottest, ripple: I * Math.sqrt(D * (1 - D)), gatePeak: d.gateDrive / Rtotal, tOn, tOff, other: Math.max(0, total - gate - hottest.rms ** 2 * R - hottest.extra) };
}
// The bus capacitors carry the inverter's input ripple, I√(D(1−D)) RMS, shared
// across the bank. Ripple voltage is the charge drawn each cycle plus the step
// across the bank's ESR.
export function busRipple(d, { phase, duty }) {
  const c = d.capacitor, D = Math.min(1, Math.max(0, duty)), rms = phase * Math.sqrt(D * (1 - D)), esr = c.esr * 1e-3 / c.count;
  return { rms, share: rms / c.count, volts: phase * D * (1 - D) / (d.pwm * c.count * c.capacitance * 1e-6) + phase * esr, loss: rms * rms * esr, mass: c.count * c.grams / 1000 };
}
export const powerAssumptions = d => [
  { id: "power.mosfet", label: `${d.mosfet.part}: datasheet typical values`, kind: "nominal", source: "RDS(on) rises with temperature per the datasheet's curve; use the maximum column for a worst case." },
  { id: "power.switching", label: "Hard-switching loss from gate charge and gate current", kind: "nominal", source: "A standard analytic estimate. A double-pulse test on the real board measures it." },
  { id: "power.capacitor", label: `${d.capacitor.count} × ${d.capacitor.part} on the battery leads`, kind: "nominal", source: "Listed values: the listing does not state the frequency of the ESR, the ripple rating is at 100 kHz, and the mass is an estimate. Lead inductance is not modelled." },
  { id: "power.block", label: "Block (120°) commutation with synchronous rectification", kind: "nominal", source: "Field-oriented control spreads the current differently; losses are of the same order." },
];
