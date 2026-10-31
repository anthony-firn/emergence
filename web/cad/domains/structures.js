// A motor arm as a cantilevered tube carrying thrust and torque at its tip.
// Euler–Bernoulli bending and Saint-Venant torsion in closed form, combined by
// von Mises; the first bending mode by Rayleigh's method with the standard 33/140
// of the arm's own mass at the tip. Resonance is a Campbell check: the mode must
// clear once-per-rev and blade-pass across the flying speed range. A clamped root
// is an upper bound on stiffness: a real joint is softer, so measure the mode with
// a tap test. Laminates are not isotropic; the allowable is the supplier's figure.
export const structuresExample = () => ({ name: "Motor arm", od: 25, wall: 1.5, length: 310, modulus: 90, density: 1.55, strength: 600, safetyFactor: 1.5, deflectionLimit: 2, separation: 15, blades: 2, tipMass: .2, loadThrust: 40, loadTorque: 1, source: "" });
const finite = (v, lo, hi) => Number.isFinite(v) && v >= lo && v <= hi;
export const STRUCTURE_FIELDS = { od: [1, 1000], wall: [.05, 100], length: [5, 5000], modulus: [.1, 1000], density: [.1, 20], strength: [1, 10000], safetyFactor: [1, 10], deflectionLimit: [.01, 50], separation: [0, 100], blades: [1, 12], tipMass: [0, 100], loadThrust: [0, 1e5], loadTorque: [0, 1e5] };
export function validateStructures(d) {
  if (!d || typeof d.name !== "string" || !d.name.trim() || d.name.length > 60) throw new Error("Name the structure.");
  for (const [k, [lo, hi]] of Object.entries(STRUCTURE_FIELDS)) if (!finite(d[k], lo, hi)) throw new Error(`${k} must be ${lo}–${hi}.`);
  if (2 * d.wall >= d.od) throw new Error("The wall is thicker than the tube allows.");
  if (d.source && (typeof d.source !== "string" || !/^https:\/\/\S{4,300}$/.test(d.source))) throw new Error("A material source is an https link.");
  return { name: d.name, ...Object.fromEntries(Object.keys(STRUCTURE_FIELDS).map(k => [k, d[k]])), source: d.source || "" };
}
// Stresses in MPa, lengths in mm, loads in N and N·m.
export function armResults(d, ov = {}) {
  const od = ov.od ?? d.od, wall = ov.wall ?? d.wall, L = (ov.length ?? d.length) / 1000, T = ov.thrust ?? d.loadThrust, Q = ov.torque ?? d.loadTorque, tip = ov.tip_mass ?? d.tipMass;
  const ro = od / 2000, ri = ro - wall / 1000, area = Math.PI * (ro * ro - ri * ri), I = Math.PI / 4 * (ro ** 4 - ri ** 4), J = 2 * I, E = d.modulus * 1e9;
  const M = T * L * d.safetyFactor, sigma = M * ro / I / 1e6, tau = Q * d.safetyFactor * ro / J / 1e6, vonMises = Math.sqrt(sigma * sigma + 3 * tau * tau), margin = d.strength / vonMises - 1;
  const deflection = T * L ** 3 / (3 * E * I) * 1000, armMass = area * L * d.density * 1000, k = 3 * E * I / L ** 3, modal = tip + 33 / 140 * armMass, frequency = Math.sqrt(k / modal) / (2 * Math.PI);
  const lo = (ov.rpm_low ?? 0) / 60, hi = (ov.rpm_high ?? 0) / 60, sep = d.separation / 100;
  const bands = hi > 0 ? [{ name: "once per rev", lo, hi }, { name: "blade pass", lo: lo * d.blades, hi: hi * d.blades }] : [];
  const conflicts = bands.filter(b => frequency >= b.lo * (1 - sep) && frequency <= b.hi * (1 + sep));
  return { area, I, J, sigma, tau, vonMises, margin, deflection, armMass, stiffness: k, modal, frequency, bands, conflicts, T, Q, tip, od, wall, length: L * 1000 };
}
export const structuresAssumptions = d => [
  { id: "structures.material", label: `Tube modulus ${d.modulus} GPa, allowable ${d.strength} MPa`, kind: "nominal", source: "Supplier figures for one layup, for reference only. Tube suppliers ask for physical testing." },
  { id: "structures.root", label: "Arm root taken as fully clamped", kind: "nominal", source: "An upper bound on stiffness and frequency; a tap test on the real joint measures the true mode." },
  { id: "structures.factor", label: `Safety factor ${d.safetyFactor} on peak thrust and torque`, kind: "nominal", source: "Covers manoeuvres but not a crash, and there is no fatigue or impact analysis." },
];
