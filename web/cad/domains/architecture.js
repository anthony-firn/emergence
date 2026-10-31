// Rooms on a plan, metres. Each room's steady heat balance is exact for its
// inputs: envelope conductance UA, solar gain through glazing, internal gains
// and supply air. A supplied room floats to
//   T = (C·T_supply + UA·T_out + Q_internal + Q_solar) / (C + UA),  C = ρ·cp·V̇,
// a conditioned room reports the cooling load at its set point. Interior
// walls carry no heat (adjacent rooms assumed at similar temperatures).
export const SIDES = ["north", "east", "south", "west"];
// Peak vertical irradiance relative to the site value by orientation,
// a nominal mid-latitude summer afternoon; replace with site data.
export const ORIENTATION = { north: .25, east: .8, south: .75, west: 1 };
export const SENSIBLE_PER_PERSON = 75; // W, seated light work (ASHRAE Fundamentals ch. 18 table range).
const wall = (exterior, U = .35, window = 0) => ({ exterior, U, window });
// The starter building is one room with the geometry and envelope of
// ANSI/ASHRAE Standard 140 BESTEST case 600: 8 × 6 × 2.7 m, 12 m² of south
// glazing, light walls U 0.514, roof U 0.318 W/m²K (as published; confirm
// against the standard before relying on them). Here it is ventilated and
// solved in steady state, not the standard's annual dynamic test.
export const BESTEST = "ANSI/ASHRAE Standard 140, BESTEST case 600 geometry and envelope";
export function architectureExample() {
  return {
    site: { name: "Site", outdoor: 32, solar: 450 },
    glazing: { U: 3, shgc: .787 },
    units: 6, limit: 35,
    rooms: [
      { id: "room", name: "Equipment room (BESTEST 600)", x: 0, y: 0, width: 8, depth: 6, height: 2.7, people: 1, lighting: 8, mode: "supplied", supply: 1400, supplyTemp: 20, setpoint: 24, equipment: true, roofU: .318,
        walls: { north: wall(true, .514, 0), east: wall(true, .514, 0), south: wall(true, .514, 12), west: wall(true, .514, 0) } },
    ],
  };
}
const finite = (n, lo, hi) => Number.isFinite(n) && n >= lo && n <= hi;
export function validateArchitecture(d) {
  if (!d || !d.site || typeof d.site.name !== "string" || d.site.name.length > 60 || !finite(d.site.outdoor, -40, 55) || !finite(d.site.solar, 0, 1200)
    || !d.glazing || !finite(d.glazing.U, .3, 7) || !finite(d.glazing.shgc, .05, .9) || !Number.isInteger(d.units) || d.units < 0 || d.units > 10000 || !finite(d.limit, 5, 80)
    || !finite(d.equipmentPower ?? 0, 0, 1e6) || !Array.isArray(d.rooms) || d.rooms.length > 32 || new Set(d.rooms.map(r => r?.id)).size !== d.rooms.length) throw new Error("Invalid building settings.");
  const rooms = d.rooms.map(r => {
    if (!r || typeof r.id !== "string" || !/^[a-zA-Z][a-zA-Z0-9_-]{0,23}$/.test(r.id) || typeof r.name !== "string" || !r.name.trim() || r.name.length > 60
      || !finite(r.x, -500, 500) || !finite(r.y, -500, 500) || !finite(r.width, .5, 200) || !finite(r.depth, .5, 200) || !finite(r.height, 2, 20)
      || !finite(r.people, 0, 2000) || !finite(r.lighting, 0, 50) || !["supplied", "conditioned"].includes(r.mode) || !finite(r.supply, 0, 1e6)
      || !finite(r.supplyTemp, 5, 40) || !finite(r.setpoint, 10, 35) || typeof r.equipment !== "boolean" || !finite(r.roofU, 0, 5) || !r.walls) throw new Error(`Invalid room ${r?.id ?? ""}.`);
    const walls = Object.fromEntries(SIDES.map(s => {
      const w = r.walls[s], length = s === "north" || s === "south" ? r.width : r.depth;
      if (!w || typeof w.exterior !== "boolean" || !finite(w.U, .05, 6) || !finite(w.window, 0, length * r.height * .95)) throw new Error(`${r.name}: ${s} wall U-value and window area (at most 95% of the wall).`);
      return [s, { exterior: w.exterior, U: w.U, window: w.exterior ? w.window : 0 }];
    }));
    return { id: r.id, name: r.name, x: r.x, y: r.y, width: r.width, depth: r.depth, height: r.height, people: r.people, lighting: r.lighting, mode: r.mode,
      supply: r.supply, supplyTemp: r.supplyTemp, setpoint: r.setpoint, equipment: r.equipment, roofU: r.roofU, walls };
  });
  return { site: { name: d.site.name, outdoor: d.site.outdoor, solar: d.site.solar }, glazing: { U: d.glazing.U, shgc: d.glazing.shgc }, units: d.units, limit: d.limit, equipmentPower: d.equipmentPower ?? 0, rooms };
}
export const AIR = { rho: 1.2, cp: 1006 }; // supply air at ~20 °C
export function roomLoads(d, r, equipmentW = 0) {
  const floor = r.width * r.depth;
  let UA = r.roofU * floor, solar = 0, glass = 0, wallArea = 0;
  for (const s of SIDES) {
    const w = r.walls[s]; if (!w.exterior) continue;
    const area = (s === "north" || s === "south" ? r.width : r.depth) * r.height;
    UA += w.U * (area - w.window) + d.glazing.U * w.window; solar += w.window * d.glazing.shgc * d.site.solar * ORIENTATION[s];
    glass += w.window; wallArea += area;
  }
  const internal = r.people * SENSIBLE_PER_PERSON + r.lighting * floor + (r.equipment ? equipmentW : 0);
  const C = AIR.rho * AIR.cp * r.supply / 3600;
  const temperature = r.mode === "supplied" ? (C * r.supplyTemp + UA * d.site.outdoor + internal + solar) / (C + UA) : r.setpoint;
  // Load the supply must remove at the set point; for a supplied room, the
  // load its actual supply removes at the floating temperature.
  const cooling = r.mode === "supplied" ? C * (temperature - r.supplyTemp) : UA * (d.site.outdoor - r.setpoint) + internal + solar;
  // ASHRAE 62.1 office breathing-zone outdoor air: 2.5 L/s per person + 0.3 L/s per m².
  const ventilation = 2.5 * r.people + .3 * floor;
  const capacity = C * (r.setpoint - r.supplyTemp);
  return { floor, UA, solar, internal, temperature, cooling, ventilation, supplyLs: r.supply / 3.6, glass, wallArea, capacity };
}
export function architectureResults(d, overrides = {}) {
  const perUnit = overrides.equipment_power ?? d.equipmentPower ?? 0, equipmentW = perUnit * d.units;
  const rooms = d.rooms.map(r => ({ id: r.id, name: r.name, mode: r.mode, ...roomLoads(d, r, equipmentW) }));
  const equipmentRoom = rooms.find((r, i) => d.rooms[i].equipment) || rooms[0];
  return { rooms, equipmentRoom, equipmentW, floor: rooms.reduce((s, r) => s + r.floor, 0), UA: rooms.reduce((s, r) => s + r.UA, 0),
    cooling: rooms.reduce((s, r) => s + Math.max(0, r.cooling), 0) };
}
export function overlaps(d) {
  const out = [];
  for (let i = 0; i < d.rooms.length; i++) for (const b of d.rooms.slice(i + 1)) {
    const a = d.rooms[i];
    if (a.x < b.x + b.width - 1e-9 && b.x < a.x + a.width - 1e-9 && a.y < b.y + b.depth - 1e-9 && b.y < a.y + a.depth - 1e-9) out.push([a, b]);
  }
  return out;
}
export function architectureAssumptions(d) {
  return [
    { id: "architecture.climate", label: `Summer design ${d.site.outdoor} °C, ${d.site.solar} W/m² peak vertical irradiance`, kind: "nominal", source: "Use the ASHRAE climatic design data for your site." },
    { id: "architecture.orientation", label: "Façade irradiance factors N 0.25 · E 0.8 · S 0.75 · W 1.0", kind: "nominal", source: "Northern mid-latitude afternoon peak; not a solar-position calculation." },
    { id: "architecture.envelope", label: `Wall U-values and glazing U ${d.glazing.U} W/m²K, SHGC ${d.glazing.shgc}`, kind: "nominal", source: "Typical modern construction; use the specified assemblies." },
    { id: "architecture.steady", label: "Steady state: no thermal mass, no interior-wall transfer", kind: "nominal", source: "Peak-hour heat balance, not a dynamic simulation." },
    { id: "architecture.people", label: `${SENSIBLE_PER_PERSON} W sensible per occupant`, kind: "nominal", source: "ASHRAE Fundamentals, seated light work." },
  ];
}
