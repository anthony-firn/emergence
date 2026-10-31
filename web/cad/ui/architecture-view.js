// Architecture studio: the plan and each room's heat balance.
import { project, onSolved } from "../app.js";
import { SIDES, ORIENTATION, BESTEST } from "../domains/architecture.js";
import { node, svg, button, format, numberField, selectField, textField, checkCards, certificates, assumptionList, section, linkedChip, live, run } from "./kit.js";
const root = document.getElementById("architecture-studio");
let solution = null, selected = "charge";
function plan(d, x) {
  const xs = d.rooms.flatMap(r => [r.x, r.x + r.width]), ys = d.rooms.flatMap(r => [r.y, r.y + r.depth]), x0 = Math.min(...xs) - 2, y0 = Math.min(...ys) - 2, w = Math.max(...xs) - x0 + 2, h = Math.max(...ys) - y0 + 2;
  const s = svg("svg", { viewBox: `${x0} ${y0} ${w} ${h}`, class: "arch-plan", role: "img", "aria-label": "Floor plan in metres, north up" });
  const max = Math.max(...(x?.rooms || []).map(r => r.temperature), 1), min = Math.min(...(x?.rooms || []).map(r => r.temperature), 0);
  for (const r of d.rooms) {
    const z = x?.rooms.find(q => q.id === r.id), t = z ? (z.temperature - min) / Math.max(1e-9, max - min) : 0;
    const g = svg("g", { class: "arch-room", role: "button", tabindex: 0, "aria-pressed": String(selected === r.id), "aria-label": `${r.name}${z ? `, ${z.temperature.toFixed(1)} °C` : ""}` });
    g.append(svg("rect", { x: r.x, y: r.y, width: r.width, height: r.depth, style: `fill:hsl(${Math.round(210 - 190 * t)} 60% 50% / .18)` }));
    for (const side of SIDES) {
      const wl = r.walls[side], [ax, ay, bx, by] = { north: [r.x, r.y, r.x + r.width, r.y], south: [r.x, r.y + r.depth, r.x + r.width, r.y + r.depth], west: [r.x, r.y, r.x, r.y + r.depth], east: [r.x + r.width, r.y, r.x + r.width, r.y + r.depth] }[side];
      g.append(svg("line", { x1: ax, y1: ay, x2: bx, y2: by, class: wl.exterior ? "wall-ext" : "wall-int" }));
      if (wl.window > 0) { const len = Math.hypot(bx - ax, by - ay), frac = Math.min(.9, wl.window / (len * r.height)), mx = (ax + bx) / 2, my = (ay + by) / 2, hx = (bx - ax) / 2 * frac, hy = (by - ay) / 2 * frac;
        g.append(svg("line", { x1: mx - hx, y1: my - hy, x2: mx + hx, y2: my + hy, class: "wall-window" })); }
    }
    g.append(svg("text", { x: r.x + r.width / 2, y: r.y + r.depth / 2 - .3, "text-anchor": "middle", class: "arch-label" }, r.name), svg("text", { x: r.x + r.width / 2, y: r.y + r.depth / 2 + .9, "text-anchor": "middle", class: "arch-sub" }, z ? `${z.temperature.toFixed(1)} °C · ${r.mode}` : r.mode));
    g.addEventListener("click", () => { selected = r.id; render(); }); g.addEventListener("keydown", e => { if (e.key === "Enter") { selected = r.id; render(); } });
    s.append(g);
  }
  s.append(svg("text", { x: x0 + 1, y: y0 + 1.4, class: "arch-sub" }, "N ↑"));
  return s;
}
function render() {
  const d = project.store.get("architecture"), r = solution?.results.architecture, x = r?.details, status = node("p", "", "studio-status"); status.id = "architecture-status";
  root.replaceChildren();
  const head = node("div", undefined, "studio-heading"), title = node("div");
  title.append(node("p", "ARCHITECTURE / " + d.site.name.toUpperCase(), "panel-label"), node("h2", "The building is part of the system."),
    node("p", `Each room's steady heat balance: envelope conductance, solar gain through glazing, people, lighting, equipment and supply air. Equipment heat per unit arrives through links from the other domains; ${d.units} unit${d.units === 1 ? "" : "s"} share the equipment room. ${d.rooms.some(r => /BESTEST/.test(r.name)) ? `Rooms marked BESTEST 600 use ${BESTEST}.` : ""}`));
  head.append(title); root.append(head);
  const grid = node("div", undefined, "domain-grid");
  const p = section(grid, "Plan", "metres · north up · click a room");
  const room = d.rooms.find(q => q.id === selected) || d.rooms[0], z = x?.rooms.find(q => q.id === room?.id);
  if (room) p.append(plan(d, x)); else p.append(node("p", "No rooms yet. Add the first one below.", "small muted"));
  const insp = section(grid, room?.name ?? "Rooms", z ? `${z.temperature.toFixed(2)} °C` : "");
  if (room) {
  if (z) { const stats = node("dl", undefined, "stat-grid"); for (const [t, v] of [["Envelope UA", format(z.UA, "W/K", 4)], ["Solar gain", format(z.solar, "W", 4)], ["Internal gains", format(z.internal, "W", 4)], ["Supply air", `${(room.supply / 3.6).toFixed(0)} L/s at ${room.supplyTemp} °C`], ["Cooling removed", format(z.cooling, "W", 4)], ["Outdoor air required (62.1)", `${z.ventilation.toFixed(1)} L/s`]]) { const div = node("div"); div.append(node("dt", t), node("dd", v)); stats.append(div); } insp.append(stats); }
  const f = node("div", undefined, "field-grid");
  textField(f, { label: "Name", id: "arch-name", value: room.name, op: "architecture.update_room", args: v => ({ id: room.id, name: v }), status });
  selectField(f, { label: "Mode", id: "arch-mode", value: room.mode, options: [["supplied", "Supplied (floats)"], ["conditioned", "Conditioned (set point)"]], op: "architecture.update_room", args: v => ({ id: room.id, mode: v }), status });
  numberField(f, { label: "Plan X", id: "arch-x", unit: "m", value: room.x, min: -500, max: 500, step: .1, op: "architecture.update_room", args: v => ({ id: room.id, at: [v, room.y] }), status });
  numberField(f, { label: "Plan Y", id: "arch-y", unit: "m", value: room.y, min: -500, max: 500, step: .1, op: "architecture.update_room", args: v => ({ id: room.id, at: [room.x, v] }), status });
  for (const [key, label, unit, min, max, step, value] of [["width", "Width", "m", .5, 200, .1, room.width], ["depth", "Depth", "m", .5, 200, .1, room.depth], ["height", "Height", "m", 2, 20, .1, room.height], ["people", "Occupants", "", 0, 2000, 1, room.people],
    ["lighting", "Lighting", "W/m2", 0, 50, .5, room.lighting], ["supply", "Supply air", "m3/h", 0, 1e6, 10, room.supply], ["supply_temp", "Supply temperature", "degC", 5, 40, .5, room.supplyTemp], ["setpoint", "Set point", "degC", 10, 35, .5, room.setpoint], ["roof_u", "Roof U", "W/m2K", 0, 5, .01, room.roofU]])
    numberField(f, { label, id: `arch-${key}`, unit, value, min, max, step, op: "architecture.update_room", args: v => ({ id: room.id, [key]: v }), status });
  insp.append(f);
  const walls = node("table", undefined, "data-table"), wh = node("thead"), wr = node("tr"); for (const h of ["Wall", "Exterior", "U (W/m²K)", "Window (m²)", "Solar factor"]) wr.append(node("th", h)); wh.append(wr); walls.append(wh);
  const wb = node("tbody");
  for (const side of SIDES) { const w = room.walls[side], tr = node("tr"); tr.append(node("td", side));
    const ex = node("input"); ex.type = "checkbox"; ex.checked = w.exterior; ex.setAttribute("aria-label", `${side} wall exterior`); ex.addEventListener("change", () => run("architecture.wall", { id: room.id, side, exterior: ex.checked }, status).catch(() => {}));
    const u = node("input"); u.type = "number"; u.step = .01; u.value = w.U; u.setAttribute("aria-label", `${side} wall U-value`); u.addEventListener("change", () => run("architecture.wall", { id: room.id, side, u: Number(u.value) }, status).catch(() => {}));
    const win = node("input"); win.type = "number"; win.step = .5; win.min = 0; win.value = w.window; win.disabled = !w.exterior; win.setAttribute("aria-label", `${side} window area`); win.id = `arch-window-${side}`; win.addEventListener("change", () => run("architecture.wall", { id: room.id, side, window: Number(win.value) }, status).catch(() => {}));
    for (const c of [ex, u, win]) { const td = node("td"); td.append(c); tr.append(td); } tr.append(node("td", String(ORIENTATION[side]))); wb.append(tr); }
  walls.append(wb); insp.append(walls);
  const eq = node("label", undefined, "check-label"), ec = node("input"); ec.type = "checkbox"; ec.checked = room.equipment; ec.id = "arch-equipment"; ec.addEventListener("change", () => run("architecture.update_room", { id: room.id, equipment: ec.checked }, status).catch(() => {}));
  eq.append(ec, document.createTextNode("Equipment room (receives the equipment heat)")); insp.append(eq, button("Remove room", () => run("architecture.remove_room", { id: room.id }, status).catch(() => {})));
  }
  insp.append(status);
  root.append(grid);
  const site = section(root, "Site, envelope and equipment");
  const sf = node("div", undefined, "field-grid");
  for (const [key, label, unit, min, max, step, value] of [["outdoor", "Outdoor design temperature", "degC", -40, 55, .5, d.site.outdoor], ["solar", "Peak vertical irradiance", "W/m2", 0, 1200, 10, d.site.solar], ["glazing_u", "Glazing U", "W/m2K", .3, 7, .1, d.glazing.U], ["shgc", "Glazing SHGC", "", .05, .9, .01, d.glazing.shgc], ["units", "Chargers running", "", 0, 10000, 1, d.units], ["limit", "Equipment room limit", "degC", 5, 80, .5, d.limit]])
    numberField(sf, { label, id: `arch-site-${key}`, unit, value, min, max, step, op: "architecture.site", args: v => ({ [key]: v }), status });
  textField(sf, { label: "Site name", id: "arch-site-name", value: d.site.name, op: "architecture.site", args: v => ({ name: v }), status });
  const eqp = numberField(sf, { label: "Heat per equipment unit", id: "arch-site-equipment", unit: "W", value: d.equipmentPower, min: 0, max: 1e6, step: 1, op: "architecture.site", args: v => ({ equipment: v }), status });
  if (solution?.applied["architecture.equipment_power"]) eqp.disabled = true;
  site.append(sf);
  const chip = linkedChip(solution, "architecture.equipment_power"); if (chip) site.append(chip);
  site.append(button("＋ Add room", () => run("architecture.add_room", { label: "New room" }, status).then(r => { selected = r.id; render(); }).catch(() => {})));
  if (x) { const table = node("table", undefined, "data-table"), th = node("thead"), tr = node("tr"); for (const h of ["Room", "Mode", "Air", "UA", "Solar", "Internal", "Cooling"]) tr.append(node("th", h)); th.append(tr); table.append(th);
    const tb = node("tbody"); for (const q of x.rooms) { const row = node("tr"); for (const v of [q.name, q.mode, format(q.temperature, "degC", 4), format(q.UA, "W/K", 3), format(q.solar, "W", 3), format(q.internal, "W", 3), format(q.cooling, "W", 4)]) row.append(node("td", v)); tb.append(row); }
    table.append(tb); site.append(table); }
  if (r) { const c = section(root, "Building checks"); checkCards(c, r.checks || []); certificates(c, r.certificates); }
  const a = section(root, "Assumptions"); assumptionList(a, "architecture");
}
onSolved(s => { solution = s; });
live(["architecture"], render);
render();
