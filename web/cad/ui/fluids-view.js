// Fluids studio: the project's air path, its elements and its operating point.
import { project, onSolved } from "../app.js";
import { FANS, fanPressure, elementDrop, air } from "../domains/fluids.js";
import { node, svg, button, format, numberField, selectField, checkCards, certificates, assumptionList, section, linkedChip, live, lineChart, run } from "./kit.js";
const root = document.getElementById("fluids-studio");
let solution = null;
function render() {
  const d = project.store.get("fluids"), r = solution?.results.fluids, x = r?.details, op = x?.op, status = node("p", "", "studio-status"); status.id = "fluids-status";
  root.replaceChildren();
  const head = node("div", undefined, "studio-heading"), title = node("div");
  title.append(node("p", `FLUIDS / ${project.store.docs.project.name.toUpperCase()}`, "panel-label"), node("h2", "Where the air goes, and why."),
    node("p", "A series air path solved at the point where the driving pressure meets the system curve. Laminar friction uses Shah and London's exact rectangular-duct fit; turbulent friction uses Colebrook. Between Re 2300 and 4000 the flow is bounded by both laws, not interpolated."));
  head.append(title); root.append(head);
  if (op && !op.error) {
    const sink = op.elements.find(e => e.kind === "channels"), hero = node("div", undefined, "system-hero");
    for (const [label, v] of [["Volume flow", `${(op.Q * 1000).toFixed(3)} L/s${op.high - op.low > 1e-9 ? ` (${(op.low * 1000).toFixed(3)}–${(op.high * 1000).toFixed(3)})` : ""}`], ["Driving pressure", format(op.pressure, "Pa", 4)],
      ["Sink approach air", format(sink?.approach, "m/s", 4)], ["Channel Reynolds", sink?.re.toFixed(0)], ["Air density", `${op.props.rho.toFixed(4)} kg/m³`], ["Air temperature", format(x.temperature, "degC", 3)],
      ["Rust graph sink Δp", format(solution.outputs["mech.sink_dp"]?.value, "Pa", 4)], ["Duct model sink Δp", format((sink.lo + sink.hi) / 2, "Pa", 4)]]) { const c = node("div", undefined, "hero-stat"); c.append(node("span", label), node("b", v ?? "—")); hero.append(c); }
    root.append(hero);
    root.append(node("p", "The two sink pressure drops come from independent models: the Rust project graph treats the fins as parallel plates (Dh = 2g) with developing-flow corrections; the duct model uses the finite-aspect channel (Dh = 2gh/(g+h)) with entry and exit losses. They are shown side by side, not reconciled.", "small muted"));
    // Schematic with per-element pressure drop.
    const sch = section(root, "Air path", `${d.elements.length} elements in series`), max = Math.max(...op.elements.map(e => (e.lo + e.hi) / 2), 1e-9), path = node("div", undefined, "air-path");
    for (const e of op.elements) {
      const card = node("article", undefined, `air-element air-${e.kind}`); card.dataset.regime = e.regime;
      const bar = node("i"); bar.style.height = `${8 + 80 * (e.lo + e.hi) / 2 / max}%`;
      card.append(node("strong", e.name), node("span", `${format((e.lo + e.hi) / 2, "Pa", 3)}${e.hi - e.lo > 1e-9 ? ` [${e.lo.toFixed(2)}–${e.hi.toFixed(2)}]` : ""}`), node("small", `${e.regime}${e.re ? ` · Re ${e.re.toFixed(0)}` : ""} · ${format(e.velocity, "m/s", 3)}`), node("div", undefined, "air-bar"));
      card.lastChild.append(bar); path.append(card);
    }
    sch.append(path);
    // Driving pressure against the system curve.
    const pts = [], sys = [], hi = [], top = d.fan.preset === "wake" ? op.Q * 2 : d.fan.free;
    for (let i = 1; i <= 60; i++) { const Q = top * i / 60, p = air(x.temperature, d.pressure); pts.push([Q * 1000, d.fan.preset === "wake" ? op.pressure : Math.max(0, fanPressure(d.fan, Q))]);
      const drops = d.elements.map(e => elementDrop(e, Q, p)); sys.push([Q * 1000, drops.reduce((s, e) => s + e.lo, 0)]); hi.push([Q * 1000, drops.reduce((s, e) => s + e.hi, 0)]); }
    const chart = section(root, "Operating point", d.fan.preset === "wake" ? "constant wake pressure source" : FANS[d.fan.preset]?.title || "custom fan curve");
    chart.append(lineChart({ series: [{ points: pts, cls: "fan" }, { points: sys, cls: "system" }, { points: hi, cls: "system-hi" }], x: { label: "flow, L/s", min: 0 }, y: { label: "pressure, Pa", min: 0 }, marks: [{ x: op.Q * 1000, y: op.pressure, label: "operating point" }] }));
  } else if (op?.error) root.append(node("p", op.error, "conflict"));
  // Elements.
  const el = section(root, "Elements");
  const table = node("table", undefined, "data-table"), thead = node("thead"), tr = node("tr"); for (const h of ["Element", "Kind", "Width × height", "Length", "Loss data", ""]) tr.append(node("th", h)); thead.append(tr); table.append(thead);
  const body = node("tbody");
  for (const e of d.elements) {
    const row = node("tr"), num = (key, step, label) => { const i = node("input"); i.type = "number"; i.step = step; i.value = e[key]; i.id = `fluids-${e.id}-${key}`; i.setAttribute("aria-label", `${e.name} ${label}`);
      i.addEventListener("change", () => run("fluids.update", { id: e.id, [key]: Number(i.value) }, status).catch(() => { i.value = e[key]; })); return i; };
    row.append(node("td", e.name), node("td", e.kind));
    const wh = node("td"); if (e.kind === "channels") wh.append(node("span", `${e.fins - 1} × ${format(e.gap, "mm", 3)} × ${format(e.height, "mm", 3)}`)); else wh.append(num("width", 1, "width"), node("span", " × "), num("height", 1, "height")); row.append(wh);
    const len = node("td"); if (e.kind === "duct") len.append(num("length", 1, "length")); else len.textContent = e.kind === "channels" ? format(e.length, "mm", 3) : "—"; row.append(len);
    const loss = node("td"); if (e.kind === "fitting") loss.append(node("span", "K "), num("K", .05, "loss coefficient")); else if (e.kind === "filter") loss.append(num("coefficient", 1, "filter coefficient"), node("small", " Pa·s/m")); else if (e.kind === "duct") loss.append(node("span", "ε "), num("roughness", .0005, "roughness")); else loss.append(node("span", "K in "), num("Kin", .05, "entry loss"), node("span", " out "), num("Kout", .05, "exit loss")); row.append(loss);
    const del = node("td"); if (e.kind !== "channels") del.append(button("Remove", () => run("fluids.remove", { id: e.id }, status).catch(() => {}))); else del.append(node("small", "⛓ geometry from Mechanical")); row.append(del);
    body.append(row);
  }
  table.append(body); el.append(table);
  const add = node("div", undefined, "rtl-actions");
  add.append(button("＋ Duct", () => run("fluids.add", { kind: "duct" }, status).catch(() => {})), button("＋ Fitting", () => run("fluids.add", { kind: "fitting" }, status).catch(() => {})), button("＋ Filter", () => run("fluids.add", { kind: "filter" }, status).catch(() => {})));
  el.append(add);
  const src = section(root, "Driving pressure and air");
  const f = node("div", undefined, "field-grid");
  selectField(f, { label: "Source", id: "fluids-fan", value: d.fan.preset, options: [...Object.entries(FANS).map(([k, v]) => [k, v.title]), ["custom", "Custom fan curve"]], op: "fluids.fan", args: v => ({ preset: v }), status });
  if (d.fan.preset === "custom") { numberField(f, { label: "Shut-off pressure", id: "fluids-shutoff", unit: "Pa", value: d.fan.shutoff, min: 1, max: 5000, op: "fluids.fan", args: v => ({ preset: "custom", shutoff: v }), status }); numberField(f, { label: "Free delivery", id: "fluids-free", unit: "m3/s", value: d.fan.free, min: 1e-5, max: 10, step: .0001, op: "fluids.fan", args: v => ({ preset: "custom", free: v }), status }); }
  const t = numberField(f, { label: "Air temperature", id: "fluids-temperature", unit: "degC", value: d.temperature, min: -40, max: 150, op: "fluids.set", args: v => ({ temperature: v }), status });
  const chip = linkedChip(solution, "fluids.air_temperature"); if (chip) { t.disabled = true; t.parentElement.append(chip); }
  numberField(f, { label: "Air temperature rise limit", id: "fluids-limit", unit: "K", value: d.heatLimit, min: .1, max: 200, op: "fluids.set", args: v => ({ heat_limit: v }), status });
  src.append(f); const w = linkedChip(solution, "fluids.drive_pressure"); if (w && d.fan.preset === "wake") src.append(w); src.append(status);
  if (r) { const c = section(root, "Fluid checks"); checkCards(c, r.checks || []); certificates(c, r.certificates); }
  const a = section(root, "Assumptions"); assumptionList(a, "fluids");
}
onSolved(s => { solution = s; });
live(["fluids"], render);
render();
