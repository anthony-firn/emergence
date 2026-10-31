// Equations studio: a project's own relations. Inputs can be driven by links,
// and each formula is an output any other domain can link to.
import { project, onSolved } from "../app.js";
import { UNITS } from "../core/units.js";
import { FUNCTION_NAMES } from "../core/expr.js";
import { node, button, format, section, run, live, linkedChip, assumptionList, miniForm, certificates, checkCards } from "./kit.js";
const root = document.getElementById("equations-studio");
let solution = null;
const units = () => { const l = node("datalist"); l.id = "equation-units"; for (const u of Object.keys(UNITS).filter(Boolean)) l.append(Object.assign(node("option"), { value: u })); return l; };
function render() {
  const d = project.store.get("equations"), r = solution?.results.equations, values = r?.details?.scope ?? {}, status = node("p", "", "studio-status"); status.id = "equations-status";
  root.replaceChildren();
  const head = node("div", undefined, "studio-heading"), title = node("div");
  title.append(node("p", `EQUATIONS / ${project.store.docs.project.name.toUpperCase()}`, "panel-label"), node("h2", "Say what no domain knows."),
    node("p", "Name your own inputs and write formulas over them. Each formula is an output port, and each input can be driven by a link, so this is how to turn a die power and a rail voltage into a current, or add a cost or a margin. Units are as you declare them; the arithmetic inside a formula is not dimension-checked."));
  head.append(title); root.append(head, units());
  const i = section(root, "Inputs", `${d.inputs.length} named values`), table = node("table", undefined, "data-table"), thead = node("thead"), tr = node("tr");
  for (const h of ["Name", "Label", "Value", "Unit", ""]) tr.append(node("th", h)); thead.append(tr); table.append(thead); const body = node("tbody");
  for (const x of d.inputs) {
    const row = node("tr"), v = node("input"), td = node("td"), del = node("td"); v.type = "number"; v.step = "any"; v.value = x.value; v.id = `equations-value-${x.id}`; v.setAttribute("aria-label", `${x.id} value`);
    v.addEventListener("change", () => run("equations.input", { ...x, value: Number(v.value) }, status).catch(() => { v.value = x.value; })); td.append(v);
    const chip = linkedChip(solution, `equations.${x.id}`); if (chip) { v.disabled = true; td.append(chip); }
    del.append(button("Remove", () => run("equations.remove", { id: x.id }, status).catch(() => {})));
    row.append(node("td", x.id), node("td", x.label), td, node("td", x.unit || "—"), del); body.append(row);
  }
  if (!d.inputs.length) { const row = node("tr"), td = node("td", "No inputs yet."); td.colSpan = 5; row.append(td); body.append(row); }
  table.append(body); i.append(table);
  miniForm(i, { id: "equations-input", action: "＋ Input", status, fields: [{ name: "id", label: "Name", placeholder: "rail_voltage" }, { name: "value", label: "Value", number: true }, { name: "unit", label: "Unit", list: "equation-units" }, { name: "label", label: "Label" }, { name: "min", label: "Lowest a link may set", number: true }, { name: "max", label: "Highest a link may set", number: true }],
    submit: a => run("equations.input", a, status) });
  const o = section(root, "Formulas", `Functions: ${FUNCTION_NAMES.join(", ")}; constants pi and e`), t2 = node("table", undefined, "data-table"), th2 = node("thead"), tr2 = node("tr");
  for (const h of ["Name", "Formula", "Value", ""]) tr2.append(node("th", h)); th2.append(tr2); t2.append(th2); const b2 = node("tbody");
  for (const x of d.outputs) {
    const row = node("tr"), e = node("input"), td = node("td"), del = node("td"); e.value = x.expr; e.id = `equations-expr-${x.id}`; e.setAttribute("aria-label", `${x.id} formula`); e.spellcheck = false;
    e.addEventListener("change", () => run("equations.output", { ...x, expr: e.value }, status).catch(() => { e.value = x.expr; })); td.append(e);
    del.append(button("Remove", () => run("equations.remove", { id: x.id }, status).catch(() => {})));
    row.append(node("td", x.label === x.id ? x.id : `${x.id} · ${x.label}`), td, node("td", format(values[x.id], x.unit, 5)), del); row.lastChild.previousSibling.id = `equations-out-${x.id}`; b2.append(row);
  }
  if (!d.outputs.length) { const row = node("tr"), td = node("td", "No formulas yet."); td.colSpan = 4; row.append(td); b2.append(row); }
  t2.append(b2); o.append(t2);
  miniForm(o, { id: "equations-output", action: "＋ Formula", status, fields: [{ name: "id", label: "Name", placeholder: "rail_current" }, { name: "expr", label: "Formula", placeholder: "power / rail_voltage", wide: true }, { name: "unit", label: "Unit", list: "equation-units" }, { name: "label", label: "Label" }],
    submit: a => run("equations.output", a, status) });
  root.append(status);
  if (r?.checks?.length) { const c = section(root, "Equation checks"); checkCards(c, r.checks); certificates(c, r.certificates); }
  const a = section(root, "Assumptions"); assumptionList(a, "equations");
}
onSolved(s => { solution = s; });
live(["equations"], render);
render();
