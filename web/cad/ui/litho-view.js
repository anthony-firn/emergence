// Lithography stage inside the silicon studio: the layout clip, its Abbe
// aerial image, the threshold print and an inverse-lithography mask.
import { project } from "../app.js";
import { rasterize, aerial, printed, patternError, ilt } from "../domains/litho.js";
import { node, button, svg, numberField, section } from "./kit.js";
let cache = null;
function heat(canvas, n, values, mode, target) {
  canvas.width = n; canvas.height = n;
  const ctx = canvas.getContext("2d"), img = ctx.createImageData(n, n), max = Math.max(...values, 1e-9);
  for (let k = 0; k < n * n; k++) {
    const v = values[k] / (mode === "mask" ? 1 : max), o = 4 * k;
    if (mode === "mask") { img.data[o] = img.data[o + 1] = img.data[o + 2] = v > .5 ? 225 : 30; }
    else { img.data[o] = 255 * Math.min(1, 2 * v); img.data[o + 1] = 255 * Math.max(0, Math.min(1, 2 * v - .6)); img.data[o + 2] = 255 * Math.max(0, 1 - 3 * v) * .7 + 40 * v; }
    if (target && target[k] && mode !== "mask") { img.data[o + 2] = Math.min(255, img.data[o + 2] + 60); }
    img.data[o + 3] = 255;
  }
  ctx.putImageData(img, 0, 0);
}
export function mountLitho(root) {
  const d = project.store.get("litho"), s = section(root, "Lithography clip", `Abbe imaging · λ ${d.wavelength} nm · NA ${d.na} · σ ${d.sigmaOuter}`);
  s.classList.add("litho"); s.id = "silicon-litho";
  s.append(node("p", "A 150 nm layout clip imaged by scalar Abbe theory: each point of the partially coherent source tilts the mask spectrum through the pupil |f| ≤ NA/λ and the intensities add. A clear mask images to exactly 1 and a coherent grating past cutoff to its DC term squared; both are tested. Inverse lithography descends the exact adjoint gradient of the printed-pattern error with respect to every mask pixel.", "small muted"));
  const grid = node("div", undefined, "litho-grid"), panels = ["Target layout", "Aerial image", "Printed (threshold)", "ILT mask"].map(t => { const f = node("figure"); const c = node("canvas"); c.setAttribute("aria-label", t); f.append(c, node("figcaption", t)); grid.append(f); return c; });
  const status = node("p", "", "studio-status"); status.id = "litho-status";
  const controls = node("div", undefined, "litho-controls");
  numberField(controls, { label: "NA", id: "litho-na", value: d.na, min: .1, max: 1.5, step: .01, op: "litho.set", args: v => ({ na: v }), status });
  numberField(controls, { label: "σ outer", id: "litho-sigma", value: d.sigmaOuter, min: 0, max: 1, step: .05, op: "litho.set", args: v => ({ sigma_outer: v }), status });
  numberField(controls, { label: "Threshold", id: "litho-threshold", value: d.threshold, min: .05, max: .95, step: .01, op: "litho.set", args: v => ({ threshold: v }), status });
  numberField(controls, { label: "Wavelength", id: "litho-wavelength", unit: "nm", value: d.wavelength, min: 13, max: 500, step: 1, op: "litho.set", args: v => ({ wavelength: v }), status });
  const run = button("Run inverse lithography", async () => {
    run.disabled = true; status.textContent = "Optimising the mask…";
    const target = rasterize(d); let last = 0;
    const r = ilt(d, target, { iterations: 36, onStep: h => { last = h.errors; } });
    heat(panels[3], d.n, r.mask, "mask"); heat(panels[2], d.n, printed(d, r.I), "mask");
    status.textContent = `ILT: printed-pattern errors ${r.history[0].errors} → ${r.errors} pixels (of ${d.n * d.n}) in ${r.history.length} adjoint steps. Serifs and assist features emerge from the gradient; this is a study model, not OPC sign-off.`;
    run.disabled = false;
  }, "accent", "litho-ilt");
  controls.append(run);
  s.append(controls, grid, status);
  queueMicrotask(() => {
    const target = rasterize(d), { I } = aerial(d, target), print = printed(d, I), err = patternError(d, target, I);
    heat(panels[0], d.n, target, "mask"); heat(panels[1], d.n, I, "heat", target); heat(panels[2], d.n, print, "mask"); heat(panels[3], d.n, target, "mask");
    status.textContent = `Uncorrected mask: ${err} of ${d.n * d.n} pixels print wrong at threshold ${d.threshold} (field ${d.n * d.pixel} nm, ${d.pixel} nm pixels, k₁ = ${(150 * d.na / d.wavelength).toFixed(2)} for 150 nm lines).`;
  });
}
