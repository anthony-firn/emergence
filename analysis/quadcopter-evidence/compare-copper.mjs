import { readFileSync } from "node:fs";
const S = (process.env.EVIDENCE_WORK ?? "/tmp/emergence-evidence") + "/";
const { kernelBoard, parseGerber } = await import(new URL("../../web/cad/domains/pcb.js", import.meta.url));
const { pads } = await import(new URL("../../web/cad/studio-model.js", import.meta.url));
const { nodeModels } = await import(new URL("../../web/cad/agent/models.mjs", import.meta.url));
const { exampleFile } = await import(new URL("../../web/cad/agent/headless.mjs", import.meta.url));
const d = (await exampleFile("drone")).documents.pcb, models = await nodeModels(), k = kernelBoard(d), kc = JSON.parse(readFileSync(S + "kicad/copper.json", "utf8"));
const { result: r } = await models.request("pcb", { input: { board: k.board, currents: k.currents, fabrication: true, cell: 5e-5 } });
const cell = .05, nx = Math.ceil(d.width / cell), ny = Math.ceil(d.height / cell);
// Even-odd scan conversion of rings given in board coordinates, y down.
function raster(rings) { const rows = Array.from({ length: ny }, () => []), g = new Uint8Array(nx * ny);
  for (const l of rings) for (let i = 0, j = l.length - 1; i < l.length; j = i++) { const [xi, yi] = l[i], [xj, yj] = l[j]; for (let r = Math.max(0, Math.floor(Math.min(yi, yj) / cell - 1)); r <= Math.min(ny - 1, Math.ceil(Math.max(yi, yj) / cell)); r++) { const y = (r + .5) * cell; if ((yi > y) !== (yj > y)) rows[r].push(xi + (xj - xi) * (y - yi) / (yj - yi)); } }
  rows.forEach((c, r) => { c.sort((a, b) => a - b); for (let i = 0; i + 1 < c.length; i += 2) for (let x = Math.max(0, Math.ceil(c[i] / cell - .5)); x <= Math.min(nx - 1, Math.floor(c[i + 1] / cell - .5)); x++) g[r * nx + x] = 1; }); return g; }
const holes = pads(d).filter(p => !p.smd), scores = {};
for (const l of r.gerber) {
  const name = l.name === "top" ? "top" : "bottom", app = raster(parseGerber(l.gerber).map(loop => loop.map(([x, y]) => [x, d.height - y]))), kicad = raster(kc[name].flat());
  let both = 0, onlyApp = 0, onlyKi = 0;
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const x = (i + .5) * cell, y = (j + .5) * cell; if (holes.some(p => Math.abs(x - p.x) < 2 && Math.abs(y - p.y) < 2 && Math.hypot(x - p.x, y - p.y) < (p.drill ?? .9) / 2)) continue; const a = app[j * nx + i], b = kicad[j * nx + i]; if (a && b) both++; else if (a) onlyApp++; else if (b) onlyKi++; }
  const mm = n => (n * cell * cell).toFixed(1);
  scores[name] = 100 * both / (both + onlyApp + onlyKi);
  console.log(name, "app Gerber", mm(both + onlyApp), "mm²; KiCad", mm(both + onlyKi), "mm²; overlap", (100 * both / (both + onlyApp + onlyKi)).toFixed(2) + "%", "only app", mm(onlyApp), "only KiCad", mm(onlyKi));
}
(await import("node:fs")).writeFileSync(S + "kicad/copper-compare.json", JSON.stringify(scores)); process.exit(0);
