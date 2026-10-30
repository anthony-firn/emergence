// The open project: one command bus, one document store and one system
// solver for every studio on the page. Which project is open comes from
// projects.js; the mechanical workspace fills in `mechanical` with accessors
// onto its revisioned state and makes the project's first revision.
import { createProject, DOMAIN_INFO } from "./domains/index.js";
import { request } from "./models.js";
export { current } from "./projects.js";
export const mechanical = { state: () => ({ params: {}, features: [], drawings: {} }), snapshot: () => ({}), restore() {}, set() {} };
const fetchJSON = async url => { const r = await fetch(url); if (!r.ok) throw new Error(`${url}: ${r.status}`); return r.json(); };
const text = name => fetch(new URL(`./data/${name}`, import.meta.url)).then(r => r.ok ? r.text() : "").catch(() => "");
const verilog = { mini8: await text("mini8.v"), "tt-motor": await text("tt_motor_controller.v") };
const { modelMesh } = await import("./ui/meshes.js");
export const project = createProject({ request, mech: mechanical, verilog, mesh: (model, overrides) => modelMesh(model, overrides), loadCorners: () => fetchJSON(new URL("./data/chip-corners.json", import.meta.url)) });
export { DOMAIN_INFO };
// Debounced, generation-checked system solves. Studios listen for results.
const solvedListeners = new Set(), changed = new Set();
let timer, generation = 0, running = null;
export let solution = null;
export const onSolved = fn => { solvedListeners.add(fn); if (solution) fn(solution); return () => solvedListeners.delete(fn); };
export const onSolving = new Set();
export function scheduleSolve(domains = []) {
  domains.forEach(d => changed.add(d)); clearTimeout(timer);
  for (const fn of onSolving) fn(true);
  timer = setTimeout(async () => {
    const gen = ++generation, edited = [...changed]; changed.clear();
    await running?.catch(() => {});
    if (gen !== generation) return;
    running = project.solve({ changed: edited });
    try { const s = await running; if (gen !== generation) return; solution = s; for (const fn of solvedListeners) fn(s); }
    catch (e) { console.error(e); for (const fn of onSolving) fn(false, e); }
  }, 140);
}
project.store.on(({ ids }) => scheduleSolve(ids));
window.emergence = Object.freeze({
  // The browser command API: the same registry the CLI and MCP server use.
  commands: domain => project.bus.list(domain),
  run: (op, args = {}) => project.bus.run(op, args, { source: "api" }),
  batch: list => project.bus.batch(list, { source: "api" }),
  parse: text => project.bus.parse(text),
  solve: () => project.solve(),
  get solution() { return solution; },
  documents: () => ({ ...project.store.all(), mech: structuredClone(mechanical.state()) }),
  // The open design as a project file: the same format the CLI and MCP server read.
  file: () => project.file(),
  undo: () => project.bus.undo(), redo: () => project.bus.redo(),
});
