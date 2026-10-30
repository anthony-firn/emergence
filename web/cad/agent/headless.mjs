// A project without a browser: the same registry, documents and solver, with
// the WebAssembly models called directly and the design kept in a project
// file. Examples are ordinary project files in web/cad/projects.
import { readFile, writeFile } from "node:fs/promises";
import { existsSync } from "node:fs";
import { nodeModels, loadCorners, loadVerilog } from "./models.mjs";
import { createProject } from "../domains/index.js";
import { MechanicalState } from "../domains/mechanical.js";
import { kicadPCB, kernelBoard } from "../domains/pcb.js";
import { gdsii, verifyGDS } from "../domains/silicon.js";
import { seal, unseal } from "../core/integrity.js";
import { DOMAINS, PROJECT_FORMAT, projectDefault, projectFile, slug } from "../core/project.js";
const LEGACY = "emergence-headless/1";
const PROJECTS = new URL("../projects/", import.meta.url);
export const examples = async () => JSON.parse(await readFile(new URL("index.json", PROJECTS), "utf8")).examples;
export async function exampleFile(id) {
  if (id === "blank") return projectFile({ project: projectDefault() });
  const all = await examples(), e = all.find(x => x.id === id);
  if (!e) throw new Error(`No example ${id}. Examples: blank, ${all.map(x => x.id).join(", ")}.`);
  return JSON.parse(await readFile(new URL(e.file, PROJECTS), "utf8"));
}
// A project file, sealed or not. Files from before projects had no manifest:
// they open with every domain in use.
export async function readFileAt(path) {
  const check = await unseal(JSON.parse(await readFile(path, "utf8")));
  if (!check.ok) throw new Error(`${path}: integrity hash does not match its content.`);
  const body = check.body;
  if (body.format === LEGACY) return projectFile({ project: { ...projectDefault(), name: "Imported project", domains: [...DOMAINS] }, ...body.docs, mech: body.mech });
  if (body.format !== PROJECT_FORMAT) throw new Error(`${path} is not an ${PROJECT_FORMAT} file.`);
  return body;
}
export async function openProject(path, { example } = {}) {
  const { request } = await nodeModels(), mech = new MechanicalState();
  const project = createProject({ request, mech, loadCorners, verilog: await loadVerilog() });
  if (path && existsSync(path)) project.load(await readFileAt(path));
  else if (example) project.load(await exampleFile(example));
  const save = async (to = path) => { if (to) await writeFile(to, JSON.stringify(await seal(project.file()), null, 1) + "\n"); };
  return { project, mech, request, save };
}
export function summary(s, project) {
  const out = Object.fromEntries(Object.entries(s.outputs).map(([k, v]) => [k, { value: Number.isFinite(v.value) ? Number(v.value.toPrecision(6)) : null, unit: v.unit }]));
  return { ...(project && { project: project.store.docs.project.name, domains: project.active() }), converged: s.converged, iterations: s.iterations, residual: s.residual, ms: Math.round(s.ms),
    checks: s.checks.map(c => ({ domain: c.domain, id: c.id, status: c.status, text: c.text, ...(c.reason && { reason: c.reason }) })),
    issues: s.issues.map(i => i.message), outputs: out };
}
export async function exportArtifact({ project, request }, kind) {
  const docs = project.store.docs, name = slug(docs.project.name), on = project.active();
  const need = domain => { if (!on.includes(domain)) throw new Error(`This project has no ${domain} domain.`); };
  if (kind === "kicad") { need("pcb"); return { name: `${name}.kicad_pcb`, text: kicadPCB(docs.pcb) }; }
  if (kind === "gerber") { need("pcb"); const k = kernelBoard(docs.pcb), { result } = await request("pcb", { input: { board: k.board, currents: k.currents, fabrication: true } });
    return { files: [...result.gerber.map(l => ({ name: `${name}-${l.name === "top" ? "F" : "B"}_Cu.gbr`, text: l.gerber })), { name: `${name}-Edge_Cuts.gbr`, text: result.profile }, { name: `${name}.drl`, text: result.excellon }] }; }
  if (kind === "gds") { need("silicon"); const bytes = gdsii(docs.silicon), v = verifyGDS(docs.silicon, bytes); if (!v.ok) throw new Error(v.problems.join(" ")); return { name: `${name}-floorplan.gds`, bytes, verified: v }; }
  if (kind === "project") return { name: `${name}.emergence.json`, text: JSON.stringify(await seal(project.file()), null, 1) + "\n" };
  throw new Error("Export kinds: kicad, gerber, gds, project.");
}
