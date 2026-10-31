// Cross-domain CAD: units, commands, merges, every domain model against its
// exact or published reference, project files, the linked solves of the
// example projects, revisions and the agent.
import test from "node:test";
import assert from "node:assert/strict";
import { spawn } from "node:child_process";
import { convert, compatible } from "./cad/core/units.js";
import { CommandBus } from "./cad/core/bus.js";
import { merge3 } from "./cad/core/merge.js";
import { System } from "./cad/core/system.js";
import { Repository, hash } from "./cad/repository.js";
import { seal, unseal, zip, crc32 } from "./cad/core/integrity.js";
import * as optics from "./cad/domains/optics.js";
import * as fluids from "./cad/domains/fluids.js";
import * as arch from "./cad/domains/architecture.js";
import * as flight from "./cad/domains/flight.js";
import * as litho from "./cad/domains/litho.js";
import { gdsii, verifyGDS, readGDS } from "./cad/domains/silicon.js";
import { autoroute, kicadPCB, kernelBoard, gerberCheck, densityIssues, pourLoss, pourSize } from "./cad/domains/pcb.js";
import { powerExample, validatePower, stageLoss, busRipple } from "./cad/domains/power.js";
import { structuresExample, validateStructures, armResults } from "./cad/domains/structures.js";
import { controlsExample, validateControls, plantOf, simulate, rateLoop, tune, flightTests } from "./cad/domains/controls.js";
import { correlate, evaluatePlan, digest } from "./cad/domains/testing.js";
import { pcbExample, siliconExample, pcbChecks, validateDocument, pads, pourIslands } from "./cad/studio-model.js";
import { nodeModels, loadCorners, loadVerilog } from "./cad/agent/models.mjs";
import { createProject } from "./cad/domains/index.js";
import { MechanicalState } from "./cad/domains/mechanical.js";
import { exampleFile, examples, readFileAt } from "./cad/agent/headless.mjs";
import { validateProject, getPath, setPath, numericPaths, DOMAINS } from "./cad/core/project.js";
import { writeFile, readFile, readdir } from "node:fs/promises";
import { buildScript } from "./cad/domains/rebuild.js";
import { parse, evaluate } from "./cad/core/expr.js";
import { validateEquations } from "./cad/domains/equations.js";
import { tmpdir } from "node:os";
import { join } from "node:path";
const close = (a, b, tol, msg) => assert(Math.abs(a - b) <= tol, `${msg ?? ""} ${a} vs ${b} (tol ${tol})`);
const models = await nodeModels();
const project = async () => createProject({ request: models.request, mech: new MechanicalState(), loadCorners, verilog: await loadVerilog() });
const opened = async id => { const p = await project(); p.load(await exampleFile(id)); return p; };

test("units convert absolutely, refuse mismatched dimensions", () => {
  close(convert(25, "degC", "K"), 298.15, 1e-12); close(convert(3, "L/s", "m3/h"), 10.8, 1e-12); close(convert(1, "kg", "g"), 1000, 1e-12);
  assert(!compatible("W", "Pa")); assert.throws(() => convert(1, "W", "A"), /Cannot convert/);
});
test("commands parse, validate, undo as one step per burst, and scripts are atomic", async () => {
  const bus = new CommandBus(); let doc = { a: 1, b: 1 };
  bus.domain("t", { snapshot: () => ({ ...doc }), restore: s => { doc = s; } });
  bus.register({ op: "t.set", domain: "t", title: "Set", coalesce: a => a.key, args: [{ name: "key", type: "enum", values: ["a", "b"], required: true }, { name: "value", type: "number", min: 0, max: 10, required: true }], run: ({ key, value }) => { doc[key] = value; } });
  bus.register({ op: "t.fail", domain: "t", title: "Fail", run: () => { doc.a = 99; throw new Error("no"); } });
  assert.deepEqual(bus.parse("t.set a 3"), { op: "t.set", args: { key: "a", value: "3" } });
  assert.deepEqual(bus.parse('set value=4 key=b', "t"), { op: "t.set", args: { value: "4", key: "b" } });
  await assert.rejects(bus.run("t.set", { key: "a", value: 11 }), /0–10/);
  for (const v of [2, 3, 4]) await bus.run("t.set", { key: "a", value: v });
  assert.equal(bus.undoStack.length, 1); bus.undo(); assert.equal(doc.a, 1); bus.redo(); assert.equal(doc.a, 4);
  await assert.rejects(bus.run("t.fail"), /no/); assert.equal(doc.a, 4);
  await assert.rejects(bus.batch([{ op: "t.set", args: { key: "b", value: 7 } }, { op: "t.fail" }]), /Step 2/); assert.deepEqual(doc, { a: 4, b: 1 });
});
test("three-way merge joins edits to different objects and reports both-sides conflicts", () => {
  const base = { w: 10, parts: [{ id: "a", x: 1 }, { id: "b", x: 2 }] };
  const ours = { w: 10, parts: [{ id: "a", x: 5 }, { id: "b", x: 2 }] }, theirs = { w: 12, parts: [{ id: "a", x: 1 }, { id: "b", x: 2 }, { id: "c", x: 3 }] };
  assert.deepEqual(merge3(base, ours, theirs).value, { w: 12, parts: [{ id: "a", x: 5 }, { id: "b", x: 2 }, { id: "c", x: 3 }] });
  const both = merge3(base, ours, { ...base, parts: [{ id: "a", x: 9 }, { id: "b", x: 2 }] });
  assert.equal(both.conflicts.length, 1); assert.equal(both.conflicts[0].path, "parts[a].x");
  assert.equal(merge3(base, ours, { ...base, parts: [{ id: "a", x: 9 }, { id: "b", x: 2 }] }, { "parts[a].x": "source" }).value.parts[0].x, 9);
  assert.deepEqual(merge3(base, base, { w: 10, parts: [{ id: "a", x: 1 }] }).value.parts, [{ id: "a", x: 1 }]);
});
test("optics: catalogue indices, lensmaker, window shift, Fresnel and the Cooke triplet", () => {
  close(optics.refractiveIndex("N-BK7", 587.56), 1.51680, 1e-5); close(optics.refractiveIndex("N-SK16", 587.56), 1.62041, 1e-5); close(optics.refractiveIndex("F2", 587.56), 1.62004, 1e-5);
  const lens = { id: "L", name: "L", kind: "lens", z: 0, r1: 60, r2: -45, thickness: 6, diameter: 25, glass: "N-BK7" }, d = { wavelength: 632.8, source: { power: 1, diameter: 4 }, elements: [lens], target: { name: "t", z: 200, diameter: 10 }, rays: 64, fields: [0] };
  close(optics.paraxial(d).efl, optics.lensmaker(lens, 632.8), 1e-10);
  const n = optics.refractiveIndex("N-BK7", 632.8), win = { id: "W", name: "W", kind: "window", z: 60, thickness: 10, diameter: 25, glass: "N-BK7" };
  close(optics.paraxial({ ...d, elements: [lens, win] }).focus - optics.paraxial(d).focus, 10 * (1 - 1 / n), 1e-9, "window focus shift");
  close(optics.fresnel(1, 1.5, 1, 1), 1 - .04, 1e-15, "normal-incidence Fresnel");
  const c = optics.validateOptics(optics.opticsExample()), r = optics.opticsResults(c, 60);
  close(r.paraxial.efl, 50, .01, "Cooke EFL"); close(optics.bestFocus(c), c.target.z, .05, "best focus at the prescription image plane");
  assert(r.trace.rms < .006); close(r.fields.find(f => f.field === 20).centroid[1] / (50 * Math.tan(20 * Math.PI / 180)), 1, .01, "20° image height");
  // Fresnel loss alone at six uncoated surfaces.
  const t = [1, 2, 3].map(() => 1).reduce(a => a, 1), on = r.trace.rays.find(ray => !ray.lost);
  assert(on.w < 1 && r.trace.transmitted > .69 && r.trace.transmitted < .72, `transmission ${r.trace.transmitted}`); void t;
});
test("fluids: Shah–London limits, Colebrook against the Moody value, balanced operating point", () => {
  close(fluids.shahLondon(0), 96, 1e-12); close(fluids.shahLondon(1), 56.91, .02);
  close(fluids.colebrook(1e5, 1e-4), .01852, 1e-4);
  const d = fluids.validateFluids({ ...fluids.fluidsExample(), fan: { preset: "axial60", shutoff: 30, free: .0085 } }), op = fluids.operatingPoint(d);
  const drop = op.elements.reduce((s, e) => s + (e.lo + e.hi) / 2, 0); close(fluids.fanPressure(d.fan, op.Q), drop, 1e-6 * drop + (op.high - op.low > 0 ? 1 : 1e-9));
  const t = fluids.friction(3000, .5, 0); assert.equal(t.regime, "transition"); assert(t.hi > t.lo);
});
test("architecture: a supplied room's heat balance closes; no gains floats to the supply", () => {
  const d = arch.validateArchitecture(arch.architectureExample()), r = arch.architectureResults(d, { equipment_power: 50 }), room = d.rooms[0], z = r.rooms[0];
  close(arch.AIR.rho * arch.AIR.cp * room.supply / 3600 * (z.temperature - room.supplyTemp), z.internal + z.solar + z.UA * (d.site.outdoor - z.temperature), 1e-9);
  const quiet = structuredClone(d); Object.assign(quiet.rooms[0], { people: 0, lighting: 0, equipment: false, roofU: 0 }); quiet.site.outdoor = quiet.rooms[0].supplyTemp; quiet.site.solar = 0;
  close(arch.architectureResults(quiet).rooms[0].temperature, quiet.rooms[0].supplyTemp, 1e-12);
});
test("flight: actuator-disc power and the mass budget", () => {
  const d = flight.validateFlight(flight.flightExample()), r = flight.flightResults(d, { esc_mass: .1, camera_mass: .01, parts_mass: .5 });
  close(r.mass, d.masses.reduce((s, m) => s + m.mass, 0) + d.battery.mass + .4 + .01 + .5, 1e-12);
  const T = r.weight / d.rotors; close(r.hover.ideal, T ** 1.5 / Math.sqrt(2 * r.rho * r.area), 1e-9);
  assert(r.design.esc > r.hover.esc && r.time > 0);
});
test("lithography: clear mask images to 1; a coherent grating past cutoff to its DC term squared; ILT improves the print", () => {
  const d = litho.validateLitho({ ...litho.lithoExample(), sigmaOuter: 0, sourcePoints: 1 }), clear = new Float64Array(d.n * d.n).fill(1);
  for (const v of litho.aerial(d, clear).I) close(v, 1, 1e-12);
  // Pitch 4 pixels × 20 nm = 80 nm < λ/NA = 331 nm: only DC passes.
  const grating = new Float64Array(d.n * d.n).map((_, k) => (k % d.n) % 4 < 2 ? 1 : 0);
  for (const v of litho.aerial(d, grating).I) close(v, .25, 1e-12);
  const p = litho.validateLitho(litho.lithoExample()), target = litho.rasterize(p), before = litho.patternError(p, target, litho.aerial(p, target).I);
  assert(litho.ilt(p, target, { iterations: 24 }).errors < before);
});
test("silicon: GDSII round-trips every block exactly", () => {
  const d = validateDocument({ format: "emergence-studio/1", kind: "silicon", design: siliconExample() }).design, bytes = gdsii(d), v = verifyGDS(d, bytes);
  assert(v.ok, v.problems.join(" ")); assert.equal(readGDS(bytes).shapes.length, d.blocks.length + 2);
  const moved = structuredClone(d); moved.blocks[1].x += .001; assert(!verifyGDS(moved, bytes).ok, "a 1 nm move is detected");
});
test("PCB: obstacle-aware routes pass the clearance checks; Gerber reads back as drawn; KiCad lists every footprint", async () => {
  const d = validateDocument({ format: "emergence-studio/1", kind: "pcb", design: pcbExample() }).design;
  d.routes.push(autoroute(d, "U1:4", "R1:1", "top", .4), autoroute(d, "J1:1", "U1:1", "top", .4), autoroute(d, "C1:1", "U1:8", "bottom", .4));
  assert.deepEqual(pcbChecks(d).issues.filter(i => i.code !== "pad-edge"), []);
  const k = kernelBoard(d), { result } = await models.request("pcb", { input: { board: k.board, currents: k.currents, fabrication: true } });
  const top = gerberCheck(d, result.gerber.find(l => l.name === "top").gerber, "top", .1);
  assert(top.agreement > .9, `Gerber agreement ${top.agreement}`);
  const pcb = kicadPCB(d); assert.equal((pcb.match(/\(footprint /g) || []).length, d.components.length); assert.equal((pcb.match(/\(segment /g) || []).length, d.routes.reduce((s, r) => s + r.points.length - 1, 0));
  d.currents.VCC = 30; assert(densityIssues(d).length > 0);
});
test("the drone example settles, refuses incompatible links and reports out-of-range drives", async () => {
  const p = await opened("drone"), s = await p.solve();
  assert(s.converged && s.cyclic, "the weight spiral is a coupled loop that settles");
  close(s.applied["mech.current"].value, s.outputs["power.fet_rms"].value, 1e-12);
  close(s.applied["flight.esc_mass"].value, s.outputs["mech.sink_mass"].value + s.outputs["pcb.board_mass"].value + s.outputs["power.capacitor_mass"].value, 1e-12);
  assert(!s.checks.some(c => c.status === "unmet" || c.status === "unavailable"), s.checks.filter(c => c.status !== "met").map(c => c.text).join("; "));
  await assert.rejects(p.bus.run("system.link", { from: "flight.total_mass", to: "mech.velocity" }), /cannot drive/);
  await p.bus.run("flight.set", { capture: .9 });
  // Either the drive leaves the mechanical model's validated range and the link
  // is refused, or the fins go turbulent and the model reports unavailable.
  const hot = await p.solve();
  assert(hot.issues.some(i => i.kind === "range" && i.link === "fluids-mech") || hot.checks.some(c => c.domain === "mech" && c.status === "unavailable"), "an out-of-range drive is refused, not extrapolated");
});
test("revisions carry every domain as content-addressed blobs; tampering is refused", async () => {
  const p = await project(), repo = await new Repository().init(); repo.blobKinds = Object.fromEntries(p.store.ids().map(id => [id, d => p.store.domains[id].validate(d)]));
  const studios = await repo.putStudios(p.store.all()), c = await repo.commit("main", [repo.root], { fins: 9, velocity: 2.5, current: 30, ambient: 45, heat: 0, opticalSize: .1 }, { author: "T", domain: "Systems" }, "All domains", [], {}, studios);
  const again = await repo.putStudios(p.store.all()); assert.deepEqual(again, studios);
  const sealed = await seal(repo.bundle()), check = await unseal(JSON.parse(JSON.stringify(sealed))); assert(check.ok);
  const copy = await new Repository().init(); copy.blobKinds = repo.blobKinds; await copy.ingest(check.body); assert.deepEqual(copy.studioDocs(copy.commits.get(c.id)).flight, p.store.get("flight"));
  const bad = structuredClone(check.body); bad.blobs.find(b => b.kind === "flight").doc.prop = 999;
  const victim = await new Repository().init(); victim.blobKinds = repo.blobKinds; await assert.rejects(victim.ingest(bad), /verification/);
  assert(!(await unseal({ ...sealed, commits: [] })).ok);
});
test("zip archives carry correct CRC-32 and sizes", () => {
  assert.equal(crc32(new TextEncoder().encode("123456789")), 0xcbf43926);
  const z = zip([{ name: "a.txt", data: "hello" }]), dv = new DataView(z.buffer); assert.equal(dv.getUint32(0, true), 0x04034b50); assert.equal(dv.getUint32(18, true), 5);
});
test("the MCP server lists tools and runs a command", async () => {
  const child = spawn(process.execPath, [new URL("./cad/agent/mcp.mjs", import.meta.url).pathname, "--example", "drone"]), lines = [];
  child.stdout.on("data", d => lines.push(...d.toString().trim().split("\n")));
  child.stdin.write(JSON.stringify({ jsonrpc: "2.0", id: 1, method: "initialize", params: {} }) + "\n" + JSON.stringify({ jsonrpc: "2.0", id: 2, method: "tools/list" }) + "\n"
    + JSON.stringify({ jsonrpc: "2.0", id: 3, method: "tools/call", params: { name: "run_command", arguments: { command: "optics.design_lens L1 60" } } }) + "\n");
  await new Promise((resolve, reject) => { const t = setTimeout(() => reject(new Error("timeout")), 60000); const i = setInterval(() => { if (lines.length >= 3) { clearInterval(i); clearTimeout(t); resolve(); } }, 50); });
  child.kill();
  const [init, list, call] = lines.map(l => JSON.parse(l));
  assert.equal(init.result.serverInfo.name, "emergence-cad"); assert(list.result.tools.some(t => t.name === "solve"));
  assert.match(JSON.parse(call.result.content[0].text).result.message, /R1 = /);
});

test("project files: every example's manifest, links and documents validate, and its ports exist", async () => {
  const list = await examples(); assert(list.length >= 2);
  for (const e of list) {
    const file = await exampleFile(e.id), m = validateProject(file.documents.project), p = await opened(e.id);
    assert.equal(m.name, e.name, `${e.id}: index name matches the file`);
    const ports = p.system.ports();
    for (const h of m.headline) assert(ports.outputs.has(h.port), `${e.id}: headline ${h.port}`);
    for (const l of p.store.docs.system.links) assert.equal(p.system.check(l, ports), null, `${e.id}: link ${l.id}`);
    const doc = d => d === "mech" ? file.documents.mech : p.store.docs[d];
    for (const x of [...m.studies, ...m.optimize?.variables ?? []]) assert(m.domains.includes(x.domain) && Number.isFinite(getPath(doc(x.domain), x.path)), `${e.id}: ${x.domain} ${x.path}`);
    for (const step of m.tour.filter(s => s.command)) assert.doesNotThrow(() => p.bus.parse(step.command), `${e.id}: tour command ${step.command}`);
  }
});
test("the edge shelter closes its loop through the building and every check holds", async () => {
  const p = await opened("edge-shelter"), s = await p.solve();
  assert(s.converged && s.cyclic, "room air → inlet → node heat → room air is a coupled loop that settles");
  assert.deepEqual(p.active(), ["mech", "parts", "pcb", "silicon", "fluids", "architecture"]);
  // Inside a loop the applied value is the last pass's: equal to the solver's tolerance.
  const settled = (a, b, msg) => close(a, b, 1e-4 * Math.abs(b), msg);
  settled(s.applied["mech.ambient"].value, s.outputs["architecture.room_temperature"].value);
  settled(s.applied["architecture.equipment_power"].value, s.outputs["silicon.power"].value + s.outputs["mech.q1_power"].value + s.outputs["pcb.dissipation"].value, "per-node heat sums its three links");
  assert(!s.checks.some(c => c.status === "unmet" || c.status === "unavailable"), s.checks.filter(c => c.status !== "met").map(c => c.text).join("; "));
  const t0 = s.outputs["architecture.room_temperature"].value;
  await p.bus.run("architecture.site", { units: 100 }); const more = await p.solve();
  assert(more.outputs["architecture.room_temperature"].value > t0 + 3 && more.outputs["silicon.junction"].value > s.outputs["silicon.junction"].value, "more nodes warm the shelter and every die in it");
});
test("a blank project grows by commands; dormant links, refused domains and multi-document undo", async () => {
  const p = await project(); p.load(await exampleFile("blank"));
  assert.deepEqual(p.active(), []); assert.equal((await p.solve()).checks.length, 0);
  await assert.rejects(p.bus.run("optics.focus"), /not in this project/);
  await p.bus.run("project.add_domain", { domain: "optics" }); await p.bus.run("project.add_domain", { domain: "fluids" });
  await p.bus.run("system.link", { from: "optics.power_on_target", to: "fluids.heat_load" });
  await p.bus.run("optics.set", { power: 2 });
  let s = await p.solve(); close(s.applied["fluids.heat_load"].value, s.outputs["optics.power_on_target"].value, 1e-12);
  await p.bus.run("project.remove_domain", { domain: "optics" }); s = await p.solve();
  assert.equal(s.dormant.length, 1); assert(!s.applied["fluids.heat_load"]); assert.equal(p.store.docs.system.links.length, 1, "the link is kept");
  await p.bus.run("project.add_domain", { domain: "optics" }); assert.equal((await p.solve()).dormant.length, 0);
  // A path edit is validated by its domain and undone as one step.
  await p.bus.run("project.set_value", { domain: "optics", path: "elements.L1.z", value: 12 }); assert.equal(getPath(p.store.docs.optics, "elements.L1.z"), 12);
  await assert.rejects(p.bus.run("project.set_value", { domain: "optics", path: "elements.L1.r1", value: 1 }), /radius/);
  p.bus.undo(); assert.equal(getPath(p.store.docs.optics, "elements.L1.z"), 10);
  await p.bus.run("project.pin", { port: "optics.efl" }); await p.bus.run("project.study", { domain: "optics", path: "target.z" });
  assert.equal(p.store.docs.project.headline[0].port, "optics.efl"); assert.equal(p.store.docs.project.studies[0].path, "target.z");
  assert(numericPaths(p.store.docs.optics).includes("elements.L2.thickness"));
  assert.throws(() => setPath(p.store.docs.optics, "elements.L1.glass", 1), /not a numeric field/);
  assert.throws(() => validateProject({ name: "x", domains: ["warp"] }), /one of/);
});
test("a saved project file reopens identically; headless files from before projects open with every domain", async () => {
  const p = await opened("drone"), path = join(tmpdir(), `emergence-${process.pid}.emergence.json`);
  await p.bus.run("flight.mass", { id: "gimbal", mass: 1.1 }); await writeFile(path, JSON.stringify(await seal(p.file())));
  const back = await project(); back.load(await readFileAt(path)); assert.deepEqual(back.store.all(), p.store.all());
  const legacy = join(tmpdir(), `emergence-legacy-${process.pid}.json`), docs = p.store.all(); delete docs.project;
  await writeFile(legacy, JSON.stringify({ format: "emergence-headless/1", docs, mech: p.file().documents.mech }));
  const old = await project(); old.load(await readFileAt(legacy)); assert.deepEqual(old.active(), DOMAINS);
});
test("review fixes: rolled-back scripts stay out of the log; one equipment room; merges see the starter document", async () => {
  const p = await project(); p.load(await exampleFile("blank"));
  await p.bus.run("project.add_domain", { domain: "architecture" });
  const before = p.bus.script().length;
  await assert.rejects(p.bus.batch([{ op: "architecture.site", args: { outdoor: 30 } }, { op: "architecture.site", args: { outdoor: 99 } }]));
  assert.equal(p.store.docs.architecture.site.outdoor, 32); assert.equal(p.bus.script().length, before, "a rolled-back step is not an exported edit");
  await p.bus.run("architecture.add_room", { label: "Lab", equipment: true });
  assert.deepEqual(p.store.docs.architecture.rooms.filter(r => r.equipment).map(r => r.name), ["Lab"]);
  await p.bus.run("architecture.add_room", { label: "Store", equipment: false });
  assert.deepEqual(p.store.docs.architecture.rooms.filter(r => r.equipment).map(r => r.name), ["Lab"], "a non-equipment room leaves the flag where it was");
  // Two branches from a revision without documents: a one-sided PCB edit is not a conflict.
  const repo = await new Repository().init(); repo.blobKinds = Object.fromEntries(p.store.ids().map(id => [id, d => p.store.domains[id].validate(d)])); repo.studioDefaults = p.store.all();
  const params = { fins: 9, velocity: 2.5, current: 30, ambient: 45, heat: 0, opticalSize: .1 }, who = { author: "T", domain: "Systems" };
  const edited = { ...p.store.all(), pcb: { ...p.store.all().pcb, width: 90 } };
  const a = await repo.commit("main", [repo.root], params, who, "Edit the board", [], {}, await repo.putStudios(edited));
  const b = await repo.commit("other", [repo.root], params, who, "Save defaults", [], {}, await repo.putStudios(p.store.all()));
  const m = repo.merge(a.id, b.id); assert.deepEqual(m.studioConflicts, []); assert.equal(m.studios.pcb.width, 90);
});
test("review fixes: chip corner libraries are preloaded again after the model Worker restarts", async () => {
  let preloads = 0;
  const request = (type, args) => { if (type === "chip" && args.input?.libs) preloads++; return models.request(type, args); };
  request.generation = 1;
  const p = createProject({ request, mech: new MechanicalState(), loadCorners, verilog: await loadVerilog() }); p.load(await exampleFile("blank"));
  await p.bus.run("project.add_domain", { domain: "silicon" }); await p.solve(); assert.equal(preloads, 1);
  await p.bus.run("silicon.chip", { frequency: 110 }); await p.solve(); assert.equal(preloads, 1, "the same Worker keeps its libraries");
  request.generation = 2; await p.bus.run("silicon.chip", { frequency: 120 }); const s = await p.solve();
  assert.equal(preloads, 2, "a restarted Worker gets them again"); assert(!s.checks.some(c => c.status === "unavailable"));
});
// A person starts with nothing. Replaying each example's build script on a
// blank project must give back every document and every number: a field no
// command can set would fail here, whatever a starter happened to contain.
test("from nothing: every example rebuilds exactly from a blank project, using only commands", async () => {
  for (const { id } of await examples()) {
    const shipped = await opened(id), file = shipped.file(), script = buildScript(file.documents), p = await project(); p.load(await exampleFile("blank"));
    await p.bus.batch(script, { source: "test", label: "Rebuild" });
    const got = p.file().documents, want = file.documents, used = ["project", "system", ...want.project.domains, ...want.project.domains.includes("silicon") ? ["litho"] : []];
    for (const d of used) assert.deepEqual(got[d], want[d], `${id}: ${d}`);
    const a = await p.solve(), b = await shipped.solve();
    assert.equal(a.checks.length, b.checks.length); assert.equal(a.checks.filter(c => c.status === "met").length, b.checks.filter(c => c.status === "met").length);
    for (const [port, o] of Object.entries(b.outputs)) if (Number.isFinite(o.value)) assert.equal(a.outputs[port].value, o.value, `${id}: ${port}`);
  }
});
test("from nothing: every domain can be empty, solves, and raises no false alarm", async () => {
  for (const domain of DOMAINS) {
    const p = await project(); p.load(await exampleFile("blank")); await p.bus.run("project.add_domain", { domain, start: "empty" });
    const s = await p.solve(), bad = s.checks.filter(c => c.status === "unmet" || c.status === "unavailable");
    assert.deepEqual(bad.map(c => c.text), [], domain);
    await p.bus.run("project.reset_domain", { domain, start: "starter" }); await p.solve();
  }
});
test("from nothing: every command has a visible control", async () => {
  const dir = new URL("./cad/", import.meta.url), names = ["studios.js", "workspace.js", "shell.js", ...(await readdir(new URL("ui/", dir))).map(f => `ui/${f}`)];
  const source = (await Promise.all(names.map(f => readFile(new URL(f, dir), "utf8")))).join("\n"), p = await project();
  assert.deepEqual(p.bus.list().map(c => c.op).filter(op => !source.includes(`"${op}"`)), []);
});
test("from nothing: a new design with its own parts, footprint, formulas and requirements", async () => {
  const p = await project(), run = (op, args) => p.bus.run(op, args); p.load(await exampleFile("blank"));
  await p.bus.batch([{ op: "project.rename", args: { name: "Handheld inspection lamp" } }, ...["parts", "pcb", "optics", "equations"].map(domain => ({ op: "project.add_domain", args: { domain, start: "empty" } }))]);
  // A body: a parametric block, bored, hollowed and patterned with vents.
  const { id } = await run("parts.add", { template: "disc", id: "body", name: "Body", material: "al-6061" });
  await run("parts.define", { id, name: "bore", value: "diameter/3" }); await run("parts.solid", { id, how: "cut", shape: "cylinder", radius: "bore/2", height: 40, z: 6 });
  await run("parts.solid", { id, how: "add", shape: "extrude", points: [[28, -6], [44, -6], [44, 6], [28, 6]], height: 8, z: 4 }); await run("parts.modify", { id, op: "polar", count: 3, axis: "z" });
  const before = (await p.solve()).outputs["parts.mass"].value; assert(before > 0);
  await run("parts.param", { id, name: "diameter", value: 80 }); assert((await p.solve()).outputs["parts.mass"].value > before, "the bore is an expression of the diameter");
  // A driver board on a footprint that is not built in.
  await run("pcb.footprint", { name: "soic14", kind: "dual", pins: 14, pitch: 1.27, span: 5.4, mass: .15, pad: .8, drill: .4 });
  await assert.rejects(run("pcb.footprint", { name: "bad", kind: "dual", pins: 8, pad: .5, drill: .6 }), /Invalid footprint/);
  await run("pcb.place", { footprint: "soic14", at: [20, 20], ref: "D1", name: "LED driver" }); await run("pcb.place", { footprint: "connector", at: [8, 20], ref: "J1" });
  await run("pcb.net", { pad: "D1:1", net: "VIN" }); await run("pcb.net", { pad: "J1:1", net: "VIN" }); await run("pcb.current", { net: "VIN", amps: .5 }); await run("pcb.route", { from: "J1:1", to: "D1:1", width: .6 });
  assert((await p.solve()).checks.find(c => c.domain === "pcb" && c.id === "geometry").reason.includes("has no net"), "an unfinished pad is still flagged");
  for (let i = 2; i <= 14; i++) await run("pcb.net", { pad: `D1:${i}`, net: "NC" }); await run("pcb.net", { pad: "J1:2", net: "NC" });
  await assert.rejects(run("pcb.connect", { from: "D1:2", to: "D1:3" }), /same assigned net/);
  // A lens, and a formula no domain knows, linked in and out.
  await run("optics.add", { kind: "lens", id: "front", z: 10, r1: 40, r2: -40, thickness: 6, diameter: 24 }); await run("optics.focus"); await run("optics.set", { power: 2 });
  await run("equations.input", { id: "power", value: 1, unit: "W" }); await run("equations.input", { id: "battery_wh", value: 7.4, unit: "" }); await run("equations.output", { id: "hours", expr: "battery_wh / (power + 0.4)", label: "Runtime", unit: "h" });
  await run("system.link", { from: "pcb.dissipation", to: "equations.power" });
  await run("project.require", { port: "equations.hours", min: 2, label: "Runs at least two hours" }); await assert.rejects(run("project.require", { port: "pcb.d1_power", max: 1 }), /No output/);
  const s = await p.solve(), req = s.checks.find(c => c.domain === "project" && c.text === "Runs at least two hours");
  assert.equal(s.checks.find(c => c.domain === "pcb" && c.id === "geometry").status, "met", "pad size and drill belong to the footprint");
  const gerber = (await p.bus.run("pcb.current", { net: "VIN", amps: .5 }), kernelBoard(p.store.docs.pcb)); assert.equal(gerber.board.vias.find(v => Math.abs(v.drill - .4e-3) < 1e-9)?.drill, .4e-3, "the drill reaches the kernel board");
  assert(!s.checks.some(c => c.status === "unavailable"), JSON.stringify(s.checks.filter(c => c.status === "unavailable")));
  assert.equal(req.status, "met"); close(s.outputs["equations.hours"].value, 7.4 / (s.outputs["pcb.dissipation"].value + .4), 1e-9);
  assert(p.system.ports().inputs.has("pcb.d1_power"), "every component has its own power port");
  await run("pcb.power", { id: "D1", watts: 30 }); assert.equal((await p.solve()).checks.find(c => c.text === "Runs at least two hours").status, "unmet");
  // It is a project file like any other: its own script rebuilds it.
  const again = await project(); again.load(await exampleFile("blank")); await again.bus.batch(buildScript(p.file().documents));
  for (const d of ["parts", "pcb", "optics", "equations", "system", "project"]) assert.deepEqual(again.file().documents[d], p.file().documents[d], d);
});
test("equations are parsed, never executed", () => {
  const scope = { a: 2, b: 3 }, value = text => evaluate(parse(text), scope);
  close(value("a + b * 2 ^ 2"), 14, 1e-12); close(value("-a ^ 2"), -4, 1e-12); close(value("max(a, b) / sqrt(4) + pi - pi"), 1.5, 1e-12); close(value("1.5e3 / (a + 1)"), 500, 1e-9);
  for (const bad of ["a +", "(a", "alert(1)", "a b", "constructor.prototype", "1 ; 2"]) assert.throws(() => value(bad), undefined, bad);
  assert.throws(() => validateEquations({ inputs: [], outputs: [{ id: "x", expr: "y + 1" }] }), /not an input/);
  assert.throws(() => validateEquations({ inputs: [{ id: "x", value: 1, unit: "furlongs" }], outputs: [] }), /unknown unit/);
});
test("a history saved before a document gained fields still opens, unchanged, and unknown keys are still refused", async () => {
  const p = await project(); await p.bus.run("project.add_domain", { domain: "pcb" });
  const kinds = Object.fromEntries(p.store.ids().map(id => [id, d => p.store.domains[id].validate(d)]));
  const saved = async doc => { const body = { kind: "pcb", doc }, blob = { ...body, id: await hash(body) }, repo = await new Repository().init(); repo.blobKinds = kinds; repo.studioDefaults = p.store.all();
    const c = await repo.make({ parents: [repo.root], branch: "main", params: { fins: 9, velocity: 2.5, current: 30, ambient: 45, heat: 0, opticalSize: .1 }, author: "T", domain: "Systems", message: "Saved earlier", time: "2026-10-30T12:00:00.000Z", studios: { pcb: blob.id } });
    await repo.ingest({ ...repo.bundle(), commits: [...repo.bundle().commits, c], blobs: [blob] }); return { repo, c, blob }; };
  const old = p.store.snapshot("pcb"); delete old.footprints; delete old.pours;
  const { repo, c, blob } = await saved(old);
  assert.deepEqual(repo.studioDocs(repo.commits.get(c.id)).pcb, p.store.get("pcb"), "reads with today's fields, so no phantom draft");
  assert.deepEqual(repo.bundle().blobs[0].doc, old, "stored as hashed"); assert.equal(await hash({ kind: "pcb", doc: repo.bundle().blobs[0].doc }), blob.id, "and still verifies");
  const copy = await new Repository().init(); copy.blobKinds = kinds; await copy.ingest(repo.bundle());
  await assert.rejects(saved({ ...p.store.snapshot("pcb"), injected: 1 }), /verification/);
});
const G = 9.80665, droneDocs = async () => (await exampleFile("drone")).documents;
test("propulsion: bench rows read back, density follows similarity, and the rotor stops where the table does", async () => {
  const P = { ...(await droneDocs()).flight.propulsion, blockage: 0 }, top = P.table.at(-1), a = P.table[4], first = P.table[0], thrust = g => g / 1000 * G;
  const row = flight.propulsionAt(P, thrust(a.grams), P.density); close(row.rpm, a.rpm, 1e-9); close(row.torque, a.torque, 1e-12); close(row.watts, a.watts, 1e-9);
  // At one rpm, thrust, torque and power all go with density.
  const thin = flight.propulsionAt(P, .9 * thrust(a.grams), .9 * P.density); close(thin.rpm, a.rpm, 1e-9); close(thin.torque, .9 * a.torque, 1e-12); close(thin.watts, .9 * a.watts, 1e-9);
  // Below the table: half the speed is a quarter of the thrust and an eighth of the power.
  const low = flight.propulsionAt(P, thrust(first.grams) / 4, P.density); close(low.rpm, first.rpm / 2, 1e-9); close(low.torque, first.torque / 4, 1e-12); close(low.watts, first.watts / 8, 1e-9);
  assert.equal(flight.propulsionAt(P, 1.01 * thrust(top.grams), P.density), null, "the rotor cannot do what was never measured");
  const blocked = flight.propulsionAt({ ...P, blockage: .05 }, .95 * thrust(a.grams), P.density); close(blocked.rpm, a.rpm, 1e-9);
  close(flight.maxThrust(P, top.volts, P.density), thrust(top.grams), 1e-9); close(flight.maxThrust(P, top.volts / 2, P.density), thrust(top.grams) / 4, 1e-9);
});
test("the published bench table is consistent with physics it was not derived from", async () => {
  const { flight: f } = await droneDocs(), P = f.propulsion, area = Math.PI * (f.prop / 2000) ** 2;
  for (const r of P.table) {
    const omega = r.rpm * 2 * Math.PI / 60, shaft = r.torque * omega, T = r.grams / 1000 * G, ideal = T ** 1.5 / Math.sqrt(2 * P.density * area);
    assert(shaft < r.watts, `${r.rpm} rpm: shaft power ${shaft} above the ${r.watts} W electrical input`);
    assert(shaft / r.watts > .7, `${r.rpm} rpm: motor efficiency ${shaft / r.watts}`);
    const fm = ideal / shaft; assert(fm > .6 && fm < .85, `${r.rpm} rpm: figure of merit ${fm} is outside momentum theory's range for a real rotor`);
    assert(r.rpm < r.volts * P.kv, `${r.rpm} rpm is above the no-load speed at ${r.volts} V`);
  }
});
test("battery: the discharge curve is monotonic and the endurance agrees with the pack's energy", async () => {
  close(flight.cellVoltage(1), 4.2, .01); for (let s = 0; s < 1; s += .05) assert(flight.cellVoltage(s + .05) >= flight.cellVoltage(s));
  const p = await opened("drone"), s = await p.solve(), r = s.results.flight.details, b = p.store.get("flight").battery, volts = Array.from({ length: 80 }, (_, i) => flight.cellVoltage(1 - b.usable * (i + .5) / 80)), mean = volts.reduce((a, v) => a + v, 0) / 80;
  const energy = b.capacity * b.usable * b.cells * mean; close(r.time / 60 * r.hover.P / energy, 1, .05, "minutes × hover watts against watt-hours");
});
test("power stage: block commutation's currents add up, and each loss scales as its physics says", () => {
  const d = validatePower(powerExample()), at = { phase: 10, bus: 44, duty: .6, tj: 25 }, s = stageLoss(d, at);
  close(s.R, 3.4e-3, 1e-12, "the datasheet's RDS(on) at 25 °C");
  close(s.hiRms ** 2 + s.loRms ** 2, 2 / 3 * 100, 1e-9, "a 120° block has RMS² of ⅔ I²");
  close(3 * (s.hiRms ** 2 + s.loRms ** 2) * s.R, s.conduction, 1e-12, "the six devices sum to two in series");
  close(s.gate, 2 * 58e-9 * 12 * 24000, 1e-12, "2·Qg·Vdrive·f");
  const fast = stageLoss({ ...d, pwm: 2 * d.pwm }, at); close(fast.conduction, s.conduction, 1e-12); close(fast.switching, 2 * s.switching, 1e-12); close(fast.gate, 2 * s.gate, 1e-12); close(fast.dead, 2 * s.dead, 1e-12);
  assert(stageLoss(d, { ...at, tj: 100 }).conduction > 1.4 * s.conduction, "RDS(on) rises with temperature");
  assert(stageLoss({ ...d, gateResistor: 20 }, at).switching > s.switching, "a slower gate spends longer switching");
  const b = busRipple(d, { phase: 10, duty: .5 }); close(b.rms, 5, 1e-12, "I/2 at half duty"); close(b.share, 5 / 4, 1e-12);
  close(b.volts, 10 * .25 / (24000 * 4 * 470e-6) + 10 * .033 / 4, 1e-12, "charge per cycle plus the ESR step"); close(b.loss, 25 * .033 / 4, 1e-12);
});
// An independent Euler–Bernoulli finite-element beam: Hermite cubics, consistent
// mass, clamped root, tip mass lumped. It shares nothing with the closed forms.
function beam({ EI, mu, L, tip, elements = 24 }) {
  const h = L / elements, N = 2 * elements, K = Array.from({ length: N }, () => new Float64Array(N)), M = Array.from({ length: N }, () => new Float64Array(N));
  const ke = [[12, 6 * h, -12, 6 * h], [6 * h, 4 * h * h, -6 * h, 2 * h * h], [-12, -6 * h, 12, -6 * h], [6 * h, 2 * h * h, -6 * h, 4 * h * h]].map(r => r.map(v => v * EI / h ** 3));
  const me = [[156, 22 * h, 54, -13 * h], [22 * h, 4 * h * h, 13 * h, -3 * h * h], [54, 13 * h, 156, -22 * h], [-13 * h, -3 * h * h, -22 * h, 4 * h * h]].map(r => r.map(v => v * mu * h / 420));
  for (let e = 0; e < elements; e++) for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) { const a = 2 * e + i - 2, b = 2 * e + j - 2; if (a >= 0 && b >= 0) { K[a][b] += ke[i][j]; M[a][b] += me[i][j]; } }
  M[N - 2][N - 2] += tip;
  const solve = rhs => { const A = K.map((r, i) => [...r, rhs[i]]); for (let c = 0; c < N; c++) { let m = c; for (let r = c; r < N; r++) if (Math.abs(A[r][c]) > Math.abs(A[m][c])) m = r; [A[c], A[m]] = [A[m], A[c]]; for (let r = c + 1; r < N; r++) { const f = A[r][c] / A[c][c]; for (let k = c; k <= N; k++) A[r][k] -= f * A[c][k]; } } const x = new Float64Array(N); for (let r = N - 1; r >= 0; r--) { let v = A[r][N]; for (let k = r + 1; k < N; k++) v -= A[r][k] * x[k]; x[r] = v / A[r][r]; } return x; };
  const unit = new Float64Array(N); unit[N - 2] = 1; const flex = solve(unit)[N - 2];
  let u = new Float64Array(N).fill(1), lambda = 0;
  for (let it = 0; it < 60; it++) { const v = solve(M.map(r => r.reduce((s, m, j) => s + m * u[j], 0))), top = Math.max(...v.map(Math.abs)); lambda = top / Math.max(...u.map(Math.abs)); u = v.map(x => x / top); }
  return { flex, frequency: Math.sqrt(1 / lambda) / (2 * Math.PI) };
}
test("structures: an independent finite-element beam agrees with the closed forms on the drone's own arm", async () => {
  const { structures: doc } = await droneDocs(), d = validateStructures(doc), r = armResults(d, { thrust: 40, tip_mass: .19 });
  const ro = d.od / 2000, ri = ro - d.wall / 1000; close(r.I, Math.PI / 64 * ((2 * ro) ** 4 - (2 * ri) ** 4), 1e-16, "the tube's second moment"); close(r.area, Math.PI * (ro * ro - ri * ri), 1e-12);
  const fe = beam({ EI: d.modulus * 1e9 * r.I, mu: r.area * d.density * 1000, L: d.length / 1000, tip: .19 });
  close(r.deflection / (40 * fe.flex * 1000), 1, 1e-6, "tip deflection, PL³/3EI");
  close(r.frequency / fe.frequency, 1, .01, "Rayleigh's first mode with 33/140 of the arm");
  // Too long, and a clamped arm walks into the blade-pass band.
  assert(armResults(d, { length: 310, tip_mass: .19, rpm_low: 4200, rpm_high: 6000 }).conflicts.length, "a long arm meets blade-pass");
  assert.equal(armResults(d, { tip_mass: .19, rpm_low: 4200, rpm_high: 6000 }).conflicts.length, 0, "the drone's arm clears it");
});
test("controls: hover trims, the rigid body conserves angular momentum, and margins match a loop with a known answer", async () => {
  const d = validateControls(controlsExample()), P = plantOf(d, { mass: 7.25, ixx: .14, iyy: .14, izz: .28, arm: 368, rotors: 4, hover_throttle: 60, max_thrust: 38, motor_tau: .05, torque_ratio: .0225 });
  close(P.n * P.Th, 7.25 * G, 1e-9, "hover thrust is the weight");
  close(P.rotors.reduce((s, r) => s + r.mr ** 2, 0), P.n / 2, 1e-9, "PX4's mix normalisation"); close(P.rotors.reduce((s, r) => s + r.mr, 0), 0, 1e-9, "roll adds no net thrust"); close(P.rotors.reduce((s, r) => s + r.my, 0), 0, 1e-9, "yaw adds no net thrust");
  const hover = simulate(d, P, () => ({ position: [0, 0, 0] }), { seconds: 3, noise: false }); assert(Math.hypot(...hover.final.pos) < .02, `drifted ${hover.final.pos}`);
  // Kick it, let the motors settle to equal thrust, and it is a free rigid body.
  const free = simulate(d, P, () => ({ roll: 0 }), { seconds: 4, noise: false, open: t => t < .2 ? [.08, .05, .2] : [0, 0, 0], record: s => s });
  const L = free.out.filter(s => s.t > 1.5).map(s => Math.hypot(P.ixx * s.w[0], P.iyy * s.w[1], P.izz * s.w[2])); assert(L[0] > .01, "it is spinning");
  close(Math.max(...L) / Math.min(...L), 1, .005, "|L| of a torque-free body");
  // An integrator with a first-order lag: crossover at K·G and 90° less the lag.
  const quiet = { ...d, loopRate: 1e6, params: { ...d.params, IMU_GYRO_CUTOFF: 1e7, IMU_DGYRO_CUTOFF: 1e7 } }, kp = .05, a = rateLoop(quiet, { ...P, tau: 1e-6 }, 0, { MC_ROLLRATE_P: kp, MC_ROLLRATE_I: 0, MC_ROLLRATE_D: 0, MC_ROLLRATE_K: 1 });
  close(a.crossover / (kp * a.authority), 1, .002); close(a.phaseMargin, 90, .5);
  const tau = .05, b = rateLoop(quiet, { ...P, tau }, 0, { MC_ROLLRATE_P: kp, MC_ROLLRATE_I: 0, MC_ROLLRATE_D: 0, MC_ROLLRATE_K: 1 }); close(b.phaseMargin, 90 - Math.atan(b.crossover * tau) * 180 / Math.PI, .5);
  // The tuner's promise holds with the integrator in the loop.
  const t = tune(d, P), tuned = { ...d, params: { ...d.params, ...Object.fromEntries(Object.entries(t).filter(([k]) => k !== "bandwidth")) } }, m = rateLoop(tuned, P);
  assert(m.phaseMargin >= 54 && (Number.isNaN(m.gainMargin) || m.gainMargin >= 8), `tuned to ${m.phaseMargin}° and ${m.gainMargin} dB`);
  assert(m.phaseMargin > rateLoop(d, P).phaseMargin + 8, "tuning buys real margin over PX4's stock gains on this airframe");
  const f = flightTests(tuned, P); assert(!f.failure.survives, "a quad does not survive losing a rotor, and the model says so");
});
test("testing: correlation predicts what it never saw, and a physical test passes only on a recorded measurement", async () => {
  const table = [1000, 2000, 3000].map(rpm => ({ rpm, grams: 1e-4 * rpm ** 2, torque: 1e-7 * rpm ** 2 })), exact = correlate({ rows: [1500, 2500].map(rpm => ({ rpm, grams: 1e-4 * rpm ** 2, torque: 1e-7 * rpm ** 2 })) }, table);
  close(exact.rms, 0, 1e-9); close(exact.torqueRms, 0, 1e-9); close(correlate({ rows: [{ rpm: 1000, grams: 50 }] }, table).rms, 100, 1e-9, "a wrong measurement shows as error");
  const docs = await droneDocs(), real = correlate(docs.testing.datasets[0], docs.flight.propulsion.table); assert(real.rms < 5 && real.worst < 5, `the KV170 fit predicts the KV340 data to ${real.rms}%`);
  const plan = { datasets: [], plan: [{ id: "tap", title: "Tap test", kind: "physical", min: 229, unit: "Hz" }, { id: "sim", title: "Mode", kind: "simulation", port: "structures.frequency", min: 229 }] };
  const before = evaluatePlan(plan, { outputs: { "structures.frequency": { value: 258, unit: "Hz" } } });
  assert.equal(before.items[0].status, "not run", "a predicted pass is not a test"); assert.equal(before.items[1].status, "passed"); assert.equal(before.summary.physical.passed, 0);
  plan.plan[0].result = { value: 200 }; assert.equal(evaluatePlan(plan, { outputs: {} }).items[0].status, "failed", "the measurement decides");
  const s = await (await opened("drone")).solve(), items = s.results.testing.details.items.filter(i => i.kind === "physical");
  assert(items.length >= 8 && items.every(i => i.status === "not run"), "nothing about this aircraft has been measured, and the file says so");
});
test("the ESC board is a real layout: datasheet pinouts, surface pads, pours, and exports that read back", async () => {
  const p = await opened("drone"), d = p.store.get("pcb"), geo = pcbChecks(d);
  assert.deepEqual(geo.issues, []); assert.equal(geo.unrouted, 0);
  const pinNets = Object.fromEntries(pads(d).map(q => [q.id, q.net]));
  // TI UCC27211 SOIC-8: VDD, HB, HO, HS, HI, LI, VSS, LO. Infineon TDSON-8: S S S G D.
  assert.deepEqual([1, 2, 3, 4, 5, 6, 7, 8].map(i => pinNets[`U1:${i}`]), ["V12", "BOOT_A", "GHD_A", "PH_A", "HI_A", "LI_A", "GND", "GLD_A"]);
  assert.deepEqual([1, 2, 3, 4, 5].map(i => pinNets[`Q2:${i}`]), ["GND", "GND", "GND", "GL_A", "PH_A"], "the low side's source is at ground and its drain on the phase");
  assert.deepEqual([1, 2, 3, 4, 5].map(i => pinNets[`Q1:${i}`]), ["PH_A", "PH_A", "PH_A", "GH_A", "VBAT"], "the high side's drain is on the bus and its source on the phase");
  const k = kernelBoard(d), { result } = await models.request("pcb", { input: { board: k.board, currents: k.currents, fabrication: true, temperature: 40 } });
  const through = pads(d).filter(q => !q.smd).length, hits = result.excellon.split("\n").filter(l => /^X/.test(l)).length; assert.equal(hits, through, "only through pads are drilled"); assert(through < pads(d).length / 2, "this is a surface-mount board");
  for (const layer of ["top", "bottom"]) { const start = performance.now(), c = gerberCheck(d, result.gerber.find(l => l.name === layer).gerber, layer); assert(c.agreement > .95, `${layer} Gerber agreement ${c.agreement}`); assert(performance.now() - start < 5000, "the read-back runs on the page's thread, so it stays quick on a poured board"); }
  const pcb = kicadPCB(d), vias = d.components.filter(c => c.footprint === "via").length; assert.equal((pcb.match(/\(zone /g) || []).length, d.pours.length); assert.equal((pcb.match(/ smd rect /g) || []).length, pads(d).filter(q => q.smd).length);
  assert.equal((pcb.match(/\(via /g) || []).length, vias, "a via is a via to KiCad"); assert.equal((pcb.match(/\(footprint /g) || []).length, d.components.length - vias);
  // Stray slivers with nothing on them do no harm; KiCad drops them. Every pad a pour serves shares one island.
  for (const q of d.pours) { const { under } = pourIslands(d, q), served = pads(d).filter(a => a.net === q.net && (!a.smd || q.layer === "top")).flatMap(under); assert(served.length && new Set(served).size === 1, `${q.id} serves its pads from ${new Set(served).size} islands`); }
  assert.equal(d.components.filter(c => /^VSA/.test(c.id)).length, 6, "six vias take the low-side source current to the plane");
  // STM32F051K8T6 LQFP-32 as AM32's F0 group B wires it: TIM1 on PA8–10 with PB1, PB0, PA7; comparators on PA0, PA4, PA5 against PA1.
  const mcu = { 4: "NRST", 6: "BEMF_A", 7: "BEMF_N", 9: "ISENSE", 10: "BEMF_B", 11: "BEMF_C", 12: "VSENSE", 13: "LI_C", 14: "LI_B", 15: "LI_A", 18: "HI_C", 19: "HI_B", 20: "HI_A", 23: "SWDIO", 24: "SWCLK", 27: "INPUT", 31: "GND", 16: "GND", 32: "GND", 1: "3V3", 5: "3V3", 17: "3V3" };
  for (const [pin, want] of Object.entries(mcu)) assert.equal(pinNets[`U4:${pin}`], want, `U4 pin ${pin}`);
  assert.deepEqual([1, 2, 3, 4, 5, 6, 7, 8, 9].map(i => pinNets[`U6:${i}`]), ["GND", "BAT_P", "BAT_P", "RON", "FB", "VCC", "BST", "SW", "GND"], "TI LM5017: RTN, VIN, UVLO, RON, FB, VCC, BST, SW and its exposed pad");
  assert.deepEqual([1, 2, 3, 4, 5, 6, 7, 8].map(i => pinNets[`U7:${i}`]), ["VBAT", "GND", "GND", "GND", "ISENSE", "3V3", "GND", "BAT_P"], "TI INA240: IN− on the load side of the shunt, IN+ on the battery side");
  assert.deepEqual([1, 2, 3, 4, 5].map(i => pinNets[`U5:${i}`]), ["V12", "GND", "V12", "BP", "3V3"], "TI LP2985: VIN, GND, ON, BP, VOUT");
  // Loss: a stub into a pour carries none of the net's current; the pour carries it, by hand.
  const bus = d.pours.find(q => q.id === "vbat"), { length, width } = pourSize(bus), R = 1.678e-8 * (1 + 3.93e-3 * 20) * length * 1e-3 / (width * 1e-3 * 2 * 34.8e-6), one = { ...d, pours: [bus], currents: { VBAT: 12 } };
  close(pourLoss(one), 144 * R, 1e-9); assert(R < .01, `the bus pour is ${R} Ω`);
  const s = await p.solve(); assert(s.outputs["pcb.copper_loss"].value < 2, `copper loss ${s.outputs["pcb.copper_loss"].value} W on a board whose phase carries ${s.outputs["flight.phase_current"].value} A`);
  const checks = async () => Object.fromEntries((await p.solve()).checks.filter(c => c.domain === "power").map(c => [c.id, c.status]));
  assert.equal((await checks())["capacitor-ripple"], "met"); await p.bus.run("power.capacitor", { count: 2 }); assert.equal((await checks())["capacitor-ripple"], "unmet", "two capacitors overrun their rating at this current"); await p.bus.run("power.capacitor", { count: 3 });
  await p.bus.run("pcb.board", { copper: 1 }); close((await p.solve()).outputs["pcb.copper_loss"].value / s.outputs["pcb.copper_loss"].value, 2, .25, "half the copper, about twice the loss");
});
test("PCB connectivity believes the copper: a cut pour is open, airwires go to the nearest pad, and parts cannot overlap", () => {
  const pad = (id, x, y, net) => ({ id, name: id, footprint: "via", x, y, rotation: 0, nets: [net] });
  const board = { width: 60, height: 30, grid: .25, rules: { clearance: .2, minWidth: .2, edge: .5, density: 50 }, stackup: { thickness: 1.6, copper: 2, epsr: 4.4 }, footprints: {}, currents: {}, powers: {}, loadNets: [],
    pours: [{ id: "bus", net: "VBAT", layer: "top", points: [[2, 2], [58, 2], [58, 28], [2, 28]] }], components: [pad("A", 8, 15, "VBAT"), pad("B", 52, 15, "VBAT"), pad("X", 30, 1.5, "GND"), pad("Y", 30, 28.5, "GND")], routes: [] };
  assert.equal(pcbChecks(board).unrouted, 1, "the ground pair is open, the bus is joined by its pour");
  assert.equal(pourIslands(board, board.pours[0]).islands, 1);
  // A ground trace right across the pour leaves two islands, and the bus is open.
  board.routes.push({ id: "cut", from: "X:1", to: "Y:1", layer: "top", width: .3, points: [[30, 1.5], [30, 28.5]] });
  const cut = pcbChecks(board); assert.equal(pourIslands(board, board.pours[0]).islands, 2);
  assert.equal(cut.unrouted, 1, "the bus now needs a connection"); assert(cut.issues.some(i => i.code === "pour-split"), "and the check says why");
  // A trace on the other layer does not cut it.
  board.routes[0].layer = "bottom"; assert.equal(pourIslands(board, board.pours[0]).islands, 1); assert(!pcbChecks(board).issues.some(i => i.code === "pour-split"));
  const open = { ...board, pours: [], routes: [], components: [pad("P", 10, 15, "N"), pad("Q", 12, 15, "N"), pad("R", 40, 15, "N"), pad("S", 42, 15, "N")] };
  open.routes.push({ id: "r1", from: "P:1", to: "Q:1", layer: "top", width: .3, points: [[10, 15], [12, 15]] }, { id: "r2", from: "R:1", to: "S:1", layer: "top", width: .3, points: [[40, 15], [42, 15]] });
  const [wire] = pcbChecks(open).nets[0].airwires; assert.deepEqual([wire[0].id, wire[1].id], ["Q:1", "R:1"], "the shortest link between the two groups, not the first pad");
  const crowded = { ...board, pours: [], routes: [], components: [{ id: "U1", name: "U1", footprint: "dip8", x: 20, y: 15, rotation: 0, nets: [] }, { id: "U2", name: "U2", footprint: "dip8", x: 21, y: 15, rotation: 0, nets: [] }] };
  assert(pcbChecks(crowded).issues.some(i => i.code === "courtyard"));
});
test("a recorded result goes stale when what it covers changes, and real software is never a physical test", async () => {
  const docs = { controls: { params: { MC_ROLLRATE_P: .07 } }, flight: { mass: 7 } }, covers = ["controls", "flight"];
  assert.equal(digest({ flight: { mass: 7 }, controls: { params: { MC_ROLLRATE_P: .07 } } }, covers), digest(docs, covers)); assert.equal(digest({ controls: { b: 1, a: 2 } }, ["controls"]), digest({ controls: { a: 2, b: 1 } }, ["controls"]), "key order is not a change");
  const plan = { datasets: [], traces: [], plan: [{ id: "fly", title: "Fly it", kind: "software", covers, min: 45, result: { value: 60, digest: digest(docs, covers) } }, { id: "hand", title: "On hardware", kind: "physical", min: 45 }] };
  const now = evaluatePlan(plan, { outputs: {}, docs }); assert.equal(now.items[0].status, "passed"); assert.equal(now.summary.software.passed, 1); assert.equal(now.summary.physical.passed, 0, "a software result does not touch the physical tests");
  const moved = evaluatePlan(plan, { outputs: {}, docs: { ...docs, flight: { mass: 8 } } }); assert.equal(moved.items[0].status, "stale", "a heavier aircraft is not the one that was flown"); assert.equal(moved.summary.software.passed, 0); assert.match(moved.items[0].detail, /Run it again/);
  // The shipped evidence is current, and each kind stays what it is.
  const p = await opened("drone"), s = await p.solve(), items = s.results.testing.details.items, soft = items.filter(i => i.kind === "software");
  assert(soft.length >= 10); assert(soft.filter(i => i.id !== "am32-build").every(i => i.status === "passed"), soft.filter(i => i.status !== "passed").map(i => `${i.id}: ${i.status}`).join(", "));
  assert.equal(soft.find(i => i.id === "am32-build").status, "not run", "the firmware has not been compiled, and the file says so");
  assert(items.filter(i => i.kind === "physical").every(i => i.status === "not run")); assert.equal(p.store.get("testing").traces.length, 2);
  await p.bus.run("flight.mass", { id: "gimbal", mass: 1.6 }); const after = (await p.solve()).results.testing.details.items.filter(i => i.kind === "software");
  assert.equal(after.find(i => i.id === "px4-phase").status, "stale", "a heavier camera invalidates the flight"); assert.equal(after.find(i => i.id === "kicad-drc").status, "passed", "but not the board");
});
test("the router goes round another net's pour when it can", () => {
  const pad = (id, x, y, net) => ({ id, name: id, footprint: "via", x, y, rotation: 0, nets: [net] });
  const board = { width: 60, height: 40, grid: .5, rules: { clearance: .2, minWidth: .2, edge: .5, density: 50 }, stackup: { thickness: 1.6, copper: 2, epsr: 4.4 }, footprints: {}, currents: {}, powers: {}, loadNets: [],
    pours: [{ id: "bus", net: "VBAT", layer: "top", points: [[20, 2], [40, 2], [40, 30], [20, 30]] }, { id: "plane", net: "GND", layer: "bottom", points: [[1, 1], [59, 1], [59, 39], [1, 39]] }], components: [pad("A", 10, 16, "SIG"), pad("B", 50, 16, "SIG")], routes: [] };
  const top = autoroute(board, "A:1", "B:1", "top", .3); assert(Math.max(...top.points.map(q => q[1])) >= 30, "it goes below the pour"); assert(top.points.every(q => !(q[0] > 20.5 && q[0] < 39.5 && q[1] > 2 && q[1] < 29.5)), "and never through it");
  board.pours[0].points = [[20, 0], [40, 0], [40, 40], [20, 40]]; const forced = autoroute(board, "A:1", "B:1", "top", .3); assert(forced.points.length >= 2, "where a pour spans the board it still crosses, and the connectivity check reports the cut");
});
test("the PX4 evidence agrees with the model it checks", async () => {
  const { testing } = await droneDocs(), item = id => testing.plan.find(t => t.id === id), r = id => item(id).result.value;
  assert(r("px4-crossover") < 10, "crossover"); assert(r("px4-roll") < 1, "roll step"); assert(r("px4-rate") < 10, "rate overshoot"); assert(r("px4-phase") >= 45);
  const rate = testing.traces.find(t => t.id === "px4-rate").rows, tail = rate.filter(q => q[0] > .1), rms = Math.sqrt(tail.reduce((s, q) => s + (q[1] - q[2]) ** 2, 0) / tail.length); assert(rms < 6, `averaged rate trace ${rms} °/s from the model`);
});
