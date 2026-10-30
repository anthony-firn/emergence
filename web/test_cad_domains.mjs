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
import { Repository } from "./cad/repository.js";
import { seal, unseal, zip, crc32 } from "./cad/core/integrity.js";
import * as optics from "./cad/domains/optics.js";
import * as fluids from "./cad/domains/fluids.js";
import * as arch from "./cad/domains/architecture.js";
import * as flight from "./cad/domains/flight.js";
import * as litho from "./cad/domains/litho.js";
import { gdsii, verifyGDS, readGDS } from "./cad/domains/silicon.js";
import { autoroute, kicadPCB, kernelBoard, gerberCheck, densityIssues } from "./cad/domains/pcb.js";
import { pcbExample, siliconExample, pcbChecks, validateDocument } from "./cad/studio-model.js";
import { nodeModels, loadCorners, loadVerilog } from "./cad/agent/models.mjs";
import { createProject } from "./cad/domains/index.js";
import { MechanicalState } from "./cad/domains/mechanical.js";
import { exampleFile, examples, readFileAt } from "./cad/agent/headless.mjs";
import { validateProject, getPath, setPath, numericPaths, DOMAINS } from "./cad/core/project.js";
import { writeFile } from "node:fs/promises";
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
  close(s.applied["mech.current"].value, s.outputs["flight.esc_current"].value, 1e-12);
  close(s.applied["flight.esc_mass"].value, s.outputs["mech.sink_mass"].value + s.outputs["pcb.board_mass"].value, 1e-12);
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
