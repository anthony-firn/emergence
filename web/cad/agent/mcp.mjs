#!/usr/bin/env node
// An MCP server (stdio, JSON-RPC 2.0) over the emergence CAD registry.
//   claude mcp add emergence -- node /path/to/web/cad/agent/mcp.mjs --project my.emergence.json [--example drone]
// Agents get the same commands, documents and coupled solve as the browser.
// The project file is created from --example (blank by default) if missing,
// and every edit is saved back to it.
import { createInterface } from "node:readline";
import { openProject, summary, exportArtifact, examples, exampleFile } from "./headless.mjs";
import { buildScript } from "../domains/rebuild.js";
const argv = process.argv.slice(2), opt = name => { const i = argv.indexOf(name); return i >= 0 ? argv[i + 1] : undefined; };
const projectPath = opt("--project"), ctx = await openProject(projectPath, { example: opt("--example") || "blank" }), { project, save } = ctx;
const DOMAINS = [...project.store.ids(), "mech"];
const TOOLS = [
  { name: "list_commands", description: "List the design commands, with argument types, units and ranges. Every edit a person can make in the browser studios is one of these.", inputSchema: { type: "object", properties: { domain: { type: "string", description: `Optional domain: ${DOMAINS.join(", ")}` } } } },
  { name: "run_command", description: "Run one command in command-line syntax, e.g. 'flight.mass gimbal 1.4' or 'pcb.route U1:4 R1:1 top 0.5'. Returns the command's message; the linked system is not re-solved until you call solve.", inputSchema: { type: "object", properties: { command: { type: "string" } }, required: ["command"] } },
  { name: "run_script", description: "Run a list of {op, args} commands atomically: if any step fails, every touched domain is restored.", inputSchema: { type: "object", properties: { steps: { type: "array", items: { type: "object", properties: { op: { type: "string" }, args: { type: "object" } }, required: ["op"] } } }, required: ["steps"] } },
  { name: "solve", description: "Solve the project's linked domains (any of flight, mechanical/electrothermal, parts, PCB, silicon chip flow, optics, fluids, architecture) and return outputs, every check with its status, link issues and convergence.", inputSchema: { type: "object", properties: {} } },
  { name: "project", description: "The open project's manifest (name, domains, headline outputs, studies, optimisation) and every input and output port of its domains, with units, for drawing links with system.link.", inputSchema: { type: "object", properties: {} } },
  { name: "build_script", description: "The open project as the ordered commands that build it from a blank project: a worked example of how to make any design, ready for run_script.", inputSchema: { type: "object", properties: {} } },
  { name: "list_examples", description: "Example project files. Domains are added with project.add_domain and linked with system.link; examples show complete projects.", inputSchema: { type: "object", properties: {} } },
  { name: "open_example", description: "Replace the open design with an example project (or blank). Saved to the project file.", inputSchema: { type: "object", properties: { id: { type: "string" } }, required: ["id"] } },
  { name: "get_document", description: "Return one domain's design document as JSON.", inputSchema: { type: "object", properties: { domain: { type: "string", enum: DOMAINS } }, required: ["domain"] } },
  { name: "assumptions", description: "List the nominal, illustrative and prescribed inputs every answer depends on, and whether the user has confirmed each.", inputSchema: { type: "object", properties: {} } },
  { name: "undo", description: "Undo the last edit.", inputSchema: { type: "object", properties: {} } },
  { name: "export", description: "Produce a design artifact: kicad (.kicad_pcb text), gerber (copper, profile, drill text), gds (GDSII as base64, verified by read-back) or project (all documents).", inputSchema: { type: "object", properties: { kind: { type: "string", enum: ["kicad", "gerber", "gds", "project"] } }, required: ["kind"] } },
];
async function call(name, a = {}) {
  if (name === "list_commands") return project.bus.list(a.domain).map(c => ({ ...c, usage: project.bus.usage(c.op) }));
  if (name === "run_command") { const { op, args } = project.bus.parse(a.command); const r = await project.bus.run(op, args, { source: "mcp" }); await save(); return { op, args, result: r }; }
  if (name === "run_script") { const r = await project.bus.batch(a.steps, { source: "mcp", label: "Agent script" }); await save(); return r; }
  if (name === "solve") return summary(await project.solve(), project);
  if (name === "project") { const ports = project.system.ports(), list = m => [...m.values()].map(p => ({ port: p.port, label: p.label, unit: p.unit, ...(p.min !== undefined && { min: p.min, max: p.max }) }));
    return { ...project.store.docs.project, tour: undefined, inputs: list(ports.inputs), outputs: list(ports.outputs), links: project.store.docs.system.links }; }
  if (name === "build_script") return buildScript(project.file().documents);
  if (name === "list_examples") return [{ id: "blank", name: "Blank project" }, ...await examples()];
  if (name === "open_example") { project.load(await exampleFile(a.id)); await save(); return { opened: project.store.docs.project.name, domains: project.active() }; }
  if (name === "get_document") return a.domain === "mech" ? ctx.mech.snapshot() : project.store.get(a.domain);
  if (name === "assumptions") return project.assumptions();
  if (name === "undo") { const ok = project.bus.undo(); await save(); return { undone: ok }; }
  if (name === "export") { const r = await exportArtifact(ctx, a.kind); return r.bytes ? { name: r.name, base64: Buffer.from(r.bytes).toString("base64"), verified: r.verified } : r; }
  throw new Error(`Unknown tool ${name}`);
}
const send = m => process.stdout.write(JSON.stringify(m) + "\n");
const rl = createInterface({ input: process.stdin });
// Requests are handled strictly in order, so a solve never overtakes an edit.
let queue = Promise.resolve();
rl.on("line", line => { queue = queue.then(() => handle(line)); });
async function handle(line) {
  let msg; try { msg = JSON.parse(line); } catch { return send({ jsonrpc: "2.0", id: null, error: { code: -32700, message: "Parse error" } }); }
  const { id, method, params } = msg;
  if (id === undefined) return; // notifications
  try {
    if (method === "initialize") return send({ jsonrpc: "2.0", id, result: { protocolVersion: params?.protocolVersion || "2025-06-18", capabilities: { tools: {} }, serverInfo: { name: "emergence-cad", version: "1.0.0" },
      instructions: "Cross-domain design and simulation. A project is a file: the domains it uses, one design document per domain, and typed links between their ports. Call project to see the manifest and ports; edit with run_command or run_script (project.add_domain, system.link, and each domain's commands); then call solve to see every linked check. Values carry units; links refuse mismatched dimensions; models report unavailable rather than extrapolate." } });
    if (method === "ping") return send({ jsonrpc: "2.0", id, result: {} });
    if (method === "tools/list") return send({ jsonrpc: "2.0", id, result: { tools: TOOLS } });
    if (method === "tools/call") {
      try { const r = await call(params.name, params.arguments); send({ jsonrpc: "2.0", id, result: { content: [{ type: "text", text: JSON.stringify(r, null, 1) }] } }); }
      catch (e) { send({ jsonrpc: "2.0", id, result: { content: [{ type: "text", text: e.message }], isError: true } }); }
      return;
    }
    send({ jsonrpc: "2.0", id, error: { code: -32601, message: `Method not found: ${method}` } });
  } catch (e) { send({ jsonrpc: "2.0", id, error: { code: -32603, message: e.message } }); }
}
