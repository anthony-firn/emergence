#!/usr/bin/env node
// emergence CAD from a terminal, on project files. The same commands as the browser:
//   node cli.mjs examples
//   node cli.mjs new my.emergence.json [--example drone|edge-shelter|blank]
//   node cli.mjs commands [domain]
//   node cli.mjs do "project.add_domain optics" --project my.emergence.json
//   node cli.mjs run script.json --project my.emergence.json      (atomic; - for stdin)
//   node cli.mjs solve --project my.emergence.json [--json]       (or --example drone)
//   node cli.mjs export kicad|gerber|gds|project <out> --project my.emergence.json
// --project saves every edit back to the file; --example works on a copy in memory.
import { readFile, writeFile, mkdir } from "node:fs/promises";
import { existsSync } from "node:fs";
import { join } from "node:path";
import { openProject, summary, exportArtifact, examples } from "./headless.mjs";
const argv = process.argv.slice(2), flag = name => { const i = argv.indexOf(name); if (i < 0) return undefined; const v = argv[i + 1]; argv.splice(i, 2); return v; }, has = name => { const i = argv.indexOf(name); if (i < 0) return false; argv.splice(i, 1); return true; };
const projectPath = flag("--project"), example = flag("--example"), asJSON = has("--json"), [verb, ...rest] = argv;
const print = s => { console.log(`${s.project}: ${s.converged ? "settled" : "NOT settled"} in ${s.iterations} passes (${s.ms} ms) · ${s.domains.join(", ") || "no domains"}`); for (const c of s.checks) console.log(`  ${c.status.padEnd(11)} ${c.domain.padEnd(12)} ${c.text}${c.reason ? " · " + c.reason.slice(0, 110) : ""}`); for (const i of s.issues) console.log(`  note: ${i}`); };
try {
  if (!verb || verb === "help") { console.log("Usage: examples | new <file> [--example id] | commands [domain] | do \"<command>\" | run <script.json|-> | solve [--json] | export <kind> <out>\n  --project <file> opens and saves a project file; --example <id> opens an example in memory."); process.exit(0); }
  if (verb === "examples") { for (const e of [{ id: "blank", name: "Blank project", summary: "No domains; add them with project.add_domain." }, ...await examples()]) console.log(`${e.id.padEnd(14)} ${e.name} · ${e.summary}`); process.exit(0); }
  if (verb === "new") {
    const out = rest[0]; if (!out) throw new Error("new <file> [--example id]"); if (existsSync(out)) throw new Error(`${out} exists; choose a new file name.`);
    const { save, project } = await openProject(undefined, { example: example || "blank" }); await save(out);
    console.log(`Wrote ${out}: ${project.store.docs.project.name} (${project.active().join(", ") || "no domains"}).`); process.exit(0);
  }
  if (["do", "run", "solve", "export"].includes(verb) && !projectPath && !example) throw new Error(`${verb} needs --project <file> or --example <id>. List examples with: cli.mjs examples`);
  if (projectPath && !existsSync(projectPath)) throw new Error(`${projectPath} does not exist. Create it with: cli.mjs new ${projectPath} --example <id>`);
  const ctx = await openProject(projectPath, { example }), { project, save } = ctx;
  if (verb === "commands") { for (const c of project.bus.list(rest[0])) console.log(`${project.bus.usage(c.op).padEnd(64)} ${c.summary || c.title}`); process.exit(0); }
  if (verb === "do") { const { op, args } = project.bus.parse(rest.join(" ")); const r = await project.bus.run(op, args, { source: "cli" }); console.log(`✓ ${op}${r?.message ? " · " + r.message : ""}`); await save(); }
  if (verb === "run") { const text = rest[0] === "-" || !rest[0] ? await new Promise(res => { let d = ""; process.stdin.on("data", c => d += c).on("end", () => res(d)); }) : await readFile(rest[0], "utf8");
    const r = await project.bus.batch(JSON.parse(text), { source: "cli", label: "Script" }); for (const x of r) console.log(`✓ ${x.op}${x.result?.message ? " · " + x.result.message : ""}`); await save(); }
  if (["do", "run", "solve"].includes(verb)) { const s = summary(await project.solve(), project); asJSON ? console.log(JSON.stringify(s, null, 1)) : print(s); }
  if (verb === "export") { const [kind, out] = rest; if (!out) throw new Error("export <kind> <output path>");
    const a = await exportArtifact(ctx, kind);
    if (a.files) { await mkdir(out, { recursive: true }); for (const f of a.files) await writeFile(join(out, f.name), f.text); console.log(`Wrote ${a.files.length} files to ${out} (not verified against a CAM tool)`); }
    else { await writeFile(out, a.bytes || a.text); console.log(`Wrote ${out}${a.verified ? ` (GDSII read back: ${a.verified.boundaries} boundaries)` : ""}`); } }
} catch (e) { console.error(`✗ ${e.message}`); process.exit(1); }
