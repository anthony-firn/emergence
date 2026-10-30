// Projects on this device. Each project is its own revision history, started
// from a file: an example, the blank template or one you open. A project keeps
// its id in every export and on every peer channel, so imports and peers join
// the right history. Switching project reloads the page on that project.
import { projectFile, projectDefault, DOMAINS, slug } from "./core/project.js";
const REGISTRY = "emergence-cad-projects-v1", SEED = "emergence-cad-seed-v1:", INCOMING = "emergence-cad-incoming-v1:", STORE = "emergence-cad-project-v1";
const read = (store, key) => { try { return JSON.parse(store.getItem(key) || "null"); } catch { return null; } };
const write = (key, value) => { try { localStorage.setItem(key, JSON.stringify(value)); return true; } catch { return false; } };
const fresh = () => `p-${crypto.randomUUID().slice(0, 12)}`;
export const EXAMPLES = await fetch(new URL("./projects/index.json", import.meta.url)).then(r => r.ok ? r.json() : { examples: [] }).then(j => j.examples).catch(() => []);
export const DEFAULT_EXAMPLE = "drone";
export async function exampleFile(id) {
  const e = EXAMPLES.find(x => x.id === id); if (!e) throw new Error(`No example ${id}.`);
  const r = await fetch(new URL(`./projects/${e.file}`, import.meta.url)); if (!r.ok) throw new Error(`${e.file}: ${r.status}`);
  return r.json();
}
export const blankFile = () => projectFile({ project: projectDefault() });

function registry() {
  let r = read(localStorage, REGISTRY);
  if (r && Array.isArray(r.projects)) return r;
  r = { projects: [], last: null };
  // Before projects, this browser held a single history under the unsuffixed
  // keys. It becomes the project "local", which keeps using those keys.
  if (localStorage.getItem(STORE)) { r.projects.push({ id: "local", name: "My project", created: new Date().toISOString(), source: { legacy: true } }); r.last = "local"; }
  write(REGISTRY, r);
  return r;
}
const reg = (() => { try { return registry(); } catch { return { projects: [], last: null }; } })();
const save = () => write(REGISTRY, reg);
export const projects = () => [...reg.projects].sort((a, b) => (b.opened || b.created).localeCompare(a.opened || a.created));
const find = id => reg.projects.find(p => p.id === id);
function add(meta) { reg.projects.push({ created: new Date().toISOString(), ...meta }); save(); return meta.id; }

// Which project this page is: ?project=<id>, ?example=<id>, ?new=blank, or the
// last one opened. A first visit opens the default example.
const query = new URLSearchParams(location.search);
async function resolve() {
  let id = query.get("project"), example = query.get("example");
  if (query.get("new") === "blank") id = add({ id: fresh(), name: "Untitled project", source: { blank: true } });
  else if (example && EXAMPLES.some(e => e.id === example)) { id = `example-${example}`; if (!find(id)) add({ id, name: EXAMPLES.find(e => e.id === example).name, source: { example } }); }
  if (!id || !find(id)) id = find(reg.last) ? reg.last : null;
  if (!id) { id = `example-${DEFAULT_EXAMPLE}`; if (!find(id)) add({ id, name: EXAMPLES.find(e => e.id === DEFAULT_EXAMPLE)?.name || "Example", source: { example: DEFAULT_EXAMPLE } }); }
  const meta = find(id);
  meta.opened = new Date().toISOString(); reg.last = id; save();
  return meta;
}
const meta = await resolve();
history.replaceState(null, "", `?project=${encodeURIComponent(meta.id)}${location.hash}`);
// Storage keys are per project; the migrated project keeps the original keys.
const keyFor = (id, base) => id === "local" ? base : `${base}:${id}`;
const key = base => keyFor(meta.id, base);
// The file a new project starts from, used once to make its first revision.
async function seed() {
  if (meta.source?.example) return exampleFile(meta.source.example);
  if (meta.source?.blank) return blankFile();
  const stored = read(localStorage, SEED + meta.id); if (stored) return stored;
  return null;
}
export const current = {
  id: meta.id, meta, key, channel: `emergence-cad-${meta.id}-v1`,
  seed, clearSeed: () => { try { localStorage.removeItem(SEED + meta.id); } catch { /* best effort */ } },
  takeIncoming: () => { const b = read(localStorage, INCOMING + meta.id); try { localStorage.removeItem(INCOMING + meta.id); } catch { /* best effort */ } return b; },
  // The documents a revision without a manifest opens with: work from before
  // projects existed used every domain.
  legacyManifest: () => ({ ...projectDefault(), name: meta.name, domains: [...DOMAINS] }),
};
export function rename(name) { const m = find(meta.id); if (m && m.name !== name) { m.name = name; save(); } }
export const open = id => { location.href = `?project=${encodeURIComponent(id)}${location.hash}`; };
export const openExample = id => { location.href = `?example=${encodeURIComponent(id)}${location.hash}`; };
export const newBlank = () => { location.href = "?new=blank#system"; };
// A project file becomes a new project; a revision bundle keeps its id, so a
// collaborator's export joins the same history on this device.
export function createFromFile(file) {
  const id = fresh(); if (!write(SEED + id, file)) throw new Error("Browser storage is unavailable, so a new project cannot be created here.");
  add({ id, name: file.documents.project?.name || "Imported project", source: { file: true } }); return id;
}
export function createFromBundle(bundle, name) {
  const id = bundle.id && /^[a-zA-Z0-9-]{1,64}$/.test(bundle.id) ? bundle.id : fresh();
  // Already here: the project merges the history when it opens.
  if (find(id)) { if (!write(INCOMING + id, bundle)) throw new Error("Browser storage is unavailable, so the history cannot be handed to the project."); return { id, exists: true }; }
  if (!write(keyFor(id, STORE), bundle)) throw new Error("Browser storage is unavailable, so a new project cannot be created here.");
  add({ id, name: name || bundle.name || "Imported project", source: { bundle: true } }); return { id, exists: false };
}
export function remove(id) {
  if (id === meta.id) throw new Error("Open another project before deleting this one.");
  reg.projects = reg.projects.filter(p => p.id !== id); save();
  try { for (const b of [STORE, "emergence-cad-autosave-v1"]) localStorage.removeItem(keyFor(id, b)); localStorage.removeItem(SEED + id); localStorage.removeItem(INCOMING + id); } catch { /* best effort */ }
}
export { slug };
