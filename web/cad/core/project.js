// A project is a file. Its manifest names the domains in use, the outputs to
// headline, the inputs to study, an optional optimisation, a guided tour and
// the 3D layout; the links and one design document per domain sit beside it.
// Nothing in the system knows what is being designed: a drone, a building
// full of servers and a single lens are all just project files.
export const PROJECT_FORMAT = "emergence-project/1";
export const DOMAINS = ["flight", "mech", "parts", "structures", "power", "pcb", "silicon", "optics", "fluids", "architecture", "controls", "testing", "equations"];
export const STUDIOS = ["system", "history", "render", "mechanical", ...DOMAINS.filter(d => d !== "mech")];
export const projectDefault = () => ({ name: "Untitled project", summary: "", domains: [], headline: [], requirements: [], studies: [], optimize: null, tour: [], scene: { placements: {} }, credits: [] });

const PORT = /^[a-z]+\.[a-zA-Z0-9_]+$/, PATH = /^[a-zA-Z0-9_-]+(\.[a-zA-Z0-9_-]+){0,7}$/, REF = /^[a-z]+:[a-zA-Z0-9_-]+(\.[a-zA-Z0-9_-]+){0,7}$/, ID = /^[a-zA-Z][a-zA-Z0-9_-]{0,31}$/;
const finite = (v, lo = -1e9, hi = 1e9) => Number.isFinite(v) && v >= lo && v <= hi;
const str = (v, max, what, optional = true) => {
  if (v === undefined || v === null || v === "") { if (optional) return undefined; throw new Error(`${what} is required.`); }
  if (typeof v !== "string" || v.length > max) throw new Error(`${what} must be text of up to ${max} characters.`);
  return v;
};
const list = (v, max, what) => { if (v === undefined) return []; if (!Array.isArray(v) || v.length > max) throw new Error(`${what} must be a list of up to ${max}.`); return v; };
const port = (v, what) => { if (typeof v !== "string" || !PORT.test(v)) throw new Error(`${what} must be a port such as flight.total_mass.`); return v; };
const path = (v, what) => { if (typeof v !== "string" || !PATH.test(v) || v.split(".").some(s => UNSAFE.has(s))) throw new Error(`${what} must be a field path such as battery.capacity.`); return v; };
const domain = (v, what) => { if (!DOMAINS.includes(v)) throw new Error(`${what} must be one of ${DOMAINS.join(", ")}.`); return v; };
const UNSAFE = new Set(["__proto__", "prototype", "constructor"]);

export function validateProject(d) {
  if (!d || typeof d !== "object" || Array.isArray(d)) throw new Error("Invalid project manifest.");
  const name = str(d.name, 80, "Project name", false).trim();
  if (!name) throw new Error("Project name is required.");
  const domains = list(d.domains, DOMAINS.length, "Domains").map(x => domain(x, "Each domain"));
  if (new Set(domains).size !== domains.length) throw new Error("Each domain appears once.");
  const headline = list(d.headline, 12, "Headline outputs").map(h => ({ port: port(h?.port, "A headline output"), ...(str(h.label, 60, "A headline label") && { label: h.label }) }));
  const requirements = list(d.requirements, 32, "Requirements").map(q => {
    if (!q || typeof q.id !== "string" || !ID.test(q.id)) throw new Error("Each requirement needs an id.");
    if (q.max === undefined && q.min === undefined || q.max !== undefined && !finite(q.max) || q.min !== undefined && !finite(q.min)) throw new Error(`Requirement ${q.id} needs a min or max.`);
    return { id: q.id, port: port(q.port, "A requirement"), ...(q.min !== undefined && { min: q.min }), ...(q.max !== undefined && { max: q.max }), ...(str(q.label, 80, "A requirement label") && { label: q.label }) };
  });
  if (new Set(requirements.map(q => q.id)).size !== requirements.length) throw new Error("Requirement ids must be unique.");
  const studies = list(d.studies, 16, "Study inputs").map(s => {
    if (!s || typeof s.id !== "string" || !ID.test(s.id)) throw new Error("Each study input needs an id.");
    const out = { id: s.id, label: str(s.label, 60, "A study label", false), domain: domain(s.domain, "A study domain"), path: path(s.path, "A study path") };
    for (const k of ["min", "max"]) if (s[k] !== undefined) { if (!finite(s[k])) throw new Error(`Study ${s.id}: ${k} must be a number.`); out[k] = s[k]; }
    return out;
  });
  if (new Set(studies.map(s => s.id)).size !== studies.length) throw new Error("Study ids must be unique.");
  let optimize = null;
  if (d.optimize) {
    const o = d.optimize;
    if (!o.objective || !["min", "max"].includes(o.objective.sense)) throw new Error("An optimisation needs an objective port and a sense (min or max).");
    const variables = list(o.variables, 6, "Optimisation variables").map(v => {
      if (!finite(v?.min) || !finite(v.max) || v.max <= v.min || !finite(v.step, 1e-9) || v.step > v.max - v.min) throw new Error("Each variable needs min < max and a step within the range.");
      return { domain: domain(v.domain, "A variable domain"), path: path(v.path, "A variable path"), label: str(v.label, 60, "A variable label", false), min: v.min, max: v.max, step: v.step, ...(v.integer && { integer: true }) };
    });
    const constraints = list(o.constraints, 8, "Constraints").map(c => {
      if (c?.max === undefined && c?.min === undefined || c.max !== undefined && !finite(c.max) || c.min !== undefined && !finite(c.min)) throw new Error("Each constraint needs a min or max.");
      return { port: port(c.port, "A constraint port"), ...(c.min !== undefined && { min: c.min }), ...(c.max !== undefined && { max: c.max }), ...(str(c.label, 80, "A constraint label") && { label: c.label }) };
    });
    if (o.trials !== undefined && (!Number.isInteger(o.trials) || o.trials < 4 || o.trials > 200)) throw new Error("Optimisation trials must be 4–200.");
    optimize = { title: str(o.title, 80, "Optimisation title", false), ...(str(o.summary, 400, "Optimisation summary") && { summary: o.summary }),
      objective: { port: port(o.objective.port, "The objective"), sense: o.objective.sense }, constraints, variables, trials: o.trials ?? 40 };
  }
  const tour = list(d.tour, 16, "Tour").map(t => {
    if (!t || !STUDIOS.includes(t.studio)) throw new Error(`Each tour step names a studio: ${STUDIOS.join(", ")}.`);
    return { studio: t.studio, title: str(t.title, 80, "A tour title", false), text: str(t.text, 600, "Tour text", false),
      ...(str(t.command, 400, "A tour command") && { command: t.command }),
      ...(t.watch?.length && { watch: list(t.watch, 6, "Watched outputs").map(p => port(p, "A watched output")) }),
      ...(t.checks?.length && { checks: list(t.checks, 4, "Watched checks").map(c => { if (typeof c !== "string" || !/^[a-z]+\/[a-zA-Z0-9_-]{1,48}$/.test(c)) throw new Error("A watched check is domain/id."); return c; }) }) };
  });
  const placements = {};
  for (const [k, p] of Object.entries(d.scene?.placements || {})) {
    domain(k, "A placement");
    if (!p || typeof p !== "object") throw new Error(`Placement ${k} must be an object.`);
    const at = p.at ?? [0, 0, 0];
    if (!Array.isArray(at) || at.length !== 3 || !at.every(v => finite(v, -1e6, 1e6))) throw new Error(`Placement ${k}: at is [x, y, z] in mm.`);
    for (const a of ["rx", "rz", "phase"]) if (p[a] !== undefined && !finite(p[a], -360, 360)) throw new Error(`Placement ${k}: ${a} is degrees.`);
    if (p.count !== undefined && !(Number.isInteger(p.count) && p.count >= 1 && p.count <= 64 || typeof p.count === "string" && REF.test(p.count))) throw new Error(`Placement ${k}: count is 1–64 or a reference such as flight:rotors.`);
    placements[k] = { at, ...Object.fromEntries(["rx", "rz", "phase", "count"].filter(a => p[a] !== undefined).map(a => [a, p[a]])) };
  }
  const site = d.scene?.site;
  if (site && (typeof site !== "object" || site.room !== undefined && (typeof site.room !== "string" || !ID.test(site.room)) || site.z !== undefined && !finite(site.z, -1e5, 1e5))) throw new Error("The scene site is { room, z }.");
  const credits = list(d.credits, 32, "Credits").map(c => {
    if (!c || typeof c.url !== "string" || !/^https:\/\/[^\s]{4,300}$/.test(c.url)) throw new Error("Each credit links an https URL.");
    return { title: str(c.title, 120, "A credit title", false), url: c.url, ...(c.domain && { domain: domain(c.domain, "A credit domain") }),
      ...(str(c.licence, 40, "A licence") && { licence: c.licence }), ...(str(c.note, 300, "A credit note") && { note: c.note }) };
  });
  return { name, summary: str(d.summary, 600, "Summary") ?? "", domains, headline, requirements, studies, optimize, tour,
    scene: { placements, ...(site && { site: { ...(site.room && { room: site.room }), ...(site.z !== undefined && { z: site.z }) } }) }, credits };
}

// Field paths: object keys, or an array element by its id (or index).
// "masses.gimbal.mass" is the mass of the budget item whose id is gimbal.
function step(node, key) {
  if (UNSAFE.has(key) || node === null || typeof node !== "object") return undefined;
  if (Array.isArray(node)) return node.find(x => x && x.id === key) ?? (/^\d+$/.test(key) ? node[Number(key)] : undefined);
  return Object.hasOwn(node, key) ? node[key] : undefined;
}
export function getPath(doc, p) { return p.split(".").reduce((n, k) => step(n, k), doc); }
export function setPath(doc, p, value) {
  const keys = p.split("."), last = keys.pop(), parent = keys.reduce((n, k) => step(n, k), doc);
  const current = step(parent, last);
  if (!Number.isFinite(current)) throw new Error(`${p} is not a numeric field of this document.`);
  if (Array.isArray(parent)) { const i = parent.findIndex(x => x && x.id === last); parent[i >= 0 ? i : Number(last)] = value; }
  else parent[last] = value;
}
// Every numeric field, for choosing study inputs and optimisation variables.
export function numericPaths(doc, limit = 400) {
  const out = [];
  const walk = (node, prefix, depth) => {
    if (out.length >= limit || depth > 6) return;
    if (Number.isFinite(node)) { out.push(prefix); return; }
    if (Array.isArray(node)) node.forEach((x, i) => walk(x, `${prefix}.${x && typeof x.id === "string" ? x.id : i}`, depth + 1));
    else if (node && typeof node === "object") for (const [k, v] of Object.entries(node)) if (k !== "model" && k !== "shape") walk(v, prefix ? `${prefix}.${k}` : k, depth + 1);
  };
  walk(doc, "", 0);
  return out.map(p => p.replace(/^\./, ""));
}
// "flight:rotors" reads a number from another domain's document.
export function reference(docs, value) {
  if (typeof value !== "string") return value;
  const [d, p] = value.split(":"), v = getPath(docs[d], p);
  return Number.isFinite(v) ? v : undefined;
}

// Reading and writing project files. Documents are validated by their own
// domains when loaded; this checks the envelope.
export function readProjectFile(file) {
  if (!file || file.format !== PROJECT_FORMAT || !file.documents || typeof file.documents !== "object" || Array.isArray(file.documents)) throw new Error(`This is not an ${PROJECT_FORMAT} file.`);
  if (!file.documents.project) throw new Error("A project file needs a project manifest.");
  return { documents: file.documents };
}
export const projectFile = documents => ({ format: PROJECT_FORMAT, documents });
export const slug = name => name.toLowerCase().normalize("NFKD").replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "").slice(0, 48) || "project";
