import { merge3 } from "./core/merge.js";
// Immutable, content-addressed revisions. Branch tips are derived from ancestry,
// so replicas never resolve concurrent revisions by silently overwriting one.
export const FORMAT = "emergence-cad/1";
export const PROJECT = "motor-driver";
export const FIELDS = {
  fins: { label: "Heat-sink fins", unit: "", min: 5, max: 14, step: 1, domain: "Mechanical", value: 9 },
  velocity: { label: "Duct air velocity", unit: "m/s", min: 0.5, max: 4, step: 0.1, domain: "Mechanical", value: 2.5 },
  current: { label: "MOSFET RMS current", unit: "A", min: 5, max: 40, step: 1, domain: "Electrical", value: 30 },
  ambient: { label: "Inlet air temperature", unit: "°C", min: 15, max: 60, step: 1, domain: "Thermal", value: 45 },
  heat: { label: "Extra heat deposited at Q1", unit: "W", min: 0, max: 10, step: 0.5, domain: "Radiation", value: 0 },
  opticalSize: { label: "Optical size parameter x", unit: "", min: 0.01, max: 12, step: 0.01, domain: "Photonics", value: 0.1 },
  width: { label: "Sink width", unit: "mm", min: 30, max: 70, step: 1, domain: "Mechanical", value: 40 },
  length: { label: "Sink length", unit: "mm", min: 30, max: 90, step: 1, domain: "Mechanical", value: 50 },
  finHeight: { label: "Fin height", unit: "mm", min: 10, max: 35, step: 1, domain: "Mechanical", value: 20 },
  base: { label: "Base thickness", unit: "mm", min: 2, max: 8, step: 0.5, domain: "Mechanical", value: 3 },
  finThickness: { label: "Fin thickness", unit: "mm", min: 0.8, max: 2, step: 0.1, domain: "Mechanical", value: 1.2 },
  traceWidth: { label: "Power trace width", unit: "mm", min: 4, max: 15, step: 0.5, domain: "Electrical", value: 10 },
  opticalPower: { label: "Incident optical power at Q1", unit: "W", min: 0, max: 5, step: 0.1, domain: "Photonics", value: 0 },
  absorption: { label: "Prescribed absorbed fraction", unit: "%", min: 0, max: 100, step: 1, domain: "Photonics", value: 50 },
};
export const DEFAULTS = Object.fromEntries(Object.entries(FIELDS).map(([key, f]) => [key, f.value]));
// Keep the first prototype's deterministic root and content hashes intact.
const LEGACY_KEYS = ["fins", "velocity", "current", "ambient", "heat", "opticalSize"];
const LEGACY_DEFAULTS = Object.fromEntries(LEGACY_KEYS.map(k => [k, DEFAULTS[k]]));
export const values = params => ({ ...DEFAULTS, ...params });
export const featureValues = commit => commit.features || [];
export const drawingValues = commit => commit.drawings || {};
export const FEATURE_DOMAINS = ["Mechanical", "Electrical", "Thermal", "Photonics", "Radiation"];
const safeID=id=>typeof id==="string"&&/^[a-zA-Z0-9_-]{1,48}$/.test(id)&&!["__proto__","constructor","prototype"].includes(id);
export function validFeatures(features) {
  return Array.isArray(features) && features.length <= 16 && new Set(features.map(f => f?.id)).size === features.length && features.every(f =>
    f && safeID(f.id) && !["sink", "board", "q1", "optic"].includes(f.id)
    && validText(f.name, 60) && ["box", "cylinder"].includes(f.shape) && FEATURE_DOMAINS.includes(f.domain)
    && ["x", "y", "z"].every(k => Number.isFinite(f[k]) && Math.abs(f[k]) <= 100)
    && ["width", "depth", "height"].every(k => Number.isFinite(f[k]) && f[k] >= 1 && f[k] <= 50)
    && (f.shape !== "cylinder" || f.width === f.depth) && Object.keys(f).length === 10);
}
export function validDrawings(ds) {
  return ds && !Array.isArray(ds) && typeof ds === "object" && Object.keys(ds).length <= 17 && Object.entries(ds).every(([id,d]) =>
    safeID(id) && d && Object.keys(d).length === 11
    && ["number","title","material","finish","revision","prepared","checked"].every(k => typeof d[k] === "string" && d[k].length <= 80)
    && typeof d.notes === "string" && d.notes.length <= 600 && typeof d.tolerance === "string" && (d.tolerance === "" || /^\d+(\.\d{1,4})?$/.test(d.tolerance) && Number(d.tolerance)>0 && Number(d.tolerance)<=10)
    && ["third","first"].includes(d.projection) && ["Draft","For review","Released"].includes(d.status)
    && (d.status !== "Released" || ["number","title","material","finish","revision","prepared","checked","tolerance"].every(k=>d[k].trim().length>0) && !/not specified|unspecified/i.test(d.material+" "+d.finish)));
}
export const canonical = value => JSON.stringify(sort(value));
function sort(value) {
  if (Array.isArray(value)) return value.map(sort);
  if (value && typeof value === "object") return Object.fromEntries(Object.keys(value).sort().map(k => [k, sort(value[k])]));
  return value;
}
export async function hash(value) {
  return [...new Uint8Array(await crypto.subtle.digest("SHA-256", new TextEncoder().encode(canonical(value))))]
    .map(n => n.toString(16).padStart(2, "0")).join("");
}
const validName = v => typeof v === "string" && /^[a-zA-Z0-9][a-zA-Z0-9/_.-]{0,47}$/.test(v);
const validText = (v, max) => typeof v === "string" && v.trim().length > 0 && v.length <= max;
const isHash = v => typeof v === "string" && /^[a-f0-9]{64}$/.test(v);
export function validParams(p) {
  return p && !Array.isArray(p) && LEGACY_KEYS.every(k => Object.hasOwn(p, k)) && Object.keys(p).every(k => Object.hasOwn(FIELDS, k))
    && Object.entries(FIELDS).every(([k, f]) => p[k] === undefined || Number.isFinite(p[k]) && p[k] >= f.min && p[k] <= f.max && (k !== "fins" || Number.isInteger(p[k])));
}
const revisionBody = c => ({ parents: c.parents, branch: c.branch, params: c.params, author: c.author, domain: c.domain, message: c.message, time: c.time, ...(c.features === undefined ? {} : { features: c.features }), ...(c.drawings === undefined ? {} : { drawings: c.drawings }), ...(c.studios === undefined ? {} : { studios: c.studios }) });
// Every other domain's document is a content-addressed blob a revision
// refers to by hash, so an unchanged board or building costs nothing to keep.
export const studiosOf = c => c.studios || {};
const blobBody = b => ({ kind: b.kind, doc: b.doc });
const noteBody = n => ({ commit: n.commit, author: n.author, domain: n.domain, text: n.text, kind: n.kind, time: n.time, nonce: n.nonce });
function withinSize(commits, notes, blobs = new Map()) {
  const bytes = new TextEncoder().encode(JSON.stringify({ format: FORMAT, project: PROJECT, commits: [...commits.values()], notes: [...notes.values()], blobs: [...blobs.values()] })).length;
  if (bytes > 3_000_000) throw new Error("The combined project exceeds this demo’s 3 MB limit.");
}
export class Repository {
  commits = new Map(); notes = new Map(); blobs = new Map(); root;
  // Set by the workspace: kind → validator returning the projected document.
  blobKinds = {};
  async putBlob(kind, doc) {
    if (!this.blobKinds[kind]) throw new Error(`Unknown design document kind ${kind}.`);
    const body = { kind, doc: this.blobKinds[kind](structuredClone(doc)) }, id = await hash(body);
    if (!this.blobs.has(id)) { withinSize(this.commits, this.notes, new Map([...this.blobs, [id, { ...body, id }]])); this.blobs.set(id, { ...body, id }); }
    return id;
  }
  async putStudios(docs) { const out = {}; for (const [kind, doc] of Object.entries(docs).sort()) out[kind] = await this.putBlob(kind, doc); return out; }
  studioDocs(commit) { return Object.fromEntries(Object.entries(studiosOf(commit)).map(([kind, id]) => [kind, structuredClone(this.blobs.get(id).doc)])); }
  async init() {
    const root = await this.make({ parents: [], branch: "main", params: LEGACY_DEFAULTS, author: "EvolveTech", domain: "Systems", message: "Motor-driver reference design", time: "2026-10-25T00:00:00.000Z" });
    this.root = root.id; this.commits.set(root.id, root); return this;
  }
  async make(body) { return { ...body, id: await hash(body) }; }
  async commit(branch, parents, params, identity, message, features, drawings, studios, time = new Date().toISOString()) {
    if (!validName(branch) || !validParams(params) || !validText(message, 240) || !validText(identity.author, 48)
      || !validText(identity.domain, 48)) throw new Error("Check the branch, author, message and parameter ranges.");
    if (!parents.length || parents.length > 2 || !parents.every(p => this.commits.has(p))) throw new Error("Revision parents are missing.");
    if (features !== undefined && !validFeatures(features)) throw new Error("Part definitions exceed the supported shape, dimension or count limits.");
    if (drawings !== undefined && !validDrawings(drawings)) throw new Error("Check drawing fields. A released sheet requires material, finish, tolerance, drawing number, revision, preparer and checker.");
    if (studios !== undefined && !Object.values(studios).every(id => this.blobs.has(id))) throw new Error("A design document is missing from this project.");
    if (this.commits.size >= 512) throw new Error("This demo supports 512 revisions. Export your project to keep a backup.");
    const c = await this.make({ parents: [...new Set(parents)], branch, params: { ...params }, author: identity.author, domain: identity.domain, message: message.trim(), time, ...(features === undefined ? {} : { features: structuredClone(features) }), ...(drawings === undefined ? {} : { drawings: structuredClone(drawings) }), ...(studios === undefined ? {} : { studios: { ...studios } }) });
    withinSize(new Map([...this.commits, [c.id, c]]), this.notes, this.blobs);
    this.commits.set(c.id, c); return c;
  }
  async note(commit, identity, text, kind = "comment") {
    if (!this.commits.has(commit) || !validText(text, 2000) || !["comment", "approve", "changes"].includes(kind)
      || !validText(identity.author, 48) || !validText(identity.domain, 48)) throw new Error("Enter a name and a review of up to 2,000 characters.");
    if (this.notes.size >= 1000) throw new Error("This demo supports 1,000 review notes.");
    const body = { commit, ...identity, text: text.trim(), kind, time: new Date().toISOString(), nonce: crypto.randomUUID() };
    const n = { ...body, id: await hash(body) }; withinSize(this.commits, new Map([...this.notes, [n.id, n]]), this.blobs); this.notes.set(n.id, n); return n;
  }
  ancestor(a, b, map = this.commits) {
    const queue = [b], seen = new Set();
    while (queue.length) { const id = queue.pop(); if (id === a) return true;
      if (!seen.has(id)) { seen.add(id); queue.push(...(map.get(id)?.parents || [])); }
    } return false;
  }
  branches() { return [...new Set([...this.commits.values()].map(c => c.branch))].sort((a, b) => a === "main" ? -1 : b === "main" ? 1 : a.localeCompare(b)); }
  heads(branch) {
    const cs = [...this.commits.values()].filter(c => c.branch === branch);
    return cs.filter(c => !cs.some(other => c.id !== other.id && this.ancestor(c.id, other.id))).map(c => c.id).sort();
  }
  mergeBase(a, b) {
    const common = [...this.commits.keys()].filter(id => this.ancestor(id, a) && this.ancestor(id, b));
    const bases = common.filter(id => !common.some(other => id !== other && this.ancestor(id, other)));
    if (!bases.length) throw new Error("These revisions have no shared history.");
    if (bases.length > 1) throw new Error("This history has multiple merge bases. This prototype requires manual reconciliation for criss-cross merges.");
    return this.commits.get(bases[0]);
  }
  studioDiff(a, b) { const sa = studiosOf(a), sb = studiosOf(b); return [...new Set([...Object.keys(sa), ...Object.keys(sb)])].filter(k => sa[k] !== sb[k]).sort(); }
  diff(a, b) { const av = values(a), bv = values(b); return Object.keys(FIELDS).filter(k => av[k] !== bv[k]); }
  featureDiff(a = [], b = []) { const am = new Map(a.map(f => [f.id, f])), bm = new Map(b.map(f => [f.id, f])); return [...new Set([...am.keys(), ...bm.keys()])].filter(id => canonical(am.get(id)) !== canonical(bm.get(id))); }
  merge(targetId, sourceId) {
    const target = this.commits.get(targetId), source = this.commits.get(sourceId), base = this.mergeBase(targetId, sourceId);
    const tv = values(target.params), sv = values(source.params), bv = values(base.params);
    const params = { ...tv }, conflicts = [];
    for (const key of Object.keys(FIELDS)) {
      if (sv[key] === bv[key] || sv[key] === tv[key]) continue;
      if (tv[key] !== bv[key]) conflicts.push(key);
      else params[key] = sv[key];
    }
    const tm = new Map(featureValues(target).map(f => [f.id, f])), sm = new Map(featureValues(source).map(f => [f.id, f])), bm = new Map(featureValues(base).map(f => [f.id, f]));
    const fm = new Map(tm), featureConflicts = [];
    for (const id of new Set([...tm.keys(), ...sm.keys(), ...bm.keys()])) {
      if (canonical(sm.get(id)) === canonical(bm.get(id)) || canonical(sm.get(id)) === canonical(tm.get(id))) continue;
      if (canonical(tm.get(id)) !== canonical(bm.get(id))) featureConflicts.push(id);
      else if (sm.has(id)) fm.set(id, sm.get(id)); else fm.delete(id);
    }
    const td = drawingValues(target), sd = drawingValues(source), bd = drawingValues(base), drawings = { ...td }, drawingConflicts = [];
    for (const id of new Set([...Object.keys(td),...Object.keys(sd),...Object.keys(bd)])) {
      if (canonical(sd[id]) === canonical(bd[id]) || canonical(sd[id]) === canonical(td[id])) continue;
      if (canonical(td[id]) !== canonical(bd[id])) drawingConflicts.push(id);
      else if (sd[id]) drawings[id] = sd[id]; else delete drawings[id];
    }
    // Domain documents merge by object id; both-sides edits become conflicts.
    const kinds = new Set([...Object.keys(studiosOf(target)), ...Object.keys(studiosOf(source))]), studios = {}, studioConflicts = [];
    const doc = (c, k) => studiosOf(c)[k] ? this.blobs.get(studiosOf(c)[k]).doc : undefined;
    for (const k of kinds) {
      const t = doc(target, k), s = doc(source, k), b = doc(base, k);
      const r = merge3(b, t ?? s, s ?? t);
      studios[k] = r.value; for (const c of r.conflicts) studioConflicts.push({ kind: k, ...c });
    }
    return { target, source, base, params, conflicts, features: [...fm.values()], featureConflicts, drawings, drawingConflicts, studios, studioConflicts };
  }
  bundle() { return { format: FORMAT, project: PROJECT, commits: [...this.commits.values()], notes: [...this.notes.values()], blobs: [...this.blobs.values()] }; }
  async ingest(bundle) {
    if (bundle?.format !== FORMAT || bundle?.project !== PROJECT || !Array.isArray(bundle.commits) || !Array.isArray(bundle.notes)
      || bundle.commits.length > 512 || bundle.notes.length > 1000) throw new Error("This is not an Emergence project history.");
    const staged = new Map(this.commits), notes = new Map(this.notes), blobs = new Map(this.blobs);
    if (bundle.blobs !== undefined && (!Array.isArray(bundle.blobs) || bundle.blobs.length > 4096)) throw new Error("This project's design documents are malformed.");
    for (const b of bundle.blobs || []) {
      if (!b || !isHash(b.id) || !this.blobKinds[b.kind]) throw new Error("A design document has an unknown kind. Nothing was imported.");
      let doc; try { doc = this.blobKinds[b.kind](structuredClone(b.doc)); } catch (e) { throw new Error(`A ${b.kind} document is invalid (${e.message}). Nothing was imported.`); }
      if (canonical(doc) !== canonical(b.doc) || await hash(blobBody(b)) !== b.id) throw new Error("A design document failed content verification. Nothing was imported.");
      blobs.set(b.id, { kind: b.kind, doc, id: b.id });
    }
    for (const c of bundle.commits) {
      if (!c || !isHash(c.id) || !Array.isArray(c.parents) || c.parents.length > 2 || !c.parents.every(isHash)
        || !validName(c.branch) || !validParams(c.params) || c.features !== undefined && !validFeatures(c.features) || c.drawings !== undefined && !validDrawings(c.drawings) || !validText(c.author, 48) || !validText(c.domain, 48)
        || !validText(c.message, 240) || !validText(c.time, 40) || !Number.isFinite(Date.parse(c.time))
        || c.studios !== undefined && (typeof c.studios !== "object" || Array.isArray(c.studios) || !Object.entries(c.studios).every(([k, id]) => blobs.get(id)?.kind === k))
        || await hash(revisionBody(c)) !== c.id) throw new Error("A revision failed content verification. Nothing was imported.");
      staged.set(c.id, { ...revisionBody(c), id: c.id });
    }
    if (staged.size > 512) throw new Error("The combined history exceeds this demo’s 512-revision limit.");
    // Validate topology before mutation. Cycles, detached roots and missing
    // parents are rejected, including malicious precomputed circular input.
    const done = new Set(), visiting = new Set();
    const visit = id => {
      if (visiting.has(id)) throw new Error("Cyclic project history.");
      if (done.has(id)) return;
      const c = staged.get(id); if (!c) throw new Error("Incomplete project history.");
      if (!c.parents.length && id !== this.root) throw new Error("Unrelated project root.");
      visiting.add(id); c.parents.forEach(visit); visiting.delete(id); done.add(id);
    };
    for (const id of staged.keys()) visit(id);
    for (const n of bundle.notes) {
      if (!n || !isHash(n.id) || !staged.has(n.commit) || !validText(n.author, 48) || !validText(n.domain, 48)
        || !validText(n.text, 2000) || !validText(n.nonce, 48) || !validText(n.time, 40) || !Number.isFinite(Date.parse(n.time))
        || !["comment", "approve", "changes"].includes(n.kind) || await hash(noteBody(n)) !== n.id) throw new Error("A review failed content verification. Nothing was imported.");
      notes.set(n.id, { ...noteBody(n), id: n.id });
    }
    if (notes.size > 1000) throw new Error("The combined project exceeds 1,000 review notes.");
    withinSize(staged, notes, blobs);
    const changed = staged.size !== this.commits.size || notes.size !== this.notes.size || blobs.size !== this.blobs.size;
    this.commits = staged; this.notes = notes; this.blobs = blobs; return changed;
  }
}
