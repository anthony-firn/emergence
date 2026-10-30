// The working copy of every non-mechanical design document. Commands mutate a
// clone, the domain validator projects it to known fields, and only then does
// it replace the document, so a failed edit never leaves a partial change.
export class Store {
  constructor(domains) { this.domains = domains; this.docs = {}; this.listeners = new Set(); for (const [id, d] of Object.entries(domains)) this.docs[id] = d.validate(d.example()); }
  get(id) { return this.docs[id]; }
  ids() { return Object.keys(this.domains); }
  on(fn) { this.listeners.add(fn); return () => this.listeners.delete(fn); }
  emit(ids, reason) { for (const fn of this.listeners) fn({ ids, reason }); }
  mutate(id, fn) {
    const draft = structuredClone(this.docs[id]), result = fn(draft);
    this.docs[id] = this.domains[id].validate(draft); this.emit([id], "edit"); return result;
  }
  replace(id, doc, reason = "replace") { this.docs[id] = this.domains[id].validate(structuredClone(doc)); this.emit([id], reason); }
  // Whole-project loads (checkout, import, peer fast-forward) validate first.
  load(docs, reason = "load") {
    const next = { ...this.docs };
    for (const [id, doc] of Object.entries(docs)) if (this.domains[id]) next[id] = this.domains[id].validate(structuredClone(doc));
    this.docs = next; this.emit(Object.keys(docs), reason);
  }
  snapshot(id) { return structuredClone(this.docs[id]); }
  restore(id, doc) { this.docs[id] = doc; this.emit([id], "restore"); }
  all() { return structuredClone(this.docs); }
}
