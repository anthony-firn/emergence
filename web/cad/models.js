// One model Worker for every studio. Requests are serialized by the Worker;
// a timeout or crash rejects everything pending and the next call restarts it.
let worker, requestId = 0;
const pending = new Map(), listeners = new Set();
export const onModelStatus = fn => listeners.add(fn);
function fail(error) {
  worker?.terminate(); worker = null;
  for (const p of pending.values()) { clearTimeout(p.timer); p.reject(error); } pending.clear();
  for (const fn of listeners) fn({ ok: false, error });
}
function start() {
  worker = new Worker(new URL("../worker.js", import.meta.url), { type: "module" });
  // State preloaded into a Worker (chip corner libraries) dies with it.
  request.generation = (request.generation ?? 0) + 1;
  worker.onmessage = ({ data }) => { const p = pending.get(data.id); if (!p) return; clearTimeout(p.timer); pending.delete(data.id);
    data.error || data.result?.error ? p.reject(new Error(data.error || data.result.error)) : p.resolve(data); };
  worker.onerror = e => fail(new Error(e.message || "The model Worker failed"));
}
// A long job gets a Worker of its own: solves carry on beside it, and their time
// limit cannot end it.
export function requestLong(type, args = {}, seconds = 300) {
  return new Promise((resolve, reject) => {
    const w = new Worker(new URL("../worker.js", import.meta.url), { type: "module" });
    const end = fn => value => { clearTimeout(timer); w.terminate(); fn(value); }, ok = end(resolve), no = end(reject), timer = setTimeout(() => no(new Error(`Model exceeded the ${seconds}-second time limit`)), seconds * 1000);
    w.onmessage = ({ data }) => { if (data.id === 1) data.error || data.result?.error ? no(new Error(data.error || data.result.error)) : ok(data); };
    w.onerror = e => no(new Error(e.message || "The model Worker failed"));
    w.postMessage({ id: 1, type, ...args });
  });
}
export function request(type, args = {}, transfer = []) {
  try { if (!worker) start(); } catch (e) { return Promise.reject(e); }
  return new Promise((resolve, reject) => {
    const id = ++requestId, timer = setTimeout(() => fail(new Error("Model exceeded the 45-second time limit")), 45000);
    pending.set(id, { resolve, reject, timer }); worker.postMessage({ id, type, ...args }, transfer);
  });
}
