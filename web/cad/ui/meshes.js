// Kernel meshes for the viewers, cached by model and parameters.
import { request } from "../models.js";
const cache = new Map();
export function modelMesh(model, overrides = {}, budget = 2.5e5) {
  const key = JSON.stringify([model, overrides, budget]);
  if (!cache.has(key)) { cache.set(key, request("model", { input: { model, overrides, budget, binary: true } }).then(r => r.result)); if (cache.size > 64) cache.delete(cache.keys().next().value); }
  const p = cache.get(key); p.catch(() => cache.delete(key)); return p;
}
export function sinkMesh(params) {
  const key = JSON.stringify(["sink", params]);
  if (!cache.has(key)) cache.set(key, request("mesh", { params, target: "sink" }).then(({ result: m }) => ({
    ...m, positions: Float32Array.from(m.vertices.flat()), indices: Uint32Array.from(m.triangles.flat()), volume: m.volume_mm3, triangle_count: m.triangles.length })));
  const p = cache.get(key); p.catch(() => cache.delete(key)); return p;
}
export function stl(mesh, name) {
  const n = mesh.indices.length / 3, buffer = new ArrayBuffer(84 + n * 50), data = new DataView(buffer), P = mesh.positions;
  new Uint8Array(buffer).set(new TextEncoder().encode(`Emergence CAD; millimetres; ${name}`).slice(0, 80)); data.setUint32(80, n, true);
  let o = 84;
  for (let t = 0; t < n; t++) {
    const [a, b, c] = [0, 1, 2].map(k => mesh.indices[3 * t + k] * 3), u = [0, 1, 2].map(k => P[b + k] - P[a + k]), v = [0, 1, 2].map(k => P[c + k] - P[a + k]);
    const nrm = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]], l = Math.hypot(...nrm) || 1;
    for (const x of [...nrm.map(x => x / l), P[a], P[a + 1], P[a + 2], P[b], P[b + 1], P[b + 2], P[c], P[c + 1], P[c + 2]]) { data.setFloat32(o, x, true); o += 4; }
    data.setUint16(o, 0, true); o += 2;
  }
  return buffer;
}
