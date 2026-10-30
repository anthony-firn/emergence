// A small WebGL2 renderer for kernel meshes: smooth normals, instancing by
// matrix, orbit/pan/zoom, a section plane, x-ray, colour-ID picking and PNG
// capture. No third-party code, under the site's script CSP.
const VS = `#version 300 es
in vec3 position; in vec3 normal; uniform mat4 model, view, proj; out vec3 vN; out vec3 vW;
void main(){ vec4 w = model * vec4(position,1.); vW = w.xyz; vN = mat3(model) * normal; gl_Position = proj * view * w; }`;
const FS = `#version 300 es
precision highp float; in vec3 vN; in vec3 vW; uniform vec4 color; uniform vec3 eye; uniform vec4 plane; uniform int picking; uniform vec3 pickColor; uniform float sky; out vec4 o;
void main(){
  if (plane.w > -1e8 && dot(vW, plane.xyz) > plane.w) discard;
  if (picking == 1) { o = vec4(pickColor, 1.); return; }
  vec3 n = normalize(gl_FrontFacing ? vN : -vN), v = normalize(eye - vW), l = normalize(vec3(.45, -.35, .82)), h = normalize(l + v);
  float hemi = .5 + .5 * n.z, diff = max(dot(n, l), 0.), spec = pow(max(dot(n, h), 0.), 48.) * .35, rim = pow(1. - max(dot(n, v), 0.), 3.) * .25;
  vec3 c = color.rgb * (mix(vec3(.30,.31,.34), vec3(.62,.66,.72), hemi) * sky + diff * .75) + spec + rim * vec3(.5,.7,.9);
  o = vec4(c, color.a);
}`;
const LVS = `#version 300 es
in vec3 position; uniform mat4 view, proj; void main(){ gl_Position = proj * view * vec4(position,1.); }`;
const LFS = `#version 300 es
precision highp float; uniform vec4 color; out vec4 o; void main(){ o = color; }`;
const mul = (a, b) => { const o = new Float32Array(16); for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) { let s = 0; for (let k = 0; k < 4; k++) s += a[k * 4 + j] * b[i * 4 + k]; o[i * 4 + j] = s; } return o; };
export const M = {
  identity: () => new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]),
  translate: (x, y, z) => new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, x, y, z, 1]),
  scale: (x, y = x, z = x) => new Float32Array([x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1]),
  rz: a => { const c = Math.cos(a), s = Math.sin(a); return new Float32Array([c, s, 0, 0, -s, c, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]); },
  rx: a => { const c = Math.cos(a), s = Math.sin(a); return new Float32Array([1, 0, 0, 0, 0, c, s, 0, 0, -s, c, 0, 0, 0, 0, 1]); },
  mul: (...ms) => ms.reduce((a, b) => mul(a, b)),
};
function perspective(fov, aspect, near, far) { const f = 1 / Math.tan(fov / 2), nf = 1 / (near - far); return new Float32Array([f / aspect, 0, 0, 0, 0, f, 0, 0, 0, 0, (far + near) * nf, -1, 0, 0, 2 * far * near * nf, 0]); }
function lookAt(e, t, u) {
  const z = norm(sub(e, t)), x = norm(cross(u, z)), y = cross(z, x);
  return new Float32Array([x[0], y[0], z[0], 0, x[1], y[1], z[1], 0, x[2], y[2], z[2], 0, -dot(x, e), -dot(y, e), -dot(z, e), 1]);
}
const sub = (a, b) => a.map((v, i) => v - b[i]), dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2], cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]], norm = a => { const l = Math.hypot(...a) || 1; return a.map(v => v / l); };
export function smoothNormals(pos, idx) {
  const n = new Float32Array(pos.length);
  for (let t = 0; t < idx.length; t += 3) {
    const [a, b, c] = [idx[t] * 3, idx[t + 1] * 3, idx[t + 2] * 3], u = [pos[b] - pos[a], pos[b + 1] - pos[a + 1], pos[b + 2] - pos[a + 2]], v = [pos[c] - pos[a], pos[c + 1] - pos[a + 1], pos[c + 2] - pos[a + 2]], f = cross(u, v);
    for (const k of [a, b, c]) { n[k] += f[0]; n[k + 1] += f[1]; n[k + 2] += f[2]; }
  }
  for (let i = 0; i < n.length; i += 3) { const l = Math.hypot(n[i], n[i + 1], n[i + 2]) || 1; n[i] /= l; n[i + 1] /= l; n[i + 2] /= l; }
  return n;
}
export function boxMesh(w, d, h) {
  const x = w / 2, y = d / 2, p = [], i = [];
  const faces = [[[0, 0, 1], [[-x, -y, h], [x, -y, h], [x, y, h], [-x, y, h]]], [[0, 0, -1], [[-x, y, 0], [x, y, 0], [x, -y, 0], [-x, -y, 0]]], [[1, 0, 0], [[x, -y, 0], [x, y, 0], [x, y, h], [x, -y, h]]],
    [[-1, 0, 0], [[-x, y, 0], [-x, -y, 0], [-x, -y, h], [-x, y, h]]], [[0, 1, 0], [[x, y, 0], [-x, y, 0], [-x, y, h], [x, y, h]]], [[0, -1, 0], [[-x, -y, 0], [x, -y, 0], [x, -y, h], [-x, -y, h]]]];
  const nrm = [];
  for (const [n, q] of faces) { const o = p.length / 3; for (const v of q) { p.push(...v); nrm.push(...n); } i.push(o, o + 1, o + 2, o, o + 2, o + 3); }
  return { positions: new Float32Array(p), indices: new Uint32Array(i), normals: new Float32Array(nrm) };
}
export function discMesh(r, segments = 64) {
  const p = [0, 0, 0], n = [0, 0, 1], i = [];
  for (let k = 0; k <= segments; k++) { const a = 2 * Math.PI * k / segments; p.push(r * Math.cos(a), r * Math.sin(a), 0); n.push(0, 0, 1); if (k) i.push(0, k, k + 1); }
  return { positions: new Float32Array(p), indices: new Uint32Array(i), normals: new Float32Array(n) };
}
export class Renderer {
  constructor(canvas) {
    this.canvas = canvas; const gl = canvas.getContext("webgl2", { antialias: true, preserveDrawingBuffer: true });
    if (!gl) throw new Error("WebGL2 is not available in this browser.");
    this.gl = gl; this.prog = this.program(VS, FS); this.lines = this.program(LVS, LFS);
    this.items = []; this.grid = null; this.yaw = -.7; this.pitch = .55; this.distance = 1200; this.target = [0, 0, 0]; this.plane = null; this.xray = false; this.sky = 1;
    gl.enable(gl.DEPTH_TEST); gl.getExtension("OES_element_index_uint");
    const pointer = { down: null };
    canvas.addEventListener("pointerdown", e => { pointer.down = [e.clientX, e.clientY, e.button, e.shiftKey]; pointer.moved = 0; canvas.setPointerCapture(e.pointerId); });
    canvas.addEventListener("pointermove", e => { if (!pointer.down) return; const dx = e.clientX - pointer.down[0], dy = e.clientY - pointer.down[1]; pointer.moved += Math.abs(dx) + Math.abs(dy); pointer.down[0] = e.clientX; pointer.down[1] = e.clientY;
      if (pointer.down[2] === 2 || pointer.down[3]) { const s = this.distance / 900, c = Math.cos(this.yaw), si = Math.sin(this.yaw); this.target[0] -= (dx * c + dy * si * Math.sin(this.pitch)) * s; this.target[1] -= (-dx * si + dy * c * Math.sin(this.pitch)) * s; this.target[2] += dy * Math.cos(this.pitch) * s; }
      else { this.yaw -= dx * .008; this.pitch = Math.max(-1.45, Math.min(1.45, this.pitch + dy * .008)); }
      this.draw(); });
    canvas.addEventListener("pointerup", e => { if (pointer.down && pointer.moved < 4) { const r = canvas.getBoundingClientRect(); this.onpick?.(this.pick(e.clientX - r.left, e.clientY - r.top)); } pointer.down = null; });
    canvas.addEventListener("contextmenu", e => e.preventDefault());
    canvas.addEventListener("wheel", e => { e.preventDefault(); this.distance = Math.max(1, this.distance * Math.exp(e.deltaY * .0012)); this.draw(); }, { passive: false });
    canvas.addEventListener("keydown", e => { const k = { ArrowLeft: () => this.yaw += .1, ArrowRight: () => this.yaw -= .1, ArrowUp: () => this.pitch = Math.min(1.45, this.pitch + .1), ArrowDown: () => this.pitch = Math.max(-1.45, this.pitch - .1), "+": () => this.distance *= .9, "=": () => this.distance *= .9, "-": () => this.distance *= 1.1 }[e.key]; if (k) { e.preventDefault(); k(); this.draw(); } });
    new ResizeObserver(() => this.draw()).observe(canvas);
  }
  program(vs, fs) {
    const gl = this.gl, p = gl.createProgram();
    for (const [type, src] of [[gl.VERTEX_SHADER, vs], [gl.FRAGMENT_SHADER, fs]]) { const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s); if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s)); gl.attachShader(p, s); }
    gl.linkProgram(p); if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p)); return p;
  }
  upload(mesh) {
    const gl = this.gl, vao = gl.createVertexArray(); gl.bindVertexArray(vao);
    const buf = (data, loc, size) => { const b = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, b); gl.bufferData(gl.ARRAY_BUFFER, data, gl.STATIC_DRAW); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, size, gl.FLOAT, false, 0, 0); };
    buf(mesh.positions, gl.getAttribLocation(this.prog, "position"), 3); buf(mesh.normals || smoothNormals(mesh.positions, mesh.indices), gl.getAttribLocation(this.prog, "normal"), 3);
    const ib = gl.createBuffer(); gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, ib); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, mesh.indices, gl.STATIC_DRAW); gl.bindVertexArray(null);
    return { vao, count: mesh.indices.length };
  }
  set(items, grid) { this.items = items; if (grid) this.setGrid(grid); this.draw(); }
  setGrid({ size, step, z = 0 }) {
    const gl = this.gl, p = []; for (let v = -size; v <= size + 1e-9; v += step) p.push(v, -size, z, v, size, z, -size, v, z, size, v, z);
    if (!this.gridVao) { this.gridVao = gl.createVertexArray(); this.gridBuf = gl.createBuffer(); }
    gl.bindVertexArray(this.gridVao); gl.bindBuffer(gl.ARRAY_BUFFER, this.gridBuf); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(p), gl.STATIC_DRAW);
    const loc = gl.getAttribLocation(this.lines, "position"); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, 3, gl.FLOAT, false, 0, 0); gl.bindVertexArray(null); this.gridCount = p.length / 3;
  }
  fit(min, max) { this.target = min.map((v, i) => (v + max[i]) / 2); this.distance = 1.6 * Math.hypot(...max.map((v, i) => v - min[i])) || 100; this.draw(); }
  camera() {
    const e = [this.target[0] + this.distance * Math.cos(this.pitch) * Math.sin(this.yaw), this.target[1] - this.distance * Math.cos(this.pitch) * Math.cos(this.yaw), this.target[2] + this.distance * Math.sin(this.pitch)];
    return { eye: e, view: lookAt(e, this.target, [0, 0, 1]) };
  }
  draw(picking = false) {
    const gl = this.gl, c = this.canvas, r = c.getBoundingClientRect(); if (!r.width || !r.height) return;
    const ratio = Math.min(devicePixelRatio || 1, 2), w = Math.round(r.width * ratio), h = Math.round(r.height * ratio);
    if (c.width !== w || c.height !== h) { c.width = w; c.height = h; }
    gl.viewport(0, 0, w, h);
    const bg = getComputedStyle(document.documentElement).getPropertyValue("--bg-2").trim() || "#0b1016", rgb = parseColor(bg);
    gl.clearColor(...(picking ? [0, 0, 0] : rgb), 1); gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    const { eye, view } = this.camera(), proj = perspective(.75, w / h, this.distance / 500, this.distance * 50);
    if (!picking && this.gridCount) { gl.useProgram(this.lines); gl.uniformMatrix4fv(gl.getUniformLocation(this.lines, "view"), false, view); gl.uniformMatrix4fv(gl.getUniformLocation(this.lines, "proj"), false, proj);
      gl.uniform4f(gl.getUniformLocation(this.lines, "color"), .45, .55, .62, .22); gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA); gl.bindVertexArray(this.gridVao); gl.drawArrays(gl.LINES, 0, this.gridCount); }
    const p = this.prog; gl.useProgram(p);
    const u = n => gl.getUniformLocation(p, n);
    gl.uniformMatrix4fv(u("view"), false, view); gl.uniformMatrix4fv(u("proj"), false, proj); gl.uniform3f(u("eye"), ...eye); gl.uniform1i(u("picking"), picking ? 1 : 0); gl.uniform1f(u("sky"), this.sky);
    gl.uniform4f(u("plane"), ...(this.plane || [0, 0, 1, -1e9]));
    const opaque = [], clear = [];
    this.items.forEach((it, k) => (picking || (!this.xray && it.color[3] >= 1) ? opaque : clear).push([it, k]));
    const one = ([it, k]) => { gl.uniformMatrix4fv(u("model"), false, it.matrix); gl.uniform4f(u("color"), it.color[0], it.color[1], it.color[2], this.xray && !picking ? Math.min(it.color[3], .28) : it.color[3]);
      gl.uniform3f(u("pickColor"), ((k + 1) & 255) / 255, (((k + 1) >> 8) & 255) / 255, (((k + 1) >> 16) & 255) / 255); gl.bindVertexArray(it.gpu.vao); gl.drawElements(gl.TRIANGLES, it.gpu.count, gl.UNSIGNED_INT, 0); };
    gl.disable(gl.BLEND); gl.depthMask(true); opaque.forEach(one);
    if (clear.length) { gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA); gl.depthMask(false); clear.forEach(one); gl.depthMask(true); }
    gl.bindVertexArray(null);
  }
  pick(x, y) {
    const gl = this.gl, ratio = this.canvas.width / this.canvas.getBoundingClientRect().width;
    this.draw(true); const px = new Uint8Array(4); gl.readPixels(Math.round(x * ratio), Math.round(this.canvas.height - y * ratio), 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, px); this.draw();
    const k = px[0] + (px[1] << 8) + (px[2] << 16) - 1; return k >= 0 ? this.items[k] : null;
  }
  png() { this.draw(); return new Promise(resolve => this.canvas.toBlob(resolve, "image/png")); }
}
function parseColor(c) {
  if (c.startsWith("#")) { const h = c.length === 4 ? c.slice(1).split("").map(x => x + x).join("") : c.slice(1, 7); return [0, 2, 4].map(i => parseInt(h.slice(i, i + 2), 16) / 255); }
  const m = c.match(/[\d.]+/g); return m ? m.slice(0, 3).map(v => Number(v) / 255) : [.04, .06, .08];
}
