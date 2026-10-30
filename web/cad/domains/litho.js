// Scalar Abbe imaging of a periodic mask clip: each source point tilts the
// mask spectrum through the circular pupil |f| ≤ NA/λ and the intensities
// add. A clear mask images to exactly 1; a coherent grating whose first order
// misses the pupil images to its DC term squared; both are tested. The resist
// is a sigmoid threshold. Inverse lithography descends the exact adjoint
// gradient of the printed-pattern error with respect to the mask pixels.
// Not modelled: thin-mask (Kirchhoff) errors at 3D mask edges, polarization,
// resist diffusion, defocus and flare. A study tool, not OPC sign-off.
export function fft(re, im, n, inverse = false) {
  for (let i = 1, j = 0; i < n; i++) { let bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit; if (i < j) { [re[i], re[j]] = [re[j], re[i]]; [im[i], im[j]] = [im[j], im[i]]; } }
  for (let len = 2; len <= n; len <<= 1) {
    const a = (inverse ? 2 : -2) * Math.PI / len, wr = Math.cos(a), wi = Math.sin(a);
    for (let i = 0; i < n; i += len) { let cr = 1, ci = 0;
      for (let k = 0; k < len / 2; k++) {
        const ur = re[i + k], ui = im[i + k], j = i + k + len / 2, vr = re[j] * cr - im[j] * ci, vi = re[j] * ci + im[j] * cr;
        re[i + k] = ur + vr; im[i + k] = ui + vi; re[j] = ur - vr; im[j] = ui - vi;
        const t = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = t;
      } }
  }
  if (inverse) for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
}
export function fft2(re, im, n, inverse = false) {
  const r = new Float64Array(n), c = new Float64Array(n);
  for (let y = 0; y < n; y++) { const o = y * n; r.set(re.subarray(o, o + n)); c.set(im.subarray(o, o + n)); fft(r, c, n, inverse); re.set(r, o); im.set(c, o); }
  for (let x = 0; x < n; x++) { for (let y = 0; y < n; y++) { r[y] = re[y * n + x]; c[y] = im[y * n + x]; } fft(r, c, n, inverse); for (let y = 0; y < n; y++) { re[y * n + x] = r[y]; im[y * n + x] = c[y]; } }
}
export function lithoExample() {
  // 150 nm features on KrF: k1 = CD·NA/λ ≈ 0.45.
  return { n: 64, pixel: 20, wavelength: 248, na: .75, sigmaOuter: .7, sigmaInner: 0, threshold: .3, steepness: 40, sourcePoints: 5,
    shapes: [{ x: 160, y: 200, w: 160, h: 880 }, { x: 480, y: 200, w: 640, h: 160 }, { x: 480, y: 360, w: 160, h: 400 }, { x: 860, y: 640, w: 180, h: 180 }] };
}
export function rasterize(d) {
  const n = d.n, m = new Float64Array(n * n);
  for (const s of d.shapes) for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) {
    const cx = (x + .5) * d.pixel, cy = (y + .5) * d.pixel;
    if (cx >= s.x && cx < s.x + s.w && cy >= s.y && cy < s.y + s.h) m[y * n + x] = 1;
  }
  return m;
}
// Source points on a square grid inside the (annular) σ disc, equal weights.
export function source(d) {
  const pts = [], k = d.sourcePoints;
  if (d.sigmaOuter === 0 || k <= 1) return [[0, 0, 1]];
  for (let i = 0; i < k; i++) for (let j = 0; j < k; j++) {
    const sx = (2 * (i + .5) / k - 1) * d.sigmaOuter, sy = (2 * (j + .5) / k - 1) * d.sigmaOuter, r = Math.hypot(sx, sy);
    if (r <= d.sigmaOuter + 1e-12 && r >= d.sigmaInner - 1e-12) pts.push([sx, sy, 1]);
  }
  const w = 1 / pts.length; return pts.map(([x, y]) => [x, y, w]);
}
function pupils(d) {
  const n = d.n, L = n * d.pixel, cutoff = d.na / d.wavelength, f = k => (k < n / 2 ? k : k - n) / L;
  return source(d).map(([sx, sy, w]) => {
    const p = new Float64Array(n * n), fx0 = sx * cutoff, fy0 = sy * cutoff;
    for (let v = 0; v < n; v++) for (let u = 0; u < n; u++) p[v * n + u] = Math.hypot(f(u) + fx0, f(v) + fy0) <= cutoff + 1e-15 ? 1 : 0;
    return { p, w };
  });
}
export function aerial(d, mask, cache = pupils(d)) {
  const n = d.n, N = n * n, mr = Float64Array.from(mask), mi = new Float64Array(N);
  fft2(mr, mi, n);
  const I = new Float64Array(N), fields = [];
  for (const { p, w } of cache) {
    const er = new Float64Array(N), ei = new Float64Array(N);
    for (let k = 0; k < N; k++) { er[k] = mr[k] * p[k]; ei[k] = mi[k] * p[k]; }
    fft2(er, ei, n, true);
    for (let k = 0; k < N; k++) I[k] += w * (er[k] * er[k] + ei[k] * ei[k]);
    fields.push({ er, ei, p, w });
  }
  return { I, fields };
}
const sig = x => 1 / (1 + Math.exp(-x));
export function resist(d, I) { return I.map(v => sig(d.steepness * (v - d.threshold))); }
export function printed(d, I) { return I.map(v => v >= d.threshold ? 1 : 0); }
export function patternError(d, target, I) { let e = 0; for (let k = 0; k < I.length; k++) if ((I[k] >= d.threshold ? 1 : 0) !== target[k]) e++; return e; }
// Gradient descent on θ with mask m = σ(βθ); exact adjoint of the Abbe model.
export function ilt(d, target, { iterations = 40, step = 4, beta = 4, onStep } = {}) {
  const n = d.n, N = n * n, cache = pupils(d), theta = target.map(t => (t - .5) * 2);
  let mask = theta.map(t => sig(beta * t)), history = [];
  for (let it = 0; it < iterations; it++) {
    const { I, fields } = aerial(d, mask, cache), R = resist(d, I), G = new Float64Array(N);
    let loss = 0;
    for (let k = 0; k < N; k++) { const e = R[k] - target[k]; loss += e * e; G[k] = 2 * e * d.steepness * R[k] * (1 - R[k]); }
    const grad = new Float64Array(N);
    for (const { er, ei, p, w } of fields) {
      const gr = new Float64Array(N), gi = new Float64Array(N);
      for (let k = 0; k < N; k++) { gr[k] = G[k] * er[k]; gi[k] = G[k] * ei[k]; }
      fft2(gr, gi, n);
      for (let k = 0; k < N; k++) { gr[k] *= p[k]; gi[k] *= p[k]; }
      fft2(gr, gi, n, true);
      for (let k = 0; k < N; k++) grad[k] += 2 * w * gr[k];
    }
    let gmax = 0; for (let k = 0; k < N; k++) { grad[k] *= beta * mask[k] * (1 - mask[k]); gmax = Math.max(gmax, Math.abs(grad[k])); }
    for (let k = 0; k < N; k++) theta[k] -= step * grad[k] / Math.max(gmax, 1e-12);
    mask = theta.map(t => sig(beta * t));
    history.push({ iteration: it, loss, errors: patternError(d, target, I) });
    onStep?.(history.at(-1));
  }
  const binary = mask.map(v => v >= .5 ? 1 : 0), final = aerial(d, binary, cache).I;
  return { mask, binary, I: final, history, errors: patternError(d, target, final) };
}
export function validateLitho(d) {
  const ok = (v, lo, hi) => Number.isFinite(v) && v >= lo && v <= hi;
  if (!d || ![32, 64, 128].includes(d.n) || !ok(d.pixel, 1, 200) || !ok(d.wavelength, 13, 500) || !ok(d.na, .1, 1.5) || !ok(d.sigmaOuter, 0, 1) || !ok(d.sigmaInner, 0, d.sigmaOuter)
    || !ok(d.threshold, .05, .95) || !ok(d.steepness, 1, 200) || !Number.isInteger(d.sourcePoints) || d.sourcePoints < 1 || d.sourcePoints > 9
    || !Array.isArray(d.shapes) || d.shapes.length > 32 || !d.shapes.every(s => s && ["x", "y", "w", "h"].every(k => ok(s[k], 0, 1e5)))) throw new Error("Invalid lithography clip settings.");
  return { n: d.n, pixel: d.pixel, wavelength: d.wavelength, na: d.na, sigmaOuter: d.sigmaOuter, sigmaInner: d.sigmaInner, threshold: d.threshold, steepness: d.steepness,
    sourcePoints: d.sourcePoints, shapes: d.shapes.map(s => ({ x: s.x, y: s.y, w: s.w, h: s.h })) };
}
