import { readFileSync, writeFileSync } from "node:fs";
const S = (process.env.EVIDENCE_WORK ?? "/tmp/emergence-evidence") + "/px4-run/";
const { controlsExample, validateControls, plantOf, simulate, rateLoop } = await import(new URL("../../web/cad/domains/controls.js", import.meta.url));
const cfg = JSON.parse(readFileSync(S + "plant.json", "utf8")), D2R = Math.PI / 180, step = cfg.step * D2R;
const rollOf = q => Math.atan2(2 * (q[0] * q[1] + q[2] * q[3]), 1 - 2 * (q[1] ** 2 + q[2] ** 2));
const at = (t, v, x) => { let lo = 0, hi = t.length - 1; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (t[m] <= x) lo = m; else hi = m; } return v[lo] + (v[hi] - v[lo]) * (x - t[lo]) / Math.max(1e-9, t[hi] - t[lo]); };
const metrics = (ts, ys, target) => { const peak = Math.max(...ys), r10 = ts[ys.findIndex(y => y >= .1 * target)], r90 = ts[ys.findIndex(y => y >= .9 * target)]; let last = -1; ys.forEach((y, i) => { if (Math.abs(y - target) > .05 * Math.abs(target)) last = i; }); return { rise: r90 - r10, overshoot: Math.max(0, 100 * (peak - target) / target), settle: last >= 0 ? ts[last] : 0, peak }; };
const cx = (a, b) => ({ re: a, im: b }), mul = (a, b) => cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re), div = (a, b) => { const m = b.re ** 2 + b.im ** 2; return cx((a.re * b.re + a.im * b.im) / m, (a.im * b.re - a.re * b.im) / m); }, add = (a, b) => cx(a.re + b.re, a.im + b.im), abs = a => Math.hypot(a.re, a.im), arg = a => Math.atan2(a.im, a.re);
const results = {};
for (const tag of ["stock", "tuned"]) {
  const L = JSON.parse(readFileSync(S + tag + ".json", "utf8")), fs = 250;
  const params = { ...controlsExample().params, ...cfg.gains[tag] }, d = validateControls({ ...controlsExample(), loopRate: fs, params });
  const P = plantOf(d, { mass: cfg.mass, ixx: cfg.ixx, iyy: cfg.iyy, izz: cfg.izz, arm: cfg.arm, rotors: 4, hover_throttle: 100 * cfg.hover, max_thrust: cfg.tmax, motor_tau: cfg.tau, torque_ratio: cfg.kappa });
  const sp = L.att_sp.t.map((t, i) => [t, rollOf([0, 1, 2, 3].map(k => L.att_sp[`q_d[${k}]`][i]))]), t0 = sp.find(([, r]) => r > step / 2)[0];
  const roll = L.att.t.map((_, i) => rollOf([0, 1, 2, 3].map(k => L.att[`q[${k}]`][i]))), grid = Array.from({ length: 400 }, (_, i) => i * .005);
  const px4 = grid.map(x => at(L.att.t, roll, t0 + x));
  const sim = simulate(d, P, t => ({ roll: t >= 1 ? step : 0, pitch: 0 }), { seconds: 3.1, noise: false, record: s => s }), o = sim.out, mine = grid.map(x => o[Math.min(o.length - 1, o.findIndex(s => s.t >= 1 + x))].euler[0]);
  const rms = Math.sqrt(grid.reduce((s, _, i) => s + (px4[i] - mine[i]) ** 2, 0) / grid.length) / D2R, m = metrics(grid, px4, step), n = metrics(grid, mine, step);
  // Rate steps: every fast departure from zero before the sweep, averaged by sign.
  const after = L.rate_sp.t.map((t, i) => [t, L.rate_sp.roll[i]]).filter(([t]) => t > t0 + 11), edges = []; for (let i = 1; i < after.length; i++) if (Math.abs(after[i - 1][1]) < 5 * D2R && Math.abs(after[i][1]) > 25 * D2R && (!edges.length || after[i][0] - edges.at(-1)[0] > .5)) edges.push([after[i][0], Math.sign(after[i][1])]);
  edges.length = 24; // The sweep follows; the flight commands exactly 24 steps.
  const rg = Array.from({ length: 76 }, (_, i) => i * .005), rate = rg.map(x => edges.reduce((s, [e, sign]) => s + sign * (at(L.rate.t, L.rate["xyz[0]"], e + x) - at(L.rate.t, L.rate["xyz[0]"], e - .02)), 0) / edges.length), rate30 = 30 * D2R, rm = metrics(rg, rate, rate30), t1 = edges.at(-1)[0];
  const rs = simulate(d, P, t => ({ rates: [t >= 1 ? rate30 : 0, 0, 0] }), { seconds: 1.45, noise: false, record: s => s }), ro = rs.out, rateSim = rg.map(x => ro[Math.min(ro.length - 1, ro.findIndex(s => s.t >= 1 + x))].w[0]), rn = metrics(rg, rateSim, rate30), rateRms = Math.sqrt(rg.reduce((s, _, i) => s + (rate[i] - rateSim[i]) ** 2, 0) / rg.length) / D2R;
  // Simulator's rate step: 30 °/s, from its rate-setpoint input (attitude P set to carry it).
  // Hover and mix, at full rate.
  const u = i => [0, 1, 2, 3].map(k => L.motors[`control[${k}]`][i]); let hov = 0, hn = 0, xs = 0, xy = 0;
  L.motors.t.forEach((t, i) => { if (t > t0 - 1 && t < t0 - .1) { const c = u(i); hov += c.reduce((s, v) => s + v, 0) / 4; hn++; } if (t > t0 - .1 && t < t0 + 1.5) { const c = u(i), x = at(L.torque.t, L.torque["xyz[0]"], t); xs += x * x; xy += x * (-c[0] + c[1] + c[2] - c[3]) / 4; } });
  // Closed-loop identification. Gyro noise feeds back into the torque command, so
  // torque and rate are both correlated with it; the setpoint is not. Cross-spectra
  // against the setpoint, averaged over overlapping segments, give an unbiased
  // plant and a coherence that says which frequencies can be trusted.
  const t2 = t1 + 7.7, dt = 1 / fs, N = Math.floor(28 * fs), R = new Float64Array(N), U = new Float64Array(N), Y = new Float64Array(N);
  for (let i = 0; i < N; i++) { const t = t2 + i * dt; R[i] = at(L.rate_sp.t, L.rate_sp.roll, t); U[i] = at(L.torque.t, L.torque["xyz[0]"], t); Y[i] = at(L.rate.t, L.rate["xyz[0]"], t); }
  const a = rateLoop(d, P), lowpass = (fc, f) => { const r = f / fc; return div(cx(1, 0), cx(1 - r * r, Math.SQRT2 * r)); }, W = 1024, hop = 128, hann = Float64Array.from({ length: W }, (_, i) => .5 - .5 * Math.cos(2 * Math.PI * i / W));
  const curve = []; let unwrap = 0, last = null;
  for (let f = 1; f < 15; f *= 1.04) {
    const w = 2 * Math.PI * f; let Sru = cx(0, 0), Sry = cx(0, 0), Srr = 0, Syy = 0;
    for (let s = 0; s + W <= N; s += hop) { let rr = 0, ri = 0, ur = 0, ui = 0, yr = 0, yi = 0; for (let i = 0; i < W; i++) { const q = w * i * dt, h = hann[i], c = Math.cos(q) * h, n = Math.sin(q) * h, j = s + i; rr += R[j] * c; ri -= R[j] * n; ur += U[j] * c; ui -= U[j] * n; yr += Y[j] * c; yi -= Y[j] * n; }
      const r = cx(rr, -ri); Sru = add(Sru, mul(cx(ur, ui), r)); Sry = add(Sry, mul(cx(yr, yi), r)); Srr += rr * rr + ri * ri; Syy += yr * yr + yi * yi; }
    const G = div(Sry, Sru), C = add(add(cx(params.MC_ROLLRATE_P, 0), cx(0, -params.MC_ROLLRATE_I / w)), mul(cx(0, w * params.MC_ROLLRATE_D), lowpass(params.IMU_DGYRO_CUTOFF, f))), Lp = mul(mul(cx(params.MC_ROLLRATE_K, 0), C), G), mag = abs(Lp), coherence = abs(Sry) ** 2 / (Srr * Syy); let ph = arg(Lp) * 180 / Math.PI + unwrap; if (last !== null) while (ph > last + 180) { ph -= 360; unwrap -= 360; } last = ph;
    curve.push({ f, mag, ph, coherence });
  }
  const crossing = (key, test, value) => { for (let i = 1; i < curve.length; i++) if (test(curve[i - 1][key]) && !test(curve[i][key])) { const x = (value(curve[i - 1][key]) - 0) / (value(curve[i - 1][key]) - value(curve[i][key])); return { i, x }; } return null; };
  const xc = crossing("mag", m => m >= 1, m => Math.log(m)), cross = xc && curve[xc.i - 1].f * (curve[xc.i].f / curve[xc.i - 1].f) ** xc.x, pm = xc && 180 + curve[xc.i - 1].ph + (curve[xc.i].ph - curve[xc.i - 1].ph) * xc.x;
  const xp = crossing("ph", p => p > -180, p => p + 180), gm = xp && -20 * Math.log10(curve[xp.i - 1].mag * (curve[xp.i].mag / curve[xp.i - 1].mag) ** xp.x), gmFreq = xp && curve[xp.i].f, gmCoh = xp && curve[xp.i].coherence, crossCoh = xc && curve[xc.i].coherence;
  results[tag] = { fs, hoverMean: hov / hn, hoverPredicted: cfg.hover, mix: xy / xs, roll: { px4: m, sim: n, rmsDegrees: rms }, rate: { px4: rm, sim: rn, rmsDegrees: rateRms, trace: { t: rg, px4: rate.map(v => v / D2R), sim: rateSim.map(v => v / D2R) } }, measured: { crossoverHz: cross, phaseMargin: pm, gainMargin: gm, gainMarginHz: gmFreq, crossCoherence: crossCoh, gainCoherence: gmCoh }, model: { crossoverHz: a.bandwidth, phaseMargin: a.phaseMargin, gainMargin: a.gainMargin }, curve, traces: { t: grid, px4: px4.map(v => v / D2R), sim: mine.map(v => v / D2R) } };
  const f3 = x => x == null ? "n/a" : x.toFixed(3), f1d = x => x == null ? "n/a" : x.toFixed(1);
  console.log(`${tag}: hover ${f3(hov / hn)} vs ${f3(cfg.hover)} | mix ${f3(xy / xs)} | roll rise ${f3(m.rise)} / ${f3(n.rise)}, overshoot ${f1d(m.overshoot)}% / ${f1d(n.overshoot)}%, settle ${f3(m.settle)} / ${f3(n.settle)}, RMS ${rms.toFixed(2)}° | rate step overshoot ${f1d(rm.overshoot)}% / ${f1d(rn.overshoot)}%, RMS ${rateRms.toFixed(2)}°/s | ${edges.length} steps averaged`);
  console.log(`   LOOP measured: crossover ${f3(cross)} Hz (coherence ${f3(crossCoh)}), PM ${f1d(pm)}°, GM ${f1d(gm)} dB at ${f1d(gmFreq)} Hz (coherence ${f3(gmCoh)}) | model: ${f3(a.bandwidth)} Hz, PM ${f1d(a.phaseMargin)}°, GM ${f1d(a.gainMargin)} dB`);
}
writeFileSync(S + "compare.json", JSON.stringify(results)); process.exit(0);
