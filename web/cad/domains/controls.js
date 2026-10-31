// Flight control: PX4's multicopter cascade flown in a 6-DoF simulation, and the
// same loop analysed in frequency. The structure and parameter names are PX4's
// (position P → velocity PID → tilt and thrust → attitude P → rate PID → mixer),
// at its 400 Hz rate loop with its 2nd-order Butterworth gyro and D-term filters.
// This is a model of that controller, not PX4 itself: the estimator is perfect
// apart from gyro noise, and there is no sensor fusion, failsafe or mission logic.
// Motors lag by a first-order time constant from real motor constants; rotor
// aerodynamics are the thrust law through the hover and full-throttle points.
export const PX4_DEFAULTS = { MC_ROLL_P: 4, MC_PITCH_P: 4, MC_YAW_P: 2.8, MC_ROLLRATE_P: .15, MC_ROLLRATE_I: .2, MC_ROLLRATE_D: .003, MC_ROLLRATE_K: 1, MC_PITCHRATE_P: .15, MC_PITCHRATE_I: .2, MC_PITCHRATE_D: .003, MC_PITCHRATE_K: 1, MC_YAWRATE_P: .2, MC_YAWRATE_I: .1, MC_YAWRATE_D: 0, MC_YAWRATE_K: 1,
  MC_ROLLRATE_MAX: 220, MC_PITCHRATE_MAX: 220, MC_YAWRATE_MAX: 200, MPC_XY_P: .95, MPC_XY_VEL_P_ACC: 1.8, MPC_XY_VEL_I_ACC: .4, MPC_XY_VEL_D_ACC: .2, MPC_Z_P: 1, MPC_Z_VEL_P_ACC: 4, MPC_Z_VEL_I_ACC: 2, MPC_Z_VEL_D_ACC: 0, MPC_TILTMAX_AIR: 45, IMU_GYRO_CUTOFF: 40, IMU_DGYRO_CUTOFF: 20 };
export const PX4_SOURCE = "https://docs.px4.io/main/en/advanced_config/parameter_reference.html";
const PLANT = { mass: 2, ixx: .03, iyy: .03, izz: .06, arm: 250, rotors: 4, hoverThrottle: 55, maxThrust: 10, tau: .04, kappa: .016, dragArea: .05 };
const finite = (v, lo, hi) => Number.isFinite(v) && v >= lo && v <= hi;
export const controlsExample = () => ({ loopRate: 400, gyroNoise: .0028, gust: 8, wind: 3, stepAngle: 15, plant: { ...PLANT }, params: { ...PX4_DEFAULTS }, limits: { phaseMargin: 45, gainMargin: 6, overshoot: 25, settle: 1.5, drift: 3 } });
export function validateControls(d) {
  if (!d || !finite(d.loopRate, 50, 8000) || !finite(d.gyroNoise, 0, 1) || !finite(d.gust, 0, 50) || !finite(d.wind, 0, 50) || !finite(d.stepAngle, 1, 45) || !d.plant || !d.params || !d.limits) throw new Error("Invalid control settings.");
  for (const k of Object.keys(PLANT)) if (!finite(d.plant[k], 0, 1e6)) throw new Error(`Plant ${k} must be a number.`);
  if (![4, 6, 8].includes(d.plant.rotors)) throw new Error("The controller flies 4, 6 or 8 rotors.");
  for (const k of Object.keys(PX4_DEFAULTS)) if (!finite(d.params[k], 0, 1e4)) throw new Error(`${k} must be a number.`);
  for (const k of ["phaseMargin", "gainMargin", "overshoot", "settle", "drift"]) if (!finite(d.limits[k], 0, 1e4)) throw new Error(`Limit ${k} must be a number.`);
  return { loopRate: d.loopRate, gyroNoise: d.gyroNoise, gust: d.gust, wind: d.wind, stepAngle: d.stepAngle, plant: Object.fromEntries(Object.keys(PLANT).map(k => [k, d.plant[k]])), params: Object.fromEntries(Object.keys(PX4_DEFAULTS).map(k => [k, d.params[k]])),
    limits: { phaseMargin: d.limits.phaseMargin, gainMargin: d.limits.gainMargin, overshoot: d.limits.overshoot, settle: d.limits.settle, drift: d.limits.drift } };
}
// Linked values replace the plant's own; the rotor layout and thrust law follow.
export function plantOf(d, ov = {}) {
  const k = { mass: "mass", ixx: "ixx", iyy: "iyy", izz: "izz", arm: "arm", rotors: "rotors", hoverThrottle: "hover_throttle", maxThrust: "max_thrust", tau: "motor_tau", kappa: "torque_ratio" }, P = { ...d.plant };
  for (const [key, port] of Object.entries(k)) if (ov[port] !== undefined) P[key] = key === "rotors" ? Math.round(ov[port]) : ov[port];
  const n = P.rotors, R = P.arm / 1000, u = P.hoverThrottle / 100, Th = P.mass * 9.80665 / n;
  // T(u) = Tmax·u^a through the hover and full-throttle points, with u the ESC
  // command, as PX4 outputs it with THR_MDL_FAC = 0.
  const a = Math.log(Math.max(1e-3, Math.min(.999, Th / P.maxThrust))) / Math.log(Math.max(.01, Math.min(.99, u)));
  const rotors = Array.from({ length: n }, (_, i) => { const t = (45 + 360 * i / n) * Math.PI / 180; return { x: R * Math.cos(t), y: R * Math.sin(t), dir: i % 2 ? 1 : -1 }; });
  // PX4 normalises the mix so roll and pitch share a scale: Σm² = N/2.
  const scale = Math.sqrt((n / 2) / rotors.reduce((s, r) => s + (r.y / R) ** 2, 0));
  for (const r of rotors) { r.mr = scale * r.y / R; r.mp = -scale * r.x / R; r.my = r.dir; }
  return { ...P, n, R, a, u0: u, Th, rotors, slope: P.maxThrust * a * u ** (a - 1) };
}
// Biquad low-pass, the Butterworth PX4 uses on the gyro and on the D term.
function butter(fc, fs) { const w = Math.tan(Math.PI * Math.min(fc, fs * .45) / fs), s = Math.SQRT2, n = 1 / (1 + s * w + w * w), b0 = w * w * n, a1 = 2 * (w * w - 1) * n, a2 = (1 - s * w + w * w) * n; let x1 = 0, x2 = 0, y1 = 0, y2 = 0, init = false; return x => { if (!init) { x1 = x2 = y1 = y2 = x; init = true; } const y = b0 * (x + 2 * x1 + x2) - a1 * y1 - a2 * y2; x2 = x1; x1 = x; y2 = y1; y1 = y; return y; }; }
function random(seed) { let s = seed >>> 0; const u = () => { s = (s + 0x6D2B79F5) >>> 0; let t = s; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; return () => Math.sqrt(-2 * Math.log(u() + 1e-12)) * Math.cos(2 * Math.PI * u()); }
const clamp = (x, lo, hi) => Math.min(hi, Math.max(lo, x)), D2R = Math.PI / 180, G = 9.80665;
const qmul = (a, b) => [a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3], a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2], a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1], a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]];
const qeuler = (r, p, y) => { const c = v => Math.cos(v / 2), s = v => Math.sin(v / 2); return [c(r) * c(p) * c(y) + s(r) * s(p) * s(y), s(r) * c(p) * c(y) - c(r) * s(p) * s(y), c(r) * s(p) * c(y) + s(r) * c(p) * s(y), c(r) * c(p) * s(y) - s(r) * s(p) * c(y)]; };
export const eulerOf = q => [Math.atan2(2 * (q[0] * q[1] + q[2] * q[3]), 1 - 2 * (q[1] ** 2 + q[2] ** 2)), Math.asin(clamp(2 * (q[0] * q[2] - q[3] * q[1]), -1, 1)), Math.atan2(2 * (q[0] * q[3] + q[1] * q[2]), 1 - 2 * (q[2] ** 2 + q[3] ** 2))];
// One flight. `plan(t, state)` returns the mode and setpoint; returns traces.
export function simulate(d, P, plan, { seconds = 8, failed = -1, failAt = Infinity, thrustScale = 1, noise = true, record = v => v, wind = () => [0, 0, 0], open = null, dtSim } = {}) {
  const g = d.params, fs = d.loopRate, sub = 3, dt = dtSim ?? 1 / (fs * sub), steps = Math.round(seconds / dt), rho = 1.2, n = P.n, I = [P.ixx, P.iyy, P.izz];
  const pos = [0, 0, 0], vel = [0, 0, 0]; let q = [1, 0, 0, 0]; const w = [0, 0, 0], speed = new Float64Array(n).fill(Math.sqrt(P.Th)), cmd = new Float64Array(n).fill(Math.sqrt(P.Th));
  const gf = [0, 1, 2].map(() => butter(g.IMU_GYRO_CUTOFF, fs)), df = [0, 1, 2].map(() => butter(g.IMU_DGYRO_CUTOFF, fs)), prev = [0, 0, 0], integ = [0, 0, 0], viz = [0, 0, 0], rnd = random(1234), gyroSigma = noise ? d.gyroNoise * D2R * Math.sqrt(fs / 2) : 0;
  const gains = [[g.MC_ROLL_P, g.MC_ROLLRATE_P, g.MC_ROLLRATE_I, g.MC_ROLLRATE_D, g.MC_ROLLRATE_K, g.MC_ROLLRATE_MAX], [g.MC_PITCH_P, g.MC_PITCHRATE_P, g.MC_PITCHRATE_I, g.MC_PITCHRATE_D, g.MC_PITCHRATE_K, g.MC_PITCHRATE_MAX], [g.MC_YAW_P, g.MC_YAWRATE_P, g.MC_YAWRATE_I, g.MC_YAWRATE_D, g.MC_YAWRATE_K, g.MC_YAWRATE_MAX]];
  let thrustSp = P.u0, qSp = [1, 0, 0, 0], saturated = 0, ticks = 0, prevVel = [0, 0, 0], u = [0, 0, 0];
  const out = [], every = Math.max(1, Math.round(.02 / dt));
  for (let step = 0; step < steps; step++) {
    const t = step * dt;
    if (step % sub === 0) {
      ticks++;
      const sp = plan(t, { pos, vel, q }) ?? {};
      if (sp.position || sp.velocity) {
        // Position → velocity → acceleration, then tilt and thrust, as in PX4.
        const vs = sp.position ? [0, 1, 2].map(i => clamp((sp.position[i] - pos[i]) * (i < 2 ? g.MPC_XY_P : g.MPC_Z_P), -12, 12)) : sp.velocity, acc = [0, 0, 0];
        for (let i = 0; i < 3; i++) { const kp = i < 2 ? g.MPC_XY_VEL_P_ACC : g.MPC_Z_VEL_P_ACC, ki = i < 2 ? g.MPC_XY_VEL_I_ACC : g.MPC_Z_VEL_I_ACC, kd = i < 2 ? g.MPC_XY_VEL_D_ACC : g.MPC_Z_VEL_D_ACC, e = vs[i] - vel[i], de = -(vel[i] - prevVel[i]) * fs;
          viz[i] = clamp(viz[i] + ki * e / fs, -6, 6); acc[i] = kp * e + viz[i] + kd * de; prevVel[i] = vel[i]; }
        const up = [acc[0], acc[1], G + acc[2]], lim = Math.tan(g.MPC_TILTMAX_AIR * D2R), h = Math.hypot(up[0], up[1]), cap = Math.max(1, up[2]) * lim; if (h > cap) { up[0] *= cap / h; up[1] *= cap / h; } if (up[2] < .2 * G) up[2] = .2 * G;
        const roll = Math.atan2(-up[1], up[2]), pitch = Math.atan2(up[0] , Math.hypot(up[1], up[2]));
        qSp = qeuler(roll, pitch, sp.yaw ?? 0); thrustSp = clamp(P.u0 * Math.hypot(...up) / G, .08, 1);
      } else { qSp = qeuler(sp.roll ?? 0, sp.pitch ?? 0, sp.yaw ?? 0); thrustSp = sp.thrust ?? P.u0; }
      const c = [q[0], -q[1], -q[2], -q[3]], qe = qmul(c, qSp), s = qe[0] >= 0 ? 2 : -2, rateSp = sp.rates ?? gains.map((k, i) => clamp(k[0] * s * qe[i + 1], -k[5] * D2R, k[5] * D2R));
      for (let i = 0; i < 3; i++) {
        const meas = gf[i](w[i] + rnd() * gyroSigma), rate = df[i]((meas - prev[i]) * fs), k = gains[i], e = rateSp[i] - meas; prev[i] = meas;
        if (Math.abs(u[i]) < .4) integ[i] = clamp(integ[i] + k[2] * e / fs, -.3, .3);
        u[i] = clamp(k[4] * (k[1] * e + integ[i] - k[3] * rate), -1, 1);
      }
      if (open) u = open(t, u);
      // Mix, then give up yaw before roll and pitch, as PX4's desaturation does.
      const raw = P.rotors.map(r => thrustSp + u[0] * r.mr + u[1] * r.mp), yawed = P.rotors.map((r, i) => raw[i] + .5 * u[2] * r.my);
      const hi = Math.max(...yawed), lo = Math.min(...yawed), shift = hi > 1 ? 1 - hi : lo < .05 ? .05 - lo : 0; if (hi > 1 || lo < .05) saturated++;
      for (let i = 0; i < n; i++) { const c = clamp(yawed[i] + shift, .05, 1); cmd[i] = Math.sqrt(P.maxThrust * thrustScale * c ** P.a); }
    }
    // Rotors: a first-order lag in speed, with thrust going as speed².
    let Fz = 0; const tau = [0, 0, 0];
    for (let i = 0; i < n; i++) { speed[i] += (cmd[i] - speed[i]) * dt / P.tau; const r = P.rotors[i], T = i === failed && t >= failAt ? 0 : speed[i] ** 2; Fz += T; tau[0] += r.y * T; tau[1] -= r.x * T; tau[2] += r.dir * P.kappa * T; }
    const [a, b, c, e] = q, zb = [2 * (b * e + a * c), 2 * (c * e - a * b), 1 - 2 * (b * b + c * c)], air = wind(t), rel = [vel[0] - air[0], vel[1] - air[1], vel[2] - air[2]], sp = Math.hypot(...rel), drag = .5 * rho * P.dragArea * sp / P.mass;
    for (let i = 0; i < 3; i++) { vel[i] += ((Fz * zb[i]) / P.mass - drag * rel[i] - (i === 2 ? G : 0)) * dt; pos[i] += vel[i] * dt; }
    const dw = [(tau[0] - (I[2] - I[1]) * w[1] * w[2]) / I[0], (tau[1] - (I[0] - I[2]) * w[2] * w[0]) / I[1], (tau[2] - (I[1] - I[0]) * w[0] * w[1]) / I[2]];
    for (let i = 0; i < 3; i++) w[i] += dw[i] * dt;
    const dq = qmul(q, [0, w[0], w[1], w[2]]); q = q.map((v, i) => v + .5 * dq[i] * dt); const m = Math.hypot(...q); q = q.map(v => v / m);
    if (step % every === 0) out.push(record({ t, pos: [...pos], vel: [...vel], euler: eulerOf(q), w: [...w] }));
  }
  return { out, saturated: saturated / Math.max(1, ticks), final: { pos: [...pos], vel: [...vel], q: [...q], w: [...w] } };
}
// Step metrics of y against target, from time t0 onward.
function stepMetrics(trace, t0, target) {
  const after = trace.filter(s => s.t >= t0), peak = after.reduce((m, s) => Math.max(m, s.y), -Infinity), last = [...after].reverse().find(s => Math.abs(s.y - target) > .05 * Math.abs(target));
  const r10 = after.find(s => s.y >= .1 * target), r90 = after.find(s => s.y >= .9 * target);
  return { rise: r10 && r90 ? r90.t - r10.t : NaN, overshoot: Math.max(0, 100 * (peak - target) / target), settle: last ? last.t - t0 : 0, peak };
}
const thin = (rows, n = 240) => rows.filter((_, i) => i % Math.max(1, Math.floor(rows.length / n)) === 0);
export function flightTests(d, P) {
  const step = d.stepAngle * D2R, roll = simulate(d, P, t => ({ roll: t >= 1 ? step : 0 }), { seconds: 5 }), rollTrace = roll.out.map(s => ({ t: s.t, y: s.euler[0], r: s.t >= 1 ? step : 0 }));
  const alt = simulate(d, P, t => ({ position: [0, 0, t >= 1 ? 2 : 0] }), { seconds: 10 }), altTrace = alt.out.map(s => ({ t: s.t, y: s.pos[2], r: s.t >= 1 ? 2 : 0 }));
  // A steady wind, then a 1−cos gust, from t = 3 s.
  const air = t => { const g = t >= 3 && t <= 5 ? .5 * d.gust * (1 - Math.cos(Math.PI * (t - 3))) : 0; return [t >= 1 ? d.wind + g : 0, 0, 0]; };
  const gust = simulate(d, P, () => ({ position: [0, 0, 0] }), { seconds: 12, wind: air }), gustTrace = gust.out.map(s => ({ t: s.t, y: Math.hypot(s.pos[0], s.pos[1]), r: 0 })), drift = gustTrace.reduce((m, s) => Math.max(m, s.y), 0);
  // One rotor stops at 2 s. Only a vehicle with six or more can hold attitude.
  const fail = simulate(d, P, () => ({ position: [0, 0, 0] }), { seconds: 5, failed: 0, failAt: 2 }), failTrace = fail.out.map(s => ({ t: s.t, y: Math.max(Math.abs(s.euler[0]), Math.abs(s.euler[1])) / D2R, r: 0 })), tilt = failTrace.reduce((m, s) => Math.max(m, s.y), 0), z = fail.out.reduce((m, s) => Math.min(m, s.pos[2]), 0);
  return { rollStep: { ...stepMetrics(rollTrace, 1, step), saturated: roll.saturated, trace: thin(rollTrace) }, altitude: { ...stepMetrics(altTrace, 1, 2), trace: thin(altTrace) }, gust: { drift, saturated: gust.saturated, trace: thin(gustTrace) }, failure: { tilt, altitudeLoss: -z, survives: tilt < 45 && z > -3, trace: thin(failTrace) } };
}
// The rate loop in frequency: controller, gyro and D-term filters, delay, motor
// lag and the rigid body. Margins come from the loop as drawn, per axis.
const cx = (a, b) => ({ re: a, im: b }), mul = (a, b) => cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re), add = (a, b) => cx(a.re + b.re, a.im + b.im), inv = a => { const m = a.re ** 2 + a.im ** 2; return cx(a.re / m, -a.im / m); }, mag = a => Math.hypot(a.re, a.im), arg = a => Math.atan2(a.im, a.re);
const low2 = (f, w) => { const r = w / (2 * Math.PI * f); return inv(cx(1 - r * r, Math.SQRT2 * r)); };
export function rateLoop(d, P, axis = 0, gains = {}) {
  const g = { ...d.params, ...gains }, name = axis ? "PITCH" : "ROLL", I = axis ? P.iyy : P.ixx, kp = g[`MC_${name}RATE_P`], ki = g[`MC_${name}RATE_I`], kd = g[`MC_${name}RATE_D`], K = g[`MC_${name}RATE_K`];
  const G = P.slope * 2 * P.rotors.reduce((s, r) => s + Math.abs((axis ? r.x : r.y) * (axis ? r.mp : r.mr)), 0) / 2 / I, delay = .5 / d.loopRate; // The sample hold: PX4 acts on each gyro sample at once.
  const L = w => { const s = cx(0, w), C = add(add(cx(kp, 0), mul(cx(ki, 0), inv(s))), mul(mul(cx(kd, 0), s), low2(g.IMU_DGYRO_CUTOFF, w))); return mul(mul(mul(mul(mul(cx(K, 0), C), low2(g.IMU_GYRO_CUTOFF, w)), inv(cx(1, P.tau * w))), cx(Math.cos(w * delay), -Math.sin(w * delay))), mul(cx(G, 0), inv(s))); };
  let crossover = NaN, pm = NaN, gm = NaN, phasePrev = null, turn = 0, magPrev = null, wPrev = 0; const curve = [];
  for (let i = 0; i <= 400; i++) {
    const w = 10 ** (0 + 3.4 * i / 400), l = L(w), m = mag(l); let ph = arg(l) + turn; if (phasePrev !== null) while (ph > phasePrev + Math.PI) { ph -= 2 * Math.PI; turn -= 2 * Math.PI; }
    if (i % 8 === 0) curve.push({ w, gain: 20 * Math.log10(m), phase: ph * 180 / Math.PI });
    if (magPrev !== null && magPrev >= 1 && m < 1) { const f = Math.log(magPrev) / Math.log(magPrev / m); crossover = wPrev * (w / wPrev) ** f; pm = 180 + (phasePrev + (ph - phasePrev) * f) * 180 / Math.PI; }
    if (phasePrev !== null && Number.isNaN(gm) && phasePrev > -Math.PI && ph <= -Math.PI) gm = -20 * Math.log10(m);
    phasePrev = ph; magPrev = m; wPrev = w;
  }
  // Closed rate loop inside the attitude P: its own margin.
  const Latt = w => { const l = L(w), T = mul(l, inv(add(cx(1, 0), l))), s = cx(0, w); return mul(mul(cx(axis ? g.MC_PITCH_P : g.MC_ROLL_P, 0), inv(s)), T); };
  let attCross = NaN, attPm = NaN, mPrev = null; for (let i = 0; i <= 300; i++) { const w = 10 ** (-.5 + 3 * i / 300), l = Latt(w), m = mag(l); if (mPrev !== null && mPrev >= 1 && m < 1) { attCross = w; attPm = 180 + arg(l) * 180 / Math.PI; if (attPm > 180) attPm -= 360; } mPrev = m; }
  return { crossover, bandwidth: crossover / (2 * Math.PI), phaseMargin: pm, gainMargin: gm, curve, authority: G, attitude: { crossover: attCross, phaseMargin: attPm } };
}
// Gains for this airframe: the fastest rate loop meeting the margins, with I at a
// fifth of the crossover and D for phase lead. Parameters stay within PX4's ranges.
export function tune(d, P) {
  let best = null;
  for (const kp of Array.from({ length: 40 }, (_, i) => .02 * 1.1 ** i).filter(v => v <= .5)) for (const ratio of [0, .01, .02, .03, .04, .06]) {
    const round = v => Number(v.toPrecision(3)), gains = I => ({ MC_ROLLRATE_P: round(kp), MC_ROLLRATE_I: I, MC_ROLLRATE_D: round(kp * ratio), MC_ROLLRATE_K: 1 }), good = a => a.phaseMargin >= 55 && (Number.isNaN(a.gainMargin) || a.gainMargin >= 8);
    // The integrator's zero sits a fifth of the way to crossover, and the margins are read with it in.
    const open = rateLoop(d, P, 0, gains(0)); if (!good(open)) continue;
    const ki = round(Math.min(.5, kp * open.crossover / 5)), a = rateLoop(d, P, 0, gains(ki)); if (!good(a)) continue;
    if (!best || a.crossover > best.a.crossover) best = { a, set: { P: round(kp), I: ki, D: round(kp * ratio), K: 1 } };
  }
  if (!best) return null;
  const { set } = best;
  return { MC_ROLLRATE_P: set.P, MC_ROLLRATE_I: set.I, MC_ROLLRATE_D: set.D, MC_ROLLRATE_K: 1, MC_PITCHRATE_P: set.P, MC_PITCHRATE_I: set.I, MC_PITCHRATE_D: set.D, MC_PITCHRATE_K: 1, bandwidth: best.a.bandwidth };
}
// QGroundControl's parameter file, ready to load onto a real autopilot.
export function paramFile(d, name) {
  const lines = [`# Onboard parameters for ${name}, from Emergence CAD`, "# Vehicle-Id Component-Id Name Value Type"];
  for (const [k, v] of Object.entries(d.params).sort()) lines.push(`1\t1\t${k}\t${v}\t9`);
  return lines.join("\n") + "\n";
}
export const controlsAssumptions = d => [
  { id: "controls.model", label: "PX4's cascade and parameters, flown in a model", kind: "nominal", source: "Not PX4 itself: no estimator, failsafes or mixer details beyond roll, pitch and yaw allocation. Software-in-the-loop with real PX4 is the next step." },
  { id: "controls.rotor", label: "Rotor thrust follows speed² with a first-order lag", kind: "nominal", source: "The lag comes from motor constants and an estimated rotor inertia. An ESC step test measures it." },
  { id: "controls.drag", label: `Airframe drag area ${d.plant.dragArea} m²`, kind: "nominal", source: "An estimate; there is no rotor flapping or ground effect." },
  { id: "controls.noise", label: `Gyro noise ${d.gyroNoise} °/s/√Hz, no vibration`, kind: "nominal", source: "IMU white noise only; rotor vibration on the frame is not in this simulation." },
];
