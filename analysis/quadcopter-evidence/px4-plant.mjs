// The aircraft as PX4's simulator takes it, from the solved project. SIH's thrust is
// linear in command, so T_MAX is the real aircraft's slope at hover: same authority.
import { writeFileSync } from "node:fs";
const { openProject } = await import(new URL("../../web/cad/agent/headless.mjs", import.meta.url));
const { plantOf, PX4_DEFAULTS } = await import(new URL("../../web/cad/domains/controls.js", import.meta.url));
const { project: p } = await openProject(undefined, { example: "drone" }); const s = await p.solve(), o = k => s.outputs[k].value, d = p.store.get("controls");
const P = plantOf(d, { mass: o("flight.total_mass"), ixx: o("flight.ixx"), iyy: o("flight.iyy"), izz: o("flight.izz"), arm: o("flight.arm_length"), rotors: 4, hover_throttle: o("flight.hover_throttle"), max_thrust: o("flight.max_thrust"), motor_tau: o("flight.motor_tau"), torque_ratio: o("flight.torque_ratio") });
const gains = k => ({ MC_ROLLRATE_P: k.MC_ROLLRATE_P, MC_ROLLRATE_I: k.MC_ROLLRATE_I, MC_ROLLRATE_D: k.MC_ROLLRATE_D, MC_ROLLRATE_K: k.MC_ROLLRATE_K, MC_PITCHRATE_P: k.MC_PITCHRATE_P, MC_PITCHRATE_I: k.MC_PITCHRATE_I, MC_PITCHRATE_D: k.MC_PITCHRATE_D, MC_PITCHRATE_K: k.MC_PITCHRATE_K });
const cfg = { mass: P.mass, ixx: P.ixx, iyy: P.iyy, izz: P.izz, tmax: P.slope, qmax: P.slope * P.kappa, lroll: P.R * Math.SQRT1_2, tau: P.tau, kappa: P.kappa, arm: P.arm, realMaxThrust: P.maxThrust, realHover: P.u0, hover: P.Th / P.slope, step: d.stepAngle, gains: { tuned: gains(d.params), stock: gains(PX4_DEFAULTS) } };
writeFileSync(process.argv[2], JSON.stringify(cfg, null, 1)); console.log(cfg); process.exit(0);
