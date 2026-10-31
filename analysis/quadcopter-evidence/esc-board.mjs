// The ESC power stage as commands: real parts, real pinouts, datasheet land patterns.
const H = await import(new URL("../../web/cad/agent/headless.mjs", import.meta.url));
export async function escBoard() {
  const { project: p } = await H.openProject(undefined, { example: "blank" }), run = (op, args) => p.bus.run(op, args);
  await run("project.add_domain", { domain: "pcb", start: "empty" });
  await run("pcb.board", { width: 78, height: 42, grid: .25, clearance: .2, min_width: .2, edge: .5, density: 50, thickness: 1.6, copper: 2, epsr: 4.4 });
  // Infineon PG-TDSON-8 recommended board pads (datasheet figure 2): pins 1–3 source, 4 gate, 5–8 one drain land.
  await run("pcb.footprint", { name: "tdson8", kind: "custom", smd: true, title: "PG-TDSON-8 (Infineon land pattern)", width: 5.15, height: 6.15, mass: .12, pads: [[-1.905, 2.4, .6, .925], [-.635, 2.4, .6, .925], [.635, 2.4, .6, .925], [1.905, 2.4, .6, .925], [0, -1.0975, 4.41, 4.455]] });
  await run("pcb.footprint", { name: "soic8", kind: "dual", pins: 8, pitch: 1.27, span: 5.4, smd: true, pad_w: 1.55, pad_h: .6, width: 6.9, height: 5, title: "SOIC-8", mass: .08 });
  await run("pcb.footprint", { name: "r0603", kind: "custom", smd: true, title: "0603", width: 2.4, height: 1.1, mass: .004, pads: [[-.8, 0, .8, .95], [.8, 0, .8, .95]] });
  await run("pcb.footprint", { name: "c1210", kind: "custom", smd: true, title: "1210", width: 4.6, height: 2.9, mass: .1, pads: [[-1.55, 0, 1.1, 2.7], [1.55, 0, 1.1, 2.7]] });
  await run("pcb.footprint", { name: "wirepad", kind: "custom", title: "Solder pad for a power wire", width: 4, height: 4, mass: .1, pad: 3.5, drill: 1.6, pads: [[0, 0]] });
  await run("pcb.footprint", { name: "header8", kind: "single", pins: 8, pitch: 2.54, title: "8-way 2.54 mm header", pad: 1.6, drill: 1, width: 21, height: 3, mass: .6 });
  // From KiCad's libraries: LQFP-32 0.8 mm, SOT-23-5, Bourns SRR1260, 2512, and the LM5017's PSOP-8 with its exposed pad.
  const side = (n, f) => Array.from({ length: n }, (_, i) => f(i));
  await run("pcb.footprint", { name: "lqfp32", kind: "custom", smd: true, title: "LQFP-32 7×7 mm 0.8 mm", width: 10, height: 10, mass: .3, pads: [...side(8, i => [-4.175, -2.8 + .8 * i, 1.5, .5]), ...side(8, i => [-2.8 + .8 * i, 4.175, .5, 1.5]), ...side(8, i => [4.175, 2.8 - .8 * i, 1.5, .5]), ...side(8, i => [2.8 - .8 * i, -4.175, .5, 1.5])] });
  await run("pcb.footprint", { name: "psop8", kind: "custom", smd: true, title: "PSOP-8 with exposed pad", width: 6.9, height: 5, mass: .15, pads: [[-2.7, -1.905, 1.55, .6], [-2.7, -.635, 1.55, .6], [-2.7, .635, 1.55, .6], [-2.7, 1.905, 1.55, .6], [2.7, 1.905, 1.55, .6], [2.7, .635, 1.55, .6], [2.7, -.635, 1.55, .6], [2.7, -1.905, 1.55, .6], [0, 0, 2.29, 3]] });
  await run("pcb.footprint", { name: "sot235", kind: "custom", smd: true, title: "SOT-23-5", width: 4.1, height: 3.1, mass: .02, pads: [[-1.1375, -.95, 1.325, .6], [-1.1375, 0, 1.325, .6], [-1.1375, .95, 1.325, .6], [1.1375, .95, 1.325, .6], [1.1375, -.95, 1.325, .6]] });
  await run("pcb.footprint", { name: "r2512", kind: "custom", smd: true, title: "2512 shunt", width: 7.4, height: 3.5, mass: .05, pads: [[-2.9625, 0, 1.225, 3.35], [2.9625, 0, 1.225, 3.35]] });
  await run("pcb.footprint", { name: "srr1260", kind: "custom", smd: true, title: "Bourns SRR1260", width: 13, height: 13, mass: 4.5, pads: [[-4.85, 0, 2.9, 5.4], [4.85, 0, 2.9, 5.4]] });
  const net = async (ref, list) => { for (const [i, n] of list.entries()) await run("pcb.net", { pad: `${ref}:${i + 1}`, net: n }); };
  const pour = (id, net, layer, x0, y0, x1, y1) => run("pcb.pour", { id, net, layer, points: [[x0, y0], [x1, y0], [x1, y1], [x0, y1]] });
  await pour("vbat", "VBAT", "top", 12.5, 1, 77, 7.35); await pour("batin", "BAT_P", "top", 1, 1, 8, 7.35); await pour("gnd", "GND", "bottom", 1, 1, 77, 41);
  let via = 0; const v = (x, y, n) => run("pcb.via", { at: [x, y], net: n, ref: `V${++via}` });
  const phases = [["A", 21], ["B", 41], ["C", 61]];
  for (const [X, c] of phases) {
    const i = "ABC".indexOf(X) + 1;
    await pour(`ph${X}`, `PH_${X}`, "top", c - 7, 7.75, c + 1.3, 17.4); await pour(`gnd${X}`, "GND", "top", c - 3, 17.75, c + 1.3, 21.4);
    await run("pcb.place", { footprint: "tdson8", at: [c, 6], ref: `Q${2 * i - 1}`, name: `${X} high side, Infineon BSC040N10NS5` }); await net(`Q${2 * i - 1}`, [`PH_${X}`, `PH_${X}`, `PH_${X}`, `GH_${X}`, "VBAT"]);
    await run("pcb.place", { footprint: "tdson8", at: [c, 16], ref: `Q${2 * i}`, name: `${X} low side, Infineon BSC040N10NS5` }); await net(`Q${2 * i}`, ["GND", "GND", "GND", `GL_${X}`, `PH_${X}`]);
    // Six vias under each low-side source: about 1.1 A each at a 10 °C rise (IPC-2221), for 7 A RMS.
    for (const [k, [dx, y]] of [[-2.4, 19.75], [-1.2, 19.75], [0, 19.75], [-2.4, 20.9], [-1.2, 20.9], [0, 20.9]].entries()) await run("pcb.via", { at: [c + dx, y], net: "GND", ref: `VS${X}${k + 1}` });
    await run("pcb.place", { footprint: "soic8", at: [c + 8.5, 12.5], ref: `U${i}`, name: `${X} gate driver, TI UCC27211` }); await net(`U${i}`, ["V12", `BOOT_${X}`, `GHD_${X}`, `PH_${X}`, `HI_${X}`, `LI_${X}`, "GND", `GLD_${X}`]);
    await run("pcb.place", { footprint: "r0603", at: [c + 4, 8.4], ref: `R${2 * i - 1}`, name: "3 Ω 0603 gate resistor" }); await net(`R${2 * i - 1}`, [`GH_${X}`, `GHD_${X}`]);
    await run("pcb.place", { footprint: "r0603", at: [c + 4, 18.4], ref: `R${2 * i}`, name: "3 Ω 0603 gate resistor" }); await net(`R${2 * i}`, [`GL_${X}`, `GLD_${X}`]);
    await run("pcb.place", { footprint: "r0603", at: [c + 3.25, 13], ref: `C${3 * i - 2}`, name: "100 nF 100 V 0603 bootstrap" }); await run("pcb.move", { id: `C${3 * i - 2}`, rotation: "90" }); await net(`C${3 * i - 2}`, [`BOOT_${X}`, `PH_${X}`]);
    await run("pcb.place", { footprint: "r0603", at: [c + 8.5, 8.75], ref: `C${3 * i - 1}`, name: "1 µF 25 V 0603 driver bypass" }); await net(`C${3 * i - 1}`, ["V12", "GND"]);
    await run("pcb.via", { at: [c + 10.25, 8.75], net: "GND", ref: `VB${X}` }); await run("pcb.via", { at: [c + 12.75, 11.75], net: "GND", ref: `VP${X}` });
    // The rail and the low-side gate cross on the bottom layer.
    await run("pcb.via", { at: [c + 6.5, 9], net: "V12", ref: `VR${X}` }); await run("pcb.via", { at: [c + 12.75, 9.75], net: `GLD_${X}`, ref: `VG${X}` }); await run("pcb.via", { at: [c + 6, 19.5], net: `GLD_${X}`, ref: `VH${X}` });
    await run("pcb.place", { footprint: "c1210", at: [c + 8.5, 4.25], ref: `C${3 * i}`, name: "4.7 µF 100 V X7S 1210 bus" }); await net(`C${3 * i}`, ["VBAT", "GND"]);
    await run("pcb.via", { at: [c + 12, 3.25], net: "GND", ref: `VC${X}1` }); await run("pcb.via", { at: [c + 12, 5.25], net: "GND", ref: `VC${X}2` });
    await run("pcb.via", { at: [c + 4.5, 13], net: `GHD_${X}`, ref: `VHO${X}` }); await run("pcb.via", { at: [c + 3.75, 10.25], net: `GHD_${X}`, ref: `VHP${X}` });
    await run("pcb.place", { footprint: "wirepad", at: [c - 4.75, 12.5], ref: `J${i + 1}`, name: `Motor phase ${X}` }); await net(`J${i + 1}`, [`PH_${X}`]);
  }
  await run("pcb.place", { footprint: "wirepad", at: [3.5, 4], ref: "J1", name: "Battery +" }); await net("J1", ["BAT_P"]);
  await run("pcb.place", { footprint: "wirepad", at: [3.75, 13], ref: "J5", name: "Battery −" }); await net("J5", ["GND"]);
  await run("pcb.place", { footprint: "r2512", at: [10, 4], ref: "RS1", name: "1 mΩ 2512 current shunt" }); await net("RS1", ["BAT_P", "VBAT"]);
  // Current sense: INA240A2, gain 50, so 50 mV/A into the ADC. Its inputs are Kelvin taps on the shunt's own pads.
  await run("pcb.place", { footprint: "soic8", at: [9.5, 13], ref: "U7", name: "Current sense, TI INA240A2" }); await net("U7", ["VBAT", "GND", "GND", "GND", "ISENSE", "3V3", "GND", "BAT_P"]);
  await run("pcb.place", { footprint: "r0603", at: [9.5, 18], ref: "C20", name: "100 nF 0603 INA240 bypass" }); await net("C20", ["GND", "3V3"]);
  const hole = (ref, x, y, n) => run("pcb.via", { at: [x, y], net: n, ref });
  await hole("VKP1", 7.5, 8.5, "BAT_P"); await hole("VKP2", 13.25, 9.25, "BAT_P"); await hole("VI1", 6.75, 16.75, "GND"); await hole("VI2", 7.25, 18, "GND"); await hole("VIS2", 10, 14.25, "ISENSE");
  // 50 V to 12 V: LM5017 with TI's Type 2 ripple feed. 604 k gives 199 kHz. UVLO is tied to VIN, which the datasheet allows.
  await run("pcb.place", { footprint: "psop8", at: [8.5, 27], ref: "U6", name: "12 V buck, TI LM5017MR" }); await net("U6", ["GND", "BAT_P", "BAT_P", "RON", "FB", "VCC", "BST", "SW", "GND"]);
  for (const [k, dx, dy] of [[1, -.6, -.8], [2, .6, -.8], [3, 0, .8]]) await hole(`VE${k}`, 8.5 + dx, 27 + dy, "GND");
  await run("pcb.place", { footprint: "c1210", at: [2.75, 26.5], ref: "C21", name: "2.2 µF 100 V X7R 1210 buck input" }); await run("pcb.move", { id: "C21", rotation: "90" }); await net("C21", ["GND", "BAT_P"]);
  await run("pcb.place", { footprint: "r0603", at: [3, 31.5], ref: "R11", name: "604 kΩ 0603 on-time" }); await net("R11", ["BAT_P", "RON"]);
  await run("pcb.place", { footprint: "r0603", at: [14, 30.75], ref: "R14", name: "1.00 kΩ 0603 feedback bottom" }); await net("R14", ["GND", "FB"]);
  await run("pcb.place", { footprint: "r0603", at: [14, 32.75], ref: "R15", name: "8.87 kΩ 0603 feedback top" }); await net("R15", ["V12", "FB"]);
  await run("pcb.place", { footprint: "r0603", at: [14, 34.75], ref: "C22", name: "33 nF 0603 ripple feed" }); await net("C22", ["V12", "FB"]);
  await run("pcb.place", { footprint: "r0603", at: [13.5, 28.25], ref: "C23", name: "1 µF 25 V 0603 VCC" }); await net("C23", ["VCC", "GND"]);
  await run("pcb.place", { footprint: "r0603", at: [13.5, 25.25], ref: "C24", name: "10 nF 0603 bootstrap" }); await net("C24", ["BST", "SW"]);
  await hole("VF1", 15.75, 26.75, "GND"); await hole("VF2", 12, 30.75, "GND");
  await run("pcb.place", { footprint: "srr1260", at: [22, 32.5], ref: "L1", name: "220 µH Bourns SRR1260-221K" }); await net("L1", ["SW", "V12"]);
  await run("pcb.place", { footprint: "c1210", at: [31.75, 24.75], ref: "C25", name: "22 µF 25 V X7R 1210 output" }); await net("C25", ["V12", "CESR"]);
  await run("pcb.place", { footprint: "r0603", at: [36, 24.75], ref: "R16", name: "0.33 Ω 0603 ripple resistor" }); await net("R16", ["CESR", "GND"]);
  await hole("VC12", 30.25, 23.6, "V12"); await hole("VF4", 37.5, 26, "GND");
  // 12 V to 3.3 V.
  await run("pcb.place", { footprint: "sot235", at: [33.5, 31.25], ref: "U5", name: "3.3 V regulator, TI LP2985-33" }); await net("U5", ["V12", "GND", "V12", "BP", "3V3"]);
  await run("pcb.place", { footprint: "r0603", at: [31, 34.5], ref: "C26", name: "1 µF 25 V 0603 LDO input" }); await net("C26", ["V12", "GND"]);
  await run("pcb.place", { footprint: "r0603", at: [37.5, 29], ref: "C27", name: "2.2 µF 10 V 0603 LDO output" }); await net("C27", ["3V3", "GND"]);
  await run("pcb.place", { footprint: "r0603", at: [37.5, 32.25], ref: "C28", name: "10 nF 0603 LDO noise" }); await net("C28", ["BP", "GND"]);
  await hole("VF5", 39.75, 29, "GND"); await hole("VF6", 39.75, 32.25, "GND"); await hole("VF7", 31, 30.5, "GND"); await hole("VF8", 33.25, 34.5, "GND"); await hole("V33A", 38.75, 27, "3V3");
  // The controller, wired as AM32's F0 group B: TIM1 on PA8–10 and PB1, PB0, PA7, comparators on PA0, PA4, PA5 against PA1,
  // input on PB4, telemetry on PA14, current on PA3, voltage on PA6.
  const mcu = ["3V3", "NC", "NC", "NRST", "3V3", "BEMF_A", "BEMF_N", "NC", "ISENSE", "BEMF_B", "BEMF_C", "VSENSE", "LI_C", "LI_B", "LI_A", "GND", "3V3", "HI_C", "HI_B", "HI_A", "NC", "NC", "SWDIO", "SWCLK", "NC", "NC", "INPUT", "NC", "NC", "NC", "GND", "GND"];
  await run("pcb.place", { footprint: "lqfp32", at: [56, 31], ref: "U4", name: "Controller, ST STM32F051K8T6, AM32" }); await run("pcb.move", { id: "U4", rotation: "270" }); await net("U4", mcu);
  await run("pcb.place", { footprint: "header8", at: [39.5, 37], ref: "J6", name: "Signal, telemetry and SWD" }); await net("J6", ["NC", "NC", "3V3", "NRST", "SWDIO", "SWCLK", "GND", "INPUT_IN"]);
  await run("pcb.place", { footprint: "r0603", at: [47.5, 33.5], ref: "R17", name: "100 Ω 0603 signal input" }); await net("R17", ["INPUT_IN", "INPUT"]);
  for (const [ref, at, nets, name] of [["C29", [61.25, 25], ["3V3", "GND"], "100 nF 0603 VDD"], ["C31", [56, 40], ["3V3", "GND"], "1 µF 0603 VDDA"]]) { await run("pcb.place", { footprint: "r0603", at, ref, name }); await net(ref, nets); }
  await run("pcb.place", { footprint: "r0603", at: [52.5, 37.5], ref: "C30", name: "100 nF 0603 VDD" }); await run("pcb.move", { id: "C30", rotation: "90" }); await net("C30", ["3V3", "GND"]);
  await hole("VM1", 61.75, 27, "GND"); await hole("VM2", 49.75, 34.75, "GND"); await hole("VM3", 51, 38.25, "GND"); await hole("VM6", 58.25, 40, "GND"); await hole("V33E", 56, 38.25, "3V3"); await hole("VM5", 63.25, 25, "GND"); await hole("V33B", 51, 35.25, "3V3"); await hole("V33D", 58, 28.5, "3V3");
  await hole("VNR", 54.75, 37.25, "NRST"); await hole("VMA", 59.25, 38.75, "BEMF_A");
  await hole("VLA1", 33.75, 13, "LI_A"); await hole("VLA2", 62.5, 28.75, "LI_A"); await hole("VLB1", 53.75, 13, "LI_B"); await hole("VLB2", 64.5, 29.75, "LI_B");
  await hole("VIS", 61.75, 35.5, "ISENSE"); await hole("VMC", 62.5, 32.25, "BEMF_C"); await hole("VMB", 62.75, 34.25, "BEMF_B");
  // Sensing: 16:1 on each phase and 18:1 on the bus, so 50.4 V stays below the 3.3 V rail.
  for (const [X, c] of phases) {
    const i = "ABC".indexOf(X);
    await run("pcb.place", { footprint: "r0603", at: [c - 5.5, 18.25], ref: `R${30 + 3 * i}`, name: "150 kΩ 0603 phase sense" }); await run("pcb.move", { id: `R${30 + 3 * i}`, rotation: "90" }); await net(`R${30 + 3 * i}`, [`PH_${X}`, `BEMF_${X}`]);
    await run("pcb.place", { footprint: "r0603", at: [X === "A" ? 11.75 : c - 8.25, 19.5], ref: `R${31 + 3 * i}`, name: "10 kΩ 0603 phase sense" }); await net(`R${31 + 3 * i}`, ["GND", `BEMF_${X}`]);
    await run("pcb.place", { footprint: "r0603", at: X === "A" ? [58, 37.75] : [66.25, 38.25 - 1.75 * (i + 1)], ref: `R${32 + 3 * i}`, name: "10 kΩ 0603 neutral" }); if (X === "A") await run("pcb.move", { id: "R32", rotation: "90" }); await net(`R${32 + 3 * i}`, X === "A" ? ["BEMF_N", "BEMF_A"] : [`BEMF_${X}`, "BEMF_N"]);
    await hole(`VN${X}`, X === "A" ? 9.5 : c - 10.5, 19.5, "GND"); if (X === "A") await hole("VBS2", 13.25, 21, "BEMF_A"); await hole(`VBS${X}`, c - 5.5, { A: 21.5, B: 20.25, C: 22.75 }[X], `BEMF_${X}`);
  }
  await run("pcb.place", { footprint: "r0603", at: [75.5, 8], ref: "R40", name: "170 kΩ 0603 bus sense" }); await run("pcb.move", { id: "R40", rotation: "90" }); await net("R40", ["VBAT", "VSENSE"]);
  await run("pcb.place", { footprint: "r0603", at: [70.5, 31.5], ref: "R41", name: "10 kΩ 0603 bus sense" }); await net("R41", ["VSENSE", "GND"]);
  await hole("VF9", 72.5, 31.5, "GND"); await hole("VV1", 1.75, 7, "BAT_P"); await hole("VV4", 1.25, 28.75, "BAT_P");
  for (const [n, a] of [["VBAT", 12], ["BAT_P", 12], ["GND", 12], ["V12", .2], ["SW", .3], ["3V3", .05], ...["A", "B", "C"].flatMap(X => [[`PH_${X}`, 11], [`GH_${X}`, .3], [`GL_${X}`, .3], [`GHD_${X}`, .3], [`GLD_${X}`, .3], [`HI_${X}`, .01], [`LI_${X}`, .01], [`BOOT_${X}`, .1]])]) await run("pcb.current", { net: n, amps: a });
  for (const X of "ABC") await run("pcb.load_net", { net: `PH_${X}` });
  const trace = (from, to, layer = "top", width = .3) => run("pcb.route", { from, to, layer, width });
  const attempt = async (...a) => { try { await trace(...a); } catch (e) { fails.push(`${a[0]}–${a[1]} ${a[2] ?? "top"}: ${e.message.slice(0, 44)}`); } }, fails = [];
  // Gates, rail and the tight loops first: everything else can go the long way round.
  for (const [X, u, i] of [["A", "U1", 1], ["B", "U2", 2], ["C", "U3", 3]]) { await attempt(`${u}:8`, `VG${X}:1`); await attempt(`VG${X}:1`, `VH${X}:1`, "bottom", .4); await attempt(`VH${X}:1`, `R${2 * i}:2`); await attempt(`${u}:1`, `VR${X}:1`); await attempt(`${u}:1`, `C${3 * i - 1}:1`); await attempt(`C${3 * i - 1}:2`, `VB${X}:1`); await attempt(`${u}:7`, `VP${X}:1`); await attempt(`C${3 * i}:2`, `VC${X}1:1`, "top", .8); await attempt(`C${3 * i}:2`, `VC${X}2:1`, "top", .8); await attempt(`${u}:2`, `C${3 * i - 2}:1`); await attempt(`${u}:3`, `VHO${X}:1`, "top", .25); await attempt(`VHO${X}:1`, `VHP${X}:1`, "bottom", .3); await attempt(`VHP${X}:1`, `R${2 * i - 1}:2`, "top", .25); }
  for (const [a, x] of [["VRA:1", 27.5], ["VRB:1", 47.5]]) { await run("pcb.trace", { from: a, to: a === "VRA:1" ? "VRB:1" : "VRC:1", layer: "bottom", width: .4, points: [[x, 9], [x + 1.25, 7.75], [x + 18.75, 7.75], [x + 20, 9]] }); }
  // Buck: the switching loop, short and wide.
  await attempt("U6:8", "L1:1", "top", .5); await attempt("U6:7", "C24:1"); await attempt("U6:8", "C24:2"); await attempt("U6:6", "C23:1"); await attempt("U6:1", "C21:1"); await attempt("U6:2", "C21:2", "top", .4);
  await attempt("RS1:1", "VKP1:1", "top", .25); await attempt("VKP1:1", "VKP2:1", "bottom", .25); await attempt("VKP2:1", "U7:8", "top", .25); await attempt("RS1:2", "U7:1", "top", .25);
  const by = (from, to, width, points, layer = "bottom") => run("pcb.trace", { from, to, layer, width, points });
  const by2 = async (...a) => { try { await by(...a); } catch (e) { fails.push(`${a[0]}–${a[1]} ${e.message.slice(0, 80)}`); } };
  // Sense taps, before logic takes the room.
  for (const X of "ABC") { const i = "ABC".indexOf(X); await attempt(`R${31 + 3 * i}:1`, `VN${X}:1`); await attempt(`R${30 + 3 * i}:2`, X === "A" ? "VBSA:1" : `R${31 + 3 * i}:2`, "top", .25); if (X !== "A") await attempt(`R${30 + 3 * i}:2`, `VBS${X}:1`, "top", .25); }
  await attempt("R31:2", "VBS2:1", "top", .25); await attempt("VBS2:1", "VBSA:1", "bottom", .25);
  // 3.3 V runs above C25 to the current-sense amplifier; then the high-side gates, which take the direct way.
  await attempt("U5:5", "C27:1", "top", .3); await attempt("C27:1", "V33A:1", "top", .3); await by2("V33A:1", "U7:6", .3, [[38.75, 27], [38.75, 22.75], [14.25, 22.75], [14, 22.5], [14, 16], [13.5, 15.5], [13.5, 14.5], [12.2, 13.635]], "top");
  await attempt("U7:6", "C20:2", "top", .25); await attempt("C20:1", "VI2:1"); await attempt("U7:5", "VIS2:1", "top", .2); await attempt("C27:1", "J6:3", "top", .3);
  await attempt("VLA2:1", "U4:15", "top", .2); await attempt("VLB2:1", "U4:14", "top", .2);
  await by2("U1:5", "U4:20", .25, [[32.2, 14.405], [32.25, 15.25], [31.75, 15.75], [29.25, 15.75], [28.75, 16.25], [28.75, 20.75], [29.75, 21.75], [39.5, 21.75], [42.5, 24.75], [55.25, 24.75], [56.4, 25.9], [56.4, 26.825]], "top");
  await attempt("U2:5", "U4:19", "top", .25); await attempt("U3:5", "U4:18", "top", .25); await by2("R41:1", "U4:12", .2, [[69.7, 31.5], [61.2, 31.5], [60.175, 31.4]], "top");
  await attempt("U4:24", "J6:6", "top", .25); await attempt("U4:23", "J6:5", "top", .25); await attempt("U3:6", "U4:13", "top", .2);
  // Bottom-layer rows, in an order that nests, with a gap in the ground plane at x = 36.
  await attempt("U1:6", "VLA1:1", "top", .25); await by2("VLA1:1", "VLA2:1", .25, [[33.75, 13], [37.2, 16.45], [37.2, 21.9], [61.5, 21.9], [62.5, 22.9], [62.5, 28.75]]);
  await attempt("U2:6", "VLB1:1", "top", .25); await by2("VLB1:1", "VLB2:1", .25, [[53.75, 13], [55.75, 15], [63.5, 15], [64.5, 16], [64.5, 29.75]]);
  await by2("VBSC:1", "VMC:1", .25, [[55.5, 22.75], [60.75, 22.75], [60.75, 30.5], [62.5, 32.25]]);
  await by2("VBSB:1", "VMB:1", .25, [[35.5, 20.25], [35.5, 22.6], [36.5, 23.6], [60.1, 23.6], [60.1, 34.25], [62.75, 34.25]]);
  await by2("VIS2:1", "VIS:1", .25, [[10, 14.25], [13.5, 17.75], [17, 21.25], [17.75, 22], [20.15, 24.4], [59.45, 24.4], [59.45, 35.5], [61.75, 35.5]]); await attempt("U4:9", "VIS:1", "top", .2);
  await by2("VBSA:1", "VMA:1", .25, [[15.5, 21.5], [15.5, 24.3], [16.3, 25.1], [58.25, 25.1], [58.75, 25.6], [58.75, 37.25], [59.25, 37.75], [59.25, 38.75]]);
  await attempt("U4:11", "VMC:1", "top", .2); await attempt("U4:10", "VMB:1", "top", .2);
  await attempt("VMC:1", "R38:1", "top", .25); await attempt("VMB:1", "R35:1", "top", .25);
  await attempt("U4:7", "R32:1", "top", .2); await attempt("U4:6", "R32:2", "top", .2); await attempt("VMA:1", "R32:2", "top", .25);
  await attempt("R38:2", "R35:2", "top", .25); await attempt("R35:2", "R32:1", "top", .25);
  await attempt("R40:2", "R41:1", "top", .25); await attempt("R41:2", "VF9:1");
  await attempt("U4:4", "VNR:1", "top", .2); await attempt("VNR:1", "J6:4", "bottom", .25);
  await attempt("U4:5", "V33E:1", "top", .2); await attempt("V33E:1", "C31:1", "top", .25); await attempt("C31:2", "VM6:1", "top", .25);
  // Power rails.
  await by2("VV1:1", "VV4:1", .4, [[1.75, 7], [1.25, 7.5], [1.25, 28.75]]); await attempt("VV4:1", "C21:2", "top", .3); await attempt("C21:2", "U6:3", "top", .3); await attempt("U6:3", "U6:2", "top", .3); await attempt("C21:2", "R11:1", "top", .25); await attempt("U6:4", "R11:2", "top", .25);
  await attempt("L1:2", "C25:1", "top", .5); await attempt("C25:1", "VC12:1", "top", .5); await attempt("VC12:1", "VRA:1", "bottom", .4); await attempt("L1:2", "U5:1", "top", .3); await attempt("U5:1", "U5:3", "top", .25); await attempt("U5:1", "C26:1", "top", .25); await attempt("U6:5", "R14:2", "top", .25); await attempt("R14:2", "R15:2", "top", .25); await attempt("R15:2", "C22:2", "top", .25); await attempt("L1:2", "C22:1", "top", .25); await attempt("C22:1", "R15:1", "top", .25);
  await attempt("V33A:1", "V33B:1", "bottom", .3); await attempt("V33B:1", "U4:1", "top", .25); await attempt("U4:1", "C30:1", "top", .25); await attempt("V33A:1", "V33D:1", "bottom", .3); await attempt("V33D:1", "U4:17", "top", .25); await attempt("V33D:1", "C29:1", "top", .25); await attempt("V33B:1", "V33E:1", "bottom", .3);
  await attempt("J6:8", "R17:1", "top", .25); await attempt("R17:2", "U4:27", "top", .2);
  if (fails.length) console.log("fails", fails.join(" | "));
  const r = await run("pcb.route_all", { width: .3 });
  return { p, r };
}
if (process.argv[2] === "run") {
  const { p, r } = await escBoard(); console.log(r.message);
  const { pcbChecks } = await import(new URL("../../web/cad/studio-model.js", import.meta.url)), { densityIssues } = await import(new URL("../../web/cad/domains/pcb.js", import.meta.url));
  const c = pcbChecks(p.store.docs.pcb); console.log("issues", c.issues.length, "unrouted", c.unrouted, "routes", p.store.docs.pcb.routes.length); const codes = {}; for (const i of c.issues) codes[i.code] = (codes[i.code] || 0) + 1; console.log(codes); for (const i of c.issues.slice(0, 14)) console.log(" ", i.message.slice(0, 110)); for (const n of c.nets.filter(n => n.remaining)) console.log(" open:", n.name, n.remaining); console.log("density", densityIssues(p.store.docs.pcb).map(i => i.message.slice(0, 90)));
}

if (process.argv[2] === "file") { const { p } = await escBoard(); const { writeFileSync } = await import("node:fs"); writeFileSync(process.argv[3], JSON.stringify(p.store.docs.pcb)); }

if (process.argv[2] === "islands") {
  const { p } = await escBoard(), d = p.store.docs.pcb, { pads, pourIslands, inPolygon, pcbChecks } = await import(new URL("../../web/cad/studio-model.js", import.meta.url)), all = pads(d), net = r => all.find(q => q.id === r.from).net;
  for (const q of d.pours) {
    const { under, islands } = pourIslands(d, q, all), by = new Map(); for (const a of all) if (a.net === q.net && (!a.smd || q.layer === "top")) for (const k of under(a)) { if (!by.has(k)) by.set(k, []); by.get(k).push(a.id); }
    if (by.size > 1) { console.log("POUR", q.id, q.net, "islands", islands, [...by.entries()].map(([k, v]) => `#${k}(${v.length}): ${v.slice(0, 7).join(" ")}`).join(" || "));
      for (const r of d.routes) if (r.layer === q.layer && net(r) !== q.net && r.points.some(v => inPolygon(q.points, v[0], v[1]))) console.log("   crossed by", net(r), r.from, "->", r.to, "w", r.width, JSON.stringify(r.points.slice(0, 8))); }
  }
  for (const r of d.routes) if (["VBAT", "BAT_P", "PH_A", "PH_B", "PH_C"].includes(net(r)) && r.layer === "bottom") console.log("bottom power route", net(r), r.from, r.to);
  const c = pcbChecks(d); for (const n of c.nets) for (const [a, b] of n.airwires) console.log("open", n.name, a.id, b.id, `(${a.x.toFixed(1)},${a.y.toFixed(1)})-(${b.x.toFixed(1)},${b.y.toFixed(1)})`);
  for (const i of c.issues) console.log(i.code, i.message.slice(0, 110));
}

if (process.argv[2] === "paths") {
  const { p } = await escBoard(), d = p.store.docs.pcb, { pads } = await import(new URL("../../web/cad/studio-model.js", import.meta.url)), all = pads(d), net = r => all.find(q => q.id === r.from).net;
  for (const r of d.routes) if (r.layer === "top" && ["HI_B", "HI_C", "LI_C", "SWDIO", "SWCLK", "BEMF_C", "VSENSE"].includes(net(r))) console.log(net(r), r.from, r.to, JSON.stringify(r.points.map(q => q.map(v => +v.toFixed(2)))));
}
