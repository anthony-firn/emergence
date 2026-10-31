// Reads AM32's target table as text and KiCad's STM32F051K8Tx pin numbers, and checks each
// pin the firmware drives against the net the board puts on it. Nothing downloaded is executed.
import { readFileSync, writeFileSync } from "node:fs";
const S = (process.env.EVIDENCE_WORK ?? "/tmp/emergence-evidence") + "/";
const { pads } = await import(new URL("../../web/cad/studio-model.js", import.meta.url));
const { exampleFile } = await import(new URL("../../web/cad/agent/headless.mjs", import.meta.url));
const h = readFileSync(process.env.AM32_TARGETS ?? S + "AM32/Inc/targets.h", "utf8"), start = h.indexOf("#ifdef HARDWARE_GROUP_F0_B"), group = h.slice(start, h.indexOf("#endif", start));
const ref = h.slice(h.indexOf("#ifdef AM32REF_F051"), h.indexOf("#endif", h.indexOf("#ifdef AM32REF_F051")));
const def = (t, k) => t.match(new RegExp(`#define\\s+${k}\\s+(\\S+)`))?.[1];
const gpio = (pin, port) => `P${def(group, port).replace("GPIO", "")}${def(group, pin).replace("LL_GPIO_PIN_", "")}`;
const comp = k => def(group, k).replace("COMP_", "");
// STM32F051K8Tx LQFP-32, from KiCad's MCU_ST_STM32F0 symbol.
const lqfp = "VDD PF0 PF1 NRST VDDA PA0 PA1 PA2 PA3 PA4 PA5 PA6 PA7 PB0 PB1 VSS VDD PA8 PA9 PA10 PA11 PA12 PA13 PA14 PA15 PB3 PB4 PB5 PB6 PB7 BOOT0 VSS".split(" ");
const pinOf = name => lqfp.indexOf(name) + 1;
const expected = [
  ["signal input", gpio("INPUT_PIN", "INPUT_PIN_PORT"), "INPUT"],
  ...["A", "B", "C"].flatMap(X => [[`phase ${X} high side`, gpio(`PHASE_${X}_GPIO_HIGH`, `PHASE_${X}_GPIO_PORT_HIGH`), `HI_${X}`], [`phase ${X} low side`, gpio(`PHASE_${X}_GPIO_LOW`, `PHASE_${X}_GPIO_PORT_LOW`), `LI_${X}`], [`phase ${X} comparator`, comp(`PHASE_${X}_COMP`), `BEMF_${X}`]]),
  ["comparator neutral (COMP1_INP)", "PA1", "BEMF_N"],
  ["bus voltage ADC", `PA${def(ref, "VOLTAGE_ADC_PIN").replace("LL_GPIO_PIN_", "")}`, "VSENSE"], ["current ADC", `PA${def(ref, "CURRENT_ADC_PIN").replace("LL_GPIO_PIN_", "")}`, "ISENSE"],
  ["SWDIO", "PA13", "SWDIO"], ["SWCLK and telemetry (USART on PA14)", "PA14", "SWCLK"], ["reset", "NRST", "NRST"], ["BOOT0 low", "BOOT0", "GND"], ["VDD", "VDD", "3V3"], ["VDDA", "VDDA", "3V3"], ["VSS", "VSS", "GND"],
];
const d = (await exampleFile("drone")).documents.pcb, nets = Object.fromEntries(pads(d).filter(p => p.component === "U4").map(p => [p.number, p.net]));
let bad = 0;
for (const [what, name, net] of expected) for (let n = 1; n <= 32; n++) if (lqfp[n - 1] === name) { const ok = nets[n] === net; if (!ok) bad++; console.log(ok ? "ok  " : "BAD ", String(n).padStart(2), name.padEnd(6), (nets[n] ?? "?").padEnd(8), what); }
// The analog ranges the firmware will read.
const divider = 180 / 10, adc = 50.4 / divider, phase = 50.4 * 10 / 160;
console.log(`bus at full charge: ${adc.toFixed(2)} V at the ADC; phase: ${phase.toFixed(2)} V at the comparator; current at 15 A: ${(15 * .05).toFixed(2)} V`);
console.log("AM32REF divider would give", (50.4 / 11).toFixed(2), "V on a 3.3 V pin");
writeFileSync(S + "px4-run/firmware-pins.json", JSON.stringify({ checked: expected.length, mismatches: bad })); console.log("mismatches", bad); process.exit(0);
