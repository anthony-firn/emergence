import { writeFileSync } from "node:fs";
const { kicadPCB } = await import(new URL("../../web/cad/domains/pcb.js", import.meta.url));
const { exampleFile } = await import(new URL("../../web/cad/agent/headless.mjs", import.meta.url));
const d = (await exampleFile("drone")).documents.pcb;
writeFileSync(process.argv[2], kicadPCB(d)); console.log("written", d.components.length, "components", d.routes.length, "routes", d.pours.length, "pours"); process.exit(0);
