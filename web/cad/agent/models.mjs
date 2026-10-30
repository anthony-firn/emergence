// The model worker's protocol without a Worker: the same WASM exports, called
// directly, for the CLI, the MCP server and the Node tests.
import { readFile } from "node:fs/promises";
import { existsSync } from "node:fs";
import { fileURLToPath } from "node:url";
const here = p => fileURLToPath(new URL(p, import.meta.url));
export const WASM_PATHS = [here("../../../fluid/target/wasm32-unknown-unknown/web/emergence_web.wasm"), here("../../dist/emergence_web.wasm")];
export async function nodeModels(path = WASM_PATHS.find(existsSync)) {
  if (!path) throw new Error("emergence_web.wasm not found. Build it: cd fluid && cargo build --profile web -p emergence-web --target wasm32-unknown-unknown");
  const { instance } = await WebAssembly.instantiate(await readFile(path), {}), api = instance.exports;
  const json = () => JSON.parse(new TextDecoder().decode(new Uint8Array(api.memory.buffer, api.result_ptr(), api.result_len())));
  const send = input => {
    const bytes = new TextEncoder().encode(JSON.stringify(input)), ptr = api.workspace_input(bytes.length);
    if (!ptr) throw new Error("Model input exceeds 1 MB");
    new Uint8Array(api.memory.buffer, ptr, bytes.length).set(bytes);
  };
  const modes = { workspace: 0, mesh: 1, model: 2, pcb: 3, line: 4, chip: 5 };
  async function request(type, args = {}) {
    const start = performance.now();
    if (!(type in modes)) throw new Error("Unknown model request");
    send(type === "workspace" || type === "mesh" ? { params: args.params || args, features: args.features || [], target: args.target || "sink" } : args.input || {});
    api.workspace_run(modes[type]);
    const result = json();
    if (result.error) throw new Error(result.error);
    if (type === "model" && args.input?.binary) {
      result.positions = new Float32Array(api.memory.buffer, api.mesh_vertices_ptr(), api.mesh_vertices_len()).slice();
      result.indices = new Uint32Array(api.memory.buffer, api.mesh_indices_ptr(), api.mesh_indices_len()).slice();
    }
    return { result, ms: performance.now() - start };
  }
  return { request, path };
}
export const loadCorners = async () => JSON.parse(await readFile(here("../data/chip-corners.json"), "utf8"));
export const loadVerilog = async () => ({ mini8: await readFile(here("../data/mini8.v"), "utf8"), "tt-motor": await readFile(here("../data/tt_motor_controller.v"), "utf8") });
