let api;
const ready = (async () => {
  const response = await fetch("emergence_web.wasm");
  if (!response.ok) throw new Error(`Could not load models (${response.status})`);
  // Byte instantiation also supports static hosts with generic WASM MIME types.
  const { instance } = await WebAssembly.instantiate(await response.arrayBuffer(), {});
  api = instance.exports;
})();
const json = () => JSON.parse(new TextDecoder().decode(new Uint8Array(api.memory.buffer, api.result_ptr(), api.result_len())));
self.onmessage = async ({ data }) => {
  try {
    await ready;
    let result;
    const start = performance.now();
    if (data.type === "optics") { api.optics_solve(data.material, data.x, data.tolerance); result = json(); }
    else if (data.type === "design") { api.design_solve(data.fins, data.velocity); result = json(); }
    else if (data.type === "workspace" || data.type === "mesh") {
      const bytes = new TextEncoder().encode(JSON.stringify({params:data.params || data, features:data.features || [], target:data.target || "sink"}));
      const ptr = api.workspace_input(bytes.length);
      if (!ptr) throw new Error("Model input exceeds 64 KB");
      new Uint8Array(api.memory.buffer,ptr,bytes.length).set(bytes);
      api.workspace_run(data.type === "mesh" ? 1 : 0); result = json();
    }
    else if (["model", "pcb", "line", "chip"].includes(data.type)) {
      // Feature-tree meshes come back as typed arrays, not JSON numbers.
      const bytes = new TextEncoder().encode(JSON.stringify(data.input || {}));
      const ptr = api.workspace_input(bytes.length);
      if (!ptr) throw new Error("Model input exceeds 1 MB");
      new Uint8Array(api.memory.buffer,ptr,bytes.length).set(bytes);
      api.workspace_run({ model: 2, pcb: 3, line: 4, chip: 5 }[data.type]); result = json();
      if (data.type === "model" && !result.error && data.input?.binary) {
        result.positions = new Float32Array(api.memory.buffer, api.mesh_vertices_ptr(), api.mesh_vertices_len()).slice();
        result.indices = new Uint32Array(api.memory.buffer, api.mesh_indices_ptr(), api.mesh_indices_len()).slice();
      }
    }
    else if (["flow-init", "flow-step", "flow-frame"].includes(data.type)) {
      if (data.type === "flow-init" && !api.flow_init(data.n, data.exact, data.re)) throw new Error("Invalid flow settings");
      if (data.type === "flow-step") api.flow_step();
      const pointer = api.flow_frame(data.z);
      if (!pointer) throw new Error("Initialize a flow first");
      const field = new Float32Array(api.memory.buffer, pointer, data.n * data.n * 4).slice();
      api.flow_report(); result = { ...json(), field };
    } else if (data.type !== "ready") throw new Error("Unknown model request");
    const transfer = [result?.field?.buffer, result?.positions?.buffer, result?.indices?.buffer].filter(Boolean);
    self.postMessage({ id: data.id, result, ms: performance.now() - start }, transfer);
  } catch (error) { self.postMessage({ id: data.id, error: error.message || String(error) }); }
};
