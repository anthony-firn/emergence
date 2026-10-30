// A project's guided tour: each step may run a real command against the live
// model and report how the outputs and checks it names moved.
import { format, label } from "../core/units.js";
export function observe(step, after, before) {
  const out = (step.watch || []).map(port => {
    const a = after?.outputs[port], b = before?.outputs[port], p = after?.ports?.outputs.get(port);
    return `${p?.label || port} ${format(b?.value, b?.unit ?? a?.unit, 4)} → ${format(a?.value, a?.unit ?? b?.unit, 4)}`;
  });
  for (const ref of step.checks || []) {
    const [domain, id] = ref.split("/"), c = after?.checks.find(x => x.domain === domain && x.id === id);
    out.push(`${c?.text || ref}: ${c?.status ?? "not evaluated"}`);
  }
  return out.length ? out.join(" · ") + "." : "Done.";
}
export { label };
