import { Repository, FIELDS, values, validParams, validFeatures, featureValues, drawingValues, validDrawings, canonical } from "./repository.js";
import { Peers } from "./peers.js";
import { Scene, assemblyParts } from "./scene.js";
import { drawingDefaults, drawingSVG, profileDXF, bomCSV, releaseProblems } from "./drawings.js";

const $ = id => document.getElementById(id);
const el = (tag, text, cls) => { const e = document.createElement(tag); if (text !== undefined) e.textContent = text; if (cls) e.className = cls; return e; };
const STORE = "emergence-cad-project-v1", VIEW = "emergence-cad-view-v1";
const domains = {
  Mechanical: { icon: "⬡", color: "var(--cyan)", subtitle: "Geometry & airflow", description: "Dimensions drive the geometry, CAD export, drawing and thermal channel model together." },
  Electrical: { icon: "ϟ", color: "var(--ember)", subtitle: "Power & losses", description: "Current changes conduction losses and the MOSFET temperature feedback loop." },
  Thermal: { icon: "≋", color: "var(--red)", subtitle: "Heat & requirements", description: "Inlet temperature feeds the shared electrical and thermal design checks." },
  Photonics: { icon: "◈", color: "var(--violet)", subtitle: "Optical heat & scattering", description: "Prescribed incident power × absorbed fraction deposits heat at Q1. The separate 550 nm Mie sphere model does not determine this fraction." },
  Radiation: { icon: "☼", color: "var(--green)", subtitle: "Prescribed heat boundary", description: "Share extra watts deposited at Q1. This is a boundary input, not a nuclear or radiation transport simulation." },
};
let toastTimer;
function toast(message) { $("toast").textContent = message; $("toast").hidden = false; clearTimeout(toastTimer); toastTimer = setTimeout(() => $("toast").hidden = true, 6500); }
const fail = e => toast(e.message || String(e));
const action = fn => async event => { try { await fn(event); } catch (e) { fail(e); } };
const repo = await new Repository().init();
let storageOK = true;
try { const stored = localStorage.getItem(STORE); if (stored) await repo.ingest(JSON.parse(stored)); }
catch (e) { storageOK = false; toast(`Saved project could not be read: ${e.message}. Export any recoverable data before clearing browser storage.`); }
let active = repo.heads("main")[0], params, features, drawings, domain = "Mechanical", view = "workspace", selected = "sink", drawingPart = "sink";
const hiddenParts = new Set(), undoStack = [], redoStack = [];
try {
  const saved = JSON.parse(sessionStorage.getItem(VIEW) || "null");
  if (repo.commits.has(saved?.head)) { active = saved.head; if (validParams(saved.params)) params = values(saved.params); if (validFeatures(saved.features)) features = saved.features;
    if (saved.drawings && validDrawings(Object.fromEntries(Object.entries(saved.drawings).map(([k,d])=>[k,{...d,status:"Draft"}])))) drawings=saved.drawings; }
  $("author").value = localStorage.getItem("emergence-cad-author") || "Guest engineer";
} catch { /* Browser storage may be unavailable; export remains usable. */ }
params ||= values(repo.commits.get(active).params); features ||= structuredClone(featureValues(repo.commits.get(active))); drawings ||= structuredClone(drawingValues(repo.commits.get(active)));
let peers, mergePlan, comparisonGeneration = 0, generation = 0, evaluateTimer;
const identity = () => ({ author: $("author").value.trim() || "Guest engineer", domain });
const current = () => repo.commits.get(active);
const dirty = () => repo.diff(current().params, params).length > 0 || repo.featureDiff(featureValues(current()), features).length > 0 || canonical(drawingValues(current())) !== canonical(drawings);
function saveDraft() {
  try { sessionStorage.setItem(VIEW, JSON.stringify({ head: active, params, features, drawings })); }
  catch { storageOK = false; }
  renderStorage();
}
function renderStorage() { $("storage-status").textContent = storageOK ? "Saved on this device" : "Storage unavailable · export to keep changes"; }
function persist(publish = true) {
  if (storageOK) try { localStorage.setItem(STORE, JSON.stringify(repo.bundle())); } catch { storageOK = false; }
  saveDraft(); if (publish) peers?.publish();
}
function requireClean() { if (dirty()) throw new Error("Save or discard your draft before switching, branching or merging."); }
function checkout(id) {
  requireClean(); active = id; params = values(current().params); features=structuredClone(featureValues(current()));drawings=structuredClone(drawingValues(current())); undoStack.length=0;redoStack.length=0; $("message").value = ""; saveDraft(); renderRepository(); renderParameters(); renderTree(); queueEvaluation();
}
const snapshot=()=>structuredClone({params,features,drawings});
function remember(){undoStack.push(snapshot());if(undoStack.length>50)undoStack.shift();redoStack.length=0;}
function reviewDrawings(id){for(const [key,d] of Object.entries(drawings))if((!id||key===id)&&d.status==="Released")drawings[key]={...d,status:"For review"};}
function restore(s){({params,features,drawings}=structuredClone(s));saveDraft();renderTree();renderParameters();renderDraft();queueEvaluation();if(view==="drawings")renderDrawing();}
$("undo").addEventListener("click",()=>{if(undoStack.length){redoStack.push(snapshot());restore(undoStack.pop());}});
$("redo").addEventListener("click",()=>{if(redoStack.length){undoStack.push(snapshot());restore(redoStack.pop());}});

// Real Rust models run in the existing lab Worker, under the same site CSP.
let requestId = 0, worker;
const pending = new Map();
function workerFailure(error) {
  worker?.terminate(); worker = null;
  for (const p of pending.values()) { clearTimeout(p.timer); p.reject(error); } pending.clear();
  $("model-status").textContent = `Models unavailable: ${error.message}. Reload to retry.`;
}
try {
  worker = new Worker(new URL("../worker.js", import.meta.url), { type: "module" });
  worker.onmessage = ({ data }) => { const p = pending.get(data.id); if (!p) return; clearTimeout(p.timer); pending.delete(data.id);
    data.error || data.result?.error ? p.reject(new Error(data.error || data.result.error)) : p.resolve(data); };
  worker.onerror = e => workerFailure(new Error(e.message || "The model Worker failed"));
} catch (e) { workerFailure(e); }
function request(type, args = {}) {
  if (!worker) return Promise.reject(new Error("The model Worker is unavailable"));
  return new Promise((resolve, reject) => {
    const id = ++requestId, timer = setTimeout(() => workerFailure(new Error("Model exceeded the 45-second time limit")), 45000);
    pending.set(id, { resolve, reject, timer }); worker.postMessage({ id, type, ...args });
  });
}
const scene = new Scene($("assembly"));
function draw() { scene.update(params, domain, $("exploded").checked,features,selected,hiddenParts); $("geometry-label").textContent = `${params.width} × ${params.length} mm · ${params.fins} fins · ${params.finHeight} mm height`; $("stl").disabled=["board","q1","optic"].includes(selected); }
scene.onselect=id=>selectPart(id);
function selectPart(id){selected=id;const part=assemblyParts(params,features).find(f=>f.id===id);if(part)domain=part.domain;renderDomains();renderTree();renderParameters();}
function renderTree(){
  if(!assemblyParts(params,features).some(f=>f.id===selected))selected="sink";
  $("feature-tree").replaceChildren();
  for(const f of assemblyParts(params,features)){
    const row=el("div",undefined,"feature-row"),b=el("button",`${f.shape==="cylinder"?"○":"▧"} ${f.name}`),toggle=el("button",hiddenParts.has(f.id)?"◌":"◉");
    b.type=toggle.type="button";b.setAttribute("aria-pressed",String(selected===f.id));b.addEventListener("click",()=>selectPart(f.id));
    toggle.setAttribute("aria-label",`${hiddenParts.has(f.id)?"Show":"Hide"} ${f.name}`);toggle.addEventListener("click",()=>{hiddenParts.has(f.id)?hiddenParts.delete(f.id):hiddenParts.add(f.id);renderTree();draw();});row.append(b,toggle);$("feature-tree").append(row);
  }
  draw();
}
$("home-view").addEventListener("click", () => scene.fit());
$("exploded").addEventListener("change", draw);
function renderDomains() {
  $("domains").replaceChildren();
  for (const [name, d] of Object.entries(domains)) {
    const b = el("button"); b.type = "button"; b.setAttribute("aria-pressed", String(name === domain)); b.style.setProperty("--domain-color", d.color);
    const text = el("span"); text.append(el("strong", name), el("small", d.subtitle)); b.append(el("span", d.icon, "domain-icon"), text);
    b.addEventListener("click", () => { domain = name; selected={Mechanical:"sink",Electrical:"board",Thermal:"sink",Photonics:"optic",Radiation:"sink"}[name];renderDomains();renderTree(); renderParameters(); draw(); if (name === "Photonics") queueEvaluation(); }); $("domains").append(b);
  }
}
function renderDraft() {
  const changed = repo.diff(current().params, params);
  $("commit").disabled = !dirty(); $("discard").disabled = !dirty();$("undo").disabled=!undoStack.length;$("redo").disabled=!redoStack.length;
  $("draft-status").textContent = dirty() ? `${changed.length} parameters, ${repo.featureDiff(featureValues(current()),features).length} parts changed · local draft` : "No uncommitted changes";
}
function renderParameters() {
  $("domain-title").textContent = domain; $("domain-description").textContent = domains[domain].description;
  $("parameters").replaceChildren();
  const custom=features.find(f=>f.id===selected);$("part-properties").hidden=!custom;
  if(custom){renderPartProperties(custom);$("domain-title").textContent=custom.name;$("domain-description").textContent="Parametric solid with revisioned placement and shop drawings. This part is not coupled to the reference thermal solve.";}
  for (const [key, f] of Object.entries(FIELDS).filter(([, f]) => !custom && f.domain === domain)) {
    const wrap = el("div", undefined, "parameter"), label = el("label", f.label), input = el("input"), out = el("output", `${params[key]} ${f.unit}`);
    input.type = "range"; input.id = `param-${key}`; input.min = f.min; input.max = f.max; input.step = f.step; input.value = params[key]; label.htmlFor = input.id; out.htmlFor = input.id;
    label.append(out); wrap.append(label, input, el("p", `${f.min}–${f.max} ${f.unit} · ${f.domain} owns this input`, "unit-note"));
    const number=el("input");number.type="number";number.min=f.min;number.max=f.max;number.step=f.step;number.value=params[key];number.setAttribute("aria-label",`${f.label} exact value`);number.className="parameter-number";
    const change=value=>{remember();params[key]=Number(value);reviewDrawings();input.value=params[key];number.value=params[key];out.value=`${params[key]} ${f.unit}`;saveDraft();renderDraft();draw();queueEvaluation();};
    input.addEventListener("input",()=>change(input.value));number.addEventListener("change",()=>{if(number.value!==""&&number.checkValidity())change(number.value);else{number.value=params[key];toast(`Enter ${f.min}–${f.max} ${f.unit}.`);}});wrap.append(number);
    $("parameters").append(wrap);
  }
  $("optical-result").hidden = domain !== "Photonics"; renderDraft(); draw();
}
function renderPartProperties(part){
  const root=$("part-properties");root.replaceChildren();
  const name=el("input");name.value=part.name;name.maxLength=60;name.setAttribute("aria-label","Part name");
  name.addEventListener("change",()=>{if(!name.value.trim()){name.value=part.name;return;}remember();part.name=name.value.trim();reviewDrawings(part.id);$("domain-title").textContent=part.name;saveDraft();renderTree();renderDraft();});root.append(name);
  for(const [key,label,min,max] of [["x","X position",-100,100],["y","Y position",-100,100],["z","Z base",-100,100],
    ["width",part.shape==="cylinder"?"Diameter":"Width",1,50],["depth","Depth",1,50],["height","Height",1,50]]){
    if(key==="depth"&&part.shape==="cylinder")continue;
    const row=el("label",`${label} / mm`),input=el("input");input.type="number";input.min=min;input.max=max;input.step=.1;input.value=part[key];input.id=`part-${key}`;row.htmlFor=input.id;
    input.addEventListener("change",()=>{if(input.value===""||!input.checkValidity()){input.value=part[key];return;}remember();part[key]=Number(input.value);reviewDrawings(part.id);if(part.shape==="cylinder"&&key==="width")part.depth=part.width;saveDraft();renderDraft();draw();});row.append(input);root.append(row);
  }
  const drawing=el("button","Prepare shop drawing");drawing.type="button";drawing.addEventListener("click",()=>{drawingPart=part.id;showView("drawings");});
  const remove=el("button","Delete part");remove.type="button";remove.addEventListener("click",()=>{remember();features=features.filter(f=>f.id!==part.id);delete drawings[part.id];selected="sink";saveDraft();renderTree();renderParameters();});root.append(drawing,remove);
}
$("add-part").addEventListener("click",()=>$("part-dialog").showModal());
$("part-form").addEventListener("submit",action(e=>{e.preventDefault();if(features.length>=16)throw new Error("The prototype supports sixteen custom solids.");remember();
  const part={id:`part-${crypto.randomUUID().slice(0,12)}`,name:$("part-name").value.trim(),shape:$("part-shape").value,domain:$("part-domain").value,x:params.width/2+20,y:0,z:0,width:12,depth:$("part-shape").value==="cylinder"?12:8,height:6};
  features.push(part);saveDraft();$("part-dialog").close();selectPart(part.id);toast("Part added to your draft. Edit dimensions, then save a revision.");
}));
$("apply-preset").addEventListener("click",()=>{remember();const preset=$("preset").value;
  params=values({});features=[];drawings={};
  if(preset==="compact")Object.assign(params,{width:32,length:35,finHeight:16,fins:10,velocity:2,current:20});
  if(preset==="optical")Object.assign(params,{width:55,length:65,fins:12,opticalPower:4,absorption:65,current:20});
  selected="sink";domain="Mechanical";saveDraft();renderDomains();renderTree();renderParameters();queueEvaluation();toast("Preset applied as an editable draft. Undo restores your previous work.");
});

function renderDrawing(){
  const parts=assemblyParts(params,features).filter(p=>p.shape!=="context");
  if(!parts.some(p=>p.id===drawingPart))drawingPart="sink";
  $("drawing-part").replaceChildren(...parts.map(p=>{const option=el("option",p.name);option.value=p.id;return option;}));$("drawing-part").value=drawingPart;
  const part=parts.find(p=>p.id===drawingPart),meta=drawings[drawingPart]||drawingDefaults(part);$("drawing-fields").replaceChildren();
  for(const [key,label] of [["number","Drawing number"],["title","Drawing title"],["revision","Drawing revision"],["material","Material / grade / condition"],["finish","Finish / surface specification"],["tolerance","General linear ± tolerance (mm)"],["projection","Projection method"],["prepared","Prepared by"],["checked","Checked by"],["status","Sheet status"],["notes","Manufacturing notes"]]){
    const labelNode=el("label",label);let input;
    if(["projection","status"].includes(key)){input=el("select");const options=key==="projection"?[["third","Third angle"],["first","First angle"]]:[["Draft","Draft"],["For review","For review"],["Released","Released by named checker"]];for(const [value,text] of options){const o=el("option",text);o.value=value;input.append(o);}}
    else if(key==="notes"){input=el("textarea");input.rows=4;input.maxLength=600;}
    else{input=el("input");input.maxLength=key==="tolerance"?8:80;if(key==="tolerance")input.inputMode="decimal";}
    input.id=`drawing-${key}`;input.value=meta[key];labelNode.htmlFor=input.id;labelNode.append(input);$("drawing-fields").append(labelNode);
    input.addEventListener("change",()=>{remember();const next={...meta,...(drawings[drawingPart]||{}),[key]:input.value};if(key!=="status"&&next.status==="Released"){next.status="For review";$("drawing-status").value=next.status;}drawings[drawingPart]=next;saveDraft();renderDraft();updateSheet();});
  }
  updateSheet();
}
function sheetState(){const part=assemblyParts(params,features).find(p=>p.id===drawingPart)||assemblyParts(params,features)[0];return{part,meta:drawings[part.id]||drawingDefaults(part)};}
function updateSheet(){
  const {part,meta}=sheetState(),svg=drawingSVG(params,part,meta,active,dirty(),"preview-arrow"),doc=new DOMParser().parseFromString(svg,"image/svg+xml");
  const print=new DOMParser().parseFromString(drawingSVG(params,part,meta,active,dirty(),"print-arrow"),"image/svg+xml");
  $("drawing-preview").replaceChildren(document.importNode(doc.documentElement,true));$("drawing-print").replaceChildren(document.importNode(print.documentElement,true));
  const problems=releaseProblems(meta,dirty());$("drawing-readiness").replaceChildren(el("strong",problems.length?"Drawing preparation checklist":"Declared release fields complete"));
  const list=el("ul");for(const problem of problems)list.append(el("li",problem));$("drawing-readiness").append(list);
  $("drawing-save").disabled=!dirty();
}
$("drawing-part").addEventListener("change",()=>{drawingPart=$("drawing-part").value;renderDrawing();});
$("drawing-form").addEventListener("submit",e=>e.preventDefault());
$("drawing-save").addEventListener("click",()=>{$("message").value=$("drawing-message").value.trim()||`Update drawing ${sheetState().meta.number}`;$("commit").click();});
$("drawing-svg").addEventListener("click",()=>{const {part,meta}=sheetState();download(drawingSVG(params,part,meta,active,dirty()),`${safeFilename(meta.number)}-${safeFilename(meta.revision)}.svg`,"image/svg+xml");});
$("drawing-dxf").addEventListener("click",()=>{const {part,meta}=sheetState();download(profileDXF(params,part),`${safeFilename(meta.number)}-profile-mm.dxf`,"application/dxf");});
$("drawing-bom").addEventListener("click",()=>download(bomCSV(params,features,drawings,dirty()?`${active}:DRAFT`:active),"emergence-parts.csv","text/csv"));
$("drawing-print-button").addEventListener("click",()=>{updateSheet();window.print();});
const safeFilename=value=>String(value).replace(/[^a-zA-Z0-9_.-]/g,"-").slice(0,80)||"drawing";
const checkNames = { "REQ-TJ": "Junction temperature ≤ 100 °C", "REQ-DROP": "PCB trace drop ≤ 50 mV", "REQ-FAN": "Heat-sink pressure drop ≤ 20 Pa" };
function renderChecks(result) {
  $("requirements").replaceChildren();
  for (const c of result.checks) {
    const card = el("article", undefined, "requirement"); card.dataset.status = c.status;
    const top = el("div", undefined, "requirement-top"); top.append(el("h3", checkNames[c.id]), el("span", c.status, "status-badge")); card.append(top);
    if (c.value) {
      const offset = c.id === "REQ-TJ" ? 273.15 : 0, scale = c.id === "REQ-DROP" ? 1000 : 1, unit = c.id === "REQ-TJ" ? "°C" : c.id === "REQ-DROP" ? "mV" : "Pa";
      const fmt = n => ((n-offset)*scale).toFixed(2);
      card.append(el("p", `${fmt(c.value.nominal)} ${unit}`, "value"), el("p", `[${fmt(c.value.lo)}, ${fmt(c.value.hi)}] ${unit} · ${c.value.rigor}`, "bounds"));
    } else card.append(el("p", c.reason || "Model unavailable for this input", "reason"));
    $("requirements").append(card);
  }
}
function queueEvaluation() {
  const gen = ++generation; clearTimeout(evaluateTimer);
  $("model-status").textContent = "Inputs changed · recalculating Rust design checks…";
  $("requirements").replaceChildren(el("p", "Recalculating checks for the current parameters…", "small muted"));
  $("handoff-values").replaceChildren(el("p","Recalculating handoffs…","small muted"));
  if (domain === "Photonics") $("optical-result").textContent = "Evaluating independent optical model…";
  evaluateTimer = setTimeout(() => evaluate(gen).catch(e => { if (gen === generation) { $("model-status").textContent = `Models unavailable: ${e.message}`; $("requirements").replaceChildren(el("p", "Checks unavailable. Reload to retry.", "small muted")); fail(e); } }), 160);
}
async function evaluate(gen) {
  const state = { ...params }, result = await request("workspace", state);
  if (gen !== generation) return;
  renderChecks(result.result); $("model-status").textContent = `Rust + WASM · ${result.ms.toFixed(0)} ms · ${result.result.computed} model runs · assumptions apply`;
  $("handoff-values").replaceChildren();
  const h=result.result.handoffs,nodes=result.result.nodes;
  for(const [name,value,detail] of [["Photonics → heat",`${h.optical_heat_w.toFixed(2)} W`,"Prescribed incident power × absorbed fraction"],["External → heat",`${h.external_heat_w.toFixed(2)} W`,"Deposited heat boundary at Q1"],["Electrical → thermal",nodes["q1.power"]?`${nodes["q1.power"].nominal.toFixed(2)} W`:"Unavailable","Temperature-dependent Q1 losses"],["Mechanical → flow",`${h.gap_mm.toFixed(2)} mm`,"Channel gap from shared geometry"]]){const card=el("article");card.append(el("strong",name),el("b",value),el("span",detail));$("handoff-values").append(card);}
  if (domain === "Photonics") {
    const optical = await request("optics", { material: 0, x: state.opticalSize, tolerance: 1e-5 });
    if (gen !== generation || domain !== "Photonics") return;
    const p = optical.result.point;
    $("optical-result").textContent = `${["Rayleigh expansion", "Folded table", "Mie series"][p.rung]} · Qext ${p.q.toFixed(5)} · reference discrepancy ${p.error.toExponential(1)} · ${p.rigor}. Mineral-dust index estimate.`;
  }
}

function renderPeers() {
  const list = peers?.list() || []; $("peer-count").textContent = list.length; $("peer-list").replaceChildren(...(list.length ? list.map(name => el("li", name)) : [el("li", "No connected sessions yet.")]));
  const heads = repo.heads(current().branch);
  $("sync-status").textContent = heads.length > 1 ? `${heads.length} divergent revisions on ${current().branch}. Compare and merge to resolve.`
    : !heads.includes(active) && !repo.ancestor(active, heads[0]) ? "Incoming revisions available. Save your draft, then compare and merge."
    : !heads.includes(active) ? "Newer revision available. Save or discard your draft to switch to it."
    : list.length ? `${list.length} connected session${list.length > 1 ? "s" : ""} · saved revisions sync directly` : "Your workspace is local. Connect a peer to exchange revisions.";
}
function renderRepository() {
  $("branch").replaceChildren();
  for (const branch of repo.branches()) for (const id of repo.heads(branch)) {
    const option = el("option", `${branch}${repo.heads(branch).length > 1 ? ` · ${id.slice(0,7)}` : ""}`); option.value = id; $("branch").append(option);
  }
  if (![...$("branch").options].some(o => o.value === active)) { const option = el("option", `${current().branch} · ${active.slice(0,7)} (working base)`); option.value = active; $("branch").append(option); }
  $("branch").value = active; $("revision-count").textContent = repo.commits.size; $("head-label").textContent = active.slice(0,7);
  renderPeers(); renderDraft(); if (view === "review") renderReview(); if (view === "history") renderHistory();
  if(view==="drawings")renderDrawing();
}
$("branch").addEventListener("change", action(() => { const id = $("branch").value; try { checkout(id); } finally { $("branch").value = active; } }));
$("new-branch").addEventListener("click", action(() => { requireClean(); $("branch-dialog").showModal(); }));
$("branch-form").addEventListener("submit", action(async e => {
  e.preventDefault(); requireClean(); const name = $("branch-name").value.trim();
  if (repo.branches().includes(name)) throw new Error("That branch already exists. Choose another name.");
  const c = await repo.commit(name, [active], params, identity(), `Create branch ${name}`,features,drawings);
  active = c.id; persist(); $("branch-dialog").close(); $("branch-name").value = ""; renderRepository(); toast(`Created ${name}.`);
}));
$("commit").addEventListener("click", action(async () => {
  if (!dirty()) return;
  const message = $("message").value.trim(); if (!message) { $("message").focus(); throw new Error("Describe the design change before saving a revision."); }
  const branch = current().branch, parent = active, values = { ...params };
  $("commit").disabled = true;
  try { const c = await repo.commit(branch, [parent], values, identity(), message,features,drawings); active = c.id;
    $("message").value = ""; persist(); renderRepository();if(view==="drawings")renderDrawing(); toast(`Saved ${c.id.slice(0,7)} to ${branch}.`);
  } finally { renderDraft(); }
}));
$("discard").addEventListener("click", () => {remember(); params=values(current().params);features=structuredClone(featureValues(current()));drawings=structuredClone(drawingValues(current())); saveDraft(); renderParameters();renderTree(); queueEvaluation(); });
$("author").addEventListener("change", () => { try { localStorage.setItem("emergence-cad-author", identity().author); } catch { /* Optional preference. */ } peers.hello(); });

function showView(name, focus = false) {
  view = name;
  for (const tab of document.querySelectorAll("[data-view]")) { const on = tab.dataset.view === name; tab.setAttribute("aria-selected", String(on)); tab.tabIndex = on ? 0 : -1; $(`${tab.dataset.view}-view`).hidden = !on; if (on && focus) tab.focus(); }
  if (name === "review") renderReview(); if (name === "history") renderHistory(); if (name === "workspace") draw();
  if(name==="drawings"){if(["sink",...features.map(f=>f.id)].includes(selected))drawingPart=selected;renderDrawing();}
}
const tabs = [...document.querySelectorAll("[data-view]")];
tabs.forEach((b, i) => {
  b.addEventListener("click", () => showView(b.dataset.view));
  b.addEventListener("keydown", e => { const n = e.key === "ArrowRight" ? (i+1)%tabs.length : e.key === "ArrowLeft" ? (i+tabs.length-1)%tabs.length : e.key === "Home" ? 0 : e.key === "End" ? tabs.length-1 : -1;
    if (n >= 0) { e.preventDefault(); showView(tabs[n].dataset.view, true); } });
});
function renderReview() {
  const previous = $("compare").value; $("compare").replaceChildren();
  for (const branch of repo.branches()) for (const id of repo.heads(branch)) if (id !== active) { const c = repo.commits.get(id), option = el("option", `${branch} · ${id.slice(0,7)} · ${c.message}`); option.value = id; $("compare").append(option); }
  if ([...$("compare").options].some(o => o.value === previous)) $("compare").value = previous;
  compare();
}
function renderNotes() {
  const id = $("compare").value || active, c = repo.commits.get(id);
  $("review-target").textContent = `${c.branch} · ${id.slice(0,7)} · ${c.message}`;
  const notes = [...repo.notes.values()].filter(n => n.commit === id).sort((a,b) => a.time.localeCompare(b.time) || a.id.localeCompare(b.id));
  $("notes").replaceChildren();
  if (!notes.length) $("notes").append(el("p", "No reviews on this revision yet.", "small muted"));
  for (const n of notes) { const card = el("article", undefined, "note"); card.append(el("div", `${n.author} · ${n.domain} · ${{ comment: "Comment", approve: "Approved", changes: "Changes requested" }[n.kind]}`, "note-top"), el("p", n.text)); $("notes").append(card); }
}
function compare() {
  ++comparisonGeneration; $("diff").replaceChildren(); $("merge").disabled = true; mergePlan = null; renderNotes();
  const source = $("compare").value;
  if (!source) { $("diff").append(el("p", "Create a branch, change a parameter and save a revision. Its differences will appear here.", "diff-summary")); return; }
  try { mergePlan = repo.merge(active, source); }
  catch (e) { $("diff").append(el("p", e.message, "diff-summary")); return; }
  const p = mergePlan, diff = repo.diff(p.target.params, p.source.params); p.choices = new Set();
  $("diff").append(el("p", `Comparing ${p.source.branch} with ${p.target.branch}. Shared base ${p.base.id.slice(0,7)}.`, "diff-summary"));
  if (!diff.length) $("diff").append(el("p", "Parameter values are identical. Histories can still be joined.", "small muted"));
  for (const key of diff) {
    const f = FIELDS[key], row = el("div", undefined, "diff-row"); row.append(el("strong", `${f.label} · ${f.domain}`), el("p", `${values(p.target.params)[key]} → ${values(p.source.params)[key]} ${f.unit}`, "delta"));
    if (p.conflicts.includes(key)) {
      row.append(el("p", `Conflict · both branches changed this input from ${values(p.base.params)[key]}. Choose a value.`, "conflict"));
      const select = el("select"); select.setAttribute("aria-label", `Resolve ${f.label}`);
      for (const [value, label] of [["", "Select a resolution…"], ["target", `Keep active: ${values(p.target.params)[key]} ${f.unit}`], ["source", `Use incoming: ${values(p.source.params)[key]} ${f.unit}`]]) { const option = el("option", label); option.value = value; select.append(option); }
      select.addEventListener("change", () => { if (!select.value) p.choices.delete(key); else { p.choices.add(key); p.params[key] = values((select.value === "source" ? p.source : p.target).params)[key]; } evaluateMerge(); }); row.append(select);
    } else if (p.params[key] === values(p.target.params)[key]) row.append(el("p", "Active value retained: the incoming branch did not change this parameter from the shared base.", "small muted"));
    $("diff").append(row);
  }
  for(const id of repo.featureDiff(featureValues(p.target),featureValues(p.source))){
    const a=featureValues(p.target).find(f=>f.id===id),b=featureValues(p.source).find(f=>f.id===id),row=el("div",undefined,"diff-row");
    row.append(el("strong",`Part: ${b?.name||a?.name}`),el("p",!a?"Added in incoming branch":!b?"Removed in incoming branch":"Dimensions, placement or description differ","small muted"));
    if(p.featureConflicts.includes(id)){
      const key=`part:${id}`;p.conflicts.push(key);row.append(el("p","Conflict · this part changed on both branches.","conflict"));
      const select=el("select");select.setAttribute("aria-label",`Resolve part ${a?.name||b?.name}`);
      for(const [value,label] of [["","Select a resolution…"],["target",a?`Keep active ${a.width} × ${a.depth} × ${a.height} at (${a.x}, ${a.y}, ${a.z})`:"Keep deletion"],["source",b?`Use incoming ${b.width} × ${b.depth} × ${b.height} at (${b.x}, ${b.y}, ${b.z})`:"Use incoming deletion"]]){const o=el("option",label);o.value=value;select.append(o);}
      select.addEventListener("change",()=>{if(!select.value)p.choices.delete(key);else{p.choices.add(key);p.features=p.features.filter(f=>f.id!==id);const f=select.value==="source"?b:a;if(f)p.features.push(structuredClone(f));}evaluateMerge();});row.append(select);
    }$("diff").append(row);
  }
  const td=drawingValues(p.target),sd=drawingValues(p.source);
  for(const id of new Set([...Object.keys(td),...Object.keys(sd)]))if(canonical(td[id])!==canonical(sd[id])){
    const row=el("div",undefined,"diff-row");row.append(el("strong",`Shop drawing: ${sd[id]?.number||td[id]?.number||id}`),el("p","Manufacturing metadata differs: material, finish, tolerances or release fields.","small muted"));
    if(p.drawingConflicts.includes(id)){
      const key=`sheet:${id}`;p.conflicts.push(key);row.append(el("p","Conflict · drawing metadata changed on both branches.","conflict"));
      const select=el("select");select.setAttribute("aria-label",`Resolve drawing ${id}`);
      for(const [value,label] of [["","Select a resolution…"],["target","Keep active drawing fields"],["source","Use incoming drawing fields"]]){const o=el("option",label);o.value=value;select.append(o);}
      select.addEventListener("change",()=>{if(!select.value)p.choices.delete(key);else{p.choices.add(key);const d=(select.value==="source"?sd:td)[id];if(d)p.drawings[id]=structuredClone(d);else delete p.drawings[id];}evaluateMerge();});row.append(select);
    }$("diff").append(row);
  }
  p.checks = el("p", "", "review-checks"); $("diff").append(p.checks); evaluateMerge();
}
async function evaluateMerge() {
  const p = mergePlan, gen = ++comparisonGeneration; $("merge").disabled = true;
  if (!p) return;
  if (p.conflicts.some(k => !p.choices.has(k))) { p.checks.textContent = "Resolve every conflict before evaluating the merge."; return; }
  if(!validFeatures(p.features)||!validDrawings(p.drawings)){p.checks.textContent="The combined parts or drawings exceed the supported limits. Reconcile them on a branch first.";return;}
  if (repo.ancestor(p.source.id, active)) { p.checks.textContent = "The incoming revision is already included in your active history."; return; }
  if (dirty()) { p.checks.textContent = "Save or discard your draft before merging."; return; }
  p.checks.textContent = "Recalculating the proposed merge in Rust…";
  try { const r = await request("workspace", { ...p.params }); if (gen !== comparisonGeneration || mergePlan !== p) return;
    p.checks.textContent = `Proposed merge · ${r.result.checks.map(c => `${c.id}: ${c.status}`).join(" · ")}. Model assumptions and uncertainty apply.`;
    $("merge").disabled = false;
  } catch (e) { if (gen === comparisonGeneration) p.checks.textContent = `Merge checks unavailable: ${e.message}`; }
}
$("compare").addEventListener("change", compare);
$("merge").addEventListener("click", action(async () => {
  requireClean(); const p = mergePlan; if (!p || $("merge").disabled || p.target.id !== active) throw new Error("Reopen the comparison before merging.");
  $("merge").disabled = true;
  const c = await repo.commit(current().branch, [active, p.source.id], p.params, identity(), `Merge ${p.source.branch} (${p.source.id.slice(0,7)})`,p.features,p.drawings);
  active = c.id; params = values(c.params);features=structuredClone(featureValues(c));drawings=structuredClone(drawingValues(c)); persist(); renderRepository(); renderParameters();renderTree(); queueEvaluation(); toast(`Merged into ${c.branch}. All source revisions are retained.`);
}));
$("post-review").addEventListener("click", action(async () => {
  await repo.note($("compare").value || active, identity(), $("review-text").value, $("review-kind").value);
  $("review-text").value = ""; persist(); renderNotes(); toast("Review saved and shared with connected peers.");
}));
function renderHistory() {
  $("history").replaceChildren();
  const cs = [...repo.commits.values()].sort((a,b) => b.time.localeCompare(a.time) || a.id.localeCompare(b.id));
  for (const c of cs) {
    const card = el("article", undefined, "history-item"); card.append(el("h3", c.message), el("p", `${c.author} · ${c.domain} · ${c.branch} · ${new Date(c.time).toLocaleString()}`), el("code", `${c.id.slice(0,12)}${c.id === active ? " · active revision" : ""}`));
    const parent = repo.commits.get(c.parents[0]);
    if (parent) { const list = el("ul"); for (const key of repo.diff(parent.params, c.params)) list.append(el("li", `${FIELDS[key].label}: ${values(parent.params)[key]} → ${values(c.params)[key]} ${FIELDS[key].unit}`));
      for(const id of repo.featureDiff(featureValues(parent),featureValues(c)))list.append(el("li",`Part ${featureValues(c).find(f=>f.id===id)?.name||featureValues(parent).find(f=>f.id===id)?.name}: geometry or placement changed`));
      if(canonical(drawingValues(parent))!==canonical(drawingValues(c)))list.append(el("li","Manufacturing drawing fields changed"));if (list.childNodes.length) card.append(list); }
    if (c.parents.length) card.append(el("p", `Parents: ${c.parents.map(id=>id.slice(0,7)).join(" + ")}`));
    if (c.id !== active) { const button = el("button", "Open revision"); button.type = "button"; button.addEventListener("click", action(() => { checkout(c.id); showView("workspace"); })); card.append(button); } $("history").append(card);
  }
}

// Bundle ingestion is serialized, including across multiple peer transports.
// An incoming fast-forward updates a clean workspace; a draft or divergence
// stays on its original base until explicitly committed or merged.
let incoming = Promise.resolve();
function receive(bundle) {
  incoming = incoming.then(async () => {
    const changed = await repo.ingest(bundle); if (!changed) return true;
    const tips = repo.heads(current().branch), wasDirty = dirty();
    if (!wasDirty && tips.length === 1 && repo.ancestor(active, tips[0])) { active = tips[0]; params=values(current().params);features=structuredClone(featureValues(current()));drawings=structuredClone(drawingValues(current())); renderParameters();renderTree(); queueEvaluation(); }
    persist(); renderRepository(); return true;
  }).catch(e => { fail(e); return false; });
  return incoming;
}
peers = new Peers(identity, receive, renderPeers, fail); peers.snapshot = () => repo.bundle();
// Ask existing tabs for a snapshot only after the callback has been assigned.
peers.hello();
$("connect").addEventListener("click", () => { renderPeers(); $("peer-dialog").showModal(); });
for (const b of document.querySelectorAll("[data-close]")) b.addEventListener("click", () => $(b.dataset.close).close());
function iceServers() {
  const servers = JSON.parse($("ice").value);
  if (!Array.isArray(servers) || servers.length > 8 || !servers.every(s => s && typeof s === "object" && (typeof s.urls === "string" || Array.isArray(s.urls)))) throw new Error("Enter an ICE server array, for example [{\"urls\":\"stun:your-server:3478\"}].");
  return servers;
}
async function connectStep(fn, success) {
  const buttons = [$("make-offer"), $("make-answer"), $("finish-peer")]; buttons.forEach(b => b.disabled = true); $("connection-status").textContent = "Discovering network candidates…";
  try { const code = await fn(); if (code) $("peer-output").value = code; $("connection-status").textContent = success; }
  catch (e) { $("connection-status").textContent = e.message; fail(e); }
  finally { buttons.forEach(b => b.disabled = false); }
}
$("make-offer").addEventListener("click", () => connectStep(() => peers.offer(iceServers()), "Offer ready. Send your code to your collaborator, then paste their answer above."));
$("make-answer").addEventListener("click", () => connectStep(() => peers.answer($("peer-input").value, iceServers()), "Answer ready. Send your code back to the person who created the offer."));
$("finish-peer").addEventListener("click", () => connectStep(() => peers.finish($("peer-input").value), "Answer accepted. Waiting for the direct connection; connected sessions appear below."));
$("copy-code").addEventListener("click", action(async () => { if (!$("peer-output").value) throw new Error("Create an offer or answer first."); try { await navigator.clipboard.writeText($("peer-output").value); toast("Connection code copied."); } catch { $("peer-output").focus(); $("peer-output").select(); toast("Select and copy the code from the text box."); } }));
$("disconnect").addEventListener("click", () => { peers.disconnect(); $("connection-status").textContent = "Direct peers disconnected. Browser-tab sync remains active."; });
function download(data, filename, type) {
  const url = URL.createObjectURL(new Blob([data], { type })), link = el("a"); link.href = url; link.download = filename; link.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
}
$("export").addEventListener("click", () => { download(JSON.stringify(repo.bundle()), "emergence-motor-driver.json", "application/json"); toast(dirty() ? "Project exported. Your draft is not committed; save it to include it in an export." : "Full project history and reviews exported."); });
$("import").addEventListener("change", action(async () => {
  const file = $("import").files[0]; $("import").value = ""; if (!file) return;
  if (file.size > 3_000_000) throw new Error("Project file exceeds this demo’s 3 MB import limit.");
  const bundle = JSON.parse(await file.text()); if (await receive(bundle)) toast("Project imported. Existing revisions and drafts are preserved.");
}));
$("stl").addEventListener("click", action(async () => {
  const b = $("stl"), target=selected,part=assemblyParts(params,features).find(f=>f.id===target); b.disabled = true; b.textContent = "Meshing in Rust…";
  try {
    const { result: mesh } = await request("mesh", { params:{...params},features:structuredClone(features),target });
    const buffer = new ArrayBuffer(84 + mesh.triangles.length*50), data = new DataView(buffer);
    const header = new TextEncoder().encode(`Emergence CAD; millimetres; ${target}; cell=${mesh.cell_mm}mm`); new Uint8Array(buffer).set(header.slice(0,80)); data.setUint32(80,mesh.triangles.length,true);
    let offset = 84;
    for (const t of mesh.triangles) {
      const [a,b,c] = t.map(i => mesh.vertices[i]), u=b.map((n,i)=>n-a[i]),v=c.map((n,i)=>n-a[i]);
      const normal=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]], size=Math.hypot(...normal)||1;
      for (const n of [...normal.map(n=>n/size),...a,...b,...c]) { data.setFloat32(offset,n,true); offset+=4; } data.setUint16(offset,0,true); offset+=2;
    }
    download(buffer, `${safeFilename(part.name)}-mm.stl`, "model/stl"); toast(`${mesh.triangles.length.toLocaleString()} triangles · ${mesh.watertight ? "watertight" : "not watertight"} · millimetres · ${mesh.cell_mm.toFixed(2)} mm mesh rounds sharp edges.`);
  } finally { b.disabled = false; b.textContent = "↓ Export selected STL"; }
}));
addEventListener("beforeunload", e => { if (dirty() && !storageOK) { e.preventDefault(); e.returnValue = ""; } });
addEventListener("pagehide", e => { if (!e.persisted) peers.close(); });
renderDomains();renderTree(); renderParameters(); renderRepository(); saveDraft(); queueEvaluation();
