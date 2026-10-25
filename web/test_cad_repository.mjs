import test from "node:test";
import assert from "node:assert/strict";
import { Repository, DEFAULTS, hash, values, validDrawings } from "./cad/repository.js";
import { drawingDefaults, profileDXF } from "./cad/drawings.js";

const me = { author: "Mechanical engineer", domain: "Mechanical" };
const thermal = { author: "Thermal engineer", domain: "Thermal" };
test("replicas converge without losing concurrent revisions; independent parameters merge", async () => {
  const a = await new Repository().init(), b = await new Repository().init();
  assert.equal(a.root, b.root);
  const ca = await a.commit("main", [a.root], { ...DEFAULTS, fins: 12 }, me, "More fins");
  const cb = await b.commit("main", [b.root], { ...DEFAULTS, ambient: 50 }, thermal, "Hotter inlet");
  await a.ingest(b.bundle()); await b.ingest(a.bundle());
  assert.deepEqual(a.heads("main"), b.heads("main")); assert.equal(a.heads("main").length, 2);
  const merge = a.merge(ca.id, cb.id); assert.deepEqual(merge.conflicts, []);
  assert.equal(merge.params.fins, 12); assert.equal(merge.params.ambient, 50);
  const c = await a.commit("main", [ca.id, cb.id], merge.params, me, "Merge thermal input");
  await b.ingest(a.bundle()); assert.deepEqual(b.heads("main"), [c.id]);
  assert.equal(b.commits.size, 4); assert.equal(await b.ingest(a.bundle()), false);
});
test("same-input edits remain explicit conflicts; same-value edits do not conflict", async () => {
  const r = await new Repository().init();
  const a = await r.commit("main", [r.root], { ...DEFAULTS, fins: 11 }, me, "Eleven fins");
  const b = await r.commit("thermal/alternative", [r.root], { ...DEFAULTS, fins: 12 }, thermal, "Twelve fins");
  assert.deepEqual(r.merge(a.id, b.id).conflicts, ["fins"]);
  const same = await r.commit("another", [r.root], { ...DEFAULTS, fins: 11 }, me, "Same choice");
  assert.deepEqual(r.merge(a.id, same.id).conflicts, []);
});
test("three-way merge retains a newer target when the source did not change that input", async () => {
  const r = await new Repository().init();
  const a = await r.commit("main", [r.root], { ...DEFAULTS, current: 20 }, me, "Derate current");
  const b = await r.commit("cooling", [r.root], { ...DEFAULTS, fins: 12 }, thermal, "Change cooling");
  const p = r.merge(a.id, b.id); assert.equal(p.params.current, 20); assert.equal(p.params.fins, 12);
  assert.equal(p.base.id, r.root);
});
test("tampered revisions, incomplete history and unrelated roots fail atomically", async () => {
  const a = await new Repository().init(), b = await new Repository().init();
  await b.commit("cooling", [b.root], { ...DEFAULTS, fins: 12 }, thermal, "Change cooling");
  const bundle = structuredClone(b.bundle()); bundle.commits[1].params.fins = 13;
  await assert.rejects(a.ingest(bundle), /verification/); assert.equal(a.commits.size, 1);
  const body = { ...bundle.commits[1], parents: ["a".repeat(64)] }; delete body.id;
  bundle.commits[1] = { ...body, id: await hash(body) };
  await assert.rejects(a.ingest(bundle), /Incomplete/); assert.equal(a.commits.size, 1);
  const root = { ...bundle.commits[0], message: "Different root" }; delete root.id;
  await assert.rejects(a.ingest({ ...bundle, commits: [{ ...root, id: await hash(root) }] }), /Unrelated/);
});
test("reviews round-trip and tampered notes do not partly import new commits", async () => {
  const a = await new Repository().init(), b = await new Repository().init();
  const c = await b.commit("cooling", [b.root], { ...DEFAULTS, fins: 12 }, thermal, "Change cooling");
  const note = await b.note(c.id, me, "Check channel pressure before merging.", "changes");
  const corrupt = structuredClone(b.bundle()); corrupt.notes[0].text = "Approved";
  await assert.rejects(a.ingest(corrupt), /review failed/); assert.equal(a.commits.size, 1);
  await a.ingest(b.bundle()); assert.equal(a.notes.get(note.id).kind, "changes");
  const restored = await new Repository().init(); await restored.ingest(JSON.parse(JSON.stringify(a.bundle())));
  assert.deepEqual(restored.heads("cooling"), [c.id]); assert.equal(restored.notes.get(note.id).text, note.text);
});
test("criss-cross merges with multiple best bases stop for explicit reconciliation", async () => {
  const r = await new Repository().init();
  const a = await r.commit("a", [r.root], { ...DEFAULTS, fins: 11 }, me, "A");
  const b = await r.commit("b", [r.root], { ...DEFAULTS, ambient: 50 }, thermal, "B");
  const merged = r.merge(a.id,b.id).params;
  const ab = await r.commit("a", [a.id,b.id], merged, me, "AB");
  const ba = await r.commit("b", [a.id,b.id], merged, thermal, "BA");
  assert.throws(() => r.merge(ab.id,ba.id), /multiple merge bases/);
});
const block={id:"mount",name:"Mount",shape:"box",domain:"Mechanical",x:30,y:0,z:0,width:12,depth:8,height:6};
test("legacy revisions keep their root and acquire expanded defaults without hash rewriting",async()=>{
  const r=await new Repository().init();
  assert.equal(r.root.slice(0,7),"ae4b7ab");assert.equal(Object.keys(r.commits.get(r.root).params).length,6);
  assert.equal(values(r.commits.get(r.root).params).width,40);
});
test("parts merge independently; edit versus delete remains a conflict",async()=>{
  const r=await new Repository().init();
  const base=await r.commit("main",[r.root],DEFAULTS,me,"Add mount",[block]);
  const a=await r.commit("main",[base.id],DEFAULTS,me,"Widen mount",[{...block,width:14}]);
  const b=await r.commit("alternate",[base.id],DEFAULTS,thermal,"Remove mount",[]);
  assert.deepEqual(r.merge(a.id,b.id).featureConflicts,["mount"]);
  const c=await r.commit("second",[base.id],DEFAULTS,me,"Add spacer",[block,{...block,id:"spacer",name:"Spacer",shape:"cylinder",width:8,depth:8}]);
  const p=r.merge(a.id,c.id);assert.deepEqual(p.featureConflicts,[]);assert.equal(p.features.find(f=>f.id==="mount").width,14);assert.equal(p.features.length,2);
  const restored=await new Repository().init();await restored.ingest(r.bundle());assert.equal(restored.commits.get(c.id).features[1].shape,"cylinder");
});
test("manufacturing sheets are hashed and concurrent sheet edits require review",async()=>{
  const r=await new Repository().init(),drawing=drawingDefaults({id:"sink",name:"Sink"});
  assert.equal(validDrawings({sink:{...drawing,status:"Released"}}),false);
  const released={...drawing,material:"Aluminium 6063-T5",finish:"Anodized",tolerance:"0.1",prepared:"Designer",checked:"Checker",status:"Released"};
  assert.equal(validDrawings({sink:released}),true);
  const base=await r.commit("main",[r.root],DEFAULTS,me,"Document sink",[],{sink:drawing});
  const a=await r.commit("main",[base.id],DEFAULTS,me,"Finish",[],{sink:{...drawing,finish:"Anodized"}});
  const b=await r.commit("shop",[base.id],DEFAULTS,me,"Tolerance",[],{sink:{...drawing,tolerance:"0.1"}});
  assert.deepEqual(r.merge(a.id,b.id).drawingConflicts,["sink"]);
  const bundle=structuredClone(r.bundle());bundle.commits.at(-1).drawings.sink.tolerance="0.2";
  await assert.rejects(r.ingest(bundle),/verification/);
});
test("DXF end-section encloses the exact analytic extrusion area",()=>{
  const dxf=profileDXF(DEFAULTS,{id:"sink"}),lines=dxf.split("\n"),points=[];
  for(let i=0;i<lines.length;i++)if(lines[i]==="LINE"){
    const item={};for(let j=i+1;j<lines.length&&lines[j]!=="0";j+=2)item[lines[j]]=Number(lines[j+1]);points.push([[item["10"],item["20"]],[item["11"],item["21"]]]);
  }
  assert.equal(points.length,4*DEFAULTS.fins);
  for(let i=0;i<points.length;i++)assert.deepEqual(points[i][1],points[(i+1)%points.length][0]);
  const area=points.reduce((sum,[a,b])=>sum+a[0]*b[1]-b[0]*a[1],0)/2;
  const exact=DEFAULTS.width*DEFAULTS.base+DEFAULTS.fins*DEFAULTS.finThickness*DEFAULTS.finHeight;
  assert.ok(Math.abs(area-exact)<.01,`${area} versus ${exact}`);
});
