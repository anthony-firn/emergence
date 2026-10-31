import test from "node:test";
import assert from "node:assert/strict";
import {pcbExample, siliconExample, documentFor, validateDocument, pads, bodyBounds, routePoints, connectPads, moveComponent, segmentDistance, pcbChecks, siliconChecks, pcbBOM} from "./cad/studio-model.js";

test("segment distances agree with exact crossing, parallel, endpoint and degenerate cases",()=>{
  assert.equal(segmentDistance([0,0],[4,4],[0,4],[4,0]),0);
  assert.equal(segmentDistance([0,0],[4,0],[0,3],[4,3]),3);
  assert.equal(segmentDistance([0,0],[1,0],[4,4],[4,6]),5);
  assert.equal(segmentDistance([0,0],[0,0],[3,4],[3,4]),5);
  assert.equal(segmentDistance([0,0],[4,0],[2,0],[6,0]),0);
});
test("rotated footprints and 45-degree routes retain their analytic dimensions",()=>{
  const board=pcbExample(),j=board.components[0],p=pads(board).filter(p=>p.component===j.id);
  assert.deepEqual(p.map(p=>[p.x,p.y]),[[10,14],[10,18]]);
  assert.deepEqual(bodyBounds(j),{x:6.5,y:12,width:7,height:8});
  assert.deepEqual(routePoints({x:2,y:3},{x:12,y:7}),[[2,3],[6,7],[12,7]]);
  assert.deepEqual(routePoints({x:12,y:7},{x:2,y:3}),[[12,7],[8,3],[2,3]]);
  assert.deepEqual(routePoints({x:2,y:3},{x:2,y:7}),[[2,3],[2,7]]);
});
test("pad connections reduce the explicit connectivity count and are removed on a footprint move",()=>{
  const b=pcbExample(),before=pcbChecks(b).unrouted,r=connectPads(b,"J1:1","U1:1","top",.5);
  b.routes.push(r);assert.equal(pcbChecks(b).unrouted,before-1);
  const saved=structuredClone(b);
  saved.routes[0].points[0][0]+=.1;
  assert.equal(pcbChecks(saved).unrouted,before);
  assert(pcbChecks(saved).issues.some(i=>i.code==="endpoint"));
  assert.equal(moveComponent(b,b.components[0],{x:11}),1);
  assert.equal(b.routes.length,0);assert.equal(pcbChecks(b).unrouted,before);
  assert.throws(()=>connectPads(b,"J1:1","J1:2","top",.5),/same assigned net/);
});
test("copper checks find crossing nets on the same layer; different layers remain separate",()=>{
  const b={width:30,height:30,grid:.5,rules:{clearance:.25,minWidth:.3,edge:1},components:[],routes:[]};
  for(const [id,x,y,net] of [["A1",7,10,"A"],["A2",23,20,"A"],["B1",7,20,"B"],["B2",23,10,"B"]])b.components.push({id,name:id,footprint:"connector",x,y,rotation:0,nets:[net,net]});
  b.routes.push(connectPads(b,"A1:1","A2:1","top",.5),connectPads(b,"B1:1","B2:1","top",.5));
  assert(pcbChecks(b).issues.some(i=>i.code==="track-clearance"));
  b.routes[1].layer="bottom";assert(!pcbChecks(b).issues.some(i=>i.code==="track-clearance"));
  b.routes[0].width=.1;assert(pcbChecks(b).issues.some(i=>i.code==="width"));
  b.components[0].x=.5;assert(pcbChecks(b).issues.some(i=>i.code==="pad-edge"));
});
test("pad clearance is checked at the exact circular copper boundary",()=>{
  const b={width:30,height:30,grid:.1,rules:{clearance:.4,minWidth:.3,edge:1},components:[
    {id:"C1",name:"Capacitor",footprint:"capacitor",x:10,y:10,rotation:0,nets:["A","B"]},
    {id:"C2",name:"Capacitor",footprint:"capacitor",x:10,y:12,rotation:0,nets:["C","D"]},
  ],routes:[]};
  const pads=()=>pcbChecks(b).issues.filter(i=>i.code==="pad-clearance");
  assert.equal(pads().length,0); // Centre spacing 2 - 2*.8 = .4 mm.
  b.components[1].y=11.999;assert(pads().length>0);
});
test("floorplan checks agree with exact rectangular area, overlap, bounds and halos",()=>{
  const d={width:100,height:100,grid:1,margin:10,spacing:5,target:70,blocks:[
    {id:"a",name:"A",kind:"logic",x:10,y:10,width:20,height:20,fixed:false},
    {id:"b",name:"B",kind:"memory",x:35,y:10,width:20,height:20,fixed:false},
  ]};
  const result=siliconChecks(d);assert.deepEqual(result.issues,[]);
  assert.equal(result.area,800);assert.equal(result.coreArea,6400);assert.equal(result.utilization,12.5);
  d.blocks[1].x=34;assert(siliconChecks(d).issues.some(i=>i.code==="spacing"));
  d.blocks[1].x=29;assert(siliconChecks(d).issues.some(i=>i.code==="overlap"));
  d.blocks[0].x=9;assert(siliconChecks(d).issues.some(i=>i.code==="boundary"));
  d.target=10;assert(siliconChecks(d).issues.some(i=>i.code==="density"));
});
test("documents round-trip, reject broken references and nonfinite values, and discard extra keys",()=>{
  for(const [kind,design] of [["pcb",pcbExample()],["silicon",siliconExample()]]) {
    const doc=documentFor(kind,design);
    assert.deepEqual(validateDocument(JSON.parse(JSON.stringify(doc))),doc);
    doc.design.width=null;assert.throws(()=>validateDocument(doc),/dimensions/);
  }
  const doc=documentFor("pcb",pcbExample());
  doc.design.routes.push({...connectPads(doc.design,"J1:1","U1:1","top",.5),to:"missing:1"});
  assert.throws(()=>validateDocument(doc),/connections/);
  doc.design.routes=[];
  doc.design.components[0].footprint="__proto__";assert.throws(()=>validateDocument(doc),/footprints/);
  const fresh=documentFor("silicon",siliconExample());fresh.design.blocks[0].injected="ignored";
  assert.equal(validateDocument(fresh).design.blocks[0].injected,undefined);
  fresh.design.margin=500;assert.throws(()=>validateDocument(fresh),/constraints/);
});
test("BOM preserves quoted user names as data",()=>{
  const b=pcbExample();b.components[0].name='A "quoted", part';
  assert(pcbBOM(b).includes('"A ""quoted"", part"'));
  assert.equal(pcbChecks(pcbExample()).issues.length,0);
  assert.equal(siliconChecks(siliconExample()).issues.length,0);
});
