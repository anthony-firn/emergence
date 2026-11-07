import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { assemblyGeometry, assemblyParts } from './cad/scene.js';
const docs=JSON.parse(await readFile(new URL('./cad/projects/inspection-drone.json',import.meta.url))).documents;
const params=docs.mech.params;
const context={board:docs.pcb,pcbPlacement:docs.project.scene.placements.pcb,mechPlacement:docs.project.scene.placements.mech};
const extent=(item,axis,fn)=>Math[fn](...item.points.map(p=>p[axis]));
const centre=(item,axis)=>(extent(item,axis,'min')+extent(item,axis,'max'))/2;
const close=(a,b)=>assert(Math.abs(a-b)<1e-9,`${a} != ${b}`);

test('the drone PCB and all six MOSFETs use their PCB positions in the heat sink frame',()=>{
  const items=assemblyGeometry(params,[],context), board=items.find(p=>p.id==='board');
  // PCB has no turn; sink is turned 90 degrees. Shared phase and translation cancel.
  close(extent(board,0,'max')-extent(board,0,'min'),docs.pcb.height);
  close(extent(board,1,'max')-extent(board,1,'min'),docs.pcb.width);
  close(extent(board,2,'min'),-3);
  for(const c of docs.pcb.components.filter(c=>/^Q[1-6]$/.test(c.id))){
    const item=items.find(p=>p.id===(c.id==='Q1'?'q1':`pcb-${c.id}`));
    close(centre(item,0),docs.pcb.height/2-c.y);
    close(centre(item,1),-(c.x-docs.pcb.width/2));
    close(extent(item,2,'min'),extent(board,2,'max'));
    assert(extent(item,2,'max')<0,'device must sit below the sink base, not on its fins');
  }
  assert(!items.some(p=>p.id==='optic'),'a separate camera is not a component on this ESC');
  assert(assemblyParts(params,[],context).some(p=>p.name.startsWith('Q6 ·')));
});

test('PCB edits move the mechanical package; fin height does not move a soldered device',()=>{
  const a=assemblyGeometry(params,[],context).find(p=>p.id==='q1');
  const changed=structuredClone(context);changed.board.components.find(p=>p.id==='Q1').x+=10;
  const b=assemblyGeometry({...params,finHeight:params.finHeight+15},[],changed).find(p=>p.id==='q1');
  close(centre(b,1)-centre(a,1),-10);close(centre(b,2),centre(a,2));
});

test('exploding the assembly keeps devices attached to their PCB',()=>{
  const a=assemblyGeometry(params,[],context),b=assemblyGeometry(params,[],context,true);
  for(let i=0;i<a.length;i++)for(let j=0;j<a[i].points.length;j++){
    close(b[i].points[j][0],a[i].points[j][0]);close(b[i].points[j][1],a[i].points[j][1]);
    close(b[i].points[j][2]-a[i].points[j][2],a[i].id==='sink'?0:-15);
  }
});

test('standalone schematic Q1 sits on its board even with tall fins',()=>{
  const items=assemblyGeometry({...params,finHeight:40}),board=items.find(p=>p.id==='board'),q=items.find(p=>p.id==='q1');
  close(extent(q,2,'min'),extent(board,2,'max'));assert(extent(q,2,'max')<0);
});
