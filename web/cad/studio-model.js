// Browser design documents. Dimensions are user inputs, not process constants.
// Footprints and floorplan blocks below are illustrative estimates, not vendor
// land patterns or PDK geometry. Geometric checks make no fabrication claim.
export const STUDIO_FORMAT = "emergence-studio/1";
export const clone = value => structuredClone(value);
const finite = (n, lo, hi) => Number.isFinite(n) && n >= lo && n <= hi;
const text = (s, max = 60) => typeof s === "string" && s.trim().length > 0 && s.length <= max;
const id = s => typeof s === "string" && /^[a-zA-Z][a-zA-Z0-9_-]{0,47}$/.test(s);
const point = p => Array.isArray(p) && p.length === 2 && p.every(n => finite(n, -1000, 1000));
const unique = xs => new Set(xs).size === xs.length;
export const snap = (n, grid) => Math.round(n / grid) * grid;

export const FOOTPRINTS = {
  connector: { title: "2-pin connector", width: 8, height: 7, pads: [[-2, 0], [2, 0]] },
  resistor: { title: "Axial resistor", width: 9, height: 3, pads: [[-5, 0], [5, 0]] },
  capacitor: { title: "Radial capacitor", width: 6, height: 6, pads: [[-1.5, 0], [1.5, 0]] },
  transistor: { title: "3-pin power device", width: 8, height: 6, pads: [[-3, 3], [0, 3], [3, 3]] },
  dip8: { title: "DIP-8", width: 8, height: 12, pads: [[-4,-4.5],[-4,-1.5],[-4,1.5],[-4,4.5],[4,4.5],[4,1.5],[4,-1.5],[4,-4.5]] },
};
export function pcbExample() {
  return {
    width: 80, height: 55, grid: .5,
    rules: { clearance: .25, minWidth: .3, edge: 1 },
    components: [
      { id:"J1", name:"Power in", footprint:"connector", x:10, y:16, rotation:90, nets:["VCC","GND"] },
      { id:"U1", name:"Gate driver", footprint:"dip8", x:32, y:26, rotation:0, nets:["VCC","PWM","GND","GATE","GND","GATE","PWM","VCC"] },
      { id:"Q1", name:"Power switch", footprint:"transistor", x:59, y:26, rotation:0, nets:["GATE","OUT","GND"] },
      { id:"R1", name:"Gate resistor", footprint:"resistor", x:46, y:16, rotation:0, nets:["GATE","GATE"] },
      { id:"C1", name:"Decoupling", footprint:"capacitor", x:19, y:32, rotation:0, nets:["VCC","GND"] },
      { id:"J2", name:"Load", footprint:"connector", x:69, y:42, rotation:90, nets:["OUT","GND"] },
    ],
    routes: [],
    stackup: { thickness: 1.6, copper: 1, epsr: 4.4 },
    currents: { VCC: .5, PWM: .05, GATE: .4, OUT: 30, GND: 30 },
    powers: { U1: .35, R1: .05 },
    loadNets: ["OUT", "GND"],
  };
}
// Typical block dissipation by kind, mW: illustrative starting values.
export const BLOCK_POWER = { logic: 300, memory: 150, io: 90, analog: 40 };
export function siliconExample() {
  return {
    width:1000, height:800, grid:10, margin:40, spacing:20, target:70,
    blocks:[
      { id:"cpu", name:"CPU core", kind:"logic", x:100, y:100, width:300, height:250, fixed:false, power:420 },
      { id:"sram", name:"SRAM macro", kind:"memory", x:500, y:100, width:350, height:280, fixed:true, power:160 },
      { id:"io", name:"I/O controller", kind:"io", x:100, y:480, width:250, height:180, fixed:false, power:90 },
      { id:"dsp", name:"Signal processing", kind:"logic", x:460, y:470, width:300, height:200, fixed:false, power:310 },
      { id:"pwm", name:"Motor PWM (Tiny Tapeout)", kind:"logic", x:790, y:470, width:150, height:200, fixed:false, power:5 },
    ],
    // RTL implementing blocks: emergence's mini8 CPU and a public Tiny Tapeout
    // motor controller (Apache-2.0), adapted to the native Verilog subset.
    chip: { rtl:[{ block:"cpu", source:"mini8", verilog:null }, { block:"pwm", source:"tt-motor", verilog:null }], frequency:100, vdd:1.2, activity:20, instances:1, theta:45, ambient:25 },
  };
}
export function pads(board) {
  return board.components.flatMap(c => {
    const angle = c.rotation * Math.PI / 180, cs = Math.round(Math.cos(angle)), sn = Math.round(Math.sin(angle));
    return FOOTPRINTS[c.footprint].pads.map(([x,y], i) => ({
      id:c.id + ":" + (i+1), component:c.id, number:i+1, net:c.nets[i],
      x:c.x+x*cs-y*sn, y:c.y+x*sn+y*cs, radius:.8,
    }));
  });
}
export function bodyBounds(c) {
  const f = FOOTPRINTS[c.footprint], rotated = c.rotation % 180 !== 0;
  const width = rotated ? f.height : f.width, height = rotated ? f.width : f.height;
  return { x:c.x-width/2, y:c.y-height/2, width, height };
}
export function routePoints(a, b) {
  const dx=b.x-a.x, dy=b.y-a.y, d=Math.min(Math.abs(dx),Math.abs(dy));
  const middle=[a.x+Math.sign(dx)*d, a.y+Math.sign(dy)*d];
  return [[a.x,a.y],middle,[b.x,b.y]].filter((p,i,all) => i === 0 || p[0] !== all[i-1][0] || p[1] !== all[i-1][1]);
}
export function connectPads(board, from, to, layer, width) {
  const all=pads(board), a=all.find(p=>p.id===from), b=all.find(p=>p.id===to);
  if (!a || !b || a.id === b.id || !a.net || a.net !== b.net) throw new Error("Choose two different pads on the same assigned net.");
  if (!["top","bottom"].includes(layer) || !finite(width,.1,5)) throw new Error("Choose a copper layer and a width from 0.1 to 5 mm.");
  if (board.routes.length >= 128) throw new Error("This editor supports 128 connections.");
  return { id:"route-"+crypto.randomUUID().slice(0,12), from, to, layer, width, points:routePoints(a,b) };
}
export function moveComponent(board, component, changes) {
  Object.assign(component, changes);
  // Existing copper stays meaningful: moving/rotating/reassigning a footprint
  // removes its connections instead of presenting stretched copper as routed.
  const prefix=component.id+":", count=board.routes.length;
  board.routes=board.routes.filter(r=>!r.from.startsWith(prefix)&&!r.to.startsWith(prefix));
  return count-board.routes.length;
}
const samePoint=(a,b)=>Math.hypot(a[0]-b[0],a[1]-b[1])<1e-9;
function pointSegment(p,a,b) {
  const dx=b[0]-a[0],dy=b[1]-a[1],length=dx*dx+dy*dy;
  const t=length ? Math.max(0,Math.min(1,((p[0]-a[0])*dx+(p[1]-a[1])*dy)/length)) : 0;
  return Math.hypot(p[0]-a[0]-t*dx,p[1]-a[1]-t*dy);
}
export function segmentDistance(a,b,c,d) {
  const cross=(u,v,w)=>(v[0]-u[0])*(w[1]-u[1])-(v[1]-u[1])*(w[0]-u[0]);
  const p=cross(a,b,c),q=cross(a,b,d),r=cross(c,d,a),s=cross(c,d,b);
  if (((p>0&&q<0)||(p<0&&q>0))&&((r>0&&s<0)||(r<0&&s>0))) return 0;
  return Math.min(pointSegment(a,c,d),pointSegment(b,c,d),pointSegment(c,a,b),pointSegment(d,a,b));
}
const segments=route=>route.points.slice(1).map((p,i)=>[route.points[i],p]);
export function pcbChecks(board) {
  const all=pads(board), byID=new Map(all.map(p=>[p.id,p])), issues=[];
  const add=(code,object,message)=>issues.push({code,object,message});
  for(const c of board.components) {
    const b=bodyBounds(c);
    if(b.x<0||b.y<0||b.x+b.width>board.width||b.y+b.height>board.height) add("body-edge",c.id,c.id+" body extends outside the board.");
  }
  for(let i=0;i<all.length;i++) {
    const a=all[i];
    if(!a.net) add("unassigned",a.component,a.id+" has no net.");
    if(Math.min(a.x,a.y,board.width-a.x,board.height-a.y)-a.radius<board.rules.edge-1e-9) add("pad-edge",a.component,a.id+" violates copper-to-edge distance.");
    for(const b of all.slice(i+1)) {
      if(a.net&&a.net===b.net) continue;
      const gap=Math.hypot(a.x-b.x,a.y-b.y)-a.radius-b.radius;
      if(gap<board.rules.clearance-1e-9) add("pad-clearance",a.component,a.id+" / "+b.id+" clearance "+gap.toFixed(3)+" mm.");
    }
  }
  for(const r of board.routes) {
    const a=byID.get(r.from), b=byID.get(r.to);
    if(r.width<board.rules.minWidth-1e-9) add("width",r.id,r.id+" is below minimum track width.");
    if(!samePoint(r.points[0],[a.x,a.y])||!samePoint(r.points.at(-1),[b.x,b.y])) add("endpoint",r.id,r.id+" does not terminate at its assigned pad centres.");
    if(r.points.some(p=>Math.min(p[0],p[1],board.width-p[0],board.height-p[1])-r.width/2<board.rules.edge-1e-9)) add("track-edge",r.id,r.id+" violates copper-to-edge distance.");
    for(const p of all) if(!p.net||p.net!==a.net) {
      const gap=Math.min(...segments(r).map(([u,v])=>pointSegment([p.x,p.y],u,v)))-p.radius-r.width/2;
      if(gap<board.rules.clearance-1e-9) add("track-pad",r.id,r.id+" / "+p.id+" clearance "+gap.toFixed(3)+" mm.");
    }
  }
  for(let i=0;i<board.routes.length;i++) for(const b of board.routes.slice(i+1)) {
    const a=board.routes[i];
    if(a.layer!==b.layer||byID.get(a.from).net===byID.get(b.from).net) continue;
    const gap=Math.min(...segments(a).flatMap(([u,v])=>segments(b).map(([w,z])=>segmentDistance(u,v,w,z))))-(a.width+b.width)/2;
    if(gap<board.rules.clearance-1e-9) add("track-clearance",a.id,a.id+" / "+b.id+" clearance "+gap.toFixed(3)+" mm on "+a.layer+".");
  }
  // Only explicit, centre-terminated pad connections establish connectivity.
  // Touching/intersecting same-net copper is deliberately not inferred.
  const nets=[...new Set(all.map(p=>p.net).filter(Boolean))].sort().map(name=>{
    const members=all.filter(p=>p.net===name), parent=new Map(members.map(p=>[p.id,p.id]));
    const root=k=>{while(parent.get(k)!==k)k=parent.get(k);return k;};
    for(const r of board.routes) if(byID.get(r.from).net===name) {
      const a=byID.get(r.from),b=byID.get(r.to);
      if(samePoint(r.points[0],[a.x,a.y])&&samePoint(r.points.at(-1),[b.x,b.y])) parent.set(root(a.id),root(b.id));
    }
    const groups=[...new Set(members.map(p=>root(p.id)))];
    const remaining=Math.max(0,groups.length-1), airwires=[];
    for(let i=1;i<groups.length;i++) airwires.push([members.find(p=>root(p.id)===groups[0]),members.find(p=>root(p.id)===groups[i])]);
    return {name,pads:members.length,remaining,airwires};
  });
  return { issues, nets, unrouted:nets.reduce((n,net)=>n+net.remaining,0) };
}
export function siliconChecks(chip) {
  const issues=[], add=(code,object,message)=>issues.push({code,object,message});
  for(const a of chip.blocks) {
    if(a.x<chip.margin||a.y<chip.margin||a.x+a.width>chip.width-chip.margin||a.y+a.height>chip.height-chip.margin) add("boundary",a.id,a.name+" crosses the core keepout.");
  }
  for(let i=0;i<chip.blocks.length;i++) for(const b of chip.blocks.slice(i+1)) {
    const a=chip.blocks[i], dx=Math.max(a.x-b.x-b.width,b.x-a.x-a.width), dy=Math.max(a.y-b.y-b.height,b.y-a.y-a.height);
    if(dx<0&&dy<0) add("overlap",a.id,a.name+" overlaps "+b.name+".");
    else if(dx<chip.spacing&&dy<chip.spacing) add("spacing",a.id,a.name+" / "+b.name+" violate the rectangular placement halo.");
  }
  const coreArea=(chip.width-2*chip.margin)*(chip.height-2*chip.margin);
  const area=chip.blocks.reduce((sum,b)=>sum+b.width*b.height,0), utilization=100*area/coreArea;
  if(utilization>chip.target) add("density","core","Requested block area exceeds the "+chip.target+"% core-area target.");
  return {issues,area,coreArea,utilization};
}
export function validateDocument(doc) {
  if(!doc||doc.format!==STUDIO_FORMAT||!["pcb","silicon"].includes(doc.kind)) throw new Error("Choose an Emergence PCB or silicon design file.");
  const d=doc.design;
  if(!d||!finite(d.width,10,doc.kind==="pcb"?300:10000)||!finite(d.height,10,doc.kind==="pcb"?300:10000)||!finite(d.grid,.1,100)) throw new Error("Invalid design dimensions or grid.");
  if(doc.kind==="pcb") {
    if(!d.rules||!finite(d.rules.clearance,.05,5)||!finite(d.rules.minWidth,.1,5)||!finite(d.rules.edge,0,10)) throw new Error("Invalid board rules.");
    if(!Array.isArray(d.components)||d.components.length>64||!unique(d.components.map(c=>c?.id))||!d.components.every(c=>c&&id(c.id)&&text(c.name)&&Object.hasOwn(FOOTPRINTS,c.footprint)&&finite(c.x,-300,600)&&finite(c.y,-300,600)&&[0,90,180,270].includes(c.rotation)&&Array.isArray(c.nets)&&c.nets.length===FOOTPRINTS[c.footprint].pads.length&&c.nets.every(n=>typeof n==="string"&&n.length<=32&&/^[a-zA-Z0-9_+./-]*$/.test(n)))) throw new Error("Invalid components, footprints or pad nets.");
    const ps=new Map(pads(d).map(p=>[p.id,p]));
    if(!Array.isArray(d.routes)||d.routes.length>128||!unique(d.routes.map(r=>r?.id))||!d.routes.every(r=>r&&id(r.id)&&ps.has(r.from)&&ps.has(r.to)&&r.from!==r.to&&ps.get(r.from).net&&ps.get(r.from).net===ps.get(r.to).net&&["top","bottom"].includes(r.layer)&&finite(r.width,.1,5)&&Array.isArray(r.points)&&r.points.length>=2&&r.points.length<=32&&r.points.every(point))) throw new Error("Invalid copper connections.");
  } else {
    if(!finite(d.margin,0,1000)||2*d.margin>=Math.min(d.width,d.height)||!finite(d.spacing,0,1000)||!finite(d.target,1,100)) throw new Error("Invalid core constraints.");
    if(!Array.isArray(d.blocks)||d.blocks.length>64||!unique(d.blocks.map(b=>b?.id))||!d.blocks.every(b=>b&&id(b.id)&&text(b.name)&&["logic","memory","io","analog"].includes(b.kind)&&finite(b.x,-10000,20000)&&finite(b.y,-10000,20000)&&finite(b.width,1,10000)&&finite(b.height,1,10000)&&typeof b.fixed==="boolean"&&(b.power===undefined||finite(b.power,0,1e6)))) throw new Error("Invalid floorplan blocks.");
  }
  // Project known fields only; imported JSON is never merged into UI objects.
  const design=doc.kind==="pcb" ? {
    width:d.width,height:d.height,grid:d.grid,rules:{clearance:d.rules.clearance,minWidth:d.rules.minWidth,edge:d.rules.edge,density:pcbExtras(d).density},
    components:d.components.map(c=>({id:c.id,name:c.name,footprint:c.footprint,x:c.x,y:c.y,rotation:c.rotation,nets:[...c.nets]})),
    routes:d.routes.map(r=>({id:r.id,from:r.from,to:r.to,layer:r.layer,width:r.width,points:r.points.map(p=>[...p])})),
    ...pcbExtras(d,true),
  } : {
    width:d.width,height:d.height,grid:d.grid,margin:d.margin,spacing:d.spacing,target:d.target,
    blocks:d.blocks.map(b=>({id:b.id,name:b.name,kind:b.kind,x:b.x,y:b.y,width:b.width,height:b.height,fixed:b.fixed,power:b.power??BLOCK_POWER[b.kind]})),
    chip:chipSettings(d),
  };
  return {format:STUDIO_FORMAT,kind:doc.kind,design};
}
// Stackup, RMS net currents, component dissipation and linked load nets.
function pcbExtras(d, full) {
  const s=d.stackup??{thickness:1.6,copper:1,epsr:4.4}, nets=new Set(d.components.flatMap(c=>c.nets).filter(Boolean));
  const density=d.rules.density??35;
  if(!finite(density,1,1000)) throw new Error("Invalid board rules.");
  if(!full) return {density};
  if(!s||!finite(s.thickness,.2,6)||!finite(s.copper,.25,6)||!finite(s.epsr,1,20)) throw new Error("Invalid board stackup.");
  const currents=d.currents??{}, powers=d.powers??{}, loadNets=d.loadNets??[];
  if(typeof currents!=="object"||Array.isArray(currents)||!Object.entries(currents).every(([n,a])=>/^[a-zA-Z0-9_+./-]{1,32}$/.test(n)&&finite(a,0,500))) throw new Error("Invalid net currents.");
  if(typeof powers!=="object"||Array.isArray(powers)||!Object.entries(powers).every(([c,w])=>d.components.some(x=>x.id===c)&&finite(w,0,500))) throw new Error("Invalid component dissipation.");
  if(!Array.isArray(loadNets)||loadNets.length>8||!loadNets.every(n=>typeof n==="string"&&n.length<=32)) throw new Error("Invalid load nets.");
  return {stackup:{thickness:s.thickness,copper:s.copper,epsr:s.epsr},
    currents:Object.fromEntries(Object.entries(currents).filter(([n])=>nets.has(n)).sort()),
    powers:Object.fromEntries(Object.entries(powers).sort()),loadNets:[...new Set(loadNets)]};
}
export const RTL_SOURCES = { mini8:"mini8: 8-bit accumulator CPU (emergence, Apache-2.0)", "tt-motor":"tt_um_hadjal001_motor_controller: dual PWM (Tiny Tapeout, Apache-2.0)", custom:"Custom Verilog" };
function chipSettings(d) {
  const c=d.chip??{rtl:[],frequency:100,vdd:1.2,activity:20,instances:1,theta:45,ambient:25};
  // Earlier documents named one block and an optional Verilog override.
  const rtl=Array.isArray(c.rtl)?c.rtl:c.block?[{block:c.block,source:c.verilog?"custom":"mini8",verilog:c.verilog??null}]:[];
  if(!c||!finite(c.frequency,.001,5000)||!finite(c.vdd,1.08,1.32)||!finite(c.activity,0,100)||!Number.isInteger(c.instances)||c.instances<1||c.instances>100000||!finite(c.theta,.1,500)||!finite(c.ambient,-40,125)
    ||rtl.length>4||!unique(rtl.map(r=>r?.block))||!rtl.every(r=>r&&d.blocks.some(b=>b.id===r.block)&&Object.hasOwn(RTL_SOURCES,r.source)&&(r.source==="custom"?typeof r.verilog==="string"&&r.verilog.length>0&&r.verilog.length<=65536:r.verilog===null)))
    throw new Error("Invalid chip flow settings.");
  return {rtl:rtl.map(r=>({block:r.block,source:r.source,verilog:r.verilog})),frequency:c.frequency,vdd:c.vdd,activity:c.activity,instances:c.instances,theta:c.theta,ambient:c.ambient};
}
export function documentFor(kind, design) { return validateDocument({format:STUDIO_FORMAT,kind,design}); }
const csv=value=>'"'+String(value).replaceAll('"','""')+'"';
export function pcbBOM(board) {
  return "Reference,Description,Illustrative footprint,X (mm),Y (mm),Rotation (deg)\r\n"+
    board.components.map(c=>[c.id,c.name,FOOTPRINTS[c.footprint].title,c.x,c.y,c.rotation].map(csv).join(",")).join("\r\n")+"\r\n";
}
