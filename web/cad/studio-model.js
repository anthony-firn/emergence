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
// A pad deliberately left unconnected: not a net, and never joined to another.
export const NO_CONNECT = "NC";

export const FOOTPRINTS = {
  connector: { title: "2-pin connector", width: 8, height: 7, pads: [[-2, 0], [2, 0]] },
  resistor: { title: "Axial resistor", width: 9, height: 3, pads: [[-5, 0], [5, 0]] },
  capacitor: { title: "Radial capacitor", width: 6, height: 6, pads: [[-1.5, 0], [1.5, 0]] },
  transistor: { title: "3-pin power device", width: 8, height: 6, pads: [[-3, 3], [0, 3], [3, 3]] },
  dip8: { title: "DIP-8", width: 8, height: 12, pads: [[-4,-4.5],[-4,-1.5],[-4,1.5],[-4,4.5],[4,4.5],[4,1.5],[4,-1.5],[4,-4.5]] },
  via: { title: "Via", width: 1, height: 1, pad: .8, drill: .4, mass: 0, pads: [[0, 0]] },
};
export function pcbExample() {
  return {
    width: 80, height: 55, grid: .5,
    rules: { clearance: .25, minWidth: .3, edge: 1 }, footprints: {}, pours: [],
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
// A board's own footprints sit beside the built-in ones.
export const footprintOf = (board, name) => board?.footprints && Object.hasOwn(board.footprints, name) ? board.footprints[name] : Object.hasOwn(FOOTPRINTS, name) ? FOOTPRINTS[name] : undefined;
// Through-hole pad patterns from a package's numbers: a single row, two rows
// numbered counter-clockwise as on a DIP, or four sides.
export function footprintPads(kind, pins, pitch, span) {
  if (kind === "single") return Array.from({ length: pins }, (_, i) => [(i - (pins - 1) / 2) * pitch, 0]);
  if (kind === "dual") { if (pins % 2) throw new Error("A dual-row footprint has an even pin count."); const n = pins / 2, y = i => (i - (n - 1) / 2) * pitch; return [...Array.from({ length: n }, (_, i) => [-span / 2, y(i)]), ...Array.from({ length: n }, (_, i) => [span / 2, y(n - 1 - i)])]; }
  if (kind === "quad") { if (pins % 4) throw new Error("A four-sided footprint has a pin count divisible by four."); const n = pins / 4, s = i => (i - (n - 1) / 2) * pitch; return [...Array.from({ length: n }, (_, i) => [-span / 2, s(i)]), ...Array.from({ length: n }, (_, i) => [s(i), span / 2]), ...Array.from({ length: n }, (_, i) => [span / 2, s(n - 1 - i)]), ...Array.from({ length: n }, (_, i) => [s(n - 1 - i), -span / 2])]; }
  throw new Error("Footprint kinds: single, dual, quad.");
}
export function pads(board) {
  return board.components.flatMap(c => {
    const angle = c.rotation * Math.PI / 180, cs = Math.round(Math.cos(angle)), sn = Math.round(Math.sin(angle));
    const f = footprintOf(board, c.footprint), turned = c.rotation % 180 !== 0;
    return f.pads.map(([x,y,w,h], i) => {
      // A surface pad is a rectangle on the top layer; a pad may carry its own size.
      const pw = w ?? f.padW, ph = h ?? f.padH, rect = f.smd && pw !== undefined;
      return { id:c.id + ":" + (i+1), component:c.id, number:i+1, net:c.nets[i],
        x:c.x+x*cs-y*sn, y:c.y+x*sn+y*cs, radius:rect ? Math.hypot(pw,ph)/2 : (f.pad ?? 1.6)/2, ...(rect && {w:turned?ph:pw,h:turned?pw:ph,smd:true}), ...(f.drill !== undefined && {drill:f.drill}) };
    });
  });
}
export function bodyBounds(c, board) {
  const f = footprintOf(board, c.footprint), rotated = c.rotation % 180 !== 0;
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
  if (!a || !b || a.id === b.id || !a.net || a.net === NO_CONNECT || a.net !== b.net) throw new Error("Choose two different pads on the same assigned net.");
  if (!["top","bottom"].includes(layer) || !finite(width,.1,5)) throw new Error("Choose a copper layer and a width from 0.1 to 5 mm.");
  if (board.routes.length >= 512) throw new Error("This editor supports 512 connections.");
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
const halves=p=>p.w?[p.w/2,p.h/2]:[p.radius,p.radius];
// Clearance between pads: circles exactly, otherwise as rectangles.
export function padGap(a,b){
  if(!a.w&&!b.w)return Math.hypot(a.x-b.x,a.y-b.y)-a.radius-b.radius;
  const [ax,ay]=halves(a),[bx,by]=halves(b),dx=Math.abs(a.x-b.x)-ax-bx,dy=Math.abs(a.y-b.y)-ay-by;
  return dx>0&&dy>0?Math.hypot(dx,dy):Math.max(dx,dy);
}
// Distance from a point to a pad's copper.
export function padDistance(p,x,y){
  if(!p.w)return Math.hypot(x-p.x,y-p.y)-p.radius;
  const dx=Math.abs(x-p.x)-p.w/2,dy=Math.abs(y-p.y)-p.h/2;
  return dx>0&&dy>0?Math.hypot(dx,dy):Math.max(dx,dy);
}
// Even-odd point in polygon, for pours.
export function inPolygon(points,x,y){
  let k=false;for(let i=0,j=points.length-1;i<points.length;j=i++){const [xi,yi]=points[i],[xj,yj]=points[j];if((yi>y)!==(yj>y)&&x<(xj-xi)*(y-yi)/(yj-yi)+xi)k=!k;}return k;
}
const padOnLayer=(p,layer)=>!p.smd||layer==="top";
// The copper a pour really fills: its outline less clearance to other nets, as
// connected islands on a raster. A trace across a pour can cut it in two, and only
// the pads on one island are joined by it.
export function pourIslands(board,pour,all=pads(board)) {
  const xs=pour.points.map(p=>p[0]),ys=pour.points.map(p=>p[1]),x0=Math.min(...xs),y0=Math.min(...ys),w=Math.max(...xs)-x0,h=Math.max(...ys)-y0;
  const g=Math.max(.1,Math.min(.25,Math.sqrt(w*h)/240)),nx=Math.ceil(w/g),ny=Math.ceil(h/g),X=i=>x0+(i+.5)*g,Y=j=>y0+(j+.5)*g,cells=new Int32Array(nx*ny),c=board.rules.clearance;
  for(let j=0;j<ny;j++)for(let i=0;i<nx;i++)if(inPolygon(pour.points,X(i),Y(j)))cells[j*nx+i]=-1;
  const clear=(ax,ay,bx,by,test)=>{for(let j=Math.max(0,Math.floor((ay-y0)/g-.5));j<=Math.min(ny-1,Math.ceil((by-y0)/g-.5));j++)for(let i=Math.max(0,Math.floor((ax-x0)/g-.5));i<=Math.min(nx-1,Math.ceil((bx-x0)/g-.5));i++)if(cells[j*nx+i]&&test(X(i),Y(j)))cells[j*nx+i]=0;};
  for(const p of all)if(p.net!==pour.net&&padOnLayer(p,pour.layer)&&p.x+p.radius+c>=x0&&p.x-p.radius-c<=x0+w&&p.y+p.radius+c>=y0&&p.y-p.radius-c<=y0+h){const r=p.radius+c;clear(p.x-r,p.y-r,p.x+r,p.y+r,(x,y)=>padDistance(p,x,y)<c);}
  const byID=new Map(all.map(p=>[p.id,p]));
  for(const r of board.routes??[])if(r.layer===pour.layer&&byID.get(r.from)?.net!==pour.net)for(const [u,v] of segments(r)){const m=r.width/2+c;clear(Math.min(u[0],v[0])-m,Math.min(u[1],v[1])-m,Math.max(u[0],v[0])+m,Math.max(u[1],v[1])+m,(x,y)=>pointSegment([x,y],u,v)<m);}
  let islands=0;const stack=[];
  for(let k=0;k<cells.length;k++)if(cells[k]===-1){cells[k]=++islands;stack.push(k);while(stack.length){const n=stack.pop(),i=n%nx;for(const q of [i>0?n-1:-1,i<nx-1?n+1:-1,n-nx,n+nx])if(q>=0&&q<cells.length&&cells[q]===-1){cells[q]=islands;stack.push(q);}}}
  // Every island a pad's copper touches.
  const under=p=>{const r=p.radius+g,found=new Set();for(let j=Math.max(0,Math.floor((p.y-r-y0)/g));j<=Math.min(ny-1,Math.ceil((p.y+r-y0)/g));j++)for(let i=Math.max(0,Math.floor((p.x-r-x0)/g));i<=Math.min(nx-1,Math.ceil((p.x+r-x0)/g));i++){const v=cells[j*nx+i];if(v&&padDistance(p,X(i),Y(j))<=g/2)found.add(v);}return [...found];};
  return {islands,under};
}
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
function segments(route){return route.points.slice(1).map((p,i)=>[route.points[i],p]);}
function bounds(route){const xs=route.points.map(p=>p[0]),ys=route.points.map(p=>p[1]);return {x0:Math.min(...xs),x1:Math.max(...xs),y0:Math.min(...ys),y1:Math.max(...ys)};}
export function pcbChecks(board) {
  const all=pads(board), byID=new Map(all.map(p=>[p.id,p])), issues=[];
  const add=(code,object,message)=>issues.push({code,object,message});
  for(const c of board.components) {
    const b=bodyBounds(c,board);
    if(b.x<0||b.y<0||b.x+b.width>board.width||b.y+b.height>board.height) add("body-edge",c.id,c.id+" body extends outside the board.");
  }
  // Bodies are courtyards: parts cannot be assembled closer than this. Vias sit anywhere.
  const parts=board.components.filter(c=>c.footprint!=="via").map(c=>[c,bodyBounds(c,board)]);
  for(let i=0;i<parts.length;i++)for(const [c,b] of parts.slice(i+1)){const [a,u]=parts[i];if(u.x<b.x+b.width-1e-6&&b.x<u.x+u.width-1e-6&&u.y<b.y+b.height-1e-6&&b.y<u.y+u.height-1e-6)add("courtyard",a.id,a.id+" and "+c.id+" overlap.");}
  for(let i=0;i<all.length;i++) {
    const a=all[i];
    if(!a.net) add("unassigned",a.component,a.id+" has no net.");
    if(Math.min(a.x,a.y,board.width-a.x,board.height-a.y)-a.radius<board.rules.edge-1e-9) add("pad-edge",a.component,a.id+" violates copper-to-edge distance.");
    for(const b of all.slice(i+1)) {
      if(a.net&&a.net!==NO_CONNECT&&a.net===b.net) continue;
      const gap=padGap(a,b);
      if(gap<board.rules.clearance-1e-9) add("pad-clearance",a.component,a.id+" / "+b.id+" clearance "+gap.toFixed(3)+" mm.");
    }
  }
  for(const r of board.routes) {
    const a=byID.get(r.from), b=byID.get(r.to);
    if(r.width<board.rules.minWidth-1e-9) add("width",r.id,r.id+" is below minimum track width.");
    if((a.smd||b.smd)&&r.layer!=="top") add("smd-layer",r.id,r.id+" reaches a surface pad from the bottom layer; use a via.");
    if(!samePoint(r.points[0],[a.x,a.y])||!samePoint(r.points.at(-1),[b.x,b.y])) add("endpoint",r.id,r.id+" does not terminate at its assigned pad centres.");
    if(r.points.some(p=>Math.min(p[0],p[1],board.width-p[0],board.height-p[1])-r.width/2<board.rules.edge-1e-9)) add("track-edge",r.id,r.id+" violates copper-to-edge distance.");
    const box=bounds(r),reach=r.width/2+board.rules.clearance;
    for(const p of all) if((!p.net||p.net!==a.net)&&padOnLayer(p,r.layer)&&p.x+p.radius+reach>=box.x0&&p.x-p.radius-reach<=box.x1&&p.y+p.radius+reach>=box.y0&&p.y-p.radius-reach<=box.y1) {
      const gap=Math.min(...segments(r).map(([u,v])=>p.w?Math.min(...Array.from({length:25},(_,k)=>padDistance(p,u[0]+(v[0]-u[0])*k/24,u[1]+(v[1]-u[1])*k/24))):pointSegment([p.x,p.y],u,v)-p.radius))-r.width/2;
      if(gap<board.rules.clearance-1e-9) add("track-pad",r.id,r.id+" / "+p.id+" clearance "+gap.toFixed(3)+" mm.");
    }
  }
  const boxes=new Map(board.routes.map(r=>[r,bounds(r)]));
  for(let i=0;i<board.routes.length;i++) for(const b of board.routes.slice(i+1)) {
    const a=board.routes[i],u=boxes.get(a),v=boxes.get(b),far=(a.width+b.width)/2+board.rules.clearance;
    if(a.layer!==b.layer||byID.get(a.from).net===byID.get(b.from).net||u.x0-far>v.x1||v.x0-far>u.x1||u.y0-far>v.y1||v.y0-far>u.y1) continue;
    const gap=Math.min(...segments(a).flatMap(([u,v])=>segments(b).map(([w,z])=>segmentDistance(u,v,w,z))))-(a.width+b.width)/2;
    if(gap<board.rules.clearance-1e-9) add("track-clearance",a.id,a.id+" / "+b.id+" clearance "+gap.toFixed(3)+" mm on "+a.layer+".");
  }
  // Only explicit, centre-terminated pad connections establish connectivity.
  // Touching/intersecting same-net copper is deliberately not inferred.
  const nets=[...new Set(all.map(p=>p.net).filter(n=>n&&n!==NO_CONNECT))].sort().map(name=>{
    const members=all.filter(p=>p.net===name), parent=new Map(members.map(p=>[p.id,p.id]));
    const root=k=>{while(parent.get(k)!==k)k=parent.get(k);return k;};
    for(const r of board.routes) if(byID.get(r.from).net===name) {
      const a=byID.get(r.from),b=byID.get(r.to);
      if(samePoint(r.points[0],[a.x,a.y])&&samePoint(r.points.at(-1),[b.x,b.y])) parent.set(root(a.id),root(b.id));
    }
    // A pour joins the pads of its net on each island of copper it really fills.
    const cut=[];
    for(const pour of board.pours??[]) if(pour.net===name) {
      const {under}=pourIslands(board,pour,all),on=new Map();
      for(const p of members)if(padOnLayer(p,pour.layer))for(const island of under(p)){if(on.has(island))parent.set(root(p.id),root(on.get(island).id));else on.set(island,p);}
      if(on.size>1)cut.push([pour,[...on.values()]]);
    }
    for(const [pour,first] of cut)if(new Set(first.map(p=>root(p.id))).size>1)add("pour-split",pour.id,pour.id+" is cut into separate islands, so "+name+" is open across it. A trace crosses the pour.");
    const groups=[...new Set(members.map(p=>root(p.id)))].map(g=>members.filter(p=>root(p.id)===g));
    // Each open group connects to the nearest copper of its net, as a ratsnest does.
    const remaining=Math.max(0,groups.length-1), airwires=[], linked=groups.slice(0,1), rest=groups.slice(1);
    while(rest.length){let best=null;for(const a of linked.flat())for(const [i,set] of rest.entries())for(const b of set){const d=Math.hypot(a.x-b.x,a.y-b.y);if(!best||d<best.d)best={a,b,i,d};}airwires.push([best.a,best.b]);linked.push(...rest.splice(best.i,1));}
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
    for(const [name,f] of Object.entries(d.footprints??{})) if(!id(name)||Object.hasOwn(FOOTPRINTS,name)||!f||!text(f.title)||!finite(f.width,.5,300)||!finite(f.height,.5,300)||!finite(f.mass??1,0,1000)||!finite(f.pad??1.6,.2,20)||!finite(f.drill??.9,.1,20)||(f.drill??.9)>=(f.pad??1.6)||!Array.isArray(f.pads)||!f.pads.length||f.pads.length>256||!f.pads.every(p=>Array.isArray(p)&&(p.length===2||p.length===4)&&p.every(n=>finite(n,-1000,1000)))||f.smd!==undefined&&typeof f.smd!=="boolean"||f.smd&&!f.pads.every(p=>p.length===4||finite(f.padW,.05,100)&&finite(f.padH,.05,100))) throw new Error(`Invalid footprint ${name}.`);
    if(Object.keys(d.footprints??{}).length>32) throw new Error("A board holds 32 footprints of its own.");
    if(!Array.isArray(d.components)||d.components.length>256||!unique(d.components.map(c=>c?.id))||!d.components.every(c=>c&&id(c.id)&&text(c.name)&&footprintOf(d,c.footprint)&&finite(c.x,-300,600)&&finite(c.y,-300,600)&&[0,90,180,270].includes(c.rotation)&&Array.isArray(c.nets)&&c.nets.length===footprintOf(d,c.footprint).pads.length&&c.nets.every(n=>typeof n==="string"&&n.length<=32&&/^[a-zA-Z0-9_+./-]*$/.test(n)))) throw new Error("Invalid components, footprints or pad nets.");
    const ps=new Map(pads(d).map(p=>[p.id,p]));
    if(d.pours!==undefined&&(!Array.isArray(d.pours)||d.pours.length>32||!unique(d.pours.map(p=>p?.id))||!d.pours.every(p=>p&&id(p.id)&&typeof p.net==="string"&&p.net.length<=32&&["top","bottom"].includes(p.layer)&&Array.isArray(p.points)&&p.points.length>=3&&p.points.length<=64&&p.points.every(point)))) throw new Error("Invalid copper pours.");
    if(!Array.isArray(d.routes)||d.routes.length>512||!unique(d.routes.map(r=>r?.id))||!d.routes.every(r=>r&&id(r.id)&&ps.has(r.from)&&ps.has(r.to)&&r.from!==r.to&&ps.get(r.from).net&&ps.get(r.from).net===ps.get(r.to).net&&["top","bottom"].includes(r.layer)&&finite(r.width,.1,5)&&Array.isArray(r.points)&&r.points.length>=2&&r.points.length<=32&&r.points.every(point))) throw new Error("Invalid copper connections.");
  } else {
    if(!finite(d.margin,0,1000)||2*d.margin>=Math.min(d.width,d.height)||!finite(d.spacing,0,1000)||!finite(d.target,1,100)) throw new Error("Invalid core constraints.");
    if(!Array.isArray(d.blocks)||d.blocks.length>64||!unique(d.blocks.map(b=>b?.id))||!d.blocks.every(b=>b&&id(b.id)&&text(b.name)&&["logic","memory","io","analog"].includes(b.kind)&&finite(b.x,-10000,20000)&&finite(b.y,-10000,20000)&&finite(b.width,1,10000)&&finite(b.height,1,10000)&&typeof b.fixed==="boolean"&&(b.power===undefined||finite(b.power,0,1e6)))) throw new Error("Invalid floorplan blocks.");
  }
  // Project known fields only; imported JSON is never merged into UI objects.
  const design=doc.kind==="pcb" ? {
    width:d.width,height:d.height,grid:d.grid,rules:{clearance:d.rules.clearance,minWidth:d.rules.minWidth,edge:d.rules.edge,density:pcbExtras(d).density},
    footprints:Object.fromEntries(Object.entries(d.footprints??{}).map(([k,f])=>[k,{title:f.title,width:f.width,height:f.height,mass:f.mass??1,pad:f.pad??1.6,drill:f.drill??.9,...(f.smd&&{smd:true}),...(f.padW!==undefined&&{padW:f.padW,padH:f.padH}),pads:f.pads.map(p=>[...p])}])),
    pours:(d.pours??[]).map(p=>({id:p.id,net:p.net,layer:p.layer,points:p.points.map(q=>[...q])})),
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
    board.components.map(c=>[c.id,c.name,footprintOf(board,c.footprint).title,c.x,c.y,c.rotation].map(csv).join(",")).join("\r\n")+"\r\n";
}
