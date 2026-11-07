import { footprintOf } from "./studio-model.js";
// Small orthographic canvas renderer. Geometry uses the same fixed dimensions
// as the Rust adapter. The exported mesh comes from cad::Solid, not this view.
const faces = [[0,3,2,1],[4,5,6,7],[0,1,5,4],[2,3,7,6],[1,2,6,5],[3,0,4,7]];
function box(x, y, z, w, d, h, color, group, id) {
  return { points: [[x-w/2,y-d/2,z],[x+w/2,y-d/2,z],[x+w/2,y+d/2,z],[x-w/2,y+d/2,z],
    [x-w/2,y-d/2,z+h],[x+w/2,y-d/2,z+h],[x+w/2,y+d/2,z+h],[x-w/2,y+d/2,z+h]], color, group, id };
}
function cylinder(f,color) {
  const points=[], polygons=[], count=24;
  for(const z of [f.z,f.z+f.height])for(let j=0;j<count;j++)points.push([f.x+Math.cos(j/count*2*Math.PI)*f.width/2,f.y+Math.sin(j/count*2*Math.PI)*f.width/2,z]);
  polygons.push(Array.from({length:count},(_,j)=>count-1-j),Array.from({length:count},(_,j)=>count+j));
  for(let j=0;j<count;j++)polygons.push([j,(j+1)%count,(j+1)%count+count,j+count]);
  return {points,faces:polygons,color,group:f.domain,id:f.id};
}
// The PCB and heat sink share the project's scene placement. Show one ESC
// instance in the heat sink's local frame; phase cancels when both repeat together.
const radians = angle => (angle || 0) * Math.PI / 180;
const rz = ([x,y,z], a) => [x*Math.cos(a)-y*Math.sin(a), x*Math.sin(a)+y*Math.cos(a), z];
const rx = ([x,y,z], a) => [x, y*Math.cos(a)-z*Math.sin(a), y*Math.sin(a)+z*Math.cos(a)];
function placed(point, p) {
  const rotated = rz(rx(point, radians(p.rx)), radians(p.rz));
  return rz(rotated.map((v,i)=>v+p.at[i]), radians(p.phase));
}
function unplaced(point, p) {
  const translated = rz(point, -radians(p.phase)).map((v,i)=>v-p.at[i]);
  return rx(rz(translated, -radians(p.rz)), -radians(p.rx));
}
const componentId = c => c.id === "Q1" ? "q1" : `pcb-${c.id}`;
export function assemblyParts(params,features=[],context=null) {
  const parts = [{id:"sink",name:"Finned heat sink",domain:"Mechanical",shape:"sink"}];
  if (context?.board) parts.push({id:"board",name:"PCB (live board dimensions)",domain:"Electrical",shape:"context"},
    ...context.board.components.filter(c=>c.footprint!=="via").map(c=>({id:componentId(c),name:`${c.id} · ${c.name}`,domain:"Electrical",shape:"context"})));
  else parts.push({id:"board",name:"Power board (schematic)",domain:"Electrical",shape:"context"},
    {id:"q1",name:"MOSFET Q1 (schematic)",domain:"Electrical",shape:"context"},{id:"optic",name:"Optical marker (schematic)",domain:"Photonics",shape:"context"});
  return [...parts,...features];
}
export function assemblyGeometry(params, features=[], context=null, exploded=false) {
  const {width:w,length:l,base:b,finHeight:h,finThickness:t}=params, e=exploded?15:0;
  const boxes=[box(0,0,0,w,l,b,[129,150,164],"Mechanical","sink")];
  for(let j=0;j<params.fins;j++) boxes.push(box(-(w-t)/2+(w-t)*j/(params.fins-1),0,b,t,l,h,[163,185,200],"Mechanical","sink"));
  if(context?.board) {
    const board=context.board, z=board.stackup.thickness;
    const boardBoxes=[box(0,0,0,board.width,board.height,z,[36,93,82],"Electrical","board")];
    for(const c of board.components) {
      if(c.footprint==="via") continue;
      const f=footprintOf(board,c.footprint), turned=c.rotation%180!==0;
      boardBoxes.push(box(c.x-board.width/2,board.height/2-c.y,z,turned?f.height:f.width,turned?f.width:f.height,f.smd?1.2:3,
        c.id.startsWith("Q")?[200,116,62]:[62,70,77],"Electrical",componentId(c)));
    }
    for(const item of boardBoxes) {
      item.points=item.points.map(point=> {
        const local=context.pcbPlacement && context.mechPlacement ? unplaced(placed(point,context.pcbPlacement),context.mechPlacement) : [point[0],point[1],point[2]-z-1.4];
        // Explode the populated PCB as one unit; devices stay soldered to it.
        local[2]-=e; return local;
      });
      boxes.push(item);
    }
  } else {
    // Legacy stand-alone thermal study: the package sits on the board,
    // below the sink base, never above its fin tips.
    boxes.push(box(0,0,-4-e,w+20,l+22,1.6,[36,93,82],"Electrical","board"),
      box(0,0,-2.4-e,10,15,2,[200,116,62],"Electrical","q1"),
      box(-w/2+1,-l/2,-2+e,8,10,5,[61,142,149],"Photonics","optic"));
  }
  const colors={Mechanical:[139,172,194],Electrical:[211,153,65],Thermal:[208,116,92],Photonics:[148,119,196],Radiation:[104,176,140]};
  for(const f of features) boxes.push(f.shape==="box"?box(f.x,f.y,f.z,f.width,f.depth,f.height,colors[f.domain],f.domain,f.id):cylinder(f,colors[f.domain]));
  return boxes;
}
export class Scene {
  constructor(canvas) {
    this.canvas = canvas; this.ctx = canvas.getContext("2d"); this.yaw = -.6; this.pitch = .52; this.zoom = 1; this.params = { fins: 9 }; this.domain = "Mechanical";
    this.resize = new ResizeObserver(() => this.draw()); this.resize.observe(canvas);
    this.features=[];this.selected="sink";this.hidden=new Set();
    canvas.addEventListener("pointerdown", e => { this.drag = [e.clientX, e.clientY];this.origin=[e.clientX,e.clientY]; canvas.setPointerCapture(e.pointerId); });
    canvas.addEventListener("pointermove", e => { if (!this.drag) return;
      this.yaw += (e.clientX - this.drag[0]) * .008; this.pitch = Math.max(0, Math.min(Math.PI/2, this.pitch + (e.clientY - this.drag[1]) * .008)); this.drag = [e.clientX, e.clientY]; this.draw(); });
    canvas.addEventListener("pointerup", e => { if(this.origin&&Math.hypot(e.clientX-this.origin[0],e.clientY-this.origin[1])<4){const r=canvas.getBoundingClientRect();this.pick(e.clientX-r.left,e.clientY-r.top);}this.drag = null;this.origin=null; }); canvas.addEventListener("pointercancel", () => this.drag = null);
    canvas.addEventListener("wheel", e => { e.preventDefault(); this.zoom = Math.max(.55, Math.min(2, this.zoom * Math.exp(-e.deltaY * .001))); this.draw(); }, { passive: false });
    canvas.addEventListener("keydown", e => {
      if (!["ArrowLeft","ArrowRight","ArrowUp","ArrowDown","+","=","-","0"].includes(e.key)) return;
      e.preventDefault();
      if (e.key === "ArrowLeft") this.yaw -= .1; if (e.key === "ArrowRight") this.yaw += .1;
      if (e.key === "ArrowUp") this.pitch = Math.min(Math.PI/2,this.pitch+.1); if (e.key === "ArrowDown") this.pitch = Math.max(0,this.pitch-.1);
      if (["+","="].includes(e.key)) this.zoom = Math.min(2,this.zoom+.1); if (e.key === "-") this.zoom = Math.max(.55,this.zoom-.1);
      if (e.key === "0") this.fit(); this.draw();
    });
    this.theme = new MutationObserver(() => this.draw()); this.theme.observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
    this.scheme = matchMedia("(prefers-color-scheme: light)"); this.scheme.addEventListener("change", () => this.draw());
  }
  fit() { this.yaw = -.6; this.pitch = .52; this.zoom = 1; this.draw(); }
  standardView(name) {
    const views={iso:[-.6,.52],top:[0,Math.PI/2],front:[0,0],right:[Math.PI/2,0]};
    if(!views[name])return;
    [this.yaw,this.pitch]=views[name];this.draw();
  }
  update(params, domain, exploded, features=[], selected="sink", hidden=new Set(), context=null) { this.context=context; this.params = params; this.domain = domain; this.exploded = exploded;this.features=features;this.selected=selected;this.hidden=hidden; this.draw(); }
  pick(x,y) {
    for(const polygon of [...(this.polys||[])].reverse()){
      let inside=false;const ps=polygon.vertices;
      for(let i=0,j=ps.length-1;i<ps.length;j=i++)if((ps[i][1]>y)!==(ps[j][1]>y)&&x<(ps[j][0]-ps[i][0])*(y-ps[i][1])/(ps[j][1]-ps[i][1])+ps[i][0])inside=!inside;
      if(inside){this.onselect?.(polygon.id);return;}
    }
  }
  project([x, y, z]) {
    const a = this.yaw, p = this.pitch, rx = x*Math.cos(a)-y*Math.sin(a), ry = x*Math.sin(a)+y*Math.cos(a);
    return [this.width/2 + rx*this.scale, this.height*.56 + (ry*Math.sin(p)-(z-10)*Math.cos(p))*this.scale, ry*Math.cos(p)+z*Math.sin(p)];
  }
  draw() {
    const c = this.canvas, ctx = this.ctx, r = c.getBoundingClientRect(); if (!r.width || !r.height) return;
    const ratio = Math.min(devicePixelRatio || 1, 2); c.width = r.width * ratio; c.height = r.height * ratio; ctx.setTransform(ratio,0,0,ratio,0,0);
    this.width = r.width; this.height = r.height;
    const boxes = assemblyGeometry(this.params,this.features,this.context,this.exploded);
    const extent=Math.max(110,...boxes.flatMap(b=>b.points.map(p=>2*Math.max(Math.abs(p[0]),Math.abs(p[1])))));
    this.scale = Math.min(r.width/extent,r.height/(extent*.76))*this.zoom;
    const css = getComputedStyle(document.documentElement), light = css.colorScheme === "light";
    ctx.fillStyle = css.getPropertyValue("--bg-2"); ctx.fillRect(0,0,r.width,r.height);
    ctx.lineWidth = .6; ctx.strokeStyle = light ? "#b7c5d24a" : "#6c86962b";
    for (let i=-100;i<=100;i+=10) { for (const line of [[[i,-100,-6],[i,100,-6]],[[-100,i,-6],[100,i,-6]]]) {
      const a=this.project(line[0]),b=this.project(line[1]);ctx.beginPath();ctx.moveTo(a[0],a[1]);ctx.lineTo(b[0],b[1]);ctx.stroke();
    } }
    const polys = [];
    for (const b of boxes) { if(this.hidden.has(b.id))continue;const ps=b.points.map(p=>this.project(p));
      for (const [i, face] of (b.faces||faces).entries()) { const vertices=face.map(k=>ps[k]);
        const cross=(vertices[1][0]-vertices[0][0])*(vertices[2][1]-vertices[0][1])-(vertices[1][1]-vertices[0][1])*(vertices[2][0]-vertices[0][0]);
        if(cross<=0)continue;
        polys.push({ id:b.id,vertices, depth: vertices.reduce((s,p)=>s+p[2],0)/vertices.length, color:b.color, shade:b.faces?i===1?1.08:.78+Math.sin(i)*.1:[.62,1.08,.76,.9,.94,.65][i], active:b.id===this.selected });
      }
    }
    polys.sort((a,b)=>a.depth-b.depth);
    this.polys=polys;
    for (const p of polys) { const rgb=p.color.map(n=>Math.round(Math.min(255,n*p.shade)));ctx.fillStyle=`rgb(${rgb.join(",")})`;
      ctx.strokeStyle=p.active?"rgba(255,179,103,.7)":"rgba(45,61,71,.4)";ctx.lineWidth=p.active?.8:.4;ctx.beginPath();p.vertices.forEach((v,i)=>i?ctx.lineTo(v[0],v[1]):ctx.moveTo(v[0],v[1]));ctx.closePath();ctx.fill();ctx.stroke(); }
    // Visual annotations do not pretend to be a CFD solution.
    if (["Mechanical","Thermal","Radiation"].includes(this.domain)) {
      ctx.strokeStyle=light?"#0b8aa3aa":"#5fd4e899";ctx.lineWidth=1.3;
      for(const x of [-12,0,12]){const a=this.project([x,-48,12]),b=this.project([x,-29,12]);ctx.beginPath();ctx.moveTo(a[0],a[1]);ctx.lineTo(b[0],b[1]);ctx.stroke();
        const angle=Math.atan2(b[1]-a[1],b[0]-a[0]);ctx.beginPath();ctx.moveTo(b[0]-6*Math.cos(angle-.5),b[1]-6*Math.sin(angle-.5));ctx.lineTo(b[0],b[1]);ctx.lineTo(b[0]-6*Math.cos(angle+.5),b[1]-6*Math.sin(angle+.5));ctx.stroke();}
    }
    const label=this.project([32,18,12]);ctx.fillStyle=css.getPropertyValue("--muted");ctx.font="9px Inter, sans-serif";ctx.fillText(this.context?.board?"PCB-linked cooling study":"Electrothermal study",Math.min(this.width-105,Math.max(12,label[0])),Math.max(50,label[1]));
  }
}
