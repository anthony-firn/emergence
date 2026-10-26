import {FOOTPRINTS, clone, snap, pcbExample, siliconExample, pads, bodyBounds, connectPads, moveComponent, pcbChecks, siliconChecks, documentFor, validateDocument, pcbBOM} from "./studio-model.js";

const node=(tag,text,cls)=>{const e=document.createElement(tag);if(text!==undefined)e.textContent=text;if(cls)e.className=cls;return e;};
const svgNode=(tag,attrs={},text)=>{const e=document.createElementNS("http://www.w3.org/2000/svg",tag);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,v);if(text!==undefined)e.textContent=text;return e;};
function button(text,id,fn,cls) {const b=node("button",text,cls);b.type="button";if(id)b.id=id;b.addEventListener("click",fn);return b;}
function download(data,name,type) {
  const url=URL.createObjectURL(new Blob([data],{type})),a=node("a");a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
}

class Studio {
  constructor(kind) {
    this.kind=kind;this.pcb=kind==="pcb";this.root=document.getElementById(kind+"-studio");
    this.design=this.pcb?pcbExample():siliconExample();this.undo=[];this.redo=[];this.selected=null;
    this.layers={top:true,bottom:true};this.airwires=true;this.routing=false;this.startPad=null;this.net=null;
    this.layer="top";this.routeWidth=.5;this.storage=true;this.readFailure=false;this.connectivity=false;
    try {const saved=localStorage.getItem("emergence-"+kind+"-v1");if(saved){const doc=validateDocument(JSON.parse(saved));if(doc.kind!==kind)throw new Error("Wrong document kind");this.design=doc.design;}}
    catch {this.storage=false;this.readFailure=true;}
    this.render();
    if(this.readFailure)this.message("The saved design could not be read. Its stored copy is preserved; export this session to keep edits.");
  }
  doc() {return documentFor(this.kind,this.design);}
  items() {return this.pcb?this.design.components:this.design.blocks;}
  message(text) {this.status.textContent=text;}
  transact(fn, message="Design updated.") {
    const before=clone(this.design);
    try {fn();this.doc();}
    catch(e) {this.design=before;this.render();this.message(e.message);return;}
    this.undo.push(before);if(this.undo.length>50)this.undo.shift();this.redo=[];
    this.save();this.render();this.message(message);
  }
  save() {
    if(this.readFailure)return;
    try {localStorage.setItem("emergence-"+this.kind+"-v1",JSON.stringify(this.doc()));this.storage=true;}
    catch {this.storage=false;}
  }
  history(back) {
    const from=back?this.undo:this.redo,to=back?this.redo:this.undo;if(!from.length)return;
    to.push(clone(this.design));this.design=from.pop();this.startPad=null;this.save();this.render();
  }
  render() {
    this.root.replaceChildren();
    const head=node("div",undefined,"studio-heading"),title=node("div");
    title.append(node("p",this.pcb?"PCB DESIGN / POWER STAGE":"SILICON DESIGN / MIXED-BLOCK EXAMPLE","panel-label"),node("h2",this.pcb?"From connections to copper.":"Give the architecture room."));
    title.append(node("p",this.pcb?"Place footprints, assign nets and connect pads across two copper layers.":"Place functional blocks, set core constraints and inspect placement conflicts."));
    const files=node("div",undefined,"studio-file-tools");
    files.append(button("↓ Design JSON",this.kind+"-export",()=>download(JSON.stringify(this.doc(),null,2)+"\n","emergence-"+this.kind+".json","application/json")));
    const upload=node("label","↑ Import design","file-button"),input=node("input");input.type="file";input.accept=".json,application/json";input.id=this.kind+"-import";upload.htmlFor=input.id;upload.append(input);
    input.addEventListener("change",async()=>{
      try {
        const file=input.files[0];if(!file)return;if(file.size>500000)throw new Error("Design files must be below 500 KB.");
        const doc=validateDocument(JSON.parse(await file.text()));if(doc.kind!==this.kind)throw new Error("Open this file in the "+doc.kind+" workspace.");
        this.transact(()=>{this.design=doc.design;this.selected=null;this.startPad=null;},"Imported design. Undo restores your previous document.");
      }catch(e){this.message("Import rejected: "+e.message);}finally{input.value="";}
    });
    files.append(upload,button(this.pcb?"↓ BOM CSV":"↓ Floorplan SVG",this.kind+"-output",()=>{
      if(this.pcb)download(pcbBOM(this.design),"emergence-pcb-bom.csv","text/csv");
      else download(this.exportSVG(),"emergence-floorplan-um.svg","image/svg+xml");
    }));
    head.append(title,files);this.root.append(head);
    const commands=node("div",undefined,"command-bar"),addLabel=node("label",this.pcb?"Footprint":"Block type"),types=node("select");types.id=this.kind+"-type";
    for(const [value,label] of this.pcb?Object.entries(FOOTPRINTS).map(([k,v])=>[k,v.title]):[["logic","Logic block"],["memory","Memory macro"],["io","I/O block"],["analog","Analog block"]]){const o=node("option",label);o.value=value;types.append(o);}
    addLabel.htmlFor=types.id;addLabel.append(types);
    commands.append(addLabel,button(this.pcb?"＋ Place component":"＋ Add block",this.kind+"-add",()=>this.add(types.value),"accent"));
    if(this.pcb) {
      const route=button(this.routing?"Cancel routing":"⌁ Connect pads",this.kind+"-route",()=>{this.routing=!this.routing;if(this.routing)this.connectivity=false;this.startPad=null;this.render();});
      route.setAttribute("aria-pressed",String(this.routing));commands.append(route);
      const netlist=button(this.connectivity?"▧ Board layout":"☷ Connectivity",this.kind+"-connectivity",()=>{this.connectivity=!this.connectivity;this.render();});
      netlist.setAttribute("aria-pressed",String(this.connectivity));commands.append(netlist);
    }
    const undo=button("↶ Undo",this.kind+"-undo",()=>this.history(true)),redo=button("↷ Redo",this.kind+"-redo",()=>this.history(false));
    undo.disabled=!this.undo.length;redo.disabled=!this.redo.length;commands.append(undo,redo);
    commands.append(node("span",this.storage?"Saved on this device":"Storage unavailable · export to keep edits","studio-save"));
    this.root.append(commands);
    if(!this.pcb) {
      const flow=node("div",undefined,"silicon-flow");
      for(const [label,state] of [["Floorplan","Editing here"],["RTL synthesis","Native Rust flow"],["Place & route","Planned"],["Timing sign-off","Planned"]]){const step=node("div");step.append(node("strong",label),node("small",state));flow.append(step);}
      this.root.append(flow);
    }
    const grid=node("div",undefined,"studio-grid"),explorer=node("aside",undefined,"studio-explorer"),centre=node("div",undefined,"studio-centre"),inspector=node("aside",undefined,"studio-inspector");
    explorer.append(node("p",this.pcb?"BOARD / mm":"FLOORPLAN / µm","panel-label"));
    this.number(explorer,"Width","width",this.design.width,n=>this.design.width=n,10,this.pcb?300:10000,this.pcb ? .5 : 10);
    this.number(explorer,"Height","height",this.design.height,n=>this.design.height=n,10,this.pcb?300:10000,this.pcb ? .5 : 10);
    this.number(explorer,"Snap grid","grid",this.design.grid,n=>this.design.grid=n,.1,100,.1);
    explorer.append(node("p",this.pcb?"COMPONENTS":"DESIGN HIERARCHY","panel-label studio-section-label"));
    this.items().forEach(item=>{
      const b=button((this.pcb?item.id+" · ":item.fixed?"▣ ":"◇ ")+item.name,null,()=>{this.selected=item.id;this.render();},"studio-tree-item");
      b.setAttribute("aria-pressed",String(this.selected===item.id));explorer.append(b);
    });
    if(this.pcb) {
      const checks=pcbChecks(this.design);explorer.append(node("p","NETS / REMAINING LINKS","panel-label studio-section-label"));
      checks.nets.forEach(net=>{
        const b=button(net.name+" · "+net.remaining,null,()=>{this.net=this.net===net.name?null:net.name;this.render();},"studio-tree-item");
        b.setAttribute("aria-pressed",String(this.net===net.name));explorer.append(b);
      });
    } else {
      const link=node("a","Native Verilog → gate netlist ↗");link.href="https://github.com/anthony-firn/emergence/blob/main/analysis/chip.md";link.className="studio-native-link";explorer.append(link);
    }
    const canvasTools=node("div",undefined,"studio-canvas-tools");
    if(this.pcb) {
      for(const [key,label] of [["top","Top copper"],["bottom","Bottom copper"],["airwires","Airwires"]]) {
        const l=node("label",undefined,"check-label"),c=node("input");c.type="checkbox";c.checked=key==="airwires"?this.airwires:this.layers[key];
        c.addEventListener("change",()=>{if(key==="airwires")this.airwires=c.checked;else this.layers[key]=c.checked;this.render();});l.append(c,document.createTextNode(label));canvasTools.append(l);
      }
    } else canvasTools.append(node("span","Block placement · schematic geometry"));
    canvasTools.append(node("span",this.design.grid+(this.pcb?" mm grid":" µm grid"),"mono"));centre.append(canvasTools);
    if(this.connectivity&&this.pcb)centre.append(this.netlist());else centre.append(this.canvas());
    centre.append(node("p",this.pcb?"Drag a component to snap its position, or select it and enter coordinates. Moving it removes its attached routes.":"Drag an unlocked block to snap its position, or select it and enter coordinates. Dark inner border marks the core keepout.","studio-hint"));
    this.status=node("p",this.routing?(this.startPad?"Choose a second pad on the same net.":"Routing: choose a start pad, then an end pad on the same net."):"Select an object to edit its properties.","studio-status");
    this.status.id=this.kind+"-status";this.status.setAttribute("role","status");centre.append(this.status);
    inspector.append(node("p","PROPERTIES","panel-label"));this.properties(inspector);
    grid.append(explorer,centre,inspector);this.root.append(grid);
    this.problems();
    this.root.append(node("p",this.pcb?
      "Layout study: footprints are illustrative through-hole geometry. Checks cover circular-pad and straight-track clearance, width, board edges and explicit pad connections. No obstacle avoidance, planes, vias, signal-integrity analysis or Gerber/drill output yet.":
      "Floorplanning study: block sizes are user inputs, with no PDK or process mapping. Checks cover core bounds, overlap, rectangular spacing halos and requested block area. Timing, power, transistor layout, DRC/LVS and GDSII export are not run here.","studio-scope"));
    this.root.append(node("p","This "+(this.pcb?"board":"floorplan")+" has its own local document and JSON file. Mechanical assembly branches and peer synchronization do not include it yet.","studio-scope"));
  }
  number(root,label,key,value,set,min,max,step) {
    const l=node("label",label),i=node("input");i.id=this.kind+"-"+key;i.type="number";i.min=min;i.max=max;i.step=step;i.value=value;l.htmlFor=i.id;l.append(i);root.append(l);
    i.addEventListener("change",()=>{
      if(i.value===""||!i.checkValidity()){i.value=value;this.message("Enter "+min+"–"+max+" for "+label+".");return;}
      this.transact(()=>set(Number(i.value)));
    });
  }
  field(root,label,key,value,set,max=60) {
    const l=node("label",label),i=node("input");i.id=this.kind+"-"+key;i.value=value;i.maxLength=max;l.htmlFor=i.id;l.append(i);root.append(l);
    i.addEventListener("change",()=>this.transact(()=>set(i.value)));
  }
  properties(root) {
    const item=this.items().find(i=>i.id===this.selected);
    if(item) {
      root.append(node("h3",this.pcb?item.id:item.name));
      this.field(root,"Name","name",item.name,v=>item.name=v);
      const shift=(key,n)=>{
        if(!this.pcb&&item.fixed)throw new Error("Unlock this block before changing its placement or size.");
        if(this.pcb)moveComponent(this.design,item,{[key]:n});else item[key]=n;
      };
      for(const [key,label] of [["x","X position"],["y","Y position"]])this.number(root,label,key,item[key],n=>shift(key,n),this.pcb?-300:-10000,this.pcb?600:20000,this.design.grid);
      if(this.pcb) {
        const l=node("label","Rotation"),select=node("select");select.id="pcb-rotation";l.htmlFor=select.id;
        for(const angle of [0,90,180,270]){const o=node("option",angle+"°");o.value=angle;select.append(o);}select.value=item.rotation;
        select.addEventListener("change",()=>this.transact(()=>moveComponent(this.design,item,{rotation:Number(select.value)})));l.append(select);root.append(l);
        root.append(node("p","PAD NET ASSIGNMENTS","panel-label studio-section-label"));
        item.nets.forEach((net,i)=>this.field(root,"Pad "+(i+1),"net-"+(i+1),net,v=>{moveComponent(this.design,item,{});item.nets[i]=v;},32));
      } else {
        for(const key of ["width","height"])this.number(root,key==="width"?"Block width":"Block height","block-"+key,item[key],n=>shift(key,n),1,10000,this.design.grid);
        const l=node("label",undefined,"check-label"),c=node("input");c.id="silicon-fixed";c.type="checkbox";c.checked=item.fixed;
        c.addEventListener("change",()=>this.transact(()=>item.fixed=c.checked));l.append(c,document.createTextNode("Lock placement and size"));root.append(l);
      }
      root.append(button("Delete "+(this.pcb?"component":"block"),this.kind+"-delete",()=>this.transact(()=>{
        if(!this.pcb&&item.fixed)throw new Error("Unlock this block before deleting it.");
        if(this.pcb){moveComponent(this.design,item,{});this.design.components=this.design.components.filter(c=>c.id!==item.id);}else this.design.blocks=this.design.blocks.filter(b=>b.id!==item.id);
        this.selected=null;
      })));
    } else root.append(node("p","Select an object in the explorer or drawing.","small muted"));
    root.append(node("p",this.pcb?"DESIGN RULES / mm":"CORE CONSTRAINTS / µm","panel-label studio-section-label"));
    if(this.pcb) {
      for(const [key,label,min,max] of [["clearance","Copper clearance",.05,5],["minWidth","Minimum width",.1,5],["edge","Copper to edge",0,10]])this.number(root,label,key,this.design.rules[key],n=>this.design.rules[key]=n,min,max,.05);
      root.append(node("p","ROUTING","panel-label studio-section-label"));
      const layerLabel=node("label","Active layer"),layer=node("select");layer.id="pcb-layer";layerLabel.htmlFor=layer.id;
      for(const [key,label] of [["top","Top copper"],["bottom","Bottom copper"]]){const o=node("option",label);o.value=key;layer.append(o);}layer.value=this.layer;
      layer.addEventListener("change",()=>this.layer=layer.value);layerLabel.append(layer);root.append(layerLabel);
      const widthLabel=node("label","Connection width"),width=node("input");width.id="pcb-route-width";width.type="number";width.min=.1;width.max=5;width.step=.1;width.value=this.routeWidth;widthLabel.htmlFor=width.id;widthLabel.append(width);root.append(widthLabel);
      width.addEventListener("change",()=>{if(width.value!==""&&width.checkValidity())this.routeWidth=Number(width.value);else{width.value=this.routeWidth;this.message("Connection width must be 0.1–5 mm.");}});
      root.append(node("p","Click two same-net pads with Connect pads active. Connections use a 45° bend; geometry checks report collisions after placement.","small muted"));
      if(this.design.routes.length) {
        const l=node("label","Connection to remove"),s=node("select");s.id="pcb-route-selection";l.htmlFor=s.id;
        for(const route of this.design.routes){const o=node("option",route.from+" → "+route.to+" / "+route.layer);o.value=route.id;s.append(o);}l.append(s);root.append(l);
        root.append(button("Remove connection","pcb-remove-route",()=>this.transact(()=>{this.design.routes=this.design.routes.filter(r=>r.id!==s.value);})));
      }
    } else {
      for(const [key,label,min,max,step] of [["margin","Edge keepout",0,1000,1],["spacing","Placement halo",0,1000,1],["target","Area target (%)",1,100,1]])this.number(root,label,key,this.design[key],n=>this.design[key]=n,min,max,step);
    }
  }
  add(type) {
    this.transact(()=>{
      if(this.items().length>=64)throw new Error("This editor supports 64 objects.");
      const key=(this.pcb?"C":"block")+"-"+crypto.randomUUID().slice(0,8);
      const x=snap(this.design.width/2,this.design.grid),y=snap(this.design.height/2,this.design.grid);
      if(this.pcb)this.design.components.push({id:key,name:FOOTPRINTS[type].title,footprint:type,x,y,rotation:0,nets:FOOTPRINTS[type].pads.map(()=>"")});
      else this.design.blocks.push({id:key,name:type+" block",kind:type,x,y,width:100,height:100,fixed:false});
      this.selected=key;
    },"Object placed at the centre. Set its position and inspect the checks.");
  }
  padClick(p) {
    if(!this.routing){this.selected=p.component;this.render();return;}
    if(!this.startPad){this.startPad=p.id;this.render();return;}
    try {
      const route=connectPads(this.design,this.startPad,p.id,this.layer,this.routeWidth);
      this.transact(()=>{this.design.routes.push(route);this.startPad=null;},"Connection added. Inspect clearance checks below the board.");
    }catch(e){this.message(e.message);}
  }
  canvas() {
    const d=this.design,pad=Math.max(d.width,d.height)*.055;
    const svg=svgNode("svg",{viewBox:[-pad,-pad,d.width+2*pad,d.height+2*pad].join(" "),"aria-label":this.pcb?"PCB layout in millimetres":"Silicon floorplan in micrometres",class:"studio-canvas",id:this.kind+"-canvas"});
    const defs=svgNode("defs"),pattern=svgNode("pattern",{id:this.kind+"-grid",width:d.grid,height:d.grid,patternUnits:"userSpaceOnUse"});
    pattern.append(svgNode("circle",{cx:0,cy:0,r:this.pcb ? .08 : .9,fill:"#8a9ead",opacity:.35}));defs.append(pattern);svg.append(defs);
    svg.append(svgNode("rect",{x:0,y:0,width:d.width,height:d.height,class:"studio-board"}),svgNode("rect",{x:0,y:0,width:d.width,height:d.height,fill:"url(#"+this.kind+"-grid)", "pointer-events":"none"}));
    if(this.pcb) {
      const all=pads(d),byID=new Map(all.map(p=>[p.id,p])),checks=pcbChecks(d);
      if(this.airwires)for(const net of checks.nets)for(const [a,b] of net.airwires)svg.append(svgNode("line",{x1:a.x,y1:a.y,x2:b.x,y2:b.y,class:"airwire",opacity:!this.net||this.net===net.name ? .7 : .12}));
      for(const r of d.routes)if(this.layers[r.layer])svg.append(svgNode("polyline",{points:r.points.map(p=>p.join(",")).join(" "),fill:"none",stroke:r.layer==="top"?"#f19177":"#71afff","stroke-width":r.width,"stroke-linecap":"round","stroke-linejoin":"round",opacity:!this.net||byID.get(r.from).net===this.net?1:.16}));
      for(const c of d.components) {
        const b=bodyBounds(c),g=this.object(svg,c,this.pcb?c.id:c.name);
        g.append(svgNode("rect",{...b,rx:.5,class:"pcb-body"}),svgNode("text",{x:c.x,y:c.y-.7,"font-size":2,"text-anchor":"middle",class:"component-label"},c.id));
      }
      for(const p of all) {
        const g=svgNode("g",{role:"button",tabindex:0,"aria-label":p.id+" "+(p.net||"unassigned"),class:"pcb-pad", "data-pad":p.id});
        const selected=this.startPad===p.id;
        g.append(svgNode("circle",{cx:p.x,cy:p.y,r:1.4,fill:"transparent"}),svgNode("circle",{cx:p.x,cy:p.y,r:p.radius,fill:selected?"#fff":"#dfba6b",stroke:selected?"#5fd4e8":"#2a3542","stroke-width":selected ? .45 : .15}),svgNode("circle",{cx:p.x,cy:p.y,r:.3,fill:"#13232d"}));
        const title=svgNode("title",{},p.id+" · "+(p.net||"no net"));g.append(title);
        g.addEventListener("click",()=>this.padClick(p));
        g.addEventListener("keydown",e=>{if(e.key==="Enter"||e.key===" "){e.preventDefault();this.padClick(p);}});
        svg.append(g);
      }
    } else {
      svg.append(svgNode("rect",{x:d.margin,y:d.margin,width:d.width-2*d.margin,height:d.height-2*d.margin,class:"core-boundary"}));
      const colors={logic:"#65b3bd",memory:"#a493db",io:"#d4ae69",analog:"#d78d91"};
      for(const b of d.blocks) {
        const g=this.object(svg,b,b.name);
        g.append(svgNode("rect",{x:b.x,y:b.y,width:b.width,height:b.height,fill:colors[b.kind],"fill-opacity":.22,stroke:colors[b.kind],"stroke-width":2}));
        g.append(svgNode("text",{x:b.x+12,y:b.y+28,"font-size":18,class:"component-label"},(b.fixed?"▣ ":"")+b.name),svgNode("text",{x:b.x+12,y:b.y+52,"font-size":12,class:"component-label"},b.width+" × "+b.height+" µm"));
      }
    }
    let drag;
    const coordinates=e=>{const p=svg.createSVGPoint();p.x=e.clientX;p.y=e.clientY;return p.matrixTransform(svg.getScreenCTM().inverse());};
    svg.addEventListener("pointerdown",e=>{
      const group=e.target.closest("[data-object]");if(!group||this.routing||e.button!==0)return;
      const item=this.items().find(i=>i.id===group.dataset.object);if(!item)return;
      const p=coordinates(e);drag={id:item.id,x:p.x,y:p.y,startX:item.x,startY:item.y,locked:item.fixed};
      svg.setPointerCapture(e.pointerId);
    });
    svg.addEventListener("pointerup",e=>{
      if(!drag)return;const p=coordinates(e),item=this.items().find(i=>i.id===drag.id),distance=Math.hypot(p.x-drag.x,p.y-drag.y),was=drag;drag=null;
      if(distance<d.grid/2){this.selected=was.id;this.render();return;}
      if(was.locked){this.message("Unlock this block before moving it.");return;}
      const x=snap(was.startX+p.x-was.x,d.grid),y=snap(was.startY+p.y-was.y,d.grid);
      this.transact(()=>{if(this.pcb)moveComponent(d,item,{x,y});else Object.assign(item,{x,y});this.selected=item.id;},"Placement updated on the snap grid.");
    });
    svg.addEventListener("pointercancel",()=>drag=null);
    return svg;
  }
  object(svg,item,label) {
    const g=svgNode("g",{"data-object":item.id,role:"button",tabindex:0,"aria-label":label,"aria-pressed":String(this.selected===item.id),class:"studio-object"});
    // Pointer selection happens on pointerup so one gesture is one undo step.
    g.addEventListener("keydown",e=>{if(e.key==="Enter"||e.key===" "){e.preventDefault();this.selected=item.id;this.render();}});
    svg.append(g);return g;
  }
  netlist() {
    const wrap=node("div",undefined,"studio-netlist");wrap.append(node("h3","Connectivity"));
    wrap.append(node("p","Logical net membership from pad assignments. This is a connection list, not a circuit simulation."));
    for(const net of pcbChecks(this.design).nets) {
      const row=node("section");row.append(node("strong",net.name),node("span",net.remaining+" explicit links remaining"));
      const members=node("div");for(const p of pads(this.design).filter(p=>p.net===net.name))members.append(button(p.id,null,()=>{this.selected=p.component;this.render();}));
      row.append(members);wrap.append(row);
    }return wrap;
  }
  problems() {
    const results=this.pcb?pcbChecks(this.design):siliconChecks(this.design),section=node("section",undefined,"studio-checks");
    section.id=this.kind+"-checks";section.setAttribute("aria-label",this.pcb?"PCB geometry checks":"Floorplan placement checks");
    const top=node("div",undefined,"section-top");top.append(node("h3",this.pcb?"Board geometry checks":"Placement checks"));
    top.append(node("span",results.issues.length+" issues"+(this.pcb?" · "+results.unrouted+" explicit links remaining":" · "+results.utilization.toFixed(1)+"% requested core area")));section.append(top);
    if(!results.issues.length)section.append(node("p","No violations in the implemented geometry checks.","studio-check-pass"));
    else {
      const list=node("ul");
      for(const issue of results.issues.slice(0,50)){const li=node("li");li.dataset.code=issue.code;li.append(button(issue.message,null,()=>{
        if(this.items().some(i=>i.id===issue.object)){this.selected=issue.object;this.render();}
        else this.message(issue.message);
      }));list.append(li);}
      section.append(list);if(results.issues.length>50)section.append(node("p","Showing the first 50 issues."));
    }
    if(this.pcb)section.append(node("p","Connectivity is counted from explicit pad-to-pad links; same-net copper intersections are not inferred.","small muted"));
    this.root.append(section);
  }
  exportSVG() {
    const svg=this.root.querySelector("svg").cloneNode(true),style=svgNode("style",{},".studio-board{fill:#12212b;stroke:#6c8998;stroke-width:2}.core-boundary{fill:none;stroke:#8998a5;stroke-width:2;stroke-dasharray:8 5}.component-label{fill:#e5edf4;font-family:sans-serif}");
    svg.prepend(style);svg.setAttribute("xmlns","http://www.w3.org/2000/svg");
    svg.append(svgNode("text",{x:0,y:-12,fill:"#e5edf4","font-size":12},"Concept floorplan · µm · no PDK / no sign-off"));
    return new XMLSerializer().serializeToString(svg);
  }
}

new Studio("pcb");new Studio("silicon");
const engines=[...document.querySelectorAll("[data-studio]")];
function selectStudio(name,focus=false) {
  if(!engines.some(b=>b.dataset.studio===name))return;
  for(const b of engines){const active=b.dataset.studio===name;b.setAttribute("aria-selected",String(active));b.tabIndex=active?0:-1;document.getElementById(b.dataset.studio+"-studio").hidden=!active;if(active&&focus)b.focus();}
}
engines.forEach((b,i)=>{
  b.addEventListener("click",()=>{selectStudio(b.dataset.studio);history.replaceState(null,"","#"+b.dataset.studio);});
  b.addEventListener("keydown",e=>{
    const n=e.key==="ArrowRight"?(i+1)%engines.length:e.key==="ArrowLeft"?(i+engines.length-1)%engines.length:e.key==="Home"?0:e.key==="End"?engines.length-1:-1;
    if(n>=0){e.preventDefault();selectStudio(engines[n].dataset.studio,true);history.replaceState(null,"","#"+engines[n].dataset.studio);}
  });
});
selectStudio(location.hash.slice(1)||"mechanical");
addEventListener("hashchange",()=>selectStudio(location.hash.slice(1)));
