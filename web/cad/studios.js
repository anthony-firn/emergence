import {FOOTPRINTS, pads, bodyBounds, pcbChecks, siliconChecks, documentFor, validateDocument, pcbBOM, RTL_SOURCES} from "./studio-model.js";
import {project, onSolved} from "./app.js";
import {request} from "./models.js";
import {kernelBoard, densityIssues, kicadPCB, gerberCheck, netOf, OZ} from "./domains/pcb.js";
import {STAGES, gdsii, verifyGDS, powerSummary} from "./domains/silicon.js";
import {zip} from "./core/integrity.js";
import {format} from "./core/units.js";
import {slug} from "./core/project.js";
const named=()=>slug(project.store.docs.project.name);
import {mountLitho} from "./ui/litho-view.js";
import {checkCards, certificates, assumptionList} from "./ui/kit.js";

const node=(tag,text,cls)=>{const e=document.createElement(tag);if(text!==undefined)e.textContent=text;if(cls)e.className=cls;return e;};
const svgNode=(tag,attrs={},text)=>{const e=document.createElementNS("http://www.w3.org/2000/svg",tag);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,v);if(text!==undefined)e.textContent=text;return e;};
function button(text,id,fn,cls) {const b=node("button",text,cls);b.type="button";if(id)b.id=id;b.addEventListener("click",fn);return b;}
function download(data,name,type) {
  const url=URL.createObjectURL(new Blob([data],{type})),a=node("a");a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
}

// The board and floorplan editors are clients of the shared command registry:
// every edit below is a pcb.*, silicon.* or system.* command with one undo.
class Studio {
  constructor(kind) {
    this.kind=kind;this.pcb=kind==="pcb";this.root=document.getElementById(kind+"-studio");this.selected=null;
    this.layers={top:true,bottom:true};this.airwires=true;this.routing=false;this.avoid=true;this.startPad=null;this.net=null;
    this.layer="top";this.routeWidth=.5;this.connectivity=false;this.heat=false;this.solution=null;this.lastMessage=null;
    project.store.on(({ids,reason})=>{if(ids.includes(kind)){if(reason!=="edit")this.startPad=null;this.render();}});
    // The undo entry is pushed after the store change that triggered render().
    project.bus.on(e=>{if(e.type==="reset")this.render();else this.syncUndo();});
    onSolved(s=>{this.solution=s;this.renderResults();});
    this.render();
  }
  get design() {return project.store.get(this.kind);}
  syncUndo() {if(!this.undoButtons)return;const [u,r]=this.undoButtons;u.disabled=!project.bus.undoStack.length;r.disabled=!project.bus.redoStack.length;u.title=project.bus.undoStack.at(-1)?"Undo: "+project.bus.undoStack.at(-1).label:"Nothing to undo";}
  doc() {return documentFor(this.kind,this.design);}
  items() {return this.pcb?this.design.components:this.design.blocks;}
  message(text) {this.lastMessage=text;if(this.status)this.status.textContent=text;}
  run(op,args,message) {
    return project.bus.run(op,args).then(r=>{this.message(r?.message||message||"Design updated.");return r;}).catch(e=>{this.render();this.message(e.message);});
  }
  render() {
    const scrollY=window.scrollY;
    this.root.replaceChildren();
    const head=node("div",undefined,"studio-heading"),title=node("div");
    title.append(node("p",`${this.pcb?"PCB DESIGN":"SILICON DESIGN"} / ${project.store.docs.project.name.toUpperCase()}`,"panel-label"),node("h2",this.pcb?"From connections to copper to fabrication.":"From RTL to timing to floorplan."));
    title.append(node("p",this.pcb?"Place footprints, route around obstacles on two layers, solve trace resistance and impedance in the kernel, then export Gerber, drill and KiCad files.":"Synthesize Verilog onto characterised cells, time it at the die's own temperature, place the blocks and image a layout clip."));
    const files=node("div",undefined,"studio-file-tools");
    files.append(button("↓ Design JSON",this.kind+"-export",()=>download(JSON.stringify(this.doc(),null,2)+"\n","emergence-"+this.kind+".json","application/json")));
    const upload=node("label","↑ Import design","file-button"),input=node("input");input.type="file";input.accept=".json,application/json";input.id=this.kind+"-import";upload.htmlFor=input.id;upload.append(input);
    input.addEventListener("change",async()=>{
      try {
        const file=input.files[0];if(!file)return;if(file.size>500000)throw new Error("Design files must be below 500 KB.");
        const doc=validateDocument(JSON.parse(await file.text()));if(doc.kind!==this.kind)throw new Error("Open this file in the "+doc.kind+" workspace.");
        this.selected=null;this.startPad=null;await project.bus.run(this.kind+".load",{design:doc});this.message("Imported design. Undo restores your previous document.");
      }catch(e){this.message("Import rejected: "+e.message);}finally{input.value="";}
    });
    files.append(upload,button(this.pcb?"↓ BOM CSV":"↓ Floorplan SVG",this.kind+"-output",()=>{
      if(this.pcb)download(pcbBOM(this.design),"emergence-pcb-bom.csv","text/csv");
      else download(this.exportSVG(),"emergence-floorplan-um.svg","image/svg+xml");
    }));
    if(this.pcb)files.append(button("↓ Fabrication package","pcb-fab",()=>this.fabrication(),"accent"),button("↓ KiCad board","pcb-kicad",()=>{download(kicadPCB(this.design),`${named()}.kicad_pcb`,"application/x-kicad-pcb");this.message("KiCad 8 board written to the documented S-expression format. Open it in KiCad to confirm; kicad-cli has not checked this file.");}));
    else files.append(button("↓ GDSII","silicon-gds",()=>this.exportGDS(),"accent"));
    head.append(title,files);this.root.append(head);
    const commands=node("div",undefined,"command-bar"),addLabel=node("label",this.pcb?"Footprint":"Block type"),types=node("select");types.id=this.kind+"-type";
    for(const [value,label] of this.pcb?Object.entries(FOOTPRINTS).map(([k,v])=>[k,v.title]):[["logic","Logic block"],["memory","Memory macro"],["io","I/O block"],["analog","Analog block"]]){const o=node("option",label);o.value=value;types.append(o);}
    addLabel.htmlFor=types.id;addLabel.append(types);
    commands.append(addLabel,button(this.pcb?"＋ Place component":"＋ Add block",this.kind+"-add",()=>this.add(types.value),"accent"));
    if(this.pcb) {
      const route=button(this.routing?"Cancel routing":"⌁ Connect pads",this.kind+"-route",()=>{this.routing=!this.routing;if(this.routing)this.connectivity=false;this.startPad=null;this.render();});
      route.setAttribute("aria-pressed",String(this.routing));commands.append(route);
      const avoid=node("label",undefined,"check-label"),a=node("input");a.type="checkbox";a.id="pcb-avoid";a.checked=this.avoid;a.addEventListener("change",()=>this.avoid=a.checked);avoid.append(a,document.createTextNode("Avoid obstacles"));commands.append(avoid);
      commands.append(button("⌁ Route all","pcb-route-all",()=>this.run("pcb.route_all",{width:this.routeWidth})));
      const netlist=button(this.connectivity?"▧ Board layout":"☷ Connectivity",this.kind+"-connectivity",()=>{this.connectivity=!this.connectivity;this.render();});
      netlist.setAttribute("aria-pressed",String(this.connectivity));commands.append(netlist);
    } else {
      const heat=button(this.heat?"▧ Block kinds":"▤ Power density",this.kind+"-heat",()=>{this.heat=!this.heat;this.render();});heat.setAttribute("aria-pressed",String(this.heat));commands.append(heat);
    }
    const undo=button("↶ Undo",this.kind+"-undo",()=>project.bus.undo()),redo=button("↷ Redo",this.kind+"-redo",()=>project.bus.redo());
    this.undoButtons=[undo,redo];this.syncUndo();commands.append(undo,redo);
    commands.append(node("span","Part of the project draft · saved with each revision","studio-save"));
    this.root.append(commands);
    if(!this.pcb){this.flow=node("div",undefined,"silicon-flow stage-flow");this.root.append(this.flow);}
    const grid=node("div",undefined,"studio-grid"),explorer=node("aside",undefined,"studio-explorer"),centre=node("div",undefined,"studio-centre"),inspector=node("aside",undefined,"studio-inspector");
    explorer.append(node("p",this.pcb?"BOARD / mm":"FLOORPLAN / µm","panel-label"));
    this.number(explorer,"Width","width",this.design.width,n=>this.pcb?["pcb.board",{width:n}]:["silicon.core",{width:n}],10,this.pcb?300:10000,this.pcb ? .5 : 10);
    this.number(explorer,"Height","height",this.design.height,n=>this.pcb?["pcb.board",{height:n}]:["silicon.core",{height:n}],10,this.pcb?300:10000,this.pcb ? .5 : 10);
    this.number(explorer,"Snap grid","grid",this.design.grid,n=>this.pcb?["pcb.board",{grid:n}]:["silicon.core",{grid:n}],.1,100,.1);
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
    } else canvasTools.append(node("span",this.heat?"Block power density · W/mm²":"Block placement · schematic geometry"));
    canvasTools.append(node("span",this.design.grid+(this.pcb?" mm grid":" µm grid"),"mono"));centre.append(canvasTools);
    if(this.connectivity&&this.pcb)centre.append(this.netlist());else centre.append(this.canvas());
    centre.append(node("p",this.pcb?"Drag a component to snap its position, or select it and enter coordinates. Moving it removes its attached routes.":"Drag an unlocked block to snap its position, or select it and enter coordinates. Dark inner border marks the core keepout.","studio-hint"));
    const hint=this.routing?(this.startPad?"Choose a second pad on the same net.":"Routing: choose a start pad, then an end pad on the same net."):"Select an object to edit its properties.";
    this.status=node("p",this.lastMessage||hint,"studio-status");
    this.status.id=this.kind+"-status";this.status.setAttribute("role","status");centre.append(this.status);this.lastMessage=null;
    inspector.append(node("p","PROPERTIES","panel-label"));this.properties(inspector);
    grid.append(explorer,centre,inspector);this.root.append(grid);
    this.problems();
    this.results=node("section",undefined,"domain-results");this.results.id=this.kind+"-results";this.root.append(this.results);this.renderResults();
    if(!this.pcb)mountLitho(this.root);
    const scope=node("details",undefined,"assumption-details");scope.append(node("summary","Assumptions, scope and what is not modelled"));
    scope.append(node("p",this.pcb?
      "Layout study: footprints are illustrative through-hole geometry. Checks cover circular-pad and straight-track clearance, width, board edges, current density and explicit pad connections. Trace resistance, Gerber and Excellon come from the Rust cad::kits::pcb kernel; the files follow their specifications but are not verified against a CAM tool.":
      "Flow study: RTL uses the native Verilog subset and the repository's 150 nm transistor-characterised cells. Block sizes and non-RTL powers are planning inputs with no PDK. Place and route, extraction, DRC/LVS and foundry GDS are not run here.","studio-scope"));
    assumptionList(scope,this.kind);this.root.append(scope);
    if(document.activeElement===document.body)window.scrollTo({top:scrollY,behavior:"instant"});
  }
  number(root,label,key,value,command,min,max,step) {
    const l=node("label",label),i=node("input");i.id=this.kind+"-"+key;i.type="number";i.min=min;i.max=max;i.step=step;i.value=value;l.htmlFor=i.id;l.append(i);root.append(l);
    i.addEventListener("change",()=>{
      if(i.value===""||!i.checkValidity()){i.value=value;this.message("Enter "+min+"–"+max+" for "+label+".");return;}
      const [op,args]=command(Number(i.value));this.run(op,args);
    });
    return i;
  }
  field(root,label,key,value,command,max=60) {
    const l=node("label",label),i=node("input");i.id=this.kind+"-"+key;i.value=value;i.maxLength=max;l.htmlFor=i.id;l.append(i);root.append(l);
    i.addEventListener("change",()=>{const [op,args]=command(i.value);this.run(op,args);});
  }
  properties(root) {
    const item=this.items().find(i=>i.id===this.selected);
    if(item) {
      root.append(node("h3",this.pcb?item.id:item.name));
      this.field(root,"Name","name",item.name,v=>this.pcb?["pcb.rename",{id:item.id,name:v}]:["silicon.update_block",{id:item.id,name:v}]);
      for(const [key,label] of [["x","X position"],["y","Y position"]])this.number(root,label,key,item[key],n=>this.pcb?["pcb.move",{id:item.id,at:key==="x"?[n,item.y]:[item.x,n]}]:["silicon.update_block",{id:item.id,[key]:n}],this.pcb?-300:-10000,this.pcb?600:20000,this.design.grid);
      if(this.pcb) {
        const l=node("label","Rotation"),select=node("select");select.id="pcb-rotation";l.htmlFor=select.id;
        for(const angle of [0,90,180,270]){const o=node("option",angle+"°");o.value=angle;select.append(o);}select.value=item.rotation;
        select.addEventListener("change",()=>this.run("pcb.move",{id:item.id,rotation:select.value}));l.append(select);root.append(l);
        this.number(root,"Dissipation (W)","power",this.design.powers[item.id]??0,n=>["pcb.power",{id:item.id,watts:n}],0,500,.01);
        const linked=item.id==="U1"&&this.solution?.applied["pcb.u1_power"];
        if(linked)root.append(node("p",`⛓ Driven by ${linked.from.join(", ")}: ${format(linked.value,"W",4)}`,"link-note"));
        root.append(node("p","PAD NET ASSIGNMENTS","panel-label studio-section-label"));
        item.nets.forEach((net,i)=>this.field(root,"Pad "+(i+1),"net-"+(i+1),net,v=>["pcb.net",{pad:item.id+":"+(i+1),net:v||"-"}],32));
      } else {
        for(const key of ["width","height"])this.number(root,key==="width"?"Block width":"Block height","block-"+key,item[key],n=>["silicon.update_block",{id:item.id,[key]:n}],1,10000,this.design.grid);
        const rtl=this.design.chip.rtl.find(r=>r.block===item.id);
        this.number(root,rtl?"Power (RTL replaces this) mW":"Power (mW)","block-power",item.power,n=>["silicon.update_block",{id:item.id,power:n}],0,1e6,1);
        const l=node("label",undefined,"check-label"),c=node("input");c.id="silicon-fixed";c.type="checkbox";c.checked=item.fixed;
        c.addEventListener("change",()=>this.run("silicon.update_block",{id:item.id,fixed:c.checked}));l.append(c,document.createTextNode("Lock placement and size"));root.append(l);
        const r=node("label","Implemented by RTL"),rs=node("select");rs.id="silicon-rtl-source";r.htmlFor=rs.id;
        for(const [k,t] of [["none","No: estimated power"],...Object.entries(RTL_SOURCES)]){const o=node("option",t);o.value=k;rs.append(o);}rs.value=rtl?.source??"none";
        rs.addEventListener("change",()=>{if(rs.value==="custom"){this.rtlBlock=item.id;this.rtlOpen=true;this.render();return;}this.run("silicon.rtl",{block:item.id,source:rs.value});});r.append(rs);root.append(r);
      }
      root.append(button("Delete "+(this.pcb?"component":"block"),this.kind+"-delete",()=>{const id=item.id;this.selected=null;this.run(this.pcb?"pcb.remove":"silicon.remove_block",{id}).then(r=>{if(!r){this.selected=id;this.render();}});}));
    } else root.append(node("p","Select an object in the explorer or drawing.","small muted"));
    root.append(node("p",this.pcb?"DESIGN RULES / mm":"CORE CONSTRAINTS / µm","panel-label studio-section-label"));
    if(this.pcb) {
      for(const [key,arg,label,min,max] of [["clearance","clearance","Copper clearance",.05,5],["minWidth","min_width","Minimum width",.1,5],["edge","edge","Copper to edge",0,10]])this.number(root,label,key,this.design.rules[key],n=>["pcb.board",{[arg]:n}],min,max,.05);
      this.number(root,"Current density limit (A/mm²)","density",this.design.rules.density,n=>["pcb.board",{density:n}],1,1000,1);
      root.append(node("p","STACKUP","panel-label studio-section-label"));
      this.number(root,"Board thickness (mm)","thickness",this.design.stackup.thickness,n=>["pcb.board",{thickness:n}],.2,6,.1);
      this.number(root,"Copper (oz/ft²)","copper",this.design.stackup.copper,n=>["pcb.board",{copper:n}],.25,6,.25);
      this.number(root,"Dielectric εr","epsr",this.design.stackup.epsr,n=>["pcb.board",{epsr:n}],1,20,.05);
      root.append(node("p","ROUTING","panel-label studio-section-label"));
      const layerLabel=node("label","Active layer"),layer=node("select");layer.id="pcb-layer";layerLabel.htmlFor=layer.id;
      for(const [key,label] of [["top","Top copper"],["bottom","Bottom copper"]]){const o=node("option",label);o.value=key;layer.append(o);}layer.value=this.layer;
      layer.addEventListener("change",()=>this.layer=layer.value);layerLabel.append(layer);root.append(layerLabel);
      const widthLabel=node("label","Connection width"),width=node("input");width.id="pcb-route-width";width.type="number";width.min=.1;width.max=5;width.step=.1;width.value=this.routeWidth;widthLabel.htmlFor=width.id;widthLabel.append(width);root.append(widthLabel);
      width.addEventListener("change",()=>{if(width.value!==""&&width.checkValidity())this.routeWidth=Number(width.value);else{width.value=this.routeWidth;this.message("Connection width must be 0.1–5 mm.");}});
      root.append(node("p","Click two same-net pads with Connect pads active. With Avoid obstacles on, the router finds a clear 45° path; otherwise a direct 45° bend is drawn and the checks report collisions.","small muted"));
      if(this.design.routes.length) {
        const l=node("label","Connection"),s=node("select");s.id="pcb-route-selection";l.htmlFor=s.id;
        for(const route of this.design.routes){const o=node("option",route.from+" → "+route.to+" / "+route.layer+" / "+route.width+" mm");o.value=route.id;s.append(o);}l.append(s);root.append(l);
        root.append(button("Remove connection","pcb-remove-route",()=>this.run("pcb.unroute",{id:s.value})),button("Solve impedance","pcb-impedance",()=>this.impedance(s.value)));
        this.impedanceOut=node("p","","impedance-out");this.impedanceOut.id="pcb-impedance-out";root.append(this.impedanceOut);
      }
    } else {
      for(const [key,label,min,max,step] of [["margin","Edge keepout",0,1000,1],["spacing","Placement halo",0,1000,1],["target","Area target (%)",1,100,1]])this.number(root,label,key,this.design[key],n=>["silicon.core",{[key]:n}],min,max,step);
    }
  }
  add(type) {
    const d=this.design,at=[d.width/2,d.height/2];
    (this.pcb?project.bus.run("pcb.place",{footprint:type,at}):project.bus.run("silicon.add_block",{kind:type,at,width:100,height:100}))
      .then(r=>{this.selected=r.id;this.render();this.message("Object placed at the centre. Set its position and inspect the checks.");}).catch(e=>this.message(e.message));
  }
  padClick(p) {
    if(!this.routing){this.selected=p.component;this.render();return;}
    if(!this.startPad){this.startPad=p.id;this.render();return;}
    const from=this.startPad;
    if(pads(this.design).find(x=>x.id===from)?.net!==p.net||!p.net||from===p.id){this.message("Choose two different pads on the same assigned net.");return;}
    this.startPad=null;
    this.run(this.avoid?"pcb.route":"pcb.connect",{from,to:p.id,layer:this.layer,width:this.routeWidth},"Connection added. Inspect clearance checks below the board.").then(r=>{if(!r){this.startPad=from;this.render();}});
  }
  canvas() {
    const d=this.design,pad=Math.max(d.width,d.height)*.055;
    const svg=svgNode("svg",{viewBox:[-pad,-pad,d.width+2*pad,d.height+2*pad].join(" "),"aria-label":this.pcb?"PCB layout in millimetres":"Silicon floorplan in micrometres",class:"studio-canvas",id:this.kind+"-canvas"});
    const defs=svgNode("defs"),pattern=svgNode("pattern",{id:this.kind+"-grid",width:d.grid,height:d.grid,patternUnits:"userSpaceOnUse"});
    pattern.append(svgNode("circle",{cx:0,cy:0,r:this.pcb ? .08 : .9,fill:"#8a9ead",opacity:.35}));defs.append(pattern);svg.append(defs);
    svg.append(svgNode("rect",{x:0,y:0,width:d.width,height:d.height,class:"studio-board"}),svgNode("rect",{x:0,y:0,width:d.width,height:d.height,fill:"url(#"+this.kind+"-grid)", "pointer-events":"none"}));
    if(this.pcb) {
      const all=pads(d),byID=new Map(all.map(p=>[p.id,p])),checks=pcbChecks(d),hot=new Set(densityIssues(d).map(i=>i.object));
      if(this.airwires)for(const net of checks.nets)for(const [a,b] of net.airwires)svg.append(svgNode("line",{x1:a.x,y1:a.y,x2:b.x,y2:b.y,class:"airwire",opacity:!this.net||this.net===net.name ? .7 : .12}));
      for(const r of d.routes)if(this.layers[r.layer]){const line=svgNode("polyline",{points:r.points.map(p=>p.join(",")).join(" "),fill:"none",stroke:r.layer==="top"?"#f19177":"#71afff","stroke-width":r.width,"stroke-linecap":"round","stroke-linejoin":"round",opacity:!this.net||byID.get(r.from).net===this.net?1:.16,class:hot.has(r.id)?"route-hot":"route"});
        line.append(svgNode("title",{},`${r.from} → ${r.to} · ${r.layer} · ${r.width} mm · ${d.currents[netOf(d,r)]??0} A`));svg.append(line);}
      for(const c of d.components) {
        const b=bodyBounds(c),g=this.object(svg,c,c.id);
        g.append(svgNode("rect",{...b,rx:.5,class:"pcb-body"}),svgNode("text",{x:c.x,y:c.y-.7,"font-size":2,"text-anchor":"middle",class:"component-label"},c.id));
      }
      for(const p of all) {
        const g=svgNode("g",{role:"button",tabindex:0,"aria-label":p.id+" "+(p.net||"unassigned"),class:"pcb-pad", "data-pad":p.id});
        const selected=this.startPad===p.id;
        g.append(svgNode("circle",{cx:p.x,cy:p.y,r:1.4,fill:"transparent"}),svgNode("circle",{cx:p.x,cy:p.y,r:p.radius,fill:selected?"#fff":"#dfba6b",stroke:selected?"#5fd4e8":"#2a3542","stroke-width":selected ? .45 : .15}),svgNode("circle",{cx:p.x,cy:p.y,r:.3,fill:"#13232d"}));
        g.append(svgNode("title",{},p.id+" · "+(p.net||"no net")));
        g.addEventListener("click",()=>this.padClick(p));
        g.addEventListener("keydown",e=>{if(e.key==="Enter"||e.key===" "){e.preventDefault();this.padClick(p);}});
        svg.append(g);
      }
    } else {
      svg.append(svgNode("rect",{x:d.margin,y:d.margin,width:d.width-2*d.margin,height:d.height-2*d.margin,class:"core-boundary"}));
      const colors={logic:"#65b3bd",memory:"#a493db",io:"#d4ae69",analog:"#d78d91"};
      const blocks=powerSummary(d,this.solution?.results.silicon?.details?.flows).blocks,peak=Math.max(1e-9,...blocks.map(b=>b.power/b.area));
      for(const b of d.blocks) {
        const g=this.object(svg,b,b.name),p=blocks.find(x=>x.id===b.id),density=p.power/p.area,heat=`hsl(${Math.round(220-200*density/peak)} 75% 55%)`;
        g.append(svgNode("rect",{x:b.x,y:b.y,width:b.width,height:b.height,fill:this.heat?heat:colors[b.kind],"fill-opacity":this.heat?.55:.22,stroke:colors[b.kind],"stroke-width":2}));
        g.append(svgNode("text",{x:b.x+12,y:b.y+28,"font-size":18,class:"component-label"},(b.fixed?"▣ ":"")+b.name),svgNode("text",{x:b.x+12,y:b.y+52,"font-size":12,class:"component-label"},this.heat?`${(p.power*1e3).toFixed(0)} mW · ${density.toFixed(2)} W/mm²${p.fromRtl?" · RTL":""}`:b.width+" × "+b.height+" µm"));
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
      if(!drag)return;const p=coordinates(e),distance=Math.hypot(p.x-drag.x,p.y-drag.y),was=drag;drag=null;
      if(distance<d.grid/2){this.selected=was.id;this.render();return;}
      if(was.locked){this.message("Unlock this block before moving it.");return;}
      const snap=v=>Math.round(v/d.grid)*d.grid,x=snap(was.startX+p.x-was.x),y=snap(was.startY+p.y-was.y);
      this.selected=was.id;this.run(this.pcb?"pcb.move":"silicon.update_block",this.pcb?{id:was.id,at:[x,y]}:{id:was.id,x,y},"Placement updated on the snap grid.");
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
    wrap.append(node("p","Logical net membership from pad assignments, with each net's RMS current for the copper-loss and current-density checks. This is a connection list, not a circuit simulation."));
    for(const net of pcbChecks(this.design).nets) {
      const row=node("section");row.append(node("strong",net.name),node("span",net.remaining+" explicit links remaining"));
      const amps=node("label","RMS current (A)"),a=node("input");a.type="number";a.min=0;a.max=500;a.step=.05;a.value=this.design.currents[net.name]??0;a.id="pcb-current-"+net.name;amps.htmlFor=a.id;amps.append(a);
      a.addEventListener("change",()=>{if(a.checkValidity()&&a.value!=="")this.run("pcb.current",{net:net.name,amps:Number(a.value)});});
      if(this.design.loadNets.includes(net.name)&&this.solution?.applied["pcb.load_current"])amps.append(node("small",`⛓ linked: ${format(this.solution.applied["pcb.load_current"].value,"A",4)}`,"link-note"));
      row.append(amps);
      const members=node("div");for(const p of pads(this.design).filter(p=>p.net===net.name))members.append(button(p.id,null,()=>{this.selected=p.component;this.connectivity=false;this.render();}));
      row.append(members);wrap.append(row);
    }return wrap;
  }
  problems() {
    const results=this.pcb?pcbChecks(this.design):siliconChecks(this.design),section=node("section",undefined,"studio-checks");
    const issues=this.pcb?[...results.issues,...densityIssues(this.design)]:results.issues;
    section.id=this.kind+"-checks";section.setAttribute("aria-label",this.pcb?"PCB geometry checks":"Floorplan placement checks");
    const top=node("div",undefined,"section-top");top.append(node("h3",this.pcb?"Board geometry checks":"Placement checks"));
    top.append(node("span",issues.length+" issues"+(this.pcb?" · "+results.unrouted+" explicit links remaining":" · "+results.utilization.toFixed(1)+"% requested core area")));section.append(top);
    if(!issues.length)section.append(node("p","No violations in the implemented geometry checks.","studio-check-pass"));
    else {
      const list=node("ul");
      for(const issue of issues.slice(0,50)){const li=node("li");li.dataset.code=issue.code;li.append(button(issue.message,null,()=>{
        if(this.items().some(i=>i.id===issue.object)){this.selected=issue.object;this.render();}
        else this.message(issue.message);
      }));list.append(li);}
      section.append(list);if(issues.length>50)section.append(node("p","Showing the first 50 issues."));
    }
    if(this.pcb)section.append(node("p","Connectivity is counted from explicit pad-to-pad links; same-net copper intersections are not inferred.","small muted"));
    this.root.append(section);
  }
  // Linked-system results: kernel numbers, checks and certificates.
  renderResults() {
    if(!this.results)return;
    const r=this.solution?.results[this.kind];this.results.replaceChildren();
    if(!this.pcb&&this.flow)this.renderFlow(r);
    if(!r){this.results.append(node("p","Solving the linked system…","small muted"));return;}
    const top=node("div",undefined,"section-top");top.append(node("h2",this.pcb?"Kernel board results":"Chip flow results"),node("span",`linked solve · ${this.solution.iterations} passes`));this.results.append(top);
    const stats=node("dl",undefined,"stat-grid"),stat=(t,v)=>{const d=node("div");d.append(node("dt",t),node("dd",v));stats.append(d);};
    if(this.pcb) {
      const k=r.details?.kernel;
      stat("Copper I²R loss",format(r.outputs.copper_loss?.value,"W",4));stat("Board + parts dissipation",format(r.outputs.dissipation?.value,"W",4));stat("Populated mass",format(r.outputs.board_mass?.value*1000,"g",3));
      stat("Load current (OUT, GND)",format(this.solution.applied["pcb.load_current"]?.value??this.design.currents.OUT,"A",4));
      if(k){stat("Sheet resistance",`${(k.sheet_resistance[0]*1e3).toFixed(4)} mΩ/□`);stat("Stack thermal resistance",`${(k.thermal_resistance*1e4).toFixed(3)} K·cm²/W`);}
      this.results.append(stats);
      if(k?.traces.length){const t=node("table",undefined,"data-table"),h=node("tr"),head=node("thead"),body=node("tbody");
        for(const c of ["Trace","Layer","Length","Width","R (40 °C)","I","Loss"])h.append(node("th",c));head.append(h);
        k.traces.forEach((x,i)=>{const route=this.design.routes[i];if(!route)return;const tr=node("tr");for(const c of [`${route.from} → ${route.to}`,route.layer,`${(x.length*1e3).toFixed(1)} mm`,`${route.width} mm`,`${(x.resistance*1e3).toFixed(2)} mΩ`,`${x.current} A`,format(x.loss,"W",3)])tr.append(node("td",c));body.append(tr);});
        t.append(head,body);this.results.append(t);}
    } else {
      const flows=r.details?.flows||{};
      stat("Die power",format(r.outputs.power?.value,"W",4));stat("Die junction",format(r.outputs.junction?.value,"degC",4));stat("Electrothermal passes",String(r.details?.iterations??"—"));stat("Slowest RTL fmax",format(r.outputs.fmax?.value,"",4)+" MHz");
      this.results.append(stats);
      const t=node("table",undefined,"data-table rtl-table"),head=node("thead"),hr=node("tr"),body=node("tbody");
      for(const c of ["Block","RTL","Cells","Transistors","Flip-flops","Equivalence","fmax","Setup slack","Power"])hr.append(node("th",c));head.append(hr);
      for(const [block,f] of Object.entries(flows)){const tr=node("tr");for(const c of [this.design.blocks.find(b=>b.id===block)?.name||block,f.module,Object.values(f.cells).reduce((a,b)=>a+b,0).toLocaleString(),f.transistors.toLocaleString(),String(f.flops),
        f.equivalence.mismatch?"MISMATCH":`${f.equivalence.points.toLocaleString()} pts${f.equivalence.exhaustive?" · proof":""}`,f.timing?`${(f.timing.fmax/1e6).toFixed(1)} MHz`:"—",f.timing?`${(f.timing.setup_slack*1e12).toFixed(0)} ps`:"—",f.power?`${(f.power.total*1e6).toFixed(3)} µW`:"—"])tr.append(node("td",c));body.append(tr);
        if(f.critical){const cr=node("tr"),td=node("td",`Critical path at ${f.critical.vdd} V, ${f.critical.temp} K to ${f.critical.endpoint}: ${f.critical.stages.map(s=>s.cell).join(" → ")}`,"critical-path");td.colSpan=9;cr.append(td);body.append(cr);}}
      t.append(head,body);this.results.append(t);
      this.rtl(this.results,flows);
    }
    checkCards(this.results,r.checks||[]);certificates(this.results,r.certificates);
  }
  renderFlow(r) {
    const flows=Object.values(r?.details?.flows||{}),f=flows[0],checks=r?.checks||[],ok=id=>checks.find(c=>c.id===id)?.status==="met"?"done":"fail";
    const state={rtl:f?"done":"wait",synth:f?"done":"wait",equiv:f?(flows.some(x=>x.equivalence.mismatch)?"fail":"done"):"wait",timing:f?(flows.every(x=>x.timing&&x.timing.setup_slack>=0&&x.timing.hold_slack>=0)?"done":"fail"):"wait",
      power:r?ok("tj"):"wait",floorplan:r?ok("floorplan"):"wait",litho:"study"};
    const detail={rtl:f?`${flows.length} design${flows.length>1?"s":""} · ${flows.map(x=>x.module.replace(/^tt_um_/,"")).join(", ")}`:"",synth:f?`${flows.reduce((s,x)=>s+x.transistors,0).toLocaleString()} transistors`:"",equiv:f?`${flows.map(x=>x.equivalence.points.toLocaleString()).join(" + ")} points`:"",
      timing:f?`fmax ${format(r.outputs.fmax?.value,"",3)} MHz`:"",power:r?format(r.outputs.junction?.value,"degC",3):"",floorplan:r?format(r.outputs.utilization?.value,"%",3):""};
    this.flow.replaceChildren();
    for(const s of STAGES){const d=node("div");d.dataset.state=s.planned?"planned":state[s.id]||"wait";d.append(node("strong",s.title),node("small",detail[s.id]||s.detail));this.flow.append(d);}
  }
  rtl(root,flows) {
    const d=this.design.chip,box=node("details",undefined,"rtl-panel");box.open=!!this.rtlOpen;box.addEventListener("toggle",()=>this.rtlOpen=box.open);
    const blocks=this.design.blocks,current=this.rtlBlock&&blocks.some(b=>b.id===this.rtlBlock)?this.rtlBlock:d.rtl[0]?.block??blocks[0]?.id,entry=d.rtl.find(r=>r.block===current);
    box.append(node("summary",`RTL sources · ${d.rtl.map(r=>`${blocks.find(b=>b.id===r.block)?.name}: ${r.source}`).join(" · ")||"none"}`));
    const pick=node("div",undefined,"rtl-actions"),bl=node("label","Block"),bs=node("select");bs.id="silicon-rtl-block";bl.htmlFor=bs.id;
    for(const b of blocks){const o=node("option",b.name);o.value=b.id;bs.append(o);}bs.value=current;bs.addEventListener("change",()=>{this.rtlBlock=bs.value;this.render();});bl.append(bs);pick.append(bl);
    pick.append(node("p",entry?RTL_SOURCES[entry.source]:"This block uses its estimated power.","small muted"));box.append(pick);
    const area=node("textarea");area.id="silicon-verilog";area.rows=14;area.spellcheck=false;area.value=entry?.source==="custom"?entry.verilog:project.rtlSources[entry?.source??"mini8"]??"";box.append(area);
    const row=node("div",undefined,"rtl-actions");
    row.append(button("Synthesize as custom RTL","silicon-synthesize",()=>this.run("silicon.rtl",{block:current,source:"custom",verilog:area.value}),"accent"),
      button("Use mini8","silicon-rtl-mini8",()=>this.run("silicon.rtl",{block:current,source:"mini8"})),button("Use Tiny Tapeout PWM","silicon-rtl-tt",()=>this.run("silicon.rtl",{block:current,source:"tt-motor"})));
    for(const [key,label,min,max,step] of [["frequency","Clock (MHz)",.001,5000,1],["vdd","Supply (V)",1.08,1.32,.01],["activity","Input activity (%)",0,100,1],["theta","Package θJA (K/W)",.1,500,1],["instances","Instances",1,100000,1]]){
      const l=node("label",label),i=node("input");i.type="number";i.min=min;i.max=max;i.step=step;i.value=d[key];i.id="silicon-chip-"+key;l.htmlFor=i.id;l.append(i);row.append(l);
      i.addEventListener("change",()=>{if(i.checkValidity()&&i.value!=="")this.run("silicon.chip",{[key]:Number(i.value)});});
    }
    const amb=this.solution?.applied["silicon.ambient"];row.append(node("p",amb?`⛓ Package ambient ${format(amb.value,"degC",3)} from ${amb.from.join(", ")}`:`Package ambient ${d.ambient} °C`,"link-note"));
    box.append(row);
    for(const [block,f] of Object.entries(flows||{}))if(f.equivalence?.mismatch)box.append(node("p",`${block}: synthesis does not match the RTL: ${f.equivalence.mismatch}`,"conflict"));
    box.append(node("p","The Tiny Tapeout controller is HADJAL001/tt-motor-controller (Apache-2.0) at commit 0f078f5, adapted to the native subset: its asynchronous reset is synchronous, and lint-only lines are removed. The header of data/tt_motor_controller.v lists every change.","small muted"));
    root.append(box);
  }
  async impedance(id) {
    const route=this.design.routes.find(r=>r.id===id),s=this.design.stackup,out=this.impedanceOut;if(!route)return;
    const h=(s.thickness-2*s.copper*OZ),t=s.copper*OZ;out.textContent="Solving the drawn cross-section in cad::laplace…";
    try{const {result:r}=await request("line",{input:{width:route.width*1e-3,height:h*1e-3,thickness:t*1e-3,eps_r:s.epsr,field:true,target:50}});
      out.textContent=`Microstrip over the opposite layer (h ${h.toFixed(3)} mm, εr ${s.epsr}): Hammerstad–Jensen ${r.formula.z0.toFixed(2)} Ω, εeff ${r.formula.eps_eff.toFixed(3)}${r.formula.valid?"":" (outside 0.05 ≤ W/h ≤ 20)"} · field solve with ${(t*1000).toFixed(0)} µm copper ${r.field.z0.toFixed(2)} Ω (${((r.field.z0/r.formula.z0-1)*100).toFixed(2)}%), εeff ${r.field.eps_eff.toFixed(3)} · 50 Ω needs ${(r.width_for_target*1e3).toFixed(3)} mm.`;
    }catch(e){out.textContent="Impedance unavailable: "+e.message;}
  }
  async fabrication() {
    const d=this.design;this.message("Writing copper, profile and drill files in the Rust kernel…");
    try{
      const k=kernelBoard(d),{result:r}=await request("pcb",{input:{board:k.board,currents:k.currents,fabrication:true,cell:5e-5}});
      const files=[],report=[];
      for(const l of r.gerber){const name=l.name==="top"?"emergence-F_Cu.gbr":"emergence-B_Cu.gbr";files.push({name,data:l.gerber});
        const c=gerberCheck(d,l.gerber,l.name==="top"?"top":"bottom");report.push(`${l.name} copper: ${c.loops} regions, read-back ${c.gerberArea.toFixed(1)} mm² vs drawn ${c.drawnArea.toFixed(1)} mm², ${(c.agreement*100).toFixed(1)}% agreement at a 0.1 mm raster`);}
      files.push({name:"emergence-Edge_Cuts.gbr",data:r.profile},{name:"emergence.drl",data:r.excellon},{name:"emergence.kicad_pcb",data:kicadPCB(d)},{name:"emergence-bom.csv",data:pcbBOM(d)},
        {name:"README.txt",data:`${project.store.docs.project.name}: board fabrication package\n\nGerber RS-274X (mm, 4.6) and Excellon written by the cad crate to their specifications. NOT verified against a CAM tool: inspect in a Gerber viewer before ordering.\nFootprints and drills are illustrative, not IPC-7351 land patterns.\n\nRead-back check:\n${report.join("\n")}\n`});
      download(zip(files),`${named()}-fabrication.zip`,"application/zip");
      this.message("Fabrication package written. "+report.join(" · ")+". Not verified against a CAM tool.");
    }catch(e){this.message("Fabrication files unavailable: "+e.message);}
  }
  exportGDS() {
    const bytes=gdsii(this.design),v=verifyGDS(this.design,bytes);
    if(!v.ok){this.message("GDSII read-back failed: "+v.problems.join(" "));return;}
    download(bytes,"emergence-floorplan.gds","application/octet-stream");
    this.message(`GDSII read back exactly: ${v.boundaries} boundaries, ${v.labels} labels, ${v.bytes} bytes. Concept floorplan on layers die 0, kinds 1–4, core 10, text 63; no PDK layer map.`);
  }
  exportSVG() {
    const svg=this.root.querySelector("svg").cloneNode(true),style=svgNode("style",{},".studio-board{fill:#12212b;stroke:#6c8998;stroke-width:2}.core-boundary{fill:none;stroke:#8998a5;stroke-width:2;stroke-dasharray:8 5}.component-label{fill:#e5edf4;font-family:sans-serif}");
    svg.prepend(style);svg.setAttribute("xmlns","http://www.w3.org/2000/svg");
    svg.append(svgNode("text",{x:0,y:-12,fill:"#e5edf4","font-size":12},"Concept floorplan · µm · no PDK / no sign-off"));
    return new XMLSerializer().serializeToString(svg);
  }
}

new Studio("pcb");new Studio("silicon");
