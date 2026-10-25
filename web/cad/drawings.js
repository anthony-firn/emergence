// Analytic orthographic views of the supported parametric parts. Dimensions
// use design values, never measurements from the rounded display/STL mesh.
const esc = v => String(v).replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&apos;"}[c]));
const n = v => Number(v.toFixed(3));
export function drawingDefaults(part) {
  return { number: part.id === "sink" ? "ET-MD-001" : `ET-${part.id.slice(0,8).toUpperCase()}`, title: part.name,
    material: part.id === "sink" ? "Aluminium 6063 / temper unspecified" : "Not specified", finish: "Not specified", tolerance: "",
    notes: "DO NOT SCALE DRAWING. Dimensions in millimetres.\nDeburring and edge requirements to be specified by designer.",
    projection: "third", revision: "A", prepared: "", checked: "", status: "Draft" };
}
export function partGeometry(params, part) {
  return part.id === "sink" ? { width:params.width,depth:params.length,height:params.base+params.finHeight,base:params.base,
    fins:params.fins,thickness:params.finThickness,shape:"sink" } : { width:part.width,depth:part.shape==="cylinder"?part.width:part.depth,height:part.height,shape:part.shape };
}
export function releaseProblems(meta, dirty) {
  const problems = [];
  if (dirty) problems.push("Save a revision to bind the drawing to committed geometry.");
  for (const [k,label] of [["number","drawing number"],["title","title"],["material","material specification"],["finish","finish"],["revision","revision"],["prepared","preparer"],["checked","checker"]]) {
    if (!meta[k]?.trim() || /not specified|unspecified/i.test(meta[k])) problems.push(`Complete ${label}.`);
  }
  if (!/^\d+(\.\d{1,4})?$/.test(meta.tolerance) || Number(meta.tolerance)<=0 || Number(meta.tolerance)>10) problems.push("Specify a general linear tolerance in mm (0 < tolerance ≤ 10).");
  if (meta.status !== "Released") problems.push("Sheet status is not Released.");
  return problems;
}
export function drawingSVG(params, part, meta, revision, dirty, markerId="drawing-arrow") {
  const g=partGeometry(params,part), fit=100/Math.max(g.width,g.depth,g.height), s=fit>=2?2:fit>=1?1:.5;
  const problems=releaseProblems(meta,dirty), state=dirty?"UNCOMMITTED DRAFT":meta.status.toUpperCase();
  const parts=[];
  const line=(x1,y1,x2,y2,cls="edge")=>parts.push(`<line class="${cls}" x1="${n(x1)}" y1="${n(y1)}" x2="${n(x2)}" y2="${n(y2)}"/>`);
  const rect=(x,y,w,h,cls="edge")=>parts.push(`<rect class="${cls}" x="${n(x)}" y="${n(y)}" width="${n(w)}" height="${n(h)}"/>`);
  const text=(x,y,t,size=3,anchor="start",cls="")=>parts.push(`<text class="${cls}" x="${n(x)}" y="${n(y)}" font-size="${size}" text-anchor="${anchor}">${esc(t)}</text>`);
  const hdim=(x1,x2,y,from,label)=>{line(x1,from,x1,y-2,"dim");line(x2,from,x2,y-2,"dim");line(x1,y,x2,y,"arrow");text((x1+x2)/2,y-1.7,label,3,"middle");};
  const vdim=(x,y1,y2,from,label)=>{line(from,y1,x-2,y1,"dim");line(from,y2,x-2,y2,"dim");line(x,y1,x,y2,"arrow");parts.push(`<text x="${n(x-2)}" y="${n((y1+y2)/2)}" text-anchor="middle" font-size="3" transform="rotate(-90 ${n(x-2)} ${n((y1+y2)/2)})">${esc(label)}</text>`);};
  const front=(x,y)=>{
    if(g.shape==="sink"){
      rect(x,y+(g.height-g.base)*s,g.width*s,g.base*s);
      for(let j=0;j<g.fins;j++){const fx=x+j*(g.width-g.thickness)/(g.fins-1)*s;rect(fx,y,g.thickness*s,(g.height-g.base)*s);}
      vdim(x+g.width*s+12,y+(g.height-g.base)*s,y+g.height*s,x+g.width*s,`${n(g.base)}`);
      line(x,y-2,x-5,y-9,"dim");line(x-5,y-9,x+g.width*s+5,y-9,"dim");
      text(x+g.width*s+6,y-8,`${g.fins}X ${n(g.thickness)} FIN`,2.8);
    }else{rect(x,y,g.width*s,g.height*s);if(g.shape==="cylinder")line(x+g.width*s/2,y-4,x+g.width*s/2,y+g.height*s+4,"center");}
    vdim(x-12,y,y+g.height*s,x,`${n(g.height)}`);text(x+g.width*s/2,y+g.height*s+9,"FRONT VIEW",2.6,"middle");
  };
  const top=(x,y)=>{
    if(g.shape==="cylinder"){
      parts.push(`<circle class="edge" cx="${n(x+g.width*s/2)}" cy="${n(y+g.width*s/2)}" r="${n(g.width*s/2)}"/>`);
      line(x-4,y+g.width*s/2,x+g.width*s+4,y+g.width*s/2,"center");line(x+g.width*s/2,y-4,x+g.width*s/2,y+g.width*s+4,"center");
    }else{
      rect(x,y,g.width*s,g.depth*s);
      if(g.shape==="sink")for(let j=0;j<g.fins;j++){const fx=x+j*(g.width-g.thickness)/(g.fins-1)*s;line(fx,y,fx,y+g.depth*s);line(fx+g.thickness*s,y,fx+g.thickness*s,y+g.depth*s);}
    }
    hdim(x,x+g.width*s,y-12,y,`${g.shape==="cylinder"?"⌀ ":""}${n(g.width)}`);
    if(g.shape!=="cylinder")vdim(x-12,y,y+g.depth*s,x,`${n(g.depth)}`);
    text(x+g.width*s/2,y+g.depth*s+9,"TOP VIEW",2.6,"middle");
  };
  const side=(x,y)=>{rect(x,y,g.depth*s,g.height*s);if(g.shape==="sink")line(x,y+(g.height-g.base)*s,x+g.depth*s,y+(g.height-g.base)*s);text(x+g.depth*s/2,y+g.height*s+9,"RIGHT VIEW",2.6,"middle");};
  rect(8,8,404,281);text(15,19,"EVOLVETECH / EMERGENCE",4.5);text(405,19,state,3,"end");
  if(meta.projection==="third"){top(60,48);front(60,175);side(240,175);}else{front(240,58);side(60,58);top(240,155);}
  const nx=meta.projection==="third"?235:25,ny=meta.projection==="third"?42:155;
  text(nx,ny,"MANUFACTURING NOTES",3.3);
  const notes=[`PART: ${part.name}`,`MATERIAL: ${meta.material}`,`FINISH: ${meta.finish}`,
    meta.tolerance?`GENERAL LINEAR TOLERANCE: ±${meta.tolerance} mm`:"GENERAL TOLERANCE: NOT SPECIFIED",
    `${meta.projection==="third"?"THIRD":"FIRST"} ANGLE PROJECTION`,"UNITS: mm. Written dimensions govern."];
  if(g.shape==="sink")notes.push(`${g.fins} FINS EQUALLY SPACED ACROSS WIDTH.`,`GAP (${n((g.width-g.fins*g.thickness)/(g.fins-1))}) REF.`,`CROSS-SECTION CONSTANT THROUGH ${n(g.depth)} LENGTH.`);
  const wrap=t=>{const result=[];for(const raw of t.split("\n")){let rest=raw;while(rest.length>66){let at=rest.lastIndexOf(" ",66);if(at<20)at=66;result.push(rest.slice(0,at));rest=rest.slice(at).trimStart();}result.push(rest);}return result;};
  notes.push(...wrap(meta.notes));
  notes.slice(0,22).forEach((t,i)=>text(nx,ny+7+i*4.4,t,2.6));
  text(16,258,problems.length?"ENGINEERING REVIEW REQUIRED":"RELEASED BY NAMED CHECKER",3.2);
  text(16,265,`CONTENT REVISION: ${revision.slice(0,16)}`,2.6);text(16,271,dirty?"GEOMETRY HAS UNSAVED CHANGES":"GEOMETRY AND DRAWING METADATA ARE REVISIONED",2.2);
  text(16,279,"Prototype drawing generator. Verify process-specific requirements.",2.2);
  rect(225,249,187,40);line(225,264,412,264);line(225,278,412,278);line(355,264,355,289);
  text(229,254,"TITLE",2);text(229,261,meta.title,3.3);
  text(229,268,`DWG NO. ${meta.number}`,2.5);text(359,269,`REV ${meta.revision}`,2.5);
  text(229,273,`PREPARED ${meta.prepared || "—"}`,2.2);text(229,277,`CHECKED ${meta.checked || "—"}`,2.2);
  text(229,283,`SCALE ${n(s)}:1 / A3 / mm`,2.5);text(359,283,"SHEET 1 OF 1",2.4);
  text(229,287,`${meta.projection.toUpperCase()} ANGLE / ${meta.status.toUpperCase()}`,2.1);
  // Projection symbol: truncated cone side view and circular end view.
  const cx=meta.projection==="third"?387:372,cy=234;
  parts.push(`<path class="edge" d="M ${cx-15} ${cy-4} L ${cx-5} ${cy-6} L ${cx-5} ${cy+6} L ${cx-15} ${cy+4} Z"/><circle class="edge" cx="${meta.projection==="third"?cx+4:cx-24}" cy="${cy}" r="5"/><circle class="edge" cx="${meta.projection==="third"?cx+4:cx-24}" cy="${cy}" r="3"/>`);
  return `<svg xmlns="http://www.w3.org/2000/svg" width="420mm" height="297mm" viewBox="0 0 420 297" id="${markerId}-sheet" role="img" aria-label="${esc(meta.title)} dimensioned manufacturing drawing"><title>${esc(meta.title)} / ${esc(meta.number)}</title><defs><marker id="${markerId}" markerWidth="6" markerHeight="5" refX="0" refY="2.5" orient="auto-start-reverse" markerUnits="strokeWidth"><path d="M6 0 L0 2.5 L6 5Z" fill="#111"/></marker></defs><style>#${markerId}-sheet text{font-family:Arial,Helvetica,sans-serif;fill:#111}#${markerId}-sheet .edge{stroke:#111;stroke-width:.35;fill:none}#${markerId}-sheet .dim{stroke:#111;stroke-width:.18;fill:none}#${markerId}-sheet .arrow{stroke:#111;stroke-width:.18;marker-start:url(#${markerId});marker-end:url(#${markerId})}#${markerId}-sheet .center{stroke:#111;stroke-width:.18;stroke-dasharray:8 2 1 2;fill:none}</style><rect width="420" height="297" fill="white"/>${parts.join("")}</svg>`;
}

export function profileDXF(params,part) {
  const g=partGeometry(params,part), entities=[];
  let handle=0x100;
  const entity=type=>`0\n${type}\n5\n${(handle++).toString(16).toUpperCase()}\n100\nAcDbEntity\n8\nPROFILE\n100\n${type==="LINE"?"AcDbLine":"AcDbCircle"}\n`;
  const line=(x1,y1,x2,y2)=>entities.push(`${entity("LINE")}10\n${n(x1)}\n20\n${n(y1)}\n30\n0\n11\n${n(x2)}\n21\n${n(y2)}\n31\n0\n`);
  if(g.shape==="cylinder")entities.push(`${entity("CIRCLE")}10\n0\n20\n0\n30\n0\n40\n${n(g.width/2)}\n`);
  else if(g.shape==="sink") {
    const pts=[[0,0],[g.width,0],[g.width,g.height]];
    for(let j=g.fins-1;j>=0;j--){const x=j*(g.width-g.thickness)/(g.fins-1);pts.push([x,g.height]);if(j>0)pts.push([x,g.base],[x-(g.width-g.fins*g.thickness)/(g.fins-1),g.base],[x-(g.width-g.fins*g.thickness)/(g.fins-1),g.height]);}
    pts.push([0,0]);for(let i=1;i<pts.length;i++)line(...pts[i-1],...pts[i]);
  } else {line(0,0,g.width,0);line(g.width,0,g.width,g.depth);line(g.width,g.depth,0,g.depth);line(0,g.depth,0,0);}
  return `0\nSECTION\n2\nHEADER\n9\n$ACADVER\n1\nAC1015\n9\n$INSUNITS\n70\n4\n9\n$MEASUREMENT\n70\n1\n0\nENDSEC\n0\nSECTION\n2\nTABLES\n0\nTABLE\n2\nLAYER\n5\nA\n100\nAcDbSymbolTable\n70\n1\n0\nLAYER\n5\nB\n100\nAcDbSymbolTableRecord\n100\nAcDbLayerTableRecord\n2\nPROFILE\n70\n0\n62\n7\n6\nCONTINUOUS\n0\nENDTAB\n0\nENDSEC\n0\nSECTION\n2\nENTITIES\n${entities.join("")}0\nENDSEC\n0\nEOF\n`;
}
export function bomCSV(params,features,drawings,revision) {
  const q=v=>{let s=String(v);if(/^[=+@-]/.test(s))s="'"+s;return`"${s.replaceAll('"','""')}"`;};
  const rows=[["Part","Quantity","Description","Material","Size mm","Drawing","Drawing revision","Content revision"],
    ["sink",1,"Finned heat sink",drawings.sink?.material || "Aluminium 6063 / temper unspecified",`${params.width} x ${params.length} x ${params.base+params.finHeight}`,drawings.sink?.number||"",drawings.sink?.revision||"",revision],
    ...features.map(f=>[f.id,1,f.name,drawings[f.id]?.material||"Not specified",f.shape==="cylinder"?`D${f.width} x ${f.height}`:`${f.width} x ${f.depth} x ${f.height}`,drawings[f.id]?.number||"",drawings[f.id]?.revision||"",revision])];
  return rows.map(r=>r.map(q).join(",")).join("\r\n");
}
