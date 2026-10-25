const $ = (id) => document.getElementById(id);
const status = $("status");
const sci = (v) => Number.isFinite(v) ? v.toExponential(2) : "—";
const ns = "http://www.w3.org/2000/svg";
const svg = (tag, attrs, content) => { const e = document.createElementNS(ns, tag); for (const [k,v] of Object.entries(attrs)) e.setAttribute(k,v); if(content !== undefined)e.textContent=content; return e; };
let worker, id = 0;
const pending = new Map();
function fail(error) {
  status.textContent = `Models could not run: ${error.message}. Reload to retry, or read the source and validation results.`;
  for(const {reject,timer} of pending.values()){clearTimeout(timer);reject(error);} pending.clear();
}
try {
  worker = new Worker(new URL("worker.js", import.meta.url), { type: "module" });
  worker.onmessage = ({data}) => { const p = pending.get(data.id); if(!p)return; clearTimeout(p.timer); pending.delete(data.id); data.error ? p.reject(new Error(data.error)) : p.resolve(data); };
  worker.onerror = (e) => fail(new Error(e.message || "Worker could not start"));
} catch(error) { fail(error); }
function request(type, args = {}) {
  if(!worker)return Promise.reject(new Error("Web Workers are unavailable"));
  return new Promise((resolve,reject) => {
    const key = ++id;
    const timer = setTimeout(() => { worker.terminate(); worker=null; fail(new Error("The model exceeded its 45-second time limit")); },45000);
    pending.set(key,{resolve,reject,timer}); worker.postMessage({id:key,type,...args});
  });
}
let tab="optics", flowReady=false, playing=false, busy=false, flowGeneration=0, renderGeneration=0;
const tabs=[...document.querySelectorAll("[data-tab]")];
function switchTab(button){
  if(tab === "flow") pause();
  tab=button.dataset.tab;
  for(const b of tabs){const on=b===button;b.setAttribute("aria-selected",String(on));b.tabIndex=on?0:-1;$(b.dataset.tab).hidden=!on;}
  if(tab==="flow"&&!flowReady)initFlow().catch(fail);
  if(tab==="design"&&!$("design-results").hasChildNodes())evaluateDesign().catch(fail);
}
for(const [j,b] of tabs.entries()){
  b.addEventListener("click",()=>switchTab(b));
  b.addEventListener("keydown",e=>{const target=e.key==="ArrowRight"?tabs[(j+1)%3]:e.key==="ArrowLeft"?tabs[(j+2)%3]:e.key==="Home"?tabs[0]:e.key==="End"?tabs[2]:null;if(target){e.preventDefault();switchTab(target);target.focus();}});
}

let opticsGeneration=0, opticsTimer;
function showOptics({point,curve}){
  const chart=$("optics-chart");chart.replaceChildren();
  const m={l:52,r:20,t:24,b:55},w=700,h=400;
  const max=Math.max(0.1,...curve.map(p=>p.reference))*1.13;
  const x=v=>m.l+(Math.log10(v)+2)/(Math.log10(12)+2)*(w-m.l-m.r);
  const y=v=>h-m.b-v/max*(h-m.t-m.b);
  for(let j=0;j<=4;j++){const v=max*j/4;chart.append(svg("line",{x1:m.l,x2:w-m.r,y1:y(v),y2:y(v),class:"grid"}),svg("text",{x:m.l-9,y:y(v)+4,"text-anchor":"end"},v.toFixed(1)));}
  for(const v of [.01,.1,1,10])chart.append(svg("text",{x:x(v),y:h-m.b+22,"text-anchor":"middle"},String(v)));
  chart.append(svg("text",{x:m.l,y:12},"↑ extinction efficiency Qext"),svg("text",{x:w-m.r,y:h-6,"text-anchor":"end"},"particle size x → (log scale)"));
  chart.append(svg("path",{class:"reference",d:curve.map((p,j)=>`${j?"L":"M"}${x(p.x)},${y(p.reference)}`).join(" ")}));
  const colors=["#ff9b4a","#5fd4e8","#a78bfa"];
  for(const p of curve)chart.append(svg("circle",{cx:x(p.x),cy:y(p.q),r:3,fill:colors[p.rung]}));
  chart.append(svg("line",{x1:x(point.x),x2:x(point.x),y1:m.t,y2:h-m.b,class:"marker"}),svg("circle",{cx:x(point.x),cy:y(point.q),r:7,fill:colors[point.rung],stroke:"var(--text)","stroke-width":2}));
  const box=$("optics-result");box.replaceChildren();
  const name=document.createElement("strong");name.textContent=["Rayleigh expansion","Folded closure table","Full Mie series"][point.rung];box.append(name);
  const dl=document.createElement("dl");
  for(const [label,value] of [["Qext",point.q.toFixed(6)],["Claimed error",sci(point.claim)],["Observed error",sci(point.error)],["Rigor",point.rigor],["Request",point.met?"Met":"Unmet"]]){
    const row=document.createElement("div"),dt=document.createElement("dt"),dd=document.createElement("dd");dt.textContent=label;dd.textContent=value;if(label==="Request"&&!point.met)dd.className="unmet";row.append(dt,dd);dl.append(row);
  }box.append(dl);
}
async function solveOptics(){
  const gen=++opticsGeneration;const x=10**Number($("size").value);$("size-value").value=x.toFixed(2);
  status.textContent="Checking the optical ladder against its reference…";
  const response=await request("optics",{material:Number($("material").value),x,tolerance:Number($("tolerance").value)});
  if(gen!==opticsGeneration)return;if(response.result.error)throw new Error(response.result.error);
  showOptics(response.result);status.textContent=`Rust models ready · last optical solve + 100-point audit: ${response.ms.toFixed(0)} ms on this device`;
}
for(const name of ["size","material","tolerance"])$(name).addEventListener("input",()=>{
  $("size-value").value=(10**Number($("size").value)).toFixed(2);clearTimeout(opticsTimer);opticsTimer=setTimeout(()=>solveOptics().catch(fail),150);
});

function flowArgs(){const n=Number($("resolution").value);return {n,exact:Number($("flow-preset").value),re:Number($("reynolds").value),z:Math.floor(Number($("slice").value)*n/100)};}
function pause(){playing=false;$("play").textContent="Start flow";$("play").setAttribute("aria-pressed","false");}
async function initFlow(){
  pause();flowReady=false;const gen=++flowGeneration;$("play").disabled=true;$("reset").disabled=true;
  const args=flowArgs();$("flow-grid").textContent=`${args.n}³ cells`;
  const r=await request("flow-init",args);if(gen!==flowGeneration)return;
  flowReady=true;drawFlow(r.result,args.n);$("play").disabled=false;$("reset").disabled=false;
}
const canvas=$("flow-canvas"),ctx=canvas.getContext("2d");
const offscreen=document.createElement("canvas"),off=offscreen.getContext("2d");
function drawFlow(r,n){
  const field=r.field;offscreen.width=n;offscreen.height=n;const im=off.createImageData(n,n);
  for(let j=0;j<n;j++)for(let i=0;i<n;i++){
    const q=4*(j*n+i),p=4*((n-1-j)*n+i),v=Math.tanh(field[q+2]*.8),a=Math.abs(v);
    const color=v>=0?[255,155,74]:[95,212,232];
    for(let c=0;c<3;c++)im.data[p+c]=[9,16,28][c]*(1-a)+color[c]*a;im.data[p+3]=255;
  }off.putImageData(im,0,0);ctx.imageSmoothingEnabled=true;ctx.drawImage(offscreen,0,0,640,640);
  ctx.strokeStyle="rgba(240,245,255,.65)";ctx.lineWidth=1.1;
  for(let j=1;j<n;j+=2)for(let i=1;i<n;i+=2){
    const q=4*(j*n+i),x=i/n*640,y=(1-j/n)*640,dx=field[q]*16,dy=-field[q+1]*16;
    if(Math.hypot(dx,dy)<1)continue;ctx.beginPath();ctx.moveTo(x-dx*.4,y-dy*.4);ctx.lineTo(x+dx*.6,y+dy*.6);
    const a=Math.atan2(dy,dx);ctx.moveTo(x+dx*.6-4*Math.cos(a-.5),y+dy*.6-4*Math.sin(a-.5));ctx.lineTo(x+dx*.6,y+dy*.6);ctx.lineTo(x+dx*.6-4*Math.cos(a+.5),y+dy*.6-4*Math.sin(a+.5));ctx.stroke();
  }
  const metrics=$("flow-result");metrics.replaceChildren();
  for(const [label,value] of [["Simulation time",r.time.toFixed(3)],["Mean kinetic energy",r.energy.toFixed(6)],["Divergence residual",sci(r.divergence)],["Outer-shell energy",(r.tail*100).toFixed(3)+"%"]]){
    const row=document.createElement("div"),dt=document.createElement("dt"),dd=document.createElement("dd");dt.textContent=label;dd.textContent=value;row.append(dt,dd);metrics.append(row);
  }
  $("flow-validity").textContent=r.relative_error!==null?`Exact planar energy: ${r.exact_energy.toFixed(8)} · measured relative error ${sci(r.relative_error)}`:r.tail>.01?"Resolution warning: more than 1% of energy is near the retained grid cutoff. Try a finer grid or more viscosity.":"Outer-shell energy monitors unresolved scales. Choose the planar flow for an exact energy check.";
  canvas.dataset.time=String(r.time);canvas.dataset.error=String(r.relative_error);
}
async function tick(){
  if(!playing||busy||document.hidden||tab!=="flow")return;busy=true;const gen=flowGeneration,args=flowArgs();
  try{const r=await request("flow-step",args);if(gen===flowGeneration){drawFlow(r.result,args.n);status.textContent=`Rust fluid step + diagnostics: ${r.ms.toFixed(0)} ms on this device`;}}
  catch(e){pause();fail(e);}finally{busy=false;if(playing)setTimeout(tick,16);}
}
$("play").addEventListener("click",()=>{if(playing){pause();return;}playing=true;$("play").textContent="Pause flow";$("play").setAttribute("aria-pressed","true");tick();});
$("reset").addEventListener("click",()=>initFlow().catch(fail));
for(const k of ["flow-preset","resolution","reynolds"])$(k).addEventListener("change",()=>initFlow().catch(fail));
$("slice").addEventListener("input",async()=>{
  $("slice-value").value=`${$("slice").value}%`;if(!flowReady)return;
  const gen=++renderGeneration,args=flowArgs();try{const r=await request("flow-frame",args);if(gen===renderGeneration)drawFlow(r.result,args.n);}catch(e){fail(e);}
});
document.addEventListener("visibilitychange",()=>{if(document.hidden)pause();});

function drawDesign(){
  const scene=$("design-scene"),n=Number($("fins").value);scene.replaceChildren();
  scene.append(svg("path",{d:"M170 245 L480 245 L550 200 L240 200 Z",fill:"#35475a",stroke:"#6a8299"}),svg("path",{d:"M170 245 L480 245 L480 260 L170 260 Z",fill:"#253449"}));
  for(let j=0;j<n;j++){const x=180+j*290/(n-1);scene.append(svg("path",{d:`M${x} 240 L${x} 105 L${x+65} 60 L${x+65} 195 Z`,fill:"#537189",stroke:"#89a5b9","stroke-width":1}),svg("path",{d:`M${x} 105 L${x+5} 105 L${x+70} 60 L${x+65} 60 Z`,fill:"#9bb1bf"}));}
  for(let j=0;j<3;j++){const y=130+j*34;scene.append(svg("path",{d:`M35 ${y} L135 ${y} M123 ${y-6} L135 ${y} L123 ${y+6}`,stroke:"#5fd4e8","stroke-width":2,fill:"none"}));}
  scene.append(svg("text",{x:35,y:110},`${$("velocity").value} m/s`),svg("text",{x:180,y:295},`${n} fins · 40 × 50 mm base · 20 mm fin height`));
  $("fins-value").value=String(n);$("velocity-value").value=`${$("velocity").value} m/s`;
}
async function evaluateDesign(){
  const button=$("design-run");button.disabled=true;button.textContent="Evaluating…";drawDesign();
  try{
    const r=await request("design",{fins:Number($("fins").value),velocity:Number($("velocity").value)});
    if(r.result.error)throw new Error(r.result.error);
    const box=$("design-results");box.replaceChildren();
    for(const c of r.result.checks){
      const card=document.createElement("article");card.className="check";card.dataset.status=c.status;
      const top=document.createElement("div");top.className="check-top";const name=document.createElement("b"),badge=document.createElement("span");name.textContent={"REQ-TJ":"Junction temperature ≤ 100 °C","REQ-DROP":"PCB trace drop ≤ 50 mV","REQ-FAN":"Heat-sink pressure drop ≤ 20 Pa"}[c.id];badge.textContent=c.status;badge.className="check-status";top.append(name,badge);card.append(top);
      const detail=document.createElement("p");
      if(c.value){const v=c.value,offset=c.id==="REQ-TJ"?273.15:0,scale=c.id==="REQ-DROP"?1000:1,unit=c.id==="REQ-TJ"?"°C":c.id==="REQ-DROP"?"mV":"Pa";detail.className="bounds";detail.textContent=`${((v.nominal-offset)*scale).toFixed(2)} ${unit} · [${((v.lo-offset)*scale).toFixed(2)}, ${((v.hi-offset)*scale).toFixed(2)}] · ${v.rigor}`;}
      else detail.textContent=c.reason||"Model unavailable for these settings";
      card.append(detail);box.append(card);
    }status.textContent=`Rust design evaluation: ${r.ms.toFixed(0)} ms on this device · ${r.result.computed} model runs`;
  }finally{button.disabled=false;button.textContent="Evaluate design";}
}
for(const k of ["fins","velocity"])$(k).addEventListener("input",()=>{drawDesign();$("design-results").replaceChildren();});
$("design-run").addEventListener("click",()=>evaluateDesign().catch(fail));
drawDesign();
request("ready").then(()=>{$("design-run").disabled=false;return solveOptics();}).catch(fail);
