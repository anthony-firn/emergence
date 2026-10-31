//! Small, dependency-free JavaScript ABI. Models run unchanged in a Worker.
//! JSON and slice pointers stay valid until the next call that replaces them.
use emergence::{Ladder, Request};
use num_complex::Complex64;
use photonics::{
    mie_ladder::{self, MieInput, MieQ},
    C,
};
use serde_json::{json, Value};
use spectral_core::{fft::Fft3, init, Grid, Solver};
use std::cell::RefCell;

thread_local! {
    static OUTPUT: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
    static FLOW: RefCell<Option<Flow>> = const { RefCell::new(None) };
    static OPTICS: RefCell<Option<(u32, Ladder<MieInput, MieQ>)>> = const { RefCell::new(None) };
    static INPUT: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
    static MESH: RefCell<(Vec<f32>, Vec<u32>)> = const { RefCell::new((Vec::new(), Vec::new())) };
    // Synthesis and the equivalence check depend only on the source; the
    // electrothermal loop re-times the same netlist at new temperatures.
    static CHIP: RefCell<Vec<ChipCache>> = const { RefCell::new(Vec::new()) };
}

fn output(v: Value) {
    OUTPUT.with(|b| *b.borrow_mut() = serde_json::to_vec(&v).unwrap());
}
#[no_mangle]
pub extern "C" fn result_ptr() -> *const u8 {
    OUTPUT.with(|b| b.borrow().as_ptr())
}
#[no_mangle]
pub extern "C" fn result_len() -> usize {
    OUTPUT.with(|b| b.borrow().len())
}

// Optical constants: estimates representative of these materials at the
// stated wavelengths, as in photonics/tests/mie_ladder.rs. These are fixed
// monochromatic sphere models, not broadband sky or cloud forecasts.
fn index(material: u32) -> C {
    match material {
        1 => C::new(1.95, 0.79),  // Bond & Bergstrom (2006), soot at 550 nm.
        2 => C::new(1.333, 2e-9), // Hale & Querry (1973), water at 550 nm.
        _ => C::new(1.53, 0.008), // Estimate: mineral dust, d'Almeida/OPAC range.
    }
}

fn optics_record(l: &Ladder<MieInput, MieQ>, material: u32, x: f64, tol: f64) -> Value {
    let i = MieInput {
        x,
        m: index(material),
    };
    let s = l.solve(&i, Request::new(tol));
    let reference = mie_ladder::exact(x, i.m);
    json!({"x":x,"q":s.value.q_ext,"reference":reference.q_ext,
        "error":mie_ladder::metric(&s.value,&reference),"claim":s.claim.error,
        "rigor":format!("{:?}",s.claim.rigor),"rung":s.rung,"met":s.met})
}

#[no_mangle]
pub extern "C" fn optics_solve(material: u32, x: f64, tol: f64) {
    if material > 2
        || !x.is_finite()
        || !(0.01..=12.0).contains(&x)
        || !tol.is_finite()
        || tol <= 0.0
    {
        output(json!({"error":"Invalid optical input"}));
        return;
    }
    OPTICS.with(|slot| {
        let mut state = slot.borrow_mut();
        if state.as_ref().is_none_or(|s| s.0 != material) {
            // Fold the existing reference into a resonance-aware table, once.
            let c = mie_ladder::closure_rung(index(material), 0.05, 12.0, 1e-5, 512);
            *state = Some((material, mie_ladder::ladder(vec![c], 1)));
        }
        let l = &state.as_ref().unwrap().1;
        let point = optics_record(l, material, x, tol);
        let curve: Vec<_> = (0..100)
            .map(|j| {
                let px = 10f64.powf(-2.0 + (12f64.log10() + 2.0) * j as f64 / 99.0);
                optics_record(l, material, px.min(12.0), tol)
            })
            .collect();
        output(json!({"point":point,"curve":curve}));
    });
}

pub struct Flow {
    pub solver: Solver,
    exact: bool,
    field: Vec<f32>,
}
impl Flow {
    pub fn new(n: usize, exact: bool, re: f64) -> Self {
        let g = Grid::new(n);
        let u = if exact {
            init::taylor_green_2d(&g)
        } else {
            init::taylor_green(&g)
        };
        Self {
            solver: Solver::new(g, 1.0 / re, u),
            exact,
            field: vec![0.0; 4 * n * n],
        }
    }
    pub fn advance(&mut self) {
        let dt = self.solver.cfl_dt(0.3).min(0.04);
        self.solver.step(dt);
    }
    fn frame(&mut self, z: usize) -> *const f32 {
        let s = &mut self.solver;
        let n = s.grid.n;
        let u = s.velocity_phys(0);
        let v = s.velocity_phys(1);
        let mut w = s.u[0].clone();
        for (p, w) in w.iter_mut().enumerate() {
            let k = s.grid.k(p);
            *w = Complex64::i() * (s.u[1][p] * k[0] - s.u[0][p] * k[1]);
        }
        let mut omega = vec![0.0; s.grid.real_len()];
        Fft3::new(&s.grid).inverse_inplace(&mut w, &mut omega);
        for j in 0..n {
            for i in 0..n {
                let p = (i * n + j) * n + z.min(n - 1);
                let q = 4 * (j * n + i);
                self.field[q] = u[p] as f32;
                self.field[q + 1] = v[p] as f32;
                self.field[q + 2] = omega[p] as f32;
                self.field[q + 3] = u[p].hypot(v[p]) as f32;
            }
        }
        self.field.as_ptr()
    }
    pub fn report(&self) -> Value {
        let d = self.solver.diagnostics();
        // Taylor & Green (1937): this planar solution has
        // E(t) = (1/4) exp(-4 nu t), on the 2pi-periodic box.
        let reference = 0.25 * (-4.0 * self.solver.nu * d.t).exp();
        json!({"time":d.t,"energy":d.energy,"enstrophy":d.enstrophy,
            "divergence":d.max_div,"tail":d.tail_energy_fraction,
            "exact_energy":if self.exact {Some(reference)}else{None},
            "relative_error":if self.exact {Some((d.energy-reference).abs()/reference)}else{None}})
    }
}

#[no_mangle]
pub extern "C" fn flow_init(n: u32, exact: u32, re: f64) -> u32 {
    if ![16, 24, 32].contains(&n) || !re.is_finite() || !(20.0..=1600.0).contains(&re) || exact > 1
    {
        return 0;
    }
    FLOW.with(|f| *f.borrow_mut() = Some(Flow::new(n as usize, exact == 1, re)));
    1
}
#[no_mangle]
pub extern "C" fn flow_step() {
    FLOW.with(|f| {
        if let Some(s) = f.borrow_mut().as_mut() {
            s.advance();
        }
    });
}
#[no_mangle]
pub extern "C" fn flow_frame(z: u32) -> *const f32 {
    FLOW.with(|f| {
        f.borrow_mut()
            .as_mut()
            .map_or(std::ptr::null(), |s| s.frame(z as usize))
    })
}
#[no_mangle]
pub extern "C" fn flow_report() {
    FLOW.with(|f| {
        output(
            f.borrow()
                .as_ref()
                .map_or(json!({"error":"Initialize a flow first"}), |s| s.report()),
        )
    });
}

fn design(fins: u32, velocity: f64) -> Result<Value, String> {
    workspace_design(fins, velocity, 30.0, 45.0, 0.0)
}

fn workspace_design(fins: u32, velocity: f64, current: f64, ambient: f64, heat: f64) -> Result<Value, String> {
    workspace_config(&json!({"fins":fins,"velocity":velocity,"current":current,"ambient":ambient,"heat":heat}))
}

fn setting(p: &Value, name: &str, default: f64, min: f64, max: f64) -> Result<f64, String> {
    let value = match p.get(name) { Some(v) => v.as_f64().ok_or_else(|| format!("Invalid {name}"))?, None => default };
    if !value.is_finite() || value < min || value > max { return Err(format!("{name} outside the supported range")); }
    Ok(value)
}

fn workspace_config(input: &Value) -> Result<Value, String> {
    let p = input.get("params").unwrap_or(input);
    let fins = setting(p,"fins",9.,5.,14.)?;
    if fins.fract() != 0. { return Err("Fin count must be an integer".into()); }
    let velocity = setting(p,"velocity",2.5,0.5,4.)?;
    let current = setting(p,"current",30.,0.5,40.)?;
    let ambient = setting(p,"ambient",45.,15.,60.)?;
    let heat = setting(p,"heat",0.,0.,10.)?;
    let width = setting(p,"width",40.,30.,70.)?;
    let length = setting(p,"length",50.,30.,90.)?;
    let height = setting(p,"finHeight",20.,10.,35.)?;
    let base = setting(p,"base",3.,2.,8.)?;
    let thickness = setting(p,"finThickness",1.2,0.8,2.)?;
    let trace = setting(p,"traceWidth",10.,4.,15.)?;
    let optical = setting(p,"opticalPower",0.,0.,5.)? * setting(p,"absorption",50.,0.,100.)? / 100.;
    let gap = (width - fins * thickness) / (fins - 1.);
    if gap <= 0. { return Err("Fins overlap: increase width or reduce count/thickness".into()); }
    let mut doc: project::Doc =
        serde_json::from_str(include_str!("../../../projects/motor-driver/project.json"))
            .map_err(|e| e.to_string())?;
    // This browser view evaluates the original thermal/electrical graph.
    // CAD, printing and signal-integrity examples remain in the native repo.
    doc.parts.clear();
    doc.analyses.remove("enc");
    doc.requirements
        .retain(|id, _| ["REQ-TJ", "REQ-DROP", "REQ-FAN"].contains(&id.as_str()));
    for (key, value, unit) in [("sink.fin_count",fins,""),("duct.velocity",velocity,"m/s"),
        ("q1.current",current,"A"),("ambient",ambient,"degC"),("sink.width",width,"mm"),
        ("sink.length",length,"mm"),("sink.fin_height",height,"mm"),("sink.base",base,"mm"),
        ("sink.fin_thickness",thickness,"mm"),("trace.width",trace,"mm")] {
        doc.params.get_mut(key).unwrap().value = format!("{value} {unit}");
    }
    // Prescribed extra watts deposited at Q1, not a radiation transport model.
    // Retains the original switching-loss interval, conservatively scaling it.
    // Optical coupling is prescribed absorbed power P*a, not inferred from
    // the independent sphere-scattering experiment. Conservation is exact
    // for this boundary model; the original loss interval is retained.
    // A linked power stage supplies its own on-resistance and switching loss.
    let switching = setting(p,"switching",3.,0.,50.)?;
    doc.params.get_mut("q1.switching_loss").unwrap().value = format!("{} W", switching + heat + optical);
    let rds = setting(p,"rdsOn",4.5,0.1,100.)?;
    if rds != 4.5 {
        let param = doc.params.get_mut("q1.rds_on_25").unwrap();
        param.value = format!("{rds} mOhm");
        param.tol = None;
    }
    let p = project::Project::new(doc, ".", Default::default())?;
    let e = p.evaluate();
    let checks: Vec<_> = e.checks.values().map(|c|json!({
        "id":c.id,"text":c.text,"status":c.status.name(),"reason":c.reason,
        "value":c.value.map(|v|json!({"nominal":v.nom,"lo":v.lo,"hi":v.hi,"rigor":format!("{:?}",v.rigor)})),
        "unit":c.unit,"owner":c.owner
    })).collect();
    let nodes: serde_json::Map<String, Value> = ["q1.power","q1.tj","sink.thermal_resistance","trace.resistance","sink.pressure_drop"]
        .iter().map(|id| ((*id).to_string(), e.value(id).map(|v| json!({"nominal":v.nom,"lo":v.lo,"hi":v.hi})).unwrap_or(Value::Null))).collect();
    Ok(json!({"checks":checks,"computed":e.computed,"cached":e.cached,"nodes":nodes,
        "handoffs":{"optical_heat_w":optical,"external_heat_w":heat,"gap_mm":gap,
        "sink_volume_mm3":width*length*base+fins*thickness*length*height}}))
}
#[no_mangle]
pub extern "C" fn design_solve(fins: u32, velocity: f64) {
    output(design(fins, velocity).unwrap_or_else(|e| json!({"error":e})));
}

#[no_mangle]
pub extern "C" fn workspace_solve(fins: u32, velocity: f64, current: f64, ambient: f64, heat: f64) {
    output(workspace_design(fins, velocity, current, ambient, heat).unwrap_or_else(|e| json!({"error":e})));
}

fn sink_mesh(fins: u32) -> Result<Value, String> {
    mesh_config(&json!({"params":{"fins":fins}}))
}

fn mesh_config(input: &Value) -> Result<Value, String> {
    let p = input.get("params").unwrap_or(input);
    let target = input.get("target").and_then(Value::as_str).unwrap_or("sink");
    let (solid, cell) = if target == "sink" {
        workspace_config(input)?;
        let n = setting(p,"fins",9.,5.,14.)? as u32;
        let w = setting(p,"width",40.,30.,70.)?;
        let l = setting(p,"length",50.,30.,90.)?;
        let h = setting(p,"finHeight",20.,10.,35.)?;
        let b = setting(p,"base",3.,2.,8.)?;
        let t = setting(p,"finThickness",1.2,0.8,2.)?;
        let mut parts = vec![cad::Solid::cuboid([0.,0.,b/2.],[w/2.,l/2.,b/2.])];
        for i in 0..n {
            let x = -(w-t)/2. + (w-t)*i as f64/(n-1) as f64;
            parts.push(cad::Solid::cuboid([x,0.,b+h/2.],[t/2.,l/2.,h/2.]));
        }
        let gap=(w-n as f64*t)/(n-1) as f64;
        let cell=(t/2.).min(gap/3.).min(0.6);
        if cell<0.2 { return Err("This fin gap requires a finer mesh than the browser export supports. Export the analytic DXF profile or widen the channels.".into()); }
        (cad::Solid::union(parts), cell)
    } else {
        let features = input.get("features").and_then(Value::as_array).ok_or("No custom parts")?;
        if features.len() > 16 { return Err("Too many custom parts".into()); }
        let f = features.iter().find(|f| f.get("id").and_then(Value::as_str) == Some(target)).ok_or("Part not found")?;
        let x = setting(f,"x",0.,-100.,100.)?;
        let y = setting(f,"y",0.,-100.,100.)?;
        let z = setting(f,"z",0.,-100.,100.)?;
        let w = setting(f,"width",10.,1.,50.)?;
        let d = setting(f,"depth",10.,1.,50.)?;
        let h = setting(f,"height",10.,1.,50.)?;
        let solid = match f.get("shape").and_then(Value::as_str) {
            Some("box") => cad::Solid::cuboid([x,y,z+h/2.],[w/2.,d/2.,h/2.]),
            Some("cylinder") => cad::Solid::cylinder_z([x,y,z+h/2.],w/2.,h),
            _ => return Err("Unsupported part shape".into())
        };
        (solid, (w.min(d).min(h)/12.).clamp(0.1,0.6))
    };
    let bounds=solid.bbox();
    let nodes=(0..3).map(|k|((bounds.1[k]-bounds.0[k])/cell).ceil()+5.).product::<f64>();
    if nodes>2_000_000. { return Err("This mesh exceeds the browser export budget. Use the analytic profile DXF or reduce part size/detail.".into()); }
    let mesh = cad::mesh_solid(&solid, cell)?;
    let r = mesh.report();
    Ok(json!({"vertices":mesh.vertices,"triangles":mesh.triangles,
        "cell_mm":cell,"watertight":r.watertight,"volume_mm3":r.volume}))
}

#[no_mangle]
pub extern "C" fn workspace_input(len: u32) -> *mut u8 {
    if len == 0 || len > 1_048_576 { return std::ptr::null_mut(); }
    INPUT.with(|b| { let mut b = b.borrow_mut(); b.resize(len as usize,0); b.as_mut_ptr() })
}

#[no_mangle]
pub extern "C" fn workspace_run(mesh: u32) {
    let result = INPUT.with(|b| {
        let input: Value = serde_json::from_slice(&b.borrow()).map_err(|e|e.to_string())?;
        match mesh { 0 => workspace_config(&input), 1 => mesh_config(&input), 2 => model_config(&input),
            3 => pcb_config(&input), 4 => line_config(&input), 5 => chip_config(&input), _ => Err("Unknown model request".into()) }
    });
    output(result.unwrap_or_else(|e: String| json!({"error":e})));
}

#[no_mangle]
pub extern "C" fn cad_mesh(fins: u32) {
    output(sink_mesh(fins).unwrap_or_else(|e| json!({"error":e})));
}

// Feature-tree models: the same JSON format and kits as `fluid cad`, built by
// cad::model. File-backed mesh features are refused: the browser has no disk.
fn refuse_files(v: &Value) -> Result<(), String> {
    match v {
        Value::Object(o) => {
            if o.get("op").and_then(Value::as_str) == Some("mesh") { return Err("File-backed mesh features are not available in the browser".into()); }
            o.values().try_for_each(refuse_files)
        }
        Value::Array(a) => a.iter().try_for_each(refuse_files),
        _ => Ok(()),
    }
}

fn model_config(input: &Value) -> Result<Value, String> {
    let m = input.get("model").ok_or("No model")?;
    refuse_files(m)?;
    let model: cad::model::Model = serde_json::from_value(m.clone()).map_err(|e| format!("model JSON: {e}"))?;
    let mut overrides = cad::Params::default();
    if let Some(o) = input.get("overrides").and_then(Value::as_object) {
        for (k, v) in o {
            let x = v.as_f64().filter(|x| x.is_finite()).ok_or_else(|| format!("override {k}: not a finite number"))?;
            overrides.set(k, x);
        }
    }
    let params = model.resolve(&overrides)?;
    let solid = model.build(&overrides, ".")?;
    let (lo, hi) = solid.bbox();
    let extent: Vec<f64> = (0..3).map(|k| hi[k] - lo[k]).collect();
    if !extent.iter().all(|e| e.is_finite() && *e > 0.) { return Err("The model is unbounded or empty; clip or intersect it with a box".into()); }
    let budget = input.get("budget").and_then(Value::as_f64).unwrap_or(1.2e6).clamp(1e4, 2.5e6);
    let requested = input.get("cell").and_then(Value::as_f64).or(model.cell)
        .unwrap_or_else(|| extent.iter().cloned().fold(f64::INFINITY, f64::min) / 24.);
    if !(requested.is_finite() && requested > 0.) { return Err("Mesh cell must be positive".into()); }
    let nodes = |c: f64| extent.iter().map(|e| (e / c).ceil() + 5.).product::<f64>();
    // The browser has a fixed budget: coarsen rather than run out of memory,
    // and say so. The analytic queries and the native CLI keep full detail.
    let mut cell = requested;
    while nodes(cell) > budget { cell *= 1.1; }
    let mesh = cad::mesh_solid(&solid, cell)?;
    let r = mesh.report();
    let binary = input.get("binary").and_then(Value::as_bool).unwrap_or(false);
    let mut out = json!({"name":model.name,"cell":cell,"requested_cell":requested,"coarsened":cell > requested * 1.000001,
        "watertight":r.watertight,"volume":r.volume,"exact":solid.exact(),"bbox":[lo,hi],
        "params":params.values,"triangle_count":mesh.triangles.len(),"vertex_count":mesh.vertices.len()});
    if binary {
        MESH.with(|b| {
            let mut b = b.borrow_mut();
            b.0 = mesh.vertices.iter().flat_map(|v| v.iter().map(|&x| x as f32)).collect();
            b.1 = mesh.triangles.iter().flat_map(|t| t.iter().map(|&i| i as u32)).collect();
        });
    } else {
        out["vertices"] = json!(mesh.vertices);
        out["triangles"] = json!(mesh.triangles);
    }
    Ok(out)
}

#[no_mangle]
pub extern "C" fn mesh_vertices_ptr() -> *const f32 { MESH.with(|b| b.borrow().0.as_ptr()) }
#[no_mangle]
pub extern "C" fn mesh_vertices_len() -> usize { MESH.with(|b| b.borrow().0.len()) }
#[no_mangle]
pub extern "C" fn mesh_indices_ptr() -> *const u32 { MESH.with(|b| b.borrow().1.as_ptr()) }
#[no_mangle]
pub extern "C" fn mesh_indices_len() -> usize { MESH.with(|b| b.borrow().1.len()) }

// A drawn board in SI units, as cad::kits::pcb::Board. Resistances come from
// the drawn length and copper section; fabrication files are written to the
// Gerber and Excellon specifications but are not verified against a CAM tool.
fn pcb_config(input: &Value) -> Result<Value, String> {
    use cad::kits::pcb::{excellon, gerber_layer, Board, Feature};
    let board: Board = serde_json::from_value(input.get("board").ok_or("No board")?.clone()).map_err(|e| format!("board JSON: {e}"))?;
    if board.features.len() > 512 || board.vias.len() > 512 || board.holes.len() > 512 { return Err("Board exceeds the browser feature budget".into()); }
    let temperature = setting(input, "temperature", 20., -55., 200.)?;
    let currents: Vec<f64> = input.get("currents").and_then(Value::as_array).map(|a| a.iter().map(|v| v.as_f64().unwrap_or(0.)).collect()).unwrap_or_default();
    let mut traces = Vec::new();
    let mut loss = 0.;
    for (i, f) in board.features.iter().enumerate() {
        if let Feature::Trace { layer, width, .. } = f {
            let r = board.trace_resistance(i, temperature)?;
            let current = currents.get(i).copied().unwrap_or(0.);
            if !current.is_finite() || current.abs() > 1000. { return Err("Trace current outside the supported range".into()); }
            loss += current * current * r;
            traces.push(json!({"feature":i,"layer":layer,"width":width,"length":f.length(),"resistance":r,"current":current,"loss":current*current*r}));
        }
    }
    let coppers = board.copper_layers();
    let sheet: Vec<f64> = (0..coppers.len()).map(|l| board.sheet_resistance(l)).collect::<Result<_, _>>()?;
    let mut out = json!({"traces":traces,"copper_loss":loss,"sheet_resistance":sheet,"thickness":board.thickness(),
        "thermal_resistance":board.thermal_resistance()?,"temperature":temperature});
    if input.get("fabrication").and_then(Value::as_bool).unwrap_or(false) {
        let cell = setting(input, "cell", 5e-5, 1e-5, 1e-3)?;
        let names = ["bottom", "inner", "top"];
        let layers: Vec<Value> = (0..coppers.len()).map(|l| {
            let name = if l == 0 { names[0] } else if l + 1 == coppers.len() { names[2] } else { names[1] };
            json!({"layer":l,"name":name,"gerber":gerber_layer(&board.copper(l), cell, &format!("{name} copper"))})
        }).collect();
        out["gerber"] = json!(layers);
        out["profile"] = json!(gerber_layer(&board.outline, cell, "board profile"));
        out["excellon"] = json!(excellon(&board));
    }
    Ok(out)
}

// Microstrip: Hammerstad-Jensen closed form, and optionally the quasi-TEM
// field solve of the same drawn section including copper thickness.
fn line_config(input: &Value) -> Result<Value, String> {
    use cad::laplace::{quasi_tem, Edge, Section};
    use cad::kits::pcb::{microstrip_width_for, microstrip_z0};
    let w = setting(input, "width", 3e-4, 1e-5, 2e-2)?;
    let h = setting(input, "height", 1.6e-3, 5e-5, 1e-2)?;
    let t = setting(input, "thickness", 3.5e-5, 0., 5e-4)?;
    let eps_r = setting(input, "eps_r", 4.4, 1., 20.)?;
    let (z0, eps_eff) = microstrip_z0(w, h, eps_r);
    let mut out = json!({"width":w,"height":h,"thickness":t,"eps_r":eps_r,"formula":{"z0":z0,"eps_eff":eps_eff,"source":"Hammerstad 1975; Hammerstad & Jensen 1980","valid":(0.05..=20.).contains(&(w/h))}});
    if let Some(target) = input.get("target").and_then(Value::as_f64) {
        if !(5. ..=200.).contains(&target) { return Err("Target impedance outside 5-200 ohm".into()); }
        out["width_for_target"] = json!(microstrip_width_for(target, h, eps_r));
    }
    if input.get("field").and_then(Value::as_bool).unwrap_or(false) {
        let span = 6. * (w + 2. * h);
        let height = 9. * h;
        let cell = (h / 16.).min(w / 6.).max(h / 48.);
        let section = Section::new(([0., 0.], [0.5 * span, height]), cell)
            .with_dielectric(cad::Region2::rect([0.25 * span, 0.5 * h], [0.5 * span, 0.5 * h]), eps_r)
            .with_conductor(cad::Region2::rect([0.25 * span, 0.], [0.5 * span, 1e-15]), 0.)
            .with_conductor(cad::Region2::rect([0., h + 0.5 * t.max(1e-9)], [0.5 * w, 0.5 * t.max(1e-9)]), 1.)
            .with_edges([Edge::Mirror, Edge::Grounded, Edge::Grounded, Edge::Grounded])
            .with_multiplicity(2.);
        let (zf, ef, sol) = quasi_tem(&section, 1., 1e-9, 40000);
        out["field"] = json!({"z0":zf,"eps_eff":ef,"cell":cell,"residual":sol.residual,"source":"cad::laplace quasi-TEM, two solves"});
    }
    Ok(out)
}

// JSON has no NaN: a table entry whose characterisation failed natively is
// written as null. Carry it through as NaN, exactly as the native flow holds
// it, rather than inventing a value. The sentinel is not a physical slew.
const NAN_SENTINEL: f64 = -1.0e300;
fn nan_in(v: Value) -> Value {
    match v {
        Value::Null => json!(NAN_SENTINEL),
        Value::Array(a) => Value::Array(a.into_iter().map(nan_in).collect()),
        Value::Object(o) => Value::Object(o.into_iter().map(|(k, x)| (k, nan_in(x))).collect()),
        x => x,
    }
}
fn restore_nan(lib: &mut semiconductor::timing::Library) {
    let fix = |x: &mut f64| if *x == NAN_SENTINEL { *x = f64::NAN };
    for c in &mut lib.cells { for o in &mut c.outs { for p in o {
        for t in [&mut p.rise, &mut p.fall, &mut p.e_rise, &mut p.e_fall] {
            t.slew.iter_mut().chain(t.load.iter_mut()).for_each(fix);
            t.delay.iter_mut().chain(t.out_slew.iter_mut()).flatten().for_each(fix);
            fix(&mut t.cert_error); fix(&mut t.cert_slew);
        }
    } } }
}

struct ChipCache { key: u64, netlist: chip::synth::Netlist, equivalence: Value, toggles: Vec<(u64, Vec<f64>)> }
fn fnv(bytes: &[u8]) -> u64 { bytes.iter().fold(0xcbf29ce484222325u64, |h, b| (h ^ *b as u64).wrapping_mul(0x100000001b3)) }

// Verilog through the native chip flow on cells characterised from the
// transistor law: synthesis, an equivalence check against the RTL, then
// timing and power interpolated across the nine (supply, temperature)
// corners. Corner libraries are supplied by the caller (see export test).
fn chip_config(input: &Value) -> Result<Value, String> {
    use chip::{corners, equiv, library, power, rtl, sta, synth};
    if let Some(libs) = input.get("libs").and_then(Value::as_array) {
        for l in libs {
            let corner = library::Corner { vdd: l["vdd"].as_f64().ok_or("corner supply")?, temp: l["temp"].as_f64().ok_or("corner temperature")? };
            let mut comb: semiconductor::timing::Library = serde_json::from_value(nan_in(l["comb"].clone())).map_err(|e| format!("cell library: {e}"))?;
            restore_nan(&mut comb);
            let dff = serde_json::from_value(l["dff"].clone()).map_err(|e| format!("flip-flop library: {e}"))?;
            library::preload(library::CellLib { corner, comb, dff });
        }
    }
    let source = input.get("verilog").and_then(Value::as_str).ok_or("No Verilog source")?;
    if source.len() > 65536 { return Err("Verilog source exceeds 64 KB".into()); }
    let samples = setting(input, "samples", 4000., 100., 50000.)? as u64;
    let key = fnv(source.as_bytes()) ^ samples.wrapping_mul(0x9e3779b97f4a7c15);
    // A few designs (CPU, PWM, a custom block) share the cache, most recent last.
    let fresh = CHIP.with(|c| { let mut c = c.borrow_mut(); match c.iter().position(|x| x.key == key) { Some(i) => { let hit = c.remove(i); c.push(hit); false } None => true } });
    if fresh {
        let module = rtl::parse(source)?;
        let netlist = synth::synthesize(&module)?;
        let eq = equiv::check(&module, &netlist, samples)?;
        let equivalence = json!({"points":eq.points,"exhaustive":eq.exhaustive,"mismatch":eq.mismatch});
        CHIP.with(|c| { let mut c = c.borrow_mut(); if c.len() >= 4 { c.remove(0); } c.push(ChipCache { key, netlist, equivalence, toggles: Vec::new() }); });
    }
    let vdd = setting(input, "vdd", 1.2, 1.08, 1.32)?;
    let temp = setting(input, "temperature", 300., 250., 425.)?;
    let freq = setting(input, "frequency", 100e6, 1e3, 1e10)?;
    let activity = setting(input, "activity", 0.2, 0., 1.)?;
    let copies = setting(input, "instances", 1., 1., 1e6)?;
    let wire = setting(input, "wire_cap", 0.5e-15, 0., 1e-12)?;
    CHIP.with(|cache| -> Result<Value, String> {
    let mut cache = cache.borrow_mut();
    let cache = cache.last_mut().unwrap();
    if !cache.toggles.iter().any(|(a, _)| *a == activity.to_bits()) {
        let t = power::activity(&cache.netlist, activity, 2000, 1);
        cache.toggles.push((activity.to_bits(), t));
    }
    let toggles = &cache.toggles.iter().find(|(a, _)| *a == activity.to_bits()).unwrap().1;
    let netlist = &cache.netlist;
    let c = sta::Constraints { period: 1. / freq, wire_cap: wire, ..Default::default() };
    let mut out = json!({"module":netlist.name,"cells":netlist.cells(),"transistors":netlist.transistors(),"flops":netlist.flops.len(),
        "inputs":netlist.inputs.iter().map(|(n,b)|json!([n,b.len()])).collect::<Vec<_>>(),"outputs":netlist.outputs.iter().map(|(n,b)|json!([n,b.len()])).collect::<Vec<_>>(),
        "equivalence":cache.equivalence.clone(),"vdd":vdd,"temperature":temp,"frequency":freq,"instances":copies,"cached":!fresh});
    if !cache.equivalence["mismatch"].is_null() { return Ok(out); }
    let (m, ind) = corners::at(netlist, &c, toggles, vdd, temp)?;
    // As in project::chip_design: the interpolation indicator is doubled.
    let period_err = m.min_period_error + 2. * ind.min_period;
    let switching = m.energy_per_cycle * freq * copies;
    let leakage = m.leakage * copies;
    out["timing"] = json!({"min_period":m.min_period,"min_period_error":period_err,"fmax":1./m.min_period,"setup_slack":c.period - m.min_period,"hold_slack":m.hold_slack});
    out["power"] = json!({"switching":switching,"leakage":leakage,"total":switching+leakage,"energy_per_cycle":m.energy_per_cycle,
        "switching_error":switching*0.15 + ind.energy_per_cycle*freq*copies,"leakage_error":leakage*(0.12+ind.leakage)});
    let near = |grid: &[f64; 3], x: f64| *grid.iter().min_by(|a, b| (*a - x).abs().total_cmp(&(*b - x).abs())).unwrap();
    let corner = library::Corner { vdd: near(&corners::VDDS, vdd), temp: near(&corners::TEMPS, temp) };
    if let Ok((_, t)) = corners::exact(netlist, &c, toggles, corner) {
        out["critical"] = json!({"vdd":corner.vdd,"temp":corner.temp,"endpoint":t.endpoint,
            "stages":t.critical.iter().map(|s|json!({"cell":s.cell,"edge":s.edge,"arrival":s.arrival})).collect::<Vec<_>>()});
    }
    Ok(out)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn small_pool(f: impl FnOnce() + Send) {
        rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap()
            .install(f);
    }
    #[test]
    fn browser_flow_matches_exact_decay_and_stays_solenoidal() {
        small_pool(|| {
            let mut f = Flow::new(16, true, 40.0);
            for _ in 0..10 {
                f.advance();
            }
            let d = f.report();
            let error = d["relative_error"].as_f64().unwrap();
            println!("Planar Taylor–Green energy relative error: {error:e}");
            assert!(error < 1e-12);
            assert!(d["divergence"].as_f64().unwrap() < 1e-14);
            f.frame(0);
            assert!(f.field.iter().all(|v| v.is_finite()));
            // Analytic curl_z = 2 sin(x) sin(y) exp(-2 nu t).
            let n = f.solver.grid.n;
            let exact = 2.0 * (-2.0 * f.solver.nu * f.solver.t).exp();
            let v = f.field[4 * ((n / 4) * n + n / 4) + 2] as f64;
            assert!((v - exact).abs() < 1e-6); // f32 display buffer, not solver state.
        });
    }
    #[test]
    fn browser_inputs_reject_invalid_values() {
        assert_eq!(flow_init(1024, 0, 100.0), 0);
        assert_eq!(flow_init(16, 0, f64::NAN), 0);
        assert!(design(0, 2.5).is_err());
        assert!(workspace_design(9, 2.5, f64::NAN, 45., 0.).is_err());
        assert!(workspace_design(9, 2.5, 30., 45., -1.).is_err());
    }
    #[test]
    fn shared_load_reaches_thermal_requirements() {
        small_pool(|| {
            let base = workspace_design(9, 2.5, 30., 45., 0.).unwrap();
            let loaded = workspace_design(9, 2.5, 30., 45., 5.).unwrap();
            let temperature = |v: &Value| v["checks"].as_array().unwrap().iter()
                .find(|c| c["id"] == "REQ-TJ").unwrap()["value"]["nominal"].as_f64().unwrap();
            assert!(temperature(&loaded) > temperature(&base));
            let mesh = sink_mesh(9).unwrap();
            assert_eq!(mesh["watertight"], true);
            // Exact volume: base plus 9 rectangular fins. Surface-net
            // rounding at 0.6 mm is bounded here by 5% for this preset.
            let exact = 40. * 50. * 3. + 9. * 1.2 * 50. * 20.;
            assert!((mesh["volume_mm3"].as_f64().unwrap() / exact - 1.).abs() < 0.05);
        });
    }
    #[test]
    fn geometry_and_absorbed_power_reach_the_shared_graph() {
        small_pool(|| {
            let base=workspace_config(&json!({})).unwrap();
            let longer=workspace_config(&json!({"length":70.})).unwrap();
            let temp=|v:&Value|v["nodes"]["q1.tj"]["nominal"].as_f64().unwrap();
            assert!(temp(&longer)<temp(&base));
            assert_eq!(longer["handoffs"]["sink_volume_mm3"],json!(40.*70.*3.+9.*1.2*70.*20.));
            let optical=workspace_config(&json!({"opticalPower":4.,"absorption":25.})).unwrap();
            let thermal=workspace_config(&json!({"heat":1.})).unwrap();
            assert_eq!(optical["handoffs"]["optical_heat_w"],json!(1.));
            assert!((temp(&optical)-temp(&thermal)).abs()<1e-12);
            assert!(workspace_config(&json!({"width":f64::INFINITY})).is_err());
            assert!(workspace_config(&json!({"fins":9.5})).is_err());
            let part=json!({"target":"spacer","features":[{"id":"spacer","shape":"cylinder","width":12.,"depth":12.,"height":10.,"x":0.,"y":0.,"z":0.}]});
            let mesh=mesh_config(&part).unwrap();
            assert_eq!(mesh["watertight"],true);
            let exact=std::f64::consts::PI*36.*10.;
            assert!((mesh["volume_mm3"].as_f64().unwrap()/exact-1.).abs()<0.03);
        });
    }
    // Characterises the nine STA corners natively and writes them for the
    // browser: cargo test --release -p emergence-web export_chip_corners -- --ignored
    #[test]
    #[ignore]
    fn export_chip_corners() {
        let mut libs = Vec::new();
        for vdd in chip::corners::VDDS { for temp in chip::corners::TEMPS {
            let l = chip::library::at(chip::library::Corner { vdd, temp }).unwrap();
            let comb = json!(l.comb);
            let failed = comb.to_string().matches("null").count();
            if failed > 0 { eprintln!("{vdd} V, {temp} K: {failed} table entries failed characterisation (NaN)"); }
            libs.push(json!({"vdd":vdd,"temp":temp,"comb":comb,"dff":l.dff,"failed_entries":failed}));
        } }
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../../web/cad/data/chip-corners.json");
        std::fs::create_dir_all(std::path::Path::new(path).parent().unwrap()).unwrap();
        std::fs::write(path, serde_json::to_string(&json!({"format":"emergence-chip-corners/1","technology":"150 nm EKV cells (analysis/gates.md)","corners":libs})).unwrap()).unwrap();
    }
    #[test]
    fn browser_chip_flow_matches_the_native_corner_flow() {
        small_pool(|| {
            let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../../web/cad/data/chip-corners.json");
            let corners: Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
            let verilog = include_str!("../../../projects/mini-cpu/mini8.v");
            let r = chip_config(&json!({"libs":corners["corners"],"verilog":verilog,"frequency":100e6,"temperature":337.5,"samples":2000})).unwrap();
            assert!(r["equivalence"]["mismatch"].is_null());
            assert!(r["transistors"].as_u64().unwrap() > 100);
            let module = chip::rtl::parse(verilog).unwrap();
            let netlist = chip::synth::synthesize(&module).unwrap();
            let toggles = chip::power::activity(&netlist, 0.2, 2000, 1);
            let c = chip::sta::Constraints { period: 1e-8, ..Default::default() };
            let (m, _) = chip::corners::at(&netlist, &c, &toggles, 1.2, 337.5).unwrap();
            assert!((r["timing"]["min_period"].as_f64().unwrap() / m.min_period - 1.).abs() < 1e-12);
            assert!(r["critical"]["stages"].as_array().unwrap().len() > 1);
        });
    }
    #[test]
    fn feature_tree_models_mesh_and_refuse_files() {
        small_pool(|| {
            let gear: Value = serde_json::from_str(include_str!("../../../models/gear.json")).unwrap();
            let r = model_config(&json!({"model":gear,"budget":3e5,"binary":true})).unwrap();
            assert_eq!(r["watertight"], true);
            assert!(MESH.with(|b| b.borrow().1.len()) == 3 * r["triangle_count"].as_u64().unwrap() as usize);
            let bad = json!({"model":{"shape":{"op":"mesh","file":"/etc/passwd"}}});
            assert!(model_config(&bad).unwrap_err().contains("browser"));
            let cube = json!({"model":{"params":{"a":2},"shape":{"op":"box","center":[0,0,0],"half":["a/2","a/2","a/2"]}},"overrides":{"a":4},"cell":0.25});
            let c = model_config(&cube).unwrap();
            assert!((c["volume"].as_f64().unwrap() / 64. - 1.).abs() < 0.02);
            assert_eq!(c["params"]["a"], json!(4.));
        });
    }
    #[test]
    fn board_resistance_and_line_impedance_use_the_kit() {
        let board = json!({"outline":{"op":"rect","center":[0.06,0.01],"half":[0.06,0.01]},
            "stackup":[{"kind":"copper","ounces":1.0},{"kind":"dielectric","thickness":0.0016},{"kind":"copper","ounces":1.0}],
            "features":[{"kind":"trace","layer":1,"points":[[0.005,0.01],[0.105,0.01]],"width":0.001}]});
        let r = pcb_config(&json!({"board":board,"currents":[2.0]})).unwrap();
        // rho L / (w t) for 100 mm of 1 mm, 1 oz copper at 20 C.
        let exact = 1.678e-8 * 0.1 / (0.001 * 34.8e-6);
        assert!((r["traces"][0]["resistance"].as_f64().unwrap() / exact - 1.).abs() < 1e-12);
        assert!((r["copper_loss"].as_f64().unwrap() - 4. * exact).abs() < 1e-12);
        let fab = pcb_config(&json!({"board":board,"fabrication":true})).unwrap();
        assert!(fab["gerber"][1]["gerber"].as_str().unwrap().contains("G36"));
        assert!(fab["excellon"].as_str().unwrap().starts_with("M48"));
        let line = line_config(&json!({"width":0.0029,"height":0.0016,"thickness":35e-6,"eps_r":4.4,"field":true,"target":50.})).unwrap();
        let (zf, z) = (line["field"]["z0"].as_f64().unwrap(), line["formula"]["z0"].as_f64().unwrap());
        assert!((zf / z - 1.).abs() < 0.05, "field {zf} formula {z}");
        assert!((line["width_for_target"].as_f64().unwrap() - 0.00306).abs() < 2e-4);
    }
    #[test]
    fn optical_selector_unfolds_and_reports_unreachable_tolerance() {
        small_pool(|| {
            let l = mie_ladder::ladder(vec![], 1);
            let small = optics_record(&l, 0, 0.02, 1e-3);
            let large = optics_record(&l, 0, 4.0, 1e-3);
            assert_eq!(small["rung"], 0);
            assert_eq!(large["rung"], 1);
            assert!(small["error"].as_f64().unwrap() <= small["claim"].as_f64().unwrap() + 2e-8);
            assert_eq!(optics_record(&l, 0, 4.0, 1e-10)["met"], false);
        });
    }
    #[test]
    fn design_preserves_undecided_and_out_of_regime_results() {
        small_pool(|| {
            let base = design(9, 2.5).unwrap();
            let checks = base["checks"].as_array().unwrap();
            assert_eq!(
                checks.iter().find(|c| c["id"] == "REQ-TJ").unwrap()["status"],
                "met"
            );
            assert_eq!(
                checks.iter().find(|c| c["id"] == "REQ-DROP").unwrap()["status"],
                "undecided"
            );
            let outside = design(6, 3.5).unwrap();
            assert_eq!(
                outside["checks"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .find(|c| c["id"] == "REQ-TJ")
                    .unwrap()["status"],
                "unavailable"
            );
        });
    }
}
