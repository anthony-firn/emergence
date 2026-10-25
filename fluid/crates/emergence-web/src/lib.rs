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
    let current = setting(p,"current",30.,5.,40.)?;
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
    doc.params.get_mut("q1.switching_loss").unwrap().value = format!("{} W", 3.0 + heat + optical);
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
    if len == 0 || len > 65536 { return std::ptr::null_mut(); }
    INPUT.with(|b| { let mut b = b.borrow_mut(); b.resize(len as usize,0); b.as_mut_ptr() })
}

#[no_mangle]
pub extern "C" fn workspace_run(mesh: u32) {
    let result = INPUT.with(|b| {
        let input: Value = serde_json::from_slice(&b.borrow()).map_err(|e|e.to_string())?;
        if mesh == 0 { workspace_config(&input) } else { mesh_config(&input) }
    });
    output(result.unwrap_or_else(|e: String| json!({"error":e})));
}

#[no_mangle]
pub extern "C" fn cad_mesh(fins: u32) {
    output(sink_mesh(fins).unwrap_or_else(|e| json!({"error":e})));
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
