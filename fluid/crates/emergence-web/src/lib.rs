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
    if !(5..=14).contains(&fins) || !velocity.is_finite() || !(0.5..=4.0).contains(&velocity) {
        return Err("Design settings outside the demo range".into());
    }
    let mut doc: project::Doc =
        serde_json::from_str(include_str!("../../../projects/motor-driver/project.json"))
            .map_err(|e| e.to_string())?;
    // This browser view evaluates the original thermal/electrical graph.
    // CAD, printing and signal-integrity examples remain in the native repo.
    doc.parts.clear();
    doc.analyses.remove("enc");
    doc.requirements
        .retain(|id, _| ["REQ-TJ", "REQ-DROP", "REQ-FAN"].contains(&id.as_str()));
    doc.params.get_mut("sink.fin_count").unwrap().value = fins.to_string();
    doc.params.get_mut("duct.velocity").unwrap().value = format!("{velocity} m/s");
    let p = project::Project::new(doc, ".", Default::default())?;
    let e = p.evaluate();
    let checks: Vec<_> = e.checks.values().map(|c|json!({
        "id":c.id,"text":c.text,"status":c.status.name(),"reason":c.reason,
        "value":c.value.map(|v|json!({"nominal":v.nom,"lo":v.lo,"hi":v.hi,"rigor":format!("{:?}",v.rigor)})),
        "unit":c.unit,"owner":c.owner
    })).collect();
    Ok(json!({"checks":checks,"computed":e.computed,"cached":e.cached}))
}
#[no_mangle]
pub extern "C" fn design_solve(fins: u32, velocity: f64) {
    output(design(fins, velocity).unwrap_or_else(|e| json!({"error":e})));
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
