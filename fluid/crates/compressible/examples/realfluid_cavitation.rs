//! Cavitation with the reference water EOS (IAPWS-95 tabulated, HEM dome)
//! vs the NASG five-equation model, same expansion tube.
//!   realfluid_cavitation N [u] [T0]
use compressible::realflow::EulerReal;
use compressible::realfluid::RealFluidEos;
use compressible::reference_data::water_iapws95;
use compressible::{Bc, Mesh, Solver};
use std::sync::Arc;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let n: usize = a.get(1).map(|v| v.parse().unwrap()).unwrap_or(400);
    let u0: f64 = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(2.0);
    let t0: f64 = a.get(3).map(|v| v.parse().unwrap()).unwrap_or(355.0);
    let eos = Arc::new(RealFluidEos::build(water_iapws95(), 0.01, 1100.0, 275.0, 900.0, 200));
    let h = eos.h().clone();
    let p0 = 101_325.0;
    let rho0 = h.density(p0, t0, 970.0).unwrap();
    let e0 = h.props(rho0, t0).u;
    let ps = h.saturation(t0).unwrap().0;
    println!("liquid water {t0} K, 1 atm: ρ {rho0:.3} kg/m³, c {:.1} m/s; IAPWS-95 p_sat({t0} K) = {ps:.1} Pa", h.props(rho0, t0).w);
    let sys = EulerReal { eos: eos.clone() };
    let mesh = Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]);
    let mut s = Solver::new(sys.clone(), mesh, |x| [rho0, if x[0] < 0.5 { -u0 } else { u0 }, 0.0, 0.0, e0]);
    s.bc[0] = [Bc::Outflow; 2];
    s.cfl = 0.4;
    let t_end = 3.2e-3;
    let t1 = std::time::Instant::now();
    if std::env::var("DEBUG").is_ok() {
        for it in 0..400 {
            let dt = s.dt_cfl();
            let wc = s.prim(n / 2, 0, 0);
            if it > 25 { println!("   centre w {:?} state {:?} two {}", wc, sys.state(&wc), eos.two_phase_region(wc[0], wc[4])); }
            let (mut cm, mut im) = (0.0f64, 0);
            for i in 0..n { let c = sys.state(&s.prim(i, 0, 0)).c; if !(c <= cm) { cm = c; im = i; } }
            let w = s.prim(im, 0, 0);
            if it % 20 == 0 || it > 30 || !(dt > 1e-12) { println!("it {it} dt {dt:e} max c {cm:.3e} at {im}: w {:?} state {:?} two {}", w, sys.state(&w), eos.two_phase_region(w[0], w[4])); }
            if !(dt > 1e-12) { break; }
            s.step(dt);
        }
        return;
    }
    let steps = s.run_to(t_end);
    let w = s.prim(n / 2, 0, 0);
    let st = sys.state(&w);
    let tex = eos.exact(w[0], w[4]).unwrap();
    println!(
        "centre after {t_end} s ({steps} steps, {:.1}s): ρ {:.3}, p {:.1} Pa, T {:.3} K, vapour quality {:.3e}; IAPWS p_sat(T) = {:.1} Pa (ratio {:.6}); exact EOS at that state p {:.1} Pa",
        t1.elapsed().as_secs_f64(),
        w[0],
        st.p,
        st.t,
        st.x,
        h.saturation(st.t).unwrap().0,
        st.p / h.saturation(st.t).unwrap().0,
        tex.p
    );
    let vapour_cells = (0..n).filter(|&i| sys.state(&s.prim(i, 0, 0)).x > 0.0).count();
    println!("cells with vapour: {vapour_cells}; shock-path faces {:.1}%, limited faces {}", 100.0 * s.stats.weno_faces as f64 / s.stats.faces as f64, s.stats.limited_faces);
}
