//! Open-boundary reflection: an acoustic pulse leaves the domain; what
//! remains afterwards is the reflection.
//!
//!   boundary_validate
//!
//! 1D simple wave in air (normal incidence), 2D cylindrical pulse (all
//! incidence angles), and a water pulse in the five-equation model.
use compressible::eos::materials::{AIR, WATER_SG};
use compressible::eos::Nasg;
use compressible::multiphase::FiveEq;
use compressible::{Bc, Euler, Mesh, Solver, System};

fn main() {
    let eos = Nasg::ideal(1.4);
    let (rho0, p0) = (1.0f64, 1.0 / 1.4);
    let c0 = (1.4 * p0 / rho0).sqrt();
    let eps = 1e-3;
    let sys = Euler { eos };
    let far = sys.to_cons(&[rho0, 0.0, 0.0, 0.0, p0]);
    println!("1D right-going pulse, relative amplitude {eps}: residual |p − p0|/(ε p0) after it has left");
    for (name, bc) in [("zero-gradient", Bc::Outflow), ("characteristic", Bc::Characteristic(far))] {
        let n = 400;
        let mesh = Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]);
        let mut s = Solver::new(sys, mesh, |x| {
            let dp = eps * p0 * (-((x[0] - 0.5) / 0.05).powi(2)).exp();
            [rho0 + dp / (c0 * c0), dp / (rho0 * c0), 0.0, 0.0, p0 + dp]
        });
        s.bc[0] = [bc; 2];
        s.cfl = 0.4;
        s.run_to(1.0 / c0);
        let res = (0..n).map(|i| (s.prim(i, 0, 0)[4] - p0).abs()).fold(0.0, f64::max) / (eps * p0);
        println!("  {name:15}: {res:.2e}");
    }
    println!("2D cylindrical pulse (all incidence angles): residual after the wave has left");
    for (name, bc) in [("zero-gradient", Bc::Outflow), ("characteristic", Bc::Characteristic(far))] {
        let n = 160;
        let mesh = Mesh::new([n, n, 1], [-1.0, -1.0, 0.0], [1.0, 1.0, 1.0]);
        let mut s = Solver::new(sys, mesh, |x| {
            let r2 = x[0] * x[0] + x[1] * x[1];
            let dp = eps * p0 * (-r2 / 0.01).exp();
            [rho0 + dp / (c0 * c0), 0.0, 0.0, 0.0, p0 + dp]
        });
        s.bc[0] = [bc; 2];
        s.bc[1] = [bc; 2];
        s.cfl = 0.4;
        // peak residual after the main front has passed the corners
        s.run_to(1.0 / c0);
        let mut mx: f64 = 0.0;
        let t1 = 3.0 / c0;
        let mut worst_late: f64 = 0.0;
        while s.t < t1 {
            let dt = s.dt_cfl().min(t1 - s.t);
            s.step(dt);
            if s.t > 2.0 / c0 {
                let m = (0..n).flat_map(|i| (0..n).map(move |j| (i, j))).map(|(i, j)| (s.prim(i, j, 0)[4] - p0).abs()).fold(0.0, f64::max);
                worst_late = worst_late.max(m);
            }
            mx = mx.max(0.0);
        }
        let _ = mx;
        println!("  {name:15}: max residual for t > 2L/c: {:.2e} of ε p0", worst_late / (eps * p0));
    }
    println!("water pulse (five-equation model, stiffened gas): residual after it has left");
    let sys5 = FiveEq { eos: [WATER_SG, AIR], kapila: true };
    let (rw, pw) = (1000.0, 1e5);
    let a = 1.0 - 1e-8;
    let w0 = [a * rw, (1.0 - a) * 1.2, 0.0, 0.0, 0.0, pw, a];
    let cw = sys5.sound_speed(&w0);
    let far5 = sys5.to_cons(&w0);
    for (name, bc) in [("zero-gradient", Bc::Outflow), ("characteristic", Bc::Characteristic(far5))] {
        let n = 400;
        let mesh = Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]);
        let amp = 1e5;
        let mut s = Solver::new(sys5, mesh, |x| {
            let dp = amp * (-((x[0] - 0.5) / 0.05).powi(2)).exp();
            let mut w = w0;
            w[0] += a * dp / (cw * cw);
            w[2] = dp / (rw * cw);
            w[5] += dp;
            w
        });
        s.bc[0] = [bc; 2];
        s.cfl = 0.4;
        s.run_to(1.0 / cw);
        let res = (0..n).map(|i| (s.prim(i, 0, 0)[5] - pw).abs()).fold(0.0, f64::max) / amp;
        println!("  {name:15}: {res:.2e}");
    }
}
