//! Reproducible CPU timing plus accuracy/conservation, no rendering.
//! cpu_validate [smooth|weno|fdweno|interface|open] [N=32] [steps=40]
use compressible::{
    eos::Nasg, multiphase::FiveEq, scheme::ShockScheme, Bc, Euler, FluxMode, Mesh, Solver, System,
};
use std::{f64::consts::TAU, time::Instant};
fn hash<const N: usize, S: System<N>>(s: &Solver<N, S>) -> u64 {
    let mut h = 0xcbf29ce484222325_u64;
    for i in 0..s.mesh.n[0] {
        for j in 0..s.mesh.n[1] {
            for k in 0..s.mesh.n[2] {
                let q = s
                    .mesh
                    .idx(i + s.mesh.g(0), j + s.mesh.g(1), k + s.mesh.g(2));
                for v in s.u[q] {
                    h = (h ^ v.to_bits()).wrapping_mul(0x100000001b3);
                }
            }
        }
    }
    h
}
fn main() {
    let a: Vec<_> = std::env::args().collect();
    let case = a.get(1).map(String::as_str).unwrap_or("smooth");
    let n = a.get(2).map(|s| s.parse::<usize>().unwrap()).unwrap_or(32);
    let steps = a.get(3).map(|s| s.parse::<usize>().unwrap()).unwrap_or(40);
    let mesh = Mesh::new([n; 3], [0.0; 3], [1.0; 3]);
    if ["smooth", "weno", "fdweno"].contains(&case) {
        let exact = |x: [f64; 3], t: f64| {
            1.0 + 0.1 * (TAU * (x[0] - 0.3 * t)).sin() * (TAU * (x[1] - 0.2 * t)).cos()
        };
        let mut s = Solver::new(
            Euler {
                eos: Nasg::ideal(1.4),
            },
            mesh,
            |x| [exact(x, 0.0), 0.3, 0.2, 0.1, 1.0],
        );
        s.mode = if case == "smooth" {
            FluxMode::Central
        } else {
            FluxMode::Weno
        };
        if case == "fdweno" {
            s.shock = ShockScheme::FdWeno;
        }
        // Fixed final time, dt small enough to resolve sixth-order space error.
        let dt = 0.02 / steps as f64;
        let start = Instant::now();
        for _ in 0..steps {
            s.step(dt);
        }
        let elapsed = start.elapsed().as_secs_f64();
        let mut err = 0.0;
        let mut mass = 0.0;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let w = s.prim(i, j, k);
                    assert!(s.sys.admissible(&w));
                    err += (w[0] - exact(s.mesh.x(i, j, k), s.t)).abs();
                    mass += w[0];
                }
            }
        }
        println!("{case} n={n} steps={steps} wall_s={elapsed:.9} l1_rho={:.9e} mass_error={:.9e} hash={:016x}",err/(n*n*n) as f64,(mass/(n*n*n) as f64-1.0).abs(),hash(&s));
    } else if case == "open" {
        let sys = Euler {
            eos: Nasg::ideal(1.4),
        };
        let far = sys.to_cons(&[1.0, 0.0, 0.0, 0.0, 1.0]);
        let mut s = Solver::new(sys, mesh, |x| {
            let r2: f64 = x.into_iter().map(|v| (v - 0.5).powi(2)).sum();
            let dp = 1e-4 * (-r2 / 0.04).exp();
            [1.0 + dp / 1.4, 0.0, 0.0, 0.0, 1.0 + dp]
        });
        s.bc = [[Bc::Characteristic(far); 2]; 3];
        s.mode = FluxMode::Central;
        let dt = 0.02 / steps as f64;
        let start = Instant::now();
        for _ in 0..steps {
            s.step(dt);
        }
        let elapsed = start.elapsed().as_secs_f64();
        let mut pressure: f64 = 0.0;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let w = s.prim(i, j, k);
                    assert!(s.sys.admissible(&w));
                    pressure = pressure.max((w[4] - 1.0).abs());
                }
            }
        }
        println!(
            "open n={n} steps={steps} wall_s={elapsed:.9} max_dp={pressure:.9e} hash={:016x}",
            hash(&s)
        );
    } else if case == "interface" {
        let sys = FiveEq {
            eos: [
                Nasg {
                    gamma: 2.0,
                    p_inf: 1e6,
                    cv: 1000.0,
                    ..Nasg::ideal(2.0)
                },
                Nasg::ideal(1.4),
            ],
            kapila: true,
        };
        let mut s = Solver::new(sys, mesh, |x| {
            let r = ((x[0] - 0.5).powi(2) + (x[1] - 0.5).powi(2) + (x[2] - 0.5).powi(2)).sqrt();
            let alpha =
                (0.5 * (1.0 + ((0.25 - r) / (1.0 / n as f64)).tanh())).clamp(1e-6, 1.0 - 1e-6);
            [
                alpha * 1000.0,
                (1.0 - alpha) * 1.2,
                0.1,
                0.0,
                0.0,
                1e5,
                alpha,
            ]
        });
        s.mode = FluxMode::Hybrid;
        let dt = 1e-5;
        let start = Instant::now();
        for _ in 0..steps {
            s.step(dt);
        }
        let elapsed = start.elapsed().as_secs_f64();
        let mut dp: f64 = 0.0;
        let mut du: f64 = 0.0;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let w = s.prim(i, j, k);
                    assert!(s.sys.admissible(&w));
                    dp = dp.max((w[5] / 1e5 - 1.0).abs());
                    du = du.max((w[2] - 0.1).abs());
                }
            }
        }
        println!("interface n={n} steps={steps} wall_s={elapsed:.9} max_relative_dp={dp:.9e} max_du={du:.9e} hash={:016x}",hash(&s));
    } else {
        panic!("smooth, weno, fdweno, interface or open expected");
    }
}
