//! Reproducible conduction-heavy fluid stepping with fixed timesteps.
//! cpu_thermal [gas|two|three] [N=32] [steps=20]
use compressible::{
    eos::Nasg, multiphase::FiveEq, thermal::Conduction, threephase::Kapila3, Euler, FluxMode, Mesh,
    Solver, System,
};
use std::{f64::consts::TAU, time::Instant};

fn probe<const N: usize, S: System<N>>(mut s: Solver<N, S>, k: Vec<f64>, case: &str, steps: usize) {
    s.mode = FluxMode::Central;
    s.conduction = Some(Conduction::new(k));
    let before = s.totals();
    let start = Instant::now();
    for _ in 0..steps {
        s.step(1e-4);
    }
    let elapsed = start.elapsed().as_secs_f64();
    let after = s.totals();
    let energy = s.sys.energy_index().unwrap();
    let energy_error = (after[energy] / before[energy] - 1.0).abs();
    assert!(energy_error < 1e-11);
    let mut hash = 0xcbf29ce484222325_u64;
    let mut min_temperature = f64::INFINITY;
    let mut max_temperature = 0.0_f64;
    for i in 0..s.mesh.n[0] {
        for j in 0..s.mesh.n[1] {
            for k in 0..s.mesh.n[2] {
                let q = s
                    .mesh
                    .idx(i + s.mesh.g(0), j + s.mesh.g(1), k + s.mesh.g(2));
                let w = s.sys.to_prim(&s.u[q]);
                assert!(s.sys.admissible(&w));
                let temperature = s.sys.temperature(&w);
                assert!(temperature.is_finite() && temperature > 0.0);
                min_temperature = min_temperature.min(temperature);
                max_temperature = max_temperature.max(temperature);
                for value in s.u[q] {
                    hash = (hash ^ value.to_bits()).wrapping_mul(0x100000001b3);
                }
            }
        }
    }
    println!("{case} n={} steps={steps} wall_s={elapsed:.9} substeps={} energy_error={energy_error:.9e} min_t={min_temperature:.9e} max_t={max_temperature:.9e} hash={hash:016x}", s.mesh.n[0], s.conduction_substeps);
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let case = args.get(1).map(String::as_str).unwrap_or("gas");
    let n = args
        .get(2)
        .map(|s| s.parse::<usize>().unwrap())
        .unwrap_or(32);
    let steps = args
        .get(3)
        .map(|s| s.parse::<usize>().unwrap())
        .unwrap_or(20);
    assert!(n > 1 && steps > 0);
    let mesh = Mesh::new([n; 3], [0.0; 3], [1.0; 3]);
    let state = |x: [f64; 3]| {
        let rho = 1.0 + 0.02 * (TAU * x[0]).sin();
        let temperature = 2.0 + 0.1 * (TAU * x[1]).cos() * (TAU * x[2]).sin();
        (rho, 0.4 * rho * temperature)
    };
    let eos = Nasg::ideal(1.4);
    match case {
        "gas" => probe(
            Solver::new(Euler { eos }, mesh, |x| {
                let (rho, p) = state(x);
                [rho, 0.01, 0.02, 0.03, p]
            }),
            vec![20.0],
            case,
            steps,
        ),
        "two" => probe(
            Solver::new(
                FiveEq {
                    eos: [eos; 2],
                    kapila: true,
                },
                mesh,
                |x| {
                    let (rho, p) = state(x);
                    let alpha = 0.4 + 0.05 * (TAU * x[0]).sin();
                    [alpha * rho, (1.0 - alpha) * rho, 0.01, 0.02, 0.03, p, alpha]
                },
            ),
            vec![10.0, 30.0],
            case,
            steps,
        ),
        "three" => probe(
            Solver::new(Kapila3 { eos: [eos; 3] }, mesh, |x| {
                let (rho, p) = state(x);
                let a0 = 0.2 + 0.02 * (TAU * x[0]).sin();
                let a1 = 0.3 + 0.03 * (TAU * x[1]).cos();
                [
                    a0 * rho,
                    a1 * rho,
                    (1.0 - a0 - a1) * rho,
                    0.01,
                    0.02,
                    0.03,
                    p,
                    a0,
                    a1,
                ]
            }),
            vec![10.0, 20.0, 30.0],
            case,
            steps,
        ),
        _ => panic!("gas, two or three expected"),
    }
}
