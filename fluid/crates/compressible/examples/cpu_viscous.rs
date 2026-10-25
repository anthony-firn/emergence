//! Reproducible viscous RHS and full-step CPU probes with fixed timesteps.
//! cpu_viscous [rhs|step|line] [N=32] [iterations=40]
use compressible::{
    eos::Nasg,
    viscous::{add_viscous, Transport},
    Euler, FluxMode, Mesh, Solver, System,
};
use std::{f64::consts::TAU, time::Instant};

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let case = args.get(1).map(String::as_str).unwrap_or("rhs");
    assert!(["rhs", "step", "line"].contains(&case));
    let n = args
        .get(2)
        .map(|s| s.parse::<usize>().unwrap())
        .unwrap_or(32);
    let iterations = args
        .get(3)
        .map(|s| s.parse::<usize>().unwrap())
        .unwrap_or(40);
    assert!(n > 1 && iterations > 0);
    let mesh = Mesh::new(
        if case == "line" { [n, 1, 1] } else { [n; 3] },
        [0.0; 3],
        [1.0; 3],
    );
    let sys = Euler {
        eos: Nasg::ideal(1.4),
    };
    let initial = |x: [f64; 3]| {
        [
            1.0 + 0.05 * (TAU * x[0]).sin(),
            0.1 * (TAU * x[1]).sin(),
            0.1 * (TAU * x[0]).cos(),
            0.1 * (TAU * x[2]).sin(),
            1.0 + 0.05 * (TAU * x[0]).cos(),
        ]
    };
    let mut solver = Solver::new(sys, mesh.clone(), initial);
    solver.mode = FluxMode::Central;
    solver.viscous = Some(Transport {
        mu: 0.02,
        power: 0.7,
        t_ref: 2.5,
        prandtl: 0.75,
        bulk: 0.01,
    });
    let mut residual = vec![[0.0; 5]; mesh.len()];
    // Analytic periodic fields also fill corners and deep ghost layers.
    if case != "step" {
        for i in 0..mesh.m(0) {
            for j in 0..mesh.m(1) {
                for k in 0..mesh.m(2) {
                    let c = [i, j, k];
                    let x = std::array::from_fn(|d| {
                        let index =
                            (c[d] as isize - mesh.g(d) as isize).rem_euclid(mesh.n[d] as isize);
                        (index as f64 + 0.5) * mesh.dx[d]
                    });
                    solver.u[mesh.idx(i, j, k)] = sys.to_cons(&initial(x));
                }
            }
        }
    }
    let start = Instant::now();
    for _ in 0..iterations {
        if case == "step" {
            solver.step(1e-4);
        } else {
            add_viscous(
                &sys,
                &mesh,
                solver.viscous.as_ref().unwrap(),
                &solver.u,
                &mut residual,
            );
        }
    }
    let elapsed = start.elapsed().as_secs_f64();
    let data = if case == "step" { &solver.u } else { &residual };
    let mut hash = 0xcbf29ce484222325_u64;
    let mut sum = [0.0; 5];
    for i in 0..mesh.n[0] {
        for j in 0..mesh.n[1] {
            for k in 0..mesh.n[2] {
                let q = mesh.idx(i + mesh.g(0), j + mesh.g(1), k + mesh.g(2));
                if case == "step" {
                    assert!(sys.admissible(&sys.to_prim(&data[q])));
                }
                for (d, value) in data[q].iter().enumerate() {
                    assert!(value.is_finite());
                    sum[d] += value;
                    hash = (hash ^ value.to_bits()).wrapping_mul(0x100000001b3);
                }
            }
        }
    }
    let balance = if case == "step" {
        (sum[0] / mesh.cells() as f64 - 1.0).abs()
    } else {
        sum.into_iter().map(f64::abs).fold(0.0, f64::max) / mesh.cells() as f64
    };
    println!("{case} n={n} steps={iterations} wall_s={elapsed:.9} balance={balance:.9e} hash={hash:016x}");
}
