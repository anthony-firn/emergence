//! CPU timestep probes. cpu_timestep [mixture|allspeed|viscous|advancing] [N=32] [queries=200]
use compressible::{
    allspeed::AllSpeedOpts, eos::Nasg, multiphase::FiveEq, viscous::Transport, Euler, Mesh, Solver,
    System,
};
use std::{f64::consts::TAU, hint::black_box, time::Instant};

fn probe<const N: usize, S: System<N>>(mut s: Solver<N, S>, count: usize, advance: bool) {
    let mut h = 0xcbf29ce484222325_u64;
    let mut smallest = f64::INFINITY;
    let start = Instant::now();
    for _ in 0..count {
        let dt = black_box(&s).dt_cfl();
        assert!(dt.is_finite() && dt > 0.0);
        smallest = smallest.min(dt);
        h = (h ^ dt.to_bits()).wrapping_mul(0x100000001b3);
        if advance {
            s.step(dt.min(1e-5));
        }
    }
    let elapsed = start.elapsed().as_secs_f64();
    for i in 0..s.mesh.n[0] {
        for j in 0..s.mesh.n[1] {
            for k in 0..s.mesh.n[2] {
                let w = s.prim(i, j, k);
                assert!(s.sys.admissible(&w));
                let q = s
                    .mesh
                    .idx(i + s.mesh.g(0), j + s.mesh.g(1), k + s.mesh.g(2));
                for v in s.u[q] {
                    h = (h ^ v.to_bits()).wrapping_mul(0x100000001b3);
                }
            }
        }
    }
    println!("n={} queries={count} wall_s={elapsed:.9} min_dt={smallest:.17e} fluid_time={:.17e} hash={h:016x}",s.mesh.n[0],s.t);
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let case = args.get(1).map(String::as_str).unwrap_or("mixture");
    let n = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(32);
    let count = args.get(3).map(|s| s.parse().unwrap()).unwrap_or(200);
    assert!(n >= 8 && count > 0);
    let mesh = Mesh::new([n; 3], [0.0; 3], [1.0; 3]);
    if case == "viscous" {
        let mut s = Solver::new(
            Euler {
                eos: Nasg::ideal(1.4),
            },
            mesh,
            |x| {
                [
                    1.0 + 0.1 * (TAU * x[0]).sin() * (TAU * x[1]).cos(),
                    0.3,
                    0.2,
                    0.1,
                    1.0,
                ]
            },
        );
        s.viscous = Some(Transport {
            mu: 1.0,
            power: 0.0,
            t_ref: 1.0,
            prandtl: 0.7,
            bulk: 0.0,
        });
        probe(s, count, false);
    } else {
        assert!(["mixture", "allspeed", "advancing"].contains(&case));
        let sys = FiveEq {
            eos: [
                Nasg {
                    p_inf: 1e6,
                    cv: 1000.0,
                    ..Nasg::ideal(2.0)
                },
                Nasg::ideal(1.4),
            ],
            kapila: true,
        };
        let mut s = Solver::new(sys, mesh, |x| {
            let a = 0.5 + 0.3 * (TAU * x[0]).sin() * (TAU * x[1]).cos();
            [
                a * 1000.0,
                (1.0 - a) * 1.2,
                0.3 + 0.1 * (TAU * x[0]).sin(),
                0.2 * (TAU * x[1]).cos(),
                0.1 * (TAU * x[2]).sin(),
                1e5,
                a,
            ]
        });
        if case == "allspeed" {
            s.allspeed = Some(AllSpeedOpts::default());
            // A completed pressure solve provides this acceleration estimate.
            // The timestep still needs advective and multiphase source limits.
            s.allspeed_stats.last_accel = 1.0;
        }
        probe(s, count, case == "advancing");
    }
}
