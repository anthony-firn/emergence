use compressible::{
    allspeed::AllSpeedOpts,
    eos::Nasg,
    viscous::{add_viscous, Transport},
    Euler, Mesh, Solver, System,
};
use std::f64::consts::TAU;

fn residual(
    n: [usize; 3],
    sys: Euler,
    tr: Transport,
    init: impl Fn([f64; 3]) -> [f64; 5],
) -> (Mesh, Vec<[f64; 5]>) {
    let m = Mesh::new(n, [0.0; 3], [1.0; 3]);
    let mut u = vec![[0.0; 5]; m.len()];
    for i in 0..m.m(0) {
        for j in 0..m.m(1) {
            for k in 0..m.m(2) {
                let c = [i, j, k];
                let x = std::array::from_fn(|d| {
                    let q = (c[d] as isize - m.g(d) as isize).rem_euclid(n[d] as isize);
                    (q as f64 + 0.5) * m.dx[d]
                });
                u[m.idx(i, j, k)] = sys.to_cons(&init(x));
            }
        }
    }
    // Poison ghosts to check that adding a residual only writes interiors.
    let mut out = vec![[f64::NAN; 5]; m.len()];
    for i in 0..n[0] {
        for j in 0..n[1] {
            for k in 0..n[2] {
                out[m.idx(i + m.g(0), j + m.g(1), k + m.g(2))] = [0.0; 5];
            }
        }
    }
    add_viscous(&sys, &m, &tr, &u, &mut out);
    for i in 0..m.m(0) {
        for j in 0..m.m(1) {
            for k in 0..m.m(2) {
                let c = [i, j, k];
                if (0..3).any(|d| c[d] < m.g(d) || c[d] >= m.g(d) + n[d]) {
                    assert!(out[m.idx(i, j, k)].iter().all(|v| v.is_nan()));
                }
            }
        }
    }
    (m, out)
}

#[test]
fn shear_stress_and_viscous_work_converge_at_sixth_order_on_every_axis() {
    let sys = Euler {
        eos: Nasg::ideal(1.4),
    };
    let (mu, amplitude) = (0.2, 0.1);
    for d in 0..3 {
        let transverse = (d + 1) % 3;
        let mut previous = None::<[f64; 2]>;
        for n in [16, 32, 64] {
            let mut shape = [1; 3];
            shape[d] = n;
            let (m, out) = residual(shape, sys, Transport::constant(mu, 0.75), |x| {
                let mut w = [1.0, 0.0, 0.0, 0.0, 1.0];
                w[1 + transverse] = amplitude * (TAU * x[d]).cos();
                w
            });
            let mut errors = [0.0; 2];
            let mut balance = [0.0; 5];
            for i in 0..n {
                let mut c = [0; 3];
                c[d] = i;
                let x = m.x(c[0], c[1], c[2])[d];
                let q = m.idx(c[0] + m.g(0), c[1] + m.g(1), c[2] + m.g(2));
                let exact = [
                    -mu * amplitude * TAU.powi(2) * (TAU * x).cos(),
                    -mu * amplitude.powi(2) * TAU.powi(2) * (2.0 * TAU * x).cos(),
                ];
                errors[0] += (out[q][1 + transverse] - exact[0]).abs() / n as f64;
                errors[1] += (out[q][4] - exact[1]).abs() / n as f64;
                for v in 0..5 {
                    balance[v] += out[q][v];
                }
                assert_eq!(out[q][0], 0.0);
            }
            assert!(balance.iter().all(|v| v.abs() < 1e-11));
            if let Some(p) = previous {
                for v in 0..2 {
                    assert!(
                        p[v] / errors[v] > 50.0,
                        "axis {d}, n {n}, errors {errors:?}, previous {p:?}"
                    );
                }
            }
            previous = Some(errors);
        }
    }
}

#[test]
fn longitudinal_bulk_stress_and_thermal_diffusion_have_the_physical_coefficients() {
    let sys = Euler {
        eos: Nasg {
            cv: 7.0,
            ..Nasg::ideal(2.3)
        },
    };
    let tr = Transport {
        bulk: 0.7,
        ..Transport::constant(0.2, 0.75)
    };
    let n = 96;
    for thermal in [false, true] {
        let (m, out) = residual([n, 1, 1], sys, tr, |x| {
            let wave = (TAU * x[0]).cos();
            let temperature = 3.0 + if thermal { 0.1 * wave } else { 0.0 };
            [
                1.0,
                if thermal { 0.0 } else { 0.1 * wave },
                0.0,
                0.0,
                (sys.eos.gamma - 1.0) * sys.eos.cv * temperature,
            ]
        });
        for i in 0..n {
            let x = m.x(i, 0, 0)[0];
            let q = m.idx(i + m.g(0), 0, 0);
            let coefficient = if thermal {
                tr.mu * sys.eos.gamma * sys.eos.cv / tr.prandtl
            } else {
                (4.0 / 3.0) * tr.mu + tr.bulk
            };
            let exact = -coefficient * 0.1 * TAU.powi(2) * (TAU * x).cos();
            assert!((out[q][if thermal { 4 } else { 1 }] - exact).abs() < 2e-7);
        }
    }
}

#[test]
fn temperature_dependent_shear_stress_converges_at_sixth_order() {
    let sys = Euler {
        eos: Nasg::ideal(1.4),
    };
    let tr = Transport {
        power: 1.0,
        t_ref: 3.0,
        ..Transport::constant(0.2, 0.75)
    };
    let mut previous = f64::INFINITY;
    for n in [24, 48, 96] {
        let (m, out) = residual([n, 1, 1], sys, tr, |x| {
            let temperature = 3.0 + 0.2 * (TAU * x[0]).sin();
            [
                1.0,
                0.0,
                0.1 * (TAU * x[0]).cos(),
                0.0,
                (sys.eos.gamma - 1.0) * temperature,
            ]
        });
        let error = (0..n)
            .map(|i| {
                let theta = TAU * m.x(i, 0, 0)[0];
                let exact = -tr.mu / tr.t_ref
                    * 0.1
                    * TAU.powi(2)
                    * (3.0 * theta.cos() + 0.4 * theta.sin() * theta.cos());
                (out[m.idx(i + m.g(0), 0, 0)][2] - exact).abs()
            })
            .sum::<f64>()
            / n as f64;
        assert!(previous / error > 50.0);
        previous = error;
    }
}

#[test]
fn viscous_timestep_uses_temperature_bulk_and_actual_heat_capacity_ratio() {
    for allspeed in [false, true] {
        for (gamma, cv, tr, temperature) in [
            (1.4, 1.0, Transport::constant(2.0, 0.75), 3.0),
            (
                1.4,
                1.0,
                Transport {
                    power: 1.0,
                    ..Transport::constant(2.0, 0.75)
                },
                50.0,
            ),
            (
                1.4,
                1.0,
                Transport {
                    bulk: 100.0,
                    ..Transport::constant(2.0, 1e30)
                },
                3.0,
            ),
            (2.3, 7.0, Transport::constant(2.0, 0.1), 3.0),
        ] {
            let rho = 2.0;
            let mut s = Solver::new(
                Euler {
                    eos: Nasg {
                        cv,
                        ..Nasg::ideal(gamma)
                    },
                },
                Mesh::new([16, 12, 1], [0.0; 3], [1.0; 3]),
                |_| [rho, 0.0, 0.0, 0.0, rho * (gamma - 1.0) * cv * temperature],
            );
            s.viscous = Some(tr);
            if allspeed {
                s.allspeed = Some(AllSpeedOpts::default());
            }
            let mu = tr.mu * (temperature / tr.t_ref).powf(tr.power);
            let diffusion = ((mu * (4.0 / 3.0) + tr.bulk) / rho)
                .max(mu * gamma / tr.prandtl * (4.0 / 3.0) / rho);
            let expected = 0.3 / (diffusion * (16.0_f64.powi(2) + 12.0_f64.powi(2)));
            assert!(
                (s.dt_cfl() / expected - 1.0).abs() < 1e-13,
                "allspeed {allspeed}, transport {tr:?}"
            );
            // This query must read current interiors, rather than unfilled ghosts.
            for value in &mut s.u {
                if value[0] == 0.0 {
                    *value = [f64::NAN; 5];
                }
            }
            assert!((s.dt_cfl() / expected - 1.0).abs() < 1e-13);
        }
    }
}

#[test]
fn viscous_timestep_finds_the_hottest_cell_even_when_it_is_denser() {
    for allspeed in [false, true] {
        let gamma = 1.4;
        let mut s = Solver::new(
            Euler {
                eos: Nasg::ideal(gamma),
            },
            Mesh::new([16, 1, 1], [0.0; 3], [1.0; 3]),
            |x| {
                let (rho, temperature) = if x[0] < 0.5 { (1.0, 1.0) } else { (2.0, 100.0) };
                [rho, 0.0, 0.0, 0.0, rho * (gamma - 1.0) * temperature]
            },
        );
        s.viscous = Some(Transport {
            power: 1.0,
            ..Transport::constant(2.0, 1e30)
        });
        if allspeed {
            s.allspeed = Some(AllSpeedOpts::default());
        }
        let expected = 0.3 / ((4.0 / 3.0) * 200.0 * 16.0_f64.powi(2));
        assert!((s.dt_cfl() / expected - 1.0).abs() < 1e-13);
    }
}

#[test]
fn cached_viscous_flux_is_reproducible_on_short_and_multidimensional_meshes() {
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let three = rayon::ThreadPoolBuilder::new()
        .num_threads(3)
        .build()
        .unwrap();
    for n in [[2, 1, 1], [1, 3, 2], [8, 9, 7], [1, 12, 16]] {
        let evaluate = || {
            residual(
                n,
                Euler {
                    eos: Nasg::ideal(1.4),
                },
                Transport {
                    power: 0.7,
                    bulk: 0.03,
                    ..Transport::constant(0.02, 0.75)
                },
                |x| {
                    [
                        1.0 + 0.05 * (TAU * x[0]).cos(),
                        0.1 * (TAU * x[1]).sin(),
                        0.1 * (TAU * x[0]).cos(),
                        0.1 * (TAU * x[2]).sin(),
                        1.0 + 0.05 * (TAU * x[2]).cos(),
                    ]
                },
            )
            .1
            .into_iter()
            .map(|v| v.map(f64::to_bits))
            .collect::<Vec<_>>()
        };
        assert_eq!(one.install(evaluate), three.install(evaluate));
    }
}

#[test]
fn viscous_timestep_includes_current_dirichlet_and_characteristic_flux_layers() {
    use compressible::Bc;
    for allspeed in [false, true] {
        for d in 0..3 {
            for characteristic in [false, true] {
                let sys = Euler {
                    eos: Nasg::ideal(1.4),
                };
                let mut shape = [1; 3];
                shape[d] = 16;
                let mut s = Solver::new(sys, Mesh::new(shape, [0.0; 3], [1.0; 3]), |x| {
                    let mut w = [1.0, 0.0, 0.0, 0.0, 0.4 * (1.0 + x[d])];
                    w[1 + d] = 3.0; // All waves leave the high-side boundary.
                    w
                });
                let mut far = [1.0, 0.0, 0.0, 0.0, if characteristic { 0.4 } else { 40.0 }];
                far[1 + d] = 3.0;
                let far = sys.to_cons(&far);
                s.bc[d] = [if characteristic {
                    Bc::Characteristic(far)
                } else {
                    Bc::Dirichlet(far)
                }; 2];
                s.viscous = Some(Transport {
                    power: 1.0,
                    ..Transport::constant(20.0, 1e30)
                });
                if allspeed {
                    s.allspeed = Some(AllSpeedOpts::default());
                }
                // Dirichlet: T=100. Supersonic characteristic: the third
                // outgoing flux layer extends the linear temperature field.
                let max_temperature = if characteristic {
                    2.0 + 2.5 * s.mesh.dx[d]
                } else {
                    100.0
                };
                let diffusion = 20.0 * max_temperature * (4.0 / 3.0);
                let expected = 0.3 * s.mesh.dx[d].powi(2) / diffusion;
                assert!(
                    (s.dt_cfl() / expected - 1.0).abs() < 1e-13,
                    "axis {d}, allspeed {allspeed}, characteristic {characteristic}"
                );
                for value in &mut s.u {
                    if value[0] == 0.0 {
                        *value = [f64::NAN; 5];
                    }
                }
                assert!((s.dt_cfl() / expected - 1.0).abs() < 1e-13);
            }
        }
    }
}

#[test]
fn hot_conductive_flux_is_bounded_by_the_light_neighbors_heat_capacity() {
    for allspeed in [false, true] {
        let mut s = Solver::new(
            Euler {
                eos: Nasg {
                    cv: 7.0,
                    ..Nasg::ideal(2.3)
                },
            },
            Mesh::new([16, 1, 1], [0.0; 3], [1.0; 3]),
            |x| {
                let (rho, temperature) = if x[0] < 0.5 {
                    (1.0, 1.0)
                } else {
                    (100.0, 1000.0)
                };
                [rho, 0.0, 0.0, 0.0, rho * 1.3 * 7.0 * temperature]
            },
        );
        s.viscous = Some(Transport {
            power: 1.0,
            ..Transport::constant(20.0, 0.1)
        });
        if allspeed {
            s.allspeed = Some(AllSpeedOpts::default());
        }
        let max_conductivity = 20.0 * 1000.0 * 2.3 * 7.0 / 0.1;
        let min_heat_capacity = 1.0 * 7.0;
        let expected =
            0.3 / ((4.0 / 3.0) * max_conductivity / min_heat_capacity * 16.0_f64.powi(2));
        assert!((s.dt_cfl() / expected - 1.0).abs() < 1e-13);
    }
}
