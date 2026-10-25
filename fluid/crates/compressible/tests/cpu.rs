use compressible::{eos::Nasg, multiphase::FiveEq, Bc, Mesh, Solver, System};
fn gas_liquid() -> Solver<7, FiveEq> {
    let sys = FiveEq {
        eos: [
            Nasg {
                p_inf: 1e6,
                ..Nasg::ideal(2.0)
            },
            Nasg::ideal(1.4),
        ],
        kapila: true,
    };
    Solver::new(sys, Mesh::new([8, 1, 1], [0.0; 3], [1.0; 3]), |_| {
        [
            0.99999 * 1000.0,
            0.00001 * 0.001,
            0.2,
            0.0,
            0.0,
            100.0,
            0.99999,
        ]
    })
}
#[test]
fn cfl_reflecting_wall_uses_current_boundary_before_first_step() {
    let mut s = gas_liquid();
    s.bc[0] = [Bc::Reflect; 2];
    let w = s.prim(0, 0, 0);
    let acoustic = s.cfl * s.mesh.dx[0] / s.sys.max_speed(&w, 0);
    let source = 0.5 / s.sys.dilatation_rate(&w, w[2] / s.mesh.dx[0]);
    let expected = acoustic.min(source);
    assert!(source < 0.1 * acoustic);
    assert!((s.dt_cfl() / expected - 1.0).abs() < 1e-13);
    // Deliberately poison only ghosts: the same physical state has the same CFL.
    for i in 0..s.mesh.m(0) {
        if i < s.mesh.g(0) || i >= s.mesh.g(0) + 8 {
            s.u[i] = [f64::NAN; 7];
        }
    }
    assert!((s.dt_cfl() / expected - 1.0).abs() < 1e-13);
}
#[test]
fn uniform_periodic_and_open_flow_do_not_get_spurious_boundary_limits() {
    for bc in [Bc::Periodic, Bc::Outflow] {
        let mut s = gas_liquid();
        s.bc[0] = [bc; 2];
        let w = s.prim(0, 0, 0);
        let expected = s.cfl * s.mesh.dx[0] / s.sys.max_speed(&w, 0);
        assert!((s.dt_cfl() / expected - 1.0).abs() < 1e-13);
    }
}

#[test]
fn characteristic_layers_extrapolate_outgoing_waves_and_reject_incoming_waves() {
    use compressible::{boundary::characteristic_ghost_layer, Euler};
    let sys = Euler {
        eos: Nasg::ideal(1.4),
    };
    let w = [1.0, 0.0, 0.0, 0.0, 1.0];
    let far = sys.to_cons(&w);
    for d in 0..3 {
        let mut l = [[0.0; 5]; 5];
        let mut r = [[0.0; 5]; 5];
        sys.eigen(&w, d, &mut l, &mut r);
        for side in 0..2 {
            let outgoing = if side == 0 { 0 } else { 4 };
            let incoming = 4 - outgoing;
            for mode in [outgoing, incoming] {
                let amplitude = 1e-4;
                let neighbor = std::array::from_fn(|k| w[k] + amplitude * r[k][mode]);
                let neighbor = sys.to_cons(&neighbor);
                for layer in [1.0, 2.0, 6.0] {
                    let actual = sys.to_prim(&characteristic_ghost_layer(
                        &sys, &far, &neighbor, &far, d, side, layer,
                    ));
                    for k in 0..5 {
                        let expected = w[k]
                            - if mode == outgoing {
                                layer * amplitude * r[k][mode]
                            } else {
                                0.0
                            };
                        assert!((actual[k] - expected).abs() < 2e-14);
                    }
                }
            }
        }
    }
}

#[test]
fn rk_updates_are_reproducible_for_every_2d_mesh_orientation() {
    use compressible::{Euler, FluxMode};
    use rayon::ThreadPoolBuilder;
    let one = ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    let three = ThreadPoolBuilder::new().num_threads(3).build().unwrap();
    for n in [[1, 96, 96], [96, 1, 96], [96, 96, 1]] {
        let evolve = || {
            let mut s = Solver::new(
                Euler {
                    eos: Nasg::ideal(1.4),
                },
                Mesh::new(n, [0.0; 3], [1.0; 3]),
                |x| {
                    [
                        1.0 + 0.01
                            * x.into_iter()
                                .map(|v| (std::f64::consts::TAU * v).sin())
                                .sum::<f64>(),
                        0.3,
                        0.2,
                        0.1,
                        1.0,
                    ]
                },
            );
            s.mode = FluxMode::Central;
            s.step(1e-4);
            s.step(1e-4);
            for i in 0..n[0] {
                for j in 0..n[1] {
                    for k in 0..n[2] {
                        assert!(s.sys.admissible(&s.prim(i, j, k)));
                    }
                }
            }
            s.u
        };
        assert_eq!(one.install(evolve), three.install(evolve));
    }
}

#[test]
fn strong_pressure_jump_reports_the_limiter_that_keeps_it_admissible() {
    use compressible::{Euler, FluxMode};
    let mut s = Solver::new(
        Euler {
            eos: Nasg::ideal(1.4),
        },
        Mesh::new([32, 1, 1], [0.0; 3], [1.0; 3]),
        |x| [1.0, 0.0, 0.0, 0.0, if x[0] < 0.5 { 1000.0 } else { 0.01 }],
    );
    s.mode = FluxMode::Central;
    s.bc[0] = [Bc::Outflow; 2];
    s.cfl = 0.2;
    s.step(s.dt_cfl());
    assert!(s.stats.limited_faces > 0);
    for i in 0..32 {
        assert!(s.sys.admissible(&s.prim(i, 0, 0)));
    }
}

#[test]
fn allspeed_first_pressure_gradient_limit_uses_current_dirichlet_boundary() {
    use compressible::{allspeed::AllSpeedOpts, Euler};
    let sys = Euler {
        eos: Nasg::ideal(1.4),
    };
    let mut s = Solver::new(sys, Mesh::new([8, 1, 1], [0.0; 3], [1.0; 3]), |_| {
        [1.0, 0.0, 0.0, 0.0, 1.0]
    });
    let opts = AllSpeedOpts::default();
    s.allspeed = Some(opts);
    s.bc[0] = [
        Bc::Dirichlet(sys.to_cons(&[1.0, 0.0, 0.0, 0.0, 100.0])),
        Bc::Outflow,
    ];
    let gradient = 99.0 / (2.0 * s.mesh.dx[0]);
    let expected = opts.cfl * (s.mesh.dx[0] / gradient).sqrt();
    let before = s.dt_cfl();
    assert!((before / expected - 1.0).abs() < 1e-13);
    for i in 0..s.mesh.m(0) {
        if i < s.mesh.g(0) || i >= s.mesh.g(0) + 8 {
            s.u[i] = [f64::NAN; 5];
        }
    }
    assert_eq!(before.to_bits(), s.dt_cfl().to_bits());
}

#[test]
fn allspeed_reflecting_wall_retains_the_multiphase_source_limit() {
    use compressible::allspeed::AllSpeedOpts;
    let mut s = gas_liquid();
    s.allspeed = Some(AllSpeedOpts::default());
    s.allspeed_stats.last_accel = 1.0;
    s.bc[0] = [Bc::Reflect; 2];
    let w = s.prim(0, 0, 0);
    let expected = 0.5 / s.sys.dilatation_rate(&w, w[2] / s.mesh.dx[0]);
    let before = s.dt_cfl();
    assert!((before / expected - 1.0).abs() < 1e-13);
    for i in 0..s.mesh.m(0) {
        if i < s.mesh.g(0) || i >= s.mesh.g(0) + 8 {
            s.u[i] = [f64::NAN; 7];
        }
    }
    assert_eq!(before.to_bits(), s.dt_cfl().to_bits());
}

#[test]
fn cached_thermal_properties_and_conserved_velocities_match_model_closures() {
    use compressible::{
        gpr::{Gpr, TauModel, NG},
        threephase::Kapila3,
        Euler,
    };
    fn check<const N: usize, S: System<N>>(sys: S, w: [f64; N]) {
        let u = sys.to_cons(&w);
        let p = sys.to_prim(&u);
        for d in 0..3 {
            assert_eq!(
                sys.velocity_conserved(&u, d).to_bits(),
                p[sys.vel_index() + d].to_bits()
            );
        }
        let ph = sys.phases(&p);
        assert_eq!(
            sys.temperature_from_phases(&p, &ph).to_bits(),
            sys.temperature(&p).to_bits()
        );
    }
    for density in [1e-6, 1.0, 1000.0] {
        check(
            Euler {
                eos: Nasg::ideal(1.4),
            },
            [density, 0.3, -0.2, 0.1, 1e5],
        );
        let sys = gas_liquid().sys;
        for a in [1e-10, 0.5, 1.0 - 1e-10] {
            check(sys, [a * density, (1.0 - a) * 1.2, 0.3, -0.2, 0.1, 1e5, a]);
        }
        check(
            Kapila3 {
                eos: [Nasg::ideal(1.4); 3],
            },
            [
                0.3 * density,
                0.2 * density,
                0.5 * density,
                0.3,
                -0.2,
                0.1,
                1e5,
                0.3,
                0.2,
            ],
        );
        let g = Gpr {
            eos: Nasg::ideal(1.4),
            cs: 1.0,
            tau: f64::INFINITY,
            tau_model: TauModel::Constant,
        };
        let mut w = [0.0; NG];
        w[0] = density;
        w[1..4].copy_from_slice(&[0.3, -0.2, 0.1]);
        w[13] = 1e5;
        for d in 0..3 {
            w[4 + 3 * d + d] = 1.0;
        }
        check(g, w);
    }
}

#[test]
fn periodic_axes_shorter_than_the_ghost_width_preserve_uniform_flow() {
    use compressible::{Euler, FluxMode};
    let expected = [1.0, 0.3, -0.2, 0.1, 1.0];
    for n in [[2, 1, 1], [1, 3, 1], [1, 1, 2], [2, 3, 4]] {
        for ng in [6, 10] {
            let mut s = Solver::new(
                Euler {
                    eos: Nasg::ideal(1.4),
                },
                Mesh::new(n, [0.0; 3], [1.0; 3]).with_ghosts(ng),
                |_| expected,
            );
            s.mode = FluxMode::Central;
            s.positivity = false;
            s.step(1e-5);
            for i in 0..n[0] {
                for j in 0..n[1] {
                    for k in 0..n[2] {
                        let w = s.prim(i, j, k);
                        assert!(
                            s.sys.admissible(&w),
                            "short periodic mesh {n:?}, ghosts={ng}"
                        );
                        for v in 0..5 {
                            assert!((w[v] - expected[v]).abs() < 1e-13);
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn parallel_pencils_preserve_uniform_scalars_for_each_mixture_carrier() {
    use compressible::{
        scalars::{Carrier, Scalar},
        FluxMode,
    };
    let sys = gas_liquid().sys;
    let mut s = Solver::new(sys, Mesh::new([12, 9, 8], [0.0; 3], [1.0; 3]), |x| {
        let a = 0.5 + 0.2 * (std::f64::consts::TAU * x[0]).sin();
        [a * 1000.0, (1.0 - a) * 1.2, 0.1, 0.05, -0.02, 1e5, a]
    });
    s.mode = FluxMode::Hybrid;
    for (carrier, value) in [
        (Carrier::Mixture, 0.25),
        (Carrier::Phase(0), 1.0),
        (Carrier::Phase(1), 0.7),
    ] {
        s.scalars
            .push(Scalar::new(&s, "uniform", carrier, 0.0, |_| value));
    }
    let totals: Vec<_> = (0..3).map(|q| s.scalar_total(q)).collect();
    for _ in 0..3 {
        s.step(1e-5);
    }
    for (q, value) in [0.25, 1.0, 0.7].into_iter().enumerate() {
        assert!((s.scalar_total(q) / totals[q] - 1.0).abs() < 1e-13);
        for i in 0..12 {
            for j in 0..9 {
                for k in 0..8 {
                    assert!((s.scalars[q].phi(&s, i, j, k) - value).abs() < 1e-12);
                }
            }
        }
    }
}

#[test]
fn short_reflecting_axes_preserve_rest_states_without_positivity_fallback() {
    use compressible::{Euler, FluxMode};
    for n in [[2, 1, 1], [1, 3, 1], [1, 1, 2], [2, 3, 4]] {
        for ng in [6, 10] {
            let mut s = Solver::new(
                Euler {
                    eos: Nasg::ideal(1.4),
                },
                Mesh::new(n, [0.0; 3], [1.0; 3]).with_ghosts(ng),
                |_| [1.0, 0.0, 0.0, 0.0, 1.0],
            );
            s.bc = [[Bc::Reflect; 2]; 3];
            s.mode = FluxMode::Central;
            s.positivity = false;
            for _ in 0..3 {
                s.step(1e-4);
            }
            for i in 0..n[0] {
                for j in 0..n[1] {
                    for k in 0..n[2] {
                        let w = s.prim(i, j, k);
                        assert!(s.sys.admissible(&w));
                        for (value, expected) in w.into_iter().zip([1.0, 0.0, 0.0, 0.0, 1.0]) {
                            assert!((value - expected).abs() < 1e-13);
                        }
                    }
                }
            }
        }
    }
}
