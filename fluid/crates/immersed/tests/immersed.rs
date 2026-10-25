//! Immersed geometry from triangle meshes against exact solutions.

use geometry::primitives::cylinder_z;
use geometry::Sdf;
use immersed::{Bc, Body, Grid, Solver};

/// Circular Couette flow between an immersed rotating inner cylinder and an
/// immersed stationary outer cylinder (an inside-out mesh: the exterior is
/// solid). Exact: `v_θ = A r + B/r`, `A = −Ω R₁²/(R₂² − R₁²)`,
/// `B = Ω R₁² R₂²/(R₂² − R₁²)`. Both walls are curved and cut the grid at
/// arbitrary positions, so this measures the immersed boundary's order.
#[test]
fn taylor_couette_between_immersed_cylinders_converges() {
    let (r1, r2, omega, nu) = (0.5, 1.0, 1.0, 0.1);
    let a = -omega * r1 * r1 / (r2 * r2 - r1 * r1);
    let b = omega * r1 * r1 * r2 * r2 / (r2 * r2 - r1 * r1);
    let vt = |r: f64| a * r + b / r;
    let half = 1.25;
    let mut rows = Vec::new();
    for n in [24usize, 48, 96] {
        let h = 2.0 * half / n as f64;
        let grid = Grid::new([n, n, 1], h, [-half, -half, 0.0]);
        // Extruded far beyond the one-cell depth; 4096 segments: sagitta
        // 3e-7 R, far below the discretisation error.
        let inner = cylinder_z([0.0, 0.0, -1.0], r1, 3.0, 4096);
        let mut outer = cylinder_z([0.0, 0.0, -1.0], r2, 3.0, 4096);
        outer.flip();
        let bodies = vec![
            Body {
                angular: [0.0, 0.0, omega],
                ..Body::fixed(Sdf::new(inner))
            },
            Body::fixed(Sdf::new_as_oriented(outer)),
        ];
        let walls = [Bc::wall(), Bc::wall()];
        let mut s = Solver::new(grid, nu, 1.0, [walls.clone(), walls, [Bc::Periodic, Bc::Periodic]], bodies);
        // Start from the exact field (inside the gap) so the run only has to
        // relax to the discrete steady state, which is what is measured.
        s.initialize(|x| {
            let r = (x[0] * x[0] + x[1] * x[1]).sqrt();
            if r <= r1 || r >= r2 {
                return [0.0; 3];
            }
            [-vt(r) * x[1] / r, vt(r) * x[0] / r, 0.0]
        });
        let (steps, rate) = s.run_until_steady(0.8, 40.0, 1e-8);
        assert!(rate < 1e-8, "n {n}: not steady ({rate:e})");
        // Error at every fluid u- and v-node.
        let (mut linf, mut l2, mut count) = (0.0f64, 0.0f64, 0usize);
        for c in 0..2 {
            let d = s.vel[c].d;
            for j in 0..d[1] as isize {
                for i in 0..d[0] as isize {
                    let x = grid.node(c, i, j, 0);
                    let r = (x[0] * x[0] + x[1] * x[1]).sqrt();
                    if r <= r1 || r >= r2 {
                        continue;
                    }
                    let exact = if c == 0 { -vt(r) * x[1] / r } else { vt(r) * x[0] / r };
                    let e = (s.vel[c].get(i, j, 0) - exact).abs();
                    linf = linf.max(e);
                    l2 += e * e;
                    count += 1;
                }
            }
        }
        let l2 = (l2 / count as f64).sqrt();
        // Torque the fluid exerts on the inner cylinder per unit length:
        // −4πρνB (it resists the rotation); the domain is one cell deep.
        let torque_exact = -4.0 * std::f64::consts::PI * nu * b * h;
        let f = s.body_force_momentum(0);
        let tq = s.body_torque_momentum(0)[2];
        println!(
            "n {n:3} (R₁/h = {:.1}): {steps} steps, L∞ {linf:.3e}, L2 {l2:.3e}, net force on inner {:.1e}; max div {:.1e}; IB mass source {:.1e}; torque error {:.2e}",
            r1 / h,
            (f[0] * f[0] + f[1] * f[1]).sqrt(),
            s.max_divergence(),
            s.ib_mass_source(),
            (tq - torque_exact).abs() / torque_exact.abs()
        );
        rows.push((linf, l2, (tq - torque_exact).abs() / torque_exact.abs()));
    }
    for w in rows.windows(2) {
        let (o_inf, o_2) = ((w[0].0 / w[1].0).log2(), (w[0].1 / w[1].1).log2());
        println!("orders: L∞ {o_inf:.2}, L2 {o_2:.2}, torque {:.2}", (w[0].2 / w[1].2).log2());
        assert!(o_2 > 1.5 && o_inf > 1.0, "orders L∞ {o_inf}, L2 {o_2}");
        assert!((w[0].2 / w[1].2).log2() > 1.5, "torque order");
    }
    assert!(rows[2].2 < 0.01, "torque error {:.2e} at R₁/h = 19", rows[2].2);
}
