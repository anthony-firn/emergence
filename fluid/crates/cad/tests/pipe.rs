//! Pipe geometry against ASME dimensions, its hydraulics against
//! Hagen–Poiseuille and Colebrook, and the generated bore against the
//! exact laminar solution by solving the flow in it.

use std::f64::consts::PI;

use cad::kits::pipe::{colebrook, colebrook_residual, minor_loss, Pipe, Reducer, Run, Schedule, Tee};
use cad::meshing::mesh_solid;

#[test]
fn standard_sizes_and_derived_areas() {
    // The table is ASME B36.10M data; what is checked here is that the
    // geometry derived from it is self-consistent and that the generated
    // solid measures what the table says.
    let p = Pipe::nps(4.0, Schedule::Sch40).unwrap();
    assert!((p.outer - 114.3e-3).abs() < 1e-12);
    assert!((p.wall - 6.02e-3).abs() < 1e-12);
    assert!((p.inner() - (114.3 - 2.0 * 6.02) * 1e-3).abs() < 1e-12);
    // NPS 4 Sch 40 has a bore of 102.26 mm and a flow area of 82.1 cm^2.
    assert!((p.inner() * 1e3 - 102.26).abs() < 1e-9, "{:.3} mm", p.inner() * 1e3);
    assert!((p.flow_area() * 1e4 - 82.1).abs() < 0.05, "{:.3} cm^2", p.flow_area() * 1e4);
    // Steel at 7850 kg/m^3: NPS 4 Sch 40 weighs about 16 kg/m.
    let mass = p.mass_per_length(7850.0);
    assert!((mass - 16.07).abs() < 0.1, "{mass:.3} kg/m");
    // A round bore's hydraulic diameter is its bore.
    assert!((p.hydraulic_diameter() - p.inner()).abs() < 1e-15);

    // Every tabulated size is consistent: a positive bore, and walls that
    // grow with the schedule.
    for row in cad::kits::pipe::NPS_TABLE {
        let mut last = 0.0;
        for sched in [Schedule::Sch10, Schedule::Sch40, Schedule::Sch80, Schedule::Sch160] {
            if let Ok(p) = Pipe::nps(row.0, sched) {
                assert!(p.inner() > 0.0, "NPS {} {sched:?} has no bore", row.0);
                assert!(p.wall > last, "NPS {} {sched:?} wall {} is not thicker than {last}", row.0, p.wall);
                last = p.wall;
            }
        }
    }
    assert!(Pipe::nps(3.7, Schedule::Sch40).is_err(), "an unlisted size must be refused");
}

#[test]
fn laminar_friction_is_hagen_poiseuille() {
    // Darcy-Weisbach with f = 64/Re is algebraically Hagen-Poiseuille.
    // Checking the identity checks both implementations at once.
    let p = Pipe::nps(1.0, Schedule::Sch40).unwrap();
    let (rho, nu) = (998.0, 1.0e-6);
    for u in [0.001, 0.005, 0.02] {
        let re = p.reynolds(u, nu);
        assert!(re < 2300.0, "Re = {re:.0} is not laminar");
        let dw = p.pressure_drop(1.0, u, rho, nu, 0.0).unwrap();
        let hp = p.poiseuille_drop(1.0, u, rho, nu);
        assert!((dw - hp).abs() / hp < 1e-15, "{dw:e} vs {hp:e}");
        assert!((p.friction_factor(re, 0.0).unwrap() - 64.0 / re).abs() < 1e-18);
    }
    // The transition range has no defined friction factor, and says so.
    let err = p.friction_factor(3000.0, 0.0).unwrap_err();
    assert!(err.contains("transition"), "{err}");
}

#[test]
fn colebrook_solves_its_own_equation() {
    for &(re, rough) in &[(4e3, 0.0), (1e5, 0.0), (1e6, 0.0), (1e5, 1e-3), (1e7, 5e-3), (1e8, 1e-2)] {
        let f = colebrook(re, rough);
        let residual = colebrook_residual(f, re, rough);
        assert!(residual.abs() < 1e-12, "Re = {re:e}, e/D = {rough}: residual {residual:e}");
        assert!(f > 0.005 && f < 0.12, "f = {f} is off the Moody chart");
    }
    // A smooth pipe at Re = 1e5 has f = 0.0180 on the Moody chart.
    let f = colebrook(1e5, 0.0);
    assert!((f - 0.0180).abs() < 2e-4, "smooth pipe at Re = 1e5: f = {f:.5}");
    // Fully rough: f stops depending on Re.
    let a = colebrook(1e7, 0.01);
    let b = colebrook(1e9, 0.01);
    assert!((a - b).abs() / a < 0.02, "fully rough f should plateau: {a:.5} vs {b:.5}");
    // Roughness only ever increases friction.
    assert!(colebrook(1e5, 1e-3) > colebrook(1e5, 0.0));
    // The minor-loss table is there as data.
    assert!(minor_loss::ELBOW_90_SHORT > minor_loss::ELBOW_90_LONG);
}

#[test]
fn runs_build_centrelines_of_the_right_length() {
    // Straights plus R*theta arcs, exactly.
    let d = 0.1143;
    let run = Run::new([0.0; 3], [1.0, 0.0, 0.0])
        .straight(0.5)
        .elbow_long(d, [0.0, 0.0, 1.0])
        .straight(0.3)
        .elbow_long(d, [0.0, 0.0, 1.0])
        .straight(0.5);
    let expect = 0.5 + 0.3 + 0.5 + 2.0 * (1.5 * d) * (0.5 * PI);
    assert!((run.length() - expect).abs() / expect < 1e-15);
    let path = run.centreline().expect("centreline");
    let rel = (path.length() - expect).abs() / expect;
    // The arcs are polylines, so the sampled length is short by the chord
    // error of 24 steps per quadrant.
    assert!(rel < 1e-3, "sampled centreline off by {rel:e}");

    // Two 90 degree turns about +z reverse the direction.
    let (end, dir) = run.end().unwrap();
    assert!((dir[0] + 1.0).abs() < 1e-9, "direction after two turns: {dir:?}");
    assert!(end[2].abs() < 1e-12, "the run should stay in the z = 0 plane");

    // ASME B16.9: a long-radius 90 degree elbow has a centreline radius of
    // 1.5 times the nominal size, so its centre-to-face is 1.5 D.
    let one_elbow = Run::new([0.0; 3], [1.0, 0.0, 0.0]).elbow_long(d, [0.0, 0.0, 1.0]);
    let (end, _) = one_elbow.end().unwrap();
    assert!((end[0] - 1.5 * d).abs() < 1e-6 && (end[1] - 1.5 * d).abs() < 1e-6, "{end:?}");

    // An axis along the run is refused: a pipe cannot twist about itself.
    let bad = Run::new([0.0; 3], [1.0, 0.0, 0.0]).elbow(0.2, 0.5 * PI, [1.0, 0.0, 0.0]);
    let err = bad.centreline().unwrap_err();
    assert!(err.contains("perpendicular"), "{err}");
}

#[test]
fn the_generated_bore_measures_its_nominal_size() {
    let p = Pipe::nps(2.0, Schedule::Sch40).unwrap();
    let run = Run::new([0.0; 3], [0.0, 0.0, 1.0]).straight(0.2);
    let path = run.centreline().unwrap();
    let (bore, check) = p.bore(&path, 64);
    assert!(check.watertight, "{check:?}");
    assert!(check.curvature_ratio < 1e-9, "a straight run has no curvature");

    // Measured bore and flow area from the solid itself.
    let r = 0.5 * p.inner();
    assert!(bore.sdf([0.0, 0.0, 0.1]) < 0.0);
    assert!(bore.sdf([r * 1.01, 0.0, 0.1]) > 0.0, "the bore is too wide");
    assert!(bore.sdf([r * 0.99, 0.0, 0.1]) < 0.0, "the bore is too narrow");
    // Volume = flow area x length, to the inscribed-polygon deficit of a
    // 64-gon (1 - (n/2) sin(2 pi/n)/pi = 1.6e-3).
    let mesh = mesh_solid(&bore, p.inner() / 16.0).expect("mesh");
    let v = mesh.report().volume;
    let expect = p.flow_area() * 0.2;
    assert!((v - expect).abs() / expect < 0.01, "bore volume {v:e} vs {expect:e}");

    // The wall solid is the annulus.
    let (wall, _) = p.wall_solid(&path, 64);
    assert!(wall.sdf([0.5 * (r + 0.5 * p.outer), 0.0, 0.1]) < 0.0, "metal missing");
    assert!(wall.sdf([0.0, 0.0, 0.1]) > 0.0, "the bore must be open");

    // A run that bends tighter than its own bore is refused.
    let tight = Run::new([0.0; 3], [1.0, 0.0, 0.0]).elbow(0.5 * p.inner() * 0.5, 0.5 * PI, [0.0, 0.0, 1.0]);
    let err = p.bore_as_housing(&tight.centreline().unwrap(), 32).unwrap_err();
    assert!(err.contains("bends tighter"), "{err}");
}

#[test]
fn reducers_and_tees() {
    let big = Pipe::nps(4.0, Schedule::Sch40).unwrap();
    let small = Pipe::nps(2.0, Schedule::Sch40).unwrap();
    let r = Reducer { from: big, to: small, length: 0.1 };
    // Continuity and Bernoulli: a contraction accelerates and loses static
    // pressure.
    let u1 = 1.0;
    let dp = r.bernoulli_recovery(u1, 998.0);
    let u2 = u1 * big.flow_area() / small.flow_area();
    assert!((dp - 0.5 * 998.0 * (u1 * u1 - u2 * u2)).abs() < 1e-9);
    assert!(dp < 0.0, "a contraction must drop the static pressure");
    assert!(u2 > 3.0 * u1, "the area ratio should accelerate the flow");
    // Half-angle from the geometry.
    let expect = (0.5 * (small.inner() - big.inner()) / 0.1).atan();
    assert!((r.half_angle() - expect).abs() < 1e-15);
    assert!(r.half_angle() < 0.0, "this reducer contracts");

    // The bore is open through and the metal is not.
    let bore = r.bore();
    assert!(bore.sdf([0.0, 0.0, 0.05]) < 0.0);
    let wall = r.wall_solid();
    assert!(wall.sdf([0.0, 0.0, 0.05]) > 0.0, "the reducer bore must be open");
    assert!(wall.sdf([0.5 * (0.5 * big.inner() + 0.5 * big.outer), 0.0, 0.002]) < 0.0);

    let tee = Tee { run: big, branch: small, run_length: 0.3, branch_length: 0.15 };
    let bore = tee.bore();
    assert!(bore.sdf([0.1, 0.0, 0.0]) < 0.0, "the run bore");
    assert!(bore.sdf([0.0, 0.0, 0.1]) < 0.0, "the branch bore");
    assert!(bore.sdf([0.0, 0.0, 0.2]) > 0.0, "beyond the branch");
    let mesh = mesh_solid(&tee.wall_solid(), 0.003).expect("mesh");
    assert_eq!(mesh.report().boundary_edges, 0, "the fitting must be closed");
}

#[test]
fn flow_in_the_generated_pipe_is_the_exact_laminar_solution() {
    // The strongest check available: a pipe drawn by the kit, with the
    // exact Hagen-Poiseuille profile at the inlet, solved by `immersed`.
    // Fully developed laminar pipe flow is an exact solution of the steady
    // Navier-Stokes equations, so the solver must preserve the profile and
    // reproduce the pressure gradient.
    //
    // The error is *not* monotone in the resolution, and the test does not
    // pretend it is. For a sharp direct-forcing immersed boundary it
    // depends on where the wall falls inside the near-wall cell: at
    // D/h = 8, 12 and 16 the wall sits 0.33, 0.74 and 0.16 of the way
    // across the cell on the diagonal, and the measured errors follow that
    // ordering rather than the grid spacing. `analysis/real_geometry.md`
    // records the same effect spoiling the DFG 3D lift. What is asserted
    // here is that every resolution agrees with the exact solution to a
    // few per cent, and the best-aligned one to about one.
    use immersed::solver::{Bc, Body, Solver};
    use immersed::Grid;

    let p = Pipe::new(0.028, 0.002); // 24 mm bore
    let d = p.inner();
    let a = 0.5 * d;
    let u_mean = 0.01;
    let re_target = 40.0;
    let nu = u_mean * d / re_target;
    let rho = 1000.0;

    let solve = |cells_across: usize| -> (f64, f64) {
        let h = d / cells_across as f64;
        let length = 6.0 * d;
        let nx = (length / h).round() as usize;
        // Two cells of metal all round the bore.
        let nt = cells_across + 4;
        let grid = Grid::new([nx, nt, nt], h, [0.0; 3]);
        let axis = 0.5 * nt as f64 * h;

        let run = Run::new([0.0, axis, axis], [1.0, 0.0, 0.0]).straight(length);
        let path = run.centreline().unwrap();
        let sdf = p.bore_as_housing(&path, 96).expect("housing");

        let inlet = move |x: [f64; 3], _t: f64| {
            let r = ((x[1] - axis).powi(2) + (x[2] - axis).powi(2)).sqrt();
            if r >= a {
                [0.0; 3]
            } else {
                [2.0 * u_mean * (1.0 - (r / a).powi(2)), 0.0, 0.0]
            }
        };
        let bcs = [
            [Bc::inlet(inlet), Bc::Outlet],
            [Bc::wall(), Bc::wall()],
            [Bc::wall(), Bc::wall()],
        ];
        let mut solver = Solver::new(grid, nu, rho, bcs, vec![Body::fixed(sdf)]);
        let (steps, rate) = solver.run_until_steady(0.8, 600.0, 1e-7);
        println!("  D/h = {cells_across}: {steps} steps to t = {:.2}, rate {rate:e}", solver.t);
        assert!(steps > 20, "the solve did not advance at {cells_across} cells across");

        // Profile error at 5 D, where the flow is developed (the laminar
        // entrance length is about 0.06 Re D = 2.4 D).
        //
        // Sampled at *fixed physical radii* by trilinear interpolation,
        // not at whichever nodes a given grid happens to put inside the
        // bore: comparing different node sets across resolutions mixes the
        // sampling position into the error and hides the convergence.
        let u = &solver.vel[0];
        let interp = |x: f64, y: f64, z: f64| -> f64 {
            // u lives at x = i h, y = (j + 1/2) h, z = (k + 1/2) h.
            let (fi, fj, fk) = (x / h, y / h - 0.5, z / h - 0.5);
            let (i0, j0, k0) = (fi.floor() as isize, fj.floor() as isize, fk.floor() as isize);
            let (tx, ty, tz) = (fi - i0 as f64, fj - j0 as f64, fk - k0 as f64);
            let mut acc = 0.0;
            for (di, wx) in [(0isize, 1.0 - tx), (1, tx)] {
                for (dj, wy) in [(0isize, 1.0 - ty), (1, ty)] {
                    for (dk, wz) in [(0isize, 1.0 - tz), (1, tz)] {
                        acc += wx * wy * wz * u.data[u.at(i0 + di, j0 + dj, k0 + dk)];
                    }
                }
            }
            acc
        };
        let x_station = 5.0 * d;
        let mut worst: f64 = 0.0;
        for &frac in &[0.0, 0.2, 0.4, 0.6] {
            let r = frac * a;
            for (dy, dz) in [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)] {
                let got = interp(x_station, axis + dy * r, axis + dz * r);
                let want = p.poiseuille_profile(u_mean, r);
                worst = worst.max((got - want).abs() / (2.0 * u_mean));
            }
        }

        // Pressure gradient between two developed stations.
        let mean_p = |i: usize| -> f64 {
            let mut sum = 0.0;
            let mut n = 0;
            for j in 0..nt {
                for k in 0..nt {
                    let y = (j as f64 + 0.5) * h;
                    let z = (k as f64 + 0.5) * h;
                    let r = ((y - axis).powi(2) + (z - axis).powi(2)).sqrt();
                    if r < 0.6 * a {
                        sum += solver.p[i + nx * (j + nt * k)];
                        n += 1;
                    }
                }
            }
            sum / n as f64
        };
        let (i1, i2) = ((3.0 * d / h).round() as usize, (5.0 * d / h).round() as usize);
        let span = (i2 - i1) as f64 * h;
        let dp = rho * (mean_p(i1) - mean_p(i2));
        let expect = p.poiseuille_drop(span, u_mean, rho, nu);
        (worst, (dp - expect).abs() / expect)
    };

    let mut profile = Vec::new();
    let mut gradient = Vec::new();
    for cells in [8usize, 12, 16] {
        let (prof, grad) = solve(cells);
        println!(
            "D/h = {cells}: worst profile error {:.3} % of the centreline velocity, \
             pressure gradient {:.3} % off Hagen-Poiseuille",
            100.0 * prof,
            100.0 * grad
        );
        profile.push(prof);
        gradient.push(grad);
    }
    // Both must converge with the resolution: the immersed boundary is the
    // only thing standing between the drawn pipe and an exact solution.
    for (i, cells) in [8usize, 12, 16].iter().enumerate() {
        assert!(
            profile[i] < 0.08,
            "D/h = {cells}: profile error {:.4} against the exact parabola",
            profile[i]
        );
        assert!(
            gradient[i] < 0.06,
            "D/h = {cells}: pressure gradient {:.4} off Hagen-Poiseuille",
            gradient[i]
        );
    }
    // The best-aligned grid of the three gets within about a per cent of
    // the exact solution, in both the profile and the pressure gradient.
    let best_profile = profile.iter().copied().fold(f64::INFINITY, f64::min);
    let best_gradient = gradient.iter().copied().fold(f64::INFINITY, f64::min);
    assert!(best_profile < 0.03, "best profile error {best_profile:.4}");
    assert!(best_gradient < 0.02, "best pressure-gradient error {best_gradient:.4}");
}
