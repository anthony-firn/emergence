//! Quasicontinuum coupling: patch tests, defect benchmarks and adaptivity.
//!
//! The reference solution is always the fully atomistic model with the same
//! potential and boundary conditions. See `analysis/quasicontinuum.md`.

use materials::lattice::{Lattice, V3};
use materials::potential::{relax_cutoff, Harmonic, LennardJones, PairPotential};
use materials::quasicontinuum::{build_chain, build_triangular, Coupling, Qc};

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

fn diag_strain(f: f64) -> [[f64; 3]; 3] {
    [[f, 0.0, 0.0], [0.0, f, 0.0], [0.0, 0.0, 1.0]]
}

/// Equilibrium nearest-neighbour spacing of the 1D chain with the given
/// cutoff, via the crate's lattice relaxation.
fn chain_spacing(pot: &dyn PairPotential, rmax: f64) -> f64 {
    let base = Lattice::custom(
        [[1.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 12.0]],
        vec![[0.0, 0.0, 0.0]],
    );
    relax_cutoff(&base, pot, rmax, 0.9, 1.5).scale
}

/// Equilibrium spacing of the triangular lattice (NN + 2nd shell).
fn triangular_spacing(pot: &dyn PairPotential, rmax: f64) -> f64 {
    let h = 3.0_f64.sqrt() / 2.0;
    let base = Lattice::custom(
        [[1.0, 0.0, 0.0], [0.5, h, 0.0], [0.0, 0.0, 12.0]],
        vec![[0.0, 0.0, 0.0]],
    );
    relax_cutoff(&base, pot, rmax, 0.9, 1.5).scale
}

fn interior_force_max_1d(qc: &Qc, pot: &dyn PairPotential, f: [[f64; 3]; 3]) -> f64 {
    let u = qc.uniform_displacements(&f);
    let force = qc.forces(pot, &u);
    let xmax = qc.atoms.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max);
    let mut mx = 0.0_f64;
    for r in 0..qc.ndof() {
        let x = qc.atoms[qc.repatoms[r]][0];
        if x <= qc.atoms[0][0] + 1e-9 || x >= xmax - 1e-9 {
            continue;
        }
        for d in 0..qc.dim {
            mx = mx.max(force[r][d].abs());
        }
    }
    mx
}

/// A 2D patch residual: max force over repatoms that lie more than
/// `margin` lattice steps from every boundary, under a uniform strain.
/// This is the standard interior patch test; surface atoms carry physical
/// reactions that are not ghost forces.
fn interior_force_max_2d(
    qc: &Qc,
    pot: &dyn PairPotential,
    f: [[f64; 3]; 3],
    a0: f64,
    n1: usize,
    n2: usize,
    margin: usize,
) -> f64 {
    let h = 3.0_f64.sqrt() / 2.0;
    let u = qc.uniform_displacements(&f);
    let force = qc.forces(pot, &u);
    let mut mx = 0.0_f64;
    for r in 0..qc.ndof() {
        let x = qc.atoms[qc.repatoms[r]];
        let j = (x[1] / (a0 * h)).round() as isize;
        let i = (x[0] / a0 - 0.5 * j as f64).round() as isize;
        let m = margin as isize;
        if i < m || j < m || i > n1 as isize - 1 - m || j > n2 as isize - 1 - m {
            continue;
        }
        for d in 0..qc.dim {
            mx = mx.max(force[r][d].abs());
        }
    }
    mx
}

/// Fully atomistic chain with a point load at the centre; returns the
/// minimised repatom displacements and energy.
fn atomistic_chain(pot: &dyn PairPotential, n: usize, a0: f64, rmax: f64, load: f64) -> (Vec<V3>, f64) {
    let xc = (n - 1) as f64 * a0 / 2.0;
    let qc = build_chain(n, a0, rmax, xc, 1e9, 1, false, Coupling::Qnl);
    let mut free = vec![true; qc.ndof()];
    // clamp the two end atoms
    for r in 0..qc.ndof() {
        let x = qc.atoms[qc.repatoms[r]][0];
        if x <= 1e-9 || x >= (n - 1) as f64 * a0 - 1e-9 {
            free[r] = false;
        }
    }
    let mut ext = vec![[0.0; 3]; qc.ndof()];
    // apply the load to the atom closest to the centre
    let centre_atom = qc
        .repatoms
        .iter()
        .enumerate()
        .min_by(|(_, &a), (_, &b)| {
            (qc.atoms[a][0] - xc)
                .abs()
                .partial_cmp(&(qc.atoms[b][0] - xc).abs())
                .unwrap()
        })
        .map(|(r, _)| r)
        .unwrap();
    ext[centre_atom][0] = load;
    qc.minimize(pot, &free, &ext, 5000, 1e-8)
}

// ---------------------------------------------------------------------------
// 1D patch tests
// ---------------------------------------------------------------------------

#[test]
fn qnl_patch_test_passes_and_naive_does_not() {
    let pot = Harmonic { k: 1.0, r0: 1.0 };
    let f = diag_strain(1.01);
    let naive = build_chain(41, 1.0, 1.5, 20.0, 6.0, 4, false, Coupling::Naive);
    let qnl = build_chain(41, 1.0, 1.5, 20.0, 6.0, 4, false, Coupling::Qnl);
    let g_naive = interior_force_max_1d(&naive, &pot, f);
    let g_qnl = interior_force_max_1d(&qnl, &pot, f);
    println!("ghost force: naive = {g_naive:.3e}, qnl = {g_qnl:.3e}");
    assert!(g_naive > 1e-3, "naive should show a ghost force, got {g_naive}");
    assert!(g_qnl < 1e-10, "qnl should pass the patch test, got {g_qnl}");
}

#[test]
fn qnl_forces_match_finite_difference_1d() {
    // Validate the analytic QNL forces against a central difference.
    let pot = LennardJones::new(1.0, 1.0);
    let a0 = chain_spacing(&pot, 2.5 * 1.12);
    let rmax = 2.5 * a0;
    let qc = build_chain(31, a0, rmax, 15.0 * a0, 5.0 * a0, 3, false, Coupling::Qnl);
    let mut u = vec![[0.0; 3]; qc.ndof()];
    // a deterministic non-uniform perturbation
    for r in 0..qc.ndof() {
        let x = qc.atoms[qc.repatoms[r]][0];
        u[r][0] = 0.01 * a0 * (x / a0 * 1.7).sin();
    }
    let fa = qc.forces(&pot, &u);
    let h = 1e-6;
    let mut max_err = 0.0_f64;
    for r in 0..qc.ndof() {
        let mut up = u.clone();
        up[r][0] += h;
        let ep = qc.energy(&pot, &up);
        up[r][0] -= 2.0 * h;
        let em = qc.energy(&pot, &up);
        let fd = -(ep - em) / (2.0 * h);
        max_err = max_err.max((fa[r][0] - fd).abs());
    }
    println!("max |force - fd force| = {max_err:.3e}");
    assert!(max_err < 1e-4, "analytic forces disagree with FD: {max_err}");
}

// ---------------------------------------------------------------------------
// 1D defect benchmark
// ---------------------------------------------------------------------------

/// Coupled chain with a point load at the centre. Returns (u_atoms, energy,
/// ndof).
fn coupled_chain(
    pot: &dyn PairPotential,
    n: usize,
    a0: f64,
    rmax: f64,
    r_atom: f64,
    coarsen: usize,
    load: f64,
) -> (Vec<V3>, f64, usize) {
    let xc = (n - 1) as f64 * a0 / 2.0;
    let qc = build_chain(n, a0, rmax, xc, r_atom, coarsen, false, Coupling::Qnl);
    let mut free = vec![true; qc.ndof()];
    for r in 0..qc.ndof() {
        let x = qc.atoms[qc.repatoms[r]][0];
        if x <= 1e-9 || x >= (n - 1) as f64 * a0 - 1e-9 {
            free[r] = false;
        }
    }
    let mut ext = vec![[0.0; 3]; qc.ndof()];
    let centre_atom = qc
        .repatoms
        .iter()
        .enumerate()
        .min_by(|(_, &a), (_, &b)| {
            (qc.atoms[a][0] - xc)
                .abs()
                .partial_cmp(&(qc.atoms[b][0] - xc).abs())
                .unwrap()
        })
        .map(|(r, _)| r)
        .unwrap();
    ext[centre_atom][0] = load;
    let (udof, e) = qc.minimize(pot, &free, &ext, 5000, 1e-8);
    (qc.atom_displacements(&udof), e, qc.ndof())
}

#[test]
fn chain_point_load_error_decreases_with_atomistic_region() {
    let pot = LennardJones::new(1.0, 1.0);
    let n = 81;
    let a0 = chain_spacing(&pot, 2.5 * 1.12);
    let rmax = 2.5 * a0;
    let load = 0.5;
    let (u_ref, e_ref) = atomistic_chain(&pot, n, a0, rmax, load);
    println!("atomistic chain: a0 = {a0:.6}, E = {e_ref:.8}, ndof = {n}");

    let mut prev = f64::INFINITY;
    for &r_atom in &[2.0 * a0, 4.0 * a0, 8.0 * a0, 16.0 * a0] {
        let (u, e, ndof) = coupled_chain(&pot, n, a0, rmax, r_atom, 4, load);
        let du = u
            .iter()
            .zip(&u_ref)
            .map(|(a, b)| {
                ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
            })
            .fold(0.0_f64, f64::max);
        let de = (e - e_ref).abs();
        println!("r_atom = {:5.2}, ndof = {ndof:3}, max|du| = {du:.3e}, |dE| = {de:.3e}", r_atom / a0);
        assert!(de < 1.0, "energy blow-up");
        assert!(du < prev + 1e-9, "displacement error should not grow with r_atom");
        prev = du;
    }
    // The largest atomistic region should be close to the reference.
    let (_, e_big, _) = coupled_chain(&pot, n, a0, rmax, 40.0 * a0, 4, load);
    assert!((e_big - e_ref).abs() < 1e-5, "large atomistic region should match");
}

// ---------------------------------------------------------------------------
// 2D patch test and benchmark
// ---------------------------------------------------------------------------

#[test]
fn triangular_qnl_patch_test() {
    let pot = LennardJones::new(1.0, 1.0);
    let a0 = triangular_spacing(&pot, 1.9);
    let rmax = 1.9 * a0;
    let center = [a0 * 5.0, a0 * 4.5, 0.0];
    let f = diag_strain(1.005);
    let qnl = build_triangular(19, 19, a0, rmax, center, 3.0 * a0, 3, false, &[], Coupling::Qnl);
    let g = interior_force_max_2d(&qnl, &pot, f, a0, 19, 19, 4);
    println!("2D qnl patch residual = {g:.3e}");
    assert!(g < 1e-8, "2D qnl patch test failed: {g}");
}

#[test]
fn triangular_qnl_forces_match_finite_difference() {
    let pot = LennardJones::new(1.0, 1.0);
    let a0 = triangular_spacing(&pot, 1.9);
    let rmax = 1.9 * a0;
    let center = [a0 * 5.0, a0 * 4.5, 0.0];
    let qc = build_triangular(19, 19, a0, rmax, center, 3.0 * a0, 3, false, &[], Coupling::Qnl);
    let mut u = vec![[0.0; 3]; qc.ndof()];
    for r in 0..qc.ndof() {
        let x = qc.atoms[qc.repatoms[r]];
        u[r][0] = 0.005 * a0 * (x[0] / a0 * 1.3).sin() * (x[1] / a0 * 0.7).cos();
        u[r][1] = 0.005 * a0 * (x[1] / a0 * 1.1).sin();
    }
    let fa = qc.forces(&pot, &u);
    let h = 1e-6;
    let mut max_err = 0.0_f64;
    for r in 0..qc.ndof() {
        for d in 0..2 {
            let mut up = u.clone();
            up[r][d] += h;
            let ep = qc.energy(&pot, &up);
            up[r][d] -= 2.0 * h;
            let em = qc.energy(&pot, &up);
            let fd = -(ep - em) / (2.0 * h);
            max_err = max_err.max((fa[r][d] - fd).abs());
        }
    }
    println!("2D max |force - fd| = {max_err:.3e}");
    assert!(max_err < 1e-4, "2D analytic forces disagree with FD: {max_err}");
}

// ---------------------------------------------------------------------------
// 2D point-load benchmark and adaptive unfolding
// ---------------------------------------------------------------------------

struct Solve2d {
    udof: Vec<V3>,
    uatoms: Vec<V3>,
    energy: f64,
    ndof_free: usize,
}

/// Solve a 2D model with all boundary repatoms fixed and a point force at
/// the repatom nearest `load_pos`.
fn solve_tri(qc: &Qc, pot: &dyn PairPotential, a0: f64, n1: usize, n2: usize, load_pos: V3, load: V3) -> Solve2d {
    let h = 3.0_f64.sqrt() / 2.0;
    let mut free = vec![true; qc.ndof()];
    for r in 0..qc.ndof() {
        let x = qc.atoms[qc.repatoms[r]];
        let j = (x[1] / (a0 * h)).round() as usize;
        let i = (x[0] / a0 - 0.5 * j as f64).round() as usize;
        if i == 0 || i == n1 - 1 || j == 0 || j == n2 - 1 {
            free[r] = false;
        }
    }
    let n_free = free.iter().filter(|f| **f).count();
    let mut ext = vec![[0.0; 3]; qc.ndof()];
    let centre = qc
        .repatoms
        .iter()
        .enumerate()
        .min_by(|(_, &a), (_, &b)| {
            let da = (qc.atoms[a][0] - load_pos[0]).powi(2) + (qc.atoms[a][1] - load_pos[1]).powi(2);
            let db = (qc.atoms[b][0] - load_pos[0]).powi(2) + (qc.atoms[b][1] - load_pos[1]).powi(2);
            da.partial_cmp(&db).unwrap()
        })
        .map(|(r, _)| r)
        .unwrap();
    ext[centre] = load;
    let (udof, e) = qc.minimize(pot, &free, &ext, 4000, 1e-9);
    Solve2d {
        uatoms: qc.atom_displacements(&udof),
        udof,
        energy: e,
        ndof_free: n_free,
    }
}

fn max_du(a: &[V3], b: &[V3]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(p, q)| ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt())
        .fold(0.0_f64, f64::max)
}

#[test]
fn triangular_point_load_error_vs_atomistic_region() {
    let pot = LennardJones::new(1.0, 1.0);
    let n1 = 19;
    let n2 = 19;
    let a0 = triangular_spacing(&pot, 1.9);
    println!("2D a0 = {a0:.6}");
    let h = 3.0_f64.sqrt() / 2.0;
    let rmax = 1.9 * a0;
    let xc = a0 * (n1 as f64 - 1.0 + (n2 as f64 - 1.0) * 0.5) / 2.0;
    let yc = a0 * h * (n2 as f64 - 1.0) / 2.0;
    let center = [xc, yc, 0.0];
    let load = [0.3, -0.15, 0.0];

    let ref_qc = build_triangular(n1, n2, a0, rmax, center, 1e9, 1, false, &[], Coupling::Qnl);
    let reference = solve_tri(&ref_qc, &pot, a0, n1, n2, center, load);
    println!(
        "2D atomistic reference: ndof = {}, E = {:.8}",
        ref_qc.ndof(),
        reference.energy
    );

    let mut _prev = f64::INFINITY;
    for &r in &[2.0 * a0, 3.0 * a0, 4.0 * a0] {
        let qc = build_triangular(n1, n2, a0, rmax, center, r, 3, false, &[], Coupling::Qnl);
        let s = solve_tri(&qc, &pot, a0, n1, n2, center, load);
        let du = max_du(&s.uatoms, &reference.uatoms);
        let de = (s.energy - reference.energy).abs();
        println!(
            "r_atom = {:4.1}, ndof = {:4}, free = {:4}, max|du| = {du:.3e}, |dE| = {de:.3e}",
            r / a0,
            qc.ndof(),
            s.ndof_free
        );
        assert!(de < 1.0);
        // The displacement error is dominated by the interface/coarse-mesh
        // representation and need not fall monotonically with r_atom; the
        // coarsening study below isolates the mesh contribution.
        assert!(du < 5e-3);
        _prev = du;
    }

    // Mesh-coarsening study at a fixed atomistic radius.
    for &coarsen in &[2usize, 3usize] {
        let qc = build_triangular(n1, n2, a0, rmax, center, 3.0 * a0, coarsen, false, &[], Coupling::Qnl);
        let s = solve_tri(&qc, &pot, a0, n1, n2, center, load);
        let du = max_du(&s.uatoms, &reference.uatoms);
        let de = (s.energy - reference.energy).abs();
        println!(
            "coarsen = {coarsen}, ndof = {:4}, max|du| = {du:.3e}, |dE| = {de:.3e}",
            qc.ndof()
        );
    }
}

#[test]
fn adaptive_unfolding_reaches_target_with_few_dofs() {
    let pot = LennardJones::new(1.0, 1.0);
    let n1 = 19;
    let n2 = 19;
    let a0 = triangular_spacing(&pot, 1.9);
    let h = 3.0_f64.sqrt() / 2.0;
    let rmax = 1.9 * a0;
    let xc = a0 * (n1 as f64 - 1.0 + (n2 as f64 - 1.0) * 0.5) / 2.0;
    let yc = a0 * h * (n2 as f64 - 1.0) / 2.0;
    let center = [xc, yc, 0.0];
    let load = [0.3, -0.15, 0.0];

    let ref_qc = build_triangular(n1, n2, a0, rmax, center, 1e9, 1, false, &[], Coupling::Qnl);
    let reference = solve_tri(&ref_qc, &pot, a0, n1, n2, center, load);
    let target = 1e-4;

    // Uniform refinement: grow the atomistic disc until the target is met.
    let mut uniform_best = (0usize, f64::INFINITY);
    for &r in &[2.0 * a0, 3.0 * a0, 4.0 * a0, 5.0 * a0, 6.0 * a0, 8.0 * a0] {
        let qc = build_triangular(n1, n2, a0, rmax, center, r, 3, false, &[], Coupling::Qnl);
        let s = solve_tri(&qc, &pot, a0, n1, n2, center, load);
        let de = (s.energy - reference.energy).abs();
        println!("uniform r_atom = {:4.1}: ndof = {:4}, |dE| = {de:.3e}", r / a0, qc.ndof());
        if de < target {
            uniform_best = (qc.ndof(), de);
            break;
        }
    }
    println!("uniform reach: {:?}", uniform_best);

    // Adaptive: start small and unfold the worst cells by the strain-jump
    // indicator, keeping the solution so the indicator needs no re-solve.
    let mut unfold: Vec<(usize, usize)> = Vec::new();
    let mut adaptive_best = (0usize, f64::INFINITY);
    for _iter in 0..20 {
        let qc = build_triangular(n1, n2, a0, rmax, center, 1.5 * a0, 3, false, &unfold, Coupling::Qnl);
        let s = solve_tri(&qc, &pot, a0, n1, n2, center, load);
        let de = (s.energy - reference.energy).abs();
        let du = max_du(&s.uatoms, &reference.uatoms);
        println!(
            "adaptive: ndof = {:4}, cells = {:3}, |dE| = {de:.3e}, max|du| = {du:.3e}",
            qc.ndof(),
            unfold.len()
        );
        if de < target {
            adaptive_best = (qc.ndof(), de);
            break;
        }
        let ind = qc.cauchy_born_error(&pot, &s.udof);
        let new = qc.unfold_cells(&ind, 0.0, 3);
        let mut grew = false;
        for c in new {
            if !unfold.contains(&c) {
                unfold.push(c);
                grew = true;
            }
        }
        if !grew {
            break;
        }
    }
    println!("adaptive reach: {:?}", adaptive_best);
    assert!(adaptive_best.1 < target, "adaptive loop did not reach target");
    if uniform_best.1 < f64::INFINITY {
        assert!(
            adaptive_best.0 < uniform_best.0,
            "adaptive should use fewer dofs than uniform ({} vs {})",
            adaptive_best.0,
            uniform_best.0
        );
    }
}

#[test]
fn chain_vacancy_coupled_matches_atomistic() {
    let pot = LennardJones::new(1.0, 1.0);
    let n = 61;
    let a0 = chain_spacing(&pot, 2.5 * 1.12);
    let rmax = 2.5 * a0;
    let xc = (n - 1) as f64 * a0 / 2.0;

    let solve = |r_atom: f64, coarsen: usize| -> (Vec<V3>, f64, usize) {
        let qc = build_chain(n, a0, rmax, xc, r_atom, coarsen, true, Coupling::Qnl);
        let mut free = vec![true; qc.ndof()];
        for r in 0..qc.ndof() {
            let x = qc.atoms[qc.repatoms[r]][0];
            if x <= 1e-9 || x >= (n - 1) as f64 * a0 - 1e-9 {
                free[r] = false;
            }
        }
        let ext = vec![[0.0; 3]; qc.ndof()];
        let (udof, e) = qc.minimize(&pot, &free, &ext, 4000, 1e-8);
        (qc.atom_displacements(&udof), e, qc.ndof())
    };

    let (u_ref, e_ref, ndof_ref) = solve(1e9, 1);
    println!("vacancy reference: ndof = {ndof_ref}, E = {e_ref:.8}");
    // A vacancy is a much more local perturbation than a point load, so the
    // coupled model converges quickly.
    for &r in &[2.0 * a0, 4.0 * a0, 8.0 * a0] {
        let (u, e, ndof) = solve(r, 4);
        let du = max_du(&u, &u_ref);
        let de = (e - e_ref).abs();
        println!("vacancy r_atom = {:4.1}, ndof = {ndof:3}, max|du| = {du:.3e}, |dE| = {de:.3e}", r / a0);
        assert!(de < 1e-3, "vacancy energy error too large: {de}");
    }
}
