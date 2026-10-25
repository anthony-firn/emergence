//! Cauchy-Born elasticity: apply a homogeneous deformation gradient to the
//! lattice and take second derivatives of the energy density.
//!
//! For a deformation gradient `F` the reference bond `R` deforms to `F R`.
//! With internal basis displacements `u` the bond becomes
//! `F R + u_j − u_i`, and the energy density per reference volume is
//!
//! ```text
//! W(F, u) = (1/(2 V0)) Σ_{bonds} φ(|F R + u_j − u_i|).
//! ```
//!
//! The frozen (affine) stiffness is the second derivative at `u = 0`; for
//! non-Bravais lattices the basis relaxes for each strain, and the physical
//! stiffness is the Schur complement
//!
//! ```text
//! C = W_FF − W_Fu K⁻¹ W_uF,   K = W_uu (reduced by the translation gauge),
//! ```
//!
//! evaluated at zero stress (the equilibrium lattice constant). The tests
//! check the analytic `C` against a finite difference of the *relaxed*
//! energy `W(F, u*(F))`.
//!
//! Benchmarks: for a pair potential on a Bravais lattice at zero pressure
//! the Cauchy relation `C12 = C44` holds exactly, cubic symmetry forces
//! `C11 = C22 = C33` (and likewise for the off-diagonals), and for LJ fcc
//! the 6x6 Voigt matrix is reported in `ε/σ³` and matched to the
//! finite-difference result.

use crate::lattice::{dot, norm, Lattice, V3};
use crate::linalg::solve;
use crate::potential::{bond_stiffness, lattice_sum, PairPotential};

/// Voigt index: 0=xx, 1=yy, 2=zz, 3=yz, 4=xz, 5=xy.
pub fn voigt(i: usize, j: usize) -> usize {
    match (i, j) {
        (0, 0) => 0,
        (1, 1) => 1,
        (2, 2) => 2,
        (1, 2) | (2, 1) => 3,
        (0, 2) | (2, 0) => 4,
        (0, 1) | (1, 0) => 5,
        _ => panic!("bad Voigt indices {i},{j}"),
    }
}

fn mat_vec(f: &[[f64; 3]; 3], r: V3) -> V3 {
    [
        f[0][0] * r[0] + f[0][1] * r[1] + f[0][2] * r[2],
        f[1][0] * r[0] + f[1][1] * r[1] + f[1][2] * r[2],
        f[2][0] * r[0] + f[2][1] * r[1] + f[2][2] * r[2],
    ]
}

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// Symmetric strain basis `D_I` in Voigt order, so that
/// `F = I + Σ_I e_I D_I` and `C_IJ = ∂²W/∂e_I∂e_J`.
pub fn strain_basis(i: usize) -> [[f64; 3]; 3] {
    let mut d = [[0.0; 3]; 3];
    match i {
        0 => d[0][0] = 1.0,
        1 => d[1][1] = 1.0,
        2 => d[2][2] = 1.0,
        3 => {
            d[1][2] = 0.5;
            d[2][1] = 0.5;
        }
        4 => {
            d[0][2] = 0.5;
            d[2][0] = 0.5;
        }
        5 => {
            d[0][1] = 0.5;
            d[1][0] = 0.5;
        }
        _ => panic!("bad strain index {i}"),
    }
    d
}

/// Energy density `W(F, u)` in units of energy per reference volume.
pub fn energy_density(
    lat: &Lattice,
    pot: &dyn PairPotential,
    rmax: f64,
    f: &[[f64; 3]; 3],
    u: &[V3],
) -> f64 {
    let bonds = lat.bond_list(rmax);
    let n = lat.natoms();
    let v0 = lat.volume();
    let mut e = 0.0;
    for i in 0..n {
        for b in &bonds[i] {
            let fr = mat_vec(f, b.r);
            let s = [
                fr[0] + u[b.j][0] - u[i][0],
                fr[1] + u[b.j][1] - u[i][1],
                fr[2] + u[b.j][2] - u[i][2],
            ];
            e += pot.energy(norm(s));
        }
    }
    e / (2.0 * v0)
}

struct Hessians {
    /// ∂²W/∂F_ab∂F_cd, row/col flat `(a*3+b)`.
    wff: [[f64; 9]; 9],
    /// ∂²W/∂F_ab∂u_r, row `(a*3+b)`, column internal index.
    wfu: Vec<[f64; 9]>,
    /// ∂²W/∂u_r∂u_s reduced (atoms 1..n fixed, atom 0 gauge).
    kuu: Vec<f64>,
    reduce: usize,
}

fn hessians(
    lat: &Lattice,
    pot: &dyn PairPotential,
    rmax: f64,
    f: &[[f64; 3]; 3],
    u: &[V3],
) -> Hessians {
    let bonds = lat.bond_list(rmax);
    let n = lat.natoms();
    let v0 = lat.volume();
    // Translation gauge: atom 0 is held fixed, atoms 1..n-1 are free, so
    // there are 3(n-1) internal coordinates (zero for a Bravais lattice).
    let reduce = 3 * n.saturating_sub(1);
    // Internal index of an atom, if it is a free coordinate (atom 0 fixed).
    let idx = |atom: usize, comp: usize| -> Option<usize> {
        if atom == 0 {
            None
        } else {
            Some((atom - 1) * 3 + comp)
        }
    };
    let mut wff = [[0.0; 9]; 9];
    let mut wfu = vec![[0.0; 9]; reduce];
    let mut kuu = vec![0.0; reduce * reduce];
    let half_v0 = 2.0 * v0;
    for i in 0..n {
        for b in &bonds[i] {
            let fr = mat_vec(f, b.r);
            let s = [
                fr[0] + u[b.j][0] - u[i][0],
                fr[1] + u[b.j][1] - u[i][1],
                fr[2] + u[b.j][2] - u[i][2],
            ];
            let m = bond_stiffness(pot, s);
            // W_FF[(a,b)][(c,d)] += R_b R_d M_{ac}.
            for a in 0..3 {
                for bb in 0..3 {
                    for c in 0..3 {
                        for dd in 0..3 {
                            wff[a * 3 + bb][c * 3 + dd] += b.r[bb] * b.r[dd] * m[a][c] / half_v0;
                        }
                    }
                }
            }
            // W_Fu += R_b M_{aβ} at the j end, − at the i end.
            for a in 0..3 {
                for bb in 0..3 {
                    let row = a * 3 + bb;
                    for beta in 0..3 {
                        if let Some(jj) = idx(b.j, beta) {
                            wfu[jj][row] += b.r[bb] * m[a][beta] / half_v0;
                        }
                        if let Some(ii) = idx(i, beta) {
                            wfu[ii][row] -= b.r[bb] * m[a][beta] / half_v0;
                        }
                    }
                }
            }
            // W_uu: +M on the diagonal, −M across the bond.
            for alpha in 0..3 {
                for beta in 0..3 {
                    let mm = m[alpha][beta] / half_v0;
                    if let Some(ii) = idx(i, alpha) {
                        if let Some(jj) = idx(i, beta) {
                            kuu[ii * reduce + jj] += mm;
                        }
                        if let Some(jj) = idx(b.j, beta) {
                            kuu[ii * reduce + jj] -= mm;
                        }
                    }
                    if let Some(ii) = idx(b.j, alpha) {
                        if let Some(jj) = idx(b.j, beta) {
                            kuu[ii * reduce + jj] += mm;
                        }
                        if let Some(jj) = idx(i, beta) {
                            kuu[ii * reduce + jj] -= mm;
                        }
                    }
                }
            }
        }
    }
    Hessians {
        wff,
        wfu,
        kuu,
        reduce,
    }
}

/// Internal gradient `∂W/∂u` at fixed `F`.
pub fn internal_gradient(
    lat: &Lattice,
    pot: &dyn PairPotential,
    rmax: f64,
    f: &[[f64; 3]; 3],
    u: &[V3],
) -> Vec<f64> {
    let bonds = lat.bond_list(rmax);
    let n = lat.natoms();
    let half_v0 = 2.0 * lat.volume();
    let mut g = vec![0.0; 3 * n];
    for i in 0..n {
        for b in &bonds[i] {
            let fr = mat_vec(f, b.r);
            let s = [
                fr[0] + u[b.j][0] - u[i][0],
                fr[1] + u[b.j][1] - u[i][1],
                fr[2] + u[b.j][2] - u[i][2],
            ];
            let rad = norm(s);
            let force = pot.d1(rad) / rad;
            for d in 0..3 {
                g[3 * b.j + d] += force * s[d] / half_v0;
                g[3 * i + d] -= force * s[d] / half_v0;
            }
        }
    }
    g
}

/// Relax the internal basis coordinates for a fixed deformation gradient
/// (atom 0 held fixed). Newton iteration on the analytic Hessian.
pub fn relax_internal(
    lat: &Lattice,
    pot: &dyn PairPotential,
    rmax: f64,
    f: &[[f64; 3]; 3],
) -> Vec<V3> {
    let n = lat.natoms();
    let mut u = vec![[0.0; 3]; n];
    if n <= 1 {
        return u;
    }
    let reduce = 3 * n - 3;
    for _ in 0..50 {
        let h = hessians(lat, pot, rmax, f, &u);
        let g = internal_gradient(lat, pot, rmax, f, &u);
        // Reduced gradient: components belonging to atoms 1..n-1.
        let gr: Vec<f64> = (3..3 * n).map(|k| g[k]).collect();
        let norm_g = gr.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm_g < 1e-13 {
            break;
        }
        let rhs: Vec<f64> = gr.iter().map(|x| -x).collect();
        let du = solve(h.kuu, rhs, reduce);
        let mut step: f64 = 0.0;
        for k in 0..reduce {
            step = step.max(du[k].abs());
        }
        for k in 0..reduce {
            let atom = 1 + k / 3;
            u[atom][k % 3] += du[k];
        }
        if step < 1e-14 {
            break;
        }
    }
    u
}

/// Fully relaxed energy density `W(F, u*(F))`.
pub fn relaxed_energy_density(
    lat: &Lattice,
    pot: &dyn PairPotential,
    rmax: f64,
    f: &[[f64; 3]; 3],
) -> f64 {
    let u = relax_internal(lat, pot, rmax, f);
    energy_density(lat, pot, rmax, f, &u)
}

/// Cauchy-Born 6x6 Voigt stiffness with internal relaxation.
pub fn elastic_tensor(lat: &Lattice, pot: &dyn PairPotential, rmax: f64) -> [[f64; 6]; 6] {
    let n = lat.natoms();
    let u0 = vec![[0.0; 3]; n];
    let h = hessians(lat, pot, rmax, &IDENTITY, &u0);
    let mut m9 = h.wff;
    if h.reduce > 0 {
        // corr = W_Fu K⁻¹ W_uF, a 9x9 correction.
        // Solve K y_q = W_Fu[q] for each column q, then corr[p][q] = W_Fu[p]·y_q.
        let mut solved: Vec<Vec<f64>> = Vec::with_capacity(9);
        for q in 0..9 {
            let rhs = h.wfu.iter().map(|row| row[q]).collect::<Vec<f64>>();
            solved.push(solve(h.kuu.clone(), rhs, h.reduce));
        }
        for p in 0..9 {
            for q in 0..9 {
                let mut s = 0.0;
                for a in 0..h.reduce {
                    s += h.wfu[a][p] * solved[q][a];
                }
                m9[p][q] -= s;
            }
        }
    }
    let mut c = [[0.0; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    c[voigt(i, j)][voigt(k, l)] = m9[i * 3 + j][k * 3 + l];
                }
            }
        }
    }
    // Numerical noise can break the exact symmetries; restore them.
    for i in 0..6 {
        for j in i + 1..6 {
            let v = 0.5 * (c[i][j] + c[j][i]);
            c[i][j] = v;
            c[j][i] = v;
        }
    }
    c
}

/// Extrapolated (infinite-cutoff) Cauchy-Born stiffness.
///
/// For a smooth short-range pair potential the second derivative of the
/// energy under strain has a truncation tail proportional to `rc^-3`: each
/// missing bond contributes `phi'' r^4` and the number of bonds in a shell
/// grows like `r^2`, so the leading omitted term is `∫ r^4 phi'' n(r) dr`
/// with `n(r) ∝ r^2` and `phi'' ∝ r^-8`, i.e. `rc^-3`. Evaluating the
/// stiffness at two cutoffs and eliminating that tail,
/// `C_inf = C2 - A rc2^-3`, `A = (C1 - C2)/(rc1^-3 - rc2^-3)`, gives the
/// converged (infinite-sum) elastic tensor without a mean-field tail whose
/// strain dependence would itself be approximate.
pub fn elastic_tensor_extrapolated(
    lat: &Lattice,
    pot: &dyn PairPotential,
    rc1: f64,
    rc2: f64,
) -> [[f64; 6]; 6] {
    assert!(rc2 > rc1 && rc1 > 0.0);
    let c1 = elastic_tensor(lat, pot, rc1);
    let c2 = elastic_tensor(lat, pot, rc2);
    let f1 = rc1.powi(-3);
    let f2 = rc2.powi(-3);
    let mut c = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let a = (c1[i][j] - c2[i][j]) / (f1 - f2);
            c[i][j] = c2[i][j] - a * f2;
        }
    }
    c
}

/// Frozen (affine) stiffness, no internal relaxation.
pub fn elastic_tensor_frozen(lat: &Lattice, pot: &dyn PairPotential, rmax: f64) -> [[f64; 6]; 6] {
    let n = lat.natoms();
    let u0 = vec![[0.0; 3]; n];
    let m9 = hessians(lat, pot, rmax, &IDENTITY, &u0).wff;
    let mut c = [[0.0; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    c[voigt(i, j)][voigt(k, l)] = m9[i * 3 + j][k * 3 + l];
                }
            }
        }
    }
    for i in 0..6 {
        for j in i + 1..6 {
            let v = 0.5 * (c[i][j] + c[j][i]);
            c[i][j] = v;
            c[j][i] = v;
        }
    }
    c
}

/// Finite-difference stiffness from the relaxed energy under symmetric
/// Voigt strains. This is the independent check of [`elastic_tensor`].
pub fn elastic_tensor_fd(
    lat: &Lattice,
    pot: &dyn PairPotential,
    rmax: f64,
    eps: f64,
) -> [[f64; 6]; 6] {
    let mut c = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in i..6 {
            let fpp = |si: f64, sj: f64| {
                let mut f = add_strain(IDENTITY, si, i);
                f = add_strain(f, sj, j);
                relaxed_energy_density(lat, pot, rmax, &f)
            };
            let v = (fpp(eps, eps) - fpp(eps, -eps) - fpp(-eps, eps) + fpp(-eps, -eps))
                / (4.0 * eps * eps);
            c[i][j] = v;
            c[j][i] = v;
        }
    }
    c
}

fn add_strain(mut f: [[f64; 3]; 3], e: f64, which: usize) -> [[f64; 3]; 3] {
    let d = strain_basis(which);
    for a in 0..3 {
        for b in 0..3 {
            f[a][b] += e * d[a][b];
        }
    }
    f
}

/// Convenience: the three independent cubic constants `(C11, C12, C44)`.
pub fn cubic_constants(c: &[[f64; 6]; 6]) -> (f64, f64, f64) {
    (c[0][0], c[0][1], c[3][3])
}

/// Bulk modulus `K = (C11 + 2 C12)/3` for a cubic crystal.
pub fn bulk_modulus(c: &[[f64; 6]; 6]) -> f64 {
    (c[0][0] + 2.0 * c[0][1]) / 3.0
}

/// Number density of atoms `n_atoms / V0`.
pub fn number_density(lat: &Lattice) -> f64 {
    lat.natoms() as f64 / lat.volume()
}

/// A lattice sum re-export useful in the elastic tests.
pub fn sum_p(lat: &Lattice, p: f64, rmax: f64) -> f64 {
    lattice_sum(lat, p, rmax)
}

/// Dot product re-export (keeps the module self-contained for tests).
pub fn dot3(a: V3, b: V3) -> f64 {
    dot(a, b)
}
