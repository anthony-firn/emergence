//! Point and planar defects on small periodic atomistic supercells.
//!
//! This is the rung immediately above the perfect-crystal Cauchy-Born rung
//! of [`crate::elastic`] and below continuum plasticity
//! ([`crate::plasticity`]). The fine model is a set of atoms interacting
//! through a central pair potential with periodic boundary conditions; the
//! coarse model is the perfect crystal. The defects are:
//!
//! - a **vacancy**, whose formation energy is `E_v = E(N-1) - (N-1) e_coh`
//!   (the removed atom is returned to a reservoir at the cohesive energy);
//! - **free surfaces**, whose energy is the excess over the bulk;
//! - the **intrinsic stacking fault** and the **generalized stacking-fault
//!   energy surface (γ-surface)** of an fcc `{111}` plane, the last
//!   continuum-relevant input for the Peierls-Nabarro model in
//!   [`crate::plasticity`].
//!
//! Minimisation uses Polak-Ribiere conjugate gradients with an Armijo
//! backtracking line search, driven by analytic forces. All sizes are
//! deliberately small (tens to a few hundred atoms) so the tests stay well
//! below 30 s.
//!
//! For a pair potential the **unrelaxed** vacancy formation energy has an
//! exact value. Removing one atom from the perfect crystal removes all its
//! bonds, a per-atom pair sum of `-2 e_coh` (because `e_coh = ½ Σ_j φ`,
//! while the total crystal energy is `½ Σ_i Σ_j φ`), so
//! `E_v^unrelaxed = -e_coh` after referring the remaining `N-1` atoms to the
//! bulk cohesive energy. Relaxation of the neighbours lowers it.
//!
//! The **intrinsic stacking-fault energy** of a pair-potential fcc crystal
//! is, to leading order, the hcp-fcc energy difference per unit area. The
//! hcp stacking `ABAB` has one reversed step per two close-packed planes
//! while fcc has none, so one fault costs half the hcp-fcc difference per
//! plane:
//! `γ_ISF = 2 (e_hcp - e_fcc) / A_{111}` with `A_{111}` the area per atom
//! on a `{111}` plane. This is the standard pair-potential result; because
//! LJ puts hcp *below* fcc by ~0.01 %, the predicted γ_ISF is slightly
//! negative, a well-known deficiency of pair potentials (real fcc metals
//! have positive stacking-fault energies).
//!
//! References: Bitzek et al. 2006 (FIRE, an alternative minimiser); Vitek,
//! *Phil. Mag.* **18**, 773
//! (1968) (γ-surface); Hirth & Lothe, *Theory of Dislocations* (stacking
//! faults); the fcc LJ cohesive energy and lattice sums are checked in
//! `analysis/materials.md`.

use crate::lattice::{cross, dot, norm, scale, Bond, Lattice, V3};
use crate::potential::{energy_per_atom_with_tail, PairPotential};

/// Replicate a lattice `reps = (n0, n1, n2)` times along its primitive
/// vectors, keeping the atoms as the basis of a single larger periodic
/// cell. Fractional coordinates are scaled so the basis is unchanged
/// physically.
pub fn replicate(base: &Lattice, reps: [usize; 3]) -> Lattice {
    let prim = [
        scale(base.prim[0], reps[0] as f64),
        scale(base.prim[1], reps[1] as f64),
        scale(base.prim[2], reps[2] as f64),
    ];
    let nb = base.basis.len();
    let mut basis = Vec::with_capacity(nb * reps[0] * reps[1] * reps[2]);
    for n0 in 0..reps[0] {
        for n1 in 0..reps[1] {
            for n2 in 0..reps[2] {
                for b in &base.basis {
                    basis.push([
                        (n0 as f64 + b[0]) / reps[0] as f64,
                        (n1 as f64 + b[1]) / reps[1] as f64,
                        (n2 as f64 + b[2]) / reps[2] as f64,
                    ]);
                }
            }
        }
    }
    Lattice {
        prim,
        basis,
        label: format!("{}x{}x{}", reps[0], reps[1], reps[2]),
    }
}

/// Raw periodic pair energy (no tail) and Cartesian forces for the
/// displacement field `u` (one three-vector per atom).
///
/// A bond of the fixed neighbour list contributes `φ(|r + u_j - u_i|)`; the
/// force on atom `i` is `Σ_b φ'(|s|) ŝ`. The bond topology is frozen, which
/// is exact for the small displacements of a relaxation.
pub fn energy_forces(
    bonds: &[Vec<Bond>],
    pot: &dyn PairPotential,
    u: &[V3],
) -> (f64, Vec<V3>) {
    let n = u.len();
    let mut e = 0.0;
    let mut f = vec![[0.0; 3]; n];
    for i in 0..n {
        for b in &bonds[i] {
            let s = [
                b.r[0] + u[b.j][0] - u[i][0],
                b.r[1] + u[b.j][1] - u[i][1],
                b.r[2] + u[b.j][2] - u[i][2],
            ];
            let r = norm(s);
            e += pot.energy(r);
            let g = pot.d1(r) / r; // φ'(r)/r
            for d in 0..3 {
                f[i][d] += g * s[d];
            }
        }
    }
    (0.5 * e, f)
}

/// Result of a relaxation.
#[derive(Clone, Debug)]
pub struct Relax {
    pub u: Vec<V3>,
    pub energy: f64,
    pub max_force: f64,
    pub steps: usize,
}

fn max_force(f: &[V3], free: [bool; 3]) -> f64 {
    f.iter()
        .flat_map(|v| (0..3).filter(|&d| free[d]).map(move |d| v[d].abs()))
        .fold(0.0, f64::max)
}

fn dot_v(a: &[V3], b: &[V3], free: [bool; 3]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (0..3).filter(|&d| free[d]).map(|d| x[d] * y[d]).sum::<f64>())
        .sum()
}

/// Conjugate-gradient minimisation of the potential energy (Polak-Ribiere
/// with an Armijo backtracking line search). `free` masks the degrees of
/// freedom that may move (`[true; 3]` is a full relaxation; `[false,
/// false, true]` relaxes only normal to the slip plane, the constrained
/// γ-surface convention of Vitek 1968).
///
/// A conjugate-gradient solver is used rather than raw FIRE because the
/// free slabs have a translational zero mode and FIRE occasionally
/// overshoots it at large cutoff; backtracking makes every step
/// energy-monotone.
pub fn minimize(
    bonds: &[Vec<Bond>],
    pot: &dyn PairPotential,
    mut u: Vec<V3>,
    free: [bool; 3],
    max_steps: usize,
    force_tol: f64,
) -> Relax {
    let n = u.len();
    let (mut e, mut f) = energy_forces(bonds, pot, &u);
    for x in f.iter_mut() {
        for d in 0..3 {
            if !free[d] {
                x[d] = 0.0;
            }
        }
    }
    // gradient g = -F and descent direction d = -g = F
    let mut g: Vec<V3> = f.iter().map(|v| [-v[0], -v[1], -v[2]]).collect();
    let mut dir: Vec<V3> = f.clone();
    let mut alpha = 0.1_f64;
    let mut steps = 0usize;
    let mut fmax = max_force(&f, free);
    while steps < max_steps && fmax > force_tol {
        let gd = dot_v(&g, &dir, free);
        if gd >= 0.0 {
            // restart with steepest descent
            dir = g.iter().map(|v| [-v[0], -v[1], -v[2]]).collect();
        }
        let gd = dot_v(&g, &dir, free);
        if gd.abs() < 1e-300 {
            break;
        }
        let u0 = u.clone();
        let mut a = alpha.max(1e-12);
        let mut ok = false;
        for _ in 0..80 {
            let mut ut = u0.clone();
            for i in 0..n {
                for d in 0..3 {
                    if free[d] {
                        ut[i][d] += a * dir[i][d];
                    }
                }
            }
            let (et, _) = energy_forces(bonds, pot, &ut);
            if et <= e + 1e-4 * a * gd {
                u = ut;
                e = et;
                alpha = a;
                ok = true;
                break;
            }
            a *= 0.5;
        }
        if !ok {
            break;
        }
        steps += 1;
        let (_, f_new) = energy_forces(bonds, pot, &u);
        let mut f_new = f_new;
        for x in f_new.iter_mut() {
            for d in 0..3 {
                if !free[d] {
                    x[d] = 0.0;
                }
            }
        }
        let g_new: Vec<V3> = f_new.iter().map(|v| [-v[0], -v[1], -v[2]]).collect();
        let gg = dot_v(&g, &g, free);
        let num = dot_v(&g_new, &g_new, free) - dot_v(&g_new, &g, free);
        let beta = if gg > 0.0 { (num / gg).max(0.0) } else { 0.0 };
        for i in 0..n {
            for d in 0..3 {
                if free[d] {
                    dir[i][d] = -g_new[i][d] + beta * dir[i][d];
                } else {
                    dir[i][d] = 0.0;
                }
            }
        }
        g = g_new;
        f = f_new;
        fmax = max_force(&f, free);
    }
    Relax {
        u,
        energy: e,
        max_force: fmax,
        steps,
    }
}

/// Exact zero-temperature surface energy of a central-force crystal in the
/// infinite-slab limit (the classical **broken-bond** lattice sum):
///
/// ```text
/// γ(n) = ½ ρ Σ_{bonds, r·n > 0} (r·n) |φ(r)|,
/// ```
///
/// where `ρ` is the number density. The reasoning: cleaving at a plane
/// removes every bond that crosses it, one bond of range `h = r·n` is
/// removed for a slab layer of thickness `h`, so the removed energy per
/// area is `ρ Σ_{r·n>0} h |φ|`; two surfaces share it, giving the `½`.
/// This reproduces the relaxed slab surface energy to within the small
/// relaxation (and exactly for the rigid lattice). Reference: the
/// broken-bond/lattice-sum construction of Lennard-Jones (1924) and
/// Shuttleworth, *Proc. Phys. Soc. A* **62**, 167 (1949).
pub fn broken_bond_surface_energy(
    base: &Lattice,
    pot: &dyn PairPotential,
    rmax: f64,
    normal: V3,
) -> f64 {
    let nrm = norm(normal);
    let n = [normal[0] / nrm, normal[1] / nrm, normal[2] / nrm];
    let rho = base.natoms() as f64 / base.volume();
    let bonds = base.bond_list(rmax);
    let nb = base.natoms();
    let mut total = 0.0;
    for i in 0..nb {
        for b in &bonds[i] {
            let h = dot(b.r, n);
            if h > 0.0 {
                total += h * pot.energy(b.radius).abs();
            }
        }
    }
    0.5 * rho * total / nb as f64
}

/// Raw periodic energy of a relaxed cell, `(E_raw, max_force, steps)`.
fn relax_energy(
    bonds: &[Vec<Bond>],
    pot: &dyn PairPotential,
    u0: Vec<V3>,
    free: [bool; 3],
    ftol: f64,
) -> (f64, f64, usize) {
    let r = minimize(bonds, pot, u0, free, 20_000, ftol);
    (r.energy, r.max_force, r.steps)
}

// ---------------------------------------------------------------------------
// Vacancy
// ---------------------------------------------------------------------------

/// Vacancy formation energy in a periodic supercell.
#[derive(Clone, Debug)]
pub struct Vacancy {
    pub n_atoms: usize,
    /// `Energy(defect) - (N-1) e_coh`, atoms held fixed.
    pub e_unrelaxed: f64,
    /// Same after relaxation.
    pub e_relaxed: f64,
    /// `e_unrelaxed - e_relaxed` (positive).
    pub relaxation_energy: f64,
    pub max_force: f64,
    pub steps: usize,
}

/// Formation energy of one vacancy in `reps` repetitions of `base`. The
/// reference is the *converged* cohesive energy
/// [`energy_per_atom_with_tail`] evaluated on the primitive cell, and the
/// defect total energy carries the same homogeneous tail at the defect
/// density, so the two are consistent.
pub fn vacancy_formation(
    base: &Lattice,
    pot: &dyn PairPotential,
    rmax: f64,
    reps: [usize; 3],
) -> Vacancy {
    let e_coh = energy_per_atom_with_tail(base, pot, rmax);
    let mut defect = replicate(base, reps);
    defect.basis.remove(0); // vacancy at the origin
    let n = defect.natoms();
    let rho = n as f64 / defect.volume();
    let tail = n as f64 * pot.tail_energy_per_atom(rho, rmax);
    let bonds = defect.bond_list(rmax);
    let u0 = vec![[0.0; 3]; n];
    let (e_raw, _) = energy_forces(&bonds, pot, &u0);
    let e_unrelaxed = e_raw + tail - n as f64 * e_coh;
    let rel = minimize(&bonds, pot, u0, [true; 3], 4_000, 1e-4);
    let e_relaxed = rel.energy + tail - n as f64 * e_coh;
    Vacancy {
        n_atoms: n,
        e_unrelaxed,
        e_relaxed,
        relaxation_energy: e_unrelaxed - e_relaxed,
        max_force: rel.max_force,
        steps: rel.steps,
    }
}

// ---------------------------------------------------------------------------
// Free surfaces
// ---------------------------------------------------------------------------

/// `(111)` slab of an fcc crystal with `layers` close-packed planes and a
/// vacuum gap of at least `vacuum` (choose `> rmax`). The in-plane cell is
/// the primitive triangular cell with nearest-neighbour distance `r_nn`,
/// one atom per layer, stacked `ABC…`.
pub fn fcc_111_slab(r_nn: f64, layers: usize, vacuum: f64) -> Lattice {
    let u = [r_nn, 0.0, 0.0];
    let v = [0.5 * r_nn, 0.5 * 3.0_f64.sqrt() * r_nn, 0.0];
    let d = r_nn * (2.0_f64 / 3.0).sqrt(); // {111} interplanar spacing
    let lz = layers as f64 * d + vacuum;
    let prim = [u, v, [0.0, 0.0, lz]];
    let offsets = [[0.0, 0.0], [1.0 / 3.0, 1.0 / 3.0], [2.0 / 3.0, 2.0 / 3.0]];
    let basis = (0..layers)
        .map(|m| {
            let o = offsets[m % 3];
            [o[0], o[1], (m as f64 * d) / lz]
        })
        .collect();
    Lattice {
        prim,
        basis,
        label: format!("fcc111-{layers}"),
    }
}

/// `(100)` slab of an fcc crystal with `layers` planes of spacing `a/2`
/// (two atoms per layer) and a vacuum gap of at least `vacuum`.
pub fn fcc_100_slab(a: f64, layers: usize, vacuum: f64) -> Lattice {
    let dz = a / 2.0;
    let lz = layers as f64 * dz + vacuum;
    let prim = [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, lz]];
    let mut basis = Vec::with_capacity(2 * layers);
    for m in 0..layers {
        let z = (m as f64 * dz) / lz;
        if m % 2 == 0 {
            basis.push([0.0, 0.0, z]);
            basis.push([0.5, 0.5, z]);
        } else {
            basis.push([0.5, 0.0, z]);
            basis.push([0.0, 0.5, z]);
        }
    }
    Lattice {
        prim,
        basis,
        label: format!("fcc100-{layers}"),
    }
}

/// Surface area of the in-plane cell of a slab.
pub fn slab_area(lat: &Lattice) -> f64 {
    norm(cross(lat.prim[0], lat.prim[1]))
}

/// Unrelaxed and relaxed surface energy of a slab:
/// `(E_slab - N e_bulk) / (2 A)`. `e_bulk` must be the per-atom energy of
/// the *same* potential and cutoff (use [`energy_per_atom`] for a
/// same-cutoff comparison, or the tail-corrected value for the converged
/// limit). `e_bulk` is passed in so the caller chooses the reference.
pub fn surface_energy(
    slab: &Lattice,
    pot: &dyn PairPotential,
    rmax: f64,
    e_bulk: f64,
    ftol: f64,
) -> (f64, f64) {
    let n = slab.natoms();
    let bonds = slab.bond_list(rmax);
    let u0 = vec![[0.0; 3]; n];
    let (e0, _) = energy_forces(&bonds, pot, &u0);
    let (e1, _, _) = relax_energy(&bonds, pot, u0.clone(), [true; 3], ftol);
    let area = slab_area(slab);
    ((e0 - n as f64 * e_bulk) / (2.0 * area), (e1 - n as f64 * e_bulk) / (2.0 * area))
}

// ---------------------------------------------------------------------------
// Stacking fault and γ-surface
// ---------------------------------------------------------------------------

/// hcp-fcc stacking-fault estimate for a pair potential:
/// `γ = 2 (e_hcp - e_fcc) / A_{111}`. `a` is the fcc cubic lattice constant.
#[derive(Clone, Copy, Debug)]
pub struct StackingFault {
    pub e_fcc: f64,
    pub e_hcp: f64,
    pub area_per_atom: f64,
    pub gamma: f64,
}

/// Evaluate the hcp-fcc relation at the same nearest-neighbour distance.
pub fn stacking_fault_hcp_fcc(a_cubic: f64, pot: &dyn PairPotential, rmax: f64) -> StackingFault {
    let r_nn = a_cubic / 2.0_f64.sqrt();
    let e_fcc = energy_per_atom_with_tail(&Lattice::fcc(a_cubic), pot, rmax);
    let e_hcp = energy_per_atom_with_tail(&Lattice::hcp(r_nn), pot, rmax);
    let area_per_atom = 3.0_f64.sqrt() / 4.0 * a_cubic * a_cubic;
    StackingFault {
        e_fcc,
        e_hcp,
        area_per_atom,
        gamma: 2.0 * (e_hcp - e_fcc) / area_per_atom,
    }
}

/// Displacement field that rigidly shifts every atom in layer `split` and
/// above (one atom per layer, as built by [`fcc_111_slab`]) by `shift`.
pub fn shift_top(slab: &Lattice, split: usize, shift: V3) -> Vec<V3> {
    let mut u = vec![[0.0; 3]; slab.natoms()];
    for x in u.iter_mut().skip(split) {
        *x = shift;
    }
    u
}

/// One point of a γ-surface scan.
#[derive(Clone, Copy, Debug)]
pub struct GammaPoint {
    /// Shift amplitude along the path direction.
    pub t: f64,
    /// Unrelaxed (rigid-shift) energy `(E - E_0)/A`.
    pub gamma_unrelaxed: f64,
    /// Energy after relaxing only the normal coordinate (Vitek 1968).
    pub gamma_relaxed: f64,
}

/// Generalized stacking-fault energy along a rigid in-plane shift
/// `shift = t * dir` of the upper half of a `(111)` slab, for `t` in
/// `0..=1` with `npts` points. `dir` should be the Shockley partial
/// `(u+v)/3` for the `⟨112⟩` path or a nearest-neighbour vector `u` for the
/// `⟨110⟩` path. The reference `E_0` is the perfect slab, relaxed for the
/// relaxed column and rigid for the unrelaxed column.
pub fn gamma_surface_111(
    r_nn: f64,
    layers: usize,
    vacuum: f64,
    split: usize,
    rmax: f64,
    pot: &dyn PairPotential,
    dir: V3,
    npts: usize,
    ftol: f64,
) -> Vec<GammaPoint> {
    let slab = fcc_111_slab(r_nn, layers, vacuum);
    let area = slab_area(&slab);
    let bonds = slab.bond_list(rmax);
    let zero = vec![[0.0; 3]; slab.natoms()];
    let (e0_raw, _) = energy_forces(&bonds, pot, &zero);
    let (e0_rel, _, _) = relax_energy(&bonds, pot, zero.clone(), [false, false, true], ftol);
    let mut out = Vec::with_capacity(npts + 1);
    for k in 0..=npts {
        let t = k as f64 / npts as f64;
        let shift = scale(dir, t);
        let u0 = shift_top(&slab, split, shift);
        let (eu, _) = energy_forces(&bonds, pot, &u0);
        let (er, _, _) = relax_energy(&bonds, pot, u0, [false, false, true], ftol);
        out.push(GammaPoint {
            t,
            gamma_unrelaxed: (eu - e0_raw) / area,
            gamma_relaxed: (er - e0_rel) / area,
        });
    }
    out
}

/// Maximum of a γ-surface scan on `(0, 1]` (the unstable stacking-fault
/// energy γ_us). The point at `t = 0` is the perfect crystal and is
/// excluded.
pub fn unstable_stacking_fault(curve: &[GammaPoint], relaxed: bool) -> f64 {
    curve
        .iter()
        .filter(|p| p.t > 1e-12)
        .map(|p| if relaxed { p.gamma_relaxed } else { p.gamma_unrelaxed })
        .fold(f64::NEG_INFINITY, f64::max)
}
