//! Adaptive atomistic/continuum coupling: the quasicontinuum (QC) method.
//!
//! This is the fold/unfold principle of
//! [`analysis/emergence_ladder.md`](../../../../../analysis/emergence_ladder.md)
//! applied to a solid. The fine model is a crystal of atoms interacting
//! through a central pair potential; the coarse model is Cauchy-Born
//! elasticity on a linear finite-element mesh; the indicator is the local
//! variation of the deformation gradient, and where it is large the
//! continuum element is *unfolded* back to atoms.
//!
//! Two couplings are implemented, so the ghost-force problem can be
//! demonstrated and fixed:
//!
//! - [`Coupling::Naive`] - the textbook energy-based QC: the atomistic
//!   region contributes its exact site energies and the continuum region
//!   its Cauchy-Born element energies. The interface atom is counted by
//!   both, and a uniform-strain patch test leaves a **ghost force** at the
//!   interface.
//! - [`Coupling::Qnl`] - the **quasi-nonlocal** energy of Shimokawa,
//!   Mortensen, Schiotz & Jacobsen, *Phys. Rev. B* **69**, 214102 (2004).
//!   Each bond is summed once at the positions of the interpolated
//!   displacement field: a bond with both ends in the atomistic region uses
//!   the exact repatom displacements, every other bond uses the local
//!   (Cauchy-Born) length obtained from the interpolation, which for a
//!   linear element is exactly `|F R|`. Under a uniform strain the two
//!   agree, the energy collapses to the exact atomistic energy, and the
//!   patch test passes to round-off. It was introduced in 1D by Shimokawa
//!   et al. and analysed by E, Ming & Yang; the review is Miller & Tadmor,
//!   *Modelling Simul. Mater. Sci. Eng.* **17**, 053001 (2009).
//!
//! The energy is a function of a reduced set of **representative atoms**
//! (repatoms): every atom in the atomistic region is a degree of freedom,
//! the continuum carries only the finite-element nodes, and every other
//! atom is obtained by linear (simplex) interpolation. The indicator
//! [`Qc::strain_jump`] returns, for each continuum element, the jump of its
//! deformation gradient to its neighbours, and [`Qc::cauchy_born_error`]
//! the local atomistic-minus-Cauchy-Born energy difference;
//! [`Qc::unfold_cells`] turns the worst continuum cells into atomistic ones.
//! The adaptive driver is a solve/indicate/unfold loop around those calls
//! (see `tests/quasicontinuum.rs`).
//!
//! All lengths are the reference lengths of the input lattice; energies are
//! in the potential's units. See `tests/quasicontinuum.rs` for the patch
//! tests and the benchmarks against the fully atomistic solution.

use crate::lattice::{norm, V3};
use crate::potential::PairPotential;

/// How the atomistic and continuum energies are coupled.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Coupling {
    /// Standard energy-based QC; interface atoms are counted twice and ghost
    /// forces appear.
    Naive,
    /// Quasi-nonlocal coupling; every bond is counted once and the patch test
    /// is exact.
    Qnl,
}

/// One undirected bond of the reference crystal. `r` is `X_j - X_i`.
#[derive(Clone, Copy, Debug)]
pub struct BondPair {
    pub i: usize,
    pub j: usize,
    pub r: V3,
    pub radius: f64,
}

/// A linear simplex element of the coarse mesh. In 2D it has three nodes,
/// in 1D two (the extra slot is unused).
#[derive(Clone, Debug)]
pub struct Element {
    pub nodes: Vec<usize>,
    /// Cartesian shape-function gradients `dN_a/dx` (only the first `dim`
    /// components are used).
    pub grads: Vec<V3>,
    /// Number of atoms represented by this element (used by the naive CB
    /// element energy).
    pub natoms: usize,
    /// True if the element lies in the continuum region (naive coupling).
    pub continuum: bool,
}

/// A 1D or 2D atomistic/continuum model.
#[derive(Clone, Debug)]
pub struct Qc {
    pub dim: usize,
    /// Reference positions of all atoms.
    pub atoms: Vec<V3>,
    /// Undirected bond list (each bond once).
    pub bonds: Vec<BondPair>,
    /// Atom indices that are degrees of freedom (repatoms).
    pub repatoms: Vec<usize>,
    /// Atom index -> degree-of-freedom index.
    pub node_of: Vec<Option<usize>>,
    /// Atom index -> true if it belongs to the atomistic region.
    pub atomistic: Vec<bool>,
    /// Atom index -> index of the coarse element that contains it.
    pub elem_of: Vec<usize>,
    pub elements: Vec<Element>,
    /// For each element, its coarse cell `(ci, cj)` (1D uses `(k, 0)`), so
    /// the adaptive driver can unfold a whole cell.
    pub element_cell: Vec<(usize, usize)>,
    /// Atom index -> (dof, coefficient) of the interpolation `u_i = sum`.
    pub interp: Vec<Vec<(usize, f64)>>,
    /// A representative bulk bond list, used by the naive Cauchy-Born
    /// element energy.
    pub cb_bonds: Vec<V3>,
    pub coupling: Coupling,
}

fn dphi(pot: &dyn PairPotential, r: f64) -> f64 {
    pot.d1(r)
}

impl Qc {
    pub fn ndof(&self) -> usize {
        self.repatoms.len()
    }

    /// Interpolate the full displacement field from the repatom values.
    pub fn atom_displacements(&self, udof: &[V3]) -> Vec<V3> {
        let mut u = vec![[0.0; 3]; self.atoms.len()];
        for i in 0..self.atoms.len() {
            for &(dof, c) in &self.interp[i] {
                for d in 0..self.dim {
                    u[i][d] += c * udof[dof][d];
                }
            }
        }
        u
    }

    /// Deformation gradient of an element: `F = I + sum_a u_a (x) grad_a`.
    pub fn element_f(&self, e: usize, udof: &[V3]) -> [[f64; 3]; 3] {
        let mut f = [[0.0; 3]; 3];
        for d in 0..3 {
            f[d][d] = 1.0;
        }
        let el = &self.elements[e];
        for (a, &node) in el.nodes.iter().enumerate() {
            let dof = self.node_of[node].expect("element node must be a repatom");
            let u = udof[dof];
            for x in 0..self.dim {
                for y in 0..self.dim {
                    f[x][y] += u[x] * el.grads[a][y];
                }
            }
        }
        f
    }

    fn apply_f(f: &[[f64; 3]; 3], r: V3, dim: usize) -> V3 {
        let mut q = [0.0; 3];
        for x in 0..dim {
            for y in 0..dim {
                q[x] += f[x][y] * r[y];
            }
        }
        q
    }

    /// Per-atom Cauchy-Born site energy `1/2 sum_bonds phi(|F R|)` at the
    /// given deformation gradient, using the representative bulk bond list.
    fn cb_site_energy(&self, pot: &dyn PairPotential, f: &[[f64; 3]; 3]) -> f64 {
        let mut e = 0.0;
        for &r in &self.cb_bonds {
            let q = Self::apply_f(f, r, self.dim);
            e += pot.energy(norm(q));
        }
        0.5 * e
    }

    /// Total energy (either coupling) at the repatom displacements `udof`.
    pub fn energy(&self, pot: &dyn PairPotential, udof: &[V3]) -> f64 {
        let u = self.atom_displacements(udof);
        match self.coupling {
            Coupling::Qnl => self.energy_qnl(pot, &u, udof),
            Coupling::Naive => self.energy_naive(pot, &u, udof),
        }
    }

    fn energy_qnl(&self, pot: &dyn PairPotential, u: &[V3], _udof: &[V3]) -> f64 {
        // Quasi-nonlocal local rule: every bond is evaluated at the positions
        // obtained from the displacement field of the interpolated
        // (continuum) atoms; for repatoms this is the exact displacement and
        // for a linear element the local bond is the Cauchy-Born bond length
        // `|F R|`. Under a uniform strain it collapses to the full atomistic
        // energy, which is why the patch test is exact.
        let mut e = 0.0;
        for b in &self.bonds {
            let s = [
                b.r[0] + u[b.j][0] - u[b.i][0],
                b.r[1] + u[b.j][1] - u[b.i][1],
                b.r[2] + u[b.j][2] - u[b.i][2],
            ];
            e += pot.energy(norm(s));
        }
        e
    }

    fn energy_naive(&self, pot: &dyn PairPotential, u: &[V3], udof: &[V3]) -> f64 {
        // Atomistic site energies for the atomistic atoms (cross bonds
        // included via the interpolated positions).
        let mut e = 0.0;
        for b in &self.bonds {
            if self.atomistic[b.i] {
                let s = [
                    b.r[0] + u[b.j][0] - u[b.i][0],
                    b.r[1] + u[b.j][1] - u[b.i][1],
                    b.r[2] + u[b.j][2] - u[b.i][2],
                ];
                e += 0.5 * pot.energy(norm(s));
            }
            if self.atomistic[b.j] {
                let s = [
                    b.r[0] + u[b.j][0] - u[b.i][0],
                    b.r[1] + u[b.j][1] - u[b.i][1],
                    b.r[2] + u[b.j][2] - u[b.i][2],
                ];
                e += 0.5 * pot.energy(norm(s));
            }
        }
        // Cauchy-Born element energies for the continuum elements.
        for (ei, el) in self.elements.iter().enumerate() {
            if !el.continuum {
                continue;
            }
            let f = self.element_f(ei, udof);
            e += el.natoms as f64 * self.cb_site_energy(pot, &f);
        }
        e
    }

    /// Force (negative energy gradient) on every repatom.
    pub fn forces(&self, pot: &dyn PairPotential, udof: &[V3]) -> Vec<V3> {
        match self.coupling {
            Coupling::Qnl => self.forces_qnl(pot, udof),
            Coupling::Naive => self.forces_fd(pot, udof),
        }
    }

    fn forces_qnl(&self, pot: &dyn PairPotential, udof: &[V3]) -> Vec<V3> {
        let u = self.atom_displacements(udof);
        let mut grad = vec![[0.0; 3]; self.ndof()];
        let mut atom_grad = vec![[0.0; 3]; self.atoms.len()];
        for b in &self.bonds {
            let s = [
                b.r[0] + u[b.j][0] - u[b.i][0],
                b.r[1] + u[b.j][1] - u[b.i][1],
                b.r[2] + u[b.j][2] - u[b.i][2],
            ];
            let rad = norm(s);
            let g = dphi(pot, rad);
            for d in 0..self.dim {
                let gg = g * s[d] / rad;
                atom_grad[b.i][d] -= gg;
                atom_grad[b.j][d] += gg;
            }
        }
        // Distribute atom-level gradient through the interpolation.
        for i in 0..self.atoms.len() {
            for &(dof, c) in &self.interp[i] {
                for d in 0..self.dim {
                    grad[dof][d] += c * atom_grad[i][d];
                }
            }
        }
        grad.iter()
            .map(|g| [-g[0], -g[1], -g[2]])
            .collect()
    }

    /// Central finite-difference forces; used for the naive coupling, which
    /// is only needed for the patch-test demonstration.
    fn forces_fd(&self, pot: &dyn PairPotential, udof: &[V3]) -> Vec<V3> {
        let h = 1e-6;
        let mut out = vec![[0.0; 3]; self.ndof()];
        let mut up = udof.to_vec();
        for r in 0..self.ndof() {
            for d in 0..self.dim {
                let old = up[r][d];
                up[r][d] = old + h;
                let ep = self.energy(pot, &up);
                up[r][d] = old - h;
                let em = self.energy(pot, &up);
                up[r][d] = old;
                out[r][d] = -(ep - em) / (2.0 * h);
            }
        }
        out
    }

    /// Minimise the energy by Polak-Ribiere conjugate gradients with an
    /// Armijo backtracking line search. `free[r]` masks the repatoms that
    /// may move; `load[r]` is an external force on each repatom.
    pub fn minimize(
        &self,
        pot: &dyn PairPotential,
        free: &[bool],
        load: &[V3],
        max_steps: usize,
        ftol: f64,
    ) -> (Vec<V3>, f64) {
        let n = self.ndof();
        assert_eq!(free.len(), n);
        let mut u = vec![[0.0; 3]; n];
        // Total potential = energy - load . u; force = -d(total)/du.
        let total = |qc: &Qc, u: &[V3]| -> f64 {
            let mut e = qc.energy(pot, u);
            for r in 0..n {
                if free[r] {
                    for d in 0..self.dim {
                        e -= load[r][d] * u[r][d];
                    }
                }
            }
            e
        };
        let force = |qc: &Qc, u: &[V3]| -> Vec<V3> {
            let mut f = qc.forces(pot, u);
            for r in 0..n {
                if free[r] {
                    for d in 0..self.dim {
                        f[r][d] += load[r][d];
                    }
                } else {
                    f[r] = [0.0; 3];
                }
            }
            f
        };
        let dot = |a: &[V3], b: &[V3]| -> f64 {
            let mut s = 0.0;
            for r in 0..n {
                if free[r] {
                    for d in 0..self.dim {
                        s += a[r][d] * b[r][d];
                    }
                }
            }
            s
        };
        let mut e = total(self, &u);
        let f = force(self, &u);
        let mut dir: Vec<V3> = f.clone();
        let mut g = f.clone();
        let mut alpha = 0.1_f64;
        let mut fmax = f
            .iter()
            .enumerate()
            .filter(|(r, _)| free[*r])
            .flat_map(|(_, v)| (0..self.dim).map(move |d| v[d].abs()))
            .fold(0.0, f64::max);
        let mut steps = 0;
        while steps < max_steps && fmax > ftol {
            let gd = dot(&g, &dir);
            if gd >= 0.0 {
                dir = g.clone();
            }
            let gd = dot(&g, &dir);
            if gd.abs() < 1e-300 {
                break;
            }
            let u0 = u.clone();
            let mut a = alpha.max(1e-12);
            let mut ok = false;
            for _ in 0..80 {
                let mut ut = u0.clone();
                for r in 0..n {
                    if free[r] {
                        for d in 0..self.dim {
                            ut[r][d] += a * dir[r][d];
                        }
                    }
                }
                let et = total(self, &ut);
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
            let f_new = force(self, &u);
            let num = dot(&f_new, &f_new) - dot(&f_new, &g);
            let den = dot(&g, &g);
            let beta = if den > 0.0 { (num / den).max(0.0) } else { 0.0 };
            for r in 0..n {
                if free[r] {
                    for d in 0..self.dim {
                        dir[r][d] = f_new[r][d] + beta * dir[r][d];
                    }
                } else {
                    dir[r] = [0.0; 3];
                }
            }
            g = f_new.clone();
            fmax = f_new
                .iter()
                .enumerate()
                .filter(|(r, _)| free[*r])
                .flat_map(|(_, v)| (0..self.dim).map(move |d| v[d].abs()))
                .fold(0.0, f64::max);
        }
        (u, e)
    }

    /// Repatom displacements for a uniform deformation gradient `F`:
    /// `u_r = (F - I) X_r`.
    pub fn uniform_displacements(&self, f: &[[f64; 3]; 3]) -> Vec<V3> {
        let mut u = vec![[0.0; 3]; self.ndof()];
        for r in 0..self.ndof() {
            let x = self.atoms[self.repatoms[r]];
            for d in 0..self.dim {
                for y in 0..self.dim {
                    u[r][d] += (f[d][y] - if d == y { 1.0 } else { 0.0 }) * x[y];
                }
            }
        }
        u
    }

    /// Deformation-gradient jump indicator: for each continuum element,
    /// `max over neighbours |F_e - F_neighbour|` (Frobenius norm of the
    /// difference). Large values mark where the continuum is straining
    /// non-uniformly and the element should unfold to atoms.
    pub fn strain_jump(&self, udof: &[V3]) -> Vec<f64> {
        let fs: Vec<[[f64; 3]; 3]> = (0..self.elements.len())
            .map(|e| self.element_f(e, udof))
            .collect();
        let mut ind = vec![0.0; self.elements.len()];
        // Two elements are neighbours if they share two nodes.
        for a in 0..self.elements.len() {
            for b in (a + 1)..self.elements.len() {
                let shared = self.elements[a]
                    .nodes
                    .iter()
                    .filter(|n| self.elements[b].nodes.contains(n))
                    .count();
                if shared < 2 {
                    continue;
                }
                let mut s = 0.0;
                for x in 0..self.dim {
                    for y in 0..self.dim {
                        s += (fs[a][x][y] - fs[b][x][y]).powi(2);
                    }
                }
                let v = s.sqrt();
                if v > ind[a] {
                    ind[a] = v;
                }
                if v > ind[b] {
                    ind[b] = v;
                }
            }
        }
        ind
    }
    /// Cauchy-Born error indicator: for each element, the total absolute
    /// difference between the per-atom atomistic site energy (at the
    /// interpolated positions) and the per-atom Cauchy-Born site energy (at
    /// the element deformation gradient). It vanishes under a uniform strain
    /// and grows with the non-affine part of the deformation, so it directly
    /// estimates the local continuum energy error.
    pub fn cauchy_born_error(&self, pot: &dyn PairPotential, udof: &[V3]) -> Vec<f64> {
        let u = self.atom_displacements(udof);
        let fs: Vec<[[f64; 3]; 3]> = (0..self.elements.len())
            .map(|e| self.element_f(e, udof))
            .collect();
        let mut err = vec![0.0; self.elements.len()];
        for i in 0..self.atoms.len() {
            let e = self.elem_of[i];
            let f = fs[e];
            let mut e_at = 0.0;
            let mut e_cb = 0.0;
            for b in &self.bonds {
                if b.i == i || b.j == i {
                    let j = if b.i == i { b.j } else { b.i };
                    let s = [
                        b.r[0] + u[j][0] - u[i][0],
                        b.r[1] + u[j][1] - u[i][1],
                        b.r[2] + u[j][2] - u[i][2],
                    ];
                    e_at += 0.5 * pot.energy(norm(s));
                    let q = Self::apply_f(&f, b.r, self.dim);
                    e_cb += 0.5 * pot.energy(norm(q));
                }
            }
            err[e] += (e_at - e_cb).abs();
        }
        err
    }

    /// Adaptive unfold: return the coarse cells of the continuum elements
    /// with the largest strain-jump indicator above `tol`, at most
    /// `max_new` of them (largest first). Passing the result to the next
    /// `build_triangular(..., unfold, ...)` makes those cells fully
    /// atomistic.
    pub fn unfold_cells(&self, ind: &[f64], tol: f64, max_new: usize) -> Vec<(usize, usize)> {
        let mut mark: Vec<(f64, (usize, usize))> = Vec::new();
        for e in 0..self.elements.len() {
            if self.elements[e].continuum && ind[e] > tol {
                mark.push((ind[e], self.element_cell[e]));
            }
        }
        mark.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        let mut out: Vec<(usize, usize)> = Vec::new();
        for (_, cell) in mark {
            if !out.contains(&cell) {
                out.push(cell);
            }
            if out.len() >= max_new {
                break;
            }
        }
        out
    }
}

// ---------------------------------------------------------------------------
// 1D chain
// ---------------------------------------------------------------------------

/// Build a 1D chain. `n` atoms at `x = i a0`. The atomistic set is the set
/// of atoms within `r_atom` of `x_center`; every `coarsen`-th atom outside
/// is a continuum repatom. A vacancy is introduced by dropping the atom
/// nearest `x_center` when `vacancy` is true.
pub fn build_chain(
    n: usize,
    a0: f64,
    rmax: f64,
    x_center: f64,
    r_atom: f64,
    coarsen: usize,
    vacancy: bool,
    coupling: Coupling,
) -> Qc {
    assert!(n >= 5 && coarsen >= 1);
    let mut atoms: Vec<V3> = Vec::new();
    let mut atomistic: Vec<bool> = Vec::new();
    let mut drop_idx: Option<usize> = None;
    if vacancy {
        let c = ((x_center / a0).round() as usize).min(n - 1);
        drop_idx = Some(c);
    }
    for i in 0..n {
        if Some(i) == drop_idx {
            continue;
        }
        let x = i as f64 * a0;
        atoms.push([x, 0.0, 0.0]);
        atomistic.push((x - x_center).abs() <= r_atom);
    }
    let na = atoms.len();

    // Bond list (undirected) with a spherical cutoff.
    let mut bonds = Vec::new();
    for i in 0..na {
        for j in (i + 1)..na {
            let dx = atoms[j][0] - atoms[i][0];
            let r = [dx, 0.0, 0.0];
            let rad = dx.abs();
            if rad > 1e-9 && rad <= rmax {
                bonds.push(BondPair {
                    i,
                    j,
                    r,
                    radius: rad,
                });
            }
        }
    }

    // Repatoms: every atomistic atom, plus coarse nodes outside. A coarse
    // node is an atom whose index in the *original* (pre-vacancy) chain is a
    // multiple of `coarsen`.
    let mut repatoms: Vec<usize> = Vec::new();
    let mut node_of = vec![None; na];
    for i in 0..na {
        if atomistic[i] {
            node_of[i] = Some(repatoms.len());
            repatoms.push(i);
        }
    }
    for i in 0..na {
        if atomistic[i] {
            continue;
        }
        // Recover the original index (vacancy shifts by one after the drop).
        let orig = match drop_idx {
            Some(c) if i >= c => i + 1,
            _ => i,
        };
        if orig % coarsen == 0 || i == 0 || i == na - 1 {
            if node_of[i].is_none() {
                node_of[i] = Some(repatoms.len());
                repatoms.push(i);
            }
        }
    }
    repatoms.sort_by(|&a, &b| atoms[a][0].partial_cmp(&atoms[b][0]).unwrap());
    // Rebuild node_of after sorting.
    node_of = vec![None; na];
    for (dof, &i) in repatoms.iter().enumerate() {
        node_of[i] = Some(dof);
    }

    // Coarse elements: consecutive repatoms.
    let mut elements: Vec<Element> = Vec::new();
    let mut element_cell: Vec<(usize, usize)> = Vec::new();
    for k in 0..repatoms.len().saturating_sub(1) {
        let p = repatoms[k];
        let q = repatoms[k + 1];
        let xp = atoms[p][0];
        let xq = atoms[q][0];
        let len = xq - xp;
        elements.push(Element {
            nodes: vec![p, q],
            grads: vec![[-1.0 / len, 0.0, 0.0], [1.0 / len, 0.0, 0.0]],
            natoms: 0, // filled below
            continuum: false,
        });
        element_cell.push((k, 0));
    }

    // Assign every atom to the element whose left node is the greatest node
    // at or below it.
    let mut elem_of = vec![0usize; na];
    let mut node_pos = vec![0usize; na];
    for (k, &i) in repatoms.iter().enumerate() {
        node_pos[i] = k;
    }
    for i in 0..na {
        // find greatest k with repatoms[k] <= i
        let mut best = 0;
        for (k, &r) in repatoms.iter().enumerate() {
            if r <= i {
                best = k;
            }
        }
        elem_of[i] = best.min(elements.len().saturating_sub(1));
    }

    // Mark continuum elements and count atoms per element.
    // The last atomistic atom's index:
    let last_a = (0..na).rev().find(|&i| atomistic[i]);
    for (ei, el) in elements.iter_mut().enumerate() {
        let p = el.nodes[0];
        el.continuum = match last_a {
            Some(la) => p >= la,
            None => true,
        };
        let q = el.nodes[1];
        let mut cnt = 0;
        for i in 0..na {
            if atoms[i][0] + 1e-12 >= atoms[p][0] && atoms[i][0] <= atoms[q][0] + 1e-12 {
                cnt += 1;
            }
        }
        el.natoms = cnt;
        let _ = ei;
    }

    // Interpolation: repatoms are identity, others linear in their element.
    let mut interp = vec![Vec::new(); na];
    for i in 0..na {
        if let Some(dof) = node_of[i] {
            interp[i].push((dof, 1.0));
        } else {
            let e = elem_of[i];
            let p = elements[e].nodes[0];
            let q = elements[e].nodes[1];
            let xp = atoms[p][0];
            let xq = atoms[q][0];
            let t = (atoms[i][0] - xp) / (xq - xp);
            interp[i].push((node_of[p].unwrap(), 1.0 - t));
            interp[i].push((node_of[q].unwrap(), t));
        }
    }

    // Representative bulk bond list (an interior atom).
    let mid = na / 2;
    let cb_bonds: Vec<V3> = bonds
        .iter()
        .filter_map(|b| {
            if b.i == mid {
                Some(b.r)
            } else if b.j == mid {
                Some([-b.r[0], -b.r[1], -b.r[2]])
            } else {
                None
            }
        })
        .collect();

    Qc {
        dim: 1,
        atoms,
        bonds,
        repatoms,
        node_of,
        atomistic,
        elem_of,
        elements,
        element_cell,
        interp,
        cb_bonds,
        coupling,
    }
}

// ---------------------------------------------------------------------------
// 2D triangular lattice
// ---------------------------------------------------------------------------

/// Build a 2D triangular-lattice sheet. Atoms sit at
/// `x = a0 (i + j/2), y = a0 (sqrt(3)/2) j` for `i in 0..n1`, `j in 0..n2`.
/// The coarse mesh has nodes every `coarsen` lattice steps; an atom is in the
/// atomistic region if its coarse cell lies within `r_atom` of `center`.
/// `unfold` lists extra coarse cells `(ci, cj)` to unfold (adaptivity).
pub fn build_triangular(
    n1: usize,
    n2: usize,
    a0: f64,
    rmax: f64,
    center: V3,
    r_atom: f64,
    coarsen: usize,
    vacancy: bool,
    unfold: &[(usize, usize)],
    coupling: Coupling,
) -> Qc {
    assert!(coarsen >= 1 && n1 > 2 && n2 > 2);
    assert!(
        (n1 - 1) % coarsen == 0 && (n2 - 1) % coarsen == 0,
        "domain must be exactly covered by whole coarse cells"
    );
    let h = (3.0_f64).sqrt() / 2.0;
    let nc1 = (n1 - 1) / coarsen;
    let nc2 = (n2 - 1) / coarsen;
    let idx = |i: usize, j: usize| i + j * n1;
    let pos = |i: usize, j: usize| -> V3 {
        [a0 * (i as f64 + 0.5 * j as f64), a0 * h * j as f64, 0.0]
    };
    let mut remove: Option<(usize, usize)> = None;
    if vacancy {
        // Choose a non-vertex interior atom so the coarse mesh stays whole.
        let mut ri = n1 / 2;
        let rj = n2 / 2;
        if ri % coarsen == 0 && rj % coarsen == 0 {
            ri += 1;
        }
        remove = Some((ri, rj));
    }
    let mut atoms = Vec::new();
    let mut cell_of = Vec::new(); // coarse cell (ci,cj) for each atom
    let mut ij = Vec::new();
    for j in 0..n2 {
        for i in 0..n1 {
            if remove == Some((i, j)) {
                continue;
            }
            atoms.push(pos(i, j));
            cell_of.push(((i / coarsen).min(nc1 - 1), (j / coarsen).min(nc2 - 1)));
            ij.push((i, j));
        }
    }
    let na = atoms.len();

    // Atomistic cells: within r_atom of the centre, plus the requested unfolds.
    let mut atom_cell = vec![false; nc1 * nc2];
    for cj in 0..nc2 {
        for ci in 0..nc1 {
            let px = a0 * (ci as f64 * coarsen as f64 + 0.5 * (cj as f64 * coarsen as f64 + 0.5));
            let py = a0 * h * (cj as f64 * coarsen as f64 + 0.5);
            let d = ((px - center[0]).powi(2) + (py - center[1]).powi(2)).sqrt();
            if d <= r_atom {
                atom_cell[ci + cj * nc1] = true;
            }
        }
    }
    for &(ci, cj) in unfold {
        if ci < nc1 && cj < nc2 {
            atom_cell[ci + cj * nc1] = true;
        }
    }
    let mut atomistic = vec![false; na];
    for (a, &(ci, cj)) in cell_of.iter().enumerate() {
        atomistic[a] = atom_cell[ci + cj * nc1];
    }

    // Bond list (undirected).
    let mut bonds = Vec::new();
    for i in 0..na {
        for j in (i + 1)..na {
            let r = [
                atoms[j][0] - atoms[i][0],
                atoms[j][1] - atoms[i][1],
                0.0,
            ];
            let rad = norm(r);
            if rad > 1e-9 && rad <= rmax {
                bonds.push(BondPair { i, j, r, radius: rad });
            }
        }
    }

    // Repatoms: atomistic atoms plus coarse vertices.
    let mut repatoms = Vec::new();
    let mut node_of = vec![None; na];
    for a in 0..na {
        let (i, j) = ij[a];
        let boundary = i == 0 || i == n1 - 1 || j == 0 || j == n2 - 1;
        if atomistic[a] || boundary || (i % coarsen == 0 && j % coarsen == 0) {
            node_of[a] = Some(repatoms.len());
            repatoms.push(a);
        }
    }

    // Coarse elements per cell: two triangles. Assumes the domain is exactly
    // covered by whole coarse cells (n1-1 and n2-1 are multiples of coarsen).
    let mut elements: Vec<Element> = Vec::new();
    let mut elem_of = vec![0usize; na];
    let mut element_cell: Vec<(usize, usize)> = Vec::new();
    for cj in 0..nc2 {
        for ci in 0..nc1 {
            let i0 = ci * coarsen;
            let j0 = cj * coarsen;
            let i1 = i0 + coarsen;
            let j1 = j0 + coarsen;
            let n00 = idx(i0, j0);
            let n10 = idx(i1, j0);
            let n01 = idx(i0, j1);
            let n11 = idx(i1, j1);
            // T1 = (00,10,01) with N_00 = 1-(li+lj)/m, N_10 = li/m, N_01 = lj/m.
            // T2 = (10,11,01) with N_10 = 1-lj/m, N_11 = li/m+lj/m-1, N_01 = 1-li/m.
            let t1 = make_tri_elem(n00, n10, n01, a0, h, coarsen as f64, false);
            let t2 = make_tri_elem(n10, n11, n01, a0, h, coarsen as f64, true);
            let e1 = elements.len();
            elements.push(t1);
            element_cell.push((ci, cj));
            let e2 = elements.len();
            elements.push(t2);
            element_cell.push((ci, cj));
            for j in j0..=j1 {
                for i in i0..=i1 {
                    if remove == Some((i, j)) {
                        continue;
                    }
                    let li = i - i0;
                    let lj = j - j0;
                    let a = idx(i, j);
                    if li + lj <= coarsen {
                        elem_of[a] = e1;
                    } else {
                        elem_of[a] = e2;
                    }
                }
            }
        }
    }

    // Mark continuum elements (naive) and count atoms.
    for el in elements.iter_mut() {
        el.continuum = !el.nodes.iter().all(|&n| atomistic[n]);
        el.natoms = 0;
    }
    for a in 0..na {
        elements[elem_of[a]].natoms += 1;
    }

    // Interpolation.
    let mut interp = vec![Vec::new(); na];
    for a in 0..na {
        if let Some(dof) = node_of[a] {
            interp[a].push((dof, 1.0));
        } else {
            let e = elem_of[a];
            let el = &elements[e];
            // Barycentric weights from the shape gradients: reconstruct the
            // affine map. Easier: solve in local integer coordinates.
            let (i, j) = ij[a];
            let (ci, cj) = cell_of[a];
            let i0 = ci * coarsen;
            let j0 = cj * coarsen;
            let li = (i - i0) as f64;
            let lj = (j - j0) as f64;
            let m = coarsen as f64;
            // Node atoms of the element in order.
            let nn = &el.nodes;
            // Determine which triangle by the cell-local integer sum.
            let is_t2 = li + lj > coarsen as f64;
            // T1 weights: N(n00)=1-(li+lj)/m, N(n10)=li/m, N(n01)=lj/m.
            // T2 nodes are (n10,n11,n01) with weights
            //   N(n10)=1-lj/m, N(n11)=li/m+lj/m-1, N(n01)=1-li/m.
            let w: [f64; 3] = if is_t2 {
                [1.0 - lj / m, li / m + lj / m - 1.0, 1.0 - li / m]
            } else {
                [1.0 - (li + lj) / m, li / m, lj / m]
            };
            for (k, &nd) in nn.iter().enumerate() {
                if w[k].abs() > 1e-14 {
                    interp[a].push((node_of[nd].unwrap(), w[k]));
                }
            }
        }
    }

    // Representative bulk bond list from an interior atom.
    let mut cb_bonds: Vec<V3> = Vec::new();
    let mid = na / 2;
    for b in &bonds {
        if b.i == mid {
            cb_bonds.push(b.r);
        } else if b.j == mid {
            cb_bonds.push([-b.r[0], -b.r[1], -b.r[2]]);
        }
    }

    Qc {
        dim: 2,
        atoms,
        bonds,
        repatoms,
        node_of,
        atomistic,
        elem_of,
        elements,
        element_cell,
        interp,
        cb_bonds,
        coupling,
    }
}

fn make_tri_elem(
    n0: usize,
    n1: usize,
    n2: usize,
    a0: f64,
    h: f64,
    m: f64,
    is_t2: bool,
) -> Element {
    // Physical shape gradients. In integer-reference coordinates with
    // x = a0 (li + lj/2), y = a0 h lj:
    //   dli/dx = 1/a0, dli/dy = -1/(2 a0 h); dlj/dx = 0, dlj/dy = 1/(a0 h).
    let dli = [1.0 / a0, -1.0 / (2.0 * a0 * h), 0.0];
    let dlj = [0.0, 1.0 / (a0 * h), 0.0];
    let mut g = [[0.0; 3]; 3];
    if is_t2 {
        // N_A = 1 - lj/m, N_B = li/m + lj/m - 1, N_C = 1 - li/m.
        for d in 0..2 {
            g[0][d] = -dlj[d] / m;
            g[1][d] = (dli[d] + dlj[d]) / m;
            g[2][d] = -dli[d] / m;
        }
    } else {
        // N_A = 1 - (li+lj)/m, N_B = li/m, N_C = lj/m.
        for d in 0..2 {
            g[0][d] = -(dli[d] + dlj[d]) / m;
            g[1][d] = dli[d] / m;
            g[2][d] = dlj[d] / m;
        }
    }
    Element {
        nodes: vec![n0, n1, n2],
        grads: g.to_vec(),
        natoms: 0,
        continuum: false,
    }
}
