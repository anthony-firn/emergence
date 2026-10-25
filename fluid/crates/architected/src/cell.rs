//! Periodic unit cells and two-scale homogenisation.
//!
//! A lattice of beams is folded to an effective continuum by solving, in one
//! unit cell, the periodic fluctuation that minimises the strain energy under
//! a prescribed macroscopic strain. The effective stiffness follows from the
//! energy:
//!
//! ```text
//!   C_ab = [ U(e_a + e_b) - U(e_a) - U(e_b) ] / V ,
//! ```
//!
//! which is exact for the quadratic energy and needs no boundary tractions.
//! Strains are engineering Voigt: `[e_xx, e_yy, g_xy]` in 2D and
//! `[e_xx, e_yy, e_zz, g_xy, g_yz, g_zx]` in 3D, with `g_ij = 2 e_ij`.
//!
//! Boundary conditions are genuinely periodic. Every member endpoint is
//! reduced modulo the lattice vectors to a canonical node; its translation
//! image is kept, and the affine displacement `eps . x` is applied at the
//! physical position. The fluctuation is the same at a node and all its
//! images, which is periodicity.

use crate::beam::*;
use crate::linalg::{solve, Dense};
use std::f64::consts::PI;

/// Floor of a fractional lattice coordinate, snapping values that are within
/// round-off of an integer. Without this a node exactly on a cell boundary
/// can be assigned to the next cell and the periodic reduction fails.
fn exact_floor(f: f64) -> f64 {
    let r = f.round();
    if (f - r).abs() < 1e-9 {
        r
    } else {
        f.floor()
    }
}

/// Sheet/beam network element kind.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// Rigid joints: axial + bending (Euler-Bernoulli, or Timoshenko if `nu`
    /// is set). Bending-dominated at low density.
    Frame,
    /// Pin joints: axial only. Stretch-dominated.
    Truss,
}

// ---------------------------------------------------------------------------
// 2D
// ---------------------------------------------------------------------------

/// A member before periodic reduction.
#[derive(Clone, Copy, Debug)]
pub struct RawMember2d {
    pub a: [f64; 2],
    pub b: [f64; 2],
    pub e: f64,
    pub area: f64,
    pub inertia: f64,
}

/// A member after periodic reduction: canonical nodes plus lattice images.
#[derive(Clone, Copy, Debug)]
pub struct Member2d {
    pub a: usize,
    pub b: usize,
    pub a_img: [i32; 2],
    pub b_img: [i32; 2],
    pub e: f64,
    pub area: f64,
    pub inertia: f64,
}

fn dedup2d(members: Vec<Member2d>) -> Vec<Member2d> {
    // Two entries are the same physical bond when their canonical node pair
    // and the relative lattice image agree; the key is translation invariant.
    let mut seen = std::collections::HashSet::new();
    let mut out = Vec::new();
    for m in members {
        let key = if m.a <= m.b {
            (m.a, m.b, [m.b_img[0] - m.a_img[0], m.b_img[1] - m.a_img[1]])
        } else {
            (m.b, m.a, [m.a_img[0] - m.b_img[0], m.a_img[1] - m.b_img[1]])
        };
        if seen.insert(key) {
            out.push(m);
        }
    }
    out
}

/// 2D periodic beam lattice.
#[derive(Clone, Debug)]
pub struct Cell2d {
    pub nodes: Vec<[f64; 2]>,
    pub members: Vec<Member2d>,
    pub a1: [f64; 2],
    pub a2: [f64; 2],
    pub kind: Kind,
    /// Poisson ratio for the Timoshenko shear correction; `None` = EB.
    pub nu: Option<f64>,
}

fn reduce2d(
    nodes: &mut Vec<[f64; 2]>,
    a1: [f64; 2],
    a2: [f64; 2],    p: [f64; 2],
) -> (usize, [i32; 2]) {
    // fractional coordinates: p = x + n1 a1 + n2 a2, x in the cell.
    let det = a1[0] * a2[1] - a1[1] * a2[0];
    assert!(det.abs() > 0.0);
    let n1f = (p[0] * a2[1] - p[1] * a2[0]) / det;
    let n2f = (a1[0] * p[1] - a1[1] * p[0]) / det;
    let n1 = exact_floor(n1f);
    let n2 = exact_floor(n2f);
    let x = [
        p[0] - n1 * a1[0] - n2 * a2[0],
        p[1] - n1 * a1[1] - n2 * a2[1],
    ];
    for (i, q) in nodes.iter().enumerate() {
        if (q[0] - x[0]).abs() < 1e-8 && (q[1] - x[1]).abs() < 1e-8 {
            return (i, [n1 as i32, n2 as i32]);
        }
    }
    nodes.push(x);
    (nodes.len() - 1, [n1 as i32, n2 as i32])
}

impl Cell2d {
    pub fn new(a1: [f64; 2], a2: [f64; 2], raw: &[RawMember2d], kind: Kind, nu: Option<f64>) -> Self {
        let mut nodes = Vec::new();
        let mut members = Vec::new();
        for r in raw {
            let (a, a_img) = reduce2d(&mut nodes, a1, a2, r.a);
            let (b, b_img) = reduce2d(&mut nodes, a1, a2, r.b);
            members.push(Member2d {
                a,
                b,
                a_img,
                b_img,
                e: r.e,
                area: r.area,
                inertia: r.inertia,
            });
        }
        Self {
            nodes,
            members: dedup2d(members),
            a1,
            a2,
            kind,
            nu,
        }
    }

    pub fn area(&self) -> f64 {
        (self.a1[0] * self.a2[1] - self.a1[1] * self.a2[0]).abs()
    }

    /// Member area fraction of the unit cell (relative density).
    pub fn relative_density(&self) -> f64 {
        let mut a = 0.0;
        for m in &self.members {
            let xa = self.physical(m.a, m.a_img);
            let xb = self.physical(m.b, m.b_img);
            a += ((xb[0] - xa[0]).powi(2) + (xb[1] - xa[1]).powi(2)).sqrt() * m.area;
        }
        a / self.area()
    }

    pub fn physical(&self, node: usize, img: [i32; 2]) -> [f64; 2] {
        let x = self.nodes[node];
        [
            x[0] + img[0] as f64 * self.a1[0] + img[1] as f64 * self.a2[0],
            x[1] + img[0] as f64 * self.a1[1] + img[1] as f64 * self.a2[1],
        ]
    }

    /// Affine (macroscopic) displacement `eps . x` at a physical point.
    pub fn affine(&self, x: [f64; 2], strain: [f64; 3]) -> [f64; 2] {
        let [exx, eyy, gxy] = strain;
        [
            exx * x[0] + 0.5 * gxy * x[1],
            0.5 * gxy * x[0] + eyy * x[1],
        ]
    }

    fn dof(&self) -> usize {
        match self.kind {
            Kind::Frame => 3,
            Kind::Truss => 2,
        }
    }

    fn element_k(&self, xa: [f64; 2], xb: [f64; 2], m: &Member2d) -> Vec<f64> {
        let dx = xb[0] - xa[0];
        let dy = xb[1] - xa[1];
        let l = (dx * dx + dy * dy).sqrt();
        let (c, s) = (dx / l, dy / l);
        match self.kind {
            Kind::Truss => truss2d_global(l, m.e, m.area, c, s).to_vec(),
            Kind::Frame => {
                let phi = match self.nu {
                    Some(nu) => timoshenko_phi(l, m.e, m.area, m.inertia, nu),
                    None => 0.0,
                };
                let kl = frame2d_local(l, m.e, m.area, m.inertia, phi);
                let kg = frame2d_global(&kl, c, s);
                kg.iter().flat_map(|r| r.iter().copied()).collect()
            }
        }
    }

    fn endpoint_g(&self, x: [f64; 2], strain: [f64; 3]) -> Vec<f64> {
        let u = self.affine(x, strain);
        match self.kind {
            Kind::Truss => vec![u[0], u[1]],
            Kind::Frame => vec![u[0], u[1], 0.0],
        }
    }

    /// Strain energy of the cell under a macroscopic strain.
    pub fn energy(&self, strain: [f64; 3]) -> f64 {
        let dof = self.dof();
        let n = self.nodes.len() * dof;
        let mut k = Dense::zeros(n);
        let mut f = vec![0.0; n];
        for m in &self.members {
            let xa = self.physical(m.a, m.a_img);
            let xb = self.physical(m.b, m.b_img);
            let ke = self.element_k(xa, xb, m);
            let mut idx = Vec::with_capacity(2 * dof);
            for t in 0..dof {
                idx.push(m.a * dof + t);
            }
            for t in 0..dof {
                idx.push(m.b * dof + t);
            }
            k.add_outer_block(&idx, &ke);
            let ga = self.endpoint_g(xa, strain);
            let gb = self.endpoint_g(xb, strain);
            let g: Vec<f64> = ga.iter().chain(gb.iter()).copied().collect();
            for (i, &ri) in idx.iter().enumerate() {
                let mut s = 0.0;
                for j in 0..2 * dof {
                    s += ke[i * 2 * dof + j] * g[j];
                }
                f[ri] -= s;
            }
        }
        // Pin node 0 completely to remove rigid-body translations/rotation.
        let pinned = dof;
        let nred = n - pinned;
        let mut map = vec![usize::MAX; n];
        for (new, old) in (pinned..n).enumerate() {
            map[old] = new;
        }
        let mut kr = Dense::zeros(nred);
        for i in pinned..n {
            for j in pinned..n {
                kr.set(map[i], map[j], k.get(i, j));
            }
        }
        let mut fr: Vec<f64> = (pinned..n).map(|i| f[i]).collect();
        solve(&mut kr, &mut fr).expect("singular unit cell");
        let mut q = vec![0.0; n];
        for i in pinned..n {
            q[i] = fr[map[i]];
        }
        let mut u = 0.0;
        for m in &self.members {
            let xa = self.physical(m.a, m.a_img);
            let xb = self.physical(m.b, m.b_img);
            let ke = self.element_k(xa, xb, m);
            let mut idx = Vec::with_capacity(2 * dof);
            for t in 0..dof {
                idx.push(m.a * dof + t);
            }
            for t in 0..dof {
                idx.push(m.b * dof + t);
            }
            let ga = self.endpoint_g(xa, strain);
            let gb = self.endpoint_g(xb, strain);
            let g: Vec<f64> = ga.iter().chain(gb.iter()).copied().collect();
            let d: Vec<f64> = (0..2 * dof).map(|i| q[idx[i]] + g[i]).collect();
            let mut e = 0.0;
            for i in 0..2 * dof {
                for j in 0..2 * dof {
                    e += d[i] * ke[i * 2 * dof + j] * d[j];
                }
            }
            u += 0.5 * e;
        }
        u
    }

    /// Fluctuation (periodic part) nodal displacement for a macroscopic
    /// strain; the pinned node's entries are zero.
    pub fn fluctuation(&self, strain: [f64; 3]) -> Vec<f64> {
        let dof = self.dof();
        let n = self.nodes.len() * dof;
        let mut k = Dense::zeros(n);
        let mut f = vec![0.0; n];
        for m in &self.members {
            let xa = self.physical(m.a, m.a_img);
            let xb = self.physical(m.b, m.b_img);
            let ke = self.element_k(xa, xb, m);
            let mut idx = Vec::with_capacity(2 * dof);
            for t in 0..dof {
                idx.push(m.a * dof + t);
            }
            for t in 0..dof {
                idx.push(m.b * dof + t);
            }
            k.add_outer_block(&idx, &ke);
            let ga = self.endpoint_g(xa, strain);
            let gb = self.endpoint_g(xb, strain);
            let g: Vec<f64> = ga.iter().chain(gb.iter()).copied().collect();
            for (i, &ri) in idx.iter().enumerate() {
                let mut s = 0.0;
                for j in 0..2 * dof {
                    s += ke[i * 2 * dof + j] * g[j];
                }
                f[ri] -= s;
            }
        }
        let pinned = dof;
        let nred = n - pinned;
        let mut map = vec![usize::MAX; n];
        for (new, old) in (pinned..n).enumerate() {
            map[old] = new;
        }
        let mut kr = Dense::zeros(nred);
        for i in pinned..n {
            for j in pinned..n {
                kr.set(map[i], map[j], k.get(i, j));
            }
        }
        let mut fr: Vec<f64> = (pinned..n).map(|i| f[i]).collect();
        solve(&mut kr, &mut fr).expect("singular unit cell");
        let mut q = vec![0.0; n];
        for i in pinned..n {
            q[i] = fr[map[i]];
        }
        q
    }

    /// Axial force (tension positive) in every member under a macroscopic
    /// strain, from the linear solution.
    pub fn axial_forces(&self, strain: [f64; 3]) -> Vec<f64> {
        let q = self.fluctuation(strain);
        let dof = self.dof();
        let mut out = Vec::with_capacity(self.members.len());
        for m in &self.members {
            let xa = self.physical(m.a, m.a_img);
            let xb = self.physical(m.b, m.b_img);
            let dx = xb[0] - xa[0];
            let dy = xb[1] - xa[1];
            let l = (dx * dx + dy * dy).sqrt();
            let (c, s) = (dx / l, dy / l);
            let mut idx = Vec::with_capacity(2 * dof);
            for t in 0..dof {
                idx.push(m.a * dof + t);
            }
            for t in 0..dof {
                idx.push(m.b * dof + t);
            }
            let ga = self.endpoint_g(xa, strain);
            let gb = self.endpoint_g(xb, strain);
            let g: Vec<f64> = ga.iter().chain(gb.iter()).copied().collect();
            let d: Vec<f64> = (0..2 * dof).map(|i| q[idx[i]] + g[i]).collect();
            match self.kind {
                Kind::Truss => {
                    let du = [d[2] - d[0], d[3] - d[1]];
                    out.push(m.e * m.area / l * (du[0] * c + du[1] * s));
                }
                Kind::Frame => {
                    let phi = match self.nu {
                        Some(nu) => timoshenko_phi(l, m.e, m.area, m.inertia, nu),
                        None => 0.0,
                    };
                    let kl = frame2d_local(l, m.e, m.area, m.inertia, phi);
                    let r = [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]];
                    let mut dl = [0.0; 6];
                    for node in 0..2 {
                        for mm in 0..3 {
                            for nn in 0..3 {
                                dl[node * 3 + mm] += r[mm][nn] * d[node * 3 + nn];
                            }
                        }
                    }
                    let f: f64 = (0..6).map(|j| kl[0][j] * dl[j]).sum();
                    // f is the local force at node a in the +x direction, so
                    // tension is -f.
                    out.push(-f);
                }
            }
        }
        out
    }

    /// Full 3x3 effective stiffness (engineering Voigt) and the cell area.
    pub fn stiffness(&self) -> [[f64; 3]; 3] {
        let basis = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let u: Vec<f64> = basis.iter().map(|&e| self.energy(e)).collect();
        let a = self.area();
        let mut c = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let mut s = [0.0; 3];
                for k in 0..3 {
                    s[k] = basis[i][k] + basis[j][k];
                }
                c[i][j] = (self.energy(s) - u[i] - u[j]) / a;
            }
        }
        c
    }
}

/// In-plane engineering moduli from a 3x3 stiffness.
#[derive(Clone, Copy, Debug)]
pub struct Moduli2d {
    pub ex: f64,
    pub ey: f64,
    pub nu_xy: f64,
    pub nu_yx: f64,
    pub gxy: f64,
}

pub fn moduli2d(c: &[[f64; 3]; 3]) -> Moduli2d {
    let m = Dense {
        n: 3,
        a: c.iter().flat_map(|r| r.iter().copied()).collect(),
    };
    let s = crate::linalg::invert(&m).expect("singular stiffness");
    Moduli2d {
        ex: 1.0 / s.get(0, 0),
        ey: 1.0 / s.get(1, 1),
        nu_xy: -s.get(1, 0) / s.get(0, 0),
        nu_yx: -s.get(0, 1) / s.get(1, 1),
        gxy: 1.0 / s.get(2, 2),
    }
}

// ---------------------------------------------------------------------------
// 3D
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
pub struct RawMember3d {
    pub a: [f64; 3],
    pub b: [f64; 3],
    pub e: f64,
    pub area: f64,
    pub iy: f64,
    pub iz: f64,
    pub j: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Member3d {
    pub a: usize,
    pub b: usize,
    pub a_img: [i32; 3],
    pub b_img: [i32; 3],
    pub e: f64,
    pub area: f64,
    pub iy: f64,
    pub iz: f64,
    pub j: f64,
}

/// 3D periodic beam lattice.
#[derive(Clone, Debug)]
pub struct Cell3d {
    pub nodes: Vec<[f64; 3]>,
    pub members: Vec<Member3d>,
    pub a1: [f64; 3],
    pub a2: [f64; 3],
    pub a3: [f64; 3],
    pub kind: Kind,
    /// Poisson ratio for the 3D frame torsion/bending shear modulus.
    pub nu: f64,
}

fn dedup3d(members: Vec<Member3d>) -> Vec<Member3d> {
    let mut seen = std::collections::HashSet::new();
    let mut out = Vec::new();
    for m in members {
        let key = if m.a <= m.b {
            (
                m.a,
                m.b,
                [
                    m.b_img[0] - m.a_img[0],
                    m.b_img[1] - m.a_img[1],
                    m.b_img[2] - m.a_img[2],
                ],
            )
        } else {
            (
                m.b,
                m.a,
                [
                    m.a_img[0] - m.b_img[0],
                    m.a_img[1] - m.b_img[1],
                    m.a_img[2] - m.b_img[2],
                ],
            )
        };
        if seen.insert(key) {
            out.push(m);
        }
    }
    out
}

fn reduce3d(
    nodes: &mut Vec<[f64; 3]>,
    a1: [f64; 3],
    a2: [f64; 3],
    a3: [f64; 3],
    p: [f64; 3],
) -> (usize, [i32; 3]) {
    // Solve the 3x3 system [a1 a2 a3] n = p by Cramer's rule.
    let det = a1[0] * (a2[1] * a3[2] - a2[2] * a3[1])
        - a1[1] * (a2[0] * a3[2] - a2[2] * a3[0])
        + a1[2] * (a2[0] * a3[1] - a2[1] * a3[0]);
    assert!(det.abs() > 0.0);
    let d1 = p[0] * (a2[1] * a3[2] - a2[2] * a3[1])
        - p[1] * (a2[0] * a3[2] - a2[2] * a3[0])
        + p[2] * (a2[0] * a3[1] - a2[1] * a3[0]);
    let d2 = a1[0] * (p[1] * a3[2] - p[2] * a3[1])
        - a1[1] * (p[0] * a3[2] - p[2] * a3[0])
        + a1[2] * (p[0] * a3[1] - p[1] * a3[0]);
    let d3 = a1[0] * (a2[1] * p[2] - a2[2] * p[1])
        - a1[1] * (a2[0] * p[2] - a2[2] * p[0])
        + a1[2] * (a2[0] * p[1] - a2[1] * p[0]);
    let n1 = exact_floor(d1 / det);
    let n2 = exact_floor(d2 / det);
    let n3 = exact_floor(d3 / det);
    let x = [
        p[0] - n1 * a1[0] - n2 * a2[0] - n3 * a3[0],
        p[1] - n1 * a1[1] - n2 * a2[1] - n3 * a3[1],
        p[2] - n1 * a1[2] - n2 * a2[2] - n3 * a3[2],
    ];
    for (i, q) in nodes.iter().enumerate() {
        if (q[0] - x[0]).abs() < 1e-8 && (q[1] - x[1]).abs() < 1e-8 && (q[2] - x[2]).abs() < 1e-8
        {
            return (i, [n1 as i32, n2 as i32, n3 as i32]);
        }
    }
    nodes.push(x);
    (
        nodes.len() - 1,
        [n1 as i32, n2 as i32, n3 as i32],
    )
}

impl Cell3d {
    pub fn new(
        a1: [f64; 3],
        a2: [f64; 3],
        a3: [f64; 3],
        raw: &[RawMember3d],
        kind: Kind,
    ) -> Self {
        let mut nodes = Vec::new();
        let mut members = Vec::new();
        for r in raw {
            let (a, ai) = reduce3d(&mut nodes, a1, a2, a3, r.a);
            let (b, bi) = reduce3d(&mut nodes, a1, a2, a3, r.b);
            members.push(Member3d {
                a,
                b,
                a_img: ai,
                b_img: bi,
                e: r.e,
                area: r.area,
                iy: r.iy,
                iz: r.iz,
                j: r.j,
            });
        }
        Self {
            nodes,
            members: dedup3d(members),
            a1,
            a2,
            a3,
            kind,
            nu: 0.3,
        }
    }

    pub fn volume(&self) -> f64 {
        let a1 = self.a1;
        let a2 = self.a2;
        let a3 = self.a3;
        (a1[0] * (a2[1] * a3[2] - a2[2] * a3[1]) - a1[1] * (a2[0] * a3[2] - a2[2] * a3[0])
            + a1[2] * (a2[0] * a3[1] - a2[1] * a3[0]))
            .abs()
    }

    /// Member volume fraction of the unit cell (relative density).
    pub fn relative_density(&self) -> f64 {
        let mut v = 0.0;
        for m in &self.members {
            let xa = self.physical(m.a, m.a_img);
            let xb = self.physical(m.b, m.b_img);
            v += ((xb[0] - xa[0]).powi(2) + (xb[1] - xa[1]).powi(2) + (xb[2] - xa[2]).powi(2))
                .sqrt()
                * m.area;
        }
        v / self.volume()
    }

    pub fn physical(&self, node: usize, img: [i32; 3]) -> [f64; 3] {
        let x = self.nodes[node];
        let mut p = x;
        for k in 0..3 {
            p[k] += img[0] as f64 * self.a1[k]
                + img[1] as f64 * self.a2[k]
                + img[2] as f64 * self.a3[k];
        }
        p
    }

    pub fn affine(&self, x: [f64; 3], strain: [f64; 6]) -> [f64; 3] {
        let [exx, eyy, ezz, gxy, gyz, gzx] = strain;
        [
            exx * x[0] + 0.5 * gxy * x[1] + 0.5 * gzx * x[2],
            0.5 * gxy * x[0] + eyy * x[1] + 0.5 * gyz * x[2],
            0.5 * gzx * x[0] + 0.5 * gyz * x[1] + ezz * x[2],
        ]
    }

    fn dof(&self) -> usize {
        match self.kind {
            Kind::Frame => 6,
            Kind::Truss => 3,
        }
    }

    fn element_k(&self, xa: [f64; 3], xb: [f64; 3], m: &Member3d) -> Vec<f64> {
        let d = [xb[0] - xa[0], xb[1] - xa[1], xb[2] - xa[2]];
        let l = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let dir = [d[0] / l, d[1] / l, d[2] / l];
        match self.kind {
            Kind::Truss => truss3d_global(l, m.e, m.area, dir).to_vec(),
            Kind::Frame => {
                let r = local_frame(dir);
                let g = m.e / (2.0 * (1.0 + self.nu));
                let kl = frame3d_local(l, m.e, g, m.area, m.iy, m.iz, m.j);
                let kg = frame3d_global(&kl, &r);
                kg.iter().flat_map(|row| row.iter().copied()).collect()
            }
        }
    }

    fn endpoint_g(&self, x: [f64; 3], strain: [f64; 6]) -> Vec<f64> {
        let u = self.affine(x, strain);
        match self.kind {
            Kind::Truss => vec![u[0], u[1], u[2]],
            Kind::Frame => vec![u[0], u[1], u[2], 0.0, 0.0, 0.0],
        }
    }

    pub fn energy(&self, strain: [f64; 6]) -> f64 {
        let dof = self.dof();
        let n = self.nodes.len() * dof;
        let mut k = Dense::zeros(n);
        let mut f = vec![0.0; n];
        for m in &self.members {
            let xa = self.physical(m.a, m.a_img);
            let xb = self.physical(m.b, m.b_img);
            let ke = self.element_k(xa, xb, m);
            let mut idx = Vec::with_capacity(2 * dof);
            for t in 0..dof {
                idx.push(m.a * dof + t);
            }
            for t in 0..dof {
                idx.push(m.b * dof + t);
            }
            k.add_outer_block(&idx, &ke);
            let ga = self.endpoint_g(xa, strain);
            let gb = self.endpoint_g(xb, strain);
            let g: Vec<f64> = ga.iter().chain(gb.iter()).copied().collect();
            for (i, &ri) in idx.iter().enumerate() {
                let mut s = 0.0;
                for j in 0..2 * dof {
                    s += ke[i * 2 * dof + j] * g[j];
                }
                f[ri] -= s;
            }
        }
        let pinned = dof;
        let nred = n - pinned;
        let mut map = vec![usize::MAX; n];
        for (new, old) in (pinned..n).enumerate() {
            map[old] = new;
        }
        let mut kr = Dense::zeros(nred);
        for i in pinned..n {
            for j in pinned..n {
                kr.set(map[i], map[j], k.get(i, j));
            }
        }
        let mut fr: Vec<f64> = (pinned..n).map(|i| f[i]).collect();
        solve(&mut kr, &mut fr).expect("singular unit cell");
        let mut q = vec![0.0; n];
        for i in pinned..n {
            q[i] = fr[map[i]];
        }
        let mut u = 0.0;
        for m in &self.members {
            let xa = self.physical(m.a, m.a_img);
            let xb = self.physical(m.b, m.b_img);
            let ke = self.element_k(xa, xb, m);
            let mut idx = Vec::with_capacity(2 * dof);
            for t in 0..dof {
                idx.push(m.a * dof + t);
            }
            for t in 0..dof {
                idx.push(m.b * dof + t);
            }
            let ga = self.endpoint_g(xa, strain);
            let gb = self.endpoint_g(xb, strain);
            let g: Vec<f64> = ga.iter().chain(gb.iter()).copied().collect();
            let d: Vec<f64> = (0..2 * dof).map(|i| q[idx[i]] + g[i]).collect();
            let mut e = 0.0;
            for i in 0..2 * dof {
                for j in 0..2 * dof {
                    e += d[i] * ke[i * 2 * dof + j] * d[j];
                }
            }
            u += 0.5 * e;
        }
        u
    }

    /// Fluctuation (periodic part) nodal displacement for a macroscopic
    /// strain; the pinned node's entries are zero.
    pub fn fluctuation(&self, strain: [f64; 6]) -> Vec<f64> {
        let dof = self.dof();
        let n = self.nodes.len() * dof;
        let mut k = Dense::zeros(n);
        let mut f = vec![0.0; n];
        for m in &self.members {
            let xa = self.physical(m.a, m.a_img);
            let xb = self.physical(m.b, m.b_img);
            let ke = self.element_k(xa, xb, m);
            let mut idx = Vec::with_capacity(2 * dof);
            for t in 0..dof {
                idx.push(m.a * dof + t);
            }
            for t in 0..dof {
                idx.push(m.b * dof + t);
            }
            k.add_outer_block(&idx, &ke);
            let ga = self.endpoint_g(xa, strain);
            let gb = self.endpoint_g(xb, strain);
            let g: Vec<f64> = ga.iter().chain(gb.iter()).copied().collect();
            for (i, &ri) in idx.iter().enumerate() {
                let mut s = 0.0;
                for j in 0..2 * dof {
                    s += ke[i * 2 * dof + j] * g[j];
                }
                f[ri] -= s;
            }
        }
        let pinned = dof;
        let nred = n - pinned;
        let mut map = vec![usize::MAX; n];
        for (new, old) in (pinned..n).enumerate() {
            map[old] = new;
        }
        let mut kr = Dense::zeros(nred);
        for i in pinned..n {
            for j in pinned..n {
                kr.set(map[i], map[j], k.get(i, j));
            }
        }
        let mut fr: Vec<f64> = (pinned..n).map(|i| f[i]).collect();
        solve(&mut kr, &mut fr).expect("singular unit cell");
        let mut q = vec![0.0; n];
        for i in pinned..n {
            q[i] = fr[map[i]];
        }
        q
    }

    /// Axial force (tension positive) in every member under a macroscopic
    /// strain.
    pub fn axial_forces(&self, strain: [f64; 6]) -> Vec<f64> {
        let q = self.fluctuation(strain);
        let dof = self.dof();
        let mut out = Vec::with_capacity(self.members.len());
        for m in &self.members {
            let xa = self.physical(m.a, m.a_img);
            let xb = self.physical(m.b, m.b_img);
            let dvec = [xb[0] - xa[0], xb[1] - xa[1], xb[2] - xa[2]];
            let l = (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]).sqrt();
            let dir = [dvec[0] / l, dvec[1] / l, dvec[2] / l];
            let mut idx = Vec::with_capacity(2 * dof);
            for t in 0..dof {
                idx.push(m.a * dof + t);
            }
            for t in 0..dof {
                idx.push(m.b * dof + t);
            }
            let ga = self.endpoint_g(xa, strain);
            let gb = self.endpoint_g(xb, strain);
            let g: Vec<f64> = ga.iter().chain(gb.iter()).copied().collect();
            let d: Vec<f64> = (0..2 * dof).map(|i| q[idx[i]] + g[i]).collect();
            match self.kind {
                Kind::Truss => {
                    let du = [d[3] - d[0], d[4] - d[1], d[5] - d[2]];
                    out.push(m.e * m.area / l * (du[0] * dir[0] + du[1] * dir[1] + du[2] * dir[2]));
                }
                Kind::Frame => {
                    let r = local_frame(dir);
                    let g = m.e / (2.0 * (1.0 + self.nu));
                    let kl = frame3d_local(l, m.e, g, m.area, m.iy, m.iz, m.j);
                    let mut dl = [0.0; 12];
                    for node in 0..4 {
                        for mm in 0..3 {
                            for nn in 0..3 {
                                dl[node * 3 + mm] += r[mm][nn] * d[node * 3 + nn];
                            }
                        }
                    }
                    let f: f64 = (0..12).map(|j| kl[0][j] * dl[j]).sum();
                    // tension is -f for the same reason as 2D.
                    out.push(-f);
                }
            }
        }
        out
    }

    /// Full 6x6 effective stiffness and the cell volume.
    pub fn stiffness(&self) -> [[f64; 6]; 6] {
        let mut basis = [[0.0f64; 6]; 6];
        for i in 0..6 {
            basis[i][i] = 1.0;
        }
        let u: Vec<f64> = basis.iter().map(|&e| self.energy(e)).collect();
        let v = self.volume();
        let mut c = [[0.0f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                let mut s = [0.0; 6];
                for k in 0..6 {
                    s[k] = basis[i][k] + basis[j][k];
                }
                c[i][j] = (self.energy(s) - u[i] - u[j]) / v;
            }
        }
        c
    }
}

/// Isotropic engineering moduli from a 6x6 stiffness.
#[derive(Clone, Copy, Debug)]
pub struct Moduli3d {
    pub e: f64,
    pub nu: f64,
    pub k: f64,
    pub g: f64,
}

pub fn moduli3d(c: &[[f64; 6]; 6]) -> Moduli3d {
    let m = Dense {
        n: 6,
        a: c.iter().flat_map(|r| r.iter().copied()).collect(),
    };
    let s = crate::linalg::invert(&m).expect("singular stiffness");
    let e = 1.0 / s.get(0, 0);
    let nu = -s.get(1, 0) / s.get(0, 0);
    // Bulk modulus from the stiffness (the 3x3 volumetric block).
    let mut mean = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            mean[i][j] = c[i][j];
        }
    }
    let km = Dense {
        n: 3,
        a: mean.iter().flat_map(|r| r.iter().copied()).collect(),
    };
    let kmi = crate::linalg::invert(&km).unwrap();
    // K = 1 / sum(S_vol) for a hydrostatic stress.
    let mut vol = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            vol += kmi.get(i, j);
        }
    }
    Moduli3d {
        e,
        nu,
        k: 1.0 / vol,
        g: 1.0 / s.get(3, 3),
    }
}

/// Relative density of a 2D lattice: member area fraction.
pub fn relative_density_2d(members: &[RawMember2d], area: f64) -> f64 {
    let mut a = 0.0;
    for m in members {
        let dx = m.b[0] - m.a[0];
        let dy = m.b[1] - m.a[1];
        let l = (dx * dx + dy * dy).sqrt();
        a += l * m.area;
    }
    a / area
}

/// Relative density of a 3D lattice: member volume fraction.
pub fn relative_density_3d(members: &[RawMember3d], volume: f64) -> f64 {
    let mut a = 0.0;
    for m in members {
        let d = [m.b[0] - m.a[0], m.b[1] - m.a[1], m.b[2] - m.a[2]];
        let l = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        a += l * m.area;
    }
    a / volume
}

/// `2 pi` for the dispersion module.
pub const TWO_PI: f64 = 2.0 * PI;
