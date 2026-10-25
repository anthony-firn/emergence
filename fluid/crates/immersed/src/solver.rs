//! Incompressible Navier–Stokes on a MAC grid with immersed geometry.
//!
//! Discretisation: staggered velocities, cell-centred (kinematic) pressure;
//! conservative second-order central convection (Harlow & Welch 1965) and
//! diffusion. Time: SSP-RK3 (Shu & Osher 1988) with a Leray projection at
//! every stage and an incremental pressure, so a steady state satisfies the
//! discrete steady equations and the immersed-boundary condition exactly,
//! independent of the time step.
//!
//! Immersed geometry: sharp-interface direct forcing (Fadlun et al., J.
//! Comput. Phys. 161, 2000; Mittal & Iaccarino, Annu. Rev. Fluid Mech. 37,
//! 2005). Nodes inside the body take the body velocity; each fluid node next
//! to it is set by linear interpolation along the surface normal between
//! the wall and a probe point `√dim · h` further out, so the no-slip
//! condition is imposed on the true (CAD) surface to second order. Because
//! the body is immersed, the pressure equation stays the box Laplacian and
//! is solved exactly by [`crate::poisson`].

use std::sync::Arc;

use geometry::{Sdf, P3};
use rayon::prelude::*;

use crate::field::{Arr, Grid};
use crate::poisson::{End, Poisson};

/// Prescribed velocity as a function of position and time.
pub type VelocityFn = Arc<dyn Fn(P3, f64) -> P3 + Send + Sync>;

/// Condition on one face of the box.
#[derive(Clone)]
pub enum Bc {
    /// Prescribed velocity: no-slip wall (zero), moving wall, or inlet.
    Velocity(VelocityFn),
    /// Pressure outlet: `p = 0`, zero normal gradient of velocity.
    Outlet,
    /// Symmetry / free slip: zero normal velocity, zero tangential stress.
    Slip,
    /// Periodic (must be paired with the opposite face).
    Periodic,
}

impl Bc {
    pub fn wall() -> Self {
        Bc::Velocity(Arc::new(|_, _| [0.0; 3]))
    }
    pub fn moving_wall(v: P3) -> Self {
        Bc::Velocity(Arc::new(move |_, _| v))
    }
    pub fn inlet(f: impl Fn(P3, f64) -> P3 + Send + Sync + 'static) -> Self {
        Bc::Velocity(Arc::new(f))
    }
    fn end(&self) -> End {
        match self {
            Bc::Outlet => End::Dirichlet,
            Bc::Periodic => End::Periodic,
            _ => End::Neumann,
        }
    }
}

/// A rigid immersed body: surface plus rigid motion
/// `u_b(x) = velocity + angular × (x − center)` (the geometry is fixed; the
/// motion enters through the no-slip condition, e.g. a rotating shaft).
pub struct Body {
    pub sdf: Sdf,
    pub velocity: P3,
    pub angular: P3,
    pub center: P3,
}

impl Body {
    /// A body at rest.
    pub fn fixed(sdf: Sdf) -> Self {
        Self {
            sdf,
            velocity: [0.0; 3],
            angular: [0.0; 3],
            center: [0.0; 3],
        }
    }

    /// Rigid-body velocity at `x`.
    pub fn velocity_at(&self, x: P3) -> P3 {
        let r = [x[0] - self.center[0], x[1] - self.center[1], x[2] - self.center[2]];
        let w = self.angular;
        [
            self.velocity[0] + w[1] * r[2] - w[2] * r[1],
            self.velocity[1] + w[2] * r[0] - w[0] * r[2],
            self.velocity[2] + w[0] * r[1] - w[1] * r[0],
        ]
    }
}

/// One forcing node at distance `d` from the wall: its value is the
/// Lagrange interpolant through the wall velocity `u_w` (at 0) and probes
/// at `d + Δ` and `d + 2Δ` along the normal, evaluated at `d`:
/// `u = l₀ u_w + l₁ u(p₁) + l₂ u(p₂)` (linear forcing: `l₂ = 0`).
#[derive(Clone, Copy, Debug)]
struct Forcing {
    node: usize,
    l: [f64; 3],
    stencil: [(usize, f64); 8],
    stencil2: [(usize, f64); 8],
    wall: f64,
    body: usize,
}

#[derive(Clone, Debug, Default)]
struct CompIb {
    /// Solid node, its body velocity component, and its body.
    solid: Vec<(usize, f64, usize)>,
    forcing: Vec<Forcing>,
}

/// One precomputed boundary assignment (flat indices into a component).
#[derive(Clone, Copy, Debug)]
enum FillOp {
    Copy(usize, usize),
    Zero(usize),
    /// Outlet normal velocity: zero gradient before a projection only.
    Extrapolate(usize, usize),
    /// Prescribed velocity of face `face` (= 2·axis + side) at `x`.
    Set { dst: usize, face: u8, x: P3 },
    /// Tangential ghost `2 u_face(x) − u_src`.
    Mirror { dst: usize, src: usize, face: u8, x: P3 },
}

/// Per-axis range of nodes updated by the momentum equation.
type Range3 = [(isize, isize); 3];

/// The solver state.
pub struct Solver {
    pub grid: Grid,
    /// Kinematic viscosity and density.
    pub nu: f64,
    pub rho: f64,
    pub bcs: [[Bc; 2]; 3],
    pub vel: [Arr; 3],
    /// Kinematic pressure `p/ρ` at cell centres (`x` fastest, no ghosts).
    pub p: Vec<f64>,
    pub t: f64,
    pub bodies: Vec<Body>,
    poisson: Poisson,
    ib: [CompIb; 3],
    /// Cells with a face held by the immersed boundary (mass-source cells).
    ib_cell: Vec<bool>,
    /// Use the mass-source treatment of [`Self::ib_cell`] (diagnostics).
    pub mass_source: bool,
    /// Quadratic (default) or linear wall reconstruction at forcing nodes
    /// (change with [`Solver::set_quadratic_ib`]).
    quadratic_ib: bool,
    fill_ops: [Vec<FillOp>; 3],
    /// Axes that carry stencil contributions (not single-cell periodic).
    active: [bool; 3],
    update: [Range3; 3],
    /// Probe offset `√dim · h` used by the immersed boundary and the
    /// surface force integration.
    pub probe: f64,
}

impl Solver {
    pub fn new(grid: Grid, nu: f64, rho: f64, bcs: [[Bc; 2]; 3], bodies: Vec<Body>) -> Self {
        for a in 0..3 {
            let per = [matches!(bcs[a][0], Bc::Periodic), matches!(bcs[a][1], Bc::Periodic)];
            assert!(per[0] == per[1], "axis {a}: periodic faces must be paired");
        }
        let ends = std::array::from_fn(|a| [bcs[a][0].end(), bcs[a][1].end()]);
        let poisson = Poisson::new(grid.n, grid.h, ends);
        let vel = std::array::from_fn(|c| Arr::new(grid.comp_dims(c)));
        let update = std::array::from_fn(|c| {
            std::array::from_fn(|a| {
                let n = grid.n[a] as isize;
                if a != c {
                    (0, n - 1)
                } else if matches!(bcs[a][0], Bc::Periodic) {
                    (0, n - 1)
                } else {
                    (1, n - 1)
                }
            })
        });
        let probe = (grid.dim().max(1) as f64).sqrt() * grid.h;
        let active = std::array::from_fn(|a| !(grid.n[a] == 1 && matches!(bcs[a][0], Bc::Periodic)));
        let mut s = Self {
            grid,
            nu,
            rho,
            bcs,
            vel,
            p: vec![0.0; grid.cells()],
            t: 0.0,
            bodies,
            poisson,
            ib: Default::default(),
            ib_cell: Vec::new(),
            mass_source: true,
            quadratic_ib: true,
            fill_ops: Default::default(),
            active,
            update,
            probe,
        };
        s.fill_ops = s.build_fill_ops();
        s.classify();
        let mut vel = std::mem::replace(&mut s.vel, std::array::from_fn(|_| Arr::new([1, 1, 1])));
        s.fill_bcs_on(&mut vel, 0.0);
        s.vel = vel;
        s
    }

    /// Choose the wall reconstruction at forcing nodes: quadratic (default;
    /// Majumdar, Iaccarino & Durbin, CTR Briefs 2001) or linear (Fadlun et
    /// al. 2000). Re-classifies the immersed boundary.
    pub fn set_quadratic_ib(&mut self, quadratic: bool) {
        self.quadratic_ib = quadratic;
        self.classify();
    }

    /// Set the initial velocity from a function, then project it.
    pub fn initialize(&mut self, f: impl Fn(P3) -> P3 + Sync) {
        for c in 0..3 {
            let g = self.grid;
            let arr = &mut self.vel[c];
            let d = arr.d;
            for k in 0..d[2] as isize {
                for j in 0..d[1] as isize {
                    for i in 0..d[0] as isize {
                        let v = f(g.node(c, i, j, k))[c];
                        arr.set(i, j, k, v);
                    }
                }
            }
        }
        let mut vel = std::mem::replace(&mut self.vel, std::array::from_fn(|_| Arr::new([1, 1, 1])));
        self.fill_bcs_on(&mut vel, self.t);
        self.apply_ib(&mut vel);
        self.project(&mut vel, 1.0, false);
        self.fill_ghosts(&mut vel, self.t);
        self.vel = vel;
    }

    // ---------------------------------------------------------------- IB

    fn classify(&mut self) {
        if self.bodies.is_empty() {
            return;
        }
        let g = self.grid;
        let probe = self.probe;
        let quadratic = self.quadratic_ib;
        let bodies = &self.bodies;
        for c in 0..3 {
            let arr = Arr::new(g.comp_dims(c));
            let d = arr.d;
            let nodes: Vec<usize> = (0..arr.data.len())
                .filter(|&q| {
                    let [i, j, k] = arr.ijk(q);
                    i >= 0 && j >= 0 && k >= 0 && i < d[0] as isize && j < d[1] as isize && k < d[2] as isize
                })
                .collect();
            // Signed distance to the union of bodies (only 8 bytes per node).
            // Far from a body (outside its box + margin) no query is needed:
            // the point is fluid, or solid for an inside-out body.
            let margin = 3.0 * probe;
            let boxes: Vec<(P3, P3, bool)> = bodies
                .iter()
                .map(|b| {
                    let (lo, hi) = b.sdf.mesh.bbox();
                    (lo.map(|v| v - margin), hi.map(|v| v + margin), b.sdf.mesh.volume() < 0.0)
                })
                .collect();
            let far = 1e3 * (g.h * g.n.iter().sum::<usize>() as f64);
            let phi_of = |x: P3| -> f64 {
                bodies
                    .iter()
                    .zip(&boxes)
                    .map(|(body, (lo, hi, inside_out))| {
                        if (0..3).all(|a| x[a] >= lo[a] && x[a] <= hi[a]) {
                            body.sdf.distance(x)
                        } else if *inside_out {
                            -far
                        } else {
                            far
                        }
                    })
                    .fold(f64::INFINITY, f64::min)
            };
            let dist: Vec<f64> = nodes
                .par_iter()
                .map(|&q| {
                    let [i, j, k] = arr.ijk(q);
                    phi_of(g.node(c, i, j, k))
                })
                .collect();
            let mut phi = vec![f64::INFINITY; arr.data.len()];
            for (&q, &d) in nodes.iter().zip(&dist) {
                phi[q] = d;
            }
            let st = arr.strides();
            let mut ib = CompIb::default();
            // Solid nodes, and fluid nodes next to the solid (forcing).
            let mut candidates = Vec::new();
            for (&q, &d0) in nodes.iter().zip(&dist) {
                if d0 < 0.0 {
                    candidates.push((q, true));
                } else if (0..3).any(|a| d[a] > 1 && (phi[q + st[a]] < 0.0 || phi[q - st[a]] < 0.0)) {
                    candidates.push((q, false));
                }
            }
            // Exact closest points only for these (wall point, normal, body).
            let info: Vec<(usize, bool, f64, P3, P3, usize)> = candidates
                .par_iter()
                .map(|&(q, solid)| {
                    let [i, j, k] = arr.ijk(q);
                    let x = g.node(c, i, j, k);
                    let (dd, n, w, b) = bodies
                        .iter()
                        .enumerate()
                        .map(|(b, body)| {
                            let cl = body.sdf.closest(x);
                            (cl.distance, cl.normal, cl.point, b)
                        })
                        .min_by(|a, b| a.0.total_cmp(&b.0))
                        .unwrap();
                    (q, solid, dd, n, w, b)
                })
                .collect();
            for (q, solid, dist, normal, wall, b) in info {
                let [i, j, k] = arr.ijk(q);
                let x = g.node(c, i, j, k);
                if solid {
                    ib.solid.push((q, bodies[b].velocity_at(x)[c], b));
                    continue;
                }
                let frac = |s: f64| -> P3 { std::array::from_fn(|a| (x[a] + s * normal[a] - g.origin[a]) / g.h - if a == c { 0.0 } else { 0.5 }) };
                let (d1, d2) = (dist + probe, dist + 2.0 * probe);
                let l = if quadratic {
                    [
                        (dist - d1) * (dist - d2) / (d1 * d2),
                        dist * (dist - d2) / (d1 * (d1 - d2)),
                        dist * (dist - d1) / (d2 * (d2 - d1)),
                    ]
                } else {
                    [1.0 - dist / d1, dist / d1, 0.0]
                };
                ib.forcing.push(Forcing {
                    node: q,
                    l,
                    stencil: arr.stencil(frac(probe)),
                    stencil2: arr.stencil(frac(2.0 * probe)),
                    wall: bodies[b].velocity_at(wall)[c],
                    body: b,
                });
            }
            self.ib[c] = ib;
        }
        // Cells whose faces are *all* held by the body (along the active
        // axes) have their divergence dictated by the immersed boundary: the
        // projection takes it as a mass source (Kim, Kim & Choi, J. Comput.
        // Phys. 171, 2001) instead of fighting the forcing, which otherwise
        // makes the incremental pressure drift. Cells with a free face keep
        // their continuity constraint (masking those lets the free face run
        // away: observed at a stagnation point).
        let held: [Vec<bool>; 3] = std::array::from_fn(|c| {
            let arr = Arr::new(g.comp_dims(c));
            let mut h = vec![false; arr.data.len()];
            for q in self.ib[c].solid.iter().map(|e| e.0).chain(self.ib[c].forcing.iter().map(|f| f.node)) {
                h[q] = true;
            }
            h
        });
        let arrs: [Arr; 3] = std::array::from_fn(|c| Arr::new(g.comp_dims(c)));
        let mut mask = vec![false; g.cells()];
        for k in 0..g.n[2] {
            for j in 0..g.n[1] {
                for i in 0..g.n[0] {
                    let (ii, jj, kk) = (i as isize, j as isize, k as isize);
                    let all = (0..3).filter(|&c| g.n[c] > 1).all(|c| {
                        let mut hi = [ii, jj, kk];
                        hi[c] += 1;
                        held[c][arrs[c].at(ii, jj, kk)] && held[c][arrs[c].at(hi[0], hi[1], hi[2])]
                    });
                    mask[i + g.n[0] * (j + g.n[1] * k)] = all;
                }
            }
        }
        self.ib_cell = mask;
    }

    /// Counts of solid and forcing nodes per component.
    pub fn ib_counts(&self) -> [(usize, usize); 3] {
        std::array::from_fn(|c| (self.ib[c].solid.len(), self.ib[c].forcing.len()))
    }

    fn apply_ib(&self, vel: &mut [Arr; 3]) {
        for c in 0..3 {
            let ib = &self.ib[c];
            let arr = &mut vel[c];
            for &(q, ub, _) in &ib.solid {
                arr.data[q] = ub;
            }
            // Forcing nodes may sit in each other's probe stencils:
            // Gauss–Seidel sweeps (weights d/(d+Δ) < 1/2 contract fast).
            for _ in 0..4 {
                for f in &ib.forcing {
                    let u1 = arr.interp(&f.stencil);
                    let u2 = if f.l[2] != 0.0 { arr.interp(&f.stencil2) } else { 0.0 };
                    arr.data[f.node] = f.l[0] * f.wall + f.l[1] * u1 + f.l[2] * u2;
                }
            }
        }
    }

    /// Set boundary-normal nodes and ghosts of every component at time `t`.
    /// Outlet faces extrapolate the normal velocity (zero gradient) before a
    /// projection; afterwards ([`Self::fill_ghosts`]) the projected outlet
    /// value is kept.
    pub fn fill_bcs_on(&self, vel: &mut [Arr; 3], t: f64) {
        self.fill(vel, t, true)
    }

    fn fill_ghosts(&self, vel: &mut [Arr; 3], t: f64) {
        self.fill(vel, t, false)
    }

    /// Build the boundary-fill program: per component, the ordered list of
    /// ghost/boundary assignments (axis by axis, so edges and corners take
    /// already-filled ghosts). Axes that are periodic with a single cell
    /// carry no stencil contributions and are skipped.
    fn build_fill_ops(&self) -> [Vec<FillOp>; 3] {
        let g = self.grid;
        std::array::from_fn(|c| {
            let arr = Arr::new(g.comp_dims(c));
            let d = arr.d.map(|v| v as isize);
            let mut ops = Vec::new();
            for a in 0..3 {
                let flat = g.n[a] == 1 && matches!(self.bcs[a][0], Bc::Periodic);
                // A single periodic cell: no stencil reaches across it, and
                // the component along it (identically zero in 2D) needs none.
                if flat || !self.active[c] {
                    continue;
                }
                let (b1, b2) = ((a + 1) % 3, (a + 2) % 3);
                for side in 0..2 {
                    let bc = &self.bcs[a][side];
                    let face = (2 * a + side) as u8;
                    for m2 in -1..=d[b2] {
                        for m1 in -1..=d[b1] {
                            let mut idx = [0isize; 3];
                            idx[b1] = m1;
                            idx[b2] = m2;
                            let at = |v: isize| {
                                let mut i = idx;
                                i[a] = v;
                                arr.at(i[0], i[1], i[2])
                            };
                            let n = d[a];
                            if a == c {
                                let (bnd, inner, ghost) = if side == 0 { (0, 1, -1) } else { (n - 1, n - 2, n) };
                                let mut ib = idx;
                                ib[a] = bnd;
                                match bc {
                                    Bc::Periodic => {
                                        if side == 1 {
                                            ops.push(FillOp::Copy(at(n - 1), at(0)));
                                            ops.push(FillOp::Copy(at(n), at(1)));
                                            ops.push(FillOp::Copy(at(-1), at(n - 2)));
                                        }
                                    }
                                    Bc::Velocity(_) => {
                                        let x = g.node(c, ib[0], ib[1], ib[2]);
                                        ops.push(FillOp::Set { dst: at(bnd), face, x });
                                        ops.push(FillOp::Copy(at(ghost), at(bnd)));
                                    }
                                    Bc::Slip => {
                                        ops.push(FillOp::Zero(at(bnd)));
                                        ops.push(FillOp::Zero(at(ghost)));
                                    }
                                    Bc::Outlet => {
                                        ops.push(FillOp::Extrapolate(at(bnd), at(inner)));
                                        ops.push(FillOp::Copy(at(ghost), at(bnd)));
                                    }
                                }
                            } else {
                                let (adj, ghost, opp) = if side == 0 { (0, -1, n - 1) } else { (n - 1, n, 0) };
                                match bc {
                                    Bc::Periodic => ops.push(FillOp::Copy(at(ghost), at(opp))),
                                    Bc::Velocity(_) => {
                                        let mut ia = idx;
                                        ia[a] = adj;
                                        let mut x = g.node(c, ia[0], ia[1], ia[2]);
                                        x[a] = if side == 0 { g.origin[a] } else { g.origin[a] + g.h * g.n[a] as f64 };
                                        ops.push(FillOp::Mirror { dst: at(ghost), src: at(adj), face, x });
                                    }
                                    Bc::Slip | Bc::Outlet => ops.push(FillOp::Copy(at(ghost), at(adj))),
                                }
                            }
                        }
                    }
                }
            }
            ops
        })
    }

    fn fill(&self, vel: &mut [Arr; 3], t: f64, extrapolate_outlet: bool) {
        for c in 0..3 {
            let data = &mut vel[c].data;
            for op in &self.fill_ops[c] {
                match *op {
                    FillOp::Copy(dst, src) => data[dst] = data[src],
                    FillOp::Zero(dst) => data[dst] = 0.0,
                    FillOp::Extrapolate(dst, src) => {
                        if extrapolate_outlet {
                            data[dst] = data[src];
                        }
                    }
                    FillOp::Set { dst, face, x } => data[dst] = self.face_velocity(face, x, t)[c],
                    FillOp::Mirror { dst, src, face, x } => data[dst] = 2.0 * self.face_velocity(face, x, t)[c] - data[src],
                }
            }
        }
    }

    #[inline]
    fn face_velocity(&self, face: u8, x: P3, t: f64) -> P3 {
        match &self.bcs[face as usize / 2][face as usize % 2] {
            Bc::Velocity(f) => f(x, t),
            _ => [0.0; 3],
        }
    }

    // ---------------------------------------------------------- operators

    /// `R = −∇·(u u) + ν∇²u − ∇p` at every updated node (zero elsewhere).
    fn rhs(&self, vel: &[Arr; 3]) -> [Arr; 3] {
        let g = self.grid;
        let h = g.h;
        let (ih, ih2, nu) = (1.0 / h, 1.0 / (h * h), self.nu);
        let n = g.n.map(|v| v as isize);
        let p = &self.p;
        let active = self.active;
        let pidx = |i: isize, j: isize, k: isize| -> usize {
            let w = |v: isize, a: usize| v.rem_euclid(n[a]) as usize;
            w(i, 0) + g.n[0] * (w(j, 1) + g.n[1] * w(k, 2))
        };
        std::array::from_fn(|c| {
            let q = &vel[c];
            let mut out = Arr::new(q.d);
            if !active[c] {
                return out;
            }
            let st = q.strides();
            let r = self.update[c];
            let sx = st[1];
            let rows_per_plane = q.d[1] + 2;
            out.data.par_chunks_mut(sx).enumerate().for_each(|(row, line)| {
                let (j, k) = ((row % rows_per_plane) as isize - 1, (row / rows_per_plane) as isize - 1);
                if k < r[2].0 || k > r[2].1 || j < r[1].0 || j > r[1].1 {
                    return;
                }
                for i in r[0].0..=r[0].1 {
                    let idx = [i, j, k];
                    let m = q.at(i, j, k);
                    let qc = q.data[m];
                    let mut conv = 0.0;
                    let mut lap = 0.0;
                    for a in 0..3 {
                        if !active[a] {
                            continue;
                        }
                        let (qp, qm) = (q.data[m + st[a]], q.data[m - st[a]]);
                        lap += qp + qm - 2.0 * qc;
                        if a == c {
                            let fp = 0.5 * (qc + qp);
                            let fm = 0.5 * (qm + qc);
                            conv += (fp * fp - fm * fm) * ih;
                        } else {
                            let arr = &vel[a];
                            // component a at the faces bounding the flux points
                            let mut ip = idx;
                            ip[a] += 1;
                            let mut ipm = ip;
                            ipm[c] -= 1;
                            let im = idx;
                            let mut imm = im;
                            imm[c] -= 1;
                            let adv_p = 0.5 * (arr.get(ip[0], ip[1], ip[2]) + arr.get(ipm[0], ipm[1], ipm[2]));
                            let adv_m = 0.5 * (arr.get(im[0], im[1], im[2]) + arr.get(imm[0], imm[1], imm[2]));
                            conv += (adv_p * 0.5 * (qc + qp) - adv_m * 0.5 * (qm + qc)) * ih;
                        }
                    }
                    let mut cm = idx;
                    cm[c] -= 1;
                    let gp = (p[pidx(idx[0], idx[1], idx[2])] - p[pidx(cm[0], cm[1], cm[2])]) * ih;
                    line[(i + 1) as usize] = -conv + nu * lap * ih2 - gp;
                }
            });
            out
        })
    }

    /// Divergence at cell centres.
    fn divergence(&self, vel: &[Arr; 3]) -> Vec<f64> {
        let g = self.grid;
        let [nx, ny, nz] = g.n;
        let ih = 1.0 / g.h;
        (0..nx * ny * nz)
            .into_par_iter()
            .map(|q| {
                let (i, j, k) = ((q % nx) as isize, ((q / nx) % ny) as isize, (q / (nx * ny)) as isize);
                let mut s = 0.0;
                for c in 0..3 {
                    let mut hi = [i, j, k];
                    hi[c] += 1;
                    s += vel[c].get(hi[0], hi[1], hi[2]) - vel[c].get(i, j, k);
                }
                s * ih
            })
            .collect()
    }

    /// Leray projection: solve `L φ = ∇·w / dt`, `w −= dt ∇φ`; returns `φ`.
    /// With `incremental`, `p += φ`.
    fn project(&mut self, vel: &mut [Arr; 3], dt: f64, incremental: bool) -> Vec<f64> {
        let mut phi = self.divergence(vel);
        if self.mass_source && !self.ib_cell.is_empty() {
            // Net leak of the held cells: in a closed box (singular
            // operator) the free cells must absorb it, so it is spread over
            // them; this keeps the right-hand side compatible and lets the
            // fixed point have φ = 0. With an outlet the outlet absorbs it.
            let (leak, free) = phi
                .iter()
                .zip(&self.ib_cell)
                .fold((0.0, 0usize), |(s, n), (v, &m)| if m { (s + v, n) } else { (s, n + 1) });
            let spread = if self.poisson.singular() && free > 0 { leak / free as f64 } else { 0.0 };
            phi.par_iter_mut().zip(self.ib_cell.par_iter()).for_each(|(v, &m)| {
                if m {
                    *v = 0.0;
                } else {
                    *v += spread;
                }
            });
        }
        phi.iter_mut().for_each(|v| *v /= dt);
        self.poisson.solve(&mut phi);
        let g = self.grid;
        let n = g.n.map(|v| v as isize);
        let coef = dt / g.h;
        let pidx = |i: [isize; 3]| i[0] as usize + g.n[0] * (i[1] as usize + g.n[1] * i[2] as usize);
        let phi_ref = &phi;
        for c in 0..3 {
            if !self.active[c] {
                continue;
            }
            let arr = &mut vel[c];
            let d = arr.d.map(|v| v as isize);
            let st = arr.strides();
            let lo_out = matches!(self.bcs[c][0], Bc::Outlet);
            let hi_out = matches!(self.bcs[c][1], Bc::Outlet);
            let periodic = matches!(self.bcs[c][0], Bc::Periodic);
            let rows_per_plane = arr.d[1] + 2;
            arr.data.par_chunks_mut(st[1]).enumerate().for_each(|(row, line)| {
                let (j, k) = ((row % rows_per_plane) as isize - 1, (row / rows_per_plane) as isize - 1);
                if k < 0 || k >= d[2] || j < 0 || j >= d[1] {
                    return;
                }
                for i in 0..d[0] {
                    let idx = [i, j, k];
                    let ic = idx[c];
                    let mut lo = idx;
                    lo[c] -= 1;
                    let grad = if ic > 0 && ic < n[c] {
                        phi_ref[pidx(idx)] - phi_ref[pidx(lo)]
                    } else if ic == 0 {
                        if periodic {
                            let mut w = idx;
                            w[c] = n[c] - 1;
                            phi_ref[pidx(idx)] - phi_ref[pidx(w)]
                        } else if lo_out {
                            2.0 * phi_ref[pidx(idx)]
                        } else {
                            0.0
                        }
                    } else if periodic {
                        let mut w = idx;
                        w[c] = 0;
                        phi_ref[pidx(w)] - phi_ref[pidx(lo)]
                    } else if hi_out {
                        -2.0 * phi_ref[pidx(lo)]
                    } else {
                        0.0
                    };
                    if grad != 0.0 {
                        line[(i + 1) as usize] -= coef * grad;
                    }
                }
            });
        }
        if incremental {
            self.p.par_iter_mut().zip(phi.par_iter()).for_each(|(p, f)| *p += f);
        }
        phi
    }

    // ---------------------------------------------------------- stepping

    /// Time step: convective CFL and the RK3 diffusion limit (real-axis
    /// stability 2.51 for eigenvalue `−4 dim ν/h²` gives `0.63 h²/(dim ν)`;
    /// 0.5 is used).
    pub fn stable_dt(&self, cfl: f64) -> f64 {
        let umax = self.vel.iter().flat_map(|a| a.data.iter()).fold(0.0f64, |m, v| m.max(v.abs())).max(1e-12);
        let dim = self.grid.dim().max(1) as f64;
        let h = self.grid.h;
        (cfl * h / umax).min(0.5 * h * h / (dim * self.nu))
    }

    /// One SSP-RK3 step of size `dt`.
    pub fn step(&mut self, dt: f64) {
        let t0 = self.t;
        let u0 = self.vel.clone();
        // SSP-RK3: w ← a u⁰ + b (w + dt R(w)).
        let stages: [(f64, f64); 3] = [(0.0, 1.0), (0.75, 0.25), (1.0 / 3.0, 2.0 / 3.0)];
        let mut w = u0.clone();
        for (s, &(a, b)) in stages.iter().enumerate() {
            let r = self.rhs(&w);
            let mut next: [Arr; 3] = std::array::from_fn(|c| {
                let mut out = u0[c].clone();
                out.data
                    .par_iter_mut()
                    .zip(w[c].data.par_iter().zip(r[c].data.par_iter()))
                    .for_each(|(o, (wv, rv))| *o = a * *o + b * (wv + dt * rv));
                out
            });
            let ts = t0 + dt * [1.0, 0.5, 1.0][s];
            self.fill_bcs_on(&mut next, ts);
            self.apply_ib(&mut next);
            self.project(&mut next, b * dt, s == 2);
            self.fill_ghosts(&mut next, ts);
            w = next;
        }
        self.vel = w;
        self.t = t0 + dt;
    }

    /// Set the (kinematic) pressure from a function of position.
    pub fn initialize_pressure(&mut self, f: impl Fn(P3) -> f64) {
        let g = self.grid;
        for k in 0..g.n[2] {
            for j in 0..g.n[1] {
                for i in 0..g.n[0] {
                    self.p[i + g.n[0] * (j + g.n[1] * k)] = f(g.node(3, i as isize, j as isize, k as isize));
                }
            }
        }
    }

    /// Step until `max |Δu|/dt < tol` (an acceleration, velocity/time) or
    /// `t ≥ t_max`. Returns `(steps, final rate)`.
    pub fn run_until_steady(&mut self, cfl: f64, t_max: f64, tol: f64) -> (usize, f64) {
        let mut steps = 0;
        let mut rate = f64::INFINITY;
        while self.t < t_max - 1e-12 * t_max.abs().max(1.0) && rate >= tol {
            let dt = self.stable_dt(cfl).min(t_max - self.t);
            let prev: Vec<Vec<f64>> = self.vel.iter().map(|a| a.data.clone()).collect();
            self.step(dt);
            rate = self
                .vel
                .iter()
                .zip(&prev)
                .flat_map(|(a, b)| a.data.iter().zip(b).map(|(x, y)| (x - y).abs()))
                .fold(0.0, f64::max)
                / dt;
            steps += 1;
        }
        (steps, rate)
    }

    /// Momentum residual `R` at the current state (zero at a steady state
    /// away from the body).
    pub fn residual(&self) -> [Arr; 3] {
        self.rhs(&self.vel)
    }

    /// `max |∇·u|` over the cells that are not mass-source (immersed-boundary)
    /// cells.
    pub fn max_divergence(&self) -> f64 {
        self.divergence(&self.vel)
            .iter()
            .enumerate()
            .filter(|(q, _)| !self.ib_cell.get(*q).copied().unwrap_or(false))
            .fold(0.0, |m, (_, v)| m.max(v.abs()))
    }

    /// Net mass source of the immersed-boundary cells, `Σ ∇·u h³`, the
    /// volume flux that leaks through the immersed surfaces (it should
    /// vanish with refinement).
    pub fn ib_mass_source(&self) -> f64 {
        let h3 = self.grid.h.powi(3);
        self.divergence(&self.vel)
            .iter()
            .zip(&self.ib_cell)
            .filter(|(_, &m)| m)
            .map(|(v, _)| v * h3)
            .sum()
    }

    /// Force of the fluid on body `b` from the momentum balance of its
    /// immersed-boundary nodes. At a steady state the forcing that holds
    /// them is `f = −R`, so the fluid pushes the body with `ρ h³ Σ R` over
    /// the body's solid and forcing nodes (exact for the discrete
    /// equations; over the whole box depth).
    pub fn body_force_momentum(&self, b: usize) -> P3 {
        let r = self.residual();
        let vol = self.grid.h.powi(3);
        std::array::from_fn(|c| {
            let ib = &self.ib[c];
            let s: f64 = ib.solid.iter().filter(|e| e.2 == b).map(|e| r[c].data[e.0]).sum::<f64>()
                + ib.forcing.iter().filter(|f| f.body == b).map(|f| r[c].data[f.node]).sum::<f64>();
            self.rho * vol * s
        })
    }

    /// Torque of the fluid on body `b` about its `center`, from the same
    /// momentum balance: `ρ h³ Σ (x − c) × R`.
    pub fn body_torque_momentum(&self, b: usize) -> P3 {
        let r = self.residual();
        let vol = self.grid.h.powi(3);
        let center = self.bodies[b].center;
        let mut t = [0.0; 3];
        for c in 0..3 {
            let arr = &self.vel[c];
            let nodes = self.ib[c].solid.iter().filter(|e| e.2 == b).map(|e| e.0).chain(self.ib[c].forcing.iter().filter(|f| f.body == b).map(|f| f.node));
            for q in nodes {
                let [i, j, k] = arr.ijk(q);
                let x = self.grid.node(c, i, j, k);
                let rel = [x[0] - center[0], x[1] - center[1], x[2] - center[2]];
                let mut f = [0.0; 3];
                f[c] = self.rho * vol * r[c].data[q];
                t[0] += rel[1] * f[2] - rel[2] * f[1];
                t[1] += rel[2] * f[0] - rel[0] * f[2];
                t[2] += rel[0] * f[1] - rel[1] * f[0];
            }
        }
        t
    }

    // ---------------------------------------------------------- sampling

    /// Velocity at a point (trilinear on each staggered component).
    pub fn velocity_at(&self, x: P3) -> P3 {
        let g = self.grid;
        std::array::from_fn(|c| {
            let f: P3 = std::array::from_fn(|a| (x[a] - g.origin[a]) / g.h - if a == c { 0.0 } else { 0.5 });
            let arr = &self.vel[c];
            arr.interp(&arr.stencil(f))
        })
    }

    /// Kinematic pressure at a point (trilinear on cell centres, clamped).
    pub fn pressure_at(&self, x: P3) -> f64 {
        let g = self.grid;
        let mut base = [0usize; 3];
        let mut t = [0.0; 3];
        for a in 0..3 {
            let f = ((x[a] - g.origin[a]) / g.h - 0.5).clamp(0.0, (g.n[a] - 1) as f64);
            let b = (f.floor() as usize).min(g.n[a].saturating_sub(2));
            base[a] = b;
            t[a] = if g.n[a] == 1 { 0.0 } else { f - b as f64 };
        }
        let mut s = 0.0;
        for m in 0..8 {
            let o = [m & 1, (m >> 1) & 1, (m >> 2) & 1];
            if (0..3).any(|a| o[a] == 1 && g.n[a] == 1) {
                continue;
            }
            let w: f64 = (0..3).map(|a| if o[a] == 1 { t[a] } else { 1.0 - t[a] }).product();
            let idx: [usize; 3] = std::array::from_fn(|a| base[a] + o[a]);
            s += w * self.p[idx[0] + g.n[0] * (idx[1] + g.n[1] * idx[2])];
        }
        s
    }
}
