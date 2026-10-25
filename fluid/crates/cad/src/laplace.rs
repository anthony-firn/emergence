//! Two-dimensional electrostatics on a drawn cross-section.
//!
//! This is the one place in the crate that *solves* something rather than
//! generating geometry, and it is here for one reason: a drawn conductor
//! has an impedance, and the only way to know that the drawing has the
//! impedance it was meant to have is to solve for it. A PCB trace width is
//! chosen for 50 ohms; the generator can check its own work.
//!
//! The equation is `div(eps grad phi) = 0` with conductors held at fixed
//! potentials, discretised on a uniform grid with the permittivity taken
//! as the mean of the two cells either side of each face. In two dimensions
//! with square cells the grid spacing cancels, so the energy per unit
//! length is exactly
//!
//! ```text
//!   W = (eps0 / 2) sum_faces eps_r,face (phi_i - phi_j)^2
//! ```
//!
//! and the capacitance per unit length of a two-conductor system at unit
//! potential difference is `C = 2W`. For a quasi-TEM line the inductance
//! follows from the same geometry in vacuum, `L = 1 / (c^2 C_air)`, so
//!
//! ```text
//!   Z0 = sqrt(L/C) = 1 / (c sqrt(C C_air)),    eps_eff = C / C_air.
//! ```
//!
//! The linear system is solved by Jacobi-preconditioned conjugate
//! gradients. The discretisation is second order where the solution is
//! smooth, but a conductor corner carries a field singularity, so the
//! capacitance converges more slowly there; [`Solution::capacitance`]
//! should be checked at two resolutions, which is what the tests do.
//!
//! A problem with a symmetry plane may be solved on half the domain with
//! an [`Edge::Mirror`]; say so with [`Section::with_multiplicity`], or the
//! capacitance comes out halved and the impedance doubled.

use crate::sketch::{Region2, V2};

/// Vacuum permittivity, F/m (CODATA 2018, exact by definition of the
/// ampere to within the fine-structure constant).
pub const EPS0: f64 = 8.8541878128e-12;
/// Speed of light in vacuum, m/s (exact).
pub const C0: f64 = 299_792_458.0;

/// What happens at an outer edge of the domain.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Edge {
    /// `phi = 0`: a far-field ground.
    Grounded,
    /// No normal flux: a mirror plane, which halves a symmetric problem
    /// exactly.
    Mirror,
}

/// A cross-section to solve.
#[derive(Clone, Debug)]
pub struct Section {
    pub bounds: (V2, V2),
    pub cell: f64,
    /// Relative permittivity outside every listed dielectric (1 for air).
    pub eps_background: f64,
    /// Dielectric regions, searched in order; the first hit wins.
    pub dielectrics: Vec<(Region2, f64)>,
    /// Conductors and their potentials, in volts.
    pub conductors: Vec<(Region2, f64)>,
    /// `[x_min, x_max, y_min, y_max]`.
    pub edges: [Edge; 4],
    /// How many copies of the solved domain make up the whole problem.
    ///
    /// A mirror plane that cuts *through* a conductor halves the problem,
    /// so the energy and the capacitance of the whole line are twice what
    /// the solved domain holds: set `multiplicity = 2`. A mirror that only
    /// truncates empty space far from the conductors is a boundary
    /// condition, not a symmetry, and leaves this at 1.
    pub multiplicity: f64,
}

/// The solved field.
#[derive(Clone, Debug)]
pub struct Solution {
    pub dims: [usize; 2],
    pub cell: f64,
    pub potential: Vec<f64>,
    /// Energy per unit length, J/m, at the potentials as given.
    pub energy: f64,
    /// Conjugate-gradient iterations used, and the final relative
    /// residual.
    pub iterations: usize,
    pub residual: f64,
}

impl Section {
    pub fn new(bounds: (V2, V2), cell: f64) -> Self {
        Self {
            bounds,
            cell,
            eps_background: 1.0,
            dielectrics: Vec::new(),
            conductors: Vec::new(),
            edges: [Edge::Grounded; 4],
            multiplicity: 1.0,
        }
    }

    /// Set how many copies of this domain make the whole problem.
    pub fn with_multiplicity(mut self, multiplicity: f64) -> Self {
        self.multiplicity = multiplicity;
        self
    }

    pub fn with_dielectric(mut self, region: Region2, eps_r: f64) -> Self {
        self.dielectrics.push((region, eps_r));
        self
    }

    pub fn with_conductor(mut self, region: Region2, volts: f64) -> Self {
        self.conductors.push((region, volts));
        self
    }

    pub fn with_edges(mut self, edges: [Edge; 4]) -> Self {
        self.edges = edges;
        self
    }

    /// Replace every dielectric with vacuum, keeping the conductors: the
    /// companion problem whose capacitance gives the line inductance.
    pub fn in_vacuum(&self) -> Self {
        Self { eps_background: 1.0, dielectrics: Vec::new(), ..self.clone() }
    }

    fn dims(&self) -> [usize; 2] {
        std::array::from_fn(|k| (((self.bounds.1[k] - self.bounds.0[k]) / self.cell).round() as usize).max(2) + 1)
    }

    fn point(&self, i: usize, j: usize) -> V2 {
        [self.bounds.0[0] + self.cell * i as f64, self.bounds.0[1] + self.cell * j as f64]
    }

    /// Solve to a relative residual of `tol` or `max_iter` iterations.
    pub fn solve(&self, tol: f64, max_iter: usize) -> Solution {
        let [nx, ny] = self.dims();
        let n = nx * ny;
        let id = |i: usize, j: usize| i + nx * j;

        // Cell-centred permittivity, and node-fixed potentials.
        let mut eps = vec![self.eps_background; (nx - 1) * (ny - 1)];
        for j in 0..ny - 1 {
            for i in 0..nx - 1 {
                let p = [
                    self.bounds.0[0] + self.cell * (i as f64 + 0.5),
                    self.bounds.0[1] + self.cell * (j as f64 + 0.5),
                ];
                for (region, e) in &self.dielectrics {
                    if region.sdf(p) < 0.0 {
                        eps[i + (nx - 1) * j] = *e;
                        break;
                    }
                }
            }
        }
        let mut fixed = vec![f64::NAN; n];
        // A node on a conductor's boundary belongs to it. The field of a
        // region evaluated exactly on its edge can come out +1e-20, and
        // testing `<= 0` then dropped a whole row of nodes: a trace drawn
        // on a dielectric surface floated one cell above it in air.
        let on_edge = 1e-9 * self.cell;
        for j in 0..ny {
            for i in 0..nx {
                let p = self.point(i, j);
                for (region, v) in &self.conductors {
                    if region.sdf(p) <= on_edge {
                        fixed[id(i, j)] = *v;
                        break;
                    }
                }
            }
        }
        for j in 0..ny {
            for i in 0..nx {
                let edge = [
                    i == 0 && self.edges[0] == Edge::Grounded,
                    i == nx - 1 && self.edges[1] == Edge::Grounded,
                    j == 0 && self.edges[2] == Edge::Grounded,
                    j == ny - 1 && self.edges[3] == Edge::Grounded,
                ];
                if edge.iter().any(|&e| e) && fixed[id(i, j)].is_nan() {
                    fixed[id(i, j)] = 0.0;
                }
            }
        }

        // Face permittivity: the mean of the cells touching the face.
        let cell_eps = |ci: isize, cj: isize| -> Option<f64> {
            if ci < 0 || cj < 0 || ci >= (nx - 1) as isize || cj >= (ny - 1) as isize {
                None
            } else {
                Some(eps[ci as usize + (nx - 1) * cj as usize])
            }
        };
        // A face's area fraction inside the domain: a node row or column
        // on a mirror plane owns only half a cell, so the faces normal to
        // the mirror carry half their area. Getting this right is what
        // makes a mirrored problem exactly equal to the whole one.
        let ax = |i: usize| -> f64 {
            if (i == 0 && self.edges[0] == Edge::Mirror) || (i == nx - 1 && self.edges[1] == Edge::Mirror) {
                0.5
            } else {
                1.0
            }
        };
        let ay = |j: usize| -> f64 {
            if (j == 0 && self.edges[2] == Edge::Mirror) || (j == ny - 1 && self.edges[3] == Edge::Mirror) {
                0.5
            } else {
                1.0
            }
        };
        // Face between (i,j) and (i+1,j): cells (i, j-1) and (i, j), with
        // the area fraction of row j.
        let eps_x = |i: usize, j: usize| -> f64 {
            let a = cell_eps(i as isize, j as isize - 1);
            let b = cell_eps(i as isize, j as isize);
            mean(a, b) * ay(j)
        };
        let eps_y = |i: usize, j: usize| -> f64 {
            let a = cell_eps(i as isize - 1, j as isize);
            let b = cell_eps(i as isize, j as isize);
            mean(a, b) * ax(i)
        };

        // Assemble the free-node operator as neighbour weights.
        let mut diag = vec![0.0f64; n];
        let mut weights: Vec<[f64; 4]> = vec![[0.0; 4]; n]; // -x, +x, -y, +y
        for j in 0..ny {
            for i in 0..nx {
                if !fixed[id(i, j)].is_nan() {
                    continue;
                }
                let mut w = [0.0; 4];
                if i > 0 {
                    w[0] = eps_x(i - 1, j);
                }
                if i + 1 < nx {
                    w[1] = eps_x(i, j);
                }
                if j > 0 {
                    w[2] = eps_y(i, j - 1);
                }
                if j + 1 < ny {
                    w[3] = eps_y(i, j);
                }
                weights[id(i, j)] = w;
                diag[id(i, j)] = w.iter().sum();
            }
        }

        // Right-hand side from the fixed nodes.
        let neighbour = |i: usize, j: usize, k: usize| -> Option<(usize, usize)> {
            match k {
                0 if i > 0 => Some((i - 1, j)),
                1 if i + 1 < nx => Some((i + 1, j)),
                2 if j > 0 => Some((i, j - 1)),
                3 if j + 1 < ny => Some((i, j + 1)),
                _ => None,
            }
        };
        let mut b = vec![0.0f64; n];
        for j in 0..ny {
            for i in 0..nx {
                let q = id(i, j);
                if !fixed[q].is_nan() {
                    continue;
                }
                for k in 0..4 {
                    if let Some((ii, jj)) = neighbour(i, j, k) {
                        let v = fixed[id(ii, jj)];
                        if !v.is_nan() {
                            b[q] += weights[q][k] * v;
                        }
                    }
                }
            }
        }

        // A x, for free nodes only.
        let apply = |x: &[f64], out: &mut [f64]| {
            for j in 0..ny {
                for i in 0..nx {
                    let q = id(i, j);
                    if !fixed[q].is_nan() {
                        out[q] = 0.0;
                        continue;
                    }
                    let mut s = diag[q] * x[q];
                    for k in 0..4 {
                        if let Some((ii, jj)) = neighbour(i, j, k) {
                            let p = id(ii, jj);
                            if fixed[p].is_nan() {
                                s -= weights[q][k] * x[p];
                            }
                        }
                    }
                    out[q] = s;
                }
            }
        };

        // Jacobi-preconditioned conjugate gradients.
        let mut x = vec![0.0f64; n];
        let mut r = b.clone();
        let mut z: Vec<f64> = r.iter().enumerate().map(|(q, &v)| if diag[q] > 0.0 { v / diag[q] } else { 0.0 }).collect();
        let mut p = z.clone();
        let mut ap = vec![0.0f64; n];
        let b_norm = dot(&b, &b).sqrt().max(1e-300);
        let mut rz = dot(&r, &z);
        let mut iterations = 0;
        let mut residual = 1.0;
        for it in 0..max_iter {
            apply(&p, &mut ap);
            let denom = dot(&p, &ap);
            if denom.abs() < 1e-300 {
                break;
            }
            let alpha = rz / denom;
            for q in 0..n {
                x[q] += alpha * p[q];
                r[q] -= alpha * ap[q];
            }
            residual = dot(&r, &r).sqrt() / b_norm;
            iterations = it + 1;
            if residual < tol {
                break;
            }
            for q in 0..n {
                z[q] = if diag[q] > 0.0 { r[q] / diag[q] } else { 0.0 };
            }
            let rz_new = dot(&r, &z);
            let beta = rz_new / rz;
            rz = rz_new;
            for q in 0..n {
                p[q] = z[q] + beta * p[q];
            }
        }

        // The full potential, fixed nodes included.
        let potential: Vec<f64> = (0..n).map(|q| if fixed[q].is_nan() { x[q] } else { fixed[q] }).collect();

        // Energy per unit length: the face sum, with the grid spacing
        // cancelling in 2D.
        let mut energy = 0.0;
        for j in 0..ny {
            for i in 0..nx {
                if i + 1 < nx {
                    let d = potential[id(i + 1, j)] - potential[id(i, j)];
                    energy += eps_x(i, j) * d * d;
                }
                if j + 1 < ny {
                    let d = potential[id(i, j + 1)] - potential[id(i, j)];
                    energy += eps_y(i, j) * d * d;
                }
            }
        }
        energy *= 0.5 * EPS0 * self.multiplicity;

        Solution {
            dims: [nx, ny],
            cell: self.cell,
            potential,
            energy,
            iterations,
            residual,
        }
    }
}

fn mean(a: Option<f64>, b: Option<f64>) -> f64 {
    match (a, b) {
        (Some(x), Some(y)) => 0.5 * (x + y),
        (Some(x), None) | (None, Some(x)) => x,
        (None, None) => 0.0,
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

impl Solution {
    /// Capacitance per unit length, F/m, for a two-conductor section at a
    /// potential difference `volts`.
    pub fn capacitance(&self, volts: f64) -> f64 {
        2.0 * self.energy / (volts * volts)
    }
}

/// Characteristic impedance and effective permittivity of a quasi-TEM line
/// from two solves of the same geometry, one with the dielectrics and one
/// in vacuum.
pub fn quasi_tem(section: &Section, volts: f64, tol: f64, max_iter: usize) -> (f64, f64, Solution) {
    let filled = section.solve(tol, max_iter);
    let vacuum = section.in_vacuum().solve(tol, max_iter);
    let c = filled.capacitance(volts);
    let c_air = vacuum.capacitance(volts);
    let z0 = 1.0 / (C0 * (c * c_air).sqrt());
    (z0, c / c_air, filled)
}
