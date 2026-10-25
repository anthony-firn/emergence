//! Post-processing on the CAD surface and the grid: surface forces, wall
//! pressure, and ParaView (VTK) output.
//!
//! Wall quantities are sampled along the outward surface normal at
//! `s₁ = √dim · h` and `s₂ = 2 s₁` (outside the immersed-boundary band):
//! the wall pressure is extrapolated linearly, `p_w = 2p(s₁) − p(s₂)`, and
//! the wall-normal velocity derivative comes from the quadratic through the
//! wall velocity and the two samples. The traction on a rigid no-slip wall
//! is `−p n + μ ∂u/∂n` (the `(∇u)ᵀ n` term vanishes there).

use geometry::{cross, norm, sub, vtk, GeomError, P3};

use crate::solver::Solver;

/// Force of the fluid on a body by surface integration.
#[derive(Clone, Copy, Debug, Default)]
pub struct SurfaceForce {
    pub pressure: P3,
    pub viscous: P3,
    /// Integrated (wetted) area inside the box.
    pub area: f64,
}

impl SurfaceForce {
    pub fn total(&self) -> P3 {
        std::array::from_fn(|a| self.pressure[a] + self.viscous[a])
    }
}

/// Wall sample at a surface point with outward normal `n`.
#[derive(Clone, Copy, Debug)]
pub struct WallSample {
    /// Pressure (not kinematic): `ρ p_w`.
    pub pressure: f64,
    /// Viscous traction `μ ∂u/∂n`.
    pub shear: P3,
}

impl Solver {
    fn inside_box(&self, x: P3) -> bool {
        let (lo, hi) = (self.grid.origin, self.grid.hi());
        (0..3).all(|a| x[a] >= lo[a] && x[a] <= hi[a])
    }

    /// Wall pressure and shear at surface point `x` (outward normal `n`) of
    /// body `b`.
    pub fn wall_sample(&self, b: usize, x: P3, n: P3) -> WallSample {
        let s1 = self.probe;
        let s2 = 2.0 * s1;
        let at = |s: f64| std::array::from_fn(|a| x[a] + s * n[a]);
        let (x1, x2): (P3, P3) = (at(s1), at(s2));
        let pw = 2.0 * self.pressure_at(x1) - self.pressure_at(x2);
        let uw = self.bodies[b].velocity_at(x);
        let (u1, u2) = (self.velocity_at(x1), self.velocity_at(x2));
        let shear = std::array::from_fn(|a| {
            let d = ((u1[a] - uw[a]) * s2 * s2 - (u2[a] - uw[a]) * s1 * s1) / (s1 * s2 * (s2 - s1));
            self.rho * self.nu * d
        });
        WallSample { pressure: self.rho * pw, shear }
    }

    /// Wall pressure at the surface point of body `b` nearest to `x`.
    pub fn wall_pressure(&self, b: usize, x: P3) -> f64 {
        let cl = self.bodies[b].sdf.closest(x);
        self.wall_sample(b, cl.point, cl.normal).pressure
    }

    /// Force of the fluid on body `b` by integrating the traction over its
    /// surface triangles, clipped to the box (so an extruded 2D body gives
    /// the force over the box depth) and split into sub-triangles of edge
    /// `≤ h/2`, one sample at each sub-triangle centroid.
    pub fn surface_force(&self, b: usize) -> SurfaceForce {
        use rayon::prelude::*;
        let mesh = &self.bodies[b].sdf.mesh;
        let h = self.grid.h;
        let (lo, hi) = (self.grid.origin, self.grid.hi());
        let parts: Vec<SurfaceForce> = (0..mesh.triangles.len())
            .into_par_iter()
            .map(|t| {
                let n = mesh.unit_normal(t);
                let poly = clip_to_box(&mesh.corners(t), lo, hi);
                let mut f = SurfaceForce::default();
                for k in 1..poly.len().saturating_sub(1) {
                    let (a, bb, c) = (poly[0], poly[k], poly[k + 1]);
                    let area = 0.5 * norm(cross(sub(bb, a), sub(c, a)));
                    if area == 0.0 {
                        continue;
                    }
                    let longest = norm(sub(a, bb)).max(norm(sub(bb, c))).max(norm(sub(c, a)));
                    let m = ((longest / (0.5 * h)).ceil() as usize).max(1);
                    let da = area / (m * m) as f64;
                    // Centroids of the m² uniform sub-triangles.
                    for i in 0..m {
                        for j in 0..m - i {
                            for up in [false, true] {
                                if up && i + j + 1 >= m {
                                    continue;
                                }
                                let o = if up { 2.0 / 3.0 } else { 1.0 / 3.0 };
                                let (u, v) = ((i as f64 + o) / m as f64, (j as f64 + o) / m as f64);
                                let x: P3 = std::array::from_fn(|q| a[q] + u * (bb[q] - a[q]) + v * (c[q] - a[q]));
                                let w = self.wall_sample(b, x, n);
                                for q in 0..3 {
                                    f.pressure[q] -= w.pressure * n[q] * da;
                                    f.viscous[q] += w.shear[q] * da;
                                }
                                f.area += da;
                            }
                        }
                    }
                }
                f
            })
            .collect();
        parts.iter().fold(SurfaceForce::default(), |mut acc, f| {
            for k in 0..3 {
                acc.pressure[k] += f.pressure[k];
                acc.viscous[k] += f.viscous[k];
            }
            acc.area += f.area;
            acc
        })
    }

    /// Cell-centred fields for ParaView: velocity, pressure, vorticity and
    /// the signed distance to the nearest body.
    pub fn write_vtk(&self, path: impl AsRef<std::path::Path>) -> Result<(), GeomError> {
        let g = self.grid;
        let [nx, ny, nz] = g.n;
        let n = nx * ny * nz;
        let mut vel = vec![[0.0; 3]; n];
        let mut p = vec![0.0; n];
        let mut vort = vec![[0.0; 3]; n];
        let mut dist = vec![f64::INFINITY; n];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let q = i + nx * (j + ny * k);
                    let (ii, jj, kk) = (i as isize, j as isize, k as isize);
                    for c in 0..3 {
                        let mut hi = [ii, jj, kk];
                        hi[c] += 1;
                        vel[q][c] = 0.5 * (self.vel[c].get(ii, jj, kk) + self.vel[c].get(hi[0], hi[1], hi[2]));
                    }
                    p[q] = self.rho * self.p[q];
                    let x = g.node(3, ii, jj, kk);
                    let d = |a: usize, c: usize| {
                        let mut xp = x;
                        let mut xm = x;
                        xp[a] += 0.5 * g.h;
                        xm[a] -= 0.5 * g.h;
                        if g.n[a] == 1 {
                            return 0.0;
                        }
                        (self.velocity_at(xp)[c] - self.velocity_at(xm)[c]) / g.h
                    };
                    vort[q] = [d(1, 2) - d(2, 1), d(2, 0) - d(0, 2), d(0, 1) - d(1, 0)];
                    for body in &self.bodies {
                        dist[q] = dist[q].min(body.sdf.distance(x));
                    }
                }
            }
        }
        vtk::write_grid(
            path,
            g.origin,
            [g.h; 3],
            g.n,
            &[
                ("velocity", vtk::Data::Vectors(&vel)),
                ("pressure", vtk::Data::Scalars(&p)),
                ("vorticity", vtk::Data::Vectors(&vort)),
                ("distance", vtk::Data::Scalars(&dist)),
            ],
            &[],
        )
    }

    /// Body surface for ParaView: per-triangle wall pressure, pressure
    /// coefficient `(p − p_ref)/(½ ρ U²)` and wall shear (triangles whose
    /// centroid is outside the box get NaN).
    pub fn write_surface_vtk(&self, b: usize, path: impl AsRef<std::path::Path>, p_ref: f64, u_ref: f64) -> Result<(), GeomError> {
        let mesh = &self.bodies[b].sdf.mesh;
        let mut p = Vec::with_capacity(mesh.triangles.len());
        let mut cp = Vec::with_capacity(mesh.triangles.len());
        let mut tau = Vec::with_capacity(mesh.triangles.len());
        for t in 0..mesh.triangles.len() {
            let [a, bb, c] = mesh.corners(t);
            let x: P3 = std::array::from_fn(|k| (a[k] + bb[k] + c[k]) / 3.0);
            if self.inside_box(x) {
                let w = self.wall_sample(b, x, mesh.unit_normal(t));
                p.push(w.pressure);
                cp.push((w.pressure - p_ref) / (0.5 * self.rho * u_ref * u_ref));
                tau.push(w.shear);
            } else {
                p.push(f64::NAN);
                cp.push(f64::NAN);
                tau.push([f64::NAN; 3]);
            }
        }
        vtk::write_surface(path, mesh, &[], &[("pressure", vtk::Data::Scalars(&p)), ("Cp", vtk::Data::Scalars(&cp)), ("wall_shear", vtk::Data::Vectors(&tau))])
    }
}

/// Clip a triangle to the box `[lo, hi]` (Sutherland–Hodgman), returning
/// the convex polygon inside.
pub fn clip_to_box(tri: &[P3; 3], lo: P3, hi: P3) -> Vec<P3> {
    let mut poly: Vec<P3> = tri.to_vec();
    for a in 0..3 {
        for (bound, keep_above) in [(lo[a], true), (hi[a], false)] {
            if poly.is_empty() {
                return poly;
            }
            let inside = |p: &P3| if keep_above { p[a] >= bound } else { p[a] <= bound };
            let mut out = Vec::with_capacity(poly.len() + 2);
            for i in 0..poly.len() {
                let (p, q) = (poly[i], poly[(i + 1) % poly.len()]);
                let (ip, iq) = (inside(&p), inside(&q));
                if ip {
                    out.push(p);
                }
                if ip != iq {
                    let t = (bound - p[a]) / (q[a] - p[a]);
                    out.push(std::array::from_fn(|k| p[k] + t * (q[k] - p[k])));
                }
            }
            poly = out;
        }
    }
    poly
}
