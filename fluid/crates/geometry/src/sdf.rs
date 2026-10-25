//! Signed distance to a triangle surface.
//!
//! The distance is exact (BVH closest point). The sign comes from the
//! angle-weighted pseudonormal of the closest feature (Bærentzen & Aanæs,
//! IEEE TVCG 11(3), 2005), which is provably correct for a closed,
//! consistently oriented, manifold surface. For a surface that is not, the
//! sign comes from the generalised winding number (Jacobson, Kavan &
//! Sorkine-Hornung, SIGGRAPH 2013): the exact sum of signed solid angles,
//! robust to holes and self-overlap but `O(triangles)` per query.

use rayon::prelude::*;

use crate::bvh::{Bvh, Feature};
use crate::mesh::TriMesh;
use crate::{add, cross, dot, norm, normalize, scale, sub, P3};

/// How the inside/outside sign is decided.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SignMode {
    PseudoNormal,
    WindingNumber,
}

/// Signed distance field of a surface (negative inside).
#[derive(Clone, Debug)]
pub struct Sdf {
    pub mesh: TriMesh,
    bvh: Bvh,
    face_n: Vec<P3>,
    /// Pseudonormal of each triangle's three edges (sum of the two faces).
    edge_n: Vec<[P3; 3]>,
    vertex_n: Vec<P3>,
    pub sign: SignMode,
}

/// Closest-point information.
#[derive(Clone, Copy, Debug)]
pub struct Closest {
    /// Signed distance (negative inside).
    pub distance: f64,
    pub point: P3,
    /// Outward unit normal of the surface at the closest point (for `p` off
    /// the surface: the unit vector from the surface towards `p`, flipped
    /// inside so it always points out of the solid).
    pub normal: P3,
    pub triangle: u32,
}

impl Sdf {
    /// Build from a mesh. Welds coincident vertices first; uses
    /// pseudonormals when the result is watertight, winding numbers
    /// otherwise.
    pub fn new(mesh: TriMesh) -> Self {
        Self::build(mesh, true)
    }

    /// As [`Sdf::new`] but keeps the given orientation: a closed surface
    /// with inward normals describes the *exterior* as solid (e.g. a pipe or
    /// housing wall around the fluid).
    pub fn new_as_oriented(mesh: TriMesh) -> Self {
        Self::build(mesh, false)
    }

    fn build(mut mesh: TriMesh, orient: bool) -> Self {
        mesh.weld_relative(1e-9);
        let report = mesh.report();
        if report.watertight && orient {
            mesh.orient_outward();
        }
        let sign = if report.watertight { SignMode::PseudoNormal } else { SignMode::WindingNumber };
        let bvh = Bvh::build(&mesh);
        let face_n: Vec<P3> = (0..mesh.triangles.len()).map(|t| mesh.unit_normal(t)).collect();
        let mut vertex_n = vec![[0.0; 3]; mesh.vertices.len()];
        for (t, tri) in mesh.triangles.iter().enumerate() {
            let c = mesh.corners(t);
            for k in 0..3 {
                let e1 = normalize(sub(c[(k + 1) % 3], c[k]));
                let e2 = normalize(sub(c[(k + 2) % 3], c[k]));
                let angle = dot(e1, e2).clamp(-1.0, 1.0).acos();
                vertex_n[tri[k] as usize] = add(vertex_n[tri[k] as usize], scale(face_n[t], angle));
            }
        }
        let mut edge_n = vec![[[0.0; 3]; 3]; mesh.triangles.len()];
        for ((a, b), uses) in mesh.edge_map() {
            let sum = uses.iter().fold([0.0; 3], |s, &(t, _)| add(s, face_n[t as usize]));
            for &(t, _) in &uses {
                let tri = mesh.triangles[t as usize];
                for k in 0..3 {
                    let (x, y) = (tri[k], tri[(k + 1) % 3]);
                    if (x.min(y), x.max(y)) == (a, b) {
                        edge_n[t as usize][k] = sum;
                    }
                }
            }
        }
        Self {
            mesh,
            bvh,
            face_n,
            edge_n,
            vertex_n,
            sign,
        }
    }

    /// Closest point, signed distance and outward normal.
    pub fn closest(&self, p: P3) -> Closest {
        let hit = self.bvh.closest(&self.mesh, p).expect("empty mesh");
        let t = hit.triangle as usize;
        let feature_n = match hit.feature {
            Feature::Face => self.face_n[t],
            Feature::Edge(k) => self.edge_n[t][k as usize],
            Feature::Vertex(k) => self.vertex_n[self.mesh.triangles[t][k as usize] as usize],
        };
        let d = hit.dist2.sqrt();
        let inside = match self.sign {
            SignMode::PseudoNormal => dot(sub(p, hit.point), feature_n) < 0.0,
            SignMode::WindingNumber => self.winding_number(p) > 0.5,
        };
        let offset = sub(p, hit.point);
        let normal = if d > 1e-12 * (1.0 + norm(p)) {
            scale(offset, if inside { -1.0 / d } else { 1.0 / d })
        } else {
            normalize(feature_n)
        };
        Closest {
            distance: if inside { -d } else { d },
            point: hit.point,
            normal,
            triangle: hit.triangle,
        }
    }

    pub fn distance(&self, p: P3) -> f64 {
        self.closest(p).distance
    }

    /// Generalised winding number: 1 inside a closed outward surface, 0
    /// outside, fractional near holes. Exact solid angles (Van Oosterom &
    /// Strackee, IEEE TBME 30(2), 1983).
    pub fn winding_number(&self, p: P3) -> f64 {
        let mut w = 0.0;
        for t in 0..self.mesh.triangles.len() {
            let [a, b, c] = self.mesh.corners(t);
            let (a, b, c) = (sub(a, p), sub(b, p), sub(c, p));
            let (la, lb, lc) = (norm(a), norm(b), norm(c));
            let num = dot(a, cross(b, c));
            let den = la * lb * lc + dot(a, b) * lc + dot(b, c) * la + dot(c, a) * lb;
            w += 2.0 * num.atan2(den);
        }
        w / (4.0 * std::f64::consts::PI)
    }

    /// Signed distances at many points, in parallel.
    pub fn sample(&self, points: &[P3]) -> Vec<f64> {
        points.par_iter().map(|&p| self.distance(p)).collect()
    }

    /// Signed distance on the nodes of a regular grid, `x` fastest:
    /// node `(i, j, k)` at `origin + h·(i, j, k)`.
    pub fn sample_grid(&self, origin: P3, h: P3, dims: [usize; 3]) -> Vec<f64> {
        (0..dims[0] * dims[1] * dims[2])
            .into_par_iter()
            .map(|q| {
                let (i, j, k) = (q % dims[0], (q / dims[0]) % dims[1], q / (dims[0] * dims[1]));
                self.distance([origin[0] + h[0] * i as f64, origin[1] + h[1] * j as f64, origin[2] + h[2] * k as f64])
            })
            .collect()
    }

    /// Face normal of triangle `t` (outward for a watertight mesh).
    pub fn face_normal(&self, t: usize) -> P3 {
        self.face_n[t]
    }
}
