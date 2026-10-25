//! Bounding-volume hierarchy over triangles: exact closest-point queries.

use crate::mesh::TriMesh;
use crate::{dot, sub, P3};

/// Which part of a triangle a closest point lies on (for pseudonormals).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Feature {
    /// Corner `k` of the triangle.
    Vertex(u8),
    /// Edge from corner `k` to corner `(k + 1) % 3`.
    Edge(u8),
    Face,
}

/// Closest point of triangle `abc` to `p` (Ericson, *Real-Time Collision
/// Detection*, 2005, §5.1.5), with the feature it lies on.
pub fn closest_on_triangle(p: P3, a: P3, b: P3, c: P3) -> (P3, Feature) {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return (a, Feature::Vertex(0));
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return (b, Feature::Vertex(1));
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return (lerp(a, ab, v), Feature::Edge(0));
    }
    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return (c, Feature::Vertex(2));
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return (lerp(a, ac, w), Feature::Edge(2));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (lerp(b, sub(c, b), w), Feature::Edge(1));
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    ([a[0] + ab[0] * v + ac[0] * w, a[1] + ab[1] * v + ac[1] * w, a[2] + ab[2] * v + ac[2] * w], Feature::Face)
}

fn lerp(a: P3, d: P3, t: f64) -> P3 {
    [a[0] + d[0] * t, a[1] + d[1] * t, a[2] + d[2] * t]
}

#[derive(Clone, Copy, Debug)]
struct Node {
    lo: P3,
    hi: P3,
    /// Leaf: first index into `order`; inner: index of the left child (the
    /// right child follows the whole left subtree, stored in `right`).
    first: u32,
    count: u32,
    right: u32,
}

/// A median-split BVH over a mesh's triangles.
#[derive(Clone, Debug)]
pub struct Bvh {
    nodes: Vec<Node>,
    order: Vec<u32>,
}

/// Result of a closest-point query.
#[derive(Clone, Copy, Debug)]
pub struct Hit {
    pub dist2: f64,
    pub point: P3,
    pub triangle: u32,
    pub feature: Feature,
}

/// Nearest ray crossing, with barycentric weights in vertex order.
#[derive(Clone, Copy, Debug)]
pub struct RayHit {
    pub distance: f64,
    pub triangle: u32,
    pub barycentric: [f64; 3],
}

fn ray_box(origin: P3, direction: P3, lo: P3, hi: P3, min: f64, max: f64) -> bool {
    let (mut a,mut b)=(min,max);
    for k in 0..3 {
        if direction[k]==0. {
            if origin[k]<lo[k]||origin[k]>hi[k] {return false;}
        } else {
            let x=(lo[k]-origin[k])/direction[k];
            let y=(hi[k]-origin[k])/direction[k];
            a=a.max(x.min(y)); b=b.min(x.max(y));
            if a>b {return false;}
        }
    }
    true
}

/// Sheared edge-function intersection, after PBRT 4e, Shapes/Triangle_Meshes.
/// Shared edges use the same transformed vertices, avoiding cracks between
/// neighbouring triangles. No barycentric epsilon expands the triangles.
fn ray_triangle(origin: P3,direction: P3,vertices: [P3;3],min: f64,max: f64) -> Option<(f64,[f64;3])> {
    let kz=(0..3).max_by(|&a,&b|direction[a].abs().total_cmp(&direction[b].abs()))?;
    if direction[kz]==0. {return None;}
    let (mut kx,mut ky)=((kz+1)%3,(kz+2)%3);
    if direction[kz]<0. {std::mem::swap(&mut kx,&mut ky);}
    let sx=-direction[kx]/direction[kz]; let sy=-direction[ky]/direction[kz];
    let p=vertices.map(|v|sub(v,origin));
    let x=p.map(|v|v[kx]+sx*v[kz]); let y=p.map(|v|v[ky]+sy*v[kz]);
    let e=[x[1]*y[2]-y[1]*x[2],x[2]*y[0]-y[2]*x[0],x[0]*y[1]-y[0]*x[1]];
    if e.iter().any(|v|*v<0.)&&e.iter().any(|v|*v>0.) {return None;}
    let det=e.iter().sum::<f64>();
    if det==0. {return None;}
    let t=(e[0]*p[0][kz]+e[1]*p[1][kz]+e[2]*p[2][kz])/(det*direction[kz]);
    if t<=min||t>=max {return None;}
    Some((t,e.map(|v|v/det)))
}

const LEAF: usize = 4;

fn dist2_box(p: P3, lo: P3, hi: P3) -> f64 {
    let mut d = 0.0;
    for k in 0..3 {
        let e = if p[k] < lo[k] {
            lo[k] - p[k]
        } else if p[k] > hi[k] {
            p[k] - hi[k]
        } else {
            0.0
        };
        d += e * e;
    }
    d
}

impl Bvh {
    /// Closest two-sided triangle crossing. `direction` need not be unit;
    /// distance is the parameter in `origin + distance * direction`.
    pub fn ray(&self,mesh:&TriMesh,origin:P3,direction:P3,min:f64,max:f64)->Option<RayHit> {
        if self.nodes.is_empty() {return None;}
        let mut best=None;
        let mut limit=max;
        // Median splits bound tree depth by log2(number of triangles).
        // u32 triangle indices imply fewer than 33 simultaneously pending nodes.
        let mut stack=[0_usize;64]; let mut len=1;
        while len>0 {
            len-=1; let node=self.nodes[stack[len]];
            if !ray_box(origin,direction,node.lo,node.hi,min,limit) {continue;}
            if node.count==0 {
                stack[len]=node.first as usize; stack[len+1]=node.right as usize; len+=2;
            } else {
                for &tri in &self.order[node.first as usize..(node.first+node.count) as usize] {
                    if let Some((distance,barycentric))=ray_triangle(origin,direction,mesh.corners(tri as usize),min,limit) {
                        limit=distance; best=Some(RayHit{distance,triangle:tri,barycentric});
                    }
                }
            }
        }
        best
    }
    pub fn build(mesh: &TriMesh) -> Self {
        let n = mesh.triangles.len();
        let mut order: Vec<u32> = (0..n as u32).collect();
        let cent: Vec<P3> = (0..n)
            .map(|t| {
                let [a, b, c] = mesh.corners(t);
                [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, (a[2] + b[2] + c[2]) / 3.0]
            })
            .collect();
        let mut nodes = Vec::with_capacity(2 * n / LEAF + 1);
        if n > 0 {
            Self::split(mesh, &cent, &mut order, 0, n, &mut nodes);
        }
        Self { nodes, order }
    }

    fn split(mesh: &TriMesh, cent: &[P3], order: &mut [u32], start: usize, end: usize, nodes: &mut Vec<Node>) -> usize {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        let mut clo = [f64::INFINITY; 3];
        let mut chi = [f64::NEG_INFINITY; 3];
        for &t in &order[start..end] {
            for v in mesh.corners(t as usize) {
                for k in 0..3 {
                    lo[k] = lo[k].min(v[k]);
                    hi[k] = hi[k].max(v[k]);
                }
            }
            for k in 0..3 {
                clo[k] = clo[k].min(cent[t as usize][k]);
                chi[k] = chi[k].max(cent[t as usize][k]);
            }
        }
        let id = nodes.len();
        nodes.push(Node {
            lo,
            hi,
            first: start as u32,
            count: (end - start) as u32,
            right: 0,
        });
        if end - start <= LEAF {
            return id;
        }
        let axis = (0..3).max_by(|&a, &b| (chi[a] - clo[a]).total_cmp(&(chi[b] - clo[b]))).unwrap();
        let mid = (start + end) / 2;
        order[start..end].select_nth_unstable_by(mid - start, |&a, &b| cent[a as usize][axis].total_cmp(&cent[b as usize][axis]));
        Self::split(mesh, cent, order, start, mid, nodes);
        let right = Self::split(mesh, cent, order, mid, end, nodes);
        nodes[id].count = 0;
        nodes[id].first = (id + 1) as u32;
        nodes[id].right = right as u32;
        id
    }

    /// Exact closest point on the mesh to `p`.
    pub fn closest(&self, mesh: &TriMesh, p: P3) -> Option<Hit> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut best = Hit {
            dist2: f64::INFINITY,
            point: [0.0; 3],
            triangle: 0,
            feature: Feature::Face,
        };
        let mut stack = Vec::with_capacity(64);
        stack.push(0usize);
        while let Some(i) = stack.pop() {
            let nd = &self.nodes[i];
            if dist2_box(p, nd.lo, nd.hi) >= best.dist2 {
                continue;
            }
            if nd.count > 0 {
                for &t in &self.order[nd.first as usize..(nd.first + nd.count) as usize] {
                    let [a, b, c] = mesh.corners(t as usize);
                    let (q, f) = closest_on_triangle(p, a, b, c);
                    let d = sub(p, q);
                    let d2 = dot(d, d);
                    if d2 < best.dist2 {
                        best = Hit {
                            dist2: d2,
                            point: q,
                            triangle: t,
                            feature: f,
                        };
                    }
                }
            } else {
                // Visit the nearer child first.
                let (l, r) = (nd.first as usize, nd.right as usize);
                let (dl, dr) = (dist2_box(p, self.nodes[l].lo, self.nodes[l].hi), dist2_box(p, self.nodes[r].lo, self.nodes[r].hi));
                if dl < dr {
                    stack.push(r);
                    stack.push(l);
                } else {
                    stack.push(l);
                    stack.push(r);
                }
            }
        }
        Some(best)
    }
}
