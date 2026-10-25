use crate::math::{Bounds, Ray, Vec3};

#[derive(Clone, Debug)]
pub enum Shape {
    Sphere { center: Vec3, radius: f64 },
    Triangle { p: [Vec3; 3] },
}
#[derive(Clone, Debug)]
pub struct Primitive {
    pub shape: Shape,
    pub material: usize,
}
#[derive(Clone, Copy, Debug)]
pub struct Hit {
    pub distance: f64,
    pub point: Vec3,
    pub normal: Vec3,
    pub material: usize,
}

struct PreparedPrimitive {
    primitive: Primitive,
    bounds: Bounds,
    normal: Vec3,
}
struct Query {
    ray: Ray,
    inverse: Vec3,
    axes: [usize; 3],
    shear: [f64; 3],
}
impl Query {
    fn new(ray: Ray) -> Self {
        let z = (0..3)
            .max_by(|&a, &b| ray.dir.0[a].abs().total_cmp(&ray.dir.0[b].abs()))
            .unwrap();
        let (mut x, mut y) = ((z + 1) % 3, (z + 2) % 3);
        if ray.dir.0[z] < 0.0 {
            std::mem::swap(&mut x, &mut y);
        }
        let inverse = Vec3(ray.dir.0.map(|v| 1.0 / v));
        Self {
            ray,
            inverse,
            axes: [x, y, z],
            shear: [
                -ray.dir.0[x] / ray.dir.0[z],
                -ray.dir.0[y] / ray.dir.0[z],
                inverse.0[z],
            ],
        }
    }
    fn bounds(&self, b: Bounds, limit: f64) -> Option<(f64, f64)> {
        b.interval_cached(self.ray, self.inverse, limit)
    }
}
impl PreparedPrimitive {
    fn new(primitive: Primitive) -> Self {
        let (bounds, normal) = match &primitive.shape {
            Shape::Sphere { center, radius } => {
                assert!(radius.is_finite() && *radius > 0.0);
                (
                    Bounds {
                        lo: *center - Vec3([*radius; 3]),
                        hi: *center + Vec3([*radius; 3]),
                    },
                    Vec3::default(),
                )
            }
            Shape::Triangle { p } => {
                assert!(p.iter().flat_map(|p| p.0).all(|v| v.is_finite()));
                let bounds = Bounds {
                    lo: Vec3(std::array::from_fn(|d| {
                        p.iter()
                            .map(|p| p.0[d])
                            .fold(f64::INFINITY, f64::min)
                            .next_down()
                    })),
                    hi: Vec3(std::array::from_fn(|d| {
                        p.iter()
                            .map(|p| p.0[d])
                            .fold(f64::NEG_INFINITY, f64::max)
                            .next_up()
                    })),
                };
                let cross = (p[1] - p[0]).cross(p[2] - p[0]);
                let normal = if cross.length() > 0.0 {
                    cross.unit()
                } else {
                    Vec3::default()
                };
                (bounds, normal)
            }
        };
        Self {
            primitive,
            bounds,
            normal,
        }
    }
    fn hit(&self, q: &Query, limit: f64) -> Option<Hit> {
        let r = q.ray;
        let (t, normal) = match &self.primitive.shape {
            Shape::Sphere { center, radius } => {
                let oc = r.origin - *center;
                let b = oc.dot(r.dir);
                let disc = b * b - oc.dot(oc) + radius * radius;
                if disc < 0.0 {
                    return None;
                }
                let root = disc.sqrt();
                let a = -b - root;
                let t = if a > 0.0 { a } else { -b + root };
                (t, (r.at(t) - *center) / *radius)
            }
            Shape::Triangle { p } => {
                // Ray-aligned shear: a shared edge uses identical transformed
                // vertices in both triangles. No absolute area/determinant cutoff.
                let v = p.map(|p| {
                    let a = p - r.origin;
                    let z = a.0[q.axes[2]];
                    [
                        a.0[q.axes[0]] + q.shear[0] * z,
                        a.0[q.axes[1]] + q.shear[1] * z,
                        z * q.shear[2],
                    ]
                });
                let e = [
                    v[1][0] * v[2][1] - v[1][1] * v[2][0],
                    v[2][0] * v[0][1] - v[2][1] * v[0][0],
                    v[0][0] * v[1][1] - v[0][1] * v[1][0],
                ];
                if e.iter().any(|e| *e < 0.0) && e.iter().any(|e| *e > 0.0) {
                    return None;
                }
                let det = e.iter().sum::<f64>();
                if det == 0.0 {
                    return None;
                }
                let ts = e[0] * v[0][2] + e[1] * v[1][2] + e[2] * v[2][2];
                if det > 0.0 {
                    if ts <= 0.0 || ts >= limit * det {
                        return None;
                    }
                } else if ts >= 0.0 || ts <= limit * det {
                    return None;
                }
                (ts / det, self.normal)
            }
        };
        if t > 0.0 && t < limit {
            Some(Hit {
                distance: t,
                point: r.at(t),
                normal,
                material: self.primitive.material,
            })
        } else {
            None
        }
    }
}
enum Node {
    Leaf {
        bounds: Bounds,
        start: usize,
        end: usize,
    },
    Branch {
        bounds: Bounds,
        left: usize,
        right: usize,
    },
}
impl Node {
    fn bounds(&self) -> Bounds {
        match self {
            Self::Leaf { bounds, .. } | Self::Branch { bounds, .. } => *bounds,
        }
    }
}
#[derive(Clone, Copy, Default)]
struct Bin {
    count: usize,
    bounds: Option<Bounds>,
}
fn area(b: Bounds) -> f64 {
    let d = b.hi - b.lo;
    2.0 * (d.0[0] * d.0[1] + d.0[1] * d.0[2] + d.0[2] * d.0[0])
}
fn union(a: Option<Bounds>, b: Option<Bounds>) -> Option<Bounds> {
    match (a, b) {
        (Some(a), Some(b)) => Some(a.union(b)),
        (Some(a), None) | (None, Some(a)) => Some(a),
        (None, None) => None,
    }
}
/// Binned surface-area-heuristic BVH with cached geometry and ray reciprocals.
/// Bounding boxes are scale-aware; shadow rays can stop at their first blocker.
pub struct Bvh {
    primitives: Vec<PreparedPrimitive>,
    nodes: Vec<Node>,
}
impl Bvh {
    pub fn new(primitives: Vec<Primitive>) -> Self {
        let mut b = Self {
            primitives: primitives.into_iter().map(PreparedPrimitive::new).collect(),
            nodes: vec![],
        };
        if !b.primitives.is_empty() {
            b.build(0, b.primitives.len(), 0);
        }
        b
    }
    fn build(&mut self, start: usize, end: usize, depth: usize) -> usize {
        const BINS: usize = 16;
        let bounds = self.primitives[start..end]
            .iter()
            .map(|p| p.bounds)
            .reduce(Bounds::union)
            .unwrap();
        let id = self.nodes.len();
        self.nodes.push(Node::Leaf { bounds, start, end });
        if end - start <= 4 || depth >= 56 {
            return id;
        }
        let centers = self.primitives[start..end]
            .iter()
            .map(|p| {
                let c = p.bounds.center();
                Bounds { lo: c, hi: c }
            })
            .reduce(Bounds::union)
            .unwrap();
        let mut best = None;
        let mut best_cost = f64::INFINITY;
        for axis in 0..3 {
            let width = centers.hi.0[axis] - centers.lo.0[axis];
            if width <= 0.0 {
                continue;
            }
            let bin = |p: &PreparedPrimitive| {
                (((p.bounds.center().0[axis] - centers.lo.0[axis]) / width * BINS as f64) as usize)
                    .min(BINS - 1)
            };
            let mut bins = [Bin::default(); BINS];
            for p in &self.primitives[start..end] {
                let b = &mut bins[bin(p)];
                b.count += 1;
                b.bounds = union(b.bounds, Some(p.bounds));
            }
            let mut left = [Bin::default(); BINS];
            let mut right = [Bin::default(); BINS];
            let mut acc = Bin::default();
            for i in 0..BINS {
                acc.count += bins[i].count;
                acc.bounds = union(acc.bounds, bins[i].bounds);
                left[i] = acc;
            }
            acc = Bin::default();
            for i in (0..BINS).rev() {
                acc.count += bins[i].count;
                acc.bounds = union(acc.bounds, bins[i].bounds);
                right[i] = acc;
            }
            for i in 0..BINS - 1 {
                if let (Some(l), Some(r)) = (left[i].bounds, right[i + 1].bounds) {
                    let cost = left[i].count as f64 * area(l) + right[i + 1].count as f64 * area(r);
                    if cost < best_cost {
                        best_cost = cost;
                        best = Some((axis, i, width));
                    }
                }
            }
        }
        let mut mid = start;
        if let Some((axis, split, width)) = best {
            // Traversal costs one node; a small leaf is cheaper than a poor split.
            if end - start <= 8 && best_cost + area(bounds) >= (end - start) as f64 * area(bounds) {
                return id;
            }
            for i in start..end {
                let b = (((self.primitives[i].bounds.center().0[axis] - centers.lo.0[axis]) / width
                    * BINS as f64) as usize)
                    .min(BINS - 1);
                if b <= split {
                    self.primitives.swap(i, mid);
                    mid += 1;
                }
            }
        }
        if mid == start || mid == end {
            let extent = centers.hi - centers.lo;
            let axis = (0..3)
                .max_by(|&a, &b| extent.0[a].total_cmp(&extent.0[b]))
                .unwrap();
            mid = (start + end) / 2;
            self.primitives[start..end].select_nth_unstable_by(mid - start, |a, b| {
                a.bounds.center().0[axis].total_cmp(&b.bounds.center().0[axis])
            });
        }
        let left = self.build(start, mid, depth + 1);
        let right = self.build(mid, end, depth + 1);
        self.nodes[id] = Node::Branch {
            bounds,
            left,
            right,
        };
        id
    }
    pub fn len(&self) -> usize {
        self.primitives.len()
    }
    pub fn is_empty(&self) -> bool {
        self.primitives.is_empty()
    }
    pub fn hit(&self, r: Ray, limit: f64) -> Option<Hit> {
        self.trace(r, limit, false)
    }
    pub fn occluded(&self, r: Ray, limit: f64) -> bool {
        self.trace(r, limit, true).is_some()
    }
    fn trace(&self, r: Ray, limit: f64, any: bool) -> Option<Hit> {
        if self.nodes.is_empty() {
            return None;
        }
        let q = Query::new(r);
        let root = q.bounds(self.nodes[0].bounds(), limit)?;
        let mut stack = [(0usize, 0.0); 64];
        stack[0] = (0, root.0);
        let mut count = 1;
        let mut best = None;
        let mut distance = limit;
        while count > 0 {
            count -= 1;
            let (id, entry) = stack[count];
            if entry > distance {
                continue;
            }
            match &self.nodes[id] {
                Node::Leaf { start, end, .. } => {
                    for p in &self.primitives[*start..*end] {
                        if let Some(h) = p.hit(&q, distance) {
                            if any {
                                return Some(h);
                            }
                            distance = h.distance;
                            best = Some(h);
                        }
                    }
                }
                Node::Branch { left, right, .. } => {
                    let a = q.bounds(self.nodes[*left].bounds(), distance);
                    let b = q.bounds(self.nodes[*right].bounds(), distance);
                    match (a, b) {
                        (Some(a), Some(b)) => {
                            assert!(count + 2 <= stack.len());
                            let (near, far) = if a.0 < b.0 {
                                ((*left, a.0), (*right, b.0))
                            } else {
                                ((*right, b.0), (*left, a.0))
                            };
                            stack[count] = far;
                            stack[count + 1] = near;
                            count += 2;
                        }
                        (Some(a), None) => {
                            stack[count] = (*left, a.0);
                            count += 1;
                        }
                        (None, Some(b)) => {
                            stack[count] = (*right, b.0);
                            count += 1;
                        }
                        (None, None) => {}
                    }
                }
            }
        }
        best
    }
}
pub fn quad(p: Vec3, u: Vec3, v: Vec3, material: usize) -> Vec<Primitive> {
    vec![
        Primitive {
            shape: Shape::Triangle {
                p: [p, p + u, p + u + v],
            },
            material,
        },
        Primitive {
            shape: Shape::Triangle {
                p: [p, p + u + v, p + v],
            },
            material,
        },
    ]
}
/// Marching tetrahedra on cell-centred volume fractions. Surface normals
/// point out of the high-fraction phase. No invented fine-scale geometry.
/// The phase must vanish on the outer cell centres (otherwise it is open).
pub fn isosurface(
    n: [usize; 3],
    lo: Vec3,
    dx: Vec3,
    alpha: &[f64],
    iso: f64,
    material: usize,
) -> Vec<Primitive> {
    assert_eq!(n.iter().product::<usize>(), alpha.len());
    assert!(n.iter().all(|&v| v >= 2));
    assert!(alpha.iter().all(|a| a.is_finite()));
    let index = |i: [usize; 3]| (i[0] * n[1] + i[1]) * n[2] + i[2];
    let corner = [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
    ];
    // Consistent face diagonals across adjacent cubes.
    let tets = [
        [0, 1, 2, 6],
        [0, 2, 3, 6],
        [0, 3, 7, 6],
        [0, 7, 4, 6],
        [0, 4, 5, 6],
        [0, 5, 1, 6],
    ];
    let mut result = vec![];
    for i in 0..n[0] - 1 {
        for j in 0..n[1] - 1 {
            for k in 0..n[2] - 1 {
                let ids = corner.map(|c| [i + c[0], j + c[1], k + c[2]]);
                let p =
                    ids.map(|c| lo + Vec3(std::array::from_fn(|d| (c[d] as f64 + 0.5) * dx.0[d])));
                let a = ids.map(|c| alpha[index(c)]);
                for tet in tets {
                    let inside: Vec<usize> = tet.iter().copied().filter(|&q| a[q] >= iso).collect();
                    let outside: Vec<usize> = tet.iter().copied().filter(|&q| a[q] < iso).collect();
                    if inside.is_empty() || outside.is_empty() {
                        continue;
                    }
                    let outward = outside
                        .iter()
                        .map(|&q| p[q])
                        .fold(Vec3::default(), |x, y| x + y)
                        / outside.len() as f64
                        - inside
                            .iter()
                            .map(|&q| p[q])
                            .fold(Vec3::default(), |x, y| x + y)
                            / inside.len() as f64;
                    let edge =
                        |q: usize, r: usize| p[q] + (p[r] - p[q]) * ((iso - a[q]) / (a[r] - a[q]));
                    let mut emit = |mut points: [Vec3; 3]| {
                        let normal = (points[1] - points[0]).cross(points[2] - points[0]);
                        if normal.length() == 0.0 {
                            return;
                        }
                        if normal.dot(outward) < 0.0 {
                            points.swap(1, 2);
                        }
                        result.push(Primitive {
                            shape: Shape::Triangle { p: points },
                            material,
                        });
                    };
                    if inside.len() == 1 {
                        let q = inside[0];
                        emit([
                            edge(q, outside[0]),
                            edge(q, outside[1]),
                            edge(q, outside[2]),
                        ]);
                    } else if outside.len() == 1 {
                        let q = outside[0];
                        emit([edge(q, inside[0]), edge(q, inside[1]), edge(q, inside[2])]);
                    } else {
                        let x = [
                            edge(inside[0], outside[0]),
                            edge(inside[0], outside[1]),
                            edge(inside[1], outside[1]),
                            edge(inside[1], outside[0]),
                        ];
                        emit([x[0], x[1], x[2]]);
                        emit([x[0], x[2], x[3]]);
                    }
                }
            }
        }
    }
    result
}
