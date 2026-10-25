//! Closed, outward-oriented test shapes with known volume and area.

use std::collections::HashMap;
use std::f64::consts::PI;

use crate::mesh::TriMesh;
use crate::{add, normalize, scale, P3};

/// Icosphere of radius `r`: an icosahedron subdivided `levels` times with
/// vertices projected to the sphere (20·4^levels triangles).
pub fn icosphere(center: P3, r: f64, levels: u32) -> TriMesh {
    let t = (1.0 + 5f64.sqrt()) / 2.0;
    let mut v: Vec<P3> = [
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ]
    .iter()
    .map(|&p| normalize(p))
    .collect();
    let mut f: Vec<[u32; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    for _ in 0..levels {
        let mut mid: HashMap<(u32, u32), u32> = HashMap::new();
        let mut nf = Vec::with_capacity(4 * f.len());
        for tri in &f {
            let mut m = [0u32; 3];
            for k in 0..3 {
                let (a, b) = (tri[k], tri[(k + 1) % 3]);
                m[k] = *mid.entry((a.min(b), a.max(b))).or_insert_with(|| {
                    v.push(normalize(add(v[a as usize], v[b as usize])));
                    v.len() as u32 - 1
                });
            }
            nf.push([tri[0], m[0], m[2]]);
            nf.push([tri[1], m[1], m[0]]);
            nf.push([tri[2], m[2], m[1]]);
            nf.push(m);
        }
        f = nf;
    }
    let mut mesh = TriMesh::new(v.iter().map(|&p| add(center, scale(p, r))).collect(), f);
    mesh.orient_outward();
    mesh
}

/// Closed cylinder along +z: base centre `c0`, radius `r`, length `len`,
/// `segments` around (with fan caps).
pub fn cylinder_z(c0: P3, r: f64, len: f64, segments: usize) -> TriMesh {
    let n = segments as u32;
    let mut v = Vec::with_capacity(2 * segments + 2);
    for z in [0.0, len] {
        for k in 0..segments {
            let a = 2.0 * PI * k as f64 / segments as f64;
            v.push([c0[0] + r * a.cos(), c0[1] + r * a.sin(), c0[2] + z]);
        }
    }
    v.push(c0);
    v.push([c0[0], c0[1], c0[2] + len]);
    let (cb, ct) = (2 * n, 2 * n + 1);
    let mut f = Vec::new();
    for k in 0..n {
        let k1 = (k + 1) % n;
        f.push([k, k1, n + k1]);
        f.push([k, n + k1, n + k]);
        f.push([cb, k1, k]);
        f.push([ct, n + k, n + k1]);
    }
    let mut mesh = TriMesh::new(v, f);
    mesh.orient_outward();
    mesh
}

/// Axis-aligned box `[lo, hi]` (12 triangles).
pub fn cuboid(lo: P3, hi: P3) -> TriMesh {
    let c = |i: usize| [if i & 1 == 0 { lo[0] } else { hi[0] }, if i & 2 == 0 { lo[1] } else { hi[1] }, if i & 4 == 0 { lo[2] } else { hi[2] }];
    let v = (0..8).map(c).collect();
    let quads = [[0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [2, 6, 7, 3], [0, 4, 6, 2], [1, 3, 7, 5]];
    let mut f = Vec::new();
    for q in quads {
        f.push([q[0], q[1], q[2]]);
        f.push([q[0], q[2], q[3]]);
    }
    let mut mesh = TriMesh::new(v, f);
    mesh.orient_outward();
    mesh
}

/// Torus around the z axis: major radius `big`, minor `small`.
pub fn torus(center: P3, big: f64, small: f64, nu: usize, nv: usize) -> TriMesh {
    let mut v = Vec::with_capacity(nu * nv);
    for i in 0..nu {
        let u = 2.0 * PI * i as f64 / nu as f64;
        for j in 0..nv {
            let w = 2.0 * PI * j as f64 / nv as f64;
            let rr = big + small * w.cos();
            v.push([center[0] + rr * u.cos(), center[1] + rr * u.sin(), center[2] + small * w.sin()]);
        }
    }
    let id = |i: usize, j: usize| ((i % nu) * nv + (j % nv)) as u32;
    let mut f = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            f.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
            f.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
        }
    }
    let mut mesh = TriMesh::new(v, f);
    mesh.orient_outward();
    mesh
}
