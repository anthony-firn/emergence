//! From the implicit kernel to triangles.
//!
//! [`surface_nets`] is the dual isosurface method (Gibson, *MICCAI* 1998):
//! one vertex per grid cell that contains a sign change, and one quad per
//! grid edge that contains one. Its two properties are what the solvers
//! need:
//!
//! * the output is **closed by construction** — every quad edge is shared
//!   by exactly two quads, because each is contributed by one of the four
//!   cells around a sign-changing edge, so there are never boundary edges.
//!   [`TriMesh::report`] checks this, and the tests assert it.
//!
//!   It is *manifold* as well whenever the surface passes through each
//!   cell once. Where two sheets of the surface meet inside one cell — two
//!   solids in contact, or a gap thinner than the cell — the single vertex
//!   of that cell is shared by both sheets and the edge through it is
//!   non-manifold. That is the correct answer for tangent solids, whose
//!   boundary really is non-manifold there, and it is a resolution
//!   artefact for a thin gap, which refining removes.
//!   `geometry::Sdf` falls back from pseudonormals to winding numbers for
//!   such a mesh, which is correct but `O(triangles)` per query;
//! * the **topology follows the sign of the field**, which the kernel
//!   invariant makes exact, so the mesh cannot have a hole that the model
//!   does not have.
//!
//! The vertex of a cell starts at the mean of the zero crossings on its
//! edges and is then projected onto the zero set by a few damped Newton
//! steps along the field gradient, clamped inside the cell. The projection
//! is what makes the volume and area converge at second order with a small
//! constant; `analysis/cad.md` has the measured rates.
//!
//! Marching cubes (Lorensen & Cline, *SIGGRAPH* 1987) would place vertices
//! on grid edges instead. It needs a 256-case table and, in its original
//! form, is not watertight for the ambiguous cases; the dual method avoids
//! both problems and loses only sharp-feature fidelity, which an isotropic
//! grid cannot represent anyway.

use rayon::prelude::*;

use geometry::{add, dot, norm, scale, sub, TriMesh, P3};
use std::collections::HashMap;

use crate::solid::Solid;

/// A regular node grid: node `(i, j, k)` sits at
/// `origin + (i, j, k) * spacing`, with the spacing allowed to differ per
/// axis.
///
/// Anisotropic cells are what make a thin feature in a wide model
/// meshable. A printed circuit board is 35 um of copper over 50 mm of
/// board: an isotropic grid fine enough for the copper needs 1e11 cells,
/// while one with 10 um in `z` and 200 um in `x` and `y` needs 1e7 and
/// resolves the copper just as well, because nothing varies quickly in the
/// plane.
#[derive(Clone, Copy, Debug)]
pub struct Grid3 {
    pub origin: P3,
    pub spacing: P3,
    /// Number of **nodes** along each axis.
    pub dims: [usize; 3],
}

impl Grid3 {
    /// A grid covering `bounds` with isotropic step `cell`, padded by
    /// `pad` cells on every side so that the outer node layer lies outside
    /// any solid contained in `bounds`.
    pub fn covering(bounds: (P3, P3), cell: f64, pad: usize) -> Result<Self, String> {
        Self::covering_with(bounds, [cell; 3], pad)
    }

    /// As [`Grid3::covering`], with a different step on each axis.
    pub fn covering_with(bounds: (P3, P3), spacing: P3, pad: usize) -> Result<Self, String> {
        if (0..3).any(|k| !(spacing[k] > 0.0)) {
            return Err("every cell size must be positive".into());
        }
        let mut origin = [0.0; 3];
        let mut dims = [0usize; 3];
        for k in 0..3 {
            let (lo, hi) = (bounds.0[k], bounds.1[k]);
            if !lo.is_finite() || !hi.is_finite() {
                return Err(format!(
                    "unbounded solid along axis {k}: clip it with Solid::clipped before meshing"
                ));
            }
            let n = ((hi - lo) / spacing[k]).ceil().max(1.0) as usize + 1 + 2 * pad;
            origin[k] = lo - pad as f64 * spacing[k];
            dims[k] = n;
        }
        Ok(Grid3 { origin, spacing, dims })
    }

    pub fn nodes(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    /// Number of cells along each axis.
    pub fn cells(&self) -> [usize; 3] {
        std::array::from_fn(|k| self.dims[k] - 1)
    }

    /// The smallest cell dimension, which is what sets the resolvable
    /// feature size.
    pub fn finest(&self) -> f64 {
        self.spacing.iter().copied().fold(f64::INFINITY, f64::min)
    }

    /// Half the cell diagonal: the distance from a cell centre to its
    /// farthest corner.
    pub fn half_diagonal(&self) -> f64 {
        0.5 * norm(self.spacing)
    }

    #[inline]
    pub fn point(&self, i: usize, j: usize, k: usize) -> P3 {
        [
            self.origin[0] + self.spacing[0] * i as f64,
            self.origin[1] + self.spacing[1] * j as f64,
            self.origin[2] + self.spacing[2] * k as f64,
        ]
    }

    #[inline]
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.dims[0] * (j + self.dims[1] * k)
    }

    /// The box spanned by the nodes.
    pub fn bounds(&self) -> (P3, P3) {
        (
            self.origin,
            std::array::from_fn(|k| self.origin[k] + self.spacing[k] * (self.dims[k] - 1) as f64),
        )
    }
}

/// Sample the kernel field on every node, in parallel.
pub fn field_on(solid: &Solid, grid: &Grid3) -> Vec<f64> {
    let [nx, ny, _] = grid.dims;
    (0..grid.nodes())
        .into_par_iter()
        .map(|q| {
            let (i, j, k) = (q % nx, (q / nx) % ny, q / (nx * ny));
            solid.sdf(grid.point(i, j, k))
        })
        .collect()
}

/// The 12 edges of a cell as pairs of corner indices (corner `c` has
/// offsets `(c & 1, (c >> 1) & 1, (c >> 2) & 1)`).
const CELL_EDGES: [(usize, usize); 12] = [
    (0, 1), (2, 3), (4, 5), (6, 7), // along x
    (0, 2), (1, 3), (4, 6), (5, 7), // along y
    (0, 4), (1, 5), (2, 6), (3, 7), // along z
];

/// Extract the zero isosurface of `solid` on `grid`.
///
/// `refine` damped Newton steps are applied to each vertex (2 is a good
/// default; 0 gives plain surface nets).
pub fn surface_nets(solid: &Solid, grid: &Grid3, refine: usize) -> TriMesh {
    let f = field_on(solid, grid);
    let [nx, ny, nz] = grid.dims;
    let (cx, cy, cz) = (nx - 1, ny - 1, nz - 1);
    let h = grid.spacing;
    let hmin = grid.finest();

    // One vertex per sign-changing cell.
    let mut cell_vertex = vec![u32::MAX; cx * cy * cz];
    let mut vertices: Vec<P3> = Vec::new();
    let cell_id = |i: usize, j: usize, k: usize| i + cx * (j + cy * k);
    for k in 0..cz {
        for j in 0..cy {
            for i in 0..cx {
                let corner = |c: usize| grid.index(i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1));
                let v: [f64; 8] = std::array::from_fn(|c| f[corner(c)]);
                if v.iter().all(|&x| x < 0.0) || v.iter().all(|&x| x >= 0.0) {
                    continue;
                }
                let base = grid.point(i, j, k);
                let mut sum = [0.0; 3];
                let mut n = 0.0;
                for (a, b) in CELL_EDGES {
                    let (fa, fb) = (v[a], v[b]);
                    if (fa < 0.0) == (fb < 0.0) {
                        continue;
                    }
                    let t = if (fb - fa).abs() > 0.0 { (-fa / (fb - fa)).clamp(0.0, 1.0) } else { 0.5 };
                    let ca = [(a & 1) as f64, ((a >> 1) & 1) as f64, ((a >> 2) & 1) as f64];
                    let cb = [(b & 1) as f64, ((b >> 1) & 1) as f64, ((b >> 2) & 1) as f64];
                    for d in 0..3 {
                        sum[d] += base[d] + h[d] * (ca[d] + t * (cb[d] - ca[d]));
                    }
                    n += 1.0;
                }
                let mut p = scale(sum, 1.0 / n);
                for _ in 0..refine {
                    let val = solid.sdf(p);
                    let g = solid.gradient(p, 0.25 * hmin);
                    let gg = dot(g, g);
                    if gg < 1e-12 {
                        break;
                    }
                    let step = scale(g, -val / gg);
                    // Stay inside the cell: the dual mesh stays manifold.
                    if norm(step) > norm(h) {
                        break;
                    }
                    let q = add(p, step);
                    let ok = (0..3).all(|d| q[d] >= base[d] - 1e-12 && q[d] <= base[d] + h[d] + 1e-12);
                    if !ok {
                        break;
                    }
                    p = q;
                }
                cell_vertex[cell_id(i, j, k)] = vertices.len() as u32;
                vertices.push(p);
            }
        }
    }

    // One quad per sign-changing node edge, from the four cells around it.
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    let quad = |a: u32, b: u32, c: u32, d: u32, flip: bool, tris: &mut Vec<[u32; 3]>| {
        if a == u32::MAX || b == u32::MAX || c == u32::MAX || d == u32::MAX {
            return;
        }
        if flip {
            tris.push([a, c, b]);
            tris.push([a, d, c]);
        } else {
            tris.push([a, b, c]);
            tris.push([a, c, d]);
        }
    };
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let here = f[grid.index(i, j, k)];
                // +x edge: cells share (y-1,y) x (z-1,z); cyclic y -> z.
                if i + 1 < nx && j >= 1 && k >= 1 && j < cy && k < cz {
                    let there = f[grid.index(i + 1, j, k)];
                    if (here < 0.0) != (there < 0.0) {
                        quad(
                            cell_vertex[cell_id(i, j - 1, k - 1)],
                            cell_vertex[cell_id(i, j, k - 1)],
                            cell_vertex[cell_id(i, j, k)],
                            cell_vertex[cell_id(i, j - 1, k)],
                            here >= 0.0,
                            &mut triangles,
                        );
                    }
                }
                // +y edge: cyclic z -> x.
                if j + 1 < ny && k >= 1 && i >= 1 && k < cz && i < cx {
                    let there = f[grid.index(i, j + 1, k)];
                    if (here < 0.0) != (there < 0.0) {
                        quad(
                            cell_vertex[cell_id(i - 1, j, k - 1)],
                            cell_vertex[cell_id(i - 1, j, k)],
                            cell_vertex[cell_id(i, j, k)],
                            cell_vertex[cell_id(i, j, k - 1)],
                            here >= 0.0,
                            &mut triangles,
                        );
                    }
                }
                // +z edge: cyclic x -> y.
                if k + 1 < nz && i >= 1 && j >= 1 && i < cx && j < cy {
                    let there = f[grid.index(i, j, k + 1)];
                    if (here < 0.0) != (there < 0.0) {
                        quad(
                            cell_vertex[cell_id(i - 1, j - 1, k)],
                            cell_vertex[cell_id(i, j - 1, k)],
                            cell_vertex[cell_id(i, j, k)],
                            cell_vertex[cell_id(i - 1, j, k)],
                            here >= 0.0,
                            &mut triangles,
                        );
                    }
                }
            }
        }
    }

    // Drop vertices no quad used (cells whose sign change touched only the
    // padded boundary).
    let mut used = vec![u32::MAX; vertices.len()];
    let mut kept: Vec<P3> = Vec::new();
    for t in &mut triangles {
        for v in t.iter_mut() {
            if used[*v as usize] == u32::MAX {
                used[*v as usize] = kept.len() as u32;
                kept.push(vertices[*v as usize]);
            }
            *v = used[*v as usize];
        }
    }
    TriMesh::new(kept, triangles)
}

/// What block culling skipped.
#[derive(Clone, Copy, Debug)]
pub struct SparseStats {
    pub blocks: usize,
    pub blocks_kept: usize,
    /// Field evaluations actually made, against the node count a dense
    /// pass would have needed.
    pub evaluations: usize,
    pub dense_evaluations: usize,
}

impl SparseStats {
    /// Evaluations made as a fraction of the dense count.
    pub fn work_fraction(&self) -> f64 {
        self.evaluations as f64 / self.dense_evaluations.max(1) as f64
    }
}

/// [`surface_nets`], skipping the parts of the grid that provably hold no
/// surface.
///
/// This is where the kernel invariant pays for itself. The field is
/// 1-Lipschitz, so if `|f|` at a block's centre exceeds the distance from
/// that centre to the block's farthest corner, `f` cannot change sign
/// anywhere in the block, and the block holds no surface. One evaluation
/// then rules out `block^3` cells.
///
/// The result is **identical** to the dense pass, not an approximation:
/// every sign-changing cell lies in a kept block (its block cannot be
/// culled), and every quad's four cells all contain the same
/// sign-changing edge, so they are all kept too. Cells are visited in the
/// same order as the dense pass, so even the vertex numbering matches.
/// The test asserts bit-identity.
///
/// Cost falls from the volume of the grid toward its surface area: for a
/// sphere at `r/h = 64` the pass makes a third of the evaluations the
/// dense one would.
///
/// It buys nothing for a **slab**, whose surfaces pass through every block
/// in the plane — a printed circuit board culls no blocks at all. That
/// case is won by anisotropic cells instead ([`Grid3::covering_with`]).
/// The two address different shapes and compose.
pub fn surface_nets_sparse(
    solid: &Solid,
    grid: &Grid3,
    refine: usize,
    block: usize,
) -> (TriMesh, SparseStats) {
    let block = block.max(1);
    let [nx, ny, nz] = grid.dims;
    let (cx, cy, cz) = (nx - 1, ny - 1, nz - 1);
    let h = grid.spacing;
    let hmin = grid.finest();
    let nb: [usize; 3] = std::array::from_fn(|k| grid.cells()[k].div_ceil(block));
    let total_blocks = nb[0] * nb[1] * nb[2];

    // One evaluation per block decides whether it can hold a surface.
    let keep: Vec<bool> = (0..total_blocks)
        .into_par_iter()
        .map(|q| {
            let (bi, bj, bk) = (q % nb[0], (q / nb[0]) % nb[1], q / (nb[0] * nb[1]));
            let lo = [bi * block, bj * block, bk * block];
            let hi: [usize; 3] = std::array::from_fn(|k| ((lo[k] + block).min(grid.cells()[k])).max(lo[k] + 1));
            let a = grid.point(lo[0], lo[1], lo[2]);
            let b = grid.point(hi[0], hi[1], hi[2]);
            let centre: P3 = std::array::from_fn(|k| 0.5 * (a[k] + b[k]));
            let reach = 0.5 * norm(sub(b, a));
            // A 1-Lipschitz field cannot reach zero inside the block if it
            // is farther than `reach` from zero at the centre. The
            // tolerance is one part in 1e12 of the reach, for round-off.
            solid.sdf(centre).abs() <= reach * (1.0 + 1e-12)
        })
        .collect();
    let kept = keep.iter().filter(|&&k| k).count();

    // Evaluate only the nodes of kept blocks, plus the one-node overlap
    // their cells need.
    let mut values: HashMap<usize, f64> = HashMap::with_capacity(kept * (block + 1).pow(3));
    let mut wanted: Vec<usize> = Vec::new();
    for q in 0..total_blocks {
        if !keep[q] {
            continue;
        }
        let (bi, bj, bk) = (q % nb[0], (q / nb[0]) % nb[1], q / (nb[0] * nb[1]));
        let lo = [bi * block, bj * block, bk * block];
        for i in lo[0]..=(lo[0] + block).min(cx) {
            for j in lo[1]..=(lo[1] + block).min(cy) {
                for k in lo[2]..=(lo[2] + block).min(cz) {
                    let idx = grid.index(i, j, k);
                    if !values.contains_key(&idx) {
                        values.insert(idx, f64::NAN);
                        wanted.push(idx);
                    }
                }
            }
        }
    }
    let evaluated: Vec<(usize, f64)> = wanted
        .par_iter()
        .map(|&idx| {
            let (i, j, k) = (idx % nx, (idx / nx) % ny, idx / (nx * ny));
            (idx, solid.sdf(grid.point(i, j, k)))
        })
        .collect();
    let evaluations = evaluated.len() + total_blocks;
    for (idx, v) in evaluated {
        values.insert(idx, v);
    }

    // Is cell (i, j, k) in a kept block?
    let cell_kept = |i: usize, j: usize, k: usize| -> bool {
        keep[(i / block) + nb[0] * ((j / block) + nb[1] * (k / block))]
    };
    let at = |i: usize, j: usize, k: usize| -> Option<f64> {
        values.get(&grid.index(i, j, k)).copied()
    };

    // The same cell and edge passes as the dense version, over kept cells
    // only and in the same order.
    let mut cell_vertex: HashMap<usize, u32> = HashMap::new();
    let mut vertices: Vec<P3> = Vec::new();
    let cell_id = |i: usize, j: usize, k: usize| i + cx * (j + cy * k);
    for k in 0..cz {
        for j in 0..cy {
            for i in 0..cx {
                if !cell_kept(i, j, k) {
                    continue;
                }
                let mut v = [0.0f64; 8];
                let mut have = true;
                for c in 0..8 {
                    match at(i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1)) {
                        Some(x) => v[c] = x,
                        None => have = false,
                    }
                }
                if !have || v.iter().all(|&x| x < 0.0) || v.iter().all(|&x| x >= 0.0) {
                    continue;
                }
                let base = grid.point(i, j, k);
                let mut sum = [0.0; 3];
                let mut n = 0.0;
                for (a, b) in CELL_EDGES {
                    let (fa, fb) = (v[a], v[b]);
                    if (fa < 0.0) == (fb < 0.0) {
                        continue;
                    }
                    let t = if (fb - fa).abs() > 0.0 { (-fa / (fb - fa)).clamp(0.0, 1.0) } else { 0.5 };
                    let ca = [(a & 1) as f64, ((a >> 1) & 1) as f64, ((a >> 2) & 1) as f64];
                    let cb = [(b & 1) as f64, ((b >> 1) & 1) as f64, ((b >> 2) & 1) as f64];
                    for d in 0..3 {
                        sum[d] += base[d] + h[d] * (ca[d] + t * (cb[d] - ca[d]));
                    }
                    n += 1.0;
                }
                let mut p = scale(sum, 1.0 / n);
                for _ in 0..refine {
                    let val = solid.sdf(p);
                    let g = solid.gradient(p, 0.25 * hmin);
                    let gg = dot(g, g);
                    if gg < 1e-12 {
                        break;
                    }
                    let step = scale(g, -val / gg);
                    if norm(step) > norm(h) {
                        break;
                    }
                    let q = add(p, step);
                    let ok = (0..3).all(|d| q[d] >= base[d] - 1e-12 && q[d] <= base[d] + h[d] + 1e-12);
                    if !ok {
                        break;
                    }
                    p = q;
                }
                cell_vertex.insert(cell_id(i, j, k), vertices.len() as u32);
                vertices.push(p);
            }
        }
    }

    let mut triangles: Vec<[u32; 3]> = Vec::new();
    let quad = |a: Option<&u32>, b: Option<&u32>, c: Option<&u32>, d: Option<&u32>, flip: bool, tris: &mut Vec<[u32; 3]>| {
        let (Some(&a), Some(&b), Some(&c), Some(&d)) = (a, b, c, d) else {
            return;
        };
        if flip {
            tris.push([a, c, b]);
            tris.push([a, d, c]);
        } else {
            tris.push([a, b, c]);
            tris.push([a, c, d]);
        }
    };
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let Some(here) = at(i, j, k) else { continue };
                if i + 1 < nx && j >= 1 && k >= 1 && j < cy && k < cz {
                    if let Some(there) = at(i + 1, j, k) {
                        if (here < 0.0) != (there < 0.0) {
                            quad(
                                cell_vertex.get(&cell_id(i, j - 1, k - 1)),
                                cell_vertex.get(&cell_id(i, j, k - 1)),
                                cell_vertex.get(&cell_id(i, j, k)),
                                cell_vertex.get(&cell_id(i, j - 1, k)),
                                here >= 0.0,
                                &mut triangles,
                            );
                        }
                    }
                }
                if j + 1 < ny && k >= 1 && i >= 1 && k < cz && i < cx {
                    if let Some(there) = at(i, j + 1, k) {
                        if (here < 0.0) != (there < 0.0) {
                            quad(
                                cell_vertex.get(&cell_id(i - 1, j, k - 1)),
                                cell_vertex.get(&cell_id(i - 1, j, k)),
                                cell_vertex.get(&cell_id(i, j, k)),
                                cell_vertex.get(&cell_id(i, j, k - 1)),
                                here >= 0.0,
                                &mut triangles,
                            );
                        }
                    }
                }
                if k + 1 < nz && i >= 1 && j >= 1 && i < cx && j < cy {
                    if let Some(there) = at(i, j, k + 1) {
                        if (here < 0.0) != (there < 0.0) {
                            quad(
                                cell_vertex.get(&cell_id(i - 1, j - 1, k)),
                                cell_vertex.get(&cell_id(i, j - 1, k)),
                                cell_vertex.get(&cell_id(i, j, k)),
                                cell_vertex.get(&cell_id(i - 1, j, k)),
                                here >= 0.0,
                                &mut triangles,
                            );
                        }
                    }
                }
            }
        }
    }

    let mut used = vec![u32::MAX; vertices.len()];
    let mut kept_vertices: Vec<P3> = Vec::new();
    for t in &mut triangles {
        for v in t.iter_mut() {
            if used[*v as usize] == u32::MAX {
                used[*v as usize] = kept_vertices.len() as u32;
                kept_vertices.push(vertices[*v as usize]);
            }
            *v = used[*v as usize];
        }
    }
    (
        TriMesh::new(kept_vertices, triangles),
        SparseStats {
            blocks: total_blocks,
            blocks_kept: kept,
            evaluations,
            dense_evaluations: grid.nodes(),
        },
    )
}

/// Mesh a solid at resolution `cell`, choosing the grid from its bounding
/// box. Fails for an unbounded solid.
pub fn mesh_solid(solid: &Solid, cell: f64) -> Result<TriMesh, String> {
    let grid = Grid3::covering(solid.bbox(), cell, 2)?;
    Ok(surface_nets(solid, &grid, 2))
}

/// Mesh a solid over explicit bounds (for clipping a periodic solid, or for
/// matching a solver's grid).
pub fn mesh_solid_in(solid: &Solid, bounds: (P3, P3), cell: f64) -> Result<TriMesh, String> {
    let grid = Grid3::covering(bounds, cell, 2)?;
    Ok(surface_nets(solid, &grid, 2))
}

/// Mesh a solid with a different cell size on each axis.
///
/// This is how a thin feature in a wide model is meshed: give the thin
/// direction a small cell and the others a large one.
pub fn mesh_solid_anisotropic(solid: &Solid, spacing: P3) -> Result<TriMesh, String> {
    let grid = Grid3::covering_with(solid.bbox(), spacing, 2)?;
    Ok(surface_nets(solid, &grid, 2))
}

/// Surface area of the zero set, from the meshed surface.
pub fn surface_area(solid: &Solid, cell: f64) -> Result<f64, String> {
    Ok(mesh_solid(solid, cell)?.area())
}

/// The outward unit normal of the field at `p` (for post-processing).
pub fn normal_at(solid: &Solid, p: P3, h: f64) -> P3 {
    let g = solid.gradient(p, h);
    let n = norm(g);
    if n > 0.0 {
        scale(g, 1.0 / n)
    } else {
        [0.0, 0.0, 1.0]
    }
}

/// Signed distance from a solid onto a solver grid, in the layout
/// `geometry::Sdf::sample_grid` uses (`x` fastest).
pub fn sample_grid(solid: &Solid, origin: P3, h: P3, dims: [usize; 3]) -> Vec<f64> {
    (0..dims[0] * dims[1] * dims[2])
        .into_par_iter()
        .map(|q| {
            let (i, j, k) = (q % dims[0], (q / dims[0]) % dims[1], q / (dims[0] * dims[1]));
            solid.sdf([origin[0] + h[0] * i as f64, origin[1] + h[1] * j as f64, origin[2] + h[2] * k as f64])
        })
        .collect()
}

/// Distance from `p` to the meshed boundary, by brute force. Used by the
/// tests to check the kernel's `|sdf| <= dist` guarantee.
pub fn distance_to_mesh(mesh: &TriMesh, p: P3) -> f64 {
    let mut best = f64::INFINITY;
    for t in 0..mesh.triangles.len() {
        let d = point_triangle_distance(p, mesh.corners(t));
        if d < best {
            best = d;
        }
    }
    best
}

fn point_triangle_distance(p: P3, c: [P3; 3]) -> f64 {
    // Ericson, Real-Time Collision Detection (2005), 5.1.5.
    let (a, b, cc) = (c[0], c[1], c[2]);
    let ab = sub(b, a);
    let ac = sub(cc, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return norm(ap);
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return norm(bp);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return norm(sub(p, add(a, scale(ab, v))));
    }
    let cp = sub(p, cc);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return norm(cp);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return norm(sub(p, add(a, scale(ac, w))));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return norm(sub(p, add(b, scale(sub(cc, b), w))));
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    norm(sub(p, add(a, add(scale(ab, v), scale(ac, w)))))
}
