//! Fused-deposition printing of a CAD part: what it weighs as printed, the
//! least time it can take, and whether its walls are printable.
//!
//! All three come from one depth field. The part is voxelised on its own
//! field's *sign*, which the kernel gets exactly right everywhere, and an
//! exact Euclidean distance transform (Felzenszwalb & Huttenlocher, *Theory
//! of Computing* 8, 2012) gives each voxel its distance to the nearest
//! outside voxel, in linear time. The depth below the surface is that minus
//! half a cell, good to within a cell. The kernel's own field could not be
//! used for depth: for a union it is only a lower bound, and would erode too
//! much and report walls too thin.
//!
//! * **Printed volume.** Walls of `perimeters × line width` all round
//!   (perimeters, top and bottom skins taken as one thickness), infill
//!   inside. Voxels straddling the wall's inner surface count by the
//!   fraction beyond it. The error is the change from a grid twice as
//!   coarse.
//! * **Least print time** = printed volume / maximum volumetric flow. No
//!   printer extrudes faster than its hotend melts, so this is a true lower
//!   bound; travel, acceleration and cooling only add to it.
//! * **Thinnest wall.** Where the depth field has a crest across some axis
//!   (a run of at most two equal voxels with shallower ones on both
//!   sides), the wall is twice the depth, plus a cell when the crest falls
//!   between two voxels. Convex edges, where the medial axis runs into a
//!   corner, give runs along the face and are not crests. Exact for walls
//!   on grid lines. Off the grid, a voxel's depth is within a half-diagonal
//!   plus a half-cell of the truth, (√3/2 + 1/2) ≈ 1.37 cells, and the
//!   wall doubles it: the claim is 3 cells. Curved walls read thin (a
//!   1 mm tube wall reads 0.69 mm at 0.2 mm cells), which errs toward
//!   flagging, not passing.

use cad::Solid;
use rayon::prelude::*;

#[derive(Clone, Copy, Debug)]
pub struct Print {
    pub printed_volume: f64,
    pub printed_error: f64,
    pub min_wall: f64,
    pub min_wall_error: f64,
}

/// One-dimensional squared distance transform of the sampled function `f`
/// (Felzenszwalb & Huttenlocher's lower envelope of parabolas).
fn dt1(f: &[f64], out: &mut [f64]) {
    let n = f.len();
    let mut v = vec![0usize; n];
    let mut z = vec![0.0f64; n + 1];
    let mut k = 0usize;
    // Skip leading infinities: no parabola is rooted there.
    let Some(first) = f.iter().position(|x| x.is_finite()) else {
        out.iter_mut().for_each(|o| *o = f64::INFINITY);
        return;
    };
    v[0] = first;
    z[0] = f64::NEG_INFINITY;
    z[1] = f64::INFINITY;
    for q in first + 1..n {
        if !f[q].is_finite() {
            continue;
        }
        loop {
            let p = v[k];
            let s = ((f[q] + (q * q) as f64) - (f[p] + (p * p) as f64)) / (2.0 * (q as f64 - p as f64));
            if s <= z[k] && k > 0 {
                k -= 1;
            } else {
                k += 1;
                v[k] = q;
                z[k] = s;
                z[k + 1] = f64::INFINITY;
                break;
            }
        }
    }
    let mut k = 0;
    for (q, o) in out.iter_mut().enumerate() {
        while z[k + 1] < q as f64 {
            k += 1;
        }
        let d = q as f64 - v[k] as f64;
        *o = d * d + f[v[k]];
    }
}

/// Depth below the surface at every voxel centre (negative or zero
/// outside), on a grid of `cell` covering the solid with one cell of
/// padding. Returns (depths, dims).
fn depth_field(solid: &Solid, cell: f64) -> (Vec<f64>, [usize; 3]) {
    let (lo, hi) = solid.bbox();
    let dims: [usize; 3] = std::array::from_fn(|a| ((hi[a] - lo[a]) / cell).ceil() as usize + 2);
    let [nx, ny, nz] = dims;
    let origin = [lo[0] - cell, lo[1] - cell, lo[2] - cell];
    let idx = |i: usize, j: usize, k: usize| i + nx * (j + ny * k);
    // Squared distance (cells²) to the nearest outside voxel centre.
    let mut d: Vec<f64> = (0..nx * ny * nz)
        .into_par_iter()
        .map(|n| {
            let (i, j, k) = (n % nx, n / nx % ny, n / (nx * ny));
            let p = [origin[0] + cell * (i as f64 + 0.5), origin[1] + cell * (j as f64 + 0.5), origin[2] + cell * (k as f64 + 0.5)];
            if solid.sdf(p) < 0.0 { f64::INFINITY } else { 0.0 }
        })
        .collect();
    // x lines are contiguous.
    d.par_chunks_mut(nx).for_each(|line| {
        let f = line.to_vec();
        dt1(&f, line);
    });
    // y and z lines, gathered and scattered.
    for axis in [1usize, 2] {
        let (len, lines): (usize, Vec<(usize, usize)>) = if axis == 1 {
            (ny, (0..nz).flat_map(|k| (0..nx).map(move |i| (i, k))).collect())
        } else {
            (nz, (0..ny).flat_map(|j| (0..nx).map(move |i| (i, j))).collect())
        };
        let results: Vec<Vec<f64>> = lines
            .par_iter()
            .map(|&(a, b)| {
                let f: Vec<f64> = (0..len).map(|t| if axis == 1 { d[idx(a, t, b)] } else { d[idx(a, b, t)] }).collect();
                let mut o = vec![0.0; len];
                dt1(&f, &mut o);
                o
            })
            .collect();
        for (&(a, b), o) in lines.iter().zip(results) {
            for (t, v) in o.into_iter().enumerate() {
                let n = if axis == 1 { idx(a, t, b) } else { idx(a, b, t) };
                d[n] = v;
            }
        }
    }
    let depth = d.into_par_iter().map(|s| if s > 0.0 { s.sqrt() * cell - 0.5 * cell } else { 0.0 }).collect();
    (depth, dims)
}

/// Volume deeper than `wall`: each voxel contributes the fraction of its
/// thickness beyond the wall (exact for a face aligned with the grid).
fn inner_volume(depth: &[f64], cell: f64, wall: f64) -> f64 {
    let v = cell * cell * cell;
    depth.iter().filter(|d| **d > 0.0).map(|&d| v * ((d - wall) / cell + 0.5).clamp(0.0, 1.0)).sum()
}

/// `volume` is the part's own volume (from its mesh); the walls and infill
/// split it.
pub fn analyse(solid: &Solid, volume: f64, cell: f64, wall: f64, infill: f64) -> Result<Print, String> {
    if !(0.0..=1.0).contains(&infill) {
        return Err(format!("infill must be between 0 and 100 %, got {:.1} %", infill * 100.0));
    }
    let (fine, [nx, ny, nz]) = depth_field(solid, cell);
    let (coarse, _) = depth_field(solid, 2.0 * cell);
    let (inf, inc) = (inner_volume(&fine, cell, wall), inner_volume(&coarse, 2.0 * cell, wall));
    let pf = volume - (1.0 - infill) * inf;
    let idx = |i: usize, j: usize, k: usize| i + nx * (j + ny * k);
    let step = |i: usize, j: usize, k: usize, a: usize, s: isize| -> (usize, usize, usize) {
        let mut c = [i as isize, j as isize, k as isize];
        c[a] += s;
        (c[0] as usize, c[1] as usize, c[2] as usize)
    };
    // A crest across axis `a`: a run of at most two equal voxels with
    // strictly shallower voxels on both sides. Convex edges of a body give
    // runs that continue along the face, so they are not crests.
    let crest = |i: usize, j: usize, k: usize, a: usize| -> Option<f64> {
        let d = fine[idx(i, j, k)];
        let n = [nx, ny, nz][a] as isize;
        let pos = [i, j, k][a] as isize;
        let mut run = 1;
        let mut bounded = [false; 2];
        for (side, s) in [(0, -1isize), (1, 1)] {
            let mut t = s;
            loop {
                if pos + t < 0 || pos + t >= n {
                    break;
                }
                let (x, y, z) = step(i, j, k, a, t);
                let e = fine[idx(x, y, z)];
                if e < d {
                    bounded[side] = true;
                    break;
                }
                if e > d {
                    break;
                }
                run += 1;
                if run > 2 {
                    return None;
                }
                t += s;
            }
        }
        (bounded[0] && bounded[1] && run <= 2).then(|| 2.0 * d + (run - 1) as f64 * cell)
    };
    let thinnest = (1..nz - 1)
        .into_par_iter()
        .map(|k| {
            let mut t = f64::INFINITY;
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    if fine[idx(i, j, k)] <= 0.0 {
                        continue;
                    }
                    for a in 0..3 {
                        if let Some(w) = crest(i, j, k, a) {
                            t = t.min(w);
                        }
                    }
                }
            }
            t
        })
        .reduce(|| f64::INFINITY, f64::min);
    if !thinnest.is_finite() {
        return Err(format!("no interior found at cell {:.3} mm: the part is thinner than the grid", cell * 1e3));
    }
    Ok(Print {
        printed_volume: pf,
        printed_error: (1.0 - infill) * (inf - inc).abs(),
        min_wall: thinnest,
        min_wall_error: 3.0 * cell,
    })
}
