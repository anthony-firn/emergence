//! Cell-linked neighbour lists for a cubic periodic box.
//!
//! For every atom the builder returns the full list of neighbours within
//! `rc + skin` (each pair appears twice, once per endpoint). The full list
//! is deliberate: the force loop is then embarrassingly parallel over
//! atoms, each atom writing only its own force slot, and Newton's third law
//! is recovered by summing both endpoints.
//!
//! Cells are enumerated in `[0, L)^3`; positions need not be wrapped (the
//! cell index is taken modulo the grid), and the minimum-image convention
//! is applied to every candidate pair. The minimum box length is required
//! to exceed `2 rc` so that a particle never interacts with its own image.

use rayon::prelude::*;

/// Minimum-image displacement `r_i - r_j` for a cubic box of side `l`.
#[inline]
pub fn min_image(mut d: [f64; 3], l: f64) -> [f64; 3] {
    d[0] -= l * (d[0] / l).round();
    d[1] -= l * (d[1] / l).round();
    d[2] -= l * (d[2] / l).round();
    d
}

/// Full neighbour lists within `rc + skin`.
///
/// `pos` may lie outside the box; min-image distances are always used.
pub fn build_neighbor_lists(pos: &[[f64; 3]], l: f64, rc: f64, skin: f64) -> Vec<Vec<usize>> {
    let n = pos.len();
    assert!(n > 0);
    let rmax = rc + skin;
    assert!(l > 2.0 * rc, "box {l} too small for min-image with rc={rc}");
    // Aim for a cell width of about rmax; clamp so there is always a grid.
    let ncells: [usize; 3] = std::array::from_fn(|_| ((l / rmax).floor() as usize).max(1));
    let nc = ncells[0] * ncells[1] * ncells[2];

    // Bin atoms. Cell index in [0, ncells).
    let mut cells: Vec<Vec<usize>> = vec![Vec::new(); nc];
    let cell_index = |p: &[f64; 3]| -> usize {
        let ix = ((p[0] / l).floor() as i64).rem_euclid(ncells[0] as i64) as usize;
        let iy = ((p[1] / l).floor() as i64).rem_euclid(ncells[1] as i64) as usize;
        let iz = ((p[2] / l).floor() as i64).rem_euclid(ncells[2] as i64) as usize;
        ix + ncells[0] * (iy + ncells[1] * iz)
    };
    for (i, p) in pos.iter().enumerate() {
        cells[cell_index(p)].push(i);
    }

    let rmax2 = rmax * rmax;
    (0..n)
        .into_par_iter()
        .map_init(
            || vec![0u32; nc],
            |stamp, i| {
                // SplitMix-ish generation counter to avoid clearing.
                let mut out: Vec<usize> = Vec::new();
                let pi = pos[i];
                let home = cell_index(&pi);
                let hx = home % ncells[0];
                let hy = (home / ncells[0]) % ncells[1];
                let hz = home / (ncells[0] * ncells[1]);
                let gen = i as u32 + 1;
                for dz in -1i64..=1 {
                    let cz = (hz as i64 + dz).rem_euclid(ncells[2] as i64) as usize;
                    for dy in -1i64..=1 {
                        let cy = (hy as i64 + dy).rem_euclid(ncells[1] as i64) as usize;
                        for dx in -1i64..=1 {
                            let cx = (hx as i64 + dx).rem_euclid(ncells[0] as i64) as usize;
                            let cid = cx + ncells[0] * (cy + ncells[1] * cz);
                            if stamp[cid] == gen {
                                continue;
                            }
                            stamp[cid] = gen;
                            for &j in &cells[cid] {
                                if j == i {
                                    continue;
                                }
                                let d = min_image(
                                    [pi[0] - pos[j][0], pi[1] - pos[j][1], pi[2] - pos[j][2]],
                                    l,
                                );
                                if d[0] * d[0] + d[1] * d[1] + d[2] * d[2] <= rmax2 {
                                    out.push(j);
                                }
                            }
                        }
                    }
                }
                out
            },
        )
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn min_image_is_bounded_by_half_the_box() {
        let l = 10.0;
        for &x in &[-24.3, -4.9, 0.0, 3.2, 9.99, 25.1] {
            let d = min_image([x, 0.0, 0.0], l);
            assert!(d[0].abs() <= l / 2.0 + 1e-12, "{x} -> {}", d[0]);
        }
    }

    #[test]
    fn neighbour_list_matches_brute_force() {
        // Deterministic pseudo-random configuration.
        let mut x: u64 = 12345;
        let mut next = || {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((x >> 11) as f64) / 9_007_199_254_740_992.0
        };
        let n = 200;
        let l = 8.0;
        let pos: Vec<[f64; 3]> = (0..n)
            .map(|_| [next() * l, next() * l, next() * l])
            .collect();
        let rc = 2.0;
        let skin = 0.3;
        let lists = build_neighbor_lists(&pos, l, rc, skin);
        let rmax2 = (rc + skin) * (rc + skin);
        for i in 0..n {
            let mut expected: Vec<usize> = (0..n)
                .filter(|&j| {
                    if j == i {
                        return false;
                    }
                    let d = min_image(
                        [
                            pos[i][0] - pos[j][0],
                            pos[i][1] - pos[j][1],
                            pos[i][2] - pos[j][2],
                        ],
                        l,
                    );
                    d[0] * d[0] + d[1] * d[1] + d[2] * d[2] <= rmax2
                })
                .collect();
            expected.sort_unstable();
            let mut got = lists[i].clone();
            let before = got.len();
            got.sort_unstable();
            got.dedup();
            assert_eq!(before, got.len(), "duplicate neighbours for atom {i}");
            assert_eq!(got, expected, "neighbour mismatch for atom {i}");
        }
    }
}
