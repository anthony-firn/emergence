//! Lattice generators: honeycomb, re-entrant honeycomb, triangular, octet
//! truss and simple-cubic beam lattices.
//!
//! Every generator emits overlapping blocks of the lattice; [`crate::cell`]
//! reduces endpoints modulo the lattice vectors and removes duplicate bonds,
//! so the generators can be written naively and stay correct.

use crate::cell::{Kind, RawMember2d, RawMember3d};

/// General hexagonal honeycomb: inclined walls of length `l` at angle
/// `theta` to the horizontal, vertical walls of length `h`, wall thickness
/// `t`. `theta = 30 deg`, `h = l` is the regular honeycomb; `theta < 0` with
/// `h/l > 0.5` is re-entrant (auxetic).
pub fn honeycomb(t: f64, l: f64, h_over_l: f64, theta: f64, _kind: Kind) -> (Vec<RawMember2d>, [f64; 2], [f64; 2], f64) {
    let h = h_over_l * l;
    let lc = l * theta.cos();
    let ls = l * theta.sin();
    let a1 = [2.0 * lc, 0.0];
    let a2 = [lc, h + ls];
    let mut raw = Vec::new();
    for n1 in -2i32..=2 {
        for n2 in -2i32..=2 {
            let cx = n1 as f64 * a1[0] + n2 as f64 * a2[0];
            let cy = n1 as f64 * a1[1] + n2 as f64 * a2[1];
            // Vertices of the (possibly re-entrant) hexagon.
            let tr = [cx + lc, cy + h / 2.0];
            let top = [cx, cy + h / 2.0 + ls];
            let tl = [cx - lc, cy + h / 2.0];
            let bl = [cx - lc, cy - h / 2.0];
            let bot = [cx, cy - h / 2.0 - ls];
            let br = [cx + lc, cy - h / 2.0];
            let verts = [tr, top, tl, bl, bot, br];
            for k in 0..6 {
                raw.push(RawMember2d {
                    a: verts[k],
                    b: verts[(k + 1) % 6],
                    e: 1.0,
                    area: t,
                    inertia: t * t * t / 12.0,
                });
            }
        }
    }
    (raw, a1, a2, h)
}

/// Regular honeycomb, `theta = 30 deg`, `h = l`.
pub fn regular_honeycomb(t: f64, l: f64, kind: Kind) -> (Vec<RawMember2d>, [f64; 2], [f64; 2]) {
    let (raw, a1, a2, _) = honeycomb(t, l, 1.0, std::f64::consts::PI / 6.0, kind);
    (raw, a1, a2)
}

/// Triangular lattice of pin-jointed bars, side `l`, bar thickness `t`.
pub fn triangular(t: f64, l: f64) -> (Vec<RawMember2d>, [f64; 2], [f64; 2]) {
    let a1 = [l, 0.0];
    let a2 = [0.5 * l, (3.0f64).sqrt() / 2.0 * l];
    let mut raw = Vec::new();
    for n1 in -2i32..=2 {
        for n2 in -2i32..=2 {
            let p = [
                n1 as f64 * a1[0] + n2 as f64 * a2[0],
                n1 as f64 * a1[1] + n2 as f64 * a2[1],
            ];
            // Three nearest-neighbour directions; images cover the rest.
            for d in [a1, a2, [a1[0] - a2[0], a1[1] - a2[1]]] {
                raw.push(RawMember2d {
                    a: p,
                    b: [p[0] + d[0], p[1] + d[1]],
                    e: 1.0,
                    area: t,
                    inertia: 0.0,
                });
            }
        }
    }
    (raw, a1, a2)
}

/// Face-centred-cubic nearest-neighbour network (the octet truss), cubic
/// lattice constant `a`, square strut side `t`. This is stretch-dominated.
pub fn octet_truss(t: f64, a: f64) -> (Vec<RawMember3d>, [f64; 3], [f64; 3], [f64; 3]) {
    let basis = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5 * a, 0.5 * a],
        [0.5 * a, 0.0, 0.5 * a],
        [0.5 * a, 0.5 * a, 0.0],
    ];
    let lat = [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]];
    let mut pts = Vec::new();
    for n1 in -1i32..=1 {
        for n2 in -1i32..=1 {
            for n3 in -1i32..=1 {
                for b in &basis {
                    pts.push([
                        b[0] + n1 as f64 * a,
                        b[1] + n2 as f64 * a,
                        b[2] + n3 as f64 * a,
                    ]);
                }
            }
        }
    }
    let l = a / 2.0f64.sqrt();
    let tol = 1e-6 * l;
    let mut raw = Vec::new();
    for i in 0..pts.len() {
        for j in (i + 1)..pts.len() {
            let d = [
                pts[j][0] - pts[i][0],
                pts[j][1] - pts[i][1],
                pts[j][2] - pts[i][2],
            ];
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            if (r - l).abs() < tol {
                raw.push(RawMember3d {
                    a: pts[i],
                    b: pts[j],
                    e: 1.0,
                    area: t * t,
                    iy: t.powi(4) / 12.0,
                    iz: t.powi(4) / 12.0,
                    j: 0.1406 * t.powi(4),
                });
            }
        }
    }
    (raw, lat[0], lat[1], lat[2])
}

/// Simple-cubic lattice of rigidly jointed beams, lattice constant `a`,
/// square strut side `t`.
pub fn cubic_frame(t: f64, a: f64) -> (Vec<RawMember3d>, [f64; 3], [f64; 3], [f64; 3]) {
    let lat = [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]];
    let mut raw = Vec::new();
    for n1 in -1i32..=1 {
        for n2 in -1i32..=1 {
            for n3 in -1i32..=1 {
                let p = [n1 as f64 * a, n2 as f64 * a, n3 as f64 * a];
                for d in lat {
                    raw.push(RawMember3d {
                        a: p,
                        b: [p[0] + d[0], p[1] + d[1], p[2] + d[2]],
                        e: 1.0,
                        area: t * t,
                        iy: t.powi(4) / 12.0,
                        iz: t.powi(4) / 12.0,
                        j: 0.1406 * t.powi(4),
                    });
                }
            }
        }
    }
    (raw, lat[0], lat[1], lat[2])
}

/// Square lattice of rigidly jointed beams, side `l`, thickness `t`.
pub fn square_frame(t: f64, l: f64) -> (Vec<RawMember2d>, [f64; 2], [f64; 2]) {
    let a1 = [l, 0.0];
    let a2 = [0.0, l];
    let mut raw = Vec::new();
    for n1 in -2i32..=2 {
        for n2 in -2i32..=2 {
            let p = [n1 as f64 * l, n2 as f64 * l];
            for d in [a1, a2] {
                raw.push(RawMember2d {
                    a: p,
                    b: [p[0] + d[0], p[1] + d[1]],
                    e: 1.0,
                    area: t,
                    inertia: t * t * t / 12.0,
                });
            }
        }
    }
    (raw, a1, a2)
}
