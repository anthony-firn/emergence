//! Elastic-wave dispersion of a 2D beam lattice by Bloch-Floquet analysis.
//!
//! For an infinite lattice, a Bloch wave has the form `u(x) = U(x) exp(i k.x)`
//! with `U` periodic. Only one unit cell is assembled. A member whose image
//! endpoint is displaced by a lattice vector `d` relative to the canonical
//! node carries the phase `exp(i k.d)`:
//!
//! ```text
//!   T = diag(1, e^{i k.d}) ,   K(k) = T^H K_e T ,
//!   K(k) phi = omega^2 M phi .
//! ```
//!
//! `K(k)` is Hermitian. Writing `K = A + i B` with `A` real symmetric and `B`
//! real antisymmetric, the Hermitian problem is equivalent to the real
//! symmetric generalised problem
//!
//! ```text
//!   [ A  -B ] [p]        [M 0] [p]
//!   [ B   A ] [q]  = w^2 [0 M] [q] ,
//! ```
//!
//! which the crate's real dense solver handles. At long wavelength the
//! acoustic branches have slope `sqrt(C*/rho*)`, the continuum sound speed,
//! and the test compares them with the homogenised stiffness.

use crate::beam::{frame2d_global, frame2d_local};
use crate::cell::{Cell2d, Kind};
use crate::linalg::{generalized_eigen, Dense};

/// Bloch wavevector `(k1, k2)` in rad / length.
#[derive(Clone, Copy, Debug)]
pub struct Wavevector {
    pub k1: f64,
    pub k2: f64,
}

/// Frequencies [rad/s] of the lowest `nbands` branches at `k`.
pub fn bands_at(cell: &Cell2d, k: Wavevector, rho: f64, nbands: usize) -> Vec<f64> {
    let (krep, mrep) = real_representation(cell, k, rho);
    let (vals, _x) = generalized_eigen(&krep, &mrep).expect("singular mass matrix");
    vals.into_iter()
        .filter(|v| *v > 1e-12)
        .take(nbands)
        .map(|v| v.sqrt())
        .collect()
}

/// Build the real symmetric representation `[[A,-B],[B,A]]` and the
/// block-diagonal mass matrix.
fn real_representation(cell: &Cell2d, k: Wavevector, rho: f64) -> (Dense, Dense) {
    assert_eq!(cell.kind, Kind::Frame, "Bloch bands implemented for frames");
    let dof = 3;
    let n = cell.nodes.len() * dof;
    let mut a = Dense::zeros(n);
    let mut b = Dense::zeros(n);
    let mut m = Dense::zeros(n);
    for mem in &cell.members {
        let xa = cell.physical(mem.a, mem.a_img);
        let xb = cell.physical(mem.b, mem.b_img);
        let dx = xb[0] - xa[0];
        let dy = xb[1] - xa[1];
        let l = (dx * dx + dy * dy).sqrt();
        let (c, s) = (dx / l, dy / l);
        let kl = frame2d_local(l, mem.e, mem.area, mem.inertia, 0.0);
        let kg = frame2d_global(&kl, c, s);
        let idx = [
            mem.a * dof,
            mem.a * dof + 1,
            mem.a * dof + 2,
            mem.b * dof,
            mem.b * dof + 1,
            mem.b * dof + 2,
        ];
        // Physical image displacement d = n1 a1 + n2 a2.
        let dimg = [
            mem.b_img[0] - mem.a_img[0],
            mem.b_img[1] - mem.a_img[1],
        ];
        let dphys = [
            dimg[0] as f64 * cell.a1[0] + dimg[1] as f64 * cell.a2[0],
            dimg[0] as f64 * cell.a1[1] + dimg[1] as f64 * cell.a2[1],
        ];
        let phi = k.k1 * dphys[0] + k.k2 * dphys[1];
        let (cf, sf) = (phi.cos(), phi.sin());
        for i in 0..6 {
            for j in 0..6 {
                let bi = i / 3;
                let bj = j / 3;
                let val = kg[i][j];
                if bi == 0 && bj == 0 {
                    a.add(idx[i], idx[j], val);
                } else if bi == 1 && bj == 1 {
                    a.add(idx[i], idx[j], val);
                } else if bi == 0 && bj == 1 {
                    // K_ab = kg_ab e^{i phi}
                    a.add(idx[i], idx[j], val * cf);
                    b.add(idx[i], idx[j], val * sf);
                } else {
                    // K_ba = kg_ab^T e^{-i phi}: Re = val cos phi,
                    // Im = -val sin phi.
                    a.add(idx[i], idx[j], val * cf);
                    b.add(idx[i], idx[j], -val * sf);
                }
            }
        }
        let mass = rho * mem.area * l;
        for t in 0..3 {
            m.add(idx[t], idx[t], 0.5 * mass);
            m.add(idx[t + 3], idx[t + 3], 0.5 * mass);
        }
    }
    // Real symmetric representation of the Hermitian pencil.
    let mut krep = Dense::zeros(2 * n);
    let mut mrep = Dense::zeros(2 * n);
    for i in 0..n {
        for j in 0..n {
            krep.set(i, j, a.get(i, j));
            krep.set(i, n + j, -b.get(i, j));
            krep.set(n + i, j, b.get(i, j));
            krep.set(n + i, n + j, a.get(i, j));
            mrep.set(i, j, m.get(i, j));
            mrep.set(n + i, n + j, m.get(i, j));
        }
    }
    (krep, mrep)
}

/// Dispersion along a straight `k` path: `(kx, ky, omega)` triples.
pub fn dispersion_path(
    cell: &Cell2d,
    from: [f64; 2],
    to: [f64; 2],
    steps: usize,
    rho: f64,
    nbands: usize,
) -> Vec<(f64, f64, f64)> {
    let mut out = Vec::new();
    for i in 0..=steps {
        let f = i as f64 / steps as f64;
        let k = Wavevector {
            k1: from[0] + f * (to[0] - from[0]),
            k2: from[1] + f * (to[1] - from[1]),
        };
        for w in bands_at(cell, k, rho, nbands) {
            out.push((k.k1, k.k2, w));
        }
    }
    out
}

/// Acoustic slope `d omega / d k` of the lowest branch along `direction`.
pub fn acoustic_slope(cell: &Cell2d, direction: [f64; 2], rho: f64, dk: f64) -> f64 {
    let w0 = bands_at(
        cell,
        Wavevector { k1: 0.0, k2: 0.0 },
        rho,
        1,
    );
    let w1 = bands_at(
        cell,
        Wavevector {
            k1: direction[0] * dk,
            k2: direction[1] * dk,
        },
        rho,
        1,
    );
    (w1[0] - w0[0]) / dk
}
