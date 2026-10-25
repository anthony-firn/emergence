//! The continuum rung: the Christoffel equation and the validity
//! indicator that decides when elasticity may replace the lattice.
//!
//! For a unit direction `n`, the acoustic sound speeds follow from the
//! eigenvalues of `ρ v² = C_{ijkl} n_j n_l`, where `C` is the Cauchy-Born
//! stiffness from [`crate::elastic`]. In the long-wavelength limit the
//! atomistic acoustic branches must approach these speeds; the fold is
//! certified by the relative dispersion error
//! `(ω_lattice − v k)/ω_lattice`, which grows with `k·a` and marks the
//! point where the continuum description must be unfolded back to phonons.

use crate::elastic::voigt;
use crate::linalg::sym_eigenvalues;

/// Normalise a direction, falling back to `(1,0,0)` for a zero vector.
pub fn unit(n: [f64; 3]) -> [f64; 3] {
    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if len < 1e-300 {
        [1.0, 0.0, 0.0]
    } else {
        [n[0] / len, n[1] / len, n[2] / len]
    }
}

/// Christoffel matrix `Γ_ik = C_ijkl n_j n_l / ρ` for a unit direction.
pub fn christoffel_matrix(c: &[[f64; 6]; 6], rho: f64, n: [f64; 3]) -> [[f64; 3]; 3] {
    let n = unit(n);
    let mut g = [[0.0; 3]; 3];
    for i in 0..3 {
        for k in 0..3 {
            let mut s = 0.0;
            for j in 0..3 {
                for l in 0..3 {
                    s += c[voigt(i, j)][voigt(k, l)] * n[j] * n[l];
                }
            }
            g[i][k] = s / rho;
        }
    }
    g
}

/// The three sound speeds along `n`, ascending, from `ρ v² = ` eigenvalues
/// of the Christoffel matrix.
pub fn christoffel_velocities(c: &[[f64; 6]; 6], rho: f64, n: [f64; 3]) -> [f64; 3] {
    let g = christoffel_matrix(c, rho, n);
    let flat = [
        g[0][0], g[0][1], g[0][2], g[1][0], g[1][1], g[1][2], g[2][0], g[2][1], g[2][2],
    ];
    let e = sym_eigenvalues(&flat, 3);
    [
        e[0].max(0.0).sqrt(),
        e[1].max(0.0).sqrt(),
        e[2].max(0.0).sqrt(),
    ]
}

/// Relative dispersion error `(ω_lattice − v k)/ω_lattice`: the fraction of
/// the lattice frequency that continuum elasticity fails to reproduce.
/// Zero as `k → 0`; order `(k a)²` for a smooth branch.
pub fn dispersion_error(omega_lattice: f64, v_continuum: f64, k: f64) -> f64 {
    if omega_lattice == 0.0 {
        return 0.0;
    }
    (omega_lattice - v_continuum * k) / omega_lattice
}

/// Magnitude of `k` for a wavevector.
pub fn k_norm(k: [f64; 3]) -> f64 {
    (k[0] * k[0] + k[1] * k[1] + k[2] * k[2]).sqrt()
}
