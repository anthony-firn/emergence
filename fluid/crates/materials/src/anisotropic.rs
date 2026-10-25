//! Anisotropic continuum rung: single-crystal elastic constants, symmetry
//! classes, orientation (Bond) rotation, the compliance, and the directional
//! elastic-wave properties they predict.
//!
//! The fine rung is the atomistic crystal ([`crate::elastic`]). Above it sits
//! the homogeneous elastic continuum of a *real* crystal, described by the
//! 6x6 Voigt stiffness `C_IJ` in the crystal frame. The rungs are linked by
//! the Cauchy-Born rule; here we take `C` as measured and propagate it:
//!
//! - symmetry classes: cubic (3 constants), hexagonal (5) and the general
//!   triclinic case (21, supplied as a full 6x6);
//! - rotation of the stiffness to any crystal orientation with the Bond
//!   matrix (equivalent to contracting the fourth-order tensor);
//! - the compliance `S = C⁻¹`;
//! - directional properties: the Christoffel wave speeds and polarisations,
//!   Young's modulus `E(n) = 1/S'_{11}` and the Zener anisotropy
//!   `A = 2C44/(C11 − C12)`.
//!
//! Constants are room-temperature single-crystal values. The primary source
//! is Simmons & Wang, *Single Crystal Elastic Constants and Calculated
//! Aggregate Properties: A Handbook* (MIT Press, 1971), with the silicon
//! values from Hopcroft, Nix & Kenny (2010) and densities from the CRC
//! Handbook. Entries that are rounded handbook estimates carry
//! `estimate = true`. Absolute units are SI (Pa, kg/m³).

use crate::elastic::voigt;
use crate::linalg::{solve, sym_eigen};

/// A 6x6 Voigt stiffness or compliance, in Pa (stiffness) or Pa⁻¹
/// (compliance).
pub type M6 = [[f64; 6]; 6];
type V3 = [f64; 3];

// ---------------------------------------------------------------------------
// Symmetry classes and stiffness builders
// ---------------------------------------------------------------------------

/// Crystal symmetry class. Only the symmetry classes needed to build the
/// stiffness from independent constants are named; a triclinic tensor is
/// given in full.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Symmetry {
    Cubic,
    Hexagonal,
    Triclinic,
}

/// Cubic stiffness (three independent constants, all in Pa).
pub fn cubic(c11: f64, c12: f64, c44: f64) -> M6 {
    let mut c = [[0.0; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = if i == j { c11 } else { c12 };
        }
    }
    for i in 3..6 {
        c[i][i] = c44;
    }
    c
}

/// Hexagonal (transversely isotropic, `z` the six-fold axis) stiffness with
/// five independent constants, all in Pa:
/// `C11, C12, C13, C33, C44` and `C66 = (C11 − C12)/2`.
pub fn hexagonal(c11: f64, c12: f64, c13: f64, c33: f64, c44: f64) -> M6 {
    let c66 = 0.5 * (c11 - c12);
    let mut c = [[0.0; 6]; 6];
    c[0][0] = c11;
    c[1][1] = c11;
    c[2][2] = c33;
    c[0][1] = c12;
    c[1][0] = c12;
    c[0][2] = c13;
    c[2][0] = c13;
    c[1][2] = c13;
    c[2][1] = c13;
    c[3][3] = c44;
    c[4][4] = c44;
    c[5][5] = c66;
    c
}

/// General anisotropic (triclinic) stiffness: symmetrise the supplied full
/// 6x6. Use this for any symmetry class not covered by the builders.
pub fn triclinic(mut c: M6) -> M6 {
    for i in 0..6 {
        for j in i + 1..6 {
            let v = 0.5 * (c[i][j] + c[j][i]);
            c[i][j] = v;
            c[j][i] = v;
        }
    }
    c
}

// ---------------------------------------------------------------------------
// Real single crystals (room temperature)
// ---------------------------------------------------------------------------

/// Room-temperature single-crystal elastic data.
#[derive(Clone, Copy, Debug)]
pub struct CrystalMaterial {
    pub name: &'static str,
    pub symmetry: Symmetry,
    /// Mass density, kg/m³.
    pub density: f64,
    /// Voigt stiffness in the crystal frame, Pa.
    pub c: M6,
    pub source: &'static str,
    /// True where the quoted constants are rounded handbook estimates.
    pub estimate: bool,
}

const SW: &str = "Simmons & Wang, Single Crystal Elastic Constants and Calculated \
Aggregate Properties: A Handbook (MIT Press, 1971)";
const HOPCROFT: &str = "Hopcroft, Nix & Kenny, J. Microelectromech. Syst. 19, 229 (2010)";

fn cubic_mat(
    name: &'static str,
    c11: f64,
    c12: f64,
    c44: f64,
    density: f64,
    source: &'static str,
    estimate: bool,
) -> CrystalMaterial {
    CrystalMaterial {
        name,
        symmetry: Symmetry::Cubic,
        density,
        c: cubic(c11, c12, c44),
        source,
        estimate,
    }
}

fn hex_mat(
    name: &'static str,
    c11: f64,
    c12: f64,
    c13: f64,
    c33: f64,
    c44: f64,
    density: f64,
    source: &'static str,
    estimate: bool,
) -> CrystalMaterial {
    CrystalMaterial {
        name,
        symmetry: Symmetry::Hexagonal,
        density,
        c: hexagonal(c11, c12, c13, c33, c44),
        source,
        estimate,
    }
}

/// The sourced room-temperature single-crystal table. Stiffnesses are in Pa
/// (the quoted handbook values are GPa); densities in kg/m³.
pub fn materials() -> Vec<CrystalMaterial> {
    vec![
        // Cubic: Si, Ge, diamond, Cu, Al, Ni, alpha-Fe, W, NaCl.
        cubic_mat("Si", 165.7e9, 63.9e9, 79.6e9, 2329.0, HOPCROFT, false),
        cubic_mat("Ge", 128.5e9, 48.3e9, 66.7e9, 5323.0, SW, false),
        cubic_mat("diamond", 1076.0e9, 125.0e9, 576.0e9, 3515.0, SW, false),
        cubic_mat("Cu", 168.4e9, 121.4e9, 75.4e9, 8960.0, SW, false),
        cubic_mat("Al", 107.3e9, 60.9e9, 28.3e9, 2700.0, SW, false),
        cubic_mat("Ni", 246.5e9, 147.3e9, 124.7e9, 8900.0, SW, false),
        cubic_mat("alpha-Fe", 231.4e9, 134.5e9, 116.4e9, 7874.0, SW, false),
        cubic_mat("W", 522.4e9, 204.4e9, 160.8e9, 19250.0, SW, false),
        cubic_mat("NaCl", 48.7e9, 12.4e9, 12.7e9, 2165.0, SW, false),
        // Hexagonal: Mg, Ti, Zn.
        hex_mat(
            "Mg", 59.7e9, 25.8e9, 21.7e9, 61.7e9, 16.4e9, 1738.0, SW, false,
        ),
        hex_mat(
            "Ti", 162.4e9, 92.0e9, 69.0e9, 180.7e9, 46.7e9, 4507.0, SW, false,
        ),
        hex_mat(
            "Zn", 162.9e9, 36.3e9, 52.7e9, 62.9e9, 39.2e9, 7140.0, SW, true,
        ),
    ]
}

/// Look up a material by name (case sensitive, e.g. `"Cu"`, `"alpha-Fe"`).
pub fn material(name: &str) -> Option<CrystalMaterial> {
    materials().into_iter().find(|m| m.name == name)
}

// ---------------------------------------------------------------------------
// Orientation: Bond matrices and stiffness rotation
// ---------------------------------------------------------------------------

/// Rotation matrix `R` (row-major) about a unit `axis` by `angle`, right-handed.
pub fn rotation_axis_angle(axis: V3, angle: f64) -> [[f64; 3]; 3] {
    let n = unit(axis);
    let (s, c) = angle.sin_cos();
    let t = 1.0 - c;
    [
        [
            t * n[0] * n[0] + c,
            t * n[0] * n[1] - s * n[2],
            t * n[0] * n[2] + s * n[1],
        ],
        [
            t * n[0] * n[1] + s * n[2],
            t * n[1] * n[1] + c,
            t * n[1] * n[2] - s * n[0],
        ],
        [
            t * n[0] * n[2] - s * n[1],
            t * n[1] * n[2] + s * n[0],
            t * n[2] * n[2] + c,
        ],
    ]
}

/// Rotation whose first row is `n` (so it maps the crystal direction `n` onto
/// the lab `x` axis). Rows are an orthonormal frame.
pub fn rotation_align_x(n: V3) -> [[f64; 3]; 3] {
    let n = unit(n);
    let helper = if n[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let t1 = unit(cross(n, helper));
    let t2 = cross(n, t1);
    [n, t1, t2]
}

/// Bond matrix for the second-rank tensor transformation `A' = R A Rᵀ`
/// (stress convention): the Voigt 6-vector of `A'` is `M` times the Voigt
/// 6-vector of `A`. `M` is the standard Bond matrix (Nye, *Physical
/// Properties of Crystals*).
pub fn bond_matrix(r: &[[f64; 3]; 3]) -> M6 {
    let a = r;
    let mut m = [[0.0; 6]; 6];
    for i in 0..3 {
        m[i][0] = a[i][0] * a[i][0];
        m[i][1] = a[i][1] * a[i][1];
        m[i][2] = a[i][2] * a[i][2];
        m[i][3] = 2.0 * a[i][1] * a[i][2];
        m[i][4] = 2.0 * a[i][0] * a[i][2];
        m[i][5] = 2.0 * a[i][0] * a[i][1];
    }
    let pair = |i: usize, k: usize, j: usize, l: usize| a[i][j] * a[k][l] + a[i][l] * a[k][j];
    m[3][0] = a[1][0] * a[2][0];
    m[3][1] = a[1][1] * a[2][1];
    m[3][2] = a[1][2] * a[2][2];
    m[3][3] = pair(1, 2, 1, 2);
    m[3][4] = pair(1, 2, 0, 2);
    m[3][5] = pair(1, 2, 0, 1);
    m[4][0] = a[2][0] * a[0][0];
    m[4][1] = a[2][1] * a[0][1];
    m[4][2] = a[2][2] * a[0][2];
    m[4][3] = pair(2, 0, 1, 2);
    m[4][4] = pair(2, 0, 0, 2);
    m[4][5] = pair(2, 0, 0, 1);
    m[5][0] = a[0][0] * a[1][0];
    m[5][1] = a[0][1] * a[1][1];
    m[5][2] = a[0][2] * a[1][2];
    m[5][3] = pair(0, 1, 1, 2);
    m[5][4] = pair(0, 1, 0, 2);
    m[5][5] = pair(0, 1, 0, 1);
    m
}

/// Rotate a stiffness to a new frame. Contracting the fourth-order tensor,
/// `C'_{ijkl} = R_{im}R_{jn}R_{ko}R_{lp}C_{mnop}`, with
/// `C_{ijkl} = c[voigt(ij)][voigt(kl)]`; the Bond-matrix form
/// [`rotate_stiffness_bond`] is equivalent and tested against it.
pub fn rotate_stiffness(c: &M6, r: &[[f64; 3]; 3]) -> M6 {
    let mut out = [[0.0; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    let mut s = 0.0;
                    for m in 0..3 {
                        for n in 0..3 {
                            for o in 0..3 {
                                for p in 0..3 {
                                    s += r[i][m]
                                        * r[j][n]
                                        * r[k][o]
                                        * r[l][p]
                                        * c[voigt(m, n)][voigt(o, p)];
                                }
                            }
                        }
                    }
                    out[voigt(i, j)][voigt(k, l)] = s;
                }
            }
        }
    }
    out
}

/// Rotate a stiffness with the Bond matrices: `C' = M C Mᵀ`. Equivalent to
/// [`rotate_stiffness`].
pub fn rotate_stiffness_bond(c: &M6, r: &[[f64; 3]; 3]) -> M6 {
    let m = bond_matrix(r);
    let mut tmp = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            tmp[i][j] = (0..6).map(|k| m[i][k] * c[k][j]).sum();
        }
    }
    let mut out = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            out[i][j] = (0..6).map(|k| tmp[i][k] * m[j][k]).sum();
        }
    }
    out
}

/// Inverse of a 6x6 matrix by Gauss elimination with six right-hand sides.
pub fn inverse(m: &M6) -> M6 {
    let mut inv = [[0.0; 6]; 6];
    for col in 0..6 {
        let mut rhs = vec![0.0; 6];
        rhs[col] = 1.0;
        let a: Vec<f64> = m.iter().flatten().copied().collect();
        let x = solve(a, rhs, 6);
        for row in 0..6 {
            inv[row][col] = x[row];
        }
    }
    inv
}

/// Compliance `S = C⁻¹`.
pub fn compliance(c: &M6) -> M6 {
    inverse(c)
}

// ---------------------------------------------------------------------------
// Directional properties
// ---------------------------------------------------------------------------

fn unit(n: V3) -> V3 {
    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if len < 1e-300 {
        [1.0, 0.0, 0.0]
    } else {
        [n[0] / len, n[1] / len, n[2] / len]
    }
}

fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Young's modulus along the crystal direction `n`: `E(n) = 1/S'_{11}`
/// where `S'` is the compliance in the frame with `n` along `x`.
pub fn youngs_modulus(c: &M6, n: V3) -> f64 {
    let r = rotation_align_x(n);
    let cr = rotate_stiffness(c, &r);
    let s = compliance(&cr);
    1.0 / s[0][0]
}

/// Zener anisotropy ratio for a cubic crystal, `A = 2C44/(C11 − C12)`.
/// `A = 1` is elastically isotropic.
pub fn zener_ratio(c11: f64, c12: f64, c44: f64) -> f64 {
    2.0 * c44 / (c11 - c12)
}

/// Christoffel matrix `Γ_{ik} = C_{ijkl} n_j n_l / ρ` in the crystal frame.
pub fn christoffel(c: &M6, rho: f64, n: V3) -> [[f64; 3]; 3] {
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

/// The three Christoffel wave speeds along `n`, ascending, and the
/// corresponding polarisation unit vectors. The polarisation returned for a
/// degenerate pair is one arbitrary orthonormal choice within the pair.
pub fn christoffel_modes(c: &M6, rho: f64, n: V3) -> ([f64; 3], [[f64; 3]; 3]) {
    let g = christoffel(c, rho, n);
    let flat = [
        g[0][0], g[0][1], g[0][2], g[1][0], g[1][1], g[1][2], g[2][0], g[2][1], g[2][2],
    ];
    let (e, v) = sym_eigen(&flat, 3);
    let speeds = [
        e[0].max(0.0).sqrt(),
        e[1].max(0.0).sqrt(),
        e[2].max(0.0).sqrt(),
    ];
    let pol = [
        [v[0][0], v[0][1], v[0][2]],
        [v[1][0], v[1][1], v[1][2]],
        [v[2][0], v[2][1], v[2][2]],
    ];
    (speeds, pol)
}

/// Longitudinal and (fast, slow) transverse speeds along `n` for a
/// polycrystal with isotropic shear/bulk moduli. Convenience for the fold in
/// [`crate::polycrystal`].
pub fn isotropic_speeds(bulk: f64, shear: f64, rho: f64) -> (f64, f64) {
    let vp = ((bulk + 4.0 * shear / 3.0) / rho).sqrt();
    let vs = (shear / rho).sqrt();
    (vp, vs)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol * b.abs().max(1.0) + 1e-6, "{a} vs {b}");
    }

    #[test]
    fn bond_and_tensor_rotation_agree() {
        let c = cubic(165.7e9, 63.9e9, 79.6e9);
        for r in [
            rotation_axis_angle([1.0, 1.0, 1.0], 0.7),
            rotation_axis_angle([0.0, 1.0, 0.0], 0.31),
            rotation_align_x([1.0, 1.0, 0.0]),
        ] {
            let a = rotate_stiffness(&c, &r);
            let b = rotate_stiffness_bond(&c, &r);
            let scale = c.iter().flatten().fold(0.0_f64, |m, &x| m.max(x.abs()));
            for i in 0..6 {
                for j in 0..6 {
                    assert!(
                        (a[i][j] - b[i][j]).abs() <= 1e-12 * scale,
                        "rotation mismatch at {i},{j}: {} vs {}",
                        a[i][j],
                        b[i][j]
                    );
                }
            }
        }
    }

    #[test]
    fn compliance_inverts_stiffness() {
        let c = cubic(107.3e9, 60.9e9, 28.3e9);
        let s = compliance(&c);
        for i in 0..6 {
            for j in 0..6 {
                let v: f64 = (0..6).map(|k| c[i][k] * s[k][j]).sum();
                let expect = if i == j { 1.0 } else { 0.0 };
                assert!((v - expect).abs() < 1e-12, "C*S != I at {i},{j}: {v}");
            }
        }
    }

    #[test]
    fn cubic_youngs_modulus_matches_closed_form() {
        // E[hkl]^{-1} = S11 - 2(S11 - S12 - S44/2)(l^2 m^2 + m^2 n^2 + n^2 l^2).
        let (c11, c12, c44) = (165.7e9, 63.9e9, 79.6e9);
        let c = cubic(c11, c12, c44);
        let s = compliance(&c);
        let (s11, s12, s44) = (s[0][0], s[0][1], s[3][3]);
        for (n, combos) in [
            ([1.0, 0.0, 0.0], 0.0),
            ([1.0, 1.0, 0.0], 0.25),
            ([1.0, 1.0, 1.0], 1.0 / 3.0),
        ] {
            let expected = 1.0 / (s11 - 2.0 * (s11 - s12 - 0.5 * s44) * combos);
            close(youngs_modulus(&c, n), expected, 1e-9);
        }
    }
}
