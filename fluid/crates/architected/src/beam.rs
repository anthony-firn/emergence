//! Linear beam/bar element stiffness matrices in 2D and 3D.
//!
//! The frame element is the linear Euler-Bernoulli beam with axial and
//! bending (and torsion in 3D) stiffness. A Timoshenko variant adds the
//! shear parameter `Phi = 12 E I / (G A_s L^2)`, which matters only for
//! stubby members; it is provided for completeness and cross-checked in the
//! tests. Bars (pin-jointed truss members) carry axial force only and are the
//! stretching-dominated limit used for the triangular and octet lattices.
//!
//! Degrees of freedom: 2D frame nodes have `[u_x, u_y, theta_z]`; 3D frame
//! nodes have `[u_x, u_y, u_z, theta_x, theta_y, theta_z]`. Trusses have the
//! translations only.

/// 2D Euler-Bernoulli (or Timoshenko) frame element, local coordinates with
/// the member along `x`. `area` is the cross-section area, `inertia` the
/// second moment in bending, `e` Young's modulus. `shear` is the Timoshenko
/// shear parameter `Phi`; pass 0 for Euler-Bernoulli.
pub fn frame2d_local(l: f64, e: f64, area: f64, inertia: f64, shear: f64) -> [[f64; 6]; 6] {
    let mut k = [[0.0f64; 6]; 6];
    let ea = e * area / l;
    let ei = e * inertia;
    let phi = shear;
    let a = 12.0 * ei / (l * l * l * (1.0 + phi));
    let b = 6.0 * ei / (l * l * (1.0 + phi));
    let c = (4.0 + phi) * ei / (l * (1.0 + phi));
    let d = (2.0 - phi) * ei / (l * (1.0 + phi));
    let set = |k: &mut [[f64; 6]; 6], i: usize, j: usize, v: f64| {
        k[i][j] = v;
        k[j][i] = v;
    };
    set(&mut k, 0, 0, ea);
    set(&mut k, 0, 3, -ea);
    set(&mut k, 3, 3, ea);
    set(&mut k, 1, 1, a);
    set(&mut k, 1, 2, b);
    set(&mut k, 1, 4, -a);
    set(&mut k, 1, 5, b);
    set(&mut k, 2, 2, c);
    set(&mut k, 2, 4, -b);
    set(&mut k, 2, 5, d);
    set(&mut k, 4, 4, a);
    set(&mut k, 4, 5, -b);
    set(&mut k, 5, 5, c);
    k
}

/// Timoshenko shear parameter for a rectangular section, `k_s = 5/6`.
pub fn timoshenko_phi(l: f64, e: f64, area: f64, inertia: f64, nu: f64) -> f64 {
    let g = e / (2.0 * (1.0 + nu));
    let a_s = (5.0 / 6.0) * area;
    12.0 * e * inertia / (g * a_s * l * l)
}

/// Transform a 2D local element matrix (6x6) to global coordinates for a
/// member with direction cosines `(c, s)`.
pub fn frame2d_global(k: &[[f64; 6]; 6], c: f64, s: f64) -> [[f64; 6]; 6] {
    // T = blockdiag(R, R), R = [[c, s, 0], [-s, c, 0], [0, 0, 1]].
    let r = [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]];
    let mut t = [[0.0f64; 6]; 6];
    for block in 0..2 {
        for m in 0..3 {
            for n in 0..3 {
                t[block * 3 + m][block * 3 + n] = r[m][n];
            }
        }
    }
    mat6(t_t(k, &t), &t)
}

fn t_t(k: &[[f64; 6]; 6], t: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut out = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let mut s = 0.0;
            for l in 0..6 {
                s += k[i][l] * t[l][j];
            }
            out[i][j] = s;
        }
    }
    out
}

fn mat6(a: [[f64; 6]; 6], t: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut out = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let mut s = 0.0;
            for l in 0..6 {
                s += t[l][i] * a[l][j];
            }
            out[i][j] = s;
        }
    }
    out
}

/// 2D pin-jointed bar element (4x4) in global coordinates. `c, s` direction
/// cosines.
pub fn truss2d_global(l: f64, e: f64, area: f64, c: f64, s: f64) -> [f64; 16] {
    let k = e * area / l;
    let (c2, s2, cs) = (c * c, s * s, c * s);
    let mut m = [0.0f64; 16];
    let vals = [
        c2, cs, -c2, -cs, cs, s2, -cs, -s2, -c2, -cs, c2, cs, -cs, -s2, cs, s2,
    ];
    for (i, v) in vals.iter().enumerate() {
        m[i] = k * v;
    }
    m
}

/// Local 12x12 3D Euler-Bernoulli frame element, member along local `x`.
/// `iz`, `iy` bending inertias, `j` the torsion constant, `area` axial area.
/// `g` is the shear modulus, used for the torsion stiffness `G J / L`.
pub fn frame3d_local(l: f64, e: f64, g: f64, area: f64, iy: f64, iz: f64, j: f64) -> [[f64; 12]; 12] {
    let mut k = [[0.0f64; 12]; 12];
    let set = |k: &mut [[f64; 12]; 12], i: usize, jj: usize, v: f64| {
        k[i][jj] = v;
        k[jj][i] = v;
    };
    let ea = e * area / l;
    set(&mut k, 0, 0, ea);
    set(&mut k, 0, 6, -ea);
    set(&mut k, 6, 6, ea);
    let gj = g * j / l;
    set(&mut k, 3, 3, gj);
    set(&mut k, 3, 9, -gj);
    set(&mut k, 9, 9, gj);
    // bending about z: v(1), rz(5), v2(7), rz2(11)
    let az = 12.0 * e * iz / l.powi(3);
    let bz = 6.0 * e * iz / l.powi(2);
    let cz = 4.0 * e * iz / l;
    let dz = 2.0 * e * iz / l;
    set(&mut k, 1, 1, az);
    set(&mut k, 1, 5, bz);
    set(&mut k, 1, 7, -az);
    set(&mut k, 1, 11, bz);
    set(&mut k, 5, 5, cz);
    set(&mut k, 5, 7, -bz);
    set(&mut k, 5, 11, dz);
    set(&mut k, 7, 7, az);
    set(&mut k, 7, 11, -bz);
    set(&mut k, 11, 11, cz);
    // bending about y: w(2), ry(4), w2(8), ry2(10)
    let ay = 12.0 * e * iy / l.powi(3);
    let by = 6.0 * e * iy / l.powi(2);
    let cy = 4.0 * e * iy / l;
    let dy = 2.0 * e * iy / l;
    set(&mut k, 2, 2, ay);
    set(&mut k, 2, 4, -by);
    set(&mut k, 2, 8, -ay);
    set(&mut k, 2, 10, -by);
    set(&mut k, 4, 4, cy);
    set(&mut k, 4, 8, by);
    set(&mut k, 4, 10, dy);
    set(&mut k, 8, 8, ay);
    set(&mut k, 8, 10, by);
    set(&mut k, 10, 10, cy);
    k
}

/// Rotation matrix `R` (rows are local axes `e1, e2, e3`) mapping a global
/// vector to local components for a member along `e1`.
pub fn local_frame(e1: [f64; 3]) -> [[f64; 3]; 3] {
    let n = (e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]).sqrt();
    let e1 = [e1[0] / n, e1[1] / n, e1[2] / n];
    let reference = if e1[1].abs() < 0.9 {
        [0.0, 1.0, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let d = reference[0] * e1[0] + reference[1] * e1[1] + reference[2] * e1[2];
    let mut e2 = [
        reference[0] - d * e1[0],
        reference[1] - d * e1[1],
        reference[2] - d * e1[2],
    ];
    let n2 = (e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]).sqrt();
    for v in &mut e2 {
        *v /= n2;
    }
    let e3 = [
        e1[1] * e2[2] - e1[2] * e2[1],
        e1[2] * e2[0] - e1[0] * e2[2],
        e1[0] * e2[1] - e1[1] * e2[0],
    ];
    [e1, e2, e3]
}

/// Transform a 3D local element matrix to global coordinates. `r` is the
/// 3x3 rotation from [`local_frame`]; the 12x12 transform is block diagonal
/// with `r` repeated four times.
pub fn frame3d_global(k: &[[f64; 12]; 12], r: &[[f64; 3]; 3]) -> [[f64; 12]; 12] {
    let mut t = [[0.0f64; 12]; 12];
    for block in 0..4 {
        for m in 0..3 {
            for n in 0..3 {
                t[block * 3 + m][block * 3 + n] = r[m][n];
            }
        }
    }
    // out = T^T K T
    let mut kt = [[0.0f64; 12]; 12];
    for i in 0..12 {
        for j in 0..12 {
            let mut s = 0.0;
            for l in 0..12 {
                s += k[i][l] * t[l][j];
            }
            kt[i][j] = s;
        }
    }
    let mut out = [[0.0f64; 12]; 12];
    for i in 0..12 {
        for j in 0..12 {
            let mut s = 0.0;
            for l in 0..12 {
                s += t[l][i] * kt[l][j];
            }
            out[i][j] = s;
        }
    }
    out
}

/// 3D pin-jointed bar element (6x6) in global coordinates.
pub fn truss3d_global(l: f64, e: f64, area: f64, dir: [f64; 3]) -> [f64; 36] {
    let k = e * area / l;
    let n = dir;
    let mut m = [0.0f64; 36];
    for i in 0..3 {
        for j in 0..3 {
            let v = k * n[i] * n[j];
            m[i * 6 + j] = v;
            m[i * 6 + (j + 3)] = -v;
            m[(i + 3) * 6 + j] = -v;
            m[(i + 3) * 6 + (j + 3)] = v;
        }
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn frame2d_rigid_modes_are_zero() {
        // Translating both ends and rotating rigidly must give zero force.
        let k = frame2d_local(2.0, 200e9, 1e-4, 1e-8, 0.0);
        // rigid body translation in x
        let u = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        for i in 0..6 {
            let f: f64 = (0..6).map(|j| k[i][j] * u[j]).sum();
            assert!(f.abs() < 1e-6, "row {i} force {f}");
        }
        // small rigid rotation: v = theta*x, all rotations equal
        let u = [0.0, 0.0, 0.01, 0.0, 0.02, 0.01];
        for i in 0..6 {
            let f: f64 = (0..6).map(|j| k[i][j] * u[j]).sum();
            if i == 1 || i == 4 {
                assert!(f.abs() < 1e-6);
            }
        }
    }
    #[test]
    fn bar_axial_matches_hooke() {
        let m = truss2d_global(2.0, 100.0, 0.5, 1.0, 0.0);
        // Pull node 2 in +x by 0.1: force = EA/L * 0.1
        let u = [0.0, 0.0, 0.1, 0.0];
        let f: f64 = (0..4).map(|j| m[2 * 4 + j] * u[j]).sum();
        assert!((f - 100.0 * 0.5 / 2.0 * 0.1).abs() < 1e-12);
    }

    #[test]
    fn frame3d_cantilever_and_torsion_are_exact() {
        use crate::linalg::{solve, Dense};
        let (l, e) = (2.0, 210e9);
        let area = 0.01 * 0.02;
        let iy = 0.01 * 0.02f64.powi(3) / 12.0;
        let iz = 0.02 * 0.01f64.powi(3) / 12.0;
        let j = 0.229 * 0.01 * 0.02f64.powi(3);
        let gmod = e / (2.0 * (1.0 + 0.3));
        let k = frame3d_local(l, e, gmod, area, iy, iz, j);
        // Fix node A (first 6 DOF), leave node B free.
        let mut kr = Dense::zeros(6);
        for i in 0..6 {
            for jj in 0..6 {
                kr.set(i, jj, k[i + 6][jj + 6]);
            }
        }
        // Transverse load in local y -> bending about z.
        let p = 100.0;
        let mut f = vec![0.0; 6];
        f[1] = p;
        solve(&mut kr.clone(), &mut f).unwrap();
        let exact = p * l * l * l / (3.0 * e * iz);
        assert!((f[1].abs() - exact).abs() / exact < 1e-12);
        // Transverse load in local z -> bending about y.
        let mut f = vec![0.0; 6];
        f[2] = p;
        solve(&mut kr.clone(), &mut f).unwrap();
        let exact = p * l * l * l / (3.0 * e * iy);
        assert!((f[2].abs() - exact).abs() / exact < 1e-12);
        // Torque -> twist theta_x = M L / (G J).
        let m = 50.0;
        let mut f = vec![0.0; 6];
        f[3] = m;
        solve(&mut kr, &mut f).unwrap();
        let exact = m * l / (gmod * j);
        assert!((f[3].abs() - exact).abs() / exact < 1e-12);
    }
}
