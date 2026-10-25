//! Benchmarks for the anisotropic continuum rung of the materials ladder
//! (`analysis/materials.md`). Directional properties of real single crystals
//! are checked against published values; rotation and polarisation identities
//! are checked exactly.

use materials::anisotropic::{
    christoffel_modes, material, rotate_stiffness, rotation_axis_angle, youngs_modulus, zener_ratio,
};

fn close_rel(actual: f64, expected: f64, tol: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tol * expected.abs().max(1e-300),
        "{what}: {actual:.6e} vs {expected:.6e}, rel tolerance {tol:.2e}"
    );
}

/// Silicon acoustic speeds (Hopcroft, Nix & Kenny 2010): [100] longitudinal
/// ~8433 m/s and transverse ~5843 m/s; [110] longitudinal ~9134 m/s.
#[test]
fn silicon_sound_speeds_match_hopcroft() {
    let si = material("Si").expect("Si in table");
    let (v100, _) = christoffel_modes(&si.c, si.density, [1.0, 0.0, 0.0]);
    // Ascending: two degenerate transverse modes, then the longitudinal one.
    close_rel(v100[0], 5843.0, 5e-3, "Si [100] T");
    close_rel(v100[1], 5843.0, 5e-3, "Si [100] T (degenerate)");
    close_rel(v100[2], 8433.0, 5e-3, "Si [100] L");

    let (v110, _) = christoffel_modes(&si.c, si.density, [1.0, 1.0, 0.0]);
    close_rel(v110[2], 9134.0, 5e-3, "Si [110] L");
    // The two transverse speeds differ: the [1-10] polarised shear and the
    // [001] polarised shear.
    close_rel(
        v110[0],
        ((si.c[0][0] - si.c[0][1]) / (2.0 * si.density)).sqrt(),
        1e-9,
        "Si [110] slow T",
    );
    close_rel(
        v110[1],
        (si.c[3][3] / si.density).sqrt(),
        1e-9,
        "Si [110] fast T",
    );
}

/// Silicon Young's moduli (Hopcroft, Nix & Kenny 2010): E[100] ~130 GPa,
/// E[110] ~169 GPa, E[111] ~188 GPa.
#[test]
fn silicon_youngs_moduli_match_hopcroft() {
    let si = material("Si").expect("Si in table");
    close_rel(
        youngs_modulus(&si.c, [1.0, 0.0, 0.0]) / 1e9,
        130.0,
        1e-2,
        "Si E[100]",
    );
    close_rel(
        youngs_modulus(&si.c, [1.0, 1.0, 0.0]) / 1e9,
        169.0,
        1e-2,
        "Si E[110]",
    );
    close_rel(
        youngs_modulus(&si.c, [1.0, 1.0, 1.0]) / 1e9,
        188.0,
        1e-2,
        "Si E[111]",
    );
}

/// Tungsten is nearly isotropic (A ~ 1.0); copper is strongly anisotropic
/// (A ~ 3.2). Zener ratio A = 2C44/(C11 - C12).
#[test]
fn zener_ratios_tungsten_and_copper() {
    let w = material("W").expect("W in table");
    let cu = material("Cu").expect("Cu in table");
    let aw = zener_ratio(w.c[0][0], w.c[0][1], w.c[3][3]);
    let acu = zener_ratio(cu.c[0][0], cu.c[0][1], cu.c[3][3]);
    println!("Zener A: W={aw:.4}, Cu={acu:.4}");
    close_rel(aw, 1.01, 3e-2, "W Zener");
    close_rel(acu, 3.21, 3e-2, "Cu Zener");
}

/// For a cubic crystal the [100] family and the [110] family must give the
/// same directional modulus under symmetry-equivalent directions.
#[test]
fn cubic_directions_are_symmetry_equivalent() {
    let si = material("Si").expect("Si in table");
    let e = |n| youngs_modulus(&si.c, n);
    for (a, b) in [
        ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        ([1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        ([1.0, 1.0, 0.0], [1.0, 0.0, 1.0]),
        ([1.0, 1.0, 1.0], [1.0, -1.0, 1.0]),
    ] {
        close_rel(e(a), e(b), 1e-9, "symmetry-equivalent E");
    }
}

/// Rotating the crystal must leave the scalar directional properties
/// unchanged when the direction is rotated with it, and the rotated
/// stiffness must remain symmetric and positive definite-sounding (real
/// Christoffel speeds).
#[test]
fn rotation_preserves_christoffel_speeds() {
    let cu = material("Cu").expect("Cu in table");
    let r = rotation_axis_angle([1.0, 2.0, 3.0], 0.9);
    let cr = rotate_stiffness(&cu.c, &r);
    for n in [
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [1.0, 1.0, 1.0],
        [0.3, 1.0, -0.5],
    ] {
        // The rotated frame direction corresponding to crystal `n` is R n.
        let rn = [
            r[0][0] * n[0] + r[0][1] * n[1] + r[0][2] * n[2],
            r[1][0] * n[0] + r[1][1] * n[1] + r[1][2] * n[2],
            r[2][0] * n[0] + r[2][1] * n[1] + r[2][2] * n[2],
        ];
        let (v0, _) = christoffel_modes(&cu.c, cu.density, n);
        let (v1, _) = christoffel_modes(&cr, cu.density, rn);
        for m in 0..3 {
            close_rel(v1[m], v0[m], 1e-9, "rotated Christoffel speed");
        }
        close_rel(
            youngs_modulus(&cr, rn),
            youngs_modulus(&cu.c, n),
            1e-9,
            "rotated E",
        );
    }
}

/// The Christoffel eigenvectors are orthonormal and satisfy
/// `Γ p = v² p`.
#[test]
fn christoffel_polarisations_are_eigenvectors() {
    for name in ["Si", "Cu", "W", "NaCl"] {
        let m = material(name).expect("material in table");
        for n in [
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
            [0.2, -0.7, 0.4],
        ] {
            let (v, p) = christoffel_modes(&m.c, m.density, n);
            for a in 0..3 {
                for b in 0..3 {
                    let dot: f64 = (0..3).map(|k| p[a][k] * p[b][k]).sum();
                    let expect = if a == b { 1.0 } else { 0.0 };
                    assert!(
                        (dot - expect).abs() < 1e-9,
                        "{name}: polarisations not orthonormal"
                    );
                }
                // Γ p should equal v² p; check the residual norm.
                let g = materials::anisotropic::christoffel(&m.c, m.density, n);
                for i in 0..3 {
                    let gp: f64 = (0..3).map(|k| g[i][k] * p[a][k]).sum();
                    let v2p = v[a] * v[a] * p[a][i];
                    assert!(
                        (gp - v2p).abs() < 1e-6 * v[a].powi(2).max(1.0),
                        "{name}: polarisation is not a Christoffel eigenvector"
                    );
                }
            }
        }
    }
}
