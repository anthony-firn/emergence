//! Benchmarks for the FFT demagnetising field (Newell-Williams-Dunlop
//! tensor, zero-padded convolution) and for the fold between the
//! local-shape closure and the full non-local field.

use magnetism::demag::DemagFft;
use magnetism::geometry::dot;
use magnetism::llg3d::Llg3d;
use magnetism::micromag::MU0;

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

fn idx(_nx: usize, ny: usize, nz: usize, i: usize, j: usize, k: usize) -> usize {
    (i * ny + j) * nz + k
}

/// Local-shape magnetostatic energy for diagonal shape factors `n`.
fn local_energy(m: &[[f64; 3]], cell: [f64; 3], n: [f64; 3], ms: f64) -> f64 {
    let v = cell[0] * cell[1] * cell[2];
    let mut acc = 0.0;
    for mi in m {
        acc += n[0] * mi[0] * mi[0] + n[1] * mi[1] * mi[1] + n[2] * mi[2] * mi[2];
    }
    0.5 * MU0 * ms * ms * v * acc
}

#[test]
fn cube_demag_factor_is_one_third() {
    // A uniformly magnetised cube has H = -M/3 on every interior cell; for
    // an odd grid the centre cell sits at the symmetry point and the
    // transverse components vanish exactly.
    let n = 5;
    let d = DemagFft::new(n, n, n, [1.0, 1.0, 1.0]);
    let m = vec![[1.0, 0.0, 0.0]; n * n * n];
    let mut h = vec![[0.0; 3]; n * n * n];
    d.field(&m, &mut h);
    let c = idx(n, n, n, n / 2, n / 2, n / 2);
    assert!(rel(h[c][0], -1.0 / 3.0) < 1e-12, "Hx = {}", h[c][0]);
    assert!(h[c][1].abs() < 1e-12 && h[c][2].abs() < 1e-12);
    // Field is linear in M.
    let m2 = vec![[2.0, 0.0, 0.0]; n * n * n];
    let mut h2 = vec![[0.0; 3]; n * n * n];
    d.field(&m2, &mut h2);
    assert!(rel(h2[c][0], 2.0 * h[c][0]) < 1e-12);
}

#[test]
fn sphere_cells_converge_to_one_third() {
    // A sphere of cubic cells on a cubic lattice: the Lorentz-sphere sum is
    // 4 pi / 3 for the field (i.e. H = -M/3) and is independent of the
    // radius once the cells form complete spherical shells.
    let n = 21;
    let c = (n / 2) as i64;
    for r in [3.0_f64, 5.0, 7.0] {
        let d = DemagFft::new(n, n, n, [1.0, 1.0, 1.0]);
        let mut m = vec![[0.0; 3]; n * n * n];
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let (di, dj, dk) = (i as i64 - c, j as i64 - c, k as i64 - c);
                    if (di * di + dj * dj + dk * dk) as f64 <= r * r {
                        m[idx(n, n, n, i, j, k)] = [1.0, 0.0, 0.0];
                    }
                }
            }
        }
        let mut h = vec![[0.0; 3]; n * n * n];
        d.field(&m, &mut h);
        let hc = h[idx(n, n, n, n / 2, n / 2, n / 2)][0];
        assert!(
            rel(hc, -1.0 / 3.0) < 1e-12,
            "sphere r={r}: centre Hx = {hc}, expected -1/3"
        );
    }
}

#[test]
fn thin_film_nz_approaches_one() {
    // A uniformly out-of-plane magnetised thin film: N_z -> 1 as the lateral
    // size grows and the thickness shrinks.
    let mut vals = Vec::new();
    for nxy in [8usize, 16, 32] {
        let nz = 1;
        let d = DemagFft::new(nxy, nxy, nz, [1.0, 1.0, 0.1]);
        let m = vec![[0.0, 0.0, 1.0]; nxy * nxy * nz];
        let mut h = vec![[0.0; 3]; nxy * nxy * nz];
        d.field(&m, &mut h);
        let c = idx(nxy, nxy, nz, nxy / 2, nxy / 2, 0);
        vals.push(h[c][2]);
    }
    for v in &vals {
        assert!(*v < 0.0 && *v > -1.0 - 1e-12, "Hz = {v}");
    }
    assert!(vals[0] > vals[1] && vals[1] > vals[2], "not monotone: {vals:?}");
    assert!(vals[2] < -0.97, "largest film Hz = {}", vals[2]);
}

#[test]
fn demag_energy_has_cubic_symmetry() {
    let n = 5;
    let d = DemagFft::new(n, n, n, [1.0, 1.0, 1.0]);
    let ms = 8.0e5;
    let ex = |dir: [f64; 3]| d.energy(&vec![dir; n * n * n], ms, MU0);
    let exx = ex([1.0, 0.0, 0.0]);
    let eyy = ex([0.0, 1.0, 0.0]);
    let ezz = ex([0.0, 0.0, 1.0]);
    assert!((exx - eyy).abs() / exx.abs() < 1e-12);
    assert!((exx - ezz).abs() / exx.abs() < 1e-12);
    // Energy is even in M and quadratic.
    assert!((d.energy(&vec![[-1.0, 0.0, 0.0]; n * n * n], ms, MU0) - exx).abs() < 1e-12);
    let e2 = d.energy(&vec![[2.0, 0.0, 0.0]; n * n * n], ms, MU0);
    assert!((e2 - 4.0 * exx).abs() / e2.abs() < 1e-12);
    // Uniform cube: E -> (mu0/2) M_s^2 V / 3.
    let v = (n * n * n) as f64; // unit cells
    let macroscopic = 0.5 * MU0 * ms * ms * v / 3.0;
    assert!(rel(exx, macroscopic) < 0.01, "energy {exx} vs {macroscopic}");
}

#[test]
fn fold_local_shape_is_enough_for_uniform_thin_film() {
    // Uniform out-of-plane thin film.  The interior FFT field converges to
    // the local shape factor N_z = 1, and the full and local energies agree
    // to within the finite-lateral-size edge correction.  The indicator is
    // the cell aspect ratio t/L: here t/L << 1.
    let (nx, ny, nz) = (32usize, 32, 1);
    let cell = [1.0, 1.0, 0.05];
    let d = DemagFft::new(nx, ny, nz, cell);
    let m = vec![[0.0, 0.0, 1.0]; nx * ny * nz];
    let mut h = vec![[0.0; 3]; nx * ny * nz];
    d.field(&m, &mut h);
    let c = idx(nx, ny, nz, nx / 2, ny / 2, 0);
    let hz = h[c][2];
    assert!(rel(hz, -1.0) < 0.02, "film interior Hz {hz} vs -1");
    let full = d.energy(&m, 8.0e5, MU0);
    let local = local_energy(&m, cell, [0.0, 0.0, 1.0], 8.0e5);
    assert!(rel(full, local) < 0.05, "film energy {full} vs local {local}");
}

#[test]
fn fold_local_shape_is_enough_for_uniform_needle() {
    // Uniformly magnetised elongated particle (a rod along x).  The local
    // shape approximation with the geometric factor captures the shape
    // anisotropy: the axial demag factor is small and the transverse one is
    // large.  The internal field is only approximately uniform along the
    // rod (a finite-length correction), which is why the local closure is
    // adequate but not exact.  The indicator is the aspect ratio, not the
    // exchange length.
    let (nx, ny, nz) = (24usize, 2, 2);
    let d = DemagFft::new(nx, ny, nz, [1.0, 1.0, 1.0]);
    let m = vec![[1.0, 0.0, 0.0]; nx * ny * nz];
    let mut h = vec![[0.0; 3]; nx * ny * nz];
    d.field(&m, &mut h);
    let mut avg = 0.0;
    for hc in &h {
        avg += hc[0];
    }
    avg /= h.len() as f64;
    // Transverse response of the same rod.
    let mt = vec![[0.0, 1.0, 0.0]; nx * ny * nz];
    let mut ht = vec![[0.0; 3]; nx * ny * nz];
    d.field(&mt, &mut ht);
    let mut avgt = 0.0;
    for hc in &ht {
        avgt += hc[1];
    }
    avgt /= ht.len() as f64;
    // Shape anisotropy: the axial factor is small and much smaller than the
    // transverse factor.
    assert!(avg > -0.1 && avg < 0.0, "needle N_x = {}", -avg);
    assert!(-avgt > 0.3, "needle N_y = {}", -avgt);
    assert!(-avgt + avg > 0.3, "shape anisotropy {}", -avgt + avg);
}

#[test]
fn fold_vortex_disc_needs_the_nonlocal_field() {
    // A vortex in a thin disc is the case where the local-shape closure
    // fails: the internal field is strongly non-uniform and the closed flux
    // lowers the energy below the uniform in-plane state in a way a single
    // shape factor cannot represent.
    let n = 25usize;
    let cell = [1.0, 1.0, 0.2];
    let d = DemagFft::new(n, n, 1, cell);
    let c = (n / 2) as f64;
    let r = 10.0;
    let mut vortex = vec![[0.0; 3]; n * n];
    for i in 0..n {
        for j in 0..n {
            let x = i as f64 - c;
            let y = j as f64 - c;
            if x * x + y * y <= r * r {
                let phi = y.atan2(x);
                vortex[idx(n, n, 1, i, j, 0)] = [-phi.sin(), phi.cos(), 0.0];
            }
        }
    }
    let mut uniform = vec![[0.0; 3]; n * n];
    for i in 0..n {
        for j in 0..n {
            let x = i as f64 - c;
            let y = j as f64 - c;
            if x * x + y * y <= r * r {
                uniform[idx(n, n, 1, i, j, 0)] = [1.0, 0.0, 0.0];
            }
        }
    }
    let ev = d.energy(&vortex, 8.0e5, MU0);
    let eu = d.energy(&uniform, 8.0e5, MU0);
    assert!(ev > 0.0 && eu > 0.0);
    assert!(ev < 0.6 * eu, "vortex {ev} vs uniform {eu}");
    // A single shape factor gives essentially the same (near-zero) energy
    // for both in-plane states, so it cannot resolve the gap.
    let lv = local_energy(&vortex, cell, [0.0, 0.0, 1.0], 8.0e5);
    let lu = local_energy(&uniform, cell, [0.0, 0.0, 1.0], 8.0e5);
    assert!(lv.abs() < 1e-6 && lu.abs() < 1e-6);
    // The internal field of the vortex is far from uniform.
    let mut hv = vec![[0.0; 3]; n * n];
    d.field(&vortex, &mut hv);
    let mut emean = [0.0; 3];
    let mut count = 0.0;
    for i in 0..n {
        for j in 0..n {
            let x = i as f64 - c;
            let y = j as f64 - c;
            if x * x + y * y <= r * r {
                let hc = hv[idx(n, n, 1, i, j, 0)];
                for a in 0..3 {
                    emean[a] += hc[a];
                }
                count += 1.0;
            }
        }
    }
    for a in 0..3 {
        emean[a] /= count;
    }
    let mut var = 0.0;
    for i in 0..n {
        for j in 0..n {
            let x = i as f64 - c;
            let y = j as f64 - c;
            if x * x + y * y <= r * r {
                let hc = hv[idx(n, n, 1, i, j, 0)];
                var += dot(
                    [hc[0] - emean[0], hc[1] - emean[1], hc[2] - emean[2]],
                    [hc[0] - emean[0], hc[1] - emean[1], hc[2] - emean[2]],
                );
            }
        }
    }
    var = (var / count).sqrt();
    // The non-local field is small (flux closure) but structured: the local
    // model predicts a near-zero in-plane field and misses it entirely.
    assert!(var > 0.005, "vortex field spread {var} (should be resolved)");
}

#[test]
fn llg3d_uniform_cube_interior_field_is_minus_m_over_three() {
    // The LLG effective field includes the physical demag field M_s * (-N m);
    // deep inside a uniformly magnetised cube this is -M_s/3.
    let n = 5;
    let ms = 8.0e5;
    let (nx, ny, nz) = (n, n, n);
    let mut s = Llg3d::new(nx, ny, nz, 1.0, 1.0, 1.0, 0.0, 0.0, ms, 0.5);
    s.demag = Some(DemagFft::new(nx, ny, nz, [1.0, 1.0, 1.0]));
    s.set_uniform([1.0, 0.0, 0.0]);
    let mut h = vec![[0.0; 3]; s.n_cells()];
    s.effective_field(&s.m, &mut h);
    let c = idx(nx, ny, nz, n / 2, n / 2, n / 2);
    assert!(rel(h[c][0], -ms / 3.0) < 1e-12, "Hx = {}", h[c][0]);
}
