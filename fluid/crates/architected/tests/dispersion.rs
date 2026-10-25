//! Bloch-Floquet dispersion: long-wavelength slopes must equal the
//! homogenised stiffness.

use architected::bloch::*;
use architected::cell::*;
use architected::lattices::square_frame;

#[test]
fn acoustic_slopes_match_homogenised_stiffness() {
    let l = 1.0;
    let t = 0.1 * l;
    let rho = 1.0;
    let (raw, a1, a2) = square_frame(t, l);
    let cell = Cell2d::new(a1, a2, &raw, Kind::Frame, None);
    let c = cell.stiffness();
    let rho_star = rho * cell.relative_density();
    // Wave along x: shear (polarisation y) and longitudinal (polarisation x).
    let c_shear = c[2][2];
    let c_long = c[0][0];
    let expected = [ (c_shear / rho_star).sqrt(), (c_long / rho_star).sqrt() ];
    // Long-wavelength limit: acoustic branches start at omega = 0 at k = 0,
    // so the slopes are just omega(k)/k for small k. The real representation
    // doubles each branch, so take the distinct values.
    let dk = 1e-4;
    let w = bands_at(&cell, Wavevector { k1: dk, k2: 0.0 }, rho, 4);
    let mut got: Vec<f64> = Vec::new();
    for wi in &w {
        let s = wi / dk;
        if got.iter().all(|g| (g - s).abs() / s > 1e-3) {
            got.push(s);
        }
    }
    assert_eq!(got.len(), 2, "expected two distinct acoustic speeds: {w:?}");
    let mut exp_sorted = expected;
    exp_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    got.sort_by(|a, b| a.partial_cmp(b).unwrap());
    for (g, e) in got.iter().zip(exp_sorted.iter()) {
        assert!(
            (g - e).abs() / e < 2e-3,
            "slope {g:.6} vs sqrt(C/rho) {e:.6}"
        );
    }
}
