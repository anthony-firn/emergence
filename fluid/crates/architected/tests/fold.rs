//! Fold/unfold: a finite honeycomb cantilever against the homogenised beam.
//!
//! The effective continuum is exact only as the sample grows. The tip
//! deflection relative to the Euler-Bernoulli prediction using `E* I` is
//! therefore the size-effect indicator; it converges to one with a boundary
//! layer at the clamp of about one cell.

use architected::cell::*;
use architected::fold::honeycomb_patch;
use architected::lattices::regular_honeycomb;
use std::f64::consts::PI;

fn continuum_ratio(nx: usize, ny: usize, t: f64, l: f64, e_star: f64) -> (f64, f64) {
    let patch = honeycomb_patch(nx, ny, t, l, 1.0, PI / 6.0);
    let (delta, length, height) = patch.cantilever_tip(1.0);
    let i_beam = height * height * height / 12.0;
    let continuum = length.powi(3) / (3.0 * e_star * i_beam);
    (delta.abs() / continuum, length / height)
}

#[test]
fn honeycomb_cantilever_converges_to_continuum() {
    let l = 1.0;
    let t = 0.02 * l;
    let (raw, a1, a2) = regular_honeycomb(t, l, Kind::Frame);
    let cell = Cell2d::new(a1, a2, &raw, Kind::Frame, None);
    let e_star = moduli2d(&cell.stiffness()).ex;

    // Fixed aspect ratio, growing cell count: L/H ~ 4.
    let samples: Vec<usize> = vec![4, 6, 8, 12];
    let mut ratios = Vec::new();
    for &n in &samples {
        let (ratio, aspect) = continuum_ratio(4 * n, n, t, l, e_star);
        assert!(
            (ratio - 1.0).abs() < 0.25,
            "N={n}, aspect {aspect:.2}: ratio {ratio:.4} far from 1"
        );
        ratios.push((n, ratio));
    }
    // Monotone approach to the continuum.
    for w in ratios.windows(2) {
        assert!(
            w[1].1 > w[0].1,
            "ratio not increasing: {:?}",
            ratios
        );
    }
    // The correction is a boundary layer: (1 - ratio) ~ c / N with c ~ 1.
    let n_last = ratios.last().unwrap().0 as f64;
    let c_last = (1.0 - ratios.last().unwrap().1) * n_last;
    assert!(
        c_last > 0.5 && c_last < 2.0,
        "boundary-layer coefficient {c_last:.3} cells, expected ~1"
    );
    assert!(
        ratios.last().unwrap().1 > 0.9,
        "largest sample ratio {:.4} should be close to 1",
        ratios.last().unwrap().1
    );
}

use architected::Kind;
