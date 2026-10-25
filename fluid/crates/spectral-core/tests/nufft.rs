//! Non-uniform FFT against direct summation.

use num_complex::Complex64;
use spectral_core::grid::wavenumber;
use spectral_core::nufft::Nufft;
use spectral_core::Grid;

fn rng(seed: &mut u64) -> f64 {
    *seed ^= *seed << 13;
    *seed ^= *seed >> 7;
    *seed ^= *seed << 17;
    (*seed >> 11) as f64 / (1u64 << 53) as f64
}

/// Random real field on the retained modes (Hermitian on the kz = 0 plane
/// enforced by building it from a real physical field).
fn random_field(g: &Grid, seed: &mut u64) -> Vec<Complex64> {
    let fft = spectral_core::fft::Fft3::new(g);
    let x: Vec<f64> = (0..g.real_len()).map(|_| rng(seed) - 0.5).collect();
    let mut s = vec![Complex64::default(); g.spec_len()];
    fft.forward(&x, &mut s);
    for i in 0..g.spec_len() {
        if !g.kept(i) {
            s[i] = Complex64::default();
        }
    }
    s
}

fn direct_eval(g: &Grid, uh: &[Complex64], x: [f64; 3]) -> f64 {
    let mut acc = 0.0;
    for idx in 0..g.spec_len() {
        if !g.kept(idx) {
            continue;
        }
        let k = g.k(idx);
        let ph = Complex64::from_polar(1.0, k[0] * x[0] + k[1] * x[1] + k[2] * x[2]);
        acc += g.weight(idx) * (uh[idx] * ph).re;
    }
    acc
}

#[test]
fn interpolation_matches_direct_sum() {
    let g = Grid::new(24);
    let mut seed = 0x1234_5678_9abc_def1u64;
    let uh = random_field(&g, &mut seed);
    let pts: Vec<[f64; 3]> = (0..20).map(|_| [0, 1, 2].map(|_| rng(&mut seed) * 2.0 * std::f64::consts::PI)).collect();
    let scale = pts.iter().map(|&x| direct_eval(&g, &uh, x).abs()).fold(0.0, f64::max);
    for hw in [4, 6, 8, 10, 12] {
        let mut nf = Nufft::new(&g, hw);
        let vals = nf.interpolate(&uh, &pts);
        let err = pts.iter().zip(&vals).map(|(&x, v)| (direct_eval(&g, &uh, x) - v).abs()).fold(0.0, f64::max) / scale;
        eprintln!("interp half_width {hw}: rel err {err:.1e}");
        if hw == 8 {
            assert!(err < 1e-6);
        }
        if hw == 12 {
            assert!(err < 1e-9);
        }
    }
}

#[test]
fn spreading_matches_direct_sum_and_is_adjoint() {
    let g = Grid::new(24);
    let mut seed = 0x0fed_cba9_8765_4321u64;
    let pts: Vec<[f64; 3]> = (0..15).map(|_| [0, 1, 2].map(|_| rng(&mut seed) * 2.0 * std::f64::consts::PI)).collect();
    let f: Vec<f64> = (0..pts.len()).map(|_| rng(&mut seed) - 0.5).collect();
    let mut nf = Nufft::new(&g, 12);
    let mut fh = vec![Complex64::default(); g.spec_len()];
    nf.spread(&f, &pts, &mut fh);
    let v = (2.0 * std::f64::consts::PI).powi(3);
    let mut err: f64 = 0.0;
    let mut scale: f64 = 0.0;
    for idx in 0..g.spec_len() {
        if !g.kept(idx) {
            continue;
        }
        let k = g.k(idx);
        let direct: Complex64 = pts
            .iter()
            .zip(&f)
            .map(|(x, &fm)| fm * Complex64::from_polar(1.0, -(k[0] * x[0] + k[1] * x[1] + k[2] * x[2])))
            .sum::<Complex64>()
            / v;
        err = err.max((direct - fh[idx]).norm());
        scale = scale.max(direct.norm());
    }
    eprintln!("spread rel err {:.1e}", err / scale);
    assert!(err / scale < 1e-8);
    // adjointness: Σ_m F_m (J u)_m = (2π)³ ⟨S F, u⟩  (weighted full-spectrum sum)
    let uh = random_field(&g, &mut seed);
    let ju = nf.interpolate(&uh, &pts);
    let lhs: f64 = f.iter().zip(&ju).map(|(a, b)| a * b).sum();
    let rhs: f64 = v * (0..g.spec_len()).map(|i| g.weight(i) * (fh[i].conj() * uh[i]).re).sum::<f64>();
    eprintln!("adjoint: {lhs:.12e} vs {rhs:.12e}");
    assert!(((lhs - rhs) / lhs.abs()).abs() < 1e-8);
    let _ = wavenumber(0, 2);
}
