//! The fast Poisson solver inverts the discrete box Laplacian exactly.

use immersed::poisson::{End, Poisson};

fn rng(seed: &mut u64) -> f64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64 / (1u64 << 53) as f64) - 0.5
}

#[test]
fn every_end_combination_is_inverted_to_round_off() {
    use End::*;
    let pairs = [[Neumann, Neumann], [Dirichlet, Dirichlet], [Neumann, Dirichlet], [Dirichlet, Neumann], [Periodic, Periodic]];
    let mut seed = 7;
    let mut worst: f64 = 0.0;
    for (sizes, h) in [([7usize, 6, 5], 0.1), ([12, 9, 8], 0.37), ([5, 4, 1], 0.2)] {
        for &ex in &pairs {
            for &ey in &pairs {
                for &ez in &[pairs[0], pairs[2], pairs[4]] {
                    if sizes[2] == 1 && ez != pairs[4] && ez != pairs[0] {
                        continue;
                    }
                    let p = Poisson::new(sizes, h, [ex, ey, ez]);
                    let n: usize = sizes.iter().product();
                    let mut phi: Vec<f64> = (0..n).map(|_| rng(&mut seed)).collect();
                    if p.singular() {
                        let m = phi.iter().sum::<f64>() / n as f64;
                        phi.iter_mut().for_each(|v| *v -= m);
                    }
                    let mut f = p.apply(&phi);
                    let incompat = p.solve(&mut f);
                    assert!(incompat.abs() < 1e-9, "compatible rhs reported incompatibility {incompat:e}");
                    let err = phi.iter().zip(&f).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
                    worst = worst.max(err);
                    assert!(err < 1e-10, "{sizes:?} {ex:?} {ey:?} {ez:?}: error {err:e}");
                }
            }
        }
    }
    println!("worst inversion error {worst:.2e}");
}
