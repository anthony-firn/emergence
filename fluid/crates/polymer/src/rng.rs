//! Small deterministic pseudo-random number generator (xoshiro256++ seeded
//! through SplitMix64; Blackman & Vigna 2018, <https://prng.di.unimi.it/>).
//!
//! It is deliberately duplicated here (rather than depending on `md`) so
//! that the polymer crate stays a leaf: the Brownian-dynamics thermostat
//! only needs uniform and standard-normal variates.

/// Seedable xoshiro256++ generator.
#[derive(Clone)]
pub struct Rng {
    s: [u64; 4],
}

#[inline]
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        let mut sm = seed;
        let s = [
            splitmix64(&mut sm),
            splitmix64(&mut sm),
            splitmix64(&mut sm),
            splitmix64(&mut sm),
        ];
        Self { s }
    }

    /// Uniform in `[0, 2^64)`.
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let result = self.s[0]
            .wrapping_add(self.s[3])
            .rotate_left(23)
            .wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    /// Uniform in `[0, 1)` with 53 bits of resolution.
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    }

    /// Standard normal via Box–Muller.
    #[inline]
    pub fn normal(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-300);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_normal_moments() {
        let mut a = Rng::new(1);
        let mut b = Rng::new(1);
        let (mut m1, mut m2) = (0.0, 0.0);
        let n = 100_000;
        for _ in 0..n {
            let x = a.normal();
            assert_eq!(x, b.normal());
            m1 += x;
            m2 += x * x;
        }
        m1 /= n as f64;
        m2 = m2 / n as f64 - m1 * m1;
        assert!(m1.abs() < 0.02, "mean {m1}");
        assert!((m2 - 1.0).abs() < 0.03, "var {m2}");
    }
}
