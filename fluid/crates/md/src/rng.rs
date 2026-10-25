//! Small deterministic pseudo-random number generator.
//!
//! The generator is xoshiro256++ seeded through SplitMix64 (Blackman &
//! Vigna 2018, <https://prng.di.unimi.it/>). It is *not* cryptographic;
//! it exists so simulations are bit-for-bit reproducible from a seed
//! without pulling a dependency into the workspace. Gaussian variates use
//! the Box–Muller transform, and Gamma variates the Marsaglia–Tsang
//! squeeze method (Marsaglia & Tsang 2000), which the Bussi thermostat
//! needs for its chi-square term.

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

    /// Uniform in `[lo, hi)`.
    #[inline]
    pub fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }

    /// Standard normal via Box–Muller.
    #[inline]
    pub fn normal(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-300);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    /// Gamma variate with given `shape` and unit scale (Marsaglia–Tsang).
    pub fn gamma(&mut self, shape: f64) -> f64 {
        debug_assert!(shape > 0.0);
        // Boost small shapes to >= 1 and correct by a power of a uniform.
        if shape < 1.0 {
            let u = self.next_f64().max(1e-300);
            return self.gamma(shape + 1.0) * u.powf(1.0 / shape);
        }
        let d = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let x = self.normal();
            let v = (1.0 + c * x).powi(3);
            if v <= 0.0 {
                continue;
            }
            let u = self.next_f64();
            let x2 = x * x;
            if u < 1.0 - 0.033_1 * x2 * x2 {
                return d * v;
            }
            if u.ln() < 0.5 * x2 + d * (1.0 - v + v.ln()) {
                return d * v;
            }
        }
    }

    /// Three Maxwell–Boltzmann velocity components with variance
    /// `temperature / mass` (reduced units, `k_B = 1`).
    pub fn maxwell_boltzmann(&mut self, temperature: f64, mass: f64) -> [f64; 3] {
        let sigma = (temperature / mass).sqrt();
        [
            sigma * self.normal(),
            sigma * self.normal(),
            sigma * self.normal(),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rng_is_deterministic_and_roughly_uniform() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        let mut mean = 0.0;
        let m = 200_000;
        for _ in 0..m {
            let x = a.next_f64();
            assert_eq!(x, b.next_f64());
            assert!((0.0..1.0).contains(&x));
            mean += x;
        }
        mean /= m as f64;
        assert!((mean - 0.5).abs() < 5e-3, "mean {mean}");

        // Gaussian moments.
        let mut r = Rng::new(7);
        let (mut m1, mut m2) = (0.0, 0.0);
        let k = 200_000;
        for _ in 0..k {
            let x = r.normal();
            m1 += x;
            m2 += x * x;
        }
        m1 /= k as f64;
        m2 /= k as f64;
        assert!(m1.abs() < 0.02, "normal mean {m1}");
        assert!((m2 - 1.0).abs() < 0.02, "normal var {m2}");

        // Gamma mean = shape, variance = shape.
        let shape = 300.0;
        let mut g = Rng::new(3);
        let (mut gm, mut gv) = (0.0, 0.0);
        let kg = 100_000;
        for _ in 0..kg {
            let x = g.gamma(shape);
            gm += x;
            gv += x * x;
        }
        gm /= kg as f64;
        gv = gv / kg as f64 - gm * gm;
        assert!((gm - shape).abs() < 3.0, "gamma mean {gm}");
        assert!((gv - shape).abs() < 20.0, "gamma var {gv}");
    }
}
