//! A small deterministic pseudo-random generator (xorshift64*), so the
//! Monte-Carlo benchmarks are reproducible without external crates.

#[derive(Clone, Debug)]
pub struct Rng {
    state: u64,
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        // Any non-zero seed works for xorshift; mix it a little.
        Self {
            state: seed.wrapping_mul(0x9E3779B97F4A7C15) | 1,
        }
    }

    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }

    /// Uniform in `[0, 1)`.
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform on the unit sphere (Marsaglia's method).
    pub fn unit_vector(&mut self) -> [f64; 3] {
        loop {
            let x = 2.0 * self.uniform() - 1.0;
            let y = 2.0 * self.uniform() - 1.0;
            let s = x * x + y * y;
            if s < 1.0 {
                let f = 2.0 * (1.0 - s).sqrt();
                return [x * f, y * f, 1.0 - 2.0 * s];
            }
        }
    }
}
