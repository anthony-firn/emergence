//! Johnson–Zollweg–Gubbins (1993) Lennard-Jones equation of state.
//!
//! The reference is a 33-parameter modified Benedict–Webb–Rubin fit to
//! computer-simulation data: J. K. Johnson, W. Zollweg and K. E. Gubbins,
//! "The Lennard-Jones equation of state revisited", *Mol. Phys.* **78**,
//! 591 (1993). The coefficients below are their Table 10 (the "full
//! potential" set, `gamma = 3`), and the pressure and residual internal
//! energy are their Eqs. (7) and (9). The fit represents the *untruncated*
//! 12-6 potential; simulations must add the homogeneous tail corrections
//! of [`crate::potential`] before comparing.
//!
//! Reduced units: `P* = P sigma^3/epsilon`, `U* = U/epsilon`, `T* = k_B T/epsilon`,
//! `rho* = rho sigma^3`.

/// The MBWR Lennard-Jones equation of state.
#[derive(Clone, Copy, Debug)]
pub struct JzgEos {
    /// Nonlinear parameter, fixed at 3 by the authors.
    pub gamma: f64,
}

impl Default for JzgEos {
    fn default() -> Self {
        Self::new()
    }
}

impl JzgEos {
    pub fn new() -> Self {
        Self { gamma: 3.0 }
    }

    #[rustfmt::skip]
    fn x(j: usize) -> f64 {
        const X: [f64; 33] = [
            0.0,
            0.862_308_509_750_742_1,
            2.976_218_765_822_098,
            -8.402_230_115_796_038,
            0.105_413_662_920_355_5,
            -0.856_458_382_817_459_8,
            1.582_759_470_107_601,
            0.763_942_194_830_545_3,
            1.753_173_414_312_048,
            2.798_291_772_190_376e3,
            -4.839_422_026_085_765_7e-2,
            0.996_326_519_772_193_5,
            -3.698_000_291_272_493e1,
            2.084_012_299_434_647e1,
            8.305_402_124_717_285e1,
            -9.574_799_715_203_068e2,
            -1.477_746_229_234_994e2,
            6.398_607_852_471_505e1,
            1.603_993_673_294_834e1,
            6.805_916_615_864_377e1,
            -2.791_293_578_795_945e3,
            -6.245_128_304_568_454,
            -8.116_836_104_958_410e3,
            1.488_735_559_561_229e1,
            -1.059_346_754_655_084e4,
            -1.131_607_632_802_822e2,
            -8.867_771_540_418_822e3,
            -3.986_982_844_450_543e1,
            -4.689_270_299_917_261e3,
            2.593_535_277_438_717e2,
            -2.694_523_589_434_903e3,
            -7.218_487_631_550_215e2,
            1.721_802_063_863_269e2,
        ];
        X[j]
    }

    /// Temperature-dependent `a_i`, Table 5.
    fn a(&self, t: f64) -> [f64; 9] {
        let x = |j| Self::x(j);
        let mut a = [0.0; 9];
        a[1] = x(1) * t + x(2) * t.sqrt() + x(3) + x(4) / t + x(5) / t.powi(2);
        a[2] = x(6) * t + x(7) + x(8) / t + x(9) / t.powi(2);
        a[3] = x(10) * t + x(11) + x(12) / t;
        a[4] = x(13);
        a[5] = x(14) / t + x(15) / t.powi(2);
        a[6] = x(16) / t;
        a[7] = x(17) / t + x(18) / t.powi(2);
        a[8] = x(19) / t.powi(2);
        a
    }

    /// Temperature-dependent `b_i`, Table 6.
    fn b(&self, t: f64) -> [f64; 7] {
        let x = |j| Self::x(j);
        let mut b = [0.0; 7];
        b[1] = x(20) / t.powi(2) + x(21) / t.powi(3);
        b[2] = x(22) / t.powi(2) + x(23) / t.powi(4);
        b[3] = x(24) / t.powi(2) + x(25) / t.powi(3);
        b[4] = x(26) / t.powi(2) + x(27) / t.powi(4);
        b[5] = x(28) / t.powi(2) + x(29) / t.powi(3);
        b[6] = x(30) / t.powi(2) + x(31) / t.powi(3) + x(32) / t.powi(4);
        b
    }

    /// Temperature-dependent `c_i`, Table 8 (internal-energy coefficients).
    fn c(&self, t: f64) -> [f64; 9] {
        let x = |j| Self::x(j);
        let mut c = [0.0; 9];
        c[1] = x(2) * t.sqrt() / 2.0 + x(3) + 2.0 * x(4) / t + 3.0 * x(5) / t.powi(2);
        c[2] = x(7) + 2.0 * x(8) / t + 3.0 * x(9) / t.powi(2);
        c[3] = x(11) + 2.0 * x(12) / t;
        c[4] = x(13);
        c[5] = 2.0 * x(14) / t + 3.0 * x(15) / t.powi(2);
        c[6] = 2.0 * x(16) / t;
        c[7] = 2.0 * x(17) / t + 3.0 * x(18) / t.powi(2);
        c[8] = 3.0 * x(19) / t.powi(2);
        c
    }

    /// Temperature-dependent `d_i`, Table 9 (internal-energy coefficients).
    fn d(&self, t: f64) -> [f64; 7] {
        let x = |j| Self::x(j);
        let mut d = [0.0; 7];
        d[1] = 3.0 * x(20) / t.powi(2) + 4.0 * x(21) / t.powi(3);
        d[2] = 3.0 * x(22) / t.powi(2) + 5.0 * x(23) / t.powi(4);
        d[3] = 3.0 * x(24) / t.powi(2) + 4.0 * x(25) / t.powi(3);
        d[4] = 3.0 * x(26) / t.powi(2) + 5.0 * x(27) / t.powi(4);
        d[5] = 3.0 * x(28) / t.powi(2) + 4.0 * x(29) / t.powi(3);
        d[6] = 3.0 * x(30) / t.powi(2) + 4.0 * x(31) / t.powi(3) + 5.0 * x(32) / t.powi(4);
        d
    }

    fn g(&self, rho: f64) -> [f64; 7] {
        let f = (-self.gamma * rho * rho).exp();
        let mut g = [0.0; 7];
        g[1] = (1.0 - f) / (2.0 * self.gamma);
        for i in 2..=6 {
            let p = rho.powi(2 * (i as i32 - 1));
            g[i] = -(f * p - 2.0 * (i as f64 - 1.0) * g[i - 1]) / (2.0 * self.gamma);
        }
        g
    }

    /// Pressure at reduced density `rho` and temperature `t`, Eq. (7).
    pub fn pressure(&self, rho: f64, t: f64) -> f64 {
        let a = self.a(t);
        let b = self.b(t);
        let f = (-self.gamma * rho * rho).exp();
        let mut p = rho * t;
        for i in 1..=8 {
            p += a[i] * rho.powi(i as i32 + 1);
        }
        for i in 1..=6 {
            p += f * b[i] * rho.powi(2 * i as i32 + 1);
        }
        p
    }

    /// Residual (configurational) internal energy per atom, Eq. (9). Add
    /// `(3/2) T` for the total energy.
    pub fn residual_energy(&self, rho: f64, t: f64) -> f64 {
        let c = self.c(t);
        let d = self.d(t);
        let g = self.g(rho);
        let mut u = 0.0;
        for i in 1..=8 {
            u += c[i] * rho.powi(i as i32) / i as f64;
        }
        for i in 1..=6 {
            u += d[i] * g[i];
        }
        u
    }

    /// Total internal energy per atom `3T/2 + U_res`.
    pub fn total_energy(&self, rho: f64, t: f64) -> f64 {
        1.5 * t + self.residual_energy(rho, t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reproduces_published_eos_values() {
        // Values evaluated from the JZG coefficient tables; these are the
        // numbers the MD tests compare against.
        let eos = JzgEos::new();
        let p = eos.pressure(0.5, 2.0);
        let u = eos.residual_energy(0.5, 2.0);
        assert!((p - 1.077_450).abs() < 1e-5, "P* {p}");
        assert!((u - (-3.144_943)).abs() < 1e-5, "U_res* {u}");

        let p2 = eos.pressure(0.8, 1.0);
        let u2 = eos.residual_energy(0.8, 1.0);
        assert!((p2 - 1.031_806).abs() < 1e-5, "P* {p2}");
        assert!((u2 - (-5.523_454)).abs() < 1e-5, "U_res* {u2}");
    }

    #[test]
    fn low_density_limit_is_ideal_gas() {
        let eos = JzgEos::new();
        let rho = 1e-6;
        let t = 2.0;
        let p = eos.pressure(rho, t);
        assert!((p - rho * t).abs() / (rho * t) < 1e-5);
    }
}
