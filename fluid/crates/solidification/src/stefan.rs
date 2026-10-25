//! The one-phase Stefan problem and its exact Neumann similarity solution.
//!
//! A pure melt, initially everywhere at its melting temperature `Tm`, is put
//! against a wall held at `Tw < Tm`. Solid grows from the wall. With constant
//! thermal diffusivity `alpha` in both phases and no density change the
//! interface position is
//!
//! ```text
//!   s(t) = 2 lambda sqrt(alpha t),
//!   lambda e^{lambda^2} erf(lambda) = St / sqrt(pi),
//! ```
//!
//! where `St = c (Tm - Tw) / L` is the Stefan number (latent heat per unit
//! volume `L`, volumetric heat capacity `c`). This is the classical Neumann
//! solution for the one-phase problem (Carslaw & Jaeger, *Conduction of Heat
//! in Solids*, 2nd ed., §11.2); the same form is quoted in
//! `analysis/compressible.md` §15. The temperature in the solid is
//!
//! ```text
//!   T(x,t) = Tw + (Tm - Tw) erf(x / (2 sqrt(alpha t))) / erf(lambda).
//! ```
//!
//! The functions here are the reference for the phase-field Stefan test in
//! [`crate::phase_field`] and for the 1D enthalpy solver [`crate::stefan::Stefan1d`].

/// `erf(x)` for `|x| < 0.5` via the Maclaurin series
/// `(2/sqrt(pi)) sum (-1)^n x^(2n+1)/(n!(2n+1))`. Beyond that the series loses
/// accuracy to cancellation and [`erfc`] switches to a continued fraction.
pub fn erf(x: f64) -> f64 {
    if x.abs() >= 1.5 {
        return if x >= 0.0 {
            1.0 - erfc_pos(x)
        } else {
            erfc_pos(-x) - 1.0
        };
    }
    let mut term = x;
    let mut sum = term;
    let mut n = 0u32;
    loop {
        n += 1;
        term *= -x * x * (2.0 * n as f64 - 1.0) / (n as f64 * (2.0 * n as f64 + 1.0));
        sum += term;
        if term.abs() < 1e-18 * sum.abs() || n > 400 {
            break;
        }
    }
    2.0 / std::f64::consts::PI.sqrt() * sum
}

/// `erfc(x) = 1 - erf(x)`.
pub fn erfc(x: f64) -> f64 {
    if x >= 1.5 {
        erfc_pos(x)
    } else if x >= 0.0 {
        1.0 - erf(x)
    } else {
        2.0 - erfc_pos(-x)
    }
}

/// `erfc(x)` for `x >= 0` by the modified-Lentz continued fraction
/// (Numerical Recipes §6.2), accurate to a few ulp after a dozen levels:
/// `erfc(x) = e^{-x^2}/sqrt(pi) / (x + (1/2)/(x + 1/(x + (3/2)/(x + ...))))`.
fn erfc_pos(x: f64) -> f64 {
    const TINY: f64 = 1e-300;
    let mut f = x;
    if f == 0.0 {
        f = TINY;
    }
    let mut c = f;
    let mut d = 0.0;
    let mut n = 0u32;
    loop {
        n += 1;
        let an = n as f64 * 0.5;
        d = x + an * d;
        if d == 0.0 {
            d = TINY;
        }
        c = x + an / c;
        if c == 0.0 {
            c = TINY;
        }
        d = 1.0 / d;
        let delta = c * d;
        f *= delta;
        if (delta - 1.0).abs() < 1e-16 || n > 300 {
            break;
        }
    }
    (-x * x).exp() / std::f64::consts::PI.sqrt() / f
}

/// Residual of the Neumann transcendental equation,
/// `f(lambda) = lambda e^{lambda^2} erf(lambda) - St/sqrt(pi)`.
pub fn neumann_residual(lambda: f64, stefan: f64) -> f64 {
    lambda * (lambda * lambda).exp() * erf(lambda) - stefan / std::f64::consts::PI.sqrt()
}

/// Solve for the similarity root `lambda` in
/// `lambda e^{lambda^2} erf(lambda) = St/sqrt(pi)`.
///
/// `f` is monotone increasing for `lambda >= 0`, so a bisection on
/// `[0, lambda_max]` with `lambda_max` grown until `f > 0` is robust. A
/// Newton polish then takes the residual to machine precision.
pub fn neumann_lambda(stefan: f64) -> f64 {
    assert!(stefan > 0.0, "Stefan number must be positive");
    let target = stefan / std::f64::consts::PI.sqrt();
    // f is monotone; grow the upper bracket.
    let mut hi = 0.5;
    let mut guard = 0;
    while neumann_residual(hi, stefan) < 0.0 {
        hi *= 1.5;
        guard += 1;
        assert!(guard < 200, "failed to bracket Neumann root");
    }
    let mut lo = 0.0;
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if neumann_residual(mid, stefan) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let mut lambda = 0.5 * (lo + hi);
    // Newton polish with numerical derivative.
    for _ in 0..8 {
        let f = neumann_residual(lambda, stefan);
        let h = 1e-7 * lambda.max(1e-3);
        let df = (neumann_residual(lambda + h, stefan) - neumann_residual(lambda - h, stefan))
            / (2.0 * h);
        if df.abs() < 1e-300 {
            break;
        }
        lambda -= f / df;
    }
    let _ = target;
    lambda
}

/// Interface position `s(t) = 2 lambda sqrt(alpha t)`.
pub fn neumann_front(lambda: f64, alpha: f64, t: f64) -> f64 {
    2.0 * lambda * (alpha * t).sqrt()
}

/// Temperature at depth `x` and time `t` in the solid, from the wall at `Tw`.
pub fn neumann_temperature(x: f64, t: f64, tw: f64, tm: f64, alpha: f64, lambda: f64) -> f64 {
    if t <= 0.0 {
        return tm;
    }
    let eta = x / (2.0 * (alpha * t).sqrt());
    tw + (tm - tw) * erf(eta) / erf(lambda)
}

/// Explicit 1D enthalpy solidification solver for a melt at rest at `Tm`
/// against a wall held at `Tw`. The enthalpy `h = c (T - Tm) + L f_s`
/// (per unit volume) is marched with face-centred conduction, then inverted
/// for `(T, f_s)` in each cell. Conductivity is the same in both phases so the
/// benchmark matches the constant-property Neumann solution.
///
/// The inversion is exact (linear interpolation between the two branches), so
/// there is no "mushy" smearing beyond the cell size.
pub struct Stefan1d {
    pub n: usize,
    pub dx: f64,
    pub alpha: f64,
    pub stefan: f64,
    pub tm: f64,
    pub tw: f64,
    /// temperature, K
    pub t: Vec<f64>,
    /// solid volume fraction in [0,1]
    pub fs: Vec<f64>,
    /// accumulated solid length (per unit area)
    solid_len: f64,
}

impl Stefan1d {
    /// `n` cells, unit length, diffusivity `alpha`, wall temperature `tw`,
    /// melting temperature `tm`. The latent heat enters only through `stefan`.
    pub fn new(n: usize, alpha: f64, stefan: f64, tm: f64, tw: f64) -> Self {
        assert!(tw < tm, "wall must be below the melting point");
        let t = vec![tm; n];
        let fs = vec![0.0; n];
        Self {
            n,
            dx: 1.0 / n as f64,
            alpha,
            stefan,
            tm,
            tw,
            t,
            fs,
            solid_len: 0.0,
        }
    }

    /// `L/c = (Tm - Tw) / St`, the enthalpy jump across the two-phase cell.
    fn inv_latent(&self) -> f64 {
        (self.tm - self.tw) / self.stefan
    }

    /// Stability time step `0.4 dx^2 / alpha`.
    pub fn stable_dt(&self) -> f64 {
        0.4 * self.dx * self.dx / self.alpha
    }

    /// Scaled enthalpy `y = (T - Tm) + (L/c) f_l`, in kelvin, where `f_l` is
    /// the liquid fraction. Liquid at `Tm` sits at `y = L/c`; removing heat
    /// walks `y` down through the two-phase interval `[0, L/c]` (solid at
    /// `Tm`), then into solid below `Tm`.
    fn enthalpy(&self, i: usize) -> f64 {
        (self.t[i] - self.tm) + self.inv_latent() * (1.0 - self.fs[i])
    }

    /// Invert the scaled enthalpy to `(T, f_s)`.
    fn invert(&self, y: f64) -> (f64, f64) {
        let l = self.inv_latent();
        if y <= 0.0 {
            (self.tm + y, 1.0)
        } else if y < l {
            (self.tm, 1.0 - y / l)
        } else {
            (self.tm + (y - l), 0.0)
        }
    }

    /// Advance one explicit step. The wall cell is pinned to `Tw`.
    pub fn step(&mut self, dt: f64) {
        let a = self.alpha;
        let n = self.n;
        let dx2 = self.dx * self.dx;
        let y: Vec<f64> = (0..n).map(|i| self.enthalpy(i)).collect();
        let mut new_t = self.t.clone();
        let mut new_fs = self.fs.clone();
        for i in 0..n {
            let t_left = if i == 0 { self.tw } else { self.t[i - 1] };
            let t_here = self.t[i];
            let t_right = if i + 1 < n { self.t[i + 1] } else { t_here };
            let lap = (t_left - 2.0 * t_here + t_right) / dx2;
            let (tt, ff) = self.invert(y[i] + a * lap * dt);
            new_t[i] = tt;
            new_fs[i] = ff;
        }
        // The wall is pinned: cell 0 is solid at Tw.
        new_t[0] = self.tw;
        new_fs[0] = 1.0;
        self.t = new_t;
        self.fs = new_fs;
        self.solid_len = self.fs.iter().sum::<f64>() * self.dx;
    }

    /// Solid length from the wall: `dx * sum(fs)`.
    pub fn front(&self) -> f64 {
        self.fs.iter().sum::<f64>() * self.dx
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn erf_matches_known_values() {
        assert!((erf(0.0)).abs() < 1e-15);
        assert!((erf(1.0) - 0.842_700_792_949_715).abs() < 1e-12);
        assert!((erf(0.5) - 0.520_499_877_813_047).abs() < 1e-12);
        assert!((erfc(2.0) - 0.004_677_734_981_047_266).abs() < 1e-12);
    }

    #[test]
    fn neumann_root_satisfies_equation() {
        for &st in &[0.05, 0.1, 0.5, 1.0, 2.0] {
            let l = neumann_lambda(st);
            assert!(neumann_residual(l, st).abs() < 1e-10, "St={st} residual");
        }
    }

    #[test]
    fn small_stefan_limit() {
        // For small St, lambda -> sqrt(St/2).
        let st = 1e-3;
        let l = neumann_lambda(st);
        assert!((l - (st / 2.0).sqrt()).abs() / l < 5e-3);
    }
}
