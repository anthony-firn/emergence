//! Static axisymmetric meniscus in a circular cup, with imposed wall height.
//!
//! Young–Laplace plus hydrostatic pressure, using far-field reference z=0:
//! h'' = h (1+h'²)^(3/2)/lc² - h'(1+h'²)/r.
//! h'(0)=0, h(R)=H. This is a prescribed-boundary static shape, not a
//! forecast of contact angle, bubbles, pouring, evaporation or coffee foam.
//! See de Gennes et al., Capillarity and Wetting Phenomena (2004), ch. 2;
//! the small-slope exact solution is H I0(r/lc)/I0(R/lc).

#[derive(Clone, Debug)]
pub struct Meniscus {
    pub radius: f64,
    pub wall_height: f64,
    pub capillary_length: f64,
    pub heights: Vec<f64>,
    pub slopes: Vec<f64>,
}

fn integrate(radius: f64, lc: f64, apex: f64, n: usize) -> Option<(Vec<f64>, Vec<f64>)> {
    let dr = radius / n as f64;
    let mut h = apex * (1. + dr * dr / (4. * lc * lc));
    let mut slope = apex * dr / (2. * lc * lc);
    let mut heights = vec![apex, h];
    let mut slopes = vec![0., slope];
    let rhs = |r: f64, z: f64, s: f64| {
        (
            s,
            z * (1. + s * s).powf(1.5) / (lc * lc) - s * (1. + s * s) / r,
        )
    };
    for i in 1..n {
        let r = i as f64 * dr;
        let a = rhs(r, h, slope);
        let b = rhs(r + dr / 2., h + dr * a.0 / 2., slope + dr * a.1 / 2.);
        let c = rhs(r + dr / 2., h + dr * b.0 / 2., slope + dr * b.1 / 2.);
        let d = rhs(r + dr, h + dr * c.0, slope + dr * c.1);
        h += dr * (a.0 + 2. * b.0 + 2. * c.0 + d.0) / 6.;
        slope += dr * (a.1 + 2. * b.1 + 2. * c.1 + d.1) / 6.;
        if !h.is_finite() || !slope.is_finite() || slope.abs() > 10. {
            return None;
        }
        heights.push(h);
        slopes.push(slope);
    }
    Some((heights, slopes))
}
impl Meniscus {
    pub fn solve(radius: f64, lc: f64, wall_height: f64, n: usize) -> Result<Self, &'static str> {
        if !radius.is_finite()
            || !lc.is_finite()
            || !wall_height.is_finite()
            || radius <= 0.
            || lc <= 0.
            || wall_height < 0.
            || n < 32
        {
            return Err("positive finite scales, nonnegative wall height, and at least 32 intervals required");
        }
        let (mut low, mut high) = (0., wall_height);
        for _ in 0..80 {
            let mid = (low + high) / 2.;
            let h = integrate(radius, lc, mid, n).map_or(f64::INFINITY, |p| *p.0.last().unwrap());
            if h > wall_height {
                high = mid;
            } else {
                low = mid;
            }
        }
        let (heights, slopes) =
            integrate(radius, lc, (low + high) / 2., n).ok_or("meniscus integration failed")?;
        if (heights[n] - wall_height).abs() > 1e-10 * lc.max(wall_height) {
            return Err("wall boundary did not converge");
        }
        Ok(Self {
            radius,
            wall_height,
            capillary_length: lc,
            heights,
            slopes,
        })
    }
    pub fn height(&self, r: f64) -> f64 {
        let n = self.heights.len() - 1;
        let x = (r / self.radius).clamp(0., 1.) * n as f64;
        let i = (x.floor() as usize).min(n - 1);
        let t = x - i as f64;
        // Cubic Hermite interpolation, using the integrated slope.
        let dr = self.radius / n as f64;
        (2. * t.powi(3) - 3. * t * t + 1.) * self.heights[i]
            + (t.powi(3) - 2. * t * t + t) * dr * self.slopes[i]
            + (-2. * t.powi(3) + 3. * t * t) * self.heights[i + 1]
            + (t.powi(3) - t * t) * dr * self.slopes[i + 1]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn i0(x: f64) -> f64 {
        let mut sum = 1.;
        let mut term = 1.;
        for k in 1..100 {
            term *= x * x / (4. * (k * k) as f64);
            sum += term;
            if term < sum * 1e-16 {
                break;
            }
        }
        sum
    }
    #[test]
    fn axisymmetric_meniscus_recovers_exact_linear_limit_and_refines() {
        let lc = crate::constants::water_capillary_length();
        let radius = 0.02;
        let height = 1e-6;
        let fine = Meniscus::solve(radius, lc, height, 1024).unwrap();
        let mut error = 0_f64;
        for i in 0..101 {
            let r = radius * i as f64 / 100.;
            let exact = height * i0(r / lc) / i0(radius / lc);
            error = error.max((fine.height(r) - exact).abs());
        }
        eprintln!("small-slope meniscus vs exact Bessel: max height error {error:e} m");
        assert!(error < 2e-13);
        let coarse = Meniscus::solve(0.044, lc, 0.0008, 256).unwrap();
        let fine = Meniscus::solve(0.044, lc, 0.0008, 512).unwrap();
        let reference = Meniscus::solve(0.044, lc, 0.0008, 2048).unwrap();
        let mut e_coarse = 0_f64;
        let mut e_fine = 0_f64;
        for i in 0..101 {
            let r = 0.044 * i as f64 / 100.;
            e_coarse = e_coarse.max((coarse.height(r) - reference.height(r)).abs());
            e_fine = e_fine.max((fine.height(r) - reference.height(r)).abs());
        }
        eprintln!("nonlinear meniscus refinement: {e_coarse:e} -> {e_fine:e} m");
        assert!(e_fine < e_coarse / 8.);
        assert!(e_fine < 1e-9);
    }
}
