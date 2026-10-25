//! Potential-field source-surface extrapolation of surface magnetograms.
//! Analytic harmonic gradients, no currents/free-energy storage/reconnection.
//! This supplies a coronal baseline, NOT a flare/CME trigger or forecast.
use crate::surface::Surface;
use gravity::{norm, Vec3};
use num_complex::Complex64;
use std::f64::consts::{PI, TAU};
pub const MU0: f64 = 1.256_637_062_12e-6;
fn index(l: usize, m: usize) -> usize {
    l * (l + 1) / 2 + m
}
/// Gauss-Legendre nodes/weights, generated once per harmonic fit. Integrating
/// in theta makes the associated-Legendre factors regular at both poles.
fn quadrature(n: usize) -> Vec<(f64, f64)> {
    (0..n)
        .map(|i| {
            let mut x = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
            let derivative = |x: f64| {
                let (mut previous, mut current) = (1.0, x);
                for k in 2..=n {
                    let next =
                        ((2 * k - 1) as f64 * x * current - (k - 1) as f64 * previous) / k as f64;
                    previous = current;
                    current = next;
                }
                (current, n as f64 * (x * current - previous) / (x * x - 1.0))
            };
            for _ in 0..32 {
                let (value, gradient) = derivative(x);
                let change = value / gradient;
                x -= change;
                if change.abs() < 2e-16 {
                    break;
                }
            }
            let (_, gradient) = derivative(x);
            (x, 2.0 / ((1.0 - x * x) * gradient * gradient))
        })
        .collect()
}
/// Orthonormal associated Legendre Q_lm with Condon-Shortley phase and its
/// theta derivative. Recurrences avoid factorial overflow and polar division.
fn angular(max: usize, mu: f64, sin: f64) -> (Vec<f64>, Vec<f64>) {
    let n = index(max, max) + 1;
    let mut q = vec![0.0; n];
    let mut d = vec![0.0; n];
    q[0] = 1.0 / (4.0 * PI).sqrt();
    for m in 0..=max {
        if m > 0 {
            let factor = -((2 * m + 1) as f64 / (2 * m) as f64).sqrt();
            let prev = index(m - 1, m - 1);
            let here = index(m, m);
            q[here] = factor * sin * q[prev];
            d[here] = factor * (mu * q[prev] + sin * d[prev]);
        }
        if m < max {
            let a = ((2 * m + 3) as f64).sqrt();
            let here = index(m + 1, m);
            let prev = index(m, m);
            q[here] = a * mu * q[prev];
            d[here] = a * (-sin * q[prev] + mu * d[prev]);
        }
        for l in m + 2..=max {
            let lf = l as f64;
            let mf = m as f64;
            let a = ((4.0 * lf * lf - 1.0) / (lf * lf - mf * mf)).sqrt();
            let b = ((2.0 * lf + 1.0) * ((lf - 1.0).powi(2) - mf * mf)
                / ((2.0 * lf - 3.0) * (lf * lf - mf * mf)))
                .sqrt();
            let here = index(l, m);
            let prev = index(l - 1, m);
            let old = index(l - 2, m);
            q[here] = a * mu * q[prev] - b * q[old];
            d[here] = a * (-sin * q[prev] + mu * d[prev]) - b * d[old];
        }
    }
    (q, d)
}
pub struct Pfss {
    radius: f64,
    source_ratio: f64,
    max: usize,
    coefficients: Vec<Complex64>,
}
impl Pfss {
    pub fn photospheric_radius(&self) -> f64 {
        self.radius
    }
    pub fn source_radius(&self) -> f64 {
        self.radius * self.source_ratio
    }
    /// Project the piecewise constant equal-area cell averages, integrating
    /// each cell with polar-regular Gauss quadrature and an analytic longitude
    /// sinc factor. This avoids interpreting averaged observations as points.
    /// Input must be explicitly flux balanced. Harmonic truncation and the
    /// piecewise constant boundary representation still require refinement.
    pub fn from_surface(
        surface: &Surface,
        max: usize,
        source_ratio: f64,
    ) -> Result<Self, &'static str> {
        let (n_mu, n_phi) = surface.dimensions();
        if max == 0
            || max > 64
            || max > n_mu / 2
            || !source_ratio.is_finite()
            || source_ratio <= 1.0
            || source_ratio > 100.0
        {
            return Err("invalid coronal harmonic/source-surface resolution");
        }
        let scale = surface.field().iter().map(|b| b.abs()).fold(0.0, f64::max);
        let mean = surface.field().iter().sum::<f64>() / surface.field().len() as f64;
        if mean.abs() > 1e-8 * scale {
            return Err("magnetogram has net flux; explicitly balance observations before PFSS");
        }
        let max_m = if n_phi == 1 { 0 } else { max };
        if n_phi > 1 && n_phi < 2 * max + 1 {
            return Err("longitude grid does not resolve coronal harmonics");
        }
        let mut coefficients = vec![Complex64::default(); index(max, max) + 1];
        let delta_phi = TAU / n_phi as f64;
        let nodes = quadrature(((max + 2) / 2).max(8));
        for i in 0..n_mu {
            let lo = (-1.0 + i as f64 * 2.0 / n_mu as f64).acos();
            let hi = (-1.0 + (i + 1) as f64 * 2.0 / n_mu as f64)
                .clamp(-1.0, 1.0)
                .acos();
            let midpoint = 0.5 * (lo + hi);
            let half = 0.5 * (lo - hi);
            let mut latitude = vec![0.0; coefficients.len()];
            for (node, weight) in &nodes {
                let theta = midpoint + half * node;
                let (sin, mu) = theta.sin_cos();
                let (q, _) = angular(max, mu, sin);
                for (sum, value) in latitude.iter_mut().zip(q) {
                    *sum += half * weight * sin * value;
                }
            }
            for m in 0..=max_m {
                let mut azimuth = Complex64::default();
                for j in 0..n_phi {
                    azimuth += Complex64::from_polar(
                        surface.field()[i * n_phi + j],
                        -(m as f64) * surface.longitude(j),
                    );
                }
                let a = m as f64 * delta_phi * 0.5;
                let longitude = delta_phi * if m == 0 { 1.0 } else { a.sin() / a };
                for l in m.max(1)..=max {
                    coefficients[index(l, m)] += azimuth * (longitude * latitude[index(l, m)]);
                }
            }
        }
        Ok(Self {
            radius: surface.parameters().radius,
            source_ratio,
            max,
            coefficients,
        })
    }
    fn radial(&self, l: usize, r: f64) -> (f64, f64) {
        let inv_source = self.source_ratio.powi(-((2 * l + 1) as i32));
        let den = l as f64 * inv_source + (l + 1) as f64;
        let grow = r.powi(l as i32 - 1) * inv_source;
        let fall = r.powi(-(l as i32) - 2);
        (
            (l as f64 * grow + (l + 1) as f64 * fall) / den,
            (grow - fall) / den,
        )
    }
    /// Cartesian magnetic field [T], radius measured from the solar center.
    /// Exact axis limits are used at either pole, including nonaxisymmetric m=1.
    pub fn magnetic_field(&self, position: Vec3) -> Result<Vec3, &'static str> {
        if position.0.iter().any(|x| !x.is_finite()) {
            return Err("invalid coronal query");
        }
        let radius = norm(position);
        let r = radius / self.radius;
        if r < 1.0 - 8.0 * f64::EPSILON || r > self.source_ratio * (1.0 + 8.0 * f64::EPSILON) {
            return Err("coronal radius outside photosphere/source surface");
        }
        let r = r.clamp(1.0, self.source_ratio);
        let mu = (position.0[2] / radius).clamp(-1.0, 1.0);
        let horizontal = position.0[0].hypot(position.0[1]);
        let sin = horizontal / radius;
        if sin < 1e-14 {
            let mut out = Vec3::default();
            for l in 1..=self.max {
                let (fr, ft) = self.radial(l, r);
                let sign = if mu < 0.0 && l % 2 == 1 { -1.0 } else { 1.0 };
                let q = ((2 * l + 1) as f64 / (4.0 * PI)).sqrt() * sign;
                out.0[2] += self.coefficients[index(l, 0)].re * fr * q * mu;
                let mut k = -0.5 * ((l * (l + 1) * (2 * l + 1)) as f64 / (4.0 * PI)).sqrt();
                if mu < 0.0 && l % 2 == 0 {
                    k = -k;
                }
                let c = self.coefficients[index(l, 1)];
                out.0[0] += 2.0 * ft * k * c.re;
                out.0[1] -= 2.0 * ft * k * c.im;
            }
            return checked(out);
        }
        // Reuse azimuthal powers across degrees instead of one trig pair per
        // (l,m). This is useful for the many queries in field-line tracing.
        let phase = Complex64::new(position.0[0] / horizontal, position.0[1] / horizontal);
        let mut phases = vec![Complex64::new(1.0, 0.0); self.max + 1];
        for m in 1..=self.max {
            phases[m] = phases[m - 1] * phase;
        }
        let (q, d) = angular(self.max, mu, sin);
        let (mut br, mut bt, mut bp) = (0.0, 0.0, 0.0);
        for l in 1..=self.max {
            let (fr, ft) = self.radial(l, r);
            for (m, phase) in phases.iter().enumerate().take(l + 1) {
                let k = index(l, m);
                let factor = if m == 0 { 1.0 } else { 2.0 };
                let c = self.coefficients[k] * phase;
                br += factor * fr * c.re * q[k];
                bt += factor * ft * c.re * d[k];
                bp -= factor * ft * c.im * m as f64 * q[k] / sin;
            }
        }
        let (cp, sp) = (position.0[0] / horizontal, position.0[1] / horizontal);
        checked(Vec3([
            br * sin * cp + bt * mu * cp - bp * sp,
            br * sin * sp + bt * mu * sp + bp * cp,
            br * mu - bt * sin,
        ]))
    }
    /// Unsigned magnetic flux through the radial source surface [Wb].
    pub fn open_flux(&self, n_mu: usize, n_phi: usize) -> Result<f64, &'static str> {
        if n_mu < 4 || n_phi < 4 {
            return Err("invalid open-flux quadrature");
        }
        let r = self.radius * self.source_ratio;
        let mut sum = 0.0;
        for i in 0..n_mu {
            let mu = -1.0 + (i as f64 + 0.5) * 2.0 / n_mu as f64;
            let sin = (1.0 - mu * mu).sqrt();
            for j in 0..n_phi {
                let phi = (j as f64 + 0.5) * TAU / n_phi as f64;
                let unit = Vec3([sin * phi.cos(), sin * phi.sin(), mu]);
                sum += self.magnetic_field(unit * r)?.dot(unit).abs();
            }
        }
        Ok(sum * (4.0 * PI * r * r / (n_mu * n_phi) as f64))
    }
}
fn checked(v: Vec3) -> Result<Vec3, &'static str> {
    if v.0.iter().any(|x| !x.is_finite()) {
        Err("unrepresentable coronal field")
    } else {
        Ok(v)
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Energy {
    pub magnetic: f64,
    pub potential_reference: f64,
    pub excess: f64,
}
/// Diagnostic only. Free magnetic energy requires a consistent reference with
/// the same boundary flux and volume. Negative excess is reported, not clamped
/// into a fabricated flare reservoir. Input cell volumes are SI m³.
pub fn energy(field: &[Vec3], reference: &[Vec3], volumes: &[f64]) -> Result<Energy, &'static str> {
    if field.len() != reference.len() || field.len() != volumes.len() {
        return Err("inconsistent coronal energy grid");
    }
    let (mut magnetic, mut potential_reference) = (0.0, 0.0);
    for ((b, p), v) in field.iter().zip(reference).zip(volumes) {
        checked(*b)?;
        checked(*p)?;
        if !v.is_finite() || *v <= 0.0 {
            return Err("invalid coronal volume");
        }
        magnetic += b.dot(*b) * v / (2.0 * MU0);
        potential_reference += p.dot(*p) * v / (2.0 * MU0);
    }
    if !magnetic.is_finite() || !potential_reference.is_finite() {
        return Err("unrepresentable coronal energy");
    }
    Ok(Energy {
        magnetic,
        potential_reference,
        excess: magnetic - potential_reference,
    })
}
