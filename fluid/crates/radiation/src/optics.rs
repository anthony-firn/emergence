use crate::math::{around, Rng, Vec3};
use crate::{C, H, KB};
use std::f64::consts::PI;

/// Tabulated spectrum, linearly interpolated; endpoints clamped.
/// Caller must keep wavelengths within the measured range for accuracy.
#[derive(Clone, Debug)]
pub struct Spectrum(pub Vec<(f64, f64)>);
impl Spectrum {
    pub fn new(points: Vec<(f64, f64)>) -> Self {
        assert!(!points.is_empty());
        assert!(points
            .iter()
            .all(|&(l, v)| l.is_finite() && l > 0.0 && v.is_finite() && v >= 0.0));
        assert!(points.windows(2).all(|p| p[0].0 < p[1].0));
        Self(points)
    }
    pub fn constant(v: f64) -> Self {
        Self::new(vec![(550e-9, v)])
    }
    pub fn at(&self, l: f64) -> f64 {
        let i = self.0.partition_point(|p| p.0 < l);
        if i == 0 {
            return self.0[0].1;
        }
        if i == self.0.len() {
            return self.0[i - 1].1;
        }
        let (a, b) = (self.0[i - 1], self.0[i]);
        a.1 + (b.1 - a.1) * (l - a.0) / (b.0 - a.0)
    }
}
#[derive(Clone, Debug)]
pub struct OpticalMaterial {
    pub name: String,
    pub n: Spectrum,
    pub k: Spectrum,
    pub source: &'static str,
}
impl OpticalMaterial {
    pub fn from_csv(name: &str, data: &str, source: &'static str) -> Self {
        let rows: Vec<[f64; 3]> = data
            .lines()
            .filter(|s| !s.starts_with('#') && !s.trim().is_empty())
            .map(|s| {
                let v: Vec<f64> = s.split(',').map(|x| x.parse().unwrap()).collect();
                assert_eq!(v.len(), 3);
                [v[0] * 1e-6, v[1], v[2]]
            })
            .collect();
        Self {
            name: name.into(),
            n: Spectrum::new(rows.iter().map(|r| (r[0], r[1])).collect()),
            k: Spectrum::new(rows.iter().map(|r| (r[0], r[2])).collect()),
            source,
        }
    }
    pub fn water() -> Self {
        Self::from_csv(
            "water (298 K)",
            include_str!("../data/water.csv"),
            "Hale & Querry 1973, doi:10.1364/AO.12.000555",
        )
    }
    pub fn iron() -> Self {
        Self::from_csv(
            "iron (room temperature)",
            include_str!("../data/iron.csv"),
            "Johnson & Christy 1974, doi:10.1103/PhysRevB.9.5056",
        )
    }
    /// Intensity absorption coefficient, m^-1.
    pub fn absorption(&self, l: f64) -> f64 {
        4.0 * PI * self.k.at(l) / l
    }
}
/// Planck spectral radiance B_lambda [W m^-2 sr^-1 m^-1], vacuum wavelength.
pub fn planck(l: f64, t: f64) -> f64 {
    assert!(l > 0.0 && l.is_finite() && t >= 0.0 && t.is_finite());
    if t == 0.0 {
        return 0.0;
    }
    let x = H * C / (l * KB * t);
    if x > 700.0 {
        return 0.0;
    }
    2.0 * H * C * C / (l.powi(5) * x.exp_m1())
}
/// Exact unpolarized dielectric Fresnel power reflectance.
pub fn fresnel_dielectric(cos: f64, ni: f64, nt: f64) -> f64 {
    assert!(ni > 0.0 && nt > 0.0);
    let c = cos.abs().clamp(0.0, 1.0);
    let st2 = (ni / nt).powi(2) * (1.0 - c * c);
    if st2 >= 1.0 {
        return 1.0;
    }
    let ct = (1.0 - st2).sqrt();
    let rs = (ni * c - nt * ct) / (ni * c + nt * ct);
    let rp = (nt * c - ni * ct) / (nt * c + ni * ct);
    0.5 * (rs * rs + rp * rp)
}
/// Snell refraction; normal points toward the incident medium.
/// None means total internal reflection.
pub fn refract(incoming: Vec3, normal: Vec3, ni: f64, nt: f64) -> Option<Vec3> {
    assert!(ni > 0.0 && nt > 0.0);
    let cos = -incoming.dot(normal);
    assert!(cos >= -1e-12);
    let eta = ni / nt;
    let st2 = eta * eta * (1.0 - cos * cos);
    if st2 >= 1.0 {
        None
    } else {
        Some((incoming * eta + normal * (eta * cos - (1.0 - st2).sqrt())).unit())
    }
}
/// Exact unpolarized Fresnel for a conductor viewed from vacuum.
pub fn fresnel_conductor(cos: f64, n: f64, k: f64) -> f64 {
    let c = cos.abs().clamp(0.0, 1.0);
    let c2 = c * c;
    let s2 = 1.0 - c2;
    let t0 = n * n - k * k - s2;
    let a2b2 = (t0 * t0 + 4.0 * n * n * k * k).sqrt();
    let a = (0.5 * (a2b2 + t0)).sqrt();
    let rs = (a2b2 + c2 - 2.0 * a * c) / (a2b2 + c2 + 2.0 * a * c);
    let rp =
        rs * (c2 * a2b2 + s2 * s2 - 2.0 * a * c * s2) / (c2 * a2b2 + s2 * s2 + 2.0 * a * c * s2);
    0.5 * (rs + rp)
}
/// HG convention: positive g scatters forward along the incoming ray.
pub fn hg(cos: f64, g: f64) -> f64 {
    (1.0 - g * g) / (4.0 * PI * (1.0 + g * g - 2.0 * g * cos).powf(1.5))
}
pub fn sample_hg(incoming: Vec3, g: f64, rng: &mut Rng) -> Vec3 {
    assert!(g.abs() < 1.0);
    let u = rng.uniform();
    let cos = if g.abs() < 1e-5 {
        2.0 * u - 1.0
    } else {
        let v = (1.0 - g * g) / (1.0 - g + 2.0 * g * u);
        (1.0 + g * g - v * v) / (2.0 * g)
    };
    around(incoming, cos.clamp(-1.0, 1.0), 2.0 * PI * rng.uniform())
}
/// Wyman, Sloan & Shirley (JCGT 2013), asymmetric Gaussian CIE 1931 fit.
pub fn cie_xyz(l: f64) -> [f64; 3] {
    let w = l * 1e9;
    let f = |m: f64, a: f64, b: f64| {
        let z = (w - m) * if w < m { a } else { b };
        (-0.5 * z * z).exp()
    };
    [
        (1.056 * f(599.8, 0.0264, 0.0323) + 0.362 * f(442.0, 0.0624, 0.0374)
            - 0.065 * f(501.1, 0.0490, 0.0382))
        .max(0.0),
        0.821 * f(568.8, 0.0213, 0.0247) + 0.286 * f(530.9, 0.0613, 0.0322),
        1.217 * f(437.0, 0.0845, 0.0278) + 0.681 * f(459.0, 0.0385, 0.0725),
    ]
}
pub fn xyz_to_rgb(v: [f64; 3]) -> [f64; 3] {
    [
        3.2406 * v[0] - 1.5372 * v[1] - 0.4986 * v[2],
        -0.9689 * v[0] + 1.8758 * v[1] + 0.0415 * v[2],
        0.0557 * v[0] - 0.2040 * v[1] + 1.0570 * v[2],
    ]
}
