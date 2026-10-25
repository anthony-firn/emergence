//! Static, zero-shift isotropic metrics: ds²=-alpha² c²dt²+beta² dx².
//! Exact Schwarzschild exterior and a declared weak-field interior closure.
use crate::{norm, Field, Vec3, C};
#[derive(Clone, Copy, Debug)]
pub struct Geometry {
    pub lapse: f64,
    pub scale: f64,
    pub grad_lapse: Vec3,
    pub grad_scale: Vec3,
}
impl Geometry {
    pub fn validate(self) -> Result<Self, &'static str> {
        if !self.lapse.is_finite()
            || self.lapse <= 0.0
            || !self.scale.is_finite()
            || self.scale <= 0.0
            || self
                .grad_lapse
                .0
                .iter()
                .chain(&self.grad_scale.0)
                .any(|v| !v.is_finite())
        {
            Err("invalid static gravitational metric")
        } else {
            Ok(self)
        }
    }
}
pub trait Metric: Sync {
    fn geometry(&self, position: Vec3) -> Result<Geometry, &'static str>;
}
#[derive(Clone, Copy, Debug)]
pub struct Minkowski;
impl Metric for Minkowski {
    fn geometry(&self, x: Vec3) -> Result<Geometry, &'static str> {
        if x.0.iter().any(|v| !v.is_finite()) {
            return Err("invalid metric coordinate");
        }
        Ok(Geometry {
            lapse: 1.0,
            scale: 1.0,
            grad_lapse: Vec3::default(),
            grad_scale: Vec3::default(),
        })
    }
}
pub struct WeakField<F>(pub F);
impl<F: Field> Metric for WeakField<F> {
    fn geometry(&self, x: Vec3) -> Result<Geometry, &'static str> {
        let s = self.0.sample(x)?.validate()?;
        let f = s.potential / (C * C);
        if f.abs() > 1e-3 {
            return Err("source exceeds declared weak-field metric domain");
        }
        let lapse = (1.0 + 2.0 * f).sqrt();
        let scale = (1.0 - 2.0 * f).sqrt();
        Geometry {
            lapse,
            scale,
            grad_lapse: -s.acceleration / (C * C * lapse),
            grad_scale: s.acceleration / (C * C * scale),
        }
        .validate()
    }
}
/// Vacuum exterior in isotropic coordinates. Radius is NOT the Schwarzschild
/// areal radius: R_areal=r*(1+GM/(2c²r))². No material interior or source spin.
#[derive(Clone, Copy, Debug)]
pub struct Schwarzschild {
    pub gm: f64,
    pub center: Vec3,
}
impl Metric for Schwarzschild {
    fn geometry(&self, x: Vec3) -> Result<Geometry, &'static str> {
        if !self.gm.is_finite()
            || self.gm < 0.0
            || x.0.iter().chain(&self.center.0).any(|v| !v.is_finite())
        {
            return Err("invalid Schwarzschild coordinates/source");
        }
        if self.gm == 0.0 {
            return Minkowski.geometry(x);
        }
        let d = x - self.center;
        let r = norm(d);
        let u = self.gm / (2.0 * C * C * r);
        if !u.is_finite() || u >= 1.0 {
            return Err("outside static Schwarzschild coordinate domain");
        }
        let p = 1.0 + u;
        Geometry {
            lapse: (1.0 - u) / p,
            scale: p * p,
            grad_lapse: (d / r) * (2.0 * u / (r * p * p)),
            grad_scale: -(d / r) * (2.0 * u * p / r),
        }
        .validate()
    }
}
