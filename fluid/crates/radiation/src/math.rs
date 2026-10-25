use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3(pub [f64; 3]);
impl From<[f64; 3]> for Vec3 {
    fn from(v: [f64; 3]) -> Self {
        Self(v)
    }
}
impl Add for Vec3 {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self(std::array::from_fn(|i| self.0[i] + b.0[i]))
    }
}
impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self(std::array::from_fn(|i| self.0[i] - b.0[i]))
    }
}
impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, b: f64) -> Self {
        Self(self.0.map(|v| v * b))
    }
}
impl Div<f64> for Vec3 {
    type Output = Self;
    fn div(self, b: f64) -> Self {
        self * (1.0 / b)
    }
}
impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self {
        self * (-1.0)
    }
}
impl Vec3 {
    pub fn dot(self, b: Self) -> f64 {
        (0..3).map(|i| self.0[i] * b.0[i]).sum()
    }
    pub fn cross(self, b: Self) -> Self {
        Self([
            self.0[1] * b.0[2] - self.0[2] * b.0[1],
            self.0[2] * b.0[0] - self.0[0] * b.0[2],
            self.0[0] * b.0[1] - self.0[1] * b.0[0],
        ])
    }
    pub fn length(self) -> f64 {
        self.dot(self).sqrt()
    }
    pub fn unit(self) -> Self {
        let n = self.length();
        assert!(n > 0.0 && n.is_finite());
        self / n
    }
    pub fn reflect(self, n: Self) -> Self {
        self - n * (2.0 * self.dot(n))
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Ray {
    pub origin: Vec3,
    pub dir: Vec3,
}
impl Ray {
    pub fn new(origin: Vec3, dir: Vec3) -> Self {
        assert!(origin.0.iter().all(|v| v.is_finite()));
        Self {
            origin,
            dir: dir.unit(),
        }
    }
    pub fn at(self, t: f64) -> Vec3 {
        self.origin + self.dir * t
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Bounds {
    pub lo: Vec3,
    pub hi: Vec3,
}
impl Bounds {
    pub fn interval(self, r: Ray, limit: f64) -> Option<(f64, f64)> {
        self.interval_cached(r, Vec3(r.dir.0.map(|v| 1.0 / v)), limit)
    }
    pub(crate) fn interval_cached(self, r: Ray, inverse: Vec3, limit: f64) -> Option<(f64, f64)> {
        let (mut a, mut b) = (0.0_f64, limit);
        for d in 0..3 {
            if r.dir.0[d] == 0.0 {
                if r.origin.0[d] < self.lo.0[d] || r.origin.0[d] > self.hi.0[d] {
                    return None;
                }
            } else {
                let x = if inverse.0[d].is_finite() {
                    (self.lo.0[d] - r.origin.0[d]) * inverse.0[d]
                } else {
                    (self.lo.0[d] - r.origin.0[d]) / r.dir.0[d]
                };
                let y = if inverse.0[d].is_finite() {
                    (self.hi.0[d] - r.origin.0[d]) * inverse.0[d]
                } else {
                    (self.hi.0[d] - r.origin.0[d]) / r.dir.0[d]
                };
                a = a.max(x.min(y));
                b = b.min(x.max(y));
                if a > b {
                    return None;
                }
            }
        }
        Some((a, b))
    }
    pub fn union(self, b: Self) -> Self {
        Self {
            lo: Vec3(std::array::from_fn(|i| self.lo.0[i].min(b.lo.0[i]))),
            hi: Vec3(std::array::from_fn(|i| self.hi.0[i].max(b.hi.0[i]))),
        }
    }
    pub fn center(self) -> Vec3 {
        (self.lo + self.hi) * 0.5
    }
}
/// SplitMix64: independent, reproducible streams per packet/pixel.
pub struct Rng(pub u64);
impl Rng {
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
    /// Open interval (0,1), safe for logarithms.
    pub fn uniform(&mut self) -> f64 {
        ((self.next_u64() >> 12) as f64 + 0.5) / 4503599627370496.0
    }
}
pub fn around(axis: Vec3, cos: f64, phi: f64) -> Vec3 {
    let a = if axis.0[2].abs() < 0.9 {
        Vec3([0.0, 0.0, 1.0])
    } else {
        Vec3([1.0, 0.0, 0.0])
    };
    let u = axis.cross(a).unit();
    let v = axis.cross(u);
    (axis * cos + (u * phi.cos() + v * phi.sin()) * (1.0 - cos * cos).max(0.0).sqrt()).unit()
}
