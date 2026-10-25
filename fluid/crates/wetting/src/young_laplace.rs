//! Axisymmetric capillary statics from the arc-length Young–Laplace equations
//! (Bashforth–Adams).
//!
//! For an axisymmetric surface of revolution with meridian arc length `s`,
//! radial coordinate `r`, axial coordinate `z` and tangent angle `φ` to the
//! horizontal (`dr/ds = cos φ`, `dz/ds = sin φ`), the Young–Laplace equation
//! with hydrostatic pressure is
//!
//! ```text
//!   dφ/ds = 2/b ± Bo z − sin φ / r,      Bo = ρ g b² / γ,
//! ```
//!
//! where `b` is the radius of curvature at the apex. The sign of the gravity
//! term is `+` for a **sessile** drop (apex at the top, `z` measured down) and
//! `−` for a **pendant** drop (apex at the bottom, `z` measured up); this is
//! the standard Bashforth–Adams convention (Bashforth & Adams 1883; see also
//! the OpenDrop / Berry et al. 2015 formulation used as a benchmark below).
//!
//! Lengths are scaled by the apex radius `b`, so the integration is
//! parameter-free in shape; the physical size enters through `Bo`.
//!
//! Volume and surface area of the body of revolution accumulate as
//! `dV/ds = π r² sin φ`, `dA/ds = 2 π r`.

use std::f64::consts::PI;

/// Which hydrostatic branch to integrate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DropKind {
    /// Apex at the top, `z` increasing downwards, gravity term `+Bo z`.
    Sessile,
    /// Apex at the bottom, `z` increasing upwards, gravity term `−Bo z`.
    Pendant,
}

/// One sample along the meridian.
#[derive(Clone, Copy, Debug)]
pub struct Sample {
    pub s: f64,
    pub r: f64,
    pub z: f64,
    pub phi: f64,
    pub volume: f64,
    pub area: f64,
}

/// A discretised axisymmetric drop profile.
#[derive(Clone, Debug)]
pub struct AxisymmetricDrop {
    pub kind: DropKind,
    pub bond: f64,
    pub samples: Vec<Sample>,
}

fn rhs(kind: DropKind, bond: f64, r: f64, z: f64, phi: f64) -> (f64, f64, f64, f64, f64) {
    // sin φ / r → 1 as r → 0 (apex); the start point avoids the singularity but
    // the guard keeps the first steps exact.
    let sph = if r < 1e-12 { 1.0 } else { phi.sin() / r };
    let gravity = match kind {
        DropKind::Sessile => bond * z,
        DropKind::Pendant => -bond * z,
    };
    (
        phi.cos(),
        phi.sin(),
        2.0 + gravity - sph,
        PI * r * r * phi.sin(),
        2.0 * PI * r,
    )
}

impl AxisymmetricDrop {
    /// Integrate from the apex with fixed RK4 step `ds` up to arclength `s_max`
    /// (or until the tangent passes `phi_stop`, if given). `s0` is the small
    /// offset used to start off the apex singularity.
    pub fn integrate(kind: DropKind, bond: f64, ds: f64, s0: f64, s_max: f64, phi_stop: Option<f64>) -> Self {
        let (mut r, mut z, mut phi, mut v, mut a) = (s0, 0.0, s0, 0.0, 0.0);
        let mut s = s0;
        let mut samples = vec![Sample { s, r, z, phi, volume: v, area: a }];
        while s < s_max {
            let k1 = rhs(kind, bond, r, z, phi);
            let k2 = rhs(kind, bond, r + 0.5 * ds * k1.0, z + 0.5 * ds * k1.1, phi + 0.5 * ds * k1.2);
            let k3 = rhs(kind, bond, r + 0.5 * ds * k2.0, z + 0.5 * ds * k2.1, phi + 0.5 * ds * k2.2);
            let k4 = rhs(kind, bond, r + ds * k3.0, z + ds * k3.1, phi + ds * k3.2);
            r += ds / 6.0 * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0);
            z += ds / 6.0 * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1);
            phi += ds / 6.0 * (k1.2 + 2.0 * k2.2 + 2.0 * k3.2 + k4.2);
            v += ds / 6.0 * (k1.3 + 2.0 * k2.3 + 2.0 * k3.3 + k4.3);
            a += ds / 6.0 * (k1.4 + 2.0 * k2.4 + 2.0 * k3.4 + k4.4);
            s += ds;
            samples.push(Sample { s, r, z, phi, volume: v, area: a });
            if let Some(ps) = phi_stop {
                if phi >= ps {
                    break;
                }
            }
            if phi >= PI - 1e-6 || r < 0.0 {
                break;
            }
        }
        AxisymmetricDrop { kind, bond, samples }
    }

    /// A sessile drop integrated to the contact angle `theta` (radians).
    pub fn sessile_to_contact_angle(bond: f64, theta: f64, ds: f64) -> Self {
        Self::integrate(DropKind::Sessile, bond, ds, 1e-8, 1e4, Some(theta))
    }

    /// A pendant drop integrated for `s_max` arclengths.
    pub fn pendant(bond: f64, s_max: f64, ds: f64) -> Self {
        Self::integrate(DropKind::Pendant, bond, ds, 1e-8, s_max, None)
    }

    /// Sample the profile at arclength `s` by linear interpolation.
    pub fn at_s(&self, s: f64) -> Sample {
        let ss = &self.samples;
        if s <= ss[0].s {
            return ss[0];
        }
        if s >= ss[ss.len() - 1].s {
            return ss[ss.len() - 1];
        }
        // samples are uniform in s, so index directly
        let ds = ss[1].s - ss[0].s;
        let idx = ((s - ss[0].s) / ds).floor() as usize;
        let t = (s - ss[idx].s) / ds;
        let a = ss[idx];
        let b = ss[idx + 1];
        Sample {
            s,
            r: a.r + t * (b.r - a.r),
            z: a.z + t * (b.z - a.z),
            phi: a.phi + t * (b.phi - a.phi),
            volume: a.volume + t * (b.volume - a.volume),
            area: a.area + t * (b.area - a.area),
        }
    }

    /// Last sample (the endpoint reached).
    pub fn last(&self) -> Sample {
        *self.samples.last().unwrap()
    }

    /// Physical height and base radius of a sessile drop with apex radius `b`.
    pub fn physical_sessile(b_phys: f64, bond: f64, theta: f64, ds: f64) -> (Self, f64, f64, f64) {
        let d = Self::sessile_to_contact_angle(bond, theta, ds);
        let e = d.last();
        (d, b_phys * e.z, b_phys * e.r, b_phys.powi(3) * e.volume)
    }
}

/// Exact geometry of a zero-gravity spherical cap (a spherical-cap drop in the
/// limit `Bo → 0`).
#[derive(Clone, Copy, Debug)]
pub struct SphericalCap {
    /// Contact angle measured through the liquid (radians).
    pub theta: f64,
    /// Sphere radius.
    pub radius: f64,
}

impl SphericalCap {
    pub fn new(theta: f64, radius: f64) -> Self {
        Self { theta, radius }
    }

    /// Solve for the sphere radius that gives volume `v` at contact angle `theta`.
    pub fn from_volume(theta: f64, v: f64) -> Self {
        let c = theta.cos();
        let shape = 2.0 - 3.0 * c + c.powi(3);
        let radius = (3.0 * v / (PI * shape)).cbrt();
        Self { theta, radius }
    }

    /// Height `h = R (1 − cos θ)`.
    pub fn height(&self) -> f64 {
        self.radius * (1.0 - self.theta.cos())
    }

    /// Base (contact) radius `a = R sin θ`.
    pub fn base_radius(&self) -> f64 {
        self.radius * self.theta.sin()
    }

    /// Volume `V = (π R³/3)(2 − 3 cos θ + cos³ θ)`.
    pub fn volume(&self) -> f64 {
        let c = self.theta.cos();
        PI * self.radius.powi(3) / 3.0 * (2.0 - 3.0 * c + c.powi(3))
    }

    /// Surface area of the spherical cap `A = 2 π R² (1 − cos θ)`.
    pub fn area(&self) -> f64 {
        2.0 * PI * self.radius.powi(2) * (1.0 - self.theta.cos())
    }
}
