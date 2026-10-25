//! Optical elements: lenses, mirrors, prisms and windows.
//!
//! An optical surface is defined by a sag equation, not by a drawing, so
//! this kit can hand the *same* numbers to two different consumers and
//! they must agree: a [`Solid`] for a mechanical or thermal model, and a
//! surface list for `waveoptics::raytrace`, which traces real rays through
//! it. The tests use that to check the geometry against the paraxial
//! formulas and against the one case where ray optics is exact — a
//! paraboloid, which brings every parallel ray to a single point.
//!
//! ```text
//!   Lens { r1, r2, thickness, semi_diameter, conics }
//!        -> meridian sag profile -> Region2 -> revolve -> Solid
//!        -> waveoptics::raytrace::Surface list -> real rays
//!        -> thick-lens focal length, edge thickness, exact volume
//! ```
//!
//! Lengths are metres. The sign convention is the standard one: a radius
//! is positive when its centre of curvature lies to the right of the
//! vertex (down the direction of propagation, `+z`), the front vertex is
//! at `z = 0`, and the back vertex at `z = thickness`.
//!
//! # The sag equation
//!
//! ```text
//!   z(r) = r^2 / (R (1 + sqrt(1 - (1 + k) r^2 / R^2)))
//!          + sum a_i r^(2i+2)
//! ```
//!
//! the standard conic-plus-even-asphere form (ISO 10110-12; Smith, *Modern
//! Optical Engineering*, 4th ed., §1.4). `k = 0` is a sphere, `k = -1` a
//! paraboloid, `-1 < k < 0` a prolate ellipsoid and `k < -1` a hyperboloid.
//!
//! # What is exact
//!
//! * **A paraboloid has no spherical aberration.** Every ray parallel to
//!   the axis reflects through `R/2` exactly, whatever its height, and the
//!   test checks that to 1e-14 — which validates the sag *and* its normal.
//!   A sphere's marginal focus falls short of `R/2`, and the test measures
//!   how far.
//! * **The thick-lens focal length** is
//!   `1/f = (n-1)[1/R1 - 1/R2 + (n-1) d/(n R1 R2)]`, and the paraxial
//!   trace in `waveoptics` must reproduce it from the generated surfaces.
//!   Two independent implementations of first-order optics, compared.
//! * **A prism's minimum deviation** is `2 asin(n sin(A/2)) - A`, and the
//!   test finds the minimum of the traced deviation numerically and
//!   compares.
//! * **Edge thickness and volume** come from the sag analytically; the
//!   test checks the volume against the closed form for a spherical cap.
//!
//! # Limits
//!
//! Centred systems only: no tilts, decentres or freeform surfaces. No
//! coatings (`photonics::multilayer` does those), no tolerancing, no
//! mounting features, and nothing here chooses a glass or corrects an
//! aberration.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::sketch::{Polygon, Region2, V2};
use crate::solid::Solid;

/// One conic-plus-asphere surface.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Sag {
    /// Vertex radius of curvature; `inf` for a plane.
    pub radius: f64,
    /// Conic constant `k`.
    #[serde(default)]
    pub conic: f64,
    /// Even aspheric coefficients `a_i` on `r^(2i+2)`, starting at
    /// `r^4`.
    #[serde(default)]
    pub aspheric: Vec<f64>,
}

impl Sag {
    pub fn sphere(radius: f64) -> Self {
        Self { radius, conic: 0.0, aspheric: Vec::new() }
    }
    pub fn plane() -> Self {
        Self { radius: f64::INFINITY, conic: 0.0, aspheric: Vec::new() }
    }
    pub fn parabola(radius: f64) -> Self {
        Self { radius, conic: -1.0, aspheric: Vec::new() }
    }
    pub fn conic(radius: f64, conic: f64) -> Self {
        Self { radius, conic, aspheric: Vec::new() }
    }

    /// Sag at radial distance `r`. `None` when the surface does not reach
    /// that far (the conic term goes imaginary).
    pub fn at(&self, r: f64) -> Option<f64> {
        if !self.radius.is_finite() {
            let mut z = 0.0;
            for (i, a) in self.aspheric.iter().enumerate() {
                z += a * r.powi(2 * i as i32 + 2);
            }
            return Some(z);
        }
        let c = 1.0 / self.radius;
        let arg = 1.0 - (1.0 + self.conic) * c * c * r * r;
        if arg < 0.0 {
            return None;
        }
        let mut z = c * r * r / (1.0 + arg.sqrt());
        for (i, a) in self.aspheric.iter().enumerate() {
            z += a * r.powi(2 * i as i32 + 2);
        }
        Some(z)
    }

    /// `dz/dr` at `r`, analytically.
    pub fn slope(&self, r: f64) -> Option<f64> {
        let mut s = 0.0;
        if self.radius.is_finite() {
            let c = 1.0 / self.radius;
            let arg = 1.0 - (1.0 + self.conic) * c * c * r * r;
            if arg < 0.0 {
                return None;
            }
            // d/dr of c r^2 / (1 + sqrt(arg)).
            s += c * r / arg.sqrt();
        }
        for (i, a) in self.aspheric.iter().enumerate() {
            let p = 2 * i as i32 + 2;
            s += a * p as f64 * r.powi(p - 1);
        }
        Some(s)
    }

    /// Outward unit normal in the meridian plane at radius `r`.
    ///
    /// It points toward `+z`, which for a surface with `radius > 0` is the
    /// open, concave side: the side light arrives from.
    pub fn normal(&self, r: f64) -> Option<V2> {
        let s = self.slope(r)?;
        let n = (1.0 + s * s).sqrt();
        Some([-s / n, 1.0 / n])
    }

    /// Largest `r` the surface reaches.
    pub fn max_radius(&self) -> f64 {
        if !self.radius.is_finite() {
            return f64::INFINITY;
        }
        let k = 1.0 + self.conic;
        if k <= 0.0 {
            f64::INFINITY
        } else {
            self.radius.abs() / k.sqrt()
        }
    }
}

/// A rotationally symmetric lens.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Lens {
    pub front: Sag,
    pub back: Sag,
    /// Axial thickness between the vertices.
    pub thickness: f64,
    pub semi_diameter: f64,
    /// Points per surface in the meridian profile.
    #[serde(default = "profile_points")]
    pub points: usize,
}

fn profile_points() -> usize {
    80
}

impl Lens {
    pub fn new(r1: f64, r2: f64, thickness: f64, semi_diameter: f64) -> Self {
        Self {
            front: Sag::sphere(r1),
            back: Sag::sphere(r2),
            thickness,
            semi_diameter,
            points: profile_points(),
        }
    }

    pub fn plano_convex(radius: f64, thickness: f64, semi_diameter: f64) -> Self {
        Self {
            front: Sag::sphere(radius),
            back: Sag::plane(),
            thickness,
            semi_diameter,
            points: profile_points(),
        }
    }

    /// Edge thickness: axial thickness plus the back sag less the front
    /// sag at the clear aperture. Negative means a knife edge.
    pub fn edge_thickness(&self) -> Result<f64, String> {
        let a = self.semi_diameter;
        let sf = self.front.at(a).ok_or("the front surface does not reach the clear aperture")?;
        let sb = self.back.at(a).ok_or("the back surface does not reach the clear aperture")?;
        Ok(self.thickness + sb - sf)
    }

    /// Thick-lens focal length in a medium of index `n_medium`
    /// (the lensmaker's equation).
    pub fn focal_length(&self, n: f64, n_medium: f64) -> f64 {
        let nr = n / n_medium;
        let (c1, c2) = (
            if self.front.radius.is_finite() { 1.0 / self.front.radius } else { 0.0 },
            if self.back.radius.is_finite() { 1.0 / self.back.radius } else { 0.0 },
        );
        let power = (nr - 1.0) * (c1 - c2 + (nr - 1.0) * self.thickness * c1 * c2 / nr);
        1.0 / power
    }

    /// Front and back focal distances from the respective vertices.
    pub fn focal_distances(&self, n: f64, n_medium: f64) -> (f64, f64) {
        let f = self.focal_length(n, n_medium);
        let nr = n / n_medium;
        let (c1, c2) = (
            if self.front.radius.is_finite() { 1.0 / self.front.radius } else { 0.0 },
            if self.back.radius.is_finite() { 1.0 / self.back.radius } else { 0.0 },
        );
        // Principal-plane offsets (Smith, Modern Optical Engineering §2.5).
        let ffd = f * (1.0 + (nr - 1.0) * self.thickness * c2 / nr);
        let bfd = f * (1.0 - (nr - 1.0) * self.thickness * c1 / nr);
        (ffd, bfd)
    }

    /// The meridian section, spanning the axis so the revolved solid's
    /// axis is interior (see [`Solid::Revolve`]).
    pub fn meridian(&self) -> Result<Polygon, String> {
        let a = self.semi_diameter;
        if a <= 0.0 || self.thickness <= 0.0 {
            return Err("a lens needs a positive thickness and aperture".into());
        }
        if self.edge_thickness()? <= 0.0 {
            return Err(format!(
                "the lens closes to a knife edge: edge thickness {:e} at semi-diameter {a}",
                self.edge_thickness()?
            ));
        }
        let n = self.points.max(8);
        let mut pts: Vec<V2> = Vec::with_capacity(2 * n + 4);
        for i in 0..=n {
            let r = -a + 2.0 * a * i as f64 / n as f64;
            let z = self.front.at(r.abs()).ok_or("front surface too steep")?;
            pts.push([r, z]);
        }
        for i in (0..=n).rev() {
            let r = -a + 2.0 * a * i as f64 / n as f64;
            let z = self.thickness + self.back.at(r.abs()).ok_or("back surface too steep")?;
            pts.push([r, z]);
        }
        Ok(Polygon::new(pts))
    }

    /// The lens as a solid of revolution about `z`.
    pub fn solid(&self) -> Result<Solid, String> {
        Ok(Solid::revolve(Region2::Polygon(self.meridian()?)))
    }

    /// Glass volume, by Gauss quadrature of `2 pi r (z_back - z_front) dr`
    /// over the clear aperture. The integrand is smooth, so 5-point
    /// quadrature on 64 intervals is good to well past single precision.
    pub fn volume(&self) -> Result<f64, String> {
        const X: [f64; 5] = [
            -0.906_179_845_938_664,
            -0.538_469_310_105_683,
            0.0,
            0.538_469_310_105_683,
            0.906_179_845_938_664,
        ];
        const W: [f64; 5] = [
            0.236_926_885_056_189,
            0.478_628_670_499_366,
            0.568_888_888_888_889,
            0.478_628_670_499_366,
            0.236_926_885_056_189,
        ];
        let a = self.semi_diameter;
        let n = 64;
        let h = a / n as f64;
        let mut total = 0.0;
        for i in 0..n {
            let (lo, hi) = (h * i as f64, h * (i + 1) as f64);
            let (mid, half) = (0.5 * (lo + hi), 0.5 * (hi - lo));
            for k in 0..5 {
                let r = mid + half * X[k];
                let zf = self.front.at(r).ok_or("front surface too steep")?;
                let zb = self.back.at(r).ok_or("back surface too steep")?;
                total += W[k] * 2.0 * PI * r * (self.thickness + zb - zf) * half;
            }
        }
        Ok(total)
    }

    pub fn mass(&self, density: f64) -> Result<f64, String> {
        Ok(self.volume()? * density)
    }

    /// The two surfaces as `(radius, conic, semi_diameter, thickness)`
    /// tuples, in the order and units `waveoptics::raytrace::Surface`
    /// takes. The caller supplies the materials.
    pub fn surface_data(&self) -> [(f64, f64, f64, f64); 2] {
        [
            (self.front.radius, self.front.conic, self.semi_diameter, self.thickness),
            (self.back.radius, self.back.conic, self.semi_diameter, 0.0),
        ]
    }
}

/// A conic mirror.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Mirror {
    pub surface: Sag,
    pub semi_diameter: f64,
    /// Substrate thickness behind the vertex.
    pub thickness: f64,
    #[serde(default = "profile_points")]
    pub points: usize,
}

impl Mirror {
    pub fn paraboloid(radius: f64, semi_diameter: f64, thickness: f64) -> Self {
        Self {
            surface: Sag::parabola(radius),
            semi_diameter,
            thickness,
            points: profile_points(),
        }
    }

    pub fn spherical(radius: f64, semi_diameter: f64, thickness: f64) -> Self {
        Self {
            surface: Sag::sphere(radius),
            semi_diameter,
            thickness,
            points: profile_points(),
        }
    }

    /// Paraxial focal length, `R/2`.
    pub fn focal_length(&self) -> f64 {
        0.5 * self.surface.radius
    }

    /// Where a ray parallel to the axis at height `h` crosses the axis
    /// after reflection, measured from the vertex.
    ///
    /// For a paraboloid this is `R/2` for every `h`; for a sphere it falls
    /// short, which is the longitudinal spherical aberration.
    ///
    /// The ray travels along `-z`, which for `radius > 0` is the direction
    /// that strikes the concave face. Getting that backwards gives the
    /// same number by a *negative* path length — a virtual focus behind
    /// the mirror — so the returned distance is checked to be forward.
    pub fn marginal_focus(&self, h: f64) -> Option<f64> {
        let z = self.surface.at(h)?;
        let n = self.surface.normal(h)?; // unit, in (r, z)
        let d = [0.0, -1.0];
        let dn = d[0] * n[0] + d[1] * n[1];
        let r = [d[0] - 2.0 * dn * n[0], d[1] - 2.0 * dn * n[1]];
        if r[0].abs() < 1e-300 {
            return None;
        }
        // Travel from (h, z) until the radial coordinate reaches zero.
        let t = -h / r[0];
        if t <= 0.0 {
            return None; // the reflected ray diverges: no real focus
        }
        Some(z + t * r[1])
    }

    /// The meridian section of the substrate.
    ///
    /// The substrate sits *behind* the reflecting face, on the side light
    /// does not come from, so that the face is the first surface a ray
    /// meets. Putting it on the other side fills the dish with glass and
    /// the mirror reflects off a flat back instead — which is how this was
    /// found, by a traced paraboloid that failed to focus.
    pub fn meridian(&self) -> Result<Polygon, String> {
        let a = self.semi_diameter;
        let n = self.points.max(8);
        let edge = self.surface.at(a).ok_or("the surface does not reach the aperture")?;
        // Uniform thickness behind the deepest point of the face: the
        // vertex for a concave surface, the rim for a convex one.
        let back = edge.min(0.0) - self.thickness;
        let mut pts: Vec<V2> = Vec::with_capacity(n + 4);
        for i in 0..=n {
            let r = -a + 2.0 * a * i as f64 / n as f64;
            pts.push([r, self.surface.at(r.abs()).ok_or("surface too steep")?]);
        }
        pts.push([a, back]);
        pts.push([-a, back]);
        Ok(Polygon::new(pts))
    }

    pub fn solid(&self) -> Result<Solid, String> {
        Ok(Solid::revolve(Region2::Polygon(self.meridian()?)))
    }
}

/// A prism of apex angle `apex`, extruded along `+z`.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct Prism {
    /// Apex angle, radians.
    pub apex: f64,
    /// Length of the two refracting faces.
    pub face: f64,
    /// Extrusion depth.
    pub depth: f64,
}

impl Prism {
    pub fn equilateral(face: f64, depth: f64) -> Self {
        Self { apex: PI / 3.0, face, depth }
    }

    pub fn right_angle(leg: f64, depth: f64) -> Self {
        Self { apex: PI / 2.0, face: leg, depth }
    }

    /// Cross-section: an isosceles triangle with its apex up, base on
    /// `y = 0`.
    pub fn region(&self) -> Region2 {
        let half_base = self.face * (0.5 * self.apex).sin();
        let height = self.face * (0.5 * self.apex).cos();
        Region2::Polygon(Polygon::new(vec![
            [-half_base, 0.0],
            [half_base, 0.0],
            [0.0, height],
        ]))
    }

    pub fn solid(&self) -> Solid {
        Solid::extrude(self.region(), 0.0, self.depth)
    }

    /// Deviation of a ray entering the first face at incidence `i1`,
    /// by Snell's law applied twice. `None` when the ray is totally
    /// internally reflected at the second face.
    pub fn deviation(&self, n: f64, i1: f64) -> Option<f64> {
        let r1 = (i1.sin() / n).asin();
        let r2 = self.apex - r1;
        let s = n * r2.sin();
        if s.abs() > 1.0 {
            return None; // total internal reflection
        }
        let i2 = s.asin();
        Some(i1 + i2 - self.apex)
    }

    /// Minimum deviation, `2 asin(n sin(A/2)) - A`.
    pub fn minimum_deviation(&self, n: f64) -> Option<f64> {
        let s = n * (0.5 * self.apex).sin();
        if s.abs() > 1.0 {
            return None;
        }
        Some(2.0 * s.asin() - self.apex)
    }

    /// Incidence that gives minimum deviation (the symmetric passage).
    pub fn minimum_deviation_incidence(&self, n: f64) -> Option<f64> {
        let s = n * (0.5 * self.apex).sin();
        if s.abs() > 1.0 {
            return None;
        }
        Some(s.asin())
    }

    /// Largest apex angle a ray can pass at all:
    /// `2 asin(1/n)`, beyond which the second face totally reflects.
    pub fn maximum_apex(n: f64) -> f64 {
        2.0 * (1.0 / n).asin()
    }
}

/// A plane-parallel window.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct Window {
    pub thickness: f64,
    pub semi_diameter: f64,
}

impl Window {
    pub fn solid(&self) -> Solid {
        Solid::cylinder(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, self.thickness],
            self.semi_diameter,
        )
    }

    /// Axial image shift a window introduces, `t (1 - 1/n)`.
    pub fn focus_shift(&self, n: f64) -> f64 {
        self.thickness * (1.0 - 1.0 / n)
    }
}
