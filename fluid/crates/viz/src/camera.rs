//! Cameras: pinhole, thin lens, and orthographic.
//!
//! The pinhole camera is exact and invertible: [`Camera::project`] is the
//! inverse of [`Camera::ray`], and the tests check the round trip to
//! 1e-15. That matters for more than tidiness — a projection that does not
//! invert its own rays cannot be used to overlay a measurement on a
//! render, or to tell which pixel a probe point falls in.
//!
//! The thin-lens camera adds a finite aperture, so the image has a depth
//! of field. Rays leave a disc of radius `aperture/2` and are aimed so
//! that every ray through a given image point meets at the focus
//! distance — exactly, by construction, which the tests also check. The
//! blur circle of an out-of-focus point then follows the textbook
//! `aperture * |1/focus - 1/distance| * focus` relation.
//!
//! For rendering *through* a real optical system — a lens with its own
//! dispersion and aberrations — use `cad::kits::optic` with
//! `waveoptics::raytrace` instead. This camera is the idealisation.

use geometry::{add, cross, dot, norm, normalize, scale, sub, P3};
use radiation::math::{Ray, Rng, Vec3};
use serde::{Deserialize, Serialize};

/// A camera looking from `eye` at `target`.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct Camera {
    pub eye: P3,
    pub target: P3,
    pub up: P3,
    /// Vertical field of view, radians. Ignored when `orthographic` is
    /// set.
    pub fov_y: f64,
    /// Entrance pupil diameter; zero for a pinhole.
    #[serde(default)]
    pub aperture: f64,
    /// Distance at which the aperture's rays converge.
    #[serde(default)]
    pub focus: f64,
    /// Vertical extent of an orthographic view; `None` for perspective.
    #[serde(default)]
    pub orthographic: Option<f64>,
}

/// The camera's own axes: `(right, up, forward)`, orthonormal.
#[derive(Clone, Copy, Debug)]
pub struct Frame {
    pub right: P3,
    pub up: P3,
    pub forward: P3,
}

impl Camera {
    pub fn looking_at(eye: P3, target: P3, fov_y: f64) -> Self {
        Self {
            eye,
            target,
            up: [0.0, 0.0, 1.0],
            fov_y,
            aperture: 0.0,
            focus: 0.0,
            orthographic: None,
        }
    }

    pub fn with_up(mut self, up: P3) -> Self {
        self.up = up;
        self
    }

    /// Give the camera a finite aperture focused at `focus`.
    pub fn with_depth_of_field(mut self, aperture: f64, focus: f64) -> Self {
        self.aperture = aperture;
        self.focus = focus;
        self
    }

    pub fn orthographic(mut self, height: f64) -> Self {
        self.orthographic = Some(height);
        self
    }

    /// Distance from the eye to the target.
    pub fn subject_distance(&self) -> f64 {
        norm(sub(self.target, self.eye))
    }

    /// The orthonormal frame. `up` is orthogonalised against the view
    /// direction; a camera whose `up` is parallel to its view direction
    /// has no frame and gets one from the least-aligned axis instead.
    pub fn frame(&self) -> Frame {
        let forward = normalize(sub(self.target, self.eye));
        let mut up = sub(self.up, scale(forward, dot(self.up, forward)));
        if norm(up) < 1e-12 * norm(self.up).max(1e-30) {
            let k = (0..3)
                .min_by(|&i, &j| forward[i].abs().partial_cmp(&forward[j].abs()).unwrap())
                .unwrap();
            let mut e = [0.0; 3];
            e[k] = 1.0;
            up = sub(e, scale(forward, dot(e, forward)));
        }
        let up = normalize(up);
        Frame {
            right: cross(forward, up),
            up,
            forward,
        }
    }

    /// The ray through normalised image coordinates `(u, v)`, both in
    /// `[-1, 1]` with `v` up and `u` right, for an image of the given
    /// aspect ratio (width over height).
    ///
    /// `rng` is used only when the aperture is finite.
    pub fn ray(&self, u: f64, v: f64, aspect: f64, rng: &mut Rng) -> Ray {
        let f = self.frame();
        if let Some(height) = self.orthographic {
            let half_h = 0.5 * height;
            let half_w = half_h * aspect;
            let origin = add(
                self.eye,
                add(scale(f.right, u * half_w), scale(f.up, v * half_h)),
            );
            return Ray::new(Vec3(origin), Vec3(f.forward));
        }
        let half_h = (0.5 * self.fov_y).tan();
        let dir = normalize(add(
            f.forward,
            add(scale(f.right, u * half_h * aspect), scale(f.up, v * half_h)),
        ));
        if self.aperture <= 0.0 || self.focus <= 0.0 {
            return Ray::new(Vec3(self.eye), Vec3(dir));
        }
        // A point on the aperture disc, and the focal point this image
        // direction must pass through.
        let (r, phi) = (
            0.5 * self.aperture * rng.uniform().sqrt(),
            std::f64::consts::TAU * rng.uniform(),
        );
        let offset = add(scale(f.right, r * phi.cos()), scale(f.up, r * phi.sin()));
        let origin = add(self.eye, offset);
        // The focus plane is perpendicular to `forward` at `focus`, so the
        // image direction reaches it at this point.
        let focal_point = add(self.eye, scale(dir, self.focus / dot(dir, f.forward)));
        Ray::new(Vec3(origin), Vec3(normalize(sub(focal_point, origin))))
    }

    /// Where a point in space falls in normalised image coordinates, and
    /// its distance along the view axis. `None` for a point behind the
    /// camera.
    ///
    /// This is the exact inverse of [`Camera::ray`] for a pinhole camera.
    pub fn project(&self, p: P3, aspect: f64) -> Option<([f64; 2], f64)> {
        let f = self.frame();
        let d = sub(p, self.eye);
        let z = dot(d, f.forward);
        if let Some(height) = self.orthographic {
            let half_h = 0.5 * height;
            return Some((
                [dot(d, f.right) / (half_h * aspect), dot(d, f.up) / half_h],
                z,
            ));
        }
        if z <= 0.0 {
            return None;
        }
        let half_h = (0.5 * self.fov_y).tan();
        Some((
            [
                dot(d, f.right) / (z * half_h * aspect),
                dot(d, f.up) / (z * half_h),
            ],
            z,
        ))
    }

    /// Pixel centre `(x, y)` in normalised coordinates, with `y`
    /// increasing downward as images are stored.
    pub fn pixel(&self, x: usize, y: usize, width: usize, height: usize, jitter: [f64; 2]) -> [f64; 2] {
        [
            2.0 * (x as f64 + jitter[0]) / width as f64 - 1.0,
            1.0 - 2.0 * (y as f64 + jitter[1]) / height as f64,
        ]
    }

    /// Diameter of the blur circle of a point at `distance`, projected
    /// back to the focus plane: `A |1/focus - 1/distance| focus`.
    pub fn blur_diameter(&self, distance: f64) -> f64 {
        if self.aperture <= 0.0 || self.focus <= 0.0 || distance <= 0.0 {
            return 0.0;
        }
        self.aperture * (1.0 / self.focus - 1.0 / distance).abs() * self.focus
    }

    /// Horizontal field of view for a given aspect ratio.
    pub fn fov_x(&self, aspect: f64) -> f64 {
        2.0 * ((0.5 * self.fov_y).tan() * aspect).atan()
    }
}
