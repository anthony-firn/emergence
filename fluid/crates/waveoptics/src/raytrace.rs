//! Exact sequential ray tracing through rotationally symmetric conic
//! surfaces, paraxial first-order properties, and wavefront aberration
//! maps on the exit-pupil reference sphere.
//!
//! A surface with curvature `c = 1/R` and conic constant `k` has the sag
//! equation `c(x² + y²) + (1 + k) c z² − 2z = 0` (vertex at the origin);
//! rays are intersected with it exactly (the near root of a quadratic in
//! the ray parameter), refracted by the vector Snell law, and their optical
//! path `Σ n·ℓ` accumulated. The wavefront aberration of a plane wave
//! entering at entrance-pupil coordinate `u` is the optical path to the
//! reference sphere (centred on the image point, through the paraxial
//! exit-pupil centre) minus that of the chief ray: zero for a perfect
//! system by Fermat's principle.

use crate::materials::Material;

/// One refracting surface; `thickness` is the axial distance to the next
/// vertex and `after` the medium behind it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Surface {
    /// Radius of curvature (m); `f64::INFINITY` for a plane.
    pub radius: f64,
    pub conic: f64,
    /// Clear semi-diameter (m).
    pub semi_diameter: f64,
    pub thickness: f64,
    pub after: Material,
}

impl Surface {
    pub fn curvature(&self) -> f64 {
        if self.radius.is_infinite() {
            0.0
        } else {
            1.0 / self.radius
        }
    }
}

/// A centred system: object-space medium, surfaces, an aperture stop, and
/// the image (sensor) plane.
#[derive(Clone, Debug, PartialEq)]
pub struct System {
    pub object_medium: Material,
    pub surfaces: Vec<Surface>,
    /// The stop is a plane at `stop_z` (from the first vertex) with radius
    /// `stop_radius`, in the medium of the surface preceding it.
    pub stop_z: f64,
    pub stop_radius: f64,
    /// Sensor plane, axial distance from the first vertex.
    pub image_z: f64,
}

/// First-order (paraxial) properties at one wavelength.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Paraxial {
    /// Effective focal length (object-space units; `1/power` in air).
    pub efl: f64,
    /// Back focal distance from the last vertex (geometric).
    pub bfd: f64,
    /// Entrance pupil position (from the first vertex) and radius.
    pub ep_z: f64,
    pub ep_radius: f64,
    /// Exit pupil position (from the first vertex).
    pub xp_z: f64,
}

/// A traced ray: position, direction (unit), optical path.
#[derive(Clone, Copy, Debug)]
pub struct Ray {
    pub p: [f64; 3],
    pub d: [f64; 3],
    pub opl: f64,
}

/// Root of `a t² + b t + c = 0` on the sheet of the conic through the
/// vertex: among the forward roots, the one whose local `z` is smallest in
/// magnitude (cancellation-free quadratic formula).
fn near_root(a: f64, b: f64, c: f64, z0: f64, dz: f64) -> Option<f64> {
    let mut roots = Vec::with_capacity(2);
    if a.abs() < 1e-300 {
        if b != 0.0 {
            roots.push(-c / b);
        }
    } else {
        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            return None;
        }
        let q = -0.5 * (b + b.signum().max(0.0).mul_add(2.0, -1.0) * disc.sqrt());
        if q != 0.0 {
            roots.push(q / a);
            roots.push(c / q);
        } else {
            roots.push(0.0);
        }
    }
    roots
        .into_iter()
        .filter(|t| t.is_finite() && *t >= -1e-9)
        .min_by(|x, y| (z0 + x * dz).abs().total_cmp(&(z0 + y * dz).abs()))
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

impl System {
    /// Vertex positions (from the first vertex).
    pub fn vertices(&self) -> Vec<f64> {
        let mut z = 0.0;
        self.surfaces
            .iter()
            .map(|s| {
                let v = z;
                z += s.thickness;
                v
            })
            .collect()
    }

    fn medium_before(&self, i: usize) -> Material {
        if i == 0 {
            self.object_medium
        } else {
            self.surfaces[i - 1].after
        }
    }

    /// Index of the space containing the stop plane (0 = object space).
    fn stop_space(&self) -> usize {
        let v = self.vertices();
        v.iter().filter(|&&z| z <= self.stop_z).count()
    }

    /// Paraxial y–u trace from surface `from` (exclusive of nothing) with
    /// height `y` and angle `u` in the space before surface `from`, at axial
    /// position `z`; returns `(y, u, z)` after the last surface.
    fn paraxial_forward(&self, lambda: f64, from: usize, mut y: f64, mut u: f64, mut z: f64) -> (f64, f64, f64) {
        let v = self.vertices();
        for i in from..self.surfaces.len() {
            y += u * (v[i] - z);
            z = v[i];
            let (n1, n2) = (self.medium_before(i).index(lambda), self.surfaces[i].after.index(lambda));
            u = (n1 * u - y * (n2 - n1) * self.surfaces[i].curvature()) / n2;
        }
        (y, u, z)
    }

    /// Backward paraxial trace from a point in space `space` (before surface
    /// `space`) to object space; returns `(y, u, z)` at the first vertex.
    fn paraxial_backward(&self, lambda: f64, space: usize, mut y: f64, mut u: f64, mut z: f64) -> (f64, f64, f64) {
        let v = self.vertices();
        for i in (0..space).rev() {
            y += u * (v[i] - z);
            z = v[i];
            let (n1, n2) = (self.medium_before(i).index(lambda), self.surfaces[i].after.index(lambda));
            // forward: n2 u2 = n1 u1 − y φ; invert for u1
            u = (n2 * u + y * (n2 - n1) * self.surfaces[i].curvature()) / n1;
        }
        (y, u, z)
    }

    pub fn paraxial(&self, lambda: f64) -> Paraxial {
        let n_obj = self.object_medium.index(lambda);
        let n_img = self.surfaces.last().unwrap().after.index(lambda);
        // parallel marginal ray, height 1
        let (y, u, z) = self.paraxial_forward(lambda, 0, 1.0, 0.0, -1.0);
        let power = -n_img * u / 1.0;
        let efl = n_obj / power;
        let bfd = -y / u;
        let _ = z;
        // entrance pupil: image of the stop centre in object space
        let s = self.stop_space();
        let (yb, ub, zb) = self.paraxial_backward(lambda, s, 0.0, 1.0, self.stop_z);
        let ep_z = zb - yb / ub;
        // pupil radius: a parallel ray of height h reaches the stop at y_s(h)
        let v = self.vertices();
        let (ys, us, zs) = if s == 0 { (1.0, 0.0, self.stop_z) } else { self.paraxial_to_space(lambda, s, 1.0) };
        let y_at_stop = ys + us * (self.stop_z - zs);
        let ep_radius = self.stop_radius / y_at_stop.abs();
        // exit pupil: image of the stop centre in image space
        let (yx, ux, zx) = self.paraxial_forward(lambda, s, 0.0, 1.0, self.stop_z);
        let xp_z = zx - yx / ux;
        let _ = v;
        Paraxial { efl, bfd, ep_z, ep_radius, xp_z }
    }

    /// Parallel ray of height `h` traced paraxially up to (but not through)
    /// surface `space`: `(y, u, z)` at the vertex of surface `space − 1`.
    fn paraxial_to_space(&self, lambda: f64, space: usize, h: f64) -> (f64, f64, f64) {
        let v = self.vertices();
        let (mut y, mut u, mut z) = (h, 0.0, -1.0);
        for i in 0..space {
            y += u * (v[i] - z);
            z = v[i];
            let (n1, n2) = (self.medium_before(i).index(lambda), self.surfaces[i].after.index(lambda));
            u = (n1 * u - y * (n2 - n1) * self.surfaces[i].curvature()) / n2;
        }
        (y, u, z)
    }

    /// Trace a ray exactly from object space through every surface;
    /// `None` if it misses a surface, exceeds a clear aperture, is blocked
    /// by the stop, or is totally internally reflected.
    pub fn trace(&self, lambda: f64, mut ray: Ray) -> Option<Ray> {
        let v = self.vertices();
        let s_space = self.stop_space();
        for (i, s) in self.surfaces.iter().enumerate() {
            if i == s_space {
                // stop plane lies in the space before surface i
                if !self.through_stop(&ray, lambda, i) {
                    return None;
                }
            }
            let n1 = self.medium_before(i).index(lambda);
            let n2 = s.after.index(lambda);
            let c = s.curvature();
            let k = s.conic;
            let p0 = [ray.p[0], ray.p[1], ray.p[2] - v[i]];
            let d = ray.d;
            let a = c * (d[0] * d[0] + d[1] * d[1]) + (1.0 + k) * c * d[2] * d[2];
            let b = 2.0 * (c * (p0[0] * d[0] + p0[1] * d[1]) + (1.0 + k) * c * p0[2] * d[2] - d[2]);
            let cc = c * (p0[0] * p0[0] + p0[1] * p0[1]) + (1.0 + k) * c * p0[2] * p0[2] - 2.0 * p0[2];
            let t = near_root(a, b, cc, p0[2], d[2])?;
            let q = [p0[0] + t * d[0], p0[1] + t * d[1], p0[2] + t * d[2]];
            if q[0] * q[0] + q[1] * q[1] > s.semi_diameter * s.semi_diameter {
                return None;
            }
            ray.opl += n1 * t;
            // normal ∇F = (2cx, 2cy, 2(1+k)cz − 2), oriented against the ray
            let mut nrm = [c * q[0], c * q[1], (1.0 + k) * c * q[2] - 1.0];
            let l = dot(nrm, nrm).sqrt();
            nrm = [nrm[0] / l, nrm[1] / l, nrm[2] / l];
            if dot(nrm, d) > 0.0 {
                nrm = [-nrm[0], -nrm[1], -nrm[2]];
            }
            let cos_i = -dot(nrm, d);
            let eta = n1 / n2;
            let kk = 1.0 - eta * eta * (1.0 - cos_i * cos_i);
            if kk < 0.0 {
                return None;
            }
            let f = eta * cos_i - kk.sqrt();
            ray.d = [eta * d[0] + f * nrm[0], eta * d[1] + f * nrm[1], eta * d[2] + f * nrm[2]];
            ray.p = [q[0], q[1], q[2] + v[i]];
        }
        if s_space == self.surfaces.len() && !self.through_stop(&ray, lambda, s_space) {
            return None;
        }
        Some(ray)
    }

    fn through_stop(&self, ray: &Ray, _lambda: f64, _space: usize) -> bool {
        if ray.d[2] <= 0.0 {
            return false;
        }
        let t = (self.stop_z - ray.p[2]) / ray.d[2];
        let x = ray.p[0] + t * ray.d[0];
        let y = ray.p[1] + t * ray.d[1];
        x * x + y * y <= self.stop_radius * self.stop_radius
    }

    /// Optical path from the plane `z = z0` (object space) to the
    /// reference sphere centred on the axial image point (sensor plane),
    /// through the paraxial exit-pupil centre, for a ray parallel to the
    /// axis at height `(x, y)`. `None` if vignetted.
    fn path_to_sphere(&self, lambda: f64, par: &Paraxial, z0: f64, x: f64, y: f64) -> Option<f64> {
        let ray = self.trace(lambda, Ray { p: [x, y, z0], d: [0.0, 0.0, 1.0], opl: 0.0 })?;
        let n_img = self.surfaces.last().unwrap().after.index(lambda);
        let centre = [0.0, 0.0, self.image_z];
        let radius = self.image_z - par.xp_z;
        // intersect p + t d with |q − centre| = radius, root nearest the exit pupil
        let w = [ray.p[0] - centre[0], ray.p[1] - centre[1], ray.p[2] - centre[2]];
        let b = dot(w, ray.d);
        let c = dot(w, w) - radius * radius;
        let disc = b * b - c;
        if disc < 0.0 {
            return None;
        }
        let r = disc.sqrt();
        // the sphere point on the exit-pupil side (z < centre for radius > 0)
        let (t1, t2) = (-b - r, -b + r);
        let pick = |t: f64| ray.p[2] + t * ray.d[2];
        let t = if radius > 0.0 {
            if pick(t1) < pick(t2) { t1 } else { t2 }
        } else if pick(t1) > pick(t2) {
            t1
        } else {
            t2
        };
        Some(ray.opl + n_img * t)
    }

    /// Wavefront aberration map `W(u)` (m) for a plane wave along the axis,
    /// sampled on an `n × n` grid over the entrance pupil (spacing
    /// `2 r_EP/n`, cell centres); `None` where the ray is vignetted.
    pub fn wavefront(&self, lambda: f64, n: usize) -> (Paraxial, Vec<Option<f64>>) {
        let par = self.paraxial(lambda);
        let z0 = par.ep_z.min(0.0) - 1e-3;
        let chief = self.path_to_sphere(lambda, &par, z0, 0.0, 0.0).expect("chief ray blocked");
        let h = 2.0 * par.ep_radius / n as f64;
        let w = (0..n * n)
            .map(|q| {
                let (i, j) = (q % n, q / n);
                let x = -par.ep_radius + h * (i as f64 + 0.5);
                let y = -par.ep_radius + h * (j as f64 + 0.5);
                if x * x + y * y > par.ep_radius * par.ep_radius * 1.0000001 {
                    return None;
                }
                self.path_to_sphere(lambda, &par, z0, x, y).map(|p| p - chief)
            })
            .collect();
        (par, w)
    }
}
