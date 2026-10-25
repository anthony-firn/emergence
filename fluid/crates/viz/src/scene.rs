//! Scenes: solids, materials, lights, and the nearest-hit query over them.
//!
//! A scene is a flat list of objects, each a [`cad::Solid`] with a
//! material. Nearest-hit marches every object's own field and keeps the
//! closest, which is correct and simple; each object's bounding box
//! rejects most rays before any field is evaluated, so the cost is the
//! objects a ray actually crosses rather than the whole list.
//!
//! Combining every object into one `Solid::Union` and marching that once
//! would be faster per step but wrong for rendering: a union's field
//! cannot say *which* object was hit, and the material lookup needs that.

use cad::solid::Solid;
use geometry::{dot, norm, sub, P3};
use radiation::math::Ray;
use serde::{Deserialize, Serialize};

use crate::trace::{march, Hit, MarchOptions, Trace};

/// How a surface scatters.
///
/// The physics is `radiation::optics`', so the Fresnel and Planck
/// behaviour here is the already-validated implementation rather than a
/// second one.
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(tag = "surface", rename_all = "snake_case")]
pub enum Material {
    /// Lambertian, with a constant spectral reflectance.
    Diffuse { albedo: f64 },
    /// Lambertian with two reflectances in a checker of the given period.
    Checker { a: f64, b: f64, period: f64 },
    /// Smooth dielectric with a real index: glass, water, a lens.
    Dielectric { index: f64 },
    /// Polished conductor with index `n` and extinction `k`.
    Conductor { index: f64, extinction: f64 },
    /// A grey emitter of the given radiance, in W/(m^2 sr).
    Emitter { radiance: f64 },
    /// A blackbody emitter at the given temperature, scaled.
    Blackbody { temperature: f64, #[serde(default = "one")] scale: f64 },
    /// A **dispersive** dielectric, by its Schott Sellmeier coefficients.
    ///
    /// This is what gives a caustic its colour: the index differs across
    /// the spectrum, so each wavelength focuses in a slightly different
    /// place and the caustic edge is fringed. The index comes from
    /// `waveoptics::materials`, so it is the same dispersion the optics
    /// crates use.
    Dispersive {
        sellmeier_b: [f64; 3],
        sellmeier_c: [f64; 3],
    },
    /// A perfect mirror: reflects everything, with no Fresnel falloff.
    ///
    /// Not physical, and there on purpose — the catacaustic of a spherical
    /// mirror is a statement about *geometry*, and checking it against a
    /// surface that also dims with angle would confound the two.
    Mirror,
}

fn one() -> f64 {
    1.0
}

impl Material {
    /// Reflectance at a point, for the diffuse kinds.
    pub fn albedo_at(&self, p: P3) -> f64 {
        match *self {
            Material::Diffuse { albedo } => albedo,
            Material::Checker { a, b, period } => {
                let f = |x: f64| (x / period).floor() as i64;
                if (f(p[0]) + f(p[1]) + f(p[2])).rem_euclid(2) == 0 {
                    a
                } else {
                    b
                }
            }
            _ => 0.0,
        }
    }

    /// Refractive index at a wavelength, for the materials that have one.
    pub fn refractive_index(&self, lambda: f64) -> Option<f64> {
        match *self {
            Material::Dielectric { index } => Some(index),
            Material::Dispersive {
                sellmeier_b,
                sellmeier_c,
            } => Some(
                waveoptics::materials::Material::Sellmeier {
                    b: sellmeier_b,
                    c: sellmeier_c,
                }
                .index(lambda),
            ),
            _ => None,
        }
    }

    /// N-BK7, the commonest optical glass (Schott data sheet, through
    /// `waveoptics::materials`).
    pub fn bk7() -> Self {
        Material::Dispersive {
            sellmeier_b: [1.039_612_12, 0.231_792_344, 1.010_469_45],
            sellmeier_c: [0.006_000_698_67, 0.020_017_914_4, 103.560_653],
        }
    }

    pub fn is_emitter(&self) -> bool {
        matches!(self, Material::Emitter { .. } | Material::Blackbody { .. })
    }

    /// Emitted radiance at wavelength `lambda`, W/(m^2 sr m).
    pub fn emission(&self, lambda: f64) -> f64 {
        match *self {
            Material::Emitter { radiance } => radiance,
            Material::Blackbody { temperature, scale } => {
                radiation::optics::planck(lambda, temperature) * scale
            }
            _ => 0.0,
        }
    }
}

/// One renderable object.
#[derive(Clone, Debug)]
pub struct Object {
    pub name: String,
    pub solid: Solid,
    pub material: usize,
}

impl Object {
    pub fn new(name: &str, solid: Solid, material: usize) -> Self {
        Self {
            name: name.to_string(),
            solid,
            material,
        }
    }
}

/// A light that can be sampled directly, which is what keeps a render from
/// needing thousands of paths to find a small source.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
#[serde(tag = "light", rename_all = "snake_case")]
pub enum Light {
    /// A direction and an irradiance on a surface facing it: the sun.
    Distant { direction: P3, irradiance: f64 },
    /// An isotropic point source of the given radiant intensity, W/sr.
    Point { position: P3, intensity: f64 },
    /// A disc of radius `radius` at `position` facing `normal`, with the
    /// given radiance.
    Disc {
        position: P3,
        normal: P3,
        radius: f64,
        radiance: f64,
    },
}

impl Light {
    /// Irradiance this light puts on a surface at `p` with normal `n`,
    /// ignoring occlusion, and the direction to sample for shadowing with
    /// the distance to go.
    pub fn sample(&self, p: P3, n: P3) -> (f64, P3, f64) {
        match *self {
            Light::Distant { direction, irradiance } => {
                let to_light = geometry::normalize(geometry::scale(direction, -1.0));
                let cos = dot(n, to_light);
                (irradiance * cos.max(0.0), to_light, f64::INFINITY)
            }
            Light::Point { position, intensity } => {
                let d = sub(position, p);
                let r = norm(d);
                if r <= 0.0 {
                    return (0.0, [0.0, 0.0, 1.0], 0.0);
                }
                let to_light = geometry::scale(d, 1.0 / r);
                let cos = dot(n, to_light);
                (intensity * cos.max(0.0) / (r * r), to_light, r)
            }
            Light::Disc {
                position,
                normal,
                radius,
                radiance,
            } => {
                // On-axis-equivalent solid angle: exact for a point on the
                // disc's axis, and the standard approximation off it.
                let d = sub(position, p);
                let r = norm(d);
                if r <= 0.0 {
                    return (0.0, [0.0, 0.0, 1.0], 0.0);
                }
                let to_light = geometry::scale(d, 1.0 / r);
                let cos_surface = dot(n, to_light).max(0.0);
                let cos_light = dot(geometry::normalize(normal), geometry::scale(to_light, -1.0)).max(0.0);
                let solid_angle =
                    std::f64::consts::PI * radius * radius * cos_light / (r * r + radius * radius);
                (radiance * solid_angle * cos_surface, to_light, r)
            }
        }
    }
}

/// A scene.
#[derive(Clone, Debug, Default)]
pub struct Scene {
    pub objects: Vec<Object>,
    pub materials: Vec<Material>,
    pub lights: Vec<Light>,
    /// Uniform background radiance, W/(m^2 sr).
    pub background: f64,
}

/// The nearest hit in a scene.
#[derive(Clone, Copy, Debug)]
pub struct SceneHit {
    pub object: usize,
    pub hit: Hit,
}

impl Scene {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_material(mut self, m: Material) -> Self {
        self.materials.push(m);
        self
    }

    pub fn with_object(mut self, o: Object) -> Self {
        self.objects.push(o);
        self
    }

    pub fn with_light(mut self, l: Light) -> Self {
        self.lights.push(l);
        self
    }

    /// Overall bounding box of the finite objects.
    pub fn bbox(&self) -> (P3, P3) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for o in &self.objects {
            let (a, b) = o.solid.bbox();
            for k in 0..3 {
                if a[k].is_finite() && b[k].is_finite() {
                    lo[k] = lo[k].min(a[k]);
                    hi[k] = hi[k].max(b[k]);
                }
            }
        }
        (lo, hi)
    }

    /// A length scale for the scene, used to pick march tolerances.
    pub fn scale(&self) -> f64 {
        let (lo, hi) = self.bbox();
        let d = norm(sub(hi, lo));
        if d.is_finite() && d > 0.0 {
            d
        } else {
            1.0
        }
    }

    /// Nearest hit along a ray, and whether any march was exhausted.
    pub fn intersect(&self, ray: Ray, opts: &MarchOptions) -> (Option<SceneHit>, bool, usize) {
        let mut best: Option<SceneHit> = None;
        let mut exhausted = false;
        let mut steps = 0usize;
        for (i, o) in self.objects.iter().enumerate() {
            let mut local = *opts;
            if let Some(b) = &best {
                local.t_max = local.t_max.min(b.hit.t);
            }
            match march(&o.solid, ray, &local) {
                Trace::Hit(h) => {
                    steps += h.steps;
                    if best.is_none_or(|b| h.t < b.hit.t) {
                        best = Some(SceneHit { object: i, hit: h });
                    }
                }
                Trace::Exhausted { steps: s, .. } => {
                    steps += s;
                    exhausted = true;
                }
                Trace::Miss => {}
            }
        }
        (best, exhausted, steps)
    }

    /// Is the path from `p` towards `dir` blocked within `limit`?
    pub fn occluded(&self, p: P3, dir: P3, limit: f64, opts: &MarchOptions) -> bool {
        self.objects
            .iter()
            .any(|o| crate::trace::occluded(&o.solid, p, dir, limit, opts))
    }

    pub fn material_of(&self, hit: &SceneHit) -> Option<&Material> {
        self.materials.get(self.objects[hit.object].material)
    }
}
