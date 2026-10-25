//! RGB surface transport for photographic CAD scenes.
//!
//! GGX visible-normal sampling: Heitz, JCGT 7(4), 2018, pp. 1–13.
//! https://jcgt.org/published/0007/04/01/
//! Smith masking, dielectric transport and MIS: PBRT 4e, chapters 9 and 13.
//! RGB albedo/indices are closures, not a spectral or polarization solver.
//! Coated diffuse neglects interreflection inside the coating; GGX is single
//! scattering. Normal maps and interpolated normals are shading approximations.

use crate::{Camera, Image};
use geometry::{bvh::Bvh, TriMesh};
use radiation::{
    math::{around, Bounds, Ray, Rng, Vec3},
    optics::{fresnel_conductor, fresnel_dielectric, refract},
};
use rayon::prelude::*;
use std::{
    f64::consts::PI,
    path::Path,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    },
};

pub type Rgb = [f64; 3];
fn mul(a: Rgb, b: Rgb) -> Rgb {
    std::array::from_fn(|k| a[k] * b[k])
}
fn add(a: &mut Rgb, b: Rgb) {
    for (x, y) in a.iter_mut().zip(b) {
        *x += y;
    }
}
pub fn luminance(p: Rgb) -> f64 {
    0.2126 * p[0] + 0.7152 * p[1] + 0.0722 * p[2]
}

/// Homogeneous medium; absorption is in inverse metres. No scattering closure.
#[derive(Clone, Copy, Debug)]
pub struct Medium {
    pub index: f64,
    pub absorption: Rgb,
}
impl Medium {
    /// Estimate: air refractive index rounded to unity.
    pub const AIR: Self = Self {
        index: 1.,
        absorption: [0.; 3],
    };
    pub fn transmission(self, distance: f64) -> Rgb {
        self.absorption.map(|s| (-s * distance).exp())
    }
}

#[derive(Clone, Copy, Debug)]
pub enum Bsdf {
    Diffuse {
        albedo: Rgb,
    },
    Coated {
        albedo: Rgb,
        index: f64,
        alpha: [f64; 2],
    },
    Conductor {
        index: Rgb,
        extinction: Rgb,
        alpha: [f64; 2],
    },
    Dielectric {
        outside: Medium,
        inside: Medium,
    },
    Emitter {
        radiance: Rgb,
    },
}

/// Trowbridge–Reitz anisotropic projected-normal density.
pub fn ggx_d(m: Vec3, alpha: [f64; 2]) -> f64 {
    if m.0[2] <= 0. {
        return 0.;
    }
    let q = (m.0[0] / alpha[0]).powi(2) + (m.0[1] / alpha[1]).powi(2) + m.0[2] * m.0[2];
    1. / (PI * alpha[0] * alpha[1] * q * q)
}
fn ggx_lambda(v: Vec3, a: [f64; 2]) -> f64 {
    let c = v.0[2].abs();
    if c == 0. {
        return f64::INFINITY;
    }
    0.5 * ((1. + ((a[0] * v.0[0]).powi(2) + (a[1] * v.0[1]).powi(2)) / (c * c)).sqrt() - 1.)
}
pub fn ggx_visible_normal(v: Vec3, a: [f64; 2], u: [f64; 2]) -> Vec3 {
    let vh = Vec3([a[0] * v.0[0], a[1] * v.0[1], v.0[2]]).unit();
    let t1 = if vh.0[2] < 0.99999 {
        Vec3([0., 0., 1.]).cross(vh).unit()
    } else {
        Vec3([1., 0., 0.])
    };
    let t2 = vh.cross(t1);
    let r = u[0].sqrt();
    let phi = 2. * PI * u[1];
    let x = r * phi.cos();
    let blend = 0.5 * (1. + vh.0[2]);
    let y = (1. - blend) * (1. - x * x).max(0.).sqrt() + blend * r * phi.sin();
    let nh = t1 * x + t2 * y + vh * (1. - x * x - y * y).max(0.).sqrt();
    Vec3([a[0] * nh.0[0], a[1] * nh.0[1], nh.0[2].max(1e-12)]).unit()
}
fn visible_pdf(v: Vec3, m: Vec3, a: [f64; 2]) -> f64 {
    ggx_d(m, a) * v.dot(m).max(0.) / ((1. + ggx_lambda(v, a)) * v.0[2])
}
impl Bsdf {
    fn glossy_probability(self) -> f64 {
        match self {
            Self::Coated { .. } => 0.5,
            Self::Conductor { .. } => 1.,
            _ => 0.,
        }
    }
    /// BRDF and sampling density in a local frame with +z along the normal.
    pub fn evaluate(self, wo: Vec3, wi: Vec3) -> (Rgb, f64) {
        if wo.0[2] <= 0. || wi.0[2] <= 0. {
            return ([0.; 3], 0.);
        }
        let (albedo, alpha) = match self {
            Self::Diffuse { albedo } => return (albedo.map(|v| v / PI), wi.0[2] / PI),
            Self::Coated { albedo, alpha, .. } => (albedo, alpha),
            Self::Conductor { alpha, .. } => ([0.; 3], alpha),
            _ => return ([0.; 3], 0.),
        };
        let m = (wo + wi).unit();
        let fresnel = match self {
            Self::Coated { index, .. } => [fresnel_dielectric(wo.dot(m).abs(), 1., index); 3],
            Self::Conductor {
                index, extinction, ..
            } => {
                std::array::from_fn(|k| fresnel_conductor(wo.dot(m).abs(), index[k], extinction[k]))
            }
            _ => unreachable!(),
        };
        let g = 1. / (1. + ggx_lambda(wo, alpha) + ggx_lambda(wi, alpha));
        let scale = ggx_d(m, alpha) * g / (4. * wo.0[2] * wi.0[2]);
        let attenuation = match self {
            Self::Coated { index, .. } => {
                (1. - fresnel_dielectric(wo.0[2], 1., index))
                    * (1. - fresnel_dielectric(wi.0[2], 1., index))
            }
            _ => 0.,
        };
        let f = std::array::from_fn(|k| fresnel[k] * scale + albedo[k] * attenuation / PI);
        let probability = self.glossy_probability();
        let pdf = probability * visible_pdf(wo, m, alpha) / (4. * wo.dot(m).abs())
            + (1. - probability) * wi.0[2] / PI;
        (f, pdf)
    }
    pub fn sample(self, wo: Vec3, rng: &mut Rng) -> Option<(Vec3, Rgb, f64)> {
        let probability = self.glossy_probability();
        let wi = if probability > 0. && rng.uniform() < probability {
            let alpha = match self {
                Self::Coated { alpha, .. } | Self::Conductor { alpha, .. } => alpha,
                _ => unreachable!(),
            };
            let m = ggx_visible_normal(wo, alpha, [rng.uniform(), rng.uniform()]);
            (-wo).reflect(m)
        } else {
            around(
                Vec3([0., 0., 1.]),
                rng.uniform().sqrt(),
                2. * PI * rng.uniform(),
            )
        };
        let (f, pdf) = self.evaluate(wo, wi);
        if pdf <= 0. || !pdf.is_finite() {
            return None;
        }
        Some((wi, f.map(|v| v * wi.0[2] / pdf), pdf))
    }
}

#[derive(Clone, Copy)]
struct Frame {
    x: Vec3,
    y: Vec3,
    z: Vec3,
}
impl Frame {
    fn new(z: Vec3) -> Self {
        let x = if z.0[2].abs() < 0.99 {
            Vec3([0., 0., 1.]).cross(z).unit()
        } else {
            Vec3([0., 1., 0.]).cross(z).unit()
        };
        Self {
            x,
            y: z.cross(x),
            z,
        }
    }
    fn local(self, v: Vec3) -> Vec3 {
        Vec3([v.dot(self.x), v.dot(self.y), v.dot(self.z)])
    }
    fn world(self, v: Vec3) -> Vec3 {
        self.x * v.0[0] + self.y * v.0[1] + self.z * v.0[2]
    }
}

/// Linear albedo pyramid: sRGB is decoded before averaging or interpolation.
pub struct Texture {
    levels: Vec<(usize, usize, Vec<Rgb>)>,
}
impl Texture {
    pub fn from_png(path: impl AsRef<Path>) -> Result<Self, Box<dyn std::error::Error>> {
        let mut decoder = png::Decoder::new(std::fs::File::open(path)?);
        decoder.set_transformations(png::Transformations::EXPAND | png::Transformations::STRIP_16);
        let mut reader = decoder.read_info()?;
        let mut data = vec![0; reader.output_buffer_size()];
        let info = reader.next_frame(&mut data)?;
        let channels = match info.color_type {
            png::ColorType::Rgb => 3,
            png::ColorType::Rgba => 4,
            _ => return Err("texture must be RGB or RGBA".into()),
        };
        let pixels = data[..info.buffer_size()]
            .chunks_exact(channels)
            .map(|p| {
                std::array::from_fn(|k| {
                    let s = p[k] as f64 / 255.;
                    if s <= 0.04045 {
                        s / 12.92
                    } else {
                        ((s + 0.055) / 1.055).powf(2.4)
                    }
                })
            })
            .collect();
        Ok(Self::from_linear(
            info.width as usize,
            info.height as usize,
            pixels,
        ))
    }
    pub fn from_linear(width: usize, height: usize, pixels: Vec<Rgb>) -> Self {
        assert!(width > 0 && height > 0 && pixels.len() == width * height);
        let mut levels = vec![(width, height, pixels)];
        while levels.last().unwrap().0 > 1 || levels.last().unwrap().1 > 1 {
            let (w, h, p) = levels.last().unwrap();
            let (nw, nh) = (w.div_ceil(2), h.div_ceil(2));
            let mut next = vec![[0.; 3]; nw * nh];
            for y in 0..nh {
                for x in 0..nw {
                    let mut count = 0.;
                    for yy in 2 * y..(2 * y + 2).min(*h) {
                        for xx in 2 * x..(2 * x + 2).min(*w) {
                            add(&mut next[y * nw + x], p[yy * w + xx]);
                            count += 1.;
                        }
                    }
                    next[y * nw + x] = next[y * nw + x].map(|v| v / count);
                }
            }
            levels.push((nw, nh, next));
        }
        Self { levels }
    }
    pub fn sample(&self, u: f64, v: f64, footprint: f64) -> Rgb {
        let level = (footprint.max(0.) * self.levels[0].0.max(self.levels[0].1) as f64)
            .max(1.)
            .log2();
        let a = (level.floor() as usize).min(self.levels.len() - 1);
        let b = (a + 1).min(self.levels.len() - 1);
        let fetch = |i: usize| {
            let (w, h, p) = &self.levels[i];
            let x = u.rem_euclid(1.) * (*w as f64) - 0.5;
            let y = v.rem_euclid(1.) * (*h as f64) - 0.5;
            let (ix, iy) = (x.floor() as isize, y.floor() as isize);
            let (fx, fy) = (x - x.floor(), y - y.floor());
            let mut out = [0.; 3];
            for (dy, wy) in [(0, 1. - fy), (1, fy)] {
                for (dx, wx) in [(0, 1. - fx), (1, fx)] {
                    let xx = (ix + dx).rem_euclid(*w as isize) as usize;
                    let yy = (iy + dy).rem_euclid(*h as isize) as usize;
                    add(&mut out, p[yy * w + xx].map(|c| c * wx * wy));
                }
            }
            out
        };
        let (x, y) = (fetch(a), fetch(b));
        let t = level.fract();
        std::array::from_fn(|k| x[k] * (1. - t) + y[k] * t)
    }
}

pub struct Mesh {
    pub geometry: TriMesh,
    bvh: Bvh,
    normals: Option<Vec<Vec3>>,
    bounds: Bounds,
}
impl Mesh {
    pub fn new(geometry: TriMesh, smooth: bool) -> Self {
        let bvh = Bvh::build(&geometry);
        let (lo, hi) = geometry.bbox();
        let normals = if smooth {
            let mut normals = vec![Vec3([0.; 3]); geometry.vertices.len()];
            for ids in &geometry.triangles {
                let [a, b, c] = ids.map(|i| Vec3(geometry.vertices[i as usize]));
                let n = (b - a).cross(c - a);
                for &i in ids {
                    normals[i as usize] = normals[i as usize] + n;
                }
            }
            Some(
                normals
                    .into_iter()
                    .map(|n| {
                        if n.length() > 0. {
                            n.unit()
                        } else {
                            Vec3([0., 0., 1.])
                        }
                    })
                    .collect(),
            )
        } else {
            None
        };
        Self {
            geometry,
            bvh,
            normals,
            bounds: Bounds {
                lo: Vec3(lo),
                hi: Vec3(hi),
            },
        }
    }
}

pub enum Shape {
    Mesh(Mesh),
    Sphere { centre: Vec3, radius: f64 },
    Quad { corner: Vec3, u: Vec3, v: Vec3 },
}
#[derive(Clone, Copy)]
pub struct Hit {
    pub distance: f64,
    pub position: Vec3,
    pub geometric: Vec3,
    pub normal: Vec3,
    pub object: usize,
}
impl Shape {
    fn bounds(&self) -> Bounds {
        match self {
            Self::Mesh(m) => m.bounds,
            Self::Sphere { centre, radius } => Bounds {
                lo: *centre - Vec3([*radius; 3]),
                hi: *centre + Vec3([*radius; 3]),
            },
            Self::Quad { corner, u, v } => {
                let p = [*corner, *corner + *u, *corner + *v, *corner + *u + *v];
                Bounds {
                    lo: Vec3(std::array::from_fn(|k| {
                        p.iter().map(|p| p.0[k]).fold(f64::INFINITY, f64::min) - 1e-12
                    })),
                    hi: Vec3(std::array::from_fn(|k| {
                        p.iter().map(|p| p.0[k]).fold(f64::NEG_INFINITY, f64::max) + 1e-12
                    })),
                }
            }
        }
    }
    fn intersect(&self, ray: Ray, min: f64, max: f64) -> Option<Hit> {
        let (distance, geometric, normal) = match self {
            Self::Mesh(m) => {
                let h = m.bvh.ray(&m.geometry, ray.origin.0, ray.dir.0, min, max)?;
                let geom = Vec3(m.geometry.unit_normal(h.triangle as usize));
                let normal = if let Some(normals) = &m.normals {
                    let ids = m.geometry.triangles[h.triangle as usize];
                    (normals[ids[0] as usize] * h.barycentric[0]
                        + normals[ids[1] as usize] * h.barycentric[1]
                        + normals[ids[2] as usize] * h.barycentric[2])
                        .unit()
                } else {
                    geom
                };
                (h.distance, geom, normal)
            }
            Self::Sphere { centre, radius } => {
                let oc = ray.origin - *centre;
                let b = oc.dot(ray.dir);
                let disc = b * b - oc.dot(oc) + radius * radius;
                if disc < 0. {
                    return None;
                }
                let root = disc.sqrt();
                let mut t = -b - root;
                if t <= min {
                    t = -b + root;
                }
                if t <= min || t >= max {
                    return None;
                }
                let n = (ray.at(t) - *centre) / *radius;
                (t, n, n)
            }
            Self::Quad { corner, u, v } => {
                let n = u.cross(*v).unit();
                let denominator = ray.dir.dot(n);
                if denominator == 0. {
                    return None;
                }
                let t = (*corner - ray.origin).dot(n) / denominator;
                if t <= min || t >= max {
                    return None;
                }
                let p = ray.at(t) - *corner;
                let (uu, uv, vv) = (u.dot(*u), u.dot(*v), v.dot(*v));
                let det = uu * vv - uv * uv;
                let a = (p.dot(*u) * vv - p.dot(*v) * uv) / det;
                let b = (p.dot(*v) * uu - p.dot(*u) * uv) / det;
                if !(0.0..=1.0).contains(&a) || !(0.0..=1.0).contains(&b) {
                    return None;
                }
                (t, n, n)
            }
        };
        Some(Hit {
            distance,
            position: ray.at(distance),
            geometric,
            normal,
            object: 0,
        })
    }
}
type Shader = Arc<dyn Fn(Vec3, Vec3, f64) -> (Bsdf, Vec3) + Send + Sync>;
type Accept = Arc<dyn Fn(Vec3) -> bool + Send + Sync>;
pub struct Object {
    pub shape: Shape,
    pub bsdf: Bsdf,
    pub shader: Option<Shader>,
    pub accept: Option<Accept>,
    light: Option<usize>,
}
impl Object {
    pub fn new(shape: Shape, bsdf: Bsdf) -> Self {
        Self {
            shape,
            bsdf,
            shader: None,
            accept: None,
            light: None,
        }
    }
    pub fn shader(
        mut self,
        f: impl Fn(Vec3, Vec3, f64) -> (Bsdf, Vec3) + Send + Sync + 'static,
    ) -> Self {
        self.shader = Some(Arc::new(f));
        self
    }
    pub fn accept(mut self, f: impl Fn(Vec3) -> bool + Send + Sync + 'static) -> Self {
        self.accept = Some(Arc::new(f));
        self
    }
    fn intersect(&self, ray: Ray, mut min: f64, max: f64) -> Option<Hit> {
        // Rejected interfaces can occur at explicit liquid/bubble unions.
        loop {
            let h = self.shape.intersect(ray, min, max)?;
            if self.accept.as_ref().is_none_or(|f| f(h.position)) {
                return Some(h);
            }
            min = h.distance + 1e-10;
        }
    }
}
#[derive(Clone, Copy)]
pub struct AreaLight {
    pub corner: Vec3,
    pub u: Vec3,
    pub v: Vec3,
    pub radiance: Rgb,
}
impl AreaLight {
    fn pdf(self, from: Vec3, to: Vec3) -> f64 {
        let d = to - from;
        let n = self.u.cross(self.v);
        let c = n.unit().dot(-d.unit()).max(0.);
        if c == 0. {
            0.
        } else {
            d.dot(d) / (c * n.length())
        }
    }
}
struct Node {
    bounds: Bounds,
    children: Option<[usize; 2]>,
    indices: Vec<usize>,
}
pub struct World {
    pub objects: Vec<Object>,
    pub lights: Vec<AreaLight>,
    pub environment: Rgb,
    nodes: Vec<Node>,
}
impl World {
    pub fn new(mut objects: Vec<Object>, lights: Vec<AreaLight>, environment: Rgb) -> Self {
        for (i, l) in lights.iter().enumerate() {
            let mut o = Object::new(
                Shape::Quad {
                    corner: l.corner,
                    u: l.u,
                    v: l.v,
                },
                Bsdf::Emitter {
                    radiance: l.radiance,
                },
            );
            o.light = Some(i);
            objects.push(o);
        }
        let mut world = Self {
            objects,
            lights,
            environment,
            nodes: Vec::new(),
        };
        if !world.objects.is_empty() {
            world.build((0..world.objects.len()).collect());
        }
        world
    }
    fn build(&mut self, mut indices: Vec<usize>) -> usize {
        let mut lo = Vec3([f64::INFINITY; 3]);
        let mut hi = Vec3([f64::NEG_INFINITY; 3]);
        for &i in &indices {
            let b = self.objects[i].shape.bounds();
            for k in 0..3 {
                lo.0[k] = lo.0[k].min(b.lo.0[k]);
                hi.0[k] = hi.0[k].max(b.hi.0[k]);
            }
        }
        let id = self.nodes.len();
        self.nodes.push(Node {
            bounds: Bounds { lo, hi },
            children: None,
            indices: Vec::new(),
        });
        if indices.len() <= 4 {
            self.nodes[id].indices = indices;
        } else {
            let axis = (0..3)
                .max_by(|&a, &b| (hi.0[a] - lo.0[a]).total_cmp(&(hi.0[b] - lo.0[b])))
                .unwrap();
            indices.sort_unstable_by(|&a, &b| {
                let x = self.objects[a].shape.bounds();
                let y = self.objects[b].shape.bounds();
                (x.lo.0[axis] + x.hi.0[axis]).total_cmp(&(y.lo.0[axis] + y.hi.0[axis]))
            });
            let right = indices.split_off(indices.len() / 2);
            let left = self.build(indices);
            let right = self.build(right);
            self.nodes[id].children = Some([left, right]);
        }
        id
    }
    pub fn intersect(&self, ray: Ray, max: f64, skip_lights: bool) -> Option<Hit> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut stack = [0_usize; 64];
        let mut len = 1;
        let mut best = None;
        let mut limit = max;
        while len > 0 {
            len -= 1;
            let node = &self.nodes[stack[len]];
            if node.bounds.interval(ray, limit).is_none() {
                continue;
            }
            if let Some([a, b]) = node.children {
                stack[len] = a;
                stack[len + 1] = b;
                len += 2;
            } else {
                for &i in &node.indices {
                    if skip_lights && self.objects[i].light.is_some() {
                        continue;
                    }
                    if let Some(mut h) = self.objects[i].intersect(ray, 1e-8, limit) {
                        limit = h.distance;
                        h.object = i;
                        best = Some(h);
                    }
                }
            }
        }
        best
    }
    pub fn radiance(&self, ray: Ray, rng: &mut Rng, max_depth: usize) -> Rgb {
        self.radiance_with_cone(ray, rng, max_depth, 0.0003)
    }
    pub fn radiance_with_cone(
        &self,
        mut ray: Ray,
        rng: &mut Rng,
        max_depth: usize,
        cone: f64,
    ) -> Rgb {
        let mut result = [0.; 3];
        let mut beta = [1.; 3];
        let mut medium = Medium::AIR;
        let mut previous: Option<(Vec3, f64)> = None;
        for depth in 0..max_depth {
            let Some(hit) = self.intersect(ray, f64::INFINITY, false) else {
                add(&mut result, mul(beta, self.environment));
                break;
            };
            beta = mul(beta, medium.transmission(hit.distance));
            let object = &self.objects[hit.object];
            let (bsdf, mut n) = object
                .shader
                .as_ref()
                .map_or((object.bsdf, hit.normal), |f| {
                    f(hit.position, hit.normal, hit.distance * cone)
                });
            let front = ray.dir.dot(hit.geometric) < 0.;
            let ng = if front { hit.geometric } else { -hit.geometric };
            if n.dot(ng) < 0. {
                n = -n;
            }
            let wo = -ray.dir;
            if n.dot(wo) <= 0. {
                n = ng;
            }
            if let Bsdf::Emitter { radiance } = bsdf {
                if front {
                    let weight = if let (Some((p, pdf)), Some(i)) = (previous, object.light) {
                        power(
                            pdf,
                            self.lights[i].pdf(p, hit.position) / self.lights.len() as f64,
                        )
                    } else {
                        1.
                    };
                    add(&mut result, mul(beta, radiance.map(|v| v * weight)));
                }
                break;
            }
            if let Bsdf::Dielectric { outside, inside } = bsdf {
                let (ni, nt, next) = if front {
                    (outside.index, inside.index, inside)
                } else {
                    (inside.index, outside.index, outside)
                };
                let r = fresnel_dielectric(wo.dot(ng).abs(), ni, nt);
                let reflected = Ray::new(hit.position + ng * 2e-7, ray.dir.reflect(ng));
                let terminal = if depth == 0 {
                    match self.intersect(reflected, f64::INFINITY, false) {
                        None => Some(self.environment),
                        Some(h) => match self.objects[h.object].bsdf {
                            Bsdf::Emitter { radiance }
                                if self.objects[h.object].shader.is_none() =>
                            {
                                Some(if reflected.dir.dot(h.geometric) < 0. {
                                    radiance
                                } else {
                                    [0.; 3]
                                })
                            }
                            _ => None,
                        },
                    }
                } else {
                    None
                };
                // Rao–Blackwellize the first reflected branch when it has an
                // exactly known terminal value. This is variance reduction,
                // not a painted highlight; indirect reflections still trace.
                if let Some(value) = terminal {
                    add(&mut result, mul(beta, value.map(|v| v * r)));
                    if r >= 1. {
                        break;
                    }
                }
                let direction = if terminal.is_none() && rng.uniform() < r {
                    ray.dir.reflect(ng)
                } else if let Some(dir) = refract(ray.dir, ng, ni, nt) {
                    medium = next;
                    beta = beta.map(|v| {
                        v * (ni / nt).powi(2) * if terminal.is_some() { 1. - r } else { 1. }
                    });
                    dir
                } else {
                    if terminal.is_some() {
                        break;
                    }
                    ray.dir.reflect(ng)
                };
                let sign = if direction.dot(hit.geometric) > 0. {
                    1.
                } else {
                    -1.
                };
                ray = Ray::new(hit.position + hit.geometric * (sign * 2e-7), direction);
                previous = None;
                continue;
            }
            let frame = Frame::new(n);
            let local_wo = frame.local(wo);
            if !self.lights.is_empty() {
                let index = ((rng.uniform() * self.lights.len() as f64) as usize)
                    .min(self.lights.len() - 1);
                let light = self.lights[index];
                let p = light.corner + light.u * rng.uniform() + light.v * rng.uniform();
                let delta = p - hit.position;
                let wi = delta.unit();
                let pdf = light.pdf(hit.position, p) / self.lights.len() as f64;
                if pdf > 0. && wi.dot(ng) > 0. {
                    let (f, bpdf) = bsdf.evaluate(local_wo, frame.local(wi));
                    let shadow = Ray::new(hit.position + ng * 2e-7, wi);
                    if self
                        .intersect(shadow, delta.length() - 4e-7, false)
                        .is_none()
                    {
                        let weight = n.dot(wi).max(0.) * power(pdf, bpdf) / pdf;
                        add(
                            &mut result,
                            mul(beta, mul(f, light.radiance.map(|v| v * weight))),
                        );
                    }
                }
            }
            let Some((wi, weight, pdf)) = bsdf.sample(local_wo, rng) else {
                break;
            };
            let wi = frame.world(wi);
            if wi.dot(ng) <= 0. {
                break;
            }
            beta = mul(beta, weight);
            previous = Some((hit.position, pdf));
            ray = Ray::new(hit.position + ng * 2e-7, wi);
            if depth >= 4 {
                let survive = beta.iter().copied().fold(0., f64::max).clamp(0.05, 0.95);
                if rng.uniform() > survive {
                    break;
                }
                beta = beta.map(|v| v / survive);
            }
        }
        result
    }
}
fn power(a: f64, b: f64) -> f64 {
    let a = a * a;
    let b = b * b;
    a / (a + b)
}

/// Deterministic range of rows, enabling bounded, restartable render jobs.
pub fn render_rows(
    world: &World,
    camera: &Camera,
    width: usize,
    height: usize,
    samples: usize,
    start: usize,
    end: usize,
) -> Image {
    assert!(width > 0 && height > 0 && samples > 0 && start < end && end <= height);
    let mut image = Image::new(width, end - start);
    let done = AtomicUsize::new(0);
    image
        .pixels
        .par_chunks_mut(width)
        .enumerate()
        .for_each(|(row, pixels)| {
            let y = start + row;
            for (x, pixel) in pixels.iter_mut().enumerate() {
                let mut sum = [0.; 3];
                for s in 0..samples {
                    let mut hash = ((y * width + x) as u64).wrapping_mul(0x9e3779b97f4a7c15)
                        ^ (s as u64).wrapping_mul(0xbf58476d1ce4e5b9)
                        ^ 0x87f11b982bdd7a35;
                    hash = (hash ^ (hash >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
                    hash = (hash ^ (hash >> 27)).wrapping_mul(0x94d049bb133111eb);
                    let mut rng = Rng(hash ^ (hash >> 31));
                    let [u, v] = camera.pixel(x, y, width, height, [rng.uniform(), rng.uniform()]);
                    let ray = camera.ray(u, v, width as f64 / height as f64, &mut rng);
                    add(
                        &mut sum,
                        world.radiance_with_cone(
                            ray,
                            &mut rng,
                            24,
                            2. * (camera.fov_y / 2.).tan() / height as f64,
                        ),
                    );
                }
                *pixel = sum.map(|v| v / samples as f64);
            }
            let n = done.fetch_add(1, Ordering::Relaxed) + 1;
            if n.is_multiple_of(32) {
                eprintln!("{n}/{} rows", end - start);
            }
        });
    image
}

pub fn read_pfm(path: impl AsRef<Path>) -> std::io::Result<Image> {
    let data = std::fs::read(path)?;
    let mut lines = data.splitn(4, |b| *b == b'\n');
    if lines.next() != Some(b"PF".as_slice()) {
        return Err(std::io::Error::other("expected RGB PFM"));
    }
    let dims = std::str::from_utf8(
        lines
            .next()
            .ok_or_else(|| std::io::Error::other("PFM dimensions"))?,
    )
    .map_err(std::io::Error::other)?;
    let sizes: Vec<usize> = dims
        .split_whitespace()
        .map(str::parse)
        .collect::<Result<_, _>>()
        .map_err(std::io::Error::other)?;
    if sizes.len() != 2 || sizes.contains(&0) || lines.next() != Some(b"-1.0".as_slice()) {
        return Err(std::io::Error::other("unsupported PFM"));
    }
    let (w, h) = (sizes[0], sizes[1]);
    let bytes = lines
        .next()
        .ok_or_else(|| std::io::Error::other("PFM pixels"))?;
    if bytes.len()
        != w.checked_mul(h)
            .and_then(|n| n.checked_mul(12))
            .ok_or_else(|| std::io::Error::other("PFM overflow"))?
    {
        return Err(std::io::Error::other("PFM size"));
    }
    let mut image = Image::new(w, h);
    for (i, p) in bytes.chunks_exact(12).enumerate() {
        image.pixels[(h - 1 - i / w) * w + i % w] = std::array::from_fn(|k| {
            f32::from_le_bytes(p[4 * k..4 * k + 4].try_into().unwrap()) as f64
        });
    }
    Ok(image)
}

/// Conventional spatial display filter. Divide out the actual material albedo
/// before filtering illumination, then restore it, so wood/fibre texture is
/// retained. Geometry/normal guides keep edges and reflections from bleeding
/// between objects. This is biased postprocessing; retain unfiltered radiance.
pub fn denoise(world: &World, camera: &Camera, image: &Image) -> Image {
    #[derive(Clone, Copy)]
    struct Guide {
        object: usize,
        normal: Vec3,
        albedo: Rgb,
    }
    let (w, h) = (image.width, image.height);
    let mut pinhole = *camera;
    pinhole.aperture = 0.;
    let guides: Vec<Guide> = (0..w * h)
        .into_par_iter()
        .map(|i| {
            let mut rng = Rng(1);
            let [u, v] = pinhole.pixel(i % w, i / w, w, h, [0.5; 2]);
            let ray = pinhole.ray(u, v, w as f64 / h as f64, &mut rng);
            if let Some(hit) = world.intersect(ray, f64::INFINITY, false) {
                let o = &world.objects[hit.object];
                let (bsdf, normal) = o.shader.as_ref().map_or((o.bsdf, hit.normal), |f| {
                    f(
                        hit.position,
                        hit.normal,
                        hit.distance * 2. * (camera.fov_y / 2.).tan() / h as f64,
                    )
                });
                let albedo = match bsdf {
                    Bsdf::Diffuse { albedo } | Bsdf::Coated { albedo, .. } => {
                        albedo.map(|v| v.max(0.08))
                    }
                    _ => [1.; 3],
                };
                Guide {
                    object: hit.object,
                    normal,
                    albedo,
                }
            } else {
                Guide {
                    object: usize::MAX,
                    normal: Vec3([0., 0., 1.]),
                    albedo: [1.; 3],
                }
            }
        })
        .collect();
    let mut light: Vec<Rgb> = image
        .pixels
        .iter()
        .zip(&guides)
        .map(|(p, g)| std::array::from_fn(|k| p[k] / g.albedo[k]))
        .collect();
    for step in [1_isize, 2, 4] {
        light = (0..w * h)
            .into_par_iter()
            .map(|i| {
                let centre = light[i];
                let g = guides[i];
                let mut sum = [0.; 3];
                let mut total = 0.;
                let sigma = 0.018 + 0.18 * luminance(centre);
                for dy in -2_isize..=2 {
                    for dx in -2_isize..=2 {
                        let x = (i % w) as isize + dx * step;
                        let y = (i / w) as isize + dy * step;
                        if x < 0 || y < 0 || x >= w as isize || y >= h as isize {
                            continue;
                        }
                        let j = y as usize * w + x as usize;
                        let other = guides[j];
                        if g.object != other.object {
                            continue;
                        }
                        let diff = luminance(light[j]) - luminance(centre);
                        let weight = (-0.5 * (dx * dx + dy * dy) as f64
                            - 48. * (1. - g.normal.dot(other.normal).clamp(0., 1.))
                            - diff * diff / (2. * sigma * sigma))
                            .exp();
                        add(&mut sum, light[j].map(|v| v * weight));
                        total += weight;
                    }
                }
                sum.map(|v| v / total)
            })
            .collect();
    }
    Image {
        width: w,
        height: h,
        pixels: light
            .into_iter()
            .zip(guides)
            .map(|(p, g)| mul(p, g.albedo))
            .collect(),
    }
}
