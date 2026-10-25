//! SI isolated gravity: spherical shells O(n) build/O(log n) query, octree
//! quadrupoles and a direct reference. No periodic image masses or cosmology.
pub use radiation::math::Vec3;
pub mod metric;
use rayon::prelude::*;
pub const G: f64 = 6.674_30e-11;
pub const C: f64 = radiation::C;
pub fn norm(v: Vec3) -> f64 {
    let squared = v.dot(v);
    if squared.is_finite() && squared >= f64::MIN_POSITIVE {
        return squared.sqrt();
    }
    v.0[0].hypot(v.0[1]).hypot(v.0[2])
}
fn finite(v: Vec3) -> bool {
    v.0.iter().all(|x| x.is_finite())
}
#[derive(Clone, Copy, Debug, Default)]
pub struct Sample {
    pub potential: f64,
    pub acceleration: Vec3,
}
impl Sample {
    pub fn validate(self) -> Result<Self, &'static str> {
        if self.potential.is_finite() && finite(self.acceleration) {
            Ok(self)
        } else {
            Err("unrepresentable gravitational field")
        }
    }
}
pub trait Field: Sync {
    fn sample(&self, position: Vec3) -> Result<Sample, &'static str>;
}
#[derive(Clone, Copy, Debug)]
pub struct PointMass {
    pub gm: f64,
    pub center: Vec3,
    pub softening: f64,
}
impl Field for PointMass {
    fn sample(&self, x: Vec3) -> Result<Sample, &'static str> {
        if !finite(x)
            || !finite(self.center)
            || !self.gm.is_finite()
            || self.gm < 0.0
            || !self.softening.is_finite()
            || self.softening < 0.0
        {
            return Err("invalid point gravity source");
        }
        point(self.center - x, self.gm, self.softening)
    }
}
fn point(d: Vec3, gm: f64, epsilon: f64) -> Result<Sample, &'static str> {
    if gm == 0.0 {
        return Ok(Sample::default());
    }
    let squared = d.dot(d) + epsilon * epsilon;
    let r = if squared.is_finite() && squared >= f64::MIN_POSITIVE {
        squared.sqrt()
    } else {
        norm(d).hypot(epsilon)
    };
    if r == 0.0 {
        return Err("unsoftened coincident gravity source");
    }
    Sample {
        potential: -gm / r,
        acceleration: (d / r) * (gm / r / r),
    }
    .validate()
}
/// Piecewise constant spherical mass density, with exact shell integrals.
#[derive(Clone, Debug)]
pub struct Spherical {
    edges: Vec<f64>,
    density: Vec<f64>,
    mass: Vec<f64>,
    outer: Vec<f64>,
    center: Vec3,
}
impl Spherical {
    pub fn new(edges: Vec<f64>, density: Vec<f64>, center: Vec3) -> Result<Self, &'static str> {
        if density.is_empty()
            || edges.len() != density.len() + 1
            || edges[0] != 0.0
            || !finite(center)
            || edges.iter().any(|x| !x.is_finite())
            || edges.windows(2).any(|x| x[1] <= x[0])
            || density.iter().any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err("invalid spherical gravity grid");
        }
        let n = density.len();
        let mut mass = vec![0.0; n + 1];
        let mut outer = vec![0.0; n + 1];
        for i in 0..n {
            let (a, b) = (edges[i], edges[i + 1]);
            mass[i + 1] = mass[i]
                + 4.0 * std::f64::consts::PI / 3.0 * density[i] * (b - a) * (b * b + a * b + a * a);
        }
        for i in (0..n).rev() {
            let (a, b) = (edges[i], edges[i + 1]);
            outer[i] = outer[i + 1] + 2.0 * std::f64::consts::PI * density[i] * (b - a) * (b + a);
        }
        if mass.iter().chain(&outer).any(|x| !x.is_finite()) {
            return Err("unrepresentable spherical source mass");
        }
        Ok(Self {
            edges,
            density,
            mass,
            outer,
            center,
        })
    }
    pub fn total_mass(&self) -> f64 {
        *self.mass.last().unwrap()
    }
    pub fn radius(&self) -> f64 {
        *self.edges.last().unwrap()
    }
    pub fn enclosed_mass(&self, r: f64) -> Result<f64, &'static str> {
        if !r.is_finite() || r < 0.0 {
            return Err("invalid spherical radius");
        }
        if r >= self.radius() {
            return Ok(self.total_mass());
        }
        let i = self.edges.partition_point(|x| *x <= r) - 1;
        let a = self.edges[i];
        Ok(self.mass[i]
            + 4.0 * std::f64::consts::PI / 3.0
                * self.density[i]
                * (r - a)
                * (r * r + r * a + a * a))
    }
}
impl Field for Spherical {
    fn sample(&self, x: Vec3) -> Result<Sample, &'static str> {
        if !finite(x) {
            return Err("invalid gravity query position");
        }
        let d = x - self.center;
        let r = norm(d);
        if r == 0.0 {
            return Sample {
                potential: -G * self.outer[0],
                acceleration: Vec3::default(),
            }
            .validate();
        }
        if r >= self.radius() {
            return point(-d, G * self.total_mass(), 0.0);
        }
        let i = self.edges.partition_point(|x| *x <= r) - 1;
        let mass = self.enclosed_mass(r)?;
        let b = self.edges[i + 1];
        let exterior =
            self.outer[i + 1] + 2.0 * std::f64::consts::PI * self.density[i] * (b - r) * (b + r);
        Sample {
            potential: -G * (mass / r + exterior),
            acceleration: -(d / r) * (G * mass / r / r),
        }
        .validate()
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Body {
    pub position: Vec3,
    pub mass: f64,
}
fn validate_bodies(b: &[Body], eps: f64) -> Result<(), &'static str> {
    if !eps.is_finite()
        || eps < 0.0
        || b.iter()
            .any(|b| !finite(b.position) || !b.mass.is_finite() || b.mass < 0.0)
    {
        Err("invalid gravitational bodies/softening")
    } else {
        Ok(())
    }
}
/// Plummer-softened reference. `exclude` skips a source index, not all bodies
/// at the same position. Softening is a resolved-scale model input in meters.
pub fn direct(
    bodies: &[Body],
    x: Vec3,
    epsilon: f64,
    exclude: Option<usize>,
) -> Result<Sample, &'static str> {
    validate_bodies(bodies, epsilon)?;
    direct_validated(bodies, x, epsilon, exclude)
}
fn direct_validated(
    bodies: &[Body],
    x: Vec3,
    epsilon: f64,
    exclude: Option<usize>,
) -> Result<Sample, &'static str> {
    if !finite(x) || exclude.is_some_and(|i| i >= bodies.len()) {
        return Err("invalid gravity query");
    }
    let mut out = Sample::default();
    for (i, b) in bodies.iter().enumerate() {
        if Some(i) == exclude {
            continue;
        }
        let f = point(b.position - x, G * b.mass, epsilon)?;
        out.potential += f.potential;
        out.acceleration = out.acceleration + f.acceleration;
    }
    out.validate()
}
/// Validates source data once; useful for fair direct/tree accuracy and timing
/// comparisons, and small source populations where a tree costs more.
pub struct Direct {
    bodies: Vec<Body>,
    epsilon: f64,
}
impl Direct {
    pub fn new(bodies: Vec<Body>, epsilon: f64) -> Result<Self, &'static str> {
        validate_bodies(&bodies, epsilon)?;
        Ok(Self { bodies, epsilon })
    }
    pub fn evaluate(&self, x: Vec3, exclude: Option<usize>) -> Result<Sample, &'static str> {
        if exclude.is_some_and(|i| i >= self.bodies.len()) {
            return Err("invalid excluded gravity source");
        }
        direct_validated(&self.bodies, x, self.epsilon, exclude)
    }
}
impl Field for Direct {
    fn sample(&self, x: Vec3) -> Result<Sample, &'static str> {
        self.evaluate(x, None)
    }
}
#[derive(Debug)]
struct Node {
    center: Vec3,
    half: f64,
    com: Vec3,
    mass: f64,
    moment: [[f64; 3]; 3],
    children: Vec<usize>,
    indices: Vec<usize>,
}
#[derive(Clone, Copy, Debug, Default)]
pub struct Work {
    pub nodes: usize,
    pub direct_pairs: usize,
    pub multipoles: usize,
}
pub struct Tree {
    bodies: Vec<Body>,
    nodes: Vec<Node>,
    epsilon: f64,
    theta: f64,
}
impl Tree {
    /// Quadrupole tree. Opening compares bounding support radius to distance,
    /// with COM displacement included. theta=0 gives the direct reference.
    pub fn new(bodies: Vec<Body>, epsilon: f64, theta: f64) -> Result<Self, &'static str> {
        validate_bodies(&bodies, epsilon)?;
        if !theta.is_finite() || !(0.0..=0.7).contains(&theta) {
            return Err("invalid gravity opening angle");
        }
        let mut s = Self {
            bodies,
            nodes: vec![],
            epsilon,
            theta,
        };
        if s.bodies.is_empty() {
            return Ok(s);
        }
        let mut lo = s.bodies[0].position;
        let mut hi = lo;
        for b in &s.bodies {
            for d in 0..3 {
                lo.0[d] = lo.0[d].min(b.position.0[d]);
                hi.0[d] = hi.0[d].max(b.position.0[d]);
            }
        }
        let center = (lo + hi) * 0.5;
        let half = (0..3)
            .map(|d| (hi.0[d] - lo.0[d]) * 0.5)
            .fold(0.0, f64::max)
            * (1.0 + 1e-12);
        s.build((0..s.bodies.len()).collect(), center, half, 0);
        if s.nodes.iter().any(|n| {
            !n.mass.is_finite()
                || !finite(n.com)
                || !n.half.is_finite()
                || n.moment.iter().flatten().any(|x| !x.is_finite())
        }) {
            return Err("unrepresentable tree moments");
        }
        Ok(s)
    }
    fn build(&mut self, indices: Vec<usize>, center: Vec3, half: f64, depth: usize) -> usize {
        let mass: f64 = indices.iter().map(|i| self.bodies[*i].mass).sum();
        let mut offset = Vec3::default();
        if mass > 0.0 {
            for i in &indices {
                let b = self.bodies[*i];
                offset = offset + (b.position - center) * (b.mass / mass);
            }
        }
        let com = center + offset;
        let mut moment = [[0.0; 3]; 3];
        for i in &indices {
            let b = self.bodies[*i];
            let d = b.position - com;
            for (a, row) in moment.iter_mut().enumerate() {
                for (c, value) in row.iter_mut().enumerate() {
                    *value += b.mass * d.0[a] * d.0[c];
                }
            }
        }
        let idx = self.nodes.len();
        self.nodes.push(Node {
            center,
            half,
            com,
            mass,
            moment,
            children: vec![],
            indices: vec![],
        });
        if indices.len() <= 8 || depth >= 48 || half == 0.0 {
            self.nodes[idx].indices = indices;
            return idx;
        }
        let mut groups: [Vec<usize>; 8] = std::array::from_fn(|_| vec![]);
        for i in indices {
            let p = self.bodies[i].position;
            let oct = (0..3).fold(0, |o, d| o | ((p.0[d] >= center.0[d]) as usize) << d);
            groups[oct].push(i);
        }
        for (oct, g) in groups.into_iter().enumerate() {
            if g.is_empty() {
                continue;
            }
            let next = Vec3(std::array::from_fn(|d| {
                center.0[d]
                    + if oct & (1 << d) != 0 {
                        half * 0.5
                    } else {
                        -half * 0.5
                    }
            }));
            let child = self.build(g, next, half * 0.5, depth + 1);
            self.nodes[idx].children.push(child);
        }
        idx
    }
    pub fn evaluate(
        &self,
        x: Vec3,
        exclude: Option<usize>,
    ) -> Result<(Sample, Work), &'static str> {
        if !finite(x) || exclude.is_some_and(|i| i >= self.bodies.len()) {
            return Err("invalid gravity tree query");
        }
        let mut out = Sample::default();
        let mut work = Work::default();
        if !self.nodes.is_empty() {
            self.walk(0, x, exclude, &mut out, &mut work)?;
        }
        Ok((out.validate()?, work))
    }
    fn walk(
        &self,
        idx: usize,
        x: Vec3,
        exclude: Option<usize>,
        out: &mut Sample,
        work: &mut Work,
    ) -> Result<(), &'static str> {
        work.nodes += 1;
        let n = &self.nodes[idx];
        if n.mass == 0.0 {
            return Ok(());
        }
        let d = n.com - x;
        let distance = norm(d);
        let inside = (0..3).all(|a| (x.0[a] - n.center.0[a]).abs() <= n.half);
        let excluded_inside = exclude.is_some_and(|i| {
            (0..3).all(|a| (self.bodies[i].position.0[a] - n.center.0[a]).abs() <= n.half)
        });
        let support = 3f64.sqrt() * n.half + norm(n.com - n.center);
        if !n.children.is_empty() && !inside && !excluded_inside && support < self.theta * distance
        {
            let q = distance.hypot(self.epsilon);
            let inv = 1.0 / q;
            let unit = d * inv;
            let sd = Vec3(std::array::from_fn(|i| {
                (0..3).map(|j| n.moment[i][j] * unit.0[j]).sum()
            }));
            let a = unit.dot(sd);
            let trace = (0..3).map(|i| n.moment[i][i]).sum::<f64>();
            out.potential -= G * inv * (n.mass + 0.5 * inv * inv * (3.0 * a - trace));
            out.acceleration = out.acceleration
                + (unit * n.mass + (unit * (7.5 * a - 1.5 * trace) - sd * 3.0) * (inv * inv))
                    * (G * inv * inv);
            work.multipoles += 1;
        } else if n.children.is_empty() {
            for i in &n.indices {
                if Some(*i) == exclude {
                    continue;
                }
                let b = self.bodies[*i];
                let s = point(b.position - x, G * b.mass, self.epsilon)?;
                out.potential += s.potential;
                out.acceleration = out.acceleration + s.acceleration;
                work.direct_pairs += 1;
            }
        } else {
            for child in &n.children {
                self.walk(*child, x, exclude, out, work)?;
            }
        }
        Ok(())
    }
    pub fn at_bodies(&self) -> Result<Vec<Sample>, &'static str> {
        let f = |i: usize| self.evaluate(self.bodies[i].position, Some(i)).map(|v| v.0);
        if self.bodies.len() >= 1024 {
            (0..self.bodies.len()).into_par_iter().map(f).collect()
        } else {
            (0..self.bodies.len()).map(f).collect()
        }
    }
}
impl Field for Tree {
    fn sample(&self, x: Vec3) -> Result<Sample, &'static str> {
        self.evaluate(x, None).map(|v| v.0)
    }
}
/// Newtonian bulk/macro bodies. KDK leapfrog rebuilds gravity at drifted
/// positions. Particle tree forces are approximate; momentum diagnostics must
/// be checked against theta=0 when cancellation accuracy matters.
pub struct NBody {
    pub bodies: Vec<Body>,
    pub velocities: Vec<Vec3>,
    pub time: f64,
    pub softening: f64,
    pub theta: f64,
}
impl NBody {
    pub fn step(&mut self, dt: f64) -> Result<(), &'static str> {
        if !dt.is_finite()
            || dt < 0.0
            || !(self.time + dt).is_finite()
            || self.velocities.len() != self.bodies.len()
            || self.velocities.iter().any(|v| !finite(*v))
        {
            return Err("invalid gravitational time step");
        }
        let field = Tree::new(self.bodies.clone(), self.softening, self.theta)?;
        let start = field.at_bodies()?;
        let mut bodies = self.bodies.clone();
        let mut velocities = self.velocities.clone();
        for i in 0..bodies.len() {
            velocities[i] = velocities[i] + start[i].acceleration * (0.5 * dt);
            bodies[i].position = bodies[i].position + velocities[i] * dt;
        }
        let end = Tree::new(bodies.clone(), self.softening, self.theta)?.at_bodies()?;
        for i in 0..bodies.len() {
            velocities[i] = velocities[i] + end[i].acceleration * (0.5 * dt);
            if !finite(velocities[i]) {
                return Err("unrepresentable body velocity");
            }
        }
        self.bodies = bodies;
        self.velocities = velocities;
        self.time += dt;
        Ok(())
    }
    /// Direct potential/kinetic energy [J] and linear momentum [kg m/s].
    pub fn invariants(&self) -> Result<(f64, Vec3), &'static str> {
        validate_bodies(&self.bodies, self.softening)?;
        if self.velocities.len() != self.bodies.len() || self.velocities.iter().any(|v| !finite(*v))
        {
            return Err("invalid body velocities");
        }
        let mut energy = 0.0;
        let mut momentum = Vec3::default();
        for (i, b) in self.bodies.iter().enumerate() {
            let v = self.velocities[i];
            energy += 0.5
                * b.mass
                * (v.dot(v) + direct(&self.bodies, b.position, self.softening, Some(i))?.potential);
            momentum = momentum + v * b.mass;
        }
        if !energy.is_finite() || !finite(momentum) {
            return Err("unrepresentable gravitational invariants");
        }
        Ok((energy, momentum))
    }
}
