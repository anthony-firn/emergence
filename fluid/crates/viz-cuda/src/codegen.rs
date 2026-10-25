//! Compiling a `cad::Solid` into a CUDA device function.
//!
//! The CAD kernel is a Rust enum tree with boxes and reference counts;
//! none of that can cross to a device. Two ways out: interpret a flattened
//! instruction tape on the GPU, or emit the field as CUDA source and let
//! nvrtc compile it. This takes the second, because the generated code is
//! *the same arithmetic* as the Rust — each node prints the formula from
//! `cad::solid`, in the same order, with the same intermediate
//! expressions — so the two agree to the last bit in `double`, and the
//! tests check exactly that. An interpreter would be as correct but slower
//! and much harder to argue about.
//!
//! Each node becomes its own `__device__` function, so loops (finite
//! arrays, polar folds, periodic lattices) come out as loops rather than
//! unrolled expressions, and nvrtc inlines what is worth inlining.
//!
//! # What cannot be compiled
//!
//! [`Solid::Mesh`] and [`Solid::Custom`] carry a triangle BVH and an
//! arbitrary Rust closure. Neither has a CUDA form, and both are reported
//! as errors naming the node rather than silently approximated. A model
//! that needs them renders on the CPU.

use std::fmt::Write as _;

use cad::sketch::Region2;
use cad::solid::{Solid, Tpms};

/// The generated source, and what it needs to be called with.
#[derive(Clone, Debug)]
pub struct Generated {
    /// CUDA source defining `__device__ double scene_sdf(double3 p)`.
    pub source: String,
    /// Number of device functions emitted, a rough measure of model size.
    pub functions: usize,
}

struct Emitter {
    out: String,
    next: usize,
}

fn v3(p: [f64; 3]) -> String {
    format!("make_double3({:.17e}, {:.17e}, {:.17e})", p[0], p[1], p[2])
}

fn v2(p: [f64; 2]) -> String {
    format!("make_double2({:.17e}, {:.17e})", p[0], p[1])
}

/// A literal, wrapped so that it takes the field's precision: in `float`
/// mode an unwrapped double literal would promote the whole expression
/// back to double and throw the speed away.
fn f(x: f64) -> String {
    if x.is_infinite() {
        if x > 0.0 { "VINF".into() } else { "(-VINF)".into() }
    } else {
        format!("R({x:.17e})")
    }
}

impl Emitter {
    fn new() -> Self {
        Self { out: String::new(), next: 0 }
    }

    fn fresh(&mut self, prefix: &str) -> String {
        let name = format!("{prefix}{}", self.next);
        self.next += 1;
        name
    }

    /// Emit a 3D node, returning the name of its device function.
    fn solid(&mut self, s: &Solid) -> Result<String, String> {
        // Children first, so their definitions precede the parent's.
        let body = match s {
            Solid::Empty => "return VINF;".to_string(),
            Solid::Full => "return -VINF;".to_string(),
            Solid::HalfSpace { normal, offset } => format!(
                "double3 n = {}; double s = v3len(n); if (s == R(0.0)) return VINF; \
                 return dot3(scale3(n, 1.0 / s), p) - {} / s;",
                v3(*normal),
                f(*offset)
            ),
            Solid::Sphere { center, radius } => {
                format!("return v3len(sub3(p, {})) - {};", v3(*center), f(*radius))
            }
            Solid::Cuboid { center, half, round } => format!(
                "double3 d = sub3(p, {c}); double3 q = make_double3(fabs(d.x) - ({hx} - {r}), \
                 fabs(d.y) - ({hy} - {r}), fabs(d.z) - ({hz} - {r})); \
                 double3 o = make_double3(fmax(q.x, 0.0), fmax(q.y, 0.0), fmax(q.z, 0.0)); \
                 return v3len(o) + fmin(fmax(fmax(q.x, q.y), q.z), 0.0) - {r};",
                c = v3(*center),
                hx = f(half[0]),
                hy = f(half[1]),
                hz = f(half[2]),
                r = f(*round)
            ),
            Solid::Cylinder { a, b, radius } => format!(
                "return sd_cylinder(p, {}, {}, {});",
                v3(*a),
                v3(*b),
                f(*radius)
            ),
            Solid::Cone { a, b, radius_a, radius_b } => format!(
                "return sd_cone(p, {}, {}, {}, {});",
                v3(*a),
                v3(*b),
                f(*radius_a),
                f(*radius_b)
            ),
            Solid::Capsule { a, b, radius } => format!(
                "return sd_capsule(p, {}, {}, {});",
                v3(*a),
                v3(*b),
                f(*radius)
            ),
            Solid::Torus { center, major, minor } => format!(
                "double3 d = sub3(p, {}); double2 q = make_double2(sqrt(d.x * d.x + d.y * d.y) - {}, d.z); \
                 return sqrt(q.x * q.x + q.y * q.y) - {};",
                v3(*center),
                f(*major),
                f(*minor)
            ),
            Solid::Extrude { region, z0, z1 } => {
                let r = self.region(region)?;
                format!(
                    "double zc = 0.5 * ({z0} + {z1}); double hh = 0.5 * fabs({z1} - ({z0})); \
                     double d = {r}(make_double2(p.x, p.y)); \
                     double2 w = make_double2(d, fabs(p.z - zc) - hh); \
                     double2 o = make_double2(fmax(w.x, 0.0), fmax(w.y, 0.0)); \
                     return sqrt(o.x * o.x + o.y * o.y) + fmin(fmax(w.x, w.y), 0.0);",
                    z0 = f(*z0),
                    z1 = f(*z1),
                    r = r
                )
            }
            Solid::Revolve { region } => {
                let r = self.region(region)?;
                format!("return {r}(make_double2(sqrt(p.x * p.x + p.y * p.y), p.z));")
            }
            Solid::Union(parts) => {
                let names = self.children(parts)?;
                let mut body = "REAL d = VINF;".to_string();
                for n in names {
                    let _ = write!(body, " d = fmin(d, {n}(p));");
                }
                body + " return d;"
            }
            Solid::Intersection(parts) => {
                let names = self.children(parts)?;
                let mut body = "REAL d = -VINF;".to_string();
                for n in names {
                    let _ = write!(body, " d = fmax(d, {n}(p));");
                }
                body + " return d;"
            }
            Solid::Difference(a, b) => {
                let (na, nb) = (self.solid(a)?, self.solid(b)?);
                format!("return fmax({na}(p), -{nb}(p));")
            }
            Solid::SmoothUnion { a, b, k } => {
                let (na, nb) = (self.solid(a)?, self.solid(b)?);
                format!(
                    "double da = {na}(p), db = {nb}(p); double k = {k}; \
                     if (k <= 0.0) return fmin(da, db); \
                     double h = fmin(fmax(0.5 + 0.5 * (db - da) / k, 0.0), 1.0); \
                     return db + (da - db) * h - k * h * (1.0 - h);",
                    k = f(*k)
                )
            }
            Solid::Offset { solid, distance } => {
                let n = self.solid(solid)?;
                format!("return {n}(p) - {};", f(*distance))
            }
            Solid::Shell { solid, thickness } => {
                let n = self.solid(solid)?;
                format!("return fabs({n}(p)) - 0.5 * {};", f(*thickness))
            }
            Solid::Placed { solid, placement } => {
                let n = self.solid(solid)?;
                let r = placement.rot;
                format!(
                    "double3 d = sub3(p, {t}); \
                     double3 q = make_double3(({r00} * d.x + {r10} * d.y + {r20} * d.z) / {s}, \
                                              ({r01} * d.x + {r11} * d.y + {r21} * d.z) / {s}, \
                                              ({r02} * d.x + {r12} * d.y + {r22} * d.z) / {s}); \
                     return {s} * {n}(q);",
                    t = v3(placement.translate),
                    s = f(placement.scale),
                    r00 = f(r[0][0]), r01 = f(r[0][1]), r02 = f(r[0][2]),
                    r10 = f(r[1][0]), r11 = f(r[1][1]), r12 = f(r[1][2]),
                    r20 = f(r[2][0]), r21 = f(r[2][1]), r22 = f(r[2][2]),
                    n = n
                )
            }
            Solid::Stretch { solid, factors } => {
                let n = self.solid(solid)?;
                format!(
                    "double m = fmin(fmin(fabs({fx}), fabs({fy})), fabs({fz})); \
                     double3 q = make_double3(p.x / {fx}, p.y / {fy}, p.z / {fz}); \
                     return m * {n}(q);",
                    fx = f(factors[0]),
                    fy = f(factors[1]),
                    fz = f(factors[2]),
                    n = n
                )
            }
            Solid::Twist { solid, rate, radius } => {
                let n = self.solid(solid)?;
                format!(
                    "double a = -({rate}) * p.z; double s = sin(a), c = cos(a); \
                     double3 q = make_double3(c * p.x - s * p.y, s * p.x + c * p.y, p.z); \
                     return {n}(q) / (1.0 + fabs({rate}) * fabs({radius}));",
                    rate = f(*rate),
                    radius = f(*radius),
                    n = n
                )
            }
            Solid::Array { solid, period, counts } => {
                let n = self.solid(solid)?;
                format!(
                    "double3 per = {per}; int cnt[3] = {{{c0}, {c1}, {c2}}}; \
                     double pc[3] = {{p.x, p.y, p.z}}; double pe[3] = {{per.x, per.y, per.z}}; \
                     int base[3]; \
                     for (int k = 0; k < 3; ++k) {{ \
                         base[k] = (cnt[k] <= 1 || pe[k] == 0.0) ? 0 \
                                 : (int)nearbyint(pc[k] / pe[k] + 0.5 * ((double)cnt[k] - 1.0)); }} \
                     REAL best = VINF; \
                     for (int di = -1; di <= 1; ++di) for (int dj = -1; dj <= 1; ++dj) \
                     for (int dk = -1; dk <= 1; ++dk) {{ \
                         int c[3] = {{base[0] + di, base[1] + dj, base[2] + dk}}; \
                         bool ok = true; \
                         for (int k = 0; k < 3; ++k) if (c[k] < 0 || c[k] >= cnt[k]) ok = false; \
                         if (!ok) continue; \
                         double3 t = make_double3( \
                             ((double)c[0] - 0.5 * ((double)cnt[0] - 1.0)) * per.x, \
                             ((double)c[1] - 0.5 * ((double)cnt[1] - 1.0)) * per.y, \
                             ((double)c[2] - 0.5 * ((double)cnt[2] - 1.0)) * per.z); \
                         best = fmin(best, {n}(sub3(p, t))); }} \
                     return best;",
                    per = v3(*period),
                    c0 = counts[0],
                    c1 = counts[1],
                    c2 = counts[2],
                    n = n
                )
            }
            Solid::Periodic { kind, cell, level, thickness } => {
                let field = match kind {
                    Tpms::Gyroid => "double fv = sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x); \
                                     double lip = 2.0 * wn;",
                    Tpms::SchwarzP => "double fv = cos(x) + cos(y) + cos(z); double lip = wn;",
                    Tpms::Diamond => {
                        "double fv = sin(x) * sin(y) * sin(z) + sin(x) * cos(y) * cos(z) \
                         + cos(x) * sin(y) * cos(z) + cos(x) * cos(y) * sin(z); \
                         double lip = 4.0 * wn;"
                    }
                };
                let sheet = match thickness {
                    None => "return g;".to_string(),
                    Some(t) => format!("return fabs(g) - 0.5 * {};", f(*t)),
                };
                format!(
                    "double wx = 6.283185307179586476925286766559e0 / {cx}; \
                     double wy = 6.283185307179586476925286766559e0 / {cy}; \
                     double wz = 6.283185307179586476925286766559e0 / {cz}; \
                     double x = wx * p.x, y = wy * p.y, z = wz * p.z; \
                     double wn = sqrt(wx * wx + wy * wy + wz * wz); \
                     {field} double g = (fv - ({level})) / lip; {sheet}",
                    cx = f(cell[0]),
                    cy = f(cell[1]),
                    cz = f(cell[2]),
                    level = f(*level),
                    field = field,
                    sheet = sheet
                )
            }
            Solid::Mesh(_) => {
                return Err(
                    "Solid::Mesh carries a triangle BVH, which has no device form; \
                     render this model on the CPU"
                        .into(),
                )
            }
            Solid::Custom(c) => {
                return Err(format!(
                    "Solid::Custom(\"{}\") is an arbitrary host closure, which has no device form; \
                     render this model on the CPU",
                    c.name
                ))
            }
        };
        let name = self.fresh("s");
        let _ = writeln!(
            self.out,
            "__device__ double {name}(double3 p) {{ {body} }}"
        );
        Ok(name)
    }

    fn children(&mut self, parts: &[Solid]) -> Result<Vec<String>, String> {
        parts.iter().map(|p| self.solid(p)).collect()
    }

    fn region_children(&mut self, parts: &[Region2]) -> Result<Vec<String>, String> {
        parts.iter().map(|p| self.region(p)).collect()
    }

    /// Emit a 2D node, returning the name of its device function.
    fn region(&mut self, r: &Region2) -> Result<String, String> {
        let body = match r {
            Region2::Empty => "return VINF;".to_string(),
            Region2::Circle { center, radius } => {
                format!("return v2len(sub2(p, {})) - {};", v2(*center), f(*radius))
            }
            Region2::Rect { center, half, round } => format!(
                "double2 q = make_double2(fabs(p.x - {cx}) - ({hx} - {r}), \
                 fabs(p.y - {cy}) - ({hy} - {r})); \
                 double2 o = make_double2(fmax(q.x, 0.0), fmax(q.y, 0.0)); \
                 return v2len(o) + fmin(fmax(q.x, q.y), 0.0) - {r};",
                cx = f(center[0]),
                cy = f(center[1]),
                hx = f(half[0]),
                hy = f(half[1]),
                r = f(*round)
            ),
            Region2::RegularPolygon { center, radius, sides, phase } => {
                let poly = Region2::regular_polygon(*center, *radius, *sides, *phase);
                return self.region(&Region2::Polygon(poly));
            }
            Region2::Polygon(poly) => {
                // Loops are flattened into one array with an index table,
                // so the device code walks exactly the Rust segment order.
                let mut points: Vec<[f64; 2]> = Vec::new();
                let mut starts: Vec<usize> = Vec::new();
                for l in &poly.loops {
                    starts.push(points.len());
                    points.extend(l.iter().copied());
                }
                starts.push(points.len());
                if points.is_empty() {
                    return self.region(&Region2::Empty);
                }
                let pts = points
                    .iter()
                    .map(|p| format!("{{{:.17e}, {:.17e}}}", p[0], p[1]))
                    .collect::<Vec<_>>()
                    .join(", ");
                let st = starts.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(", ");
                format!(
                    "const double pts[{np}][2] = {{{pts}}}; \
                     const int st[{ns}] = {{{st}}}; \
                     REAL best = VINF; bool inside = false; \
                     for (int l = 0; l + 1 < {ns}; ++l) {{ \
                         int lo = st[l], hi = st[l + 1]; int n = hi - lo; \
                         if (n <= 0) continue; \
                         for (int i = 0; i < n; ++i) {{ \
                             double2 a = make_double2(pts[lo + i][0], pts[lo + i][1]); \
                             double2 b = make_double2(pts[lo + (i + 1) % n][0], \
                                                      pts[lo + (i + 1) % n][1]); \
                             double2 e = sub2(b, a); double ee = e.x * e.x + e.y * e.y; \
                             double2 pa = sub2(p, a); \
                             double t = ee > 0.0 ? fmin(fmax((pa.x * e.x + pa.y * e.y) / ee, 0.0), 1.0) : 0.0; \
                             double2 c = make_double2(a.x + e.x * t, a.y + e.y * t); \
                             double d = v2len(sub2(p, c)); if (d < best) best = d; \
                             if ((a.y > p.y) != (b.y > p.y)) {{ \
                                 double xx = a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x); \
                                 if (p.x < xx) inside = !inside; }} }} }} \
                     return inside ? -best : best;",
                    np = points.len(),
                    ns = starts.len(),
                    pts = pts,
                    st = st
                )
            }
            Region2::HalfPlane { normal, offset } => format!(
                "double2 n = {}; double s = v2len(n); if (s == R(0.0)) return VINF; \
                 return (n.x * p.x + n.y * p.y) / s - {} / s;",
                v2(*normal),
                f(*offset)
            ),
            Region2::Slot { a, b, radius } => format!(
                "double2 a = {a}, b = {b}; double2 e = sub2(b, a); \
                 double ee = e.x * e.x + e.y * e.y; double2 pa = sub2(p, a); \
                 double t = ee > 0.0 ? fmin(fmax((pa.x * e.x + pa.y * e.y) / ee, 0.0), 1.0) : 0.0; \
                 double2 c = make_double2(a.x + e.x * t, a.y + e.y * t); \
                 return v2len(sub2(p, c)) - {r};",
                a = v2(*a),
                b = v2(*b),
                r = f(*radius)
            ),
            Region2::Ring { center, inner, outer } => format!(
                "double r = v2len(sub2(p, {})); return fmax(r - {}, {} - r);",
                v2(*center),
                f(*outer),
                f(*inner)
            ),
            Region2::Union { parts } => {
                let names = self.region_children(parts)?;
                let mut body = "REAL d = VINF;".to_string();
                for n in names {
                    let _ = write!(body, " d = fmin(d, {n}(p));");
                }
                body + " return d;"
            }
            Region2::Intersection { parts } => {
                let names = self.region_children(parts)?;
                let mut body = "REAL d = -VINF;".to_string();
                for n in names {
                    let _ = write!(body, " d = fmax(d, {n}(p));");
                }
                body + " return d;"
            }
            Region2::Difference { a, b } => {
                let (na, nb) = (self.region(a)?, self.region(b)?);
                format!("return fmax({na}(p), -{nb}(p));")
            }
            Region2::Offset { region, distance } => {
                let n = self.region(region)?;
                format!("return {n}(p) - {};", f(*distance))
            }
            Region2::Shell { region, thickness } => {
                let n = self.region(region)?;
                format!("return fabs({n}(p)) - 0.5 * {};", f(*thickness))
            }
            Region2::Placed { region, angle, translate } => {
                let n = self.region(region)?;
                format!(
                    "double2 d = sub2(p, {t}); double s = sin(-({a})), c = cos(-({a})); \
                     double2 q = make_double2(c * d.x - s * d.y, s * d.x + c * d.y); \
                     return {n}(q);",
                    t = v2(*translate),
                    a = f(*angle),
                    n = n
                )
            }
            Region2::Array { region, period, counts } => {
                let n = self.region(region)?;
                format!(
                    "double2 per = {per}; int cnt[2] = {{{c0}, {c1}}}; \
                     double pc[2] = {{p.x, p.y}}; double pe[2] = {{per.x, per.y}}; \
                     int base[2]; \
                     for (int k = 0; k < 2; ++k) {{ \
                         base[k] = (cnt[k] <= 1 || pe[k] == 0.0) ? 0 \
                                 : (int)nearbyint(pc[k] / pe[k] + 0.5 * ((double)cnt[k] - 1.0)); }} \
                     REAL best = VINF; \
                     for (int di = -1; di <= 1; ++di) for (int dj = -1; dj <= 1; ++dj) {{ \
                         int c0i = base[0] + di, c1i = base[1] + dj; \
                         if (c0i < 0 || c1i < 0 || c0i >= cnt[0] || c1i >= cnt[1]) continue; \
                         double2 t = make_double2( \
                             ((double)c0i - 0.5 * ((double)cnt[0] - 1.0)) * per.x, \
                             ((double)c1i - 0.5 * ((double)cnt[1] - 1.0)) * per.y); \
                         best = fmin(best, {n}(sub2(p, t))); }} \
                     return best;",
                    per = v2(*period),
                    c0 = counts[0],
                    c1 = counts[1],
                    n = n
                )
            }
            Region2::PolarRepeat { region, count } => {
                let n = self.region(region)?;
                format!(
                    "int count = {c}; if (count < 1) count = 1; \
                     double step = 6.283185307179586476925286766559e0 / (double)count; \
                     double k = nearbyint(atan2(p.y, p.x) / step); \
                     REAL best = VINF; \
                     for (int d = -1; d <= 1; ++d) {{ \
                         double a = -(k + (double)d) * step; \
                         double s = sin(a), c = cos(a); \
                         double2 q = make_double2(c * p.x - s * p.y, s * p.x + c * p.y); \
                         best = fmin(best, {n}(q)); }} \
                     return best;",
                    c = (*count).max(1),
                    n = n
                )
            }
        };
        let name = self.fresh("r");
        let _ = writeln!(
            self.out,
            "__device__ double {name}(double2 p) {{ {body} }}"
        );
        Ok(name)
    }
}

/// The fixed preamble: the vector helpers and the exact primitive
/// distances, transcribed from `cad::solid` so the arithmetic matches.
pub const PREAMBLE: &str = include_str!("preamble.cu");

/// What precision the device field is evaluated in.
///
/// `Double` is the one the tests check bit for bit against the host, and
/// the default for that reason. `Single` is about four times faster on a
/// consumer card, which runs `double` at a fraction of its `float` rate,
/// and its deviation from the exact path is measured rather than assumed
/// (`tests/gpu.rs`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Precision {
    Double,
    Single,
}

impl Precision {
    fn prelude(&self) -> &'static str {
        match self {
            Precision::Double => {
                "typedef double REAL;\ntypedef double2 REAL2;\ntypedef double3 REAL3;\n\
                 #define mk2 make_double2\n#define mk3 make_double3\n\
                 #define R(x) (x)\n\
                 #define VINF (__longlong_as_double(0x7ff0000000000000LL))\n"
            }
            Precision::Single => {
                "typedef float REAL;\ntypedef float2 REAL2;\ntypedef float3 REAL3;\n\
                 #define mk2 make_float2\n#define mk3 make_float3\n\
                 #define R(x) ((float)(x))\n\
                 #define VINF (__int_as_float(0x7f800000))\n"
            }
        }
    }
}

/// Rewrite the `double`-spelled source into the precision-agnostic one.
///
/// The generator and the hand-written preamble are both written in
/// `double`, which keeps them readable next to the Rust they mirror; this
/// is the one place that knows about `float`. `__longlong_as_double` is
/// left alone because it has no word boundary before `double`.
///
/// Floating literals are wrapped in `R(...)` as well. That is not
/// cosmetic: a bare `0.0` next to a `float` makes `fmin`/`fmax`
/// ambiguous between their two overloads, and in plain arithmetic it
/// promotes the whole expression back to `double` and throws away the
/// speed the single-precision mode exists for.
fn to_real(src: &str) -> String {
    let s = src
        .replace("make_double3(", "mk3(")
        .replace("make_double2(", "mk2(");
    wrap_literals(&regexless_double(&s))
}

/// Wrap every floating-point literal in `R(...)`, leaving integers and hex
/// alone. Already-wrapped literals become `R(R(x))`, which is harmless.
fn wrap_literals(s: &str) -> String {
    let b = s.as_bytes();
    let mut out = String::with_capacity(s.len() + s.len() / 8);
    let word = |c: u8| c.is_ascii_alphanumeric() || c == b'_';
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        // Hex literals pass through untouched.
        if c == b'0' && i + 1 < b.len() && (b[i + 1] | 32) == b'x' {
            let start = i;
            i += 2;
            while i < b.len() && (b[i].is_ascii_hexdigit() || (b[i] | 32) == b'l' || (b[i] | 32) == b'u') {
                i += 1;
            }
            out.push_str(&s[start..i]);
            continue;
        }
        let fresh_number = c.is_ascii_digit() && (i == 0 || (!word(b[i - 1]) && b[i - 1] != b'.'));
        if fresh_number {
            let start = i;
            while i < b.len() && b[i].is_ascii_digit() {
                i += 1;
            }
            let mut is_float = false;
            if i < b.len() && b[i] == b'.' {
                is_float = true;
                i += 1;
                while i < b.len() && b[i].is_ascii_digit() {
                    i += 1;
                }
            }
            if i < b.len() && (b[i] | 32) == b'e' {
                let save = i;
                i += 1;
                if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
                    i += 1;
                }
                if i < b.len() && b[i].is_ascii_digit() {
                    is_float = true;
                    while i < b.len() && b[i].is_ascii_digit() {
                        i += 1;
                    }
                } else {
                    i = save;
                }
            }
            let lit = &s[start..i];
            if is_float {
                out.push_str("R(");
                out.push_str(lit);
                out.push(')');
            } else {
                out.push_str(lit);
            }
            continue;
        }
        let ch = s[i..].chars().next().unwrap();
        out.push(ch);
        i += ch.len_utf8();
    }
    out
}

/// Replace whole-word `double`/`double2`/`double3` without a regex crate.
fn regexless_double(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let bytes = s.as_bytes();
    let mut i = 0;
    let word = |b: u8| b.is_ascii_alphanumeric() || b == b'_';
    while i < bytes.len() {
        let rest = &s[i..];
        let before_ok = i == 0 || !word(bytes[i - 1]);
        let mut matched = None;
        if before_ok {
            for (pat, rep) in [("double3", "REAL3"), ("double2", "REAL2"), ("double", "REAL")] {
                if rest.starts_with(pat) {
                    let after = bytes.get(i + pat.len()).copied().unwrap_or(b' ');
                    if !word(after) {
                        matched = Some((pat.len(), rep));
                        break;
                    }
                }
            }
        }
        match matched {
            Some((len, rep)) => {
                out.push_str(rep);
                i += len;
            }
            None => {
                let ch = rest.chars().next().unwrap();
                out.push(ch);
                i += ch.len_utf8();
            }
        }
    }
    out
}

/// Compile a solid into CUDA source defining `scene_sdf`, in double
/// precision.
pub fn emit(solid: &Solid) -> Result<Generated, String> {
    emit_with(solid, Precision::Double)
}

/// Compile a solid into CUDA source defining `scene_sdf`.
pub fn emit_with(solid: &Solid, precision: Precision) -> Result<Generated, String> {
    let mut e = Emitter::new();
    let root = e.solid(solid)?;
    let mut body = String::with_capacity(PREAMBLE.len() + e.out.len() + 256);
    body.push_str(PREAMBLE);
    body.push_str(&e.out);
    let _ = writeln!(
        body,
        "__device__ double scene_sdf(double3 p) {{ return {root}(p); }}"
    );
    let mut source = String::with_capacity(body.len() + 512);
    source.push_str(precision.prelude());
    source.push_str(&to_real(&body));
    Ok(Generated { source, functions: e.next })
}
