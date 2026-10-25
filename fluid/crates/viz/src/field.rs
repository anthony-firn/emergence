//! Simulation fields as something to look at.
//!
//! A solver produces numbers on a grid; a render needs to show them. Three
//! ways, each answering a different question:
//!
//! * **surface colouring** — sample the field where a ray hits geometry,
//!   which is what [`crate::Pass::Field`] records and what shows a wall
//!   temperature or a surface pressure;
//! * **isosurfaces** — hand the field to the CAD kernel as a
//!   [`cad::Solid`] and let it be shaded like any other solid, so a
//!   temperature contour can be cut, shelled and rendered with the part
//!   it came from;
//! * **colour maps** — turn a scalar into RGB, with the perceptual
//!   caveats stated below.
//!
//! # Interpolation
//!
//! A grid field is sampled trilinearly, which is continuous but has a
//! discontinuous gradient at cell faces. For a *rendered isosurface* that
//! matters: the normal, which comes from the gradient, is faceted at cell
//! boundaries. [`Field::as_solid`] therefore reports the Lipschitz bound
//! it can honestly claim, and a field isosurface is marked inexact so the
//! kernel does not treat its value as a distance.

use cad::solid::Solid;
use geometry::P3;
use serde::{Deserialize, Serialize};

/// A scalar field on a uniform grid, `x` fastest.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Field {
    pub origin: P3,
    pub spacing: P3,
    /// Node counts.
    pub dims: [usize; 3],
    pub values: Vec<f64>,
    /// Value returned outside the grid.
    #[serde(default)]
    pub outside: f64,
}

impl Field {
    pub fn new(origin: P3, spacing: P3, dims: [usize; 3], values: Vec<f64>) -> Result<Self, String> {
        let n = dims[0] * dims[1] * dims[2];
        if values.len() != n {
            return Err(format!("{} values for a {dims:?} grid, expected {n}", values.len()));
        }
        if (0..3).any(|k| dims[k] < 2 || !(spacing[k] > 0.0)) {
            return Err("every axis needs at least two nodes and a positive spacing".into());
        }
        Ok(Self {
            origin,
            spacing,
            dims,
            values,
            outside: 0.0,
        })
    }

    /// Sample a field from a closure, which is how a solver's state or an
    /// analytic test case gets in.
    pub fn from_fn(
        origin: P3,
        spacing: P3,
        dims: [usize; 3],
        f: impl Fn(P3) -> f64,
    ) -> Result<Self, String> {
        let mut values = Vec::with_capacity(dims[0] * dims[1] * dims[2]);
        for k in 0..dims[2] {
            for j in 0..dims[1] {
                for i in 0..dims[0] {
                    values.push(f([
                        origin[0] + spacing[0] * i as f64,
                        origin[1] + spacing[1] * j as f64,
                        origin[2] + spacing[2] * k as f64,
                    ]));
                }
            }
        }
        Self::new(origin, spacing, dims, values)
    }

    pub fn bounds(&self) -> (P3, P3) {
        (
            self.origin,
            std::array::from_fn(|k| self.origin[k] + self.spacing[k] * (self.dims[k] - 1) as f64),
        )
    }

    #[inline]
    fn node(&self, i: usize, j: usize, k: usize) -> f64 {
        self.values[i + self.dims[0] * (j + self.dims[1] * k)]
    }

    /// Trilinear sample. Outside the grid, `outside`.
    pub fn at(&self, p: P3) -> f64 {
        let mut f = [0.0f64; 3];
        let mut i0 = [0usize; 3];
        for k in 0..3 {
            let x = (p[k] - self.origin[k]) / self.spacing[k];
            if x < 0.0 || x > (self.dims[k] - 1) as f64 {
                return self.outside;
            }
            let i = (x.floor() as usize).min(self.dims[k] - 2);
            i0[k] = i;
            f[k] = x - i as f64;
        }
        let mut acc = 0.0;
        for dk in 0..2 {
            for dj in 0..2 {
                for di in 0..2 {
                    let w = (if di == 0 { 1.0 - f[0] } else { f[0] })
                        * (if dj == 0 { 1.0 - f[1] } else { f[1] })
                        * (if dk == 0 { 1.0 - f[2] } else { f[2] });
                    acc += w * self.node(i0[0] + di, i0[1] + dj, i0[2] + dk);
                }
            }
        }
        acc
    }

    pub fn range(&self) -> (f64, f64) {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for &v in &self.values {
            if v.is_finite() {
                lo = lo.min(v);
                hi = hi.max(v);
            }
        }
        (lo, hi)
    }

    /// The largest `|value difference| / distance` between adjacent nodes:
    /// a lower bound on the field's Lipschitz constant, and what an
    /// isosurface's march has to be divided by to be safe.
    pub fn gradient_bound(&self) -> f64 {
        let mut worst: f64 = 0.0;
        for k in 0..self.dims[2] {
            for j in 0..self.dims[1] {
                for i in 0..self.dims[0] {
                    let v = self.node(i, j, k);
                    if i + 1 < self.dims[0] {
                        worst = worst.max((self.node(i + 1, j, k) - v).abs() / self.spacing[0]);
                    }
                    if j + 1 < self.dims[1] {
                        worst = worst.max((self.node(i, j + 1, k) - v).abs() / self.spacing[1]);
                    }
                    if k + 1 < self.dims[2] {
                        worst = worst.max((self.node(i, j, k + 1) - v).abs() / self.spacing[2]);
                    }
                }
            }
        }
        worst
    }

    /// The isosurface `value = level` as a solid the kernel can render,
    /// cut, shell and boolean like any other.
    ///
    /// The field is divided by its gradient bound so the result satisfies
    /// the kernel's 1-Lipschitz invariant, which is what makes sphere
    /// tracing it sound. It is marked **inexact**: the quotient is a
    /// conservative distance, not the true one, so the march is correct
    /// but takes more steps than it would on an exact field.
    pub fn as_solid(&self, level: f64, inside_is_below: bool) -> Result<Solid, String> {
        let bound = self.gradient_bound();
        if !(bound > 0.0) {
            return Err("the field is constant: it has no isosurface".into());
        }
        let field = self.clone();
        let (lo, hi) = self.bounds();
        let sign = if inside_is_below { 1.0 } else { -1.0 };
        let outside = self.outside;
        Ok(Solid::custom(
            "field isosurface",
            (lo, hi),
            false,
            move |p| {
                let v = if (0..3).all(|k| p[k] >= lo[k] && p[k] <= hi[k]) {
                    field.at(p)
                } else {
                    outside
                };
                sign * (v - level) / bound
            },
        ))
    }
}

/// A colour map from a scalar in `[0, 1]` to linear RGB.
///
/// `Viridis` is the perceptually uniform default: lightness increases
/// monotonically, so the map does not invent structure the data does not
/// have. `Turbo` is included because it is familiar, with the warning that
/// its lightness is not monotone and it therefore does create apparent
/// banding. `Grey` is honest and unambiguous. These are the standard
/// control points of each map, interpolated linearly in linear-light RGB.
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum ColourMap {
    Grey,
    Viridis,
    Turbo,
    /// Blue-white-red, centred: for a signed quantity where zero matters.
    Diverging,
}

const VIRIDIS: [[f64; 3]; 9] = [
    [0.267, 0.005, 0.329],
    [0.283, 0.141, 0.458],
    [0.254, 0.265, 0.530],
    [0.207, 0.372, 0.553],
    [0.164, 0.471, 0.558],
    [0.128, 0.567, 0.551],
    [0.135, 0.659, 0.518],
    [0.267, 0.749, 0.441],
    [0.993, 0.906, 0.144],
];

const TURBO: [[f64; 3]; 9] = [
    [0.190, 0.072, 0.232],
    [0.254, 0.395, 0.904],
    [0.147, 0.665, 0.981],
    [0.096, 0.879, 0.734],
    [0.459, 0.985, 0.406],
    [0.806, 0.919, 0.228],
    [0.985, 0.713, 0.169],
    [0.939, 0.357, 0.086],
    [0.479, 0.012, 0.011],
];

impl ColourMap {
    /// Map `t`, clamped to `[0, 1]`, to linear RGB.
    pub fn at(&self, t: f64) -> [f64; 3] {
        let t = t.clamp(0.0, 1.0);
        match self {
            ColourMap::Grey => [t, t, t],
            ColourMap::Viridis => interpolate(&VIRIDIS, t),
            ColourMap::Turbo => interpolate(&TURBO, t),
            ColourMap::Diverging => {
                if t < 0.5 {
                    let u = 2.0 * t;
                    [u, u, 1.0]
                } else {
                    let u = 2.0 * (1.0 - t);
                    [1.0, u, u]
                }
            }
        }
    }

    /// Whether the map's luminance increases monotonically, which is what
    /// decides if it can be read without inventing structure.
    pub fn perceptually_ordered(&self) -> bool {
        matches!(self, ColourMap::Grey | ColourMap::Viridis)
    }

    /// Relative luminance of the map at `t` (Rec. 709 weights).
    pub fn luminance(&self, t: f64) -> f64 {
        let c = self.at(t);
        0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]
    }
}

fn interpolate(table: &[[f64; 3]], t: f64) -> [f64; 3] {
    let n = table.len();
    let x = t * (n - 1) as f64;
    let i = (x.floor() as usize).min(n - 2);
    let f = x - i as f64;
    std::array::from_fn(|k| table[i][k] * (1.0 - f) + table[i + 1][k] * f)
}

/// Colour a channel of values into an image, with the range given or taken
/// from the data.
pub fn colourise(
    channel: &crate::image::Channel,
    map: ColourMap,
    range: Option<(f64, f64)>,
) -> crate::image::Image {
    let (lo, hi) = range.unwrap_or_else(|| channel.range());
    let span = if (hi - lo).abs() > 0.0 { hi - lo } else { 1.0 };
    let mut img = crate::image::Image::new(channel.width, channel.height);
    for (i, &v) in channel.values.iter().enumerate() {
        img.pixels[i] = if v.is_finite() {
            map.at((v - lo) / span)
        } else {
            [0.0; 3]
        };
    }
    img
}
