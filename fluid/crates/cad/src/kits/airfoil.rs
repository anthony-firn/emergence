//! Aerofoils and blades.
//!
//! The NACA four-digit series is the one shape family in this crate that
//! is *defined* by a published polynomial rather than by a drawing, which
//! makes it unusually well suited to validation: the thickness
//! distribution, the maximum-thickness position, the leading-edge radius
//! and the enclosed area are all exact closed forms, and the tests check
//! the generated section against every one of them.
//!
//! ```text
//!   NACA MPXX -> camber line + thickness -> Region2 -> Solid (extrude)
//!                                        -> blade: stacked, twisted,
//!                                           tapered sections -> loft
//! ```
//!
//! # The definition
//!
//! For a section `NACA mpxx` of chord `c`, with `m = M/100` the maximum
//! camber, `p = P/10` its position and `t = XX/100` the thickness
//! (Jacobs, Ward & Pinkerton, *NACA Report 460*, 1933; Abbott & von
//! Doenhoff, *Theory of Wing Sections*, 1959, §6.3):
//!
//! ```text
//!   y_t(x) = 5 t (0.2969 sqrt(x) - 0.1260 x - 0.3516 x^2
//!                 + 0.2843 x^3 - a4 x^4)
//!   y_c(x) = m/p^2 (2 p x - x^2)                        for x <= p
//!          = m/(1-p)^2 ((1 - 2p) + 2 p x - x^2)         for x > p
//! ```
//!
//! with `x` in chords and `a4 = 0.1015` for the original open trailing
//! edge or `0.1036` for the closed one. The surfaces are offset
//! perpendicular to the camber line, which is the standard construction
//! and differs from a vertical offset by `O(dy_c/dx)^2`.
//!
//! Three consequences used as tests:
//!
//! * the thickness is greatest at `x = 0.2998` — the series' nominal
//!   "30 % chord" — where it equals `t c`;
//! * the leading edge is a parabola `y^2 = 2 r x` with
//!   `r = (5 t * 0.2969)^2 / 2 = 1.1019 t^2` chords, the published
//!   leading-edge radius;
//! * the enclosed area is `2 integral y_t dx = 0.685083 t c^2` open, or
//!   `0.680883 t c^2` closed, from the exact integral of the polynomial.
//!
//! # Limits
//!
//! No five-digit, six-series or modified four-digit sections; no
//! aerodynamic prediction of any kind — this draws shapes. A blade is a
//! loft through stacked sections, with no account of how a real one is
//! manufactured or how its sections should be stacked for a given duty.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use geometry::{add, scale, TriMesh};

use crate::path::{loft, Path3, Profile};
use crate::sketch::{Polygon, Region2, V2};
use crate::solid::Solid;

/// The four-digit thickness polynomial's coefficients.
const A: [f64; 5] = [0.2969, -0.1260, -0.3516, 0.2843, 0.0];
/// Last coefficient: open trailing edge (NACA Report 460).
pub const A4_OPEN: f64 = 0.1015;
/// Last coefficient that closes the trailing edge exactly.
pub const A4_CLOSED: f64 = 0.1036;

/// A NACA four-digit section.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct Naca4 {
    /// Maximum camber in chords (the first digit over 100).
    pub camber: f64,
    /// Position of maximum camber in chords (the second digit over 10).
    pub camber_position: f64,
    /// Maximum thickness in chords (the last two digits over 100).
    pub thickness: f64,
    pub chord: f64,
    /// `true` to close the trailing edge exactly.
    #[serde(default)]
    pub closed: bool,
    /// Points per surface.
    #[serde(default = "points")]
    pub points: usize,
}

fn points() -> usize {
    80
}

impl Naca4 {
    /// From the four digits, e.g. `Naca4::designation(2412, 1.0)`.
    pub fn designation(code: u32, chord: f64) -> Result<Self, String> {
        if code > 9999 {
            return Err(format!("{code} is not a four-digit designation"));
        }
        let m = (code / 1000) % 10;
        let p = (code / 100) % 10;
        let t = code % 100;
        if t == 0 {
            return Err("a section needs a non-zero thickness".into());
        }
        if m > 0 && p == 0 {
            return Err("a cambered section needs a camber position".into());
        }
        Ok(Self {
            camber: m as f64 / 100.0,
            camber_position: p as f64 / 10.0,
            thickness: t as f64 / 100.0,
            chord,
            closed: false,
            points: points(),
        })
    }

    pub fn symmetric(thickness: f64, chord: f64) -> Self {
        Self {
            camber: 0.0,
            camber_position: 0.0,
            thickness,
            chord,
            closed: false,
            points: points(),
        }
    }

    fn a4(&self) -> f64 {
        if self.closed {
            A4_CLOSED
        } else {
            A4_OPEN
        }
    }

    /// Half-thickness at `x` (in chords), in chords.
    pub fn half_thickness(&self, x: f64) -> f64 {
        let x = x.clamp(0.0, 1.0);
        5.0 * self.thickness
            * (A[0] * x.sqrt() + A[1] * x + A[2] * x * x + A[3] * x.powi(3) - self.a4() * x.powi(4))
    }

    /// Camber-line ordinate at `x` (in chords), in chords.
    pub fn camber_at(&self, x: f64) -> f64 {
        let (m, p) = (self.camber, self.camber_position);
        if m == 0.0 || p <= 0.0 {
            return 0.0;
        }
        let x = x.clamp(0.0, 1.0);
        if x <= p {
            m / (p * p) * (2.0 * p * x - x * x)
        } else {
            m / ((1.0 - p) * (1.0 - p)) * ((1.0 - 2.0 * p) + 2.0 * p * x - x * x)
        }
    }

    /// Camber-line slope at `x`.
    pub fn camber_slope(&self, x: f64) -> f64 {
        let (m, p) = (self.camber, self.camber_position);
        if m == 0.0 || p <= 0.0 {
            return 0.0;
        }
        let x = x.clamp(0.0, 1.0);
        if x <= p {
            2.0 * m / (p * p) * (p - x)
        } else {
            2.0 * m / ((1.0 - p) * (1.0 - p)) * (p - x)
        }
    }

    /// Upper and lower surface points at `x` (in chords), in metres, with
    /// the thickness applied perpendicular to the camber line.
    pub fn surfaces(&self, x: f64) -> (V2, V2) {
        let yt = self.half_thickness(x);
        let yc = self.camber_at(x);
        let theta = self.camber_slope(x).atan();
        let (s, co) = theta.sin_cos();
        (
            [self.chord * (x - yt * s), self.chord * (yc + yt * co)],
            [self.chord * (x + yt * s), self.chord * (yc - yt * co)],
        )
    }

    /// Leading-edge radius in metres: `1.1019 t^2` chords, exactly
    /// `(5 t a0)^2 / 2`.
    pub fn leading_edge_radius(&self) -> f64 {
        let k = 5.0 * self.thickness * A[0];
        0.5 * k * k * self.chord
    }

    /// Trailing-edge gap in metres; zero for a closed section.
    pub fn trailing_edge_gap(&self) -> f64 {
        2.0 * self.half_thickness(1.0).abs() * self.chord
    }

    /// `2 integral y_t dx` times `c^2`, in closed form: the area a
    /// *vertically* offset section would enclose.
    ///
    /// For the four-digit series this is `0.685083 t c^2` with the open
    /// trailing edge and `0.680883 t c^2` with the closed one. It is the
    /// section's true area only when the camber is zero — see
    /// [`Naca4::area`].
    pub fn area_vertical(&self) -> f64 {
        let i = A[0] * 2.0 / 3.0 + A[1] / 2.0 + A[2] / 3.0 + A[3] / 4.0 - self.a4() / 5.0;
        2.0 * 5.0 * self.thickness * i * self.chord * self.chord
    }

    /// Enclosed area of the section as actually constructed.
    ///
    /// The thickness is laid off perpendicular to the camber line, so the
    /// area is `2 integral y_t ds = 2 integral y_t sqrt(1 + y_c'^2) dx`,
    /// which **exceeds** [`Naca4::area_vertical`] for any cambered
    /// section — by 0.8 % for NACA 4415. There is no elementary closed
    /// form, so it is integrated: substituting `x = u^2` removes the
    /// leading edge's square-root singularity, and the range is split at
    /// the camber position, where the camber slope has a kink. With
    /// 5-point Gauss-Legendre on 32 intervals per piece the result is good
    /// to well under 1e-12 relative.
    pub fn area(&self) -> f64 {
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
        // In u with x = u^2: integrand 2 y_t(u^2) sqrt(1 + y_c'(u^2)^2) 2u.
        let f = |u: f64| {
            let x = u * u;
            let slope = self.camber_slope(x);
            2.0 * self.half_thickness(x) * (1.0 + slope * slope).sqrt() * 2.0 * u
        };
        let mut breaks = vec![0.0];
        if self.camber > 0.0 && self.camber_position > 0.0 && self.camber_position < 1.0 {
            breaks.push(self.camber_position.sqrt());
        }
        breaks.push(1.0);
        let mut total = 0.0;
        for w in breaks.windows(2) {
            let n = 32;
            let h = (w[1] - w[0]) / n as f64;
            for i in 0..n {
                let (a, b) = (w[0] + h * i as f64, w[0] + h * (i + 1) as f64);
                let (mid, half) = (0.5 * (a + b), 0.5 * (b - a));
                for k in 0..5 {
                    total += W[k] * f(mid + half * X[k]) * half;
                }
            }
        }
        total * self.chord * self.chord
    }

    /// Where the section is thickest, in chords, by bisection on the
    /// polynomial's derivative. The four-digit series gives 0.2998.
    pub fn max_thickness_position(&self) -> f64 {
        let dt = |x: f64| {
            0.5 * A[0] / x.sqrt()
                + A[1]
                + 2.0 * A[2] * x
                + 3.0 * A[3] * x * x
                - 4.0 * self.a4() * x.powi(3)
        };
        let (mut lo, mut hi) = (0.05, 0.8);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if dt(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Cosine-spaced sample abscissae, which cluster points at the
    /// leading edge where the curvature is highest.
    fn abscissae(&self) -> Vec<f64> {
        let n = self.points.max(8);
        (0..=n).map(|i| 0.5 * (1.0 - (PI * i as f64 / n as f64).cos())).collect()
    }

    /// The closed section outline: upper surface from nose to tail, lower
    /// surface back.
    pub fn polygon(&self) -> Polygon {
        let xs = self.abscissae();
        let mut pts: Vec<V2> = Vec::with_capacity(2 * xs.len());
        for &x in &xs {
            pts.push(self.surfaces(x).0);
        }
        // Skip the duplicated trailing edge when it is closed.
        let skip_tail = self.trailing_edge_gap() < 1e-12 * self.chord;
        for &x in xs.iter().rev().skip(usize::from(skip_tail)) {
            if x <= 0.0 {
                continue; // the nose is already in
            }
            pts.push(self.surfaces(x).1);
        }
        Polygon::new(pts)
    }

    pub fn region(&self) -> Region2 {
        Region2::Polygon(self.polygon())
    }

    /// The section extruded into a constant-chord wing of `span`.
    pub fn wing(&self, span: f64) -> Solid {
        Solid::extrude(self.region(), 0.0, span)
    }

    /// The outline as a sweep profile, for lofting a blade.
    pub fn profile(&self, twist: f64, scale_chord: f64) -> Profile {
        let (s, c) = twist.sin_cos();
        Profile::new(
            self.polygon()
                .loops[0]
                .iter()
                .map(|p| {
                    let (x, y) = (p[0] * scale_chord / self.chord, p[1] * scale_chord / self.chord);
                    [c * x - s * y, s * x + c * y]
                })
                .collect(),
        )
    }
}

/// A blade or wing lofted through stacked sections.
#[derive(Clone, Debug)]
pub struct Blade {
    /// Section at the root; the same shape is used at every station.
    pub section: Naca4,
    pub span: f64,
    /// Chord at root and tip.
    pub root_chord: f64,
    pub tip_chord: f64,
    /// Twist at root and tip, radians (positive leading edge up).
    pub root_twist: f64,
    pub tip_twist: f64,
    /// Sweep of the quarter-chord line, radians.
    pub sweep: f64,
    pub stations: usize,
}

impl Blade {
    pub fn new(section: Naca4, span: f64, root_chord: f64, tip_chord: f64) -> Self {
        Self {
            section,
            span,
            root_chord,
            tip_chord,
            root_twist: 0.0,
            tip_twist: 0.0,
            sweep: 0.0,
            stations: 12,
        }
    }

    pub fn twisted(mut self, root: f64, tip: f64) -> Self {
        self.root_twist = root;
        self.tip_twist = tip;
        self
    }

    pub fn swept(mut self, sweep: f64) -> Self {
        self.sweep = sweep;
        self
    }

    /// Chord at a span fraction (linear taper).
    pub fn chord_at(&self, t: f64) -> f64 {
        self.root_chord + (self.tip_chord - self.root_chord) * t
    }

    /// Planform area, exactly the trapezoid `span (c_root + c_tip)/2`.
    pub fn planform_area(&self) -> f64 {
        0.5 * self.span * (self.root_chord + self.tip_chord)
    }

    /// Aspect ratio `b^2 / S` for a full wing of this half-span pair.
    pub fn aspect_ratio(&self) -> f64 {
        self.span * self.span / self.planform_area()
    }

    /// Taper ratio.
    pub fn taper(&self) -> f64 {
        self.tip_chord / self.root_chord
    }

    /// Mean aerodynamic chord of a linearly tapered planform:
    /// `(2/3) c_r (1 + l + l^2)/(1 + l)`.
    pub fn mean_aerodynamic_chord(&self) -> f64 {
        let l = self.taper();
        2.0 / 3.0 * self.root_chord * (1.0 + l + l * l) / (1.0 + l)
    }

    /// The lofted surface. Sections are stacked along `+z` with the
    /// quarter-chord line swept in `x`.
    pub fn mesh(&self) -> TriMesh {
        let n = self.stations.max(2);
        let spine: Vec<[f64; 3]> = (0..=n)
            .map(|i| {
                let t = i as f64 / n as f64;
                let z = self.span * t;
                [z * self.sweep.tan(), 0.0, z]
            })
            .collect();
        let path = Path3::new(spine, false);
        let profiles: Vec<Profile> = (0..=n)
            .map(|i| {
                let t = i as f64 / n as f64;
                let twist = self.root_twist + (self.tip_twist - self.root_twist) * t;
                // The spine runs along +z, so the frame's u axis is the
                // chordwise direction and v the thickness direction.
                let p = self.section.profile(twist, self.chord_at(t));
                // Place the quarter-chord at the spine.
                Profile::new(
                    p.points
                        .iter()
                        .map(|q| [q[0] - 0.25 * self.chord_at(t), q[1]])
                        .collect(),
                )
            })
            .collect();
        loft(&path, &profiles)
    }

    pub fn solid(&self) -> Solid {
        Solid::from_trimesh(self.mesh())
    }
}

/// Offset a 2D point by a vector, for callers assembling sections.
pub fn offset(p: V2, by: V2) -> V2 {
    let a = add([p[0], p[1], 0.0], [by[0], by[1], 0.0]);
    [a[0], a[1]]
}

/// Scale a 2D point about the origin.
pub fn scaled(p: V2, s: f64) -> V2 {
    let a = scale([p[0], p[1], 0.0], s);
    [a[0], a[1]]
}
