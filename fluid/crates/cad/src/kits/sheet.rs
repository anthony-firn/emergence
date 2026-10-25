//! Sheet metal: formed parts and the flat patterns they are cut from.
//!
//! Sheet metal is the one domain here where the *manufacturing* geometry is
//! a different shape from the part: a bent bracket is cut flat, and the
//! blank is longer than the sum of its flanges by the amount the bend
//! stretches. Getting that bookkeeping right is the whole job, and it has
//! a standard, published form.
//!
//! ```text
//!   flange lengths to the outside apex + bend angles, radii, K-factor
//!        -> neutral-axis path (straights + arcs)
//!        -> formed cross-section -> Region2 -> Solid
//!        -> flat pattern (the blank to cut)
//! ```
//!
//! # The bend formulas
//!
//! For a bend of included angle `A` (the angle turned through), inside
//! radius `R`, material thickness `T` and K-factor `K`:
//!
//! ```text
//!   bend allowance   BA   = A (R + K T)      the neutral axis's arc length
//!   outside setback  OSSB = (R + T) tan(A/2)
//!   bend deduction   BD   = 2 OSSB - BA
//!   flat length           = sum(flange to apex) - sum(BD)
//!                         = sum(straight) + sum(BA)
//! ```
//!
//! These are the standard shop formulas (Boljanovic, *Sheet Metal Forming
//! Processes and Die Design*, 2nd ed., 2014, ch. 2). The two expressions
//! for the flat length are algebraically equivalent but reach it by
//! different routes — one from the apex dimensions a drawing carries, the
//! other from the neutral-axis geometry — so the tests check them against
//! each other and against the generated path.
//!
//! # The K-factor is shop data, not a prediction
//!
//! `K` is where the neutral axis sits, as a fraction of the thickness from
//! the inside surface. It depends on the material, the tooling and the
//! ratio `R/T`, and published tables disagree by several per cent.
//! [`k_factor`] offers the common rule-of-thumb values with that caveat
//! attached; nothing here predicts `K`, and a production flat pattern
//! should use the value measured on the press that will form it.
//!
//! One consequence worth knowing: bending conserves volume only when
//! `K = 0.5`. At any other value the model stretches one surface and
//! compresses the other, so the formed solid's volume differs from the
//! blank's by a known amount, and [`SheetPart::volume_change`] reports it.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::sketch::{Polygon, Region2, V2};
use crate::solid::Solid;

/// Rule-of-thumb K-factors by bend radius over thickness.
///
/// These are the values commonly tabulated for mild steel (Boljanovic
/// 2014, table 2.3, and the equivalent table in most press-brake
/// handbooks). **Published tables differ by several per cent** and the
/// value depends on the tooling, so this is a starting point for a flat
/// pattern, not an answer.
pub fn k_factor(radius_over_thickness: f64) -> f64 {
    match radius_over_thickness {
        r if r < 1.0 => 0.33,
        r if r < 3.0 => 0.40,
        _ => 0.50,
    }
}

/// One bend.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct Bend {
    /// Angle turned through, radians. Positive turns the flange to the
    /// left (counter-clockwise in the section), negative to the right; a
    /// part with one of each is a Z-section. Only the magnitude affects
    /// the bend allowance.
    pub angle: f64,
    /// Inside bend radius.
    pub inside_radius: f64,
}

impl Bend {
    pub fn at_right_angle(inside_radius: f64) -> Self {
        Self { angle: 0.5 * PI, inside_radius }
    }
}

/// A part formed by bending along parallel lines: a section bent in the
/// plane and extruded along the bend axis. This is the case where the flat
/// pattern is exact.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct SheetPart {
    pub thickness: f64,
    /// Where the neutral axis sits, as a fraction of the thickness from
    /// the inside surface.
    pub k_factor: f64,
    /// Extrusion depth along the bend axis.
    pub width: f64,
    /// Flange lengths **to the outside apex**, as a drawing dimensions
    /// them: one more than the number of bends.
    pub flanges: Vec<f64>,
    pub bends: Vec<Bend>,
}

/// What one bend costs in material.
#[derive(Clone, Copy, Debug)]
pub struct BendNumbers {
    pub allowance: f64,
    pub outside_setback: f64,
    pub deduction: f64,
    /// Radius of the neutral axis through the bend.
    pub neutral_radius: f64,
}

impl SheetPart {
    pub fn new(thickness: f64, width: f64, flanges: Vec<f64>, bends: Vec<Bend>) -> Result<Self, String> {
        if flanges.len() != bends.len() + 1 {
            return Err(format!(
                "{} flanges need {} bends, got {}",
                flanges.len(),
                flanges.len() - 1,
                bends.len()
            ));
        }
        if thickness <= 0.0 || width <= 0.0 {
            return Err("thickness and width must be positive".into());
        }
        let k = k_factor(bends.first().map_or(1.0, |b| b.inside_radius / thickness));
        Ok(Self { thickness, k_factor: k, width, flanges, bends })
    }

    pub fn with_k_factor(mut self, k: f64) -> Self {
        self.k_factor = k;
        self
    }

    /// The three bend quantities for bend `i`.
    pub fn numbers(&self, i: usize) -> Result<BendNumbers, String> {
        let b = *self.bends.get(i).ok_or("no such bend")?;
        // The sign of the angle only says which way the flange turns; the
        // material cost depends on its magnitude.
        let angle = b.angle.abs();
        if angle <= 0.0 || angle >= PI {
            return Err(format!(
                "bend {i}: the angle turned through must be between 0 and pi, got {}",
                b.angle
            ));
        }
        if b.inside_radius <= 0.0 {
            return Err(format!("bend {i}: the inside radius must be positive"));
        }
        let neutral_radius = b.inside_radius + self.k_factor * self.thickness;
        let allowance = angle * neutral_radius;
        let outside_setback = (b.inside_radius + self.thickness) * (0.5 * angle).tan();
        Ok(BendNumbers {
            allowance,
            outside_setback,
            deduction: 2.0 * outside_setback - allowance,
            neutral_radius,
        })
    }

    /// Straight length of flange `i` along the neutral axis: the apex
    /// dimension less the setback of each bend it touches.
    pub fn straight(&self, i: usize) -> Result<f64, String> {
        let apex = *self.flanges.get(i).ok_or("no such flange")?;
        let mut length = apex;
        if i > 0 {
            length -= self.numbers(i - 1)?.outside_setback;
        }
        if i < self.bends.len() {
            length -= self.numbers(i)?.outside_setback;
        }
        if length <= 0.0 {
            return Err(format!(
                "flange {i} is shorter than its setbacks: a flange of {apex} cannot be formed with \
                 these radii (the shop minimum is about R + T per bend)"
            ));
        }
        Ok(length)
    }

    /// Flat blank length from the apex dimensions:
    /// `sum(flange) - sum(bend deduction)`.
    pub fn flat_length(&self) -> Result<f64, String> {
        let mut total: f64 = self.flanges.iter().sum();
        for i in 0..self.bends.len() {
            total -= self.numbers(i)?.deduction;
        }
        Ok(total)
    }

    /// The same length by the other route:
    /// `sum(straight) + sum(bend allowance)`.
    pub fn flat_length_from_geometry(&self) -> Result<f64, String> {
        let mut total = 0.0;
        for i in 0..self.flanges.len() {
            total += self.straight(i)?;
        }
        for i in 0..self.bends.len() {
            total += self.numbers(i)?.allowance;
        }
        Ok(total)
    }

    /// The neutral-axis path of the formed section, as a polyline with
    /// `per_bend` points through each arc.
    pub fn neutral_axis(&self, per_bend: usize) -> Result<Vec<V2>, String> {
        let per_bend = per_bend.max(2);
        let mut pts = vec![[0.0, 0.0]];
        let mut p = [0.0, 0.0];
        let mut dir = 0.0f64; // heading, radians
        for i in 0..self.flanges.len() {
            let s = self.straight(i)?;
            p = [p[0] + s * dir.cos(), p[1] + s * dir.sin()];
            pts.push(p);
            if i < self.bends.len() {
                let n = self.numbers(i)?;
                let b = self.bends[i];
                let sign = b.angle.signum();
                // Centre of the arc, to the left for a positive bend.
                let centre = [
                    p[0] - sign * n.neutral_radius * dir.sin(),
                    p[1] + sign * n.neutral_radius * dir.cos(),
                ];
                let start = (p[1] - centre[1]).atan2(p[0] - centre[0]);
                for k in 1..=per_bend {
                    let t = start + sign * b.angle.abs() * k as f64 / per_bend as f64;
                    p = [
                        centre[0] + n.neutral_radius * t.cos(),
                        centre[1] + n.neutral_radius * t.sin(),
                    ];
                    pts.push(p);
                }
                dir += b.angle;
            }
        }
        Ok(pts)
    }

    /// Length of the generated neutral-axis path, measured from the
    /// geometry (arcs exactly, by summing `R theta` rather than chords).
    pub fn neutral_axis_length(&self) -> Result<f64, String> {
        self.flat_length_from_geometry()
    }

    /// The formed cross-section: the band between the inside surface, at
    /// `K T` inside the neutral axis, and the outside, at `(1 - K) T`.
    pub fn section(&self, per_bend: usize) -> Result<Polygon, String> {
        let axis = self.neutral_axis(per_bend)?;
        if axis.len() < 2 {
            return Err("the section has no length".into());
        }
        let normal_at = |i: usize| -> V2 {
            // Left normal of the local tangent.
            let (a, b) = if i == 0 {
                (axis[0], axis[1])
            } else if i == axis.len() - 1 {
                (axis[axis.len() - 2], axis[axis.len() - 1])
            } else {
                (axis[i - 1], axis[i + 1])
            };
            let t = [b[0] - a[0], b[1] - a[1]];
            let n = (t[0] * t[0] + t[1] * t[1]).sqrt();
            if n > 0.0 {
                [-t[1] / n, t[0] / n]
            } else {
                [0.0, 1.0]
            }
        };
        let inside = self.k_factor * self.thickness;
        let outside = (1.0 - self.k_factor) * self.thickness;
        let mut pts: Vec<V2> = Vec::with_capacity(2 * axis.len());
        // The inside surface is on the left of a positive (left-turning)
        // bend, so it is offset along the left normal.
        for (i, p) in axis.iter().enumerate() {
            let n = normal_at(i);
            pts.push([p[0] + n[0] * inside, p[1] + n[1] * inside]);
        }
        for (i, p) in axis.iter().enumerate().rev() {
            let n = normal_at(i);
            pts.push([p[0] - n[0] * outside, p[1] - n[1] * outside]);
        }
        Ok(Polygon::new(pts))
    }

    pub fn region(&self, per_bend: usize) -> Result<Region2, String> {
        Ok(Region2::Polygon(self.section(per_bend)?))
    }

    /// The formed part.
    pub fn solid(&self, per_bend: usize) -> Result<Solid, String> {
        Ok(Solid::extrude(self.region(per_bend)?, 0.0, self.width))
    }

    /// The flat blank: a rectangle of the developed length by the width.
    pub fn flat_pattern(&self) -> Result<Region2, String> {
        let l = self.flat_length()?;
        Ok(Region2::rect([0.5 * l, 0.5 * self.width], [0.5 * l, 0.5 * self.width]))
    }

    /// The blank as a solid of the material thickness.
    pub fn blank(&self) -> Result<Solid, String> {
        Ok(Solid::extrude(self.flat_pattern()?, 0.0, self.thickness))
    }

    /// Blank volume.
    pub fn blank_volume(&self) -> Result<f64, String> {
        Ok(self.flat_length()? * self.width * self.thickness)
    }

    /// Formed-section area times width, from the generated polygon.
    pub fn formed_volume(&self, per_bend: usize) -> Result<f64, String> {
        Ok(self.section(per_bend)?.area().abs() * self.width)
    }

    /// Relative volume change from blank to formed part.
    ///
    /// Zero only at `K = 0.5`: the neutral-axis model stretches the
    /// outside and compresses the inside, and the two cancel exactly only
    /// when the neutral axis is at mid-thickness. The number is a property
    /// of the model, not of the metal, which really does conserve volume.
    pub fn volume_change(&self, per_bend: usize) -> Result<f64, String> {
        let blank = self.blank_volume()?;
        Ok((self.formed_volume(per_bend)? - blank) / blank)
    }

    /// Shop minimum flange length, `R + T` per bend it touches: shorter
    /// than this and the flange is inside the die.
    pub fn minimum_flange(&self, i: usize) -> Result<f64, String> {
        let mut m = 0.0;
        if i > 0 {
            m += self.bends[i - 1].inside_radius + self.thickness;
        }
        if i < self.bends.len() {
            m += self.bends[i].inside_radius + self.thickness;
        }
        Ok(m)
    }

    /// Overall size of the formed part, from the generated section.
    pub fn formed_bbox(&self, per_bend: usize) -> Result<(V2, V2), String> {
        Ok(self.section(per_bend)?.bbox())
    }
}
