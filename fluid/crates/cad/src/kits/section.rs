//! Structural sections: the profiles a beam is made from, and the section
//! properties a beam model needs.
//!
//! A beam element is a line plus four numbers — area, two second moments
//! and a torsion constant. This kit draws the cross-section as a
//! [`Region2`] and computes those numbers from the drawing, so the solid
//! that gets meshed and the stiffness that gets solved come from the same
//! geometry. [`Properties`] is what `architected`'s frame and beam
//! elements take.
//!
//! ```text
//!   Section::i_beam { d, bf, tw, tf } -> Region2 -> Solid (extrude)
//!                                     -> Properties { a, ixx, iyy, sx, rx, j, zx }
//!                                     -> architected::frame::Member
//! ```
//!
//! All lengths are metres, consistently; the shapes are parametric, so a
//! standard size is entered from its own table.
//!
//! # What is exact, and what is not
//!
//! * **Area, centroid and the second moments are exact** for the drawn
//!   polygon, by Green's theorem ([`crate::Polygon::second_moments`]). The
//!   tests check them against the textbook composite-rectangle formulas,
//!   which is an independent computation path, to 1e-15.
//! * **The shapes are idealised**: square corners, no root fillets, no
//!   rolling tolerance. A rolled wide-flange shape's tabulated `Ix`
//!   includes its fillets and so is a per cent or two above the idealised
//!   value; a tabulated area likewise. Nothing here reproduces a mill
//!   table — the dimensions are the caller's input.
//! * **The torsion constant of an open section is approximate.** For thin
//!   open shapes `J ~ sum(b t^3 / 3)` (Timoshenko & Goodier, *Theory of
//!   Elasticity*, 3rd ed., §109), which ignores the junction effect and is
//!   low by 10-20 % for a wide flange. It is marked as an estimate in the
//!   returned struct. The closed sections' `J` is the exact thin-wall
//!   Bredt formula where the wall is thin, and that too is stated.
//! * **The plastic modulus** assumes a doubly symmetric shape bending
//!   about its own axis of symmetry.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::sketch::{Polygon, Region2, V2};
use crate::solid::Solid;

/// Section properties about the centroid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Properties {
    /// Area.
    pub area: f64,
    pub centroid: V2,
    /// Second moment about the centroidal `x` axis (bending in the `y`
    /// direction).
    pub ixx: f64,
    /// Second moment about the centroidal `y` axis.
    pub iyy: f64,
    /// Product of inertia; zero for a symmetric shape.
    pub ixy: f64,
    /// Elastic section modulus about `x`, `Ixx / y_max`.
    pub sx: f64,
    pub sy: f64,
    /// Radii of gyration.
    pub rx: f64,
    pub ry: f64,
    /// Torsion constant. `exact_j` says whether it is the thin-wall
    /// formula's own value or an estimate for an open shape.
    pub j: f64,
    pub exact_j: bool,
    /// Perimeter of the drawn outline.
    pub perimeter: f64,
}

/// A parametric structural shape. Every dimension is an outside dimension
/// except the wall thicknesses.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
#[serde(tag = "shape", rename_all = "snake_case")]
pub enum Section {
    /// Solid rectangle, `width` by `height`.
    Rectangle { width: f64, height: f64 },
    /// Solid round bar.
    Round { diameter: f64 },
    /// Wide flange / I-beam: overall depth `d`, flange width `bf`, web
    /// thickness `tw`, flange thickness `tf`.
    IBeam { d: f64, bf: f64, tw: f64, tf: f64 },
    /// Channel, opening toward `+x`.
    Channel { d: f64, bf: f64, tw: f64, tf: f64 },
    /// Equal or unequal angle, legs along `+x` and `+y`.
    Angle { leg_x: f64, leg_y: f64, thickness: f64 },
    /// Tee: flange on top.
    Tee { d: f64, bf: f64, tw: f64, tf: f64 },
    /// Rectangular hollow section.
    Box { width: f64, height: f64, thickness: f64 },
    /// Circular hollow section.
    Tube { diameter: f64, thickness: f64 },
}

impl Section {
    /// The outline as a sketch region, centred on its own centroid.
    pub fn region(&self) -> Result<Region2, String> {
        let poly = self.polygon()?;
        Ok(match *self {
            Section::Round { diameter } => Region2::circle([0.0, 0.0], 0.5 * diameter),
            Section::Tube { diameter, thickness } => {
                if thickness <= 0.0 || 2.0 * thickness >= diameter {
                    return Err("a tube's wall must be positive and thinner than its radius".into());
                }
                Region2::ring([0.0, 0.0], 0.5 * diameter - thickness, 0.5 * diameter)
            }
            _ => Region2::Polygon(poly),
        })
    }

    /// The outline as a polygon (with a hole for a box section). Round
    /// shapes are approximated here; [`Section::region`] keeps them exact.
    pub fn polygon(&self) -> Result<Polygon, String> {
        let positive = |v: f64, name: &str| -> Result<(), String> {
            if v > 0.0 {
                Ok(())
            } else {
                Err(format!("{name} must be positive"))
            }
        };
        let poly = match *self {
            Section::Rectangle { width, height } => {
                positive(width, "width")?;
                positive(height, "height")?;
                rect(width, height)
            }
            Section::Round { diameter } => {
                positive(diameter, "diameter")?;
                Polygon::new(
                    (0..256)
                        .map(|i| {
                            let a = 2.0 * PI * i as f64 / 256.0;
                            [0.5 * diameter * a.cos(), 0.5 * diameter * a.sin()]
                        })
                        .collect(),
                )
            }
            Section::IBeam { d, bf, tw, tf } => {
                check_web(d, bf, tw, tf)?;
                let (hw, ht) = (0.5 * tw, 0.5 * d - tf);
                Polygon::new(vec![
                    [-0.5 * bf, -0.5 * d],
                    [0.5 * bf, -0.5 * d],
                    [0.5 * bf, -ht],
                    [hw, -ht],
                    [hw, ht],
                    [0.5 * bf, ht],
                    [0.5 * bf, 0.5 * d],
                    [-0.5 * bf, 0.5 * d],
                    [-0.5 * bf, ht],
                    [-hw, ht],
                    [-hw, -ht],
                    [-0.5 * bf, -ht],
                ])
            }
            Section::Channel { d, bf, tw, tf } => {
                check_web(d, bf, tw, tf)?;
                let ht = 0.5 * d - tf;
                Polygon::new(vec![
                    [0.0, -0.5 * d],
                    [bf, -0.5 * d],
                    [bf, -ht],
                    [tw, -ht],
                    [tw, ht],
                    [bf, ht],
                    [bf, 0.5 * d],
                    [0.0, 0.5 * d],
                ])
            }
            Section::Angle { leg_x, leg_y, thickness } => {
                positive(thickness, "thickness")?;
                if thickness >= leg_x || thickness >= leg_y {
                    return Err("an angle's thickness must be less than both legs".into());
                }
                Polygon::new(vec![
                    [0.0, 0.0],
                    [leg_x, 0.0],
                    [leg_x, thickness],
                    [thickness, thickness],
                    [thickness, leg_y],
                    [0.0, leg_y],
                ])
            }
            Section::Tee { d, bf, tw, tf } => {
                check_web(d, bf, tw, tf)?;
                let ht = d - tf;
                Polygon::new(vec![
                    [-0.5 * tw, 0.0],
                    [0.5 * tw, 0.0],
                    [0.5 * tw, ht],
                    [0.5 * bf, ht],
                    [0.5 * bf, d],
                    [-0.5 * bf, d],
                    [-0.5 * bf, ht],
                    [-0.5 * tw, ht],
                ])
            }
            Section::Box { width, height, thickness } => {
                positive(thickness, "thickness")?;
                if 2.0 * thickness >= width.min(height) {
                    return Err("a box section's walls must not meet".into());
                }
                let mut p = rect(width, height);
                let (iw, ih) = (width - 2.0 * thickness, height - 2.0 * thickness);
                p.loops.push(rect(iw, ih).loops.remove(0));
                p
            }
            Section::Tube { diameter, thickness } => {
                positive(thickness, "thickness")?;
                if 2.0 * thickness >= diameter {
                    return Err("a tube's walls must not meet".into());
                }
                let ring = |r: f64| -> Vec<V2> {
                    (0..256)
                        .map(|i| {
                            let a = 2.0 * PI * i as f64 / 256.0;
                            [r * a.cos(), r * a.sin()]
                        })
                        .collect()
                };
                let mut p = Polygon::new(ring(0.5 * diameter));
                p.loops.push(ring(0.5 * diameter - thickness));
                p
            }
        };
        Ok(poly)
    }

    /// Section properties from the drawn outline.
    ///
    /// Area, centroid and the second moments come from Green's theorem on
    /// the polygon and are exact for it. Closed-form areas and moments are
    /// used for the two round shapes, where a polygon would only
    /// approximate them.
    pub fn properties(&self) -> Result<Properties, String> {
        let poly = self.polygon()?;
        let (area, centroid, ixx, iyy, ixy) = match *self {
            Section::Round { diameter } => {
                let r = 0.5 * diameter;
                let i = PI * r.powi(4) / 4.0;
                (PI * r * r, [0.0, 0.0], i, i, 0.0)
            }
            Section::Tube { diameter, thickness } => {
                let (ro, ri) = (0.5 * diameter, 0.5 * diameter - thickness);
                let i = PI * (ro.powi(4) - ri.powi(4)) / 4.0;
                (PI * (ro * ro - ri * ri), [0.0, 0.0], i, i, 0.0)
            }
            _ => {
                let c = poly.centroid();
                let (ixx, iyy, ixy) = poly.second_moments();
                (poly.area(), c, ixx, iyy, ixy)
            }
        };
        if area <= 0.0 {
            return Err("the section has no area".into());
        }
        // Extreme fibre distances from the centroid.
        let (lo, hi) = poly.bbox();
        let y_max = (hi[1] - centroid[1]).abs().max((centroid[1] - lo[1]).abs());
        let x_max = (hi[0] - centroid[0]).abs().max((centroid[0] - lo[0]).abs());
        let (j, exact_j) = self.torsion_constant()?;
        Ok(Properties {
            area,
            centroid,
            ixx,
            iyy,
            ixy,
            sx: if y_max > 0.0 { ixx / y_max } else { 0.0 },
            sy: if x_max > 0.0 { iyy / x_max } else { 0.0 },
            rx: (ixx / area).sqrt(),
            ry: (iyy / area).sqrt(),
            j,
            exact_j,
            perimeter: poly.perimeter(),
        })
    }

    /// Torsion constant, and whether it is a closed form or an estimate.
    fn torsion_constant(&self) -> Result<(f64, bool), String> {
        Ok(match *self {
            Section::Round { diameter } => (PI * (0.5 * diameter).powi(4) / 2.0, true),
            Section::Tube { diameter, thickness } => {
                let (ro, ri) = (0.5 * diameter, 0.5 * diameter - thickness);
                (PI * (ro.powi(4) - ri.powi(4)) / 2.0, true)
            }
            Section::Rectangle { width, height } => {
                // Saint-Venant for a rectangle, to within 0.5 % for any
                // aspect ratio (Timoshenko & Goodier §109, the series
                // truncated at one term plus the thin-strip limit).
                let (a, b) = (0.5 * width.max(height), 0.5 * width.min(height));
                let k = 1.0 - 0.63 * (b / a) * (1.0 - (b / a).powi(4) / 12.0);
                (k * (2.0 * a) * (2.0 * b).powi(3) / 3.0, false)
            }
            Section::Box { width, height, thickness } => {
                // Bredt's thin-wall formula: J = 4 A_m^2 / integral(ds/t).
                let (bm, hm) = (width - thickness, height - thickness);
                let am = bm * hm;
                let path = 2.0 * (bm + hm) / thickness;
                (4.0 * am * am / path, true)
            }
            Section::IBeam { d, bf, tw, tf } => {
                let web = (d - 2.0 * tf) * tw.powi(3);
                ((web + 2.0 * bf * tf.powi(3)) / 3.0, false)
            }
            Section::Channel { d, bf, tw, tf } => {
                let web = (d - 2.0 * tf) * tw.powi(3);
                ((web + 2.0 * bf * tf.powi(3)) / 3.0, false)
            }
            Section::Tee { d, bf, tw, tf } => (((d - tf) * tw.powi(3) + bf * tf.powi(3)) / 3.0, false),
            Section::Angle { leg_x, leg_y, thickness } => {
                ((leg_x + leg_y - thickness) * thickness.powi(3) / 3.0, false)
            }
        })
    }

    /// Plastic section modulus about the centroidal `x` axis, for a shape
    /// symmetric about that axis.
    pub fn plastic_modulus_x(&self) -> Result<f64, String> {
        Ok(match *self {
            Section::Rectangle { width, height } => width * height * height / 4.0,
            Section::Round { diameter } => (0.5 * diameter).powi(3) * 4.0 / 3.0,
            Section::IBeam { d, bf, tw, tf } => {
                let hw = d - 2.0 * tf;
                bf * tf * (d - tf) + tw * hw * hw / 4.0
            }
            Section::Box { width, height, thickness } => {
                let (iw, ih) = (width - 2.0 * thickness, height - 2.0 * thickness);
                width * height * height / 4.0 - iw * ih * ih / 4.0
            }
            Section::Tube { diameter, thickness } => {
                let (ro, ri) = (0.5 * diameter, 0.5 * diameter - thickness);
                4.0 / 3.0 * (ro.powi(3) - ri.powi(3))
            }
            _ => return Err("the plastic modulus here assumes symmetry about the bending axis".into()),
        })
    }

    /// The section extruded into a beam of `length` along `+z`.
    pub fn beam(&self, length: f64) -> Result<Solid, String> {
        Ok(Solid::extrude(self.region()?, 0.0, length))
    }
}

fn rect(width: f64, height: f64) -> Polygon {
    Polygon::new(vec![
        [-0.5 * width, -0.5 * height],
        [0.5 * width, -0.5 * height],
        [0.5 * width, 0.5 * height],
        [-0.5 * width, 0.5 * height],
    ])
}

fn check_web(d: f64, bf: f64, tw: f64, tf: f64) -> Result<(), String> {
    if d <= 0.0 || bf <= 0.0 || tw <= 0.0 || tf <= 0.0 {
        return Err("every dimension must be positive".into());
    }
    if 2.0 * tf >= d {
        return Err("the flanges would meet: 2 tf >= d".into());
    }
    if tw >= bf {
        return Err("the web is wider than the flange".into());
    }
    Ok(())
}

/// Tip deflection of a cantilever of length `L` under a tip load `P`:
/// `P L^3 / (3 E I)`, the Euler-Bernoulli closed form.
pub fn cantilever_tip_deflection(load: f64, length: f64, e: f64, inertia: f64) -> f64 {
    load * length.powi(3) / (3.0 * e * inertia)
}

/// Euler buckling load of a pinned column, `pi^2 E I / L^2`.
pub fn euler_buckling(e: f64, inertia: f64, length: f64) -> f64 {
    PI * PI * e * inertia / (length * length)
}
