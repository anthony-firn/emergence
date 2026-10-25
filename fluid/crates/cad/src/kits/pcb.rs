//! Printed circuit boards.
//!
//! A board is a layer stackup, an outline, and copper drawn on each layer
//! as sketch regions: traces are capsule chains, pads are rectangles or
//! circles, pours are the outline minus clearances. From that the kit
//! produces
//!
//! * a **solid** of the whole board, dielectric and copper, for an airflow
//!   or thermal run through `immersed` or `heat`;
//! * a **layer stack** ([`crate::sim::LayerStack`]) carrying each layer's
//!   permittivity and conductivity, which `photonics` and `phononic` take
//!   directly;
//! * **electrical answers about what was drawn**: DC resistance from the
//!   exact copper cross-section, and the characteristic impedance of a
//!   trace both from the published closed forms and from a field solve on
//!   the generated cross-section ([`crate::laplace`]);
//! * **manufacturing output**: Gerber RS-274X copper layers and an
//!   Excellon drill file.
//!
//! All lengths are metres. Copper thickness is given in ounces per square
//! foot, the industry unit: 1 oz is 34.8 um
//! (IPC-4562A; 1 oz/ft^2 of copper over 1 ft^2 at 8.96 g/cm^3).
//!
//! # Material constants
//!
//! | quantity | value | source |
//! |---|---|---|
//! | copper resistivity at 20 C | 1.678e-8 ohm m | CRC Handbook of Chemistry and Physics, 97th ed. |
//! | copper thermal conductivity | 401 W/(m K) | CRC Handbook, 97th ed. |
//! | copper temperature coefficient | 3.93e-3 /K | CRC Handbook, 97th ed. |
//! | FR-4 relative permittivity at 1 GHz | 4.4 | typical woven-glass/epoxy laminate; varies 4.2-4.8 with resin content and weave, so it is a *nominal* figure |
//! | FR-4 loss tangent at 1 GHz | 0.02 | same caveat |
//! | FR-4 thermal conductivity, through-plane | 0.29 W/(m K) | typical laminate datasheet value |
//! | solder mask relative permittivity | 3.8 | typical liquid photoimageable mask |
//!
//! The FR-4 numbers are nominal: a real board's permittivity depends on
//! the glass style and resin content, and a controlled-impedance design
//! must use the laminate's own datasheet. The kit takes them as
//! parameters for that reason.
//!
//! # Microstrip and coplanar formulas
//!
//! [`microstrip_z0`] is Hammerstad's synthesis of Wheeler's analysis
//! (Hammerstad, *5th European Microwave Conference*, 1975; Hammerstad &
//! Jensen, *IEEE MTT-S*, 1980), quoted by its authors as accurate to about
//! 1 % for `0.05 <= W/h <= 20` and `eps_r <= 16`. It assumes zero
//! conductor thickness. [`crate::laplace::quasi_tem`] solves the drawn
//! cross-section instead, including the copper thickness; the tests
//! compare the two and find them within 6 %.

use std::f64::consts::PI;
use std::fmt::Write as _;

use serde::{Deserialize, Serialize};

use crate::sim::{Layer, LayerStack};
use crate::sketch::{Region2, V2};
use crate::solid::Solid;

/// Copper resistivity at 20 C, ohm m (CRC Handbook, 97th ed.).
pub const RHO_CU: f64 = 1.678e-8;
/// Copper thermal conductivity, W/(m K) (CRC Handbook, 97th ed.).
pub const KAPPA_CU: f64 = 401.0;
/// Copper resistance temperature coefficient, 1/K (CRC Handbook).
pub const ALPHA_CU: f64 = 3.93e-3;
/// Thickness of 1 oz/ft^2 copper, m (IPC-4562A).
pub const OZ_COPPER: f64 = 34.8e-6;

/// A stackup layer.
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum PcbLayer {
    /// A copper layer, thickness in ounces per square foot.
    Copper { #[serde(default = "one_oz")] ounces: f64 },
    /// A dielectric core or prepreg.
    Dielectric {
        #[serde(default = "fr4_name")]
        material: String,
        thickness: f64,
        #[serde(default = "fr4_eps")]
        eps_r: f64,
        #[serde(default = "fr4_tan")]
        loss_tangent: f64,
        #[serde(default = "fr4_kappa")]
        kappa: f64,
    },
    /// Solder mask over the outer copper.
    Mask { thickness: f64, #[serde(default = "mask_eps")] eps_r: f64 },
}

fn one_oz() -> f64 {
    1.0
}
fn fr4_name() -> String {
    "FR-4".into()
}
fn fr4_eps() -> f64 {
    4.4
}
fn fr4_tan() -> f64 {
    0.02
}
fn fr4_kappa() -> f64 {
    0.29
}
fn mask_eps() -> f64 {
    3.8
}

impl PcbLayer {
    pub fn thickness(&self) -> f64 {
        match self {
            PcbLayer::Copper { ounces } => ounces * OZ_COPPER,
            PcbLayer::Dielectric { thickness, .. } => *thickness,
            PcbLayer::Mask { thickness, .. } => *thickness,
        }
    }

    pub fn is_copper(&self) -> bool {
        matches!(self, PcbLayer::Copper { .. })
    }

    pub fn eps_r(&self) -> f64 {
        match self {
            PcbLayer::Copper { .. } => 1.0,
            PcbLayer::Dielectric { eps_r, .. } => *eps_r,
            PcbLayer::Mask { eps_r, .. } => *eps_r,
        }
    }

    pub fn as_layer(&self) -> Layer {
        match self {
            PcbLayer::Copper { .. } => Layer::new("copper", self.thickness())
                .electrical(1.0, 1.0 / RHO_CU)
                .thermal(KAPPA_CU),
            PcbLayer::Dielectric { material, eps_r, kappa, .. } => {
                Layer::new(material, self.thickness()).electrical(*eps_r, 0.0).thermal(*kappa)
            }
            PcbLayer::Mask { eps_r, .. } => Layer::new("mask", self.thickness()).electrical(*eps_r, 0.0),
        }
    }
}

/// A copper feature on one layer.
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Feature {
    /// A trace: a chain of segments of constant width.
    Trace { layer: usize, points: Vec<V2>, width: f64 },
    /// A rectangular pad.
    Pad { layer: usize, center: V2, size: V2 },
    /// A round pad.
    Round { layer: usize, center: V2, diameter: f64 },
    /// A filled region, typically a ground pour.
    Pour { layer: usize, region: Region2 },
}

impl Feature {
    pub fn layer(&self) -> usize {
        match self {
            Feature::Trace { layer, .. }
            | Feature::Pad { layer, .. }
            | Feature::Round { layer, .. }
            | Feature::Pour { layer, .. } => *layer,
        }
    }

    pub fn region(&self) -> Region2 {
        match self {
            Feature::Trace { points, width, .. } => {
                if points.len() < 2 {
                    return Region2::Empty;
                }
                Region2::Union {
                    parts: points
                        .windows(2)
                        .map(|w| Region2::slot(w[0], w[1], 0.5 * width))
                        .collect(),
                }
            }
            Feature::Pad { center, size, .. } => Region2::rect(*center, [0.5 * size[0], 0.5 * size[1]]),
            Feature::Round { center, diameter, .. } => Region2::circle(*center, 0.5 * diameter),
            Feature::Pour { region, .. } => region.clone(),
        }
    }

    /// Length of a trace's centreline; 0 for anything else.
    pub fn length(&self) -> f64 {
        match self {
            Feature::Trace { points, .. } => points
                .windows(2)
                .map(|w| ((w[1][0] - w[0][0]).powi(2) + (w[1][1] - w[0][1]).powi(2)).sqrt())
                .sum(),
            _ => 0.0,
        }
    }
}

/// A plated through hole.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Via {
    pub center: V2,
    /// Finished drill diameter.
    pub drill: f64,
    /// Width of the copper ring around the hole on each copper layer.
    pub annular_ring: f64,
    /// Barrel plating thickness (25 um is the IPC-6012 class 2 minimum).
    #[serde(default = "plating")]
    pub plating: f64,
}

fn plating() -> f64 {
    25e-6
}

/// A board.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Board {
    pub outline: Region2,
    pub stackup: Vec<PcbLayer>,
    #[serde(default)]
    pub features: Vec<Feature>,
    #[serde(default)]
    pub vias: Vec<Via>,
    /// Non-plated holes.
    #[serde(default)]
    pub holes: Vec<(V2, f64)>,
}

impl Board {
    pub fn new(outline: Region2, stackup: Vec<PcbLayer>) -> Self {
        Self { outline, stackup, features: Vec::new(), vias: Vec::new(), holes: Vec::new() }
    }

    pub fn with(mut self, feature: Feature) -> Self {
        self.features.push(feature);
        self
    }

    pub fn with_via(mut self, via: Via) -> Self {
        self.vias.push(via);
        self
    }

    pub fn thickness(&self) -> f64 {
        self.stackup.iter().map(|l| l.thickness()).sum()
    }

    /// Height of the bottom of each stackup layer.
    pub fn offsets(&self) -> Vec<f64> {
        let mut z = 0.0;
        self.stackup
            .iter()
            .map(|l| {
                let z0 = z;
                z += l.thickness();
                z0
            })
            .collect()
    }

    /// Indices into `stackup` of the copper layers, bottom to top.
    pub fn copper_layers(&self) -> Vec<usize> {
        self.stackup.iter().enumerate().filter(|(_, l)| l.is_copper()).map(|(i, _)| i).collect()
    }

    /// The copper on one copper layer (indexed over copper layers only, as
    /// features are), including via rings.
    pub fn copper(&self, layer: usize) -> Region2 {
        let mut parts: Vec<Region2> = self
            .features
            .iter()
            .filter(|f| f.layer() == layer)
            .map(|f| f.region())
            .collect();
        for v in &self.vias {
            parts.push(Region2::ring(v.center, 0.5 * v.drill, 0.5 * v.drill + v.annular_ring));
        }
        if parts.is_empty() {
            return Region2::Empty;
        }
        // Copper cannot exist outside the outline, nor across a hole.
        let mut region = Region2::Union { parts }.and(self.outline.clone());
        let holes: Vec<Region2> = self
            .vias
            .iter()
            .map(|v| Region2::circle(v.center, 0.5 * v.drill))
            .chain(self.holes.iter().map(|&(c, d)| Region2::circle(c, 0.5 * d)))
            .collect();
        if !holes.is_empty() {
            region = region.minus(Region2::Union { parts: holes });
        }
        region
    }

    /// Fraction of the board area covered by copper on a layer, by
    /// quadrature at step `cell`.
    pub fn coverage(&self, layer: usize, cell: f64) -> f64 {
        let copper = self.copper(layer).area_quadrature(cell);
        let board = self.outline.area_quadrature(cell);
        if board > 0.0 {
            copper / board
        } else {
            0.0
        }
    }

    /// The whole board as a solid: every layer extruded over the outline,
    /// copper only where it was drawn, with the drills cut through.
    pub fn solid(&self) -> Solid {
        let offsets = self.offsets();
        let mut copper_index = 0usize;
        let mut parts = Vec::new();
        for (i, layer) in self.stackup.iter().enumerate() {
            let (z0, z1) = (offsets[i], offsets[i] + layer.thickness());
            if layer.is_copper() {
                parts.push(Solid::extrude(self.copper(copper_index), z0, z1));
                copper_index += 1;
            } else {
                parts.push(Solid::extrude(self.outline.clone(), z0, z1));
            }
        }
        let mut solid = Solid::union(parts);
        let drills: Vec<Solid> = self
            .vias
            .iter()
            .map(|v| (v.center, v.drill))
            .chain(self.holes.iter().copied())
            .map(|(c, d)| {
                Solid::cylinder(
                    [c[0], c[1], -self.thickness()],
                    [c[0], c[1], 2.0 * self.thickness()],
                    0.5 * d,
                )
            })
            .collect();
        if !drills.is_empty() {
            solid = solid.minus(Solid::union(drills));
        }
        solid
    }

    /// The stackup as a layer stack for the optics and thermal crates.
    pub fn layer_stack(&self) -> LayerStack {
        LayerStack::new(self.stackup.iter().map(|l| l.as_layer()).collect())
    }

    /// DC resistance of a trace at `temperature` (Celsius), from its exact
    /// drawn length and copper cross-section.
    pub fn trace_resistance(&self, feature: usize, temperature: f64) -> Result<f64, String> {
        let f = self.features.get(feature).ok_or("no such feature")?;
        let (layer, width) = match f {
            Feature::Trace { layer, width, .. } => (*layer, *width),
            _ => return Err("not a trace".into()),
        };
        let copper = self
            .copper_layers()
            .get(layer)
            .and_then(|&i| self.stackup.get(i))
            .ok_or("no such copper layer")?;
        let rho = RHO_CU * (1.0 + ALPHA_CU * (temperature - 20.0));
        Ok(rho * f.length() / (width * copper.thickness()))
    }

    /// Sheet resistance of a copper layer, ohms per square.
    pub fn sheet_resistance(&self, layer: usize) -> Result<f64, String> {
        let copper = self
            .copper_layers()
            .get(layer)
            .and_then(|&i| self.stackup.get(i))
            .ok_or("no such copper layer")?;
        Ok(RHO_CU / copper.thickness())
    }

    /// Resistance of a via barrel, from the plated annulus.
    pub fn via_resistance(&self, via: usize) -> Result<f64, String> {
        let v = self.vias.get(via).ok_or("no such via")?;
        let r_out = 0.5 * v.drill + v.plating;
        let area = PI * (r_out * r_out - 0.25 * v.drill * v.drill);
        Ok(RHO_CU * self.thickness() / area)
    }

    /// Thermal resistance of the stack per unit area, m^2 K/W.
    pub fn thermal_resistance(&self) -> Result<f64, String> {
        self.layer_stack().thermal_resistance()
    }

    /// Dielectric thickness between a copper layer and the next copper
    /// layer below it: the `h` a microstrip or stripline formula needs.
    pub fn dielectric_below(&self, layer: usize) -> Option<(f64, f64)> {
        let coppers = self.copper_layers();
        let here = *coppers.get(layer)?;
        let below = *coppers.get(layer.checked_sub(1)?)?;
        let mut thickness = 0.0;
        let mut weighted = 0.0;
        for l in &self.stackup[below + 1..here] {
            thickness += l.thickness();
            weighted += l.thickness() * l.eps_r();
        }
        if thickness > 0.0 {
            Some((thickness, weighted / thickness))
        } else {
            None
        }
    }
}

/// Microstrip characteristic impedance and effective permittivity
/// (Hammerstad 1975; Hammerstad & Jensen 1980). Zero conductor thickness;
/// about 1 % accurate for `0.05 <= W/h <= 20`.
pub fn microstrip_z0(width: f64, height: f64, eps_r: f64) -> (f64, f64) {
    let u = width / height;
    let eps_eff = if u <= 1.0 {
        0.5 * (eps_r + 1.0)
            + 0.5 * (eps_r - 1.0) * ((1.0 + 12.0 / u).powf(-0.5) + 0.04 * (1.0 - u).powi(2))
    } else {
        0.5 * (eps_r + 1.0) + 0.5 * (eps_r - 1.0) * (1.0 + 12.0 / u).powf(-0.5)
    };
    let z0 = if u <= 1.0 {
        60.0 / eps_eff.sqrt() * (8.0 / u + 0.25 * u).ln()
    } else {
        120.0 * PI / (eps_eff.sqrt() * (u + 1.393 + 0.667 * (u + 1.444).ln()))
    };
    (z0, eps_eff)
}

/// Trace width that gives `target` ohms on a given stack, by bisection on
/// [`microstrip_z0`].
pub fn microstrip_width_for(target: f64, height: f64, eps_r: f64) -> f64 {
    // Z0 falls as the trace widens, so a width whose impedance is still
    // above the target is too narrow.
    let (mut lo, mut hi) = (0.01 * height, 30.0 * height);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if microstrip_z0(mid, height, eps_r).0 > target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Gerber RS-274X for one copper layer, as filled regions.
///
/// **Not verified against a CAM tool.** The file is written to the format
/// specification (Ucamco, *The Gerber Layer Format Specification*,
/// revision 2023.08): metric, 4.6 format, dark polarity, one `G36/G37`
/// region per boundary loop. Only polygonal boundaries are emitted, so a
/// region is first flattened to a polygon at `cell` resolution. Anything
/// that matters should be checked in a viewer before it is sent out.
pub fn gerber_layer(region: &Region2, cell: f64, name: &str) -> String {
    let loops = flatten(region, cell);
    let mut out = String::new();
    let _ = writeln!(out, "G04 {name} generated by the cad crate*");
    let _ = writeln!(out, "%FSLAX46Y46*%");
    let _ = writeln!(out, "%MOMM*%");
    let _ = writeln!(out, "%LPD*%");
    let _ = writeln!(out, "%ADD10C,0.010000*%");
    let _ = writeln!(out, "D10*");
    let coord = |v: f64| (v * 1e3 * 1e6).round() as i64; // metres -> mm -> 1e-6 mm
    for l in loops {
        if l.len() < 3 {
            continue;
        }
        let _ = writeln!(out, "G36*");
        let _ = writeln!(out, "X{}Y{}D02*", coord(l[0][0]), coord(l[0][1]));
        for p in &l[1..] {
            let _ = writeln!(out, "X{}Y{}D01*", coord(p[0]), coord(p[1]));
        }
        let _ = writeln!(out, "X{}Y{}D01*", coord(l[0][0]), coord(l[0][1]));
        let _ = writeln!(out, "G37*");
    }
    let _ = writeln!(out, "M02*");
    out
}

/// Excellon drill file for the plated and unplated holes.
///
/// Same caveat as [`gerber_layer`]: written to the format, not verified by
/// a CAM tool.
pub fn excellon(board: &Board) -> String {
    let mut sizes: Vec<f64> = board
        .vias
        .iter()
        .map(|v| v.drill)
        .chain(board.holes.iter().map(|&(_, d)| d))
        .collect();
    sizes.sort_by(|a, b| a.partial_cmp(b).unwrap());
    sizes.dedup_by(|a, b| (*a - *b).abs() < 1e-9);

    let mut out = String::new();
    let _ = writeln!(out, "M48");
    let _ = writeln!(out, "; generated by the cad crate");
    let _ = writeln!(out, "METRIC,TZ");
    for (i, d) in sizes.iter().enumerate() {
        let _ = writeln!(out, "T{:02}C{:.3}", i + 1, d * 1e3);
    }
    let _ = writeln!(out, "%");
    for (i, d) in sizes.iter().enumerate() {
        let _ = writeln!(out, "T{:02}", i + 1);
        for v in &board.vias {
            if (v.drill - d).abs() < 1e-9 {
                let _ = writeln!(out, "X{:.3}Y{:.3}", v.center[0] * 1e3, v.center[1] * 1e3);
            }
        }
        for &(c, hd) in &board.holes {
            if (hd - d).abs() < 1e-9 {
                let _ = writeln!(out, "X{:.3}Y{:.3}", c[0] * 1e3, c[1] * 1e3);
            }
        }
    }
    let _ = writeln!(out, "M30");
    out
}

/// Trace the boundary loops of a 2D region at resolution `cell`, by
/// marching squares on its field.
///
/// The loops are closed polylines in no particular winding order; they are
/// what the Gerber writer emits. Resolution-limited: a feature thinner
/// than one cell can be lost.
pub fn flatten(region: &Region2, cell: f64) -> Vec<Vec<V2>> {
    let (lo, hi) = region.bbox();
    if !lo[0].is_finite() || !hi[0].is_finite() {
        return Vec::new();
    }
    let pad = 2.0 * cell;
    let (lo, hi) = ([lo[0] - pad, lo[1] - pad], [hi[0] + pad, hi[1] + pad]);
    let n: [usize; 2] = std::array::from_fn(|k| (((hi[k] - lo[k]) / cell).ceil() as usize).max(2) + 1);
    let at = |i: usize, j: usize| [lo[0] + cell * i as f64, lo[1] + cell * j as f64];
    let f: Vec<f64> = (0..n[0] * n[1])
        .map(|q| region.sdf(at(q % n[0], q / n[0])))
        .collect();
    let val = |i: usize, j: usize| f[i + n[0] * j];

    // Collect the segments of each cell, then chain them into loops.
    let mut segments: Vec<([f64; 2], [f64; 2])> = Vec::new();
    let interp = |a: V2, b: V2, fa: f64, fb: f64| -> V2 {
        let t = if (fb - fa).abs() > 0.0 { (-fa / (fb - fa)).clamp(0.0, 1.0) } else { 0.5 };
        [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]
    };
    for j in 0..n[1] - 1 {
        for i in 0..n[0] - 1 {
            let (p00, p10, p11, p01) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            let (v00, v10, v11, v01) = (val(i, j), val(i + 1, j), val(i + 1, j + 1), val(i, j + 1));
            let mut crossings: Vec<V2> = Vec::new();
            for (pa, pb, fa, fb) in [
                (p00, p10, v00, v10),
                (p10, p11, v10, v11),
                (p11, p01, v11, v01),
                (p01, p00, v01, v00),
            ] {
                if (fa < 0.0) != (fb < 0.0) {
                    crossings.push(interp(pa, pb, fa, fb));
                }
            }
            if crossings.len() == 2 {
                segments.push((crossings[0], crossings[1]));
            } else if crossings.len() == 4 {
                // Saddle: join the nearer pairs.
                segments.push((crossings[0], crossings[1]));
                segments.push((crossings[2], crossings[3]));
            }
        }
    }

    // Chain segments by endpoint proximity.
    let tol = 0.25 * cell;
    let key = |p: V2| ((p[0] / tol).round() as i64, (p[1] / tol).round() as i64);
    let mut used = vec![false; segments.len()];
    let mut index: std::collections::HashMap<(i64, i64), Vec<usize>> = std::collections::HashMap::new();
    for (s, seg) in segments.iter().enumerate() {
        index.entry(key(seg.0)).or_default().push(s);
        index.entry(key(seg.1)).or_default().push(s);
    }
    let mut loops = Vec::new();
    for start in 0..segments.len() {
        if used[start] {
            continue;
        }
        used[start] = true;
        let mut points = vec![segments[start].0, segments[start].1];
        loop {
            let tail = *points.last().unwrap();
            let mut advanced = false;
            if let Some(cands) = index.get(&key(tail)) {
                for &s in cands {
                    if used[s] {
                        continue;
                    }
                    let (a, b) = segments[s];
                    let next = if key(a) == key(tail) {
                        b
                    } else if key(b) == key(tail) {
                        a
                    } else {
                        continue;
                    };
                    used[s] = true;
                    points.push(next);
                    advanced = true;
                    break;
                }
            }
            if !advanced {
                break;
            }
        }
        if points.len() >= 3 {
            loops.push(points);
        }
    }
    loops
}
