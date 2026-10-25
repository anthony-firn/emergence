//! JSON models: a parametric feature tree on disk.
//!
//! A model is a parameter table and a feature tree. Any number in the tree
//! may be written as a string expression over the parameters, so one file
//! describes a family of parts:
//!
//! ```json
//! {
//!   "name": "flange",
//!   "params": { "bolts": 6, "bolt_circle": 0.08, "hole": "0.0033",
//!               "thickness": 0.006, "od": "bolt_circle + 3*hole" },
//!   "shape": {
//!     "op": "difference",
//!     "a": { "op": "cylinder_z", "center": [0,0,0],
//!            "radius": "od/2", "height": "thickness" },
//!     "b": { "op": "union", "parts": [
//!              { "op": "cylinder_z", "center": [0,0,0],
//!                "radius": "od/4", "height": "2*thickness" },
//!              { "op": "polar_array", "axis": [0,0,1], "count": "bolts",
//!                "shape": { "op": "cylinder_z",
//!                           "center": ["bolt_circle/2",0,0],
//!                           "radius": "hole", "height": "2*thickness" } } ] }
//!   }
//! }
//! ```
//!
//! Parameters may refer to one another in any order; they are resolved by
//! repeated passes and a cycle is reported as an error. Expressions are
//! substituted into the tree *before* it is deserialised: every JSON string
//! that evaluates is replaced by its value, except under the keys in
//! [`NAME_KEYS`], which are always names. A string that does not evaluate
//! is left alone, so a parameter typo surfaces as an unknown-field or type
//! error naming the offending key.
//!
//! Units are the caller's; the generators in [`crate::kits`] document
//! theirs (SI metres, except the PCB kit, which also accepts millimetres on
//! input).

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use serde_json::Value as Json;

use crate::expr::Params;
use crate::sketch::Region2;
use crate::solid::{Placement, Solid, Tpms};

/// Keys whose string values are names, never expressions.
///
/// Every serde tag and string-valued discriminant in the crate belongs
/// here. A parameter whose name collides with a variant under some *other*
/// key would be substituted into it — but the result is a loud
/// deserialisation error naming the key ("expected variant identifier"),
/// not silent corruption, so the failure is diagnosable. Extend the list
/// when a new tagged enum is added.
pub const NAME_KEYS: [&str; 18] = [
    "op",
    "kind",
    "file",
    "name",
    "material",
    "kit",
    "leg",
    "part",
    "shape",
    "profile",
    "pattern",
    "structure",
    "topology",
    "species",
    "schedule",
    "state",
    "dopant",
    "relaxation",
];

/// A parametric model as read from disk.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Model {
    #[serde(default)]
    pub name: String,
    /// Parameter definitions: numbers or expressions over other
    /// parameters, in any order.
    #[serde(default)]
    pub params: BTreeMap<String, Json>,
    /// The feature tree, with expressions still unresolved.
    pub shape: Json,
    /// Default meshing resolution, if the file gives one.
    #[serde(default)]
    pub cell: Option<f64>,
}

impl Model {
    pub fn from_json(text: &str) -> Result<Self, String> {
        serde_json::from_str(text).map_err(|e| format!("model JSON: {e}"))
    }

    pub fn read(path: impl AsRef<Path>) -> Result<Self, String> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        Self::from_json(&text)
    }

    /// Resolve the parameter table, with `overrides` applied on top (this is
    /// how a design table or a CLI flag sweeps a dimension).
    pub fn resolve(&self, overrides: &Params) -> Result<Params, String> {
        let mut out = overrides.clone();
        let mut pending: Vec<(String, String)> = Vec::new();
        for (k, v) in &self.params {
            if out.values.contains_key(k) {
                continue; // overridden
            }
            match v {
                Json::Number(n) => {
                    out.set(k, n.as_f64().ok_or_else(|| format!("param {k}: not a number"))?);
                }
                Json::String(s) => pending.push((k.clone(), s.clone())),
                other => return Err(format!("param {k}: expected number or expression, found {other}")),
            }
        }
        for _ in 0..pending.len() + 1 {
            if pending.is_empty() {
                break;
            }
            let before = pending.len();
            pending.retain(|(k, s)| match out.eval(s) {
                Ok(v) => {
                    out.set(k, v);
                    false
                }
                Err(_) => true,
            });
            if pending.len() == before {
                let names: Vec<&str> = pending.iter().map(|(k, _)| k.as_str()).collect();
                return Err(format!(
                    "parameters {names:?} cannot be resolved (unknown name or a cycle)"
                ));
            }
        }
        Ok(out)
    }

    /// Build the solid. `base` is the directory that `file` references in
    /// mesh features are resolved against.
    pub fn build(&self, overrides: &Params, base: impl AsRef<Path>) -> Result<Solid, String> {
        let params = self.resolve(overrides)?;
        let mut tree = self.shape.clone();
        substitute(&mut tree, &params, false)?;
        let feature: Feature = serde_json::from_value(tree).map_err(|e| format!("shape: {e}"))?;
        feature.build(base.as_ref())
    }
}

/// Replace every expression string in a JSON tree by its value.
fn substitute(v: &mut Json, params: &Params, is_name: bool) -> Result<(), String> {
    match v {
        Json::String(s) => {
            if is_name {
                return Ok(());
            }
            if let Ok(x) = params.eval(s) {
                // Emit an integer when the value is one, so counts and
                // tooth numbers deserialize into integer fields.
                *v = if x.fract() == 0.0 && x.abs() < 9.0e15 {
                    serde_json::json!(x as i64)
                } else {
                    serde_json::json!(x)
                };
            }
            Ok(())
        }
        Json::Array(a) => {
            for x in a.iter_mut() {
                substitute(x, params, false)?;
            }
            Ok(())
        }
        Json::Object(o) => {
            for (k, x) in o.iter_mut() {
                substitute(x, params, NAME_KEYS.contains(&k.as_str()))?;
            }
            Ok(())
        }
        _ => Ok(()),
    }
}

fn default_one() -> f64 {
    1.0
}

/// One node of a JSON feature tree.
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(tag = "op", rename_all = "snake_case", deny_unknown_fields)]
pub enum Feature {
    Empty,
    Sphere { center: P3, radius: f64 },
    Box { center: P3, half: P3, #[serde(default)] round: f64 },
    BoxBetween { min: P3, max: P3 },
    Cylinder { a: P3, b: P3, radius: f64 },
    CylinderZ { center: P3, radius: f64, height: f64 },
    Cone { a: P3, b: P3, radius_a: f64, radius_b: f64 },
    Capsule { a: P3, b: P3, radius: f64 },
    Torus { center: P3, major: f64, minor: f64 },
    HalfSpace { normal: P3, offset: f64 },
    Extrude { sketch: Region2, z0: f64, z1: f64 },
    Revolve { sketch: Region2, #[serde(default)] angle: Option<f64> },
    /// An STL/OBJ/STEP file through `geometry`.
    Mesh {
        file: String,
        #[serde(default = "default_one")]
        scale: f64,
        #[serde(default)]
        translate: P3,
        /// Treat the exterior as solid (a housing or pipe wall).
        #[serde(default)]
        inside_out: bool,
    },
    Union { parts: Vec<Feature> },
    Intersection { parts: Vec<Feature> },
    Difference { a: Box<Feature>, b: Box<Feature> },
    Blend { a: Box<Feature>, b: Box<Feature>, k: f64 },
    Offset { shape: Box<Feature>, distance: f64 },
    /// Morphological rounding; `cell` is the remesh resolution it needs
    /// (see [`Solid::round_convex`]).
    RoundConvex { shape: Box<Feature>, radius: f64, cell: f64 },
    RoundConcave { shape: Box<Feature>, radius: f64, cell: f64 },
    Shell { shape: Box<Feature>, thickness: f64 },
    Translate { shape: Box<Feature>, by: P3 },
    Rotate { shape: Box<Feature>, axis: P3, angle: f64 },
    Scale { shape: Box<Feature>, factor: f64 },
    Stretch { shape: Box<Feature>, factors: P3 },
    Mirror { shape: Box<Feature>, normal: P3 },
    Twist { shape: Box<Feature>, rate: f64, radius: f64 },
    Array { shape: Box<Feature>, period: P3, counts: [usize; 3] },
    PolarArray { shape: Box<Feature>, axis: P3, count: usize },
    Periodic {
        kind: String,
        cell: P3,
        #[serde(default)]
        level: f64,
        #[serde(default)]
        thickness: Option<f64>,
    },
    Clip { shape: Box<Feature>, min: P3, max: P3 },
    /// A generator from [`crate::kits`], by name, with its own arguments.
    Kit { kit: String, args: Json },
}

type P3 = [f64; 3];

impl Feature {
    pub fn build(&self, base: &Path) -> Result<Solid, String> {
        Ok(match self {
            Feature::Empty => Solid::Empty,
            Feature::Sphere { center, radius } => Solid::sphere(*center, *radius),
            Feature::Box { center, half, round } => Solid::rounded_box(*center, *half, *round),
            Feature::BoxBetween { min, max } => Solid::box_between(*min, *max),
            Feature::Cylinder { a, b, radius } => Solid::cylinder(*a, *b, *radius),
            Feature::CylinderZ { center, radius, height } => Solid::cylinder_z(*center, *radius, *height),
            Feature::Cone { a, b, radius_a, radius_b } => Solid::cone(*a, *b, *radius_a, *radius_b),
            Feature::Capsule { a, b, radius } => Solid::capsule(*a, *b, *radius),
            Feature::Torus { center, major, minor } => Solid::torus(*center, *major, *minor),
            Feature::HalfSpace { normal, offset } => Solid::half_space(*normal, *offset),
            Feature::Extrude { sketch, z0, z1 } => Solid::extrude(sketch.clone(), *z0, *z1),
            Feature::Revolve { sketch, angle } => match angle {
                None => Solid::revolve(sketch.clone()),
                Some(a) => Solid::revolve_arc(sketch.clone(), *a),
            },
            Feature::Mesh { file, scale, translate, inside_out } => {
                let path: PathBuf = if Path::new(file).is_absolute() {
                    file.into()
                } else {
                    base.join(file)
                };
                let mut mesh = read_mesh(&path)?;
                if *scale != 1.0 {
                    mesh.scale_by(*scale);
                }
                if *translate != [0.0; 3] {
                    mesh.translate(*translate);
                }
                if *inside_out {
                    mesh.orient_outward();
                    mesh.flip();
                    Solid::mesh(geometry::Sdf::new_as_oriented(mesh))
                } else {
                    Solid::from_trimesh(mesh)
                }
            }
            Feature::Union { parts } => Solid::union(build_all(parts, base)?),
            Feature::Intersection { parts } => Solid::intersection(build_all(parts, base)?),
            Feature::Difference { a, b } => a.build(base)?.minus(b.build(base)?),
            Feature::Blend { a, b, k } => a.build(base)?.blend(b.build(base)?, *k),
            Feature::Offset { shape, distance } => shape.build(base)?.offset(*distance),
            Feature::RoundConvex { shape, radius, cell } => shape.build(base)?.round_convex(*radius, *cell)?,
            Feature::RoundConcave { shape, radius, cell } => shape.build(base)?.round_concave(*radius, *cell)?,
            Feature::Shell { shape, thickness } => shape.build(base)?.shell(*thickness),
            Feature::Translate { shape, by } => shape.build(base)?.translated(*by),
            Feature::Rotate { shape, axis, angle } => shape.build(base)?.rotated(*axis, *angle),
            Feature::Scale { shape, factor } => shape.build(base)?.scaled(*factor),
            Feature::Stretch { shape, factors } => shape.build(base)?.stretched(*factors),
            Feature::Mirror { shape, normal } => shape.build(base)?.mirrored(*normal),
            Feature::Twist { shape, rate, radius } => shape.build(base)?.twisted(*rate, *radius),
            Feature::Array { shape, period, counts } => shape.build(base)?.array(*period, *counts),
            Feature::PolarArray { shape, axis, count } => shape.build(base)?.polar_array(*axis, *count),
            Feature::Periodic { kind, cell, level, thickness } => Solid::Periodic {
                kind: match kind.as_str() {
                    "gyroid" => Tpms::Gyroid,
                    "schwarz_p" | "primitive" => Tpms::SchwarzP,
                    "diamond" => Tpms::Diamond,
                    other => return Err(format!("unknown periodic surface {other:?}")),
                },
                cell: *cell,
                level: *level,
                thickness: *thickness,
            },
            Feature::Clip { shape, min, max } => shape.build(base)?.clipped(*min, *max),
            Feature::Kit { kit, args } => crate::kits::from_json(kit, args)?,
        })
    }

    /// A placement helper for callers assembling features in Rust.
    pub fn placed(self, base: &Path, placement: Placement) -> Result<Solid, String> {
        Ok(self.build(base)?.placed(placement))
    }
}

fn build_all(parts: &[Feature], base: &Path) -> Result<Vec<Solid>, String> {
    parts.iter().map(|f| f.build(base)).collect()
}

fn read_mesh(path: &Path) -> Result<geometry::TriMesh, String> {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
    match ext.as_str() {
        "stl" | "obj" => geometry::TriMesh::read(path).map_err(|e| format!("{}: {e:?}", path.display())),
        "step" | "stp" | "iges" | "igs" | "brep" => {
            geometry::cad::tessellate(path, &geometry::cad::Tessellation::default())
                .map_err(|e| format!("{}: {e:?}", path.display()))
        }
        other => Err(format!("{}: unknown geometry extension {other:?}", path.display())),
    }
}
