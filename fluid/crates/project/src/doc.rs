//! The project file.
//!
//! Plain JSON with sorted keys, so it diffs and merges in git like code.
//! Node ids are stable: a label can change freely, and every reference
//! (expressions, analysis inputs, part bindings, requirements) uses the id.
//!
//! ```json
//! {
//!   "format": "fluid-project/1",
//!   "name": "Motor driver",
//!   "params":   { "trace.width": { "value": "3 mm", "tol": "0.1 mm", "owner": "electrical" } },
//!   "derived":  { "q1.tj": { "expr": "ambient + q1.power * q1.rja", "owner": "thermal", "unit": "degC" } },
//!   "analyses": { "sink": { "kind": "plate_fin_heatsink", "owner": "thermal",
//!                           "inputs": { "width": "sink.width" }, "materials": { "fin": "al6063" } } },
//!   "parts":    { "sink_body": { "model": "heatsink.json", "params": { "width": "sink.width" },
//!                              "material": "al6063", "owner": "mechanical" } },
//!   "requirements": { "REQ-TJ": { "text": "...", "of": "q1.tj", "max": "125 degC", "owner": "systems" } }
//! }
//! ```

use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

pub const FORMAT: &str = "fluid-project/1";

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Doc {
    #[serde(default)]
    pub format: String,
    pub name: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub description: String,
    /// Team name → what it owns (free text).
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub teams: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub materials: BTreeMap<String, MaterialDef>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub params: BTreeMap<String, ParamDef>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub derived: BTreeMap<String, DerivedDef>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub analyses: BTreeMap<String, AnalysisDef>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub parts: BTreeMap<String, PartDef>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub requirements: BTreeMap<String, RequirementDef>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaterialDef {
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub description: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub owner: String,
    pub props: BTreeMap<String, PropDef>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PropDef {
    pub value: String,
    /// Relative uncertainty, e.g. `"20%"`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tol: Option<String>,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub source: String,
}

/// A tolerance: `"10%"` (relative), `"0.1 mm"` (absolute, ±), or an explicit
/// `["35 degC", "45 degC"]` range.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum Tol {
    Sym(String),
    Range([String; 2]),
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ParamDef {
    pub value: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tol: Option<Tol>,
    pub owner: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub label: String,
    /// Display unit (default: the unit `value` is written in).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unit: Option<String>,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub note: String,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DerivedDef {
    pub expr: String,
    pub owner: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub label: String,
    /// Display unit; also checked against the expression's dimension.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unit: Option<String>,
    /// Makes this node the unknown of a feedback loop it sits on, solved as
    /// a fixed point inside `[lo, hi]` (e.g. `["25 degC", "175 degC"]`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub range: Option<[String; 2]>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AnalysisDef {
    pub kind: String,
    pub owner: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub label: String,
    /// Input name → expression.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub inputs: BTreeMap<String, String>,
    /// Material slot → material name.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub materials: BTreeMap<String, String>,
    /// Non-numeric settings, e.g. `{"rtl": "cpu.v"}` for a chip design.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub config: BTreeMap<String, String>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PartDef {
    /// Parametric CAD model (`cad` JSON), relative to the project file.
    pub model: String,
    /// CAD parameter → expression (lengths in any unit; converted to SI).
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub params: BTreeMap<String, String>,
    pub material: String,
    pub owner: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub label: String,
    /// Mesh cell, e.g. `"0.25 mm"` (default: the model's own).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cell: Option<String>,
    /// Printed by fused deposition: bindings for `nozzle`, `perimeters`,
    /// `infill` and `max_flow`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fdm: Option<BTreeMap<String, String>>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RequirementDef {
    pub text: String,
    /// Expression the limits apply to.
    pub of: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min: Option<String>,
    pub owner: String,
}

impl Doc {
    pub fn from_json(text: &str) -> Result<Doc, String> {
        let doc: Doc = serde_json::from_str(text).map_err(|e| format!("project JSON: {e}"))?;
        if !doc.format.is_empty() && doc.format != FORMAT {
            return Err(format!("unsupported format `{}` (this build reads {FORMAT})", doc.format));
        }
        Ok(doc)
    }

    pub fn read(path: impl AsRef<Path>) -> Result<Doc, String> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        Doc::from_json(&text).map_err(|e| format!("{}: {e}", path.display()))
    }

    pub fn to_json(&self) -> String {
        let mut s = serde_json::to_string_pretty(self).expect("serialisable");
        s.push('\n');
        s
    }

    /// The owner of a node id (params, derived, analyses, parts and their
    /// outputs, requirements).
    pub fn owner(&self, id: &str) -> Option<&str> {
        if let Some(p) = self.params.get(id) {
            return Some(&p.owner);
        }
        if let Some(d) = self.derived.get(id) {
            return Some(&d.owner);
        }
        if let Some(r) = self.requirements.get(id) {
            return Some(&r.owner);
        }
        let head = id.rsplit_once('.').map(|(h, _)| h).unwrap_or(id);
        for key in [id, head] {
            if let Some(a) = self.analyses.get(key) {
                return Some(&a.owner);
            }
            if let Some(p) = self.parts.get(key) {
                return Some(&p.owner);
            }
        }
        if let Some((m, _)) = crate::materials::Library::split(id) {
            return Some(self.materials.get(m).map(|d| d.owner.as_str()).filter(|o| !o.is_empty()).unwrap_or("materials library"));
        }
        None
    }
}
