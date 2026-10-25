//! Triangle surface meshes: STL/OBJ input and output, welding, topology
//! checks and mass properties.

use std::collections::HashMap;
use std::fmt;
use std::path::Path;

use crate::{add, cross, dot, norm, scale, sub, P3};

/// Errors from reading or converting geometry.
#[derive(Debug)]
pub enum GeomError {
    Io(std::io::Error),
    Parse(String),
    Unsupported(String),
    Tool(String),
}

impl fmt::Display for GeomError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            GeomError::Io(e) => write!(f, "i/o: {e}"),
            GeomError::Parse(s) => write!(f, "parse: {s}"),
            GeomError::Unsupported(s) => write!(f, "unsupported: {s}"),
            GeomError::Tool(s) => write!(f, "external tool: {s}"),
        }
    }
}

impl std::error::Error for GeomError {}

impl From<std::io::Error> for GeomError {
    fn from(e: std::io::Error) -> Self {
        GeomError::Io(e)
    }
}

/// An indexed triangle mesh. Triangles are counter-clockwise seen from
/// outside (outward normals) once [`TriMesh::orient_outward`] has run.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TriMesh {
    pub vertices: Vec<P3>,
    pub triangles: Vec<[u32; 3]>,
}

/// Topology and quality report.
#[derive(Clone, Debug, PartialEq)]
pub struct MeshReport {
    pub vertices: usize,
    pub triangles: usize,
    /// Edges used by one triangle only (holes).
    pub boundary_edges: usize,
    /// Edges used by more than two triangles.
    pub nonmanifold_edges: usize,
    /// Edges whose two triangles traverse them in the same direction.
    pub inconsistent_edges: usize,
    /// Triangles with (relative) zero area.
    pub degenerate: usize,
    /// Closed, manifold and consistently oriented.
    pub watertight: bool,
    /// Signed enclosed volume (positive for outward normals).
    pub volume: f64,
    pub area: f64,
    pub bbox: (P3, P3),
}

impl TriMesh {
    pub fn new(vertices: Vec<P3>, triangles: Vec<[u32; 3]>) -> Self {
        Self { vertices, triangles }
    }

    /// Read a surface by file extension: `.stl`, `.obj`, or a CAD format
    /// (`.step`, `.stp`, `.iges`, `.igs`, `.brep`) through [`crate::cad`].
    /// STL soups are welded so the result has connectivity.
    pub fn read(path: impl AsRef<Path>) -> Result<Self, GeomError> {
        let path = path.as_ref();
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
        let mut mesh = match ext.as_str() {
            "stl" => Self::from_stl_bytes(&std::fs::read(path)?)?,
            "obj" => Self::from_obj(&std::fs::read_to_string(path)?)?,
            "step" | "stp" | "iges" | "igs" | "brep" => crate::cad::tessellate(path, &crate::cad::Tessellation::default())?,
            other => return Err(GeomError::Unsupported(format!("surface format `.{other}`"))),
        };
        mesh.weld_relative(1e-9);
        Ok(mesh)
    }

    /// Parse STL, binary or ASCII. The result is a triangle soup (three
    /// vertices per triangle); call [`TriMesh::weld`] for connectivity.
    pub fn from_stl_bytes(bytes: &[u8]) -> Result<Self, GeomError> {
        // Binary: 80-byte header, u32 count, 50 bytes per triangle. An ASCII
        // file can start with "solid" too, so trust the size check first.
        if bytes.len() >= 84 {
            let n = u32::from_le_bytes(bytes[80..84].try_into().unwrap()) as usize;
            if bytes.len() == 84 + 50 * n {
                return Ok(Self::from_binary_stl(&bytes[84..], n));
            }
        }
        let text = std::str::from_utf8(bytes).map_err(|_| GeomError::Parse("STL is neither binary (size mismatch) nor UTF-8 ASCII".into()))?;
        Self::from_ascii_stl(text)
    }

    fn from_binary_stl(data: &[u8], n: usize) -> Self {
        let mut vertices = Vec::with_capacity(3 * n);
        let mut triangles = Vec::with_capacity(n);
        let f = |o: usize| f32::from_le_bytes(data[o..o + 4].try_into().unwrap()) as f64;
        for t in 0..n {
            let base = 50 * t + 12; // skip the stored normal
            for v in 0..3 {
                let o = base + 12 * v;
                vertices.push([f(o), f(o + 4), f(o + 8)]);
            }
            let i = 3 * t as u32;
            triangles.push([i, i + 1, i + 2]);
        }
        Self { vertices, triangles }
    }

    fn from_ascii_stl(text: &str) -> Result<Self, GeomError> {
        let mut vertices = Vec::new();
        for (line_no, line) in text.lines().enumerate() {
            let mut it = line.split_whitespace();
            if it.next() == Some("vertex") {
                let mut v = [0.0; 3];
                for c in v.iter_mut() {
                    *c = it
                        .next()
                        .and_then(|s| s.parse().ok())
                        .ok_or_else(|| GeomError::Parse(format!("bad vertex on line {}", line_no + 1)))?;
                }
                vertices.push(v);
            }
        }
        if vertices.is_empty() || vertices.len() % 3 != 0 {
            return Err(GeomError::Parse(format!("ASCII STL has {} vertices (not a positive multiple of 3)", vertices.len())));
        }
        let triangles = (0..vertices.len() as u32 / 3).map(|t| [3 * t, 3 * t + 1, 3 * t + 2]).collect();
        Ok(Self { vertices, triangles })
    }

    /// Parse Wavefront OBJ (`v` and `f` records; polygons are fanned into
    /// triangles; `v/vt/vn` and negative indices are accepted).
    pub fn from_obj(text: &str) -> Result<Self, GeomError> {
        let mut vertices = Vec::new();
        let mut triangles = Vec::new();
        for (line_no, line) in text.lines().enumerate() {
            let mut it = line.split_whitespace();
            match it.next() {
                Some("v") => {
                    let mut v = [0.0; 3];
                    for c in v.iter_mut() {
                        *c = it
                            .next()
                            .and_then(|s| s.parse().ok())
                            .ok_or_else(|| GeomError::Parse(format!("bad vertex on line {}", line_no + 1)))?;
                    }
                    vertices.push(v);
                }
                Some("f") => {
                    let idx: Result<Vec<u32>, GeomError> = it
                        .map(|tok| {
                            let k: i64 = tok
                                .split('/')
                                .next()
                                .and_then(|s| s.parse().ok())
                                .ok_or_else(|| GeomError::Parse(format!("bad face on line {}", line_no + 1)))?;
                            let k = if k < 0 { vertices.len() as i64 + k } else { k - 1 };
                            if k < 0 || k >= vertices.len() as i64 {
                                return Err(GeomError::Parse(format!("face index out of range on line {}", line_no + 1)));
                            }
                            Ok(k as u32)
                        })
                        .collect();
                    let idx = idx?;
                    for w in 1..idx.len().saturating_sub(1) {
                        triangles.push([idx[0], idx[w], idx[w + 1]]);
                    }
                }
                _ => {}
            }
        }
        Ok(Self { vertices, triangles })
    }

    /// Binary STL bytes (normals computed from the vertices).
    pub fn to_stl_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(84 + 50 * self.triangles.len());
        let mut header = [0u8; 80];
        let tag = b"fluid geometry crate binary STL";
        header[..tag.len()].copy_from_slice(tag);
        out.extend_from_slice(&header);
        out.extend_from_slice(&(self.triangles.len() as u32).to_le_bytes());
        for t in 0..self.triangles.len() {
            let n = self.unit_normal(t);
            for c in n {
                out.extend_from_slice(&(c as f32).to_le_bytes());
            }
            for &v in &self.triangles[t] {
                for c in self.vertices[v as usize] {
                    out.extend_from_slice(&(c as f32).to_le_bytes());
                }
            }
            out.extend_from_slice(&[0, 0]);
        }
        out
    }

    pub fn write_stl(&self, path: impl AsRef<Path>) -> Result<(), GeomError> {
        Ok(std::fs::write(path, self.to_stl_bytes())?)
    }

    /// Merge vertices closer than `tol` (absolute), drop triangles that
    /// collapse, and drop unused vertices.
    pub fn weld(&mut self, tol: f64) {
        let inv = 1.0 / tol.max(f64::MIN_POSITIVE);
        let key = |p: P3| [(p[0] * inv).round() as i64, (p[1] * inv).round() as i64, (p[2] * inv).round() as i64];
        let mut cells: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
        let mut remap = vec![0u32; self.vertices.len()];
        let mut kept: Vec<P3> = Vec::new();
        for (i, &p) in self.vertices.iter().enumerate() {
            let k = key(p);
            let mut found = None;
            // A neighbour within tol can sit in an adjacent cell.
            'search: for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(list) = cells.get(&[k[0] + dx, k[1] + dy, k[2] + dz]) {
                            for &j in list {
                                if norm(sub(kept[j as usize], p)) <= tol {
                                    found = Some(j);
                                    break 'search;
                                }
                            }
                        }
                    }
                }
            }
            remap[i] = match found {
                Some(j) => j,
                None => {
                    kept.push(p);
                    let j = kept.len() as u32 - 1;
                    cells.entry(k).or_default().push(j);
                    j
                }
            };
        }
        self.vertices = kept;
        self.triangles = self
            .triangles
            .iter()
            .map(|t| t.map(|v| remap[v as usize]))
            .filter(|t| t[0] != t[1] && t[1] != t[2] && t[0] != t[2])
            .collect();
    }

    /// [`TriMesh::weld`] with `tol` relative to the bounding-box diagonal.
    pub fn weld_relative(&mut self, rel: f64) {
        let (lo, hi) = self.bbox();
        let diag = norm(sub(hi, lo));
        self.weld(rel * diag.max(f64::MIN_POSITIVE));
    }

    pub fn bbox(&self) -> (P3, P3) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for v in &self.vertices {
            for d in 0..3 {
                lo[d] = lo[d].min(v[d]);
                hi[d] = hi[d].max(v[d]);
            }
        }
        (lo, hi)
    }

    pub fn corners(&self, t: usize) -> [P3; 3] {
        self.triangles[t].map(|v| self.vertices[v as usize])
    }

    /// Area-weighted normal (length = 2 × area).
    pub fn normal(&self, t: usize) -> P3 {
        let [a, b, c] = self.corners(t);
        cross(sub(b, a), sub(c, a))
    }

    pub fn unit_normal(&self, t: usize) -> P3 {
        let n = self.normal(t);
        let l = norm(n);
        if l > 0.0 {
            scale(n, 1.0 / l)
        } else {
            [0.0; 3]
        }
    }

    pub fn triangle_area(&self, t: usize) -> f64 {
        0.5 * norm(self.normal(t))
    }

    pub fn area(&self) -> f64 {
        (0..self.triangles.len()).map(|t| self.triangle_area(t)).sum()
    }

    /// Signed volume by the divergence theorem (exact for the polyhedron).
    pub fn volume(&self) -> f64 {
        (0..self.triangles.len())
            .map(|t| {
                let [a, b, c] = self.corners(t);
                dot(a, cross(b, c)) / 6.0
            })
            .sum()
    }

    /// Centroid of the enclosed solid.
    pub fn centroid(&self) -> P3 {
        let mut s = [0.0; 3];
        let mut v = 0.0;
        for t in 0..self.triangles.len() {
            let [a, b, c] = self.corners(t);
            let dv = dot(a, cross(b, c)) / 6.0;
            v += dv;
            s = add(s, scale(add(add(a, b), c), dv / 4.0));
        }
        scale(s, 1.0 / v)
    }

    /// Reverse every triangle.
    pub fn flip(&mut self) {
        for t in &mut self.triangles {
            t.swap(1, 2);
        }
    }

    /// Make normals point outward (positive enclosed volume). Assumes a
    /// consistently oriented closed surface.
    pub fn orient_outward(&mut self) {
        if self.volume() < 0.0 {
            self.flip();
        }
    }

    pub fn transform(&mut self, f: impl Fn(P3) -> P3) {
        for v in &mut self.vertices {
            *v = f(*v);
        }
    }

    pub fn translate(&mut self, d: P3) {
        self.transform(|p| add(p, d));
    }

    pub fn scale_by(&mut self, s: f64) {
        self.transform(|p| scale(p, s));
    }

    /// Append another mesh.
    pub fn merge(&mut self, other: &TriMesh) {
        let off = self.vertices.len() as u32;
        self.vertices.extend_from_slice(&other.vertices);
        self.triangles.extend(other.triangles.iter().map(|t| t.map(|v| v + off)));
    }

    /// Undirected edges with the triangles using them, keyed `(min, max)`;
    /// the bool is whether the triangle traverses the edge min → max.
    pub fn edge_map(&self) -> HashMap<(u32, u32), Vec<(u32, bool)>> {
        let mut m: HashMap<(u32, u32), Vec<(u32, bool)>> = HashMap::with_capacity(3 * self.triangles.len() / 2);
        for (t, tri) in self.triangles.iter().enumerate() {
            for k in 0..3 {
                let (a, b) = (tri[k], tri[(k + 1) % 3]);
                m.entry((a.min(b), a.max(b))).or_default().push((t as u32, a < b));
            }
        }
        m
    }

    pub fn report(&self) -> MeshReport {
        let edges = self.edge_map();
        let (mut boundary, mut nonmanifold, mut inconsistent) = (0, 0, 0);
        for uses in edges.values() {
            match uses.len() {
                1 => boundary += 1,
                2 => {
                    if uses[0].1 == uses[1].1 {
                        inconsistent += 1;
                    }
                }
                _ => nonmanifold += 1,
            }
        }
        let (lo, hi) = self.bbox();
        let diag = norm(sub(hi, lo));
        let degenerate = (0..self.triangles.len())
            .filter(|&t| self.triangle_area(t) <= 1e-14 * diag * diag)
            .count();
        MeshReport {
            vertices: self.vertices.len(),
            triangles: self.triangles.len(),
            boundary_edges: boundary,
            nonmanifold_edges: nonmanifold,
            inconsistent_edges: inconsistent,
            degenerate,
            watertight: boundary == 0 && nonmanifold == 0 && inconsistent == 0 && !self.triangles.is_empty(),
            volume: self.volume(),
            area: self.area(),
            bbox: (lo, hi),
        }
    }
}
