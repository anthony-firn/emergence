//! Real-world geometry for the solvers.
//!
//! * [`mesh`]: triangle surfaces: STL (binary/ASCII) and OBJ input, binary
//!   STL output, welding, watertightness/orientation checks, exact volume,
//!   area and centroid;
//! * [`cad`]: STEP/IGES/BREP through gmsh's OpenCASCADE importer, with a
//!   controlled chordal error;
//! * [`bvh`] and [`sdf`]: exact closest points and signed distance
//!   (pseudonormal sign for watertight surfaces, winding number otherwise);
//! * [`primitives`]: closed test shapes;
//! * [`vtk`]: ParaView-readable output for grids and surfaces.

pub mod bvh;
pub mod cad;
pub mod mesh;
pub mod primitives;
pub mod sdf;
pub mod vtk;

pub use mesh::{GeomError, MeshReport, TriMesh};
pub use sdf::{Closest, Sdf, SignMode};

/// A point or vector in 3D.
pub type P3 = [f64; 3];

#[inline]
pub fn add(a: P3, b: P3) -> P3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub fn sub(a: P3, b: P3) -> P3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub fn scale(a: P3, s: f64) -> P3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub fn dot(a: P3, b: P3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub fn cross(a: P3, b: P3) -> P3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
#[inline]
pub fn norm(a: P3) -> f64 {
    dot(a, a).sqrt()
}
#[inline]
pub fn normalize(a: P3) -> P3 {
    let n = norm(a);
    if n > 0.0 {
        scale(a, 1.0 / n)
    } else {
        a
    }
}
