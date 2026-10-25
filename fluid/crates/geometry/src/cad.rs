//! CAD import (STEP, IGES, BREP) by tessellation through gmsh/OpenCASCADE.
//!
//! A B-rep kernel is out of scope for this crate, so exact CAD is read by
//! gmsh's OpenCASCADE importer and tessellated to a surface mesh with a
//! controlled chordal error. The tool lives in a project-local virtualenv
//! (`fluid/scripts/setup-cad.sh` creates `external/tools/cad`); set
//! `FLUID_GMSH_PYTHON` to use another Python with the `gmsh` module.
//!
//! Chordal error: with `n` elements per `2π` of curvature, an arc of radius
//! `R` is approximated by chords whose sagitta is `R (1 − cos(π/n))`; the
//! default `n = 64` gives `1.2e-3 R`.

use std::path::{Path, PathBuf};
use std::process::Command;

use crate::mesh::{GeomError, TriMesh};

/// Tessellation controls.
#[derive(Clone, Copy, Debug)]
pub struct Tessellation {
    /// Elements per `2π` of curvature (chordal accuracy).
    pub elements_per_2pi: u32,
    /// Largest edge length (model units); `0` = no limit.
    pub max_size: f64,
}

impl Default for Tessellation {
    fn default() -> Self {
        Self {
            elements_per_2pi: 64,
            max_size: 0.0,
        }
    }
}

impl Tessellation {
    /// Relative sagitta of the chords on a circle, `1 − cos(π/n)`.
    pub fn relative_chordal_error(&self) -> f64 {
        1.0 - (std::f64::consts::PI / self.elements_per_2pi as f64).cos()
    }
}

/// The Python interpreter that has `gmsh`, if any.
pub fn gmsh_python() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("FLUID_GMSH_PYTHON") {
        return Some(PathBuf::from(p));
    }
    let local = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../external/tools/cad/bin/python");
    local.exists().then_some(local)
}

const SCRIPT: &str = r#"
import sys, gmsh
src, dst, n, hmax = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4])
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.occ.importShapes(src)
gmsh.model.occ.synchronize()
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", n)
if hmax > 0:
    gmsh.option.setNumber("Mesh.MeshSizeMax", hmax)
gmsh.option.setNumber("Mesh.Binary", 1)
gmsh.model.mesh.generate(2)
gmsh.write(dst)
gmsh.finalize()
"#;

/// Tessellate a CAD file (STEP/IGES/BREP) into a welded triangle mesh.
pub fn tessellate(path: &Path, opt: &Tessellation) -> Result<TriMesh, GeomError> {
    let python = gmsh_python().ok_or_else(|| {
        GeomError::Tool("no gmsh available: run fluid/scripts/setup-cad.sh or set FLUID_GMSH_PYTHON (or export STL from your CAD tool)".into())
    })?;
    let out = std::env::temp_dir().join(format!("fluid-cad-{}-{}.stl", std::process::id(), path.file_stem().and_then(|s| s.to_str()).unwrap_or("part")));
    let status = Command::new(&python)
        .arg("-c")
        .arg(SCRIPT)
        .arg(path)
        .arg(&out)
        .arg(opt.elements_per_2pi.to_string())
        .arg(opt.max_size.to_string())
        .output()
        .map_err(|e| GeomError::Tool(format!("could not run {}: {e}", python.display())))?;
    if !status.status.success() {
        return Err(GeomError::Tool(format!("gmsh failed on {}: {}", path.display(), String::from_utf8_lossy(&status.stderr))));
    }
    let bytes = std::fs::read(&out)?;
    let _ = std::fs::remove_file(&out);
    let mut mesh = TriMesh::from_stl_bytes(&bytes)?;
    mesh.weld_relative(1e-9);
    Ok(mesh)
}
