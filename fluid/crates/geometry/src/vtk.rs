//! VTK legacy-format writers (binary), readable by ParaView and VisIt:
//! uniform grids with cell or point data, and triangle surfaces with
//! per-vertex and per-triangle data (e.g. pressure coefficient and wall
//! shear on the CAD surface).

use std::io::Write;
use std::path::Path;

use crate::mesh::{GeomError, TriMesh};
use crate::P3;

/// One named data array.
pub enum Data<'a> {
    Scalars(&'a [f64]),
    Vectors(&'a [P3]),
}

impl Data<'_> {
    fn len(&self) -> usize {
        match self {
            Data::Scalars(s) => s.len(),
            Data::Vectors(v) => v.len(),
        }
    }
}

fn put_f64(out: &mut Vec<u8>, v: f64) {
    out.extend_from_slice(&v.to_be_bytes());
}

fn write_arrays(out: &mut Vec<u8>, kind: &str, n: usize, arrays: &[(&str, Data)]) -> Result<(), GeomError> {
    if arrays.is_empty() {
        return Ok(());
    }
    writeln!(out, "{kind} {n}")?;
    for (name, data) in arrays {
        if data.len() != n {
            return Err(GeomError::Parse(format!("VTK array `{name}` has {} values, expected {n}", data.len())));
        }
        let name = name.replace(' ', "_");
        match data {
            Data::Scalars(s) => {
                write!(out, "SCALARS {name} double 1\nLOOKUP_TABLE default\n")?;
                for &v in s.iter() {
                    put_f64(out, v);
                }
            }
            Data::Vectors(v) => {
                writeln!(out, "VECTORS {name} double")?;
                for p in v.iter() {
                    for &c in p {
                        put_f64(out, c);
                    }
                }
            }
        }
        out.push(b'\n');
    }
    Ok(())
}

/// A uniform grid of `cells[0] × cells[1] × cells[2]` cells (`x` fastest),
/// origin at the low corner, with cell-centred and node data.
pub fn write_grid(
    path: impl AsRef<Path>,
    origin: P3,
    spacing: P3,
    cells: [usize; 3],
    cell_data: &[(&str, Data)],
    point_data: &[(&str, Data)],
) -> Result<(), GeomError> {
    let mut out = Vec::new();
    write!(
        out,
        "# vtk DataFile Version 3.0\nfluid immersed-boundary output\nBINARY\nDATASET STRUCTURED_POINTS\nDIMENSIONS {} {} {}\nORIGIN {} {} {}\nSPACING {} {} {}\n",
        cells[0] + 1,
        cells[1] + 1,
        cells[2] + 1,
        origin[0],
        origin[1],
        origin[2],
        spacing[0],
        spacing[1],
        spacing[2]
    )?;
    write_arrays(&mut out, "CELL_DATA", cells.iter().product(), cell_data)?;
    write_arrays(&mut out, "POINT_DATA", cells.iter().map(|c| c + 1).product(), point_data)?;
    Ok(std::fs::write(path, out)?)
}

/// A triangle surface with per-vertex and per-triangle data.
pub fn write_surface(path: impl AsRef<Path>, mesh: &TriMesh, point_data: &[(&str, Data)], cell_data: &[(&str, Data)]) -> Result<(), GeomError> {
    let mut out = Vec::new();
    write!(out, "# vtk DataFile Version 3.0\nfluid surface output\nBINARY\nDATASET POLYDATA\nPOINTS {} double\n", mesh.vertices.len())?;
    for p in &mesh.vertices {
        for &c in p {
            put_f64(&mut out, c);
        }
    }
    write!(out, "\nPOLYGONS {} {}\n", mesh.triangles.len(), 4 * mesh.triangles.len())?;
    for t in &mesh.triangles {
        out.extend_from_slice(&3i32.to_be_bytes());
        for &v in t {
            out.extend_from_slice(&(v as i32).to_be_bytes());
        }
    }
    out.push(b'\n');
    write_arrays(&mut out, "POINT_DATA", mesh.vertices.len(), point_data)?;
    write_arrays(&mut out, "CELL_DATA", mesh.triangles.len(), cell_data)?;
    Ok(std::fs::write(path, out)?)
}
