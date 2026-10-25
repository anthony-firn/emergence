//! Handing geometry to the solvers.
//!
//! Each solver in the workspace takes geometry in its own form. This module
//! is the set of closures between the kernel and those forms, so a model is
//! authored once and used by every rung:
//!
//! | solver | what it needs | what this module gives it |
//! |---|---|---|
//! | `immersed` (Navier–Stokes) | watertight surface + exact signed distance | [`write_stl`], [`body_sdf`], [`case_body`] |
//! | any grid solver | material label per cell | [`LabelGrid`], [`label_grid`] |
//! | `photonics`, `phononic` | a layer stack | [`LayerStack`] |
//! | `semiconductor` | doping and region at a point | [`kits::silicon`](crate::kits::silicon) |
//! | `architected` | a strut graph in a periodic cell | [`StrutGraph`] |
//! | `heat`, `materials` | per-material properties on the label grid | [`LabelGrid::fractions`] |
//! | ParaView | fields on a grid, data on a surface | [`write_vtk_field`], [`write_vtk_surface`] |

use std::collections::BTreeMap;
use std::path::Path;

use geometry::vtk;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::meshing::{mesh_solid, Grid3};
use crate::solid::Solid;
use geometry::{P3, TriMesh};

/// Write a solid as a binary STL at resolution `cell`.
pub fn write_stl(solid: &Solid, cell: f64, path: impl AsRef<Path>) -> Result<TriMesh, String> {
    let mesh = mesh_solid(solid, cell)?;
    mesh.write_stl(path).map_err(|e| format!("{e:?}"))?;
    Ok(mesh)
}

/// Mesh a solid and wrap it in the exact signed distance the flow solver
/// uses for its immersed boundary.
pub fn body_sdf(solid: &Solid, cell: f64) -> Result<geometry::Sdf, String> {
    Ok(geometry::Sdf::new(mesh_solid(solid, cell)?))
}

/// One body entry for an `immersed` case file: the STL is written next to
/// the case and referenced by name.
pub fn case_body(
    solid: &Solid,
    cell: f64,
    dir: impl AsRef<Path>,
    name: &str,
) -> Result<serde_json::Value, String> {
    let dir = dir.as_ref();
    std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    let file = format!("{name}.stl");
    let mesh = write_stl(solid, cell, dir.join(&file))?;
    let report = mesh.report();
    if !report.watertight {
        return Err(format!(
            "{name}: meshed body is not watertight ({} boundary, {} non-manifold edges)",
            report.boundary_edges, report.nonmanifold_edges
        ));
    }
    Ok(serde_json::json!({ "file": file }))
}

/// A material label per grid cell, with `0` reserved for the background.
#[derive(Clone, Debug)]
pub struct LabelGrid {
    pub grid: Grid3,
    /// Label of each **cell** (`x` fastest), indexed `i + nx*(j + ny*k)`
    /// with `nx = grid.dims[0] - 1`.
    pub labels: Vec<u16>,
    pub names: Vec<String>,
}

impl LabelGrid {
    pub fn cells(&self) -> [usize; 3] {
        self.grid.cells()
    }

    /// Volume fraction of every label.
    pub fn fractions(&self) -> Vec<f64> {
        let mut count = vec![0usize; self.names.len()];
        for &l in &self.labels {
            if (l as usize) < count.len() {
                count[l as usize] += 1;
            }
        }
        let n = self.labels.len() as f64;
        count.iter().map(|&c| c as f64 / n).collect()
    }

    /// A per-cell field built from one value per material, ready for a
    /// finite-volume solver or for VTK.
    pub fn property(&self, values: &[f64]) -> Vec<f64> {
        self.labels
            .iter()
            .map(|&l| values.get(l as usize).copied().unwrap_or(0.0))
            .collect()
    }

    pub fn write_vtk(&self, path: impl AsRef<Path>) -> Result<(), String> {
        let f: Vec<f64> = self.labels.iter().map(|&l| l as f64).collect();
        vtk::write_grid(
            path,
            self.grid.origin,
            self.grid.spacing,
            self.cells(),
            &[("material", vtk::Data::Scalars(&f))],
            &[],
        )
        .map_err(|e| format!("{e:?}"))
    }
}

/// Label every cell centre by the first solid in `materials` that contains
/// it; cells in none get label `0`.
///
/// Later entries do not overwrite earlier ones, so the list reads like a
/// process: substrate first, then what is added on top.
pub fn label_grid(grid: Grid3, background: &str, materials: &[(&str, &Solid)]) -> LabelGrid {
    let cells = grid.cells();
    let total = cells[0] * cells[1] * cells[2];
    let h = grid.spacing;
    let labels: Vec<u16> = (0..total)
        .into_par_iter()
        .map(|q| {
            let (i, j, k) = (q % cells[0], (q / cells[0]) % cells[1], q / (cells[0] * cells[1]));
            let p = [
                grid.origin[0] + h[0] * (i as f64 + 0.5),
                grid.origin[1] + h[1] * (j as f64 + 0.5),
                grid.origin[2] + h[2] * (k as f64 + 0.5),
            ];
            for (n, (_, s)) in materials.iter().enumerate() {
                if s.sdf(p) < 0.0 {
                    return n as u16 + 1;
                }
            }
            0
        })
        .collect();
    let mut names = vec![background.to_string()];
    names.extend(materials.iter().map(|(n, _)| n.to_string()));
    LabelGrid { grid, labels, names }
}

/// One layer of a planar stack: what `photonics::multilayer` and
/// `phononic::layered` consume, and what a PCB or a process flow produces.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Layer {
    pub material: String,
    /// Thickness in metres.
    pub thickness: f64,
    /// Refractive index and extinction coefficient at the wavelength of
    /// interest, when known.
    #[serde(default)]
    pub n: Option<f64>,
    #[serde(default)]
    pub k: Option<f64>,
    /// Relative permittivity at low frequency, when known.
    #[serde(default)]
    pub eps_r: Option<f64>,
    /// Electrical conductivity, S/m.
    #[serde(default)]
    pub sigma: Option<f64>,
    /// Thermal conductivity, W/(m K).
    #[serde(default)]
    pub kappa: Option<f64>,
}

impl Layer {
    pub fn new(material: &str, thickness: f64) -> Self {
        Self {
            material: material.to_string(),
            thickness,
            n: None,
            k: None,
            eps_r: None,
            sigma: None,
            kappa: None,
        }
    }
    pub fn optical(mut self, n: f64, k: f64) -> Self {
        self.n = Some(n);
        self.k = Some(k);
        self
    }
    pub fn electrical(mut self, eps_r: f64, sigma: f64) -> Self {
        self.eps_r = Some(eps_r);
        self.sigma = Some(sigma);
        self
    }
    pub fn thermal(mut self, kappa: f64) -> Self {
        self.kappa = Some(kappa);
        self
    }
}

/// A planar stack, bottom to top.
#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct LayerStack {
    pub layers: Vec<Layer>,
}

impl LayerStack {
    pub fn new(layers: Vec<Layer>) -> Self {
        Self { layers }
    }

    pub fn total_thickness(&self) -> f64 {
        self.layers.iter().map(|l| l.thickness).sum()
    }

    /// `(thickness, n, k)` per layer, the argument order the optics crates
    /// use. Fails when a layer has no optical data.
    pub fn optical(&self) -> Result<Vec<(f64, f64, f64)>, String> {
        self.layers
            .iter()
            .map(|l| match (l.n, l.k) {
                (Some(n), Some(k)) => Ok((l.thickness, n, k)),
                _ => Err(format!("layer {:?} has no optical constants", l.material)),
            })
            .collect()
    }

    /// Series thermal resistance per unit area, m²K/W.
    pub fn thermal_resistance(&self) -> Result<f64, String> {
        self.layers
            .iter()
            .map(|l| l.kappa.map(|k| l.thickness / k).ok_or_else(|| format!("layer {:?} has no kappa", l.material)))
            .sum()
    }

    /// Height of the bottom of each layer above the stack base.
    pub fn offsets(&self) -> Vec<f64> {
        let mut z = 0.0;
        self.layers
            .iter()
            .map(|l| {
                let z0 = z;
                z += l.thickness;
                z0
            })
            .collect()
    }

    /// The stack as a solid: each layer a slab over the footprint `region`.
    pub fn solid(&self, region: &crate::sketch::Region2) -> Solid {
        let offsets = self.offsets();
        Solid::union(
            self.layers
                .iter()
                .zip(offsets)
                .map(|(l, z0)| Solid::extrude(region.clone(), z0, z0 + l.thickness))
                .collect(),
        )
    }
}

/// A periodic strut network: the form `architected::Cell3d` takes.
#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct StrutGraph {
    /// Node positions in the unit cell, in cell-fraction coordinates.
    pub nodes: Vec<P3>,
    /// `(from, to, image)`: a strut from `from` to `to` translated by
    /// `image` whole cells, so a periodic connection is `image != [0,0,0]`.
    pub struts: Vec<(usize, usize, [i32; 3])>,
    /// Strut radius, in the same length units as `cell`.
    pub radius: f64,
    pub cell: P3,
}

impl StrutGraph {
    /// Physical position of node `i` in image `img`.
    pub fn position(&self, i: usize, img: [i32; 3]) -> P3 {
        std::array::from_fn(|k| (self.nodes[i][k] + img[k] as f64) * self.cell[k])
    }

    /// Total strut length inside one cell (each periodic strut counted
    /// once).
    pub fn length(&self) -> f64 {
        self.struts
            .iter()
            .map(|&(a, b, img)| {
                let (p, q) = (self.position(a, [0; 3]), self.position(b, img));
                geometry::norm(geometry::sub(q, p))
            })
            .sum()
    }

    /// Relative density if the struts did not overlap: `sum(pi r^2 L) / V`.
    /// This is an upper bound; the solid's measured volume fraction is
    /// lower by the double-counted material at the nodes.
    pub fn relative_density_upper(&self) -> f64 {
        let v = self.cell[0] * self.cell[1] * self.cell[2];
        std::f64::consts::PI * self.radius * self.radius * self.length() / v
    }

    /// The struts as a solid: a union of capsules, clipped to `counts`
    /// whole cells.
    ///
    /// Cells one deep outside the block are generated too, because their
    /// struts reach inside it; without that halo the material near the
    /// block's faces is missing and the volume fraction comes out low.
    pub fn solid(&self, counts: [usize; 3]) -> Solid {
        let mut parts = Vec::new();
        for i in -1..=counts[0] as i32 {
            for j in -1..=counts[1] as i32 {
                for k in -1..=counts[2] as i32 {
                    for &(a, b, img) in &self.struts {
                        let p = self.position(a, [i, j, k]);
                        let q = self.position(b, [i + img[0], j + img[1], k + img[2]]);
                        parts.push(Solid::capsule(p, q, self.radius));
                    }
                }
            }
        }
        let lo = [0.0; 3];
        let hi: P3 = std::array::from_fn(|k| counts[k] as f64 * self.cell[k]);
        Solid::union(parts).clipped(lo, hi)
    }

    /// `(a, b, image)` struts as physical endpoint pairs, for a frame
    /// solver that wants explicit coordinates.
    pub fn members(&self) -> Vec<(P3, P3)> {
        self.struts
            .iter()
            .map(|&(a, b, img)| (self.position(a, [0; 3]), self.position(b, img)))
            .collect()
    }
}

/// Write the kernel field on a grid for ParaView.
pub fn write_vtk_field(solid: &Solid, grid: &Grid3, path: impl AsRef<Path>) -> Result<(), String> {
    let f = crate::meshing::field_on(solid, grid);
    let cells = grid.cells();
    vtk::write_grid(
        path,
        grid.origin,
        grid.spacing,
        cells,
        &[],
        &[("distance", vtk::Data::Scalars(&f))],
    )
    .map_err(|e| format!("{e:?}"))
}

/// Write a surface with optional per-vertex scalars.
pub fn write_vtk_surface(mesh: &TriMesh, scalars: &BTreeMap<String, Vec<f64>>, path: impl AsRef<Path>) -> Result<(), String> {
    let data: Vec<(&str, vtk::Data)> = scalars
        .iter()
        .map(|(k, v)| (k.as_str(), vtk::Data::Scalars(v)))
        .collect();
    vtk::write_surface(path, mesh, &data, &[]).map_err(|e| format!("{e:?}"))
}
