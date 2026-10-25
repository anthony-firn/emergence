//! Periodic lattices: strut networks and triply periodic surfaces.
//!
//! Two families, both of which `architected` and `heat` can consume:
//!
//! * **Strut lattices** — a [`StrutGraph`] of nodes and periodic struts.
//!   The same graph becomes a solid (a union of capsules) for a flow or
//!   thermal solver, and a member list for the beam homogenisation in
//!   `architected::cell::Cell3d`. The topologies are the standard ones:
//!   simple cubic, body-centred, the octet truss (the face-centred-cubic
//!   nearest-neighbour network; Deshpande, Fleck & Ashby, *J. Mech. Phys.
//!   Solids* 49, 2001) and simple cubic with body diagonals.
//! * **Triply periodic minimal surfaces** — the nodal approximations to
//!   Schoen's gyroid and Schwarz's primitive and diamond surfaces, as
//!   either a solid network (`f <= level`) or a sheet (`|f - level| <= t/2`).
//!
//! Two exact properties are used as tests. First, each nodal function is
//! odd under a symmetry of its own lattice, so at `level = 0` the solid
//! fills exactly half the cell; the measured fraction converges to 0.5.
//! Second, the strut graph's length is exact, so its relative density
//! `pi r^2 L / V` is an exact upper bound on the solid's volume fraction,
//! and the gap between them is the material double-counted at the nodes.
//!
//! No claim is made here about matching published *minimal* surface areas:
//! the trigonometric level sets are not the minimal surfaces, only close to
//! them, and the difference was not measured.

use serde::{Deserialize, Serialize};

use crate::sim::StrutGraph;
use crate::solid::{Solid, Tpms};

/// Which strut topology. Each is a Bravais lattice with a basis, wired up
/// along its nearest-neighbour distance.
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Topology {
    /// Simple cubic: struts along the three cell edges (3 per cell).
    SimpleCubic,
    /// Body-centred cubic: corner-to-centre struts (8 per cell).
    BodyCentred,
    /// The octet truss: the face-centred-cubic nearest-neighbour network,
    /// 24 struts of length `a/sqrt(2)` per cell. Stretch-dominated
    /// (Deshpande, Fleck & Ashby 2001).
    Octet,
    /// Simple cubic plus the body diagonals (3 + 8 per cell).
    CubicDiagonal,
}

impl Topology {
    /// Basis sites in cell-fraction coordinates.
    pub fn basis(&self) -> Vec<[f64; 3]> {
        match self {
            Topology::SimpleCubic => vec![[0.0, 0.0, 0.0]],
            Topology::BodyCentred | Topology::CubicDiagonal => vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
            Topology::Octet => vec![
                [0.0, 0.0, 0.0],
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [0.5, 0.5, 0.0],
            ],
        }
    }

    /// Bonds to wire up: a length in units of the cell edge, and
    /// optionally the one pair of basis sites it may join. The restriction
    /// is what separates "simple cubic plus body diagonals", where only
    /// the corner sublattice carries edges, from body-centred cubic with
    /// edges on both sublattices.
    pub fn bonds(&self) -> Vec<(f64, Option<(usize, usize)>)> {
        let d3 = 3.0_f64.sqrt() / 2.0;
        match self {
            Topology::SimpleCubic => vec![(1.0, None)],
            Topology::BodyCentred => vec![(d3, None)],
            Topology::CubicDiagonal => vec![(d3, None), (1.0, Some((0, 0)))],
            Topology::Octet => vec![(1.0 / 2.0_f64.sqrt(), None)],
        }
    }

    /// Expected number of struts per cell, for a test to pin down.
    pub fn struts_per_cell(&self) -> usize {
        match self {
            Topology::SimpleCubic => 3,
            Topology::BodyCentred => 8,
            Topology::CubicDiagonal => 11,
            Topology::Octet => 24,
        }
    }
}

/// Wire up every pair of basis sites whose separation is one of
/// `bonds` (in cell-edge units), counting each strut once.
///
/// A strut is the pair `(i, j, image)`; the same strut seen from the other
/// end is `(j, i, -image)`, so one of the two is kept by a canonical
/// ordering. Cubic cells only: the search covers images in `-1..=1`, which
/// is enough for any bond shorter than one cell edge diagonal.
pub fn nearest_struts(basis: &[[f64; 3]], bonds: &[(f64, Option<(usize, usize)>)]) -> Vec<(usize, usize, [i32; 3])> {
    let mut out = Vec::new();
    let tol = 1e-9;
    for i in 0..basis.len() {
        for j in 0..basis.len() {
            for ix in -1i32..=1 {
                for iy in -1i32..=1 {
                    for iz in -1i32..=1 {
                        let img = [ix, iy, iz];
                        if i == j && img == [0, 0, 0] {
                            continue;
                        }
                        let d: f64 = (0..3)
                            .map(|k| (basis[j][k] + img[k] as f64 - basis[i][k]).powi(2))
                            .sum::<f64>()
                            .sqrt();
                        let matched = bonds.iter().any(|&(l, pair)| {
                            (d - l).abs() < tol && pair.is_none_or(|(x, y)| (x, y) == (i, j) || (x, y) == (j, i))
                        });
                        if !matched {
                            continue;
                        }
                        // Keep one of the two directions.
                        let keep = match i.cmp(&j) {
                            std::cmp::Ordering::Less => true,
                            std::cmp::Ordering::Greater => false,
                            std::cmp::Ordering::Equal => img > [0, 0, 0],
                        };
                        if keep {
                            out.push((i, j, img));
                        }
                    }
                }
            }
        }
    }
    out
}

/// A strut lattice.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct LatticeSpec {
    pub topology: Topology,
    /// Cell size (cubic unless all three are given).
    pub cell: [f64; 3],
    pub radius: f64,
    #[serde(default = "one_cell")]
    pub counts: [usize; 3],
}

fn one_cell() -> [usize; 3] {
    [1, 1, 1]
}

impl LatticeSpec {
    pub fn new(topology: Topology, cell: f64, radius: f64) -> Self {
        Self { topology, cell: [cell; 3], radius, counts: [1, 1, 1] }
    }

    pub fn with_counts(mut self, counts: [usize; 3]) -> Self {
        self.counts = counts;
        self
    }

    /// Nodes and periodic struts in cell-fraction coordinates.
    pub fn graph(&self) -> StrutGraph {
        let nodes = self.topology.basis();
        let struts = nearest_struts(&nodes, &self.topology.bonds());
        StrutGraph { nodes, struts, radius: self.radius, cell: self.cell }
    }

    pub fn solid(&self) -> Solid {
        self.graph().solid(self.counts)
    }
}

/// A triply periodic level set.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct TpmsSpec {
    /// `gyroid`, `schwarz_p` or `diamond`.
    pub kind: String,
    pub cell: [f64; 3],
    /// Level offset; 0 gives the balanced surface.
    #[serde(default)]
    pub level: f64,
    /// Sheet thickness; omit for the solid network `f <= level`.
    #[serde(default)]
    pub thickness: Option<f64>,
    /// Block to clip to, in cells. Omit to leave the solid unbounded.
    #[serde(default)]
    pub counts: Option<[usize; 3]>,
}

impl TpmsSpec {
    pub fn new(kind: &str, cell: f64) -> Self {
        Self { kind: kind.to_string(), cell: [cell; 3], level: 0.0, thickness: None, counts: None }
    }

    pub fn sheet(mut self, thickness: f64) -> Self {
        self.thickness = Some(thickness);
        self
    }

    pub fn block(mut self, counts: [usize; 3]) -> Self {
        self.counts = Some(counts);
        self
    }

    pub fn tpms(&self) -> Result<Tpms, String> {
        Ok(match self.kind.as_str() {
            "gyroid" => Tpms::Gyroid,
            "schwarz_p" | "primitive" => Tpms::SchwarzP,
            "diamond" => Tpms::Diamond,
            other => return Err(format!("unknown periodic surface {other:?}")),
        })
    }

    /// The sheet's achieved geometric wall thickness, measured as
    /// `2 V / A` from the meshed sheet (exact for a thin sheet of uniform
    /// thickness). It exceeds the nominal `thickness` because the nodal
    /// field is normalised by a conservative gradient bound.
    pub fn wall_thickness(&self, cell: f64) -> Result<f64, String> {
        let block = self.counts.unwrap_or([1, 1, 1]);
        let spec = Self { counts: Some(block), ..self.clone() };
        if spec.thickness.is_none() {
            return Err("wall_thickness needs a sheet (set thickness)".into());
        }
        let mesh = crate::meshing::mesh_solid(&spec.solid(), cell)?;
        let report = mesh.report();
        Ok(2.0 * report.volume / report.area)
    }

    pub fn solid(&self) -> Solid {
        let kind = self.tpms().unwrap_or(Tpms::Gyroid);
        let body = Solid::Periodic {
            kind,
            cell: self.cell,
            level: self.level,
            thickness: self.thickness,
        };
        match self.counts {
            None => body,
            Some(c) => {
                let hi: [f64; 3] = std::array::from_fn(|k| c[k] as f64 * self.cell[k]);
                body.clipped([0.0; 3], hi)
            }
        }
    }
}

/// Octet-truss relative density, counting each strut once and ignoring the
/// material double-counted where struts meet.
///
/// The cubic cell of side `a` holds 24 struts of length `l = a/sqrt(2)`, so
/// `rho = 24 pi r^2 l / a^3 = 6 sqrt(2) pi (r/l)^2`, which is the published
/// relation (Deshpande, Fleck & Ashby, *J. Mech. Phys. Solids* 49(8), 2001,
/// eq. 1). It is an upper bound on the solid's volume fraction.
pub fn octet_relative_density(radius: f64, cell: f64) -> f64 {
    let l = cell / 2.0_f64.sqrt();
    6.0 * 2.0_f64.sqrt() * std::f64::consts::PI * (radius / l).powi(2)
}
