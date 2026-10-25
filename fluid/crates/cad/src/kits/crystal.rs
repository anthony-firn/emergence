//! Atomistic structures: crystals, nanoparticles, nanowires, graphene and
//! nanotubes.
//!
//! "Anything that can be simulated" includes things that are not solids at
//! all. A molecular-dynamics run or a tight-binding calculation takes a
//! list of atom positions and a periodic cell, and building that list is
//! as much a geometry problem as drawing a flange: pick a lattice, cut a
//! shape out of it, roll a sheet into a tube. This kit does that, and
//! hands the result to `materials` and `md` in the form each wants.
//!
//! ```text
//!   Crystal { lattice, basis, a } -> Cell  --supercell--> atoms + box
//!                                         --cut by a Solid--> nanoparticle,
//!                                             nanowire, any shape at all
//!   Sheet (graphene) --roll (n, m)--> nanotube
//! ```
//!
//! Lengths are metres, as everywhere else in the crate; lattice constants
//! are therefore in the 1e-10 range. Atom species are carried as strings.
//!
//! # What is checked against published values
//!
//! * **Coordination numbers**: 12 for close packing, 8 for body-centred
//!   cubic, 6 for simple cubic, 4 for diamond — counted from the generated
//!   positions, not asserted.
//! * **Packing fractions**: `pi/(3 sqrt 2) = 0.740480` for face-centred
//!   cubic and ideal hexagonal close packing, `pi sqrt 3 / 8 = 0.680175`
//!   for body-centred cubic, `pi/6` for simple cubic and
//!   `pi sqrt 3 / 16 = 0.340087` for diamond (Ashcroft & Mermin, *Solid
//!   State Physics*, 1976, ch. 4).
//! * **Nearest-neighbour distances**: `a/sqrt 2`, `a sqrt 3 / 2`, `a` and
//!   `a sqrt 3 / 4` for the same four.
//! * **Nanotube geometry**: an `(n, m)` tube rolled from graphene of
//!   lattice constant `a` has circumference `a sqrt(n^2 + nm + m^2)`, so
//!   diameter `a sqrt(n^2 + nm + m^2) / pi`, and chiral angle
//!   `atan(sqrt 3 m / (2n + m))` (Saito, Dresselhaus & Dresselhaus,
//!   *Physical Properties of Carbon Nanotubes*, 1998, ch. 3). The test
//!   measures all three off the generated atoms.
//!
//! # Limits
//!
//! Ideal structures only: no defects, no vacancies, no substitutional
//! disorder, no surface relaxation or reconstruction, and no equilibration
//! — a cut nanoparticle has the bulk lattice right up to its surface,
//! which is not what a relaxed one looks like. Rolling a sheet into a tube
//! places the atoms on a cylinder exactly, which is the standard
//! construction but ignores the curvature-induced bond-length changes a
//! relaxation would find.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use geometry::{add, norm, scale, sub, P3};

use crate::solid::Solid;

/// A Bravais lattice with a basis.
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Structure {
    /// Simple cubic.
    Sc,
    /// Body-centred cubic.
    Bcc,
    /// Face-centred cubic.
    Fcc,
    /// Hexagonal close packed with the ideal `c/a = sqrt(8/3)`.
    Hcp,
    /// Diamond cubic.
    Diamond,
    /// Zincblende (diamond with two species).
    Zincblende,
    /// Rock salt.
    RockSalt,
    /// Cubic perovskite `ABX3`.
    Perovskite,
}

impl Structure {
    /// Primitive (here: conventional) cell vectors for lattice constant
    /// `a`.
    pub fn cell_vectors(&self, a: f64) -> [P3; 3] {
        match self {
            Structure::Hcp => {
                let c = a * (8.0f64 / 3.0).sqrt();
                [[a, 0.0, 0.0], [-0.5 * a, a * 3.0f64.sqrt() / 2.0, 0.0], [0.0, 0.0, c]]
            }
            _ => [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]],
        }
    }

    /// Basis sites in cell-fraction coordinates, with a species index.
    pub fn basis(&self) -> Vec<(P3, usize)> {
        match self {
            Structure::Sc => vec![([0.0, 0.0, 0.0], 0)],
            Structure::Bcc => vec![([0.0, 0.0, 0.0], 0), ([0.5, 0.5, 0.5], 0)],
            Structure::Fcc => vec![
                ([0.0, 0.0, 0.0], 0),
                ([0.0, 0.5, 0.5], 0),
                ([0.5, 0.0, 0.5], 0),
                ([0.5, 0.5, 0.0], 0),
            ],
            // (2/3, 1/3, 1/2) in this cell's fractional coordinates puts
            // the second atom exactly `a` from its three neighbours in the
            // layer below, which is what makes the packing close.
            Structure::Hcp => vec![
                ([0.0, 0.0, 0.0], 0),
                ([2.0 / 3.0, 1.0 / 3.0, 0.5], 0),
            ],
            Structure::Diamond => {
                let mut v: Vec<(P3, usize)> = Structure::Fcc.basis();
                for (p, _) in Structure::Fcc.basis() {
                    v.push(([p[0] + 0.25, p[1] + 0.25, p[2] + 0.25], 0));
                }
                v
            }
            Structure::Zincblende => {
                let mut v: Vec<(P3, usize)> = Structure::Fcc.basis();
                for (p, _) in Structure::Fcc.basis() {
                    v.push(([p[0] + 0.25, p[1] + 0.25, p[2] + 0.25], 1));
                }
                v
            }
            Structure::RockSalt => {
                let mut v: Vec<(P3, usize)> = Structure::Fcc.basis();
                for (p, _) in Structure::Fcc.basis() {
                    v.push(([p[0] + 0.5, p[1], p[2]], 1));
                }
                v
            }
            Structure::Perovskite => vec![
                ([0.0, 0.0, 0.0], 0),              // A at the corner
                ([0.5, 0.5, 0.5], 1),              // B at the centre
                ([0.5, 0.5, 0.0], 2),              // X on the faces
                ([0.5, 0.0, 0.5], 2),
                ([0.0, 0.5, 0.5], 2),
            ],
        }
    }

    pub fn atoms_per_cell(&self) -> usize {
        self.basis().len()
    }

    /// Nearest-neighbour distance in units of `a`.
    pub fn nn_over_a(&self) -> f64 {
        match self {
            Structure::Sc => 1.0,
            Structure::Bcc => 3.0f64.sqrt() / 2.0,
            Structure::Fcc => 1.0 / 2.0f64.sqrt(),
            // For the hexagonal cell, `a` *is* the nearest-neighbour
            // distance; the cubic a/sqrt(2) of face-centred cubic does not
            // apply.
            Structure::Hcp => 1.0,
            Structure::Diamond | Structure::Zincblende => 3.0f64.sqrt() / 4.0,
            Structure::RockSalt => 0.5,
            Structure::Perovskite => 0.5,
        }
    }

    /// Published coordination number of the first shell.
    pub fn coordination(&self) -> usize {
        match self {
            Structure::Sc => 6,
            Structure::Bcc => 8,
            Structure::Fcc | Structure::Hcp => 12,
            Structure::Diamond | Structure::Zincblende => 4,
            Structure::RockSalt => 6,
            Structure::Perovskite => 6,
        }
    }

    /// Hard-sphere packing fraction of the monatomic structures, taking
    /// the sphere radius as half the nearest-neighbour distance.
    pub fn packing_fraction(&self) -> Option<f64> {
        let n = match self {
            Structure::Sc => 1.0,
            Structure::Bcc => 2.0,
            Structure::Fcc => 4.0,
            Structure::Diamond => 8.0,
            Structure::Hcp => 2.0,
            _ => return None,
        };
        let r = 0.5 * self.nn_over_a();
        let cell_volume = match self {
            // The hexagonal cell holds two atoms in a volume
            // a^2 sqrt(3)/2 * c with c = a sqrt(8/3).
            Structure::Hcp => 3.0f64.sqrt() / 2.0 * (8.0f64 / 3.0).sqrt(),
            _ => 1.0,
        };
        Some(n * 4.0 / 3.0 * PI * r.powi(3) / cell_volume)
    }
}

/// One atom.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Atom {
    pub species: String,
    pub position: P3,
}

/// A generated structure: atoms plus the periodic cell they sit in.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Atoms {
    pub atoms: Vec<Atom>,
    /// Periodic cell vectors; `None` for a finite cluster.
    pub cell: Option<[P3; 3]>,
}

impl Atoms {
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Count of each species.
    pub fn composition(&self) -> BTreeMap<String, usize> {
        let mut out = BTreeMap::new();
        for a in &self.atoms {
            *out.entry(a.species.clone()).or_insert(0) += 1;
        }
        out
    }

    /// Cell volume; `None` for a cluster.
    pub fn volume(&self) -> Option<f64> {
        let c = self.cell?;
        Some(geometry::dot(c[0], geometry::cross(c[1], c[2])).abs())
    }

    /// Number density, atoms per cubic metre.
    pub fn number_density(&self) -> Option<f64> {
        Some(self.len() as f64 / self.volume()?)
    }

    /// Mass density for a map of species to atomic mass in kg.
    pub fn mass_density(&self, masses: &BTreeMap<String, f64>) -> Option<f64> {
        let mass: f64 = self.atoms.iter().map(|a| masses.get(&a.species).copied().unwrap_or(0.0)).sum();
        Some(mass / self.volume()?)
    }

    /// Smallest separation between any two atoms, by brute force.
    pub fn min_separation(&self) -> f64 {
        let mut best = f64::INFINITY;
        for i in 0..self.atoms.len() {
            for j in i + 1..self.atoms.len() {
                let d = norm(sub(self.atoms[j].position, self.atoms[i].position));
                if d < best {
                    best = d;
                }
            }
        }
        best
    }

    /// How many neighbours each atom has within `cutoff`, with the
    /// periodic cell applied when there is one.
    pub fn coordination(&self, cutoff: f64) -> Vec<usize> {
        let images: Vec<P3> = match self.cell {
            None => vec![[0.0; 3]],
            Some(c) => {
                let mut v = Vec::with_capacity(27);
                for i in -1..=1i32 {
                    for j in -1..=1i32 {
                        for k in -1..=1i32 {
                            v.push(add(
                                add(scale(c[0], i as f64), scale(c[1], j as f64)),
                                scale(c[2], k as f64),
                            ));
                        }
                    }
                }
                v
            }
        };
        self.atoms
            .iter()
            .map(|a| {
                let mut count = 0;
                for b in &self.atoms {
                    for img in &images {
                        let d = norm(sub(add(b.position, *img), a.position));
                        if d > 1e-15 && d <= cutoff {
                            count += 1;
                        }
                    }
                }
                count
            })
            .collect()
    }

    /// Bounding box of the atom positions.
    pub fn bbox(&self) -> (P3, P3) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for a in &self.atoms {
            for k in 0..3 {
                lo[k] = lo[k].min(a.position[k]);
                hi[k] = hi[k].max(a.position[k]);
            }
        }
        (lo, hi)
    }

    /// Centre the structure on its own centroid.
    pub fn centred(&self) -> Atoms {
        let n = self.len() as f64;
        if n == 0.0 {
            return self.clone();
        }
        let mut c = [0.0; 3];
        for a in &self.atoms {
            c = add(c, a.position);
        }
        c = scale(c, 1.0 / n);
        Atoms {
            atoms: self
                .atoms
                .iter()
                .map(|a| Atom { species: a.species.clone(), position: sub(a.position, c) })
                .collect(),
            cell: self.cell,
        }
    }

    /// Keep only the atoms inside `solid`: this is how a nanoparticle, a
    /// nanowire or any other shape is cut out of a crystal.
    pub fn cut_to(&self, solid: &Solid) -> Atoms {
        Atoms {
            atoms: self
                .atoms
                .iter()
                .filter(|a| solid.sdf(a.position) < 0.0)
                .cloned()
                .collect(),
            cell: None,
        }
    }

    /// Positions and species as the two parallel arrays most solvers take.
    pub fn arrays(&self) -> (Vec<P3>, Vec<String>) {
        (
            self.atoms.iter().map(|a| a.position).collect(),
            self.atoms.iter().map(|a| a.species.clone()).collect(),
        )
    }

    /// Extended XYZ, which every atomistic tool reads.
    pub fn to_xyz(&self) -> String {
        use std::fmt::Write as _;
        let mut out = String::new();
        let _ = writeln!(out, "{}", self.len());
        match self.cell {
            Some(c) => {
                // Angstroms, as the format expects.
                let _ = writeln!(
                    out,
                    "Lattice=\"{:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6}\" Properties=species:S:1:pos:R:3",
                    c[0][0] * 1e10, c[0][1] * 1e10, c[0][2] * 1e10,
                    c[1][0] * 1e10, c[1][1] * 1e10, c[1][2] * 1e10,
                    c[2][0] * 1e10, c[2][1] * 1e10, c[2][2] * 1e10
                );
            }
            None => {
                let _ = writeln!(out, "Properties=species:S:1:pos:R:3");
            }
        }
        for a in &self.atoms {
            let _ = writeln!(
                out,
                "{} {:.6} {:.6} {:.6}",
                a.species,
                a.position[0] * 1e10,
                a.position[1] * 1e10,
                a.position[2] * 1e10
            );
        }
        out
    }
}

/// A crystal to build from.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Crystal {
    pub structure: Structure,
    /// Cubic (or hexagonal `a`) lattice constant, metres.
    pub a: f64,
    /// One species name per basis species index.
    pub species: Vec<String>,
}

impl Crystal {
    pub fn new(structure: Structure, a: f64, species: &[&str]) -> Self {
        Self {
            structure,
            a,
            species: species.iter().map(|s| s.to_string()).collect(),
        }
    }

    /// Copper: face-centred cubic, a = 3.615 A (CRC Handbook, 97th ed.).
    pub fn copper() -> Self {
        Self::new(Structure::Fcc, 3.615e-10, &["Cu"])
    }

    /// Iron: body-centred cubic, a = 2.866 A (CRC Handbook).
    pub fn iron() -> Self {
        Self::new(Structure::Bcc, 2.866e-10, &["Fe"])
    }

    /// Silicon: diamond cubic, a = 5.431 A (CRC Handbook).
    pub fn silicon() -> Self {
        Self::new(Structure::Diamond, 5.431e-10, &["Si"])
    }

    /// Nearest-neighbour distance.
    pub fn nn_distance(&self) -> f64 {
        self.a * self.structure.nn_over_a()
    }

    fn species_of(&self, index: usize) -> String {
        self.species
            .get(index)
            .cloned()
            .unwrap_or_else(|| format!("X{index}"))
    }

    /// A supercell of `counts` conventional cells, with the periodic cell
    /// set so the structure tiles.
    pub fn supercell(&self, counts: [usize; 3]) -> Atoms {
        let v = self.structure.cell_vectors(self.a);
        let basis = self.structure.basis();
        let mut atoms = Vec::with_capacity(counts.iter().product::<usize>() * basis.len());
        for i in 0..counts[0] {
            for j in 0..counts[1] {
                for k in 0..counts[2] {
                    let origin = add(
                        add(scale(v[0], i as f64), scale(v[1], j as f64)),
                        scale(v[2], k as f64),
                    );
                    for (frac, species) in &basis {
                        let p = add(
                            origin,
                            add(
                                add(scale(v[0], frac[0]), scale(v[1], frac[1])),
                                scale(v[2], frac[2]),
                            ),
                        );
                        atoms.push(Atom { species: self.species_of(*species), position: p });
                    }
                }
            }
        }
        Atoms {
            atoms,
            cell: Some([
                scale(v[0], counts[0] as f64),
                scale(v[1], counts[1] as f64),
                scale(v[2], counts[2] as f64),
            ]),
        }
    }

    /// Cut a shape out of the crystal. The supercell is sized from the
    /// solid's bounding box with one cell of margin, so any solid the
    /// kernel can describe becomes a cluster.
    pub fn carve(&self, solid: &Solid) -> Result<Atoms, String> {
        let (lo, hi) = solid.bbox();
        if (0..3).any(|k| !lo[k].is_finite() || !hi[k].is_finite()) {
            return Err("cannot carve from an unbounded solid".into());
        }
        let v = self.structure.cell_vectors(self.a);
        let extent: P3 = std::array::from_fn(|k| hi[k] - lo[k]);
        let span: P3 = std::array::from_fn(|k| v[0][k].abs() + v[1][k].abs() + v[2][k].abs());
        let counts: [usize; 3] =
            std::array::from_fn(|k| (extent[k] / span[k].max(1e-30)).ceil() as usize + 3);
        let big = self.supercell(counts);
        // Shift so the supercell covers the solid's box.
        let (blo, _) = big.bbox();
        let shift: P3 = std::array::from_fn(|k| lo[k] - blo[k] - span[k]);
        let shifted = Atoms {
            atoms: big
                .atoms
                .into_iter()
                .map(|a| Atom { species: a.species, position: add(a.position, shift) })
                .collect(),
            cell: None,
        };
        Ok(shifted.cut_to(solid))
    }

    /// A spherical nanoparticle of the given diameter.
    pub fn nanoparticle(&self, diameter: f64) -> Result<Atoms, String> {
        Ok(self
            .carve(&Solid::sphere([0.0; 3], 0.5 * diameter))?
            .centred())
    }

    /// A nanowire: a cylinder along `+z` of the given diameter and length.
    pub fn nanowire(&self, diameter: f64, length: f64) -> Result<Atoms, String> {
        self.carve(&Solid::cylinder(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, length],
            0.5 * diameter,
        ))
    }
}

/// A honeycomb sheet: graphene, or any two-site hexagonal monolayer.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Sheet {
    /// Hexagonal lattice constant (2.46 A for graphene; Saito et al.
    /// 1998, ch. 3).
    pub a: f64,
    pub species: String,
}

impl Sheet {
    pub fn graphene() -> Self {
        Self { a: 2.46e-10, species: "C".into() }
    }

    /// Carbon-carbon bond length, `a / sqrt 3` (1.42 A for graphene).
    pub fn bond_length(&self) -> f64 {
        self.a / 3.0f64.sqrt()
    }

    /// Lattice vectors of the hexagonal cell.
    pub fn vectors(&self) -> [[f64; 2]; 2] {
        [
            [self.a * 3.0f64.sqrt() / 2.0, self.a / 2.0],
            [self.a * 3.0f64.sqrt() / 2.0, -self.a / 2.0],
        ]
    }

    /// A flat sheet of `counts` cells, two atoms each, in the `xy` plane.
    pub fn flat(&self, counts: [usize; 2]) -> Atoms {
        let v = self.vectors();
        let b = self.bond_length();
        let mut atoms = Vec::with_capacity(2 * counts[0] * counts[1]);
        for i in 0..counts[0] {
            for j in 0..counts[1] {
                let o = [
                    v[0][0] * i as f64 + v[1][0] * j as f64,
                    v[0][1] * i as f64 + v[1][1] * j as f64,
                ];
                atoms.push(Atom { species: self.species.clone(), position: [o[0], o[1], 0.0] });
                atoms.push(Atom {
                    species: self.species.clone(),
                    position: [o[0] + b, o[1], 0.0],
                });
            }
        }
        Atoms {
            atoms,
            cell: Some([
                [v[0][0] * counts[0] as f64, v[0][1] * counts[0] as f64, 0.0],
                [v[1][0] * counts[1] as f64, v[1][1] * counts[1] as f64, 0.0],
                [0.0, 0.0, 1e-9],
            ]),
        }
    }

    /// Diameter of the `(n, m)` nanotube rolled from this sheet:
    /// `a sqrt(n^2 + n m + m^2) / pi`.
    pub fn tube_diameter(&self, n: usize, m: usize) -> f64 {
        self.a * ((n * n + n * m + m * m) as f64).sqrt() / PI
    }

    /// Chiral angle of the `(n, m)` tube, `atan(sqrt 3 m / (2n + m))`.
    pub fn chiral_angle(&self, n: usize, m: usize) -> f64 {
        (3.0f64.sqrt() * m as f64).atan2((2 * n + m) as f64)
    }

    /// Atoms in one translational unit cell of the `(n, m)` tube:
    /// `4 (n^2 + nm + m^2) / d_R` with `d_R = gcd(2n + m, 2m + n)`.
    pub fn tube_atoms_per_cell(&self, n: usize, m: usize) -> usize {
        let dr = gcd(2 * n + m, 2 * m + n);
        4 * (n * n + n * m + m * m) / dr
    }

    /// Length of that translational unit cell,
    /// `sqrt(3) a sqrt(n^2 + nm + m^2) / d_R`.
    pub fn tube_cell_length(&self, n: usize, m: usize) -> f64 {
        let dr = gcd(2 * n + m, 2 * m + n) as f64;
        3.0f64.sqrt() * self.a * ((n * n + n * m + m * m) as f64).sqrt() / dr
    }

    /// Roll the sheet into an `(n, m)` nanotube of `cells` translational
    /// unit cells along its axis.
    ///
    /// The chiral vector `C = n a1 + m a2` becomes the circumference and
    /// the translation vector `T` the axis. Atoms are placed on the
    /// cylinder exactly: the construction is the standard one and does not
    /// relax the curvature-induced bond changes.
    pub fn nanotube(&self, n: usize, m: usize, cells: usize) -> Result<Atoms, String> {
        if n == 0 && m == 0 {
            return Err("a nanotube needs a non-zero chiral index".into());
        }
        let v = self.vectors();
        let chiral = [
            v[0][0] * n as f64 + v[1][0] * m as f64,
            v[0][1] * n as f64 + v[1][1] * m as f64,
        ];
        let circumference = (chiral[0] * chiral[0] + chiral[1] * chiral[1]).sqrt();
        let radius = circumference / (2.0 * PI);
        // The translation vector: the shortest lattice vector normal to C.
        let dr = gcd(2 * n + m, 2 * m + n) as i64;
        let (t1, t2) = ((2 * m as i64 + n as i64) / dr, -(2 * n as i64 + m as i64) / dr);
        let trans = [
            v[0][0] * t1 as f64 + v[1][0] * t2 as f64,
            v[0][1] * t1 as f64 + v[1][1] * t2 as f64,
        ];
        let t_len = (trans[0] * trans[0] + trans[1] * trans[1]).sqrt();
        // Unit vectors along the chiral and translation directions.
        let ch = [chiral[0] / circumference, chiral[1] / circumference];
        let th = [trans[0] / t_len, trans[1] / t_len];

        // Enough sheet cells to cover the (C, T) parallelogram, then keep
        // the atoms whose (u, v) coordinates fall inside it.
        let reach = ((n + m + t1.unsigned_abs() as usize + t2.unsigned_abs() as usize) + 2) as i64;
        let b = self.bond_length();
        let mut atoms = Vec::new();
        let mut seen: Vec<[f64; 2]> = Vec::new();
        for i in -reach..=reach {
            for j in -reach..=reach {
                let o = [
                    v[0][0] * i as f64 + v[1][0] * j as f64,
                    v[0][1] * i as f64 + v[1][1] * j as f64,
                ];
                for site in [[o[0], o[1]], [o[0] + b, o[1]]] {
                    let u = site[0] * ch[0] + site[1] * ch[1];
                    let w = site[0] * th[0] + site[1] * th[1];
                    // Fold into the cell, with a tolerance that keeps one
                    // copy of each boundary atom.
                    let tol = 1e-9 * self.a;
                    if u < -tol || u >= circumference - tol || w < -tol || w >= t_len - tol {
                        continue;
                    }
                    if seen
                        .iter()
                        .any(|s| (s[0] - u).abs() < 1e-6 * self.a && (s[1] - w).abs() < 1e-6 * self.a)
                    {
                        continue;
                    }
                    seen.push([u, w]);
                    for c in 0..cells.max(1) {
                        let angle = 2.0 * PI * u / circumference;
                        atoms.push(Atom {
                            species: self.species.clone(),
                            position: [
                                radius * angle.cos(),
                                radius * angle.sin(),
                                w + t_len * c as f64,
                            ],
                        });
                    }
                }
            }
        }
        Ok(Atoms {
            atoms,
            cell: Some([
                [4.0 * radius, 0.0, 0.0],
                [0.0, 4.0 * radius, 0.0],
                [0.0, 0.0, t_len * cells.max(1) as f64],
            ]),
        })
    }
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 {
        a.max(1)
    } else {
        gcd(b, a % b)
    }
}
