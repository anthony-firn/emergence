//! Case files: a whole engineering run described in JSON.
//!
//! ```json
//! {
//!   "name": "cylinder-in-channel",
//!   "fluid": { "density": 1.0, "kinematic_viscosity": 0.001 },
//!   "domain": { "min": [0, 0, 0], "max": [2.2, 0.41, 0.005], "cell": 0.005 },
//!   "boundaries": {
//!     "x_min": { "type": "inlet", "velocity": [0.3, 0, 0], "profile": "parabolic" },
//!     "x_max": { "type": "outlet" },
//!     "y_min": { "type": "wall" }, "y_max": { "type": "wall" },
//!     "z_min": { "type": "periodic" }, "z_max": { "type": "periodic" }
//!   },
//!   "bodies": [ { "file": "cylinder.step" } ],
//!   "run": { "steady": true, "t_end": 60, "cfl": 0.8, "tolerance": 1e-6, "initial": "inlet" },
//!   "output": { "directory": "out", "reference": { "velocity": 0.2, "area": 0.0005 },
//!               "probes": [[0.15, 0.2, 0.0025]], "vtk": true }
//! }
//! ```
//!
//! Boundary types: `inlet` (uniform or `parabolic` profile, `velocity` is
//! the peak), `wall` (no-slip; `velocity` makes it a moving wall), `outlet`
//! (pressure outlet, `p = 0`), `slip` (symmetry), `periodic`. A parabolic
//! profile is the product of `4s(1 − s)` over the face's non-periodic
//! tangential directions (`s` the normalised coordinate), as in the DFG
//! channel benchmarks. Bodies are read by extension (`.stl`, `.obj`, `.step`,
//! `.stp`, `.iges`, `.igs`, `.brep`), relative to the case file.

use std::path::{Path, PathBuf};
use std::time::Instant;

use geometry::cad::Tessellation;
use geometry::{Sdf, TriMesh, P3};
use serde::{Deserialize, Serialize};

use crate::field::Grid;
use crate::solver::{Bc, Body, Solver};

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Case {
    #[serde(default = "default_name")]
    pub name: String,
    pub fluid: Fluid,
    pub domain: Domain,
    pub boundaries: Boundaries,
    #[serde(default)]
    pub bodies: Vec<BodySpec>,
    #[serde(default)]
    pub run: Run,
    #[serde(default)]
    pub output: Output,
}

fn default_name() -> String {
    "case".into()
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Fluid {
    pub density: f64,
    pub kinematic_viscosity: f64,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Domain {
    pub min: P3,
    pub max: P3,
    /// Cubic cell size; every extent must be a whole number of cells.
    pub cell: f64,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Boundaries {
    pub x_min: Face,
    pub x_max: Face,
    pub y_min: Face,
    pub y_max: Face,
    pub z_min: Face,
    pub z_max: Face,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Face {
    Inlet {
        velocity: P3,
        #[serde(default)]
        profile: Profile,
    },
    Wall {
        #[serde(default)]
        velocity: P3,
    },
    Outlet,
    Slip,
    Periodic,
}

#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Profile {
    #[default]
    Uniform,
    Parabolic,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct BodySpec {
    pub file: PathBuf,
    #[serde(default = "one")]
    pub scale: f64,
    #[serde(default)]
    pub translate: P3,
    #[serde(default)]
    pub velocity: P3,
    #[serde(default)]
    pub angular_velocity: P3,
    #[serde(default)]
    pub center: P3,
    /// The surface encloses the fluid (housing, pipe): keep its
    /// orientation and treat the outside as solid. The file's normals must
    /// then point into the fluid.
    #[serde(default)]
    pub inside_out: bool,
    /// CAD tessellation: elements per 2π of curvature.
    #[serde(default = "default_curvature")]
    pub elements_per_2pi: u32,
}

fn one() -> f64 {
    1.0
}
fn default_curvature() -> u32 {
    64
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Run {
    /// Stop when `max |Δu|/dt < tolerance` (steady) or at `t_end`.
    #[serde(default = "yes")]
    pub steady: bool,
    #[serde(default = "default_t_end")]
    pub t_end: f64,
    #[serde(default = "default_cfl")]
    pub cfl: f64,
    #[serde(default = "default_tol")]
    pub tolerance: f64,
    /// `"inlet"`: start from the inlet profile extended through the domain;
    /// `"rest"`: start from rest.
    #[serde(default = "default_initial")]
    pub initial: String,
}

fn yes() -> bool {
    true
}
fn default_t_end() -> f64 {
    50.0
}
fn default_cfl() -> f64 {
    0.8
}
fn default_tol() -> f64 {
    1e-6
}
fn default_initial() -> String {
    "inlet".into()
}

impl Default for Run {
    fn default() -> Self {
        Self {
            steady: true,
            t_end: default_t_end(),
            cfl: default_cfl(),
            tolerance: default_tol(),
            initial: default_initial(),
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Output {
    #[serde(default = "default_dir")]
    pub directory: PathBuf,
    #[serde(default)]
    pub reference: Option<Reference>,
    /// Points where the wall pressure of the nearest body is reported.
    #[serde(default)]
    pub probes: Vec<P3>,
    #[serde(default = "yes")]
    pub vtk: bool,
    /// Interval (time units) between rows of `forces.csv`.
    #[serde(default = "default_every")]
    pub every: f64,
}

fn default_dir() -> PathBuf {
    "out".into()
}
fn default_every() -> f64 {
    1.0
}

impl Default for Output {
    fn default() -> Self {
        Self {
            directory: default_dir(),
            reference: None,
            probes: Vec::new(),
            vtk: true,
            every: default_every(),
        }
    }
}

/// Force coefficients `C = 2F/(ρ U² A)`.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct Reference {
    pub velocity: f64,
    pub area: f64,
}

/// What a run produced.
#[derive(Clone, Debug, Serialize)]
pub struct Summary {
    pub name: String,
    pub cells: [usize; 3],
    pub cell: f64,
    pub steps: usize,
    pub time: f64,
    pub converged: bool,
    pub final_rate: f64,
    pub wall_seconds: f64,
    pub bodies: Vec<BodySummary>,
    /// Wall pressure at each probe (nearest body).
    pub probe_pressure: Vec<f64>,
    pub max_divergence: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct BodySummary {
    pub file: String,
    pub triangles: usize,
    pub watertight: bool,
    /// Force of the fluid on the body from the momentum balance (primary).
    pub force: P3,
    /// The same by surface integration (pressure + viscous), for comparison.
    pub force_surface: P3,
    pub torque: P3,
    pub coefficients: Option<P3>,
}

#[derive(Debug)]
pub enum CaseError {
    Parse(String),
    Geometry(geometry::GeomError),
    Setup(String),
    Io(std::io::Error),
}

impl std::fmt::Display for CaseError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            CaseError::Parse(s) => write!(f, "case file: {s}"),
            CaseError::Geometry(e) => write!(f, "geometry: {e}"),
            CaseError::Setup(s) => write!(f, "setup: {s}"),
            CaseError::Io(e) => write!(f, "i/o: {e}"),
        }
    }
}

impl std::error::Error for CaseError {}

impl From<std::io::Error> for CaseError {
    fn from(e: std::io::Error) -> Self {
        CaseError::Io(e)
    }
}
impl From<geometry::GeomError> for CaseError {
    fn from(e: geometry::GeomError) -> Self {
        CaseError::Geometry(e)
    }
}

impl Case {
    pub fn from_json(text: &str) -> Result<Self, CaseError> {
        serde_json::from_str(text).map_err(|e| CaseError::Parse(e.to_string()))
    }

    pub fn load(path: impl AsRef<Path>) -> Result<(Self, PathBuf), CaseError> {
        let path = path.as_ref();
        let case = Self::from_json(&std::fs::read_to_string(path)?)?;
        let base = path.parent().map(Path::to_path_buf).unwrap_or_default();
        Ok((case, base))
    }

    /// The grid implied by the domain (every extent a whole number of cells).
    pub fn grid(&self) -> Result<Grid, CaseError> {
        let d = &self.domain;
        let mut n = [0usize; 3];
        for a in 0..3 {
            let len = d.max[a] - d.min[a];
            let cells = len / d.cell;
            let r = cells.round();
            if r < 1.0 || (cells - r).abs() > 1e-6 * r.max(1.0) {
                return Err(CaseError::Setup(format!("domain extent {len} along axis {a} is not a whole number of cells of {}", d.cell)));
            }
            n[a] = r as usize;
        }
        Ok(Grid::new(n, d.cell, d.min))
    }

    fn faces(&self) -> [[&Face; 2]; 3] {
        let b = &self.boundaries;
        [[&b.x_min, &b.x_max], [&b.y_min, &b.y_max], [&b.z_min, &b.z_max]]
    }

    fn bcs(&self) -> Result<[[Bc; 2]; 3], CaseError> {
        let (lo, hi) = (self.domain.min, self.domain.max);
        let faces = self.faces();
        let periodic: [bool; 3] = std::array::from_fn(|a| matches!(faces[a][0], Face::Periodic));
        for a in 0..3 {
            if matches!(faces[a][0], Face::Periodic) != matches!(faces[a][1], Face::Periodic) {
                return Err(CaseError::Setup(format!("axis {a}: periodic faces must come in pairs")));
            }
        }
        Ok(std::array::from_fn(|a| {
            std::array::from_fn(|side| match faces[a][side] {
                Face::Periodic => Bc::Periodic,
                Face::Outlet => Bc::Outlet,
                Face::Slip => Bc::Slip,
                Face::Wall { velocity } => Bc::moving_wall(*velocity),
                Face::Inlet { velocity, profile } => {
                    let (v, profile) = (*velocity, *profile);
                    Bc::inlet(move |x: P3, _t| {
                        let mut s = 1.0;
                        if profile == Profile::Parabolic {
                            for b in 0..3 {
                                if b != a && !periodic[b] {
                                    let u = (x[b] - lo[b]) / (hi[b] - lo[b]);
                                    s *= 4.0 * u * (1.0 - u);
                                }
                            }
                        }
                        [v[0] * s, v[1] * s, v[2] * s]
                    })
                }
            })
        }))
    }

    fn load_bodies(&self, base: &Path) -> Result<(Vec<Body>, Vec<(String, usize, bool)>), CaseError> {
        let mut bodies = Vec::new();
        let mut info = Vec::new();
        for spec in &self.bodies {
            let path = base.join(&spec.file);
            let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
            let mut mesh = if matches!(ext.as_str(), "step" | "stp" | "iges" | "igs" | "brep") {
                let opt = Tessellation {
                    elements_per_2pi: spec.elements_per_2pi,
                    max_size: 0.0,
                };
                let mut m = geometry::cad::tessellate(&path, &opt)?;
                m.weld_relative(1e-9);
                m
            } else {
                TriMesh::read(&path)?
            };
            let (s, t) = (spec.scale, spec.translate);
            mesh.transform(|p| [p[0] * s + t[0], p[1] * s + t[1], p[2] * s + t[2]]);
            let r = mesh.report();
            let sdf = if spec.inside_out { Sdf::new_as_oriented(mesh) } else { Sdf::new(mesh) };
            info.push((spec.file.display().to_string(), r.triangles, r.watertight));
            bodies.push(Body {
                sdf,
                velocity: spec.velocity,
                angular: spec.angular_velocity,
                center: spec.center,
            });
        }
        Ok((bodies, info))
    }

    /// Build the solver (grid, boundary conditions, bodies, initial state).
    pub fn build(&self, base: &Path) -> Result<(Solver, Vec<(String, usize, bool)>), CaseError> {
        let grid = self.grid()?;
        let (bodies, info) = self.load_bodies(base)?;
        let mut s = Solver::new(grid, self.fluid.kinematic_viscosity, self.fluid.density, self.bcs()?, bodies);
        if self.run.initial == "inlet" {
            // Extend the first inlet's profile through the domain.
            let faces = self.faces();
            if let Some((a, side)) = (0..3).flat_map(|a| (0..2).map(move |s| (a, s))).find(|&(a, s)| matches!(faces[a][s], Face::Inlet { .. })) {
                let bc = s.bcs[a][side].clone();
                if let Bc::Velocity(f) = bc {
                    let (lo, hi) = (self.domain.min, self.domain.max);
                    s.initialize(move |mut x| {
                        x[a] = if side == 0 { lo[a] } else { hi[a] };
                        f(x, 0.0)
                    });
                }
            }
        } else if self.run.initial != "rest" {
            return Err(CaseError::Setup(format!("unknown initial condition `{}` (inlet|rest)", self.run.initial)));
        }
        Ok((s, info))
    }

    /// Run the case, writing `forces.csv`, `summary.json` and (optionally)
    /// VTK files to the output directory (relative to `base`).
    pub fn run(&self, base: &Path, mut log: impl FnMut(&str)) -> Result<Summary, CaseError> {
        let start = Instant::now();
        let (mut s, info) = self.build(base)?;
        let out = base.join(&self.output.directory);
        std::fs::create_dir_all(&out)?;
        log(&format!(
            "{}: {}×{}×{} cells of {}, {} bodies, setup {:.2?}",
            self.name,
            s.grid.n[0],
            s.grid.n[1],
            s.grid.n[2],
            s.grid.h,
            info.len(),
            start.elapsed()
        ));
        let mut csv = String::from("t");
        for b in 0..s.bodies.len() {
            csv += &format!(",fx{b},fy{b},fz{b}");
        }
        csv.push('\n');
        let mut steps = 0;
        let mut rate = f64::INFINITY;
        let mut next = self.output.every;
        while s.t < self.run.t_end - 1e-12 {
            let target = next.min(self.run.t_end);
            let tol = if self.run.steady { self.run.tolerance } else { 0.0 };
            let (n, r) = s.run_until_steady(self.run.cfl, target, tol);
            steps += n;
            rate = r;
            let mut row = format!("{:.6}", s.t);
            for b in 0..s.bodies.len() {
                let f = s.body_force_momentum(b);
                row += &format!(",{:.9e},{:.9e},{:.9e}", f[0], f[1], f[2]);
            }
            csv += &row;
            csv.push('\n');
            log(&format!("t {:.3}: {row} (rate {r:.2e})", s.t));
            next += self.output.every;
            if self.run.steady && r < self.run.tolerance {
                break;
            }
        }
        std::fs::write(out.join("forces.csv"), csv)?;
        let mut bodies = Vec::new();
        for (b, (file, triangles, watertight)) in info.into_iter().enumerate() {
            let force = s.body_force_momentum(b);
            let fs = s.surface_force(b).total();
            let coefficients = self.output.reference.map(|r| {
                let q = 0.5 * self.fluid.density * r.velocity * r.velocity * r.area;
                [force[0] / q, force[1] / q, force[2] / q]
            });
            if self.output.vtk {
                let u = self.output.reference.map_or(1.0, |r| r.velocity);
                s.write_surface_vtk(b, out.join(format!("body{b}.vtk")), 0.0, u)?;
            }
            bodies.push(BodySummary {
                file,
                triangles,
                watertight,
                force,
                force_surface: fs,
                torque: s.body_torque_momentum(b),
                coefficients,
            });
        }
        if self.output.vtk {
            s.write_vtk(out.join("fields.vtk"))?;
        }
        let probe_pressure = self
            .output
            .probes
            .iter()
            .map(|&x| {
                let b = (0..s.bodies.len())
                    .min_by(|&i, &j| s.bodies[i].sdf.distance(x).abs().total_cmp(&s.bodies[j].sdf.distance(x).abs()))
                    .expect("probes need a body");
                s.wall_pressure(b, x)
            })
            .collect();
        let summary = Summary {
            name: self.name.clone(),
            cells: s.grid.n,
            cell: s.grid.h,
            steps,
            time: s.t,
            converged: rate < self.run.tolerance,
            final_rate: rate,
            wall_seconds: start.elapsed().as_secs_f64(),
            bodies,
            probe_pressure,
            max_divergence: s.max_divergence(),
        };
        std::fs::write(out.join("summary.json"), serde_json::to_string_pretty(&summary).unwrap())?;
        Ok(summary)
    }
}
