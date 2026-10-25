//! Pipe and duct runs: the fluid domain, the metal around it, and the
//! hydraulics of what was drawn.
//!
//! This is the kit that feeds the flow solver most directly. A [`Run`] is a
//! chain of straights and elbows that builds a centreline; a [`Pipe`]
//! carries it as either the **bore** (the fluid domain, which is what an
//! internal-flow case meshes) or the **wall** (the metal, which an
//! external-flow or thermal case wants). The same geometry then answers
//! the hydraulic questions, so the drawn run and the pressure drop quoted
//! for it cannot disagree.
//!
//! ```text
//!   NPS + schedule -> OD, wall  ---> Pipe
//!   Run: straight / elbow chain ---> centreline Path3
//!        Pipe::bore(run)        ---> fluid domain for `immersed`
//!        Pipe::wall(run)        ---> metal, for thermal or external flow
//!        Pipe::pressure_drop()  ---> Darcy-Weisbach on the same numbers
//! ```
//!
//! All lengths are metres.
//!
//! # References and what is exact
//!
//! * **Dimensions** are as tabulated in ASME B36.10M (welded and seamless
//!   wrought steel pipe) for schedules 10, 40, 80 and 160, and ASME B16.9
//!   for the long- and short-radius elbow bend radii (`1.5 D` and `1.0 D`
//!   of the nominal size). Only the sizes in [`NPS_TABLE`] are built in;
//!   anything else is given directly as `Pipe::new(od, wall)`. For a
//!   pressure-containing application, check the figures against the
//!   current edition of the standard: they are reproduced here, not
//!   derived.
//! * **Laminar friction** is `f = 64/Re`, which with Darcy-Weisbach is
//!   algebraically identical to Hagen-Poiseuille `dp = 32 mu L u / D^2`.
//!   The test checks the identity to 1e-15 and then checks the *generated
//!   pipe* against the exact parabolic profile by solving the flow in it.
//! * **Turbulent friction** is the Colebrook-White implicit law
//!   (Colebrook, *J. Inst. Civ. Eng.* 11, 1939), solved by fixed-point
//!   iteration to a stated residual. It is a correlation, not an exact
//!   result: its own scatter against measurement is a few per cent.
//! * **Minor losses** are tabulated `K` factors (Crane Technical Paper
//!   410) and are the crudest numbers here; they are offered as data with
//!   their source, not as predictions.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use geometry::{add, cross, dot, norm, normalize, scale, sub, P3};

use crate::path::{sweep_checked, Path3, Profile, SweepCheck};
use crate::solid::{Placement, Solid};

/// Wall thickness schedule.
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Schedule {
    Sch10,
    Sch40,
    Sch80,
    Sch160,
}

/// `(NPS in inches, OD mm, Sch 10, Sch 40, Sch 80, Sch 160)` wall
/// thicknesses in millimetres, as tabulated in ASME B36.10M. A `None`
/// means the standard does not list that schedule at that size.
pub const NPS_TABLE: &[(f64, f64, Option<f64>, Option<f64>, Option<f64>, Option<f64>)] = &[
    (0.5, 21.3, Some(2.11), Some(2.77), Some(3.73), Some(4.78)),
    (0.75, 26.7, Some(2.11), Some(2.87), Some(3.91), Some(5.56)),
    (1.0, 33.4, Some(2.77), Some(3.38), Some(4.55), Some(6.35)),
    (1.25, 42.2, Some(2.77), Some(3.56), Some(4.85), Some(6.35)),
    (1.5, 48.3, Some(2.77), Some(3.68), Some(5.08), Some(7.14)),
    (2.0, 60.3, Some(2.77), Some(3.91), Some(5.54), Some(8.74)),
    (2.5, 73.0, Some(3.05), Some(5.16), Some(7.01), Some(9.53)),
    (3.0, 88.9, Some(3.05), Some(5.49), Some(7.62), Some(11.13)),
    (4.0, 114.3, Some(3.05), Some(6.02), Some(8.56), Some(13.49)),
    (5.0, 141.3, Some(3.40), Some(6.55), Some(9.53), Some(15.88)),
    (6.0, 168.3, Some(3.40), Some(7.11), Some(10.97), Some(18.26)),
    (8.0, 219.1, Some(3.76), Some(8.18), Some(12.70), Some(23.01)),
    (10.0, 273.0, Some(4.19), Some(9.27), Some(15.09), Some(28.58)),
    (12.0, 323.8, Some(4.57), Some(10.31), Some(17.48), Some(33.32)),
];

/// A pipe cross-section.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct Pipe {
    /// Outside diameter, m.
    pub outer: f64,
    /// Wall thickness, m.
    pub wall: f64,
}

impl Pipe {
    pub fn new(outer: f64, wall: f64) -> Self {
        Self { outer, wall }
    }

    /// A standard size: nominal pipe size in inches and a schedule.
    pub fn nps(nps: f64, schedule: Schedule) -> Result<Self, String> {
        let row = NPS_TABLE
            .iter()
            .find(|r| (r.0 - nps).abs() < 1e-9)
            .ok_or_else(|| format!("NPS {nps} is not in the built-in ASME B36.10M table"))?;
        let wall = match schedule {
            Schedule::Sch10 => row.2,
            Schedule::Sch40 => row.3,
            Schedule::Sch80 => row.4,
            Schedule::Sch160 => row.5,
        }
        .ok_or_else(|| format!("ASME B36.10M lists no {schedule:?} at NPS {nps}"))?;
        Ok(Self { outer: row.1 * 1e-3, wall: wall * 1e-3 })
    }

    /// Inside diameter, m.
    pub fn inner(&self) -> f64 {
        self.outer - 2.0 * self.wall
    }

    /// Bore cross-sectional area, m².
    pub fn flow_area(&self) -> f64 {
        0.25 * PI * self.inner().powi(2)
    }

    /// Metal cross-sectional area, m².
    pub fn metal_area(&self) -> f64 {
        0.25 * PI * (self.outer.powi(2) - self.inner().powi(2))
    }

    /// Mass per unit length for a given material density, kg/m.
    pub fn mass_per_length(&self, density: f64) -> f64 {
        self.metal_area() * density
    }

    /// Hydraulic diameter `4A/P`, which for a round bore is the bore.
    pub fn hydraulic_diameter(&self) -> f64 {
        4.0 * self.flow_area() / (PI * self.inner())
    }

    /// Reynolds number from a bulk velocity and a kinematic viscosity.
    pub fn reynolds(&self, velocity: f64, kinematic_viscosity: f64) -> f64 {
        velocity * self.inner() / kinematic_viscosity
    }

    /// Bulk velocity for a volumetric flow rate, m³/s.
    pub fn velocity(&self, flow_rate: f64) -> f64 {
        flow_rate / self.flow_area()
    }

    /// Darcy friction factor. Laminar below `Re = 2300`, Colebrook-White
    /// above it, and nothing sensible in between — the transition is
    /// reported rather than interpolated.
    pub fn friction_factor(&self, re: f64, roughness: f64) -> Result<f64, String> {
        if re <= 0.0 {
            return Err("Reynolds number must be positive".into());
        }
        if re < 2300.0 {
            return Ok(64.0 / re);
        }
        if re < 4000.0 {
            return Err(format!(
                "Re = {re:.0} is in the transition range 2300-4000, where neither the laminar law \
                 nor Colebrook applies; no friction factor is defined"
            ));
        }
        Ok(colebrook(re, roughness / self.inner()))
    }

    /// Pressure drop over `length` by Darcy-Weisbach,
    /// `dp = f (L/D) (rho u^2 / 2)`.
    pub fn pressure_drop(
        &self,
        length: f64,
        velocity: f64,
        density: f64,
        kinematic_viscosity: f64,
        roughness: f64,
    ) -> Result<f64, String> {
        let re = self.reynolds(velocity, kinematic_viscosity);
        let f = self.friction_factor(re, roughness)?;
        Ok(f * length / self.inner() * 0.5 * density * velocity * velocity)
    }

    /// The laminar pressure drop straight from Hagen-Poiseuille,
    /// `32 mu L u / D^2`. Identical to [`Pipe::pressure_drop`] below
    /// `Re = 2300`; kept separate so the test can check the identity.
    pub fn poiseuille_drop(&self, length: f64, velocity: f64, density: f64, kinematic_viscosity: f64) -> f64 {
        32.0 * density * kinematic_viscosity * length * velocity / self.inner().powi(2)
    }

    /// Centreline velocity of fully developed laminar flow, exactly twice
    /// the bulk velocity.
    pub fn poiseuille_centreline(&self, velocity: f64) -> f64 {
        2.0 * velocity
    }

    /// The fully developed laminar profile at radius `r`.
    pub fn poiseuille_profile(&self, velocity: f64, r: f64) -> f64 {
        let a = 0.5 * self.inner();
        2.0 * velocity * (1.0 - (r / a).powi(2)).max(0.0)
    }

    /// The bore as a solid: the fluid domain of this run.
    pub fn bore(&self, run: &Path3, sides: usize) -> (Solid, SweepCheck) {
        let (mesh, check) = sweep_checked(run, &Profile::circle(0.5 * self.inner(), sides));
        (Solid::from_trimesh(mesh), check)
    }

    /// The metal as a solid: the swept outside diameter less the bore.
    pub fn wall_solid(&self, run: &Path3, sides: usize) -> (Solid, SweepCheck) {
        let (outer, check) = sweep_checked(run, &Profile::circle(0.5 * self.outer, sides));
        let (inner, _) = sweep_checked(run, &Profile::circle(0.5 * self.inner(), sides));
        (
            Solid::from_trimesh(outer).minus(Solid::from_trimesh(inner)),
            check,
        )
    }

    /// The bore as an inside-out solid, so a flow solver treats the
    /// exterior as the wall: what `immersed` wants for internal flow.
    pub fn bore_as_housing(&self, run: &Path3, sides: usize) -> Result<geometry::Sdf, String> {
        let (mesh, check) = sweep_checked(run, &Profile::circle(0.5 * self.inner(), sides));
        if !check.watertight {
            return Err("the swept bore is not watertight".into());
        }
        if check.curvature_ratio >= 1.0 {
            return Err(format!(
                "the run bends tighter than the bore: the section reaches {:.2} of the way past the \
                 centre of curvature at station {}",
                check.curvature_ratio, check.worst_station
            ));
        }
        let mut mesh = mesh;
        mesh.orient_outward();
        mesh.flip();
        Ok(geometry::Sdf::new_as_oriented(mesh))
    }
}

/// Colebrook-White, solved by fixed-point iteration on `1/sqrt(f)`.
pub fn colebrook(re: f64, relative_roughness: f64) -> f64 {
    // Swamee-Jain as the starting guess, then Colebrook to convergence.
    let mut x = -2.0 * (relative_roughness / 3.7 + 5.74 / re.powf(0.9)).log10();
    for _ in 0..200 {
        let next = -2.0 * (relative_roughness / 3.7 + 2.51 * x / re).log10();
        if (next - x).abs() < 1e-15 * next.abs().max(1.0) {
            x = next;
            break;
        }
        x = next;
    }
    1.0 / (x * x)
}

/// Residual of the Colebrook equation, for checking a solution.
pub fn colebrook_residual(f: f64, re: f64, relative_roughness: f64) -> f64 {
    let x = 1.0 / f.sqrt();
    x + 2.0 * (relative_roughness / 3.7 + 2.51 * x / re).log10()
}

/// Minor-loss coefficients `K`, as tabulated in Crane Technical Paper 410.
/// These are handbook data with a spread of tens of per cent between
/// sources, not predictions.
pub mod minor_loss {
    /// 90° long-radius elbow (r/D = 1.5).
    pub const ELBOW_90_LONG: f64 = 0.3;
    /// 90° short-radius (standard) elbow.
    pub const ELBOW_90_SHORT: f64 = 0.9;
    /// 45° elbow.
    pub const ELBOW_45: f64 = 0.4;
    /// Tee, flow through the run.
    pub const TEE_RUN: f64 = 0.2;
    /// Tee, flow turning into the branch.
    pub const TEE_BRANCH: f64 = 1.0;
    /// Sudden entry from a large vessel.
    pub const ENTRY_SHARP: f64 = 0.5;
    /// Discharge to a large vessel.
    pub const EXIT: f64 = 1.0;
    /// Fully open gate valve.
    pub const GATE_VALVE_OPEN: f64 = 0.15;
    /// Fully open globe valve.
    pub const GLOBE_VALVE_OPEN: f64 = 10.0;
}

/// One leg of a run.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
#[serde(tag = "leg", rename_all = "snake_case")]
pub enum Leg {
    Straight { length: f64 },
    /// Turn `angle` radians about `axis`, which must be perpendicular to
    /// the current direction, on a centreline radius `radius`.
    Elbow { radius: f64, angle: f64, axis: P3 },
}

/// A chain of straights and elbows, built into a centreline.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Run {
    pub start: P3,
    pub direction: P3,
    pub legs: Vec<Leg>,
    /// Points per metre of straight leg, and per radian of elbow.
    #[serde(default = "straight_steps")]
    pub straight_steps: usize,
    #[serde(default = "elbow_steps")]
    pub elbow_steps: usize,
}

fn straight_steps() -> usize {
    4
}
fn elbow_steps() -> usize {
    24
}

impl Run {
    pub fn new(start: P3, direction: P3) -> Self {
        Self {
            start,
            direction: normalize(direction),
            legs: Vec::new(),
            straight_steps: straight_steps(),
            elbow_steps: elbow_steps(),
        }
    }

    pub fn straight(mut self, length: f64) -> Self {
        self.legs.push(Leg::Straight { length });
        self
    }

    pub fn elbow(mut self, radius: f64, angle: f64, axis: P3) -> Self {
        self.legs.push(Leg::Elbow { radius, angle, axis });
        self
    }

    /// A long-radius 90° elbow, ASME B16.9: centreline radius `1.5 D` of
    /// the nominal size.
    pub fn elbow_long(self, nominal: f64, axis: P3) -> Self {
        self.elbow(1.5 * nominal, 0.5 * PI, axis)
    }

    /// A short-radius 90° elbow, ASME B16.9: centreline radius `1.0 D`.
    pub fn elbow_short(self, nominal: f64, axis: P3) -> Self {
        self.elbow(nominal, 0.5 * PI, axis)
    }

    /// Total centreline length, exactly (straights plus `R θ` arcs).
    pub fn length(&self) -> f64 {
        self.legs
            .iter()
            .map(|leg| match *leg {
                Leg::Straight { length } => length,
                Leg::Elbow { radius, angle, .. } => radius * angle.abs(),
            })
            .sum()
    }

    /// Build the centreline, checking that every elbow axis is
    /// perpendicular to the direction it turns.
    pub fn centreline(&self) -> Result<Path3, String> {
        let mut points = vec![self.start];
        let mut p = self.start;
        let mut d = normalize(self.direction);
        for (n, leg) in self.legs.iter().enumerate() {
            match *leg {
                Leg::Straight { length } => {
                    if length <= 0.0 {
                        return Err(format!("leg {n}: a straight needs a positive length"));
                    }
                    let steps = ((self.straight_steps as f64 * length / length.max(1e-12)) as usize)
                        .max(self.straight_steps)
                        .max(1);
                    for i in 1..=steps {
                        points.push(add(p, scale(d, length * i as f64 / steps as f64)));
                    }
                    p = add(p, scale(d, length));
                }
                Leg::Elbow { radius, angle, axis } => {
                    if radius <= 0.0 {
                        return Err(format!("leg {n}: an elbow needs a positive radius"));
                    }
                    let a = normalize(axis);
                    if norm(axis) < 1e-12 {
                        return Err(format!("leg {n}: the elbow axis is zero"));
                    }
                    if dot(a, d).abs() > 1e-9 {
                        return Err(format!(
                            "leg {n}: the elbow axis is not perpendicular to the run \
                             (dot = {:.3e}); a pipe cannot twist about its own axis",
                            dot(a, d)
                        ));
                    }
                    // Turning about +a with the centre offset along a x d
                    // keeps the motion going forward along d.
                    let centre = add(p, scale(cross(a, d), radius));
                    let steps = ((self.elbow_steps as f64 * angle.abs() / (0.5 * PI)).ceil() as usize).max(2);
                    let r0 = sub(p, centre);
                    for i in 1..=steps {
                        let t = angle * i as f64 / steps as f64;
                        let rot = Placement::rotation(a, t);
                        points.push(add(centre, rot.apply(r0)));
                    }
                    let rot = Placement::rotation(a, angle);
                    p = add(centre, rot.apply(r0));
                    d = normalize(rot.apply(d));
                }
            }
        }
        if points.len() < 2 {
            return Err("a run needs at least one leg".into());
        }
        Ok(Path3::new(points, false))
    }

    /// Where the run ends, and which way it is pointing.
    pub fn end(&self) -> Result<(P3, P3), String> {
        let path = self.centreline()?;
        let n = path.points.len();
        Ok((path.points[n - 1], normalize(sub(path.points[n - 1], path.points[n - 2]))))
    }
}

/// A conical reducer between two pipes.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct Reducer {
    pub from: Pipe,
    pub to: Pipe,
    pub length: f64,
}

impl Reducer {
    /// The bore: a cone between the two inside diameters, along `+z` from
    /// the origin.
    pub fn bore(&self) -> Solid {
        Solid::cone(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, self.length],
            0.5 * self.from.inner(),
            0.5 * self.to.inner(),
        )
    }

    /// The metal.
    pub fn wall_solid(&self) -> Solid {
        Solid::cone(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, self.length],
            0.5 * self.from.outer,
            0.5 * self.to.outer,
        )
        .minus(Solid::cone(
            [0.0, 0.0, -self.length],
            [0.0, 0.0, 2.0 * self.length],
            0.5 * self.from.inner(),
            0.5 * self.to.inner(),
        ))
    }

    /// Half-angle of the taper, radians. Above about 0.12 rad (7°) a
    /// diffuser separates, so the sign of this number decides whether a
    /// Bernoulli estimate of the pressure recovery means anything.
    pub fn half_angle(&self) -> f64 {
        (0.5 * (self.to.inner() - self.from.inner()) / self.length).atan()
    }

    /// Ideal (loss-free) pressure change from continuity and Bernoulli:
    /// `dp = rho/2 (u1^2 - u2^2)` with `u2 = u1 A1/A2`.
    pub fn bernoulli_recovery(&self, velocity_in: f64, density: f64) -> f64 {
        let u2 = velocity_in * self.from.flow_area() / self.to.flow_area();
        0.5 * density * (velocity_in * velocity_in - u2 * u2)
    }
}

/// An equal or reducing tee: a run with a branch at `angle` to it.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct Tee {
    pub run: Pipe,
    pub branch: Pipe,
    /// Length of the run, centred on the branch.
    pub run_length: f64,
    /// Length of the branch from the run's centreline.
    pub branch_length: f64,
}

impl Tee {
    /// The bore: the run's bore united with the branch's, which is the
    /// fluid domain of the fitting.
    pub fn bore(&self) -> Solid {
        let half = 0.5 * self.run_length;
        Solid::union(vec![
            Solid::cylinder([-half, 0.0, 0.0], [half, 0.0, 0.0], 0.5 * self.run.inner()),
            Solid::cylinder([0.0, 0.0, 0.0], [0.0, 0.0, self.branch_length], 0.5 * self.branch.inner()),
        ])
    }

    /// The metal: the two outside diameters less the two bores.
    pub fn wall_solid(&self) -> Solid {
        let half = 0.5 * self.run_length;
        let outer = Solid::union(vec![
            Solid::cylinder([-half, 0.0, 0.0], [half, 0.0, 0.0], 0.5 * self.run.outer),
            Solid::cylinder([0.0, 0.0, 0.0], [0.0, 0.0, self.branch_length], 0.5 * self.branch.outer),
        ]);
        outer.minus(self.bore())
    }
}
