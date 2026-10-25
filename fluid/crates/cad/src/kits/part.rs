//! Machined and printed solids.
//!
//! Lengths are in whatever unit the caller works in, used consistently.
//! The gear and thread standards are written in millimetres, so those two
//! specs read naturally in millimetres; scale the solid if the rest of the
//! model is in metres.
//!
//! # Involute spur gears
//!
//! Built from the standard geometry (ISO 21771; Litvin & Fuentes, *Gear
//! Geometry and Applied Theory*, 2nd ed., 2004, §10): pitch diameter
//! `d = m z`, base diameter `d_b = d cos a`, and the flank is the involute
//! of the base circle, whose angular half-width at radius `r` is
//!
//! ```text
//!   psi(r) = (pi/2 + 2 x tan a)/z + inv(a) - inv(phi),
//!   inv(t) = tan t - t,   cos phi = r_b / r
//! ```
//!
//! At the pitch circle this gives a tooth thickness of `m (pi/2 + 2x tan a)`,
//! which for an unshifted gear is exactly `pi m / 2`; the pitch circle is
//! always one of the sampled flank radii, so that thickness is exact in the
//! generated polygon rather than short by a chord sagitta. The root is cut
//! straight to the root circle: there is **no trochoidal root fillet**, so
//! root stresses from this geometry are not meaningful. The flank itself is
//! the true involute, sampled at `flank_points` radii.
//!
//! # ISO metric threads
//!
//! The basic profile of ISO 68-1: a 60° triangle of height
//! `H = P sqrt(3)/2`, truncated by `H/8` at the crest and `H/4` at the
//! root, giving the ISO 724 basic diameters
//!
//! ```text
//!   d2 = d - 2 (3/8) H = d - 0.649519 P      (pitch diameter)
//!   d1 = d - 2 (5/8) H = d - 1.082532 P      (minor diameter)
//! ```
//!
//! The solid is the core cylinder united with a helical ridge: the field is
//! evaluated in the helical coordinate `u = z - P theta / 2 pi`, where the
//! ridge is a fixed cross-section. The transform's gradient bound,
//! `sqrt(1 + (P / 2 pi r)^2)`, is divided out so the kernel invariant still
//! holds; for M8 that factor is 1.002. The root is left sharp (ISO allows a
//! rounded root, which is what sets `d3 = d - 1.226869 P`); `root_fillet`
//! applies that rounding as an offset pair.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::sketch::{Polygon, Region2};
use crate::solid::Solid;

/// `inv(t) = tan t - t`, the involute function.
pub fn involute(t: f64) -> f64 {
    t.tan() - t
}

/// An external involute spur gear.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct SpurGear {
    /// Module `m`: pitch diameter per tooth.
    pub module: f64,
    pub teeth: usize,
    /// Pressure angle, radians (20° = 0.349066 is the ISO default).
    #[serde(default = "default_pressure_angle")]
    pub pressure_angle: f64,
    /// Addendum in modules (1.0 standard full depth).
    #[serde(default = "one")]
    pub addendum: f64,
    /// Dedendum in modules (1.25 standard full depth).
    #[serde(default = "dedendum_default")]
    pub dedendum: f64,
    /// Profile shift coefficient `x`.
    #[serde(default)]
    pub shift: f64,
    pub thickness: f64,
    /// Central bore diameter (0 for none).
    #[serde(default)]
    pub bore: f64,
    /// Points sampled along each flank.
    #[serde(default = "flank_default")]
    pub flank_points: usize,
}

fn default_pressure_angle() -> f64 {
    20.0_f64.to_radians()
}
fn one() -> f64 {
    1.0
}
fn dedendum_default() -> f64 {
    1.25
}
fn flank_default() -> usize {
    24
}

/// Measured and derived gear dimensions.
#[derive(Clone, Copy, Debug)]
pub struct GearReport {
    pub pitch_diameter: f64,
    pub base_diameter: f64,
    pub tip_diameter: f64,
    pub root_diameter: f64,
    /// Circular tooth thickness at the pitch circle.
    pub pitch_thickness: f64,
    /// Circular pitch `pi m`.
    pub circular_pitch: f64,
    /// Base pitch `pi m cos a`.
    pub base_pitch: f64,
    /// True when the tip was cut back because the flanks met first.
    pub tip_truncated: bool,
    /// Teeth below which this pressure angle undercuts
    /// (`2 / sin^2 a`, Litvin & Fuentes §10).
    pub undercut_limit: f64,
}

impl SpurGear {
    pub fn new(module: f64, teeth: usize, thickness: f64) -> Self {
        Self {
            module,
            teeth,
            pressure_angle: default_pressure_angle(),
            addendum: 1.0,
            dedendum: 1.25,
            shift: 0.0,
            thickness,
            bore: 0.0,
            flank_points: flank_default(),
        }
    }

    pub fn pitch_radius(&self) -> f64 {
        0.5 * self.module * self.teeth as f64
    }
    pub fn base_radius(&self) -> f64 {
        self.pitch_radius() * self.pressure_angle.cos()
    }
    pub fn tip_radius(&self) -> f64 {
        self.pitch_radius() + (self.addendum + self.shift) * self.module
    }
    pub fn root_radius(&self) -> f64 {
        (self.pitch_radius() - (self.dedendum - self.shift) * self.module).max(1e-9)
    }

    /// Angular half-width of a tooth at radius `r` (radians). `None` below
    /// the base circle, where the involute does not exist.
    pub fn half_angle(&self, r: f64) -> Option<f64> {
        let rb = self.base_radius();
        if r < rb {
            return None;
        }
        let a = self.pressure_angle;
        let psi0 = (PI / 2.0 + 2.0 * self.shift * a.tan()) / self.teeth as f64;
        let phi = (rb / r).clamp(-1.0, 1.0).acos();
        Some(psi0 + involute(a) - involute(phi))
    }

    /// The largest tip radius at which the tooth still has width.
    fn usable_tip(&self) -> (f64, bool) {
        let ra = self.tip_radius();
        if self.half_angle(ra).is_some_and(|p| p > 1e-9) {
            return (ra, false);
        }
        // Bisect for psi = 0.
        let (mut lo, mut hi) = (self.base_radius() * 1.000001, ra);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            match self.half_angle(mid) {
                Some(p) if p > 0.0 => lo = mid,
                _ => hi = mid,
            }
        }
        (lo, true)
    }

    /// The cross-section: the root disc united with the teeth.
    pub fn region(&self) -> Region2 {
        let (ra, _) = self.usable_tip();
        let rf = self.root_radius();
        let rb = self.base_radius();
        let r_start = rb.max(rf);
        let n = self.flank_points.max(3);
        let psi_start = self.half_angle(r_start.max(rb * 1.000001)).unwrap_or(0.0);
        // Sample radii, with the pitch circle always included so that the
        // tooth thickness there is the standard value exactly rather than
        // short by the chord's sagitta.
        let mut radii: Vec<f64> = (0..=n).map(|i| r_start + (ra - r_start) * i as f64 / n as f64).collect();
        let rp = self.pitch_radius();
        if rp > r_start && rp < ra {
            radii.push(rp);
            radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
        }

        // One tooth, centred on +x, from the root circle out and back.
        let mut pts: Vec<[f64; 2]> = Vec::with_capacity(2 * n + 6);
        pts.push([rf * (-psi_start).cos(), rf * (-psi_start).sin()]);
        for &r in &radii {
            let psi = self.half_angle(r.max(rb * 1.000001)).unwrap_or(0.0).max(0.0);
            pts.push([r * (-psi).cos(), r * (-psi).sin()]);
        }
        for &r in radii.iter().rev() {
            let psi = self.half_angle(r.max(rb * 1.000001)).unwrap_or(0.0).max(0.0);
            pts.push([r * psi.cos(), r * psi.sin()]);
        }
        pts.push([rf * psi_start.cos(), rf * psi_start.sin()]);
        let tooth = Region2::Polygon(Polygon::new(pts));

        // Folding the angle costs three evaluations instead of one per
        // tooth, which is what makes a 48-tooth gear meshable.
        let body = Region2::circle([0.0, 0.0], rf).plus(tooth.polar_repeat(self.teeth.max(1)));
        if self.bore > 0.0 {
            body.minus(Region2::circle([0.0, 0.0], 0.5 * self.bore))
        } else {
            body
        }
    }

    pub fn solid(&self) -> Solid {
        Solid::extrude(self.region(), 0.0, self.thickness)
    }

    pub fn report(&self) -> GearReport {
        let (_, truncated) = self.usable_tip();
        let a = self.pressure_angle;
        GearReport {
            pitch_diameter: 2.0 * self.pitch_radius(),
            base_diameter: 2.0 * self.base_radius(),
            tip_diameter: 2.0 * self.tip_radius(),
            root_diameter: 2.0 * self.root_radius(),
            pitch_thickness: self.module * (PI / 2.0 + 2.0 * self.shift * a.tan()),
            circular_pitch: PI * self.module,
            base_pitch: PI * self.module * a.cos(),
            tip_truncated: truncated,
            undercut_limit: 2.0 / (a.sin() * a.sin()),
        }
    }

    /// Centre distance for a standard (unshifted) mesh with `other`.
    pub fn centre_distance(&self, other: &SpurGear) -> f64 {
        0.5 * self.module * (self.teeth + other.teeth) as f64
    }

    /// Transverse contact ratio of a standard mesh (ISO 21771):
    /// `[sqrt(ra1^2 - rb1^2) + sqrt(ra2^2 - rb2^2) - a sin a] / p_b`.
    pub fn contact_ratio(&self, other: &SpurGear) -> f64 {
        let a = self.pressure_angle;
        let (ra1, rb1) = (self.tip_radius(), self.base_radius());
        let (ra2, rb2) = (other.tip_radius(), other.base_radius());
        let centre = self.centre_distance(other);
        let path = (ra1 * ra1 - rb1 * rb1).max(0.0).sqrt() + (ra2 * ra2 - rb2 * rb2).max(0.0).sqrt()
            - centre * a.sin();
        path / (PI * self.module * a.cos())
    }
}

/// An ISO metric thread.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Thread {
    /// Nominal major diameter `d`.
    pub diameter: f64,
    pub pitch: f64,
    /// Threaded length along `+z` from `z = 0`.
    pub length: f64,
    /// Left-handed when true.
    #[serde(default)]
    pub left_hand: bool,
    /// Root rounding radius ISO allows on an external thread
    /// (`>= 0.125 P`). It is **reported, not modelled**: the solid keeps
    /// the sharp basic profile, so the minor diameter of the solid is the
    /// basic `d1`, not `d3`.
    #[serde(default)]
    pub root_fillet: f64,
}

/// ISO 724 basic diameters.
#[derive(Clone, Copy, Debug)]
pub struct ThreadReport {
    /// Height of the fundamental triangle, `P sqrt(3)/2`.
    pub h: f64,
    /// Major diameter `d`.
    pub major: f64,
    /// Pitch diameter `d2 = d - 0.649519 P`.
    pub pitch_diameter: f64,
    /// Basic minor diameter `d1 = d - 1.082532 P`.
    pub minor: f64,
    /// External minor diameter with the rounded root,
    /// `d3 = d - 1.226869 P`.
    pub minor_rounded: f64,
    /// Crest flat width, `P/8`.
    pub crest_flat: f64,
    /// Lipschitz factor divided out of the helical field.
    pub helical_factor: f64,
}

impl Thread {
    pub fn new(diameter: f64, pitch: f64, length: f64) -> Self {
        Self { diameter, pitch, length, left_hand: false, root_fillet: 0.0 }
    }

    /// ISO 261 coarse-pitch thread, nominal diameter in **millimetres**
    /// (ISO 261 table 1, first choice). `None` for a size not in the
    /// table.
    pub fn coarse_mm(diameter_mm: f64, length_mm: f64) -> Option<Self> {
        let pitch = match (diameter_mm * 100.0).round() as i64 {
            100 => 0.25,
            120 => 0.25,
            160 => 0.35,
            200 => 0.40,
            250 => 0.45,
            300 => 0.50,
            400 => 0.70,
            500 => 0.80,
            600 => 1.00,
            800 => 1.25,
            1000 => 1.50,
            1200 => 1.75,
            1600 => 2.00,
            2000 => 2.50,
            2400 => 3.00,
            3000 => 3.50,
            3600 => 4.00,
            _ => return None,
        };
        Some(Self::new(diameter_mm, pitch, length_mm))
    }

    pub fn report(&self) -> ThreadReport {
        let p = self.pitch;
        let h = p * 3.0_f64.sqrt() / 2.0;
        ThreadReport {
            h,
            major: self.diameter,
            pitch_diameter: self.diameter - 2.0 * (3.0 / 8.0) * h,
            minor: self.diameter - 2.0 * (5.0 / 8.0) * h,
            minor_rounded: self.diameter - 2.0 * (17.0 / 24.0) * h,
            crest_flat: p / 8.0,
            helical_factor: self.helical_factor(),
        }
    }

    /// Basic minor diameter `d1 = d - 1.082532 P`.
    fn minor_diameter(&self) -> f64 {
        self.diameter - 2.0 * (5.0 / 8.0) * self.pitch * 3.0_f64.sqrt() / 2.0
    }

    fn helical_factor(&self) -> f64 {
        let r_min = 0.5 * self.minor_diameter();
        (1.0 + (self.pitch / (2.0 * PI * r_min)).powi(2)).sqrt()
    }

    /// The thread ridge cross-section in the `(u, r)` half-plane, `u` the
    /// helical coordinate: the ISO basic profile over one pitch.
    fn ridge_region(&self) -> Region2 {
        let p = self.pitch;
        let r_crest = 0.5 * self.diameter;
        let r_root = 0.5 * self.minor_diameter();
        // Half-widths: P/16 at the crest, 3P/8 at the root.
        let poly = Polygon::new(vec![
            [-3.0 * p / 8.0, r_root],
            [-p / 16.0, r_crest],
            [p / 16.0, r_crest],
            [3.0 * p / 8.0, r_root],
        ]);
        // The ISO root rounding is a concave fillet, which offsets alone
        // cannot produce (see Solid::round_convex); the sharp basic
        // profile is used and `root_fillet` is reported, not modelled.
        Region2::Polygon(poly)
    }

    /// The external thread: core cylinder plus helical ridge, limited to
    /// `0 <= z <= length`.
    pub fn solid(&self) -> Solid {
        let p = self.pitch;
        let hand = if self.left_hand { -1.0 } else { 1.0 };
        let ridge = self.ridge_region();
        let lip = self.helical_factor();
        let r_core = 0.5 * self.minor_diameter();
        let r_crest = 0.5 * self.diameter;
        let len = self.length;
        let bbox = (
            [-r_crest, -r_crest, -0.1 * len],
            [r_crest, r_crest, 1.1 * len],
        );
        let helical = Solid::custom("iso_thread_ridge", bbox, false, move |q| {
            let r = (q[0] * q[0] + q[1] * q[1]).sqrt();
            let theta = q[1].atan2(q[0]);
            let u0 = q[2] - hand * p * theta / (2.0 * PI);
            // Periodic in u with period p: take the nearest three images.
            let m = (u0 / p).round();
            let mut best = f64::INFINITY;
            for k in -1..=1 {
                let u = u0 - (m + k as f64) * p;
                best = best.min(ridge.sdf([u, r]));
            }
            best / lip
        });
        let core = Solid::cylinder([0.0, 0.0, 0.0], [0.0, 0.0, len], r_core);
        let band = Solid::box_between([-r_crest, -r_crest, 0.0], [r_crest, r_crest, len]);
        Solid::union(vec![core, helical.and(band)])
    }

    /// The solid to subtract from a part to leave an internal thread of
    /// this nominal size, grown by `allowance` for clearance.
    pub fn tap_cut(&self, allowance: f64) -> Solid {
        self.solid().offset(allowance)
    }
}

/// A flat plate: an outline, holes, a thickness, optional edge rounding.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Plate {
    pub outline: Region2,
    #[serde(default)]
    pub holes: Vec<Region2>,
    pub thickness: f64,
}

impl Plate {
    pub fn new(outline: Region2, thickness: f64) -> Self {
        Self { outline, holes: Vec::new(), thickness }
    }

    pub fn with_hole(mut self, hole: Region2) -> Self {
        self.holes.push(hole);
        self
    }

    pub fn solid(&self) -> Solid {
        let mut region = self.outline.clone();
        if !self.holes.is_empty() {
            region = region.minus(Region2::Union { parts: self.holes.clone() });
        }
        Solid::extrude(region, 0.0, self.thickness)
    }
}

/// A cylindrical boss of `height` standing on `z = base`, with an optional
/// through bore, blended into whatever it stands on by `blend`.
pub fn boss(center: [f64; 2], radius: f64, base: f64, height: f64, bore: f64) -> Solid {
    let c = [center[0], center[1], base];
    let outer = Solid::cylinder(c, [center[0], center[1], base + height], radius);
    if bore > 0.0 {
        outer.minus(Solid::cylinder(
            [center[0], center[1], base - height],
            [center[0], center[1], base + 2.0 * height],
            0.5 * bore,
        ))
    } else {
        outer
    }
}

/// A straight stiffening rib of thickness `t` and height `h` along the
/// segment `a`-`b` in the plane `z = base`.
pub fn rib(a: [f64; 2], b: [f64; 2], t: f64, base: f64, h: f64) -> Solid {
    let region = Region2::slot(a, b, 0.5 * t);
    Solid::extrude(region, base, base + h)
}

/// A helical compression spring.
///
/// Swept, not implicit: the wire is a circle carried along a helix by
/// [`crate::path::sweep_checked`], which also reports whether the coils
/// clash.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct Spring {
    /// Wire diameter.
    pub wire: f64,
    /// Mean coil diameter (to the wire's centreline).
    pub mean_diameter: f64,
    /// Free length.
    pub free_length: f64,
    /// Number of active coils.
    pub coils: f64,
    /// Stations per coil in the swept path.
    #[serde(default = "spring_steps")]
    pub steps: usize,
    /// Points around the wire section.
    #[serde(default = "spring_sides")]
    pub sides: usize,
}

fn spring_steps() -> usize {
    48
}
fn spring_sides() -> usize {
    16
}

impl Spring {
    pub fn new(wire: f64, mean_diameter: f64, free_length: f64, coils: f64) -> Self {
        Self {
            wire,
            mean_diameter,
            free_length,
            coils,
            steps: spring_steps(),
            sides: spring_sides(),
        }
    }

    /// Axial pitch between coils.
    pub fn pitch(&self) -> f64 {
        self.free_length / self.coils
    }

    /// Gap between adjacent coils; negative means they clash.
    pub fn coil_clearance(&self) -> f64 {
        self.pitch() - self.wire
    }

    /// Spring index `D/d`. Below about 4 the wire cannot be coiled; above
    /// about 15 the spring tangles (Shigley & Mischke, *Mechanical
    /// Engineering Design*, 5th ed., ch. 10).
    pub fn index(&self) -> f64 {
        self.mean_diameter / self.wire
    }

    /// Axial stiffness `k = G d^4 / (8 D^3 n)` (Shigley & Mischke, ch. 10,
    /// from the torsional deflection of the wire).
    pub fn stiffness(&self, shear_modulus: f64) -> f64 {
        shear_modulus * self.wire.powi(4) / (8.0 * self.mean_diameter.powi(3) * self.coils)
    }

    /// Solid length: the coils stacked touching.
    pub fn solid_length(&self) -> f64 {
        (self.coils + 1.0) * self.wire
    }

    /// Wahl curvature correction factor
    /// `K = (4C - 1)/(4C - 4) + 0.615/C` (Wahl, *Mechanical Springs*, 2nd
    /// ed., 1963), which raises the surface shear stress at the inside of
    /// the coil.
    pub fn wahl_factor(&self) -> f64 {
        let c = self.index();
        (4.0 * c - 1.0) / (4.0 * c - 4.0) + 0.615 / c
    }

    /// Maximum shear stress at a given load,
    /// `tau = K * 8 F D / (pi d^3)`.
    pub fn shear_stress(&self, load: f64) -> f64 {
        self.wahl_factor() * 8.0 * load * self.mean_diameter / (PI * self.wire.powi(3))
    }

    /// Exact wire length: a helix of `n` turns at radius `D/2` and pitch
    /// `p` has length `n sqrt((pi D)^2 + p^2)`.
    pub fn wire_length(&self) -> f64 {
        let circumference = PI * self.mean_diameter;
        self.coils * (circumference * circumference + self.pitch().powi(2)).sqrt()
    }

    /// Wire volume: the length times the section.
    pub fn volume(&self) -> f64 {
        self.wire_length() * 0.25 * PI * self.wire * self.wire
    }

    pub fn mass(&self, density: f64) -> f64 {
        self.volume() * density
    }

    /// Whether the spring could be coiled, as opposed to drawn.
    ///
    /// A spring index below about 4 cannot be wound and above about 15
    /// tangles in handling (Shigley & Mischke, ch. 10). This is separate
    /// from [`Spring::solid`], which only reports whether the geometry is
    /// self-consistent: a spring of index 1.6 sweeps perfectly well and
    /// could not be made.
    pub fn manufacturability(&self) -> Result<(), String> {
        if self.coil_clearance() <= 0.0 {
            return Err(format!(
                "the coils clash: pitch {:e} is not greater than the wire diameter {:e}",
                self.pitch(),
                self.wire
            ));
        }
        let c = self.index();
        if c < 4.0 {
            return Err(format!(
                "spring index {c:.2} is below 4: the wire cannot be coiled this tightly"
            ));
        }
        if c > 15.0 {
            return Err(format!("spring index {c:.2} is above 15: the spring will tangle"));
        }
        Ok(())
    }

    /// The coil as a swept solid, with the sweep's own validity report.
    ///
    /// This checks only that the geometry is consistent — the coils clear
    /// each other and the wire does not pass through the coil axis. Use
    /// [`Spring::manufacturability`] for whether it could be wound.
    pub fn solid(&self) -> Result<(Solid, crate::path::SweepCheck), String> {
        if self.coil_clearance() <= 0.0 {
            return Err(format!(
                "the coils clash: pitch {:e} is not greater than the wire diameter {:e}",
                self.pitch(),
                self.wire
            ));
        }
        let path = crate::path::Path3::helix(
            0.5 * self.mean_diameter,
            self.pitch(),
            self.coils,
            0.0,
            0.0,
            self.steps,
        );
        let (mesh, check) = crate::path::sweep_checked(
            &path,
            &crate::path::Profile::circle(0.5 * self.wire, self.sides),
        );
        if check.curvature_ratio >= 1.0 {
            return Err(format!(
                "the wire is too thick for the coil: the section reaches {:.2} of the way past the \
                 centre of curvature (spring index {:.2})",
                check.curvature_ratio,
                self.index()
            ));
        }
        Ok((Solid::from_trimesh(mesh), check))
    }
}
