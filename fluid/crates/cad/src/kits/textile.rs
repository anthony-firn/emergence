//! Woven and weft-knitted fabrics.
//!
//! A fabric is a set of yarn centrelines and a cross-section carried along
//! them, so it is built with [`crate::path::sweep`] rather than with the
//! implicit algebra; the swept yarns then become a [`crate::Solid`] through
//! their exact signed distance, and compose with everything else.
//!
//! ```text
//!   yarn count (tex) + fibre density + packing -> yarn diameter
//!   weave pattern + thread spacings            -> centreline paths
//!   path + elliptical cross-section            -> swept yarn solids
//!                                              -> cover factor, crimp,
//!                                                 areal density, fibre
//!                                                 volume fraction
//! ```
//!
//! # The published relations this is checked against
//!
//! * **Yarn diameter from count**: `d = sqrt(4 T / (pi rho phi))` with `T`
//!   the linear density in kg/m, `rho` the fibre density and `phi` the
//!   packing factor (the fraction of the yarn's circle that is fibre, about
//!   0.6 for a ring-spun staple yarn). This is the standard definition
//!   (Hearle, Grosberg & Backer, *Structural Mechanics of Fibers, Yarns and
//!   Fabrics*, 1969, ch. 1).
//! * **Peirce's plain-weave geometry** (Peirce, *J. Text. Inst.* 28, 1937):
//!   the flexible-thread model whose constraint `h1 + h2 = d1 + d2` fixes
//!   the fabric thickness, with
//!   `p1 = (l2 - D theta2) cos theta2 + D sin theta2` relating thread
//!   spacing, yarn length and weave angle. [`peirce_plain`] solves it.
//! * **Cover factor** (Peirce 1937): `K = d/p` in each direction, and the
//!   fabric cover `K1 + K2 - K1 K2`.
//! * **Munden's dimensional relations for plain knitting** (Munden,
//!   *J. Text. Inst.* 50, 1959): in a relaxed fabric the loop takes a shape
//!   independent of the yarn, so courses and wales per unit length go as
//!   `c = Kc/l`, `w = Kw/l` with `l` the loop length. The published
//!   constants are `Kc = 5.0, Kw = 3.8` dry-relaxed and
//!   `Kc = 5.3, Kw = 4.1` wet-relaxed, giving a stitch density
//!   `Ks = Kc Kw` of 19.0 and 21.7.
//!
//! # Honest limits
//!
//! * The centreline is a **sinusoid through the interlacing points**, not
//!   Peirce's circular-arc-and-straight-line path. [`peirce_plain`] is
//!   provided and tested, but the geometry that gets swept is the sinusoid,
//!   which is what most textile models use and which is smooth enough to
//!   sweep. Crimp from the two differs by a few per cent at high cover.
//! * **No mechanics.** Nothing here relaxes the structure, so a jammed or
//!   over-set weave produces interpenetrating yarns.
//!   [`Weave::interference`] reports how far into that state the geometry
//!   is, and [`Weave::check`] refuses to claim a watertight solid when the
//!   yarns overlap.
//! * **Munden's constants are empirical** and were measured on relaxed
//!   wool, cotton, nylon and acrylic jersey. They are used here only to
//!   *choose* the course and wale spacing of a relaxed fabric; the geometry
//!   is then exact for those spacings. Nothing here predicts them.
//! * Yarns are smooth elliptical tubes: no fibre-level structure, no
//!   hairiness, no twist.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::path::{sweep_checked_oriented, Path3, Profile, SweepCheck};
use crate::solid::Solid;
use geometry::TriMesh;

/// Fibre densities, kg/m^3 (Morton & Hearle, *Physical Properties of
/// Textile Fibres*, 4th ed., 2008, table 1.1).
pub mod fibre_density {
    pub const COTTON: f64 = 1540.0;
    pub const WOOL: f64 = 1310.0;
    pub const POLYESTER: f64 = 1380.0;
    pub const NYLON_66: f64 = 1140.0;
    pub const E_GLASS: f64 = 2540.0;
    pub const CARBON: f64 = 1780.0;
}

/// A yarn: its count, what it is made of, and the cross-section it takes
/// in the fabric.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Yarn {
    /// Linear density in tex (g/km).
    pub tex: f64,
    /// Fibre density, kg/m^3.
    pub fibre_density: f64,
    /// Packing factor: the fibre fraction of the yarn's own section.
    #[serde(default = "packing")]
    pub packing: f64,
    /// Flattening in the fabric: the section is an ellipse of width
    /// `flattening * d` and thickness `d / flattening`, so its area is
    /// unchanged.
    #[serde(default = "one")]
    pub flattening: f64,
}

fn packing() -> f64 {
    0.6
}
fn one() -> f64 {
    1.0
}

impl Yarn {
    pub fn new(tex: f64, fibre_density: f64) -> Self {
        Self { tex, fibre_density, packing: packing(), flattening: 1.0 }
    }

    pub fn flattened(mut self, flattening: f64) -> Self {
        self.flattening = flattening;
        self
    }

    /// Linear density in kg/m.
    pub fn linear_density(&self) -> f64 {
        self.tex * 1e-6
    }

    /// Circular-equivalent diameter, m.
    pub fn diameter(&self) -> f64 {
        (4.0 * self.linear_density() / (PI * self.fibre_density * self.packing)).sqrt()
    }

    /// Width and thickness of the elliptical section in the fabric.
    pub fn section(&self) -> (f64, f64) {
        let d = self.diameter();
        (d * self.flattening, d / self.flattening)
    }

    /// Cross-sectional area of the yarn, m^2 (unchanged by flattening).
    pub fn area(&self) -> f64 {
        let (w, t) = self.section();
        0.25 * PI * w * t
    }

    /// Fibre area inside the yarn, m^2.
    pub fn fibre_area(&self) -> f64 {
        self.linear_density() / self.fibre_density
    }
}

/// A weave pattern as a repeat of warp-over-weft decisions.
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(tag = "pattern", rename_all = "snake_case")]
pub enum Pattern {
    /// 1/1 plain weave.
    Plain,
    /// `over/under` twill advancing by `step` each warp end.
    Twill { over: usize, under: usize, #[serde(default = "one_step")] step: usize },
    /// Satin of `repeat` ends with the given interlacing `step`.
    Satin { repeat: usize, step: usize },
    /// Explicit repeat: `true` where the warp passes over the weft.
    Custom { repeat: Vec<Vec<bool>> },
}

fn one_step() -> usize {
    1
}

impl Pattern {
    /// Size of the repeat in (warp ends, weft picks).
    pub fn repeat_size(&self) -> (usize, usize) {
        match self {
            Pattern::Plain => (2, 2),
            Pattern::Twill { over, under, .. } => {
                let n = over + under;
                (n, n)
            }
            Pattern::Satin { repeat, .. } => (*repeat, *repeat),
            Pattern::Custom { repeat } => (repeat.len(), repeat.first().map_or(0, |r| r.len())),
        }
    }

    /// Does warp end `i` pass over weft pick `j`?
    pub fn over(&self, i: usize, j: usize) -> bool {
        match self {
            Pattern::Plain => (i + j) % 2 == 0,
            Pattern::Twill { over, under, step } => {
                let n = over + under;
                ((j + n - (i * step) % n) % n) < *over
            }
            Pattern::Satin { repeat, step } => (j + repeat - (i * step) % repeat) % repeat == 0,
            Pattern::Custom { repeat } => {
                if repeat.is_empty() {
                    false
                } else {
                    let row = &repeat[i % repeat.len()];
                    !row.is_empty() && row[j % row.len()]
                }
            }
        }
    }

    /// Floats per repeat: the mean number of consecutive picks a warp end
    /// stays on one side. A plain weave is 1.
    pub fn mean_float(&self) -> f64 {
        let (ni, nj) = self.repeat_size();
        if ni == 0 || nj == 0 {
            return 0.0;
        }
        let mut changes = 0usize;
        for i in 0..ni {
            for j in 0..nj {
                if self.over(i, j) != self.over(i, (j + 1) % nj) {
                    changes += 1;
                }
            }
        }
        if changes == 0 {
            nj as f64
        } else {
            (ni * nj) as f64 / changes as f64
        }
    }
}

/// A woven fabric.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Weave {
    pub pattern: Pattern,
    pub warp: Yarn,
    pub weft: Yarn,
    /// Spacing between warp ends (the pitch across the fabric), m.
    pub warp_spacing: f64,
    /// Spacing between weft picks along the fabric, m.
    pub weft_spacing: f64,
    /// Repeats generated in each direction.
    #[serde(default = "two")]
    pub repeats: [usize; 2],
    /// Stations per weft spacing along a yarn path.
    #[serde(default = "stations")]
    pub stations: usize,
    /// Points around a yarn section. An `n`-gon inscribed in the section
    /// has area `(n/2) sin(2 pi/n)` times the ellipse's, so 16 sides
    /// under-fills by 2.6 %, 24 by 1.1 % and 32 by 0.6 %: a swept yarn's
    /// volume is that much below [`Yarn::area`] times its length.
    #[serde(default = "sides")]
    pub sides: usize,
}

fn two() -> [usize; 2] {
    [1, 1]
}
fn stations() -> usize {
    12
}
fn sides() -> usize {
    24
}

impl Weave {
    pub fn new(pattern: Pattern, warp: Yarn, weft: Yarn, warp_spacing: f64, weft_spacing: f64) -> Self {
        Self {
            pattern,
            warp,
            weft,
            warp_spacing,
            weft_spacing,
            repeats: [1, 1],
            stations: stations(),
            sides: sides(),
        }
    }

    pub fn with_repeats(mut self, repeats: [usize; 2]) -> Self {
        self.repeats = repeats;
        self
    }

    /// Crimp amplitude of each yarn's centreline: Peirce's constraint
    /// `h1 + h2 = d1 + d2` shared equally, so each centreline swings by
    /// `+-(d1 + d2)/4`.
    pub fn crimp_height(&self) -> f64 {
        0.5 * (self.warp.section().1 + self.weft.section().1)
    }

    /// Fabric thickness: the crimp swing plus half a yarn either side,
    /// which for the shared constraint is exactly `d1 + d2`.
    pub fn thickness(&self) -> f64 {
        self.warp.section().1 + self.weft.section().1
    }

    /// Ends and picks per metre.
    pub fn density(&self) -> (f64, f64) {
        (1.0 / self.warp_spacing, 1.0 / self.weft_spacing)
    }

    /// Cover factor in each direction, `d/p`, and the fabric cover
    /// `K1 + K2 - K1 K2` (Peirce 1937). A cover above 1 in one direction
    /// means the yarns cannot fit side by side.
    pub fn cover(&self) -> (f64, f64, f64) {
        let k1 = self.warp.section().0 / self.warp_spacing;
        let k2 = self.weft.section().0 / self.weft_spacing;
        (k1, k2, k1 + k2 - k1 * k2)
    }

    /// How much the yarns overlap, as a fraction of a yarn width: zero
    /// when the fabric is open enough for the geometry to be consistent,
    /// positive when yarns interpenetrate.
    pub fn interference(&self) -> f64 {
        let (k1, k2, _) = self.cover();
        k1.max(k2) - 1.0
    }

    /// The range of yarn flattening for which this setting is
    /// geometrically consistent, as `(minimum, maximum)`.
    ///
    /// Two conditions bound it, and they pull in opposite directions:
    ///
    /// * neighbouring yarns must fit side by side, `d f / p <= 1`, which
    ///   caps the flattening at `p/d`;
    /// * the centreline must not curve more tightly than the section
    ///   reaches. The raised-cosine path of amplitude `t/2` over a half
    ///   period `p` has a radius of curvature `p^2/(pi^2 t/2)` at its
    ///   crest, and the section reaches `t/2` inward, so
    ///   `d pi / (2 f p) <= 1` and the flattening must exceed
    ///   `d pi / (2 p)`.
    ///
    /// The interval is empty unless `p > d sqrt(pi/2) = 1.2533 d`, that
    /// is, unless the cover factor of the round yarn is below 0.8. A
    /// close-set fabric cannot be drawn with a sinusoidal centreline at
    /// all, whatever the flattening: that is a property of the path model,
    /// not of real cloth, and Peirce's arc path (whose radius is
    /// `d1 + d2`) does not share it.
    pub fn consistent_flattening(&self) -> (f64, f64) {
        let pairs = [
            (self.warp.diameter(), self.warp_spacing),
            (self.weft.diameter(), self.weft_spacing),
        ];
        let lo = pairs
            .iter()
            .map(|&(d, p)| d * PI / (2.0 * p))
            .fold(0.0, f64::max);
        let hi = pairs.iter().map(|&(d, p)| p / d).fold(f64::INFINITY, f64::min);
        (lo, hi)
    }

    /// Centreline of warp end `i` over `picks` weft crossings.
    pub fn warp_path(&self, i: usize, picks: usize) -> Path3 {
        let amp = 0.5 * self.crimp_height();
        let y = i as f64 * self.warp_spacing;
        let p = self.weft_spacing;
        let n = (self.stations * picks).max(2);
        let length = picks as f64 * p;
        Path3::from_fn(n, false, |t| {
            let x = t * length;
            // The sign at each crossing comes from the pattern; between
            // crossings the path is a sinusoid through those points.
            let z = amp * self.interlace_wave(i, x / p, true);
            [x, y, z]
        })
    }

    /// Centreline of weft pick `j` over `ends` warp crossings.
    pub fn weft_path(&self, j: usize, ends: usize) -> Path3 {
        let amp = 0.5 * self.crimp_height();
        let x = j as f64 * self.weft_spacing;
        let p = self.warp_spacing;
        let n = (self.stations * ends).max(2);
        let length = ends as f64 * p;
        Path3::from_fn(n, false, |t| {
            let y = t * length;
            let z = -amp * self.interlace_wave(j, y / p, false);
            [x, y, z]
        })
    }

    /// A smooth wave that is `+1` at crossings where this yarn is on top
    /// and `-1` where it is underneath, interpolated by a raised cosine.
    fn interlace_wave(&self, index: usize, u: f64, warp: bool) -> f64 {
        let k = u.floor();
        let frac = u - k;
        let side = |m: i64| -> f64 {
            let m = m.rem_euclid(1 << 20) as usize;
            let over = if warp { self.pattern.over(index, m) } else { !self.pattern.over(m, index) };
            if over {
                1.0
            } else {
                -1.0
            }
        };
        let (a, b) = (side(k as i64), side(k as i64 + 1));
        // Raised cosine from a to b: value a at frac = 0, b at frac = 1,
        // with zero slope at both ends, so the path is smooth.
        0.5 * (a + b) + 0.5 * (a - b) * (PI * frac).cos()
    }

    /// Every yarn in the generated block, as a swept mesh with its own
    /// validity check.
    pub fn yarns(&self) -> Vec<(TriMesh, SweepCheck)> {
        let (ri, rj) = self.pattern.repeat_size();
        let ends = (ri * self.repeats[0]).max(1);
        let picks = (rj * self.repeats[1]).max(1);
        let (ww, wt) = self.warp.section();
        let (fw, ft) = self.weft.section();
        // The profile's u axis is swept along +z, so its first half-extent
        // is the yarn's thickness and its second its width. Stating the up
        // direction is what keeps a flattened yarn lying flat.
        let up = [0.0, 0.0, 1.0];
        let warp_profile = Profile::ellipse(0.5 * wt, 0.5 * ww, self.sides);
        let weft_profile = Profile::ellipse(0.5 * ft, 0.5 * fw, self.sides);
        let mut out = Vec::new();
        for i in 0..ends {
            out.push(sweep_checked_oriented(&self.warp_path(i, picks), &warp_profile, up));
        }
        for j in 0..picks {
            out.push(sweep_checked_oriented(&self.weft_path(j, ends), &weft_profile, up));
        }
        out
    }

    /// The fabric as a solid: the union of the yarns.
    pub fn solid(&self) -> Solid {
        Solid::union(
            self.yarns()
                .into_iter()
                .map(|(mesh, _)| Solid::from_trimesh(mesh))
                .collect(),
        )
    }

    /// Total centreline length of all generated yarns.
    pub fn yarn_length(&self) -> f64 {
        let (ri, rj) = self.pattern.repeat_size();
        let ends = (ri * self.repeats[0]).max(1);
        let picks = (rj * self.repeats[1]).max(1);
        (0..ends).map(|i| self.warp_path(i, picks).length()).sum::<f64>()
            + (0..picks).map(|j| self.weft_path(j, ends).length()).sum::<f64>()
    }

    /// Crimp of each set: `path length / straight length - 1`.
    pub fn crimp(&self) -> (f64, f64) {
        let (ri, rj) = self.pattern.repeat_size();
        let ends = (ri * self.repeats[0]).max(1);
        let picks = (rj * self.repeats[1]).max(1);
        let warp = (0..ends)
            .map(|i| self.warp_path(i, picks).length() / (picks as f64 * self.weft_spacing) - 1.0)
            .sum::<f64>()
            / ends as f64;
        let weft = (0..picks)
            .map(|j| self.weft_path(j, ends).length() / (ends as f64 * self.warp_spacing) - 1.0)
            .sum::<f64>()
            / picks as f64;
        (warp, weft)
    }

    /// Fabric mass per unit area, kg/m^2, from the yarn counts, the thread
    /// densities and the crimp: the standard textile calculation
    /// `(ends/m * tex_warp * (1 + c1) + picks/m * tex_weft * (1 + c2)) / 1e6`.
    pub fn areal_density(&self) -> f64 {
        let (ends, picks) = self.density();
        let (c1, c2) = self.crimp();
        ends * self.warp.linear_density() * (1.0 + c1) + picks * self.weft.linear_density() * (1.0 + c2)
    }

    /// Fibre volume fraction of the fabric: fibre area times yarn length
    /// per unit cell volume.
    pub fn fibre_volume_fraction(&self) -> f64 {
        let (ri, rj) = self.pattern.repeat_size();
        let ends = (ri * self.repeats[0]).max(1);
        let picks = (rj * self.repeats[1]).max(1);
        let cell = (ends as f64 * self.warp_spacing) * (picks as f64 * self.weft_spacing) * self.thickness();
        let (c1, c2) = self.crimp();
        let warp_len = ends as f64 * picks as f64 * self.weft_spacing * (1.0 + c1);
        let weft_len = picks as f64 * ends as f64 * self.warp_spacing * (1.0 + c2);
        (self.warp.fibre_area() * warp_len + self.weft.fibre_area() * weft_len) / cell
    }

    /// Whether the generated yarns are geometrically consistent: no
    /// interference between neighbours, and every swept yarn watertight
    /// and free of self-intersection.
    pub fn check(&self) -> Result<(), String> {
        if self.interference() > 0.0 {
            return Err(format!(
                "yarns interfere: cover factors are {:?}, so neighbouring yarns overlap by {:.1} % of a width",
                self.cover(),
                100.0 * self.interference()
            ));
        }
        for (n, (_, check)) in self.yarns().iter().enumerate() {
            if !check.watertight {
                return Err(format!("yarn {n} did not sweep watertight"));
            }
            if check.curvature_ratio >= 1.0 {
                return Err(format!(
                    "yarn {n} self-intersects at station {}: the section reaches {:.2} of the way \
                     past the centre of curvature. Flatten the yarn or open the fabric out.",
                    check.worst_station, check.curvature_ratio
                ));
            }
        }
        Ok(())
    }
}

/// Peirce's plain-weave geometry (Peirce, *J. Text. Inst.* 28, 1937).
///
/// Given the yarn diameters and the two thread spacings, solve his
/// equations for the weave angles and crimp heights:
///
/// ```text
///   p1 = (l2 - D th2) cos th2 + D sin th2
///   h2 = (l2 - D th2) sin th2 + D (1 - cos th2)
///   h1 + h2 = D = d1 + d2
/// ```
///
/// Here the unknowns are the two angles; the yarn lengths follow. Returns
/// `(theta1, theta2, h1, h2, l1, l2)` with lengths per repeat.
pub fn peirce_plain(d1: f64, d2: f64, p1: f64, p2: f64) -> Result<(f64, f64, f64, f64, f64, f64), String> {
    let dd = d1 + d2;
    if p1 <= 0.0 || p2 <= 0.0 || dd <= 0.0 {
        return Err("spacings and diameters must be positive".into());
    }
    // Peirce's relations, with h1 + h2 = D, reduce to one equation per
    // direction once h is known; solve the pair by iterating on h1.
    // Given h1, theta1 and l1 follow from the first two relations with
    // (p2, h1); then h2 = D - h1 gives theta2 and l2 from (p1, h2).
    let solve_dir = |p: f64, h: f64| -> Option<(f64, f64)> {
        // Unknown theta in (0, pi/2): eliminate l from
        //   p = (l - D th) cos th + D sin th
        //   h = (l - D th) sin th + D (1 - cos th)
        // so (l - D th) = (h - D(1 - cos th)) / sin th and
        //   p sin th = (h - D(1 - cos th)) cos th + D sin^2 th.
        let f = |th: f64| {
            let (s, c) = th.sin_cos();
            p * s - ((h - dd * (1.0 - c)) * c + dd * s * s)
        };
        let (mut lo, mut hi) = (1e-9, 0.5 * PI - 1e-9);
        if f(lo).signum() == f(hi).signum() {
            return None;
        }
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if f(mid).signum() == f(lo).signum() {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let th = 0.5 * (lo + hi);
        let (s, c) = th.sin_cos();
        let l = (h - dd * (1.0 - c)) / s + dd * th;
        Some((th, l))
    };
    // h1 is set by requiring both directions to be consistent with the
    // same D; bisect on the split.
    let residual = |h1: f64| -> Option<f64> {
        let (th1, _) = solve_dir(p2, h1)?;
        let (th2, _) = solve_dir(p1, dd - h1)?;
        // Peirce's symmetric solution balances the two weave angles
        // against the spacings: p1 sin th1 = p2 sin th2 at the solution of
        // the simultaneous pair.
        Some(p1 * th1.sin() - p2 * th2.sin())
    };
    let (mut lo, mut hi) = (1e-6 * dd, dd * (1.0 - 1e-6));
    let (mut flo, fhi) = (residual(lo), residual(hi));
    if flo.is_none() || fhi.is_none() {
        return Err("Peirce's equations have no solution for these spacings (the weave is jammed)".into());
    }
    if flo.unwrap().signum() == fhi.unwrap().signum() {
        return Err("Peirce's equations have no bracketed solution for these spacings".into());
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        let fm = residual(mid);
        match fm {
            None => return Err("Peirce iteration left the valid range".into()),
            Some(v) => {
                if v.signum() == flo.unwrap().signum() {
                    lo = mid;
                    flo = Some(v);
                } else {
                    hi = mid;
                }
            }
        }
    }
    let _ = fhi;
    let h1 = 0.5 * (lo + hi);
    let h2 = dd - h1;
    let (th1, l1) = solve_dir(p2, h1).ok_or("direction 1 failed")?;
    let (th2, l2) = solve_dir(p1, h2).ok_or("direction 2 failed")?;
    Ok((th1, th2, h1, h2, l1, l2))
}

/// How relaxed a knitted fabric is, which is what sets Munden's constants.
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Relaxation {
    /// Dry-relaxed: `Kc = 5.0`, `Kw = 3.8`.
    Dry,
    /// Wet-relaxed: `Kc = 5.3`, `Kw = 4.1`.
    Wet,
}

impl Relaxation {
    /// Munden's `(Kc, Kw)` (Munden, *J. Text. Inst.* 50, 1959).
    pub fn constants(&self) -> (f64, f64) {
        match self {
            Relaxation::Dry => (5.0, 3.8),
            Relaxation::Wet => (5.3, 4.1),
        }
    }

    /// Stitch density constant `Ks = Kc Kw`.
    pub fn stitch_constant(&self) -> f64 {
        let (kc, kw) = self.constants();
        kc * kw
    }
}

/// A plain weft-knitted (jersey) fabric.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Knit {
    pub yarn: Yarn,
    /// Loop length, m: the yarn consumed by one stitch.
    pub loop_length: f64,
    /// Course spacing (between rows), m.
    pub course_spacing: f64,
    /// Wale spacing (between columns), m.
    pub wale_spacing: f64,
    #[serde(default = "two")]
    pub repeats: [usize; 2],
    #[serde(default = "knit_stations")]
    pub stations: usize,
    #[serde(default = "sides")]
    pub sides: usize,
}

fn knit_stations() -> usize {
    48
}

impl Knit {
    /// A relaxed fabric of the given loop length, with the course and wale
    /// spacing taken from Munden's constants. The geometry is then exact
    /// for those spacings; Munden's law only chooses them.
    pub fn relaxed(yarn: Yarn, loop_length: f64, state: Relaxation) -> Self {
        let (kc, kw) = state.constants();
        Self {
            yarn,
            loop_length,
            course_spacing: loop_length / kc,
            wale_spacing: loop_length / kw,
            repeats: [1, 1],
            stations: knit_stations(),
            sides: sides(),
        }
    }

    pub fn with_repeats(mut self, repeats: [usize; 2]) -> Self {
        self.repeats = repeats;
        self
    }

    /// Stitch density, loops per square metre.
    pub fn stitch_density(&self) -> f64 {
        1.0 / (self.course_spacing * self.wale_spacing)
    }

    /// Mass per unit area, kg/m^2: loop length times stitch density times
    /// linear density.
    pub fn areal_density(&self) -> f64 {
        self.loop_length * self.stitch_density() * self.yarn.linear_density()
    }

    /// Head and sinker arc radius: a quarter of the wale spacing, so the
    /// two semicircles of a stitch span one wale.
    pub fn arc_radius(&self) -> f64 {
        0.25 * self.wale_spacing
    }

    /// Length of each straight leg.
    ///
    /// In the fabric plane a stitch is two legs plus a head semicircle and
    /// a sinker semicircle, so its planar length is `2 h + pi w / 2`. The
    /// yarn also crosses from one face of the fabric to the other, and
    /// that excursion adds length, so `h` is solved against the **three
    /// dimensional** length: the generated stitch consumes exactly
    /// `loop_length`, which is the quantity a knitter sets. The correction
    /// is about half a per cent for a normal yarn, and grows as the yarn
    /// thickens relative to the loop.
    pub fn leg_length(&self) -> Result<f64, String> {
        let arcs = PI * self.wale_spacing / 2.0;
        if self.loop_length <= arcs {
            return Err(format!(
                "loop length {:.3} mm is shorter than the head and sinker arcs alone \
                 ({:.3} mm at this wale spacing): no loop of this shape exists",
                self.loop_length * 1e3,
                arcs * 1e3
            ));
        }
        let planar = 0.5 * (self.loop_length - arcs);
        if self.spatial_length(planar) <= self.loop_length {
            return Ok(planar);
        }
        let (mut lo, mut hi) = (0.0, planar);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if self.spatial_length(mid) < self.loop_length {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Ok(0.5 * (lo + hi))
    }

    /// Planar length of a stitch with legs of length `h`.
    pub fn planar_length(&self, h: f64) -> f64 {
        2.0 * h + PI * self.wale_spacing / 2.0
    }

    /// The true three-dimensional length of a stitch with legs `h`, by
    /// quadrature of the generated centreline.
    pub fn spatial_length(&self, h: f64) -> f64 {
        let planar = self.planar_length(h);
        let n = 4000;
        let mut total = 0.0;
        let mut prev = self.stitch_point(0.0, h);
        for i in 1..=n {
            let p = self.stitch_point(planar * i as f64 / n as f64, h);
            total += geometry::norm(geometry::sub(p, prev));
            prev = p;
        }
        total
    }

    /// How far the loop reaches in the course direction, as a multiple of
    /// the course spacing. A plain-knit loop interlooks with the course
    /// below, so this is above one.
    pub fn reach_in_courses(&self) -> f64 {
        let h = self.leg_length().unwrap_or(0.0);
        (h + self.arc_radius()) / self.course_spacing
    }

    /// A point on one stitch at arc length `s` along it.
    ///
    /// The stitch is two straight legs joined by a head semicircle and
    /// closed by a sinker semicircle, which is the standard plain-knit
    /// loop shape (Leaf & Glaskin, *J. Text. Inst.* 46, 1955). It is
    /// parameterised by arc length, so a stitch consumes exactly
    /// `loop_length` of yarn by construction. The head faces one side of
    /// the fabric and the sinker loop the other, which is the `z`
    /// excursion.
    fn stitch_point(&self, s: f64, h: f64) -> [f64; 3] {
        let a = self.arc_radius();
        let r = 0.5 * self.yarn.diameter();
        let (s1, arc) = (h, PI * a);
        let s2 = s1 + arc;
        let s3 = s2 + h;
        let total = s3 + arc;
        let head_centre = s1 + 0.5 * arc;
        let z = r * (2.0 * PI * (s - head_centre) / total).cos();
        if s <= s1 {
            [0.0, s, z]
        } else if s <= s2 {
            let f = (s - s1) / arc;
            let th = PI * (1.0 - f);
            [a + a * th.cos(), h + a * th.sin(), z]
        } else if s <= s3 {
            [2.0 * a, h - (s - s2), z]
        } else {
            let f = ((s - s3) / arc).min(1.0);
            let th = PI + PI * f;
            [3.0 * a + a * th.cos(), a * th.sin(), z]
        }
    }

    /// The centreline of one course, spanning `wales` stitches.
    pub fn course_path(&self, course: usize, wales: usize) -> Path3 {
        let h = self.leg_length().unwrap_or(0.0);
        let w = self.wale_spacing;
        let y0 = course as f64 * self.course_spacing;
        let wales = wales.max(1);
        let n = (self.stations * wales).max(8);
        let total = self.loop_length;
        Path3::from_fn(n, false, move |t| {
            let s_all = t * wales as f64 * total;
            let k = (s_all / total).floor().min(wales as f64 - 1.0);
            let p = self.stitch_point(s_all - k * total, h);
            [p[0] + k * w, y0 + p[1], p[2]]
        })
    }

    /// The generated courses as swept meshes with their validity checks.
    pub fn loops(&self) -> Vec<(TriMesh, SweepCheck)> {
        let profile = Profile::circle(0.5 * self.yarn.diameter(), self.sides);
        (0..self.repeats[1].max(1))
            .map(|course| {
                sweep_checked_oriented(
                    &self.course_path(course, self.repeats[0].max(1)),
                    &profile,
                    [0.0, 0.0, 1.0],
                )
            })
            .collect()
    }

    pub fn solid(&self) -> Solid {
        Solid::union(self.loops().into_iter().map(|(m, _)| Solid::from_trimesh(m)).collect())
    }

    /// The measured yarn length per stitch of the generated geometry,
    /// which must equal `loop_length` up to the polyline's chord error.
    pub fn measured_loop_length(&self) -> f64 {
        let wales = 4;
        self.course_path(0, wales).length() / wales as f64
    }

    /// Whether the generated loops sweep without self-intersecting.
    pub fn check(&self) -> Result<(), String> {
        self.leg_length()?;
        for (n, (_, check)) in self.loops().iter().enumerate() {
            if check.curvature_ratio >= 1.0 {
                return Err(format!(
                    "loop {n} self-intersects at station {}: the yarn reaches {:.2} of the way past \
                     the centre of curvature; the loop is too tight for this yarn",
                    check.worst_station, check.curvature_ratio
                ));
            }
        }
        Ok(())
    }
}
