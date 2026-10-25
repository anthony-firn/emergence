//! Silicon process flows.
//!
//! A planar process is geometry in one direction and a stack in the other:
//! masks are intervals in `x`, films stack in `y`, and dopants are profiles
//! in `y` blended laterally at the mask edges. That is exactly the shape of
//! the input `semiconductor::device::Device` takes — a material function
//! and a doping function of `(x, y)` on a tensor mesh — so a [`Process`]
//! here is a drop-in replacement for a hand-written device structure, and
//! the test suite builds the validated MOSFET of `semiconductor::mosfet`
//! through this path.
//!
//! Conventions match that crate: `y = 0` is the silicon surface, `y > 0` is
//! into the substrate, and films sit at `y < 0`. Lengths are metres,
//! concentrations per cubic metre, doses per square metre.
//!
//! ```text
//!   substrate doping
//!        + grow/deposit films (each over its own mask windows)
//!        + implants (depth profile x lateral erfc blend at the windows)
//!        = material(x, y) and doping(x, y) -> semiconductor::Device
//!                                         -> Solid, LayerStack
//! ```
//!
//! # What is modelled, and what is not
//!
//! * **Implant depth profiles** are Gaussian (ion implantation), surface
//!   Gaussian (an implant driven in from the surface) or complementary
//!   error function (diffusion from a constant source). All three are the
//!   textbook forms (Sze & Ng, *Physics of Semiconductor Devices*, 3rd
//!   ed., ch. 1; Plummer, Deal & Griffin, *Silicon VLSI Technology*, 2000,
//!   ch. 7-8). Dose is conserved exactly, which the tests check against the
//!   analytic integral.
//! * **Lateral spread** at a mask edge is the standard
//!   `0.5 erfc((x - edge)/sigma_lat)` blend. Real lateral straggle is about
//!   0.7 of the vertical for a crystalline target and depends on the mask
//!   edge profile; it is a parameter here, not a prediction.
//! * **Oxidation** is the Deal-Grove law [`deal_grove`], solved exactly.
//!   **No Arrhenius constants are built in**: `A`, `B` and `tau` are the
//!   caller's, because the tabulated values differ between compilations and
//!   depend on orientation, ambient and doping. The two asymptotic limits
//!   are what the tests check.
//! * **Diffusion during anneal** broadens a profile as
//!   `sigma^2 -> sigma^2 + 2 D t` ([`Implant::annealed`]), the exact
//!   solution of Fick's equation for a Gaussian. [`arrhenius`] evaluates
//!   `D = D0 exp(-Ea/kT)` but the crate ships **no species constants**, for
//!   the same reason.
//! * **Not modelled**: segregation at the oxide interface, transient
//!   enhanced diffusion, clustering and solubility limits, channelling
//!   tails (the Gaussian misses these badly in the deep tail), stress, or
//!   anything three-dimensional. This is a structure generator, not a
//!   process simulator.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::sim::{Layer, LayerStack};
use crate::sketch::Region2;
use crate::solid::Solid;

/// Boltzmann constant, J/K (CODATA 2018, exact).
pub const K_B: f64 = 1.380649e-23;
/// Elementary charge, C (CODATA 2018, exact).
pub const Q_E: f64 = 1.602176634e-19;

/// Complementary error function.
///
/// Chebyshev fit with relative error below 1.2e-7 (Press et al.,
/// *Numerical Recipes*, 2nd ed., §6.2). Good enough for doping profiles,
/// which is all it is used for; it is the same fit `semiconductor` uses, so
/// profiles built here and there agree exactly.
pub fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let r = t
        * (-z * z - 1.265_512_23
            + t * (1.000_023_68
                + t * (0.374_091_96
                    + t * (0.096_784_18
                        + t * (-0.186_288_06
                            + t * (0.278_868_07
                                + t * (-1.135_203_98 + t * (1.488_515_87 + t * (-0.822_152_23 + t * 0.170_872_77)))))))))
            .exp();
    if x >= 0.0 {
        r
    } else {
        2.0 - r
    }
}

/// Whether a dopant donates electrons or accepts them.
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum DopantType {
    Donor,
    Acceptor,
}

/// A dopant species. Only its type affects the structure; the name is
/// carried through for reporting.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Dopant {
    pub name: String,
    pub kind: DopantType,
}

impl Dopant {
    pub fn donor(name: &str) -> Self {
        Self { name: name.into(), kind: DopantType::Donor }
    }
    pub fn acceptor(name: &str) -> Self {
        Self { name: name.into(), kind: DopantType::Acceptor }
    }
    pub fn boron() -> Self {
        Self::acceptor("boron")
    }
    pub fn phosphorus() -> Self {
        Self::donor("phosphorus")
    }
    pub fn arsenic() -> Self {
        Self::donor("arsenic")
    }
}

/// A depth profile, in concentration per cubic metre as a function of `y`.
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
#[serde(tag = "profile", rename_all = "snake_case")]
pub enum Depth {
    /// Ion implantation: `N(y) = dose / (sqrt(2 pi) straggle) *
    /// exp(-(y - range)^2 / (2 straggle^2))`. Dose is conserved.
    Implanted { dose: f64, range: f64, straggle: f64 },
    /// A surface source driven in: `N(y) = peak exp(-(y/length)^2)`.
    Surface { peak: f64, length: f64 },
    /// Diffusion from a constant surface source:
    /// `N(y) = surface erfc(y / length)`.
    Erfc { surface: f64, length: f64 },
}

impl Depth {
    /// A surface profile whose concentration falls to `background` at depth
    /// `junction`: `length = junction / sqrt(ln(peak/background))`. This is
    /// how a junction depth is specified in practice.
    pub fn for_junction(peak: f64, background: f64, junction: f64) -> Self {
        Depth::Surface { peak, length: junction / (peak / background).ln().sqrt() }
    }

    pub fn at(&self, y: f64) -> f64 {
        match *self {
            Depth::Implanted { dose, range, straggle } => {
                dose / ((2.0 * PI).sqrt() * straggle) * (-((y - range) / straggle).powi(2) / 2.0).exp()
            }
            Depth::Surface { peak, length } => peak * (-(y / length).powi(2)).exp(),
            Depth::Erfc { surface, length } => surface * erfc(y / length),
        }
    }

    /// Total dose per unit area over `y >= 0`, analytically.
    pub fn dose(&self) -> f64 {
        match *self {
            // The full Gaussian integral; the part at y < 0 is lost into
            // the film above, and is negligible when range > 3 straggle.
            Depth::Implanted { dose, range, straggle } => {
                0.5 * dose * erfc(-range / (straggle * 2.0_f64.sqrt()))
            }
            Depth::Surface { peak, length } => 0.5 * peak * length * PI.sqrt(),
            Depth::Erfc { surface, length } => surface * length / PI.sqrt(),
        }
    }

    /// Depth at which the profile crosses `level`, by bisection; `None` if
    /// it never does.
    pub fn junction_depth(&self, level: f64) -> Option<f64> {
        let peak_depth = match *self {
            Depth::Implanted { range, .. } => range,
            _ => 0.0,
        };
        if self.at(peak_depth) < level {
            return None;
        }
        let (mut lo, mut hi) = (peak_depth, peak_depth + 1e-3);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if self.at(mid) > level {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some(0.5 * (lo + hi))
    }

    /// The profile after an anneal of diffusion length `dt = D t` (m^2):
    /// a Gaussian broadens as `sigma^2 -> sigma^2 + 2 D t`, the exact
    /// solution of Fick's equation.
    pub fn annealed(&self, dt: f64) -> Depth {
        match *self {
            Depth::Implanted { dose, range, straggle } => Depth::Implanted {
                dose,
                range,
                straggle: (straggle * straggle + 2.0 * dt).sqrt(),
            },
            // A surface Gaussian N0 exp(-(y/L)^2) has sigma = L/sqrt(2),
            // and conserves dose, so the peak falls as L grows.
            Depth::Surface { peak, length } => {
                let new_length = (length * length + 4.0 * dt).sqrt();
                Depth::Surface { peak: peak * length / new_length, length: new_length }
            }
            // A constant-source profile keeps its surface value.
            Depth::Erfc { surface, length } => Depth::Erfc {
                surface,
                length: (length * length + 4.0 * dt).sqrt(),
            },
        }
    }
}

/// `D = D0 exp(-Ea / kT)`, with `d0` in m^2/s, `ea` in eV and `t` in
/// kelvin. The crate ships no species constants: supply the ones that
/// match the process.
pub fn arrhenius(d0: f64, ea_ev: f64, t: f64) -> f64 {
    d0 * (-ea_ev * Q_E / (K_B * t)).exp()
}

/// Deal-Grove oxide thickness: the root of `x^2 + A x = B (t + tau)`
/// (Deal & Grove, *J. Appl. Phys.* 36(12), 1965).
///
/// `a` (m) and `b` (m^2/s) are the linear and parabolic rate constants and
/// `tau` (s) the offset for any oxide already present. The short-time
/// limit is `x -> (B/A)(t + tau)` and the long-time limit `x -> sqrt(B t)`.
pub fn deal_grove(a: f64, b: f64, tau: f64, time: f64) -> f64 {
    0.5 * a * ((1.0 + 4.0 * b * (time + tau) / (a * a)).sqrt() - 1.0)
}

/// The `tau` that accounts for an oxide of thickness `x0` already grown.
pub fn deal_grove_offset(a: f64, b: f64, x0: f64) -> f64 {
    (x0 * x0 + a * x0) / b
}

/// A deposited or grown film, present only over its mask windows.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Film {
    pub material: String,
    pub thickness: f64,
    /// `x` intervals where the film exists; empty means everywhere.
    #[serde(default)]
    pub windows: Vec<[f64; 2]>,
    /// Relative permittivity, if known (3.9 for thermal SiO2).
    #[serde(default)]
    pub eps_r: Option<f64>,
}

impl Film {
    pub fn new(material: &str, thickness: f64) -> Self {
        Self { material: material.into(), thickness, windows: Vec::new(), eps_r: None }
    }

    pub fn over(mut self, windows: Vec<[f64; 2]>) -> Self {
        self.windows = windows;
        self
    }

    pub fn with_eps(mut self, eps_r: f64) -> Self {
        self.eps_r = Some(eps_r);
        self
    }

    pub fn covers(&self, x: f64) -> bool {
        self.windows.is_empty() || self.windows.iter().any(|w| x >= w[0] && x <= w[1])
    }
}

/// A doping step: a depth profile let in through mask windows.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Implant {
    pub dopant: Dopant,
    pub depth: Depth,
    /// Lateral straggle at a window edge.
    #[serde(default)]
    pub lateral: f64,
    /// `x` intervals open to the implant; empty means a blanket implant.
    /// Infinite bounds are allowed and are how a half-open window (a
    /// source or drain outside a gate) is written.
    #[serde(default)]
    pub windows: Vec<[f64; 2]>,
}

impl Implant {
    pub fn blanket(dopant: Dopant, depth: Depth) -> Self {
        Self { dopant, depth, lateral: 0.0, windows: Vec::new() }
    }

    pub fn through(dopant: Dopant, depth: Depth, windows: Vec<[f64; 2]>, lateral: f64) -> Self {
        Self { dopant, depth, lateral, windows }
    }

    /// Fraction of the implant reaching `x`: the sum over windows of the
    /// standard two-sided erfc blend, capped at one.
    pub fn lateral_factor(&self, x: f64) -> f64 {
        if self.windows.is_empty() {
            return 1.0;
        }
        if self.lateral <= 0.0 {
            return if self.windows.iter().any(|w| x >= w[0] && x <= w[1]) { 1.0 } else { 0.0 };
        }
        let mut total = 0.0;
        for w in &self.windows {
            let left = if w[0].is_infinite() { 1.0 } else { 0.5 * erfc((w[0] - x) / self.lateral) };
            let right = if w[1].is_infinite() { 1.0 } else { 0.5 * erfc((x - w[1]) / self.lateral) };
            total += (left + right - 1.0).max(0.0);
        }
        total.min(1.0)
    }

    pub fn at(&self, x: f64, y: f64) -> f64 {
        self.depth.at(y) * self.lateral_factor(x)
    }

    /// This implant after an anneal of diffusion length `dt = D t`.
    pub fn annealed(&self, dt: f64) -> Implant {
        Implant {
            depth: self.depth.annealed(dt),
            // Lateral diffusion broadens the edge by the same amount.
            lateral: (self.lateral * self.lateral + 4.0 * dt).sqrt(),
            ..self.clone()
        }
    }
}

/// A planar process flow.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Process {
    /// Background doping of the wafer.
    pub substrate: (Dopant, f64),
    /// Films above the silicon surface, in deposition order.
    #[serde(default)]
    pub films: Vec<Film>,
    #[serde(default)]
    pub implants: Vec<Implant>,
    /// Silicon depth modelled below the surface.
    pub depth: f64,
    /// `x` extent of the cross-section.
    pub extent: [f64; 2],
}

impl Process {
    pub fn new(substrate: Dopant, concentration: f64, extent: [f64; 2], depth: f64) -> Self {
        Self {
            substrate: (substrate, concentration),
            films: Vec::new(),
            implants: Vec::new(),
            depth,
            extent,
        }
    }

    pub fn film(mut self, film: Film) -> Self {
        self.films.push(film);
        self
    }

    /// Grow a thermal oxide of the Deal-Grove thickness for the given rate
    /// constants and time. The permittivity is thermal SiO2's 3.9 (Sze &
    /// Ng, 3rd ed., appendix G).
    pub fn oxidize(self, a: f64, b: f64, tau: f64, time: f64) -> Self {
        let t = deal_grove(a, b, tau, time);
        self.film(Film::new("oxide", t).with_eps(3.9))
    }

    pub fn implant(mut self, implant: Implant) -> Self {
        self.implants.push(implant);
        self
    }

    /// Anneal every implant placed so far by a diffusion length `dt = D t`.
    pub fn anneal(mut self, dt: f64) -> Self {
        self.implants = self.implants.iter().map(|i| i.annealed(dt)).collect();
        self
    }

    /// Total film thickness at `x` (the height of the surface above
    /// `y = 0`).
    pub fn surface_height(&self, x: f64) -> f64 {
        self.films.iter().filter(|f| f.covers(x)).map(|f| f.thickness).sum()
    }

    /// Which material is at `(x, y)`: a film name, `"silicon"`, or `None`
    /// for ambient above the stack.
    pub fn material_at(&self, x: f64, y: f64) -> Option<&str> {
        if y >= 0.0 {
            return if y <= self.depth { Some("silicon") } else { None };
        }
        // Films stack upward from the surface in deposition order.
        let mut top = 0.0;
        for f in self.films.iter().filter(|f| f.covers(x)) {
            let bottom = top;
            top += f.thickness;
            if -y > bottom && -y <= top {
                return Some(&f.material);
            }
        }
        None
    }

    /// Net donor and acceptor concentration at `(x, y)`, the pair
    /// `semiconductor::device::Device::new` takes.
    ///
    /// Implants are summed per type and the substrate added, which is the
    /// usual convention: compensation is handled by the device solver
    /// through `N_D - N_A`.
    pub fn doping_at(&self, x: f64, y: f64) -> (f64, f64) {
        let (mut donors, mut acceptors) = (0.0, 0.0);
        match self.substrate.0.kind {
            DopantType::Donor => donors += self.substrate.1,
            DopantType::Acceptor => acceptors += self.substrate.1,
        }
        for i in &self.implants {
            let n = i.at(x, y.max(0.0));
            match i.dopant.kind {
                DopantType::Donor => donors += n,
                DopantType::Acceptor => acceptors += n,
            }
        }
        (donors, acceptors)
    }

    /// A closure of the same signature the device solver wants.
    pub fn doping(&self) -> impl Fn(f64, f64) -> (f64, f64) + '_ {
        move |x, y| self.doping_at(x, y)
    }

    /// Net doping `N_D - N_A` at `(x, y)`.
    pub fn net_doping(&self, x: f64, y: f64) -> f64 {
        let (d, a) = self.doping_at(x, y);
        d - a
    }

    /// Every depth at which the net doping changes sign under `x`,
    /// shallowest first.
    ///
    /// A stacked flow has several: an emitter-base junction, a base-well
    /// junction and a well-substrate junction all appear under the same
    /// `x`. The scan uses `samples` intervals, so two junctions closer
    /// than `depth/samples` can be missed together.
    pub fn junctions(&self, x: f64, samples: usize) -> Vec<f64> {
        let samples = samples.max(2);
        let mut out = Vec::new();
        let mut y0 = 0.0;
        let mut n0 = self.net_doping(x, y0);
        for i in 1..=samples {
            let y1 = self.depth * i as f64 / samples as f64;
            let n1 = self.net_doping(x, y1);
            if n0 != 0.0 && n1 != 0.0 && n0.signum() != n1.signum() {
                let (mut lo, mut hi) = (y0, y1);
                let s = n0.signum();
                for _ in 0..100 {
                    let mid = 0.5 * (lo + hi);
                    if self.net_doping(x, mid).signum() == s {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                out.push(0.5 * (lo + hi));
            }
            y0 = y1;
            n0 = n1;
        }
        out
    }

    /// Depth of the shallowest metallurgical junction under `x`; `None`
    /// if the net doping never changes sign.
    pub fn junction_depth(&self, x: f64) -> Option<f64> {
        self.junctions(x, 2000).first().copied()
    }

    /// The cross-section as a solid, extruded over `width` in `z`: the
    /// films where they exist plus the silicon body.
    pub fn solid(&self, width: f64) -> Solid {
        let mut parts = Vec::new();
        let body = Region2::rect(
            [0.5 * (self.extent[0] + self.extent[1]), 0.5 * self.depth],
            [0.5 * (self.extent[1] - self.extent[0]), 0.5 * self.depth],
        );
        parts.push(Solid::extrude(body, 0.0, width));
        let mut top = 0.0;
        for f in &self.films {
            let bottom = top;
            top += f.thickness;
            let spans: Vec<[f64; 2]> = if f.windows.is_empty() {
                vec![self.extent]
            } else {
                f.windows
                    .iter()
                    .map(|w| [w[0].max(self.extent[0]), w[1].min(self.extent[1])])
                    .filter(|w| w[1] > w[0])
                    .collect()
            };
            for s in spans {
                let region = Region2::rect(
                    [0.5 * (s[0] + s[1]), -0.5 * (bottom + top)],
                    [0.5 * (s[1] - s[0]), 0.5 * (top - bottom)],
                );
                parts.push(Solid::extrude(region, 0.0, width));
            }
        }
        Solid::union(parts)
    }

    /// The blanket part of the stack as a layer stack, silicon first.
    /// Films that are masked (not full-area) are left out, with their
    /// names reported, since a 1D stack cannot represent them.
    pub fn layer_stack(&self) -> (LayerStack, Vec<String>) {
        let mut layers = vec![Layer::new("silicon", self.depth).electrical(11.7, 0.0)];
        let mut skipped = Vec::new();
        for f in &self.films {
            if f.windows.is_empty() {
                let mut l = Layer::new(&f.material, f.thickness);
                if let Some(e) = f.eps_r {
                    l = l.electrical(e, 0.0);
                }
                layers.push(l);
            } else {
                skipped.push(f.material.clone());
            }
        }
        layers.reverse();
        (LayerStack::new(layers), skipped)
    }
}

/// A planar MOSFET process: substrate, gate oxide, gate, and a
/// self-aligned source/drain implant.
///
/// The parameters are the ones a device structure is quoted with, and the
/// resulting [`Process`] reproduces the structure of
/// `semiconductor::mosfet::MosfetSpec` exactly, which is what the tests
/// check.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct MosProcess {
    pub gate_length: f64,
    pub oxide: f64,
    /// Substrate concentration, m^-3.
    pub substrate: f64,
    /// Source/drain peak concentration, m^-3.
    pub source_drain: f64,
    /// Lateral straggle of the source/drain implant.
    pub lateral: f64,
    /// Source/drain junction depth.
    pub junction: f64,
    /// Length of source and drain outside the gate.
    pub access: f64,
    pub depth: f64,
    /// p-channel when true: the substrate is n-type and the implant p.
    #[serde(default)]
    pub pmos: bool,
    /// Gate film thickness (polysilicon).
    #[serde(default = "gate_thickness")]
    pub gate_thickness: f64,
}

fn gate_thickness() -> f64 {
    100e-9
}

impl MosProcess {
    pub fn nmos(gate_length: f64) -> Self {
        Self {
            gate_length,
            oxide: 10e-9,
            substrate: 1e23,
            source_drain: 1e26,
            lateral: 20e-9,
            junction: 0.2e-6,
            access: 1e-6,
            depth: 1e-6,
            pmos: false,
            gate_thickness: gate_thickness(),
        }
    }

    /// Mask edges of the gate.
    pub fn edges(&self) -> (f64, f64) {
        (self.access, self.access + self.gate_length)
    }

    pub fn total_length(&self) -> f64 {
        2.0 * self.access + self.gate_length
    }

    pub fn process(&self) -> Process {
        let (xs, xd) = self.edges();
        let (substrate_dopant, sd_dopant) = if self.pmos {
            (Dopant::phosphorus(), Dopant::boron())
        } else {
            (Dopant::boron(), Dopant::arsenic())
        };
        Process::new(substrate_dopant, self.substrate, [0.0, self.total_length()], self.depth)
            .film(Film::new("oxide", self.oxide).with_eps(3.9))
            .film(Film::new("poly", self.gate_thickness).over(vec![[xs, xd]]))
            .implant(Implant::through(
                sd_dopant,
                Depth::for_junction(self.source_drain, self.substrate, self.junction),
                vec![[f64::NEG_INFINITY, xs], [xd, f64::INFINITY]],
                self.lateral,
            ))
    }

    /// Lateral distance from the mask edge to the surface junction, where
    /// the implant falls to the substrate level.
    pub fn lateral_junction(&self) -> f64 {
        let target = 2.0 * self.substrate / self.source_drain;
        let (mut lo, mut hi) = (-5.0, 8.0);
        for _ in 0..200 {
            let m = 0.5 * (lo + hi);
            if erfc(m) > target {
                lo = m;
            } else {
                hi = m;
            }
        }
        0.5 * (lo + hi) * self.lateral
    }

    /// Metallurgical channel length at the surface.
    pub fn effective_length(&self) -> f64 {
        self.gate_length - 2.0 * self.lateral_junction()
    }
}
