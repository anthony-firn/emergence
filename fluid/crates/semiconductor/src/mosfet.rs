//! Planar bulk MOSFET structures for the drift–diffusion solver.
//!
//! Cross-section (x along the channel, y down into the silicon):
//! gate on top of the oxide between the source and drain mask edges;
//! n⁺ (p⁺ for pMOS) source/drain with Gaussian vertical and erfc lateral
//! profiles; source/drain contacts on the silicon surface at the outer ends
//! (the oxide above them is inert); body contact at the bottom.

use crate::device::{Contact, ContactKind, Device};
use crate::mesh::{graded, Material, Mesh};
use crate::silicon::{band_gap, Silicon};

#[derive(Clone, Debug)]
pub struct MosfetSpec {
    /// gate (mask) length (m)
    pub l_gate: f64,
    pub tox: f64,
    /// substrate doping (m⁻³)
    pub n_sub: f64,
    /// source/drain peak doping (m⁻³)
    pub n_sd: f64,
    /// lateral straggle of the source/drain profile (m)
    pub sigma_x: f64,
    /// source/drain junction depth (m)
    pub y_j: f64,
    /// source/drain region length outside the gate (m)
    pub l_sd: f64,
    /// silicon depth (m)
    pub depth: f64,
    /// device width (m)
    pub width: f64,
    pub pmos: bool,
    /// mesh refinement factor (1 = default spacing)
    pub refine: f64,
}

impl MosfetSpec {
    pub fn nmos(l_gate: f64) -> Self {
        MosfetSpec { l_gate, tox: 10e-9, n_sub: 1e23, n_sd: 1e26, sigma_x: 20e-9, y_j: 0.2e-6, l_sd: 1e-6, depth: 1e-6, width: 1e-6, pmos: false, refine: 1.0 }
    }
    /// A scaled technology for short-channel studies: 4 nm oxide,
    /// substrate 5e17 cm⁻³, 60 nm junctions, abrupt lateral profiles.
    pub fn scaled(l_gate: f64, pmos: bool) -> Self {
        MosfetSpec { l_gate, tox: 4e-9, n_sub: 5e23, n_sd: 1e26, sigma_x: 8e-9, y_j: 60e-9, l_sd: 0.3e-6, depth: 0.5e-6, width: 1e-6, pmos, refine: 1.0 }
    }
    /// Source and drain mask edges (x).
    pub fn edges(&self) -> (f64, f64) {
        (self.l_sd, self.l_sd + self.l_gate)
    }
    /// Lateral distance from the mask edge to the metallurgical junction
    /// at the surface.
    pub fn lateral_junction(&self) -> f64 {
        // 0.5 erfc(ξ) N_sd = N_sub
        let target = 2.0 * self.n_sub / self.n_sd;
        let (mut lo, mut hi) = (-5.0, 8.0);
        for _ in 0..200 {
            let m = 0.5 * (lo + hi);
            if erfc(m) > target {
                lo = m;
            } else {
                hi = m;
            }
        }
        0.5 * (lo + hi) * self.sigma_x
    }
    /// Metallurgical channel length at the surface.
    pub fn l_eff(&self) -> f64 {
        self.l_gate - 2.0 * self.lateral_junction()
    }
    /// Gate work-function offset (n⁺ poly for nMOS, p⁺ poly for pMOS).
    pub fn gate_offset(&self, t: f64) -> f64 {
        let o = 0.5 * band_gap(t);
        if self.pmos {
            -o
        } else {
            o
        }
    }

    pub fn build(&self, si: Silicon) -> Device {
        let s = self.clone();
        let (xs, xd) = self.edges();
        let lt = xd + self.l_sd;
        let lj = self.lateral_junction();
        let r = self.refine;
        // x: fine at the junctions, growing quickly to L/25 mid-channel
        let hj = 2e-9 / r;
        let hmid = self.l_gate / 25.0 / r;
        let hc = (100e-9f64).min(self.l_sd / 6.0) / r;
        let lch = self.l_eff();
        let tl = (0.15e-6f64).min(0.25 * lch);
        let h1 = (15e-9 / r).min(hmid);
        let outer = (0.3e-6f64).min(0.6 * self.l_sd);
        let mut bp = vec![(0.0, hc), (xs - outer, hc), (xs + lj, hj)];
        if 2.0 * tl < 0.8 * lch {
            bp.push((xs + lj + tl, h1));
            bp.push((0.5 * (xs + xd), hmid.max(h1)));
            bp.push((xd - lj - tl, h1));
        } else {
            bp.push((0.5 * (xs + xd), hmid.min(hj.max(lch / 20.0))));
        }
        bp.extend([(xd - lj, hj), (xd + outer, hc), (lt, hc)]);
        let x = graded(&bp);
        // y: oxide, then silicon fine at the interface (inversion layer)
        let hs = 0.2e-9 / r;
        let mut y = graded(&[(-self.tox, self.tox / 4.0 / r), (0.0, hs)]);
        y.pop();
        y.extend(graded(&[(0.0, hs), (5e-9, 1.5e-9 / r), (self.y_j * 1.3, self.y_j / 8.0 / r), (self.depth, self.depth / 6.0 / r)]));
        // x breakpoints may collide for short gates: keep strictly increasing
        let mut xx: Vec<f64> = Vec::with_capacity(x.len());
        for v in x {
            if xx.last().is_none_or(|&l| v > l + 1e-13) {
                xx.push(v);
            }
        }
        let mesh = Mesh::new(xx, y, |_, y| if y < 0.0 { Material::Oxide } else { Material::Silicon });
        let j0 = mesh.y.iter().position(|&v| v == 0.0).unwrap();
        let ny = mesh.ny();
        let lc = (0.5e-6f64).min(0.5 * self.l_sd);
        let mut src = vec![];
        let mut drn = vec![];
        let mut gate = vec![];
        let mut body = vec![];
        for i in 0..mesh.nx() {
            let xi = mesh.x[i];
            if xi <= lc + 1e-12 {
                src.push(mesh.node(i, j0));
            }
            if xi >= lt - lc - 1e-12 {
                drn.push(mesh.node(i, j0));
            }
            if xi >= xs - 1e-12 && xi <= xd + 1e-12 {
                gate.push(mesh.node(i, 0));
            }
            body.push(mesh.node(i, ny - 1));
        }
        let offset = self.gate_offset(si.t);
        let contacts = vec![
            Contact { name: "source".into(), nodes: src, kind: ContactKind::Ohmic, v: 0.0 },
            Contact { name: "drain".into(), nodes: drn, kind: ContactKind::Ohmic, v: 0.0 },
            Contact { name: "gate".into(), nodes: gate, kind: ContactKind::Gate { offset }, v: 0.0 },
            Contact { name: "body".into(), nodes: body, kind: ContactKind::Ohmic, v: 0.0 },
        ];
        let sigma_y = self.y_j / (self.n_sd / self.n_sub).ln().sqrt();
        let doping = move |x: f64, y: f64| -> (f64, f64) {
            let vert = (-(y / sigma_y).powi(2)).exp();
            let lat = 0.5 * erfc((x - xs) / s.sigma_x) + 0.5 * erfc((xd - x) / s.sigma_x);
            let sd = s.n_sd * vert * (lat.min(1.0));
            if s.pmos {
                (s.n_sub, sd)
            } else {
                (sd, s.n_sub)
            }
        };
        Device::new(mesh, si, doping, contacts, self.width)
    }
}

/// Complementary error function (Chebyshev fit, relative error < 1.2e-7:
/// used only for doping profiles).
pub fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let r = t * (-z * z - 1.265_512_23 + t * (1.000_023_68 + t * (0.374_091_96 + t * (0.096_784_18 + t * (-0.186_288_06 + t * (0.278_868_07 + t * (-1.135_203_98 + t * (1.488_515_87 + t * (-0.822_152_23 + t * 0.170_872_77))))))))).exp();
    if x >= 0.0 {
        r
    } else {
        2.0 - r
    }
}
