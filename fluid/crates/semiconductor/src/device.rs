//! Drift–diffusion device solver (the TCAD rung of the electronics ladder).
//!
//!   −∇·(ε∇ψ) = q(p − n + N_d − N_a)
//!   ∂n/∂t = ∇·J_n/q − R + G,  J_n = −qμ_n n ∇φ_n
//!   ∂p/∂t = −∇·J_p/q − R + G, J_p = −qμ_p p ∇φ_p
//! with n = n_i e^{ψ + θ_n − φ_n}, p = n_i e^{φ_p − ψ + θ_p} (potentials in
//! units of V_t). The effective potentials θ carry the physics beyond
//! Boltzmann statistics: band-gap narrowing (ΔE_g/2 per band) and Fermi–
//! Dirac degeneracy, θ = ln(F_{1/2}(η)/e^η), evaluated only where η > −12
//! (below, the Boltzmann fold is exact to 2e-6).
//!
//! Discretization: box method on a tensor mesh, materials per element (the
//! oxide carries only the Poisson equation; Si/SiO₂ interface nodes have
//! partial semiconductor boxes). Currents along edges by Scharfetter–Gummel
//! on the effective potentials ψ + θ_n and ψ − θ_p: exact for constant
//! current and effective field along the edge. Unknowns ψ, φ_n, φ_p (units
//! of V_t; densities stay positive by construction), coupled Newton with a
//! blocked banded LU (`band.rs`), exact derivatives of θ, damped steps and
//! bias continuation.
//!
//! Physics switches, each a fold with an indicator:
//! - statistics: Boltzmann / Fermi–Dirac (per node, by η);
//! - band-gap narrowing (Slotboom–de Graaff);
//! - inversion-layer mobility (Lombardi: phonon and surface-roughness
//!   scattering in the perpendicular field), velocity saturation;
//! - impact ionization (van Overstraeten–de Man), evaluated only where the
//!   field exceeds 5e4 V/cm;
//! - transients (BDF1/BDF2): carrier storage and displacement current, with
//!   total current conserved exactly by the discrete equations.
//!
//! Contacts: ohmic (charge-neutral equilibrium values at the applied
//! voltage), gates (on oxide nodes: ψ = V_g + work-function offset), and
//! fixed values (manufactured-solution tests). Free boundaries are
//! insulating (zero normal field and current).

use crate::band::Band;
use crate::mesh::{Material, Mesh};
use crate::silicon::{alpha_n, alpha_p, band_gap, bgn_slotboom, degeneracy, nc, ni_measured, nv, vsat_n, vsat_p, Silicon, EPS_OX, EPS_SI, LOMBARDI_N, LOMBARDI_P, Q};

#[derive(Clone, Debug)]
pub enum ContactKind {
    Ohmic,
    /// gate on oxide: Fermi level of the gate relative to the intrinsic
    /// level of the silicon (V): n⁺ poly ≈ +E_g/2, p⁺ poly ≈ −E_g/2
    Gate { offset: f64 },
    /// prescribed (ψ, φ_n, φ_p) per node, in units of V_t
    Fixed(Vec<[f64; 3]>),
}

#[derive(Clone, Debug)]
pub struct Contact {
    pub name: String,
    pub nodes: Vec<usize>,
    pub kind: ContactKind,
    pub v: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Statistics {
    Boltzmann,
    FermiDirac,
}

#[derive(Clone, Copy, Debug)]
struct Edge {
    a: usize,
    b: usize,
    h: f64,
    /// Σ ε·(half element width) across the edge (Poisson coupling)
    eps_w: f64,
    /// semiconductor part of the coupling width (carrier fluxes)
    si_w: f64,
    /// low-field (doping-dependent) mobilities
    mu_n: f64,
    mu_p: f64,
    /// along x (true) or y
    along_x: bool,
}

/// One stored time level: (n̂, p̂, ψ) per node.
type Level = (Vec<f64>, Vec<f64>, Vec<f64>);

/// Time-integration state: BDF coefficients and history levels.
#[derive(Clone, Debug)]
struct Transient {
    dt: f64,
    /// ∂x/∂t ≈ (a₀ x + a₁ x₁ + a₂ x₂)/dt
    a: [f64; 3],
    /// previous levels (n̂, p̂, ψ) — [0] is t_n, [1] is t_{n−1}
    hist: Vec<Level>,
    /// full state (ψ, φ_n, φ_p) of the last committed level
    committed: Level,
    /// length of the last committed step (variable-step BDF2)
    dt_prev: f64,
    /// after `accept_step`: the history and coefficients that step used, so
    /// terminal currents of the committed step can still be evaluated
    used: Option<(Vec<Level>, [f64; 3], f64)>,
}

#[derive(Clone, Debug)]
pub struct Device {
    pub mesh: Mesh,
    pub si: Silicon,
    pub nd: Vec<f64>,
    pub na: Vec<f64>,
    pub contacts: Vec<Contact>,
    /// device extent in z (m); 1D meshes: area (m²)
    pub depth: f64,
    pub psi: Vec<f64>,
    pub phin: Vec<f64>,
    pub phip: Vec<f64>,
    /// constant mobilities (m²/(V s)) instead of the doping-dependent model
    pub mu_const: Option<(f64, f64)>,
    /// field-dependent mobility (Caughey–Thomas velocity saturation, driven
    /// by the electric field along each edge)
    pub velocity_saturation: bool,
    /// include ∂μ/∂ψ in the Newton matrix. Off by default: the exact terms
    /// (negative differential mobility) spoil the M-matrix-like structure
    /// and Newton diverges from equilibrium on MOSFETs, while re-evaluating
    /// the mobility at every iterate converges in 5–7 steps to the same
    /// solution (the residual always uses the current mobility)
    pub vsat_exact_jacobian: bool,
    /// Lombardi inversion-layer mobility (perpendicular-field degradation)
    pub lombardi: bool,
    /// re-screen impurity scattering by the local carrier density
    /// (Brooks–Herring; Masetti is its neutral-bulk fold)
    pub screening: bool,
    pub statistics: Statistics,
    /// band-gap narrowing (Slotboom–de Graaff)
    pub bgn: bool,
    /// quantum confinement by density gradient along y (mass in m₀): the
    /// quantum potential Λ (units of V_t) adds to θ_n; see
    /// `update_quantum`/`solve_quantum`
    pub dg_mass: Option<f64>,
    /// quantum potential per node (units of V_t)
    pub lambda: Vec<f64>,
    /// impact ionization (van Overstraeten–de Man). The generation is
    /// re-evaluated at every Newton iterate but not linearized: its exact
    /// derivatives with respect to the majority-carrier fluxes near
    /// depletion edges dwarf the minority-carrier rows they enter, which
    /// then nearly duplicate the majority equations (a singular matrix).
    /// The fixed point is the same; convergence is linear with rate
    /// (M − 1)/M, M the multiplication factor.
    pub impact_ionization: bool,
    /// recombination switch (SRH + Auger)
    pub recombination: bool,
    /// per-node recombination multiplier (None: 1 everywhere)
    pub recomb_scale: Option<Vec<f64>>,
    /// extra residual per node and equation (manufactured solutions)
    pub source: Option<Vec<[f64; 3]>>,
    pub newton_tol: f64,
    pub max_step: f64,
    edges: Vec<Edge>,
    vol_si: Vec<f64>,
    semi: Vec<bool>,
    contact_of: Vec<Option<usize>>,
    /// band-gap narrowing per node, units of V_t
    dg: Vec<f64>,
    /// charge-neutral ψ at φ = 0 per node (statistics-consistent)
    psi0: Vec<f64>,
    tr: Option<Transient>,
    /// nodes where Fermi–Dirac was evaluated in the last assembly
    pub fd_nodes: std::cell::Cell<usize>,
    /// nodes where impact ionization was evaluated in the last assembly
    pub ii_nodes: std::cell::Cell<usize>,
}

/// Bernoulli function B(x) = x/(eˣ − 1) and its derivative.
#[inline]
pub fn bernoulli(x: f64) -> (f64, f64) {
    if x.abs() < 1e-3 {
        let x2 = x * x;
        (1.0 - 0.5 * x + x2 / 12.0 - x2 * x2 / 720.0, -0.5 + x / 6.0 - x2 * x / 180.0)
    } else {
        let b = if x > 0.0 { x * (-x).exp() / -(-x).exp_m1() } else { x / x.exp_m1() };
        (b, b * (1.0 - b - x) / x)
    }
}

/// Carrier state at the nodes for one assembly: densities (units of n_i),
/// effective potentials, and dθ/dη.
struct Carriers {
    n: Vec<f64>,
    p: Vec<f64>,
    /// ψ + θ_n and ψ − θ_p
    sn: Vec<f64>,
    tp: Vec<f64>,
    dthn: Vec<f64>,
    dthp: Vec<f64>,
    /// e^{θ_n + θ_p}: equilibrium n̂p̂
    eq: Vec<f64>,
}

/// One edge's fluxes (units of q n_i V_t per depth) and their partials
/// with respect to (ψ_a, ψ_b, φ_a, φ_b) of the carrier's own quasi-Fermi
/// potential.
#[derive(Clone, Copy, Default)]
struct Flux {
    g: f64,
    d: [f64; 4],
    /// magnitude of the separate drift/diffusion terms (round-off scale)
    s: f64,
}

impl Device {
    pub fn new(mesh: Mesh, si: Silicon, doping: impl Fn(f64, f64) -> (f64, f64), contacts: Vec<Contact>, depth: f64) -> Self {
        let n = mesh.nodes();
        let (nx, ny) = (mesh.nx(), mesh.ny());
        let mut nd = vec![0.0; n];
        let mut na = vec![0.0; n];
        let mut vol_si = vec![0.0; n];
        let mut semi = vec![false; n];
        for i in 0..nx {
            for j in 0..ny {
                let k = mesh.node(i, j);
                let (d, a) = doping(mesh.x[i], mesh.y[j]);
                nd[k] = d;
                na[k] = a;
                for (ei, ej) in [(i.wrapping_sub(1), j.wrapping_sub(1)), (i.wrapping_sub(1), j), (i, j.wrapping_sub(1)), (i, j)] {
                    if ei < nx - 1 && ej < ny - 1 && mesh.material(ei, ej) == Material::Silicon {
                        vol_si[k] += 0.25 * (mesh.x[ei + 1] - mesh.x[ei]) * (mesh.y[ej + 1] - mesh.y[ej]);
                        semi[k] = true;
                    }
                }
            }
        }
        let eps = |m: Material| if m == Material::Silicon { EPS_SI } else { EPS_OX };
        let mut edges = Vec::new();
        for i in 0..nx {
            for j in 0..ny {
                if i + 1 < nx {
                    let (mut ew, mut sw) = (0.0, 0.0);
                    for ej in [j.wrapping_sub(1), j] {
                        if ej < ny - 1 {
                            let hy = 0.5 * (mesh.y[ej + 1] - mesh.y[ej]);
                            let m = mesh.material(i, ej);
                            ew += eps(m) * hy;
                            if m == Material::Silicon {
                                sw += hy;
                            }
                        }
                    }
                    edges.push(Edge { a: mesh.node(i, j), b: mesh.node(i + 1, j), h: mesh.x[i + 1] - mesh.x[i], eps_w: ew, si_w: sw, mu_n: 0.0, mu_p: 0.0, along_x: true });
                }
                if j + 1 < ny {
                    let (mut ew, mut sw) = (0.0, 0.0);
                    for ei in [i.wrapping_sub(1), i] {
                        if ei < nx - 1 {
                            let hx = 0.5 * (mesh.x[ei + 1] - mesh.x[ei]);
                            let m = mesh.material(ei, j);
                            ew += eps(m) * hx;
                            if m == Material::Silicon {
                                sw += hx;
                            }
                        }
                    }
                    edges.push(Edge { a: mesh.node(i, j), b: mesh.node(i, j + 1), h: mesh.y[j + 1] - mesh.y[j], eps_w: ew, si_w: sw, mu_n: 0.0, mu_p: 0.0, along_x: false });
                }
            }
        }
        let mut contact_of = vec![None; n];
        for (c, ct) in contacts.iter().enumerate() {
            for &k in &ct.nodes {
                contact_of[k] = Some(c);
            }
        }
        let mut d = Device {
            mesh,
            si,
            nd,
            na,
            contacts,
            depth,
            psi: vec![0.0; n],
            phin: vec![0.0; n],
            phip: vec![0.0; n],
            mu_const: None,
            velocity_saturation: false,
            vsat_exact_jacobian: false,
            lombardi: false,
            screening: false,
            dg_mass: None,
            lambda: vec![0.0; n],
            statistics: Statistics::Boltzmann,
            bgn: false,
            impact_ionization: false,
            recombination: true,
            recomb_scale: None,
            source: None,
            newton_tol: 1e-10,
            max_step: 5.0,
            edges,
            vol_si,
            semi,
            contact_of,
            dg: vec![0.0; n],
            psi0: vec![0.0; n],
            tr: None,
            fd_nodes: std::cell::Cell::new(0),
            ii_nodes: std::cell::Cell::new(0),
        };
        d.update_mobility();
        d.update_statistics();
        d.initial_guess();
        d
    }

    /// Recompute edge low-field mobilities (after changing `mu_const` or
    /// temperature).
    pub fn update_mobility(&mut self) {
        let t = self.si.t;
        for e in self.edges.iter_mut() {
            let ntot = 0.5 * (self.nd[e.a] + self.na[e.a] + self.nd[e.b] + self.na[e.b]);
            let (mn, mp) = self.mu_const.unwrap_or((self.si.mob_n.mobility(ntot, t), self.si.mob_p.mobility(ntot, t)));
            e.mu_n = mn;
            e.mu_p = mp;
        }
    }

    /// Recompute band-gap narrowing and the neutral contact potentials
    /// (after changing `statistics`, `bgn` or temperature).
    pub fn update_statistics(&mut self) {
        let vt = self.si.vt();
        for k in 0..self.psi.len() {
            self.dg[k] = if self.bgn { bgn_slotboom(self.nd[k] + self.na[k]) / vt } else { 0.0 };
        }
        for k in 0..self.psi.len() {
            self.psi0[k] = if self.semi[k] { self.solve_neutral(k) } else { 0.0 };
        }
    }

    /// Set the lattice temperature: intrinsic density, mobilities, band
    /// gap (poly gates follow ±E_g/2), statistics. The state is rescaled to
    /// the new thermal voltage.
    pub fn set_temperature(&mut self, t: f64) {
        let vt_old = self.si.vt();
        self.si.t = t;
        self.si.ni = ni_measured(t);
        let vt = self.si.vt();
        for c in self.contacts.iter_mut() {
            if let ContactKind::Gate { offset } = &mut c.kind {
                *offset = offset.signum() * 0.5 * band_gap(t);
            }
        }
        let r = vt_old / vt;
        for k in 0..self.psi.len() {
            self.psi[k] *= r;
            self.phin[k] *= r;
            self.phip[k] *= r;
        }
        self.update_mobility();
        self.update_statistics();
        self.apply_dirichlet();
    }

    /// ln(N_c/n_i), ln(N_v/n_i): band edges relative to the intrinsic level.
    fn chi(&self) -> (f64, f64) {
        let t = self.si.t;
        ((nc(t) / self.si.ni).ln(), (nv(t) / self.si.ni).ln())
    }

    /// θ_n, dθ_n/dη at node k for ψ − φ_n = u.
    #[inline]
    fn theta_n(&self, k: usize, u: f64, chi_c: f64) -> (f64, f64) {
        let half = 0.5 * self.dg[k] + self.lambda[k];
        if self.statistics == Statistics::Boltzmann {
            return (half, 0.0);
        }
        let (t, dt) = degeneracy(u + half - chi_c);
        (half + t, dt)
    }
    #[inline]
    fn theta_p(&self, k: usize, w: f64, chi_v: f64) -> (f64, f64) {
        let half = 0.5 * self.dg[k];
        if self.statistics == Statistics::Boltzmann {
            return (half, 0.0);
        }
        let (t, dt) = degeneracy(w + half - chi_v);
        (half + t, dt)
    }

    /// Charge-neutral ψ at φ_n = φ_p = 0 (Newton on p̂ − n̂ + Ĉ, bisection
    /// safeguard).
    fn solve_neutral(&self, k: usize) -> f64 {
        let (cc, cv) = self.chi();
        let c = self.c_hat(k);
        let f = |psi: f64| -> f64 {
            let (tn, _) = self.theta_n(k, psi, cc);
            let (tp, _) = self.theta_p(k, -psi, cv);
            (tp - psi).exp() - (psi + tn).exp() + c
        };
        let (mut lo, mut hi) = (-80.0, 80.0);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if f(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Net doping (N_d − N_a)/n_i at node k.
    #[inline]
    fn c_hat(&self, k: usize) -> f64 {
        (self.nd[k] - self.na[k]) / self.si.ni
    }

    pub fn initial_guess(&mut self) {
        for k in 0..self.psi.len() {
            self.psi[k] = self.psi0[k];
            self.phin[k] = 0.0;
            self.phip[k] = 0.0;
        }
        self.apply_dirichlet();
    }

    fn dirichlet(&self, k: usize) -> Option<[f64; 3]> {
        let c = self.contact_of[k]?;
        let ct = &self.contacts[c];
        let vt = self.si.vt();
        let v = ct.v / vt;
        Some(match &ct.kind {
            ContactKind::Ohmic => [v + self.psi0[k], v, v],
            ContactKind::Gate { offset } => [v + offset / vt, self.phin[k], self.phip[k]],
            ContactKind::Fixed(vals) => vals[ct.nodes.iter().position(|&q| q == k).unwrap()],
        })
    }

    fn apply_dirichlet(&mut self) {
        for k in 0..self.psi.len() {
            if let Some(d) = self.dirichlet(k) {
                self.psi[k] = d[0];
                self.phin[k] = d[1];
                self.phip[k] = d[2];
            }
        }
    }

    fn carriers(&self) -> Carriers {
        let nn = self.psi.len();
        let (cc, cv) = self.chi();
        let mut c = Carriers { n: vec![0.0; nn], p: vec![0.0; nn], sn: vec![0.0; nn], tp: vec![0.0; nn], dthn: vec![0.0; nn], dthp: vec![0.0; nn], eq: vec![1.0; nn] };
        let mut fd = 0;
        for k in 0..nn {
            let (tn, dtn) = self.theta_n(k, self.psi[k] - self.phin[k], cc);
            let (tp, dtp) = self.theta_p(k, self.phip[k] - self.psi[k], cv);
            if self.statistics == Statistics::FermiDirac && (dtn != 0.0 || dtp != 0.0) {
                fd += 1;
            }
            c.sn[k] = self.psi[k] + tn;
            c.tp[k] = self.psi[k] - tp;
            c.n[k] = (c.sn[k] - self.phin[k]).exp();
            c.p[k] = (self.phip[k] - c.tp[k]).exp();
            c.dthn[k] = dtn;
            c.dthp[k] = dtp;
            c.eq[k] = (tn + tp).exp();
        }
        self.fd_nodes.set(fd);
        c
    }

    /// Electron density (m⁻³) at node k.
    pub fn n(&self, k: usize) -> f64 {
        let (cc, _) = self.chi();
        let (tn, _) = self.theta_n(k, self.psi[k] - self.phin[k], cc);
        self.si.ni * (self.psi[k] + tn - self.phin[k]).exp()
    }
    /// Hole density (m⁻³) at node k.
    pub fn p(&self, k: usize) -> f64 {
        let (_, cv) = self.chi();
        let (tp, _) = self.theta_p(k, self.phip[k] - self.psi[k], cv);
        self.si.ni * (self.phip[k] - self.psi[k] + tp).exp()
    }

    /// Recombination R/n_i (1/s) and derivatives with respect to n̂, p̂ and
    /// the equilibrium product E = e^{θ_n+θ_p}.
    fn recomb(&self, k: usize, n: f64, p: f64, e: f64) -> (f64, f64, f64, f64) {
        if !self.recombination {
            return (0.0, 0.0, 0.0, 0.0);
        }
        let ntot = self.nd[k] + self.na[k];
        let (tn, tp) = (self.si.tau_n(ntot), self.si.tau_p(ntot));
        let se = e.sqrt();
        let num = n * p - e;
        let den = tp * (n + se) + tn * (p + se);
        let r = num / den;
        let rn = (p * den - num * tp) / (den * den);
        let rp = (n * den - num * tn) / (den * den);
        let re = (-den - num * (tp + tn) * 0.5 / se) / (den * den);
        let ni2 = self.si.ni * self.si.ni;
        let a = ni2 * (self.si.c_n * n + self.si.c_p * p);
        let m = self.recomb_scale.as_ref().map_or(1.0, |s| s[k]);
        (m * (r + a * num), m * (rn + ni2 * self.si.c_n * num + a * p), m * (rp + ni2 * self.si.c_p * num + a * n), m * (re - a))
    }

    /// Node electric-field components (V/m) from semiconductor neighbours
    /// (one-sided at the oxide interface and boundaries).
    fn node_fields(&self) -> Vec<(f64, f64)> {
        let m = &self.mesh;
        let (nx, ny) = (m.nx(), m.ny());
        let vt = self.si.vt();
        let grad = |k: usize, lo: Option<usize>, hi: Option<usize>, x: &dyn Fn(usize) -> f64| -> f64 {
            let ok = |q: Option<usize>| q.filter(|&q| self.semi[q]);
            match (ok(lo), ok(hi)) {
                (Some(a), Some(b)) => -(self.psi[b] - self.psi[a]) * vt / (x(b) - x(a)),
                (Some(a), None) => -(self.psi[k] - self.psi[a]) * vt / (x(k) - x(a)),
                (None, Some(b)) => -(self.psi[b] - self.psi[k]) * vt / (x(b) - x(k)),
                (None, None) => 0.0,
            }
        };
        (0..self.psi.len())
            .map(|k| {
                if !self.semi[k] {
                    return (0.0, 0.0);
                }
                let (i, j) = m.ij(k);
                let xs = |q: usize| m.x[m.ij(q).0];
                let ys = |q: usize| m.y[m.ij(q).1];
                let ex = grad(k, (i > 0).then(|| m.node(i - 1, j)), (i + 1 < nx).then(|| m.node(i + 1, j)), &xs);
                let ey = if ny > 2 { grad(k, (j > 0).then(|| m.node(i, j - 1)), (j + 1 < ny).then(|| m.node(i, j + 1)), &ys) } else { 0.0 };
                (ex, ey)
            })
            .collect()
    }

    /// Edge mobilities (Lombardi on the perpendicular field, then velocity
    /// saturation on the parallel field) and d ln μ/dψ_b (exact-Jacobian
    /// option only).
    fn edge_mobility_d(&self, e: &Edge, fields: Option<&[(f64, f64)]>) -> (f64, f64, f64, f64) {
        let t = self.si.t;
        let (mut mun, mut mup) = (e.mu_n, e.mu_p);
        if self.screening && self.mu_const.is_none() {
            let ntot = 0.5 * (self.nd[e.a] + self.na[e.a] + self.nd[e.b] + self.na[e.b]);
            let ns = 0.5 * (self.n(e.a) + self.p(e.a) + self.n(e.b) + self.p(e.b));
            mun = self.si.mob_n.mobility_screened(ntot, ns, t, true);
            mup = self.si.mob_p.mobility_screened(ntot, ns, t, false);
        }
        if self.lombardi && self.mu_const.is_none() {
            if let Some(f) = fields {
                let perp = |k: usize| if e.along_x { f[k].1.abs() } else { f[k].0.abs() };
                let ep = 0.5 * (perp(e.a) + perp(e.b));
                let ntot = 0.5 * (self.nd[e.a] + self.na[e.a] + self.nd[e.b] + self.na[e.b]);
                mun = LOMBARDI_N.combine(mun, ep, ntot, t);
                mup = LOMBARDI_P.combine(mup, ep, ntot, t);
            }
        }
        if !self.velocity_saturation {
            return (mun, mup, 0.0, 0.0);
        }
        let vt = self.si.vt();
        // smoothed driving force: F_s = √(F² + F₀²) − F₀ (F₀ = 1e4 V/m,
        // far below saturation) removes the |F|^β cusp of Caughey–Thomas at
        // zero field
        let f0 = 1e4;
        let ct = |mu0: f64, dpsi: f64, vs: f64, beta: f64| -> (f64, f64) {
            let fsig = dpsi * vt / e.h;
            let root = (fsig * fsig + f0 * f0).sqrt();
            let fs = root - f0;
            let xb = (mu0 * fs / vs).powf(beta);
            let mu = mu0 / (1.0 + xb).powf(1.0 / beta);
            let dl = if !self.vsat_exact_jacobian || fs <= 0.0 { 0.0 } else { -(xb / fs) / (1.0 + xb) * fsig / root * vt / e.h };
            (mu, dl)
        };
        let (vn, bn) = vsat_n(t);
        let (vp, bp) = vsat_p(t);
        let dpsi = self.psi[e.b] - self.psi[e.a];
        let (mn, dn) = ct(mun, dpsi, vn, bn);
        let (mp, dp) = ct(mup, dpsi, vp, bp);
        (mn, mp, dn, dp)
    }

    /// Electron and hole fluxes on an edge with partials (chain rule
    /// through the effective potentials).
    fn edge_flux(&self, e: &Edge, c: &Carriers, mob: (f64, f64, f64, f64)) -> (Flux, Flux) {
        let (a, b) = (e.a, e.b);
        let (mun, mup, dlmn, dlmp) = mob;
        // electrons, effective potential s = ψ + θ_n
        let ds = c.sn[b] - c.sn[a];
        let (bp, dbp) = bernoulli(ds);
        let (bm, dbm) = bernoulli(-ds);
        let cn = mun * e.si_w / e.h;
        let g = cn * (c.n[b] * bp - c.n[a] * bm);
        let dg_dsb = cn * (c.n[b] * (bp + dbp) + c.n[a] * dbm) + g * dlmn;
        let dg_dsa = cn * (-c.n[b] * dbp - c.n[a] * (bm + dbm)) - g * dlmn;
        let dg_dfb = -cn * c.n[b] * bp;
        let dg_dfa = cn * c.n[a] * bm;
        let fn_ = Flux {
            g,
            s: cn * (c.n[b] * bp + c.n[a] * bm),
            d: [dg_dsa * (1.0 + c.dthn[a]), dg_dsb * (1.0 + c.dthn[b]), dg_dfa - dg_dsa * c.dthn[a], dg_dfb - dg_dsb * c.dthn[b]],
        };
        // holes, effective potential t = ψ − θ_p
        let dt = c.tp[b] - c.tp[a];
        let (bp, dbp) = bernoulli(dt);
        let (bm, dbm) = bernoulli(-dt);
        let cp = mup * e.si_w / e.h;
        let h = cp * (c.p[a] * bp - c.p[b] * bm);
        let dh_dta = cp * (-c.p[a] * bp - c.p[a] * dbp - c.p[b] * dbm) - h * dlmp;
        let dh_dtb = cp * (c.p[a] * dbp + c.p[b] * bm + c.p[b] * dbm) + h * dlmp;
        let dh_dfa = cp * c.p[a] * bp;
        let dh_dfb = -cp * c.p[b] * bm;
        let fp = Flux {
            g: h,
            s: cp * (c.p[a] * bp + c.p[b] * bm),
            d: [dh_dta * (1.0 + c.dthp[a]), dh_dtb * (1.0 + c.dthp[b]), dh_dfa - dh_dta * c.dthp[a], dh_dfb - dh_dtb * c.dthp[b]],
        };
        (fn_, fp)
    }

    /// Fluxes on every edge for the present state.
    fn all_fluxes(&self, car: &Carriers, fields: Option<&[(f64, f64)]>) -> Vec<(Flux, Flux)> {
        self.edges.iter().map(|e| if e.si_w == 0.0 { (Flux::default(), Flux::default()) } else { self.edge_flux(e, car, self.edge_mobility_d(e, fields)) }).collect()
    }

    /// Impact-ionization generation per node in continuity-row units
    /// (V·α|Ĵ| summed over carriers) with its linearization (column,
    /// derivative), α frozen at the node field. Nodes below the
    /// ionization threshold field are skipped (the fold).
    #[allow(clippy::type_complexity)]
    fn generation(&self, car: &Carriers, fluxes: &[(Flux, Flux)], fields: &[(f64, f64)], fixed: &[bool]) -> Vec<Option<(f64, Vec<(usize, f64)>)>> {
        let nn = self.psi.len();
        let mut eidx = vec![[usize::MAX; 4]; nn]; // x−, x+, y−, y+
        for (ie, e) in self.edges.iter().enumerate() {
            if e.si_w == 0.0 {
                continue;
            }
            if e.along_x {
                eidx[e.a][1] = ie;
                eidx[e.b][0] = ie;
            } else {
                eidx[e.a][3] = ie;
                eidx[e.b][2] = ie;
            }
        }
        let mut out = vec![None; nn];
        let mut count = 0;
        for k in 0..nn {
            if fixed[k] || !self.semi[k] {
                continue;
            }
            let fmag = (fields[k].0.powi(2) + fields[k].1.powi(2)).sqrt();
            let (an, ap) = (alpha_n(fmag), alpha_p(fmag));
            if an == 0.0 && ap == 0.0 {
                continue;
            }
            count += 1;
            let v = self.vol_si[k];
            let mut total = 0.0;
            let mut lin: Vec<(usize, f64)> = Vec::new();
            let nloc = (self.nd[k] + self.na[k]) / self.si.ni;
            for (carrier, alpha) in [(0usize, an), (1usize, ap)] {
                if alpha == 0.0 {
                    continue;
                }
                // a carrier above 10% of the local doping is majority (Debye
                // tail of a depletion edge): its current there sits below
                // double-precision resolution of φ, and the ionization over
                // that ~one Debye length is negligible
                let dens = if carrier == 0 { car.n[k] } else { car.p[k] };
                if dens > 0.1 * nloc && nloc > 0.0 {
                    continue;
                }
                // flux-per-width components along x and y, averaged over the
                // node's edges in that direction, with linear forms
                let comp = |dirs: [usize; 2]| -> (f64, Vec<(usize, f64)>) {
                    let used: Vec<usize> = dirs.iter().map(|&d| eidx[k][d]).filter(|&ie| ie != usize::MAX).collect();
                    let wgt = 1.0 / used.len().max(1) as f64;
                    let mut val = 0.0;
                    let mut l = Vec::new();
                    for ie in used {
                        let e = &self.edges[ie];
                        let fl = if carrier == 0 { fluxes[ie].0 } else { fluxes[ie].1 };
                        val += wgt * fl.g / e.si_w;
                        let phi = 1 + carrier;
                        for (q, col) in [(0, 3 * e.a), (1, 3 * e.b), (2, 3 * e.a + phi), (3, 3 * e.b + phi)] {
                            l.push((col, wgt * fl.d[q] / e.si_w));
                        }
                    }
                    (val, l)
                };
                let (jx, lx) = comp([0, 1]);
                let (jy, ly) = comp([2, 3]);
                let mag = (jx * jx + jy * jy).sqrt();
                if mag == 0.0 {
                    continue;
                }
                total += v * alpha * mag;
                for (l, jc) in [(lx, jx), (ly, jy)] {
                    for (col, d) in l {
                        lin.push((col, v * alpha * jc / mag * d));
                    }
                }
            }
            out[k] = Some((total, lin));
        }
        self.ii_nodes.set(count);
        out
    }

    fn nvar(&self) -> usize {
        3 * self.psi.len()
    }
    fn band(&self) -> Band {
        let w = 3 * self.mesh.ny() + 2;
        Band::new(self.nvar(), w, w)
    }

    /// Residual F and (optionally) Jacobian.
    fn assemble(&self, jac: Option<&mut Band>) -> Vec<f64> {
        self.assemble_scaled(jac).0
    }

    /// Residual, its per-row term magnitudes (Σ|terms|, the round-off
    /// scale of each row), and optionally the Jacobian.
    fn assemble_scaled(&self, jac: Option<&mut Band>) -> (Vec<f64>, Vec<f64>) {
        let nn = self.psi.len();
        let mut f = vec![0.0; 3 * nn];
        let mut fs = vec![0.0; 3 * nn];
        let ni = self.si.ni;
        let vt = self.si.vt();
        let kappa = Q * ni / vt;
        let fixed: Vec<bool> = (0..nn).map(|k| self.contact_of[k].is_some()).collect();
        let mut jac = jac;
        if let Some(j) = jac.as_deref_mut() {
            j.clear();
        }
        macro_rules! add {
            ($r:expr, $c:expr, $v:expr) => {
                if let Some(j) = jac.as_deref_mut() {
                    j.add($r, $c, $v);
                }
            };
        }
        let car = self.carriers();
        let fields = if self.lombardi || self.impact_ionization { Some(self.node_fields()) } else { None };
        let fluxes = self.all_fluxes(&car, fields.as_deref());
        for (ie, e) in self.edges.iter().enumerate() {
            let (a, b) = (e.a, e.b);
            // Poisson
            let g = e.eps_w / e.h;
            let pf = g * (self.psi[a] - self.psi[b]);
            for (r, s, sg) in [(a, b, 1.0), (b, a, -1.0)] {
                if !fixed[r] {
                    f[3 * r] += sg * pf;
                    fs[3 * r] += g * (self.psi[a].abs() + self.psi[b].abs());
                    add!(3 * r, 3 * r, g);
                    add!(3 * r, 3 * s, -g);
                }
            }
            if e.si_w == 0.0 {
                continue;
            }
            let (fnx, fpx) = fluxes[ie];
            for (r, sg) in [(a, 1.0), (b, -1.0)] {
                if fixed[r] {
                    continue;
                }
                f[3 * r + 1] += sg * fnx.g;
                fs[3 * r + 1] += fnx.s;
                fs[3 * r + 2] += fpx.s;
                add!(3 * r + 1, 3 * a, sg * fnx.d[0]);
                add!(3 * r + 1, 3 * b, sg * fnx.d[1]);
                add!(3 * r + 1, 3 * a + 1, sg * fnx.d[2]);
                add!(3 * r + 1, 3 * b + 1, sg * fnx.d[3]);
                f[3 * r + 2] += sg * fpx.g;
                add!(3 * r + 2, 3 * a, sg * fpx.d[0]);
                add!(3 * r + 2, 3 * b, sg * fpx.d[1]);
                add!(3 * r + 2, 3 * a + 2, sg * fpx.d[2]);
                add!(3 * r + 2, 3 * b + 2, sg * fpx.d[3]);
            }
        }
        for k in 0..nn {
            if let Some(dv) = self.dirichlet(k) {
                f[3 * k] = self.psi[k] - dv[0];
                f[3 * k + 1] = self.phin[k] - dv[1];
                f[3 * k + 2] = self.phip[k] - dv[2];
                fs[3 * k] = self.psi[k].abs() + dv[0].abs() + 1.0;
                fs[3 * k + 1] = self.phin[k].abs() + dv[1].abs() + 1.0;
                fs[3 * k + 2] = self.phip[k].abs() + dv[2].abs() + 1.0;
                for v in 0..3 {
                    add!(3 * k + v, 3 * k + v, 1.0);
                }
                continue;
            }
            if !self.semi[k] {
                f[3 * k + 1] = 0.0;
                f[3 * k + 2] = 0.0;
                fs[3 * k + 1] = 1.0;
                fs[3 * k + 2] = 1.0;
                add!(3 * k + 1, 3 * k + 1, 1.0);
                add!(3 * k + 2, 3 * k + 2, 1.0);
            } else {
                let v = self.vol_si[k];
                let (n, p) = (car.n[k], car.p[k]);
                // ∂n̂/∂(ψ, φ_n) and ∂p̂/∂(ψ, φ_p)
                let (dn_psi, dn_phi) = (n * (1.0 + car.dthn[k]), -n * (1.0 + car.dthn[k]));
                let (dp_psi, dp_phi) = (-p * (1.0 + car.dthp[k]), p * (1.0 + car.dthp[k]));
                // Poisson charge
                f[3 * k] -= kappa * v * (p - n + self.c_hat(k));
                fs[3 * k] += kappa * v * (p + n + self.c_hat(k).abs());
                add!(3 * k, 3 * k, -kappa * v * (dp_psi - dn_psi));
                add!(3 * k, 3 * k + 1, kappa * v * dn_phi);
                add!(3 * k, 3 * k + 2, -kappa * v * dp_phi);
                // recombination (E = e^{θ_n+θ_p} depends on ψ, φ through θ)
                let (r, rn, rp, re) = self.recomb(k, n, p, car.eq[k]);
                let e = car.eq[k];
                let (de_psi, de_n, de_p) = (e * (car.dthn[k] - car.dthp[k]), -e * car.dthn[k], e * car.dthp[k]);
                let s = v / vt;
                let r_psi = rn * dn_psi + rp * dp_psi + re * de_psi;
                let r_n = rn * dn_phi + re * de_n;
                let r_p = rp * dp_phi + re * de_p;
                f[3 * k + 1] -= s * r;
                fs[3 * k + 1] += (s * r).abs();
                fs[3 * k + 2] += (s * r).abs();
                add!(3 * k + 1, 3 * k, -s * r_psi);
                add!(3 * k + 1, 3 * k + 1, -s * r_n);
                add!(3 * k + 1, 3 * k + 2, -s * r_p);
                f[3 * k + 2] += s * r;
                add!(3 * k + 2, 3 * k, s * r_psi);
                add!(3 * k + 2, 3 * k + 1, s * r_n);
                add!(3 * k + 2, 3 * k + 2, s * r_p);
                // carrier storage (transients)
                if let Some(tr) = &self.tr {
                    let w = s * tr.a[0] / tr.dt;
                    let hn: f64 = tr.hist.iter().enumerate().map(|(l, h)| tr.a[l + 1] * h.0[k]).sum();
                    let hp: f64 = tr.hist.iter().enumerate().map(|(l, h)| tr.a[l + 1] * h.1[k]).sum();
                    f[3 * k + 1] -= s * (tr.a[0] * n + hn) / tr.dt;
                    fs[3 * k + 1] += s * (tr.a[0] * n + hn.abs()) / tr.dt;
                    fs[3 * k + 2] += s * (tr.a[0] * p + hp.abs()) / tr.dt;
                    add!(3 * k + 1, 3 * k, -w * dn_psi);
                    add!(3 * k + 1, 3 * k + 1, -w * dn_phi);
                    f[3 * k + 2] += s * (tr.a[0] * p + hp) / tr.dt;
                    add!(3 * k + 2, 3 * k, w * dp_psi);
                    add!(3 * k + 2, 3 * k + 2, w * dp_phi);
                }
            }
            if let Some(src) = &self.source {
                for v in 0..3 {
                    f[3 * k + v] += src[k][v];
                }
            }
        }
        // impact ionization: generation α_n|J_n| + α_p|J_p| at the nodes
        if self.impact_ionization {
            let gen = self.generation(&car, &fluxes, fields.as_ref().unwrap(), &fixed);
            for (k, gk) in gen.iter().enumerate() {
                let Some((val, lin)) = gk else { continue };
                let _ = lin;
                f[3 * k + 1] += val;
                f[3 * k + 2] -= val;
                fs[3 * k + 1] += val.abs();
                fs[3 * k + 2] += val.abs();
            }
        }
        (f, fs)
    }

    /// Newton solve at the present contact voltages. Returns the iteration
    /// count, or None if it did not converge (state restored).
    pub fn solve(&mut self) -> Option<usize> {
        let saved = (self.psi.clone(), self.phin.clone(), self.phip.clone());
        // impact ionization unfolds from the converged solution without it:
        // generation built from the spurious currents of early iterates far
        // from the solution would destabilize Newton
        let mut pre = 0;
        if self.impact_ionization {
            self.impact_ionization = false;
            let r = self.solve();
            self.impact_ionization = true;
            pre = r?;
        }
        self.apply_dirichlet();
        let max_it = if self.impact_ionization { 2000 } else { 100 };
        for it in 1..=max_it {
            let mut j = self.band();
            let (mut f, fs) = self.assemble_scaled(Some(&mut j));
            // backward-error test: every residual at the round-off level of
            // its own terms means the state solves a problem perturbed at
            // machine precision — converged, whatever the next update says
            // (slow transients balance storage terms below the round-off of
            // near-cancelling drift and diffusion fluxes)
            let berr = f.iter().zip(&fs).map(|(a, b)| if *b > 0.0 { a.abs() / b } else { a.abs() }).fold(0.0, f64::max);
            if it > 1 && berr < 1e-14 {
                return Some(it - 1 + pre);
            }
            f.iter_mut().for_each(|v| *v = -*v);
            let ok = j.solve(&mut f).is_some();
            if !ok || f.iter().any(|v| !v.is_finite()) {
                break;
            }
            let dmax = f.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            let s = (self.max_step / dmax).min(1.0);
            for k in 0..self.psi.len() {
                self.psi[k] += s * f[3 * k];
                self.phin[k] += s * f[3 * k + 1];
                self.phip[k] += s * f[3 * k + 2];
            }
            // convergence: with impact ionization, quasi-Fermi updates of
            // carriers below 1e-14 of the local doping are not waited for —
            // majority-carrier currents in heavily doped regions sit below
            // double-precision resolution of φ, so the generation built from
            // them is noise that only moves such negligible carriers
            let dconv = if self.impact_ionization {
                let (cc, cv) = self.chi();
                (0..self.psi.len())
                    .map(|k| {
                        let floor = 1e-14 * ((self.nd[k] + self.na[k]) / self.si.ni + 1.0);
                        let (tn, _) = self.theta_n(k, self.psi[k] - self.phin[k], cc);
                        let (tp, _) = self.theta_p(k, self.phip[k] - self.psi[k], cv);
                        let n = (self.psi[k] + tn - self.phin[k]).exp();
                        let p = (self.phip[k] - self.psi[k] + tp).exp();
                        let wn = if n > floor { 1.0 } else { 0.0 };
                        let wp = if p > floor { 1.0 } else { 0.0 };
                        f[3 * k].abs().max(wn * f[3 * k + 1].abs()).max(wp * f[3 * k + 2].abs())
                    })
                    .fold(0.0, f64::max)
            } else {
                dmax
            };
            if dconv * s < self.newton_tol && s == 1.0 {
                return Some(it + pre);
            }
        }
        (self.psi, self.phin, self.phip) = saved;
        None
    }

    /// Ramp contact `c` to voltage `v` in steps of at most `dv`, halving on
    /// failure. Returns the total Newton iterations.
    pub fn ramp(&mut self, c: usize, v: f64, dv: f64) -> Option<usize> {
        let mut total = 0;
        let mut step = dv;
        while (self.contacts[c].v - v).abs() > 1e-12 {
            let v0 = self.contacts[c].v;
            let next = if (v - v0).abs() <= step { v } else { v0 + step * (v - v0).signum() };
            self.contacts[c].v = next;
            match self.solve() {
                Some(it) => {
                    total += it;
                    step = (step * 1.5).min(dv);
                }
                None => {
                    self.contacts[c].v = v0;
                    step *= 0.5;
                    if step < 1e-6 {
                        return None;
                    }
                }
            }
        }
        if total == 0 {
            total = self.solve()?;
        }
        Some(total)
    }

    pub fn contact(&self, name: &str) -> usize {
        self.contacts.iter().position(|c| c.name == name).unwrap_or_else(|| panic!("no contact {name}"))
    }

    // -----------------------------------------------------------------
    // quantum confinement (density gradient along y, per column)

    /// Recompute the quantum potential: in every mesh column, solve the 1D
    /// density-gradient equation 2b u'' = u (2 ln u − s) for u = √n̂ given
    /// the classical exponent s = ψ + θ_n,classical − φ_n, with u small at
    /// the Si/SiO₂ interface (e^{−15} of its classical value: no
    /// penetration) and classical in the bulk; Λ = 2 ln u − s. Returns the
    /// largest change of Λ (units of V_t).
    pub fn update_quantum(&mut self) -> f64 {
        let Some(m_dg) = self.dg_mass else { return 0.0 };
        let vt = self.si.vt();
        let b = crate::tunnel::HBAR.powi(2) / (12.0 * Q * m_dg * crate::tunnel::M0 * vt);
        let (cc, _) = self.chi();
        let m = self.mesh.clone();
        let (nx, ny) = (m.nx(), m.ny());
        let mut change: f64 = 0.0;
        for i in 0..nx {
            // semiconductor nodes of this column, top to bottom
            let col: Vec<usize> = (0..ny).map(|j| m.node(i, j)).filter(|&k| self.semi[k]).collect();
            if col.len() < 3 {
                continue;
            }
            // interface at the top only if the node above is oxide
            let top = col[0];
            let (_, jt) = m.ij(top);
            // (contacts are metal: no confining oxide wall at them)
            let wall = jt > 0 && !self.semi[m.node(i, jt - 1)] && self.contact_of[top].is_none();
            if !wall {
                continue;
            }
            let ys: Vec<f64> = col.iter().map(|&k| m.y[m.ij(k).1]).collect();
            let nc_ = col.len();
            // classical exponent (without the present Λ)
            let sv: Vec<f64> = col
                .iter()
                .map(|&k| {
                    let lam = self.lambda[k];
                    let (tn, _) = self.theta_n(k, self.psi[k] - self.phin[k], cc);
                    self.psi[k] + tn - lam - self.phin[k]
                })
                .collect();
            let w: Vec<f64> = (0..nc_).map(|q| 0.5 * (if q > 0 { ys[q] - ys[q - 1] } else { 0.0 } + if q + 1 < nc_ { ys[q + 1] - ys[q] } else { 0.0 })).collect();
            let mut u: Vec<f64> = col.iter().zip(&sv).map(|(&k, &s0)| (0.5 * (s0 + self.lambda[k])).exp()).collect();
            u[0] = (-15.0 + 0.5 * sv[0]).exp();
            u[nc_ - 1] = (0.5 * sv[nc_ - 1]).exp();
            for _ in 0..60 {
                // Newton on interior nodes (tridiagonal)
                let mut a = Band::new(nc_, 1, 1);
                let mut f = vec![0.0; nc_];
                for q in 0..nc_ {
                    if q == 0 || q == nc_ - 1 {
                        a.add(q, q, 1.0);
                        continue;
                    }
                    let mut lap = 0.0;
                    for r in [q - 1, q + 1] {
                        let h = (ys[r] - ys[q]).abs();
                        lap += 2.0 * b * (u[r] - u[q]) / h;
                        a.add(q, r, 2.0 * b / h);
                        a.add(q, q, -2.0 * b / h);
                    }
                    let uq = u[q].max(1e-300);
                    f[q] = -(lap - w[q] * uq * (2.0 * uq.ln() - sv[q]));
                    a.add(q, q, -w[q] * (2.0 * uq.ln() - sv[q] + 2.0));
                }
                if a.solve(&mut f).is_none() {
                    break;
                }
                let mut rel: f64 = 0.0;
                for q in 1..nc_ - 1 {
                    let un = (u[q] + f[q]).max(0.1 * u[q]).min(10.0 * u[q]);
                    rel = rel.max(((un - u[q]) / u[q]).abs());
                    u[q] = un;
                }
                if rel < 1e-12 {
                    break;
                }
            }
            for (q, &k) in col.iter().enumerate() {
                let lam = 2.0 * u[q].ln() - sv[q];
                change = change.max((lam - self.lambda[k]).abs());
                self.lambda[k] = lam;
            }
        }
        change
    }

    /// Self-consistent solve with the density-gradient quantum potential
    /// (outer iteration: device solve, then column density-gradient
    /// solves, until Λ changes by less than 1e-5 V_t; the Λ step is halved
    /// when a device solve fails).
    pub fn solve_quantum(&mut self) -> Option<usize> {
        let mut total = self.solve()?;
        if self.dg_mass.is_none() {
            return Some(total);
        }
        for _ in 0..300 {
            let old = self.lambda.clone();
            let ch = self.update_quantum();
            let target = self.lambda.clone();
            // full outer step; halved if the device solve fails
            let mut frac = 1.0;
            let mut ok = None;
            for _ in 0..12 {
                for k in 0..old.len() {
                    self.lambda[k] = old[k] + frac * (target[k] - old[k]);
                }
                ok = self.solve();
                if ok.is_some() {
                    break;
                }
                frac *= 0.5;
            }
            total += ok?;
            if ch < 1e-5 {
                return Some(total);
            }
        }
        None
    }

    // -----------------------------------------------------------------
    // transients

    /// Start time stepping from the present (steady) state.
    pub fn begin_transient(&mut self) {
        let c = self.carriers();
        self.tr = Some(Transient { dt: 0.0, a: [1.0, -1.0, 0.0], hist: vec![(c.n, c.p, self.psi.clone())], committed: (self.psi.clone(), self.phin.clone(), self.phip.clone()), dt_prev: 0.0, used: None });
    }
    /// Leave transient mode (back to steady-state solves).
    pub fn end_transient(&mut self) {
        self.tr = None;
    }
    /// Solve the time step of length dt to the present contact voltages
    /// (variable-step BDF2 once two levels of history exist, BDF1 before).
    /// The state is not committed: call `accept_step` after a successful
    /// solve, or `reject_step` to retry.
    pub fn step(&mut self, dt: f64) -> Option<usize> {
        let tr = self.tr.as_mut().expect("begin_transient first");
        tr.dt = dt;
        tr.used = None;
        tr.a = if tr.hist.len() >= 2 && tr.dt_prev > 0.0 {
            let w = dt / tr.dt_prev;
            [(1.0 + 2.0 * w) / (1.0 + w), -(1.0 + w), w * w / (1.0 + w)]
        } else {
            [1.0, -1.0, 0.0]
        };
        self.solve()
    }
    /// Commit the solved step to the history.
    pub fn accept_step(&mut self) {
        let c = self.carriers();
        let state = (self.psi.clone(), self.phin.clone(), self.phip.clone());
        let tr = self.tr.as_mut().expect("begin_transient first");
        tr.used = Some((tr.hist.clone(), tr.a, tr.dt));
        tr.hist.insert(0, (c.n, c.p, state.0.clone()));
        tr.hist.truncate(2);
        tr.committed = state;
        tr.dt_prev = tr.dt;
    }
    /// Restore the state of the last committed level (after a rejected
    /// step).
    pub fn reject_step(&mut self) {
        let tr = self.tr.as_ref().expect("begin_transient first");
        (self.psi, self.phin, self.phip) = tr.committed.clone();
    }

    /// Advance from the committed level by dt, with contact voltages set by
    /// `bias(dev, s)` at the fraction s ∈ (0, 1] of the interval; halves
    /// the step on Newton failure (down to dt/2¹⁶). Returns the sub-steps.
    pub fn advance(&mut self, dt: f64, bias: &dyn Fn(&mut Device, f64)) -> Option<usize> {
        let mut done = 0.0f64;
        let mut h = 1.0f64;
        let mut count = 0;
        while done < 1.0 - 1e-12 {
            let s = (done + h).min(1.0);
            bias(self, s);
            if self.step((s - done) * dt).is_some() {
                self.accept_step();
                done = s;
                count += 1;
                h = (h * 2.0).min(1.0);
            } else {
                self.reject_step();
                h *= 0.5;
                if h < 1.0 / 65536.0 {
                    return None;
                }
            }
        }
        Some(count)
    }

    /// Total terminal current during a transient (A, into the device):
    /// conduction plus displacement, weighted by the contact's Laplace
    /// function. The discrete equations conserve total current exactly,
    /// so this equals the current through the contact.
    pub fn transient_current(&self, c: usize) -> f64 {
        let tr = self.tr.as_ref().expect("begin_transient first");
        // the step just solved (not yet committed) or the last committed one
        let (hist, a, dt) = match &tr.used {
            Some((h, a, dt)) => (h, *a, *dt),
            None => (&tr.hist, tr.a, tr.dt),
        };
        let g = self.weight(c);
        let car = self.carriers();
        let fields = if self.lombardi { Some(self.node_fields()) } else { None };
        let vt = self.si.vt();
        let (mut cond, mut disp) = (0.0, 0.0);
        for e in &self.edges {
            let dg = g[e.a] - g[e.b];
            if dg == 0.0 {
                continue;
            }
            if e.si_w > 0.0 {
                let (fnx, fpx) = self.edge_flux(e, &car, self.edge_mobility_d(e, fields.as_deref()));
                cond += dg * (fnx.g + fpx.g);
            }
            let d_now = self.psi[e.a] - self.psi[e.b];
            let hsum: f64 = hist.iter().enumerate().map(|(l, h)| a[l + 1] * (h.2[e.a] - h.2[e.b])).sum();
            disp += dg * e.eps_w / e.h * (a[0] * d_now + hsum) / dt;
        }
        (Q * self.si.ni * vt * cond + vt * disp) * self.depth
    }

    // -----------------------------------------------------------------
    // terminal quantities

    /// Weighting function of contact c: Laplace's equation on the mesh
    /// with g = 1 on c and 0 on the other contacts.
    fn weight(&self, c: usize) -> Vec<f64> {
        let nn = self.psi.len();
        let w = self.mesh.ny();
        let mut a = Band::new(nn, w, w);
        let mut b = vec![0.0; nn];
        for e in &self.edges {
            let g = e.eps_w / e.h;
            for (r, s) in [(e.a, e.b), (e.b, e.a)] {
                if self.contact_of[r].is_none() {
                    a.add(r, r, g);
                    a.add(r, s, -g);
                }
            }
        }
        for k in 0..nn {
            if let Some(ck) = self.contact_of[k] {
                a.add(k, k, 1.0);
                b[k] = (ck == c) as u8 as f64;
            }
        }
        a.solve(&mut b).expect("weighting function");
        b
    }

    /// Current flowing from contact c into the device (A), electrons and
    /// holes (steady state). Evaluated as I = Σ_edges J·(g_a − g_b) ∓
    /// Σ_nodes g q (R − G) V with the contact's weighting function g:
    /// identical to the flux through the contact when current is
    /// conserved, but insensitive to the round-off of majority-carrier
    /// fluxes (tiny differences of huge drift and diffusion terms).
    pub fn current(&self, c: usize) -> (f64, f64) {
        let g = self.weight(c);
        let ni = self.si.ni;
        let vt = self.si.vt();
        let car = self.carriers();
        let fields = if self.lombardi || self.impact_ionization { Some(self.node_fields()) } else { None };
        let (mut jn, mut jp) = (0.0, 0.0);
        for e in &self.edges {
            let dg = g[e.a] - g[e.b];
            if e.si_w == 0.0 || dg == 0.0 {
                continue;
            }
            let (fnx, fpx) = self.edge_flux(e, &car, self.edge_mobility_d(e, fields.as_deref()));
            jn += dg * fnx.g;
            jp += dg * fpx.g;
        }
        // net recombination (and impact generation) at the non-contact nodes
        let gen = if self.impact_ionization {
            let fixed: Vec<bool> = (0..self.psi.len()).map(|k| self.contact_of[k].is_some()).collect();
            self.generation(&car, &self.all_fluxes(&car, fields.as_deref()), fields.as_ref().unwrap(), &fixed)
        } else {
            vec![None; self.psi.len()]
        };
        for k in 0..self.psi.len() {
            if self.contact_of[k].is_some() || !self.semi[k] || g[k] == 0.0 {
                continue;
            }
            let (r, _, _, _) = self.recomb(k, car.n[k], car.p[k], car.eq[k]);
            let net = self.vol_si[k] * r / vt - gen[k].as_ref().map_or(0.0, |x| x.0);
            jn -= g[k] * net;
            jp += g[k] * net;
        }
        let s = Q * ni * vt * self.depth;
        (s * jn, s * jp)
    }

    /// Electron and hole current densities (A/m², +x direction) on the x
    /// edges of node row j: (x of the edge midpoint, J_n, J_p).
    pub fn x_currents(&self, j: usize) -> Vec<(f64, f64, f64)> {
        let ni = self.si.ni;
        let vt = self.si.vt();
        let m = &self.mesh;
        let car = self.carriers();
        let fields = if self.lombardi { Some(self.node_fields()) } else { None };
        (0..m.nx() - 1)
            .map(|i| {
                let (a, b) = (m.node(i, j), m.node(i + 1, j));
                let e = self.edges.iter().find(|e| e.a == a && e.b == b).unwrap();
                let (fnx, fpx) = self.edge_flux(e, &car, self.edge_mobility_d(e, fields.as_deref()));
                let s = Q * ni * vt / e.si_w.max(1e-300);
                (0.5 * (m.x[i] + m.x[i + 1]), s * fnx.g, s * fpx.g)
            })
            .collect()
    }

    /// Charge on contact c (C): Gauss's law over its nodes.
    pub fn charge(&self, c: usize) -> f64 {
        let vt = self.si.vt();
        let mut q = 0.0;
        for e in &self.edges {
            let (ca, cb) = (self.contact_of[e.a], self.contact_of[e.b]);
            let g = e.eps_w / e.h * (self.psi[e.a] - self.psi[e.b]);
            if ca == Some(c) && cb != Some(c) {
                q += g;
            } else if cb == Some(c) && ca != Some(c) {
                q -= g;
            }
        }
        q * vt * self.depth
    }

    /// Quasi-static terminal charges (C), one per contact, summing to zero:
    /// Q_c = (Gauss charge of electrode c) + Σ_k g^σ_c(k) ρ_k V_k, where
    /// g^σ_c solves the conduction-weighted Laplace equation (edge
    /// conductances from the Scharfetter–Gummel fluxes): the share of a
    /// change of space charge that is supplied through contact c by the
    /// charging currents (Ward–Dutton, generalized to 2D).
    pub fn qs_charges(&self) -> Vec<f64> {
        let nn = self.psi.len();
        let car = self.carriers();
        let fields = if self.lombardi { Some(self.node_fields()) } else { None };
        let w = self.mesh.ny();
        // conductances
        let cond: Vec<f64> = self
            .edges
            .iter()
            .map(|e| {
                if e.si_w == 0.0 {
                    return 0.0;
                }
                let (fnx, fpx) = self.edge_flux(e, &car, self.edge_mobility_d(e, fields.as_deref()));
                // ∂G/∂φ at fixed ψ (symmetric part)
                0.5 * (fnx.d[2].abs() + fnx.d[3].abs()) + 0.5 * (fpx.d[2].abs() + fpx.d[3].abs())
            })
            .collect();
        let gmax = cond.iter().cloned().fold(0.0, f64::max);
        let floor = 1e-12 * gmax;
        let rho: Vec<f64> = (0..nn).map(|k| if self.semi[k] { Q * self.si.ni * self.vol_si[k] * (car.p[k] - car.n[k] + self.c_hat(k)) * self.depth } else { 0.0 }).collect();
        let mut out = Vec::with_capacity(self.contacts.len());
        for (ci, ct) in self.contacts.iter().enumerate() {
            let mut q = self.charge(ci);
            if !matches!(ct.kind, ContactKind::Gate { .. }) {
                let mut a = Band::new(nn, w, w);
                let mut b = vec![0.0; nn];
                for (e, &gc) in self.edges.iter().zip(&cond) {
                    if e.si_w == 0.0 {
                        continue;
                    }
                    let g = gc + floor;
                    for (r, s) in [(e.a, e.b), (e.b, e.a)] {
                        if self.contact_of[r].is_none() && self.semi[r] {
                            a.add(r, r, g);
                            a.add(r, s, -g);
                        }
                    }
                }
                for k in 0..nn {
                    if let Some(ck) = self.contact_of[k] {
                        a.add(k, k, 1.0);
                        let is_carrier_contact = !matches!(self.contacts[ck].kind, ContactKind::Gate { .. });
                        b[k] = (ck == ci && is_carrier_contact) as u8 as f64;
                    } else if !self.semi[k] {
                        a.add(k, k, 1.0);
                    }
                }
                a.solve(&mut b).expect("conduction weighting");
                q += (0..nn).map(|k| b[k] * rho[k]).sum::<f64>();
            }
            out.push(q);
        }
        out
    }

    /// Joule and recombination heat per node (W): H = −J_n·∇φ_n − J_p·∇φ_p
    /// on the edges (half to each end) plus q(R − G)(φ_p − φ_n) at the
    /// nodes. Its sum equals Σ_c I_c V_c exactly (discrete energy balance).
    pub fn heat(&self) -> Vec<f64> {
        let nn = self.psi.len();
        let car = self.carriers();
        let fields = if self.lombardi || self.impact_ionization { Some(self.node_fields()) } else { None };
        let vt = self.si.vt();
        let s = Q * self.si.ni * vt * vt * self.depth;
        let mut h = vec![0.0; nn];
        for e in &self.edges {
            if e.si_w == 0.0 {
                continue;
            }
            let (fnx, fpx) = self.edge_flux(e, &car, self.edge_mobility_d(e, fields.as_deref()));
            let p = s * (fnx.g * (self.phin[e.a] - self.phin[e.b]) + fpx.g * (self.phip[e.a] - self.phip[e.b]));
            h[e.a] += 0.5 * p;
            h[e.b] += 0.5 * p;
        }
        let gen = if self.impact_ionization {
            let fixed: Vec<bool> = (0..nn).map(|k| self.contact_of[k].is_some()).collect();
            self.generation(&car, &self.all_fluxes(&car, fields.as_deref()), fields.as_ref().unwrap(), &fixed)
        } else {
            vec![None; nn]
        };
        for k in 0..nn {
            if !self.semi[k] || self.contact_of[k].is_some() {
                continue;
            }
            // net recombination in continuity-row units (V R̂/V_t − V ĝ/V_t)
            let (r, _, _, _) = self.recomb(k, car.n[k], car.p[k], car.eq[k]);
            let net = self.vol_si[k] * r / vt - gen[k].as_ref().map_or(0.0, |x| x.0);
            h[k] += s * net * (self.phip[k] - self.phin[k]);
        }
        h
    }

    /// Total mobile + fixed charge in the semiconductor (C).
    pub fn semiconductor_charge(&self) -> f64 {
        (0..self.psi.len()).filter(|&k| self.semi[k]).map(|k| Q * self.vol_si[k] * (self.p(k) - self.n(k) + self.nd[k] - self.na[k])).sum::<f64>() * self.depth
    }

    /// Full residual vector (diagnostic).
    pub fn residual_vec(&self) -> Vec<f64> {
        self.assemble(None)
    }

    /// Residual norm (diagnostic).
    pub fn residual(&self) -> f64 {
        self.assemble(None).iter().fold(0.0f64, |m, v| m.max(v.abs()))
    }

    pub fn semiconductor_volume(&self, k: usize) -> f64 {
        self.vol_si[k]
    }
    pub fn is_semiconductor(&self, k: usize) -> bool {
        self.semi[k]
    }
    pub fn contact_of(&self, k: usize) -> Option<usize> {
        self.contact_of[k]
    }
    /// Node electric field magnitude (V/m).
    pub fn field_magnitude(&self) -> Vec<f64> {
        self.node_fields().iter().map(|f| (f.0 * f.0 + f.1 * f.1).sqrt()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::{uniform, Mesh};

    /// Jacobian against central differences of the residual at a non-
    /// equilibrium state: Boltzmann, Fermi–Dirac + band-gap narrowing on
    /// degenerate doping, velocity saturation (exact option), transients.
    #[test]
    fn jacobian_matches_finite_differences() {
        for case in 0..4 {
            let mut y = vec![-4e-9, -2e-9];
            y.extend(uniform(0.0, 0.1e-6, 6));
            let mesh = Mesh::new(uniform(0.0, 0.2e-6, 12), y, |_, y| if y < 0.0 { Material::Oxide } else { Material::Silicon });
            let ny = mesh.ny();
            let left: Vec<usize> = (2..ny).map(|j| mesh.node(0, j)).collect();
            let right: Vec<usize> = (2..ny).map(|j| mesh.node(mesh.nx() - 1, j)).collect();
            let gate: Vec<usize> = (3..9).map(|i| mesh.node(i, 0)).collect();
            let contacts = vec![
                Contact { name: "s".into(), nodes: left, kind: ContactKind::Ohmic, v: 0.0 },
                Contact { name: "d".into(), nodes: right, kind: ContactKind::Ohmic, v: 0.3 },
                Contact { name: "g".into(), nodes: gate, kind: ContactKind::Gate { offset: 0.56 }, v: 1.0 },
            ];
            let nsd = if case >= 2 { 2e26 } else { 1e25 };
            let mut dev = Device::new(mesh, Silicon::new(300.0), move |x, _| if !(0.05e-6..=0.15e-6).contains(&x) { (nsd, 1e23) } else { (0.0, 1e23) }, contacts, 1e-6);
            match case {
                1 => {
                    dev.velocity_saturation = true;
                    dev.vsat_exact_jacobian = true;
                }
                2 => {
                    dev.statistics = Statistics::FermiDirac;
                    dev.bgn = true;
                    dev.update_statistics();
                }
                3 => {
                    dev.statistics = Statistics::FermiDirac;
                    dev.bgn = true;
                    dev.update_statistics();
                    dev.begin_transient();
                    dev.tr.as_mut().unwrap().dt = 1e-12;
                }
                _ => {}
            }
            dev.apply_dirichlet();
            for k in 0..dev.psi.len() {
                let (i, j) = dev.mesh.ij(k);
                if dev.contact_of[k].is_none() {
                    dev.psi[k] += 0.7 * ((i as f64) * 0.9 + j as f64 * 0.4).sin();
                    dev.phin[k] += 2.0 + 0.3 * (i as f64 * 0.5).cos();
                    dev.phip[k] -= 0.4 * (j as f64 * 0.7).sin();
                }
            }
            let mut jac = dev.band();
            let _ = dev.assemble(Some(&mut jac));
            assert!(case < 2 || dev.fd_nodes.get() > 0, "Fermi–Dirac not exercised");
            let nv = dev.nvar();
            let mut worst = 0.0f64;
            for col in 0..nv {
                let (k, v) = (col / 3, col % 3);
                let h = 1e-6;
                let bump = |d: &mut Device, delta: f64| match v {
                    0 => d.psi[k] += delta,
                    1 => d.phin[k] += delta,
                    _ => d.phip[k] += delta,
                };
                bump(&mut dev, h);
                let fp = dev.assemble(None);
                bump(&mut dev, -2.0 * h);
                let fm = dev.assemble(None);
                bump(&mut dev, h);
                for row in 0..nv {
                    let gate_node = matches!(dev.contact_of[row / 3].map(|c| &dev.contacts[c].kind), Some(ContactKind::Gate { .. }));
                    if row % 3 != 0 && (!dev.semi[row / 3] || gate_node) {
                        continue;
                    }
                    let fd = (fp[row] - fm[row]) / (2.0 * h);
                    let an = jac.get(row, col);
                    let rowmax = (0..nv).map(|c| jac.get(row, c).abs()).fold(0.0, f64::max).max(1e-300);
                    let err = (fd - an).abs() / rowmax;
                    worst = worst.max(err);
                    assert!(err < 1e-5, "case {case}: row {row} col {col} (node {k} var {v}): analytic {an:e}, FD {fd:e}");
                }
            }
            eprintln!("case {case}: worst relative Jacobian error {worst:.2e} (Fermi–Dirac nodes {})", dev.fd_nodes.get());
        }
    }
}
