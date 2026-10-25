//! Pseudo-two-dimensional Doyle–Fuller–Newman (DFN) model: the reference rung
//! above the electrolyte-enhanced SPM.
//!
//! The cell is resolved through the thickness (`x`) into negative electrode |
//! separator | positive electrode, and every electrode node carries a spherical
//! particle resolved in the radial coordinate (`r`) — hence *pseudo-2D*. The
//! model couples
//!
//! ```text
//!   solid diffusion      ∂c_s/∂t = (D_s/r²) ∂_r(r² ∂_r c_s),
//!                        −D_s ∂_r c_s|_R = j/F
//!   electrolyte mass     ε_e ∂c_e/∂t = ∂_x(D_e^eff ∂_x c_e) + (1−t₊) a_s j/F
//!   electrolyte charge   ∂_x i_e = a_s j,   i_e = −κ^eff ∂_x φ_e + κ_D^eff ∂_x ln c_e
//!   solid charge         ∂_x(σ^eff ∂_x φ_s) = a_s j
//!   Butler–Volmer        j = i₀[exp(α_a F η/RT) − exp(−α_c F η/RT)],
//!                        η = φ_s − φ_e − U(c_s^surf),   i₀ = k F √(c_e c_s(c_max−c_s))
//! ```
//!
//! with the concentrated-solution diffusion-potential coefficient
//! `κ_D^eff/κ^eff = 2RT(1−t₊)/F` and Bruggeman corrections
//! `D_e^eff = D_e ε_e^b`, `κ^eff = κ ε_e^b`. The electrolyte current is
//! `i_e`, the solid-phase current `i_s = −σ^eff ∂_x φ_s`, and `i_s + i_e = I`
//! with `I = i_app` the applied (discharge-positive) current density.
//!
//! The nonlinear system is discretised with **finite volumes** on a variable
//! through-thickness grid and **backward (implicit) Euler** in time, then
//! solved by a **Newton iteration** on the whole coupled state
//! `(φ_e, φ_s, j, c_e, c_s)`; the Jacobian is formed by forward differences and
//! the dense linear system is row/column equilibrated before Gaussian
//! elimination. The electrolyte potential is pinned at the negative current
//! collector, and the reaction kinetics fix the solid potential through the
//! per-electrode charge balance.
//!
//! # References
//! * M. Doyle, T. F. Fuller, J. Newman, *J. Electrochem. Soc.* **140**, 1526
//!   (1993) — the DFN model.
//! * K. E. Thomas, J. Newman, R. M. Darling, *Advances in Lithium-Ion
//!   Batteries* (2002) — concentrated-solution transport.
//! * C.-H. Chen et al., *J. Electrochem. Soc.* **167**, 080534 (2020) — LG M50
//!   (`Chen2020`) parameter set.
//! * F. Marquis et al., *J. Electrochem. Soc.* **166**, A3693 (2019) — SPMe and
//!   its error relative to the DFN.
//!
//! # Honesty note
//! This crate only carries the published graphite OCP (negative) and a
//! **surrogate** NMC sigmoid (positive); the electrolyte transport properties,
//! porosities, Bruggeman exponent, ionic conductivity and solid conductivities
//! are flagged `ESTIMATE`. The DFN therefore reproduces the placement of the
//! fold (reaction non-uniformity, electrolyte limitation) and the
//! conservation/charge-balance laws, not the published LG M50 voltage digits.

use crate::spm::{ocp_graphite, ocp_nmc, CellParams};
use crate::{F_FARADAY, R_GAS};

const NEG: usize = 0;
const SEP: usize = 1;
const POS: usize = 2;

/// DFN material and cell parameters. The electrochemical half-cell values
/// (particle radii, solid diffusivities, maximum/initial concentrations and
/// kinetics) come from [`CellParams`]; the electrolyte and matrix properties
/// are added here.
#[derive(Clone, Debug)]
pub struct DfnParams {
    /// Half-cell parameters (Chen 2020 LG M50, as used by the SPM/SPMe).
    pub cell: CellParams,
    /// Separator thickness (m). `ESTIMATE` (matches `spme`).
    pub l_sep: f64,
    /// Electrolyte volume fraction `[negative, separator, positive]`.
    /// `ESTIMATE`, chosen compatible with the crate's active fractions.
    pub eps_e: [f64; 3],
    /// Bruggeman exponent `b` in `ε_e^b`. `ESTIMATE` (1.5).
    pub bruggeman: f64,
    /// Electrolyte salt diffusivity (m² s⁻¹). `ESTIMATE`.
    pub d_e: f64,
    /// Ionic conductivity at `c_e` (S m⁻¹). `ESTIMATE`.
    pub kappa: f64,
    /// Cation transference number `t₊^0`. Chen 2020.
    pub t_plus: f64,
    /// Solid matrix electronic conductivity `[negative, positive]` (S m⁻¹).
    /// `ESTIMATE`.
    pub sigma: [f64; 2],
}

impl DfnParams {
    /// LG M50 (`Chen2020`)-based DFN parameter set.
    pub fn lg_m50() -> Self {
        Self {
            cell: CellParams::lg_m50(),
            l_sep: 12.0e-6, // ESTIMATE separator thickness (as in spme.rs)
            // ESTIMATE porosities; separator from PyBaMM Chen2020, electrodes
            // set so that eps_s + eps_e <= 1.
            eps_e: [0.25, 0.47, 0.335],
            bruggeman: 1.5, // ESTIMATE
            d_e: 2.0e-10,   // ESTIMATE LiPF6/carbonate, as in spme.rs
            kappa: 1.0,     // ESTIMATE ~1 S/m at 1 M
            t_plus: 0.2594, // Chen 2020
            sigma: [100.0, 10.0], // ESTIMATE matrix conductivities
        }
    }
}

/// A running DFN cell.
#[derive(Clone, Debug)]
pub struct Dfn {
    pub params: DfnParams,

    // ---- through-thickness grid ----
    n: usize,
    n_el: usize,
    x: Vec<f64>,
    dx: Vec<f64>,
    region: Vec<usize>,
    cell_el: Vec<Option<usize>>,
    el_global: Vec<usize>,
    el_region: Vec<usize>,
    el_local: Vec<usize>,
    el_ne: usize,
    el_np: usize,

    // ---- radial grid per electrode (index 0 = negative, 1 = positive) ----
    m: usize,
    rvol: [Vec<f64>; 2],
    rarea: [Vec<f64>; 2],
    rdr: [f64; 2],

    // ---- state ----
    /// Electrolyte potential (V) on every through-thickness node; `phi_e[0]`
    /// is the pinned reference and is held at 0.
    pub phi_e: Vec<f64>,
    /// Solid-phase potential (V) on the electrode nodes only.
    pub phi_s: Vec<f64>,
    /// Interfacial transfer current density `j` (A m⁻²), anodic positive.
    pub j: Vec<f64>,
    /// Electrolyte concentration (mol m⁻³) on every node.
    pub c_e: Vec<f64>,
    /// Solid concentration (mol m⁻³), flat `[el cell][shell]`.
    pub c_s: Vec<f64>,
    /// Elapsed time (s).
    pub t: f64,

    // ---- linear layout ----
    off_pe: usize,
    off_ps: usize,
    off_j: usize,
    off_ce: usize,
    off_cs: usize,
    ntot: usize,
}

/// Effective material coefficient at through-thickness node `i`.
impl Dfn {
    fn eps_e(&self, i: usize) -> f64 {
        self.params.eps_e[self.region[i]]
    }

    /// `D_e^eff = D_e ε_e^b` (m² s⁻¹).
    fn d_eff(&self, i: usize) -> f64 {
        self.params.d_e * self.eps_e(i).powf(self.params.bruggeman)
    }

    /// `κ^eff = κ ε_e^b` (S m⁻¹).
    fn kappa_eff(&self, i: usize) -> f64 {
        self.params.kappa * self.eps_e(i).powf(self.params.bruggeman)
    }

    /// `κ_D^eff / κ^eff = 2RT(1−t₊)/F` (V).
    fn kd_over_kappa(&self) -> f64 {
        2.0 * R_GAS * self.params.cell.t / F_FARADAY * (1.0 - self.params.t_plus)
    }

    /// Specific active surface area `a_s = 3 ε_s / R` at node `i` (m² m⁻³),
    /// zero in the separator.
    fn a_s(&self, i: usize) -> f64 {
        match self.region[i] {
            NEG => 3.0 * self.params.cell.neg.eps_s / self.params.cell.neg.radius,
            POS => 3.0 * self.params.cell.pos.eps_s / self.params.cell.pos.radius,
            _ => 0.0,
        }
    }

    /// Series (harmonic) face diffusivity between nodes `i-1` and `i`.
    fn d_face(&self, i: usize) -> f64 {
        let h = self.x[i] - self.x[i - 1];
        let dl = self.d_eff(i - 1);
        let dr = self.d_eff(i);
        h / (0.5 * self.dx[i - 1] / dl + 0.5 * self.dx[i] / dr)
    }

    /// Series (harmonic) face ionic conductivity between nodes `i-1` and `i`.
    fn k_face(&self, i: usize) -> f64 {
        let h = self.x[i] - self.x[i - 1];
        let kl = self.kappa_eff(i - 1);
        let kr = self.kappa_eff(i);
        h / (0.5 * self.dx[i - 1] / kl + 0.5 * self.dx[i] / kr)
    }

    /// Electrolyte current density at the interior face `i` (between nodes
    /// `i-1` and `i`), positive in `+x`.
    fn i_e_face(&self, i: usize, phi_e: &[f64], c_e: &[f64]) -> f64 {
        let h = self.x[i] - self.x[i - 1];
        let k = self.k_face(i);
        let kd = k * self.kd_over_kappa();
        let dphi = phi_e[i] - phi_e[i - 1];
        let dln = safe_ln(c_e[i]) - safe_ln(c_e[i - 1]);
        -k * dphi / h + kd * dln / h
    }

    /// Diffusive electrolyte flux `D^eff ∂_x c_e` at face `i`.
    fn d_flux_face(&self, i: usize, c_e: &[f64]) -> f64 {
        let h = self.x[i] - self.x[i - 1];
        self.d_face(i) * (c_e[i] - c_e[i - 1]) / h
    }

    /// Open-circuit potential (V) of the electrode owning node `i` at solid
    /// surface stoichiometry `x`.
    fn ocp(&self, region: usize, x: f64) -> f64 {
        match region {
            NEG => ocp_graphite(x),
            _ => ocp_nmc(x),
        }
    }

    /// Charge-transfer (exchange) current density `i₀` (A m⁻²) for an
    /// electrode with surface concentration `c_s^surf` and electrolyte `c_e`.
    fn i0(&self, region: usize, c_s_surf: f64, c_e: f64) -> f64 {
        let el = match region {
            NEG => &self.params.cell.neg,
            _ => &self.params.cell.pos,
        };
        let cs = c_s_surf.clamp(1e-6 * el.c_max, el.c_max - 1e-6 * el.c_max);
        let arg = (c_e.max(1.0) * cs * (el.c_max - cs)).max(0.0).sqrt();
        el.k * F_FARADAY * arg
    }

    /// Terminal voltage (V): solid potential at the positive current collector
    /// minus the negative current collector.
    pub fn voltage(&self) -> f64 {
        self.phi_s[self.n_el - 1] - self.phi_s[0]
    }

    /// Total solid lithium per unit area (mol m⁻²).
    pub fn solid_lithium(&self) -> f64 {
        let mut s = 0.0;
        for e in 0..self.n_el {
            let g = self.el_global[e];
            let eps_s = if self.el_region[e] == NEG {
                self.params.cell.neg.eps_s
            } else {
                self.params.cell.pos.eps_s
            };
            let avg = self.particle_average(e);
            s += eps_s * self.dx[g] * avg;
        }
        s
    }

    /// Total lithium in the electrolyte per unit area (mol m⁻²).
    pub fn electrolyte_lithium(&self) -> f64 {
        (0..self.n)
            .map(|i| self.eps_e(i) * self.dx[i] * self.c_e[i])
            .sum()
    }

    /// Integrated reaction current density `∫ a_s j dx` (A m⁻²) in the
    /// negative and positive electrode. Charge balance requires it to equal
    /// `(+i_app, −i_app)` for the discharge-positive convention.
    pub fn integrated_reaction_current(&self) -> (f64, f64) {
        let mut neg = 0.0;
        let mut pos = 0.0;
        for e in 0..self.n_el {
            let g = self.el_global[e];
            let v = self.a_s(g) * self.j[e] * self.dx[g];
            if self.el_region[e] == NEG {
                neg += v;
            } else {
                pos += v;
            }
        }
        (neg, pos)
    }

    /// Volume average of the solid concentration in electrode-cell `e`.
    pub fn particle_average(&self, e: usize) -> f64 {
        let er = self.el_region[e];
        let vol = &self.rvol[if er == NEG { 0 } else { 1 }];
        let num: f64 = (0..self.m)
            .map(|k| self.c_s[e * self.m + k] * vol[k])
            .sum();
        let den: f64 = vol.iter().sum();
        num / den
    }

    /// Surface stoichiometry `c_s^surf / c_max` of electrode-cell `e`.
    pub fn surface_stoichiometry(&self, e: usize) -> f64 {
        let er = self.el_region[e];
        let cmax = if er == NEG {
            self.params.cell.neg.c_max
        } else {
            self.params.cell.pos.c_max
        };
        self.c_s[e * self.m + self.m - 1] / cmax
    }

    /// Reaction non-uniformity indicator: the largest per-electrode standard
    /// deviation of `a_s j` divided by that electrode's mean magnitude (zero
    /// when the reaction is uniform). On discharge the two electrodes carry
    /// opposite currents, so the ratio is taken separately.
    pub fn reaction_nonuniformity(&self) -> f64 {
        let mut worst = 0.0f64;
        for er in [NEG, POS] {
            let vals: Vec<f64> = (0..self.n_el)
                .filter(|&e| self.el_region[e] == er)
                .map(|e| self.a_s(self.el_global[e]) * self.j[e])
                .collect();
            let mean = vals.iter().sum::<f64>() / vals.len() as f64;
            if mean.abs() < 1e-30 {
                continue;
            }
            let var =
                vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / vals.len() as f64;
            worst = worst.max(var.sqrt() / mean.abs());
        }
        worst
    }
}

/// `ln c_e` with a floor so the Newton step cannot feed a non-positive
/// concentration into the diffusion potential.
fn safe_ln(c: f64) -> f64 {
    c.max(1.0).ln()
}

impl Dfn {
    /// Build a DFN cell with `n_n`, `n_sep`, `n_p` finite volumes through the
    /// negative, separator and positive regions and `m` radial shells per
    /// particle.
    pub fn new(params: DfnParams, n_n: usize, n_sep: usize, n_p: usize, m: usize) -> Self {
        assert!(n_n >= 1 && n_sep >= 1 && n_p >= 1 && m >= 1);
        let neg = &params.cell.neg;
        let pos = &params.cell.pos;
        let n = n_n + n_sep + n_p;
        let n_el = n_n + n_p;

        // Faces and cell centres (cell centres are face midpoints, so
        // neighbouring regions may have different widths).
        let dx_n = neg.thickness / n_n as f64;
        let dx_s = params.l_sep / n_sep as f64;
        let dx_p = pos.thickness / n_p as f64;
        let mut xf = Vec::with_capacity(n + 1);
        xf.push(0.0);
        for k in 1..=n_n {
            xf.push(k as f64 * dx_n);
        }
        for k in 1..=n_sep {
            xf.push(neg.thickness + k as f64 * dx_s);
        }
        for k in 1..=n_p {
            xf.push(neg.thickness + params.l_sep + k as f64 * dx_p);
        }
        let x: Vec<f64> = (0..n).map(|i| 0.5 * (xf[i] + xf[i + 1])).collect();
        let dx: Vec<f64> = (0..n).map(|i| xf[i + 1] - xf[i]).collect();
        let region: Vec<usize> = (0..n)
            .map(|i| {
                if i < n_n {
                    NEG
                } else if i < n_n + n_sep {
                    SEP
                } else {
                    POS
                }
            })
            .collect();
        let mut cell_el = vec![None; n];
        let mut el_global = Vec::with_capacity(n_el);
        let mut el_region = Vec::with_capacity(n_el);
        let mut el_local = Vec::with_capacity(n_el);
        for i in 0..n_n {
            cell_el[i] = Some(el_global.len());
            el_global.push(i);
            el_region.push(NEG);
            el_local.push(i);
        }
        for i in n_n + n_sep..n {
            cell_el[i] = Some(el_global.len());
            el_global.push(i);
            el_region.push(POS);
            el_local.push(i - (n_n + n_sep));
        }

        // Radial grids.
        let mut rvol = [Vec::new(), Vec::new()];
        let mut rarea = [Vec::new(), Vec::new()];
        let mut rdr = [0.0, 0.0];
        for (e, radius) in [(0usize, neg.radius), (1usize, pos.radius)] {
            let dr = radius / m as f64;
            rdr[e] = dr;
            rvol[e] = (0..m)
                .map(|k| {
                    let ri = k as f64 * dr;
                    let ro = (k as f64 + 1.0) * dr;
                    (4.0 / 3.0) * std::f64::consts::PI * (ro.powi(3) - ri.powi(3))
                })
                .collect();
            rarea[e] = (0..=m)
                .map(|k| 4.0 * std::f64::consts::PI * (k as f64 * dr).powi(2))
                .collect();
        }

        // State: uniform initial concentrations, zero electrolyte potential,
        // solid potential at its equilibrium value (zero overpotential).
        let c_e = vec![params.cell.c_e; n];
        let mut phi_s = Vec::with_capacity(n_el);
        let mut c_s = vec![0.0; n_el * m];
        for e in 0..n_el {
            let (el, er) = if el_region[e] == NEG {
                (&neg, 0usize)
            } else {
                (&pos, 1usize)
            };
            let x0 = el.c0 / el.c_max;
            phi_s.push(params.cell.c_e * 0.0 + 0.0); // phi_e = 0
            for k in 0..m {
                c_s[e * m + k] = el.c0;
            }
            phi_s[e] = if er == 0 {
                ocp_graphite(x0)
            } else {
                ocp_nmc(x0)
            };
        }
        let j = vec![0.0; n_el];
        let phi_e = vec![0.0; n];

        // Linear layout.
        let off_pe = 0;
        let off_ps = off_pe + (n - 1);
        let off_j = off_ps + n_el;
        let off_ce = off_j + n_el;
        let off_cs = off_ce + n;
        let ntot = off_cs + n_el * m;

        Self {
            params,
            n,
            n_el,
            x,
            dx,
            region,
            cell_el,
            el_global,
            el_region,
            el_local,
            el_ne: n_n,
            el_np: n_p,
            m,
            rvol,
            rarea,
            rdr,
            phi_e,
            phi_s,
            j,
            c_e,
            c_s,
            t: 0.0,
            off_pe,
            off_ps,
            off_j,
            off_ce,
            off_cs,
            ntot,
        }
    }

    /// Pack the state into the Newton vector.
    fn pack(&self) -> Vec<f64> {
        let mut z = vec![0.0; self.ntot];
        for i in 1..self.n {
            z[self.off_pe + i - 1] = self.phi_e[i];
        }
        for e in 0..self.n_el {
            z[self.off_ps + e] = self.phi_s[e];
            z[self.off_j + e] = self.j[e];
        }
        for i in 0..self.n {
            z[self.off_ce + i] = self.c_e[i];
        }
        for k in 0..self.n_el * self.m {
            z[self.off_cs + k] = self.c_s[k];
        }
        z
    }

    /// Unpack the Newton vector into the state.
    fn unpack(&mut self, z: &[f64]) {
        for i in 1..self.n {
            self.phi_e[i] = z[self.off_pe + i - 1];
        }
        for e in 0..self.n_el {
            self.phi_s[e] = z[self.off_ps + e];
            self.j[e] = z[self.off_j + e];
        }
        for i in 0..self.n {
            self.c_e[i] = z[self.off_ce + i];
        }
        for k in 0..self.n_el * self.m {
            self.c_s[k] = z[self.off_cs + k];
        }
    }

    /// Per-row residual scales, so the Newton convergence test is a relative
    /// one across terms with very different physical units.
    fn row_scale(&self, dt: f64, i_app: f64) -> Vec<f64> {
        let mut s = vec![1.0; self.ntot];
        let i_scale = i_app.abs().max(1e-3 * self.params.cell.i_ref()).max(1e-12);
        for k in self.off_pe..self.off_ce {
            s[k] = i_scale;
        }
        for i in 0..self.n {
            let eps = self.eps_e(i);
            s[self.off_ce + i] = (eps * self.dx[i] * self.params.cell.c_e / dt).max(1e-30);
        }
        for e in 0..self.n_el {
            let er = self.el_region[e];
            let cmax = if er == NEG {
                self.params.cell.neg.c_max
            } else {
                self.params.cell.pos.c_max
            };
            let vol = &self.rvol[if er == NEG { 0 } else { 1 }];
            for k in 0..self.m {
                s[self.off_cs + e * self.m + k] = (vol[k] * cmax / dt).max(1e-30);
            }
        }
        s
    }

    /// The full residual vector `F(z)` at time level `n+1`.
    ///
    /// * rows `off_pe..`   — electrolyte charge `∂_x i_e − a_s j = 0` (cells 1..)
    /// * rows `off_ps..`   — solid charge `∂_x(σ ∂_x φ_s) − a_s j = 0`
    /// * rows `off_j..`    — Butler–Volmer `j − j_BV(η) = 0`
    /// * rows `off_ce..`   — electrolyte mass (backward Euler)
    /// * rows `off_cs..`   — solid diffusion (backward Euler)
    fn residual(
        &self,
        z: &[f64],
        ce_old: &[f64],
        cs_old: &[f64],
        dt: f64,
        i_app: f64,
    ) -> Vec<f64> {
        let mut f = vec![0.0; self.ntot];
        // Unpack without mutating self.
        let mut phi_e = vec![0.0; self.n];
        for i in 1..self.n {
            phi_e[i] = z[self.off_pe + i - 1];
        }
        let mut phi_s = vec![0.0; self.n_el];
        let mut j = vec![0.0; self.n_el];
        for e in 0..self.n_el {
            phi_s[e] = z[self.off_ps + e];
            j[e] = z[self.off_j + e];
        }
        let c_e = &z[self.off_ce..self.off_ce + self.n];
        let c_s = &z[self.off_cs..self.off_cs + self.n_el * self.m];

        // ---- electrolyte charge on cells 1..n-1 ----
        // Cell 0's equation is replaced by the reference phi_e[0] = 0.
        for i in 1..self.n {
            let ie_r = if i + 1 < self.n {
                self.i_e_face(i + 1, &phi_e, c_e)
            } else {
                0.0
            };
            let ie_l = self.i_e_face(i, &phi_e, c_e);
            let src = match self.cell_el[i] {
                Some(e) => self.a_s(i) * j[e],
                None => 0.0,
            };
            f[self.off_pe + i - 1] = ie_r - ie_l - src * self.dx[i];
        }

        // ---- solid charge and Butler–Volmer on each electrode node ----
        for e in 0..self.n_el {
            let g = self.el_global[e];
            let er = self.el_region[e];
            let (lflux, rflux) = self.solid_face_fluxes(e, &phi_s, i_app);
            let src = self.a_s(g) * j[e];
            f[self.off_ps + e] = rflux - lflux - src * self.dx[g];

            // Butler–Volmer.
            let cs_surf = c_s[e * self.m + self.m - 1];
            let cmax = if er == NEG {
                self.params.cell.neg.c_max
            } else {
                self.params.cell.pos.c_max
            };
            let i0 = self.i0(er, cs_surf, c_e[g]);
            let eta = phi_s[e] - phi_e[g] - self.ocp(er, (cs_surf / cmax).clamp(0.0, 1.0));
            let ffe = F_FARADAY * eta / (R_GAS * self.params.cell.t);
            let (aa, ac) = if er == NEG {
                (self.params.cell.neg.alpha_a, self.params.cell.neg.alpha_c)
            } else {
                (self.params.cell.pos.alpha_a, self.params.cell.pos.alpha_c)
            };
            let jbv = i0 * (clamped_exp(aa * ffe) - clamped_exp(-ac * ffe));
            f[self.off_j + e] = j[e] - jbv;
        }

        // ---- electrolyte mass on all cells ----
        for i in 0..self.n {
            let eps = self.eps_e(i);
            let flux_r = if i + 1 < self.n {
                self.d_flux_face(i + 1, c_e)
            } else {
                0.0
            };
            let flux_l = if i > 0 { self.d_flux_face(i, c_e) } else { 0.0 };
            let src = match self.cell_el[i] {
                Some(e) => (1.0 - self.params.t_plus) * self.a_s(i) * j[e] / F_FARADAY,
                None => 0.0,
            };
            f[self.off_ce + i] = eps * self.dx[i] * (c_e[i] - ce_old[i]) / dt
                - (flux_r - flux_l)
                - src * self.dx[i];
        }

        // ---- solid diffusion on all electrode nodes ----
        for e in 0..self.n_el {
            let er = self.el_region[e];
            let (d_s, vol, area, dr) = if er == NEG {
                (
                    self.params.cell.neg.d_s,
                    &self.rvol[0],
                    &self.rarea[0],
                    self.rdr[0],
                )
            } else {
                (
                    self.params.cell.pos.d_s,
                    &self.rvol[1],
                    &self.rarea[1],
                    self.rdr[1],
                )
            };
            for k in 0..self.m {
                let idx = e * self.m + k;
                let c = c_s[idx];
                let c_old = cs_old[idx];
                let inner = if k == 0 {
                    0.0
                } else {
                    area[k] * d_s * (c_s[idx - 1] - c) / dr
                };
                let outer = if k + 1 < self.m {
                    area[k + 1] * d_s * (c_s[idx + 1] - c) / dr
                } else {
                    -area[self.m] * j[e] / F_FARADAY
                };
                f[self.off_cs + idx] =
                    vol[k] * (c - c_old) / dt - inner - outer;
            }
        }
        f
    }

    /// Solid-phase flux `σ ∂_x φ_s` at the left and right faces of electrode
    /// node `e` (positive in `+x`). The collector face carries `−i_app`; the
    /// electrode/separator face is insulating.
    fn solid_face_fluxes(
        &self,
        e: usize,
        phi_s: &[f64],
        i_app: f64,
    ) -> (f64, f64) {
        let sigma = self.params.sigma[if self.el_region[e] == NEG { 0 } else { 1 }];
        let local = self.el_local[e];
        let ne = if self.el_region[e] == NEG {
            self.el_ne
        } else {
            self.el_np
        };
        let left = if local == 0 {
            if self.el_region[e] == NEG {
                -i_app // collector at x = 0
            } else {
                0.0 // separator interface
            }
        } else {
            let g = self.el_global[e];
            let gl = self.el_global[e - 1];
            let h = self.x[g] - self.x[gl];
            sigma * (phi_s[e] - phi_s[e - 1]) / h
        };
        let right = if local + 1 == ne {
            if self.el_region[e] == POS {
                -i_app // collector at x = L
            } else {
                0.0 // negative/separator interface
            }
        } else {
            let g = self.el_global[e];
            let gr = self.el_global[e + 1];
            let h = self.x[gr] - self.x[g];
            sigma * (phi_s[e + 1] - phi_s[e]) / h
        };
        (left, right)
    }

    /// One implicit (backward Euler) step of size `dt` at applied current
    /// density `i_app` (A m⁻², discharge positive).
    pub fn step(&mut self, dt: f64, i_app: f64) {
        let ce_old = self.c_e.clone();
        let cs_old = self.c_s.clone();
        let z0 = self.pack();
        let result = self.newton(z0, &ce_old, &cs_old, dt, i_app);
        self.unpack(&result);
        self.t += dt;
    }

    /// Damped Newton solve of the coupled implicit step. Returns the updated
    /// state vector.
    fn newton(
        &self,
        z0: Vec<f64>,
        ce_old: &[f64],
        cs_old: &[f64],
        dt: f64,
        i_app: f64,
    ) -> Vec<f64> {
        let scale = self.row_scale(dt, i_app);
        let norm = |f: &[f64]| -> f64 {
            f.iter()
                .zip(&scale)
                .map(|(a, s)| (a / s).abs())
                .fold(0.0f64, f64::max)
        };
        let mut z = z0;
        let mut f = self.residual(&z, ce_old, cs_old, dt, i_app);
        let mut rnorm = norm(&f);
        let tol = 1e-9;
        for _ in 0..50 {
            if rnorm < tol {
                break;
            }
            let jac = self.fd_jacobian(&z, &f, ce_old, cs_old, dt, i_app);
            let mut rhs: Vec<f64> = f.iter().map(|v| -v).collect();
            let dz = solve_equilibrated(&jac, self.ntot, &mut rhs);
            // Backtracking line search on the scaled residual.
            let mut alpha = 1.0;
            let mut accepted = false;
            for _ in 0..16 {
                let znew: Vec<f64> = z
                    .iter()
                    .zip(&dz)
                    .map(|(a, b)| a + alpha * b)
                    .collect();
                let fnew = self.residual(&znew, ce_old, cs_old, dt, i_app);
                let nn = norm(&fnew);
                if nn.is_finite() && nn < rnorm {
                    z = znew;
                    f = fnew;
                    rnorm = nn;
                    accepted = true;
                    break;
                }
                alpha *= 0.5;
            }
            if !accepted {
                // Take a small step and stop; the outer time step will carry on.
                let znew: Vec<f64> = z
                    .iter()
                    .zip(&dz)
                    .map(|(a, b)| a + 1e-4 * b)
                    .collect();
                z = znew;
                break;
            }
        }
        z
    }

    /// Forward-difference Jacobian of the residual.
    fn fd_jacobian(
        &self,
        z: &[f64],
        f0: &[f64],
        ce_old: &[f64],
        cs_old: &[f64],
        dt: f64,
        i_app: f64,
    ) -> Vec<f64> {
        let n = self.ntot;
        let mut jac = vec![0.0; n * n];
        let mut zp = z.to_vec();
        for k in 0..n {
            let dzk = 1e-7 * (1.0 + z[k].abs());
            zp[k] = z[k] + dzk;
            let f1 = self.residual(&zp, ce_old, cs_old, dt, i_app);
            for i in 0..n {
                jac[i * n + k] = (f1[i] - f0[i]) / dzk;
            }
            zp[k] = z[k];
        }
        jac
    }

    /// Constant-current discharge from the current state to `v_cutoff`.
    ///
    /// Returns `(capacity_mAh_m2, trace)` with trace entries
    /// `(t, V, reaction_nonuniformity)`.
    pub fn discharge(
        &self,
        c_rate: f64,
        v_cutoff: f64,
        dt: f64,
        max_steps: usize,
    ) -> (f64, Vec<(f64, f64, f64)>) {
        let mut cell = self.clone();
        let i_app = c_rate * self.params.cell.i_ref();
        let mut trace = Vec::new();
        for _ in 0..max_steps {
            let v = cell.voltage();
            trace.push((cell.t, v, cell.reaction_nonuniformity()));
            if v < v_cutoff {
                break;
            }
            cell.step(dt, i_app);
        }
        (i_app * cell.t / 3600.0 * 1000.0, trace)
    }
}

/// Clamped exponential for Butler–Volmer.
#[inline]
fn clamped_exp(x: f64) -> f64 {
    x.clamp(-700.0, 700.0).exp()
}

/// Solve `A x = b` after row and column equilibration, using dense Gaussian
/// elimination with partial pivoting. `a` is row-major `n×n`; `b` is consumed.
pub fn solve_equilibrated(a: &[f64], n: usize, b: &mut [f64]) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    // Row/column max norms.
    let mut row = vec![0.0f64; n];
    let mut col = vec![0.0f64; n];
    for i in 0..n {
        for k in 0..n {
            let v = a[i * n + k].abs();
            if v > row[i] {
                row[i] = v;
            }
            if v > col[k] {
                col[k] = v;
            }
        }
    }
    for v in row.iter_mut() {
        if *v <= 0.0 {
            *v = 1.0;
        }
    }
    for v in col.iter_mut() {
        if *v <= 0.0 {
            *v = 1.0;
        }
    }
    let mut m = vec![0.0f64; n * n];
    for i in 0..n {
        for k in 0..n {
            m[i * n + k] = a[i * n + k] / (row[i] * col[k]);
        }
    }
    let mut rhs: Vec<f64> = (0..n).map(|i| b[i] / row[i]).collect();

    // LU with partial pivoting (Doolittle, in place).
    for c in 0..n {
        let mut piv = c;
        let mut mx = m[c * n + c].abs();
        for r in c + 1..n {
            let v = m[r * n + c].abs();
            if v > mx {
                mx = v;
                piv = r;
            }
        }
        if mx < 1e-300 {
            continue;
        }
        if piv != c {
            for k in 0..n {
                m.swap(piv * n + k, c * n + k);
            }
            rhs.swap(piv, c);
        }
        let d = m[c * n + c];
        for r in c + 1..n {
            let factor = m[r * n + c] / d;
            if factor == 0.0 {
                continue;
            }
            m[r * n + c] = 0.0;
            for k in c + 1..n {
                m[r * n + k] -= factor * m[c * n + k];
            }
            rhs[r] -= factor * rhs[c];
        }
    }
    // Back substitution.
    let mut y = vec![0.0f64; n];
    for i in (0..n).rev() {
        let mut s = rhs[i];
        for k in i + 1..n {
            s -= m[i * n + k] * y[k];
        }
        let d = m[i * n + i];
        y[i] = if d.abs() > 1e-300 { s / d } else { 0.0 };
    }
    // Undo the column scaling: x_k = y_k / col_k.
    (0..n).map(|k| y[k] / col[k]).collect()
}
