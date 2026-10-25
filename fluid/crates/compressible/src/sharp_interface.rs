//! Sharp solid–fluid interface in 1D, coupled by a real-ghost-fluid /
//! exact-Riemann interface flux.
//!
//! The crate's multi-material model is the *diffuse-interface* five-equation
//! model ([`crate::multiphase::FiveEq`]); at a water/steel contrast it mixes
//! the two equations of state over a few cells and the transmitted amplitude
//! carries an O(10%) error at a 1024-cell run, with `|R| > 1` (see
//! `analysis/fsi.md`). This module is the **unfolded** rung: each cell holds
//! one pure material, the material interface is tracked (a level set), and
//! the face flux on the interface is the *exact* two-material Riemann
//! solution. The interface is one cell face wide, so there is no EOS mixing.
//!
//! # Materials
//!
//! * A **fluid** is a stiffened gas (NASG, `b = q = 0`), with normal stress
//!   equal to the pressure, `τ = p`.
//! * A **solid** is the 1D reduction of the GPR unified-continuum model
//!   (`crate::gpr::Gpr`). In 1D the transverse distortion is constant and the
//!   longitudinal stretch is slaved to the density by `det A = ρ/ρ₀`
//!   (the A-equation contracts to `D(a)/Dt = −a ∂ₓu`), so the solid is a
//!   hyperelastic material with
//!   ```text
//!     a      = ρ/ρ₀,
//!     E₂(ρ)  = (c_s²/6)(a² − 1)²,          (GPR elastic energy, ρ0-fixed)
//!     f(ρ)   = (2/3) c_s² ρ a²(a² − 1),    (elastic normal-stress part)
//!     τ      = p_h + f(ρ),   e_tot = e_h(ρ,p_h) + E₂(ρ),
//!     c_L²   = γ(p_h+p∞)/ρ + f'(ρ),
//!   ```
//!   which is exactly `τ = ρ² ∂e_tot/∂ρ|_s` and gives the GPR longitudinal
//!   speed `c_L = sqrt(c₀² + 4c_s²/3)` in the linear limit. Shear decouples
//!   at normal incidence, so this is the whole longitudinal response.
//!
//! # Exact solid–fluid Riemann problem
//!
//! At a material interface the normal velocity and the normal traction are
//! continuous and the inviscid fluid carries no shear traction, so the
//! interface state `(u*, τ*)` solves
//! `u_L + F_L(τ*) = u* = u_R − F_R(τ*)`, where `F(τ*)` is the velocity jump
//! across the wave connecting a side to the interface. For a rarefaction
//! `F = ∫ c dρ/ρ` along the isentrope; for a shock `F = (τ*−τ)/J` with the
//! mass flux `J = sqrt((τ*−τ)/(v−v*))` from the Rankine–Hugoniot relations.
//! Both are evaluated with the material EOS above, so a stiffened gas and an
//! elastic solid are handled by the same code. The exact solver is used as
//! the interface face flux and independently as the validation reference.
//!
//! Two conservative discretisations are provided: an Eulerian level-set
//! [`SharpSolver`] and a Lagrangian moving-mesh [`LagrangianSolver`] in which
//! the interface is a mesh node. Both use one single-valued flux per face.
//!
//! # What is conserved exactly
//!
//! The update is `dU_i/dt = −(F_{i+½} − F_{i−½})/Δx` with **one** flux per
//! face (the exact interface flux at an interface face, HLLC elsewhere). The
//! flux telescopes, so total mass, momentum and energy over the domain are
//! conserved to machine precision for any level-set motion, because the
//! conserved variables are never relabelled. The solid's elastic energy is
//! part of the conserved total energy through the EOS, so no separate strain
//! ledger is needed.
//!
//! # References
//!
//! * Wang, Liu & Khoo, *A real ghost fluid method for the simulation of
//!   compressible multifluid flow*, SIAM J. Sci. Comput. 28 (2006).
//! * Ndanou, Favrie & Gavrilyuk, *Multi-solid and multi-fluid diffuse
//!   interface Godunov method*, J. Comput. Phys. 295 (2015) (exact
//!   solid–fluid Riemann structure).
//! * Peshkov & Romenski, *A hyperbolic model for viscous Newtonian flows*,
//!   Continuum Mech. Thermodyn. 28 (2016) (GPR EOS and `c_L`).
//! * Toro, *Riemann Solvers and Numerical Methods for Fluid Dynamics*
//!   (shock/rarefaction functions).

use crate::eos::Nasg;
use crate::scheme::weno5z;

/// Atmospheric reference pressure [Pa].
pub const P0: f64 = 101_325.0;
/// Stiffened-gas exponent used to fit the equation of state at `P0`; it only
/// controls the (weak) nonlinearity of the reference fit. Water uses the
/// usual 4.4 shock-fitting value (Le Métayer & Saurel 2016 give 1.19 for NASG
/// liquid water; 4.4 is the coarser stiffened-gas fit used throughout this
/// crate, `crate::eos::materials::WATER_SG`).
pub const GAMMA: f64 = 4.4;

/// Stiffened-gas NASG fitted to `(ρ, c)` at pressure `p0` exactly for `b=0`:
/// `c² = γ(p0 + p∞)/ρ`, so `p∞ = ρ c²/γ − p0`.
pub fn nasg_at(rho: f64, c: f64, p0: f64, gamma: f64) -> Nasg {
    Nasg::stiffened(gamma, rho * c * c / gamma - p0)
}

/// What kind of material a cell holds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Fluid,
    Solid,
}

/// A pure 1D material: stiffened-gas hydrostatic EOS plus optional elastic
/// normal-stress stiffening (GPR 1D reduction).
#[derive(Clone, Copy, Debug)]
pub struct Material {
    pub kind: Kind,
    /// hydrostatic (NASG) equation of state
    pub eos: Nasg,
    /// reference density ρ₀ (elastic reference state)
    pub rho0: f64,
    /// GPR shear sound speed c_s (0 for a fluid)
    pub cs: f64,
}

impl Material {
    /// A stiffened-gas fluid with density `rho` and sound speed `c` at `P0`.
    pub fn fluid(rho: f64, c: f64) -> Self {
        Self { kind: Kind::Fluid, eos: nasg_at(rho, c, P0, GAMMA), rho0: rho, cs: 0.0 }
    }

    /// A GPR elastic solid with density `rho`, P-wave speed `c_p` and shear
    /// speed `c_s`. The NASG bulk speed is set to `c₀ = sqrt(c_p² − 4c_s²/3)`
    /// so the longitudinal speed at the reference state is exactly `c_p`.
    pub fn solid(rho: f64, c_p: f64, c_s: f64) -> Self {
        let c0sq = (c_p * c_p - 4.0 / 3.0 * c_s * c_s).max(1.0);
        Self { kind: Kind::Solid, eos: nasg_at(rho, c0sq.sqrt(), P0, GAMMA), rho0: rho, cs: c_s }
    }

    /// Longitudinal linear speed at the reference state: `sqrt(c₀² + 4c_s²/3)`.
    pub fn c_longitudinal(&self) -> f64 {
        (self.eos.gamma * (P0 + self.eos.p_inf) / self.rho0 + 4.0 / 3.0 * self.cs * self.cs).sqrt()
    }

    #[inline]
    fn stretch(&self, rho: f64) -> f64 {
        rho / self.rho0
    }

    /// `f(ρ) = (2/3)c_s² ρ a²(a²−1)`, the elastic part of the normal stress.
    #[inline]
    pub fn elastic_stress(&self, rho: f64) -> f64 {
        if self.kind == Kind::Fluid {
            return 0.0;
        }
        let a = self.stretch(rho);
        (2.0 / 3.0) * self.cs * self.cs * rho * a * a * (a * a - 1.0)
    }

    /// `f'(ρ) = (2/3)c_s²(5a⁴ − 3a²)`, the derivative of the elastic stress.
    #[inline]
    pub fn elastic_stress_deriv(&self, rho: f64) -> f64 {
        if self.kind == Kind::Fluid {
            return 0.0;
        }
        let a = self.stretch(rho);
        (2.0 / 3.0) * self.cs * self.cs * (5.0 * a * a * a * a - 3.0 * a * a)
    }

    /// GPR elastic energy per unit mass `E₂(ρ) = (c_s²/6)(a²−1)²`.
    #[inline]
    pub fn elastic_energy(&self, rho: f64) -> f64 {
        if self.kind == Kind::Fluid {
            return 0.0;
        }
        let a = self.stretch(rho);
        self.cs * self.cs / 6.0 * (a * a - 1.0) * (a * a - 1.0)
    }

    /// Hydrostatic pressure from density and *total* specific internal energy.
    #[inline]
    pub fn p_hydro(&self, rho: f64, e_tot: f64) -> f64 {
        self.eos.pressure(rho, e_tot - self.elastic_energy(rho))
    }

    /// Total normal traction `τ = p_h + f(ρ)`.
    #[inline]
    pub fn tau(&self, rho: f64, p_h: f64) -> f64 {
        p_h + self.elastic_stress(rho)
    }

    /// Total specific internal energy.
    #[inline]
    pub fn energy(&self, rho: f64, p_h: f64) -> f64 {
        self.eos.energy(rho, p_h) + self.elastic_energy(rho)
    }

    /// Isentropic longitudinal sound speed `c² = γ(p_h+p∞)/ρ + f'(ρ)`.
    #[inline]
    pub fn sound_speed(&self, rho: f64, p_h: f64) -> f64 {
        let bulk = self.eos.sound_speed_sq(rho, p_h);
        (bulk + self.elastic_stress_deriv(rho)).max(0.0).sqrt()
    }

    /// A finite, positive-density, positive-(p+p∞) state.
    #[inline]
    pub fn admissible(&self, w: &State) -> bool {
        w.rho.is_finite()
            && w.u.is_finite()
            && w.p.is_finite()
            && w.rho > 0.0
            && w.p + self.eos.p_inf > 0.0
    }

    /// Isentropic sound speed along the isentrope with constant
    /// `A = (p_h+p∞) ρ^{−γ}`.
    #[inline]
    fn sound_speed_isentropic(&self, a0: f64, rho: f64) -> f64 {
        let bulk = a0 * self.eos.gamma * rho.powf(self.eos.gamma - 1.0);
        (bulk + self.elastic_stress_deriv(rho)).max(0.0).sqrt()
    }

    /// Isentrope constant `A = (p_h+p∞) ρ^{−γ}` of a state.
    #[inline]
    fn isentrope_constant(&self, rho: f64, p_h: f64) -> f64 {
        (p_h + self.eos.p_inf) * rho.powf(-self.eos.gamma)
    }
}

/// Directional state at a material interface: density, velocity, hydrostatic
/// pressure and the material.
#[derive(Clone, Copy, Debug)]
pub struct State {
    pub rho: f64,
    pub u: f64,
    pub p: f64,
}

impl Material {
    /// Conserved state `[ρ, ρu, E]` of a primitive state.
    pub fn to_cons(&self, w: State) -> Cons {
        [w.rho, w.rho * w.u, w.rho * (self.energy(w.rho, w.p) + 0.5 * w.u * w.u)]
    }

    /// Primitive hydrostatic state of a conserved state.
    pub fn to_prim(&self, c: Cons) -> State {
        let rho = c[0].max(1e-300);
        let u = c[1] / rho;
        let e_tot = c[2] / rho - 0.5 * u * u;
        State { rho, u, p: self.p_hydro(rho, e_tot) }
    }
}

/// Interface state of an exact two-material Riemann problem.
#[derive(Clone, Copy, Debug)]
pub struct Star {
    /// common interface normal traction
    pub tau: f64,
    /// common interface normal velocity
    pub u: f64,
    /// density on the left of the contact
    pub rho_l: f64,
    /// density on the right of the contact
    pub rho_r: f64,
}

/// Velocity jump `F(τ*)` across the wave joining a side to the interface,
/// its star density, and `dF/dτ*`.
fn wave_jump(m: &Material, rho0: f64, p0: f64, tau_star: f64) -> (f64, f64, f64) {
    let tau0 = m.tau(rho0, p0);
    if tau_star <= tau0 {
        // isentropic rarefaction: solve a0 ρ^γ − p∞ + f(ρ) = τ*.
        let a0 = m.isentrope_constant(rho0, p0);
        let rho_s = solve_isentrope(m, a0, rho0, tau_star);
        let du = integrate_sound_speed(m, a0, rho0, rho_s);
        let c_s = m.sound_speed_isentropic(a0, rho_s);
        (du, rho_s, 1.0 / (rho_s * c_s).max(1e-300))
    } else {
        // shock: Rankine–Hugoniot determines the star density.
        let rho_s = solve_hugoniot(m, rho0, p0, tau_star);
        let (v0, vs) = (1.0 / rho0, 1.0 / rho_s);
        let j = ((tau_star - tau0) / (v0 - vs).max(1e-300)).max(0.0).sqrt();
        ((tau_star - tau0) / j.max(1e-300), rho_s, 1.0 / j.max(1e-300))
    }
}

/// Solve `a0 ρ^γ − p∞ + f(ρ) = τ*` for ρ on the isentrope (monotone
/// increasing since `c² > 0`). For the rarefaction branch `τ* < τ0` the root
/// lies below `ρ0`, so a fixed bracket `[ρ0·1e−10, ρ0]` is valid; a
/// safeguarded Newton iteration (derivative `dτ/dρ = c²`) converges in a few
/// steps for the acoustic case.
fn solve_isentrope(m: &Material, a0: f64, rho0: f64, tau_star: f64) -> f64 {
    let g = |rho: f64| a0 * rho.powf(m.eos.gamma) - m.eos.p_inf + m.elastic_stress(rho) - tau_star;
    // dg/drho = d(tau)/drho = c(rho)^2 along the isentrope.
    let dg = |rho: f64| m.sound_speed_isentropic(a0, rho).powi(2).max(1e-300);
    let mut lo = rho0 * 1e-10;
    let mut hi = rho0;
    // g(lo) < 0 < g(hi) for tau_star < tau0; safeguarded Newton on a
    // strictly increasing g (the acoustic case converges in two steps).
    let mut rho = rho0;
    for _ in 0..50 {
        let gv = g(rho);
        if gv > 0.0 {
            hi = rho;
        } else {
            lo = rho;
        }
        let mut next = rho - gv / dg(rho);
        if !next.is_finite() || next <= lo || next >= hi {
            next = 0.5 * (lo + hi);
        }
        if (next - rho).abs() < 1e-15 * (1.0 + rho.abs()) {
            rho = next;
            break;
        }
        rho = next;
    }
    rho
}

/// `∫_{ρ0}^{ρ1} (c(ρ)/ρ) dρ` along the isentrope, composite 8-point
/// Gauss–Legendre with the number of panels set by how many decades the
/// rarefaction spans (a near-vacuum fan gets 64 panels, an acoustic one four),
/// so a near-vacuum rarefaction is integrated to well below 1e-6 relative
/// while the small-amplitude solve stays cheap. This is `−∫ dτ/(ρc)` (since
/// `dτ = ρc² dρ`), the velocity jump across a rarefaction.
fn integrate_sound_speed(m: &Material, a0: f64, rho0: f64, rho1: f64) -> f64 {
    // 8-point Gauss–Legendre on [-1, 1]
    const X: [f64; 8] = [
        -0.9602898564975363,
        -0.7966664774136267,
        -0.5255324099163290,
        -0.1834346424956498,
        0.1834346424956498,
        0.5255324099163290,
        0.7966664774136267,
        0.9602898564975363,
    ];
    const W: [f64; 8] = [
        0.1012285362903763,
        0.2223810344533745,
        0.3137066458778873,
        0.3626837833783620,
        0.3626837833783620,
        0.3137066458778873,
        0.2223810344533745,
        0.1012285362903763,
    ];
    let half = 0.5 * (rho1 - rho0);
    let mid = 0.5 * (rho1 + rho0);
    // Composite Gauss–Legendre 8-point rule. The integrand c(rho)/rho varies
    // on the scale of log(rho), so spend panels in proportion to the number
    // of decades crossed: a near-vacuum fan gets up to 64 panels, an acoustic
    // one only a handful. This is what makes the small-amplitude interface
    // Riemann solve cheap without degrading the exact reference.
    let ratio = (rho1 / rho0).abs().max(1e-300);
    let decades = ratio.ln().abs() / std::f64::consts::LN_10;
    let sub = (4 + (32.0 * decades) as usize).clamp(4, 64);
    let ph = half / sub as f64;
    let mut s = 0.0;
    for p in 0..sub {
        let c = mid + (2.0 * p as f64 + 1.0 - sub as f64) * ph;
        for k in 0..8 {
            let rho = c + ph * X[k];
            s += W[k] * m.sound_speed_isentropic(a0, rho) / rho;
        }
    }
    ph * s
}

/// Solve the Hugoniot
/// `e_tot(ρ*) − e₀ − ½(τ*+τ₀)(1/ρ₀ − 1/ρ*) = 0` for the post-shock density.
/// `H(ρ₀) > 0` and `H` decreases monotonically with compression, so the
/// bracket is `[ρ₀, ρ_max]` with `ρ_max` expanded until `H < 0`.
fn solve_hugoniot(m: &Material, rho0: f64, p0: f64, tau_star: f64) -> f64 {
    let tau0 = m.tau(rho0, p0);
    let e0 = m.energy(rho0, p0);
    let h = |rho: f64| {
        let p_h = tau_star - m.elastic_stress(rho);
        let e = m.eos.energy(rho, p_h) + m.elastic_energy(rho);
        e - e0 - 0.5 * (tau_star + tau0) * (1.0 / rho0 - 1.0 / rho)
    };
    let mut hi = rho0 * 1.01;
    let mut it = 0;
    while h(hi) > 0.0 && it < 100 {
        hi = rho0 + (hi - rho0) * 2.0;
        it += 1;
    }
    let mut lo = rho0;
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if h(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Exact two-material Riemann solver: solve
/// `u_L − F_L(τ*) = u_R + F_R(τ*)` for the interface traction `τ*`, where
/// `F_K` is the velocity jump across the wave (positive for a shock, negative
/// for a rarefaction) with the sign appropriate to the wave direction.
pub fn exact_star(ml: &Material, wl: &State, mr: &Material, wr: &State) -> Star {
    let tau_l = ml.tau(wl.rho, wl.p);
    let tau_r = mr.tau(wr.rho, wr.p);
    let fl = |tau: f64| wave_jump(ml, wl.rho, wl.p, tau).0;
    let fr = |tau: f64| wave_jump(mr, wr.rho, wr.p, tau).0;
    // g is decreasing in tau (F_K is increasing).
    let g = |tau: f64| wl.u - fl(tau) - (wr.u + fr(tau));

    // acoustic initial guess: F_K ≈ (τ−τ_K)/Z_K
    let zl = wl.rho * ml.sound_speed(wl.rho, wl.p);
    let zr = wr.rho * mr.sound_speed(wr.rho, wr.p);

    // The traction cannot fall below the vacuum limit τ > −min(p∞).
    let tau_min = -ml.eos.p_inf.min(mr.eos.p_inf)
        + 1e-9 * (1.0 + tau_l.abs() + tau_r.abs());
    let tau_ac = ((wl.u - wr.u + tau_l / zl + tau_r / zr) / (1.0 / zl + 1.0 / zr)).max(tau_min);

    // bracket by expansion: below the root g > 0, above it g < 0
    let mut lo = tau_ac;
    let mut hi = tau_ac;
    let step0 = 0.25 * (zl + zr) * (wl.u.abs() + wr.u.abs() + 100.0) + 1e3;
    let mut step = step0;
    let mut n = 0;
    while g(lo) < 0.0 && lo > tau_min && n < 200 {
        lo = (lo - step).max(tau_min);
        step *= 1.4;
        n += 1;
    }
    step = step0;
    n = 0;
    while g(hi) > 0.0 && n < 200 {
        hi += step;
        step *= 1.4;
        n += 1;
    }
    // Bisection on the traction: `g` is strictly decreasing, and 55
    // iterations take the bracket to well below 1e-9 relative. The derivative
    // in `wave_jump` is only exact on the rarefaction branch, so Newton can
    // stagnate on strong shocks; plain bisection is robust and, after the
    // adaptive rarefaction quadrature, cheap.
    for _ in 0..55 {
        let mid = 0.5 * (lo + hi);
        let gv = wl.u - fl(mid) - (wr.u + fr(mid));
        if gv > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let tau = 0.5 * (lo + hi);
    let u = wl.u - fl(tau);
    let rho_l = wave_jump(ml, wl.rho, wl.p, tau).1;
    let rho_r = wave_jump(mr, wr.rho, wr.p, tau).1;
    Star { tau, u, rho_l, rho_r }
}

/// Acoustic (linearised) two-material interface state. Exact for small
/// amplitudes; the limit of [`exact_star`] as the wave strength goes to zero.
pub fn acoustic_star(ml: &Material, wl: &State, mr: &Material, wr: &State) -> Star {
    let zl = wl.rho * ml.sound_speed(wl.rho, wl.p);
    let zr = wr.rho * mr.sound_speed(wr.rho, wr.p);
    let tau_l = ml.tau(wl.rho, wl.p);
    let tau_r = mr.tau(wr.rho, wr.p);
    let tau = (wl.u - wr.u + tau_l / zl + tau_r / zr) / (1.0 / zl + 1.0 / zr);
    let u = wl.u - (tau - tau_l) / zl;
    Star { tau, u, rho_l: wl.rho, rho_r: wr.rho }
}

/// State at similarity coordinate `ξ = x/t` of the exact two-material
/// Riemann solution, for comparison with the finite-volume run.
pub fn sample_exact(ml: &Material, wl: &State, mr: &Material, wr: &State, star: &Star, xi: f64) -> State {
    let left = xi <= star.u;
    let (m, w, rho_s) = if left { (ml, wl, star.rho_l) } else { (mr, wr, star.rho_r) };
    let tau0 = m.tau(w.rho, w.p);
    let star_state = |rho_s: f64| {
        let p = star.tau - m.elastic_stress(rho_s);
        State { rho: rho_s, u: star.u, p }
    };
    if star.tau > tau0 {
        // shock: initial state is unswept outside the shock position
        let rho_shock = solve_hugoniot(m, w.rho, w.p, star.tau);
        let j = ((star.tau - tau0) / (1.0 / w.rho - 1.0 / rho_shock)).max(0.0).sqrt();
        let s = if left { w.u - j / w.rho } else { w.u + j / w.rho };
        let unswept = if left { xi < s } else { xi > s };
        if unswept {
            *w
        } else {
            star_state(rho_s)
        }
    } else {
        // rarefaction: fan bounded by the head and the star-side tail
        let a0 = m.isentrope_constant(w.rho, w.p);
        let c0 = m.sound_speed(w.rho, w.p);
        let c_star = m.sound_speed_isentropic(a0, rho_s);
        let head = if left { w.u - c0 } else { w.u + c0 };
        let tail = if left { star.u - c_star } else { star.u + c_star };
        let outside = if left { xi < head } else { xi > head };
        let beyond = if left { xi > tail } else { xi < tail };
        if outside {
            return *w;
        }
        if beyond {
            return star_state(rho_s);
        }
        // solve h(τ) = ξ for the local traction, h = u(τ) ∓ c(τ)
        let h = |tau: f64| -> f64 {
            let (du, rho, _) = wave_jump(m, w.rho, w.p, tau);
            let c = m.sound_speed_isentropic(a0, rho);
            if left {
                w.u - du - c
            } else {
                w.u + du + c
            }
        };
        let (mut tau_lo, mut tau_hi) = (star.tau, tau0);
        // right side: h decreases in τ; left side: h increases in τ
        let decreasing = !left;
        for _ in 0..80 {
            let mid = 0.5 * (tau_lo + tau_hi);
            let v = h(mid);
            if decreasing {
                if v > xi {
                    tau_lo = mid;
                } else {
                    tau_hi = mid;
                }
            } else if v < xi {
                tau_lo = mid;
            } else {
                tau_hi = mid;
            }
        }
        let tau = 0.5 * (tau_lo + tau_hi);
        let (du, rho, _) = wave_jump(m, w.rho, w.p, tau);
        let u = if left { w.u - du } else { w.u + du };
        State { rho, u, p: tau - m.elastic_stress(rho) }
    }
}

// =====================================================================
// Conservative 1D finite-volume solver with a tracked material interface
// =====================================================================

/// Per-cell conserved state `[ρ, ρu, E]`.
pub type Cons = [f64; 3];

/// A 1D sharp-interface solver on a uniform fixed grid.
///
/// The material of each cell is given by the sign of a level set `φ`
/// (`φ > 0` is solid material `mats[1]`, `φ ≤ 0` is the fluid `mats[0]`).
/// The level set is advected with the local velocity; the interface is the
/// zero crossing.
pub struct SharpSolver {
    pub n: usize,
    pub ng: usize,
    pub lo: f64,
    pub dx: f64,
    /// conserved state, ghost-padded (`n + 2 ng` cells)
    pub u: Vec<Cons>,
    /// signed level set at interior cell centres
    pub phi: Vec<f64>,
    /// material id per cell, ghost-padded
    pub mat: Vec<usize>,
    pub mats: [Material; 2],
    pub cfl: f64,
    pub t: f64,
}

impl SharpSolver {
    /// Build from an initial conservation state and a level set. `phi > 0`
    /// selects `solid`.
    pub fn new(n: usize, ng: usize, lo: f64, hi: f64, mats: [Material; 2], init: impl Fn(f64) -> (Cons, f64)) -> Self {
        let dx = (hi - lo) / n as f64;
        let mut u = vec![[0.0; 3]; n + 2 * ng];
        let mut phi = vec![0.0; n];
        let mut mat = vec![0usize; n + 2 * ng];
        for i in 0..n {
            let x = lo + (i as f64 + 0.5) * dx;
            let (c, p) = init(x);
            u[i + ng] = c;
            phi[i] = p;
            mat[i + ng] = if p > 0.0 { 1 } else { 0 };
        }
        for g in 0..ng {
            mat[g] = mat[ng];
            mat[n + ng + g] = mat[n + ng - 1];
        }
        SharpSolver { n, ng, lo, dx, u, phi, mat, mats, cfl: 0.4, t: 0.0 }
    }

    #[inline]
    fn mat_at(&self, i: isize) -> usize {
        let idx = i.clamp(0, (self.n - 1) as isize) as usize;
        self.mat[idx + self.ng]
    }

    /// Primitive `(ρ, u, p_h)` of cell `i` (interior index).
    #[inline]
    pub fn prim(&self, i: usize) -> State {
        self.prim_from(&self.u[i + self.ng], self.mat[i + self.ng])
    }

    #[inline]
    fn prim_from(&self, c: &Cons, m: usize) -> State {
        let rho = c[0].max(1e-300);
        let u = c[1] / rho;
        let e_tot = c[2] / rho - 0.5 * u * u;
        let mat = &self.mats[m];
        let p = mat.p_hydro(rho, e_tot);
        State { rho, u, p }
    }

    #[inline]
    fn cons_from(&self, w: &State, m: usize) -> Cons {
        let mat = &self.mats[m];
        [w.rho, w.rho * w.u, w.rho * (mat.energy(w.rho, w.p) + 0.5 * w.u * w.u)]
    }

    /// Time step from the CFL condition using the local longitudinal speed.
    pub fn dt_cfl(&self) -> f64 {
        let mut rate = 0.0f64;
        for i in 0..self.n {
            let w = self.prim(i);
            let m = &self.mats[self.mat[i + self.ng]];
            let c = m.sound_speed(w.rho, w.p);
            rate = rate.max((c + w.u.abs()) / self.dx);
        }
        self.cfl / rate.max(1e-300)
    }

    /// Fill ghost cells by zero-gradient extrapolation.
    fn fill_ghosts(&mut self) {
        for g in 0..self.ng {
            let lo_src = self.u[self.ng];
            let hi_src = self.u[self.ng + self.n - 1];
            self.u[g] = lo_src;
            self.u[self.ng + self.n + g] = hi_src;
            self.mat[g] = self.mat[self.ng];
            self.mat[self.ng + self.n + g] = self.mat[self.ng + self.n - 1];
        }
    }

    /// Physical flux of a pure material state.
    #[inline]
    fn physical_flux(&self, w: &State, m: usize) -> Cons {
        let mat = &self.mats[m];
        let tau = mat.tau(w.rho, w.p);
        let e = mat.energy(w.rho, w.p) + 0.5 * w.u * w.u;
        [w.rho * w.u, w.rho * w.u * w.u + tau, (w.rho * e + tau) * w.u]
    }

    /// HLLC flux for a same-material face (Toro). `τ` plays the role of `p`.
    fn hllc_flux(&self, wl: &State, wr: &State, m: usize) -> Cons {
        let mat = &self.mats[m];
        let cl = mat.sound_speed(wl.rho, wl.p);
        let cr = mat.sound_speed(wr.rho, wr.p);
        let sl = (wl.u - cl).min(wr.u - cr);
        let sr = (wl.u + cl).max(wr.u + cr);
        if sl >= 0.0 {
            return self.physical_flux(wl, m);
        }
        if sr <= 0.0 {
            return self.physical_flux(wr, m);
        }
        let (taul, taur) = (mat.tau(wl.rho, wl.p), mat.tau(wr.rho, wr.p));
        let den = wl.rho * (sl - wl.u) - wr.rho * (sr - wr.u);
        let ss = (taur - taul + wl.rho * wl.u * (sl - wl.u) - wr.rho * wr.u * (sr - wr.u)) / den;
        let star = |w: &State, tau: f64, s: f64| -> (f64, Cons) {
            let fac = w.rho * (s - w.u) / (s - ss);
            let e = mat.energy(w.rho, w.p) + 0.5 * w.u * w.u;
            let ustar = fac * ss;
            let estar = fac * (e + (ss - w.u) * (ss + tau / (w.rho * (s - w.u))));
            let c = [fac, ustar, estar];
            (fac, c)
        };
        if ss >= 0.0 {
            let (_, cstar) = star(wl, taul, sl);
            let f = self.physical_flux(wl, m);
            std::array::from_fn(|k| f[k] + sl * (cstar[k] - self.cons_from(wl, m)[k]))
        } else {
            let (_, cstar) = star(wr, taur, sr);
            let f = self.physical_flux(wr, m);
            std::array::from_fn(|k| f[k] + sr * (cstar[k] - self.cons_from(wr, m)[k]))
        }
    }

    /// Exact two-material interface flux. The conserved flux is taken from
    /// the star state on the upwind side of the contact, so the face flux is
    /// single-valued and the update conserves exactly.
    fn interface_flux(&self, wl: &State, ml: usize, wr: &State, mr: usize) -> Cons {
        let st = exact_star(&self.mats[ml], wl, &self.mats[mr], wr);
        let (rho_s, m) = if st.u >= 0.0 { (st.rho_l, ml) } else { (st.rho_r, mr) };
        let mat = &self.mats[m];
        let p_h = st.tau - mat.elastic_stress(rho_s);
        let e = mat.energy(rho_s, p_h) + 0.5 * st.u * st.u;
        [rho_s * st.u, rho_s * st.u * st.u + st.tau, (rho_s * e + st.tau) * st.u]
    }

    /// Right-hand side `−(F_{i+½} − F_{i−½})/Δx` for the interior cells.
    fn rhs(&self, u: &[Cons]) -> Vec<Cons> {
        // Work on a copy whose ghosts are refreshed for this stage: the RK
        // stage states are otherwise left with the previous stage's ghosts,
        // which breaks the boundary flux balance.
        let mut u = u.to_vec();
        for g in 0..self.ng {
            u[g] = u[self.ng];
            u[self.ng + self.n + g] = u[self.ng + self.n - 1];
        }
        let mut out = vec![[0.0; 3]; self.n];
        // build primitive states from the passed conserved array
        let prim_of = |c: &Cons, m: usize| -> State {
            let rho = c[0].max(1e-300);
            let uu = c[1] / rho;
            let e_tot = c[2] / rho - 0.5 * uu * uu;
            State { rho, u: uu, p: self.mats[m].p_hydro(rho, e_tot) }
        };
        let prim = |c: isize| -> State {
            let idx = (c + self.ng as isize).clamp(0, (u.len() - 1) as isize) as usize;
            prim_of(&u[idx], self.mat_at(c))
        };
        let flux = |f: usize| -> Cons {
            // face f separates cells f-1 and f, so cell i has faces i and i+1
            let i = f as isize - 1;
            let ml = self.mat_at(i);
            let mr = self.mat_at(i + 1);
            let all = |m: usize, c: isize| (-2..=2).all(|k| self.mat_at(c + k) == m);
            let wl = if all(ml, i) {
                let (p0, p1, p2, p3, p4) = (prim(i - 2), prim(i - 1), prim(i), prim(i + 1), prim(i + 2));
                State { rho: weno5z(p0.rho, p1.rho, p2.rho, p3.rho, p4.rho), u: weno5z(p0.u, p1.u, p2.u, p3.u, p4.u), p: weno5z(p0.p, p1.p, p2.p, p3.p, p4.p) }
            } else {
                prim(i)
            };
            let wr = if all(mr, i + 1) {
                let (p0, p1, p2, p3, p4) = (prim(i + 3), prim(i + 2), prim(i + 1), prim(i), prim(i - 1));
                State { rho: weno5z(p0.rho, p1.rho, p2.rho, p3.rho, p4.rho), u: weno5z(p0.u, p1.u, p2.u, p3.u, p4.u), p: weno5z(p0.p, p1.p, p2.p, p3.p, p4.p) }
            } else {
                prim(i + 1)
            };
            // A WENO overshoot must not feed an inadmissible state into the
            // Riemann solver; fall back to the cell value there.
            let wl = if self.mats[ml].admissible(&wl) { wl } else { prim(i) };
            let wr = if self.mats[mr].admissible(&wr) { wr } else { prim(i + 1) };
            if ml == mr {
                self.hllc_flux(&wl, &wr, ml)
            } else {
                self.interface_flux(&wl, ml, &wr, mr)
            }
        };
        let fl: Vec<Cons> = (0..=self.n).map(flux).collect();
        for i in 0..self.n {
            for k in 0..3 {
                out[i][k] = -(fl[i + 1][k] - fl[i][k]) / self.dx;
            }
        }
        out
    }

    /// One SSP-RK3 step of size `dt` (advecting the level set once).
    pub fn step(&mut self, dt: f64) {
        self.fill_ghosts();
        let l1 = self.rhs(&self.u);
        let mut u1 = self.u.clone();
        for i in self.ng..self.ng + self.n {
            for k in 0..3 {
                u1[i][k] = self.u[i][k] + dt * l1[i - self.ng][k];
            }
        }
        let l2 = self.rhs(&u1);
        let mut u2 = self.u.clone();
        for i in self.ng..self.ng + self.n {
            for k in 0..3 {
                u2[i][k] = 0.75 * self.u[i][k] + 0.25 * (u1[i][k] + dt * l2[i - self.ng][k]);
            }
        }
        let l3 = self.rhs(&u2);
        for i in self.ng..self.ng + self.n {
            for k in 0..3 {
                self.u[i][k] = self.u[i][k] / 3.0 + 2.0 / 3.0 * (u2[i][k] + dt * l3[i - self.ng][k]);
            }
        }
        self.advect_level_set(dt);
        self.t += dt;
    }

    /// Advect the level set with the local velocity and re-tag the materials.
    fn advect_level_set(&mut self, dt: f64) {
        let n = self.n;
        let dx = self.dx;
        let mut dphi = vec![0.0; n];
        for i in 0..n {
            let u = self.u[i + self.ng][1] / self.u[i + self.ng][0];
            let (lo, hi) = (self.phi[i.saturating_sub(1)], self.phi[(i + 1).min(n - 1)]);
            dphi[i] = if u > 0.0 { (self.phi[i] - lo) / dx } else { (hi - self.phi[i]) / dx };
            // note: dphi already the spatial derivative; multiply below
            dphi[i] *= u;
        }
        for i in 0..n {
            self.phi[i] -= dt * dphi[i];
            self.mat[i + self.ng] = if self.phi[i] > 0.0 { 1 } else { 0 };
        }
    }

    /// Advance to `t_end` with CFL steps; returns the number of steps.
    pub fn run_to(&mut self, t_end: f64) -> usize {
        let mut steps = 0;
        while self.t < t_end - 1e-14 * t_end.abs() {
            let dt = self.dt_cfl().min(t_end - self.t);
            self.step(dt);
            steps += 1;
        }
        steps
    }

    /// Totals of the conserved variables over the interior.
    pub fn totals(&self) -> Cons {
        let mut s = [0.0; 3];
        for i in 0..self.n {
            for k in 0..3 {
                s[k] += self.u[i + self.ng][k] * self.dx;
            }
        }
        s
    }
}

// =====================================================================
// Lagrangian moving-mesh solver (conservative sharp interface)
// =====================================================================

/// Interface traction and velocity of an HLLC Riemann solver for one
/// material (used at same-material mesh nodes; cheap and nonlinear).
fn hllc_star(m: &Material, wl: &State, wr: &State) -> (f64, f64) {
    let cl = m.sound_speed(wl.rho, wl.p);
    let cr = m.sound_speed(wr.rho, wr.p);
    let sl = (wl.u - cl).min(wr.u - cr);
    let sr = (wl.u + cl).max(wr.u + cr);
    if sl >= 0.0 {
        return (m.tau(wl.rho, wl.p), wl.u);
    }
    if sr <= 0.0 {
        return (m.tau(wr.rho, wr.p), wr.u);
    }
    let (taul, taur) = (m.tau(wl.rho, wl.p), m.tau(wr.rho, wr.p));
    let den = wl.rho * (sl - wl.u) - wr.rho * (sr - wr.u);
    let ss = (taur - taul + wl.rho * wl.u * (sl - wl.u) - wr.rho * wr.u * (sr - wr.u)) / den;
    if ss >= 0.0 {
        (taul + wl.rho * (sl - wl.u) * (ss - wl.u), ss)
    } else {
        (taur + wr.rho * (sr - wr.u) * (ss - wr.u), ss)
    }
}

/// 1D Lagrangian (moving-mesh) finite-volume solver with a sharp material
/// interface at a mesh node.
///
/// Cells keep their mass; nodes move with the interface Riemann velocity.
/// The updates `m du/dt = −(τ*_{j+1} − τ*_j)` and
/// `m de/dt = −(τ*_{j+1}u*_{j+1} − τ*_j u*_j)` make total mass constant and
/// change momentum and energy only through boundary tractions/work: the
/// scheme is **exactly conservative**. The material interface is a mesh node,
/// so it stays sharp and no mass crosses it. Same-material nodes use an HLLC
/// star; a material interface uses the exact two-material Riemann solver.
/// Piecewise-linear (MUSCL) reconstruction makes it second order in smooth
/// regions.
pub struct LagrangianSolver {
    pub n: usize,
    /// cell masses (constant in time)
    pub m: Vec<f64>,
    /// cell velocities
    pub u: Vec<f64>,
    /// cell specific total energy
    pub e: Vec<f64>,
    /// material id per cell
    pub mat: Vec<usize>,
    /// node positions, length `n + 1`
    pub x: Vec<f64>,
    pub mats: [Material; 2],
    pub cfl: f64,
    pub t: f64,
    pub second_order: bool,
    /// Use the linearised (acoustic) two-material Riemann solver at the
    /// material interface instead of the exact nonlinear one. This is exact to
    /// O(amplitude²) and makes the small-amplitude acoustic benchmark cheap;
    /// leave it `false` for finite-amplitude/shock problems.
    pub acoustic_interface: bool,
    /// Accumulated boundary traction impulse `[left, right]`; the exact
    /// momentum ledger for the reflecting walls.
    pub wall_impulse: [f64; 2],
}

impl LagrangianSolver {
    /// Build a uniform mesh from a cell-centre initial state and material.
    pub fn new(n: usize, lo: f64, hi: f64, mats: [Material; 2], init: impl Fn(f64) -> (State, usize)) -> Self {
        let dx = (hi - lo) / n as f64;
        let mut m = vec![0.0; n];
        let mut u = vec![0.0; n];
        let mut e = vec![0.0; n];
        let mut mat = vec![0usize; n];
        let mut x = vec![0.0; n + 1];
        for k in 0..=n {
            x[k] = lo + k as f64 * dx;
        }
        for j in 0..n {
            let (w, mj) = init(0.5 * (x[j] + x[j + 1]));
            let mm = &mats[mj];
            m[j] = w.rho * dx;
            u[j] = w.u;
            e[j] = mm.energy(w.rho, w.p) + 0.5 * w.u * w.u;
            mat[j] = mj;
        }
        Self { n, m, u, e, mat, x, mats, cfl: 0.4, t: 0.0, second_order: true, acoustic_interface: false, wall_impulse: [0.0; 2] }
    }

    #[inline]
    pub fn density(&self, j: usize) -> f64 {
        self.m[j] / (self.x[j + 1] - self.x[j])
    }

    #[inline]
    pub fn prim(&self, j: usize) -> State {
        let rho = self.density(j);
        State {
            rho,
            u: self.u[j],
            p: self.mats[self.mat[j]].p_hydro(rho, self.e[j] - 0.5 * self.u[j] * self.u[j]),
        }
    }

    /// CFL time step from the local longitudinal speed.
    pub fn dt_cfl(&self) -> f64 {
        let mut rate = 0.0f64;
        for j in 0..self.n {
            let w = self.prim(j);
            let c = self.mats[self.mat[j]].sound_speed(w.rho, w.p);
            rate = rate.max((c + w.u.abs()) / (self.x[j + 1] - self.x[j]));
        }
        self.cfl / rate.max(1e-300)
    }

    /// Minmod-limited slope of a cell variable for MUSCL. Reconstruction is
    /// kept strictly inside one material so the interface stays sharp.
    fn slope(&self, j: usize, f: impl Fn(&State) -> f64) -> f64 {
        if j == 0 || j + 1 == self.n {
            return 0.0;
        }
        if self.mat[j - 1] != self.mat[j] || self.mat[j + 1] != self.mat[j] {
            return 0.0;
        }
        let d1 = f(&self.prim(j)) - f(&self.prim(j - 1));
        let d2 = f(&self.prim(j + 1)) - f(&self.prim(j));
        if d1 * d2 <= 0.0 {
            0.0
        } else if d1.abs() < d2.abs() {
            d1
        } else {
            d2
        }
    }

    /// Left/right node states at node `j` (between cells `j-1` and `j`).
    fn node_states(&self, j: usize) -> (State, State) {
        let wl = self.prim(j - 1);
        let wr = self.prim(j);
        if !self.second_order {
            return (wl, wr);
        }
        let sl = [self.slope(j - 1, |w| w.rho), self.slope(j - 1, |w| w.u), self.slope(j - 1, |w| w.p)];
        let sr = [self.slope(j, |w| w.rho), self.slope(j, |w| w.u), self.slope(j, |w| w.p)];
        (
            State { rho: wl.rho + 0.5 * sl[0], u: wl.u + 0.5 * sl[1], p: wl.p + 0.5 * sl[2] },
            State { rho: wr.rho - 0.5 * sr[0], u: wr.u - 0.5 * sr[1], p: wr.p - 0.5 * sr[2] },
        )
    }

    /// RHS of the Lagrangian update `(du, de, dx, boundary_traction)`.
    fn rhs(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>, [f64; 2]) {
        let mut du = vec![0.0; self.n];
        let mut de = vec![0.0; self.n];
        let mut dxn = vec![0.0; self.n + 1];
        let mut pstar = vec![0.0; self.n + 1];
        let mut ustar = vec![0.0; self.n + 1];
        for k in 1..self.n {
            let (wl, wr) = self.node_states(k);
            let (ml, mr) = (self.mat[k - 1], self.mat[k]);
            if ml == mr {
                let (p, uu) = hllc_star(&self.mats[ml], &wl, &wr);
                pstar[k] = p;
                ustar[k] = uu;
            } else if self.acoustic_interface {
                let st = acoustic_star(&self.mats[ml], &wl, &self.mats[mr], &wr);
                pstar[k] = st.tau;
                ustar[k] = st.u;
            } else {
                let st = exact_star(&self.mats[ml], &wl, &self.mats[mr], &wr);
                pstar[k] = st.tau;
                ustar[k] = st.u;
            }
        }
        // fixed reflecting walls
        pstar[0] = self.prim(0).p;
        ustar[0] = 0.0;
        pstar[self.n] = self.prim(self.n - 1).p;
        ustar[self.n] = 0.0;
        for j in 0..self.n {
            du[j] = -(pstar[j + 1] - pstar[j]) / self.m[j];
            de[j] = -(pstar[j + 1] * ustar[j + 1] - pstar[j] * ustar[j]) / self.m[j];
        }
        for k in 1..self.n {
            dxn[k] = ustar[k];
        }
        (du, de, dxn, [pstar[0], pstar[self.n]])
    }

    /// One explicit RK2 (Heun) step.
    pub fn step(&mut self, dt: f64) {
        let (k1u, k1e, k1x, k1bc) = self.rhs();
        let (u0, e0, x0) = (self.u.clone(), self.e.clone(), self.x.clone());
        for j in 0..self.n {
            self.u[j] = u0[j] + dt * k1u[j];
            self.e[j] = e0[j] + dt * k1e[j];
        }
        for k in 1..self.n {
            self.x[k] = x0[k] + dt * k1x[k];
        }
        let (k2u, k2e, k2x, k2bc) = self.rhs();
        for j in 0..self.n {
            self.u[j] = u0[j] + 0.5 * dt * (k1u[j] + k2u[j]);
            self.e[j] = e0[j] + 0.5 * dt * (k1e[j] + k2e[j]);
        }
        for k in 1..self.n {
            self.x[k] = x0[k] + 0.5 * dt * (k1x[k] + k2x[k]);
        }
        // Heun quadrature of the wall tractions: the same average that is
        // applied to the interior momentum.
        for side in 0..2 {
            self.wall_impulse[side] += 0.5 * dt * (k1bc[side] + k2bc[side]);
        }
        self.t += dt;
    }

    pub fn run_to(&mut self, t_end: f64) -> usize {
        let mut steps = 0;
        while self.t < t_end - 1e-14 * t_end.abs() {
            let dt = self.dt_cfl().min(t_end - self.t);
            self.step(dt);
            steps += 1;
        }
        steps
    }

    /// Totals `[mass, momentum, energy]`.
    pub fn totals(&self) -> [f64; 3] {
        let mut s = [0.0; 3];
        for j in 0..self.n {
            s[0] += self.m[j];
            s[1] += self.m[j] * self.u[j];
            s[2] += self.m[j] * self.e[j];
        }
        s
    }

    /// Cell index whose centre is nearest `x` (for probes).
    pub fn cell_at(&self, x: f64) -> usize {
        let mut best = 0;
        let mut d = f64::INFINITY;
        for j in 0..self.n {
            let c = 0.5 * (self.x[j] + self.x[j + 1]);
            if (c - x).abs() < d {
                d = (c - x).abs();
                best = j;
            }
        }
        best
    }
}

// =====================================================================
// Reference materials (same sources as phononic::materials)
// =====================================================================

/// Water at 20 C (CRC Handbook: ρ = 998.2 kg/m³, c = 1482 m/s).
pub fn water() -> Material {
    Material::fluid(998.2, 1482.0)
}

/// Structural steel (Kinsler et al. / NDT handbook: ρ = 7850 kg/m³,
/// c_p = 5960 m/s, c_s = 3230 m/s).
pub fn steel() -> Material {
    Material::solid(7850.0, 5960.0, 3230.0)
}
