//! Sharp solid–fluid interface in **2D**, coupled by a real-ghost-fluid /
//! exact 1D solid–fluid Riemann flux in the face-normal direction.
//!
//! This is the 2D unfold of [`crate::sharp_interface`] (see
//! `analysis/sharp_interface.md` and `analysis/fsi.md`). The interface is a
//! tracked level set; every cell holds one pure material, so there is no
//! equation-of-state mixing. The coupled face flux uses the exact one
//! dimensional solid–fluid Riemann solution along the face normal, which is
//! exactly the 1D rung's [`sharp_interface::exact_star`] solver.
//!
//! # Physics
//!
//! * A **fluid** is a stiffened gas (`Nasg`, `b = q = 0`) in 2D Euler form
//!   `[ρ, ρu, ρv, E]`, with normal stress `τ = p`.
//! * A **solid** is linear isotropic elasticity in plane strain, in the
//!   velocity–stress formulation with **compression-positive** stresses
//!   `[ρu, ρv, σxx, σyy, σxy]`:
//!   ```text
//!     ∂u/∂t   = −(∂σxx/∂x + ∂σxy/∂y)/ρ
//!     ∂v/∂t   = −(∂σxy/∂x + ∂σyy/∂y)/ρ
//!     ∂σxx/∂t = −(λ+2μ)∂u/∂x − λ∂v/∂y
//!     ∂σyy/∂t = −λ∂u/∂x − (λ+2μ)∂v/∂y
//!     ∂σxy/∂t = −μ(∂u/∂y + ∂v/∂x)
//!   ```
//!   with `λ = ρ(c_p² − 2c_s²)`, `μ = ρ c_s²` (the same handbook constants as
//!   [`phononic::materials`] and [`sharp_interface::steel`]). The plane-strain
//!   P modulus is `λ+2μ = ρ c_p²`, which is exactly the 1D reduction's
//!   longitudinal modulus, so the interface-normal P wave is the same as in
//!   1D.
//!
//! # Interface coupling
//!
//! At a face whose two cells are different materials the face state is built
//! from **one-sided (upwind) reconstructions** of both sides (a 3-point
//! one-sided interpolation at the face, second order; it falls back to the
//! cell-centre state when a layer is thinner than three cells), and the
//! face-normal state is then solved by the exact 1D solid–fluid Riemann
//! problem:
//!
//! * normal velocity `u_n` and normal traction `τ` are **continuous** and come
//!   from [`sharp_interface::exact_star`];
//! * the inviscid fluid carries no shear traction, so the solid's tangential
//!   traction is set to zero (`σ_nt = 0`) through the outgoing shear
//!   characteristic `u_t^* = u_t ∓ σ_nt/Z_s` (`Z_s = ρ c_s`), i.e. a
//!   free-surface (traction-free) condition for the shear subsystem;
//! * the fluid tangential velocity is free (slip).
//!
//! The 1D normal Riemann solve cannot, by itself, represent the tangential
//! (P↔SV) coupling that makes the solid's normal response depend on the
//! angle. That coupling is restored **inside the flux** by propagating the
//! solid's outgoing P and SV Riemann invariants from the cell centre to the
//! face over half a cell with the transverse source terms of the 2D
//! characteristic equations, `S_P = ∂_∥(λ u_t − c_pσ_nt)` and
//! `S_SV = ∂_∥(μ u_n − c_sσ_tt)` (the tangential derivative `∂_∥` along the
//! face, `σ_tt` the tangential normal stress).
//! This is the standard multidimensional characteristic (NSCBC-type)
//! correction to a normal Riemann solve and is what lets the solid's own P
//! and SV waves respond to the interface traction.
//!
//! Same-material faces use an HLLC flux for the fluid and an exact upwind
//! (P/S) flux for the solid. The update is the unsplit
//! `dU/dt = −(F_{i+½,j}−F_{i−½,j})/Δx − (G_{i,j+½}−G_{i,j−½})/Δy` with one
//! single-valued flux per face; the CFL condition is
//! `Δt = CFL / max_cells[(c+|u|)/Δx + (c+|v|)/Δy]` with `CFL = 0.4` and `c`
//! the material sound speed (`c_p` for the solid). The interface is held
//! fixed: at the small amplitudes here the interface displacement is second
//! order.
//!
//! # What is conserved
//!
//! For a closed domain the fluid mass, momentum and energy and the solid
//! momentum telescope exactly (upwind/HLLC fluxes are conservative). The
//! solid total energy `∫(½ρ|u|² + ½σ:η) dV` is conserved by the compatible
//! velocity–stress flux to round-off except for the O(amplitude²) mismatch
//! in the normal-momentum flux at a *moving* mixed face.
//!
//! # References
//!
//! * Wang, Liu & Khoo, *A real ghost fluid method …*, SIAM J. Sci. Comput. 28
//!   (2006) (real-ghost-fluid interface flux).
//! * Ndanou, Favrie & Gavrilyuk, *Multi-solid and multi-fluid diffuse
//!   interface Godunov method*, J. Comput. Phys. 295 (2015) (exact
//!   solid–fluid Riemann structure used in 1D by `sharp_interface`).
//! * Brekhovskikh & Godin, *Acoustics of Layered Media* (oblique fluid–solid
//!   reflection, mode conversion and the critical angles).
//! * Faran, *Sound scattering by solid cylinders and spheres*, J. Acoust.
//!   Soc. Am. 23, 405 (1951) (curved-interface reference).

use crate::scheme::weno5z;
use crate::sharp_interface::{acoustic_star, exact_star, Material, State};

/// Atmospheric reference pressure [Pa].
pub const P0: f64 = crate::sharp_interface::P0;

/// Default CFL number for the unsplit 2D update.
pub const CFL: f64 = 0.4;

/// A 2D primitive fluid state.
#[derive(Clone, Copy, Debug)]
pub struct FPrim {
    pub rho: f64,
    pub u: f64,
    pub v: f64,
    pub p: f64,
}

/// A 2D primitive solid state (compression-positive stresses).
#[derive(Clone, Copy, Debug)]
pub struct SPrim {
    pub rho: f64,
    pub u: f64,
    pub v: f64,
    pub sxx: f64,
    pub syy: f64,
    pub sxy: f64,
}

/// One face flux: the fluid part is acted on by a fluid cell and the solid
/// part by a solid cell. At a same-material face the unused half is zero.
#[derive(Clone, Copy, Debug, Default)]
struct FF {
    f: [f64; 4],
    s: [f64; 5],
}

/// Rigid-wall mirror in `y` (normal velocity and shear stress flip).
#[inline]
fn reflect_y(mat: u8, q: [f64; 4], s: [f64; 5]) -> ([f64; 4], [f64; 5]) {
    if mat == 0 {
        ([q[0], q[1], -q[2], q[3]], s)
    } else {
        (q, [s[0], -s[1], s[2], s[3], -s[4]])
    }
}

/// Rigid-wall mirror in `x` (normal velocity and shear stress flip).
#[inline]
fn reflect_x(mat: u8, q: [f64; 4], s: [f64; 5]) -> ([f64; 4], [f64; 5]) {
    if mat == 0 {
        ([q[0], -q[1], q[2], q[3]], s)
    } else {
        (q, [-s[0], s[1], s[2], s[3], -s[4]])
    }
}

/// The 2D sharp solid–fluid solver.
///
/// Cell `(i,j)` (interior `0..nx`, `0..ny`) is solid when the level set
/// `phi(i,j) > 0` and fluid otherwise. The level set is not advected; the
/// interface is static at the linear amplitudes this module targets.
pub struct Sharp2D {
    pub nx: usize,
    pub ny: usize,
    pub ng: usize,
    pub dx: f64,
    pub dy: f64,
    /// signed level set at interior cell centres
    pub phi: Vec<f64>,
    /// padded material id, `0` fluid `1` solid
    pub mat: Vec<u8>,
    /// padded fluid conserved state `[ρ, ρu, ρv, E]`
    pub q: Vec<[f64; 4]>,
    /// padded solid state `[ρu, ρv, σxx, σyy, σxy]`
    pub s: Vec<[f64; 5]>,
    /// `[fluid, solid]`
    pub mats: [Material; 2],
    /// Lamé parameters of `mats[1]`
    pub lam: f64,
    pub mu: f64,
    pub periodic_x: bool,
    /// Reflecting (rigid) `y` walls. When false the `y` ghosts copy the
    /// interior (radiating/zero-gradient outflow), which is used for the
    /// scattering benchmarks.
    pub wall_y: bool,
    /// Reflecting `x` walls when not periodic.
    pub wall_x: bool,
    pub acoustic_interface: bool,
    pub cfl: f64,
    pub t: f64,
    w: usize,
}

impl Sharp2D {
    /// Build a grid on `[0,lx] × [0,ly]` from a level set, a fluid primitive
    /// closure and a solid primitive closure.
    pub fn new(
        nx: usize,
        ny: usize,
        ng: usize,
        lx: f64,
        ly: f64,
        mats: [Material; 2],
        phi: impl Fn(f64, f64) -> f64,
        fluid: impl Fn(f64, f64) -> FPrim,
        solid: impl Fn(f64, f64) -> SPrim,
    ) -> Self {
        let dx = lx / nx as f64;
        let dy = ly / ny as f64;
        let w = nx + 2 * ng;
        let h = ny + 2 * ng;
        let mut phi_v = vec![0.0; nx * ny];
        let mut mat = vec![0u8; w * h];
        let mut q = vec![[0.0; 4]; w * h];
        let mut s = vec![[0.0; 5]; w * h];
        for j in 0..ny {
            for i in 0..nx {
                let x = (i as f64 + 0.5) * dx;
                let y = (j as f64 + 0.5) * dy;
                let p = phi(x, y);
                phi_v[j * nx + i] = p;
                let m = if p > 0.0 { 1 } else { 0 };
                mat[(j + ng) * w + (i + ng)] = m;
                if m == 0 {
                    let f = fluid(x, y);
                    let e = mats[0].eos.energy(f.rho, f.p) + 0.5 * (f.u * f.u + f.v * f.v);
                    q[(j + ng) * w + (i + ng)] = [f.rho, f.rho * f.u, f.rho * f.v, f.rho * e];
                } else {
                    let d = solid(x, y);
                    s[(j + ng) * w + (i + ng)] =
                        [d.rho * d.u, d.rho * d.v, d.sxx, d.syy, d.sxy];
                }
            }
        }
        let c_p = mats[1].c_longitudinal();
        let lam = mats[1].rho0 * (c_p * c_p - 2.0 * mats[1].cs * mats[1].cs);
        let mu = mats[1].rho0 * mats[1].cs * mats[1].cs;
        let mut sol = Self {
            nx,
            ny,
            ng,
            dx,
            dy,
            phi: phi_v,
            mat,
            q,
            s,
            mats,
            lam,
            mu,
            periodic_x: true,
            wall_y: false,
            wall_x: false,
            acoustic_interface: false,
            cfl: CFL,
            t: 0.0,
            w,
        };
        sol.fill_ghosts_inplace();
        sol
    }

    #[inline]
    fn at(&self, i: isize, j: isize) -> usize {
        ((j + self.ng as isize) as usize) * self.w + (i + self.ng as isize) as usize
    }

    /// Primitive fluid state of a padded cell.
    #[inline]
    fn fprim_at(&self, q: &[[f64; 4]], i: isize, j: isize) -> FPrim {
        let c = &q[self.at(i, j)];
        let rho = c[0].max(1e-300);
        let u = c[1] / rho;
        let v = c[2] / rho;
        let e_tot = c[3] / rho - 0.5 * (u * u + v * v);
        FPrim { rho, u, v, p: self.mats[0].eos.pressure(rho, e_tot) }
    }

    /// Primitive solid state of a padded cell.
    #[inline]
    fn sprim_at(&self, s: &[[f64; 5]], i: isize, j: isize) -> SPrim {
        let c = &s[self.at(i, j)];
        let rho = self.mats[1].rho0;
        SPrim {
            rho,
            u: c[0] / rho,
            v: c[1] / rho,
            sxx: c[2],
            syy: c[3],
            sxy: c[4],
        }
    }

    /// Admissibility of a reconstructed fluid state.
    #[inline]
    fn fadm(&self, w: &FPrim) -> bool {
        w.rho > 0.0 && w.p + self.mats[0].eos.p_inf > 0.0 && w.rho.is_finite() && w.p.is_finite()
    }

    /// Public fluid primitive accessor (interior indices).
    pub fn fluid_prim(&self, i: usize, j: usize) -> FPrim {
        self.fprim_at(&self.q, i as isize, j as isize)
    }

    /// Public solid primitive accessor (interior indices).
    pub fn solid_prim(&self, i: usize, j: usize) -> SPrim {
        self.sprim_at(&self.s, i as isize, j as isize)
    }

    pub fn is_solid(&self, i: usize, j: usize) -> bool {
        self.mat[self.at(i as isize, j as isize)] == 1
    }

    /// CFL time step.
    pub fn dt_cfl(&self) -> f64 {
        let mut rate = 0.0f64;
        for j in 0..self.ny {
            for i in 0..self.nx {
                let m = &self.mats[self.mat[self.at(i as isize, j as isize)] as usize];
                let (c, ua, va) = if self.mat[self.at(i as isize, j as isize)] == 0 {
                    let w = self.fluid_prim(i, j);
                    (m.sound_speed(w.rho, w.p), w.u.abs(), w.v.abs())
                } else {
                    let w = self.solid_prim(i, j);
                    (m.c_longitudinal(), w.u.abs(), w.v.abs())
                };
                rate = rate.max((c + ua) / self.dx + (c + va) / self.dy);
            }
        }
        self.cfl / rate.max(1e-300)
    }

    /// Fill the padded ghosts: periodic in `x`, zero-gradient in `y`.
    pub fn fill_ghosts_inplace(&mut self) {
        let (ng, nx, ny, w) = (self.ng, self.nx, self.ny, self.w);
        // y: reflecting wall or radiating outflow
        for i in 0..w {
            let lo_m = self.mat[ng * w + i];
            let hi_m = self.mat[(ng + ny - 1) * w + i];
            let qlo = self.q[ng * w + i];
            let qhi = self.q[(ng + ny - 1) * w + i];
            let slo = self.s[ng * w + i];
            let shi = self.s[(ng + ny - 1) * w + i];
            let (qlo, slo) = if self.wall_y { reflect_y(lo_m, qlo, slo) } else { (qlo, slo) };
            let (qhi, shi) = if self.wall_y { reflect_y(hi_m, qhi, shi) } else { (qhi, shi) };
            for g in 0..ng {
                self.mat[g * w + i] = lo_m;
                self.q[g * w + i] = qlo;
                self.s[g * w + i] = slo;
                self.mat[(ng + ny + g) * w + i] = hi_m;
                self.q[(ng + ny + g) * w + i] = qhi;
                self.s[(ng + ny + g) * w + i] = shi;
            }
        }
        // x: periodic wrap, or rigid wall when not periodic
        if self.periodic_x {
            for jr in 0..(ny + 2 * ng) {
                for g in 0..ng {
                    let src_lo = jr * w + (ng + nx - 1 - g);
                    let src_hi = jr * w + (ng + g);
                    let dst_lo = jr * w + g;
                    let dst_hi = jr * w + (ng + nx + g);
                    self.mat[dst_lo] = self.mat[src_lo];
                    self.q[dst_lo] = self.q[src_lo];
                    self.s[dst_lo] = self.s[src_lo];
                    self.mat[dst_hi] = self.mat[src_hi];
                    self.q[dst_hi] = self.q[src_hi];
                    self.s[dst_hi] = self.s[src_hi];
                }
            }
        } else {
            for jr in 0..(ny + 2 * ng) {
                let lo = jr * w + ng;
                let hi = jr * w + (ng + nx - 1);
                let (ql, sl) = if self.wall_x {
                    reflect_x(self.mat[lo], self.q[lo], self.s[lo])
                } else {
                    (self.q[lo], self.s[lo])
                };
                let (qh, sh) = if self.wall_x {
                    reflect_x(self.mat[hi], self.q[hi], self.s[hi])
                } else {
                    (self.q[hi], self.s[hi])
                };
                for g in 0..ng {
                    self.mat[jr * w + g] = self.mat[lo];
                    self.q[jr * w + g] = ql;
                    self.s[jr * w + g] = sl;
                    self.mat[jr * w + (ng + nx + g)] = self.mat[hi];
                    self.q[jr * w + (ng + nx + g)] = qh;
                    self.s[jr * w + (ng + nx + g)] = sh;
                }
            }
        }
    }

    fn filled(&self, q: &[[f64; 4]], s: &[[f64; 5]]) -> (Vec<[f64; 4]>, Vec<[f64; 5]>) {
        let (ng, nx, ny, w) = (self.ng, self.nx, self.ny, self.w);
        let mut q = q.to_vec();
        let mut s = s.to_vec();
        for i in 0..w {
            let qlo = q[ng * w + i];
            let qhi = q[(ng + ny - 1) * w + i];
            let slo = s[ng * w + i];
            let shi = s[(ng + ny - 1) * w + i];
            let (qlo, slo) = if self.wall_y {
                reflect_y(self.mat[ng * w + i], qlo, slo)
            } else {
                (qlo, slo)
            };
            let (qhi, shi) = if self.wall_y {
                reflect_y(self.mat[(ng + ny - 1) * w + i], qhi, shi)
            } else {
                (qhi, shi)
            };
            for g in 0..ng {
                q[g * w + i] = qlo;
                s[g * w + i] = slo;
                q[(ng + ny + g) * w + i] = qhi;
                s[(ng + ny + g) * w + i] = shi;
            }
        }
        if self.periodic_x {
            for jr in 0..(ny + 2 * ng) {
                for g in 0..ng {
                    let src_lo = jr * w + (ng + nx - 1 - g);
                    let src_hi = jr * w + (ng + g);
                    let dst_lo = jr * w + g;
                    let dst_hi = jr * w + (ng + nx + g);
                    q[dst_lo] = q[src_lo];
                    s[dst_lo] = s[src_lo];
                    q[dst_hi] = q[src_hi];
                    s[dst_hi] = s[src_hi];
                }
            }
        } else {
            for jr in 0..(ny + 2 * ng) {
                let lo = jr * w + ng;
                let hi = jr * w + (ng + nx - 1);
                let (ql, sl) = if self.wall_x {
                    reflect_x(self.mat[lo], q[lo], s[lo])
                } else {
                    (q[lo], s[lo])
                };
                let (qh, sh) = if self.wall_x {
                    reflect_x(self.mat[hi], q[hi], s[hi])
                } else {
                    (q[hi], s[hi])
                };
                for g in 0..ng {
                    q[jr * w + g] = ql;
                    s[jr * w + g] = sl;
                    q[jr * w + (ng + nx + g)] = qh;
                    s[jr * w + (ng + nx + g)] = sh;
                }
            }
        }
        (q, s)
    }

    // -----------------------------------------------------------------
    // Second-order (MUSCL/minmod) reconstruction
    // -----------------------------------------------------------------

    /// Fluid primitive states at the low and high face of cell `(i,j)`,
    /// WENO5. Reconstruction is disabled on a stencil that touches another
    /// material, so the interface stays sharp.
    fn recon_fluid(&self, q: &[[f64; 4]], i: isize, j: isize, dir: u8) -> (FPrim, FPrim) {
        let w1 = self.fprim_at(q, i, j);
        let off = |k: isize| if dir == 0 { (i + k, j) } else { (i, j + k) };
        let m = self.mat[self.at(i, j)];
        let ok = (-2..=2).all(|k| {
            let (a, b) = off(k);
            self.mat[self.at(a, b)] == m
        });
        if !ok {
            return (w1, w1);
        }
        let p = |k: isize| self.fprim_at(q, off(k).0, off(k).1);
        let c = |f: fn(&FPrim) -> f64, k: isize| f(&p(k));
        let low = |f: fn(&FPrim) -> f64| weno5z(c(f, 2), c(f, 1), c(f, 0), c(f, -1), c(f, -2));
        let high = |f: fn(&FPrim) -> f64| weno5z(c(f, -2), c(f, -1), c(f, 0), c(f, 1), c(f, 2));
        (
            FPrim {
                rho: low(|w| w.rho),
                u: low(|w| w.u),
                v: low(|w| w.v),
                p: low(|w| w.p),
            },
            FPrim {
                rho: high(|w| w.rho),
                u: high(|w| w.u),
                v: high(|w| w.v),
                p: high(|w| w.p),
            },
        )
    }

    /// Solid primitive states at the low and high face of cell `(i,j)`, WENO5.
    fn recon_solid(&self, s: &[[f64; 5]], i: isize, j: isize, dir: u8) -> (SPrim, SPrim) {
        let w1 = self.sprim_at(s, i, j);
        let off = |k: isize| if dir == 0 { (i + k, j) } else { (i, j + k) };
        let m = self.mat[self.at(i, j)];
        let ok = (-2..=2).all(|k| {
            let (a, b) = off(k);
            self.mat[self.at(a, b)] == m
        });
        if !ok {
            return (w1, w1);
        }
        let p = |k: isize| self.sprim_at(s, off(k).0, off(k).1);
        let c = |f: fn(&SPrim) -> f64, k: isize| f(&p(k));
        let low = |f: fn(&SPrim) -> f64| weno5z(c(f, 2), c(f, 1), c(f, 0), c(f, -1), c(f, -2));
        let high = |f: fn(&SPrim) -> f64| weno5z(c(f, -2), c(f, -1), c(f, 0), c(f, 1), c(f, 2));
        let build = |ux: f64, vx: f64, sxx: f64, syy: f64, sxy: f64| SPrim {
            rho: w1.rho,
            u: ux,
            v: vx,
            sxx,
            syy,
            sxy,
        };
        (
            build(low(|w| w.u), low(|w| w.v), low(|w| w.sxx), low(|w| w.syy), low(|w| w.sxy)),
            build(high(|w| w.u), high(|w| w.v), high(|w| w.sxx), high(|w| w.syy), high(|w| w.sxy)),
        )
    }

    /// One-sided (upwind) fluid primitive at a face of cell `(i,j)`: `upper`
    /// selects the high face, and the stencil is taken from inside the same
    /// material only. Falls back to the cell-centre state when the material
    /// layer is thinner than three cells.
    fn fprim_face(&self, q: &[[f64; 4]], i: isize, j: isize, dir: u8, upper: bool) -> FPrim {
        let m = self.mat[self.at(i, j)];
        let sgn: isize = if upper { -1 } else { 1 };
        let off = |k: isize| if dir == 0 { (i + sgn * k, j) } else { (i, j + sgn * k) };
        let ok = |k: isize| {
            let (a, b) = off(k);
            self.mat[self.at(a, b)] == m
        };
        if !ok(1) || !ok(2) {
            return self.fprim_at(q, i, j);
        }
        let p = |k: isize| {
            let (a, b) = off(k);
            self.fprim_at(q, a, b)
        };
        let (w0, w1, w2) = (1.875, -1.25, 0.375);
        let (a, b, c) = (p(0), p(1), p(2));
        FPrim {
            rho: w0 * a.rho + w1 * b.rho + w2 * c.rho,
            u: w0 * a.u + w1 * b.u + w2 * c.u,
            v: w0 * a.v + w1 * b.v + w2 * c.v,
            p: w0 * a.p + w1 * b.p + w2 * c.p,
        }
    }

    /// One-sided (upwind) solid primitive at a face of cell `(i,j)`; see
    /// [`Self::fprim_face`].
    fn sprim_face(&self, s: &[[f64; 5]], i: isize, j: isize, dir: u8, upper: bool) -> SPrim {
        let m = self.mat[self.at(i, j)];
        let sgn: isize = if upper { -1 } else { 1 };
        let off = |k: isize| if dir == 0 { (i + sgn * k, j) } else { (i, j + sgn * k) };
        let ok = |k: isize| {
            let (a, b) = off(k);
            self.mat[self.at(a, b)] == m
        };
        if !ok(1) || !ok(2) {
            return self.sprim_at(s, i, j);
        }
        let p = |k: isize| {
            let (a, b) = off(k);
            self.sprim_at(s, a, b)
        };
        let (w0, w1, w2) = (1.875, -1.25, 0.375);
        let (a, b, c) = (p(0), p(1), p(2));
        SPrim {
            rho: self.mats[1].rho0,
            u: w0 * a.u + w1 * b.u + w2 * c.u,
            v: w0 * a.v + w1 * b.v + w2 * c.v,
            sxx: w0 * a.sxx + w1 * b.sxx + w2 * c.sxx,
            syy: w0 * a.syy + w1 * b.syy + w2 * c.syy,
            sxy: w0 * a.sxy + w1 * b.sxy + w2 * c.sxy,
        }
    }

    // -----------------------------------------------------------------
    // Face fluxes
    // -----------------------------------------------------------------

    /// HLLC flux for a fluid–fluid face. `dir = 0` is x, `dir = 1` is y.
    fn hllc_fluid(&self, wl: &FPrim, wr: &FPrim, dir: u8) -> [f64; 4] {
        let m = &self.mats[0];
        let (unl, utl) = if dir == 0 { (wl.u, wl.v) } else { (wl.v, wl.u) };
        let (unr, utr) = if dir == 0 { (wr.u, wr.v) } else { (wr.v, wr.u) };
        let el = m.eos.energy(wl.rho, wl.p) + 0.5 * (wl.u * wl.u + wl.v * wl.v);
        let er = m.eos.energy(wr.rho, wr.p) + 0.5 * (wr.u * wr.u + wr.v * wr.v);
        let flux = |rho: f64, un: f64, ut: f64, p: f64, e: f64| -> [f64; 4] {
            [rho * un, rho * un * un + p, rho * ut * un, (rho * e + p) * un]
        };
        let cl = m.sound_speed(wl.rho, wl.p);
        let cr = m.sound_speed(wr.rho, wr.p);
        let sl = (unl - cl).min(unr - cr);
        let sr = (unl + cl).max(unr + cr);
        let f = if sl >= 0.0 {
            flux(wl.rho, unl, utl, wl.p, el)
        } else if sr <= 0.0 {
            flux(wr.rho, unr, utr, wr.p, er)
        } else {
            let den = wl.rho * (sl - unl) - wr.rho * (sr - unr);
            let ss = (wr.p - wl.p + wl.rho * unl * (sl - unl) - wr.rho * unr * (sr - unr)) / den;
            let star = |rho: f64, un: f64, ut: f64, p: f64, e: f64, s: f64| -> [f64; 4] {
                let fac = rho * (s - un) / (s - ss);
                let e_star = fac * (e + (ss - un) * (ss + p / (rho * (s - un))));
                [fac, fac * ss, fac * ut, e_star]
            };
            if ss >= 0.0 {
                let fl = flux(wl.rho, unl, utl, wl.p, el);
                let w = [wl.rho, wl.rho * unl, wl.rho * utl, wl.rho * el];
                let st = star(wl.rho, unl, utl, wl.p, el, sl);
                std::array::from_fn(|k| fl[k] + sl * (st[k] - w[k]))
            } else {
                let fr = flux(wr.rho, unr, utr, wr.p, er);
                let w = [wr.rho, wr.rho * unr, wr.rho * utr, wr.rho * er];
                let st = star(wr.rho, unr, utr, wr.p, er, sr);
                std::array::from_fn(|k| fr[k] + sr * (st[k] - w[k]))
            }
        };
        // map [ρ, n-mom, t-mom, E] back to [ρ, ρu, ρv, E]
        if dir == 0 {
            f
        } else {
            [f[0], f[2], f[1], f[3]]
        }
    }

    /// Exact upwind (P and S) flux for a solid–solid face.
    fn upwind_solid(&self, wl: &SPrim, wr: &SPrim, dir: u8) -> [f64; 5] {
        let rho = self.mats[1].rho0;
        let c_p = self.mats[1].c_longitudinal();
        let c_s = self.mats[1].cs;
        let zp = rho * c_p;
        let zs = rho * c_s;
        if dir == 0 {
            // P: (u, σxx); S: (v, σxy)
            let un = 0.5 * (wl.u + wr.u) + 0.5 * (wl.sxx - wr.sxx) / zp;
            let snn = 0.5 * (wl.sxx + wr.sxx) + 0.5 * zp * (wl.u - wr.u);
            let ut = 0.5 * (wl.v + wr.v) + 0.5 * (wl.sxy - wr.sxy) / zs;
            let snt = 0.5 * (wl.sxy + wr.sxy) + 0.5 * zs * (wl.v - wr.v);
            [snn, snt, (self.lam + 2.0 * self.mu) * un, self.lam * un, self.mu * ut]
        } else {
            // P: (v, σyy); S: (u, σxy)
            let un = 0.5 * (wl.v + wr.v) + 0.5 * (wl.syy - wr.syy) / zp;
            let snn = 0.5 * (wl.syy + wr.syy) + 0.5 * zp * (wl.v - wr.v);
            let ut = 0.5 * (wl.u + wr.u) + 0.5 * (wl.sxy - wr.sxy) / zs;
            let snt = 0.5 * (wl.sxy + wr.sxy) + 0.5 * zs * (wl.u - wr.u);
            [snt, snn, self.lam * un, (self.lam + 2.0 * self.mu) * un, self.mu * ut]
        }
    }

    /// Tangential source terms `(S_P, S_SV)` for the solid's outgoing `P`
    /// and `SV` Riemann invariants at cell `(i,j)`, for an interface whose
    /// normal is `dir` (`0` = x, `1` = y). They are the components of the
    /// 2D characteristic equations that a purely normal Riemann solve drops:
    /// with normal `y`, `S_P = ∂ₓ(λu − c_pσ_xy)`,
    /// `S_SV = ∂ₓ(μv − c_sσ_xx)`; with normal `x` the roles of `(u,v)` and
    /// `(σ_xx,σ_yy)` swap. These terms are exactly what converts a normal
    /// traction into a mix of P and SV motion.
    fn char_sources(&self, s: &[[f64; 5]], i: isize, j: isize, dir: u8) -> (f64, f64) {
        let c_p = self.mats[1].c_longitudinal();
        let c_s = self.mats[1].cs;
        let (a, b, inv) = if dir == 0 {
            ((i, j - 1), (i, j + 1), 1.0 / (2.0 * self.dy))
        } else {
            ((i - 1, j), (i + 1, j), 1.0 / (2.0 * self.dx))
        };
        let wl = self.sprim_at(s, a.0, a.1);
        let wr = self.sprim_at(s, b.0, b.1);
        if dir == 0 {
            // normal x, tangential y
            let fp = self.lam * wl.v - c_p * wl.sxy;
            let fr = self.lam * wr.v - c_p * wr.sxy;
            let gp = self.mu * wl.u - c_s * wl.syy;
            let gr = self.mu * wr.u - c_s * wr.syy;
            ((fr - fp) * inv, (gr - gp) * inv)
        } else {
            // normal y, tangential x
            let fp = self.lam * wl.u - c_p * wl.sxy;
            let fr = self.lam * wr.u - c_p * wr.sxy;
            let gp = self.mu * wl.v - c_s * wl.sxx;
            let gr = self.mu * wr.v - c_s * wr.sxx;
            ((fr - fp) * inv, (gr - gp) * inv)
        }
    }

    /// Mixed fluid–solid face flux. `fl` is the fluid cell, `sl` the solid
    /// cell; `fluid_lower` says which side of the face the fluid is on.
    /// `sp_src`/`ss_src` are [`Self::char_sources`] at the solid cell.
    fn interface_flux(
        &self,
        fl: &FPrim,
        sl: &SPrim,
        dir: u8,
        fluid_lower: bool,
        sp_src: f64,
        ss_src: f64,
    ) -> FF {
        let fmat = &self.mats[0];
        let smat = &self.mats[1];
        let c_p = smat.c_longitudinal();
        let c_s = smat.cs;
        let (un_f, ut_f, sn_f) = if dir == 0 { (fl.u, fl.v, fl.p) } else { (fl.v, fl.u, fl.p) };
        let (un_s, ut_s, sn_s) = if dir == 0 {
            (sl.u, sl.v, sl.sxx)
        } else {
            (sl.v, sl.u, sl.syy)
        };
        // Propagate the solid's outgoing characteristics from its cell centre
        // to the face over half a cell, including the tangential (P<->SV)
        // source terms that the 1D normal solve omits. For a solid on the
        // lower/left side the outgoing invariant is ρc v + σ and the traction
        // correction has the opposite sign to the upper/right case.
        let dt_p = if dir == 0 { self.dx / (2.0 * c_p) } else { self.dy / (2.0 * c_p) };
        let dt_s = if dir == 0 { self.dx / (2.0 * c_s) } else { self.dy / (2.0 * c_s) };
        let sign = if fluid_lower { -1.0 } else { 1.0 };
        let sn_s = sn_s + sign * sp_src * dt_p;
        // 1D Riemann in the normal direction (lower side first).
        let (wl, ml, wr, mr) = if fluid_lower {
            (
                State { rho: fl.rho, u: un_f, p: sn_f },
                fmat,
                State { rho: sl.rho, u: un_s, p: sn_s },
                smat,
            )
        } else {
            (
                State { rho: sl.rho, u: un_s, p: sn_s },
                smat,
                State { rho: fl.rho, u: un_f, p: sn_f },
                fmat,
            )
        };
        let st = if self.acoustic_interface {
            acoustic_star(ml, &wl, mr, &wr)
        } else {
            exact_star(ml, &wl, mr, &wr)
        };
        let tau = st.tau;
        let ustar = st.u;
        let rho_star = if fluid_lower { st.rho_l } else { st.rho_r };
        // fluid flux
        let e_star = fmat.eos.energy(rho_star, tau) + 0.5 * (ustar * ustar + ut_f * ut_f);
        let fflux = if dir == 0 {
            [
                rho_star * ustar,
                rho_star * ustar * ustar + tau,
                rho_star * ut_f * ustar,
                (rho_star * e_star + tau) * ustar,
            ]
        } else {
            [
                rho_star * ustar,
                rho_star * ut_f * ustar,
                rho_star * ustar * ustar + tau,
                (rho_star * e_star + tau) * ustar,
            ]
        };
        // solid flux: free-surface shear characteristic, corrected by the
        // tangential P<->SV source so the solid's own shear motion responds to
        // the interface normal traction.
        let zs = sl.rho * c_s;
        let ut_star = if fluid_lower {
            ut_s - sl.sxy / zs + ss_src * dt_s / zs
        } else {
            ut_s + sl.sxy / zs + ss_src * dt_s / zs
        };
        let sflux = if dir == 0 {
            [
                tau,                      // normal momentum σxx
                0.0,                      // tangential momentum, σxy* = 0
                (self.lam + 2.0 * self.mu) * ustar,
                self.lam * ustar,
                self.mu * ut_star,
            ]
        } else {
            [
                0.0,                      // tangential momentum, σxy* = 0
                tau,                      // normal momentum σyy
                self.lam * ustar,
                (self.lam + 2.0 * self.mu) * ustar,
                self.mu * ut_star,
            ]
        };
        FF { f: fflux, s: sflux }
    }

    /// Accumulate the conservative RHS over all faces.
    fn rhs(&self, q: &[[f64; 4]], s: &[[f64; 5]]) -> (Vec<[f64; 4]>, Vec<[f64; 5]>) {
        let (q, s) = self.filled(q, s);
        let (nx, ny) = (self.nx, self.ny);
        // x faces
        let mut xf = vec![FF::default(); (nx + 1) * ny];
        for j in 0..ny {
            for i in 0..=nx {
                let il = i as isize - 1;
                let ir = i as isize;
                let ml = self.mat[self.at(il, j as isize)];
                let mr = self.mat[self.at(ir, j as isize)];
                let idx = j * (nx + 1) + i;
                xf[idx] = if ml == mr {
                    if ml == 0 {
                        let wl = self.recon_fluid(&q, il, j as isize, 0).1;
                        let wr = self.recon_fluid(&q, ir, j as isize, 0).0;
                        let wl = if self.fadm(&wl) { wl } else { self.fprim_at(&q, il, j as isize) };
                        let wr = if self.fadm(&wr) { wr } else { self.fprim_at(&q, ir, j as isize) };
                        FF { f: self.hllc_fluid(&wl, &wr, 0), s: [0.0; 5] }
                    } else {
                        let wl = self.recon_solid(&s, il, j as isize, 0).1;
                        let wr = self.recon_solid(&s, ir, j as isize, 0).0;
                        FF { f: [0.0; 4], s: self.upwind_solid(&wl, &wr, 0) }
                    }
                } else {
                    let (wl_f, wr_s, wr_f, wl_s) = if ml == 0 {
                        (
                            self.fprim_face(&q, il, j as isize, 0, true),
                            self.sprim_face(&s, ir, j as isize, 0, false),
                            self.fprim_at(&q, ir, j as isize),
                            self.sprim_at(&s, il, j as isize),
                        )
                    } else {
                        (
                            self.fprim_at(&q, il, j as isize),
                            self.sprim_at(&s, ir, j as isize),
                            self.fprim_face(&q, ir, j as isize, 0, false),
                            self.sprim_face(&s, il, j as isize, 0, true),
                        )
                    };
                    if ml == 0 {
                        let (sp, ss) = self.char_sources(&s, ir, j as isize, 0);
                        self.interface_flux(&wl_f, &wr_s, 0, true, sp, ss)
                    } else {
                        let (sp, ss) = self.char_sources(&s, il, j as isize, 0);
                        self.interface_flux(&wr_f, &wl_s, 0, false, sp, ss)
                    }
                };
            }
        }
        // y faces
        let mut yf = vec![FF::default(); nx * (ny + 1)];
        for j in 0..=ny {
            for i in 0..nx {
                let jl = j as isize - 1;
                let jr = j as isize;
                let ml = self.mat[self.at(i as isize, jl)];
                let mr = self.mat[self.at(i as isize, jr)];
                let idx = j * nx + i;
                yf[idx] = if ml == mr {
                    if ml == 0 {
                        let wl = self.recon_fluid(&q, i as isize, jl, 1).1;
                        let wr = self.recon_fluid(&q, i as isize, jr, 1).0;
                        let wl = if self.fadm(&wl) { wl } else { self.fprim_at(&q, i as isize, jl) };
                        let wr = if self.fadm(&wr) { wr } else { self.fprim_at(&q, i as isize, jr) };
                        FF { f: self.hllc_fluid(&wl, &wr, 1), s: [0.0; 5] }
                    } else {
                        let wl = self.recon_solid(&s, i as isize, jl, 1).1;
                        let wr = self.recon_solid(&s, i as isize, jr, 1).0;
                        FF { f: [0.0; 4], s: self.upwind_solid(&wl, &wr, 1) }
                    }
                } else {
                    let (wl_f, wr_s, wr_f, wl_s) = if ml == 0 {
                        (
                            self.fprim_face(&q, i as isize, jl, 1, true),
                            self.sprim_face(&s, i as isize, jr, 1, false),
                            self.fprim_at(&q, i as isize, jr),
                            self.sprim_at(&s, i as isize, jl),
                        )
                    } else {
                        (
                            self.fprim_at(&q, i as isize, jl),
                            self.sprim_at(&s, i as isize, jr),
                            self.fprim_face(&q, i as isize, jr, 1, false),
                            self.sprim_face(&s, i as isize, jl, 1, true),
                        )
                    };
                    if ml == 0 {
                        let (sp, ss) = self.char_sources(&s, i as isize, jr, 1);
                        self.interface_flux(&wl_f, &wr_s, 1, true, sp, ss)
                    } else {
                        let (sp, ss) = self.char_sources(&s, i as isize, jl, 1);
                        self.interface_flux(&wr_f, &wl_s, 1, false, sp, ss)
                    }
                };
            }
        }
        let mut dq = vec![[0.0; 4]; nx * ny];
        let mut ds = vec![[0.0; 5]; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;
                let xl = &xf[j * (nx + 1) + i];
                let xr = &xf[j * (nx + 1) + i + 1];
                let yl = &yf[j * nx + i];
                let yr = &yf[(j + 1) * nx + i];
                if self.mat[self.at(i as isize, j as isize)] == 0 {
                    for k in 0..4 {
                        dq[idx][k] = -(xr.f[k] - xl.f[k]) / self.dx - (yr.f[k] - yl.f[k]) / self.dy;
                    }
                } else {
                    for k in 0..5 {
                        ds[idx][k] = -(xr.s[k] - xl.s[k]) / self.dx - (yr.s[k] - yl.s[k]) / self.dy;
                    }
                }
            }
        }
        (dq, ds)
    }

    /// One SSP-RK3 step.
    pub fn step(&mut self, dt: f64) {
        let (k1q, k1s) = self.rhs(&self.q, &self.s);
        let mut q1 = self.q.clone();
        let mut s1 = self.s.clone();
        self.combine(&mut q1, &mut s1, &k1q, &k1s, dt, 1.0, 0.0);
        let (k2q, k2s) = self.rhs(&q1, &s1);
        let mut q2 = q1.clone();
        let mut s2 = s1.clone();
        self.combine(&mut q2, &mut s2, &k2q, &k2s, dt, 0.25, 0.75);
        let (k3q, k3s) = self.rhs(&q2, &s2);
        for j in 0..self.ny {
            for i in 0..self.nx {
                let p = self.at(i as isize, j as isize);
                let r = j * self.nx + i;
                for k in 0..4 {
                    self.q[p][k] = self.q[p][k] / 3.0
                        + (2.0 / 3.0) * (q2[p][k] + dt * k3q[r][k]);
                }
                for k in 0..5 {
                    self.s[p][k] = self.s[p][k] / 3.0
                        + (2.0 / 3.0) * (s2[p][k] + dt * k3s[r][k]);
                }
            }
        }
        self.t += dt;
    }

    /// `out = a*base + b*(stage + dt*k)` over interior cells.
    fn combine(
        &self,
        q_out: &mut [[f64; 4]],
        s_out: &mut [[f64; 5]],
        kq: &[[f64; 4]],
        ks: &[[f64; 5]],
        dt: f64,
        b: f64,
        a: f64,
    ) {
        for j in 0..self.ny {
            for i in 0..self.nx {
                let p = self.at(i as isize, j as isize);
                let r = j * self.nx + i;
                for k in 0..4 {
                    q_out[p][k] = a * self.q[p][k] + b * (q_out[p][k] + dt * kq[r][k]);
                }
                for k in 0..5 {
                    s_out[p][k] = a * self.s[p][k] + b * (s_out[p][k] + dt * ks[r][k]);
                }
            }
        }
    }

    /// Advance to `t_end`; returns the number of steps.
    pub fn run_to(&mut self, t_end: f64) -> usize {
        let mut steps = 0;
        while self.t < t_end - 1e-14 * t_end.abs() {
            let dt = self.dt_cfl().min(t_end - self.t);
            self.step(dt);
            steps += 1;
        }
        steps
    }

    // -----------------------------------------------------------------
    // Diagnostics
    // -----------------------------------------------------------------

    /// Fluid totals `[mass, x-momentum, y-momentum, energy]` over the domain.
    pub fn fluid_totals(&self) -> [f64; 4] {
        let mut acc = [0.0; 4];
        for j in 0..self.ny {
            for i in 0..self.nx {
                if self.mat[self.at(i as isize, j as isize)] != 0 {
                    continue;
                }
                let c = self.q[self.at(i as isize, j as isize)];
                for k in 0..4 {
                    acc[k] += c[k] * self.dx * self.dy;
                }
            }
        }
        acc
    }

    /// Solid totals `[x-momentum, y-momentum]` over the domain.
    pub fn solid_momentum(&self) -> [f64; 2] {
        let mut acc = [0.0; 2];
        for j in 0..self.ny {
            for i in 0..self.nx {
                if self.mat[self.at(i as isize, j as isize)] != 1 {
                    continue;
                }
                let c = self.s[self.at(i as isize, j as isize)];
                acc[0] += c[0] * self.dx * self.dy;
                acc[1] += c[1] * self.dx * self.dy;
            }
        }
        acc
    }

    /// Elastic strain energy density of a compression-positive stress state.
    fn elastic_energy_density(&self, w: &SPrim) -> f64 {
        let (l, m) = (self.lam, self.mu);
        let inv = 1.0 / (4.0 * m * (l + m));
        let exx = inv * ((l + 2.0 * m) * w.sxx - l * w.syy);
        let eyy = inv * ((l + 2.0 * m) * w.syy - l * w.sxx);
        let exy = w.sxy / (2.0 * m);
        0.5 * (w.sxx * exx + w.syy * eyy + 2.0 * w.sxy * exy)
    }

    /// Solid kinetic + elastic energy over the domain (excluding the constant
    /// rest-state elastic energy).
    pub fn solid_energy(&self) -> f64 {
        let mut e = 0.0;
        let rho = self.mats[1].rho0;
        for j in 0..self.ny {
            for i in 0..self.nx {
                if self.mat[self.at(i as isize, j as isize)] != 1 {
                    continue;
                }
                let w = self.solid_prim(i, j);
                e += (0.5 * rho * (w.u * w.u + w.v * w.v) + self.elastic_energy_density(&w))
                    * self.dx
                    * self.dy;
            }
        }
        e
    }

    /// Total fluid energy `∫E dV` over the domain.
    pub fn fluid_energy(&self) -> f64 {
        let mut e = 0.0;
        for j in 0..self.ny {
            for i in 0..self.nx {
                if self.mat[self.at(i as isize, j as isize)] != 0 {
                    continue;
                }
                let c = self.q[self.at(i as isize, j as isize)];
                e += c[3] * self.dx * self.dy;
            }
        }
        e
    }

    /// Total (fluid + solid) perturbation energy.
    pub fn total_energy(&self) -> f64 {
        self.fluid_energy() + self.solid_energy()
    }
}

/// The 2D plane-strain Lamé parameters of `sharp_interface::steel`-style
/// handbook constants, exposed for tests.
pub fn lame(rho: f64, c_p: f64, c_s: f64) -> (f64, f64) {
    (rho * (c_p * c_p - 2.0 * c_s * c_s), rho * c_s * c_s)
}

/// Oblique reflection coefficient of a plane acoustic wave at a fluid–solid
/// **half-space** (the exact P/SV solution; Brekhovskikh & Godin). `theta` is
/// the incidence angle from the normal. The solid's input normal impedance is
/// `Z_in = Z_L cos²(2θ_T) + Z_T sin²(2θ_T)` with the oblique normal
/// impedances `Z_L = ρ c_p/cosθ_L`, `Z_T = ρ c_s/cosθ_T` and the fluid
/// `Z_f = ρ_f c_f/cosθ`; `R = (Z_in − Z_f)/(Z_in + Z_f)`. Past a critical
/// angle the corresponding cosine continues as `i|cos|`, so `|R| = 1`.
/// This is the same P/SV physics as `phononic::fluid_coupling`; it agrees
/// with that crate's layered Fresnel result at normal incidence.
pub fn halfspace_reflection(rho_f: f64, c_f: f64, rho_s: f64, c_p: f64, c_s: f64, theta: f64) -> (f64, f64) {
    // complex square root sqrt(x), branch with non-negative imaginary part
    let csqrt = |x: f64| -> (f64, f64) {
        if x >= 0.0 {
            (x.sqrt(), 0.0)
        } else {
            (0.0, (-x).sqrt())
        }
    };
    // multiply and divide (a+ib)(c+id)
    let cmul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    let cdiv = |a: (f64, f64), b: (f64, f64)| {
        let d = b.0 * b.0 + b.1 * b.1;
        ((a.0 * b.0 + a.1 * b.1) / d, (a.1 * b.0 - a.0 * b.1) / d)
    };
    // oblique normal impedance Z = rho c / cos(theta), complex past critical
    let zimp = |rho: f64, cc: f64, s: f64| -> (f64, f64) {
        let (cr, ci) = csqrt(1.0 - s * s);
        cmul((rho * cc, 0.0), cdiv((1.0, 0.0), (cr, ci)))
    };
    let sin_t = theta.sin();
    let cos_f = theta.cos();
    let s_l = c_p / c_f * sin_t;
    let s_t = c_s / c_f * sin_t;
    let zl = zimp(rho_s, c_p, s_l);
    let zt = zimp(rho_s, c_s, s_t);
    // cos(2θ_T) real, sin(2θ_T) = 2 sinθ_T cosθ_T complex
    let c2 = 1.0 - 2.0 * s_t * s_t;
    let (ctr, cti) = csqrt(1.0 - s_t * s_t);
    let s2 = (2.0 * s_t * ctr, 2.0 * s_t * cti);
    let s2sq = cmul(s2, s2);
    let zl_term = (zl.0 * c2 * c2, zl.1 * c2 * c2);
    let zt_term = cmul(zt, s2sq);
    let zin = (zl_term.0 + zt_term.0, zl_term.1 + zt_term.1);
    let zf = (rho_f * c_f / cos_f, 0.0);
    let num = (zin.0 - zf.0, zin.1 - zf.1);
    let den = (zin.0 + zf.0, zin.1 + zf.1);
    cdiv(num, den)
}

/// Convenience: a stiffened-gas fluid `Material` at `(rho, c)`.
pub fn fluid_material(rho: f64, c: f64) -> Material {
    Material::fluid(rho, c)
}

/// Convenience: a GPR 1D-reduction solid `Material` at `(rho, c_p, c_s)`.
pub fn solid_material(rho: f64, c_p: f64, c_s: f64) -> Material {
    Material::solid(rho, c_p, c_s)
}

/// Re-export the 1D material helpers so callers need only this module.
pub use crate::sharp_interface::{steel, water};
