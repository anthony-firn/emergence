//! Natural convection: the unfold from conduction to flow, and its coupling to
//! the isobaric phase-change rung.
//!
//! The isobaric rung (`isobaric2d`) solves conduction and the phase-change
//! dilatation. It is exact while the layer is stably stratified; the cheap
//! indicator `Ra = g β ΔT L³/(ν α)` (with `Ra_c = 1707.762` between rigid,
//! isothermal plates, Drazin & Reid 1981 / Chandrasekhar 1961) says when it
//! fails. This module is the unfold it points to: a two-dimensional,
//! incompressible **Boussinesq** solver on the same uniform grid, with the
//! buoyancy source `g β (T − T_ref)`, coupled to temperature
//! advection–diffusion and, through a Darcy / Carman–Kozeny drag, to the
//! isobaric enthalpy field.
//!
//! # Method
//!
//! A staggered (MAC) grid carries the `x`-face velocity `u`, the `y`-face
//! velocity `v`, and the cell-centred pressure `p` and temperature `T`.
//! Advection and diffusion are second order in space (central differences,
//! with a minmod-limited upwind face value in the energy coupling); time
//! stepping is explicit with a CFL/diffusion-limited step. Each step is a
//! standard projection:
//!
//! ```text
//!   u* = u + Δt[ −(u·∇)u + ν∇²u + gβ(T−T_ref) ŷ − A u ],
//!   ∇²p = ∇·u*/Δt,   u = u* − Δt ∇p,   ∇·u = 0,
//! ```
//!
//! with `∂p/∂n = 0` on the no-slip walls. The pressure Poisson problem is
//! matrix-free conjugate gradients with a Jacobi preconditioner, warm-started
//! from the previous step; with no Dirichlet pressure boundary the constant
//! null space is removed (the enclosure has constant volume). The Darcy term
//! `−A u` uses the Voller & Prakash (1987) Carman–Kozeny form
//! `A = C (1−f_l)²/(f_l³ + q)`, treated semi-implicitly so the solid limit is
//! stable; `f_l` is the liquid mass fraction of the isobaric rung, so the solid
//! is held rigid and only the liquid advects the enthalpy.
//!
//! # Benchmarks
//!
//! [`Boussinesq2D::cavity`] runs the de Vahl Davis (1983) differentially
//! heated square cavity (`Pr = 0.71`); the average Nusselt number and the
//! mid-plane velocity maxima reproduce the benchmark table. [`benard_growth`]
//! measures the perturbation growth rate for Rayleigh–Bénard convection
//! between rigid plates and changes sign at the Drazin–Reid value.
//! [`EnthalpyConvection2D`] couples the solver to `Isobaric2D` for melting in
//! the presence of convection (Gau & Viskanta 1986 gallium; Brent, Voller &
//! Reid 1988 enthalpy–porosity).
//!
//! All physical constants in the tests carry their source; the cavity and
//! Bénard cases are dimensionless (only `Ra` and `Pr` enter).

use crate::isobaric2d::{convection_regime, Isobaric2D, ConvectionRegime};

/// Critical Rayleigh number for the onset of steady convection between rigid,
/// isothermal plates (Drazin & Reid 1981, Chandrasekhar 1961).
pub const RA_CRIT_RIGID: f64 = 1707.762;

/// Critical wavenumber (non-dimensionalised by the layer depth) for the same
/// onset. Chandrasekhar (1961), Table II: `k_c = 3.116`.
pub const K_CRIT_RIGID: f64 = 3.116;

/// Diffusivity of the de Vahl Davis benchmark (Prandtl number `0.71`).
pub const PR_DVAHL_DAVIS: f64 = 0.71;

/// Two-dimensional Boussinesq solver on a uniform staggered grid.
pub struct Boussinesq2D {
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
    /// kinematic viscosity, m²/s
    pub nu: f64,
    /// thermal diffusivity, m²/s (only used when `temperature_active`)
    pub alpha: f64,
    /// thermal expansion coefficient, 1/K
    pub beta: f64,
    /// gravitational acceleration, m/s²
    pub g: f64,
    /// reference temperature of the buoyancy, K
    pub t_ref: f64,
    /// cell-centred temperature, K (row-major `j*nx+i`)
    pub t: Vec<f64>,
    /// `x`-face velocity, row-major `j*(nx+1)+i` at `x = i dx`
    pub u: Vec<f64>,
    /// `y`-face velocity, row-major `j*nx+i` at `y = j dy`
    pub v: Vec<f64>,
    /// cell-centred pressure (up to a constant)
    pub p: Vec<f64>,
    /// cell-centred Darcy / Carman–Kozeny drag `A`, 1/s (0 in pure fluid)
    pub drag: Vec<f64>,
    /// periodic in `x` (a Bénard layer); otherwise no-slip side walls
    pub periodic_x: bool,
    /// west/east/south/north wall temperature (`None`: adiabatic)
    pub t_wall: [Option<f64>; 4],
    /// drive the temperature advection–diffusion inside [`Self::advance`]
    /// (`false` when the enthalpy rung owns the energy equation)
    pub temperature_active: bool,
    /// use first-order upwind advection in the momentum equation. Central
    /// differences are second order but unstable once the cell Reynolds
    /// number `|u|Δx/ν` exceeds ~2 (liquid metals do); upwind stabilises it at
    /// the cost of some numerical diffusion.
    pub upwind_momentum: bool,
    /// physical time
    pub time: f64,
    /// last CG iteration count
    pub last_cg_iters: usize,
}

impl Boussinesq2D {
    /// `nx × ny` grid of size `dx × dy` with `nu, alpha, beta, g` and a
    /// buoyancy reference `t_ref`. Wall temperatures start adiabatic.
    #[allow(clippy::too_many_arguments)]
    pub fn new(nx: usize, ny: usize, dx: f64, dy: f64, nu: f64, alpha: f64, beta: f64, g: f64, t_ref: f64) -> Self {
        let n = nx * ny;
        Boussinesq2D {
            nx,
            ny,
            dx,
            dy,
            nu,
            alpha,
            beta,
            g,
            t_ref,
            t: vec![t_ref; n],
            u: vec![0.0; (nx + 1) * ny],
            v: vec![0.0; nx * (ny + 1)],
            p: vec![0.0; n],
            drag: vec![0.0; n],
            periodic_x: false,
            t_wall: [None; 4],
            temperature_active: true,
            upwind_momentum: false,
            time: 0.0,
            last_cg_iters: 0,
        }
    }

    /// `u` on a possibly out-of-range `(i, j)` (`0` on a wall).
    #[inline]
    fn u_cmp(&self, i: isize, j: isize) -> f64 {
        if j < 0 || j >= self.ny as isize {
            return 0.0;
        }
        if self.periodic_x {
            let kk = ((i % self.nx as isize) + self.nx as isize) % self.nx as isize;
            self.u[self.uf(kk as usize, j as usize)]
        } else if i < 0 || i > self.nx as isize {
            0.0
        } else {
            self.u[self.uf(i as usize, j as usize)]
        }
    }

    /// `v` on a possibly out-of-range `(i, j)` (`0` on a wall).
    #[inline]
    fn v_cmp(&self, i: isize, j: isize) -> f64 {
        if j < 0 || j > self.ny as isize {
            return 0.0;
        }
        if self.periodic_x {
            let kk = ((i % self.nx as isize) + self.nx as isize) % self.nx as isize;
            self.v[self.vf(kk as usize, j as usize)]
        } else if i < 0 || i >= self.nx as isize {
            0.0
        } else {
            self.v[self.vf(i as usize, j as usize)]
        }
    }

    #[inline]
    pub fn ci(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }
    #[inline]
    pub fn uf(&self, i: usize, j: usize) -> usize {
        j * (self.nx + 1) + i
    }
    #[inline]
    pub fn vf(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    /// `u` at a possibly wrapped `x`-face index (periodic seam).
    #[inline]
    fn u_at(&self, k: isize, j: usize) -> f64 {
        let n = self.nx as isize;
        let kk = ((k % n) + n) % n;
        self.u[self.uf(kk as usize, j)]
    }

    /// `v` at a possibly wrapped `x`-face index; `0` outside a non-periodic
    /// domain (a side wall).
    #[inline]
    fn v_at(&self, k: isize, j: usize) -> f64 {
        if self.periodic_x {
            let n = self.nx as isize;
            let kk = ((k % n) + n) % n;
            self.v[self.vf(kk as usize, j)]
        } else if k < 0 || k >= self.nx as isize {
            0.0
        } else {
            self.v[self.vf(k as usize, j)]
        }
    }

    /// Cell-centred velocities from the face values (for diagnostics).
    pub fn cell_velocity(&self) -> (Vec<f64>, Vec<f64>) {
        let (nx, ny) = (self.nx, self.ny);
        let mut uc = vec![0.0; nx * ny];
        let mut vc = vec![0.0; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                uc[j * nx + i] = 0.5 * (self.u[self.uf(i, j)] + self.u[self.uf(i + 1, j)]);
                vc[j * nx + i] = 0.5 * (self.v[self.vf(i, j)] + self.v[self.vf(i, j + 1)]);
            }
        }
        (uc, vc)
    }

    /// Largest face velocity.
    pub fn max_speed(&self) -> f64 {
        self.u.iter().chain(self.v.iter()).fold(0.0f64, |m, a| m.max(a.abs()))
    }

    /// Mid-plane velocity maxima `(max|u| on x = Lx/2, max|v| on y = Ly/2)`,
    /// the quantities tabulated by de Vahl Davis (1983).
    pub fn midline_maxima(&self) -> (f64, f64) {
        let (nx, ny) = (self.nx, self.ny);
        let (uc, vc) = self.cell_velocity();
        let mut um = 0.0f64;
        let mut vm = 0.0f64;
        for j in 0..ny {
            for i in 0..nx {
                if i == nx / 2 {
                    um = um.max(uc[j * nx + i].abs());
                }
                if j == ny / 2 {
                    vm = vm.max(vc[j * nx + i].abs());
                }
            }
        }
        (um, vm)
    }

    /// Average Nusselt number over the hot (`west`) wall,
    /// `Nu = ⟨−∂T/∂x⟩ / (ΔT/L)` with the half-cell Dirichlet gradient.
    pub fn nusselt_west(&self) -> f64 {
        let tw = match self.t_wall[0] {
            Some(t) => t,
            None => return f64::NAN,
        };
        let mut s = 0.0;
        for j in 0..self.ny {
            s += 2.0 * (tw - self.t[self.ci(0, j)]) / self.dx;
        }
        s / self.ny as f64
    }

    /// Average Nusselt number over the hot (`south`) wall (Bénard layer).
    pub fn nusselt_south(&self) -> f64 {
        let tw = match self.t_wall[2] {
            Some(t) => t,
            None => return f64::NAN,
        };
        let mut s = 0.0;
        for i in 0..self.nx {
            s += 2.0 * (tw - self.t[self.ci(i, 0)]) / self.dy;
        }
        s / self.nx as f64
    }

    /// Stable explicit step from CFL (`0.4 Δx/|u|`) and diffusion
    /// (`0.4 Δx²/4 max(ν,α)`).
    pub fn stable_dt(&self) -> f64 {
        let dmin = self.dx.min(self.dy);
        let diff = 0.4 * dmin * dmin / (4.0 * self.nu.max(self.alpha));
        let umax = self.max_speed();
        let cfl = if umax > 0.0 { 0.4 * dmin / umax } else { f64::INFINITY };
        diff.min(cfl).max(f64::MIN_POSITIVE)
    }

    /// Advance by `dt`, sub-stepping to the stable limit; returns the number
    /// of sub-steps taken.
    pub fn advance(&mut self, dt: f64) -> usize {
        let mut done = 0.0;
        let mut steps = 0;
        while done < dt * (1.0 - 1e-12) {
            let h = self.stable_dt().min(dt - done).max(f64::MIN_POSITIVE);
            self.step(h);
            done += h;
            steps += 1;
        }
        steps
    }

    /// Momentum predictor, pressure projection and (optionally) the
    /// temperature update, one explicit step of size `h`.
    pub fn step(&mut self, h: f64) {
        self.momentum(h);
        self.project(h);
        if self.temperature_active {
            self.advect_temperature(h);
        }
        self.time += h;
    }

    // ------------------------------------------------------------ momentum

    fn momentum(&mut self, h: f64) {
        let (nx, ny, dx, dy) = (self.nx, self.ny, self.dx, self.dy);
        let nu = self.nu;
        let mut us = self.u.clone();
        let mut vs = self.v.clone();

        // --- u faces (x-momentum)
        let uilo = if self.periodic_x { 0 } else { 1 };
        let uihi = if self.periodic_x { nx } else { nx - 1 };
        for j in 0..ny {
            for i in uilo..uihi {
                let a = self.uf(i, j);
                let uc = self.u[a];
                let vf = 0.25
                    * (self.v_at(i as isize - 1, j)
                        + self.v_at(i as isize - 1, j + 1)
                        + self.v_at(i as isize, j)
                        + self.v_at(i as isize, j + 1));
                let (dudx, dudy) = if self.upwind_momentum {
                    // second-order upwind (unlimited): O(Δx²) numerical
                    // diffusion, stable at large cell Reynolds number
                    let dxu = if uc >= 0.0 {
                        (3.0 * uc - 4.0 * self.u_cmp(i as isize - 1, j as isize) + self.u_cmp(i as isize - 2, j as isize))
                            / (2.0 * dx)
                    } else {
                        (-3.0 * uc + 4.0 * self.u_cmp(i as isize + 1, j as isize) - self.u_cmp(i as isize + 2, j as isize))
                            / (2.0 * dx)
                    };
                    let dyu = if vf >= 0.0 {
                        (3.0 * uc - 4.0 * self.u_cmp(i as isize, j as isize - 1) + self.u_cmp(i as isize, j as isize - 2))
                            / (2.0 * dy)
                    } else {
                        (-3.0 * uc + 4.0 * self.u_cmp(i as isize, j as isize + 1) - self.u_cmp(i as isize, j as isize + 2))
                            / (2.0 * dy)
                    };
                    (dxu, dyu)
                } else {
                    let dudx = if self.periodic_x {
                        (self.u_at(i as isize + 1, j) - self.u_at(i as isize - 1, j)) / (2.0 * dx)
                    } else {
                        (self.u[self.uf(i + 1, j)] - self.u[self.uf(i - 1, j)]) / (2.0 * dx)
                    };
                    let dudy = if ny == 1 {
                        0.0
                    } else if j == 0 {
                        // no-slip ghost u(-dy/2) = -u(0)
                        (self.u[self.uf(i, 1)] + uc) / (2.0 * dy)
                    } else if j == ny - 1 {
                        -uc / dy
                    } else {
                        (self.u[self.uf(i, j + 1)] - self.u[self.uf(i, j - 1)]) / (2.0 * dy)
                    };
                    (dudx, dudy)
                };
                let adv = uc * dudx + vf * dudy;
                let lap_x = if self.periodic_x {
                    (self.u_at(i as isize + 1, j) - 2.0 * uc + self.u_at(i as isize - 1, j)) / (dx * dx)
                } else {
                    (self.u[self.uf(i + 1, j)] - 2.0 * uc + self.u[self.uf(i - 1, j)]) / (dx * dx)
                };
                let lap_y = if ny == 1 {
                    0.0
                } else if j == 0 {
                    (self.u[self.uf(i, 1)] - 3.0 * uc) / (dy * dy)
                } else if j == ny - 1 {
                    (self.u[self.uf(i, j - 1)] - 3.0 * uc) / (dy * dy)
                } else {
                    (self.u[self.uf(i, j + 1)] - 2.0 * uc + self.u[self.uf(i, j - 1)]) / (dy * dy)
                };
                // Darcy drag, averaged to the face, semi-implicit.
                let il = if self.periodic_x {
                    (i + nx - 1) % nx
                } else {
                    i - 1
                };
                let cl = self.drag[j * nx + il];
                let cr = self.drag[j * nx + i];
                let af = 0.5 * (cl + cr);
                us[a] = (uc + h * (-adv + nu * (lap_x + lap_y))) / (1.0 + h * af.max(0.0));
            }
        }

        // --- v faces (y-momentum)
        for j in 1..ny {
            for i in 0..nx {
                let a = self.vf(i, j);
                let vc = self.v[a];
                let uf = 0.25
                    * (self.u[self.uf(i, j - 1)]
                        + self.u[self.uf(i + 1, j - 1)]
                        + self.u[self.uf(i, j)]
                        + self.u[self.uf(i + 1, j)]);
                let (dvdx, dvdy) = if self.upwind_momentum {
                    let dxv = if uf >= 0.0 {
                        (3.0 * vc - 4.0 * self.v_cmp(i as isize - 1, j as isize) + self.v_cmp(i as isize - 2, j as isize))
                            / (2.0 * dx)
                    } else {
                        (-3.0 * vc + 4.0 * self.v_cmp(i as isize + 1, j as isize) - self.v_cmp(i as isize + 2, j as isize))
                            / (2.0 * dx)
                    };
                    let dyv = if vc >= 0.0 {
                        (3.0 * vc - 4.0 * self.v_cmp(i as isize, j as isize - 1) + self.v_cmp(i as isize, j as isize - 2))
                            / (2.0 * dy)
                    } else {
                        (-3.0 * vc + 4.0 * self.v_cmp(i as isize, j as isize + 1) - self.v_cmp(i as isize, j as isize + 2))
                            / (2.0 * dy)
                    };
                    (dxv, dyv)
                } else {
                    let dvdy = (self.v[self.vf(i, j + 1)] - self.v[self.vf(i, j - 1)]) / (2.0 * dy);
                    let dvdx = if self.periodic_x {
                        (self.v_at(i as isize + 1, j) - self.v_at(i as isize - 1, j)) / (2.0 * dx)
                    } else if i == 0 {
                        (self.v[self.vf(1, j)] + vc) / (2.0 * dx)
                    } else if i == nx - 1 {
                        -vc / dx
                    } else {
                        (self.v[self.vf(i + 1, j)] - self.v[self.vf(i - 1, j)]) / (2.0 * dx)
                    };
                    (dvdx, dvdy)
                };
                let adv = uf * dvdx + vc * dvdy;
                let lap_y = (self.v[self.vf(i, j + 1)] - 2.0 * vc + self.v[self.vf(i, j - 1)]) / (dy * dy);
                let lap_x = if self.periodic_x {
                    (self.v_at(i as isize + 1, j) - 2.0 * vc + self.v_at(i as isize - 1, j)) / (dx * dx)
                } else if i == 0 {
                    (self.v[self.vf(1, j)] - 3.0 * vc) / (dx * dx)
                } else if i == nx - 1 {
                    (self.v[self.vf(i - 1, j)] - 3.0 * vc) / (dx * dx)
                } else {
                    (self.v[self.vf(i + 1, j)] - 2.0 * vc + self.v[self.vf(i - 1, j)]) / (dx * dx)
                };
                let tcell = 0.5 * (self.t[self.ci(i, j - 1)] + self.t[self.ci(i, j)]);
                let buoy = self.g * self.beta * (tcell - self.t_ref);
                let ab = 0.5
                    * (self.drag[j * nx + i] + self.drag[(j - 1) * nx + i]);
                vs[a] = (vc + h * (-adv + nu * (lap_x + lap_y) + buoy)) / (1.0 + h * ab.max(0.0));
            }
        }

        // walls fixed at no-slip
        if !self.periodic_x {
            for j in 0..ny {
                us[self.uf(0, j)] = 0.0;
                us[self.uf(nx, j)] = 0.0;
            }
        }
        for i in 0..nx {
            vs[self.vf(i, 0)] = 0.0;
            vs[self.vf(i, ny)] = 0.0;
        }
        self.u = us;
        self.v = vs;
    }

    // ------------------------------------------------------------ projection

    fn project(&mut self, h: f64) {
        let (nx, ny, dx, dy) = (self.nx, self.ny, self.dx, self.dy);
        let n = nx * ny;
        let mut div = vec![0.0; n];
        for j in 0..ny {
            for i in 0..nx {
                let a = self.ci(i, j);
                let ur = if self.periodic_x && i == nx - 1 {
                    self.u[self.uf(0, j)]
                } else {
                    self.u[self.uf(i + 1, j)]
                };
                div[a] = (ur - self.u[self.uf(i, j)]) / dx
                    + (self.v[self.vf(i, j + 1)] - self.v[self.vf(i, j)]) / dy;
            }
        }
        let scale = 1.0 / h;
        let p = self.solve_poisson(&div, scale);
        let mut us = self.u.clone();
        let mut vs = self.v.clone();
        // u faces
        let uilo = if self.periodic_x { 0 } else { 1 };
        let uihi = if self.periodic_x { nx } else { nx - 1 };
        for j in 0..ny {
            for i in uilo..uihi {
                let pm = if self.periodic_x && i == 0 {
                    p[self.ci(nx - 1, j)]
                } else {
                    p[self.ci(i - 1, j)]
                };
                let gx = (p[self.ci(i, j)] - pm) / dx;
                us[self.uf(i, j)] -= h * gx;
            }
        }
        for j in 1..ny {
            for i in 0..nx {
                let gy = (p[self.ci(i, j)] - p[self.ci(i, j - 1)]) / dy;
                vs[self.vf(i, j)] -= h * gy;
            }
        }
        if !self.periodic_x {
            for j in 0..ny {
                us[self.uf(0, j)] = 0.0;
                us[self.uf(nx, j)] = 0.0;
            }
        }
        for i in 0..nx {
            vs[self.vf(i, 0)] = 0.0;
            vs[self.vf(i, ny)] = 0.0;
        }
        self.u = us;
        self.v = vs;
        self.p = p;
    }

    fn apply_lap(&self, p: &[f64], out: &mut [f64]) {
        let (nx, ny, dx, dy) = (self.nx, self.ny, self.dx, self.dy);
        for j in 0..ny {
            for i in 0..nx {
                let a = self.ci(i, j);
                let xm = if self.periodic_x {
                    p[self.ci((i + nx - 1) % nx, j)]
                } else if i == 0 {
                    p[a]
                } else {
                    p[self.ci(i - 1, j)]
                };
                let xp = if self.periodic_x {
                    p[self.ci((i + 1) % nx, j)]
                } else if i == nx - 1 {
                    p[a]
                } else {
                    p[self.ci(i + 1, j)]
                };
                let ym = if j == 0 { p[a] } else { p[self.ci(i, j - 1)] };
                let yp = if j == ny - 1 { p[a] } else { p[self.ci(i, j + 1)] };
                out[a] = (xp - 2.0 * p[a] + xm) / (dx * dx) + (yp - 2.0 * p[a] + ym) / (dy * dy);
            }
        }
    }

    /// Solve `∇²p = div · scale` (Neumann in `y`, periodic or Neumann in `x`).
    fn solve_poisson(&mut self, div: &[f64], scale: f64) -> Vec<f64> {
        let (nx, ny, dx, dy) = (self.nx, self.ny, self.dx, self.dy);
        let n = nx * ny;
        let mut rhs: Vec<f64> = div.iter().map(|d| d * scale).collect();
        let mean = rhs.iter().sum::<f64>() / n as f64;
        for r in rhs.iter_mut() {
            *r -= mean;
        }
        let diag: Vec<f64> = (0..n)
            .map(|a| {
                let i = a % nx;
                let j = a / nx;
                let mut d = 0.0;
                if self.periodic_x {
                    d += 2.0 / (dx * dx);
                } else {
                    if i > 0 {
                        d += 1.0 / (dx * dx);
                    }
                    if i < nx - 1 {
                        d += 1.0 / (dx * dx);
                    }
                }
                if j > 0 {
                    d += 1.0 / (dy * dy);
                }
                if j < ny - 1 {
                    d += 1.0 / (dy * dy);
                }
                d
            })
            .collect();
        let mut x = self.p.clone();
        let mut r = vec![0.0; n];
        self.apply_lap(&x, &mut r);
        for i in 0..n {
            r[i] = rhs[i] - r[i];
        }
        let mr = r.iter().sum::<f64>() / n as f64;
        for v in r.iter_mut() {
            *v -= mr;
        }
        let mut z = vec![0.0; n];
        let mut p = vec![0.0; n];
        let mut ap = vec![0.0; n];
        for i in 0..n {
            z[i] = r[i] / diag[i].max(1e-300);
            p[i] = z[i];
        }
        let mut rz: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
        let rhs_norm = rhs.iter().map(|a| a * a).sum::<f64>().sqrt().max(1e-300);
        let max_iter = (8 * n).min(8000).max(80);
        let mut iters = 0;
        for _ in 0..max_iter {
            iters += 1;
            self.apply_lap(&p, &mut ap);
            let pap: f64 = p.iter().zip(&ap).map(|(a, b)| a * b).sum();
            if pap.abs() < 1e-300 {
                break;
            }
            let alpha = rz / pap;
            for i in 0..n {
                x[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
            }
            let mr = r.iter().sum::<f64>() / n as f64;
            for v in r.iter_mut() {
                *v -= mr;
            }
            let rnorm = r.iter().map(|a| a * a).sum::<f64>().sqrt();
            if rnorm < 1e-8 * rhs_norm {
                break;
            }
            for i in 0..n {
                z[i] = r[i] / diag[i].max(1e-300);
            }
            let rz_new: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
            let beta = rz_new / rz;
            for i in 0..n {
                p[i] = z[i] + beta * p[i];
            }
            rz = rz_new;
        }
        self.last_cg_iters = iters;
        x
    }

    // ------------------------------------------------------------ temperature

    fn advect_temperature(&mut self, h: f64) {
        let (nx, ny, dx, dy) = (self.nx, self.ny, self.dx, self.dy);
        let mut tn = self.t.clone();
        for j in 0..ny {
            for i in 0..nx {
                let a = self.ci(i, j);
                let ue = if self.periodic_x && i == nx - 1 {
                    self.u[self.uf(0, j)]
                } else {
                    self.u[self.uf(i + 1, j)]
                };
                let uw = self.u[self.uf(i, j)];
                let vn = self.v[self.vf(i, j + 1)];
                let vs = self.v[self.vf(i, j)];
                let dtdx = if self.periodic_x {
                    let l = if i == 0 { nx - 1 } else { i - 1 };
                    let r = if i == nx - 1 { 0 } else { i + 1 };
                    (self.t[self.ci(r, j)] - self.t[self.ci(l, j)]) / (2.0 * dx)
                } else if i == 0 {
                    (self.t[self.ci(1, j)] - self.t[a]) / dx
                } else if i == nx - 1 {
                    (self.t[a] - self.t[self.ci(i - 1, j)]) / dx
                } else {
                    (self.t[self.ci(i + 1, j)] - self.t[self.ci(i - 1, j)]) / (2.0 * dx)
                };
                let dtdy = if j == 0 {
                    (self.t[self.ci(i, 1)] - self.t[a]) / dy
                } else if j == ny - 1 {
                    (self.t[a] - self.t[self.ci(i, j - 1)]) / dy
                } else {
                    (self.t[self.ci(i, j + 1)] - self.t[self.ci(i, j - 1)]) / (2.0 * dy)
                };
                let adv = 0.5 * (ue + uw) * dtdx + 0.5 * (vn + vs) * dtdy;
                // diffusion
                let lap_x = if self.periodic_x {
                    let l = if i == 0 { nx - 1 } else { i - 1 };
                    let r = if i == nx - 1 { 0 } else { i + 1 };
                    (self.t[self.ci(r, j)] - 2.0 * self.t[a] + self.t[self.ci(l, j)]) / (dx * dx)
                } else if i == 0 {
                    match self.t_wall[0] {
                        Some(tw) => (self.t[self.ci(1, j)] - 3.0 * self.t[a] + 2.0 * tw) / (dx * dx),
                        None => (self.t[self.ci(1, j)] - self.t[a]) / (dx * dx),
                    }
                } else if i == nx - 1 {
                    match self.t_wall[1] {
                        Some(tw) => (self.t[self.ci(i - 1, j)] - 3.0 * self.t[a] + 2.0 * tw) / (dx * dx),
                        None => (self.t[self.ci(i - 1, j)] - self.t[a]) / (dx * dx),
                    }
                } else {
                    (self.t[self.ci(i + 1, j)] - 2.0 * self.t[a] + self.t[self.ci(i - 1, j)]) / (dx * dx)
                };
                let lap_y = if j == 0 {
                    match self.t_wall[2] {
                        Some(tw) => (self.t[self.ci(i, 1)] - 3.0 * self.t[a] + 2.0 * tw) / (dy * dy),
                        None => (self.t[self.ci(i, 1)] - self.t[a]) / (dy * dy),
                    }
                } else if j == ny - 1 {
                    match self.t_wall[3] {
                        Some(tw) => (self.t[self.ci(i, j - 1)] - 3.0 * self.t[a] + 2.0 * tw) / (dy * dy),
                        None => (self.t[self.ci(i, j - 1)] - self.t[a]) / (dy * dy),
                    }
                } else {
                    (self.t[self.ci(i, j + 1)] - 2.0 * self.t[a] + self.t[self.ci(i, j - 1)]) / (dy * dy)
                };
                tn[a] = self.t[a] + h * (-adv + self.alpha * (lap_x + lap_y));
            }
        }
        self.t = tn;
    }

    // ------------------------------------------------------------ benchmarks

    /// Set up the de Vahl Davis (1983) differentially heated square cavity:
    /// unit box, west wall `T = 1`, east wall `T = 0`, adiabatic top/bottom,
    /// no-slip, `Pr = 0.71`. The properties are chosen so that
    /// `Ra = g β ΔT L³/(ν α)` with `L = ΔT = 1`. Returns the solver with a
    /// linear temperature initial condition.
    pub fn cavity(ra: f64, n: usize) -> Self {
        let alpha = 1.0;
        let nu = PR_DVAHL_DAVIS;
        let g = 1.0;
        let beta = ra * nu * alpha; // ΔT = L = ρ = 1
        let mut s = Boussinesq2D::new(n, n, 1.0 / n as f64, 1.0 / n as f64, nu, alpha, beta, g, 0.5);
        s.t_wall = [Some(1.0), Some(0.0), None, None];
        for j in 0..n {
            for i in 0..n {
                let x = (i as f64 + 0.5) / n as f64;
                s.t[j * n + i] = 1.0 - x;
            }
        }
        s
    }

    /// Integrate the cavity to a steady state (max `ΔT` per outer step below
    /// `tol`) or `max_steps`, whichever comes first. Returns `(steps, Nu)`.
    pub fn run_cavity_steady(&mut self, tol: f64, max_steps: usize) -> (usize, f64) {
        let mut prev = self.t.clone();
        let mut steps = 0;
        while steps < max_steps {
            self.advance(0.02);
            steps += 1;
            let d = self
                .t
                .iter()
                .zip(&prev)
                .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
            prev.clone_from(&self.t);
            if d < tol {
                break;
            }
        }
        (steps, self.nusselt_west())
    }

    /// Number of cells above the critical Rayleigh number, from the cheap
    /// indicator. Convenience wrapper over [`crate::isobaric2d`].
    pub fn regime(ra: f64) -> ConvectionRegime {
        convection_regime(ra)
    }
}

/// Rayleigh–Bénard growth rate between rigid, isothermal plates at wavenumber
/// `k` (periodic box of width `2π/k`, depth `1`, `Pr`). The temperature is
/// seeded with `T = 1 − y + ε cos(kx) sin(πy)`; the perturbation energy
/// `E = Σ(u²+v²)` is tracked and its logarithmic slope over the second half of
/// the run is returned together with `(E_end/E_start)`. A positive rate is the
/// onset of convection. `n` is the number of cells in each direction.
pub fn benard_growth(ra: f64, pr: f64, k: f64, n: usize, t_end: f64) -> (f64, f64) {
    let lx = 2.0 * std::f64::consts::PI / k;
    let dx = lx / n as f64;
    let dy = 1.0 / n as f64;
    let alpha = 1.0;
    let nu = pr;
    let g = 1.0;
    let beta = ra * nu * alpha; // ΔT = L = 1
    let mut s = Boussinesq2D::new(n, n, dx, dy, nu, alpha, beta, g, 0.5);
    s.periodic_x = true;
    s.t_wall = [None, None, Some(1.0), Some(0.0)];
    let eps = 1e-4;
    for j in 0..n {
        for i in 0..n {
            let x = (i as f64 + 0.5) * dx;
            let y = (j as f64 + 0.5) * dy;
            let a = j * n + i;
            s.t[a] = (1.0 - y) + eps * (k * x).cos() * (std::f64::consts::PI * y).sin();
        }
    }
    let energy = |s: &Boussinesq2D| -> f64 {
        s.u.iter().map(|a| a * a).sum::<f64>() + s.v.iter().map(|a| a * a).sum::<f64>()
    };
    let e0 = energy(&s).max(1e-300);
    // two windows of equal length to fit the exponential slope
    let mut e_mid = 0.0;
    let nstep = 400;
    let h = t_end / nstep as f64;
    for step in 0..nstep {
        s.advance(h);
        if step == nstep / 2 - 1 {
            e_mid = energy(&s).max(1e-300);
        }
    }
    let e1 = energy(&s).max(1e-300);
    let rate = (e1 / e_mid.max(1e-300)).ln() / (0.5 * t_end);
    (rate, e1 / e0)
}

/// A phase-change / natural-convection system: the isobaric enthalpy rung
/// carries the energy, the Boussinesq solver carries the momentum, and the
/// solid is held rigid by a Carman–Kozeny Darcy drag in the mushy zone
/// (Voller & Prakash 1987).
pub struct EnthalpyConvection2D {
    pub pc: Isobaric2D,
    pub flow: Boussinesq2D,
    /// Carman–Kozeny constant `C`, 1/s. Larger pins the solid harder.
    pub mushy_c: f64,
    /// cached liquid mass fraction per cell
    pub liquid: Vec<f64>,
}

impl EnthalpyConvection2D {
    /// Couple an `Isobaric2D` (`flow.expansion = false` recommended) to a
    /// `Boussinesq2D` on the same grid. The flow's temperature is driven by
    /// the enthalpy rung, so its own diffusion is switched off here.
    pub fn new(pc: Isobaric2D, mut flow: Boussinesq2D, mushy_c: f64) -> Self {
        assert_eq!(pc.nx, flow.nx);
        assert_eq!(pc.ny, flow.ny);
        assert!((pc.dx - flow.dx).abs() < 1e-15 && (pc.dy - flow.dy).abs() < 1e-15, "grids must match");
        flow.temperature_active = false;
        let n = pc.len();
        EnthalpyConvection2D { pc, flow, mushy_c, liquid: vec![0.0; n] }
    }

    /// Carman–Kozeny drag of a cell with liquid fraction `f`.
    #[inline]
    fn drag_of(&self, f: f64) -> f64 {
        let fl = f.clamp(0.0, 1.0);
        let a = self.mushy_c * (1.0 - fl).powi(2) / (fl.powi(3) + 1e-3);
        a.min(1e12)
    }

    /// One coupled step of size `dt`: momentum + projection with the current
    /// temperature and liquid fraction, conduction/phase change in the
    /// enthalpy rung, then conservative advection of the enthalpy by the
    /// flow.
    pub fn step(&mut self, dt: f64) {
        let n = self.pc.len();
        let temps = self.pc.temperatures();
        for c in 0..n {
            self.flow.t[c] = temps[c];
            let fl = self.pc.liquid_fraction(c);
            self.liquid[c] = fl;
            self.flow.drag[c] = self.drag_of(fl);
        }
        // momentum and incompressibility first, so the energy is advected by
        // the new (divergence-free) velocity
        self.flow.advance(dt);
        self.pc.advance(dt);
        self.pc.advect_enthalpy(&self.flow.u, &self.flow.v, dt);
    }

    /// Integrate for `t_end` in steps of `dt`.
    pub fn run(&mut self, t_end: f64, dt: f64) {
        let mut t = 0.0;
        while t < t_end * (1.0 - 1e-12) {
            let h = dt.min(t_end - t).max(f64::MIN_POSITIVE);
            self.step(h);
            t += h;
        }
    }

    /// The `x` coordinate of the solid–liquid interface at height `j`,
    /// interpolated where the liquid fraction crosses `1/2`. Returns the box
    /// width if the row is fully liquid.
    pub fn front_x(&self, j: usize) -> f64 {
        let nx = self.pc.nx;
        let mut prev = self.liquid[self.pc.idx(0, j)];
        for i in 1..nx {
            let cur = self.liquid[self.pc.idx(i, j)];
            if (prev - 0.5) * (cur - 0.5) <= 0.0 && prev != cur {
                let f = (0.5 - prev) / (cur - prev);
                return (i as f64 - 1.0 + f + 0.5) * self.pc.dx;
            }
            prev = cur;
        }
        nx as f64 * self.pc.dx
    }

    /// The `y` coordinate of the solid–liquid interface at column `i`,
    /// interpolated where the liquid fraction crosses `1/2` (the Bénard
    /// melting geometry, hot below). Returns the box height if fully liquid.
    pub fn front_y(&self, i: usize) -> f64 {
        let ny = self.pc.ny;
        let mut prev = self.liquid[self.pc.idx(i, 0)];
        for j in 1..ny {
            let cur = self.liquid[self.pc.idx(i, j)];
            if (prev - 0.5) * (cur - 0.5) <= 0.0 && prev != cur {
                let f = (0.5 - prev) / (cur - prev);
                return (j as f64 - 1.0 + f + 0.5) * self.pc.dy;
            }
            prev = cur;
        }
        ny as f64 * self.pc.dy
    }

    /// Melted volume fraction (liquid area over the box).
    pub fn melt_fraction(&self) -> f64 {
        let cell = self.pc.dx * self.pc.dy;
        let area: f64 = self.liquid.iter().sum::<f64>() * cell;
        area / (self.pc.nx as f64 * self.pc.dx * self.pc.ny as f64 * self.pc.dy)
    }
}

#[cfg(test)]
mod unit {
    use super::*;

    #[test]
    fn critical_constants() {
        assert!((RA_CRIT_RIGID - 1707.762).abs() < 1e-9);
        assert!((K_CRIT_RIGID - 3.116).abs() < 1e-9);
        assert!((2.0 * std::f64::consts::PI / K_CRIT_RIGID - 2.016).abs() < 0.01);
    }
}
