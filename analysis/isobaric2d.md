# The isobaric phase-change model in 2D

`compressible::isobaric2d` (`fluid/crates/compressible/src/isobaric2d.rs`,
tests in `fluid/crates/compressible/tests/isobaric2d.rs`).

This is the second rung of the ladder in `analysis/compressible.md` §15. The
1D isobaric model (`isobaric.rs`) is a Lagrangian enthalpy balance in mass
coordinates; the 2D model below is its Eulerian (fixed uniform grid) form, for
melting and freezing where the front is not planar and, eventually, for melt
pools coupled to flow. It now carries the **expansion (contraction) flow**
driven by the density change on melting/freezing, which the earlier fixed-grid
rung neglected.

## 1. Field equations and state recovery

At fixed pressure `p₀`, the model carries the **volumetric enthalpy**

```text
  U = ρ h,        ∂U/∂t + ∇·(U u) = ∇·(k ∇T) + q,
  ∂ρ/∂t + ∇·(ρ u) = 0,
```

with `ρ`, the temperature `T` and the phase mass fractions `y_k` from the
NASG equation of state at `p₀` and the phase equilibrium (`isobaric::PhaseSet`,
`Ladder` on the Gibbs energies). The temperature and the phase state are
recovered from the **specific** enthalpy `h = U/ρ` by an **exact
enthalpy–temperature inversion**:

* Inside one phase `h_k(T)` is monotone in `T`. For a NASG phase without the
  magnetic or tabulated-thermal term, `h = c_p T + c₀` is linear, so
  `T = (h − c₀)/c_p`. The bcc phase of iron carries the CALPHAD magnetic Gibbs
  energy (`isobaric::Magnetic`), so its `h(T)` is nonlinear and the inversion
  uses a safeguarded Newton iteration.
* At a transition `T = T_tr` and the mixture moves along the tie line,
  `h(φ) = (1−φ) h_a + φ h_b`, giving the fraction `φ` in closed form.
* `h` (and hence `U = h/v`) is strictly monotone along the condensed-phase
  equilibrium path (checked in `volumetric_enthalpy_is_monotone`), so the
  state is single valued.

Across a **liquid–vapour** transition the density collapses by ~10³ and
`U = ρh` is no longer single valued; the vapour rung is the mass-coordinate 1D
model, not this one. Tests here keep the ladder in the condensed range.

Conduction is explicit and stable, limited by each cell's *sensible* volumetric
heat capacity `ρ min_k c_{p,k}` against its face conductances (safety factor
0.4), exactly as in `isobaric.rs`. Faces use harmonic conductivities; Dirichlet
sides sit half a cell outside the last centre; periodic pairings wrap the last
row/column. The expansion flow adds a CFL limit `Δt < 0.4 Δx / max|u|`.

## 2. The expansion flow

At constant pressure density is a function of the local equilibrium state, so
mass conservation requires a dilatation `∇·u = −(1/ρ) Dρ/Dt`. In the slow limit
the flow has no vorticity source and is written as `u = ∇ψ`, with the
variable-coefficient Poisson problem

```text
  ∇·(ρ ∇ψ) = −Dρ/Dt = (ρ − ρ^eq)/Δt =: s,     u = ∇ψ,
```

solved each sub-step. The right-hand side is the density change produced by
conduction at a fixed material point: after the conduction update the EOS gives
a provisional `ρ^eq = 1/v(h)`, and the difference from the carried `ρ` is the
volume source. Solving with the *same* face coefficients that update `ρ` makes
the projection exact up to round-off — the cell ends the step at `ρ^eq` — and
the same conservative mass flux advects `U` with upwind specific enthalpy.
Both `∫ρ dV` and `∫U dV` are then conserved by telescoping.

* **Boundary conditions.** A fixed wall is `∂ψ/∂n = 0` (`u = 0`); an open,
  pressure-fixed side is `ψ = 0` (normal velocity free), like the free end of
  the 1D Lagrangian mesh. `ExpansionBc` sets each side to `Wall` or `Open`;
  periodic directions (from `Boundaries2D`) override. With no open side the
  problem is a pure Neumann/periodic enclosure and the source is made
  mean-zero (constant volume).
* **Solver.** Matrix-free conjugate gradients with a Jacobi preconditioner,
  warm-started from the previous sub-step; 1D grids use a direct tridiagonal
  (Thomas) solve. The closed enclosures (all walls) are solved on the
  mean-zero-subtracted source.
* **Negligible flow.** If the relative density change over the step is below
  `1e-6`, the expansion velocity is far below the conduction velocity scale;
  the step then applies the zero-flow projection (`ρ ← ρ^eq`, no advection).
  This is what makes the equal-density verification material — whose only
  residual density change is thermal, `Δρ/ρ ~ 10⁻⁷` per step — run at the old
  speed while the phase-change cases solve the flow.
* `expansion = false` restores the old zero-flow closure, for measuring the
  error it removes.

## 3. Benchmarks and the measured error

### 3.1 1D Stefan/Neumann in a 2D domain

The one-phase plane Stefan problem (solid/liquid at `T_t`, wall at
`T_t + ΔT`, `s = 2λ√(αt)`, `λ e^{λ²} erf λ = St/√π`; Carslaw & Jaeger,
*Conduction of Heat in Solids*, 2nd ed., §11.2) is run on a 2D grid with the
front aligned to the grid (64 cells, `t = 1 s`, closed enclosure):

| case | transition | `Δρ/ρ` on change | zero-flow (before) | expansion (after) |
|---|---|---|---|---|
| model material (equal densities) | 1811 K | 0 | **−0.087 %** | **−0.087 %** |
| iron (`elements::fit`, melting) | 1811 K | +4.9 % | **+10.65 %** | **+1.23 %** |
| water/ice (freezing) | 273.15 K | +9.0 % | **−14.53 %** | **−0.71 %** |

The equal-density case is the method's verification and is unchanged. The
expansion flow closes the iron and water gaps by roughly an order of
magnitude. With a fixed wall at the cold side and open outflow at the far side
(`FlowBc::Open`, the physical free-boundary arrangement, 64 cells) the water
front error is +2.23 %.

### 3.2 Consistency with the 1D Lagrangian rung

The 1D Lagrangian model carries the expansion exactly. Water freezing, 96
cells, `t = 1 s`:

| quantity | front |
|---|---|
| exact Neumann | 0.377112 mm |
| 1D Lagrangian | 0.376939 mm (−0.046 %) |
| 2D Eulerian, expansion | 0.376625 mm (**2D/1D −0.083 %**) |
| 2D Eulerian, zero-flow | 0.325032 mm (2D/1D **−13.77 %**) |

Equal-density material gives 1D vs 2D **+0.029 %**, unchanged.

### 3.3 Radially symmetric freezing with density change

Cylindrical freezing from a cold line at the origin: water at `T_m` freezes
outward, the solid (ice) inside carries the similarity profile
`T = T_m − (Q/4πk_s)[E1(η²) − E1(λ²)]`, the front is `s = 2λ√(α_s t)` with
`λ² e^{λ²} = Q c_s/(4π k_s L)`, and the liquid is isothermal so the expansion
flow is the only density effect. `St_Q = 0.1` (so the near-origin temperature
stays above the ladder floor), `t = 1 s`, closed enclosure, the line sink
regularised over four central cells:

| `nx` | zero-flow error | expansion error |
|---|---|---|
| 24 | +5.7 % | +15.6 % |
| 48 | −8.0 % | **+3.1 %** |
| 64 | −9.7 % | **+1.6 %** |

With the expansion flow the error is first order in the grid (the line sink is
regularised over one cell) and converges; the zero-flow closure is biased and
does not. The coarse `nx = 24` expansion overshoots only because the
source-regularisation error dominates there. The line-source **melting**
(equal-density, gas-free) benchmark is unchanged because its expansion is
skipped as thermal: errors +15.96 % / +3.40 % / +0.82 % at `nx = 24/48/96`,
roughly second order in the front position.

### 3.4 Conservation to round-off

With adiabatic walls and an active expansion flow (iron hot spot, 32 × 32,
density jump 4.9 %):

| quantity | drift over the run |
|---|---|
| energy `∫U dV` | **+1.3 × 10⁻¹⁵ relative** |
| mass `∫ρ dV` | **+1.4 × 10⁻¹⁴ relative** |

The flux updates telescope, so both are exact up to summation round-off.

### 3.5 Grid orientation

The manufactured one-phase solution `T(ξ, t)` with `ξ = (x+y)/√2` is imposed
on all four boundaries (per-cell Dirichlet, `set_dirichlet_row`) and used as the
initial condition; the front is then a diagonal staircase. Equal-density
material, a front resolved over ~10 cells and advanced to `ξ = 0.40` in a unit
box at `nx = 32`:

| orientation | front error |
|---|---|
| aligned (`ξ = x`) | −1.35 % |
| 45° (`ξ = (x+y)/√2`) | +8.58 % |
| **grid-orientation error** | **+9.9 percentage points** |

The diagonal front's staircase remains the dominant error at this resolution.

### 3.6 Conduction vs convection: the unfold indicator

`rayleigh_benard`, `rayleigh_internal`, `convection_regime` and
`Isobaric2D::buoyancy` are the cheap interface to a flow solver. The critical
values are the same as `photothermal::indicators`:
`Ra_c = 1707.762` (Rayleigh–Bénard; Drazin & Reid 1981, Chandrasekhar 1961)
and `Ra_I,c = 2772` (internal heating; Sparrow, Goldstein & Jonsson, J. Fluid
Mech. 18, 513, 1964). For a 1 cm water layer (`ν = 1e-6 m²/s`,
`α = 1.4e-7 m²/s`, `β = 2.07e-4 /K`):

| case | number | regime |
|---|---|---|
| `ΔT = 0.01 K` | `Ra = 145` | conduction (isobaric rung) |
| `ΔT = 1 K` | `Ra = 14 505` | convection → unfold to flow |
| `q = 1e4 W/m³` | `Ra_I = 2.42e4` | convection → unfold to flow |

The flow coupling itself is now implemented, as the next rung:
`compressible::convection2d` (`analysis/convection2d.md`) is a 2D Boussinesq
MAC-projection solver that takes `buoyancy` as its momentum source, advects
the enthalpy field with the resulting divergence-free velocity, and holds the
solid rigid with a Carman–Kozeny Darcy drag built from `liquid_fraction`. The
isobaric rung supplies the source, the switch and the energy update; the
**expansion flow** added here is the *phase-change* dilatation, a separate
irrotational flow, not thermal convection.

Two methods were added to `Isobaric2D` for the coupling.
`liquid_fraction(i)` returns the liquid mass fraction of a cell (0 = rigid
solid) used by the Darcy drag. `advect_enthalpy(ux, uy, dt)` conservatively
advects the volumetric enthalpy `U` with a prescribed staggered velocity,
sub-stepped by CFL 0.4 and using a minmod-limited second-order upwind face
value to keep a moving front monotone. Its enthalpy flux is the **mass** flux
`ρ u` times the upwind specific enthalpy, so for a constant-density Boussinesq
flow with `∇·u = 0` the update reduces to `ρ Dh/Dt = 0`; omitting the density
factor makes the coupling ~`ρ` times too weak. The flux telescopes, so with
adiabatic walls `∫U dV` is conserved to round-off even with an active flow.

## 4. API summary

* `Isobaric2D::new(model, p, k, nx, ny, dx, dy, phase, T0, (T_lo, T_hi))`.
* `advance(dt)` — stable explicit sub-stepping (conduction + expansion flow),
  returns the number of steps.
* `state_h(h) → (T, y)`, `state_at(i) → (T, y, ρ_eq)`, `temperatures()`,
  `states()`, `energy()`, `mass()`, `phase_area(k)`. The old volumetric
  inversion `state(U)` and `u_pure` remain for the zero-flow closure and
  initial conditions.
* `set_cell`, `set_temperatures`, `set_dirichlet_row(side, per-cell T)`.
* `liquid_fraction(i)`, `advect_enthalpy(ux, uy, dt)` — the convection coupling
  (see `analysis/convection2d.md`).
* `source: Vec<f64>` (W/m³), `bc: Boundaries2D`
  (`Dirichlet` / `Adiabatic` / periodic), `flow_bc: ExpansionBc`
  (`Wall` / `Open`), `expansion: bool`, 1D and 2D grids.
* `rayleigh_benard`, `rayleigh_internal`, `convection_regime`,
  `RA_CRIT_RAYLEIGH_BENARD`, `RA_CRIT_INTERNAL`, `buoyancy`.

## 5. Limits and open items

* **Condensed phases only**: the volumetric `U = ρh` inversion is single
  valued only while the density does not collapse. A general boiling or
  condensation problem needs a mass-tracking (specific-enthalpy) formulation
  in 2D, as the 1D rung has.
* **Expansion flow, not convection**: the Marangoni surface stress and the
  buoyancy-driven 2D velocity field are a separate unfold, now implemented in
  `convection2d` on top of this rung. The expansion Poisson here is
  irrotational and constant in time within a step; it does not include a
  momentum equation.
* **Convection coupling is Boussinesq constant density**: `convection2d`
  advects `U` with a divergence-free velocity and neglects the density change
  on melting (small for gallium). A variable-density projection is the next
  step.
* **Interface conductivity**: a two-phase cell uses the volume-fraction
  weighted conductivity. When the phases' conductivities differ strongly
  (water: ice 2.22, liquid 0.561 W/m/K) this under-resolves the flux at the
  front; the 1D rung's layered resistances would be the model to port in each
  grid direction.
* **Radial line-source freezing** has a regularised line sink, so its error is
  first order in the grid; a resolved-source or boundary-driven cylindrical
  set-up would isolate the expansion flow more cleanly.
* Kinetics (nucleation/Avrami) and solid–solid transitions are available in
  the 1D rung but not yet exposed in 2D.

## References

* H. S. Carslaw and J. C. Jaeger, *Conduction of Heat in Solids*, 2nd ed.
  (plane Stefan §11.2; cylindrical/line source §10.4–11.7).
* P. G. Drazin and W. H. Reid, *Hydrodynamic Stability* (1981), and
  S. Chandrasekhar, *Hydrodynamic and Hydromagnetic Stability* (1961):
  `Ra_c = 1707.762`.
* E. M. Sparrow, R. J. Goldstein and V. K. Jonsson, "Thermal instability in a
  horizontal fluid layer…", J. Fluid Mech. **18**, 513 (1964): `Ra_I,c = 2772`.
* Water/ice and iron phase data and conductivities as in
  `analysis/compressible.md` §14–15 and `examples/phase_change.rs`.
