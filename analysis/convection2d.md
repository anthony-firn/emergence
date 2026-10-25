# Natural convection in 2D: the unfold from conduction to flow

`compressible::convection2d`
(`fluid/crates/compressible/src/convection2d.rs`,
tests in `fluid/crates/compressible/tests/convection2d.rs`).

This is the next rung of the ladder in `analysis/emergence_ladder.md` after
the isobaric phase-change model (`analysis/isobaric2d.md`). The isobaric rung
solves conduction and the phase-change dilatation exactly while the layer is
stably stratified; the cheap indicator
`Ra = g β ΔT L³/(ν α)` (with `Ra_c = 1707.762` between rigid, isothermal
plates — Drazin & Reid 1981, Chandrasekhar 1961) says when it fails. This
module is the unfold it points to: a two-dimensional incompressible
Boussinesq solver, coupled to temperature advection–diffusion and, through a
Darcy / Carman–Kozeny drag, to the isobaric enthalpy field.

## 1. Method

A staggered (MAC) grid carries the `x`-face velocity `u`, the `y`-face
velocity `v`, and cell-centred pressure `p` and temperature `T`. Second-order
central differences are used for the viscous and (by default) the advective
terms; time stepping is explicit with a CFL / diffusion limit. Each step is a
standard projection,

```text
  u* = u + Δt[ −(u·∇)u + ν∇²u + gβ(T−T_ref) ŷ − A u ],
  ∇²p = ∇·u*/Δt,   u = u* − Δt ∇p,   ∇·u = 0,
```

with `∂p/∂n = 0` on the no-slip walls. The pressure Poisson problem is
solved matrix-free by conjugate gradients with a Jacobi preconditioner,
warm-started from the previous step; the constant null space of the
all-Neumann enclosure is removed each iteration. The domain can be periodic
in `x` (`Boussinesq2D::periodic_x`, the Bénard layer) or closed with no-slip
side walls (the cavity). `Boussinesq2D::cavity(ra, n)` builds the de Vahl
Davis cell; `nusselt_west` / `nusselt_south` give the half-cell Dirichlet wall
flux and `midline_maxima` the benchmark velocities.

**Darcy / Carman–Kozeny coupling.** For phase change the cell drag is the
Voller & Prakash (1987) form

```text
  A = C (1−f_l)² / (f_l³ + q),       q = 1e-3, C = 1e6 1/s,
```

treated semi-implicitly in the momentum predictor,
`u* = (u + Δt F)/(1 + Δt A)` so the solid limit is unconditionally stable.
`f_l` is the liquid mass fraction from the isobaric rung
(`Isobaric2D::liquid_fraction`). The energy is carried by the enthalpy rung,
whose `expansion` is switched off (constant-density Boussinesq), and the
convective velocity advects the volumetric enthalpy conservatively through
`Isobaric2D::advect_enthalpy`. The enthalpy flux is the **mass** flux times
the upwind specific enthalpy, `ρ u h` (not `u h`); the density factor is what
makes the coupling physical — with it the gallium interface leans strongly
(§3), without it the advection is ~6000× too weak and the interface stays
vertical.

**Momentum advection at large cell Reynolds number.** Liquid metals have
`Re_cell = |u|Δx/ν ~ 10²`, where central advection is unstable. Setting
`upwind_momentum = true` switches to a second-order upwind discretisation
(`O(Δx²)` numerical diffusion), which is stable and, unlike first-order
upwind, not too diffusive to permit the convection.

## 2. Benchmark: the de Vahl Davis cavity

The de Vahl Davis (1983) differentially heated square cavity: unit box, west
wall `T = 1`, east wall `T = 0`, adiabatic top/bottom, no-slip, `Pr = 0.71`.
The properties are chosen so that `Ra = gβΔT L³/(νρα)` with `L = ΔT = 1`; the
velocities are in units of `α/L`, exactly the benchmark normalisation.

| `Ra` | `n` | Nu (computed) | Nu (de Vahl Davis) | error | u_max | v_max |
|---|---|---|---|---|---|---|
| `10³` | 24 | 1.1197 | 1.118 | +0.15 % | 3.638 (3.649) | 3.693 (3.697) |
| `10⁴` | 40 | 2.2584 | 2.243 | +0.69 % | 16.105 (16.178) | 19.609 (19.617) |
| `10⁵` | 48 | 4.6042 | 4.519 | +1.89 % | 34.963 (34.73) | 68.031 (68.59) |
| `10⁵` | 64 | 4.5687 | 4.519 | +1.10 % | 34.978 (34.73) | 68.409 (68.59) |

`u_max` and `v_max` are the tabulated mid-plane maxima (`|u|` on `x = L/2`,
`|v|` on `y = L/2`); the global `|u|` maximum is larger, so the mid-plane
definition matters. The grid study is in the `#[ignore]` test
`cavity_ra_1e5_converged_grid_study` (`--ignored --nocapture`). The `Ra = 10⁴`
Nusselt number converges towards the benchmark with resolution
(`Nu(24) = 2.2782`, `Nu(48) = 2.2546`, reference `2.243`), verified by
`cavity_grid_convergence`. All errors are under 2 %, without tuning; the
tests hold 2–5 % tolerances.

## 3. Benchmark: Rayleigh–Bénard onset

A layer of depth 1 in a periodic box of width `2π/k_c` with rigid, isothermal
plates is seeded with `T = 1 − y + ε cos(kx) sin(πy)` (`ε = 1e-4`). The
perturbation energy grows or decays; the logarithmic slope over the second
half of a `t = 6` run is

| `Ra` | growth rate |
|---|---|
| 2000 | **+1.50** |
| 1400 | **−4.12** |

so the rate changes sign between them, bracketing `Ra_c = 1707.762` at the
critical wavenumber `k_c = 3.116`. This is exactly the sign change the
indicator predicts; `rayleigh_benard_growth_changes_sign_at_ra_c` asserts it.

## 4. Coupling to phase change: gallium melting

The Gau & Viskanta (1986) cell — gallium in a rectangular cavity, heated
vertical wall on the left, the right wall at the melting point, adiabatic
top/bottom — is run with the physical properties of Brent, Voller & Reid
(1988): `ρ = 6093 kg/m³`, `ν = 3.0e-7 m²/s`, `c_p = 381.5 J/(kg K)`,
`k = 32 W/(m K)`, `β = 1.2e-4 /K`, `L = 80160 J/kg`. With `ΔT = 8.91 K`
(`Ste = c_p ΔT/L = 0.0424`) and the cavity 6.35 cm wide × 8.94 cm high
(aspect ratio 0.71) the Rayleigh number is `Ra = 1.8 × 10⁶`,
`Pr = ν/α = 0.0216`.

On a 17 × 24 grid at `t = 300 s` (`cargo test`, fast) the coupled run gives

| quantity | value |
|---|---|
| flow speed `max\|u\|` | `1.36 × 10⁻² m/s` |
| melt fraction | 0.587 |
| interface (mid-height) | 0.0483 m |
| interface top − bottom | **+0.0588 m** (top ahead) |
| liquid temperature top − bottom | **+6.85 K** (top hotter) |

The top of the interface has melted through to the far wall while the bottom
is still near the heated wall: the interface **leans strongly with the top
ahead**, the qualitative signature of the Gau–Viskanta experiment. Buoyant
melt rises along the heated wall, flows right along the top, and melts the
upper solid faster.

**Honest comparison.** The published Gau–Viskanta / Brent–Voller–Reid
interface positions could not be digitised here, so a cell-by-cell
quantitative comparison is not claimed. The `Ra` we compute (`1.8 × 10⁶`) is
also higher than some papers quote for the same cell because of the length
and property conventions; the lean is therefore likely exaggerated relative
to the experiment. What is verified is the qualitative behaviour and the
mechanism, not the digitised data.

**Energy conservation.** `Isobaric2D::advect_enthalpy` uses the same flux on
both sides of every face and zero flux at the wall faces, so the total
enthalpy telescopes. With a closed (adiabatic) box, an off-centre hot spot
drives a flow yet the enthalpy is conserved to round-off:

| quantity | drift over the run |
|---|---|
| energy `∫U dV` | **0.0** (below double precision) |

`coupled_energy_is_conserved_to_round_off` asserts `< 1e-11`; the open
Dirichlet-melting case is not energy conserving (heat enters and leaves) and
no assertion is made there.

## 5. The indicator in action

`isobaric2d::{rayleigh_benard, convection_regime, RA_CRIT_RAYLEIGH_BENARD}`
supply the cheap switch. For a 1 cm water layer (`ν = 1e-6`, `α = 1.4e-7`,
`β = 2.07e-4`) `ΔT = 0.01 K` gives `Ra = 145` (conduction) and `ΔT = 1 K`
gives `Ra = 14505` (convection).

On the gallium cell, scaling the wall superheat to move across the critical
value:

| case | `Ra` | `max\|u\|` | melt fraction | top − bottom (lean) |
|---|---|---|---|---|
| super-critical, `ΔT = 8.91 K` | `1.81 × 10⁶` | `1.36 × 10⁻²` | 0.587 | `+5.9 × 10⁻² m` |
| sub-critical, `Ra = 854` | `8.5 × 10²` | `0` | 0.0007 | `0` |

Below `Ra_c` the coupled flow is numerically zero, no convection develops and
the conduction-only isobaric rung is exact; above it the flow is large, the
melt fraction is hundreds of times larger and the melt front leans. This is
the unfold the ladder asks for.

## 6. Limits and open items

* **Boussinesq, constant density.** The energy is carried by the enthalpy
  rung with `expansion = false`, so the ~3 % gallium density change on melting
  is neglected. A variable-density flow coupling needs the momentum equation
  to carry `ρ` and the pressure projection to be variable-coefficient.
* **Darcy drag is a model.** `C = 1e6 1/s` and `q = 1e-3` are Voller–Prakash
  values; the mushy-zone width and the interface velocity depend on them.
* **No interface thermodynamics.** The gallium solid is held at the melting
  point; a subcooled solid, contact resistance or a curved-interface Gibbs–
  Thomson term would refine the comparison.
* **Quantitative gallium data missing.** See §4: the qualitative lean is
  reproduced, the published interface positions are not digitised here.
* **Advection schemes.** The temperature/momentum advection is second order;
  the enthalpy advection is minmod-limited second-order upwind, first order
  at extrema. Shock-free flows only.
* **Bénard from below** is not included: an equilibrium-melt Bénard layer
  seeded with a single roll did not sustain convection in this
  implementation, because the latent heat smears the front over a broad mushy
  band whose Darcy drag pins the thin liquid layer. The side-heated cell
  (§4), where the energy is advected the same way, does convect.

## References

* G. de Vahl Davis, "Natural convection of air in a square cavity: a
  benchmark numerical solution", Int. J. Numer. Methods Fluids **3**, 249
  (1983) (Nusselt numbers `1.118 / 2.243 / 4.519` and mid-plane velocity
  maxima at `Ra = 10³/10⁴/10⁵`, `Pr = 0.71`).
* P. G. Drazin and W. H. Reid, *Hydrodynamic Stability* (1981); S.
  Chandrasekhar, *Hydrodynamic and Hydromagnetic Stability* (1961)
  (`Ra_c = 1707.762`, `k_c = 3.116` for rigid, isothermal plates).
* C. Gau and R. Viskanta, "Melting and solidification of a pure metal on a
  vertical wall", J. Heat Transfer **108**, 174 (1986) (gallium experiment).
* A. D. Brent, V. R. Voller and K. J. Reid, "Enthalpy–porosity technique for
  modeling convection–diffusion phase change: application to the melting of a
  pure metal", Numer. Heat Transfer **13**, 297 (1988) (gallium properties,
  enthalpy–porosity method, melt-front data).
* V. R. Voller and C. Prakash, "A fixed grid numerical modelling methodology
  for convection–diffusion mushy region phase-change problems", Int. J. Heat
  Mass Transfer **30**, 1709 (1987) (Carman–Kozeny Darcy drag).
* Sparrow, Goldstein & Jonsson, J. Fluid Mech. **18**, 513 (1964)
  (`Ra_I,c = 2772`, the volumetric critical value re-exported by the
  isobaric rung).
