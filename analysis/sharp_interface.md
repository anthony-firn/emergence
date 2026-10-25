# A sharp solid–fluid interface in 1D

Crate: `fluid/crates/compressible`. Module
`compressible::sharp_interface`, validated by `tests/sharp_interface.rs`
(run `cargo test -p compressible --test sharp_interface`, ~16 s total).
The motivation is in `analysis/fsi.md`: the crate's only multi-material model
is the diffuse-interface five-equation model, and at a water/steel interface
it gives ~6% absolute `|T|` error at 1024 cells, `|R| > 1`, and a
shear-bearing GPR mixture that is unstable at that contrast. This note is the
**unfolded** rung: two pure materials, a tracked interface, and the exact
solid–fluid Riemann flux.

## 1. What is modelled

* A **fluid** is a stiffened gas (NASG, `b = q = 0`) with normal stress
  `τ = p`. Water: `ρ = 998.2 kg/m³`, `c = 1482 m/s` (CRC Handbook).
* A **solid** is the 1D reduction of the GPR unified-continuum model
  (`gpr::Gpr`). In 1D the transverse distortion is constant and the
  longitudinal stretch is slaved to the density by `det A = ρ/ρ₀`, so the
  solid is hyperelastic with

  ```
    a      = ρ/ρ₀,
    E₂(ρ)  = (c_s²/6)(a² − 1)²,
    f(ρ)   = (2/3) c_s² ρ a²(a² − 1),     τ = p_h + f(ρ),
    c_L²   = γ(p_h+p∞)/ρ + f'(ρ),
  ```

  which is exactly `τ = ρ² ∂e_tot/∂ρ|_s` and gives the GPR longitudinal speed
  `c_L = sqrt(c₀² + 4c_s²/3)`. Steel: `ρ = 7850 kg/m³`, `c_p = 5960 m/s`,
  `c_s = 3230 m/s` (Kinsler et al. / NDT handbook, the same constants as
  `phononic::materials`), so `c₀ = 4648 m/s` and `c_L = 5960 m/s`. The
  module reproduces `c_L` to 1e-12. The solid is **elastic**; no plastic
  yield surface is in this 1D rung (stated as a limitation in §6).

At normal incidence the shear degrees of freedom decouple, so this
longitudinal stretch is the whole response. The fluid and solid are handled
by the *same* two-material Riemann solver, which is the point of the rung.

## 2. Exact solid–fluid Riemann problem

At a material interface the normal velocity and normal traction are
continuous and the inviscid fluid carries no shear traction. The interface
state `(u*, τ*)` solves

```
u_L + F_L(τ*) = u* = u_R − F_R(τ*),
```

where `F(τ*)` is the velocity jump across the wave joining a side to the
interface: `F = ∫ c dρ/ρ` (evaluated along the isentrope) across a
rarefaction, and `F = (τ*−τ)/J`, `J = sqrt((τ*−τ)/(v−v*))`, across a shock
(Rankine–Hugoniot). The outer root is a bisection on `τ*`; the isentropic
density is a safeguarded Newton solve and the rarefaction integral is
composite Gauss–Legendre with the panel count set by how many decades the fan
spans. `exact_star` returns `(τ*, u*, ρ*_L, ρ*_R)`; `sample_exact` evaluates
the similarity solution at `ξ = x/t` for comparison with a run. An
`acoustic_star` linearisation is also exposed (the exact small-amplitude
limit); the plate benchmark uses it as the interface flux, while the shock
benchmark uses `exact_star`.

**Independent check.** On two ideal gases the module's solver agrees with the
crate's unrelated `exact::star` (Toro's exact Riemann solver) to machine
precision over five standard cases including the strong double-rarefaction
case:

| quantity | worst relative/absolute error |
|---|---|
| interface traction `τ*` | 2.1e-12 |
| interface velocity `u*` | 5.3e-13 m/s |
| star densities `ρ*_L, ρ*_R` | 1.5e-12 |

The solver also reduces to the acoustic (impedance) result as the wave
amplitude goes to zero (`exact_star_reduces_to_acoustic`).

## 3. Two conservative sharp schemes

**Eulerian, level-set tracked** (`SharpSolver`). Each cell holds one pure
material (`φ > 0` solid), WENO5 on same-material stencils, the exact
two-material flux on the interface face, HLLC elsewhere. The interface face
flux is taken from the contact-upwind star state, so it is single-valued and
the update telescopes:

```
dU_i/dt = −(F_{i+½} − F_{i−½})/Δx.
```

Total mass, momentum and total energy over the domain are therefore conserved
to round-off **for any level-set motion**, because the conserved variables
are never relabelled (only the EOS used to interpret `U` changes when a cell
changes tag). The solid elastic energy lives inside the conserved total
energy through the EOS, so there is no separate strain ledger. Verified by
`eulerian_conservation_is_exact` (<1e-8 relative in all three totals) and by
a stationary contact that stays uniform to `|Δp|/p < 1e-6` over 400 steps.

**Lagrangian, moving mesh** (`LagrangianSolver`). The material interface is a
mesh **node**, so it can never diffuse and no mass crosses it. Cell masses
are constant, and the updates

```
m du/dt = −(τ*_{j+1} − τ*_j),   m de/dt = −(τ*_{j+1}u*_{j+1} − τ*_j u*_j)
```

make mass exactly constant, energy exactly constant with rigid walls (the
walls do no work), and momentum change exactly equal to the accumulated wall
traction impulse (`wall_impulse`). Same-material nodes use an HLLC star; a
material interface uses `exact_star` (or `acoustic_star` when
`acoustic_interface` is set). MUSCL/minmod reconstruction is second order
inside one material and disabled next to a material node to keep the
interface sharp. Verified:

| conservation law | residual |
|---|---|
| mass (with walls) | < 1e-13 |
| energy (with walls, `u*_{wall}=0`) | < 1e-11 |
| momentum vs `∫(τ_L − τ_R)dt` | < 1e-10 |

## 4. Benchmarks

### (a) Exact solid–fluid Riemann problem

Covered by §2: machine-precision agreement with the crate's independent exact
Riemann solver, plus the acoustic limit.

### (b) Small-amplitude pulse in water through a steel plate

Setup: 1D water, a 20 mm steel plate, a one-way right-going
Gaussian-modulated 80 kHz pressure pulse (100 Pa, ~1e-7 Mach), probes one cell
either side of the plate; the exact reference is the `phononic::layered::Stack`
transfer matrix of the actual probe–plate–probe geometry (the real probe cell
centres, so there is no fitted phase). `|R|` is `(p_with − p_without)/p_without`
at the front probe, `|T|` is `p_with(back)/p_without(front)`. Errors are
**incident-power-weighted** RMS over 40–120 kHz, where the pulse has its
energy. The solver uses the `acoustic_star` interface flux (exact to
O(amplitude²)); `exact_star` is used independently for the reference.

| N | cells/plate | \|T\| abs err | \|R\| abs err | energy identity \|R\|²+\|T\|²−1 | steps | Δt (s) |
|---|---|---|---|---|---|---|
| 200 | 16 | 5.70e-3 | 1.06e-1 | 1.77e-1 | 1908 | 8.39e-8 |
| 400 | 32 | 2.63e-3 | 4.32e-2 | 7.23e-2 | 3815 | 4.20e-8 |
| 800 | 64 | **1.36e-3** | **2.84e-3** | 4.83e-3 | 7629 | 2.10e-8 |

The absolute `|T|` error falls `5.70e-3 → 2.63e-3 → 1.36e-3`, i.e. order
`log2(5.70/2.63) = 1.12` then `log2(2.63/1.36) = 0.95` — first order in the
plate resolution here (the plate is only 16–64 cells thick and the interface
node is first order). At **64 cells per plate thickness the `|T|` error is
1.4e-3, two orders below the 1% target**, and the water-on-both-sides energy
identity `|R|² + |T|² = 1` holds to 4.8e-3. The at-40-and-80-kHz values at
64 cells are `|T| = 0.0622` vs exact `0.0632`, `|R| = 0.9971` vs exact
`0.9980`.

Run times in the test profile: 0.9 s (N=200), 3.3 s (N=400), 13.4 s (N=800).
The whole benchmark stays under 30 s.

### (c) Strong water shock onto steel

A 1 GPa shock in water hits steel at `p = 1 atm`. The exact two-material
Riemann solution: `τ* = 9.48828e8 Pa`, `u* = 20.1617 m/s`, `ρ*_L = 990.347`
(water, compressed), `ρ*_R = 7876.49` (steel, shocked). The reflected wave in
water is a rarefaction (`τ* < τ_L`), the transmitted wave in steel is a
shock. Run to `t = 5e-5 s` with the Lagrangian solver at N = 1600 (interface
stress uses `exact_star`):

| quantity | error vs exact similarity solution |
|---|---|
| pressure, incident-weighted L1 | 4.65e-3 |
| velocity, mean absolute | 0.193 m/s |
| interface stress `τ*` | < 1e-2 relative |

## 5. Cost and accuracy against the diffuse-interface model

`analysis/fsi.md` measured the five-equation diffuse model on the same
water/steel plate: at N = 1024 (domain 0.6 m, ≈34 cells across the plate) its
absolute `|T|` error was **6.18e-2** and its `|R|` relative error **0.519**;
a single half-space interface over-reflected to peak `|r| = 1.000` against the
Fresnel value 0.939; and the shear-bearing GPR mixture gave peak `|r| = 4.05`
at the steel contrast (unstable).

The sharp interface at a comparable ≈32 cells/plate gives absolute `|T|`
error **2.6e-3** (≈24× smaller) with `|R| = 0.98` (no over-reflection), and
it stays stable at the full contrast. The comparison is not perfectly
apples-to-apples — `fsi.md` bands 20–250 kHz while this note weights 40–120
kHz — but the order-of-magnitude gap and the absence of `|R| > 1` are the
point. Two structural differences explain it:

* the sharp model never mixes the two NASG equations of state, so the
  `p∞` gap (~1e10 Pa) that makes the Kapila mixture speed negative and the
  diffuse interface linearly unstable never enters;
* the interface flux is the exact (or acoustic) two-material Riemann solution,
  so the impedance contrast is imposed rather than approached as the diffuse
  zone is refined.

**When to use which (the fold).** The sharp interface is the resolved rung:
use it when the interface impedance/transmission matters and the layer can be
given a few tens of cells — here 64 cells per plate already gives `|T|` error
1.4e-3. The diffuse five-equation model is the folded/coarse rung: use it when
the layer is **sub-grid**, where many interfaces must be carried without
tracking each one. Its validity indicator is the number of cells across each
material layer: as that goes to zero the physically correct closure is the
effective medium / transfer-matrix homogenisation (`analysis/fsi.md` §2 shows
the diffuse solver following the exact Bragg gap while homogenisation fails as
the indicator becomes O(1)). A sub-grid water/steel layer should therefore be
folded through the exact transfer matrix, not through the diffuse EOS mixing.

## 6. Limits

* **1D normal incidence only.** At oblique incidence the shear modes matter
  and the 1D GPR reduction is not enough; the multi-D extension is the next
  rung (a level-set/GFM in 2D with the full GPR `A` tensor).
* **The solid is elastic.** The 1D reduction has no plastic yield surface;
  "elastic(-plastic)" here means elastic only. Adding a von Mises/lip yield
  (or the GPR relaxation `F_A`) is a separate rung.
* The Eulerian scheme sharpens the interface through retagging; it does not
  carry a separate strain variable, so the elastic energy is only the
  equilibrium part inside `E`. The Lagrangian scheme has no such ambiguity but
  needs a moving mesh.
* The `acoustic_interface` shortcut is only valid at small amplitude; the
  shock benchmark and the conservation tests use the exact Riemann flux.
* The transfer-matrix benchmark is a thin-plate, high-Q measurement, so the
  window must outlast the plate ring-down; the geometry here was chosen so
  the walls do not contaminate the first arrival. Wall parcels of the signal
  are excluded by the one-way source and the short `t_end`.
