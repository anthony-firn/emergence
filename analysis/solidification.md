# Crystal growth from the melt: phase-field solidification and its folds

This note fills the item left open in [compressible.md](compressible.md) §15,
"dendritic ice and solidification microstructure", and is an instance of the
ladder in [emergence_ladder.md](emergence_ladder.md):

```
  isobaric enthalpy model        macroscopic cell, no interface       (coarse)
       ^  indicator: tip diameter / dendrite spacing vs dx_macro, and
       |             W0/d0 (interface resolution)
  Karma-Rappel phase field       anisotropic, W0/d0 -> 0               (fine)
       ^  thin-interface limit
  sharp-interface Stefan + Gibbs-Thomson                              (exact)
       ^
  1D Neumann similarity solution                                      (exact)
```

The crate is `fluid/crates/solidification/`. Everything below is reproducible
with `cargo test -p solidification` (unit + benchmark tests) and the examples
`dendrite`, `gibbs_thomson`, `isobaric_closure`.

## 1. The model

`phase_field.rs` integrates the pure-substance Karma–Rappel model
(Karma & Rappel, *Phys. Rev. E* **57**, 4323 (1998); the 1996 paper is
*Phys. Rev. E* **54**, 3065). Phase field `p = +1` solid, `-1` liquid;
dimensionless temperature `u = (T - T_m)/(L/c_p)`:

```
  tau(theta) dp/dt = div( W(theta)^2 grad p ) + (p - lambda u (1-p^2))(1-p^2)
  du/dt            = D lap u + (1/2) dp/dt
```

with `W(theta) = W0 (1 + eps4 cos 4 theta)`, `theta` the interface-normal
angle. The Karma–Rappel kinetic-free choice is used:

```
  lambda = a1 W0 / d0,   tau0 = a2 lambda W0^2 / D,
  a1 = 5 sqrt(2)/8 = 0.883883..., a2 = 0.6267.
```

`tau(theta) = tau0 (W(theta)/W0)^2` is the anisotropic extension that keeps the
*angular* kinetic coefficient zero as well: with `W(theta) = W0 a(theta)`,
`tau = a2 lambda W0^2 a^2 / D`, so the interface sits on the local
Gibbs–Thomson temperature `u_i = -d0(theta) kappa` up to `O(W0)` products of
curvature and the diffusion field.

Implementation: uniform grid, cell-centred, explicit forward Euler,
second-order central differences, ghost-cell boundary conditions
(symmetry / far field), rows advanced in parallel with rayon. The anisotropic
Laplacian is expanded as `W^2 lap p + 2 W W' (grad theta . grad p)`, with
`theta_x = (p_x p_xy - p_y p_xx)/|grad p|^2` and the symmetric form for
`theta_y`. Two operators are exposed: `PhaseField::step` is this plain
operator, and `PhaseField::step_kr` is the corrected Karma–Rappel operator
(KR 1998 Eq. 88): KR's 9-point isotropic Laplacian plus the variational
anisotropy term `p_x d_y Q - p_y d_x Q`, `Q = W dW/dtheta`. At `eps4 = 0` the
variational term vanishes and `step_kr` is bit-for-bit `step`, so the
isotropic benchmarks below are unchanged.

`stefan.rs` holds the exact Neumann similarity solution and a 1D enthalpy
solver. `closure.rs` holds the isobaric fold. `materials.rs` holds sourced
properties.

## 2. Benchmarks

### (a) 1D planar Stefan problem vs the exact Neumann solution

A pure melt at `T_m` against a wall at `T_w < T_m`. The exact front is
`s = 2 lambda sqrt(alpha t)` with `lambda e^{lambda^2} erf(lambda) = St/sqrt(pi)`
(Carslaw & Jaeger, *Conduction of Heat in Solids*, 2nd ed., §11.2). The 1D
enthalpy solver (`Stefan1d`) marches the scaled enthalpy
`y = (T - T_m) + (L/c) f_liquid` and inverts it exactly in each cell.

| quantity | model | exact | error |
|---|---|---|---|
| front at `t = 0.5`, `St = 0.063`, `n = 400` | 0.249683 | 0.248424 | **0.51 %** |
| `T` a quarter into the new phase | 269.3688 K | 269.3947 K | 0.026 K |

At `n = 800` and `t = 0.1, 0.2, 0.5, 1.0` the front error falls
`+0.56, +0.40, +0.25, +0.18 %`: the transient offset is the first cell being
fully solid while the film is thinner than a cell, and vanishes with `dx`.
This reproduces the same law quoted for the isobaric model in
`analysis/compressible.md` §15 (which used 400 cells and got 0.0–1.9 % for
materials with density change; this constant-property reference is 0.5 %).

### (b) Gibbs–Thomson: critical radius of a circular seed

A seed of radius `R` in a melt at `u_far = -Delta` grows above
`R_c = d0/Delta` and shrinks below it. Bisection on the sign of the front
displacement (`gibbs_thomson` example, `d0 = 1`, `D = 2`, `Delta = 0.3`,
`d0/Delta = 3.3333`):

| `W0/d0` | `dx` | `W0/dx` | measured `R_c` | `R_c/(d0/Delta)` |
|---|---|---|---|---|
| 0.8 | 0.25 | 3.2 | 3.3354 | 1.0006 |
| 0.4 | 0.25 | 1.6 | 3.3094 | 0.9928 |
| 0.2 | 0.25 | 0.8 | 3.1620 | 0.9486 |
| 0.1 | 0.25 | 0.4 | 1.8891 | 0.567 |

For `W0/dx >= 1.6` the critical radius is recovered to **< 1 %**. The two
finest cases are *not* window-in-`W` convergence: at `W0/dx < 1` the interface
is under-resolved and the measured value is meaningless. A clean
`W -> 0` sweep requires `dx ~ W0/2`, whose cost (roughly `W0^{-5}` in 2D for
fixed physical box and time) was not affordable in this pass; the resolved
points already show the thin-interface limit is reached and are the honest
statement of the fold.

The in-crate test uses `W0 = 0.4`, `dx = 0.4` on a `64^2` box, 9 bisection
iterations, and reports `R_c = 3.3164` vs `3.3333`, **0.51 %**.

### (c) The kinetic-coefficient-free condition

Two checks:

* The critical radius above is the Gibbs–Thomson value with no kinetic offset.
  With the KR choice `beta = a1 tau/(lambda W0) - a2 W0/D = 0` to leading
  order. Multiplying `tau0` by 5 switches on a large kinetic term but leaves
  `R_c` unchanged (the term is proportional to the velocity and vanishes at
  threshold); it changes the *growth rate*: for a seed at `1.5 R_c` the front
  moves `0.792` in the KR case and `1.452` with `5 tau0` (test
  `kinetic_choice_changes_the_interface`).
* Enthalpy is conserved: with insulated (symmetry) boundaries the discrete
  identity `dt(u - f_s) = D lap u` telescopes, and `sum(u - f_s)` drifts by
  `4.4e-15` over 500 steps on a `32^2` grid.

### (d) 2D free dendrite — resolved

Benchmark: `Delta = 0.55`, `eps4 = 0.05`, four-fold anisotropy, at the
Karma–Rappel Table II point `W0 = 1`, `D = 2`, `d0 = a1 a2 W0/D = 0.277`
(which makes `tau0 = 1`). The `dendrite` example (quarter domain, moving
frame, `N = 400 x 200`, `dx = 0.4`, `W0/dx = 2.5`) now uses the corrected
operator `PhaseField::step_kr` and gives a tip velocity that plateaus:

| `t` | tip `x` | `V` | `V d0/D` | `R/W0` | `V R^2/(2 D d0)` |
|---|---|---|---|---|---|
| 320 | 61.6 | 0.1416 | 0.01961 | 3.59 | 1.65 |
| 512 | 87.3 | 0.1303 | 0.01804 | 3.73 | 1.63 |
| 640 | 103.5 | 0.1255 | 0.01738 | 3.71 | 1.56 |
| 768 | 118.9 | 0.1178 | 0.01632 | 3.59 | 1.37 |

`V d0/D` settles near 0.017 (0.0163–0.0174 over `t = 640–768`) instead of
decaying; the slow drop after that is the finite-box moving-frame artifact
described in the GPU note. This reproduces the GPU f32 run on the same point
(`0.0164–0.0172`, KR 0.0168, Green function 0.0170) to within the
finite-box/run-to-run scatter, so the CPU and GPU operators agree. The
`R/W0` and `V R^2/(2 D d0)` columns in the table above are the *old* cap fit
(a `p = 0.5` parabola over a fixed 8 rows, `2.8 W0` at `dx = 0.4`), which sits
inside the diffuse-interface cap and biases `R` low by ~2.4x; the corrected
sharp-interface values are `R/W0 = 2 D p/V ~ 9` and `sigma_star =
(V d0/D)/(2 p^2) ~ 0.12`. See §3(ii) and
[solidification_cuda.md](solidification_cuda.md) §5.

A CPU run at the older thin-interface default (`W0/d0 = 0.4`, `W0/dx = 0.8`,
no moving frame) still shows a falling velocity, but that is a different,
under-resolved parameter point, not the KR benchmark. The earlier reported
miss — `V d0/D` passing through ~0.017 and then decaying while the `p = 0.5`
contour lost its parabolic tip — came from the plain `step` omitting the
variational anisotropy term of KR Eq. (88), not from the explicit scheme or
the interface width; with the term omitted the least-squares tip-radius fit
returns `NaN`. See §4 of [solidification_cuda.md](solidification_cuda.md).

## 3. The folds

### (i) Phase field -> sharp interface

The Gibbs–Thomson critical radius and the kinetic-free growth are the
thin-interface limit. The measured `R_c` agrees with `d0/Delta` to < 1 % for
resolved interfaces and converges as `W0 -> 0` with `dx`. The 1D enthalpy
solver is the sharp-interface Stefan reference that the phase field reduces
to; its front matches the exact Neumann solution to 0.5 %.

### (ii) Microstructure -> isobaric closure

When the dendrite tip is finer than the macroscopic cell, the isobaric model
of `analysis/compressible.md` §15 cannot represent it and needs a closure.
`closure.rs` provides:

* **Ivantsov + selection.** The 2D Ivantsov relation
  `Delta = sqrt(pi p) e^p erfc(sqrt(p))`, `p = v R/(2 D)` (Ivantsov 1947),
  fixes only `v R`. The microscopic-solvability condition
  `R^2 v = 2 D d0 / sigma_star` (Langer & Müller-Krumbhaar, *Acta Metall.*
  **26**, 1681 (1978)) supplies the missing relation, giving
  `v = 2 D sigma_star p^2 / d0`. `sigma_star` is a closure parameter; it is
  calibrated once against the phase-field / published tip and then reused
  across undercoolings. It is never tuned inside a validation.
* **Avrami mushy zone.** `X(t) = 1 - exp(-g N (v t)^d)` with `g = pi` (2D) or
  `4 pi/3` (3D), and its rate `dX/dt` for the latent-heat source
  `L dX/dt`. This is the same sub-grid model used for supercooled droplets in
  `analysis/compressible.md` §15, now with a tip-velocity-based growth rate.
* **Indicator.** `microstructure_indicator` compares the tip diameter
  `2 R = 2 d0/(sigma_star p)` against `cells_per_feature * dx_macro`. Below
  that the macroscopic cell averages the mush (use the closure); above it the
  cell can resolve the needle and should unfold to the phase field.

The dimensional closure (example `isobaric_closure`, `sigma_star = 0.1288`
from the quoted 2D benchmark at `Delta = 0.55`) gives for succinonitrile:

| `Delta T` (K) | `p` | `v` (m/s) | `R` (m) |
|---|---|---|---|
| 1.0 | 0.0006 | 4.0e-6 | 3.5e-5 |
| 2.0 | 0.0027 | 7.2e-5 | 8.2e-6 |
| 5.0 | 0.0203 | 4.1e-3 | 1.1e-6 |

At `Delta T = 2 K` this is `v ~ 72 um/s`, `R ~ 8 um`, of the same order as
Glicksman's SCN dendrite measurements (tens of um/s and ~10 um). At
`Delta T = 1 K` the 2D Ivantsov Péclet number is so small that `R` is
overestimated by roughly an order of magnitude; a 3D Ivantsov relation
(`Delta = p e^p E_1(p)`) would be the physical choice for bulk 3D dendrites,
and is a clear next step. For a 1 mm isobaric cell every case above is
`Closure` (tip-diameter/cell `~1e-2` or less), which is why the macroscopic
model cannot show dendrites.

The value `sigma_star = 0.1288` is fixed, not fitted. Karma & Rappel's 2D
benchmark is `Delta = 0.55`, `eps4 = 0.05`, and their sharp-interface
(Green-function) tip speed is `V d0/D = 0.0170` (Karma & Rappel, *Phys. Rev.
E* **57**, 4323 (1998), Table II; the phase-field column is `0.0168`). Their
table tabulates `V` only; it does not list `R` or `sigma_star` (their shape
figure compares the interface with the Ivantsov paraboloid in units of `d0`,
and the tip region has `R ~ 30-40 d0`). The operating point is Ivantsov to
within the measured scatter, so with the 2D Ivantsov Péclet number
`p(0.55) = 0.25693` the selection constant is
`sigma_star = (V d0/D)/(2 p^2) = 0.0170/(2*0.25693^2) = 0.1288`
(`closure::sigma_star_2d_kr_benchmark`). The corresponding sharp-interface tip
radius is `R = 2 D p/V = d0/(sigma_star p) = 30.2 d0`, i.e. `R/W0 = 8.4` at
the `D = 2` point and `5.6` at the `D = 3` point.

**Why an earlier pass measured `sigma_star ~ 0.5-0.7`.** The phase-field tip
radius was fitted from the `p = 0.5` contour over a fixed **8 rows**, which at
`dx/W0 = 0.4` is only `2.8 W0`. That window lies entirely inside the
diffuse-interface cap, where the local radius of curvature is a few `W0` and
has nothing to do with the macroscopic tip radius; it returns `R/W0 = 3.8`
(`D = 2`) and `2.95` (`D = 3`) instead of `8.4` and `5.6`. The bias is a factor
`2.2-2.4`, exactly the factor between `sigma_star = 0.13` and the reported
`0.5-0.7` (since `sigma_star` scales as `1/R^2`). Fitting the `p = 0` contour
over the single-valued monotone needle (this is now what
`PhaseField::tip_radius` does), excluding the first `1.5 W0` of the rounded cap
and capping at `40 W0`, gives `R/W0 = 9.3 +/- 0.5` for the steady `D = 2`,
`dx/W0 = 0.4` run and `5.46 +/- 0.4` for `D = 3`, `dx/W0 = 0.2`, matching
`2 D p/V` and `sigma_star = 0.12-0.13` (the residual is the finite `dx` and
finite `W0/d0` scatter, see §2(d)). **The closure's `sigma_star = 0.1288` is
therefore correct; the inconsistency reported earlier was in the tip-radius
estimator, not in the closure.** A whole-needle fit that extends beyond
`~40 W0` overshoots `R` by ~25 % because the selected shape is only
approximately parabolic away from the cap.

### Materials (`materials.rs`)

* **Succinonitrile** (the standard dendrite material): `T_m = 331.24 K`,
  `L = 4.55e7 J/m^3`, `c_p,vol = 1.97e6 J/(m^3 K)`, `gamma = 8.9e-3 J/m^2`,
  `alpha_l = 1.1e-7 m^2/s`, from Glicksman's SCN review and NASA microgravity
  dendrite experiments (Glicksman, Winsa, Hahn, Lograsso & Tirmizi, 1976–93).
  Derived: `Gamma = gamma T_m/L = 6.5e-8 K m`, `d0 = 2.82e-9 m`. Measured
  surface anisotropy `eps4 = 0.0055`.
* **Ice/water**: `T_m = 273.15 K`, `L = 3.34e8 J/m^3`, `gamma = 2.91e-2 J/m^2`
  (Hardy 1977; IAPWS), `c_p,vol = 1.93e6 J/(m^3 K)`, `alpha = 1.14e-6 m^2/s`
  (k = 2.2 W/(m K)). Derived `d0 = 1.4e-10 m`, the smallest capillary length
  here because water's volumetric latent heat is large. `eps4 = 0.05` is an
  *estimate*; real ice anisotropy is six-fold and face-dependent. For a 1 mm
  lake/casting cell the tip diameter/cell is `~5e-2`, so ice microstructure
  is always a closure at that resolution.

## 4. What is established, what is not

Established, with errors:

* Neumann Stefan: 0.5 % front, 0.03 K temperature (1D enthalpy, `n = 400`).
* Gibbs–Thomson `R_c = d0/Delta`: < 1 % for `W0/dx >= 1.6`.
* KR kinetic-free condition: no kinetic offset in `R_c`; `tau0` controls the
  growth rate; enthalpy conserved to round-off.
* 2D free dendrite (KR Table II): the corrected CPU operator gives a plateau
  `V d0/D = 0.0163–0.0174`, matching the GPU (`0.0164–0.0172`) and KR's 0.0168
  (Green function 0.0170) to ~1–2 %. The plain operator decays and loses the
  parabolic tip.

Not established / open:

* 3D Ivantsov closure and the 3D dendrite.
* The GPU `D = 3` plateau is 13 % below KR (0.0150 vs 0.0175), an unresolved
  deficit; `D = 4` was not run to steady state. The CPU `step_kr` is
  cell-for-cell the GPU operator, but repeating those long moving-frame runs on
  the CPU was out of scope here.
* A phase-field version of the cold-wall Stefan test (the enthalpy solver is
  the sharp-interface reference here rather than the phase field itself).
* Coupling the closure to the actual isobaric 2D/3D solver
  (`analysis/compressible.md` §15 leaves that solver one-dimensional).

## Sources

* G. Caginalp / J. B. Collins & H. Levine — early phase-field; not used.
* A. Karma & W.-J. Rappel, *Phys. Rev. E* **54**, 3065 (1996); **57**, 4323
  (1998) — thin-interface model and constants `a1, a2`.
* G. P. Ivantsov, Dokl. Akad. Nauk SSSR **58**, 567 (1947) — diffusion solution.
* J. S. Langer & H. Müller-Krumbhaar, *Acta Metall.* **26**, 1681 (1978) —
  microscopic solvability.
* H. S. Carslaw & J. C. Jaeger, *Conduction of Heat in Solids*, 2nd ed. §11.2 —
  Neumann solution.
* M. E. Glicksman et al., NASA microgravity dendritic-growth experiments and
  SCN reviews (1976–1993) — SCN properties and anisotropy.
* S. C. Hardy, *J. Cryst. Growth* **38**, 29 (1977) — ice/water interfacial
  free energy.
* T. Koop et al., *Nature* **406**, 611 (2000) — homogeneous ice nucleation
  (used only by the closure example, via `analysis/compressible.md` §15).

## GPU convergence (appended note)

The CUDA port `fluid/crates/solidification-cuda/` (see
[solidification_cuda.md](solidification_cuda.md)) resolves the §2(d) miss and
identifies its cause. The plain `solidification` operator
(`PhaseField::step`) omits the variational anisotropy term of Karma–Rappel
Eq. (88), `p_x d_y Q − p_y d_x Q` with `Q = W dW/dθ`; with it omitted the tip
velocity passes through ~0.017 and then
decays while the `p = 0.5` contour loses its parabolic tip (the radius fit
fails). Adding the term restores a stable tip; the 9-point Laplacian of KR
then lowers the velocity toward the published value. On the exact KR Table II
point (`W0 = 1`, `Δ = 0.55`, `ε4 = 0.05`, `D = 2`, `d0 = 0.277`,
`dx = 0.4`), GPU f32 gives a quasi-steady `V d0/D ≈ 0.0164–0.0172`
(sharp-interface `R = 2 D p/V ≈ 9 W0`), versus KR's 0.0168 (Green function
0.0170) — **~1–2 %**. The
same operator at `D = 3` plateaus at 0.0150 versus KR's 0.0175 (**13 % low**),
an unresolved deficit; `D = 4` was not run to steady state. The GPU f64 step
reproduces the CPU to `3e-16`, and GPU f32 runs **~2.7×** the multi-core CPU
f64 solver (GPU f64 is slower: consumer FP64 is throttled). The CPU crate now
exposes the same corrected operator as `PhaseField::step_kr`, and its
`dendrite` example on this point gives `V d0/D = 0.0163–0.0174` over
`t = 640–768`, matching the GPU table above (see §2(d)). The plain
`PhaseField::step` (no variational term) is retained for the GPU = CPU check
and keeps the old isotropic results bit-for-bit.

## GPU convergence: the D = 3 deficit resolved (appended note)

Follow-up work in `solidification-cuda` (`examples/dendrite_convergence.rs`,
new `SHIFT_AT`/`SHIFT_TO`/`DIAG`/`TAU_FACTOR` switches, gentle moving frame
and a generic `r_tip = 2 W0` parabolic needle seed) shows the 13 % D = 3
deficit reported above is **grid resolution**, and that the earlier D = 2
"agreement" was a transient. True steady states (tip speed and radius
constant for thousands of time units, far field exactly `u_far = −Δ`) at
`Δ = 0.55`, `ε4 = 0.05`, `W0 = 1`:

| `D` | `d0/W0` | `dx/W0` | `V d0/D` | `R/W0` | `σ* = 2Dd0/(R²V)` |
|---|---|---|---|---|---|
| 2 | 0.277 | 0.4 | 0.01536 | 9.27 | 0.116 |
| 2 | 0.277 | 0.3 | 0.01574 | 9.04 | 0.119 |
| 2 | 0.277 | 0.2 | 0.01663 | 8.56 | 0.126 |
| 3 | 0.1846 | 0.4 | 0.01490 | 6.37 | 0.113 |
| 3 | 0.1846 | 0.3 | 0.01536 | 6.18 | 0.116 |
| 3 | 0.1846 | 0.2 | 0.01600 | 5.93 | 0.121 |
| 4 | 0.1385 | 0.4 | 0.01450 | 4.91 | 0.110 |

`R` is the sharp-interface tip radius `2 D p/V` (`p = 0.25693`); because the
operating point is Ivantsov it is fixed by `V d0/D`, and the whole-needle
`p = 0` measurements at two of these points agree with it. The `R` values in
earlier versions of this note (`3.8`, `2.95`, `2.66`) were a `p = 0.5`
parabola fit over a fixed 8 rows (`2.8 W0` at `dx/W0 = 0.4`), i.e. the
diffuse-interface cap, and biased `R` low by ~2.4x. For D = 3 the velocity
rises `0.0149 → 0.0154 → 0.0160` as `dx/W0` falls `0.4 → 0.3 → 0.2` and
extrapolates to ≈0.0173, i.e. KR's phase-field 0.0175 and the Green-function
0.0170. The old D = 2 point at `t = 768` (`0.0164`) was still relaxing from
the seed; run to steady it converges to 0.01536 at `dx/W0 = 0.4` and 0.01663
at 0.2, so D = 2 agrees with KR only after refinement, not before. The
residual fixed-`dx` ordering (D = 2 > 3 > 4) is the finite-`W0/d0` kinetic
correction: `τ0` is still a spectator only asymptotically, and the phase
field approaches the sharp-interface result as `W0/d0 → 0`. Domain, moving
frame and boundary treatment are ruled out (diffusion length `2D/V` vs box;
`SHIFT = 0` control; lateral-size control), as are `f32` round-off and an
operator error. The `σ*` values (0.11–0.13) follow from the measured `(V, R)` and are
consistent with the `σ* = 0.1288` used to calibrate the sub-grid closure in
§3(ii); the earlier `0.5–0.7` was the tip-radius estimator sitting inside the
diffuse-interface cap, and the closure has been left at `0.1288` (now fixed
by `closure::sigma_star_2d_kr_benchmark`). Full details and reproducibility
are in [solidification_cuda.md](solidification_cuda.md) §5.
