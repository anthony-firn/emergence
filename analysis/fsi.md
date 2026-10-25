# Fluid–solid coupling against exact layered acoustics

Crate: `fluid/crates/compressible`. Everything below comes from
`examples/fsi_layered.rs` (run `cargo run --release -p compressible --example
fsi_layered -- all`) and is asserted by `tests/fsi_layered.rs`. The exact
reference is `fluid/crates/phononic` (`layered::Stack` and `homogenize`).

The organising idea is in `emergence_ladder.md`: the compressible solver is the
resolved model; the exact transfer matrix is the closure that a sub-wavelength
laminate should fold onto. This note measures how well the solver's
fluid–solid interface reproduces the closure, and where it does not.

## 0. What the crate represents, and the missing rung

* **Solids** are the GPR unified-continuum model (`gpr::Gpr`): a
  *single-material* `System<14>` with one NASG EOS and one shear speed `c_s`.
  Its linear longitudinal speed is `c_L = sqrt(c_0^2 + 4 c_s^2/3)`.
* **Multi-material flow** is the Kapila/Allaire five-equation model
  (`multiphase::FiveEq`): a `System<7>` with two partial densities and a
  volume fraction, one NASG EOS per phase, a **diffuse interface** and Wood's
  (or Allaire's) mixture sound speed.

There is **no shear-bearing multi-material system**; `threephase.rs` states
that GPR elasticity per phase "is a separate rung". To test solid–fluid
coupling directly we implemented that rung in the example file as `GprMix`
(GPR distortion carried with the Kapila partial densities). It is correct for
single-material elastic waves but **unstable at the water/steel contrast**
(§4). The quantitative validation therefore uses the reduction that is exact
at normal incidence:

> At normal incidence the shear degrees of freedom decouple, and a solid's
> P-wave reflection is fixed by the constrained modulus `M = K + 4G/3`, i.e.
> by `(rho, c_L)`. A two-phase stiffened-gas fluid with the same `(rho, c_L)`
> has the same impedance `rho c_L` and speed, so `FiveEq` with steel at `c_L`
> reproduces the longitudinal solid response exactly.

`phononic` is given the same `rho` and `c_L = sqrt((K+4G/3)/rho)`; the shear
speed enters only through `c_L`. Material constants (water `998.2 kg/m3`,
`1482 m/s`; steel `7850`, `5960`; aluminium `2700`, `6420`) are from
`phononic::materials`, which carries the handbook sources.

## 1. Single plate: a pulse in water onto 20 mm of steel

Setup: 1D, water, a 20 mm steel plate at `x = 0.25 m`, a Gaussian-modulated
80 kHz pressure pulse (amplitude 100 Pa, ~1e-7 Mach) launched at `x = 0.1 m`.
Probes sit half a cell in front of and behind the plate; the exact reference
stack includes those water spacers, so no propagation de-phasing is fitted.
The incident is measured by an identical run with no plate, so the reflected
spectrum is `(p_with - p_without)/p_without` and the transmitted is
`p_with(probe right)/p_without(probe left)`.

Exact values at 40/80 kHz: `|R| = 0.996/0.998`, `|T| = 0.0870/0.0637`. The
transmitted amplitude is small, so the meaningful accuracy is the **absolute**
incident-weighted error `sqrt(<(|T|meas-|T|exact)^2>)`.

| N | steps | `\|T\|` abs err (20–250 kHz) | `\|R\|` rel err | `\|T\|` at 40 kHz | `\|T\|` at 80 kHz | phase rms (R/T) |
|---|---|---|---|---|---|---|
| 256 | 1653 | 7.77e-1 | 2.96 | 0.246 (exact 0.087) | 0.537 (0.064) | 2.77 / 2.38 rad |
| 512 | 3306 | 2.98e-1 | 1.11 | 0.149 | 0.231 | 1.54 / 2.52 |
| 1024 | 6612 | **6.18e-2** | 0.519 | 0.0955 | 0.0959 | 1.40 / 1.61 |

The absolute `|T|` error falls `0.777 -> 0.298 -> 0.0618`, i.e. order
`log2(0.78/0.30) = 1.4` then `log2(0.30/0.062) = 2.3`: between first and
second order, tending to second. At `N = 1024` the 40 kHz transmission is
within 10% of exact and the 80 kHz phase is within about 1.6 rad (the phase
is dominated by the residual delay of the diffuse interface).

**Error source.** The interface is diffuse over a fixed number of cells
(`1.5 Delta x`). At `N = 256` the 20 mm plate is only ~9 cells thick, so the
two interfaces overlap and the plate never reaches its full impedance; the
computed plate is nearly transparent (`|T| = 0.25` vs exact 0.087). As the
plate gains cells the impedance contrast is recovered and the error falls.
The remaining error is the EOS mixing inside the diffuse zone; the exact
reference is a sharp interface.

**Half-space check.** `FiveEq` on a single water/steel interface gives a peak
reflection ratio `1.000` against the Fresnel value `(Z2-Z1)/(Z2+Z1) = 0.939`
at both `N = 256` and `N = 512`, so the interface does carry the correct
impedance when it does not overlap a second interface.

## 2. Periodic laminate: the homogenisation fold

Setup: 8 aluminium/steel quarter-wave periods at 400 kHz in a steel matrix
(phase 0 aluminium `6420 m/s`, phase 1 steel `5960 m/s`), a broadband 350 kHz
pulse. The exact first Bragg gap from `phononic::quarter_wave_band_edges` is
**278–522 kHz**. Reuss effective medium: `rho_eff = 5179 kg/m3`,
`c_eff = 5498 m/s`.

`N = 1024` (each layer ~6–7 cells); measured vs exact `|T|` band error 0.308,
homogenised vs exact 0.984.

| f (kHz) | `\|T\|` meas | `\|T\|` exact | `\|T\|` homog | indicator `(omega Lambda/c_eff)^2/12` |
|---|---|---|---|---|
| 150 (below gap) | 0.960 | 0.838 | 0.902 | 0.147 |
| 250 (edge) | 0.763 | 0.802 | 0.903 | 0.407 |
| 350 (gap) | 0.014 | 0.0013 | 0.983 | 0.798 |
| 450 (gap) | 0.0046 | 0.0013 | 0.979 | 1.32 |
| 650 (pass band) | 0.796 | 0.838 | 0.905 | 2.75 |

The fold is visible exactly as the ladder predicts:

* **below the gap** the indicator is 0.15, and the homogenised medium is
  within 0.06 of the exact transmission (0.902 vs 0.838), while the solver is
  within 0.12 (0.960 vs 0.838);
* **inside the gap** the indicator is O(1), the exact transmission collapses
  to `1.3e-3`, the homogenised medium still predicts `0.98`, and the solver
  follows the exact one (`1.4e-2`).

`tests/fsi_layered.rs::laminate_fold_and_gap` asserts these four inequalities,
so the fold cannot silently regress.

## 3. Oblique 2D case

**Not done.** The only shear-bearing system is single-material, so an oblique
fluid–solid test would need the multi-material GPR rung extended to 2D (with
the shear potential), which is outside this task. The normal-incidence
reduction above is the resolved model that is currently available.

## 4. Negative result: the shear-bearing two-material extension

`GprMix` in `examples/fsi_layered.rs` carries the GPR distortion together with
the Kapila partial densities and the volume fraction, with the shear speed
Voigt-averaged. Against single-material references it is correct:

| check | result |
|---|---|
| pure-steel longitudinal pulse speed, N = 800 | 5971 m/s vs `c_L = 5960` (+0.18%) |
| ideal-gas contact (all `p_inf = 0`), 500 steps | max `\|dp\|/p = 1.8e-14` |
| two-phase stiffened-gas contact at rest | `\|dp\|/p = 1.6e-8` after 200 steps, growing |

It fails at the water/steel contrast in two ways.

1. **Negative mixture sound speed.** The Kapila coefficient
   `K = alpha1 alpha2 (rho2 c2^2 - rho1 c1^2) / (...)` times
   `d(rho e)/d(alpha)` makes `rho c^2 < 0` once a little stiff material is
   mixed into water (e.g. `alpha_steel = 0.1`), because the two NASG
   `p_inf` differ by ~1e10. The solver's clamp hides the sign but the signal
   speeds are wrong. Switching to the Allaire closure (`K = 0`, as the
   validation above uses) removes the negative speed but not the next point.
2. **The diffuse interface is linearly unstable.** A stationary smooth
   water/steel interface, with no pulse, grows from `|dp|/p = 7.5e-10` to
   `1.6e-8` over 200 steps of the non-dissipative central path; extrapolated
   over a plate run this reaches O(1) and dominates the reflection. The
   over-reflection grows with contrast:

   | second medium | peak `\|r\|` | Fresnel |
   |---|---|---|
   | 2000 kg/m3, 2000 m/s | 0.32 | 0.46 |
   | 4000, 3000 | 1.69 | 0.78 |
   | 7850, 5960 | 4.05 | 0.94 |

   Forcing the shock path (`FluxMode::Weno`) does not change it, because the
   smooth interface never trips the sensor; the growth is intrinsic to the
   diffuse-interface energy/EOS mixing with a large `p_inf`.

The honest reading: **the crate cannot currently couple a shear-bearing solid
to a fluid at the steel contrast.** The five-equation model is stable on the
same contrast and is used for the quantitative results above; a stable
multi-material GPR would need a pressure-equilibrium-preserving or
internal-energy formulation rather than the total-energy one.

## 5. Limits

* All results are **normal incidence**. The solid is represented by its
  constrained modulus, so shear, mode conversion, coincidence and Lamb waves
  are not tested (the phononic fluid–plate model exists; the solver path does
  not).
* The interface is diffuse over `1.5 Delta x`; its physical width therefore
  shrinks with the grid. The exact reference is a sharp interface, so the
  errors above are the diffuse-interface/EOS-mixing error, not a grid
  truncation error of the smooth solver.
* The laminate demonstration is aluminium/steel in a **steel matrix**, and the
  plate is steel in water; the two-phase model cannot hold three constituents,
  so a steel/epoxy/water stack is not available.
* The example default is `FSI_N=256,512` (~40 s); `FSI_N=1024,2048` gives the
  converged runs (up to ~2 min at 2048, release only).
