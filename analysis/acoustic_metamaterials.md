# Acoustic metamaterials from fluid-filled resonators

Implemented in `fluid/crates/phononic/src/helmholtz.rs` (benchmarks in
`fluid/crates/phononic/tests/helmholtz.rs`). This is the Helmholtz-resonator
rung of the acoustics ladder in [phononic.md](phononic.md): a fluid cavity
behind a neck is a lumped mass–spring oscillator, and an array of them loads a
waveguide so strongly that it acquires a **negative effective bulk modulus**.
The organising idea is [emergence_ladder.md](emergence_ladder.md): a cheap
lumped model, a cheap validity indicator (length over wavelength), and an
exact duct transfer matrix to unfold to where the lumped model fails.

## The models

**Lumped Helmholtz resonator.** A cavity of volume `V` behind a neck of
cross-section `S = pi a^2` and geometric length `L` has acoustic mass
`M_a = rho L_eff / S` and acoustic compliance `C_a = V / B0`, with
`B0 = rho c^2`. The series impedance is
`Z = i omega M_a + 1/(i omega C_a)` and the resonance is

\[
  \omega_0 = c\sqrt{\frac{S}{V L_{\mathrm{eff}}}}.
\]

The effective neck length carries the fluid that moves outside the geometric
neck:

\[
  L_{\mathrm{eff}} = L + \delta_{\mathrm{mouth}} + \delta_{\mathrm{junction}}.
\]

For a circular neck the low-frequency end corrections are `delta/a = 8/(3 pi)
= 0.84883` for a flanged (baffled) opening and `delta/a = 0.6133` for an
unflanged (free-space) opening. The flanged value is the low-`ka` reactance of
a baffled piston; the unflanged value is Rayleigh's open-pipe correction
(Rayleigh, *Theory of Sound*, Vol. II; Levine & Schwinger, J. Acoust. Soc. Am.
**20**, 637 (1948); Kinsler, Frey, Coppens & Sanders, *Fundamentals of
Acoustics*, 4th ed., ch. 9). Both coefficients are returned exactly by
`EndCondition::coefficient` and checked in the tests.

**Thermoviscous narrow-tube losses.** The Zwikker–Kosten / Stinson
low-reduced-frequency (LRF) model (Zwikker & Kosten, *Sound Absorbing
Materials*, 1949; Stinson, J. Acoust. Soc. Am. **89**, 550 (1991)) gives a
complex density and compressibility for a circular tube of radius `a`:

\[
  \rho_c = \frac{\rho_0}{1-\dfrac{2J_1(\alpha a)}{\alpha a J_0(\alpha a)}},
  \qquad
  C_c = \frac{1}{B_0}\!\left[1+(\gamma-1)\frac{2J_1(\alpha_t a)}{\alpha_t a J_0(\alpha_t a)}\right],
\]
\[
  \alpha^2 = -i\omega\rho_0/\mu,\qquad \alpha_t = \sqrt{\mathrm{Pr}}\,\alpha .
\]

The characteristic impedance is `Z_c = sqrt(rho_c/C_c)`, the wavenumber
`k_c = omega sqrt(rho_c C_c)`, and the neck impedance is
`Z_neck = i omega rho_c L_eff / S`. This is what `lumped_impedance_lossy` and
`input_impedance_exact_lossy` use. In the Poiseuille limit it reduces to the
classical resistance `8 mu L/(pi a^4)`.

**Exact duct model.** [`Helmholtz::input_impedance_exact`] concatenates the
neck and cavity as transmission-line ducts and terminates the cavity rigidly;
the input impedance is `-m22/m21` of the product matrix. The same machinery,
with complex `Z_c, k_c`, gives the lossy LRF input impedance.

**Side-branch waveguide.** A dead-end branch of impedance `Z_b` on a guide of
impedance `Z0 = rho c / A` transmits the pressure amplitude

\[
  t = \frac{2 Z_b}{2 Z_b + Z_0},
\]

so at resonance (`Z_b -> 0`) the guide is short-circuited (`t -> 0`,
total reflection) and off resonance it is transparent. A periodic array,
spacing `a`, has the exact Bloch relation `cos(q a) = tr(M_cell)/2` with the
cell the duct segment times the shunt. At long wavelength it folds to an
effective medium with

\[
  \frac{1}{B_{\mathrm{eff}}}
   = \frac{1}{B_0}\!\left[1-\frac{F\,\omega_0^2}{\omega^2-\omega_0^2+i\Gamma\omega}\right],
  \qquad F = \frac{V}{A a},
\]

which is the Fang et al. form (Nature Materials **5**, 452 (2006)); the
filling factor `F = V/(Aa)` is fixed by the geometry, not fitted. `Re B_eff < 0`
on `(omega_0, omega_0 sqrt(1+F))`: the negative-bulk-modulus band.

*Sign convention.* A dead-end branch removes volume velocity from the
through-flow, `U_2 = U_1 - Y_b p`, so the shunt matrix has off-diagonal
`m21 = -Y_b`. With the opposite sign the folded medium reproduces the
mirror-image band *below* `omega_0`, which contradicts the exact transfer
matrix. This was found while benchmarking and is the reason the shunt matrix is
written with the minus sign.

**Elastic-walled cavity.** A spherical cavity in a solid of shear modulus `G`
bulges under internal pressure; the Lamé solution gives `u(R) = pR/(4G)`, so
`Delta V/V = 3p/(4G)` and the wall adds `3/(4G)` to the cavity compliance
`1/B0`. A long cylindrical cavity in plane strain adds `1/G`. The folded
cavity stiffness is `1/B_cav = 1/B0 + 1/B_wall`. The wall term is negligible
unless `B0/G` is not small — true for water in a polymer, irrelevant for air.

## Benchmarks (measured)

Run `cd fluid && cargo test -p phononic --test helmholtz`. Every value is
printed by the tests; "error" is the deviation from the exact or published
result.

| Benchmark | Exact / published | Computed | Error |
|---|---|---|---|
| Flanged end correction `delta/a` | `8/(3 pi) = 0.848826` | 0.848826 | round-off |
| Unflanged end correction `delta/a` | 0.6133 | 0.613300 | round-off |
| Lumped vs exact resonance, `L/lambda=0.016`, `Lc/lambda=0.043` | 416.56 Hz lumped | 410.88 Hz exact | 1.36% |
| Stinson neck resistance vs Poiseuille `8 mu L/(pi a^4)`, `s=0.32` | 3.71786e6 Pa s/m³ | 3.71789e6 | 9.3e-6 |
| Stinson wide tube `|Z_c-rho c|/rho c`, `s~100` | 0 | 5.17e-4 | — |
| Stinson wide tube `|Im k|/Re k` | 0 | 1.03e-3 | — |
| Lumped lossy resonance vs exact LRF (short neck) | 432.25 Hz | 436.95 Hz | 1.1% |
| Single branch at resonance `|t|` | 0 | 0 | round-off |
| Single branch off resonance `|t|` | 1 | 1.000000 | <1e-6 |
| Lossless branch energy `|r|²+|t|²`, `r=t-1` | 1 | 1 | <1e-9 |
| Hybridisation gap (exact Bloch), `omega_0/2pi=416.6 Hz` | — | 413.4–435.8 Hz | contains `omega_0` |
| Bragg gap (exact Bloch), `c/2a=1715.0 Hz` | — | 1715.0–1725.7 Hz | centre 0.3% high |
| Fang `B_eff` vs exact Bloch `q`, `ka=0.153` | equal | equal | 9.6e-6 |
| Fang `B_eff` vs exact Bloch `q`, `ka=0.305` | equal | equal | 4.95e-5 |
| Fang `B_eff` vs exact Bloch `q`, `ka=0.458` | equal | equal | 1.88e-4 |
| Fang `B_eff` vs exact Bloch `q`, `ka=0.611` | equal | equal | 9.74e-4 |
| Negative-B sign at `0.5 omega_0` / `0.5(omega_0+w_top)` | `B>0` / `B<0` | holds | exact |
| 20-cell TL dip inside Bloch gap | in gap | 416.6 Hz, 1295 dB | in gap |
| Water/steel wall-to-fluid compliance `3B0/(4G)` | — | 2.01e-2 | — |
| Water/PMMA wall-to-fluid compliance | — | 0.788 | — |
| Air/steel wall-to-fluid compliance | — | 1.30e-6 | — |
| Water/PMMA cavity `f0` shift, folding wall compliance | 1799.85 Hz rigid | 1346.11 Hz elastic | -25.2% |
| Water/air resonance ratio (same geometry) | `c_w/c_a = 4.3207` | 4.3207 | round-off |
| Double-negative overlap `B<0`, `rho<0`, `n<0`, `omega=1.05` | — | `B=-0.519, rho=-1.530, n=-1.717` | — |

### Where the folded model fails

Sweeping the cavity length with the neck fixed (so `L_c/lambda` grows and the
lumped resonance drops), the folded mass+compliance error grows with the
length-over-wavelength indicator:

| `L_c` [m] | lumped `f0` [Hz] | exact `f` [Hz] | `L_c/lambda` | relative error |
|---|---|---|---|---|
| 0.03 | 452.30 | 446.80 | 0.0396 | 1.22e-2 |
| 0.10 | 247.73 | 239.41 | 0.0722 | 3.36e-2 |
| 0.30 | 143.03 | 129.83 | 0.1251 | 9.23e-2 |
| 0.60 | 101.14 | 84.35 | 0.1769 | 1.66e-1 |
| 1.00 | 78.34 | 59.26 | 0.2284 | 2.44e-1 |
| 1.20 | 71.51 | 51.81 | 0.2502 | 2.76e-1 |

The exact cavity becomes a quarter-wave resonator as `L_c -> lambda/4`, which a
pure compliance cannot represent; the error reaches tens of percent. The
indicator to watch is `max(L_eff, L_c)/lambda`; the fold is safe below about
`lambda/20` (≈1% error here) and fails by `lambda/8`–`lambda/4`.

### Effective-medium validity

The Fang form with the geometric `F` matches the exact Bloch wavevector to
`9.6e-6` at `ka = 0.153` and the error grows quadratically with the
sub-wavelength indicator `ka = omega a/c` (to `9.7e-4` at `ka = 0.611`). This
is the analogous result to the laminate fold in [phononic.md](phononic.md):
the error is `O(ka^2)` at small `ka` (the fitted `err/ka^2` rises from
`4.1e-4` at `ka=0.15` to `2.6e-3` at `ka=0.61` as the hybridisation gap is
approached). The same effective medium
correctly places the negative-`B` band on `(omega_0, omega_0 sqrt(1+F))`, which
the exact array realises as its hybridisation gap. Unfold to the exact cell
matrix once `ka` is no longer small or the drive nears `omega_0`.

### Losses

The LRF model reproduces the Poiseuille resistance to `9e-6` at `s = 0.32`.
The inertial term carries an unsteady boundary-layer added mass: at `s = 0.32`
the imaginary part is `1.33x` the ideal `omega rho L/S`, a small correction to
the total impedance (the resistance is `~58x` larger). The wide-tube limit
returns `Z_c -> rho c` and `k -> omega/c` to `~1e-3`. Loss shifts the
resonance by `~1%` for a short neck and turns the lumped input impedance at
`omega_0` from zero into `5.2e4 + i5.1e4` Pa s/m³.

## Water vs air; wall compliance

Filling the same resonator with water instead of air multiplies `omega_0` by
`c_w/c_a = 4.32` and raises the bulk modulus by four orders of magnitude. The
wall-to-fluid compliance ratio `3 B0/(4G)` (spherical cavity) is:

| fluid / wall | `B0` [Pa] | `3B0/(4G)` | `B_elastic/B0` |
|---|---|---|---|
| air / steel | 1.416e5 | 1.30e-6 | 1.0000 |
| air / PMMA | 1.416e5 | 5.09e-5 | 0.9999 |
| water / steel | 2.192e9 | 2.01e-2 | 0.9803 |
| water / PMMA | 2.192e9 | 0.788 | 0.5594 |
| water / epoxy | 2.192e9 | 1.04 | 0.4913 |

Wall compliance is irrelevant for air in any solid, a 2% correction for water
in steel, and order-one for water in a polymer: a water/PMMA cavity's folded
resonance drops 25.2% (1799.85 -> 1346.11 Hz). The two effects are folded into
`cavity_bulk_modulus_elastic`.

*Approximations.* The wall model is quasi-static: it keeps only the breathing
(monopole) compliance of an infinite elastic medium and neglects wall inertia,
wall resonances, radiation from the wall, finite-thickness and finite-extent
effects, and viscoelastic damping. It also assumes the cavity pressure is
hydrostatic; a real cavity has a spatial pressure distribution that couples to
wall modes. The fluid is treated as an inviscid, lossless bulk for the
compliance; the neck carries the thermoviscous loss.

## Double-negative band (optional)

Combining a Helmholtz array (negative `B_eff`, band `(omega_0, omega_0
sqrt(1+F))`) with a membrane-loaded guide (negative `rho_eff` above the
membrane resonance) gives a double-negative medium, following Lee, Park, Seo,
Wang, Kim, Kim, Rho, Kim & Lee, Phys. Rev. Lett. **104**, 054301 (2010). At
`omega = 1.05 omega_0`, `B_eff = -0.519` and `rho_eff = -1.530` in scaled
units, and the effective index `n = -1.717`. The overlap is real but the
current membrane model is a single-DOF effective-density law with a
phenomenological damping `gamma`; the exact membrane–fluid coupling
(tension, mass per area, fluid loading) is not implemented.

## Fold/unfold map

| coarse model (fold) | cheap indicator | unfold when | where |
|---|---|---|---|
| Lumped Helmholtz mass + compliance | `max(L_eff, L_c)/lambda` | `lambda/20`–`lambda/8` | `input_impedance_exact` |
| Ideal (lossless) neck | `s = a sqrt(omega rho/mu)`, `ka` | boundary layers matter | Stinson LRF (`lumped_impedance_lossy`) |
| Fang effective bulk modulus (one number per frequency) | `ka = omega a/c` | `ka` no longer small, or near `omega_0` | `cell_matrix`, `bloch_wavevector` |
| Rigid-wall cavity compliance | `3 B0/(4G)` | comparable to 1 | `cavity_bulk_modulus_elastic` |
| Single-DOF membrane density | `omega/omega_m` | near the membrane resonance | (not implemented) |

Hand-off in practice:

- Fold the resonator to `M_a, C_a` while both the neck and the cavity are
  shorter than about `lambda/20`; unfold to the duct matrix when the indicator
  rises, and never use the compliance past `lambda/8`.
- Fold the array to `B_eff` while `ka` is small; the measured prefactor is
  `~4e-4` per `ka^2`, far below the laminate envelope, so `ka < 0.3` is already
  a `1e-4` phase error. Unfold past the hybridization gap.
- Fold the elastic wall only when `3 B0/(4G)` is not small; for stiff metal
  walls and soft fluids it is a sub-percent correction.

## Limits and honest misses

* The exact model is 1D (plane waves, pressure uniform across the duct). It
  ignores higher-order duct modes, which matter above the first cut-off
  `k a ~ 1.84` (the first non-planar mode of a circular duct).
* The Stinson Bessel model is the LRF approximation: pressure uniform across
  the tube. It is benchmarked against the Poiseuille resistance and the
  lossless wide-tube limit, not against a full 3D thermoviscous solution.
* The side-branch array uses lossless branch impedances in the exact Bloch
  calculation; the effective medium carries a phenomenological `Gamma`. A
  thermoviscous `Gamma` derived from the Stinson resistance is not implemented.
* The Bragg gap is only approximately at `c/2a` because the branches load the
  guide; at the tested loading the centre is 0.3% high.
* The wall-compliance model is static and single-shape (spherical/cylindrical
  cavity in an infinite medium); finite geometry, wall inertia and wall modes
  are neglected.
* The double-negative membrane model is a single-DOF effective density, not a
  membrane–fluid coupled eigenproblem.
* No nonlinear acoustics, no flow, no temperature gradients, no radiation
  damping at the neck mouth beyond the low-`ka` end correction.

## References

* J. W. Strutt (Lord Rayleigh), *The Theory of Sound*, Vol. II (open-pipe end
  correction `0.6133 a`; baffled-piston reactance `8/(3 pi) ka`).
* H. Levine and J. Schwinger, "On the radiation of sound from an unflanged
  circular pipe", Phys. Rev. **73**, 383 (1948); J. Acoust. Soc. Am. **20**,
  637 (1948) (unflanged end correction).
* L. E. Kinsler, A. R. Frey, A. B. Coppens and J. V. Sanders, *Fundamentals of
  Acoustics*, 4th ed., ch. 9 (Helmholtz resonator, end corrections).
* C. Zwikker and C. W. Kosten, *Sound Absorbing Materials* (Elsevier, 1949);
  M. R. Stinson, "The propagation of plane sound waves in narrow and wide
  circular tubes, and generalization of the bulk viscosity coefficient",
  J. Acoust. Soc. Am. **89**, 550 (1991) (thermoviscous narrow-tube model).
* N. Fang, D. Xi, J. Xu, M. Ambati, W. Srituravanich, C. Sun and X. Zhang,
  "Ultrasonic metamaterials with negative modulus", Nature Materials **5**, 452
  (2006) (negative bulk modulus, filling factor).
* S. H. Lee, C. M. Park, Y. M. Seo, Z. G. Wang and C. K. Kim, "Acoustic
  metamaterial with negative density", Phys. Lett. A **350**, 146 (2005);
  S. H. Lee et al., "Composite acoustic medium with simultaneously negative
  density and bulk modulus", Phys. Rev. Lett. **104**, 054301 (2010)
  (double-negative acoustic medium).
* *CRC Handbook of Chemistry and Physics* (air and water density, sound speed,
  viscosity, Prandtl number at 20 C).
