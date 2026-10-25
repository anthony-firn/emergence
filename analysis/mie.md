# Mie scattering: the exact closure for a spherical particle, and its folds

Implemented in `fluid/crates/photonics/src/mie.rs`, integration tests in
`fluid/crates/photonics/tests/mie.rs`. This rung extends the wave-optics crate
(`analysis/photonics.md`) from planar and periodic matter to a single sphere of
arbitrary size, which is exactly the closure the geometrical-optics `radiation`
crate needs when a participating medium is made of particles, droplets or
bubbles. Nothing in `radiation` is changed; the two share the convention
`N = n + i kappa` (positive imaginary part = absorption) and SI units.

The organising idea is `analysis/emergence_ladder.md`: a cheap model (Rayleigh or
geometrical optics), a cheap validity indicator (`x`, `|m|x`, the phase delay
`rho`), and an unfold to the exact Mie series where the indicator fails.

## The exact model

For a sphere of radius `a`, relative index `m = m_particle / n_medium`, and
ambient index `n_medium`, the size parameter is

\[
  x = \frac{2\pi n_{\rm medium} a}{\lambda},
\]

and the Mie coefficients are (Bohren & Huffman 1983, eqs. 4.52, 4.88)

\[
  a_n = \frac{\left(\frac{D_n(mx)}{m}+\frac{n}{x}\right)\psi_n(x)-\psi_{n-1}(x)}
             {\left(\frac{D_n(mx)}{m}+\frac{n}{x}\right)\xi_n(x)-\xi_{n-1}(x)},\qquad
  b_n = \frac{\left(mD_n(mx)+\frac{n}{x}\right)\psi_n(x)-\psi_{n-1}(x)}
             {\left(mD_n(mx)+\frac{n}{x}\right)\xi_n(x)-\xi_{n-1}(x)},
\]

with the Riccati–Bessel functions `psi_n = x j_n`, `chi_n = -x y_n`,
`xi_n = psi_n + i chi_n`, and the logarithmic derivative
`D_n(z) = psi_n'(z)/psi_n(z)`. The integrals are

\[
  Q_{\rm ext} = \frac{2}{x^2}\sum (2n+1)\Re(a_n+b_n),\qquad
  Q_{\rm sca} = \frac{2}{x^2}\sum (2n+1)(|a_n|^2+|b_n|^2),\qquad
  Q_{\rm abs} = Q_{\rm ext}-Q_{\rm sca},
\]

\[
  Q_{\rm back} = \frac{1}{x^2}\left|\sum (2n+1)(-1)^n(a_n-b_n)\right|^2,
\]

\[
  g = \frac{4}{x^2Q_{\rm sca}}\sum\left[
      \frac{n(n+2)}{n+1}\Re(a_na^*_{n+1}+b_nb^*_{n+1})
      +\frac{2n+1}{n(n+1)}\Re(a_nb^*_n)\right].
\]

The angular amplitudes `S1`, `S2` follow from the `pi_n`, `tau_n` recurrence
(Bohren & Huffman eq. 4.47). The unpolarised phase function

\[
  p(\theta) = \frac{2(|S_1|^2+|S_2|^2)}{x^2 Q_{\rm sca}}
\]

is normalised to `(1/4 pi) integral p dOmega = 1`, and the unpolarised
differential cross section is `d sigma/dOmega = (|S1|^2+|S2|^2)/(2 k^2)`.
`MieScattering::phase_function` returns `p`; `s1_s2` returns the amplitudes.

### Numerical stability (Wiscombe 1980)

* `D_n(mx)` by the **downward** recurrence, seeded by the Lentz continued
  fraction. Downward is the minimal solution and does not overflow for complex
  arguments; upward would be contaminated once `n` passes `|mx|`.
* `psi_n(x)` by **Miller's downward** recurrence with rescaling, then normalised
  against the analytic `psi_0 = sin x`, `psi_1 = sin x/x - cos x`.
* `chi_n(x)` by the **upward** recurrence, which is the growing solution.
* The series is truncated at `n_stop = x + 4.05 x^{1/3} + 2` terms, giving about
  `1e-6` accuracy in the summed efficiencies.

The `x = 10^4` cases below (about `10^4` terms) and the strongly absorbing
`m = 10-10i` cases exercise all of this.

### Perfect-conductor limit

`MieScattering::perfect_conductor(x)` implements the `|m| -> infinity` limit
(`a_n = (n/x psi_n - psi_{n-1})/(n/x xi_n - xi_{n-1})`, `b_n = psi_n/xi_n`).
This is *not* the same as `new(x, C::new(0, 100))`: for a finite index the
general series differs from the limit by `O(1/|m|)`. At `x = 100` the finite
`m = 100i` value is `Q_sca = 2.01379` while the conductor limit is `2.008102`
(the `m = 100i` value was confirmed against a 50-digit mpmath evaluation of the
general formula, `2.013792663`). This is the reason the published MIEV0
"perfectly conducting" cases are checked through the dedicated constructor.

## Benchmarks

### Published Mie efficiencies

The following are the Wiscombe MIEV0 test cases (Wiscombe 1979, NCAR/TN-140+STR)
and the Bohren & Huffman (1983) Appendix test sphere (`m = 1.55`,
`x = 2 pi * 0.525/0.6328 = 5.212819668567135`). They are the same values used by
the `miepython` reference suite. "err" is the absolute deviation from the
published value.

| m | x | Q_sca published | Q_sca computed | err | g published | err g |
|---|---|---|---|---|---|---|
| 1.55 | 5.21282 | 3.1054300 | 3.1054255 | 4.5e-6 | 0.633140 | 3.2e-6 |
| 1.55−0.1i | 5.21282 | 1.6642491 | 1.6642491 | 2.0e-12 | 0.801290 | 4.6e-12 |
| 0.75 | 0.099 | 7.4e-6 | 7.4e-6 | 1.8e-8 | 0.001448 | 2.3e-7 |
| 0.75 | 10 | 2.2322650 | 2.2322648 | 1.6e-7 | 0.896473 | 4.5e-7 |
| 0.75 | 1000 | 1.9979080 | 1.9979082 | 1.8e-7 | 0.844944 | 2.9e-7 |
| 1.33−1e-5i | 100 | 2.0965940 | 2.0965935 | 4.9e-7 | 0.868959 | 2.7e-7 |
| 1.33−1e-5i | 10000 | 1.7238570 | 1.7238572 | 2.2e-7 | 0.907840 | 3.7e-7 |
| 1.5−1i | 1 | 0.6634538 | 0.6634538 | 3.8e-8 | 0.192136 | 4.0e-7 |
| 1.5−1i | 100 | 1.2836970 | 1.2836970 | 4.9e-8 | 0.850252 | 2.3e-9 |
| 1.5−1i | 10000 | 1.2365750 | 1.2365743 | 6.9e-7 | 0.846309 | 9.6e-7 |
| 10−10i | 1 | 2.0494050 | 2.0494050 | 6.9e-9 | −0.110664 | 3.6e-7 |
| 10−10i | 100 | 1.8367850 | 1.8367854 | 4.0e-7 | 0.556215 | 4.8e-7 |
| 10−10i | 10000 | 1.7953930 | 1.7953930 | 3.0e-8 | 0.548194 | 3.9e-8 |

The published values are quoted to 6 significant figures, so most errors are at
the rounding floor; the largest is `4.5e-6` on the `m = 1.55` case, whose
published `Q_ext` is rounded from `3.1054255`. The absorption convention here is
`N = n + i kappa`; the `m = 1.55−0.1i` published entry is reproduced by passing
`C::new(1.55, 0.1)`.

The backscatter efficiency is checked separately:

| m | x | Q_back published | computed | err |
|---|---|---|---|---|
| 1.55 | 5.21282 | 2.92534 | 2.92421 (at rounded x=5.213) / 2.92534 (exact x) | <1e-5 |
| 1.55−0.1i | 5.21282 | 0.20599534 | 0.20599534 | <1e-8 |

The Bohren & Huffman `m = 1.55` table is tabulated at the exact wavelength and
radius `lambda0 = 0.6328`, `radius = 0.525`, not at the rounded `x = 5.213`; using
the rounded value moves `Q_ext` by `4e-4` and the test uses the exact ratio.

### Perfect conductor (MIEV0)

| x | Q_sca published | computed | g published | err g |
|---|---|---|---|---|
| 0.101 | 0.000348 | 0.000348 | −0.397262 | <1e-6 |
| 100 | 2.008102 | 2.008102 | 0.500926 | <1e-6 |
| 10000 | 2.000289 | 2.000289 | 0.500070 | <1e-6 |

### Angular amplitudes

`m = 1.5−1i`, `x = 1`, at `theta = 0, 30, ..., 180°` (the `miepython`
"wiscombe" normalisation, equal to Bohren & Huffman `S1`, `S2`). Example rows:

| theta | Re S1 | Im S1 | Re S2 | Im S2 |
|---|---|---|---|---|
| 0° | 0.584080 | 0.190515 | 0.584080 | 0.190515 |
| 60° | 0.517525 | 0.178443 | 0.287964 | 0.041054 |
| 90° | 0.456340 | 0.167167 | 0.0362285 | −0.0618265 |
| 180° | 0.348844 | 0.146829 | −0.348844 | −0.146829 |

Matched to `<1e-5`. (This test also checks the sign convention: the crate's
positive-imaginary absorption is conjugated internally to the `N = n - i kappa`
Bohren–Huffman kernel.)

### Exact identities (errors are round-off)

* **Optical theorem**: `Q_ext = (4/x^2) Re S(0)` for real and absorbing spheres,
  matched to `<1e-12` relative over `x = 2..50`, `m = 0.75..1.5`.
* **Energy conservation**: `Q_ext = Q_sca` and `Q_abs = 0` to `<1e-11` for real
  `m = 0.75, 1.05, 1.2, 1.33, 1.5, 2.0` and `x = 0.1..500`.
* **Phase-function normalisation**: `(1/4 pi) integral p dOmega = 1` to `<1e-6`
  for `m = 1.5−0.1i`, `x = 5`.
* **Perfect conductor**: `Q_ext = Q_sca` to round-off (no absorption).

## Folds

### Rayleigh limit (`x << 1`, `|m|x << 1`)

\[
  Q_{\rm sca} = \frac{8}{3}x^4\left|\frac{m^2-1}{m^2+2}\right|^2,\qquad
  p(\theta) = \frac{3}{4}(1+\cos^2\theta),
\]

and for an absorbing particle `Q_ext = 4x Im[(m^2-1)/(m^2+2)] + Q_sca`. The
**validity indicator** is `max(x, |m| x)` (`rayleigh_indicator`); the relative
error of the scattering fold is `O(x^2)`. Measured relative error of `Q_sca`:

| x | m = 1.05 | m = 1.33 | m = 1.5 |
|---|---|---|---|
| 0.5 | 9.089e-2 | 2.416e-2 | 1.024e-2 |
| 0.25 | 2.195e-2 | 4.941e-3 | 3.977e-3 |
| 0.125 | 5.440e-3 | 1.171e-3 | 1.076e-3 |
| 0.0625 | 1.357e-3 | 2.887e-4 | 2.741e-4 |
| 0.03125 | 3.391e-4 | 7.194e-5 | 6.883e-5 |

Halving `x` cuts the relative error by 4 (clean `O(x^2)`), and the leading
coefficient grows with contrast, as expected. The dipole phase function is
matched to `<1e-4` for `x = 1e-3`.

### Large-particle limit (`x >> 1`)

The extinction paradox, `Q_ext -> 2` (a large opaque sphere removes twice its
geometric shadow because of diffraction). Measured for `m = 1.5`:

| x | Q_ext | Q_ext − 2 |
|---|---|---|
| 10 | 2.88199895 | 8.820e-1 |
| 100 | 2.09438781 | 9.439e-2 |
| 1000 | 2.01394465 | 1.394e-2 |
| 10000 | 2.00461747 | 4.617e-3 |
| 100000 | 2.00094201 | 9.420e-4 |

The deviation falls as about `1/x` after the oscillations average out; the
indicator for the fold is `1/x`.

### Anomalous diffraction (intermediate fold, `x >> 1`, `|m-1| << 1`)

van de Hulst's anomalous-diffraction extinction for a non-absorbing sphere,

\[
  Q_{\rm ext} = 2 - \frac{4}{\rho}\sin\rho + \frac{4}{\rho^2}(1-\cos\rho),
  \qquad \rho = 2x|m-1|,
\]

keeps the phase delay but ignores refraction. Indicator `rho`: the fold is
accurate once `rho >> 1`. At fixed `m = 1.05` the absolute error against the
exact series is oscillatory in `rho`:

| x | rho | \|err Q_ext\| |
|---|---|---|
| 10 | 1.0 | 1.062e-2 |
| 20 | 2.0 | 6.723e-2 |
| 50 | 5.0 | 1.309e-1 |
| 100 | 10.0 | 8.835e-2 |
| 200 | 20.0 | 4.394e-2 |
| 400 | 40.0 | 2.821e-2 |

The peak near `rho ~ 5` is the first diffraction oscillation; beyond it the
envelope falls roughly as `1/x`. The cleaner convergence statement is to hold
`rho` fixed while `x` grows and `m -> 1`, which is the fold's asymptotic regime:

| x | m = 1 + 5/x | \|err Q_ext\| |
|---|---|---|
| 10 | 1.5 | 5.908e-1 |
| 20 | 1.25 | 3.769e-1 |
| 40 | 1.125 | 1.765e-1 |
| 80 | 1.0625 | 1.065e-1 |
| 160 | 1.03125 | 5.915e-2 |
| 320 | 1.015625 | 3.193e-2 |

The error is cleanly `O(1/x)` (a factor ~1.85 per doubling), as expected.

### Fold map

| Coarse model (fold) | Cheap indicator | Unfold when | Where |
|---|---|---|---|
| Rayleigh dipole `Q_sca = (8/3)x^4\|K\|^2` | `max(x, \|m\|x)` | indicator exceeds the required accuracy (relative error `~indicator^2`) | `mie::rayleigh_qsca` |
| Rayleigh phase `(3/4)(1+cos^2)` | `max(x, \|m\|x)` | same | `mie::rayleigh_phase` |
| Anomalous diffraction (van de Hulst) | phase delay `rho = 2x\|m-1\|` | `rho` not large, or `\|m-1\|` not small | `mie::anomalous_diffraction_qext_real` |
| Extinction paradox `Q_ext -> 2` | `1/x` | `1/x` above the tolerance | large-`x` Mie |
| Henyey–Greenstein phase with matched `g` | deviation of `p(theta)` from HG, or number of scatterings | phase shape matters (small optical depth, strong forward peak) | exact `s1_s2` / `phase_function` |
| Perfect conductor `\|m\| -> inf` | `1/\|m\|` | finite index matters | `perfect_conductor` vs `new` |

## Applications that connect the repo

### (a) Stratocumulus-like water cloud -> `radiation` participating medium

A gamma size distribution with effective radius `r_eff = 10 µm` and effective
variance `v_eff = 0.1` (Hansen & Travis 1974), number density `N = 100 cm^-3`, at
`lambda = 550 nm`. Water index from the crate's Daimon–Masumura fit
(`dispersion::Sellmeier::water`), `n = 1.3347`. The distribution has
`alpha = 1/v_eff - 3 = 7` and `beta = 1/(v_eff r_eff)`; the integration is
composite Simpson over `[0, r_eff(1+8 sqrt(v_eff))]`.

| quantity | value |
|---|---|
| liquid water content `LWC` | 0.3016 g/m^3 |
| `sigma_ext` | 0.04725 m^-1 (47.25 km^-1) |
| `sigma_sca` | 0.04725 m^-1 |
| `sigma_abs` | 2.31e-19 m^-1 |
| single-scattering albedo `omega` | 1.0000000 |
| asymmetry `g` | 0.8638 |
| vertical optical depth over 300 m | 14.2 |
| large-particle law `3 LWC/(2 rho_w r_eff)` | 0.04524 m^-1 |

The large-`x` result is within 4.4% of the `Q_ext -> 2` law, as it should be at
`x ~ 150`; the residual is the real `Q_ext` oscillating about 2. For a typical
stratocumulus (`LWC ~ 0.3 g/m^3`, `r_eff ~ 10 µm`, `N ~ 100-300 cm^-3`) the
extinction coefficient is `~30-60 km^-1` (e.g. Stephens 1978; Petty 2006) and
`g ~ 0.85`, so the computed `47.25 km^-1`, `g = 0.864` and conservative albedo
agree with the textbook range. The optical constants are the crate's own, not a
tuned number.

**Exactly the `radiation` inputs.** `radiation::Medium` stores `absorption`,
`scattering` (both per metre), `g`, and `n`, and its
`coefficients(lambda, scale)` returns `(sigma_a * scale, sigma_s * scale)`
where `scale` is a volume fraction (so a cell holding droplet fraction `f`
carries `f sigma`). `ParticipatingMedium` therefore maps as

```text
Medium { absorption = sigma_abs, scattering = sigma_sca, g = g, n = n_medium }
```

`ParticipatingMedium::albedo()` is the single-scattering albedo the transport
solver uses; `g` is the asymmetry the Henyey–Greenstein fold consumes. The exact
`p(theta)` from `phase_function` is available when the HG shape is not good
enough (the unfold).

### (b) Optical forces on a 1 µm silica bead in water (1064 nm)

Exact Mie `Q_pr = Q_ext - g Q_sca` for a silica bead (`n = 1.4496` from
Malitson 1965) in water (`n = 1.3245` from Daimon–Masumura), `m_rel = 1.0944`,
`x = 3.911`:

| quantity | value |
|---|---|
| `Q_ext = Q_sca` | 0.2533 |
| `g` | 0.8681 |
| `Q_pr` | 0.0334 |
| scattering force at `P = 100 mW` | 12.8 pN |
| Rayleigh gradient force (max slope) | 219 pN |
| Rayleigh stiffness `dF/dr` | 0.875 pN/nm |
| potential depth `(1/2) k w0^2` | 2.6e4 k_B T |

Assumptions: Gaussian waist `w0 = 0.5 µm`, `I(r) = I0 exp(-2r^2/w0^2)`,
`I0 = 2P/(pi w0^2)`; the scattering force is `(n_medium/c) Q_pr` times the
intercepted power `P(1 - exp(-2a^2/w0^2))` (for `a << w0` this is the usual
`I0 pi a^2`); the gradient force and stiffness are the standard Rayleigh result
`F_grad = (2 pi n_medium a^3/c) K grad I` with `K = (m^2-1)/(m^2+2) = 0.0619`
(Ashkin 1986; Neuman & Block 2004).

**Indicator.** The trap is deep if `U / k_B T >> 1`; here the Rayleigh estimate
gives `2.6e4`. The Rayleigh gradient force at `x ~ 4` is only
order-of-magnitude: measured trap stiffnesses for a 1 µm bead are typically
`0.01-0.1 pN/nm` per 100 mW (Neuman & Block 2004), about 10-100× below this
point-dipole estimate, which would give `U/k_B T ~ 3e2 - 3e3` — still an
order of magnitude or more above unity, so the bead is trapped. The exact Mie
`Q_pr` is the reliable part of the calculation; the gradient estimate is
labelled as an estimate. Both forces scale linearly with power, which the test
checks.

### (c) Virion-sized particle in water (detection cross-section)

A weakly scattering particle with `m_rel = 1.1` (protein/lipid envelope,
estimate) in water at 550 nm. `C_sca = Q_sca pi a^2`:

| diameter | x | Q_sca | `C_sca` (m^2) |
|---|---|---|---|
| 20 nm | 0.1525 | 6.126e-6 | 1.925e-21 |
| 50 nm | 0.3812 | 2.307e-4 | 4.530e-19 |
| 100 nm | 0.7624 | 3.225e-3 | 2.533e-17 |
| 200 nm | 1.5247 | 2.953e-2 | 9.276e-16 |
| 500 nm | 3.8118 | 2.690e-1 | 5.282e-14 |

A 100 nm virion scatters `2.5e-17 m^2`, about `3e-3` of its geometric area
`pi a^2 = 7.85e-15 m^2`, and 64× more than a 50 nm particle (Rayleigh `d^6`
scaling, `x < 1`). Detection is therefore photon-starved: the cross-section is
of order `10^-17 m^2`, so a modest number density and a bright source are needed,
while the strong forward peak (`g > 0`) keeps most of the scattered light near
the beam.

## Coated sphere: not implemented

A coated sphere (Aden–Kerker) would give the shell/core refractive-index contrast
that homogeneous Mie cannot. It is **not** implemented here. No independent
reference (published table or a second implementation) was available in the
offline/limited-network environment, and the constraint to validate every model
against an exact or published result before claiming it works made a
self-consistent-only implementation the wrong call within the time budget. The
homogeneous sphere, both folds and all three applications are complete and
validated instead. This is a scope limit, stated rather than hidden.

## Limits and honest misses

* The model is exact for a single homogeneous, isotropic, non-magnetic sphere in
  a non-magnetic medium. It does not include coating, anisotropy, roughness,
  a complex ambient index, or orientation (no non-spherical particles).
* The Rayleigh and anomalous-diffraction folds are documented as switches with
  measured error envelopes, not certified bounds. Anomalous diffraction's error
  is oscillatory in `rho` and can peak near `rho ~ 5` before decaying as `1/x`.
* The finite-index versus perfect-conductor distinction matters: `m = 0 + 100i`
  is not the conductor limit. The published MIEV0 conductor cases require the
  dedicated `perfect_conductor` constructor.
* The cloud case uses a fixed `N = 100 cm^-3` and `v_eff = 0.1`; the gamma
  distribution is only valid for `v_eff <= 1/3` (a non-negative shape exponent).
* Water absorption at 550 nm is negligible, so the computed single-scattering
  albedo is 1 to seven digits; real clouds have trace absorbers. The optical
  constants are the crate's Daimon–Masumura fit, which is the `n` in the
  visible, not a cloud-droplet absorption model.
* The optical-trap gradient force uses the point-dipole (Rayleigh) formula at
  `x ~ 4`, where it overestimates the true stiffness by roughly an order of
  magnitude. Only `Q_pr` is exact. A ray-optics or full-field computation would
  be the unfold.
* The virion index `m_rel = 1.1` is an estimate, not a measurement; the
  cross-section scales strongly with it. The `d^6` statement is checked only for
  `x < 1`.
* Large `x` is exercised to `10^4` (about `10^4` terms); `10^5` is used only for
  the extinction-paradox check and is not in the published benchmark set.

## Reproduction

```sh
cd fluid
cargo test -p photonics --lib mie                 # 14 unit tests, ~0.5 s
cargo test -p photonics --test mie                # 6 integration tests, ~0.2 s
cargo test -p photonics --test mie -- --nocapture # fold and application tables
```

The full `cargo test -p photonics` suite passes; the Mie tests add well under a
second.

## References

* C. F. Bohren and D. R. Huffman, *Absorption and Scattering of Light by Small
  Particles*, Wiley (1983). Mie coefficients eqs. 4.52, 4.88; efficiencies
  eqs. 4.62–4.64, 4.71; angular functions eq. 4.47; Appendix A/B test sphere
  `m = 1.55`, `x = 2 pi * 0.525/0.6328`; Appendix B MIEV0 program.
* W. J. Wiscombe, "Improved Mie scattering algorithms", Appl. Opt. **19**, 1505
  (1980), doi:10.1364/AO.19.001505. Downward `D_n` recurrence, Miller's `psi_n`,
  the truncation `n_stop = x + 4.05 x^{1/3} + 2`.
* W. J. Wiscombe, "MIEV0: A portable program for Mie scattering calculations",
  NCAR Technical Note NCAR/TN-140+STR (1979). Source of the MIEV0 test cases
  (`m = 0.75`, `m = 1.33-1e-5i`, `m = 1.5-1i`, `m = 10-10i`, perfect
  conductor).
* W. J. Lentz, "Generating Bessel functions in Mie scattering calculations using
  continued fractions", Appl. Opt. **15**, 668 (1976),
  doi:10.1364/AO.15.000668. Continued fraction seeding `D_n`.
* H. C. van de Hulst, *Light Scattering by Small Particles*, Dover (1981).
  Anomalous diffraction and the extinction paradox.
* J. E. Hansen and L. D. Travis, "Light scattering in planetary atmospheres",
  Space Sci. Rev. **16**, 527 (1974), doi:10.1007/BF00168069. Gamma size
  distribution and effective radius/variance.
* A. Ashkin, J. M. Dziedzic, J. E. Bjorkholm, S. Chu, "Observation of a single-beam
  gradient force optical trap for dielectric particles", Opt. Lett. **11**, 288
  (1986), doi:10.1364/OL.11.000288.
* K. C. Neuman and S. M. Block, "Optical trapping", Rev. Sci. Instrum. **75**,
  2787 (2004), doi:10.1063/1.1785844. Typical trap stiffnesses.
* G. L. Stephens, "Radiation profiles in extended water clouds. II: Parameterization
  schemes", J. Atmos. Sci. **35**, 2123 (1978), doi:10.1175/1520-0469(1978)035<2123:RPIBWC>2.0.CO;2.
  Cloud extinction and asymmetry.
* I. H. Malitson, JOSA **55**, 1205 (1965) (fused silica); M. Daimon and
  A. Masumura, Appl. Opt. **46**, 3811 (2007) (water). Both via
  `photonics::dispersion`.
* `miepython` (S. Prahl), reference implementation and test suite used to
  cross-check the published values and the finite-index vs conductor distinction.
* mpmath 50-digit evaluation of the general Mie series, used to confirm the
  finite `m = 100i` value `Q_ext = 2.013792663` at `x = 100`.
