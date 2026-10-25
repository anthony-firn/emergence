# Photothermal coupling: light to heat to sound, and where to fold

Implemented in `fluid/crates/photothermal`. This crate is the rung that turns
the optical absorbed-power density from `fluid/crates/photonics` into heat and,
when the pulse is short enough, into an acoustic source. It also exposes the
three dimensionless numbers that decide where each coarse model stops being
valid, following the pattern of [emergence_ladder.md](emergence_ladder.md).

```
  optics: field / Beer-Lambert      --p(z)-->   conduction (this crate)
                                                        |
                       stress confinement               v
  photoacoustics:  p0 = Gamma mu_a F  --d'Alembert-->  layered acoustics
                                   (phononic transfer matrices / FDTD)
```

## 1. The models

**Conduction.** The heat equation

\[
  \rho c_p\,\partial_t T = \partial_z\!\left(k\,\partial_z T\right) + p(z,t)
\]

is solved by a conservative finite-volume scheme on a grid aligned with the
layer interfaces. Cell temperature is a cell average and the face flux is
`G (T_L - T_R)` with the series conductance
`G = 1/(d_L/2 k_L + d_R/2 k_R)`. The interface flux is continuous by
construction and the scheme reduces to `k/dz` inside a layer. Time stepping is
backward Euler (unconditionally stable for diffusion), with a Thomas solve.
Boundaries are adiabatic, fixed-temperature, or fixed surface flux. The absorbed
power can be sampled from `photonics::multilayer::Solution::power_density`
([`source_from_optics`]) or taken analytically as Beer–Lambert
`p = mu_a I_0 e^{-mu_a z}`.

**Photoacoustics.** Under stress confinement an absorbed fluence `F` launches
the initial pressure

\[
  p_0(z) = \Gamma\,\mu_a F\,e^{-\mu_a z},
  \qquad \Gamma = \frac{\beta c^2}{c_p}.
\]

The 1D wave equation with `p(z,0)=p_0`, `p_t(z,0)=0` has the exact d'Alembert
solution `p = ½[P_0(z-ct) + P_0(z+ct)]`, where `P_0` is the **odd** extension of
`p_0` for a free (pressure-release) surface and the **even** extension for a
rigid surface. The free surface therefore produces a bipolar signal and the
rigid surface a monopolar one.

For a layered stack the initial condition is propagated by a staggered
finite-difference solver (pressure at nodes, velocity at faces, harmonic
interface modulus and exact half-step initialisation). Independently, a pulse
can be launched in the frequency domain: FFT the incident pressure, multiply by
the pressure transmission coefficient of a `phononic::layered::Stack`, inverse
FFT (with the exact DC limit `2 Z_exit/(Z_in+Z_exit)`).

**Indicators.** The hand-offs are cheap dimensionless numbers: the thermal
confinement ratio `tau_p alpha / delta^2`, the stress confinement ratio
`tau_p c / delta`, the coherence indicator `n d / L_c`, and the Rayleigh number
`Ra = g beta Delta T L^3/(nu alpha)`. They are described in section 4.

## 2. Benchmarks (measured)

Run `cd fluid && cargo test -p photothermal -- --nocapture`. Every value below
is printed by the tests; "error" is the deviation from the exact or published
result. No tolerance is loosened to make a case pass.

| Benchmark | Exact / published | Computed | Error |
|---|---|---|---|
| Complementary error function `erfc(0.5,1,2,3)` | 0.4795001221869535, 0.15729920705028513, 4.677734981047266e-3, 2.209049699858544e-5 | same | `<=6.8e-21` |
| Semi-infinite surface flux, `T(z=dz/2)` at 1 s | erfc solution | — | `4.4e-4` K (rel. L2 `6.0e-4`) |
| Steady Beer–Lambert slab, 400 cells | analytic | — | `1.53e-5` K, order `1.998` |
| Transient Beer–Lambert slab (`dt=0.005 s`, 5 s) | Fourier series (300 modes) | — | rel. L2 `1.66e-5`, time order `1.10` |
| Adiabatic layered energy, absorbed = stored | equality | ratio `1.0` | `6.0e-9` |
| Energy ledger with outflow, deposited = stored + conducted | equality | residual | `1.7e-11` |
| Interface flux continuity, front flux `q` | `q` on every face | — | `1.2e-8` |
| `photonics` gold film source `∫p dz` vs `I0 A` | equal | — | `5.5e-6` |
| ARC window (MgF2 / 2 µm glass) integrated source | `I0 A` | — | `9.7e-5` (A = `4.5e-3`) |
| Water 10.6 µm, `mu_a = 4 pi kappa/lambda` | 5.9275e4 m⁻¹ (Hale–Querry `kappa=0.05`) | 5.9275e4 m⁻¹ | round-off; `delta=16.87` µm |
| Water 10.6 µm Beer–Lambert absorbed power | exact | — | `2.9e-16` |
| Grüneisen water from `beta c^2/c_p` | 0.11 (quoted) | 0.108713 | `1.3e-3` |
| Free vs rigid surface shape | bipolar / monopolar | range `[-10.8,10.7]` Pa / min `1.1e-17` Pa | sign structure exact |
| FDTD Gaussian vs d'Alembert (rigid) | exact | — | rel. L2 `1.37e-3` |
| FDTD exponential vs d'Alembert (rigid) | exact | — | rel. L2 `9.94e-3` (first order, moving kink) |
| FDTD smooth profile vs d'Alembert (free) | exact | — | rel. L2 `2.47e-3` |
| Spectral phononic layer, pure delay `d/c` | exact shift | — | rel. L2 `3.0e-14` |
| Spectral single interface, `t = 2Z2/(Z1+Z2)` | 1.938700002 | — | `1.2e-13` |

The erfc, energy, flux-continuity, pure-delay and Fresnel rows are exact
identities, so their errors are round-off. The d'Alembert and Fourier-series
rows are numerical-scheme errors.

### Spatial convergence (steady Beer–Lambert slab)

| cells | max error [K] | observed order |
|---|---|---|
| 50 | 9.703e-4 | — |
| 100 | 2.436e-4 | 1.994 |
| 200 | 6.104e-5 | 1.997 |
| 400 | 1.528e-5 | 1.998 |

The conservative interface conductance is second-order.

### Time convergence (transient slab vs Fourier series)

| `dt` [s] | relative L2 | observed order |
|---|---|---|
| 0.02 | 7.439e-5 | — |
| 0.01 | 3.556e-5 | 1.065 |
| 0.005 | 1.664e-5 | 1.096 |

Backward Euler is first order, as expected.

### Photoacoustic solver convergence

The staggered scheme is second order for smooth fields: the Gaussian rigid
case has relative L2 `1.37e-3` at `sigma/dz = 15`, and the smooth free-surface
profile falls `7.6e-3 -> 2.5e-3 -> 7.9e-4 -> 2.5e-4` when the grid is refined
`10 -> 20 -> 40 -> 80` cells per absorption depth (order `~1.7`). The
*exponential* profile even-extended at a rigid wall has a kink in its third
derivative, so the exact field carries a moving kink and the error falls only
first order (`8.1e-2 -> 4.0e-2 -> 2.0e-2 -> 9.9e-3`); this is reported rather
than hidden.

## 3. The example stacks

* **Gold film (50 nm) on glass under water**, 633 nm: the photonics
  Drude–Lorentz field deposits `I0 A` in the film; the thermal source integrates
  to it to `5.5e-6` and reaches a 10.2 K steady rise at `I0 = 1e5 W/m²` with a
  fixed glass back.
* **Anti-reflection-coated window** (MgF2 quarter wave at 550 nm on 2 µm of
  faintly absorbing glass): `R = 0.094`, `A = 4.5e-3`; the coherent transfer
  matrix is used only because the pane is a few wavelengths thick (see the
  coherence indicator below).
* **Strongly absorbing water at 10.6 µm** (`kappa = 0.05`, Hale & Querry 1973):
  `mu_a = 5.9275e4 m⁻¹`, absorption depth `16.87 µm`, Beer–Lambert heating stores
  exactly the absorbed energy, and the stress-confined photoacoustic peak for
  `F = 0.01 J/m²` is `p_0 = 64.5 Pa`.

## 4. Fold / unfold map

| coarse model (fold) | cheap indicator | unfold when | where |
|---|---|---|---|
| Beer–Lambert source from `I0 A` | layer thickness vs `delta = 1/mu_a` | `delta <` layer thickness | `photonics::multilayer::power_density` |
| Coherent field / transfer matrix | `n d / L_c` | `n d / L_c > 1`: use incoherent intensity sum | `photonics::incoherent`, `photothermal::indicators::optical_regime` |
| Instantaneous heat deposition (no heat flow during the pulse) | thermal confinement `tau_p alpha / delta^2` | `> 1`: solve conduction during the pulse | `thermal::ThermalModel` |
| Instantaneous pressure deposition (`p0 = Gamma mu_a F`) | stress confinement `tau_p c / delta` | `> 1`: solve the coupled acoustic problem during the pulse | `photoacoustic` |
| Conduction-only transport | Rayleigh number `Ra` | `Ra > 1708` (surface-like) or `Ra_I > 2772` (volumetric): unfold to the flow solver | `indicators` |
| FDTD layered acoustics | grid points per absorption depth | too few points: first-order kink smearing | `photoacoustic::AcousticModel` |

### Thermal vs stress confinement

For water at 20 °C (`alpha = 1.47e-7 m²/s`, `c = 1482 m/s`) and an absorption
depth `delta = 10 µm`, the thermal relaxation time is `delta²/alpha = 0.68 ms`
and the acoustic transit time is `delta/c = 6.7 ns`. A 1 ns pulse is therefore
thermally confined (ratio `1.5e-6`) but only marginally stress confined (ratio
`0.148`); 1 ns is short enough that the `p0 = Gamma mu_a F` initial condition is
appropriate, while a 1 µs pulse would not be stress confined.

### Optics: coherent vs incoherent

The coherence indicator is `n d / L_c`. The tests confirm the photonics
crate's own finding: a 1 mm glass pane (optical thickness `~2700` wavelengths)
is far beyond the coherent transfer matrix's numerical range, and the coherent
solver returns an absorptance of `0.209` that is round-off noise, not physics.
The correct unfold is the incoherent intensity sum (`photonics::incoherent`) or
Beer–Lambert. The photothermal crate therefore samples the coherent field only
for optically thin films (the gold film, the few-wavelength ARC pane).

### Conduction vs convection

For a 1 cm water layer with `nu = 1e-6 m²/s`:

* `Delta T = 0.01 K` gives `Ra = 145 < 1708`: conduction only;
* `Delta T = 1 K` gives `Ra = 14505 > 1708`: unfold to the flow solver.

For volumetric absorbed heating (absorption depth not small compared with the
layer) the internal Rayleigh number `Ra_I = g beta q L^5/(k nu alpha)` is used
with the internal critical value `2772`; `q = 1e4 W/m³` over 1 cm gives
`Ra_I = 2.4e4 > 2772`. `absorbed_heating_rayleigh` switches on `mu_a L >= 3`.

## 5. Limits and honest misses

* The thermal solver is 1D and purely conductive: no radiation, no
  temperature-dependent properties, no phase change, no contact resistance and
  no convection. It does not feed temperature back into the optical constants.
* The source is treated as a time-independent volumetric term; a moving or
  pulsed source needs a step-by-step source update (the `step` API allows this).
* The photoacoustic initial condition assumes stress confinement. Outside it,
  the coupled thermo-acoustic problem must be solved during the pulse; the
  indicator says when.
* The acoustic solver is lossless and normal incidence only; the fluid–plate,
  oblique-incidence and anisotropic cases belong to `phononic`.
* The spectral launch uses only the pressure transmission coefficient, so it
  gives the transmitted waveform at the back face and cannot represent internal
  reflections returning to the source. The time-domain solver is the tool for
  that.
* The exponential initial pressure at a free surface has a pressure jump; the
  exact d'Alembert solution is discontinuous and the FDTD error converges only
  first order. We report this instead of tuning the grid.
* The Grüneisen value is a single room-temperature number; `Gamma` is
  temperature dependent, and the beta used (`2.07e-4 K⁻¹`) is the volumetric
  coefficient of water at 20 °C.
* The internal-heating critical Rayleigh number `2772` is inherited from the
  literature (Sparrow et al. 1964); the test only exercises clear-cut cases
  (`Ra_I` well above and below it).

## References

* F. P. Incropera and D. P. DeWitt, *Fundamentals of Heat and Mass Transfer*
  (thermal properties of water, fused silica, Pyrex, gold; cp and k tables).
* *CRC Handbook of Chemistry and Physics* (water density and volumetric thermal
  expansion, gold density and heat capacity).
* G. M. Hale and M. R. Querry, "Optical constants of water in the 200-nm to
  200-µm wavelength region", Appl. Opt. **12**, 555 (1973), doi:10.1364/AO.12.000555
  (water `n`, `kappa` at 10.6 µm).
* M. Daimon and A. Masumura, Appl. Opt. **46**, 3811 (2007);
  I. H. Malitson, JOSA **55**, 1205 (1965) (water and silica optics, via
  `photonics`).
* A. D. Rakić et al., Appl. Opt. **37**, 5271 (1998) (gold dispersion, via
  `photonics`).
* P. G. Drazin and W. H. Reid, *Hydrodynamic Stability* (1981); S. Chandrasekhar,
  *Hydrodynamic and Hydromagnetic Stability* (1961) (`Ra_c = 1707.762` for
  rigid, isothermal Rayleigh–Bénard convection).
* E. M. Sparrow, R. J. Goldstein and V. K. Jonsson, "Thermal instability in a
  horizontal fluid layer: effect of boundary conditions and non-linear
  temperature profile", J. Fluid Mech. **18**, 513 (1964) (internal-heating
  critical Rayleigh number `2772`).
* L. M. Brekhovskikh and O. A. Godin, *Acoustics of Layered Media* (layered
  transfer matrices, via `phononic`).
