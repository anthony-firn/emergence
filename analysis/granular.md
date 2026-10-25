# Granular matter and suspensions: particles folding to a continuum rheology

Implemented in `fluid/crates/granular` (new crate, no dependencies on the
other crates). This is the particulate rung of
[emergence_ladder.md](emergence_ladder.md). Two fine models live here and
two closures fold them to the continuum:

| coarse model | validity indicator | unfold / closure | where |
|---|---|---|---|
| Newtonian fluid with a viscosity `eta(phi)` | volume fraction `phi` vs `phi_m`; particle Reynolds `Re_p`; Péclet `Pe` | exact single-sphere Stokes stresslet (dilute) → Batchelor–Green → Krieger–Dougherty | `stokes.rs`, `suspension.rs` |
| fluid-immersed dense suspension | viscous number `J = eta_f gdot / P` | `mu(J)`, `phi(J)` (Boyer et al. 2011) | `rheology.rs` |
| dry granular flow | inertial number `I = gdot d / sqrt(P/rho_p)` | `mu(I)` (Jop et al. 2006) | `dem.rs`, `rheology.rs` |

Everything runs in seconds per test. The DEM uses reduced units fixed by
the particle diameter `d = 1` and density `rho_p = 1`; the shear cell also
works in `d = 1` (radius 0.5). Heaps and the larger shear sweep live in
`examples/large_runs.rs` (run with `--release`).

## 1. Suspension rung

### 1.1 Stokes flow and the drag

`stokes::translating_sphere_velocity` is the exact solution for a rigid
sphere of radius `a` translating with velocity `V` in a fluid at rest at
infinity (Landau & Lifshitz, *Fluid Mechanics* 2nd ed. §20; Happel &
Brenner 1965 §3-3):

```text
u = V A(r) + (V.n) n B(r),   n = x/r,
A = (3/4)(a/r) + (1/4)(a/r)^3,   B = (3/4)(a/r)(1 - (a/r)^2),
p = (3/2) mu a (V.x)/r^3.
```

The gradient is analytic and is checked against central finite differences
to `4e-11`. The surface traction is uniform, `t = -(3 mu / 2a) V`, and its
surface integral is the Stokes drag `F = -6 pi mu a V` (Stokes 1851).

| benchmark | measured | exact | error |
|---|---|---|---|
| drag `\|F\|` from a 24×48 Gauss–Legendre/trapezoid surface quadrature of the stress tensor | `1.785967323083e1` | `6 pi mu a \|V\| = 1.785967323083e1` | **3.4e-15** relative |

The force direction is antiparallel to `V` to `1e-9`.

### 1.2 Einstein's viscosity from the stresslet

For a rigid, force- and torque-free sphere in an imposed symmetric
traceless rate of strain `E`, the stresslet is exactly
`S = (20/3) pi mu a^3 E` (Batchelor 1970; Kim & Karrila 2005 §3.3). The
surface traction that realises it is `t = 5 mu E.n`; integrating its
symmetric first moment over the sphere numerically returns `S`.

| benchmark | measured | exact | error |
|---|---|---|---|
| stresslet from a 24×48 surface quadrature of `t = 5 mu E.n` | matches | `(20/3) pi mu a^3 E` | **3.2e-15** relative |
| dilute bulk stress `2 mu E + n S` for `phi = 0.1` | `eta/eta0 = 1.25` | `1 + 5/2 phi = 1.25` | **0** (machine) |

Summing `N` exact single-sphere stresslets, with
`(N/V)(20/3) pi a^3 = (5/2) phi`, reproduces `eta = eta0 (1 + 5/2 phi)`
exactly. So the Einstein coefficient `[eta] = 5/2` is obtained from the
exact stresslet rather than fitted. **Which method:** the task's third
option — the exact single-sphere stresslet summed in the dilute limit — is
the one used; the numerics that accompany it are genuine surface
quadratures of the exact stress field (drag) and of the strain traction
(stresslet).

### 1.3 Batchelor–Green and Krieger–Dougherty

Batchelor & Green (1972) computed the `O(phi^2)` coefficient for a random
suspension of spheres **without Brownian motion** (the high-frequency
viscosity): `k2 = 5.2`. Batchelor (1977) added the Brownian stresses and
found `k2 = 6.17`; the pair-Smoluchowski result of Cichocki & Felderhof
(1988) is `k2 = 5.91` (all three as compiled by Mendoza &
Santamaría-Holek 2009, arXiv:0810.5111, Eqs. 1–2). This crate provides all
three. **The one compared against the `O(phi^2)` series is the
non-Brownian `5.2`; the Brownian `6.17` is reported but not used as the
primary closure.**

| quantity | value |
|---|---|
| Batchelor & Green (1972) `k2`, non-Brownian | **5.2** (numeric series exactly 5.200000) |
| Batchelor (1977) `k2`, Brownian | 6.17 |
| Krieger–Dougherty (1959) implied `k2` at `phi_m = 0.63`, `[eta] = 5/2` | **5.109** |
| KD divergence exponent `[eta] phi_m` at `phi_m = 0.63` | 1.575 |

The KD closure is `eta/eta0 = (1 - phi/phi_m)^(-[eta] phi_m)`. Hard-sphere
maximum packing runs from random close packing `phi_m ~ 0.63` (low-shear) to
fcc `phi_m ~ 0.7404` (high-shear); Mendoza & Santamaría-Holek (2009) report
`phi_m = 0.63` fits the low-shear data of **de Kruif et al. (1985)**
(*J. Chem. Phys.* **83**, 4717) and **Krieger (1972)** (*Adv. Colloid
Interface Sci.* **3**, 111). The implied `k2 = 5.11` is within 3% of
Batchelor–Green's 5.2, which is a genuine consistency check of the closure.

### 1.4 Dense reference points, and where KD misses

Comparison against published low-shear hard-sphere relative viscosities
(Phillies 2004 reanalysis, arXiv:cond-mat/0403305, Figure data):

| `phi` | published `eta/eta0` | Krieger–Dougherty (`phi_m=0.63`) | relative error |
|---|---|---|---|
| 0.41 (dynamic crossover) | 11 | 5.24 | **-52%** |
| 0.494 (hard-sphere melting) | 49 | 11.19 | **-77%** |

The KD closure is a non-Brownian / high-frequency form. It **misses** the
Brownian low-shear viscosity near freezing, where the zero-shear viscosity
diverges much faster than KD's `(1 - phi/phi_m)^(-1.575)`; the observed
divergence near jamming is closer to `(phi_c - phi)^(-2)` (Boyer et al.
2011). This is reported rather than tuned away: KD is left as the
high-frequency/dense closure it is, and the Brownian regime is handed to
the Brownian form (Batchelor 1977 / Cichocki–Felderhof 1988).

### 1.5 Indicators: when the Newtonian suspension fold fails

`particle_reynolds = rho_f gdot a^2 / eta_f`, `Pe = 6 pi eta_f a^3 gdot /
(k_B T)`, and `phi/phi_m` classify a suspension
(`suspension::suspension_regime`):

| regime | condition |
|---|---|
| dilute Newtonian (Einstein/Stokes) | `phi < 0.1`, `Re_p << 1`, `Pe >> 1` |
| dense Newtonian (BG / KD) | `phi` up to `phi_m`, `Re_p << 1` |
| Brownian | `Pe <~ 1` |
| inertial shear thickening | `Re_p >~ 1` — the Stokes closure fails |
| jamming | `phi/phi_m >~ 0.95` |

## 2. Granular rung

### 2.1 Contact law

`contact::hertz_mindlin` implements the standard soft-sphere DEM contact:
Hertz normal elasticity `F_n = (4/3) E* sqrt(R*) delta^(3/2)`, Tsuji,
Tanaka & Ishida (1992) viscous damping with
`beta = -ln(e)/sqrt(ln^2 e + pi^2)`, a Cundall–Strack tangential spring
`k_t = S_t = 8 G* sqrt(R*) sqrt(delta)` with the Coulomb limit
`|F_t| <= mu F_n`, and an optional rolling-resistance torque
(Zhou et al. 1999; Ai et al. 2011). `1/E* = (1-nu_1^2)/E_1 +
(1-nu_2^2)/E_2`, `1/R* = 1/R_1 + 1/R_2`.

### 2.2 Hertz contact time

For a head-on elastic Hertz collision the contact duration is exactly

```text
t_c = 2.868 (m*^2 / (R* E*^2 v))^(1/5),
```

where `2.868 = 2 I (15/16)^(2/5)` and
`I = ∫_0^1 (1-u^(5/2))^(-1/2) du = 1.4714` (Johnson, *Contact Mechanics*,
1985 §11.7). The integration in `dem::two_sphere_collision` uses
`dt = t_c/400`:

| impact speed `v` | simulated `t_c` | exact `t_c` | relative error |
|---|---|---|---|
| 0.01 | `4.458689e-1` | `4.458689e-1` | **< 1e-3** (printed 0.000) |
| 0.05 | `3.231567e-1` | `3.231567e-1` | **< 1e-3** |
| 0.2 | `2.449070e-1` | `2.449070e-1` | **< 1e-3** |

### 2.3 Coefficient of restitution

With the Tsuji damping set from a target `e`, the measured restitution
follows the target closely (head-on, no friction, impact `v = 0.1`):

| target `e` | measured | error |
|---|---|---|
| 0.3 | 0.2997 | 0.0003 |
| 0.5 | 0.4999 | 0.0001 |
| 0.7 | 0.7001 | 0.0001 |
| 0.9 | 0.9000 | 0.0000 |

### 2.4 Angle of repose (a limited, non-robust result)

A quasi-2D heap (periodic in `z`, floor and two side walls) of 176 spheres
is collapsed from a dense block (`examples/large_runs.rs`, sliding friction
0.5, 40 000 steps). For spheres the heap angle is controlled mainly by
**rolling** friction, since sliding friction alone lets the spheres roll
apart (Zhou et al. 1999; Ai et al. 2011; for the DEM reference values see
the review of Ai et al. 2011, *Powder Technol.* **206**, 269):

| rolling friction `mu_r` | run A angle | run B angle |
|---|---|---|
| 0.0 | 2.0° | 6.5° |
| 0.2 | 2.8° | 2.2° |
| 0.4 | 5.8° | 5.9° |

This is reported as a **miss**: with only 176 grains the heap does not
reach a reproducible steady profile, and the measured angles scatter over
2–7° with no clean monotonic trend (the two runs use different seeds).
The example prints the same measurement, so the scatter is reproducible,
not a transcription error. A quantitative angle of repose needs several
thousand grains, slow pouring and non-spherical grains or strong rolling
resistance; the small DEM here is not that benchmark. The rolling-friction
trend is visible only in the larger of the two runs (2.0° → 5.8°).

### 2.5 Simple shear and the inertial number

`dem::ShearCell` imposes a mean shear `v_x = gdot y` with the SLLOD
equations for the peculiar velocity,

```text
xdot = u + gdot y e_x,    udot = F/m - gdot u_y e_x,
```

with Lees–Edwards periodic images, a Berendsen constant-pressure barostat
(target `P0 = 0.05`) and the virial stress. It measures
`mu = |sigma_xy|/P` and `I = gdot d / sqrt(P/rho_p)`.

Test configuration: 108 spheres, radius 0.5, density 1, `E = 1000`,
`nu = 0.3`, `mu_particle = 0.4`, target restitution 0.5, `dt = 0.005`,
5000 equilibration + 5000 production steps at each rate.

| `gdot` | `I` | `mu` | `phi` |
|---|---|---|---|
| 0.005 | 0.0224 | 0.3298 | 0.5914 |
| 0.010 | 0.0447 | 0.3887 | 0.5896 |
| 0.020 | 0.0894 | 0.3956 | 0.5921 |
| 0.040 | 0.1788 | 0.4344 | 0.5674 |
| 0.080 | 0.3578 | 0.5359 | 0.5454 |
| 0.160 | 0.7168 | 0.5903 | 0.5232 |

The dilatancy (`phi` falls with `gdot`) and the rise of `mu` are the
expected behaviour. The `sign` of `sigma_xy` is a convention; `mu` uses the
magnitude.

## 3. Fold to the continuum: `mu(I)`

### 3.1 Jop, Forterre & Pouliquen (2006)

`mu(I) = mu_s + (mu_2 - mu_s)/(I0/I + 1)`, with glass-bead constants
`mu_s = tan 20.9° = 0.3819`, `mu_2 = tan 32.8° = 0.6445`, `I0 = 0.279`
(Nature **441**, 727, 2006).

Fitting the DEM points above (`rheology::Jop2006::fit`, which scans `I0`
and solves the two-parameter linear problem) gives

| parameter | DEM fit | Jop glass beads (2006) |
|---|---|---|
| `mu_s` | 0.322 | 0.382 |
| `mu_2` | 0.751 | 0.644 |
| `I0` | 0.412 | 0.279 |
| RMS residual | 0.015 | — |

The DEM reproduces the **shape** of the glass-bead `mu(I)` law: the
quasi-static friction (0.32) sits near the interparticle friction
`mu_particle = 0.4`, and the large-`I` friction is within 0.11 of Jop's
`mu_2`. The cross-over `I0` is 48% larger. These differences are expected:
DEM friction, restitution and contact stiffness (here `E = 1000`, which
sets a finite overlap of order 1%) shift the fitted numbers. The fit is
not tuned to the Jop constants.

### 3.2 Regimes and the fluid-immersed fold

`rheology::granular_regime(I)`:

| `I` | regime |
|---|---|
| `< 1e-3` | quasi-static (Coulomb friction; no rate dependence) |
| `1e-3 – 0.1` | dense flow, `mu(I)` |
| `0.1 – 0.3` | transition |
| `> 0.3` | collisional (kinetic theory) |

When the grains are immersed in a fluid, inertia is replaced by viscous
stress and the relevant number is the viscous number
`J = eta_f gdot / P` (Boyer, Guazzelli & Pouliquen 2011, *Phys. Rev. Lett.*
**107**, 188301). `rheology::Boyer2011` implements their closures

```text
mu(J)  = mu1 + (mu2 - mu1)/(1 + I0/J) + J + (5/2) phi_m sqrt(J),
phi(J) = phi_m / (1 + sqrt(J)),
```

with frictional particles `phi_m = 0.575`, `mu1 = 0.32`, `mu2 = 0.7`,
`I0 = 0.005`, and frictionless particles `phi_m = 0.644`, `mu1 = 0.1` (as
quoted for Boyer et al. 2011 in Guazzelli & Pouliquen 2018,
arXiv:1801.03805). This is the bridge between the suspension rung and the
granular rung: `J -> 0` recovers the jamming packing `phi_m`, while
`J -> inf` recovers the high-`J` friction `mu2`.

## 4. Honest misses and limitations

- **Krieger–Dougherty at high `phi` and low shear.** KD with `phi_m = 0.63`
  underestimates the Brownian zero-shear hard-sphere viscosity by 52% at
  `phi = 0.41` and 77% at `phi = 0.494`. KD's divergence exponent 1.575 is
  weaker than the observed ~2. This is a genuine limitation of using a
  non-Brownian closure in the Brownian regime, not a numerical error.
- **Brownian `O(phi^2)`.** The Brownian coefficient is quoted as
  Batchelor's 6.17, but the exact pair-Smoluchowski value is 5.91
  (Cichocki & Felderhof 1988); the crate provides both and flags the
  difference.
- **Angle of repose.** Not validated. A 176-grain quasi-2D heap gives
  scattered angles (2–7°) with no clean trend across seeds; see §2.4. A
  quantitative measurement needs a much larger heap and slow pouring.
- **`mu(I)` cross-over.** DEM `I0 = 0.41` vs the glass-bead 0.279. Contact
  stiffness and interparticle friction set this; it was not tuned.
- **Shear-cell geometry.** A single small cell (108–180 particles) at one
  pressure cannot resolve the pressure dependence of `mu_s` that the full
  `mu(I)` rheology has; the fit is a single pressure slice.

## 5. Sources

- G. G. Stokes (1851), drag `6 pi mu a U`; Landau & Lifshitz, *Fluid
  Mechanics* 2nd ed. §20; Happel & Brenner (1965) §3-3.
- A. Einstein (1906) *Ann. Phys.* **19**, 289 — `[eta] = 5/2`.
- G. K. Batchelor (1970) *J. Fluid Mech.* **41**, 545 — stresslet.
- G. K. Batchelor & J. T. Green (1972) *J. Fluid Mech.* **56**, 401 —
  non-Brownian `k2 = 5.2`.
- G. K. Batchelor (1977) *J. Fluid Mech.* **83**, 97 — Brownian `k2 = 6.17`.
- B. Cichocki & B. U. Felderhof (1988) *J. Chem. Phys.* **89**, 3705 —
  exact `O(phi^2)`, `k2 = 5.91`.
- I. M. Krieger & T. J. Dougherty (1959) *Trans. Soc. Rheol.* **3**, 137.
- C. G. de Kruif et al. (1985) *J. Chem. Phys.* **83**, 4717; I. M. Krieger
  (1972) *Adv. Colloid Interface Sci.* **3**, 111.
- C. I. Mendoza & I. Santamaría-Holek (2009) arXiv:0810.5111 — compilation
  of the coefficients and KD fits.
- G. D. J. Phillies (2004) arXiv:cond-mat/0403305 — hard-sphere
  relative-viscosity reference points.
- P. A. Cundall & O. D. L. Strack (1979) *Géotechnique* **29**, 47.
- Y. Tsuji, T. Tanaka & T. Ishida (1992) *Powder Technol.* **71**, 239.
- R. D. Mindlin (1949) *J. Appl. Mech.* **16**, 259.
- K. L. Johnson, *Contact Mechanics* (1985) §11.7 — Hertz contact time.
- P. Jop, Y. Forterre & O. Pouliquen (2006) *Nature* **441**, 727 — `mu(I)`.
- F. Boyer, É. Guazzelli & O. Pouliquen (2011) *Phys. Rev. Lett.* **107**,
  188301 — `mu(J)`, `phi(J)`.
- É. Guazzelli & O. Pouliquen (2018) arXiv:1801.03805 — review and
  compilation of the Boyer et al. constants.
- J. Ai, J.-F. Chen, J. M. Rotter & J. Y. Ooi (2011) *Powder Technol.*
  **206**, 269 — rolling friction and the angle of repose.
