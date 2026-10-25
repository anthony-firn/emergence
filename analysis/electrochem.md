# Electrochemistry and batteries: the interface, the cell, and where each model folds

Implemented in `fluid/crates/electrochem` (new crate; added to the workspace
`members`). This crate walks the electrochemistry rung of
[emergence_ladder.md](emergence_ladder.md): a cheap interface/cell model, a
cheap validity indicator, and an unfold where the indicator fails.

```
  equivalent circuit (Thevenin RC)         BMS / system level
        ^  indicator: pulse length / tau_RC
  single-particle model (SPM)              cell control level
        ^  indicator: t_electrolyte / t_discharge
  electrolyte-enhanced SPM (SPMe)          finer cell model
        ^  (the true reference is the pseudo-2D DFN; implemented as `dfn`)
  interface kinetics + 1D transport        electrochemical rung
```

Module map: `interface` (Nernst, Butler–Volmer, Tafel/linear),
`diffusion` (Crank–Nicolson, Cottrell), `cv` (reversible cyclic voltammetry,
Randles–Ševčík), `particle` (spherical solid diffusion, exact constant-flux
series), `spm` (single-particle cell), `spme` (electrolyte-enhanced SPM),
`dfn` (the pseudo-2D Doyle–Fuller–Newman reference), `ecm` (Thevenin RC), `thermal` (Bernardi heat and the `heat` interface).

---

## 1. Interface rung: Nernst, Butler–Volmer and its limits

For `O + n e⁻ ⇌ R` (Bard & Faulkner, *Electrochemical Methods*, 2nd ed. 2001):

```
Nernst          E = E° + (RT/nF) ln(a_O/a_R)
Butler–Volmer   i = i₀ [exp(α_a nFη/RT) − exp(−α_c nFη/RT)]
Tafel   (η≫RT/nF)   i ≈ i₀ exp(α_a nFη/RT)
linear  (η≪RT/nF)   i ≈ i₀ nFη/RT
```

The Tafel and linear branches are the two folds of Butler–Volmer, and
`butler_volmer_overpotential` inverts the full law by Newton (warm-started from
`2 asinh(i/2i₀)/f`) so the cell model can use the exact kinetics.

| benchmark | result | error |
|---|---|---|
| Tafel slope `2.303RT/(αnF)`, α=0.5, n=1, 298.15 K | 118.319 mV/dec | exact (numeric vs closed form, test) |
| Nernst decade `RT/F ln10` at 298.15 K | 59.159 mV | exact (test) |
| Butler–Volmer → linear as η→0 | `i/i_lin − 1 = O(η²)` | < 1e-9 in test |
| Butler–Volmer → Tafel for large η | relative residual | < 1e-5 at η = 0.3 V |
| BV inverse round-trip | `i → η → i` | < 1e-9 relative |

Mass transport is folded through a Nernst diffusion layer:
`i_L = nF D c_bulk/δ`, surface depletion `c_s = c_bulk(1 − i/i_L)`, and a
Koutecký–Levich series `i = i_k i_L/(i_k+i_L)` between kinetics and the plateau.

---

## 2. Transport benchmarks

### Cottrell transient (`diffusion`)

A potential step to a transport-limited value clamps the surface concentration,
and the current is the exact
`i(t) = nFAc√(D/(πt))` (Bard & Faulkner eq. 5.2.11). The finite-volume
Crank–Nicolson solver (`Diffusion1D`, Thomas solve) is validated on the exact
erf solution of the same problem:

| t (D=1, dimensionless) | numerical surface flux | exact `√(D/πt)` | rel. error |
|---|---|---|---|
| 0.5 | 0.797959 | 0.797885 | 0.0093 % |
| 1.0 | 0.564216 | 0.564190 | 0.0046 % |

The scheme's **order** is checked on a smooth Neumann eigenmode
(`cos(πx/L) e^{−Dπ²t/L²}`), where it converges at **2nd order in both Δt and
Δx** (test `eigenmode_second_order_in_time/_space`; the differences isolate the
O(Δt²), O(Δx²) terms). The potential-step problem itself has a √t singularity
at t=0, which is why the Cottrell error above is ~1e-4 relative rather than at
round-off.

### Reversible cyclic voltammetry (`cv`)

For equal diffusion coefficients the total concentration is uniform, so the
Nernst relation at the electrode is a time-dependent Dirichlet boundary
`a(0,τ) = 1/(1+e^{−σ})`, `σ = nF(E−E°)/RT`. The dimensionless current peaks at
the **Randles–Ševčík constant** (Nicholson & Shain, *Anal. Chem.* **36**, 706
(1964); Bard & Faulkner eq. 6.2.12):

```
i_p = 0.4463 n F A c √(n F v D/(R T)).
```

Simulated peak: **0.446292** vs 0.4463 → **−0.0018 %**. A second-order
one-sided surface gradient is required for this accuracy; the naive
`(a₁−a₀)/Δx` is only first order and gave +4 % at Δx=0.01. The peak scales as
`√v` and as `n^{3/2}` (tests).

### Spherical solid diffusion (`particle`)

A finite-volume sphere with an imposed outward surface flux `J₀` conserves moles
exactly. The exact constant-flux surface-concentration series is derived by
eigenfunction expansion (Crank, *The Mathematics of Diffusion*, 2nd ed. 1975,
ch. 6):

```
c_s(t) = c₀ − 3J₀t/R − (2RJ₀/D) Σ_{n≥1} (1−e^{−D μ_n² t/R²})/μ_n²,   tan μ_n = μ_n.
```

| quantity | result | error |
|---|---|---|
| volume average `c̄ = c₀ − 3J₀t/R` | exact | < 1e-11 (test `volume_average_is_exact_linear`) |
| surface vs series (n=50, t=0.3R²/D) | 0.989101 vs 0.989112 | 0.10 % |
| surface vs series (n=100) | 0.989052 | 0.55 % |
| surface vs series (n=200) | 0.989027 | 0.78 % |

The surface error is set by the thin outer shell resolving a flux boundary; the
volume average is exact to round-off. The zero-flux state is a fixed point and
the total moles track `A_R J₀ t` exactly (tests).

---

## 3. Cell rung: single-particle model (SPM)

Each porous electrode is folded to one spherical particle with a **uniform**
reaction current `i_loc = i_app/(a_s L)`, `a_s = 3ε_s/R`. Solid diffusion is
resolved (radially), the surface reaction is Butler–Volmer, and

```
V = U_p(x_p^s) − U_n(x_n^s) + η_p − η_n − i_app R_Ω.
```

Parameters: C.-H. Chen et al., *J. Electrochem. Soc.* **167**, 080534 (2020),
"LG M50" (the PyBaMM `Chen2020` set). Geometry, diffusivities, maximum
concentrations and initial concentrations are from that set:

| parameter | value | source |
|---|---|---|
| L_n, L_p | 85.2, 75.6 µm | Chen 2020 |
| R_n, R_p | 5.86, 5.22 µm | Chen 2020 |
| c_max,n, c_max,p | 33133, 63104 mol m⁻³ | Chen 2020 |
| c_n0, c_p0 | 29866, 17038 mol m⁻³ | Chen 2020 |
| D_s,n, D_s,p | 3.3e-14, 4.0e-15 m² s⁻¹ | Chen 2020 |
| ε_s,n, ε_s,p | 0.75, 0.665 | Chen 2020 (approx.) |
| k_n, k_p | 2.334e-11, 6.48e-11 m^2.5 mol^-0.5 s^-1 | **ESTIMATE** |
| R_Ω | 2.0e-3 Ω m² | **ESTIMATE** |

**Honesty note.** The negative-electrode OCP is the published Chen 2020 graphite
fit. The **positive-electrode OCP is a surrogate sigmoid**
(`ocp_nmc = 3.85 + 0.40 tanh(−4(x−0.45))`); the published NMC fit was not
available offline. The model therefore reproduces the *shape*, the conservation
laws and the fold structure, **not the published voltage digits**. The
pseudo-2D DFN is now implemented (`dfn`, Section 7); the intermediate rung used
for the fold comparison below is the electrolyte-enhanced SPM (SPMe).

Benchmarks (all tests; run < 1 s):

* **charge/mass conservation**: the negative-electrode volume-average drop is
  `i_app t/(F ε_s L)` and solid Li in the two particles is conserved, both to
  < 1e-9.
* **low-rate limit**: at C/1000 the terminal voltage equals the OCV to < 5 mV.
* `i_ref` (1C) = **51.02 A m⁻²** from the positive usable range Δx_p = 0.6;
  1C·1h = 51 020 mAh m⁻².

Capacity vs C-rate (discharge to 2.5 V), from the report run:

| C-rate | capacity (mAh m⁻²) | fraction of 1C·1h | V(0) | t_end (s) |
|---|---|---|---|---|
| 0.5 | 48 678 | 0.954 | 3.923 | 6870 |
| 1 | 47 899 | 0.939 | 3.846 | 3380 |
| 2 | 46 482 | 0.911 | 3.706 | 1640 |
| 3 | 43 789 | 0.858 | 3.577 | 1030 |

The fall with C-rate is the expected polarisation (ohmic + kinetics + solid
diffusion lag).

---

## 4. Folds and their indicators

### 4.1 SPM ← SPMe (electrolyte concentration)

SPMe restores a 1D electrolyte concentration across negative | separator |
positive with a uniform reaction current, plus the concentration overpotential
`(2RT/F)(1−t₊) ln(c_e⁺/c_e⁻)` (`spme`). The SPM validity indicator is the ratio
of the electrolyte diffusion time to the discharge time,

```
t_e / t_d,   t_e = L_e²/D_e ≈ 55.5 s for LG M50,   t_d = 3600/C.
```

Indices and errors (SPMe minus SPM, worst over the discharge):

| C-rate | indicator t_e/t_d | max |ΔV| SPMe−SPM | c_e at −/sep | c_e at sep/+ |
|---|---|---|---|---|---|
| 0.5 | 0.0077 | 0.68 mV | 1010.8 | 993.1 |
| 1 | 0.0154 | 1.37 mV | 1013.9 | 978.4 |
| 2 | 0.0309 | 2.76 mV | 1020.0 | 949.1 |
| 5 | 0.0771 | 7.16 mV | 1037.5 | 860.1 |

The correction is **linear in the indicator**, roughly `|ΔV| ≈ 0.09 V × (t_e/t_d)`.
With a 1 % voltage tolerance (~37 mV) the SPM fold is valid up to
`t_e/t_d ≈ 0.4`, i.e. tens of C for this parameter set: the LG M50 electrolyte
equilibrates in ~1 min, so **SPM does not need to unfold for terminal-voltage
accuracy at practical rates**. What SPMe still cannot represent — and where the
true DFN is required — is the **position dependence of the reaction**, whose
indicator is the reaction-nonuniformity number `a_s i_loc L²/(F D_e c_e)`; SPMe
assumes it is zero. That rung is now implemented as `dfn` (Section 7).

### 4.2 Thevenin RC ← SPM/DFN (`ecm`)

`V = U(z) − i R₀ − V_RC`, `V̇_RC = i/C₁ − V_RC/(R₁C₁)`, with the RC branch
integrated exactly. Parameters are extracted from a 1C/60 s SPM pulse response:
`R₀` from the instantaneous offset and `R₁, τ` from the remaining steady drop.

| C-rate | R₀ (Ω m²) | R₁ (Ω m²) | τ (s) |
|---|---|---|---|
| 1 | 0.003123 | 0.000747 | 30.17 |
| 3 | 0.002799 | 0.000875 | 32.50 |

Error of the ECM against the SPM that generated the pulse:

| C-rate | pulse 60 s | pulse 600 s |
|---|---|---|
| 1 | 5.2 mV | 170 mV |
| 3 | 21.1 mV | 474 mV |

The single RC branch (time constant `τ ≈ 30 s`) cannot follow the `√t`
solid-diffusion tail. The indicator is `t_pulse/τ`: the ECM is good to a few mV
for pulses ≲ τ, and **unfolds when the pulse is long compared with τ** — a
600 s (10τ) pulse already gives 4.6 % voltage error at 1C. This is the expected
boundary for a BMS-level model: use it for short transients and state
estimation, not for long galvanostatic discharges.

---

## 5. Thermal coupling

Bernardi's heat rate per unit electrode area
(D. Bernardi, E. Pawlikowski, J. Newman, *J. Electrochem. Soc.* **132**, 5
(1985)):

```
q = I(U − V) − I T dU/dT.
```

`thermal::total_heat` implements it; at 1C for the LG M50 SPM above,
`OCV = 4.005 V`, `V = 3.846 V`, so the irreversible term is
`I(U−V) = 8.13 W m⁻²`. With the sign-changing estimate
`dU/dT ∈ [1.35, 3.0]e-4 V K⁻¹` the reversible term is `−2.1 … −4.6 W m⁻²`
(cooling during this discharge), giving `q ≈ 6.1 W m⁻²` at z = 0.8. The
coefficient is an **ESTIMATE** and should be replaced by a measured `dU/dT(z)`.

**Interface to `heat`** (read-only, not modified): a continuum solver assembles

```
ρ c_p ∂T/∂t = ∇·(κ ∇T) + q_vol,    q_vol = q_areal / L   (`volumetric_source`),
```

using `heat::effective_kappa(T, L)` for `κ(T,L)` and its Fourier-validity flag.
The exchange is just the source density; no `heat` code changes are needed.

---

## 6. What is built, what is not, and where it unfolds

| rung | model | indicator | unfolds when | status |
|---|---|---|---|---|
| interface | linear / Tafel | `nF|η|/RT` | `|η| ≳ RT/nF` | done, exact limits tested |
| interface | Butler–Volmer | — | transport limit `i → i_L` | done, Koutecký–Levich |
| transport | Crank–Nicolson diffusion | `√(Dt)/L` | boundary reaches the domain edge | done, Cottrell + order tested |
| transport | Randles–Ševčík CV | `t` vs peak | beyond the peak / finite domain | done, −0.0018 % |
| solid | constant-flux sphere | shell thickness | thin boundary shell | done, series < 0.1 % |
| cell | SPM | `t_e/t_d` | `≈ 0.4` (tens of C here) | done, conservation + capacity |
| cell | SPMe | reaction nonuniformity | reaction is nonuniform | done (electrolyte only) |
| cell | DFN | reaction nonuniformity | (reference) | done (Section 7) |
| system | Thevenin RC | `t_pulse/τ` | `t_pulse ≳ few τ` | done, parameters extracted |
| thermal | Bernardi heat | — | — | done, `heat` interface documented |

**Negative results / limitations, stated plainly.**

* The positive-electrode OCP is a surrogate, not the published Chen 2020 fit, so
  the absolute voltage curve is not a validated cell prediction.
* The DFN is now implemented; it resolves the reaction distribution and the
  electrolyte potential, so it can certify SPM/SPMe (Section 7). Its positive
  OCP is still the surrogate, so its absolute voltage is not a published
  prediction.
* The kinetic rate constants, ohmic resistance, separator thickness, porosity,
  Bruggeman exponent and electrolyte diffusivity are flagged **ESTIMATE**;
  they set the SPMe magnitude, which is correspondingly uncertain.
* For the LG M50 parameter set the electrolyte is fast (`t_e ≈ 56 s`), so the
  SPM voltage error is only ~0.09 V × (t_e/t_d) — small through 5C. The fold
  demonstration is therefore conservative; cells with thicker electrodes or
  lower `D_e` unfold at lower rates.

All 44 crate tests pass in ~16 s (`cd fluid && cargo test -p electrochem`),
each well under the 30 s limit; the fine-grid full discharge is `#[ignore]`d.

---

## 7. The pseudo-2D DFN (reference rung)

`dfn` is the true reference rung: the full Doyle–Fuller–Newman model. The
through-thickness direction is a finite-volume grid over negative electrode |
separator | positive electrode; every electrode node carries a spherical
particle resolved radially (pseudo-2D). Concentrated-solution electrolyte
transport includes diffusion and migration through the transference number,
Bruggeman effective properties, and the concentrated-solution
diffusion-potential term; solid- and electrolyte-phase potentials come from
their own charge balances, and Butler–Volmer couples them at every node.

Governing equations (M. Doyle, T. F. Fuller, J. Newman, *J. Electrochem. Soc.*
**140**, 1526 (1993); F. Marquis et al., *J. Electrochem. Soc.* **166**, A3693
(2019)):

```
∂c_s/∂t = (D_s/r²) ∂_r(r² ∂_r c_s),      −D_s ∂_r c_s|_R = j/F
ε_e ∂c_e/∂t = ∂_x(D_e^eff ∂_x c_e) + (1−t₊) a_s j/F
∂_x i_e = a_s j,    i_e = −κ^eff ∂_x φ_e + κ_D^eff ∂_x ln c_e
∂_x(σ^eff ∂_x φ_s) = a_s j
j = i₀[exp(α_a Fη/RT) − exp(−α_c Fη/RT)],   η = φ_s − φ_e − U(c_s^surf)
i₀ = k F √(c_e c_s(c_max−c_s)),   D_e^eff = D_e ε_e^b,   κ^eff = κ ε_e^b,
κ_D^eff/κ^eff = 2RT(1−t₊)/F,   V = φ_s(L) − φ_s(0)
```

The electrolyte current `i_e` and the solid current `i_s = −σ^eff ∂_x φ_s`
satisfy `i_s + i_e = i_app`. **Time stepping is backward (implicit) Euler.
The coupled state `(φ_e, φ_s, j, c_e, c_s)` is solved by damped Newton**; the
Jacobian is formed by forward differences and the dense linear system is
row/column equilibrated before Gaussian elimination. `φ_e` is pinned at the
negative current collector, and the reaction kinetics fix the solid potential
through the per-electrode charge balance. Grid sizes are `(n_n, n_sep, n_p)`
through-thickness volumes and `m` radial shells per particle.

### 7.1 Laws the model must satisfy

| benchmark | result | error |
|---|---|---|
| solid lithium conservation (1C, 30 × 10 s) | invariant | ~1e-16 relative |
| electrolyte lithium conservation | invariant | ~1e-16 relative |
| charge balance `∫a_s j dx = i_app` per electrode | satisfied | < 1e-8 relative |
| reduction to SPM / SPMe at C/100 | 0.25 / 0.24 mV | linear in C-rate (ratio 1.99) |
| reduction at C/50 | 0.50 / 0.47 mV | linear in C-rate |
| spatial refinement (h ∝ 1/4, 1/6, 1/8) | 25.3, 18.7, 13.2 mV | ~ first order |
| temporal refinement (Δt = 40, 20, 10 s) | 1.28, 0.60, 0.26 mV | ~ first order (2.13, 2.33 per halving) |

The low-rate limit: at C/100 the DFN agrees with the crate's SPM/SPMe to
0.25 mV, and the difference scales linearly with C-rate (0.25 → 0.50 mV from
C/100 to C/50). The comparison removes the SPM/SPMe lumped `ESTIMATE`
`r_ohmic`, because the DFN resolves the solid- and electrolyte-phase ohmic
drops internally; keeping the lumped term would double-count it.

### 7.2 Fold table: SPM, SPMe, ECM against the DFN

Maximum terminal-voltage error over a **600 s constant-current window**
(at the same instant the same charge has passed in every model), and the ECM
error over the 60 s DFN pulse it was fitted to. `unif` is the per-electrode
reaction non-uniformity `std(a_s j)/|mean(a_s j)|` at `t = 100 s`. The SPM/SPMe
comparison uses `r_ohmic = 0` (see above).

| C-rate | V_DFN (600 s) | SPM err (mV) | SPMe err (mV) | unif | ECM 60 s err (mV) |
|---|---|---|---|---|---|
| 0.1 | 3.9794 | 2.74 | 2.60 | 0.188 | 0.61 |
| 0.5 | 3.8720 | 11.55 | 10.89 | 0.204 | 3.00 |
| 1 | 3.7234 | 20.41 | 19.13 | 0.247 | 6.28 |
| 2 | 3.3733 | 68.41 | 65.64 | 0.352 | 15.72 |
| 3 | 3.3349 | 135.86 | 131.78 | 0.598 | 36.35 |
| 5 | 2.2842 | 509.26 | 502.15 | 0.363 | 81.69 |

Reading the table against a 1 % voltage tolerance (~37 mV): SPM and SPMe hold
to ~1C (20 mV, 0.5 %); by 2C the error is 68 mV (1.8 %) and the fold has
failed; at 5C it is 0.51 V (13 %). SPMe is consistently a little better than
SPM but by only a few per cent, because for LG M50 the electrolyte is fast
(`t_e ≈ 56 s`) and the missing physics is the **reaction distribution**, which
SPMe does not resolve. The reaction non-uniformity grows with rate from 0.19
(C/10) to 0.60 (3C), consistent with the dimensionless indicator
`a_s i_loc L²/(F D_e c_e)`; at 5C the value is lower (0.36) because the state
is already deep into the discharge, but the *reaction distribution* is by then
strongly rate-limited. The ECM is best for short pulses — 0.6 mV at C/10,
6 mV at 1C, 82 mV at 5C — and unfolds as the pulse approaches the solid-
diffusion time, matching Section 4.2.

### 7.3 Validation against published data, and limitations

**No published DFN voltage curve was available offline, so the absolute
voltage is not validated against Chen et al. (2020) or a PyBaMM `Chen2020`
run.** The obstacle is not the DFN solver but the inputs: the crate's
positive-electrode OCP is a surrogate sigmoid, and the electrolyte
conductivity/diffusivity, porosities, Bruggeman exponent and solid
conductivities are flagged `ESTIMATE`. What *is* validated is the model's
structure and its laws: lithium conservation to round-off, exact charge
balance, the linear low-rate reduction to SPM/SPMe, first-order spatial and
temporal convergence, and the placement of the fold. As an independent
consistency check, a fine-grid (`12|4|12`, `m = 10`) 1C DFN discharge to 2.5 V
gives **48 324 mAh m⁻²**, within **0.9 %** of the crate's SPM (47 899 mAh m⁻²,
Section 2), the difference being the DFN's resolved ohmic drop versus the
SPM's lumped `R_Ω`.

Limitations: the dense Newton Jacobian is `O(n_unk²)` in memory and `O(n_unk³)`
to factor, so the fine-grid run is `#[ignore]`d and the tests use coarse grids;
a banded/block solver would be the next efficiency step. The solid potential is
resolved but contact resistance is not. Temperature is constant (isothermal).
The DFN is the reference only up to the material parameters it is fed.
