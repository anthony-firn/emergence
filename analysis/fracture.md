# The fracture rung: bond breaking → continuum toughness

This rung sits between the defect rung ([`defects.md`](defects.md)) and the
continuum elasto-plastic solid, and it is the fold/unfold pattern of
[`emergence_ladder.md`](emergence_ladder.md) applied to a crack:

- **coarse model** — linear-elastic fracture mechanics (LEFM): a universal
  near-tip field, an energy release rate `G = K²/E'`, and the Griffith
  criterion `G = 2 γ_s` (Griffith 1921);
- **validity indicator** — the ratio of the process-zone radius to the
  K-dominance / sample size, and the Rice (1992) cleavage-vs-emission
  competition `γ_s/γ_us`;
- **fold / unfold** — a 2D triangular LJ sheet with a pre-crack, loaded by
  the exact finite-crack K-field on the outer boundary; the closure for the
  coarse model is a cohesive-zone law with work of separation `2 γ_s`.

The module is `fluid/crates/materials/src/fracture.rs`; the benchmarks are
`fluid/crates/materials/tests/fracture.rs` (12 active tests, ≈24 s total,
plus one `#[ignore]`d size-sweep test run on request).

## 1. The Griffith fold

For a brittle crystal the critical energy release rate is `G_c = 2 γ_s`,
with `γ_s` the surface energy of the defects module. In plane strain the
mode-I toughness is

```text
K_IC = sqrt(2 γ_s E / (1 − ν²));
```

in plane stress and in a 2D monolayer, `K_IC = sqrt(2 γ_s E')` with the
appropriate modulus.

### 1.1 fcc Lennard-Jones crystal, 3D

The infinite-cutoff Cauchy–Born constants of LJ fcc at zero pressure are
`C11 = 105.182`, `C12 = C44 = 60.187 ε/σ³` (extrapolated, `analysis/materials.md`;
central-force Cauchy relation `C12 = C44`). The Voigt–Reuss–Hill average
gives `E = 103.165 ε/σ³`, `ν = 0.27131`; the Zener anisotropy is

```text
A = 2 C44 / (C11 − C12) = 2.675.
```

With the relaxed surface energies `γ(111) ≈ 3.09`, `γ(100) ≈ 3.19 ε/σ²`
([`defects.md`](defects.md)):

| crack | modulus | `K_IC` (`ε/σ^{3/2}`) | source of factor |
|---|---|---|---|
| `{100}`, isotropic (Hill) | `E, ν` | 26.66 | `√(2γ_s E/(1−ν²))` |
| `{111}`, isotropic (Hill) | `E, ν` | 26.23 | `√(2γ_s E/(1−ν²))` |
| `{100}` anisotropic | Stroh `H₂₂ = 0.020145` | **25.17** | `√(4γ_s/H₂₂)` |
| `{111}[1−10]` anisotropic | Stroh `H₂₂ = 0.017538`, `H₁₂ = 0.003045` | **26.55** | `√(4γ_s/H₂₂)` |

The general anisotropic factor is the Stroh Irwin matrix
`H = Re(i A B⁻¹)` (Stroh 1958; Sih, Paris & Irwin 1965; Ting 1996), so that
with the near-tip field `σ_{2i} = K_i/√(2π r)`

```text
G = ½ Kᵀ H K,     K_Ic = sqrt(4 γ_s / H₂₂).
```

The implementation solves the sextic `det[Q + p(R+Rᵀ) + p²T] = 0` by
Durand–Kerner iteration and builds the Stroh eigenvectors. It is validated
in the tests two independent ways:

- the antiplane entry is exact, `H₃₃ = 1/C44` (machine precision);
- for the cubic `{100}` crack the in-plane factor reproduces the
  Sih–Paris–Irwin orthotropic closed form
  `G/K² = (s₁₁s₂₂/2)^{1/2} [(s₂₂/s₁₁)^{1/2} + (2s₁₂+s₆₆)/(2s₁₁)]^{1/2}`
  (with the plane-strain reduced compliances) to `< 1e-12` relative.

The Cauchy relation makes LJ fcc strongly anisotropic (`A = 2.68`), yet the
`{100}` toughness is only **5.6 % below** the Hill-isotropic estimate
(`K_aniso/K_iso = 0.944`). The `{111}` crack has a non-zero Irwin coupling
`H₁₂`, i.e. pure mode-I loading induces a mode-II component; the mode-I
toughness is within 1.2 % of isotropic.

*Sources.* Griffith, *Phil. Trans. R. Soc. A* **221**, 163 (1921); Stroh,
*Phil. Mag.* **3**, 625 (1958); Sih, Paris & Irwin, *Int. J. Fract. Mech.*
**1**, 189 (1965); Ting, *Anisotropic Elasticity* (Oxford, 1996); Rice,
*J. Mech. Phys. Solids* **40**, 239 (1992).

## 2. Atomistic crack and lattice trapping

### 2.1 Model

A finite 2D triangular lattice of 12–6 Lennard-Jones atoms with **nearest
neighbour bonds only** (`r_max = 1.3 a₀`). At this cutoff the zero-stress
spacing is exactly the pair minimum `a₀ = 2^{1/6} σ`, so every bond is at
its force-free length and the `K = 0` state is stress-free to round-off
(forces `< 2 × 10⁻¹³`).

A half-row of `n_c` atoms is removed on `x ≤ 0` to make a mode-I pre-crack
of length `2a = n_c a₀`. The outer boundary is held at the **exact finite
centre-crack Westergaard field** (plane stress),

```text
Z(z) = σ z / √(z²−a²),   Z̄ = σ √(z²−a²),   K = σ √(π a),
u_x = [(1−ν) Re Z̄ − (1+ν) y Im Z]/E,
u_y = [2 Im Z̄ − (1+ν) y Re Z]/E,
```

evaluated with the Cauchy–Born in-plane moduli. This is the correct
boundary condition for a crack of the modelled length (the near-tip
asymptotic field alone is inconsistent with the finite crack mouth). Near
the tip (`r/a ≲ 10⁻³`) the field reduces to the universal mode-I K-field to
`3.2 × 10⁻⁴` relative (test `center_crack_field_reduces_to_the_near_tip_k_field`).

The 2D Cauchy–Born properties of the NN triangular LJ crystal are exact:

| quantity | value |
|---|---|
| `a₀` | `2^{1/6} σ = 1.1224620 σ` |
| `E_2D` | `65.9870 ε/σ²` |
| `ν_2D` | `1/3` (exact) |
| `γ_2D` (broken bond) | `|φ(a₀)|/a₀ = 1/a₀ = 0.890899 ε/σ` |
| `K_G = √(2 γ_2D E_2D)` | `10.8432 ε/σ^{3/2}` |

The broken-bond surface energy `1/a₀` is confirmed by the direct
crack-formation energy: extending the crack by one site raises the excess
energy by `2 γ_2D a₀` to 1e-10 (test
`two_d_surface_energy_and_griffith_prediction`).

### 2.2 Energy release rate by virtual crack extension

The physically meaningful atomistic test is the **energy release rate at
fixed grip**. Relax the reference crack (length `2a`) and the extended
crack (length `2a + a₀`) under the same boundary displacement field, and
form

```text
G = [E(a) − E(a + a₀)] / a₀.
```

`G < 0` means the pristine crack is stable; `G > 0` means extending it
lowers the total energy (surface energy included), i.e. the crack advances.
The crossing `G = 0` is the Griffith threshold.

#### Size extrapolation of the applied threshold

The applied K is larger than K_G because the boundary is loaded with the
sharp-crack field while the modelled crack is a one-row **blunt notch**, and
the imposed field uses the Cauchy–Born isotropic modulus. To separate a
finite-size effect from the notch conversion we repeated the measurement at
fixed crack length `2a = 10 a₀` for several sample half-widths `R = nx a₀`
(`R/a = 2 nx / 10`). Each `K_c` is the linear interpolation of `G(K)`
between the two bracketing loads on a grid of `G` computed by independent
relaxations at each K (a warm-started curve can hop basins and is not used):

| `nx` (boundary half-width, sites) | `R/a` | atoms | `K_c(applied)/K_G` |
|---|---|---|---|
| 18 | 3.6 | 1358 | 1.549 |
| 25 | 5.0 | 2590 | 1.497 |
| 37 | 7.4 | 5614 | 1.497 |

Above `K/K_G ≈ 1.6` the tip jumps by more than one site and the
one-site VCE is no longer the relevant measure, so only the crossing is
quoted. `K_c/K_G` converges to ≈1.50 by `R/a = 5` and is flat to
`R/a = 7.4`. A least-squares `A + B/(R/a)` fit over the three points gives
`A = 1.44`, but the smallest sample is not K-dominated and the two largest
already agree to `3 × 10⁻⁴`, so the plateau **1.497** is the more defensible
`R → ∞` estimate; the 1/R form is not the correct correction law. The
applied threshold therefore does **not** tend to `K_G`: the ≈1.5 factor is a
notch / boundary-condition conversion that is essentially independent of
sample size, not a finite-sample artefact and not a thermodynamic lattice
trap.

At the crossing the **effective tip K is exactly `K_G`**, by construction of
the energy balance: `G_net = 0` means the elastic driving force is exactly
the surface work `2 γ_s`, and `K_tip = √(2 γ_s E_2D)`. This is the atomistic
Griffith result — the energy balance closes — so there is no evidence here of
a thermodynamic shift of the advance threshold above `K_G`; the discrete
lattice enters only through the notch conversion and the dynamic (barrier)
question of §2.3. The fast test
`atomistic_crack_advances_above_griffith` asserts the sign change and
`K_c > K_G`; the size sweep is the ignored test `experiment_size_sweep`
(≈2.5 min).

### 2.3 Lattice trapping: local bond rupture

Loading and unloading quasistatically and recording the number of bonds
stretched past the LJ spinodal `r_s = (26/7)^{1/6} σ = 1.244455 σ` (where
`φ'' = 0`, so the bond has lost its stiffness) gives a hysteresis window.
For a 5-site crack in a 26-site-wide sheet (counting only bonds more than
two sites from the outer free surface):

| bound | `K/K_G` | definition |
|---|---|---|
| `K_+` | 0.964 | first tip bond reaches the spinodal on loading |
| `K_-` | 0.827 | that bond heals on unloading |
| `K_+/K_-` | 1.17 | metastable **local bond-rupture** window |

**This is not lattice trapping.** The tip bond snaps at `0.96 K_G`, far
below the thermodynamic advance threshold `1.497 K_G` of §2.2; the snapped
bond is a **metastable precursor**, not a crack advance, and it heals on
unloading below `0.83 K_G`. The real load-controlled trapping bounds
`K_-` and `K_+` are defined by the vanishing of the forward/backward energy
**barrier** for one-site advance at fixed load (Thomson, Hsieh & Rana 1971;
Bernstein & Hess 2003). A full saddle search along the tip reaction
coordinate was attempted here (harmonic-restraint path scan and
nudged-elastic-band refinement) but did **not** converge: the crack-tip
minimum-energy path in the chosen coordinate (`u_y` of the tip atom) was
rugged and basin-hopping, and no reproducible barrier-vs-K curve was
obtained. Those failed diagnostics were removed from the module. What the
present calculation does establish is only the *precursor* window: the
discrete crack neither heals nor advances over a finite band of loads, with
ratio `K_+/K_- ≈ 1.17`. The genuine trapping inequality
`K_- < K_G < K_+` remains an open question; it is **not** claimed here.

One thermodynamic ingredient that would have produced a *static* trapping
range has been ruled out: the zero-load relaxed surface energy is **exactly
linear in crack length** — successive one-site increments equal `2 γ_s a₀`
to `< 10⁻⁸` (test `relaxed_surface_energy_is_linear_in_crack_length`). So
the relaxed minima carry no lattice-periodic Peierls amplitude, and any
`K_±` window must come from the dynamic barrier, not from a difference of
the relaxed energies at successive tip positions.

The test suite exposes the whole `G(K)` history, so the thermodynamic
threshold (`critical_k_vce`) and the bond-rupture precursor
(`lattice_trapping`) are reported separately and never conflated.

### 2.4 K-field validation: a caveat

Comparing the relaxed interior displacement against the imposed field
globally is **not** a good test: the crack faces are a free surface with a
nearly zero-energy rigid-sliding mode, so the conjugate-gradient
displacement error stays large (`≈ 0.8` of the field) even after the
residual force falls to `10⁻⁵`. The field is instead validated
analytically (Westergaard reduces to the K-field, §2.1) and by the
energy balance (§2.2). This is recorded as an honest limitation of
displacement-based checks on free cracks.

## 3. Indicators and the cohesive-zone closure

### 3.1 Process zone and K-dominance

The Irwin plane-stress process-zone radius is

```text
r_p = K² / (2 π σ_y²),
```

and LEFM is valid only when `r_p` is small compared with the K-dominance
zone (a fraction of the smallest sample dimension). For the LJ sheet the
cohesive zone is only a few lattice sites, so the K-dominant zone is small;
this is exactly why the applied `K_c` in §2.2 converges to `≈1.50 K_G`
(a blunt-notch / boundary-field conversion) while the effective tip K is
`K_G`. The module provides `process_zone_radius` and
`lefm_valid(K, σ_y, L, threshold)`.

### 3.2 Rice brittle-vs-ductile criterion

Rice (1992) compares the energy for cleavage `G_cleave = 2 γ_s` with the
energy to emit a dislocation from the tip `G_emit ≈ 8 γ_us`
(`γ_us` the unstable stacking-fault energy of [`defects.md`](defects.md)).
The crystal cleaves when `2 γ_s < 8 γ_us`, i.e. `γ_us/γ_s < 1/4`.

For LJ fcc, with the relaxed `⟨112⟩` value `γ_us = 0.4226 ε/σ²` and
`γ_s(111) = 3.09 ε/σ²`:

```text
γ_us / γ_s = 0.137,     γ_s / γ_us = 7.31,
G_cleave = 6.18,        G_emit = 3.38   (ε/σ²),
```

so the criterion predicts **dislocation emission — ductile behaviour** for
LJ fcc. This is a genuine prediction of the pair potential, not a tuned
agreement: `γ_us` is only 14 % of `γ_s`. It is also *marginal*: the
often-quoted Rice threshold is `γ_s/γ_us ≈ 8–10`, and LJ sits at 7.3, just
on the emission side. Given that LJ fcc also has a **negative** intrinsic
stacking-fault energy (§3 of [`defects.md`](defects.md)), the pair potential
is not a good quantitative model of fcc ductility; the criterion is
reported as a consistent order-of-magnitude statement.

### 3.3 Cohesive-zone closure

The closure the coarse model uses where the atomistic rung is too fine is
the Rose–Smith–Ferrante universal binding-energy relation
(Rose, Smith & Ferrante, *Phys. Rev. Lett.* **47**, 675 (1981)):

```text
E(δ) = −E₀ (1 + δ/λ) e^{−δ/λ},     T(δ) = dE/dδ = E₀ (δ/λ²) e^{−δ/λ},
```

with the work of separation `E(∞) − E(0) = E₀ = 2 γ_s`, peak traction

```text
T_max = (2 γ_s/λ) e^{−1}   at   δ = λ.
```

The single parameter `λ` is the screening length (external closure). The
module provides `universal_binding`, `traction_from_energy`,
`cohesive_zone_universal` and `cohesive_work`; the tests integrate the
UBER traction numerically and recover `2 γ_s` and `T_max` to `10⁻⁶`.
By construction the cohesive law reproduces the Griffith toughness exactly,
so it is the correct continuation of the atomistic crack into the
continuum.

## Limits and honest misses

- **Stroh degeneracy.** The sextic root finder is ill-conditioned in the
  elastically isotropic limit (repeated roots `p = ±i`), where the two
  in-plane Stroh eigenvectors coalesce. The implementation is validated on
  anisotropic cubic and orthotropic cases; the isotropic limit is checked
  through the Sih–Paris–Irwin closed form, which matches to machine
  precision.
- **Atomistic K-dominance.** The applied thermodynamic threshold converges
to `K_c/K_G ≈ 1.497` as the sample radius grows (`R/a = 3.6, 5.0, 7.4`);
the residual `≈1.5` factor is the blunt-notch / boundary-field conversion,
not a finite-size artefact. The effective tip K at the crossing is exactly
`K_G` by the energy balance.
- **True lattice trapping.** The load-controlled barrier that defines
`K_-` and `K_+` (Thomson, Hsieh & Rana 1971) was **not** obtained: the
restrained-path and NEB saddle searches did not converge, and those
diagnostics were removed. The reported `K_+/K_- = 1.17` is a local
bond-rupture precursor, not the barrier-endpoint ratio, and does **not**
establish `K_- < K_G < K_+`.
- **2D model.** Nearest-neighbour LJ only, one atomic layer, plane stress
  (no out-of-plane coupling). The 3D fcc anisotropic toughness is computed
  separately from the Cauchy–Born constants and the defects surface
  energies; the two are not coupled in one simulation.
- **Notch.** The pre-crack is a removed row, i.e. a notch of width
  `2h = √3 a₀`; a lattice with zero-thickness crack faces would require
  bond cutting rather than atom removal.
- **{111} coupling.** The cubic `{111}[1−10]` crack has `H₁₂ ≠ 0`, so
  mode-I loading couples to mode II; only the diagonal `K_Ic` is quoted.
- **Rice factor.** The factor 8 in `G_emit` is Rice's sharp-crack estimate
  for the `{111}` slip geometry; the exact factor depends on the inclination
  and is not derived here.

## Constants and sources

- Griffith criterion: A. A. Griffith, *Phil. Trans. R. Soc. A* **221**, 163
  (1921).
- Stroh formalism and the Irwin matrix: D. M. Barnett & J. Lothe,
  *Phys. Status Solidi B* **67**, 105 (1975); T. C. T. Ting, *Anisotropic
  Elasticity: Theory and Applications* (Oxford, 1996).
- Anisotropic energy release rate: G. C. Sih, P. C. Paris & G. R. Irwin,
  *Int. J. Fract. Mech.* **1**, 189 (1965).
- Lattice trapping: R. Thomson, C. Hsieh & R. Rana, *J. Appl. Phys.* **42**,
  3154 (1971); N. Bernstein & D. W. Hess, *Phys. Rev. Lett.* **91**, 025501
  (2003) (brittle-lattice trapping via barrier calculations).
- Rice brittle-vs-ductile criterion: J. R. Rice, *J. Mech. Phys. Solids*
  **40**, 239 (1992).
- Universal binding energy: J. H. Rose, J. R. Smith & J. Ferrante,
  *Phys. Rev. Lett.* **47**, 675 (1981).
- LJ fcc equilibrium, cohesive energy, lattice sums and extrapolated
  Cauchy–Born constants `C11 = 105.182`, `C12 = C44 = 60.187`: Jones (1924);
  Lennard-Jones & Ingham (1925), as tabulated in
  [`materials.md`](materials.md).
- Surface energies `γ(111) ≈ 3.09`, `γ(100) ≈ 3.19`, unstable stacking-fault
  energy `γ_us(⟨112⟩) = 0.4226`: [`defects.md`](defects.md).
- Cauchy–Born in-plane moduli of the 2D triangular LJ lattice and the
  broken-bond surface energy: `materials::elastic` and
  `materials::defects` (this work; see the tests).
- Westergaard finite-crack field: H. M. Westergaard, *J. Appl. Mech.* **6**,
  A49 (1939); T. L. Anderson, *Fracture Mechanics* (3rd ed., CRC, 2005).
