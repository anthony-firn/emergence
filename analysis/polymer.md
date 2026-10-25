# Polymers and viscoelasticity: gloopiness from chains to continuum

Implemented in `fluid/crates/polymer` (new crate, no dependencies). This is
the viscoelastic rung of [emergence_ladder.md](emergence_ladder.md): the
fine model is a bead-spring chain, the coarse model is a viscoelastic
constitutive law (Maxwell, standard linear solid, Oldroyd-B), the fold is
the set of normal-mode relaxation times and the stress relaxation modulus
`G(t)` (plus its continuum sum, a Prony series), and the validity
indicator is the Deborah / Weissenberg number. The chain model itself is
unfolded again — from Hookean to FENE springs — when the flow is strong
enough to stretch the coil.

| module | content |
|---|---|
| `rouse.rs` | Rouse normal modes `τ_p`, `G(t)`, `η₀`, `G*(ω)`, and the fold of the Rouse spectrum onto a few-mode Maxwell/Prony model |
| `zimm.rs` | pre-averaged hydrodynamic interaction as a generalized eigenproblem, giving the Zimm `τ₁ ∝ N^{3/2}` |
| `bd.rs` | Brownian dynamics of a Rouse chain (Euler–Maruyama) and its Green–Kubo `G(t)` |
| `fene.rs` | Hookean dumbbell in extension (coil–stretch at `Wi = ½`) and the FENE-P / FENE chain fold |
| `continuum.rs` | Maxwell, standard linear solid, multi-mode Maxwell (Prony), Oldroyd-B, with exact start-up and oscillatory benchmarks |
| `rheology.rs` | virtual DMA/rheometer and time–temperature superposition (WLF) |
| `indicator.rs` | Deborah / Weissenberg regime selection and the Hookean-extension validity check |
| `gpr_fold.rs` | map a Maxwell mode onto the GPR strain-relaxation time and shear speed |

All numbers below come from `cargo test -p polymer` or
`cargo run -p polymer --example benchmarks`; reduced units
(`k_B T = b = ζ = 1`, one chain per unit volume) unless stated.

## 1. Rouse chain: normal modes and the stress modulus

A Rouse chain of `N` beads joined by `N-1` entropic springs of stiffness
`k = 3k_BT/b²` has Rouse-matrix eigenvalues
`λ_p = 4 sin²(pπ/2N)` and relaxation times
`τ_p = ζ/(kλ_p)`, `p = 1..N-1`. The stress modulus is
`G(t) = n k_B T Σ_p e^{-2t/τ_p}`. The factor 2 is not a convention: the
stress is bilinear in the Gaussian mode amplitudes `X_p`, whose
autocorrelation decays as `e^{-t/τ_p}`, so the Wick contraction gives
`e^{-2t/τ_p}`. Integrating gives the Rouse viscosity
`η₀ = n k_B T Σ_p τ_p/2 = n ζ N² b² / 36` for large `N`.

| check, `N = 50` | measured | reference | error |
|---|---|---|---|
| `τ₁` exact discrete | 84.4621 | `ζ N² b²/(3π² k_BT)` = 84.4343 | 0.03% |
| `η₀` exact finite sum | 69.4167 | `n ζ N² b²/36` = 69.4444 | 0.04% |
| `G(0) = n k_BT(N-1)` | 49 | 49 | exact |
| Rouse high-frequency slope `d ln G''/d ln ω` | 0.4996 | ½ | 0.1% |

`G(t) = n k_BT Σ e^{-2t/τ_p}` and `G''(ω)/ω → η₀` are asserted in the
tests. The high-frequency slope is measured between the slowest and
fastest mode of an `N = 200` chain.

### Folding the Rouse spectrum onto a few Maxwell modes

`Rouse::prony_lumped(m)` keeps the `m` slowest Rouse modes exactly and
lumps all faster modes into one Maxwell mode that preserves both the
viscosity `η = Σ g_i λ_i` and the first moment `Σ g_i λ_i²`. Viscosity and
the terminal moduli are therefore exact for every `m`; the error is
concentrated where the Rouse spectrum is non-Maxwellian.

| modes kept | `η₀` error | max `|ΔG|/G` on `[0.1,10]τ₁` | max `|ΔG''|/G''`, `ωτ₁∈[10⁻³,10⁻¹]` | max `|ΔG''|/G''`, `ωτ₁∈[10⁻³,10³]` |
|---|---|---|---|---|
| 2 (m=1) | 0 | 26.4% | 1e-9 | 99.5% |
| 3 (m=2) | 0 | 5.46% | 1e-10 | 98.2% |
| 4 (m=3) | 0 | 2.88% | 2e-11 | 95.5% |
| 6 (m=5) | 0 | 0.059% | 2e-12 | 84.7% |

The lesson is the one the ladder predicts: a few relaxation modes capture
the terminal and time-domain response (a 4-mode Maxwell reproduces the
stress relaxation to 3%), but the Rouse `ω^{1/2}` high-frequency regime is
intrinsically *not* a sum of a few Maxwell modes, so the high-band `G''`
error stays tens of per cent until many modes are kept. The fold is
valid where its indicator — a broad separation between the slowest mode
and the band edge — holds.

## 2. Brownian dynamics and its Green–Kubo modulus

`bd.rs` integrates the overdamped SDE
`ζ dR_i = F_i dt + √(2k_BTζ) dW_i` with Euler–Maruyama and computes the
Kramers shear stress `S_xy = Σ_bonds k r_x r_y`. The estimator is
`G(t)/(n k_BT) = ⟨S_xy(t)S_xy(0)⟩/(k_BT)²`, one `k_BT` per spring.

| check, `N = 8`, 200 chains, `dt = 0.005`, `T = 180` | measured | reference |
|---|---|---|
| `G(0)` | 7.208 | `N-1 = 7` (3.0%) |
| max rel. error of `G(t)` where `G_an > 0.3 k_BT n` | 9.8% at `t = 1.30` | analytic Rouse sum |
| end-to-end `⟨R²⟩` | 0.9979 `(N-1)b²` | `(N-1)b²` |
| end-to-end long-time log-slope | 0.465 | `1/τ₁ = 0.457` (1.8%) |

The `G(t)` curve matches the analytic Rouse sum to 10% over the
well-resolved range; beyond `t ≈ 2τ₁` the finite window `T` leaves a
positive bias of order `Σ_p τ_p/T ≈ 0.02` that we do not tune away.

**A pitfall worth recording.** Subtracting the per-trajectory sample mean
from `S_xy` biases the long-time tail strongly negative (the bilinear
stress of a slowly relaxing mode has a poorly determined sample mean over
a few `τ₁`). The same effect appears if trajectories are shorter than a
few `τ₁`. The test uses long trajectories and no mean subtraction, and
states the residual bias.

## 3. Zimm: pre-averaged hydrodynamic interaction

`zimm.rs` builds the Rouse matrix `A` and the dimensionless friction
matrix `H_ij = δ_ij + h*/√|i-j|`, and solves the mobility eigenproblem
`(H A)φ = μφ` (symmetrised as `Lᵀ A L`, `H = LLᵀ`). The relaxation times
are `τ_p = ζ/(k μ_p)`. For `h* = 0` the code reproduces the Rouse times
to `1e-9` (test).

| `h*` | `τ₁(40)/τ₁(20)` exponent | `τ₁(80)/τ₁(40)` | `τ₁(160)/τ₁(80)` |
|---|---|---|---|
| 0 (Rouse) | 2.0 | — | 2.00 |
| 0.25 | 1.58 | 1.56 | 1.54 |

The longest time crosses over from `N²` to the Zimm `N^{3/2}` (i.e.
`N^{3ν}` with `ν = ½`); in a good solvent, excluded-volume swelling gives
`ν ≈ 0.588` and `N^{1.76}`. The measured exponent approaches 1.5 from
above as `N` grows, the expected finite-size drift.

## 4. Continuum folds: Maxwell to Oldroyd-B

Oldroyd-B is the *exact* continuum fold of a dilute Hookean-dumbbell
solution: one Maxwell mode with modulus `G = n k_BT`, time `λ`, polymer
viscosity `η_p = Gλ`, plus solvent `η_s`.

| benchmark | measured | exact |
|---|---|---|
| UCM start-up `σ₁₂(t) = η γ̇(1-e^{-t/λ})` | — | asserted |
| UCM start-up `N₁(t) = 2η_pλγ̇²[1-(1+t/λ)e^{-t/λ}]` | short-time `∝ t²`, steady | asserted |
| Oldroyd-B `N₁` steady, `η_s=0.1, η_p=0.9, λ=2, γ̇=0.3` | 0.324000 | `2η_pλγ̇²` = 0.324000 |
| Oldroyd-B `σ₁₂` steady | 0.300000 | `η₀γ̇` = 0.300000 |
| Maxwell terminal `G'/ω²` | 0.500000 | `Gλ²` = 0.5 |
| Maxwell terminal `G''/ω` | 1.000000 | `η₀ = Gλ` = 1 |

The Oldroyd-B identity `N₁ = 2λσ₁₂,p²/η_p` is asserted at steady state.
The start-up `N₁` has the exact form `2η_pλγ̇²[1-(1+t/λ)e^{-t/λ}]`, not
the commonly mis-quoted `(1-e^{-t/λ})²`.

## 5. GPR fold and the factor of six

In the GPR unified model (`compressible/src/gpr.rs`) the elastic shear
modulus is `G = ρc_s²` and the `De ≪ 1` emergent viscosity is
`μ = ρc_s²τ/6`. To represent a Maxwell mode `(G, λ)` set
`ρc_s² = G` and `τ_GPR = 6λ`; then `μ = ρc_s²τ/6 = Gλ = η`. The linear
GPR shear strain decays as `e^{-6t/τ}`, i.e. with Maxwell time `τ/6` — the
origin of the factor 6. A naive `τ = λ` mapping would put GPR's viscosity
at `η/6`. Example: `G = 2.5×10⁵ Pa`, `λ = 10⁻³ s`, `ρ = 1000 kg/m³` gives
`c_s = 15.81 m/s`, `τ = 6×10⁻³ s`, `μ = 250 Pa s = η`.

## 6. Finite extensibility: the coil–stretch transition

A Hookean dumbbell in steady uniaxial extension has
`c_xx = ⟨Q_x²⟩/⟨Q_x²⟩_eq = 1/(1-2Wi)`, diverging at `Wi = λε̇ = ½`; the
Oldroyd-B elongational viscosity
`η_E = 3η_s + 3η_p/[(1-2Wi)(1+Wi)]` diverges there. Brownian dynamics of
the dumbbell reproduces the Hookean stretch (`Wi = 0.3`: measured 2.5 in
`c_xx`, analytic 2.5, test tolerance 8%). Above the transition a FENE
spring is required.

| FENE, `b = H Q₀²/k_BT = 50`, `Wi = 0.8` | result |
|---|---|
| FENE-P closure `c_xx` | 20.07 |
| FENE-P trace `tr(A)/A_eq` | 20.93 `< b` |
| FENE Brownian dynamics `c_xx` | 14.84 |
| FENE BD `⟨Q²⟩` | 5.24 (finite, below `Q₀² = b/H = 16.7`) |

The FENE-P Peterlin closure over-predicts the stretch by ~25% at this
`Wi`; the BD is the reference. The essential point — no divergence above
`Wi = ½`, and the extension saturates below `Q₀²` — holds.

## 7. Virtual DMA, WLF and time–temperature superposition

`Rheometer::oscillatory_from_ode` integrates the diagonal Maxwell ODEs
under `γ(t) = γ₀ sin ωt` with an exact (piecewise-linear) integrating
factor and extracts `G', G''` by projecting the steady stress on `sin` and
`cos`.

| check | measured | reference |
|---|---|---|
| virtual DMA `G', G''` vs analytic, ω = 0.1, 1, 10 | max rel. error 5.9e-4 | `G*(ω)` |
| step stress relaxation `σ = γ₀G(t)` | exact | asserted |
| WLF `log₁₀ a_T` at `T_g + 50 K` | −8.5827 | universal constants |
| TTS collapse `G*(ω;T) = G*(a_Tω;T_ref)` | 3.8e-16 | exact by construction |

The WLF shift is `log₁₀ a_T = −C₁(T−T_ref)/(C₂+T−T_ref)` with the
universal `C₁ = 17.44`, `C₂ = 51.6 K` relative to `T_g` (Williams,
Landel & Ferry 1955). At `T_g + 50 K` this gives `a_T = 2.6×10⁻⁹`.

Polystyrene-like melt parameters used for the DMA/TTS tests:
`T_g = 373 K`, plateau modulus `G_N⁰ ≈ 2.0×10⁵ Pa`, entanglement molar
mass `M_e ≈ 13.3 kg/mol` (Ferry 1980; Brandrup & Immergut; Fetters et
al. 1994). The spectrum is Rouse-like (`λ_p = λ_1/p²`, `λ_1 = 0.1 s` is
an estimate at ~453 K); its absolute time scale is *not* a sourced number
and appears only in the shape/`G_N⁰` and in the WLF collapse, which is
exact. The resulting zero-shear viscosity `≈ 5.9×10³ Pa s` is an estimate.

## 8. Indicator: where the coarse model is valid

`indicator.rs` computes `De = λ/t_obs` (or `λω`) and `Wi = λγ̇` and
selects a regime: `De < 0.1` Newtonian; `0.1 < De < 10` viscoelastic;
`De > 10` elastic solid. For steady shear `t_obs = 1/γ̇`, so `De = Wi`.

The chain-level validity indicator is independent: a Hookean
(Oldroyd-B) closure has no steady state in strong extension beyond
`Wi = ½`, so the indicator unfolds to FENE there. This is the rung's
"cheap check": it uses only `λ`, `t_obs` and `γ̇`, already known to the
continuum solver.

## Honest limitations

- The Rouse and Zimm models are free-draining / pre-averaged; real melts
  need entanglements (reptation), which are not modelled.
- The Zimm exponent is measured at small `N` (20–160) and approaches 1.5
  slowly; excluded-volume swelling (`N^{1.76}`) is stated but not
  simulated.
- The BD Green–Kubo test is stochastic; the tolerance (10% where
  `G > 0.3 n k_BT`) is a statement about the well-resolved part of the
  curve. The long-time tail carries a finite-window bias that is reported,
  not tuned away.
- FENE-P over-predicts the stretch by ~25% at `Wi = 0.8`, `b = 50`
  relative to FENE Brownian dynamics. No closure is exact.
- The polystyrene spectrum is a shape surrogate with estimated time
  constants; only the plateau modulus, `T_g` and the exact WLF/TTS
  relation are quantitative.

## Sources

- P. E. Rouse, *J. Chem. Phys.* **21**, 1272 (1953) — Rouse model.
- B. H. Zimm, *J. Chem. Phys.* **24**, 269 (1956) — hydrodynamic
  interaction, `N^{3/2}`.
- M. Doi & S. F. Edwards, *The Theory of Polymer Dynamics* (1986),
  ch. 4 — normal modes, `G(t)`.
- M. Rubinstein & R. H. Colby, *Polymer Physics* (2003), ch. 8 — Rouse
  and Zimm spectra, plateau modulus.
- R. B. Bird, R. C. Armstrong & O. Hassager, *Dynamics of Polymeric
  Liquids*, Vol. 1, 2nd ed. (1987), ch. 8 — Maxwell/UCM/Oldroyd-B
  start-up and oscillatory solutions.
- R. B. Bird, C. F. Curtiss, R. C. Armstrong & O. Hassager, *Dynamics of
  Polymeric Liquids*, Vol. 2, 2nd ed. (1987), ch. 13 — dumbbells and
  FENE-P.
- J. G. Oldroyd, *Proc. R. Soc. A* **200**, 523 (1950) — Oldroyd-B.
- P. G. de Gennes, *J. Chem. Phys.* **60**, 5030 (1974) — coil–stretch
  transition at `Wi = ½`.
- H. R. Warner, *Ind. Eng. Chem. Fundam.* **11**, 379 (1972) — FENE
  spring; A. Peterlin, *Makromol. Chem.* **44**, 338 (1961) — Peterlin
  closure.
- M. L. Williams, R. F. Landel & J. D. Ferry, *J. Am. Chem. Soc.* **77**,
  3701 (1955) — WLF equation, `C₁ = 17.44`, `C₂ = 51.6 K`.
- J. D. Ferry, *Viscoelastic Properties of Polymers*, 3rd ed. (1980) —
  polystyrene plateau modulus and shift factors.
- J. Brandrup & E. H. Immergut, *Polymer Handbook* — polystyrene `T_g`.
- L. J. Fetters, D. J. Lohse, D. Richter, T. A. Witten & A. Zirkel,
  *Macromolecules* **27**, 4639 (1994) — entanglement molar mass.
- I. Peshkov & E. Romenski, *J. Comput. Phys.* (2016); M. Dumbser,
  I. Peshkov, E. Romenski & O. Zanotti (2016) — GPR unified continuum
  model; `compressible/src/gpr.rs`.
