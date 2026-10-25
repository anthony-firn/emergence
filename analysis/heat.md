# Phonon heat transport: the diffusive → ballistic fold

Implemented in `fluid/crates/heat`. This is the phonon rung of the heat
ladder. The coarse model is the continuum heat equation with a conductivity
`κ(T, L)`; the fine model is the phonon Boltzmann equation; the fold between
them is controlled by a single dimensionless indicator, the **phonon Knudsen
number** `Kn = Λ/L`, where `Λ` is the spectral mean free path and `L` the
sample size.

```text
  continuum heat equation  ρc ∂T/∂t = ∇·(κ(T,L) ∇T)   <-- coarse model
           ^   indicator: Kn = Λ/L  (phonon Knudsen number)
  kinetic (Callaway/Holland) κ(T), boundary scattering, ballistic Landauer
           ^   validation: grey BTE between two diffuse walls
  discrete ordinates / exact integral Boltzmann equation
```

- **coarse model** — a `∇·(κ∇T)` heat equation with a size-limited `κ(T, L)`,
  as already used by `semiconductor/src/thermal.rs` (which today takes
  `k_si(T)` from a power-law fit);
- **validity indicator** — `Kn = Λ/L` from the kinetic model; `Kn ≪ 1` is
  Fourier, `Kn ≫ 1` is ballistic;
- **fold / unfold** — fold the boundary scattering into the relaxation time
  (`τ_B⁻¹ = v/L`, Casimir) while `Kn` is small; when `Kn ≳ 1` unfold to the
  grey BTE, whose exact solution gives the effective conductivity in between.

The crate is `fluid/crates/heat`, with no dependencies. Modules: `phonon`,
`scattering`, `kinetic`, `bte1d`, `materials`, `effective`.

## 1. Kinetic (Callaway/Holland) conductivity

The relaxation-time expression used is exactly

\[
  \kappa = \frac{1}{3}\sum_{\text{branches}} g \int C(\omega)\,v_g(\omega)^2\,
           \tau(\omega)\,g_{\text{DOS}}(\omega)\,d\omega,
\]

with `C(ω) = k_B x² e^x/(e^x−1)²`, `x = ħω/k_BT`, and the Debye density of
states `g_DOS = q²/(2π² v_g)`. Two per-branch dispersions are provided:

- **Debye**: `ω = v q`, branch cutoff `ω_max = k_B Θ/ħ`;
- **Born–von Kármán sine**: `ω = ω_max sin(πq/2q_max)` with `q_max` chosen so
  `dω/dq|₀ = v`; the group velocity falls to zero at the zone boundary.

The scattering rates (`scattering.rs`) are the Callaway/Holland set:

| process | rate | source |
|---|---|---|
| umklapp | `τ_U⁻¹ = ħγ²ω²T /(M v²Θ) · e^{−Θ/3T}` | Morelli et al. 2002 |
| normal | `τ_N⁻¹ = β τ_U⁻¹` | Callaway 1959 |
| point defect (isotope) | `τ_I⁻¹ = V₀ Γ ω⁴ /(4π v³)` | Klemens 1958 |
| boundary (Casimir) | `τ_B⁻¹ = v/L` | Casimir 1938 |

### Silicon parameters

Branch-resolved set of **Morelli, Heremans & Slack, *Phys. Rev. B* 66, 195304
(2002)**, as tabulated by **Kazan et al., *J. Appl. Phys.* 107, 083503
(2010)**, Table II:

| quantity | value |
|---|---|
| `v_T`, `v_L` | 5840, 8430 m/s |
| `Θ_T`, `Θ_L` | 216, 566 K |
| `γ_T`, `γ_L` | 0.45, 1.09 |
| `M` | 4.66×10⁻²⁶ kg |
| `β` | 1.5 |
| `V₀` | 2.0×10⁻²⁹ m³ (a³/8, a = 0.5431 nm) |
| reference sample | 14 mm (Kazan et al.) |

These are **fitted** parameters. Agreement with the measured κ(T) is a
consistency check, not a prediction. The isotope term is *not* fitted: Γ is
computed from the isotopic abundances,
`Γ = Σ fᵢ(1 − mᵢ/M̄)² = 2.0088×10⁻⁴` for natural Si. This matches the `g² =
200.7×10⁻⁶` quoted by Inyushkin et al., *J. Appl. Phys.* 123, 095104 (2018),
to 0.09 %.

### Benchmarks

All numbers from the checked-in tests (`tests/benchmarks.rs`).

| quantity | reference | computed | error |
|---|---|---|---|
| `Γ`(natural Si) | 200.7×10⁻⁶ (Inyushkin 2018) | 2.0088×10⁻⁴ | 0.09 % |
| κ(300 K), RTA | 148 W/(m K) (Glassbrenner & Slack 1964) | 126.96 | −14.2 % |
| κ(300 K), Callaway | 148 | 218.00 | +47.3 % |
| natural-Si peak (14 mm) | few ×10³ W/(m K) near 20–30 K | 6901 W/(m K) at 17.25 K | location low, value high |
| ²⁸Si/nat at 300 K | +8 ± 1 % (Inyushkin 2018) | +2.3 % | model 3–4× too small |
| ²⁸Si/nat at the peak | ≈ 10 at 24 K (99.995 %) | 4.18 at 23.8 K | model ≈2.4× too small |
| κ(300 K), sine vs Debye | — | 113.6 vs 127.0 | sine −10.5 % |

The RTA value uses the combined rate including normal processes; this is the
formula quoted in the task. The Callaway normal-drift correction is provided
as well (`kinetic::kappa_callaway`) but it **overshoots** (218 vs 148). The
isotope effect is a genuine *prediction* (only `Γ` changes, and `Γ` is
computed, not fitted) and it is under-predicted by a factor 3–4. Both misses
are honest and are not tuned.

## 2. Size effects: nanowires and the Casimir fold

Setting `τ_B⁻¹ = v/d` (diameter `d`) folds boundary scattering into the
kinetic integral. Measured Li et al., *Appl. Phys. Lett.* 83, 2934 (2003)
room-temperature values for individual smooth Si nanowires (read from
Fig. 3; approximate):

| diameter | Li et al. 2003 | model | error |
|---|---|---|---|
| 22 nm | ≈ 9 W/(m K) | 8.78 | −2.4 % |
| 37 nm | ≈ 15 | 13.43 | −10.5 % |
| 56 nm | ≈ 20 | 18.33 | −8.3 % |
| 115 nm | ≈ 35 | 29.46 | −15.8 % |

The bulk 300 K mean free path is `Λ = 246 nm`, so these wires are in the
**transitional** regime (`Kn = 0.09…0.47`). The agreement is quantitative for
the two thinnest wires and within 16 % for the largest; the simple Casimir
`v/d` rate is a grey approximation and does not capture the full geometry.

## 3. The grey Boltzmann equation: exact fold

### Set-up and exact solution

For a grey gas (`C`, `v`, `τ` frequency-independent) between two diffuse
walls at `T₁, T₂`, the steady linearised RTA equation is

\[
  \mu\,\frac{\partial i}{\partial x} = \frac{i_0 - i}{\Lambda},
  \qquad i_0(x) = \frac{1}{2}\int_{-1}^{1} i(x,\mu)\,d\mu,
\]

with `Λ = vτ` and wall emission `i(0,μ>0)=CT₁/2`, `i(L,μ<0)=CT₂/2`. Energy
conservation fixes the local equilibrium `i₀`. Integrating along
characteristics gives the exact Fredholm equation

\[
  i_0(x) = \frac{1}{2}\Big[i_1 E_2(x/\Lambda) + i_2 E_2((L-x)/\Lambda)\Big]
         + \frac{1}{2\Lambda}\int_0^L i_0(x')\,E_1(|x-x'|/\Lambda)\,dx',
\]

with `E_n(z) = ∫₀¹ μ^{n-2} e^{−z/μ} dμ` the standard exponential integrals.
The flux follows from the same solution,
`q = i₁/2 − i₂E₃(L/Λ) − (1/Λ)∫₀^L i₀E₂(x/Λ)dx`. The crate solves this
equation by a Nyström method with piecewise-linear `i₀` and **exact** cell
integrals (using `dE₂/dz = −E₁` etc.), which is the independent benchmark.

The dimensionless flux is `q/(C v ΔT)` and the effective conductivity ratio is

\[
  \frac{\kappa_{\text{eff}}}{\kappa_0} = \frac{3q}{2\,\mathrm{Kn}},
  \qquad \kappa_0 = \tfrac13 C v \Lambda,
\]

with `Kn = Λ/L`. The exact limits are

- **Fourier** (`Kn → 0`): `κ_eff/κ₀ → 1`, with the leading temperature-jump
  correction `1/(1 + 2ζ Kn)`, `ζ = 0.7104` the Milne extrapolation length;
- **ballistic** (`Kn → ∞`): `κ_eff/κ₀ → (3/4)/Kn`, i.e. the Casimir/Landauer
  flux `q = ¼ C v ΔT`.

### Benchmarks of the fold

| `Kn` | exact integral eq. | DOM (24×160) | relative | slip/ballistic |
|---|---|---|---|---|
| 0.01 | 0.97562 | — | — | — |
| 0.05 | 0.93303 | — | — | 0.93367 (slip, 0.07 %) |
| 0.1 | 0.87542 | 0.87088 | 4.7×10⁻³ | 0.87559 (slip, 0.02 %) |
| 0.2 | 0.77858 | 0.77500 | 4.6×10⁻³ | — |
| 0.5 | 0.58508 | 0.58302 | 3.5×10⁻³ | — |
| 1 | 0.41505 | 0.41405 | 2.4×10⁻³ | — |
| 2 | 0.26406 | 0.26370 | 1.4×10⁻³ | — |
| 5 | 0.12738 | 0.12732 | 4.7×10⁻⁴ | — |
| 10 | 0.06868 | 0.06868 | 4×10⁻⁵ | — |
| 100 | 0.00743 | 0.00744 | 1.2×10⁻³ | 0.007500 (ballistic, 1.0 %) |
| 1000 | 0.000749 | — | — | 0.000750 (ballistic, 0.1 %) |

The discrete-ordinates solver (first-order upwind, block-tridiagonal
elimination) agrees with the exact integral equation to better than 0.6 %
over `0.1 ≤ Kn ≤ 100`. The exact solution interpolates between the two limits
and is **not** a simple rational function; the interpolation
`1/(1 + (4/3)Kn)` happens to be within ≈3 % and is used only as a diagnostic,
not as the benchmark.

### The fold in numbers (silicon)

At 300 K the bulk mean free path is `Λ = 246 nm`:

| sample `L` | `Kn` | regime | Fourier valid? |
|---|---|---|---|
| 1 mm | 2.5×10⁻⁴ | Fourier | yes |
| 100 nm | 2.46 | transitional | no |
| 10 nm | 24.6 | ballistic | no |

For the 10 nm case the Landauer value `¼ C v L = 3.56 W/(m K)` and the
boundary-limited Callaway model gives 4.36 (22 % high — the relaxation-time
fold is a grey approximation).

## 4. `κ(T, L)` for a continuum heat solver

`effective::effective_kappa(branches, params, T, L)` returns

```text
EffectiveKappa {
    kappa_bulk, kappa_effective, mean_free_path,
    knudsen, regime, fourier_valid,
}
```

`conductivity()` is the coefficient to hand to `∇·(κ∇T)`. The rule for a
solver is:

- `fourier_valid == true` (`Kn < 0.1`): use `kappa_effective`; the Fourier
  closure is accurate to the boundary-scattering correction already folded in;
- `fourier_valid == false`: the diffusion closure is **not** certified; unfold
  to the grey BTE (`bte1d`) or apply the ballistic closure
  `κ ≈ ¼ C v L` (`kappa_ballistic`).

`silicon_effective_kappa(T, L)` is the convenience wrapper for silicon. The
existing `semiconductor/src/thermal.rs` uses `k_si(T)` (a power-law fit to
Glassbrenner–Slack); it can be replaced by this call without changing the
solver, but that crate is deliberately not modified here.

## Limits and honest misses

- **300 K silicon is not predicted.** The published Morelli/Kazan parameters
  give 127 W/(m K) in the RTA (−14 % vs 148) and 218 with the Callaway
  correction (+47 %). The parameters are a fit to a particular data set and
  branch-resolved Debye model; the discrepancy is reported, not tuned away.
- **The isotope effect is under-predicted** (a genuine prediction, since only
  the computed Γ changes): +2.3 % at 300 K against +8 ± 1 % measured, and a
  peak ratio 4.2 against ≈10. The Klemens point-defect rate with the physical
  Γ is too weak at the thermal frequencies; a larger effective Γ (or the full
  Callaway drift with the true spectral weights) would be needed.
- **The natural-Si peak location is low** (17 K vs 25–30 K measured) and the
  peak height is high (6900 vs ≈5000 W/(m K)). The low-temperature peak is
  very sensitive to the chosen Casimir length (14 mm here) and to the
  umklapp cutoff temperature.
- **The grey BTE benchmark is a model problem.** It assumes a
  frequency-independent `τ`, isotropic dispersion and perfectly diffuse walls.
  Real nanowires have specular/rough surfaces and a frequency-dependent `τ`;
  the Li et al. comparison is therefore only quantitative to ≈ 16 %.
- **The exact integral equation is converged only while the grid resolves the
  `O(Kn)` boundary layer** (`m ≳ 20/Kn`); for `Kn < 0.05` a boundary-layer
  grid or a spectral method would be needed. The tests use `m = 200–400` and
  `Kn ≥ 0.05`.
- **Germanium is an estimate.** Its sound speeds and branch split are
  plausible but not a published fit; it is included for completeness and is
  flagged in the source.

## Constants and sources

- Phonon umklapp form `τ_U⁻¹ ∝ γ²ω²T e^{−Θ/3T}`: Morelli, Heremans & Slack,
  *Phys. Rev. B* **66**, 195304 (2002); Kazan, Guisbiers, Pereira, Correia,
  Masri, Bruyant, Volz & Royer, *J. Appl. Phys.* **107**, 083503 (2010).
- Callaway model with normal processes: Callaway, *Phys. Rev.* **113**, 1046
  (1959).
- Point-defect (isotope) rate: Klemens, *Solid State Physics* **7**, 1
  (1958).
- Boundary (Casimir) limit: Casimir, *Physica* **5**, 495 (1938); Landauer,
  *J. Phys. Chem. Solids* **49**, 115 (1988).
- Measured silicon κ(T): Glassbrenner & Slack, *Phys. Rev.* **134**, A1058
  (1964) (≈148 W/(m K) at 300 K).
- Si nanowire κ(d): Li, Wu, Kim, Shi, Yang & Majumdar, *Appl. Phys. Lett.*
  **83**, 2934 (2003).
- Isotope effect: Ruf, Henn, Asen-Palmer, Gmelin, Cardona, Pohl, Devyatych &
  Sennikov, *Solid State Commun.* **115**, 243 (2000); revised by
  Inyushkin, Taldenkov, Gibin, Gusev & Pohl, *Phys. Rev. B* **70**, 045205
  (2004) and Inyushkin, Taldenkov, Ager, Haller, ... *J. Appl. Phys.* **123**,
  095104 (2018); mass fluctuation `g² = 200.7×10⁻⁶`.
- Milne extrapolation length `ζ = 0.710446…`: Case & Zweifel, *Linear
  Transport Theory* (1967); Chandrasekhar, *Radiative Transfer* (1950).
- Exponential integrals `E_n`: Abramowitz & Stegun, *Handbook of Mathematical
  Functions*, ch. 5.
- Silicon sound speeds / Debye temperature cross-check: Hopcroft, Nix & Kenny,
  *J. Microelectromech. Syst.* **19**, 229 (2010) (in `materials` crate).
