# The virtual materials-testing lab

Implemented in `fluid/crates/testlab` (depends read-only on `materials`,
`md`, `compressible`, `heat`, `polymer`, `magnetism`, `photonics`, `wetting`
and `granular`). This is the standard material test programme run on one
material at every rung of the materials ladder in
[emergence_ladder.md](emergence_ladder.md), with the fold between rungs
checked where it exists and flagged where it does not.

```sh
cd fluid
cargo test -p testlab                              # 18 tests, ~11 s total
cargo run -p testlab --example report              # full standard-test battery
cargo run -p testlab --example report -- --full    # plus MD and plate impact
```

The battery covers mechanical (fracture, yield), thermal (heat capacity,
thermal expansion, conductivity), viscoelastic (DMA, WLF), magnetic
(hysteresis, domain walls), optical (reflectance, Reststrahlen) and
interfacial/granular (contact angle, `μ(I)`) tests. The report prints one
table per material and a cross-material summary of which folds were
certified. Running the example in the default (unoptimised) profile takes
about 75 s because of the molecular-dynamics-fold and DEM paths; the test
profile is opt-level 3 and the whole suite runs in about 11 s.

The crate is SI throughout (m, kg, s, Pa). Every `Measurement` carries the
material, the test, the `Level` it ran at, the value, a published or exact
reference, the error and a `fold_valid` flag; `report::format_report` prints
one table per material. Nothing is tuned to the number it predicts: the
dynamic EOS is closed from the measured `(c0, s)` and the residual
curvature is reported as an error.

## Levels and folds

| level | model | validity indicator | unfold when |
|---|---|---|---|
| atomistic | LJ pair potential (Ar); measured single-crystal `C_ij` as the bottom rung (Cu, Al) | — | pair potentials inadequate for metals |
| single crystal | anisotropic elasticity / Cauchy–Born | `(ka)²` dispersion error | short wavelength → full phonons |
| polycrystal | Voigt/Reuss/Hill isotropic average | grain size / wavelength `d/λ` | stochastic/geometric scattering |
| continuum | GPR elastic–plastic solid + NASG EOS | `Us(up)` and `σ_HEL` | strength or EOS form fail |

## Materials

| material | density (kg/m³) | atomic model | `C11,C12,C44` (GPa) | source |
|---|---|---|---|---|
| Cu | 8960 | **no** pair potential; measured single-crystal constants are the bottom rung | 168.4, 121.4, 75.4 | Simmons & Wang, *Single Crystal Elastic Constants* (1971) |
| Al | 2700 | as Cu | 107.3, 60.9, 28.3 | Simmons & Wang (1971) |
| Ar | 1416 (LJ, 0 K) | 12–6 Lennard-Jones, self-consistent | 4.42, 2.52, 2.52 (from Cauchy–Born) | Rowley & Nicholson 1972 LJ parameters (`σ = 0.3405 nm`, `ε/k_B = 119.8 K`) |

Polycrystal handbook references: Cu `E = 128 GPa, ν = 0.34`; Al
`E = 70 GPa, ν = 0.35` (Ashby & Jones, *Engineering Materials*; CRC).
Yield stresses (annealed polycrystal, order-of-magnitude): Cu 70 MPa,
Al 35 MPa (Callister, *Materials Science and Engineering*, Table 6.1).
Shock Hugoniots `Us = c0 + s·up`: Cu `c0 = 3.94 km/s, s = 1.49`; Al
`c0 = 5.35 km/s, s = 1.34` (Marsh, *LASL Shock Hugoniot Data*, 1980).

## Static tensile tests

Virtual tensile tests on the single-crystal compliance along
`[100]/[110]/[111]` (reference: the closed-form cubic formula, an
independent code path) and the Hill polycrystal average (reference:
handbook). The single-crystal rows reproduce the closed forms to
`~1e-15` relative; the physically interesting numbers are:

| material | `E[100]` | `E[110]` | `E[111]` | `ν[100]` | `ν[110]` | `ν[111]` |
|---|---|---|---|---|---|---|
| Cu | 66.69 | 130.34 | 191.15 GPa | 0.419 | 0.342 | 0.268 |
| Al | 63.20 | 72.04 | 75.57 GPa | 0.362 | 0.343 | 0.335 |
| Ar | 2.568 | 4.483 | 5.966 GPa | 0.364 | 0.263 | 0.184 |

The `[110]`/`[111]` Poisson ratios are the mean over the two transverse
directions (for `[110]` the two differ, one expanding and one contracting).
The single-crystal `E` values agree with the handbook directional moduli
of Simmons & Wang (1971) to `≲1 %`; the code checks them against the closed
form rather than against a rounded handbook figure. The Hill fold:

| material | `E_Hill` (GPa) | ref | error | `ν_Hill` | ref | error |
|---|---|---|---|---|---|---|
| Cu | 127.35 | 128 | −0.51 % | 0.3452 | 0.34 | +1.5 % |
| Al | 70.38 | 70 | +0.54 % | 0.3464 | 0.35 | −1.0 % |

This is the classic result: the orientation average of the strongly
anisotropic Cu crystal (Zener `A = 3.21`) lands within a percent of the
isotropic handbook values.

## Ultrasonic tests

Christoffel speeds along the three directions (reference: the independent
`materials::continuum` solver, agreeing to `~1e-12` relative) and the Hill
isotropic speeds against handbook values.

| material | `v_L[100]` | `v_L[110]` | `v_L[111]` | `v_T[100]` | `v_T[110]` | `v_T[111]` (m/s) |
|---|---|---|---|---|---|---|
| Cu | 4335 | 4959 | 5150 | 2901 | 2260 | 2134 |
| Al | 6304 | 6452 | 6501 | 3238 | 3084 | 3037 |
| Ar | 1549 | 1806 | 1884 | 1172 | 944 | 895 |

| material | `v_L` Hill | handbook | error | `v_T` Hill | handbook | error | fold |
|---|---|---|---|---|---|---|---|
| Cu | 4727 | 4689 m/s | +0.80 % | 2299 | 2309 m/s | −0.45 % | ok |
| Al | 6418 | 6451 m/s | −0.50 % | 3111 | 3099 m/s | +0.41 % | ok |
| Ar | 1718 | 1716 m/s (mean) | +0.08 % | 962 | 1033 m/s (mean) | −6.9 % | ok (bracket) |

The polycrystal fold is **valid only when the isotropic speeds lie inside
the single-crystal directional bracket** (they do for all three materials)
**and the grains are in the Rayleigh regime**. `scattering_fold(d, f)`
returns the regime and the leading Rayleigh attenuation. At 5 MHz:

| material | `d` = 10 µm | `d` = 10 mm |
|---|---|---|
| Cu | `d/λ = 0.011`, Rayleigh, fold ok | `d/λ = 10.6`, geometric, **unfold** |
| Al | `d/λ = 0.0078`, Rayleigh, fold ok | `d/λ = 7.8`, stochastic, **unfold** |
| Ar | `d/λ = 0.029`, Rayleigh, fold ok | `d/λ = 29`, geometric, **unfold** |

For coarse grains the constant-moduli isotropic average must be unfolded
grain by grain; the crate says so (`fold_valid = false`) instead of
pretending the fold holds. The Rayleigh attenuation keeps only the leading
Stanke–Kino `d³f⁴` scaling and an approximate anisotropy prefactor
(materials.md); it is not a precise attenuation.

## Argon: the atomistic → continuum chain

The only fully self-consistent material. All three folds close:

| fold | measured | reference | error | level |
|---|---|---|---|---|
| Cauchy–Born `C11` | 105.064 | 105.182 `ε/σ³` | −0.11 % | atomistic |
| Cauchy–Born `C12` | 60.119 | 60.187 | −0.11 % | atomistic |
| Cauchy–Born `C44` | 60.131 | 60.187 | −0.09 % | atomistic |
| phonon slope vs Christoffel speed | — | — | max `4.0e-8` | continuum |
| MD `T*` | 2.0448 | 2.0 | +2.2 % | atomistic |
| MD `P*` | 1.0926 | 1.0775 (JZG) | +1.4 % | atomistic |
| MD `U_res*` | −3.1418 | −3.1449 (JZG) | +0.10 % | atomistic |

The Cauchy–Born constants use `rc = 12σ, 20σ` and the `rc⁻³`
extrapolation of `materials.md`; the reference is the converged
infinite-lattice sum (Jones 1924; Lennard-Jones & Ingham 1925). The
central-force Cauchy relation `C12 = C44` holds to `1.1e-3`. The
phonon/Christoffel comparison uses the same raw-cutoff configuration and
`q = 1e-3`, so the `4.0e-8` residual is the `(ka)²` dispersion error, not a
model discrepancy. The MD state point is 256 atoms, tail-corrected LJ, NVT
with a Bussi thermostat in reduced units; the reference is the
Johnson–Zollweg–Gubbins (1993) MBWR equation of state. The `+2 %`
temperature and `+1 %` pressure residuals are the statistical scatter of a
short run, not a tolerance.

## Dynamic tests: 1D symmetric plate impact

The flyer and target are the same material; in the target frame the flyer
moves at `V = 200 m/s` and the target is at rest, so symmetry fixes the
particle velocity behind the plastic shock to exactly `V/2`. The solid is
the GPR elastic–plastic model with the **polycrystal** shear speed
`c_s = √(G_Hill/ρ)`, the polycrystal bulk speed `c0 = √(K/ρ)` and the
handbook yield stress; the EOS is a NASG stiffened gas with
`γ = 4s − 1` (so `c0` and the weak-shock slope are matched) and
`p∞ = ρc0²/γ`. Grid `n = 600`, domain `±0.4 m`, CFL 0.4, Strang splitting,
`n = 400` hardening.

| material | quantity | measured | reference | error | fold |
|---|---|---|---|---|---|
| Cu | particle velocity | 100.01 | 100.00 m/s | +0.015 % | ok |
| Cu | `σ_HEL` | 143.1 | 148.0 MPa | −3.3 % | ok |
| Cu | `u_HEL` | 3.383 | 3.495 m/s | −3.2 % | ok |
| Cu | precursor speed | 4765 | 4727 m/s | +0.82 % | ok |
| Cu | shock speed | 4065 | 4089 m/s | −0.58 % | ok |
| Cu | NASG–RH shock speed | 4063 | 4089 m/s | −0.64 % | ok |
| Al | particle velocity | 100.00 | 100.00 m/s | +0.008 % | ok |
| Al | `σ_HEL` | 72.2 | 74.5 MPa | −3.1 % | ok |
| Al | `u_HEL` | 4.162 | 4.297 m/s | −3.2 % | ok |
| Al | precursor speed | 6461 | 6418 m/s | +0.66 % | ok |
| Al | shock speed | 5453 | 5484 m/s | −0.57 % | ok |
| Al | NASG–RH shock speed | 5454 | 5484 m/s | −0.55 % | ok |

Here `σ_HEL = Y(1−ν)/(1−2ν)` (the uniaxial-strain von Mises yield stress, a
few percent below the ideal-plasticity value because the GPR power law
yields around `0.93 σ_Y`). The particle velocity is exact to the grid. The
elastic precursor is `0.6–0.8 %` above `c_L` — finite-amplitude steepening
plus a little grid diffusion at `n = 600`; the residual falls with
resolution. The plastic shock is within `0.6 %` of the Marsh Hugoniot.

### The EOS-form error, reported separately

`rh_us` solves the stiffened-gas Rankine–Hugoniot exactly. With `γ` set
from the weak-shock slope, the closure reproduces the linear Marsh fit well
at weak shocks and drifts at strong ones:

| material | `γ` | `Us(100)` NASG | Marsh | error | `Us(500)` NASG | Marsh | error |
|---|---|---|---|---|---|---|---|
| Cu | 4.960 | 4092 | 4089 m/s | +0.07 % | 4755 | 4685 m/s | +1.49 % |
| Al | 4.360 | 5486 | 5484 m/s | +0.03 % | 6062 | 6020 m/s | +0.69 % |

A NASG/stiffened-gas EOS gives a **linear `Us–up` only asymptotically**;
the curvature is `+1.5 %` at `up = 500 m/s` for Cu and `+0.7 %` for Al. The
solver's own shock speed follows its closure to `~1e-3` (see the table), so
the `~0.6 %` gap to the Marsh line is the EOS form, not the scheme.

## Honest misses and limits

- **Pair potentials for metals.** Cu and Al have no atomistic rung here.
  A central pair potential gives the Cauchy relation `C12 = C44` and a
  negative stacking-fault energy, neither of which is true for a metal;
  `analysis/defects.md` documents this failure. The measured single-crystal
  constants are used as the bottom rung and the report says so.
- **Solid argon is not a classical LJ crystal.** The LJ Cauchy–Born modulus
  is a 0 K classical result; quantum solid argon has a larger lattice
  constant and a lower modulus, so there is no handbook comparison at the
  single-crystal or polycrystal rungs. Only the internal folds (Cauchy–Born
  ↔ phonons, MD ↔ JZG) are checked.
- **The polycrystal fold is order-of-magnitude.** Voigt/Reuss/Hill is an
  estimate; Hashin–Shtrikman is tighter but still not a rigorous bound for
  a finite body. The Rayleigh attenuation uses a leading anisotropy
  prefactor, not the full Stanke–Kino angular integrals.
- **The dynamic model is 1D and uses idealised strength.** The GPR
  closure uses a single yield stress and a power-law hardening exponent,
  with no rate dependence, no dislocation density and a homogeneous
  polycrystal. The `σ_HEL` residual (`−3 %`) is the power law, not a fit.
- **The NASG Hugoniot is curved.** Matching `c0` and the weak-shock slope
  leaves a `0.7–1.5 %` curvature error at `up = 500 m/s`. The report keeps
  the closure error (`EOS shock speed`) separate from the solver error
  (`shock speed`).
- **Grid resolution.** The plate-impact precursor speed and `σ_HEL` carry
  a `~1–3 %` grid/nonlinearity error at `n = 600`; the particle velocity
  and shock speed are essentially grid-converged. Tests use small 1D grids
  so the whole suite runs in about 11 s.

## Constants and sources

- Single-crystal elastic constants and densities: Simmons & Wang,
  *Single Crystal Elastic Constants and Calculated Aggregate Properties:
  A Handbook* (MIT Press, 1971); densities from the CRC Handbook.
- Polycrystal `E`, `ν` and ultrasonic speeds: Ashby & Jones,
  *Engineering Materials*; CRC Handbook.
- Yield stresses: Callister, *Materials Science and Engineering*,
  Table 6.1 (annealed polycrystal, order of magnitude).
- Shock Hugoniot `Us = c0 + s·up`: Marsh, *LASL Shock Hugoniot Data*
  (1980).
- Argon LJ parameters: Rowley & Nicholson (1972), as carried by
  `md::argon`; lattice sums and converged Cauchy–Born constants:
  `analysis/materials.md` (Jones 1924; Lennard-Jones & Ingham 1925).
- Lennard-Jones reference equation of state: Johnson, Zollweg & Gubbins,
  *Mol. Phys.* **78**, 591 (1993).
- Cauchy–Born / Christoffel / harmonic phonons, Voigt–Reuss–Hill,
  Hashin–Shtrikman and the grain-size/wavelength indicator: the
  `materials` crate and `analysis/materials.md`.
- GPR elastic–plastic solid: Peshkov & Romenski (2016); Dumbser et al.
  (2016), implemented in `compressible/gpr.rs`; NASG EOS: Le Métayer &
  Saurel (2016), `compressible/eos.rs`.

# Extended battery: the other emergent rungs

The battery now runs every material named in the tables below through the
standard tests of the mechanical, thermal, viscoelastic, magnetic, optical
and interfacial/granular rungs, reusing the `Measurement`/`Level` structure
and the same `fold_valid` convention. The report example prints one table
per material and then a cross-material summary of certified and unfolded
folds; the whole battery is 100 measurements, **84 certified and 16
unfolded**. The heavy paths (the quasi-harmonic argon spectrum, the small
DEM shear cell) are cached with `OnceLock`, so repeated calls are free and
the test suite stays at about 11 s.

| material | model | reference |
|---|---|---|
| Si | diamond, `a = 0.5431 nm`; phonons from the Morelli–Kazan set | Glassbrenner & Slack 1964 |
| polystyrene | Rouse/Prony melt, `T_g = 373 K`, `G_N⁰ = 2.0×10⁵ Pa` | Ferry 1980; Williams–Landel–Ferry 1955 |
| Fe (bcc) | `A = 2.1e-11 J/m`, `M_s = 1.71e6 A/m`, `K₁ = 4.8e4 J/m³` | Coey 2010 |
| Ni₈₀Fe₂₀ | `A = 1.3e-11 J/m`, `M_s = 8.0e5 A/m`, `K ≈ 0` | Coey 2010 |
| Au | Rakić Lorentz–Drude fit | Rakić et al. 1998 |
| 3C-SiC | `ω_TO = 796.2`, `ω_LO = 972.7 cm⁻¹`, `ε_∞ = 6.52` | Palik 1985; Spitzer et al. 1959 |
| water | `γ = 72.8 mN/m`, `ρ = 998.21 kg/m³` | Vargaftik et al. 1983; CRC |
| glass beads | `μ_s = tan 20.9°`, `μ_2 = tan 32.8°`, `I₀ = 0.279` | Jop et al. 2006 |

## Mechanical: the fracture and yield folds

The fracture test is the Griffith energy balance
`K_IC = √(2 γ_s E/(1−ν²))`. For LJ argon both inputs are atomistic: `γ_s`
is the broken-bond lattice sum (`materials::defects`) and `E, ν` the
Cauchy–Born Hill moduli. The reference is the anisotropic Stroh
Irwin-matrix toughness from the *same* atomistic inputs, so the residual is
the Zener anisotropy of LJ fcc. For the metals the handbook surface energy
and measured moduli are used and compared with the measured `K_IC`: the
fold fails by two orders of magnitude because LEFM describes cleavage, not
the dislocation-mediated tearing of a ductile metal.

| material | `K_IC` (Griffith) | reference | error | fold |
|---|---|---|---|---|
| Cu | 0.721 MPa√m | 100 MPa√m (ductile) | −99.3 % | **no**; Rice ratio `γ_us/γ_s = 0.089 < 1/4` predicts emission |
| Al | 0.438 MPa√m | 35 MPa√m (ductile) | −98.7 % | **no**; `γ_us/γ_s = 0.133 < 1/4` |
| Ar | 0.02030 MPa√m | 0.01917 MPa√m (Stroh) | +5.9 % | ok (anisotropy `A = 2.68`) |
| Ar, 2D sheet | 10.84322 `ε/σ^{3/2}` | 10.84320 (exact) | +2×10⁻⁶ | ok; the finite atomistic crack advances at `K_c/K_G ≈ 1.48` |

The fold is simply **invalid** for Cu and Al because LEFM describes cleavage,
not dislocation-mediated tearing.

The yield test uses the Peierls–Nabarro fold where a pair potential supplies
the γ-surface. For LJ argon the classical Peierls stress reproduces the
tabulated `3.296 ε/σ³ = 0.081 G` to −0.1 %, but the core width implied by
the measured atomistic `γ_us = 0.4226 ε/σ²` is `3.96 σ` against the
classical `0.61 σ`, so the fold is flagged invalid for a quantitative yield
stress. The ideal fault-plane strength `π γ_us/b` matches `1.22 ε/σ³` to
−0.2 %. For Cu and Al only the handbook polycrystal yield exists
(70 and 35 MPa); a pair potential cannot supply the CRSS, and the report
says so.

## Thermal: quasi-harmonic argon and phonon conductivity of silicon

The argon rung minimises the quantum quasi-harmonic free energy over volume
on a `4³` Monkhorst–Pack grid (`materials::quasiharmonic`). The thermal
expansion is compared with the experimental lattice constants of Peterson,
Batchelder & Simmons (1966); the heat capacity with the Debye model at the
calorimetric `θ_D = 93.5 K`.

| `T` (K) | `a_QHA` (Å) | `a_exp` (Å) | error | `C_V/k_B` | Debye | error |
|---|---|---|---|---|---|---|
| 4.2 | 5.3096 | 5.311 | −0.03 % | 0.079 | 0.021 | +275 % |
| 10 | 5.3117 | 5.312 | −0.01 % | 0.617 | 0.274 | +125 % |
| 20 | 5.3225 | 5.317 | +0.10 % | 1.737 | 1.227 | +42 % |
| 40 | 5.3729 | 5.342 | +0.58 % | 2.630 | 2.317 | +13.5 % |

The expansion is accurate at the low-temperature anchors (the quantum
zero-point expansion is captured) and degrades on heating as the harmonic
spinodal approaches. The heat-capacity comparison is an **honest miss** at
low `T`: the full fcc spectrum is not Debye-like, and the finite `4³` grid
departs from the Debye continuum, so the Debye reference underestimates
`C_V` by 42 % at 20 K. The well-resolved 40 K point is within 13.5 % and is
the only heat-capacity row certified. This is reported, not hidden.

For silicon the `heat` crate gives the Callaway/Holland kinetic
conductivity with the fitted Morelli–Kazan parameters and a Casimir
boundary rate `τ_B⁻¹ = v/L`. The phonon Knudsen number `Kn = Λ/L` flags the
fold.

| sample | `κ` (W/(m K)) | reference | `Kn` | regime | Fourier valid? |
|---|---|---|---|---|---|
| bulk (14 mm) | 126.96 | 148 (Glassbrenner–Slack) | 1.8×10⁻⁵ | Fourier | yes; the −14 % is the fitted parameter set |
| 100 nm film | 27.05 | 127.57 (bulk) | 2.46 | transitional | no |
| 10 nm film | 4.36 | 127.57 (bulk) | 24.6 | ballistic | no |

The size effect lowers `κ` by factors of 4.7 and 29; the diffusion closure
is certified only for the bulk sample.

## Viscoelastic: DMA and WLF on a polystyrene melt

The polymer crate's virtual rheometer integrates the Maxwell ODEs under
`γ = γ₀ sin ωt` and projects the steady stress onto `sin`/`cos`; the
reference is the analytic Prony `G*(ω)`.

| `ω` (rad/s) | `G'` ODE (Pa) | `G'` analytic (Pa) | `G''` ODE (Pa) | `G''` analytic (Pa) |
|---|---|---|---|---|
| 1 | 428.16 | 428.16 | 5814.9 | 5814.1 |
| 10 | 23075.9 | 23060.3 | 37873.1 | 37889.7 |
| 100 | 114341 | 112939 | 69071 | 69414 |

The largest relative error is 1.2 % at `ω = 100 rad/s` (`De = 10`), where
the slowest mode is only marginally resolved by the 40-period window. The
plateau modulus is the sum of the Prony weights, `G_N⁰ = 2.0×10⁵ Pa`, exact
against Ferry (1980). The WLF shift `log₁₀ a_T = −8.58268` at
`T_g + 50 K` matches the universal `C₁ = 17.44, C₂ = 51.6 K` to 3×10⁻⁶,
and the time–temperature superposition `G*(ω;T) = G*(a_Tω;T_ref)` holds to
machine precision by construction. The Deborah number `De = λ₁ω` certifies
the linear-viscoelastic closure only in the window `0.1 < De < 10`; at
`De = 10` the fold is at its edge.

## Magnetic: domain walls and Stoner–Wohlfarth coercivity

| material | quantity | value | reference | error | fold |
|---|---|---|---|---|---|
| Fe | wall width `δ` | 20.92 nm | 20–40 nm (Coey; Cullity) | +4.6 % | ok |
| Fe | `H_c` macrospin | 22337 A/m | 80 A/m (soft iron) | +2.8×10² | **no** (`D ≫ D_c = 9.8 nm`) |
| Fe | `h_sw(45°)/H_K` | 0.50000 | 0.5 (astroid) | −3×10⁻⁶ | ok |
| Ni₈₀Fe₂₀ | exchange length | 5.69 nm | ≈5 nm | +13.7 % | ok |
| Ni₈₀Fe₂₀ | `H_c` macrospin | 0 A/m | 1 A/m | −100 % | **no** (`K ≈ 0`, shape-dominated) |

The wall width and exchange length are the micromagnetics rung; the LLG
macrospin reproduces the exact Stoner–Wohlfarth astroid minimum at 45°.
The macrospin coercivity is flagged invalid for both bulk magnets because
they are far larger than the single-domain diameter and reverse by
domain-wall motion and pinning, not coherent rotation — the two differ by
orders of magnitude, which is the point of the validity indicator.

## Optical: gold reflectance and the SiC Reststrahlen band

| material | quantity | value | reference | error | fold |
|---|---|---|---|---|---|
| Au | `R(633 nm)` | 0.8927 | 0.933 (Johnson–Christy/Palik) | −4.3 % | ok |
| SiC | `R` mid-band | 0.9635 | 1 (lossless mirror) | −3.7 % | ok |
| SiC | `ω_TO` | 796.2 cm⁻¹ | 790 | +0.8 % | ok |
| SiC | `ω_LO` | 972.7 cm⁻¹ | 970 | +0.3 % | ok |
| SiC | `R(0.5 ω_TO)`, 2 µm film | 0.57045 | 0.5704 (TMM) | +9×10⁻⁵ | ok |

The gold reflectance uses the causal Rakić Lorentz–Drude fit; the −4.3 %
miss is the interband loss that a pure Drude form would omit. The SiC
Reststrahlen band is a perfect mirror where `Re ε < 0`, with the TO/LO
edges reproduced to a few wavenumbers; the in-band damped reflectance is
0.9635 rather than 1 because of the finite linewidth.

## Interfaces and granular: textured wetting and `μ(I)`

| material | quantity | value | reference | error | fold |
|---|---|---|---|---|---|
| water (textured) | Wenzel `θ*` (fine model) | 120.85° | 123.18° (fold) | −1.9 % | ok (`a/λ = 17.2`) |
| water (textured) | Cassie `θ*` (fine model) | 137.11° | 137.46° (fold) | −0.3 % | ok |
| water (textured) | sag/height | 0.176 | 0 (no impalement) | — | ok |
| glass beads | `μ(γ̇=0.005)` | 0.330 | 0.401 (Jop `μ(I)`) | −17.8 % | ok (dense flow) |
| glass beads | `μ(γ̇=0.02)` | 0.425 | 0.446 (Jop `μ(I)`) | −4.6 % | ok |
| glass beads | repose angle | 18.25° | 20.9° (Jop `μ_s`) | −12.7 % | **no** (dense, not quasi-static) |

The wetting fine model is a direct free-energy minimisation on an explicit
2D square wave (`r = 1.60`, `φ_s = 0.40`, `θ_Y = 110°`); at `a/λ = 17.2`
it lands within 2.3° of the Wenzel fold, and the residual is contact-line
pinning. The granular `(I, μ)` points come from a small (108-sphere)
constant-pressure DEM shear cell; they reproduce the Jop glass-bead `μ(I)`
law to 5–18 %, with the largest miss at the lowest inertial number where the
DEM cross-over `I₀` differs from the glass-bead value. The angle of repose
is flagged invalid because the small cell is in the dense-flow regime, not
quasi-static — the same honest conclusion as `analysis/granular.md`.

## Cross-material fold summary

The report prints, for every material, the number of tests and the folds
that must be unfolded. The 16 unfolded folds are the physically expected
ones: the ductile-metal Griffith toughness (Cu, Al), the coarse-grain
ultrasonic scattering fold (all three elastic materials), the argon Peierls
core width, the three low-temperature argon heat-capacity points on the
coarse QHA grid, the ballistic silicon films (100 nm and 10 nm), the
macrospin coercivity of both bulk magnets, the granular angle of repose,
and the collisional `I` indicator. Everything else — single-crystal
elasticity, ultrasonic speeds, the atomistic Cauchy–Born/Christoffel folds,
the polycrystal Hill bracket, the argon MD state point, the polystyrene DMA
and WLF, the iron wall width and astroid, the gold and SiC reflectances, and
the wetting Wenzel/Cassie folds — is certified.

## Additional constants and sources

- Solid-argon `a(T)`: Peterson, Batchelder & Simmons, *Phys. Rev.* **150**,
  703 (1966) (rounded anchors); calorimetric `θ_D ≈ 93.5 K` (standard
  rare-gas table). See [quasiharmonic.md](quasiharmonic.md).
- Silicon conductivity `κ(300 K) = 148 W/(m K)`: Glassbrenner & Slack,
  *Phys. Rev.* **134**, A1058 (1964); parameters from Morelli, Heremans &
  Slack (2002) via Kazan et al. (2010). See [heat.md](heat.md).
- Polystyrene `T_g`, `G_N⁰`: Ferry, *Viscoelastic Properties of Polymers*,
  3rd ed. (1980); WLF shift: Williams, Landel & Ferry, *J. Am. Chem. Soc.*
  **77**, 3701 (1955). See [polymer.md](polymer.md).
- Magnetic constants: Coey, *Magnetism and Magnetic Materials* (2010);
  Cullity & Graham, *Introduction to Magnetic Materials* (2nd ed.);
  wall width `δ = √(A/K)`; Stoner–Wohlfarth astroid (1948). See
  [magnetism.md](magnetism.md).
- Gold optical constants: Rakić et al., *Appl. Opt.* **37**, 5271 (1998);
  measured reflectance: Johnson & Christy (1972) / Palik, *Handbook of
  Optical Constants of Solids*. SiC phonons: Spitzer, Kleinman & Walsh,
  *Phys. Rev.* **113**, 127 (1959); Palik (1985). See
  [photonics.md](photonics.md) and [polariton.md](polariton.md).
- Water surface tension and density: Vargaftik, Volkov & Voljak, *J. Phys.
  Chem. Ref. Data* **12**, 817 (1983); IAPWS-95/CRC. Texture geometry and
  Wenzel/Cassie balances: Wenzel (1936); Cassie & Baxter (1944). See
  [wetting.md](wetting.md).
- Glass-bead `μ(I)` constants: Jop, Forterre & Pouliquen, *Nature* **441**,
  727 (2006). See [granular.md](granular.md).
- Griffith and Rice criteria: Griffith (1921); Rice, *J. Mech. Phys.
  Solids* **40**, 239 (1992). Peierls–Nabarro: Peierls (1940); Nabarro
  (1947). Metal `K_IC` and yields: Callister, *Materials Science and
  Engineering*; the handbook values are typical annealed polycrystal
toughnesses and the surface energies are estimates in the ranges 1.3–2.2
J/m² (Cu) and 1.1–1.2 J/m² (Al), flagged in the code. See
[fracture.md](fracture.md) and [defects.md](defects.md).
