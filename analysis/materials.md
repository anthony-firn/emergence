# Materials ladder: atoms → continuum elasticity

Implemented in `fluid/crates/materials`. This is the first rung of the
materials ladder in [emergence_ladder.md](emergence_ladder.md): the fine
model is an atomistic crystal with a central pair potential, the coarse
model is continuum elasticity, and the fold between them is the
Cauchy–Born rule. The same structure is then repeated at the next rungs:
real single crystals (anisotropic continuum) fold to isotropic
polycrystals, with the grain-size/wavelength ratio as the indicator.

- **coarse model** — a homogeneous strain field, the 6×6 Voigt stiffness
  `C_ij`, and the acoustic sound speeds from the Christoffel equation;
- **validity indicator** — the relative dispersion error
  `(ω_lattice − v k)/ω_lattice` as a function of `k·a`, zero at long
  wavelength and growing like `(k a)²`; at the polycrystal rung, the
  grain-size/wavelength ratio `d/λ`;
- **fold / unfold** — apply the Cauchy–Born deformation to the lattice and
  relax the basis (non-Bravais only); when the indicator is not small,
  unfold to the full harmonic phonon dispersion. A single crystal folds to
  an isotropic polycrystal by orientation averaging (Voigt/Reuss/Hill and
  Hashin–Shtrikman); when `d ~ λ` or `d ≫ λ` the isotropic fold fails and
  the crystal must be unfolded grain by grain.

## The rung

| module | content |
|---|---|
| `lattice.rs` | primitive vectors + fractional basis for sc, bcc, fcc, hcp (ideal and given `c/a`), diamond, and a custom cell; neighbour shells and the spherical-cutoff bond list (periodic images included, no minimum-image approximation) |
| `potential.rs` | Lennard-Jones and Morse with exact `φ′`, `φ″`; a harmonic spring for exact tests; true lattice sums, a homogeneous LJ tail correction, and equilibrium by minimising over the lattice constant |
| `lattice/linalg` | a cyclic Jacobi eigensolver and a pivoted dense solve, the only numerical dependencies |
| `elastic.rs` | Cauchy–Born energy density `W(F,u)`, analytic Hessians, internal-basis relaxation, the relaxed (Schur-complement) and frozen stiffnesses, and a finite-difference check |
| `phonon.rs` | the harmonic dynamical matrix `D(k)` for any lattice/potential, phonon frequencies, the exact 1D-chain reference, and an fcc high-symmetry path Γ–X–W–L–Γ |
| `continuum.rs` | the Christoffel matrix `Γ_ik = C_ijkl n_j n_l / ρ`, its three sound speeds, and the dispersion-error indicator |
| `anisotropic.rs` | real single-crystal elastic constants and densities (cubic, hexagonal, general triclinic), Bond-matrix orientation rotation, the compliance `S = C⁻¹`, directional Young's modulus and the Zener ratio, and the Christoffel speeds/polarisations |
| `polycrystal.rs` | Voigt/Reuss/Hill orientation averages, Hashin–Shtrikman shear bounds for cubic polycrystals, isotropic moduli and sound speeds, the GPR solid-model adapter, and the grain-size/wavelength validity indicator with the Rayleigh attenuation |

The pair energy per atom of a crystal with `n` basis atoms is

\[
  E=\frac{1}{2n}\sum_i\sum_{R}\varphi(|R+\tau_j-\tau_i|),
\]

the `1/2` because every bond is shared by two atoms; the sum is a true
periodic lattice sum. For Lennard-Jones,

\[
  E=2\varepsilon\left[A_{12}\left(\frac{\sigma}{r_{nn}}\right)^{12}
      -A_6\left(\frac{\sigma}{r_{nn}}\right)^{6}\right],
\]

so the equilibrium and cohesive energy are fixed by the fcc sums
`A12 = 12.13188`, `A6 = 14.45392` (classic fcc lattice sums; Jones 1924,
Lennard-Jones & Ingham 1925). The r⁻⁶ sum converges slowly; the
homogeneous tail `Σ_{r>rc}(r_nn/r)^6 ≈ 4πρ r_nn^6/(3 rc³)` (standard
mean-field molecular-simulation tail) is applied so a modest cutoff
already reproduces it.

For a deformation gradient `F` and internal displacements `u`, the energy
density per reference volume is
`W(F,u) = (1/2V0) Σ φ(|F R + u_j − u_i|)`. The physical (relaxed)
stiffness is the Schur complement

\[
  C = W_{FF} - W_{Fu}\,K^{-1}\,W_{uF},\qquad K=W_{uu},
\]

with the translation gauge fixed by holding atom 0. For a Bravais lattice
`K` is empty and `C = W_{FF}`. The dynamical matrix enforces the acoustic
sum rule by construction: each bond adds `M` on site and `−M e^{ik·R}`
across the bond, with
`M_ab = φ″ r̂_a r̂_b + (φ′/r)(δ_ab − r̂_a r̂_b)`.

## Benchmarks and measured errors

LJ runs use `σ = ε = 1` (so `C` is in `ε/σ³`, `ρ` in `σ⁻³`, and with unit
masses the sound speeds are in `(ε/m_σ)^{1/2}`). All numbers below are
produced by the checked-in tests in `tests/atomistic.rs`; the full suite
runs in about 1 s.

### Neighbour shells (exact counts)

| structure | first shells | expected | measured |
|---|---|---|---|
| fcc | 12, 6, 24, 12 | 12, 6, 24, 12 | exact |
| bcc | 8, 6, 12 | 8, 6, 12 | exact |
| diamond | 4, 12 | 4, 12 | exact |
| ideal hcp | 12 | 12 | exact, `r = a` to 1e-12 |

### Static lattice sums and cohesion

| quantity | reference | computed | error |
|---|---|---|---|
| fcc `A12` (`rmax=3`) | 12.13188 | 12.131876 | 4e-6 |
| fcc `A6` tail-corrected `rmax=4` | 14.45392 | 14.453605 | 3.1e-4 |
| fcc `A6` tail-corrected `rmax=6` | 14.45392 | 14.453839 | 8.1e-5 |
| fcc `A6` tail-corrected `rmax=8` | 14.45392 | 14.453927 | 7.4e-6 |
| fcc `A6` tail-corrected `rmax=10` | 14.45392 | 14.453921 | 1.0e-6 |
| fcc LJ `r0` | 1.09017297 σ | 1.09017298 σ | 9e-9 (8e-9 rel) |
| fcc LJ cohesive energy | −8.61019905 ε | −8.61020012 ε | 1.1e-6 abs (1.2e-7 rel) |

Phase ordering at the same (tail-corrected) cutoff: `E_hcp = −8.611097 <
E_fcc = −8.610200 < E_bcc = −8.237300`. The hcp result is lower than fcc
by `9.0e-4 ε`, i.e. **0.0104 %**, consistent with the known ~0.01 %.

### Cauchy–Born elasticity (LJ fcc at zero pressure)

At the zero-pressure (raw-cutoff) fcc lattice constant the analytic and
finite-difference stiffnesses agree; these are the constants of the
*truncated* potential and are deliberately not converged in the cutoff:

| constant | analytic | finite difference | error |
|---|---|---|---|
| `C11` | 104.92733 | 104.92735 | 2e-7 rel |
| `C12` | 60.01864 | 60.01866 | 3e-7 rel |
| `C44` | 60.01864 | 60.01866 | 3e-7 rel |

Evaluated at the exact `r0` and extrapolated to infinite cutoff (the
leading truncation tail of the second derivative is `∝ rc⁻³`; two
cutoffs, `rc = 12σ` and `20σ`, are enough), the constants converge to the
independent infinite-sum reference:

| constant | reference | extrapolated | error |
|---|---|---|---|
| `C11` | 105.19 ± 0.02 | 105.182 | 0.008 |
| `C12` | 60.19 ± 0.02 | 60.187 | 0.003 |
| `C44` | 60.19 ± 0.02 | 60.187 | 0.003 |

The reference values are the converged (infinite-sum) static values;
`C11` pushes to `105.182` by `rc = 120σ` (the raw `8σ` value is
`105.331`, `C12 = C44 = 60.211`). The central-force Cauchy relation holds
to `1e-5` relative after extrapolation.

Cubic symmetry holds to <1e-9 and the central-force Cauchy relation
`C12 = C44` holds to <1e-6 relative (at zero pressure, a Bravais lattice,
no internal relaxation). Internal relaxation is essential for the
non-Bravais structures:

| structure / potential | frozen | relaxed | finite difference |
|---|---|---|---|
| diamond, NN spring, `C44` | 0.288675 | 0.000000 | 0.000000 |
| hcp, NN spring, `C11` | 3.535534 | 3.417683 | 3.417683 |
| hcp, NN spring, `C12` | 1.178511 | 1.296362 | 1.296362 |
| hcp, LJ, `C11` | 142.48845 | 133.54035 | 133.54050 |

Diamond with only nearest-neighbour central forces is isostatic
(`z = 4`, `d = 3`); the affine `C44` is finite but the sublattice
relaxation is a floppy mode and removes it completely. The relaxed
values agree with the independent finite difference of the *relaxed*
energy to <1e-4 relative.

### Harmonic lattice dynamics

| benchmark | reference | error |
|---|---|---|
| 1D monatomic chain `ω = 2√(K/m)|sin(ka/2)|` | exact | <1e-15 (machine precision) |
| acoustic sum rule, `ω(Γ)` | 0 | all `|ω²| < 1e-10` |
| stability at equilibrium (fcc LJ, several k) | `ω² ≥ 0` | all `ω² > −1e-9`, branches rise from Γ |

### The fold: Christoffel sound speeds vs phonons

At the fcc zero-stress equilibrium, with `k = q n` and `q = 1e-3`:

| direction | mode | `v_lattice` | `v_Christoffel` | relative error |
|---|---|---|---|---|
| [100] | T (×2) | 7.41724418 | 7.41724437 | 2.6e-8 |
| [100] | L | 9.80717273 | 9.80717265 | 8.1e-9 |
| [110] | T | 4.53680056 | 4.53680063 | 1.7e-8 |
| [110] | T | 7.41724430 | 7.41724437 | 9.5e-9 |
| [110] | L | 11.42863001 | 11.42863025 | 2.1e-8 |
| [111] | T (×2) | 5.66217350 | 5.66217373 | 4.0e-8 |
| [111] | L | 11.92020292 | 11.92020305 | 1.1e-8 |

The largest relative error over all directions and modes at `q = 1e-3`
is `4.0e-8`. Doubling `q` multiplies the error by 4.000, confirming the
`(k a)²` approach to the continuum limit. The comparison must use the
same (zero-stress) reference for the phonons and for `C`: at a
configuration with residual stress the long-wavelength dynamical matrix
gives the Huang tensor, which differs from the energy second derivative
by stress terms.

### Validity indicator

The relative dispersion error of the longitudinal [100] branch:

| `q a` | `(ω_lat − v k)/ω_lat` |
|---|---|
| 5e-3 | 2.031e-7 |
| 1e-2 | 8.124e-7 |

The ratio is 4.000 for a factor-2 change in `q`, i.e. the elasticity
fold is valid to `O((ka)²)`. Operationally the indicator is evaluated at
the local wavenumber; the continuum rung is used while it stays below the
required accuracy, and the full `D(k)` is unfolded where it is not.

## The anisotropic continuum: real single crystals

The Cauchy–Born rung fixes the *form* of the stiffness but not its
numbers for a real material. `anisotropic.rs` takes measured
room-temperature single-crystal constants and propagates them. Cubic
crystals need `(C11, C12, C44)`, hexagonal ones `(C11, C12, C13, C33,
C44)` with `C66 = (C11−C12)/2`, and the general triclinic case a full
symmetric 6×6.

| material | class | ρ (kg/m³) | source |
|---|---|---|---|
| Si, Ge, diamond, Cu, Al, Ni, α-Fe, W, NaCl | cubic | 2329, 5323, 3515, 8960, 2700, 8900, 7874, 19250, 2165 | Si: Hopcroft et al. 2010; rest: Simmons & Wang 1971; densities CRC |
| Mg, Ti, Zn | hexagonal | 1738, 4507, 7140 | Simmons & Wang 1971 (Zn rounded, flagged `estimate`) |

The Bond matrix rotates the stiffness to any crystal orientation,
`C' = M C Mᵀ`; it is checked against the direct fourth-order contraction
`C'_{ijkl} = R_{im}R_{jn}R_{ko}R_{lp}C_{mnop}` to `1e-12` relative. The
compliance is `S = C⁻¹`; `youngs_modulus(C, n) = 1/S'₁₁`. Directional
properties are invariants, so rotating the crystal and the direction
together leaves the speeds and `E` unchanged (checked to `1e-9`).

### Directional benchmarks

Silicon (Hopcroft, Nix & Kenny, *J. Microelectromech. Syst.* **19**, 229
(2010)):

| quantity | reference | computed | error |
|---|---|---|---|
| `v_L` [100] | 8433 m/s | 8434.8 m/s | 2.1e-4 |
| `v_T` [100] (×2) | 5843 m/s | 5846.2 m/s | 5.5e-4 |
| `v_L` [110] | 9134 m/s | 9136.2 m/s | 2.4e-4 |
| `E` [100] | 130 GPa | 130.13 GPa | 0.10 % |
| `E` [110] | 169 GPa | 169.10 GPa | 0.06 % |
| `E` [111] | 188 GPa | 187.85 GPa | 0.08 % |

Anisotropy `A = 2C44/(C11−C12)`: tungsten `A = 1.011` (nearly
isotropic), copper `A = 3.209` (strongly anisotropic), matching the
literature values. The Christoffel eigenproblem returns the polarisations
as well as the speeds; they are orthonormal and satisfy `Γp = v²p` to
`1e-6` in every material and direction tested.

## The fold to isotropic polycrystals

A random-orientation polycrystal is statistically isotropic at scales
large compared with a grain. Orientation averaging of the stiffness
(Voigt, uniform strain) and of the compliance (Reuss, uniform stress)
brackets the effective response; the Hill average is their mean. For a
cubic crystal the bulk modulus is exactly orientation-independent,
`K = (C11+2C12)/3`, so only the shear modulus is nontrivial.

| material | `G_V` | `G_R` | `G_Hill` | measured `G` | `E_Hill` | measured `E` |
|---|---|---|---|---|---|---|
| Cu | 54.64 | 40.03 | 47.34 | ≈48 GPa | 127.4 | ≈128 GPa |
| Al | 26.26 | 26.01 | 26.14 | ≈26 GPa | 70.4 | ≈70 GPa |
| W | 160.08 | 160.08 | 160.08 | — | 409.8 | — |
| Si | 68.12 | 64.95 | 66.54 | — | 162.7 | — |

The measured polycrystalline values are the usual room-temperature
quasi-static handbook numbers (Ashby & Jones, *Engineering Materials*;
CRC Handbook). The Hill average is within a few percent in both cases.

### Hashin–Shtrikman bounds

The Hashin–Shtrikman shear bounds use the standard two-phase formula with
the cubic shear modes as the two phases: `C' = (C11−C12)/2` (weight 2/5)
and `C44` (weight 3/5), with the bulk modulus `K` common to both (Watt,
Davies & O'Connell 1976; Peselnick & Meister 1975). They are tighter than
Voigt–Reuss and contain the Hill average to 0.1 %:

| material | `G_HS⁻` | `G_HS⁺` | width |
|---|---|---|---|
| Cu | 45.96 | 49.45 | 3.49 GPa |
| Al | 26.142 | 26.152 | 0.010 GPa |
| Ni | 84.95 | 88.10 | 3.15 GPa |
| W | 160.078 | 160.078 | 0.011 GPa |

This construction uses the orientation average of the cubic shear
projector; it is the leading bound and is not the full
Peselnick–Meister angular integral (the Hill average can fall a hundredth
of a percent outside it, as for Al).

### Dropping a polycrystal into the GPR solid solver

`polycrystal::isotropic_parameters_stiffness` returns the density, the
shear sound speed `c_s = √(G/ρ)`, the bulk sound speed `c_b = √(K/ρ)` and
the longitudinal speed `√((K+4G/3)/ρ)`. To drive the existing GPR model
in `fluid/crates/compressible` (not modified here) without adding a
dependency, build a NASG with `b = q = 0`, `p_inf = K/γ` so that at
`p = 0` its sound speed is `c0 = √(γ p_inf/ρ) = √(K/ρ)`, then set
`Gpr { eos, cs: c_s, tau: f64::INFINITY, tau_model: TauModel::Constant }`
for a lossless elastic solid. For copper this gives `ρ = 8960 kg/m³`,
`c_s = 2299 m/s`, `c_b = 3911 m/s`, `v_p = 4727 m/s`.

### Emergent boundary: grain size vs wavelength

`scattering_regime(d, λ)` returns Rayleigh (`d/λ < 0.1`), stochastic
(`0.1 < d/λ < 10`) or geometric (`d/λ > 10`). The isotropic polycrystal
is valid (up to a scattering loss) only in the Rayleigh regime; in the
stochastic and geometric regimes the single-crystal anisotropy must be
unfolded grain by grain. The thresholds are conventional
order-of-magnitude boundaries, not sharp transitions.

The leading Rayleigh longitudinal attenuation is
`α_L = (8π³/375) (d³ f⁴ / v_L⁴) · Var[M(n)]/⟨M⟩²` with the cubic
longitudinal modulus `M(n) = C_{ijkl}n_i n_j n_k n_l` and
`Var[M]/⟨M⟩² = (16/525)(C11−C12−2C44)²/⟨M⟩²`. It scales as `d³ f⁴`
(checked exactly) and grows with anisotropy (Cu `α_L = 6.2e-3 m⁻¹`,
Al `6.2e-5 m⁻¹` at `d = 0.1 mm`, `f = 5 MHz`; an isotropic crystal has
zero). **This is approximate**: it uses the leading Born/scattering
scaling with the orientation variance of the longitudinal modulus as the
anisotropy factor, but the order-unity prefactor and the full angular
integrals of Stanke & Kino (1984) / Weaver (1990) are not reproduced. It
is only valid for `d ≪ λ`.

## Limits and honest misses

- The elastic tensor is the second derivative of the energy under a
  prescribed deformation gradient, i.e. the Born tensor. At nonzero
  stress it must be used together with the stress (Huang) terms to
  predict sound speeds; the tests therefore compare phonons and
  Christoffel speeds at zero stress. Adding the stress correction would
  let the fold be used at pressure.
- The hcp LJ ideal `c/a` carries a small deviatoric stress (ideal is not
  the exact LJ minimiser at fixed `a`), which shifts the shear constant;
  the zero-stress tests use either the isotropic fcc equilibrium or an
  ideal hcp built from relaxed nearest-neighbour springs.
- Pair potentials miss covalent and angular forces: diamond is modelled
  as central-force floppy, and the Morse/LJ diamond equilibrium is not
  physical. That is the next rung (bond-order / angular potentials), not
  a defect of this one.
- The lattice sums use a spherical cutoff; the tail is a homogeneous
  mean-field estimate, whose residual is the `rmax`-convergence table
  above. The bond list is the true periodic sum but is not minimum-image,
  by design.
- The polycrystal Hashin–Shtrikman shear bound uses the two-phase formula
  with the cubic shear-mode weights; it is not the full
  Peselnick–Meister angular integral, and the Hill average can fall a
  hundredth of a percent outside it (Al). The Rayleigh attenuation keeps
  the Stanke–Kino scaling but only a leading-order anisotropy factor, so
  its prefactor is approximate to an order-unity constant.
- The single-crystal table quotes literature constants without propagating
  their measurement uncertainty; the `estimate` flag only marks values
  that are rounded handbook figures.

## Constants and sources

- fcc lattice sums `A6 = 14.45392`, `A12 = 12.13188`: standard tabulated
  values (Jones 1924; Lennard-Jones & Ingham 1925).
- ideal hcp `c/a = sqrt(8/3)`.
- Cauchy relation for central-force Bravais crystals at zero pressure:
  Born & Huang, *Dynamical Theory of Crystal Lattices*.
- Christoffel equation and the 1D chain dispersion: standard crystal
  acoustics / lattice-dynamics texts.
- Argon LJ parameters (`σ = 0.3405 nm`, `ε/k_B = 119.8 K`) are included
  in `potential.rs` with the Rowley–Nicholson source; the benchmarks above
  use dimensionless `σ = ε = 1` to match the exact lattice sums.
- Single-crystal elastic constants: Simmons & Wang, *Single Crystal
  Elastic Constants and Calculated Aggregate Properties: A Handbook*
  (MIT Press, 1971); silicon from Hopcroft, Nix & Kenny, *J.
  Microelectromech. Syst.* **19**, 229 (2010); densities from the CRC
  Handbook. The zinc entry is a rounded handbook estimate (`estimate`
  flag).
- Bond matrices and crystal-orientation rotation: Nye, *Physical
  Properties of Crystals*.
- Voigt/Reuss/Hill averages: Voigt (1910); Reuss (1929); Hill (1952).
- Hashin–Shtrikman bounds: Hashin & Shtrikman, *J. Mech. Phys. Solids*
  **11**, 127 (1963); for cubic polycrystals, Watt, Davies & O'Connell,
  *Rev. Geophys. Space Phys.* **14**, 541 (1976); Peselnick & Meister,
  *J. Geophys. Res.* **80**, 4017 (1975).
- Rayleigh attenuation scaling and cubic anisotropy: Stanke & Kino,
  *J. Acoust. Soc. Am.* **75**, 665 (1984); Weaver, *J. Mech. Phys.
  Solids* **38**, 55 (1990).
