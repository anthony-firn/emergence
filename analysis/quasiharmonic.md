# Quasi-harmonic ladder: atoms → thermal properties

Implemented in `fluid/crates/materials/src/quasiharmonic.rs`, with the
benchmarks in `fluid/crates/materials/tests/quasiharmonic.rs`. This is the
rung above the static crystal of [materials.md](materials.md): the fine
model is the full harmonic phonon spectrum on a Brillouin-zone grid, the
coarse model is a temperature-dependent equation of state obtained by
minimising the free energy over volume, and the fold between them is the
quasi-harmonic approximation (QHA), in which the only volume dependence of
the vibrations is through the harmonic frequencies. The two validity
indicators are the quantum/classical crossover `T/θ_D` and the Lindemann
ratio `√⟨u²⟩/d_nn`; where they fail, the model must unfold to molecular
dynamics (or another anharmonic/closure rung).

## The rung

| piece | content |
|---|---|
| Monkhorst–Pack grid | `(n1,n2,n3)` mesh in the primitive reciprocal cell, equal weights, no symmetry reduction |
| density of states | histogram of `3N_q` frequencies, normalised to 3 per atom |
| harmonic thermodynamics | `F_vib(V,T) = (1/N_q) Σ_{qν} [ħω/2 + k_B T ln(1 − e^{−ħω/k_BT})]`, internal energy, `C_V`, entropy; the classical limit `k_B T ln(ħω/k_BT)` is a switch |
| QHA | minimise `F(V,T) = E_static(V) + F_vib(V,T)` over uniform volume at zero pressure |
| derived quantities | `α(T)`, `B_T(T)`, `C_V(T)`, `C_P − C_V = α² B_T V T`, `γ = α B_T V / C_V`, mode Grüneisen `−∂lnω/∂lnV` |
| indicators | `T/θ_D`, Lindemann ratio, distance to the harmonic spinodal |

The dynamical matrix and bond list are those of the existing
`phonon.rs`/`lattice.rs`; the Brillouin-zone sums use a closed-form 3×3
Hermitian eigensolver (a 12-digit match to the generic Jacobi solver is
checked in a unit test). One practical constraint inherited from the
shared bond-list builder is that it discards bonds closer than an absolute
`1e-9`, so the argon benchmark runs in reduced LJ units (`σ = ε = m = 1`)
and converts with `σ = 3.405 Å`, `ε/k_B = 119.8 K`, `m = 39.948 u`.

## Benchmarks and measured errors

All numbers are produced by `cargo test -p materials` (the
quasi-harmonic suite runs in ≈10 s). Reduced LJ runs use `σ = ε = m = 1`
and `ħ = k_B = 1`, so `T` is in units of `ε/k_B`.

### Brillouin-zone grid and DOS

| quantity | reference | measured |
|---|---|---|
| reciprocal orthogonality `a_i·b_j = 2π δ_ij` | exact | < 1e-10 |
| even mesh contains Γ | no | no |
| odd mesh contains Γ | yes | yes |
| `∫ g(ω) dω` per atom (grid 8³) | `3 − 3/N_q` | 3 to < 0.6 % |

### Dulong–Petit (high-temperature limit)

At the static volume, grid 12³, dimensionless LJ fcc:

| `T*(=ε/k_B)` | `C_V/k_B` | error from 3 |
|---|---|---|
| 80 | 2.98516 | 0.49 % |
| 160 | 2.99628 | 0.124 % |
| 320 | 2.99907 | 0.031 % |
| any | classical branch | exact (`< 1e-12`) |

### Debye `T³` law (dimensionless LJ fcc)

`C_V = (12π⁴/5) k_B (T/θ_D)³`. Extracting `θ_D` from the low-temperature
coefficient (grid 32³, `T/θ_D ≈ 0.1`, where the zone sum resolves the
acoustic branches) and comparing with the independent Debye temperature
from the elastic tensor and the directional-average Christoffel sound
speed at the same lattice:

| route | `θ_D` (natural units) |
|---|---|
| low-`T` `C_V/T³` fit (`T* = 2.5, 3.0, 3.5`) | 26.75 |
| elasticity + average sound speed | 26.34 |
| relative error | **+1.6 %** |

The `C_V/T³` plateau is flat to 4 % over the fit range. Below
`T* ≈ 2` the discrete zone sum under-resolves the acoustic modes and the
extracted `θ_D` drifts upward, which is the usual finite-grid limitation
of a direct low-`T` Debye extraction.

### Classical-limit QHA thermal expansion

Classical statistics (`ħ → 0`), dimensionless LJ fcc, grid 6³. At low
temperature the equilibrium volume is linear in `T`, so `α` is a
constant:

| `T*` | `s` | `α` |
|---|---|---|
| 0.02 | 1.5432 | 0.1345 |
| 0.05 | 1.5453 | 0.1386 |
| 0.10 | 1.5490 | 0.1461 |
| 0.15 | 1.5529 | 0.1547 |

(The static minimum of the QHA energy, with the cutoff scaled with the
lattice, is `s ≈ 1.5432`, 0.13 % below the `static_scale` found by
`potential::relax`, which uses a fixed cutoff. This is a cutoff-consistency
detail of the existing API, not a QHA result.) The two independent
Grüneisen estimates `α B_T V / C_V` and the mode average of
`−∂lnω/∂lnV` agree to < 0.05 %.

### Solid argon

Reduced-unit LJ fcc argon (`σ = 3.405 Å`, `ε/k_B = 119.8 K`,
`m = 39.948 u`; Rowley & Nicholson 1972 for the potential; grid 6³).

**Zero-point lattice expansion.** The static minimum is `a = 5.2501 Å`
(the converged LJ value is 5.2496 Å). The quantum QHA at `T → 0` gives
`a = 5.3141 Å`: the zero-point energy expands the lattice by **0.064 Å**
(1.2 %), and the value agrees with the experimental `a(4.2 K) = 5.311 Å` to
0.06 %. This is a genuine quantum effect that the classical fold does not
contain. (An earlier version quoted a 0.053 Å expansion from a static
reference of 5.2611 Å; that reference came from the tail-at-fixed-cutoff
bug in `potential::relax`, fixed in the materials crate.)

**Lattice constant vs temperature.** Experimental values are rounded from
Peterson, Batchelder & Simmons (1966); they are approximate to the last
digit and are treated as anchors, not as a fit target.

| `T` (K) | `a_QHA` (Å) | `a_exp` (Å) | error |
|---|---|---|---|
| 4.2 | 5.3142 | 5.311 | +0.06 % |
| 10 | 5.3152 | 5.312 | +0.06 % |
| 20 | 5.3239 | 5.317 | +0.13 % |
| 40 | 5.3643 | 5.342 | +0.42 % |
| 50 | 5.3944 | 5.362 | +0.60 % |
| 60 | 5.4090 | 5.389 | +0.37 % |
| 80 | 5.4090 | 5.462 | −0.97 % |

The error grows monotonically with `T` (QHA over-expands), then the
equilibrium is pinned at the harmonic spinodal and the `a(T)` curve stops
responding. The columns beyond 50 K are inside the breakdown region and
are quoted only to show the failure.

**Debye temperature.** Elasticity plus the directional-average sound speed
give `θ_D(Ar) = 93.29 K`; the low-temperature heat-capacity value for
solid argon is ≈ 93.5 K (standard rare-gas table). Error **−0.2 %**.
The direct low-`T` `C_V` extraction is not useful for argon: `θ_D` is only
`0.78 ε/k_B`, so the required `T ≪ θ_D` lies below what a feasible
Brillouin-zone grid can resolve. This is reported as a negative result
for the direct route, and the elastic route (which uses the same
long-wavelength limit) is the reliable one.

**Bulk modulus (isothermal, zero pressure).** `B_T` falls from 2.65 GPa
at 4 K to 2.48 GPa at 20 K and 1.96 GPa at 40 K, a 26 % softening over
that range; the low-`T` value is the right order for solid argon
(≈1.6–2.7 GPa). Near the spinodal the one-sided finite difference is no
longer physically meaningful.

**Thermodynamic Grüneisen parameter.** `γ = α B_T V / C_V` in the quantum
QHA runs 3.29 (4 K), 3.10 (20 K) and 3.16 (40 K). The mode
heat-capacity average agrees with it. The commonly quoted zero-pressure
value for rare-gas/LJ fcc is ≈2.6–3.0; the calculation lands
**+4 % to +20 % above** that band, consistently with the LJ pair
potential overestimating the volume dependence (real rare gases have
many-body contributions). This is a miss on the upper edge, not a
disguised success.

## Emergent boundaries

1. **Quantum → classical crossover.** The indicator is `T/θ_D`. For
   `T/θ_D ≳ 1` the classical free energy `k_B T ln(ħω/k_BT)` and
   `C_V = 3 k_B` are valid; the largest deviation from Dulong–Petit in
   the table above is 0.49 % at `T/θ_D ≈ 3`. For `T/θ_D ≪ 1` the
   zero-point term dominates and the classical fold is qualitatively
   wrong (it has no zero-point expansion and no `T³` law). A
   molecular-dynamics rung would take over exactly at the upper edge of
   this window, where the classical potential-energy surface is the only
   input needed.
2. **Harmonic → anharmonic.** The indicator is the Lindemann ratio
   `√⟨u²⟩/d_nn`. For argon it is 0.051 at 4 K, 0.082 at 40 K, 0.094 at
   50 K and crosses the traditional ≈0.1 melting threshold between 50 and
   60 K. At the same point the equilibrium volume reaches the harmonic
   spinodal (`a = 5.404 Å`): the soft acoustic mode makes the QHA free
   energy diverge and the crystal is mechanically unstable. Real argon
   melts at 83.8 K, so the QHA has already declared defeat ~25–30 K early;
   that gap is the anharmonic stabilization a molecular-dynamics rung
   would supply.
3. **Spinodal indicator.** The module computes the spinodal once (the
   first volume where `min ω² ≤ 1e-4 max ω²` on a coarse grid); beyond it
   `F(V,T)` is set to `+∞` so the volume minimiser cannot enter the
   unstable region. `ThermalProperties::near_instability` marks a
   temperature whose equilibrium sits on that boundary.

## Limits and honest misses

- The **direct low-`T` Debye extraction is grid-limited**: for a high
  Debye temperature (`θ_D ≈ 26` in reduced units) the zone sum needs
  `T ≫ ω_min` and the interesting `T ≲ θ_D/10` is only reachable with a
  few tens of points per axis. The extracted `θ_D` is good to ≈2 % at
  `T ≈ θ_D/9` and degrades below that. For argon the gap is
  unfixable at present, and only the elastic route is quoted.
- The **Grüneisen parameter is ≈3.1**, at the high end of / above the
  rare-gas band 2.6–3.0. A pair potential lacks the three-body terms that
  soften the volume dependence; the discrepancy is intrinsic to this
  rung.
- The **harmonic spinodal precedes melting** by ~25 K in LJ argon, so
  the QHA cannot be used through the melting line. This is expected: the
  harmonic model has no mechanism to stabilize the expanding crystal.
- **Multi-atom bases** use the generic Jacobi eigensolver and the
  mean-square displacement uses an averaged mass; the fast 3×3 path and
  the quantitative Lindemann benchmark are for monatomic Bravais cells.
- **The bond-list cutoff is scaled with volume** in the free energy
  (`rmax·s`) so `F(V,T)` is smooth, but the shared bond-list builder has
  an absolute `1e-9` exclusion that forces atomic-scale SI models to be
  run in reduced units. Fixing that belongs to `lattice.rs`, which is
  outside this rung.
- The experimental argon lattice constants are rounded and should be
  treated as anchors, not as a precision measurement; the trend (QHA
  over-expanding with rising `T`, then failing) is the robust conclusion.

## Constants and sources

| constant | value | source |
|---|---|---|
| argon `σ` | 0.3405 nm | Rowley & Nicholson 1972 (via `potential.rs`) |
| argon `ε/k_B` | 119.8 K | Rowley & Nicholson 1972 |
| argon mass | 39.948 u | standard atomic weight |
| `k_B` | 1.380649e-23 J/K | SI 2019 |
| `ħ` | 1.054571817e-34 J·s | SI 2019 |
| Monkhorst–Pack mesh | `u_i = (2m−n_i−1)/(2n_i)` | Monkhorst & Pack, Phys. Rev. B 13, 5188 (1976) |
| Debye `C_V` coefficient | `12π⁴/5` | standard Debye model |
| Debye temperature from elasticity | `θ_D = (ħ v_D/k_B)(6π²n)^{1/3}` | standard; `3/v_D³ = ⟨Σ_s v_s⁻³⟩` |
| solid-argon `a(T)` | table above | Peterson, Batchelder & Simmons 1966 (rounded) |
| solid-argon `θ_D` | ≈93.5 K | standard rare-gas Debye table |
| Lindemann ratio ≈0.1 at melting | empirical | Lindemann 1910 |
| corrected `C_P−C_V = α²B_TVT` | exact | standard thermodynamic identity |
