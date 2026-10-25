# The defect rung: point and planar defects, and the fold to plasticity

This rung sits between the perfect-crystal Cauchy–Born rung
([materials.md](materials.md)) and continuum plasticity. The fine model is
a small periodic atomistic supercell with a central pair potential and
relaxed atom positions ([`materials::defects`](../fluid/crates/materials/src/defects.rs));
the coarse model is the homogeneous continuum, and above it the
elasto-plastic GPR solid of `compressible/gpr.rs`. The fold is
Peierls–Nabarro: a dislocation is a soliton of the γ-surface, its Peierls
stress is the single-crystal CRSS, and a Taylor factor maps that onto the
GPR yield stress ([`materials::plasticity`](../fluid/crates/materials/src/plasticity.rs)).

The organizing pattern is the one in
[emergence_ladder.md](emergence_ladder.md): a cheap continuum yield law, a
validity indicator, and an unfold where the indicator fails.

- **coarse model** — the GPR solid with `TauModel::Plastic { tau0,
  sigma_y, n }` (a yield stress above which the shear strain relaxes);
- **validity indicator** — the Peierls stress `σ_P` and the sample-size /
  dislocation-spacing ratio `L / (1/√ρ)`;
- **fold / unfold** — below `σ_P` the material is elastic; above it and
  when `L ≫ 1/√ρ` the dislocation ensemble is statistical and the GPR
  yield law applies; when `L ≲ 1/√ρ` individual dislocation cores and
  source-limited "smaller is stronger" strengthening must be unfolded.

All atomistic numbers below are for the 12–6 Lennard-Jones potential with
`σ = ε = 1`, at its zero-pressure fcc equilibrium
`a = 1.5417373 σ`, `r_nn = 1.0901729 σ`, `e_coh = −8.6102077 ε`. The
lattice sums and the cohesive energy are the validated values of
`analysis/materials.md` (Jones 1924; Lennard-Jones & Ingham 1925).

The minimiser is Polak–Ribière conjugate gradients with an Armijo
backtracking line search, driven by analytic forces. The neighbour topology
is frozen (exact for the small relaxation displacements). The whole defect
suite (`fluid/crates/materials/tests/defects.rs`) runs in about 1.7 s.

## 1. Vacancy

For a **pair potential** the *unrelaxed* vacancy formation energy is an
exact identity. The cohesive energy per atom is `e_coh = ½ Σ_j φ`, while
the total crystal energy is `½ Σ_i Σ_j φ`. Removing one atom removes all
its bonds, a per-atom sum `Σ_j φ = 2 e_coh`; referring the remaining `N−1`
atoms to the bulk cohesive energy gives

```text
E_v^unrelaxed = − e_coh.
```

For LJ fcc this is exactly `8.6102 ε`; the supercell calculation reproduces
it to the finite-size/cutoff residual. Full relaxation of the neighbours
lowers it.

| quantity | reference | computed | note |
|---|---|---|---|
| unrelaxed `E_v` | `8.6102 ε` (exact identity) | `8.6050 ε` | 4×4×4 supercell, `rc = 6σ` (0.06 % finite-size) |
| relaxed `E_v` | — | `8.4894 ε` | 6×6×6, `rc = 6σ` |
| relaxation energy | — | `0.120 ε` (1.4 %) | finite-size converged |
| truncated-LJ `E_v` | `0.0690 eV = 7.6 ε` (Beeler & Beeler 1981) | `8.49 ε` (full range) | their `rc = 1.94 r0` removes the attractive tail, lowering `E_v` |

The relaxed value is stable against supercell size: 4×4×4 → `8.4942 ε`,
5×5×5 → `8.4906 ε`, 6×6×6 → `8.4894 ε`. The published truncated-potential
calculation of Beeler & Beeler, *Vacancy Behavior in a Compressed fcc
Lennard-Jones Crystal* (OSTI 5669586, 1981), quotes `E_f(SCV) = 0.0690 eV`
with their `ε = 9.1×10⁻³ eV`; the difference from our full-range value is
the missing `r⁻⁶` tail, i.e. a **cutoff** effect, not a physical
disagreement. Pair potentials are known to overestimate vacancy energies
because they lack directional bonding.

## 2. Free surfaces

The rigid-lattice surface energy of a central-force crystal has an exact
**broken-bond** form. Cleaving at a plane removes every bond with
`h = r·n > 0`; one bond of range `h` is removed for a slab layer of
thickness `h`, so the removed energy per area is `ρ Σ_{r·n>0} h |φ(r)|`,
shared by the two surfaces:

```text
γ(n) = ½ ρ Σ_{r·n>0} (r·n) |φ(r)|.
```

`materials::defects::broken_bond_surface_energy` evaluates it on the true
periodic bond list. The relaxed slab
`γ = (E_slab − N e_bulk)/(2A)` reproduces the rigid value to **1×10⁻⁶
relative** (the bond lists are identical), which validates the slab
construction and the vacuum handling.

| face | `rc` | rigid broken-bond `γ` | relaxed slab `γ` | relaxation |
|---|---|---|---|---|
| (111) | 5 σ | 2.9820 | 2.9709 | 0.37 % |
| (111) | 6 σ | 3.0207 | 3.0100 | 0.35 % |
| (111) | 12 σ | 3.1003 | — | — |
| (100) | 5 σ | 3.1017 | 3.0716 | 0.97 % |
| (100) | 6 σ | 3.1399 | 3.1106 | 0.93 % |
| (100) | 12 σ | 3.2197 | — | — |

All energies in `ε/σ²`. The exact rigid sums converge slowly (the `r⁻⁶`
tail contributes `∝ rc⁻²`); at `rc = 12 σ` the (111) value is `3.100` and
the (100) value `3.220`, so the infinite-range relaxed values are about
`γ(111) ≈ 3.09` and `γ(100) ≈ 3.19 ε/σ²`. The close-packed (111) face is
the lower-energy one, `γ(100)/γ(111) = 1.03–1.04`.

This agrees with the classical result of Shuttleworth, *Proc. Phys. Soc.
A* **62**, 167 (1949), who found the inert-gas surface energy
proportional to the sublimation energy, **independent of face to within
15 %**, and a `{100}` relaxation of `+2.54 %` outward making a **negative
`1 %`** contribution to the surface energy. We find the (100) relaxation
lowers `γ` by 0.93 %, matching that statement. The modern atomistic
reference for the temperature dependence and the full curvature is
Broughton & Gilmer, *Acta Metall.* **31**, 845 (1983).

## 3. Stacking faults and the γ-surface

For a pair potential the **intrinsic stacking-fault energy** follows from
the hcp–fcc energy difference. The hcp stacking `ABAB` contains one
reversed step per two close-packed planes while fcc contains none, so one
fault costs half the hcp–fcc difference per plane:

```text
γ_ISF = 2 (e_hcp − e_fcc) / A_{111},
```

with `A_{111} = √3 a²/4 = 1.0293 σ²` the area per atom. The materials
crate gives `e_hcp − e_fcc` (tail-corrected) approaching the
`analysis/materials.md` value `−9.0×10⁻⁴ ε`:

| `rc` | `e_hcp − e_fcc` (`ε`) | `γ_ISF` (`ε/σ²`) |
|---|---|---|
| 6 σ | `−4.91×10⁻⁴` | `−9.53×10⁻⁴` |
| 8 σ | `−1.09×10⁻³` | `−2.12×10⁻³` |
| 10 σ | `−8.22×10⁻⁴` | `−1.60×10⁻³` |
| 12 σ | `−8.63×10⁻⁴` | `−1.68×10⁻³` |

The sign is **negative**: LJ puts hcp marginally below fcc, so a stacking
fault lowers the energy. This is the well-known failure of pair potentials
(real fcc metals have positive stacking-fault energies of tens of
mJ/m², i.e. `O(0.01–0.1) ε/σ²`). The direct constrained (normal-only)
fault at `t = 1` on the `⟨112⟩` path gives `+0.0016 ε/σ²` at `rc = 6σ`,
i.e. the two independent routes agree that `|γ_ISF| ≲ 0.002 ε/σ²`.

The **generalized stacking-fault energy** `γ(s)` is computed by rigidly
shifting the upper half of a `(111)` slab by `s = t·dir` and relaxing only
the normal coordinate (Vitek's constrained γ-surface convention; Vitek,
*Phil. Mag.* **18**, 773 (1968)). The unstable stacking-fault energy
`γ_us` is the maximum for `t ∈ (0, 1]`.

| path | `γ_us` rigid | `γ_us` relaxed | peak at |
|---|---|---|---|
| `⟨112⟩` (to the Shockley partial) | `0.5514` | `0.4226` | `t = 0.5` |
| `⟨110⟩` (full lattice period) | `3.226` | `1.251` | `t = 0.5` |

Energies in `ε/σ²`, `rc = 6σ`, 9-layer slab, 12 sample points. The
`⟨110⟩` barrier is about 3× the `⟨112⟩` barrier, as expected. Relaxation
lowers `γ_us` by ~23 % on `⟨112⟩`.

## 4. Peierls–Nabarro: a crude closure

For an edge dislocation with Burgers vector `b` on a slip plane of spacing
`d`, the classical core width and Peierls stress are

```text
ζ = d / (2 (1 − ν)),      σ_P ≈ (2 G / (1 − ν)) exp(−2 π ζ / b).
```

With a **sinusoidal** γ-surface `γ(φ) = (γ_us/2)(1 − cos(2πφ/b))` the
continuum PN equation has the exact soliton

```text
φ(x) = b/2 + (b/π) arctan((x − x0)/ζ),
```

with `ζ = G b² / (4π²(1−ν) γ_us)`; the elastic amplitude
`γ_us = G b²/(2π²d)` reproduces the classical width. The continuum
soliton is translation invariant, so the exponential Peierls stress comes
from the discreteness of the lattice, not from the continuum solution.

Using the extrapolated LJ fcc elastic constants
(`C11 = 105.182`, `C12 = C44 = 60.187 ε/σ³`; `analysis/materials.md`) and
the Hill average `G = 40.575 ε/σ³`, `ν = 0.2713`, `b = 1.0902 σ`,
`d = 0.8901 σ`:

| quantity | value |
|---|---|
| core width `ζ = d/(2(1−ν))` | `0.6108 σ` |
| Peierls stress `σ_P` | `3.296 ε/σ³ = 0.0812 G` |
| elastic sinusoidal `γ_us = G b²/(2π²d)` | `2.744 ε/σ²` |
| fault-plane ideal strength `π γ_us / b` (elastic) | `7.91 ε/σ³ = 0.195 G` |
| **measured** `γ_us` (relaxed `⟨112⟩`) | `0.4226 ε/σ²` |
| width implied by the measured `γ_us` | `3.97 σ` |
| fault-plane ideal strength (measured `γ_us`) | `1.22 ε/σ³ = 0.030 G` |

**PN is not a predictive closure.** The atomistic `γ_us` is 6.5× smaller
than the elastic/sinusoidal amplitude that would reproduce the classical
core width, and the width it implies (`3.97 σ`) is 6.5× the classical
`ζ`. The exponential estimate `0.081 G` and the measured-`γ_us` ideal
strength `0.030 G` differ by 2.7×, and the exponential is exponentially
sensitive to `ζ` (a 10 % change in `ζ` changes `σ_P` by ~30 %). The
classical formula should be read as an order-of-magnitude estimate only
(Schoeck, *Mater. Sci. Eng. A* **400–401**, 7 (2005)).

## 5. The fold to continuum plasticity and the GPR closure

The polycrystal yield stress is `σ_y = M_T τ_CRSS` with the random-fcc
Taylor factor `M_T = 3.06` (Taylor 1938; Kocks, *Metall. Trans.* **1**,
1121 (1970)). Taking the PN stress as the single-crystal CRSS gives
`σ_y = 3.06 × 3.296 = 10.09 ε/σ³ = 0.249 G`. Alternatively the Taylor
hardening law `τ = α G b √ρ` (Taylor 1934) supplies a density-dependent
CRSS; with `α = 0.3`, `b = 1.09 σ`, `G = 40.57 ε/σ³` and
`ρ = 10⁻³ σ⁻³` it gives `τ = 0.420 ε/σ³`.

`materials::plasticity` returns a `GprPlastic { tau0, sigma_y, n }` whose
fields map one-to-one onto `compressible::gpr::TauModel::Plastic`
(**read, not modified**):

```text
GprPlastic { tau0, sigma_y, n }  <->  TauModel::Plastic { tau0, sigma_y, n }
sigma_y = M_T tau_CRSS
tau0    = M_T B / (G rho b^2)
```

The `tau0` mapping: above yield the GPR deviatoric stress relaxes as
`dσ/dt ≈ −σ/τ0`, so the plastic strain rate is `γ̇_GPR = σ/(G τ0)`;
Orowan glide gives `γ̇ = ρ b² τ_res / B` with drag `B` and resolved stress
`τ_res = σ/M_T`. Equating the two eliminates `σ` and gives the relation.
The exponent `n` is a hardening parameter (large `n` is nearly ideal
plasticity; `n → ∞` makes the GPR radial return project exactly onto the
yield surface).

### Validity indicator

The continuum-plasticity rung is valid when both:
1. `σ ≥ σ_P` — below the Peierls stress the dislocations are pinned and the
   response is elastic, so the plastic closure is inactive;
2. `L ≥ 1/√ρ` — the smallest stressed dimension exceeds the mean
   dislocation spacing, so the discrete cores average into a density.

`materials::plasticity::plasticity_rung(L, ρ, σ, σ_P)` returns `Elastic`,
`Discrete` or `ContinuumPlastic`. When `L ≲ 1/√ρ` the continuum density is
meaningless and the **discrete** rung applies: source-limited strengthening
(`σ ≈ G b / L`, the Orowan/Stokes–Nabarro estimate; Uchic et al., *Science*
**305**, 986 (2004); Greer & Nix, *Phys. Rev. B* **73**, 245410 (2006))
gives the "smaller is stronger" size effect, and polycrystals add the
Hall–Petch term `k/√d` (Hall 1951; Petch 1953). Both are provided as
estimates, not calibrated laws. For LJ fcc, `σ_P/G = 0.081`, so the
elastic→plastic boundary sits at about `8 %` of the shear modulus — an
order of magnitude above real fcc metals (where `σ_P/G ~ 10⁻³–10⁻⁴`),
another symptom of the pair-potential closure.

## Limits and honest misses

- **Vacancy**: the full-range relaxed value is `8.49 ε`; pair potentials
  overestimate vacancy energies because they lack bond directionality. The
  finite supercell residual is a few `10⁻³ ε`; the published truncated-LJ
  result (`7.6 ε` in its own units) is dominated by the missing tail.
- **Surfaces**: the rigid broken-bond sum is exact and validated against
  the slab to `10⁻⁶` relative; the relaxed values are quoted at the
  same-cutoff references and converge slowly (`∝ rc⁻²`). The classical
  comparison is the face-independence and the `~1 %` (100) relaxation of
  Shuttleworth (1949), reproduced here.
- **Stacking fault**: the pair-potential `γ_ISF` is **negative**, a
  qualitative failure. The hcp–fcc and direct-slab routes agree only in
  magnitude (`|γ_ISF| ≲ 0.002 ε/σ²`); the direct constrained fault is
  cutoff-sensitive at exactly this scale.
- **γ-surface**: only the normal coordinate is relaxed (Vitek's
  constrained convention); full in-plane relaxation would lower `γ_us`
  further and is not attempted.
- **PN**: crude, as documented above. The exponential form is the leading
  discrete-lattice result and its prefactor is not universally agreed.
- **Isotropic moduli**: the PN fold uses the Voigt–Reuss–Hill average of
  the cubic crystal, whereas a real dislocation sees the anisotropic
  single-crystal moduli. The exponent makes this a large uncertainty.
- **GPR mapping**: `tau0` needs a dislocation drag coefficient `B`, which
  is not computed here; it is a labelled external closure. The mapping
  assumes dilute, non-interacting glide and a single mobile density.

## Constants and sources

- fcc LJ equilibrium, cohesive energy, lattice sums, extrapolated elastic
  constants: `analysis/materials.md` (Jones 1924; Lennard-Jones & Ingham
  1925).
- Unrelaxed vacancy identity `E_v = −e_coh`: standard pair-potential
  defect counting; the same recipe is written in Beeler & Beeler (1981).
- Vacancy in a compressed fcc LJ crystal: R. J. Beeler & M. F. Beeler,
  *Vacancy Behavior in a Compressed fcc Lennard-Jones Crystal*,
  OSTI 5669586 (1981) — `E_f(SCV) = 0.0690 eV` at `rc = 1.94 r0`.
- Surface energy of inert-gas crystals: Shuttleworth, *Proc. Phys. Soc. A*
  **62**, 167 (1949) — face-independence to 15 %, `{100}` relaxation
  `+2.54 %`, `1 %` energy contribution. Modern atomistic study:
  Broughton & Gilmer, *Acta Metall.* **31**, 845 (1983).
- Generalized stacking-fault / γ-surface and constrained relaxation:
  Vitek, *Phil. Mag.* **18**, 773 (1968).
- Peierls–Nabarro model: Peierls, *Proc. Phys. Soc.* **52**, 34 (1940);
  Nabarro, *Proc. Phys. Soc.* **59**, 256 (1947); Hirth & Lothe, *Theory of
  Dislocations*; critical review in Schoeck, *Mater. Sci. Eng. A*
  **400–401**, 7 (2005).
- Taylor factor and hardening: Taylor, *Proc. R. Soc. A* **145**, 362
  (1934); Kocks, *Metall. Trans.* **1**, 1121 (1970) (fcc `M ≈ 3.06`).
- Size effects: Uchic et al., *Science* **305**, 986 (2004); Greer & Nix,
  *Phys. Rev. B* **73**, 245410 (2006); Hall, *Proc. Phys. Soc. B* **64**,
  747 (1951); Petch, *J. Iron Steel Inst.* **174**, 25 (1953).
- GPR elasto-plastic closure: Peshkov & Romenski (2016); Dumbser et al.
  (2016), implemented in `compressible/gpr.rs`.
- Conjugate-gradient minimiser: Polak & Ribière (1969); Armijo (1966).
