# The quasicontinuum rung: adaptive atomistic/continuum coupling

This rung is the fold/unfold principle of
[emergence_ladder.md](emergence_ladder.md) applied to a solid.
The fine model is a crystal of atoms with a central pair potential, the
coarse model is Cauchy-Born elasticity on a linear finite-element mesh,
and the fold between them is the **quasicontinuum (QC)** method of Tadmor,
Ortiz & Phillips, *Phil. Mag. A* **73**, 1529 (1996); the modern review is
Miller & Tadmor, *Modelling Simul. Mater. Sci. Eng.* **17**, 053001 (2009).

- **coarse model** - a piecewise-linear displacement field on a coarse
  simplex mesh, with the Cauchy-Born energy density `W(F)` of the crystal;
- **validity indicator** - the local variation of the deformation gradient:
  the jump of `F` between neighbouring elements, or (sharper) the
  Cauchy-Born error `|E_i^atom - E_i^CB|` per element;
- **unfold** - where the indicator is large, a coarse cell is promoted to
  fully atomistic (every atom becomes a degree of freedom); where it is
  small the cell keeps only the finite-element nodes as degrees of freedom.

The implementation is `fluid/crates/materials/src/quasicontinuum.rs`; the
benchmarks are `fluid/crates/materials/tests/quasicontinuum.rs`. The whole
QC suite runs in about 20 s.

## 1. The model

A **repatom** is a degree of freedom. Every atom inside the atomistic
region is a repatom; the continuum carries only the finite-element nodes
of the coarse mesh; every other atom's displacement is the linear
interpolation `u_i = Σ_r N_{ir} u_r` of the element that contains it. The
number of degrees of freedom is therefore far smaller than the number of
atoms.

For a pair potential the per-atom site energy is
`E_i = ½ Σ_{j∈N(i)} φ(|r_ij|)`. Two couplings are provided.

### Naive energy-based QC (the ghost-force problem)

The textbook energy is

```text
E_naive = Σ_{i ∈ A} E_i^atom  +  Σ_{e ∈ C} V_e W(F_e),
```

with `A` the atomistic atoms and `C` the continuum elements. An interface
atom is counted by both the atomistic site sum and the adjacent continuum
element. Under a uniform deformation the two agree, but the *derivatives*
do not, and a uniform-strain patch test leaves a spurious **ghost force**
at the interface. This is the well-known failure of sharp energy-based
coupling (Miller & Tadmor 2009, §4).

### Quasi-nonlocal (QNL) coupling (the fix)

The QNL energy of Shimokawa, Mortensen, Schiøtz & Jacobsen,
*Phys. Rev. B* **69**, 214102 (2004) sums every bond exactly once, at the
positions of the interpolated displacement field:

```text
E_QNL = Σ_bonds φ(|R + u_j − u_i|),
```

where `u` is exact for repatoms and interpolated for continuum atoms. For
a linear element the interpolated difference of two atoms inside the
element is exactly the Cauchy-Born bond length `F R`, so the continuum
energy is the Cauchy-Born energy; a near-interface bond between an
atomistic and a continuum atom uses the exact position of the first and
the interpolated position of the second. Under a uniform strain the whole
energy collapses to the exact atomistic energy, so the patch test passes
to round-off. This is the QNL "local rule"; in 1D it is exactly
ghost-force free and in higher dimensions it is patch-test consistent
once the cutoff is not placed exactly on a neighbour shell (see the
limits below).

## 2. Patch tests: ghost force and its removal

The 1D test is a harmonic nearest-neighbour chain at `a0 = 1`,
`φ = ½ k (r−a0)²`, with the atomistic region `[−6a0, 6a0]`, coarsening 4,
and a uniform strain `F = 1.01`. For a harmonic chain the exact atomistic
force on every interior atom vanishes for **any** uniform `F`, so any
interior force in the coupled model is a ghost force. All boundary
repatoms are held fixed.

| coupling | max interior ghost force |
|---|---|
| naive | `7.50 × 10⁻³` |
| QNL | `0` (below `10⁻¹⁰`, machine round-off) |

The 2D test is a triangular LJ sheet (NN + second shell, `a0 = 1.096966 σ`,
`rmax = 1.9 a0`), 19×19 atoms, coarsening 3, atomistic disc radius `3a0`,
uniform strain `F = 1.005`, checked on repatoms more than four lattice
steps from every free surface. The QNL residual is `3.3 × 10⁻¹³`.

The analytic QNL forces are validated against a central difference of the
energy: the largest disagreement is `5.2 × 10⁻⁹` (1D) and
`1.3 × 10⁻⁶` (2D) for a non-uniform test displacement.

> The 2D free-surface patch test must exclude atoms within a few lattice
> steps of the surface: a surface atom has missing bonds, so its atomistic
> force is a physical reaction, not a ghost force, and the finite-element
> interpolation spreads it to the neighbouring nodes. The cutoff must also
> avoid landing exactly on a neighbour shell; with `rmax = 2.0 a0` the
> `2a0` shell is included asymmetrically by floating-point rounding and a
> spurious `O(10⁻¹)` residual appears. Both are numerical, not modelling,
> artefacts and are documented here rather than hidden.

## 3. 1D defect benchmarks

All 1D runs use the 12–6 Lennard-Jones potential with `σ = ε = 1`, at its
zero-stress chain spacing `a0 = 1.119611 σ` (NN + next-nearest neighbours,
`rmax = 2.5 a0`), found by minimising the raw cut-off energy. The two end
atoms are clamped; the reference is the fully atomistic model with the
same atoms and boundary conditions.

### Point load at the centre (`n = 81`, load `0.5`)

Reference: `N = 81` degrees of freedom, `E = −82.51008312`.

| `r_atom / a0` | ndof | `max |Δu|` | `|ΔE|` |
|---|---|---|---|
| 2 | 24 | `7.82 × 10⁻³` | `2.43 × 10⁻⁵` |
| 4 | 27 | `7.70 × 10⁻³` | `2.43 × 10⁻⁵` |
| 8 | 33 | `6.93 × 10⁻³` | `2.30 × 10⁻⁵` |
| 16 | 45 | `3.91 × 10⁻³` | `2.26 × 10⁻⁶` |
| 40 (all) | 81 | — | `< 10⁻⁵` |

Displacement and energy errors fall as the atomistic region grows, and the
all-atomistic coupled model reproduces the reference energy to better than
`10⁻⁵`.

### Vacancy at the centre (`n = 61`)

Reference: `N = 60` atoms, `E = −59.78645059`.

| `r_atom / a0` | ndof | `max |Δu|` | `|ΔE|` |
|---|---|---|---|
| 2 | 18 | `4.23 × 10⁻³` | `6.67 × 10⁻⁵` |
| 4 | 21 | `4.02 × 10⁻³` | `6.68 × 10⁻⁵` |
| 8 | 27 | `3.13 × 10⁻³` | `6.83 × 10⁻⁵` |

The vacancy perturbation is much more local than the point load, so the
coupled model converges quickly; the residual `≈ 7 × 10⁻⁵` is the QNL
interface error at coarsening 4, which does not shrink with `r_atom`
because the interface simply moves to a region of similar (small) strain.

## 4. 2D triangular sheet

19×19 LJ triangular sheet, `a0 = 1.096966 σ`, `rmax = 1.9 a0`, all free
surface atoms clamped, a point force `(0.3, −0.15)` at the centre.
Reference: 361 atoms, `E = −1063.87029944`.

### Error vs atomistic-region radius (`coarsen = 3`)

| `r_atom / a0` | ndof | `max |Δu|` | `|ΔE|` |
|---|---|---|---|
| 2 | 105 | `2.78 × 10⁻³` | `3.15 × 10⁻³` |
| 3 | 129 | `2.78 × 10⁻³` | `3.56 × 10⁻³` |
| 4 | 153 | `2.18 × 10⁻³` | `1.35 × 10⁻³` |

### Error vs mesh coarsening (`r_atom = 3 a0`)

| coarsen | ndof | `max |Δu|` | `|ΔE|` |
|---|---|---|---|
| 2 | 160 | `2.78 × 10⁻³` | `2.24 × 10⁻³` |
| 3 | 129 | `2.78 × 10⁻³` | `3.56 × 10⁻³` |

The displacement error falls with a larger atomistic disc and the energy
error falls when the continuum mesh is refined. The total-energy error is
**not monotone** in `r_atom` at a fixed `coarsen`: shrinking the atomistic
disc leaves a differently shaped coarse mesh whose Cauchy-Born error
changes non-monotonically, and the residual `O(10⁻³)` is dominated by that
coarse-mesh error, not by the interface. This is an honest negative: for a
point load, growing the atomistic region alone is not an efficient error
reduction strategy. The next section shows the indicator does much better.

## 5. Adaptivity: reaching a target error with few degrees of freedom

The adaptive driver starts from a small atomistic disc (`1.5 a0`) and, at
each step, unfolds the three continuum cells with the largest Cauchy-Born
error `Σ_i |E_i^atom − E_i^CB(F_e)|`, then re-solves. Target energy error
`10⁻⁴`.

### Uniform refinement (grow the disc)

| `r_atom / a0` | ndof | `|ΔE|` |
|---|---|---|
| 2 | 105 | `3.15 × 10⁻³` |
| 3 | 129 | `3.56 × 10⁻³` |
| 4 | 153 | `1.35 × 10⁻³` |
| 5 | 169 | `1.65 × 10⁻³` |
| 6 | 225 | `9.77 × 10⁻⁴` |
| 8 | 283 | `2.06 × 10⁻³` |

Uniform refinement **never reaches** the `10⁻⁴` target within the domain
(the smallest error is `9.8 × 10⁻⁴` at 225 ndof).

### Adaptive unfolding

| step | ndof | unfolded cells | `|ΔE|` | `max |Δu|` |
|---|---|---|---|---|
| 0 | 105 | 0 | `3.15 × 10⁻³` | `2.78 × 10⁻³` |
| 1 | 129 | 3 | `7.07 × 10⁻⁴` | `2.57 × 10⁻³` |
| 2 | 143 | 6 | `3.69 × 10⁻³` | `2.78 × 10⁻³` |
| 3 | **151** | 7 | **`3.91 × 10⁻⁵`** | `2.12 × 10⁻³` |

The Cauchy-Born-error indicator reaches the `10⁻⁴` target at **151
degrees of freedom**, versus **≥ 283** for uniform refinement (which does
not reach the target). The adaptive model uses about **18 %** of the
361 fully atomistic degrees of freedom and about **53 %** of the best
uniform model. The decisive advantage comes from marking the cells where
the continuum *energy* is actually in error, not where the strain is
largest: the strain-jump indicator stalls at `≈ 2.3 × 10⁻³` because the
coarse-graining error accumulates over low-strain regions and is not
localised at strain concentrations.

## 6. Limits and honest misses

- **QNL is exact in 1D only in the sense of the patch test.** In 2D/3D the
  local rule used here is patch-test consistent, but the interface still
  carries an `O(a0)` modelling error that decays only as the interface
  moves into a lower-strain region. More sophisticated ghost-force
  corrections (Shenoy et al. 1999; the "ghost-force correction" of the
  force-based variant) are not implemented.
- **Total-energy error is non-monotone in the atomistic radius** at fixed
  coarsening (section 4). The coarse-mesh Cauchy-Born error is the
  dominant contribution for a point load, and it does not fall cleanly
  when the atomistic region grows.
- The continuum energy uses the **interpolated atom positions**, which for
  a linear element is the Cauchy-Born bond length, rather than the
  element-aggregated `V_e W(F_e)` form. Inter-element bonds are therefore
  treated with the local rule of their endpoint, not with a single `W(F)`;
  in the mesh-refinement limit the two agree. The naive coupling *does*
  use the element-aggregated Cauchy-Born energy, which is why it exhibits
  the ghost force.
- The implementation is **Bravais only** (no internal basis relaxation in
  the continuum), so non-Bravais crystals and angular/many-body potentials
  are out of scope; the Cauchy-Born `W(F)` for those would require the
  internal-basis Schur complement of `elastic.rs`.
- The indicator uses the **energy** difference (Cauchy-Born error). It is
  an empirical estimator, not a certified bound; the patch test and the
  comparison to the fully atomistic solution are the validation.
- The minimiser is Polak-Ribière conjugate gradients with backtracking.
  The 2D reference and coupled solves use a few thousand iterations and
  are the bulk of the run time; larger systems would need a proper sparse
  Newton/quasi-Newton solve.
- The triangular sheet uses a **structured nested coarse mesh** whose
  cells are unions of lattice cells. Unstructured adaptivity (arbitrary
  element sets, hanging nodes) is not implemented.

## 7. Constants and sources

- Quasicontinuum method: Tadmor, Ortiz & Phillips, *Phil. Mag. A* **73**,
  1529 (1996); review: Miller & Tadmor, *Modelling Simul. Mater. Sci.
  Eng.* **17**, 053001 (2009).
- Quasi-nonlocal coupling: Shimokawa, Mortensen, Schiøtz & Jacobsen,
  *Phys. Rev. B* **69**, 214102 (2004); analysis: E, Ming & Yang,
  *J. Comput. Phys.* (2006-2010).
- Ghost-force problem and force-based corrections: Shenoy, Miller,
  Tadmor, Rodney, Phillips & Ortiz, *J. Mech. Phys. Solids* **47**, 611
  (1999).
- Lennard-Jones 12-6 potential and the fcc lattice sums (`A6 = 14.45392`,
  `A12 = 12.13188`): Jones (1924); Lennard-Jones & Ingham (1925), as used
  in [materials.md](materials.md).
- Cauchy-Born rule for central-force crystals at zero pressure: Born &
  Huang, *Dynamical Theory of Crystal Lattices*.
- Conjugate-gradient minimiser (Polak-Ribière 1969) and Armijo
  backtracking (1966).
- The equilibrium spacings (`a0 = 1.119611 σ` for the 1D chain with
  `rmax = 2.5 a0`; `a0 = 1.096966 σ` for the triangular sheet with
  `rmax = 1.9 a0`) are computed in the tests by `relax_cutoff`, the raw
  cut-off zero-stress configuration of `potential.rs`, so the reference
  atomistic model is stress free at `F = I`.
