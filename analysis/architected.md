# Architected mechanical metamaterials: beam lattices folded to a continuum

Implemented in `fluid/crates/architected` (new crate; no dependencies). This
is the mechanical instance of the ladder in
[`emergence_ladder.md`](emergence_ladder.md). The fine model is a linear frame
(pin-jointed bar or rigid-jointed Euler-Bernoulli/Timoshenko beam) lattice of
struts; the coarse model is an effective elastic stiffness `C*`; the fold is
periodic homogenisation of one unit cell. The two validity indicators are the
**sample size in cells** (boundary layer and shear) and the **strut stress
against Euler buckling**.

```sh
cd fluid
cargo test -p architected          # all tests, each well under 30 s
```

## The models

* `beam`: the 6-DOF 2D frame element (axial + bending, Timoshenko via
  `Phi = 12 E I/(G A_s L^2)`) and the 12-DOF 3D frame element (axial, bending,
  torsion `G J/L`), plus the 4- and 6-DOF pin-jointed bar elements. The 2D and
  3D cantilever/torsion elements are exact: a single element reproduces
  `P L^3/(3 E I)`, `P L^3/(192 E I)`, `M L/(G J)` to `10^{-13}`.
* `cell`: periodic unit cells. Every member endpoint is reduced modulo the
  lattice vectors to a canonical node plus a lattice image `d`; images receive
  the same periodic fluctuation, and the affine displacement `eps.x` is applied
  at the physical position. The fluctuation is the energy minimiser, and the
  effective stiffness is read off the strain energy,
  `C_ab = [U(e_a+e_b) - U(e_a) - U(e_b)]/V`, which needs no boundary tractions
  (engineering Voigt strains; `g_ij = 2 e_ij`).
* `lattices`: regular and re-entrant hexagonal honeycombs of any wall angle,
  triangular trusses, the octet truss (FCC nearest-neighbour network), a
  simple-cubic frame and a square frame. Overlapping blocks are generated and
  duplicate bonds are removed automatically, so the generators are written
  naively.
* `frame`: a small direct 2D frame solver for the benchmark beams.
* `fold`: a finite honeycomb cantilever compared with the homogenised beam.
* `buckling`: the Euler load of the most compressed strut under a macroscopic
  strain, giving the strain at which the linear fold stops holding.
* `bloch`: elastic-wave bands of a 2D frame cell by Bloch-Floquet. `K(k)` is
  Hermitian; writing `K = A + iB`, the real symmetric representation
  `[[A,-B],[B,A]] phi = omega^2 [[M,0],[0,M]] phi` is diagonalised.

## Benchmarks (measured)

Every entry is printed by the tests in `tests/benchmarks.rs`,
`tests/fold.rs` and `tests/dispersion.rs`. "Error" is the deviation from the
exact or published value, never a tuned tolerance.

| Benchmark | Exact / published | Computed | Error |
|---|---|---|---|
| 2D cantilever tip, end load | `P L^3/(3 E I)` | one frame element | `<1e-12` |
| 2D fixed-fixed centre, point load | `P L^3/(192 E I)` | two frame elements | `<1e-12` |
| 3D cantilever, bending about y/z | `P L^3/(3 E I)` | one 3D frame element | `<1e-12` |
| 3D torsion | `M L/(G J)` | one 3D frame element | `<1e-12` |
| 2D Timoshenko cantilever, stubby | `P L^3/(3 E I) + P L/(G A_s)` | one frame element, `L/h=3.3` | `<1e-12` |
| Regular hexagon, bending only `t/l=0.005` | `(4/sqrt3)(t/l)^3 = 2.8868e-7` | `2.88654e-7` | `9e-5` (0.009%) |
| Regular hexagon, bending only `t/l=0.01` | `2.3094e-6` | `2.30871e-6` | `3e-7` (0.03%) |
| Regular hexagon, bending only `t/l=0.04` | `1.47802e-4` | `1.47096e-4` | `7.1e-7` (0.48%) |
| Regular hexagon, Poisson `t/l -> 0` | `1` | `0.99990` at `t/l=0.005` | `1e-4` |
| Regular hexagon, stretching (Masters & Evans) | see below | frame cell | `<1e-13` |
| General hexagon bending `E1/Es` (3 cells) | `(t/l)^3 cosTh/((h/l+sinTh) sin^2 Th)` | frame cell | `<0.1%` at `t/l=0.005` |
| General hexagon bending `nu12` | `cos^2Th/((h/l+sinTh) sinTh)` | frame cell | `<1e-3` absolute |
| Re-entrant `h/l=2`, `Th=-30 deg`, bending | `nu12 = -1` | `-0.99960` at `t/l=0.01` | `4e-4` (`O((t/l)^2)`) |
| Triangular truss | `E*/Es = rho/3`, `nu=1/3`, isotropic | exact | `<1e-15` |
| Octet truss (FCC truss), `t/l=0.02..0.05` | `E*/Es = rho/9` (Deshpande 2001) | exact | `<1e-9` |
| Octet truss symmetry | cubic `C11 = 2 C12 = 2 C44`, `nu=1/3` | exact | `<1e-9` |
| Simple-cubic frame shear | `G*/Es = rho^2/18` | exact, slope 2.000 | `<1e-9` |
| Honeycomb buckling strain | `eps_cr ~ (t/l)^2` | slope 2.00, `eps_cr = 1.644 (t/l)^2` | `<0.2%` in coefficient |
| Square-frame acoustic slopes | `sqrt(C11/rho*)`, `sqrt(C44/rho*)` | Bloch `k -> 0` | `<2e-3` |
| Honeycomb cantilever vs continuum | ratio `-> 1` as cells grow | see fold below | boundary layer ~1 cell |

### Regular honeycomb with stretching

The unit-cell beam analysis (Gibson & Ashby's one-strut cell with axial
extension retained; the stretching-inclusive model of Masters & Evans 1996)
gives, for the regular hexagon (`h = l`, `Th = 30 deg`),

\[
  \frac{E^*}{E_s} = \frac{4}{\sqrt3}\,\frac{(t/l)^3}{1 + 3 (t/l)^2},
  \qquad
  \nu^* = \frac{1 - (t/l)^2}{1 + 3 (t/l)^2}.
\]

The frame unit cell reproduces both to machine precision (`<1e-13`) for
`t/l = 0.05 .. 0.25`, so this is an exact check of the element, the periodic
reduction and the energy homogenisation at once. At `t/l=0.2`,
`E*/Es = 1.6496e-2` and `nu = 0.857143 = 6/7`, both matching the closed form.
The bending-only Gibson-Ashby value `2.309 (t/l)^3` is the `t/l -> 0` limit
and is high by 10.7% at `t/l=0.2`, as expected.

The octet truss deserves a caveat: the FCC nearest-neighbour network is
**not** elastically isotropic. Its elastic constants satisfy
`C11 = 2 C12 = 2 C44` (Zener factor 2), i.e. `nu = 1/3` along `<100>` but
direction-dependent. The published `E*/Es = (1/9) rho` is the `<100>` Young's
modulus and is reproduced exactly.

## Fold / unfold map

The homogenised `C*` is exact only when the deformation wavelength is large
compared with the cell and when no strut has buckled. The two boundaries:

| coarse model | indicator | value / threshold | unfold |
|---|---|---|---|
| homogenised `C*` | sample length `N` in cells | relative error `~ 1.13/N` (one-cell clamp boundary layer) | discrete frame lattice |
| homogenised `C*` | `h_cell / lambda_deformation` | small `k h_cell`; Bloch slope error `O((k h)^2)` | Bloch bands per cell |
| linear `C*` | strut stress / Euler load | `eps_cr = P_euler / N(eps=1)` | post-buckling / nonlinear strut |
| effective beam `E* I` | shear `E*/(G*) (H/L)^2` | grows for short, thick beams | Timoshenko beam |

**Size effect.** For a honeycomb cantilever with fixed aspect ratio (`L/H ~ 4`)
and growing cell count `N`, the tip deflection relative to the
Euler-Bernoulli continuum `P L^3/(3 E^* I)` is

| `N` (cells along length) | 4 | 6 | 8 | 12 | 16 |
|---|---|---|---|---|---|
| `delta_FEM / delta_cont` | 0.794 | 0.838 | 0.868 | 0.906 | 0.930 |
| `(1 - ratio) N` | 0.82 | 0.97 | 1.06 | 1.13 | 1.13 |

The ratio rises monotonically to one; the residual is a clamped-end boundary
layer of about one cell, `(1 - ratio) ~ 1.13/N`. The size-effect indicator is
therefore `h_cell / L`: the fold is certified while it is below the tolerance
one needs, and the discrete lattice is unfolded otherwise.

**Buckling.** Under uniaxial compression the axial force in each strut is
linear in the macroscopic strain, so the first strut buckles at
`eps_cr = P_euler/|N(eps=1)|`. For the regular honeycomb
`eps_cr = 1.644 (t/l)^2` (Euler, pinned ends) — a quadratic, not cubic,
scaling, so buckling overtakes the linear modulus as the limiting process at
moderate `t/l`. The homogenised `C*` is valid only below this strain.

**Dispersion.** The `k -> 0` slopes of the square-frame acoustic branches are
`sqrt(C11/rho*)` (longitudinal) and `sqrt(C44/rho*)` (shear) to `2e-3`,
confirming that the same `C*` controls both the static fold and the
long-wavelength dynamics.

## Limitations / honest notes

* Only linear elasticity and small strains. No geometric nonlinearity, joint
  compliance, or post-buckling path; the buckling result is the *onset*
  boundary, not the collapse load.
* Homogenisation uses one unit cell and the energy method; it gives the
  effective stiffness, not a rigorous error bound on `C*` for a finite body.
  The fold test measures that error empirically.
* The 2D Timoshenko option exists and is exact for the end-load cantilever
  (validated above); the 3D frame uses `G = E/(2(1+nu))` for torsion. The
  benchmark lattices are all Euler-Bernoulli (`nu` default 0.3 only enters 3D
  torsion).
* The octet truss is cubic, not isotropic; only its `<100>` modulus matches
  the published `rho/9`.
* Dense linear algebra is used throughout (fine for unit cells and the fold
  blocks, which are up to a few thousand DOF).
