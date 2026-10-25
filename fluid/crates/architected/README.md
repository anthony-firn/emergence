# Architected: beam lattices and their effective continuum

Mechanical metamaterials in the ladder of
[`analysis/emergence_ladder.md`](../../../analysis/emergence_ladder.md): a
linear frame (or pin-jointed bar) lattice is *folded* to an effective elastic
stiffness `C*` by periodic homogenisation of one unit cell, and *unfolded* to
the struts when the size-effect or buckling indicator fails.

- `beam`: 2D and 3D linear frame and bar element stiffness matrices. Exact
  single-element cantilever, fixed-fixed and torsion checks.
- `cell`: periodic unit cells, periodic boundary conditions, macroscopic strain
  loading and effective stiffness by the energy method (frame or truss).
- `lattices`: regular and re-entrant honeycombs, triangular, octet truss,
  simple-cubic and square beam lattices.
- `frame`: a small direct 2D frame solver for benchmark beams.
- `fold`: a finite honeycomb cantilever against the homogenised beam; the
  size-effect indicator and the clamped-end boundary layer.
- `buckling`: the Euler load of the most compressed strut under a macroscopic
  strain (the instability boundary of the linear fold).
- `bloch`: Bloch-Floquet elastic-wave bands of a 2D frame cell; long-wavelength
  slopes equal the homogenised sound speeds.
- `linalg`: the small dense solver (LU, inverse, Cholesky, symmetric and
  generalised eigenproblems) used throughout. No external dependencies.

```sh
cd fluid
cargo test -p architected
```

Benchmarks and the fold/unfold map are in
[`analysis/architected.md`](../../../analysis/architected.md).
