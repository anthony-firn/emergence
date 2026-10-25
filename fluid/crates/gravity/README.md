# Isolated gravity and static spacetime

All source coordinates, softening lengths, masses and clocks use SI. `Sample`
returns gravitational potential in m²/s², with zero at infinity, and acceleration
in m/s². `G=6.67430e-11` is pinned to CODATA 2022; `c=299792458` is exact.

## Source representation

`Spherical::new(edges, densities, center)` integrates concentric constant-density
shells exactly, including their interior potential. Construction takes O(n)
and queries O(log n). It has no artificial central singularity. Use this for a
spherical stellar background rather than replacing the interior with a point mass.

`Direct::new(bodies, softening)` validates sources once and sums the Plummer law
`phi=-Gm/sqrt(r²+epsilon²)` exactly up to floating-point reduction. `Tree::new`
uses an isolated octree with softened monopole and quadrupole moments. Its
acceptance test bounds the source support relative to the distance from the
center of mass and opens nodes containing either the target or an excluded
source. `theta=0` opens every node; permitted theta is 0..=0.7. `evaluate(x,
Some(i))` excludes precisely source i. `Work` reports traversal and direct pair
counts. `at_bodies()` evaluates histories in parallel from 1024 bodies.

Softening represents unresolved source size; choose it and theta by spatial
convergence and direct-force comparisons. The tree has no universal error
certificate. It is not a periodic Poisson solver. The implementation follows
the multipole approach described in [Springel's GADGET-2 paper](https://wwwmpa.mpa-garching.mpg.de/~volker/gadget/gadget2-paper.pdf),
with its own geometric acceptance criterion. Timing example:

```sh
cd fluid
RAYON_NUM_THREADS=8 cargo run -p gravity --release --example benchmark
cargo test -p gravity --release
```

`NBody` provides kick-drift-kick leapfrog, rebuilding the field at drifted
positions. `invariants()` independently sums direct mechanical energy and
momentum, at O(n²) cost. Approximate tree forces do not enforce pairwise action
and reaction; monitor drift and refine theta when necessary. A suitable time
step still needs orbital/encounter convergence. This is Newtonian bulk matter,
not a relativistic N-body Einstein solver.

`compressible::gravity::{sources,self_field,kick}` uses cell mass `rho*dV` and
applies a frozen-field impulse with exact kinetic work, preserving internal
energy. It prepares all updates before committing. The caller owns temporal
splitting, field rebuilds and boundary conditions. This source alone is not a
well-balanced hydrostatic solver. Conserving self-gravitating total energy
through mass advection requires consistent gravitational potential work in the
hydrodynamic flux update; that closure is not implemented here.

## Relativistic particles

`metric::Metric` supplies static isotropic geometry
`ds²=-alpha² c² dt²+beta² dx²`. `Minkowski`, `WeakField<F>` and `Schwarzschild`
are implemented. WeakField uses `alpha²=1+2phi/c²`, `beta²=1-2phi/c²` and rejects
`abs(phi)/c²>1e-3`. It retains the leading gravitational redshift and spatial
curvature; it is not a complete first-post-Newtonian stellar solution.
Schwarzschild is the exact spherical vacuum metric in isotropic coordinates,
outside the horizon. Its coordinate radius differs from areal radius:
`R=r*(1+GM/(2*c²*r))²`. It cannot describe a material stellar interior.

The [Standard Model particle API](../standard-model/README.md) integrates the
Hamiltonian `H=alpha*sqrt(m²+P²/beta²)` with canonical `P=beta*p_local` and
proper-time decays. Photons and neutrinos follow null geodesics, and massive
species follow timelike geodesics. Tests include solar-limb light bending,
redshift and the stationary Killing energy. Geometry is prescribed and frozen;
source backreaction, dynamical spacetime and simultaneous Lorentz-force plus
geodesic integration need additional coupling.
