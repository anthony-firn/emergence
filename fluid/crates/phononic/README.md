# phononic

Elastic and acoustic waves in layered and structured materials, and their
homogenised (folded) descriptions. This is the mechanical counterpart of
[`photonics`](../photonics/README.md), the rung below a ray/geometrical
description: it keeps the phase of multiply reflected sound and the
sub-wavelength structure of the matter, and provides the cheap indicator and
the fold/unfold for each.

## Modules

- `complex`: small `C` complex type, principal `sqrt` (decay for a passive
  complex modulus), and a rank-safe small complex linear solver.
- `materials`: sourced handbook properties of water, steel, aluminium and PMMA
  (density, P and S speeds, derived Lamé/bulk/Young moduli).
- `layered`: transfer matrices for 1D elastic/acoustic stacks at normal
  incidence, loss through a complex modulus, `r, t, R, T, A`, and the exact
  Rytov Bloch dispersion `cos(qΛ) = cos(k1d1)cos(k2d2) − ½(Z1/Z2+Z2/Z1)
  sin(k1d1)sin(k2d2)` of the infinite periodic medium.
- `homogenize`: the exact long-wavelength laminate modulus (Reuss) and density,
  Voronoi-symmetric Voigt/Reuss/Hashin–Shtrikman bounds for isotropic two-phase
  composites, and the phase-velocity validity indicator.
- `resonant`: the mass-in-mass locally resonant lattice, its exact and
  effective-mass dispersions, the negative-mass interval and the sub-wavelength
  band gap.
- `fluid_coupling`: oblique fluid–solid–fluid transmission through an elastic
  plate with P and SV waves, showing coincidence and Lamb-wave peaks.

## Example

```rust
use phononic::layered::{Layer, Stack};
use phononic::materials::{STEEL, WATER};
use std::f64::consts::PI;

// Steel plate in water: transparent at the half-wave thickness resonances.
let d = 5.0e-3;
let f = STEEL.c_p / (2.0 * d);           // first resonance
let plate = Layer::from_solid(&STEEL, d);
let s = Stack::from_media(&WATER, vec![plate], &WATER).solve(2.0 * PI * f);
assert!((s.transmittance - 1.0).abs() < 1e-8);
```

```sh
cd fluid
cargo test -p phononic
```

## Validation

17 tests compare against exact or published results: half-wave transmission
resonances, the lossless energy identity `|r|² + (Z_out/Z_in)|t_v|² = 1`, the
analytic quarter-wave band edges and their agreement with the exact Rytov gap,
the Bloch gap against the finite-stack transmission dip, quadratic convergence
of the Bloch phase velocity to the Reuss speed, the Hashin–Shtrikman ordering,
the mass-in-mass effective mass against the exact 2-DOF eigenproblem, the
sub-wavelength gap, the fluid–plate normal-incidence limit, the coincidence
peaks and energy conservation at oblique incidence. Measured errors are
tabulated in [`analysis/phononic.md`](../../../analysis/phononic.md). The suite
runs in under two seconds.

## Limits

- Isotropic, planar, perfectly bonded layers; no anisotropy, roughness or
  interface slip.
- The fluid–plate transfer matrix mixes growing and decaying evanescent
  solutions, so it loses precision past the shear critical angle when
  `|Im k_zS| d` is large; a scattering-matrix reformulation is the unfold.
- Hashin–Shtrikman bounds are for isotropic two-phase composites without
  percolation or interfacial resistance.
- The locally resonant model is a discrete mass-spring chain, not a continuum
  and not viscoelastic.
- The homogenisation indicator is an empirical quadratic envelope, a switch
  rather than a proven bound.
