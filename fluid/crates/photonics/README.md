# photonics

Wave optics of layered and structured materials, the rung below
[`radiation`](../radiation/README.md). Geometrical optics is what remains when
the phase of multiply-reflected partial waves and the sub-wavelength structure
of the matter are both unresolved; this crate keeps both and provides the
indicator and the unfold/fold for each.

## Modules

- `complex`: small `C` complex type; the principal `sqrt` branch gives decay
  for absorbing media and evanescent waves.
- `multilayer`: transfer-matrix method for planar isotropic stacks, any angle,
  TE/TM, complex `N = n + i kappa`. Returns `r, t, R, T, A`, per-layer
  absorptance, the internal field, and the absorbed power density (W m⁻³).
- `bloch`: exact 1D photonic-crystal bands from one unit cell,
  `cos(K Lambda) = (1/2) tr M`; stop-band decay length.
- `dispersion`: causal Sellmeier / Drude / Lorentz models. Sources: fused
  silica (Malitson 1965), water (Daimon–Masumura 2007), gold and silver
  (Rakić et al. 1998).
- `effective_medium`: Maxwell Garnett, Bruggeman, and the exact long-wavelength
  laminate limit (`eps_par`, `1/eps_perp`); anisotropic-slab comparison and a
  `Lambda/lambda` error indicator.
- `incoherent`: the geometrical-optics intensity sum and the partial-coherence
  crossover from coherent Airy to incoherent (indicator `n d / L_c`).
- `plasmon`: optional Kretschmann surface-plasmon reflectance and dip search.

## Example

```rust
use photonics::{Layer, Polarization, Stack};

// Quarter-wave antireflection coating: n_c = sqrt(n_s).
let n_s = 1.5;
let stack = Stack::single(1.0, Layer::quarter_wave(n_s.sqrt(), 550e-9), n_s);
let sol = stack.solve(550e-9, 0.0, Polarization::TE);
assert!(sol.reflectance < 1e-20);
// Heat sink for a thermal solver on a lossy film:
let lossy = Stack::single(1.0, Layer::new(0.18, 3.0, 80e-9), 1.5);
let s = lossy.solve(633e-9, 0.0, Polarization::TE);
let profile = s.power_density_profile(1000.0, 8); // (z [m], W m^-3)
```

```sh
cd fluid
cargo test -p photonics
```

## Validation

17 tests compare against exact or published results: Fresnel, Brewster, total
internal reflection, quarter-wave ARC, the quarter-wave Bragg mirror formula,
Bloch band edges, energy conservation, reciprocity, the Rakić gold/silver
fits, a numerical Kramers–Kronig check, effective-medium convergence, the
coherence crossover, the heat-source energy balance and the Kretschmann dip.
Measured errors are tabulated in [`analysis/photonics.md`](../../../analysis/photonics.md).
Tests run in well under a second.

## Limits

Planar, isotropic, static layers only. The stop-band formula is the
leading-order one (the exact gap is 2.1% wider at `n_H/n_L = 1.586`). The
effective-medium error indicator is an empirical quadratic envelope, not a
proven bound. The incoherent partial-coherence sum assumes lossless
interfaces. No roughness, anisotropy, nonlinear optics, waveguiding or 2D/3D
photonic crystals yet.
