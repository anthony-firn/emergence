# heat

Phonon (lattice) heat transport and the **diffusive → ballistic** emergent
boundary. This is the phonon rung of the heat ladder: a continuum heat
equation with an effective conductivity, a cheap validity indicator, and a
kinetic unfold.

```text
  continuum heat equation  ρc ∂T/∂t = ∇·(κ(T,L) ∇T)   <-- coarse model
           ^   indicator: Kn = Λ/L  (phonon Knudsen number)
  kinetic (Callaway/Holland) κ(T), boundary scattering, ballistic Landauer
           ^   exact test case: grey BTE between two diffuse walls
  discrete ordinates / integral Boltzmann equation
```

## Modules

- `phonon`: per-branch Debye and Born–von Kármán sine dispersions, group
  velocities and the mode heat capacity.
- `scattering`: umklapp `τ_U⁻¹ ∝ γ² ω² T e^{−Θ/3T}`, normal, Klemens
  point-defect (isotope) `τ_I⁻¹ ∝ ω⁴` and Casimir boundary `τ_B⁻¹ = v/L`.
- `kinetic`: relaxation-time conductivity
  `κ = (1/3) Σ ∫ C(ω) v_g² τ g dω`, the Callaway normal-process correction,
  and the heat-capacity-weighted mean free path.
- `bte1d`: the grey Boltzmann equation between two diffuse walls, solved by
  discrete ordinates and by the exact integral equation; the Fourier and
  ballistic limits.
- `materials`: silicon (natural and ²⁸Si-enriched) and germanium parameters
  with sources.
- `effective`: `κ(T, L)` and the Fourier-validity flag for a continuum solver,
  plus the Landauer/Casimir ballistic limit.

## Example

```rust
use heat::{silicon_kappa, silicon_effective_kappa};
use heat::effective::TransportRegime;

// Bulk Si at 300 K (Morelli/Kazan fit; measured ≈148 W/(m K)).
let k = silicon_kappa(300.0);

// 50 nm nanowire: boundary scattering dominates.
let e = silicon_effective_kappa(300.0, 50e-9);
assert_eq!(e.regime, TransportRegime::Transitional);
assert!(e.kappa_effective < e.kappa_bulk);
```

```sh
cd fluid
cargo test -p heat
```

## Validation

15 tests compare against exact or published results: the silicon κ(300 K)
consistency check against Glassbrenner & Slack (1964), the low-temperature
peak, the isotope mass-fluctuation parameter against Inyushkin et al. (2018),
the isotope-enrichment effect, Si nanowire κ vs diameter against Li et al.
(2003), the exact quadrature and exponential integrals, the grey-BTE
discrete-ordinates solution against the independently derived integral
equation, the Fourier (Milne temperature-jump) and ballistic (Landauer)
limits, and the effective-conductivity regime flags. Measured errors are
tabulated in [`analysis/heat.md`](../../../analysis/heat.md). The suite runs
in about 30 seconds.

## The honest misses

The published Morelli/Kazan silicon parameters, used in the RTA as specified,
give κ(300 K) = 127 W/(m K) (the measured value is 148, −14 %). The Callaway
normal-drift correction overshoots to 218 W/(m K). The model also
under-predicts the isotope effect: +2.3 % at 300 K and a peak ratio of 4.2,
against the measured +8 % and ≈10. These are stated as consistency checks, not
predictions; the size effect and the BTE fold are the parts that are validated
independently.
