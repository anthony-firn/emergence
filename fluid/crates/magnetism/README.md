# magnetism

Magnetic materials across the emergence ladder: atomistic Heisenberg spins,
their fold to micromagnetics, and the collapse of a small particle to a
Stoner-Wohlfarth macrospin. This is the magnetic counterpart of the
materials ladder in [`materials`](../materials/README.md); the same
cheap-model / validity-indicator / unfold structure is used, and the
benchmarks are collected in [`analysis/magnetism.md`](../../../analysis/magnetism.md).

## Modules

- `geometry`: conventional sc/bcc/fcc cells, nearest-neighbour tables and
  the atomistic-to-continuum exchange-stiffness fold
  `A = J S^2 sum_delta |delta|^2 / (12 V_atom)`.
- `spin`: classical Heisenberg Hamiltonian, Metropolis Monte Carlo on a
  finite periodic cluster, the Binder cumulant and Curie temperature,
  linear spin-wave theory, and a spin-chain Landau-Lifshitz-Gilbert (LLG)
  integrator.
- `micromag`: the 1D finite-difference micromagnetic LLG solver (exchange,
  uniaxial anisotropy, Zeeman, local-shape demagnetising field), the exact
  Bloch wall, the reduced 1D Walker model, Walker breakdown and the Kittel
  thin-film resonance.
- `hysteresis`: the Stoner-Wohlfarth astroid, the LLG macrospin, the
  random-orientation assembly coercivity, and the single-domain indicator.
- `demag`: the Newell-Williams-Dunlop cuboid demagnetising tensor and a
  zero-padded FFT demagnetising-field operator (the non-local unfold of the
  local-shape closure).
- `llg3d`: a 3-D finite-difference LLG solver built on the FFT demag field.
- `sp4`: the NIST µMAG standard problem #4 reversal as a reusable run.
- `materials`: sourced room-temperature constants for iron, nickel,
  permalloy and NdFeB.
- `rng`: a small deterministic xorshift generator, so the Monte-Carlo runs
  are reproducible without external crates.

## Example

```rust
use magnetism::geometry::{exchange_stiffness, LatticeKind};
use magnetism::micromag::{wall_delta, wall_sigma, Llg1d};

// Fold the atomistic exchange to the micromagnetic stiffness.
let (j, s, a) = (1.0, 0.5, 0.2866e-9);
let a_ex = exchange_stiffness(LatticeKind::Bcc, j, s, a);

// Exact Bloch wall: delta = sqrt(A/K), sigma = 4 sqrt(A K).
let k = 4.8e4;
let _ = (wall_delta(a_ex, k), wall_sigma(a_ex, k));

// 1D LLG relaxation of the same wall on a finite-difference grid.
let mut g = Llg1d::new(200, wall_delta(a_ex, k) / 8.0, a_ex, k, 1.71e6, 1.0);
g.seed_wall();
g.relax(2.0e-14, 40_000);
```

## Benchmarks (numbers in `analysis/magnetism.md`)

- classical sc Heisenberg Monte-Carlo Curie temperature: Binder crossing
  `1.43 ± 0.04` versus the published `1.44299` (Chen, Ferrenberg & Landau
  1993); mean field gives `2.0`.
- spin-wave dispersion `omega = 2 J S (z - sum cos k.delta)` agrees with
  the LLG spin chain to `2.5e-5`.
- discrete domain-wall fold: continuum `sigma = 4 sqrt(A K)` and
  `delta = sqrt(A/K)` reproduced to `0.04 %` at `delta/a = 10`, with the
  discrete correction growing to `4.5 %` at `delta/a = 1`.
- micromagnetic Bloch wall: `sigma` to `6.5e-4`, `delta` to `1.3e-3`.
- Walker breakdown: `H_W = alpha M_s/2` and `v = gamma delta H/alpha` are
  exact in the reduced 1D model (below `H_W`), with the average velocity
  falling to `0.26` and `0.14` of the linear law at `1.5 H_W` and `2 H_W`.
- Stoner-Wohlfarth astroid reproduced by the LLG macrospin to <4 %; the
  random-orientation coercivity is `0.492` versus the exact `0.479 H_K`.
- FFT demagnetising field (Newell-Williams-Dunlop): uniform cube `H_x =
  -M_s/3` and energy `(mu0/2) M_s^2 V/3` to machine precision; the
  Lorentz-sphere sum is `-M_s/3` at every radius; thin-film `N_z` tends to
  1; the full and local demag agree for a uniformly magnetised thin film
  and rod but not for a vortex disc.
- NIST µMAG standard problem 4 (coarse): first `<m_x>` / `<m_y>` zero
  crossings within `2.5 %` / `6.7 %` of the downloaded OOMMF reference.

## Examples

- `mumag_sp4`: the standard-problem-4 reversal (`cargo run --release
  --example mumag_sp4 -- [nx ny nz]`, default `64 16 1`).
