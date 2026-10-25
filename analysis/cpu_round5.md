# Fifth CPU pass: viscous flux reuse and transport timestep corrections

The viscous RHS is 1.71–1.86× faster on the measured 32³ and 48³ grids,
with bit-identical residuals and fixed-step conserved states. The timestep
estimate now uses temperature-dependent transport and the configured
thermodynamics, including fluxes from current boundary states.
The [raw runs and validation outputs](cpu_round5.json) accompany this report.
The next pass adds [conduction and reflecting-boundary improvements](cpu_round6.md).

## Kernel changes

Each gradient is computed once, used to form compact momentum/energy fluxes,
and discarded. Those fluxes replace the full gradient buffer; each flux is
reused by neighboring divergence stencils. The primitive and flux buffers
occupy the same padded-grid storage as the old primitive and gradient buffers.
The interior-index vector and residual staging vector are removed: disjoint
output rows accumulate directly in parallel. For Euler, this removes another
48 bytes per interior cell, or 1.5 MiB on a 32³ grid, per RHS evaluation.
The sixth-order derivative and its per-component arithmetic order are retained.

## Corrected timestep estimates

Previously, explicit and all-speed CFL queries used reference viscosity,
minimum density, and a hardcoded heat-capacity ratio of 1.4. They ignored
`Transport::power`, `Transport::bulk`, and the EOS heat capacities.
A hot region or a large bulk viscosity could therefore exceed the assumed
diffusion coefficient. Using only the largest local μ/ρ also misses a hot
flux node driving diffusion into a lighter neighboring cell.

The shared estimator now collects maximum μ and κ = μ cp/Pr, minimum
interior density, and minimum interior volumetric heat capacity ρ cv.
It applies

```
nu = max((4 mu_max/3 + bulk)/rho_min, 4 kappa_max/(3 capacity_min))
dt_viscous = 0.3 / (nu * sum_active_axes(1/dx²))
```

This retains the previous thermal safety factor. `System::cv` defaults to
the phase-derived heat capacity; Euler supplies its EOS value without a
phase-vector allocation. Dirichlet states and all three required
characteristic ghost flux layers contribute transport coefficients.
These states are reconstructed from current interiors, including before
the first step and with NaN-filled stored ghosts. Ghosts have no evolved
heat capacity, so capacity minima come from interiors.

The estimate freezes coefficients at the current state. Coupled nonlinear
heating and rapidly changing transport can require smaller steps.

## Measured CPU performance

Baseline `3cfd80e`, final runtime `dcef97a`; Intel Core i7-13620H,
eight Rayon threads, System76 Battery profile, CPU 10–50%, no turbo.
Release builds use thin LTO. Four alternating before/after pairs per case,
with compilation completed and no concurrent heavy task. Times include
first-use allocation and the complete measured loop.

| Workload | Previous median | New median | Previous/new |
|---|---:|---:|---:|
| Viscous RHS 32³, 100 evaluations | 0.395 s | 0.230 s | 1.71× |
| Viscous RHS 48³, 60 evaluations | 0.717 s | 0.386 s | 1.86× |
| Complete viscous flow 24³, 60 steps | 0.924 s | 0.765 s | 1.21× |
| Viscous RHS 256×1×1, 1,000 evaluations | 0.079 s | 0.072 s | 1.10× |
| Viscous RHS 4096×1×1, 200 evaluations | 0.066 s | 0.037 s | 1.77× |
| Complete viscous flow 32³, 100 steps | 2.564 s | 2.291 s | 1.12× |

Timing ranges are disjoint for the 32³/48³ RHS cases and the 4096-cell
line. They overlap for the 256-cell line and both complete-step cases;
the complete-step median gains are 12–21%, with visible run-to-run variance.
Every pair has the same result hash and printed conservation metric.
The probes use temperature-dependent viscosity and nonzero bulk viscosity.
Full stepping uses fixed dt = 1e-4 to measure the operator independently
of the deliberately corrected CFL rule.

The 32³ mixture and all-speed CFL probes retain identical timestep/state
hashes and overlapping timing ranges. The constant-viscosity Euler CFL
probe also retains an identical timestep hash, but its richer coefficient
scan increases median query time from 0.187 s to 0.209 s for 1,000 queries
(about 12%). The measured cost is approximately 22 μs per query.

## Accuracy and regressions

63 tests pass across compressible flow, radiation, and semiconductors.
Eight new viscous tests check:

- Sixth-order convergence of shear momentum and viscous energy work on
  every axis, with periodic momentum/energy balance.
- Longitudinal bulk stress and thermal conductivity against analytic
  derivatives, including gamma=2.3 and cv=7.
- Sixth-order convergence with temperature-dependent viscosity.
- Explicit and all-speed diffusion bounds for hot states, bulk viscosity,
  EOS heat capacities, hot dense cells next to light cells, and prescribed
  or characteristic boundary states in every orientation.
- Identical residual bits across one and three Rayon threads on short,
  two-dimensional, and three-dimensional meshes; untouched output ghosts.

The same test file on the previous runtime has four failures, all in the
new timestep regressions. Its four spatial/concurrency tests pass.
Central Euler, WCNS, FD-WENO, hybrid water/air, and characteristic open-flow
24³ probes retain identical accuracy metrics and conserved-state hashes.
The hybrid Becker viscous-shock validation at Mach 2 and μ=0.02,
with 50, 100, 200, and 400 cells, retains identical printed accuracy
metrics and WENO fractions.

## Reproduce

From the repository root, build before timing:

```sh
git worktree add --detach /tmp/ns-viscous-baseline 3cfd80e
cargo build --manifest-path /tmp/ns-viscous-baseline/fluid/Cargo.toml -p compressible --release --example cpu_viscous --example cpu_timestep --example cpu_validate --example becker
cargo build --manifest-path fluid/Cargo.toml -p compressible --release --example cpu_viscous --example cpu_timestep --example cpu_validate --example becker
python3 fluid/scripts/cpu-bench.py /tmp/ns-viscous-baseline/fluid/target/release/examples/cpu_viscous fluid/target/release/examples/cpu_viscous --threads 8 --repeats 4 --case rhs:32:100 --case rhs:48:60 --case step:24:60 --case step:32:100 --case line:256:1000 --case line:4096:200 --output /tmp/ns-viscous.json
python3 fluid/scripts/cpu-bench.py /tmp/ns-viscous-baseline/fluid/target/release/examples/cpu_timestep fluid/target/release/examples/cpu_timestep --threads 8 --repeats 4 --case mixture:32:1000 --case allspeed:32:1000 --case viscous:32:1000 --output /tmp/ns-viscous-cfl.json
env -u DEBUG RAYON_NUM_THREADS=8 MODE=hybrid /tmp/ns-viscous-baseline/fluid/target/release/examples/becker
env -u DEBUG RAYON_NUM_THREADS=8 MODE=hybrid fluid/target/release/examples/becker
cargo test --manifest-path fluid/Cargo.toml -p compressible -p radiation -p semiconductor --release --lib --tests
git worktree remove /tmp/ns-viscous-baseline
```
