# Fourth CPU pass: periodic boundary accuracy

Periodic axes shorter than the ghost width now wrap every ghost layer back
to an interior cell. The correction applies to both fluid state and
transported scalar fields, including corners in multiple dimensions.
The [previous timestep and conduction improvements](cpu_round3.md) remain.
The next pass adds [viscous CPU improvements and transport timestep corrections](cpu_round5.md).

## Corrected behavior

The old periodic copy used source indices `n + gi` and `ng + gi`. Those are
interior indices only when `n >= ng`. On shorter axes, deep layers could
read uninitialized or stale ghosts. In a two-cell axis with six ghost
layers, a uniform Euler flow became non-finite in its first step with
positivity fallback disabled. The limiter could hide the bad stencil when
enabled, but correct boundary data should not depend on that fallback.

The new shared mapping reduces the physical ghost coordinate modulo the
interior axis length. A direct-copy branch keeps the existing source indices
when `n >= ng`. Inactive one-cell axes continue to have no allocated ghosts.
This repairs boundary data; it does not resolve waves below the grid scale
or raise the scheme's formal order.

55 tests pass across compressible flow, radiation and semiconductors.
Four new checks verify:

- Uniform flow on short periodic axes in all three orientations and a
  2×3×4 grid, using six and ten ghost layers, with no positivity fallback.
- Every fluid ghost, including corners, against an independent periodic
  coordinate mapping, starting from NaN-filled ghost storage.
- The same corner and deep-layer property for scalar fields.
- Three uniform scalars carried by the mixture and by each individual
  phase of a moving mixture: scalar mass conservation below 1e-13 relative
  error and specific-value deviations below 1e-12.

The existing scalar-advection, phase-carried scalar, scalar-diffusion, Sod
and acoustic-pulse examples retain identical printed accuracy metrics and
step counts. A 32³ evolved water/air render at 128×80, 64 spectral
paths/pixel, after 12 fluid steps retains byte-identical HDR output.
The [raw results](cpu_round4.json) include those comparisons and the HDR hash.

## CPU verification

The allocation, batching and FD-WENO cache experiments did not demonstrate
reliable gains across the tested workloads, so they were removed. The
existing flux kernels are retained. The CPU probe now also supports `weno`
and `fdweno`, providing separate WCNS and FD-WENO regression measurements.

Baseline runtime `2338d79`, retained candidate `bf21e43`: Intel Core
i7-13620H, eight Rayon threads, System76 Battery profile, CPU 10–50%, no
turbo. Release builds use thin LTO. Four workloads each ran four alternating
before/after pairs after compilation, with no overlapping heavy task.
Timings include the complete stepping loop and first-use buffers.

| Workload | Previous | Retained version | Previous/new |
|---|---:|---:|---:|
| Central Euler 32³, 100 steps | 1.487 s | 1.450 s | 1.03× |
| Hybrid water/air 32³, 40 steps | 1.514 s | 1.474 s | 1.03× |
| WCNS Euler 24³, 30 steps | 0.496 s | 0.431 s | 1.15× |
| FD-WENO Euler 24³, 30 steps | 0.600 s | 0.607 s | 0.99× |

All timing ranges overlap; this pass does not establish a new CPU speedup.
Every repeated before/after pair has the same conserved-state hash and
printed accuracy metrics on these normal-sized grids. The boundary fix
intentionally changes the previously invalid short-axis stencils.

## Reproduce

From the repository root, build both executables before timing:

```sh
git worktree add --detach /tmp/ns-cpu-baseline 2338d79
cargo build --manifest-path /tmp/ns-cpu-baseline/fluid/Cargo.toml -p compressible --release --example cpu_validate
cargo build --manifest-path fluid/Cargo.toml -p compressible --release --example cpu_validate
python3 fluid/scripts/cpu-bench.py /tmp/ns-cpu-baseline/fluid/target/release/examples/cpu_validate fluid/target/release/examples/cpu_validate --threads 8 --repeats 4 --case smooth:32:100 --case interface:32:40 --case weno:24:30 --case fdweno:24:30 --output /tmp/ns-cpu-round4.json
cargo test --manifest-path fluid/Cargo.toml -p compressible -p radiation -p semiconductor --release --lib --tests
git worktree remove /tmp/ns-cpu-baseline
```
