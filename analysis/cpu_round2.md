# Second CPU improvement pass

The [third pass](cpu_round3.md) adds faster timestep queries and corrects
all-speed boundary and heat-conduction stability estimates.

This pass improves the fluid solver on top of the previous CPU optimizations.
It reuses two-point split fluxes between adjacent faces, reuses cell states
and physical fluxes inside the positivity limiter, prepares characteristic
boundaries once per pencil rather than per ghost layer, and parallelizes
independent RK cell updates on larger multidimensional grids.

It also repairs two accuracy/diagnostic issues: the positivity limiter now
applies its blend to the primitive face states used by general
non-conservative products (the GPR solid/viscoelastic model), and the solver
retains the number of limited faces when collecting pencil statistics.
Previously it could report zero even when limiting was active.

## Timings against the previous version

Baseline `572c566` contains the previous solver plus the new open-boundary
probe; candidate `34ce425` contains the changes. Intel Core i7-13620H,
eight Rayon threads, Battery profile, CPU limited to 10–50%, no turbo.
Release builds use thin LTO. Six workloads each ran four alternating pairs;
no build or other heavy task overlapped the timed runs. The table reports
median wall time inside the stepping loop, including its first-use buffers.
[All measurements, ranges and state hashes](cpu_round2.json) are retained.

| Workload | Previous | New | Previous/new |
|---|---:|---:|---:|
| smooth 32³, 100 steps | 1.562 s | 1.464 s | 1.07× |
| smooth 48³, 60 steps | 2.578 s | 2.229 s | 1.16× |
| interface 24³, 60 steps | 1.343 s | 1.254 s | 1.07× |
| interface 32³, 40 steps | 1.574 s | 1.485 s | 1.06× |
| open 24³, 60 steps | 1.193 s | 0.672 s | 1.77× |
| open 32³, 40 steps | 1.257 s | 0.875 s | 1.44× |

Every repeated before/after pair has the same interior conserved-state hash
and printed accuracy metrics. Gains are workload-specific; short timings
retain CPU scheduling/frequency noise. The separated ranges in the two
open-boundary cases support a clearer gain than the smaller periodic gains.
These improvements are incremental to the prior pass; they are not a new
renderer speedup. The renderer's sampling and radiance estimator are unchanged.

## Work saved without changing the stage arithmetic

The sixth-order central flux combines six two-point evaluations per face.
A lag-L pair contributes to L neighboring faces. A small tagged stack ring
reuses lag-2 and lag-3 pairs, reducing a fully central run to approximately
three new pair evaluations per face. Tags handle gaps where the hybrid
scheme selects its shock path. Fluxes accumulate in their original order.
The unused WCNS midpoint array is no longer allocated in pure central or
FD-WENO mode.

The positivity limiter carries the preceding cell's conserved state and
physical flux into the next face, instead of evaluating each twice. Its
admissibility tests and 30 bisection iterations are unchanged. Limited
non-conservative face states now use the same blend and first-order mean
partner as the flux and normal velocity; failed non-finite reconstructions
recover that mean too.

Each characteristic boundary pencil computes primitive states, eigenvectors
and outgoing-wave flags once per stage, shared by its six ghost layers.
Each layer retains the original projection and accumulation order. This
saves decompositions rather than changing the boundary closure.

The three SSP-RK stages update disjoint contiguous rows in parallel for
multidimensional meshes with at least 8,192 interior cells. Smaller meshes
and 1D meshes keep a serial path. Scalar formulas, stage order and ghost
refreshes are preserved. This uses safe slice partitioning and adds no
shared mutable pointers or persistent state buffer.

## Checks

45 tests pass across compressible flow, radiation and semiconductors.
New checks cover:

- A deliberately unsafe GPR energy flux: limiting makes both one-sided
  conservative states admissible, and blends every non-conservative face
  component consistently, including recovery from a non-finite component.
- Characteristic acoustic eigenmodes in all three axes, both boundary sides
  and multiple ghost layers: outgoing modes extrapolate, incoming modes
  take the far field.
- Two steps on all three orientations of a 96² mesh: one and three threads
  produce identical states, and every interior cell remains admissible.
- A strong pressure jump: limiting remains active and admissible, and the
  solver now reports its activation count.

Existing Sod, acoustic-pulse and GPR Stokes validation examples retain
identical printed accuracy metrics and step counts against the baseline.
The forced-all-speed acoustic example still has its existing coarse-grid
error (14.6% of pulse amplitude at 128 cells, versus 0.30% for explicit);
this pass does not improve that time discretization. An end-to-end 32³
water/air render at 128×80, 64 spectral paths/pixel and 12 fluid steps
produces byte-identical HDR output before and after. The raw JSON includes
these comparisons and the HDR SHA-256.

The benchmark open case advances a small 3D acoustic disturbance to t=0.02.
Its timing probes boundary work; it does not establish long-time reflection
quality or improve the existing normal-incidence boundary approximation.
The changes retain the existing RK3 temporal order, sixth-order central
spatial method, shock closures and EOS. Caching preserves the probe results;
the GPR limited-face correction intentionally changes formerly inconsistent
fallback states when that limiter activates.

## Reproduce

From the repository root, build both versions completely before timing:

```sh
git worktree add --detach /tmp/ns-cpu-baseline 572c566
cargo build --manifest-path /tmp/ns-cpu-baseline/fluid/Cargo.toml -p compressible --release --example cpu_validate
cargo build --manifest-path fluid/Cargo.toml -p compressible --release --example cpu_validate
python3 fluid/scripts/cpu-bench.py /tmp/ns-cpu-baseline/fluid/target/release/examples/cpu_validate fluid/target/release/examples/cpu_validate --threads 8 --repeats 4 --output /tmp/ns-cpu-round2.json
cargo test --manifest-path fluid/Cargo.toml -p compressible -p radiation -p semiconductor --release --lib --tests
git worktree remove /tmp/ns-cpu-baseline
```

The benchmark tool alternates executables and fails if any workload's
conserved-state hash differs. Use matching thread counts and power settings,
and keep other heavy work out of the timing interval.
