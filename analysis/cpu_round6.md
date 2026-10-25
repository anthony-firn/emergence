# Sixth CPU pass: heat conduction and reflecting boundaries

Conduction-heavy fluid stepping is 1.16–1.30× faster on the measured 32³
Euler, two-phase, and three-phase fixtures. Every repeated pair retains
identical conserved-state hashes, temperature extrema, energy errors, and
conduction substep counts. Short reflecting axes now fill deep ghost layers
from current physical or prescribed boundary states.
The [raw runs](cpu_round6.json) include timings, regression outputs, and
preliminary measurements.

## Heat conduction

Euler, FiveEq, and Kapila3 evaluate conduction properties with fixed-size
phase arrays, avoiding a phase-vector construction for every cell and
substep. `System::conduction_properties` has a default implementation that
uses existing phase and temperature closures for custom systems.
FiveEq and Kapila3 also use these arrays for temperature queries.

Contiguous interior rows map directly to the padded state, removing the
interior-index vector (256 KiB on a 32³ grid). Computing face rates also
reduces the timestep limit in that traversal. Cells then receive energy
and optional phase-relaxation updates directly, avoiding a collected list
of changed states and a serial scatter. The updates run in parallel at
8,192 cells or more; smaller grids use the serial update path.

The numerical operator uses the same harmonic face conductivities,
half-cell thermal walls, phase-capacity cutoff, and conductance-based
substeps. Relaxation receives the heated state and its substep duration;
its mutations are committed when it returns true. The new tests verify
that mutations returned with false are discarded and that cooling hooks
produce reproducible, physically consistent energy changes.

## Deep reflecting ghosts

The previous mirror indices could land in unfilled opposite ghosts when
an active axis had fewer cells than the ghost width. A uniform resting
Euler state on a two-cell reflecting axis became inadmissible with the
positivity fallback disabled. NaN-filled storage also exposed invalid
sources at mixed boundaries and in multidimensional corners.

The new image resolver reflects coordinates until they reach an interior
cell or the opposite boundary condition. Each reflection reverses normal
momentum. Periodic, outflow, fixed-state, and characteristic data are
resolved from current inputs. Two reflecting walls therefore support
repeated images, including their alternating momentum parity.

A shared implementation has separate compile-time ordinary and deep-image
variants. A single mesh/boundary check selects the uncommon path; the
ordinary ghost-copy kernel contains no deep-image traversal calls.
Inactive one-cell axes continue to have no allocated ghosts.

## Measured performance

Baseline `50664c9`, retained runtime `035602d`; Intel Core i7-13620H,
eight Rayon threads, System76 Battery profile, CPU 10–50%, no turbo.
Release builds use thin LTO. Each case has four alternating before/after
pairs, after compilation and with no concurrent heavy task. Timings cover
the full stepping loop and its allocations.

| Conduction-heavy workload | Previous median | New median | Previous/new |
|---|---:|---:|---:|
| Euler gas 24³, 40 steps | 0.451 s | 0.448 s | 1.01× |
| Euler gas 32³, 30 steps | 0.858 s | 0.737 s | 1.16× |
| Two-phase mixture 32³, 30 steps | 1.256 s | 0.965 s | 1.30× |
| Three-phase mixture 32³, 30 steps | 1.409 s | 1.125 s | 1.25× |

The two- and three-phase 32³ timing ranges are disjoint. The gas ranges
overlap, including an outlier in the 24³ candidate. The measured gains
apply to these conduction-heavy fixtures. All phases in the timing fixtures
use the same ideal EOS, with different phase conductivities and evolving
volume fractions; the accuracy tests separately cover distinct phase heat
capacities. Fixed dt is 1e-4. At 32³ the gas case takes 480 conduction
substeps, and the two- and three-phase cases each take 570, across 30 flow
steps. Relative total-energy errors remain below 6e-15.

| Ordinary flow check | Previous median | New median | Previous/new |
|---|---:|---:|---:|
| Central Euler 32³, 100 steps | 1.349 s | 1.369 s | 0.99× |
| Hybrid water/air 32³, 40 steps | 1.471 s | 1.497 s | 0.98× |
| Characteristic open flow 24³, 60 steps | 0.680 s | 0.717 s | 0.95× |

Ordinary-flow medians are 1–6% slower in this sample, with overlapping
ranges in all three cases. These runs establish identical physical outputs;
they do not establish an ordinary-flow speedup. Earlier boundary layouts
showed larger median slowdowns, which prompted separating the deep-image
path. Their raw results are included for comparison.

## Accuracy and regression checks

68 tests pass across compressible flow, radiation, and semiconductors.
Five new tests verify:

- Every deep reflecting ghost, including corners, against an independent
  mirrored-coordinate map with six and ten ghost layers.
- Reflections into opposite periodic, outflow, fixed-state, and
  characteristic boundaries, in both orientations on every axis.
- Rest-state preservation on short reflecting meshes, including 2×3×4,
  with positivity fallback disabled.
- Second-order cosine heat-diffusion convergence and conservation in two-
  and three-phase systems with different phase heat capacities.
- Thermal relaxation call counts, analytic total-energy decay from a
  cooling source, commit-on-true behavior, and identical one-/three-thread
  states on every 2D orientation and a small 1D grid.

The three new reflecting-boundary tests fail on the previous runtime.
Both new thermal tests pass there, confirming preserved numerical behavior.
WCNS and FD-WENO flow probes, viscous RHS and full stepping, and explicit,
all-speed, and viscous timestep queries also retain identical result hashes
and printed physical metrics. Those comparisons are stored in the raw data.

## Reproduce

From the repository root:

```sh
git worktree add --detach /tmp/ns-thermal-baseline 50664c9
cargo build --manifest-path /tmp/ns-thermal-baseline/fluid/Cargo.toml -p compressible --release --example cpu_thermal --example cpu_validate --example cpu_viscous --example cpu_timestep
cargo build --manifest-path fluid/Cargo.toml -p compressible --release --example cpu_thermal --example cpu_validate --example cpu_viscous --example cpu_timestep
python3 fluid/scripts/cpu-bench.py /tmp/ns-thermal-baseline/fluid/target/release/examples/cpu_thermal fluid/target/release/examples/cpu_thermal --threads 8 --repeats 4 --case gas:24:40 --case gas:32:30 --case two:32:30 --case three:32:30 --output /tmp/ns-thermal.json
python3 fluid/scripts/cpu-bench.py /tmp/ns-thermal-baseline/fluid/target/release/examples/cpu_validate fluid/target/release/examples/cpu_validate --threads 8 --repeats 4 --case smooth:32:100 --case interface:32:40 --case open:24:60 --output /tmp/ns-thermal-flow.json
cargo test --manifest-path fluid/Cargo.toml -p compressible -p radiation -p semiconductor --release --lib --tests
git worktree remove /tmp/ns-thermal-baseline
```
