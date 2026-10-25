# Third CPU improvement pass

The [fourth pass](cpu_round4.md) corrects deep periodic ghost layers on
short axes and adds scalar conservation checks.

This pass reduces timestep-query work and fixes two stability estimates:
all-speed CFL queries now evaluate current boundary states, and heat
conduction chooses substeps from the actual face conductances. It builds
on the [previous fluid CPU improvements](cpu_round2.md).

## Measured performance

Baseline `efb293b` contains the previous runtime and the new `cpu_timestep`
probe; candidate `ec94d2f` includes this pass. Intel Core i7-13620H, eight
Rayon threads, System76 Battery profile, CPU limited to 10–50%, no turbo.
Release builds use thin LTO. Each workload ran four alternating before/after
pairs, after all builds completed and without overlapping heavy work.
[Raw runs, ranges, hashes and accuracy comparisons](cpu_round3.json) are retained.

| Workload | Previous | New | Previous/new |
|---|---:|---:|---:|
| Multiphase explicit, 32³, 300 CFL queries | 0.780 s | 0.280 s | 2.79× |
| Multiphase explicit, 48³, 200 CFL queries | 1.670 s | 0.658 s | 2.54× |
| Multiphase all-speed, 32³, 300 CFL queries | 0.330 s | 0.142 s | 2.33× |
| Viscous Euler, 32³, 500 CFL queries | 0.107 s | 0.095 s | 1.12× |
| Multiphase explicit, 32³, 60 CFL queries + RK steps | 2.367 s | 2.186 s | 1.08× |

The first four rows time repeated CFL queries on a fixed state, including
Rayon initialization. They measure timestep selection rather than a complete
simulation. The all-speed probe uses a saved acceleration of 1 and an
unlimited acoustic Courant cap, as after a pressure solve; it does not time
the pressure solve itself. Its first-query boundary correction does not
change the global limiting timestep for this particular periodic field.
The viscous case uses a constant-viscosity diffusion limit.

The final row queries CFL and advances one full explicit RK step each time,
with physical step capped at 1e-5. Its 8% median improvement is modest and
its ranges overlap (previous 2.136–2.373 s, new 2.038–2.367 s), so the
whole-step gain remains less certain than the separated query timings.
These measurements do not establish a new rendering speedup.

Every repeated pair has the same hash of all queried timestep bit patterns
and every final interior conserved component. The JSON's existing `steps`
field is the query count for this probe; `state_hash` includes the timestep
bits too. The benchmark script retains its original default fluid workloads
and now accepts repeatable `--case NAME:N:COUNT` overrides.

## Work removed

Explicit CFL selection combines acoustic, non-conservative source and
viscous density bounds in one grid scan. Each center cell converts to
primitive variables once. Velocity-only neighbor queries read momentum and
mass directly for Euler, FiveEq, Kapila3 and GPR, using the same division as
the primitive conversion. Interior pairs use direct references with the
mesh stride; only boundary cells construct boundary states. Custom `System`
implementations retain a primitive-conversion default.

All-speed selection similarly collects the viscous minimum density in its
existing scan. It skips the pressure gradient when the previous solve's
acceleration replaces that estimate, and skips the sound-speed calculation
when the acoustic cap is infinite. Both its pressure-gradient and source
bounds evaluate the current periodic, reflecting, outflow, prescribed or
characteristic boundary closure instead of stale ghost storage.

Conduction evaluates phase properties once per cell and reuses them to
compute temperature. Its temperature/conductivity/capacity and heat-rate
arrays are allocated once per call and reused across substeps. Custom
temperature closures keep their existing default behavior.

## Corrected stability and accuracy checks

51 tests pass across compressible flow, radiation and semiconductors.
Six new tests cover the following:

- A resting Euler grid with a prescribed high-pressure boundary receives
  the analytic first-step acceleration bound. Poisoning ghosts with NaNs
  leaves the timestep unchanged.
- A reflecting wall in a stiff mixture retains its analytic source limit
  in all-speed mode even with saved acceleration; poisoned ghosts again
  leave the result unchanged.
- Direct conserved velocities match primitive velocities bit for bit across
  the four supported flow systems, density ranges and nearly absent phases.
  Cached Euler and multiphase temperatures also match their closures.
- A hot cell beside a cold wall previously fell from T=100 to T=-18.8 in
  one accepted conduction substep. The corrected bound takes multiple
  substeps and keeps all temperatures within the initial/wall range.
- A low-capacity cell between more conductive neighbors previously also
  went negative. The corrected update keeps temperature bounded and
  conserves closed-domain energy to absolute error below 1e-13.
- Periodic cosine conduction retains second-order convergence and
  conserves energy to absolute error below 1e-13:

| Cells | Temperature RMS error at t=0.02 |
|---:|---:|
| 16 | 9.1351e-4 |
| 32 | 2.2813e-4 |
| 64 | 5.7015e-5 |

The heat substep now uses
`h <= 2*safety * min(Cv_volume / sum(face_conductance))`, including harmonic
conductivities and half-cell wall distances. With `0 < safety <= 0.5`, a
frozen positive-capacity temperature update is a convex combination.
The old cell-only estimate missed both the wall contribution and conductivity
contrast. Invalid conductivities, capacities or non-finite face rates now
fail explicitly. This row bound does not prove a maximum principle for
every nonlinear mixture/latent-heat closure; phase relaxation remains
applied after each substep.

Separate old/new smooth-flow and water/air-interface probes retain identical
state hashes and printed accuracy metrics. Sod and forced-all-speed acoustic
pulse examples at 128 cells retain identical errors and step counts. The
existing all-speed pulse error remains 14.6% of amplitude, versus 0.30% for
explicit stepping; this pass improves boundary timestep estimates rather
than that time discretization. A 32³ evolved water/air snapshot rendered at
128×80, 64 spectral paths/pixel, after 12 fluid steps produces byte-identical
HDR output. Its metadata and SHA-256 are in the JSON.

## Reproduce

From the repository root:

```sh
git worktree add --detach /tmp/ns-cpu-baseline efb293b
cargo build --manifest-path /tmp/ns-cpu-baseline/fluid/Cargo.toml -p compressible --release --example cpu_timestep
cargo build --manifest-path fluid/Cargo.toml -p compressible --release --example cpu_timestep
python3 fluid/scripts/cpu-bench.py /tmp/ns-cpu-baseline/fluid/target/release/examples/cpu_timestep fluid/target/release/examples/cpu_timestep --threads 8 --repeats 4 --case mixture:32:300 --case mixture:48:200 --case allspeed:32:300 --case viscous:32:500 --case advancing:32:60 --output /tmp/ns-cpu-round3.json
cargo test --manifest-path fluid/Cargo.toml -p compressible -p radiation -p semiconductor --release --lib --tests
git worktree remove /tmp/ns-cpu-baseline
```

Complete both builds before timing, keep matching thread counts and power
settings, and exclude other heavy work during the comparisons.
