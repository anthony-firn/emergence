# GPU molecular dynamics: keeping the atoms resident

Implemented in `fluid/crates/md-cuda` (new crate; depends on `md` read-only
and on `materials` for the Lennard-Jones parameters). This is the
GPU backend of the molecular rung of [emergence_ladder.md](emergence_ladder.md),
whose CPU engine and closures are documented in [md.md](md.md). The point
is the same fold as there — atoms to a continuum with an equation of state
and transport coefficients — but with the fine model held on the device so
that the cheap model can be sampled at 10⁵–10⁶ atoms.

The design follows `compressible-cuda`: CUDA kernels are NVRTC source
strings (`src/kernels.cu`) compiled at runtime for `sm_86`, instantiated for
both `f64` and `f32` from one templated source, and the CPU engine is the
reference for every kernel.

## What is on the device

`GpuMd<R>` owns the whole state across steps:

| buffer | content |
|---|---|
| `pos`, `vel`, `force` | `[n][3]` atom-major |
| `pe` | per-atom pair energy (already halved) |
| `stress` | per-atom unsymmetrised pair virial `[n][6]` |
| `counts`, `starts`, `order`, `cell_of` | cell list (`i32`) |
| `diag` | 13 reduction slots (`pe`, `stress[6]`, kinetic tensor `[6]`) |

Only diagnostics (`positions()`, `forces()`, `potential_energy()`,
`kinetic_energy()`, `pressure_tensor()`, `internal_pair_virial()`,
`cell_stats()`) copy anything back; `step()` is `kick`–`drift`–`force`–`kick`
entirely on the device, exactly mirroring `md::System::step_nve` (unit mass,
positions not wrapped, minimum image applied to every candidate pair).

The neighbour scheme is a **plain cell list rebuilt on every force
evaluation** (no Verlet skin), which is the right trade on a device: the
build is cheap and it keeps the CPU and GPU trajectories on the same
algorithm. Kernels:

- `cell_counts`: bin atoms, `atomicAdd` per cell;
- `cell_scan`: exclusive scan of the counts into offsets. Written as a
  *blocked* scan — one 256-thread block loops over the cells in chunks with
  a shared-memory Hillis–Steele scan and a running carry. A serial scan is
  `O(ncells)` on one thread and dominated the whole step at 256 k atoms
  (~20 000 cells); the blocked scan removed that.
- `cell_fill`: scatter atom indices using the offsets as a cursor;
- `force`: visit the up-to-27 neighbouring cells, `Cutoff` schemes 0/1 on
  the device, scheme 2 (bare) with the analytic tail correction added on
  the host;
- `kick`, `drift`: velocity-Verlet;
- `diag`: `atomicAdd` reductions for the diagnostics.

The grid is chosen as `nc = floor(L/rc)` per axis, so the cell width is
`≥ rc` and `< 2 rc`; the neighbour loop therefore visits adjacent cells
only. `nc ≤ 2` is special-cased (the `−1` and `+1` images coincide), which
is what makes a tiny `fcc(3)` box (2×2×2 = 8 cells) exact. Like the CPU engine
this is a *full* list, so each pair is visited from both endpoints and the
per-atom energy/virial are halved when totalled.

`gpu_available()` is panic-safe: under cudarc's `dynamic-loading` a missing
driver makes `CudaContext::new` panic rather than return `Err`, so it probes
for the CUDA library and wraps the context creation in `catch_unwind`. Every
test returns early with a printed `[skip]` note when it is false.

## Correctness

`cargo test -p md-cuda` (7 tests, ~22 s). All values are `f64` GPU against
the `md` CPU engine on the *same* configuration; relative errors are
`max|a−b| / max|a|` for forces and `|·|/|reference|` for energies/virials.

| configuration | N | cells (nc) | force rel | PE rel | virial rel |
|---|---|---|---|---|---|
| fcc + 0.05 σ displacement | 108 | 8 (2) | 3.27e-15 | 1.79e-16 | 1.93e-16 |
| liquid snapshot (melted fcc) | 108 | 8 (2) | 1.96e-15 | 3.97e-16 | 4.92e-16 |
| fcc k=5 (boundary case) | 500 | 27 (3) | 3.44e-15 | 6.19e-16 | 3.33e-16 |
| fcc k=11 (multi-chunk scan) | 5324 | 343 (7) | 3.20e-15 | 2.79e-15 | 3.14e-15 |
| tail-corrected fcc k=4 | 256 | 8 (2) | 1.66e-15 | 1.32e-16 | 0 (pair virial excludes tail) |

The `fcc(3)` perfect lattice is deliberately displaced by ±0.05 σ before the
force comparison: a perfect lattice has cancelling forces, so the
relative-force metric would be measuring the round-off of a zero. The
k=11 case exercises the scan carry across chunks of 256 cells; the k=5 case
exercises `nc = 3` images. Everything agrees to a few 10⁻¹⁵. The scan is
additionally pinned down on its own: a dilute random configuration on a
24³ = 13 824-cell grid (54 scan chunks) reproduces the CPU per-cell
occupancy **exactly, integer for integer** (`cell_list_scan_matches_cpu_binning`);
a carry bug would have moved atoms to the wrong cells and shown up as a
mismatch.

**NVE trajectory, 108-atom fcc at ρ = 0.8, dt = 0.002, 2000 steps.** The
initial total energies differ by 2.4e-16 relative; the first-step position
error is 1.1e-16. Chaotic divergence (this is a round-off separation, not a
bug) grows as

| step | 1 | 10 | 100 | 200 | 500 | 1000 | 2000 |
|---|---|---|---|---|---|---|---|
| max │x_GPU − x_CPU│ | 1.1e-16 | 1.1e-16 | 8.9e-16 | 1.8e-15 | 1.2e-14 | 3.3e-13 | 5.7e-11 |

The growth rate is set by the Lyapunov exponent of the LJ liquid: each
factor of ~10 in time costs roughly one digit, so a round-off match can
only be expected for a few hundred steps, not indefinitely. Run-to-run the 2000-step
value varies at the 10⁻¹¹ level (4.3–5.7e-11 across runs) because the cell
list order comes from `atomicAdd` in `cell_fill`, so the order of terms in
the force sum is not reproducible; this is exactly the round-off-scale
nondeterminism the method is designed to have. The CPU trajectory's own
energy deviation over the run is 1.5e-5, and the GPU `f64` deviation is
identical to the digits shown — the drift is integrator truncation, not
precision.

**Single precision.** Against the CPU `f64` force on a displaced fcc, the
`f32` force is off by max 2.9e-5 absolute, 3.6e-6 relative — the `f32`
machine epsilon (`6e-8`) amplified by cancellation in the pair sum, and the
level at which `f32` can be used safely for force-limited work. Over the
same 2000-step NVE run the `f32` maximum `|ΔE|/|E|` is 1.7e-5 and the `f64`
value is 1.2e-5; the end-state energy `f32` vs the initial CPU energy is
1.3–1.8e-5. So over a few thousand steps the `dt` truncation error
(~10⁻⁵) dominates the single-precision round-off, and `f32` looks as good
as `f64`.

A 50 000-step run (dt = 0.002, 100 reduced time units) makes the split
between *energy* and *trajectory* explicit, and is the honest statement of
what `f32` is good for. The initial `f32` energy differs from `f64` by
1.1e-7, and the maximum energy deviations are

| horizon (steps) | 5 000 | 20 000 | 50 000 |
|---|---|---|---|
| f64 max │ΔE│/│E│ | 2.2e-6 | 5.5e-6 | 1.3e-5 |
| f32 max │ΔE│/│E│ | 2.5e-6 | 1.7e-5 | 3.3e-5 |
| max │x_f32 − x_f64│ (σ) | 2.3 | 5.6 | 9.6 |

Energy is conserved only ~2.5× worse in `f32` even after 50 000 steps —
`f32` is fine for thermodynamic averages, which do not need the trajectory.
But the `f32` and `f64` trajectories separate to ~2σ within 5 000 steps and
to ~10σ by 50 000, so `f32` cannot be used to continue a specific `f64`
trajectory: it is a different (statistically equivalent) realisation. This
is the same message as the `f64` GPU-vs-CPU divergence above, just at a much
larger amplitude.

## Performance

`cargo run -p md-cuda --release --example perf -- [N ...]` (defaults
4000, 32 000, 256 000). Sizes are fcc at ρ = 0.8, rc = 2.5, dt = 0.005
reduced = 0.0108 ps; the CPU engine is `md` with rayon, the GPU rebuilds
the cell list every step. `SKIP_CPU=1` / `SKIP_GPU=1` restrict the run.
Each step is timed separately and the *minimum* and *mean* are reported:
minimum filters scheduler jitter, mean is the amortised throughput.

Argon time unit `σ√(m/ε) = 2.156e-12 s`; `ns/day` is the standard
molecular-dynamics unit and is `steps/s × dt × 2.156e-12 × 86400 / 1e-9`.

### The environment matters

The test GPU is an RTX 3050 Laptop (GA107, 6 GB), driver 580.173.02, CUDA
13.0, kernels built for `sm_86`. Two things make the absolute numbers
pessimistic and must be kept in mind:

- **The card is software power-capped to 15 W** (`nvidia-smi`: current and
  requested limit 15 W, default 35 W, max 65 W) and during the sweeps the
  SM clock sat at **210 MHz under 100 % load** against a 2100 MHz maximum.
  So the GPU is running at ~1/10 of its nominal clock, with the memory at
  405/5501 MHz.
- **Other CUDA agents share the card.** A co-resident long-running kernel
  inflates a step from 4 ms to 200 ms for the whole measurement window, so
  the best-of-step filtering cannot rescue a run that overlaps one. The
  numbers below were taken with `nvidia-smi --query-compute-apps` reporting
  no co-resident process; a run started while the card was busy is ~30–65×
  slower and is discarded.

The numbers are therefore honest **lower bounds** for this hardware; the
relative CPU/GPU crossover is the robust part.

### Table (idle GPU, uncontended)

`steps/s` and `ns/day` from the best step; parenthesised values are the
mean-step values.

| N | engine | steps/s (best) | ns/day | steps/s (mean) | best s | mean s |
|---|---|---|---|---|---|---|
| 4 000 | CPU f64 | 741 | 690 | 232 | 1.35e-3 | 4.32e-3 |
| 4 000 | GPU f64 | 44.7 | 41.6 | 43.1 | 2.24e-2 | 2.32e-2 |
| 4 000 | GPU f32 | 641 | 597 | 491 | 1.56e-3 | 2.04e-3 |
| 32 000 | CPU f64 | 101 | 94.1 | 3.3 | 9.90e-3 | 3.04e-1 |
| 32 000 | GPU f64 | 9.1 | 8.5 | 8.1 | 1.09e-1 | 1.24e-1 |
| 32 000 | GPU f32 | 241 | 224 | 210 | 4.15e-3 | 4.77e-3 |
| 62 500 | CPU f64 | 44.7 | 41.6 | 0.75 | 2.24e-2 | 1.34 |
| 62 500 | GPU f64 | 5.5 | 5.2 | 4.9 | 1.80e-1 | 2.04e-1 |
| 62 500 | GPU f32 | 159 | 149 | 139 | 6.28e-3 | 7.21e-3 |
| 256 000 | CPU f64 | — | — | — | >150 (build) | — |
| 256 000 | GPU f64 | 1.5 | 1.4 | 1.0 | 6.73e-1 | 1.02 |
| 256 000 | GPU f32 | 39.2 | 36.5 | 28.0 | 2.55e-2 | 3.57e-2 |

The CPU mean is the amortised throughput including neighbour-list rebuilds
(`md` rebuilds when an atom has moved more than `skin/2 = 0.15 σ`, i.e.
every few tens of steps at T ≈ 1), so it is the fair end-to-end number. Its
spread between best and mean is large and grows with N — the rebuild is
the CPU bottleneck. `fcc(40)` = 256 000 atoms did not finish the CPU
neighbour-list build inside a 150 s timeout on the shared machine, so the
CPU 256 k entry is missing rather than guessed; extrapolating the measured
best-step scaling (741 → 101 → 44.7 steps/s) gives roughly 10 steps/s,
labelled an extrapolation.

### Crossover and speedup

- **GPU f32 vs CPU f64.** At 4 000 atoms the GPU is *slower*
  (641 vs 741 steps/s, 0.86×) — kernel-launch and host round-trip overhead
  on a 200 µs CPU step. By 32 000 atoms it is 2.4× faster and by 62 500 it
  is 3.6×. The crossover is around **5 000–8 000 atoms**. On the amortised
  (mean) throughput the speedup is larger still (≈ 64× at 32 000, order
  10²× at 62 500) because the CPU neighbour-list rebuild dominates and the
  GPU rebuilds cheaply every step — but the mean depends on where a rebuild
  happens to fall in the sampled window, so the best-step ratio is the
  stable measure and the amortised ratio is order-of-magnitude only.
- **GPU f64 vs CPU f64.** The GPU `f64` path *never* wins on this consumer
  card: 0.06× at 4 k, 0.09× at 32 k, 0.12× at 62.5 k, and even at 256 k the
  extrapolated CPU (~10 steps/s) beats 1.5 steps/s. This is the 1/64
  double-precision rate of GA10x; the gap narrows slightly with N as the
  kernel becomes memory-bound but the CPU always leads. `f64` is for
  checking, `f32` is for production.
- **GPU f32 vs GPU f64.** 14× at 4 k rising to ~26× at large N, consistent
  with the nominal 1/64 throughput being partly hidden by memory access.

The 256 k GPU run also shows the cell-list cost is under control: the
first version paid a serial `O(ncells)` scan and stalled at 256 k; the
blocked scan brings the 256 k `f32` step to 25.5 ms. A per-kernel
breakdown (`cell_stats` vs `compute_forces`, min of several reps) gives,
for the clean `f32` paths: 4 000 atoms 0.06/1.6 ms, 32 000
0.18/3.8 ms, 62 500 0.38/5.6 ms, 256 000 1.6/20.5 ms, where the first
number is the cell-list build and the second the whole force evaluation —
so the cell list is ~4–8 % of the force evaluation and a Verlet skin would
buy little. The
remaining cost is the neighbour traversal itself, which is gather-bound
(positions fetched from ~13 atoms per cell × 27 cells).

## Files

- `fluid/crates/md-cuda/src/kernels.cu` — all kernels, templated `f64`/`f32`.
- `fluid/crates/md-cuda/src/lib.rs` — `GpuMd<R>`, `gpu_available()`, diagnostics.
- `fluid/crates/md-cuda/tests/md_cuda.rs` — the correctness tests above.
- `fluid/crates/md-cuda/examples/perf.rs` — the performance sweep.


## Re-measured with the Performance power profile (coordinator, 2026-10-07)

The tables above were taken while the laptop's system76 profile had
reverted to Battery (GPU pinned at 210 MHz, CPU capped at 50% with no
turbo). Re-run of `examples/perf.rs 32000 62500` with the Performance
profile (GPU at 990 MHz, 9.6 W during the run; the CPU was shared with
three other agents' builds, so the CPU column is somewhat low):

| N | CPU f64 (steps/s) | GPU f64 | GPU f32 | GPU f32 speedup |
|---|---|---|---|---|
| 32 000 | 87.1 | 72.5 | 1875.8 | 21.5× |
| 62 500 | 47.0 | 45.0 | 1294.2 | 27.5× |

GPU f32 reaches 1747 / 1206 ns/day (argon units, dt = 10.8 fs); GPU f64
is at parity with the CPU on this consumer card.

## Green–Kubo transport on the device (size and cutoff study)

This section extends the crate with the GPU pressure-tensor sampling and a
Green–Kubo (GK) / Einstein–Helfand (EH) estimator, and uses them to test
whether the CPU miss of the Lennard-Jones triple-point shear viscosity
(`η* = 2.45 ± 0.41` from `md` at `N = 256`, `rc = 2.5`, versus Meier,
Laesecke & Kabelac 2004 I: `3.258 ± 0.033`) closes with longer runs, larger
systems and a larger cutoff. The reference state is `T* = 0.722`,
`ρ* = 0.8442`, `dt* = 0.004`, `Cutoff::TailCorrected` (bare truncation plus
the analytic isotropic tail). Reduced units, `k_B = 1`.

### Per-step pressure tensor, f64 accumulation

`md_press_f64/f32` in `kernels.cu` reduces the instantaneous numerator
`out[k] = Σ_i [m v_i v_j + ½ stress_i[k]]` for
`k = xx, yy, zz, xy, xz, yz`, where `stress_i` is the per-atom pair virial
already summed over neighbours. The reduction is an `atomicAdd` on
`double`, so the **f32 force path still accumulates the virial sum in f64**
— the f32 round-off enters through the trajectory, not through the
diagnostic. `GpuMd::reserve_pressure_samples` allocates a device buffer and
`sample_pressure_tensor` appends one 6-vector per step (no host
synchronisation in the hot loop); `take_pressure_samples` divides by `V`
and adds the isotropic tail pressure on the diagonal.

Validation against the CPU engine (`cargo test -p md_cuda`, GPU required):

- `pressure_tensor_matches_cpu`: on a liquid configuration the f64 GPU
  tensor matches `md::System::pressure_tensor` to `< 1e-11` per component;
  the f32 path to `< 2e-4`.
- `gk_running_integral_matches_cpu_md_crate`: over a 400-step f64 NVE
  trajectory the GPU pressure samples match the CPU to `3.1e-13`
  component-wise, and the GPU block estimator (one block) reproduces
  `md::transport::running_gk` applied to the same samples to `< 1e-9`.
- The full `md-cuda` suite (10 integration + 4 unit tests) runs in ~6 s.

### Estimator and error bars

`gk.rs` splits the trajectory into contiguous blocks, forms the
autocorrelation of each off-diagonal component inside the block, and
returns the per-lag standard error of the mean of the running integral
across blocks. The three off-diagonal components are averaged. Because a
single block gives no error bar, production runs use independent seeds and
the quoted error is the standard error across seeds (each seed's curve is
also block-analysable). The plateau is located as the peak of the running
integral smoothed over 0.4 reduced time units (the stress correlation time
is `≈ 0.14`); the estimator also reports a fixed integration window
`t ∈ [1, 4]`, where all curves are flat, as a cross-check. The EH viscosity
is the slope of `〈[∫₀ᵗ P_xy]²〉` over the same window.

### Size and cutoff study

The `examples/gk_study.rs` example melts the metastable fcc at `T* = 1.5`,
equilibrates at `T* = 0.722` by velocity rescaling, then runs microcanonical
production. Pressure tensors are sampled every 2 steps and unwrapped
positions every `steps/400` steps for the MSD. Runs are `40 000` steps
(`30 000` for the `N = 32 000` rc = 2.5 case), except the `15 000`-step
self-diffusion sweep. η plateau, `t ∈ [1, 4]`, with the seed SEM at `t* = 2`
(flat region) in brackets; `η_EH` over the plateau:

| `N` | seeds | `rc` | `η*` (win.) | `η*`(`t*=2`) | `η*_EH` |
|---|---|---|---|---|---|
| 500 | 16 | 2.5 | 3.11 | 3.184 ± 0.089 | 3.243 ± 0.118 |
| 2048 | 12 | 2.5 | 3.16 | 3.235 ± 0.097 | 3.214 ± 0.120 |
| 8788 | 4 | 2.5 | 3.33 | 3.321 ± 0.126 | 3.648 ± 0.315 |
| 32000 | 4 | 2.5 | 3.03 | 3.001 ± 0.156 | 3.289 ± 0.166 |
| 2048 | 3 | 4.0 | 3.07 | 2.949 ± 0.166 | 3.413 ± 0.160 |

All rc = 2.5 values are consistent with Meier's `3.258` within ~1.5σ; the
seed-to-seed scatter is the limiting error, not the size. **The miss
closes.** The original `2.45 ± 0.41` came from three 6 000-step, 256-atom
seeds; with 16 seeds of 40 000 steps at `N = 500` the same estimator gives
`3.18 ± 0.09`, and the value is not size-sensitive. This is a
statistics-limited miss, not a finite-size or cutoff one.

The `rc = 4.0` case is more subtle. The `N = 2048` running integral keeps
rising slowly to `t* = 12` (reaching `3.69` at the end of the window)
instead of flattening, which is the known slow convergence of the GK
integral through the long-time stress tail. The fixed `t ∈ [1, 4]` value is
`3.07 ± 0.13`, consistent with `rc = 2.5`; the late-time value is not a
converged plateau and should not be quoted as η. A larger cutoff therefore
does not move the plateau estimate outside the rc = 2.5 band on the
trajectory lengths affordable here; it mainly lengthens the unresolved
tail. (The earlier claim that `rc = 2.5` sits low relative to `3.258` is not
supported: Meier's own `rc = 2.5` simulations give `η* ≈ 3.0–3.1`, and the
GPU rc = 2.5 plateau is `3.0–3.3`.)

### Self-diffusion and the Yeh–Hummer correction

`D` is measured from the MSD of unwrapped positions with the CPU-`md`
compatible window `t ∈ [0.3, 6.0]` (`D_short`); a whole-trajectory fit
(`D_long`) is biased high in the f32 runs at larger cutoff and is not used.
Values are seed means ± SEM; `D0` applies Yeh & Hummer's
`D0 = D_pbc − T ξ/(6 π η L)`, `ξ = 2.837297`, with the measured `η`.

| `N` | `L/σ` | `D_pbc` rc = 2.5 | `D0` rc = 2.5 | `D_pbc` rc = 4.0 | `D0` rc = 4.0 |
|---|---|---|---|---|---|
| 500 | 8.40 | 0.03162 ± 0.00039 | 0.0276 | 0.03310 ± 0.00054 | 0.0289 |
| 2048 | 13.44 | 0.03291 ± 0.00009 | 0.0304 | 0.03307 ± 0.00019 | 0.0305 |
| 8788 | 21.39 | 0.03335 ± 0.00016 | 0.0318 | 0.03383 ± 0.00034 | 0.0322 |
| 32000 | 33.60 | 0.03350 ± 0.00005 | 0.0325 | 0.03366 ± 0.00002 | 0.0324 |

For comparison, the CPU `md` crate at `N = 500`, `rc = 3.0` (the
`long_runs` example, short window) gives `D_pbc* = 0.0323` and `D0* =
0.0264`, matching the GPU `rc = 3.0` `D_short = 0.0324` — the GPU and CPU
diffusion agree once the system and fit window match.

**Negative result on Yeh–Hummer.** The correction is the right order of
magnitude (`0.004` at `N = 500` falling to `0.001` at `N = 32 000`), but it
does *not* flatten the data: the measured `D_pbc` **increases** weakly with
`L` (from `0.0316` to `0.0335` at rc = 2.5), whereas the Yeh–Hummer sign
(`D_pbc = D0 + T ξ/(6 π η L)`) predicts a decrease. The corrected `D0`
therefore still rises with `L` (`0.0276 → 0.0325`). Either the correction is
not asymptotic at `L ≲ 8σ`, or the small-box measurements carry a residual
systematic (residual structure left by the melt, or the f32 trajectory).
The `rc = 4.0` data are flatter (`0.0331 → 0.0337`), which points at the
small-box rc = 2.5 `N = 500` point being the anomalous one. `D0 ≈ 0.028–0.032`
lies inside the Stokes–Einstein band `0.023–0.035` computed from Meier's
viscosity; a tabulated Meier self-diffusion is not available in this
repository for a direct comparison.

### What caused the miss

Statistics. The LJ triple-point stress autocorrelation has a long tail, so
a single short trajectory can give η anywhere in `2.0–4.3`; three
6 000-step seeds gave `2.45 ± 0.41`, while 16 seeds of 40 000 steps give
`3.18 ± 0.09` at the same `N = 500`, `rc = 2.5`. Finite size at rc = 2.5 is
weak (the four sizes agree within the seed errors), and a larger cutoff
lengthens the unresolved GK tail rather than moving the plateau. The
per-step f64 pressure reduction and the block/seed error bars are the
prerequisites for seeing this.

### Files

- `fluid/crates/md-cuda/src/gk.rs` — block-averaged GK/EH estimators.
- `fluid/crates/md-cuda/examples/gk_study.rs` — size/cutoff production runs.
- `kernels.cu` — `md_press_*` (f64 reduction), `md_scale_*`.
- `lib.rs` — `sample_pressure_tensor`, `take_pressure_samples`,
  `pressure_tensor_f64`.
