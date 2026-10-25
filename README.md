# Emergence

Simulation and engineering design with explicit limits on every shortcut.

**[Try the browser lab](https://www.evolvetech.dev/play/emergence/)** ·
**[Project overview](https://www.evolvetech.dev/projects/emergence/)** ·
[Model ladder](analysis/emergence_ladder.md) · [Measured results](analysis/emergence_engine.md)

The organising idea is small: use a cheap model while its validity indicator
meets the requested error, then unfold to a finer model where it fails. Fold
expensive calculations into validated closures when they will be reused. Carry
the error and its rigor through to the result someone is designing against.

This is a Rust research workspace, with CPU reference implementations and
optional CUDA backends. It brings fluid simulation, light and rendering,
materials, electronics, CAD and engineering requirements into one set of
models. Domains have different maturity; the links below state their checked
results and limits. It is not yet a generally coupled world simulator.

## Start here

Install stable Rust, then run the CPU examples without a CUDA toolkit:

```sh
cd fluid
cargo test -p emergence
cargo test -p emergence-web -- --nocapture
cargo run --release -p fluid-cli -- tg --backend cpu --n 24 --re 200 --t-end 0.5 --out /tmp/tg.csv
cargo run --release -p fluid-cli -- project report projects/motor-driver/project.json --why
```

`cargo test` defaults to the engine, browser adapters, CPU spectral solver,
photonics, project graph and CLI. Each crate's own tests can be run separately:
`cargo test -p <crate>`. CUDA crates are opt-in; see
[the CUDA setup](fluid/scripts/cuda-env.sh), then build the CLI with
`cargo build --release -p fluid-cli --features cuda`.

## What is here

| System | Models and design tools | Validation and current limits |
| --- | --- | --- |
| Emergence engine | Typed claims, cheapest admissible rung, reference audits, adaptive closure tables, goal-oriented error budgets | [Engine results](analysis/emergence_engine.md); only some domain ladders use the shared engine today |
| Air and fluid motion | f64 spectral incompressible flow; compressible flow, buoyancy/convection, shocks, two-phase interfaces, wetting, fluid–structure coupling | [Regimes](analysis/regime_map.md), [compressible flow](analysis/compressible.md), [convection](analysis/convection2d.md), [water/air interfaces](analysis/twophase_and_mhd.md). These are building blocks, not a validated weather forecast model |
| Light, atmosphere and water optics | Rayleigh/Mie sphere scattering, optical stacks, wave optics, photons and radiation transport | [Mie](analysis/mie.md), [radiation](analysis/radiation.md), [wave optics](analysis/double_slit.md). Optical constants and independent-particle assumptions are stated in the models |
| Rendering | CAD sphere tracing, participating media, spectral emission and evolving fluid snapshots | [Renderer](analysis/viz.md), [fluid rendering](analysis/radiation.md), [photographic scene](analysis/photographic_scene/README.md). Appearance estimates and unresolved realism are documented |
| Chip and CPU design | Device/cell models, a synthesizable Verilog subset, equivalence checks, timing and power, a small CPU running Fibonacci | [Chip flow](analysis/chip.md), [mini-CPU project](fluid/projects/mini-cpu). One technology, ideal clock, no physical layout |
| CAD and geometry | Parametric parts and assemblies, signed distance, STEP/IGES via gmsh, STL/OBJ, immersed boundaries, forces, probes and VTK | [CAD](analysis/cad.md), [engineering flow](analysis/real_geometry.md). Uniform-grid laminar immersed flow |
| Engineering projects | Unit-checked parameters, tolerances, shared materials, team-owned requirements, electrothermal feedback, HVAC and printable parts | [Project graph](analysis/project_model.md), [HVAC](analysis/hvac.md), [motor driver](fluid/projects/motor-driver), [office](fluid/projects/office-hvac) |
| Materials and heat | Atomistic → continuum elasticity, transport, fracture, solidification, polymers, magnetism, batteries and phonon heat transport | [Materials overview](analysis/materials_overview.md); every domain retains its benchmark and assumptions |
| Gravity and solar models | Isolated gravity, particle transport, published solar background, surface flux and continuum rendering | [Gravity](analysis/gravity_solar.md), [solar rendering](analysis/solar_render.md). Prescribed solar flows/sources, no cycle or flare forecast |

## Why the indicator matters

These are previously measured native results, retained with their original
benchmarks and hardware context. They are not browser performance claims.

| Task | Observed result | Evidence |
| --- | --- | --- |
| Cloud scattering field | 10.7× faster than the full Mie series after a one-off closure build; every sampled node agreed to 4.7 × 10⁻¹² | [Engine results](analysis/emergence_engine.md#mie-scattering-photonicstestsmie_ladderrs) |
| Certified fluid resolution | At Re = 200, t = 0.5, a 10⁻³ requested L² tolerance selected 16³ with a 3.5 × 10⁻⁴ bound; distance from the 48³ comparison was 8.7 × 10⁻⁵ | [Certified resolution](analysis/emergence_engine.md#rigorous-certified-resolution-spectral-coretestscertified_ladderrs) |
| Heat transport budget | For the stated independent-grey-mode film model, goal-weighted error budgets cost 3.6 × 10⁵ less than a uniform per-mode tolerance | [Heat model and limitations](analysis/emergence_engine.md) |
| Engineering hand-off | The motor driver's trace-drop nominal passes 50 mV, but its interval crosses the limit: the requirement is reported **undecided** | [Project model](analysis/project_model.md) |

Claims are tagged **asymptotic**, **validated**, or **rigorous** (up to the
floating-point qualifications in the relevant model). An unmet request stays
unmet. A claim about numerical error inside a model is not a bound on its
physical modelling error. The native adaptive fluid ladder also has an honest
negative result: on Taylor–Green it did not beat the best fixed grid.

## Browser lab

The lab compiles the existing Rust engine and models to WebAssembly. It runs
inside a Web Worker, with no simulation server or external runtime downloads:

1. **Light:** Rayleigh expansion → resonance-aware folded table → Mie series.
   Tighten tolerance or change size, then inspect the selected rung, claimed
   error, actual reference discrepancy and an intentionally unreachable request.
2. **Fluid:** the actual f64 spectral solver on a small 16³–32³ periodic grid.
   The planar Taylor–Green preset checks energy against
   `E(t) = ¼ exp(−4νt)`; the 3D preset shows resolution warnings separately.
3. **Design:** the motor-driver thermal/electrical graph, with its original
   tolerances and MOSFET feedback. Change fins and airflow; see requirements
   become met, undecided or unavailable outside the laminar regime.

Build the static files:

```sh
rustup target add wasm32-unknown-unknown
python3 web/build.py
```

`web/dist/` is intended for `/play/emergence/` on evolvetech.dev and uses that
site's shared CSS, theme and fonts. To integrate a local checkout:

```sh
python3 /path/to/evolvetech.dev/tools/build_demos.py emergence --local-source /path/to/emergence
python3 /path/to/evolvetech.dev/tools/check_site.py
```

[Browser adapter and validation](analysis/browser_lab.md) describes the ABI,
exact checks and display precision. The browser demo is illustrative at small
grids; it does not claim automatic rigorous fluid refinement or a fully coupled
atmosphere, water and light simulation.

## Reproducibility and sources

Code is in `fluid/crates/`, example designs in `fluid/projects/`, and the
benchmark write-ups and retained results in `analysis/`. Material constants
cite a publication, handbook or an explicit estimate in their source comments.
See [third-party data and assets](SOURCES.md) and the per-domain documentation.

This repository starts with a clean snapshot of the simulation and design
work previously developed in a local Navier–Stokes audit checkout. The proof
audit, downloaded proof repositories, papers and their git history are not
part of this project. No history of the original checkout was rewritten.

Code is licensed under [Apache 2.0](LICENSE). Third-party reference data and
generated appearance assets retain the provenance described in `SOURCES.md`.
