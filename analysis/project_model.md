# The shared project model: one graph for every team

Crate `fluid/crates/project`, command `fluid project`, example
`fluid/projects/motor-driver/`.

## Why

Cross-functional hardware work fails at the hand-offs, not inside the tools.
Electrical sends "about 9 W" to thermal. Thermal reads a datasheet copper
that differs from the copper electrical used. Mechanical trims two fins to
save weight without knowing that the junction-temperature limit sits on that
sink. Every team is right in its own tool, and the product still misses its
requirements.

This crate makes the hand-off itself the product. One plain JSON file per
product, versioned in git, holds every team's parameters, the analyses and
CAD parts that use them, and the requirements they must meet. Every value
that crosses from one team to another carries four things:

| carried | what it prevents |
|---|---|
| **dimension** (SI, strict parsing) | unit mistakes between teams; `degC` is only legal alone, so 5 degC can never be read as a 5 K tolerance |
| **interval** containing the true value | false confidence: a requirement is *met*, *failed* or honestly *undecided* |
| **rigor** (the emergence engine's rigorous / validated / asymptotic) | treating an estimate as a bound |
| **owner** | changes nobody can trace; every result names the teams its uncertainty comes from |

## The model

* **params**: a value, a tolerance (`"10%"`, `"0.1 mm"`, or a range) and an
  owner. A tolerance is a stated bound, so a parameter's interval is rigorous
  by construction.
* **materials**: one sourced library (Incropera Tables A.1/A.4/A.6, ASM
  Vol. 2, IEC 60028; estimates marked) shared by all teams, extendable per
  project. Temperature-dependent properties are tables. A query outside the
  table is an error, not an extrapolation.
* **derived**: unit-checked interval expressions (`ambient + q1.power *
  q1.rth_total`). A declared display unit is also checked against the
  expression.
* **analyses**: typed models with dimensioned inputs and outputs, declared
  monotonicity, and stated assumptions. Built in: `pcb_trace`,
  `pcb_microstrip`, `conduction_slab`, `plate_fin_heatsink`, and
  `chip_design` (Verilog synthesized onto transistor-characterised cells;
  see `analysis/chip.md` and `projects/mini-cpu`), and the HVAC kinds
  `zone_loads`, `refrigeration_cycle`, `duct_run`, `fan`, `cooling_coil`,
  `thermal_comfort` (`analysis/hvac.md`, `projects/office-hvac`).
* **parts**: parametric CAD models (`cad` JSON) driven by project
  parameters, giving volume, area, mass and envelope. They are exported as
  STL and rendered (`viz` sphere tracer) from the same values. A part with
  an `fdm` block also reports what it weighs as printed, the least time it
  can print in, and its thinnest wall.
* **requirements**: limits on any expression, owned by a team.
* **feedback loops**: a cycle through a derived node that declares a
  `range` is solved as a fixed point (below); any other cycle is refused,
  with its path named.

Analyses run at the corners of their toleranced inputs. When a model is
monotone in each toleranced input separately, its extremes over the
tolerance box lie on corners, so the enclosure is a bound. Each model
declares which inputs that holds for, and a test checks every declaration
against the model on random designs. When a toleranced input is not declared
monotone, the result is labelled asymptotic and the report says which input
caused it. Every model adds its own error claim on top. Runs are memoised by
exact input bits, so after an edit only what the edit reaches runs again,
and the result is bit-identical to a fresh evaluation.

## Commands

```sh
fluid project report projects/motor-driver/project.json --why --html report.html
fluid project check  project.json          # CI gate: 0 all met, 1 failed, 3 undecided/unavailable
fluid project diff   project.json --rev main
fluid project set    project.json sink.fin_count 7
fluid project export project.json sink_body sink.stl
fluid project render project.json sink_body sink.png
fluid project serve  project.json --port 8080   # live shared editor
fluid project kinds
```

`serve` turns the project into a shared web page. Every parameter is an
input. An edit is checked (unit, tolerance, a project that still loads),
evaluated incrementally, written to the project file, and pushed to every
viewer, with an activity feed of what it did to each team's requirements.
Edits made outside the editor (a text editor, a `git pull`) are picked up
the same way, so the file stays the source of truth for git and CI.
Requests are handled one at a time, so two edits can never interleave.
There is no login: it binds to loopback unless told otherwise.

`diff` is the design review. It lists changed definitions with their owners,
then every requirement that changed status, with the upstream changes that
caused it. When a requirement can no longer be judged, it gives the root
cause and the path that cause travelled. A real example from the motor
driver: mechanical cuts the fins from 9 to 6 and raises the duct velocity to
3.5 m/s.

```
  – REQ-TJ     met → unavailable   [systems] Q1 junction at or below 100 degC ...
      74.54 degC [63.96, 89.97] → unavailable
      why: channel Reynolds number 3161 exceeds 2300: turbulent channel flow has no rung yet
           [at nominal] (via sink.thermal_resistance) (via q1.rth_total) (via q1.tj)
      caused by: duct.velocity (mechanical), sink.fin_count (mechanical)
```

The tool refuses to give a laminar answer outside the laminar regime, and it
tells systems which mechanical change took their requirement out of reach.

## The first seam: electrical → thermal → fluids → mechanical

The motor-driver example joins six teams: MOSFET losses and a power trace
(electrical), a gap pad and a junction stack-up (thermal), a ducted
plate-fin heat sink (thermal and fluids), the sink's CAD geometry and mass
(mechanical), and an encoder line (signal integrity) on a stack-up the PCB
fabricator owns.

Result (`fluid project report --why`):

| requirement | value | status | uncertainty mostly from |
|---|---|---|---|
| Q1 junction ≤ 100 °C (electro-thermal loop) | 73.4 °C [62.7, 90.7] | met, 9.3 K worst-case margin | q1.switching_loss (electrical) 25 %, gap-pad k (thermal) 19 %, q1.rds_on_25 (electrical) 17 % |
| encoder line 50 Ω ± 10 % | 47.6 Ω [40.7, 53.6] | **undecided** | etch width (signal integrity) 25 %, prepreg thickness (fab) 19 %, permittivity (fab) 15 % |
| trace drop ≤ 50 mV | 46.4 mV [38.4, 56.4] | **undecided** | copper thickness tolerance (electrical) 53 % |
| sink mass ≤ 60 g | 45.3 g [43.3, 47.3] | met | fin thickness tolerance (mechanical) 61 % |
| sink Δp ≤ 20 Pa | 6.7 Pa [5.1, 8.7] | met | duct velocity (mechanical) 79 % |

"Undecided" is the honest answer for the trace. The nominal design passes,
but 2 oz copper at ±10 % can push the drop past the limit. The report names
the input to tighten and the team that owns it.

### The heat-sink model

The fins form channels with simultaneously developing laminar flow. The
analysis marches the parabolized (boundary-layer) equations between
parallel plates, the method of Bodoia & Osterle (1961): implicit
streamwise, with the pressure gradient set at each station by mass
conservation. The mean Nusselt number gives an exact ε-NTU heat balance
(`θ_b = exp(−NTU)`), 1-D fin efficiency, and base conduction. Its
discretisation error is the change from a grid twice as coarse.

| check | reference | result |
|---|---|---|
| Blasius wall shear f''(0) | 0.332057336 | 0.3320573362 |
| Pohlhausen θ'(0), Pr = 0.7 | 0.2927 | 0.29268 |
| Pr = 1: θ'(0) = f''(0) | exact identity | to 1e-8 |
| fully developed Nu_Dh, uniform wall T | 7.5407 (Shah & London 1978) | 7.5401 (8e-5) |
| fully developed dP/dX, mid-plane u | −3, 1.5 (Poiseuille) | −2.9999, 1.49995 |
| leading edge: Nu·√X / 4θ'(0), extrapolated to X → 0 | 1 | 1.0005 |
| example sink Nu_mean, fine vs coarse grid | — | 11.961 ± 0.012 |

### Feedback loops

The MOSFET's on-resistance rises with its junction temperature, and that
temperature depends on the resistance: `Rds(Tj) = Rds25 (1 + tc (Tj − 25
°C))`. Giving `q1.tj` a `range` makes the cycle a fixed-point problem
`X = F(X)`, solved by `X ← F(X) ∩ X`. Every fixed point inside the range
equals its own image, so it survives each step. If the intersection
becomes empty, there provably is no operating point in the range, which
is how thermal runaway is reported. When the final interval maps into
itself, an operating point provably exists in it (Brouwer). The nominal is
iterated from the bottom of the range, so an increasing loop settles on
its lowest operating point, the stable one.

| check | result |
|---|---|
| linear loop `t = base + gain·t`, base ∈ [10, 12] K, gain ∈ [0.2, 0.3] | [12.5, 120/7] K and nominal 44/3 K, all to 1e-11; rigorous; existence proven |
| loop gain 1.2 | "no operating point … (runaway, or the range is too small)" |
| motor-driver Tj loop | settles in 16 iterations; nominal equals the closed-form solution to 1e-9; existence proven |

### Spreading resistance

The MOSFET tab (10 × 15 mm) is smaller than the sink base (40 × 50 mm),
so heat spreads through the base. The analysis uses the exact cosine-series
solution for a centred uniform-flux rectangle on a plate with adiabatic
edges and a uniform h beneath; the fin array is smeared into an effective
h, `1/(R_conv · base area)`. Mode `(m, n)` responds at the heated face with
`G = (1 + ε tanh βt)/(kβ (tanh βt + ε))`, `ε = h/(kβ)`.

| check | result |
|---|---|
| source covering the whole plate | 0 (to 1e-15) |
| small square source → half-space, `R k c = (4/3)(1 − √2 + 3 asinh 1)/(2π)` = 0.473201 (Yovanovich: 0.4732) | extrapolated in c/L: 0.473180 (4e-5) |
| thin, strongly cooled plate → footprint's 1-D resistance | within 1.2 % (an O(t/c) edge effect) |
| example | 0.23 K/W, 13 % of the sink's resistance |

A thicker base conducts worse but spreads better, so `base_thickness` is
no longer declared monotone.

### Microstrip impedance (signal integrity)

`pcb_microstrip` is a two-rung ladder. Copper thinner than h/12 cannot be
drawn on the uniform grid, so it takes Hammerstad–Jensen (quoted 1 %) and
counts Wheeler's thickness correction as error rather than applying it.
Otherwise it runs `cad`'s quasi-TEM field solve of the drawn cross-section,
on a grid chosen to put the copper's thickness and half-width on grid
lines. The error is the change from halving the cell plus the change from
enlarging the box. It is monotone in width, thickness and permittivity by
Dirichlet's principle, and the test checks that.

| check | result |
|---|---|
| claim vs an independent solve, finer (h/36) and in a larger box | claim ±0.63 Ω, miss 0.38 Ω |
| 70 µm copper on 0.2 mm prepreg | field 47.6 Ω; zero-thickness formula 51.5 Ω (8 % high) |
| 35 µm on 1.6 mm FR-4 (thin rung) | the 50 Ω rule of thumb, inside the stated error |

**Wrapping the solver found a bug in it.** `cad::laplace` classed a node as
conductor when the region's field was `≤ 0`, but exactly on the region's
edge the field evaluated to +1.4e-20. A trace drawn on its dielectric lost
its whole bottom row of nodes and floated one cell above the substrate in
air. 70 µm copper read 64.1 Ω (ε_eff 2.14) at h/8 and crept down only as
the gap shrank. With a round-off tolerance the same solve gives 49.2,
47.3, 47.3 Ω at h/8, h/16, h/32. Fixed at the source, with a regression
test: growing a conductor by 1 pm must change nothing.

### FDM printing

Printed mass, filament length, least print time and thinnest wall all
come from one depth field. The part is voxelised on the kernel's exact
sign, then an exact Euclidean distance transform (Felzenszwalb &
Huttenlocher) runs in linear time. The kernel's own field could not be
used, because for a union it is only a lower bound on depth.

* Walls of `perimeters × line width` all round, infill inside. Voxels
  straddling the wall's inner surface count by the fraction beyond it.
* The least print time, printed volume / hotend flow, is a true lower bound:
  no printer extrudes faster than it melts.
* A wall is a crest of the depth field across an axis, a run of at most
  two voxels. Convex edges, where the medial axis runs into a corner, are
  not mistaken for walls. The claim is 3 cells: a voxel's depth is within
  (√3/2 + 1/2) cells of the truth, and the wall doubles that.

| check | result |
|---|---|
| 40 × 30 × 20 mm box, 0.8 mm walls, 20 % infill | 7.9442 cm³ vs 7.9470 exact (claim ±0.0137) |
| 1.2 mm plate on a block | 1.200 mm |
| 1 mm tube wall (curved) | 0.69 mm: thin, inside the claim, erring toward flagging |
| example duct (2 mm walls, 3 × 0.4 mm perimeters) | prints solid, so printed mass = PLA mass; ~30 min least time; thinnest wall 1.81 mm, short at the reflex inner corners |

The duct is sized from the sink's parameters, so it follows mechanical's
changes, and REQ-PRINT (owned by the FDM team) checks its walls against
two extrusion lines.

### Other checks

| check | result |
|---|---|
| CAD part volume vs exact `W L t_b + N t_f H L` = 16.80 cm³ | 16.76 cm³ (−0.24 %), exact value inside the claimed interval |
| 5000 random points in the trace's tolerance box | all inside the rigorous interval |
| declared monotone inputs (heat sink, trace) | checked on random designs: no reversal |
| edit one parameter | unaffected analyses reused; result identical to a fresh run, bit for bit |
| contribution of `x` (±10 %) to `x + y` with `y` ±1 % | 10/11 exactly |
| air table (Incropera A.4) | ν = μ/ρ and Pr = μc_p/k hold to 0.2–0.3 % at every row (catches transcription errors) |
| mistakes refused | temperature bound to a length, misspelt input, wrong declared unit, dependency cycle (named), unknown id, limit in the wrong unit, Celsius tolerance, duplicate ids |

Timing for the example on the i7-13620H, with tolerance corners, field
solves and per-input re-evaluations run in parallel: a full evaluation
takes 2.9 s (it was 12.6 s serially) and `report --why` takes 7.6 s (was
33 s). An edit in the live editor reruns only what it reaches.

## Honest limits

* **Laminar, ducted heat sinks only.** Above Re_Dh = 2300 the analysis
  refuses. There is no turbulent rung, no bypass flow, and no entrance or
  exit pressure losses; spreading treats the fins as a uniform h. Each is
  listed under the model's assumptions in every report. The error claim
  covers the numerics within those assumptions and nothing outside them.
* **Loops are solved by contraction.** A loop whose gain is near 1
  converges slowly (400 iterations at most), and its interval
  iteration does not track correlation, so a strongly coupled loop's
  enclosure can be wide.
* **Corner enumeration costs 2^k runs.** Up to ten toleranced inputs are
  enumerated exhaustively. Beyond that, the ends are sampled one input at a
  time and the result is labelled asymptotic. Intervals do not track
  correlation between inputs (`x − x` is not 0), so they are conservative.
* **Contributions do not add up to 100 %.** Each is the share of the width
  that vanishes when that one input is pinned, which is the question a team
  asks ("what if we tighten this?"). It is not an additive budget.
* **Collaboration has no accounts.** The live editor trusts whoever can
  reach it. Ownership is advisory (editing another team's parameter asks
  for confirmation), and the activity feed lives only as long as the
  server. The durable record is the project file's git history.
* **CAD parts give mass properties only.** Faces have no persistent names,
  so a boundary condition cannot yet attach to "this face" of a part. That
  is the prerequisite for driving the flow solvers from the project.
* **FDM is geometric.** Mass, time and wall checks only: no warping,
  residual stress or anisotropic strength, and support material and
  overhangs are not considered.
* **Five analysis kinds plus parts.** The optics and flow solvers
  elsewhere in the workspace are not yet exposed as project analyses. Each
  needs a thin wrapper that declares ports, monotonicity and assumptions.
