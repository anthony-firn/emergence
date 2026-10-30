# Cross-domain CAD: any project, eight domains, one solve

`/play/emergence/cad/` is a design environment, not a design. A project is a
file that names the domains it uses, carries one design document per domain,
the typed, unit-checked links between their ports, and what to watch. One
solver runs the project's domains in dependency order and iterates the loops
between them until the linked values stop changing, or says that they did not.
Nothing in the environment knows what is being designed.

## Projects

A project file (`emergence-project/1`, `core/project.js`) is JSON:

```
{ "format": "emergence-project/1",
  "documents": {
    "project": { name, summary, domains, headline, studies, optimize, tour, scene, credits },
    "system":  { links, confirmed },
    "flight" | "parts" | "pcb" | "silicon" | "litho" | "optics" | "fluids" | "architecture": { … },
    "mech":    { params, features, drawings } } }
```

| manifest field | what it does |
|---|---|
| `domains` | which studios the project uses; the others keep their documents but are not solved, and links to them go dormant |
| `headline` | the outputs the System studio leads with |
| `studies` | numeric fields, by path (`masses.gimbal.mass`: list items by id), raised 5% for the coupled sensitivity study |
| `optimize` | objective port and sense, constraint ports, variables by path with range and step; a compass search where feasible beats infeasible |
| `tour` | steps that run real commands and report the outputs and checks they name |
| `scene` | where each domain's geometry sits in the 3D view (`at`, `rx`, `rz`, `count` copies spaced by `phase`; `count` may read another document, such as `flight:rotors`) |
| `credits` | the public reference designs the documents come from |

Every field is edited by commands (`project.add_domain`, `project.pin`,
`project.study`, `project.set_value`, …) and versioned like any other document.
A domain's commands are refused while it is not in the project. In the browser
each project is its own revision history on the device: opening an example or a
file makes a first revision with a fixed author and time, so an example makes the
same revision everywhere and peers who open it share history. Exports, peer
channels and imports carry the project's id; another project's history is
refused rather than merged. The CLI and MCP server read and write the same files.

### Example projects (`web/cad/projects/`)

**Inspection quadcopter** (`inspection-drone.json`): vehicle scale, all eight domains.

```
             parts.mass ─┐            optics.mass ─┐
 mech.sink_mass ─┐       ├──► flight ──► esc_current ──► mech.current ──► Q1 losses
 pcb.board_mass ─┴───────┘     │  │                         │
                               │  └ wake_pressure ──► fluids ──► duct_velocity ──► mech.velocity
                               └ altitude ──► optics (ground sample distance)
 architecture.outdoor ──► mech.ambient, fluids, silicon, flight air density
 silicon.power ──► pcb U1     pcb.copper_loss ──► mech.heat     flight.charger_heat ──► hangar
```

The weight spiral (mass → thrust → ESC current → heat → heat-sink size → mass)
is a genuine coupled loop. The reference design settles in four passes from a
cold start and two to three after an edit (warm start, edited domain first):
8.27 kg take-off, 24.6 min hover endurance, Q1 at 51 °C, 0.414 cm ground sample
distance at 60 m. Its optimisation finds the lightest ESC heat sink that keeps Q1
at or below 100 °C.

**Edge compute shelter** (`edge-compute-shelter.json`): building scale, no flight or optics.

```
 silicon.power ─┐
 mech.q1_power ─┼──► architecture.equipment_power (per node × nodes) ──► room_temperature
 pcb.dissipation┘                                                          │
        ▲                                                                  ▼
 die leakage, VRM losses, fan duct  ◄────────────── silicon.ambient, mech.ambient, fluids.air_temperature
```

Accelerator nodes (a 48 mm² die whose control CPU and fan PWM blocks are
synthesised from public RTL, a VRM power stage on the node board, its heat sink
in a 60 mm fan duct) sit in a BESTEST 600 shelter. Each node's heat goes into the
shelter and the shelter's air is every node's inlet, so the loop runs through five
domains. With 24 nodes the shelter air is 29.3 °C and every check is met or
honestly undecided; the project's optimisation ("Fill the shelter") finds 127
nodes before the air reaches its 35 °C limit. The node board ships unrouted: at
25 A a plain trace fails the current-density rule (real VRMs use copper pours,
which this board model does not have), and the PCB studio says so once routed.

A **blank project** has no domains; adding one opens its studio with a starter
design (a motor-driver board, a five-block SoC, the Cooke triplet, a fan duct,
one BESTEST room, a quadcopter, a bolted flange).

## Architecture

| layer | file | what it guarantees |
|---|---|---|
| projects | `web/cad/core/project.js`, `web/cad/projects.js` | a project is a validated file; field paths address any numeric field; each project is its own history on the device |
| units | `web/cad/core/units.js` | SI dimension vectors; `°C`↔`K` with offset; links refuse mismatched dimensions |
| commands | `web/cad/core/bus.js` | every edit is a registered command with typed, ranged arguments; one undo stack across domains; slider bursts coalesce; scripts are atomic across domains |
| store | `web/cad/core/store.js` | a failed edit never leaves a partial document: commands edit a clone, validators project it |
| links and solve | `web/cad/core/system.js` | typed ports; topological order; Gauss–Seidel with warm start and relaxation on divergence; range refusal; residual history |
| merge | `web/cad/core/merge.js` | three-way merge by object id, so different components, rooms or lenses edited on two branches merge cleanly |
| revisions | `web/cad/repository.js` | each domain document is a content-addressed blob a revision refers to; legacy revision hashes are unchanged |
| integrity | `web/cad/core/integrity.js` | exports carry SHA-256 of their canonical JSON and are read back before download |
| domains | `web/cad/domains/*.js` | each domain's model, validator, ports, checks and live certificates |

The browser studios, the command line, the guided tour, scripts,
`window.emergence`, the CLI (`web/cad/agent/cli.mjs`) and the MCP server
(`web/cad/agent/mcp.mjs`) all call the same registry. The CLI and MCP server
load the same WebAssembly as the browser, so an agent gets the same physics.

## The kernels behind each domain

| domain | model | source |
|---|---|---|
| Mechanical | the motor-driver electrothermal project graph, heat-sink mesh | Rust `project`, `cad` (unchanged) |
| Parts | JSON feature trees: primitives, booleans, shells, arrays and every kit | Rust `cad::model` in WebAssembly, file features refused |
| PCB | trace resistance from drawn length and section, Gerber RS-274X, Excellon | Rust `cad::kits::pcb` |
| PCB | microstrip Z₀: Hammerstad–Jensen and the quasi-TEM field solve of the drawn section | Rust `cad::laplace` |
| Silicon | Verilog → synthesis → equivalence → nine-corner STA and power | Rust `chip` with pre-characterised corner libraries |
| Silicon | Abbe imaging, threshold resist, inverse lithography | `domains/litho.js` |
| Optics | 3D real-ray trace, Sellmeier glasses, unpolarized Fresnel, paraxial ABCD | `domains/optics.js` |
| Fluids | Shah–London laminar, Colebrook turbulent, bounded transition, operating point | `domains/fluids.js` |
| Architecture | steady room heat balance, ASHRAE 62.1 outdoor-air rate | `domains/architecture.js` |
| Flight | actuator-disc momentum theory, pack sag, endurance | `domains/flight.js` |

## Public reference designs

| domain | design | licence | adaptation |
|---|---|---|---|
| Silicon | `tt_um_hadjal001_motor_controller`, [HADJAL001/tt-motor-controller](https://github.com/HADJAL001/tt-motor-controller) at `0f078f5`, a Tiny Tapeout dual-channel PWM motor controller | Apache-2.0 | asynchronous reset made synchronous, lint-only lines removed; every change is listed in the file header |
| Silicon | `mini8`, `fluid/projects/mini-cpu/mini8.v` | Apache-2.0 | none |
| Optics | Cooke triplet from the [Optiland gallery](https://optiland.readthedocs.io/en/latest/gallery/basic_lenses/Cooke_Triplet.html), after Kingslake | MIT | the stop is a plane just past surface 4; clear diameters chosen here |
| Architecture | ANSI/ASHRAE Standard 140 BESTEST case 600 geometry and envelope for the charging room | published test case | solved as a ventilated steady-state room, not the standard's annual test; U-values to be confirmed against the standard |
| Parts | the repository's flange, gear, bracket, elbow, pin-fin, gyroid and wing models | Apache-2.0 | rescaled to millimetres |

## What each number was checked against

`node --test web/test_cad_domains.mjs` (18 tests, including every example file's
manifest, links and ports), the Rust adapter tests in `emergence-web` (9), and the
browser suites `web/test_cad_{browser,studios_browser,system_browser,manufacturing,projects_browser}.py`.

| check | reference | result |
|---|---|---|
| N-BK7, N-SK16, F2 at 587.56 nm | Schott catalogue n_d | 1.51680, 1.62041, 1.62004, within 1e-5 |
| single-lens paraxial EFL | thick-lens lensmaker equation (as `cad::kits::optic::Lens::focal_length`) | 1e-10 relative |
| flat window focus shift | t(1 − 1/n) | 1e-9 mm |
| normal-incidence Fresnel | 1 − ((n−1)/(n+1))² | 1e-15 |
| Cooke triplet | Optiland: f = 50 mm at f/5; image plane 42.208 mm behind the last vertex | EFL 50.000 mm; least-RMS focus within 0.03 mm of the prescription's image plane; 20° image height within 0.4% of f·tanθ (the design's distortion) |
| Shah–London fRe | 96 (parallel plates), 56.91 (square duct) | exact, 0.02 |
| Colebrook | Moody value at Re 1e5, ε/D 1e-4 | 0.01852 ± 1e-4 |
| fan / system balance | the operating point | residual below 1e-6 of the pressure |
| supplied-room heat balance | supply removes internal + solar + envelope | 1e-9 W |
| actuator-disc hover power | T^1.5 / √(2ρA) | 1e-9 relative |
| Abbe imaging, clear mask | intensity exactly 1 | 1e-12 |
| coherent 80 nm-pitch grating past cutoff | DC term squared, 0.25 | 1e-12 |
| inverse lithography | printed-pattern error falls | 192 → fewer pixels in 24–36 adjoint steps |
| GDSII | every block boundary read back | exact to 1 nm; a 1 nm move is detected |
| Gerber read-back | drawn copper rasterised at 0.1 mm | 100.0% agreement on the reference board |
| trace resistance (kernel) | ρL/(wt) with the copper temperature coefficient | 1e-12 relative |
| microstrip field solve | Hammerstad–Jensen, W/h 1.8 | within 5% (the field solve includes copper thickness) |
| browser chip flow | native `chip::corners::at` | min period identical to 1e-12 |
| drone example | coupled loop | settles; ESC current and per-ESC mass carried exactly |
| shelter example | coupled loop through the building | settles; inlet air and per-node heat equal their sources to the solver tolerance; more nodes warm the shelter and every die |
| project files | every example | manifest validates, every link checks, every headline, study, variable and tour command resolves; a saved file reopens identically |

Every domain also reports live certificates in the UI: the lensmaker identity
for the current lens, the closed room heat balance, the balanced operating
point, the actuator-disc power identity, kernel trace resistance, intervals
containing their nominal and netlist ≡ RTL.

## Honest refusals the system keeps

* Fin channels above Re 2300 have no rung in the mechanical graph: Q1 is
  reported **unavailable**, not extrapolated. Above the mechanical model's
  validated 0.5–4 m/s approach velocity the fluids link is **refused** and the
  mechanical design keeps its own value; the System studio says so.
* Ducts between Re 2300 and 4000 are bounded by both friction laws and the flow
  is reported as an interval.
* A chip junction outside the 250–425 K characterised corners stops the chip
  flow instead of extrapolating the cell tables.
* Links between incompatible dimensions cannot be created.

## A finding in the native flow

Exporting the nine corner libraries showed that the 1.08 V / 425 K corner has
one table entry that failed characterisation: the full adder's fall output slew
at the smallest input slew and load is NaN. JSON writes NaN as `null`, so the
native disk cache for that corner never reloads and the corner is silently
re-characterised on every run. The browser carries the entry as NaN, exactly as
the native flow holds it; the characterisation itself is not changed here.

## Assumptions every answer inherits

Each studio lists its own, and the System studio counts how many the user has
confirmed against their own data. Among them: the rotor figure of merit and
motor and ESC efficiencies; the fan or wake source; fitting K values; nominal
FR-4 permittivity; illustrative footprints and drills; planning block powers and
the package θJA; Sellmeier glass data and air n = 1; climate, façade irradiance
factors and envelope U-values; handbook material densities.

## Inspiration

Three of the open CAD projects studied — [CADCraft](https://github.com/storytold/cadcraft),
[SolidCraft](https://github.com/tangency-labs/solidcraft) and
[Open CAD Studio](https://github.com/HakanSeven12/OpenCADStudio) — share one
command set across UI, CLI and MCP. [Purrametric](https://purrametric.deepmagic.de/)
shows a fully-defined-sketch view, which became the confirmed-assumptions view here.
Open CAD Studio's verified saves became read-back exports. [ZiroEDA](https://github.com/ZiroEDA/ziro-designer)
suggested KiCad interop and a Gerber read-back. [OpenROAD-flow-scripts](https://github.com/The-OpenROAD-Project/OpenROAD-flow-scripts),
[KLayout](https://github.com/KLayout/klayout), [OpenIC](https://github.com/AhmedLilah/OpenIC)
and [vibeic](https://github.com/vibeic/vibeic-eda) shaped the staged silicon flow, GDSII
and the rule never to claim sign-off. [lithosim](https://github.com/VLSIDA/lithosim) and
[OpenILT](https://github.com/OpenOPC/OpenILT) shaped the lithography stage,
[Optiland](https://github.com/optiland/optiland) the optics analyses, and
[OpenMDAO](https://openmdao.org/) and preCICE the coupled solve, sensitivities and
the weight-spiral optimiser.

## Limits and next steps

* No KiCad *import* yet: the PCB model's footprints are five illustrative
  through-hole packages, so a public KiCad board cannot be represented faithfully.
  The `.kicad_pcb` export is written to the documented format but has not been
  opened by `kicad-cli`.
* No STEP import or export in the browser. OCCT's import build is vendored but
  not wired in.
* The renderer draws kernel meshes; it does not yet sphere-trace the distance
  field directly, which the 1-Lipschitz invariant would allow.
* Optics are monochromatic and uncoated, with no tolerancing; lithography is
  scalar and in focus; the building is steady state.
* Sensitivities are finite differences through the full solve (one solve per
  study input), and the optimiser is a compass search, which finds a local
  optimum over a few variables.
* The 3D scene places each domain's geometry where the project file says; it
  does not check that placed geometry clears.
