# Cross-domain CAD: any project, one solve

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
| `requirements` | your own limits on any output, checked on every solve; an interval that straddles a limit reads undecided |
| `studies` | numeric fields, by path (`masses.gimbal.mass`: list items by id), raised 5% for the coupled sensitivity study |
| `optimize` | objective port and sense, constraint ports, variables by path with range and step; a compass search where feasible beats infeasible |
| `tour` | steps that run real commands and report the outputs and checks they name |
| `scene` | where each domain's geometry sits in the 3D view (`at`, `rx`, `rz`, `count` copies spaced by `phase`; `count` may read another document, such as `flight:rotors`) |
| `credits` | the public reference designs the documents come from |

Every field is edited by commands (`project.add_domain`, `project.pin`,
`project.require`, `project.study`, `project.optimize`, `project.tour_step`,
`project.place`, `project.credit`, …), each with a visible control, and versioned
like any other document.
A domain's commands are refused while it is not in the project. In the browser
each project is its own revision history on the device: opening an example or a
file makes a first revision with a fixed author and time, so an example makes the
same revision everywhere and peers who open it share history. Exports, peer
channels and imports carry the project's id; another project's history is
refused rather than merged. The CLI and MCP server read and write the same files.

### Example projects (`web/cad/projects/`)

**Inspection quadcopter** (`inspection-drone.json`): vehicle scale, twelve domains,
built from parts you can buy and every value traceable to its source.

```
             parts.mass ─┐            optics.mass ─┐
 mech.sink_mass ─┐       ├──► flight ──► phase current, bus, duty ──► power ──► FET RMS, switching ──► mech ──► Q1 Tj ─┐
 pcb.board_mass ─┴───────┘     │  │ │                                                                  ▲               │
                               │  │ └ thrust, torque, rpm ──► structures (stress, mode against blade-pass)└───────────────┘
                               │  └ wake_pressure ──► fluids ──► duct_velocity ──► mech.velocity
                               ├ altitude ──► optics (ground sample distance)
                               └ mass, inertia, thrust, motor lag ──► controls (PX4 cascade, 6-DoF)
 architecture.outdoor ──► mech.ambient, fluids, silicon, flight air density
 pcb.copper_loss ──► mech.heat     flight.charger_heat ──► hangar     every output ──► testing
```

What is real in it, and where it came from:

| part | what the model uses | source |
|---|---|---|
| Motor and propeller | T-MOTOR MN5008 KV170 with P18×6.1: the 11-row bench table (throttle, volts, amps, rpm, torque, grams), 128 g, 270 mΩ, 15 A | authorised distributor's page, which marks it for reference only |
| Propeller mass | 31.5 g per blade | same distributor; other listings say 13 g, so the heavier figure is used |
| Battery | Tattu Plus 1.0 Compact 16000 mAh 12S 15C, 4250 g, with a LiPo voltage curve | maker's page |
| Autopilot | Holybro Pixhawk 6X on a mini baseboard, 57.8 g | retailer listing |
| MOSFET | Infineon BSC040N10NS5: RDS(on), gate charges, Qoss, Qrr, land pattern | data sheet Rev. 2.4 |
| Gate driver | TI UCC27211 SOIC-8 pinout | data sheet |
| Bus capacitors | three Rubycon 63ZLJ470M12.5X20 (470 µF, 63 V, 33 mΩ, 2.4 A at 100 kHz) on the battery leads | distributor listing, which does not give the ESR's frequency; 4 g each is an estimate |
| ESC board | 78 × 42 mm, two layers, 2 oz: six MOSFETs, three drivers, an STM32F051K8T6, an LM5017 buck and SRR1260 inductor, LP2985 regulator, INA240 on a 1 mΩ shunt, phase and bus dividers, poured bus and phase copper, ground plane | laid out here by commands; pinouts and land patterns from TI, Infineon and KiCad's libraries |
| ESC firmware | AM32's F0 group B pin map, with this board's own target: 18:1 bus divider, 50 mV/A, 500 ns dead time | [am32-target.h](quadcopter-evidence/am32-target.h); written and pin-checked, **not yet compiled** |
| Arm | 35 × 2 mm carbon tube, 240 mm free span from the plate corner to the motor | supplier's modulus and strength for one layup |
| Controller | PX4's cascade, parameter names and defaults, then gains tuned to this airframe; exports a QGroundControl parameter file | PX4 parameter reference |

The gimbal, wiring, gear and fasteners are still estimates, and the file says so
line by line. The flight controller is the Pixhawk, not a die: the silicon domain
stays in the project as a study of moving motor control onto a small chip, not
as a part on the aircraft.

The weight spiral (mass → thrust → motor current → MOSFET heat → heat-sink size
→ mass) is a genuine coupled loop, and it settles in five passes from a cold start:
7.35 kg take-off, 39.7 min hover to 20% reserve, 61% hover throttle at 4200 rpm,
thrust to weight 2.04 full and 1.81 at end of charge, 11.0 A design phase current,
7.1 A RMS in the hottest MOSFET, 1.45 W inverter loss, 0.18 V bus ripple, Q1 at 45 °C. The optimisation
finds the lightest heat sink; the shipped sink is the lightest the laminar channel
model can still vouch for, and anything smaller is refused rather than extrapolated.

Realism found things the old efficiencies hid. PX4's stock gains leave this airframe
about 40° of phase margin, under the 45° required, so the file carries tuned gains (59°,
8.1 dB, 5.6 Hz at 400 Hz). A 30 mm tube passes on strength but puts the arm's first mode in the
rotors' blade-pass band, so the arm is 35 mm (259 Hz, against a band that ends at 229 Hz
with 15% separation). A quadcopter does not survive losing a rotor, and the simulation shows that.

Testing keeps four kinds of evidence apart. Thirteen simulation items pass; they are
predictions. One correlation item fits T = c·n² to the KV170 table and predicts
the same propeller on the KV340 motor's published data, which it never saw: 1.3% RMS,
2.2% worst, against a 5% limit set before the result was known. Ten real-software results
are recorded from PX4 and KiCad, and each goes stale if a document it covers changes. Nine
physical tests (weigh it, thrust stand on the arm, double pulse, thermal soak, tap test, motor
step, tethered vibration, roll step in flight, measured endurance) stay "not run" until a
person records a measurement, as does compiling the firmware. **Nothing in this project has
been built, and it has only ever flown on a simulated aircraft.**

#### What real tools said

Real PX4 v1.17.0 firmware, built from source, flew this aircraft's numbers in its SIH
simulator with the shipped gains, against this project's own simulator on the same plant at
the same 250 Hz loop:

| | PX4 | model |
|---|---|---|
| hover command | 0.419 | 0.419 |
| roll mix gain | 0.707 | 0.707 |
| 15° roll step: rise / settle | 0.275 s / 0.92 s | 0.275 s / 0.88 s |
| rate step overshoot, 24 steps averaged | 27.2% | 25.6% |
| stock gains: rate step overshoot | 7.4% and 14.4%, two flights | 15.0% |
| crossover | 5.90 Hz | 5.65 Hz |
| phase margin | 64.9° | 57.1° |
| gain margin (coherence 0.24, weak) | 9.3 dB | 7.7 dB |
| stock gains: phase margin | 43.9° | 38.7° |

All 29 parameters exist in PX4 with the same defaults, and the tuned gains are in range.
KiCad 9.0.9 reports 0 errors and 0 unconnected items; its only warnings say the footprint
library is not installed. The app's Gerbers overlap the copper KiCad computes by 99.8% on top
and 99.5% on the bottom on a 50 µm grid. All 22 pin assignments match AM32's source. The
scripts are in [quadcopter-evidence](quadcopter-evidence/README.md).

These tools corrected the work, which is the point of running them:

* KiCad found the bus pour cut in two, so phases B and C were not on the battery. The app had
  called the board fully connected because it treated a pour as one piece. Connectivity now
  follows the copper a pour really fills, and the router avoids other nets' pours.
* KiCad found courtyard overlaps; the app now checks them, and writes real vias.
* The flight log showed the model carried about 6 ms of delay that PX4 does not have. Removing
  it moved phase margin 8° closer on the same flights. Overshoot agrees to within the scatter
  between repeat flights, which is several points even after averaging 24 steps.
* AM32's reference divider would put 4.6 V on a 3.3 V ADC pin at 12S, and its dead time is
  ten times what the power model assumed. The board has its own divider, and the model now
  uses the firmware's 500 ns.
* A long-standing tuner claim failed: it read margins before adding the integrator.

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

## Building from nothing

A project starts empty: no domains, no links. Each domain can be added empty or
with a starter design to edit, and can be emptied again. What makes that claim
checkable is `domains/rebuild.js`, which writes any project as the ordered
commands that build it from a blank one. A field no command can set could not be
reproduced, so replaying the script is the proof:

| check | result |
|---|---|
| inspection quadcopter, 173 commands from blank | every document identical; all 53 outputs equal exactly; 29 checks equal |
| edge compute shelter, 133 commands from blank | every document identical; all 35 outputs equal exactly; 20 checks equal |
| every domain added empty | solves with no false alarm |
| every command | has a visible control (tested against the UI source) |
| a new design through visible controls alone (`web/test_cad_scratch_browser.py`) | a housing from a disc, a custom footprint, a lens, a fin duct, a room, a formula, a link, a requirement, a study, an optimisation, a tour step, a credit and a 3D layout; saved and reloaded identically |

The build script of any project is one button in the System studio and
`cli.mjs script --project f`, so a design can be read step by step.

The tools for a design that is not an example:

* **Parts**: start from a block or disc, then `parts.solid` adds, cuts or keeps a
  box, cylinder, sphere, cone, torus, or a sketched extrude or revolve;
  `parts.modify` shells, offsets, rounds, patterns, mirrors, moves, rotates, scales
  or clips what exists; `parts.define` makes parameters, and any dimension may be
  an expression of them. Each step wraps the tree, so it reads as the order you
  worked in.
* **PCB**: `pcb.footprint` makes a through-hole footprint from a package's numbers
  (one row, two rows, four sides) or a pad list, with its own pad and drill.
  `NC` marks a pin deliberately left unconnected. Every component has its own
  power port (`pcb.d1_power`), and `pcb.trace` draws copper through your own points.
* **Equations**: named inputs and formula outputs, parsed by hand and never
  executed. Inputs can be driven by links and formulas are ordinary ports, so a
  rail current from a die power and a voltage, a cost or a margin can join the solve.
* **Requirements**: `project.require` puts a limit on any output of any domain.

### What this does not make possible

The system can now author everything in its example projects, and a new design
within what its domains can express. It is a design and analysis environment, and
each domain is one kind of model:

| domain | what it can express | what it cannot |
|---|---|---|
| Flight | multirotors of 3, 4, 6 or 8 rotors, on efficiencies or a measured bench table | fixed wing, VTOL, forward flight |
| Mechanical | one power stage: a MOSFET on a fin heat sink | any other electrothermal circuit; stress, vibration or fatigue anywhere in the system |
| Structures | a cantilevered tube arm: von Mises stress, deflection, first mode against rotor excitation | joints, laminate analysis, fatigue, impact, a whole-frame model |
| Power | a three-phase MOSFET inverter from data sheet values, loss by loss, and its bus capacitor bank | waveforms, ringing, EMI, current sensing, lead inductance |
| Parts | up to 24 solids, each meshed | assemblies with constraints; STEP or other exact geometry out |
| PCB | two layers, through-hole and surface pads, pours, vias, 256 components, 512 traces | schematic capture, more layers, thermal relief, creepage rules |
| Silicon | a Verilog subset, 150 nm study cells, a block floorplan | place and route, DRC/LVS, any foundry PDK |
| Optics | rotationally symmetric spherical lenses, one wavelength | mirrors, aspheres, coatings, tolerancing, MTF |
| Fluids | one series air path | branches, liquids, CFD |
| Architecture | steady room heat balances, one equipment room | thermal mass, annual energy, structure |
| Controls | PX4's position, attitude and rate cascade in a 6-DoF simulation, with loop margins and flight tests | PX4 itself: no estimator, no firmware in the loop, no rotor flapping or vibration |
| Testing | a plan separating simulated, correlated, real-software and physically measured evidence, with results that go stale | it records results; it does not measure anything |

So the quadcopter example is a sourced design with real software checks behind it, not a
build package for an aircraft that flies. The ESC board is complete on paper, but its firmware
has not been compiled, nothing has been fabricated, and there is no schematic, wiring harness or
frame drawing. Two copper layers put more inductance in the hot loop than a production ESC would
accept; only the double-pulse test will say how much. The 0.2 mm clearance relies on solder mask
at 50 V. Estimated masses, rotor inertia, the 1.5 safety factor and the rigid arm root all need
the physical tests in the plan, and PX4's SIH shares none of a real aircraft's vibration,
rotor aerodynamics or estimator error.

Revisions outlive schema growth: a stored document must still validate, and
validation may add fields it has since gained but never drop or change one, so
histories saved by earlier versions open unchanged and keep their hashes.

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
| Flight | a measured bench table with rotor similarity, or actuator-disc momentum theory; LiPo discharge under sag; motor constants | `domains/flight.js` |
| Structures | Euler–Bernoulli bending, Saint-Venant torsion, von Mises, Rayleigh's first mode, a Campbell check | `domains/structures.js` |
| Power | block commutation with synchronous rectification: conduction, switching from gate charge, Qoss, Qrr, gate and dead-time loss; input ripple in the bus capacitors | `domains/power.js` |
| Controls | PX4's cascade at the loop rate, Butterworth gyro and D-term filters, mixer with yaw desaturation, quaternion 6-DoF, frequency-domain margins, a tuner | `domains/controls.js` |
| Testing | rotor similarity fitted to one table and asked to predict another; plan status by evidence kind | `domains/testing.js` |

## Public reference designs

| domain | design | licence | adaptation |
|---|---|---|---|
| Silicon | `tt_um_hadjal001_motor_controller`, [HADJAL001/tt-motor-controller](https://github.com/HADJAL001/tt-motor-controller) at `0f078f5`, a Tiny Tapeout dual-channel PWM motor controller | Apache-2.0 | asynchronous reset made synchronous, lint-only lines removed; every change is listed in the file header |
| Silicon | `mini8`, `fluid/projects/mini-cpu/mini8.v` | Apache-2.0 | none |
| Optics | Cooke triplet from the [Optiland gallery](https://optiland.readthedocs.io/en/latest/gallery/basic_lenses/Cooke_Triplet.html), after Kingslake | MIT | the stop is a plane just past surface 4; clear diameters chosen here |
| Architecture | ANSI/ASHRAE Standard 140 BESTEST case 600 geometry and envelope for the charging room | published test case | solved as a ventilated steady-state room, not the standard's annual test; U-values to be confirmed against the standard |
| Parts | the repository's flange, gear, bracket, elbow, pin-fin, gyroid and wing models | Apache-2.0 | rescaled to millimetres |

## What each number was checked against

`node --test web/test_cad_domains.mjs` (38 tests, including every example file's
manifest, links and ports and every from-nothing rebuild), the Rust adapter tests
in `emergence-web` (9, one more ignored), and the browser suites
`web/test_cad_{browser,studios_browser,system_browser,manufacturing,projects_browser,scratch_browser}.py`.

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
| published bench table | power and momentum: shaft below electrical, implied figure of merit 0.6–0.85, speed below no-load, on every row | all 11 rows hold |
| propulsion lookup | a row reads back; at one rpm T, Q, P ∝ ρ; below the table T ∝ n², P ∝ n³; no thrust past it | 1e-9 |
| hover endurance | pack watt-hours against minutes × hover watts | within 5% |
| connectivity | a trace across a pour opens the net, and the check says why | exact, on a test board and on the shipped ESC |
| evidence | a recorded result goes stale when a document it covers changes | a heavier camera stales the flights, not the board |
| inverter currents | ⅔ I² for a 120° block, and six devices summing to two in series | 1e-9 |
| bus ripple | I/2 at half duty; charge per cycle plus the ESR step | 1e-12; two capacitors overrun their rating, three do not |
| arm deflection and mode | an independent Euler–Bernoulli finite-element beam in the test | deflection 1e-6; Rayleigh's mode within 1% |
| flight simulator | hover trim, mix normalisation, and \|L\| of a torque-free rigid body | \|L\| within 0.5% over 2.5 s |
| rate-loop margins | an integrator (K·G crossover, 90°), then one with a known lag | crossover 0.2%, phase 0.5° |
| thrust model | the KV340 motor's published data, not used in the fit | 1.3% RMS, 2.2% worst |
| ESC board | UCC27211 and TDSON-8 pinouts, Gerber read-back on both layers, drilled pads, KiCad zones | exact pin nets; 97.4% top and 99.5% bottom at 0.1 mm |
| physical tests | the shipped file | none recorded; every one reads "not run" |
| drone example | coupled loop | settles; MOSFET current and per-ESC mass carried exactly |
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

* No KiCad *import* yet, though surface pads and pours mean a public board could now
  be represented. The export is checked by KiCad, but only through `quadcopter-evidence`.
* Nothing has been physically tested. The bench table is the maker's, copied from a
  distributor; no result in this project has been measured here.
* AM32 for this board has not been compiled, and PX4 has run only in SIH. The browser's own
  simulator is still a model of PX4, now checked against it. An in-flight step response is
  the next test of the controls.
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
