# Real-world geometry: CAD in, engineering answers out

The spectral solver works in a periodic box with immersed walls. Engineering
work needs CAD import, real boundary conditions (inlets, outlets, walls) and
post-processing. This document covers the new path from a CAD file to forces,
pressure distributions and ParaView output. The new pieces are two crates,
`fluid/crates/geometry` and `fluid/crates/immersed`, and the `fluid case`
command.

```
STEP / IGES / BREP ──gmsh+OpenCASCADE──┐
STL (binary/ASCII) / OBJ ──────────────┴─> TriMesh ─> checks/repair ─> exact signed distance (BVH)
                                                                               │
JSON case: fluid, box, per-face BCs, bodies, run, outputs ─────────────────────┤
                                                                               v
             Cartesian MAC grid, sharp immersed boundary (no meshing), exact FFT pressure solve
                                                                               │
          forces (momentum balance + surface integral), Cp, wall shear, probes, VTK, CSV, JSON
```

## Why an immersed Cartesian method

* **No meshing step.** The CAD surface is used as given. Its signed
  distance and normals drive the boundary condition, so the same run works
  for a cylinder or for a 900k-triangle tessellation. Commercial tools charge
  separately for the meshing that this skips.
* **The pressure equation stays the box Laplacian.** It is inverted exactly
  (to 3e-15) by one real trigonometric transform per axis, for any grid size
  and any mix of wall, inlet, outlet and periodic faces.
* **Costs.** The grid is uniform, so boundary layers must be resolved by
  global refinement; local refinement is not built yet. The cut cells need
  care for mass conservation (below).

## Pieces and how each was checked

### Geometry (`geometry` crate)

| piece | check | result |
|---|---|---|
| STL (binary/ASCII), OBJ (polygons, `v/vt/vn`), binary STL out | round trip; ASCII and OBJ parsing | watertight after weld; f32 volume change < 1e-6 |
| weld, watertight / manifold / orientation report, volume/area/centroid | cube, icosphere (O(h²) volume convergence), torus, cylinder | exact for polyhedra; analytic limits within 1e-3 |
| BVH closest point (Ericson 2005) | brute force on 3000 points | identical to 1e-14 |
| signed distance: pseudonormal sign (Bærentzen & Aanæs 2005); winding numbers for open meshes | box (exact, incl. edges/corners), sphere and torus vs analytic | box 1e-12; curved within the computed Hausdorff/sagitta bound; open-mesh signs correct |
| STEP/IGES/BREP via gmsh + OpenCASCADE (project-local venv, `fluid/scripts/setup-cad.sh`) | STEP cylinder made by OpenCASCADE, imported | watertight; volume error 1.2e-3 vs inscribed-polygon deficit 1.6e-3 |
| VTK writers (grid cell/point data, surfaces) | files read back by the reference VTK library (9.5.2) | every array and the geometry match |

### Flow solver (`immersed` crate)

The solver is a staggered (MAC) grid with conservative second-order central
differences and SSP-RK3 time stepping. Each stage ends with a Leray
projection, and the pressure is incremental, so a steady state satisfies the
discrete steady equations exactly at any time step. The boundary conditions
are velocity inlets (uniform or parabolic), no-slip and moving walls,
pressure outlets, slip and periodic faces.

The immersed boundary is sharp-interface direct forcing:
* nodes inside a body take the rigid-body velocity;
* each fluid node next to the body takes the quadratic through the wall
  velocity and two probes along the CAD normal (Majumdar, Iaccarino & Durbin
  2001).

Cells whose faces are all held by the body take their divergence as a mass
source (Kim, Kim & Choi 2001); in a closed box their net leak is spread over
the free cells. Bodies can translate and rotate. A body can also be
"inside-out", with the exterior solid, for housings and pipes.

| check | reference | result |
|---|---|---|
| pressure solve, all end combinations | discrete Laplacian | inverted to 3.1e-15 |
| channel with inlet/outlet/walls | closed-form *discrete* solution `A(y(1−y) + h²/4)` | to 1e-12 (error vs the parabola is `h² − 3h³`) |
| Taylor–Green vortex | exact | second order (ratios 3.93, 3.98) |
| lid-driven cavity, Re = 100, 64² | Ghia, Ghia & Shin 1982 | within 0.0086 of the lid speed on both centrelines |
| Taylor–Couette between two **immersed** cylinders (one rotating, one inside-out) | exact `Ar + B/r` | L2 orders 2.0 then 2.9, L∞ 1.6 then 3.0; torque error 0.47 % → 0.055 % (R₁/h = 9.6 → 19.2) |
| wall shear by surface sampling (outer Couette wall) | exact `2μB/R²` | 8.6 % → 1.5 % → 0.43 % (linear reconstruction: 13 % → 5.9 % → 3.3 %) |

Two problems found on the way, kept because they matter to anyone building
this:
* **Mass-source cells.** Without them, the projection fights the forcing and
  the incremental pressure drifts. The observed residual grew 1.5× per 55
  time units. Exempting too many cells, every cell with a single held face,
  lets the free face at a stagnation point run away.
* **Round-off in 2D slabs.** A probe normal with a z component of −1e-17
  picked the ghost layer of a single-cell axis. It was latent until an
  optimisation stopped filling that layer.

### DFG benchmark 2D-1 (Re = 20), cylinder read from a file

The reference is Schäfer & Turek (1996), with the FeatFlow values
`C_D = 5.57953523384`, `C_L = 0.010618948146`, `Δp = 0.11752016697`.
Admissible intervals: C_D ∈ [5.57, 5.59], C_L ∈ [0.0104, 0.0110],
Δp ∈ [0.1172, 0.1176].

| D/h | wall reconstruction | C_D (momentum) | C_D (surface) | C_L | Δp |
|---|---|---|---|---|---|
| 10 | linear | +5.21 % | −10.3 % | +71 % | +1.39 % |
| 20 | linear | +1.39 % | −6.6 % | +20.6 % | +0.26 % |
| 40 | linear | +0.48 % | −2.6 % | +1.6 % | +0.24 % |
| 10 | quadratic | +1.99 % | −7.3 % | +51 % | −0.77 % |
| 20 | quadratic | +0.39 % | −3.75 % | +9.3 % | −0.60 % |
| 40 | quadratic | **+0.16 %** | −1.16 % | **−0.3 %** | **−0.23 %** |

At D/h = 40 with the quadratic reconstruction, all three values lie inside
the DFG admissible intervals: C_D 5.5885, C_L 0.01059, Δp 0.11724. The run
took 586 s on 880×164 cells with the CPU shared. The momentum-balance drag
converges at about order 2.3 (+1.99 → +0.39 %) and then 1.3, as it
approaches the level of the reference intervals.

Reading the table:
* **The momentum-balance force is the primary estimator.** Its error is
  exact for the discrete equations and converges at about second order.
* **The surface integral gives the pressure and shear distribution** on the
  CAD triangles. Its total converges more slowly.
* **Δp uses linear extrapolation of the wall pressure**, which is first
  order; it is the weakest number.
* **C_L ≈ 0.0106 is a tiny difference of large forces.** It is the hardest
  quantity for any immersed method.

### DFG benchmark 3D-1Z (Re = 20), cylinder from STEP

Reference values from Braack & Richter (2006): `C_D = 6.185333`,
`C_L = 0.009401`, `Δp = 0.171007`. The channel is 2.4984 m long instead of
2.5 m, so that 0.41 m is a whole number of cells.

| D/h | cells | C_D (momentum) | C_D (surface) | C_L | Δp | wall time |
|---|---|---|---|---|---|---|
| 15.6 | 390×64×64 = 1.6M | 6.2955 (+1.78 %) | 5.961 (−3.6 %) | −0.043 (ref +0.0094) | 0.17089 (−0.07 %) | 15 min, 1894 steps |

The run used 16 CPU threads and the cylinder came from the STEP file, so
the whole case file path was exercised.
* **Drag and Δp** are where the 2D convergence study puts them at this
  resolution: the 2D drag error was +1.99 % at D/h = 10 and +0.39 % at 20.
* **The lift is wrong in sign and size.** The physical lift is 0.15 % of the
  drag. At D/h = 15.6 the cylinder centre is not on grid lines, and the
  up/down asymmetry the immersed boundary introduces is larger than that.
* **No convergence study in 3D.** The next grid (D/h = 23.4, 5.4M cells,
  about 1.2 GB) was not run: the machine had about 1 GB free alongside the
  desktop applications. `fluid case` runs it unchanged on a larger machine.

**Setup cost.** Exact closest-point queries are made only near the body, and
the case uses a sensible tessellation (68k triangles). The case is ready to
step in 4 s with a 354 MB peak. The first version took 74 s and 832 MB,
because 512 elements per 2π gave 1.08M triangles.

## Using it

```sh
fluid/scripts/setup-cad.sh                 # once, for STEP/IGES import
cd fluid && cargo build --release -p fluid-cli
./target/release/fluid case cases/schafer-turek-2d/case.json
```

The case file (see `immersed::case` docs and `fluid/cases/*/case.json`)
lists:
* the fluid properties;
* the box and the cell size;
* one condition per face;
* the bodies (file, scale/translate, rigid motion, `inside_out`,
  tessellation);
* the run (steady tolerance or `t_end`);
* the outputs.

Results go to `out/`:
* `forces.csv`: the force history;
* `summary.json`: forces (momentum and surface), coefficients, torque, probe
  wall pressures and convergence;
* `fields.vtk`: velocity, pressure, vorticity and distance;
* `body*.vtk`: pressure, Cp and wall shear on the CAD surface.

Setup errors, such as a non-integer extent or unpaired periodic faces, are
reported as messages, not panics.

## Honest limits (what an engineering buyer will hit)

* **Laminar only.** There is no turbulence model (RANS, LES) and no wall
  function. Typical engineering Reynolds numbers (1e5–1e7) are out of reach
  except as under-resolved DNS. This is the largest gap against Ansys,
  COMSOL and Siemens.
* **Uniform grid.** There is no local refinement near bodies. Large 3D
  cases cost cells everywhere. The compressible crate has AMR, but this
  solver does not.
* **CPU only.** The spectral and compressible solvers have CUDA backends;
  this one does not yet.
* **Explicit time stepping.** Diffusion limits the step at low Reynolds
  number, so steady states take thousands of steps. Implicit diffusion or a
  steady solver would help.
* **Geometry is static.** Rigid motion enters through the wall velocity, so
  rotating shafts and moving walls work. Bodies that actually move, and
  fluid–structure interaction, are not supported.
* **Physics.** There is no heat transfer, no species, and a single phase.
* **Wall pressure is extrapolated linearly from samples at √dim·h and
  2√dim·h.** That works at stationary walls: in Couette the outer-wall
  pressure is within 0.1 %, the 2D Δp within 0.6 %, and the 3D Δp within
  0.07 %. It fails at rotating walls, where the centrifugal pressure
  curvature (`p'' = −4.3` at the inner Couette wall) makes any polynomial
  extrapolation over 2–4 cells miss by 5–20 %. Those walls need the
  wall-normal momentum equation instead.
* **The surface-integrated drag is low by 2.6–3.8 %** at D/h = 20–40. The
  momentum-balance force is the one to quote.
* **CAD import uses an external tool** (gmsh/OpenCASCADE, 39 MB wheel,
  project-local). The repository has no B-rep kernel. STL/OBJ need nothing
  extra.
* **Tessellation is isotropic.** Over-fine curvature settings explode the
  triangle count: 512 elements per 2π gives 1.08M triangles on the 3D
  cylinder, 128 gives 68k at a sagitta still 400× below the cell size.
