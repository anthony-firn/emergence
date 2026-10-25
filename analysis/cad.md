# Parametric CAD: one model, every solver

The solvers in this workspace all need geometry and each needs it in a
different form. The flow solver wants a watertight surface and an exact
signed distance; the drift–diffusion solver wants doped regions on a tensor
mesh; the optics crates want a layer stack; the lattice crate wants a strut
graph. Before `fluid/crates/cad` the only way in was `geometry`, which
*imports* CAD (STL/OBJ, and STEP through gmsh) but cannot make any. This
document covers the generator: what it models, what every number was
checked against, and where it falls short.

```
parameters + expressions (expr.rs)
      |
sketches (sketch.rs) ----extrude/revolve----+
NURBS curves + surfaces (nurbs.rs) ---------+
paths (path.rs) ----------sweep/loft--------+--> Solid (solid.rs)
imported CAD (geometry::Sdf) ---------------+      | booleans, offsets,
                                                   | shells, arrays, polar
generators (kits/): part, lattice, pcb,            | folds, periodic level
silicon, textile, pipe, section, airfoil,          | sets
optic, crystal, sheet  --------------------------->|
                                                   v
              meshing.rs: watertight triangles, anisotropic cells,
                          provably safe block culling
                                                   |
 sim.rs: distance grids, material label grids, layer stacks, strut graphs,
         doping closures, VTK, STL, immersed case bodies
 laplace.rs: 2D electrostatics on a drawn cross-section
```

Eleven generators, 107 tests, each suite under 9 seconds. The strongest
checks are the five that put generated geometry through another crate's
solver and compare against an exact solution — a pipe in the
Navier–Stokes solver against Hagen–Poiseuille, a process flow's device in
the drift–diffusion solver, lenses in the ray tracer, an I-beam in the
frame solver, a lattice against the homogenisation — because those test
the whole path rather than any one formula.

## The kernel invariant, and why there is no B-rep

A `Solid` is a scalar field, not a boundary representation. Every node of
the tree — primitive, boolean, offset, shell, placement, array, periodic
level set — keeps one invariant:

> `sdf` is 1-Lipschitz on R³, negative exactly on the interior, positive
> exactly outside, and zero exactly on the boundary.

Everything else follows from it. A 1-Lipschitz function vanishing on the
boundary satisfies `|sdf(p)| ≤ dist(p, boundary)`, so

* **the sign is exact everywhere.** Membership, meshing and the
  inside/outside test a solver needs are never approximate, however deep
  the tree;
* **the magnitude is a safe lower bound** on the clearance, so the field can
  size a grid, march a ray or offset a surface and can only be
  conservative;
* `Solid::exact` reports when the bound is tight. Primitives, extrusions,
  full revolutions and imported meshes are exact; booleans and arrays keep
  the sign but may fall short near a concave seam.

This is standard function representation (Pasko et al. 1995). The cost is
real and is listed under **Limits**: no NURBS, no exact fillet topology, no
feature references to faces and edges.

The triply periodic level sets are brought inside the invariant rather than
left outside it: each nodal function is divided by an analytic bound on its
gradient (`2|w|` for the gyroid, `|w|` for Schwarz P, `4|w|` for the
diamond, with `w` the per-axis wavenumber). The field is then 1-Lipschitz
and conservative — and *measurably* conservative, which is why a periodic
sheet's wall comes out thicker than its nominal thickness (below).

### A negative result worth keeping

**A morphological fillet cannot be done with offsets.** Eroding and
dilating a distance field cancel exactly: `(f + r) − r = f`. The first
version of `Solid::fillet` was therefore a no-op that type-checked, built
and looked plausible. `Solid::round_convex(r, cell)` re-distances through a
remesh instead, and says in its own documentation that it costs one remesh
and carries that mesh's `O(cell²)` error. The test asserts both halves: that
the offset pair cancels to 1e-15, and that the real operation moves a cube's
corner by `r(√3 − 1)`.

## Exact smooth geometry

Everything above is an implicit field; everything in `path` is a polyline.
Real CAD is neither — an aerofoil, a cam profile, a duct transition or a
lens edge is a rational spline — so `nurbs` adds NURBS curves and surface
patches: Cox–de Boor basis, de Boor evaluation in homogeneous coordinates,
derivatives as a lower-degree spline, knot insertion, chord-length
interpolation, exact rational-quadratic arcs and exact surfaces of
revolution (Piegl & Tiller, *The NURBS Book*, 2nd ed.). Curves feed sweeps
and sketches; surfaces tessellate back into the implicit kernel, so a
smooth surface can still be cut and shelled.

| check | reference | result |
|---|---|---|
| basis is a partition of unity, degrees 1–5 | `sum N_i = 1` | 1e-14, all non-negative |
| **the rational quadratic circle is exact** | `|p − c| = r` | 1e-14 over 2000 points, tilted planes too |
| its circumference | `2 pi r` | 1e-12 |
| clamped endpoints and end derivative | `3(P1 − P0)` for a cubic Bézier | exact |
| analytic derivative | finite differences | 1e-7 polynomial, 1e-6 rational |
| **knot insertion does not move the curve** | the curve itself | 1e-13 over 1000 points |
| interpolation, degrees 2–4 | the data points | 1e-10 |
| **a revolved arc is an exact sphere and torus** | the implicit equations | 1e-13 |
| tessellated area | `4 pi r^2`, `4 pi^2 R r` | second order, 5e-4 at 128×64 |

Two parameterisation traps, both now in the API rather than latent. **The
derivative of a polynomial curve has every weight component zero**, so
carrying the derivative spline as (point, weight) pairs and projecting
through the weight silently returned a zero derivative for every
non-rational curve; it is now carried as raw homogeneous 4-vectors. And
**arc length must be integrated span by span**: the speed kinks at every
knot, so a three-quarter arc was off by 5e-6 where a full circle was exact,
because the full circle's knots happened to align with the uniform
quadrature intervals.

A rational circle's uniform parameter is not uniform angle, so
`to_polygon` is inscribed but is *not* the equal-angle polygon;
`to_polygon_by_length` is, and the tests pin both.

## Meshing

`meshing::surface_nets` is the dual isosurface method (Gibson 1998): one
vertex per sign-changing cell, one quad per sign-changing grid edge, then a
few damped Newton steps projecting each vertex onto the zero set. It is
**closed by construction** — every quad edge is contributed by one of the
four cells around an edge, so there are never boundary edges. It is
*manifold* as well wherever the surface passes through each cell once;
where two sheets meet in one cell the edge through that cell's single
vertex is non-manifold. For two solids in genuine contact that is the right
answer (the boundary really is non-manifold there); for a gap thinner than
a cell it is a resolution artefact that refining removes.

| check | reference | result |
|---|---|---|
| box distance | `geometry::Sdf` on the polyhedron, itself validated analytically | < 1e-14 over 4000 points |
| capped cylinder | its closed form | < 1e-14 |
| capped cone | a 2048-gon loft, brute force | < 1e-6 (the loft's own sagitta is 2.4e-7 r) |
| boolean tree is 1-Lipschitz | finite differences over 4000 random pairs | ratio ≤ 1 + 1e-12 |
| `|sdf| ≤ dist` | brute force to the meshed boundary | excess < 5h², i.e. the mesh's own position error |
| revolution volume | Pappus, exactly `2πR·πa²` | < 2 % by quadrature at h = R/60 |
| polygon area, with a hole | Green's theorem | 1e-15 |
| rectangle second moments | `bh³/12` | 1e-20 |
| sphere volume, area | analytic | 3.9e-4 and 1.8e-4 at h = r/32, order 1.6–2.0 |
| 12-sphere array | analytic | 5e-3 at h = r/10, as the sphere study predicts |
| watertightness | `TriMesh::report` | asserted on every meshed model |

### Making a thin feature in a wide model meshable

A uniform grid was the limit this document used to list as blocking real
work: 35 um of copper on a 1.5 mm board could not be meshed whole. Two
fixes, for different shapes, that compose.

**Anisotropic cells.** Nothing varies quickly in the plane of a board, so
12 um through the thickness and 250 um in the plane resolves the copper
exactly as well as an isotropic 12 um grid, with 3.7M nodes instead of
1.5e9. A 40 mm board with a trace and a pour meshes closed, 264 756
triangles, with the copper layers present in the mesh.

**Block culling**, which is where the kernel invariant pays for itself.
The field is 1-Lipschitz, so if `|f|` at a block's centre exceeds the
distance to the block's farthest corner, `f` cannot change sign anywhere in
the block — one evaluation rules out `block^3` cells. The result is
*identical* to the dense pass, not an approximation: every sign-changing
cell lies in a kept block, every quad's four cells contain the same
sign-changing edge so they are all kept, and cells are visited in the dense
order so the vertex numbering matches. The test asserts bit-identity at
four block sizes, including every vertex position.

| check | result |
|---|---|
| sparse mesh against dense, blocks of 2, 4, 8, 16 | identical vertices and triangle indices |
| a sphere at r/h = 16, 32, 64 | 83 %, 55 %, 32 % of the dense evaluations |
| a plate 1000× wider than thick, anisotropic | 0.5 % volume error on 121 275 nodes |
| the same plate, isotropic at the coarse spacing | surface lost entirely, as it must be |

Culling buys nothing on a board, and the test asserts that rather than
reporting a win it did not get: a slab's surfaces pass through every block
in the plane. Anisotropy wins that case; culling wins the compact one.

## Generators

Every generator is checked against an exact identity or a published
standard. Nothing was tuned to the number it is meant to predict.

### Mechanical parts

| check | reference | result |
|---|---|---|
| thread `d2`, `d1` for M6, M8, M10, M16, M24 | ISO 724 basic dimensions | 5e-4 mm |
| the thread solid's own profile | ISO 68-1: half-pitch metal at `d2`, `P/8` crest flat | 0.500 ± 0.010, 0.125 ± 0.010 |
| lead | the helix closes over one turn | exact, 200 sample points |
| gear tooth thickness at the pitch circle | ISO 21771, `m(π/2 + 2x tan α)` | 1e-6 over four (module, teeth, shift) cases |
| the generated flank really is an involute | its normal line is tangent to the base circle | lever arm = `r_b` to 2e-3 |
| flank chord error | refinement | 7.0e-5 at 6 flank points, 2.1e-7 at 200 |
| contact ratio, m = 2, 20 + 20 teeth | ISO 21771 formula | 1.557 |

The thread is a helical field: the ISO basic profile in the helical
coordinate `u = z − Pθ/2π`, periodic in `u`, divided by the transform's
gradient bound `√(1 + (P/2πr)²)` — 1.002 for M8 — so the kernel invariant
survives. The root is left sharp: ISO's rounded root (which is what defines
`d3`) is **reported, not modelled**, so the solid's minor diameter is the
basic `d1`.

### Lattices

| check | reference | result |
|---|---|---|
| octet truss topology | 24 struts of length `a/√2` per cubic cell | exact |
| relative density | Deshpande, Fleck & Ashby 2001, `6√2 π (r/l)²` | 1e-14 |
| homogenised 6×6 stiffness | `architected::lattices::octet_truss` through independent code | 1e-10 of `C₁₁` |
| level-zero volume fraction, gyroid / Schwarz P / diamond | 0.5, by the odd symmetry of each nodal function | 5e-3 |

Two measured numbers rather than predictions. At `r/l = 0.0707` the
meshed strut lattice's volume fraction is **16.8 % below** the no-overlap
bound: that gap is the material shared where twelve struts meet at a node,
and `StrutGraph::relative_density_upper` is documented as an upper bound for
exactly this reason. And a gyroid sheet's geometric wall is **1.95–2.09×**
its nominal field thickness, because the nodal field is normalised by a
conservative gradient bound; `TpmsSpec::wall_thickness` measures the
achieved value instead of assuming it.

No claim is made about matching published *minimal* surface areas. The
trigonometric level sets are not the minimal surfaces, and the difference
was not measured.

### Printed circuit boards

The kit draws a board and then answers electrical questions about what was
drawn, which is the point: a trace width is chosen for 50 Ω, and the
generator can check its own work. `laplace` solves 2D electrostatics on the
drawn cross-section — permittivity averaged across each face, so the energy
per unit length is a plain face sum — by preconditioned conjugate
gradients, and `quasi_tem` gets `Z₀ = 1/(c√(C C_air))` from two solves.

| check | reference | result |
|---|---|---|
| parallel-plate capacitance | `ε₀ε_r W/d`, exact with mirror walls | 1e-9 |
| two dielectrics in series | the series formula | 1e-9 |
| microstrip `Z₀`, W/h = 0.5, 1.0, 1.8 | Hammerstad 1975 / Hammerstad & Jensen 1980 | −3.70 %, −2.81 %, −0.59 % |
| the same, refined | h/cell = 8, 16, 32 | 67.7, 69.1, 69.8 Ω → −0.8 % of the formula |
| effective permittivity | the same formulas | 3.043 vs 3.057 (0.5 %) |
| 1 oz copper sheet resistance | `ρ_Cu/t` | 0.4822 mΩ/square |
| 1 mm × 100 mm trace | `ρL/wt`, with the temperature coefficient | 48.2 mΩ |
| copper coverage | the analytic capsule area | 1 % |
| boundary tracing (what Gerber emits) | a disc's area and perimeter | 1e-3 and 2 %, converging |

The field solve comes out **below** the formula at every width, with the
same sign, which is what the formula's zero-thickness assumption requires:
the real 35 µm copper raises the capacitance. That is a point in the field
solve's favour, not an error in it.

Gerber RS-274X and Excellon output are written to their specifications and
structurally self-tested, but are **not verified against a CAM tool**. The
FR-4 constants are nominal — real laminate permittivity runs 4.2–4.8 with
resin content and weave — so a controlled-impedance design must supply its
own.

### Silicon process flows

A planar process is masks in `x` and a stack in `y`, which is exactly the
shape `semiconductor::device::Device::new` takes: a material function and a
doping function of `(x, y)`.

| check | reference | result |
|---|---|---|
| **a device built through the CAD path** | the same device from `semiconductor::mosfet::MosfetSpec`, after an equilibrium solve and gate and drain ramps | drain current to **1.2e-15** relative |
| its doping closure | the reference closure over the whole cross-section | 1e-14 of the peak |
| lateral junction, effective channel length | the reference implementation | 1e-15 |
| implant dose, three profile families | the analytic integral | 1e-6 |
| anneal | Fick: `σ² → σ² + 2Dt`, dose fixed | exact; dose to 1e-15 |
| Deal–Grove thickness | `x² + Ax = B(t + τ)` | 1e-12 |
| its short-time limit | `x → (B/A)t` | 1e-6 |
| its long-time limit | `x → √(Bt) − A/2`, residual = `A²/(8√(Bt))` | the next asymptotic term, to 1 % |

A CAD-authored device therefore **inherits** that crate's published
validation rather than claiming its own. No process constants ship with the
crate: the Deal–Grove rate constants and the diffusivities are the
caller's, because tabulated values differ between compilations and depend on
orientation, ambient and doping. Not modelled at all: segregation at the
oxide interface, transient enhanced diffusion, clustering and solubility,
channelling tails (the Gaussian misses these badly in the deep tail),
stress, and anything three-dimensional.

`Process::junction_depth` returns the *shallowest* junction and
`Process::junctions` all of them. A single bisection silently returned an
arbitrary one in a stacked flow, which the three-implant bipolar test
caught.

### Textiles

| check | reference | result |
|---|---|---|
| yarn diameter from count | `d = √(4T/πρφ)` | exact; flattening preserves area to 1e-15 |
| plain, 2/2 twill, 5-end satin | the interlacing definitions and mean floats | exact, including the twill diagonal as a shift identity |
| fabric thickness | Peirce 1937, `h₁ + h₂ = d₁ + d₂` | the generated solid's z extent to 1e-9 |
| Peirce's plain-weave equations | Peirce 1937 | both relations to 1e-9, `h₁+h₂ = D` to 1e-12; a jammed setting reported as having no solution |
| crimp | the exact arc length of the centreline | 0.1237 vs 0.1244 (5.3e-3) |
| swept yarn volume | section area × path length | 1.6e-3, with the inscribed-polygon deficit (a 24-gon holds 0.9886 of the ellipse) accounted for |
| knitted course and wale spacing | Munden 1959, `K_c = 5.0, K_w = 3.8` dry and `5.3, 4.1` wet | exact; `K_s = K_cK_w` = 19.0 and 21.7 |
| the generated stitch's yarn consumption | the loop length it was given | second order: 1.7e-3 at 48 stations, 1.1e-4 at 192 |

Three findings, all kept:

* **The sweep's end caps were wound with the adjacent side triangles**, so
  every swept solid had inconsistent orientation and reported as not
  watertight. The fabrics found it; the threads and the cone reference did
  not, because neither test looked at orientation.
* **The self-intersection indicator was the wrong test.** Comparing the
  section's largest radius with the path's curvature radius forbids a
  flattened yarn bending about its thin axis — exactly the case in cloth.
  It now measures how far the section reaches *toward the centre of
  curvature* at each station, which is the actual condition.
* **A sinusoidal centreline cannot draw a close-set weave.**
  `Weave::consistent_flattening` gives the flattening interval that works,
  and it is empty unless the thread spacing exceeds `d√(π/2)` — the cosine
  crest curves more tightly than the yarn is thick. Peirce's arc path,
  whose radius is `d₁ + d₂`, does not have this problem. Likewise a single
  sinusoid cannot reach a knitted loop's 3.8× length-to-advance ratio
  without bending tighter than the yarn, so the loop is built as two legs
  plus head and sinker semicircles, solved against its true
  three-dimensional length. The resulting loops reach 1.78 (dry) and 1.94
  (wet) course spacings, which is how a real interlooped fabric behaves.

A fabric's yarns are exactly tangent at every crossing, so the union's
boundary is genuinely non-manifold there: the mesh is closed but not
manifold, and nothing in the crate pretends otherwise.

### Pipe and duct runs

The kit that feeds the flow solver most directly. A `Run` is a chain of
straights and elbows; a `Pipe` carries it as the **bore** (the fluid
domain an internal-flow case meshes) or the **wall** (the metal a thermal
or external-flow case wants). Dimensions are ASME B36.10M schedules
10/40/80/160 and ASME B16.9 elbow radii, reproduced with their source and
flagged as reproduced rather than derived.

| check | reference | result |
|---|---|---|
| **flow in the generated pipe**, D/h = 8, 12, 16 | exact Hagen–Poiseuille profile, solved by `immersed` | 2.6 %, 6.7 %, 2.2 % on the velocity profile at 5 D |
| the same run's pressure gradient | `32 mu L u / D^2` | 1.5 %, 4.9 %, 1.2 % |
| NPS 4 Sch 40 | the standard's bore, area, mass | 102.26 mm, 82.1 cm², 16.07 kg/m |
| meshed bore volume | flow area × length | 1 % |
| Darcy–Weisbach at `f = 64/Re` | Hagen–Poiseuille | 1e-15 |
| Colebrook–White, six (Re, roughness) pairs | its own equation | residual 1e-12 |
| smooth pipe at Re = 1e5 | the Moody chart | `f = 0.0180` |
| run length | straights plus `R theta` | exact |
| long-radius elbow centre-to-face | ASME B16.9, `1.5 D` | exact |

**The flow error is not monotone in the resolution**, and the test says so
instead of asserting an order it does not have. For a sharp direct-forcing
immersed boundary the error depends on where the wall falls inside the
near-wall cell: at these three grids it sits 0.33, 0.74 and 0.16 of the way
across the diagonal cell, and the errors follow that ordering rather than
the grid spacing. `real_geometry.md` records the same effect spoiling the
DFG 3D lift.

Two things the kit refuses rather than draws: an elbow whose axis is not
perpendicular to the run, and a run that bends tighter than its own bore.

### Structural sections and aerofoils

`section` draws rectangle, round bar, I-beam, channel, angle, tee, box and
tube, and computes the section properties a beam element needs from the
drawing. `airfoil` draws NACA four-digit sections — the one family here
*defined* by a published polynomial — and lofts blades from them.

| check | reference | result |
|---|---|---|
| area and both second moments | the composite-rectangle formulas | 1e-14 |
| a rectangle, a circle | `bh^3/12`, `pi r^4/4`, `J = 2I` | 1e-15 |
| channel and tee centroids | their first-moment sums | 1e-14 |
| **a CAD section in the beam solver** | `P L^3 / 3EI`, through `architected::frame` | 2.8e-13 |
| every shape extruded | area × length | 3 % |
| NACA maximum-thickness position | 0.2998 chords | 1e-3 |
| NACA leading-edge radius | `1.1019 t^2` chords | 1e-3 |
| NACA area, open and closed trailing edge | `0.685083 t c^2`, `0.680883 t c^2` | 1e-6 |
| surfaces from the camber line | the half-thickness | 1e-14 |
| lofted planform, taper, aspect ratio, MAC | the trapezoid formulas | exact |

Two statements the tests pin rather than tolerate. **The NACA polynomial's
maximum is 1.00027 t, not t**: the series is nominally "t at 30 % chord"
and the published ordinate tables are, but the closed form overshoots by
0.027 % at x = 0.2998. And **a cambered section encloses more area than
`2 integral y_t dx`**, because the thickness is laid off perpendicular to
the camber line — 0.8 % more for NACA 4415 — so `Naca4::area` integrates
the real thing and `area_vertical` keeps the closed form it is.

An open section's torsion constant is marked `exact_j: false`: the
thin-walled sum is low by 10–20 % for a wide flange, and the shapes have no
root fillets, so nothing here reproduces a mill table.

### Optical elements

Lenses, conic mirrors, prisms and windows from the standard
conic-plus-even-asphere sag equation. The same numbers go to two consumers
that must agree: a solid of revolution, and a surface list for
`waveoptics::raytrace`.

| check | reference | result |
|---|---|---|
| **a paraboloid has no spherical aberration** | every parallel ray through `R/2` | 1e-14 of `f`, 50 heights |
| a sphere's marginal focus | falls short as `h^2` | ratio 4.09 for doubled height |
| **lensmaker's equation** | an independent paraxial trace in `waveoptics` | 2e-16 or exact, four lens forms |
| back focal distance | the same two routes | 1e-10 |
| a real ray at 0.1 mm | the paraxial focus | 1e-5 |
| lens volume | the spherical-cap closed form | 1e-12 |
| **prism minimum deviation** | `2 asin(n sin(A/2)) − A`, against 200 001 traced incidences | 1e-6, at the symmetric passage |
| an equilateral N-BK7 prism | 38.647 degrees | 0.01 degrees |
| a right-angle prism of the same glass | passes nothing: apex > `2 asin(1/n)` | total internal reflection reported |

A lens that would close to a knife edge is refused rather than drawn. Two
of my own reference numbers were wrong and the tests caught both: the
spherical-cap volume, and the prism's minimum deviation, which is 38.65
degrees and not 37.2.

### Atomistic structures

What a molecular-dynamics or tight-binding run takes: atom lists and
periodic cells. Simple cubic, body- and face-centred cubic, hexagonal close
packed, diamond, zincblende, rock salt and cubic perovskite; supercells;
graphene; and `(n, m)` nanotubes. **Any solid the implicit kernel can
describe cuts a cluster out of a crystal**, so a nanoparticle, a nanowire
or a hollow shell is `Crystal::carve` with the shape — which is what makes
the atomistic kit part of the same CAD system.

| check | reference | result |
|---|---|---|
| coordination numbers, **counted** from the positions | 6, 8, 12, 12, 4 | exact |
| packing fractions | Ashcroft & Mermin ch. 4 | 1e-14 |
| mass densities of Cu, Si, Fe | handbook values | 8930, 2329, 7878 kg/m³ |
| a 3 nm copper particle | the bulk number density | 1193 atoms against 1197 |
| **nanotube diameter, atom count, chiral angle** | Saito, Dresselhaus & Dresselhaus ch. 3 | 1e-12, exact, exact |
| nearest-neighbour distances | `materials::lattice`, independently | 1e-12 |

That last cross-check found a real bug in `materials::lattice`, now fixed:
`nn_distance` and `bond_list` both used an **absolute** epsilon of 1e-9 to
skip the zero vector. That is sensible for the dimensionless lattice
constants the crate's benchmarks use, but it exceeds every bond in a
lattice measured in metres, so `nn_distance` silently returned a distant
shell and `bond_list` returned no bonds at all — any lattice sum on an
SI-unit lattice would have been silently zero. Both now scale the epsilon
by the primitive vectors; every existing `materials` test is unchanged, and
a new SI-units test pins the first four face-centred shells at 12, 6, 24,
12 in metres as they are at `a = 1`.

### Sheet metal

The one domain here where the manufacturing geometry is a different shape
from the part. The standard bend bookkeeping: allowance `A(R + KT)`,
outside setback `(R+T)tan(A/2)`, deduction `2 OSSB − BA`.

| check | reference | result |
|---|---|---|
| the handbook worked example (1 mm, R = 1 mm, 90°, K = 0.44) | BA 2.262, OSSB 2.000, BD 1.738 mm | 1e-3 mm |
| flat length two ways | apex less deductions; straights plus allowances | 1e-15 |
| the **generated** neutral-axis path | the developed length | second order: 8.8e-4 → 3.4e-6 |
| a U-channel's formed box | its drawing dimensions | 0.2 mm |
| spring stiffness | `G d^4 / 8 D^3 n` | exact, 4699 N/m |
| spring wire length | `n sqrt((pi D)^2 + p^2)` | 1e-15 |

**Volume is conserved only at K = 0.5.** The neutral-axis model stretches
the outside of a bend and compresses the inside, and the two cancel exactly
only at mid-thickness; at K = 0.33, 0.40 and 0.44 the formed solid holds
0.46 %, 0.27 % and 0.16 % more material than the blank. Real metal
conserves volume, so `volume_change` reports the model's discrepancy
instead of hiding it.

**Drawable is not makeable**, and the two are separate calls. A spring of
index 1.6 sweeps perfectly well — the wire clears the coil axis — and
cannot be wound; `solid` reports geometric consistency and
`manufacturability` the index limits.

The K-factor is shop data, not a prediction: published tables disagree by
several per cent and the value depends on the press.

## Handing geometry to the solvers

| check | reference | result |
|---|---|---|
| **a CAD body in the flow solver** | the same cylinder from `geometry::primitives`, `immersed` at Re = 20 | drag to **6.2e-4** relative |
| `Solid` → mesh → `geometry::Sdf` | the analytic distance | second order; 2e-3 of the radius at r/16 |
| **a CAD layer stack in the optics solver** | a quarter-wave coating must not reflect; a half-wave layer must reproduce the bare interface | R < 1e-12; exact to 1e-12 |
| PCB stack thermal resistance | the series sum | 1e-12 |
| label-grid volume fractions | analytic | 1e-12 |

## Using it

```sh
cd fluid && cargo build --release -p fluid-cli
./target/release/fluid cad models/flange.json --cell 0.0006 --stl flange.stl
./target/release/fluid cad models/flange.json --set bolts=12   # a design sweep
```

A model is a parameter table and a feature tree, with any number written as
an expression over the parameters, so one file is a family of parts.
Parameters may refer to one another in any order; a cycle is reported.
`models/` holds four worked examples:

| model | result |
|---|---|
| bolted flange | 113 948 triangles in 0.30 s, watertight |
| spur gear, m = 2, 24 teeth | 952 928 triangles in 14.6 s, watertight |
| gyroid sheet lattice, 3×3×2 cells | 331 916 triangles in 0.38 s, watertight |
| pin-fin heat sink | 135 612 triangles in 0.52 s, watertight |
| NPS 2 Sch 40 elbow run | 251 264 triangles in 4.4 s, watertight; volume within 0.4 % of area × centreline, the 48-gon profile's own deficit |
| formed sheet-metal bracket | 99 904 triangles in 0.64 s, watertight |
| 25 mm biconvex lens | 76 064 triangles in 1.1 s, watertight |
| NACA 2412 wing, 0.6 m span | 240 384 triangles in 3.2 s, watertight |

The kits reachable from a JSON model are `spur_gear`, `thread`, `plate`,
`strut_lattice`, `tpms`, `woven`, `knit`, `spring`, `beam`, `wing`, `lens`,
`mirror`, `prism`, `window`, `sheet_part` and `pipe_run`.

One sharp edge in the expression substitution, found by writing the elbow
model: a parameter named `straight` was substituted into
`{"leg": "straight"}`, because the deny-list of keys that are always names
did not yet include `leg`. It now lists every serde tag in the crate, and
the failure mode is a loud deserialisation error naming the key, not silent
corruption.

Two performance findings came out of running them. A capped cylinder's
bounding box used to be a ball at each end, which inflates it along the
axis by the radius — a sixth of the flange's meshing grid was empty space.
And a 24-tooth gear as a union of 24 tooth polygons costs 24 evaluations
per query; `Region2::PolarRepeat` folds the query's angle into one sector
and checks three copies instead, which is what makes a many-toothed gear
meshable at all.

## Limits

* **No topological B-rep kernel.** `nurbs` has rational curves and single
  untrimmed surface patches, but no trimmed faces, no surface–surface
  intersection, no exact fillet topology and no feature history referring
  to faces and edges — so no STEP export. Exact CAD still comes in through
  `geometry::cad` (STEP/IGES by tessellation) and goes out as triangles.
* **Sharp edges are rounded by the grid** when a solid is meshed: the dual
  isosurface cannot represent a crease sharper than one cell. Analytic
  queries keep the edge exactly.
* **No constraint solver.** Sketches are explicit; `expr` gives
  dimension-driven parametrics, not geometric constraints.
* **No assemblies, mates, tolerances, GD&T or drawings.**
* **Manufacturing output is unverified** against a CAM tool (Gerber,
  Excellon). No tool paths, no post-processors, no DFM checks beyond the
  per-kit refusals listed above.
* **The generators are structure generators, not process simulators.**
  Nothing relaxes a fabric or a nanoparticle, oxidises with real rate
  constants, predicts a K-factor, or chooses a glass; each kit says which
  of its inputs are the caller's measurements.
* **Meshing is grid-based, not adaptive.** Anisotropic cells and block
  culling between them handle a thin feature in a wide model and a compact
  solid in a large box, but there is no octree: a model with features at
  three very different scales *and* no preferred direction still costs the
  finest cell everywhere.
