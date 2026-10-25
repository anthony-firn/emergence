# Rendering: the same geometry the solvers use, and the same physics

A renderer in this repository has to answer to the rest of it. It draws
`cad::Solid` directly, with no export step and no second representation to
keep in sync; its light transport is `radiation::optics`', its colour is
`waveoptics::color`'s, its diffraction is `waveoptics::diffraction`'s and
its layered surfaces are `photonics::multilayer`'s. What the `viz` and
`viz-cuda` crates add is the geometry path, the choice of optical model,
and the device backend.

```
cad::Solid ----sphere trace (trace.rs)----+
simulation field (field.rs) --------------+--> Render (image.rs)
camera (camera.rs) -----------------------+     beauty + depth, normal,
     |                                          object id, position,
scene: materials, lights (scene.rs)              march steps, field
     |
integrator (shade.rs): direct lighting, spectral path tracing
     |
wave.rs: the optical ladder — rays, Fraunhofer, Fresnel, Rayleigh–Sommerfeld
     |
anim.rs: parameter curves driving cad::model, NURBS camera paths
     |
viz-cuda: the field compiled to CUDA, sphere-traced on the device
```

## Why sphere tracing is sound here

The CAD kernel guarantees `|sdf(p)| ≤ dist(p, boundary)` at every node of a
model (`analysis/cad.md`). A ray step of `|sdf(p)|` therefore cannot reach
the boundary, let alone cross it, so the march converges on the **first**
crossing and cannot tunnel through thin geometry — the standard failure of
a hand-written distance estimator. The test that matters: a 10 µm shell at
1 m range is found by all 200 rays aimed at it, to 1e-9.

What the guarantee does not cover is cost. A ray that grazes a surface
tangentially takes unboundedly many steps, because the field falls without
the ray crossing. So the march has three outcomes, not two, and
`Trace::Exhausted` is distinct from `Trace::Miss`: a renderer that
conflates them puts holes in grazing geometry and offers no way to tell.
The exhausted case converges linearly in the step count — 64 steps reach
4e-4 of the surface, 4096 get ten times closer — which is the evidence
that it is a cost problem and not a correctness one.

| check | reference | result |
|---|---|---|
| ray–sphere and ray–box hits, 2000+ rays each | the analytic root | zero analytic hits missed, none invented |
| hit position | `|sdf|` at the hit | within the tolerance |
| hit distance | `tolerance / cos(incidence)` | the right bound for a grazing hit, not a loose constant |
| a 10 µm shell at 1 m | its outer surface | 1e-9, no tunnelling |
| normals from the field gradient | the analytic sphere normal | 1e-8; a box's face normals exactly |
| a four-operation boolean tree | no step inside the solid | checked at 40 points along every ray |
| camera projection | its own rays | inverse to 1e-14 |
| thin-lens rays | meeting at the focus distance | 1e-12 |

The projection being the exact inverse of the rays is not tidiness: it is
what lets a measurement be overlaid on a render, or a probe point be
located in the image.

## Rays are rung zero

A ray tracer is a model of light with a domain of validity, not a
description of it. `wave.rs` makes that explicit by putting it on an
`emergence::Ladder` (`analysis/emergence_ladder.md`):

```
rung 0  geometrical optics (the shadow)             cost 1
rung 1  Fraunhofer                                  cost ~1
rung 2  Fresnel (exact quadratic phase)             cost ~1
rung 3  Rayleigh–Sommerfeld (reference)             cost ~N²
```

Ask for a tolerance and the cheapest rung whose claim meets it is
evaluated. Rungs 1 and 2 carry the **rigorous** bounds `waveoptics` proves
for them, so where they apply the ladder answers with a proof.

The interesting part is where they do not apply. Fraunhofer and Fresnel are
*far-field* models, and their own bounds say so: at a Fresnel number of
4545 they claim errors of 8e8 and 7e15, so the ladder refuses them and a
tight request falls through to Rayleigh–Sommerfeld. At F = 0.0091 the
picture inverts — rays sit at their ceiling and Fresnel serves 1e-3 with a
rigorous 1.9e-7. The machinery earns its keep by rejecting rungs, not only
by selecting them.

| check | reference | result |
|---|---|---|
| the exact knife edge at the geometric edge | `1/4` of the open irradiance (Born & Wolf §8.7.1) | 0.250000000 |
| its first maximum | 1.370 at `v = 1.217` (table 8.1) | 1.370443 at 1.2172 |
| a slit's first minimum | `sin θ = λ/a` | 0.00e0 |
| two-slit fringes | `m λ/d`, with midpoints vanishing | exact |
| Airy's first zero | `1.2197 λ/D` | 1e-9 |
| a thin film | the two-beam Airy closed form | 1e-12 at four thicknesses |
| a non-absorbing film | `R + T = 1` | 1e-12 |
| Brewster's angle | `R_p = 0` at `atan(n₂/n₁)` | below 1e-6 at 56.31° |
| photon statistics | Poisson: variance = mean | four exposures |

### Two findings about shadows

**Diffraction misplaces no net irradiance.** The signed error of a sharp
shadow integrates to zero, halving as the range doubles (2.6e-3 to 3.2e-4
over ±40 to ±320 Fresnel units). Light moves across the edge — 37 %
overshoot on the lit side, a tail into the shadow — but none is created.

**A sharp shadow is far better than it looks, once averaged over a pixel.**
Pointwise at the edge geometrical optics is wrong by exactly 3/4 (it says
1, the exact solution says 1/4), and that never shrinks however small the
blur. But over a pixel the step error cancels against the symmetric part of
the fringe pattern, leaving only the fringe asymmetry:

| penumbra width (pixels) | worst pixel error |
|---|---|
| 0.0023 | 5.7e-5 |
| 0.023 | 8.9e-4 |
| 0.235 | 0.027 |
| 0.742 | 0.119 |
| 7.42 | 0.696 |

So the indicator that matters to a renderer is the **penumbra width in
pixels**, `√(λz)/pixel`, not the Fresnel number, and the claim is
`min(0.25·blur, 3/4)`. The ceiling is derived — it is the pointwise edge
error. The slope is *fitted*, and the reason is worth recording: the
obvious derivation ("a pixel of width `dv` containing the edge misplaces a
fixed amount `K`, so its error is `K/dv`") fails precisely because `K = 0`
by the conservation above. The error therefore falls faster than `1/dv` for
a centred pixel, and the worst case is an offset one.

## Caustics, and what they measure

A caustic is where a family of rays folds: the map from the source to the
image plane becomes singular and the ray density goes to infinity. It is
the hardest thing for a renderer to get right and the easiest thing to
check, because unlike a diffuse render its *position* has a closed form.

`caustic.rs` traces light forward — from the emitter, through the specular
chain, onto a receiver (Arvo 1986; Jensen 2001) — and splats each photon
with a bilinear tent, which conserves flux exactly rather than nearly. The
bookkeeping is part of the interface: every photon's flux ends up
`deposited`, `spilled` (landed outside the receiver) or `lost` (absorbed or
escaped), and `PhotonStats::accounted()` is their sum over what was
emitted. Across the tests it is 1 to within 1e-12, which is the rounding of
summing a few hundred thousand terms and not a tolerance anybody chose.

### Selection, not approximation

`PhotonOptions::select_events` keeps only photons with exactly `k`
specular events. This is how one caustic is separated from another and it
biases nothing: for a sphere in parallel light, two events is the bright
forward beam, three is the primary rainbow, four is the secondary. Without
it the forward beam buries the bows — which is also why a rainbow needs
rain and not a window pane. The rejected flux is counted as lost, so the
books still balance.

### What the references are

| reference | closed form | measured | error |
|---|---|---|---|
| primary rainbow deviation, n = 1.333 | 137.9219° | 137.9739° | **0.052°** |
| secondary rainbow deviation | 129.1092° (folded) | 129.1073° | **0.0019°** |
| spherical catacaustic vs adjacent-ray intersections | — | — | **< 7e-7** at a 1e-6 step |
| paraboloid marginal focus at any height | R/2 | R/2 | **< 1e-12** |
| flux accounted for | 1 | 1.000000000001 | **1e-12** |
| noise scaling, `sd·√N` over a 16× range of N | constant | 0.272–0.336 | **1.24×** spread |

The rainbow numbers are the sharpest statement here. Descartes'
construction gives the least-deviation angle in closed form; the photons
have to find it by piling up there, with no knowledge of the answer, after
refracting in, reflecting internally and refracting out. 0.052° on the
primary from 400 000 photons (17 442 of which are three-event) is a
measurement of the parabolic fit's resolution, not of the optics. Water's
dispersion (Daimon & Masumura 2007, via `waveoptics`) puts the red bow at
42.295° and the violet at 40.629°, a width of 1.67° against a published
1.8°.

### Caustics as an instrument: three things they found

This is the part worth keeping. Each of these was a real defect that the
caustic test exposed and nothing else had.

**`cad::kits::optic::Mirror` put its substrate on the wrong side.** The
back face was placed at `sag(aperture) + thickness`, which fills the dish
with glass; a traced paraboloid then reflects off a flat back and focused
0.24% of its flux into a 40 mm aperture. Nothing in the CAD tests could see
this — the solid was watertight, the right volume and the right sag. Only
shining light at it showed the material was on the lit side. The fix is a
back at `min(0, sag(aperture)) − thickness`, and `marginal_focus` now takes
the incoming direction as `−z` and refuses a negative path length, because
getting it backwards returns the same focal length as a *virtual* focus
behind the mirror.

**The Airy asymptotic was missing its `sin` series.** `airy()` evaluates a
Taylor series near the origin and an asymptotic expansion beyond. The
oscillatory branch was only the leading `cos(ζ − π/4)` term, and near a
zero of that cosine the omitted `u₁/ζ` term is not a small correction: at
`x = −8` it is 3% of the answer. The test that caught it compares the two
branches *at the same argument* — comparing across the crossover instead
measures `Ai'`, which is O(1) there and hides a 3% error completely. With
both series to five orders (DLMF 9.7.9) the branches agree to 2e-10 at
`x = −7` and the worst relative error anywhere is about 1e-7 just above
`x = 6`, where `Ai` is already 1e-6.

**A reference constant of mine was wrong, again.** The asserted `Ai(5)` was
`1.0832779e-4`; the true value is `1.0834442808e-4`. The series disagreed,
I assumed the series was broken, and arbitrating it in 60-digit arithmetic
showed the series had been right. That is the third time in this work that
a number I supplied as the reference was the thing in error, which is an
argument for deriving references rather than recalling them — the
catacaustic below is now derived in its doc comment for that reason.

### Where the model's limits show, and how they are stated

Two of these tests fail if written as a single threshold, and in both cases
the threshold would have been wrong rather than the code.

**A spherical mirror's caustic cusp is a narrow-beam limit.** The cusp of
the catacaustic is at the paraxial focus `R/2`, but the brightest *plane*
of a wide beam is not: rays out to `h = a` cross the axis anywhere from
`R/2` down to `R − R/(2 cos asin(a/R))`, about `a²/4R` short of it. So the
test narrows the beam and watches the peak walk onto the cusp:

| beam half-width | peak-flux plane | distance from the cusp | `a²/4R` |
|---|---|---|---|
| 0.40 R | 0.4523 R | 0.0477 R | 0.0400 R |
| 0.20 R | 0.4897 R | 0.0103 R | 0.0100 R |
| 0.10 R | 0.4990 R | 0.0010 R | 0.0025 R |

The closed form is checked separately and more strictly, against what an
envelope actually is: the limit of the intersection of neighbouring
reflected rays, which it matches to 7e-7 at a 1e-6 parameter step.

**A traced paraboloid's spot is its mesh, not its optics.**
`Mirror::meridian` revolves a polygon of `points` segments, so every facet
is a flat mirror whose normal is wrong by about `dz/dr` across half a facet
and whose focus is displaced proportionally to the facet width. The exact
surface has no aberration at all — `marginal_focus` is `R/2` to 1e-12 at
every height — so the honest test is the scaling. Shrink the aperture in
step with `1/points` and the caught fraction holds still (0.9973, 0.9972,
0.9973 at 40, 80 and 160 facets); hold the aperture at 5 mm and refining
the mesh tightens the focus onto it (0.17, 0.69, 1.00 at 40, 160 and 640).
That is a statement that the residual is discretisation, which is stronger
than any single pass mark.

A third case is a near-miss worth recording. The convergence test first
used eight runs per photon count and reported the noise falling 7.9× for
16× the photons where `1/√N` predicts 4×. That was not a violation: a
sample standard deviation from `n` runs carries about `1/√(2(n−1))`, so
eight runs give 27% on each and 38% on a ratio. With 24 runs and
well-separated seeds, `sd·√N` is constant to 1.24× over a 16× range, and
the test now asserts that constancy rather than a single ratio.

### The wave caustic

Geometrical optics says the irradiance at a fold goes as `1/√s` on the lit
side and zero on the dark — infinite at the caustic and discontinuous
across it, which is the signature of a model being asked a question outside
its range. The wave answer is finite (Airy 1838; Berry & Upstill 1980):

```text
  I(s) = (2π/√w) Ai²(−s/w),    w = (λ²ρ/2)^(1/3)
```

with `ρ` the caustic's radius of curvature. `fold_geometric(0.0)` returns
infinity rather than a convenient finite number, because that is where the
cheap rung fails and hiding it would defeat the point of having an
indicator. The tests check that the wave peak is finite where the ray model
diverges, that the two agree to 0.3% from forty wave scales out, and that
the peak brightens as `λ^(−1/3)` while the fringe spacing widens as
`λ^(2/3)` — measured 1.2599 and 1.5874 against `2^(1/3)` and `2^(2/3)`.
Whether this matters in an image is the same indicator as everywhere else
in this crate: the wave scale **in pixels**. At `ρ = 1 mm` it is 0.27 px
and invisible; at `ρ = 10 m` it is 5.7 px and the fringes are the picture.

## The GPU

`viz-cuda` compiles the field rather than interpreting it: each node prints
its formula from `cad::solid` into a `__device__` function, and nvrtc builds
it with a fixed marching kernel. The claim is identity, and it is tested as
such.

| check | result |
|---|---|
| device field vs host field, 22 node kinds × 4000 points | **6.0e-16**, about one ulp |
| GPU trace vs CPU trace, per pixel | identical classification; depth 2.2e-15, normals 4.4e-16, **identical step counts** |
| single precision vs the exact field | 4.4e-8 absolute, measured |
| `Solid::Mesh`, `Solid::Custom` | refused by name, including inside a tree |

The performance result on an RTX 3050 laptop card, measured rather than
hoped:

| comparison | ratio |
|---|---|
| GPU f64 vs one CPU core | 8–10× |
| GPU f64 vs all CPU cores | 0.9–1.1×, a tie |
| GPU f32 vs all CPU cores | 1.0–1.8×, inside the noise |
| GPU f32 vs GPU f64 | 1.3× |

**This card does not beat a multi-core CPU at double-precision sphere
tracing.** Two reasons. Its f64 throughput is a thirty-second of its f32,
which puts it in the same class as a modern CPU. And the kernel is
**divergence-bound, not throughput-bound**: neighbouring pixels take wildly
different numbers of steps, a warp waits for its slowest lane, and that is
why single precision buys 1.3× where the arithmetic ratio is 32×. Precision
is the wrong knob; per-warp work balancing is the right one, and is not
done. The test asserts only the robust single-core margin and that f32 is
the faster GPU mode, and prints the all-cores ratio rather than asserting a
win that straddles one.

## Determinism

Every sample's generator is seeded from the pixel and sample indices, so a
render is reproducible to the bit regardless of thread scheduling, and two
renders of the same scene differ only where the scene differs.

## Limits

* **No global illumination acceleration.** The path tracer samples lights
  directly and scatters, with no bidirectional or many-light machinery, so
  a scene lit only by a small emitter converges slowly.
* **No participating media in the implicit path.** `radiation::render` has
  volumes with Henyey–Greenstein scattering; this crate's field
  visualisation samples a field at surfaces and builds isosurfaces, and
  does not yet integrate along the ray.
* **Wave rungs are per-aperture, not per-image.** The ladder answers "what
  is the irradiance behind this rectangular opening"; it does not propagate
  a whole image plane through an arbitrary CAD scene. The indicator says
  when that would matter, which is the useful half.
* **One precision per compiled scene** on the GPU, and no mesh or custom
  nodes there.
* **No denoiser, no adaptive sampling** in this crate (`radiation::render`
  has an adaptive sampler for its own geometry).
* **Caustics are forward-traced and unbiased, but not fast.** There is no
  photon map, no kd-tree, no final gather and no importance sampling of the
  specular chain: a photon that misses the receiver is simply lost. That is
  the right trade for a *reference* — nothing to tune, nothing to bias —
  and the wrong one for interactive work.
* **The wave caustic is the fold alone.** `fold` is the generic
  codimension-one caustic; the cusp and the higher catastrophes are not
  implemented, so a cusp is rendered as rays with no wave correction and no
  indicator saying so.
* The thin-lens camera is an idealisation. Rendering *through* a real
  optical system means `cad::kits::optic` with `waveoptics::raytrace`,
  which has its own dispersion and aberrations.
