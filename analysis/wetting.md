# Wetting of textured surfaces

The wetting instance of the ladder in
[emergence_ladder.md](emergence_ladder.md):

| rung | model | validity indicator | unfold / closure |
|---|---|---|---|
| coarse | an **effective contact angle** `θ*` at the flow-solver wall | `a/λ`, the drop base radius over the texture period | — |
| fold | Wenzel `cos θ* = r cos θ_Y`; Cassie–Baxter `cos θ* = φ_s(cos θ_Y+1) − 1` | the Cassie→Wenzel energy criterion, the impalement pressure, `Δcos θ` hysteresis | resolve the texture |
| fine | a 2D drop on an explicit square-wave substrate, free-energy minimised | local meniscus curvature and contact-line pinning | — |

The coarse model is a **closure for the wall boundary condition**: the
compressible multiphase solver represents surface tension with a
continuum (Brackbill) force (`fluid/crates/compressible/src/capillary.rs`,
§9 of [compressible.md](compressible.md)); at a solid wall the missing
ingredient is the equilibrium angle that the diffuse interface makes
with the wall. A phase-field wall condition imposes
`κ ∂_n α = ±γ cosθ*` (Jacqmin, *J. Comput. Phys.* **155**, 96 (1999);
**161**, 463 (2000); Ding & Spelt, *J. Comput. Phys.* **226**, 2078
(2007)); the texture supplies `θ*` through the fold, and `a/λ` says
whether the fold may be trusted. This crate computes the fine model,
folds it, and checks the fold.

The crate has no external dependencies; the water constants live in
`constants.rs`, each with its source.

## 1. Exact capillary statics

Axisymmetric sessile/pendant drops are integrated from the apex with the
arc-length Young–Laplace equations in the Bashforth–Adams form
(`young_laplace.rs`),

```text
dr/ds = cos φ,  dz/ds = sin φ,  dφ/ds = 2/b ± Bo z − sin φ / r,
Bo = ρ g b²/γ,
```

with `+` for sessile and `−` for pendant. Lengths are scaled by the apex
radius `b`, so the shape is parameter-free. Volume and area accumulate as
`dV/ds = π r² sin φ`, `dA/ds = 2π r`.

The planar results (`capillary.rs`) are closed-form:

* Jurin rise in a tube of radius `r`: `h = 2 γ cos θ /(ρ g r)` (Jurin
  1718; de Gennes, Brochard-Wyart & Quéré, *Capillarity and Wetting
  Phenomena*, 2004, ch. 2).
* Wall meniscus: `γ dφ/ds = −ρ g z` has the first integral
  `z = 2 l_c sin(φ/2)`, the wall tangent is `φ_w = π/2 − θ`, so the
  height is `z_w = l_c √(2(1−sin θ))` and
  `x(z) = l_c[F(u_w) − F(u)]`, `F(u) = 2√(1−u²) − artanh√(1−u²)`.
* Large 2D puddle: `h = 2 l_c sin(θ/2)`.

`l_c = √(γ/ρg) = 2.7271 mm` for water at 20 °C.

### Benchmarks (water at 20 °C)

| check | exact / reference | computed | error |
|---|---|---|---|
| zero-gravity sessile vs spherical cap `h = R(1−cosθ)`, `a = R sinθ` | exact, θ = 30…120° | Bashforth–Adams RK4, `ds=1e-4` | `max‖Δh,Δa‖ = 4.3e-6` (apex-radius units) |
| zero-gravity cap volume `πR³/3(2−3cosθ+cos³θ)` | exact | integrated `dV/ds` | `< 3e-3` |
| pendant profile, `Bo = 0.21` | OpenDrop reference (Berry, Neeson, Dagastine, Chan & Tabor, *J. Colloid Interface Sci.* **454**, 226 (2015)) | integrated | `max|Δr| = 5.6e-6`, `max|Δz| = 3.8e-5` |
| pendant volume 5.53648, area 15.9890 | OpenDrop | integrated | `ΔV = 1.5e-4`, `ΔA = 2.7e-4` |
| Jurin, `r = 0.5 mm`, θ = 0 | 29.7474 mm | closed form | exact |
| planar puddle `2 l_c sin(θ/2)` | exact | RK4 from the first integral | `0` (integration ends on θ) |
| wall meniscus `x(z)` | closed form | RK4 `dφ/ds = −(2/l_c) sin(φ/2)` | `< 2e-3 x` |

The OpenDrop reference numbers are at `Bo = U²... = 0.21` in its
`test_younglaplace.cpp`, quoted to the digits used there; the agreement
is at the level of the reference's own tabulation. The zero-gravity
spherical-cap limit is the sharpest test: with `Bo = 0` the Bashforth–
Adams system reduces to the exact cap.

## 2. Texture homogenisation (the fold)

Geometry is exact (`texture.rs`); the hysteresis formulas are empirical
and labelled.

### Geometry

Square array of square posts, pitch `λ`, post width `a`, height `h`:

```text
φ_s = a²/λ²,    r = 1 + 4 a h/λ².
```

Hexagonal array of regular hexagonal posts, side `s`, lattice constant
`λ`: unit-cell area `(√3/2)λ²`, hexagon area `(3√3/2)s²`, perimeter `6s`,

```text
φ_s = 3 s²/λ²,   r = 1 + 6 s h /((√3/2)λ²) = 1 + 4√3 s h/λ².
```

For a 2D (cylindrical) square wave the same quantities are `φ_s = a/λ`
and `r = 1 + 2h/λ` (the factors differ because the area/line measure
does; the crate keeps the 2D versions for the fold check and the
3D post versions for the flow-solver boundary condition).

### Wenzel and Cassie–Baxter

Wenzel (*Ind. Eng. Chem.* **28**, 988 (1936)): the liquid penetrates the
texture, so the wetted area is `r ×` projected,

```text
cos θ*_W = r cos θ_Y.
```

Cassie & Baxter (*Trans. Faraday Soc.* **40**, 546 (1944)): the liquid
rests on the post tops with vapour in the gaps,

```text
cos θ*_CB = φ_s (cos θ_Y + 1) − 1.
```

The Cassie state has the lower surface energy when
`cos θ_CB > cos θ_W`, i.e.

```text
cos θ_Y < (φ_s − 1)/(r − φ_s).
```

Both are exact area/energy balances; no fitting is involved.

### Impalement of the Cassie state

A bridging meniscus over a gap of width `w` can hold a Laplace pressure
at most `p_imp = −2γ cos θ_Y/w` (Bico, Thiele & Quéré, *Europhys. Lett.*
**61**, 348 (2003); Lafuma & Quéré, *Nature Mater.* **2**, 457 (2003)).
For a drop of radius `R` the base pressure is `2γ/R`, so a Cassie drop is
lost below

```text
R_c = w/|cos θ_Y|.
```

The sag criterion: at the impalement pressure the pinned meniscus has
radius `R_c = γ/|p_imp| = w/(2|cos θ_Y|)` and sagitta
`s = R_c − √(R_c² − w²/4)`; impalement by sag needs `s ≥ h`. For the
worked example `λ=10 µm`, `a=3 µm`, `h=8 µm`, `θ_Y=110°`:
`p_imp = 7.11×10³ Pa`, `R_c = 20 µm`, `s = 0.61 µm ≪ h`, so sag is not
the limiting mechanism — the Laplace-pressure criterion `R_c` is.

### Hysteresis (empirical)

Joanny & de Gennes (*J. Chem. Phys.* **81**, 552 (1984)) show that a
dilute array of strong pinning defects gives a hysteresis proportional to
the defect density; Reyssat & Quéré (*J. Phys. Chem. B* **113**, 7996
(2009)) confirm the scaling on fakir surfaces and find it controlled by
`φ_s`. The crate provides

```text
Δcos θ = −C φ_s,
```

with `C` a **supplied pinning strength**; the scaling is empirical and
the coefficient is **not predicted** by the crate. This is stated in the
doc comment of `dilute_defect_hysteresis_delta_cos`.

## 3. Validity of the fold

The fine model (`lattice.rs`) resolves an explicit square-wave substrate
and minimises the true capillary energy of a 2D drop over meniscus
shapes. The free surface is parametrised as a circular arc (the exact
gravity-free 2D shape) with the **exact** wetted arc length of the square
wave under the footprint (vertical walls and groove floors included), and
the **exact** spherical-cap area relation `A = R²(ψ − sinψ cosψ)`, valid
for obtuse contact angles. The only degree of freedom left is the contact
line position `x_c`; the energy

```text
E(x_c) = γ L_arc(x_c) + (γ_SL − γ_SV) L_wet(x_c)
```

is minimised directly over `x_c`. For the Cassie branch `L_wet` is the
projected post-top length and the liquid–vapour menisci over the grooves
add `γ(1−φ_s)·2x_c`. The apparent angle is read from the fitted arc at
the mean plane. The validity indicator is `a/λ = x_c/λ`.

This is the sharp-interface variant the task suggests; it resolves the
texture geometry and the contact-line pinning, but **not** the local
meniscus curvature inside a groove and **not** the sub-λ penetration of
the free surface. The convergence below therefore comes from the
contact-line degree of freedom and the exact wetted length, not from
resolving the fine meniscus. That is the honest limitation.

### Flat-wall benchmark

On a flat substrate the exact minimum is the circular segment with
Young's angle for every size. Measured apparent angle minus `θ_Y`,
`A = 1000`, `30° ≤ θ_Y ≤ 140°`:

```text
max |θ* − θ_Y| = 0.031°.
```

### Fold convergence

Square wave `λ=10`, `a=4`, `h=3` (so `r=1.6`, `φ_s=0.4`), `θ_Y=110°`.
Predicted `θ*_W = 123.18°`, `θ*_CB = 137.46°`.

| area | `x_c/λ` | Wenzel measured | error | Cassie measured | error |
|---|---|---|---|---|---|
| 400 | 1.18 | 111.14° | −12.04° | 137.16° | −0.30° |
| 1600 | 2.20 | 117.86° | −5.32° | 131.32° | −6.14° |
| 6400 | 4.20 | 121.68° | −1.50° | 142.89° | +5.43° |
| 25600 | 8.20 | 123.55° | +0.37° | 134.25° | −3.21° |
| 102400 | 16.20 | 124.46° | +1.28° | 137.14° | −0.32° |

The error is **not monotonic**: the contact line pins at the post
edges/corners, and the energy landscape has many local minima
(`E(x_c)` has cusps wherever the contact crosses a wall). A given drop
can land on a pinned branch whose apparent angle is a degree or two off
the fold, while the *same* drop on another pinning phase is closer. That
scatter is the physical metastability of a Cassie/Wenzel state on a
sharp texture, not numerical noise; it is exactly why a hysteresis width
exists. The large-drop limit is approached by the envelope, and the
largest Cassie drop is within `0.32°` and the Wenzel drop within a few
degrees. The tests check the physically meaningful statement — the error
at `x_c/λ ≳ 4–16` is smaller than at `x_c/λ ≈ 1–2` and below a stated
bound — rather than monotonicity.

Two complementary branches appear in the same energy at the same size
(Wenzel ≈ 124°, Cassie ≈ 137°), separated by an energy barrier; the
lower one is the stable state and the other is metastable, in agreement
with the analytic criterion `cos θ_Y < (φ_s−1)/(r−φ_s)` (here the
threshold is `θ_Y = 119.1°`, so at `θ_Y = 110°` the Wenzel state is
favoured).

## 4. How the fold plugs into the flow solver

For the compressible diffuse-interface solver (`capillary.rs`, §9 of
[compressible.md](compressible.md)) the wall condition on the reference
volume fraction is `κ ∂_n α = γ cos θ*` (Jacqmin 2000; Ding & Spelt
2007). The wetting crate supplies:

* `texture::{SquarePosts, HexPosts}::{wenzel_cos, cassie_cos}` and the
  transition/impalement tests for the local textured patch;
* the indicator `size_ratio(base_radius, period) = a/λ`: use the fold
  when `a/λ ≫ 1`; otherwise resolve the texture or use the hysteresis
  band;
* `lattice::minimise` to produce the fine apparent angle when the
  indicator is small.

## 5. What is exact, what is empirical, what is unresolved

* **Exact**: all of §1; the Wenzel/Cassie area balances and the
  geometry of §2; the Cassie→Wenzel energy criterion; the impalement
  pressure and sag geometry.
* **Empirical**: the hysteresis scaling `Δcos θ = −C φ_s` (the
  coefficient is not predicted).
* **Unresolved / negative**: the fold validity model resolves the
  texture and the pinning but approximates the free surface as an arc
  and its wetted length as the homogeneous square-wave value; it
  therefore cannot see the sub-λ penetration profile or a finite
  impalement energy barrier. A full Cahn–Hilliard or lattice-gas
  resolution of the contact line was attempted during development but
  was not used: the symmetric double-well diffuses a fixed-domain drop
  into a uniform film, and the eight-neighbour lattice-gas carries large
  anisotropic surface-tension artefacts (digital-circle perimeter about
  2.5× the Euclidean value) that dominate the contact angle. Those
  attempts are not shipped. The shipped minimiser is deterministic, has
  no lattice artefacts, and its flat-wall benchmark is exact to 0.03°;
  its fold convergence is reported above with the pinning scatter that
  the model genuinely produces.

## Sources

* water `γ = 72.8 mN/m` at 20 °C — Vargaftik, Volkov & Voljak, *J. Phys.
  Chem. Ref. Data* **12**, 817 (1983); CRC Handbook.
* water `ρ = 998.21 kg/m³` at 20 °C — IAPWS-95 / CRC Handbook.
* `g = 9.80665 m/s²` — CGPM 1901.
* Young, *Phil. Trans. R. Soc.* **95**, 65 (1805).
* Wenzel, *Ind. Eng. Chem.* **28**, 988 (1936).
* Cassie & Baxter, *Trans. Faraday Soc.* **40**, 546 (1944).
* Jurin, *Phil. Trans. R. Soc.* **30**, 739 (1718).
* Bashforth & Adams, *An Attempt to Test the Theories of Capillary
  Action*, Cambridge (1883).
* Bico, Thiele & Quéré, *Europhys. Lett.* **61**, 348 (2003).
* Lafuma & Quéré, *Nature Mater.* **2**, 457 (2003).
* Joanny & de Gennes, *J. Chem. Phys.* **81**, 552 (1984).
* Reyssat & Quéré, *J. Phys. Chem. B* **113**, 7996 (2009).
* Berry, Neeson, Dagastine, Chan & Tabor, *J. Colloid Interface Sci.*
  **454**, 226 (2015) (OpenDrop).
* de Gennes, Brochard-Wyart & Quéré, *Capillarity and Wetting
  Phenomena*, Springer (2004).
* Jacqmin, *J. Comput. Phys.* **155**, 96 (1999); **161**, 463 (2000).
* Ding & Spelt, *J. Comput. Phys.* **226**, 2078 (2007).
