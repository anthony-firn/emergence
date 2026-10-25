# A sharp solid–fluid interface in 2D

Crate: `fluid/crates/compressible`. Module `compressible::sharp_interface2d`,
validated by `tests/sharp_interface2d.rs`
(`cargo test -p compressible --test sharp_interface2d`, ~15–25 s: ten fast
tests, including an oblique manufactured regression, an oblique-stability
regression and a closed-box energy balance; the grid-convergence sweep and
diagnostic probes are `#[ignore]`d).
This note is the 2D unfold of `analysis/sharp_interface.md`. The motivation
is `analysis/fsi.md`: the crate's only multi-material model is the
diffuse-interface five-equation model, and at a water/steel interface it
gives ~6% absolute `|T|` error at 1024 cells, `|R| > 1`, and an unstable
shear-bearing GPR mixture at that contrast.

## 1. What is modelled

* A **fluid** is a stiffened gas (`Nasg`, `b = q = 0`) in 2D Euler form
  `[ρ, ρu, ρv, E]`, normal stress `τ = p`. Water `ρ = 998.2 kg/m³`,
  `c = 1482 m/s` (CRC Handbook).
* A **solid** is linear isotropic elasticity in **plane strain**, in the
  velocity–stress formulation with compression-positive stresses
  `[ρu, ρv, σxx, σyy, σxy]`:
  ```
    ∂u/∂t   = −(∂σxx/∂x + ∂σxy/∂y)/ρ
    ∂v/∂t   = −(∂σxy/∂x + ∂σyy/∂y)/ρ
    ∂σxx/∂t = −(λ+2μ)∂u/∂x − λ∂v/∂y
    ∂σyy/∂t = −λ∂u/∂x − (λ+2μ)∂v/∂y
    ∂σxy/∂t = −μ(∂u/∂y + ∂v/∂x)
  ```
  with `λ = ρ(c_p² − 2c_s²)`, `μ = ρc_s²`. Steel `ρ = 7850 kg/m³`,
  `c_p = 5960 m/s`, `c_s = 3230 m/s` (Kinsler et al. / NDT handbook, the
  same constants as `phononic::materials` and `sharp_interface::steel`).
  The plane-strain P modulus is `λ+2μ = ρc_p²`, exactly the 1D reduction's
  longitudinal modulus, so the interface-normal P wave coincides with the
  1D rung.

The level set is **tracked but not advected**: the interface is held fixed,
which is second-order accurate at the small amplitudes here.

## 2. Interface coupling

At a face whose two cells are different materials the face state is built
from **one-sided (upwind) reconstructions** of both sides (a 3-point
one-sided interpolation at the face, second order; it falls back to the
cell-centre state when a material layer is thinner than three cells), and
the face-normal state is then solved by the exact 1D solid–fluid Riemann
problem:

* normal velocity `u_n` and normal traction `τ` are continuous and come
  from the star state;
* the inviscid fluid carries no shear traction, so the solid's tangential
  traction is set to zero (`σ_nt = 0`) through the outgoing shear
  characteristic `u_t* = u_t ∓ σ_nt/Z_s`, `Z_s = ρc_s` (a traction-free
  condition for the shear subsystem);
* the fluid tangential velocity is free (slip).

The 1D normal Riemann solve alone cannot represent the tangential (P↔SV)
coupling that makes the solid's normal response angle-dependent. The module
restores that coupling **inside the flux** by propagating the solid's
going P and SV Riemann invariants from the cell centre to the face over
half a cell with the transverse source terms of the 2D characteristic
equations, `S_P = ∂_∥(λu_t − c_pσ_nt)` and `S_SV = ∂_∥(μu_n − c_sσ_tt)`,
where `∂_∥` is the derivative along the face and `σ_tt` the tangential
normal stress (the standard multidimensional characteristic / NSCBC-type
correction). This lets the solid's own P and SV
waves respond to the interface traction instead of imposing a purely normal
P impedance.

Same-material faces use an HLLC flux for the fluid and the exact upwind
(P/S) characteristic flux for the solid. The update is the **unsplit**
`dU/dt = −(F_{i+½}−F_{i−½})/Δx − (G_{j+½}−G_{j−½})/Δy` with one
single-valued flux per face, SSP-RK3 in time and WENO5 primitive
reconstruction (reconstruction is switched off where a stencil crosses a
material, replaced by the one-sided face reconstruction above). CFL:
`Δt = CFL / max_cells[(c+|u|)/Δx + (c+|v|)/Δy]`, `CFL = 0.4`, `c` the
material sound speed (`c_p` for the solid).

## 3. What is conserved and what it costs

For a closed domain the fluid mass, momentum and energy and the solid
momentum telescope exactly (upwind/HLLC fluxes are conservative); the
uniform-flow and stationary-contact tests hold to `< 1e-10` and `< 1e-6`
relative over 50 and 300 steps. The solid total energy is conserved to
round-off apart from the O(amplitude²) normal-momentum mismatch at a moving
mixed face. The `closed_box_energy_balance` fast test (fluid pulse in a
water/steel box with rigid walls, four periods) measures a relative drift of
the total perturbation energy of `2.8e-5`.

Cost: at equal resolution the sharp run uses two independent state vectors
(4-vector fluid, 5-vector solid) and the exact Riemann solve only on the
interface faces; the diffuse five-equation model additionally carries
volume fraction, mixture EOS and (for shear) the full GPR matrix. The
relevant accuracy comparison is the diffuse-interface error itself: in
`analysis/fsi.md` the diffuse model needs `N = 1024` (40 kHz plate, ~6–7
cells/layer) to reach `|T|` absolute error `6.2e-2`, while the sharp 1D rung
reaches `< 1e-2` at 200–800 cells. The 2D module inherits the sharp
interface; its normal-incidence error is `1.3–3.3%` at 37–59 cells per
wavelength (below).

## 4. Exact P/SV reference

`halfspace_reflection(rho_f, c_f, rho_s, c_p, c_s, theta)` returns the
complex reflection coefficient of a plane wave at a fluid–solid half-space
(the exact P/SV solution, Brekhovskikh & Godin):
`Z_in = Z_L cos²(2θ_T) + Z_T sin²(2θ_T)`, `Z_L = ρ_s c_p/cosθ_L`,
`Z_T = ρ_s c_s/cosθ_T`, `Z_f = ρ_f c_f/cosθ`, `R = (Z_in − Z_f)/(Z_in+Z_f)`,
with the cosines continued as `i√(sin²θ−1)` past a critical angle. The
implementation uses complex arithmetic so `|R| = 1` exactly past the shear
critical angle; `tests/sharp_interface2d.rs::halfspace_reflection_reference`
checks

| θ (deg) | \|R\| exact | note |
|--------:|------------:|------|
| 0   | 0.93870 | normal |
| 8   | 0.93773 | subcritical |
| 14.4| —       | longitudinal critical angle |
| 20  | 0.91958 | only S transmitted |
| 27.3| —       | shear critical angle |
| 40  | 1.0     | total reflection |

## 5. Benchmarks

### (b) Normal incidence recovers the 1D result — passes

`normal_incidence_recovers_1d_star`: a uniform normally-incident fluid
state (`dp = 100 Pa`) is driven into a steel half-space; the steady interface
traction is compared with `sharp_interface::exact_star` for the same
one-dimensional states. Measured relative error `5e-5` (fluid) and `5e-5`
(solid). The fluid–solid Riemann solver is also verified against
`acoustic_star` at `dp = ±100 Pa` to `< 1e-9` relative
(`fluid_solid_star_matches_acoustic`).

A stronger, pulse-free normal-incidence benchmark
(`normal_incidence_manufactured`) initialises the exact normal-incidence
standing wave in the fluid and the transmitted P wave in the solid (the
analytic reflection `R = (Z_s−Z_f)/(Z_s+Z_f)`, transmission `1+R`) and
evolves it for one period, measuring the maximum error in a band
`|y−y_if| < 0.5λ` around the interface (the fluid outflow at the far
boundary is excluded):

| cells/λ (water) | rel. error (fluid+ solid) |
|----------------:|--------------------------:|
| 30 | 5.1e-2 |
| 37 | 3.0e-2  (fast test, bound `< 4%`) |
| 59 | 1.3e-2 |

The error falls with resolution; with the one-sided face reconstruction the
normal-incidence manufactured error falls nearly second order (see the
0° row of the convergence table below, order ≈ 1.95).

### (a) Oblique incidence — **manufactured P/SV benchmark now converges**

The old note reported that oblique incidence did not match the exact P/SV
reflection. The cause was the interface closure: the 1D normal Riemann
solve sees only the normal traction, so the solid's input impedance is the
P impedance `ρ_s c_p` independent of the tangential wavenumber, while the
exact impedance `Z_in(θ)` depends on `θ` through mode conversion (§2). The
flux now carries the transverse P/SV characteristic correction, so the
solid's own P and SV wave families respond to the interface traction.

The benchmark is a **manufactured** plane wave: the exact fluid
(incident+reflected) and solid (P+SV) half-space fields from
`halfspace_reflection` are initialised on a periodic-`x` box and evolved,
and the L∞ error in a band `±0.5λ` around the interface is measured at one
interior column (`oblique_field_error`). The exact fields satisfy
normal-velocity continuity, normal-traction continuity and zero shear
traction to machine precision (`oblique_exact_interface_check`, `<1e-12` /
`<1e-9`), so this measures the interface operator directly with no pulse or
FFT extraction. The fast regression `oblique_manufactured` (20°, half a
period) gives `rel = 1.0e-1` at 20 cells/λ in ~14 s. The ignored
`oblique_convergence_sweep` refines each angle over one period:

| θ (deg) | 8 cells/λ | 12 cells/λ | 16 cells/λ | order (12→16) |
|--------:|----------:|-----------:|-----------:|--------------:|
| 0 (normal)             | 4.7e-1 | 2.5e-1 | 1.4e-1 | 1.95 |
| 8 (subcritical)        | 5.2e-1 | 2.7e-1 | 1.5e-1 | 2.04 |
| 20 (P evanescent, S propagating) | 5.0e-1 | 2.6e-1 | 1.7e-1 | 1.37 |
| 32 (beyond S critical, `\|R\|=1`) | 1.21 | 5.9e-1 | 3.8e-1 | 1.50 |

All four angles converge at roughly first-to-second order and — unlike the
old result — the oblique error has the same character as the normal one, so
the P/SV mode conversion is being represented rather than omitted. The
oblique constant is larger than the normal one and grows near the shear
critical angle (evanescent shear, largest sensitivity to the interface
closure).

What is still *not* exact: the mixed-face flux is a normal Riemann solve
plus a **single** transverse-characteristic source step. That correction is
only first-order accurate in `Δx`, so the scheme is consistent but the
converged constant is large at coarse resolution (38% at 16 cells/λ at
32°). A tangential-wavenumber (Dirichlet-to-Neumann) interface operator
would recover the exact oblique impedance; that is the natural next rung.
The fast regression is therefore a manufactured field benchmark rather than
a converged `|R|(θ)` pulse benchmark.

### (c) Curved interface — not attempted

The cylindrical steel shell/rod benchmark against the Faran (1951)
partial-wave series was not attempted in this run.

## 6. Limitations

* Interface held fixed (no advection); valid only at linear amplitude.
* Solid is linear elastic, plane strain, no plasticity or anisotropy.
* The interface closure is a normal Riemann solve plus a first-order
  transverse characteristic correction. It represents P↔SV conversion and
  the manufactured oblique field converges under refinement, but it is not
  an exact tangential-wavenumber (DtN) interface operator, so coarse-grid
  oblique errors are large (up to 38% at 16 cells/λ beyond the shear
  critical angle).
* The oblique benchmark is a manufactured exact field, not a converged
  `|R|(θ)` pulse measurement (the old reference-subtraction sweeps were
  dominated by the interface-band discretisation error and are not used).
* Curved-interface Faran scattering not attempted.

## 7. Files

* `fluid/crates/compressible/src/sharp_interface2d.rs`
* `fluid/crates/compressible/tests/sharp_interface2d.rs`
* `pub mod sharp_interface2d;` in `fluid/crates/compressible/src/lib.rs`
