# Electrohydrodynamics: charged interfaces, drops and electrorheology

Crate: `fluid/crates/electrofluid`. All numbers below come from
`cargo run -p electrofluid --example electrofluid` (release) and are
re-checked by `cargo test -p electrofluid` (9 tests, ~12 s).

This is the electromagnetic branch of `emergence_ladder.md`:

| coarse model | indicator | unfold / closure | where |
|---|---|---|---|
| quasi-static potential (∇×E=0, ∇·J=0) | charge-relaxation time τ_e = ε/σ vs flow time | full Maxwell | this crate's electric solution |
| Gouy–Chapman / Poisson–Boltzmann | u = zeψ/kT; λ_D/h | Poisson–Nernst–Planck (resolved) | `pnp.rs` |
| Debye–Hückel (linear) | \|u\| ≪ 1 | Gouy–Chapman | `pnp.rs` |
| Helmholtz–Smoluchowski slip | λ_D/h | resolved Stokes + PNP body force | `eo.rs` |
| leaky-dielectric drop (Taylor–Melcher) | Ca_E; D ≪ 1 | free-boundary Stokes / full Maxwell | `leaky.rs` |
| point-dipole ER chain | Mason number Mn; φ | particle simulation / Bingham | `er.rs` |

Material constants (in `pnp.rs`): ε₀ = 8.854 187 8128×10⁻¹² F/m,
e = 1.602 176 634×10⁻¹⁹ C, k_B = 1.380 649×10⁻²³ J/K (CODATA 2018);
water ε_r = 78.5 at 298 K (Malmberg & Maryott 1956 give 78.36 at 25 °C;
78.5 is the usual electrokinetic rounding), μ = 8.90×10⁻⁴ Pa·s (CRC).

## 1. Electrokinetics rung: PNP, Gouy–Chapman, Debye layer

At equilibrium the ion fluxes J± vanish, so PNP collapses to the
Poisson–Boltzmann equation for a symmetric z:z electrolyte,
ε ε₀ ψ″ = 2 z e n₀ sinh(zeψ/kT), i.e. u″ = λ_D⁻² sinh u. This is the exact
no-flux steady state of PNP. λ_D = √(ε ε₀ kT / 2 z² e² n₀).

### Debye length (water, 298 K, 1:1)

| c (mM) | λ_D computed | published (1:1, 25 °C) |
|---|---|---|
| 1 | 9.620 nm | ≈ 9.61 nm |
| 10 | 3.042 nm | ≈ 3.04 nm |
| 100 | 0.962 nm | ≈ 0.962 nm |

(The published values are the standard 1:1 table; the computed values
agree to better than 0.1 %.)

### Gouy–Chapman, Debye–Hückel, Grahame

The PB equation integrates once to the Gouy–Chapman profile
tanh(u/4) = tanh(u₀/4) e^{−x/λ_D}, and Gauss's law gives the Grahame
surface charge σ = √(8 n₀ ε ε₀ kT) sinh(zeψ₀/2kT).

Numerical Newton solve of the nonlinear PB equation (the PNP steady state)
on [0, 8λ_D], ψ₀ = −100 mV, c = 1 mM:

| check | value |
|---|---|
| max \|PB_num − Gouy–Chapman\| | 2.59×10⁻⁵ V |
| max \|PB_num − Debye–Hückel\| | 8.89×10⁻³ V (u₀ = 4, so DH is wrong by ~0.5 %) |
| Grahame σ vs σ from the profile | −1.2732×10⁻² vs −1.2731×10⁻² C/m², rel. 3.8×10⁻⁵ |
| small-signal ψ₀ = −5 mV: max \|PB − DH\|/\|ψ₀\| | 3.4×10⁻⁴ |

The indicator of the Gouy–Chapman rung is u = zeψ/kT; of the Debye–Hückel
rung, \|u\| ≪ 1. The small-signal test is the linear-limit confirmation.

## 2. Electro-osmotic fold: Helmholtz–Smoluchowski slip

Resolved rung: steady Stokes with the PB charge as a body force,
−μ u″ = ρ_e E, u = 0 on both walls of a channel of width h, both walls at
ζ. Fold: for λ_D/h → 0 the body force is a boundary layer and the core is a
plug flow with the Helmholtz–Smoluchowski slip
u_s = −ε ε₀ ζ E/μ. The indicator is λ_D/h.

Channel h = 60 nm, ζ = −25 mV, E = 10⁴ V/m, 600 cells:

| c (mM) | λ_D/h | u_center | u_HS | rel. error |
|---|---|---|---|---|
| 0.1 | 0.507 | 7.107×10⁻⁵ | 1.952×10⁻⁴ | 6.4×10⁻¹ |
| 1 | 0.160 | 1.784×10⁻⁴ | 1.952×10⁻⁴ | 8.6×10⁻² |
| 10 | 0.051 | 1.95220×10⁻⁴ | 1.95240×10⁻⁴ | 1.0×10⁻⁴ |
| 100 | 0.016 | 1.952398×10⁻⁴ | 1.952398×10⁻⁴ | 1.5×10⁻¹⁴ |

The fold is clean: the error collapses as (λ_D/h)² once the double layers
no longer overlap, and fails by 64 % at λ_D/h ≈ 0.5. This is the
certificate for replacing the Debye layer by a slip boundary condition.

## 3. Leaky-dielectric (Taylor–Melcher) rung

Quasi-static electrics plus Stokes. With R = σ_i/σ_o, S = ε_i/ε_o the
potential is analytic,
φ_o = (−E₀ r + A/r²)cosθ, φ_i = B r cosθ,
A = E₀ a³(R−1)/(R+2), B = −3E₀/(R+2),
and the interface free charge produces local electric tractions
f_t^E = q_s E_t = −T sinθ cosθ (tangential), T = 9 ε_o ε₀ E₀²(R−S)/(R+2)²,
and [[p^E]] = ½[[ε(E_n²−E_t²)]] (normal). The tangential traction shears the
interface into a Hadamard recirculation; the normal Maxwell-stress jump
enters the normal-stress balance. The two are coupled through the same P₂
Stokes mode, which is implemented in full in `leaky.rs`:

    u_r = U(r) P₂,  u_θ = V(r) P₂',  ee = −T a⁵/[5(μ_o+μ_i)],
    inside : U = ee(r/a⁵ − r³/a⁷),  V = ee(r/(2a⁵) − 5r³/(6a⁷)),
    outside: U = ee(r⁻⁴ − a⁻²r⁻²),  V = −ee r⁻⁴/3,

The sign of `ee` is fixed by the physical viscous traction jump
[[σ_rθ]] = −f_t^E = +T sinθ cosθ. The interfacial velocity is the standard
Taylor circulation

    u_θ(a,θ) = U_T sin2θ,
    U_T = −(9/10) ε_o ε₀ E₀² a /[μ_o(1+λ)] · (R−S)/(R+2)²,  λ = μ_i/μ_o,

i.e. equator→pole for R > S (prolate A) and pole→equator for R < S,
matching the Vlahovska (2019) review (whose printed formula carries a
spurious extra 1/S).

Small-deformation Taylor (1966) result, written in the requested form
D = (9/16) Ca_E Φ/(2+R)², Ca_E = ε_o ε₀ E₀² a/γ:

    Φ = (R² + 1 − 2S) + (3/5) λ_v (R − S),   λ_v = (2μ_o + 3μ_i)/(μ_o + μ_i).

D > 0 is prolate, D < 0 oblate; the sign rule is R > S → prolate,
R < S → oblate. The conductor limit R → ∞ gives D = (9/16)Ca_E. The P₂
normal-stress balance 4γ s/a = [[p^E]] + [[σ_rr^visc]] returns the *shape
amplitude* s in r = a(1 + s P₂):

    s = (3/4) Ca_E Φ/(R+2)²,   so that s/D = 4/3 exactly.

| case | R | S | Ca_E | Φ | D | s = 4D/3 |
|---|---|---|---|---|---|---|
| water in oil | 10⁶ | 40 | 4.25×10⁻³ | ≫1 | +2.39×10⁻³ | +3.19×10⁻³ (prolate) |
| R = 1, S > 1 | 1 | 4 | 4.43×10⁻³ | −10.5 | −2.91×10⁻³ | −3.87×10⁻³ (oblate) |
| R = 1, S < 1 | 1 | 0.5 | 4.43×10⁻³ | +1.75 | +4.84×10⁻⁴ | +6.46×10⁻⁴ (prolate) |
| R = 0.1, S = 1.5 | 0.1 | 1.5 | 4.43×10⁻³ | −4.09 | −2.31×10⁻³ | −3.08×10⁻³ (oblate) |

### The 4/3 is geometry, not missing physics

The earlier version of this crate reported that the P₂ balance came out
exactly 4/3 × Taylor's D and attributed the factor to the *flow driven by
the normal Maxwell stress*. That diagnosis is wrong. For r = a(1 + s P₂),
L/a = 1 + s (pole) and B/a = 1 − s/2 (equator), so

    D = (L−B)/(L+B) = (3s/2)/(2 + s/2) = 3s/(4 + s) → 3s/4.

Taylor's deformation parameter is D = 3s/4 and the stress balance returns s;
the 4/3 is the exact conversion. This is stated explicitly in Vlahovska
(2019, sidebar: "D = 3s/4") and Das & Saintillan (2020, Eq. 5.5). The full
small-deformation problem is already complete: the flow is driven by the
tangential electric traction and the normal Maxwell stress enters the normal
balance. A mean flow driven directly by the normal traction would violate
the steady kinematic condition u_r(a) = 0 and gives the wrong viscosity
structure (e.g. 5μ_o+4μ_i in the denominator instead of 5(μ_o+μ_i)).

### What is checked and what is analytic only

- **Analytic:** potential, circulation, shape amplitude and D.
- **Numerical check:** finite differences of the analytic fields confirm
  (a) velocity continuity across r = a to 2–4×10⁻⁴ of the velocity scale,
  and (b) the physical tangential traction jump [[σ_rθ]] = +T sinθ cosθ to
  1.7–2.3×10⁻³ of T. A test also checks the circulation amplitude against
the published formula U_T above to 10⁻¹² relative, and the flow direction.
- **Finite-Ca / published comparison.** We did not build an independent
  free-boundary boundary-integral solver within the budget. Instead we
  compare with the published axisymmetric BEM of Sorgentone & Vlahovska
  (J. Fluid Mech. 2021, arXiv:2004.05493): for an isolated drop with
  R = 0.1, S = 1, λ = 1, Ca = 0.1 their figure 4(a) reaches
  D ≈ −0.030, while the small-deformation prediction is D_T = −0.0298
  (1% below the BEM value). The experiments of Salipante & Vlahovska
  (Phys. Fluids 22, 112110 (2010)) are closed access here; their data (as
  reproduced in Das & Saintillan, J. Fluid Mech. 880, A22 (2020),
  Table 1 and figure 2) agree with this small-deformation theory at low
  fields and show the expected prolate/oblate transitions. The conductor
  limit and the R=1 cases in the table are exact analytic checks.

## 4. Electrorheology rung: chains, yield stress, Bingham fold

Dielectric spheres in an insulating oil polarise and attract along the
field. We use the exact point-dipole tensor force (equivalently
p = 4π ε₀ ε_f a³ β E, β = (ε_p−ε_f)/(ε_p+2ε_f)),

    F = (3p²/4π ε₀ ε_f r⁴)[ (1−5cos²θ) r̂ + 2 cosθ ŷ ].

The aligned contact force is 1.5π ε₀ ε_f β² E² a². (The frequently quoted
contact form (3/16)π ε₀ ε_f d² β² E², d=2a, is exactly half of this
normalisation; factors of two between conventions are common. In the
conventional half normalisation our measured coefficient below doubles,
to ≈ 0.43, matching Parthasarathy & Klingenberg 1996.) The electrostatic
yield-stress scale is τ_y = C φ ε₀ ε_f β² E².

### Simulation

2D overdamped particles (48 spheres, φ = 0.311), dipole pair force,
steep short-range repulsion plus hard-core projection, Stokes drag,
Lees–Edwards shear. Units: a = 1, F₀ = ε₀ ε_f β² E² a², mobility 1, so
the dimensionless shear γ̇ is related to the standard Mason number
Mn = μ γ̇/(ε₀ ε_f β² E²) = γ̇/(6π). τ is the dipole virial (1/V)Σ F_x Δy.

Zero shear: chain order 1.00, interaction energy U/N = −2.97 (strong
chains; the perfect initial lattice is a stationary saddle, so a small
symmetry-breaking jitter is required).

| Mn | γ̇ | τ_E | U/N | chain order |
|---|---|---|---|---|
| 0.02 | 0.38 | +0.266 | −2.56 | 0.82 |
| 0.05 | 0.94 | +0.272 | −2.05 | 0.52 |
| 0.10 | 1.88 | +0.214 | −1.70 | 0.42 |
| 0.20 | 3.77 | +0.115 | −1.17 | 0.27 |
| 0.50 | 9.42 | +0.063 | −0.75 | 0.29 |
| 1.00 | 18.8 | +0.083 | −0.80 | 0.29 |
| 2.00 | 37.7 | +0.076 | −0.65 | 0.28 |
| 5.00 | 94.2 | +0.075 | −0.78 | 0.26 |

- The low-shear plateau is τ_y ≈ 0.27 (τ_y/φ ≈ 0.87 in our normalisation,
  ≈ 0.43 in the conventional half normalisation) — consistent with the
  point-dipole structure factor C ∼ 0.1–1 (Klingenberg 1993;
  Parthasarathy & Klingenberg 1996; Halsey 1992).
- The chain energy rises from −2.56 to −0.8 as Mn grows from 0.02 to 0.5:
  the chains break and the shear stress falls to the viscous-dominated
  value. The onset is at Mn ≈ 0.1–0.2 in our convention, i.e. within an
  order-one factor of Mn ≈ 1; different authors put the 6π and β
  conventions of the Mason number in different places, so we quote both.
- Continuum Bingham fold: τ = τ_y + η γ̇ for τ > τ_y, and γ̇ = 0 below.
  The indicator is Mn: chains are the coarse model and are valid for
  Mn ≲ 0.1; above that the particle simulation (the unfold) is the honest
  model. The measured minimum separation stayed at 2a throughout, so the
  hard-core projection is well resolved.

## Sources

- G. I. Taylor, *Studies in electrohydrodynamics I*, Proc. R. Soc. A 291,
  159 (1966); J. R. Melcher & G. I. Taylor, Annu. Rev. Fluid Mech. 1, 111
  (1969) — leaky-dielectric potential and deformation.
- D. A. Saville, *Electrohydrodynamics: the Taylor–Melcher leaky
  dielectric model*, Annu. Rev. Fluid Mech. 29, 27 (1997) — discriminating
  function and viscosity factor.
- P. M. Vlahovska, *Electrohydrodynamics of drops and vesicles*, Annu.
  Rev. Fluid Mech. 51, 305 (2019) — sidebar "Drop deformation in a uniform
  DC electric field": D = 3s/4, the discriminating function and the
  interfacial (Taylor circulation) velocity U_T.
- D. Das & D. Saintillan, *Three-dimensional electrohydrodynamic drop
  theory*, J. Fluid Mech. 880, A22 (2020) — steady shape coefficient and
  the normal/tangential stress balance for the complete Melcher–Taylor
  model.
- C. Sorgentone & P. M. Vlahovska, *Electrohydrodynamic interaction of
  drop pairs*, J. Fluid Mech. (2021), arXiv:2004.05493 — axisymmetric BEM
  used for the finite-Ca comparison (R = 0.1, S = 1, λ = 1, Ca = 0.1).
- P. F. Salipante & P. M. Vlahovska, *Electrohydrodynamics of drops in
  strong uniform dc electric fields*, Phys. Fluids 22, 112110 (2010) —
  experiment; closed access here, data as reproduced in Das & Saintillan
  (2020).
- Gouy (1910), Chapman (1913); Grahame (1947) — diffuse double layer and
  surface charge. Debye & Hückel (1923) — linear limit.
- Helmholtz (1879), Smoluchowski (1903) — electro-osmotic slip.
- D. J. Klingenberg, J. Rheol. 37, 199 (1993); M. Parthasarathy &
  D. J. Klingenberg, Mater. Sci. Eng. R 17, 57 (1996); T. C. Halsey,
  Science 258, 761 (1992) — point-dipole force and yield stress.
