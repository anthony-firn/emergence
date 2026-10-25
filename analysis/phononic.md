# Phononic waves: layered transfer matrices, Bloch bands, homogenisation and locally resonant gaps

Implemented in `fluid/crates/phononic`. This is the mechanical counterpart of
`fluid/crates/photonics` and sits below a ray/geometrical-acoustics
description. Geometrical acoustics is exact only once two things are
unresolved: the phase of multiply reflected sound, and the sub-wavelength
structure of the matter. This crate keeps both, then provides the cheap
indicator and the unfold/fold for each. Nothing in `radiation` or `photonics`
is changed.

## The models

**Layered transfer matrices.** Longitudinal waves at normal incidence in a
layer are described by continuous pressure `p` and particle velocity `v`. With
impedance `Z = rho c = sqrt(rho M)`, wavenumber `k = omega sqrt(rho/M)`, and
`theta = k d`, the characteristic matrix and stack solution are

\[
  M_j = \begin{pmatrix}\cos\theta & -i Z \sin\theta\\
  -i \sin\theta/Z & \cos\theta\end{pmatrix},
  \qquad
  Z_{\mathrm{in}} = \frac{m_{22} - m_{12}/Z_s}{-m_{21} + m_{11}/Z_s},
  \qquad
  r=\frac{Z_{\mathrm{in}}-Z_0}{Z_{\mathrm{in}}+Z_0}.
\]

The power transmittance is `T = (Z_0/Z_s)|t|^2`; the particle-velocity
amplitude `t_v = t Z_0/Z_s` makes the lossless identity symmetric,
`|r|^2 + (Z_s/Z_0)|t_v|^2 = 1`. Loss enters through `M = M_0(1 + i eta)`.

**Bloch dispersion.** The product of one unit cell's matrices gives the exact
band structure of the infinite periodic medium through the Rytov relation

\[
  \cos(q\Lambda) = \cos(k_1 d_1)\cos(k_2 d_2)
    - \tfrac12\!\left(\frac{Z_1}{Z_2}+\frac{Z_2}{Z_1}\right)
      \sin(k_1 d_1)\sin(k_2 d_2).
\]

For a quarter-wave stack (`k_i d_i = pi/2` at `f_0`) the edges are at
`sin theta = 2 sqrt(Z_1 Z_2)/(Z_1+Z_2)`, giving the closed-form relative width
`(4/pi) asin((Z_1-Z_2)/(Z_1+Z_2))`.

**Effective medium (fold).** A laminate loaded in series is an isostress
problem: `1/M_eff = sum f_i/M_i`, `rho_eff = sum f_i rho_i`, and the exact
phase velocity of the Bloch wave tends to `c_eff = sqrt(M_eff/rho_eff)` as
`q Lambda -> 0`. The relative error is `O((omega Lambda/c_eff)^2)`.

**Bounds.** For an isotropic two-phase composite the effective bulk and shear
moduli are bracketed by Voigt (arithmetic) and Reuss (harmonic) averages and,
more tightly, by the Hashin–Shtrikman bounds. The stiffer phase is identified
internally, and the exact ordering `Reuss <= HS_lower <= HS_upper <= Voigt` is
tested for both `K` and `G`.

**Locally resonant metamaterial.** A chain of outer masses `m1` on springs
`k1`, each carrying an inner resonator `m2` on a spring `k2`, has the exact
dynamic mass

\[
  m_{\mathrm{eff}}(\omega) = m_1 + \frac{m_2\omega_0^2}{\omega_0^2-\omega^2},
  \qquad \omega_0^2 = k_2/m_2,
\]

which is negative on `(omega_0, omega_0 sqrt(1+m_2/m_1))`. The lattice
dispersion is `cos(qa) = 1 - omega^2 m_eff(omega)/(2 k_1)`, whose gap extends
slightly below `omega_0`: the acoustic branch tops out where
`omega^2 m_eff = 4 k_1`.

**Fluid–solid–fluid.** At oblique incidence the fluid exerts only a normal
traction and the tangential traction vanishes. The plate carries P and SV
potentials with `k_zP = sqrt(k_P^2-k_x^2)`, `k_zS = sqrt(k_S^2-k_x^2)`. The
fluid imposes `Z_f = rho_f c_f/cos theta`, `p + Z_f v_z` fixed by the incident
wave at the front, `p - Z_f v_z = 0` at the back, and `sigma_xz = 0` on both
faces.

Two solvers share that physics. `solve` uses an unconditionally stable
**global (stiffness) matrix** formulation: the four partial-wave amplitudes
(`pf, pb, sf, sb`) are normalised at the face they originate from, so every
exponential factor is `exp(i k_z d)` with `|exp(i k_z d)| <= 1`. This is the
single-layer instance of the Rokhlin–Wang global/stiffness-matrix method; it
conserves energy to round-off past the shear critical angle at any `f d`.
`solve_transfer_matrix` is the historical 4x4 propagator. It advances
`[u_x, u_z, sigma_zz, sigma_xz]` from `z=0` to `z=d` by `E(d) E(0)^{-1}`; that
mix of growing and decaying evanescent partial waves loses precision once
`|Im k_zS| d` is large, so it is kept only as the well-conditioned reference.

**Rayleigh–Lamb free plate.** Setting the four face stresses to zero gives a
4x4 homogeneous system `free_plate_determinant(solid, d, omega, c) = 0` on the
guided modes; splitting by mirror symmetry recovers the classical
Rayleigh–Lamb equations (with `h = d/2`, `delta = k^2 - beta^2`)

\[
  S:\; \delta^2 \sin(\beta h)\cos(\alpha h) + 4\alpha\beta k^2 \sin(\alpha h)\cos(\beta h)=0,
\]
\[
  A:\; 4\alpha\beta k^2 \sin(\beta h)\cos(\alpha h) + \delta^2 \sin(\alpha h)\cos(\beta h)=0.
\]

Below the P and S body-wave speeds `alpha` and `beta` are continued as `i|alpha|`,
`i|beta|`, so `lamb::residual` evaluates the appropriate real combination in
each of the three speed intervals. Cross-multiplying the `tan` ratio forms
introduces removable roots at `c_s` (the `beta -> 0` degeneracy) and at `c_p`
on the antisymmetric branch; those speeds are excluded. A free mode with
`c > c_f` has coincidence angle `asin(c_f/c)`, and the immersed-plate
transmission peaks sit there (shifted by fluid loading).

**2D phononic crystal (plane-wave expansion).** For a square lattice (period
`a`) of circular cylinders (radius `r`, filling fraction `f = pi r^2/a^2`) the
Fourier coefficients of a material function `P(r) = P_m + (P_i-P_m) chi(r)` are
`P_0 = P_m + f (P_i-P_m)` and `P_G = (P_i-P_m) F(G)`, with the disk structure
factor `F(G) = 2 f J_1(Gr)/(Gr)` (Kushwaha et al. 1993). Two polarisations are
implemented.

* **Anti-plane / SH** (`u = w z_hat`): `div(mu grad w) + rho omega^2 w = 0`
  becomes `M_mu (K.K') w = omega^2 M_rho w`. This is **mathematically
  identical to TM photonics** `div((1/eps) grad H) + (omega^2/c^2) H = 0`
  under `mu <-> 1/eps`, `rho <-> 1/c^2`; `tm_frequencies` is the same solver
  with `A = eps_mat/eps_inc`, `B = 1`.
* **In-plane mixed P-SV**: `rho omega^2 u_i = K_j (C_ijkl)_{G-G'} K'_k u_l`
  with the isotropic stiffness, i.e.
  `M_il = mu (K.K') delta_il + lambda K_i K'_l + mu K_l K'_i` and mass
  `rho_{G-G'} delta_il`.

Internally the problem is dimensionless (lengths in `a`, moduli in `mu_matrix`,
densities in `rho_matrix`), so eigenvalues are `Omega^2 = omega^2 rho_m a^2/mu_m`
of order one. The generalised symmetric eigenproblem is reduced by Cholesky and
solved by Householder tridiagonalisation plus the implicit QL algorithm
(`symmetric_eigenvalues`), with a Jacobi implementation used as the reference
in the unit tests.

## Benchmarks (measured)

Run `cd fluid && cargo test -p phononic`. Every value is printed by the tests;
"error" is the deviation from the exact or published result.

| Benchmark | Exact / published | Computed | Error |
|---|---|---|---|
| Steel plate in water, `T=1` at `d=n lambda/2`, n=1..5 | 1 | 1.000000000000 | < 1e-12 |
| Off-resonance plate transmission | < 1 | 0.00795 | — |
| Lossless energy `R + (Z_out/Z_in)|t_v|² = 1`, 4-layer stack, 200 freqs | 1 | 1 | 4.4e-16 |
| Lossy stack `R+T+A=1` (steel, loss 0.01) | 1 | 1 | round-off; `A=4.996e-3` |
| Quarter-wave edges `(4/π)asin((Z1−Z2)/(Z1+Z2))` vs Rytov scan (Z=2:1) | 0.783653–1.216347 MHz | 0.783654–1.216346 MHz | 1e-6 MHz |
| `cos(qΛ)` at the analytic edges | −1 | −1 | 1e-10 |
| Bloch gap vs finite dip, steel/Al quarter-wave, 20 periods | 0.695542–1.304458 MHz | 0.711880–1.288120 MHz | centre 0; width 5.4% narrower |
| Bloch `v_phase -> c_eff` (Reuss), `ωΛ/c_eff = 0.4 -> 0.0125` | `c_eff=2598.076 m/s` | quadratic | `err/x² -> 1.042e-2` |
| Hashin–Shtrikman ordering, K and G, 5 fractions | Reuss ≤ HS≤ Voigt | holds | exact |
| HS end members `f=0,1` | `K_2`, `K_1` | equal | < 1e-6 |
| Mass-in-mass effective mass vs exact 2-DOF eigenproblem, 201 q | equal | equal | 4.1e-15 |
| Negative-mass interval inside true gap | `(1,1.4142)` ⊂ `(0.8740,1.4142)` | holds | exact |
| Sub-wavelength gap, `k2=0.01` | `lambda/a` large | 44.4 | — |
| Fluid–plate at `θ=0` vs longitudinal layered, 0.1–10 MHz | equal | equal | 2.3e-14 |
| Fluid–plate energy `R+T=1`, `θ in [0,40°]`, 8001 angles, 2 MHz | 1 | 1 | 2.6e-14 |
| Fluid–plate energy beyond shear critical, 1, 3, 5, 10 MHz / 5 mm steel | 1 | 1 | 4.8e-14 worst (stable) |
| Old transfer matrix, same 3, 5, 10 MHz scans | 1 | 1 | 2.6e2, 6.8e23, 1.2e79 (fails) |
| Stable vs transfer matrix, 0.5 MHz / 5 mm, `θ in [0,40°]` | equal | equal | 5.7e-13 |
| P coincidence, water–steel | 14.398° | peak 14.19° | 0.21° |
| S0 (extensional) coincidence | 15.841° | peak 16.36° | 0.52° |
| S0 thin plate `sqrt(E/(rho(1-nu²)))`, fd 0.10 / 0.25 MHz·mm | 5429.1 m/s | 5428.6 / 5425.8 | 0.01% / 0.06% |
| A0 flexural limit, fd 0.10 / 0.25 MHz·mm | 0 (as fd->0) | 958.8 / 1444.2 | shrinking |
| `lamb::residual` vs `free_plate_determinant`, 12 modes at 3 MHz / 5 mm | both zero | worst ratio 1.5e-10 | 1.5e-10 |
| Leaky Lamb peaks at `asin(c_f/c)`, f=2,3 MHz, modes `c>c_s` | coincidence | local max within 1.4° | T>0.999 |
| SH vs TM photonic, `mu=1/eps`, `rho=1`, m=3 | identical | identical | 0 |
| PWE empty lattice vs `c|k+G|` (folded free phonon) | exact | exact | 0 |
| In-plane PWE empty lattice vs folded Navier `c_p,c_s\|k+G\|` | exact | exact | 4.3e-10 rad/s |
| SH low-frequency `c` vs Reuss–Voigt bracket | inside | inside | — |
| SH `c` vs 2D Maxwell–Garnett/Hashin, m=5->11 | 963.7 m/s | 1068.8 -> 1009.8 | 10.9% -> 4.8% |
| SH complete gap, steel/epoxy f=0.4, extrapolated | see text | `Ω=2.30–5.85` | width/mid 0.87 |
| In-plane complete gap, steel/epoxy f=0.4, m=3 / 5 | see text | `Ω=4.80–10.68 / 4.31–9.42` | — |

The half-wave, energy, edge, ordering, eigenproblem and normal-incidence rows
are exact identities, so their errors are round-off. The quarter-wave edge
formula is the exact all-orders result for a lossless quarter-wave stack (not a
leading-order approximation), so the 1e-6 MHz is scan resolution. The finite
dip is 5.4% narrower than the infinite-stack gap because 20 periods is finite
(the field is not yet fully evanescent at the gap edge); the centres coincide.

### Effective medium: measured convergence

Laminate layers `rho1=3000, c1=4000` and `rho2=2000, c2=2000` with quarter-wave
thicknesses, `Lambda` the period. Exact Bloch phase velocity from the Rytov
relation against `c_eff`:

| `x = ωΛ/c_eff` | relative error | error / x² |
|---|---|---|
| 0.4 | 1.69e-3 | 1.057e-2 |
| 0.2 | 4.18e-4 | 1.046e-2 |
| 0.1 | 1.04e-4 | 1.043e-2 |
| 0.05 | 2.60e-5 | 1.042e-2 |
| 0.025 | 6.51e-6 | 1.042e-2 |
| 0.0125 | 1.63e-6 | 1.042e-2 |

The error falls by 4x each time `x` halves: clean `O(x²)`. The measured
prefactor is `1.04e-2` for this contrast; `homogenize::homogenisation_error`
uses the conservative envelope `(ωΛ/c_eff)²/12 = 8.3e-2 x²`, which bounds the
measured error by a factor 8. It is a switch, not a certificate.

### The sub-wavelength gap

For `m1=m2=1`, `k1=1`, `k2=0.01` (units set the frequencies) the inner
resonance is `omega_0 = 0.1`, `m_eff < 0` on `(0.1, 0.1414)`, and the true
lattice gap is `(0.099875, 0.141421)`. The gap wavelength is `44.4` lattice
constants, while the bare-chain Bragg frequency is `20x` the gap bottom: the
gap is set by the resonator, not by the lattice period. At zero frequency
`m_eff -> m1+m2` to 1e-6, which is the folded static chain.

### Lamb modes and leaky-wave coincidence

The free-plate Rayleigh–Lamb modes are found from `lamb::residual` (the
pole-free symmetric/antisymmetric equations) and cross-checked against
`fluid_coupling::free_plate_determinant`, which comes from the independent
global (stiffness) matrix. At 3 MHz / 5 mm steel every one of the 12 modes has
`|det|` a factor `1.5e-10` smaller than at `c(1 +- 0.005)`. Thin-plate limits:
`S0 -> 5429 m/s` (`sqrt(E/(rho(1-nu^2)))`) to within 0.06%, and `A0` falls
from 1444 m/s at `f d = 0.25 MHz·mm` to 959 m/s at 0.10 MHz·mm, toward the
flexural zero.

The water-loaded plate radiates at the coincidence angle `asin(c_f/c)` of each
mode faster than `c_s` (below the shear critical angle). At 2 and 3 MHz every
such free mode has a transmission maximum `T > 0.999` within 1.4° of its
coincidence angle (typically 0.0–0.2°), the residual gap being the fluid-loading
shift. Modes slower than `c_s` lie past the shear critical angle and produce no
radiating peak, which is the expected near-total reflection.

### 2D phononic crystals: PWE

The scalar (SH) solver is verified exactly against three independent results:
(1) the algebraic identity with TM photonics under `mu <-> 1/eps`,
`rho <-> 1/c^2` (`max|SH - TM| = 0`); (2) the homogeneous limit, where the
folded free-phonon spectrum `c|k+G|` is recovered to machine precision; and
(3) the low-frequency effective medium.

For the anti-plane effective medium, the long-wavelength slope must approach
the published 2D Maxwell–Garnett / Hashin result
`mu_eff = mu_m + f/(1/(mu_i-mu_m) + (1-f)/(2 mu_m))`. Square steel/epoxy,
`f = 0.4` (`c_MG = 963.74 m/s`, `rho_eff = 3848 kg/m^3`):

| `m` | `c` [m/s] | `mu_eff` [Pa] | `|c - c_MG|/c_MG` |
|---|---|---|---|
| 5 | 1068.80 | 4.3957e9 | 1.090e-1 |
| 7 | 1036.63 | 4.1351e9 | 7.564e-2 |
| 9 | 1021.19 | 4.0128e9 | 5.961e-2 |
| 11 | 1009.82 | 3.9240e9 | 4.782e-2 |

Every value lies between the rigorous Reuss (`824.0 m/s`) and Voigt
(`2959.9 m/s`) bounds of `homogenize`; at the tested plane-wave counts the 3D
Hashin–Shtrikman shear bracket (`987.5–2422.6 m/s`) also contains them. The
error to the 2D value falls monotonically with `m`, but only as the inverse
power of the plane-wave count: this is the documented slow PWE convergence for
high contrast.

**First SH gap, square steel/epoxy, `f = 0.4`.** The complete gap of the lowest
two bands (both edges at the `X` point) converges from above as the plane-wave
count increases. `Omega = omega a / c_T(epoxy)` and reduced frequency
`f a/c_T`:

| `m` | `Omega` lower | `Omega` upper | `f a/c_T` lower | `f a/c_T` upper |
|---|---|---|---|---|
| 3 | 2.9014 | 6.9390 | 0.4618 | 1.1044 |
| 5 | 2.6042 | 6.4489 | 0.4145 | 1.0264 |
| 7 | 2.5182 | 6.2786 | 0.4008 | 0.9993 |
| 9 | 2.4719 | 6.1800 | 0.3934 | 0.9836 |
| 11 | 2.4415 | 6.1209 | 0.3886 | 0.9742 |
| extrap. `1/m` | 2.3043 | 5.8549 | 0.3668 | 0.9319 |

The extrapolated `1/m` fit from the two finest runs agrees with the fit from
`m=7,9` to 0.3%, giving a gap width-to-midgap ratio of 0.87. The mixed in-plane
P-SV operator has a gap between bands 2 and 3 over the high-symmetry points
`Gamma`, `X`, `M` at `Omega = 4.8010–10.6784` (`m=3`) and `4.3128–9.4164`
(`m=5`), again moving downward with `m`; a full-BZ scan would be needed to
certify it as complete.

*What we compare and what we do not.* Vasseur et al., *J. Phys.: Condens.
Matter* **6**, 8759 (1994), and Kushwaha et al., *Phys. Rev. Lett.* **71**,
2022 (1993), established complete elastic band gaps in square fibre composites
including steel/epoxy. We reproduce that system at filling fraction 0.4 rather
than a digitised point from a figure. The original paper is not available in
this sandbox, so no published edge is quoted as a hard number; the benchmark is
the existence and breadth of the gap and its convergence. The slow convergence
is the headline: between `m=9` and `m=11` the SH edges still move by ~1%, so any
comparison to a converged band structure must state the plane-wave count. The
reliable, cheap check is the folded effective medium, which needs only a few
waves.

## Fold/unfold map

The organising idea is `emergence_ladder.md`. The acoustics instance:

| Coarse model (fold) | Cheap indicator | Unfold when | Where |
|---|---|---|---|
| Homogeneous laminate (Reuss) | `(ωΛ/c_eff)²`; measured prefactor ~1e-2, envelope `x²/12` | exceeds target phase accuracy | `homogenize` |
| Static heavy chain `m1+m2` | `ω/ω₀` | drive approaches the inner resonance | `resonant` |
| Infinite-crystal band structure (one unit cell) | `N |Im q| Λ` (periods × decay/period) | small: finite size dominates | `layered` |
| Plane-wave crystal band structure | number of plane waves `m` (edge drift) | edges still move: refine `m` or extrapolate | `pwe` |
| Near-total reflection beyond the shear critical angle | `|Im k_zS| d` | never: the global (stiffness) solver is unconditionally stable | `fluid_coupling` (fixed) |

### Hand-off in practice

- **Fold a laminate up** when `ωΛ/c_eff` is below the required accuracy. The
  exact long-wavelength modulus is the Reuss average; the Bloch calculation is
  only needed once the error indicator rises.
- **Unfold a resonant chain** once `ω` nears `omega_0`: the static mass is a
  cheap model only well below the gap.
- **Fold an infinite stack down** to one unit cell: the Rytov relation costs
  one cell, while a finite stack needs `N` cells; the decay indicator
  `N|Im q|Λ` says when the finite stack behaves like the infinite one.
- **Unfold to the global (stiffness) matrix** for every fluid–plate angle. The
  bounded-exponential normalisation removes the conditioning limit; the old
  4x4 transfer matrix is retained only as a well-conditioned cross-check.
- **Refine the plane-wave count** until the gap edges stop moving (or
  extrapolate in `1/m`) before trusting a phononic band structure; the
  folded effective medium needs only a few waves.

## Limits and honest misses

* The transfer matrix is exact only for planar, isotropic, perfectly bonded
  layers; no anisotropy, roughness, interface slip or moving media.
* The historical fluid–plate 4x4 transfer matrix still mixes growing and
  decaying evanescent solutions, but it is no longer used: `solve` is the
  global (stiffness) matrix, whose energy error is round-off at every angle
  and `f d` tested (worst 4.8e-14 up to 10 MHz / 5 mm), while the old solver
  reaches O(1e2)–O(1e79). `solve_transfer_matrix` exists only as the
  well-conditioned reference (agreement 5.7e-13 at 0.5 MHz).
* The homogenisation indicator is an empirical quadratic envelope from the
  benchmark, not a proven bound; it can be exceeded near a constituent
  resonance or at very high contrast.
* Hashin–Shtrikman bounds assume isotropic two-phase composites with no
  percolation, no interfacial resistance and no size effects; the 3D
  `homogenize` bounds are not the 2D plane-strain bounds, so for the SH
  anti-plane slope the comparison uses the rigorous Reuss–Voigt bracket and
  the published 2D Maxwell–Garnett/Hashin value (the measured slope sits just
  above the 3D lower bound because it is still converging from above).
* The PWE converges slowly for high elastic contrast (discontinuous
  displacement, `1/G` Fourier decay): steel/epoxy at `f=0.4` still moves ~1%
  between `m=9` and `m=11`, so the reported gap is a bounded extrapolation, not
  a converged number. Only a square lattice of circular isotropic inclusions
  is implemented; no Fourier factorisation / inverse-rule acceleration.
* The locally resonant model is a discrete mass-spring chain; it is not a
  viscoelastic continuum and has no damping.
* No anisotropy (hexagonal/cubic Christoffel), no 3D crystals, no guided-wave
  solver, no nonlinear acoustics and no thermoelastic coupling.

## References

* L. M. Brekhovskikh and O. A. Godin, *Acoustics of Layered Media* (transfer
  matrices, Rytov dispersion).
* W. T. Thomson, "Transmission of elastic waves through a stratified solid",
  J. Appl. Phys. **21**, 89 (1950); N. A. Haskell, "The dispersion of surface
  waves on multilayered media", Bull. Seismol. Soc. Am. **43**, 17 (1953).
* Z. Hashin and S. Shtrikman, "A variational approach to the theory of the
  elastic behaviour of multiphase materials", J. Mech. Phys. Solids **11**, 127
  (1963).
* Z. Liu, X. Zhang, Y. Mao, Y. Y. Zhu, Z. Yang, C. T. Chan and P. Sheng,
  "Locally resonant sonic materials", Science **289**, 1734 (2000).
* L. E. Kinsler, A. R. Frey, A. B. Coppens and J. V. Sanders, *Fundamentals of
  Acoustics*, 4th ed. (steel, aluminium and PMMA P/S speeds).
* *CRC Handbook of Chemistry and Physics* (water density and sound speed at
  20 C).
* R. D. Mindlin and H. H. Fu, "On the propagation of flexural waves in an
  elastic plate" (Rayleigh–Lamb branches).
* K. F. Graff, *Wave Motion in Elastic Solids* (Dover, 1975), ch. 6, and
  J. L. Rose, *Ultrasonic Waves in Solid Media* (Cambridge, 1999), ch. 8
  (Rayleigh–Lamb symmetric/antisymmetric equations).
* S. I. Rokhlin and L. Wang, "Stable recursive algorithm for elastic wave
  propagation in layered anisotropic media: stiffness matrix method",
  J. Acoust. Soc. Am. **112**, 822 (2002) (global/stiffness matrix).
* M. S. Kushwaha, P. Halevi, L. Dobrzynski and B. Djafari-Rouhani, "Acoustic
  band structure of periodic elastic composites", Phys. Rev. Lett. **71**, 2022
  (1993); M. S. Kushwaha et al., "Theory of acoustic band structure of periodic
  elastic composites", Phys. Rev. B **49**, 2313 (1994) (PWE, disk structure
  factor).
* J. O. Vasseur, B. Djafari-Rouhani, L. Dobrzynski, M. S. Kushwaha and
  P. Halevi, "Complete acoustic band gaps in periodic fibre reinforced
  composite materials", J. Phys.: Condens. Matter **6**, 8759 (1994)
  (square steel/epoxy complete gaps).
* Z. Hashin, "On the elastic behaviour of fibre-reinforced materials",
  J. Appl. Mech. **32**, 143 (1965), and J. C. Maxwell Garnett, Phil. Trans. R.
  Soc. A **203**, 385 (1904) (2D anti-plane effective shear modulus).
* M. Abramowitz and I. A. Stegun, *Handbook of Mathematical Functions*,
  eq. 9.4.7 (large-argument `J_1`).
