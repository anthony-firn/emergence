# Anisotropic optics: the Berreman 4x4 rung, liquid crystals, and where it folds

Implemented in `fluid/crates/photonics/src/berreman.rs`, benchmarked in
`fluid/crates/photonics/tests/berreman.rs` (also short unit tests in the
module). Run `cd fluid && cargo test -p photonics`.

The wave-optics rung in [photonics.md](photonics.md) assumes isotropic
layers: a 2x2 characteristic matrix per polarisation, `TE` and `TM`
independent. Real crystals and liquid crystals are anisotropic, so the two
polarisations couple and the layered mode functions are the four components
of the tangential field. The Berreman 4x4 matrix is that exact local model.
It is the fine rung; the two coarse closures below it are the
sub-wavelength periodic medium and the slowly-varying director, each with a
cheap validity indicator and a measured convergence.

## The model

Maxwell's equations in a uniform medium are written for the tangential state

```text
  psi = (E_x, E_y, H_x, H_y)^T          (H scaled by the vacuum impedance)
  d psi / d z = i k0 Delta psi,          k0 = 2 pi / lambda,
```

with conserved tangential wavevector `s = n0 sin(theta0)` and complex
relative permittivity tensor `eps_ij` (a scalar `mu` is also carried). With
`eps = e`, the Berreman matrix is

```text
  row for dE_x/dz : [ -s e_zx/e_zz,  -s e_zy/e_zz,  0,  mu - s^2/e_zz   ]
  row for dE_y/dz : [  0,             0,           -mu, 0               ]
  row for dH_x/dz : [ -(e_yx - e_yz e_zx/e_zz), -(e_yy - e_yz e_zy/e_zz - s^2/mu), 0, s e_yz/e_zz ]
  row for dH_y/dz : [  e_xx - e_xz e_zx/e_zz,   e_xy - e_xz e_zy/e_zz, 0, -s e_xz/e_zz ]
```

The layer matrix is `P = expm(i k0 d Delta)` (scaling-and-squaring with a
truncated Taylor series), the stack matrix is the ordered product in the
front-to-back sense `psi(d) = M psi(0)`, and the stack is matched to the
incident, reflected and exit plane-wave modes. The unknowns are the complex
reflection and transmission amplitudes for an incident `(E_x, E_y)`; the
result is a full 2x2 complex **Jones** matrix. Reflectance and transmittance
use the normal Poynting flux weights `Re(q)` (`s`) and `Re(n^2/q)` (`p`), so
`R + T = 1` for every lossless stack. `eps` may be uniaxial, biaxial, or
gyrotropic (`eps_xy = -i g`), arbitrarily oriented.

The isotropic path is recovered exactly: the `eps = n^2 I` matrix is block
diagonal, and the 4x4 reduces to the crate's 2x2 TE/TM calculation to
round-off.

## Benchmarks (measured)

| Benchmark | Exact / published | Computed | Error |
|---|---|---|---|
| Isotropic stacks vs `multilayer` TE/TM, 0–69°, 3 index sets | 2x2 transfer matrix | 4x4 Berreman | \(7.9\times10^{-15}\) |
| Layer eigen-phase (in-plane axis, normal incidence) | `Gamma = 2 pi dn d / lambda` = 1.827835726 rad | 1.827835726 rad | \(<10^{-9}\) |
| Quarter-wave plate, 45° linear in | \(\|S_3\|/S_0 = 1\) | 0.991351 | 0.86% (Fresnel ellipticity) |
| Half-wave plate, 45° linear in | 135° | 135.95° | 0.95° |
| Uniaxial interface TE, `n_e` (axis normal to the plane) | Fresnel | Berreman | \(2.7\times10^{-15}\) |
| Uniaxial interface TM, `n_o` | Fresnel | Berreman | \(5.6\times10^{-16}\) |
| Uniaxial Brewster `tan theta_B = n_o/n0` | `R_p = 0` | \(8.9\times10^{-33}\) | round-off |
| Uniaxial TIR beyond `theta_c` (`n_o` and `n_e`) | \(R = 1\) | 1.000000 | \(<10^{-13}\) |
| Energy conservation, lossless uniaxial/biaxial/TN/cholesteric | `R + T = 1` | 1 | \(1.5\times10^{-14}\) |
| Faraday rotation, `eps_xy = -i g`, index-matched | \((k_0/2)(n_+ - n_-) d\) = 2.481513e-2 rad | 2.481512e-2 rad | \(3.6\times10^{-7}\) rel |
| TN vs Gooch–Tarry, parallel polarisers | 0.039439 | 0.039304 | \(1.3\times10^{-4}\) |
| TN vs Gooch–Tarry, crossed polarisers | 0.960561 | 0.960555 | \(6.2\times10^{-6}\) |
| Cholesteric band centre `lambda_0 = n_bar p` (E7) | 500.0 nm | 500.0 nm | design |
| Cholesteric band width `Delta lambda = dn p` (de Vries) | 69.0 nm | 73.0 nm | 5.8% |
| Cholesteric in-band, co-handed / opposite-handed | 1 / 0 | 0.9988 / 0.0014 | — |
| Cholesteric out-of-band (`lambda_0 +/- 2.5 Delta lambda`) | 0 | 0.038 / 0.0005 | — |

The isotropic-reduction, Fresnel, Brewster, TIR, Faraday and energy checks
are exact identities, so their errors are round-off. The quarter/half-wave
ellipticity is the physical ordinary/extraordinary Fresnel amplitude
difference at the air interfaces (the intrinsic retardance test is exact);
it is reported rather than removed by an idealised model. The TN and
cholesteric finite-size deviations are discussed below.

### Brewster and total internal reflection for a uniaxial medium

With the optic axis normal to the plane of incidence the two polarisations
decouple into two effective isotropic media: `TE` (E along the optic axis)
sees `n_e`, `TM` (E in the plane) sees `n_o`. Setting the exit index equal to
the index seen by the polarisation under test makes the back interface
transparent, so the total reflectance equals the front Fresnel reflectance
exactly. For a `1.8 -> 1.5/1.7` interface that gives the analytic
`theta_B = atan(n_o/n0) = 39.8056°` with `R_p = 8.9e-33`, and the critical
angles `asin(n_o/n0) = 56.4427°` and `asin(n_e/n0) = 70.8119°`, with
`R = 1` beyond.

### Twisted nematic and Gooch–Tarry

A `90°` twisted nematic of 5CB (`n_o = 1.532`, `n_e = 1.717` at 589 nm,
S. T. Wu, Phys. Rev. A **33**, 1270 (1986)), `d = 5 µm`, at 550 nm, is
discretised into 400 thin uniaxial layers and compared with the analytic
coupled-mode (Mauguin) result

```text
  J = R(phi) [ cos x - i (Gamma/(2x)) sin x ,  (phi/x) sin x
               -(phi/x) sin x,  cos x + i (Gamma/(2x)) sin x ]
  x = sqrt(phi^2 + (Gamma/2)^2),   Gamma = 2 pi dn d / lambda,
```

whose polariser transmissions are the Gooch–Tarry formula
(Gooch & Tarry, J. Phys. D **8**, 1575 (1975)). At these parameters the
Mauguin number `u = pi dn d/(phi lambda) = 3.64`, and the crossed
transmission is `0.9606` (the polarisation is guided by the twist) against
`0.0394` parallel. The exact Berreman stack reproduces the formula to
`6e-6` (crossed) and `1.3e-4` (parallel); the residual is the
slowly-varying-envelope approximation in the analytic formula.

### Cholesteric selective reflection (de Vries)

A cholesteric with pitch `p` and E7 indices (`n_o = 1.5217`,
`n_e = 1.7472`; S. T. Wu, U. Efron, L. D. Hess, Appl. Opt. **23**, 3911
(1984)) at normal incidence reflects one circular handedness in the band

```text
  lambda_0 = n_bar p,   Delta lambda = dn p,   n_bar = (n_e + n_o)/2.
```

For `p = 305.9 nm` the de Vries band is `lambda_0 = 500.0 nm`,
`Delta lambda = 69.0 nm`. The 20-turn Berreman stack gives in-band
`R = 0.9988` for the co-handed circular input and `0.0014` for the other
handedness, and the measured full width at `R = 0.5` is `73.0 nm`, 5.8%
above de Vries. The excess is the finite-length (Fabry–Perot) narrowing and
the staircase discretisation of the helix (40 sublayers per turn); it
decreases with more turns. The opposite handedness is transmitted
(`R = 0.0014`) across the band.

## The fold/unfold map

The organising idea is [emergence_ladder.md](emergence_ladder.md). The
anisotropic instance:

| Coarse model (fold) | Cheap indicator | Unfold when | Where |
|---|---|---|---|
| Homogeneous uniaxial/effective medium of a short-pitch cholesteric | `p/lambda`; reflectance error `O((p/lambda)^2)` | indicator exceeds target reflectance accuracy | `berreman::cholesteric` vs `effective_index_indicator` |
| Adiabatic (Mauguin) waveguiding by a slowly varying director | inverse Mauguin number `phi lambda/(pi dn d)`; loss `1 - abs(<E_out,n_out>)^2` | indicator not small: polarisation no longer follows the director | `twisted_nematic` vs `twisted_nematic_jones` |
| Isotropic 2x2 TE/TM | anisotropy `abs(eps_par - eps_perp)/abs(eps)` | non-zero: polarisations couple | `multilayer` vs `berreman` |

### Cholesteric homogenisation: measured convergence

Fix the total thickness at `h = 20 µm` and shrink the pitch `p`, so the
number of periods grows while the folded medium stays the same. At normal
incidence the helical director averages to an in-plane permittivity
`(n_o^2 + n_e^2)/2` and an axial `n_o^2`, so the effective index is
`n_eff = sqrt((n_o^2 + n_e^2)/2) = 1.63833`; the incident and exit media are
matched to it. The reflected power (the validity indicator being `p/lambda`)
is:

| `p` (nm) | `p/lambda` | `R` | rotatory power `rho` (rad/m) |
|---|---|---|---|
| 160 | 0.291 | \(1.19\times10^{-3}\) | \(2.02\times10^{-10}\) |
| 80 | 0.146 | \(8.38\times10^{-5}\) | \(2.53\times10^{-11}\) |
| 40 | 0.0727 | \(1.58\times10^{-5}\) | \(3.16\times10^{-12}\) |
| 20 | 0.0364 | \(3.68\times10^{-6}\) | \(3.95\times10^{-13}\) |

From 80 nm to 40 nm the ratio is `0.188` and from 40 nm to 20 nm `0.234`,
i.e. `R ~ (p/lambda)^{2.1-2.4}`: clean quadratic convergence to the folded
effective medium. The folded model used for comparison is the uniaxial slab
`eps_par = (n_o^2+n_e^2)/2`, `eps_perp = n_o^2`; at `p = 100 nm` the fine
stack differs from it by at most `4.6e-4` in reflectance over 0–23°.

The residual optical activity (chirality) survives the fold: the
de Vries/Belyakov rotatory power `rho ~ pi dn^2 p^3/(4 n_bar^2 lambda^2)`
falls as `p^3` and is included as a diagnostic, but the leading reflection
error is the even `O((p/lambda)^2)` term. `effective_index_indicator =
p/lambda` is an empirical switch, not a proven bound, exactly as the
isotropic laminate indicator in [photonics.md](photonics.md).

### The other fold: the adiabatic (Mauguin) limit

A slowly varying director guides the polarisation: the transmitted power on
the output director approaches one as the Mauguin number grows. Measured
loss `1 - |<E_out|n_out>|^2` for a `90°` 5CB cell:

| `u = pi dn d/(phi lambda)` | measured loss | analytic (coupled mode) |
|---|---|---|
| 1.35 | \(8.38\times10^{-2}\) | \(8.43\times10^{-2}\) |
| 2.69 | \(1.17\times10^{-1}\) | \(1.16\times10^{-1}\) |
| 5.38 | \(1.82\times10^{-2}\) | \(1.81\times10^{-2}\) |
| 10.76 | \(7.89\times10^{-3}\) | \(7.82\times10^{-3}\) |
| 21.53 | \(9.11\times10^{-4}\) | \(9.07\times10^{-4}\) |

The exact Berreman stack tracks the analytic coupled-mode loss to ~1%; the
inverse Mauguin number `1/u` is the cheap indicator. When it is not small
the adiabatic closure fails and the exact 4x4 stack must be used — this is
the same "fold where the indicator certifies accuracy" rule as everywhere
else in this repository.

## Limits and honest misses

* The model is planar and static: no roughness, no moving media, no
  nonlinear response, no spatial dispersion (the cholesteric is resolved as
  a stack, not by an effective optically active constitutive law).
* The cholesteric band width is 5.8% high at 20 turns; the discrepancy is
  finite-length plus staircase discretisation, not a model constant. More
  turns reduce it, at linear cost in the layer count.
* The quarter/half-wave tests use a real plate in air, so the Fresnel
  amplitude difference leaves a 0.86% ellipticity for the quarter-wave and a
  0.95° rotation error for the half-wave. The intrinsic retardance
  `Gamma = 2 pi dn d/lambda` is checked exactly via the layer eigenmodes.
* The Mauguin/Gooch–Tarry analytic matrix is the slowly-varying-envelope
  solution; the Berreman stack is the exact Maxwell solution. Where they
  differ, the 4x4 result is the reference.
* The exit medium of the general stack is isotropic. The Brewster/TIR
  benchmarks exploit the decoupling of the two polarisations with the optic
  axis normal to the plane of incidence, which makes the backward interface
  transparent without needing a general anisotropic half-space.
* The gyrotropic form is the standard `eps = eps0 I - i g [a]_x`
  (Faraday/Voigt); no natural optical activity tensor and no magnetic
  permeability tensor are carried.

## References

* D. W. Berreman, "Optics in stratified and anisotropic media: 4x4-matrix
  formulation", J. Opt. Soc. Am. **62**, 502 (1972),
  doi:10.1364/JOSA.62.000502.
* Hl. de Vries, "Rotatory power and other optical properties of certain
  liquid crystals", Acta Cryst. **4**, 219 (1951),
  doi:10.1107/S0365110X51000751 (selective reflection `lambda_0 = n_bar p`,
  `Delta lambda = dn p`, rotatory power).
* C. H. Gooch and H. A. Tarry, "The optical properties of twisted nematic
  liquid crystal structures with twist angles <= 90°", J. Phys. D **8**, 1575
  (1975), doi:10.1088/0022-3727/8/13/007 (transmission formula).
* S. T. Wu, U. Efron, L. D. Hess, "Birefringence measurements of liquid
  crystals", Appl. Opt. **23**, 3911 (1984), doi:10.1364/AO.23.003911 (E7 at
  589 nm: `n_e = 1.7472`, `n_o = 1.5217`).
* S. T. Wu, "Birefringence dispersions of liquid crystals", Phys. Rev. A
  **33**, 1270 (1986), doi:10.1103/PhysRevA.33.1270 (5CB at 589 nm:
  `n_e = 1.717`, `n_o = 1.532`).
* P. G. de Gennes and J. Prost, *The Physics of Liquid Crystals*, 2nd ed.,
  Oxford University Press (1993) (Mauguin/adiabatic waveguiding, cholesteric
  optics).
* P. Yeh and C. Gu, *Optics of Liquid Crystal Displays*, Wiley (1999)
  (twisted-nematic Jones matrix).
* M. Born and E. Wolf, *Principles of Optics*, 7th ed., Cambridge University
  Press (1999) (Stokes parameters and circular-polarisation convention).
* A. Yariv and P. Yeh, *Optical Waves in Crystals*, Wiley (1984)
  (anisotropic plane waves, uniaxial Brewster and total reflection).
