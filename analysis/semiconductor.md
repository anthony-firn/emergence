# Silicon and transistors: results

Crate: `fluid/crates/semiconductor`. Every number below comes from
`examples/semiconductor.rs` and can be rerun (`semiconductor <case>`).

This is the lower part of the electronics ladder in
`emergence_ladder.md`:
- **device (TCAD):** drift–diffusion;
- **circuit:** compact models folded from TCAD, each with an error
  certificate;
- **gate level:** inverters, with noise margins as the validity indicator
  of the Boolean abstraction.

Silicon as a material (melting) is in `compressible.md` §15.

## 1. Silicon (`semiconductor silicon`)

Material model, with sources in `silicon.rs`:

| quantity | source |
|---|---|
| band gap | Thurmond (Varshni form) |
| densities of states | Green 1990 |
| intrinsic density used by the device solver | measured, Misiakos & Tsamakis 1993 |
| mobility vs doping | Masetti et al. 1983 |
| velocity saturation | Canali et al. 1975 |
| Auger recombination | Dziewior & Schmid 1977 |
| SRH lifetimes | Scharfetter doping dependence |

| T (K) | 200 | 250 | 300 | 350 | 400 |
|---|---|---|---|---|---|
| n_i from bands vs measured | +10.9% | +10.5% | +10.1% | +9.8% | +9.5% |

- The band-structure n_i is 10% above the measured values, uniformly
  across temperature. This is the known gap between Green's densities of
  states and the measured n_i.
- The device solver therefore uses the measured n_i.
- Mobilities reproduce the familiar resistivity table: 1e15 cm⁻³ gives
  4.6 Ω·cm (n-type) and 13.5 Ω·cm (p-type). The Masetti fit is built on
  these data, so this is not an independent check.

## 2. Drift–diffusion solver (`device.rs`)

Model: Poisson's equation plus electron and hole continuity, with
Boltzmann statistics, SRH and Auger recombination, and optional velocity
saturation.

**Discretization.**
- Box method on tensor meshes, with a material per element. Oxide
  carries only Poisson's equation; Si/SiO₂ interface nodes have partial
  boxes.
- Scharfetter–Gummel edge currents: exact for constant current and field
  along an edge.
- Unknowns ψ, φ_n, φ_p (units of V_t), so carrier densities stay positive.
- Fully coupled Newton with damped steps and bias continuation.

**Linear solver: blocked banded LU with partial pivoting.**
- Panel factorization plus a parallel rank-32 trailing update, with an
  AVX2 kernel: 0.15 s per factorization at 32 k unknowns, 6× faster than
  the unblocked version.
- faer's supernodal sparse LU was tried and rejected. Partial pivoting
  after a fill-reducing column order lost too much accuracy on these
  Jacobians: iterative refinement diverged, and Newton failed on the MOS
  capacitor.

**Terminal currents: weighting function.**
- I = Σ_edges J·(g_a − g_b), where g solves Laplace's equation with
  g = 1 on the contact and 0 on the others.
- This is exactly equivalent to the flux through the contact when current
  is conserved.
- It also suppresses round-off: majority-carrier fluxes are tiny
  differences of drift and diffusion terms up to 10⁷ A/m². Read directly
  at the contact edge, the diode's reverse current came out 4.7% low.

**Manufactured solutions (`semiconductor mms`).** Doping up to
±3e20 m⁻³, ψ amplitude 8 V_t, SRH and Auger on.

| mesh | max error ψ / φ_n / φ_p (V_t) | observed order |
|---|---|---|
| 1D, 20 cells | 1.6e-2 / 3.0e-2 / 4.7e-3 | |
| 1D, 320 cells | 6.4e-5 / 1.2e-4 / 1.9e-5 | 2.00 / 2.00 / 1.98 |
| 2D graded (spacing ratio 4), 149×112 | 1.3e-4 / 5.3e-4 / 1.3e-4 | 2.00 / 1.98 / 1.99 |

The manufactured sources are computed by finite differences of the
continuous operator. Their step must be small: at 1e-3·L, the d²
truncation error put a floor of 1e-3 on φ_p.

## 3. pn junction (`semiconductor diode`)

Abrupt junction, 1e17 / 1e16 cm⁻³, both bases 10 diffusion lengths long
(L_n 41 µm, L_p 41 µm).

- **Equilibrium:** quasi-Fermi levels flat to 1.7e-16 V (detailed
  balance), terminal current 5e-16 A/cm².
- **Depletion capacitance** vs abrupt-junction theory with the 2V_t tail
  correction:

| V | 0 | −0.5 | −1 | −2 | −5 |
|---|---|---|---|---|---|
| error | −3.19% | −0.50% | −0.08% | −0.00% | +0.00% |

  The −3% at zero bias is the depletion approximation's own error, not
  the solver's.

- **Shockley's law, exact test.** Shockley assumes no recombination in the
  depletion region, so it is removed there (bias-dependent edges):

| V | −1.0 | +0.1 | +0.3 | +0.5 | +0.55 | +0.6 |
|---|---|---|---|---|---|---|
| J / J_Shockley | 0.99964 | 1.00031 | 1.00028 | 0.99845 | 0.98903 | 0.93087 |
| corrected for neutral-region I·R drop | 0.99964 | 1.00031 | 1.00028 | 0.99981 | 0.99827 | 0.98901 |

  The last 1% at 0.6 V is the onset of high injection (p/N_d ≈ 1%).

- **With full recombination:** ideality 1.66 at 0.2 V falling to 1.02
  near 0.5 V (depletion-region recombination, then diffusion), then
  rising again from series resistance and high injection.

## 4. MOS capacitor (`semiconductor moscap`)

n⁺ poly / 10 nm SiO₂ / p-Si 1e17 cm⁻³. The reference is the exact
Poisson–Boltzmann solution, which integrates once in closed form:
E_s = √(2qp₀V_t/ε) F(u).

| interface spacing (nm) | 0.1 | 0.05 | 0.025 |
|---|---|---|---|
| max gate-charge error / max charge | 1.6e-4 | 4.0e-5 | 1.0e-5 |

- Convergence is second order.
- The low-frequency C/C_ox agrees to 5 digits through accumulation,
  depletion and inversion, e.g. 0.24471 at the minimum, 0.97635 in
  strong inversion.

## 5. Long-channel MOSFET vs Pao–Sah (`semiconductor mosfet L`)

Device: 2D nMOSFET, 10 nm oxide, 1e17 cm⁻³ substrate, 1e20 cm⁻³
source/drain, constant mobility 500 cm²/(V s).

Reference: the Pao–Sah double integral, the exact long-channel
(gradual-channel) limit of drift–diffusion. It uses the exact
Poisson–Boltzmann charge with the channel quasi-Fermi split.

Results at V_d = 0.05 V:

| gate length | subthreshold slope (TCAD / Pao–Sah) | subthreshold I_d / Pao–Sah | strong inversion I_d / Pao–Sah | implied channel shortening ΔL |
|---|---|---|---|---|
| 10 µm | 80.4 / 80.4 mV/dec | 1.037 | 1.002–0.999 | 0.35 µm |
| 3 µm | 80.4 / 80.4 | 1.117 | 1.004–0.994 | 0.31 µm |
| 1 µm | 80.3 / 80.4 | 1.47 | 1.002–0.980 | 0.29 µm |

- **Strong inversion** (V_g ≥ 1.6 V) agrees to within 0.6% at 10 and
  3 µm, and to 2% at 1 µm.
- **Subthreshold slope** agrees to 0.1 mV/dec at every length.
- **Subthreshold current** runs higher than Pao–Sah, by a fixed effective
  channel shortening ΔL ≈ 0.3 µm (ratio = L/(L − ΔL)) at every length.
  These are the source/drain depletion regions, which the gradual-channel
  approximation ignores. The 2D solver contains this effect; Pao–Sah
  cannot.
- At V_d = 1 V the 1 µm device also departs in saturation (+20%):
  channel-length modulation, where the gradual-channel approximation
  fails at pinch-off.

This marks the boundary between the two rungs. The compact long-channel
description is valid while ΔL/L is small. The excess current is the
indicator, and it grows exactly as 1/L.

## 6. Short-channel effects (`semiconductor short`)

Scaled nMOS technology: 4 nm oxide, 5e17 cm⁻³ substrate, 60 nm
junctions, V_dd 1.2 V. Masetti mobility, plus velocity saturation driven
by the edge field. Threshold by constant current (1e-7 A·W/L); DIBL =
ΔV_t / ΔV_d between 0.05 V and 1.2 V. Nothing below is fitted: the
effects emerge from the same drift–diffusion equations as the long
channel.

| L_gate | L_met | V_t (0.05 V) | V_t (1.2 V) | DIBL (mV/V) | I_off (A/µm) | I_on (A/µm) | SS (mV/dec) |
|---|---|---|---|---|---|---|---|
| 1000 nm | 971 nm | 0.369 | 0.361 | 7 | 2.4e-12 | 6.6e-5 | 77.9 |
| 300 nm | 271 nm | 0.353 | 0.340 | 11 | 1.4e-11 | 1.9e-4 | 77.5 |
| 200 nm | 171 nm | 0.335 | 0.311 | 21 | 4.9e-11 | 2.6e-4 | 76.8 |
| 150 nm | 121 nm | 0.303 | 0.256 | 41 | 3.7e-10 | 3.3e-4 | 77.0 |
| 120 nm | 91 nm | 0.257 | 0.171 | 74 | 7.1e-9 | 4.1e-4 | 82.0 |
| 100 nm | 71 nm | 0.199 | 0.068 | 114 | 1.9e-7 | 4.8e-4 | 96.6 |
| 80 nm | 51 nm | 0.105 | −0.152 | 223 | 9.4e-6 | 6.0e-4 | 174 |
| 65 nm | 36 nm | −0.008 | −1.05 | 904 | 9.1e-5 | 7.6e-4 | 592 |

What emerges as the channel shortens:
- **Threshold roll-off and DIBL.** The drain's field reaches the source
  barrier.
- **Punch-through below about 80 nm.** The subthreshold slope degrades
  from 77 to 592 mV/dec, and the off-current rises by 7 orders of
  magnitude: the gate has lost control.

The design rule follows: this technology is usable down to about 120 nm.
Shorter gates need a thinner oxide, higher doping or shallower junctions.
That is the scaling of the natural length l = √(ε_si t_ox W_d/ε_ox),
24.2 nm here (Liu et al. 1993, with η = 1).

**Roll-off decay length.** The roll-off follows ΔV_t ∝ exp(−L/2l)
asymptotically. Local decay lengths between successive devices:

| L_met interval (nm) | 271–171 | 171–121 | 121–91 | 91–71 | 71–51 | 51–36 |
|---|---|---|---|---|---|---|
| local decay length (nm) | 64.7 | 38.1 | 28.1 | 24.1 | 22.6 | 21.1 |

- At short L these approach the natural length, 24.2 nm, within about
  10%.
- At longer L the roll-off decays more slowly, because Liu's form has a
  prefactor and a second exponential.

**Solver note.** Velocity saturation uses a frozen mobility at each
Newton iterate. The exact ∂μ/∂ψ terms, checked against finite
differences to 3e-9, made Newton diverge from equilibrium. Such terms
carry negative differential mobility, which spoils the M-matrix-like
structure of the drift–diffusion matrix. Re-evaluating the mobility at
every iterate converges in 5–7 steps to the same solution.

## 7. Folding TCAD into circuits (`semiconductor cmos`)

nMOS and pMOS from section 6, 150 nm gate, V_dd 1.2 V.

**The fold** (`compact.rs`).
- Each device is sampled by TCAD on a 13 × 16 grid in (V_gs, V_ds); this
  takes 951 s for both.
- The table model factors the current: I_d = I_sat(V_gs) · s(V_gs, V_ds),
  with monotone PCHIP in asinh(I_sat) and in the bounded shape s.

**The certificate.** The table is compared with fresh TCAD solutions at
the 180 cell centres, none of them used to build it:

| device | worst relative error (above 1e-3 I_max) | worst error below that |
|---|---|---|
| nMOS | 5.35% | 0.015 decades |
| pMOS | 5.12% | 0.015 decades |

- The first version interpolated asinh(I_d) along V_ds as well. The
  certificate caught it: 99.9% error at cell centres. asinh(I_d) jumps
  from 0 to about 25 within the first V_ds cell, a log singularity at
  V_ds = 0. That is why the current is now factored.
- A unit test against the EKV all-region law gives a 3.3% worst error at
  cell midpoints.
- Outside the sampled domain the model is not trusted. That is where the
  circuit would unfold to TCAD.

**Devices.**
- I_on: nMOS 330.5 µA/µm, pMOS 230.2 µA/µm.
- The gate capacitances are identical (1.031 fF). The pMOS is the
  electrostatic mirror image of the nMOS (complementary doping, ±E_g/2
  gates), and mobility does not enter.

**Inverter** (W_p/W_n = 1.56 for a switching threshold V_M = V_dd/2), on
the tables:
- peak gain 18.8;
- V_IL 0.512 V, V_IH 0.684 V, V_OL 0.078 V, V_OH 1.120 V;
- noise margins NM_L 0.434 V, NM_H 0.436 V.

Positive noise margins are the indicator that the Boolean gate-level
abstraction holds for this inverter.

**Check against mixed-mode TCAD.** Both drift–diffusion devices are
coupled by Kirchhoff's current law at the output (regula falsi on V_out,
5–12 TCAD evaluations per point):

| V_in | 0.2 | 0.4 | 0.5 | 0.55 | 0.6 | 0.65 | 0.7 | 0.8 | 1.0 |
|---|---|---|---|---|---|---|---|---|---|
| V_out, table (V) | 1.1999 | 1.1817 | 1.1303 | 1.0641 | 0.6000 | 0.1255 | 0.0630 | 0.0161 | 0.0001 |
| table − TCAD (mV) | 0.0 | 0.0 | −0.2 | +0.7 | 0.0 | −0.6 | +0.2 | 0.0 | 0.0 |

The table-based transfer curve agrees with direct TCAD to within 0.7 mV.
That holds through the high-gain region, where a 5% current error moves
V_out least, because both devices' errors largely cancel in the KCL
balance.

**Ring oscillator** (5 stages).
- Node capacitance 3.86 fF:
  - next-stage gates from the TCAD gate charge, averaged over 0 → V_dd;
  - drain junctions from the depletion capacitance at V_dd/2 over the
    source/drain area.
- Period 109.6 ps (9.1 GHz), delay 11.0 ps per stage. The simple estimate
  C·V_dd/(2 I_on) gives 6.7 ps; the average switching current is below
  I_on.
- Limitation: the capacitances are lumped and bias-averaged, not the
  charge-conserving terminal charges of a transient TCAD model.

## 8. Folding and unfolding: the map

Each piece of physics below is a fold/unfold pair:
- a cheap model is used wherever an indicator says it is valid;
- the heavier model switches on only where the indicator fails;
- every switch is checked against an exact or reference solution.

| physics | fold (cheap) | unfold (heavy) | indicator / switch | check |
|---|---|---|---|---|
| carrier statistics | Boltzmann | Fermi–Dirac (per node) | η > −12 (error of the fold < 2.2e-6) | degenerate MMS, 2nd order (§9) |
| band gap | rigid | Slotboom narrowing | doping | in the same MMS |
| impurity scattering | Masetti (neutral bulk) | Brooks–Herring re-screened by local carriers | carrier density ≠ doping | mobility universality (§10) |
| inversion-layer mobility | bulk | Lombardi (phonon + roughness in E⊥) | E⊥ | universal curve (§10) |
| impact ionization | off | van Overstraeten–de Man | field > 5e4 V/cm (α < 1e-5 /cm below) | breakdown vs ionization integral (§11) |
| time dependence | quasi-static tables (I, Q), 10⁴–10⁵× cheaper | mixed-mode transient TCAD, BDF2 | input edge / channel transit time ≳ 10 (error 4% at 69, 11% at 7) | inverter switching (§12) |
| minority-carrier supply | low-frequency C–V | deep depletion | ramp time / generation time | MOS ramps vs both exact limits (§12) |
| lattice temperature | isothermal device at the mean T | (non-isothermal DD: not built) | temperature spread inside the device | heat solver exact, energy balance (§13) |
| coolant | Nu = 140/17 (fully developed) | Graetz march (entrance) | x* = x/(D_h Re Pr) | Nu → 140/17, entrance asymptote (§13) |
| gate leakage | post-processing | self-consistent gate current (not built) | I_g/I_d | Fowler–Nordheim limit (§14) |
| confinement | classical | density gradient (quasi-2D, inside the device solver), then Schrödinger–Poisson | Q_cl/Q_SP − 1, DG certificate | SP reference; 2D DG = 1D DG (§15) |

The same holds inside the numerics. Where double precision cannot
resolve a quantity, the solver's convergence test is stated in terms of
what it can resolve, instead of chasing noise:
- majority-carrier currents in heavy doping;
- slow transients balanced below the round-off of near-cancelling
  fluxes.

## 9. Heavy doping: Fermi–Dirac statistics, band-gap narrowing

- Carriers are n = n_i e^{ψ+θ_n−φ_n}. The effective potential θ carries
  ΔE_g/2 (Slotboom–de Graaff) and the degeneracy correction
  ln(F_{1/2}(η)/e^η).
- Scharfetter–Gummel runs on the effective potentials, so it stays exact
  for constant current and effective field.
- The Newton matrix includes dθ/dη exactly (finite-difference Jacobian
  test: 7e-9 relative error with 97 degenerate nodes, transients
  included).
- F_{1/2} and F_{−1/2} come from Gauss–Legendre panels, checked against
  the non-degenerate series and the Sommerfeld expansion to series
  truncation (5e-8).

Manufactured solution, degenerate (`semiconductor mms`):
- potential ±23.5 V_t, η ≈ +2 on both sides;
- compensated 3e20 cm⁻³ doping, N_d + N_a constant, N_d − N_a = n − p of
  the manufactured state.

| cells | 20 | 40 | 80 | 160 | 320 |
|---|---|---|---|---|---|
| max error ψ (V_t) | 1.6e-1 | 3.7e-2 | 9.1e-3 | 2.3e-3 | 5.7e-4 |
| order | | 2.11 | 2.02 | 2.00 | 2.00 |

My first degenerate test was ill-posed. Its potential was a step
(asinh of the doping), or carried space charge that was unphysical by
orders of magnitude. Newton failed on it even for plain Boltzmann
statistics. The test now uses a neutral, smooth construction.

## 10. Inversion-layer mobility (`semiconductor mobility`)

Method: long-channel devices (2 µm, 4 nm oxide) at three substrate
dopings. μ_eff = I_d L/(W Q_inv V_d) at mid-channel, against E_eff =
(|Q_dep| + |Q_inv|/2)/ε_si.

| E_eff | Masetti + Lombardi: spread across 3e16–3e17 | + carrier screening: spread | μ_eff |
|---|---|---|---|
| 0.5 MV/cm | 35.7% | 2.5% | 454–465 cm²/Vs |
| 0.7 MV/cm | 32.1% | 1.7% | 385–391 cm²/Vs |
| 1.0 MV/cm | 28.9% | 1.1% | 311–315 cm²/Vs |

- Universality, one curve independent of doping, is how measured
  inversion-layer mobility behaves. It emerges only once impurity
  scattering is re-screened by the inversion carriers.
- Masetti is fitted to neutral silicon, where the screening density
  equals the doping. Brooks–Herring screening is unfolded where the two
  differ, and gives Masetti back exactly in neutral bulk.
- The values are close to the measured universal electron mobility,
  about 300 cm²/Vs near 1 MV/cm (Takagi et al. 1994).
- On the 150 nm device the full mobility model gives I_on 349.5 µA/µm.

## 11. Avalanche breakdown (`semiconductor breakdown`)

One-sided p⁺n junction, 1e19/1e16 cm⁻³, reverse bias, impact ionization
on vs off:

| V | −20 | −30 | −40 | −45 | −48 | −50 | −51 |
|---|---|---|---|---|---|---|---|
| multiplication M | 1.024 | 1.144 | 1.590 | 2.358 | 3.865 | 7.543 | 15.67 |
| ionization integral | 0.009 | 0.066 | 0.265 | 0.480 | 0.671 | 0.834 | 0.928 |

Breakdown voltage:

| estimate | V |
|---|---|
| drift–diffusion, extrapolating 1/M → 0 | 51.93 |
| ionization integral = 1 on the same field, with the same majority-tail fold | 51.88 |
| ionization integral, full | 51.70 |
| Sze–Gibbons empirical law | 62 |

The 16% gap to Sze–Gibbons is the ionization-coefficient set: that law
predates the van Overstraeten–de Man coefficients.

Three numerical lessons:
- **Generation is not linearized.** Exact derivatives of α|J| with
  respect to majority-carrier fluxes at depletion edges dwarf the
  minority rows they enter, and the matrix goes singular. The generation
  is instead evaluated at each Newton iterate, which converges at rate
  (M−1)/M.
- **It unfolds from the converged solution without it.** Generation
  built from early iterates' spurious currents destabilizes Newton.
- **Majority-carrier Debye tails are excluded.** In 1e19 material, the
  majority quasi-Fermi variation per edge (4e-15 V_t) is below double
  precision, so the local majority current there is noise. Over one
  Debye length this changes the breakdown voltage by 0.18 V (reported
  above), and the convergence test ignores carriers below 1e-14 of the
  doping.

## 12. Transients (`semiconductor transient`, `semiconductor switching`)

**Method.**
- Time integration: variable-step BDF2 with carrier storage, and
  displacement current through the box-method Poisson equation.
- Total current (conduction plus displacement) is conserved by the
  discrete equations, so terminal currents use the weighting function
  exactly.
- Steps that fail are subdivided.

**MOS capacitor under gate ramps** (n⁺ poly / 10 nm / 1e17): C =
I_g/(dV/dt). Two exact references: deep depletion (no inversion
electrons) and low frequency (full equilibrium).

| V_g | 1 V/µs | deep depletion (exact) | 0.1 V/s | quasi-static fold | low frequency (exact) |
|---|---|---|---|---|---|
| −1.5 | 0.90443 | 0.90436 | 0.90792 | 0.90438 | 0.90436 |
| −0.6 | 0.40489 | 0.40511 | 0.40488 | 0.40511 | 0.40511 |
| 0.2 | 0.23849 | 0.23850 | 0.24113 | 0.24226 | 0.24226 |
| 0.5 | 0.21373 | 0.21373 | 0.25603 | 0.76761 | 0.76758 |
| 1.0 | 0.18546 | 0.18546 | 0.29600 | 0.93766 | 0.93763 |
| 1.9 | 0.15429 | 0.15429 | 0.36888 | 0.97131 | 0.97126 |

- **Fast ramp:** reproduces deep depletion to 4–5 digits, with total
  current conserved to 1e-15.
- **Quasi-static fold** (steady solves): reproduces the low-frequency
  curve to 4–5 digits.
- **0.1 V/s:** the ramp takes 40 s, comparable to the ~40 s thermal
  generation time of the inversion charge, and the capacitance lies
  between the two limits. That time ratio is the indicator between the
  folds.
- **Slower ramps:** go to the quasi-static fold by necessity as well as
  by economy. At 1 mV/s the displacement current (~3.5e-6 A/m²) is
  ~1e-14 of the accumulation-layer drift and diffusion fluxes, below
  double precision. A transient run there gave noise even after 12,936
  adaptive steps.

Two numerical lessons:
- **Backward-error convergence.** In slow transients the electron
  equation balances storage terms below the round-off of
  drift/diffusion fluxes up to 1e9 times larger. Newton then stagnates
  at 1e-5–1e-3 V_t. A step is now converged when every residual is below
  1e-14 of its own row's terms. At 1e-12 the accepted residuals exceeded
  the true displacement current of a 1 mV/s ramp.
- **Currents of a committed step.** Evaluating the transient current
  after committing a step used the new level as history and read
  capacitances 5–15% low. Fixed by keeping the history the step used.

**Inverter switching** (`semiconductor switching`):
- Quasi-static table model: currents plus charge-conserving terminal
  charges from TCAD. Gate charge from Gauss's law; drain charge as
  Ward–Dutton, generalized to 2D through the conduction-weighted Laplace
  function.
- Reference: mixed-mode transient TCAD, two drift–diffusion devices plus
  the load coupled by KCL, BDF2 at the device and circuit level alike.
- Channel transit time L²/(μ_eff(V_dd − V_t)) ≈ 0.73 ps: the quasi-static
  fold needs input edges much longer than this.

Inverter with a 2 fF load, output falling:

| input edge | edge / transit | delay t_pHL: table | TCAD | fall 90→10%: table | TCAD | cost: table | TCAD |
|---|---|---|---|---|---|---|---|
| 50 ps | 69 | 21.23 ps | 22.14 ps (−4.1%) | 19.87 ps | 20.75 ps | 0.04 s | 2815 s* |
| 5 ps | 6.9 | 13.41 ps | 15.04 ps (−10.8%) | 14.58 ps | 16.02 ps | 0.03 s | 2447 s |
| 1 ps | 1.4 | 13.40 ps | 14.81 ps (−9.6%) | 14.58 ps | 16.02 ps | 0.05 s | 3042 s |

\* With the cores shared with another run; about half that alone.

- The quasi-static fold is 4% fast at slow edges, consistent with its 5%
  current certificate. It is 10⁴–10⁵ times cheaper than mixed-mode TCAD.
- **Miller overshoot** (V_out above V_dd as the input rises):
  - 50 ps edge: table 1.253 V, TCAD 1.246 V;
  - 5 ps edge: table 1.663 V, TCAD 1.430 V;
  - 1 ps edge: table 2.091 V, TCAD 1.490 V.
- At fast edges the quasi-static model assumes the channel charge
  responds instantly. It therefore overstates the gate-to-drain coupling
  and the switching speed.
- The fold's validity boundary emerges where expected:
  - at 7 transit times per edge the delay error (11%) has outgrown the 5%
    certificate;
  - at 1.4 the overshoot error reaches 0.6 V, while the delay error stays
    near 10%, because the waveform is distorted rather than shifted.
- Below an edge/transit ratio of about 10, the circuit should unfold to
  transient TCAD, or to a non-quasi-static compact model (charges
  relaxing over the transit time).

## 13. Temperature (`semiconductor thermal`)

**Heat solver** (2D box method, k_Si = 148(300/T)^1.3 by Picard,
oxide 1.4 W/(m K)). Reference: a silicon slab with uniform heating, exact
through the Kirchhoff transform even with nonlinear k:

| cells | 10 | 20 | 40 | 80 |
|---|---|---|---|---|
| max error (K) | 1.2e-2 | 3.0e-3 | 7.6e-4 | 1.9e-4 |
| order | | 2.00 | 2.00 | 2.00 |

**Energy bookkeeping.** H = −J_n·∇φ_n − J_p·∇φ_p + q(R−G)(φ_p−φ_n). On
the 150 nm device at V_g = V_d = 1.2 V, Σ H equals I_d·V_d to 6e-14.

**Temperature-dependent physics.** Diode saturation current against
Shockley with n_i(T), μ(T), V_t(T):

| T (K) | 250 | 300 | 350 | 400 |
|---|---|---|---|---|
| J₀ / Shockley | 1.00031 | 1.00031 | 1.00014 | 0.99606 |

That is over a range in which n_i² changes by 10 orders of magnitude.

**Leakage of the 150 nm nMOS** (V_g = 0, V_d = 1.2 V):

| T (K) | 300 | 325 | 350 | 375 | 400 |
|---|---|---|---|---|---|
| I_off (A/µm) | 3.1e-10 | 1.0e-9 | 2.9e-9 | 7.4e-9 | 1.7e-8 |

It doubles every 17.4 K.

**Self-heating** (electro-thermal loop, bottom of the 0.5 µm domain at
the die temperature):
- converges in 2 iterations: mean device temperature +3.3 K, I_on 349.5 →
  347.2 µA/µm;
- temperature spread inside the device 2.5 K, so the isothermal-device
  fold holds; it would unfold to non-isothermal drift–diffusion near
  ~10 K.

**Coolant.** Laminar flow between parallel plates, both walls at constant
heat flux (the Graetz problem). The energy equation is marched in the
flow direction with axial conduction folded away (Pe ≫ 1), BDF2.

| x* = x/(D_h Re Pr) | 1e-5 | 1e-4 | 1e-3 | 1e-2 | ≥ 5e-2 |
|---|---|---|---|---|---|
| Nu_x | 68.99 | 32.15 | 15.43 | 8.80 | 8.2353 |
| entrance asymptote 1.490 x*^−1/3 | 69.14 | 32.09 | 14.90 | — | — |

- Fully developed: 8.23529, the exact 140/17.
- Downstream of x* ≈ 0.05 the fold Nu = 140/17 is exact to 4 digits.
- Crank–Nicolson rang around the limit (8.19–8.28) once the steps grew.
  BDF2 is L-stable and doesn't.

**Chip and thermal runaway.**
- Chip: 1 cm², 2e9 µm of off transistors leaking I_off(T)·V_dd, dynamic
  power 40 W, coolant at 300 K.
- Fixed points of T = T_c + R_th(P_dyn + P_leak(T)):
  - microchannels (h = 140/17·k/D_h, R_th 0.455 K/W): 318.9 K stable,
    449.7 K unstable;
  - air sink (0.6 K/W): 325.2 / 439.5 K;
  - poor sink (1.2 K/W): 356.5 / 403.7 K.
- With worse cooling the pair merges and disappears: runaway. The
  unstable point is the margin; a transient excursion past it never
  returns.

## 14. Gate tunnelling (`semiconductor tunnel`)

Tsu–Esaki current with closed-form WKB transmission through the actual
trapezoidal oxide barrier (Φ_B 3.15 eV, m_ox 0.5 m₀), as post-processing
of the device solution.

**Fowler–Nordheim limit** (metal emitter, 10 nm oxide):

| F (MV/cm) | 6 | 8 | 10 | 12 | 14 |
|---|---|---|---|---|---|
| J / J_FN | 0.979 | 0.972 | 0.966 | 0.959 | 0.953 |

The ratio falls linearly with F and extrapolates to 1.000 at F → 0. That
is the size of the barrier linearization in the Fowler–Nordheim
derivation itself.

**Direct tunnelling** in the inversion layer (n⁺ poly / p 1e18, V_g 1 V):

| t_ox (nm) | 3.0 | 2.0 | 1.6 | 1.4 | 1.2 | 1.0 |
|---|---|---|---|---|---|---|
| J_g (A/cm²) | 6.1e-9 | 1.3e-3 | 0.19 | 2.3 | 27 | 340 |

- 1.07–1.09 decades per 0.2 nm, matching the well-known one decade per
  2 Å of SiO₂.
- The 4 nm oxide of the scaled technology leaks negligibly. Below
  ~1.5 nm the fold indicator I_g/I_d grows past 1e-3 and the gate current
  would have to enter the device equations.

## 15. Quantum confinement (`semiconductor quantum`)

Three models of the inversion layer (n⁺ poly / SiO₂ / p-Si), 1D in
depth:
- **Schrödinger–Poisson:** Δ2/Δ4 valley subbands, 2D Fermi–Dirac
  occupation, self-consistent. Its eigen-solver reproduces the infinite
  square well within 0.2% on the graded mesh.
- **Density gradient:** one mass parameter, calibrated to Schrödinger–
  Poisson on stack A at V_g = 1 V, giving m_dg = 0.201 m₀.
- **Classical:** the drift–diffusion limit.

| stack | V_g | Q_inv classical | Q_inv SP | DG error | centroid classical / SP / DG |
|---|---|---|---|---|---|
| A: 1e18, 2 nm (calibration) | 0.6 | 4.08e-3 C/m² | 2.31e-3 | −1.3% | 0.40 / 1.53 / 1.73 nm |
| | 1.8 | 2.36e-2 | 1.90e-2 | +0.9% | 0.25 / 1.23 / 1.35 nm |
| B: 3e17, 3 nm | 0.6 | 3.67e-3 | 2.76e-3 | +0.6% | 0.69 / 2.02 / 2.14 nm |
| | 1.8 | 1.68e-2 | 1.45e-2 | +0.4% | 0.42 / 1.53 / 1.62 nm |
| C: 3e18, 1.5 nm | 0.6 | 2.90e-3 | 5.04e-4 | −22.8% | 0.25 / 1.22 / 1.45 nm |
| | 1.0 | 1.10e-2 | 6.45e-3 | −2.4% | 0.21 / 1.13 / 1.34 nm |
| | 1.8 | 2.85e-2 | 2.09e-2 | +1.3% | 0.16 / 1.02 / 1.18 nm |

| stack | threshold shift SP − classical | DG − SP | effective oxide penalty | fold indicator at 1.2 V |
|---|---|---|---|---|
| A | +115 mV | +2.7 mV | +0.31 nm | 31% |
| B | +77 mV | −1.9 mV | +0.32 nm | 19% |
| C | +173 mV | +14 mV | +0.31 nm | 54% |

- The centroid moves from ~0.2–0.7 nm (classical) to 1–2 nm. The
  inversion capacitance acts like ~0.3 nm of extra oxide.
- For these stacks classical drift–diffusion is invalid: its indicator is
  19–54% in charge.
- Density gradient certifies within 2.4% above threshold on stacks it
  was not calibrated on. It fails near threshold at 3e18 / 1.5 nm
  (−23%), and there the certificate says to unfold to
  Schrödinger–Poisson.
**Density gradient inside the 2D solver** (`semiconductor quantum2d`).
Quasi-2D construction:
- in every mesh column under oxide, solve the 1D density-gradient
  equation for √n given the classical exponent ψ + θ_n − φ_n;
- the quantum potential Λ = 2 ln √n − (ψ + θ_n − φ_n) enters θ_n;
- an outer iteration alternates the 2D solve and the column solves,
  damped, and backs off when an inner solve fails.

This folds full 2D density gradient for inversion layers, whose
confinement is along y. The indicator is channel length against
inversion thickness.

2D MOS capacitor (1e18 / 2 nm), inversion charge against the 1D
solutions:

| V_g | 0.6 | 1.0 | 1.4 | 1.8 |
|---|---|---|---|---|
| 2D classical vs 1D classical | +0.00% | −0.00% | −0.00% | +0.00% |
| 2D DG vs 1D DG | +0.00% | +0.00% | +0.00% | +0.00% |
| 2D DG vs Schrödinger–Poisson | −1.34% | +0.00% | +0.53% | +0.88% |

**The 150 nm nMOS (4 nm oxide, 5e17) with the quantum potential:**

| model | V_t (constant current, V_d 0.05 V) | I_on (V_g = V_d = 1.2 V) |
|---|---|---|
| classical | 0.3122 V | 349.5 µA/µm |
| density gradient (2D) | 0.3841 V | 281.0 µA/µm |

- The threshold rises by 72 mV and the on-current falls by 20%.
- The 1D references for this stack are threshold shifts of +96.2 mV
  (Schrödinger–Poisson) and +95.4 mV (density gradient), at a charge
  criterion of 1e12 cm⁻². The 0.8 mV gap certifies the fold on this
  stack as well.
- The 2D shift is smaller. It is defined by constant current, at lower
  inversion, with 2D electrostatics included.
- Confinement is the largest single correction to the classical on-
  current in this technology. Classical drift–diffusion is not adequate
  for oxides of a few nanometres, which the fold indicator says
  directly.

## 16. Still open

- Non-isothermal drift–diffusion: needed when the temperature spread
  inside a device reaches ~10 K.
- Self-consistent gate current: needed when I_g/I_d > 1e-3, oxides below
  ~1.5 nm.
- Full 2D density gradient (lateral confinement, source–drain
  tunnelling) for channels below ~20 nm.
- Ballistic and quasi-ballistic transport below ~20 nm channels, beyond
  drift–diffusion (the next rung down).
- The gate and architecture rungs: timing from these tables, then logic.
