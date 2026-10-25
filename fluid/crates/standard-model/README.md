# Standard Model kinetic module

CPU transport for the complete elementary species set: six quark flavors and
six lepton flavors, their antiparticles, photon, eight adjoint gluon components,
W+, W−, Z and Higgs. There are 37 catalog entries. Quark colors are averaged,
with their factor of three included in thermal degeneracies. Spin, mass,
charge, baryon/lepton numbers, PDG identifiers and available widths are exposed.

This implements **kinetic closures** within the simulator's emergence ladder.
It does not solve quantum fields, confinement, parton showers, electroweak
production or arbitrary Standard Model scattering. In particular, listing a
species does not supply every interaction it participates in.

## Run

From `fluid/`:

```sh
RAYON_NUM_THREADS=8 cargo run -p standard-model --release --example particles -- 1000 partonic
RAYON_NUM_THREADS=8 cargo run -p standard-model --release --example particles -- 1000 vacuum
cargo run -p standard-model --release --example relativistic
cargo test -p standard-model --release
```

The demo injects copies of every allowed species, evolves decay cascades for
100 microseconds, checks total energy/momentum and charges, and prints a
Compton encounter, photon-packet conversion and ideal-gas thermal states.

## Particle transport

```rust
use standard_model::{Domain, Particle, Simulation, Species, Vec3};

fn main() -> Result<(), &'static str> {
let mut sim = Simulation::new(Domain::Vacuum);
sim.insert(Particle::new(
    Species::Muon,
    Vec3([0.0, 0.0, 0.0]), // meters
    Vec3([0.5, 0.0, 0.0]), // GeV/c, numerically GeV with c=1
    42,                    // reproducible per-history seed
    Domain::Vacuum,
)?)?;
sim.advance(10e-6)?;        // seconds; proper-time lifetime with Lorentz dilation
let (four, charges) = sim.totals();
Ok(())
}
```

`advance` moves particles at c p/E, processes decay at its sampled event time,
and moves each daughter only for the remainder of the interval. No fluid grid
CFL or substeps at the shortest boson lifetime are required. Histories are O(n),
parallel for at least 1024 particles, and reduce in input order. Ordinary
ballistic histories avoid individual queues and heap allocations. The finite
cascade depth is bounded by the implemented channel graph. Thread count does
not change a history's RNG stream. A rejected step leaves the simulation intact.

`advance_uniform_field(dt, UniformField { electric, magnetic })` integrates the
covariant Lorentz force and spacetime worldline together in prescribed SI fields.
Pure electric/magnetic fields use analytic solutions; mixed fields use the
constant four-generator's matrix exponential and its time integral. A monotone
proper-time solve finds the laboratory endpoint or the sampled decay event.
Daughters follow their own field orbits for the remaining interval. Magnetic
motion needs no gyro-period substeps, even across many turns. Time-varying or
spatially varying fields require an external field sampling/integration closure.

The returned `FourMomentum` is the energy and momentum supplied by the field:
`totals_after = totals_before + field_exchange`. Charged histories are not
closed systems. `proper_age()` exposes elapsed proper time. Known-mass
light-cone boosts and pair invariant calculations retain mass information that
subtracting E²−p² loses for nearly parallel, highly boosted particles. Use
`pair_mass`/`dot_on_shell` with the actual species masses for such calculations;
the generic `mass2()` cannot recover information already rounded away in E.

`Particle::push_uniform(E, B, dt)` remains available as a Boris momentum source
for externally managed drift/force splitting. It returns field work in GeV.
`set_momentum` preserves the sampled proper decay clock, including after recoil.
`particles_mut()` exposes these operations for live particles.

Prescribed-field trajectories omit self-consistent Maxwell evolution, radiation
reaction, synchrotron/bremsstrahlung emission, spin precession and strong-field
QED. `UniformField::quantum_chi` reports a rest-field nonlinearity indicator;
small chi alone does not establish negligible radiative losses. The EM solver
uses special relativity in flat spacetime.

`advance_gravity(dt, metric, geodesic::Options::default())` transports all live
species in a prescribed static isotropic [gravitational metric](../gravity/README.md),
including null photon/neutrino paths and local proper-time decays. Particle
momenta are measured in local orthonormal frames, while positions and simulation
time use the metric coordinates. Canonical momentum is `beta*p_local`. Adaptive
RK4 with step doubling integrates the Hamiltonian trajectory and decay clock;
daughters follow the remaining coordinate interval. Rejected updates are atomic.
`Options` controls relative error, absolute position error, maximum coordinate
step and maximum integration attempts. This integrator is not symplectic;
monitor tolerances and energy drift for long orbits.

`killing_energy(metric)` sums `alpha(x)*E_local`, the stationary conserved energy,
including the unresolved event ledger at recorded positions. The return from
`advance_gravity` is the change in the sum of local four-momenta for bookkeeping,
not a globally conserved GR four-vector: vectors at different positions inhabit
different frames. Geometry is frozen; no source backreaction or simultaneous
EM-plus-gravity force integrator is supplied. Schwarzschild is vacuum exterior
only, and the weak-field closure is leading order, with an explicit domain guard.

`Simulation::compton(photon_index, lepton_index)` applies one scheduled
Klein–Nishina encounter and retains the target's proper decay clock. The
unpolarized free-lepton kernel accepts moving electrons/positrons, muons and
taus, with recoil and Lorentz boosts. Its cross section is in m².
Forward scattering preserves the sampled energy fraction and 1−cos(theta)
separately, retaining recoil when cos(theta) rounds to one. Numerical stress
tests at extreme energies check this identity within the tree-level model;
they do not establish the adequacy of tree-level physics at those energies.
`annihilate(electron_index, positron_index)` replaces any supported moving pair
with two photons using the finite-mass tree-level differential QED kernel.
`annihilate_at_threshold` remains available for explicitly comoving pairs.
`pair_produce(first_photon_index, second_photon_index, flavor)` implements the
reverse Breit–Wheeler process for `Electron`, `Muon` or `Tau` and its conjugate,
with invariant thresholds, finite-mass angular weights and new proper lifetimes.
Total pair/annihilation cross sections are in m². The annihilation event cross
section includes the identical-final-photon factor 1/2 over the full sphere:
sigma_pair = 2 beta_CM² sigma_ann for the electron channel. At threshold sigma_ann
diverges, but `annihilation_sigma_v` gives the finite CM Møller flux product in
m³/s; its limit is pi r_e² c. The Møller flux speed 2 beta c is not a particle speed.
The pair sampler rejects CM collinear peaks with 4 m²/s < 1e-14 rather than
fabricating angles below double-precision resolution.

`scatter_neutrino(neutrino_index, electron_index, couplings)` adds elastic
neutrino/antineutrino–electron recoil for all three flavors. Use
`weak::Couplings::tree(species, sin2_theta_W)` for explicit tree-level chiral
coefficients (electron flavor includes charged current, antiparticles swap
left/right), or supply calibrated effective coefficients. The angle and its
renormalization convention are caller inputs. The differential spectrum retains
the electron mass interference term; total cross sections are in m² and
differential recoil cross sections in m²/GeV. The low-transfer Fermi closure
rejects 2 m_e E_nu / M_W² > 1e-3 in the electron rest frame. It has no radiative
corrections, nuclear targets, inelastic weak channels or absorption.

These methods apply **caller-scheduled encounters**; the caller must supply
collision rates, spatial eligibility and event scheduling. Indices can change
after decays or annihilation. `qed::elastic` supplies two-body kinematics for an
external angular/rate closure; it does not supply a strong interaction law.

## Decays and retained products

The inclusive channel graph covers muon, tau, top, W±, Z and Higgs:

- Muon and leptonic tau decays use unpolarized tree-level V−A phase space with
  finite charged-daughter mass and massless neutrinos. The electron Michel
  spectrum and both neutrino energy marginals are checked independently.
  Radiative corrections and polarization/spin correlations are omitted.
- Top uses the dominant bW channel. Rare sW/dW modes are omitted.
- W has leptonic and CKM-weighted two-parton modes; Z has charged-lepton,
  neutrino and five-flavor quark pairs. Inclusive PDG weights are normalized
  together to resolve rounding/inconsistency of independently measured values.
  Approximate CKM weights are not a precision electroweak calculation.
- Higgs has bb, cc, tau pairs, muon pairs, gamma pairs, gluon pairs and Z gamma.
  WW* and ZZ* remain explicit unresolved events; two on-shell W/Z bosons are
  forbidden at its mass. Small unlisted modes are omitted. The reference
  branching weights are at 125.09 GeV, used as a near-mass closure at 125.20 GeV.
- Hadronic tau modes remain unresolved in both domains, including their
  neutrino energy: they require a hadronic spectral-function decay closure.

`Domain::Vacuum` rejects independent light quarks and gluons. Colored daughters
become **hadronization records**, each retaining four-momentum, additive charges,
position, time and origin. Process these partons jointly with a QCD closure.
`Domain::Partonic` transports them as color-averaged, noninteracting on-shell
partons. Light quark masses are running mass parameters used as kinetic model
masses. No isolated light-quark lifetime is fabricated: these species are
transported without a weak/QCD decay closure, rather than declared stable in
nature. Gluon labels are adjoint basis components, not separately observed
particle species; all share PDG ID 21.

`unresolved()` shows the retained event ledger; `totals()` includes its energy,
momentum and charges. No unresolved probability is renormalized away or
silently turned into fluid heat. `take_unresolved()` transfers ownership to an
external closure; that caller then owns the corresponding conservation ledger.
These records do not themselves propagate or thermalize.

Neutrinos use massless flavor transport plus the elastic electron kernel above.
Finite masses, oscillations and weak absorption need additional kernels.

## Thermal and fluid/photon adapters

`thermal::species(s, kBT_GeV)` and `thermal::equilibrium(kBT_GeV, domain)` give
zero-chemical-potential ideal Fermi–Dirac/Bose–Einstein number density, pressure
and energy. Massive species use momentum quadrature; massless ones use analytic
integrals. SI conversion methods give m⁻³, Pa and J/m³. All antiparticles and
color degeneracies are counted once, with a formal relativistic limit
**g*=106.75**. Vacuum thermal sums omit all colored species and all hadrons.
Neither domain is an interacting EOS for the QCD/electroweak transition, dense
matter, a chemically frozen gas, or a fluid-cell matter EOS replacement.

`coupling::photon_packet(&particle, seed)` maps one physical photon to the
existing SI radiation packet (photon_count=1). Rendering still requires a
scene and optical data; high-energy photons are outside visible camera bands.
Do not interpret a packet containing many optical photons as one energetic
particle in a QED collision.

`coupling::deposit(cell, transferred_four)` converts GeV energy/momentum to an
SI `radiation::transport::Deposit`. The compressible crate's
`standard_model::deposit(&mut solver, transfers)` applies these sources using
the existing atomic admissibility check. The interaction closure must decide
what transfers: photons/neutrinos/hadrons are not automatically absorbed.
Transferred energy is **total** energy; decide explicitly whether rest energy
belongs in the fluid source. Existing fluid mass/species fields do not change.

## Data provenance

Values are pinned for reproducibility, rather than fetched during simulation.

- [PDG 2025 mass/width/ID data](https://pdg.lbl.gov/2025/mcdata/mass_width_2025.txt),
  generated April 22, 2025: the central particle masses and widths.
- [PDG 2025 gauge/Higgs summary](https://pdg.lbl.gov/2025/tables/rpp2025-sum-gauge-higgs-bosons.pdf):
  W and Z inclusive branching weights. Experimental Higgs fractions are not
  substituted for theoretical SM branching weights.
- [PDG 2025 lepton summary](https://pdg.lbl.gov/2025/tables/rpp2025-sum-leptons.pdf):
  tau leptonic branching weights.
- [PDG muon decay review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-muon-decay-params.pdf):
  unpolarized Michel spectrum for validation of the V−A kernel.
- [CERN Report 3 Higgs branching table at 125.09 GeV](https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAtMH12509_2014):
  pinned approximate Higgs branching weights.
- [PDG passage of particles through matter](https://pdg.lbl.gov/2025/reviews/rpp2024-rev-passage-particles-matter.pdf):
  radiation interaction context. The Compton total kernel is independently
  tested against numerical integration of its differential kernel.
- [Chin, relativistic motion in constant electromagnetic fields](https://arxiv.org/html/0809.0859):
  proper-time field-generator solution and its Lorentz-covariant formulation.
- [Geant4 two-photon annihilation reference](https://geant4.web.cern.ch/documentation/pipelines/master/prm_html/PhysicsReferenceManual/electromagnetic/electron_incident/annil.html):
  independent Heitler total cross section for the moving annihilation test.
- [WarpX linear Breit–Wheeler reference implementation](https://warpx.readthedocs.io/en/26.09/_static/doxyhtml/_linear_breit_wheeler_cross_section_8_h_source.html):
  total cross section and angular kernel; our finite-mass QED law is also checked
  by angular integration, detailed balance and sampled moments.
- [Tomalak and Hill, elastic neutrino–electron scattering](https://arxiv.org/html/1907.03379):
  finite-recoil chiral Fermi kernel and G_F. Our tree-level closure does not
  implement the paper's radiative corrections or its calibrated precision rates.

Mass/width uncertainties, in-medium shifts, running couplings and scheme
conversions are not propagated. Catalog values alone are not precision QCD data.
