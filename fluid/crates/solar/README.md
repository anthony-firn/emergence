# Solar background and magnetic transport

This crate supplies a published stellar background, surface magnetic transport,
and a potential coronal extrapolation. It does not yet evolve a self-consistent
Sun or provide validated solar-activity forecasts.

```sh
cd fluid
RAYON_NUM_THREADS=8 cargo run -p solar --release --example sun
RAYON_NUM_THREADS=8 cargo run -p solar --release --example hindcast
RAYON_NUM_THREADS=8 cargo run -p solar --release --example render -- /tmp/solar-render 28
# Explicit camera orbit for inspection, separately from simulated evolution:
RAYON_NUM_THREADS=8 cargo run -p solar --release --example render -- /tmp/solar-orbit 28 orbit
# Optional first argument writes the forecast grid in tesla:
cargo run -p solar --release --example hindcast -- /tmp/solar_forecast.csv
cargo test -p solar --release
python3 scripts/fetch_solar_reference.py  # reproduce bundled GONG maps
```

## Stellar structure and particle gravity

`structure::Profile::model_s()` loads 2482 published knots of sound speed,
density, pressure, Gamma1 and temperature, converts cgs to SI, and retains the
original Model S radius and mass convention. `gravity(center)` uses the exact
volume-average density of each linearly interpolated shell to build the
regular interior [gravity field](../gravity/README.md). It uses modern pinned G;
the historical model's G is separately exposed. `at(r)` interpolates the data,
without an EOS inversion or nuclear evolution.

The source is [Christensen-Dalsgaard's Model S archive](https://users-phys.au.dk/jcd/solar_models/),
specifically its [basic six-column table](https://users-phys.au.dk/jcd/solar_models/cptrho.l5bi.d.15c).
The archive notes a corrected Ledoux discriminant in full FGONG files. We use
the basic thermodynamic table and do not read that discriminant. SHA-256:
`65ecb920ed81b6b41f733cb8ab6f8c30941f7c743b0b6fec831de30e9a7322cc`.

Nominal radius, GM and luminosity are the exact conversion constants from
[IAU Resolution B3](https://arxiv.org/html/1510.07674), separately from the
historical Model S parameters. The `sun` example demonstrates exterior vacuum
photon redshift with Schwarzschild geodesics. Matter opacity is not applied in
that example. Gravity can act on every transported elementary species, but
solar ions are composite matter represented by bulk density, not free quarks.

## Surface flux and polar fields

On equal-area cells `mu=sin(latitude)` and periodic longitude `phi`, `Surface`
evolves radial field B in tesla:

```text
dB/dt = -Omega(mu) dB/dphi - d(u_mu B)/dmu
        + eta/R² [d((1-mu²)dB/dmu)/dmu + d²B/dphi²/(1-mu²)] + S
u_mu = 2 v_peak/R mu(1-mu²)
```

This is the classical [surface flux transport equation](https://arxiv.org/html/2303.01209).
Defaults are poleward peak flow 15 m/s, effective diffusivity 250 km²/s, and
the Snodgrass synodic rotation profile minus Carrington's 13.2 degrees/day.
These are prescribed model choices, not an interior dynamo solution. Diffusivity
models unresolved supergranular mixing, not microscopic resistivity.

Fourier phase shifts integrate differential rotation exactly for resolved
longitude modes. Latitude fluxes are conservative second-order linear MUSCL
with RK2; diffusion uses a second-order [L-stable SDIRK method](https://mooseframework.inl.gov/source/timeintegrators/LStableDirk2.html). Its two stages
reuse a tridiagonal factorization, cached across steps and shared by positive
and negative Fourier modes. Symmetric splitting is second-order in time. Polar face fluxes vanish. The meridional
CFL and configurable maximum step control integration; refine time and space
for accuracy. Unlimited slopes can ring on unresolved sharp
features. Even-grid longitude Nyquist power is explicitly removed on input.

`add_bipole` supplies equal and opposite spherical magnetic lobes using specified
location, flux, size, separation and tilt. It rejects severely underresolved
sources. Emergence must come from observations or an explicitly uncertain source
model; no artificial cycle or flare events are injected automatically.
`balance_flux` explicitly removes and reports a magnetogram's uniform bias.
`net_flux`, `unsigned_flux`, `axial_dipole` and `mean_field` are diagnostics.

`axisymmetric()` folds the full model into its longitude mean. For the same
prescribed linear coefficients and source averages, this reproduces the full
discrete mean and axial dipole to roundoff with much less work. It discards
local active-region information. Nonlinear or longitude-dependent flow feedback
would invalidate this exact fold.

## Coronal field

`corona::Pfss::from_surface` projects explicitly flux-balanced cell averages
onto resolved orthonormal spherical harmonics. Latitude integrals use Gaussian
quadrature in polar angle, avoiding polar endpoint singularities; longitude
cell integrals use an analytic sinc factor. Averaged observations are not
silently treated as point samples. The potential-field source-surface solution has
analytic gradients and regular polar limits, including nonaxisymmetric modes.
The inner radial field is a finite harmonic projection; the outer source
surface is radial. Source radius and harmonic degree are caller choices.
The boundary representation and harmonic truncation still require resolution checks.
`magnetic_field` returns Cartesian tesla and `open_flux` integrates unsigned
outer flux. [Mackay and Yeates](https://link.springer.com/article/10.12942/lrsp-2012-6)
describe the potential-field approximation and its physical limits.

`corona::energy` compares fields to a supplied potential reference over SI
volumes. Boundaries and integration domain must agree. Negative excess is
reported, not silently clamped. PFSS itself contains no electric currents,
free-energy reservoir, heating or reconnection, and cannot predict flares/CMEs.

## Physical rendering and magnetic geometry

The `render` example writes 900×900 PNG frames, a local HTML viewer with a time
slider, per-frame diagnostic CSVs, raw field maps and SI trace vertices. A
28-frame sequence advances the initial GONG map over 27 days without new
emergence. The default camera stays fixed in Carrington coordinates; `orbit`
adds a labeled inspection turntable and does not represent solar rotation.

`photosphere` uses the Eddington grey LTE relation
`T^4=3/4 Teff^4(tau+2/3)`. Instead of assigning a constant surface blackbody,
it integrates the formal outward Planck transfer with 32-point Gauss-Laguerre
quadrature. Visible radiance is integrated through the existing CIE response.
Nominal solar luminosity/radius set Teff≈5772 K. The analytic bolometric law
integrates to `sigma*Teff^4`, with a limb/center ratio of 0.4. This is an
approximate continuum atmosphere, not a line-blanketed measured solar spectrum.
[MESA's atmosphere documentation](https://docs.mesastar.org/en/24.03.1/atm/t-tau.html)
states the same grey temperature closure.

`trace::trace` integrates `dx/ds=±B/|B|` in both directions using adaptive RK4
step doubling. Coordinates and exported vertices are SI. The default boundary
and local position-error tolerance is 1e-6 solar radii. Only two photospheric endpoints
classify a closed curve; a photospheric/source-surface pair classifies an open
curve. Near-null, length, step and underflow limits remain explicit. Curves are
spatial field geometry, not particle trajectories or plasma time evolution.

The renderer hides rear geometry behind the opaque photosphere. Red/blue disk
colors encode actual signed field values on an explicit asinh scale; gold and
green lines mark closed/open curves. These overlays have no radiometric meaning.
Seed counts are not magnetic flux fractions. The separate continuum PFM retains
linear RGB radiance; PNG exposure, Reinhard and sRGB are display transforms.
Angle interpolation and pixel quadrature have discretization error; no Monte
Carlo standard error is claimed for deterministic rendering.

There are no fabricated spots, granulation, coronal emission or heating. The
render uses straight orthographic rays, independently of the particle geodesic
solver; GR light propagation corrections are omitted. The coronal boundary is
a modeling choice, not a self-consistent solar-wind solution. The
[preview, movies and verification report](../../../analysis/solar_render.md)
document the generated artifacts and limitations.

## Observational comparison and remaining targets

The bundled GONG radial synoptic maps are dated 2019-01-08 and 2019-02-04.
[The provenance manifest](data/gong_provenance.json) records source URLs, hashes,
UTC dates and preparation. The stdlib-only fetcher checks FITS metadata,
reindexes absolute Carrington longitude, averages 3x3 equal-area cells to
60x120, and converts Gauss to tesla. The demo applies identical Nyquist filtering
and explicit mean removal to both maps, reporting those removed means.

`hindcast` fixes defaults before evaluation, advances only the initial map,
then loads the target to report area-weighted RMSE against persistence and
`1-MSE_forecast/MSE_persistence`. It assumes no new flux emergence; it never
fits parameters or sources to the target. [NSO's product documentation](https://magmap.nso.edu/)
explains that synoptic maps combine observations over a rotation, use a radial
field approximation and include fitted polar fields. This single dated test is
a reproducible smoke comparison, not a calibrated multi-case forecast claim.

The requested prediction targets have the following implementation status:

| Target | Available foundation | Required work before predictive claims |
|---|---|---|
| Cycle strength/timing | Efficient polar-field and axial-dipole evolution | Emergence ensembles, cycle calibration, interior dynamo, held-out cycles and uncertainty |
| Magnetic-region evolution | Full surface transport, explicit bipoles, GONG comparison | Observational source assimilation, measured flows, multi-case spatial/time convergence and forecast scores |
| Flares/CMEs | Potential coronal topology, open flux, energy comparison API | Vector boundary data, evolving nonpotential induction, energy/helicity injection, reconnection, eruption physics and held-out event validation |

Full stellar evolution additionally needs a composition-dependent EOS, nuclear
reaction network, radiative/convective energy transport and a hydrostatically
balanced self-gravity update. The existing low magnetic Reynolds number imposed
field module is not solar dynamo MHD. Simulating individual elementary particles
across the entire Sun is not the reduction used here: bulk gravity and resolved
magnetic transport retain the relevant macroscopic variables efficiently.
