# Radiation and simulation cameras

Shared SI optical data power two different estimators: energy packets for
matter coupling, and radiance paths for cameras. A packet represents many
incoherent photons, with count Eλ/(hc). Every photon is not stored separately.

- `transport::fly_batch`: parallel static gas/diffuse-volume transport,
  exact piecewise-constant voxel traversal, implicit capture, HG scattering,
  energy and momentum deposition. Finite flight times persist live packets.
  Event-budget exhaustion returns pending time that MUST be resumed.
- `transport::emit_band`: isotropic Kirchhoff emission and opposite matter
  recoil in a narrow wavelength bin. Apply its signed withdrawals before
  flight. Returns an error if a step exhausts the supplied thermal budget.
- `thermal::exchange`: positive, conservative implicit grey LTE exchange;
  no absorption-time CFL. Uses n=1 radiation, a supplied Planck-mean opacity
  and fixed Cv. The fluid's `System::phases` must supply a valid fixed-volume
  heat capacity rather than a generic approximation.
- `thermal::diffuse_closed`: conservative O(n) implicit 1D diffusion with
  harmonic face diffusivity. Closed boundaries; requires transport optical
  depth ≥10 in every cell. This is a thick-medium operator, not a 3D hybrid.
- `thermal::Indicators`: chooses an eligible rung using optical depth,
  gradients, reduced flux, directional structure, equilibration and feature
  size. It does not automatically transfer fields between rungs.
- `render::Scene`: CPU spectral path tracing, BVH surfaces, exact Fresnel
  reflection/refraction, polished conductors, volume absorption/emission/HG
  scattering, area-light sampling, Russian roulette and deterministic
  parallel pixels. Four correlated wavelength paths form each sample bundle to reduce color
  noise, with independent paths when dispersion changes them. The sample
  count is wavelength paths/pixel (thus 256 gives 64 bundles). Stratified
  wavelengths use a full-support CIE importance PDF over 380–780 nm; CIE 1931 is the
  Wyman–Sloan–Shirley analytic fit. PNG is exposed/Reinhard/sRGB; PFM retains
  unexposed floating-point RGB. Per-pixel luminance standard errors and
  truncated-path counts are retained in `Image`.

`compressible::radiation` reads arbitrary solver snapshots, applies signed
packet energy/momentum sources atomically, and exchanges a grey radiation
field with thermally equilibrated cells. Radiation is an explicit split
operator, not silently enabled for existing simulations. EOS fits cannot
determine optical properties: supply measured data for each material/phase.

```sh
cd fluid
cargo test -p radiation --release
cargo test -p compressible --test radiation --release
cargo run -p compressible --release --example photon_render -- renders 40 256 640 12
```

The output path is relative to the current directory. Demo arguments are output
directory, fluid grid size, samples/pixel (multiple of 16), image width,
fluid steps, and optional mode (`fixed`, `exact`, or `adaptive`). Height is 5/8 of width. An evolved water/air droplet with an
internal bubble supplies the interface and opacity snapshot. The polished
hot iron object, floor and studio lights are illustrative scene geometry.
The renderer samples the fluid state without modifying it.

CPU tracing uses a binned SAH BVH, cached ray reciprocals and triangle normals,
shared-edge triangle intersections, and wavelength-local medium coefficients.
The fixed renderer defaults to a maximum omitted thermal source of 1e-6
W m⁻² sr⁻¹ m⁻¹ per weighted history, including roulette weights. Attenuation
and scattering remain unchanged. `Scene::radiance` and `exact` mode evaluate
all thermal emission; `render_with_budget(..., 0.0)` selects exact emission
for an image. This bound applies to thermal omission, not Monte Carlo error,
mesh error or depth truncation. Packet transport retains its full source physics.

`adaptive` discards an independent pilot and chooses a fixed main sample
count per pixel. Only the main samples contribute to the image; this avoids
bias from stopping on the samples being averaged. `AdaptiveOptions` controls
pilot/minimum counts, absolute/relative luminance targets and thermal budget.
`Image` records actual main counts and total pilot paths. A pilot can miss rare
caustics, and costs extra work when most pixels reach the maximum. Fixed
sampling remains the default. [Measured CPU results](../../../analysis/cpu_optimization.md).

Limits: CPU only; no GPU tracer yet. The packet physics has n=1, static,
straight-ray transport with open boundaries; refractive interfaces are in
the camera estimator. No automatic 3D diffusion/packet handoff, M1 solver,
moving-medium frequency shifts, polarization, coherence, diffraction,
fluorescence, photoionization, semiconductor photocarriers or gamma-ray
cross sections. No measured ice/liquid-nitrogen presets yet. The camera has
ideal smooth dielectric interfaces and polished/diffuse surfaces, not a
roughness/microfacet or thin-film model. Caustics converge slowly. The
depth cap and ray offset are reported numerical limitations. Meshes follow
the resolved phase fraction, so a coarse fluid grid remains visible.

References: [PBRT volume transport](https://pbr-book.org/4ed/Volume_Scattering/Transmittance),
[Fresnel](https://pbr-book.org/4ed/Reflection_Models/Specular_Reflection_and_Transmission),
[CIE fit](https://jcgt.org/published/0002/02/01/), and [optical data](data/SOURCES.md).

Validation: 24 physics, camera and fluid-coupling tests. Planck integrates
to Stefan–Boltzmann within 2e-9 relative; closed cosine diffusion RMS errors
are 3.243e-4, 8.123e-5 and 2.032e-5 at 16/32/64 cells (second order). Packet
energy/momentum and thermal emission/recoil balance to roundoff. Camera
tests include emitting slabs, Lambert reflection, Kirchhoff emissivity,
closed dielectric radiance, and thread-count reproducibility.
