# Young's double slit, recorded by physical devices

Crate `fluid/crates/waveoptics`. The goal is a double-slit experiment in
which every step from emitter to recorded count is a physical model with a
stated error. That means a partially coherent, polarized source; certified
scalar diffraction; and a real recording device. The devices are a bare
sensor, a camera with a vendor lens prescription, a pinhole camera, or the
human eye with its cone mosaic. Photons are detected with Poisson
statistics and silicon quantum efficiency. Every claim below comes from a
test in `cargo test -p waveoptics` (each test under 15 s) or from the scene
program, and each states its error.

Run the scenes with `cargo run --release -p waveoptics --example
double_slit_scenes` (about 4 minutes on CPU, under 0.4 GB). It writes PNGs
and `report.txt` to `fluid/renders/double_slit/`.

## The ladder of models

| Rung | What it computes | Error control |
|---|---|---|
| Tilted Fraunhofer / tilted Fresnel | Field of one rectangular opening at one point. Exact phase, gradient and diagonal Hessian at the opening centre; closed-form integrals | Rigorous bound on everything dropped (cross term, cubic remainder, amplitude curvature) |
| Line integral | Elongated openings: exact phase along the long axis, local expansion across | Validated bound |
| Rayleigh–Sommerfeld quadrature | The reference | Refinement estimate |
| **Angular spectrum, band-limited** (new) | Whole planes: each opening's closed-form angular spectrum, propagated exactly, restricted to the band that can reach the target (a sensor grid, or the directions a device's detector window accepts) | Out-of-band power is reported. Incident-wave model bound for point sources (2.9e-6 rad in the test geometry). Audited against the certified ladder |
| Device pupil + Fourier imaging | Traced wavefront of the real optics per wavelength; matrix DFT onto the detector window | Pupil step keeps every replica out of the window; the device's own ray-angle errors widen the band |
| Isoplanatic coherent imaging (FFT) | A lens or pinhole camera photographing a diffuse screen. The pinhole is the same imager with pupil phase `|u|²(1/D + 1/f)/2` | Paraxial; validated on speckle statistics and the coherent-edge response |

Which rung to use is decided by cost per point. The pointwise ladder is
certified but costs about 1 ms per point in the near field, where it
escalates to the line integral. The angular-spectrum rung computes a
whole grid with two separable sums, and a strided sample of its points is
re-solved by the certified ladder as an audit. On a 640 × 480 sensor 0.3 m
behind the slits, the same frame took **396 s through the pointwise ladder
and 0.8 s through the grid rung** (audit 3.1e-5 of the field scale).

## Validation

**Diffraction and coherence**: `tests/diffraction.rs`, `tests/double_slit.rs`.

* Absolute centre irradiance: 1.2032 vs 1.2036 W/m² in closed form.
* Pattern within 3.2e-4 of cos²·sinc², and fringe spacing within 3e-5.
* Finite-slit van Cittert–Zernike visibility within 3.3e-4.
* Bandwidth fringes within 2.6e-4 of an independent two-beam integral.
* Which-path marking and the quantum eraser with polarizers.
* **New:** the grid rung matches the certified ladder at every point of a
  sensor-sized grid to **1.9e-5 of peak irradiance**, for partially
  coherent light in the near field. It is 11× faster even at 656 points.

**Angular spectrum vs ladder** (`tests/devices.rs`), at 20 mm, Fresnel
number 0.8:

* Tilted plane wave: agrees to 1.1e-9 of the field scale.
* Off-axis point source: agrees to 7.5e-8.

Writing this test exposed a sign error in the diffraction module. Its code
used `e^{−ik s·ξ}` for a tilted plane wave, contrary to its own
documentation and inconsistent with point sources. It was invisible
because sources only ever produced untilted plane waves. It is now fixed.

**Devices** (`tests/devices.rs`):

| Check | Result |
|---|---|
| Ideal camera, 20 mm behind the slits, records the exact Fraunhofer pattern `|a h sinc sinc 2cos|²/(λf)²` | max error 9.9e-5 of peak |
| Fringe period on the detector vs `f λ/d` | 63.281 vs 63.280 µm |
| Camera fringe visibility vs the slits' mutual coherence (extended source) | 0.8433 vs 0.8440, 0.3535 vs 0.3554 |
| Rough screen: fully developed speckle, contrast | 0.985, with ~8900 speckles (σ ≈ 0.011) |
| Speckle mean vs the pupil's share `π f_c² δ²` | −1.4 % |
| Coherent straight edge, Gibbs overshoot `(½ + Si(π)/π)²` | 1.1871 vs 1.1873 |
| Cone catches vs `E λ/(hc) · A · s_Q(λ) · t` | exact (to 1e-9) |
| Variance/mean of cone catches (Poisson) | 1.00 |
| Mosaic reconstruction of a uniform broadband field vs its CIE 2006 XYZ | within 0.05 % |

Earlier tests for the optics stay as they were:

* The glasses match the Schott catalog.
* The AC254-050-A matches its vendor EFL/BFD.
* The model eye is 60.68 D and emmetropic to 11 µm.
* With physiological dispersion, its chromatic aberration is within 0.08 D
  of Atchison & Smith and 0.06 D of Thibos.
* The Airy pattern, encircled energy, Maréchal limit and pinhole limits
  hold.

## The scenes

Slits are 50 µm wide, 250 µm apart and 3 mm tall. The laser is a He–Ne:
1 mW expanded to a 5 mm beam (estimate). The white LED is a
phosphor-converted spectrum: a 450 nm pump plus a 560 nm phosphor, 1 : 3,
from a 1 mm emitter 2 m away (estimate). The full numbers for each scene
are in `report.txt`.

**a) Bare CMOS at the screen plane (L = 0.3 m)**

* Raw 12-bit frame with photon shot noise, dark current, read noise and
  full-well clipping.
* Single-photon build-up frames at 20 … 2 × 10⁶ detections: the pattern
  emerges from individual events.
* Fringe period measured from dark-fringe positions: +0.05 % from λL/d.

**b) Eye looking through the slits** (held 20 mm in front of the cornea)
at the LED. The pupil is 3 mm; the model is Escudero-Sanz & Navarro.

* Outputs: the retinal image (colorimetric), the photoisomerisations of
  every cone in 0.1 s on a 2.5 µm triangular mosaic, and the colour that
  mosaic encodes.
* An aberration-free lens with the eye's pupil and focal length
  reproduces the closed-form two-beam white-light pattern to **2.4e-3 of
  peak**, with centre visibility 0.868 vs 0.872.
* The eye's own chromatic aberration then changes the pattern by 4.4 %
  of peak with physiological dispersion and 6.0 % with the published
  (more dispersive) Navarro media.
* The cone-encoded image has no blue signal in the central fovea, because
  the foveola has no S cones. No perceptual filling-in is modelled.

**c) Camera: AC254-050-A at f/8 and a CMOS sensor**, in the eye's
position.

* The expected-charge profile is within 4.8e-3 of the closed form,
  including the lens's residual aberrations and silicon QE.
* Centre visibility: 0.868 in the raw frame vs 0.875 in closed form.

**d) Laser fringes on a matte screen** (L = 1 m, albedo 0.9, estimate),
photographed by two cameras:

* **f/8 lens camera at 0.5 m:**
  * Fringe period +0.09 % from m λL/d (speckle-free mean).
  * Speckle contrast 0.87 ± 0.02 vs 0.865 for fully developed speckle
    integrated over the 3.45 µm pixels (Goodman).
* **0.25 mm pinhole camera at 0.15 m:**
  * Fringe period +0.28 %.
  * Speckle contrast ≈ 1 ± 0.1 (only ~70 speckles in view).
  * The fringes are largely washed out by its 1 mm blur on the screen,
    as a real pinhole would show.

## Limits and estimates

* **Scalar diffraction only.** Openings must be at least 10 λ wide;
  Kirchhoff boundary values are used.
* **Device optics:**
  * Each device is focused at infinity, or isoplanatic and paraxial for
    the screen scene.
  * The eye's retina is its tangent plane.
  * The traced pupil is used on axis for the whole window, i.e. the
    window is assumed isoplanatic. That holds for the ≲1° windows used
    here, but is not checked.
* **Estimates, not measurements:**
  * Retina: cone collecting area, proportions, the S-free radius and the
    0.1 s integration time.
  * Sensor parameters.
  * The LED spectrum and screen albedo.
  * Absolute photoisomerisation rates are therefore estimates; only the
    ratios and statistics are validated.
* **Rough screen:** modelled as uniformly random phase per sample, which
  is fully developed speckle. Surface-height statistics, wavelength
  decorrelation and subsurface scattering are not modelled.
* **White-light fringe period:** it is not one number. A dark-fringe
  regression drifts toward the narrow blue LED line, whose fringes
  survive more orders. For white light the scenes therefore compare whole
  profiles with the closed form instead.
* **Sensor colour:** colour images from cameras are ideal colorimetric
  (CIE XYZ), noise-free. Only the monochrome raw frame carries sensor
  noise; no colour filter array is modelled.
