//! The incoherent (geometrical-optics) limit of a thin film.
//!
//! A plane-parallel film reflects a train of partial waves: the front-surface
//! reflection, then waves that have made one, two, ... round trips. Whether
//! they interfere depends on the source coherence length `L_c` relative to the
//! round-trip optical path `2 n d`. For `L_c >> n d` the film shows sharp
//! interference fringes (the `multilayer` rung); for `L_c << n d` the partial
//! waves are mutually incoherent and only intensities add, which is the model
//! used by geometrical optics.
//!
//! The crossover is made explicit by weighting each pair of paths by the
//! degree of coherence `gamma((m-m') 2 n d)`:
//!
//! ```text
//!   R = sum_{m,m'} A_m A_{m'}^* gamma((m-m') 2 n d).
//! ```
//!
//! As `L_c -> 0`, `gamma` becomes a Kronecker delta and `R -> sum_m |A_m|^2`,
//! the incoherent intensity sum. As `L_c -> infinity`, `R -> |sum_m A_m|^2`,
//! the coherent Airy result. The indicator that decides which rung to use is
//! `n d / L_c`: below about one the coherent model is required, above it the
//! incoherent closure is accurate.
use crate::complex::C;

/// Incoherent reflectance of a film with interface power reflectances `r1`,
/// `r2` and single-pass internal intensity factor `a = exp(-alpha d)`.
pub fn incoherent_film_reflectance(r1: f64, r2: f64, a: f64) -> f64 {
    let t1 = 1.0 - r1;
    let denom = 1.0 - r1 * r2 * a * a;
    r1 + t1 * t1 * r2 * a * a / denom
}

/// Incoherent transmittance of the same film.
pub fn incoherent_film_transmittance(r1: f64, r2: f64, a: f64) -> f64 {
    (1.0 - r1) * (1.0 - r2) * a / (1.0 - r1 * r2 * a * a)
}

/// Reflected partial-wave amplitudes for a lossless film. `r_front` is the
/// front interface coefficient, `r_back` the internal reflection at the rear,
/// `delta = 2 pi n d / lambda` the single-pass phase.
fn partial_amplitudes(r_front: f64, r_back: f64, delta: f64, max_terms: usize) -> Vec<C> {
    let t = 1.0 - r_front * r_front; // t1 t1' (Stokes, lossless)
    let mut out = Vec::with_capacity(max_terms);
    out.push(C::real(r_front));
    let rr = r_front * r_back;
    let mut amp = t * r_back;
    for m in 1..max_terms {
        if amp.abs() < 1e-18 {
            break;
        }
        // A_m = (1-r1^2) r2 (r1 r2)^{m-1} e^{2 i m delta}
        out.push(C::new(
            amp * (2.0 * m as f64 * delta).cos(),
            amp * (2.0 * m as f64 * delta).sin(),
        ));
        amp *= rr;
    }
    out
}

/// Incoherent intensity sum `sum_m |A_m|^2`.
pub fn incoherent_slab_reflectance(r_front: f64, r_back: f64) -> f64 {
    let amps = partial_amplitudes(r_front, r_back, 0.0, 4000);
    amps.iter().map(|a| a.abs2()).sum()
}

/// Fully coherent Airy reflectance.
pub fn coherent_slab_reflectance(r_front: f64, r_back: f64, delta: f64) -> f64 {
    let amps = partial_amplitudes(r_front, r_back, delta, 4000);
    let sum: C = amps.iter().fold(C::real(0.0), |acc, a| acc + *a);
    sum.abs2()
}

/// Partially coherent reflectance for a Gaussian degree of coherence
/// `gamma(L) = exp(-(L/L_c)^2)`. `optical_path = 2 n d` sets the round-trip
/// path difference; `coherence_length` is `L_c` in the same units.
pub fn partially_coherent_slab_reflectance(
    r_front: f64,
    r_back: f64,
    delta: f64,
    optical_path: f64,
    coherence_length: f64,
) -> f64 {
    assert!(coherence_length > 0.0);
    let amps = partial_amplitudes(r_front, r_back, delta, 4000);
    let mut total = 0.0;
    for (m, a) in amps.iter().enumerate() {
        for (p, b) in amps.iter().enumerate() {
            let path = (m as f64 - p as f64) * optical_path;
            let g = (-(path / coherence_length).powi(2)).exp();
            total += (*a * b.conj()).re * g;
        }
    }
    total
}

/// `n d / L_c`, the coherence indicator. Values below one require the coherent
/// rung; much above one certify the incoherent closure.
pub fn coherence_indicator(n: f64, thickness: f64, coherence_length: f64) -> f64 {
    n * thickness / coherence_length
}

/// A thin-film `Layer`-free convenience: Fresnel amplitude coefficient at a
/// real dielectric interface.
pub fn amplitude_r(n0: f64, n1: f64) -> f64 {
    (n0 - n1) / (n0 + n1)
}
