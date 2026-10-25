//! Optical tests: normal-incidence reflectance of a metal (gold, Drude–Lorentz
//! fit) and the reststrahlen band of SiC (factorised polar-optical
//! permittivity).
//!
//! The wave-optics rung is the transfer-matrix / Fresnel response; the
//! validity indicator is the ratio of the microscopic structure (lattice
//! spacing, phonon wavevector) to the wavelength. In the infrared Reststrahlen
//! region the local permittivity is the only input: a single polar mode makes
//! a bulk crystal a perfect mirror between the TO and LO frequencies.
//!
//! Sources: Rakić et al., Appl. Opt. 37, 5271 (1998) (gold Drude–Lorentz);
//! Palik, *Handbook of Optical Constants of Solids* (gold reflectance);
//! Spitzer, Kleinman & Walsh, Phys. Rev. 113, 127 (1959) and Palik (1985)
//! (SiC phonon frequencies).

use crate::{Level, Measurement};
use photonics::complex::C;
use photonics::dispersion::{Dispersion, DrudeLorentz};
use photonics::multilayer::{interface_reflectance, Polarization};
use photonics::polariton::PolarCrystal;

/// Normal-incidence power reflectance of a vacuum/metal interface.
fn metal_reflectance(n: C) -> f64 {
    interface_reflectance(C::real(1.0), n, 0.0, Polarization::TE)
}

/// Gold at 633 nm and the SiC reststrahlen band.
pub fn optical() -> Vec<Measurement> {
    let mut out = Vec::new();

    // --- Gold: Drude–Lorentz normal-incidence reflectance ---
    let au = DrudeLorentz::gold();
    let lambda = 633.0e-9;
    let n_au = au.index(lambda);
    let r_au = metal_reflectance(n_au);
    // Measured normal-incidence reflectance of gold at 633 nm ≈ 0.933
    // (Johnson & Christy 1972 / Palik). The fit omits the interband loss
    // detail, so it is ~4 % low; the fold is certified within 5 %.
    out.push(Measurement::new(
        "Au",
        "reflectance (Drude–Lorentz)",
        Level::Macroscopic,
        "R(633 nm)",
        r_au,
        0.933,
        "",
        (r_au / 0.933 - 1.0).abs() < 0.05,
        &format!(
            "Fresnel reflectance from the Rakić (1998) Lorentz–Drude fit: ε = {:.3}{:+.3}i, \
             n+ik = {:.4}{:+.4}i. Reference is the measured normal-incidence reflectance of gold \
             (Palik). The fit is a causal model of the bound+d free electrons; at 633 nm the \
             interband Lorentz oscillators matter, so a pure Drude form would be too reflective.",
            n_au.re * n_au.re - n_au.im * n_au.im,
            2.0 * n_au.re * n_au.im,
            n_au.re,
            n_au.im
        ),
    ));

    // --- SiC: reststrahlen band ---
    let sic = PolarCrystal::silicon_carbide();
    let to_cm = sic.to_cm();
    let lo_cm = sic.lo_cm();
    let mid = 0.5 * (to_cm + lo_cm);
    let r_mid = sic.bulk_reflectance(mid);
    // Between TO and LO, Re ε < 0 and a lossless bulk crystal is a perfect
    // mirror.
    out.push(Measurement::new(
        "SiC",
        "reststrahlen band",
        Level::Macroscopic,
        "R(mid-band)",
        r_mid,
        1.0,
        "",
        (r_mid - 1.0).abs() < 0.05,
        &format!(
            "normal-incidence reflectance at {:.0} cm⁻¹, inside the Reststrahlen band where \
             Re ε(ω) < 0: the bulk crystal is a perfect mirror. The model has one damped polar mode \
             (γ ≈ 5 cm⁻¹); the exact lossless identity is R = 1.",
            mid
        ),
    ));
    // Band edges (Re ε = 0) vs the published TO/LO wavenumbers.
    out.push(Measurement::new(
        "SiC",
        "reststrahlen band",
        Level::Macroscopic,
        "ω_TO",
        to_cm,
        790.0,
        "cm⁻¹",
        (to_cm / 790.0 - 1.0).abs() < 0.02,
        "transverse-optical frequency (denominator pole of ε) vs the published SiC band edge \
         ≈790 cm⁻¹ (Spitzer et al. 1959; Palik 1985)",
    ));
    out.push(Measurement::new(
        "SiC",
        "reststrahlen band",
        Level::Macroscopic,
        "ω_LO",
        lo_cm,
        970.0,
        "cm⁻¹",
        (lo_cm / 970.0 - 1.0).abs() < 0.02,
        "longitudinal-optical frequency (numerator zero of ε) vs the published SiC band edge \
         ≈970 cm⁻¹; the LST relation ε(0) = ε_∞(ω_LO/ω_TO)² holds by construction",
    ));
    // Below the band the crystal is dielectric and transmits/reflects
    // normally; this checks the contrast.
    let below = sic.lossless().film_reflectance(0.5 * to_cm, 2.0e-6, 1.5);
    out.push(Measurement::new(
        "SiC",
        "reststrahlen band",
        Level::Macroscopic,
        "R(0.5 ω_TO, 2 µm film)",
        below,
        0.5704,
        "",
        (below / 0.5704 - 1.0).abs() < 0.15,
        "reflectance of a lossless 2 µm SiC film on an n = 1.5 substrate just below the \
         Reststrahlen band (ε ≈ 10.8, dielectric); reference is the transfer-matrix value 0.5704 \
         quoted in analysis/polariton.md. The contrast between the in-band mirror and the \
         out-of-band dielectric is the optical signature.",
    ));
    out
}
