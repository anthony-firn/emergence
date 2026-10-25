//! Light → heat → sound: the coupling across the emergent boundaries of the
//! photothermal chain.
//!
//! [`photonics`] resolves the optical field and returns an absorbed power
//! density; [`phononic`] resolves acoustic propagation in a layered stack. This
//! crate is the rung that connects them to heat and to a photoacoustic source,
//! and it makes the hand-offs explicit as indicator functions.
//!
//! # The chain
//!
//! ```text
//!   optics: field (photonics::multilayer)  --p(z)-->  thermal: conduction
//!                                                          |
//!                            stress confinement             v
//!   photoacoustics:  p0 = Gamma mu_a F  --d'Alembert-->  acoustic stack
//!                                        (phononic transfer matrix / FDTD)
//! ```
//!
//! * [`thermal`]: conservative, layer-aligned, implicit 1D heat conduction
//!   with a volumetric source. Benchmarks against the semi-infinite
//!   surface-flux erfc solution, the Beer–Lambert steady state and a Fourier
//!   series transient.
//! * [`photoacoustic`]: the exact d'Alembert split of the initial pressure for
//!   free (bipolar) and rigid (monopolar) surfaces, a staggered-grid time-domain
//!   solver for layered media, and a frequency-domain launch through
//!   [`phononic::layered`] transfer matrices.
//! * [`indicators`]: thermal and stress confinement, the coherent/incoherent
//!   optics boundary, and the Rayleigh number for the onset of convection.
//!
//! # Folding and unfolding
//!
//! The organising idea is in `analysis/emergence_ladder.md`; the photothermal
//! instance is written up in `analysis/photothermal.md`. In short:
//!
//! | folded model | indicator | unfold when |
//! |---|---|---|
//! | Beer–Lambert source from the optical absorptance | layer thickness vs absorption length | absorption length < layer thickness ([`photonics::multilayer::Solution::power_density`]) |
//! | instantaneous heat deposition | thermal confinement `tau_p alpha / delta^2` | > 1 ([`thermal::ThermalModel`]) |
//! | instantaneous pressure deposition | stress confinement `tau_p c / delta` | > 1 ([`photoacoustic`]) |
//! | conduction-only transport | Rayleigh number | > 1708 ([`indicators`]) |
pub mod fft;
pub mod indicators;
pub mod photoacoustic;
pub mod thermal;

pub use indicators::{
    absorption_depth, convects_internal, convects_rayleigh_benard, is_stress_confined,
    is_thermally_confined, rayleigh_benard, rayleigh_from_flux, rayleigh_internal,
    stress_confinement_ratio, thermal_confinement_ratio, OpticalRegime, RA_CRIT_INTERNAL,
    RA_CRIT_RAYLEIGH_BENARD,
};
pub use photoacoustic::{
    dalembert, gruneisen, initial_pressure, spectral_transmit, AcousticLayer, AcousticModel,
    FarBoundary, Surface, GRUNEISEN_WATER,
};
pub use thermal::{
    beer_lambert_cell_average, beer_lambert_source, semi_infinite_surface_flux, source_from_optics,
    Boundary, Material, ThermalLayer, ThermalModel, FUSED_SILICA, GOLD, MGF2, PYREX, WATER,
};
