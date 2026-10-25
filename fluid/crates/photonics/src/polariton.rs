//! Phonon polaritons in polar crystals: the electromagnetic closure of the
//! lattice vibrations.
//!
//! A polar diatomic crystal supports a transverse-optical (TO) phonon at
//! `omega_TO` and a longitudinal-optical (LO) phonon at `omega_LO`, split by the
//! macroscopic electric field the ions produce when they move against each
//! other. At long wavelength the crystal's infrared response is captured by a
//! frequency-dependent permittivity
//!
//! ```text
//!   eps(omega) = eps_inf (omega_LO^2 - omega^2 - i gamma_LO omega)
//!                        / (omega_TO^2 - omega^2 - i gamma_TO omega).
//! ```
//!
//! With `gamma_TO = gamma_LO = gamma` this is the single Lorentz oscillator;
//! with independent damping of the two poles it is the factorised
//! (Gervais-Piriou) form. Its static limit is the Lyddane-Sachs-Teller (LST)
//! relation `eps(0)/eps_inf = (omega_LO/omega_TO)^2`, which the tests check
//! exactly. The denominator pole sits at the TO frequency the lattice
//! dynamical matrix predicts; the numerator zero sits at the LO frequency the
//! Born effective charge produces through the non-analytic (depolarisation)
//! term. That matching is the fold: the photons see a continuum permittivity
//! whose poles are the zone-centre lattice modes.
//!
//! The rungs this module exposes:
//!
//! * [`PolarCrystal`] — the sourced infrared permittivity of SiC, NaCl, GaAs
//!   and the in-plane response of hBN, in the factorised form.
//! * Reststrahlen reflectance — bulk (exact Fresnel) and thin films through
//!   [`crate::multilayer`]; `R -> 1` between `omega_TO` and `omega_LO`.
//! * Bulk polariton dispersion `k = (omega/c) sqrt(eps(omega))` (two branches
//!   with a gap between TO and LO) and the surface phonon polariton
//!   `k_sp = (omega/c) sqrt(eps/(eps+1))`, which exists where `eps < -1` and
//!   diverges at the surface-mode frequency `omega_s` where `eps = -1`. The
//!   Otto/Kretschmann attenuated-total-reflection dip computed with the
//!   transfer-matrix code recovers `omega_s` and the analytic dip angle.
//! * The fold from lattice dynamics: [`DiatomicChain`] gives the zone-centre
//!   TO frequency from a 2x2 dynamical matrix, and the Born effective charge
//!   turns it into the LO frequency and the static permittivity. The cheap
//!   indicators `q a` (spatial dispersion) and `a/d` (film thickness in lattice
//!   spacings) say when the local `eps(omega)` closure may be used.
//!
//! Frequencies are handled in angular units rad/s internally and in
//! spectroscopic wavenumbers cm^-1 in the public entry points, since the
//! literature parameters are quoted in cm^-1.

use crate::complex::C;
use crate::dispersion::Dispersion;
use crate::multilayer::{Layer, Polarization, Stack};
use std::f64::consts::PI;

/// Converts a spectroscopic wavenumber `nu_tilde` (cm^-1) to angular
/// frequency `omega = 2 pi c nu_tilde` with `nu_tilde` converted to m^-1.
pub const CM_TO_RAD: f64 = 2.0 * PI * crate::C0 * 100.0;

/// Elementary charge [C], CODATA 2018 (exact by SI definition).
pub const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19;
/// Atomic mass unit [kg], CODATA 2018.
pub const AMU: f64 = 1.660_539_066_60e-27;

/// Spectroscopic wavenumber (cm^-1) to angular frequency (rad/s).
pub fn cm_to_rad(nu_cm: f64) -> f64 {
    nu_cm * CM_TO_RAD
}
/// Angular frequency (rad/s) to spectroscopic wavenumber (cm^-1).
pub fn rad_to_cm(omega: f64) -> f64 {
    omega / CM_TO_RAD
}
/// Vacuum wavelength (m) of a spectroscopic wavenumber (cm^-1).
pub fn wavenumber_to_lambda(nu_cm: f64) -> f64 {
    1.0e-2 / nu_cm
}
/// Spectroscopic wavenumber (cm^-1) of a vacuum wavelength (m).
pub fn lambda_to_wavenumber(lambda_m: f64) -> f64 {
    1.0e-2 / lambda_m
}

// ---------------------------------------------------------------------------
// The infrared dielectric function
// ---------------------------------------------------------------------------

/// One polar phonon mode: the transverse and longitudinal frequencies and the
/// damping of the two factorised poles, all in cm^-1.
#[derive(Clone, Copy, Debug)]
pub struct PhononMode {
    pub to_cm: f64,
    pub lo_cm: f64,
    pub gamma_to_cm: f64,
    pub gamma_lo_cm: f64,
}

impl PhononMode {
    /// `to_cm` and `lo_cm` in cm^-1; the LO frequency must exceed the TO one.
    pub fn new(to_cm: f64, lo_cm: f64, gamma_to_cm: f64, gamma_lo_cm: f64) -> Self {
        assert!(to_cm > 0.0 && lo_cm > to_cm, "need 0 < TO < LO");
        assert!(gamma_to_cm >= 0.0 && gamma_lo_cm >= 0.0);
        Self {
            to_cm,
            lo_cm,
            gamma_to_cm,
            gamma_lo_cm,
        }
    }
    /// Lossless mode (`gamma = 0`); the LST relation and the pole positions are
    /// unchanged by damping, so the exact identities use this.
    pub fn lossless(to_cm: f64, lo_cm: f64) -> Self {
        Self::new(to_cm, lo_cm, 0.0, 0.0)
    }
    /// The factorised factor `(omega_LO^2 - omega^2 - i gamma_LO omega) /
    /// (omega_TO^2 - omega^2 - i gamma_TO omega)` at angular frequency omega.
    pub fn factor(&self, omega: f64) -> C {
        let wto = cm_to_rad(self.to_cm);
        let wlo = cm_to_rad(self.lo_cm);
        let gto = cm_to_rad(self.gamma_to_cm);
        let glo = cm_to_rad(self.gamma_lo_cm);
        let num = C::real(wlo * wlo - omega * omega) - C::i() * (glo * omega);
        let den = C::real(wto * wto - omega * omega) - C::i() * (gto * omega);
        num / den
    }
}

/// A polar crystal described by its high-frequency permittivity and one or
/// more factorised polar modes. Sources are cited on the constructors.
#[derive(Clone, Debug)]
pub struct PolarCrystal {
    pub name: &'static str,
    pub eps_inf: f64,
    pub modes: Vec<PhononMode>,
    /// Independently measured static permittivity (handbook), for the LST
    /// consistency check. `f64::NAN` if none is quoted.
    pub eps0_measured: f64,
    pub source: &'static str,
}

impl PolarCrystal {
    /// Complex relative permittivity at angular frequency `omega` (rad/s),
    /// factorised (Gervais-Piriou) over all modes.
    pub fn epsilon(&self, omega: f64) -> C {
        let mut eps = C::real(self.eps_inf);
        for m in &self.modes {
            eps = eps * m.factor(omega);
        }
        eps
    }
    /// Complex relative permittivity at spectroscopic wavenumber `nu_cm`.
    pub fn epsilon_wavenumber(&self, nu_cm: f64) -> C {
        self.epsilon(cm_to_rad(nu_cm))
    }
    /// Complex refractive index `sqrt(eps)` at vacuum wavelength `lambda_m`.
    /// The principal square root gives `Im(n) >= 0` for loss and `n` imaginary
    /// inside a Reststrahlen band, which the transfer-matrix code needs.
    pub fn index_at_wavelength(&self, lambda_m: f64) -> C {
        let omega = 2.0 * PI * crate::C0 / lambda_m;
        self.epsilon(omega).sqrt()
    }
    /// Static permittivity `eps(0) = eps_inf prod (omega_LO/omega_TO)^2`.
    pub fn static_epsilon(&self) -> f64 {
        self.eps_inf
            * self
                .modes
                .iter()
                .map(|m| (m.lo_cm / m.to_cm).powi(2))
                .product::<f64>()
    }
    /// LST ratio `eps(0)/eps_inf`.
    pub fn lst_ratio(&self) -> f64 {
        self.static_epsilon() // already eps_inf * ratio
            / self.eps_inf
    }
    /// Lossless copy; the exact identities (LST, perfect Reststrahlen,
    /// undamped polariton branches) use it.
    pub fn lossless(&self) -> Self {
        let mut c = self.clone();
        for m in &mut c.modes {
            m.gamma_to_cm = 0.0;
            m.gamma_lo_cm = 0.0;
        }
        c
    }
    /// Number of polar modes.
    pub fn n_modes(&self) -> usize {
        self.modes.len()
    }
    /// Convenience for a single-mode crystal: the TO frequency in cm^-1.
    pub fn to_cm(&self) -> f64 {
        self.modes[0].to_cm
    }
    /// Convenience for a single-mode crystal: the LO frequency in cm^-1.
    pub fn lo_cm(&self) -> f64 {
        self.modes[0].lo_cm
    }

    // -- sourced materials ---------------------------------------------------

    /// Cubic (3C) silicon carbide. Reststrahlen band ~796-973 cm^-1.
    /// TO/LO and eps_inf from the infrared lattice bands; eps0 is the
    /// independent static measurement. Source: E. D. Palik, *Handbook of
    /// Optical Constants of Solids* (Academic Press, 1985), SiC entry, and
    /// W. G. Spitzer, D. Kleinman, D. Walsh, Phys. Rev. **113**, 127 (1959).
    pub fn silicon_carbide() -> Self {
        Self {
            name: "3C-SiC",
            eps_inf: 6.52,
            modes: vec![PhononMode::new(796.2, 972.7, 5.0, 5.0)],
            eps0_measured: 9.72,
            source: "Palik, Handbook of Optical Constants of Solids (1985); \
                     Spitzer, Kleinman & Walsh, Phys. Rev. 113, 127 (1959)",
        }
    }

    /// Rock-salt sodium chloride. Source: E. D. Palik, *Handbook of Optical
    /// Constants of Solids*, Vol. 2 (Academic Press, 1991), NaCl entry; the
    /// mode frequencies are the classic Born-Huang values.
    pub fn sodium_chloride() -> Self {
        Self {
            name: "NaCl",
            eps_inf: 2.34,
            modes: vec![PhononMode::new(164.0, 264.0, 5.0, 5.0)],
            eps0_measured: 5.90,
            source: "Palik, Handbook of Optical Constants of Solids II (1991); \
                     Born & Huang, Dynamical Theory of Crystal Lattices (1954)",
        }
    }

    /// Gallium arsenide. Source: J. S. Blakemore, J. Appl. Phys. **53**, R123
    /// (1982), and Palik, *Handbook of Optical Constants of Solids* (1985).
    pub fn gallium_arsenide() -> Self {
        Self {
            name: "GaAs",
            eps_inf: 10.9,
            modes: vec![PhononMode::new(268.3, 292.3, 2.0, 2.0)],
            eps0_measured: 12.9,
            source: "Blakemore, J. Appl. Phys. 53, R123 (1982); Palik (1985)",
        }
    }

    /// Hexagonal boron nitride, in-plane (E perpendicular to c) response.
    /// Source: R. Geick, C. H. Perry, G. Rupprecht, Phys. Rev. **146**, 543
    /// (1966).
    pub fn hbn_inplane() -> Self {
        Self {
            name: "hBN (in-plane)",
            eps_inf: 4.95,
            modes: vec![PhononMode::new(1367.0, 1580.0, 10.0, 10.0)],
            eps0_measured: 6.70,
            source: "Geick, Perry & Rupprecht, Phys. Rev. 146, 543 (1966)",
        }
    }

    // -- Reststrahlen reflectance -------------------------------------------

    /// Normal-incidence power reflectance of a bulk (semi-infinite) crystal,
    /// exact Fresnel from the complex permittivity. `R -> 1` where `eps < 0`.
    pub fn bulk_reflectance(&self, nu_cm: f64) -> f64 {
        let eps = self.epsilon_wavenumber(nu_cm);
        crate::multilayer::interface_reflectance(C::real(1.0), eps.sqrt(), 0.0, Polarization::TE)
    }

    /// Normal-incidence reflectance of a film of `thickness_m` on a
    /// non-absorbing substrate of index `substrate_index`, through the
    /// transfer-matrix code. The film is described by `eps(omega)`.
    pub fn film_reflectance(&self, nu_cm: f64, thickness_m: f64, substrate_index: f64) -> f64 {
        let eps = self.epsilon_wavenumber(nu_cm);
        let layer = Layer::from_eps_mu(eps, C::real(1.0), thickness_m);
        Stack::new(1.0, vec![layer], substrate_index)
            .solve(wavenumber_to_lambda(nu_cm), 0.0, Polarization::TE)
            .reflectance
    }

    // -- polariton dispersion ------------------------------------------------

    /// Bulk (transverse) polariton wavevector `k = (omega/c) sqrt(eps)`.
    /// Real `k` on the two propagating branches (`eps > 0`); imaginary inside
    /// the Reststrahlen gap between TO and LO.
    pub fn bulk_polariton_wavevector(&self, omega: f64) -> C {
        C::real(omega / crate::C0) * self.epsilon(omega).sqrt()
    }

    /// Surface phonon polariton wavevector at a crystal/vacuum interface,
    /// `k_sp = (omega/c) sqrt(eps/(eps+1))`. Real and `> omega/c` where
    /// `eps < -1`; imaginary between `omega_s` (`eps = -1`) and `omega_LO`.
    pub fn surface_polariton_wavevector(&self, omega: f64) -> C {
        let eps = self.epsilon(omega);
        C::real(omega / crate::C0) * (eps / (eps + C::real(1.0))).sqrt()
    }

    /// Surface-mode frequency where the lossless single-mode permittivity
    /// equals -1, in cm^-1:
    /// `omega_s^2 = (eps_inf omega_LO^2 + omega_TO^2)/(1 + eps_inf)`.
    pub fn surface_mode_wavenumber(&self) -> f64 {
        assert_eq!(self.modes.len(), 1, "analytic form is single-mode");
        let m = &self.modes[0];
        let wto = cm_to_rad(m.to_cm);
        let wlo = cm_to_rad(m.lo_cm);
        let ws2 = (self.eps_inf * wlo * wlo + wto * wto) / (1.0 + self.eps_inf);
        rad_to_cm(ws2.sqrt())
    }

    /// Numeric surface-mode frequency from `Re eps(omega) = -1`, used as an
    /// independent check of [`Self::surface_mode_wavenumber`] and valid with
    /// damping. Assumes a single Reststrahlen band.
    pub fn surface_mode_wavenumber_numeric(&self) -> f64 {
        assert_eq!(self.modes.len(), 1, "bisection assumes one band");
        let m = &self.modes[0];
        let mut lo = m.to_cm * (1.0 + 1e-12);
        let mut hi = m.lo_cm * (1.0 - 1e-12);
        let f = |x: f64| self.epsilon_wavenumber(x).re + 1.0;
        // Re eps -> -inf just above TO, -> 0- just below LO.
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if f(mid) < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Wavenumber at which the surface polariton has fixed tangential
    /// wavevector `s = k_x/k0 > 1`, i.e. `Re eps(omega) = s^2/(1 - s^2)`.
    /// This is the frequency of the Otto/Kretschmann reflectance dip at that
    /// angle. Single Reststrahlen band.
    pub fn surface_mode_wavenumber_at_s(&self, s: f64) -> f64 {
        assert!(s > 1.0, "only evanescent coupling has s > 1");
        assert_eq!(self.modes.len(), 1);
        let target = s * s / (1.0 - s * s); // negative
        let m = &self.modes[0];
        let mut lo = m.to_cm * (1.0 + 1e-12);
        let mut hi = m.lo_cm * (1.0 - 1e-12);
        let f = |x: f64| self.epsilon_wavenumber(x).re - target;
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if f(mid) < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Analytic Otto/Kretschmann dip sine, `sin(theta) =
    /// sqrt(eps/(eps+1))/n_prism`, from the complex permittivity.
    pub fn surface_polariton_sine(&self, prism_index: f64, nu_cm: f64) -> f64 {
        let eps = self.epsilon_wavenumber(nu_cm);
        (eps / (eps + C::real(1.0))).sqrt().re / prism_index
    }

    /// p-polarised reflectance of prism / air gap / bulk crystal (Otto
    /// configuration) at wavenumber `nu_cm` and incidence angle `theta`.
    pub fn otto_reflectance(
        &self,
        prism_index: f64,
        gap_m: f64,
        nu_cm: f64,
        theta: f64,
    ) -> f64 {
        let eps = self.epsilon_wavenumber(nu_cm);
        let lambda = wavenumber_to_lambda(nu_cm);
        let s = prism_index * theta.sin();
        let stack = Stack::with_complex(
            C::real(prism_index),
            vec![Layer::new(1.0, 0.0, gap_m)],
            eps.sqrt(),
        );
        stack.solve_s(lambda, s, Polarization::TM).reflectance
    }

    /// Search the Otto reflectance minimum over `[theta_lo, theta_hi]` by
    /// ternary search. Returns `(angle, reflectance)`.
    pub fn otto_dip_angle(
        &self,
        prism_index: f64,
        gap_m: f64,
        nu_cm: f64,
        theta_lo: f64,
        theta_hi: f64,
    ) -> (f64, f64) {
        let f = |t: f64| self.otto_reflectance(prism_index, gap_m, nu_cm, t);
        let (mut a, mut b) = (theta_lo, theta_hi);
        for _ in 0..200 {
            let m1 = a + (b - a) / 3.0;
            let m2 = b - (b - a) / 3.0;
            if f(m1) < f(m2) {
                b = m2;
            } else {
                a = m1;
            }
        }
        let t = 0.5 * (a + b);
        (t, f(t))
    }

    /// Search the Otto reflectance minimum in wavenumber at fixed angle
    /// `theta` over `[nu_lo, nu_hi]`. Returns `(nu_cm, reflectance)`.
    pub fn otto_dip_wavenumber(
        &self,
        prism_index: f64,
        gap_m: f64,
        theta: f64,
        nu_lo: f64,
        nu_hi: f64,
    ) -> (f64, f64) {
        let f = |x: f64| self.otto_reflectance(prism_index, gap_m, x, theta);
        let (mut a, mut b) = (nu_lo, nu_hi);
        for _ in 0..200 {
            let m1 = a + (b - a) / 3.0;
            let m2 = b - (b - a) / 3.0;
            if f(m1) < f(m2) {
                b = m2;
            } else {
                a = m1;
            }
        }
        let x = 0.5 * (a + b);
        (x, f(x))
    }
}

impl Dispersion for PolarCrystal {
    fn index(&self, lambda_m: f64) -> C {
        self.index_at_wavelength(lambda_m)
    }
}

// ---------------------------------------------------------------------------
// The fold: lattice dynamics -> the dielectric function
// ---------------------------------------------------------------------------

/// A one-dimensional diatomic ionic chain: two masses per cell separated by
/// `a/2`, nearest-neighbour harmonic spring `k`, cell length `a`. This is the
/// exactly solvable model whose zone-centre optical mode becomes the pole of
/// the dielectric function. It mirrors the `Harmonic` spring and `chain_lattice`
/// construction of the `materials` crate, reimplemented here so this crate has
/// no dependency on it.
#[derive(Clone, Copy, Debug)]
pub struct DiatomicChain {
    /// Unit-cell length [m].
    pub a: f64,
    pub mass1: f64,
    pub mass2: f64,
    /// Nearest-neighbour spring constant [N/m].
    pub spring: f64,
}

impl DiatomicChain {
    pub fn new(a: f64, mass1: f64, mass2: f64, spring: f64) -> Self {
        assert!(a > 0.0 && mass1 > 0.0 && mass2 > 0.0 && spring > 0.0);
        Self {
            a,
            mass1,
            mass2,
            spring,
        }
    }

    /// Squared phonon frequencies `(omega_acoustic^2, omega_optical^2)` at
    /// wavevector `q` (rad/m) along the chain. Derived from the 2x2 dynamical
    /// matrix with nearest-neighbour bonds at `+a/2` and `-a/2`.
    pub fn omega_sq(&self, q: f64) -> (f64, f64) {
        let k = self.spring;
        let (m1, m2) = (self.mass1, self.mass2);
        let mean = k * (1.0 / m1 + 1.0 / m2);
        let split =
            ((1.0 / m1 - 1.0 / m2).powi(2) + 4.0 * (q * self.a / 2.0).cos().powi(2) / (m1 * m2))
                .sqrt();
        let root = k * split;
        ((mean - root).max(0.0), mean + root)
    }

    /// Zone-centre optical (TO) frequency [rad/s]:
    /// `omega_TO^2 = 2 k (1/m1 + 1/m2)`.
    pub fn omega_to(&self) -> f64 {
        (2.0 * self.spring * (1.0 / self.mass1 + 1.0 / self.mass2)).sqrt()
    }

    /// Zone-centre TO wavenumber [cm^-1].
    pub fn to_wavenumber(&self) -> f64 {
        rad_to_cm(self.omega_to())
    }

    /// Optical frequency at the Brillouin-zone boundary `q = pi/a` [rad/s].
    pub fn omega_optical_zone_boundary(&self) -> f64 {
        self.omega_sq(PI / self.a).1.sqrt()
    }

    /// Reduced mass `m1 m2/(m1+m2)`.
    pub fn reduced_mass(&self) -> f64 {
        self.mass1 * self.mass2 / (self.mass1 + self.mass2)
    }
}

/// Oscillator strength of an ionic mode with Born effective charge `z_star`
/// (in units of `e`), ion-pair number density `n_pairs` [m^-3] and reduced
/// mass `mu` [kg]:
///
/// ```text
///   S = z_star^2 e^2 n_pairs / (eps0 mu)     [rad^2/s^2].
/// ```
///
/// This is the coefficient of the `(Z^* e) E` driving term in the relative
/// displacement equation of motion, so `eps(0) = eps_inf + S/omega_TO^2` and
/// the depolarisation field shifts the longitudinal mode to
/// `omega_LO^2 = omega_TO^2 + S/eps_inf` — LST follows identically.
pub fn born_oscillator_strength(z_star: f64, n_pairs: f64, reduced_mass: f64) -> f64 {
    z_star * z_star * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE * n_pairs
        / (crate::EPS0 * reduced_mass)
}

/// Builds the long-wavelength dielectric closure of a diatomic ionic chain:
/// the TO pole from the dynamical matrix, the LO frequency from the Born
/// charge through the non-analytic term, and the static permittivity that
/// follows. `n_pairs` is the ion-pair density and `z_star` the Born effective
/// charge; `gamma_cm` is a phenomenological damping (0 for the exact test).
pub fn polar_crystal_from_chain(
    chain: &DiatomicChain,
    eps_inf: f64,
    z_star: f64,
    n_pairs: f64,
    gamma_cm: f64,
) -> PolarCrystal {
    let wto = chain.omega_to();
    let s = born_oscillator_strength(z_star, n_pairs, chain.reduced_mass());
    let wlo = (wto * wto + s / eps_inf).sqrt();
    PolarCrystal {
        name: "model diatomic ionic chain",
        eps_inf,
        modes: vec![PhononMode::new(
            rad_to_cm(wto),
            rad_to_cm(wlo),
            gamma_cm,
            gamma_cm,
        )],
        eps0_measured: f64::NAN,
        source: "Born effective charge (non-analytic term); see polariton.rs",
    }
}

/// Cheap spatial-dispersion indicator: the local permittivity `eps(omega)`
/// assumes `q -> 0`. The chain's optical frequency departs from its zone-centre
/// value by `O((q a)^2)`, so `q a` is the switch.
pub fn spatial_dispersion_indicator(q: f64, a: f64) -> f64 {
    (q * a).abs()
}

/// Cheap film-thickness indicator: a film `d` thick is homogeneous only when it
/// spans many lattice spacings `a`; the discrete atomic layers (and the
/// thickness-dependent image modes) matter when `a/d` is not small.
pub fn film_thickness_indicator(a: f64, d: f64) -> f64 {
    a / d
}

/// Single Lorentz-oscillator permittivity
/// `eps_inf (omega_LO^2 - omega^2 - i gamma omega)/(omega_TO^2 - omega^2 - i gamma omega)`.
/// Identical to the factorised form with `gamma_LO = gamma_TO = gamma`.
pub fn lorentz_single(eps_inf: f64, omega_to: f64, omega_lo: f64, gamma: f64, omega: f64) -> C {
    let num = C::real(omega_lo * omega_lo - omega * omega) - C::i() * (gamma * omega);
    let den = C::real(omega_to * omega_to - omega * omega) - C::i() * (gamma * omega);
    C::real(eps_inf) * num / den
}
