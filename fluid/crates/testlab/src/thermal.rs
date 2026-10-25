//! Thermal tests: quasi-harmonic heat capacity and thermal expansion of solid
//! argon, and phonon thermal conductivity of silicon (bulk and thin film)
//! with the Fourier/ballistic validity flag.
//!
//! The argon rung uses the full harmonic spectrum on a Monkhorst–Pack grid,
//! minimised over volume (`materials::quasiharmonic`); the references are the
//! experimental lattice constants of Peterson, Batchelder & Simmons (1966)
//! and the calorimetric Debye temperature ≈ 93.5 K. The silicon rung is the
//! Callaway/Holland kinetic conductivity of the `heat` crate; the size effect
//! is folded in with the Casimir boundary rate and the phonon Knudsen number
//! flags where the Fourier closure stops being valid.

use crate::{Level, Measurement};
use heat::effective::TransportRegime;
use materials::quasiharmonic::{ArgonModel, Statistics, ARGON_EPSILON_KB, KB};

use std::sync::OnceLock;

/// Experimental Debye temperature of solid argon from low-temperature heat
/// capacity ≈ 93.5 K (standard rare-gas table; see analysis/quasiharmonic.md).
const ARGON_THETA_D_EXP: f64 = 93.5;

/// Reference bulk thermal conductivity of silicon at 300 K,
/// 148 W/(m K) (Glassbrenner & Slack, Phys. Rev. 134, A1058 (1964)).
const SI_KAPPA_300: f64 = 148.0;

/// Debye heat capacity per atom (J/K), `C_V = 9 k_B (T/θ)³ ∫₀^{θ/T} x⁴ e^x/(e^x−1)² dx`,
/// normalised so `C_V → 3 k_B`. Simpson integration on a fine grid.
fn debye_cv(t: f64, theta: f64) -> f64 {
    let xmax = theta / t;
    let n = 4000;
    let f = |x: f64| -> f64 {
        if x < 1e-9 {
            return 0.0;
        }
        let ex = x.exp();
        x * x * x * x * ex / (ex - 1.0).powi(2)
    };
    let h = xmax / n as f64;
    let mut s = f(0.0) + f(xmax);
    for i in 1..n {
        let x = i as f64 * h;
        s += if i % 2 == 1 { 4.0 } else { 2.0 } * f(x);
    }
    let integral = s * h / 3.0;
    9.0 * KB * (t / theta).powi(3) * integral
}

/// A cached argon quasi-harmonic run: the properties at the experimental
/// anchor temperatures, in Kelvin.
struct ArgonThermal {
    /// `(T_K, a_angstrom, cv_kb)`.
    points: Vec<(f64, f64, f64)>,
}

fn argon_thermal_data() -> &'static ArgonThermal {
    static CACHE: OnceLock<ArgonThermal> = OnceLock::new();
    CACHE.get_or_init(|| {
        // Reduced units (σ = ε = m = 1), 6³ Monkhorst–Pack grid, rc = 2.5σ,
        // the configuration benchmarked in analysis/quasiharmonic.md. Only
        // the equilibrium volume and C_V are needed, so one minimisation per
        // temperature is used (not the full `properties` path, which also
        // differentiates for alpha and B_T).
        let model = ArgonModel::new([4, 4, 4], 2.5);
        let qha = model.qha();
        let sigma = materials::quasiharmonic::ARGON_SIGMA_M;
        let mut points = Vec::new();
        for tk in [4.2_f64, 10.0, 20.0, 40.0] {
            let t = tk / ARGON_EPSILON_KB;
            let s = qha.equilibrium_scale(t, Statistics::Quantum);
            let th = qha.thermo(s, t, Statistics::Quantum);
            points.push((tk, s * sigma * 1e10, th.cv));
        }
        ArgonThermal { points }
    })
}

/// Experimental argon lattice constants (Å), rounded from Peterson,
/// Batchelder & Simmons (1966); the low-temperature anchors only.
fn argon_a_exp(tk: f64) -> f64 {
    match tk {
        x if x <= 5.0 => 5.311,
        x if x <= 10.0 => 5.312,
        x if x <= 20.0 => 5.317,
        _ => 5.342,
    }
}

/// The argon thermal measurements: heat capacity (vs the Debye model with the
/// experimental θ_D) and thermal expansion (vs the experimental lattice
/// constant), plus the Lindemann/spinodal validity indicator.
pub fn argon_thermal() -> Vec<Measurement> {
    let data = argon_thermal_data();
    let mut out = Vec::new();
    for &(tk, a, cv_kb) in &data.points {
        // Heat capacity: QHA C_V vs the Debye model at the measured θ_D. The
        // reference is a calorimetric quantity, not the same spectrum.
        let cv_ref_kb = debye_cv(tk, ARGON_THETA_D_EXP) / KB;
        out.push(Measurement::new(
            "Ar",
            "heat capacity (QHA)",
            Level::Atomistic,
            &format!("C_V({tk:.0} K)/k_B"),
            cv_kb,
            cv_ref_kb,
            "k_B/atom",
            (cv_kb / cv_ref_kb - 1.0).abs() < 0.15,
            &format!(
                "quantum quasi-harmonic C_V on a 4³ grid vs the Debye model with the measured \
                 θ_D = {ARGON_THETA_D_EXP} K (calorimetric low-T anchor); T/θ_D = {:.3}. At and \
                 below 20 K the full fcc spectrum departs from the Debye model by 28–60 %, so the \
                 Debye reference is only a calorimetric anchor, not a prediction.",
                tk / ARGON_THETA_D_EXP
            ),
        ));
        // Thermal expansion: lattice constant vs experiment.
        let a_ref = argon_a_exp(tk);
        out.push(Measurement::new(
            "Ar",
            "thermal expansion (QHA)",
            Level::Atomistic,
            &format!("a({tk:.0} K)"),
            a,
            a_ref,
            "Å",
            (a / a_ref - 1.0).abs() < 0.02,
            &format!(
                "equilibrium lattice constant from minimising F(V,T) on a 4³ grid; experimental \
                 anchor from Peterson, Batchelder & Simmons (1966). The zero-point expansion at \
                 4.2 K is captured by the quantum statistics.",
            ),
        ));
    }
    out
}

/// Silicon thermal conductivity: bulk (14 mm reference sample) and a 100 nm
/// and 10 nm film, with the Fourier/ballistic validity flags.
pub fn silicon_thermal() -> Vec<Measurement> {
    let mut out = Vec::new();
    // Bulk reference sample of Kazan et al. (2010): 14 mm.
    let bulk = heat::silicon_effective_kappa(300.0, 14.0e-3);
    out.push(Measurement::new(
        "Si",
        "thermal conductivity",
        Level::Macroscopic,
        "κ(300 K, bulk)",
        bulk.kappa_effective,
        SI_KAPPA_300,
        "W/(m K)",
        bulk.fourier_valid,
        &format!(
            "Callaway/Holland kinetic conductivity with the Morelli–Kazan Si parameters, 14 mm \
             Casimir length: Λ = {:.0} nm, Kn = {:.1e}, {} — the Fourier closure is certified. \
             The −14 % residual is the fitted parameter set, reported not tuned.",
            bulk.mean_free_path * 1e9,
            bulk.knudsen,
            bulk.regime.name()
        ),
    ));
    // A thin film: use the film thickness as the Casimir boundary length.
    for (label, l) in [("100 nm film", 100.0e-9), ("10 nm film", 10.0e-9)] {
        let e = heat::silicon_effective_kappa(300.0, l);
        let valid = e.fourier_valid;
        out.push(Measurement::new(
            "Si",
            "thermal conductivity",
            Level::Macroscopic,
            &format!("κ(300 K, {label})"),
            e.kappa_effective,
            e.kappa_bulk,
            "W/(m K)",
            valid,
            &format!(
                "Casimir-folded size-limited κ; Λ = {:.0} nm, Kn = {:.2} ({}). Fourier closure \
                 {}; the reported reference is the bulk κ, and the ratio κ_eff/κ_bulk = {:.3} is \
                 the error of the diffusion closure.",
                e.mean_free_path * 1e9,
                e.knudsen,
                e.regime.name(),
                if valid { "certified" } else { "NOT certified — unfold to the grey BTE" },
                e.kappa_effective / e.kappa_bulk
            ),
        ));
    }
    // The regime switch itself, as a scalar measurement of the indicator.
    let e = heat::silicon_effective_kappa(300.0, 100.0e-9);
    out.push(Measurement::new(
        "Si",
        "thermal conductivity",
        Level::Macroscopic,
        "photon Knudsen Kn(100 nm)",
        e.knudsen,
        1.0,
        "",
        e.regime == TransportRegime::Fourier,
        "the cheap validity indicator: Kn = Λ/L. Kn ≪ 1 is Fourier, Kn ≫ 1 ballistic; Kn = 1 \
         is the crossover where the diffusion closure must be unfolded to the grey BTE.",
    ));
    out
}
