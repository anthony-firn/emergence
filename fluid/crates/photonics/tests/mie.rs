//! Mie scattering integration tests: published benchmarks, the Rayleigh and
//! anomalous-diffraction folds with their measured convergence, and the three
//! repo applications (cloud radiative properties, optical forces, virion
//! detection). Every published number is cited in `analysis/mie.md`.
//!
//! Run with `cargo test -p photonics --test mie -- --nocapture` to see the
//! convergence and application tables.

use photonics::mie::*;
use photonics::{Dispersion, C};
use std::f64::consts::PI;

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a} vs {b} (tol {tol})");
}

/// Bohren & Huffman 1983 Appendix / Wiscombe MIEV0 published efficiencies,
/// reached through the public `MieScattering` API.
#[test]
fn published_wiscombe_and_bohren_huffman() {
    // m = 1.55, x = 2 pi * 0.525 / 0.6328 (Bohren & Huffman dielectric test).
    let x = 2.0 * PI * 0.525 / 0.6328;
    let e = MieScattering::new(x, C::new(1.55, 0.0)).efficiencies();
    close(e.q_ext, 3.10543, 1e-5);
    close(e.q_sca, 3.10543, 1e-5);
    close(e.q_back, 2.92534, 1e-5);
    close(e.g, 0.63314, 1e-5);

    // m = 1.55 - 0.1i (absorbing): Q_ext, Q_sca, Q_back, g.
    let e = MieScattering::new(x, C::new(1.55, 0.1)).efficiencies();
    close(e.q_ext, 2.86165188243, 1e-7);
    close(e.q_sca, 1.66424911991, 1e-7);
    close(e.q_back, 0.20599534080, 1e-7);
    close(e.g, 0.80128972639, 1e-7);

    // MIEV0 non-absorbing m = 0.75 at x = 10 and 1000.
    let e = MieScattering::new(10.0, C::new(0.75, 0.0)).efficiencies();
    close(e.q_sca, 2.232265, 1e-6);
    close(e.g, 0.896473, 1e-6);
    let e = MieScattering::new(1000.0, C::new(0.75, 0.0)).efficiencies();
    close(e.q_sca, 1.997908, 1e-6);
    close(e.g, 0.844944, 1e-6);

    // MIEV0 absorbing m = 1.5 - 1i, x = 100 and 10000. The last case spans
    // ~10^4 terms and exercises the downward D recurrence at scale.
    let e = MieScattering::new(100.0, C::new(1.5, 1.0)).efficiencies();
    close(e.q_ext, 2.097502, 1e-6);
    close(e.q_sca, 1.283697, 1e-3);
    close(e.g, 0.850252, 1e-3);
    let e = MieScattering::new(10000.0, C::new(1.5, 1.0)).efficiencies();
    close(e.q_ext, 2.004368, 1e-6);
    close(e.q_sca, 1.236575, 1e-6);
    close(e.g, 0.846309, 1e-6);
}

/// Rayleigh fold: relative error versus `x` for two contrasts, showing the O(x^2)
/// rate. `Q_sca` is checked against the exact series.
#[test]
fn fold_convergence_table() {
    println!("\nRayleigh fold, Q_sca relative error (exact vs (8/3)x^4|K|^2):");
    println!(
        "{:>9} {:>14} {:>14} {:>14}",
        "x", "err m=1.05", "err m=1.33", "err m=1.5"
    );
    for &x in &[0.5, 0.25, 0.125, 0.0625, 0.03125] {
        let mut row = [0.0f64; 3];
        for (i, &mr) in [1.05f64, 1.33, 1.5].iter().enumerate() {
            let m = C::new(mr, 0.0);
            let exact = efficiencies(x, m).q_sca;
            let approx = rayleigh_qsca(x, m);
            row[i] = (exact - approx).abs() / exact;
        }
        println!(
            "{x:>9.5} {:>14.3e} {:>14.3e} {:>14.3e}",
            row[0], row[1], row[2]
        );
    }

    println!("\nAnomalous-diffraction fold, Q_ext absolute error (exact vs van de Hulst):");
    println!("{:>9} {:>14} {:>14}", "x", "rho(m=1.05)", "|err|");
    for &x in &[10.0, 20.0, 50.0, 100.0, 200.0, 400.0] {
        let m = 1.05;
        let exact = efficiencies(x, C::new(m, 0.0)).q_ext;
        let adt = anomalous_diffraction_qext_real(x, m);
        println!(
            "{x:>9.1} {:>14.3} {:>14.3e}",
            2.0 * x * (m - 1.0),
            (exact - adt).abs()
        );
    }

    // Anomalous-diffraction fold at fixed phase delay rho = 2x(m-1): x grows,
    // m -> 1, and the error should fall as O(1/x).
    println!("\nAnomalous-diffraction fold at fixed rho = 10 (m -> 1 as x grows):");
    println!("{:>9} {:>14} {:>14}", "x", "m", "|err Q_ext|");
    for &x in &[10.0, 20.0, 40.0, 80.0, 160.0, 320.0] {
        let m = 1.0 + 10.0 / (2.0 * x);
        let exact = efficiencies(x, C::new(m, 0.0)).q_ext;
        let adt = anomalous_diffraction_qext_real(x, m);
        println!("{x:>9.1} {:>14.5} {:>14.3e}", m, (exact - adt).abs());
    }

    println!("\nExtinction paradox: Q_ext -> 2 as x grows (m = 1.5):");
    for &x in &[10.0, 100.0, 1000.0, 10000.0, 100000.0] {
        let e = efficiencies(x, C::new(1.5, 0.0));
        println!(
            "{x:>9.0} Q_ext = {:.8}  (Q_ext-2) = {:.3e}",
            e.q_ext,
            e.q_ext - 2.0
        );
    }
    println!();
}

/// Print the published-vs-computed table used in `analysis/mie.md`.
#[test]
fn benchmark_error_table() {
    // (m_re, m_im, x, published Q_sca, published g)
    #[rustfmt::skip]
    let cases: &[(f64, f64, f64, f64, f64)] = &[
        (1.55, 0.0, 5.212819668567135, 3.10543, 0.63314),
        (1.55, 0.1, 5.212819668567135, 1.66424911991, 0.80128972639),
        (0.75, 0.0, 0.099, 7.4e-6, 0.001448),
        (0.75, 0.0, 10.0, 2.232265, 0.896473),
        (0.75, 0.0, 1000.0, 1.997908, 0.844944),
        (1.33, 1e-5, 100.0, 2.096594, 0.868959),
        (1.33, 1e-5, 10000.0, 1.723857, 0.907840),
        (1.5, 1.0, 1.0, 0.6634538, 0.192136),
        (1.5, 1.0, 100.0, 1.283697, 0.850252),
        (1.5, 1.0, 10000.0, 1.236575, 0.846309),
        (10.0, 10.0, 1.0, 2.049405, -0.110664),
        (10.0, 10.0, 100.0, 1.836785, 0.556215),
        (10.0, 10.0, 10000.0, 1.795393, 0.548194),
    ];
    println!(
        "\n{:>5} {:>7} {:>8} {:>12} {:>12} {:>9} {:>10} {:>9}",
        "m_re", "m_im", "x", "Qsca pub", "Qsca comp", "err", "g pub", "err g"
    );
    for &(mr, mi, x, qs, gp) in cases {
        let e = MieScattering::new(x, C::new(mr, mi)).efficiencies();
        println!(
            "{mr:>5.2} {mi:>7.0e} {x:>8.3} {qs:>12.7} {:>12.7} {:>9.1e} {:>10.6} {:>9.1e}",
            e.q_sca,
            (e.q_sca - qs).abs(),
            gp,
            (e.g - gp).abs()
        );
    }
    println!();
}

/// Stratocumulus-like cloud of water droplets: the participating-medium inputs
/// `radiation` consumes, compared with the large-particle law.
#[test]
fn cloud_radiative_properties() {
    let dist = GammaSizeDistribution::stratocumulus();
    // Water index from the crate's Daimon-Masumura (2007) Sellmeier fit.
    let n_water = photonics::dispersion::Sellmeier::water().index(550e-9).re;
    let m = C::new(n_water, 0.0);
    let wavelength = 550e-9;
    let n_number = 100e6; // 100 cm^-3, a stratocumulus number density (estimate)
    let med = cloud_medium(&dist, m, wavelength, 1.0, n_number, 4000);

    // Mass and the large-particle (Q_ext -> 2) law sigma = 3 LWC / (2 rho r_eff).
    let rho = 1000.0;
    let lwc = n_number * rho * 4.0 / 3.0 * PI * dist.moment(3.0);
    let sigma_large = 3.0 * lwc / (2.0 * rho * dist.r_eff);
    println!("\nStratocumulus-like cloud (r_eff = 10 µm, v_eff = 0.1, N = 100 cm^-3):");
    println!("  LWC               = {:.4} g/m^3", lwc * 1e3);
    println!(
        "  sigma_ext         = {:.5} m^-1 = {:.2} km^-1",
        med.sigma_ext,
        med.sigma_ext * 1e3
    );
    println!("  sigma_sca         = {:.5} m^-1", med.sigma_sca);
    println!("  sigma_abs         = {:.3e} m^-1", med.sigma_abs);
    println!("  single-scat albedo= {:.7}", med.albedo());
    println!("  g                 = {:.4}", med.g);
    println!(
        "  3 LWC/(2 rho r_e) = {:.5} m^-1 (Q_ext -> 2 law)",
        sigma_large
    );
    println!("  vertical tau over 300 m = {:.1}", med.sigma_ext * 300.0);

    // Published/typical stratocumulus values: sigma_ext ~ 30-60 km^-1,
    // single-scattering albedo > 0.999, g ~ 0.85.
    assert!(
        med.sigma_ext > 0.03 && med.sigma_ext < 0.06,
        "sigma_ext {}",
        med.sigma_ext
    );
    assert!(med.albedo() > 0.999);
    assert!((med.g - 0.85).abs() < 0.02, "g {}", med.g);
    assert!((med.sigma_ext - sigma_large).abs() / sigma_large < 0.05);
}

/// 1 µm silica bead in a 1064 nm Gaussian trap: exact Mie `Q_pr`, the Rayleigh
/// gradient force, and the depth in units of `k_B T`.
#[test]
fn optical_force_on_silica_bead() {
    let radius = 0.5e-6;
    // Fused silica (Malitson 1965) and water (Daimon-Masumura 2007) from the crate.
    let n_silica = photonics::dispersion::Sellmeier::fused_silica()
        .index(1064e-9)
        .re;
    let n_medium = photonics::dispersion::Sellmeier::water().index(1064e-9).re;
    let m = C::new(n_silica / n_medium, 0.0);
    let power = 0.1; // 100 mW
    let w0 = 0.5e-6;
    let f = optical_force(radius, 1064e-9, n_medium, m, power, w0, 300.0);

    println!("\nOptical trap: 1 µm silica bead in water, 1064 nm, {power} W, w0 = 0.5 µm:");
    println!("  Q_pr              = {:.4}", f.q_pr);
    println!(
        "  scattering force  = {:.3e} N = {:.3} pN",
        f.scattering_force,
        f.scattering_force * 1e12
    );
    println!(
        "  gradient force    = {:.3e} N = {:.3} pN",
        f.gradient_force,
        f.gradient_force * 1e12
    );
    println!(
        "  stiffness (Rayl.) = {:.3e} N/m = {:.3} pN/nm",
        f.stiffness,
        f.stiffness * 1e12 / 1e9
    );
    println!("  depth / k_B T     = {:.0}", f.depth_over_kbt);
    println!();

    assert!(f.q_pr > 0.01 && f.q_pr < 0.5, "q_pr {}", f.q_pr);
    // The trap is deep: thousands of k_B T.
    assert!(f.depth_over_kbt > 100.0, "U/kT {}", f.depth_over_kbt);
    // Exact Mie gives a radiation-pressure force of tens of pN at 100 mW.
    assert!(f.scattering_force > 1e-12 && f.scattering_force < 1e-9);
}

/// Virion-sized particle: cross-section versus size for detection in water.
#[test]
fn virion_detection_cross_section() {
    // Protein/lipid virion in water; relative index taken as ~1.1 (estimate).
    let n_medium = photonics::dispersion::Sellmeier::water().index(550e-9).re;
    let wavelength = 550e-9;
    let m_rel = 1.1;
    println!("\nVirion in water (m_rel = {m_rel}, lambda = 550 nm):");
    println!(
        "{:>8} {:>10} {:>12} {:>12}",
        "d (nm)", "x", "Q_sca", "C_sca (m^2)"
    );
    for &d_nm in &[20.0, 50.0, 100.0, 200.0, 500.0] {
        let radius = d_nm * 0.5e-9;
        let s =
            MieScattering::from_radius(radius, wavelength, C::new(m_rel * n_medium, 0.0), n_medium);
        let e = s.efficiencies();
        let c_sca = e.q_sca * PI * radius * radius;
        println!(
            "{d_nm:>8.0} {:>10.4} {:>12.5e} {:>12.5e}",
            s.x, e.q_sca, c_sca
        );
        assert!(c_sca.is_finite() && c_sca >= 0.0);
    }
    println!();

    // A 100 nm virion is a weak Rayleigh-like scatterer: cross-section well below
    // its geometric area pi a^2 = 7.85e-15 m^2.
    let radius = 50e-9;
    let s = MieScattering::from_radius(radius, wavelength, C::new(m_rel * n_medium, 0.0), n_medium);
    let e = s.efficiencies();
    let c_sca = e.q_sca * PI * radius * radius;
    let geometric = PI * radius * radius;
    assert!(
        c_sca < 0.1 * geometric,
        "C_sca/geom = {}",
        c_sca / geometric
    );
    // Rayleigh scaling: C_sca ~ d^6 (ratio over a factor-2 size change is 64).
    let c2 = MieScattering::from_radius(25e-9, wavelength, C::new(m_rel * n_medium, 0.0), n_medium)
        .efficiencies()
        .q_sca
        * PI
        * (25e-9f64).powi(2);
    let ratio = c_sca / c2;
    assert!(ratio > 40.0 && ratio < 80.0, "C(100)/C(50) = {ratio}");
}
