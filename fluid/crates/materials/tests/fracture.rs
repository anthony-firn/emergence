//! Fracture benchmarks: the Griffith fold, the atomistic crack and
//! lattice trapping, the indicators and the cohesive-zone closure.
//!
//! References are collected in `analysis/fracture.md`. Tolerances are set
//! by the physics (machine precision for the Stroh algebra, a few percent
//! for the finite atomistic crack), never loosened to force agreement.

use materials::fracture::*;
use materials::lattice::Lattice;
use materials::potential::{LennardJones, PairPotential};

fn lj() -> LennardJones {
    LennardJones::new(1.0, 1.0)
}

/// 2D triangular LJ properties: `(a0, E_2D, nu_2D, gamma_2D)`.
fn props() -> (f64, f64, f64, f64) {
    let pot = lj();
    let a0 = lj_2d_equilibrium();
    let (e2d, nu2d) = triangular_2d_moduli(&pot, a0, 1.3 * a0);
    let h = 3.0_f64.sqrt() / 2.0;
    let lat = Lattice::custom(
        [[a0, 0.0, 0.0], [0.5 * a0, h * a0, 0.0], [0.0, 0.0, 1.0]],
        vec![[0.0, 0.0, 0.0]],
    );
    let gamma = surface_energy_broken_bond(&lat, &pot, 1.3 * a0, [0.0, 1.0, 0.0]);
    (a0, e2d, nu2d, gamma)
}

fn cubic(c11: f64, c12: f64, c44: f64) -> [[f64; 6]; 6] {
    let mut c = [[0.0; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = if i == j { c11 } else { c12 };
        }
    }
    for i in 3..6 {
        c[i][i] = c44;
    }
    c
}

// ---------------------------------------------------------------------------
// 1. Griffith fold: the anisotropic Stroh/Sih-Paris-Irwin factor
// ---------------------------------------------------------------------------

#[test]
fn stroh_irwin_matches_sih_paris_irwin_and_antiplane() {
    // For a crack on a cubic symmetry plane the general Stroh Irwin matrix
    // must reproduce (a) the antiplane result H33 = 1/C44 and (b) the
    // Sih-Paris-Irwin orthotropic energy-release factor computed from the
    // plane-strain reduced compliances. This is the independent check.
    for c44 in [2.0_f64, 1.3, 0.5] {
        let c = cubic(3.0, 1.0, c44);
        let h = stroh_irwin_matrix(&c);
        assert!((h[2][2] - 1.0 / c44).abs() < 1e-9, "antiplane H33");
        let s = materials::anisotropic::compliance(&c);
        let s11 = s[0][0] - s[0][2] * s[0][2] / s[2][2];
        let s12 = s[0][1] - s[0][2] * s[1][2] / s[2][2];
        let s22 = s[1][1] - s[1][2] * s[1][2] / s[2][2];
        let s66 = s[5][5];
        let spi = (s11 * s22 / 2.0).sqrt()
            * ((s22 / s11).sqrt() + (2.0 * s12 + s66) / (2.0 * s11)).sqrt();
        let rel = (0.5 * h[1][1] - spi).abs() / spi;
        assert!(rel < 1e-12, "c44={c44}: G/K^2 Stroh {} vs SPI {spi}", 0.5 * h[1][1]);
    }
}

#[test]
fn isotropic_griffith_limits_are_exact() {
    // Plane strain K_IC = sqrt(2 gamma E/(1-nu^2)); plane stress
    // sqrt(2 gamma E); 2D sqrt(2 gamma E2D). Also the Sih-Paris-Irwin
    // formula must collapse to 1/E' for isotropic compliances.
    let (gamma, e, nu) = (1.0, 200.0, 0.3);
    let kp = griffith_plane_strain(gamma, e, nu);
    let expected = (2.0 * gamma * e / (1.0 - nu * nu)).sqrt();
    assert!((kp - expected).abs() < 1e-12);
    assert!((griffith_plane_stress(gamma, e) - (2.0 * gamma * e).sqrt()).abs() < 1e-12);
    assert!((griffith_2d(gamma, e) - (2.0 * gamma * e).sqrt()).abs() < 1e-12);
    // isotropic reduced compliances
    let c = cubic(3.0, 1.0, 1.0); // lambda = mu = 1
    let s = materials::anisotropic::compliance(&c);
    let s11 = s[0][0] - s[0][2] * s[0][2] / s[2][2];
    let s12 = s[0][1] - s[0][2] * s[1][2] / s[2][2];
    let s66 = s[5][5];
    let gk2 = (s11 * s11 / 2.0).sqrt() * (1.0 + (2.0 * s12 + s66) / (2.0 * s11)).sqrt();
    let eprime = 2.5 / (1.0 - 0.25 * 0.25);
    assert!((gk2 - 1.0 / eprime).abs() < 1e-12, "{gk2} vs {}", 1.0 / eprime);
}

// ---------------------------------------------------------------------------
// 2. Atomistic crack: the K-field and the energy release rate
// ---------------------------------------------------------------------------

#[test]
fn center_crack_field_reduces_to_the_near_tip_k_field() {
    // For r/a << 1 the exact finite-crack Westergaard field must approach
    // the universal mode-I K-field with the same K.
    let (_, e, nu, _) = props();
    let a = 100.0;
    let k = 1.0;
    let mut max_rel = 0.0_f64;
    for (x, y) in [(0.05, 0.02), (0.1, 0.0), (0.1, 0.1), (0.2, -0.1)] {
        let c = center_crack_field(x, y, a, k, e, nu);
        let m = mode_i_field(x, y, k, e, nu);
        let sc = m[0].abs().max(m[1].abs()).max(1e-12);
        max_rel = max_rel.max((c[0] - m[0]).abs().max((c[1] - m[1]).abs()) / sc);
    }
    // error is O(r/a) ~ 2e-3; require 1e-2.
    assert!(max_rel < 1e-2, "near-tip relative error {max_rel}");
}

#[test]
fn two_d_surface_energy_and_griffith_prediction() {
    // For the NN-only triangular LJ sheet the broken-bond surface energy is
    // exactly |phi(a0)|/a0 = 1/a0 (a0 = 2^{1/6} sigma), matching the direct
    // crack-formation energy per unit advance.
    let pot = lj();
    let (a0, e2d, nu2d, gamma) = props();
    assert!((gamma - 1.0 / a0).abs() < 1e-12, "broken-bond gamma {gamma}");
    assert!((nu2d - 1.0 / 3.0).abs() < 1e-9, "2D Poisson ratio {nu2d}");
    let kg = griffith_2d(gamma, e2d);
    assert!(kg > 0.0 && kg.is_finite());
    // direct formation energy: E(cs+1)-E(cs) per unit advance = 2 gamma a0
    let ebulk = 3.0 * pot.energy(a0);
    let mut diff = 0.0_f64;
    for cs in [6_i32, 7] {
        let sheet = build_cracked_sheet(20, 20, a0, 1.3 * a0, cs, &pot);
        let u = vec![[0.0; 3]; sheet.atoms.len()];
        let (e, _) = materials::defects::energy_forces(&sheet.bonds, &pot, &u);
        let excess = e - sheet.atoms.len() as f64 * ebulk;
        if cs == 6 {
            diff -= excess;
        } else {
            diff += excess;
        }
    }
    let gamma_direct = diff / (2.0 * a0);
    assert!((gamma_direct - gamma).abs() < 1e-10, "gamma direct {gamma_direct} vs {gamma}");
}

#[test]
fn atomistic_crack_advances_above_griffith() {
    // Virtual crack extension at fixed grip: G(K) is negative (stable)
    // below the Griffith load and positive (advancing) above. The crossing
    // is the atomistic critical K; lattice trapping puts it above K_G.
    let pot = lj();
    let (a0, e2d, nu2d, gamma) = props();
    let kg = griffith_2d(gamma, e2d);
    // A crack of 10 sites in a 50x50-site sheet: the cohesive zone is a few
    // sites, so this is the smallest sample that is K-dominated.
    let (nx, ny, cs) = (25, 25, 10);
    // Locate the crossing of G(K) through zero on a bracket that stays in
    // the single-advance regime (K <= ~1.8 K_G).
    let kc = critical_k_vce(
        &pot, nx, ny, a0, 1.3 * a0, cs, cs, e2d, nu2d, 1.0 * kg, 1.8 * kg, 0.3, 2500, 1e-5,
    );
    assert!(kc > kg, "K_c = {kc} should exceed K_G = {kg}");
    assert!(kc < 1.8 * kg, "K_c = {kc} out of range");
}

// ---------------------------------------------------------------------------
// 3. Indicators and the cohesive-zone closure
// ---------------------------------------------------------------------------

#[test]
fn process_zone_and_lefm_indicator() {
    // r_p = K^2/(2 pi sigma_y^2); the LEFM flag fires when r_p < threshold L.
    let k = 10.0;
    let sy = 5.0;
    let rp = process_zone_radius(k, sy);
    assert!((rp - 100.0 / (2.0 * std::f64::consts::PI * 25.0)).abs() < 1e-12);
    assert!(lefm_valid(k, sy, 100.0, 0.1));
    assert!(!lefm_valid(k, sy, 0.001, 0.1));
}

#[test]
fn rice_brittle_ductile_threshold_and_ratio() {
    // Rice (1992): cleavage G = 2 gamma_s competes with emission
    // G = 8 gamma_us; the threshold is gamma_us/gamma_s = 1/4.
    let s = 1.0;
    let brittle = rice_ductile_brittle(0.5, s); // ratio 0.5 > 0.25
    assert!(brittle.brittle);
    let ductile = rice_ductile_brittle(0.1, s); // ratio 0.1 < 0.25
    assert!(!ductile.brittle);
    assert!((brittle.g_cleave - 2.0).abs() < 1e-12);
    assert!((brittle.g_emit - 4.0).abs() < 1e-12);
    assert!((brittle.ratio - 0.5).abs() < 1e-12);
}

#[test]
fn universal_binding_work_of_separation_is_two_gamma() {
    // Integrate the UBER traction and check the work equals 2 gamma_s and
    // the peak traction is (2 gamma_s/lambda) e^{-1} at delta = lambda.
    let gamma = 0.9;
    let lambda = 0.5;
    let cz = cohesive_zone_universal(gamma, lambda);
    assert!((cz.work_of_separation - 2.0 * gamma).abs() < 1e-12);
    let e0 = 2.0 * gamma;
    let n = 20_000;
    let dmax = 20.0 * lambda;
    let mut tractions = Vec::with_capacity(n + 1);
    let mut deltas = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let d = dmax * i as f64 / n as f64;
        deltas.push(d);
        tractions.push(traction_from_energy(|x| universal_binding(x, e0, lambda), d, 1e-6));
    }
    let work = cohesive_work(&tractions, &deltas);
    assert!((work - 2.0 * gamma).abs() < 1e-6, "work {work}");
    let tmax = tractions.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert!((tmax - cz.t_max).abs() < 1e-6, "tmax {tmax} vs {}", cz.t_max);
    assert!((cz.t_max - 2.0 * gamma / lambda * (-1.0_f64).exp()).abs() < 1e-12);
}

#[test]
fn local_tip_bond_rupture_is_a_precursor_not_trapping() {
    // The original `lattice_trapping` routine reports the load window over
    // which a tip bond exceeds the LJ spinodal. This is a *local bond
    // instability* (a precursor), not the load-controlled crack-advance
    // bound: it lies below the thermodynamic advance threshold. The test
    // pins that down by comparing with the VCE threshold, and deliberately
    // does NOT claim K_- < K_G < K_+.
    let pot = lj();
    let (a0, e2d, nu2d, gamma) = props();
    let kg = griffith_2d(gamma, e2d);
    let sheet = build_cracked_sheet(13, 13, a0, 1.3 * a0, 5, &pot);
    let r_break = (26.0_f64 / 7.0).powf(1.0 / 6.0);
    let res = lattice_trapping(&sheet, &pot, e2d, nu2d, gamma, 1.5 * kg, 14, r_break, 2);
    println!(
        "precursor bond-rupture window: K+/KG={:.4} K-/KG={:.4}",
        res.k_plus / kg,
        res.k_minus / kg
    );
    assert!(res.k_plus.is_finite() && res.k_minus > 0.0);
    assert!(res.k_minus < res.k_plus, "K_- = {}, K_+ = {}", res.k_minus, res.k_plus);
    // Precursor window is below the symmetric Griffith load.
    assert!(res.k_plus < kg, "precursor K_+ = {} should be below K_G", res.k_plus);
}

#[test]
fn relaxed_surface_energy_is_linear_in_crack_length() {
    // Relaxing the pristine sheet at zero load gives a surface energy that
    // is *exactly* proportional to the removed-row length: successive
    // increments equal 2 gamma_s a0 to machine precision. Hence the
    // relaxed energies of successive crack positions carry no lattice-
    // periodic signal; any trapping barrier is dynamic (a Peierls barrier
    // of the tip), not a variation of these minima.
    let pot = lj();
    let (a0, e2d, nu2d, gamma) = props();
    let ebulk = 3.0 * pot.energy(a0);
    let mut prev: Option<f64> = None;
    for cs in 4..=10 {
        let sheet = build_cracked_sheet(30, 30, a0, 1.3 * a0, cs, &pot);
        let u0 = k_field_displacements(&sheet, 0.0, e2d, nu2d);
        let load = vec![[0.0; 3]; sheet.atoms.len()];
        let (_, e, _) = minimize_atoms(&sheet, &pot, u0, &load, 6000, 1e-8);
        let esurf = e - sheet.atoms.len() as f64 * ebulk;
        if let Some(p) = prev {
            let inc = esurf - p;
            assert!(
                (inc - 2.0 * gamma * a0).abs() < 1e-8,
                "surface-energy increment {inc} vs {}",
                2.0 * gamma * a0
            );
        }
        prev = Some(esurf);
    }
}

#[test]
fn effective_tip_k_at_vce_crossing_is_griffith() {
    // The VCE threshold is where the net energy release rate is zero; the
    // elastic driving force is then exactly 2 gamma_s and the effective
    // tip K is K_G. Check the helper agrees.
    let (_, e2d, _, gamma) = props();
    let kg = griffith_2d(gamma, e2d);
    let keff = effective_tip_k(0.0, gamma, e2d);
    assert!((keff - kg).abs() < 1e-12);
    // A negative net G (stable crack) implies a smaller effective tip K.
    assert!(effective_tip_k(-0.5, gamma, e2d) < kg);
}

#[test]
#[ignore]
fn experiment_size_sweep() {
    let pot = lj();
    let (a0, e2d, nu2d, gamma) = props();
    let kg = griffith_2d(gamma, e2d);
    let cs = 10;
    // Several sample radii with the crack length held fixed at 10 sites.
    let sizes: Vec<i32> = std::env::var("SWEEP_NX")
        .ok()
        .and_then(|s| s.parse().ok())
        .map(|n| vec![n])
        .unwrap_or_else(|| vec![18, 25, 37]);
    // Independent relaxations per K (a warm-started curve can hop basins).
    let ks: Vec<f64> = (0..=5).map(|i| (1.10 + 0.125 * i as f64) * kg).collect();
    let mut table = Vec::new();
    for nx in sizes {
        let mut gs = Vec::new();
        for &k in &ks {
            gs.push(energy_release_rate_vce(
                &pot, nx, nx, a0, 1.3 * a0, cs, cs, k, e2d, nu2d, 4000, 1e-7,
            ));
        }
        let mut kc = f64::NAN;
        for i in 1..ks.len() {
            if gs[i - 1] <= 0.0 && gs[i] >= 0.0 {
                let w = -gs[i - 1] / (gs[i] - gs[i - 1]);
                kc = ks[i - 1] + w * (ks[i] - ks[i - 1]);
            }
        }
        let r_over_a = nx as f64 / (0.5 * cs as f64);
        println!("nx={nx} cs={cs} R/a={r_over_a:.1}: Kc/KG={:.4}", kc / kg);
        for (k, g) in ks.iter().zip(&gs) {
            println!("   K/KG={:.2}  G={:+.5}", k / kg, g);
        }
        table.push((r_over_a, kc / kg));
    }
    // Extrapolate Kc/KG = A + B/(R/a) by least squares over the three
    // points. The two largest already agree to 1e-3, so the fit is bounded.
    let n = table.len() as f64;
    let sx: f64 = table.iter().map(|t| 1.0 / t.0).sum();
    let sy: f64 = table.iter().map(|t| t.1).sum();
    let sxx: f64 = table.iter().map(|t| 1.0 / (t.0 * t.0)).sum();
    let sxy: f64 = table.iter().map(|t| t.1 / t.0).sum();
    let slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    let inter = (sy - slope * sx) / n;
    let plateau = table.last().unwrap().1;
    println!(
        "extrapolated R->inf: 1/R fit Kc/KG = {inter:.4} (slope {slope:.4}); \
         largest-R plateau Kc/KG = {plateau:.4}"
    );
    assert!(inter > 1.0 && inter < 1.6, "extrapolated Kc/KG {inter}");
    assert!(plateau > 1.0 && plateau < 1.6, "plateau Kc/KG {plateau}");
}

#[test]
fn cubic_lj_toughness_anisotropy() {
    // LJ fcc infinite-cutoff Cauchy-Born constants from analysis/materials.md
    // (C12 = C44 by the central-force Cauchy relation): C11 = 105.182,
    // C12 = C44 = 60.187 eps/sigma^3. The Zener ratio is 2.68, so the
    // {100} crack toughness is 5.6% below the Hill-isotropic estimate.
    let c = cubic(105.182, 60.187, 60.187);
    let zener = 2.0 * 60.187 / (105.182 - 60.187);
    assert!((zener - 2.6753_f64).abs() < 1e-3, "Zener {zener}");
    let iso = materials::polycrystal::hill(&c);
    let (e, nu) = (iso.young(), iso.poisson());
    let g100 = 3.19;
    let k_iso = griffith_plane_strain(g100, e, nu);
    assert!((k_iso - 26.655).abs() < 0.05, "isotropic K_IC {k_iso}");
    let h = stroh_irwin_matrix(&c);
    let k_aniso = (4.0 * g100 / h[1][1]).sqrt();
    assert!((k_aniso - 25.167).abs() < 0.05, "anisotropic K_IC {k_aniso}");
    let ratio = k_aniso / k_iso;
    assert!((ratio - 0.9442).abs() < 2e-3, "anisotropy ratio {ratio}");
    // {111}[1-10] crack: same construction, coupled mode I/II.
    let n = [1.0 / 3.0_f64.sqrt(); 3];
    let inv2 = 1.0 / 2.0_f64.sqrt();
    let p = [inv2, -inv2, 0.0];
    let t = [
        n[1] * p[2] - n[2] * p[1],
        n[2] * p[0] - n[0] * p[2],
        n[0] * p[1] - n[1] * p[0],
    ];
    let c111 = materials::anisotropic::rotate_stiffness(&c, &[p, n, t]);
    let h111 = stroh_irwin_matrix(&c111);
    let k111 = (4.0 * 3.09 / h111[1][1]).sqrt();
    assert!(k111 > 20.0 && k111 < 35.0, "K_IC(111) {k111}");
}

