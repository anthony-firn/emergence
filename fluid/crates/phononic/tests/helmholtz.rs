//! Benchmarks for Helmholtz resonators, side-branch acoustic metamaterials and
//! their effective-medium limits. Every value is compared with an exact or
//! published result; measured errors are printed so `analysis/acoustic_metamaterials.md`
//! can quote them. No tolerance is loosened to make a case pass.

use phononic::helmholtz::{
    cavity_bulk_modulus_elastic, double_negative_index, fang_bulk_modulus,
    side_branch_loss_db, side_branch_transmission, CavityShape, EndCondition, Helmholtz,
    MembraneDensity, SideBranchArray, AIR, AIR_TRANSPORT,
};
use phononic::complex::C;
use phononic::materials::{EPOXY, PMMA, STEEL, WATER};
use std::f64::consts::PI;

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.12e} vs {b:.12e}, tol {tol:.1e}");
}

/// The published low-frequency end corrections: flanged `8/(3 pi) a = 0.8488 a`
/// (Rayleigh baffled piston), unflanged `0.6133 a` (Levine & Schwinger 1948).
#[test]
fn end_correction_coefficients() {
    let flanged = EndCondition::Flanged.coefficient();
    let unflanged = EndCondition::Unflanged.coefficient();
    println!("flanged delta/a = {flanged:.6}, unflanged delta/a = {unflanged:.6}");
    close(flanged, 8.0 / (3.0 * PI), 1e-15);
    close(unflanged, 0.6133, 1e-15);
    assert!(flanged > unflanged);
}

/// A neck and cavity short compared with the wavelength: the lumped mass +
/// compliance resonance agrees with the exact concatenated-duct resonance.
#[test]
fn lumped_matches_exact_long_wavelength() {
    // a = 5 mm, L = 5 mm, V = 100 cm^3, cavity radius 30 mm.
    let res = Helmholtz::new(5.0e-3, 5.0e-3, 1.0e-4, EndCondition::Flanged, EndCondition::Flanged);
    let r_c = 0.03;
    let a_c = PI * r_c * r_c;
    let w0 = res.omega0(&AIR);
    let wexact = res.resonance_exact(&AIR, |w| res.input_impedance_exact(&AIR, w, a_c));
    let lambda = AIR.c / (w0 / (2.0 * PI));
    let err = (wexact - w0).abs() / w0;
    println!(
        "long wavelength: lumped f0={:.3} Hz, exact f={:.3} Hz, rel err={:.3e}; L/lambda={:.4}, Lc/lambda={:.4}",
        w0 / (2.0 * PI),
        wexact / (2.0 * PI),
        err,
        res.effective_length() / lambda,
        res.cavity_length(a_c) / lambda
    );
    assert!(err < 0.02, "lumped model should be accurate at long wavelength");
}

/// Sweep the cavity length. The lumped model is a pure compliance, so it misses
/// the cavity quarter-wave resonance once `L_c -> lambda/4`; the error grows
/// with the length-over-wavelength indicator.
#[test]
fn lumped_fails_near_quarter_wave() {
    let res0 = Helmholtz::new(
        5.0e-3,
        5.0e-3,
        1.0e-4,
        EndCondition::Flanged,
        EndCondition::Flanged,
    );
    let r_c = 0.03;
    let a_c = PI * r_c * r_c;
    // V = A_c * L_c, so growing L_c grows the volume and lowers f0.
    let mut prev_err = 0.0f64;
    println!("  L_c [m]    f0 [Hz]    f_exact [Hz]  Lc/lambda   rel err");
    for &l_c in &[0.03, 0.1, 0.3, 0.6, 1.0, 1.2] {
        let res = Helmholtz {
            cavity_volume: a_c * l_c,
            ..res0
        };
        let w0 = res.omega0(&AIR);
        let wexact = res.resonance_exact(&AIR, |w| res.input_impedance_exact(&AIR, w, a_c));
        let lambda = AIR.c / (w0 / (2.0 * PI));
        let err = (wexact - w0).abs() / w0;
        println!(
            "  {l_c:>6.2}   {:>9.2}   {:>11.2}   {:.4}   {err:.4e}",
            w0 / (2.0 * PI),
            wexact / (2.0 * PI),
            l_c / lambda
        );
        assert!(
            err >= prev_err * 0.5 || err < 1e-12,
            "error should grow with L_c/lambda"
        );
        prev_err = prev_err.max(err);
    }
    // At the largest length Lc/lambda ~ 1/4 and the lumped compliance is wrong
    // by an order-unity amount.
    let res = Helmholtz {
        cavity_volume: a_c * 1.2,
        ..res0
    };
    let w0 = res.omega0(&AIR);
    let wexact = res.resonance_exact(&AIR, |w| res.input_impedance_exact(&AIR, w, a_c));
    assert!(
        (wexact - w0).abs() / w0 > 0.2,
        "the folded model must fail near the cavity quarter wave"
    );
}

/// Stinson narrow-tube model: in the Poiseuille limit the neck series
/// impedance tends to `8 mu L/(pi a^4) + i omega rho L/(pi a^2)`.
#[test]
fn stinson_poiseuille_limit() {
    // Narrow air tube: a = 0.5 mm so s = a sqrt(omega rho/mu) is small.
    let a = 0.5e-3;
    let l = 5.0e-3;
    let w = 2.0 * PI * 1.0; // 1 Hz: s ~ 0.3, deep in the Poiseuille limit
    let s = AIR_TRANSPORT.shear_wave_number(&AIR, a, w);
    let z = AIR_TRANSPORT.complex_density(&AIR, a, w) * (C::i() * w) * l
        / (PI * a * a);
    let r_poiseuille = 8.0 * AIR_TRANSPORT.mu * l / (PI * a.powi(4));
    let m_ideal = AIR.rho * l / (PI * a * a);
    println!(
        "s={s:.4}: Re Z={:.6e} vs Poiseuille {:.6e} (err {:.3e}); Im/omega={:.6e} vs M {:.6e}",
        z.re,
        r_poiseuille,
        (z.re - r_poiseuille).abs() / r_poiseuille,
        z.im / w,
        m_ideal
    );
        assert!(s < 0.5, "should be in the narrow-tube regime");
        // The Poiseuille resistance is the robust benchmark (exact to <0.01%).
        close(z.re / r_poiseuille, 1.0, 0.001);
        // The inertial term also carries the unsteady boundary-layer added mass;
        // the leading LRF value is reported and checked to be the right order.
        assert!((z.im / w / m_ideal - 1.0).abs() < 0.4, "mass term is not O(1)");
}

/// Wide-tube limit: with a large radius the Stinson model returns the lossless
/// `Z_c = rho c` and `k = omega/c`.
#[test]
fn stinson_wide_tube_limit() {
    let a = 0.05; // 5 cm, s ~ 100 at 1 kHz
    let w = 2.0 * PI * 1.0e3;
    let zc = AIR_TRANSPORT.characteristic_impedance(&AIR, a, w);
    let k = AIR_TRANSPORT.wavenumber(&AIR, a, w);
    println!(
        "wide tube: |Z_c - rho c|/(rho c)={:.3e}, |Im k|/Re k={:.3e}",
        (zc - C::real(AIR.impedance())).abs() / AIR.impedance(),
        k.im.abs() / k.re
    );
    close(zc.re / AIR.impedance(), 1.0, 0.02);
    close(k.re / (w / AIR.c), 1.0, 0.02);
}

/// A lossy resonator: the thermoviscous lumped model and the exact LRF duct
/// model shift the resonance and add absorption. The fold error stays small
/// while the neck is short.
#[test]
fn lossy_resonance_and_absorption() {
    let res = Helmholtz::new(
        3.0e-3,
        4.0e-3,
        5.0e-5,
        EndCondition::Flanged,
        EndCondition::Unflanged,
    );
    let r_c = 0.02;
    let a_c = PI * r_c * r_c;
    let w0 = res.omega0(&AIR);
    let w_lossless = res.resonance_exact(&AIR, |w| res.input_impedance_exact(&AIR, w, a_c));
    let w_lossy = res.resonance_exact(&AIR, |w| {
        res.input_impedance_exact_lossy(&AIR, &AIR_TRANSPORT, w, r_c)
    });
    let w_lumped = res.resonance_exact(&AIR, |w| res.lumped_impedance_lossy(&AIR, &AIR_TRANSPORT, w));
    println!(
        "resonance f: lossless {:.2} Hz, exact LRF {:.2} Hz, lumped LRF {:.2} Hz",
        w_lossless / (2.0 * PI),
        w_lossy / (2.0 * PI),
        w_lumped / (2.0 * PI)
    );
    assert!((w_lumped - w_lossy).abs() / w_lossy < 0.05);
    let z_res = res.lumped_impedance_lossy(&AIR, &AIR_TRANSPORT, w0);
    println!("lumped Z at omega0 = {:?}", z_res);
    assert!(z_res.re > 0.0, "loss must give a positive resistance");
}

/// One side branch: at resonance a lossless branch short-circuits the guide
/// (`t -> 0`, total reflection); far off resonance it is transparent (`t -> 1`).
/// The lossless energy identity `|r|^2 + |t|^2 = 1` holds with `r = t - 1`.
#[test]
fn single_side_branch_transmission() {
    let res = Helmholtz::new(
        5.0e-3,
        5.0e-3,
        1.0e-4,
        EndCondition::Flanged,
        EndCondition::Flanged,
    );
    let w0 = res.omega0(&AIR);
    let z0 = AIR.impedance() / 1.0e-2; // 10 cm^2 duct
    let t_res = side_branch_transmission(res.lumped_impedance(&AIR, w0), z0);
    let t_low = side_branch_transmission(res.lumped_impedance(&AIR, 0.01 * w0), z0);
    let r_low = t_low - 1.0;
    println!(
        "at resonance |t|={:.3e}, TL={:.2} dB; low frequency |t|={:.6}",
        t_res.abs(),
        side_branch_loss_db(res.lumped_impedance(&AIR, w0), z0),
        t_low.abs()
    );
    close(t_res.abs(), 0.0, 1e-9);
    assert!(t_low.abs() > 0.999);
    close(r_low.abs2() + t_low.abs2(), 1.0, 1e-9);
}

/// The periodic side-branch array has a sub-wavelength hybridisation gap at
/// `omega0` and a Bragg gap at `omega_B = pi c/a`.
#[test]
fn periodic_array_hybridization_and_bragg() {
    let res = Helmholtz::new(
        5.0e-3,
        5.0e-3,
        1.0e-4,
        EndCondition::Flanged,
        EndCondition::Flanged,
    );
    let a = 0.1;
    let arr = SideBranchArray::new(1.0e-2, a, res, AIR_TRANSPORT);
    let w0 = res.omega0(&AIR);
    let wb = PI * AIR.c / a;
    println!("omega0/2pi={:.1} Hz, Bragg/2pi={:.1} Hz", w0 / (2.0 * PI), wb / (2.0 * PI));
    // Find the gaps by scanning.
    let mut gaps: Vec<(f64, f64)> = Vec::new();
    let n = 200_000;
    let (lo, hi) = (0.2 * w0, 1.4 * wb);
    let mut start = f64::NAN;
    for i in 0..=n {
        let w = lo + (hi - lo) * i as f64 / n as f64;
        let g = arr.in_band_gap(&AIR, w);
        if g && start.is_nan() {
            start = w;
        }
        if !g && !start.is_nan() {
            gaps.push((start, w));
            start = f64::NAN;
        }
    }
    for (g0, g1) in &gaps {
        println!("  gap {:.1} .. {:.1} Hz", g0 / (2.0 * PI), g1 / (2.0 * PI));
    }
    let near = |target: f64| {
        gaps.iter().any(|&(g0, g1)| {
            let centre = 0.5 * (g0 + g1);
            (centre - target).abs() / target < 0.05
        })
    };
    assert!(near(w0), "no hybridisation gap near omega0");
    assert!(near(wb), "no Bragg gap near pi c / a");
}

/// Folded effective bulk modulus vs the exact Bloch wavevector. At long
/// wavelength (`k a << 1`) the Fang form with `F = V/(A a)` reproduces the
/// exact dispersion; the error grows with the validity indicator `k a`.
#[test]
fn effective_medium_vs_exact_bloch() {
    let res = Helmholtz::new(
        5.0e-3,
        5.0e-3,
        1.0e-4,
        EndCondition::Flanged,
        EndCondition::Flanged,
    );
    let a = 0.1;
    let arr = SideBranchArray::new(1.0e-2, a, res, AIR_TRANSPORT);
    let w0 = res.omega0(&AIR);
    for &frac in &[0.2, 0.4, 0.6, 0.8] {
        let w = frac * w0;
        let q = arr.bloch_wavevector(&AIR, w).re;
        let q_eff = arr.effective_medium_wavenumber(&AIR, w).re;
        let ka = arr.ka(&AIR, w);
        let err = (q - q_eff).abs() / q.abs().max(1e-30);
        println!("  {frac:.2}   {ka:.4}   {q:.6}   {q_eff:.6}   {err:.3e}");
    }
    let err_low = (arr.bloch_wavevector(&AIR, 0.2 * w0).re
        - arr.effective_medium_wavenumber(&AIR, 0.2 * w0).re)
        .abs()
        / arr.bloch_wavevector(&AIR, 0.2 * w0).re.abs();
    let w_hi = 0.8 * w0;
    let err_hi = (arr.bloch_wavevector(&AIR, w_hi).re
        - arr.effective_medium_wavenumber(&AIR, w_hi).re)
        .abs()
        / arr.bloch_wavevector(&AIR, w_hi).re.abs().max(1e-30);
    println!("low-ka error={err_low:.3e}, higher-ka error={err_hi:.3e}");
    assert!(err_low < 0.15, "effective medium should work at low ka");
    assert!(err_hi > err_low, "fold error should grow with ka");
}

/// Fang negative-bulk-modulus band: `Re B_eff < 0` on
/// `(omega0, omega0 sqrt(1+F))`, and the exact array is in a gap there.
#[test]
fn negative_bulk_modulus_band() {
    let res = Helmholtz::new(
        5.0e-3,
        5.0e-3,
        1.0e-4,
        EndCondition::Flanged,
        EndCondition::Flanged,
    );
    let a = 0.1;
    let arr = SideBranchArray::new(1.0e-2, a, res, AIR_TRANSPORT);
    let w0 = res.omega0(&AIR);
    let f = arr.filling_factor();
    let wtop = w0 * (1.0 + f).sqrt();
    let wmid = 0.5 * (w0 + wtop);
    let b_mid = arr.effective_bulk_modulus(&AIR, wmid);
    let b_below = arr.effective_bulk_modulus(&AIR, 0.5 * w0);
    println!(
        "F={f:.4}, negative-B band {:.1}..{:.1} Hz; B_eff(mid)={:?}, B_eff(0.5 w0)={:?}",
        w0 / (2.0 * PI),
        wtop / (2.0 * PI),
        b_mid,
        b_below
    );
    assert!(b_mid.re < 0.0, "B_eff must be negative in the band");
    assert!(b_below.re > 0.0, "B_eff must be positive below omega0");
    assert!(arr.in_band_gap(&AIR, wmid), "exact array should be in a gap");
    // At the top of the band 1/B_eff passes through zero (B_eff -> -inf).
    let inv_b = C::real(1.0) / arr.effective_bulk_modulus(&AIR, wtop);
    assert!(inv_b.abs() < 1e-9 * (1.0 / AIR.bulk()));
}

/// Water vs air and rigid vs elastic walls. The wall-to-fluid compliance ratio
/// `(3 B0)/(4 G)` (spherical) or `B0/G` (cylindrical) controls whether the
/// elastic wall matters: negligible for air in any solid, ~2% for water in
/// steel, and order one for water in a polymer.
#[test]
fn water_vs_air_and_wall_compliance() {
    for fluid in [&AIR, &WATER] {
        for solid in [&STEEL, &PMMA, &EPOXY] {
            let ratio = phononic::helmholtz::wall_to_fluid_compliance(
                fluid,
                solid,
                CavityShape::Spherical,
            );
            let b_rigid = fluid.bulk();
            let b_elastic = cavity_bulk_modulus_elastic(fluid, solid, CavityShape::Spherical);
            println!(
                "{} / {}: B0={:.3e} Pa, wall/fluid={:.3e}, B_elastic/B0={:.4}",
                fluid.name,
                solid.name,
                b_rigid,
                ratio,
                b_elastic / b_rigid
            );
        }
    }
    let air_steel =
        phononic::helmholtz::wall_to_fluid_compliance(&AIR, &STEEL, CavityShape::Spherical);
    let water_steel =
        phononic::helmholtz::wall_to_fluid_compliance(&WATER, &STEEL, CavityShape::Spherical);
    let water_pmma =
        phononic::helmholtz::wall_to_fluid_compliance(&WATER, &PMMA, CavityShape::Spherical);
    assert!(air_steel < 1e-4, "air is too soft for wall compliance to matter");
    assert!(water_steel < 0.1, "steel wall is a small correction for water");
    assert!(water_pmma > 0.5, "polymer wall materially softens a water cavity");

    // Folding the wall compliance shifts omega0 down by 1/sqrt(1+ratio).
    let res = Helmholtz::new(
        5.0e-3,
        5.0e-3,
        1.0e-4,
        EndCondition::Flanged,
        EndCondition::Flanged,
    );
    let w_rigid = res.omega0(&WATER);
    let b_elastic = cavity_bulk_modulus_elastic(&WATER, &PMMA, CavityShape::Spherical);
    // omega0 = c sqrt(S/(V Leff)); replacing B0 by B_elastic is equivalent to
    // c^2 -> B_elastic/rho.
    let c_eff = (b_elastic / WATER.rho).sqrt();
    let w_elastic = c_eff * (res.neck_area() / (res.cavity_volume * res.effective_length())).sqrt();
    println!(
        "water/PMMA cavity: f0 rigid={:.2} Hz, elastic={:.2} Hz, shift={:.1}%",
        w_rigid / (2.0 * PI),
        w_elastic / (2.0 * PI),
        100.0 * (w_elastic - w_rigid) / w_rigid
    );
    assert!(w_elastic < w_rigid);
}

/// The Fang helper and the array method agree, including the lossy form.
#[test]
fn fang_helper_consistency() {
    let res = Helmholtz::new(
        5.0e-3,
        5.0e-3,
        1.0e-4,
        EndCondition::Flanged,
        EndCondition::Flanged,
    );
    let mut arr = SideBranchArray::new(1.0e-2, 0.1, res, AIR_TRANSPORT);
    arr.gamma_loss = 25.0;
    let w = 1.1 * res.omega0(&AIR);
    let a = arr.effective_bulk_modulus(&AIR, w);
    let b = fang_bulk_modulus(
        AIR.bulk(),
        w,
        res.omega0(&AIR),
        arr.filling_factor(),
        arr.gamma_loss,
    );
    close((a - b).abs(), 0.0, 1e-9 * b.abs());
    assert!(a.im.abs() > 0.0, "loss must make B_eff complex");
}

/// Double-negative band: a membrane (negative effective density) overlapping
/// the Helmholtz negative-bulk-modulus band gives `Re rho < 0`, `Re B < 0` and
/// a negative effective index.
#[test]
fn double_negative_band() {
    let w0 = 1.0_f64;
    let filling = 0.3;
    let b_neg = fang_bulk_modulus(1.0, 1.05, w0, filling, 0.0);
    let mem = MembraneDensity::new(0.95, 0.3, 0.0);
    let rho_neg = mem.effective_density(1.0, 1.05);
    println!("at w=1.05: B_eff={b_neg:?}, rho_eff={rho_neg:?}");
    assert!(b_neg.re < 0.0);
    assert!(rho_neg.re < 0.0);
    let n = double_negative_index(rho_neg, b_neg);
    println!("double-negative index n={n:?}");
    assert!(n.re < 0.0, "index should be negative in the double-negative band");
    // At w=1.138 only the bulk modulus is negative; the membrane density has
    // already turned positive again (rho band top 1.135).
    let b_only = fang_bulk_modulus(1.0, 1.138, w0, filling, 0.0);
    let rho_pos = mem.effective_density(1.0, 1.138);
    println!("at w=1.138: B_eff={b_only:?}, rho_eff={rho_pos:?}");
    assert!(b_only.re < 0.0 && rho_pos.re > 0.0);
}

/// A finite array of resonators has a deep transmission-loss dip inside the
/// hybridisation band gap predicted by the exact Bloch relation. This validates
/// the periodic transfer matrix against a directly measurable TL.
#[test]
fn finite_array_transmission_dip_in_gap() {
    let res = Helmholtz::new(
        5.0e-3,
        5.0e-3,
        1.0e-4,
        EndCondition::Flanged,
        EndCondition::Flanged,
    );
    let a = 0.1;
    let arr = SideBranchArray::new(1.0e-2, a, res, AIR_TRANSPORT);
    let w0 = res.omega0(&AIR);
    let cells = 20;
    // Scan for the largest TL near the hybridisation gap.
    let n = 20_000;
    let (lo, hi) = (0.85 * w0, 1.15 * w0);
    let (mut best_w, mut best_tl) = (f64::NAN, 0.0f64);
    let mut in_gap_lo = f64::NAN;
    let mut in_gap_hi = f64::NAN;
    for i in 0..=n {
        let w = lo + (hi - lo) * i as f64 / n as f64;
        let tl = arr.stack_loss_db(&AIR, w, cells);
        if tl > best_tl {
            best_tl = tl;
            best_w = w;
        }
        if arr.in_band_gap(&AIR, w) {
            if in_gap_lo.is_nan() {
                in_gap_lo = w;
            }
            in_gap_hi = w;
        }
    }
    println!(
        "20-cell dip: f={:.1} Hz, TL={:.1} dB; Bloch gap {:.1}..{:.1} Hz",
        best_w / (2.0 * PI),
        best_tl,
        in_gap_lo / (2.0 * PI),
        in_gap_hi / (2.0 * PI)
    );
    assert!(best_tl > 40.0, "the resonance array must give a deep TL dip");
    assert!(
        best_w > in_gap_lo && best_w < in_gap_hi,
        "TL dip must lie inside the Bloch gap"
    );
}

/// The same small resonator filled with water instead of air resonates at
/// `c_water/c_air` times the air frequency (the geometry is unchanged).
#[test]
fn water_vs_air_resonance_scaling() {
    let res = Helmholtz::new(
        5.0e-3,
        5.0e-3,
        1.0e-4,
        EndCondition::Flanged,
        EndCondition::Flanged,
    );
    let f_air = res.frequency0(&AIR);
    let f_water = res.frequency0(&WATER);
    let ratio = f_water / f_air;
    println!("f0 air={f_air:.1} Hz, water={f_water:.1} Hz, ratio={ratio:.6}");
    close(ratio, WATER.c / AIR.c, 1e-12);
}
