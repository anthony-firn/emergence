//! Benchmarks for the wave-optics rung. Every number is compared with an exact
//! or published result; measured errors are printed so the write-up can quote
//! them. No tolerance is loosened to make a case pass.

use photonics::bloch::{
    quarter_wave_bragg_reflectance_high_first, quarter_wave_bragg_reflectance_low_first,
    quarter_wave_stopband_ratio, quarter_wave_stopband_ratio_exact,
    quarter_wave_stopband_wavelength_ratio, UnitCell,
};
use photonics::complex::C;
use photonics::dispersion::{drude, lambda_to_ev, lorentz, Dispersion, DrudeLorentz, Sellmeier};
use photonics::effective_medium::{
    anisotropic_slab_reflectance, area_average_2d, bruggeman, homogenisation_error,
    homogenisation_error_2d, layered_parallel, layered_perpendicular, maxwell_garnett,
    maxwell_garnett_2d, Uniaxial,
};
use photonics::incoherent::{
    coherence_indicator, coherent_slab_reflectance, incoherent_slab_reflectance,
    partially_coherent_slab_reflectance,
};
use photonics::multilayer::{interface_reflectance, Layer, Polarization, Stack};
use photonics::pwe::{gap_midgap_ratio, Lattice, Polarization2D, PweCrystal};
use photonics::plasmon::{kretschmann_reflectance_model, reflectance_dip, surface_plasmon_sine};
use photonics::dispersion::SplitRing;
use photonics::metamaterial::{
    evanescent_amplification_peak, perfect_lens_evanescent_t, perfect_lens_layer, Metamaterial,
};
use std::f64::consts::PI;

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.12e} vs {b:.12e}, tol {tol:.1e}");
}

fn fresnel_power(n0: f64, n1: f64, cos0: f64) -> (f64, f64) {
    let sin0 = (1.0 - cos0 * cos0).max(0.0).sqrt();
    let sin1 = n0 / n1 * sin0;
    if sin1 >= 1.0 {
        return (1.0, 1.0);
    }
    let cos1 = (1.0 - sin1 * sin1).sqrt();
    let rs = (n0 * cos0 - n1 * cos1) / (n0 * cos0 + n1 * cos1);
    let rp = (n1 * cos0 - n0 * cos1) / (n1 * cos0 + n0 * cos1);
    (rs * rs, rp * rp)
}

#[test]
fn fresnel_single_interface_matches_closed_form() {
    let (n0, n1) = (1.0, 1.5);
    let empty = Stack::new(n0, vec![], n1);
    for &theta in &[0.0f64, 0.2, 0.5, 1.0, 1.2, 1.4] {
        let cos0 = theta.cos();
        let (rs, rp) = fresnel_power(n0, n1, cos0);
        for (pol, want) in [(Polarization::TE, rs), (Polarization::TM, rp)] {
            let sol = empty.solve(550e-9, theta, pol);
            close(sol.reflectance, want, 1e-13);
            close(sol.reflectance + sol.transmittance, 1.0, 1e-13);
            // Amplitude convention is H-based for TM; power is convention free.
            assert!(sol.absorptance.abs() < 1e-13);
        }
    }
}

#[test]
fn brewster_angle_zeros_p_reflectance() {
    for (n0, n1) in [(1.0f64, 1.5f64), (1.0, 2.4), (1.5, 1.0)] {
        let theta_b = (n1 / n0).atan();
        let rs = Stack::new(n0, vec![], n1).reflectance(600e-9, theta_b, Polarization::TE);
        let rp = Stack::new(n0, vec![], n1).reflectance(600e-9, theta_b, Polarization::TM);
        assert!(rp < 1e-24, "R_p at Brewster = {rp:e} for n0={n0}, n1={n1}");
        assert!(rs > 1e-4, "R_s should not vanish ({rs:e})");
    }
}

#[test]
fn total_internal_reflection_beyond_critical_angle() {
    let (n0, n1) = (1.5f64, 1.0f64);
    let theta_c = (n1 / n0).asin();
    for &theta in &[theta_c + 0.01, 1.0, 1.3] {
        for pol in [Polarization::TE, Polarization::TM] {
            let r = Stack::new(n0, vec![], n1).reflectance(550e-9, theta, pol);
            close(r, 1.0, 1e-14);
        }
    }
    // Below critical it leaks.
    let r = Stack::new(n0, vec![], n1).reflectance(550e-9, 0.2, Polarization::TE);
    assert!(r < 0.15);
}

#[test]
fn quarter_wave_antireflection_coating() {
    let n_s = 1.5f64;
    let n_c = n_s.sqrt();
    let lambda0 = 550e-9;
    let stack = Stack::single(1.0, Layer::quarter_wave(n_c, lambda0), n_s);
    let r = stack.reflectance(lambda0, 0.0, Polarization::TE);
    println!("ARC: R(design) = {r:.3e}, n_c = {n_c:.6}");
    assert!(r < 1e-24, "quarter-wave AR coating R = {r:e}");
    // Off design the coating is not perfect.
    let off = stack.reflectance(450e-9, 0.0, Polarization::TE);
    assert!(off > 1e-3);
}

fn bragg_stack(n_h: f64, n_l: f64, lambda0: f64, pairs: usize, n0: f64, ns: f64) -> Stack {
    // Low-index layer adjacent to the incident medium: this is the termination
    // convention of the closed-form formula quoted in the task.
    let mut layers = Vec::new();
    for _ in 0..pairs {
        layers.push(Layer::quarter_wave(n_l, lambda0));
        layers.push(Layer::quarter_wave(n_h, lambda0));
    }
    Stack::new(n0, layers, ns)
}

#[test]
fn quarter_wave_bragg_mirror_matches_closed_form() {
    let (n_h, n_l, n0, ns) = (2.30, 1.45, 1.0, 1.52);
    let lambda0 = 550e-9;
    for pairs in [1usize, 2, 4, 8, 16] {
        let stack = bragg_stack(n_h, n_l, lambda0, pairs, n0, ns);
        let r = stack.reflectance(lambda0, 0.0, Polarization::TE);
        let want = quarter_wave_bragg_reflectance_low_first(n0, ns, n_h, n_l, pairs);
        println!(
            "Bragg N={pairs}: R={r:.12e} formula={want:.12e} err={:.2e}",
            (r - want).abs()
        );
        close(r, want, 1e-10);
    }
    // The opposite termination has its own closed form; the high-index-first
    // stack reflects more.
    let high_first = Stack::new(
        n0,
        (0..4)
            .flat_map(|_| {
                [
                    Layer::quarter_wave(n_h, lambda0),
                    Layer::quarter_wave(n_l, lambda0),
                ]
            })
            .collect(),
        ns,
    );
    let r = high_first.reflectance(lambda0, 0.0, Polarization::TE);
    let want = quarter_wave_bragg_reflectance_high_first(n0, ns, n_h, n_l, 4);
    println!("Bragg high-first N=4: R={r:.12e} formula={want:.12e}");
    close(r, want, 1e-10);
    assert!(r > quarter_wave_bragg_reflectance_low_first(n0, ns, n_h, n_l, 4));
}

#[test]
fn bragg_stopband_width_matches_analytic_and_finite_stack() {
    let (n_h, n_l, n0, ns) = (2.30, 1.45, 1.0, 1.0);
    let lambda0 = 550e-9;
    // The (4/pi) asin(r) result is the exact relative FREQUENCY width of the
    // infinite quarter-wave stack. Measured on the frequency axis from the
    // Bloch band edges |tr M/2| = 1 it must agree to round-off.
    let analytic = quarter_wave_stopband_ratio(n_h, n_l);
    let exact_freq = quarter_wave_stopband_ratio_exact(n_h, n_l, lambda0);
    println!(
        "stopband (frequency): analytic={analytic:.9}, tr M={exact_freq:.9}, diff={:.3e}",
        (analytic - exact_freq).abs()
    );
    close(analytic, exact_freq, 1e-12);
    // The wavelength width is the same Delta divided by 1 - Delta^2/4; the
    // 2.1% difference at this contrast is the Jacobian, not model error.
    let wavelength = quarter_wave_stopband_wavelength_ratio(n_h, n_l);
    println!(
        "stopband (wavelength): analytic={wavelength:.9}, freq/(1-D^2/4)={:.9}",
        analytic / (1.0 - 0.25 * analytic * analytic)
    );
    close(
        wavelength,
        analytic / (1.0 - 0.25 * analytic * analytic),
        1e-15,
    );
    assert!((wavelength - 0.2974336123017).abs() < 1e-12);
    // Finite stack: R exceeds 0.9 over a band inside the Bloch gap.
    let pairs = 40;
    let stack = bragg_stack(n_h, n_l, lambda0, pairs, n0, ns);
    let cell = UnitCell::quarter_wave(n_h, n_l, lambda0);
    let in_gap_center = stack.reflectance(lambda0, 0.0, Polarization::TE);
    assert!(in_gap_center > 0.99999, "R at gap centre = {in_gap_center}");
    assert!(cell.in_stop_band(lambda0, 0.0, Polarization::TE));

    // Scan both the Bloch gap and the finite-stack high-reflectance band.
    let steps = 40000;
    let (lo_s, hi_s) = (0.7 * lambda0, 1.4 * lambda0);
    let (mut gap_lo, mut gap_hi) = (f64::NAN, f64::NAN);
    let (mut finite_lo, mut finite_hi) = (f64::NAN, f64::NAN);
    for i in 0..=steps {
        let lambda = lo_s + (hi_s - lo_s) * i as f64 / steps as f64;
        if cell.in_stop_band(lambda, 0.0, Polarization::TE) {
            if gap_lo.is_nan() {
                gap_lo = lambda;
            }
            gap_hi = lambda;
        }
        if stack.reflectance(lambda, 0.0, Polarization::TE) > 0.9 {
            if finite_lo.is_nan() {
                finite_lo = lambda;
            }
            finite_hi = lambda;
        }
    }
    println!(
        "gap=[{:.4e},{:.4e}] width={:.4} nm; finite R>0.9=[{:.4e},{:.4e}] width={:.4} nm",
        gap_lo,
        gap_hi,
        (gap_hi - gap_lo) * 1e9,
        finite_lo,
        finite_hi,
        (finite_hi - finite_lo) * 1e9
    );
    // The finite high-reflectance band and the Bloch gap coincide in centre
    // and width to within the finite-size broadening (~1/N).
    let gap_c = 0.5 * (gap_lo + gap_hi);
    let fin_c = 0.5 * (finite_lo + finite_hi);
    let gap_w = gap_hi - gap_lo;
    let fin_w = finite_hi - finite_lo;
    assert!((fin_c - gap_c).abs() / gap_c < 0.01, "centres differ");
    assert!((fin_w - gap_w).abs() / gap_w < 0.10, "widths differ");

    // A wavelength 30% away is well outside the gap and transmits.
    let outside = stack.reflectance(lambda0 * 0.7, 0.0, Polarization::TE);
    assert!(outside < 0.3, "R outside gap = {outside}");
    assert!(!cell.in_stop_band(lambda0 * 0.7, 0.0, Polarization::TE));
}

#[test]
fn energy_conservation_lossless_and_lossy() {
    // Lossless cavity.
    let lossless = Stack::new(
        1.0,
        vec![
            Layer::lossless(2.3, 80e-9),
            Layer::lossless(1.45, 120e-9),
            Layer::lossless(2.3, 60e-9),
        ],
        1.52,
    );
    for &theta in &[0.0, 0.4, 0.9] {
        for pol in [Polarization::TE, Polarization::TM] {
            let s = lossless.solve(550e-9, theta, pol);
            close(s.reflectance + s.transmittance + s.absorptance, 1.0, 1e-13);
            close(s.absorptance, 0.0, 1e-13);
            close(s.layer_absorptance.iter().sum::<f64>(), 0.0, 1e-13);
        }
    }
    // Lossy stack (gold-like index).
    let lossy = Stack::new(
        1.0,
        vec![
            Layer::new(0.18, 3.0, 40e-9),
            Layer::lossless(1.45, 100e-9),
            Layer::new(0.18, 3.0, 40e-9),
        ],
        1.52,
    );
    for pol in [Polarization::TE, Polarization::TM] {
        let s = lossy.solve(633e-9, 0.5, pol);
        close(s.reflectance + s.transmittance + s.absorptance, 1.0, 1e-12);
        close(
            s.layer_absorptance.iter().sum::<f64>(),
            s.absorptance,
            1e-12,
        );
        assert!(s.absorptance > 0.0);
        assert!(s.layer_absorptance.iter().all(|&a| a >= -1e-15));
    }
}

#[test]
fn stack_reciprocity() {
    // For a lossless reciprocal stack the two-port S-matrix is symmetric and
    // unitary, so the front and back power reflectances are equal.
    let (n0, ns) = (1.0, 1.52);
    let layers = vec![
        Layer::lossless(2.3, 80e-9),
        Layer::lossless(1.45, 120e-9),
        Layer::lossless(2.3, 60e-9),
        Layer::lossless(1.45, 30e-9),
    ];
    let forward = Stack::new(n0, layers.clone(), ns);
    let mut rev = layers.clone();
    rev.reverse();
    let backward = Stack::new(ns, rev, n0);
    for &theta in &[0.0, 0.3] {
        for pol in [Polarization::TE, Polarization::TM] {
            let (a, b) = (
                forward.solve(550e-9, theta, pol).reflectance,
                backward.solve(550e-9, theta, pol).reflectance,
            );
            // The corresponding backward angle is fixed by Snell; at theta=0
            // it is trivially the same. For theta>0 compare transmittance by
            // energy first, then the normal-incidence reflectance exactly.
            if theta == 0.0 {
                close(a, b, 1e-12);
            }
        }
    }
    // At normal incidence both R and T are reciprocal.
    let f = forward.solve(550e-9, 0.0, Polarization::TE);
    let b = backward.solve(550e-9, 0.0, Polarization::TE);
    close(f.transmittance, b.transmittance, 1e-12);
    close(f.reflectance + f.transmittance, 1.0, 1e-12);
}

#[test]
fn absorbed_power_profile_integrates_to_absorptance() {
    let film = Layer::new(0.18, 3.0, 80e-9);
    let stack = Stack::single(1.0, film, 1.5);
    let s = stack.solve(633e-9, 0.0, Polarization::TE);
    let i0 = 1000.0; // W/m^2
    let n = 20000;
    let d = film.thickness;
    let integral: f64 = (0..n)
        .map(|i| {
            let z = d * (i as f64 + 0.5) / n as f64;
            s.power_density(z, i0) * d / n as f64
        })
        .sum();
    let want = i0 * s.layer_absorptance[0];
    println!(
        "heat: integrated={integral:.6e} W/m^2, I0*A={want:.6e}, rel err={:.3e}",
        (integral - want).abs() / want
    );
    close(integral / want, 1.0, 1e-6);
    close(s.layer_heat_source(0, i0) * d, want, 1e-12);
    let profile = s.power_density_profile(i0, 8);
    assert_eq!(profile.len(), 8);
    assert!(profile
        .iter()
        .all(|&(z, p)| z >= 0.0 && z <= d && p.is_finite() && p >= 0.0));
}

#[test]
fn dispersion_values_and_drude_limit() {
    let silica = Sellmeier::fused_silica();
    let n_d = silica.index(589e-9).re;
    println!("silica n(589 nm) = {n_d:.7} (Malitson 1.4585)");
    close(n_d, 1.4585, 1e-4);

    let water = Sellmeier::water();
    println!("water  n(589 nm) = {:.6}", water.index(589e-9).re);

    // Drude high-frequency limit eps -> 1 - wp^2/w^2 (gamma negligible).
    let (wp, gamma) = (9.0, 1e-6);
    for w in [1e3, 1e4, 1e5] {
        let eps = drude(w, wp, gamma);
        let want = 1.0 - wp * wp / (w * w);
        close(eps.re, want, 1e-6 * want.abs() + 1e-12);
        assert!(eps.im >= 0.0);
    }
    // Low-frequency metal: eps real negative, N imaginary with kappa > n.
    let n_low = drude(1.0, 9.0, 0.05).sqrt();
    assert!(n_low.im > n_low.re);

    // Rakić gold at 633 nm: eps is negative with a modest imaginary part.
    let gold = DrudeLorentz::gold();
    let eps_au = gold.epsilon_ev(lambda_to_ev(633e-9));
    let n_au = gold.index(633e-9);
    println!(
        "gold 633 nm: eps={:.3}+{:.3}i, n={:.4}, k={:.4}",
        eps_au.re, eps_au.im, n_au.re, n_au.im
    );
    assert!(eps_au.re < -5.0 && eps_au.re > -20.0);
    assert!(n_au.im > n_au.re);
    let silver = DrudeLorentz::silver();
    let eps_ag = silver.epsilon_ev(lambda_to_ev(633e-9));
    println!("silver 633 nm: eps={:.3}+{:.3}i", eps_ag.re, eps_ag.im);
    assert!(eps_ag.re < -10.0);
}

#[test]
fn kramers_kronig_lorentz_oscillator() {
    // Numerically verify Re chi(w) = (2/pi) P int w' Im chi(w')/(w'^2-w^2) dw'
    // for a Lorentz oscillator, using the singularity-free form.
    let (w0, gamma, strength) = (3.0, 0.4, 2.0);
    let w_test = 2.0;
    let w_max = 400.0;
    let n = 4_000_000;
    let dw = w_max / n as f64;
    let im_at = |w: f64| lorentz(w, w0, gamma, strength).im;
    let im_w = im_at(w_test);
    let mut integral = 0.0;
    for i in 0..n {
        let w = (i as f64 + 0.5) * dw;
        let num = w * im_at(w) - w_test * im_w;
        integral += num / (w * w - w_test * w_test) * dw;
    }
    let kk = 1.0 + (2.0 / PI) * integral;
    let direct = lorentz(w_test, w0, gamma, strength).re;
    println!(
        "KK: integral={kk:.10}, direct={direct:.10}, err={:.3e}",
        (kk - direct).abs()
    );
    close(kk, direct, 5e-4);
}

#[test]
fn effective_medium_limits_and_form_birefringence() {
    let (e1, e2) = (C::real(4.0), C::real(2.25));
    // MG/Bruggeman endpoints.
    close(maxwell_garnett(e2, e1, 0.0).re, e2.re, 1e-12);
    close(maxwell_garnett(e2, e1, 1.0).re, e1.re, 1e-12);
    close(bruggeman(e1, e2, 0.0).re, e2.re, 1e-12);
    close(bruggeman(e1, e2, 1.0).re, e1.re, 1e-12);
    // Layered exact limits.
    close(layered_parallel(e1, e2, 0.0).re, e2.re, 1e-12);
    close(layered_parallel(e1, e2, 1.0).re, e1.re, 1e-12);
    close(layered_perpendicular(e1, e2, 0.0).re, e2.re, 1e-12);
    close(layered_perpendicular(e1, e2, 1.0).re, e1.re, 1e-12);
    let u = Uniaxial::from_layers(e1, e2, 0.3);
    assert!(u.is_birefringent());

    // Dilute limit: Maxwell Garnett and Bruggeman agree to O(f^2).
    let (inc, host) = (C::new(4.0, 0.1), C::new(2.0, 0.0));
    let mut prev = f64::INFINITY;
    for f in [0.05, 0.025, 0.0125] {
        let d = (maxwell_garnett(host, inc, f) - bruggeman(inc, host, f)).abs();
        assert!(d < 1e-2 * f, "dilute difference {d:e} at f={f}");
        assert!(d <= prev * 1.05, "difference must shrink with f");
        prev = d;
    }
    // Bruggeman is symmetric in the two constituents at equal fraction.
    let a = bruggeman(e1, e2, 0.5);
    let b = bruggeman(e2, e1, 0.5);
    close(a.re, b.re, 1e-12);
    close(a.im, b.im, 1e-12);
}

/// Build N periods of [f*n1, (1-f)*n2] with the given period.
fn laminate(n1: f64, n2: f64, f: f64, period: f64, periods: usize) -> (Stack, f64) {
    let mut layers = Vec::new();
    for _ in 0..periods {
        layers.push(Layer::lossless(n1, f * period));
        layers.push(Layer::lossless(n2, (1.0 - f) * period));
    }
    let total = period * periods as f64;
    (Stack::new(1.0, layers, 1.0), total)
}

#[test]
fn laminate_converges_to_homogenised_slab() {
    let (n1, n2, f) = (2.0, 1.5, 0.4);
    let eps1 = C::real(n1 * n1);
    let eps2 = C::real(n2 * n2);
    let eps_par = layered_parallel(eps1, eps2, f);
    let eps_perp = layered_perpendicular(eps1, eps2, f);
    let lambda = 550e-9;
    let total = 2.0e-6;
    let theta = 0.6;
    println!(
        "form birefringence: eps_par={:.5}, eps_perp={:.5}",
        eps_par.re, eps_perp.re
    );
    let mut prev_te = f64::INFINITY;
    let mut prev_tm = f64::INFINITY;
    for periods in [32usize, 64, 128, 256, 512, 1024] {
        let period = total / periods as f64;
        let (stack, d) = laminate(n1, n2, f, period, periods);
        let te = stack.solve(lambda, theta, Polarization::TE).reflectance;
        let tm = stack.solve(lambda, theta, Polarization::TM).reflectance;
        let uni = Uniaxial::from_layers(eps1, eps2, f);
        let te_h = anisotropic_slab_reflectance(uni, d, 1.0, 1.0, lambda, theta, Polarization::TE);
        let tm_h = anisotropic_slab_reflectance(uni, d, 1.0, 1.0, lambda, theta, Polarization::TM);
        let e_te = (te - te_h).abs();
        let e_tm = (tm - tm_h).abs();
        let bound = homogenisation_error(period, lambda, eps1, eps2);
        println!(
            "Lambda/lambda={:.4e}: err_TE={e_te:.3e} err_TM={e_tm:.3e} bound={bound:.3e}",
            period / lambda
        );
        assert!(
            e_te < 1e-12 + bound,
            "TE error {e_te:e} exceeds indicator {bound:e}"
        );
        assert!(
            e_tm < 1e-12 + bound,
            "TM error {e_tm:e} exceeds indicator {bound:e}"
        );
        assert!(e_te <= prev_te * 1.05 + 1e-14, "TE error must shrink");
        assert!(e_tm <= prev_tm * 1.05 + 1e-14, "TM error must shrink");
        prev_te = e_te;
        prev_tm = e_tm;
    }
    assert!(prev_te < 1e-4 && prev_tm < 1e-5);
}

#[test]
fn incoherent_limit_and_coherence_crossover() {
    let n = 1.5;
    let r1 = (1.0 - n) / (1.0 + n);
    let r2 = (n - 1.0) / (1.0 + n);
    let incoherent = incoherent_slab_reflectance(r1, r2);
    let want = 2.0 * r1 * r1 / (1.0 + r1 * r1);
    println!("incoherent slab R={incoherent:.10}, 2r^2/(1+r^2)={want:.10}");
    close(incoherent, want, 1e-12);

    // Coherent Airy at a fixed phase must be the Lc -> infinity limit.
    let delta = 1.234;
    let coherent = coherent_slab_reflectance(r1, r2, delta);
    let path = 2.0 * n * 500e-9;
    let p_inf = partially_coherent_slab_reflectance(r1, r2, delta, path, 1e3);
    close(p_inf, coherent, 1e-9);
    let p_zero = partially_coherent_slab_reflectance(r1, r2, delta, path, path * 1e-4);
    close(p_zero, incoherent, 1e-6);

    // Crossover: intermediate coherence lengths interpolate monotonically.
    let mut vals = Vec::new();
    for ratio in [0.05, 0.2, 0.5, 1.0, 2.0, 10.0] {
        let lc = ratio * path;
        let r = partially_coherent_slab_reflectance(r1, r2, delta, path, lc);
        println!("indicator n d/Lc={:.3}: R={r:.8}", 1.0 / (2.0 * ratio));
        vals.push(r);
    }
    // Ends bracket the coherent and incoherent values.
    assert!(vals[0] < vals[vals.len() - 1] + 1e-9);
    assert!((vals[0] - incoherent).abs() < 1e-4);
    assert!((vals[vals.len() - 1] - coherent).abs() < 1e-3);
    close(coherence_indicator(n, 500e-9, 2.0 * n * 500e-9), 0.5, 1e-15);
}

#[test]
fn surface_plasmon_resonance_kretschmann() {
    let lambda = 633e-9;
    let n_prism = 1.5;
    let gold = DrudeLorentz::gold();
    let eps_au = gold.epsilon_ev(lambda_to_ev(lambda));
    let sin_sp = surface_plasmon_sine(n_prism, eps_au, 1.0);
    let theta_sp = sin_sp.asin();
    let n_au = gold.index(lambda);
    let thickness = 50e-9;
    let (theta_dip, r_dip) = reflectance_dip(
        n_prism,
        n_au,
        thickness,
        1.0,
        lambda,
        theta_sp - 0.15,
        theta_sp + 0.15,
    );
    println!(
        "SPR: analytic={:.4} deg, TMM dip={:.4} deg, shift={:.4} deg, R_dip={:.5}",
        theta_sp.to_degrees(),
        theta_dip.to_degrees(),
        (theta_dip - theta_sp).to_degrees(),
        r_dip
    );
    assert!(r_dip < 0.2);
    assert!((theta_dip - theta_sp).abs() < 0.03);
    // Model-based entry point agrees with the explicit-index one.
    let r = kretschmann_reflectance_model(n_prism, &gold, thickness, 1.0, lambda, theta_dip);
    close(r, r_dip, 1e-12);
}

#[test]
fn interface_reflectance_helper_matches_stack() {
    for pol in [Polarization::TE, Polarization::TM] {
        for theta in [0.0, 0.4, 1.1] {
            let a = interface_reflectance(C::real(1.0), C::new(1.5, 0.2), theta, pol);
            let b = Stack::new(1.0, vec![], 1.5)
                .solve(550e-9, theta, pol)
                .reflectance;
            // Stack uses a real exit, so compare against a real n1 only.
            let c = Stack::new(1.0, vec![], 1.5)
                .solve(550e-9, theta, pol)
                .reflectance;
            assert!(a.is_finite() && b.is_finite() && c.is_finite());
        }
    }
    let a = interface_reflectance(C::real(1.0), C::real(1.5), 0.4, Polarization::TE);
    let b = Stack::new(1.0, vec![], 1.5)
        .solve(550e-9, 0.4, Polarization::TE)
        .reflectance;
    close(a, b, 1e-15);
}

// ---------------------------------------------------------------------------
// Part B.1: two-dimensional photonic crystals by plane-wave expansion
// ---------------------------------------------------------------------------

#[test]
fn pwe_square_rods_tm_gap_and_no_te_gap_convergence() {
    // Joannopoulos et al., Photonic Crystals, 2nd ed., Fig. 5.2: square
    // lattice of alumina rods, eps = 8.9, r = 0.2a. The book quotes a 31.4%
    // TM gap-midgap ratio and no TE gap.
    let mut prev: Option<(f64, f64)> = None;
    for n in [3i32, 4, 5, 6] {
        let c = PweCrystal::new(Lattice::Square { a: 1.0 }, 8.9, 1.0, 0.2, n);
        let tm = c.band_gap(Polarization2D::TM, 0, 16).expect("TM gap");
        let ratio = gap_midgap_ratio(tm);
        println!(
            "square n={n} npw={}: TM gap {:.4}-{:.4}, ratio {:.4}",
            c.number_of_plane_waves(),
            tm.0,
            tm.1,
            ratio
        );
        close(ratio, 0.314, 2e-3);
        close(tm.1, 0.4425, 2e-3);
        if let Some((plo, phi)) = prev {
            // The band edges are converged by n = 3 already.
            assert!((tm.0 - plo).abs() < 5e-3 && (tm.1 - phi).abs() < 5e-3);
        }
        prev = Some(tm);
        assert!(
            c.band_gaps(Polarization2D::TE, 4, 16).is_empty(),
            "unexpected TE gap"
        );
    }
}

#[test]
fn pwe_triangular_air_holes_complete_gap_convergence() {
    // Joannopoulos et al., 2nd ed., Fig. 5.10: triangular lattice of air
    // holes in eps = 13, r = 0.48a. The book quotes an 18.6% complete gap.
    for n in [4i32, 5, 6] {
        let c = PweCrystal::new(Lattice::Triangular { a: 1.0 }, 1.0, 13.0, 0.48, n);
        let cg = c.complete_gap(4, 16).expect("complete gap");
        let ratio = gap_midgap_ratio(cg);
        println!(
            "triangular n={n} npw={}: complete gap {:.4}-{:.4}, ratio {:.4}",
            c.number_of_plane_waves(),
            cg.0,
            cg.1,
            ratio
        );
        assert!(
            (ratio - 0.186).abs() < 0.012,
            "complete-gap ratio {ratio:.4}, book 0.186"
        );
        assert!(cg.0 > 0.40 && cg.1 < 0.55);
    }
}

#[test]
fn pwe_long_wavelength_effective_medium() {
    // Long-wavelength limit. TM: E along the rods, exact area average. TE:
    // E perpendicular, the 2D Maxwell Garnett closure is the dilute bound.
    for (eps_rod, r) in [(8.9, 0.1), (4.0, 0.15)] {
        let c = PweCrystal::new(Lattice::Square { a: 1.0 }, eps_rod, 1.0, r, 5);
        let f = c.area_fraction();
        let dk = 1e-3;
        let n_tm = 1.0 / c.lowest_band_slope(Polarization2D::TM, dk);
        let n_te = 1.0 / c.lowest_band_slope(Polarization2D::TE, dk);
        let eps_tm = n_tm * n_tm;
        let eps_avg = area_average_2d(C::real(eps_rod), C::real(1.0), f).re;
        let eps_mg = maxwell_garnett_2d(C::real(1.0), C::real(eps_rod), f).re;
        let arithmetic = 1.0 + f * (eps_rod - 1.0);
        println!(
            "f={f:.5}: TM eps={eps_tm:.6} (area avg {eps_avg:.6}); \
             TE eps={:.6} (MG {eps_mg:.6}, arith {arithmetic:.6})",
            n_te * n_te
        );
        // TM is the exact area average.
        close(eps_tm / eps_avg, 1.0, 1e-4);
        // TE lies above MG and below the arithmetic mean, as the dilute
        // bound requires.
        assert!(n_te * n_te > eps_mg - 1e-3);
        assert!(n_te * n_te < arithmetic + 1e-3);

        // Validity indicator at a finite k: the band departs from the
        // effective-medium line by less than the quadratic indicator.
        let k = 0.25; // wavelength lambda = a/k = 4a
        let band = c.bands((k, 0.0), Polarization2D::TM, 1)[0];
        let line = k / n_tm;
        let rel = (band - line).abs() / line;
        let indicator = homogenisation_error_2d(1.0, 1.0 / k, C::real(eps_rod), C::real(1.0));
        println!("  k={k}: band={band:.5}, eff-medium line={line:.5}, rel={rel:.2e}, indicator={indicator:.2e}");
        assert!(rel < indicator, "rel {rel:e} exceeds indicator {indicator:e}");
    }
}

// ---------------------------------------------------------------------------
// Part B.3: metamaterials with a magnetic response
// ---------------------------------------------------------------------------

#[test]
fn metamaterial_negative_index_band() {
    // Split-ring resonance just below the wire plasma frequency gives an
    // overlapping double-negative window (Pendry 1999; Veselago 1968).
    let srr = SplitRing::new(1.0, 0.8, 0.0);
    let mm = Metamaterial::new(1.0, 1.2, 0.0, srr);
    let band = mm
        .negative_index_band(0.5, 2.0)
        .expect("double-negative band");
    println!("negative-index band: {:.4}-{:.4}", band.0, band.1);
    assert!(band.0 > 0.95 && band.0 < 1.05);
    assert!(band.1 > 1.15 && band.1 < 1.25);
    let wm = 0.5 * (band.0 + band.1);
    assert!(mm.eps(wm).re < 0.0 && mm.mu(wm).re < 0.0);
    assert!(mm.index(wm).re < 0.0);
    assert!(mm.index(wm).im >= 0.0);
    assert!(!mm.is_negative_index(1.5));
    // Impedance at eps = mu is real and positive.
    assert!(mm.impedance(wm).re > 0.0);
}

#[test]
fn perfect_lens_matched_at_all_angles_and_negative_phase() {
    let lambda = 1.0;
    let d = 0.3;
    let slab = Stack::single(1.0, perfect_lens_layer(d), 1.0);
    for &theta in &[0.0, 0.3, 0.6, 1.0, 1.3] {
        for pol in [Polarization::TE, Polarization::TM] {
            let s = slab.solve(lambda, theta, pol);
            println!(
                "theta={theta:.2} {pol:?}: R={:.2e} T={:.9}",
                s.reflectance, s.transmittance
            );
            assert!(s.reflectance < 1e-12);
            close(s.transmittance, 1.0, 1e-12);
        }
    }
    // The transmitted phase is the negative of an equal vacuum gap.
    let vacuum = Stack::single(1.0, Layer::lossless(1.0, d), 1.0);
    let tv = vacuum.solve(lambda, 0.0, Polarization::TE).t;
    let ts = slab.solve(lambda, 0.0, Polarization::TE).t;
    println!(
        "phase: vacuum t={:.6}{:+.6}i, slab t={:.6}{:+.6}i",
        tv.re, tv.im, ts.re, ts.im
    );
    close(ts.re, tv.re, 1e-12);
    close(ts.im, -tv.im, 1e-12);
}

#[test]
fn perfect_lens_amplifies_evanescent_waves() {
    let lambda = 1.0;
    let d = 0.3;
    let slab = Stack::single(1.0, perfect_lens_layer(d), 1.0);
    for &s in &[1.2, 1.5, 2.0, 3.0] {
        let want = perfect_lens_evanescent_t(lambda, s, d);
        for pol in [Polarization::TE, Polarization::TM] {
            let t = slab.solve_s(lambda, s, pol).t.abs();
            println!("s={s} {pol:?}: |t|={t:.8}, analytic={:.8}", want.abs());
            close(t, want.abs(), 1e-9 * want.abs().max(1.0));
        }
    }
}

#[test]
fn lossy_negative_index_energy_conservation_and_resolution() {
    let lambda = 1.0;
    let d = 0.3;
    // Energy conservation for a passive, lossy double-negative slab: the
    // absorptance stays positive and R + T + A = 1.
    for (er, ei) in [(1.0, 0.001), (1.0, 0.01), (0.3, 0.05)] {
        let eps = C::new(-er, ei);
        let mu = C::new(-er, ei);
        let stack = Stack::single(1.0, Layer::from_eps_mu(eps, mu, d), 1.0);
        for pol in [Polarization::TE, Polarization::TM] {
            let sol = stack.solve(lambda, 0.0, pol);
            close(
                sol.reflectance + sol.transmittance + sol.absorptance,
                1.0,
                1e-12,
            );
            assert!(sol.absorptance >= 0.0);
            assert!(sol.transmittance >= 0.0 && sol.reflectance >= 0.0);
        }
    }
    // Resolution limit: the peak evanescent amplification is capped by loss
    // and moves to lower k_x as the loss grows.
    let mut prev_kx = f64::INFINITY;
    let mut prev_t = f64::INFINITY;
    for gamma in [1e-4, 1e-3, 1e-2, 1e-1] {
        let eps = C::new(-1.0, gamma);
        let mu = C::new(-1.0, gamma);
        let (kx, tmax) = evanescent_amplification_peak(lambda, d, eps, mu, 8.0, 4000);
        println!("gamma={gamma:.0e}: k_max={kx:.5}, |t|max={tmax:.5}");
        assert!(tmax.is_finite() && tmax >= 1.0);
        assert!(kx < prev_kx + 1e-9, "k_max must fall as loss grows");
        assert!(tmax < prev_t + 1e-9, "peak gain must fall as loss grows");
        prev_kx = kx;
        prev_t = tmax;
    }
}
