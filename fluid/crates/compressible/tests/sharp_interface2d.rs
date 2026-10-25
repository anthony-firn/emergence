//! Validation of the 2D sharp solid–fluid interface module.
//!
//! Fast tests here; the converged angle sweep is `#[ignore]`d and can be run
//! with `cargo test -p compressible --test sharp_interface2d -- --ignored
//! --nocapture`. Numbers are collected in `analysis/sharp_interface2d.md`.

use compressible::sharp_interface::P0;
use compressible::sharp_interface2d::*;
use std::f64::consts::PI;

const F0: f64 = 100_000.0;
const NU: usize = 4096;

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

/// Minimal complex arithmetic for the two-probe decomposition.
#[derive(Clone, Copy, Debug)]
struct C {
    re: f64,
    im: f64,
}
impl C {
    fn new(re: f64, im: f64) -> Self {
        C { re, im }
    }
    fn expi(t: f64) -> Self {
        C { re: t.cos(), im: t.sin() }
    }
    fn abs2(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
}
impl std::ops::Add for C {
    type Output = C;
    fn add(self, o: C) -> C {
        C::new(self.re + o.re, self.im + o.im)
    }
}
impl std::ops::Sub for C {
    type Output = C;
    fn sub(self, o: C) -> C {
        C::new(self.re - o.re, self.im - o.im)
    }
}
impl std::ops::Mul for C {
    type Output = C;
    fn mul(self, o: C) -> C {
        C::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    }
}
impl std::ops::Div for C {
    type Output = C;
    fn div(self, o: C) -> C {
        let d = o.abs2();
        C::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }
}
impl std::ops::Neg for C {
    type Output = C;
    fn neg(self) -> C {
        C::new(-self.re, -self.im)
    }
}

// ---------------------------------------------------------------------
// A plane-wave packet in a water/steel half-space
// ---------------------------------------------------------------------

/// Geometry: fluid `[0, y_if]`, solid `[y_if, ly]`, two fluid probes.
struct Setup {
    nx: usize,
    ny: usize,
    lx: f64,
    ly: f64,
    y_if: f64,
    y0: f64,
    width: f64,
    y_probe: [f64; 2],
    amp: f64,
    theta: f64,
}

impl Setup {
    fn c_f() -> f64 {
        water().c_longitudinal()
    }

    /// `theta` is the incidence angle; the x-wavenumber is the single periodic
    /// mode `m = 1`, so `lx = lambda / sin(theta)`.
    fn new(theta: f64, cells_per_lambda: f64) -> Self {
        let c_f = Self::c_f();
        let lam = c_f / F0;
        let lx = if theta.abs() < 1e-12 { 2.0 * lam } else { lam / theta.sin() };
        let dx = lam / cells_per_lambda;
        let nx = (lx / dx).round().max(4.0) as usize;
        // Fluid 7.5 wavelengths deep, solid 12.5 wavelengths deep. The
        // 2-wavelength-wide beam has an angular spread of ~4.5 deg, small
        // enough that no ray leaks backwards at the smallest incidence
        // angle; its centre is 4.5 wavelengths above the floor so the
        // tail at the outflow is < 1%.
        let ly = 20.0 * lam;
        let y_if = 7.5 * lam;
        let ny = (ly / dx).round() as usize;
        Setup {
            nx,
            ny,
            lx: nx as f64 * dx,
            ly,
            y_if,
            y0: 4.5 * lam,
            width: 2.0 * lam,
            y_probe: [y_if - 0.8 * lam, y_if - 2.0 * lam],
            amp: 100.0,
            theta,
        }
    }

    fn build(&self, with_solid: bool) -> Sharp2D {
        let w = water();
        let s = steel();
        let k0 = 2.0 * PI * F0 / Self::c_f();
        let kx = k0 * self.theta.sin();
        let ky = k0 * self.theta.cos();
        let (y_if, y0, width, amp) = (self.y_if, self.y0, self.width, self.amp);
        let phi = move |_x: f64, y: f64| if with_solid { y - y_if } else { -1.0 };
        let fluid = move |x: f64, y: f64| {
            let env = (-(((y - y0) / width).powi(2))).exp();
            let ph = amp * env * (kx * x + ky * y).cos();
            let c = Setup::c_f();
            let om = 2.0 * PI * F0;
            FPrim {
                rho: w.rho0 + ph / (c * c),
                u: kx / (om * w.rho0) * ph,
                v: ky / (om * w.rho0) * ph,
                p: P0 + ph,
            }
        };
        let solid = |_x: f64, _y: f64| SPrim {
            rho: s.rho0,
            u: 0.0,
            v: 0.0,
            sxx: P0,
            syy: P0,
            sxy: 0.0,
        };
        let mut sol = Sharp2D::new(self.nx, self.ny, 3, self.lx, self.ly, [w, s], phi, fluid, solid);
        sol.periodic_x = true;
        sol.acoustic_interface = std::env::var("SHARP2D_ACOUSTIC").is_ok();
        sol
    }

    fn probe_cells(&self) -> [(usize, usize); 2] {
        let dy = self.ly / self.ny as f64;
        let i = (0.31 * self.nx as f64) as usize;
        self.y_probe.map(|y| {
            let j = (y / dy) as usize;
            (i.min(self.nx - 1), j.min(self.ny - 1))
        })
    }
}

/// Run, returning time and the pressure perturbation at the two probes.
fn run_two(setup: &Setup, with_solid: bool, t_end: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut sol = setup.build(with_solid);
    let pr = setup.probe_cells();
    let mut t = Vec::new();
    let (mut p1, mut p2) = (Vec::new(), Vec::new());
    while sol.t < t_end - 1e-15 * t_end.abs() {
        let dt = sol.dt_cfl().min(t_end - sol.t);
        sol.step(dt);
        t.push(sol.t);
        p1.push(sol.fluid_prim(pr[0].0, pr[0].1).p - P0);
        p2.push(sol.fluid_prim(pr[1].0, pr[1].1).p - P0);
    }
    (t, p1, p2)
}

// ---------------------------------------------------------------------
// FFT
// ---------------------------------------------------------------------

fn fft(re: &mut [f64], im: &mut [f64]) {
    let n = re.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            re.swap(i, j);
            im.swap(i, j);
        }
    }
    let mut len = 2;
    while len <= n {
        let ang = -2.0 * PI / len as f64;
        let (wr, wi) = (ang.cos(), ang.sin());
        let mut i = 0;
        while i < n {
            let (mut wlr, mut wli) = (1.0f64, 0.0f64);
            for k in 0..len / 2 {
                let ur = re[i + k];
                let ui = im[i + k];
                let vr = re[i + k + len / 2] * wlr - im[i + k + len / 2] * wli;
                let vi = re[i + k + len / 2] * wli + im[i + k + len / 2] * wlr;
                re[i + k] = ur + vr;
                im[i + k] = ui + vi;
                re[i + k + len / 2] = ur - vr;
                im[i + k + len / 2] = ui - vi;
                let nwr = wlr * wr - wli * wi;
                wli = wlr * wi + wli * wr;
                wlr = nwr;
            }
            i += len;
        }
        len <<= 1;
    }
}

fn resample(t: &[f64], v: &[f64], nu: usize, t_end: f64) -> Vec<f64> {
    let dtu = t_end / nu as f64;
    let mut out = vec![0.0; nu];
    let mut j = 0usize;
    for (m, o) in out.iter_mut().enumerate() {
        let tt = m as f64 * dtu;
        while j + 1 < t.len() && t[j + 1] < tt {
            j += 1;
        }
        *o = if t.is_empty() {
            0.0
        } else if j + 1 >= t.len() {
            *v.last().unwrap()
        } else {
            let f = (tt - t[j]) / (t[j + 1] - t[j]);
            v[j] * (1.0 - f) + v[j + 1] * f
        };
    }
    out
}

/// Run recording the pressure perturbation at the first probe, with a fixed
/// time step so the with/without-solid runs share the same sampling.
fn run_probe(
    setup: &Setup,
    with_solid: bool,
    t_end: f64,
    dt_override: Option<f64>,
) -> (Vec<f64>, Vec<f64>, f64) {
    let mut sol = setup.build(with_solid);
    let pr = setup.probe_cells();
    let dt = dt_override.unwrap_or_else(|| sol.dt_cfl());
    let mut t = Vec::new();
    let mut p = Vec::new();
    while sol.t < t_end - 1e-15 * t_end.abs() {
        let d = dt.min(t_end - sol.t);
        sol.step(d);
        t.push(sol.t);
        p.push(sol.fluid_prim(pr[0].0, pr[0].1).p - P0);
    }
    (t, p, dt)
}

/// Incident-weighted RMS `|R|` error via reference subtraction: run with and
/// without the solid, subtract the two pressure histories at one fluid probe,
/// and form `|reflected|/|incident|` per frequency against the exact
/// half-space P/SV reflection. This mirrors the 1D plate benchmark in
/// `tests/sharp_interface.rs`. Returns `(rms_err, |R|@f0, exact@f0)`.
fn measure_r(setup: &Setup, t_end: f64) -> (f64, f64, f64) {
    let (t1, p1, dt) = run_probe(setup, true, t_end, None);
    let (t2, p2, _) = run_probe(setup, false, t_end, Some(dt));
    let s1 = resample(&t1, &p1, NU, t_end);
    let s2 = resample(&t2, &p2, NU, t_end);
    let mut tr = s1.clone();
    let mut ti = vec![0.0; NU];
    let mut ir = s2.clone();
    let mut ii = vec![0.0; NU];
    fft(&mut tr, &mut ti);
    fft(&mut ir, &mut ii);
    let c = Setup::c_f();
    let f_lo = 0.4 * F0;
    let f_hi = 1.6 * F0;
    let df = 1.0 / t_end;
    let k0 = (f_lo / df).round() as usize;
    let k1 = (f_hi / df).round() as usize;
    let (mut num, mut den) = (0.0f64, 0.0f64);
    let (mut r0m, mut r0e) = (0.0f64, 0.0f64);
    let kx = 2.0 * PI * F0 / c * setup.theta.sin();
    for k in k0..=k1 {
        let f = k as f64 * df;
        let om = 2.0 * PI * f;
        if kx >= om / c {
            continue;
        }
        let im2 = ir[k] * ir[k] + ii[k] * ii[k];
        if im2 < 1e-30 {
            continue;
        }
        let rr = tr[k] - ir[k];
        let ri = ti[k] - ii[k];
        let rmag = ((rr * rr + ri * ri) / im2).sqrt();
        let sin_th = (kx * c / om).clamp(0.0, 1.0);
        let th = sin_th.asin();
        let (er, ei) = halfspace_reflection(
            water().rho0,
            c,
            steel().rho0,
            steel().c_longitudinal(),
            steel().cs,
            th,
        );
        let ee = (er * er + ei * ei).sqrt();
        num += im2 * (rmag - ee).powi(2);
        den += im2;
        if (f - F0).abs() < df {
            r0m = rmag;
            r0e = ee;
            if std::env::var("SHARP2D_REPORT").is_ok() {
                eprintln!(
                    "f0: |inc|={:.4e} |refl|={:.4e} rmag={:.4} exact={:.4}",
                    im2.sqrt(),
                    (rr * rr + ri * ri).sqrt(),
                    rmag,
                    ee
                );
            }
        }
    }
    ((num / den.max(1e-300)).sqrt(), r0m, r0e)
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

/// A uniform fluid at rest stays uniform and its totals are unchanged.
#[test]
fn uniform_fluid_stays_uniform() {
    let w = water();
    let s = steel();
    let mut sol = Sharp2D::new(
        32,
        32,
        3,
        1.0,
        1.0,
        [w, s],
        |_x, _y| -1.0,
        |_x, _y| FPrim { rho: w.rho0, u: 2.0, v: 0.0, p: P0 },
        |_x, _y| SPrim { rho: s.rho0, u: 0.0, v: 0.0, sxx: P0, syy: P0, sxy: 0.0 },
    );
    sol.wall_y = true;
    let before = sol.fluid_totals();
    let dt = sol.dt_cfl();
    for _ in 0..50 {
        sol.step(dt);
    }
    let dev = (0..sol.ny)
        .flat_map(|j| (0..sol.nx).map(move |i| (i, j)))
        .map(|(i, j)| {
            let w = sol.fluid_prim(i, j);
            (w.rho - water().rho0).abs() / water().rho0
                + (w.u - 2.0).abs()
                + w.v.abs()
                + (w.p - P0).abs() / P0
        })
        .fold(0.0f64, f64::max);
    assert!(dev < 1e-10, "deviation {dev:.3e}");
    let after = sol.fluid_totals();
    for k in 0..4 {
        assert!(rel(after[k], before[k]) < 1e-10, "component {k}: {} vs {}", after[k], before[k]);
    }
}

/// Stationary water/steel interface at uniform pressure stays uniform.
#[test]
fn stationary_contact_stays_uniform() {
    let w = water();
    let s = steel();
    let mut sol = Sharp2D::new(
        48,
        48,
        3,
        1.0,
        1.0,
        [w, s],
        |_x, y| y - 0.5,
        |_x, _y| FPrim { rho: w.rho0, u: 0.0, v: 0.0, p: P0 },
        |_x, _y| SPrim { rho: s.rho0, u: 0.0, v: 0.0, sxx: P0, syy: P0, sxy: 0.0 },
    );
    sol.wall_y = true;
    let dt = sol.dt_cfl();
    for _ in 0..300 {
        sol.step(dt);
    }
    let mut dev = 0.0f64;
    for j in 0..sol.ny {
        for i in 0..sol.nx {
            if sol.is_solid(i, j) {
                dev = dev.max((sol.solid_prim(i, j).syy - P0).abs() / P0);
            } else {
                dev = dev.max((sol.fluid_prim(i, j).p - P0).abs() / P0);
            }
        }
    }
    assert!(dev < 1e-6, "max |dp|/p = {dev:.3e}");
}

/// Normal incidence: the 2D solver must recover the 1D Fresnel reflection
/// (the limit of `sharp_interface`).
/// The fluid–solid Riemann solver must reduce to the linear acoustic interface
/// at small amplitude (the solid half of `exact_star_reduces_to_acoustic`).
#[test]
fn fluid_solid_star_matches_acoustic() {
    use compressible::sharp_interface::{acoustic_star, exact_star, steel, water, State};
    let w = water();
    let s = steel();
    let c = w.c_longitudinal();
    for dp in [1.0f64, 100.0, -100.0] {
        let wl = State { rho: w.rho0 + dp / (c * c), u: dp / (w.rho0 * c), p: P0 + dp };
        let wr = State { rho: s.rho0, u: 0.0, p: P0 };
        let ex = exact_star(&w, &wl, &s, &wr);
        let ac = acoustic_star(&w, &wl, &s, &wr);
        assert!((ex.tau - ac.tau).abs() / ac.tau.abs() < 1e-9, "tau {} vs {}", ex.tau, ac.tau);
        assert!((ex.u - ac.u).abs() / ac.u.abs().max(1e-30) < 1e-5, "u {} vs {}", ex.u, ac.u);
    }
}

/// Normal-incidence P/SV reflection coefficient and the two critical angles
/// of a water/steel half-space. Below the longitudinal critical angle both P
/// and S are transmitted; between the two critical angles only S; beyond the
/// shear critical angle the reflection is total (`|R| = 1`).
#[test]
fn halfspace_reflection_reference() {
    let mag = |c: (f64, f64)| (c.0 * c.0 + c.1 * c.1).sqrt();
    let r0 = halfspace_reflection(water().rho0, water().c_longitudinal(), steel().rho0, steel().c_longitudinal(), steel().cs, 0.0);
    let zr = halfspace_reflection(water().rho0, water().c_longitudinal(), steel().rho0, steel().c_longitudinal(), steel().cs, 20.0f64.to_radians());
    let rs = halfspace_reflection(water().rho0, water().c_longitudinal(), steel().rho0, steel().c_longitudinal(), steel().cs, 40.0f64.to_radians());
    assert!((r0.0 - 0.9387).abs() < 5e-4, "normal R = {}", r0.0);
    assert!(mag(zr) < 1.0 && mag(zr) > 0.5, "between critical angles |R| = {}", mag(zr));
    assert!((mag(rs) - 1.0).abs() < 1e-9, "beyond shear critical |R| = {}", mag(rs));
}

/// Converged benchmark (slow): normal-incidence pulse reflection vs the
/// exact Fresnel value. Run with `--ignored --nocapture`.
#[ignore]
#[test]
fn normal_incidence_matches_fresnel() {
    let setup = Setup::new(0.0, 8.0);
    let t_end = 8.0 / (F0);
    let (err, rm, re) = measure_r(&setup, t_end);
    if std::env::var("SHARP2D_REPORT").is_ok() {
        eprintln!("normal incidence: rms|R| err={err:.4e} |R|@f0={rm:.4} exact={re:.4}");
    }
    assert!(err < 0.05, "normal-incidence rms |R| error {err:.4e}");
}

/// Oblique incidence below the longitudinal critical angle.
/// Converged benchmark (slow): oblique reflection below the longitudinal
/// critical angle. Run with `--ignored --nocapture`.
#[ignore]
#[test]
fn oblique_subcritical_matches_psv() {
    let setup = Setup::new(8.0f64.to_radians(), 8.0);
    let t_end = 8.0 / F0;
    let (err, rm, re) = measure_r(&setup, t_end);
    if std::env::var("SHARP2D_REPORT").is_ok() {
        eprintln!("oblique 8 deg: rms|R| err={err:.4e} |R|@f0={rm:.4} exact={re:.4}");
    }
    assert!(err < 0.05, "oblique-subcritical rms |R| error {err:.4e}");
}

/// Oblique incidence between the P and S critical angles, and beyond the
/// shear critical angle where |R| -> 1.
/// Converged benchmark (slow): between the P and S critical angles and beyond
/// the shear one. Run with `--ignored --nocapture`.
#[ignore]
#[test]
fn oblique_supercritical_matches_psv() {
    for theta in [20.0f64, 40.0] {
        let setup = Setup::new(theta.to_radians(), 8.0);
        let t_end = 8.0 / F0;
        let (err, rm, re) = measure_r(&setup, t_end);
        if std::env::var("SHARP2D_REPORT").is_ok() {
            eprintln!("oblique {theta} deg: rms|R| err={err:.4e} |R|@f0={rm:.4} exact={re:.4}");
        }
        assert!(err < 0.06, "oblique {theta} rms |R| error {err:.4e}");
    }
}

/// Manufactured normal-incidence solution: an exact standing P wave in the
/// fluid plus the transmitted P wave in the solid must be preserved by the
/// interface. This bypasses all pulse/FFT extraction and directly measures
/// the relative error in a band around the interface. It is the reliable
/// 1D-recovery / normal-incidence benchmark.
#[test]
fn normal_incidence_manufactured() {
    // Exact normal-incidence standing-wave solution, checked against the 2D
    // solver after a short evolution. This bypasses all pulse/FFT extraction.
    let w = water();
    let s = steel();
    let c = w.c_longitudinal();
    let cp = s.c_longitudinal();
    let lam = s.rho0 * (cp * cp - 2.0 * s.cs * s.cs);
    let mu = s.rho0 * s.cs * s.cs;
    let k = 2.0 * PI * F0 / c;
    let kp = 2.0 * PI * F0 / cp;
    let om = 2.0 * PI * F0;
    let zf = w.rho0 * c;
    let zs = s.rho0 * cp;
    let r = (zs - zf) / (zs + zf);
    let tcoef = 1.0 + r;
    let lam_over = lam / (lam + 2.0 * mu);
    let amp = 100.0;
    let lam_f = c / F0;
    let y_if = 0.05;
    let ly = 0.16;
    let nx = 4;
    let ny = 400;
    let lx = lam_f / 8.0; // kx = 0 field is uniform in x
    let dy = ly / ny as f64;
    let yp = |y: f64| y - y_if;
    let pf = |y: f64, t: f64| {
        let yp = yp(y);
        amp * ((k * yp - om * t).cos() + r * (k * yp + om * t).cos())
    };
    let vf = |y: f64, t: f64| {
        let yp = yp(y);
        amp / zf * ((k * yp - om * t).cos() - r * (k * yp + om * t).cos())
    };
    let pt = |y: f64, t: f64| tcoef * amp * (kp * yp(y) - om * t).cos();
    let vt = |y: f64, t: f64| tcoef * amp / zs * (kp * yp(y) - om * t).cos();
    let mut sol = Sharp2D::new(
        nx,
        ny,
        3,
        lx,
        ly,
        [w, s],
        |_x, y| y - y_if,
        |_x, y| FPrim { rho: w.rho0 + pf(y, 0.0) / (c * c), u: 0.0, v: vf(y, 0.0), p: P0 + pf(y, 0.0) },
        |_x, y| SPrim {
            rho: s.rho0,
            u: 0.0,
            v: vt(y, 0.0),
            sxx: P0 + lam_over * pt(y, 0.0),
            syy: P0 + pt(y, 0.0),
            sxy: 0.0,
        },
    );
    sol.periodic_x = true;
    let t_end = 1.0 / F0;
    while sol.t < t_end {
        let dt = sol.dt_cfl().min(t_end - sol.t);
        sol.step(dt);
    }
    let t = sol.t;
    let mut err_f = 0.0f64;
    let mut err_s = 0.0f64;
    let mut scale = 0.0f64;
    for j in 0..sol.ny {
        let y = (j as f64 + 0.5) * dy;
        if (y - y_if).abs() > 0.5 * lam_f {
            continue; // exclude the outflow boundary region
        }
        if sol.is_solid(2, j) {
            let num = sol.solid_prim(2, j).syy - P0;
            let ex = pt(y, t);
            err_s = err_s.max((num - ex).abs());
            scale = scale.max(ex.abs());
        } else {
            let num = sol.fluid_prim(2, j).p - P0;
            let ex = pf(y, t);
            err_f = err_f.max((num - ex).abs());
            scale = scale.max(ex.abs());
        }
    }
    eprintln!("manufactured normal tF0={:.3} max|err_f|={:.3e} max|err_s|={:.3e} scale={:.1} rel={:.3e}", t * F0, err_f, err_s, scale, err_f.max(err_s) / scale);
    assert!(err_f.max(err_s) / scale < 0.04, "manufactured rel error {:.3e}", err_f.max(err_s) / scale);
}

#[ignore]
#[test]
fn debug_fast_normal() {
    let c = Setup::c_f();
    let lam = c / F0;
    let setup = Setup {
        nx: 4,
        ny: 272,
        lx: lam,
        ly: 17.0 * lam,
        y_if: 5.0 * lam,
        y0: 3.0 * lam,
        width: 1.2 * lam,
        y_probe: [5.0 * lam - 0.5 * lam, 5.0 * lam - 1.5 * lam],
        amp: 100.0,
        theta: 0.0,
    };
    let (err, rm, re) = measure_r(&setup, 6.0 / F0);
    eprintln!("fast normal: err={err:.4e} rmag={rm:.4} exact={re:.4}");
}

#[ignore]
#[test]
fn debug_solid_plane_wave() {
    let s = steel();
    let c = Setup::c_f();
    let lam = s.rho0 * (s.c_longitudinal().powi(2) - 2.0 * s.cs * s.cs);
    let mu = s.rho0 * s.cs * s.cs;
    let z_p = s.rho0 * s.c_longitudinal();
    let lam_over = lam / (lam + 2.0 * mu);
    let lam_w = water();
    let progress = |y: f64| {
        let env = (-(((y - 4.5 * (c / F0)) / (2.0 * c / F0)).powi(2))).exp();
        env * (2.0 * PI * F0 / c * (y - 4.5 * (c / F0))).cos()
    };
    let mut sol = Sharp2D::new(
        4,
        160,
        3,
        2.0 * (c / F0),
        20.0 * (c / F0),
        [lam_w, s],
        |_x, _y| 1.0,
        |_x, _y| FPrim { rho: lam_w.rho0, u: 0.0, v: 0.0, p: P0 },
        |_x, y| {
            let g = 100.0 * progress(y);
            SPrim {
                rho: s.rho0,
                u: 0.0,
                v: g / z_p,
                sxx: P0 + g * lam_over,
                syy: P0 + g,
                sxy: 0.0,
            }
        },
    );
    sol.periodic_x = true;
    let dy = sol.dy;
    let j1 = (6.0 * (c / F0) / dy).round() as usize;
    let j2 = (9.0 * (c / F0) / dy).round() as usize;
    let mut t = Vec::new();
    let mut p1 = Vec::new();
    let mut p2 = Vec::new();
    let t_end = 8.0 / F0;
    while sol.t < t_end {
        let dt = sol.dt_cfl().min(t_end - sol.t);
        sol.step(dt);
        t.push(sol.t);
        p1.push(sol.solid_prim(2, j1).syy - P0);
        p2.push(sol.solid_prim(2, j2).syy - P0);
    }
    let peak = |v: &[f64], a: f64, b: f64| {
        v.iter()
            .enumerate()
            .filter(|(m, _)| t[*m] * F0 >= a && t[*m] * F0 <= b)
            .map(|(_, x)| x.abs())
            .fold(0.0f64, f64::max)
    };
    eprintln!("solid P at y=6: peak={:.3}", peak(&p1, 0.0, 8.0));
    eprintln!("solid P at y=9: peak={:.3}", peak(&p2, 0.0, 8.0));
}

#[ignore]
#[test]
fn debug_fluid_solid_star() {
    use compressible::sharp_interface::{acoustic_star, exact_star, steel, water, State};
    let w = water();
    let s = steel();
    let c = w.c_longitudinal();
    for dp in [1.0f64, 100.0, -100.0, 1.0e6, -1.0e6] {
        let wl = State { rho: w.rho0 + dp / (c * c), u: dp / (w.rho0 * c), p: P0 + dp };
        let wr = State { rho: s.rho0, u: 0.0, p: P0 };
        let ex = exact_star(&w, &wl, &s, &wr);
        let ac = acoustic_star(&w, &wl, &s, &wr);
        eprintln!(
            "dp={dp:e}: exact tau={:.6} u={:.6e} | acoustic tau={:.6} u={:.6e} | rel {:.3e} {:.3e}",
            ex.tau, ex.u, ac.tau, ac.u,
            (ex.tau - ac.tau).abs() / ac.tau.abs(),
            (ex.u - ac.u).abs() / ac.u.abs().max(1e-30)
        );
    }
}

#[ignore]
#[test]
fn debug_interface_pressure() {
    let setup = Setup::new(0.0, 8.0);
    let t_end = 8.0 / F0;
    let dy = setup.ly / setup.ny as f64;
    let i = setup.nx / 2;
    let j1 = (6.0 * (Setup::c_f() / F0) / dy).round() as usize;
    let j2 = (5.0 * (Setup::c_f() / F0) / dy).round() as usize;
    let mut sol = setup.build(true);
    let mut t = Vec::new();
    let mut p1 = Vec::new();
    let mut p2 = Vec::new();
    while sol.t < t_end {
        let dt = sol.dt_cfl().min(t_end - sol.t);
        sol.step(dt);
        t.push(sol.t);
        p1.push(sol.fluid_prim(i, j1).p - P0);
        p2.push(sol.fluid_prim(i, j2).p - P0);
    }
    // peak in incident window and reflected window
    let peak = |v: &[f64], a: f64, b: f64| {
        v.iter()
            .enumerate()
            .filter(|(m, _)| t[*m] * F0 >= a && t[*m] * F0 <= b)
            .map(|(_, x)| x.abs())
            .fold(0.0f64, f64::max)
    };
    for (name, v) in [("y=6.0", &p1), ("y=5.0", &p2)] {
        eprintln!(
            "{name}: inc max={:.3} refl max={:.3} ratio={:.3}",
            peak(v, 0.0, 3.0),
            peak(v, 3.0, 7.5),
            peak(v, 3.0, 7.5) / peak(v, 0.0, 3.0)
        );
    }
}

#[ignore]
#[test]
fn debug_two_probe_reference() {
    let setup = Setup::new(0.0, 8.0);
    let t_end = 8.0 / F0;
    let (t, p1, p2) = run_two(&setup, false, t_end);
    let s1 = resample(&t, &p1, NU, t_end);
    let s2 = resample(&t, &p2, NU, t_end);
    let mut r1 = s1.clone(); let mut i1 = vec![0.0; NU];
    let mut r2 = s2.clone(); let mut i2 = vec![0.0; NU];
    fft(&mut r1, &mut i1); fft(&mut r2, &mut i2);
    let pr = setup.probe_cells();
    let dy = setup.ly / setup.ny as f64;
    let c = Setup::c_f();
    let y1 = (pr[0].1 as f64 + 0.5) * dy;
    let y2 = (pr[1].1 as f64 + 0.5) * dy;
    let df = 1.0 / t_end;
    let k = (F0 / df).round() as usize;
    let om = 2.0 * PI * F0;
    let ky = om / c;
    let a = C::expi(-ky*y1); let b = C::expi(ky*y1);
    let cc = C::expi(-ky*y2); let d = C::expi(ky*y2);
    let det = a*d - b*cc;
    let p1c = C::new(r1[k]/NU as f64, i1[k]/NU as f64);
    let p2c = C::new(r2[k]/NU as f64, i2[k]/NU as f64);
    let pi = (p1c*d - b*p2c)/det;
    let prr = (a*p2c - cc*p1c)/det;
    eprintln!("no-solid f0: |Pi|={:.4} |Pr|={:.4} ratio={:.4}", pi.abs2().sqrt(), prr.abs2().sqrt(), (prr.abs2()/pi.abs2()).sqrt());
    eprintln!("probe y1={:.5} y2={:.5} dy={:.5}", y1, y2, dy);
}

#[ignore]
#[test]
fn debug_normal_time() {
    let mut setup = Setup::new(0.0, 8.0);
    let lam = Setup::c_f() / F0;
    setup.y_probe = [setup.y_if - 0.25 * lam, setup.y_if - 2.0 * lam];
    let t_end = 8.0 / F0;
    let (t, p1w, _p2w) = run_two(&setup, true, t_end);
    let (t2, p1r, _p2r) = run_two(&setup, false, t_end);
    let tg: Vec<f64> = (0..4000).map(|m| m as f64 * t_end / 3999.0).collect();
    let w = resample(&t, &p1w, 4000, t_end);
    let r = resample(&t2, &p1r, 4000, t_end);
    let d: Vec<f64> = w.iter().zip(r.iter()).map(|(a,b)| a-b).collect();
    let mx = |v:&[f64]| v.iter().cloned().fold(0.0f64,|a,b| a.max(b.abs()));
    let amax = |v:&[f64]| { let mut bi=0; let mut bv=0.0; for (i,x) in v.iter().enumerate(){if x.abs()>bv{bv=x.abs();bi=i;}} (tg[bi],bv,v[bi]) };
    eprintln!("inc peak {:?}", amax(&r));
    eprintln!("refl peak {:?}", amax(&d));
    eprintln!("ratio {}", mx(&d)/mx(&r));
    let bi = {
        let mut bi = 0usize; let mut bv=0.0f64;
        for (i,x) in d.iter().enumerate(){ if x.abs()>bv {bv=x.abs();bi=i;} }
        bi
    };
    for m in bi.saturating_sub(20)..(bi+21) {
        eprintln!("  t={:.3e} inc={:+.3} refl={:+.3}", tg[m], r[m], d[m]);
    }
    let i1: f64 = r.iter().map(|x| x*x).sum::<f64>()*t_end/3999.0;
    let i2: f64 = d.iter().map(|x| x*x).sum::<f64>()*t_end/3999.0;
    eprintln!("energy ratio {}", i2/i1);
}

/// Normal incidence must recover the exact 1D solid–fluid Riemann star: the
/// steady interface traction equals `sharp_interface::exact_star`, i.e. the
/// 1D rung's result. This is the 1D-recovery benchmark in 2D.
#[test]
fn normal_incidence_recovers_1d_star() {
    use compressible::sharp_interface::{exact_star, State};
    let w = water();
    let s = steel();
    let ny = 200;
    let ly = 1.0;
    let y_if = 0.5;
    let dp = 100.0;
    let c = w.c_longitudinal();
    let v0 = dp / (w.rho0 * c);
    // exact 1D star for the incident fluid state hitting the solid at rest
    let wl = State { rho: w.rho0 + dp / (c * c), u: v0, p: P0 + dp };
    let wr = State { rho: s.rho0, u: 0.0, p: P0 };
    let star = exact_star(&w, &wl, &s, &wr);
    let mut sol = Sharp2D::new(
        4, ny, 3, 2.0, ly, [w, s],
        |_x, y| y - y_if,
        |_x, _y| FPrim { rho: w.rho0 + dp / (c * c), u: 0.0, v: v0, p: P0 + dp },
        |_x, _y| SPrim { rho: s.rho0, u: 0.0, v: 0.0, sxx: P0, syy: P0, sxy: 0.0 },
    );
    sol.periodic_x = true;
    sol.wall_y = true;
    let jif = (y_if / (ly / ny as f64)) as usize; // first solid cell
    let jf = jif - 2; // fluid cell just below the interface
    let dt = sol.dt_cfl();
    for _ in 0..120 {
        sol.step(dt);
    }
    let f = sol.fluid_prim(1, jf);
    let sd = sol.solid_prim(1, jif);
    let rel = |a: f64, b: f64| (a - b).abs() / b.abs().max(1e-300);
    assert!(rel(f.p, star.tau) < 2e-3, "fluid p {} vs star tau {}", f.p, star.tau);
    assert!(rel(sd.syy, star.tau) < 2e-3, "solid syy {} vs star tau {}", sd.syy, star.tau);
}

#[ignore]
#[test]
fn debug_1d_halfspace() {
    use compressible::sharp_interface::{steel, water, SharpSolver, State};
    let w = water();
    let s = steel();
    let c = w.c_longitudinal();
    let lam = c / F0;
    let n = 1600;
    let (lo, hi) = (0.0, 1.0);
    let x_if = 0.5;
    let x0 = 0.25;
    let width = 1.4 * lam;
    let mk = |with_solid: bool| {
        SharpSolver::new(n, 3, lo, hi, [w, s], move |x| {
            let env = (-(((x - x0) / width).powi(2))).exp();
            let ph = 100.0 * env * ((2.0 * PI * F0 / c) * (x - x0)).cos();
            let m = if with_solid && x >= x_if { s } else { w };
            let matid = if with_solid && x >= x_if { 1 } else { 0 };
            (
                m.to_cons(State { rho: w.rho0 + ph / (c * c), u: ph / (w.rho0 * c), p: P0 + ph }),
                if matid == 1 { x - x_if } else { -1.0 },
            )
        })
    };
    let pi = ((0.42 - lo) / (hi - lo) * n as f64) as usize - 1;
    let run = |mut sol: SharpSolver| -> (Vec<f64>, Vec<f64>) {
        let mut ts = Vec::new();
        let mut ps = Vec::new();
        let mut t = 0.0;
        let t_end = 10.0 / F0;
        while t < t_end {
            let dt = sol.dt_cfl().min(t_end - t);
            sol.step(dt);
            t = sol.t;
            ts.push(t);
            ps.push(sol.prim(pi).p - P0);
        }
        (ts, ps)
    };
    let (t1, p1) = run(mk(true));
    let (t2, p2) = run(mk(false));
    let a = resample(&t1, &p1, 6000, 10.0 / F0);
    let b = resample(&t2, &p2, 6000, 10.0 / F0);
    let d: Vec<f64> = a.iter().zip(b.iter()).map(|(x, y)| x - y).collect();
    let mx = |v: &[f64]| v.iter().cloned().fold(0.0f64, |p, q| p.max(q.abs()));
    eprintln!("1d inc={:.4} refl={:.4} ratio={:.4}", mx(&b), mx(&d), mx(&d) / mx(&b));
}

/// Oblique incidence must be stable: a small 8 deg packet on a coarse grid
/// runs for several periods with every primitive state finite and bounded.
/// This is the regression test for the interface shear-characteristic sign
/// (a flipped sign makes the oblique run diverge within a few periods).
#[test]
fn oblique_incidence_is_stable() {
    let setup = Setup::new(8.0f64.to_radians(), 4.0);
    let mut sol = setup.build(true);
    let t_end = 6.0 / F0;
    let mut maxp = 0.0f64;
    while sol.t < t_end {
        let dt = sol.dt_cfl().min(t_end - sol.t);
        sol.step(dt);
        for j in 0..sol.ny {
            for i in 0..sol.nx {
                if sol.is_solid(i, j) {
                    let w = sol.solid_prim(i, j);
                    assert!(
                        w.u.is_finite() && w.v.is_finite() && w.syy.is_finite(),
                        "solid state diverged at ({i},{j})"
                    );
                    maxp = maxp.max((w.syy - P0).abs());
                } else {
                    let w = sol.fluid_prim(i, j);
                    assert!(
                        w.rho.is_finite() && w.rho > 0.0 && w.p.is_finite(),
                        "fluid state diverged at ({i},{j})"
                    );
                    maxp = maxp.max((w.p - P0).abs());
                }
            }
        }
    }
    assert!(maxp < 1.0e6, "oblique run grew to {maxp:.3e}");
}

/// Closed box energy balance: fluid + solid perturbation energy is conserved
/// by the unsplit flux (the interface is held fixed, so no work is done on
/// it). Reports the fractional drift.
#[test]
fn closed_box_energy_balance() {
    let w = water();
    let s = steel();
    let (lx, ly, nx, ny) = (1.0, 1.0, 48, 96);
    let y_if = 0.5;
    let c = w.c_longitudinal();
    let amp = 100.0;
    let mut sol = Sharp2D::new(
        nx,
        ny,
        3,
        lx,
        ly,
        [w, s],
        |_x, y| y - y_if,
        |_x, y| {
            let env = (-(((y - 0.25) / 0.08).powi(2))).exp();
            let dp = amp * env;
            FPrim { rho: w.rho0 + dp / (c * c), u: 0.0, v: 0.0, p: P0 + dp }
        },
        |_x, _y| SPrim { rho: s.rho0, u: 0.0, v: 0.0, sxx: P0, syy: P0, sxy: 0.0 },
    );
    sol.wall_x = true;
    sol.wall_y = true;
    let (mut nf, mut ns) = (0.0, 0.0);
    for j in 0..ny {
        for i in 0..nx {
            if sol.is_solid(i, j) {
                ns += 1.0;
            } else {
                nf += 1.0;
            }
        }
    }
    let area = sol.dx * sol.dy;
    let fluid_rest = w.rho0 * w.eos.energy(w.rho0, P0) * nf * area;
    let solid_rest = P0 * P0 / (2.0 * (sol.lam + sol.mu)) * ns * area;
    let pert = |sol: &Sharp2D| sol.total_energy() - fluid_rest - solid_rest;
    let e0 = pert(&sol);
    sol.run_to(4.0 / F0);
    let e1 = pert(&sol);
    let drift = (e1 - e0).abs() / e0.abs().max(1e-300);
    eprintln!("closed box: E0={e0:.6e} E1={e1:.6e} rel drift={drift:.3e}");
    assert!(drift < 5e-3, "closed-box energy drift {drift:.3e}");
}

/// Complex exponential `exp(z)`.
fn cexp(z: C) -> C {
    let e = z.re.exp();
    C::new(e * z.im.cos(), e * z.im.sin())
}

/// Exact oblique P/SV half-space fields (interface at `y = 0`, fluid `y < 0`,
/// solid `y > 0`), real part at time `t`. The reflection amplitude `R` and
/// hence the solid P/SV amplitudes `B`, `C_sv` come from the exact
/// Brekhovskikh–Godin impedance `halfspace_reflection`; the unit incident
/// pressure amplitude is `amp`. Used by the manufactured oblique test.
struct ObliqueExact {
    kx: f64,
    ky: f64,
    alpha: C,
    beta: C,
    b: C,
    csv: C,
    r: C,
    kp2: f64,
    lam: f64,
    mu: f64,
    omega: f64,
    rho_f: f64,
    c_f: f64,
    amp: f64,
}

impl ObliqueExact {
    fn new(theta: f64, amp: f64) -> Self {
        let om = 2.0 * PI * F0;
        let c = water().c_longitudinal();
        let kx = om / c * theta.sin();
        Self::with_kx(kx, amp)
    }

    /// Exact fields for a **specific** tangential wavenumber `kx` (used by the
    /// periodic manufactured benchmark so the x-discretisation is exactly
    /// periodic: `kx = 2π m / lx`).
    fn with_kx(kx: f64, amp: f64) -> Self {
        let w = water();
        let s = steel();
        let om = 2.0 * PI * F0;
        let c = w.c_longitudinal();
        let cp = s.c_longitudinal();
        let cs = s.cs;
        let ky = (om * om / (c * c) - kx * kx).max(0.0).sqrt();
        let kp = om / cp;
        let ks = om / cs;
        let csqrt = |x: f64| if x >= 0.0 { C::new(x.sqrt(), 0.0) } else { C::new(0.0, (-x).sqrt()) };
        let alpha = csqrt(kp * kp - kx * kx);
        let beta = csqrt(ks * ks - kx * kx);
        let (rr, ri) = halfspace_reflection(w.rho0, c, s.rho0, cp, cs, (kx * c / om).clamp(0.0, 1.0).asin());
        let r = C::new(rr, ri);
        let lam = s.rho0 * (cp * cp - 2.0 * cs * cs);
        let mu = s.rho0 * cs * cs;
        // sigma_xy(0)=0  =>  C_sv = 2 alpha kx/(kx^2 - beta^2) B
        let c0 = (C::new(2.0 * kx, 0.0) * alpha) / (C::new(kx * kx, 0.0) - beta * beta);
        // sigma_yy(0) = -A(1+R) in the tension-positive convention (the
        // solver stores compression-positive stresses, so the sign is moved
        // here and the returned stresses are negated below)
        let d = C::new(-(lam * kp * kp), 0.0) - C::new(2.0 * mu, 0.0) * alpha * alpha
            + C::new(2.0 * mu * kx, 0.0) * beta * c0;
        let b = -C::new(amp, 0.0) * (C::new(1.0, 0.0) + r) / d;
        let csv = c0 * b;
        Self { kx, ky, alpha, beta, b, csv, r, kp2: kp * kp, lam, mu, omega: om, rho_f: w.rho0, c_f: c, amp }
    }

    fn fluid(&self, x: f64, y: f64, t: f64) -> FPrim {
        let pi = cexp(C::new(0.0, 1.0) * C::new(self.kx * x + self.ky * y - self.omega * t, 0.0));
        let pr = cexp(C::new(0.0, 1.0) * C::new(self.kx * x - self.ky * y - self.omega * t, 0.0));
        let a = C::new(self.amp, 0.0);
        let pt = a * (pi + self.r * pr);
        let iw = C::new(self.rho_f * self.omega, 0.0);
        let vx = (C::new(self.kx, 0.0) * a * (pi + self.r * pr) / iw).re;
        let vy = (C::new(self.ky, 0.0) * a * (pi - self.r * pr) / iw).re;
        let p = pt.re;
        FPrim { rho: self.rho_f + p / (self.c_f * self.c_f), u: vx, v: vy, p: P0 + p }
    }

    fn solid(&self, x: f64, y: f64, t: f64) -> SPrim {
        let i = C::new(0.0, 1.0);
        let tau = C::new(self.kx * x - self.omega * t, 0.0);
        let phi = self.b * cexp(i * (self.alpha * C::new(y, 0.0) + tau));
        let psi = self.csv * cexp(i * (self.beta * C::new(y, 0.0) + tau));
        let kxc = C::new(self.kx, 0.0);
        // displacement potentials; velocities are -i*omega times them
        let ux_c = i * kxc * phi + i * self.beta * psi;
        let uy_c = i * self.alpha * phi - i * kxc * psi;
        let vf = C::new(0.0, -self.omega);
        let ux = (vf * ux_c).re;
        let uy = (vf * uy_c).re;
        let sxx = -(-(C::new(self.lam * self.kp2, 0.0) + C::new(2.0 * self.mu * self.kx * self.kx, 0.0)) * phi
            - C::new(2.0 * self.mu * self.kx, 0.0) * self.beta * psi)
            .re;
        let syy = -(-(C::new(self.lam * self.kp2, 0.0) + C::new(2.0 * self.mu, 0.0) * self.alpha * self.alpha) * phi
            + C::new(2.0 * self.mu * self.kx, 0.0) * self.beta * psi)
            .re;
        let sxy = -(C::new(self.mu, 0.0)
            * (-C::new(2.0 * self.kx, 0.0) * self.alpha * phi + (C::new(self.kx * self.kx, 0.0) - self.beta * self.beta) * psi))
            .re;
        SPrim { rho: steel().rho0, u: ux, v: uy, sxx: P0 + sxx, syy: P0 + syy, sxy }
    }
}

/// Manufactured **oblique** field error: the exact fluid (incident+reflected)
/// and solid (P+SV) plane-wave fields of a water/steel half-space are
/// initialised and evolved for `periods` acoustic periods. Returns
/// `(rel_err, max_abs_err, amplitude_scale, nx, ny)` with the L∞ error taken
/// over a band of `±0.5` fluid wavelengths around the interface, at one
/// interior column. `kx` must be an exact periodic mode `2π m / lx` so the
/// fluid incident wave is periodic in `x`; the box is sized so that outgoing
/// waves cannot return to the measurement band within `periods`. This
/// directly exercises the tangential P↔SV coupling.
fn oblique_field_error(
    kx: f64,
    cells_per_lambda: f64,
    periods: f64,
    amp: f64,
) -> (f64, f64, f64, usize, usize) {
    let w = water();
    let s = steel();
    let lam_f = w.c_longitudinal() / F0;
    let dx = lam_f / cells_per_lambda;
    // Exact periodic tangential mode m = 1: lx = 2π/kx. Normal incidence
    // (kx = 0) is uniform in x and uses a one-wavelength box.
    let lx0 = if kx.abs() < 1e-9 { lam_f } else { 2.0 * PI / kx };
    let nx = (lx0 / dx).round().max(8.0) as usize;
    let lx = nx as f64 * dx;
    let kx = if kx.abs() < 1e-9 { 0.0 } else { 2.0 * PI / lx }; // exactly periodic
    let ex = ObliqueExact::with_kx(kx, amp);
    // Interface with 2 wavelengths of fluid below; solid deep enough that the
    // transmitted P wave (speed cp) cannot reach the top and return.
    let y_if = 2.0 * lam_f;
    let solid_depth = (periods * s.c_longitudinal() / w.c_longitudinal() + 1.5) * lam_f;
    let ly = y_if + solid_depth;
    let ny = (ly / dx).round() as usize;
    let mut sol = Sharp2D::new(
        nx,
        ny,
        3,
        lx,
        ly,
        [w, s],
        |_x, y| y - y_if,
        |x, y| ex.fluid(x, y - y_if, 0.0),
        |x, y| ex.solid(x, y - y_if, 0.0),
    );
    sol.periodic_x = true;
    let t_end = periods / F0;
    while sol.t < t_end - 1e-15 * t_end {
        let dt = sol.dt_cfl().min(t_end - sol.t);
        sol.step(dt);
    }
    let t = sol.t;
    let mut err = 0.0f64;
    let mut scale = 0.0f64;
    let xc = 0.5 * lx;
    let i = ((xc / dx) as usize).min(nx - 1);
    for j in 0..ny {
        let y = (j as f64 + 0.5) * dx - y_if;
        if y.abs() > 0.5 * lam_f {
            continue;
        }
        let (num, sc) = if sol.is_solid(i, j) {
            let w = sol.solid_prim(i, j);
            let e = ex.solid(xc, y, t);
            ((w.syy - e.syy).abs(), e.syy.abs())
        } else {
            let w = sol.fluid_prim(i, j);
            let e = ex.fluid(xc, y, t);
            ((w.p - e.p).abs(), e.p.abs())
        };
        err = err.max(num);
        scale = scale.max(sc - P0);
    }
    let rel = err / scale.max(1e-300);
    (rel, err, scale, nx, ny)
}

/// Fast manufactured-oblique regression at 20 deg (between the P and S
/// critical angles, so the solid carries both a propagating P and SV wave):
/// the exact P/SV field must be preserved to better than 15% after half a
/// period on a 20 cells/λ grid (~10% measured; run time ~14 s). The bound
/// comes from the measured resolution trend (see the `#[ignore]`d
/// `oblique_convergence_sweep`, order ~1.5), not from tuning to one number.
#[test]
fn oblique_manufactured() {
    let kx = 2.0 * PI * F0 / water().c_longitudinal() * 20.0f64.to_radians().sin();
    let (rel, err, scale, nx, ny) = oblique_field_error(kx, 20.0, 0.5, 100.0);
    eprintln!(
        "oblique manufactured 20 deg: {nx}x{ny} max err={err:.3e} scale={scale:.1} rel={rel:.3e}"
    );
    assert!(rel < 0.15, "oblique manufactured rel error {rel:.3e}");
}

/// Grid-convergence sweep of the manufactured oblique field at angles below
/// the longitudinal critical angle (14.4 deg), between the two critical angles
/// (20 deg), beyond the shear critical angle (32 deg, where the reflection is
/// total) and at normal incidence. Slow; run with `--ignored --nocapture`.
#[ignore]
#[test]
fn oblique_convergence_sweep() {
    for &theta_deg in &[0.0f64, 8.0, 20.0, 32.0] {
        let theta = theta_deg.to_radians();
        let kx = 2.0 * PI * F0 / water().c_longitudinal() * theta.sin();
        eprintln!("\n== theta = {theta_deg:.1} deg ==");
        let mut prev: Option<(f64, f64)> = None;
        for &cells in &[8.0f64, 12.0, 16.0] {
            let (rel, err, scale, nx, ny) = oblique_field_error(kx, cells, 1.0, 100.0);
            let order = prev.map(|(pc, pe)| (pe / err.max(1e-300)).ln() / (cells / pc).ln());
            eprintln!(
                "  {cells:>4.0} cells/lam ({nx}x{ny}): rel={rel:.3e} err={err:.3e} scale={scale:.1} order={}",
                order.map(|o| format!("{o:.2}")).unwrap_or_else(|| "-".to_string())
            );
            prev = Some((cells, err));
        }
    }
}

#[test]
fn oblique_exact_interface_check() {
    for theta_deg in [0.0f64, 8.0, 20.0, 40.0] {
        let ex = ObliqueExact::new(theta_deg.to_radians(), 100.0);
        let ef = ex.fluid(0.0, 0.0, 0.0);
        let es = ex.solid(0.0, 0.0, 0.0);
        eprintln!(
            "theta={:.1}: duy={:.3e} dsigma_yy={:.3e} sigma_xy={:.3e} (scale {:.1})",
            theta_deg,
            es.v - ef.v,
            es.syy - ef.p,
            es.sxy,
            (ef.p - P0).abs().max((es.syy - P0).abs()).max(1.0)
        );
        assert!((es.v - ef.v).abs() < 1e-12, "normal velocity jump at {theta_deg} deg");
        assert!((es.syy - ef.p).abs() < 1e-9, "normal traction jump at {theta_deg} deg");
        assert!(es.sxy.abs() < 1e-9, "shear traction at {theta_deg} deg");
    }
}
