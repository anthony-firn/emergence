//! Static elastic and ultrasonic tests, at the single-crystal and
//! polycrystal levels.
//!
//! The single-crystal test is a virtual tensile test: the rotated compliance
//! gives the axial and transverse strains, so `E(n)` and `ν(n)` follow. The
//! reference is an independent tensor contraction of the same `C_ij`, and
//! the handbook polycrystal values are used at the polycrystal rung. The
//! ultrasonic test is the Christoffel eigenproblem for the single crystal
//! and the Hill average for the isotropic polycrystal; the fold is valid
//! only while the polycrystal speeds lie inside the single-crystal bracket
//! and the grain size is in the Rayleigh regime.

use crate::Measurement;
use materials::anisotropic::{
    christoffel_modes, compliance, rotate_stiffness, rotation_align_x, M6,
};
use materials::continuum::christoffel_velocities;
use materials::polycrystal::{
    hill, isotropic_parameters_stiffness, rayleigh_longitudinal_attenuation, scattering_regime,
    ScatteringRegime,
};
use materials::V3;

use crate::Level;

/// The three crystal directions of the single-crystal tensile/ultrasonic tests.
pub const DIRECTIONS: [( &str, V3); 3] = [
    ("[100]", [1.0, 0.0, 0.0]),
    ("[110]", [1.0, 1.0, 0.0]),
    ("[111]", [1.0, 1.0, 1.0]),
];

fn unit(a: V3) -> V3 {
    let n = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    [a[0] / n, a[1] / n, a[2] / n]
}

/// Closed-form cubic directional Young's modulus,
/// `1/E = S11 − 2(S11 − S12 − S44/2)(l²m² + m²n² + n²l²)`.
pub fn cubic_e_closed(c: &M6, combos: f64) -> f64 {
    let s = compliance(c);
    let a = s[0][0] - s[0][1] - 0.5 * s[3][3];
    1.0 / (s[0][0] - 2.0 * a * combos)
}

/// Mean-transverse Poisson ratio via the Bond-rotated compliance: align `n`
/// with x and read `−(S'₁₂ + S'₁₃)/(2 S'₁₁)`.
pub fn poisson_bond(c: &M6, n: V3) -> f64 {
    let r = rotation_align_x(n);
    let s = compliance(&rotate_stiffness(c, &r));
    -(s[0][1] + s[0][2]) / (2.0 * s[0][0])
}

/// Closed-form cubic mean-transverse Poisson ratio. The trace of the strain
/// under uniaxial stress is `S11 + 2 S12` in any direction, so with
/// `S'11 = S11 − 2(S11 − S12 − S44/2)·combos`,
/// `ν = (S'11 − S11 − 2 S12)/(2 S'11)`.
pub fn cubic_poisson_closed(c: &M6, combos: f64) -> f64 {
    let s = compliance(c);
    let a = s[0][0] - s[0][1] - 0.5 * s[3][3];
    let s11 = s[0][0] - 2.0 * a * combos;
    (s11 - s[0][0] - 2.0 * s[0][1]) / (2.0 * s11)
}

/// Static tensile test at every elastic level. Returns, for each direction,
/// the Young's modulus and Poisson ratio (reference: independent tensor
/// contraction), then the polycrystal Hill values (reference: handbook).
pub fn static_elastic(mat: &crate::Material) -> Vec<Measurement> {
    let mut out = Vec::new();
    for (label, n) in DIRECTIONS {
        let n = unit(n);
        let e = mat.young(n);
        // Closed form combos: [100] 0, [110] 1/4, [111] 1/3.
        let combos = match label {
            "[100]" => 0.0,
            "[110]" => 0.25,
            _ => 1.0 / 3.0,
        };
        let e_ref = cubic_e_closed(&mat.c, combos);
        out.push(Measurement::new(
            mat.name,
            "static tensile",
            Level::SingleCrystal,
            &format!("E {label}"),
            e / 1e9,
            e_ref / 1e9,
            "GPa",
            (e / e_ref - 1.0).abs() < 1e-9,
            "virtual tensile test on the single-crystal compliance; reference \
             is the closed-form cubic formula",
        ));
        let nu = poisson_bond(&mat.c, n);
        let nu_ref = cubic_poisson_closed(&mat.c, combos);
        out.push(Measurement::new(
            mat.name,
            "static tensile",
            Level::SingleCrystal,
            &format!("ν {label}"),
            nu,
            nu_ref,
            "",
            (nu - nu_ref).abs() < 1e-9,
            "mean-transverse Poisson ratio via the Bond-rotated compliance; \
             reference is the closed-form cubic formula",
        ));
    }
    // Polycrystal Hill average against the handbook (skipped when no
    // handbook value exists, as for the classical LJ crystal).
    if mat.young_ref.is_finite() && mat.poisson_ref.is_finite() {
        let h = hill(&mat.c);
        out.push(Measurement::new(
            mat.name,
            "static tensile",
            Level::Polycrystal,
            "E (Hill)",
            h.young() / 1e9,
            mat.young_ref / 1e9,
            "GPa",
            (h.young() / mat.young_ref - 1.0).abs() < 0.05,
            "Voigt/Reuss/Hill average of the single-crystal C_ij vs handbook \
             polycrystal value",
        ));
        out.push(Measurement::new(
            mat.name,
            "static tensile",
            Level::Polycrystal,
            "ν (Hill)",
            h.poisson(),
            mat.poisson_ref,
            "",
            (h.poisson() - mat.poisson_ref).abs() < 0.05,
            "Hill average vs handbook",
        ));
    }
    out
}

/// A single ultrasonic result: speeds along one direction.
#[derive(Clone, Copy, Debug)]
pub struct WaveSpeeds {
    pub v_long: f64,
    pub v_shear: f64,
}

/// Christoffel speeds (longitudinal, mean transverse) along a direction.
pub fn christoffel_speeds(mat: &crate::Material, n: V3) -> WaveSpeeds {
    let (v, _) = christoffel_modes(&mat.c, mat.density, n);
    WaveSpeeds {
        v_long: v[2],
        v_shear: 0.5 * (v[0] + v[1]),
    }
}

/// Ultrasonic tests: single-crystal Christoffel speeds, the isotropic
/// polycrystal speeds (Hill), and the single-crystal bracket that certifies
/// the polycrystal fold.
pub fn ultrasonic(mat: &crate::Material) -> Vec<Measurement> {
    let mut out = Vec::new();
    let mut l_min = f64::INFINITY;
    let mut l_max: f64 = 0.0;
    let mut t_min = f64::INFINITY;
    let mut t_max: f64 = 0.0;
    for (label, n) in DIRECTIONS {
        let n = unit(n);
        let ws = christoffel_speeds(mat, n);
        // Independent implementation (`materials::continuum`) as reference.
        let vref = christoffel_velocities(&mat.c, mat.density, n);
        let l_ref = vref[2];
        let t_ref = 0.5 * (vref[0] + vref[1]);
        out.push(Measurement::new(
            mat.name,
            "ultrasonic",
            Level::SingleCrystal,
            &format!("v_L {label}"),
            ws.v_long,
            l_ref,
            "m/s",
            (ws.v_long / l_ref - 1.0).abs() < 1e-9,
            "Christoffel eigenproblem; reference is the independent \
             `materials::continuum` solver",
        ));
        out.push(Measurement::new(
            mat.name,
            "ultrasonic",
            Level::SingleCrystal,
            &format!("v_T {label}"),
            ws.v_shear,
            t_ref,
            "m/s",
            (ws.v_shear / t_ref - 1.0).abs() < 1e-9,
            "mean of the two transverse Christoffel speeds",
        ));
        l_min = l_min.min(ws.v_long);
        l_max = l_max.max(ws.v_long);
        t_min = t_min.min(ws.v_shear);
        t_max = t_max.max(ws.v_shear);
    }

    let ip = isotropic_parameters_stiffness(&mat.c, mat.density);
    let v_l = ip.longitudinal_speed;
    let v_t = ip.shear_sound_speed;
    // Handbook isotropic speeds from E and ν, when a handbook value exists;
    // otherwise the mean of the single-crystal directional speeds.
    let (l_hand, t_hand) = if mat.young_ref.is_finite() && mat.poisson_ref.is_finite() {
        let e = mat.young_ref;
        let nu = mat.poisson_ref;
        (
            (e * (1.0 - nu) / (mat.density * (1.0 + nu) * (1.0 - 2.0 * nu))).sqrt(),
            (e / (2.0 * mat.density * (1.0 + nu))).sqrt(),
        )
    } else {
        (0.5 * (l_min + l_max), 0.5 * (t_min + t_max))
    };
    let l_note = if mat.young_ref.is_finite() {
        "Hill isotropic longitudinal speed vs handbook; fold certified only \
         if it lies inside the single-crystal directional bracket"
    } else {
        "Hill isotropic speed vs the mean single-crystal directional speed \
         (no handbook value for the classical LJ crystal)"
    };
    let t_note = if mat.young_ref.is_finite() {
        "Hill isotropic shear speed vs handbook; fold bracket checked"
    } else {
        "Hill isotropic shear speed vs the mean directional value"
    };
    let l_in = v_l >= l_min && v_l <= l_max;
    let t_in = v_t >= t_min && v_t <= t_max;
    // With a handbook value, require 5 % agreement as well; for the LJ
    // crystal only the directional bracket is meaningful.
    let l_fold = l_in && (mat.young_ref.is_finite().then(|| (v_l / l_hand - 1.0).abs() < 0.05).unwrap_or(true));
    let t_fold = t_in && (mat.young_ref.is_finite().then(|| (v_t / t_hand - 1.0).abs() < 0.05).unwrap_or(true));
    out.push(Measurement::new(
        mat.name,
        "ultrasonic",
        Level::Polycrystal,
        "v_L (Hill)",
        v_l,
        l_hand,
        "m/s",
        l_fold,
        l_note,
    ));
    out.push(Measurement::new(
        mat.name,
        "ultrasonic",
        Level::Polycrystal,
        "v_T (Hill)",
        v_t,
        t_hand,
        "m/s",
        t_fold,
        t_note,
    ));
    out
}

/// The grain-size/wavelength validity indicator of the polycrystal fold.
/// `fold_valid` is true only in the Rayleigh regime, where the isotropic
/// constant-moduli model is valid up to a scattering attenuation.
pub fn scattering_fold(mat: &crate::Material, grain_size: f64, frequency: f64) -> Measurement {
    let ip = isotropic_parameters_stiffness(&mat.c, mat.density);
    let v_l = ip.longitudinal_speed;
    let lambda = v_l / frequency;
    let regime = scattering_regime(grain_size, lambda);
    let alpha = rayleigh_longitudinal_attenuation(
        mat.c11, mat.c12, mat.c44, mat.density, grain_size, frequency,
    );
    let valid = regime == ScatteringRegime::Rayleigh;
    Measurement::new(
        mat.name,
        "ultrasonic",
        Level::Polycrystal,
        "fold regime d/λ",
        grain_size / lambda,
        0.0,
        "",
        valid,
        &format!(
            "grain d = {:.3e} m, f = {:.3e} Hz, λ = {:.3e} m → {:?}; \
             α_L = {:.3e} 1/m (leading Stanke–Kino scaling). The isotropic \
             fold is certified only in the Rayleigh regime.",
            grain_size, frequency, lambda, regime, alpha
        ),
    )
}
