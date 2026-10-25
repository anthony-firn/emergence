//! Sun position and clear-sky irradiance (ASHRAE Handbook — Fundamentals
//! 2017, Ch. 14).
//!
//! Equation of time and declination (eqs. 5, 8), apparent solar time
//! (eq. 7), altitude and azimuth (eqs. 12–14), the ASHRAE clear-sky model
//! with site optical depths τ_b, τ_d (eqs. 17–20, Kasten–Young air mass),
//! and irradiance on a surface: beam, diffuse (the ratio Y for vertical
//! surfaces, eq. 28) and ground-reflected (eq. 31). Angles in degrees at the
//! interface; azimuths from south, west positive; tilt from horizontal.

fn rad(d: f64) -> f64 {
    d.to_radians()
}

#[derive(Clone, Copy, Debug)]
pub struct Site {
    /// Latitude (° north), longitude (° east), standard meridian (° east,
    /// 15 × time zone).
    pub latitude: f64,
    pub longitude: f64,
    pub meridian: f64,
    /// Clear-sky optical depths for the month (ASHRAE climatic data).
    pub tau_b: f64,
    pub tau_d: f64,
    /// Ground reflectance.
    pub ground_reflectance: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Sun {
    /// Altitude above the horizon and azimuth (from south, west +), °.
    pub altitude: f64,
    pub azimuth: f64,
    /// Beam normal and diffuse horizontal irradiance (W/m²).
    pub beam: f64,
    pub diffuse: f64,
}

/// Equation of time (minutes) for day of year `n` (eq. 5).
pub fn equation_of_time(n: f64) -> f64 {
    let g = rad(360.0 * (n - 1.0) / 365.0);
    2.2918 * (0.0075 + 0.1868 * g.cos() - 3.2077 * g.sin() - 1.4615 * (2.0 * g).cos() - 4.089 * (2.0 * g).sin())
}

/// Declination (°), eq. 8.
pub fn declination(n: f64) -> f64 {
    23.45 * rad(360.0 * (n + 284.0) / 365.0).sin()
}

/// Extraterrestrial normal irradiance (W/m²), eq. 4.
pub fn extraterrestrial(n: f64) -> f64 {
    1367.0 * (1.0 + 0.033 * rad(360.0 * (n - 3.0) / 365.0).cos())
}

/// The sun at local standard time `lst` (hours) on day `n`.
pub fn sun(site: &Site, n: f64, lst: f64) -> Sun {
    let ast = lst + equation_of_time(n) / 60.0 + (site.longitude - site.meridian) / 15.0;
    let h = rad(15.0 * (ast - 12.0));
    let (l, d) = (rad(site.latitude), rad(declination(n)));
    let sb = (l.cos() * d.cos() * h.cos() + l.sin() * d.sin()).clamp(-1.0, 1.0);
    let beta = sb.asin();
    // Azimuth from its sine and cosine (eqs. 13–14 over cos β), so it stays
    // well conditioned at noon where arccos would amplify round-off.
    let phi = (h.sin() * d.cos()).atan2(h.cos() * d.cos() * l.sin() - d.sin() * l.cos());
    let alt = beta.to_degrees();
    let (beam, diffuse) = if alt > 0.0 {
        let m = 1.0 / (sb + 0.50572 * (6.07995 + alt).powf(-1.6364));
        let (tb, td) = (site.tau_b, site.tau_d);
        let ab = 1.454 - 0.406 * tb - 0.268 * td + 0.021 * tb * td;
        let ad = 0.507 + 0.205 * tb - 0.080 * td - 0.190 * tb * td;
        let e0 = extraterrestrial(n);
        (e0 * (-tb * m.powf(ab)).exp(), e0 * (-td * m.powf(ad)).exp())
    } else {
        (0.0, 0.0)
    };
    Sun { altitude: alt, azimuth: phi.to_degrees(), beam, diffuse }
}

/// Irradiance on a surface of `azimuth` (from south, west +) and `tilt`
/// (° from horizontal): (beam, diffuse, ground-reflected) W/m², and the
/// incidence angle (°).
pub fn on_surface(site: &Site, s: &Sun, azimuth: f64, tilt: f64) -> (f64, f64, f64, f64) {
    if s.altitude <= 0.0 {
        return (0.0, 0.0, 0.0, 90.0);
    }
    let (b, g, t) = (rad(s.altitude), rad(s.azimuth - azimuth), rad(tilt));
    let cos_theta = b.cos() * g.cos() * t.sin() + b.sin() * t.cos();
    let beam = s.beam * cos_theta.max(0.0);
    // Diffuse: ASHRAE's Y for vertical surfaces, blended for tilt.
    let y = (0.55 + 0.437 * cos_theta + 0.313 * cos_theta * cos_theta).max(0.45);
    let diffuse = if tilt <= 90.0 { s.diffuse * (y * t.sin() + t.cos()) } else { s.diffuse * y * t.sin() };
    let ground = (s.beam * b.sin() + s.diffuse) * site.ground_reflectance * (1.0 - t.cos()) / 2.0;
    (beam, diffuse, ground, cos_theta.clamp(-1.0, 1.0).acos().to_degrees())
}

/// Sol-air temperature (°C): outdoor air plus absorbed solar, less the
/// long-wave deficit for surfaces seeing the sky (ASHRAE: ΔR = 63 W/m² for
/// horizontal, 0 for vertical; blended by tilt), with h_o = 1/R_se.
pub fn sol_air(t_out: f64, absorbed: f64, tilt: f64, h_o: f64) -> f64 {
    let lw = 63.0 * rad(tilt).cos().max(0.0);
    t_out + (absorbed - lw) / h_o
}

/// ASHRAE fraction of the daily range below the maximum, hours 1–24
/// (Fundamentals 2017 Ch. 14 Table 6).
pub const DAILY_RANGE_FRACTION: [f64; 24] = [
    0.88, 0.92, 0.95, 0.98, 1.00, 0.98, 0.91, 0.74, 0.55, 0.38, 0.23, 0.13, 0.05, 0.00, 0.00, 0.06, 0.14, 0.24, 0.39, 0.50, 0.59,
    0.68, 0.75, 0.82,
];

/// Design-day dry bulb (°C) at hour `h` (1–24, fractional allowed).
pub fn design_temperature(t_max: f64, range: f64, h: f64) -> f64 {
    let h = h.rem_euclid(24.0);
    let i = (h.floor() as usize + 23) % 24; // hour h is index h−1
    let j = (i + 1) % 24;
    let f = h - h.floor();
    let frac = DAILY_RANGE_FRACTION[i] * (1.0 - f) + DAILY_RANGE_FRACTION[j] * f;
    t_max - frac * range
}
