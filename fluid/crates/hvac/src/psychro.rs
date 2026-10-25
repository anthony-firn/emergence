//! Moist air (ASHRAE Handbook — Fundamentals 2017, Ch. 1).
//!
//! Saturation pressure of water vapour over liquid water (0–200 °C) and over
//! ice (−100–0 °C) by Hyland & Wexler (1983), ASHRAE eqs. (5) and (6); the
//! ideal-gas relations for humidity ratio (eq. 20), enthalpy (eq. 32) and
//! specific volume (eq. 28); the psychrometric equation for the wet bulb
//! (eqs. 33, 35). The enhancement factor (≈ 1.004 at 1 atm) is neglected, as
//! in those equations. Temperatures in °C where the Handbook uses °C;
//! pressures in Pa; humidity ratio in kg water / kg dry air; enthalpy in
//! J/kg dry air.

/// Standard atmosphere at sea level (Pa).
pub const P_ATM: f64 = 101_325.0;
/// Ratio of molar masses, water / dry air (ASHRAE 2017 eq. 20).
const EPS: f64 = 0.621945;

/// Saturation pressure (Pa) at `t` °C, over ice below 0 °C.
pub fn p_ws(t: f64) -> f64 {
    let tk = t + 273.15;
    if t < 0.0 {
        let c = [-5.6745359e3, 6.3925247, -9.6778430e-3, 6.2215701e-7, 2.0747825e-9, -9.4840240e-13, 4.1635019];
        (c[0] / tk + c[1] + c[2] * tk + c[3] * tk * tk + c[4] * tk.powi(3) + c[5] * tk.powi(4) + c[6] * tk.ln()).exp()
    } else {
        let c = [-5.8002206e3, 1.3914993, -4.8640239e-2, 4.1764768e-5, -1.4452093e-8, 6.5459673];
        (c[0] / tk + c[1] + c[2] * tk + c[3] * tk * tk + c[4] * tk.powi(3) + c[5] * tk.ln()).exp()
    }
}

/// Humidity ratio from vapour partial pressure.
pub fn w_from_pw(pw: f64, p: f64) -> f64 {
    EPS * pw / (p - pw)
}

pub fn pw_from_w(w: f64, p: f64) -> f64 {
    p * w / (EPS + w)
}

/// Humidity ratio at saturation.
pub fn w_s(t: f64, p: f64) -> f64 {
    w_from_pw(p_ws(t), p)
}

/// Enthalpy (J/kg dry air) of moist air at `t` °C, humidity ratio `w`.
pub fn enthalpy(t: f64, w: f64) -> f64 {
    1000.0 * (1.006 * t + w * (2501.0 + 1.86 * t))
}

/// Dry-bulb temperature (°C) from enthalpy and humidity ratio.
pub fn t_from_h(h: f64, w: f64) -> f64 {
    (h / 1000.0 - 2501.0 * w) / (1.006 + 1.86 * w)
}

/// Specific volume (m³/kg dry air).
pub fn volume(t: f64, w: f64, p: f64) -> f64 {
    0.287042 * (t + 273.15) * (1.0 + 1.607858 * w) / (p / 1000.0)
}

/// Humidity ratio from dry bulb and relative humidity (0–1).
pub fn w_from_rh(t: f64, rh: f64, p: f64) -> f64 {
    w_from_pw(rh * p_ws(t), p)
}

pub fn rh(t: f64, w: f64, p: f64) -> f64 {
    pw_from_w(w, p) / p_ws(t)
}

fn bisect(f: impl Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> f64 {
    let flo = f(lo);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if (f(mid) > 0.0) == (flo > 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-12 {
            break;
        }
    }
    0.5 * (lo + hi)
}

/// Dew point (°C): where saturation pressure equals the vapour pressure
/// (frost point below 0 °C).
pub fn dew_point(w: f64, p: f64) -> f64 {
    let pw = pw_from_w(w, p);
    bisect(|t| p_ws(t) - pw, -100.0, 200.0)
}

/// Humidity ratio from dry bulb and thermodynamic wet bulb (ASHRAE eqs. 33
/// and 35).
pub fn w_from_wet_bulb(t: f64, twb: f64, p: f64) -> f64 {
    let ws = w_s(twb, p);
    if twb >= 0.0 {
        ((2501.0 - 2.326 * twb) * ws - 1.006 * (t - twb)) / (2501.0 + 1.86 * t - 4.186 * twb)
    } else {
        ((2830.0 - 0.24 * twb) * ws - 1.006 * (t - twb)) / (2830.0 + 1.86 * t - 2.1 * twb)
    }
}

/// Thermodynamic wet bulb (°C) by inverting the psychrometric equation.
pub fn wet_bulb(t: f64, w: f64, p: f64) -> f64 {
    bisect(|twb| w_from_wet_bulb(t, twb, p) - w, -100.0, t)
}

/// Standard atmospheric pressure at altitude `z` m (ASHRAE 2017 eq. 3).
pub fn pressure_at(z: f64) -> f64 {
    101_325.0 * (1.0 - 2.25577e-5 * z).powf(5.2559)
}
