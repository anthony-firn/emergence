//! Plane-parallel Eddington grey LTE atmosphere; no spectral lines or spots.
use radiation::optics::{cie_xyz, planck, xyz_to_rgb};
use std::f64::consts::PI;
pub const SIGMA: f64 = 5.670_374_419e-8;
// 32-point Gauss-Laguerre quadrature of integral exp(-u) f(u) du.
// Generated with numpy.polynomial.laguerre.laggauss(32).
const QUADRATURE: [(f64, f64); 32] = [
    (0.04448936583326739, 0.10921834195242515),
    (0.23452610951961825, 0.2104431079387924),
    (0.5768846293018863, 0.23521322966984298),
    (1.072448753817818, 0.19590333597287354),
    (1.7224087764446456, 0.1299837862860684),
    (2.5283367064257942, 0.07057862386571556),
    (3.492213273021994, 0.03176091250917397),
    (4.616456769749767, 0.011918214834838204),
    (5.903958504174244, 0.003738816294611417),
    (7.358126733186241, 0.0009808033066149246),
    (8.982940924212597, 0.0002148649188013578),
    (10.78301863253997, 3.920341967987801e-05),
    (12.763697986742725, 5.934541612868448e-06),
    (14.931139755522558, 7.41640457866734e-07),
    (17.292454336715316, 7.604567879120552e-08),
    (19.855860940336054, 6.350602226625618e-09),
    (22.630889013196775, 4.281382971040738e-10),
    (25.628636022459247, 2.305899491891237e-11),
    (28.862101816323474, 9.799379288726772e-13),
    (32.346629153964734, 3.2378016577291567e-14),
    (36.10049480575197, 8.171823443420549e-16),
    (40.14571977153944, 1.542133833393811e-17),
    (44.50920799575494, 2.119792290163547e-19),
    (49.22439498730864, 2.0544296737880036e-21),
    (54.33372133339691, 1.3469825866373617e-23),
    (59.89250916213402, 5.661294130397462e-26),
    (65.97537728793505, 1.4185605454629684e-28),
    (72.68762809066271, 1.9133754944539943e-31),
    (80.18744697791352, 1.1922487600982012e-34),
    (88.7353404178924, 2.6715112192396625e-38),
    (98.82954286828398, 1.3386169421064243e-42),
    (111.7513980979377, 4.5105361938987784e-48),
];

fn validate(mu: f64, effective_temperature: f64) -> Result<(), &'static str> {
    if !(0.0..=1.0).contains(&mu)
        || !effective_temperature.is_finite()
        || effective_temperature <= 0.0
    {
        return Err("invalid photosphere direction/temperature");
    }
    Ok(())
}
/// Local source temperature T^4 = 3/4 Teff^4 (tau + 2/3).
pub fn temperature(effective_temperature: f64, optical_depth: f64) -> Result<f64, &'static str> {
    validate(1.0, effective_temperature)?;
    if !optical_depth.is_finite() || optical_depth < 0.0 {
        return Err("invalid atmosphere optical depth");
    }
    let t = effective_temperature * (0.75 * optical_depth + 0.5).sqrt().sqrt();
    if !t.is_finite() {
        return Err("unrepresentable atmosphere temperature");
    }
    Ok(t)
}
/// Analytic bolometric intensity [W/m²/sr]. The limiting limb ray has mu=0.
/// Integrating 2*pi*I(mu)*mu dmu returns sigma*Teff^4.
pub fn bolometric(mu: f64, effective_temperature: f64) -> Result<f64, &'static str> {
    validate(mu, effective_temperature)?;
    let value = 3.0 * SIGMA * effective_temperature.powi(4) / (4.0 * PI) * (mu + 2.0 / 3.0);
    if !value.is_finite() {
        return Err("unrepresentable photosphere radiance");
    }
    Ok(value)
}
/// Formal outward transfer integral [W/m²/sr/m], with tau=mu*u:
/// I_lambda(mu) = integral exp(-u) B_lambda(T(mu*u)) du.
/// Fixed quadrature targets the visible continuum; it is not a general extreme
/// UV/temperature accuracy guarantee. There is no invented coronal emission.
pub fn spectral(mu: f64, wavelength: f64, effective_temperature: f64) -> Result<f64, &'static str> {
    validate(mu, effective_temperature)?;
    if !wavelength.is_finite() || wavelength <= 0.0 {
        return Err("invalid wavelength");
    }
    let mut intensity = 0.0;
    for (u, w) in QUADRATURE {
        intensity += w * planck(wavelength, temperature(effective_temperature, mu * u)?);
    }
    if !intensity.is_finite() {
        return Err("unrepresentable spectral radiance");
    }
    Ok(intensity)
}
/// Visible spectral radiance integrated through the existing CIE response.
/// Linear sRGB primaries; exposure/tone mapping is a later display operation.
pub fn rgb(mu: f64, effective_temperature: f64) -> Result<[f64; 3], &'static str> {
    let mut xyz = [0.0; 3];
    let h = 5e-9;
    for k in 0..=80 {
        let wavelength = 380e-9 + k as f64 * h;
        let weight = if k == 0 || k == 80 { 0.5 } else { 1.0 };
        let intensity = spectral(mu, wavelength, effective_temperature)?;
        for (sum, response) in xyz.iter_mut().zip(cie_xyz(wavelength)) {
            *sum += weight * h * intensity * response;
        }
    }
    Ok(xyz_to_rgb(xyz))
}
