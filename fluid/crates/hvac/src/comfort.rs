//! Thermal comfort: Fanger's predicted mean vote and predicted percentage
//! dissatisfied, exactly as ISO 7730:2005 Annex D gives the computation
//! (also ASHRAE 55).

/// PMV and PPD (%).
///
/// `ta` air and `tr` mean radiant temperature (°C), `var` relative air
/// speed (m/s), `rh` relative humidity (0–1), `met` metabolic rate (met),
/// `clo` clothing (clo), `wme` external work (met).
pub fn pmv_ppd(ta: f64, tr: f64, var: f64, rh: f64, met: f64, clo: f64, wme: f64) -> (f64, f64) {
    // Water vapour pressure (Pa), the standard's own saturation fit.
    let pa = rh * 1000.0 * (16.6536 - 4030.183 / (ta + 235.0)).exp();
    let icl = 0.155 * clo;
    let m = met * 58.15;
    let w = wme * 58.15;
    let mw = m - w;
    let fcl = if icl <= 0.078 { 1.0 + 1.29 * icl } else { 1.05 + 0.645 * icl };
    let hcf = 12.1 * var.sqrt();
    let taa = ta + 273.0;
    let tra = tr + 273.0;
    let tcla = taa + (35.5 - ta) / (3.5 * icl + 0.1);
    let p1 = icl * fcl;
    let p2 = p1 * 3.96;
    let p3 = p1 * 100.0;
    let p4 = p1 * taa;
    let p5 = 308.7 - 0.028 * mw + p2 * (tra / 100.0).powi(4);
    let mut xn = tcla / 100.0;
    let mut xf = tcla / 50.0;
    let mut hc = hcf;
    for _ in 0..150 {
        xf = (xf + xn) / 2.0;
        let hcn = 2.38 * (100.0 * xf - taa).abs().powf(0.25);
        hc = hcf.max(hcn);
        xn = (p5 + p4 * hc - p2 * xf.powi(4)) / (100.0 + p3 * hc);
        if (xn - xf).abs() <= 0.00015 {
            break;
        }
    }
    let tcl = 100.0 * xn - 273.0;
    let hl1 = 3.05 * 0.001 * (5733.0 - 6.99 * mw - pa);
    let hl2 = if mw > 58.15 { 0.42 * (mw - 58.15) } else { 0.0 };
    let hl3 = 1.7e-5 * m * (5867.0 - pa);
    let hl4 = 0.0014 * m * (34.0 - ta);
    let hl5 = 3.96 * fcl * (xn.powi(4) - (tra / 100.0).powi(4));
    let hl6 = fcl * hc * (tcl - ta);
    let ts = 0.303 * (-0.036 * m).exp() + 0.028;
    let pmv = ts * (mw - hl1 - hl2 - hl3 - hl4 - hl5 - hl6);
    let ppd = 100.0 - 95.0 * (-0.03353 * pmv.powi(4) - 0.2179 * pmv * pmv).exp();
    (pmv, ppd)
}

/// Operative temperature (°C): the mean of air and radiant temperatures
/// weighted by the convective and radiative coefficients; ISO 7730 allows
/// the plain mean below 0.2 m/s.
pub fn operative(ta: f64, tr: f64, var: f64) -> f64 {
    let a = if var < 0.2 { 0.5 } else if var < 0.6 { 0.6 } else { 0.7 };
    a * ta + (1.0 - a) * tr
}
