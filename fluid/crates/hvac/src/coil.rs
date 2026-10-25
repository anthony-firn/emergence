//! Coils and heat exchangers.
//!
//! * **Cooling and dehumidifying coil** by the apparatus dew point and
//!   bypass factor (ASHRAE Fundamentals; Carrier): the leaving air is the
//!   entering air mixed with saturated air at the ADP in the proportion of
//!   the bypass factor, in temperature and humidity ratio. When the ADP is
//!   above the entering dew point the coil is dry and the humidity ratio
//!   passes unchanged.
//! * **Sensible exchangers** by ε-NTU: counterflow and parallel exact,
//!   crossflow with both streams unmixed by the standard approximation
//!   (Incropera et al., 7th ed., Table 11.3).

use crate::psychro;

#[derive(Clone, Copy, Debug)]
pub struct CoilResult {
    pub t_out: f64,
    pub w_out: f64,
    /// Total, sensible and latent capacity (W), condensate (kg/s).
    pub total: f64,
    pub sensible: f64,
    pub latent: f64,
    pub condensate: f64,
    pub wet: bool,
}

/// Air at `t_in` °C and `w_in` through a coil with apparatus dew point
/// `adp` °C and bypass factor `bf`, dry-air mass flow `m_dot` kg/s.
pub fn cooling_coil(t_in: f64, w_in: f64, adp: f64, bf: f64, m_dot: f64, p: f64) -> CoilResult {
    let w_adp = psychro::w_s(adp, p);
    let wet = w_in > w_adp;
    let t_out = adp + bf * (t_in - adp);
    let w_out = if wet { w_adp + bf * (w_in - w_adp) } else { w_in };
    let (h_in, h_out) = (psychro::enthalpy(t_in, w_in), psychro::enthalpy(t_out, w_out));
    let total = m_dot * (h_in - h_out);
    // Sensible: the enthalpy change at the entering humidity ratio.
    let sensible = m_dot * (psychro::enthalpy(t_in, w_in) - psychro::enthalpy(t_out, w_in));
    CoilResult { t_out, w_out, total, sensible, latent: total - sensible, condensate: m_dot * (w_in - w_out), wet }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Arrangement {
    Counterflow,
    Parallel,
    CrossflowUnmixed,
}

/// Effectiveness at `ntu` and capacity ratio `cr` = C_min/C_max.
pub fn effectiveness(arr: Arrangement, ntu: f64, cr: f64) -> f64 {
    if cr <= 0.0 {
        return 1.0 - (-ntu).exp();
    }
    match arr {
        Arrangement::Counterflow => {
            if (cr - 1.0).abs() < 1e-12 {
                ntu / (1.0 + ntu)
            } else {
                let e = (-ntu * (1.0 - cr)).exp();
                (1.0 - e) / (1.0 - cr * e)
            }
        }
        Arrangement::Parallel => (1.0 - (-ntu * (1.0 + cr)).exp()) / (1.0 + cr),
        Arrangement::CrossflowUnmixed => 1.0 - ((ntu.powf(0.22) / cr) * ((-cr * ntu.powf(0.78)).exp() - 1.0)).exp(),
    }
}
