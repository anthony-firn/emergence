//! When does the continuum (Navier–Stokes with these closures) hold?
//!
//! Two independent dimensionless indicators decide:
//!
//! 1. **Knudsen number** `Kn = lambda / L`, with the mean free path
//!    estimated from the self-diffusion coefficient, `lambda = 3 D / v_mean`
//!    and `v_mean = sqrt(8 T/(pi m))` (reduced units). Navier–Stokes with
//!    no-slip boundary conditions requires `Kn <~ 0.1`; below that the
//!    constitutive law is a small correction to the stress.
//! 2. **Scale separation** `L / (c_s tau)`, the ratio of the flow length to
//!    the distance over which a stress/heat-flux correlation persists,
//!    `ell = c_s tau`. This is the ratio the autocorrelation decay
//!    measures directly. Continuum behaviour needs `L >> c_s tau`.
//!
//! The sound speed is estimated from the ideal-gas value with `gamma = 5/3`,
//! `c_s = sqrt(5 T/(3 m))`, which is within a few per cent of the LJ value
//! over the fluid range and is deliberately cheap (it is an indicator, not
//! a prediction). The validity test reports both numbers and a combined
//! boolean.

/// Flow scales supplied to the indicator.
#[derive(Clone, Copy, Debug)]
pub struct FlowScales {
    /// Characteristic flow length.
    pub length: f64,
    /// Characteristic flow velocity.
    pub velocity: f64,
    /// Temperature.
    pub temperature: f64,
    /// Molecular mass.
    pub mass: f64,
}

/// Result of the continuum-validity check.
#[derive(Clone, Copy, Debug)]
pub struct ContinuumIndicator {
    /// `Kn = lambda / L`.
    pub knudsen: f64,
    /// Mean free path `lambda = 3 D / v_mean`.
    pub mean_free_path: f64,
    /// Sound-speed estimate `c_s`.
    pub sound_speed: f64,
    /// Correlation length `ell = c_s tau`.
    pub correlation_length: f64,
    /// `L / (c_s tau)`.
    pub length_ratio: f64,
    /// `(L/U) / tau`, the same separation expressed as a time ratio.
    pub time_ratio: f64,
}

impl ContinuumIndicator {
    /// Heuristic continuum criterion: `Kn < 0.1` and at least a factor 10
    /// of scale separation.
    pub fn holds(&self) -> bool {
        self.knudsen < 0.1 && self.length_ratio > 10.0
    }
}

/// Build the indicator from the self-diffusion coefficient `d`, a stress or
/// heat-flux correlation time `tau`, and the flow scales.
pub fn continuum_indicator(flow: FlowScales, d: f64, tau: f64) -> ContinuumIndicator {
    let v_mean = (8.0 * flow.temperature / (std::f64::consts::PI * flow.mass)).sqrt();
    let mean_free_path = 3.0 * d / v_mean;
    let sound_speed = (5.0 * flow.temperature / (3.0 * flow.mass)).sqrt();
    let correlation_length = sound_speed * tau;
    let length_ratio = if correlation_length > 0.0 {
        flow.length / correlation_length
    } else {
        f64::INFINITY
    };
    let time_ratio = if tau > 0.0 {
        (flow.length / flow.velocity.max(f64::MIN_POSITIVE)) / tau
    } else {
        f64::INFINITY
    };
    ContinuumIndicator {
        knudsen: mean_free_path / flow.length,
        mean_free_path,
        sound_speed,
        correlation_length,
        length_ratio,
        time_ratio,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dense_liquid_is_continuum_and_rarefied_gas_is_not() {
        let flow = FlowScales {
            length: 100.0,
            velocity: 1.0,
            temperature: 1.0,
            mass: 1.0,
        };
        // Dense liquid: short mean free path, short correlation time.
        let liq = continuum_indicator(flow, 0.03, 1.0);
        assert!(liq.holds(), "liquid indicator {liq:?}");
        // Dilute gas: D is large, so the mean free path is long.
        let gas = continuum_indicator(flow, 20.0, 1.0);
        assert!(!gas.holds(), "gas indicator {gas:?}");
        assert!(gas.knudsen > 0.1);
    }
}
