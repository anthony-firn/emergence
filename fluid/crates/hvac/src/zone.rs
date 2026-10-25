//! Design loads of a room (a thermal zone) by heat balance.
//!
//! The room air is held at its set point; the load is what the HVAC system
//! must remove (cooling) or add (heating) to hold it. Every opaque element
//! is a transient finite-volume wall (`envelope::Wall`, combined convective
//! and radiative surface coefficient R_si to the room), driven outside by
//! the sol-air temperature or, for internal partitions and floors on an
//! identical storey, by the room itself. Windows conduct with no storage and
//! transmit SHGC × incident solar. Radiant gains (transmitted solar, the
//! radiant part of people, lights and equipment) are absorbed by the room's
//! surfaces in proportion to their area and reach the air only through
//! them, which is what delays and lowers a cooling peak; convective gains
//! reach the air at once. Outdoor air brings sensible and latent heat.
//!
//! Cooling: the ASHRAE design day (Fundamentals Ch. 14 profile, clear sky)
//! repeated until the response is periodic. Heating: steady state at the
//! winter design temperature, no sun, no internal gains (ASHRAE practice).

use crate::envelope::{Assembly, Wall};
use crate::psychro;
use crate::solar::{self, Site};

/// Latent heat of vaporisation at 0 °C used with humidity ratios (J/kg),
/// as in ASHRAE's enthalpy (eq. 32).
const H_FG: f64 = 2.501e6;
/// Air density and specific heat for airflow loads (ASHRAE standard air:
/// 1.204 kg/m³, 1.006 kJ/(kg K) dry air).
const RHO_AIR: f64 = 1.204;
const CP_AIR: f64 = 1006.0;

#[derive(Clone, Debug)]
pub struct Opaque {
    pub name: String,
    pub assembly: Assembly,
    pub area: f64,
    /// Outside: Some((azimuth, tilt, solar absorptance)) for an exterior
    /// surface; None for an internal one (the same conditions on both sides).
    pub exterior: Option<(f64, f64, f64)>,
}

#[derive(Clone, Debug)]
pub struct Window {
    pub name: String,
    pub area: f64,
    pub u: f64,
    pub shgc: f64,
    pub azimuth: f64,
    pub tilt: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Gains {
    pub people: f64,
    /// Per person (W): sensible and latent.
    pub person_sensible: f64,
    pub person_latent: f64,
    /// Lighting and equipment (W).
    pub lighting: f64,
    pub equipment: f64,
    /// Radiant fractions of people, lighting, equipment.
    pub radiant: (f64, f64, f64),
    /// Fraction of each day the gains are on, as (start hour, end hour).
    pub schedule: (f64, f64),
}

#[derive(Clone, Debug)]
pub struct Zone {
    pub floor_area: f64,
    pub volume: f64,
    pub opaque: Vec<Opaque>,
    pub windows: Vec<Window>,
    pub gains: Gains,
    /// Outdoor air (m³/s) from infiltration and ventilation together.
    pub outdoor_air: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct DesignDay {
    pub day_of_year: f64,
    pub t_max: f64,
    pub range: f64,
    /// Humidity ratio of the outdoor air (kg/kg), held over the day.
    pub w_out: f64,
    pub pressure: f64,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Hour {
    pub hour: f64,
    pub t_out: f64,
    pub sensible: f64,
    pub latent: f64,
    pub solar_transmitted: f64,
    pub envelope: f64,
    pub outdoor_air: f64,
    pub internal: f64,
}

#[derive(Clone, Debug)]
pub struct CoolingLoad {
    pub hours: Vec<Hour>,
    /// Peak sensible and the latent load at that hour (W).
    pub peak_sensible: f64,
    pub peak_latent: f64,
    pub peak_hour: f64,
    /// Days simulated until the response repeated.
    pub days: usize,
}

impl Zone {
    /// Sum of U·A of the exterior envelope, windows included (W/K).
    pub fn ua(&self) -> f64 {
        let walls: f64 = self.opaque.iter().filter(|o| o.exterior.is_some()).map(|o| o.assembly.u_value() * o.area).sum();
        walls + self.windows.iter().map(|w| w.u * w.area).sum::<f64>()
    }

    /// Steady heating load (W) at outdoor `t_out` and indoor `t_in`: no sun
    /// and no internal gains (ASHRAE practice for heating design).
    pub fn heating_load(&self, t_in: f64, t_out: f64) -> f64 {
        (self.ua() + RHO_AIR * CP_AIR * self.outdoor_air) * (t_in - t_out)
    }

    fn interior_area(&self) -> f64 {
        self.opaque.iter().map(|o| o.area).sum()
    }

    /// Cooling load over the design day at set point `t_in` °C and indoor
    /// humidity ratio `w_in`; `steps_per_hour` time steps.
    pub fn cooling_load(&self, site: &Site, day: &DesignDay, t_in: f64, w_in: f64, steps_per_hour: usize) -> CoolingLoad {
        let dt = 3600.0 / steps_per_hour as f64;
        let mut walls: Vec<Wall> = self.opaque.iter().map(|o| if o.exterior.is_some() { Wall::new(&o.assembly, 4, t_in) } else { Wall::new(&half(&o.assembly), 4, t_in) }).collect();
        let area = self.interior_area();
        let g = &self.gains;
        let on = |h: f64| -> bool { h >= g.schedule.0 && h < g.schedule.1 };
        let mut last_day: Vec<Hour> = Vec::new();
        let mut days = 0;
        for d in 0..40 {
            days = d + 1;
            let mut today = Vec::new();
            for k in 1..=24 * steps_per_hour {
                let hour = k as f64 / steps_per_hour as f64; // end of the step, hours 0–24
                let lst = hour - 0.5 / steps_per_hour as f64; // solar at the step's middle
                let t_out = solar::design_temperature(day.t_max, day.range, hour);
                let sun = solar::sun(site, day.day_of_year, lst);
                // Transmitted solar and window conduction.
                let mut solar_in = 0.0;
                let mut window_q = 0.0;
                for w in &self.windows {
                    let (b, df, gr, _) = solar::on_surface(site, &sun, w.azimuth, w.tilt);
                    solar_in += w.shgc * w.area * (b + df + gr);
                    window_q += w.u * w.area * (t_out - t_in);
                }
                let (people, light, equip) = if on(hour) { (g.people * g.person_sensible, g.lighting, g.equipment) } else { (0.0, 0.0, 0.0) };
                let radiant = solar_in + people * g.radiant.0 + light * g.radiant.1 + equip * g.radiant.2;
                let convective = people * (1.0 - g.radiant.0) + light * (1.0 - g.radiant.1) + equip * (1.0 - g.radiant.2);
                let q_abs = radiant / area; // W/m² on every interior surface
                let mut envelope = 0.0;
                for (o, wall) in self.opaque.iter().zip(walls.iter_mut()) {
                    let q = match o.exterior {
                        Some((az, tilt, alpha)) => {
                            let (b, df, gr, _) = solar::on_surface(site, &sun, az, tilt);
                            let t_sa = solar::sol_air(t_out, alpha * (b + df + gr), tilt, 1.0 / o.assembly.r_se);
                            wall.step(dt, t_in, t_sa, q_abs)
                        }
                        // Internal: half the element, adiabatic at its
                        // mid-plane (identical rooms on both sides).
                        None => wall.step(dt, t_in, t_in, q_abs),
                    };
                    envelope += q * o.area;
                }
                let oa_sensible = RHO_AIR * CP_AIR * self.outdoor_air * (t_out - t_in);
                let latent = RHO_AIR * H_FG * self.outdoor_air * (day.w_out - w_in) + if on(hour) { g.people * g.person_latent } else { 0.0 };
                // The surfaces deliver the absorbed radiant heat through
                // `envelope`; the air gets the rest directly.
                let sensible = envelope + window_q + convective + oa_sensible;
                today.push(Hour {
                    hour,
                    t_out,
                    sensible,
                    latent,
                    solar_transmitted: solar_in,
                    envelope,
                    outdoor_air: oa_sensible,
                    internal: convective,
                });
            }
            let converged = !last_day.is_empty()
                && today.iter().zip(&last_day).all(|(a, b)| (a.sensible - b.sensible).abs() < 1e-6 * (1.0 + a.sensible.abs()));
            last_day = today;
            if converged {
                break;
            }
        }
        let peak = last_day.iter().cloned().fold(Hour { sensible: f64::NEG_INFINITY, ..Default::default() }, |a, b| if b.sensible > a.sensible { b } else { a });
        CoolingLoad { hours: last_day, peak_sensible: peak.sensible, peak_latent: peak.latent, peak_hour: peak.hour, days }
    }
}

/// The inner half of a symmetric element, adiabatic at the cut (its far
/// surface resistance is effectively infinite).
pub fn half(a: &Assembly) -> Assembly {
    let total: f64 = a.layers.iter().map(|l| l.thickness).sum();
    let mut left = 0.5 * total;
    let mut layers = Vec::new();
    for l in &a.layers {
        if left <= 0.0 {
            break;
        }
        let t = l.thickness.min(left);
        layers.push(crate::envelope::Layer { thickness: t, ..*l });
        left -= t;
    }
    Assembly { layers, r_si: a.r_si, r_se: 1e12 }
}

/// Supply airflow (m³/s) to meet a sensible load with supply air
/// `delta_t` K below the room.
pub fn supply_airflow(sensible: f64, delta_t: f64) -> f64 {
    sensible / (RHO_AIR * CP_AIR * delta_t)
}

/// Outdoor humidity ratio from design dry bulb and coincident wet bulb.
pub fn humidity_from_wet_bulb(t: f64, twb: f64, p: f64) -> f64 {
    psychro::w_from_wet_bulb(t, twb, p)
}
