//! Sun, sky and a room's design loads.

use hvac::envelope::{Assembly, Flow, Layer};
use hvac::solar::{self, Site};
use hvac::zone::*;

fn site() -> Site {
    // 40° N, standard meridian 75° W, a clear summer sky (optical depths of
    // the order ASHRAE tabulates for July at mid-latitude sites).
    Site { latitude: 40.0, longitude: -75.0, meridian: -75.0, tau_b: 0.45, tau_d: 2.2, ground_reflectance: 0.2 }
}

#[test]
fn solar_geometry_identities() {
    let s = site();
    for n in [80.0, 172.0, 264.0, 355.0] {
        let dec = solar::declination(n);
        // At apparent solar noon the altitude is 90° − L + δ.
        let noon = 12.0 - solar::equation_of_time(n) / 60.0;
        let sun = solar::sun(&s, n, noon);
        assert!((sun.altitude - (90.0 - s.latitude + dec)).abs() < 1e-9, "day {n}: {} vs {}", sun.altitude, 90.0 - s.latitude + dec);
        assert!(sun.azimuth.abs() < 1e-6);
        // Sunrise where cos H = −tan L tan δ: altitude zero there.
        let hs = (-(s.latitude.to_radians().tan()) * dec.to_radians().tan()).acos().to_degrees();
        let rise = noon - hs / 15.0;
        assert!(solar::sun(&s, n, rise).altitude.abs() < 1e-6);
    }
    // Cooper's declination (eq. 8, what ASHRAE uses) against Spencer's
    // more accurate Fourier series: a statement of the formula's accuracy.
    let mut worst: f64 = 0.0;
    for n in 1..=365 {
        let g = 2.0 * std::f64::consts::PI * (n as f64 - 1.0) / 365.0;
        let spencer = (0.006918 - 0.399912 * g.cos() + 0.070257 * g.sin() - 0.006758 * (2.0 * g).cos() + 0.000907 * (2.0 * g).sin()
            - 0.002697 * (3.0 * g).cos()
            + 0.00148 * (3.0 * g).sin())
        .to_degrees();
        worst = worst.max((solar::declination(n as f64) - spencer).abs());
    }
    println!("Cooper vs Spencer declination: worst {worst:.2}°");
    assert!(worst < 2.0, "{worst}");
    // A south wall at noon in summer gets less beam than the roof.
    let sun = solar::sun(&s, 172.0, 12.0);
    let (roof, _, _, _) = solar::on_surface(&s, &sun, 0.0, 0.0);
    let (wall, _, _, _) = solar::on_surface(&s, &sun, 0.0, 90.0);
    assert!(roof > wall && wall > 0.0);
    // Clear-sky beam normal of a few hundred to ~900 W/m² at that altitude.
    println!("21 June noon at 40° N: beam normal {:.0} W/m², diffuse horizontal {:.0} W/m²", sun.beam, sun.diffuse);
    assert!(sun.beam > 600.0 && sun.beam < 1000.0);
}

fn polar_night() -> Site {
    Site { latitude: -89.0, ..site() }
}

fn room() -> Zone {
    let wall = Assembly::new(
        vec![
            Layer { thickness: 0.013, k: 0.57, rho: 1300.0, cp: 1000.0 },
            Layer { thickness: 0.10, k: 0.77, rho: 1700.0, cp: 1000.0 },
            Layer { thickness: 0.08, k: 0.035, rho: 30.0, cp: 1450.0 },
            Layer { thickness: 0.10, k: 0.77, rho: 1700.0, cp: 1000.0 },
        ],
        Flow::Horizontal,
    );
    let slab = Assembly::new(vec![Layer { thickness: 0.2, k: 1.4, rho: 2300.0, cp: 900.0 }], Flow::Downward);
    let partition = Assembly::new(vec![Layer { thickness: 0.1, k: 0.25, rho: 800.0, cp: 1000.0 }], Flow::Horizontal);
    Zone {
        floor_area: 30.0,
        volume: 81.0,
        opaque: vec![
            Opaque { name: "south wall".into(), assembly: wall.clone(), area: 13.5, exterior: Some((0.0, 90.0, 0.6)) },
            Opaque { name: "west wall".into(), assembly: wall, area: 15.0, exterior: Some((90.0, 90.0, 0.6)) },
            Opaque { name: "floor".into(), assembly: slab.clone(), area: 30.0, exterior: None },
            Opaque { name: "ceiling".into(), assembly: slab, area: 30.0, exterior: None },
            Opaque { name: "partitions".into(), assembly: partition, area: 27.0, exterior: None },
        ],
        windows: vec![Window { name: "south glazing".into(), area: 6.0, u: 1.6, shgc: 0.4, azimuth: 0.0, tilt: 90.0 }],
        gains: Gains {
            people: 3.0,
            person_sensible: 75.0,
            person_latent: 55.0,
            lighting: 300.0,
            equipment: 450.0,
            radiant: (0.6, 0.67, 0.3),
            schedule: (8.0, 18.0),
        },
        outdoor_air: 0.03,
    }
}

#[test]
fn steady_state_is_ua_delta_t() {
    // No sun (polar night: 89° S in July), no daily range, no gains: the
    // load is exactly the steady heat balance. (A huge optical depth does
    // not darken the ASHRAE sky: the diffuse exponent turns negative.)
    let s = polar_night();
    let mut z = room();
    z.gains.schedule = (0.0, 0.0);
    let day = DesignDay { day_of_year: 200.0, t_max: 32.0, range: 0.0, w_out: 0.012, pressure: 101_325.0 };
    let c = z.cooling_load(&s, &day, 24.0, 0.009, 4);
    let steady = (z.ua() + 1.204 * 1006.0 * z.outdoor_air) * (32.0 - 24.0);
    // Sol-air of the vertical walls equals the air (no sky deficit, no sun).
    println!("steady: {:.4} W vs UA·ΔT {:.4} W ({} days)", c.peak_sensible, steady, c.days);
    assert!((c.peak_sensible / steady - 1.0).abs() < 1e-6);
    assert!((z.heating_load(21.0, -10.0) - (z.ua() + 1.204 * 1006.0 * z.outdoor_air) * 31.0).abs() < 1e-9);
}

#[test]
fn daily_mean_is_the_steady_balance_of_daily_means() {
    // The model is linear and periodic: storage cancels over a day, so the
    // mean load equals the steady balance of the mean drivers, including
    // the share of absorbed radiant heat that leaks out through exterior
    // walls (1 − U·R_si of it returns to the room).
    let s = site();
    let z = room();
    let day = DesignDay { day_of_year: 200.0, t_max: 33.0, range: 11.0, w_out: 0.013, pressure: 101_325.0 };
    let steps = 4;
    let c = z.cooling_load(&s, &day, 24.0, 0.0093, steps);
    let n = c.hours.len() as f64;
    let mean = |f: &dyn Fn(&Hour) -> f64| c.hours.iter().map(f).sum::<f64>() / n;
    let area: f64 = z.opaque.iter().map(|o| o.area).sum();
    // Mean absorbed radiant flux and mean sol-air temperatures, recomputed.
    let g = &z.gains;
    let mut q_rad = 0.0;
    let mut ext = 0.0;
    for k in 1..=24 * steps {
        let hour = k as f64 / steps as f64;
        let lst = hour - 0.5 / steps as f64;
        let t_out = hvac::solar::design_temperature(day.t_max, day.range, hour);
        let sun = solar::sun(&s, day.day_of_year, lst);
        let on = hour >= g.schedule.0 && hour < g.schedule.1;
        let mut rad = if on { g.people * g.person_sensible * g.radiant.0 + g.lighting * g.radiant.1 + g.equipment * g.radiant.2 } else { 0.0 };
        for w in &z.windows {
            let (b, d, gr, _) = solar::on_surface(&s, &sun, w.azimuth, w.tilt);
            rad += w.shgc * w.area * (b + d + gr);
        }
        q_rad += rad / area;
        for o in &z.opaque {
            if let Some((az, tilt, alpha)) = o.exterior {
                let (b, d, gr, _) = solar::on_surface(&s, &sun, az, tilt);
                let tsa = solar::sol_air(t_out, alpha * (b + d + gr), tilt, 1.0 / o.assembly.r_se);
                ext += o.assembly.u_value() * o.area * (tsa - 24.0);
            }
        }
    }
    let (q_rad, ext) = (q_rad / n, ext / n);
    let returned: f64 = z
        .opaque
        .iter()
        .map(|o| if o.exterior.is_some() { o.area * q_rad * (1.0 - o.assembly.u_value() * o.assembly.r_si) } else { o.area * q_rad })
        .sum();
    let predicted = ext + returned + mean(&|h| h.sensible - h.envelope);
    let got = mean(&|h| h.sensible);
    println!("daily mean sensible {got:.4} W, steady balance of means {predicted:.4} W ({:+.1e}); peak {:.0} W sensible + {:.0} W latent at {:.2} h", got / predicted - 1.0, c.peak_sensible, c.peak_latent, c.peak_hour);
    assert!((got / predicted - 1.0).abs() < 1e-6);
}

#[test]
fn radiant_gains_arrive_later_and_lower() {
    // The same 1 kW of equipment, all radiant against all convective: the
    // radiant case peaks lower and later, because the room's mass takes it
    // first.
    let s = polar_night();
    let day = DesignDay { day_of_year: 200.0, t_max: 24.0, range: 0.0, w_out: 0.009, pressure: 101_325.0 };
    let run = |radiant: f64| {
        let mut z = room();
        z.windows.clear();
        z.outdoor_air = 0.0;
        z.gains = Gains { people: 0.0, person_sensible: 0.0, person_latent: 0.0, lighting: 0.0, equipment: 1000.0, radiant: (0.0, 0.0, radiant), schedule: (8.0, 18.0) };
        z.cooling_load(&s, &day, 24.0, 0.009, 4)
    };
    let (conv, rad) = (run(0.0), run(1.0));
    println!("1 kW 8–18 h: convective peak {:.0} W at {:.2} h, radiant peak {:.0} W at {:.2} h", conv.peak_sensible, conv.peak_hour, rad.peak_sensible, rad.peak_hour);
    assert!((conv.peak_sensible - 1000.0).abs() < 1.0);
    assert!(rad.peak_sensible < 0.9 * conv.peak_sensible && rad.peak_hour >= conv.peak_hour);
}
