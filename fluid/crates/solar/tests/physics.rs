use gravity::{norm, Field, Vec3, G};
use solar::{
    structure::{Profile, MODEL_S_MASS, MODEL_S_RADIUS},
    surface::{self, Bipole, Parameters, Surface},
};
use std::f64::consts::{PI, TAU};
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t, "{a:e} != {b:e}, tol {t:e}");
}
fn surface(n: usize, k: usize, p: Parameters, f: impl Fn(f64, f64) -> f64) -> Surface {
    let mut b = vec![];
    for i in 0..n {
        let mu = -1.0 + (i as f64 + 0.5) * 2.0 / n as f64;
        for j in 0..k {
            b.push(f(mu, (j as f64 + 0.5) * TAU / k as f64));
        }
    }
    Surface::new(n, k, p, b).unwrap()
}
#[test]
fn model_s_units_enclosed_mass_and_hydrostatic_balance_are_physical() {
    let p = Profile::model_s().unwrap();
    assert_eq!(p.layers.len(), 2482);
    close(p.layers[0].temperature, 1.5668470e7, 0.1);
    close(p.layers[0].density, 153889.36, 0.01);
    let g = p.gravity(Vec3::default()).unwrap();
    close(g.total_mass(), MODEL_S_MASS, 1e-4 * MODEL_S_MASS);
    assert_eq!(
        g.sample(Vec3::default()).unwrap().acceleration,
        Vec3::default()
    );
    let mut worst: f64 = 0.0;
    for f in [0.1, 0.25, 0.5, 0.75, 0.9] {
        let r = f * MODEL_S_RADIUS;
        let h = 0.001 * MODEL_S_RADIUS;
        let center = p.at(r).unwrap();
        let dp = (p.at(r + h).unwrap().pressure - p.at(r - h).unwrap().pressure) / (2.0 * h);
        let expected = -center.density * norm(g.sample(Vec3([r, 0.0, 0.0])).unwrap().acceleration);
        worst = worst.max((dp / expected - 1.0).abs());
    }
    assert!(worst < 0.005, "hydrostatic residual {worst}");
    let exterior = g.sample(Vec3([2.0 * MODEL_S_RADIUS, 0.0, 0.0])).unwrap();
    close(
        exterior.potential,
        -G * g.total_mass() / (2.0 * MODEL_S_RADIUS),
        1e-8,
    );
}

#[test]
fn meridional_transport_converges_to_exact_characteristic_cell_averages() {
    // u_mu=k*mu*(1-mu²). A starting B=mu has an analytic characteristic
    // primitive mu_initial²/2, allowing an independent cell-average reference.
    let mut errors = Vec::new();
    for n in [64, 128, 256] {
        let parameters = Parameters {
            radius: 1.0,
            meridional_peak: 0.2,
            diffusivity: 0.0,
            rotation: [0.0; 3],
            max_step: 1e-3,
        };
        let mut s = surface(n, 1, parameters, |mu, _| mu);
        s.advance(1.0).unwrap();
        let exponential = 0.8_f64.exp();
        let primitive = |mu: f64| 0.5 * mu * mu / (exponential + (1.0 - exponential) * mu * mu);
        let h = 2.0 / n as f64;
        let error: f64 = s
            .field()
            .iter()
            .enumerate()
            .map(|(i, b)| {
                let lo = -1.0 + i as f64 * h;
                (b - (primitive(lo + h) - primitive(lo)) / h).abs() * h
            })
            .sum();
        errors.push(error);
    }
    assert!(
        errors[0] / errors[1] > 3.0 && errors[1] / errors[2] > 3.0,
        "meridional spatial errors {errors:?}"
    );
}
#[test]
fn spherical_dipole_diffusion_matches_analytic_decay_and_second_order_time() {
    let mut errors = vec![];
    for dt in [0.1, 0.05, 0.025] {
        let params = Parameters {
            radius: 1.0,
            diffusivity: 1.0,
            meridional_peak: 0.0,
            rotation: [0.0; 3],
            max_step: dt,
        };
        let mut s = surface(64, 1, params, |mu, _| mu);
        let initial = s.axial_dipole();
        s.advance(1.0).unwrap();
        let ratio = s.axial_dipole() / initial;
        errors.push((ratio - (-2.0f64).exp()).abs());
    }
    assert!(
        errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8,
        "SDIRK errors {errors:?}"
    );
}

#[test]
fn stiff_diffusion_damps_unresolved_modes_and_preserves_uniform_flux() {
    let parameters = Parameters {
        radius: 1.0,
        diffusivity: 1.0,
        meridional_peak: 0.0,
        rotation: [0.0; 3],
        max_step: 100.0,
    };
    let mut s = surface(64, 1, parameters, |mu, _| 2.0 + mu);
    let initial = s.net_flux();
    s.advance(100.0).unwrap();
    assert!(s.field().iter().all(|b| (b - 2.0).abs() < 0.003));
    close(s.net_flux(), initial, 1e-10);
}
#[test]
fn spectral_rotation_shears_each_latitude_at_its_own_rate_without_energy_loss() {
    let params = Parameters {
        radius: 1.0,
        diffusivity: 0.0,
        meridional_peak: 0.0,
        rotation: [0.2, -0.1, 0.0],
        max_step: 0.3,
    };
    let mut s = surface(32, 64, params, |mu, phi| (1.0 - mu * mu).sqrt() * phi.cos());
    let before = s.field().iter().map(|b| b * b).sum::<f64>();
    s.advance(10.0).unwrap();
    for i in 0..32 {
        let mu = s.mu(i);
        let omega = 0.2 - 0.1 * mu * mu;
        for j in 0..64 {
            let expected = (1.0 - mu * mu).sqrt() * (s.longitude(j) - omega * 10.0).cos();
            close(s.field()[i * 64 + j], expected, 3e-14);
        }
    }
    close(s.field().iter().map(|b| b * b).sum::<f64>(), before, 1e-10);
}
#[test]
fn longitude_average_is_an_exact_emergent_fold_with_sources_and_transport() {
    let params = Parameters::default();
    let mut full = surface(64, 128, params, |mu, phi| {
        1e-4 * mu + 3e-5 * (1.0 - mu * mu) * phi.cos()
    });
    full.add_bipole(Bipole {
        latitude: 0.3,
        longitude: 1.0,
        separation: 0.15,
        tilt: 0.4,
        width: 0.1,
        flux: 1e13,
    })
    .unwrap();
    let net = full.net_flux();
    let mut folded = full.axisymmetric().unwrap();
    full.advance(90.0 * 86400.0).unwrap();
    folded.advance(90.0 * 86400.0).unwrap();
    for (a, b) in full.mean_field().iter().zip(folded.field()) {
        close(*a, *b, 2e-17);
    }
    close(full.axial_dipole(), folded.axial_dipole(), 1e-17);
    close(full.net_flux(), net, 5e-13 * full.unsigned_flux());
}
#[test]
fn bipole_flux_balance_boundaries_and_rejection_keep_state() {
    let mut s = surface(64, 512, Parameters::default(), |_, _| 0.0);
    let event = Bipole {
        latitude: 0.4,
        longitude: PI,
        separation: 0.3,
        tilt: 0.7,
        width: 0.04,
        flux: 1e13,
    };
    s.add_bipole(event).unwrap();
    close(s.net_flux(), 0.0, 1e-13 * s.unsigned_flux());
    close(s.unsigned_flux(), 2e13, 0.002 * 2e13);
    let before = s.field().to_vec();
    assert!(s.advance(f64::INFINITY).is_err());
    assert_eq!(s.field(), before);
    let bad = Bipole {
        flux: f64::NAN,
        ..event
    };
    assert!(s.add_bipole(bad).is_err());
    assert_eq!(s.field(), before);
}
#[test]
fn held_out_forecast_metrics_compare_persistence_without_fabricating_skill() {
    let observation = [1.0, -1.0];
    let prediction = [0.9, -0.9];
    let initial = [0.0, 0.0];
    let skill = surface::score(&observation, &prediction, &initial).unwrap();
    close(skill.rmse, 0.1, 1e-15);
    close(skill.mse_skill.unwrap(), 0.99, 1e-15);
    assert!(surface::score(&observation, &prediction, &observation)
        .unwrap()
        .mse_skill
        .is_none());
    assert!(surface::score(&[], &[], &[]).is_err());
}
