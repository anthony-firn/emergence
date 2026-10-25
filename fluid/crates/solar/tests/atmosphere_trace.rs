use gravity::{norm, Vec3};
use solar::{
    corona::Pfss,
    photosphere,
    render::{Options as RenderOptions, View},
    surface::{Parameters, Surface},
    trace::{self, End, Options},
};
use std::f64::consts::PI;

fn dipole() -> Pfss {
    let n = 256;
    let values = (0..n)
        .map(|i| -1.0 + (i as f64 + 0.5) * 2.0 / n as f64)
        .collect();
    let surface = Surface::new(
        n,
        1,
        Parameters {
            radius: 1.0,
            ..Parameters::default()
        },
        values,
    )
    .unwrap();
    Pfss::from_surface(&surface, 1, 2.5).unwrap()
}
#[test]
fn grey_atmosphere_has_normalized_flux_and_physical_limb_darkening() {
    let t = 5772.0;
    assert!((photosphere::temperature(t, 2.0 / 3.0).unwrap() / t - 1.0).abs() < 1e-15);
    let n = 128;
    let flux = (0..=n)
        .map(|i| {
            let mu = i as f64 / n as f64;
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            weight * photosphere::bolometric(mu, t).unwrap() * mu
        })
        .sum::<f64>()
        * 2.0
        * PI
        / (3.0 * n as f64);
    assert!((flux / (photosphere::SIGMA * t.powi(4)) - 1.0).abs() < 1e-14);
    let limb = photosphere::bolometric(0.0, t).unwrap();
    let center = photosphere::bolometric(1.0, t).unwrap();
    assert!((limb / center - 0.4).abs() < 1e-15);
    assert!(photosphere::bolometric(f64::NAN, t).is_err());
    assert!(photosphere::temperature(t, -1.0).is_err());
}
#[test]
fn spectral_transfer_matches_independent_integral_and_bolometric_energy() {
    let t = 5772.0;
    for mu in [0.0, 0.1, 1.0] {
        for wavelength in [380e-9, 550e-9, 780e-9] {
            let n = 4096;
            let h = 40.0 / n as f64;
            let integral = (0..=n)
                .map(|i| {
                    let u = i as f64 * h;
                    let weight = if i == 0 || i == n {
                        1.0
                    } else if i % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    let temperature = photosphere::temperature(t, mu * u).unwrap();
                    weight * (-u).exp() * radiation::optics::planck(wavelength, temperature)
                })
                .sum::<f64>()
                * h
                / 3.0;
            let spectral = photosphere::spectral(mu, wavelength, t).unwrap();
            assert!(
                (spectral / integral - 1.0).abs() < 5e-6,
                "mu={mu}, lambda={wavelength}, ratio={}",
                spectral / integral
            );
        }
        // Independent spectral integration on a log-wavelength grid.
        let n = 1200;
        let lo = 1e-9_f64.ln();
        let h = (1e-2_f64.ln() - lo) / n as f64;
        let integral = (0..=n)
            .map(|i| {
                let wavelength = (lo + i as f64 * h).exp();
                let weight = if i == 0 || i == n {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                weight * wavelength * photosphere::spectral(mu, wavelength, t).unwrap()
            })
            .sum::<f64>()
            * h
            / 3.0;
        assert!((integral / photosphere::bolometric(mu, t).unwrap() - 1.0).abs() < 1e-7);
    }
}
#[test]
fn traced_dipole_preserves_flux_function_and_locates_open_closed_boundaries() {
    let field = dipole();
    for (theta, is_open) in [
        (20.0_f64.to_radians(), true),
        (60.0_f64.to_radians(), false),
    ] {
        let seed = Vec3([theta.sin(), 0.0, theta.cos()]) * 1.0001;
        let line = trace::trace(
            &field,
            seed,
            Options {
                tolerance: 1e-7,
                ..Options::default()
            },
        )
        .unwrap();
        assert_eq!(line.open(), is_open);
        assert_eq!(line.closed(), !is_open);
        let flux = |p: Vec3| {
            let r = norm(p);
            let mu = p.0[2] / r;
            (r * r / 2.5_f64.powi(3) + 2.0 / r) * (1.0 - mu * mu)
        };
        let initial = flux(seed);
        for p in &line.points {
            assert!(norm(*p) >= 1.0 - 1e-14 && norm(*p) <= 2.5 + 1e-14);
            assert!((flux(*p) / initial - 1.0).abs() < 3e-6);
        }
        for (end, point) in line
            .ends
            .iter()
            .zip([line.points.first().unwrap(), line.points.last().unwrap()])
        {
            let expected = if *end == End::Photosphere { 1.0 } else { 2.5 };
            assert!((norm(*point) - expected).abs() < 1e-14);
        }
    }
}
#[test]
fn unresolved_traces_are_not_mislabeled_as_connectivity() {
    let s = Surface::new(
        32,
        1,
        Parameters {
            radius: 1.0,
            ..Parameters::default()
        },
        vec![0.0; 32],
    )
    .unwrap();
    let zero = Pfss::from_surface(&s, 1, 2.5).unwrap();
    let line = trace::trace(&zero, Vec3([1.1, 0.0, 0.0]), Options::default()).unwrap();
    assert_eq!(line.ends, [End::Null; 2]);
    assert!(!line.open() && !line.closed());
    let incomplete = trace::trace(
        &dipole(),
        Vec3([1.1, 0.0, 0.0]),
        Options {
            max_length: 0.01,
            ..Options::default()
        },
    )
    .unwrap();
    assert_eq!(incomplete.ends, [End::LengthLimit; 2]);
    assert!(!incomplete.open() && !incomplete.closed());
    assert!(trace::trace(&dipole(), Vec3([0.1, 0.0, 0.0]), Options::default()).is_err());
}
#[test]
fn rendering_occults_rear_geometry_and_retains_continuum_radiance() {
    let view = View {
        longitude: 0.0,
        latitude: 0.0,
        extent: 1.2,
    };
    assert!(view.visible(Vec3([2.0, 0.0, 0.0])).unwrap());
    assert!(!view.visible(Vec3([-2.0, 0.0, 0.0])).unwrap());
    assert!(view.visible(Vec3([-2.0, 2.0, 0.0])).unwrap());
    let image = solar::render::continuum(
        5772.0,
        RenderOptions {
            width: 32,
            height: 32,
            samples_per_axis: 2,
            view,
        },
    )
    .unwrap();
    assert!(image.rgb.iter().flatten().all(|v| v.is_finite()));
    assert_eq!(image.rgb[0], [0.0; 3]);
    let center = image.rgb[16 * 32 + 16];
    assert!(center.iter().all(|v| *v > 0.0));
    assert!(center.iter().sum::<f64>() > image.rgb[16 * 32 + 27].iter().sum::<f64>());
}
