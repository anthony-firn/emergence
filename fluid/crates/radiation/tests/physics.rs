use radiation::{
    math::{Bounds, Ray, Rng, Vec3},
    medium::{Medium, Volume, Voxel},
    optics::*,
    thermal::*,
    transport::*,
    A_RAD, C, SIGMA_SB,
};
use std::f64::consts::PI;

fn slab(n: usize, a: f64, s: f64, g: f64) -> Volume {
    Volume::new(
        [n, 1, 1],
        Bounds {
            lo: Vec3([0.0, -10.0, -10.0]),
            hi: Vec3([1.0, 10.0, 10.0]),
        },
        vec![Medium::grey(a, s, g)],
        vec![
            Voxel {
                medium: 0,
                scale: 1.0,
                temperature: 300.0
            };
            n
        ],
    )
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} vs {b:e}, tolerance {tol:e}");
}
#[test]
fn planck_integrates_to_stefan_boltzmann() {
    for t in [100.0, 300.0, 1800.0, 6500.0] {
        let n = 40000;
        let lo = (1e-9_f64).ln();
        let hi = (0.1_f64).ln();
        let h = (hi - lo) / n as f64;
        let flux = PI
            * h
            * (0..=n)
                .map(|i| {
                    let l = (lo + h * i as f64).exp();
                    let w = if i == 0 || i == n { 0.5 } else { 1.0 };
                    w * planck(l, t) * l
                })
                .sum::<f64>();
        close(flux / (SIGMA_SB * t.powi(4)), 1.0, 2e-9);
    }
}
#[test]
fn fresnel_normal_incidence_and_total_internal_reflection() {
    let ray = Vec3([0.5, 0.0, -(0.75_f64).sqrt()]);
    let transmitted = refract(ray, Vec3([0.0, 0.0, 1.0]), 1.0, 1.5).unwrap();
    close(transmitted.0[0], 1.0 / 3.0, 1e-15);
    assert!(refract(
        Vec3([(0.75_f64).sqrt(), 0.0, -0.5]),
        Vec3([0.0, 0.0, 1.0]),
        1.5,
        1.0
    )
    .is_none());
    for n in [1.1, 1.333, 1.5, 2.4] {
        close(
            fresnel_dielectric(1.0, 1.0, n),
            ((n - 1.0) / (n + 1.0)).powi(2),
            1e-15,
        );
    }
    close(fresnel_dielectric(0.5, 1.5, 1.0), 1.0, 0.0);
    let (n, k) = (2.4, 3.1);
    close(
        fresnel_conductor(1.0, n, k),
        ((n - 1.0).powi(2) + k * k) / ((n + 1.0).powi(2) + k * k),
        1e-15,
    );
    let water = OpticalMaterial::water();
    assert!(water.n.at(550e-9) > 1.33 && water.n.at(550e-9) < 1.34);
    assert!(water.absorption(700e-9) > water.absorption(450e-9));
}
#[test]
fn hg_normalization_and_sampled_asymmetry() {
    let axis = Vec3([0.0, 0.0, 1.0]);
    for g in [-0.6, 0.0, 0.8] {
        let n = 100000;
        let integral = (0..n)
            .map(|i| hg(-1.0 + 2.0 * (i as f64 + 0.5) / n as f64, g))
            .sum::<f64>()
            * 4.0
            * PI
            / n as f64;
        close(integral, 1.0, 1e-7);
        let mut rng = Rng(12);
        let mean = (0..n)
            .map(|_| sample_hg(axis, g, &mut rng).dot(axis))
            .sum::<f64>()
            / n as f64;
        close(mean, g, 0.006);
    }
}
#[test]
fn exact_voxel_transmittance_faces_corners_and_reverse_rays() {
    let mut v = slab(10, 2.0, 0.0, 0.0);
    for (i, c) in v.cells.iter_mut().enumerate() {
        c.scale = (i + 1) as f64;
    }
    for r in [
        Ray::new(Vec3([-1.0, 0.0, 0.0]), Vec3([1.0, 0.0, 0.0])),
        Ray::new(Vec3([2.0, 0.0, 0.0]), Vec3([-1.0, 0.0, 0.0])),
    ] {
        let seg: Vec<_> = v.segments(r, 100.0).collect();
        assert_eq!(seg.len(), 10);
        close(seg.iter().map(|s| s.end - s.start).sum(), 1.0, 1e-14);
        close(v.transmittance(r, 100.0, 550e-9), (-11.0_f64).exp(), 1e-17);
    }
    let cube = Volume::new(
        [8; 3],
        Bounds {
            lo: Vec3([0.0; 3]),
            hi: Vec3([1.0; 3]),
        },
        vec![Medium::grey(2.0, 0.0, 0.0)],
        vec![
            Voxel {
                medium: 0,
                scale: 1.0,
                temperature: 300.0
            };
            512
        ],
    );
    let r = Ray::new(Vec3([-1.0; 3]), Vec3([1.0; 3]));
    close(
        cube.transmittance(r, 100.0, 550e-9),
        (-2.0 * 3.0_f64.sqrt()).exp(),
        1e-15,
    );
    let reverse = Ray::new(Vec3([2.0; 3]), Vec3([-1.0; 3]));
    close(
        cube.transmittance(reverse, 100.0, 550e-9),
        cube.transmittance(r, 100.0, 550e-9),
        1e-15,
    );
    assert_eq!(cube.segments(r, 100.0).count(), 8);
}
#[test]
fn packet_slab_matches_beer_lambert_and_conserves_energy_momentum() {
    let v = slab(12, 2.0, 0.0, 0.0);
    let n = 40000;
    let packets = (0..n)
        .map(|i| Packet {
            ray: Ray::new(Vec3([0.0; 3]), Vec3([1.0, 0.0, 0.0])),
            wavelength: 550e-9,
            energy: 1.0 / n as f64,
            rng_state: i as u64,
        })
        .collect();
    let flights = fly_batch(&v, packets, 2.0 / C, 100);
    let escaped: f64 = flights.iter().map(|f| f.escaped_energy).sum();
    close(escaped, (-2.0_f64).exp(), 0.004);
    let absorbed: f64 = flights
        .iter()
        .flat_map(|f| &f.deposits)
        .map(|d| d.energy)
        .sum();
    close(escaped + absorbed, 1.0, 1e-12);
    let px: f64 = flights
        .iter()
        .map(|f| f.escaped_momentum.0[0] + f.deposits.iter().map(|d| d.momentum.0[0]).sum::<f64>())
        .sum();
    close(px * C, 1.0, 1e-12);
}
#[test]
fn finite_time_and_event_budget_retain_live_packets() {
    let v = slab(4, 0.0, 0.0, 0.0);
    let p = Packet {
        ray: Ray::new(Vec3([0.0; 3]), Vec3([1.0, 0.0, 0.0])),
        wavelength: 500e-9,
        energy: 1.0,
        rng_state: 42,
    };
    let f = fly(&v, p.clone(), 0.25 / C, 1);
    close(f.live.unwrap().ray.origin.0[0], 0.25, 1e-15);
    let v = slab(8, 0.3, 100.0, 0.8);
    let mut f = fly(&v, p, 0.25 / C, 1);
    assert!(f.pending_time > 0.0);
    let mut energy = 0.0;
    let mut momentum = Vec3::default();
    loop {
        energy += f.escaped_energy + f.deposits.iter().map(|d| d.energy).sum::<f64>();
        momentum = momentum
            + f.escaped_momentum
            + f.deposits
                .iter()
                .map(|d| d.momentum)
                .fold(Vec3::default(), |a, b| a + b);
        if f.pending_time == 0.0 {
            if let Some(p) = f.live {
                energy += p.energy;
                momentum = momentum + p.ray.dir * (p.energy / C);
            }
            break;
        }
        f = fly(&v, f.live.unwrap(), f.pending_time, 1);
    }
    close(energy, 1.0, 1e-14);
    close((momentum * C - Vec3([1.0, 0.0, 0.0])).length(), 0.0, 1e-14);
}
#[test]
fn thermal_emission_withdraws_equal_energy_and_recoil() {
    let v = slab(3, 2.0, 0.0, 0.0);
    let available = vec![1e6; 3];
    let emitted = emit_band(&v, 1e-3, 10e-6, 1e-7, 4, &available, 19).unwrap();
    let e: f64 = emitted.packets.iter().map(|p| p.energy).sum();
    let removed: f64 = emitted.withdrawals.iter().map(|d| d.energy).sum();
    close(e + removed, 0.0, e * 1e-14);
    let p = emitted
        .packets
        .iter()
        .map(|p| p.ray.dir * (p.energy / C))
        .chain(emitted.withdrawals.iter().map(|d| d.momentum))
        .fold(Vec3::default(), |a, b| a + b);
    close(p.length(), 0.0, e / C * 1e-14);
    assert!(emit_band(&v, 1.0, 10e-6, 1e-7, 4, &[0.0; 3], 1).is_err());
}
#[test]
fn stiff_grey_exchange_conserves_and_reaches_lte() {
    for cv in [1.0, 1e3, 4e6] {
        for multiplier in [1e-12, 1.0, 1e12] {
            let t = 1000.0;
            let e = A_RAD * 300.0_f64.powi(4);
            let r = exchange(t, cv, e, 1.0, multiplier / C);
            close(
                cv * r.temperature + r.radiation,
                cv * t + e,
                (cv * t + e) * 2e-15,
            );
            assert!(r.radiation >= 0.0 && r.temperature >= 300.0 && r.temperature <= 1000.0);
            if multiplier == 1e12 {
                close(r.radiation / (A_RAD * r.temperature.powi(4)), 1.0, 1e-9);
            }
        }
    }
}
#[test]
fn thick_diffusion_converges_to_closed_cosine_mode() {
    let k = 1e6;
    let d = C / (3.0 * k);
    let end = 0.03 / d;
    let mut errors = vec![];
    for n in [16, 32, 64] {
        let dx = 1.0 / n as f64;
        let mut e: Vec<f64> = (0..n)
            .map(|i| 1.0 + 0.5 * (PI * (i as f64 + 0.5) * dx).cos())
            .collect();
        let count = (end / (0.05 * dx * dx / d)).ceil() as usize;
        let dt = end / count as f64;
        for _ in 0..count {
            diffuse_closed(&mut e, &vec![k; n], dx, dt).unwrap();
        }
        close(e.iter().sum::<f64>() / n as f64, 1.0, 2e-12);
        assert!(e.iter().all(|&v| v > 0.0));
        let err = e
            .iter()
            .enumerate()
            .map(|(i, &e)| {
                (e - (1.0 + 0.5 * (PI * (i as f64 + 0.5) * dx).cos() * (-d * PI * PI * end).exp()))
                    .powi(2)
            })
            .sum::<f64>()
            .sqrt()
            / (n as f64).sqrt();
        println!("diffusion n={n}: rms={err:.6e}");
        errors.push(err);
    }
    assert!(errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8);
    assert!(diffuse_closed(&mut [1.0; 2], &[0.0; 2], 1.0, 1.0).is_err());
}
#[test]
fn anisotropy_blocks_a_false_diffusion_fold() {
    let mut i = Indicators {
        transport_depth: 100.0,
        mean_free_path_over_gradient_length: 0.01,
        reduced_flux: 0.01,
        absorption_relaxations: 1000.0,
        lte_relative_error: 0.001,
        directional_structure: false,
        wavelength_over_feature: 1e-5,
    };
    assert_eq!(i.regime(), Regime::LocalEquilibrium);
    i.directional_structure = true;
    assert_eq!(i.regime(), Regime::AngularPackets);
    i.directional_structure = false;
    i.wavelength_over_feature = 0.1;
    assert_eq!(i.regime(), Regime::WaveOpticsRequired);
}
