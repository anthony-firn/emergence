use gravity::{norm, Vec3};
use solar::{
    corona::{self, Pfss},
    surface::{Parameters, Surface},
};
use std::f64::consts::TAU;
fn map(n: usize, k: usize, f: impl Fn(f64, f64) -> f64) -> Surface {
    let mut b = vec![];
    for i in 0..n {
        let mu = -1.0 + (i as f64 + 0.5) * 2.0 / n as f64;
        for j in 0..k {
            b.push(f(mu, (j as f64 + 0.5) * TAU / k as f64));
        }
    }
    Surface::new(
        n,
        k,
        Parameters {
            radius: 1.0,
            ..Parameters::default()
        },
        b,
    )
    .unwrap()
}
#[test]
fn dipole_pfss_matches_analytic_field_and_radial_source_boundary() {
    let s = map(256, 1, |mu, _| mu);
    let field = Pfss::from_surface(&s, 1, 2.5).unwrap();
    let polar = field.magnetic_field(Vec3([0.0, 0.0, 1.0])).unwrap();
    assert!((polar.0[2] - 1.0).abs() < 2e-5);
    let ratio = 2.5f64.powi(-3);
    let equator = field.magnetic_field(Vec3([1.0, 0.0, 0.0])).unwrap();
    assert!((equator.0[2] - (ratio - 1.0) / (ratio + 2.0)).abs() < 2e-5);
    let point = Vec3([1.5, 0.0, 2.0]);
    let out = field.magnetic_field(point).unwrap();
    assert!(norm(out.cross(point)) < 1e-14);
    assert!(field.open_flux(64, 128).unwrap() > 0.0);
}

#[test]
fn harmonic_fit_integrates_cell_averages_instead_of_sampling_centers() {
    // Br=sign(mu) is represented exactly by any even equal-area latitude
    // grid. Analytic Legendre coefficients are B1=3/2 and B3=-7/8.
    // At the pole the l<=3 projection is therefore 5/8, even on 8 cells.
    let s = map(8, 1, |mu, _| mu.signum());
    let field = Pfss::from_surface(&s, 3, 2.5).unwrap();
    let pole = field.magnetic_field(Vec3([0.0, 0.0, 1.0])).unwrap();
    assert!((pole.0[2] - 0.625).abs() < 1e-12, "pole={pole:?}");
}
#[test]
fn nonaxisymmetric_coronal_field_is_curl_free_divergence_free_and_pole_regular() {
    let s = map(256, 64, |mu, phi| (1.0 - mu * mu).sqrt() * phi.cos());
    let field = Pfss::from_surface(&s, 1, 2.5).unwrap();
    let pole = field.magnetic_field(Vec3([0.0, 0.0, 1.5])).unwrap();
    let near = field.magnetic_field(Vec3([1e-10, 1e-10, 1.5])).unwrap();
    assert!(norm(pole - near) < 1e-9);
    let p = Vec3([1.1, 0.4, 0.8]);
    let h = 1e-5;
    let mut gradient = [[0.0; 3]; 3];
    for (axis, row) in gradient.iter_mut().enumerate() {
        let mut lo = p;
        let mut hi = p;
        lo.0[axis] -= h;
        hi.0[axis] += h;
        let d = (field.magnetic_field(hi).unwrap() - field.magnetic_field(lo).unwrap()) / (2.0 * h);
        row.copy_from_slice(&d.0);
    }
    let divergence = gradient[0][0] + gradient[1][1] + gradient[2][2];
    assert!(divergence.abs() < 1e-9);
    for (i, row) in gradient.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            assert!((value - gradient[j][i]).abs() < 1e-9);
        }
    }
}
#[test]
fn coronal_energy_is_a_diagnostic_and_incompatible_boundaries_are_rejected() {
    let s = map(32, 32, |_, _| 1.0);
    assert!(Pfss::from_surface(&s, 4, 2.5).is_err());
    let b = [Vec3([1.0, 0.0, 0.0])];
    let p = [Vec3([2.0, 0.0, 0.0])];
    let e = corona::energy(&b, &p, &[1.0]).unwrap();
    assert!(e.excess < 0.0);
    assert!(corona::energy(&b, &p, &[-1.0]).is_err());
}
