use compressible::{eos::Nasg, Euler, Mesh, Solver};
use standard_model::{FourMomentum, Vec3, C, GEV_J};
#[test]
fn kinetic_sources_apply_gevs_in_si_and_reject_atomically() {
    let mut s = Solver::new(
        Euler {
            eos: Nasg::ideal(1.4),
        },
        Mesh::new([4, 1, 1], [0.0; 3], [1.0; 3]),
        |_| [1.0, 0.0, 0.0, 0.0, 1.0],
    );
    let q = s.mesh.idx(s.mesh.g(0), 0, 0);
    let old = s.u[q];
    // 1 J in GeV, including the momentum of absorption of a physical photon.
    let transfer = FourMomentum {
        energy: 1.0 / GEV_J,
        momentum: Vec3([1.0 / GEV_J, 0.0, 0.0]),
    };
    let amount = compressible::standard_model::deposit(&mut s, &[(0, transfer)]).unwrap();
    assert!((amount - 1.0).abs() < 1e-15);
    assert!((s.u[q][4] - old[4] - 4.0).abs() < 1e-14);
    assert!((s.u[q][1] - 4.0 / C).abs() < 1e-22);
    let before = s.u.clone();
    let cooling = FourMomentum {
        energy: -1e30,
        momentum: Vec3::default(),
    };
    assert!(compressible::standard_model::deposit(&mut s, &[(0, transfer), (1, cooling)]).is_err());
    assert_eq!(s.u, before);
    assert!(compressible::standard_model::deposit(&mut s, &[(4, transfer)]).is_err());
    assert_eq!(s.u, before);
    let nan = FourMomentum {
        energy: f64::NAN,
        momentum: Vec3::default(),
    };
    assert!(compressible::standard_model::deposit(&mut s, &[(0, transfer), (1, nan)]).is_err());
    assert_eq!(s.u, before);
}
