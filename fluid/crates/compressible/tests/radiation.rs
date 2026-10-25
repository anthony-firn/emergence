use compressible::{eos::Nasg, Euler, Mesh, Solver, System};
use radiation::{
    math::Vec3,
    medium::{Medium, Voxel},
    transport::Deposit,
    A_RAD, C,
};
fn gas() -> Solver<5, Euler> {
    let eos = Nasg {
        gamma: 5.0 / 3.0,
        p_inf: 0.0,
        b: 0.0,
        q: 0.0,
        cv: 1000.0,
        q_prime: 0.0,
    };
    Solver::new(
        Euler { eos },
        Mesh::new([4, 1, 1], [0.0; 3], [1.0; 3]),
        |_| [1.0, 0.0, 0.0, 0.0, (eos.gamma - 1.0) * eos.cv * 1000.0],
    )
}
#[test]
fn grey_source_preserves_fluid_plus_radiation_and_correct_cv() {
    let mut s = gas();
    let mut e = vec![A_RAD * 300.0_f64.powi(4); 4];
    let old = s.prim(0, 0, 0);
    let cv = s.sys.phases(&old)[0].2;
    assert_eq!(cv, 1000.0);
    let total = |s: &Solver<5, Euler>, e: &[f64]| {
        (0..4)
            .map(|i| s.u[s.mesh.idx(i + s.mesh.g(0), 0, 0)][4] + e[i])
            .sum::<f64>()
    };
    let before = total(&s, &e);
    compressible::radiation::exchange_grey(&mut s, &mut e, &[1.0; 4], 1e12 / C).unwrap();
    assert!((total(&s, &e) / before - 1.0).abs() < 1e-15);
    let t = s.sys.temperature(&s.prim(0, 0, 0));
    assert!((e[0] / (A_RAD * t.powi(4)) - 1.0).abs() < 1e-9);
}
#[test]
fn snapshot_and_photon_sources_use_interior_indices_and_si_volumes() {
    let mut s = gas();
    let v =
        compressible::radiation::snapshot(&s, vec![Medium::grey(2.0, 0.0, 0.0)], |sys, w| Voxel {
            medium: 0,
            scale: w[0],
            temperature: sys.temperature(w),
        });
    assert_eq!(v.cells.len(), 4);
    assert_eq!(v.cells[0].temperature, 1000.0);
    assert_eq!(v.cell_volume(), 0.25);
    let idx = s.mesh.idx(s.mesh.g(0), 0, 0);
    let old = s.u[idx];
    let d = Deposit {
        cell: 0,
        energy: 1.0,
        momentum: Vec3([1.0 / C, 0.0, 0.0]),
    };
    assert_eq!(compressible::radiation::deposit(&mut s, &[d]).unwrap(), 1.0);
    assert_eq!(s.u[idx][4] - old[4], 4.0);
    assert_eq!(s.u[idx][1], 4.0 / C);
    let before = s.u.clone();
    let bad = [
        Deposit {
            cell: 0,
            energy: 1.0,
            momentum: Vec3::default(),
        },
        Deposit {
            cell: 1,
            energy: -1e20,
            momentum: Vec3::default(),
        },
    ];
    assert!(compressible::radiation::deposit(&mut s, &bad).is_err());
    assert_eq!(s.u, before);
}

#[test]
fn grey_exchange_rejects_loss_of_internal_thermal_equilibrium_atomically() {
    use compressible::multiphase::FiveEq;
    let liquid = Nasg {
        gamma: 2.0,
        p_inf: 1e6,
        b: 0.0,
        q: 0.0,
        cv: 1000.0,
        q_prime: 0.0,
    };
    let gas = Nasg {
        p_inf: 0.0,
        gamma: 1.4,
        ..liquid
    };
    let sys = FiveEq {
        eos: [liquid, gas],
        kapila: true,
    };
    let mut s = Solver::new(sys, Mesh::new([1; 3], [0.0; 3], [1.0; 3]), |_| {
        [0.55, 0.125, 0.0, 0.0, 0.0, 1e5, 0.5]
    });
    let before = s.u.clone();
    let mut e = [1e8];
    let before_e = e;
    assert!(compressible::radiation::exchange_grey(&mut s, &mut e, &[1.0], 1e12 / C).is_err());
    assert_eq!(s.u, before);
    assert_eq!(e, before_e);
}
