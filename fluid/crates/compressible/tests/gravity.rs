use compressible::{eos::Nasg, Euler, Mesh, Solver};
use gravity::{Field, Sample, Vec3};
struct Uniform(Vec3);
impl Field for Uniform {
    fn sample(&self, _: Vec3) -> Result<Sample, &'static str> {
        Ok(Sample {
            potential: 0.0,
            acceleration: self.0,
        })
    }
}
#[test]
fn gravity_kick_has_exact_kinetic_work_and_rejects_atomically() {
    let mut s = Solver::new(
        Euler {
            eos: Nasg::ideal(1.4),
        },
        Mesh::new([4, 1, 1], [0.0; 3], [1.0; 3]),
        |_| [2.0, 3.0, 0.0, 0.0, 1e5],
    );
    let work = compressible::gravity::kick(&mut s, &Uniform(Vec3([0.0, -10.0, 0.0])), 0.2).unwrap();
    assert!((work - 4.0).abs() < 1e-12);
    for i in 0..4 {
        let w = s.prim(i, 0, 0);
        assert!((w[2] + 2.0).abs() < 1e-12);
        assert!((w[4] - 1e5).abs() < 1e-8);
    }
    let before = s.u.clone();
    assert!(
        compressible::gravity::kick(&mut s, &Uniform(Vec3([f64::NAN, 0.0, 0.0])), 1.0).is_err()
    );
    assert_eq!(s.u, before);
    let sources = compressible::gravity::sources(&s).unwrap();
    assert!((sources.iter().map(|x| x.mass).sum::<f64>() - 2.0).abs() < 1e-14);
}
