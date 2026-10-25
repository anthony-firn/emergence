//! Explicit Newtonian bulk/fluid adapter. Rebuild sources after fluid advection;
//! a source kick alone is not a well-balanced hydrostatic discretization.
use crate::{Solver, System};
use ::gravity::{Body, Field, Tree, Vec3};
/// Interior cell masses [kg] at cell centers, including composite nuclei and
/// ordinary matter through the bulk density; no free-quark surrogate for ions.
pub fn sources<const N: usize, S: System<N>>(s: &Solver<N, S>) -> Result<Vec<Body>, &'static str> {
    let m = &s.mesh;
    let volume = m.dx.iter().product::<f64>();
    if !volume.is_finite() || volume <= 0.0 {
        return Err("invalid gravitational cell volume");
    }
    let mut bodies = Vec::with_capacity(m.cells());
    for i in 0..m.n[0] {
        for j in 0..m.n[1] {
            for k in 0..m.n[2] {
                let rho = s.sys.density(&s.prim(i, j, k));
                if !rho.is_finite() || rho < 0.0 {
                    return Err("invalid gravitational mass density");
                }
                bodies.push(Body {
                    position: Vec3(m.x(i, j, k)),
                    mass: rho * volume,
                });
            }
        }
    }
    Ok(bodies)
}
pub fn self_field<const N: usize, S: System<N>>(
    s: &Solver<N, S>,
    softening: f64,
    theta: f64,
) -> Result<Tree, &'static str> {
    Tree::new(sources(s)?, softening, theta)
}
/// Frozen field impulse with exact kinetic work, preserving internal energy.
/// Prepares all cells before committing. Returns gravitational work [J]. The
/// caller supplies temporal splitting/field rebuilds and a hydrostatic balance
/// closure. Self-gravity total-energy conservation additionally needs consistent
/// gravitational potential work during mass advection.
pub fn kick<const N: usize, S: System<N>, F: Field>(
    s: &mut Solver<N, S>,
    field: &F,
    dt: f64,
) -> Result<f64, &'static str> {
    if !dt.is_finite() || dt < 0.0 {
        return Err("invalid fluid gravity interval");
    }
    let ie = s.sys.energy_index().ok_or("fluid has no energy equation")?;
    let m = &s.mesh;
    let volume = m.dx.iter().product::<f64>();
    let mut updates = Vec::with_capacity(m.cells());
    let mut work = 0.0;
    for i in 0..m.n[0] {
        for j in 0..m.n[1] {
            for k in 0..m.n[2] {
                let index = m.idx(i + m.g(0), j + m.g(1), k + m.g(2));
                let mut u = s.u[index];
                let w = s.sys.to_prim(&u);
                let rho = s.sys.density(&w);
                if !rho.is_finite() || rho <= 0.0 {
                    return Err("invalid gravitational cell density");
                }
                let dv = field.sample(Vec3(m.x(i, j, k)))?.validate()?.acceleration * dt;
                let mut energy = 0.0;
                for d in 0..3 {
                    let momentum = s.sys.momentum_index(d);
                    energy += u[momentum] * dv.0[d] + 0.5 * rho * dv.0[d] * dv.0[d];
                    u[momentum] += rho * dv.0[d];
                }
                u[ie] += energy;
                if u.iter().any(|x| !x.is_finite()) || !s.sys.admissible(&s.sys.to_prim(&u)) {
                    return Err("gravity source produces inadmissible fluid state");
                }
                updates.push((index, u));
                work += energy * volume;
            }
        }
    }
    if !work.is_finite() {
        return Err("unrepresentable gravitational work");
    }
    for (index, u) in updates {
        s.u[index] = u;
    }
    // Density/internal energy are unchanged, so carried thermodynamic pressure
    // remains valid in the low-Mach path.
    Ok(work)
}
