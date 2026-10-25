//! Euler buckling inside the macroscopic state: the instability boundary of
//! the linear homogenised model.
//!
//! The homogenised stiffness is linear and says nothing about when a strut
//! stops carrying its share of the load. Under compression the axial force
//! `N` in the most loaded strut grows linearly with the macroscopic strain,
//! and the strut buckles when `|N|` reaches its Euler load. The critical
//! macroscopic strain is therefore just the ratio of the two:
//!
//! ```text
//!   eps_cr = P_euler / |N(eps = 1)|.
//! ```
//!
//! `P_euler = pi^2 E I / (K L)^2` with `K` the effective-length factor
//! (`K = 1` pinned-pinned, `0.5` built-in). Real lattice joints lie between;
//! the pinned value is the conservative lower bound reported here.

use crate::cell::{Cell2d, Cell3d};
use std::f64::consts::PI;

/// Buckling result for one macroscopic strain path.
#[derive(Clone, Copy, Debug)]
pub struct Buckling {
    /// The member that buckles first.
    pub member: usize,
    /// Axial force per unit macroscopic strain in that member (compression
    /// positive).
    pub force_per_strain: f64,
    /// Euler load of that member.
    pub p_euler: f64,
    /// Macroscopic strain at which the strut buckles.
    pub critical_strain: f64,
}

/// Uniaxial compression along `direction` (0 = x, 1 = y) of a 2D cell.
/// `end_factor` is the effective-length factor `K`.
pub fn uniaxial_2d(cell: &Cell2d, direction: usize, end_factor: f64) -> Buckling {
    assert!(direction < 2);
    let mut strain = [0.0; 3];
    strain[direction] = -1.0;
    // The fluctuation (and hence the forces) are linear in strain, so a unit
    // strain gives the per-unit-strain forces directly.
    let forces = cell.axial_forces(strain);
    let mut best: Option<Buckling> = None;
    for (i, m) in cell.members.iter().enumerate() {
        let dx = cell.physical(m.b, m.b_img)[0] - cell.physical(m.a, m.a_img)[0];
        let dy = cell.physical(m.b, m.b_img)[1] - cell.physical(m.a, m.a_img)[1];
        let l = (dx * dx + dy * dy).sqrt();
        let p_euler = PI * PI * m.e * m.inertia / (end_factor * l).powi(2);
        let compression = -forces[i];
        if compression <= 0.0 {
            continue;
        }
        let crit = p_euler / compression;
        if best.map(|b| crit < b.critical_strain).unwrap_or(true) {
            best = Some(Buckling {
                member: i,
                force_per_strain: compression,
                p_euler,
                critical_strain: crit,
            });
        }
    }
    best.expect("no compressed strut under uniaxial compression")
}

/// Uniaxial compression along `direction` (0 = x, 1 = y, 2 = z) of a 3D cell.
pub fn uniaxial_3d(cell: &Cell3d, direction: usize, end_factor: f64) -> Buckling {
    assert!(direction < 3);
    let mut strain = [0.0; 6];
    strain[direction] = -1.0;
    let forces = cell.axial_forces(strain);
    let mut best: Option<Buckling> = None;
    for (i, m) in cell.members.iter().enumerate() {
        let xa = cell.physical(m.a, m.a_img);
        let xb = cell.physical(m.b, m.b_img);
        let d = [xb[0] - xa[0], xb[1] - xa[1], xb[2] - xa[2]];
        let l = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        // Weakest bending axis.
        let inertia = m.iy.min(m.iz);
        let p_euler = PI * PI * m.e * inertia / (end_factor * l).powi(2);
        let compression = -forces[i];
        if compression <= 0.0 {
            continue;
        }
        let crit = p_euler / compression;
        if best.map(|b| crit < b.critical_strain).unwrap_or(true) {
            best = Some(Buckling {
                member: i,
                force_per_strain: compression,
                p_euler,
                critical_strain: crit,
            });
        }
    }
    best.expect("no compressed strut under uniaxial compression")
}

/// Relative-density scaling of the buckling strain. For bending-dominated
/// lattices the axial force at unit strain scales as `E t`, so
/// `eps_cr ~ E t^3 / (E t l^2) = (t/l)^2`. Returning the log-log slope of
/// `eps_cr` against `t/l` lets the tests confirm this exponent.
pub fn buckling_exponent(mut sample: impl FnMut(f64) -> f64, rho1: f64, rho2: f64) -> f64 {
    let (e1, e2) = (sample(rho1), sample(rho2));
    (e2.ln() - e1.ln()) / (rho2.ln() - rho1.ln())
}
