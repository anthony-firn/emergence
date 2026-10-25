//! Crystal geometry for the atomistic rung and the atomistic-to-continuum
//! exchange-stiffness fold.
//!
//! Conventional cubic cells are used throughout: an sc cell holds one atom,
//! bcc two and fcc four, all on a cube of side `a`.  Nearest neighbours are
//! found by scanning the 27 surrounding cells and comparing the bond length
//! with the exact nearest-neighbour distance, so the coordination, the bond
//! vectors and the exchange stiffness are derived rather than tabulated.
//!
//! # The fold: exchange stiffness
//!
//! With the classical Heisenberg Hamiltonian `H = -J sum_<ij> S_i . S_j`
//! (unit spins, magnitude `S`) the small-angle expansion of the bond energy
//! gives, per atom,
//! `e = J S^2/(4 V_atom) sum_delta (delta . grad m)^2` and hence, for a
//! cubic crystal, the continuum energy density `A (grad m)^2` with
//!
//! ```text
//! A = J S^2 sum_delta |delta|^2 / (12 V_atom) = J S^2 z d_nn^2 / (12 V_atom),
//! ```
//!
//! where the sum runs over the `z` nearest-neighbour bonds of one atom,
//! `V_atom` is the atomic volume and `d_nn` the nearest-neighbour distance.
//! For the three conventional cubic cells `sum_delta |delta|^2 = 6 a^2`, so
//! `A = J S^2 n_atoms / (2 a)` (`sc`: `J S^2/(2a)`, `bcc`: `J S^2/a`,
//! `fcc`: `2 J S^2/a`).  The often-quoted form `J S^2 z a^2/(2V)` carries a
//! different prefactor because it assumes the bond sum is `z a^2 (grad m)^2`
//! rather than the isotropic average `z d_nn^2 (grad m)^2/3`; the convention
//! above is the one used consistently here.

pub type V3 = [f64; 3];

pub fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LatticeKind {
    Sc,
    Bcc,
    Fcc,
}

impl LatticeKind {
    /// Fractional coordinates of the basis atoms in the conventional cell.
    pub fn basis(&self) -> &'static [V3] {
        match self {
            LatticeKind::Sc => &[[0.0, 0.0, 0.0]],
            LatticeKind::Bcc => &[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
            LatticeKind::Fcc => &[
                [0.0, 0.0, 0.0],
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [0.5, 0.5, 0.0],
            ],
        }
    }

    pub fn n_atoms(&self) -> usize {
        self.basis().len()
    }

    /// Conventional-cell volume.
    pub fn cell_volume(&self, a: f64) -> f64 {
        a * a * a
    }

    /// Volume per atom.
    pub fn atomic_volume(&self, a: f64) -> f64 {
        self.cell_volume(a) / self.n_atoms() as f64
    }

    /// Exact nearest-neighbour distance `d_nn(a)`.
    pub fn nn_distance(&self, a: f64) -> f64 {
        match self {
            LatticeKind::Sc => a,
            LatticeKind::Bcc => 0.5 * 3.0_f64.sqrt() * a,
            LatticeKind::Fcc => a / 2.0_f64.sqrt(),
        }
    }

    pub fn coordination(&self) -> usize {
        match self {
            LatticeKind::Sc => 6,
            LatticeKind::Bcc => 8,
            LatticeKind::Fcc => 12,
        }
    }

    /// Nearest-neighbour table: for each basis atom, the list of
    /// `(cell_offset, destination_basis, delta_vector)`.  The bond vector is
    /// `a * (cell_offset + basis[j] - basis[i])`.
    pub fn neighbors(&self, a: f64) -> Vec<Vec<Neighbor>> {
        let basis = self.basis();
        let nn = self.nn_distance(a);
        let tol = 1e-9 * nn.max(1.0);
        let mut out = vec![Vec::new(); basis.len()];
        for i in 0..basis.len() {
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let cell = [dx, dy, dz];
                        for (j, bj) in basis.iter().enumerate() {
                            let frac = [
                                cell[0] as f64 + bj[0] - basis[i][0],
                                cell[1] as f64 + bj[1] - basis[i][1],
                                cell[2] as f64 + bj[2] - basis[i][2],
                            ];
                            let delta = [a * frac[0], a * frac[1], a * frac[2]];
                            let r = norm(delta);
                            if (r - nn).abs() < tol {
                                out[i].push(Neighbor {
                                    cell,
                                    basis: j,
                                    delta,
                                });
                            }
                        }
                    }
                }
            }
        }
        out
    }

    /// `sum_delta |delta|^2` over the nearest-neighbour bonds of one atom.
    pub fn sum_delta_squared(&self, a: f64) -> f64 {
        self.neighbors(a)[0]
            .iter()
            .map(|n| dot(n.delta, n.delta))
            .sum()
    }

    /// The isotropic bond tensor `sum_delta delta_a delta_b = Lambda delta_ab`,
    /// with `Lambda = (1/3) sum |delta|^2` for a cubic crystal.
    pub fn bond_tensor_lambda(&self, a: f64) -> f64 {
        self.sum_delta_squared(a) / 3.0
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Neighbor {
    pub cell: [i32; 3],
    pub basis: usize,
    pub delta: V3,
}

/// Atomistic-to-continuum exchange stiffness `A` (J/m) for a cubic lattice
/// with exchange `J` (J), spin magnitude `S` and lattice constant `a` (m).
///
/// `A = J S^2 sum_delta |delta|^2 / (12 V_atom)`; see the module docs.
pub fn exchange_stiffness(kind: LatticeKind, j: f64, s: f64, a: f64) -> f64 {
    let lambda = kind.bond_tensor_lambda(a);
    j * s * s * 3.0 * lambda / (12.0 * kind.atomic_volume(a))
    // 3*Lambda = sum |delta|^2
}

/// The same relation written directly from `z` and `d_nn`.
pub fn exchange_stiffness_z(kind: LatticeKind, j: f64, s: f64, a: f64) -> f64 {
    let z = kind.coordination() as f64;
    let d = kind.nn_distance(a);
    j * s * s * z * d * d / (12.0 * kind.atomic_volume(a))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coordination_and_distances_match_known_lattices() {
        for (kind, z) in [
            (LatticeKind::Sc, 6usize),
            (LatticeKind::Bcc, 8),
            (LatticeKind::Fcc, 12),
        ] {
            let nb = kind.neighbors(1.0);
            assert_eq!(nb[0].len(), z, "{kind:?}");
            for n in &nb[0] {
                assert!((norm(n.delta) - kind.nn_distance(1.0)).abs() < 1e-12);
            }
        }
        assert!((LatticeKind::Bcc.nn_distance(1.0) - 3.0_f64.sqrt() / 2.0).abs() < 1e-12);
        assert!((LatticeKind::Fcc.nn_distance(1.0) - 1.0 / 2.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn sum_delta_squared_is_six_a_squared() {
        for kind in [LatticeKind::Sc, LatticeKind::Bcc, LatticeKind::Fcc] {
            let s = kind.sum_delta_squared(1.0);
            assert!((s - 6.0).abs() < 1e-12, "{kind:?} {s}");
        }
    }

    #[test]
    fn exchange_stiffness_matches_closed_forms() {
        // sc: J S^2/(2a), bcc: J S^2/a, fcc: 2 J S^2/a.
        let (j, s, a) = (1.0, 0.5, 1.0);
        assert!((exchange_stiffness(LatticeKind::Sc, j, s, a) - 0.125).abs() < 1e-12);
        assert!((exchange_stiffness(LatticeKind::Bcc, j, s, a) - 0.25).abs() < 1e-12);
        assert!((exchange_stiffness(LatticeKind::Fcc, j, s, a) - 0.5).abs() < 1e-12);
    }
}
