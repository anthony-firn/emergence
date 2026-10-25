//! Crystal structures as primitive vectors plus a basis, and the neighbour
//! shells and bonds that every higher rung is built from.
//!
//! A [`Lattice`] stores Cartesian primitive vectors `a1, a2, a3` (rows) and
//! fractional basis coordinates. An atom of basis `i` in cell `n` sits at
//! `Σ_c (n_c + basis_i[c]) a_c`. The bond list enumerates every atom inside
//! a spherical cutoff around a reference atom, including periodic images;
//! no minimum-image convention is used, so the sums are the true lattice
//! sums (cut off at `rmax`).
//!
//! Structures provided: simple cubic (sc), body-centred (bcc),
//! face-centred (fcc), hexagonal close-packed (hcp) with ideal or given
//! c/a, and diamond (two-atom fcc basis). These cover the exact shell
//! counts used as benchmarks (fcc 12, 6, 24, 12, ...; bcc 8, 6, 12;
//! diamond 4, 12).

pub type V3 = [f64; 3];

pub fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

/// Ideal hcp axial ratio c/a = sqrt(8/3), exact (not the rounded 1.633).
pub const IDEAL_C_OVER_A: f64 = 1.632_993_161_855_452;

/// A crystal structure, independent of the lattice constant. Hexagonal
/// close packing carries its own axial ratio so the two independent
/// parameters (a and c) can be varied.
#[derive(Clone, Debug, PartialEq)]
pub enum Crystal {
    Sc,
    Bcc,
    Fcc,
    Hcp {
        c_over_a: f64,
    },
    Diamond,
    /// Arbitrary primitive vectors and fractional basis, already in
    /// Cartesian units. Used for the 1D-chain benchmark.
    Custom {
        prim: [V3; 3],
        basis: Vec<V3>,
    },
}

impl Crystal {
    pub fn ideal_hcp() -> Self {
        Crystal::Hcp {
            c_over_a: IDEAL_C_OVER_A,
        }
    }

    /// Primitive vectors for lattice constant `a` (rows).
    pub fn prim(&self, a: f64) -> [V3; 3] {
        match self {
            Crystal::Sc => [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]],
            Crystal::Bcc => [
                [-0.5 * a, 0.5 * a, 0.5 * a],
                [0.5 * a, -0.5 * a, 0.5 * a],
                [0.5 * a, 0.5 * a, -0.5 * a],
            ],
            Crystal::Fcc | Crystal::Diamond => [
                [0.0, 0.5 * a, 0.5 * a],
                [0.5 * a, 0.0, 0.5 * a],
                [0.5 * a, 0.5 * a, 0.0],
            ],
            Crystal::Hcp { c_over_a } => {
                let c = c_over_a * a;
                [
                    [a, 0.0, 0.0],
                    [-0.5 * a, 0.5 * 3.0_f64.sqrt() * a, 0.0],
                    [0.0, 0.0, c],
                ]
            }
            Crystal::Custom { prim, .. } => *prim,
        }
    }

    /// Fractional basis coordinates.
    pub fn basis(&self) -> Vec<V3> {
        match self {
            Crystal::Sc | Crystal::Bcc | Crystal::Fcc => vec![[0.0, 0.0, 0.0]],
            Crystal::Hcp { .. } => vec![[0.0, 0.0, 0.0], [2.0 / 3.0, 1.0 / 3.0, 0.5]],
            Crystal::Diamond => vec![[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
            Crystal::Custom { basis, .. } => basis.clone(),
        }
    }
}

/// A lattice at a fixed size: primitive vectors and basis.
#[derive(Clone, Debug)]
pub struct Lattice {
    pub prim: [V3; 3],
    pub basis: Vec<V3>,
    pub label: String,
}

/// A neighbour of a reference atom: Cartesian bond vector `r`, its length,
/// the basis indices at each end, and the lattice translation of the
/// neighbour's cell.
#[derive(Clone, Copy, Debug)]
pub struct Bond {
    pub r: V3,
    pub radius: f64,
    pub i: usize,
    pub j: usize,
    pub cell: [i32; 3],
}

/// A group of bonds at the same radius.
#[derive(Clone, Debug, PartialEq)]
pub struct Shell {
    pub radius: f64,
    pub count: usize,
}

impl Lattice {
    pub fn new(crystal: &Crystal, a: f64) -> Self {
        let label = match crystal {
            Crystal::Sc => "sc".to_string(),
            Crystal::Bcc => "bcc".to_string(),
            Crystal::Fcc => "fcc".to_string(),
            Crystal::Hcp { c_over_a } => format!("hcp(c/a={c_over_a:.6})"),
            Crystal::Diamond => "diamond".to_string(),
            Crystal::Custom { .. } => "custom".to_string(),
        };
        Self {
            prim: crystal.prim(a),
            basis: crystal.basis(),
            label,
        }
    }

    pub fn sc(a: f64) -> Self {
        Self::new(&Crystal::Sc, a)
    }
    pub fn bcc(a: f64) -> Self {
        Self::new(&Crystal::Bcc, a)
    }
    pub fn fcc(a: f64) -> Self {
        Self::new(&Crystal::Fcc, a)
    }
    pub fn hcp(a: f64) -> Self {
        Self::new(&Crystal::ideal_hcp(), a)
    }
    pub fn hcp_ca(a: f64, c_over_a: f64) -> Self {
        Self::new(&Crystal::Hcp { c_over_a }, a)
    }
    pub fn diamond(a: f64) -> Self {
        Self::new(&Crystal::Diamond, a)
    }
    pub fn custom(prim: [V3; 3], basis: Vec<V3>) -> Self {
        Self::new(&Crystal::Custom { prim, basis }, 1.0)
    }

    /// Scale every primitive vector by `s`; fractional basis is unchanged,
    /// so the whole structure scales.
    pub fn scaled(&self, s: f64) -> Self {
        Self {
            prim: self.prim.map(|v| scale(v, s)),
            basis: self.basis.clone(),
            label: self.label.clone(),
        }
    }

    pub fn natoms(&self) -> usize {
        self.basis.len()
    }

    /// Primitive-cell (reference) volume.
    pub fn volume(&self) -> f64 {
        dot(self.prim[0], cross(self.prim[1], self.prim[2])).abs()
    }

    /// Cartesian position of `Σ n_c a_c + basis_offset`.
    pub fn cartesian(&self, n: [i32; 3], basis_offset: V3) -> V3 {
        let mut r = [0.0; 3];
        for c in 0..3 {
            let coeff = n[c] as f64 + basis_offset[c];
            for a in 0..3 {
                r[a] += coeff * self.prim[c][a];
            }
        }
        r
    }

    /// Nearest-neighbour distance of the structure: the shortest non-zero
    /// bond over basis pairs and cells within two periodic repeats. A local
    /// search is exact for every structure here and avoids sizing a box
    /// from the primitive vectors (which is pathological for the
    /// deliberately anisotropic 1D-chain custom lattice).
    ///
    /// The "non-zero" test is **relative to the primitive vectors**. It
    /// used to be an absolute `1e-9`, which is a sensible epsilon for the
    /// dimensionless lattice constants this crate's own benchmarks use but
    /// rejects *every* neighbour of a lattice whose constant is given in
    /// metres (about `1e-10`): the function then silently returned the
    /// shortest bond longer than a nanometre, which for copper is the
    /// fourth shell.
    pub fn nn_distance(&self) -> f64 {
        let n = self.basis.len();
        let scale = (0..3).map(|k| norm(self.prim[k])).fold(0.0, f64::max);
        let eps = 1e-9 * scale;
        let mut best = f64::INFINITY;
        for i in 0..n {
            for j in 0..n {
                for n0 in -2..=2 {
                    for n1 in -2..=2 {
                        for n2 in -2..=2 {
                            let r = self.cartesian([n0, n1, n2], self.basis_diff(i, j));
                            let radius = norm(r);
                            if radius > eps && radius < best {
                                best = radius;
                            }
                        }
                    }
                }
            }
        }
        best
    }

    /// All bonds within `rmax` of each basis atom, grouped by reference
    /// basis index.
    ///
    /// As in [`Lattice::nn_distance`], the test that skips the zero vector
    /// is relative to the primitive vectors. An absolute `1e-9` exceeds
    /// every bond of a lattice measured in metres, so the list came back
    /// empty and every lattice sum built on it was silently zero.
    pub fn bond_list(&self, rmax: f64) -> Vec<Vec<Bond>> {
        let n = self.basis.len();
        let scale = (0..3).map(|k| norm(self.prim[k])).fold(0.0, f64::max);
        let eps = 1e-9 * scale;
        // Largest inter-basis offset, used to size the image search box.
        let mut base_max = 0.0_f64;
        for i in 0..n {
            for j in 0..n {
                base_max = base_max.max(norm(self.cartesian([0, 0, 0], self.basis_diff(i, j))));
            }
        }
        let reach = rmax + base_max;
        let v0 = self.volume();
        // Perpendicular spacing between lattice planes: any lattice vector
        // with component n_c along the reciprocal direction has length at
        // least |n_c| * spacing_c, which bounds the search.
        let mut spacing = [0.0; 3];
        for c in 0..3 {
            let a1 = self.prim[(c + 1) % 3];
            let a2 = self.prim[(c + 2) % 3];
            spacing[c] = v0 / norm(cross(a1, a2));
        }
        let bounds: [i32; 3] = std::array::from_fn(|c| (reach / spacing[c]).ceil() as i32 + 1);
        let mut out = vec![Vec::new(); n];
        for i in 0..n {
            for n0 in -bounds[0]..=bounds[0] {
                for n1 in -bounds[1]..=bounds[1] {
                    for n2 in -bounds[2]..=bounds[2] {
                        let cell = [n0, n1, n2];
                        for j in 0..n {
                            let r = self.cartesian(cell, self.basis_diff(i, j));
                            let radius = norm(r);
                            if radius > eps && radius <= rmax {
                                out[i].push(Bond {
                                    r,
                                    radius,
                                    i,
                                    j,
                                    cell,
                                });
                            }
                        }
                    }
                }
            }
        }
        out
    }

    fn basis_diff(&self, i: usize, j: usize) -> V3 {
        sub(self.basis[j], self.basis[i])
    }

    /// Neighbour shells around basis atom `i`, grouped with relative
    /// tolerance `rel_tol` on the radius.
    pub fn shells(&self, i: usize, rmax: f64, rel_tol: f64) -> Vec<Shell> {
        let mut radii: Vec<f64> = self.bond_list(rmax)[i].iter().map(|b| b.radius).collect();
        radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut shells: Vec<Shell> = Vec::new();
        for r in radii {
            match shells.last_mut() {
                Some(last) if (r - last.radius).abs() <= rel_tol * last.radius => {
                    last.count += 1;
                }
                _ => shells.push(Shell {
                    radius: r,
                    count: 1,
                }),
            }
        }
        shells
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primitive_cell_volumes_are_the_standard_ones() {
        let a = 1.0;
        assert!((Lattice::sc(a).volume() - a.powi(3)).abs() < 1e-12);
        assert!((Lattice::bcc(a).volume() - a.powi(3) / 2.0).abs() < 1e-12);
        assert!((Lattice::fcc(a).volume() - a.powi(3) / 4.0).abs() < 1e-12);
        assert!((Lattice::diamond(a).volume() - a.powi(3) / 4.0).abs() < 1e-12);
        let hcp = Lattice::hcp(a);
        let expected = 3.0_f64.sqrt() / 2.0 * a * a * IDEAL_C_OVER_A * a;
        assert!((hcp.volume() - expected).abs() < 1e-9);
    }

    #[test]
    fn nearest_neighbour_distances() {
        let a = 1.0;
        assert!((Lattice::sc(a).nn_distance() - a).abs() < 1e-9);
        assert!((Lattice::fcc(a).nn_distance() - a / 2.0_f64.sqrt()).abs() < 1e-9);
        assert!((Lattice::bcc(a).nn_distance() - 3.0_f64.sqrt() * a / 2.0).abs() < 1e-9);
        assert!((Lattice::diamond(a).nn_distance() - 3.0_f64.sqrt() * a / 4.0).abs() < 1e-9);
        assert!((Lattice::hcp(a).nn_distance() - a).abs() < 1e-9);
    }

    /// The same, with the lattice constant in metres. The zero-vector test
    /// has to be relative: an absolute epsilon of 1e-9 is larger than
    /// every bond in the structure and used to skip all of them, so
    /// `nn_distance` returned a distant shell and `bond_list` returned
    /// nothing at all.
    #[test]
    fn neighbours_in_si_units() {
        let a = 3.615e-10; // copper
        let rel = |got: f64, want: f64| (got - want).abs() / want;
        assert!(rel(Lattice::sc(a).nn_distance(), a) < 1e-12);
        assert!(rel(Lattice::fcc(a).nn_distance(), a / 2.0_f64.sqrt()) < 1e-12);
        assert!(rel(Lattice::bcc(a).nn_distance(), 3.0_f64.sqrt() * a / 2.0) < 1e-12);
        assert!(rel(Lattice::diamond(a).nn_distance(), 3.0_f64.sqrt() * a / 4.0) < 1e-12);
        assert!(rel(Lattice::hcp(a).nn_distance(), a) < 1e-12);

        // The first four shells of face-centred cubic, in metres, must be
        // the same 12, 6, 24, 12 as at a = 1.
        let shells = Lattice::fcc(a).shells(0, 1.45 * a, 1e-6);
        let counts: Vec<usize> = shells.iter().map(|s| s.count).collect();
        assert_eq!(counts, vec![12, 6, 24, 12]);
        assert!(rel(shells[0].radius, a / 2.0_f64.sqrt()) < 1e-12);
        assert_eq!(Lattice::bcc(a).shells(0, 1.2 * a, 1e-6)[0].count, 8);
        assert_eq!(Lattice::diamond(a).shells(0, 0.5 * a, 1e-6)[0].count, 4);
    }
}
