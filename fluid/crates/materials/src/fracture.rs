//! Fracture across scales: from bond breaking to continuum fracture
//! toughness.
//!
//! This module is the fracture rung of the materials ladder
//! ([`analysis/emergence_ladder.md`](../../../../../analysis/emergence_ladder.md)).
//! The organising pattern is the same everywhere:
//!
//! - **coarse model** - linear-elastic fracture mechanics (LEFM): the
//!   crack-tip field is universal, the energy release rate is
//!   `G = K^2/E'`, and the crack grows when `G = 2 γ_s` (Griffith 1921).
//! - **validity indicator** - the ratio of the process-zone size to the
//!   sample (K-dominance) size, and the Rice (1992) competition between
//!   cleavage and dislocation emission,
//!   `γ_us / γ_s`.
//! - **fold / unfold** - the atomistic rung is a 2D triangular LJ sheet
//!   with a pre-crack loaded by the exact finite-crack Westergaard field;
//!   the closure for the coarse model is the cohesive-zone
//!   traction-separation law derived from the atomistic binding curve,
//!   with peak stress and work of separation `2 γ_s`.
//!
//! ## 1. The Griffith fold
//!
//! For a brittle crystal `G_c = 2 γ_s` where `γ_s` is the surface energy
//! ([`crate::defects`]); in plane strain the mode-I toughness is
//!
//! ```text
//! K_IC = sqrt(2 γ_s E / (1 - ν²)).
//! ```
//!
//! For an anisotropic crystal the Irwin matrix of Stroh elasticity
//! replaces `E/(1-ν²)`: with the near-tip field
//! `σ_{2i} = K_i / sqrt(2 π r)` the energy release rate is
//! `G = ½ Kᵀ H K` and [`stroh_irwin_matrix`] returns `H`.
//!
//! ## 2. The atomistic crack and lattice trapping
//!
//! A 2D triangular LJ sheet (nearest-neighbour bonds) with a row of atoms
//! removed on `x < 0` is loaded by the exact finite-crack mode-I
//! Westergaard field imposed on the outer boundary. The thermodynamic
//! advance threshold follows from [`critical_k_vce`] (the applied K at
//! which the virtual-crack-extension energy release rate crosses zero);
//! the discrete tip bond ruptures and heals at the local spinodal over a
//! finite window of loads, the lattice-trap precursor measured by
//! [`lattice_trapping`]. A full load-controlled barrier calculation is not
//! implemented (see `analysis/fracture.md`).
//!
//! ## 3. Indicators and closures
//!
//! [`process_zone_radius`] is the Irwin estimate `r_p = K²/(2π σ_y²)`;
//! [`rice_ductile_brittle`] is the Rice (1992) competition
//! `G_cleave = 2 γ_s` vs `G_emit ≈ 8 γ_us`; [`universal_binding`] is the
//! Rose-Smith-Ferrante universal binding energy, and
//! [`cohesive_zone_universal`] gives a traction-separation law whose work
//! equals `2 γ_s`.

use crate::lattice::{norm, Bond, Lattice, V3};
use crate::potential::PairPotential;

// ---------------------------------------------------------------------------
// 1. Griffith fold
// ---------------------------------------------------------------------------

/// Plane-strain Griffith toughness `K_IC = sqrt(2 γ_s E / (1 - ν²))`.
pub fn griffith_plane_strain(gamma_s: f64, e: f64, nu: f64) -> f64 {
    (2.0 * gamma_s * e / (1.0 - nu * nu)).sqrt()
}

/// Plane-stress Griffith toughness `K_IC = sqrt(2 γ_s E)`.
pub fn griffith_plane_stress(gamma_s: f64, e: f64) -> f64 {
    (2.0 * gamma_s * e).sqrt()
}

/// Two-dimensional (monolayer) Griffith toughness, same as plane stress
/// with the 2D Young modulus: `K_IC = sqrt(2 γ_s E_2D)`, where `γ_s` is an
/// energy per unit crack length and `E_2D` an energy per area.
pub fn griffith_2d(gamma_s: f64, e_2d: f64) -> f64 {
    (2.0 * gamma_s * e_2d).sqrt()
}

// ---------------------------------------------------------------------------
// Complex arithmetic (self-contained; the crate has no dependencies)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq)]
struct Cx {
    re: f64,
    im: f64,
}

impl Cx {
    fn new(re: f64, im: f64) -> Self {
        Cx { re, im }
    }
    fn zero() -> Self {
        Cx::new(0.0, 0.0)
    }
    fn add(self, o: Cx) -> Cx {
        Cx::new(self.re + o.re, self.im + o.im)
    }
    fn sub(self, o: Cx) -> Cx {
        Cx::new(self.re - o.re, self.im - o.im)
    }
    fn mul(self, o: Cx) -> Cx {
        Cx::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
    fn scale(self, s: f64) -> Cx {
        Cx::new(self.re * s, self.im * s)
    }
    fn div(self, o: Cx) -> Cx {
        let d = o.re * o.re + o.im * o.im;
        Cx::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }
    fn abs(self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }
}

type C3 = [[Cx; 3]; 3];

fn c3_mul(a: &C3, b: &C3) -> C3 {
    let mut c = [[Cx::zero(); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut s = Cx::zero();
            for k in 0..3 {
                s = s.add(a[i][k].mul(b[k][j]));
            }
            c[i][j] = s;
        }
    }
    c
}

fn c3_inverse(a: &C3) -> C3 {
    // Adjugate / determinant.
    let mut inv = [[Cx::zero(); 3]; 3];
    let c00 = a[1][1].mul(a[2][2]).sub(a[1][2].mul(a[2][1]));
    let c01 = a[1][2].mul(a[2][0]).sub(a[1][0].mul(a[2][2]));
    let c02 = a[1][0].mul(a[2][1]).sub(a[1][1].mul(a[2][0]));
    let det = a[0][0].mul(c00).add(a[0][1].mul(c01)).add(a[0][2].mul(c02));
    // cofactor matrix, then transpose
    let cof = [
        [
            c00,
            c01,
            c02,
        ],
        [
            a[0][2].mul(a[2][1]).sub(a[0][1].mul(a[2][2])),
            a[0][0].mul(a[2][2]).sub(a[0][2].mul(a[2][0])),
            a[0][1].mul(a[2][0]).sub(a[0][0].mul(a[2][1])),
        ],
        [
            a[0][1].mul(a[1][2]).sub(a[0][2].mul(a[1][1])),
            a[0][2].mul(a[1][0]).sub(a[0][0].mul(a[1][2])),
            a[0][0].mul(a[1][1]).sub(a[0][1].mul(a[1][0])),
        ],
    ];
    for i in 0..3 {
        for j in 0..3 {
            inv[i][j] = cof[j][i].div(det);
        }
    }
    inv
}

/// A real degree-`n` polynomial, coefficients `c[0] + c[1] p + ...`.
fn poly_mul(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut c = vec![0.0; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            c[i + j] += x * y;
        }
    }
    c
}

fn poly_add(a: &[f64], b: &[f64]) -> Vec<f64> {
    let n = a.len().max(b.len());
    let mut c = vec![0.0; n];
    for i in 0..n {
        c[i] = a.get(i).copied().unwrap_or(0.0) + b.get(i).copied().unwrap_or(0.0);
    }
    c
}

fn poly_scale(a: &[f64], s: f64) -> Vec<f64> {
    a.iter().map(|x| x * s).collect()
}

/// Roots of a real polynomial by Durand-Kerner (Weierstrass) iteration.
/// Returns the complex roots. Repeated roots (the isotropic limit) are the
/// hard case; the caller should validate away from degeneracy.
fn poly_roots(coeffs: &[f64]) -> Vec<Cx> {
    // Strip leading near-zero coefficients (highest order).
    let mut c = coeffs.to_vec();
    while c.len() > 1 && c.last().unwrap().abs() < 1e-14 {
        c.pop();
    }
    let n = c.len() - 1;
    if n == 0 {
        return Vec::new();
    }
    let lead = *c.last().unwrap();
    let norm: Vec<f64> = c.iter().map(|x| x / lead).collect();
    // Initial roots on a circle (0.4 + 0.9i)^k, avoiding symmetry axes.
    let seed = Cx::new(0.4, 0.9);
    let mut roots = Vec::with_capacity(n);
    let mut z = Cx::new(1.0, 0.0);
    for _ in 0..n {
        roots.push(z);
        z = z.mul(seed);
    }
    let eval = |x: Cx| -> Cx {
        let mut acc = Cx::zero();
        for &k in norm.iter().rev() {
            acc = acc.mul(x).add(Cx::new(k, 0.0));
        }
        acc
    };
    for _ in 0..2000 {
        let mut max_step = 0.0_f64;
        let old = roots.clone();
        for i in 0..n {
            let mut denom = Cx::new(1.0, 0.0);
            for j in 0..n {
                if j != i {
                    denom = denom.mul(old[i].sub(old[j]));
                }
            }
            if denom.abs() < 1e-300 {
                continue;
            }
            let step = eval(old[i]).div(denom);
            roots[i] = old[i].sub(step);
            max_step = max_step.max(step.abs());
        }
        if max_step < 1e-13 {
            break;
        }
    }
    roots
}

/// Compute the 3x3 matrices `A = [a_α]`, `B = [b_α]` of the Stroh
/// formalism for a crack with faces on `x2 = 0` and front along `x3`.
/// Only eigenvectors with `Im(p) > 0` are used. The eigenvectors are
/// normalised so that `2 a_α · b_α = 1`.
fn stroh_ab(c: &[[f64; 6]; 6]) -> (C3, C3) {
    use crate::elastic::voigt;
    // Q_{ik} = C_{i1k1}, R_{ik} = C_{i1k2}, T_{ik} = C_{i2k2} (0-based
    // coordinates: 0 = x1, 1 = x2, 2 = x3).
    let mut q = [[0.0; 3]; 3];
    let mut r = [[0.0; 3]; 3];
    let mut t = [[0.0; 3]; 3];
    for i in 0..3 {
        for k in 0..3 {
            q[i][k] = c[voigt(i, 0)][voigt(k, 0)];
            r[i][k] = c[voigt(i, 0)][voigt(k, 1)];
            t[i][k] = c[voigt(i, 1)][voigt(k, 1)];
        }
    }
    // Sextic det[Q + p (R + Rᵀ) + p² T] = 0, built exactly with degree-2
    // polynomial entries.
    let entry = |m0: f64, m1: f64, m2: f64| -> Vec<f64> { vec![m0, m1, m2] };
    let mut a: [[Vec<f64>; 3]; 3] = Default::default();
    for i in 0..3 {
        for k in 0..3 {
            a[i][k] = entry(q[i][k], r[i][k] + r[k][i], t[i][k]);
        }
    }
    let det = poly_add(
        &poly_add(
            &poly_mul(
                &a[0][0],
                &poly_add(&poly_mul(&a[1][1], &a[2][2]), &poly_scale(&poly_mul(&a[1][2], &a[2][1]), -1.0)),
            ),
            &poly_scale(
                &poly_mul(&a[0][1], &poly_add(&poly_mul(&a[1][0], &a[2][2]), &poly_scale(&poly_mul(&a[1][2], &a[2][0]), -1.0))),
                -1.0,
            ),
        ),
        &poly_mul(
            &a[0][2],
            &poly_add(&poly_mul(&a[1][0], &a[2][1]), &poly_scale(&poly_mul(&a[1][1], &a[2][0]), -1.0)),
        ),
    );
    let roots = poly_roots(&det);
    let mut upper: Vec<Cx> = roots.into_iter().filter(|p| p.im > 1e-9).collect();
    // A cubic-symmetric crystal on a symmetry cut can put a root on the
    // real axis; keep at most three and pad if degenerate.
    upper.truncate(3);
    let mut a_mat: C3 = [[Cx::zero(); 3]; 3];
    let mut b_mat: C3 = [[Cx::zero(); 3]; 3];
    for (col, &p) in upper.iter().enumerate() {
        // M = Q + p(R+Rᵀ) + p²T.
        let mut m: C3 = [[Cx::zero(); 3]; 3];
        for i in 0..3 {
            for k in 0..3 {
                let m0 = Cx::new(q[i][k], 0.0);
                let m1 = Cx::new(r[i][k] + r[k][i], 0.0).mul(p);
                let m2 = Cx::new(t[i][k], 0.0).mul(p).mul(p);
                m[i][k] = m0.add(m1).add(m2);
            }
        }
        // Null vector from cross products of rows.
        let cross = |u: [Cx; 3], v: [Cx; 3]| -> [Cx; 3] {
            [
                u[1].mul(v[2]).sub(u[2].mul(v[1])),
                u[2].mul(v[0]).sub(u[0].mul(v[2])),
                u[0].mul(v[1]).sub(u[1].mul(v[0])),
            ]
        };
        let cands = [
            cross(m[0], m[1]),
            cross(m[1], m[2]),
            cross(m[0], m[2]),
        ];
        let mut best = cands[0];
        let mut best_n = 0.0;
        for cand in cands {
            let nn = cand.iter().map(|x| x.abs() * x.abs()).sum::<f64>().sqrt();
            if nn > best_n {
                best_n = nn;
                best = cand;
            }
        }
        if best_n < 1e-14 {
            continue;
        }
        let aa = [
            best[0].div(Cx::new(best_n, 0.0)),
            best[1].div(Cx::new(best_n, 0.0)),
            best[2].div(Cx::new(best_n, 0.0)),
        ];
        // b = (Rᵀ + p T) a.
        let mut bb = [Cx::zero(); 3];
        for i in 0..3 {
            let mut s = Cx::zero();
            for k in 0..3 {
                let rt = Cx::new(r[k][i], 0.0);
                let tt = Cx::new(t[i][k], 0.0).mul(p);
                s = s.add(rt.add(tt).mul(aa[k]));
            }
            bb[i] = s;
        }
        // Normalise so that 2 (a · b) = 1 (no conjugation).
        let ab = aa[0].mul(bb[0]).add(aa[1].mul(bb[1])).add(aa[2].mul(bb[2]));
        if ab.abs() < 1e-300 {
            continue;
        }
        let s = Cx::new(1.0, 0.0).div(ab.mul(Cx::new(2.0, 0.0))).abs().sqrt();
        for k in 0..3 {
            a_mat[k][col] = aa[k].scale(s);
            b_mat[k][col] = bb[k].scale(s);
        }
    }
    (a_mat, b_mat)
}

/// Stroh Irwin matrix `H` for a crack on `x2 = 0` with front along `x3`:
/// the energy release rate for a stress-intensity vector
/// `K = (K_II, K_I, K_III)` (components of `σ_{2i} = K_i/sqrt(2πr)`) is
/// `G = ½ Kᵀ H K`. `Re(i A B⁻¹)` under the Stroh normalisation
/// `2 a_α·b_α = 1` (Stroh 1958; Sih, Paris & Irwin 1965; Ting 1996).
pub fn stroh_irwin_matrix(c: &[[f64; 6]; 6]) -> [[f64; 3]; 3] {
    let (a, b) = stroh_ab(c);
    let binv = c3_inverse(&b);
    let ab = c3_mul(&a, &binv);
    // i * ab
    let mut h = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let z = ab[i][j].mul(Cx::new(0.0, 1.0));
            h[i][j] = z.re;
        }
    }
    // Symmetrise (the construction is symmetric up to round-off).
    for i in 0..3 {
        for j in i + 1..3 {
            let v = 0.5 * (h[i][j] + h[j][i]);
            h[i][j] = v;
            h[j][i] = v;
        }
    }
    h
}

/// Mode-I toughness from the Stroh Irwin matrix: `K_Ic = sqrt(4 γ_s / H_22)`.
pub fn griffith_anisotropic_kic(c: &[[f64; 6]; 6], gamma_s: f64) -> f64 {
    let h = stroh_irwin_matrix(c);
    assert!(h[1][1] > 0.0, "Irwin matrix must be positive definite");
    (4.0 * gamma_s / h[1][1]).sqrt()
}

// ---------------------------------------------------------------------------
// 2. Atomistic 2D triangular sheet with a pre-crack
// ---------------------------------------------------------------------------

/// A finite 2D triangular-lattice sheet with nearest-neighbour bonds and a
/// half-row of atoms removed to make a mode-I pre-crack.
#[derive(Clone, Debug)]
pub struct CrackedSheet {
    /// Reference positions (the crack tip is at the origin).
    pub atoms: Vec<V3>,
    /// Per-atom neighbour bonds (both directions present).
    pub bonds: Vec<Vec<Bond>>,
    /// True where the atom is free to relax (false = K-field boundary).
    pub free: Vec<bool>,
    /// Nearest-neighbour spacing.
    pub a0: f64,
    /// Number of atoms removed on the crack plane (crack length in sites).
    pub crack_sites: i32,
    /// Crack half-length `a` (used by the centre-crack field).
    pub crack_half: f64,
    /// Half-widths of the atomic block in lattice steps.
    pub nx: i32,
    pub ny: i32,
    /// Index of the tip atom (first intact atom on the crack plane).
    pub tip_atom: usize,
}

fn triangle_pos(i: i32, j: i32, a0: f64) -> V3 {
    let h = 3.0_f64.sqrt() / 2.0;
    [a0 * (i as f64 + 0.5 * j as f64), a0 * h * j as f64, 0.0]
}

/// Build a triangular sheet with `i in -nx..=nx`, `j in -ny..=ny` and
/// remove the crack-plane row `j = 0` for `i in -crack_sites..=0`.
/// Nearest-neighbour bonds only (`rmax` must be below `sqrt(3) a0`).
pub fn build_cracked_sheet(
    nx: i32,
    ny: i32,
    a0: f64,
    rmax: f64,
    crack_sites: i32,
    _pot: &dyn PairPotential,
) -> CrackedSheet {
    let mut atoms = Vec::new();
    let mut index = std::collections::HashMap::new();
    for j in -ny..=ny {
        for i in -nx..=nx {
            let remove = j == 0 && i >= -crack_sites && i <= 0;
            if remove {
                continue;
            }
            index.insert((i, j), atoms.len());
            atoms.push(triangle_pos(i, j, a0));
        }
    }
    let n = atoms.len();
    let mut bonds = vec![Vec::new(); n];
    for j in -ny..=ny {
        for i in -nx..=nx {
            let Some(&a) = index.get(&(i, j)) else {
                continue;
            };
            // Six nearest neighbours.
            let neigh = [
                (i + 1, j),
                (i - 1, j),
                (i, j + 1),
                (i - 1, j + 1),
                (i, j - 1),
                (i + 1, j - 1),
            ];
            for &(ii, jj) in &neigh {
                if let Some(&b) = index.get(&(ii, jj)) {
                    let r = [
                        atoms[b][0] - atoms[a][0],
                        atoms[b][1] - atoms[a][1],
                        atoms[b][2] - atoms[a][2],
                    ];
                    let rad = norm(r);
                    if rad > 1e-9 && rad <= rmax {
                        bonds[a].push(Bond {
                            r,
                            radius: rad,
                            i: a,
                            j: b,
                            cell: [0, 0, 0],
                        });
                    }
                }
            }
        }
    }
    // Boundary atoms are held at the K-field displacement.
    let mut free = vec![true; n];
    for j in -ny..=ny {
        for i in -nx..=nx {
            if let Some(&a) = index.get(&(i, j)) {
                if i == -nx || i == nx || j == -ny || j == ny {
                    free[a] = false;
                }
            }
        }
    }
    let tip_atom = *index.get(&(1, 0)).expect("tip atom must exist");
    CrackedSheet {
        atoms,
        bonds,
        free,
        a0,
        crack_sites,
        crack_half: 0.5 * crack_sites as f64 * a0,
        nx,
        ny,
        tip_atom,
    }
}

/// Isotropic mode-I plane-stress K-field displacement at `(x, y)` relative
/// to the crack tip, with Young modulus `e` and Poisson ratio `nu`.
///
/// `u_x = (K/2μ) √(r/2π) cos(θ/2)(κ − cosθ)`,
/// `u_y = (K/2μ) √(r/2π) sin(θ/2)(κ − cosθ)`,
/// `κ = (3 − ν)/(1 + ν)`. This is the near-tip asymptotic field.
pub fn mode_i_field(x: f64, y: f64, k: f64, e: f64, nu: f64) -> V3 {
    let mu = e / (2.0 * (1.0 + nu));
    let kappa = (3.0 - nu) / (1.0 + nu);
    let r = (x * x + y * y).sqrt();
    let theta = y.atan2(x);
    let amp = k / (2.0 * mu) * (r / (2.0 * std::f64::consts::PI)).sqrt();
    let c = (theta * 0.5).cos();
    let s = (theta * 0.5).sin();
    let cth = theta.cos();
    [
        amp * c * (kappa - cth),
        amp * s * (kappa - cth),
        0.0,
    ]
}

/// Principal complex square root.
fn csqrt(w: Cx) -> Cx {
    let r = w.abs();
    let re = ((r + w.re) * 0.5).max(0.0).sqrt();
    let im_mag = ((r - w.re) * 0.5).max(0.0).sqrt();
    Cx::new(re, if w.im < 0.0 { -im_mag } else { im_mag })
}

/// Exact plane-stress mode-I displacement field of a **finite** centre
/// crack of half-length `a` under a remote tension `sigma`, written in
/// terms of the tip stress intensity factor `K = sigma sqrt(pi a)`.
///
/// Westergaard:
/// `u_x = (1/E)[(1−ν) Re Z̄ − (1+ν) y Im Z]`,
/// `u_y = (1/E)[2 Im Z̄ − (1+ν) y Re Z]`,
/// with `Z = sigma z/√(z²−a²)`, `Z̄ = ∫Z dz = sigma √(z²−a²)`, and `z`
/// measured from the crack centre. The crack occupies `|x'| < a` and the
/// right tip is at `z' = a`. See Westergaard (1939); Anderson, *Fracture
/// Mechanics*.
pub fn center_crack_field(x: f64, y: f64, a: f64, k: f64, e: f64, nu: f64) -> V3 {
    // Shift so the crack centre is at the origin; the right tip is at x' = a.
    let z = Cx::new(x + a, y);
    let sigma = k / (std::f64::consts::PI * a).sqrt();
    let root = csqrt(z.mul(z).sub(Cx::new(a * a, 0.0)));
    let zbar = root.scale(sigma);
    let zfunc = z.mul(Cx::new(sigma, 0.0)).div(root);
    let ux = ((1.0 - nu) * zbar.re - (1.0 + nu) * y * zfunc.im) / e;
    let uy = (2.0 * zbar.im - (1.0 + nu) * y * zfunc.re) / e;
    [ux, uy, 0.0]
}

fn max_force_masked(f: &[V3], free: &[[bool; 3]]) -> f64 {
    let mut m = 0.0_f64;
    for (v, fr) in f.iter().zip(free) {
        for d in 0..3 {
            if fr[d] {
                m = m.max(v[d].abs());
            }
        }
    }
    m
}

fn dot_masked(a: &[V3], b: &[V3], free: &[[bool; 3]]) -> f64 {
    let mut s = 0.0;
    for ((x, y), fr) in a.iter().zip(b).zip(free) {
        for d in 0..3 {
            if fr[d] {
                s += x[d] * y[d];
            }
        }
    }
    s
}

/// Polak-Ribiere conjugate gradients with an Armijo backtracking line
/// search, with a per-atom, per-component free mask. Returns
/// `(displacements, energy, max free force)`.
pub fn minimize_masked(
    sheet: &CrackedSheet,
    pot: &dyn PairPotential,
    mut u: Vec<V3>,
    load: &[V3],
    free: &[[bool; 3]],
    max_steps: usize,
    ftol: f64,
) -> (Vec<V3>, f64, f64) {
    use crate::defects::energy_forces;
    let n = u.len();
    let total = |u: &[V3]| -> f64 {
        let (e, _) = energy_forces(&sheet.bonds, pot, u);
        let mut t = e;
        for i in 0..n {
            for d in 0..3 {
                if free[i][d] {
                    t -= load[i][d] * u[i][d];
                }
            }
        }
        t
    };
    let force = |u: &[V3]| -> Vec<V3> {
        let (_, mut f) = energy_forces(&sheet.bonds, pot, u);
        for i in 0..n {
            for d in 0..3 {
                if free[i][d] {
                    f[i][d] += load[i][d];
                } else {
                    f[i][d] = 0.0;
                }
            }
        }
        f
    };
    let mut e = total(&u);
    let f = force(&u);
    // g is the gradient (-force); dir is the descent direction (-g = force).
    let mut g: Vec<V3> = f.iter().map(|v| [-v[0], -v[1], -v[2]]).collect();
    let mut dir = f.clone();
    let mut alpha = 0.1_f64;
    let mut fmax = max_force_masked(&f, free);
    let mut steps = 0;
    while steps < max_steps && fmax > ftol {
        let gd = dot_masked(&g, &dir, free);
        if gd >= 0.0 {
            dir = g.iter().map(|v| [-v[0], -v[1], -v[2]]).collect();
        }
        let gd = dot_masked(&g, &dir, free);
        if gd.abs() < 1e-300 {
            break;
        }
        let u0 = u.clone();
        let mut a = alpha.max(1e-12);
        let mut ok = false;
        for _ in 0..80 {
            let mut ut = u0.clone();
            for i in 0..n {
                for d in 0..3 {
                    if free[i][d] {
                        ut[i][d] += a * dir[i][d];
                    }
                }
            }
            let et = total(&ut);
            if et <= e + 1e-4 * a * gd {
                u = ut;
                e = et;
                alpha = a;
                ok = true;
                break;
            }
            a *= 0.5;
        }
        if !ok {
            break;
        }
        steps += 1;
        let f_new = force(&u);
        let g_new: Vec<V3> = f_new.iter().map(|v| [-v[0], -v[1], -v[2]]).collect();
        let num = dot_masked(&g_new, &g_new, free) - dot_masked(&g_new, &g, free);
        let gg = dot_masked(&g, &g, free);
        let beta = if gg > 0.0 { (num / gg).max(0.0) } else { 0.0 };
        for i in 0..n {
            for d in 0..3 {
                if free[i][d] {
                    dir[i][d] = -g_new[i][d] + beta * dir[i][d];
                } else {
                    dir[i][d] = 0.0;
                }
            }
        }
        g = g_new;
        fmax = max_force_masked(&f_new, free);
    }
    (u, e, fmax)
}

/// Polak-Ribiere conjugate gradients with an Armijo backtracking line
/// search, with a per-atom free mask. Returns
/// `(displacements, energy, max free-atom force)`.
pub fn minimize_atoms(
    sheet: &CrackedSheet,
    pot: &dyn PairPotential,
    u: Vec<V3>,
    load: &[V3],
    max_steps: usize,
    ftol: f64,
) -> (Vec<V3>, f64, f64) {
    let free: Vec<[bool; 3]> = sheet.free.iter().map(|&fr| [fr; 3]).collect();
    minimize_masked(sheet, pot, u, load, &free, max_steps, ftol)
}

/// Initial displacement field = K-field imposed on the boundary, and the
/// same field used as the interior starting guess (helps the relaxation).
pub fn k_field_displacements(
    sheet: &CrackedSheet,
    k: f64,
    e: f64,
    nu: f64,
) -> Vec<V3> {
    sheet
        .atoms
        .iter()
        .map(|x| center_crack_field(x[0], x[1], sheet.crack_half, k, e, nu))
        .collect()
}

/// Warm-start displacements for a K increment: `u_prev` relaxed at
/// `k_prev`, plus the elastic increment `(k - k_prev)` times the
/// normalised K-field. Boundary atoms are set exactly to the new field.
pub fn k_field_increment(
    sheet: &CrackedSheet,
    u_prev: &[V3],
    k_prev: f64,
    k: f64,
    e: f64,
    nu: f64,
) -> Vec<V3> {
    let mut u = Vec::with_capacity(sheet.atoms.len());
    for (i, x) in sheet.atoms.iter().enumerate() {
        let f = center_crack_field(x[0], x[1], sheet.crack_half, 1.0, e, nu);
        if sheet.free[i] {
            u.push([
                u_prev[i][0] + (k - k_prev) * f[0],
                u_prev[i][1] + (k - k_prev) * f[1],
                0.0,
            ]);
        } else {
            u.push([k * f[0], k * f[1], 0.0]);
        }
    }
    u
}

/// Deformed bond lengths and the number of bonds longer than the LJ
/// spinodal `r_s = (26/7)^(1/6) σ` (where `φ''` vanishes), counting only
/// bonds whose two ends are at least `margin` sites from the **outer free
/// surface** (so boundary-layer bond stretching is not mistaken for a
/// crack-tip event).
pub fn broken_bond_count(
    sheet: &CrackedSheet,
    _pot: &dyn PairPotential,
    u: &[V3],
    r_break: f64,
    margin: i32,
) -> usize {
    let a0 = sheet.a0;
    let h = 3.0_f64.sqrt() / 2.0;
    let m = margin;
    let mut count = 0;
    for i in 0..sheet.atoms.len() {
        if !sheet.free[i] {
            continue;
        }
        let xi = sheet.atoms[i];
        let jj = (xi[1] / (a0 * h)).round() as i32;
        let ii = (xi[0] / a0 - 0.5 * jj as f64).round() as i32;
        if ii <= -sheet.nx + m || ii >= sheet.nx - m || jj <= -sheet.ny + m || jj >= sheet.ny - m {
            continue;
        }
        for b in &sheet.bonds[i] {
            if b.j <= i {
                continue;
            }
            let xj = sheet.atoms[b.j];
            let j2 = (xj[1] / (a0 * h)).round() as i32;
            let i2 = (xj[0] / a0 - 0.5 * j2 as f64).round() as i32;
            if i2 <= -sheet.nx + m || i2 >= sheet.nx - m || j2 <= -sheet.ny + m || j2 >= sheet.ny - m {
                continue;
            }
            let s = [
                b.r[0] + u[b.j][0] - u[i][0],
                b.r[1] + u[b.j][1] - u[i][1],
                b.r[2] + u[b.j][2] - u[i][2],
            ];
            if norm(s) > r_break {
                count += 1;
            }
        }
    }
    count
}

/// Result of a quasistatic loading/unloading sweep.
#[derive(Clone, Debug)]
pub struct LatticeTrapping {
    /// Smallest `K` at which the tip bond ruptures during loading
    /// (exceeds the LJ spinodal). This is the **local bond instability**,
    /// not the thermodynamic advance threshold (see
    /// [`critical_k_vce`]).
    pub k_plus: f64,
    /// Largest `K` at which the ruptured tip bond heals during unloading.
    pub k_minus: f64,
    /// Griffith prediction `sqrt(2 gamma_s E_2D)`.
    pub k_griffith: f64,
    /// `k_plus / k_minus` (the metastable bond-rupture window).
    pub ratio: f64,
    /// Broken-bond counts over the loading history: `(K, count)`.
    pub loading: Vec<(f64, usize)>,
    /// Broken-bond counts over the unloading history.
    pub unloading: Vec<(f64, usize)>,
}

/// Load/unload a cracked sheet quasistatically and record the **tip-bond
/// rupture** window, an operational proxy for lattice trapping:
///
/// 1. increase `K` in `n_steps` increments, relaxing at each step; record
///    the broken-bond count. `k_plus` is the first `K` at which the count
///    exceeds the pristine value (the tip bond stretched past the LJ
///    spinodal `(26/7)^{1/6} sigma`);
/// 2. starting from the loaded state, decrease `K`; `k_minus` is the
///    largest `K` at which the count returns to the pristine value.
///
/// The local bond instability can lie well **below** the thermodynamic
/// advance threshold `critical_k_vce`, where the total energy still
/// favours the pristine crack; the bond is then metastable and heals on
/// unloading. Both are reported in `analysis/fracture.md`.
#[allow(clippy::too_many_arguments)]
pub fn lattice_trapping(
    sheet: &CrackedSheet,
    pot: &dyn PairPotential,
    e_2d: f64,
    nu: f64,
    gamma_s: f64,
    k_max: f64,
    n_steps: usize,
    r_break: f64,
    margin: i32,
) -> LatticeTrapping {
    let pristine = broken_bond_count(sheet, pot, &vec![[0.0; 3]; sheet.atoms.len()], r_break, margin);
    // --- loading ---
    let mut u = k_field_displacements(sheet, 0.0, e_2d, nu);
    let load = vec![[0.0; 3]; u.len()];
    let mut loading = Vec::new();
    let mut k_plus = f64::INFINITY;
    let mut u_plus = u.clone();
    let mut k_prev = 0.0;
    for step in 0..=n_steps {
        let k = k_max * step as f64 / n_steps as f64;
        u = k_field_increment(sheet, &u, k_prev, k, e_2d, nu);
        let (un, _, _) = minimize_atoms(sheet, pot, u, &load, 3000, 1e-6);
        u = un;
        k_prev = k;
        let count = broken_bond_count(sheet, pot, &u, r_break, margin);
        loading.push((k, count));
        if count > pristine && k_plus.is_infinite() {
            k_plus = k;
            u_plus = u.clone();
        }
    }
    // --- unloading from the last loaded state ---
    let mut u = u_plus.clone();
    let mut unloading = Vec::new();
    let mut k_minus = 0.0_f64;
    let n_unload = n_steps;
    let mut k_prev = k_plus;
    for step in 0..=n_unload {
        let k = k_plus * (n_unload - step) as f64 / n_unload as f64;
        u = k_field_increment(sheet, &u, k_prev, k, e_2d, nu);
        let (un, _, _) = minimize_atoms(sheet, pot, u, &load, 3000, 1e-6);
        u = un;
        k_prev = k;
        let count = broken_bond_count(sheet, pot, &u, r_break, margin);
        unloading.push((k, count));
        if count <= pristine {
            k_minus = k;
            break;
        }
    }
    let ratio = if k_minus > 0.0 { k_plus / k_minus } else { f64::INFINITY };
    LatticeTrapping {
        k_plus,
        k_minus,
        k_griffith: griffith_2d(gamma_s, e_2d),
        ratio,
        loading,
        unloading,
    }
}

// ---------------------------------------------------------------------------
// 2b. Atomistic energy release rate by virtual crack extension
// ---------------------------------------------------------------------------

/// Relax a cracked sheet with the boundary displacement field of a
/// **fixed-grip** test: the centre-crack field of half-length `a_field`
/// and tip K-factor `k` is imposed on the outer boundary; the interior
/// (including the crack faces) relaxes. Returns the potential energy.
#[allow(clippy::too_many_arguments)]
pub fn relaxed_crack_energy(
    pot: &dyn PairPotential,
    nx: i32,
    ny: i32,
    a0: f64,
    rmax: f64,
    crack_sites: i32,
    a_field: f64,
    k: f64,
    e_2d: f64,
    nu: f64,
    max_steps: usize,
    ftol: f64,
) -> f64 {
    let sheet = build_cracked_sheet(nx, ny, a0, rmax, crack_sites, pot);
    let load = vec![[0.0; 3]; sheet.atoms.len()];
    let u0: Vec<V3> = sheet
        .atoms
        .iter()
        .map(|x| center_crack_field(x[0], x[1], a_field, k, e_2d, nu))
        .collect();
    let (_, en, _) = minimize_atoms(&sheet, pot, u0, &load, max_steps, ftol);
    en
}

/// Atomistic energy release rate by **virtual crack extension** at fixed
/// grip: advance the crack by one lattice site `a0` (crack length
/// `+a0`, two new surfaces of length `a0`) while keeping the boundary
/// displacement field of the reference crack `a_ref`,
/// `G = (E(a) - E(a + a0))/a0`. Positive `G` means the extended crack
/// lowers the total energy (surface energy included in `E`), so the crack
/// advances spontaneously; the Griffith threshold is `G = 0`, i.e.
/// elastic driving force `= 2 γ_s`.
#[allow(clippy::too_many_arguments)]
pub fn energy_release_rate_vce(
    pot: &dyn PairPotential,
    nx: i32,
    ny: i32,
    a0: f64,
    rmax: f64,
    crack_sites: i32,
    anchor_sites: i32,
    k: f64,
    e_2d: f64,
    nu: f64,
    max_steps: usize,
    ftol: f64,
) -> f64 {
    let a_ref = anchor_sites as f64 * a0 / 2.0;
    let e0 = relaxed_crack_energy(pot, nx, ny, a0, rmax, crack_sites, a_ref, k, e_2d, nu, max_steps, ftol);
    let e1 = relaxed_crack_energy(pot, nx, ny, a0, rmax, crack_sites + 1, a_ref, k, e_2d, nu, max_steps, ftol);
    (e0 - e1) / a0
}

/// Critical applied K at which the atomistic energy release rate crosses
/// zero (the crack advances). Bisection on `[lo, hi]`; `tol` is the
/// absolute tolerance in K. `anchor_sites` fixes the crack half-length
/// used for the boundary displacement field (the pristine tip for the
/// forward bound, the advanced tip for the backward bound).
#[allow(clippy::too_many_arguments)]
pub fn critical_k_vce(
    pot: &dyn PairPotential,
    nx: i32,
    ny: i32,
    a0: f64,
    rmax: f64,
    crack_sites: i32,
    anchor_sites: i32,
    e_2d: f64,
    nu: f64,
    lo: f64,
    hi: f64,
    tol: f64,
    max_steps: usize,
    ftol: f64,
) -> f64 {
    let mut lo = lo;
    let mut hi = hi;
    let g_lo = energy_release_rate_vce(pot, nx, ny, a0, rmax, crack_sites, anchor_sites, lo, e_2d, nu, max_steps, ftol);
    let g_hi = energy_release_rate_vce(pot, nx, ny, a0, rmax, crack_sites, anchor_sites, hi, e_2d, nu, max_steps, ftol);
    assert!(g_lo <= 0.0 && g_hi >= 0.0, "bracket must straddle G=0");
    while hi - lo > tol {
        let mid = 0.5 * (lo + hi);
        let g = energy_release_rate_vce(pot, nx, ny, a0, rmax, crack_sites, anchor_sites, mid, e_2d, nu, max_steps, ftol);
        if g > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Effective tip stress-intensity factor implied by a *net* virtual-crack-
/// extension energy release rate `G_net`: the elastic driving force is
/// `G_el = G_net + 2 γ_s` (the surface energy is already included in
/// `G_net`), and in plane stress `G_el = K_tip² / E_2D`. The value at
/// `G_net = 0` is `sqrt(2 γ_s E_2D) = K_G` by construction.
pub fn effective_tip_k(g_net: f64, gamma_s: f64, e_2d: f64) -> f64 {
    let g_el = (g_net + 2.0 * gamma_s).max(0.0);
    (e_2d * g_el).sqrt()
}

// ---------------------------------------------------------------------------
// 3. Indicators and closure
// ---------------------------------------------------------------------------

/// Irwin plane-stress process-zone radius `r_p = K²/(2 π σ_y²)`.
/// `sigma_y` is the tensile yield stress.
pub fn process_zone_radius(k: f64, sigma_y: f64) -> f64 {
    k * k / (2.0 * std::f64::consts::PI * sigma_y * sigma_y)
}

/// K-dominance / LEFM validity indicator: the coarse model is valid when
/// the K-dominance zone (a fraction of the smallest dimension `L`) is much
/// larger than the process zone, `r_p / L ≪ 1`.
pub fn lefm_valid(k: f64, sigma_y: f64, sample_size: f64, threshold: f64) -> bool {
    process_zone_radius(k, sigma_y) < threshold * sample_size
}

/// Rice (1992) brittle-vs-ductile competition at a sharp crack tip.
///
/// Cleavage requires `G = 2 γ_s`; emitting a dislocation from the tip
/// requires `G = 8 γ_us` (Rice's estimate for emission on the slip plane
/// through the tip). The crystal cleaves when
///
/// ```text
/// 2 γ_s < 8 γ_us   i.e.   γ_us / γ_s > 1/4,
/// ```
///
/// and emits a dislocation (ductile) otherwise.
#[derive(Clone, Copy, Debug)]
pub struct RiceCriterion {
    pub gamma_us: f64,
    pub gamma_s: f64,
    pub g_cleave: f64,
    pub g_emit: f64,
    pub ratio: f64,
    pub brittle: bool,
}

/// Evaluate the Rice (1992) criterion. `gamma_us` is the unstable
/// stacking-fault energy from [`crate::defects`].
pub fn rice_ductile_brittle(gamma_us: f64, gamma_s: f64) -> RiceCriterion {
    let g_cleave = 2.0 * gamma_s;
    let g_emit = 8.0 * gamma_us;
    RiceCriterion {
        gamma_us,
        gamma_s,
        g_cleave,
        g_emit,
        ratio: gamma_us / gamma_s,
        brittle: g_cleave < g_emit,
    }
}

/// Rose-Smith-Ferrante universal binding-energy relation (Rose, Smith &
/// Ferrante, *Phys. Rev. Lett.* **47**, 675 (1981)):
/// `E(δ) = -E_0 (1 + δ/λ) e^{−δ/λ}`, `δ` the separation, `E_0` the binding
/// energy and `λ` the scale (Thomas-Fermi screening length).
pub fn universal_binding(delta: f64, e0: f64, lambda: f64) -> f64 {
    -e0 * (1.0 + delta / lambda) * (-delta / lambda).exp()
}

/// Traction-separation law derived from a binding energy `E(δ)`:
/// `T(δ) = dE/dδ` (with the convention that tension is positive when
/// `E` is the energy relative to the separated state). The work of
/// separation is `∫_0^∞ T dδ = E(∞) - E(0)`.
pub fn traction_from_energy<F: Fn(f64) -> f64>(energy: F, delta: f64, h: f64) -> f64 {
    (energy(delta + h) - energy(delta - h)) / (2.0 * h)
}

/// Peak traction and work of separation of the universal binding law,
/// expressed through the surface energy `γ_s`. The atomistic
/// energy-separation curve has `E(∞) - E(0) = 2 γ_s` (two surfaces).
/// The UBER peak stress is `T_max = (2 γ_s/λ) e^{−1}` at `δ = λ`.
#[derive(Clone, Copy, Debug)]
pub struct CohesiveZone {
    pub gamma_s: f64,
    pub lambda: f64,
    pub work_of_separation: f64,
    pub t_max: f64,
    pub delta_max: f64,
}

/// Universal-binding cohesive-zone parameters for a given surface energy
/// and screening length: work of separation `2 γ_s`, peak traction
/// `(2 γ_s/λ) e^{−1}` at `δ = λ`.
pub fn cohesive_zone_universal(gamma_s: f64, lambda: f64) -> CohesiveZone {
    let work = 2.0 * gamma_s;
    let t_max = work / lambda * (-1.0_f64).exp();
    CohesiveZone {
        gamma_s,
        lambda,
        work_of_separation: work,
        t_max,
        delta_max: lambda,
    }
}

/// Work of separation of a numerically tabulated traction-separation law
/// by the trapezoidal rule.
pub fn cohesive_work(tractions: &[f64], deltas: &[f64]) -> f64 {
    assert_eq!(tractions.len(), deltas.len());
    let mut w = 0.0;
    for i in 1..tractions.len() {
        w += 0.5 * (tractions[i] + tractions[i - 1]) * (deltas[i] - deltas[i - 1]);
    }
    w
}

/// Surface energy (energy per unit area in 3D, per unit length in 2D) of a
/// central-force lattice by the broken-bond sum
/// `γ(n) = ½ ρ Σ_{r·n>0} (r·n) |φ(r)|` ([`crate::defects`]). This is a
/// thin wrapper so the 2D sheet and the 3D crystal use the same routine.
pub fn surface_energy_broken_bond(lat: &Lattice, pot: &dyn PairPotential, rmax: f64, normal: V3) -> f64 {
    crate::defects::broken_bond_surface_energy(lat, pot, rmax, normal)
}

/// 2D triangular nearest-neighbour lattice parameter at zero stress:
/// minimising `3 φ(r)` gives `r = 2^{1/6} σ` exactly for a 12-6
/// Lennard-Jones potential.
pub fn lj_nn_triangular_spacing() -> f64 {
    2.0_f64.powf(1.0 / 6.0)
}

/// Cauchy-Born in-plane elastic constants of a 2D triangular lattice
/// (no `z` bonds), returned as `(E_2D, ν_2D)` with `E_2D = (C11²−C12²)/C11`
/// and `ν_2D = C12/C11` (plane stress).
pub fn triangular_2d_moduli(pot: &dyn PairPotential, a0: f64, rmax: f64) -> (f64, f64) {
    // A 3D cell with unit thickness in z and only in-plane bonds: the
    // Cauchy-Born tensor entries are energies per unit area.
    let h = 3.0_f64.sqrt() / 2.0;
    let lat = Lattice::custom(
        [[a0, 0.0, 0.0], [0.5 * a0, h * a0, 0.0], [0.0, 0.0, 1.0]],
        vec![[0.0, 0.0, 0.0]],
    );
    let c = crate::elastic::elastic_tensor(&lat, pot, rmax);
    let c11 = c[0][0];
    let c12 = c[0][1];
    let e = (c11 * c11 - c12 * c12) / c11;
    (e, c12 / c11)
}

/// The lattice constant of the 2D sheet for a Lennard-Jones potential with
/// a nearest-neighbour cutoff: `2^{1/6} σ`. Provided as the exact
/// reference for [`triangular_2d_moduli`].
pub fn lj_2d_equilibrium() -> f64 {
    lj_nn_triangular_spacing()
}
