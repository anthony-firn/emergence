//! Poisson–Nernst–Planck (PNP) for a symmetric z:z electrolyte next to a
//! charged wall, and the Gouy–Chapman / Debye–Hückel benchmarks.
//!
//! The emergence ladder for charged interfaces has two coarse rungs above the
//! resolved PNP equations:
//! 1. **Gouy–Chapman**: at equilibrium the ion fluxes J± = 0, so the ions obey
//!    Boltzmann statistics c± = n₀ e^{∓u} with u = zeψ/kT, and PNP collapses to
//!    the Poisson–Boltzmann (PB) equation
//!        ε ε₀ ψ'' = 2 z e n₀ sinh(zeψ/kT),   u'' = λ_D⁻² sinh u.
//!    This integrates once in closed form: the Gouy–Chapman profile
//!        tanh(u/4) = tanh(u₀/4) e^{−x/λ_D}.
//! 2. **Debye–Hückel**: for |u| ≪ 1 the PB equation linearises,
//!    u'' = λ_D⁻² u, so u = u₀ e^{−x/λ_D} and the surface charge is linear in
//!    the potential.
//!
//! The validity indicator of rung 1 is |zeψ/kT|; of rung 2 as well. The fold
//! back to PNP (rung 0) is a Newton solve of the nonlinear PB equation, which
//! is the exact no-flux steady state of PNP.
//!
//! Material constants:
//! - ε₀ = 8.854 187 8128e-12 F/m (CODATA 2018);
//! - e = 1.602 176 634e-19 C (CODATA 2018);
//! - k_B = 1.380 649e-23 J/K (CODATA 2018);
//! - water ε_r = 78.5 at 298 K (Malmberg & Maryott 1956 give 78.36 at 25 °C;
//!   78.5 is the usual rounded value used in electrokinetics).

pub const EPS0: f64 = 8.854_187_812_8e-12;
pub const E_CHARGE: f64 = 1.602_176_634e-19;
pub const KB: f64 = 1.380_649e-23;
pub const NA: f64 = 6.022_140_76e23;
pub const EPS_WATER: f64 = 78.5;
pub const MU_WATER: f64 = 8.90e-4; // Pa s, 298 K (IAPWS/CRC)
pub const T0: f64 = 298.15;

/// Debye length λ_D = sqrt(ε ε₀ kT / (2 z² e² n₀)), n₀ the number density of
/// each sign species. `c_molar` in mol/L.
pub fn debye_length(eps_r: f64, t: f64, c_molar: f64, z: f64) -> f64 {
    let n0 = c_molar * 1.0e3 * NA; // mol/m³ → 1/m³
    (eps_r * EPS0 * KB * t / (2.0 * n0 * z * z * E_CHARGE * E_CHARGE)).sqrt()
}

/// Dimensionless potential u = z e ψ / (k_B T).
pub fn u_of_psi(psi: f64, t: f64, z: f64) -> f64 {
    z * E_CHARGE * psi / (KB * t)
}

pub fn psi_of_u(u: f64, t: f64, z: f64) -> f64 {
    u * KB * t / (z * E_CHARGE)
}

/// Gouy–Chapman equilibrium profile ψ(x) for a single wall at x = 0 held at
/// ψ(0) = ψ₀, decaying to zero at +∞.
pub fn gouy_chapman(x: f64, lambda: f64, psi0: f64, t: f64, z: f64) -> f64 {
    let u0 = u_of_psi(psi0, t, z);
    let u = 4.0 * ((u0 / 4.0).tanh() * (-x / lambda).exp()).atanh();
    psi_of_u(u, t, z)
}

/// Debye–Hückel (linear) profile ψ(x) = ψ₀ e^{−x/λ_D}.
pub fn debye_huckel(x: f64, lambda: f64, psi0: f64) -> f64 {
    psi0 * (-x / lambda).exp()
}

/// Grahame equation: surface charge density σ (C/m²) of a symmetric z:z
/// electrolyte at wall potential ψ₀:
///     σ = √(8 n₀ ε ε₀ k_B T) sinh( z e ψ₀ / 2 k_B T ).
pub fn grahame_sigma(eps_r: f64, t: f64, c_molar: f64, z: f64, psi0: f64) -> f64 {
    let n0 = c_molar * 1.0e3 * NA;
    (8.0 * n0 * eps_r * EPS0 * KB * t).sqrt() * (z * E_CHARGE * psi0 / (2.0 * KB * t)).sinh()
}

/// Surface charge from the resolved profile: σ = −ε ε₀ dψ/dx|₀ (Gauss).
pub fn sigma_from_profile(psi: &[f64], x: &[f64], eps_r: f64) -> f64 {
    // one-sided second-order derivative at the wall
    let h = x[1] - x[0];
    let d = (-3.0 * psi[0] + 4.0 * psi[1] - psi[2]) / (2.0 * h);
    -eps_r * EPS0 * d
}

/// Newton solve of the dimensionless PB equation u'' = λ⁻² sinh u on [0, L]
/// with u(0) = u0 and u(L) = 0 (the PNP no-flux steady state next to a single
/// charged wall). Returns (x, ψ). Small tridiagonal Newton system.
pub fn pb_solve(eps_r: f64, t: f64, c_molar: f64, z: f64, psi0: f64, l: f64, n: usize) -> (Vec<f64>, Vec<f64>) {
    pb_solve_dirichlet(eps_r, t, c_molar, z, psi0, 0.0, l, n)
}

/// Newton solve of the dimensionless PB equation with fixed potentials at both
/// ends, u(0)=ψ_left, u(L)=ψ_right (voltages, not yet dimensionless).
pub fn pb_solve_dirichlet(
    eps_r: f64,
    t: f64,
    c_molar: f64,
    z: f64,
    psi_left: f64,
    psi_right: f64,
    l: f64,
    n: usize,
) -> (Vec<f64>, Vec<f64>) {
    let lambda = debye_length(eps_r, t, c_molar, z);
    let h = l / n as f64;
    let ul = u_of_psi(psi_left, t, z);
    let ur = u_of_psi(psi_right, t, z);
    let invl2 = 1.0 / (lambda * lambda);
    let x: Vec<f64> = (0..=n).map(|i| i as f64 * h).collect();
    let mut u: Vec<f64> = x
        .iter()
        .map(|&xi| ul + (ur - ul) * xi / l)
        .collect();
    u[0] = ul;
    u[n] = ur;
    let h2 = h * h;
    for _ in 0..80 {
        // r_i = (u_{i-1} - 2 u_i + u_{i+1})/h² - invl2 sinh u_i
        // J diag = -2/h² - invl2 cosh u_i, off = 1/h²
        let m = n - 1; // interior unknowns u_1..u_{n-1}
        let mut a = vec![0.0; m]; // sub
        let mut b = vec![0.0; m]; // diag
        let mut c = vec![0.0; m]; // super
        let mut rhs = vec![0.0; m];
        for i in 1..n {
            let k = i - 1;
            let r = (u[i - 1] - 2.0 * u[i] + u[i + 1]) / h2 - invl2 * u[i].sinh();
            b[k] = -2.0 / h2 - invl2 * u[i].cosh();
            if k > 0 {
                a[k] = 1.0 / h2;
            }
            if k < m - 1 {
                c[k] = 1.0 / h2;
            }
            rhs[k] = -r;
        }
        let du = thomas(&a, &b, &c, &rhs);
        let mut maxd = 0.0f64;
        for i in 1..n {
            u[i] += du[i - 1];
            maxd = maxd.max(du[i - 1].abs());
        }
        if maxd < 1e-13 {
            break;
        }
    }
    let psi: Vec<f64> = u.iter().map(|&ui| psi_of_u(ui, t, z)).collect();
    (x, psi)
}

/// Thomas algorithm for a tridiagonal system (a sub, b diag, c super).
pub fn thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = d.len();
    let mut cp = vec![0.0; n];
    let mut dp = vec![0.0; n];
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    for i in 1..n {
        let m = b[i] - a[i] * cp[i - 1];
        cp[i] = if i < n - 1 { c[i] / m } else { 0.0 };
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m;
    }
    let mut x = vec![0.0; n];
    x[n - 1] = dp[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }
    x
}

/// Charge density ρ_e = −2 z e n₀ sinh u for a symmetric z:z electrolyte.
pub fn charge_density(u: f64, c_molar: f64, z: f64) -> f64 {
    let n0 = c_molar * 1.0e3 * NA;
    -2.0 * z * E_CHARGE * n0 * u.sinh()
}
