//! Explicit functional-inequality constants on the (2π)³ torus, in the
//! coefficient norms the solver uses (box means: ‖f‖² = Σ_k |f̂_k|²,
//! ‖∇f‖² = Σ |k|²|f̂_k|², ‖Δf‖² = Σ |k|⁴|f̂_k|²), for mean-zero fields.
//!
//! Agmon: ‖f‖_∞ ≤ Σ_k |f̂_k| ≤ C_A ‖∇f‖^{1/2} ‖Δf‖^{1/2}.
//! Proof: split at |k| = K and use Cauchy–Schwarz on each part,
//!   Σ|f̂| ≤ √A(K) ‖∇f‖ + √B(K) ‖Δf‖,  A(K) = Σ_{0<|k|≤K} |k|⁻²,
//!   B(K) = Σ_{|k|>K} |k|⁻⁴.
//! With A(K) ≤ αK and B(K) ≤ β/K for all K > 0, choosing K optimally gives
//! C_A = 2 (αβ)^{1/4}. α and β are computed from exact lattice sums up to a
//! radius K_max plus rigorous integral bounds for the remainder.

/// Lattice sums S_p(K) = Σ_{0<|k|≤K} |k|^{-p} on a grid of radii, by exact
/// enumeration of k ∈ Z³ with |k| ≤ k_max. Returns sorted (|k|², count).
fn shells(k_max: i64) -> Vec<(i64, u64)> {
    let mut counts = std::collections::BTreeMap::new();
    for a in -k_max..=k_max {
        for b in -k_max..=k_max {
            for c in -k_max..=k_max {
                let r2 = a * a + b * b + c * c;
                if r2 > 0 && r2 <= k_max * k_max {
                    *counts.entry(r2).or_insert(0u64) += 1;
                }
            }
        }
    }
    counts.into_iter().collect()
}

/// Rigorous (α, β) and C_A.
///
/// Lattice points are compared with integrals over their unit cubes: for x
/// in the cube of k, | |x| − |k| | ≤ √3/2, hence for |k| ≥ R
/// |k|⁻ᵖ ≤ (1 + √3/(2R))ᵖ |x|⁻ᵖ, and the cubes of points with
/// R < |k| ≤ K lie in the shell R − √3/2 < |x| < K + √3/2.
///
/// * α = sup_K A(K)/K. A is a step function, so on each interval between
///   shell radii the sup is at the left end; those are enumerated exactly
///   up to R = K_max. Beyond: A(K) ≤ A(R) + c₂·4π(K − R + √3) with
///   c₂ = (1 + √3/(2R))², so A(K)/K ≤ max((A(R) + c₂4π(√3 − R))/R + c₂4π,
///   c₂4π) — the bound is monotone in K, its sup is at K = R or K → ∞.
/// * β = sup_K K·B(K). B is a decreasing step function; on each interval
///   the sup is approached at the right end. The enumerated part uses the
///   exact partial sums plus the tail bound below. For K ≥ R:
///   K·B(K) ≤ c₄(K)·4πK/(K − √3/2), c₄ = (1 + √3/(2K))⁴, decreasing in K,
///   so its sup is at K = R.
pub fn agmon() -> (f64, f64, f64) {
    let k_max: i64 = 60;
    let sh = shells(k_max);
    let r = k_max as f64;
    let h = 3f64.sqrt() / 2.0;
    let four_pi = 4.0 * std::f64::consts::PI;
    let c4 = (1.0 + h / r).powi(4);
    let c2 = (1.0 + h / r).powi(2);
    // tail Σ_{|k|>R} |k|⁻⁴ ≤ c₄ ∫_{|x|>R−√3/2} |x|⁻⁴ dx
    let tail4 = c4 * four_pi / (r - h);
    let total4: f64 = sh.iter().map(|&(r2, n)| n as f64 / (r2 * r2) as f64).sum::<f64>() + tail4;
    let (mut a_sum, mut b_rem) = (0.0, total4);
    let (mut alpha, mut beta): (f64, f64) = (0.0, total4); // K < 1: A = 0, K·B ≤ B(0)
    for (j, &(r2, n)) in sh.iter().enumerate() {
        let rj = (r2 as f64).sqrt();
        a_sum += n as f64 / r2 as f64;
        b_rem -= n as f64 / (r2 * r2) as f64;
        alpha = alpha.max(a_sum / rj);
        let r_next = sh.get(j + 1).map(|&(q, _)| (q as f64).sqrt()).unwrap_or(r);
        beta = beta.max(r_next * b_rem);
    }
    alpha = alpha.max((a_sum + c2 * four_pi * (h * 2.0 - r)) / r + c2 * four_pi).max(c2 * four_pi);
    beta = beta.max(c4 * four_pi * r / (r - h));
    (alpha, beta, 2.0 * (alpha * beta).powf(0.25))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The inequality must hold for random coefficient vectors.
    #[test]
    fn agmon_constant_holds_on_random_fields() {
        let (alpha, beta, ca) = agmon();
        assert!(alpha > 4.0 * std::f64::consts::PI * 0.9 && beta > 4.0 * std::f64::consts::PI * 0.9);
        let mut s = 0x1234_5678u64;
        let mut rnd = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s >> 11) as f64 / (1u64 << 53) as f64
        };
        for trial in 0..200 {
            let kc = 1 + trial % 12;
            let (mut l1, mut g2, mut d2) = (0.0, 0.0, 0.0);
            for a in -16i64..=16 {
                for b in -16i64..=16 {
                    for c in -16i64..=16 {
                        let k2 = (a * a + b * b + c * c) as f64;
                        if k2 == 0.0 {
                            continue;
                        }
                        let amp = rnd() * (-(k2.sqrt() / kc as f64)).exp();
                        l1 += amp;
                        g2 += k2 * amp * amp;
                        d2 += k2 * k2 * amp * amp;
                    }
                }
            }
            assert!(l1 <= ca * g2.powf(0.25) * d2.powf(0.25), "trial {trial}");
        }
    }
}
