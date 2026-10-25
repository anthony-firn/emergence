//! Optical materials with dispersion.
//!
//! Glasses use the Schott Sellmeier form `n² − 1 = Σ Bᵢλ²/(λ² − Cᵢ)` (λ in
//! µm). Ocular media use the indices of the Escudero-Sanz & Navarro (JOSA A
//! 16, 1881 (1999)) wide-angle eye at 458, 543, 589.3 and 632.8 nm (as
//! tabulated in the visisipy eye-model package), interpolated by a
//! least-squares Cauchy law `n = A + B/λ² + C/λ⁴`. The resulting
//! longitudinal chromatic aberration is a prediction, checked against the
//! measured human curve of Atchison & Smith (JOSA A 22, 29 (2005)).
//! (A first attempt with water dispersion in every medium, after the
//! reduced "chromatic eye" of Thibos et al. 1992, gave only ~75 % of the
//! measured LCA in this four-surface eye: the crystalline lens disperses
//! more than water.)

/// An optical medium.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Material {
    Vacuum,
    /// Air at 15 °C, 101 325 Pa (Ciddor 1996 simplified: n − 1 ≈ 2.77e-4
    /// across the visible; dispersion neglected, estimate).
    Air,
    Sellmeier { b: [f64; 3], c: [f64; 3] },
    /// Constant index plus the water dispersion relative to a reference
    /// wavelength.
    WaterLike { n_ref: f64, lambda_ref: f64 },
    /// Cauchy law `A + B/λ² + C/λ⁴` (λ in µm).
    Cauchy { a: f64, b: f64, c: f64 },
}

impl Material {
    pub fn index(&self, lambda: f64) -> f64 {
        match *self {
            Material::Vacuum => 1.0,
            Material::Air => 1.000_277,
            Material::Sellmeier { b, c } => sellmeier(b, c, lambda),
            Material::WaterLike { n_ref, lambda_ref } => n_ref + water(lambda) - water(lambda_ref),
            Material::Cauchy { a, b, c } => {
                let l2 = (lambda * 1e6).powi(2);
                a + b / l2 + c / (l2 * l2)
            }
        }
    }

    /// Schott N-BK7 (Schott optical glass data sheet).
    pub fn n_bk7() -> Self {
        Material::Sellmeier { b: [1.039_612_12, 0.231_792_344, 1.010_469_45], c: [0.006_000_698_67, 0.020_017_914_4, 103.560_653] }
    }
    /// Schott N-BAF10.
    pub fn n_baf10() -> Self {
        Material::Sellmeier { b: [1.585_149_5, 0.143_559_385, 1.085_212_69], c: [0.009_266_812_82, 0.042_448_980_5, 105.613_573] }
    }
    /// Schott SF10 (the lead-containing glass named in the Thorlabs
    /// AC254-050-A Zemax file).
    pub fn sf10() -> Self {
        Material::Sellmeier { b: [1.616_259_77, 0.259_229_334, 1.077_623_17], c: [0.012_753_455_9, 0.058_198_395_4, 116.607_68] }
    }

    /// Ocular media of the Escudero-Sanz & Navarro (1999) eye.
    pub fn cornea() -> Self {
        cauchy_fit(&[1.3828, 1.3777, 1.3760, 1.3747])
    }
    pub fn aqueous() -> Self {
        cauchy_fit(&[1.3445, 1.3391, 1.3374, 1.3360])
    }
    pub fn crystalline_lens() -> Self {
        cauchy_fit(&[1.4292, 1.4222, 1.4200, 1.4183])
    }
    pub fn vitreous() -> Self {
        cauchy_fit(&[1.3428, 1.3377, 1.3360, 1.3347])
    }
}

/// Wavelengths (µm) of the Navarro ocular index table.
const NAVARRO_UM: [f64; 4] = [0.458, 0.543, 0.5893, 0.6328];

/// Least-squares Cauchy coefficients through indices at [`NAVARRO_UM`].
fn cauchy_fit(n: &[f64; 4]) -> Material {
    // normal equations for basis (1, λ⁻², λ⁻⁴)
    let basis = |l: f64| [1.0, l.powi(-2), l.powi(-4)];
    let mut m = [[0.0; 3]; 3];
    let mut r = [0.0; 3];
    for (l, ni) in NAVARRO_UM.iter().zip(n) {
        let b = basis(*l);
        for i in 0..3 {
            r[i] += b[i] * ni;
            for j in 0..3 {
                m[i][j] += b[i] * b[j];
            }
        }
    }
    // Gaussian elimination (3×3, well conditioned after scaling is fine here)
    for col in 0..3 {
        let piv = (col..3).max_by(|&a, &b| m[a][col].abs().total_cmp(&m[b][col].abs())).unwrap();
        m.swap(col, piv);
        r.swap(col, piv);
        for row in col + 1..3 {
            let f = m[row][col] / m[col][col];
            for k in col..3 {
                m[row][k] -= f * m[col][k];
            }
            r[row] -= f * r[col];
        }
    }
    let mut x = [0.0; 3];
    for row in (0..3).rev() {
        x[row] = (r[row] - (row + 1..3).map(|k| m[row][k] * x[k]).sum::<f64>()) / m[row][row];
    }
    Material::Cauchy { a: x[0], b: x[1], c: x[2] }
}

fn sellmeier(b: [f64; 3], c: [f64; 3], lambda: f64) -> f64 {
    let l2 = (lambda * 1e6).powi(2);
    (1.0 + (0..3).map(|i| b[i] * l2 / (l2 - c[i])).sum::<f64>()).sqrt()
}

/// Distilled water at 20 °C, Daimon & Masumura (2007) four-term Sellmeier
/// (refractiveindex.info "Daimon-20.0C"), valid 0.182–1.129 µm.
pub fn water(lambda: f64) -> f64 {
    sellmeier4(
        [5.684_027_565e-1, 1.726_177_391e-1, 2.086_189_578e-2, 1.130_748_688e-1],
        [5.101_829_712e-3, 1.821_153_936e-2, 2.620_722_293e-2, 1.069_792_721e1],
        lambda,
    )
}

fn sellmeier4(b: [f64; 4], c: [f64; 4], lambda: f64) -> f64 {
    let l2 = (lambda * 1e6).powi(2);
    (1.0 + (0..4).map(|i| b[i] * l2 / (l2 - c[i])).sum::<f64>()).sqrt()
}

/// Abbe number `V_d = (n_d − 1)/(n_F − n_C)` (d 587.56 nm, F 486.13 nm,
/// C 656.27 nm).
pub fn abbe(m: &Material) -> f64 {
    (m.index(587.56e-9) - 1.0) / (m.index(486.13e-9) - m.index(656.27e-9))
}
