//! Photon-counting image sensors.
//!
//! The expected number of photo-electrons in a pixel is
//! `N̄ = Σ_λ (t/(hc/λ)) ∫_pixel E_λ(x) dA · QE(λ)`. Photon arrivals of a
//! (coherent or thermal-but-multimode) beam are Poisson; detection with
//! probability `QE` thins a Poisson process into a Poisson process, so the
//! collected charge is exactly `Poisson(N̄ + dark·t)`. Charge saturates at
//! the full well, the readout adds Gaussian noise, and the ADC quantises
//! (`gain` electrons per digital number).
//!
//! **Pixel integration (a fold).** A pixel's signal is the irradiance
//! integrated over its active area. Where the pattern varies slowly on the
//! pixel scale, the centre value times the area is accurate; where it does
//! not, more points are needed. With the shortest variation scale `P` (for
//! fringes, the period), an `n`-point Gauss–Legendre rule per axis
//! integrates `e^{2πix/P}` over a width `ℓ` with relative error below
//! `(πℓ/P)^{2n}/(2n)!`; the smallest `n` meeting the tolerance is used.
//!
//! **Silicon quantum efficiency (physical model, estimates).** Light
//! entering through an anti-reflection coating (transmission `T`) is
//! absorbed with the absorption coefficient `α(λ)` of crystalline silicon
//! (Green, Sol. Energy Mater. Sol. Cells 92, 1305 (2008), 300 K table);
//! charge from a surface dead layer `d₀` recombines and only the depth
//! `d₀ … d₀ + W` is collected:
//! `QE = T [e^{−α d₀} − e^{−α (d₀+W)}]`.

use std::f64::consts::PI;

/// Planck constant × speed of light (J m), CODATA 2018 exact.
pub const HC: f64 = 6.626_070_15e-34 * 2.997_924_58e8;

/// Green (2008) silicon absorption coefficient at 300 K, (nm, cm⁻¹).
const SI_ALPHA: [(f64, f64); 39] = [
    (380.0, 2.93e5), (390.0, 1.50e5), (400.0, 9.52e4), (410.0, 6.74e4), (420.0, 5.00e4), (430.0, 3.92e4), (440.0, 3.11e4), (450.0, 2.55e4),
    (460.0, 2.10e4), (470.0, 1.72e4), (480.0, 1.48e4), (490.0, 1.27e4), (500.0, 1.11e4), (510.0, 9.70e3), (520.0, 8.80e3), (540.0, 7.05e3),
    (560.0, 5.78e3), (580.0, 4.88e3), (600.0, 4.14e3), (620.0, 3.52e3), (640.0, 3.04e3), (660.0, 2.58e3), (680.0, 2.21e3), (700.0, 1.90e3),
    (720.0, 1.66e3), (740.0, 1.42e3), (760.0, 1.19e3), (780.0, 1.01e3), (800.0, 8.50e2), (850.0, 5.35e2), (900.0, 3.06e2), (950.0, 1.57e2),
    (1000.0, 6.40e1), (1050.0, 1.63e1), (1100.0, 3.50e0), (1110.0, 2.0), (1120.0, 1.0), (1150.0, 0.1), (1200.0, 0.01),
];

/// Silicon absorption coefficient (1/m), log-linear interpolation of Green
/// (2008); the last four points beyond 1100 nm are an extrapolation
/// (estimate) towards the indirect band edge.
pub fn silicon_alpha(lambda: f64) -> f64 {
    let nm = lambda * 1e9;
    let t = &SI_ALPHA;
    if nm <= t[0].0 {
        return t[0].1 * 100.0;
    }
    for w in t.windows(2) {
        if nm <= w[1].0 {
            let f = (nm - w[0].0) / (w[1].0 - w[0].0);
            return (w[0].1.ln() * (1.0 - f) + w[1].1.ln() * f).exp() * 100.0;
        }
    }
    t[t.len() - 1].1 * 100.0
}

/// Quantum efficiency model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Qe {
    Constant(f64),
    /// Silicon photodiode: AR transmission, dead layer and collection depth (m).
    Silicon { transmission: f64, dead_layer: f64, depth: f64 },
}

impl Qe {
    /// A back-side illuminated CMOS sensor (estimates: T = 0.95,
    /// d₀ = 50 nm, W = 3 µm).
    pub fn bsi_cmos() -> Self {
        Qe::Silicon { transmission: 0.95, dead_layer: 50e-9, depth: 3e-6 }
    }

    pub fn at(&self, lambda: f64) -> f64 {
        match *self {
            Qe::Constant(q) => q,
            Qe::Silicon { transmission, dead_layer, depth } => {
                let a = silicon_alpha(lambda);
                transmission * ((-a * dead_layer).exp() - (-a * (dead_layer + depth)).exp())
            }
        }
    }
}

/// A monochrome image sensor in the observation plane.
#[derive(Clone, Debug, PartialEq)]
pub struct Sensor {
    pub pixels: [usize; 2],
    /// Pixel pitch (m).
    pub pitch: f64,
    /// Fraction of the pixel area that collects light (a centred square).
    pub fill_factor: f64,
    pub qe: Qe,
    /// Full-well capacity (e⁻).
    pub full_well: f64,
    /// Read noise (e⁻ rms).
    pub read_noise: f64,
    /// Dark current (e⁻/s per pixel).
    pub dark_current: f64,
    /// Conversion gain (e⁻ per digital number).
    pub gain: f64,
    pub bits: u32,
    /// Exposure time (s).
    pub exposure: f64,
    /// Sensor centre in the observation plane (m).
    pub center: [f64; 2],
}

impl Sensor {
    /// A typical BSI CMOS sensor (values are estimates of a mid-range
    /// scientific/machine-vision part): 3.45 µm pixels, 10 000 e⁻ full
    /// well, 2.5 e⁻ read noise, 12 bit.
    pub fn cmos(pixels: [usize; 2], pitch: f64, exposure: f64) -> Self {
        Self {
            pixels,
            pitch,
            fill_factor: 1.0,
            qe: Qe::bsi_cmos(),
            full_well: 10_000.0,
            read_noise: 2.5,
            dark_current: 2.0,
            gain: 10_000.0 / 4095.0,
            bits: 12,
            exposure,
            center: [0.0, 0.0],
        }
    }

    /// Centre of pixel `(i, j)` (`i` along x).
    pub fn pixel_center(&self, i: usize, j: usize) -> [f64; 2] {
        [
            self.center[0] + self.pitch * (i as f64 + 0.5 - 0.5 * self.pixels[0] as f64),
            self.center[1] + self.pitch * (j as f64 + 0.5 - 0.5 * self.pixels[1] as f64),
        ]
    }

    /// Gauss points per axis for a pattern whose shortest variation scale
    /// is `scale`, so the pixel integral is accurate to `tol` (relative).
    pub fn points_per_axis(&self, scale: f64, tol: f64) -> usize {
        let l = self.pitch * self.fill_factor.sqrt();
        let x = PI * l / scale.max(1e-300);
        let mut fact = 1.0f64; // (2n)!
        for n in 1..=16usize {
            fact *= ((2 * n - 1) * (2 * n)) as f64;
            if x.powi(2 * n as i32) / fact < tol {
                return n;
            }
        }
        16
    }

    /// Quadrature points (all pixels) and weights (area fractions) for an
    /// `n × n` Gauss rule over each pixel's active area.
    pub fn quadrature(&self, n: usize) -> (Vec<[f64; 2]>, Vec<f64>) {
        let (xg, wg) = gauss_legendre(n);
        let side = self.pitch * self.fill_factor.sqrt();
        let mut pts = Vec::with_capacity(self.pixels[0] * self.pixels[1] * n * n);
        let mut w = Vec::with_capacity(pts.capacity());
        for j in 0..self.pixels[1] {
            for i in 0..self.pixels[0] {
                let c = self.pixel_center(i, j);
                for (a, wa) in xg.iter().zip(&wg) {
                    for (b, wb) in xg.iter().zip(&wg) {
                        pts.push([c[0] + 0.5 * side * a, c[1] + 0.5 * side * b]);
                        w.push(0.25 * wa * wb * side * side);
                    }
                }
            }
        }
        (pts, w)
    }

    /// Tensor-grid form of [`Sensor::quadrature`]: the Gauss abscissae of
    /// every pixel along each axis (`pixels[a]·n` values) and the 1-D
    /// weights (fractions of the active side, summing to 1 per pixel).
    pub fn axis_quadrature(&self, n: usize) -> ([Vec<f64>; 2], Vec<f64>) {
        let (xg, wg) = gauss_legendre(n);
        let side = self.pitch * self.fill_factor.sqrt();
        let axis = |a: usize| -> Vec<f64> {
            (0..self.pixels[a])
                .flat_map(|i| {
                    let c = self.center[a] + self.pitch * (i as f64 + 0.5 - 0.5 * self.pixels[a] as f64);
                    xg.iter().map(move |x| c + 0.5 * side * x).collect::<Vec<_>>()
                })
                .collect()
        };
        ([axis(0), axis(1)], wg.iter().map(|w| 0.5 * w).collect())
    }

    /// Expected photo-electrons per pixel from spectral irradiance on the
    /// tensor grid of [`Sensor::axis_quadrature`] (index `i + nx j`).
    pub fn expected_electrons_grid(&self, n: usize, spectral: &[(f64, Vec<f64>)]) -> Vec<f64> {
        let (_, w) = self.axis_quadrature(n);
        let [px, py] = self.pixels;
        let nx = px * n;
        let area = self.pitch * self.pitch * self.fill_factor;
        let mut out = vec![0.0; px * py];
        for (lambda, e) in spectral {
            let k = self.exposure * lambda / HC * self.qe.at(*lambda) * area;
            for (q, v) in e.iter().enumerate() {
                let (gi, gj) = (q % nx, q / nx);
                out[gi / n + px * (gj / n)] += k * v * w[gi % n] * w[gj % n];
            }
        }
        out
    }

    /// Expected photo-electrons per pixel from spectral irradiance samples
    /// `(λ, E_λ at the quadrature points)` of an `n × n` rule.
    pub fn expected_electrons(&self, n: usize, spectral: &[(f64, Vec<f64>)]) -> Vec<f64> {
        let (_, w) = self.quadrature(n);
        let per = n * n;
        let npix = self.pixels[0] * self.pixels[1];
        let mut out = vec![0.0; npix];
        for (lambda, e) in spectral {
            let k = self.exposure * lambda / HC * self.qe.at(*lambda);
            for p in 0..npix {
                let s: f64 = (0..per).map(|q| e[p * per + q] * w[p * per + q]).sum();
                out[p] += k * s;
            }
        }
        out
    }

    /// One exposure: Poisson charge (+ dark), full-well clipping, read
    /// noise, quantisation. Returns digital numbers.
    pub fn expose(&self, expected: &[f64], rng: &mut Rng) -> Vec<u32> {
        let max = (1u64 << self.bits) as f64 - 1.0;
        expected
            .iter()
            .map(|&m| {
                let e = rng.poisson(m + self.dark_current * self.exposure).min(self.full_well);
                let v = (e + self.read_noise * rng.normal()) / self.gain;
                v.round().clamp(0.0, max) as u32
            })
            .collect()
    }
}

/// `n`-point Gauss–Legendre nodes and weights on [−1, 1] (Newton on P_n).
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    for i in 0..n.div_ceil(2) {
        let mut z = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut pp = 1.0;
        for _ in 0..100 {
            let (mut p1, mut p2) = (1.0, 0.0);
            for k in 1..=n {
                let p3 = p2;
                p2 = p1;
                p1 = ((2 * k - 1) as f64 * z * p2 - (k - 1) as f64 * p3) / k as f64;
            }
            pp = n as f64 * (z * p1 - p2) / (z * z - 1.0);
            let dz = p1 / pp;
            z -= dz;
            if dz.abs() < 1e-15 {
                break;
            }
        }
        x[i] = -z;
        x[n - 1 - i] = z;
        w[i] = 2.0 / ((1.0 - z * z) * pp * pp);
        w[n - 1 - i] = w[i];
    }
    (x, w)
}

/// Deterministic random numbers (SplitMix64), with exact Poisson and
/// normal deviates.
pub struct Rng(pub u64);

impl Rng {
    pub fn uniform(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        (((z ^ (z >> 31)) >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }

    pub fn normal(&mut self) -> f64 {
        let (u, v) = (self.uniform(), self.uniform());
        (-2.0 * u.ln()).sqrt() * (2.0 * PI * v).cos()
    }

    /// Exact Poisson deviate: multiplication method below 30, PTRS
    /// (Hörmann, Insurance Math. Econ. 12, 1993) above.
    pub fn poisson(&mut self, mean: f64) -> f64 {
        if mean <= 0.0 {
            return 0.0;
        }
        if mean < 30.0 {
            let l = (-mean).exp();
            let (mut k, mut p) = (0.0, 1.0);
            loop {
                p *= self.uniform();
                if p <= l {
                    return k;
                }
                k += 1.0;
            }
        }
        let slam = mean.sqrt();
        let loglam = mean.ln();
        let b = 0.931 + 2.53 * slam;
        let a = -0.059 + 0.02483 * b;
        let invalpha = 1.1239 + 1.1328 / (b - 3.4);
        let vr = 0.9277 - 3.6224 / (b - 2.0);
        loop {
            let u = self.uniform() - 0.5;
            let v = self.uniform();
            let us = 0.5 - u.abs();
            let k = ((2.0 * a / us + b) * u + mean + 0.43).floor();
            if us >= 0.07 && v <= vr {
                return k;
            }
            if k < 0.0 || (us < 0.013 && v > us) {
                continue;
            }
            if (v * invalpha / (a / (us * us) + b)).ln() <= -mean + k * loglam - ln_gamma(k + 1.0) {
                return k;
            }
        }
    }
}

/// `ln Γ(x)` (Lanczos, g = 7, n = 9; ~1e-15 relative).
pub fn ln_gamma(x: f64) -> f64 {
    const G: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        return (PI / (PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = G[0];
    let t = x + 7.5;
    for (i, g) in G.iter().enumerate().skip(1) {
        a += g / (x + i as f64);
    }
    0.5 * (2.0 * PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// Draw `n` single-photon detection positions (pixel indices) from the
/// expected-electron map (the detection probability density is ∝ |ψ|²
/// integrated over the pixel, times QE).
pub fn single_photon_events(expected: &[f64], n: usize, rng: &mut Rng) -> Vec<usize> {
    let mut cdf = Vec::with_capacity(expected.len());
    let mut s = 0.0;
    for &e in expected {
        s += e;
        cdf.push(s);
    }
    (0..n)
        .map(|_| {
            let u = rng.uniform() * s;
            cdf.partition_point(|&c| c < u).min(expected.len() - 1)
        })
        .collect()
}
