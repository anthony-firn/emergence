//! Full-support spectral importance sampling of the analytic CIE observer.
use crate::optics::cie_xyz;
pub struct WavelengthSampler {
    cdf: [f64; 400],
    weights: [f64; 400],
    sum: f64,
}
impl Default for WavelengthSampler {
    fn default() -> Self {
        let weights = std::array::from_fn(|i| {
            let xyz = cie_xyz((380.0 + i as f64 + 0.5) * 1e-9);
            // A positive floor retains support even in the observer's tails.
            xyz.into_iter().sum::<f64>() + 0.05
        });
        let mut sum = 0.0;
        let cdf = weights.map(|w| {
            sum += w;
            sum
        });
        Self { cdf, weights, sum }
    }
}
impl WavelengthSampler {
    /// Piecewise constant PDF [m^-1]. No approximation to the rendered
    /// spectrum: lambda and its exact sampling density are returned together.
    pub fn sample(&self, u: f64) -> (f64, f64) {
        assert!(u > 0.0 && u < 1.0);
        let target = u * self.sum;
        let i = self.cdf.partition_point(|&c| c < target).min(399);
        let prior = if i == 0 { 0.0 } else { self.cdf[i - 1] };
        (
            (380.0 + i as f64 + (target - prior) / self.weights[i]) * 1e-9,
            self.weights[i] / (self.sum * 1e-9),
        )
    }
}
