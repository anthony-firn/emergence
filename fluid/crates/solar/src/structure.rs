//! Published Model S background, not an evolved/calibrated solar model.
use gravity::{Spherical, Vec3};
pub const MODEL_S_RADIUS: f64 = 6.959_894_677e8;
pub const MODEL_S_MASS: f64 = 1.989e30;
pub const MODEL_S_G: f64 = 6.67232e-11;
pub const MODEL_S_SHA256: &str = "65ecb920ed81b6b41f733cb8ab6f8c30941f7c743b0b6fec831de30e9a7322cc";
#[derive(Clone, Copy, Debug)]
pub struct Layer {
    pub radius: f64,
    pub sound_speed: f64,
    pub density: f64,
    pub pressure: f64,
    pub gamma1: f64,
    pub temperature: f64,
}
pub struct Profile {
    pub layers: Vec<Layer>,
}
impl Profile {
    /// Convert the basic public Model S table from cgs into SI, retaining its
    /// original radius/mass convention rather than silently fitting nominal R.
    pub fn model_s() -> Result<Self, &'static str> {
        let mut layers = vec![];
        for line in include_str!("../data/model_s.txt").lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let values: Result<Vec<f64>, _> = line.split_whitespace().map(str::parse).collect();
            let v = values.map_err(|_| "invalid Model S table")?;
            if v.len() != 6 || v.iter().any(|x| !x.is_finite()) {
                return Err("invalid Model S row");
            }
            layers.push(Layer {
                radius: v[0] * MODEL_S_RADIUS,
                sound_speed: v[1] * 0.01,
                density: v[2] * 1000.0,
                pressure: v[3] * 0.1,
                gamma1: v[4],
                temperature: v[5],
            });
        }
        layers.reverse();
        if layers.len() < 2
            || layers[0].radius != 0.0
            || layers.windows(2).any(|p| p[1].radius <= p[0].radius)
        {
            return Err("invalid solar background grid");
        }
        Ok(Self { layers })
    }
    /// Linear interpolation of published thermodynamic knots; no EOS inversion.
    pub fn at(&self, r: f64) -> Result<Layer, &'static str> {
        if !r.is_finite()
            || r < 0.0
            || r > self.layers.last().ok_or("empty stellar profile")?.radius
        {
            return Err("radius outside stellar profile");
        }
        let i = self.layers.partition_point(|p| p.radius < r);
        if i == 0 {
            return Ok(self.layers[0]);
        }
        if i >= self.layers.len() {
            return Ok(*self.layers.last().unwrap());
        }
        let (a, b) = (self.layers[i - 1], self.layers[i]);
        let t = (r - a.radius) / (b.radius - a.radius);
        let mix = |x: f64, y: f64| x * (1.0 - t) + y * t;
        Ok(Layer {
            radius: r,
            sound_speed: mix(a.sound_speed, b.sound_speed),
            density: mix(a.density, b.density),
            pressure: mix(a.pressure, b.pressure),
            gamma1: mix(a.gamma1, b.gamma1),
            temperature: mix(a.temperature, b.temperature),
        })
    }
    /// Exact volume-average shell densities of the linearly interpolated data.
    /// Gravity uses modern pinned G, not Model S's historical G value.
    pub fn gravity(&self, center: Vec3) -> Result<Spherical, &'static str> {
        if self.layers.len() < 2 {
            return Err("empty stellar gravity profile");
        }
        let edges = self.layers.iter().map(|p| p.radius).collect();
        let density = self
            .layers
            .windows(2)
            .map(|p| {
                let (a, b) = (p[0].radius, p[1].radius);
                let mid = 0.5 * (a + b);
                let rho_mid = 0.5 * (p[0].density + p[1].density);
                (p[0].density * a * a + 4.0 * rho_mid * mid * mid + p[1].density * b * b)
                    / (2.0 * (a * a + a * b + b * b))
            })
            .collect();
        Spherical::new(edges, density, center)
    }
}
