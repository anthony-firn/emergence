//! Jones calculus for fully polarized fields (paraxial: the field is
//! transverse, `(E_x, E_y)`).

use num_complex::Complex64 as Cx;

/// A Jones vector `(E_x, E_y)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Jones(pub [Cx; 2]);

impl Jones {
    pub fn horizontal() -> Self {
        Jones([Cx::new(1.0, 0.0), Cx::new(0.0, 0.0)])
    }
    pub fn vertical() -> Self {
        Jones([Cx::new(0.0, 0.0), Cx::new(1.0, 0.0)])
    }
    /// Linear polarization at angle `theta` from x.
    pub fn linear(theta: f64) -> Self {
        Jones([Cx::new(theta.cos(), 0.0), Cx::new(theta.sin(), 0.0)])
    }
    /// Circular (`+1` left, `−1` right, in the source's convention).
    pub fn circular(handedness: f64) -> Self {
        let s = std::f64::consts::FRAC_1_SQRT_2;
        Jones([Cx::new(s, 0.0), Cx::new(0.0, s * handedness.signum())])
    }
    pub fn intensity(&self) -> f64 {
        self.0[0].norm_sqr() + self.0[1].norm_sqr()
    }
    pub fn normalized(&self) -> Self {
        let n = self.intensity().sqrt();
        Jones([self.0[0] / n, self.0[1] / n])
    }
    pub fn scale(&self, c: Cx) -> Self {
        Jones([self.0[0] * c, self.0[1] * c])
    }
    pub fn add(&self, o: &Jones) -> Self {
        Jones([self.0[0] + o.0[0], self.0[1] + o.0[1]])
    }
}

/// A Jones matrix.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JonesMatrix(pub [[Cx; 2]; 2]);

impl JonesMatrix {
    pub fn identity() -> Self {
        let (o, z) = (Cx::new(1.0, 0.0), Cx::new(0.0, 0.0));
        JonesMatrix([[o, z], [z, o]])
    }
    /// Ideal linear polarizer with transmission axis at `theta` from x.
    pub fn polarizer(theta: f64) -> Self {
        let (c, s) = (theta.cos(), theta.sin());
        JonesMatrix([[Cx::new(c * c, 0.0), Cx::new(c * s, 0.0)], [Cx::new(c * s, 0.0), Cx::new(s * s, 0.0)]])
    }
    /// Linear retarder (fast axis at `theta`) of retardance `delta`
    /// (π/2: quarter-wave plate, π: half-wave plate).
    pub fn retarder(theta: f64, delta: f64) -> Self {
        let (c, s) = (theta.cos(), theta.sin());
        let e = Cx::from_polar(1.0, delta);
        let o = Cx::new(1.0, 0.0);
        // R(−θ) diag(1, e^{iδ}) R(θ)
        JonesMatrix([[o * c * c + e * s * s, (o - e) * c * s], [(o - e) * c * s, o * s * s + e * c * c]])
    }
    pub fn apply(&self, j: &Jones) -> Jones {
        let m = &self.0;
        Jones([m[0][0] * j.0[0] + m[0][1] * j.0[1], m[1][0] * j.0[0] + m[1][1] * j.0[1]])
    }
}
