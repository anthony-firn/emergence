//! Tensor-product meshes with per-element materials.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Material {
    Silicon,
    Oxide,
}

#[derive(Clone, Debug)]
pub struct Mesh {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// material of element (i, j) between nodes i..i+1, j..j+1, at
    /// index i·(ny − 1) + j
    pub elem: Vec<Material>,
}

impl Mesh {
    pub fn new(x: Vec<f64>, y: Vec<f64>, material: impl Fn(f64, f64) -> Material) -> Self {
        let mut elem = Vec::with_capacity((x.len() - 1) * (y.len() - 1));
        for i in 0..x.len() - 1 {
            for j in 0..y.len() - 1 {
                elem.push(material(0.5 * (x[i] + x[i + 1]), 0.5 * (y[j] + y[j + 1])));
            }
        }
        Mesh { x, y, elem }
    }
    /// One-dimensional mesh along x (a unit-height strip of two node rows:
    /// densities are per m², currents per m²).
    pub fn line(x: Vec<f64>, material: impl Fn(f64) -> Material) -> Self {
        Mesh::new(x, vec![0.0, 1.0], |x, _| material(x))
    }
    pub fn nx(&self) -> usize {
        self.x.len()
    }
    pub fn ny(&self) -> usize {
        self.y.len()
    }
    pub fn nodes(&self) -> usize {
        self.x.len() * self.y.len()
    }
    #[inline]
    pub fn node(&self, i: usize, j: usize) -> usize {
        i * self.y.len() + j
    }
    #[inline]
    pub fn ij(&self, k: usize) -> (usize, usize) {
        (k / self.y.len(), k % self.y.len())
    }
    pub fn material(&self, i: usize, j: usize) -> Material {
        self.elem[i * (self.y.len() - 1) + j]
    }
}

/// Node positions with local spacing interpolated linearly between
/// breakpoints (position, spacing); breakpoints are nodes.
pub fn graded(points: &[(f64, f64)]) -> Vec<f64> {
    let mut x = vec![points[0].0];
    for s in points.windows(2) {
        let ((a, ha), (b, hb)) = (s[0], s[1]);
        let l = b - a;
        // ∫ dx / h(x), h linear in x
        let total = if (hb - ha).abs() < 1e-12 * ha { l / ha } else { l / (hb - ha) * (hb / ha).ln() };
        let n = total.ceil().max(1.0) as usize;
        for k in 1..=n {
            let u = total * k as f64 / n as f64;
            let xx = if (hb - ha).abs() < 1e-12 * ha { a + u * ha } else { a + ha * ((u * (hb - ha) / l).exp() - 1.0) * l / (hb - ha) };
            x.push(if k == n { b } else { xx });
        }
    }
    x
}

/// Uniform nodes on [a, b].
pub fn uniform(a: f64, b: f64, n: usize) -> Vec<f64> {
    (0..=n).map(|i| a + (b - a) * i as f64 / n as f64).collect()
}
