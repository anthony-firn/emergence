//! Spherical solid diffusion in an electrode particle, and the exact series
//! for a constant surface flux.
//!
//! ```text
//!   ∂c/∂t = (D/r²) ∂_r(r² ∂_r c),    ∂_r c(0,t)=0,   −D ∂_r c(R,t)=J₀
//! ```
//!
//! with `J₀` the outward surface flux (mol m⁻² s⁻¹). A finite-volume
//! discretisation on node-centred shells conserves moles exactly; the outer
//! face carries the applied flux and the inner face is a symmetry plane. The
//! exact benchmark is the series for the **constant-flux surface
//! concentration** (Crank, *The Mathematics of Diffusion*, 2nd ed., ch. 6),
//!
//! ```text
//!   c_s(t) = c₀ − 3J₀t/R − (2RJ₀/D) Σ_{n≥1} (1−e^{−D μ_n² t/R²})/μ_n²,
//!   tan μ_n = μ_n.
//! ```
//!
//! The volume average is the exact linear result `c̄ = c₀ − 3J₀t/R`.

/// Finite-volume spherical particle.
pub struct SphereParticle {
    pub nodes: usize,
    /// Particle radius (m).
    pub radius: f64,
    /// Solid diffusivity (m² s⁻¹).
    pub d: f64,
    /// Node radial positions `r_i = (i+½)Δr`.
    pub r: Vec<f64>,
    /// Node concentrations (mol m⁻³).
    pub c: Vec<f64>,
    dr: f64,
    vol: Vec<f64>,
    area: Vec<f64>,
}

impl SphereParticle {
    /// Build a particle with `nodes` shells, radius `r`, diffusivity `d` and
    /// uniform initial concentration `c0`.
    pub fn new(nodes: usize, radius: f64, d: f64, c0: f64) -> Self {
        let dr = radius / nodes as f64;
        let r: Vec<f64> = (0..nodes).map(|i| (i as f64 + 0.5) * dr).collect();
        let vol: Vec<f64> = (0..nodes)
            .map(|i| {
                let ri = i as f64 * dr;
                let ro = (i as f64 + 1.0) * dr;
                (4.0 / 3.0) * std::f64::consts::PI * (ro.powi(3) - ri.powi(3))
            })
            .collect();
        let area: Vec<f64> = (0..=nodes)
            .map(|i| 4.0 * std::f64::consts::PI * (i as f64 * dr).powi(2))
            .collect();
        Self { nodes, radius, d, r, c: vec![c0; nodes], dr, vol, area }
    }

    /// Surface concentration (outermost node, mol m⁻³).
    pub fn surface(&self) -> f64 {
        self.c[self.nodes - 1]
    }

    /// Volume-averaged concentration (mol m⁻³).
    pub fn average(&self) -> f64 {
        let num: f64 = self.c.iter().zip(&self.vol).map(|(c, v)| c * v).sum();
        let den: f64 = self.vol.iter().sum();
        num / den
    }

    /// Total moles in the particle (mol).
    pub fn moles(&self) -> f64 {
        self.c.iter().zip(&self.vol).map(|(c, v)| c * v).sum()
    }

    /// Advance with Crank–Nicolson under outward surface flux `j0`
    /// (mol m⁻² s⁻¹) held over the step.
    pub fn step(&mut self, dt: f64, j0: f64) {
        let n = self.nodes;
        let k = self.d / self.dr;
        let cap: Vec<f64> = (0..n).map(|i| 0.5 * dt / self.vol[i]).collect();
        let mut a = vec![0.0; n];
        let mut b = vec![0.0; n];
        let mut cc = vec![0.0; n];
        let mut rhs = vec![0.0; n];
        for i in 0..n {
            let ai = self.area[i];
            // No outer node beyond the last shell: its outer face carries the
            // imposed flux (added to the RHS below), not a diffusive coupling.
            let ao = if i + 1 < n { self.area[i + 1] } else { 0.0 };
            a[i] = -cap[i] * ai * k;
            b[i] = 1.0 + cap[i] * (ai + ao) * k;
            cc[i] = -cap[i] * ao * k;
            let explicit = ai * (if i == 0 { 0.0 } else { self.c[i - 1] - self.c[i] }) * k
                + ao * (if i + 1 < n { self.c[i + 1] - self.c[i] } else { 0.0 }) * k;
            rhs[i] = self.c[i] + cap[i] * explicit;
        }
        // Outer-face flux leaves the last node.
        let ar = self.area[n];
        rhs[n - 1] -= dt * ar * j0 / self.vol[n - 1];
        self.c = thomas(&mut a, &mut b, &mut cc, &mut rhs);
    }

    /// Advance `steps` steps of size `dt` under constant flux `j0`.
    pub fn run(&mut self, dt: f64, steps: usize, j0: f64) {
        for _ in 0..steps {
            self.step(dt, j0);
        }
    }
}

fn thomas(a: &mut [f64], b: &mut [f64], cc: &mut [f64], d: &mut [f64]) -> Vec<f64> {
    let n = b.len();
    for i in 1..n {
        let m = a[i] / b[i - 1];
        b[i] -= m * cc[i - 1];
        d[i] -= m * d[i - 1];
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d[n - 1] / b[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (d[i] - cc[i] * x[i + 1]) / b[i];
    }
    x
}

/// The roots `μ_n` of `tan μ = μ` (n = 1, 2, …), each in
/// `((n−½)π, (n+½)π)`, found by Newton from the midpoint.
pub fn neumann_roots(count: usize) -> Vec<f64> {
    let mut roots = Vec::with_capacity(count);
    for n in 1..=count {
        let mut mu = (n as f64 + 0.5) * std::f64::consts::PI;
        for _ in 0..60 {
            let f = mu.tan() - mu;
            let fp = 1.0 / mu.cos().powi(2) - 1.0;
            let step = f / fp;
            mu -= step;
            if step.abs() < 1e-14 {
                break;
            }
        }
        roots.push(mu);
    }
    roots
}

/// Exact surface concentration (mol m⁻³) under a constant outward flux `j0`
/// for time `t`, from uniform `c0`. Uses `terms` roots of `tan μ = μ`.
pub fn constant_flux_surface(t: f64, radius: f64, d: f64, j0: f64, c0: f64, terms: usize) -> f64 {
    let mut sum = 0.0;
    for mu in neumann_roots(terms) {
        let lam2 = d * mu * mu * t / (radius * radius);
        sum += (1.0 - (-lam2).exp()) / (mu * mu);
    }
    c0 - 3.0 * j0 * t / radius - (2.0 * radius * j0 / d) * sum
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn volume_average_is_exact_linear() {
        // Outward flux removes moles at the rate A J0: c̄ = c0 − 3 J0 t/R.
        let (r, d, j0, c0) = (1.0, 1.0, 0.01, 1.0);
        let mut p = SphereParticle::new(200, r, d, c0);
        let (dt, steps) = (2e-4, 250);
        p.run(dt, steps, j0);
        let t = dt * steps as f64;
        let exact = c0 - 3.0 * j0 * t / r;
        assert!((p.average() - exact).abs() < 1e-10, "avg {} vs {}", p.average(), exact);
    }

    #[test]
    fn surface_matches_constant_flux_series() {
        let (r, d, j0, c0) = (1.0, 1.0, 0.01, 1.0);
        let (dt, steps) = (2e-4, 1500);
        let mut p = SphereParticle::new(100, r, d, c0);
        p.run(dt, steps, j0);
        let t = dt * steps as f64;
        let exact = constant_flux_surface(t, r, d, j0, c0, 400);
        let rel = (p.surface() - exact).abs() / (c0 - exact).abs();
        assert!(rel < 1e-2, "surface {} vs series {} (rel {rel})", p.surface(), exact);
    }

    #[test]
    fn zero_flux_holds_uniform_state() {
        let mut p = SphereParticle::new(100, 1.0, 1.0, 0.7);
        p.run(1e-3, 100, 0.0);
        assert!((p.surface() - 0.7).abs() < 1e-12);
        assert!((p.average() - 0.7).abs() < 1e-12);
    }
}
