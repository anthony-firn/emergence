//! Walls, roofs and floors: layered assemblies.
//!
//! * **Steady**: U = 1/(R_si + Σ d/λ + R_se), surface resistances from
//!   ISO 6946:2017 Table 7 (R_si 0.13, 0.10, 0.17 m²K/W for horizontal,
//!   upward and downward heat flow; R_se 0.04).
//! * **Periodic** (ISO 13786:2017): each layer's exact heat-transfer matrix
//!   for a sinusoidal temperature of period P,
//!   `[θ₁; q₁] = [[cosh γd, sinh γd/(λγ)], [λγ sinh γd, cosh γd]] [θ₂; q₂]`,
//!   `γ = √(iω/a)`, giving the decrement factor (how much of an outdoor
//!   swing reaches indoors) and the time lag.
//! * **Transient**: a finite-volume wall (several nodes per layer,
//!   BDF2, which is L-stable: hourly steps damp the stiff modes of thin
//!   layers, where Crank–Nicolson leaves them ringing) for hour-by-hour
//!   load simulation. It is checked
//!   against the exact periodic solution in `tests/envelope.rs`.

use std::f64::consts::PI;

#[derive(Clone, Copy, Debug)]
pub struct Layer {
    pub thickness: f64,
    /// Thermal conductivity (W/(m·K)).
    pub k: f64,
    pub rho: f64,
    pub cp: f64,
}

/// Direction of heat flow through the element, for the surface resistances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Flow {
    Horizontal,
    Upward,
    Downward,
}

/// Internal and external surface resistances (ISO 6946:2017 Table 7).
pub fn surface_resistances(flow: Flow) -> (f64, f64) {
    let rsi = match flow {
        Flow::Horizontal => 0.13,
        Flow::Upward => 0.10,
        Flow::Downward => 0.17,
    };
    (rsi, 0.04)
}

#[derive(Clone, Debug)]
pub struct Assembly {
    /// Layers from inside to outside.
    pub layers: Vec<Layer>,
    pub r_si: f64,
    pub r_se: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct C(f64, f64);

impl C {
    fn add(self, o: C) -> C {
        C(self.0 + o.0, self.1 + o.1)
    }
    fn mul(self, o: C) -> C {
        C(self.0 * o.0 - self.1 * o.1, self.0 * o.1 + self.1 * o.0)
    }
    fn scale(self, s: f64) -> C {
        C(self.0 * s, self.1 * s)
    }
    fn div(self, o: C) -> C {
        let d = o.0 * o.0 + o.1 * o.1;
        C((self.0 * o.0 + self.1 * o.1) / d, (self.1 * o.0 - self.0 * o.1) / d)
    }
    fn abs(self) -> f64 {
        self.0.hypot(self.1)
    }
    fn arg(self) -> f64 {
        self.1.atan2(self.0)
    }
    fn cosh(self) -> C {
        C(self.0.cosh() * self.1.cos(), self.0.sinh() * self.1.sin())
    }
    fn sinh(self) -> C {
        C(self.0.sinh() * self.1.cos(), self.0.cosh() * self.1.sin())
    }
}

type M = [[C; 2]; 2];

fn mm(a: M, b: M) -> M {
    let mut r = [[C(0.0, 0.0); 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            r[i][j] = a[i][0].mul(b[0][j]).add(a[i][1].mul(b[1][j]));
        }
    }
    r
}

/// Periodic response of an assembly.
#[derive(Clone, Copy, Debug)]
pub struct Periodic {
    /// |Y₁₂|: periodic thermal transmittance (W/(m²K)).
    pub transmittance: f64,
    /// |Y₁₂| / U.
    pub decrement: f64,
    /// Delay of the indoor heat-flow peak after the outdoor peak (s).
    pub lag: f64,
    /// |Y₁₁|: internal areal admittance (W/(m²K)), how strongly the inner
    /// face stores heat.
    pub admittance: f64,
}

impl Assembly {
    pub fn new(layers: Vec<Layer>, flow: Flow) -> Self {
        let (r_si, r_se) = surface_resistances(flow);
        Assembly { layers, r_si, r_se }
    }

    pub fn r_total(&self) -> f64 {
        self.r_si + self.layers.iter().map(|l| l.thickness / l.k).sum::<f64>() + self.r_se
    }

    pub fn u_value(&self) -> f64 {
        1.0 / self.r_total()
    }

    /// Heat capacity per area (J/(m²K)).
    pub fn heat_capacity(&self) -> f64 {
        self.layers.iter().map(|l| l.thickness * l.rho * l.cp).sum()
    }

    /// ISO 13786 periodic response at period `period` (s), indoor to
    /// outdoor, surface resistances included.
    pub fn periodic(&self, period: f64) -> Periodic {
        let w = 2.0 * PI / period;
        let res = |r: f64| -> M { [[C(1.0, 0.0), C(-r, 0.0)], [C(0.0, 0.0), C(1.0, 0.0)]] };
        // ISO 13786 §6: Z = Z_se · Z_N ⋯ Z_1 · Z_si, layer 1 inside.
        let mut z = res(self.r_si);
        for l in &self.layers {
            let a = l.k / (l.rho * l.cp);
            // γ = √(iω/a) = (1 + i)/δ, δ = √(2a/ω)
            let g = (w / (2.0 * a)).sqrt();
            let gd = C(g * l.thickness, g * l.thickness);
            let kg = C(l.k * g, l.k * g);
            let m: M = [[gd.cosh(), gd.sinh().div(kg).scale(-1.0)], [kg.mul(gd.sinh()).scale(-1.0), gd.cosh()]];
            z = mm(m, z);
        }
        z = mm(res(self.r_se), z);
        // Y₁₂ = −1/Z₁₂ (inside heat flow per outside temperature swing),
        // Y₁₁ = −Z₁₁/Z₁₂ (inside admittance).
        let y12 = C(-1.0, 0.0).div(z[0][1]);
        let y11 = z[0][0].scale(-1.0).div(z[0][1]);
        let u = self.u_value();
        let mut lag = -y12.arg() / w;
        while lag < 0.0 {
            lag += period;
        }
        Periodic { transmittance: y12.abs(), decrement: y12.abs() / u, lag, admittance: y11.abs() }
    }
}

/// A finite-volume wall for time stepping: `n` nodes per layer, node 0 at
/// the inside face, interfaces exact, BDF2 in time (backward Euler for the
/// first step).
#[derive(Clone, Debug)]
pub struct Wall {
    /// Node temperatures (°C), inside to outside.
    pub t: Vec<f64>,
    prev: Option<Vec<f64>>,
    cap: Vec<f64>,
    /// Conductances between consecutive nodes (W/(m²K)).
    g: Vec<f64>,
    /// Conductance from the inside face node to the inside surface, and the
    /// outside face node to the outside surface (half cells).
    g_in: f64,
    g_out: f64,
    pub r_si: f64,
    pub r_se: f64,
}

impl Wall {
    pub fn new(a: &Assembly, nodes_per_layer: usize, t0: f64) -> Wall {
        let mut cap = Vec::new();
        let mut half_r = Vec::new(); // (left half resistance, right half resistance) of each cell
        for l in &a.layers {
            let dx = l.thickness / nodes_per_layer as f64;
            for _ in 0..nodes_per_layer {
                cap.push(l.rho * l.cp * dx);
                half_r.push(0.5 * dx / l.k);
            }
        }
        let g = (0..cap.len() - 1).map(|i| 1.0 / (half_r[i] + half_r[i + 1])).collect();
        let n = cap.len();
        Wall { t: vec![t0; n], prev: None, cap, g, g_in: 1.0 / half_r[0], g_out: 1.0 / half_r[n - 1], r_si: a.r_si, r_se: a.r_se }
    }

    /// Advance `dt` (constant from step to step) to inside air `t_in`, the
    /// outside surface driven by `t_out` (outdoor or sol-air temperature)
    /// and a flux `q_in` (W/m², e.g. radiant gains) absorbed at the inside
    /// surface, all at the end of the step. Returns the heat flux into the
    /// room then.
    pub fn step(&mut self, dt: f64, t_in: f64, t_out: f64, q_in: f64) -> f64 {
        let n = self.t.len();
        let gi = 1.0 / (self.r_si + 1.0 / self.g_in);
        let go = 1.0 / (self.r_se + 1.0 / self.g_out);
        // c·C T¹ − A T¹ = rhs: BDF2 c = 3/(2dt), rhs = C(4T⁰ − T⁻¹)/(2dt);
        // backward Euler c = 1/dt, rhs = C T⁰/dt.
        let (c, rhs): (f64, Vec<f64>) = match &self.prev {
            Some(p) => (1.5 / dt, (0..n).map(|i| self.cap[i] * (4.0 * self.t[i] - p[i]) / (2.0 * dt)).collect()),
            None => (1.0 / dt, (0..n).map(|i| self.cap[i] * self.t[i] / dt).collect()),
        };
        let mut lo = vec![0.0; n];
        let mut di: Vec<f64> = (0..n).map(|i| c * self.cap[i]).collect();
        let mut up = vec![0.0; n];
        let mut r = rhs;
        for i in 0..n - 1 {
            di[i] += self.g[i];
            di[i + 1] += self.g[i];
            up[i] -= self.g[i];
            lo[i + 1] -= self.g[i];
        }
        // The absorbed inside flux reaches the node through the half cell
        // in series with R_si: its share is R_si·gi.
        di[0] += gi;
        r[0] += gi * t_in + q_in * self.r_si * gi;
        di[n - 1] += go;
        r[n - 1] += go * t_out;
        for i in 1..n {
            let m = lo[i] / di[i - 1];
            di[i] -= m * up[i - 1];
            r[i] -= m * r[i - 1];
        }
        let mut t = vec![0.0; n];
        t[n - 1] = r[n - 1] / di[n - 1];
        for i in (0..n - 1).rev() {
            t[i] = (r[i] - up[i] * t[i + 1]) / di[i];
        }
        self.prev = Some(std::mem::replace(&mut self.t, t));
        self.flux_in(t_in, q_in)
    }

    /// Heat flux from the wall into the room air (W/m²): convective and
    /// radiative exchange through R_si, at the current state.
    pub fn flux_in(&self, t_in: f64, q_abs: f64) -> f64 {
        // Surface temperature from the series of the half cell and R_si,
        // with the absorbed flux arriving at the surface.
        let ts = (self.t[0] * self.g_in + t_in / self.r_si + q_abs) / (self.g_in + 1.0 / self.r_si);
        (ts - t_in) / self.r_si
    }

    /// Inside surface temperature (°C).
    pub fn surface_in(&self, t_in: f64, q_abs: f64) -> f64 {
        (self.t[0] * self.g_in + t_in / self.r_si + q_abs) / (self.g_in + 1.0 / self.r_si)
    }
}
