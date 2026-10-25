//! 1D transient heat conduction through layered media, with a volumetric
//! optical heat source.
//!
//! This is the rung between the optical heat source (`photonics::multilayer`
//! returns an absorbed power density `p(z)` in W m^-3) and the flow solver: it
//! answers "how hot does the stack get, and when". The PDE is
//!
//! \[
//!   \rho c_p \frac{\partial T}{\partial t}
//!     = \frac{\partial}{\partial z}\!\left(k \frac{\partial T}{\partial z}\right)
//!     + p(z,t),
//! \]
//!
//! solved by a conservative finite-volume method on a grid aligned with the
//! layer interfaces. Cell temperature is the cell average; the face flux is
//! `G (T_L - T_R)` with the series conductance
//! `G = 1 / (d_L/(2 k_L) + d_R/(2 k_R))`. This makes the flux continuous across
//! a material interface by construction, and reduces to `k/dz` inside a layer.
//! Time stepping is backward Euler (unconditionally stable for diffusion), so
//! the only stability limit is the physics, not the step.
//!
//! # Emergent boundary
//!
//! This conduction model is the coarse (folded) rung. Its validity indicator is
//! the Rayleigh number of the heated layer ([`crate::indicators`]): below the
//! critical value the fluid is motionless and conduction is exact; above it,
//! the model must unfold to the flow solver.

use std::f64::consts::PI;

/// Thermal properties of a homogeneous material.
#[derive(Clone, Copy, Debug)]
pub struct Material {
    pub name: &'static str,
    /// Mass density [kg m^-3].
    pub rho: f64,
    /// Specific heat capacity [J kg^-1 K^-1].
    pub cp: f64,
    /// Thermal conductivity [W m^-1 K^-1].
    pub k: f64,
    pub source: &'static str,
}

impl Material {
    /// Volumetric heat capacity [J m^-3 K^-1].
    pub fn rho_cp(&self) -> f64 {
        self.rho * self.cp
    }
    /// Thermal diffusivity `alpha = k / (rho cp)` [m^2 s^-1].
    pub fn alpha(&self) -> f64 {
        self.k / self.rho_cp()
    }
    /// Thermal effusivity `sqrt(rho cp k)` [J m^-2 K^-1 s^-1/2].
    pub fn effusivity(&self) -> f64 {
        (self.rho_cp() * self.k).sqrt()
    }
}

/// Water at 20 C, 1 atm. `rho` from the CRC Handbook; `cp`, `k` from
/// Incropera & DeWitt, *Fundamentals of Heat and Mass Transfer*, Table A.6
/// (saturated liquid at 300 K).
pub const WATER: Material = Material {
    name: "water",
    rho: 997.0,
    cp: 4179.0,
    k: 0.613,
    source: "CRC Handbook (rho, 20 C); Incropera & DeWitt Table A.6, 300 K (cp=4179, k=0.613)",
};

/// Fused silica (amorphous SiO2), representative room-temperature values from
/// Incropera & DeWitt, Table A.3.
pub const FUSED_SILICA: Material = Material {
    name: "fused silica",
    rho: 2200.0,
    cp: 745.0,
    k: 1.38,
    source: "Incropera & DeWitt, Fundamentals of Heat and Mass Transfer, Table A.3 (fused silica)",
};

/// Borosilicate glass (Pyrex 7740), representative room-temperature values
/// from Incropera & DeWitt, Table A.3.
pub const PYREX: Material = Material {
    name: "borosilicate glass",
    rho: 2225.0,
    cp: 835.0,
    k: 1.4,
    source: "Incropera & DeWitt, Table A.3 (Pyrex 7740: rho=2225, cp=835, k=1.4)",
};

/// Gold. Density and `cp` from the CRC Handbook; `k` from Incropera &
/// DeWitt, Table A.1 (pure gold at 300 K). Thin-film values are lower.
pub const GOLD: Material = Material {
    name: "gold",
    rho: 19300.0,
    cp: 129.0,
    k: 317.0,
    source: "CRC Handbook (rho, cp); Incropera & DeWitt Table A.1 (k=317 W/m/K, 300 K)",
};

/// Magnesium fluoride MgF2, a common single-layer antireflection coating.
/// Bulk thermal values are approximate thin-film estimates; the optical
/// coating is what matters here, and the thermal mass is negligible.
pub const MGF2: Material = Material {
    name: "MgF2 (coating)",
    rho: 3180.0,
    cp: 1010.0,
    k: 11.0,
    source: "estimate: bulk MgF2 handbook values (rho=3180, cp=1010, k=11)",
};

/// Boundary condition at one end of the 1D domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Boundary {
    /// Zero heat flux (perfect insulation).
    Adiabatic,
    /// Fixed temperature [K].
    Fixed(f64),
    /// Inward surface heat flux [W m^-2]. Only meaningful at the front face,
    /// where it is the absorbed optical irradiance.
    Flux(f64),
}

/// One homogeneous thermal layer.
#[derive(Clone, Copy, Debug)]
pub struct ThermalLayer {
    pub material: Material,
    pub thickness: f64,
}

/// A layered 1D conduction domain on a cell-centered, layer-aligned grid.
///
/// Cell `i` has centre `z[i]` and width `dz[i]`; `g[i]` is the conductance
/// [W m^-2 K^-1] between cell `i-1` and cell `i` (`g[0]` the front boundary,
/// `g[n]` the back boundary), so that the heat flux from `i-1` into `i` is
/// `g[i] (T[i-1] - T[i])`.
pub struct ThermalModel {
    pub layers: Vec<ThermalLayer>,
    pub z: Vec<f64>,
    pub dz: Vec<f64>,
    pub mat: Vec<Material>,
    pub g: Vec<f64>,
    pub front: Boundary,
    pub back: Boundary,
}

impl ThermalModel {
    /// Build a grid with `cells_per_layer` uniform cells in every layer.
    pub fn new(
        layers: Vec<ThermalLayer>,
        cells_per_layer: usize,
        front: Boundary,
        back: Boundary,
    ) -> Self {
        assert!(!layers.is_empty() && cells_per_layer > 0);
        assert!(
            !matches!(back, Boundary::Flux(_)),
            "Flux is a front-face condition"
        );
        let mut z = Vec::new();
        let mut dz = Vec::new();
        let mut mat = Vec::new();
        let mut z0 = 0.0;
        for l in &layers {
            assert!(l.thickness > 0.0 && l.material.k > 0.0);
            let h = l.thickness / cells_per_layer as f64;
            for i in 0..cells_per_layer {
                z.push(z0 + h * (i as f64 + 0.5));
                dz.push(h);
                mat.push(l.material);
            }
            z0 += l.thickness;
        }
        let mut m = Self {
            layers,
            z,
            dz,
            mat,
            g: Vec::new(),
            front,
            back,
        };
        m.rebuild_conductances();
        m
    }

    /// Number of cells.
    pub fn n(&self) -> usize {
        self.z.len()
    }

    /// Total stack thickness [m].
    pub fn length(&self) -> f64 {
        self.dz.iter().sum()
    }

    fn rebuild_conductances(&mut self) {
        let n = self.n();
        let mut g = vec![0.0; n + 1];
        // Front boundary.
        g[0] = match self.front {
            Boundary::Fixed(_) => 2.0 * self.mat[0].k / self.dz[0],
            _ => 0.0,
        };
        // Interior faces: series resistance to the common interface.
        for i in 1..n {
            let a = &self.mat[i - 1];
            let b = &self.mat[i];
            g[i] = 1.0 / (0.5 * self.dz[i - 1] / a.k + 0.5 * self.dz[i] / b.k);
        }
        // Back boundary.
        g[n] = match self.back {
            Boundary::Fixed(_) => 2.0 * self.mat[n - 1].k / self.dz[n - 1],
            _ => 0.0,
        };
        self.g = g;
    }

    /// One backward-Euler step. `source[i]` is the volumetric absorbed power
    /// density [W m^-3] in cell `i`.
    pub fn step(&self, t: &mut [f64], source: &[f64], dt: f64) {
        let n = self.n();
        assert_eq!(t.len(), n);
        assert_eq!(source.len(), n);
        assert!(dt > 0.0);
        let mut sub = vec![0.0; n];
        let mut diag = vec![0.0; n];
        let mut sup = vec![0.0; n];
        let mut rhs = vec![0.0; n];
        for i in 0..n {
            let cap = self.mat[i].rho_cp() * self.dz[i];
            diag[i] = cap / dt + self.g[i] + self.g[i + 1];
            rhs[i] = cap / dt * t[i] + source[i] * self.dz[i];
            if i > 0 {
                sub[i] = -self.g[i];
            }
            if i + 1 < n {
                sup[i] = -self.g[i + 1];
            }
        }
        match self.front {
            Boundary::Fixed(tb) => rhs[0] += self.g[0] * tb,
            Boundary::Flux(q) => rhs[0] += q,
            Boundary::Adiabatic => {}
        }
        if let Boundary::Fixed(tb) = self.back {
            rhs[n - 1] += self.g[n] * tb;
        }
        let sol = thomas(&sub, &diag, &sup, &rhs);
        t.copy_from_slice(&sol);
    }

    /// Advance for `steps` steps of `dt`, reusing one source vector.
    pub fn run(&self, t: &mut [f64], source: &[f64], dt: f64, steps: usize) {
        for _ in 0..steps {
            self.step(t, source, dt);
        }
    }

    /// Steady-state temperature for a time-independent source.
    pub fn steady_state(&self, source: &[f64]) -> Vec<f64> {
        let n = self.n();
        assert_eq!(source.len(), n);
        let mut sub = vec![0.0; n];
        let mut diag = vec![0.0; n];
        let mut sup = vec![0.0; n];
        let mut rhs = vec![0.0; n];
        for i in 0..n {
            diag[i] = self.g[i] + self.g[i + 1];
            rhs[i] = source[i] * self.dz[i];
            if i > 0 {
                sub[i] = -self.g[i];
            }
            if i + 1 < n {
                sup[i] = -self.g[i + 1];
            }
        }
        match self.front {
            Boundary::Fixed(tb) => rhs[0] += self.g[0] * tb,
            Boundary::Flux(q) => rhs[0] += q,
            Boundary::Adiabatic => {}
        }
        if let Boundary::Fixed(tb) = self.back {
            rhs[n - 1] += self.g[n] * tb;
        }
        thomas(&sub, &diag, &sup, &rhs)
    }

    /// Total stored sensible energy per unit area [J m^-2] relative to 0 K.
    pub fn stored_energy(&self, t: &[f64]) -> f64 {
        t.iter()
            .zip(&self.mat)
            .zip(&self.dz)
            .map(|((&ti, m), &d)| m.rho_cp() * d * ti)
            .sum()
    }

    /// Energy per unit area [J m^-2] deposited by a constant source up to `time`.
    pub fn source_energy(&self, source: &[f64], time: f64) -> f64 {
        source
            .iter()
            .zip(&self.dz)
            .map(|(&p, &d)| p * d)
            .sum::<f64>()
            * time
    }

    /// Heat flux [W m^-2] entering the domain at the front face for the
    /// current state. Positive heats the stack.
    pub fn front_flux(&self, t: &[f64]) -> f64 {
        match self.front {
            Boundary::Fixed(tb) => self.g[0] * (tb - t[0]),
            Boundary::Flux(q) => q,
            Boundary::Adiabatic => 0.0,
        }
    }

    /// Heat flux [W m^-2] entering the domain at the back face.
    pub fn back_flux(&self, t: &[f64]) -> f64 {
        let n = self.n();
        match self.back {
            Boundary::Fixed(tb) => self.g[n] * (tb - t[n - 1]),
            Boundary::Adiabatic => 0.0,
            Boundary::Flux(_) => unreachable!("Flux is a front-face condition"),
        }
    }
}

/// Solve a tridiagonal system `sub[i] x[i-1] + diag[i] x[i] + sup[i] x[i+1] = rhs[i]`.
pub fn thomas(sub: &[f64], diag: &[f64], sup: &[f64], rhs: &[f64]) -> Vec<f64> {
    let n = diag.len();
    let mut c = vec![0.0; n];
    let mut d = vec![0.0; n];
    c[0] = sup[0] / diag[0];
    d[0] = rhs[0] / diag[0];
    for i in 1..n {
        let m = diag[i] - sub[i] * c[i - 1];
        c[i] = if i + 1 < n { sup[i] / m } else { 0.0 };
        d[i] = (rhs[i] - sub[i] * d[i - 1]) / m;
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = d[i] - c[i] * x[i + 1];
    }
    x
}

/// Beer-Lambert volumetric source in a semi-infinite absorber,
/// `p(z) = mu_a I0 exp(-mu_a z)` [W m^-3]. `mu_a` is the absorption
/// coefficient [m^-1], `I0` the incident irradiance [W m^-2].
pub fn beer_lambert_source(mu_a: f64, irradiance: f64, z: f64) -> f64 {
    mu_a * irradiance * (-mu_a * z).exp()
}

/// Cell-average of the Beer-Lambert source over `[z0, z0+dz]` [W m^-3], exact:
/// `I0 (exp(-mu z0) - exp(-mu (z0+dz))) / dz`.
pub fn beer_lambert_cell_average(mu_a: f64, irradiance: f64, z0: f64, dz: f64) -> f64 {
    irradiance * ((-mu_a * z0).exp() - (-mu_a * (z0 + dz)).exp()) / dz
}

/// Sample the photonic absorbed power density `p(z)` at the thermal cell
/// centres. Feed the result to [`ThermalModel::step`]. Sampling at centres is
/// second-order accurate for the cell average on a smooth profile; where the
/// field varies on the cell scale, refine the grid.
pub fn source_from_optics(
    solution: &photonics::multilayer::Solution,
    irradiance: f64,
    z: &[f64],
) -> Vec<f64> {
    z.iter()
        .map(|&zi| solution.power_density(zi, irradiance))
        .collect()
}

/// Exact semi-infinite solution for a constant surface heat flux `q` [W m^-2]
/// switched on at `t=0`, insulated surface, uniform material:
///
/// \[
///   T(z,t) = \frac{2 q}{k}\sqrt{\alpha t}\,
///            \left[\frac{1}{\sqrt\pi} e^{-z^2/(4\alpha t)}
///              - \frac{z}{2\sqrt{\alpha t}}\,\mathrm{erfc}\!\left(\frac{z}{2\sqrt{\alpha t}}\right)\right].
/// \]
///
/// At `z=0` this is `2 q sqrt(t/(pi rho cp k))`.
pub fn semi_infinite_surface_flux(mat: &Material, z: f64, t: f64, q: f64) -> f64 {
    if t <= 0.0 {
        return 0.0;
    }
    let a = mat.alpha();
    let s = (a * t).sqrt();
    let x = z / (2.0 * s);
    (2.0 * q / mat.k) * s * ((-x * x).exp() / PI.sqrt() - x * erfc(x))
}

/// Complementary error function. Uses a Taylor series below `x = 1` and the
/// continued fraction for the upper incomplete gamma `Gamma(1/2, x^2)` above
/// it, so it is accurate to about `1e-14` over the range the heat benchmarks
/// need. This is a numerical helper, not a physical model.
pub fn erfc(x: f64) -> f64 {
    erfc_high(x)
}

/// High-accuracy complementary error function (see [`erfc`]).
pub fn erfc_high(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc_high(-x);
    }
    if x < 1.0 {
        return 1.0 - erf_taylor(x);
    }
    // erfc(x) = exp(-x^2) x/sqrt(pi) / D, where
    // D = x^2 + 1/2 - (1/2)/(x^2 + 5/2 - 3/(x^2 + 9/2 - ...)).
    // Evaluate with the modified Lentz algorithm.
    let x2 = x * x;
    let tiny = 1e-300;
    let mut f = x2 + 0.5;
    let mut c = f;
    let mut d = 0.0;
    for n in 1..400 {
        let nf = n as f64;
        let a = -nf * (nf - 0.5);
        let b = x2 + 2.0 * nf + 0.5;
        d = b + a * d;
        if d.abs() < tiny {
            d = tiny;
        }
        c = b + a / c;
        if c.abs() < tiny {
            c = tiny;
        }
        d = 1.0 / d;
        let delta = c * d;
        f *= delta;
        if (delta - 1.0).abs() < 1e-17 {
            break;
        }
    }
    (-x2).exp() * x / PI.sqrt() / f
}

fn erf_taylor(x: f64) -> f64 {
    // erf(x) = 2/sqrt(pi) sum_{n>=0} (-1)^n x^{2n+1}/(n!(2n+1))
    let x2 = x * x;
    let mut term = x;
    let mut sum = x;
    let mut n = 1.0f64;
    loop {
        term *= -x2 / n;
        let add = term / (2.0 * n + 1.0);
        sum += add;
        if add.abs() < 1e-18 * sum.abs().max(1e-300) {
            break;
        }
        n += 1.0;
        if n > 200.0 {
            break;
        }
    }
    2.0 / PI.sqrt() * sum
}
