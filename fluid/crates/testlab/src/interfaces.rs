//! Interface and granular tests: the apparent contact angle on a textured
//! surface (wetting fold) and the effective friction / angle of repose of
//! glass beads (the granular μ(I) fold).
//!
//! Wetting: the coarse model is the Wenzel/Cassie–Baxter effective angle;
//! the validity indicator is the drop size over the texture period `a/λ`.
//! The fine model is a direct free-energy minimisation on an explicit
//! square-wave substrate (`wetting::lattice`), which has to converge to the
//! fold as the drop grows.
//!
//! Granular: the coarse model is `μ(I)` with the Jop et al. (2006) glass-bead
//! constants; the validity indicator is the inertial-number regime. A small
//! constant-pressure DEM shear cell supplies an independent `(I, μ)` point,
//! which is compared with the glass-bead closure.
//!
//! Sources: Young (1805); Wenzel (1936); Cassie & Baxter (1944); Bico,
//! Thiele & Quéré (2003); Jop, Forterre & Pouliquen, Nature 441, 727 (2006);
//! water γ = 72.8 mN/m and ρ = 998.21 kg/m³ at 20 °C (Vargaftik et al. 1983;
//! IAPWS-95/CRC).

use crate::{Level, Measurement};
use std::sync::OnceLock;
use wetting::lattice::{minimise, size_ratio, Branch, SquareWave};
use wetting::texture::{size_ratio as texture_size_ratio, SquarePosts};

/// The worked textured example of analysis/wetting.md: pitch λ = 10 µm,
/// post width a = 3 µm, height h = 8 µm, Young angle θ_Y = 110° (water on a
/// hydrophobic coating).
const PITCH: f64 = 10.0e-6;
const POST: f64 = 4.0e-6;
const HEIGHT: f64 = 3.0e-6;
const THETA_Y_DEG: f64 = 110.0;

/// A cached fine-model minimisation on the 2D square wave.
struct TextureFold {
    wenzel: (f64, f64), // (angle_deg, x_c/λ)
    cassie: (f64, f64),
}

fn texture_fold() -> &'static TextureFold {
    static CACHE: OnceLock<TextureFold> = OnceLock::new();
    CACHE.get_or_init(|| {
        let theta = THETA_Y_DEG.to_radians();
        // 2D cylindrical square wave (period, post width, height) in µm.
        let wave = SquareWave::new(PITCH * 1e6, POST * 1e6, HEIGHT * 1e6);
        let area = 102_400.0; // µm²
        let x_max = 200.0; // µm
        let w = minimise(Branch::Wenzel, &wave, theta.cos(), area, x_max, 400);
        let c = minimise(Branch::Cassie, &wave, theta.cos(), area, x_max, 400);
        TextureFold {
            wenzel: (w.angle.to_degrees(), size_ratio(&w, &wave)),
            cassie: (c.angle.to_degrees(), size_ratio(&c, &wave)),
        }
    })
}

/// The wetting-fold measurements.
pub fn wetting_contact_angle() -> Vec<Measurement> {
    let mut out = Vec::new();
    let theta = THETA_Y_DEG.to_radians();
    let posts = SquarePosts::new(PITCH, POST, HEIGHT);
    // The fine model is a 2D (cylindrical) square wave, so its reference is
    // the 2D Wenzel/Cassie fold, not the 3D post formula.
    let wave2 = SquareWave::new(PITCH * 1e6, POST * 1e6, HEIGHT * 1e6);
    let w_fold = wetting::texture::apparent_angle(wetting::lattice::wenzel_cos(&wave2, theta.cos()))
        .to_degrees();
    let c_fold = wetting::texture::apparent_angle(wetting::lattice::cassie_cos(&wave2, theta.cos()))
        .to_degrees();
    let r2 = wave2.roughness();
    let phis2 = wave2.solid_fraction();
    let transition2 = (phis2 - 1.0) / (r2 - phis2);
    let fold = texture_fold();

    // The fine model minimises the true capillary energy on the explicit
    // texture; the reference is the Wenzel homogenisation.
    out.push(Measurement::new(
        "water (textured)",
        "contact angle (Wenzel)",
        Level::Macroscopic,
        "θ* (fine model)",
        fold.wenzel.0,
        w_fold,
        "°",
        fold.wenzel.1 > 3.0,
        &format!(
            "2D drop free-energy minimisation on a square wave with r = {r2:.2}, φ_s = {phis2:.2}; \
             Wenzel fold θ* = {w_fold:.2}° at x_c/λ = {:.1}. The residual ({:+.2}°) is contact-line \
             pinning, the physical origin of hysteresis; the fold is certified once the drop \
             spans several texture periods.",
            fold.wenzel.1,
            fold.wenzel.0 - w_fold
        ),
    ));
    out.push(Measurement::new(
        "water (textured)",
        "contact angle (Cassie)",
        Level::Macroscopic,
        "θ* (fine model)",
        fold.cassie.0,
        c_fold,
        "°",
        fold.cassie.1 > 3.0,
        &format!(
            "same minimiser on the Cassie branch (liquid bridging the grooves). Cassie–Baxter \
             fold θ* = {c_fold:.2}°; the two branches are separated by an energy barrier, the \
             lower one stable (Wenzel here, since cos θ_Y = {:.3} < (φ_s−1)/(r−φ_s) = {transition2:.3}).",
            theta.cos(),
        ),
    ));
    // The cheap indicator itself.
    out.push(Measurement::new(
        "water (textured)",
        "contact angle (fold validity)",
        Level::Macroscopic,
        "a/λ",
        texture_size_ratio(0.5 * PITCH * fold.wenzel.1 * 2.0, PITCH),
        0.0,
        "",
        fold.wenzel.1 > 3.0,
        "the validity indicator: drop base radius over texture period. The homogenised angle is \
         certified only for a/λ ≫ 1; below that the texture must be resolved.",
    ));
    // Impalement check.
    let p_imp = posts.impalement_pressure(wetting::constants::WATER_GAMMA, theta);
    let r_c = posts.critical_drop_radius(theta);
    let sag = posts.sag_before_impalement(wetting::constants::WATER_GAMMA, theta);
    out.push(Measurement::new(
        "water (textured)",
        "contact angle (Cassie stability)",
        Level::Macroscopic,
        "sag/height",
        sag / HEIGHT,
        0.0,
        "",
        sag < HEIGHT,
        &format!(
            "bridging-meniscus sag before impalement = {:.3} µm vs the post height {:.1} µm \
             (sag/height = {:.3}), so sag is not the failure mechanism; the Laplace-pressure \
             criterion gives p_imp = {:.3e} Pa and a critical drop radius R_c = {:.1} µm. For a \
             1 mm drop the base pressure 2γ/R = {:.0} Pa ≪ p_imp, so Cassie is stable.",
            sag * 1e6,
            HEIGHT * 1e6,
            sag / HEIGHT,
            p_imp,
            r_c * 1e6,
            wetting::texture::drop_laplace_pressure(wetting::constants::WATER_GAMMA, 1e-3),
        ),
    ));
    out
}

/// A cached small DEM shear cell: `(I, mu)` at two shear rates.
struct DemMuI {
    points: Vec<(f64, f64, f64)>, // (gdot, I, mu)
}

fn dem_mu_i() -> &'static DemMuI {
    static CACHE: OnceLock<DemMuI> = OnceLock::new();
    CACHE.get_or_init(|| {
        use granular::contact::ContactParams;
        use granular::dem::ShearCell;
        let contact = ContactParams {
            young: 1000.0,
            poisson: 0.3,
            mu: 0.4,
            mu_rolling: 0.0,
            restitution: 0.5,
        };
        let mut points = Vec::new();
        for (k, &g) in [0.005_f64, 0.02].iter().enumerate() {
            let mut cell = ShearCell::new_cubic(108, 0.58, g, contact, 0.5, 1.0, 0.005, 1000 + k as u64);
            cell.set_pressure(0.05, 2.0);
            // 10000 steps: half equilibration, half production (the crate
            // averages over the second half), the configuration of the
            // granular-crate benchmark.
            let (p, tau, mu, i, _tg) = cell.run(10_000);
            let _ = (p, tau);
            points.push((g, i, mu));
        }
        DemMuI { points }
    })
}

/// The granular-fold measurements: μ(I) and the angle of repose.
pub fn glass_beads() -> Vec<Measurement> {
    use granular::rheology::{granular_regime, mu_i_glass_beads, GranularRegime, Jop2006};
    let mut out = Vec::new();
    let jop = Jop2006::glass_beads();
    let dem = dem_mu_i();
    for &(g, i, mu) in &dem.points {
        let mu_ref = mu_i_glass_beads(i);
        let regime = granular_regime(i);
        out.push(Measurement::new(
            "glass beads",
            "μ(I) shear cell",
            Level::Macroscopic,
            &format!("μ(γ̇ = {g:.3})"),
            mu,
            mu_ref,
            "",
            regime == GranularRegime::DenseFlow,
            &format!(
                "small constant-pressure DEM shear cell (108 spheres, P₀ = 0.05) at inertial \
                 number I = {i:.4}; reference is the Jop et al. (2006) glass-bead closure \
                 μ(I) = μ_s + (μ_2−μ_s)/(I₀/I+1) with μ_s = {:.3}, μ_2 = {:.3}, I₀ = {:.3}. The \
                 regime is {:?}. The DEM cross-over I₀ is larger than the glass-bead value \
                 (contact stiffness and interparticle friction), so the fit is not tuned.",
                jop.mu_s, jop.mu_2, jop.i0, regime
            ),
        ));
    }
    // The angle of repose from the low-inertial DEM point: tan θ = μ.
    let low = dem
        .points
        .iter()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .copied()
        .unwrap_or((0.0, 0.0, 0.33));
    let theta_repose = low.2.atan().to_degrees();
    let theta_jop = jop.mu_s.atan().to_degrees();
    out.push(Measurement::new(
        "glass beads",
        "angle of repose (DEM)",
        Level::Macroscopic,
        "θ_repose",
        theta_repose,
        theta_jop,
        "°",
        false,
        &format!(
            "tan θ = μ at the lowest resolved inertial number in the small cell (I = {:.3}); the \
             reference is the Jop glass-bead quasi-static friction angle atan({:.3}) = {:.1}°. \
             The DEM is in the dense-flow regime, not quasi-static, and its cross-over I₀ \
             differs, so the angle-of-repose fold is reported invalid; a quantitative repose \
             angle needs a much larger heap and slow pouring.",
            low.1, jop.mu_s, theta_jop
        ),
    ));
    // The regime indicator.
    out.push(Measurement::new(
        "glass beads",
        "μ(I) shear cell",
        Level::Macroscopic,
        "I at γ̇ = 0.3, P = 0.05",
        granular::rheology::inertial_number(0.3, 1.0, 0.05, 1.0),
        0.1,
        "",
        granular_regime(0.1) == GranularRegime::DenseFlow,
        "the cheap validity indicator I = γ̇ d/√(P/ρ_p). The μ(I) closure is valid for \
         I ∈ (10⁻³, 0.1); above it the collisional/kinetic-theory unfold is required.",
    ));
    out
}
