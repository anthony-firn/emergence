//! A model human eye: the Escudero-Sanz & Navarro (JOSA A 16, 1881 (1999))
//! wide-angle schematic eye (relaxed), with ocular-media dispersion.
//!
//! Two dispersion variants:
//!
//! * [`EyeDispersion::Navarro`]: the published media indices at 458, 543,
//!   589.3 and 632.8 nm (Cauchy-interpolated). In this four-surface eye it
//!   gives about 1.3× the measured human longitudinal chromatic aberration
//!   (LCA) over 458–633 nm.
//! * [`EyeDispersion::Physiological`]: the same media with every
//!   dispersion `n(λ) − n(589 nm)` scaled by one common factor, chosen by
//!   least squares so that the eye's LCA matches the physiological LCA of
//!   Atchison & Smith (JOSA A 22, 29 (2005)) over 450–650 nm — the
//!   procedure those authors used to make schematic eyes physiological.
//!   The shape of the curve over 450–700 nm and the independent
//!   Thibos et al. (Appl. Opt. 31, 3594 (1992)) chromatic-eye formula are
//!   then checks, not inputs.
//!
//! The retina is represented by its tangent plane (adequate within a few
//! degrees of the fovea), and the 3 mm pupil sits at the front of the
//! crystalline lens.

use std::sync::OnceLock;

use crate::materials::Material;
use crate::raytrace::{Surface, System};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EyeDispersion {
    Navarro,
    Physiological,
}

/// Physiological LCA, Atchison & Smith (2005): chromatic difference of
/// refraction (D), λ in metres (fit in µm).
pub fn lca_atchison_smith(lambda: f64) -> f64 {
    let l = lambda * 1e6;
    1.60911 - 0.670941 * l.powi(-2) + 0.0555334 * l.powi(-4) - 0.00559998 * l.powi(-6)
}

/// Thibos et al. (1992) chromatic eye: refractive error (D), λ in metres
/// (fit in nm).
pub fn lca_thibos(lambda: f64) -> f64 {
    1.68524 - 633.46 / (lambda * 1e9 - 214.102)
}

const REF: f64 = 589.3e-9;

fn scaled(m: Material, s: f64) -> Material {
    match m {
        Material::Cauchy { a, b, c } => {
            let l2 = (REF * 1e6).powi(2);
            let disp_ref = b / l2 + c / (l2 * l2);
            Material::Cauchy { a: a + disp_ref - s * disp_ref, b: s * b, c: s * c }
        }
        other => other,
    }
}

fn build(pupil_diameter: f64, s: f64) -> System {
    let mm = 1e-3;
    System {
        object_medium: Material::Air,
        surfaces: vec![
            Surface { radius: 7.72 * mm, conic: -0.26, semi_diameter: 5.0 * mm, thickness: 0.55 * mm, after: scaled(Material::cornea(), s) },
            Surface { radius: 6.50 * mm, conic: 0.0, semi_diameter: 5.0 * mm, thickness: 3.05 * mm, after: scaled(Material::aqueous(), s) },
            Surface { radius: 10.2 * mm, conic: -3.1316, semi_diameter: 4.5 * mm, thickness: 4.0 * mm, after: scaled(Material::crystalline_lens(), s) },
            Surface { radius: -6.0 * mm, conic: -1.0, semi_diameter: 4.5 * mm, thickness: 16.3203 * mm, after: scaled(Material::vitreous(), s) },
        ],
        stop_z: 3.60 * mm,
        stop_radius: 0.5 * pupil_diameter,
        image_z: 23.9203 * mm,
    }
}

/// The eye's chromatic difference of refraction (D) relative to 589.3 nm:
/// minus the change of power.
pub fn model_lca(eye: &System, lambda: f64) -> f64 {
    -(1.0 / eye.paraxial(lambda).efl - 1.0 / eye.paraxial(REF).efl)
}

/// Common dispersion scale making the eye physiological (cached).
pub fn physiological_scale() -> f64 {
    static S: OnceLock<f64> = OnceLock::new();
    *S.get_or_init(|| {
        let eye = build(3e-3, 1.0);
        let (mut num, mut den) = (0.0, 0.0);
        for nm in (450..=650).step_by(10) {
            let l = nm as f64 * 1e-9;
            let m = model_lca(&eye, l);
            let h = lca_atchison_smith(l) - lca_atchison_smith(REF);
            num += m * h;
            den += m * m;
        }
        num / den
    })
}

/// The model eye with a pupil of `pupil_diameter` (m).
pub fn eye(pupil_diameter: f64, dispersion: EyeDispersion) -> System {
    let s = match dispersion {
        EyeDispersion::Navarro => 1.0,
        EyeDispersion::Physiological => physiological_scale(),
    };
    build(pupil_diameter, s)
}
