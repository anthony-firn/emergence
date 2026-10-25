//! Which models hold for some real flows (the model-ladder selector).
use spectral_core::regime::{FlowNumbers, Material};

fn main() {
    let cases = [
        ("mercury homopolar motor: 5 cm dish, 5 cm/s, 0.3 T", Material::MERCURY, 0.05, 0.05, 0.3, 101_325.0),
        ("tap-water stream: 1 cm, 1 m/s", Material::WATER, 0.01, 1.0, 0.0, 101_325.0),
        ("ship propeller tip: 0.5 m, 30 m/s, 1 m deep", Material::WATER, 0.5, 30.0, 0.0, 111_000.0),
        ("nitrogen jet: 1 cm, 500 m/s", Material::NITROGEN_GAS, 0.01, 500.0, 0.0, 101_325.0),
        ("liquid nitrogen pour: 1 cm, 0.5 m/s", Material::LIQUID_NITROGEN, 0.01, 0.5, 0.0, 101_325.0),
    ];
    for (name, mat, l, u, b, p) in cases {
        let f = FlowNumbers::new(&mat, l, u, b, p);
        println!("{name}");
        println!(
            "  Re {:.2e}  Ma {:.2e}  We {:.2e}  σ_c {:.2e}  Rm {:.2e}  Ha {:.2e}  τ_e·U/L {:.2e}",
            f.reynolds, f.mach, f.weber, f.cavitation, f.magnetic_reynolds, f.hartmann, f.charge_relaxation
        );
        let (ok, unfold): (Vec<_>, Vec<_>) = f.regimes().into_iter().partition(|r| r.1);
        println!("  holds:  {}", ok.iter().map(|r| r.0).collect::<Vec<_>>().join("; "));
        println!("  breaks: {}", unfold.iter().map(|r| r.0).collect::<Vec<_>>().join("; "));
    }
}
