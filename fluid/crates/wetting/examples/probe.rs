//! Print the Wenzel/Cassie fold and the finite-size convergence of the direct
//! free-energy minimiser on an explicit square-wave texture.
use wetting::lattice::{self, Branch, SquareWave};

fn main() {
    let wave = SquareWave::new(10.0, 4.0, 3.0);
    let costh = 110.0f64.to_radians().cos();
    let wenzel = lattice::wenzel_cos(&wave, costh).acos().to_degrees();
    let cassie = lattice::cassie_cos(&wave, costh).acos().to_degrees();
    println!("square posts r={:.4} phi_s={:.3}", wave.roughness(), wave.solid_fraction());
    println!("θ_Y = 110°, Wenzel = {wenzel:.2}°, Cassie = {cassie:.2}°");
    println!("  area    x_c/λ   Wenzel(measured)  err    Cassie(measured)  err");
    for &area in &[400.0f64, 1600.0, 6400.0, 25600.0, 102400.0] {
        let w = lattice::minimise(Branch::Wenzel, &wave, costh, area, 400.0, 8000);
        let c = lattice::minimise(Branch::Cassie, &wave, costh, area, 400.0, 8000);
        println!(
            "  {area:7.0} {:7.2}   {:8.2}°  {:6.2}   {:8.2}°  {:6.2}",
            lattice::size_ratio(&w, &wave),
            w.angle.to_degrees(),
            w.angle.to_degrees() - wenzel,
            c.angle.to_degrees(),
            c.angle.to_degrees() - cassie
        );
    }
}
