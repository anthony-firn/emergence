fn main() {
    let (a, b, c) = spectral_core::constants::agmon();
    println!("Agmon on the 2π-torus (coefficient norms): alpha = {a:.6}, beta = {b:.6}, C_A = {c:.6}  (4π = {:.6})", 4.0 * std::f64::consts::PI);
}
