//! Molecular substances (`substances.rs`) against data the fits never saw:
//! standard state, liquid vapour pressure, sublimation pressure, triple
//! point and (water) the melting curve of ice Ih.
//!
//!   substances
use compressible::elements::{Phase, P_ATM};
use compressible::substances::all;

fn main() {
    for s in all() {
        let m = s.model();
        let st = match m.stable(P_ATM, 298.15) {
            Phase::Solid => "solid",
            Phase::Liquid => "liquid",
            Phase::Vapor => "gas",
        };
        println!("{} ({}): grade {:?}, stable at 298 K / 1 atm: {st}", s.name, s.formula, m.grade);
        let rel = |t_model: Option<f64>, t: f64| t_model.map_or("—".to_string(), |tm| format!("{:+.2}%", 100.0 * (tm / t - 1.0)));
        if !s.vapour_pressure.is_empty() {
            let row: Vec<String> = s
                .vapour_pressure
                .iter()
                .map(|&(t, p)| format!("{:.0} kPa: T {}", p / 1e3, rel(m.coexistence_t(Phase::Liquid, Phase::Vapor, p, 0.3 * t, 3.0 * t), t)))
                .collect();
            println!("  boiling curve (T_model/T_data − 1): {}", row.join(", "));
        }
        if !s.sublimation_pressure.is_empty() {
            let row: Vec<String> = s
                .sublimation_pressure
                .iter()
                .map(|&(t, p)| format!("{p:.1} Pa: T {}", rel(m.coexistence_t(Phase::Solid, Phase::Vapor, p, 0.3 * t, 3.0 * t), t)))
                .collect();
            println!("  sublimation curve: {}", row.join(", "));
        }
        if let Some((tt, pt)) = s.triple {
            match m.triple_point(tt) {
                Some((t, p)) => println!("  triple point: model {t:.2} K, {:.4} kPa; data {tt:.2} K, {:.4} kPa", p / 1e3, pt / 1e3),
                None => println!("  triple point: not found"),
            }
        }
        if !s.melting_curve.is_empty() {
            let row: Vec<String> = s
                .melting_curve
                .iter()
                .map(|&(p, t)| {
                    let tm = m.coexistence_t(Phase::Solid, Phase::Liquid, p, 200.0, 300.0);
                    format!("{:.0} MPa: {} K (data {t:.1})", p / 1e6, tm.map_or("—".into(), |x| format!("{x:.2}")))
                })
                .collect();
            println!("  melting curve of ice (IAPWS R14-08): {}", row.join(", "));
            let dtdp = (m.coexistence_t(Phase::Solid, Phase::Liquid, P_ATM + 1e6, 200.0, 300.0).unwrap() - m.coexistence_t(Phase::Solid, Phase::Liquid, P_ATM, 200.0, 300.0).unwrap()) / 1e6;
            println!("  melting-point slope at 1 atm: {:.4} K/MPa (measured −0.074 K/MPa: ice is lighter than water)", dtdp * 1e6);
        }
        println!("  estimated inputs: {}", s.inputs.estimated.join("; "));
    }
}
