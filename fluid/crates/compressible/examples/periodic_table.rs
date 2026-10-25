//! Sweep the periodic table: fit each element's three-phase model and check
//! its predictions against data that was not used in the fit.
//!
//!   periodic_table [--verbose]
//!
//! Checks: standard state at 298.15 K / 1 atm (vs PubChem), triple point
//! (vs CRC, 17 elements), vapour-pressure curve 1 Pa … 100 kPa (vs the
//! measurements compiled on Wikipedia, 92 elements; reported as the error
//! in the temperature at which each pressure is reached).
use compressible::elements::{all, fit, vapor_pressure, Grade, Phase, P_ATM};

fn main() {
    let verbose = std::env::args().any(|a| a == "--verbose");
    let mut counts = [0usize; 5];
    let mut state_ok = 0;
    let mut state_n = 0;
    let mut state_bad = Vec::new();
    let mut vp_err: Vec<(Grade, f64, String)> = Vec::new();
    let mut tp_rows = Vec::new();
    println!("{:>3} {:<3} {:<14} grade  state(model/data)   vapour curve: rms/max |ΔT|/T   triple point model vs data", "Z", "", "name");
    for el in all() {
        let model = match fit(el) {
            Ok(m) => m,
            Err(e) => {
                counts[4] += 1;
                if verbose {
                    println!("{:>3} {:<3} {:<14} D      ({e})", el.z, el.symbol, el.name);
                }
                continue;
            }
        };
        counts[model.grade as usize] += 1;
        // standard state
        let ph = model.stable(P_ATM, 298.15);
        let ph_s = match ph {
            Phase::Solid => "Solid",
            Phase::Liquid => "Liquid",
            Phase::Vapor => "Gas",
        };
        let data_state = el.standard_state.clone();
        let measured_state = !data_state.contains("Expected");
        let ok = data_state.contains(ph_s);
        if measured_state {
            state_n += 1;
            if ok {
                state_ok += 1;
            } else {
                state_bad.push(format!("{} (model {ph_s}, data {data_state})", el.symbol));
            }
        }
        // vapour-pressure curve
        let vp = vapor_pressure(el.z);
        let mut errs = Vec::new();
        for &(p, t) in &vp {
            if let Some(tm) = model.vaporization_t(p, 0.1 * t, 5.0 * t) {
                errs.push((tm - t) / t);
            }
        }
        let vp_txt = if errs.is_empty() {
            "      —".to_string()
        } else {
            let rms = (errs.iter().map(|e| e * e).sum::<f64>() / errs.len() as f64).sqrt();
            let mx = errs.iter().map(|e| e.abs()).fold(0.0, f64::max);
            vp_err.push((model.grade, rms, el.symbol.clone()));
            format!("{:5.2}% / {:5.2}% ({} pts)", 100.0 * rms, 100.0 * mx, errs.len())
        };
        // triple point
        let tp_txt = match (el.t_triple, el.p_triple, model.triple_point(el.t_melt.unwrap())) {
            (Some(td), Some(pd), Some((tm, pm))) => {
                tp_rows.push((el.symbol.clone(), td, pd, tm, pm));
                format!("T {tm:8.2} vs {td:8.2} K, p {pm:9.3e} vs {pd:9.3e} Pa")
            }
            (Some(td), Some(pd), None) => format!("not found (data T {td:.2} K, p {pd:.3e} Pa)"),
            _ => String::new(),
        };
        println!(
            "{:>3} {:<3} {:<14} {:?}      {:<6}/{:<10} {:<28} {}",
            el.z,
            el.symbol,
            el.name,
            model.grade,
            ph_s,
            data_state.replace("Expected to be a ", "~"),
            vp_txt,
            tp_txt
        );
        if verbose {
            println!("        estimated: {}", model.estimated.join("; "));
        }
    }
    println!();
    println!("grades: A {} · B {} · C {} · P (predicted inputs) {} · D (not modelled) {}", counts[0], counts[1], counts[2], counts[3], counts[4]);
    println!("standard state at 298 K, 1 atm: {state_ok}/{state_n} correct{}", if state_bad.is_empty() { String::new() } else { format!("; wrong: {}", state_bad.join(", ")) });
    for g in [Grade::A, Grade::B, Grade::C, Grade::P] {
        let mut e: Vec<f64> = vp_err.iter().filter(|x| x.0 == g).map(|x| x.1).collect();
        if e.is_empty() {
            continue;
        }
        e.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!(
            "vapour curve, grade {:?}: {} elements, median rms |ΔT|/T {:.2}%, worst {:.2}%",
            g,
            e.len(),
            100.0 * e[e.len() / 2],
            100.0 * e[e.len() - 1]
        );
    }
    if !tp_rows.is_empty() {
        let et: Vec<f64> = tp_rows.iter().map(|r| ((r.3 - r.1) / r.1).abs()).collect();
        let ep: Vec<f64> = tp_rows.iter().map(|r| (r.4 / r.2).ln().abs()).collect();
        println!(
            "triple points ({}): median |ΔT|/T {:.2}%, median |ln(p_model/p_data)| {:.2} (factor {:.2})",
            tp_rows.len(),
            100.0 * median(&et),
            median(&ep),
            median(&ep).exp()
        );
    }
}

fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    s[s.len() / 2]
}
