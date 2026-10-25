//! One dated held-out synoptic-map comparison, with no future emergence.
mod common;
use solar::surface::{score, Parameters};
use std::time::Instant;
const TARGET: &str = include_str!("../data/gong_20190204.csv");
const INTERVAL: f64 = 27.0 * 86400.0 - 600.0;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // These defaults are fixed independently of the held-out target.
    let parameters = Parameters::default();
    let mut forecast = common::map(common::INITIAL, parameters)?;
    let initial_bias = forecast.balance_flux();
    let persistence = forecast.field().to_vec();
    let initial_dipole = forecast.axial_dipole();
    let initial_flux = forecast.net_flux();
    let start = Instant::now();
    let steps = forecast.advance(INTERVAL)?;
    let elapsed = start.elapsed().as_secs_f64();

    // Target values first enter the calculation after the forecast is finished.
    let mut observed = common::map(TARGET, parameters)?;
    let target_bias = observed.balance_flux();
    let skill = score(observed.field(), forecast.field(), &persistence)?;
    println!("GONG synoptic maps: 2019-01-08T11:14 -> 2019-02-04T11:04 UTC");
    println!(
        "grid=60x120 interval_days={:.9} steps={steps} runtime_s={elapsed:.6}",
        INTERVAL / 86400.0
    );
    println!("assumption=no new flux emergence; prescribed default flows; no target fitting");
    println!("preprocessing=aligned equal-area 3x3 means, Nyquist removal, explicit flux balance");
    println!(
        "removed_initial_mean_G={:.9} removed_target_mean_G={:.9}",
        initial_bias * 1e4,
        target_bias * 1e4
    );
    println!(
        "forecast_rmse_G={:.9} persistence_rmse_G={:.9} mse_skill={:?}",
        skill.rmse * 1e4,
        skill.persistence_rmse * 1e4,
        skill.mse_skill
    );
    println!(
        "axial_dipole_G: initial={:.9} forecast={:.9} observed={:.9}",
        initial_dipole * 1e4,
        forecast.axial_dipole() * 1e4,
        observed.axial_dipole() * 1e4
    );
    println!(
        "signed_flux_drift_Wb={:.9e}",
        forecast.net_flux() - initial_flux
    );
    println!("Synoptic maps span a rotation and include fitted polar fields. One case is not validated forecast skill.");
    if let Some(path) = std::env::args().nth(1) {
        let csv = forecast
            .field()
            .chunks(120)
            .map(|row| {
                row.iter()
                    .map(|v| format!("{v:.12e}"))
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .collect::<Vec<_>>()
            .join("\n");
        std::fs::write(path, csv + "\n")?;
    }
    Ok(())
}
